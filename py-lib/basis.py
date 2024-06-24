import yaml 
import re
import sys
import get_atomic_data as libatomic
import orbitals as liborb
import os
import collections.abc
from pathlib import Path
from subprocess import run
from gen_job_script import write_job_script

def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def count_total_orbitals(core_orbitals, valence_orbitals):
    """ Counts the total number of orbitals given lists of core and valence shells """
    num_core_orbitals = count_orbitals(core_orbitals)
    num_val_orbitals = count_orbitals(valence_orbitals)
    
    return num_core_orbitals + num_val_orbitals, num_core_orbitals, num_val_orbitals

def count_orbitals(orbitals):
    """ Counts the total number of orbitals from a list of orbitals """
    orbitals = orbitals.split()
    num_orbitals = 0
    for orbital in orbitals:
        try:
            l = re.findall('[spdfghikl]+', orbital)[0]
        except IndexError:
            print(orbital + ' is not a valid orbital')
        if l == 's':
            num_orbitals += 1
        else:
            num_orbitals += 2

    return num_orbitals

def get_key_breit(breit_bool):
    """ Returns key for breit from input whether to use breit or not """
    key_breit = int(re.sub('(no|No|n|N|false|False)', '0', re.sub('(yes|Yes|y|Y|true|True)', '2', str(breit_bool))))
    return key_breit

def get_key_vw(kvw_str):
    """ Returns key kvw used in inf.vw """
    kvw = -1

    if kvw_str == 'ci+all-order' or kvw_str == 'all-order' or kvw_str == 'all':
        kvw = 1
    elif kvw_str == 'ci+second-order' or kvw_str == 'second-order' or kvw_str == 'second':
        kvw = 0

    return kvw

def gen_lists_orbitals(core_orbitals, valence_orbitals):
    """ Generates a list of orbitals with possible values of J for HFD.INP"""
    NL, J, QQ, KP, NC = [], [], [], [], []
    num_core_electrons = 0
    for orbital in core_orbitals.split():
        NL.append(orbital[0] + orbital[1].upper())
        NC.append('0')
        KP.append('0')
        if orbital[1] == 's':
            J.append('1/2')
            QQ.append('2.0000')
            num_core_electrons += 2
        else:
            NL.append(orbital[0] + orbital[1].upper())
            NC.append('0')
            KP.append('0')
            if orbital[1] == 'p':
                J.append('1/2')
                J.append('3/2')
                QQ.append('2.0000')
                QQ.append('4.0000')
                num_core_electrons += 6
            if orbital[1] == 'd':
                J.append('3/2')
                J.append('5/2')
                QQ.append('4.0000')
                QQ.append('6.0000')
                num_core_electrons += 10
            if orbital[1] == 'f':
                J.append('5/2')
                J.append('7/2')
                QQ.append('6.0000')
                QQ.append('8.0000')
                num_core_electrons += 14

    count = 0
    nval = 0
    for orbital in valence_orbitals.split():
        count += 1
        nval += 1
        NL.append(orbital[0] + orbital[1].upper())
        NC.append(str(count))
        KP.append('0')
        if orbital[1] == 's':
            J.append('1/2')
            QQ.append('1.0000')
        else:
            nval += 1
            NL.append(orbital[0] + orbital[1].upper())
            NC.append(str(count))
            KP.append('0')
            QQ.append('1.0000')
            QQ.append('0.0000')
            if orbital[1] == 'p':
                J.append('1/2')
                J.append('3/2')
            if orbital[1] == 'd':
                J.append('3/2')
                J.append('5/2')
            if orbital[1] == 'f':
                J.append('5/2')
                J.append('7/2')
            if orbital[1] == 'g':
                J.append('7/2')
                J.append('9/2')

    return NL, J, QQ, KP, NC, num_core_electrons, nval

def gen_lists_kappa(Z, num_core_electrons, core_orbitals, valence_orbitals):
    """ Generates a list of relativistic quantum numbers kappa and energy guesses for bas_wj.in """
    N, kappa, iters, energies = [], [], [], []
    lowest_n = {}
    
    for orbital in core_orbitals.split():
        N.append(orbital[0])
        iters.append('0')
        energies.append(0.00)
        if orbital[1] == 's':
            kappa.append('-1')
        else:
            N.append(orbital[0])
            iters.append('0')
            energies.append(0.00)
            if orbital[1] == 'p':
                kappa.append(' 1')
                kappa.append('-2')
            if orbital[1] == 'd':
                kappa.append(' 2')
                kappa.append('-3')
            if orbital[1] == 'f':
                kappa.append(' 3')
                kappa.append('-4')
            if orbital[1] == 'g':
                kappa.append(' 4')
                kappa.append('-5')
    
    # Retrieve lowest n for each l for valence orbitals
    for orbital in valence_orbitals.split():
        if orbital[1] not in lowest_n:
            lowest_n[orbital[1]] = int(orbital[0])

    for orbital in valence_orbitals.split():
        N.append(orbital[0])
        iters.append('1')
        energies.append(get_energy_guess(Z-num_core_electrons, int(orbital[0]), lowest_n[orbital[1]], orbital[1]))
        if orbital[1] == 's':
            kappa.append('-1')
        else:
            N.append(orbital[0])
            iters.append('1')
            energies.append(get_energy_guess(Z-num_core_electrons, int(orbital[0]), lowest_n[orbital[1]], orbital[1]))
            if orbital[1] == 'p':
                kappa.append(' 1')
                kappa.append('-2')
            if orbital[1] == 'd':
                kappa.append(' 2')
                kappa.append('-3')
            if orbital[1] == 'f':
                kappa.append(' 3')
                kappa.append('-4')
            if orbital[1] == 'g':
                kappa.append(' 4')
                kappa.append('-5')    

    return N, kappa, iters, energies

def get_energy_guess(deg_ion, n, n0, l):
    """ 
    Returns the prediction for DHF energy 

    Inputs:
    deg_ion = degree of ionization (1 for Cs I, 2 for Ba II and so on)
    n = principal quantum number
    l = orbital quantum number
    n0 = lowest principal quantum number
    """
    # Inital guess for neutral atom and lowest principal quantum number for spdfg
    initial_guess = {'s': -0.13, 'p': -0.09, 'd': -0.064, 'f': -0.03, 'g': -0.02}
    # Scaling factor to use to account for higher principal quantum numbers
    scaling_factor = {'s': 2.0, 'p': 1.8, 'd': 1.7, 'f': 1.45, 'g': 1.3}

    # Difference between actual and lowest non-core quantum number
    qn_diff = n - n0

    # Empirical formula to account for higher n and degree of ionization
    energy = initial_guess[l] * deg_ion**2 / (scaling_factor[l]**qn_diff)

    # Use different variation of this formula for different degrees of ionization
    if l == 's' or l == 'p':
        if deg_ion > 1 and deg_ion < 10:
            energy = initial_guess[l]*((deg_ion*0.8)**2)/(scaling_factor[l]**qn_diff)
        elif deg_ion > 10:
            energy = initial_guess[l]*((deg_ion*0.5)**2)/(scaling_factor[l]**qn_diff)
        elif deg_ion > 20:
            energy = initial_guess[l]*((deg_ion*0.4)**2)/(scaling_factor[l]**qn_diff)
    elif l == 'd':
        energy = initial_guess[l]*((deg_ion*0.7)**2)/(scaling_factor[l]**qn_diff)

    # Do not allow very small energy predictions, set default to -0.001
    if -energy < 0.001: energy = -0.001

    return energy


def get_ao_valence(core_orbitals, valence_orbitals, val_aov):
    """ Returns the lowest spdf orbitals in kappa designation """
    N, kappa = [], []
    nmin = [1, 2, 3, 4, 5]
    nmax = [0, 0, 0, 0, 0]
    for orbital in core_orbitals.split():
        n = int(re.findall('[0-9]+', orbital)[0])
        if orbital[-1] == 's':
            nmin[0] = n + 1
        elif orbital[-1] == 'p':
            nmin[1] = n + 1
        elif orbital[-1] == 'd':
            nmin[2] = n + 1
        elif orbital[-1] == 'f':
            nmin[3] = n + 1
        elif orbital[-1] == 'g':
            nmin[4] = n + 1
        else:
            pass

    # Set maximum n for each l (4 lowest s, p, d; 3 lowest f) by default
    if val_aov == {}:
        nmax = [n + 4 for n in nmin]
        nmax[3] -= 1
        nmax[4] -= 4
    else:
        n=0
        for l in val_aov:
            nmax[n] = nmin[n] + int(val_aov[l])
            n += 1

    # Set kappas for valence orbitals
    for n in range(len(nmax)):
        if n == 0:
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append('-1')
        elif n == 1:
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append(' 1')
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append('-2')
        elif n == 2:
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append(' 2')
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append('-3')
        elif n == 3:
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append(' 3')
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append('-4')
        elif n == 4:
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append(' 4')
            for i in range(nmin[n],nmax[n]):
                N.append(i)
                kappa.append('-5')        

    """
    # Find orbital # for 23f7/2 
    nmax = [23, 23, 23, 23]
    Nmax = nmax[0] - nmin[0]            # number of s orbitals
    print(Nmax)
    Nmax += (nmax[1] - nmin[1]) * 2     # number of p orbitals 
    print(Nmax)
    Nmax += (nmax[2] - nmin[2]) * 2     # number of d orbitals
    print(Nmax)
    Nmax += (nmax[3] - nmin[3]) * 2     # number of f orbitals
    print(Nmax)
    Nmax -= NSO                         # remove core orbitals
    print(Nmax)
    """

    return N, kappa

def write_hfd_inp(filename, system, NS, NSO, Z, AM, kbr, NL, J, QQ, KP, NC, rnuc, K_is, C_is):
    """ Writes HFD.INP """
    # Define default values
    KL = 0
    JM = -2.0

    with open(filename, 'w') as f:
        f.write(' ' + system['system']['name'] + '\n')
        f.write(' KL =  ' + str(KL) + '\n')
        f.write(' NS =  ' + str(NS) + '\n')
        f.write(' NSO=  ' + str(NSO) + '\n')
        f.write(' Z  =  ' + str(Z) + '\n')
        try:
            f.write(' AM =  ' + '{:.2f}'.format(round(AM)) + '\n')
        except:
            print('AM could not be found')
            sys.exit()
        f.write(' JM =  ' + str(JM) + '\n')
        f.write(' R2 =  ' + str(float(system['basis']['cavity_radius'])) + '\n')
        f.write(' kbr= ' + str(kbr) + '\n')
        if C_is != 0:
            f.write('K_is= ' + str(K_is) + '\n')
            f.write('C_is= ' + str(C_is) + '\n')
        f.write('rnuc= ' + '{:.4f}'.format(rnuc) + '\n\n')
        f.write('        NL   J       QQ     KP   NC\n\n')
        for i in range(len(NL)):
            f.write(str(i+1).rjust(3," ") + '     ' + NL[i] + ' (' + J[i] + ')' + '   ' 
                        + QQ[i] + '    ' + KP[i] + '   ' + NC[i].rjust(2,' ') + '\n')
    f.close()
    print(filename + ' has been written')

def construct_vvorbs(core, valence, codename, nmax, lmax):
# Construct list of valence and virtual orbitals
    vorbs = []
    norbs = []

    ## First we need to construct list of orbitals from nmax + lmax
    nmin = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    core_val = core + ' ' + valence
    l = ['s','p','d','f','g','h','i','k','l']
    for orbital in core_val.split():
        n = int(re.findall('[0-9]+', orbital)[0])
        if orbital[-1] == 's' and n > nmin[0]:
            nmin[0] = n
        elif orbital[-1] == 'p' and n > nmin[1]:
            nmin[1] = n
        elif orbital[-1] == 'd' and n > nmin[2]:
            nmin[2] = n
        elif orbital[-1] == 'f' and n > nmin[3]:
            nmin[3] = n
        elif orbital[-1] == 'g' and n > nmin[4]:
            nmin[4] = n
        elif orbital[-1] == 'h' and n > nmin[5]:
            nmin[5] = n
        elif orbital[-1] == 'i' and n > nmin[6]:
            nmin[6] = n
        elif orbital[-1] == 'k' and n > nmin[7]:
            nmin[7] = n
        elif orbital[-1] == 'l' and n > nmin[8]:
            nmin[8] = n
        else:
            pass

    ### Set minimum n for each l, effectively removing core orbitals
    for n in range(len(nmin)):
        if nmin[n] <= n:
            nmin[n] = n
        nmin[n] += 1

    # Write valence orbitals
    if codename == 'ci':
        if len(valence.split()) != 0:
            count = 0
            for i in range(0, len(valence.split())):
                orbs = liborb.convert_char_to_digital(valence.split()[i])
                for orb in orbs:
                    count += 1
                    vorbs += [orb]
                    norbs += [orb[-6]+orb[-4]]
    else:
        if len(valence.split()) != 0:
            count = 0
            for i in range(0, len(valence.split())):
                orbs = liborb.convert_char_to_digital(valence.split()[i])
                for orb in orbs:
                    count += 1
                    vorbs += [orb]
                    norbs += [orb[-6]+orb[-4]]
    nval = count

    ### Loop over n to nmax, l to lmax
    for n in range(nmax+1):
        for i in range(lmax+1):
            if (nmin[i] <= nmin[lmax] and n == nmin[i]): 
                strorb = str(nmin[i]) + str(l[i])
                orbs = liborb.convert_char_to_digital(strorb)
                for orb in orbs:
                    count += 1
                    vorbs += [orb]
                    norbs += [orb[-6]+orb[-4]]
                nmin[i] += 1
    nvvorbs = count

    return vorbs, norbs, nval, nvvorbs


def write_bass_inp(filename, system, NSO, Z, AM, kbr, vorbs, norbs, nmax, lmax, codename, core, valence, K_is, C_is):
    """ Writes BASS.INP """
    # Define default values
    Nv = len(valence.split())
    Ksg = 1
    Kdg = 1
    Kkin = 1
    kout = 0

    vorbs, norbs, nval, nvvorbs = construct_vvorbs(core, valence, codename, nmax, lmax)

    # Set first orbital for diagonalization to be first valence orbital
    fvalorb = vorbs[0]
    fvalorb_str = liborb.convert_digital_to_char(fvalorb)
    fvalorb = fvalorb_str[0] + " " + fvalorb_str[1][0]

    # Set first orbital to apply kinetic balance to be first virtual orbital
    fvirorb = vorbs[nval]
    fvirorb_str = liborb.convert_digital_to_char(fvirorb)
    fvirorb = fvirorb_str[0] + " " + fvirorb_str[1][0]

    # Set last frozen orbital to be last orbital in core
    frorb_str = liborb.convert_digital_to_char(liborb.convert_char_to_digital(core.split()[-1])[-1])
    frorb = frorb_str[0] + " " + frorb_str[1][0]    

    with open(filename, 'w') as f:
        f.write(' ' + system['system']['name'] + '\n')
        f.write(' Z  =  ' + str(Z) + '\n')
        f.write(' Am =  ' + '{:.1f}'.format(round(AM)) + '\n')
        f.write(' Nso=' + str(NSO).rjust(5," ") + '# number of core orbitals (defines DF operator)\n')
        f.write(' Nv =' + str(nvvorbs).rjust(5," ") + '# number of valence & virtual orbitals\n')
        f.write(' Ksg=' + str(Ksg).rjust(5," ") + '# defines Hamiltonian: 1-DF, 3-DF+Breit\n')
        f.write(' Kdg=' + str(Kdg).rjust(5," ") + '# diagonalization of Hamiltonian (0=no,1,2=yes)\n')
        f.write(' orb=' + fvalorb.rjust(5," ") + '# first orbital for diagonalization\n')
        f.write(' Kkin' + str(Kkin).rjust(5," ") + '# kinetic balance (0,1,or 2)\n')
        f.write(' orb=' + fvirorb.rjust(5," ") + '# first orbital to apply kin.bal.\n')
        f.write(' orb=' + frorb.rjust(5," ") + '# last frozen orbital\n')
        f.write('kout=' + str(0).rjust(5," ") + '# detail rate in the output\n')
        f.write('kbrt=' + str(kbr).rjust(5," ") + '# 0,1,2 - Coulomb, Gaunt, Breit\n')
        if C_is != 0:
            f.write('K_is= ' + str(K_is) + '\n')
            f.write('C_is= ' + str(C_is) + '\n')
        f.write('----------------------------------------------------------\n')
        # Write core orbitals
        if len(core.split()) != 0:
            count = 0
            for i in range(0, len(core.split())):
                orbs = liborb.convert_char_to_digital(core.split()[i])
                for orb in orbs:
                    if orb[0] == '-':
                        f.write('    ' + orb)
                    else:
                        f.write('     ' + orb)
                    count += 1
                    if count == 6:
                        f.write('\n')
                        count = 0
            f.write('\n\n')

        # Write valence and virtual orbitals
        for i, orb in enumerate(vorbs):
            orb = orb[:-2] + '01'
            if (i < nval):
                f.write(str(i+1).rjust(3," ") + orb.rjust(8," ") + "\n")
            else:
                f.write(str(i+1).rjust(3," ") + orb.rjust(8," ") + "  3 " + orb.rjust(7," ") + "\n")
            
            # Write a new line in between different n
            if (i == nval-1):
                f.write(" \n")
            try:
                if (i >= nval and norbs[i+1] != norbs[i]):
                    f.write(" \n")
            except IndexError:
                f.write('>>>>>>>>>> END <<<<<<<<<<')
    f.close()
    print(filename + ' has been written')

def write_bas_wj_in(filename, symbol, Z, AM, NS, NSO, N, kappa, iters, energies, cfermi):
    """ Writes bas_wj.in """
    with open(filename,'w') as f: 
        f.write(' ' + symbol + ' ' + str(NS).rjust(4, ' ') + str(int(Z)).rjust(4, ' ') 
                    + str(int(AM)).rjust(4,' ') + '   0   0   9' + str(NSO+1).rjust(4, ' ')
                    + '   1\n')
        for i in range(len(N)):
            f.write(N[i].rjust(4, ' ') + kappa[i].rjust(4, ' ') + iters[i].rjust(4, ' ') +
                         "{:.2f}".format(energies[i]).rjust(7, ' ') + '\n')
        if symbol == 'Na':
            f.write('   1.5\n')
        elif symbol == 'Yb':
            f.write('   2.0\n')
        else:
            f.write('   1.0\n')
        grid = 0.0002
        f.write(str(grid).rjust(9, ' ') + '  0.03  500\n')
        f.write('   1\n')
        f.write('   0.0000' + str(round(cfermi, 4)).rjust(10, ' ') + '    2.3\n')
        f.write('   0.0')
    f.close()
    print('bas_wj.in has been written')

def write_inf_aov(filename, val_N, val_kappa, NSO, nmax, lmax, kval, energies):
    """ Writes inf.aov """
    with open(filename,'w') as f: 
        f.write(str(NSO) + '\n')
        for i in range(NSO):
            f.write(N[i] + ' ' + kappa[i] + '\n')
        f.write(str(nmax) + '  ' + str(lmax) + '\n')
        f.write('0   0\n') # internal parameters
        f.write('30\n') # max iterations
        f.write('2 7 1\n') # stabilizer code parameters
        f.write('0.d0\n') # damping factor
        f.write(str(kval) + '\n') # kval
        # If kval = 2, write energies 
        if kval == 2:
            f.write(str(len(energies)-1) + '\n')
            n = 0
            for l in energies:
                energy = energies[l]
                if isinstance(energy, list):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy[0]) + '  ' + "{:.5f}".format(energy[1]) + '\n')
                elif isinstance(energy, float):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy) + '\n')
            n += 1
        f.write(str(len(val_kappa)) + '\n')
        for i in range(len(val_kappa)):
            f.write(str(val_N[i]).rjust(2, ' ') + '  ' + val_kappa[i] + ' 30\n' )
    f.close()
    print('inf.aov has been written')
        

def write_inf_vw(filename, val_N, val_kappa, NSO, nmax, lmax, kvw, kval, energies):
    """ Writes inf.vw """
    l_array = ['s', 'p', 'd', 'f']
    with open(filename,'w') as f: 
        f.write(str(NSO) + '\n')
        for i in range(NSO):
            f.write(N[i] + ' ' + kappa[i] + '\n')
        f.write(str(nmax) + ' ' + str(lmax) + '\n')
        f.write('250\n')
        f.write('4\n')
        f.write('9\n')
        f.write('250 100\n')
        f.write(str(kvw) + '\n')
        f.write(str(kval) + '\n')
        # If kval = 2, write energies 
        if kval == 2:
            f.write(str(len(energies)-1) + '\n')
            n = 0
            for l in l_array:
                energy = energies[l]
                if isinstance(energy, list):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy[0]) + '  ' + "{:.5f}".format(energy[1]) + '\n')
                elif isinstance(energy, float):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy) + '\n')
            n += 1
    f.close()
    print('inf.vw has been written')

def write_spl_in(filename, radius, spl_params):
    """ Writes spl.in """
    
    with open(filename,'w') as f: 
        f.write(str(spl_params['lmax']) + '\n') 
        f.write(str(radius) + '\n')  # Cavity radius
        f.write(str(spl_params['nmax']) + ' ' + str(spl_params['k']) + '\n')  # Number of splines and order of splines
        f.write('0.0 0.00 500')

    f.close()
    print('spl.in has been written')

def write_inputs(system, C_is, kvw):
    # Write HFD.INP
    write_hfd_inp('HFD.INP', system, NS, NSO, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, system['optional']['isotope_shifts']['K_is'], C_is)
    
    # Write BASS.INP
    write_bass_inp('BASS.INP', system, NSO, Z, AM, kbrt, vorbs, norbs, system['basis']['orbitals']['nmax'], system['basis']['orbitals']['lmax'], system['optional']['code_method'], system['basis']['orbitals']['core'], system['basis']['orbitals']['valence'], system['optional']['isotope_shifts']['K_is'], C_is)

    # Write bas_wj.in
    write_bas_wj_in('bas_wj.in', symbol, Z, AM, NS, NSO, N, kappa, iters, energies, cfermi)

    # Write spl.in
    write_spl_in('spl.in', system['basis']['cavity_radius'], system['basis']['b_splines'])
    
    # Write inf.aov
    write_inf_aov('inf.aov', val_N, val_kappa, NSO, system['basis']['orbitals']['nmax'], system['basis']['orbitals']['lmax'], kval, system['basis']['val_energies']['energies'])

    # Write inf.vw
    write_inf_vw('inf.vw', val_N, val_kappa, NSO, system['basis']['orbitals']['nmax'], system['basis']['orbitals']['lmax'], kvw, kval, system['basis']['val_energies']['energies'])

def generate_batch_qed(kqed, krot, kbrt):
    """ Writes batch.qed """
    with open('q.in','w') as f:
        f.write("1 \n")
        if kqed == True:
            f.write('2 \n')
        else:
            f.write('1 \n')
        f.write(str(kbrt) + '\n')
    f.close()

    with open('batch.qed','w') as f: 
        f.write("#! /bin/bash -fe \n")
        f.write("vpkg_require pci \n")
        f.write("kvar=1  # variant of QED potential \n")
        f.write("iter=25 # max number of iterations \n")
        f.write("##################################### \n")
        f.write("cat >qedpot.inp <<EndofFile\n")
        f.write(" $kvar \n")
        f.write(" HFD.DAT \n")
        f.write("EndofFile \n")
        f.write("##################################### \n")
        f.write("cat >q.in <<EndofFile\n")
        f.write("1 \n")
        if kqed == True:
            f.write('2 \n')
        else:
            f.write('1 \n')
        f.write(str(kbrt) + '\n')
        f.write(". \n")
        f.write(". \n")
        f.write("  \n")
        f.write("EndofFile\n")
        f.write("##################################### \n")
        f.write("n=1 \n")
        f.write("while [ $n -lt $iter ]; do \n")
        f.write("echo 'Iteration '$n \n")
        f.write("qedpot_conf <q.in >qp.res \n")
        if krot == True:
            f.write("qed_rot <q.in >qr.res \n")
        f.write("grep 'changed' \"QED_ROT.RES\" \n")
        f.write("  if grep -q reached \"QED_ROT.RES\"; then \n")
        f.write("  echo 'Converged in '$n' iterations' \n")
        f.write("  break \n")
        f.write("  fi \n")
        f.write("  let n=n+1 \n")
        f.write("done \n")
        f.write("##################################### \n")
    print('batch.qed has been written')
    
def check_errors(filename):
    # This function checks output files for errors and returns the number of errors
    # Currently only supports output of program bass
    num_errors = 0
    warning_msgs = ['failed', 'error', 'warning']
    if filename == 'bass.out':
        try:
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if any(msg in line for msg in warning_msgs):
                    num_errors += 1
        except FileNotFoundError:
            num_errors = 1   
        
        return num_errors
    else:
        print(filename + "not currently supported")

def run_executables(K_is, C_is):
    # Run hfd
    run('hfd > hfd.out', shell=True)
    print("hfd complete")

    # Produce B-splines
    if kbrt == 0:
        run('tdhf < bas_wj.in > tdhf.out', shell=True)
        print("tdhf complete")
        run('nspl40 < spl.in > nspl40.out', shell=True)
        print("nspl40 complete")
    else:
        run('bdhf < bas_wj.in > bdhf.out', shell=True)
        print("bdhf complete")
        run('bspl40 < spl.in > bspl40.out', shell=True)
        print("bspl40 complete")

    with open('bwj.in','w') as f: 
        f.write(str(basis_lmax) + '\n')
        for i in range(basis_lmax+1):
            f.write(str(basis_nmax) + '\n')
        f.write('\n')
        f.write('\n')
        f.write('1')
    f.close()

    run('bas_wj < bwj.in > bas_wj.out', shell=True)
    run(['rm','bwj.in'])
    print("bas_wj complete")

    # Edit BASS.INP
    f = open('BASS.INP', 'r')
    lines = f.readlines()
    f.close()

    # Set first orbital for diagonalization to be first valence orbital
    fvalorb = vorbs[0]
    fvalorb_str = liborb.convert_digital_to_char(fvalorb)
    fvalorb = fvalorb_str[0] + " " + fvalorb_str[1][0]

    # Set first orbital to apply kinetic balance to be first virtual orbital
    fvirorb = vorbs[nvalb]
    fvirorb_str = liborb.convert_digital_to_char(fvirorb)
    fvirorb = fvirorb_str[0] + " " + fvirorb_str[1][0]

    # Set last frozen orbital to be last orbital in core
    frorb_str = liborb.convert_digital_to_char(liborb.convert_char_to_digital(core_orbitals.split()[-1])[-1])
    frorb = frorb_str[0] + " " + frorb_str[1][0]   
    
    f = open('BASS.INP','w')
    for line in lines:
        if line[:4] == ' lst':
            continue
        elif 'first' in line.split() and 'diagonalization' in line.split():
            f.write(' orb=' + fvalorb.rjust(5," ") + '# first orbital for diagonalization\n')
        elif 'kin.bal.' in line.split():
            f.write(' orb=' + fvirorb.rjust(5," ") + '# first orbital to apply kin.bal.\n')
        elif 'frozen' in line.split():
            f.write(' orb=' + frorb.rjust(5," ") + '# last frozen orbital\n')
        elif line[:4] == 'kbrt':
            f.write('kbrt=' + str(kbrt).rjust(2," ") + '   # 0,1,2 - Coulomb, Gaunt, Breit\n')
            if K_is != 0:
                f.write('K_is=' + str(K_is).rjust(2," ") + '\n')
                f.write('C_is=' + str(C_is) + '\n')
        elif line[:3].strip().isdigit() and int(line[:3].strip()) <= nval:
            f.write(line[:12] + "\n")
        else:
            f.write(line)
    f.close()

    # run bass until there are no errors in output
    with open('bass.in','w') as f: 
        f.write('WJ.DAT')
    f.close()

    # check if bass.out exists and remove if it does
    if os.path.isfile('bass.out'):
        run(['rm','bass.out'])

    maxNumTries = 5
    nTry = 1

    while check_errors('bass.out') > 0:
        print('bass attempt', nTry)
        run('bass < bass.in > bass.out', shell=True)
        if (nTry >= maxNumTries):
            print("bass did not converge after", nTry, "attempts")
            break
        nTry += 1
    else:
        print("bass completed with no errors after", nTry, "attempts")

    run(['rm', 'bass.in'])
    run(['rm','hfspl.1','hfspl.2'])

    # Run qed
    #if system['rotate_basis'] == True or system['include_qed'] == True:
    #    run('cp HFD.DAT HFD-noQED.DAT', shell=True)
    #    generate_batch_qed(system['include_qed'],system['rotate_basis'],kbrt)
    #    run('chmod +x batch.qed', shell=True)
    #    run('./batch.qed > qed.out', shell=True)
    #    print("qed complete")

    # Run bas_x
    run('bas_x > bas_x.out', shell=True)
    print("bas_x complete")

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)

    # Set parameters from config
    name = config['system']['name']
    isotope = config['system']['isotope']
    core_orbitals = config['basis']['orbitals']['core']
    valence_orbitals = config['basis']['orbitals']['valence']
    include_breit = config['system']['include_breit']
    basis_nmax = config['basis']['orbitals']['nmax']
    basis_lmax = config['basis']['orbitals']['lmax']
    kval = config['basis']['val_energies']['kval']
    val_aov = config['basis']['val_aov']

    include_isotope_shifts = config['optional']['isotope_shifts']['include']
    if include_isotope_shifts:
        K_is = config['optional']['isotope_shifts']['K_is']
        C_is = config['optional']['isotope_shifts']['C_is']
        c_list = [-C_is,-C_is/2,0,C_is/2,C_is]
        K_is_dict = {0: '', 1: 'FS', 2: 'SMS', 3: 'NMS', 4: 'MS'}

    code_method = config['optional']['code_method']
    run_ao_codes = config['optional']['run_ao_codes']
    pci_version = config['optional']['pci_version']

    # Get atomic data
    Z, AM, symbol, cfermi, rnuc, num_rem_ele = libatomic.get_atomic_data(name, isotope)

    # Get orbital information
    NS, NSO, num_val_orbitals = count_total_orbitals(core_orbitals, valence_orbitals)

    # Get breit key
    kbrt = get_key_breit(include_breit)

    # Generate body of HFD.INP including orbitals, values of J, and occupation numbers
    NL, J, QQ, KP, NC, num_core_electrons, nval = gen_lists_orbitals(core_orbitals, valence_orbitals)
    
    # Assign keys for which inputs to generate
    kvw = get_key_vw(code_method)

    vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(core_orbitals, valence_orbitals, code_method, basis_nmax, basis_lmax)

    # Generate body of bas_wj.in including orbitals, values of kappa, and energy guesses
    N, kappa, iters, energies = gen_lists_kappa(Z, num_core_electrons, core_orbitals, valence_orbitals)

    # Get valence orbitals for all-order calculations
    val_N, val_kappa = get_ao_valence(core_orbitals, valence_orbitals, val_aov)

    # Write input files to basis directory
    if include_isotope_shifts and K_is > 0:
        for method in code_method:
            dir_path = os.getcwd()
            is_dir = method + '/' + K_is_dict[K_is]
            Path(dir_path+'/'+is_dir).mkdir(parents=True, exist_ok=True)
            os.chdir(dir_path+'/'+is_dir)
            for c in c_list:
                dir_path = os.getcwd()
                if c < 0:
                    dir_prefix = 'minus' 
                elif c > 0:
                    dir_prefix = 'plus'
                else:
                    dir_prefix = ''
                dir_name = dir_prefix+str(abs(c))+'/basis'
                Path(dir_path+'/'+dir_name).mkdir(parents=True, exist_ok=True)
                os.chdir(dir_name)
                run('pwd', shell=True)
                write_inputs(config,c,get_key_vw(method))
                os.chdir('../../')
            if K_is_dict[K_is]:
                os.chdir('../../')
            else:
                os.chdir('../')
    else:
        if isinstance(code_method, list):
            dir_path = os.getcwd()
            for method in code_method:
                Path(dir_path+'/'+method+'/basis').mkdir(parents=True, exist_ok=True)
                os.chdir(method+'/basis')
                run('pwd', shell=True)
                write_inputs(config, 0, get_key_vw(method))
                os.chdir('../../')
        else:
            dir_path = os.getcwd()
            Path(dir_path+'/basis').mkdir(parents=True, exist_ok=True)
            os.chdir('basis')
            run('pwd', shell=True)
            write_inputs(config, 0, kvw)
            os.chdir('../')

    # Construct basis set by running sequence of programs if desired
    if run_ao_codes:
        print("Running codes...")
        if include_isotope_shifts and K_is > 0:
            for method in code_method:
                dir_path = os.getcwd()
                is_dir = method + '/' + K_is_dict[K_is]
                Path(dir_path+'/'+is_dir).mkdir(parents=True, exist_ok=True)
                os.chdir(dir_path+'/'+is_dir)
                for c in c_list:
                    dir_path = os.getcwd()
                    if c < 0:
                        dir_prefix = 'minus' 
                    elif c > 0:
                        dir_prefix = 'plus'
                    else:
                        dir_prefix = ''
                    dir_name = dir_prefix+str(abs(c))+'/basis'
                    os.chdir(dir_name)
                    run('pwd', shell=True)
                    run_executables(K_is, c)
                    script_name = write_job_script('.', method, 1, 1, True, 0, 'standard', pci_version)
                    run('sbatch ' + script_name, shell=True)
                    os.chdir('../../')
                if K_is_dict[K_is]:
                    os.chdir('../../')
                else:
                    os.chdir('../')
        else:
            if isinstance(code_method, list):
                for method in code_method:
                    dir_path = os.getcwd()
                    Path(dir_path+'/'+method+'/basis').mkdir(parents=True, exist_ok=True)
                    os.chdir(method+'/basis')
                    run('pwd', shell=True)
                    run_executables(0, 0)
                    script_name = write_job_script('.', method, 1, 1, True, 0, 'standard', pci_version)
                    run('sbatch ' + script_name, shell=True)
                    os.chdir('../../')
            else:
                dir_path = os.getcwd()
                Path(dir_path+'/basis').mkdir(parents=True, exist_ok=True)
                os.chdir('basis')
                run('pwd', shell=True)
                run_executables(0, 0)
                script_name = write_job_script('.', code_method, 1, 1, True, 0, 'standard', pci_version)
                run('sbatch ' + script_name, shell=True)
                os.chdir('../')
                
