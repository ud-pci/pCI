import yaml
import re
import math
import sys
import get_atomic_data as libatomic
import orbitals as liborb
from subprocess import Popen, PIPE, STDOUT, run

def read_yaml(filename):
    """ Reads yaml input file and returns contents """ 

    isotope = 0
    val_aov = []
    # read yaml file with inputs
    with open(filename,'r') as f:
        config = yaml.safe_load(f)

        # Check to see if isotope is specified. If not, set to 0
        try:
            isotope = config['isotope']
        except KeyError as e:
            config[e.args[0]] = 0
        
        # Check to see if val_aov is specified. If not, set to []
        try:
            val_aov = config['val_aov']
        except KeyError as e:
            config[e.args[0]] = []

        # Check to see if energies are specified. If not, set to []
        try:
            energies = config['energies']
        except KeyError as e:
            config[e.args[0]] = []

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
            l = re.findall('[spdfgh]+', orbital)[0]
        except IndexError:
            print(orbital + ' is not a valid orbital')
        if l == 's':
            num_orbitals += 1
        else:
            num_orbitals += 2

    return num_orbitals

def get_atomic_data(name, isotope):
    """ Gets atomic number Z and atomic mass AM from periodic table """
    ptable = libatomic.get_periodic_table()
    rtable = libatomic.get_radii()
    symbol = re.findall('[a-zA-Z]+', name)[0]
    rnuc = 0
    cfermi = 0

    # Check if atomic symbol in periodic table
    if ptable['Symbol'].eq(symbol).any():
        Z = ptable[ptable['Symbol'] == symbol]['Z'].values.astype(float)[0]
        AM = ptable[ptable['Symbol'] == symbol]['Mass'].values.astype(float)[0]

    # Check if isotope was specified
    if isotope != 0:
        AM = isotope
    else:
        print('Isotope was not specified, so default mass is used')

    # Check if atomic symbol in radii table
    if rtable['Elem.'].eq(symbol).any():
        try: 
            rnuc = rtable[(rtable['Elem.'] == symbol) & (rtable['Mass'] == str(round(AM)))]['R_av(fm)'].values.astype(float)[0]
        except IndexError:
            pass
    
    
    # If rnuc could not be found, search manual table for rnuc
    if rnuc == 0:
        rtable2 = libatomic.get_extra_radii()
        if rtable2['Elem.'].eq(symbol).any():
            rnuc = rtable2[(rtable2['Elem.'] == symbol)]['R_av(fm)'].values.astype(float)[0]
            print('rnuc was not found in nuclear charge radii table, so a guess is used.')
        else:
            print('Could not find rnuc in any tables. Enter parameters in HFD.INP and bas_wj.in manually.')
            pass
    else:
        # Calculate cfermi
        try:
            cfermi = calc_c_fermi(rnuc)
        except UnboundLocalError as e:
            print('ERROR: Nuclear charge radius could not be found for', name)
            sys.exit()
        except ValueError as e:
            print('ERROR: Nuclear charge radius could not be found for', name)
            sys.exit()

    return Z, AM, symbol, cfermi, rnuc

def calc_c_fermi(rnuc):
    """ Calculates c fermi from rms radius """
    at = 2.3/(4*math.log(3.0))
    a1 = (5/3)*(rnuc**2)
    a2 = (7/3)*((math.pi*at)**2)
    cfermi = math.sqrt(a1-a2)

    return cfermi

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
    nmin = [0, 0, 0, 0]
    nmax = [0, 0, 0, 0]
    for orbital in core_orbitals.split():
        n = int(re.findall('[0-9]+', orbital)[0])
        if orbital[-1] == 's' and n > nmin[0]:
            nmin[0] = n
        elif orbital[-1] == 'p' and n > nmin[1]:
            nmin[1] = n
        elif orbital[-1] == 'd' and n > nmin[2]:
            nmin[2] = n
        elif orbital[-1] == 'f' and n > nmin[3]:
            nmin[3] = n
        else:
            pass

    # Set minimum n for each l
    for n in range(len(nmin)):
        if nmin[n] < n:
            nmin[n] = n + 1
        
    # Set maximum n for each l (4 lowest s, p, d; 3 lowest f) by default
    if val_aov == []:
        nmax = [n + 4 for n in nmin]
        nmax[3] -= 1
    else:
        for n in range(len(val_aov)):
            nmax[n] = nmin[n] + int(list(val_aov[n].values())[0])

    # Set kappas for valence orbitals
    for n in range(len(nmax)):
        if n == 0:
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append('-1')
        elif n == 1:
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append(' 1')
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append('-2')
        elif n == 2:
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append(' 2')
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append('-3')
        elif n == 3:
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append(' 3')
            for i in range(nmin[n],nmax[n]):
                N.append(i+1)
                kappa.append('-4')

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

def write_hfd_inp(filename, system, NS, NSO, Z, AM, kbr, NL, J, QQ, KP, NC, rnuc):
    """ Writes HFD.INP """
    # Define default values
    KL = 0
    JM = -2.0

    with open(filename, 'w') as f:
        f.write(' ' + system['name'] + '\n')
        f.write(' KL =  ' + str(KL) + '\n')
        f.write(' NS =  ' + str(NS) + '\n')
        f.write(' NSO=  ' + str(NSO) + '\n')
        f.write(' Z  =  ' + str(Z) + '\n')
        f.write(' AM =  ' + '{:.2f}'.format(AM) + '\n')
        f.write(' JM =  ' + str(JM) + '\n')
        f.write(' R2 =  ' + str(float(system['radius'])) + '\n')
        f.write(' kbr= ' + str(kbr) + '\n')
        f.write('rnuc= ' + '{:.4f}'.format(rnuc) + '\n\n')
        f.write('        NL   J       QQ     KP   NC\n\n')
        for i in range(len(NL)):
            f.write(str(i+1).rjust(3," ") + '     ' + NL[i] + ' (' + J[i] + ')' + '   ' 
                        + QQ[i] + '    ' + KP[i] + '   ' + NC[i].rjust(2,' ') + '\n')
    f.close()
    print('HFD.INP has been written')

def construct_vvorbs(core, valence, codename, nmax, lmax):
# Construct list of valence and virtual orbitals
    vorbs = []
    norbs = []

    ## First we need to construct list of orbitals from nmax + lmax
    nmin = [0, 0, 0, 0, 0, 0, 0]
    core_val = core + ' ' + valence
    l = ['s','p','d','f','g','h','i']
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


def write_bass_inp(filename, system, NSO, Z, AM, kbr, vorbs, norbs, nmax, lmax, codename, core, valence):
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
        f.write(' ' + system['name'] + '\n')
        f.write(' Z  =  ' + str(Z) + '\n')
        f.write(' Am =  ' + '{:.1f}'.format(AM) + '\n')
        f.write(' Nso=' + str(NSO).rjust(5," ") + '# number of core orbitals (defines DF operator)\n')
        f.write(' Nv =' + str(Nv).rjust(5," ") + '# number of valence & virtual orbitals\n')
        f.write(' Ksg=' + str(Ksg).rjust(5," ") + '# defines Hamiltonian: 1-DF, 3-DF+Breit\n')
        f.write(' Kdg=' + str(Kdg).rjust(5," ") + '# diagonalization of Hamiltonian (0=no,1,2=yes)\n')
        f.write(' orb=' + fvalorb.rjust(5," ") + '# first orbital for diagonalization\n')
        f.write(' Kkin' + str(Kkin).rjust(5," ") + '# kinetic balance (0,1,or 2)\n')
        f.write(' orb=' + fvirorb.rjust(5," ") + '# first orbital to apply kin.bal.\n')
        f.write(' orb=' + frorb.rjust(5," ") + '# last frozen orbital\n')
        f.write('kout=' + str(0).rjust(5," ") + '# detail rate in the output\n')
        f.write('kbrt=' + str(kbr).rjust(5," ") + '# 0,1,2 - Coulomb, Gaunt, Breit\n')
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
            f.write('\n')

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

def write_bas_wj_in(filename, symbol, Z, AM, NS, NSO, N, kappa, iters, energies, cfermi):
    """ Writes bas_wj.in """
    with open(filename,'w') as f: 
        f.write(' ' + symbol + ' ' + str(NS).rjust(4, ' ') + str(int(Z)).rjust(4, ' ') 
                    + str(round(AM)).rjust(4,' ') + '   0   0   9' + str(NSO+1).rjust(4, ' ')
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
            for n in range(len(energies)):
                energy = list(energies[n].values())[0]
                if isinstance(energy, list):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy[0]) + '  ' + "{:.5f}".format(energy[1]) + '\n')
                elif isinstance(energy, float):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy) + '\n')
        f.write(str(len(val_kappa)) + '\n')
        for i in range(len(val_kappa)):
            f.write(str(val_N[i]).rjust(2, ' ') + '  ' + val_kappa[i] + ' 30\n' )
    f.close()
    print('inf.aov has been written')
        

def write_inf_vw(filename, val_N, val_kappa, NSO, nmax, lmax, kvw, kval, energies):
    """ Writes inf.vw """
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
            for n in range(len(energies)):
                energy = list(energies[n].values())[0]
                if isinstance(energy, list):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy[0]) + '  ' + "{:.5f}".format(energy[1]) + '\n')
                elif isinstance(energy, float):
                    f.write(str(n).rjust(2,' ') + '  ' + "{:.5f}".format(energy) + '\n')
    f.close()
    print('inf.vw has been written')

def write_spl_in(filename, radius):
    """ Writes spl.in """
    with open(filename,'w') as f: 
        f.write('6\n')
        f.write(str(radius) + '\n')
        f.write('40 7\n')
        f.write('0.0 0.00 500')
    print('spl.in has been written')

if __name__ == "__main__":
    # Read yaml file for system configurations
    system = read_yaml('basis_Sr.yml')

    # Get atomic data
    Z, AM, symbol, cfermi, rnuc = get_atomic_data(system['name'], system['isotope'])

    # Get orbital information
    NS, NSO, num_val_orbitals = count_total_orbitals(system['core'], system['valence'])

    # Get breit key
    kbrt = get_key_breit(system['include_breit'])

    # Generate body of HFD.INP including orbitals, values of J, and occupation numbers
    NL, J, QQ, KP, NC, num_core_electrons, nval = gen_lists_orbitals(system['core'], system['valence'])
    
    # Assign keys for which inputs to generate
    kvw = get_key_vw(system['codes'])

    # Write HFD.INP
    write_hfd_inp('HFD.INP', system, NS, NSO, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc)
    
    # Write BASS.INP
    vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(system['core'], system['valence'], system['codes'], system['nmax'], system['lmax'])
    write_bass_inp('BASS.INP', system, NSO, Z, AM, kbrt, vorbs, norbs, system['nmax'], system['lmax'], system['codes'], system['core'], system['valence'])

    # Generate body of bas_wj.in including orbitals, values of kappa, and energy guesses
    N, kappa, iters, energies = gen_lists_kappa(Z, num_core_electrons, system['core'], system['valence'])

    # Write bas_wj.in
    write_bas_wj_in('bas_wj.in', symbol, Z, AM, NS, NSO, N, kappa, iters, energies, cfermi)

    # Write spl.in
    write_spl_in('spl.in', system['radius'])

    # Get valence orbitals for all-order calculations
    val_N, val_kappa = get_ao_valence(system['core'], system['valence'], system['val_aov'])

    # Set key for energies
    kval = system['kval']

    # Write inf.aov
    write_inf_aov('inf.aov', val_N, val_kappa, NSO, system['nmax'], system['lmax'], kval, system['energies'])

    # Write inf.vw
    write_inf_vw('inf.vw', val_N, val_kappa, NSO, system['nmax'], system['lmax'], kvw, kval, system['energies'])

    # Construct basis set by running sequence of programs
    ## Pure CI
    if system['codes'] == 'ci':
        # Run hfd
        run(["./hfd", ">", "hfd.out"])
        # Run bass
        run(["./bass", ">", "bass.out"])
    elif system['codes'] == 'ci+all-order':
        # Run hfd
        run('./hfd > hfd.out', shell=True)
        print("hfd complete")

        # Produce B-splines
        if kbrt == 0:
            run('./tdhf < bas_wj.in > tdhf.out', shell=True)
            print("tdhf complete")
            run('./nspl40 < spl.in > nspl40.out', shell=True)
            print("nspl40 complete")
        else:
            run('./bdhf < bas_wj.in > bdhf.out', shell=True)
            print("bdhf complete")
            run('./bspl40 < spl.in > bspl40.out', shell=True)
            print("bspl40 complete")

        with open('bwj.in','w') as f: 
            f.write(str(system['lmax']) + '\n')
            for i in range(system['lmax']+1):
                f.write(str(system['nmax']) + '\n')
            f.write('\n')
            f.write('\n')
            f.write('1')
        f.close()

        run('./bas_wj < bwj.in > bas_wj.out', shell=True)
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
        frorb_str = liborb.convert_digital_to_char(liborb.convert_char_to_digital(system['core'].split()[-1])[-1])
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
            elif line[:3].strip().isdigit() and int(line[:3].strip()) <= nval:
                f.write(line[:12] + "\n")
            else:
                f.write(line)
        f.close()

        # run bass
        with open('bass.in','w') as f: 
            f.write('WJ.DAT')
        f.close()
        run('./bass < bass.in > bass.out', shell=True)
        run(['rm', 'bass.in'])
        run(['rm','hfspl.1','hfspl.2'])
        print("bass complete")

        # Run bas_x
        run('./bas_x > bas_x.out', shell=True)
        print("bas_x complete")

        print("CI+all-order basis set construction complete")
    