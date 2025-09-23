""" Basis

This script allows the user to automate the basis set construction from inputted parameters in a "config.yml" file.
The "config.yml" file should have the following blocks:

    * system - general parameters (name of atomic system, isotope number, inclusion of breit)
    * basis - parameters used by basis programs (cavity radius, core and valence orbitals, b-splines)
    * optional - optional parameters (isotope shifts, code methods, running all-order codes, pci versions)

From these parameters, this script will create all input files required for the various basis codes.
After the input files are created, the sequence of basis set codes will be executed if the parameter run_codes is set to "True".

This python script has 2 main capabilities for basis set construction:
1. Construction of basis for isotope shift calculations 
2. Construction of basis for multiple code methods

"""
import yaml 
import re
import sys
import get_atomic_data as libatomic
import orbitals as liborb
import os
from pathlib import Path
from utils import run_shell, get_dict_value, check_slurm_installed
from gen_job_script import write_job_script

def read_yaml(filename):
    """ Reads a configuration file in YAML format and returns a dictionary of config parameters """ 
    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def validate_config(config):
    """ Validates the configuration file for required fields """
    required_fields = ['system', 'atom', 'basis', 'optional']
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required field '{field}' in configuration file.")

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
        f.write(' ' + system['atom']['name'] + '\n')
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
        if rnuc:
            f.write('rnuc= ' + '{:.4f}'.format(rnuc) + '\n\n')
        else:
            f.write('\n')
        f.write('        NL   J       QQ     KP   NC\n\n')
        for i in range(len(NL)):
            f.write(str(i+1).rjust(3," ") + '     ' + NL[i] + ' (' + J[i] + ')' + '   ' 
                        + QQ[i] + '    ' + str(KP[i]) + '   ' + str(NC[i]).rjust(2,' ') + '\n')
    f.close()
    print(filename + ' has been written')

def write_hfd_inp_ci(filename, system, num_electrons, Z, AM, kbrt, NL_base, J_base, QQ_base, KP_base, NC_base, rnuc, K_is, C_is):
    """ Write multiple HFD.INP files for case of pure CI"""

    basis = system['basis']
    cavity_radius = basis['cavity_radius']
    core_orbitals = basis['orbitals']['core']
    valence_orbitals = basis['orbitals']['valence']
    order = [shell.strip() for shell in basis['orbitals']['order'].split('/')]
    
    num_hfd_inps = len(order)
    
    # Obtain base list of shells
    NL_base, J_base, QQ_base, KP_base, NC_base, num_core_electrons, nval = gen_lists_orbitals(core_orbitals, valence_orbitals)
    
    index = 1
    base_hfd_index = 1
    shells_found = {}
    frozen_found = { }
    frozen_shells = []
    
    core_shells = []        
    for orbital in core_orbitals.split(" "):
        core_shells.append(orbital[0] + orbital[1].capitalize())
        
    if num_core_electrons == num_electrons:
        core_shells.pop()
    
    # Loop through list of shells to form orbitals from
    for shell_list in order:
        NL, J, QQ, KP, NC = [], [], [], [], []

        shells = shell_list.split(" ")
        
        # Count the number of electrons when forming orbitals
        electron_cnt = 0
        electron_removed = False
        
        # Reset shells found
        shells_found = shells_found.fromkeys(shells_found, False)
        frozen_found = frozen_found.fromkeys(frozen_found, False)
        
        # Add new shells to shells_found
        for shell in shells:
            # Get the shell in the same format as NL
            shell_fmt = shell[0] + shell[1].capitalize()        
            shells_found[shell_fmt] = False
        
        nc_it = 0
        
        # Loop through base list of shells formed from all core and valence orbitals
        for i in range(len(NL_base)):
            NL.append(NL_base[i])
            J.append(J_base[i])
            
            QQ_num = int(J_base[i].split('/')[0]) + 1
            
            electron_cnt += QQ_num
            
            # If NL is in the list of found shells, set shells_found to 'True'
            if NL[i] in list(shells_found.keys()): 
                shells_found[NL[i]] = True
                
            # If NL is in the list of frozen shells, set frozen_found to 'True
            if NL[i] in list(frozen_found.keys()):
                frozen_found[NL[i]] = True
            
            # Exit the loop if the following conditions are met:
            # 1. all shells in the "order" list have been found
            # 2. NL is not in the list of "order" shells
            all_shells_found = all(found == True for found in shells_found.values())
            if all_shells_found and NL[i] not in shells_found.keys():
                NL.pop()
                J.pop()
                break
            
            # Finding when to remove an electron distinguishes core and valence shells  
            # If an electron hasn't been removed yet, check if the electron count has surpassed the N-1 threshold
            if not electron_removed:
                # If it has, remove an electron from QQ
                if electron_cnt > num_electrons - 1:
                    QQ.append(f"{QQ_num - 1:.4f}")
                    electron_removed = True
                else:
                    QQ.append(f"{QQ_num:.4f}")
                    
            # If an electron has already been removed, set QQ to 0.0000 if any of the following conditions match:
            # 1. if the current shell is frozen
            # 2. current shell is the same shell as last iteration
            # 3. current shell is not in the shells_found list
            else:
                if NL[i] in frozen_shells or NL[i] == NL[i-1] or NL[i] not in shells_found.keys():
                    QQ.append(f"{0:.4f}")
                else:
                    QQ.append(f"{1:.4f}")
            
            # By default, first list has KP = 0 and NC = 0 for all shells
            if index == 1:
                KP.append(0)
                NC.append(0)
            else:                
                # Set KP = 1 if the current shell is in core or is frozen, else 0
                if NL[i] in core_shells or NL[i] in frozen_found.keys():
                    KP.append(1)
                else:
                    KP.append(0)
                
                # Set NC = 0 if the current shell is in core or is frozen, else iterate per shell
                if NL[i] in core_shells or NL[i] in frozen_found.keys() and NC[i-1] == 0:
                    NC.append(0)
                else:
                    if NL[i] != NL[i-1]:
                        nc_it += 1
                    NC.append(nc_it)
        
        # Freeze the shells for next iteration
        for shell in shells:
            shell_fmt = shell[0] + shell[1].capitalize() 
            if shell_fmt not in frozen_shells:
                frozen_shells.append(shell_fmt)
                frozen_found[shell_fmt] = False

        # Write the HFD.INP to construct orbitals for this shell
        NS = len(NL)
        NSO = count_orbitals_in_list(core_shells)
        write_hfd_inp('HFD'+str(index)+'.INP', system, NS, NSO, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, K_is, C_is)
        
        index += 1

def count_orbitals_in_list(shell_list):
    # count number of orbitals in list of shells
    num_orbitals = 0
    for shell in shell_list:
        if 'S' in shell or 's' in shell:
            num_orbitals += 1
        else:
            num_orbitals += 2
    
    return num_orbitals

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


def write_bass_inp(filename, system, NSO, Z, AM, kbr, vorbs, norbs, K_is, C_is):
    """ Writes BASS.INP """
    # Define default values
    basis = system['basis']
    nmax = basis['orbitals']['nmax']
    lmax = basis['orbitals']['lmax']
    codename = atom['code_method']
    method = basis['method']
    core = basis['orbitals']['core']
    valence = basis['orbitals']['valence']
    
    Nv = len(valence.split())
    Ksg = 1
    
    if system['basis']['diagonalized']:
        Kdg = 1
    else:
        Kdg = 0
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

    # Handle custom orbitals
    custom_orbs = {}
    custom_vorbs = {}
    try:
        custom = system['basis']['orbitals']['custom']
    except KeyError:
        custom = ""
    
    if custom:
        for orbital in custom:
            orb = orbital.split(' ')[0]
            from_orb = orbital.split(' ')[-1]
            custom_orbs[orb] = from_orb
            for vorb in liborb.convert_char_to_digital(orb):
                vorb_fmt = vorb[:-2] + '01'
                orb_fmt = orb[:-2] + '01'
                from_orb_fmt = from_orb[:-2] + '01'
                # If virtual orbital constructed from hfd, set value same as key
                if from_orb == 'hfd':
                    custom_vorbs[vorb_fmt] = vorb_fmt
                else:
                    for from_orb_fmt in liborb.convert_char_to_digital(from_orb):
                        if vorb_fmt[0] == from_orb_fmt[0]:
                            from_orb_fmt = from_orb_fmt[:-2] + '01'
                            custom_vorbs[vorb_fmt] = from_orb_fmt    
    
    with open(filename, 'w') as f:
        f.write(' ' + system['atom']['name'] + '\n')
        f.write(' Z  = ' + str(Z) + '\n')
        f.write(' Am = ' + '{:.1f}'.format(round(AM)) + '\n')
        f.write(' Nso=' + str(NSO).rjust(5," ") + ' # number of core orbitals (defines DF operator)\n')
        f.write(' Nv =' + str(nvvorbs).rjust(5," ") + ' # number of valence & virtual orbitals\n')
        f.write(' Ksg=' + str(Ksg).rjust(5," ") + ' # defines Hamiltonian: 1-DF, 3-DF+Breit\n')
        f.write(' Kdg=' + str(Kdg).rjust(5," ") + ' # diagonalization of Hamiltonian (0=no,1,2=yes)\n')
        f.write(' orb=' + fvalorb.rjust(5," ") + ' # first orbital for diagonalization\n')
        f.write(' Kkin' + str(Kkin).rjust(5," ") + ' # kinetic balance (0,1,or 2)\n')
        f.write(' orb=' + fvirorb.rjust(5," ") + ' # first orbital to apply kin.bal.\n')
        f.write(' orb=' + frorb.rjust(5," ") + ' # last frozen orbital\n')
        f.write('kout=' + str(0).rjust(5," ") + ' # detail rate in the output\n')
        f.write('kbrt=' + str(kbr).rjust(5," ") + ' # 0,1,2 - Coulomb, Gaunt, Breit\n')
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
            if (i < nval) or method == 'dirac-fock':
                f.write(str(i+1).rjust(3," ") + orb.rjust(8," "))
                # Check if orbital is in the list of custom virtual orbitals
                if orb in list(custom_vorbs.keys()):
                    # Check if key and value for the custom orbital is the same 
                    if orb == custom_vorbs[orb]:
                        f.write("  3 ")
                    else:
                        f.write("    ")
                    f.write(custom_vorbs[orb].rjust(7," ") + "\n")
                else:
                    f.write("\n")
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
        grid = '0.00004'
        f.write(' ' + str(grid).rjust(9, ' ') + '  0.00  500\n')
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

def write_ao_inputs(system, K_is, C_is, kvw, basis_method):
    if basis_method == 'b-splines':
        write_hfd_inp('HFD.INP', system, NS, NSO, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, K_is, C_is)
        vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(core_orbitals, valence_orbitals, code_method, basis_nmax, basis_lmax)
        write_bass_inp('BASS.INP', system, NSO, Z, AM, kbrt, vorbs, norbs, K_is, C_is)
        write_bas_wj_in('bas_wj.in', symbol, Z, AM, NS, NSO, N, kappa, iters, energies, cfermi)
        write_spl_in('spl.in', system['basis']['cavity_radius'], system['basis']['b_splines'])
    elif basis_method == 'dirac-fock':
        write_hfd_inp_ci('HFD.INP', config, num_electrons, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, K_is, C_is)
        vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(core_orbitals, valence_orbitals, code_method, basis_nmax, basis_lmax)
        write_bass_inp('BASS.INP', config, NSO, Z, AM, kbrt, vorbs, norbs, K_is, C_is)
    else:
        raise ValueError('basis_method not valid. Expecting "b-splines" or "dirac-fock".')
    write_inf_aov('inf.aov', val_N, val_kappa, NSO, system['basis']['orbitals']['nmax'], system['basis']['orbitals']['lmax'], kval, system['basis']['val_energies']['energies'])
    write_inf_vw('inf.vw', val_N, val_kappa, NSO, system['basis']['orbitals']['nmax'], system['basis']['orbitals']['lmax'], kvw, kval, system['basis']['val_energies']['energies'])

def generate_batch_qed(bin_dir, kqed, kbrt):
    """ Writes batch.qed """
    
    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
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
        f.write(bin_dir + "qedpot_conf <q.in >qp.res \n")
        f.write(bin_dir + "qed_rot <q.in >qr.res \n")
        f.write("grep 'changed' \"QED_ROT.RES\" \n")
        f.write("  if grep -q reached \"QED_ROT.RES\"; then \n")
        f.write("  echo 'Converged in '$n' iterations' \n")
        f.write("  break \n")
        f.write("  fi \n")
        f.write("  let n=n+1 \n")
        f.write("done \n")
        f.write("##################################### \n")
    print('batch.qed has been written')
    run_shell('chmod +x batch.qed')
    
def check_errors(filename):
    # This function checks output files for errors and returns the number of errors
    # Currently only supports output of program bass
    num_errors = 0
    
    if filename == 'bass.out':
        warning_msgs = ['failed', 'error', 'warning']
    elif filename == 'HFD.RES':
        warning_msgs = ['NaN']
    else:
        print(filename + "not currently supported")

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

def run_ci_executables(bin_dir, order, custom):
    # Remove old HFD.INP
    if os.path.isfile('HFD.INP'):
        run_shell('rm HFD.INP')
    file_list = os.listdir(".")
    hfd_list = []
    for file in file_list:
        if file[:3] == 'HFD' and file[-3:] == 'INP':
            hfd_list.append(file)
    
    hfd_list.sort()
            
    # Run hfd for HFD.INP files
    print('Found the following HFD.INP files:', ', '.join(hfd_list))
    for file in hfd_list:
        run_shell('cp ' + file + ' HFD.INP')
        run_shell(bin_dir + 'hfd > hfd' + file[-5] + '.out')
        run_shell('cp HFD.DAT ' + file[:-3] + 'DAT')
        run_shell('cp HFD.RES ' + file[:-3] + 'RES')
        print('hfd completed with ' + file)
        
    # Clean up
    run_shell('rm HFD.INP HFD.RES')
    
    # Find base HFD.DAT to construct basis set
    # Check if key and value is the same in custom_vorbs (from hfd)
    # We assume corresponding HFD.DAT was constructed to make orbitals, and use previous as the base
    from_hfd = []
    for shell_origin in custom:
        shell = shell_origin.split(' ')[0]
        origin = shell_origin.split(' ')[-1]
        if origin == 'hfd':
            from_hfd.append(shell)

    # Figure out which HFD.INP created those orbitals - HFD.DAT should correspond to the latest HFD.INP that does not contain those orbitals
    order_list = order.split('/')
    for base_hfd_index in range(len(order_list)):
        order_list_fmt = order_list[base_hfd_index].strip()
        shells = order_list_fmt.split(' ')
        if not any(shell in from_hfd for shell in shells):
            base_hfd_index += 1
    
    # Prepare inputs for bass
    run_shell('cp HFD' + str(base_hfd_index) + '.DAT HFD.DAT')
    with open('bass.in', 'w') as f:
        f.write('HFD' + str(base_hfd_index + 1) + '.DAT')
    
    # Run bass
    # check if bass.out exists and remove if it does
    if os.path.isfile('bass.out'):
        run_shell('rm bass.out')
        
    maxNumTries = 5
    nTry = 1

    while check_errors('bass.out') > 0:
        print('bass attempt', nTry)
        run_shell(bin_dir + 'bass < bass.in > bass.out')
        
        run_shell('cp bass.out ' + 'bass' + str(nTry) + '.out')
            
        if (nTry >= maxNumTries):
            print("bass did not converge after", nTry, "attempts")
            break
        nTry += 1
        
    else:
        print("bass completed with no errors after", nTry - 1, "attempts")
    
def run_ao_executables(diag_basis, K_is, C_is, bin_dir, order, custom, basis_method):
    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
    # Produce B-splines
    if basis_method == 'b-splines':
        # Run hfd
        run_shell(bin_dir + 'hfd > hfd.out')
        print("hfd complete")

        if check_errors('HFD.RES'):
            print("Error found in HFD.RES. Please check.")
            sys.exit()

        if kbrt == 0:
            run_shell(bin_dir + 'tdhf < bas_wj.in > tdhf.out')
            print("tdhf complete")
            run_shell(bin_dir + 'nspl < spl.in > nspl.out')
            print("nspl complete")
        else:
            run_shell(bin_dir + 'bdhf < bas_wj.in > bdhf.out')
            print("bdhf complete")
            run_shell(bin_dir + 'bspl < spl.in > bspl.out')
            print("bspl complete")
    
        with open('bwj.in','w') as f: 
            f.write(str(basis_lmax) + '\n')
            for i in range(basis_lmax+1):
                f.write(str(basis_nmax) + '\n')
            f.write('\n')
            f.write('\n')
            f.write('1')
    
        run_shell(bin_dir + 'bas_wj < bwj.in > bas_wj.out')
        run_shell('rm bwj.in')
        print("bas_wj complete")
    
        # Edit BASS.INP
        with open('BASS.INP', 'r') as f:
            lines = f.readlines()
    
        # Set diagonalization
        Kdg = '0' if not diag_basis else '1'
    
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
        
        with open('BASS.INP','w') as f:
            for line in lines:
                if line[:4] == ' lst':
                    continue
                elif 'diagonalization' in line.split() and 'Hamiltonian' in line.split():
                    f.write(' Kdg=' + Kdg.rjust(5," ") + '# diagonalization of Hamiltonian (0=no,1,2=yes)\n')
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
    
        # run bass until there are no errors in output
        with open('bass.in','w') as f: 
            f.write('WJ.DAT')

        # check if bass.out exists and remove if it does
        if os.path.isfile('bass.out'):
            run_shell('rm bass.out')
    
        maxNumTries = 5
        nTry = 1
        needs_hfd_dat = False
        for orb_build in custom:
            if 'from hfd' in orb_build:
                needs_hfd_dat = True
    
        while check_errors('bass.out') > 0:
            print('bass attempt', nTry)
            if needs_hfd_dat:
                run_shell(bin_dir + 'bass < bass.in > bass.out')
            else:
                run_shell(bin_dir + 'bass > bass.out')
            if (nTry >= maxNumTries):
                print("bass did not converge after", nTry, "attempts")
                break
            nTry += 1
        else:
            print("bass completed with no errors after", nTry, "attempts")
    elif basis_method == 'dirac-fock':
        run_ci_executables(bin_dir, order, custom)

    if os.path.isfile('bass.in'):
        run_shell('rm bass.in')
    if os.path.isfile('hfspl.1'):
        run_shell('rm hfspl.1 hfspl.2')

    # Run bas_x
    run_shell(bin_dir + 'bas_x > bas_x.out')
    print("bas_x complete")

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    validate_config(config)
    
    system = get_dict_value(config, 'system')
    atom = get_dict_value(config, 'atom')
    basis = get_dict_value(config, 'basis')
    optional = get_dict_value(config, 'optional')
    
    # system parameters
    bin_dir = get_dict_value(system, 'bin_directory')
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
        
    on_hpc = get_dict_value(system, 'on_hpc')
    run_codes = get_dict_value(system, 'run_codes')
    pci_version = get_dict_value(system, 'pci_version')
    
    # hpc parameters
    on_slurm = check_slurm_installed()
    if on_hpc and on_slurm:
        hpc = get_dict_value(config, 'hpc')
        submit_job = get_dict_value(hpc, 'submit_job')
        if hpc:
            partition = get_dict_value(hpc, 'partition')
            nodes = get_dict_value(hpc, 'nodes')
            tasks_per_node = get_dict_value(hpc, 'tasks_per_node')
        else:
            print('hpc block was not found in', yml_file)
            partition, nodes, tasks_per_node = None, 1, 1
    else:
        on_hpc = False
        submit_job = False
        
    # atom parameters
    name = atom['name']
    try:
        isotope = atom['isotope']
    except KeyError:
        isotope = ""
    include_breit = atom['include_breit']
    code_method = atom['code_method']
    
    # basis parameters
    orbitals = get_dict_value(basis, 'orbitals')
    core_orbitals = get_dict_value(orbitals, 'core')
    valence_orbitals = get_dict_value(orbitals, 'valence')
    diagonalize_basis = get_dict_value(basis, 'diagonalized')
    basis_nmax = get_dict_value(orbitals, 'nmax')
    basis_lmax = get_dict_value(orbitals, 'lmax')
    val_energies = get_dict_value(basis, 'val_energies')
    if val_energies: kval = get_dict_value(val_energies, 'kval')
    val_aov = get_dict_value(basis, 'val_aov')
    order = get_dict_value(orbitals, 'order')
    custom = get_dict_value(orbitals, 'custom')
    basis_method = get_dict_value(basis, 'method')

    # isotope shift parameters
    isotope_shifts = get_dict_value(optional, 'isotope_shifts')
    include_isotope_shifts = get_dict_value(isotope_shifts, 'include')
    if include_isotope_shifts:
        K_is = get_dict_value(isotope_shifts, 'K_is')
        C_is = get_dict_value(isotope_shifts, 'C_is')
        c_list = [-C_is,-C_is/2,0,C_is/2,C_is]
        K_is_dict = {0: '', 1: 'FS', 2: 'SMS', 3: 'NMS', 4: 'MS'}
    
    # qed parameters
    qed = get_dict_value(optional, 'qed')
    include_qed = get_dict_value(qed, 'include')

    # Get atomic data
    Z, AM, symbol, cfermi, rnuc, num_electrons = libatomic.get_atomic_data(name, isotope)

    # Get orbital information
    NS, NSO, num_val_orbitals = count_total_orbitals(core_orbitals, valence_orbitals)

    # Get breit key
    kbrt = get_key_breit(include_breit)

    # Generate body of HFD.INP including orbitals, values of J, and occupation numbers
    NL, J, QQ, KP, NC, num_core_electrons, nval = gen_lists_orbitals(core_orbitals, valence_orbitals)
    
    ao_methods = ['all-order','second-order','ci+all-order','ci+second-order']
    ao_method_in_code_methods = any(x in code_method for x in ao_methods)

    if ao_method_in_code_methods:
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
                    run_shell('pwd')
                    write_ao_inputs(config, K_is, c, get_key_vw(method), basis_method)
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
                    run_shell('pwd')
                    write_ao_inputs(config, 0, 0, get_key_vw(method), basis_method)
                    os.chdir('../../')
            else:
                dir_path = os.getcwd()
                Path(dir_path+'/basis').mkdir(parents=True, exist_ok=True)
                os.chdir('basis')
                run_shell('pwd')
                write_ao_inputs(config, 0, 0, kvw, basis_method)
                os.chdir('../')

        # Construct basis set by running sequence of programs if desired
        if run_codes:
            if not on_hpc:
                print('run_codes option is only available with HPC access')
            else:
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
                            run_shell('pwd')
                            run_ao_executables(diagonalize_basis, K_is, c, bin_dir, order, custom, basis_method)
                            script_name = write_job_script('.', method, nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
                            if script_name and submit_job:
                                run_shell('sbatch ' + script_name)
                            else:
                                print('job script was not submitted. check job script and submit manually.')
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
                            run_shell('pwd')
                            run_ao_executables(diagonalize_basis, 0, 0, bin_dir, order, custom, basis_method)
                            script_name = write_job_script('.', method, nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
                            if script_name and submit_job:
                                run_shell('sbatch ' + script_name)
                            else:
                                print('job script was not submitted - check job script and submit manually.')
                            os.chdir('../../')
                    else:
                        dir_path = os.getcwd()
                        Path(dir_path+'/basis').mkdir(parents=True, exist_ok=True)
                        os.chdir('basis')
                        run_shell('pwd')
                        run_ao_executables(diagonalize_basis, 0, 0, bin_dir, order, custom, basis_method)
                        script_name = write_job_script('.', code_method, nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
                        if script_name and submit_job:
                            run_shell('sbatch ' + script_name)
                        else:
                            print('job script was not submitted - check job script and submit manually.')
                        os.chdir('../')
            
    elif code_method == 'ci':
        if include_isotope_shifts and K_is > 0:
            dir_path = os.getcwd()
            is_dir = K_is_dict[K_is]
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
                Path(dir_name).mkdir(parents=True, exist_ok=True)
                os.chdir(dir_name)
                print(dir_name)
                write_hfd_inp_ci('HFD.INP', config, num_electrons, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, K_is, c)
                vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(core_orbitals, valence_orbitals, code_method,   basis_nmax, basis_lmax)
                write_bass_inp('BASS.INP', config, NSO, Z, AM, kbrt, vorbs, norbs, K_is, c)
                
                if run_codes:
                    run_ci_executables(bin_dir, order, custom)
                    
                if K_is_dict[K_is]:
                    os.chdir('../../')
                else:
                    os.chdir('../')
                
        else:
            dir_path = os.getcwd()
            Path(dir_path+'/basis').mkdir(parents=True, exist_ok=True)
            os.chdir('basis')
            run_shell('pwd')
            write_hfd_inp_ci('HFD.INP', config, num_electrons, Z, AM, kbrt, NL, J, QQ, KP, NC, rnuc, 0, 0)
            vorbs, norbs, nvalb, nvvorbs = construct_vvorbs(core_orbitals, valence_orbitals, code_method, basis_nmax, basis_lmax)
            write_bass_inp('BASS.INP', config, NSO, Z, AM, kbrt, vorbs, norbs, 0, 0)
            
            if run_codes:
                run_ci_executables(bin_dir, order, custom)
                        
        if include_qed:
            generate_batch_qed(bin_dir, include_qed, kbrt)
            if run_codes: 
                run_shell(bin_dir+'sgc0')
                run_shell('./batch.qed > qed.out')
        
        os.chdir('../')

    else:
        print('code method', code_method, ' is not supported')