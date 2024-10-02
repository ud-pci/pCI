""" ci

This script allows the user to automate the configuration list construction from inputted parameters in a "config.yml" file. 
The "config.yml" file should have the following blocks:

    * system - system parameter (bin_directory, run_codes, on_hpc)
    * basis.orbitals.core - core orbitals are read to ensure consistent basis set
    * ci - parameters used by ci programs (reference configurations, basis set, orbitals, excitations)
    * optional - optional parameters (isotope shifts, code methods, running all-order codes, pci versions)

From these parameters, this script will create all input files required for execution of the add program.
After the input files are created, the add program will be executed to create the list of configurations CONF.INP.
If optional.generate_directories is set to "True", the script will generate respective directories for CI calculations (e.g. /even and /odd)
If optional.run_codes is set to "True", the script will then submit the job script in the respective directories. 

This python script has 3 main capabilities for configuration list construction:
1. General construction of configuration lists in root directory
2. Construction of configuration lists for isotope shift calculations 
3. Construction of configuration lists for multiple code methods

"""
import yaml
import os
import sys
import re
import math
import orbitals as orb_lib
import get_atomic_data as libatomic
from utils import run_shell, get_dict_value
from pathlib import Path
from gen_job_script import write_job_script

def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def write_add_inp(filename, Z, AM, config, multiplicity, num_val, orb_occ, parity):
    """Write ADD.INP file for specified parity"""

    configurations = config['add']['ref_configs']
    orbitals = config['add']['orbitals']
    # Check if there are configurations 
    if configurations[parity] == []:
        print(parity.capitalize(), 'parity configurations were not specified')
        return

    if parity == 'odd':
        parity_config = config['conf']['odd']
    elif parity == 'even':
        parity_config = config['conf']['even']
    else:
        print('parity ' + parity + ' is not valid')
        return
    
    # Set number of energy levels
    J = parity_config['J']
    JM = parity_config['JM']
    J_selection = parity_config['J_selection']
    num_energy_levels = parity_config['num_energy_levels']
    num_dvdsn_iterations = parity_config['num_dvdsn_iterations']

    # Define name of file to export
    filename = filename[0:-4] + parity + filename[-4:]
    f = open(filename,'w')

    # Write header
    f.write('Ncor= {:d}\n'.format(len(configurations[parity])))
    f.write('NsvNR{:d}\n'.format(orb_occ['num_orb']))
    f.write('mult= {:d}\n'.format(multiplicity))
    f.write(' NE = {:d}\n\n'.format(num_val))

    # Reformat configurations
    configurations_formatted = format_configurations(configurations[parity])
    max_orb_per_line = 6

    # Write list of basic configurations
    for configuration in configurations_formatted:
        num_lines = math.ceil(len(configuration)/max_orb_per_line)
        for line in range(num_lines):
            f.write('L:  ' + ' '.join(configuration[line*max_orb_per_line:(line+1)*max_orb_per_line]) + '\n')
    f.write('\n')

    # Write list of orbitals and occupation numbers
    nmin = orb_occ['nmin']
    nmax = orb_occ['nmax']
    count = 0

    for orb in orb_occ:
        orb_occ_formatted = format_orb_occ(orb, orb_occ[orb])
        f.write('  ' + orb_occ_formatted)
        count += 1
        if count == 6:
            f.write('\n')
            count = 0

    f.write('\n')

    num_core_orb = 0
    if orb_occ['core'] != 0:
        try:
            core_formatted = orb_lib.convert_char_to_digital(orb_occ['core'])
            num_core_orb = len(core_formatted)
        except Exception as e:
            print(e)
            pass
    else:
        core_formatted = []
        
    atom = config['atom']
    conf = config['conf']
    # Write head of CONF.INP
    f.write('>>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>\n')
    f.write('  ' + atom['name'] + ' ' + parity + '\n')
    f.write('  Z = ' + str(Z) + '\n')
    f.write(' Am = ' + "{:.1f}".format(round(AM)) + '\n')
    f.write('  J = ' + str(J) + '\n')
    f.write(' Jm = ' + str(JM) + '\n')
    f.write(' Nso=  ' + str(num_core_orb) + '\n')
    f.write(' Nc =   10 \n')
    if J_selection:
        f.write(' Kv =  3 \n')
    else:
        f.write(' Kv =  4 \n')
    f.write(' Nlv=  ' + str(num_energy_levels) + '\n')
    f.write(' Ne =  ' + str(num_val) + '\n')
    f.write(' Kl4=  1 \n')
    f.write(' Nc4=999 \n')
    f.write(' Gj = 0.0000 \n')
    f.write('Crt4= 0.0001 \n')
    f.write('kout= 0 \n')
    f.write('Ncpt= 0 \n')
    f.write('Cut0= 0.0001 \n')
    f.write('N_it= ' + str(num_dvdsn_iterations) + '\n')
    if atom['include_breit']:
        kbrt = 2
    else:
        kbrt = 0
    f.write('Kbrt= ' + str(kbrt) + '\n')
    
    # Write core orbitals
    if num_core_orb != 0:
        count = 0
        for i in range(0, num_core_orb):
            orb = core_formatted[i]
            if orb[0] == '-':
                f.write('    ' + orb)
            else:
                f.write('     ' + orb)
            count += 1
            if count == 6:
                f.write('\n')
                count = 0
    
    f.write('\n    ')

    print(filename + ' has been written')
    f.close()

def format_configurations(configurations):
    """ Format configurations as 'nnlqq',
        where nn represents principal quantum number
               l represents angular quantum number
              qq represents occupation number """
    configurations_formatted = []
    for configuration in configurations:
        orbitals_formatted = []
        for orbital in configuration.split():
            num_digits = len(re.findall('[0-9]+', orbital)[0])
            if num_digits == 1:
                orbital_formatted = orbital.center(5, ' ')
            elif num_digits == 2:
                orbital_formatted = orbital.ljust(5, ' ')
            else:
                print('ERROR: orbital ' + orbital + ' exceeds n=99')     
                sys.exit()         
            orbitals_formatted.append(orbital_formatted)
        configurations_formatted.append(orbitals_formatted)
    return configurations_formatted

def format_orb_occ(orb, occ):
    """ Format configurations as 'nnl qq qq',
        where nn represents principal quantum number
               l represents angular quantum number
              qq represents occupation number """
    ignore_list = ['num_orb', 'core', 'nmin', 'nmax']
    if orb in ignore_list:
        orb_occ_formatted = ''
    else:
        orb_occ_formatted = orb.rjust(3, ' ') + occ[0].rjust(3, ' ') + occ[1].rjust(3, ' ')

    return orb_occ_formatted

def create_add_inp(config):
    # Get atomic data
    atom = config['atom']
    
    try:
        Z, AM, symbol, cfermi, rnuc, num_rem_ele = libatomic.get_atomic_data(atom['name'], atom['isotope'])
    except KeyError as e:
        Z, AM, symbol, cfermi, rnuc, num_rem_ele = libatomic.get_atomic_data(atom['name'], "")
    
    add = config['add']
    num_val = orb_lib.count_valence(add['ref_configs'])
    orb_occ = orb_lib.expand_orbitals(add['basis_set'], add['ref_configs'], add['orbitals'])
    multiplicity = orb_lib.count_excitations(add['excitations'])

    if add['ref_configs']['even']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'even')
    else:
        print('no even reference configurations specified')
        
    if add['ref_configs']['odd']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'odd')
    else:
        print('no odd reference configurations specified')

def form_conf_inp(parity, bin_dir):
    if bin_dir and bin_dir[-1] == '/':
        bin_dir = bin_dir[:-1]
    run_shell("cp ADD" + parity + ".INP ADD.INP")
    if on_hpc:
        run_shell("add < add.in > add" + parity + ".out")
    else:
        run_shell(bin_dir + "/add < add.in > add" + parity + ".out")
    print("output of add saved to add" + parity + ".out")
    run_shell("cp CONF.INP CONF" + parity + ".INP")
    print("CONF" + parity + ".INP created")

def move_conf_inp(root_dir, parity, run_codes, include_lsj, write_hij):
    if not os.path.isdir(parity):
        run_shell("mkdir " + parity)
    if os.path.isfile('basis/HFD.DAT'):
        run_shell("cp basis/HFD.DAT " + parity)
    if os.path.isfile(root_dir + '/CONF' + parity + '.INP'):
        run_shell("cp " + root_dir + "/CONF" + parity + ".INP " + parity + "/CONF.INP")
    if os.path.isfile(root_dir + '/ADD' + parity + '.INP'):
        run_shell("cp " + root_dir + "/ADD" + parity + ".INP " + parity + "/ADD.INP")
    if os.path.isfile('basis/SGC.CON'):
        run_shell("cp basis/SGC.CON "  + parity)
    if os.path.isfile('basis/SCRC.CON'):
        run_shell("cp basis/SCRC.CON "  + parity)
    if run_codes and os.path.isfile(root_dir + '/ci.qs'):
        run_shell("cp " + root_dir + "/ci.qs " + parity)
        
    Kw = '1' if write_hij else '0'
    KLSJ = '1' if include_lsj else '0'

    if os.path.isfile('basis/SGC.CON'):
        with open(parity + '/ci.in', 'w') as f:
            f.write('Kl = 2 \n')
            f.write('Ksig = 2 \n')
            f.write('Kdsig = 0 \n')
            f.write('Kw = ' + Kw + '\n')
            f.write('KLSJ = ' + KLSJ)
    else:
        with open(parity + '/ci.in', 'w') as f:
            f.write('Kl = 0 \n')
            f.write('Ksig = 0 \n')
            f.write('Kdsig = 0 \n')
            f.write('Kw = ' + Kw + '\n')
            f.write('KLSJ = ' + KLSJ)

def check_conf_inp_exists(dir):
    if not os.path.isfile('HFD.DAT'):
        print('HFD.DAT is missing from ' + dir + ' directory')
        sys.exit()
    if not os.path.isfile('CONF.INP'):
        print('CONF.INP is missing from ' + dir + ' directory')
        sys.exit()
    if not os.path.isfile('ci.in'):
        print('ci.in is missing from ' + dir + ' directory')
        sys.exit()

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    
    include_isotope_shifts = config['optional']['isotope_shifts']['include']
    if include_isotope_shifts:
        K_is = config['optional']['isotope_shifts']['K_is']
        C_is = config['optional']['isotope_shifts']['C_is']
        c_list = [-C_is,-C_is/2,0,C_is/2,C_is]
        K_is_dict = {0: '', 1: 'FS', 2: 'SMS', 3: 'NMS', 4: 'MS'}
        
    code_method = config['atom']['code_method']
    run_codes = config['system']['run_codes']
    include_lsj = config['conf']['include_lsj']
    pci_version = config['system']['pci_version']
    write_hij = config['conf']['write_hij']
    bin_dir = config['system']['bin_directory']
    on_hpc = config['system']['on_hpc']
    
    # hpc parameters
    if on_hpc:
        hpc = get_dict_value(config, 'hpc')
        if hpc:
            partition = get_dict_value(hpc, 'partition')
            nodes = get_dict_value(hpc, 'nodes')
            tasks_per_node = get_dict_value(hpc, 'tasks_per_node')
        else:
            print('hpc block was not found in', yml_file)
            partition, nodes, tasks_per_node = None, 1, 1
    
    # Ensure basis and add core orbitals match
    basis_core = config['basis']['orbitals']['core']
    add_core = config['add']['orbitals']['core']
    
    # Check if user wants to generate directories for CI computations
    gen_dir = run_codes if run_codes else config['system']['generate_directories']
    
    # CONF.INP should only need to be constructed once, then copied to respective directories
    # Read input to add from add.in if it exists, otherwise create it
    try:
        open("add.in", 'rb').read()
    except:
        f = open("add.in", "w")
        f.write("0 \n 0")
        f.close()
    
    # BASS.INP will dictate order of orbitals in ADD.INP, so find a BASS.INP file from the basis directory
    bass_inp_paths = []
    if not os.path.isfile('BASS.INP'):
        for root, dirs, files in os.walk('.'):
            if 'BASS.INP' in files:
                bass_inp_paths.append(root)
    if code_method == 'ci':
        run_shell('cp ./basis/BASS.INP .')
        print('BASS.INP taken from basis/')
    else:
        if os.path.isdir('ci+all-order'):
            run_shell('cp ./ci+all-order/basis/BASS.INP .')
            print('BASS.INP taken from ci+all-order/basis/')
        elif os.path.isdir('ci+second-order'):
            run_shell('cp ./ci+second-order/basis/BASS.INP .')
            print('BASS.INP taken from ci+all-order/basis/')

    # Run script to generate ADD.INP for even and odd configurations
    create_add_inp(config)
    
    # Cleanup BASS.INP
    if os.path.isfile('BASS.INP'):
        run_shell('rm BASS.INP')
    
    # Check if ADDeven.INP and ADDodd.INP exist
    even_exists = os.path.isfile('ADDeven.INP')
    odd_exists = os.path.isfile('ADDodd.INP')
    
    # Run add and form CONF.INP from respective ADD.INP files
    parities = []
    if even_exists: 
        form_conf_inp('even', bin_dir)
        parities.append('even')
    if odd_exists: 
        form_conf_inp('odd', bin_dir)
        parities.append('odd')
    
    # Cleanup - remove add.in, ADD.INP, CONF.INP and CONF_.INP
    run_shell("rm add.in ADD.INP CONF.INP CONF_.INP")
    
    # Create a ci.qs job script if it doesn't exist yet
    if on_hpc:
        print('generating new ci.qs in ' + os.getcwd() + ' directory')
        script_name = write_job_script('.','ci', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
    
    # Copy ADD.INP and CONF.INP to all directories if gen_dir == True
    if gen_dir:
        root_dir = os.getcwd()
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
                    dir_name = dir_prefix+str(abs(c))
                    Path(dir_path+'/'+dir_name).mkdir(parents=True, exist_ok=True)
                    os.chdir(dir_name)
                    run_shell('pwd')
                    for parity in parities:
                        move_conf_inp(root_dir, parity, run_codes, include_lsj, write_hij)
                        # Submit CI job if run_codes == True
                        if on_hpc and run_codes: 
                            os.chdir(parity)
                            if script_name:
                                run_shell("sbatch " + script_name)
                            else:
                                print('job script was not submitted. check job script and submit manually.')
                            os.chdir('../')
                        else:
                            print("run_codes option is only available with HPC access")
                    os.chdir('../')
                if K_is_dict[K_is]:
                    os.chdir('../../')
                else:
                    os.chdir('../')
        else:
            if isinstance(code_method, list):
                dir_path = os.getcwd()
                for method in code_method:
                    Path(dir_path+'/'+method+'/basis').mkdir(parents=True, exist_ok=True)
                    os.chdir(method)
                    run_shell('pwd')
                    for parity in parities:
                        move_conf_inp(root_dir, parity, run_codes, include_lsj, write_hij)
                        # Submit CI job if run_codes == True
                        if on_hpc and run_codes: 
                            os.chdir(parity)
                            if script_name:
                                run_shell("sbatch " + script_name)
                            else:
                                print('job script was not submitted. check job script and submit manually.')
                            os.chdir('../')
                        else:
                            print("run_codes option is only available with HPC access - please run ci codes manually")
                    os.chdir('../')
            else:
                dir_path = os.getcwd()
                run_shell('pwd')
                for parity in parities:
                    move_conf_inp(root_dir, parity, run_codes, include_lsj, write_hij)
                    # Submit CI job if run_codes == True
                    if on_hpc and run_codes: 
                        os.chdir(parity)
                        if script_name:
                            run_shell("sbatch " + script_name)
                        else:
                            print('job script was not submitted. check job script and submit manually.')
                        os.chdir('../')
                    else:
                        print("run_codes option is only available with HPC access - please run ci codes manually")

    print('add script completed')

