import yaml
import collections.abc
from subprocess import run
import os
import sys
import re
import math
import orbitals as orb_lib
import get_atomic_data as libatomic
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
        
    config_sys = config['system']
    config_conf = config['conf']
    # Write head of CONF.INP
    f.write('>>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>\n')
    f.write('  ' + config_sys['name'] + ' ' + parity + '\n')
    f.write('  Z = ' + str(Z) + '\n')
    f.write(' Am = ' + "{:.1f}".format(round(AM)) + '\n')
    f.write('  J = ' + str(config_conf['J']) + '\n')
    f.write(' Jm = ' + str(config_conf['JM']) + '\n')
    f.write(' Nso=  ' + str(num_core_orb) + '\n')
    f.write(' Nc =   10 \n')
    if config_conf['J_selection']:
        f.write(' Kv =  3 \n')
    else:
        f.write(' Kv =  4 \n')
    f.write(' Nlv=  ' + str(config_conf['num_energy_levels']) + '\n')
    f.write(' Ne =  ' + str(num_val) + '\n')
    f.write(' Kl4=  1 \n')
    f.write(' Nc4=999 \n')
    f.write(' Gj = 0.0000 \n')
    f.write('Crt4= 0.0001 \n')
    f.write('kout= 0 \n')
    f.write('Ncpt= 0 \n')
    f.write('Cut0= 0.0001 \n')
    f.write('N_it= ' + str(config_conf['num_dvdsn_iterations']) + '\n')
    if config_sys['include_breit']:
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
    config_sys = config['system']
    Z, AM, symbol, cfermi, rnuc, num_rem_ele = libatomic.get_atomic_data(config_sys['name'], config_sys['isotope'])

    config_add = config['add']
    num_val = orb_lib.count_valence(config_add['ref_configs'])
    orb_occ = orb_lib.expand_orbitals(config_add['basis_set'], config_add['orbitals'])
    multiplicity = orb_lib.count_excitations(config_add['excitations'])

    if config_add['ref_configs']['even']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'even')
    else:
        print('no even reference configurations specified')
        
    if config_add['ref_configs']['odd']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'odd')
    else:
        print('no odd reference configurations specified')

def form_conf_inp(parity):
    run("cp ADD" + parity + ".INP ADD.INP", shell=True)
    run("add < add.in", shell=True)
    run("cp CONF.INP CONF" + parity + ".INP", shell=True)
    print("CONF" + parity + ".INP created")

def move_conf_inp(root_dir, parity, run_ci, include_lsj, write_hij):
    if not os.path.isdir(parity):
        run("mkdir " + parity, shell=True)
    if os.path.isfile('basis/HFD.DAT'):
        run("cp basis/HFD.DAT " + parity, shell=True)
    if os.path.isfile(root_dir + '/CONF' + parity + '.INP'):
        run("cp " + root_dir + "/CONF" + parity + ".INP " + parity + "/CONF.INP", shell=True )
    if os.path.isfile('basis/SGC.CON'):
        run("cp basis/SGC.CON "  + parity, shell=True)
    if os.path.isfile('basis/SCRC.CON'):
        run("cp basis/SCRC.CON "  + parity, shell=True)
    if run_ci and os.path.isfile(root_dir + '/ci.qs'):
        run("cp " + root_dir + "/ci.qs " + parity, shell=True)
        
    Kw = '1' if write_hij else '0'
    kLSJ = '1' if include_lsj else '0'

    if os.path.isfile('basis/SGC.CON') and os.path.isfile('basis/SCRC.CON'):
        with open(parity + '/c.in', 'w') as f:
            f.write('2, 2, 0, ' + Kw + ', ' + kLSJ)
    else:
        with open(parity + '/c.in', 'w') as f:
            f.write('0, 0, 0, ' + Kw + ', ' + kLSJ)

def check_conf_inp_exists(dir):
    if not os.path.isfile('HFD.DAT'):
        print('HFD.DAT is missing from ' + dir + ' directory')
        sys.exit()
    if not os.path.isfile('CONF.INP'):
        print('CONF.INP is missing from ' + dir + ' directory')
        sys.exit()
    if not os.path.isfile('c.in'):
        print('c.in is missing from ' + dir + ' directory')
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
        
    code_method = config['optional']['code_method']
    run_ci = config['optional']['run_ci']
    include_lsj = config['conf']['include_lsj']
    pci_version = config['optional']['pci_version']
    write_hij = config['conf']['write_hij']
    
    # Ensure basis and add core orbitals match
    basis_core = config['basis']['orbitals']['core']
    add_core = config['add']['orbitals']['core']
    
    if basis_core != add_core: 
        print('ERROR: core orbitals of basis and add do not match in config.yml')
        sys.exit()
    
    # Check if user wants to generate directories for CI computations
    gen_dir = run_ci if run_ci else config['optional']['generate_directories']
    
    # CONF.INP should only need to be constructed once, then copied to respective directories
    # Read input to add from add.in if it exists, otherwise create it
    try:
        open("add.in", 'rb').read()
    except:
        f = open("add.in", "w")
        f.write("0 \n 0")
        f.close()
    
    # BASS.INP will dictate order of orbitals in ADD.INP, so find a BASS.INP file from any basis directory
    if not os.path.isfile('BASS.INP'):
        for root, dirs, files in os.walk('.'):
            if 'BASS.INP' in files:
                run('cp ' + root + '/BASS.INP .', shell=True)
                break
    
    # Run script to generate ADD.INP for even and odd configurations
    create_add_inp(config)
    
    # Cleanup BASS.INP
    run('rm BASS.INP', shell=True)
    
    # Check if ADDeven.INP and ADDodd.INP exist
    even_exists = os.path.isfile('ADDeven.INP')
    odd_exists = os.path.isfile('ADDodd.INP')
    
    # Run add and form CONF.INP from respective ADD.INP files
    parities = []
    if even_exists: 
        form_conf_inp('even')
        parities.append('even')
    if odd_exists: 
        form_conf_inp('odd')
        parities.append('odd')
    
    # Cleanup - remove add.in, ADD.INP, CONF.INP and CONF_.INP
    run("rm add.in ADD.INP CONF.INP CONF_.INP", shell=True)
    
    # Create a ci.qs job script if it doesn't exist yet
    if not os.path.isfile('ci.qs'):
        print('generating new ci.qs in ' + os.getcwd() + ' directory')
        script_name = write_job_script('.','ci', 2, 64, True, 0, 'standard', pci_version)
    
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
                    run('pwd', shell=True)
                    for parity in parities:
                        move_conf_inp(root_dir, parity, run_ci, include_lsj, write_hij)
                        # Submit CI job if run_ci == True
                        if run_ci: 
                            os.chdir(parity)
                            run("sbatch " + script_name, shell=True)
                            os.chdir('../')
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
                    run('pwd', shell=True)
                    for parity in parities:
                        move_conf_inp(root_dir, parity, run_ci, include_lsj, write_hij)
                        # Submit CI job if run_ci == True
                        if run_ci: 
                            os.chdir(parity)
                            run("sbatch " + script_name, shell=True)
                            os.chdir('../')
                    os.chdir('../')
            else:
                dir_path = os.getcwd()
                run('pwd', shell=True)
                for parity in parities:
                    move_conf_inp(root_dir, parity, run_ci, include_lsj, write_hij)
                    # Submit CI job if run_ci == True
                    if run_ci: 
                        os.chdir(parity)
                        run("sbatch " + script_name, shell=True)
                        os.chdir('../')

    print('add script completed')

