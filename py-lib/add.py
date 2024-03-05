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

def move_conf_inp(parity, run_ci, include_lsj, pci_version):
    if os.path.isfile('../basis/HFD.DAT'):
        run("cp ../basis/HFD.DAT .", shell=True)
    if os.path.isfile('../CONF' + parity + '.INP'):
        run("cp ../CONF" + parity + ".INP CONF.INP", shell=True)
    if os.path.isfile('../basis/SGC.CON'):
        run("cp ../basis/SGC.CON .", shell=True)
    if os.path.isfile('../basis/SCRC.CON'):
        run("cp ../basis/SCRC.CON .", shell=True)
    if os.path.isfile('../basis/SGC.CON') and os.path.isfile('../basis/SCRC.CON'):
        with open('c.in', 'w') as f:
            if include_lsj:
                f.write('2, 2, 0, 0, 1')
            else:
                f.write('2, 2, 0, 0, 0')
            f.close()
    else:
        with open('c.in', 'w') as f:
            if include_lsj:
                f.write('0, 0, 0, 0, 1')
            else:
                f.write('0, 0, 0, 0, 0')
            f.close()
    if run_ci: 
        if not os.path.isfile('ci.qs'):
            print('generating new ci.qs in ' + os.getcwd() + ' directory')
        write_job_script('.','ci', 2, 64, True, 0, 'standard', pci_version)
        run("sbatch ci.qs", shell=True)
        
    os.chdir('../')

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
    
    code_method = config['optional']['code_method']
    run_ci = config['optional']['run_ci']
    include_lsj = config['conf']['include_lsj']
    pci_version = config['optional']['pci_version']
    
    # Ensure basis and add core orbitals match
    basis_core = config['basis']['orbitals']['core']
    add_core = config['add']['orbitals']['core']
    
    if basis_core != add_core: 
        print('ERROR: core orbitals of basis and add do not match in config.yml')
        sys.exit()
    
    if run_ci:
        gen_dir = run_ci
    else:
        gen_dir = config['optional']['generate_directories']
        
    if gen_dir:
        if isinstance(code_method, collections.abc.Sequence):
            dir_path = os.getcwd()
            for method in code_method:
                full_path = dir_path+'/'+method
                Path(full_path).mkdir(parents=True, exist_ok=True)
                os.chdir(full_path)
                
                # Read input to add from add.in if it exists, otherwise create it
                try:
                    open("add.in", 'rb').read()
                except:
                    f = open("add.in", "w")
                    f.write("0 \n 0")
                    f.close()
                
                # Pull BASS.INP from basis directory if it exists
                if os.path.isfile('basis/BASS.INP'):
                    run('cp basis/BASS.INP .', shell=True)
            
                # Run script to generate ADD.INP for even and odd configurations
                create_add_inp(config)
                
                # Cleanup BASS.INP
                run('rm BASS.INP', shell=True)
            
                # Check if ADDeven.INP and ADDodd.INP exist
                even_exists = os.path.isfile('ADDeven.INP')
                odd_exists = os.path.isfile('ADDodd.INP')
            
                if even_exists: form_conf_inp('even')
                if odd_exists: form_conf_inp('odd')
            
                # Cleanup - remove add.in
                run("rm add.in", shell=True)
                
                if even_exists: 
                    even_path = full_path + '/even'
                    Path(even_path).mkdir(parents=True, exist_ok=True)
                    os.chdir(even_path)
                    move_conf_inp('even', run_ci, include_lsj, pci_version)
                    os.chdir('../')
                if odd_exists: 
                    odd_path = full_path + '/odd'
                    Path(odd_path).mkdir(parents=True, exist_ok=True)
                    os.chdir(odd_path)
                    move_conf_inp('odd', run_ci, include_lsj, pci_version)
                    os.chdir('../')
                os.chdir('../')
        else:
            # Read input to add from add.in if it exists, otherwise create it
            try:
                open("add.in", 'rb').read()
            except:
                f = open("add.in", "w")
                f.write("0 \n 0")
                f.close()
            
            # Run script to generate ADD.INP for even and odd configurations
            create_add_inp(config)
            
            # Check if ADDeven.INP and ADDodd.INP exist
            even_exists = os.path.isfile('ADDeven.INP')
            odd_exists = os.path.isfile('ADDodd.INP')
            
            if even_exists: form_conf_inp('even')
            if odd_exists: form_conf_inp('odd')
            
            # Cleanup - remove add.in
            run("rm add.in", shell=True)
                
            dir_path = os.getcwd()
            for dirs in ['even','odd']:
                Path(dir_path+'/'+dirs).mkdir(parents=True, exist_ok=True)
                os.chdir(dirs)
                move_conf_inp(dirs)
                os.chdir('../')

    print('add script completed')

