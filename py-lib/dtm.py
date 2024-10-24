""" dtm

This script allows the user to automate the transition matrix element calculations from inputted parameters in a "config.yml" file. 
The "config.yml" file should have the following blocks:

    * system - system parameter (bin_directory, run_codes, on_hpc)
    * atom.code_method - list of methods (CI, CI+all-order, CI+MBPT)
    * basis.core_orbitals - list of core orbitals in basis set (used to form MBPT.INP)
    * dtm - parameters used by dtm program (include_rpa, DM, TM)

From these parameters, this script will create all input files required for execution of the add program.
After the input files are created, the add program will be executed to create the list of configurations CONF.INP.
The script will generate respective directories for DTM calculations (TM: /tm, DM: /dm_even, /dm_odd)
If optional.run_codes is set to "True", the script will then submit the slurm job in the generated directories. 

This python script has 2 main capabilities for polarizability calculations:
1. Set up file directory for matrix elements calculations.
2. Submit job script for matrix elements calculations.

"""
import yaml
import os
import sys
from pathlib import Path
from utils import run_shell, get_dict_value, check_slurm_installed
from gen_job_script import write_job_script


def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def write_mbpt_inp(basis, matrix_elements):
    core_orbitals = basis['orbitals']['core']
    Nso = 0
    for orbital in core_orbitals.split(' '):
        if orbital[-1] == 's':
            Nso += 1
        else:
            Nso += 2
    
    key_dict = {
        'A_hf': '0',
        'B_hf': '0',
        'E1_L': '0',
        'EDM': '0',
        'PNC': '0',
        'E1_V': '0',
        'AM': '0',
        'MQM': '0',
        'M1': '0',
        'E2': '0',
        'E3': '0',
        'M2': '0',
        'M3': '0'
    }
    for matrix_element in matrix_elements:
        key_dict[matrix_element] = '1'
        if matrix_element == 'E1':
            key_dict['E1_L'] = '1'
            key_dict['E1_V'] = '1'

    with open('MBPT.INP','w') as f:
        f.write('MBPT>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Used by ALL MBPT programs \n')
        f.write('Nso = ' + str(Nso) + ' - CI core \n')
        f.write('Nsh = ' + str(Nso) + ' - defines SCF field (For MS calculations Nsh=Nso) \n')
        f.write('Nss =999 \n')
        f.write('Nsv = ' + str(Nso+1) + ' - =Nso+1 \n')
        f.write('Nmax=210 \n')
        f.write('Lmax=  4 - max (l_i,l_k) for valence radial integrals \n')
        f.write('Kmax=  9 - max multipolarity of two-electron valence integrals \n')
        f.write('Kt  =  1 - Keep this fixed \n')
        f.write('Kbrt=  2 - Breit interaction \n')
        f.write('Kout=  0 - Details in output file \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')
        f.write('SMS:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> (p_i dot p_k) interaction \n')
        f.write('C_sms= 0.00000   - scaling of SMS interaction \n')
        f.write('Klow=  1         - lower component included/ignored \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')
        f.write('RPA >>>>>>>>>>>>>>>>>>Used by programs which deal with RPA MEs \n')
        f.write('A_hf   ' + key_dict['A_hf'] + ' | \n')
        f.write('B_hf   ' + key_dict['B_hf'] + ' | \n')
        f.write('E1_L   ' + key_dict['E1_L'] + ' | \n')
        f.write('EDM    ' + key_dict['EDM']  + ' | - Right hand side operators \n')
        f.write('PNC    ' + key_dict['PNC']  + ' | \n')
        f.write('E1_V   ' + key_dict['E1_V'] + ' | \n')
        f.write('AM     ' + key_dict['AM']   + ' | \n')
        f.write('MQM    ' + key_dict['MQM']  + ' | \n')
        f.write('M1     ' + key_dict['M1']   + ' | \n')
        f.write('E2     ' + key_dict['E2']   + ' | \n')
        f.write('E3     ' + key_dict['E3']   + ' | \n')
        f.write('M2     ' + key_dict['M2']   + ' | \n')
        f.write('M3     ' + key_dict['M3']   + ' | \n')
        f.write('========================= \n')
        f.write('Nhf = ' + str(Nso) + ' - SCF procedure includes Nhf upper shells (Nsh,Nsh-1,...) \n')
        f.write('Kmg =  0 - if not zero, Omega gives frequency of external field \n')
        f.write('Omega= 0.057580(a.u.) \n')
        f.write('Kex =  1 - key for exchange (0 - skip, 1 - include) \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')

    return

def write_dtm_in(mode, levels, operators):
    """ Write dtm.in """
    with open('dtm.in','w') as f:
        f.write('Mode = ' + mode + '\n')
        if levels:
            f.write('Levels = ' + levels + '\n')
        if operators:
            f.write('Operators = ' + operators)

def gen_key_list(matrix_elements):
    key_list = []
    if matrix_elements:
        if isinstance(matrix_elements, list):
            key_list = matrix_elements
        else:
            if len(matrix_elements.split(' ')) == 1:
                key_list = [matrix_elements]
            else:
                for matrix_element in matrix_elements.split(' '):
                    key_list.append(matrix_element.replace('[','').replace(']','').replace(',',''))
    return key_list

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    
    atom = get_dict_value(config, 'atom')
    code_method = get_dict_value(atom, 'code_method')
    if not isinstance(code_method, list): code_method = [code_method]
    
    system = get_dict_value(config, 'system')
    pci_version = get_dict_value(system, 'pci_version')
    on_hpc = get_dict_value(system, 'on_hpc')
    bin_dir = get_dict_value(system, 'bin_directory')
    run_codes = get_dict_value(system, 'run_codes')
    
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
    
    basis = get_dict_value(config, 'basis')
    conf = get_dict_value(config, 'conf')
    odd_J = get_dict_value(conf['odd'], 'J')
    even_J = get_dict_value(conf['even'], 'J')
    conf_odd_path = 'odd' + str(odd_J)[0]
    conf_even_path = 'even' + str(even_J)[0]

    dtm = get_dict_value(config, 'dtm')
    include_rpa = get_dict_value(dtm, 'include_rpa')
    
    # Create an array of directories for dtm
    dtm_dirs = []
    
    # DM parameters
    dm = get_dict_value(dtm, 'DM')
    dm_parity = get_dict_value(dm, 'parity')
    dm_matrix_elements = get_dict_value(dm, 'matrix_elements')
    if not dm_parity: dm_parity = 'both'
    dm_range = get_dict_value(dm, 'level_range')
    include_dm = True if dm_range else False
    if include_dm:
        dm_odd = get_dict_value(dm_range, 'odd')
        dm_even = get_dict_value(dm_range, 'even')
        if dm_odd:
            from_level_odd = dm_odd.split(' ')[0]
            to_level_odd = dm_odd.split(' ')[1]
            dtm_dirs.append('dm_odd')
        if dm_even:
            from_level_even = dm_even.split(' ')[0]
            to_level_even = dm_even.split(' ')[1]
            dtm_dirs.append('dm_even')
    
    # TM parameters
    tm = get_dict_value(dtm, 'TM')
    tm_matrix_elements = get_dict_value(tm, 'matrix_elements')
    tm_from = get_dict_value(tm, 'from')
    tm_from_parity = get_dict_value(tm_from, 'parity')
    tm_from_range = get_dict_value(tm_from, 'level_range')
    tm_to = get_dict_value(tm, 'to')
    tm_to_parity = get_dict_value(tm_to, 'parity')
    tm_to_range = get_dict_value(tm_to, 'level_range')
    include_tm = True if tm_from_range and tm_to_range else False
    if include_tm:
        from_level_initial = tm_from_range.split(' ')[0]
        from_level_final = tm_from_range.split(' ')[1]
        to_level_initial = tm_to_range.split(' ')[0]
        to_level_final = tm_to_range.split(' ')[1]
        dtm_dirs.append('tm')
    
    # list of keys for dtm operators
    dm_key_list = gen_key_list(dm_matrix_elements)
    tm_key_list = gen_key_list(tm_matrix_elements)

    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
    dir_path = os.getcwd()
    for method in code_method:
        if len(code_method) > 1:
            full_path = dir_path+'/'+method
        else:
            full_path = dir_path
        os.chdir(full_path)
        
        even_exists, odd_exists = False, False
        if tm_from_parity == 'even' or tm_to_parity == 'even' or dm_parity == 'even' or dm_parity == 'both':
            if os.path.isfile(conf_even_path + '/CONF.RES') or os.path.isfile(conf_even_path + '/CONFFINAL.RES'):
                    even_exists = True
        if tm_from_parity == 'odd' or tm_to_parity == 'odd' or dm_parity == 'odd' or dm_parity == 'both':
            if os.path.isfile(conf_odd_path + '/CONF.RES') or os.path.isfile(conf_odd_path + '/CONFFINAL.RES'):
                    odd_exists = True
        
        # Create dtm directories with dtm input and job script
        for dtm_dir in dtm_dirs:
            Path(full_path+'/'+dtm_dir).mkdir(parents=True, exist_ok=True)
            if include_rpa:
                write_dtm_in('Init','','')
                if 'tm' in dtm_dir:
                    write_mbpt_inp(basis, tm_key_list)
                elif 'dm' in dtm_dir:
                    write_mbpt_inp(basis, dm_key_list)
                run_shell('cp dtm.in ' + dtm_dir + '/dtm.in')
                run_shell('cp MBPT.INP ' + dtm_dir + '/MBPT.INP')
                run_shell('cp basis/HFD.DAT ' + dtm_dir + '/HFD.DAT')
                with open('rpa.in', 'w') as f:
                    f.write('2')
                run_shell('cp rpa.in ' + dtm_dir + '/rpa.in')
                
                if on_hpc:
                    script_name = write_job_script('.','dtm_rpa', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
                if script_name:
                    run_shell('cp dtm_rpa.qs ' + dtm_dir + '/dtm_rpa.qs')
            else:
                if dtm_dir == 'tm':
                    levels = from_level_initial + ' ' + from_level_final + ', ' + to_level_initial + ' ' + to_level_final
                    write_dtm_in('TM', levels, ', '.join(tm_key_list))
                elif dtm_dir == 'dm_even':
                    write_dtm_in('DM', from_level_even + ' ' + to_level_even, ', '.join(dm_key_list))
                elif dtm_dir =='dm_odd':
                    write_dtm_in('DM', from_level_odd + ' ' + from_level_odd, ', '.join(dm_key_list))
                if on_hpc:
                    script_name = write_job_script('.','dtm', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
                if script_name:
                    run_shell('cp dtm.qs ' + dtm_dir + '/dtm.qs')
            run_shell('cp dtm.in ' + dtm_dir + '/dtm.in')

            if not even_exists and not odd_exists:
                print('ci directories could not be found')
                sys.exit()
            else:
                # if including rpa, generate an initial DTM.INT for rpa programs             
                if dtm_dir == 'tm':
                    if tm_from_parity == 'even': 
                        from_path = conf_even_path
                    elif tm_from_parity == 'odd':
                        from_path = conf_odd_path
                    run_shell('cp ' + from_path  + '/CONF.INP ' + dtm_dir + '/CONF.INP')
                    run_shell('cp ' + from_path  + '/CONF.DET ' + dtm_dir + '/CONF.DET')
                    run_shell('cp ' + from_path  + '/CONF.XIJ ' + dtm_dir + '/CONF.XIJ')
                    run_shell('cp ' + from_path  + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR.RES')
                    run_shell('cp ' + from_path  + '/CONF.DAT ' + dtm_dir + '/CONF.DAT')
                    run_shell('cp ' + from_path  + '/CONF.INT ' + dtm_dir + '/CONF.INT')

                    if tm_to_parity == 'even': 
                        to_path = conf_even_path
                    elif tm_to_parity == 'odd':
                        to_path = conf_odd_path
                    run_shell('cp ' + to_path + '/CONF.INP ' + dtm_dir + '/CONF1.INP')
                    run_shell('cp ' + to_path + '/CONF.DET ' + dtm_dir + '/CONF1.DET')
                    run_shell('cp ' + to_path + '/CONF.XIJ ' + dtm_dir + '/CONF1.XIJ')
                    run_shell('cp ' + to_path + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR1.RES')
                elif dtm_dir == 'dm_even':
                    run_shell('cp ' + conf_even_path + '/CONF.INP ' + dtm_dir + '/CONF.INP')
                    run_shell('cp ' + conf_even_path + '/CONF.DET ' + dtm_dir + '/CONF.DET')
                    run_shell('cp ' + conf_even_path + '/CONF.XIJ ' + dtm_dir + '/CONF.XIJ')
                    run_shell('cp ' + conf_even_path + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR.RES')
                    run_shell('cp ' + conf_even_path + '/CONF.DAT ' + dtm_dir + '/CONF.DAT')
                    run_shell('cp ' + conf_even_path + '/CONF.INT ' + dtm_dir + '/CONF.INT')
                elif dtm_dir == 'dm_odd':
                    run_shell('cp ' + conf_odd_path + '/CONF.INP ' + dtm_dir + '/CONF.INP')
                    run_shell('cp ' + conf_odd_path + '/CONF.DET ' + dtm_dir + '/CONF.DET')
                    run_shell('cp ' + conf_odd_path + '/CONF.XIJ ' + dtm_dir + '/CONF.XIJ')
                    run_shell('cp ' + conf_odd_path + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR.RES')
                    run_shell('cp ' + conf_odd_path + '/CONF.DAT ' + dtm_dir + '/CONF.DAT')
                    run_shell('cp ' + conf_odd_path + '/CONF.INT ' + dtm_dir + '/CONF.INT')
        
        # cd into new dtm directories and submit job scripts if run_codes
        if run_codes:
            for dtm_dir in dtm_dirs:
                dtm_path = full_path+'/'+dtm_dir
                os.chdir(dtm_path)
                if include_rpa:
                    if on_hpc: 
                        run_shell('mpirun -n 1 ' + bin_dir + 'pdtm > dtm0.out')
                        if dtm_dir == 'tm':
                            levels = from_level_initial + ' ' + from_level_final + ', ' + to_level_initial + ' ' + to_level_final
                            write_dtm_in('TM', levels, ', '.join(tm_key_list))
                        elif dtm_dir == 'dm_even':
                            write_dtm_in('DM', from_level_even + ' ' + to_level_even, ', '.join(dm_key_list))
                        elif dtm_dir =='dm_odd':
                            write_dtm_in('DM', from_level_odd + ' ' + from_level_odd, ', '.join(dm_key_list))
                        if submit_job:
                            run_shell('sbatch dtm_rpa.qs')
                else:    
                    if submit_job:
                        run_shell('sbatch dtm.qs')
                os.chdir('..')
            os.chdir('..')
