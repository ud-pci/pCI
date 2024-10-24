""" pol

This script allows the user to automate polarizability calculation from inputted parameters in a "config.yml" file. 
The "config.yml" file should have the following blocks:

    * system - system parameter (bin_directory, run_codes, on_hpc)
    * atom.code_method - list of methods (CI, CI+all-order, CI+MBPT)
    * pol - parameters used by pol program (parity, level, method, wavelength_range, step_size)

From these parameters, this script will create move the required input files for execution of the pol program.
After the input files are moved, a job script will be generated to run the job if working on HPC. 
If optional.run_codes is set to "True", the script will then submit the slurm job in the generated directory. 

This python script has 2 main capabilities for polarizability calculations:
1. Set up file directory for polarizability calculation.
2. Submit job script for polarizability calculation.

"""
import yaml
import os
import sys
from pathlib import Path
from utils import run_shell, get_dict_value, convert_params_to_list, check_slurm_installed
from gen_job_script import write_job_script


def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def write_pol_in(method, level, include_static, include_dynamic, wavelength_range, step_size):
    """ Write pol.in """
    with open('pol.in','w') as f:
        f.write('Mode = 0 \n')
        f.write('Method = ' + str(method) + '\n')
        f.write('Level = ' + str(level) + '\n')
        f.write('Ranges = ')
        if include_static:
            f.write('0 0 0, ')
        if include_dynamic:
            ranges = ''
            for irange in range(len(wavelength_range)):
                ranges += wavelength_range[irange] + ' ' + step_size[irange] + ', '
            ranges = ranges[:-2]
        f.write(ranges)

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    
    # system parameters
    system = get_dict_value(config, 'system')
    pci_version = get_dict_value(system, 'pci_version')
    on_hpc = get_dict_value(system, 'on_hpc')
    bin_dir = get_dict_value(system, 'bin_directory')
    run_codes = get_dict_value(system, 'run_codes')
    
    atom = get_dict_value(config, 'atom')
    code_method = get_dict_value(atom, 'code_method')
    if not isinstance(code_method, list): code_method = [code_method]
    
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
    
    conf = get_dict_value(config, 'conf')
    odd_J = get_dict_value(conf['odd'], 'J')
    even_J = get_dict_value(conf['even'], 'J')
    conf_odd_path = 'odd' + str(odd_J)[0]
    conf_even_path = 'even' + str(even_J)[0]
    
    # pol parameters
    pol = get_dict_value(config, 'pol')
    parity = get_dict_value(pol, 'parity')
    level = get_dict_value(pol, 'level')
    pol_method = get_dict_value(pol, 'method')
    field_type = convert_params_to_list(get_dict_value(pol, 'field_type'))
    include_static = True if 'static' in field_type else False
    include_dynamic = True if 'dynamic' in field_type else False
    
    if include_dynamic:
        wavelength_range = convert_params_to_list(get_dict_value(pol, 'wavelength_range'))
        step_size = convert_params_to_list(get_dict_value(pol, 'step_size'))
    
    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
    # Each code_method will have its own directory
    dir_path = os.getcwd()
    for method in code_method:
        if len(code_method) > 1:
            full_path = dir_path+'/'+method
        else:
            full_path = dir_path
        Path(full_path).mkdir(parents=True, exist_ok=True)
        os.chdir(full_path)

        # Make pol directory with pol input and job script
        pol_path = full_path+'/pol_'+parity+str(level)
        Path('pol_'+parity+str(level)).mkdir(parents=True, exist_ok=True)
        write_pol_in(pol_method, level, include_static, include_dynamic, wavelength_range, step_size)

        run_shell('mv pol.in ' + pol_path + '/pol.in')

        script_name = write_job_script('.','pol', nodes, tasks_per_node, True, 0, partition, pci_version,bin_dir)
        if script_name:
            run_shell('mv pol.qs ' + pol_path + '/pol.qs')
        else:
            print('job script was not submitted. check job script and submit manually.')

        # Find even and odd directories with completed ci runs
        even_exists, odd_exists = False, False
        if os.path.isfile(conf_even_path + '/CONF.RES') or os.path.isfile(conf_even_path + '/CONFFINAL.RES'):
            even_exists = True
            if parity == 'even':
                run_shell('cp ' + conf_even_path + '/CONF.DET ' + pol_path + '/CONF0.DET')
                run_shell('cp ' + conf_even_path + '/CONF.XIJ ' + pol_path + '/CONF0.XIJ')
            elif parity == 'odd':
                run_shell('cp ' + conf_even_path + '/CONF.DET ' + pol_path + '/CONF.DET')
                run_shell('cp ' + conf_even_path + '/CONF.XIJ ' + pol_path + '/CONF.XIJ')
                run_shell('cp ' + conf_even_path + '/CONF.INP ' + pol_path + '/CONF.INP')
                run_shell('cp ' + conf_even_path + '/CONF.DAT ' + pol_path + '/CONF.DAT')
                run_shell('cp ' + conf_even_path + '/nprocs.conf ' + pol_path + '/nprocs.conf')
                run_shell('cp ' + conf_even_path + '/CONFp.HIJ ' + pol_path + '/CONFp.HIJ')
                run_shell('cp ' + conf_even_path + '/CONFp.JJJ ' + pol_path + '/CONFp.JJJ')
                
        if os.path.isfile(conf_odd_path + '/CONF.RES') or os.path.isfile(conf_odd_path + '/CONFFINAL.RES'):
            odd_exists = True
            if parity == 'odd':
                run_shell('cp ' + conf_odd_path + '/CONF.DET ' + pol_path + '/CONF0.DET')
                run_shell('cp ' + conf_odd_path + '/CONF.XIJ ' + pol_path + '/CONF0.XIJ')
            elif parity == 'even':
                run_shell('cp ' + conf_odd_path + '/CONF.DET ' + pol_path + '/CONF.DET')
                run_shell('cp ' + conf_odd_path + '/CONF.XIJ ' + pol_path + '/CONF.XIJ')
                run_shell('cp ' + conf_odd_path + '/CONF.INP ' + pol_path + '/CONF.INP')
                run_shell('cp ' + conf_odd_path + '/CONF.DAT ' + pol_path + '/CONF.DAT')
                run_shell('cp ' + conf_odd_path + '/nprocs.conf ' + pol_path + '/nprocs.conf')
                run_shell('cp ' + conf_odd_path + '/CONFp.HIJ ' + pol_path + '/CONFp.HIJ')
                run_shell('cp ' + conf_odd_path + '/CONFp.JJJ ' + pol_path + '/CONFp.JJJ')

        if not even_exists and not odd_exists:
            print('ci directories could not be found')
            sys.exit()
        
        # Move DTM.INT to pol directory
        run_shell('cp tm/DTM.INT ' + pol_path + '/DTM.INT')

        # cd into new pol directory and submit job script
        os.chdir(pol_path)
        if submit_job:
            run_shell('sbatch pol.qs')
