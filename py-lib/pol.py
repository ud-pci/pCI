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
from utils import run_shell, get_dict_value
from gen_job_script import write_job_script


def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

def write_pol_in(method, level, wavelength_range, step_size):
    """ Write pol.in """
    with open('pol.in','w') as f:
        f.write('Mode = 0 \n')
        f.write('Method = ' + str(method) + '\n')
        f.write('Level = ' + str(level) + '\n')
        f.write('Ranges = ' + wavelength_range + ' ' + str(step_size))

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    code_method = config['atom']['code_method']
    pci_version = config['system']['pci_version']
    on_hpc = config['system']['on_hpc']
    bin_dir = config['system']['bin_directory']
    run_codes = config['system']['run_codes']
    
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
    
    # pol parameters
    pol = config['pol']
    parity = pol['parity']
    level = pol['level']
    pol_method = pol['method']
    wavelength_range = pol['wavelength_range']
    step_size = pol['step_size']
    
    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
    # Each code_method will have its own directory
    if isinstance(code_method, list):
        dir_path = os.getcwd()
        for method in code_method:
            full_path = dir_path+'/'+method
            Path(full_path).mkdir(parents=True, exist_ok=True)
            os.chdir(full_path)

            # Make pol directory with pol input and job script
            pol_path = full_path+'/pol_'+parity+str(level)
            Path('pol_'+parity+str(level)).mkdir(parents=True, exist_ok=True)
            write_pol_in(pol_method, level, wavelength_range, step_size)

            run_shell('mv pol.in ' + pol_path + '/pol.in')

            script_name = write_job_script('.','pol', nodes, tasks_per_node, True, 0, partition, pci_version,bin_dir)
            if script_name:
                run_shell('mv pol.qs ' + pol_path + '/pol.qs')
            else:
                print('job script was not submitted. check job script and submit manually.')
            
            # Check if pol_path already has CONFp.HIJ and CONFp.JJJ
            hij_exists = True if os.path.isfile(pol_path + '/CONF.HIJ') else False
            jjj_exists = True if os.path.isfile(pol_path + '/CONF.JJJ') else False
            phij_exists = True if os.path.isfile(pol_path + '/CONFp.HIJ') else False
            pjjj_exists = True if os.path.isfile(pol_path + '/CONFp.JJJ') else False

            # Find even and odd directories with completed ci runs
            even_exists, odd_exists = False, False
            if os.path.isfile('even/CONF.RES') or os.path.isfile('even/CONFFINAL.RES'):
                even_exists = True
                if parity == 'even':
                    run_shell('cp even/CONF.DET ' + pol_path + '/CONF0.DET')
                    run_shell('cp even/CONF.XIJ ' + pol_path + '/CONF0.XIJ')
                elif parity == 'odd':
                    run_shell('cp even/CONF.DET ' + pol_path + '/CONF.DET')
                    run_shell('cp even/CONF.XIJ ' + pol_path + '/CONF.XIJ')
                    run_shell('cp even/CONF.INP ' + pol_path + '/CONF.INP')
                    run_shell('cp even/CONF.DAT ' + pol_path + '/CONF.DAT')
                    if not hij_exists or not jjj_exists:
                        run_shell('cp even/nprocs.conf ' + pol_path + '/nprocs.conf')
                    if not hij_exists and not phij_exists: run_shell('cp even/CONFp.HIJ ' + pol_path + '/CONFp.HIJ')
                    if not jjj_exists and not pjjj_exists: run_shell('cp even/CONFp.JJJ ' + pol_path + '/CONFp.JJJ')
            if os.path.isfile('odd/CONF.RES') or os.path.isfile('odd/CONFFINAL.RES'):
                odd_exists = True
                if parity == 'odd':
                    run_shell('cp odd/CONF.DET ' + pol_path + '/CONF0.DET')
                    run_shell('cp odd/CONF.XIJ ' + pol_path + '/CONF0.XIJ')
                elif parity == 'even':
                    run_shell('cp odd/CONF.DET ' + pol_path + '/CONF.DET')
                    run_shell('cp odd/CONF.XIJ ' + pol_path + '/CONF.XIJ')
                    run_shell('cp odd/CONF.INP ' + pol_path + '/CONF.INP')
                    run_shell('cp odd/CONF.DAT ' + pol_path + '/CONF.DAT')
                    if not hij_exists or not jjj_exists:
                        run_shell('cp even/nprocs.conf ' + pol_path + '/nprocs.conf')
                    if not hij_exists and not phij_exists: run_shell('cp odd/CONFp.HIJ ' + pol_path + '/CONFp.HIJ')
                    if not jjj_exists and not pjjj_exists: run_shell('cp odd/CONFp.JJJ ' + pol_path + '/CONFp.JJJ')

            if not even_exists and not odd_exists:
                print('ci directories could not be found')
                sys.exit()
            
            # Move DTM.INT to pol directory
            run_shell('cp dtm/DTM.INT ' + pol_path + '/DTM.INT')

            # cd into new pol directory and submit job script
            os.chdir(pol_path)
            if on_hpc and run_codes:
                run_shell('sbatch pol.qs')
    else:
        dir_path = os.getcwd()
        