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


if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    code_method = config['optional']['code_method']
    calc_E1 = config['dtm']['E1']
    num_levels = config['conf']['num_energy_levels']
    
    if isinstance(code_method, collections.abc.Sequence):
        dir_path = os.getcwd()
        for method in code_method:
            full_path = dir_path+'/'+method
            Path(full_path).mkdir(parents=True, exist_ok=True)
            os.chdir(full_path)
            
            # Make dtm directory with dtm input and job script
            Path(full_path+'/dtm').mkdir(parents=True, exist_ok=True)
            with open('dtm.in', 'w') as f:
                f.write('2\n')
                f.write('1 ' + str(num_levels) + ' 1 ' + str(num_levels) + '\n')
                f.write('E1')
            f.close()
            write_job_script('.','dtm', 2, 64, True, 0, 'standard')
            
            run('mv dtm.in dtm/dtm.in', shell=True)
            run('mv dtm.qs dtm/dtm.qs', shell=True)
            
            # Find even and odd directories with completed ci runs
            if os.path.isfile('even/CONFFINAL.RES'):
                run('cp even/CONF.INP dtm/CONF.INP', shell=True)
                run('cp even/CONF.DET dtm/CONF.DET', shell=True)
                run('cp even/CONF.XIJ dtm/CONF.XIJ', shell=True)
                run('cp even/CONFSTR.RES dtm/CONFSTR.RES', shell=True)
                run('cp even/CONF.DAT dtm/CONF.DAT', shell=True)
                run('cp even/CONF.INT dtm/CONF.INT', shell=True)
            if os.path.isfile('odd/CONFFINAL.RES'):
                run('cp odd/CONF.INP dtm/CONF1.INP', shell=True)
                run('cp odd/CONF.DET dtm/CONF1.DET', shell=True)
                run('cp odd/CONF.XIJ dtm/CONF1.XIJ', shell=True)
                run('cp odd/CONFSTR.RES dtm/CONFSTR1.RES', shell=True)
                        
            # cd into new dtm directory and submit job script
            dtm_path = full_path+'/dtm'
            os.chdir(dtm_path)
            run('sbatch dtm.qs', shell=True)
            
            # TODO - If include_rpa, run dtm to make DTM.INT, create MBPT.INP, then run rpa > rpa_dtm > dtm, else run dtm
    else:
        # TODO - run dtm for even + odd directories
        print('not supported yet')
        sys.exit()