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

def write_mbpt_inp(basis):
    core_orbitals = basis['orbitals']['core']
    Nso = 0
    for orbital in core_orbitals.split(' '):
        if orbital[-1] == 's':
            Nso += 1
        else:
            Nso += 2

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
        f.write('A_hf   0 | \n')
        f.write('B_hf   0 | \n')
        f.write('E1_L   1 | \n')
        f.write('EDM    0 | - Right hand side operators \n')
        f.write('PNC    0 | \n')
        f.write('E1_V   0 | \n')
        f.write('AM     0 | \n')
        f.write('MQM    0 | \n')
        f.write('M1     0 | \n')
        f.write('E2     0 | \n')
        f.write('E3     0 | \n')
        f.write('M2     0 | \n')
        f.write('M3     0 | \n')
        f.write('========================= \n')
        f.write('Nhf = 12 - SCF procedure includes Nhf upper shells (Nsh,Nsh-1,...) \n')
        f.write('Kmg =  0 - if not zero, Omega gives frequency of external field \n')
        f.write('Omega= 0.057580(a.u.) \n')
        f.write('Kex =  1 - key for exchange (0 - skip, 1 - include) \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')

    return


if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    code_method = config['optional']['code_method']
    calc_E1 = config['dtm']['E1']
    num_levels = config['conf']['num_energy_levels']
    include_rpa = config['dtm']['include_rpa']
    pci_version = config['optional']['pci_version']
    basis = config['basis']
    
    if isinstance(code_method, collections.abc.Sequence):
        dir_path = os.getcwd()
        for method in code_method:
            full_path = dir_path+'/'+method
            Path(full_path).mkdir(parents=True, exist_ok=True)
            os.chdir(full_path)
            
            # Make dtm directory with dtm input and job script
            Path(full_path+'/dtm').mkdir(parents=True, exist_ok=True)
            with open('dtm.in', 'w') as f:
                if include_rpa:
                    f.write('3\n')
                else:
                    f.write('2\n')
                f.write('1 ' + str(num_levels) + ' 1 ' + str(num_levels) + '\n')
                f.write('E1')

            run('mv dtm.in dtm/dtm.in', shell=True)
            if include_rpa:
                write_mbpt_inp(basis)
                write_job_script('.','dtm_rpa', 2, 64, True, 0, 'large-mem', pci_version)
                run('mv dtm_rpa.qs dtm/dtm_rpa.qs', shell=True)
            else:
                write_job_script('.','dtm', 2, 64, True, 0, 'large-mem', pci_version)
                run('mv dtm.qs dtm/dtm.qs', shell=True)
            
            # Find even and odd directories with completed ci runs
            if os.path.isfile('even/CONFFINAL.RES'):
                if include_rpa: 
                    run('cp even/HFD.DAT dtm/HFD.DAT', shell=True)
                    run('cp MBPT.INP dtm/MBPT.INP', shell=True)
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
            if include_rpa:
                with open('rpa.in', 'w') as f:
                    f.write('2')
                run('mpirun -n 1 dtm', shell=True)
                with open('dtm.in', 'w') as f:
                    f.write('2\n')
                    f.write('1 ' + str(num_levels) + ' 1 ' + str(num_levels) + '\n')
                    f.write('E1')
                run('sbatch dtm_rpa.qs', shell=True)
            else:    
                run('sbatch dtm.qs', shell=True)
            
    else:
        # TODO - run dtm for even + odd directories
        print('not supported yet')
        sys.exit()