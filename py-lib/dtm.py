import yaml
import os
import sys
from pathlib import Path
from utils import run_shell
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
        f.write('Lvls = ' + levels + '\n')
        f.write('Ops = ' + operators)

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)
    code_method = config['atom']['code_method']
    matrix_elements = config['dtm']['matrix_elements']
    num_levels = config['conf']['num_energy_levels']
    include_rpa = config['dtm']['include_rpa']
    pci_version = config['system']['pci_version']
    on_hpc = config['system']['on_hpc']
    bin_dir = config['system']['bin_directory']
    basis = config['basis']
    
    key_list = []
    if isinstance(matrix_elements, list):
        key_list = matrix_elements
    else:
        if len(matrix_elements.split(' ')) == 1:
            key_list = [matrix_elements]
        else:
            for matrix_element in matrix_elements.split(' '):
                key_list.append(matrix_element.replace('[','').replace(']','').replace(',',''))

    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    
    if isinstance(code_method, list):
        dir_path = os.getcwd()
        for method in code_method:
            full_path = dir_path+'/'+method
            Path(full_path).mkdir(parents=True, exist_ok=True)
            os.chdir(full_path)
            
            # Make dtm directory with dtm input and job script
            Path(full_path+'/dtm').mkdir(parents=True, exist_ok=True)
            if include_rpa:
                write_dtm_in('Init',
                             '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                             ', '.join(key_list))
            else:
                write_dtm_in('TM',
                             '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                             ', '.join(key_list))

            run_shell('mv dtm.in dtm/dtm.in')
            if include_rpa:
                write_mbpt_inp(basis, key_list)
                if on_hpc:
                    write_job_script('.','dtm_rpa', 2, 64, True, 0, 'large-mem', pci_version, bin_dir)
                    run_shell('mv dtm_rpa.qs dtm/dtm_rpa.qs')
            else:
                write_job_script('.','dtm', 2, 64, True, 0, 'large-mem', pci_version, bin_dir)
                run_shell('mv dtm.qs dtm/dtm.qs')
            
            # Find even and odd directories with completed ci runs
            even_exists, odd_exists = False, False
            if os.path.isfile('even/CONF.RES') or os.path.isfile('even/CONFFINAL.RES'):
                if include_rpa: 
                    run_shell('cp even/HFD.DAT dtm/HFD.DAT')
                    run_shell('cp MBPT.INP dtm/MBPT.INP')
                run_shell('cp even/CONF.INP dtm/CONF.INP')
                run_shell('cp even/CONF.DET dtm/CONF.DET')
                run_shell('cp even/CONF.XIJ dtm/CONF.XIJ')
                run_shell('cp even/CONFSTR.RES dtm/CONFSTR.RES')
                run_shell('cp even/CONF.DAT dtm/CONF.DAT')
                run_shell('cp even/CONF.INT dtm/CONF.INT')
            if os.path.isfile('odd/CONF.RES') or os.path.isfile('odd/CONFFINAL.RES'):
                odd_exists = True
                run_shell('cp odd/CONF.INP dtm/CONF1.INP')
                run_shell('cp odd/CONF.DET dtm/CONF1.DET')
                run_shell('cp odd/CONF.XIJ dtm/CONF1.XIJ')
                run_shell('cp odd/CONFSTR.RES dtm/CONFSTR1.RES')
                
            if not even_exists and not odd_exists:
                print('ci directories could not be found')
                sys.exit()
                        
            # cd into new dtm directory and submit job script
            dtm_path = full_path+'/dtm'
            os.chdir(dtm_path)
            if include_rpa:
                with open('rpa.in', 'w') as f:
                    f.write('2')
                    
                run_shell('mpirun -n 1 ' + bin_dir + 'pdtm')
                write_dtm_in('TM',
                             '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                             ', '.join(key_list))
                
                if on_hpc: 
                    run_shell('sbatch dtm_rpa.qs')
            else:    
                if on_hpc:
                    run_shell('sbatch dtm.qs')
            
    else:
        # Make dtm directory with dtm input and job script
        dir_path = os.getcwd()
        dtm_path = dir_path+'/dtm'
        Path(dir_path+'/dtm').mkdir(parents=True, exist_ok=True)
        if include_rpa:
            write_dtm_in('Init',
                         '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                         ', '.join(key_list))
        else:
            write_dtm_in('TM',
                         '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                         ', '.join(key_list))

        run_shell('mv dtm.in dtm/dtm.in')
        if include_rpa:
            write_mbpt_inp(basis, key_list)
            if on_hpc:
                write_job_script('.','dtm_rpa', 2, 64, True, 0, 'large-mem', pci_version, bin_dir)
                run_shell('mv dtm_rpa.qs dtm/dtm_rpa.qs')
        else:
            if on_hpc:
                write_job_script('.','dtm', 2, 64, True, 0, 'large-mem', pci_version, bin_dir)
                run_shell('mv dtm.qs dtm/dtm.qs')
        
        # Find even and odd directories with completed ci runs
        even_exists, odd_exists = False, False
        if os.path.isfile('even/CONF.RES') or os.path.isfile('even/CONFFINAL.RES'):
            even_exists = True
            if include_rpa: 
                run_shell('cp even/HFD.DAT dtm/HFD.DAT')
                run_shell('cp MBPT.INP dtm/MBPT.INP')
            run_shell('cp even/CONF.INP dtm/CONF.INP')
            run_shell('cp even/CONF.DET dtm/CONF.DET')
            run_shell('cp even/CONF.XIJ dtm/CONF.XIJ')
            run_shell('cp even/CONFSTR.RES dtm/CONFSTR.RES')
            run_shell('cp even/CONF.DAT dtm/CONF.DAT')
            run_shell('cp even/CONF.INT dtm/CONF.INT')
        if os.path.isfile('odd/CONF.RES') or os.path.isfile('odd/CONFFINAL.RES'):
            odd_exists = True
            run_shell('cp odd/CONF.INP dtm/CONF1.INP')
            run_shell('cp odd/CONF.DET dtm/CONF1.DET')
            run_shell('cp odd/CONF.XIJ dtm/CONF1.XIJ')
            run_shell('cp odd/CONFSTR.RES dtm/CONFSTR1.RES')
            
        if not even_exists and not odd_exists:
            print('ci directories could not be found')
            sys.exit()
                    
        # cd into new dtm directory and submit job script
        os.chdir(dtm_path)
        if include_rpa:
            with open('rpa.in', 'w') as f:
                f.write('2')
                
            run_shell('mpirun -n 1 ' + bin_dir + 'pdtm')
            
            write_dtm_in('TM',
                         '1 ' + str(num_levels) + ' 1 ' + str(num_levels),
                         ', '.join(key_list))
            
            if on_hpc:
                run_shell('sbatch dtm_rpa.qs')
        else:    
            if on_hpc:
                run_shell('sbatch dtm.qs')