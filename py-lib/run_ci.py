import yaml
import re
import math
import sys
import os
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT, run
from functools import reduce


'''
This python script has 3 main capabilities for isotope shift calculations:
1. Run basc in IS directories
2. Run conf in IS directories
3. Analyze and compile the outputs from IS directories

'''
def read_yaml(filename):
    """ Reads yaml input file and returns contents """ 

    isotope = 0
    val_aov = []
    # read yaml file with inputs
    with open(filename,'r') as f:
        config = yaml.safe_load(f)

        # Check to see if isotope is specified. If not, set to 0
        try:
            isotope = config['isotope']
        except KeyError as e:
            config[e.args[0]] = 0
        
        # Check to see if val_aov is specified. If not, set to []
        try:
            val_aov = config['val_aov']
        except KeyError as e:
            config[e.args[0]] = []

        # Check to see if energies are specified. If not, set to []
        try:
            energies = config['energies']
        except KeyError as e:
            config[e.args[0]] = []

        # Check to see if isotope shift key is specified. If not, set to 0
        try:
            isotope = config['K_is']
        except KeyError as e:
            config[e.args[0]] = 0
        try:
            C_is = config['C_is']
        except KeyError as e:
            config[e.args[0]] = 0 

        # Check whether or not to run all-order codes after basis set construction
        try:
            run_ao_codes = config['run_ao_codes']
        except KeyError as e:
            config[e.args[0]] = 0 

    return config

def write_job_script(filename, code, num_nodes, num_procs_per_node, exclusive, mem, partition):
    if code == 'conf':
        with open('c.in', 'w') as f:
            f.write('2, 2, 0, 0, 1')
        f.close()

    with open(filename, 'w') as f:
        f.write('#!/bin/bash -l \n')
        f.write(' \n')
        f.write('#SBATCH --nodes=' + str(num_nodes) + ' \n')
        f.write('#SBATCH --tasks-per-node=' + str(num_procs_per_node) + ' \n')
        if exclusive: 
            f.write('#SBATCH --exclusive=user \n')
        f.write('#SBATCH --cpus-per-task=1 \n')
        f.write('#SBATCH --mem=' + str(mem) + ' \n')
        f.write('#SBATCH --job-name=' + code + ' \n')
        f.write('#SBATCH --partition=' + partition + ' \n')
        f.write('#SBATCH --time=05-00:00:00 \n')
        f.write('#SBATCH --export=NONE \n')
        f.write(' \n')
        f.write('vpkg_require pci \n')
        f.write(' \n')
        f.write('UD_PREFER_MEM_PER_CPU=YES \n')
        f.write('UD_REQUIRE_MEM_PER_CPU=YES \n')
        f.write(' \n')
        f.write('. /opt/shared/slurm/templates/libexec/openmpi.sh \n')
        f.write(' \n')
        if code == 'conf':
            f.write('ulimit -s unlimited \n')
            f.write('CONF_MAX_BYTES_PER_CPU=$((SLURM_MEM_PER_CPU*1024*1024)) \n')
            f.write('export CONF_MAX_BYTES_PER_CPU \n')
        f.write('${UD_MPIRUN} ' + code + ' \n')
        f.write('mpi_rc=$? \n')
        f.write(' \n')
        f.write('exit $mpi_rc \n')
    f.close()

def parse_final_res():
    f = open('CONFFINAL.RES', 'r')
    lines = f.readlines()
    f.close()

    energies_ev = []
    main_confs = []
    terms_list = []
    for line in lines[1:]:
        ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]
        energies_ev.append(nums[0])
        main_confs.append(confs[0])
        terms_list.append(terms[0])

    return main_confs, terms_list, energies_ev

def get_radius():
    f = open('0/HFD.RES', 'r') # TODO - default 0 directory, change to be dynamic and check all radii in HFD.RES
    lines = f.readlines()
    f.close()

    for line in lines:
        if 'R1' in line:
            radius = float(line.split('=')[-1])
    
    return radius

if __name__ == "__main__":
    # Read yaml file for system configurations
    yml_file = input("Input yml-file: ")
    system = read_yaml(yml_file)

    program = input("Enter program you want to run: ")

    if program == 'analysis':  
        C_is = system['C_is']
        step_size = C_is/2
        c_list = [-C_is, -C_is/2, 0, C_is/2, C_is] 

        energies = []
        for c in c_list:
            dir_path = os.getcwd()
            dir_prefix = ''
            if c < 0:
                dir_prefix = 'minus'
            elif c > 0:
                dir_prefix = 'plus'
            dir_name = dir_prefix+str(abs(c))
            os.chdir(dir_name)
            main_confs, terms_list, energies_ev = parse_final_res()
            energies.append(energies_ev)
            os.chdir('../')
            
        
        # calculate 4-point derivatives and field shift coefficients K_fs
        points = [1, -8, 0, 8, -1] 
        radius = get_radius()
        au_to_si_conversion_factor_fs = 2.3497*10**-3  # 1 a.u. = 2.3497 x 10^-3 GHz/fm^2
        derivatives = [] 
        K_fs = [] 
        for num_energy in range(len(energies[0])): # runs over the total number of energy levels
            derivative = 0
            for num_energy2 in range(len(energies)): # runs over the 5 c_is values
                derivative += points[num_energy2] * float(energies[num_energy2][num_energy]) 
            derivative = derivative / 0.12
            derivatives.append(derivative)
            coefficient = (5/6) * au_to_si_conversion_factor_fs * derivative / (radius**2)
            K_fs.append(coefficient)
 

        # calculate uncertainties
        zero_plus1 = []
        zero_minus1 = []
        plus1_minus1 = []
        uncertainty1 = []
        uncertainty2 = []
        for num_energy in range(len(energies[0])): # runs over the total number of energy levels
            zero_plus1 = (float(energies[3][num_energy]) - float(energies[2][num_energy])) / step_size
            zero_minus1 = (float(energies[2][num_energy]) - float(energies[1][num_energy])) / step_size
            plus1_minus1 = (float(energies[3][num_energy]) - float(energies[1][num_energy])) / (2*step_size)
            uncertainty1 = plus1_minus1 - zero_plus1
            uncertainty2 = zero_minus1 - plus1_minus1

        with open('energies.csv','w') as f:
            conf_len = len(reduce(lambda x, y: x if len(x) > len(y) else y, main_confs))
            f.write('conf.' + ',' + 'term' + ',')
            for c in c_list:
                f.write(str(c).rjust(len(energies[0][0]), ' ') + ',')
            f.write(' 4pt deriv,' + 'K_fs (GHz/fm2),' + '0-p1,' + '0-m1,' + 'p1-m1,' + 'uncertainty1,' + 'uncertainty2')
            f.write('\n')

            for num_energy in range(len(energies[0])): # runs over the total number of energy levels
                f.write(main_confs[num_energy].rjust(conf_len, ' ') + ',' + terms_list[num_energy].strip() + ',')
                for num_energy2 in range(len(energies)): # runs over the 5 c_is values
                    f.write(str(energies[num_energy2][num_energy]) + ',')
                f.write("%.8f" % derivatives[num_energy])
                f.write(',')
                f.write("%.2f" % K_fs[num_energy])
                f.write(',' + str(zero_plus1) + ',' + str(zero_minus1) + ',' + str(plus1_minus1) + ',' + str(uncertainty1) + ',' + str(uncertainty2) + ',')
                f.write('\n')

        with open('energies.res','w') as f:
            conf_len = len(reduce(lambda x, y: x if len(x) > len(y) else y, main_confs))
            f.write('conf.'.rjust(conf_len, ' ') + ' ' + 'term' + '  ')
            for c in c_list:
                f.write(str(c).rjust(len(energies[0][0]), ' ') + '  ')
            f.write(' 4pt deriv')
            f.write('\n')

            for num_energy in range(len(energies[0])): # runs over the total number of energy levels
                f.write(main_confs[num_energy].rjust(conf_len, ' ') + '  ' + terms_list[num_energy].strip() + '  ')
                for num_energy2 in range(len(energies)): # runs over the 5 c_is values
                    f.write(str(energies[num_energy2][num_energy]) + '  ')
                f.write("%.8f" % derivatives[num_energy])
                f.write('\n')
                
        f.close()
        print('analysis complete')
        sys.exit()

    if system['codes'] == 'ci':
        print('This feature has not been implemented yet.')
        sys.exit()
    elif system['codes'] == 'ci+all-order':
        # Run executables
        if system['K_is'] == 0:
            write_job_script(program + '.qs', program, 2, 64, True, 0, 'standard')
            run('sbatch ' + program + '.qs', shell=True)
        else:
            C_is = system['C_is']
            c_list = [-C_is,-C_is/2,0,C_is/2,C_is]
            for c in c_list:
                dir_path = os.getcwd()
                dir_prefix = ''
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                dir_name = dir_prefix+str(abs(c))
                os.chdir(dir_name)
                write_job_script(program + '.qs', program, 5, 64, True, 0, 'large-mem')
                run('sbatch ' + program + '.qs', shell=True)
                os.chdir('../')
    else:
        print(system['codes'] + 'not supported')
        sys.exit()