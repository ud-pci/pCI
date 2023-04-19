import yaml
import re
import math
import sys
import os
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
    if code == 'conf' or code == 'ci':
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
        f.write('ulimit -s unlimited \n')
        f.write('CONF_MAX_BYTES_PER_CPU=$((SLURM_MEM_PER_CPU*1024*1024)) \n')
        f.write('export CONF_MAX_BYTES_PER_CPU \n')
        if code == 'basc' or code == 'ci':
            f.write('${UD_MPIRUN} basc \n')
        if code == 'conf' or code == 'ci':
            f.write('${UD_MPIRUN} conf \n')

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
        K_is = system['K_is']
        if K_is != 0:
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
            zero_plus1 = [0]*len(energies[0])
            zero_minus1 = [0]*len(energies[0])
            plus1_minus1 = [0]*len(energies[0])
            uncertainty1 = [0]*len(energies[0])
            uncertainty2 = [0]*len(energies[0])
            for num_energy in range(len(energies[0])): # runs over the total number of energy levels
                zero_plus1[num_energy] = (float(energies[3][num_energy]) - float(energies[2][num_energy])) / step_size
                zero_minus1[num_energy] = (float(energies[2][num_energy]) - float(energies[1][num_energy])) / step_size
                plus1_minus1[num_energy] = (float(energies[3][num_energy]) - float(energies[1][num_energy])) / (2*step_size)
                uncertainty1[num_energy] = plus1_minus1[num_energy] - zero_plus1[num_energy]
                uncertainty2[num_energy] = zero_minus1[num_energy] - plus1_minus1[num_energy]

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
                    f.write(',' + str(zero_plus1[num_energy]) + ',' + str(zero_minus1[num_energy]) + ',' + \
                                str(plus1_minus1[num_energy]) + ',' + str(uncertainty1[num_energy]) + ',' + \
                                str(uncertainty2[num_energy]) + ',')
                    f.write('\n')

            with open('energies.res','w') as f:
                conf_len = len(reduce(lambda x, y: x if len(x) > len(y) else y, main_confs))
                f.write('conf.'.rjust(conf_len, ' ') + ' ' + 'term' + '  ')
                for c in c_list:
                    f.write(str(c).rjust(len(energies[0][0]), ' ') + '  ')
                f.write('4pt deriv'.rjust(11, ' '))
                f.write('K_fs'.rjust(13, ' '))
                f.write('0-p1'.rjust(13, ' '))
                f.write('0-m1'.rjust(13, ' '))
                f.write('p1-m1'.rjust(13, ' '))
                f.write('uncert1'.rjust(13, ' '))
                f.write('uncert2'.rjust(13, ' '))
                f.write('\n')

                for num_energy in range(len(energies[0])): # runs over the total number of energy levels
                    f.write(main_confs[num_energy].rjust(conf_len, ' ') + '  ' + terms_list[num_energy].strip() + '  ')
                    for num_energy2 in range(len(energies)): # runs over the 5 c_is values
                        f.write(str(energies[num_energy2][num_energy]) + '  ')
                    f.write('{: .8f}'.format(derivatives[num_energy]) + '  ')
                    f.write('{: .8f}'.format(K_fs[num_energy]) + '  ')
                    f.write('{: .8f}'.format(zero_plus1[num_energy]) + '  ')
                    f.write('{: .8f}'.format(zero_minus1[num_energy]) + '  ')
                    f.write('{: .8f}'.format(plus1_minus1[num_energy]) + '  ')
                    f.write('{: .8f}'.format(uncertainty1[num_energy]) + '  ')
                    f.write('{: .8f}'.format(uncertainty2[num_energy]))
                    f.write('\n')

            f.close()
            print('analysis complete')
            sys.exit()
        else:
            confs = []
            terms = []
            energies = []
            
            dir_path = os.getcwd()
            for dir_name in ['EVEN','ODD']:
                os.chdir(dir_name)
                conf, term, energy_ev = parse_final_res()
                for nlevel in range(len(conf)):
                    confs.append(conf[nlevel])
                    terms.append(term[nlevel])
                    energies.append(float(energy_ev[nlevel]))
                os.chdir('../')
            for nlevel in range(len(confs)):
                print(confs[nlevel], terms[nlevel], energies[nlevel])

            # find ground state corresponding to lowest energy
            lowest_energy = max(energies)
            ind_lowest_energy = energies.index(max(energies))
            print('ground state has energy ' + str(lowest_energy) + ' corresponding to ' + confs[ind_lowest_energy] + terms[ind_lowest_energy])
            ht_to_cm = 219474.63 # hartree to cm-1
            energies_cm = [(-energy + lowest_energy) * ht_to_cm for energy in energies]

            confs_terms_energies = []
            for nlevel in range(len(confs)):
                confs_terms_energies.append((confs[nlevel], terms[nlevel], energies_cm[nlevel]))

            sorted_confs_terms_energies = sorted(confs_terms_energies, key=lambda x: x[2])

            with open('analysis.csv', 'w') as f:
                f.write('configuration, term, energy (cm-1), \n')
                for nlevel in range(len(confs)):
                    f.write(",".join([str(item) for item in sorted_confs_terms_energies[nlevel]]) + '\n')

            print('analysis complete')
            sys.exit()
    elif program == 'basc':
        ## TODO - check if CONFFINAL.RES exists - then ask user if they REALLY want to restart calculations
        # Run executables
        if system['K_is'] == 0:
            dir_path = os.getcwd()
            for dir_name in ['EVEN','ODD']:
                os.chdir(dir_name)
                write_job_script('basc.qs', program, 5, 64, True, 0, 'large-mem')
                run('sbatch basc.qs', shell=True)
                os.chdir('../')
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
    elif program == 'conf':
        ## TODO - check if CONFFINAL.RES exists - then ask user if they REALLY want to restart calculations
        # Run executables
        if system['K_is'] == 0:
            dir_path = os.getcwd()
            for dir_name in ['EVEN','ODD']:
                os.chdir(dir_name)
                write_job_script('conf.qs', program, 5, 64, True, 0, 'large-mem')
                run('sbatch conf.qs', shell=True)
                os.chdir('../')
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
    elif program == 'ci':
        ## TODO - check if CONFFINAL.RES exists - then ask user if they REALLY want to restart calculations
        # Run executables
        if system['K_is'] == 0:
            dir_path = os.getcwd()
            for dir_name in ['EVEN','ODD']:
                os.chdir(dir_name)
                write_job_script('ci.qs', program, 5, 64, True, 0, 'xlarge-mem')
                run('sbatch ci.qs', shell=True)
                os.chdir('../')
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
                write_job_script(program + '.qs', program, 5, 64, True, 0, 'xlarge-mem')
                run('sbatch ' + program + '.qs', shell=True)
                os.chdir('../')
    else:
        print(program + ' is not supported')
        sys.exit()