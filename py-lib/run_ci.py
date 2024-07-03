""" Isotope shifts

This script allows the user to automate isotope shift calculations from an already completed non-IS calculation.
The root/non-IS calculation should be in a directory called "0", and should contain the following files:

    * basis.sh - shell script to construct the basis set
    * HFD.INP - input file for hfd
    * BASS.INP - input file for bass
    * inf.aov, inf.vw (optional) - input files for all-order 
    * CONF.INP, c.in - input file for CI procedure
    * ci.qs - job script for parallel ci 

This python script has 3 main capabilities for isotope shift calculations:
1. Reconstruct basis in IS directories
2. Run CI procedure in IS directories
3. Analyze and compile the outputs from IS directories

If the user needs to start isotope shift calculations from scratch, e.g. do not have a "0" directory already completed,
the script 'basis.py" can be used instead.

"""
import re
import yaml
import sys
import os
from subprocess import run
from functools import reduce

def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

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
    # Check if config.yml file exists - if so, take K_is and C_is from there
    if os.path.isfile('config.yml'):
        config = read_yaml('config.yml')
        K_is = int(config['optional']['isotope_shifts']['K_is'])
        C_is = float(config['optional']['isotope_shifts']['C_is'])
    
    # If K_is in config.yml is 0, ask user for K_is
    if K_is == 0:
        K_is = int(input("Enter K_is (1 - FS, 2 - SMS, 3 - NMS, 4 - MS): "))
        C_is = float(input("Enter C_is: "))
    
    # Ask user for program to run
    program = input("Enter program you want to run (basis, ci, analysis): ")

    if program == 'analysis':  
        if K_is != 0:
            step_size = C_is/2
            c_list = [-C_is, -C_is/2, 0, C_is/2, C_is] 

            # Parse CONFFINAL.RES files in IS directories for energies 
            energies = []
            for c in c_list:
                dir_path = os.getcwd()
                dir_prefix = ''
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                dir_name = dir_prefix+'{:.5f}'.format(abs(c))
                if c == 0: dir_name = dir_prefix + '0'
                os.chdir(dir_name)
                main_confs, terms_list, energies_ev = parse_final_res()
                energies.append(energies_ev)
                os.chdir('../')

            # calculate 2-point and 4-point derivatives and isotope coefficients (K_fs, K_sms, K_nms, K_ms)
            points = [1, -8, 0, 8, -1] 
            radius = get_radius()
            au_to_si_conversion_factor_fs = 2.3497*10**-3  # 1 a.u. = 2.3497 x 10^-3 GHz/fm^2
            au_to_si_conversion_factor_ms = 3609.48        # 1 a.u. = 3609.48 GHz*amu
            derivatives = [] 
            K_is = [] 
            for num_energy in range(len(energies[0])): # runs over the total number of energy levels
                derivative  = 0
                for num_energy2 in range(len(energies)): # runs over the 5 c_is values
                    derivative += points[num_energy2] * float(energies[num_energy2][num_energy]) 
                derivative = derivative / (12*step_size) 
                derivatives.append(derivative)
                if K_is == 1:
                    coefficient = (5/6) * au_to_si_conversion_factor_fs * derivative / (radius**2)
                else:    
                    coefficient = au_to_si_conversion_factor_ms * derivative
                
                K_is.append(coefficient)

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

            filename = 'is' + C_is + '.out'
            with open(filename,'w') as f:
                conf_len = len(reduce(lambda x, y: x if len(x) > len(y) else y, main_confs))
                f.write('conf.'.rjust(conf_len, ' ') + ' ' + 'term' + '  ')
                for c in c_list:
                    f.write(str(c).rjust(len(energies[0][0]), ' ') + '  ')
                f.write('4pt deriv'.rjust(11, ' '))
                if K_is == 1: 
                    f.write('K_fs'.rjust(13, ' '))
                elif K_is == 2:
                    f.write('K_sms'.rjust(13, ' '))
                elif K_is == 3:
                    f.write('K_nms'.rjust(13, ' '))
                elif K_is == 4:
                    f.write('K_ms'.rjust(13, ' '))
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
                    f.write('{: .8f}'.format(K_is[num_energy]) + '  ')
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

    elif program == 'ci':
        # Run executables
        if K_is == 0:
            print('nothing to do')
            sys.exit()
        else:
            c_list = [-C_is,-C_is/2,C_is/2,C_is]
            
            # Check if directories already exist, else exit and have them run basis first
            dir_exists = True
            for c in c_list:
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                if not os.path.isdir(dir_prefix + '{:.5f}'.format(abs(c))):
                    print(dir_prefix + '{:.5f}'.format(abs(c)) + ' does not exist')
                    dir_exists = False
            if not dir_exists:
                print('Please re-run with "basis" option to generate basis for IS directories first')
                sys.exit()
            
            for c in c_list:
                dir_path = os.getcwd()
                dir_prefix = ''
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                dir_name = dir_prefix+'{:.5f}'.format(abs(c))
                print(dir_name)
                os.chdir(dir_name)
                run('sbatch ci.qs', shell=True)
                os.chdir('../')  
        
    elif program == 'basis':
        # Creating basis sets for isotopes shifts from an already completed /0 directory
        if K_is == 0:
            print('nothing to do')
            sys.exit()
        else:
            c_list = [-C_is,-C_is/2,C_is/2,C_is]
            
            # Check if directories already exist and prompt user to delete them
            dir_exists = False
            for c in c_list:
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                if os.path.isdir(dir_prefix + '{:.5f}'.format(abs(c))):
                    print(dir_prefix + '{:.5f}'.format(abs(c)) + ' already exists')
                    dir_exists = True
            if dir_exists:
                del_dirs =  eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input('Would you like to delete these directory? ')))))
                if del_dirs:
                    for c in c_list:
                        if c < 0:
                            dir_prefix = 'minus'
                        elif c > 0:
                            dir_prefix = 'plus'
                        run('rm -r ' + dir_prefix + '{:.5f}'.format(abs(c)), shell=True)
                else:
                    print('Please move directories before continuing..')
                    sys.exit()
            
            for c in c_list:
                dir_path = os.getcwd()
                dir_prefix = ''
                if c < 0:
                    dir_prefix = 'minus'
                elif c > 0:
                    dir_prefix = 'plus'
                dir_name = dir_prefix+'{:.5f}'.format(abs(c))
                print(dir_name)
                run("pwd", shell=True)
                run("cp -r 0 " + dir_name, shell=True)
                os.chdir(dir_name)
                if abs(c) > 0.0:
                    run("find . -type f -name \"*.INP\" -exec sed -i 's/K_is= 0/K_is= " + '{:.5f}'.format(K_is) + "/g' \{\} \;", shell=True)
                    run("find . -type f -name \"*.INP\" -exec sed -i 's/C_is= 0/C_is= " + '{:.5f}'.format(c) + "/g' \{\} \;", shell=True)
                    if K_is > 1:
                        run("find . -type f -name \"*.INP\" -exec sed -i 's/Klow= 0/Klow= 2" + "/g' \{\} \;", shell=True)
                run('./basis.sh', shell=True)
                os.chdir('../')  
              
    else:
        print(program + ' is not supported')
        sys.exit()