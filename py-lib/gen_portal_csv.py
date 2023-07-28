##
# 1. create energy csv file (NIST if it exists, our data if NIST does not exist)
    # - e.g. Sr1_Energies.csv
# 2. create matrix element file
    # - read outputs from dtm (all order and MBPT) and produce csv file [Sr1_Matrix_Elements.csv]
    # - all theory data
    # - create uncertainty from abs(all-order - MBPT)
# 3. option to hide for display G levels

import re
import sys
import pandas as pd

def parse_final_res(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    energies_au = []
    energies_cm = []
    main_confs = []
    terms_list = []
    for line in lines[1:]:
        ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        confs = [conf.replace(' ', '.') for conf in confs]
        Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]
      
        energies_au.append(float(nums[0]))
        energies_cm.append(float(nums[1]))
        main_confs.append(confs[0])
        terms_list.append(terms[0].replace(' ', ''))

    return main_confs, terms_list, energies_au, energies_cm

def edit_vipul(filename, resfilename):

    # Read vipul's data file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Parse CONFFINAL.RES file to obtain energies in a.u.
    _, _, energies_au, _ = parse_final_res(resfilename)

    lowest_energy = max(energies_au)
    ht_to_cm = 219474.63 # hartree to cm-1
    energies_cm = [str(round((-energy + lowest_energy) * ht_to_cm, 1)) for energy in energies_au]

    f = open(filename.split('.')[0] + '_fixed.txt', 'w')

    f.write(lines[0][:-1] + '  Level (a.u.) \n')
    for line in lines[1:]:
        NIST_config = line.split()[0]
        NIST_term = line.split()[1]
        NIST_J = line.split()[2]
        NIST_energy = line.split()[3]
        # Manually adjust odd parity energies until Vipul's update
        odd_energy = 14317.507
        if NIST_config[0].isnumeric() and filename == 'Odd.txt': 
            NIST_energy = str(round(float(NIST_energy) + odd_energy, 3))
        NIST_uncertainty = line.split()[4]
        theory_config = line.split()[5]
        theory_config2 = line.split()[6]
        theory_term = line.split()[7]
        theory_J = line.split()[8]
        theory_energy = line.split()[9]
        theory_delta_cm = line.split()[10]
        theory_delta = line.split()[11]
        if theory_energy in energies_cm:
            theory_energy_au = str(energies_au[energies_cm.index(theory_energy)])
            min_energy = 0.61723101
            if theory_config[0].isnumeric() and filename == 'Odd.txt': 
                theory_energy = str(round((-float(theory_energy_au) + min_energy) * ht_to_cm, 1))
        else:
            theory_energy_au = '-'
            
        #new_line = line[:-1] + str(theory_energy_au).rjust(14, ' ') + '\n'
        new_line = NIST_config.rjust(7, ' ') + NIST_term.rjust(7, ' ') + NIST_J.rjust(4, ' ') + \
                    NIST_energy.rjust(15, ' ') + NIST_uncertainty.rjust(19, ' ') + theory_config.rjust(10, ' ') + \
                    theory_config2.rjust(10, ' ') + theory_term.rjust(7, ' ') + theory_J.rjust(4, ' ') + \
                    theory_energy.rjust(15, ' ') + theory_delta_cm.rjust(7, ' ') + theory_delta.rjust(8, ' ') + \
                    theory_energy_au.rjust(14, ' ') + '\n'

        f.write(new_line)

    f.close()

        
    return
    

def read_vipul():
    f = open('Even_fixed.txt', 'r')
    lines = f.readlines()
    f.close()

    f = open('Odd_fixed.txt', 'r')
    lines = lines + f.readlines()
    f.close()
    
    # name columns of dataframe
    df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J',
                               'energy', 'energy_uncertainty', 'is_from_theory'])
    
    min_energy = 0
    ht_to_cm = 219474.63 # hartree to cm-1

    # find minimum energy in a.u. to determine ground state
    for line in lines:
        theory_energy_au = line.split()[12]
        if theory_energy_au.replace('.', '').isnumeric():
            if float(theory_energy_au) > min_energy:
                min_energy = float(theory_energy_au)
    if not min_energy:
        print('minimum energy could not be found')
        sys.exit()
        
    for line in lines:
        NIST_config = line.split()[0]
        NIST_term = line.split()[1]
        NIST_J = line.split()[2]
        NIST_energy = line.split()[3]
        NIST_uncertainty = line.split()[4]
        theory_config = line.split()[5]
        theory_config2 = line.split()[6]
        theory_term = line.split()[7]
        theory_J = line.split()[8]
        theory_energy = line.split()[9]
        theory_delta_cm = line.split()[10]
        theory_delta = line.split()[11]
        theory_energy_au = line.split()[12]

        # select relevant data for portal database 
        config = ''
        if NIST_config[0].isnumeric():
            config = NIST_config
            term = NIST_term
            J = NIST_J
            energy = NIST_energy
            energy_uncertainty = NIST_uncertainty
            is_from_theory = False
        elif theory_config[0].isnumeric():
            config = theory_config
            term = theory_term
            J = theory_J
            energy = round((-float(theory_energy_au) + min_energy) * ht_to_cm, 1)
            energy_uncertainty = '-'
            is_from_theory = True
        else:
            continue

        row = {'state_configuration': config, 'state_term': term, 'state_J': J,
                'energy': energy, 'energy_uncertainty': energy_uncertainty, 'is_from_theory': is_from_theory}

        df.loc[len(df.index)] = row
    
    df.to_csv('test.csv', index=False)
    print('test.csv' + ' has been written')

    return df

def write_energy_csv(element):
    '''
    This function writes the energy csv file
    '''
    filename = element + '_Energies.csv'

    df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J',
                               'energy', 'energy_uncertainty', 'is_from_theory'])
    
    confs, terms, energies_au, energies_cm = merge_res('CONFFINALeven.RES', 'CONFFINALodd.RES')
    confs_MBPT, terms_MBPT, energies_au_MBPT, energies_cm_MBPT = merge_res('CONFFINALevenMBPT.RES', 'CONFFINALoddMBPT.RES')

    # Determine minimum energies
    min_energy = max(energies_au)
    min_energy_MBPT = max(energies_au_MBPT)

    # Calculate energies in cm-1 from ground state energy
    ht_to_cm = 219474.63 # hartree to cm-1
    energies_cm = [(-energy + min_energy) * ht_to_cm for energy in energies_au]
    energies_cm_MBPT = [(-energy + min_energy_MBPT) * ht_to_cm for energy in energies_au_MBPT]

    # Calculate uncertainties
    uncertainties = calc_uncertainties(energies_cm, energies_cm_MBPT)

    confs_terms_energies = []
    for nlevel in range(len(confs)):
        confs_terms_energies.append((confs[nlevel], terms[nlevel], energies_cm[nlevel], uncertainties[nlevel]))

    sorted_confs_terms_energies = sorted(confs_terms_energies, key=lambda x: x[2])
    for state in sorted_confs_terms_energies:
        row = {'state_configuration': state[0], 'state_term': state[1][:2], 'state_J': state[1][-1],
               'energy': "{:.3f}".format(state[2]), 'energy_uncertainty': state[3], 'is_from_theory': 'TRUE'}
        df.loc[len(df.index)] = row

    df['temp_conf'] = df['state_configuration'].str.replace('.','')
        
    # TODO - read in Vipul's table and update values
    vipul_df = read_vipul()
#
    #merged_df = df.merge(vipul_df, how='outer', on=['NIST_conf','state_term', 'state_J'])
    #print(merged_df.to_string())

    ###

    df.to_csv(filename, index=False)
    print(filename + ' has been written')

    return

def write_matrix_csv(element):
    '''
    This function writes the matrix element csv file
    '''
    filename = element + '_Matrix_Elements_Theory.csv'
    f = open('E1.RES', 'r') 
    lines = f.readlines()
    f.close()

    df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty'])

    for line in lines:
        matrix_element = re.findall(r'\<.*?\>', line)[0]
        Tk = re.findall(r'\|.*?\|', matrix_element)[0].replace('|', '').replace(' ', '')
        state1 = re.findall(r'\|.*?\>', matrix_element)[0].replace(Tk, '').replace('|', '').replace('>', '')
        state2 = re.findall(r'\<.*?\|', matrix_element)[0].replace('<', '').replace('|', '')
        conf1 = '.'.join(state1.split()[:-1])
        conf2 = '.'.join(state2.split()[:-1])
        term1 = state1.split()[-1][0:2]
        term2 = state2.split()[-1][0:2]
        J1 = state1.split()[-1][-1]
        J2 = state2.split()[-1][-1]

        matrix_element_value = re.findall("\d+\.\d+", line)[:1]
        matrix_element_value = matrix_element_value[0] if matrix_element_value else None

        row = {'state_one_configuration': conf1, 'state_one_term': term1, 'state_one_J': J1,
               'state_two_configuration': conf2, 'state_two_term': term2, 'state_two_J': J2,
               'matrix_element': matrix_element_value, 'matrix_element_uncertainty': None}
        df.loc[len(df.index)] = row
    
    df = df.tail(-1)

    df.to_csv(filename, index=False)
    print(filename + ' has been written')

    return

def calc_uncertainties(energies_cm1, energies_cm2):
    '''
    This function parses calculates energy uncertainties between them
    '''

    # Check if number of levels are the same
    if len(energies_cm1) != len(energies_cm2):
        print('Number of levels are not the same')
        sys.exit()
    nlv = len(energies_cm1)

    # Calculate uncertainties
    uncertainties = []
    for n in range(nlv):
        if energies_cm1[n] > 0 or energies_cm2[n] > 0:
            uncertainty = round(abs(energies_cm2[n] - energies_cm1[n])/(abs(energies_cm1[n] + energies_cm2[n])/2),4)
        else:
            uncertainty = 0
        uncertainties.append(uncertainty)
    
    return uncertainties

def reorder_levels(confs1, terms1, confs2, terms2, energies_au2, energies_cm2):
    '''
    This function swaps energies of any re-ordered levels
    '''

    # Ensure the number of energy levels are the same
    nlv = len(confs1)
    nlv2 = len(confs2)
    if nlv != nlv2:
        print('number of levels are not the same')
        sys.exit()

    # Check and see if any levels are re-ordered
    if confs1 == confs2 and terms1 == terms2:
        print('levels are the same')
    else:
        # Fix any level re-orderings
        for n in range(nlv):
            if confs1[n] != confs2[n]:
                if confs1[n] == confs2[n-1] and terms1[n] == terms2[n-1]:
                    conf = confs2[n]
                    confs2[n] = confs2[n-1]
                    confs2[n-1] = conf

                    term = terms2[n]
                    terms2[n] = terms2[n-1]
                    terms2[n-1] = term

                    energy = energies_au2[n]
                    energies_au2[n] = energies_au2[n-1]
                    energies_au2[n-1] = energy

                    energy_cm = energies_cm2[n]
                    energies_cm2[n] = energies_cm2[n-1]
                    energies_cm2[n-1] = energy_cm

        # Final check of different levels
        for n in range(nlv):
            if confs1[n] != confs2[n] or terms1[n] != terms2[n]: 
                if terms1[n] != terms2[n-1]:
                    print('Different configurations that are not re-ordered found for Level #' + str(n), confs1[n], terms1[n], '->', confs2[n], terms2[n])

    return confs1, terms1, energies_au2, energies_cm2


def merge_res(res1, res2):
    '''
    This function merges the results of two CONFFINAL.RES files and returns relevant data

    confs - configurations
    terms - terms ^(2S+1)L_J
    energies_au - energies in a.u.
    energies_cm - energies in cm-1
    '''
    confs1, terms1, energies_au1, energies_cm1 = parse_final_res(res1)
    confs2, terms2, energies_au2, energies_cm2 = parse_final_res(res2)

    confs = confs1 + confs2
    terms = terms1 + terms2
    energies_au = energies_au1 + energies_au2
    energies_cm = energies_cm1 + energies_cm2

    return confs, terms, energies_au, energies_cm



if __name__ == "__main__":
    element = 'Sr1'

    confs, terms, energies_au, energies_cm = merge_res('CONFFINALeven.RES', 'CONFFINALodd.RES')
    confs_MBPT, terms_MBPT, energies_au_MBPT, energies_cm_MBPT = merge_res('CONFFINALevenMBPT.RES', 'CONFFINALoddMBPT.RES')

    #reorder_levels(confs, terms, confs_MBPT, terms_MBPT, energies_au_MBPT, energies_cm_MBPT)
    #even_confs, even_terms, even_energies, even_energies_MBPT, even_uncertainties = calc_theory_uncertainties('CONFFINALeven.RES', 'CONFFINALevenMBPT.RES')

    #edit_vipul('Even.txt', 'CONFFINALeven.RES')
    #edit_vipul('Odd.txt', 'CONFFINALodd.RES')
    write_energy_csv(element)
    #write_matrix_csv(element)