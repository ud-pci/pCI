##
# 1. create energy csv file (NIST if it exists, our data if NIST does not exist)
    # - e.g. Sr1_Energies.csv
# 2. create matrix element file
    # - read outputs from dtm (all order and MBPT) and produce csv file [Sr1_Matrix_Elements.csv]
    # - all theory data
    # - create uncertainty from abs(all-order - MBPT)
# 3. option to hide for display G levels

import re
import pandas as pd

def parse_final_res(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    energies_ev = []
    main_confs = []
    terms_list = []
    for line in lines[1:]:
        ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        confs = [conf.replace(' ', '.') for conf in confs]
        Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]
      
        energies_ev.append(float(nums[0]))
        main_confs.append(confs[0])
        terms_list.append(terms[0].replace(' ', ''))

    return main_confs, terms_list, energies_ev


def read_vipul():
    f = open('Even.txt', 'r')
    lines = f.readlines()
    f.close()

    f = open('Odd.txt', 'r')
    #lines = lines + f.readlines()
    f.close()
    
    df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J',
                               'energy', 'energy_uncertainty', 'is_from_theory'])
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

        # TODO - select data 
        config = ''
        if NIST_config[0].isnumeric():
            config = NIST_config
            term = NIST_term
            J = NIST_J
            energy = NIST_energy
            energy_uncertainty = NIST_uncertainty
            is_from_theory = False
        else:
            config = theory_config
            term = theory_term
            J = theory_J
            energy = theory_energy
            energy_uncertainty = '-'
            is_from_theory = True

        row = {'state_configuration': config, 'state_term': term, 'state_J': J,
                'energy': energy, 'energy_uncertainty': energy_uncertainty, 'is_from_theory': is_from_theory}

        df.loc[len(df.index)] = row
    
    df = df.tail(-1)

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
    
    even_confs, even_terms, even_energies = parse_final_res('CONFFINALeven.RES')
    odd_confs, odd_terms, odd_energies = parse_final_res('CONFFINALodd.RES')

    confs = even_confs + odd_confs
    terms = even_terms + odd_terms
    energies = even_energies + odd_energies

    lowest_energy = max(energies)
    ind_lowest_energy = energies.index(lowest_energy)
    #print('ground state has energy ' + str(lowest_energy) + ' corresponding to ' + confs[ind_lowest_energy] + ' ' + terms[ind_lowest_energy])
    ht_to_cm = 219474.63 # hartree to cm-1
    energies_cm = [(-energy + lowest_energy) * ht_to_cm for energy in energies]

    confs_terms_energies = []
    for nlevel in range(len(confs)):
        confs_terms_energies.append((confs[nlevel], terms[nlevel], energies_cm[nlevel]))

    sorted_confs_terms_energies = sorted(confs_terms_energies, key=lambda x: x[2])
    for state in sorted_confs_terms_energies:
        row = {'state_configuration': state[0], 'state_term': state[1][:2], 'state_J': state[1][-1],
               'energy': "{:.3f}".format(state[2]), 'energy_uncertainty': None, 'is_from_theory': 'TRUE'}
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

if __name__ == "__main__":
    element = 'Sr1'
    write_energy_csv(element)
    write_matrix_csv(element)