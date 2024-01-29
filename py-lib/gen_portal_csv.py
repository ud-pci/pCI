import re
import sys
import os
import pandas as pd
from fractions import Fraction
from UDRead import *
from parse_asd import *
from get_atomic_term import *

def parse_final_res(filename):
    try:
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
    except:
        print(filename + ' not found')
        return [], [], [], []

    energies_au = []
    energies_cm = []
    main_confs = []
    terms_list = []
    ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
    Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    for line in lines[1:]:
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        confs = [conf.replace(' ', '.') for conf in confs]
        
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]
      
        energies_au.append(float(nums[0]))
        energies_cm.append(float(nums[1]))
        main_confs.append(confs[0])
        terms_list.append(terms[0].replace(' ', ''))

    # Separate terms
    i = 0
    for term in terms_list:
        term = term[0:2] + ',' + term[2:]
        terms_list[i] = term
        i += 1

    return main_confs, terms_list, energies_au, energies_cm

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
            uncertainty = round(abs(energies_cm2[n] - energies_cm1[n]))
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

def create_mapping():
    '''
    This function reads Vipul's energy level table and creates the mapping between experimental and theory data
    Data to map: Config, Term, J, Energy(cm-1)
    '''
    filepath = "DATA_Output/"+name+"_Even.txt"
    f = open(filepath, 'r')
    lines = f.readlines()
    f.close()

    filepath = "DATA_Output/"+name+"_Odd.txt"
    f = open(filepath, 'r')
    lines = lines + f.readlines()
    f.close()
    
    mapping = []
    for line in lines[1:]:
        NIST_config = line.split()[0]
        NIST_term = line.split()[1]
        NIST_J = line.split()[2]
        NIST_energy = line.split()[3]
        NIST_uncertainty = line.split()[4]
        final_config = line.split()[5]
        theory_config = line.split()[6]
        corrected_config = line.split()[5]
        theory_config2 = line.split()[7]
        theory_term = line.split()[8]
        try:
            theory_J = str(Fraction(line.split()[9]))
        except:
            theory_J = line.split()[9]
        theory_energy_cm = line.split()[10]
        theory_uncertainty = line.split()[11]
        theory_energy_au = line.split()[12]
        #theory_delta_cm = line.split()[13]
        #theory_delta = line.split()[14]

        # select relevant data for portal database 
        if NIST_config != 'Config':
            mapping.append([[NIST_config, NIST_term, NIST_J, NIST_energy, NIST_uncertainty],
                    [theory_config, theory_term, theory_J, theory_energy_cm, theory_uncertainty, corrected_config, theory_energy_au]])
    
    return mapping

def write_new_conf_res(name, filepath):
    '''
    This function creates a new CONF.RES file with theory uncertainties
    '''
    
    matrix_file_exists = True
    # Check if raw files exist
    if not os.path.exists(filepath + 'CONFFINALeven.RES'):
        print('CONFFINALeven.RES not found in', filepath)
        sys.exit()
    if not os.path.exists(filepath + 'CONFFINALodd.RES'):
        print('CONFFINALodd.RES not found in', filepath)
        sys.exit()
    if not os.path.exists(filepath + 'CONFFINALevenMBPT.RES'):
        print('CONFFINALevenMBPT.RES not found in', filepath)
    if not os.path.exists(filepath + 'CONFFINALoddMBPT.RES'):
        print('CONFFINALoddMBPT.RES not found in', filepath)
    if not os.path.exists(filepath + 'E1.RES'):
        print('E1.RES not found in', filepath)
        matrix_file_exists = False
    if not os.path.exists(filepath + 'E1MBPT.RES'):
        print('E1MBPT.RES not found in', filepath)
        matrix_file_exists = False

    # Read CONF.RES files
    _, _, energies_au_even, energies_cm_even = parse_final_res(filepath + 'CONFFINALeven.RES')
    _, _, energies_au_even_MBPT, _ = parse_final_res(filepath + 'CONFFINALevenMBPT.RES')

    _, _, energies_au_odd, energies_cm_odd = parse_final_res(filepath + 'CONFFINALodd.RES')
    _, _, energies_au_odd_MBPT, _ = parse_final_res(filepath + 'CONFFINALoddMBPT.RES')

    # Merge even and odd parity CONF.RES files
    confs, terms, energies_au, energies_cm = merge_res(filepath + 'CONFFINALeven.RES', filepath + 'CONFFINALodd.RES')
    _, _, energies_au_MBPT, energies_cm_MBPT = merge_res(filepath + 'CONFFINALevenMBPT.RES', filepath + 'CONFFINALoddMBPT.RES')

    # Determine minimum energies
    min_energy = max(energies_au)
    try:
        min_energy_MBPT = max(energies_au_MBPT)
    except:
        min_energy_MBPT = 0

    # Calculate energies in cm-1 from ground state energy
    ht_to_cm = 219474.63 # hartree to cm-1
    energies_cm_even = [(-energy + min_energy) * ht_to_cm for energy in energies_au_even]
    if min_energy_MBPT: energies_cm_MBPT_even = [(-energy + min_energy_MBPT) * ht_to_cm for energy in energies_au_even_MBPT]
    energies_cm_odd = [(-energy + min_energy) * ht_to_cm for energy in energies_au_odd]
    if min_energy_MBPT: energies_cm_MBPT_odd = [(-energy + min_energy_MBPT) * ht_to_cm for energy in energies_au_odd_MBPT]
    energies_cm = [(-energy + min_energy) * ht_to_cm for energy in energies_au]
    if min_energy_MBPT: energies_cm_MBPT = [(-energy + min_energy_MBPT) * ht_to_cm for energy in energies_au_MBPT]

    # Calculate uncertainties
    uncertainties = []
    if min_energy_MBPT:
        uncertainties_even = calc_uncertainties(energies_cm_even, energies_cm_MBPT_even)
        uncertainties_odd = calc_uncertainties(energies_cm_odd, energies_cm_MBPT_odd)
        uncertainties = calc_uncertainties(energies_cm, energies_cm_MBPT)
    else:
        uncertainties_even = [0]*len(energies_cm_even)
        uncertainties_odd = [0]*len(energies_cm_odd)

    # Write csv-formatted CONF.RES files with uncertainties
    convert_res_to_csv(filepath + 'CONFFINALeven.RES', uncertainties_even, name)
    convert_res_to_csv(filepath + 'CONFFINALodd.RES', uncertainties_odd, name)

    # Determine energy shift between odd and even parity lowest energy levels
    energy_shift = abs(ht_to_cm * (energies_au_even[0] - energies_au_odd[0]))
    
    # Determine parity of ground state level
    gs_parity = ''
    if energies_au_even[0] > energies_au_odd[0]:
        gs_parity = 'even'
    else:
        gs_parity = 'odd'

    return confs, terms, energies_cm, uncertainties, energy_shift, gs_parity, matrix_file_exists

def convert_type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    s = s.replace(" ", "")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        return s 
    
def convert_res_to_csv(filename, uncertainties, name):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    if 'odd' in filename:
        parity = 'Odd'
    if 'even' in filename:
        parity = 'Even'
    csvfile = "DATA_Filtered/UD/"+name+'_UD_' + parity + '.csv'

    os.makedirs(os.path.dirname(csvfile), exist_ok=True)
    f = open(csvfile, 'w')
    f.write('n, conf, term, E_n (a.u.), DEL (cm^-1), S, L, J, gf, conf%, converged, conf2, conf2%, uncertainty \n')

    num_cols = 13
    i = 0
    for line in lines[1:]:
        # reformatting
        newline = re.sub('\s{3,}', '  ', line).replace('  ', ',') 
        if newline[0] == ',':
            newline = newline[1:]
        newline = newline.lstrip()[:-1]

        if newline.count(',') < num_cols-1:
            num_extra_cols = num_cols - newline.count(',') - 1
            newline += num_extra_cols*',' + str(uncertainties[i]) + '\n'
        i += 1
        string_to_list = newline.split(",")
        modified_list = [str(Convert_Type(i)) for i in string_to_list]
        newline = ",".join(modified_list) + '\n'

        f.write(newline)
    f.close()
    print(csvfile + ' has been written')

    return

def find_parity(configuration):
    '''
    This function finds parity of specified configuration
    '''
    p = 0
    parity = ''
    ldict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
    
    orbitals = configuration.split('.')
    for orbital in orbitals:
        nq = re.findall('[0-9]+', orbital)
        if len(nq) <= 1: 
            q = 1
        else:
            q = int(nq[1])
        
        l_str = re.findall('[spdfghi]+', orbital)[0]
        l = ldict[l_str]
        p += l*q
    
    if p%2 == 0:
        parity = 'even'
    else:
        parity= 'odd'
        
    return parity

def write_energy_csv(name, mapping, NIST_shift, theory_shift, gs_parity):
    '''
    This function writes the energy csv file
    '''
    filename = name + '_Energies.csv'

    # select columns for portal dataframe
    portal_df = pd.DataFrame(columns = ['state_configuration', 'state_term', 'state_J', 'energy', 
                                     'energy_uncertainty', 'is_from_theory'])

    for level in mapping:
        # if experimental data does not exist, use theory values
        if level[0][3] == '-':
            is_from_theory = True
            state_config = level[1][5]
            state_term = level[1][1]
            state_J = level[1][2]
            if find_parity(state_config) != gs_parity:
                state_energy = "{:.1f}".format(float(level[1][3]) + float(theory_shift))
            else:
                state_energy = level[1][3]
            state_uncertainty = level[1][4]
        else:
            is_from_theory = False
            state_config = level[0][0]
            state_term = level[0][1]
            state_J = level[0][2]
            if find_parity(state_config) != gs_parity:
                state_energy = "{:.3f}".format(float(level[0][3]) + float(NIST_shift))
            else:
                state_energy = level[0][3]
            state_uncertainty = level[0][4]
        row = {'state_configuration': state_config, 'state_term': state_term, 'state_J': state_J, 
               'energy': state_energy, 'energy_uncertainty': state_uncertainty,'is_from_theory': is_from_theory}
        portal_df.loc[len(portal_df.index)] = row

    portal_df.to_csv(filename, index=False)

    print(filename + ' has been written')

    return

def write_matrix_csv(element, filepath, mapping, gs_parity, theory_shift, expt_shift):
    '''
    This function writes the matrix element csv file
    '''
    filename = element + '_Matrix_Elements_Theory.csv'
    f = open(filepath + 'E1.RES', 'r') 
    lines = f.readlines()
    f.close()

    f = open(filepath + 'E1MBPT.RES', 'r') 
    lines_MBPT = f.readlines()
    f.close()

    # Obtain energies from CI+MBPT
    energies_MBPT = []
    for line in lines_MBPT[1:]:
        matrix_element_value = re.findall("\d+\.\d+", line)[:1]
        matrix_element_value = matrix_element_value[0] if matrix_element_value else None
        energies_MBPT.append(matrix_element_value)
    
    df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty'])
    
    tr_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty', 
                               'energy1(cm-1)', 'energy2(cm-1)',
                               'wavelength(nm)','transition_rate(s-1)'])
    
    i = 0
    for line in lines[1:]:
        # E1.RES format: < conf2 || E1 || conf1 > E1_L  E1_V  E2  E1  E2-E1  WL  Tr. Rate
        matrix_element = re.findall(r'\<.*?\>', line)[0]
        Tk = re.findall(r'\|\|.*?\|\|', matrix_element)[0].replace('||', '').replace(' ', '')
        state1 = re.findall(r'\|\|.*?\>', matrix_element)[0].replace(Tk, '').replace('||', '').replace('>', '')
        state2 = re.findall(r'\<.*?\|\|', matrix_element)[0].replace('<', '').replace('||', '')
        conf1 = '.'.join(state1.split()[:-1])
        conf2 = '.'.join(state2.split()[:-1])
        term1 = state1.split()[-1][0:2]
        term2 = state2.split()[-1][0:2]
        J1 = state1.split()[-1][-1]
        J2 = state2.split()[-1][-1]
        energy1 = re.findall("\d+\.\d+", line)[2:4][1]
        energy2 = re.findall("\d+\.\d+", line)[2:4][0]
        wavelength = float(re.findall("\d+\.\d+", line)[5])

        matrix_element_value = re.findall("\d+\.\d+", line)[:1]
        matrix_element_value = matrix_element_value[0] if matrix_element_value else None
        uncertainty = round(abs(float(matrix_element_value) - float(energies_MBPT[i])),5) if matrix_element_value else None
        i += 1

        # Use mapping to correct confs and terms and use experimental energies
        c1, c2 = False, False
        energy1cm, energy2cm = 0.0, 0.0
        for line_theory in mapping:
            if line_theory[1][4] == '-': continue
            if abs(float(line_theory[1][4]) - float(energy1)) < 1e-7:
                conf1 = line_theory[1][5]
                term1 = line_theory[1][1]
                J1 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-': 
                    energy1cm = float(line_theory[0][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(expt_shift)
                else:
                    energy1cm = float(line_theory[1][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(theory_shift)
                    
                c1 = True
            if abs(float(line_theory[1][4]) - float(energy2)) < 1e-7:
                conf2 = line_theory[1][5]
                term2 = line_theory[1][1]
                J2 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-': 
                    energy2cm = float(line_theory[0][3])
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(expt_shift)
                else:
                    energy2cm = float(line_theory[1][3])
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(theory_shift)
                c2 = True

        if c1 and c2:
            wavelength = 1e7/abs(energy2cm-energy1cm)
            if energy2 > energy1:
                trate = (2.02613*10**18)/((2*int(J1)+1)*(wavelength*10)**3)*float(matrix_element_value)**2
            else:
                trate = (2.02613*10**18)/((2*int(J2)+1)*(wavelength*10)**3)*float(matrix_element_value)**2
                
            row = {'state_one_configuration': conf1, 'state_one_term': term1, 'state_one_J': J1,
                   'state_two_configuration': conf2, 'state_two_term': term2, 'state_two_J': J2,
                   'matrix_element': matrix_element_value, 'matrix_element_uncertainty': uncertainty}
            df.loc[len(df.index)] = row
            
            trrow = {'state_one_configuration': conf1, 'state_one_term': term1, 'state_one_J': J1,
                   'state_two_configuration': conf2, 'state_two_term': term2, 'state_two_J': J2,
                   'matrix_element': matrix_element_value, 'matrix_element_uncertainty': uncertainty,
                   'energy1(cm-1)': f"{energy1cm:.2f}", 'energy2(cm-1)': f"{energy2cm:.2f}",
                   'wavelength(nm)': f"{wavelength:.2f}", 'transition_rate(s-1)': f"{trate:.4e}"}
            tr_df.loc[len(df.index)] = trrow
    
    tr_df.to_csv('tr_test.csv', index=False)
    df.to_csv(filename, index=False)
    print(filename + ' has been written')

    return

def nist_parity(term):
    ''' 
    this function determines parity of energy level based on term from NIST designated with a '*' 
    '''
    if '*' in term:
        parity = 'odd'
    else:
        parity = 'even'
    
    return parity

def find_energy_shift(df):
    '''
    this function finds the energy shift between lowest odd and even parity energy levels
    '''
    ground_parity = nist_parity(df['state_term'].values[:1][0])

    for index, row in df.iterrows():
        parity = nist_parity(row['state_term'])
        if parity != ground_parity:
            energy_shift = row['energy']
            break
    
    return energy_shift

if __name__ == "__main__":
    atom = 'Sr I'
    ri = False
    fac = 2 # maximum energy difference (in percent) for comparison
    name = atom.replace(" ","_")

    '''
    Method:
    1. Read all CONF*.RES files from the /DATA_RAW directory and add uncertainties to CONF.RES file in csv format
        - CONFFINALeven.RES
        - CONFFINALodd.RES
        - CONFFINALevenMBPT.RES
        - CONFFINALoddMBPT.RES
        - E1.RES
        - E1MBPT.RES
    2. Use Vipul's code to correct misidentified configurations
    3. Reformat data for use on Atom portal
    '''
    
    # 1. Write new CONF.RES (CONFFINAL.csv) with uncertainties
    # TODO - modularize this part to export and import output files
    raw_path = "DATA_RAW/"
    if os.path.isdir(raw_path):
        print('Reading raw files from ' + raw_path)
    else:
        os.makedirs(os.path.dirname(raw_path), exist_ok=True)
        print('Please put raw files in ' + raw_path)
        print('The files should be named: CONFFINALeven.RES, CONFFINALodd.RES, CONFFINALevenMBPT.RES, CONFFINALoddMBPT.RES, E1.RES, E1MBPT.RES')
        sys.exit()
    confs, terms, energies_cm, uncertainties, theory_shift, gs_parity, matrix_file_exists = write_new_conf_res(name, raw_path)

    # 2. Use Vipul's code to correct misidentified configurations
    # Store filtered data of even or odd parity in DATA_Filtered/NIST/ 
    url_nist = generate_asd_url(atom)
    data_nist = generate_df_from_asd(url_nist)
    data_nist = reformat_df_to_atomdb(data_nist)
    NIST_shift = find_energy_shift(data_nist)

    # Create directories for filtered data if it doesn't exist
    path_filtered_nist = "DATA_Filtered/NIST/"
    path_filtered_theory = "DATA_Filtered/UD/"
    os.makedirs(os.path.dirname(path_filtered_nist), exist_ok=True)
    os.makedirs(os.path.dirname(path_filtered_theory), exist_ok=True)

    df_to_csv(data_nist,"DATA_Filtered/NIST/"+atom,'odd')
    df_to_csv(data_nist,"DATA_Filtered/NIST/"+atom,'even')
    df_to_csv(data_nist,"DATA_Filtered/NIST/"+atom)
    
    # Filter NIST and UD data
    path_nist_even = "DATA_Filtered/NIST/"+name+"_NIST_Even.csv"
    path_ud_even = "DATA_Filtered/UD/"+name+"_UD_Even.csv"

    path_nist_odd = "DATA_Filtered/NIST/"+name+"_NIST_Odd.csv"
    path_ud_odd = "DATA_Filtered/UD/"+name+"_UD_Odd.csv"
    
    # TODO - automatically set the maximum number of levels
    # Set maximum number of levels to be read from NIST for each parity
    nist_max_odd = 40
    nist_max_even = 30
    
    # Export filtered data to output directory
    path_output = "DATA_Output/"
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
    
    # Filtering
    data_final_even = MainCode(path_nist_even, path_ud_even, nist_max_even)
    data_final_odd = MainCode(path_nist_odd, path_ud_odd, nist_max_odd)
    
    path = "DATA_Output/"+name+"_Even.txt" 
    ConvertToTXT(data_final_even, path)
    path = "DATA_Output/"+name+"_Odd.txt" 
    ConvertToTXT(data_final_odd, path)
    
    ## Finding Missing Levels
    data_final_even_missing = Missing_Levels(data_final_even)
    data_final_odd_missing = Missing_Levels(data_final_odd)    

    path = "DATA_Output/"+name+"_Even+missing.txt" 
    ConvertToTXT(data_final_even_missing, path)
    path = "DATA_Output/"+name+"_Odd+missing.txt" 
    ConvertToTXT(data_final_odd_missing, path)
    
    # 3. Create mapping of NIST data to theory data and reformat data for use on Atom portal
    mapping = create_mapping()

    write_energy_csv(name, mapping, NIST_shift, theory_shift, gs_parity)
    if matrix_file_exists: 
        print('Writing matrix elements...')
        write_matrix_csv(name, raw_path, mapping, gs_parity, theory_shift, NIST_shift)
    else:
        print('E1.RES files were not found, so matrix csv file was not generated')