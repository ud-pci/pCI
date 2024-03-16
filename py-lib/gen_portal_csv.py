import yaml
import re
import sys
import os
import pandas as pd
from fractions import Fraction
from UDRead import *
from parse_asd import *
from get_atomic_term import *
from pathlib import Path
from subprocess import run
from compare_res import *

def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

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

def write_new_conf_res(name, filepath, data_nist):
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

    # Read CONF.RES files
    conf_res_odd, full_res_odd, swaps_odd, fixes_odd = cmp_res(filepath + 'CONFFINALodd.RES', filepath + 'CONFFINALoddMBPT.RES')
    conf_res_even, full_res_even, swaps_even, fixes_even = cmp_res(filepath + 'CONFFINALeven.RES', filepath + 'CONFFINALevenMBPT.RES')
    swaps = swaps_odd + swaps_even
    fixes = fixes_odd + fixes_even
    
    # Merge even and odd parity CONF.RES files and obtain uncertainties
    gs_parity, merged_res = merge_res(conf_res_even, conf_res_odd)

    confs, terms, energies_au, energies_cm, energies_au_MBPT, energies_cm_MBPT, uncertainties = [], [], [], [], [], [], []
    f = open('final_res.csv','w')
    f.write('conf,term,J,energy(a.u.),energy(cm-1),energy_MBPT(a.u.),energy_MBPT(cm-1),uncertainty\n')
    for conf in merged_res:
        confs.append(conf[0])
        terms.append(conf[1])
        energies_au.append(conf[2])
        energies_cm.append(conf[3])
        energies_au_MBPT.append(conf[4])
        energies_cm_MBPT.append(conf[5])
        uncertainties.append(conf[6])
        f.write(','.join(str(i) for i in conf) + '\n')
    f.close()
    
    even_res = [conf for conf in merged_res if find_parity(conf[0]) == 'even']
    odd_res = [conf for conf in merged_res if find_parity(conf[0]) == 'odd']
    
    # Update uncertainties after merging even and odd parity CONF.RES files
    for ilvl in range(len(even_res)):
        full_res_even[ilvl][13] = even_res[ilvl][6]
    for ilvl in range(len(odd_res)):
        full_res_odd[ilvl][13] = odd_res[ilvl][6]

    # Determine if ground state level exists in theory results
    gs_exists = False
    nist_conf = data_nist['Configuration'].iloc[0]
    nist_term = data_nist['Term'].iloc[0]
    nist_J = data_nist['J'].iloc[0]

    for i in range(len(confs)): 
        conf = confs[i]
        term = terms[i].split(',')[0]
        J = terms[i].split(',')[1]
        
        str_diff, num_diff = SubtractStr(nist_conf, conf)
        if num_diff > 0:
            if nist_conf.replace(str_diff, '') == conf and nist_term == term and nist_J == J:
                gs_exists = True
                print('ground state found:', confs[i])
    
    # If ground state level not found in theory results, ask user for energy (a.u.) of ground state level
    if not gs_exists:
        print(data_nist['Configuration'].iloc[0], ' not in', confs)
        th_gs_au = float(input('Ground state level was not found in theory results. Enter energy (a.u.) of ground state level: '))
        gs_parity = find_parity(data_nist['Configuration'].iloc[0])

    # TODO - reimplement shift by user-inputted ground state
    
    # Determine energy shift between odd and even parity lowest energy levels
    ht_to_cm = 219474.63 # hartree to cm-1
    min_energy_even = even_res[0][2]
    min_energy_odd = odd_res[0][2]
    energy_shift = abs(ht_to_cm * (min_energy_even - min_energy_odd))
        
    # List of J values in theory results
    theory_J = {}
    theory_J['even'] = [conf[1].split(',')[1] for conf in even_res]
    theory_J['odd'] = [conf[1].split(',')[1] for conf in odd_res]
    
    # Write csv-formatted CONF.RES files with uncertainties
    convert_res_to_csv(filepath + 'CONFFINALeven.RES', full_res_even, gs_exists, name)
    convert_res_to_csv(filepath + 'CONFFINALodd.RES', full_res_odd, gs_exists, name)

    return confs, terms, energies_au, energies_cm, uncertainties, energy_shift, theory_J, gs_parity, matrix_file_exists, gs_exists, swaps, fixes

def convert_type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    s = s.replace(" ", "")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        return s 
    
def convert_res_to_csv(filename, full_res, gs_exists, name):

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

    for row in full_res:
        f.write(','.join([str(item) for item in row[0:2]]) + ',' + row[2].replace(',','') + ',' + ','.join([str(item) for item in row[3:]]) + '\n')

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

def write_matrix_csv(element, filepath, mapping, gs_parity, theory_shift, expt_shift, swaps, fixes, ignore_g):
    '''
    This function writes the matrix element csv file
    '''
    matrix_element_filename = element + '_Matrix_Elements.csv'
    transition_rate_filename = element + '_Transition_Rates.csv'

    # Read E1.RES and E1MBPT.RES and return E1.RES table with uncertainties
    e1_res = cmp_matrix_res(filepath + 'E1.RES', filepath + 'E1MBPT.RES', swaps, fixes)
    
    df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty'])
    
    tr_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty', 
                               'energy1(cm-1)', 'energy2(cm-1)',
                               'wavelength(nm)','transition_rate(s-1)'])

    if ignore_g:
        print('IGNORING G STATES')
    
    for line in e1_res:
        # E1.RES format: [conf11, term11, conf12, term12, me1, uncertainty, energy1, energy2, wavelength]
        conf1 = line[0]
        conf2 = line[2]
        term1 = line[1][0:2]
        term2 = line[3][0:2]
        J1 = line[1][2]
        J2 = line[3][2]
        matrix_element_value = line[4]
        uncertainty = line[5]
        energy1 = line[6]
        energy2 = line[7]
        wavelength = line[8]

        # Use mapping to correct confs and terms and use experimental energies
        c1, c2 = False, False
        energy1cm, energy2cm = 0.0, 0.0
        for line_theory in mapping:
            # mapping structure:
            # [expt=[conf, term, J, energy, unc], theory=[conf, term, J, energy, unc, final_conf, energy_au]]
            if line_theory[1][6] == '-': continue
            if abs(float(line_theory[1][6]) - float(energy1)) < 1e-7:
                conf1 = line_theory[1][5]
                term1 = line_theory[1][1]
                if ignore_g:
                    if 'g' in conf1 or 'G' in term1:
                        continue
                J1 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-': 
                    conf1 = line_theory[0][0]
                    energy1cm = float(line_theory[0][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(expt_shift)
                else:
                    energy1cm = float(line_theory[1][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(theory_shift)
                c1 = True
                
            if abs(float(line_theory[1][6]) - float(energy2)) < 1e-7:
                conf2 = line_theory[1][5]
                term2 = line_theory[1][1]
                if ignore_g:
                    if 'g' in conf2 or 'G' in term2:
                        continue
                J2 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-': 
                    conf2 = line_theory[0][0]
                    energy2cm = float(line_theory[0][3])
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(expt_shift)
                else:
                    energy2cm = float(line_theory[1][3])
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(theory_shift)
                c2 = True

        if c1 and c2:
            wavelength = 1e7/(energy2cm-energy1cm)
            if energy2 > energy1:
                trate = (2.02613*10**18)/((2*int(J1)+1)*(abs(wavelength)*10)**3)*float(matrix_element_value)**2
            else:
                trate = (2.02613*10**18)/((2*int(J2)+1)*(abs(wavelength)*10)**3)*float(matrix_element_value)**2
                
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
    
    print('TOTAL MATRIX ELEMENTS:',len(df))
    
    df.to_csv(matrix_element_filename, index=False)
    print(matrix_element_filename + ' has been written')
    
    tr_df.to_csv(transition_rate_filename, index=False)
    print(transition_rate_filename + ' has been written')

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
    # Read atom name from config.yml if it exists
    config_exists = os.path.isfile('config.yml')

    atom = ''
    if config_exists:
        config = read_yaml('config.yml')
        ignore_g = config['portal']['ignore_g']
        config_name = config['system']['name']
        if len(config_name.split()) == 1:
            atom = config_name + ' I'
    else:
        atom = input('Input name of atom: ')
    name = atom.replace(" ","_")
    
    ri = False # 
    fac = 2 # maximum energy difference (in percent) for comparison
    
    # Find input files from directories if they exist and put into DATA_RAW directory
    dir_path = os.getcwd()
    data_raw_path = 'DATA_RAW'
    if not os.path.isdir(data_raw_path):
        Path(data_raw_path).mkdir(parents=True, exist_ok=True)
    
    all_order_path = 'ci+all-order'
    if os.path.isdir(all_order_path):
        use_path = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input(all_order_path + ' directory was found - use data from this directory? '))))
        if use_path:
            if os.path.isdir(all_order_path + '/even'):
                run('cp ci+all-order/even/CONFFINAL.RES DATA_RAW/CONFFINALeven.RES', shell=True)
            if os.path.isdir(all_order_path + '/odd'):   
                run('cp ci+all-order/odd/CONFFINAL.RES DATA_RAW/CONFFINALodd.RES', shell=True)
            if os.path.isdir(all_order_path + '/dtm'):   
                run('cp ci+all-order/dtm/E1.RES DATA_RAW/E1.RES', shell=True)
            print('data from ' + all_order_path + ' moved to DATA_RAW directory')
    
    second_order_path = 'ci+second-order'
    if os.path.isdir(second_order_path):
        use_path = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input(second_order_path + ' directory was found - use data from this directory? '))))
        if use_path:
            if os.path.isdir(second_order_path + '/even'):
                run('cp ci+second-order/even/CONFFINAL.RES DATA_RAW/CONFFINALevenMBPT.RES', shell=True)
            if os.path.isdir(second_order_path + '/odd'):   
                run('cp ci+second-order/odd/CONFFINAL.RES DATA_RAW/CONFFINALoddMBPT.RES', shell=True)
            if os.path.isdir(second_order_path + '/dtm'):   
                run('cp ci+second-order/dtm/E1.RES DATA_RAW/E1MBPT.RES', shell=True)
            print('data from ' + second_order_path + ' moved to DATA_RAW directory')
        
    # Parse NIST Atomic Spectral Database for full list of energy levels
    url_nist = generate_asd_url(atom)
    data_nist = generate_df_from_asd(url_nist)
    
    # Write new CONF.RES (CONFFINAL.csv) with uncertainties
    raw_path = "DATA_RAW/"
    if os.path.isdir(raw_path):
        print('Reading raw files from ' + raw_path)
    else:
        os.makedirs(os.path.dirname(raw_path), exist_ok=True)
        print('Please put raw files in ' + raw_path)
        print('The files should be named: CONFFINALeven.RES, CONFFINALodd.RES, CONFFINALevenMBPT.RES, CONFFINALoddMBPT.RES, E1.RES, E1MBPT.RES')
        sys.exit()
    confs, terms, energies_au, energies_cm, uncertainties, theory_shift, theory_J, gs_parity, matrix_file_exists, gs_exists, swaps, fixes = write_new_conf_res(name, raw_path, data_nist)

    data_nist = reformat_df_to_atomdb(data_nist, theory_J)
    if gs_exists:
        NIST_shift = find_energy_shift(data_nist)
    else:
        NIST_shift = 0
        theory_shift = 0

    # Store filtered data of even or odd parity in DATA_Filtered/NIST/ 
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
    
    # TODO - automatically set the maximum number of levels based on last matching configuration between theory and NIST
    # Set maximum number of levels to be read from NIST for each parity
    nist_max_odd = 44
    nist_max_even = 40
    
    # Export filtered data to output directory
    path_output = "DATA_Output/"
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
    
    # Use Vipul's code to correct misidentified configurations
    data_final_even = MainCode(path_nist_even, path_ud_even, nist_max_even, gs_exists)
    data_final_odd = MainCode(path_nist_odd, path_ud_odd, nist_max_odd, gs_exists)
    
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
        write_matrix_csv(name, raw_path, mapping, gs_parity, theory_shift, NIST_shift, swaps, fixes, ignore_g)
    else:
        print('E1.RES files were not found, so matrix csv file was not generated')