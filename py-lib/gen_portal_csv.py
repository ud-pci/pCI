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
from glob import glob
from utils import get_dict_value

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

def create_mapping(num_levels_even, num_levels_odd):
    '''
    This function reads Vipul's energy level table and creates the mapping between experimental and theory data
    Data to map: Config, Term, J, Energy(cm-1)
    '''
    filepath = "DATA_Output/"+name+"_Even.txt"
    f = open(filepath, 'r')
    lines = f.readlines()[:num_levels_even + 1]
    f.close()

    filepath = "DATA_Output/"+name+"_Odd.txt"
    f = open(filepath, 'r')
    lines = lines + f.readlines()[:num_levels_odd + 1]
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

def normalize_config(config):
    '''
    Normalize an electronic configuration string by sorting the subshells.
    Example: "5s.4d" -> "4d.5s
    '''
    if not isinstance(config, str):
        return config
    
    parts = [part.strip() for part in config.split('.') if part.strip()]
    parts.sort(key=lambda x: (
        int(re.match(r'(\d+)', x).group(1)) if re.match(r'(\d+)', x) else 0,
        x
    ))
    
    return ".".join(parts)

def write_new_conf_res(name, filepath, data_nist):
    '''
    This function creates a new CONF.RES file with theory uncertainties
    '''
    
    second_order_exists = True
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
        second_order_exists = False
    if not os.path.exists(filepath + 'CONFFINALoddMBPT.RES'):
        print('CONFFINALoddMBPT.RES not found in', filepath)
        second_order_exists = False
    if not os.path.exists(filepath + 'E1.RES'):
        print('E1.RES not found in', filepath)
        matrix_file_exists = False
    if not os.path.exists(filepath + 'E1MBPT.RES'):
        print('E1MBPT.RES not found in', filepath)
        second_order_exists = False

    # Read CONF.RES files
    conf_res_odd, full_res_odd, swaps_odd, fixes_odd = cmp_res(filepath + 'CONFFINALodd.RES', filepath + 'CONFFINALoddMBPT.RES')
    conf_res_even, full_res_even, swaps_even, fixes_even = cmp_res(filepath + 'CONFFINALeven.RES', filepath + 'CONFFINALevenMBPT.RES')
    swaps = swaps_odd + swaps_even
    fixes = fixes_odd + fixes_even
    
    # Merge even and odd parity CONF.RES files and obtain uncertainties
    gs_parity, merged_res = merge_res(conf_res_even, conf_res_odd, second_order_exists)

    confs, terms, energies_au, energies_cm, energies_au_MBPT, energies_cm_MBPT, uncertainties = [], [], [], [], [], [], []
    with open('final_res.csv','w') as f:
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
    
    even_res = [conf for conf in merged_res if find_parity(conf[0]) == 'even']
    odd_res = [conf for conf in merged_res if find_parity(conf[0]) == 'odd']
    
    # Update uncertainties after merging even and odd parity CONF.RES files
    if second_order_exists:
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
        
        if normalize_config(conf) == normalize_config(nist_conf):
                gs_exists = True
                print('ground state found:', conf)
        else:
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

def write_matrix_csv(element, filepath, mapping, gs_parity, theory_shift, expt_shift, swaps, fixes, ignore_g, min_unc_per):
    '''
    This function writes the matrix element csv file
    '''
    matrix_element_filename = element + '_Matrix_Elements_Theory.csv'
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
        
    # List of default minimum uncertainties for different systems
    default_min_uncertainties = {
        'Mg1': 0.3,
        'Ca1': 1.3,
        'Sr1': 1.5
    }
    
    # Use default uncertainty if defined
    if element in default_min_uncertainties:
        min_unc_per = default_min_uncertainties[element]
        print('DEFAULT MINIMUM MATRIX ELEMENT UNCERTAINTY USED FOR', element + ':', str(min_unc_per) + '%')

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

        # Set minimum matrix element value
        if matrix_element_value == '0.00000':
            matrix_element_value = '0.00001'
            
        # Set minimum uncertainty
        try:
            uncertainty = '{:,.5f}'.format(math.sqrt(float(uncertainty)**2 + (min_unc_per/100)**2))
        except TypeError:
            uncertainty = '-'
            print('Type Error for: ' + conf1 + ' ' + term1+J1 + ', ' + conf2 + ' ' + term2+J2)
            continue
        except ValueError as ve:
            print('Uncertainty not found for: ' + conf1 + ' ' + term1+J1 + ', ' + conf2 + ' ' + term2+J2)
            continue
        
        # Use mapping to correct confs and terms and use experimental energies
        c1, c2 = False, False
        energy1cm, energy2cm = 0.0, 0.0
        for line_theory in mapping:
            # mapping structure:
            # [expt=[conf, term, J, energy, unc], theory=[conf, term, J, energy, unc, final_conf, energy_au]]
            if line_theory[1][6] == '-': continue
            if line_theory[1][6] == energy1:
                conf1 = line_theory[1][5]
                term1 = line_theory[1][1]
                if ignore_g:
                    if 'g' in conf1 or 'G' in term1:
                        continue
                J1 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy, configuration and term
                if line_theory[0][3] != '-': 
                    conf1 = line_theory[0][0]
                    term1 = line_theory[0][1]
                    J1 = line_theory[0][2]
                    energy1cm = float(line_theory[0][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(expt_shift)
                else:
                    energy1cm = float(line_theory[1][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(theory_shift)
                c1 = True
                
            if line_theory[1][6] == energy2:
                conf2 = line_theory[1][5]
                term2 = line_theory[1][1]
                if ignore_g:
                    if 'g' in conf2 or 'G' in term2:
                        continue
                J2 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-': 
                    conf2 = line_theory[0][0]
                    term2 = line_theory[0][1]
                    J2 = line_theory[0][2]
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
    
    num_E1 = len(df)
    print('TOTAL MATRIX ELEMENTS:', num_E1)
    
    df.to_csv(matrix_element_filename, index=False)
    print(matrix_element_filename + ' has been written')
    
    tr_df.to_csv(transition_rate_filename, index=False)
    print(transition_rate_filename + ' has been written')

    return num_E1

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

def convert_roman_to_num(roman):
    roman_val = {'I':1, 'V':5, 'X':10, 'L':50}
    num = 0
    i = 0
    while (i < len(roman)):
        c1 = roman_val[roman[i]]
        if (i + 1 < len(roman)):
            c2 = roman_val[roman[i + 1]]
            if (c1 >= c2):
                num = num + c1
                i = i + 1
            else:
                num = num + c2 - c1
                i = i + 2
        else:
            num = num + c1
            i = i + 1
    
    return num

def find_ci_dirs(ci_path):
    os.chdir(ci_path)
    
    even_dirs = glob('even*')
    odd_dirs = glob('odd*')
    dtm_dirs = glob('tm*')
    
    if even_dirs and odd_dirs:
        use_path = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input(ci_path + ' directory was found - use data from this directory? ')))))
    
    if not use_path:
        return None, None, None, None
    
    if len(even_dirs) > 1:
        even_dir = input(str(len(even_dirs)) + ' even directories were found: (' + ', '.join(even_dirs) + ') - which one would you like to use data from? ')
    else:
        even_dir = even_dirs[0]
    if len(odd_dirs) > 1:
        odd_dir = input(str(len(odd_dirs)) + ' odd directories were found: (' + ', '.join(odd_dirs) + ') - which one would you like to use data from? ')
    else:
        odd_dir = odd_dirs[0]
    if len(dtm_dirs) > 1:
        dtm_dir1 = input(str(len(dtm_dirs)) + ' dtm directories were found: (' + ', '.join(dtm_dirs) + ') - which one would you like to use data from? ')
        dtm_dir2 = input('Select another dtm directory if desired: ')
        if dtm_dir2 not in dtm_dirs:
            if dtm_dir2 != '':
                print(dtm_dir2 + ' was not found.')
            dtm_dir2 = None
    elif len(dtm_dirs) == 1:
        dtm_dir1 = dtm_dirs[0]
        dtm_dir2 = None
    else:
        dtm_dir1 = None
        dtm_dir2 = None
        
    os.chdir('..')
        
    return even_dir, odd_dir, dtm_dir1, dtm_dir2

def combine_tm(raw_path):
    e1_res = []
    
    with open(raw_path + '/E1a.RES', 'r') as f:
        lines = f.readlines()
    with open(raw_path + '/E1.RES', 'w') as f:
        for line in lines:
            f.write(line)

    for line in lines[1:]:
        matrix_element = re.findall(r'\<.*?\>', line)[0]
        e1_res.append(matrix_element)

    with open(raw_path + '/E1b.RES', 'r') as f:
        lines2 = f.readlines()
    
    with open(raw_path + '/E1.RES', 'a') as f:
        for line in lines2[1:]:
            matrix_element = re.findall(r'\<.*?\>', line)[0]
            if (matrix_element) not in e1_res:
                f.write(line)

    with open(raw_path + '/E1MBPTa.RES', 'r') as f:
        lines = f.readlines()
    with open(raw_path + '/E1MBPT.RES', 'w') as f:
        for line in lines:
            f.write(line)

    for line in lines[1:]:
        matrix_element = re.findall(r'\<.*?\>', line)[0]
        e1_res.append(matrix_element)

    with open(raw_path + '/E1MBPTb.RES', 'r') as f:
        lines2 = f.readlines()
    
    with open(raw_path + '/E1MBPT.RES', 'a') as f:
        for line in lines2[1:]:
            matrix_element = re.findall(r'\<.*?\>', line)[0]
            if (matrix_element) not in e1_res:
                f.write(line)

if __name__ == "__main__":
    use_config_yml = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input('Using a config.yml file? ')))))
    
    if use_config_yml:
        config_yml = input("Input yml-file: ")
        config = read_yaml(config_yml)
        atom_name = get_dict_value(config['atom'],'name')
        
        # Parse atom name
        if len(atom_name.split()) == 1:
            atom = atom_name + ' I'
        elif len(atom_name.split()) == 2:
            if '+' in atom_name.split()[1]: 
                atom = atom_name
            elif atom_name.split()[1].isnumeric():
                atom = atom_name.split()[0] + ' ' + str(int(atom_name.split()[1])-1) + '+'
            else:
                atom = atom_name.split()[0] + ' ' + str(convert_roman_to_num(atom_name.split()[1])-1) + '+'
        else:
            print('ERROR: atom name not supported')
            sys.exit()
        
        # conf parameters
        conf = get_dict_value(config, 'conf')
        even = get_dict_value(conf, 'even')
        even_J = get_dict_value(even, 'J')
        odd = get_dict_value(conf, 'odd')
        odd_J = get_dict_value(odd, 'J')
        even_dir = 'even' + str(even_J)[0] if even_J else None
        odd_dir = 'odd' + str(odd_J)[0] if odd_J else None
        tm_dir = 'tm' if even_J and odd_J else None
        tm_dir1 = None
        tm_dir2 = None
        
        # portal parameters
        portal = get_dict_value(config, 'portal')
        
        # set default ignore configurations with 'g' and terms with 'G'
        ignore_g = get_dict_value(portal, 'ignore_g') if portal else True
        
        # set default minimum uncertainty as percentage of value to 1.5
        min_uncertainty = float(get_dict_value(portal, 'min_uncertainty')) if portal else 1.5
    else:
        atom = input('Input name of atom: ')
        even_dir = None
        odd_dir = None
        tm_dir = None
        tm_dir1 = None
        tm_dir2 = None
        ignore_g = True
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
        if not even_dir or not odd_dir or not tm_dir:
            even_dir, odd_dir, tm_dir1, tm_dir2 = find_ci_dirs(dir_path + '/' + all_order_path)
        if even_dir and odd_dir:
            if os.path.isdir(all_order_path + '/' + even_dir):
                run('cp ci+all-order/' + even_dir + '/FINAL.RES DATA_RAW/CONFFINALeven.RES', shell=True)
            if os.path.isdir(all_order_path + '/' + odd_dir):   
                run('cp ci+all-order/' + odd_dir + '/FINAL.RES DATA_RAW/CONFFINALodd.RES', shell=True)
        if tm_dir1 and os.path.isdir(all_order_path + '/' + tm_dir1):   
            run('cp ci+all-order/' + tm_dir1 + '/E1.RES DATA_RAW/E1a.RES', shell=True)
            run('cp ci+all-order/' + tm_dir1 + '/E1.RES DATA_RAW/E1.RES', shell=True)
        if tm_dir2 and os.path.isdir(all_order_path + '/' + tm_dir2):   
            run('cp ci+all-order/' + tm_dir2 + '/E1.RES DATA_RAW/E1b.RES', shell=True)
        if even_dir and odd_dir or tm_dir1 or tm_dir2:
            print('data from ' + all_order_path + ' moved to DATA_RAW directory')
        os.chdir(dir_path)
    
    second_order_path = 'ci+second-order'
    if os.path.isdir(second_order_path):
        if not even_dir or not odd_dir or not tm_dir:
            even_dir, odd_dir, tm_dir1, tm_dir2 = find_ci_dirs(dir_path + '/' + second_order_path)
        if even_dir and odd_dir:
            if os.path.isdir(second_order_path + '/' + even_dir):
                run('cp ci+second-order/' + even_dir + '/FINAL.RES DATA_RAW/CONFFINALevenMBPT.RES', shell=True)
            if os.path.isdir(second_order_path + '/' + odd_dir):   
                run('cp ci+second-order/' + odd_dir + '/FINAL.RES DATA_RAW/CONFFINALoddMBPT.RES', shell=True)
        if tm_dir1 and os.path.isdir(second_order_path + '/' + tm_dir1):   
            run('cp ci+second-order/' + tm_dir1 + '/E1.RES DATA_RAW/E1MBPTa.RES', shell=True)
            run('cp ci+second-order/' + tm_dir1 + '/E1.RES DATA_RAW/E1MBPT.RES', shell=True)
        if tm_dir2 and os.path.isdir(second_order_path + '/' + tm_dir2):   
            run('cp ci+second-order/' + tm_dir2 + '/E1.RES DATA_RAW/E1MBPTb.RES', shell=True)
        if even_dir and odd_dir or tm_dir1 or tm_dir2:
            print('data from ' + second_order_path + ' moved to DATA_RAW directory')
        os.chdir(dir_path)
    
    if tm_dir1 and tm_dir2:
        combine_tm(data_raw_path)
        
    # Parse NIST Atomic Spectral Database for full list of energy levels
    url_nist = generate_asd_url(atom)
    print(url_nist)
    data_nist = generate_df_from_asd(url_nist)
    
    # Write new CONF.RES (CONFFINAL.csv) with uncertainties
    raw_path = dir_path + "/DATA_RAW/"
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
    
    # Set maximum number of levels to be read from NIST equal to number of levels in CONFFINAL.RES
    with open(raw_path + 'CONFFINALeven.RES','r') as f:
        lines = f.readlines()
        num_levels_theory_even = len(lines) - 1
    nist_max_even = num_levels_theory_even
    
    with open(raw_path + 'CONFFINALodd.RES','r') as f:
        lines = f.readlines()
        num_levels_theory_odd = len(lines) - 1
    nist_max_odd = num_levels_theory_odd
    
    # Set maximum number of levels to be outputted in csv files
    num_levels_output_even = 1 if uncertainties[0] == '-' else 0
    for i in range(1, num_levels_theory_even):
        if uncertainties[i] == '-':
            break
        else:
            num_levels_output_even += 1
    num_levels_output_odd = num_levels_theory_odd
    #num_levels_output_odd = 0
    #for i in range(num_levels_theory_even, num_levels_theory_even + num_levels_theory_odd):
    #    num_levels_output_odd = i + 1 - num_levels_theory_even
    #    if uncertainties[i] == '-':
    #        break
        
    print('Number of even parity levels: ', num_levels_output_even)
    print('Number of odd parity levels: ', num_levels_output_odd)

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
    mapping = create_mapping(num_levels_output_even, num_levels_output_odd)

    # Reformat atom name to numerical value for filenames
    if name.split('_')[1].isalpha():
        name = name.split('_')[0] + str(convert_roman_to_num(name.split('_')[1]))
    else:
        name = name.split('_')[0] + name.split('_')[1]
        
    write_energy_csv(name, mapping, NIST_shift, theory_shift, gs_parity)
    
    # Create a list of all possible transitions between states
    print('even parity configurations:')
    even_confs = []
    for line in mapping[:num_levels_output_even]:
        if line[0][0] == '-': 
            continue
        even_confs.append([line[0][0],line[0][1],line[0][2]])
    print('odd parity configurations:')
    odd_confs = []
    for line in mapping[num_levels_output_even:]:
        if line[0][0] == '-': 
            continue
        odd_confs.append([line[0][0],line[0][1],line[0][2]])

    possible_E1 = []
    for conf_odd in odd_confs:
        J_odd = int(conf_odd[2])
        for conf_even in even_confs:
            J_even = int(conf_even[2])
            if J_even == 0 and J_odd == 0: continue
            if abs(J_even - J_odd) <= 1:
                possible_E1.append([conf_odd, conf_even])
    num_possible_E1 = len(possible_E1)
    print('Number of possible E1: ', len(possible_E1))
    
    if matrix_file_exists: 
        print('Writing matrix elements...')
        num_E1 = write_matrix_csv(name, raw_path, mapping, gs_parity, theory_shift, NIST_shift, swaps, fixes, ignore_g, min_uncertainty)
        print(str(round(num_E1/num_possible_E1*100,2)) + "% of possible E1 accounted for")
    else:
        print('E1.RES files were not found, so matrix csv file was not generated')