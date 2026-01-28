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
        try:
            # Extract energy difference percentage (last column, remove '%' sign)
            energy_diff_pct = float(line.split()[14].rstrip('%'))
        except:
            # If parsing fails, set to 0 (perfect match)
            energy_diff_pct = 0.0

        # select relevant data for portal database
        if NIST_config != 'Config':
            mapping.append([[NIST_config, NIST_term, NIST_J, NIST_energy, NIST_uncertainty],
                    [theory_config, theory_term, theory_J, theory_energy_cm, theory_uncertainty, corrected_config, theory_energy_au],
                    energy_diff_pct])
    
    return mapping

def generate_mapping_fixes(mapping):
    """
    Generate fixes for E1.RES based on differences between theory_config and corrected_config.

    When NIST matching determines a different label (corrected_config) than the original theory label (theory_config), we need to fix E1.RES to use the corrected label.

    Args:
        mapping: list of [[NIST_data], [theory_data], energy_diff_pct] where theory_data = [theory_config, theory_term, theory_J, energy_cm, uncertainty, corrected_config, energy_au]

    Returns:
        list of fixes: [[old_conf, old_term, energy, new_conf, new_term], ...]
    """
    fixes = []
    seen = set()

    for entry in mapping:
        theory_data = entry[1]

        if theory_data[0] == '-':
            continue

        theory_config = theory_data[0]
        theory_term = theory_data[1]
        theory_J = theory_data[2]
        corrected_config = theory_data[5]
        theory_energy_au = theory_data[6]

        # If configs differ, generate a fix
        if theory_config != corrected_config and theory_config != '-' and corrected_config != '-':
            # Format term with comma for fix format
            term_with_comma = theory_term + ',' + str(theory_J)

            # Create fix key to avoid duplicates
            fix_key = (theory_config, term_with_comma, theory_energy_au)
            if fix_key not in seen:
                seen.add(fix_key)
                # Fix format: [old_conf, old_term, energy, new_conf, new_term]
                fixes.append([theory_config, term_with_comma, theory_energy_au, corrected_config, term_with_comma])

    if fixes:
        print(f'Generated {len(fixes)} mapping-based fixes:')
        for fix in fixes:
            print(f'  {fix[0]} {fix[1].replace(",","")} -> {fix[3]} {fix[4].replace(",","")}')

    return fixes

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
    # E1.RES and E1MBPT.RES are in DATA_Processed/, check there
    if not os.path.exists('DATA_Processed/E1.RES'):
        print('E1.RES not found in DATA_Processed/')
        matrix_file_exists = False
    if not os.path.exists('DATA_Processed/E1MBPT.RES'):
        print('E1MBPT.RES not found in DATA_Processed/')
        second_order_exists = False

    # Read CONF.RES files
    conf_res_odd, full_res_odd, swaps_odd, fixes_odd, unmatched_odd, allorder_e2l_odd, mbpt_e2l_odd = cmp_res(filepath + 'CONFFINALodd.RES', filepath + 'CONFFINALoddMBPT.RES')
    conf_res_even, full_res_even, swaps_even, fixes_even, unmatched_even, allorder_e2l_even, mbpt_e2l_even = cmp_res(filepath + 'CONFFINALeven.RES', filepath + 'CONFFINALevenMBPT.RES')
    swaps = swaps_odd + swaps_even
    # Combine energy to level mappings from both parities (all-order mapping for fixes, MBPT mapping for swaps)
    energy_to_level = {**allorder_e2l_odd, **allorder_e2l_even}
    mbpt_energy_to_level = {**mbpt_e2l_odd, **mbpt_e2l_even}
    # fixes_odd and fixes_even are tuples (fixes1, fixes2)
    # Combine all-order fixes together and MBPT fixes together
    fixes1_odd, fixes2_odd = fixes_odd
    fixes1_even, fixes2_even = fixes_even
    fixes = (fixes1_odd + fixes1_even, fixes2_odd + fixes2_even)
    unmatched_energies = {'odd': unmatched_odd, 'even': unmatched_even}
    
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
                break
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

    return confs, terms, energies_au, energies_cm, uncertainties, energy_shift, theory_J, gs_parity, matrix_file_exists, gs_exists, swaps, fixes, unmatched_energies, energy_to_level, mbpt_energy_to_level

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

def write_energy_csv(name, mapping, NIST_shift, theory_shift, gs_parity, min_energy_diff_percent):
    '''
    This function writes the energy csv file
    Mapping should already be filtered to only include levels within min_energy_diff_percent
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
            state_energy = level[1][3]
            state_uncertainty = level[1][4]
        else:
            is_from_theory = False
            state_config = level[0][0]
            state_term = level[0][1]
            state_J = level[0][2]
            state_energy = level[0][3]
            state_uncertainty = level[0][4]
        row = {'state_configuration': state_config, 'state_term': state_term, 'state_J': state_J,
               'energy': state_energy, 'energy_uncertainty': state_uncertainty,'is_from_theory': is_from_theory}
        portal_df.loc[len(portal_df.index)] = row

    # Sort by energy
    portal_df['energy'] = portal_df['energy'].astype(float)
    portal_df = portal_df.sort_values(by='energy').reset_index(drop=True)

    portal_df.to_csv(filename, index=False)

    print(f'{filename} has been written with {len(portal_df)} levels (min energy diff: {min_energy_diff_percent}%)')

    return

def write_matrix_csv(element, filepath, mapping, gs_parity, theory_shift, expt_shift, swaps, fixes, ignore_g, min_unc_per, min_energy_diff_percent, energy_to_level, mbpt_energy_to_level):
    '''
    This function writes the matrix element csv file
    Note: Transition rates are now calculated separately using generate_transition_rates.py
    '''
    matrix_element_filename = element + '_Matrix_Elements_Theory.csv'

    # Read E1.RES and E1MBPT.RES and return E1.RES table with uncertainties
    e1_res, unmatched_matrix = cmp_matrix_res('DATA_Processed/E1.RES', 'DATA_Processed/E1MBPT.RES', swaps, fixes, energy_to_level, mbpt_energy_to_level)

    df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'matrix_element', 'matrix_element_uncertainty'])
    
    # Track added transitions to avoid duplicates (E1.RES contains both A->B and B->A)
    added_transitions = set()

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

    print('MIN ENERGY DIFF PERCENT FOR NIST-THEORY:', str(min_energy_diff_percent) + '%')

    for line in e1_res:
        # E1.RES format: [conf11, term11, conf12, term12, me1, uncertainty, energy1, energy2, wavelength]
        conf1 = line[0]
        conf2 = line[2]
        term1 = line[1][0:2]
        term2 = line[3][0:2]
        J1_e1res = line[1][2]  # J from E1.RES (for matching)
        J2_e1res = line[3][2]
        J1 = J1_e1res  # Will be updated with matched state's J
        J2 = J2_e1res
        try:
            matrix_element_value = float(line[4])
            uncertainty = float(line[5])
        except (ValueError, TypeError):
            continue
        energy1 = line[6]
        energy2 = line[7]
        wavelength = line[8]
        
        # Set minimum uncertainty
        extra_uncertainty = matrix_element_value * min_unc_per / 100
        uncertainty = np.sqrt(uncertainty**2 + extra_uncertainty**2)
        if uncertainty == 0:
            uncertainty = 0.00001
        
        # Use mapping to correct confs and terms and use experimental energies
        c1, c2 = False, False
        energy1cm, energy2cm = 0.0, 0.0
        energy1_float = abs(float(energy1))  # Use absolute value (E1.RES uses negative binding energies)
        energy2_float = abs(float(energy2))
        energy_tolerance = 1e-6  # Tolerance for floating point comparison (in a.u.)

        # Track best matches (smallest energy difference)
        best_match1_diff = float('inf')
        best_match2_diff = float('inf')

        for line_theory in mapping:
            # mapping structure:
            # [expt=[conf, term, J, energy, unc], theory=[conf, term, J, energy, unc, final_conf, energy_au]]
            if line_theory[1][6] == '-': continue

            # Compare energies as floats with tolerance instead of exact string match
            theory_energy_au = abs(float(line_theory[1][6]))  # Use absolute value for comparison
            energy_diff1 = abs(theory_energy_au - energy1_float)

            # Check if this is a better match than previous best
            # Prefer exact energy matches, and use J as tiebreaker for equal energies
            is_better_match1 = False
            if energy_diff1 < energy_tolerance:
                if energy_diff1 < best_match1_diff:
                    is_better_match1 = True
                elif abs(energy_diff1 - best_match1_diff) < 1e-12:
                    # Same energy (within numerical precision) - prefer J match as tiebreaker
                    if line_theory[1][2] == J1_e1res:
                        is_better_match1 = True

            if is_better_match1:
                conf1 = line_theory[1][5]  # corrected_config for output
                term1 = line_theory[1][1]
                
                if ignore_g:
                    if 'g' in conf1 or 'G' in term1:
                        continue
                J1 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy, configuration and term
                if line_theory[0][3] != '-':
                    # Check energy difference between NIST and theory
                    nist_energy = float(line_theory[0][3])
                    theory_energy = float(line_theory[1][3])

                    # Calculate percentage difference
                    if nist_energy != 0:
                        energy_diff_percent = abs((nist_energy - theory_energy) / nist_energy * 100)
                    else:
                        energy_diff_percent = 0.0

                    # Skip if energy difference exceeds minimum threshold
                    if energy_diff_percent > min_energy_diff_percent:
                        continue

                    conf1 = line_theory[0][0]
                    term1 = line_theory[0][1]
                    J1 = line_theory[0][2]
                    energy1cm = nist_energy
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(expt_shift)
                else:
                    energy1cm = float(line_theory[1][3])
                    if find_parity(conf1) != gs_parity:
                        energy1cm = energy1cm + float(theory_shift)
                c1 = True
                best_match1_diff = energy_diff1

            energy_diff2 = abs(theory_energy_au - energy2_float)

            # Check if this is a better match than previous best for state 2
            is_better_match2 = False
            if energy_diff2 < energy_tolerance:
                if energy_diff2 < best_match2_diff:
                    is_better_match2 = True
                elif abs(energy_diff2 - best_match2_diff) < 1e-12:
                    # Same energy (within numerical precision) - prefer J match as tiebreaker
                    if line_theory[1][2] == J2_e1res:
                        is_better_match2 = True

            if is_better_match2:
                conf2 = line_theory[1][5]  # corrected_config for output
                term2 = line_theory[1][1]

                if ignore_g:
                    if 'g' in conf2 or 'G' in term2:
                        continue
                J2 = line_theory[1][2]
                # Check if NIST energy exists - if it does, overwrite theory energy
                if line_theory[0][3] != '-':
                    # Check energy difference between NIST and theory
                    nist_energy = float(line_theory[0][3])
                    theory_energy = float(line_theory[1][3])

                    # Calculate percentage difference
                    if nist_energy != 0:
                        energy_diff_percent = abs((nist_energy - theory_energy) / nist_energy * 100)
                    else:
                        energy_diff_percent = 0.0

                    # Skip if energy difference exceeds minimum threshold
                    if energy_diff_percent > min_energy_diff_percent:
                        continue

                    conf2 = line_theory[0][0]
                    term2 = line_theory[0][1]
                    J2 = line_theory[0][2]
                    energy2cm = nist_energy
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(expt_shift)
                else:
                    energy2cm = float(line_theory[1][3])
                    if find_parity(conf2) != gs_parity:
                        energy2cm = energy2cm + float(theory_shift)
                c2 = True
                best_match2_diff = energy_diff2

        if c1 and c2:
            # Apply E1 selection rules
            try:
                J1_float = float(J1) if '/' not in str(J1) else float(eval(J1))
                J2_float = float(J2) if '/' not in str(J2) else float(eval(J2))
                delta_J = abs(J1_float - J2_float)

                # Skip forbidden transitions
                if (J1_float == 0 and J2_float == 0) or delta_J > 1:
                    continue
            except:
                continue

            # Create unique transition identifier to avoid duplicates (both A->B and B->A)
            state1 = f"{conf1} {term1}{J1}"
            state2 = f"{conf2} {term2}{J2}"
            trans_id = tuple(sorted([state1, state2]))

            # Skip if we've already added this transition
            if trans_id in added_transitions:
                continue

            added_transitions.add(trans_id)

            row = {'state_one_configuration': conf1, 'state_one_term': term1, 'state_one_J': J1,
                   'state_two_configuration': conf2, 'state_two_term': term2, 'state_two_J': J2,
                   'matrix_element': matrix_element_value, 'matrix_element_uncertainty': uncertainty}
            df.loc[len(df.index)] = row
    
    num_E1 = len(df)
    print(f'TOTAL MATRIX ELEMENTS: {num_E1}')

    df.to_csv(matrix_element_filename, index=False)
    print(matrix_element_filename + ' has been written')

    return num_E1, unmatched_matrix

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

    energy_shift = 0.0
    for index, row in df.iterrows():
        parity = nist_parity(row['state_term'])
        if parity != ground_parity:
            # Skip if energy is missing (marked as '-')
            if row['energy'] != '-':
                energy_shift = float(row['energy'])
                break

    return energy_shift

def convert_num_to_roman(num):
    """
    Convert an integer to Roman numerals.

    Args:
        num: Integer to convert (1-3999)

    Returns:
        Roman numeral string
    """
    val = [
        100, 90, 50, 40,
        10, 9, 5, 4,
        1
    ]
    syms = [
        "C", "XC", "L", "XL",
        "X", "IX", "V", "IV",
        "I"
    ]
    roman_num = ''
    i = 0
    while num > 0:
        for _ in range(num // val[i]):
            roman_num += syms[i]
            num -= val[i]
        i += 1
    return roman_num

def convert_roman_to_num(roman):
    """
    Convert Roman numerals to an integer.

    Args:
        roman: Roman numeral string

    Returns:
        Integer value
    """
    roman_val = {
        'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100
    }
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

def atom_name_to_filename(atom):
    """
    Convert atom name to filename format.

    Examples:
    - Ba I → Ba1 (neutral)
    - Ba II or Ba+ → Ba2 (singly ionized)
    - Ba III or Ba++ → Ba3 (doubly ionized)

    Args:
        atom: Atom name in format "Element Roman" or "Element Charge"

    Returns:
        Filename format (e.g., "Ba1", "Ba2")
    """
    parts = atom.split()
    if len(parts) == 1:
        # Just element, assume neutral (I)
        return parts[0] + '1'

    element = parts[0]
    suffix = parts[1]

    if '+' in suffix:
        # Count number of + signs and add 1 (e.g., Ba+ → Ba2, Ba++ → Ba3)
        ionization = suffix.count('+') + 1
        return element + str(ionization)
    elif suffix.isnumeric():
        # Already in numeric format
        return atom.replace(' ', '')
    else:
        # Roman numeral (e.g., Ba I → Ba1, Ba II → Ba2)
        ionization = convert_roman_to_num(suffix)
        return element + str(ionization)

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

def combine_tm(raw_path, filtered_path):
    """
    Combine E1a/E1b and E1MBPTa/E1MBPTb files into E1.RES and E1MBPT.RES

    Args:
        raw_path: Path to raw data (E1a.RES, E1b.RES, etc.)
        filtered_path: Path to filtered data (output E1.RES, E1MBPT.RES)
    """
    # Track E1 matrix elements separately from E1MBPT
    e1_res = []
    e1_mbpt_res = []

    # Ensure filtered path exists
    os.makedirs(filtered_path, exist_ok=True)

    # Combine E1a.RES and E1b.RES into E1.RES
    with open(raw_path + '/E1a.RES', 'r') as f:
        lines = f.readlines()
    with open(filtered_path + '/E1.RES', 'w') as f:
        for line in lines:
            f.write(line)

    for line in lines[1:]:
        matrix_element = re.findall(r'\<.*?\>', line)[0]
        e1_res.append(matrix_element)

    with open(raw_path + '/E1b.RES', 'r') as f:
        lines2 = f.readlines()

    with open(filtered_path + '/E1.RES', 'a') as f:
        for line in lines2[1:]:
            matrix_element = re.findall(r'\<.*?\>', line)[0]
            if (matrix_element) not in e1_res:
                f.write(line)

    # Combine E1MBPTa.RES and E1MBPTb.RES into E1MBPT.RES (if they exist)
    # Use separate list to avoid mixing with E1 data
    if os.path.exists(raw_path + '/E1MBPTa.RES') and os.path.exists(raw_path + '/E1MBPTb.RES'):
        with open(raw_path + '/E1MBPTa.RES', 'r') as f:
            lines = f.readlines()
        with open(filtered_path + '/E1MBPT.RES', 'w') as f:
            for line in lines:
                f.write(line)

        for line in lines[1:]:
            matrix_element = re.findall(r'\<.*?\>', line)[0]
            e1_mbpt_res.append(matrix_element)

        with open(raw_path + '/E1MBPTb.RES', 'r') as f:
            lines2 = f.readlines()

        with open(filtered_path + '/E1MBPT.RES', 'a') as f:
            for line in lines2[1:]:
                matrix_element = re.findall(r'\<.*?\>', line)[0]
                if (matrix_element) not in e1_mbpt_res:
                    f.write(line)
        print(f'E1.RES and E1MBPT.RES written to {filtered_path}')
    else:
        print(f'E1.RES written to {filtered_path}')
        print('E1MBPTa.RES and/or E1MBPTb.RES not found, skipping E1MBPT.RES creation')

if __name__ == "__main__":
    use_config_yml = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input('Using a config.yml file? ')))))
    
    if use_config_yml:
        config_yml = input("Input yml-file: ")
        config = read_yaml(config_yml)
        atom_name = get_dict_value(config['atom'],'name')
        
        # Parse atom name - keep in Roman numeral format for NIST
        # NIST uses: Ba I (neutral), Ba II (singly ionized), etc.
        if len(atom_name.split()) == 1:
            # Just element symbol, assume neutral (I)
            atom = atom_name + ' I'
        elif len(atom_name.split()) == 2:
            element = atom_name.split()[0]
            suffix = atom_name.split()[1]

            if '+' in suffix:
                # Convert + notation to Roman numerals for NIST
                # Ba+ → Ba II, Ba++ → Ba III
                ionization = suffix.count('+') + 1
                atom = element + ' ' + convert_num_to_roman(ionization)
            elif suffix.isnumeric():
                # Convert numeric to Roman numerals for NIST
                # Ba1 → Ba I, Ba2 → Ba II
                ionization = int(suffix)
                atom = element + ' ' + convert_num_to_roman(ionization)
            else:
                # Already in Roman numeral format, keep as-is
                atom = atom_name
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
        
        # set default ignore configurations with 'g' and terms with 'G' to True
        ignore_g_value = get_dict_value(portal, 'ignore_g') if portal else None
        ignore_g = ignore_g_value if ignore_g_value is not None else True

        # minimum uncertainty as percentage of value
        min_unc_value = get_dict_value(portal, 'min_uncertainty') if portal else None
        if min_unc_value is not None:
            min_uncertainty = float(min_unc_value)
        else:
            min_uncertainty = float(input('min_uncertainty not found in config.yml. Enter minimum uncertainty (as % of value): '))

        # set default minimum energy difference percentage between NIST and theory to 3.0
        min_diff_value = get_dict_value(portal, 'min_energy_diff_percent') if portal else None
        min_energy_diff_percent = float(min_diff_value) if min_diff_value is not None else 3.0

        # optional global energy cutoff in cm^-1 (if not specified, use min_energy_diff_percent logic)
        cutoff_value = get_dict_value(portal, 'energy_cutoff') if portal else None
        energy_cutoff = float(cutoff_value) if cutoff_value is not None else None
    else:
        atom = input('Input name of atom: ')
        even_dir = None
        odd_dir = None
        tm_dir = None
        tm_dir1 = None
        tm_dir2 = None
        ignore_g = True
        min_uncertainty = float(input('Enter minimum uncertainty for matrix elements (as % of value): '))
        min_energy_diff_percent = 3.0
        energy_cutoff = None
    name = atom_name_to_filename(atom)
    
    ri = False # 
    fac = 2 # maximum energy difference (in percent) for comparison
    
    # Find input files from directories if they exist and put into DATA_RAW directory
    dir_path = os.getcwd()
    data_raw_path = 'DATA_RAW'
    data_filtered_path = 'DATA_Filtered/UD/'
    data_processed_path = 'DATA_Processed/'
    if not os.path.isdir(data_raw_path):
        Path(data_raw_path).mkdir(parents=True, exist_ok=True)
    if not os.path.isdir(data_filtered_path):
        Path(data_filtered_path).mkdir(parents=True, exist_ok=True)
    if not os.path.isdir(data_processed_path):
        Path(data_processed_path).mkdir(parents=True, exist_ok=True)
    
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

    # Combine E1a/E1b and E1MBPTa/E1MBPTb files if they exist in DATA_RAW/
    if os.path.exists(data_raw_path + '/E1a.RES') and os.path.exists(data_raw_path + '/E1b.RES'):
        combine_tm(data_raw_path, data_processed_path)
    else:
        print('E1a.RES and/or E1b.RES not found in DATA_RAW/, skipping combine_tm')
        
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
        print('The files should be named: CONFFINALeven.RES, CONFFINALodd.RES, CONFFINALevenMBPT.RES, CONFFINALoddMBPT.RES, E1a.RES, E1b.RES, E1MBPTa.RES, E1MBPTb.RES')
        print('Note: E1.RES and E1MBPT.RES will be generated and placed in ' + data_processed_path)
        sys.exit()
    confs, terms, energies_au, energies_cm, uncertainties, theory_shift, theory_J, gs_parity, matrix_file_exists, gs_exists, swaps, fixes, unmatched_energies, energy_to_level, mbpt_energy_to_level = write_new_conf_res(name, raw_path, data_nist)

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

    df_to_csv(data_nist,"DATA_Filtered/NIST/"+name,'odd')
    df_to_csv(data_nist,"DATA_Filtered/NIST/"+name,'even')
    df_to_csv(data_nist,"DATA_Filtered/NIST/"+name)
    
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

    # Calculate energy offsets to reference all energies from ground state
    # NIST energies use NIST_shift, theory energies use theory_shift
    # The parity opposite to ground state needs to be shifted
    nist_even_offset = 0.0 if gs_parity == 'even' else NIST_shift
    nist_odd_offset = 0.0 if gs_parity == 'odd' else NIST_shift
    theory_even_offset = 0.0 if gs_parity == 'even' else theory_shift
    theory_odd_offset = 0.0 if gs_parity == 'odd' else theory_shift

    print(f'Ground state parity: {gs_parity}')
    print(f'NIST energy offset for even parity: {nist_even_offset} cm^-1')
    print(f'NIST energy offset for odd parity: {nist_odd_offset} cm^-1')
    print(f'Theory energy offset for even parity: {theory_even_offset} cm^-1')
    print(f'Theory energy offset for odd parity: {theory_odd_offset} cm^-1')

    # Use Vipul's code to correct misidentified configurations
    data_final_even = MainCode(path_nist_even, path_ud_even, nist_max_even, gs_exists, nist_offset=nist_even_offset, theory_offset=theory_even_offset)
    data_final_odd = MainCode(path_nist_odd, path_ud_odd, nist_max_odd, gs_exists, nist_offset=nist_odd_offset, theory_offset=theory_odd_offset)
    
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
    
    # Filter mapping: truncate at first level exceeding energy cutoff
    # Apply separately to even and odd parity since they're ordered independently
    even_mapping = mapping[:num_levels_output_even]
    odd_mapping = mapping[num_levels_output_even:]

    # Determine global cutoff energy
    if energy_cutoff is not None:
        # User specified an explicit energy cutoff in cm^-1
        global_cutoff_energy = energy_cutoff
        print(f'Using user-specified global energy cutoff: {global_cutoff_energy:.2f} cm^-1')
    else:
        # Calculate cutoff from min_energy_diff_percent
        # Find truncation point for even parity
        even_truncate_idx = len(even_mapping)
        even_cutoff_energy = float('inf')
        for i, level in enumerate(even_mapping):
            if level[2] > min_energy_diff_percent:
                even_truncate_idx = i
                # Get the energy from the LAST GOOD level (i-1), not the first bad one
                if i > 0:
                    last_good_level = even_mapping[i-1]
                    if last_good_level[0][3] != '-':
                        even_cutoff_energy = float(last_good_level[0][3])
                    else:
                        even_cutoff_energy = float(last_good_level[1][3])
                    print(f'Even parity: truncating at level {i} (energy diff: {level[2]:.2f}% > {min_energy_diff_percent}%, last good energy: {even_cutoff_energy:.2f} cm^-1)')
                else:
                    # First level already exceeds threshold
                    even_cutoff_energy = 0.0
                    print(f'Even parity: first level already exceeds threshold (energy diff: {level[2]:.2f}% > {min_energy_diff_percent}%)')
                break

        # Find truncation point for odd parity
        odd_truncate_idx = len(odd_mapping)
        odd_cutoff_energy = float('inf')
        for i, level in enumerate(odd_mapping):
            if level[2] > min_energy_diff_percent:
                odd_truncate_idx = i
                # Get the energy from the LAST GOOD level (i-1), not the first bad one
                if i > 0:
                    last_good_level = odd_mapping[i-1]
                    if last_good_level[0][3] != '-':
                        odd_cutoff_energy = float(last_good_level[0][3])
                    else:
                        odd_cutoff_energy = float(last_good_level[1][3])
                    print(f'Odd parity: truncating at level {i} (energy diff: {level[2]:.2f}% > {min_energy_diff_percent}%, last good energy: {odd_cutoff_energy:.2f} cm^-1)')
                else:
                    # First level already exceeds threshold
                    odd_cutoff_energy = 0.0
                    print(f'Odd parity: first level already exceeds threshold (energy diff: {level[2]:.2f}% > {min_energy_diff_percent}%)')
                break

        # Apply global cutoff: use the lower of the two cutoff energies to prevent gaps
        global_cutoff_energy = min(even_cutoff_energy, odd_cutoff_energy)

    if global_cutoff_energy < float('inf'):
        if energy_cutoff is None:
            print(f'\nApplying global cutoff at {global_cutoff_energy:.2f} cm^-1 (lower of even/odd cutoffs)')
        else:
            print(f'\nApplying global cutoff at {global_cutoff_energy:.2f} cm^-1')

        # Re-truncate even parity based on global cutoff
        even_truncate_idx = len(even_mapping)
        for i, level in enumerate(even_mapping):
            if level[0][3] != '-':
                energy = float(level[0][3])
            else:
                energy = float(level[1][3])

            if energy > global_cutoff_energy:
                even_truncate_idx = i
                print(f'  Even parity: truncating at level {i} (energy {energy:.2f} > {global_cutoff_energy:.2f} cm^-1)')
                break

        # Re-truncate odd parity based on global cutoff
        odd_truncate_idx = len(odd_mapping)
        for i, level in enumerate(odd_mapping):
            if level[0][3] != '-':
                energy = float(level[0][3])
            else:
                energy = float(level[1][3])

            if energy > global_cutoff_energy:
                odd_truncate_idx = i
                print(f'  Odd parity: truncating at level {i} (energy {energy:.2f} > {global_cutoff_energy:.2f} cm^-1)')
                break

    # Create filtered mapping with truncated levels
    filtered_mapping = even_mapping[:even_truncate_idx] + odd_mapping[:odd_truncate_idx]
    if energy_cutoff is not None:
        print(f'Filtered {len(mapping)} levels to {len(filtered_mapping)} levels using energy cutoff: {energy_cutoff:.2f} cm^-1')
    else:
        print(f'Filtered {len(mapping)} levels to {len(filtered_mapping)} levels using min energy diff: {min_energy_diff_percent}%')
    print(f'  Even parity: {len(even_mapping)} -> {even_truncate_idx}')
    print(f'  Odd parity: {len(odd_mapping)} -> {odd_truncate_idx}')

    # Filter out states with 'g' in configuration or 'G' in term if ignore_g is True
    if ignore_g:
        before_g_filter = len(filtered_mapping)
        filtered_mapping_no_g = []
        for level in filtered_mapping:
            # Check both NIST and theory configurations/terms
            nist_config = level[0][0]
            nist_term = level[0][1]
            theory_config = level[1][5]  # corrected config
            theory_term = level[1][1]

            # Skip if 'g' in any configuration or 'G' in any term
            if 'g' in nist_config or 'g' in theory_config or 'G' in nist_term or 'G' in theory_term:
                continue
            filtered_mapping_no_g.append(level)

        filtered_mapping = filtered_mapping_no_g
        print(f'Filtered out {before_g_filter - len(filtered_mapping)} states with g/G (ignore_g=True)')
        print(f'Final mapping has {len(filtered_mapping)} levels')

    # Generate mapping-based fixes: when theory_config differs from corrected_config
    # These fixes change E1.RES entries to match the NIST-determined labels
    mapping_fixes = generate_mapping_fixes(filtered_mapping)
    if mapping_fixes:
        # Add mapping fixes to both fixes1 (for E1.RES) and fixes2 (for E1MBPT.RES)
        fixes1, fixes2 = fixes
        fixes = (fixes1 + mapping_fixes, fixes2 + mapping_fixes)

    write_energy_csv(name, filtered_mapping, NIST_shift, theory_shift, gs_parity, min_energy_diff_percent)

    # Create a list of all possible transitions between states
    print('even parity configurations:')
    even_confs = []
    odd_confs = []
    for line in filtered_mapping:
        # Use NIST config if available, otherwise use theory config
        if line[0][0] != '-':
            configuration = line[0][0]
            term = line[0][1]
            J = line[0][2]
        else:
            configuration = line[1][5]
            term = line[1][1]
            J = line[1][2]

        if find_parity(configuration) == 'even':
            even_confs.append([configuration, term, J])
        elif find_parity(configuration) == 'odd':
            odd_confs.append([configuration, term, J])
        else:
            raise ValueError(f'Configuration {configuration} is not valid.')

    possible_E1 = []
    for conf_odd in odd_confs:
        try:
            float(conf_odd[2])
            J_odd = int(conf_odd[2])
        except ValueError:
            J_odd = int(conf_odd[1])
        for conf_even in even_confs:
            try:
                float(conf_even[2])
                J_even = int(conf_even[2])
            except ValueError:
                J_even = int(conf_even[1])
            if J_even == 0 and J_odd == 0: continue
            if abs(J_even - J_odd) <= 1:
                possible_E1.append([conf_odd, conf_even])
    num_possible_E1 = len(possible_E1)
    print('Number of possible E1: ', len(possible_E1))

    unmatched_matrix = []
    if matrix_file_exists:
        print('Writing matrix elements...')
        num_E1, unmatched_matrix = write_matrix_csv(name, path_filtered_theory, filtered_mapping, gs_parity, theory_shift, NIST_shift, swaps, fixes, ignore_g, min_uncertainty, min_energy_diff_percent, energy_to_level, mbpt_energy_to_level)
        coverage = round(num_E1/num_possible_E1*100, 2)
        print(f'{coverage}% of possible E1 transitions accounted for ({num_E1}/{num_possible_E1})')
    else:
        print('E1.RES files were not found, so matrix csv file was not generated')

    # Write unmatched items to file
    write_unmatched_file(unmatched_energies, unmatched_matrix)