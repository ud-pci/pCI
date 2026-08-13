"""
gen_portal_csv.py -- generate atomic data CSV files for the ATOM Portal from pCI calculation results and available NIST data.

Pipeline overview
-----------------
1. Read config YAML: atom name, J values for even/odd parity, and portal options (min_uncertainty, min_energy_diff_percent, energy_cutoff, gauge).

2. Collect all input files:
    a. Copy pconf.csv and tm.csv from ci+all-order/ and ci+second-order/ to DATA_RAW/ ("_MBPT" is appended to CI+MBPT results).
    b. Fetch NIST ASD data live from the NIST website using the atom name.

3. Process all input data:
    a. process_pconf_levels: read pconf.csv for both CI+all-order and CI+MBPT.
       Match levels by (conf, term) and compute level energy uncertainty as |CI+all-order - CI+MBPT| in cm^-1,
       both as absolute excitation from the ground state. Write merged CSV files to DATA_Filtered/UD/.
    b. combine_tm: merge the CI+all-order and CI+MBPT tm.csv files and write DATA_Processed/tm.csv with matrix element uncertainty.
    c. Reformat NIST data and write parity-filtered CSV files to DATA_Filtered/NIST/.

4. Run MainCode (compare_res.py) for each parity to match theory levels to NIST ASD levels.
    Results are written to DATA_Output/{atom}_Even.txt and DATA_Output/{atom}_Odd.txt.

5. create_mapping: read the txt files and build a list of mapping tuples:
    [ [NIST_config, NIST_term, NIST_J, NIST_energy, NIST_uncertainty],
      [theory_config, theory_term, theory_J, theory_energy_cm, theory_uncertainty, corrected_config, theory_energy_au],
      energy_diff_pct,   # |NIST - theory| / NIST * 100 (%)
      has_nist ]         # True if a NIST match was found

6. Filter mapping:
    - Drop levels with no MBPT match.
    - Apply energy ceiling.
    - Drop levels where has_nist=True but NIST-theory energy agreement exceeds the threshold.
    - Optionally drop g-orbital levels.

7. write_energy_csv: write {atom}_Energies.csv.
    - If has_nist=True: use NIST config/term/J and NIST energy.
    - If has_nist=False (is_from_theory=True): use corrected theory labels and theory energy.

8. write_matrix_csv: write {atom}_Matrix_Elements_Theory.csv.
    - Match each tm.csv row to the filtered mapping via theory_energy_au (tolerance 1e-6 au).
    - Substitute NIST labels for matched levels (has_nist=True).
    - Apply minimum percentage uncertainty floor in quadrature.
    - Round matrix element and uncertainty to 5 decimal places; floor uncertainty at 1e-05 if rounding produces zero.
"""
import yaml
import re
import sys
import os
import numpy as np
import pandas as pd
from fractions import Fraction
from UDRead import *
from parse_asd import *
from get_atomic_term import *
from compare_res import *
from utils import get_dict_value
import shutil

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

def create_mapping(name, num_levels_even, num_levels_odd):
    '''
    This function reads Vipul's energy level table and creates the mapping between experimental and theory data
    Data to map: Config, Term, J, Energy(cm-1)
    '''
    with open(f'DATA_Output/{name}_Even.txt', 'r') as f:
        lines = f.readlines()[:num_levels_even + 1]

    with open(f'DATA_Output/{name}_Odd.txt', 'r') as f:
        lines = lines + f.readlines()[:num_levels_odd + 1]
    
    mapping = []
    for line in lines[1:]:
        NIST_config = line.split()[0]
        NIST_term = line.split()[1]
        NIST_J = line.split()[2]
        NIST_energy = line.split()[3]
        NIST_uncertainty = line.split()[4]
        theory_config = line.split()[6]
        corrected_config = line.split()[5]
        theory_term = line.split()[8]
        try:
            theory_J = str(Fraction(line.split()[9]))
        except ValueError:
            theory_J = line.split()[9]
        theory_energy_cm = line.split()[10]
        theory_uncertainty = line.split()[11]
        theory_energy_au = line.split()[12]
        try:
            # Extract energy difference percentage (last column, remove '%' sign)
            energy_diff_pct = float(line.split()[14].rstrip('%'))
        except (ValueError, IndexError):
            energy_diff_pct = 0.0

        # select relevant data for portal database
        if NIST_config != 'Config':
            has_nist = NIST_config != '-'
            mapping.append([[NIST_config, NIST_term, NIST_J, NIST_energy, NIST_uncertainty],
                    [theory_config, theory_term, theory_J, theory_energy_cm, theory_uncertainty, corrected_config, theory_energy_au],
                    energy_diff_pct, has_nist])
    
    return mapping

def generate_mapping_fixes(mapping):
    """
    Find levels where NIST matching assigned a different label than the original theory label.
    Returns list of [old_conf, old_term, energy_au, new_conf, new_term] entries.
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


def read_pconf_csv(path):
    '''
    Read a pconf.csv file and return a list of rows.
    Each row is: [conf, term, energy_au, energy_cm, S, L, J_str, gf, conf%, converged, conf2, conf2%]
    Convert spaces in conf and conf2 to dots for consistency with NIST.
    '''
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        J_str    = _j_suffix(float(r['J']))
        term_raw = '' if pd.isna(r['term']) else str(r['term'])
        S_val    = '' if pd.isna(r['S'])      else str(r['S'])
        L_val    = '' if pd.isna(r['L'])      else str(r['L'])
        gf_val   = '' if pd.isna(r['gf'])     else str(r['gf'])
        conf2    = '' if pd.isna(r['conf2'])   else str(r['conf2'])
        conf2pct = '' if pd.isna(r['conf2%'])  else str(r['conf2%'])
        conf_str  = str(r['configuration']).replace(' ', '.')
        conf2_str = conf2.replace(' ', '.')
        rows.append([conf_str, term_raw,
                     float(r['energy_au']), float(r['energy_cm']),
                     S_val, L_val, J_str, gf_val,
                     str(r['conf%']), str(r['converged']), conf2_str, conf2pct])
    return rows


def write_pconf_csv(name, parity, ao_rows, level_rows):
    '''
    Write a merged pconf.csv file with uncertainties for one parity.
    ao_rows: output of read_pconf_csv
    level_rows: output of _build_levels
    '''
    parity_cap = parity.capitalize()
    csvfile = f'DATA_Filtered/UD/{name}_UD_{parity_cap}.csv'
    os.makedirs(os.path.dirname(csvfile), exist_ok=True)
    with open(csvfile, 'w') as f:
        f.write('n, conf, term, E_n (a.u.), DEL (cm^-1), S, L, J, gf, conf%, converged, conf2, conf2%, uncertainty \n')
        for idx, (r, lv) in enumerate(zip(ao_rows, level_rows), 1):
            conf, term = r[0], r[1]
            energy_au = r[2]
            S, L, J_str = r[4], r[5], r[6]
            gf, conf_pct, converged = r[7], r[8], r[9]
            conf2, conf2pct = r[10], r[11]
            energy_cm = lv[3]
            uncertainty = lv[4]
            S_fmt = '{:.3f}'.format(float(S)) if S else ''
            L_fmt = '{:.3f}'.format(float(L)) if L else ''
            gf_fmt = '{:.5f}'.format(float(gf)) if gf else ''
            cpct_fmt = '{:.2f}'.format(float(conf_pct)) if conf_pct else ''
            c2pct_fmt = '{:.2f}'.format(float(conf2pct)) if conf2pct else ''
            f.write(','.join([str(idx), conf, term, str(energy_au), str(energy_cm),
                              S_fmt, L_fmt, J_str, gf_fmt, cpct_fmt, str(converged),
                              conf2, c2pct_fmt, str(uncertainty)]) + '\n')
    print(f'{csvfile} has been written')


_L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6, 'K': 7}
_term_re = re.compile(r'^(\d+)([A-Z])(\d+)$')
_term_cache = {}

def fix_term(conf, term, S_val, L_val):
    """
    If term is not a valid LS term for conf, replace it with the valid term 
    whose S and L are closest to the CI expectation values S_val/L_val.
    Returns the original term unchanged if no better match is found.
    """
    conf_space = conf.replace('.', ' ')
    if conf_space not in _term_cache:
        _term_cache[conf_space] = scrape_term(conf_space)
    valid = _term_cache[conf_space]
    if not valid or term in valid:
        return term
    m = _term_re.match(term)
    if not m:
        return term
    J_target = int(m.group(3))
    valid_same_J = [(t, mt) for t in valid for mt in [_term_re.match(t)] if mt and int(mt.group(3)) == J_target]
    if not valid_same_J:
        return term
    try:
        S_ci = float(S_val)
        L_ci = float(L_val)
    except (ValueError, TypeError):
        return term
    best, best_dist = term, float('inf')
    for t, mt in valid_same_J:
        S_t = (int(mt.group(1)) - 1) / 2.0
        L_t = _L_map.get(mt.group(2), -1)
        if L_t < 0:
            continue
        dist = abs(S_ci - S_t) + abs(L_ci - L_t)
        if dist < best_dist:
            best_dist = dist
            best = t
    if best != term:
        print(f'  fix_term: {conf} {term} -> {best} (S={float(S_val):.3f}, L={float(L_val):.3f})')
    return best

def process_pconf_levels(name, filepath):
    '''
    Read pconf_even.csv / pconf_odd.csv (and optional MBPT variants) from filepath,
    compute level energies globally referenced to the ground state, 
    estimate theory uncertainties as |CI+all-order - CI+MBPT| in cm^-1, 
    and write the final CSV files for both parities.
    '''
    
    second_order_exists = True

    path_even = filepath + 'pconf_even.csv'
    path_odd = filepath + 'pconf_odd.csv'
    path_even_mbpt = filepath + 'pconf_even_MBPT.csv'
    path_odd_mbpt = filepath + 'pconf_odd_MBPT.csv'
    if not os.path.exists(path_even):
        print(f'{path_even} not found')
        sys.exit()
    if not os.path.exists(path_odd):
        print(f'{path_odd} not found')
        sys.exit()
    if not os.path.exists(path_even_mbpt):
        print(f'{path_even_mbpt} not found')
        second_order_exists = False
    if not os.path.exists(path_odd_mbpt):
        print(f'{path_odd_mbpt} not found')
        second_order_exists = False

    ao_even   = read_pconf_csv(path_even)
    ao_odd    = read_pconf_csv(path_odd)
    mbpt_even = read_pconf_csv(path_even_mbpt) if second_order_exists else []
    mbpt_odd  = read_pconf_csv(path_odd_mbpt)  if second_order_exists else []

    # Fix terms in all files before matching
    for rows in (ao_even, ao_odd, mbpt_even, mbpt_odd):
        for r in rows:
            r[1] = fix_term(r[0], r[1], r[4], r[5])

    # Larger energy_au = more tightly bound (valence energy convention in pconf)
    if ao_even[0][2] > ao_odd[0][2]:
        gs_parity    = 'even'
        gs_energy_au = ao_even[0][2]
    else:
        gs_parity    = 'odd'
        gs_energy_au = ao_odd[0][2]

    ht_to_cm = 219474.63

    # MBPT lookup: (conf, term) -> [(n_mbpt, row), ...] by row index order.
    # Using a list per key rather than a flat dict lets duplicates get mapped to each other.
    def _build_mbpt_groups(mbpt_rows):
        groups = {}
        for n, mr in enumerate(mbpt_rows, 1):
            key = (mr[0], mr[1])
            if key not in groups:
                groups[key] = []
            groups[key].append((n, mr))
        return groups

    mbpt_groups_even = _build_mbpt_groups(mbpt_even)
    mbpt_groups_odd = _build_mbpt_groups(mbpt_odd)
    mbpt_gs_au = (mbpt_even[0][2] if gs_parity == 'even' else mbpt_odd[0][2]) if second_order_exists else None

    def _build_levels(ao_rows, mbpt_groups, parity):
        levels = []
        occurrence = {}
        for r in ao_rows:
            conf, term, energy_au, energy_cm = r[0], r[1], r[2], r[3]
            J_str = r[6]
            if parity == gs_parity:
                energy_cm_global = energy_cm
            else:
                energy_cm_global = round((gs_energy_au - energy_au) * ht_to_cm, 2)
            key = (conf, term)
            occ = occurrence.get(key, 0)
            occurrence[key] = occ + 1
            group = mbpt_groups.get(key, [])
            if occ < len(group) and mbpt_gs_au is not None:
                _, mr = group[occ]
                ao_abs = (gs_energy_au - energy_au) * ht_to_cm
                mbpt_abs = (mbpt_gs_au - mr[2]) * ht_to_cm
                unc = round(abs(ao_abs - mbpt_abs))
            else:
                unc = '-'
            levels.append([conf, term, energy_au, energy_cm_global, unc, J_str])
        return levels

    even_levels = _build_levels(ao_even, mbpt_groups_even, 'even')
    odd_levels  = _build_levels(ao_odd,  mbpt_groups_odd,  'odd')
    all_levels  = even_levels + odd_levels

    confs         = [r[0] for r in all_levels]
    terms         = [r[1] for r in all_levels]
    energies_au   = [r[2] for r in all_levels]
    energies_cm   = [r[3] for r in all_levels]
    uncertainties = [r[4] for r in all_levels]

    with open('final_res.csv', 'w') as f:
        f.write('conf,term,J,energy(a.u.),energy(cm-1),uncertainty\n')
        for r in all_levels:
            f.write(','.join([r[0], r[1], r[5],
                              str(r[2]), str(r[3]), str(r[4])]) + '\n')

    theory_J = {
        'even': [r[6] for r in ao_even],
        'odd':  [r[6] for r in ao_odd],
    }

    def _write_pconf_comp(ao_rows, mbpt_groups, level_rows, parity):
        csvfile = f'DATA_Filtered/UD/pconf_comp_{parity}.csv'
        os.makedirs(os.path.dirname(csvfile), exist_ok=True)
        occurrence = {}
        with open(csvfile, 'w') as f:
            f.write('n_ao,n_mbpt,conf_ao,term_ao,conf_mbpt,term_mbpt,'
                    'energy_cm_ao,energy_cm_mbpt,uncertainty,conf2_ao,conf2_mbpt\n')
            for n_ao, (r, lv) in enumerate(zip(ao_rows, level_rows), 1):
                conf_ao = r[0]
                term_ao = r[1]
                conf2_ao = r[10]
                energy_cm_ao = lv[3]
                unc = lv[4]
                key = (conf_ao, term_ao)
                occ = occurrence.get(key, 0)
                occurrence[key] = occ + 1
                group = mbpt_groups.get(key, [])
                if occ < len(group):
                    n_mbpt, mr    = group[occ]
                    conf_mbpt     = mr[0]
                    term_mbpt     = mr[1]
                    conf2_mbpt    = mr[10]
                    mbpt_energy_cm = round((mbpt_gs_au - mr[2]) * ht_to_cm, 2)
                else:
                    n_mbpt = conf_mbpt = term_mbpt = conf2_mbpt = '-'
                    mbpt_energy_cm = '-'
                f.write(','.join([str(n_ao), str(n_mbpt),
                                  conf_ao, term_ao,
                                  str(conf_mbpt), str(term_mbpt),
                                  str(energy_cm_ao), str(mbpt_energy_cm),
                                  str(unc), conf2_ao, str(conf2_mbpt)]) + '\n')
        print(f'{csvfile} has been written')
    
    print(f'Ground state parity: {gs_parity}')
    
    write_pconf_csv(name, 'even', ao_even, even_levels)
    write_pconf_csv(name, 'odd',  ao_odd,  odd_levels)
    _write_pconf_comp(ao_even, mbpt_groups_even, even_levels, 'even')
    _write_pconf_comp(ao_odd,  mbpt_groups_odd,  odd_levels,  'odd')

    return (confs, terms, energies_au, energies_cm, uncertainties, theory_J)

def convert_type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    s = s.replace(" ", "")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        return s 
    
def write_final_res(full_res, outfile):
    '''
    Write FINAL.RES file with an additional uncertainty column calculated between CI+all-order and CI+MBPT.
    '''
    ht_to_cm = 219474.63
    energy0_au = full_res[0][3]

    with open(outfile, 'w') as f:
        header = f'{"n":>3s}  {"conf":>8s}   {"term":<4s}  {"E_n (a.u.)":>14s}  {"DEL (cm^-1)":>13s}  {"S":>4s}  {"L":>4s}  {"J":>4s}  {"gf":>5s}  {"conf%":>8s}  {"converged":>9s}  {"conf2":>8s}  {"conf2%":>7s}  uncertainty\n'
        f.write(header)
        for row in full_res:
            idx = row[0]
            conf = row[1].replace('.', ' ')
            term = row[2].replace(',', '')
            energy_au = row[3]
            energy_cm = ht_to_cm * (energy0_au - energy_au)
            s = row[5]
            l = row[6]
            j = row[7]
            gf = row[8]
            cntrb1 = row[9]
            converged = row[10]
            sec_conf = row[11].replace('.', ' ') if row[11] else ''
            cntrb2 = row[12] if row[12] else ''
            uncertainty = row[13] if len(row) > 13 else '-'

            # Format gf
            if gf == '-----' or gf == 0:
                gf_str = '-----'
            else:
                gf_str = f'{float(gf):.3f}'

            # Build line - always pad conf2/conf2% columns to fixed width
            line = f'{idx:3d}  {conf:>8s}   {term:<4s}  {energy_au:14.8f}  {energy_cm:13.1f}  {float(s):.2f}  {float(l):.2f}  {float(j):.2f}  {gf_str:>5s}  {cntrb1:>8s}  {str(converged):>9s}  {sec_conf:>8s}  {cntrb2:>7s}  {uncertainty}'
            f.write(line + '\n')

    print(f'{outfile} has been written')


def write_energy_csv(name, mapping, min_energy_diff_percent):
    '''
    This function writes the energy csv file
    Mapping should already be filtered to only include levels within min_energy_diff_percent
    '''
    filename = name + '_Energies.csv'

    # select columns for portal dataframe
    portal_df = pd.DataFrame(columns = ['state_configuration', 'state_term', 'state_J', 'energy',
                                     'energy_uncertainty', 'is_from_theory'])

    for level in mapping:
        is_from_theory = not level[3]
        if is_from_theory:
            state_config = level[1][5]
            state_term = level[1][1]
            state_J = level[1][2]
            state_energy = level[1][3]
            state_uncertainty = level[1][4]
        else:
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

def write_matrix_csv(element, mapping, ignore_g, min_unc_per, min_energy_diff_percent, gauge='L'):
    '''
    Read DATA_Processed/tm.csv and write {element}_Matrix_Elements_Theory.csv.
    All operators (E1, E2, E3, M1, M2, M3) are included. 
    For E1, gauge selects which form to keep: 'L' (length, default) or 'V' (velocity); 
    the redundant E1 gauge column is dropped and the kept column is renamed to 'E1'. 
    States are matched to the NIST/theory mapping by energy_au so that experimental configurations and energies are substituted where available.
    tm.csv already has calculated uncertainty from combine_tm, and zero-energy transitions already filtered.
    '''
    tm_path = 'DATA_Processed/tm.csv'
    if not os.path.isfile(tm_path):
        print(f'{tm_path} not found; skipping matrix elements')
        return

    df_tm = pd.read_csv(tm_path)
    
    # Drop the E1 gauge not selected; keep all other operators unchanged
    drop_gauge = 'E1_V' if gauge == 'L' else 'E1_L'
    df_tm = df_tm[df_tm['operator'] != drop_gauge].reset_index(drop=True)
    df_tm['operator'] = df_tm['operator'].replace('E1_' + gauge, 'E1')

    # Drop zero matrix elements
    df_tm = df_tm[df_tm['matrix_element_value'] != 0].reset_index(drop=True)
    
    matrix_element_filename = element + '_Matrix_Elements_Theory.csv'

    df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                               'state_two_configuration', 'state_two_term', 'state_two_J',
                               'operator', 'matrix_element', 'matrix_element_uncertainty'])

    added_transitions = set()

    if ignore_g:
        print('IGNORING G STATES')

    default_min_uncertainties = {
        'Mg1': 0.3,
        'Ca1': 1.3,
        'Sr1': 1.5
    }
    if element in default_min_uncertainties:
        min_unc_per = default_min_uncertainties[element]
        print(f'DEFAULT MINIMUM MATRIX ELEMENT UNCERTAINTY USED FOR {element}: {min_unc_per}%')

    print(f'MIN ENERGY DIFF PERCENT FOR NIST-THEORY: {min_energy_diff_percent}%')

    energy_tolerance = 1e-6

    for _, row in df_tm.iterrows():
        operator = str(row['operator'])
        try:
            matrix_element_value = float(row['matrix_element_value'])
            unc_raw = row['matrix_element_uncertainty']
            uncertainty = 0.0 if pd.isna(unc_raw) else float(unc_raw)
        except (ValueError, TypeError):
            continue

        energy1_au = abs(float(row['state_one_energy_au']))
        energy2_au = abs(float(row['state_two_energy_au']))

        extra_unc = abs(matrix_element_value) * min_unc_per / 100
        uncertainty = np.sqrt(uncertainty**2 + extra_unc**2)

        c1, c2 = False, False
        out_conf1, out_term1, out_J1 = '', '', ''
        out_conf2, out_term2, out_J2 = '', '', ''
        best1, best2 = float('inf'), float('inf')

        for lt in mapping:
            if lt[1][6] == '-':
                continue
            th_au = abs(float(lt[1][6]))

            diff1 = abs(th_au - energy1_au)
            if diff1 < energy_tolerance and diff1 < best1:
                cand_conf = lt[1][5]
                cand_term = lt[1][1]
                cand_J    = lt[1][2]
                if ignore_g and ('g' in cand_conf or 'G' in cand_term):
                    continue
                if lt[3]:
                    nist_e = float(lt[0][3])
                    th_e = float(lt[1][3])
                    pct = abs((nist_e - th_e) / nist_e * 100) if nist_e != 0 else 0.0
                    if pct > min_energy_diff_percent:
                        continue
                    cand_conf = lt[0][0]
                    cand_term = lt[0][1]
                    cand_J = lt[0][2]
                out_conf1, out_term1, out_J1 = cand_conf, cand_term, cand_J
                c1 = True
                best1 = diff1

            diff2 = abs(th_au - energy2_au)
            if diff2 < energy_tolerance and diff2 < best2:
                cand_conf = lt[1][5]
                cand_term = lt[1][1]
                cand_J = lt[1][2]
                if ignore_g and ('g' in cand_conf or 'G' in cand_term):
                    continue
                if lt[3]:
                    nist_e = float(lt[0][3])
                    th_e = float(lt[1][3])
                    pct = abs((nist_e - th_e) / nist_e * 100) if nist_e != 0 else 0.0
                    if pct > min_energy_diff_percent:
                        continue
                    cand_conf = lt[0][0]
                    cand_term = lt[0][1]
                    cand_J    = lt[0][2]
                out_conf2, out_term2, out_J2 = cand_conf, cand_term, cand_J
                c2 = True
                best2 = diff2

        if not (c1 and c2):
            continue

        if (out_conf1, out_term1, out_J1) == (out_conf2, out_term2, out_J2):
            continue

        # Include operator so different operators for the same state pair each get their own row
        trans_id = (tuple(sorted([out_conf1 + ' ' + out_term1 + str(out_J1),
                                  out_conf2 + ' ' + out_term2 + str(out_J2)])),
                    operator)
        if trans_id in added_transitions:
            continue
        added_transitions.add(trans_id)

        df.loc[len(df.index)] = {
            'state_one_configuration': out_conf1, 'state_one_term': out_term1, 'state_one_J': out_J1,
            'state_two_configuration': out_conf2, 'state_two_term': out_term2, 'state_two_J': out_J2,
            'operator': operator,
            'matrix_element': round(matrix_element_value, 5), 'matrix_element_uncertainty': round(uncertainty, 5) or 1e-05,
        }

    print(f'TOTAL MATRIX ELEMENTS: {len(df)}')
    df.to_csv(matrix_element_filename, index=False)
    print(f'{matrix_element_filename} has been written')

def nist_parity(term):
    ''' 
    this function determines parity of energy level based on term from NIST designated with a '*' 
    '''
    if '*' in term:
        parity = 'odd'
    else:
        parity = 'even'
    
    return parity


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

def _j_suffix(j):
    f = float(j)
    i = int(f)
    return str(i) if f == i else str(f)

def collect_portal_files(method_dir, j0, j1, data_raw_path, res_suffix=''):
    if not os.path.isdir(method_dir):
        print(f'{method_dir} directory not found, skipping')
        return
    
    s0 = str(int(float(j0)))
    s1 = str(int(float(j1)))
    sfx = ('_' + res_suffix) if res_suffix else ''
    for parity in ('even', 'odd'):
        src = method_dir + '/' + parity + s0 + '/pconf.csv'
        dst = data_raw_path + '/pconf_' + parity + sfx + '.csv'
        if os.path.isfile(src):
            shutil.copy(src, dst)
            print(f'copied {src} to {dst}')
        else:
            print(f'{src} not found, skipping')
    tm_dirs = [
        'tm_even' + s0 + '_even' + s1,
        'tm_even' + s0 + '_odd'  + s1,
        'tm_odd'  + s0 + '_odd'  + s1,
        'tm_odd'  + s0 + '_even' + s1,
    ]
    for tm_label in tm_dirs:
        src = method_dir + '/' + tm_label + '/tm.csv'
        dst = data_raw_path + '/' + tm_label + sfx + '.csv'
        if os.path.isfile(src):
            shutil.copy(src, dst)
            print(f'copied {src} to {dst}')
        else:
            print(f'{src} not found, skipping')

def combine_tm(j0, j1, data_raw_path, data_processed_path, filtered_path):
    s0 = str(int(float(j0)))
    s1 = str(int(float(j1)))
    tm_labels = [
        'tm_even' + s0 + '_even' + s1,
        'tm_even' + s0 + '_odd'  + s1,
        'tm_odd'  + s0 + '_odd'  + s1,
        'tm_odd'  + s0 + '_even' + s1,
    ]
    match_keys = ['state_one_configuration', 'state_one_term', 'state_two_configuration', 'state_two_term', 'operator']
    os.makedirs(filtered_path, exist_ok=True)
    
    frames = []
    unmatched_frames = []
    for label in tm_labels:
        path_ao   = data_raw_path + '/' + label + '.csv'
        path_mbpt = data_raw_path + '/' + label + '_MBPT.csv'
        if not os.path.isfile(path_ao):
            print(f'{path_ao} not found, skipping')
            continue
        df_ao = pd.read_csv(path_ao)
        n_unmatched = 0
        if os.path.isfile(path_mbpt):
            df_mbpt = pd.read_csv(path_mbpt)
            merged = df_ao.merge(df_mbpt, on=match_keys, suffixes=('', '_mbpt'), how='left')
            no_mbpt = merged[merged['matrix_element_value_mbpt'].isna()]
            n_unmatched = len(no_mbpt)
            if n_unmatched > 0:
                tag = no_mbpt[match_keys].copy()
                tag.insert(0, 'tm_file', label)
                unmatched_frames.append(tag)
            merged = merged[merged['matrix_element_value_mbpt'].notna()].reset_index(drop=True)
            merged['matrix_element_uncertainty'] = (
                merged['matrix_element_value'].abs() - merged['matrix_element_value_mbpt'].abs()
            ).abs()
            merged['matrix_element_value'] = merged['matrix_element_value'].abs()
            mbpt_cols = [c for c in merged.columns if c.endswith('_mbpt')]
            merged = merged.drop(columns=mbpt_cols)
            me_idx = merged.columns.get_loc('matrix_element_value') + 1
            cols = list(merged.columns)
            cols.remove('matrix_element_uncertainty')
            cols.insert(me_idx, 'matrix_element_uncertainty')
            merged = merged[cols]
        else:
            print(f'{path_mbpt} not found; no uncertainty for {label}')
            merged = df_ao.copy()
            merged['matrix_element_uncertainty'] = None
        filtered_out = filtered_path + '/' + label + '.csv'
        merged.to_csv(filtered_out, index=False)
        unmatched_str = f' ({n_unmatched} transitions unmatched)' if n_unmatched > 0 else ''
        print(f'{label}: {len(merged)} matched transitions{unmatched_str} -> {filtered_out}')
        frames.append(merged)
        
    if unmatched_frames:
        unmatched_path = filtered_path + '/tm_unmatched.csv'
        pd.concat(unmatched_frames, ignore_index=True).to_csv(unmatched_path, index=False)
        print(f'unmatched transitions written to {unmatched_path}')

    if not frames:
        print('no tm.csv files found; skipping combine_tm')
        return
    
    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined['transition_energy_cm'] != 0]
    os.makedirs(data_processed_path, exist_ok=True)
    out_path = data_processed_path + '/tm.csv'
    combined.to_csv(out_path, index=False)
    
    print(f'combined tm.csv written to {out_path} ({len(combined)} transitions)')

if __name__ == "__main__":
    use_config_yml = input('Using a config.yml file? ').strip().lower() in ('yes', 'y', 'true')
    
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
            min_uncertainty = float(input('portal.min_uncertainty not set in config.yml. Enter minimum matrix element uncertainty as a percentage of the value (e.g. 1.5 for 1.5%): '))

        # set default minimum energy difference percentage between NIST and theory to 3.0
        min_diff_value = get_dict_value(portal, 'min_energy_diff_percent') if portal else None
        min_energy_diff_percent = float(min_diff_value) if min_diff_value is not None else 3.0

        # optional global energy cutoff in cm^-1 (if not specified, use min_energy_diff_percent logic)
        cutoff_value = get_dict_value(portal, 'energy_cutoff') if portal else None
        energy_cutoff = float(cutoff_value) if cutoff_value is not None else None
        
        # E1 gauge: 'L' (length, default) or 'V' (velocity)
        gauge_value = get_dict_value(portal, 'gauge') if portal else None
        gauge = gauge_value if gauge_value in ('L', 'V') else 'L'
    else:
        atom = input('Input name of atom: ')
        even_J = None
        odd_J = None
        ignore_g = True
        min_uncertainty = float(input('Enter minimum uncertainty for matrix elements (as % of value): '))
        min_energy_diff_percent = 3.0
        energy_cutoff = None
        gauge = 'L'
    name = atom_name_to_filename(atom)
    
    # Find input files from directories if they exist and put into DATA_RAW directory
    dir_path = os.getcwd()
    data_raw_path = 'DATA_RAW'
    data_filtered_theory_path = 'DATA_Filtered/UD'
    data_filtered_nist_path = 'DATA_Filtered/NIST'
    data_processed_path = 'DATA_Processed'
    for d in [data_raw_path, data_filtered_theory_path, data_filtered_nist_path, data_processed_path]:
        os.makedirs(d, exist_ok=True)
    
    j0, j1 = None, None
    if even_J is not None and odd_J is not None:
        j_values = sorted(set([float(even_J), float(odd_J)]))
        if len(j_values) >= 2:
            j0, j1 = j_values[0], j_values[1]
            collect_portal_files('ci+all-order', j0, j1, data_raw_path)
            collect_portal_files('ci+second-order', j0, j1, data_raw_path, 'MBPT')
        else:
            print(f'even and odd J are the same ({j_values[0]}); cannot determine TM directory pairs')
    else:
        print('J values not available; skipping portal file collection')

    # Parse NIST Atomic Spectral Database for full list of energy levels
    url_nist = generate_asd_url(atom)
    print(url_nist)
    data_nist = generate_df_from_asd(url_nist)
    
    # Write new pconf.csv with uncertainties
    raw_path = dir_path + "/DATA_RAW/"
    if os.path.isdir(raw_path):
        print(f'Reading raw files from {raw_path}')
    else:
        os.makedirs(os.path.dirname(raw_path), exist_ok=True)
        print(f'Please put raw files in {raw_path}')
        print('The files should be named: pconf_even.csv, pconf_odd.csv, pconf_even_MBPT.csv, pconf_odd_MBPT.csv')
        sys.exit()
    confs, terms, energies_au, energies_cm, uncertainties, theory_J = process_pconf_levels(name, raw_path)

    if j0 is not None:
        combine_tm(j0, j1, data_raw_path, data_processed_path, data_filtered_theory_path)
    matrix_file_exists = os.path.exists(f'{data_processed_path}/tm.csv')
    if not matrix_file_exists:
        print(f'tm.csv not found in {data_processed_path}/')

    data_nist = reformat_df_to_atomdb(data_nist, theory_J)

    # Store filtered data of even or odd parity in DATA_Filtered/NIST/
    df_to_csv(data_nist, f'{data_filtered_nist_path}/{name}', 'odd')
    df_to_csv(data_nist, f'{data_filtered_nist_path}/{name}', 'even')
    df_to_csv(data_nist, f'{data_filtered_nist_path}/{name}')
    
    # Filter NIST and UD data
    path_nist_even = f'{data_filtered_nist_path}/{name}_NIST_Even.csv'
    path_ud_even   = f'{data_filtered_theory_path}/{name}_UD_Even.csv'

    path_nist_odd = f'{data_filtered_nist_path}/{name}_NIST_Odd.csv'
    path_ud_odd   = f'{data_filtered_theory_path}/{name}_UD_Odd.csv'
    
    # Set maximum number of levels to be read from NIST equal to number of levels in pconf.csv
    num_levels_theory_even = len(pd.read_csv(raw_path + 'pconf_even.csv'))
    nist_max_even = num_levels_theory_even

    num_levels_theory_odd = len(pd.read_csv(raw_path + 'pconf_odd.csv'))
    nist_max_odd = num_levels_theory_odd
    
    print(f'Number of even parity levels: {num_levels_theory_even}')
    print(f'Number of odd parity levels: {num_levels_theory_odd}')

    # Export filtered data to output directory
    path_output = "DATA_Output"
    os.makedirs(path_output, exist_ok=True)

    data_final_even = MainCode(path_nist_even, path_ud_even, nist_max_even)
    data_final_odd = MainCode(path_nist_odd, path_ud_odd, nist_max_odd)

    ConvertToTXT(data_final_even, f'{path_output}/{name}_Even.txt')
    ConvertToTXT(data_final_odd,  f'{path_output}/{name}_Odd.txt')

    ## Finding Missing Levels
    data_final_even_missing = Missing_Levels(data_final_even)
    data_final_odd_missing = Missing_Levels(data_final_odd)

    ConvertToTXT(data_final_even_missing, f'{path_output}/{name}_Even+missing.txt')
    ConvertToTXT(data_final_odd_missing,  f'{path_output}/{name}_Odd+missing.txt')
    
    mapping = create_mapping(name, num_levels_theory_even, num_levels_theory_odd)

    # Split by parity then drop levels with no MBPT uncertainty
    even_mapping = [lv for lv in mapping[:num_levels_theory_even] if lv[1][4] != '-']
    odd_mapping  = [lv for lv in mapping[num_levels_theory_even:] if lv[1][4] != '-']

    # Determine global energy cutoff in cm^-1.
    # Exclude levels whose energy exceeds the cutoff from the outputs.
    # If energy_cutoff is not set, find the first level in each parity whose NIST-theory percentage difference exceeds min_energy_diff_percent, 
    # then place the ceiling 1 cm^-1 below that level's energy.
    # The lower of the two parities is used so even and odd coverage stay consistent.
    # global_cutoff_energy = float('inf') when no ceiling applies.
    if energy_cutoff is not None:
        global_cutoff_energy = energy_cutoff
        print(f'Using user-specified global energy cutoff: {global_cutoff_energy:.2f} cm^-1')
    else:
        global_cutoff_energy = float('inf')
        for parity_name, pmap in (('Even', even_mapping), ('Odd', odd_mapping)):
            for i, level in enumerate(pmap):
                if level[2] > min_energy_diff_percent:
                    bad_e = float(level[0][3]) if level[3] else float(level[1][3])
                    cutoff = bad_e - 1.0
                    print(f'{parity_name} parity: first bad match at level {i} (energy diff: {level[2]:.2f}% > {min_energy_diff_percent}%, cutoff: {cutoff:.2f} cm^-1)')
                    global_cutoff_energy = min(global_cutoff_energy, cutoff)
                    break
        if global_cutoff_energy < float('inf'):
            print(f'\nApplying global cutoff at {global_cutoff_energy:.2f} cm^-1 (lower of even/odd cutoffs)')

    # Keep a level if: 
    # (1) its energy is at or below the cutoff, and
    # (2) it either has no NIST match or its NIST-theory agreement is good enough.
    filtered_even, filtered_odd = [], []
    for level in even_mapping:
        energy = float(level[0][3]) if level[3] else float(level[1][3])
        if energy <= global_cutoff_energy and not (level[3] and level[2] > min_energy_diff_percent):
            filtered_even.append(level)
    for level in odd_mapping:
        energy = float(level[0][3]) if level[3] else float(level[1][3])
        if energy <= global_cutoff_energy and not (level[3] and level[2] > min_energy_diff_percent):
            filtered_odd.append(level)
    excluded_even = len(even_mapping) - len(filtered_even)
    excluded_odd = len(odd_mapping) - len(filtered_odd)

    # Create filtered mapping
    filtered_mapping = filtered_even + filtered_odd
    print(f'Filtered {len(mapping)} levels to {len(filtered_mapping)} levels')
    print(f'  Even parity: {len(even_mapping)} -> {len(filtered_even)} (excluded {excluded_even})')
    print(f'  Odd parity: {len(odd_mapping)} -> {len(filtered_odd)} (excluded {excluded_odd})')

    # Filter out states with 'g' in configuration if ignore_g is True
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
            if 'g' in nist_config or 'g' in theory_config:
                continue
            filtered_mapping_no_g.append(level)

        filtered_mapping = filtered_mapping_no_g
        print(f'Filtered out {before_g_filter - len(filtered_mapping)} states with g orbital (ignore_g=True)')
        print(f'Final mapping has {len(filtered_mapping)} levels')

    write_energy_csv(name, filtered_mapping, min_energy_diff_percent)
    
    if matrix_file_exists:
        print('Writing matrix elements...')
        write_matrix_csv(name, filtered_mapping, ignore_g, min_uncertainty, min_energy_diff_percent, gauge=gauge)
    else:
        print('tm.csv not found; matrix csv file was not generated')