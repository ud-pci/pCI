import re
import math
import os
import sys

'''
This script combines two CONFFINAL.RES files, e.g. ci+all-order and ci+mbpt results
'''

def fix_skipped_config_levels(conf_res, confs_terms_old=None):
    """
    Fix configuration labels when theory skips a principal quantum number.

    For example, if theory has:
      - 4s.5g (primary), 4s.6g (secondary)
      - 4s.7g (primary), 4s.6g (secondary)  <- should be relabeled as 4s.6g
      - 4s.8g (primary)                     <- should be relabeled as 4s.7g

    This happens because the secondary config has higher percentage than expected,
    indicating the levels are mislabeled. This can occur for any orbital type,
    typically with higher principal quantum numbers.

    Configurations can have 2 or more orbitals, e.g., '4s.7g' or '4s.4p.7g'.
    The last orbital is the one we check for skipped levels.

    Returns:
        tuple: (conf_res, skipped_fixes) where skipped_fixes is a list of
               [old_conf, old_term, energy, new_conf, new_term] for fixing matrix elements
    """
    skipped_fixes = []
    # Group levels by orbital type (last character of config, e.g., 'g' in '4s.7g')
    # and by the inner orbitals (e.g., '4s' in '4s.7g' or '4s.4p' in '4s.4p.7g')
    orbital_groups = {}

    for idx, level in enumerate(conf_res):
        config = level[1]  # primary configuration
        sec_config = level[11] if level[11] else ''  # secondary configuration

        # Parse the configuration and split into parts
        parts = config.split('.')
        if len(parts) < 2:
            continue

        # Inner orbitals are all orbitals up to the last, outer orbital is the last
        inner_orbitals = '.'.join(parts[:-1])
        outer_orbital = parts[-1]

        # Extract the principal quantum number and orbital type from outer orbital
        match = re.match(r'(\d+)([a-z])', outer_orbital)
        if not match:
            continue

        n = int(match.group(1))
        orbital_type = match.group(2)

        key = (inner_orbitals, orbital_type)
        if key not in orbital_groups:
            orbital_groups[key] = []

        # Parse secondary config if present - check if it has same inner orbitals and orbital type
        sec_n = None
        if sec_config:
            sec_parts = sec_config.split('.')
            if len(sec_parts) >= 2:
                sec_inner = '.'.join(sec_parts[:-1])
                sec_outer = sec_parts[-1]
                sec_match = re.match(r'(\d+)([a-z])', sec_outer)
                if sec_match and sec_inner == inner_orbitals and sec_match.group(2) == orbital_type:
                    sec_n = int(sec_match.group(1))

        orbital_groups[key].append({
            'idx': idx,
            'n': n,
            'sec_n': sec_n,
            'config': config,
            'sec_config': sec_config
        })

    # For each orbital group, detect and fix skipped levels
    for (inner_orbitals, orbital_type), levels in orbital_groups.items():
        # Sort by principal quantum number
        levels.sort(key=lambda x: x['n'])

        if not levels:
            continue

        # Get unique n values in sorted order
        n_values = sorted(set(level['n'] for level in levels))

        # Check for gaps where a level's secondary config fills it
        skip_detected = False
        skip_at_n = None

        for i in range(len(n_values) - 1):
            current_n = n_values[i]
            next_n = n_values[i + 1]

            # If there's a gap (next_n > current_n + 1)
            if next_n > current_n + 1:
                gap_n = current_n + 1  # The skipped n value

                # Check if any level with next_n has secondary config at gap_n
                for level_info in levels:
                    if level_info['n'] == next_n and level_info['sec_n'] == gap_n:
                        skip_detected = True
                        skip_at_n = gap_n
                        break

                if skip_detected:
                    break

        if not skip_detected:
            continue

        # Calculate how much to shift each level
        # All levels with n > skip_at_n should be decremented by 1
        for level_info in levels:
            if level_info['n'] > skip_at_n:
                idx = level_info['idx']
                old_n = level_info['n']
                new_n = old_n - 1
                old_config = conf_res[idx][1]
                term_with_comma = conf_res[idx][2]
                new_config = inner_orbitals + '.' + str(new_n) + orbital_type

                # Use original term for matrix element matching (E1.RES has original terms)
                if confs_terms_old and idx < len(confs_terms_old):
                    original_term = confs_terms_old[idx][1]  # Original term before fixes
                else:
                    original_term = term_with_comma

                # Get energy for exact matching
                energy = conf_res[idx][3]

                # Add to fixes list for matrix element correction
                # Format: [old_conf, old_term, energy, new_conf, new_term]
                # old_term uses original term to match E1.RES, new_term also uses original since term doesn't change
                skipped_fixes.append([old_config, original_term, energy, new_config, original_term])

                # Update the configuration
                old_sec = conf_res[idx][11]
                conf_res[idx][1] = new_config

                # If the old secondary matches the new primary, swap
                if old_sec == new_config:
                    conf_res[idx][11] = old_config

    return conf_res, skipped_fixes

def parse_final_res(filename):
    try:
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
    except:
        print(filename + ' not found')
        return []

    print('========== PARSING', filename,'==========')
    conf_res = []

    ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
    Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    
    # Get 1st energy level
    ht_to_cm = 219474.63
    energy0_au = float([num for num in lines[1].split('  ') if '.' in num][0])
    
    for line in lines[1:]:
        index = int(re.findall(r"\d+",line)[0])
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        confs = [conf.strip() for conf in confs]
        confs = [conf.replace(' ', '.') for conf in confs]
        
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]

        main_conf = confs[0]
        term = terms[0].replace(' ', '')
        
        energies_au = float(nums[0])
        energies_cm = ht_to_cm*(energy0_au-energies_au)
        
        s = float(nums[2])
        l = float(nums[3])
        j = int(float(nums[4]))
        if '%' in nums[5]:
            gf = '-----'
        else:
            gf = float(nums[5])
        
        cntrb1 = re.findall(r"\d+.\d+%",line)[0]
        try:
            cntrb2 = re.findall(r"\d+.\d+%",line)[1]
        except:
            cntrb2 = ''
        
        if 'True' in line:
            converged = 'True'
        else:
            converged = 'False'
        
        try:
            sec_conf = confs[1]
        except:
            sec_conf = ''
                
        conf_res.append([index, main_conf, term, energies_au, energies_cm, s, l, j, gf, cntrb1, converged, sec_conf, cntrb2])

    # Count number of electrons
    even = False
    num_electrons = 0
    for orbital in conf_res[0][1].split('.'):
        if orbital[-1].isnumeric():
            num_electrons += int(orbital[-1])
        else:
            num_electrons += 1
    if num_electrons % 2 == 0:
        even = True
    
    # Save copy of old confs and terms
    confs_terms_old = [[level[1],level[2][0:2] + ',' + level[2][2]] for level in conf_res]
    
    # Separate terms and J
    fixes = []
    i = 1
    for level in conf_res:
        # Fix term if 2 appears for even number of electrons
        conf = level[1]
        term = level[2]
        energy = level[3]
        s = level[5]
        old_term = term[0:2] + ',' + term[2]
        if even and term[0] == '2':
            if float(s) > 0.5:
                new_s = '3'
            else:
                new_s = '1'
            new_term = new_s + term[1] + ',' + term[2:]

            fixes.append([conf, old_term, energy, conf, new_term])
        else:
            new_term = term[0:2] + ',' + term[2:]

        level[2] = new_term
        i += 1

    # Check for duplicates
    confs_terms = [[level[1],level[2]] for level in conf_res]
    for ilvl in range(1, len(confs_terms)):
        if confs_terms[ilvl] in confs_terms[:ilvl]:
            existing_ilvl = confs_terms[:ilvl].index(confs_terms[ilvl])

            # Check duplicates of term
            old_term = conf_res[ilvl][2]
            s = conf_res[ilvl][5]
            existing_s = conf_res[existing_ilvl][5]
            if s > existing_s:
                new_term_existing = str(round(2*math.floor(float(existing_s))+1)) + conf_res[existing_ilvl][2][1:]
                new_term_current = str(round(2*math.ceil(float(s))+1)) + conf_res[ilvl][2][1:]
                conf_res[existing_ilvl][2] = new_term_existing
                conf_res[ilvl][2] = new_term_current
            else:
                new_term_existing = str(round(2*math.ceil(float(existing_s))+1)) + conf_res[existing_ilvl][2][1:]
                new_term_current = str(round(2*math.floor(float(s))+1)) + conf_res[ilvl][2][1:]
                conf_res[existing_ilvl][2] = new_term_existing
                conf_res[ilvl][2] = new_term_current

            # Print warning about duplicates to alert user to check results
            print(f'WARNING: DUPLICATE {confs_terms[ilvl][0]} {confs_terms[ilvl][1].replace(",","")} found at levels {existing_ilvl+1} and {ilvl+1}')
            print(f'  Level {existing_ilvl+1}: S={existing_s:.2f} -> {new_term_existing.replace(",","")}')
            print(f'  Level {ilvl+1}: S={s:.2f} -> {new_term_current.replace(",","")}')
            print(f'  Please verify these levels in the original calculation.')
                
            # If term remains same, check duplicates of configuration
            if conf_res[existing_ilvl][2] == confs_terms[ilvl][1]:
                if conf_res[existing_ilvl][11] and conf_res[existing_ilvl][1] == confs_terms[ilvl][0]:
                    # Check secondary config already in list
                    if [conf_res[existing_ilvl][11], conf_res[existing_ilvl][2]] not in confs_terms[:ilvl]:
                        main_conf = conf_res[existing_ilvl][1]
                        new_conf = conf_res[existing_ilvl][11]
                        conf_res[existing_ilvl][1] = new_conf
                        conf_res[existing_ilvl][11] = main_conf
                        # Generate fix for E1.RES: change primary to secondary config
                        # Format: [old_conf, old_term, energy, new_conf, new_term]
                        # Use original term (from confs_terms_old) for matching E1.RES
                        original_term = confs_terms_old[existing_ilvl][1] if existing_ilvl < len(confs_terms_old) else conf_res[existing_ilvl][2]
                        fixes.append([main_conf, original_term, conf_res[existing_ilvl][3], new_conf, original_term])

            # Check if there's currently a list of fixes, update the one for existing_ilvl
            if fixes:
                # Match by energy to find the correct fix for the existing level
                existing_energy = conf_res[existing_ilvl][3]
                for fix in fixes:
                    if fix[2] == existing_energy:  # Match by energy
                        fix[4] = conf_res[existing_ilvl][2]  # Update new_term

    # Fix skipped configuration levels (e.g., 5g, 7g -> 5g, 6g when 6g is secondary)
    # Pass original terms so fixes can match against matrix element file (which has original terms)
    conf_res, skipped_fixes = fix_skipped_config_levels(conf_res, confs_terms_old)

    # Merge term fixes and skipped fixes by (old_conf, old_term, energy)
    # If both exist for the same original state, combine config change and term change
    # Fix format: [old_conf, old_term, energy, new_conf, new_term]
    merged_fixes = {}
    for fix in fixes:  # term fixes: [old_conf, old_term, energy, conf, new_term]
        key = (fix[0], fix[1], fix[2])  # (old_conf, old_term, energy)
        merged_fixes[key] = [fix[0], fix[1], fix[2], fix[3], fix[4]]

    for fix in skipped_fixes:  # config fixes: [old_conf, old_term, energy, new_conf, new_term]
        key = (fix[0], fix[1], fix[2])  # (old_conf, old_term, energy)
        if key in merged_fixes:
            # Combine: keep new_term from term fix, use new_conf from skipped fix
            merged_fixes[key][3] = fix[3]  # update new_conf
        else:
            merged_fixes[key] = [fix[0], fix[1], fix[2], fix[3], fix[4]]

    fixes = list(merged_fixes.values())

    # Print merged fixes
    if fixes:
        print('FIXES for', filename + ':')
        for fix in fixes:
            old_conf, old_term, energy, new_conf, new_term = fix
            old_term_print = old_term.replace(',', '')
            new_term_print = new_term.replace(',', '')
            print(f'  {old_conf} {old_term_print} -> {new_conf} {new_term_print}')

    return conf_res, fixes

def parse_matrix_res(filename):
    try:
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
    except:
        print(filename + ' not found')
        return []
    
    matrix_res = []
    for line in lines[1:]:
        # E1.RES format: < conf2 || E1 || conf1 > E1_L  E1_V  E2  E1  E2-E1  WL  Tr. Rate
        # M2.RES format: < conf2 || M2 || conf1 > M2   E2  E1  E2-E1  WL
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
        
        if Tk == 'E1':
            energy1 = re.findall(r"\d+\.\d+", line)[2:4][1]
            energy2 = re.findall(r"\d+\.\d+", line)[2:4][0]
            wavelength = float(re.findall(r"\d+\.\d+", line)[5])
        else:
            energy1 = re.findall(r"\d+\.\d+", line)[1:3][1]
            energy2 = re.findall(r"\d+\.\d+", line)[1:3][0]
            wavelength = float(re.findall(r"\d+\.\d+", line)[4])
        index1 = int(re.findall(r"\d+",line)[0])
        index2 = int(re.findall(r"\d+",line)[1])

        matrix_element_value = re.findall(r"\d+\.\d+", line)[:1]
        matrix_element_value = matrix_element_value[0] if matrix_element_value else None
        
        matrix_res.append([conf1, term1, J1, conf2, term2, J2, matrix_element_value, energy1, energy2, wavelength, index1, index2])
        
    return matrix_res

def fix_matrix_res(fixes, res):
    """
    Apply fixes to matrix element results.
    fixes is a list of [old_conf, old_term, energy, new_conf, new_term].
    Uses energy to disambiguate when multiple fixes exist for the same (config, term).
    """
    # Group fixes by (config, term) for initial lookup
    # Then use energy to find the best match
    fix_groups = {}
    for fix in fixes:
        key = (fix[0], fix[1])  # (old_conf, old_term)
        if key not in fix_groups:
            fix_groups[key] = []
        fix_groups[key].append(fix)

    def find_best_fix(conf, term, energy):
        """Find the fix that best matches the given state by energy."""
        key = (conf, term)
        if key not in fix_groups:
            return None
        candidates = fix_groups[key]
        if len(candidates) == 1:
            return candidates[0]
        # Multiple candidates - find closest energy match
        # Energy in E1.RES is negative, in CONFFINAL is positive
        energy_abs = abs(float(energy))
        best_fix = None
        best_diff = float('inf')
        for fix in candidates:
            fix_energy = abs(float(fix[2]))
            diff = abs(energy_abs - fix_energy)
            if diff < best_diff:
                best_diff = diff
                best_fix = fix
        return best_fix

    for row in res:
        orig_conf1 = row[0]
        orig_term1 = row[1] + ',' + row[2]
        energy1 = row[7]
        orig_conf2 = row[3]
        orig_term2 = row[4] + ',' + row[5]
        energy2 = row[8]

        # Check if conf1 needs fix
        fix1 = find_best_fix(orig_conf1, orig_term1, energy1)
        if fix1:
            new_conf, new_term = fix1[3], fix1[4]
            term = new_term.split(',')[0]
            J = new_term.split(',')[1]
            row[0] = new_conf
            row[1] = term
            row[2] = J

        # Check if conf2 needs fix
        fix2 = find_best_fix(orig_conf2, orig_term2, energy2)
        if fix2:
            new_conf, new_term = fix2[3], fix2[4]
            term = new_term.split(',')[0]
            J = new_term.split(',')[1]
            row[3] = new_conf
            row[4] = term
            row[5] = J

    return res

def cmp_matrix_res(res1, res2, swaps, fixes):
    # fixes is a tuple (fixes1, fixes2) where:
    # - fixes1: fixes for CI+all-order: E1.RES (res1)
    # - fixes2: fixes for CI+MBPT: E1MBPT.RES (res2)
    fixes1, fixes2 = fixes

    matrix_res1 = parse_matrix_res(res1)
    matrix_res2 = parse_matrix_res(res2)
    second_order_exists = True
    if not matrix_res2:
        second_order_exists = False

    # Build a set of (config, term) that have fixes2 - these should not be swapped
    # Fixes take precedence over swaps to ensure NIST-matched labels are preserved
    fix2_keys = set()
    if fixes2:
        for fix in fixes2:
            fix2_keys.add((fix[0], fix[1]))  # (old_conf, old_term)

    # Check if a swap has to be made in res2 (second-order)
    # Skip swaps when there's a fix for the same (config, term) - fixes take precedence
    if swaps:
        for row in matrix_res2:
            for swap in swaps:
                # Check state 1
                term = row[1] + ',' + row[2]
                if row[0] == swap[1] and term == swap[2]:
                    # Skip if there's a fix for this (config, term)
                    if (row[0], term) in fix2_keys:
                        continue
                    row[0] = swap[0]

                # Check state 2
                term = row[4] + ',' + row[5]
                if row[3] == swap[1] and term == swap[2]:
                    # Skip if there's a fix for this (config, term)
                    if (row[3], term) in fix2_keys:
                        continue
                    row[3] = swap[0]

    # Apply fixes to the correct E1.RES files
    # fixes1 (from all-order CONFFINAL) -> E1.RES (matrix_res1)
    # fixes2 (from MBPT CONFFINAL) -> E1MBPT.RES (matrix_res2)
    if fixes1:
        matrix_res1 = fix_matrix_res(fixes1, matrix_res1)
    if second_order_exists and fixes2:
        matrix_res2 = fix_matrix_res(fixes2, matrix_res2)

    # Make a E1.RES array starting with results from res1
    matrix_res = []
    unmatched_matrix = []
    num_matches = 0

    for row1 in matrix_res1:
        matched = False
        conf11 = row1[0]
        term11 = row1[1] + row1[2]
        conf12 = row1[3]
        term12 = row1[4] + row1[5]
        me1 = row1[6]
        energy1 = row1[7]
        energy2 = row1[8]
        wavelength = row1[9]
        if second_order_exists:
            for row2 in matrix_res2:
                conf21 = row2[0]
                term21 = row2[1] + row2[2]
                conf22 = row2[3]
                term22 = row2[4] + row2[5]
                me2 = row2[6]
                if conf11 == conf21 and term11 == term21 and conf12 == conf22 and term12 == term22:
                    uncertainty = round(abs(float(me2) - float(me1)),5)
                    matrix_res.append([conf11, term11, conf12, term12, me1, uncertainty, energy1, energy2, wavelength])
                    matched = True
                    num_matches += 1
                    break
            if not matched:
                unmatched_matrix.append(row1)
                matrix_res.append([conf11, term11, conf12, term12, me1, '-', energy1, energy2, wavelength])
        else:
            matrix_res.append([conf11, term11, conf12, term12, me1, '-', energy1, energy2, wavelength])

    print(num_matches,'/',len(matrix_res), 'MATRIX ELEMENTS MATCHED')

    return matrix_res, unmatched_matrix

def cmp_res(res1, res2):
    conf_res1, fixes1 = parse_final_res(res1)
    second_order_exists = True
    try:
        conf_res2, fixes2 = parse_final_res(res2)
    except:
        second_order_exists = False
    
    confs1 = [level[1] for level in conf_res1]
    terms1 = [level[2] for level in conf_res1]
    energies_au1 = [level[3] for level in conf_res1]
    energies_cm1 = [level[4] for level in conf_res1]
    sec_confs1 = [level[11] for level in conf_res1]
    
    if second_order_exists:
        confs2 = [level[1] for level in conf_res2]
        terms2 = [level[2] for level in conf_res2]
        energies_au2 = [level[3] for level in conf_res2]
        energies_cm2 = [level[4] for level in conf_res2]
        sec_confs2 = [level[11] for level in conf_res2]
    
    # Make a CONF.RES array starting with results from res1
    conf_res = []
    for ilvl in range(len(confs1)):
        conf_res.append([confs1[ilvl], terms1[ilvl], energies_au1[ilvl], energies_cm1[ilvl]])
    
    # Loop through res2 results and find matches to res1 results
    if second_order_exists:
        not_matched1 = []
        not_matched2 = []
        matched_with_sec_confs = []
        for ilvl in range(len(confs2)):
            matched = False
            for ilvl2 in range(len(conf_res)):
                # Skip all-order states that have already been matched (have MBPT energies appended)
                if len(conf_res[ilvl2]) > 4:
                    continue
                if confs2[ilvl] == conf_res[ilvl2][0] and terms2[ilvl] == conf_res[ilvl2][1]:
                    conf_res[ilvl2].append(energies_au2[ilvl])
                    conf_res[ilvl2].append(energies_cm2[ilvl])
                    matched = True
                    continue
            # If primary configurations were not matched, check secondary configurations
            if not matched:
                # Check if MBPT primary matches all-order secondary
                for ilvl2 in range(len(sec_confs1)):
                    # Skip all-order states that have already been matched
                    if len(conf_res[ilvl2]) > 4:
                        continue
                    if confs2[ilvl] == sec_confs1[ilvl2] and terms2[ilvl] == conf_res[ilvl2][1]:
                        conf_res[ilvl2].append(energies_au2[ilvl])
                        conf_res[ilvl2].append(energies_cm2[ilvl])
                        matched_with_sec_confs.append([conf_res[ilvl2][0],confs2[ilvl],conf_res[ilvl2][1]]) # ind1 = res2, ind2 = res1
                        matched = True
                        break
            # If still not matched, check if MBPT secondary matches all-order primary
            if not matched:
                for ilvl2 in range(len(conf_res)):
                    # Skip all-order states that have already been matched
                    if len(conf_res[ilvl2]) > 4:
                        continue
                    if sec_confs2[ilvl] and sec_confs2[ilvl] == conf_res[ilvl2][0] and terms2[ilvl] == conf_res[ilvl2][1]:
                        conf_res[ilvl2].append(energies_au2[ilvl])
                        conf_res[ilvl2].append(energies_cm2[ilvl])
                        matched_with_sec_confs.append([sec_confs2[ilvl],conf_res[ilvl2][0],conf_res[ilvl2][1]]) # ind1 = res2, ind2 = res1
                        matched = True
                        break
            # Add res2 data that were not matched to array not_matched2
            if not matched:
                not_matched2.append([ilvl+1,confs2[ilvl], terms2[ilvl]])

        # Add res1 data that were not matched to array not_matched1
        for ilvl in range(len(confs1)):
            if len(conf_res[ilvl]) <= 4:
                not_matched1.append([ilvl+1,conf_res[ilvl][0],conf_res[ilvl][1]])
            
        # Print matched configurations from comparing secondary configurations
        for match in matched_with_sec_confs:
            print('SECONDARY MATCH FOUND: ', match)

        # Store unmatched configurations for external file output
        unmatched_energies = {
            'res1': res1,
            'res2': res2,
            'not_matched1': not_matched1,
            'not_matched2': not_matched2
        }
    
        # Assign uncertainty of 0 if no match was found
        for row in conf_res:
            if len(row) == 4:
                row.append(0)
                row.append(0)
        
        # Add uncertainties
        conf_res = add_uncertainties(conf_res)
        for i in range(len(conf_res1)):
            conf_res1[i].append(conf_res[i][-1])
        
        # Deduplicate fixes1 and fixes2 separately
        def deduplicate_fixes(fix_list):
            seen = set()
            result = []
            for fix in fix_list:
                fix_key = (fix[0], fix[1], fix[2], fix[3], fix[4])
                if fix_key not in seen:
                    seen.add(fix_key)
                    result.append(fix)
            return result

        fixes1 = deduplicate_fixes(fixes1)
        fixes2 = deduplicate_fixes(fixes2)

    else:
        for level in conf_res:
            level.append('-')
            level.append('-')
            level.append('-')
        for level in conf_res1:
            level.append('-')
        fixes2 = []
        matched_with_sec_confs = []
        unmatched_energies = {}

    return conf_res, conf_res1, matched_with_sec_confs, (fixes1, fixes2), unmatched_energies
    
def add_uncertainties(combined_conf_res):
    
    for row in combined_conf_res:
        if row[5] == 0:
            uncertainty = '-'
        else:
            uncertainty = abs(round(row[3]-row[5]))
        row.append(uncertainty)
    
    return combined_conf_res

def merge_res(res_even, res_odd, second_order_exists):
    ht_to_cm = 219474.63 # hartree to cm-1
    gs_parity = ''
    merged_res = []
    energy_cm1_shifted2 = '-'
    energy_cm2_shifted2 = '-'
    uncertainty = '-'
    if res_even[0][2] > res_odd[0][2]:
        gs_parity = 'even'
        merged_res = res_even
        for conf in res_odd:
            energy_cm2_shifted = round((res_even[0][2]-conf[2])*ht_to_cm, 2)
            if second_order_exists: 
                energy_cm2_shifted2 = round((res_even[0][4]-conf[4])*ht_to_cm, 2)
                uncertainty = round(abs(energy_cm2_shifted2-energy_cm2_shifted))
            if conf[4] == 0:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm2_shifted,conf[4],0,0])
            else:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm2_shifted,conf[4],energy_cm2_shifted2,uncertainty])
    else:
        gs_parity = 'odd'
        merged_res = res_odd
        for conf in res_even:
            energy_cm1_shifted = round((res_odd[0][2]-conf[2])*ht_to_cm, 2)
            if second_order_exists:
                energy_cm1_shifted2 = round((res_odd[0][4]-conf[4])*ht_to_cm, 2)
                uncertainty = round(abs(energy_cm1_shifted2-energy_cm1_shifted))
            if conf[4] == 0:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],0,0])
            else:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],energy_cm1_shifted2,uncertainty])
    
    return gs_parity, merged_res

def write_unmatched_file(unmatched_energies, unmatched_matrix, filename='unmatched.txt'):
    """
    Write unmatched energy levels and matrix elements to a file.

    Args:
        unmatched_energies: dict with 'odd' and 'even' keys, each containing
                           {'res1': filename, 'res2': filename, 'not_matched1': [...], 'not_matched2': [...]}
        unmatched_matrix: list of unmatched matrix element rows
        filename: output filename (default: 'unmatched.txt')

    Returns:
        True if file was written, False if no unmatched items
    """
    has_unmatched = False

    # Check for unmatched energies
    unmatched_energy_list = []
    for label, key in [('ODD PARITY', 'odd'), ('EVEN PARITY', 'even')]:
        unmatched = unmatched_energies.get(key, {})
        if unmatched and (unmatched.get('not_matched1') or unmatched.get('not_matched2')):
            has_unmatched = True
            unmatched_energy_list.append((label, unmatched))

    # Check for unmatched matrix elements
    if unmatched_matrix:
        has_unmatched = True

    # Write to file if there are unmatched items
    if has_unmatched:
        with open(filename, 'w') as f:
            # Write unmatched energy levels
            if unmatched_energy_list:
                f.write('='*60 + '\n')
                f.write('UNMATCHED ENERGY LEVELS\n')
                f.write('='*60 + '\n\n')
                for label, unmatched in unmatched_energy_list:
                    f.write(f'--- {label} ---\n')
                    if unmatched.get('not_matched1'):
                        f.write(f"From {unmatched['res1']}:\n")
                        for item in unmatched['not_matched1']:
                            f.write(f"  Level {item[0]}: {item[1]} {item[2]}\n")
                    if unmatched.get('not_matched2'):
                        f.write(f"From {unmatched['res2']}:\n")
                        for item in unmatched['not_matched2']:
                            f.write(f"  Level {item[0]}: {item[1]} {item[2]}\n")
                    f.write('\n')

            # Write unmatched matrix elements
            if unmatched_matrix:
                f.write('='*60 + '\n')
                f.write('UNMATCHED MATRIX ELEMENTS\n')
                f.write('='*60 + '\n\n')
                f.write(f"Total: {len(unmatched_matrix)} unmatched\n\n")
                for row in unmatched_matrix:
                    # row format: [conf1, term1, J1, conf2, term2, J2, me, energy1, energy2, wavelength, idx1, idx2]
                    conf1, term1, J1, conf2, term2, J2 = row[0], row[1], row[2], row[3], row[4], row[5]
                    f.write(f"  <{conf2} {term2}{J2} || E1 || {conf1} {term1}{J1}>\n")

        print(f"Unmatched items written to {filename}")
        return True

    return False

if __name__ == "__main__":
    conf_res_odd, full_res_odd, swaps_odd, fixes_odd, unmatched_odd = cmp_res('DATA_RAW/CONFFINALodd.RES','DATA_RAW/CONFFINALoddMBPT.RES')
    conf_res_even, full_res_even, swaps_even, fixes_even, unmatched_even = cmp_res('DATA_RAW/CONFFINALeven.RES','DATA_RAW/CONFFINALevenMBPT.RES')

    swaps = swaps_odd + swaps_even
    # fixes_odd and fixes_even are tuples (fixes1, fixes2)
    # Combine all-order fixes together and MBPT fixes together
    fixes1_odd, fixes2_odd = fixes_odd
    fixes1_even, fixes2_even = fixes_even
    fixes = (fixes1_odd + fixes1_even, fixes2_odd + fixes2_even)

    e1_res, unmatched_matrix = cmp_matrix_res('DATA_RAW/E1.RES','DATA_RAW/E1MBPT.RES',swaps,fixes)

    # Write unmatched items to file
    unmatched_energies = {'odd': unmatched_odd, 'even': unmatched_even}
    write_unmatched_file(unmatched_energies, unmatched_matrix)
    
