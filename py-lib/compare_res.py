import re
import math
import os
import sys

'''
This script combines two CONFFINAL.RES files, e.g. ci+all-order and ci+mbpt results
'''

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
        index = int(re.findall("\d+",line)[0])
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
        
        cntrb1 = re.findall("\d+.\d+%",line)[0]
        try:
            cntrb2 = re.findall("\d+.\d+%",line)[1]
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
        s = level[5]
        old_term = term[0:2] + ',' + term[2]
        if even and term[0] == '2':
            if float(s) > 0.5:
                new_s = '3'
            else:
                new_s = '1'
            new_term = new_s + term[1] + ',' + term[2:]
            
            fixes.append([i,conf,old_term,conf,new_term])
        else:
            new_term = term[0:2] + ',' + term[2:]
            
        level[2] = new_term
        i += 1

    # Check for duplicates
    confs_terms = [[level[1],level[2]] for level in conf_res]
    for ilvl in range(1, len(confs_terms)):
        if confs_terms[ilvl] in confs_terms[:ilvl-1]:
            print('DUPLICATE FOUND:', confs_terms[ilvl], 'AT INDICES', confs_terms[:ilvl-1].index(confs_terms[ilvl])+1, 'AND', ilvl+1)
            existing_ilvl = confs_terms[:ilvl-1].index(confs_terms[ilvl])
            
            # Check duplicates of term
            old_term = conf_res[ilvl][2]
            s = conf_res[ilvl][5]
            existing_s = conf_res[existing_ilvl][5]
            if s > existing_s:
                conf_res[existing_ilvl][2] = str(round(2*math.floor(float(existing_s))+1)) + conf_res[existing_ilvl][2][1:]
                conf_res[ilvl][2] = str(round(2*math.ceil(float(s))+1)) + conf_res[ilvl][2][1:]
            else:
                conf_res[existing_ilvl][2] = str(round(2*math.ceil(float(existing_s))+1)) + conf_res[existing_ilvl][2][1:]
                conf_res[ilvl][2] = str(round(2*math.floor(float(s))+1)) + conf_res[ilvl][2][1:]
                
            # If term is remains same, check duplicates of configuration
            if conf_res[existing_ilvl][2] == confs_terms[ilvl][1]:
                if conf_res[existing_ilvl][1] == confs_terms[ilvl][0]:
                    # Check secondary config already in list
                    if [conf_res[existing_ilvl][11],conf_res[existing_ilvl][2]] not in confs_terms[:ilvl-1]:
                        main_conf = conf_res[existing_ilvl][1]
                        conf_res[existing_ilvl][1] = conf_res[existing_ilvl][11]
                        conf_res[existing_ilvl][11] = main_conf
                                    
            # Check if there's currently a list of fixes
            if fixes:
                # Check if duplicate is in list of fixes
                for fix in fixes:
                    # If index to fix is in list of fixes, update the configuration and term in fixes
                    if fix[0] == existing_ilvl+1:
                        fix[3] = conf_res[existing_ilvl][1]
                        fix[4] = conf_res[existing_ilvl][2]
                    if [ilvl+1,conf_res[ilvl][1],confs_terms_old[ilvl][1],conf_res[ilvl][1],conf_res[ilvl][2]] not in fixes:
                        fixes.append([existing_ilvl+1,confs_terms[ilvl][0],confs_terms[ilvl][1],conf_res[existing_ilvl][1],conf_res[existing_ilvl][2]])
                        fixes.append([ilvl+1,conf_res[ilvl][1],old_term,conf_res[ilvl][1],conf_res[ilvl][2]])
                
            print('     TERM OF LEVEL',existing_ilvl+1,'HAS BEEN UPDATED TO',conf_res[existing_ilvl][1],conf_res[existing_ilvl][2],'(WAS PREVIOUSLY',confs_terms[ilvl][0],confs_terms[ilvl][1] + ')')
            print('     TERM OF LEVEL',ilvl+1,'HAS BEEN UPDATED TO',conf_res[ilvl][1],conf_res[ilvl][2],'(WAS PREVIOUSLY',confs_terms[ilvl][0],confs_terms[ilvl][1] + ')')

    # Print fixes
    if fixes:
        print('FIXES for', filename + ':')
        for fix in fixes:
            print('#' + str(fix[0]) + ':',fix[1],fix[2],'->',fix[3],fix[4])

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
            energy1 = re.findall("\d+\.\d+", line)[2:4][1]
            energy2 = re.findall("\d+\.\d+", line)[2:4][0]
            wavelength = float(re.findall("\d+\.\d+", line)[5])
        else:
            energy1 = re.findall("\d+\.\d+", line)[1:3][1]
            energy2 = re.findall("\d+\.\d+", line)[1:3][0]
            wavelength = float(re.findall("\d+\.\d+", line)[4])
        index1 = int(re.findall("\d+",line)[0])
        index2 = int(re.findall("\d+",line)[1])

        matrix_element_value = re.findall("\d+\.\d+", line)[:1]
        matrix_element_value = matrix_element_value[0] if matrix_element_value else None
        
        matrix_res.append([conf1, term1, J1, conf2, term2, J2, matrix_element_value, energy1, energy2, wavelength, index1, index2])
        
    return matrix_res

def fix_matrix_res(fixes, res):
    for row in res:
        index1 = row[11]
        index2 = row[10]
        conf1 = row[3]
        term1 = row[4] + ',' + row[5]
        conf2 = row[0]
        term2 = row[1] + ',' + row[2]
        for fix in fixes:
            fix_index = fix[0]
            old_conf = fix[1]
            old_term = fix[2]
            new_conf = fix[3]
            new_term = fix[4]
            if index1 == fix_index and conf2 == old_conf and term2 == old_term:
                term = new_term.split(',')[0]
                J = new_term.split(',')[1]
                row[0] = new_conf
                row[1] = term
                row[2] = J
            if index2 == fix_index and conf1 == old_conf and term1 == old_term:
                term = fix[4].split(',')[0]
                J = fix[4].split(',')[1]
                row[3] = new_conf
                row[4] = term
                row[5] = J

    return res

def cmp_matrix_res(res1, res2, swaps, fixes):
    matrix_res1 = parse_matrix_res(res1)
    matrix_res2 = parse_matrix_res(res2)
    second_order_exists = True
    if not matrix_res2:
        second_order_exists = False
    
    # Check if a swap has to be made in res2 (second-order)
    if swaps:
        for row in matrix_res2:
            for swap in swaps:
                # Check state 1
                term = row[1] + ',' + row[2]
                if row[0] == swap[1] and term == swap[2]:
                    #print('SWAP FOUND:', row[0], '-->', swap[0])
                    row[0] = swap[0]
                                       
                # Check state 2
                term = row[4] + ',' + row[5]
                if row[3] == swap[1] and term == swap[2]:
                    #print('SWAP FOUND:', row[3], '-->', swap[0])
                    row[3] = swap[0]
    
    # Check if any levels have to be fixed
    if fixes:
        matrix_res1 = fix_matrix_res(fixes, matrix_res1)
        if second_order_exists: 
            matrix_res2 = fix_matrix_res(fixes, matrix_res2)

    # Make a E1.RES array starting with results from res1
    matrix_res = []
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
                print('NOT MATCHED:', row1)
                matrix_res.append([conf11, term11, conf12, term12, me1, '-', energy1, energy2, wavelength])
        else:
            matrix_res.append([conf11, term11, conf12, term12, me1, '-', energy1, energy2, wavelength])
    
    print(num_matches,'/',len(matrix_res), 'MATRIX ELEMENTS MATCHED')
    
    return matrix_res

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
                if confs2[ilvl] == conf_res[ilvl2][0] and terms2[ilvl] == conf_res[ilvl2][1]:
                    conf_res[ilvl2].append(energies_au2[ilvl])
                    conf_res[ilvl2].append(energies_cm2[ilvl])
                    matched = True
                    continue
            # If primary configurations were not matched, check secondary configurations
            if not matched:
                for ilvl2 in range(len(sec_confs1)):
                    if confs2[ilvl] == sec_confs1[ilvl2] and terms2[ilvl] == conf_res[ilvl2][1]:
                        conf_res[ilvl2].append(energies_au2[ilvl2])
                        conf_res[ilvl2].append(energies_cm2[ilvl2])
                        matched_with_sec_confs.append([conf_res[ilvl2][0],confs2[ilvl],conf_res[ilvl2][1]]) # ind1 = res2, ind2 = res1
                        matched = True
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

        # Print unmatched configurations from each CONFFINAL.RES file
        if not_matched1 or not_matched2:
            print('UNMATCHED CONFIGURATIONS FOUND:')
            print('from ' + res1 + ': ', not_matched1)
            print('from ' + res2 + ': ', not_matched2)
    
        # Assign uncertainty of 0 if no match was found
        for row in conf_res:
            if len(row) == 4:
                row.append(0)
                row.append(0)
        
        # Add uncertainties
        conf_res = add_uncertainties(conf_res)
        for i in range(len(conf_res1)):
            conf_res1[i].append(conf_res[i][-1])
        
        # Combine list of fixes
        fixes = fixes1 + fixes2
        
    else:
        for level in conf_res:
            level.append('-')
            level.append('-')
            level.append('-')
        for level in conf_res1:
            level.append('-')
        fixes = fixes1
        matched_with_sec_confs = []
        
    return conf_res, conf_res1, matched_with_sec_confs, fixes
    
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
                uncertainty = round(abs(energy_cm2_shifted2-energy_cm2_shifted))
            if conf[4] == 0:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],0,0])
            else:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],energy_cm1_shifted2,uncertainty])
    
    return gs_parity, merged_res

if __name__ == "__main__":
    conf_res_odd, full_res_odd, swaps_odd, fixes_odd = cmp_res('DATA_RAW/CONFFINALodd.RES','DATA_RAW/CONFFINALoddMBPT.RES')
    conf_res_even, full_res_even, swaps_even, fixes_even = cmp_res('DATA_RAW/CONFFINALeven.RES','DATA_RAW/CONFFINALevenMBPT.RES')
        
    swaps = swaps_odd + swaps_even
    fixes = fixes_odd + fixes_even

    e1_res = cmp_matrix_res('DATA_RAW/E1.RES','DATA_RAW/E1MBPT.RES',swaps,fixes)
    
