import re

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
        return [], [], [], [], []

    energies_au = []
    energies_cm = []
    main_confs = []
    sec_confs = []
    terms_list = []
    ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
    Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    for line in lines[1:]:
        confs = [conf for conf in line.split('  ') if any(l in conf for l in ls)]
        confs = [conf.strip() for conf in confs]
        confs = [conf.replace(' ', '.') for conf in confs]
        
        terms = [term for term in line.split('  ') if any(L in term for L in Ls)]
        nums = [num for num in line.split('  ') if '.' in num]
      
        energies_au.append(float(nums[0]))
        energies_cm.append(float(nums[1]))
        main_confs.append(confs[0])
        terms_list.append(terms[0].replace(' ', ''))
        
        try:
            sec_confs.append(confs[1])
        except:
            sec_confs.append('')

    # Separate terms
    i = 0
    for term in terms_list:
        term = term[0:2] + ',' + term[2:]
        terms_list[i] = term
        i += 1

    return main_confs, terms_list, energies_au, energies_cm, sec_confs

def parse_matrix_res(filename):
    try:
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
    except:
        print(filename + ' not found')
        return [], [], [], []
    
    matrix_res = []
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
        
        matrix_res.append([conf1, term1, J1, conf2, term2, J2, matrix_element_value, energy1, energy2, wavelength])
        
    return matrix_res
        
def cmp_matrix_res(res1, res2, swaps):
    matrix_res1 = parse_matrix_res(res1)
    matrix_res2 = parse_matrix_res(res2)
    
    # Check if a swap has to be made in res2 (second-order)
    for row in matrix_res2:
        if swaps:
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
        
    # Make a E1.RES array starting with results from res1
    matrix_res = []
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
                continue
        if not matched:
            print('NOT MATCHED:', row1)
            matrix_res.append([conf11, term11, conf12, term12, me1, '-', energy1, energy2, wavelength])
    
    print(len(matrix_res), 'MATRIX ELEMENTS MATCHED')
            
    
    return matrix_res

def cmp_res(res1, res2):
    confs1, terms1, energies_au1, energies_cm1, sec_confs1 = parse_final_res(res1)
    confs2, terms2, energies_au2, energies_cm2, sec_confs2 = parse_final_res(res2)
    
    # Make a CONF.RES array starting with results from res1
    conf_res = []
    for ilvl in range(len(confs1)):
        conf_res.append([confs1[ilvl], terms1[ilvl], energies_au1[ilvl], energies_cm1[ilvl]])
    
    # Loop through res2 results and find matches to res1 results
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
    
    return conf_res, matched_with_sec_confs
    
def add_uncertainties(combined_conf_res):
    
    for row in combined_conf_res:
        if row[5] == 0:
            uncertainty = '-'
        else:
            uncertainty = abs(round(row[3]-row[5]))
        row.append(uncertainty)
    
    return combined_conf_res

def merge_res(res_even, res_odd):
    ht_to_cm = 219474.63 # hartree to cm-1
    gs_parity = ''
    merged_res = []
    if res_even[0][2] > res_odd[0][2]:
        gs_parity = 'even'
        merged_res = res_even
        for conf in res_odd:
            energy_cm2_shifted = round((res_even[0][2]-conf[2])*ht_to_cm, 2)
            energy_cm2_shifted2 = round((res_even[0][4]-conf[4])*ht_to_cm, 2)
            if conf[4] == 0:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm2_shifted,conf[4],0])
            else:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm2_shifted,conf[4],energy_cm2_shifted2])
    else:
        gs_parity = 'odd'
        merged_res = res_odd
        for conf in res_even:
            energy_cm1_shifted = round((res_odd[0][2]-conf[2])*ht_to_cm, 2)
            energy_cm1_shifted2 = round((res_odd[0][4]-conf[4])*ht_to_cm, 2)
            if conf[4] == 0:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],0])
            else:
                merged_res.append([conf[0],conf[1],conf[2],energy_cm1_shifted,conf[4],energy_cm1_shifted2])
    
    return gs_parity, merged_res

if __name__ == "__main__":
    conf_res_odd, swaps_odd = cmp_res('ci+all-order/odd/CONFFINAL.RES','ci+second-order/odd/CONFFINAL.RES')
    conf_res_even, swaps_even = cmp_res('ci+all-order/even/CONFFINAL.RES','ci+second-order/even/CONFFINAL.RES')
    
    swaps = swaps_odd + swaps_even
    e1_res = cmp_matrix_res('ci+all-order/dtm/E1.RES','ci+second-order/dtm/E1.RES',swaps)
