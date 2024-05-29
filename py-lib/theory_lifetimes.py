import pandas as pd
from compare_res import parse_matrix_res
import datetime

def format_time_in_seconds(time):
    str_time = ''
    time_cut = 0
    unit = ''
    
    # if time less than 1 us
    if time < 10e-6:
        time_cut = time*1e9
        unit = 'ns'
    # if time between 1 us and 1 ms
    elif time >= 10e-6 and time < 10e-3:
        time_cut = time*1e6
        unit = 'us'
    # if time between 3 min and 3 hrs
    elif time >= 3*60 and time < 3*3600:
        time_cut = time/60
        unit = 'min'
    # if time between 3 hrs and 3 days
    elif time >= 3*3600 and time < 3*24*3600:
        time_cut = time/3600
        unit = 'hr'
    # if time between 3 days and 3 years
    elif time >= 3*24*3600 and time < 3*365*24*3600:
        time_cut = time/(24*3600)
        unit = 'days'
    # else time is in seconds
    else:
        time_cut = time
        unit = 's'
    
    str_time = f'{time_cut:.0f}' + ' ' + unit
        
    return str_time
        

def parse_dtm_res():
    '''
    
    '''
    # Open files to write lifetimes and branching ratios
    filename_lifetimes = 'lifetimes.csv'    
    filename_br_ratios = 'transition_rates.csv'
    
    lifetimes_df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J', 'lifetime'])
    br_ratios_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                                         'state_two_configuration', 'state_two_term', 'state_two_J', 'type',
                                         'wavelength (nm)', 'reduced_matrix_element', 
                                         'branching_ratio', 'transition_rate'])
    # Parse the RES files
    e1_res = parse_matrix_res('E1.RES')
    e2_res = parse_matrix_res('E2.RES')
    e3_res = parse_matrix_res('E3.RES')
    
    m1_res = parse_matrix_res('M1.RES')
    m2_res = parse_matrix_res('M2.RES')
    m3_res = parse_matrix_res('M3.RES')
    
    # Generate dictionary of matrix elements
    matrix_elements = {}
    if e1_res: add_res_to_dict('E1', e1_res, matrix_elements)
    if e2_res: add_res_to_dict('E2', e2_res, matrix_elements)
    if e3_res: add_res_to_dict('E3', e3_res, matrix_elements)
    
    if m1_res: add_res_to_dict('M1', m1_res, matrix_elements)
    if m2_res: add_res_to_dict('M2', m2_res, matrix_elements)
    if m3_res: add_res_to_dict('M3', m3_res, matrix_elements)
    
    
    # Calculate lifetimes and branching ratios
    ignore_g = True
    lifetime_threshold = 1e10 # set lifetime threshold to 317.1 yrs
    lifetimes = []
    f = open('transitions.txt','w')
    
    # Write format
    f.write('upper state (lifetime) -> \n')
    f.write(' '*6 + 'Tk lower_state transition_rate reduced_matrix_element_value wavelength branching_ratio \n\n')
    
    for config, rates in matrix_elements.items():
        configuration = config.split(' ')[0]
        term = config.split(' ')[1][0:2]
        J = config.split(' ')[1][-1]
        
        total_rates = 0
        for rate in rates:
            total_rates += rate[2]
        
        lifetime = 1/total_rates*1e9
        lifetimes.append([configuration, term, J, round(lifetime,3)])
        
        # skip to next configuration if conditions are met
        if lifetime >= lifetime_threshold: 
            print(configuration + ' ' + term + J + ' lifetime exceeded threshold of 317.1 yrs, so it is skipped')
            continue
        if ignore_g and 'G' in term:
            print(configuration + ' ' + term + J + ' skipped due to presence of G term')
            continue
        
        f.write(config + ' (lifetime = ' + format_time_in_seconds(lifetime) + ') ->\n')
        for rate in rates:
            matrix_element_type = rate[0]
            state2 = rate[1]
            configuration2 = state2.split(' ')[0]
            term2 = state2.split(' ')[1][0:2]
            J2 = state2.split(' ')[1][-1]
            matrix_element = rate[3]
            wavelength = rate[4]
            tr_rate = rate[2]
            branching_ratio = tr_rate/total_rates
            
            # create row for csv
            row = {'state_one_configuration': configuration, 'state_one_term': term, 'state_one_J': J,
                   'state_two_configuration': configuration2, 'state_two_term': term2, 'state_two_J': J2, 'type': matrix_element_type,
                   'wavelength (nm)': f"{wavelength:.2f}", 'reduced_matrix_element': matrix_element, 
                   'branching_ratio': f"{branching_ratio:.3e}", 'transition_rate': f"{tr_rate:.3e}"}
            br_ratios_df.loc[len(br_ratios_df.index)] = row
            
            # write to transitions file
            f.write(' '*6 + matrix_element_type + ' '*2 + state2 + ' '*4 + 
                    f"{tr_rate:.3e}" + ' '*4 + f"{matrix_element:.6e}" + ' '*4 + 
                    f"{wavelength:.2e}" + ' '*4 + f"{branching_ratio:.3e}" + '\n')
    f.close()
    
    # Write lifetimes to file
    for lifetime in lifetimes:
        if lifetime[3] >= lifetime_threshold: 
            continue
        if ignore_g and 'G' in lifetime[0]:
            continue
        row = {'state_configuration': lifetime[0], 'state_term': lifetime[1], 'state_J': lifetime[2], 'lifetime': f"{lifetime[3]:.3e}"}
        lifetimes_df.loc[len(lifetimes_df.index)] = row
        
    br_ratios_df.to_csv(filename_br_ratios, index=False)
    lifetimes_df.to_csv(filename_lifetimes, index=False)

def add_res_to_dict(matrix_element_type, res, dict):
    
    # [conf1, term1, J1, conf2, term2, J2, matrix_element_value, energy1, energy2, wavelength, index1, index2])
    for line in res:
        # configuration 1
        conf1 = line[0]
        term1 = line[1]
        J1 = line[2]
        termJ1 = term1 + J1
        energy1 = float(line[7])
        config1 = conf1 + ' ' + termJ1
        
        # configuration 2
        conf2 = line[3]
        term2 = line[4]
        J2 = line[5]
        termJ2 = term2 + J2
        energy2 = float(line[8])
        config2 = conf2 + ' ' + termJ2
        
        # matrix elementa
        matrix_element_val = float(line[6])
        
        # recalculated wavelength for more precision
        ht_to_cm = 219474.63 # hartree to cm-1
        wavelength_nm = 1e7/(abs(energy2-energy1)*ht_to_cm) # in nm
        wavelength = wavelength_nm*10 # wavelength in angstroms
        
        # calculate transition rate from formulae (https://www1.udel.edu/atom/about.html)
        J = float(J1) if energy2 > energy1 else float(J2)
        
        line_strength = matrix_element_val**2

        if matrix_element_type == 'E1':
            tr_rate = ((2.02613*10**18)/((2*J+1)*wavelength**3))*line_strength
        elif matrix_element_type == 'M1':
            tr_rate = ((2.69735*10**13)/((2*J+1)*wavelength**3))*line_strength
        elif matrix_element_type == 'E2':
            tr_rate = ((1.11995*10**18)/((2*J+1)*wavelength**5))*line_strength
        elif matrix_element_type == 'M2':
            tr_rate = ((1.49097*10**13)/((2*J+1)*wavelength**5))*line_strength
        elif matrix_element_type == 'E3':
            tr_rate = ((3.14441*10**17)/((2*J+1)*wavelength**7))*line_strength
        elif matrix_element_type == 'M3':
            tr_rate = ((4.18610*10**12)/((2*J+1)*wavelength**7))*line_strength
        else:
            tr_rate = 0

        # add to dictionary
        if (energy2 < energy1): 
            try:
                dict[config1].append([matrix_element_type, config2, tr_rate, matrix_element_val, wavelength_nm])
            except KeyError:
                dict[config1] = []
                dict[config1].append([matrix_element_type, config2, tr_rate, matrix_element_val, wavelength_nm])
        else:
            try:
                dict[config2].append([matrix_element_type, config1, tr_rate, matrix_element_val, wavelength_nm])
            except KeyError:
                dict[config2] = []
                dict[config2].append([matrix_element_type, config1, tr_rate, matrix_element_val, wavelength_nm])
    
    

if __name__ == "__main__":
    parse_dtm_res()