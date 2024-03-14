import pandas as pd

def format_time_in_seconds(time):
    time_in_seconds = 0
    unit = ''
    if time < 10e-6:
        time_in_seconds = time*1e9
        unit = 'ns'
    elif time >= 10e-6 and time < 10e-3:
        time_in_seconds = time*1e6
        unit = 'us'
    else:
        time_in_seconds = time
        unit = 's'
        
    return f'{time_in_seconds:.3f}' + ' ' + unit
        

def calc_lifetimes(tr_file):
    
    # Read transition rates from file tr_file
    f = open(tr_file, 'r')
    lines = f.readlines()
    f.close()
    
    atom = tr_file.split('_')[0] + '_' + tr_file.split('_')[1]
    
    # Open files to write lifetimes and branching ratios
    filename_lifetimes = atom + '_Lifetimes_Error_Check.csv'    
    filename_br_ratios = atom + '_Transition_Rates_Error_Check.csv'
    
    lifetimes_df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J', 'lifetime_display'])
    br_ratios_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                                         'state_two_configuration', 'state_two_term', 'state_two_J',
                                         'wavelength_display', 'matrix_element_display', 
                                         'branching_ratio_display', 'transition_rate_display'])
    
    # Parse transition rates file for energy levels and sum up transition rates
    tr_rates = {}
    for line in lines[1:]:
        energy1 = float(line.split(',')[8])
        energy2 = float(line.split(',')[9])
        configuration1 = line.split(',')[0]
        term1 = line.split(',')[1]
        J1 = line.split(',')[2]
        termJ1 = term1 + J1
        configuration2 = line.split(',')[3]
        term2 = line.split(',')[4]
        J2 = line.split(',')[5]
        termJ2 = term2 + J2
        matrix_element = line.split(',')[6]
        wavelength = float(line.split(',')[10])
        tr_rate = float(line.split(',')[11][:-1])
        config1 = configuration1 + ' ' + termJ1
        config2 = configuration2 + ' ' + termJ2
        if (energy2 < energy1): 
            try:
                tr_rates[config1].append([config2, tr_rate, matrix_element, wavelength])
            except KeyError:
                tr_rates[config1] = []
                tr_rates[config1].append([config2, tr_rate, matrix_element, wavelength])
        else:
            try:
                tr_rates[config2].append([config1, tr_rate, matrix_element, wavelength])
            except KeyError:
                tr_rates[config2] = []
                tr_rates[config2].append([config1, tr_rate, matrix_element, wavelength])
    
    # Calculate lifetimes and branching ratios
    lifetimes = []
    f = open('transitions.txt','w')
    for config, rates in tr_rates.items():
        configuration = config.split(' ')[0]
        term = config.split(' ')[1][0:2]
        J = config.split(' ')[1][-1]
        
        total_rates = 0
        for rate in rates:
            total_rates += rate[1]
        
        lifetime = 1/total_rates*1e9
        lifetimes.append([configuration, term, J, round(lifetime,3)])
        
        f.write(config+' ->\n')
        for rate in rates:
            f.write('      ' + rate[0] + ' ' + str(rate[1]) + ' ' + rate[2] + ' ' +  str(rate[3]) + '\n')
            configuration2 = rate[0].split(' ')[0]
            term2 = rate[0].split(' ')[1][0:2]
            J2 = rate[0].split(' ')[1][-1]
            matrix_element = rate[2]
            wavelength = rate[3]
            tr_rate = rate[1]
            branching_ratio = tr_rate/total_rates
            row = {'state_one_configuration': configuration, 'state_one_term': term, 'state_one_J': J,
                   'state_two_configuration': configuration2, 'state_two_term': term2, 'state_two_J': J2,
                   'wavelength_display': f"{wavelength:.2f}", 'matrix_element_display': matrix_element, 
                   'branching_ratio_display': f"{branching_ratio:.3e}", 'transition_rate_display': f"{tr_rate:.3e}"}
            br_ratios_df.loc[len(br_ratios_df.index)] = row
    f.close()

    # Write lifetimes to file
    for lifetime in lifetimes:
        row = {'state_configuration': lifetime[0], 'state_term': lifetime[1], 'state_J': lifetime[2], 'lifetime_display': lifetime[3]}
        lifetimes_df.loc[len(lifetimes_df.index)] = row
        
    br_ratios_df.to_csv(filename_br_ratios, index=False)
    lifetimes_df.to_csv(filename_lifetimes, index=False)

if __name__ == "__main__":
    filename = input("Input csv file: ")
    calc_lifetimes(filename)