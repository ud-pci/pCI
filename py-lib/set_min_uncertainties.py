import os 
import sys
import yaml

def read_yaml(filename):
    """ 
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """ 

    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return config

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

def set_min_uncertainties(name, min_percentage):
    
    filename = name + '_Matrix_Elements_Theory.csv'
    try:
        f = open(filename,'r')
        lines = f.readlines()
        f.close()
    except FileNotFoundError as e:
        print('ERROR:', filename, 'NOT FOUND')
        sys.exit()
            
    matrix_res = []
    for line in lines:
        matrix_res.append(line.replace('\n','').split(','))
    
    cnt = 0
    for line in matrix_res[1:]:
        val = float(line[6]) # matrix element value
        try:
            unc = float(line[7]) # matrix element uncertainty
        except ValueError:
            continue
        if val*min_percentage/100 >= unc: cnt += 1
        line[7] = '{:,.5f}'.format(max(val*min_percentage/100, unc))
        
    new_filename = filename.split('.')[0] + '_adjusted.' + filename.split('.')[1]
    print('# of matrix element uncertainties adjusted:', cnt)
    f = open(new_filename, 'w')
    for line in matrix_res:
        f.write(','.join(line) + '\n')
    f.close()
    print(new_filename + ' has been written')

if __name__ == '__main__':
    
    # Read atom name from config.yml if it exists
    config_exists = os.path.isfile('config.yml')

    atom = ''
    if config_exists:
        config = read_yaml('config.yml')
        config_name = config['system']['name']
        if len(config_name.split()) == 1:
            atom = config_name + ' I'
        elif len(config_name.split()) == 2:
            if '+' in config_name.split()[1]: 
                atom = config_name
            elif config_name.split()[1].isnumeric():
                atom = config_name.split()[0] + ' ' + str(int(config_name.split()[1])-1) + '+'
            else:
                atom = config_name.split()[0] + ' ' + str(convert_roman_to_num(config_name.split()[1])-1) + '+'
        else:
            print('ERROR: atom name not supported')
            sys.exit()
    else:
        atom = input('Input name of atom: ')
    name = atom.replace(" ","_")
    
    min_percentages = {
        'Mg_I': 0.3,
        'Ca_I': 1.3,
        'Sr_I': 1.5
    }
    
    if name in min_percentages:
        min_percent = min_percentages[name]
    else:
        #min_percent = float(input('Input minimum uncertainty to set in percentage: '))
        min_percent = 1.5
    print('A minimum uncertainty of ' + str(min_percent) + '% has been set for all matrix elements.')
        
    set_min_uncertainties(name, min_percent)
    
    