import get_atomic_term
import sys

'''
This python script determines missing terms from NIST database/conf calculations/portal database
'''

def find_missing_terms(filename):
    '''
    This function finds missing terms from the specified input file
    NIST database and portal files should be csv-formatted
    conf calculation file can be csv-formatted or RES-formatted
    '''
    # Import csv data
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Organize existing configurations and terms
    config_terms = {}
    duplicate_terms = {}
    if filename[-3:] == 'csv':
        for line in lines[1:]:
            config = line.split(',')[0]
            term = line.split(',')[1] + line.split(',')[2]
            if config not in config_terms:
                config_terms[config] = []
            else:
                if term in config_terms[config]:
                    if config not in duplicate_terms:
                        duplicate_terms[config] = []
                    duplicate_terms[config] += [term]
            config_terms[config] += [term]
    elif filename[-3:] == 'RES':
        ls = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        Ls = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
        for line in lines[1:]:
            config = [conf for conf in line.split('  ') if any(l in conf for l in ls)][0].replace(' ', '.')
            term = [term for term in line.split('  ') if any(L in term for L in Ls)][0].strip()
            if config not in config_terms:
                config_terms[config] = []
            else:
                if term in config_terms[config]:
                    if config not in duplicate_terms:
                        duplicate_terms[config] = []
                    duplicate_terms[config] += [term]
            config_terms[config] += [term]
    else:
        print(filename + ' is not compatible')
        sys.exit()

    # Obtain all possible terms for each configuration
    possible_config_terms = {}
    for config in config_terms:
        possible_config_terms[config] = get_atomic_term.scrape_term(config)

    # Determine missing terms for each configurations
    missing_terms = {}
    for config in config_terms:
        for term in possible_config_terms[config]:
            if term not in config_terms[config]:
                if config not in missing_terms:
                    missing_terms[config] = []
                missing_terms[config] += [term]
    
    # Find extra terms or terms that can't exist
    extra_terms = {}
    for config in config_terms:
        for term in config_terms[config]:
            if term not in possible_config_terms[config]:
                if config not in extra_terms:
                    extra_terms[config] = []
                extra_terms[config] += [term]
    
    return missing_terms, extra_terms, duplicate_terms

if __name__ == '__main__':
    filename = input('Name of input file: ')
    
    missing_terms, extra_terms, duplicate_terms = find_missing_terms(filename)
    
    # Print out missing terms
    print('MISSING TERMS:')
    for config in missing_terms:
        print(config + ':', missing_terms[config])
    
    # Print out extra terms
    print('EXTRA TERMS:')
    for config in extra_terms:
        print(config + ':', extra_terms[config])
        
    # Print out duplicate terms
    print('DUPLICATE TERMS:')
    for config in duplicate_terms:
        print(config + ':', duplicate_terms[config])   