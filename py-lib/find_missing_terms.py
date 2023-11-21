import get_atomic_term

'''
This python script determines missing terms from portal csv
'''

def find_missing_terms(filename):
    # Import csv data
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Organize existing configurations and terms
    config_terms = {}
    for line in lines[1:]:
        config = line.split(',')[0]
        term = line.split(',')[1] + line.split(',')[2]
        if config not in config_terms:
            config_terms[config] = []
        config_terms[config] += [term]

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
    
    return missing_terms

if __name__ == '__main__':
    # TODO - add implementation with NIST/conf/portal
    # TODO - flag extra terms/terms that can't exist (from conf code)
    filename = input('Name of portal-csv file: ')
    
    missing_terms = find_missing_terms(filename)
    
    # Print out missing terms
    print('MISSING TERMS:')
    for config in missing_terms:
        print(config + ':', missing_terms[config])
    