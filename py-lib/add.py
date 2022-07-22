import yaml
import re
import sys

def read_config(filename):
    """Read yaml input file""" 
    system = []
    # set default parameters
    num_dvdsn_iterations = 50
    num_energy_levels = 5
    system = system + [{'num_dvdsn_iterations':num_dvdsn_iterations}, {'num_energy_levels':num_energy_levels}]
    
    # read yml file with inputs
    with open(filename,'r') as f:
        config = yaml.safe_load(f)
        try:
            system = system + config['system']
            configurations = config['configurations']
            orbitals = config['orbitals']
            excitations = config['excitations']
            basis = config['basis']
        except KeyError as e:
            print('ERROR: ' + str(e) + ' is missing')
        
    return system, configurations, basis, orbitals, excitations

def count_valence(configurations):
    """ Count the number of valence electrons in each basic configuration """
    configurations_even = configurations['even']
    configurations_odd = configurations['odd']
    count = 0
    num_val = sum(int(re.findall('[0-9]+', orbital)[-1]) for orbital in configurations_even[0].split())

    configurations = configurations_even + configurations_odd

    # Check that each configuration has the same number of valence electrons
    for configuration in configurations:
        count = sum(int(re.findall('[0-9]+', orbital)[-1]) for orbital in configuration.split())
        if count != num_val:
            print('ERROR: ' + configuration + ' has ' + str(count) + ' electrons, but ' 
                    + configurations[0] + ' has ' + str(num_val))
            sys.exit()

    return num_val

def count_excitations(excitations):
    """Return the multiplicity of excitations given list of excitations"""
    count = 0
    for excitation in excitations:
        if list(excitation.values()) == [True]:
            count = count + 1
    return count

def expand_orbitals(basis, orbital_list):
    """Return a table of orbitals and occupation numbers"""

    orb_occ = {}
    nmin = 1
    nmax = 0

    core_orbs = 0
    active_orbs = orbital_list['active']

    # Add orbitals from basis to table
    orb_occ, nmax = expand_basis(basis)

    # Add orbitals from active orbital list to table
    for orbital in active_orbs:
        min_occ = list(orbital.values())[0].split()[0]
        max_occ = list(orbital.values())[0].split()[1]
        for key in orbital:
            l = re.findall('[spdfgh]+', key)
            n1 = int(re.findall('[0-9]+', key)[0])
            try:
                n2 = int(re.findall('[0-9]+', key)[1])
            except IndexError:
                n2 = n1
            if n1 < nmin:
                nmin = n1
            if n2 > nmax:
                nmax = n2
            for n in range(n1, n2 + 1):
                orb_occ[str(n) + l[0]] = min_occ, max_occ

    # Strip core orbitals from table
    try:
        core_orbs = re.findall('[1-9][spdfgh]+', orbital_list['core'])
        for orbital in core_orbs:
            if orbital in orb_occ:
                orb_occ.pop(orbital)
    except KeyError as e:
        print('No orbitals in core')
    except TypeError as e:
        print('No orbitals in core')

    orb_occ['num_orb'] = len(orb_occ)
    orb_occ['core'] = core_orbs
    orb_occ['nmin'] = nmin
    orb_occ['nmax'] = nmax

    return orb_occ

def expand_basis(basis):
    orb_occ = {}
    min_occ = 0
    max_occ = 2

    nmin_dict = {'s': 1, 'p': 2, 'd': 3, 'f': 4, 'g': 5, 'h':6}
    n_array = re.findall('[0-9]+', basis)
    l_array = re.findall('[spdfghi]+', basis)

    for i in range(len(n_array)):
        for l in l_array[i]:
            if l in nmin_dict:
                for n in range(nmin_dict[l], int(n_array[i]) + 1):
                    orb_occ[str(n) + l] = str(min_occ), str(max_occ)

    return orb_occ, int(max(n_array))
    
def write_add_inp(filename, system, configurations, orbitals, multiplicity, num_val, orb_occ, parity):
    """Write ADD.INP file for specified parity"""
    filename = filename[0:-4] + parity + filename[-4:]
    f = open(filename,'w')

    # Write header
    f.write('Ncor= {:d}\n'.format(len(configurations[parity])))
    f.write('NsvNR {:d}\n'.format(orb_occ['num_orb']))
    f.write('mult= {:d}\n'.format(multiplicity))
    f.write(' NE = {:d}\n\n'.format(num_val))

    # Reformat configurations
    configurations_formatted = format_configurations(configurations[parity])
    max_orb_per_line = 6

    # Write list of basic configurations
    for configuration in configurations_formatted:
        num_lines = len(configuration)//max_orb_per_line + 1
        for line in range(num_lines):
            f.write('L:  ' + ' '.join(configuration[line*max_orb_per_line:(line+1)*max_orb_per_line]) + '\n')
    f.write('\n')

    # Write list of orbitals and occupation numbers
    nmin = orb_occ['nmin']
    nmax = orb_occ['nmax']
    count = 0

    for i in range(nmin, nmax+1):
        for l in 'spdfgh':
            orb = str(i)+l
            if orb in orb_occ:
                orb_occ_formatted = format_orb_occ(orb, orb_occ[orb])
                f.write('  ' + orb_occ_formatted)
                count += 1
                if count == 6:
                    f.write('\n')
                    count = 0

    f.write('\n')

    num_core_orb = 0
    try:
        core_formatted = convert_to_digits(orb_occ['core'])
        num_core_orb = len(core_formatted)
    except:
        pass
        
    # Write head of CONF.INP
    f.write('>>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>\n')
    f.write('  ' + system['name'] + ' ' + parity + '\n')
    f.write('  Z = ' + str(system['atomic_number']) + '\n')
    f.write(' Am = ' + str(system['atomic_mass']) + '\n')
    f.write('  J =  0.0 \n')
    f.write(' Jm =  0.0 \n')
    f.write(' Nso=  ' + str(num_core_orb) + '\n')
    f.write(' Nc =   10 \n')
    f.write(' Kv =  4 \n')
    f.write(' Nlv=  ' + str(system['num_energy_levels']) + '\n')
    f.write(' Ne =  ' + str(num_val) + '\n')
    f.write(' Kl4=  1 \n')
    f.write(' Nc4= 20 \n')
    f.write(' Gj = 0.0000 \n')
    f.write('Crt4= 0.0001 \n')
    f.write('kout= 0 \n')
    f.write('Ncpt= 0 \n')
    f.write('Cut0= 0.0001 \n')
    f.write('N_it= ' + str(system['num_dvdsn_iterations']) + '\n')
    kbrt = sum((system['include_gaunt'], system['include_breit']))
    f.write('Kbrt= ' + str(kbrt) + '\n')
    f.write('Gnuc= 1.07 \n')
    
    # Write core orbitals
    if num_core_orb != 0:
        count = 0
        for i in range(0, num_core_orb):
            orb = core_formatted[i]
            f.write('    ' + orb)
            count += 1
            if count == 6:
                f.write('\n')
                count = 0

    print(filename + ' has been written')

def format_configurations(configurations):
    """ Format configurations as 'nnlqq',
        where nn represents principal quantum number
               l represents angular quantum number
              qq represents occupation number """
    configurations_formatted = []
    for configuration in configurations:
        orbitals_formatted = []
        for orbital in configuration.split():
            num_digits = len(re.findall('[0-9]+', orbital)[0])
            if num_digits == 1:
                orbital_formatted = orbital.center(5, ' ')
            elif num_digits == 2:
                orbital_formatted = orbital.ljust(5, ' ')
            else:
                print('ERROR: orbital ' + orbital + ' exceeds n=99')     
                sys.exit()         
            orbitals_formatted.append(orbital_formatted)
        configurations_formatted.append(orbitals_formatted)
    return configurations_formatted

def format_orb_occ(orb, occ):
    """ Format configurations as 'nnl qq qq',
        where nn represents principal quantum number
               l represents angular quantum number
              qq represents occupation number """
    orb_occ_formatted = orb.rjust(3, ' ') + occ[0].rjust(3, ' ') + occ[1].rjust(3, ' ')

    return orb_occ_formatted

def parse_system(system):
    """ Create dictionary with system parameters """ 
    params = {}
    for param in system:
        params[list(param.keys())[0]] = param[list(param.keys())[0]]

    return params

def convert_to_digits(orbitals):
    """ Converts a list of non-relativistic subshells to digital format """
    occ_dict = {'s': 2, 'p': (2, 4), 'd': (4, 6), 'f': (6, 8), 'g': (8, 10)}
    l_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}
    digital_orbitals = []

    if type(orbitals) == str:
        orbitals = [orbitals]

    for orbital in orbitals:
        nn = re.findall('[0-9]+', orbital.split()[0])[0].rjust(2, '0')
        l = re.findall('[spdfg]+', orbital.split()[0])[0]
        if l == 's':
            qq = str(occ_dict[l]).rjust(2, '0')
            nnlqq = nn[:1] + '.' + nn[1:] + str(l_dict[l]) + qq
        else:
            qq = str(occ_dict[l][0]).rjust(2, '0'), str(occ_dict[l][1]).rjust(2, '0')
            nnlqq = ['-' + nn[:1] + '.' + nn[1:] + str(l_dict[l]) + qq[0], nn[:1] + '.' + nn[1:] + str(l_dict[l]) + qq[1]]
        digital_orbitals.append(nnlqq)

    digital_orbitals = flatten_list(digital_orbitals)

    return digital_orbitals

def flatten_list(nested_list):
    """ Flatten a list of lists into a single list """
    flat_list = []
    for sublist in nested_list:
        if isinstance(sublist, str):
            flat_list.append(sublist)
        else:
            for item in sublist:
                flat_list.append(item)

    return flat_list

if __name__ == "__main__":
    system, configurations, basis, orbitals, excitations = read_config('add_config.yml')

    params = parse_system(system)
    num_val = count_valence(configurations)
    orb_occ = expand_orbitals(basis, orbitals)
    multiplicity = count_excitations(excitations)

    write_add_inp('ADD.INP', params, configurations, orbitals, multiplicity, num_val, orb_occ, 'even')
    write_add_inp('ADD.INP', params, configurations, orbitals, multiplicity, num_val, orb_occ, 'odd')


    