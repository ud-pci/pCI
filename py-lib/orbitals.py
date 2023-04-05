import re
import itertools
from collections import OrderedDict
import copy
import sys

j_dict = {'s': ('1/2','1/2'), 'p': ('1/2', '3/2'), 'd': ('3/2', '5/2'), 'f': ('5/2', '7/2'), 'g': ('7/2', '9/2'), 'h': ('9/2', '11/2'), 'i': ('11/2','13/2'), 'k': ('13/2','15/2'), 'l': ('15/2, 17/2')}
l_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'k': 7, 'l': 8}
q_dict = {'s': (2), 'p': (2, 4), 'd': (4, 6), 'f': (6, 8), 'g': (8, 10), 'h': (10, 12), 'i': (12,14), 'k': (14,16), 'l': (16,18)}
maxq_dict = {'s': 2, 'p': 6, 'd': 10, 'f': 14, 'g': 18, 'h': 22, 'i': 26, 'k': 30, 'l': 34}

def read_bass(filename):
    """Read bass input file""" 
    orbitals = []
    with open(filename,'r') as f:
        for line in f:
            line_strip = re.sub(' +', ' ', line.strip()).split(" ")
            if is_integer(line_strip[0]): 
                nl, j = convert_digital_to_char(line_strip[1])
                orbitals.append([nl, j])
    return orbitals

def find_key(my_dict, value):
    for key, val in my_dict.items():
        if val == int(value): 
            return key
    return None

def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

def flatten_list(nested_list):
    """ Flattens a nested list into a single list """
    flat_list = []
    for sublist in nested_list:
        if isinstance(sublist, str):
            flat_list.append(sublist)
        else:
            for item in sublist:
                flat_list.append(item)

    return flat_list

def convert_digital_to_char(orbital):
    """ Converts an orbital from digital format to relativistic character format """

    if orbital[0] == '-':
        kj = 0
        n = int(orbital[1]+orbital[3])
        l = find_key(l_dict, orbital[4])
        j = j_dict[l][0]
    elif is_integer(orbital[0]) == True:
        kj = 1
        n = int(orbital[0]+orbital[2])
        l = find_key(l_dict, orbital[3])
        j = j_dict[l][1]

    return str(n)+l, j

def get_nn(orbital):
    """ Gets principal quantum number nn from an orbital string of format nnlqq """
    try:
        nn = re.findall('[0-9]+', orbital.split()[0])[0].rjust(2, '0')
    except IndexError:
        return None
    return nn

def get_l(orbital):
    """ Gets angular momentum quantum number l from an orbital string of format nnlqq """
    try:
        l = re.findall('[spdfghikl]+', orbital.split()[0])[0]
    except IndexError:
        return None
    return l

def get_qq(orbital):
    """ Gets occupation number qq from an orbital string of format nnlqq """
    try:
        qq = re.findall('[0-9]+', orbital.split()[0])[1]
        if int(qq) < 10:
            qq = '0' + str(qq)
    except IndexError:
        return q_dict[get_l(orbital)]
    return qq

def get_q2(l):
    """ Gets occupation number qq from an orbital string of format nnlqq """
    try:
        q2 = q_dict[l][0]
    except TypeError:
        q2 = q_dict[l]
    return q2

def convert_char_to_digital(orbitals):
    """ Converts a list of non-relativistic subshells to digital format """
    occ_dict = {'s': 2, 'p': (2, 4), 'd': (4, 6), 'f': (6, 8), 'g': (8, 10), 'h': (10, 12), 'i': (12,14), 'k': (14,16), 'l': (16,18)}
    digital_orbitals = []

    if type(orbitals) == str:
        orbitals = [orbitals]

    for orbital in orbitals:
        nn = get_nn(orbital)
        l = get_l(orbital)
    
        if l == 's':
            qq = str(occ_dict[l]).rjust(2, '0')
            qq = get_qq(orbital)
            nnlqq = nn[:1] + '.' + nn[1:] + str(l_dict[l]) + str(qq).rjust(2, '0')
        else:
            qq = str(occ_dict[l][0]).rjust(2, '0'), str(occ_dict[l][1]).rjust(2, '0')
            qq = get_qq(orbital)
            nnlqq = ['-' + nn[:1] + '.' + nn[1:] + str(l_dict[l]) + str(qq[0]).rjust(2, '0'), nn[:1] + '.' + nn[1:] + str(l_dict[l]) + str(qq[1]).rjust(2, '0')]
        digital_orbitals.append(nnlqq)

    digital_orbitals = flatten_list(digital_orbitals)

    return digital_orbitals

def assign_occ_to_config(config, occupation, sign):
    """ Updates the occupation number of a configuration """
    if occupation < 10:
        qq = '0' + str(occupation)
    else:
        qq = str(occupation)

    nn = get_nn(config)
    l = get_l(config)
    
    new_config = sign + nn[:1] + '.' + nn[1:] + str(l_dict[l]) + qq

    return new_config

def convert_nr_to_rel_config(nr_config):
    """ Expands a non-relativistic configuration into a list of relativistic configurations"""
    rel_configs = []
    nr_configs = nr_config.split()
    
    configs = []
    occupancies = []
    for config in nr_configs:
        occupancies.append(expand_occupancies(config))

    config = []
    nconfs = []
    nconf = 0
    for cindex in range(len(nr_configs)):
        l = get_l(nr_configs[cindex])
        if l == 's':
            config.append(assign_occ_to_config(nr_configs[cindex], occupancies[cindex][0][0], ' '))
            nconf += 1
        else:
            for occs in occupancies[cindex]:
                config1 = assign_occ_to_config(nr_configs[cindex], occs[0], '-')
                config2 = assign_occ_to_config(nr_configs[cindex], occs[1], ' ')
                config.append([config1, config2])
                nconf += 1
        configs.append(config)
        nconfs.append(nconf)
        nconf, config = 0, []
        
    rel_configs = list(itertools.product(*configs))

    return rel_configs

def expand_occupancies(nr_config):
    """ Expands a configuration into relativistic form """
    occupancies = []
    nn = re.findall('[0-9]+', nr_config.split()[0])[0].rjust(2, '0')
    l = re.findall('[spdfghikl]+', nr_config.split()[0])[0]
    qq = re.findall('[0-9]+', nr_config.split()[0])[1]

    if l == 's':
        return [[int(qq)]]
    
    diff = get_q2(l) - int(qq)
    if diff > 0:
        norb1 = int(qq)
        norb2 = 0
    else:
        norb1 = get_q2(l)
        norb2 = abs(diff)

    while type(q_dict[l]) == tuple: 
        occupancies.append([norb1, norb2])
        if norb1 == 0 or norb2 == q_dict[l][1]:
            break
        else:
            norb1 -= 1
            norb2 += 1
    
    return occupancies


def expand_nl(basis):
    """ Expands a basis of type string and returns a dictionary of maximum principle quantum number value for each partial wave key"""
    nlmax = {'s': 0, 'p': 0, 'd': 0, 'f': 0, 'g': 0, 'h': 0, 'i': 0, 'k': 0, 'l': 0}
    n_array = re.findall('[0-9]+', basis)
    l_array = re.findall('[spdfghikl]+', basis)

    for l in range(len(l_array)):
        for l2 in l_array[l]:
            nlmax[l2] = n_array[l]

    return nlmax

def write_rel_configs(rel_configs, fh):
    nnr = 0
    nr = 0
    for nr_config in rel_configs:
        nnr += 1
        fh.write(str(nnr).rjust(4, ' ')+'\n')
        for rel_config in nr_config:
            nr += 1
            count = 0
            if nr <= 9999:
                line = str(nr).rjust(4, ' ')
            else:
                line = str(nr)[-4:].rjust(4, ' ')
            for config in flatten_list(rel_config):  
                if count == 6:
                    line += '\n    '
                    count = 0
                if config[-1] != '0': 
                    line += config + '    '
                    count += 1
            fh.write(line[:-4])
            fh.write('\n')

def convert_list_nr_to_rel_configs(nr_configs):
    rel_configs = []
    for config in nr_configs:
        rel_configs.append(convert_nr_to_rel_config(config))
    return len(flatten_list(rel_configs)), rel_configs


def count_valence(configurations):
    """ Count the number of valence electrons in each basic configuration """
    try:
        configurations_even = configurations['even']
    except KeyError:
        configurations_even = []
    try:
        configurations_odd = configurations['odd']
    except KeyError:
        configurations_odd = []
    configurations = configurations_even + configurations_odd
    if configurations == []:
        print('No configurations were specified')
        sys.exit()

    num_val = sum(int(re.findall('[0-9]+', orbital)[-1]) for orbital in configurations[0].split())

    # Check that each configuration has the same number of valence electrons
    count = 0
    for configuration in configurations:
        count = sum(int(re.findall('[0-9]+', orbital)[-1]) for orbital in configuration.split())
        if count != num_val:
            print('ERROR: ' + configuration + ' has ' + str(count) + ' electrons, but ' 
                    + configurations[0] + ' has ' + str(num_val))
            sys.exit()

    return num_val

def count_excitations(excitations):
    """ Return the multiplicity of excitations given list of excitations"""
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

    # Add orbitals from BASS.INP to table
    try:
        orbitals = read_bass('BASS.INP')
        orb_occ = OrderedDict()
        for orbital in orbitals:
            orb_occ[orbital[0]] = '0', '2'
        nmax = int(re.findall('[0-9]+', next(reversed(orb_occ)))[0])
    # Add orbitals from basis if BASS.INP not found
    except FileNotFoundError as e:
        orb_occ, nmax = expand_basis(basis)

    # Remove orbitals not in basis set
    nlmax = expand_nl(basis)
    orb_occ_copy = copy.deepcopy(orb_occ)
    for orbital in orb_occ_copy:
        n = int(re.findall('[0-9]+', orbital)[0])
        l = re.findall('[spdfghikl]+', orbital)[0]
        nmax = int(nlmax[l])
        if n > nmax:
            orb_occ.pop(orbital)

    # Add orbitals from active orbital list to table
    orb_occ2 = {}
    for orbital in active_orbs:
        min_occ = list(orbital.values())[0].split()[0]
        max_occ = list(orbital.values())[0].split()[1]
        for key in orbital:
            l = re.findall('[spdfghikl]+', key)
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
                if str(n) + l[0] not in list(orb_occ.keys()):
                    orb_occ2[str(n) + l[0]] = min_occ, max_occ
                else:
                    orb_occ[str(n) + l[0]] = min_occ, max_occ
    orb_occ2.update(orb_occ)
    orb_occ = orb_occ2

    # Strip core orbitals from table
    try:
        core_orbs = re.findall('[1-9][spdfghikl]+', orbital_list['core'])
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

    nmin_dict = {'s': 1, 'p': 2, 'd': 3, 'f': 4, 'g': 5, 'h': 6}
    n_array = re.findall('[0-9]+', basis)
    l_array = re.findall('[spdfghikl]+', basis)

    for i in range(len(n_array)):
        for l in l_array[i]:
            if l in nmin_dict:
                for n in range(nmin_dict[l], int(n_array[i]) + 1):
                    orb_occ[str(n) + l] = str(min_occ), str(max_occ)

    return orb_occ, int(max(n_array))

if __name__ == "__main__":
    #print(convert_char_to_digital("22g"))
    #print(convert_digital_to_char("2.2408"))

    #print(convert_nr_to_rel_config('21p1 21d1'))
    #print(convert_nr_to_rel_config('4d10 4f8 6p1'))
    nr_configs = ['4d10 4f9', '4d9 4f9 5s1', '4d9 4f9 5d1', '4d9 4f9 5g1', '4d9 4f9 6s1']
    nr, rel_configs = convert_list_nr_to_rel_configs(nr_configs)
    filename = 'rel_configs.txt'
    f = open(filename, 'w')
    write_rel_configs(rel_configs, f)
    f.close()