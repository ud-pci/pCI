import re
from contextlib import suppress
from copy import deepcopy
from itertools import combinations
import read_bass
import orbitals
import copy

l_val_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
max_occ_dict = {'s': 2, 'p': 6, 'd': 10, 'f': 14, 'g': 18, 'h': 22}

def check_parity(orbitals):
    """ Checks parity of a dictionary of subshells and electron occupancies """
    parity = 0
    for orbital in orbitals:
        l_str = re.findall('[spdfgh]+', orbital)[0]
        if orbitals[orbital] > 0:
            parity += (l_val_dict[l_str] * orbitals[orbital]) % 2
    parity %= 2
    return parity

class Basis:
    """ Class for basis set """
    def __init__(self, basis_str, core_str, active_dict, filename):
        self.name = basis_str
        self.core = core_str.split(" ")
        self.active = []

        try:
            self.read_orbitals_from_bass(filename)
            print('Basis orbitals built from', filename)
        except:
            self.expand_basis()
            print('Basis orbitals built from', basis_str)
        self.remove_core()
        self.remove_unwanted_orbitals(basis_str)
        self.add_active(active_dict)
        self.generate_basis_orbitals()

    def add_active(self, active):
        # Adds active orbitals to basis set
        orb_occ_dict2 = {}
        for orbital in active:
            nmin = int(re.findall('[0-9]+', str(list(orbital.keys())[0]))[0])
            try:
                nmax = int(re.findall('[0-9]+', str(list(orbital.keys())[0]))[1])
            except IndexError:
                nmax = nmin
            min_occ = int(re.findall('[0-9]+', str(list(orbital.values())[0]))[0])
            max_occ = int(re.findall('[0-9]+', str(list(orbital.values())[0]))[1])
            l = re.findall('[spdfgh]+', str(list(orbital.keys())[0]))[0]
            for i in range(nmin,nmax+1):
                if str(i) + l not in list(self.orb_occ_dict.keys()):
                    orb_occ_dict2[(str(i) + l)] = str(min_occ), str(max_occ)
        orb_occ_dict2.update(self.orb_occ_dict)
        self.orb_occ_dict = orb_occ_dict2

    def remove_unwanted_orbitals(self, basis):
        # Remove orbitals not in basis set
        nlmax = orbitals.expand_nl(basis)
        orb_occ_copy = copy.deepcopy(self.orb_occ_dict)
        for orbital in orb_occ_copy:
            n = int(re.findall('[0-9]+', orbital)[0])
            l = re.findall('[spdfgh]+', orbital)[0]
            nmax = int(nlmax[l])
            if n > nmax:
                self.orb_occ_dict.pop(orbital)

    def remove_core(self):
        for orbital in self.core:
            if orbital in self.orb_occ_dict:
                self.orb_occ_dict.pop(orbital)

    def generate_basis_orbitals(self):
        for orbital in list(self.orb_occ_dict.keys()):
            self.active.append(orbital)
    
    def expand_basis(self):
        orb_occ = {}
        min_occ = 0
        max_occ = 2

        nmin_dict = {'s': 1, 'p': 2, 'd': 3, 'f': 4, 'g': 5, 'h': 6}
        n_array = re.findall('[0-9]+', self.name)
        l_array = re.findall('[spdfgh]+', self.name)

        for i in range(len(n_array)):
            for l in l_array[i]:
                if l in nmin_dict:
                    for n in range(nmin_dict[l], int(n_array[i]) + 1):
                        orb_occ[str(n) + l] = str(min_occ), str(max_occ)

        self.orb_occ_dict = orb_occ

    def read_orbitals_from_bass(self, filename):
        """ reads orbitals from BASS.INP to create a dictionary of {key:value} = {subshells: electron occupancies} """
        bass_orbitals = read_bass.read_bass(filename)
        orb_occ = {}
        min_occ = 0
        max_occ = 2
        for orbital in bass_orbitals:
            orb_occ[orbital[0]] = str(min_occ), str(max_occ)

        self.orb_occ_dict = orb_occ

    def print_basis(self):
        print('Basis set:', self.basis_str)
        print('List of orbitals:', [self.active[i].str for i in range(0,len(self.active))])

class Configuration:
    """ Class for configuration """
    def __init__(self, config_str, basis):
        self.name = config_str
        self.orbitals = {}
        self.set_orbitals(basis)
        self.set_parity()

    def print(self):
        print(self.orbitals)
        print('Parity:',  re.sub('(0)', 'Even', re.sub('(1)', 'Odd', str(self.parity))))

    def set_name(self):
        self.name = ''
        for orbital in self.orbitals:
            if self.orbitals[orbital] != 0:
                self.name += orbital + str(self.orbitals[orbital]) + ' '
        self.name = self.name.strip()

    def set_orbitals(self, basis):
        """ create dictionary of {key:value} = {subshells: electron occupancies} """
        for subshell in basis.active:
            self.orbitals[subshell] = 0
        for config in self.name.split(" "):
            self.orbitals[re.findall('[0-9]+', config)[0] + re.findall('[spdfgh]+', config)[0]] += int(re.findall('[0-9]+', config)[1])
        #for orbital in basis.orb_occ_dict:
        #    if basis.orb_occ_dict[orbital][0] == basis.orb_occ_dict[orbital][1]:
        #        self.orbitals.pop(orbital)
        #print(self.orbitals)

    def set_parity(self):
        """ set parity of configuration """
        self.parity = check_parity(self.orbitals)
    
    def add_electron(self, subshell):
        """ add an electron to a subshell """
        # check to see if subshell is filled
        l = re.findall('[spdfgh]+', subshell)[0]
        filled = self.orbitals[subshell] == max_occ_dict[l]
        # If not filled, add electron to that subshell
        if not filled:
            self.orbitals[subshell] += 1
        else:
            raise AttributeError('Subshell ' + subshell + ' is full.')

    def remove_electron(self, subshell):
        """ remove an electron from a subshell """
        # check to see if subshell is empty
        empty = self.orbitals[subshell] == 0
        # If not filled, add electron to that subshell
        if not empty:
            self.orbitals[subshell] -= 1
        else:
            raise AttributeError('Subshell ' + subshell + ' is empty.')

    def excite_electron(self, initial, final):
        """ excite an electron from an initial subshell to a final subshell """
        try:
            self.remove_electron(initial)
            self.add_electron(final)
            self.set_name()
            self.set_parity()
        except Exception as e:
            return e

def get_single_excitations(configuration, basis, basics, singles):
    active = [k for k, v in configuration.orbitals.items() if v > 0]
    for orbital in active:
        for basis_orbital in configuration.orbitals:
            if basis.orb_occ_dict[orbital][0] == basis.orb_occ_dict[orbital][1]:
                continue
            test_config = deepcopy(configuration)
            e = test_config.excite_electron(orbital, basis_orbital)
            if (e != None): continue
            if configuration.parity == test_config.parity and test_config.name not in basics and test_config.name not in singles:
                singles.append(test_config.name)

def get_double_excitations(configuration, basis, basics, singles, doubles):
    orbitals = [k for k, v in configuration.orbitals.items() if v > 0]
    num_electrons = [v for k, v in configuration.orbitals.items() if v > 0]
    active = []
    configurations = list(set(basics + singles))
    i = 0
    for orbital in orbitals:
        for num in range(num_electrons[i]):
            active.append(orbital)
        i += 1
    active = list(set(combinations(active, 2)))

    for orbital1, orbital2 in active:
        for basis_orbital in configuration.orbitals:
            for basis_orbital2 in configuration.orbitals:
                if basis.orb_occ_dict[orbital1][0] == basis.orb_occ_dict[orbital1][1]: continue
                if basis.orb_occ_dict[orbital2][0] == basis.orb_occ_dict[orbital2][1]: continue
                test_config = deepcopy(configuration)
                e = test_config.excite_electron(orbital1, basis_orbital)
                f = test_config.excite_electron(orbital2, basis_orbital2)
                if (e != None or f != None): continue
                if test_config.parity == configuration.parity and test_config.name not in configurations and test_config.name not in doubles:
                    doubles.append(test_config.name)    

def get_triple_excitations(configuration, basis, basics, singles, doubles, triples):
    orbitals = [k for k, v in configuration.orbitals.items() if v > 0]
    num_electrons = [v for k, v in configuration.orbitals.items() if v > 0]
    active = []
    configurations = list(set(basics + singles))
    i = 0
    for orbital in orbitals:
        for num in range(num_electrons[i]):
            active.append(orbital)
        i += 1
    active = list(set(combinations(active, 3)))

    for orbital1, orbital2, orbital3 in active:
        for basis_orbital in configuration.orbitals:
            for basis_orbital2 in configuration.orbitals:
                for basis_orbital3 in configuration.orbitals:
                    if basis.orb_occ_dict[orbital1][0] == basis.orb_occ_dict[orbital1][1]: continue
                    if basis.orb_occ_dict[orbital2][0] == basis.orb_occ_dict[orbital2][1]: continue
                    if basis.orb_occ_dict[orbital3][0] == basis.orb_occ_dict[orbital3][1]: continue
                    test_config = deepcopy(configuration)
                    e = test_config.excite_electron(orbital1, basis_orbital)
                    f = test_config.excite_electron(orbital2, basis_orbital2)
                    g = test_config.excite_electron(orbital3, basis_orbital3)
                    if (e != None or f != None or g != None): continue
                    if test_config.parity == configuration.parity and test_config.name not in configurations and test_config.name not in triples:
                        triples.append(test_config.name)    

def get_excitations(basics, basis):
    configurations = []
    singles = []
    doubles = []
    triples = []
    configurations.append(basics)
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_single_excitations(config, basis, basics, singles)
    configurations.append(singles)
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_double_excitations(config, basis, basics, singles, doubles)
    configurations.append(doubles)
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_triple_excitations(config, basis, basics, singles, doubles, triples)
    configurations.append(triples)

    if config.parity == 0:
        write_nr_configs('nr_configs_even.txt', basis, basics, singles, doubles, triples)
    else:
        write_nr_configs('nr_configs_odd.txt', basis, basics, singles, doubles, triples)
    
    return configurations

def write_nr_configs(fh, basis, basics, singles, doubles, triples):
    f = open(fh, 'w')
    f.write("Basics:\n")
    for configuration in range(len(basics)):
        f.write(str(basics[configuration])+"\n")
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_single_excitations(config, basis, basics, singles)
    f.write("\nSingles:\n")
    for configuration in range(len(singles)):
        f.write(str(singles[configuration])+"\n")
    f.write("\nDoubles:\n")
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_double_excitations(config, basis, basics, singles, doubles)
    for configuration in range(len(doubles)):
        f.write(str(doubles[configuration])+"\n")
    f.write("\nTriples:\n")
    for configuration in basics:
        config = Configuration(configuration, basis)
        get_triple_excitations(config, basis, basics, singles, doubles, triples)
    for configuration in range(len(triples)):
        f.write(str(triples[configuration])+"\n")
    f.write('\nSUMMARY: \n')
    f.write(str(len(basics)) + ' basic configurations\n')
    f.write(str(len(singles)) + ' singly excited configurations\n')
    f.write(str(len(doubles)) + ' doubly excited configurations\n')
    f.write(str(len(triples)) + ' triply excited configurations\n')
    f.write('Total number of configurations: ' + str(len(basics) + len(singles) + len(doubles) + len(triples)))
    f.close()

if __name__ == "__main__":
    basis_set = '6sp5d4f'
    core = '1s 2s 2p 3s 3p 3d 4s 4p 4d'
    basics = ['5s2 5p6 4f12 6s2 5d1']
    active = [{'5s': '2  2'}, {'5p':  '6  6'}, {'4f': '9 14'}, {'6s':  '0  2'}, {'5d': '0  4'}, {'6p':  '0  4'}]
    singles = []
    doubles = []
    triples = []
    basis = Basis(basis_set, core, active, None)
    get_excitations(basics, basis)
