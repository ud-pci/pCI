import yaml
import re
import sys
import orbitals as orb_lib
import get_excitations
import math
import get_atomic_data as libatomic

def read_yaml(filename):
    """
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """
    with open(filename,'r') as f:
        config = yaml.safe_load(f)

    return {**config['system'], **config['add'], **config['conf']}
    
def write_add_inp(filename, Z, AM, system, multiplicity, num_val, orb_occ, parity):
    """Write ADD.INP file for specified parity"""

    configurations = system['ref_configs']
    orbitals = system['orbitals']
    # Check if there are configurations 
    if configurations[parity] == []:
        print(parity.capitalize(), 'parity configurations were not specified')
        return

    # Define name of file to export
    filename = filename[0:-4] + parity + filename[-4:]
    f = open(filename,'w')

    # Write header
    f.write('Ncor= {:d}\n'.format(len(configurations[parity])))
    f.write('NsvNR{:d}\n'.format(orb_occ['num_orb']))
    f.write('mult= {:d}\n'.format(multiplicity))
    f.write(' NE = {:d}\n\n'.format(num_val))

    # Reformat configurations
    configurations_formatted = format_configurations(configurations[parity])
    max_orb_per_line = 6

    # Write list of basic configurations
    for configuration in configurations_formatted:
        num_lines = math.ceil(len(configuration)/max_orb_per_line)
        for line in range(num_lines):
            f.write('L:  ' + ' '.join(configuration[line*max_orb_per_line:(line+1)*max_orb_per_line]) + '\n')
    f.write('\n')

    # Write list of orbitals and occupation numbers
    nmin = orb_occ['nmin']
    nmax = orb_occ['nmax']
    count = 0

    for orb in orb_occ:
        orb_occ_formatted = format_orb_occ(orb, orb_occ[orb])
        f.write('  ' + orb_occ_formatted)
        count += 1
        if count == 6:
            f.write('\n')
            count = 0

    f.write('\n')

    num_core_orb = 0
    if orb_occ['core'] != 0:
        try:
            core_formatted = orb_lib.convert_char_to_digital(orb_occ['core'])
            num_core_orb = len(core_formatted)
        except Exception as e:
            print(e)
            pass
    else:
        core_formatted = []
        
    # Write head of CONF.INP
    f.write('>>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>\n')
    f.write('  ' + system['name'] + ' ' + parity + '\n')
    f.write('  Z = ' + str(Z) + '\n')
    f.write(' Am = ' + "{:.1f}".format(round(AM)) + '\n')
    f.write('  J = ' + str(system['J']) + '\n')
    f.write(' Jm = ' + str(system['JM']) + '\n')
    f.write(' Nso=  ' + str(num_core_orb) + '\n')
    f.write(' Nc =   10 \n')
    f.write(' Kv =  4 \n')
    f.write(' Nlv=  ' + str(system['num_energy_levels']) + '\n')
    f.write(' Ne =  ' + str(num_val) + '\n')
    f.write(' Kl4=  1 \n')
    f.write(' Nc4=999 \n')
    f.write(' Gj = 0.0000 \n')
    f.write('Crt4= 0.0001 \n')
    f.write('kout= 0 \n')
    f.write('Ncpt= 0 \n')
    f.write('Cut0= 0.0001 \n')
    f.write('N_it= ' + str(system['num_dvdsn_iterations']) + '\n')
    if system['include_breit']:
        kbrt = 2
    else:
        kbrt = 0
    f.write('Kbrt= ' + str(kbrt) + '\n')
    
    # Write core orbitals
    if num_core_orb != 0:
        count = 0
        for i in range(0, num_core_orb):
            orb = core_formatted[i]
            if orb[0] == '-':
                f.write('    ' + orb)
            else:
                f.write('     ' + orb)
            count += 1
            if count == 6:
                f.write('\n')
                count = 0
    
    f.write('\n    ')

    print(filename + ' has been written')
    f.close()

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
    ignore_list = ['num_orb', 'core', 'nmin', 'nmax']
    if orb in ignore_list:
        orb_occ_formatted = ''
    else:
        orb_occ_formatted = orb.rjust(3, ' ') + occ[0].rjust(3, ' ') + occ[1].rjust(3, ' ')

    return orb_occ_formatted

if __name__ == "__main__":
    filename = input('Enter input file: ')
    try:
        config = read_yaml(filename)
    except FileNotFoundError:
        print('ERROR: The file', filename, 'was not found')
        sys.exit()
    except yaml.scanner.ScannerError:
        print('ERROR: The file', filename, 'was not in YAML format')
        sys.exit()

    # Get atomic data
    Z, AM, symbol, cfermi, rnuc, num_rem_ele = libatomic.get_atomic_data(config['name'], config['isotope'])

    num_val = orb_lib.count_valence(config['ref_configs'])
    orb_occ = orb_lib.expand_orbitals(config['basis_set'], config['orbitals'])
    multiplicity = orb_lib.count_excitations(config['excitations'])

    if config['ref_configs']['even']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'even')
    else:
        print('no even reference configurations specified')
    if config['ref_configs']['odd']:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, 'odd')
    else:
        print('no odd reference configurations specified')
