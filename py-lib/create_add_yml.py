import sys
import re

if __name__ == "__main__":
    system = {}

    f = open('add_test.yml', 'w')

    # Read system parameters
    system['name'] = str(input("Enter the name of the system: "))
    system['atomic_number'] = str(input("Enter the atomic number: "))
    system['atomic_mass'] = str(input("Enter the atomic mass: "))
    system['include_gaunt'] = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Include gaunt? "))))
    system['include_breit'] = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input ("Include breit? "))))

    # Write system parameters
    f.write('system:\n  [\n')
    f.write('    name: ' + str(system['name'] + ',\n'))
    f.write('    atomic_number: ' + str(system['atomic_number'] + ',\n'))
    f.write('    atomic_mass: ' + str(system['atomic_mass'] + ',\n'))
    f.write('    include_gaunt: ' + str(system['include_gaunt'] + ',\n'))
    f.write('    include_breit: ' + str(system['include_breit'] + ',\n'))
    f.write('  ]\n\n')

    # Read configurations
    configurations = {}
    configurations['odd'] = filter(None,str(input("Enter your list of odd configurations separated by commas: ")).split(','))
    configurations['even'] = filter(None,str(input("Enter your list of even configurations separated by commas: ")).split(','))

    # Write configurations
    f.write('configurations:\n  odd:\n    [\n')
    for configuration in configurations['odd']:
        f.write('      ' + configuration.strip() + ',\n')
    f.write('    ]\n  even:\n    [\n')
    for configuration in configurations['even']:
        f.write('      ' + configuration.strip() + ',\n')
    f.write('    ]\n\n')

    # Read basis
    basis = str(input("Enter your basis set (e.g. 12spd8fg): "))

    # Write basis
    f.write('basis: ' + basis + '\n\n')

    # Read orbital information
    orbitals = {}
    orbitals['core'] = str(input("Enter core subshells: "))
    orbitals['active'] = filter(None,str(input("Enter your list of active subshells along with min/max occupation numbers (e.g. 4-5s: 0 2, 4-5p: 0 4): ")).split(','))

    # Write orbital information
    f.write('orbitals:\n')
    f.write('  core: ' + orbitals['core'] + '\n')
    f.write('  active:\n    [\n')
    for orbital in orbitals['active']:
        f.write('      ' + orbital.strip() + ',\n')
    f.write('    ]\n\n')

    # Read excitation information
    excitations = {}
    excitations['single'] = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Include single excitations? "))))
    excitations['double'] = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Include double excitations? "))))
    excitations['triple'] = re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Include triple excitations? "))))

    # Write excitation information
    f.write('excitations:\n  [\n')
    f.write('    single: ' + excitations['single'] + ',\n')
    f.write('    double: ' + excitations['double'] + ',\n')
    f.write('    triple: ' + excitations['triple'] + ',\n')
    f.write('  ]\n')

    