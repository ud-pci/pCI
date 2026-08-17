""" ci

This script allows the user to automate the configuration list construction from inputted parameters in a "config.yml" file.
The "config.yml" file should have the following blocks:

    * system - system parameter (bin_directory, run_codes, on_hpc)
    * basis.orbitals.core - core orbitals are read to ensure consistent basis set
    * ci - parameters used by ci programs (reference configurations, basis set, orbitals, excitations)
    * optional - optional parameters (isotope shifts, code methods, running all-order codes, pci versions)

From these parameters, this script will create all input files required for execution of the add program.
After the input files are created, the add program will be executed to create the list of configurations CONF.INP.
If optional.generate_directories is set to "True", the script will generate respective directories for CI calculations (e.g. /even and /odd)
If optional.run_codes is set to "True", the script will then submit the job script in the respective directories.

This python script has 3 main capabilities for configuration list construction:
1. General construction of configuration lists in root directory
2. Construction of configuration lists for isotope shift calculations
3. Construction of configuration lists for multiple code methods

"""
import yaml
import os
import sys
import re
import math
import orbitals as orb_lib
import get_atomic_data as libatomic
from utils import run_shell, get_dict_value, check_slurm_installed
from pathlib import Path
from gen_job_script import write_job_script

def read_yaml(filename):
    """
    This function reads a configuration file in YAML format and returns a dictionary of config parameters
    """

    with open(filename, 'r') as f:
        config = yaml.safe_load(f)

    return config

def write_add_inp(filename, Z, AM, config, multiplicity, num_val, orb_occ, parity, J, JM, J_selection, num_energy_levels, num_dvdsn_iterations, K_is, C_is):
    """Write ADD.INP file for specified parity"""

    configurations = config['add']['ref_configs']
    # Check if there are configurations
    if configurations[parity] == []:
        print(parity.capitalize(), 'parity configurations were not specified')
        return

    # Define name of file to export
    filename = f'{filename[:-4]}{parity}{filename[-4:]}'

    with open(filename, 'w') as f:
        # Write header
        f.write(f'Ncor= {len(configurations[parity]):d}\n')
        f.write(f'NsvNR{orb_occ["num_orb"]:d}\n')
        f.write(f'mult= {multiplicity:d}\n')
        f.write(f' NE = {num_val:d}\n\n')

        # Reformat configurations
        configurations_formatted = format_configurations(configurations[parity])
        max_orb_per_line = 6

        # Write list of basic configurations
        for configuration in configurations_formatted:
            num_lines = math.ceil(len(configuration)/max_orb_per_line)
            for line in range(num_lines):
                f.write(f"L:  {' '.join(configuration[line*max_orb_per_line:(line+1)*max_orb_per_line])}\n")
        f.write('\n')

        # Write list of orbitals and occupation numbers
        count = 0

        for orb in orb_occ:
            orb_occ_formatted = format_orb_occ(orb, orb_occ[orb])
            f.write(f'  {orb_occ_formatted}')
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
        else:
            core_formatted = []

        atom = config['atom']
        # Write head of CONF.INP
        f.write('>>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>\n')
        f.write(f'  {atom["name"]} {parity}\n')
        f.write(f'  Z = {Z}\n')
        f.write(f' Am = {round(AM):.1f}\n')
        f.write(f'  J = {J}\n')
        f.write(f' Jm = {JM}\n')
        f.write(f' Nso=  {num_core_orb}\n')
        f.write(' Nc =   10 \n')
        if J_selection:
            f.write(' Kv =  3 \n')
        else:
            f.write(' Kv =  4 \n')
        f.write(f' Nlv=  {num_energy_levels}\n')
        f.write(f' Ne =  {num_val}\n')
        f.write(' Kl4=  1 \n')
        f.write(' Nc4=999 \n')
        f.write('Crt4= 0.0001 \n')
        f.write('kout= 0 \n')
        f.write('Ncpt= 0 \n')
        f.write('Cut0= 0.0001 \n')
        f.write(f'N_it= {num_dvdsn_iterations}\n')
        if atom['include_breit']:
            kbrt = 2
        else:
            kbrt = 0
        f.write(f'Kbrt= {kbrt}\n')
        if K_is > 0:
            f.write(f'K_is= {K_is}\n')
            f.write(f'C_is= {C_is}\n')

        # Write core orbitals
        if num_core_orb != 0:
            count = 0
            for i in range(num_core_orb):
                orb = core_formatted[i]
                if orb[0] == '-':
                    f.write(f'    {orb}')
                else:
                    f.write(f'     {orb}')
                count += 1
                if count == 6:
                    f.write('\n')
                    count = 0

        f.write('\n    ')

    print(f'{filename} has been written')

def check_orbital_exists(orb):
    """Exit with an error if orbital n < l+1 (e.g. 3f, 4g do not exist)."""
    m = re.match(r'(\d+)([a-z])', orb)
    if not m:
        return
    n, l_char = int(m.group(1)), m.group(2)
    if l_char not in orb_lib.l_dict:
        return
    l = orb_lib.l_dict[l_char]
    if n < l + 1:
        print(f'ERROR: orbital {orb} does not exist (min n for {l_char} is {l + 1})')
        sys.exit()

def format_configurations(configurations):
    """ Format configurations as 'nnlqq',
        where nn represents principal quantum number
               l represents angular quantum number
              qq represents occupation number """
    configurations_formatted = []
    for configuration in configurations:
        orbitals_formatted = []
        for orbital in configuration.split():
            check_orbital_exists(orbital)
            num_digits = len(re.findall('[0-9]+', orbital)[0])
            if num_digits == 1:
                orbital_formatted = orbital.center(5, ' ')
            elif num_digits == 2:
                orbital_formatted = orbital.ljust(5, ' ')
            else:
                print(f'ERROR: orbital {orbital} exceeds n=99')
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
        orb_occ_formatted = f'{orb.rjust(3)}{occ[0].rjust(3)}{occ[1].rjust(3)}'

    return orb_occ_formatted

def create_add_inp(config, parity, J, JM, J_selection, num_energy_levels, num_dvdsn_iterations, K_is, c, j_suffix=None):
    # Get atomic data
    atom = config['atom']

    try:
        Z, AM, *_ = libatomic.get_atomic_data(atom['name'], atom['isotope'])
    except KeyError:
        Z, AM, *_ = libatomic.get_atomic_data(atom['name'], '')

    add = config['add']
    num_val = orb_lib.count_valence(add['ref_configs'])
    multiplicity = orb_lib.count_excitations(add['excitations'])
    orb_occ = orb_lib.expand_orbitals(add['basis_set'], add['ref_configs'], add['orbitals'], add['excitations'])

    if add['ref_configs'][parity]:
        write_add_inp('ADD.INP', Z, AM, config, multiplicity, num_val, orb_occ, parity, J, JM, J_selection, num_energy_levels, num_dvdsn_iterations, K_is, c)
        if j_suffix is not None:
            old_name = f'ADD{parity}.INP'
            new_name = f'ADD{parity}{j_suffix}.INP'
            if os.path.isfile(old_name):
                os.rename(old_name, new_name)
    else:
        print(f'no {parity} reference configurations specified')

def form_conf_inp(parity, bin_dir, j_suffix=None):
    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'
    suffix = '' if j_suffix is None else str(j_suffix)
    run_shell(f'cp ADD{parity}{suffix}.INP ADD.INP')
    run_shell(f'{bin_dir}add < add.in > add{parity}{suffix}.out')
    print(f'output of add saved to add{parity}{suffix}.out')
    run_shell(f'cp CONF.INP CONF{parity}{suffix}.INP')
    print(f'CONF{parity}{suffix}.INP created')
    if j_suffix is not None:
        run_shell('rm ADD.INP')

def write_ci_in(conf_dir, write_hij, include_lsj):
    Kw = '1' if write_hij else '0'
    KLSJ = '1' if include_lsj else '0'
    with open(f'{conf_dir}/ci.in', 'w') as f:
        if os.path.isfile('basis/SGC.CON'):
            f.write('Kl = 2 \n')
            f.write('Ksig = 2 \n' if os.path.isfile('basis/SCRC.CON') else 'Ksig = 1 \n')
        else:
            f.write('Kl = 0 \n')
            f.write('Ksig = 0 \n')
        f.write('Kdsig = 0 \n')
        f.write(f'Kw = {Kw}\n')
        f.write(f'KLSJ = {KLSJ}')

def move_conf_inp(conf_dir, root_dir, parity, run_codes, include_lsj, write_hij, K_is, C_is, j_suffix=None):
    if not os.path.isdir(conf_dir):
        run_shell(f'mkdir {conf_dir}')
    if os.path.isfile('basis/HFD.DAT'):
        run_shell(f'cp basis/HFD.DAT {conf_dir}')
    else:
        print('basis/HFD.DAT was not found.. try running basis.py first')
        sys.exit()

    if j_suffix is not None:
        conf_src = f'{root_dir}/CONF{parity}{j_suffix}.INP'
        add_src  = f'{root_dir}/ADD{parity}{j_suffix}.INP'
    elif K_is > 0:
        conf_src = f'{root_dir}/CONF{parity}{C_is}.INP'
        add_src  = f'{root_dir}/ADD{parity}{C_is}.INP'
    else:
        conf_src = f'{root_dir}/CONF{parity}.INP'
        add_src  = f'{root_dir}/ADD{parity}.INP'

    if os.path.isfile(conf_src):
        run_shell(f'cp {conf_src} {conf_dir}/CONF.INP')
    if os.path.isfile(add_src):
        run_shell(f'cp {add_src} {conf_dir}/ADD.INP')
    if os.path.isfile('basis/SGC.CON'):
        run_shell(f'cp basis/SGC.CON {conf_dir}')
    if os.path.isfile('basis/SCRC.CON'):
        run_shell(f'cp basis/SCRC.CON {conf_dir}')
    if run_codes and os.path.isfile(f'{root_dir}/ci.qs'):
        run_shell(f'cp {root_dir}/ci.qs {conf_dir}')

    write_ci_in(conf_dir, write_hij, include_lsj)

def submit_ci_job(conf_path, script_name, submit_job):
    os.chdir(conf_path)
    if script_name and submit_job:
        run_shell(f'sbatch {script_name}')
    else:
        print('job script was not submitted. check job script and submit manually.')
    os.chdir('../')

if __name__ == '__main__':
    # Read yaml file for system configurations
    yml_file = input('Input yml-file: ')
    config = read_yaml(yml_file)

    system = get_dict_value(config, 'system')
    gen_dir = get_dict_value(system, 'generate_directories')
    run_codes = get_dict_value(system, 'run_codes')
    pci_version = get_dict_value(system, 'pci_version')
    bin_dir = get_dict_value(system, 'bin_directory')
    on_hpc = get_dict_value(system, 'on_hpc')

    # hpc parameters
    on_slurm = check_slurm_installed()
    if on_hpc and on_slurm:
        hpc = get_dict_value(config, 'hpc')
        if hpc:
            submit_job = get_dict_value(hpc, 'submit_job')
            partition = get_dict_value(hpc, 'partition')
            nodes = get_dict_value(hpc, 'nodes')
            tasks_per_node = get_dict_value(hpc, 'tasks_per_node')
        else:
            print('hpc block was not found in', yml_file)
            submit_job = False
            partition, nodes, tasks_per_node = None, 1, 1
    else:
        on_hpc = False
        submit_job = False

    optional = get_dict_value(config, 'optional')
    isotope_shifts = get_dict_value(optional, 'isotope_shifts')
    include_isotope_shifts = get_dict_value(isotope_shifts, 'include')
    if include_isotope_shifts:
        K_is = get_dict_value(isotope_shifts, 'K_is')
        C_is = get_dict_value(isotope_shifts, 'C_is')
        c_list = [-C_is, -C_is/2, 0, C_is/2, C_is]
        K_is_dict = {0: '', 1: 'FS', 2: 'SMS', 3: 'NMS', 4: 'MS'}

    atom = get_dict_value(config, 'atom')
    code_method = get_dict_value(atom, 'code_method')
    conf = get_dict_value(config, 'conf')
    for_portal = get_dict_value(system, 'for_portal')

    include_lsj = get_dict_value(conf, 'include_lsj')
    write_hij = get_dict_value(conf, 'write_hij')

    # Check if user wants to generate directories for CI computations
    gen_dir = gen_dir or run_codes

    # CONF.INP should only need to be constructed once, then copied to respective directories
    # Read input to add from add.in if it exists, otherwise create it
    if not os.path.isfile('add.in'):
        with open('add.in', 'w') as f:
            f.write('0 \n 0')

    # BASS.INP dictates the order of orbitals in ADD.INP.
    # Find a BASS.INP file from the basis directory.
    if not isinstance(code_method, list):
        if include_isotope_shifts and K_is > 0:
            run_shell(f'cp ./{K_is_dict[K_is]}/0/basis/BASS.INP .')
            code_method = ['']
        else:
            run_shell('cp ./basis/BASS.INP .')
            print('BASS.INP taken from basis/')
    else:
        if os.path.isdir('ci+all-order'):
            run_shell('cp ./ci+all-order/basis/BASS.INP .')
            print('BASS.INP taken from ci+all-order/basis/')
        elif os.path.isdir('ci+second-order'):
            run_shell('cp ./ci+second-order/basis/BASS.INP .')
            print('BASS.INP taken from ci+second-order/basis/')

    # Run script to generate ADD.INP for even and odd configurations
    J_odd = conf['odd']['J']
    JM_odd = conf['odd']['JM']
    J_even = conf['even']['J']
    JM_even = conf['even']['JM']

    J_selection_odd = conf['odd']['J_selection']
    num_energy_levels_odd = conf['odd']['num_energy_levels']
    num_dvdsn_iterations_odd = conf['odd']['num_dvdsn_iterations']

    J_selection_even = conf['even']['J_selection']
    num_energy_levels_even = conf['even']['num_energy_levels']
    num_dvdsn_iterations_even = conf['even']['num_dvdsn_iterations']

    if include_isotope_shifts:
        for c in c_list:
            create_add_inp(config, 'even', J_even, JM_even, J_selection_even, num_energy_levels_even, num_dvdsn_iterations_even, K_is, c)
            create_add_inp(config, 'odd', J_odd, JM_odd, J_selection_odd, num_energy_levels_odd, num_dvdsn_iterations_odd, K_is, c)
            even_exists = os.path.isfile('ADDeven.INP')
            odd_exists = os.path.isfile('ADDodd.INP')
            if even_exists:
                run_shell(f'cp ADDeven.INP ADDeven{c}.INP')
            if odd_exists:
                run_shell(f'cp ADDodd.INP ADDodd{c}.INP')
            # Run add and form CONF.INP from respective ADD.INP files
            parities = []
            if even_exists:
                form_conf_inp('even', bin_dir)
                run_shell(f'cp CONFeven.INP CONFeven{c}.INP')
                parities.append('even')
            if odd_exists:
                form_conf_inp('odd', bin_dir)
                run_shell(f'cp CONFodd.INP CONFodd{c}.INP')
                parities.append('odd')
    else:
        if for_portal:
            j_values = sorted({int(J) for J in [J_even, J_odd]})
            print(f'Portal mode enabled - generating configurations for J={j_values} for both parities')
            parities = []

            # Generate ADD files for all combinations
            for parity in ['even', 'odd']:
                parity_config = conf[parity]
                if config['add']['ref_configs'][parity]:
                    parities.append(parity)
                    for J in j_values:
                        create_add_inp(config, parity, J, J,
                                       parity_config['J_selection'],
                                       parity_config['num_energy_levels'],
                                       parity_config['num_dvdsn_iterations'],
                                       0, 0, j_suffix=J)

            # Generate CONF files from ADD files
            for parity in parities:
                for J in j_values:
                    if os.path.isfile(f'ADD{parity}{J}.INP'):
                        form_conf_inp(parity, bin_dir, j_suffix=J)
        else:
            # Regular mode: use specified J values for each parity
            create_add_inp(config, 'even', J_even, JM_even, J_selection_even, num_energy_levels_even, num_dvdsn_iterations_even, 0, 0)
            create_add_inp(config, 'odd', J_odd, JM_odd, J_selection_odd, num_energy_levels_odd, num_dvdsn_iterations_odd, 0, 0)
            # Check if ADDeven.INP and ADDodd.INP exist
            even_exists = os.path.isfile('ADDeven.INP')
            odd_exists = os.path.isfile('ADDodd.INP')
            # Run add and form CONF.INP from respective ADD.INP files
            parities = []
            if even_exists:
                form_conf_inp('even', bin_dir)
                parities.append('even')
            if odd_exists:
                form_conf_inp('odd', bin_dir)
                parities.append('odd')

    # Create a ci.qs job script if it doesn't exist yet
    if on_hpc:
        print(f'generating new ci.qs in {os.getcwd()} directory')
    script_name = write_job_script('.', 'ci', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir) if on_hpc else None

    # Copy ADD.INP and CONF.INP to all directories if gen_dir == True
    if gen_dir:
        root_dir = os.getcwd()
        if include_isotope_shifts and K_is > 0:
            for method in code_method:
                dir_path = os.getcwd()
                is_dir = f'{method}/{K_is_dict[K_is]}'
                Path(f'{dir_path}/{is_dir}').mkdir(parents=True, exist_ok=True)
                os.chdir(f'{dir_path}/{is_dir}')
                for c in c_list:
                    dir_path = os.getcwd()
                    if c < 0:
                        dir_prefix = 'minus'
                    elif c > 0:
                        dir_prefix = 'plus'
                    else:
                        dir_prefix = ''
                    dir_name = f'{dir_prefix}{abs(c)}'
                    Path(f'{dir_path}/{dir_name}').mkdir(parents=True, exist_ok=True)
                    os.chdir(dir_name)
                    print(dir_name)
                    for parity in parities:
                        J = get_dict_value(conf[parity], 'J')
                        conf_path = f'{parity}{str(J)[0]}'
                        # check if run has already been completed
                        if os.path.isfile(f'{conf_path}/FINAL.RES'):
                            print(f'Completed run detected in {conf_path}')
                            print('Skipping..')
                            continue
                        move_conf_inp(conf_path, root_dir, parity, run_codes, include_lsj, write_hij, K_is, c)
                        if on_hpc and run_codes:
                            submit_ci_job(conf_path, script_name, submit_job)
                        else:
                            print('run_codes option is only available with HPC access')
                    os.chdir('../')
                if K_is_dict[K_is]:
                    os.chdir('../')
        else:
            if not isinstance(code_method, list):
                code_method = [None]
            dir_path = os.getcwd()
            for method in code_method:
                if method is not None:
                    Path(f'{dir_path}/{method}/basis').mkdir(parents=True, exist_ok=True)
                    os.chdir(method)
                if for_portal:
                    for parity in parities:
                        for J in j_values:
                            conf_path = f'{parity}{J}'
                            if os.path.isfile(f'{root_dir}/CONF{parity}{J}.INP'):
                                move_conf_inp(conf_path, root_dir, parity, run_codes, include_lsj, write_hij, 0, 0, j_suffix=J)
                                if on_hpc and run_codes:
                                    submit_ci_job(conf_path, script_name, submit_job)
                                else:
                                    print('run_codes option is only available with HPC access - please run ci codes manually')
                else:
                    for parity in parities:
                        J = get_dict_value(conf[parity], 'J')
                        conf_path = f'{parity}{str(J)[0]}'
                        move_conf_inp(conf_path, root_dir, parity, run_codes, include_lsj, write_hij, 0, 0)
                        if on_hpc and run_codes:
                            submit_ci_job(conf_path, script_name, submit_job)
                        else:
                            print('run_codes option is only available with HPC access')
                if method is not None:
                    os.chdir('../')

    # Cleanup - remove add.in, ADD.INP, CONF.INP and CONF_.INP from root directory
    run_shell('rm add.in ADD*.INP CONF*.INP add*.out BASS.INP')

    print('add script completed')
