""" dtm

Automates setup and optional execution of pdtm transition matrix element calculations from a config.yml file.

Required config.yml blocks:

    system:
        bin_directory   - path to directory containing pdtm binary
        run_codes       - if True, run or submit pdtm after setup
        on_hpc          - if True, use Slurm for job submission
    atom:
        code_method     - one or more methods (CI, CI+MBPT, CI+all-order)
    basis:
        orbitals.core   - space-separated list of core orbitals (used to set Nso in MBPT.INP)
    dtm:
        include_rpa     - if True, run pdtm in Init mode first to generate DTM.INT, then run rpa -> rpa_dtm -> pdtm
        DM:
            matrix_elements - operators for diagonal matrix elements
            level_range:
                even / odd  - "from_level to_level" range strings;
                              a dm_even or dm_odd directory is created only if the corresponding range is present
        TM:
            matrix_elements - operators (e.g. E1, M1, E2, ...) for non-portal mode
            from:
                parity      - parity of the bra CI set
                level_range - "from_level to_level" range string
            to:
                parity      - parity of the ket CI set
                level_range - "from_level to_level" range string

Portal mode (system.for_portal = True):
    Ignores DM and TM level_range settings. Instead creates four TM directories
    covering all parity combinations, each with all Nlv levels read from CONF.INP:

        tm_even0_odd1  : E1, M2, E3
        tm_odd0_even1  : E1, M2, E3
        tm_even0_even1 : M1, E2, M3
        tm_odd0_odd1   : M1, E2, M3

    Requires even0, even1, odd0, odd1 CI directories with CONF.RES or
    CONFFINAL.RES present.

For each dtm directory, the script writes dtm.in (Mode = TM or DM, level ranges, operator list), 
MBPT.INP (RPA flags), and copies CI files (CONF.INP, CONF.DET, CONF.XIJ, CONFSTR.RES, CONF.DAT, CONF.INT). 
If on_hpc is True, a Slurm job script is also written and optionally submitted.
"""
import yaml
import os
import sys
from pathlib import Path
from utils import run_shell, get_dict_value, check_slurm_installed
from gen_job_script import write_job_script


TM_DIR_PATHS = {
    'tm_even0_odd1':  ('even0', 'odd1'),
    'tm_odd0_even1':  ('odd0', 'even1'),
    'tm_even0_even1': ('even0', 'even1'),
    'tm_odd0_odd1':   ('odd0', 'odd1'),
}

TM_DIR_OPERATORS = {
    'tm_even0_odd1':  ['E1', 'M2', 'E3'],
    'tm_odd0_even1':  ['E1', 'M2', 'E3'],
    'tm_even0_even1': ['M1', 'E2', 'M3'],
    'tm_odd0_odd1':   ['M1', 'E2', 'M3'],
}


def read_yaml(filename):
    with open(filename, 'r') as f:
        config = yaml.safe_load(f)
    return config

def write_mbpt_inp(basis, matrix_elements):
    core_orbitals = basis['orbitals']['core']
    Nso = 0
    for orbital in core_orbitals.split(' '):
        if orbital[-1] == 's':
            Nso += 1
        else:
            Nso += 2

    key_dict = {
        'A_hf': '0', 'B_hf': '0', 'E1_L': '0', 'EDM': '0', 'PNC': '0',
        'E1_V': '0', 'AM':  '0', 'MQM':  '0', 'M1':  '0', 'E2':  '0',
        'E3':   '0', 'M2':  '0', 'M3':   '0',
    }
    for matrix_element in matrix_elements:
        key_dict[matrix_element] = '1'
        if matrix_element == 'E1':
            key_dict['E1_L'] = '1'
            key_dict['E1_V'] = '1'

    with open('MBPT.INP', 'w') as f:
        f.write('MBPT>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Used by ALL MBPT programs \n')
        f.write('Nso = ' + str(Nso) + ' - CI core \n')
        f.write('Nsh = ' + str(Nso) + ' - defines SCF field (For MS calculations Nsh=Nso) \n')
        f.write('Nss =999 \n')
        f.write('Nsv = ' + str(Nso+1) + ' - =Nso+1 \n')
        f.write('Nmax=210 \n')
        f.write('Lmax=  4 - max (l_i,l_k) for valence radial integrals \n')
        f.write('Kmax=  9 - max multipolarity of two-electron valence integrals \n')
        f.write('Kt  =  1 - Keep this fixed \n')
        f.write('Kbrt=  2 - Breit interaction \n')
        f.write('Kout=  0 - Details in output file \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')
        f.write('SMS:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> (p_i dot p_k) interaction \n')
        f.write('C_sms= 0.00000   - scaling of SMS interaction \n')
        f.write('Klow=  1         - lower component included/ignored \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')
        f.write('RPA >>>>>>>>>>>>>>>>>>Used by programs which deal with RPA MEs \n')
        f.write('A_hf   ' + key_dict['A_hf'] + ' | \n')
        f.write('B_hf   ' + key_dict['B_hf'] + ' | \n')
        f.write('E1_L   ' + key_dict['E1_L'] + ' | \n')
        f.write('EDM    ' + key_dict['EDM']  + ' | - Right hand side operators \n')
        f.write('PNC    ' + key_dict['PNC']  + ' | \n')
        f.write('E1_V   ' + key_dict['E1_V'] + ' | \n')
        f.write('AM     ' + key_dict['AM']   + ' | \n')
        f.write('MQM    ' + key_dict['MQM']  + ' | \n')
        f.write('M1     ' + key_dict['M1']   + ' | \n')
        f.write('E2     ' + key_dict['E2']   + ' | \n')
        f.write('E3     ' + key_dict['E3']   + ' | \n')
        f.write('M2     ' + key_dict['M2']   + ' | \n')
        f.write('M3     ' + key_dict['M3']   + ' | \n')
        f.write('========================= \n')
        f.write('Nhf = ' + str(Nso) + ' - SCF procedure includes Nhf upper shells (Nsh,Nsh-1,...) \n')
        f.write('Kmg =  0 - if not zero, Omega gives frequency of external field \n')
        f.write('Omega= 0.057580(a.u.) \n')
        f.write('Kex =  1 - key for exchange (0 - skip, 1 - include) \n')
        f.write('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n')

def write_dtm_in(mode, levels, operators):
    """ Write dtm.in """
    with open('dtm.in', 'w') as f:
        f.write('Mode = ' + mode + '\n')
        if levels:
            f.write('Levels = ' + levels + '\n')
        if operators:
            f.write('Operators = ' + operators)

def gen_key_list(matrix_elements):
    key_list = []
    if matrix_elements:
        if isinstance(matrix_elements, list):
            key_list = matrix_elements
        else:
            if len(matrix_elements.split(' ')) == 1:
                key_list = [matrix_elements]
            else:
                for matrix_element in matrix_elements.split(' '):
                    key_list.append(matrix_element.replace('[','').replace(']','').replace(',',''))
    return key_list

def read_nlv_from_conf(conf_path):
    """Read the number of levels (Nlv) from CONF.INP header"""
    conf_file = os.path.join(conf_path, 'CONF.INP')
    if not os.path.isfile(conf_file):
        return None
    try:
        with open(conf_file, 'r') as f:
            for line in f:
                if 'Nlv=' in line:
                    nlv_str = line.split('Nlv=')[1].strip().split()[0]
                    return int(nlv_str)
    except (IOError, ValueError, IndexError):
        return None
    return None

def copy_ci_files(src, dtm_dir, conf_prefix='CONF'):
    """Copy the standard set of CI files from src into dtm_dir."""
    for fname in ['CONF.INP', 'CONF.DET', 'CONF.XIJ', 'CONFSTR.RES', 'CONF.DAT', 'CONF.INT']:
        dest = fname.replace('CONF', conf_prefix, 1)
        run_shell('cp ' + src + '/' + fname + ' ' + dtm_dir + '/' + dest)

if __name__ == "__main__":
    yml_file = input("Input yml-file: ")
    config = read_yaml(yml_file)

    atom = get_dict_value(config, 'atom')
    code_method = get_dict_value(atom, 'code_method')
    if not isinstance(code_method, list): code_method = [code_method]

    system = get_dict_value(config, 'system')
    pci_version = get_dict_value(system, 'pci_version')
    on_hpc = get_dict_value(system, 'on_hpc')
    bin_dir = get_dict_value(system, 'bin_directory')
    run_codes = get_dict_value(system, 'run_codes')
    for_portal = get_dict_value(system, 'for_portal')

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

    basis = get_dict_value(config, 'basis')
    conf = get_dict_value(config, 'conf')

    if for_portal:
        print('Portal mode enabled for DTM')
    else:
        odd_J = get_dict_value(conf['odd'], 'J')
        even_J = get_dict_value(conf['even'], 'J')
        conf_odd_path = 'odd' + str(odd_J)[0]
        conf_even_path = 'even' + str(even_J)[0]

    dtm = get_dict_value(config, 'dtm')
    include_rpa = get_dict_value(dtm, 'include_rpa')

    # Create list of dtm directories
    dtm_dirs = []

    # DM parameters
    dm = get_dict_value(dtm, 'DM')
    dm_parity = get_dict_value(dm, 'parity') or 'both'
    dm_matrix_elements = get_dict_value(dm, 'matrix_elements')
    dm_range = get_dict_value(dm, 'level_range')
    include_dm = bool(dm_range)
    if include_dm:
        dm_odd = get_dict_value(dm_range, 'odd')
        dm_even = get_dict_value(dm_range, 'even')
        if dm_odd:
            from_level_odd = dm_odd.split(' ')[0]
            to_level_odd = dm_odd.split(' ')[1]
            dtm_dirs.append('dm_odd')
        if dm_even:
            from_level_even = dm_even.split(' ')[0]
            to_level_even = dm_even.split(' ')[1]
            dtm_dirs.append('dm_even')

    # TM parameters
    tm = get_dict_value(dtm, 'TM')
    tm_matrix_elements = get_dict_value(tm, 'matrix_elements')
    tm_from = get_dict_value(tm, 'from')
    tm_from_parity = get_dict_value(tm_from, 'parity')
    tm_from_range = get_dict_value(tm_from, 'level_range')
    tm_to = get_dict_value(tm, 'to')
    tm_to_parity = get_dict_value(tm_to, 'parity')
    tm_to_range = get_dict_value(tm_to, 'level_range')

    if for_portal:
        dtm_dirs.extend(list(TM_DIR_PATHS))
    else:
        include_tm = bool(tm_from_range and tm_to_range)
        if include_tm:
            from_level_initial = tm_from_range.split(' ')[0]
            from_level_final = tm_from_range.split(' ')[1]
            to_level_initial = tm_to_range.split(' ')[0]
            to_level_final = tm_to_range.split(' ')[1]
            dtm_dirs.append('tm')

    dm_key_list = gen_key_list(dm_matrix_elements)
    tm_key_list = gen_key_list(tm_matrix_elements)

    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'

    dir_path = os.getcwd()
    for method in code_method:
        full_path = dir_path + '/' + method if len(code_method) > 1 else dir_path
        os.chdir(full_path)

        if for_portal:
            current_dir = os.getcwd()
            all_ci_dirs = set(p for pair in TM_DIR_PATHS.values() for p in pair)
            dir_exists = {d: os.path.isfile(os.path.join(current_dir, d, 'CONF.RES')) or
                             os.path.isfile(os.path.join(current_dir, d, 'CONFFINAL.RES'))
                          for d in all_ci_dirs}
            all_ci_dirs_exist = all(dir_exists.values())
            if not all_ci_dirs_exist:
                missing = [d for d, ok in dir_exists.items() if not ok]
                print('Portal mode in ' + current_dir + ': missing CI directories: ' + ', '.join(missing))
                sys.exit()
        else:
            even_exists, odd_exists = False, False
            if tm_from_parity == 'even' or tm_to_parity == 'even' or dm_parity in ('even', 'both'):
                if os.path.isfile(conf_even_path + '/CONF.RES') or os.path.isfile(conf_even_path + '/CONFFINAL.RES'):
                    even_exists = True
            if tm_from_parity == 'odd' or tm_to_parity == 'odd' or dm_parity in ('odd', 'both'):
                if os.path.isfile(conf_odd_path + '/CONF.RES') or os.path.isfile(conf_odd_path + '/CONFFINAL.RES'):
                    odd_exists = True
            if not even_exists and not odd_exists:
                print('ci directories could not be found')
                sys.exit()

        # Create dtm directories with input files and job scripts
        for dtm_dir in dtm_dirs:
            Path(full_path + '/' + dtm_dir).mkdir(parents=True, exist_ok=True)
            if include_rpa:
                write_dtm_in('Init', '', '')
                if 'tm' in dtm_dir:
                    ops = TM_DIR_OPERATORS[dtm_dir] if dtm_dir in TM_DIR_OPERATORS else tm_key_list
                    write_mbpt_inp(basis, ops)
                elif 'dm' in dtm_dir:
                    write_mbpt_inp(basis, dm_key_list)
                run_shell('cp dtm.in ' + dtm_dir + '/dtm.in')
                run_shell('cp MBPT.INP ' + dtm_dir + '/MBPT.INP')
                run_shell('cp basis/HFD.DAT ' + dtm_dir + '/HFD.DAT')
                with open('rpa.in', 'w') as f:
                    f.write('2')
                run_shell('cp rpa.in ' + dtm_dir + '/rpa.in')
                script_name = write_job_script('.', 'dtm_rpa', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir) if on_hpc else None
                if script_name:
                    run_shell('cp dtm_rpa.qs ' + dtm_dir + '/dtm_rpa.qs')
            else:
                if dtm_dir == 'tm':
                    levels = from_level_initial + ' ' + from_level_final + ', ' + to_level_initial + ' ' + to_level_final
                    write_dtm_in('TM', levels, ', '.join(tm_key_list))
                elif dtm_dir in TM_DIR_PATHS:
                    from_key, to_key = TM_DIR_PATHS[dtm_dir]
                    from_nlv = read_nlv_from_conf(os.path.join(full_path, from_key))
                    to_nlv = read_nlv_from_conf(os.path.join(full_path, to_key))
                    if from_nlv and to_nlv:
                        levels = '1 ' + str(from_nlv) + ', 1 ' + str(to_nlv)
                        write_dtm_in('TM', levels, ', '.join(TM_DIR_OPERATORS[dtm_dir]))
                        print(dtm_dir + ': using levels 1-' + str(from_nlv) + ' (' + from_key + ') -> 1-' + str(to_nlv) + ' (' + to_key + ')')
                elif dtm_dir == 'dm_even':
                    write_dtm_in('DM', from_level_even + ' ' + to_level_even, ', '.join(dm_key_list))
                elif dtm_dir == 'dm_odd':
                    write_dtm_in('DM', from_level_odd + ' ' + to_level_odd, ', '.join(dm_key_list))
                run_shell('cp dtm.in ' + dtm_dir + '/dtm.in')
                script_name = write_job_script('.', 'dtm', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir) if on_hpc else None
                if script_name:
                    run_shell('cp dtm.qs ' + dtm_dir + '/dtm.qs')

            # Copy CI files to dtm directory
            if dtm_dir == 'tm':
                from_path = conf_even_path if tm_from_parity == 'even' else conf_odd_path
                to_path   = conf_even_path if tm_to_parity   == 'even' else conf_odd_path
                copy_ci_files(from_path, dtm_dir, 'CONF')
                run_shell('cp ' + to_path + '/CONF.INP '    + dtm_dir + '/CONF1.INP')
                run_shell('cp ' + to_path + '/CONF.DET '    + dtm_dir + '/CONF1.DET')
                run_shell('cp ' + to_path + '/CONF.XIJ '    + dtm_dir + '/CONF1.XIJ')
                run_shell('cp ' + to_path + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR1.RES')
            elif dtm_dir in TM_DIR_PATHS:
                from_key, to_key = TM_DIR_PATHS[dtm_dir]
                copy_ci_files(from_key, dtm_dir, 'CONF')
                run_shell('cp ' + to_key + '/CONF.INP '    + dtm_dir + '/CONF1.INP')
                run_shell('cp ' + to_key + '/CONF.DET '    + dtm_dir + '/CONF1.DET')
                run_shell('cp ' + to_key + '/CONF.XIJ '    + dtm_dir + '/CONF1.XIJ')
                run_shell('cp ' + to_key + '/CONFSTR.RES ' + dtm_dir + '/CONFSTR1.RES')
            elif dtm_dir == 'dm_even':
                copy_ci_files(conf_even_path, dtm_dir, 'CONF')
            elif dtm_dir == 'dm_odd':
                copy_ci_files(conf_odd_path, dtm_dir, 'CONF')

        # Submit job scripts
        if run_codes:
            for dtm_dir in dtm_dirs:
                dtm_path = full_path + '/' + dtm_dir
                os.chdir(dtm_path)
                if include_rpa:
                    if on_hpc:
                        run_shell('mpirun -n 1 ' + bin_dir + 'pdtm > dtm0.out')
                        if dtm_dir == 'tm':
                            levels = from_level_initial + ' ' + from_level_final + ', ' + to_level_initial + ' ' + to_level_final
                            write_dtm_in('TM', levels, ', '.join(tm_key_list))
                        elif dtm_dir in TM_DIR_PATHS:
                            from_key, to_key = TM_DIR_PATHS[dtm_dir]
                            from_nlv = read_nlv_from_conf(os.path.join(full_path, from_key))
                            to_nlv = read_nlv_from_conf(os.path.join(full_path, to_key))
                            if from_nlv and to_nlv:
                                levels = '1 ' + str(from_nlv) + ', 1 ' + str(to_nlv)
                                write_dtm_in('TM', levels, ', '.join(TM_DIR_OPERATORS[dtm_dir]))
                        elif dtm_dir == 'dm_even':
                            write_dtm_in('DM', from_level_even + ' ' + to_level_even, ', '.join(dm_key_list))
                        elif dtm_dir == 'dm_odd':
                            write_dtm_in('DM', from_level_odd + ' ' + to_level_odd, ', '.join(dm_key_list))
                        if submit_job:
                            run_shell('sbatch dtm_rpa.qs')
                else:
                    if submit_job:
                        run_shell('sbatch dtm.qs')
                os.chdir(full_path)
            os.chdir(dir_path)