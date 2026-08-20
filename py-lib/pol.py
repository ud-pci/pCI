""" ine

This script sets up and submits a polarizability calculation using the pine
program, driven by parameters in a config.yml file.

The config.yml file should have the following blocks:

    * system - system parameters (bin_directory, run_codes, on_hpc)
    * atom.code_method - list of methods (CI, CI+all-order, CI+MBPT)
    * ine - parameters for pine (parity, N0, N2, rhs, lhs, field_type, wavelength_range, step_size,
            and optionally tm_dir to override the TM directory name)

This script has 2 main capabilities:
1. Set up the directory for a polarizability calculation.
2. Submit a job script to run pine on HPC.

"""
import yaml
import os
import sys
from pathlib import Path
from utils import run_shell, get_dict_value, convert_params_to_list, check_slurm_installed
from gen_job_script import write_job_script


def read_yaml(filename):
    with open(filename, 'r') as f:
        config = yaml.safe_load(f)
    return config


def write_ine_in(rhs, lhs, n0, n2, include_static, include_dynamic, wavelength_range, step_size):
    """ Write ine.in for pine. """
    ranges_parts = []
    if include_static:
        ranges_parts.append('0 0 0')
    if include_dynamic:
        for i in range(len(wavelength_range)):
            ranges_parts.append(f'{wavelength_range[i]} {step_size[i]}')
    with open('ine.in', 'w') as f:
        f.write(f'RHS={rhs}\n')
        f.write(f'LHS={lhs}\n')
        f.write(f'N0={n0}\n')
        f.write(f'N2={n2}\n')
        f.write(f'Ranges= {", ".join(ranges_parts)}\n')


if __name__ == "__main__":
    yml_file = input('Input yml-file: ')
    config = read_yaml(yml_file)

    # system parameters
    system = get_dict_value(config, 'system')
    pci_version = get_dict_value(system, 'pci_version')
    on_hpc = get_dict_value(system, 'on_hpc')
    bin_dir = get_dict_value(system, 'bin_directory')

    atom = get_dict_value(config, 'atom')
    code_method = get_dict_value(atom, 'code_method')
    if not isinstance(code_method, list):
        code_method = [code_method]

    # hpc parameters
    on_slurm = check_slurm_installed()
    if on_hpc and on_slurm:
        hpc = get_dict_value(config, 'hpc')
        submit_job = get_dict_value(hpc, 'submit_job')
        if hpc:
            partition = get_dict_value(hpc, 'partition')
            nodes = get_dict_value(hpc, 'nodes')
            tasks_per_node = get_dict_value(hpc, 'tasks_per_node')
        else:
            print(f'hpc block was not found in {yml_file}')
            partition, nodes, tasks_per_node = None, 1, 1
    else:
        on_hpc = False
        submit_job = False

    for_portal = get_dict_value(system, 'for_portal')

    conf = get_dict_value(config, 'conf')
    odd_J = get_dict_value(conf['odd'], 'J')
    even_J = get_dict_value(conf['even'], 'J')
    conf_odd_path = f'odd{str(odd_J)[0]}'
    conf_even_path = f'even{str(even_J)[0]}'

    # pine parameters
    ine = get_dict_value(config, 'ine')
    parity = get_dict_value(ine, 'parity')
    n0 = get_dict_value(ine, 'N0')
    n2 = get_dict_value(ine, 'N2')
    rhs = get_dict_value(ine, 'rhs')
    lhs = get_dict_value(ine, 'lhs')
    field_type = convert_params_to_list(get_dict_value(ine, 'field_type'))

    # derive tm_dir from config, or compute a default
    tm_dir = get_dict_value(ine, 'tm_dir')
    if not tm_dir:
        if not for_portal:
            tm_dir = 'tm'
        else:
            # opposite-parity operators (E1, H_pnc) connect even<->odd
            same_parity_rhs = [5]  # E2
            if rhs in same_parity_rhs:
                raise ValueError(f'Cannot derive tm_dir for same-parity operator rhs={rhs}. Set tm_dir in the ine config block.')
            other_parity = 'odd' if parity == 'even' else 'even'
            J_init = int(even_J) if parity == 'even' else int(odd_J)
            J_other = int(odd_J) if parity == 'even' else int(even_J)
            tm_dir = f'tm_{parity}{J_init}_{other_parity}{J_other}'
    include_static = 'static' in field_type
    include_dynamic = 'dynamic' in field_type

    if include_dynamic:
        wavelength_range = convert_params_to_list(get_dict_value(ine, 'wavelength_range'))
        step_size = convert_params_to_list(get_dict_value(ine, 'step_size'))
    else:
        wavelength_range, step_size = [], []

    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'

    # each code_method gets its own directory
    dir_path = os.getcwd()
    for method in code_method:
        full_path = f'{dir_path}/{method}' if len(code_method) > 1 else dir_path
        Path(full_path).mkdir(parents=True, exist_ok=True)
        os.chdir(full_path)

        ine_path = f'{full_path}/pol_{parity}{n0}'
        Path(f'pol_{parity}{n0}').mkdir(parents=True, exist_ok=True)
        write_ine_in(rhs, lhs, n0, n2, include_static, include_dynamic, wavelength_range, step_size)
        run_shell(f'mv ine.in {ine_path}/ine.in')

        script_name = write_job_script('.', 'ine', nodes, tasks_per_node, True, 0, partition, pci_version, bin_dir)
        if script_name:
            run_shell(f'mv ine.qs {ine_path}/ine.qs')
        else:
            print('job script was not generated. check parameters and generate manually.')

        # copy CI files from even and odd directories
        even_exists, odd_exists = False, False
        if os.path.isfile(f'{conf_even_path}/CONF.RES') or os.path.isfile(f'{conf_even_path}/CONFFINAL.RES'):
            even_exists = True
            if parity == 'even':
                run_shell(f'cp {conf_even_path}/CONF.DET {ine_path}/CONF0.DET')
                run_shell(f'cp {conf_even_path}/CONF.XIJ {ine_path}/CONF0.XIJ')
                run_shell(f'cp {conf_even_path}/FINAL.RES {ine_path}/FINAL.RES')
            elif parity == 'odd':
                run_shell(f'cp {conf_even_path}/CONF.DET {ine_path}/CONF.DET')
                run_shell(f'cp {conf_even_path}/CONF.XIJ {ine_path}/CONF.XIJ')
                run_shell(f'cp {conf_even_path}/CONF.INP {ine_path}/CONF.INP')
                run_shell(f'cp {conf_even_path}/CONF.DAT {ine_path}/CONF.DAT')
                run_shell(f'cp {conf_even_path}/pCONF.HIJ {ine_path}/pCONF.HIJ')
                run_shell(f'cp {conf_even_path}/pCONF.JJJ {ine_path}/pCONF.JJJ')

        if os.path.isfile(f'{conf_odd_path}/CONF.RES') or os.path.isfile(f'{conf_odd_path}/CONFFINAL.RES'):
            odd_exists = True
            if parity == 'odd':
                run_shell(f'cp {conf_odd_path}/CONF.DET {ine_path}/CONF0.DET')
                run_shell(f'cp {conf_odd_path}/CONF.XIJ {ine_path}/CONF0.XIJ')
                run_shell(f'cp {conf_odd_path}/FINAL.RES {ine_path}/FINAL.RES')
            elif parity == 'even':
                run_shell(f'cp {conf_odd_path}/CONF.DET {ine_path}/CONF.DET')
                run_shell(f'cp {conf_odd_path}/CONF.XIJ {ine_path}/CONF.XIJ')
                run_shell(f'cp {conf_odd_path}/CONF.INP {ine_path}/CONF.INP')
                run_shell(f'cp {conf_odd_path}/CONF.DAT {ine_path}/CONF.DAT')
                run_shell(f'cp {conf_odd_path}/pCONF.HIJ {ine_path}/pCONF.HIJ')
                run_shell(f'cp {conf_odd_path}/pCONF.JJJ {ine_path}/pCONF.JJJ')

        if not even_exists and not odd_exists:
            print('CI directories could not be found')
            sys.exit()

        if os.path.isfile(f'{tm_dir}/DTM_RPA.INT'):
            print(f'copying DTM_RPA.INT from {tm_dir}')
            run_shell(f'cp {tm_dir}/DTM_RPA.INT {ine_path}/DTM.INT')
        else:
            print(f'copying DTM.INT from {tm_dir}')
            run_shell(f'cp {tm_dir}/DTM.INT {ine_path}/DTM.INT')

        os.chdir(ine_path)
        if submit_job:
            run_shell('sbatch ine.qs')
