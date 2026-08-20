import sys
from subprocess import run, CalledProcessError
import os

def write_job_script(path, code, num_nodes, num_procs_per_node, exclusive, mem, partition, pci_version, bin_dir):
    """
    This function writes a SLURM job script and returns the name of the job script.
    """

    # Determine cluster information from 'scontrol show config'
    try:
        scontrol_output = run(['scontrol', 'show', 'config'], capture_output=True, text=True, check=True)
        for line in scontrol_output.stdout.splitlines():
            if line.startswith('ClusterName'):
                cluster = line.split('=')[1].strip()
                print(f'Cluster found: {cluster}')
    except (CalledProcessError, FileNotFoundError):
        print('SLURM_CLUSTER_NAME environment variable not found. Job script will not be written.')
        return None

    # Obtain list of partitions
    try:
        sinfo_output = run(['sinfo', '--format=%P'], capture_output=True, text=True, check=True)
        partitions = set(line.rstrip('*') for line in sinfo_output.stdout.strip().splitlines())
        if partition in partitions:
            print(f'partition {partition} found on {cluster}')
        else:
            print(f'partition {partition} not found on {cluster}')
            return None
    except CalledProcessError as e:
        print('Failed to retrieve partitions.')
        print(e)
        return None

    # Set default job parameters
    if not num_nodes:
        num_nodes = 1
    if not num_procs_per_node:
        num_procs_per_node = 1
    if not partition:
        if cluster == 'caviness':
            partition = 'standard'
        elif cluster == 'darwin':
            partition = 'idle'
        else:
            print('partition name was not specified')
            return None

    if bin_dir and bin_dir[-1] != '/':
        bin_dir += '/'

    filenames = {'ci+all-order': 'ao.qs',
                 'all-order': 'ao.qs',
                 'ci+second-order': 'cis.qs',
                 'second-order': 'cis.qs',
                 'conf': 'conf.qs',
                 'basc': 'basc.qs',
                 'ci': 'ci.qs',
                 'dtm': 'dtm.qs',
                 'dtm_rpa': 'dtm_rpa.qs',
                 'ine': 'ine.qs'}

    is_serial = {'ci+all-order': True,
                 'all-order': True,
                 'ci+second-order': True,
                 'second-order': True,
                 'conf': False,
                 'basc': False,
                 'ci': False,
                 'dtm': False,
                 'dtm_rpa': False,
                 'ine': False}

    os.chdir(path)

    filename = filenames[code]
    with open(filename, 'w') as f:
        f.write('#!/bin/bash -l\n')
        f.write('\n')
        if is_serial[code]:
            f.write('#SBATCH --ntasks=1\n')
        else:
            f.write(f'#SBATCH --nodes={num_nodes}\n')
            if num_procs_per_node == 1:
                print('tasks-per-node is set to minimum of 2 for parallel programs')
                num_procs_per_node = 2
            f.write(f'#SBATCH --tasks-per-node={num_procs_per_node}\n')
        if exclusive:
            f.write('#SBATCH --exclusive=user\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write(f'#SBATCH --mem={mem}\n')
        f.write(f'#SBATCH --job-name={code}\n')
        f.write(f'#SBATCH --partition={partition}\n')
        f.write('#SBATCH --time=05-00:00:00\n')
        f.write('#SBATCH --export=NONE\n')
        f.write('\n')
        if pci_version != 'default':
            f.write(f'vpkg_require pci/{pci_version}\n')
        else:
            f.write('vpkg_require pci\n')
        if cluster == 'darwin':
            f.write('\n')
            f.write('UD_PREFER_MEM_PER_CPU=YES\n')
            f.write('UD_REQUIRE_MEM_PER_CPU=YES\n')
        f.write('\n')
        f.write('. /opt/shared/slurm/templates/libexec/openmpi.sh\n')
        f.write('\n')

        if not is_serial[code]:
            f.write('CONF_MAX_BYTES_PER_CPU=$((SLURM_MEM_PER_CPU*1024*1024))\n')
            f.write('export CONF_MAX_BYTES_PER_CPU\n\n')

        # Executables
        if code == 'ci':
            f.write(f'${{UD_MPIRUN}} {bin_dir}pbasc\n')
            f.write('if [ $? -eq 0 ]; then\n')
            f.write(f'    ${{UD_MPIRUN}} {bin_dir}pconf\n')
            f.write('else\n')
            f.write('    echo "pbasc failed, so pconf will not run."\n    exit 1\nfi\n')
        elif code == 'dtm':
            f.write(f'${{UD_MPIRUN}} {bin_dir}pdtm\n')
        elif code == 'dtm_rpa':
            f.write(f'{bin_dir}rpa < rpa.in\n')
            f.write(f'{bin_dir}rpa_dtm\n')
            f.write(f'${{UD_MPIRUN}} {bin_dir}pdtm\n')
        elif code == 'ine':
            f.write(f'${{UD_MPIRUN}} {bin_dir}pine\n')
        elif code == 'all-order' or code == 'ci+all-order':
            f.write(f'{bin_dir}allcore-ci <inf.aov >out.core\n')
            f.write(f'{bin_dir}valsd-ci <inf.aov >out.val\n')
            f.write(f'{bin_dir}sdvw-ci <inf.aov >out.vw\n')
            f.write(f'{bin_dir}second-ci <inf.vw >out.second.vw\n')
        elif code == 'second-order' or code == 'ci+second-order':
            f.write(f'{bin_dir}second-ci <inf.vw >out.second.vw\n')
        else:
            print(f'{code} is not supported')
            sys.exit()

        f.write('\n')
        f.write('mpi_rc=$?\n')
        f.write('exit $mpi_rc\n')

    print(f'{filename} has been generated')
    return filename


if __name__ == '__main__':
    output_path = '.'
    code = 'ci'
    num_nodes = 5
    num_procs_per_node = 64
    exclusive = True
    mem = 0
    partition = 'large-mem'
    pci_version = 'default'
    bin_dir = '/lustre/safrono/users/1813/pCI-dev/bin'

    print(write_job_script(output_path, code, num_nodes, num_procs_per_node, exclusive, mem, partition, pci_version, bin_dir))
