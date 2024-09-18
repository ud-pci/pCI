import sys
from subprocess import run
import os

def find_cluster():
    '''
    This function finds the name of the cluster you are working on - caviness or darwin
    '''
    cluster = ''
    
    # Fetch slurm node and core information
    run('sinfo -o "%n %e %m %a %c %C" > sinfo.out', shell=True)
    
    f = open('sinfo.out', 'r')
    lines = f.readlines()
    f.close()

    if lines[0].split()[0] == 'HOSTNAMES':
        hostname = lines[1].split()[0]
    else:
        hostname = ''
    
    # Current logic for determining cluster:
    # Caviness hostname has length 6: rxxnxx
    # Darwin hostname has length 5: rxnxx
    if len(hostname) == 5:
        cluster = 'darwin'
    elif len(hostname) == 6:
        cluster = 'caviness'
    else:
        cluster = 'cluster could not be identified'
    
    run('rm sinfo.out', shell=True)
    
    return cluster

def write_job_script(path, code, num_nodes, num_procs_per_node, exclusive, mem, partition, pci_version, bin_dir):
    
    # Determine cluster information
    cluster = find_cluster()
    
    # Ensure partition is in cluster
    darwin_partitions = ['standard', 'large-mem', 'xlarge-mem', 'idle']
    if partition in darwin_partitions and cluster == 'darwin':
        print('partition ' + partition + ' found on darwin')
    else:
        print('partition ' + partition + ' was not found on darwin')
        partition = 'standard'
        print('partition has been defaulted to the standard partition')
    
    # Specify directory of executables
    if bin_dir and bin_dir[-1] != '/':
        bin_dir = bin_dir + '/'
    
    # Default filenames
    filenames = {'ci+all-order': 'ao.qs', 
                 'all-order': 'ao.qs', 
                 'ci+second-order': 'cis.qs',
                 'second-order': 'cis.qs',
                 'conf': 'conf.qs',
                 'basc': 'basc.qs',
                 'ci': 'ci.qs',
                 'dtm': 'dtm.qs',
                 'dtm_rpa': 'dtm_rpa.qs'}
    
    # Serial codes
    is_serial = {'ci+all-order': True, 
                 'all-order': True, 
                 'ci+second-order': True,
                 'second-order': True,
                 'conf': False,
                 'basc': False,
                 'ci': False,
                 'dtm': False,
                 'dtm_rpa': False,}
    
    os.chdir(path)
    
    # Write job script
    filename = filenames[code]
    with open(filename, 'w') as f:
        f.write('#!/bin/bash -l \n')
        f.write(' \n')
        if is_serial[code]:
            f.write('#SBATCH --ntasks=1 \n')
        else:
            f.write('#SBATCH --nodes=' + str(num_nodes) + ' \n')
            f.write('#SBATCH --tasks-per-node=' + str(num_procs_per_node) + ' \n')
        if exclusive: 
            f.write('#SBATCH --exclusive=user \n')
        f.write('#SBATCH --cpus-per-task=1 \n')
        f.write('#SBATCH --mem=' + str(mem) + ' \n')
        f.write('#SBATCH --job-name=' + code + ' \n')
        f.write('#SBATCH --partition=' + partition + ' \n')
        f.write('#SBATCH --time=05-00:00:00 \n')
        f.write('#SBATCH --export=NONE \n')
        f.write(' \n')
        if pci_version != 'default':
            f.write('vpkg_require pci' + '/' + pci_version + ' \n')
        else:
            f.write('vpkg_require pci \n')
        if cluster == 'darwin':
            f.write(' \n')
            f.write('UD_PREFER_MEM_PER_CPU=YES \n')
            f.write('UD_REQUIRE_MEM_PER_CPU=YES \n')
        f.write(' \n')
        f.write('. /opt/shared/slurm/templates/libexec/openmpi.sh \n')
        f.write(' \n')
        
        if not is_serial[code]:
            f.write('CONF_MAX_BYTES_PER_CPU=$((SLURM_MEM_PER_CPU*1024*1024)) \n')
            f.write('export CONF_MAX_BYTES_PER_CPU \n')
        
        # Executables
        if code == 'ci':
            f.write('${UD_MPIRUN} ' + bin_dir + 'pbasc \n')
            f.write('${UD_MPIRUN} ' + bin_dir + 'pconf \n')
        elif code == 'dtm':
            f.write('${UD_MPIRUN} ' + bin_dir + 'pdtm \n')
        elif code == 'dtm_rpa':
            f.write(bin_dir + 'rpa < rpa.in \n')
            f.write(bin_dir + 'rpa_dtm \n')
            f.write('${UD_MPIRUN} ' + bin_dir + 'pdtm \n')
        elif code == 'all-order' or code == 'ci+all-order':
            f.write('time allcore-rle-ci <inf.aov >out.core \n')
            f.write('time valsd-rle-cis <inf.aov >out.val \n')
            f.write('time sdvw-rle-cis <inf.aov >out.vw \n')
            f.write('time second-cis <inf.vw >out.second.vw \n')
        elif code == 'second-order' or code == 'ci+second-order':
            f.write('time second-cis <inf.vw >out.second.vw \n')
        else:
            print(code + ' is not supported')
            sys.exit()
            
        f.write(' \n')
        f.write('mpi_rc=$? \n')
        f.write('exit $mpi_rc \n')
    f.close()  
    print(filename + ' has been generated')
    
    return filename    
    
if __name__ == '__main__':
    code = 'ci+all-order'
    num_nodes = 5
    num_procs_per_node = 64
    exclusive = True
    mem = 0
    partition = 'safrono'
    
    write_job_script(code, num_nodes, num_procs_per_node, exclusive, mem, partition)