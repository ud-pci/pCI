import re
import math
import sys
import os
from subprocess import Popen, PIPE, STDOUT, run
from datetime import datetime

'''
This python script upscales ci calculations in the following way:

Example:
17g -> 22g -> 26g -> 30g
 |      |      |      |
17h -> 22h -> 26h -> 30h
 |      |      |      | 
17i -> 22i -> 26i -> 30i

This script will generate a job script requesting a specified amount of resources
then run CI calculations starting from the top left to the bottom right. 
All calculations are done in their own designated directories.
After every CI calculation, con_cut is run to strip off configurations with weights lower than 10^-9
17g -> 22g > 17h -> 22h
26g + 

This is done in 3 main parts:
1. first run (upper-left block, e.g. 17g)
2. second run (right of first run, e.g. 22g) 
3. subsequent runs (bottom of second run, but requires bottom of first run as well, e.g. 22h requires 22g and 17h)
'''

def run_executable(code, nl):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    # This function 
    if code == 'add_nr':
        # add_nr 
        run('./add_nr < add.in > add' + nl + '.out', shell=True)
    elif code == 'con_cut':
        run('./con_cut < ../list/cut9.in > con_cut' + nl + '.out', shell=True)
    elif code == 'concmp':
        # concmp compares two configuration lists C_A.INP and C_B.INP
        # and writes a list of configurations only in C_B.INP into file C_B_new.INP
        run('time ./concmp < c3.in > concmp' + nl + '.out', shell=True)
    elif code == 'merge_ci':
        run('./merge_ci > merge_ci' + nl + '.out', shell=True)
    elif code == 'ci':
        run('sh ci.qs', shell=True)
    else:
        print('code not implemented')
        sys.exit()
    
    print(code + ' completed at' + current_time)

def mkdir(dir_name):
    # This function makes a directory specified by dir_name and throws an exception if it already exists
    try:
        os.mkdir(dir_name)
    except OSError as e:
        print(e)

def cp(path_from, path_to):
    # This function imitates the shell command cp to copy and paste a file from one directory to another

    # Check to see if copying multiple files or not by checking whether '*' is specified or if path leads to a directory
    isdir = os.path.isdir(path_from) or path_from[-1] == '*'
    
    # If isdir is true, include the -r recursive flag
    if isdir == True:
        run('cp -r ' + path_from + ' ' + path_to, shell=True)
    else:
        run('cp ' + path_from + ' ' + path_to, shell=True)

def initial_run(n, l):
    # This function runs the first CI calculation for nl
    nl = str(n) + str(l)
    add_inp_path = '../list/E/'
    add_nl_inp = 'ADD' + nl + '.INP'
    add_nl_inp_path = add_inp_path + add_nl_inp

    # Make a directory for nl and throw an exception if it already exists
    mkdir(nl)

    # Move into nl directory
    os.chdir(nl)

    # Copy files from list directory to working directory
    cp('../list/*', '.')

    # Copy relevant ADD.INP to working directory
    cp(add_nl_inp_path, '.')

    # Run add_nr to construct CONF.INP
    cp(add_nl_inp, 'ADD.INP')
    run_executable('add_nr', nl)

    # Run CI codes
    run_executable('ci', nl)

    # Run con_cut to strip all configurations with weights lower than 10^-9
    run_executable('con_cut', nl)

    # Return to root directory
    os.chdir('../')


def upscale_n_or_l(n_or_l, n0, l0):
    # This function upscales the n0l0 calculation to nl
    # e.g. n0l0 = 17g, nl = 22g; or n0l0 = 17g, nl = 17h
    n0l0 = str(n0) + str(l0)

    # Determine if n0 or l0 is being upscaled
    # if n_or_l is a number, n0 is being upscaled; e.g. 17g -> 22g
    # if n_or_l is a letter, s < p < d < f < g < h < i < k ...
    if isinstance(n_or_l, int):
        # check that upscaled n is larger than n0
        if int(n_or_l) > int(n0):
            nl = str(n_or_l) + str(l0)
        else:
            print(n0l0 + ' cannot be upscaled to ' + nl)
            sys.exit()
    elif n_or_l.isalpha():
        # TODO - check upscaling h > g > f > d > p > s
        nl = str(n0) + str(n_or_l)

    add_inp_path = '../list/E/'
    add_n0l0_inp = 'ADD' + n0l0 + '.INP'
    add_nl_inp = 'ADD' + nl + '.INP'
    add_n0l0_inp_path = add_inp_path + add_n0l0_inp
    add_nl_inp_path = add_inp_path + add_nl_inp

    # Make a directory for nl and throw an exception if it already exists
    mkdir(nl)

    # Move into nl directory
    os.chdir(nl)

    # Copy files from list directory to working directory
    cp('../list/*', '.')

    # Copy relevant ADD.INP to working directory
    cp(add_n0l0_inp_path, '.')
    cp(add_nl_inp_path, '.')
    
    # Run add_nr for n0l0 to construct C_A.INP 
    cp(add_n0l0_inp, 'ADD.INP')
    run_executable('add_nr', nl)
    cp('CONF.INP', 'C_A.INP')

    # Run add_nr for nl to construct C_B.INP
    cp(add_nl_inp, 'ADD.INP')
    run_executable('add_nr', nl)
    cp('CONF.INP', 'C_B.INP')

    # Run concmp to produce list of configurations that are in C_B.INP but not in C_A.INP
    run_executable('concmp', nl)

    # Merge configurations only in C_B.INP with select weighted configurations from n0l0 list
    con_cut_n0l0_path = '../' + n0l0 + '/CON_CUT.RES'
    cp(con_cut_n0l0_path, 'C_A.INP')
    cp('C_B_new.INP', 'CONF.INP')
    run_executable('merge_ci', nl)
    cp('C_M.INP', 'CONF.INP')

    # Run CI codes
    run_executable('ci', nl)

    # Run con_cut to strip all configurations with weights lower than 10^-9
    run_executable('con_cut', nl)

    # Return to root directory
    os.chdir('../')

def upscale_nl(n, l, n0, l0):
    # This function upscales to nl by merging previous n0l and nl0 runs
    # e.g. nl0 = 22g, n0l = 17h -> nl = 22h
    nl = str(n) + str(l) 
    nl0 = str(n) + str(l0)
    n0l = str(n0) + str(l)

    add_inp_path = '../list/E/'
    add_nl0_inp = 'ADD' + nl0 + '.INP'
    add_n0l_inp = 'ADD' + n0l + '.INP'
    add_nl_inp = 'ADD' + nl + '.INP'
    add_nl0_inp_path = add_inp_path + add_nl0_inp
    add_n0l_inp_path = add_inp_path + add_n0l_inp
    add_nl_inp_path = add_inp_path + add_nl_inp

    # Make a directory for nl and throw an exception if it already exists
    mkdir(nl)

    # Move into nl directory
    os.chdir(nl)

    # Copy files from list directory to working directory
    cp('../list/*', '.')

    # Copy relevant ADD.INP to working directory
    cp(add_nl0_inp_path, '.')
    cp(add_n0l_inp_path, '.')
    cp(add_nl_inp_path, '.')

    # Subtract n0l configurations from nl configurations, e.g. 17h from 22h
    cp(add_n0l_inp, 'ADD.INP')
    run_executable('add_nr', n0l)
    cp('CONF.INP', 'C_A.INP')

    cp(add_nl_inp, 'ADD.INP')
    run_executable('add_nr', nl)
    cp('CONF.INP', 'C_B.INP')

    # Obtain list of configurations unique to nl between nl and n0l
    run_executable('concmp', n0l)

    # Subtract nl0 configurations from nl configurations, e.g. 22g from 22h
    cp(add_nl0_inp, 'ADD.INP')
    run_executable('add_nr', nl0)
    cp('CONF.INP', 'C_A.INP')

    cp('C_B_new.INP', 'C_B.INP')

    # Obtain list of configurations unique to nl between nl and n1l1
    run_executable('concmp', nl0) 

    # Merge con_cut -9 17h and 22g runs
    con_cut_nl0_path = '../' + nl0 + '/CON_CUT.RES'  # 22g
    con_cut_n0l_path = '../' + n0l + '/CON_CUT.RES'  # 17h

    cp(con_cut_nl0_path, 'CONF.INP') # 22g
    cp(con_cut_n0l_path, 'C_A.INP')  # 17h

    run_executable('merge_ci',nl0 + '+' + n0l)
    cp('C_M.INP', 'C_A.INP')

    # Merge with new configurations
    cp('C_B_new.INP', 'CONF.INP')
    run_executable('merge_ci', nl)
    cp('C_M.INP', 'CONF.INP')

    # Run CI codes
    run_executable('ci', nl)

    # Run con_cut to strip all configurations with weights lower than 10^-9
    run_executable('con_cut', nl)

    # Return to root directory
    os.chdir('../')

def read_env():
    # This function reads slurm environment variables
    nnodes = os.environ['SLURM_NNODES']
    ntasks = os.environ['SLURM_NTASKS']
    ntasks_per_node = os.environ['SLURM_TASKS_PER_NODE']
    mem_per_cpu = os.environ['SLURM_MEM_PER_CPU']
    partition = os.environ['SLURM_JOB_PARTITION']
    
    return nnodes, ntasks, partition
    
def write_job_script(filename, code, num_nodes, num_procs_per_node, exclusive, mem, partition):
    if code == 'conf':
        with open('c.in', 'w') as f:
            f.write('2, 2, 0, 0, 1')
        f.close()

    with open(filename, 'w') as f:
        f.write('#!/bin/bash -l \n')
        f.write(' \n')
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
        f.write('vpkg_require pci \n')
        f.write(' \n')
        f.write('UD_PREFER_MEM_PER_CPU=YES \n')
        f.write('UD_REQUIRE_MEM_PER_CPU=YES \n')
        f.write(' \n')
        f.write('. /opt/shared/slurm/templates/libexec/openmpi.sh \n')
        f.write(' \n')
        f.write('ulimit -s unlimited \n')
        f.write('CONF_MAX_BYTES_PER_CPU=$((SLURM_MEM_PER_CPU*1024*1024)) \n')
        f.write('export CONF_MAX_BYTES_PER_CPU \n')
        if code == 'basc' or code == 'ci':
            f.write('${UD_MPIRUN} basc > basc.out \n')
        if code == 'conf' or code == 'ci':
            f.write('${UD_MPIRUN} conf > conf.out \n')

        f.write('mpi_rc=$? \n')
        f.write(' \n')
        f.write('exit $mpi_rc \n')
    f.close()

if __name__ == "__main__":
    '''
    arguments
    n_array - array of principal quantum numbers
    l_array - array of angular momentum quantum numbers

    example:
    n = [17, 22, 26, 30]
    l = [g, h, i]

    loop n then l:
        17g, 22g, 26g, 30g
        17h, 22h, 26h, 30h
        17i, 22i, 26i, 30i

    function upscale_n_or_l() upscales in one direction to the right or downward (higher n or l)
    function upscale_nl() upscales diagonally (higher n and l) by merging previous off-diagonal calculation

    '''
    '''
    Required bin files in list directory:
    add_nr
    con_cut
    '''

    print('upscaling procedure starting at ' + current_time)

    # Read slurm environment variables to generate a job script for CI with identical parameters
    nnodes, ntasks, partition = read_env()
    write_job_script('ci.qs', 'ci', nnodes, ntasks, False, 0, partition)

    # Move job script to master list directory
    cp('c.in', 'list/c.in')
    cp('ci.qs', 'list/ci.qs')
    
    # test arrays
    n = [17, 22, 26, 30]
    l = ['g', 'h', 'i']

    # Remove test directories
    #run('rm -r 17g', shell=True)
    #run('rm -r 17h', shell=True)
    #run('rm -r 22g', shell=True)
    #run('rm -r 22h', shell=True)

    # Find current working directory
    dir_path = os.getcwd()

    '''
    # Run initial run with n[0], l[0] - e.g. 17g
    initial_run(n[0],l[0])

    # Upscale n0l0 to nl0, e.g. 17g -> 22g
    upscale_n_or_l(n[1],n[0],l[0])

    # Upscale n0l0 to n0l, e.g. 17g -> 17h
    upscale_n_or_l(l[1],n[0],l[0])
    '''
    # Upscale n0l0 to nl by merging n0l and nl0 runs, e.g. 22g + 17h -> 22h
    upscale_nl(n[1],l[1],n[0],l[0])

    # Upscale n1l0 to n2l0, e.g. 22g -> 26g
    upscale_n_or_l(n[2],n[1],l[0])

    # 26h with 26g + 22h
    upscale_nl(n[2],l[1],n[1],l[0])
     
    print('upscaling')