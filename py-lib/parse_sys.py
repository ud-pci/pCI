import os
from datetime import datetime
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT, run
import re
import sys

'''
This script reads the output of sinfo and predicts the most efficient distribution of resources for a job

'''

def read_file_lines(filename):
    '''
    This function returns the lines of the file 'filename'
    '''
    f = open(filename, 'r') 
    lines = f.readlines()
    f.close()
    return lines

def count_num_idle_nodes():
    '''
    This function takes in a partition name and returns the number of nodes running 'idle'
    '''
    standard = 'r1n'
    large = 'r2l'
    xlarge = 'r2x'

    num_idle_nodes = 0
    running_idle_nodes = []

    run('squeue > s2.out', shell=True)
    lines = read_file_lines('s2.out')

    num_running_idle_standard, num_running_idle_large, num_running_idle_xlarge = 0, 0, 0
    running_idle_standard, running_idle_large, running_idle_xlarge = [], [], []

    # Obtain an array of nodelists containing respective partitions
    for line in lines:
        job_id = line.split()[0]
        partition = line.split()[1]
        jobname = line.split()[2]
        username = line.split()[3]
        state = line.split()[4]
        time = line.split()[5]
        num_nodes = line.split()[6]
        nodelist = line.split()[7]

        if partition == 'idle' and state == 'R':
            if standard in nodelist:
                running_idle_standard.append(nodelist)
            if large in nodelist:
                running_idle_large.append(nodelist)
            if xlarge in nodelist:
                running_idle_xlarge.append(nodelist)

    num_running_idle_standard = count_num_running_idle_nodes('r1n', running_idle_standard)
    num_running_idle_large = count_num_running_idle_nodes('r2l', running_idle_large)
    num_running_idle_xlarge = count_num_running_idle_nodes('r2x', running_idle_xlarge)

    return num_running_idle_standard, num_running_idle_large, num_running_idle_xlarge

def count_num_running_idle_nodes(partition, running_idle_partition):
    '''
    This function takes in a partition id and nodelist obtained from squeue and returns the number of nodes running on the partition in idle
    '''
    list_running_idle = []
    num_running_between = 0
    num_running_idle = 0

    for nodes in running_idle_partition:
        # Check two places after "r2l" for either node number or list of nodes in large partition
        if partition in nodes:
            # count number of nodes found
            node_split = re.findall(r'\d+', nodes.split(partition,1)[1])
            for node in node_split:
                list_running_idle.append(node)
            # count number of nodes "between" in nodelist
            if re.search(r'-', nodes):
                num_running_between += int(nodes[7:9]) - int(nodes[4:6]) - 1
    
    list_running_idle = [*set(list_running_idle)]
    num_running_idle = len(list_running_idle) + num_running_between
    return num_running_idle

if __name__ == "__main__":
    dir_path = os.getcwd()
    run('sinfo > s.out', shell=True)
    lines = read_file_lines('s.out')
    num_free_nodes = {'standard' : 0, 'large-mem' : 0, 'xlarge-mem' : 0}
    num_idle_nodes = {'standard' : 0, 'large-mem' : 0, 'xlarge-mem' : 0}
    partition = ['standard', 'large-mem', 'xlarge-mem']
    for line in lines:
        if line.split()[0].replace('*','') in partition and line.split()[4] == 'idle':
            num_free_nodes[line.split()[0].replace('*','')] = line.split()[3]
    
    current_time = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    print('TIME: ' + current_time)
    print('FREE NODES:')
    for item in num_free_nodes:
        print(item, ':', num_free_nodes[item])

    num_idle_standard, num_idle_large, num_idle_xlarge = count_num_idle_nodes()

    print('==================')
    print('IDLE NODES:')
    print('standard:', num_idle_standard)
    print('large-mem:', num_idle_large)
    print('xlarge-mem:', num_idle_xlarge)

    print('==================')
    print('TOTAL FREE NODES:')
    print('standard:', num_idle_standard + int(num_free_nodes['standard']))
    print('large-mem:', num_idle_large + int(num_free_nodes['large-mem']))
    print('xlarge-mem:', num_idle_xlarge + int(num_free_nodes['xlarge-mem']))

    run('rm s.out s2.out', shell=True)
    #num_req_nodes = input('How many nodes would you like to use? ')
    #req_mem = input('How much memory is required per core? ')