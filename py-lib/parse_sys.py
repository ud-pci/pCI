from datetime import datetime
from subprocess import run
import re
import sys

'''
This script parses the SLURM scheduling queue and prints the number of available nodes for each partition 

'''

def read_file_lines(filename):
    '''
    This function returns the lines of the file 'filename'
    '''
    f = open(filename, 'r') 
    lines = f.readlines()
    f.close()
    return lines

def get_num_available_nodes():
    '''
    This function takes in a partition id and nodelist obtained from squeue and returns the number of nodes running on the partition
    '''
    
    # Initialize arrays of partition nodes
    standard, large, xlarge = {}, {}, {}

    [standard.setdefault(str(i).zfill(2), False) for i in range(48)]
    [large.setdefault(str(i).zfill(2), False) for i in range(32)]
    [xlarge.setdefault(str(i).zfill(2), False) for i in range(11)]

    # Get running job information from squeue
    run('squeue > s.out', shell=True)
    lines = read_file_lines('s.out')

    # Parse lines from squeue
    for line in lines:
        job_id = line.split()[0]
        partition = line.split()[1]
        jobname = line.split()[2]
        username = line.split()[3]
        state = line.split()[4]
        time = line.split()[5]
        num_nodes = line.split()[6]
        nodelist = line.split()[7]
        if partition == 'standard':
            standard = set_true(nodelist, standard)
        elif partition == 'large-mem':
            large = set_true(nodelist, large)
        elif partition == 'xlarge-me':
            xlarge = set_true(nodelist, xlarge)
        else:
            pass

    num_standard = sum(1 for v in standard.values() if v == False)
    num_large = sum(1 for v in large.values() if v == False)
    num_xlarge = sum(1 for v in xlarge.values() if v == False)
    
    return num_standard, num_large, num_xlarge

def set_true(nodelist, ptn_dict):
    '''
    This function sets values of partition node dictionary to be true if a job is running on it
    '''

    # if running on idle partition, job could span different partition nodes, e.g. r1n, r1t, r2t etc
    #n = re.findall(r'' + pid + '(?:\[.*?\]|\d+)', nodelist)[0]
    nodes = re.findall(r'(?:\d+-\d+|\d+)', nodelist[3:])
    for node in nodes:
        if '-' in node:
            n1 = int(node[0:2])
            n2 = int(node[3:5])
            for n in range(n1,n2+1):
                ptn_dict[str(n).zfill(2)] = True
        else:
            ptn_dict[node] = True

    return ptn_dict

if __name__ == "__main__":
    # Print time
    current_time = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    print('TIME: ' + current_time)
    
    # Get number of available nodes
    num_standard, num_large, num_xlarge = get_num_available_nodes()

    # Print number of available nodes for each partition
    print('AVAILABLE NODES:')
    print('STANDARD:', num_standard)
    print('LARGE:', num_large)
    print('XLARGE:', num_xlarge)

    # Clean up
    run('rm s.out', shell=True)
