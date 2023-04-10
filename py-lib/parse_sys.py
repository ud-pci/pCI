import os
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT, run

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

if __name__ == "__main__":
    dir_path = os.getcwd()
    run('sinfo > s.out', shell=True)
    lines = read_file_lines('s.out')
    num_free_nodes = {'standard' : 0, 'large-mem' : 0, 'xlarge-mem' : 0}
    for line in lines:
        partition = ['standard', 'large-mem', 'xlarge-mem']
        
        if line.split()[0].replace('*','') in partition and line.split()[4] == 'idle':
            num_free_nodes[line.split()[0].replace('*','')] = line.split()[3]
    
    print('FREE NODES:')
    for item in num_free_nodes:
        print(item, ':', num_free_nodes[item])

    #num_req_nodes = input('How many nodes would you like to use? ')
    #req_mem = input('How much memory is required per core? ')