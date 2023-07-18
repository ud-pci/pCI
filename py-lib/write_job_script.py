
'''
This script writes job scripts for use on the University of Delaware DARWIN or Caviness clusters

'''

def write_job_script(file_name, job_name, mem, partition, time, env_setup, commands, parallel):
    '''
    This function writes the serial job script to submit via SLURM command sbatch
    
    Inputs:
    file_name : str
        the name of the job script file
    job_name : str
        the name of the job
    mem : int or str
        the amount of memory to allocate for the job
    partition : str
        the name of the partition to submit the job to
    time : int
        the number of days to allow the job to run for
    env_setup : str or list(str)
        commands to set up environment
    commands : str or list(str)
        commands to execute in batch job
    parallel : dict
        dictionary of parameters to include for a parallel job
    '''

    # Check if parallel job
    if parallel['include']:
        num_nodes = parallel['num_nodes']
        num_tasks_per_node = parallel['num_tasks_per_node']

    # Convert lists of commands to a string
    if isinstance(env_setup, list):
        env_setup = '\n'.join(line for line in env_setup)
        
    if isinstance(commands, list):
        commands = '\n'.join(line for line in commands)

    # Write job script file
    with open(file_name,'w') as f:
        f.write('#!/bin/bash -l \n\n')

        if parallel['include']:
            f.write('#SBATCH --nodes=' + str(num_nodes) + '\n')
            f.write('#SBATCH --tasks_per_node=' + str(num_tasks_per_node) + '\n')

        else:
            f.write('#SBATCH --ntasks=1 \n')

        f.write('#SBATCH --mem=' + str(mem) + ' \n')
        f.write('#SBATCH --job-name=' + job_name + ' \n')
        f.write('#SBATCH --partition=' + partition + ' \n')
        f.write('#SBATCH --time=' + str(time) + '-00:00:00 \n')
        f.write('#SBATCH --export=NONE \n\n')
        f.write('. /opt/shared/slurm/templates/libexec/common.sh \n')
        f.write(env_setup + '\n\n')
        f.write(commands + '\n')

    f.close()


if __name__ == "__main__":
    # General parameters for a serial job
    script_name = 'test.qs'
    mem = 0
    job_name = 'test'
    partition = 'standard'
    time = 5

    # Extra parameters for a parallel job
    pparams = {
        'include' : False,
        'num_nodes' : 1,
        'num_tasks_per_node' : 64
    }

    # Set up environment for job and run commands
    env_setup = 'vpkg_require pci'
    commands = ['time allcore-rle-ci <inf.aov >out.core',
                'time valsd-rle-cis <inf.aov >out.val',
                'time sdvw-rle-cis <inf.aov >out.vw',
                'time second-cis <inf.vw >out.second.vw']
    
    write_job_script(script_name, job_name, mem, partition, time, env_setup, commands, pparams)