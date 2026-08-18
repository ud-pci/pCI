from subprocess import run, PIPE, CalledProcessError
import re
import sys

def run_shell(command):
    """ This function runs an inputted shell command and throws an error if is not successful """
    try:
        run(command, shell=True, check=True, stdout=PIPE)
    except CalledProcessError as e:
        print(e)
        sys.exit()
        
def get_dict_value(param_dict, key):
    """ This function obtains returns the value from a dictionary given a key and returns None if the key doesn't exist """
    try:
        return param_dict[key]
    except KeyError:
        return None
    except TypeError:
        return None
    
def convert_params_to_list(param):
    """ Converts a string parameter into a list if it's not already a list"""
    
    if ',' in str(param): param = str(param).split(',')
    
    if not isinstance(param, list): 
        param = [str(param)]
    else:
        for istr in range(len(param)):
            param[istr] = param[istr].strip()

    return param

def get_basis_dir_name(basis):
    """Get directory name from basis string.
    Single-n basis (e.g. '17spdfg') -> '17g'. Mixed-n basis (e.g. '8spdf6g') -> full string."""
    n_array = re.findall('[0-9]+', basis)
    l_array = re.findall('[spdfghikl]+', basis)
    if not n_array or not l_array:
        return basis
    if len(n_array) == 1:
        return f'{n_array[0]}{l_array[0][-1]}'
    return basis

def check_slurm_installed():
    """ Check for the scontrol command """
    
    try:
        result = run(['scontrol', '--version'], stdout=PIPE)
        return result.returncode == 0
    except FileNotFoundError:
        return False
