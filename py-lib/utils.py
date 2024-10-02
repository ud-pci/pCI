from subprocess import run, PIPE, CalledProcessError
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