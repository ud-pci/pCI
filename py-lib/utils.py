from subprocess import run, PIPE, CalledProcessError
import sys

def run_shell(command):
    """ This function runs an inputted shell command and throws an error if is not successful """
    try:
        run(command, shell=True, check=True, stdout=PIPE)
    except CalledProcessError as e:
        print(e)
        sys.exit()