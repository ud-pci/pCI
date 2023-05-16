from subprocess import run
import sys
import os
import re
from pathlib import Path

# Read input to add from add.in if it exists, otherwise create it
try:
    open("add.in", 'rb').read()
except:
    f = open("add.in", "w")
    f.write("0 \n 0")
    f.close()

# Run script to generate ADD.INP for even and odd configurations
run("python3 create_add_inp.py", shell=True)

# Check if ADDeven.INP and ADDodd.INP exist
even_exists = os.path.isfile('ADDeven.INP')
odd_exists = os.path.isfile('ADDodd.INP')

if even_exists:
    # Generate CONFeven.INP from ADDeven.INP
    run("cp ADDeven.INP ADD.INP", shell=True)
    run("add_nr < add.in", shell=True)
    run("cp CONF.INP CONFeven.INP", shell=True)
    print("CONFeven.INP created")

if odd_exists:
    # Generate CONFodd.INP from ADDodd.INP
    run("cp ADDodd.INP ADD.INP", shell=True)
    run("add_nr < add.in", shell=True)
    run("cp CONF.INP CONFodd.INP", shell=True)
    print("CONFodd.INP created")

# Cleanup - remove add.in
run("rm add.in", shell=True)

# Create even and odd directories
gen_dir = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Generate directories? ")))))

if gen_dir:
    dir_path = os.getcwd()
    for dirs in ['EVEN','ODD']:
        Path(dir_path+'/'+dirs).mkdir(parents=True, exist_ok=True)
        os.chdir(dirs)
        run("cp ../HFD.DAT .", shell=True)
        run("cp ../CONF" + dirs.lower() + ".INP CONF.INP", shell=True)
        if os.path.isfile('../SGC.CON'):
            run("cp ../SGC.CON .", shell=True)
        if os.path.isfile('../SCRC.CON'):
            run("cp ../SCRC.CON .", shell=True)
        os.chdir('../')

print('add script completed')

