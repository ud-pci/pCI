from subprocess import run
import os
import sys
import re
from pathlib import Path
from gen_job_script import write_job_script
    
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
    run("add < add.in", shell=True)
    run("cp CONF.INP CONFeven.INP", shell=True)
    print("CONFeven.INP created")

if odd_exists:
    # Generate CONFodd.INP from ADDodd.INP
    run("cp ADDodd.INP ADD.INP", shell=True)
    run("add < add.in", shell=True)
    run("cp CONF.INP CONFodd.INP", shell=True)
    print("CONFodd.INP created")

# Cleanup - remove add.in
run("rm add.in", shell=True)

# Ask if user wants even and odd directories to be generated
gen_dir = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Generate directories? ")))))

if gen_dir:
    dir_path = os.getcwd()
    for dirs in ['EVEN','ODD']:
        Path(dir_path+'/'+dirs).mkdir(parents=True, exist_ok=True)
        os.chdir(dirs)
        if os.path.isfile('../HFD.DAT'):
            run("cp ../HFD.DAT .", shell=True)
        if os.path.isfile('../CONF' + dirs.lower() + '.INP'):
            run("cp ../CONF" + dirs.lower() + ".INP CONF.INP", shell=True)
        if os.path.isfile('../SGC.CON'):
            run("cp ../SGC.CON .", shell=True)
        if os.path.isfile('../SCRC.CON'):
            run("cp ../SCRC.CON .", shell=True)
        if os.path.isfile('../SGC.CON') and os.path.isfile('../SCRC.CON'):
            with open('c.in', 'w') as f:
                f.write('2, 2, 0, 0, 1')
                f.close()
        else:
            with open('c.in', 'w') as f:
                f.write('0, 0, 0, 0, 1')
                f.close()
        os.chdir('../')
    
    # Ask if user wants to submit ci jobs
    run_ci = eval(re.sub('(no|No|n|N|false)', 'False', re.sub('(yes|Yes|y|Y|true)', 'True', str(input("Submit CI jobs? ")))))
    if run_ci:
        for dirs in ['EVEN','ODD']:
            # Change to working directory
            os.chdir(dirs)
            
            # Check if files required for CI exist
            if not os.path.isfile('HFD.DAT'):
                print('HFD.DAT is missing from ' + dirs + ' directory')
                sys.exit()
            if not os.path.isfile('CONF.INP'):
                print('CONF.INP is missing from ' + dirs + ' directory')
                sys.exit()
            if not os.path.isfile('c.in'):
                print('c.in is missing from ' + dirs + ' directory')
                sys.exit()
                
            # Generate new job script if it doesn't exist yet
            if not os.path.isfile('ci.qs'):
                print('generating new ci.qs in ' + dirs + ' directory')
                write_job_script('ci', 5, 64, True, 0, 'safrono')
            
            # Submit job script
            run("sbatch ci.qs", shell=True)
                
            # Change to root directory
            os.chdir('../')


print('add script completed')

