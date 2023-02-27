from subprocess import run
import sys

# Read input to add from add.in if it exists, otherwise create it
try:
    open("add.in", 'rb').read()
except:
    f = open("add.in", "w")
    f.write("0")
    f.close()

# Run script to generate ADD.INP for even and odd configurations
run("python3 create_add_inp.py", shell=True)

# Check to see if ADDeven.INP and ADDodd.INP were created
try:
    open("ADDeven.INP", 'r')
except Exception as e:
    print(e)
    sys.exit()
try:
    open("ADDodd.INP", 'r')
except Exception as e:
    print(e)
    sys.exit()

# Generate CONFeven.INP from ADDeven.INP
run("cp ADDeven.INP ADD.INP", shell=True)
run("add < add.in", shell=True)
run("cp CONF.INP CONFeven.INP", shell=True)
print("CONFeven.INP created")

# Generate CONFodd.INP from ADDodd.INP
run("cp ADDodd.INP ADD.INP", shell=True)
run("add < add.in", shell=True)
run("cp CONF.INP CONFodd.INP", shell=True)
print("CONFodd.INP created")

# Cleanup - remove add.in
run("rm add.in", shell=True)