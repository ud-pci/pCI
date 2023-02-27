from subprocess import Popen, PIPE, STDOUT, run

add_inp = open("add.in", 'rb').read()

run(["python3", "create_add_inp.py"])

run(["cp","ADDeven.INP","ADD.INP"])
proc = Popen(['add'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
proc.communicate(input=add_inp)[0]
run(["cp","CONF.INP","CONFeven.INP"])
print("CONFeven.INP created")
proc.kill()

run(["cp","ADDodd.INP","ADD.INP"])
proc = Popen(['add'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
proc.communicate(input=add_inp)[0]
run(["cp","CONF.INP","CONFodd.INP"])
print("CONFodd.INP created")
proc.kill()
