import struct 
import time

filename = input("Enter the name of the file you would like to convert: ")
fnproc = open("nprocs.conf", "rb")

try:
    f = open(filename, "rb")
except OSError:
    print("Could not open/read file:" + filename)
    sys.exit()
except FileNotFoundError:
    print(f"File {filename} not found.  Aborting")
    sys.exit(1)

with f:
	num_cores = struct.unpack('i', f.read(4))[0]
	num_cores = struct.unpack('i', fnproc.read(4))[0]
	list_elements = []
	old_matrix = [[],[],[]]
	new_matrix = []
	print("Reading " + filename + "...")
	tic = time.perf_counter()
	for num in range(num_cores):
		num_elements = struct.unpack('i', f.read(4))[0]
		list_elements.append(num_elements)
	for num in range(num_cores):
		for i in range(list_elements[num]):
			n = struct.unpack('i', f.read(4))[0]
			old_matrix[0].append(n)
	for num in range(num_cores):
		for i in range(list_elements[num]):
			k = struct.unpack('i', f.read(4))[0]
			old_matrix[1].append(k)
	for num in range(num_cores):
		for i in range(list_elements[num]):
			t = struct.unpack('d', f.read(8))[0]
			old_matrix[2].append(t)
	toc = time.perf_counter()
	print("Reading " + filename + f" took {toc - tic:0.4f} seconds ")
	tic = time.perf_counter()
	numzero=0
	for num in range(len(old_matrix[0])):
		if old_matrix[2][num] != 0:
			new_matrix.append([old_matrix[1][num],old_matrix[0][num],old_matrix[2][num]])
		else:
			numzero=numzero+1
	toc = time.perf_counter()
	print("Removed " + str(numzero) + " zero-valued matrix elements")
	print("Re-packaging " + filename + f" took {toc - tic:0.4f} seconds ")
	tic = time.perf_counter()
	new_matrix = sorted(new_matrix, key=lambda x: x[1])
	toc = time.perf_counter()
	print("Sorting " + filename + f" took {toc - tic:0.4f} seconds ")

with open(filename[0:4] + filename[5:9], "wb") as f:
	f.write(struct.pack('i',len(new_matrix)))
	print("# of matrix elements=" + str(len(new_matrix)))
	for num in range(len(new_matrix)):
		f.write(struct.pack('iid', new_matrix[num][0], new_matrix[num][1], new_matrix[num][2]))

print(filename + " has been sorted")
