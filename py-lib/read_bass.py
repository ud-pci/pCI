import struct
import re
import orbitals as orb_lib

def read_bass(filename):
    """Read bass input file""" 
    orbitals = []
    with open(filename,'r') as f:
        for line in f:
            line_strip = re.sub(' +', ' ', line.strip()).split(" ")
            if orb_lib.is_integer(line_strip[0]): 
                nl, j = orb_lib.convert_digital_to_char(line_strip[1])
                orbitals.append([nl, j])
    return orbitals

if __name__ == "__main__":
    orbitals = read_bass('BASS.INP')
    print(orbitals)
