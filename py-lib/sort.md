# sort.py

This is a python script that takes in a CONF.JJJ file written in parallel by parallel conf and sorts it in order of serial CONF.JJJ file.

CHANGES: any files that read CONF.JJJ before has to be modified
	1. read NumJ first before reading matrix elements
		e.g. con_cut has READ(11) numjj 
	2. remove reading index from file
		e.g. con_cut has READ(11) idum, K, N, T
			- remove idum from the list