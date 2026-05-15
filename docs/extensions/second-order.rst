MBPT package
------------

The MBPT part of the package calculates corrections to the bare Hamiltonian due to the core shells using second order MBPT for the CI code ``conf``, but for a much larger part of the Hamiltonian than the all-order code since high accuracy is not required for corrections associated with higher orbitals.  
The MBPT package consists of a single code ``second-cis``, which can be omitted, but then ``conf`` will not include electronic correlations associated with any of the core shells. If the all-order calculation was carried out, it will overwrite the second-order results with the all-order results where available. Such overlay of the MBPT and the all-order parts drastically improves the efficiency of the method.

The following code execution reads the files ``hfspl.1``, ``hfspl.2`` and ``HFD.DAT`` (to read the list of orbitals), and writes the files ``SGC.CON`` and ``SCRC.CON`` used by ``conf``. It also takes in the input file ``inf.vw`` and writes the results to the respective ``out.sdvw`` files.

.. code-block:: 

 	./second-cis <inf.vw >out.second.vw   # writes SGC.CON and SCRC.CON files

.. collapse:: Click here to see a description of inf.vw

	.. code-block:: 

		12             	 # ncore - number of core shells 
		1 -1           	 # n and kappa of the core shells
		2 -1
		2  1
		2 -2
		3 -1
		3  1
		3 -2
		3  2
		3 -3
		4 -1
		4  1
		4 -2
		35  5   		 # nmax and lmax in correlation diagram summations   
		195				 # Nmax  (max_orb) # of orbitals (from BASS.INPT) for how many Sigmas to calculate
		3				 # Lmax (lvmax) for how many Sigmas to calculate, supersedes Nmax. Maybe best to set to 4 when 4f is 	important
		 9				 # Kmax (kvmax) for higest partial wave for Sigma2 S_Kmax(ijkl). May increase for f shell cases
		 195 100 		 # nav1 nav - maximum average of orbitals to calculate Sigma for nav1=(i+j)/2 for Sigma(ij) and nav=(i+j	+k+l)/4 for Sigma_K(ijkl)
		0      			 # if just second order, 1 to read all-order input from PAIR.VW and sigma1
		1     			 # kval, set the same as in the all-order, see explanations of inf.aov

		==============================================================
		Additional explanations about Sigma_1 and Sigma_2 restrictions:
		(1) Put the first number (Sigma 1) to be the last orbital you plan to list in ADD.INP. 
		Look up the number in BASS.INP and add the number of core shells. 
		Here, the last orbital to be included in conf is 21d5/2, so the number is 195. 
		Example: Last orbital in BASS.INP:
		183  2.1201  3  2.1201      # Here, the last orbital to be included in conf is 21d5/2, so the number is 183+12 = 195
		Note: this is inputted twice, keep it the same. 
		(2) Set 100 as the second value. This was tested a while ago, can increase.

|