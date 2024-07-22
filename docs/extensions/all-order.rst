All-order package
-----------------

The all-order part of the package calculates corrections to the bare Hamiltonian due to the core shells for the CI code ``conf``. The *all-order* label refers to inclusion of the large number of terms (second-, third-, fourth-order, etc.) in order-by-order many-body perturbation theory expansion using the iterative solutions until the sufficient numerical convergence is achieved. This code version implements a variant of the linearized single double coupled-cluster (CC) method. This CC version has been developed specifically for atoms fully utilizing atomic symmetries and capable of being efficiently ran with very large basis sets (over 1000 orbitals), reaching negligible numerical uncertainty associated with the choice of basis set.

The all-order package consists of three codes: ``allcore-ci``, ``valsd-ci``, and ``sdvw-ci``, which calculates core, core-valence and valence-valence excitations, respectively. These codes store resulting data in ``SGC.CON`` (small file) and ``SCRC.CON`` (up to a few GB file). The all-order package can be omitted if high precision is not required, leading to a method referred to as CI+MBPT.

The following code executions read the files ``hfspl.1`` and ``hfspl.2`` and generate the files ``SGC.CON`` and ``SCRC.CON``. It also takes in the input file ``inf.aov`` and writes the results to the respective ``out.`` files.

.. code-block:: 

	./allcore-rle-ci <inf.aov >out.core      # computes Sigma_ma, Sigma_mnab -> writes pair.3 file
	./valsd-rle-cis <inf.aov >out.val        # computes Sigma_mv (and thus Sigma_1, the one-body correcton to the Hamiltonian), Sigma_mnva -> writes val2 and sigma files
	./sdvw-rle-cis <inf.aov >out.sdvw        # computes Sigma_mnvw (Sigma_2, the two-body correction to the Hamiltonian) -> writes pair.vw and sigma1 files


.. raw:: html

	<details>
	<summary>Click here to see a description of inf.aov</summary>

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
	35  5    	     # nmax and lmax in correlation diagram summations   
	0   0    	     # internal parameters, do not change
	30       	     # max number of iterations for core
	2 7 1  	         # Stabilizer code parameters, see below 
	0.d0    	     # damping factor, not active not zero
	1      	         # kval (key_en) - key for energies,  see explanation below                          
	24    	         # nval - number of valence orbitals in list to follow   
	 5  -1 30        # n and kappa of the valence orbitals, max number of iteration 
	 6  -1 30
	 7  -1 30
	 8  -1 30
	 5   1 30
	 6   1 30
	 7   1 30
	 8   1 30
	 5  -2 30
	 6  -2 30
	 7  -2 30
	 8  -2 30
	 4   2 30
	 5   2 30
	 6   2 30
	 7   2 30
	 4  -3 30
	 5  -3 30
	 6  -3 30
	 7  -3 30
	 4   3 30
	 5   3 30
	 4  -4 30
	 5  -4 30 

	==============================================================
	Additional explanations:
	2 7 1  # Stabilizer code parameters 
	2 - DIIS method (change to 1 for RLE method)
	7 - number in iterations before stabilizer code runs
	1 - do not change (controls the type of linear algebra code)

.. raw:: html

	</details>

|

.. note::

	It is important to always check all output files. In ``out.core`` and ``out.val``, check that the core and valence orbitals have converged, respectively. If any cases completed 30 iterations (max) with values slowly drifting up, you can reduce the number of interactions to 3 or 5, and rerun the respective codes. Here, it's important to look at how much energies fluctuated during the first few iterations. In cases of severe divergence, check all inputs for errors. If nothing is found, and the atom/ion with closed d shell, but not p shell, set ``kval=2``.

.. note::
	
	``kval`` controls the values of :math:`\tilde{\epsilon}_v` in denominators (see `Phys. Rev. A 80, 012516 (2009) <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.80.012516>`_ for formulas).  
	``kval=1`` is the default choice, where :math:`\tilde{\epsilon}_v` for all :math:`nl_j` is set to the DHF energy of the lowest valence :math:`n` for the particular partial wave. For example for Sr, :math:`\tilde{\epsilon}_v` for :math:`v=ns` is set with the DHF energy of the :math:`5s` state, :math:`v=np_{1/2}` is set with the DHF energy of the :math:`5p_{1/2}` state, :math:`v=nd_{3/2}` is set with the DHF energy of the :math:`5d_{3/2}` state, and so on.  
	``kval=2`` is only used when the all-order valence energies are severely divergent. So far, this was observed with highly-charged ions with filling :math:`p`-shell (e.g. Sn-like). By "severe", we mean that the energies begin to diverge after the first or second iteration, immediately driving the correlation energy to be very large. Such a divergence cannot be fixed by a stabilizer. In this case, we have to identify which partial waves diverge and manually change the energies for these orbitals - we set them to the lowest DHF energy of the partial wave for which the all-order converged. For instance, if :math:`s` diverges, but :math:`p` does not, then set the :math:`ns` energies to the lowest :math:`np_{1/2}` DHF energies. The rest are left as DHF as in ``kval=1``.  
	The format would be as follows:  

	.. code-block:: 

		2                      		# kval                                         
		3                      		# lmax for the input to follow                
		0  -0.28000        	   	    # l=0 energies                                 
		1  -0.22000  -0.22000 	 	# l=1 energies p1/2 p3/2                       
		2  -0.31000  -0.31000     	# l=2 energies d3/2 d5/2		
		3  -0.13000  -0.13000      	# l=3 energies f5/2 f7/2