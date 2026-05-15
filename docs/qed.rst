How to include QED corrections
==============================

*The following instructions assume familiarity with the main programs of the pCI package.*

Steps to run a QED calculation
------------------------------

In this section, we will introduce calculations including QED corrections. 

1. Construct basis set by running ``hfd`` and ``bass`` to obtain ``HFD.DAT``
2. Run ``sgc0`` to form an empty ``SGC.CON`` file
3. ``cp HFD.DAT HFD-noQED.DAT`` - save a copy of ``HFD.DAT`` without QED
4. Create a file ``qedpot.inp`` with number corresponding to the variant of the QED potential and the name of the ``HFD.DAT`` file.  

    .. collapse:: Click here to see a description of qedpot.inp.

        .. code-block:: 

            1            # kvar=1-5, Variant of QED potential
            HFD.DAT      # name of file holding basis sets
    
5. Create a file ``qed.in`` selecting options.
   
    .. collapse:: Click here to see a description of qed.in.

        .. code-block::

            1	# 1 for general diagonalization, 2 for first-order
            1	# 1 for no QED, 2 for QED
            2	# 0 for Coulomb, 1 for Gaunt, 2 for Full Breit

6. Loop until convergence:
   
	* Run ``qedpot_conf`` to construct selected QED potential
	* Run ``qed_rot`` to rotate orbitals to diagonalize Hamiltonian with QED corrections  

7. ``cp SGC.CON SGC-noQED.CON`` - save a copy of ``SGC.CON`` without QED
8. Run ``qedpot_conf``
9. Rename ``SGCqed.CON`` to ``SGC.CON``
10.  Run ``conf``

.. collapse:: You can also use the following batch.qed script for steps 3-9, making sure to change inputs relevant to your job.

    .. code-block:: 

        #! /bin/bash -fe
        kvar=1  # variant of QED potential
        bin='./'
        qedpot=$bin'qedpot_conf'
        qedrot=$bin'qed_rot'
        iter=25 # max number of iterations
        #####################################
        cat >qedpot.inp <<EndofFile
         $kvar
         HFD.DAT
        EndofFile
        #####################################
        cat >qed.in <<EndofFile
        1   # diagonalization
        1   # 1 for noQED, 2 for QED
        2   # 0 for Coulomb, 1 for Gaunt, 2 for Full Breit

        EndofFile
        #####################################
        n=1
        while [ $n -lt $iter ]; do
        echo 'Iteration '$n
        $qedpot >qp.res
        $qedrot <qed.in >qr.res
        grep 'changed' "QED_ROT.RES"
          if grep -q reached "QED_ROT.RES"; then
          echo 'Converged in '$n' iterations'
          break
          fi
          let n=n+1
        done
        cp SGC.CON SCG-noQED.CON
        ./qedpot_conf >qp.res 
        cp SGCqed.CON SGC.CON
        #####################################

|