conf - configuration interaction 
--------------------------------

In this section, we will introduce how to run the CI program ``conf``. 

Here is a summary of the input and output files used in ``conf``.

Input Files:

* ``HFD.DAT`` - basis set radial orbitals :math:`\phi_{nlj}` and radial derivatives of the orbitals :math:`\partial_r\phi_{nlj}`
* ``CONF.DAT`` - basis set radial orbitals :math:`\phi_{nlj}` and functions :math:`\chi_{nlj}=h_\text{HF}^r\phi_{nlj}`, where :math:`h_\text{HF}^r` is the radial part of the Dirac-Fock operator
* ``CONF.GNT`` - relativistic Gaunt coefficients produced by ``basc``
* ``CONF.INT`` - relativistic Coulomb coefficients produced by ``basc``
* ``SGC.CON`` (optional) - one-electron effective radial integrals of the MBPT/all-order corrections
* ``SCRC.CON`` (optional) - two-electron effective radial integrals of the MBPT/all-order corrections
* ``CONF.INP`` - list of relativistic configurations and user defined parameters
* ``c.in`` - input file with keys ``Kl``, ``Ksig``, and ``Kdsig``
  
	* ``Kl`` = (0 - start, 1 - continue calculation, 2 - include corrections, 3 - add configurations)
	* ``Ksig`` = (0 - pure CI, 1 - include one-electron corrections, 2 - include one- and two-electron corrections)
	* ``Kdsig`` = (0 - no energy dependence on Sigma, 1 - energy dependence on Sigma)
	* ``Kw`` = (0 - do not write CONF.HIJ, 1 - write CONF.HIJ)
	* ``kLSJ`` = (0 - do not calculate expectation values of :math:`S^2` and :math:`L^2` and form approximate terms for each energy level, 1 - calculate :math:`LSJ`)

Output Files:

* ``CONF.DET`` - basis set of determinants
* ``CONF.HIJ`` - indices and values of the Hamiltonian matrix elements
* ``CONF.JJJ`` - indices and values of the matrix elements of the operator :math:`J^2`
* ``CONF.XIJ`` - quantum numbers, eigenvalues and eigenvectors of the Hamiltonian
* ``CONF.RES`` - final table of energy eigenvalues and the weights of all configurations contributing to each term

The following is a sample of the head of a ``CONF.INP`` for calculating the even-parity states of :math:`\text{Ir}^{17+}`.

.. code-block:: 

      Ir17+_even            # ion_parity                                               
      Z = 77.0              # atomic number
     Am =193.0              # atomic weight
      J =  2.0              # total angular momentum
     Jm =  2.0              # angular momentum projection
     Nso=  14               # number of closed core shells
     Nc = 481               # number of relativistic configurations
     Kv =   4               # Kv = (3 - use projections, 4 - no projections)           
     Nlv= 5                 # number of energy levels
     Ne =  14               # number of valence electrons                          
     Kl4=   1               # Kl4 = (1 - initial approx. from energy matrix, 2 - initial approx. from CONF.XIJ file)
     Nc4=  28               # number of relativistic configurations in initial approximation
     Gj =  0.0000           # 
    Crt4=  0.0001           # cutoff criteria for davidson convergence
    kout= 0                 # key for level of output (0 - low detail output, 1 - detailed output)
    Ncpt=  1132             # number of relativistic configurations in PT block (used in conf_pt)
    Cut0= 0.0001            # cutoff criteria for davidson convergence
    N_it=  100              # number of davidson iterations
    Kbrt= 1                 # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)
    Gnuc= 1.07              # gyromagnetic ratio
         0.1002     0.2002    -0.2102     0.2104     0.3002    -0.3102
         0.3104    -0.3204     0.3206     0.4002    -0.4102     0.4104
        -0.4204     0.4206
       1
       1-0.4306     0.4308
       2
       2-0.4306     0.4306     0.5002
       3-0.4305     0.4307     0.5002
       4-0.4304     0.4308     0.5002
       3
       5-0.4305     0.4308    -0.5101
       6-0.4305     0.4308     0.5101
       7-0.4306     0.4307    -0.5101
       8-0.4306     0.4307     0.5101
       :
    

.. note::

	The first 5 columns up to the '=' sign are fixed, and the program will give an error if there are any discrepancies here. 

.. note::

	The list of core shells are fixed to have a maximum of 6 shells per row.

In the ``CONF.INP`` file shown above, we include 481 relativistic configurations in the CI space, and 1132 relativistic configurations in the PT space (if ``conf_pt`` is to be used after ``conf``).


Running conf
~~~~~~~~~~~~
Next, you must create a file named ``c.in`` with the following parameters:

.. code-block:: 

    Kl         ! (0 - start, 1 - continue calculation, 2 - include corrections, 3 - add configurations)
    Ksig       ! (0 - pure CI, 1 - include one-electron corrections, 2 - include one- and two-electron corrections)
    Kdsig      ! (0 - no energy dependence on Sigma, 1 - energy dependence on Sigma)
    Kw         ! (0 - do not write CONF.HIJ, 1 - write CONF.HIJ)
    kLSJ       ! (0 - do not calculate expectation values of S^2 and L^2 and form approximate terms for each energy level, 1 - calculate LSJ)


To run parallel ``conf``, run the command:

.. code-block:: 

    mpirun -n <nprocs> ./conf