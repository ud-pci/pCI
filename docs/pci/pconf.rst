pconf - configuration interaction 
---------------------------------

The parallel program ``pconf`` performs the configuration interaction method in the CI space defined by the list of configurations contained in the ``CONF.INP`` input file created by the ``add`` program. It takes the ``CONF.DAT``, ``CONF.GNT``, ``CONF.INT``, and ``CONF.INP`` files as input. In addition, ``pconf`` can also take in the ``SGC.CON`` and ``SCRC.CON`` files, which contain one- and two-electron effective radial integrals, respectively. 

Here is a summary of the input and output files used in ``pconf``.

Input Files:

* ``HFD.DAT`` - basis set radial orbitals :math:`\phi_{nlj}` and radial derivatives of the orbitals :math:`\partial_r\phi_{nlj}`
* ``CONF.DAT`` - basis set radial orbitals :math:`\phi_{nlj}` and functions :math:`\chi_{nlj}=h_\text{HF}^r\phi_{nlj}`, where :math:`h_\text{HF}^r` is the radial part of the Dirac-Fock operator
* ``CONF.GNT`` - relativistic Gaunt coefficients produced by ``basc``
* ``CONF.INT`` - relativistic Coulomb coefficients produced by ``basc``
* ``SGC.CON`` (optional) - one-electron effective radial integrals of the MBPT/all-order corrections
* ``SCRC.CON`` (optional) - two-electron effective radial integrals of the MBPT/all-order corrections
* ``CONF.INP`` - list of relativistic configurations and user defined parameters

Output Files:

* ``CONF.DET`` - basis set of determinants
* ``CONFp.HIJ`` - list of matrix elements of the Hamiltonian (generated only if ``Kw = 1`` in ``ci.in``)
* ``CONFp.JJJ`` - list of matrix elements of the operator :math:`J^2` (generated only if ``Kw = 1`` in ``ci.in``)
* ``CONF.XIJ`` - quantum numbers, eigenvalues and eigenvectors of the Hamiltonian
* ``CONF.ENG`` - tables of calculated quantum numbers :math:`J_i` and energy eigenvalues :math:`E_i` calculated each time ``CONF.XIJ`` is constructed during the Davidson procedure
* ``CONF.LVL`` - tables of top contributing configurations for each energy level calculated each time ``CONF.XIJ`` is constructed during the Davidson procedure
* ``CONF.RES`` - results of ``pconf`` program
* ``FINAL.RES`` - final table of quantum numbers :math:`J_i` and energy eigenvalues :math:`E_i` 
* ``LEVELS.RES`` - final table of top contributing configurations for each energy level
* ``CONFSTR.RES`` - list of top contributing configurations along with their atomic term symbol for each energy level

The following is a sample of the head of a ``CONF.INP`` for calculating the even-parity states of :math:`\text{Ir}^{17+}`. Here, we include 481 relativistic configurations in the CI space, and 1132 relativistic configurations in the PT space (if ``conf_pt`` is to be used after ``pconf``).

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


In addition to the input files listed above, ``pconf`` also requires the file ``ci.in``, which contains a list of key-value pairs defining the CI computation:

.. code-block::
  
	Kl = (0,1,2,3) - (0 - start, 1 - continue calculation, 2 - include corrections, 3 - add configurations)
	Ksig = (0,1,2) - (0 - pure CI, 1 - include one-electron corrections, 2 - include one- and two-electron corrections)
	Kdsig = (0,1) - (0 - no energy dependence on Sigma, 1 - energy dependence on Sigma)
	Kw = (0,1) - (0 - do not write matrices, 1 - write matrices)
	KLSJ = (0,1) - (0 - do not calculate LSJ, 1 - calculate LSJ)

The value of ``Kl`` can take the following values:

* ``Kl = 0`` - start a new CI calculation
* ``Kl = 1`` - continue a previous CI calculation
* ``Kl = 2`` - start a new CI calculation, including corrections from ``SGC.CON``
* ``Kl = 3`` - continue a previous CI calculation, adding new configurations

.. note::

    Note that for ``Kl=1,3`` to work, the file ``CONFp.HIJ`` must have been successfully written with ``Kw=1`` in the original run. ``CONFp.HIJ`` stores the previous Hamiltonian matrix, which must be read to continue a CI calculation. However, this file can be as big as 1 TB or more for large systems, so it should not be used unless the user has exceptional available computational resources. 

The value of ``Ksig`` can take the following values:

* ``Ksig = 0`` - pure CI
* ``Ksig = 1`` - include one-electron corrections from ``SCRC.CON``
* ``Ksig = 2`` - include one- and two-electron corrections from ``SCRC.CON``

The value of ``Kdsig`` can take the following values:

* ``Kdsig = 0`` - automatic approximation of the energy dependence of the operator :math:`\Sigma(E)`
* ``Kdsig = 1`` - manually specify the energy :math:`E_\mathrm{val}` to treat the energy dependence of :math:`\Sigma(E)`

The value of ``Kw`` can take the following values:

* ``Kw = 0`` - do not write the files ``CONFp.HIJ`` and ``CONFp.JJJ``
* ``Kw = 1`` - write the files ``CONFp.HIJ`` and ``CONFp.JJJ``

The value of ``KLSJ`` can take the following values:

* ``KLSJ = 0`` - do not calculate :math:`\langle S^2 \rangle`, :math:`\langle L^2 \rangle`, or form approximate terms for each energy level
* ``KLSJ = 1`` - calculate :math:`\langle S^2 \rangle`, :math:`\langle L^2 \rangle`, and form approximate term symbols for each energy level

After reading the keys from ``ci.in``, the ``pconf`` program reads the general parameters and the list of configurations from the ``CONF.INP`` file. Next, the information about the basis set is read from ``CONF.DAT``, relativistic Gaunt coefficients are read from ``CONF.GNT``, radial integrals are read from ``CONF.INT``, and optionally, effective radial integrals are read from ``SGC.CON`` and ``SCRC.CON``. 

Having read all required input files, the ``pconf`` program forms a list of determinants from the list of relativistic configurations and writes them to the file ``CONF.DET``. With the list of determinants, it forms the Hamiltonian matrix and then the matrix of the operator :math:`J^2`. The construction of the Hamiltonian matrix is the most time-consuming part of the ``pconf`` program. 

After the matrices are constructed, the ``pconf`` program enters the Davidson iterative procedure, where the Hamiltonian is diagonalized to obtain a specified number of low-lying energy eigenvalues and eigenvectors. The progress of the Davidson iterative procedure is written in the ``CONF.PRG`` file at each iteration. At selected intervals before the final iteration, the eigenvalues and eigenvectors are saved to ``CONF.XIJ``. Each time ``CONF.XIJ`` is written on disk, a table of the energy levels is appended to the file ``CONF.ENG``, and tables of the top contributing configurations for each level are appended to ``CONF.LVL``. 

Once the Davidson iterative procedure has converged, the final eigenvalues and eigenvectors are saved once again to ``CONF.XIJ``, the final energy table is saved to ``FINAL.RES``, the final list of the top contributing configurations for each level is saved to ``LEVELS.RES``, and the configurations along with their atomic term symbol is written to ``CONFSTR.RES``. 

If polarizability calculations are required, the Hamiltonian matrix elements will have to be written to the file ``CONFp.HIJ``by setting the key ``Kw=1`` in ``ci.in``. Note that depending on the size of the Hamiltonian, this file could take up hundreds of GB to over 1 TB. By default, ``Kw=0`` is set to not write the Hamiltonian to disk. As a reference, polarizability calculations have typically been done with Hamiltonian sizes of up to 4 million determinants. 


Running conf
~~~~~~~~~~~~

To run ``pconf``, run the command:

.. code-block:: 

    mpirun -n <nprocs> pconf