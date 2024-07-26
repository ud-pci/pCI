dtm - density transition matrix 
-------------------------------

The ``dtm`` program calculates matrix elements of one-electron operators between many-electron states, under the density (or transition) matrix formalism. This formalism allows us to express the matrix elements between many-electron states via one-electron matrix elements. The ``dtm`` program forms these reduced density (or transition) matrices and calculates the reduced matrix elements. The following quantities can be calculated from this program:  

* electron g-factors  
* magnetic dipole and electronic quadrupole hyperfine structure constants :math:`A` and :math:`B`  
* electric :math:`Ek` and magnetic :math:`Mk` multipole transition amplitudes, where :math:`k = 1,2,3` corresponds to the dipole, quadrupole, and octupole transitions  
* nuclear spin independent parity nonconserving PNC amplitude  
* amplitude of the electron interaction with the P-odd nuclear anapole moment AM 
* P, T-odd interaction of the dipole electric dipole moment  
* nucleus magnetic quadrupole moment

This program begins by reading the file ``CONF.INP`` for system parameters and the list of configurations. Next, basis radial orbitals are read from the file ``CONF.DAT``, and radial integrals for all operators are calculated and written to the file ``DTM.INT``. If this file already exists, ``dtm`` uses it and does not recalculate the radial integrals. 

For the diagonal matrix elements, the list of determinants and eigenvectors corresponding to the state of interest are read from the files ``CONF.DET`` and ``CONF.XIJ``, respectively. For the non-diagonal matrix elements, the initial state is read from the file ``CONF.DET`` and ``CONF.XIJ``, and the final state is read from the files ``CONF1.DET`` and ``CONF1.XIJ``. The results of the diagonal and non-diagonal matrix elements are written to the files ``DM.RES`` and ``TM.RES``, respectively. 

``dtm`` takes in as input the input file ``dtm.in``:

.. code-block:: 

    2            # 1 for DM (as g-factor or hyperfine), 2 for transitions
    1  1 12      # from 1st even level to 1st-12th odd levels

.. note::

    Some :math:`3J`-coefficients might be zero in some cases, such as trying to compute :math:`E1` matrix element for  :math:`5s^2\, {}^1S_0 \rightarrow 5s5p\,{}^3P_1`. This will fail if the odd run had :math:`J=0`, :math:`M_J=0`. You would need to have an odd run with :math:`J=1`, :math:`M_J=1`.


Running dtm
~~~~~~~~~~~
Next, you must create a file named ``dtm.in`` with the following parameters:

.. code-block:: 

    Kl1                 ! (1 - density matrices, 2 - transition matrices, 3 - form DTM.INT)
    l01 l11 l10 l11     ! range of initial and final levels for desired matrix elements
                        ! l01 = initial level of first file, l11 - final level of first file
                        ! l10 - final level of second file, l11 - final level of second file
                        ! e.g. 1 3 1 4 corresponds to matrix elements between 1st to 3rd level of file 1 and 1st to 4th level of file 2 
    operator1           ! (optional) first matrix operator to summarize
    operator2           ! (optional) second matrix operator to summarize
    ...                 ! matrix operators: E1, E2, E3, M1, M2, M3, EDM, PNC, AM, MQM, A_hf, B_hf, GF

To run parallel ``dtm``, run the command:

.. code-block:: 

    mpirun -n <nprocs> ./dtm