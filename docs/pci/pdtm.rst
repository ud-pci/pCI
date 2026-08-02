pdtm - density transition matrix
--------------------------------

The ``pdtm`` program calculates matrix elements of one-electron operators between many-electron states, under the density (or transition) matrix formalism. This formalism allows us to express the matrix elements between many-electron states via one-electron matrix elements. The ``pdtm`` program forms these reduced density (or transition) matrices and calculates the reduced matrix elements. The following quantities can be calculated from this program:

* electron g-factors
* magnetic dipole and electronic quadrupole hyperfine structure constants :math:`A` and :math:`B`
* electric :math:`Ek` and magnetic :math:`Mk` multipole transition amplitudes, where :math:`k = 1,2,3` corresponds to the dipole, quadrupole, and octupole transitions
* nuclear spin independent parity nonconserving PNC amplitude
* amplitude of the electron interaction with the P-odd nuclear anapole moment AM
* P, T-odd interaction of the dipole electric dipole moment
* nucleus magnetic quadrupole moment

This program begins by reading the file ``CONF.INP`` for system parameters and the list of configurations. Next, basis radial orbitals are read from the file ``CONF.DAT``, and radial integrals for all operators are calculated and written to the file ``DTM.INT``. If this file already exists, ``pdtm`` uses it and does not recalculate the radial integrals.

For the diagonal matrix elements, the list of determinants and eigenvectors corresponding to the state of interest are read from the files ``CONF.DET`` and ``CONF.XIJ``, respectively. For the non-diagonal matrix elements, the initial state is read from the files ``CONF.DET`` and ``CONF.XIJ``, and the final state is read from the files ``CONF1.DET`` and ``CONF1.XIJ``. The results of the diagonal and non-diagonal matrix elements are written to the files ``DM.RES`` and ``TM.RES``, respectively.

Here is a summary of the input and output files used in ``pdtm``.

Input Files:

* ``CONF.INP`` - list of relativistic configurations and system parameters
* ``CONF.DAT`` - basis set radial orbitals
* ``CONF.DET`` - list of determinants for the initial (or only) state
* ``CONF.XIJ`` - eigenvectors and eigenvalues for the initial (or only) state
* ``CONF1.DET`` (TM mode only) - list of determinants for the final state
* ``CONF1.XIJ`` (TM mode only) - eigenvectors and eigenvalues for the final state
* ``CONFSTR.RES`` (optional) - term symbols for the initial state, produced by ``pconf``; used to label output
* ``CONFSTR1.RES`` (optional, TM mode only) - term symbols for the final state
* ``DTM.INT`` (optional) - pre-computed radial integrals; if present, skips recalculation

Output Files:

* ``DTM.INT`` - radial integrals for all operators (written on first run, reused on subsequent runs)
* ``DM.RES`` - full density matrix results (DM mode)
* ``TM.RES`` - full transition matrix results (TM mode)
* ``DTM.RES`` - program log recording parameters and confirmation that ``DTM.INT`` was formed (Init mode)
* ``DM0.RES`` - intermediate density matrix elements (written in both DM and TM modes)

Each operator listed in ``Operators`` produces an additional output file:

*DM mode operators:*

* ``GF.RES`` - electron g-factors
* ``gQED.RES`` - QED correction to the g-factor
* ``A_hf.RES`` - magnetic dipole hyperfine constant :math:`A`
* ``B_hf.RES`` - electronic quadrupole hyperfine constant :math:`B`

*TM mode operators:*

* ``E1.RES`` - electric dipole transition amplitudes (both length and velocity gauges)
* ``E1_L.RES`` - electric dipole amplitude, length gauge only
* ``E1_V.RES`` - electric dipole amplitude, velocity gauge only
* ``E2.RES`` - electric quadrupole transition amplitudes
* ``E3.RES`` - electric octupole transition amplitudes
* ``M1.RES`` - magnetic dipole transition amplitudes
* ``M2.RES`` - magnetic quadrupole transition amplitudes
* ``M3.RES`` - magnetic octupole transition amplitudes
* ``PNC.RES`` - nuclear spin independent parity nonconserving amplitude
* ``AM.RES`` - electron interaction with the P-odd nuclear anapole moment
* ``EDM.RES`` - P, T-odd interaction of the electric dipole moment
* ``MQM.RES`` - nuclear magnetic quadrupole moment


``pdtm`` takes in the input file ``dtm.in``:

.. code-block::

    Mode = (DM, TM, Init)                           # DM - density matrix, TM - transition matrix, Init - form DTM.INT
    Levels = level_range_from (level_range_to)      # from level_range_from to level_range_to
    Operators = E1                                  # form final tables of listed operators

The value of ``Mode`` can take the following values:

* ``DM`` - form density matrix and calculate diagonal matrix elements
* ``TM`` - form transition matrix and calculate non-diagonal matrix elements
* ``Init`` - create the file ``DTM.INT``, which is used for ``rpa_dtm`` and ``pol`` programs

The value of ``Levels`` can take in a range of levels depending on whether ``Mode`` was specified as ``DM`` or ``TM``. The following are some examples.

To calculate DM for energy levels 1 to 5, with a table of g-factors:

.. code-block::

    Mode = DM
    Levels = 1 5
    Operators = GF

To calculate TM for energy levels 1 to 4 from ``CONF.XIJ`` to energy levels 4 to 9 from ``CONF1.XIJ``, with tables of :math:`E1` and :math:`M2` reduced matrix elements:

.. code-block::

    Mode = TM
    Levels = 1 5, 4 9
    Operators = E1, M2

The value of ``Operators`` takes in a comma-separated list of operators. The available operators depend on the ``Mode``:

*DM mode:* ``GF``, ``gQED``, ``A_hf``, ``B_hf``

*TM mode:* ``E1``, ``E1_L``, ``E1_V``, ``E2``, ``E3``, ``M1``, ``M2``, ``M3``, ``PNC``, ``AM``, ``EDM``, ``MQM``

.. note::

    ``Levels`` and ``Operators`` are not read if ``Mode == Init``.

.. note::

    ``E1`` computes both the length and velocity gauge amplitudes and writes them to a single ``E1.RES`` file. Use ``E1_L`` or ``E1_V`` to compute only one gauge.

.. note::

    Some :math:`3J`-coefficients might be zero in some cases, such as trying to compute :math:`E1` matrix element for :math:`5s^2\, {}^1S_0 \rightarrow 5s5p\,{}^3P_1`. This will fail if the odd run had :math:`J=0`, :math:`M_J=0`. One needs to re-run CI with :math:`J=1`, :math:`M_J=1`.


Running pdtm
~~~~~~~~~~~~

To run parallel ``pdtm``, run the command:

.. code-block::

    mpirun -n <nprocs> pdtm
