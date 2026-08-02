pol - polarizabilities
----------------------

The ``pol`` program calculates the dc and ac polarizabilities of the specified atomic states. The expression for electric-dipole ac polarizability at the frequency :math:`\omega` of the state :math:`|JM\rangle` can be written (in a.u.) as a sum over unperturbed intermediate states :math:`n`,

.. math::

    \alpha(\omega) = 2 \sum_n \frac{(E_n-E) |\langle JM|D_z|n\rangle|^2}{(E_n-E)^2 - \omega^2} ,

where :math:`D` is an electric dipole moment operator and :math:`E` and :math:`E_n` are the energies of the initial and intermediate states, respectively.

To find :math:`\alpha`, we can rewrite this expression as

.. math::

    \alpha(\omega) = \sum_n \langle JM|D_z|n\rangle \,\langle n|D_z| JM\rangle \left[ \frac{1}{E_n-E+\omega} + \frac{1}{E_n-E-\omega} \right] .  

Then we use the Sternheimer\ :footcite:p:`Ste50` or Dalgarno-Lewis\ :footcite:p:`DalLew55` method and solve the inhomogeneous equations

.. math::

    (H - E \pm \omega)\, |\delta \phi_{\pm} \rangle = D_z\, |JM \rangle,

to find :math:`|\delta \phi_{\pm} \rangle`,

.. math::

    |\delta \phi_{\pm} \rangle &= \frac{1}{H - E \pm \omega} D_z |JM \rangle \\
    &= \sum_n \frac{1}{H - E \pm \omega}|n\rangle \langle n |D_z |JM \rangle,

where :math:`H` is the Hamiltonian, and we use the closure relation :math:`\sum_n | n \rangle \langle n | = 1`. After that, the polarizability can be found as the sum of two matrix elements

.. math::

    \alpha(\omega) =  \langle JM| D_z |\delta \phi_{+} \rangle + \langle JM| D_z |\delta \phi_{-} \rangle .

If electrons in an atomic system are divided into valence and core electrons, the polarizability can be divided accordingly as

.. math::

    \alpha \equiv \alpha_v + \alpha_c,

where :math:`\alpha_v` and :math:`\alpha_c` are the valence and core contributions.

The program ``pol`` calculates only the valence polarizability :math:`\alpha_v`. The core polarizability needs to be computed separately with a different program. 

Disregarding the vector polarizability, we can present the expression for :math:`\alpha(\omega)` as the sum of the scalar and tensor parts,

.. math::

    \alpha(\omega) = \alpha_0 + \alpha_2 \, \frac{3M^2-J(J+1)}{J(2J-1)}.

The ``pol`` program gives both scalar and tensor polarizabilities if the latter is not zero. 

This program requires several input files from previously run ``pconf`` and ``pdtm`` programs, including ``CONF.DET`` and ``CONF.XIJ`` of the parity of the level of interest (renamed to ``CONF0.DET`` and ``CONF0.XIJ``), ``CONF.INP``, ``CONF.XIJ``, ``CONF.HIJ``, and ``CONF.JJJ`` of the opposite parity, and the file ``DTM.INT`` from ``pdtm``.

Here is a summary of the input and output files used in ``pol``.

Input Files:

* ``CONF.INP`` - list of relativistic configurations and system parameters
* ``CONF.DAT`` - basis set radial orbitals
* ``CONF.HIJ`` - sorted Hamiltonian matrix elements (from ``sort`` applied to ``CONFp.HIJ``)
* ``CONF.JJJ`` - sorted :math:`J^2` matrix elements (from ``sort`` applied to ``CONFp.JJJ``)
* ``CONF.XIJ`` - eigenvectors and eigenvalues of the opposite-parity state
* ``CONF0.DET`` - list of determinants for the state whose polarizability is computed
* ``CONF0.JJJ`` - :math:`J^2` matrix elements for the state whose polarizability is computed
* ``CONF0.XIJ`` - eigenvectors and eigenvalues for the state whose polarizability is computed
* ``DTM.INT`` - radial integrals from ``pdtm``

Output Files:

* ``POL.RES`` - scalar and tensor polarizabilities for the requested level and frequency ranges
* ``POL_E1.RES`` - intermediate E1 matrix elements used in the polarizability sum
* ``POL_p.XIJ`` - solution vector :math:`|\delta\phi_+\rangle`; written during calculation and re-read when ``Mode = 1``
* ``POL_m.XIJ`` - solution vector :math:`|\delta\phi_-\rangle`; written during calculation and re-read when ``Mode = 1``
* ``POL_J_p.XIJ`` - :math:`J`-decomposition of :math:`|\delta\phi_+\rangle`
* ``POL_J_m.XIJ`` - :math:`J`-decomposition of :math:`|\delta\phi_-\rangle`

For example, if we want to calculate polarizabilities for an even state:

.. code-block:: 

    cp CONFeven.DET CONF0.DET
    cp CONFeven.XIJ CONF0.XIJ
    
    cp CONFodd.INP CONF.INP
    cp CONFodd.XIJ CONF.XIJ
    cp CONFodd.HIJ CONF.HIJ
    cp CONFodd.JJJ CONF.JJJ

Note that ``pconf`` writes ``CONFp.HIJ`` and ``CONFp.JJJ`` in unsorted parallel format. Before running ``pol``, these must be sorted using the ``sort`` program, which produces the sorted serial files ``CONF.HIJ`` and ``CONF.JJJ`` that ``pol`` expects.

The ``pol`` program requires a list of key-value parameters in a file ``pol.in``:

.. code-block::

    Mode   = (0, 1)
    Method = (0, 1, 2)
    Level  = N
    Ranges = initial_wavelength final_wavelength step_size [, ...]
    IP1    = N

The value of ``Mode`` can take the following values:

* ``Mode = 0`` - start a new calculation
* ``Mode = 1`` - continue a previous calculation, reading ``POL_p.XIJ`` and ``POL_m.XIJ`` from disk

The value of ``Method`` can take the following values:

* ``Method = 0`` - invert the matrix and iterate until convergence; restart with iteration if inversion diverges
* ``Method = 1`` - invert the matrix only, without iterative refinement
* ``Method = 2`` - modified iteration that restarts every 2 iterations while retaining solution vectors; use when ``Method = 0`` diverges

The value of ``Level`` is the index of the eigenvector in ``CONF0.XIJ`` for which the polarizability is calculated.

The ``Ranges`` field takes a comma-separated list of wavelength ranges, each specified as ``initial_wavelength final_wavelength step_size`` (in nm). The special range ``0 0 0`` computes the dc (static) polarizability.

The value of ``IP1`` sets the dimension of the initial approximation subspace (default: ``IP1 = 15000``).

For example, to calculate dc polarizabilities, as well as ac polarizabilities from :math:`\lambda=500` nm to :math:`\lambda=505` nm in steps of :math:`0.5` nm, for the first state in ``CONF0.XIJ``, we can use the following ``pol.in``:

.. code-block::

    Mode = 0
    Method = 0
    Level = 1
    Ranges = 0 0 0, 500 505 0.5

The range ``0 0 0`` corresponds to calculations of dc polarizabilities, while ``500 505 0.5`` corresponds to calculating ac polarizabilities from :math:`\lambda=500` nm to :math:`\lambda=505` nm in steps of :math:`0.5` nm.

Before parallel ``pol`` can be run, the parallel matrix files ``CONFp.HIJ`` and ``CONFp.JJJ`` have to be sorted. This can be done using the ``sort`` program, which sorts the matrix elements of the operators :math:`H` and :math:`J^2` in order of ascending index (as done in the serial version of ``pconf``). It takes in ``CONFp.JJJ`` or ``CONFp.HIJ`` and writes a sorted serial version of the file, with an additional integer at the start storing the total number of matrix elements.

Running pol
~~~~~~~~~~~

To run ``pol``, run the command:

    .. code-block::

        pol


**References**

.. footbibliography::