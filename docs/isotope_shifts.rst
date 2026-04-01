Isotope shift calculations
==========================

The transition frequency shift of an isotope :math:`A^\prime` compared to an isotope :math:`A` can be written as

.. math::

    \delta \nu_{AA^\prime} = (K_\mathrm{NMS} + K_\mathrm{SMS})(\frac{1}{A}-\frac{1}{A^\prime}) + K_\mathrm{FS} \delta \langle r^2 \rangle^{AA^\prime},

where :math:`A` and :math:`A^\prime` are mass numbers of the two isotopes and :math:`\langle r^2 \rangle` is the mean-square nuclear radius. The first term on the right represents the mass shift and the second term represents the field shift. The mass shift consists of two parts: a normal mass shift (NMS) and a specific mass shift (SMS). For more detail, see references below :footcite:p:`KorKoz07, SafPorKoz18`.

Normal mass shift
~~~~~~~~~~~~~~~~~
The normal mass shift coefficient is given by

.. math::

    K_\mathrm{NMS} = \frac{\nu_\mathrm{expt}}{1822.888},

where :math:`\nu_\mathrm{expt}` is the experimental frequency and the factor :math:`1/1822.888` is the electron mass expressed in atomic mass unit.

Field shift
~~~~~~~~~~~
We use the finite-field method to calculate field shift coefficients, in which the initial Hamiltonian is modified with an arbitrary perturbation coefficient :math:`\lambda`:

.. math:: 
    
    H\rightarrow H_\lambda = H + \lambda H_\textrm{FS},

where the field shift operator :math:`H_\textrm{FS}` modifies the Coulomb potential inside the nucleus. Here, the coefficient :math:`\lambda` has to be chosen to be sufficiently large enough to make the effect of the field shift significantly larger than the numerical uncertainty of the calculations, while also small enough to keep the change in the energy linear with :math:`\lambda`. Typically :math:`\lambda=0.01` meets both criteria. The energy eigenvalues :math:`E` are obtained by diagonalizing :math:`H_\lambda`, and then the field shift coefficient :math:`K_\mathrm{FS}` are found as

.. math::

    K_\mathrm{FS} = \frac{5}{6R^2} \frac{\partial E}{\partial \lambda}.

The conversion factor from atomic units to SI units used for the coefficient :math:`K_\mathrm{FS}` is :math:`1\,\mathrm{a.u.}=2.3497\times 10^{-3}\, \mathrm{GHz/fm}^2`.

Specific mass shift
~~~~~~~~~~~~~~~~~~~
The finite-field approach is also used to calculate the specific mass shift, in which the initial Hamiltonian is modified with the SMS operator:

.. math:: 
    
    H\rightarrow H_\lambda = H + \lambda H_\textrm{SMS}.

The SMS coefficient is given by the corresponding derivative 

.. math::
    
    K_\mathrm{SMS}=\partial E /\partial \lambda.

The conversion factor from atomic units to SI units for the coefficient :math:`K_\mathrm{SMS}` is :math:`1\,\mathrm{a.u.}=3609.46\, \mathrm{GHz/amu}`.

Running isotope shift calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To run isotope-shift (IS) calculations, a few additional keys must be added to the ``HFD.INP``, ``BASS.INP`` and ``CONF.INP`` files before running the CI procedure as usual. The relevant isotope-shift keys are:  

* ``K_is`` - selects the type of isotope shift to calculate  
  
    * ``K_is = 0`` - no isotope shift
    * ``K_is = 1`` - field shift (FS)
    * ``K_is = 2`` - specific mass shift (SMS)
    * ``K_is = 3`` - normal mass shift (NMS)
    * ``K_is = 4`` - total mass shift (SMS + NMS)
  
* ``C_is`` - scaling parameter :math:`\lambda` used to evaluate IS coefficients in the finite-field method (C_is = dR_nucl/R_nucl)
* ``Klow`` - lower component key 

.. note::

	For isotope shifts, it is important to recompile ``pbasc`` and ``pconf`` with arrays for 2-electron integrals set to ``double precision`` type. This can be done by setting the flag ``-DUSE_DP_INTEGRALS`` during compilation, or in the ``params.f90`` file by changing the value of ``Integer, Parameter :: type2_real = sp`` from ``sp``  to ``dp``. You can confirm this change in the title of the ``BASC.RES`` and ``CONF.RES`` output files, which should mention ``double precision for 2e integrals``.

via pCI-py scripts
------------------

Isotope shift calculations can also be automated using the pCI-py scripts. They generate basis sets and performs the corresponding CI calculations at different values of :math:`\lambda` to obtain :math:`E(\lambda=0)`, :math:`E(\lambda=\pm \delta \lambda)` and :math:`E(\lambda=\pm 2\delta\lambda)`. The field shift coefficient is then calculated using a five-point stencil finite difference formula using these energies.

To do this, set the keys in the ``config.yml`` file. For example, we can set our field shift calculation with ``C_is=0.01``:

.. code-block:: 

    optional:
        isotope_shifts: 
            include: True
            K_is: 1
            C_is: 0.01

References
----------

.. footbibliography::