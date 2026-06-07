Isotope shifts
==============

The transition frequency shift of an isotope :math:`A^\prime` compared to an isotope :math:`A` can be written as

.. math::

    \delta \nu_{AA^\prime} = (K_\mathrm{NMS} + K_\mathrm{SMS})\left(\frac{1}{A}-\frac{1}{A^\prime}\right) + K_\mathrm{FS} \delta \langle r^2 \rangle^{AA^\prime},

where :math:`A` and :math:`A^\prime` are mass numbers of the two isotopes and :math:`\langle r^2 \rangle` is the mean-square nuclear radius. The first term represents the mass shift and the second term represents the field shift. The mass shift consists of two parts: a normal mass shift (NMS) and a specific mass shift (SMS). For more detail, see references below :footcite:p:`KorKoz07, SafPorKoz18`.

Normal mass shift
~~~~~~~~~~~~~~~~~
The normal mass shift coefficient is given by

.. math::

    K_\mathrm{NMS} = \frac{\nu_\mathrm{expt}}{1822.888},

where :math:`\nu_\mathrm{expt}` is the experimental frequency and the factor :math:`1/1822.888` is the electron mass in atomic mass units.

Field shift
~~~~~~~~~~~
We use the finite-field method to calculate field shift coefficients, in which the initial Hamiltonian is modified with an arbitrary perturbation coefficient :math:`\lambda`:

.. math::

    H\rightarrow H_\lambda = H + \lambda H_\textrm{FS},

where the field shift operator :math:`H_\textrm{FS}` modifies the Coulomb potential inside the nucleus. The coefficient :math:`\lambda` must be large enough to make the field shift effect significantly larger than the numerical uncertainty, while small enough to keep the energy change linear in :math:`\lambda`. Typically :math:`\lambda=0.01` meets both criteria. The field shift coefficient :math:`K_\mathrm{FS}` is then found as

.. math::

    K_\mathrm{FS} = \frac{5}{6R^2} \frac{\partial E}{\partial \lambda}.

The conversion factor from atomic units to SI units for :math:`K_\mathrm{FS}` is :math:`1\,\mathrm{a.u.}=2.3497\times 10^{-3}\, \mathrm{GHz/fm}^2.`

Specific mass shift
~~~~~~~~~~~~~~~~~~~
The finite-field approach is also used to calculate the specific mass shift, in which the initial Hamiltonian is modified with the SMS operator:

.. math::

    H\rightarrow H_\lambda = H + \lambda H_\textrm{SMS}.

The SMS coefficient is given by the corresponding derivative

.. math::

    K_\mathrm{SMS}=\frac{\partial E}{\partial \lambda}.

The conversion factor from atomic units to SI units for :math:`K_\mathrm{SMS}` is :math:`1\,\mathrm{a.u.}=3609.46\, \mathrm{GHz/amu}`.

Running isotope shift calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To run isotope-shift calculations, the following keys must be added to ``HFD.INP``, ``BASS.INP``, and ``CONF.INP`` before running the CI procedure as usual:

* ``k_is`` - selects the type of isotope shift:

    * ``k_is = 0`` - no isotope shift (default)
    * ``k_is = 1`` - volume/field shift (VS)
    * ``k_is = 2`` - specific mass shift (SMS)
    * ``k_is = 3`` - normal mass shift (NMS)
    * ``k_is = 4`` - total mass shift (SMS + NMS)

* ``c_is`` - the finite-field scaling parameter :math:`\lambda`. For field shift, this corresponds to :math:`\delta R_\mathrm{nucl}/R_\mathrm{nucl}`.

* ``klow`` - lower component key (controls inclusion of small-component contributions)

.. note::

    Isotope shift calculations require ``pbasc`` and ``pconf`` to be compiled with double-precision two-electron integrals. Build with:

    .. code-block::

        cmake -DUSE_DP_INTEGRALS=ON ..

    You can confirm this in the header of ``BASC.RES`` and ``CONF.RES``, which will read ``double precision for 2e integrals``.

via pCI-py scripts
------------------

Isotope shift calculations can also be automated using the ``isotope_shifts.py`` script. It generates basis sets and performs CI calculations at different values of :math:`\lambda` to obtain :math:`E(\lambda=0)`, :math:`E(\lambda=\pm\delta\lambda)`, and :math:`E(\lambda=\pm 2\delta\lambda)`. The field shift coefficient is then calculated using a five-point stencil finite difference formula.

The script requires a completed non-IS calculation in a directory called ``0``, which must contain the relevant ``HFD.INP``, ``BASS.INP``, ``CONF.INP``, ``ci.in``, and job scripts. It reads ``K_is`` and ``C_is`` from the ``optional.isotope_shifts`` block of ``config.yml``:

.. code-block::

    optional:
        isotope_shifts:
            include: True
            K_is: 1
            C_is: 0.01

The script prompts the user for which mode to run:

* ``basis`` — copies the basis from the ``0`` directory into four IS subdirectories (``minus{c}``, ``plus{c}`` for :math:`c \in \{C\_is/2,\, C\_is\}`) and updates ``K_is`` and ``C_is`` in the input files for each.
* ``ci`` — submits the CI job script in each IS directory. Requires ``basis`` mode to have been run first.
* ``analysis`` — reads ``FINAL.RES`` from all five IS directories (including ``0`` for :math:`\lambda = 0`) and computes the IS coefficients using a five-point stencil finite difference formula.

References
----------

.. footbibliography::
