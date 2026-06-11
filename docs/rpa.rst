RPA corrections
===============

*The following instructions assume familiarity with the main programs of the pCI package.*

The Random Phase Approximation (RPA) accounts for core-polarization corrections to radial matrix elements of external-field operators. When an external perturbation (e.g., an electromagnetic field) acts on the atom, it distorts the core electron cloud in addition to acting directly on the valence electrons. RPA captures this effect self-consistently by solving a set of linear equations for the induced change in each core orbital.

The corrections are computed for radial matrix elements of the following operators:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Label
     - Operator
   * - ``A_hf``
     - Hyperfine constant A
   * - ``B_hf``
     - Hyperfine constant B
   * - ``E1_L``
     - Electric dipole (length form)
   * - ``E1_V``
     - Electric dipole (velocity form)
   * - ``E2``
     - Electric quadrupole
   * - ``E3``
     - Electric octupole
   * - ``M1``
     - Magnetic dipole
   * - ``M2``
     - Magnetic quadrupole
   * - ``M3``
     - Magnetic octupole
   * - ``EDM``
     - Electric dipole moment
   * - ``PNC``
     - Parity non-conserving amplitude
   * - ``ANM``
     - Anapole moment
   * - ``MQM``
     - Magnetic quadrupole moment

Programs
--------

rpa
~~~

Solves the RPA equations and writes corrected radial integrals for each operator. It runs in two modes selected interactively:

* **Core mode** (``1``) — solves the RPA equations self-consistently for the core shells and writes core-core matrix elements.
* **Valence mode** (``2``) — computes RPA corrections for valence shells. This is the mode used in the typical workflow.

Each run produces one ``RPA_N.INT`` file per active operator, along with a summary in ``RPA.RES``.

Input files:

* ``MBPT.INP`` — defines the shell structure, active operators, and RPA options
* ``HFD.DAT`` — orbitals from ``hfd``
* ``CONF.DAT`` — CI wave functions (valence mode only)

Output files:

* ``RPA.RES`` — printed summary of RPA corrections
* ``RPA_N.INT`` — corrected radial integrals for each operator (one file per operator)

rpa_dtm
~~~~~~~

Substitutes the zeroth-order radial integrals in ``DTM.INT`` with the RPA-corrected integrals from the ``RPA_N.INT`` files produced by ``rpa``. Run this after ``pdtm`` has generated ``DTM.INT``.

Input files:

* ``DTM.INT`` — transition matrix integrals from ``pdtm``
* ``RPA_N.INT`` — corrected integrals from ``rpa``

Output files:

* ``RPA_DTM.RES`` — summary of substituted integrals
* ``DTM.INT`` — updated in place with RPA-corrected integrals

Workflow
--------

1. Run ``pdtm`` to compute non-RPA radial integrals and produce ``DTM.INT``.
2. Run ``rpa`` in valence mode (enter ``2`` when prompted) to compute RPA-corrected radial integrals.
3. Run ``rpa_dtm`` to substitute the integrals in ``DTM.INT`` with the RPA-corrected ones.

``MBPT.INP`` configuration
--------------------------

The ``rpa`` program reads shell and operator settings from the ``MBPT`` and ``RPA`` blocks in ``MBPT.INP``:

.. code-block::

    MBPT
     Nso=  n     # number of core shells (CI core)
     Nsh=  n     # last shell included in the self-consistent RPA field
     Nss=  n     # last virtual shell for the RPA sum
     Nsv=  n     # first valence shell (usually Nso+1)
     Nmax= n     # last valence shell for effective radial integrals
     Lmax= n     # maximum angular momentum for valence integrals
           n     # max multipolarity (not used)
     Kt2=  n     # accuracy for Coulomb radial integrals
     Kbrt= n     # Breit corrections (0 = off, 1 = on)

    RPA
     n           # A_hf  (1 = compute, 0 = skip)
     n           # B_hf
     n           # E1_L
     n           # EDM
     n           # PNC
     n           # E1_V
     n           # ANM
     n           # MQM
     n           # M1
     n           # E2
     n           # E3
     n           # M2
     n           # M3
          n      # Nhf: number of core shells in the SCF procedure
          n      # Kmg: frequency type (0 = static)
      n.nnnnnn   # Omega: perturbation frequency (ignored if Kmg=0)
          n      # Kex: include exchange (0 = skip, 1 = include)
