Overview of the All-order code package
======================================

.. warning::

   The all-order code package is legacy code and is no longer actively maintained.

The all-order code package provides corrections to the bare CI Hamiltonian arising from core-electron correlations :footcite:p:`2009AO, 2016Pb`. These corrections are stored as effective one-electron integrals in ``SGC.CON`` and two-electron screening integrals in ``SCRC.CON``, which are read by ``pconf`` when ``Ksig > 0``.

The package consists of four programs: ``allcore-ci``, ``valsd-ci``, ``sdvw-ci``, and ``second-cis``. The first three compute core, core-valence, and valence-valence all-order excitations, respectively. The fourth, ``second-cis``, computes second-order MBPT corrections over a larger set of orbitals and can be run on its own or after the all-order codes to extend their coverage.

Two levels of approximation are available:

* **CI+MBPT** — runs only ``second-cis`` to compute MBPT corrections over a large set of orbitals. This is the baseline approach when all-order corrections are not needed.

* **CI+all-order** — runs all four programs: the all-order codes first for core and core-valence correlations at high accuracy, then ``second-cis`` which reads the all-order results and extends MBPT coverage to higher orbitals. The combined output provides high-accuracy corrections where all-order is available and MBPT coverage elsewhere.

Both methods require a combined B-spline + DHF basis set described in :doc:`basis`. The typical workflow is to first run the all-order codes (``allcore-ci``, ``valsd-ci``, ``sdvw-ci``), then run ``second-cis`` last to extend MBPT coverage to higher orbitals, optionally reading the all-order results where available, as described in :doc:`all-order`.

**References**

.. footbibliography::
