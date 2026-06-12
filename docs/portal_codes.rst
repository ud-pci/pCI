Portal codes
============

In this section, we describe supplementary python scripts used to generate csv-formatted databases of energy levels and matrix elements for the University of Delaware's `Portal for High-Precision Atomic Data and Computation <https://www1.udel.edu/atom/>`_. The required input files can be obtained from ``pconf`` for energies and ``pdtm`` for matrix elements. This process has been automated using the :doc:`pCI-py scripts <pci-py>`. 

Method
------
The method of generating the csv-formatted atomic data files is summarized in the following steps:

1. *Reading and reformatting input files*

    First, we read all ``FINAL.RES`` outputs from ``pconf`` from each calculation done. The script will first check the ``ci+all-order/`` and ``ci+second-order/`` top-level directories for ``even*`` and ``odd*`` parity subdirectories (e.g. ``even0``, ``odd1``). If a config.yml is provided the J values in ``conf.even.J`` and ``conf.odd.J`` determine which subdirectory is selected automatically; otherwise the script prompts interactively. If these top-level directories are not detected, then it will check the ``DATA_RAW/`` directory, where users should place the files manually. The all-order output files should be named ``CONFFINALeven.RES`` and ``CONFFINALodd.RES``, and the second-order output files should be named ``CONFFINALevenMBPT.RES`` and ``CONFFINALoddMBPT.RES``. If the all-order and second-order directories are detected, the script will copy and rename the files into ``DATA_RAW/`` itself. Additionally, matrix element output files such as ``E1.RES`` from ``pdtm`` can be read. The script will check ``tm*`` subdirectories within ``ci+all-order/`` and ``ci+second-order/`` first before ``DATA_RAW/``. If two ``tm*`` directories are found (e.g. for two J values), their ``E1.RES`` files are merged into a single ``E1.RES``; the second-order counterpart is handled identically to produce ``E1MBPT.RES``.
    Once all input files are read, if second-order/MBPT data exist, energy uncertainties are estimated as the absolute difference between the CI+all-order and CI+second-order energies and a summary file ``final_res.csv`` is written. Otherwise uncertainties are set to 0.
    Experimental data are then acquired by parsing the `NIST Atomic Spectral Database <https://physics.nist.gov/PhysRefData/ASD/levels_form.html>`_ for the full list of energy levels and stored for comparison with the data generated through pCI.

2. *Filtering and correcting misidentified configurations*

    Sometimes the configurations and terms of energy levels outputted from the pCI codes are misidentified due to strong configuration mixing or numerical precision. At this stage, we attempt to correct misidentified configurations and establish a one-to-one correspondence between NIST and theory levels. The process proceeds in several sub-steps:

    a. *Preprocessing.* Both the NIST and theory datasets are loaded. A common ground-state core prefix is subtracted from configuration strings so that comparisons are made on the valence part only. Theory levels are filtered to only those flagged as converged. For NIST levels that have no term label (e.g. levels with :math:`jj`-coupled designations), a term is inferred from the best-matching theory state.

    b. *Configuration correction.* The theory sometimes assigns the wrong principal quantum number to an orbital (e.g., predicts ``2s.7p`` when the physical level is ``2s.6p``). The corrected configurations are found by tracking the expected sequence of principal quantum numbers for each :math:`(L, 2S+1)` group and detecting gaps where the theory has jumped ahead by one.

    c. *NIST-theory matching.* For each NIST level, candidate theory counterparts are found through a series of passes:

       - **Pass 1**: :math:`J` must match.
       - **Pass 2**: Configuration must match the primary theory configuration (or its inverted orbital ordering).
       - **Pass 3**: Term symbol must match exactly; if not, a correction of ±1 to the multiplicity is tried (Pass 4).
       - **Pass 5**: If the term matches but configurations differ, the theory principal quantum number of the last orbital is decremented by one and re-checked against the NIST configuration.

       When multiple candidates survive, they are ranked by the percentage energy difference :math:`\Delta E / E_\mathrm{NIST}`. Secondary-configuration matches receive a 2% additive penalty to prefer primary-configuration matches. Assignments are then made greedily in order of score, with each theory level allowed to match at most one NIST level. Theory levels that remain unmatched receive their corrected (sequential) configuration label, provided it does not conflict with an already-assigned NIST designation.

    d. *Missing level detection.* After the mapping is complete, the expected set of term symbols for each configuration is computed from the electron occupation and compared against the combined NIST+theory set. Any term/J combinations that are theoretically expected but absent from both datasets are appended to the output and flagged with a ``~`` prefix.

    The correspondence is written to ``DATA_Output/{Element}_Even.txt`` and ``DATA_Output/{Element}_Odd.txt``, with a companion ``+missing`` file listing the absent levels.

3. *Outputting data for portal*

    The final part of the portal codes reformats the mapping of NIST and pCI energy levels for use on the UD Atom portal. The output is a csv-formatted file of the energies of the system, with a preference for NIST data over pCI-calculated data, i.e. for each configuration, NIST energies and identifications are chosen over the theory data. If NIST data is not available, then the theory values are used. The final column of the csv file ``is_from_theory`` is set to ``True`` if theory values are used. The energy levels are stored in ``Element_Energies.csv`` and the matrix elements are stored in ``Element_Matrix_Elements_Theory.csv``. Matrix element uncertainties are combined in quadrature from the theory uncertainty (difference between all-order and second-order matrix elements) and a minimum floor based on the ``min_uncertainty`` parameter. The total matrix element uncertainty is given by:

    .. math::

        \Delta_\mathrm{total} = \sqrt{\Delta_\mathrm{theory}^2 + \left(\frac{f_\mathrm{min}}{100} \cdot |d|\right)^2}

    where :math:`\Delta_\mathrm{theory}` is the uncertainty from comparing CI+all-order and CI+second-order matrix elements, :math:`f_\mathrm{min}` is the minimum uncertainty percentage (``portal.min_uncertainty`` in the config), and :math:`d` is the matrix element value. Additionally, the script ``gen_rates_lifetimes.py`` can be used to generate csv-formatted files of transition rates and lifetimes from these two output files. These are outputted as ``Element_Lifetimes_Error_Check.csv`` and ``Element_Transition_Rates_Error_Check.csv``.

Programs
--------

gen_portal_csv.py
~~~~~~~~~~~~~~~~~

Reads pCI energy and matrix element output files, compares them against the NIST Atomic Spectra Database, corrects misidentified configurations, and writes the final csv files used by the portal.

**Usage**

.. code-block:: bash

    python gen_portal_csv.py

The script is interactive. It first asks whether a ``config.yml`` file will be used. If yes, it prompts for the filename; if no, it prompts for the atom name and minimum matrix element uncertainty.

**Configuration (config.yml)**

When a config.yml is provided, the following keys are read:

- ``atom.name`` — element and ionization stage in any of the formats ``Ba``, ``Ba I``, ``Ba+``, or ``Ba2``
- ``conf.even.J``, ``conf.odd.J`` — J values used to locate the ``even{J}`` and ``odd{J}`` output directories
- ``portal.ignore_g`` (optional, default ``True``) — exclude states with a ``g`` orbital from the portal output
- ``portal.min_uncertainty`` (optional) — minimum matrix element uncertainty as a percentage of the value; prompted interactively if absent
- ``portal.min_energy_diff_percent`` (optional, default ``3.0``) — levels whose NIST-theory energy difference exceeds this percentage are excluded
- ``portal.energy_cutoff`` (optional) — explicit global energy cutoff in cm\ :sup:`-1`; overrides the ``min_energy_diff_percent`` cutoff when set

**Input files**

The script first checks for ``ci+all-order/`` and ``ci+second-order/`` subdirectories and copies their output into ``DATA_RAW/`` automatically. If those directories are absent, files must be placed in ``DATA_RAW/`` manually with the following names:

+----------------------------------+-------------------------------------------------+
| File                             | Description                                     |
+==================================+=================================================+
| ``CONFFINALeven.RES``            | CI+all-order even parity energies               |
+----------------------------------+-------------------------------------------------+
| ``CONFFINALodd.RES``             | CI+all-order odd parity energies                |
+----------------------------------+-------------------------------------------------+
| ``CONFFINALevenMBPT.RES``        | CI+second-order even parity energies            |
+----------------------------------+-------------------------------------------------+
| ``CONFFINALoddMBPT.RES``         | CI+second-order odd parity energies             |
+----------------------------------+-------------------------------------------------+
| ``E1a.RES``, ``E1b.RES``         | E1 matrix elements (combined into E1.RES)       |
+----------------------------------+-------------------------------------------------+
| ``E1MBPTa.RES``, ``E1MBPTb.RES`` | MBPT matrix elements (combined into E1MBPT.RES) |
+----------------------------------+-------------------------------------------------+

The second-order files are optional; if absent, uncertainties are set to 0.

**Output files**

+-----------------------------------------------+---------------------------------------------------+
| File                                          | Description                                       |
+===============================================+===================================================+
| ``{Element}_Energies.csv``                    | Energy levels preferring NIST over theory         |
+-----------------------------------------------+---------------------------------------------------+
| ``{Element}_Matrix_Elements_Theory.csv``      | E1 matrix elements with uncertainties             |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Processed/FINALeven.RES``              | Even parity energies with uncertainty column      |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Processed/FINALodd.RES``               | Odd parity energies with uncertainty column       |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Output/{Element}_Even.txt``            | NIST-theory mapping table, even parity            |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Output/{Element}_Odd.txt``             | NIST-theory mapping table, odd parity             |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Filtered/UD/{Element}_UD_Even.csv``    | Filtered theory energies, even parity             |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Filtered/UD/{Element}_UD_Odd.csv``     | Filtered theory energies, odd parity              |
+-----------------------------------------------+---------------------------------------------------+
| ``DATA_Filtered/NIST/{Element}_NIST_*.csv``   | Filtered NIST energies by parity                  |
+-----------------------------------------------+---------------------------------------------------+
| ``final_res.csv``                             | Merged even+odd energies with uncertainties       |
+-----------------------------------------------+---------------------------------------------------+

gen_rates_lifetimes.py
~~~~~~~~~~~~~~~~~~~~~~

Reads the energy and matrix element csv files produced by ``gen_portal_csv.py`` and computes transition rates, lifetimes, and branching ratios in three steps:

1. Calculate E1 transition rates from matrix elements and energy differences.
2. Sum rates per upper state to obtain lifetimes and branching ratios.
3. Format all quantities with parenthetical uncertainty notation for portal display.

**Usage**

.. code-block:: bash

    python gen_rates_lifetimes.py <atom_name>

where ``<atom_name>`` is the filename prefix used by ``gen_portal_csv.py`` (e.g. ``Ba1``, ``Ca2``).

**Input files**

+--------------------------------------------+---------------------------------------------------+
| File                                       | Description                                       |
+============================================+===================================================+
| ``{atom}_Energies.csv``                    | Required. Output of ``gen_portal_csv.py``         |
+--------------------------------------------+---------------------------------------------------+
| ``{atom}_Matrix_Elements_Theory.csv``      | Required. Output of ``gen_portal_csv.py``         |
+--------------------------------------------+---------------------------------------------------+
| ``{atom}_Lifetimes_Exp.csv``               | Optional. Experimental lifetimes to replace       |
|                                            | theory values. Columns: ``state_configuration``,  |
|                                            | ``state_term``, ``state_J``,                      |
|                                            | ``lifetime_display``, ``lifetime_ref``            |
+--------------------------------------------+---------------------------------------------------+

**Output files**

+------------------------------------------------+--------------------------------------------------+
| File                                           | Description                                      |
+================================================+==================================================+
| ``{atom}_Transition_Rates.csv``                | Step 1: raw rates per matrix element row         |
+------------------------------------------------+--------------------------------------------------+
| ``{atom}_Transition_Properties.csv``           | Step 2: upper to lower transitions with display  |
|                                                | columns for wavelength, matrix element,          |
|                                                | branching ratio, and rate                        |
+------------------------------------------------+--------------------------------------------------+
| ``{atom}_Lifetimes.csv``                       | Step 2: per-state lifetimes and uncertainties    |
+------------------------------------------------+--------------------------------------------------+
| ``{atom}_Transition_Rates_Error_Check.csv``    | Step 3: transition properties with parenthetical |
|                                                | uncertainty notation                             |
+------------------------------------------------+--------------------------------------------------+
| ``{atom}_Lifetimes_Error_Check.csv``           | Step 3: lifetimes with parenthetical uncertainty |
|                                                | notation; experimental values substituted where  |
|                                                | ``{atom}_Lifetimes_Exp.csv`` is provided         |
+------------------------------------------------+--------------------------------------------------+
| ``transitions.txt``                            | Human-readable decay summary per upper state     |
+------------------------------------------------+--------------------------------------------------+
