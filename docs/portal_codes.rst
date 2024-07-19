Portal codes
============

In this section, we describe supplementary python scripts used to generate csv-formatted databases of energy levels and matrix elements for the University of Delaware's `Portal for High-Precision Atomic Data and Computation <https://www1.udel.edu/atom/>`_. The required input files can be obtained from ``conf`` for energies and ``dtm`` for matrix elements. This process has been heavily automated using the `pCI-py scripts <pci-py.md>`_. 

Method
------
The method of generating the csv-formatted atomic data files is summarized in the following steps:

1. *Reading and reformatting input files*  
   
    First, we read all ``CONFFINAL.RES`` outputs from ``conf`` from each calculation done. The script will first check directories ``/ci+all-order/even/``, ``/ci+all-order/odd/``, ``/ci+second-order/even/`` and ``/ci+second-order/odd/`` for this file. If these directories are not detected, then it will check the ``/DATA_RAW/`` directory, where users should place them. If using this directory, the all-order output files should be named ``CONFFINALeven.RES`` and ``CONFFINALodd.RES``, and the second-order output files should be named ``CONFFINALevenMBPT.RES`` and ``CONFFINALoddMBPT.RES``. If the all-order and second-order directories are detected, the script will place them in ``/DATA_RAW/`` itself. Additionally, matrix element output files, such as ``E1.RES`` from ``dtm`` can be read. The script will check directories ``/ci+all-order/dtm/`` and ``/ci+second-order/dtm/`` first before ``/DATA_RAW/``. The all-order and second-order output files should be named, for example, ``E1.RES`` and ``E1MBPT.RES``.  
    Once all input files are read, if second-order/MBPT exist, then uncertainties are calculated and a new csv-formatted ``CONFFINAL.csv`` is written with uncertainties. Otherwise, they are all set to 0 in this file.  
    Experimental data are then acquired by parsing the `NIST Atomic Spectral Database <https://physics.nist.gov/PhysRefData/ASD/levels_form.html>`_ for full list of energy levels and stored for comparison with the data we generated through pCI. 

2. *Filtering and correcting misidentified configurations*  
   
    Sometimes the configurations and terms of energy levels outputted from the pCI codes might be be misidentified due to discrepancies in basis or numerical precision. At this stage, we attempt to correct misidentified configurations using the data from NIST. Here we create a correspondence that maps the parsed NIST identifications to the pCI theory data. This correspondence or mapping is written to the files ``/DATA_Output/Element_Even.txt`` and  ``/DATA_Output/Element_Odd.txt``. 

3. *Outputting data for portal*  
   
    The final part of the portal codes reformats the mapping of NIST and pCI energy levels for use on the UD Atom portal. The output is a csv-formatted file of the energies of the system, with a preference for NIST data over pCI-calculated data, i.e. for each configuration, NIST energies and identifications are chosen over the theory data. If NIST data is not available, then the theory values are used. The final column of the csv file ``is_from_theory`` is set to ``True`` if theory values are used. The energy levels are stored in ``Element_Energies.csv`` and the matrix elements are stored in ``Element_Matrix_Elements_Theory.csv``. Additionally, the script ``calc_lifetimes.py`` can be used to generate csv-formatted files of transition rates and lifetimes. These are outputted as ``Element_Lifetimes_Error_Check.csv`` and ``Element_Transition_Rates_Error_Check.csv``.
