Supplementary programs
======================

In this section, we describe supplementary programs used to generate basis sets or optimized configuration lists.

bdhfA - Dirac-Hartree-Fock code
-------------------------------

``bdhfA`` can be used in place of ``bdhf`` to make a B-spline basis set. This program reads in an input file ``xx.dat`` with a much simpler input than for ``bdhf``, and will write ``bdhf.in``, which can be renamed to ``bas_wj.in`` for subsequent use in all-order codes. 

.. raw:: html

    <details>
    <summary>
        Click here to see a description of a <a href="../src/cs_min.dat" download>minimum input file for Cs</a>.
    </summary>
    
.. code-block:: 

    1
    1
    Cs 
    1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 
    ==========
    7 70 9 220.0

.. raw:: html

    </details>  
    <details>
    <summary>
        Click here to see a description of a <a href="../src/cs_min.dat" download>complete input file for Cs</a>.
    </summary>


.. code-block:: 

    1
    1
    Cs 137 0.0001
    1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 
    6 9 s1/2
    6 9 p1/2
    5d3/2
    5d5/2
    ==========
    7 70 9 220.0

.. raw:: html

    </details>  
    <details>
    <summary>
        Click here to see a description of a <a href="../src/cs_min.dat" download>complete input file for Yb</a>.
    </summary>

.. code-block:: 

    1
    1
    Yb 176 0.0001
    1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f
    6s1/2
    5d3/2
    5d5/2
    ==========
    7 50 9 60.0

.. raw:: html

    </details>
    <br>

The input files have the following format:

.. code-block::

    1                   # 0 for old input, 1 for new input
    1                   # convergence parameter: 1 - most cases, 5 - in case of convergence problems
    Cs                  # name of element - e.g. Cs, Yb, etc. nuclear charge Z is set using this name
                        # Isotope can also be entered. If no isotope is entered, most abundant or 
                        # most stable isotope will be used. 
                        # First grid point can be entered if finer grid is used inside the nuclei.
                        # The default is set to 0.0005, which will produce a few points inside the nucleus.
                        # Use 0.0001, for example, for more points (~25).
                        # Examples (either will work):
                        # Cs
                        # Cs 137
                        # Cs 0.0001
    1s 2s ... 5p        # List of filled core shells. This code has to use fully filled subshells for core potential.
                        # Valence electron lines can be omitted if the code is used for bspl.
                        # To include valence electrons, input nmin nmax l j, or nlj:
                        # 6 8 s1/2   # for 6, 7, 8 s1/2
                        # 6 s1/2     # for just 6s1/2
    ==========          # required line for formatting
    7 70 9 220.0        # highest partial wave, number of splines, order of splines, and cavity radius in a.u.


con_cut - truncating configuration lists
----------------------------------------

``con_cut`` is used to truncate a configuration list to only hold all configurations with weights above a user-specified cutoff. This program outputs the file ``CONF_CUT.RES``, which can be used as a new ``CONF.INP`` file. The configurations in ``CONF_CUT.RES`` are also listed in order of descending weights, so the most important configurations are at the top of the list. You can copy the top configurations to your ``ADD.INP`` file and reconstruct a new configuration list with the top configurations as basic configurations, or use ``merge_ci`` to put the top configurations in another ``CONF.INP`` file.

merge_ci - merging configuration lists 
--------------------------------------

``merge_ci`` takes in two ``CONF.INP`` files named ``C_A.INP`` and ``CONF.INP`` and outputs the file ``C_M.INP``. This output file can be renamed ``CONF.INP`` again to be used in ``conf``. 

In summary:

1. Run ``con_cut``, inputting a cutoff threshold for the weight, to obtain ``CON_CUT.RES``.
2. Rename ``CON_CUT.RES`` to ``C_A.INP``.
3. Replace the basic configurations in ``ADD.INP`` with the top configurations from ``C_A.INP`` and run ``add``.
4. Run ``merge_ci``, combining the configurations in ``C_A.INP`` and ``CONF.INP`` to obtain ``C_M.INP``.
5. Rename ``CONF.INP`` to ``C_B.INP``, and ``C_M.INP`` to ``CONF.INP``.
6. Run ``conf``. 

conf_pt - valence perturbation theory
-------------------------------------

When running CI calculations with a large number of configurations relative to the number of work resources, it is often times necessary to determine the most important configurations in the CI space, and truncate the list of configurations to make successive calculations feasible. 

The ``conf_pt`` program begins the same way ``conf`` begins, reading in several input parameters and the list of configurations from the file ``CONF.INP``, the basis set from the files ``HFD.DAT`` and ``CONF.DAT``, and the radial integrals from ``CONF.INT`` and ``CONF.GNT``. In addition, it also reads the CI eigenvectors from the file ``CONF.XIJ``. 

qed_pot_conf - QED potentials
-----------------------------

The pCI package includes 5 variants (``kvar``) of QED potentials:

1. QEDMOD
2. Flambaum local potential + QEDMOD non-local correction
3. Flambaum local potential
4. QEDPOT
5. Semi-empirical approach

sort.py - converts parallel matrix element files to serial format
-----------------------------------------------------------------

``sort.py`` is a python program that sorts the matrix elements of the operator :math:`J^2` in order of ascending index (as done in the serial version of the ``conf`` program). This program takes in the parallel file ``CONF.JJJ`` or ``CONF.HIJ`` and returns a serial version of the inputted file. However, there is one change made to the file. There is an additional integer at the start of the file that stores the total number of matrix elements. 