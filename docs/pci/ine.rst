ine - polarizabilities
----------------------

The ``ine`` program calculates static and dynamic polarizabilities of specified atomic levels. ``ine`` only gives the valence polarizability. Core polarizability needs to be computed separately with a different code. The code will give both scalar and tensor polarizabilities if the tensor polarizability is not zero, but not the vector polarizability. There are several version of the code, but for now we will use ``ine_dyn_E28``. This program requires several input files from previously ran ``conf`` and ``dtm`` programs, including ``CONF.DET`` and ``CONF.XIJ`` of the parity of the level of interest (renamed to ``CONF0.DET`` and ``CONF0.XIJ``), ``CONF.INP``, ``CONF.XIJ``, ``CONF.HIJ``, and ``CONF.JJJ`` of the opposite parity, and the file ``DTM.INT`` from ``dtm``. 

For example, if we want to calculate polarizabilities for an even state:

.. code-block:: 

    cp CONFeven.DET CONF0.DET
    cp CONFeven.XIJ CONF0.XIJ
    
    cp CONFodd.INP CONF.INP
    cp CONFodd.XIJ CONF.XIJ
    cp CONFodd.HIJ CONF.HIJ
    cp CONFodd.JJJ CONF.JJJ

Running ine
~~~~~~~~~~~

Before parallel ``ine`` can be run, the parallel matrix files ``CONFp.HIJ`` and ``CONFp.JJJ`` have to be sorted. This can be done using the python script ``/pCI/src/auxiliary/sort.py``. 

1. First load the python3
   
    .. code-block:: 

        vpkg_require python/3
    

2. Run python script ``sort.py``
   
    .. code-block:: 

        python3 sort.py
    

3. The program will ask which file you want sorted, so enter ``CONFp.HIJ`` or ``CONFp.JJJ``
   
    .. code-block:: 

        CONFp.HIJ


4. Now re-run the program for the other matrix
   
    .. code-block:: 

        python3 sort.py  
        CONFp.JJJ

5. After these files are obtained, the input file ``ine.in`` must be constructed:

    .. code-block:: 

        0                          | kIters= (0- invert and iterate if diverged, 1-invert only, 2-2-step iteration)
        0                          | Kl = (0-new, 1-use X1, 2-use X1,Y1,Y2)
        2                          | Kli = (1 - H_p, 2 - E1(L), 5 - E2) for RHS of equation
        2                          | Klf = (1 - H_p, 2 - E1(L), 3 - H_am, 4 - E1(V), 5 - E2) for LHS of equation
        5                          | N0 = record number of X0
        5                          | N2 = record number of X2
        3                          | nlambda - number of ranges of wavelengths to include (in this case, 3)
        813.01 813.02 0.01         | range 1   (lambda1, lambda2, step_size)
        813.02 813.025 0.001       | range 2
        813.025 813.03 0.01        | range 3

6. Finally, ``ine`` can be run:

    .. code-block::

        ./ine

