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


``ine_dyn_E28`` can either solve the inhomogeneous equation iteratively by solving for a smaller matrix first, or by direct matrix inversion via the LAPACK library. It is controlled by the parameter ``IP1`` in ``conf.par``:

.. code-block:: 

    PARAMETER(IP1   =  15000,  ! Nd1    - number of determinants for direct diagonalization

This parameter can be set to be larger than the number of determinants in your problem if you don't want the program to iterate at all. The problem with iterations is that they diverge in many cases for dynamic polarizabilities. However, the problem with the direct solution is that it takes a long time to run (about 20-30 min even for 10, 000 determinants).

The program can be executed via the command:

.. code-block:: 

    ./ine_dyn_E28 <inf.ine

.. code-block:: 

    0  # start new computation
    1  # calculate polarizability of the first level
    0  # 0 for static, omega for dynamic

For dynamic polarizability, you need to run ``ine_dyn_E28`` twice: once with :math:`+\omega` and once with :math:`-\omega`. For example, if one needs the polarizability for :math:`\lambda=800 \text{ nm}`, compute :math:`\omega` in a.u.:  
:math:`\omega=1.0\times 10^7 / (au\times\lambda) = 0.056954191`, where \( au=219474.6313705 \).  
Run ``ine_dyn_E28`` twice, once with 

.. code-block:: 

    0
    1
    0.056954191 

and once with 

.. code-block:: 

    0
    1
    -0.056954191 

then average the two results. 