add - creating the configuration list
-------------------------------------

The ``add`` program constructs a list of configurations to define the CI space by exciting electrons from a set of reference configurations to a set of active non-relativistic shells. It takes in the input file ``ADD.INP``, which specifies the reference configurations, active non-relativistic shells, and minimum and maximum occupation numbers of each shell. It writes the file ``CONF.INP``, which includes a list of user-defined parameters and the list of configurations constructed by exciting electrons from a list of basic configurations.

The following is a sample input ``ADD.INP`` file. Each line has a description of the respective variable. The third block starting with ``4f  9 14`` is a list of the orbitals and minimum and maximum occupation numbers. For example, ``4f  9 14`` refers to having a minimum of 9 electrons or a maximum of 14 electrons for the 4f orbital. 

.. code-block:: 

    Ncor=  4               !# number of basic configurations. Must match the list below.
    NsvNR 16               !# number of active NR shells. The list below may be longer. 
    mult=  2               !# multiplicity of excitations. For full CI use mult=Ne   
     NE = 14               !# number of valence electrons

    L:   4f14              !# list of basic configurations
    L:   4f13  5p1         !# from which electrons are excited from.
    L:   4f12  5s2         !# the number of configurations listed here
    L:   4f11  5s2   5p1   !# must match the number on the first line 'Ncor= 4'
     nnlee nnlee nnlee  !# formatting of configurations
                           !# the numbers nn refer to the principal quantum number
                           !# the letters l refer to the angular momentum quantum number
                           !# the numbers ee refer to the occupation of that orbital
       4f  9 14   5s  0  2   5p  0  3   5d  0  2   5f  0  2   5g  0  2  
       6s  0  2   6p  0  2   6d  0  2   6f  0  2   6g  0  2   7s  0  2  
       7p  0  2   7d  0  2   7f  0  2   7g  0  2          
    ##nnl ee ee  nnl ee ee  nnl ee ee  nnl ee ee  nnl ee ee  nnl ee ee
    >>>>>>>>>>>>> Head of the file CONF.INP >>>>>>>>>>>>>>>>>>>>>>>>
      Ir17+_even            #   ion_parity                                                                                                      
      Z = 77.0              # atomic number    
     Am = 192.0             # atomic weight       
      J =  4.0              # total angular momentum      
     Jm =  4.0              # angular momentum projection      
     Nso=  14               # number of closed core shells     
     Nc =   10              # number of relativistic configurations (ignored in add program)         
     Kv =   4               # Kv = (3 - use projections, 4 - no projections)                 
     Nlv=  5                # number of energy levels    
     Ne =  14               # number of valence electrons                               
     Kl4=   1               # Kl4 = (1 - initial approx. from energy matrix, 2 - initial approx. from CONF.XIJ file)
     Nc4=  20               # number of relativistic configurations in initial approximation
     Gj =  0.0000           #                                                           
    Crt4=  0.0001           # cutoff criteria for davidson  convergence                                                          
    kout= 0                 # key for level of output (0 - low detail output, 1 - detailed output)
    Ncpt=   0               # number of relativistic configurations in PT block (ignored in add program)
    Cut0= 0.0001            # cutoff criteria for weights of PT configurations
    N_it=  100              # number of davidson iterations
    Kbrt= 1                 # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)
    Gnuc= 1.07              # gyromagnetic ratio             
         0.1002     0.2002    -0.2102     0.2104     0.3002    -0.3102
         0.3104    -0.3204     0.3206     0.4002    -0.4102     0.4104
        -0.4204     0.4206


.. note::

    The second block listing the basic configurations has a specific formatting ``__nnlee__``, where ``__`` indicate spaces, ``nn`` is the principal quantum number, ``l`` is the angular momentum quantum number as a letter (``s=0``, ``l=1``, ``d=2``, ...), and ``ee`` is the number of electrons in that orbital. 

.. note::

    The order in which the configurations and basis orbitals must be listed identically with those from ``BASS.INP``.