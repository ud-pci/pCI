Fe XVII and Ni XIX
==================

*The following instructions assume familiarity with the main programs of the pCI package.*

In this section, we describe a method used to construct basis sets for the cases of Fe XVII and Ni XIX. In this example, we run the ``hfd`` program in 3 sequential stages to construct the DF orbitals, then run the ``bass`` program to form virtual orbitals to account for correlations. The following is a list of the input files used in this example.

1. ``h_m_1.inp`` - In the first stage, we construct the :math:`1s, 2s, 2p, 3s, 3p, 3d` orbitals with the :math:`1s^2 2s^2 2p^5 3d` configuration. ``QQ`` is ``0`` for :math:`3s`` and :math:`3p`, but the orbital is still formed. DF is solved with 0 occupation number for these orbitals, but they won't be very good orbitals. We can think of :math:`3s0` and :math:`3p0` as placeholders to keep the order of orbitals, and re-construct them in the next step.

2. ``h_m_2.inp`` - In the second stage, we freeze the :math:`1s, 2s, 2p, 3d` orbitals and re-construct the :math:`3s, 3p` orbitals from the :math:`1s^2 2s^2 2p^5 3s` and :math:`1s^2 2s^2 2p^5 3p` configurations. Note that ``NS`` is unchanged here, since we do not add any additional orbitals.

    The following code block shows ``h_m_1.inp`` on the left and ``h_m_2.inp`` on the right of the partition. The head of both input files are identical. 

    .. code-block:: 

        Fe XVII & Ni XIX
        KL =   0                # (0 - new calculation, 1 - continue)
        NS =   9                # number of orbitals
        NSO=   2                # number of closed orbitals (in this case only 1s2, 2s2)
        Z  =  26.0              # atomic number
        AM =  56.000            # atomic mass
        JM =  -2.0              # default parameter (do not change)
        R2 =  20.0              # radius of cavity
        kbr= 2                  # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)
        
               NL   J       QQ     KP   NC          |                NL   J       QQ     KP   NC
                                                    |
         1     1S (1/2)   2.0000    0    0          |          1     1S (1/2)   2.0000    1    0
         2     2S (1/2)   2.0000    0    0          |          2     2S (1/2)   2.0000    1    0
         3     2P (1/2)   2.0000    0    0          |          3     2P (1/2)   2.0000    1    0
         4     2P (3/2)   3.0000    0    0          |          4     2P (3/2)   3.0000    1    0
         5     3S (1/2)   0.0000    0    0          |          5     3S (1/2)   1.0000    0    1
         6     3P (1/2)   0.0000    0    0          |          6     3P (1/2)   1.0000    0    2
         7     3P (3/2)   0.0000    0    0          |          7     3P (3/2)   0.0000    0    2
         8     3D (3/2)   1.0000    0    0          |          8     3D (3/2)   0.0000    1    3
         9     3D (5/2)   0.0000    0    0          |          9     3D (5/2)   0.0000    1    3
        
        
3. ``h_m_3.inp`` - At the third stage, we freeze the :math:`1s, 2s, 2p, 3s, 3p, 3d` orbitals, then construct the :math:`4s, 4p, 4d, 4f, 5g` orbitals from the :math:`1s^2 2s^2 2p^5 4s, 1s^2 2s^2 2p^5 4p, \dots, 1s^2 2s^2 2p^5 5g` configurations. Note that the number of orbitals ``Ns`` has changed from 9 to 18, since we are adding 9 additional orbitals. 

    .. code-block:: 

        Fe XVII & Ni XIX
        KL =   0                #    
        NS =  18                # number of orbitals   
        NSO=   2                # number of closed orbitals (in this case only 1s2, 2s2)   
        Z  =  26.0              # atomic number     
        AM =  56.000            # atomic mass       
        JM =  -2.0              #     
        R2 =  20.0              # radius of cavity     
        kbr= 2                  # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit) 

               NL   J       QQ     KP   NC

         1     1S (1/2)   2.0000    1    0
         2     2S (1/2)   2.0000    1    0
         3     2P (1/2)   2.0000    1    0
         4     2P (3/2)   3.0000    1    0
         5     3S (1/2)   0.0000    1    0
         6     3P (1/2)   0.0000    1    0
         7     3P (3/2)   0.0000    1    0
         8     3D (3/2)   0.0000    1    0
         9     3D (5/2)   0.0000    1    0
        10     4S (1/2)   1.0000    0    1
        11     4P (1/2)   1.0000    0    2
        12     4P (3/2)   0.0000    0    2
        13     4D (3/2)   1.0000    0    3
        14     4D (5/2)   0.0000    0    3     
        15     4F (5/2)   1.0000    0    4
        16     4F (7/2)   0.0000    0    4
        17     5G (7/2)   1.0000    0    5
        18     5G (9/2)   0.0000    0    5


4. ``b_m_2.inp`` - Finally, we run the ``bass`` program to add virtual orbitals to ``HFD.DAT`` to account for correlations. Here, we are constructing the :math:`24spdfg` basis set, where the designation :math:`24spdfg` means that all orbitals up to :math:`n=24` are included for the :math:`spdfg` partial waves. Therefore, the list of orbitals included here extends to :math:`24g`, designated by ``-2.4401`` and ``2.4401``. The digital format here is represented as ``sn.nlqq``, where ``nn`` represents the principal quantum number, ``l`` is the orbital angular momentum quantum number, ``l`` is the orbital angular momentum quantum number, and ``qq`` is the occupation number of the orbital. The sign ``s`` corresponds to the total angular momentum, represeented by ``-`` for :math:`j=l-1/2`, or an empty space for :math:`j=l+1/2`. For brevity, not all orbitals are displayed. 

    .. code-block:: 

         Fe XVII & Ni XIX
          Z = 26.0
         Am = 52.0
         Nso=    4         # number of core orbitals (defines DF operator)
         Nv =  194         # number of valence & virtual orbitals
         Ksg=    1         # defines Hamiltonian: 1-DF, 3-DF+Breit
         Kdg=    0         # diagonalization of Hamiltonian (0=no,1,2=yes)
         orb= 4s 1         # first orbital for diagonalization
         Kkin    1         # kinetic balance (0,1,or 2)
         orb= 5s 1         # first orbital to apply kin.bal.
         orb= 2p 3         # last frozen orbital
         orb= 0p 3         # last orbital in basis set
        kout= 0            # detail rate in the output
        kbrt= 2            # 0,1,2 - Coulomb, Gaunt, Breit
        ----------------------------------------------------------
             0.1002     0.2002    -0.2102     0.2104
        
          1  0.3001               # 
          2 -0.3101               # These orbitals are in HFD.DAT already run by hfd
          3  0.3101               # 
          4 -0.3201               # 
          5  0.3201               # 
                              
          6  0.4001  3  0.4001    # reading 4s from 4s from HFD.DAT
          7 -0.4101  3 -0.4101    # key '3' means 'read in from HFD.DAT'
          8  0.4101  3  0.4101    # HFD.DAT is h_m_3.inp in this case
          9 -0.4201  3 -0.4201    
         10  0.4201  3  0.4201    
         11 -0.4301  3 -0.4301
         12  0.4301  3  0.4301
                              
         13  0.5001               # key '0' or ' ' means 'build nl from (n-1)l'
         14 -0.5101               # e.g. 5s is built from 4s, 5p from 4p
         15  0.5101               #      5d from 4d, ...
         16 -0.5201    
         17  0.5201     
         18 -0.5301    
         19  0.5301     
         20 -0.5401  3 -0.5401    # since key '3' is present, 5f is read in from HFD.DAT
         21  0.5401  3  0.5401 
         22 -0.6401    
         23  0.6401     
         :
         :
         :
        193 -2.4401
        194  2.4401

The following bash script utilizes the above input files and forms the final :math:`24spdfg` basis set for Fe XVII and Ni XIX.

.. code-block:: 

    #! /bin/bash 
    #####################################################################
    # script to form basis set for Fe 16+ and Ni 18+
    cp h_m_1.inp HFD.INP    
    ./hfd                   
    cp h_m_2.inp HFD.INP     
    ./hfd                     
    cp HFD.DAT h0.dat          
    cp h_m_3.inp HFD.INP       
    ./hfd                       
    mv HFD.DAT h_m.dat        
    mv h0.dat HFD.DAT          
    cp b_m_2.inp BASS.INP    
    ./bass <b.in              
    ./bass                    
    echo "    End of script"
    
