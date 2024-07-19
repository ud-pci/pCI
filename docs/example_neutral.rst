Example: Neutral atoms
======================

*The following instructions assume familiarity with the main programs of the pCI package.*

In this section, we describe a method used to construct basis sets for the case of neutral atoms. In these example, we utilize a single ``HFD.INP`` to construct the orbitals. 

Ac
--

.. code-block::

      Ac            
       KL =   0       # (0 - new calculation, 1 - continue)  
       NS =  36       # number of orbitals     
       NSO=  24       # number of closed orbitals (in this case only 1s2, 2s2)     
       Z  =  89.0     # atomic number       
       AM = 227.00    # atomic mass        
       JM =  -2.0     # (-2 - average for non-relativistic configuration)       
       R2 =  60.0     # radius of cavity       
       kbr= 2         # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)   
      rnuc= 5.7350    # (optional) rms nuclear radius (https://www-nds.iaea.org/radii/)   

              NL   J       QQ     KP   NC

        1     1S (1/2)   2.0000    0    0
        2     2S (1/2)   2.0000    0    0
        3     2P (1/2)   2.0000    0    0
        4     2P (3/2)   4.0000    0    0
        5     3S (1/2)   2.0000    0    0
        6     3P (1/2)   2.0000    0    0
        7     3P (3/2)   4.0000    0    0
        8     3D (3/2)   4.0000    0    0
        9     3D (5/2)   6.0000    0    0
       10     4S (1/2)   2.0000    0    0
       11     4P (1/2)   2.0000    0    0
       12     4P (3/2)   4.0000    0    0
       13     4D (3/2)   4.0000    0    0
       14     4D (5/2)   6.0000    0    0
       15     4F (5/2)   6.0000    0    0
       16     4F (7/2)   8.0000    0    0
       17     5S (1/2)   2.0000    0    0
       18     5P (1/2)   2.0000    0    0
       19     5P (3/2)   4.0000    0    0
       20     5D (3/2)   4.0000    0    0
       21     5D (5/2)   6.0000    0    0
       22     6S (1/2)   2.0000    0    0
       23     6P (1/2)   2.0000    0    0
       24     6P (3/2)   4.0000    0    0
       25     7S (1/2)   1.0000    0    1
       26     6D (3/2)   1.0000    0    2
       27     6D (5/2)   0.0000    0    2
       28     7P (1/2)   1.0000    0    3
       29     7P (3/2)   0.0000    0    3
       30     5F (5/2)   1.0000    0    4
       31     5F (7/2)   0.0000    0    4
       32     8S (1/2)   1.0000    0    5
       33     7D (3/2)   1.0000    0    6
       34     7D (5/2)   0.0000    0    6
       35     8P (1/2)   1.0000    0    7
       36     8P (3/2)   0.0000    0    7
      

Sr
--

.. code-block:: 

      Sr III
       KL =   0       # (0 - new calculation, 1 - continue)             
       NS =  20       # number of orbitals         
       NSO=  12       # number of closed orbitals (in this case only 1s2, 2s2)        
       Z  =  38.0     # atomic number             
       AM =  90.000   # atomic mass                
       JM =  -2.0     # (-2 - average for non-relativistic configuration)             
       R2 =  60.0     # radius of cavity             
       kbr= 0         # key for Breit (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)     

              NL   J       QQ     KP   NC

        1     1S (1/2)   2.0000    0    0
        2     2S (1/2)   2.0000    0    0
        3     2P (1/2)   2.0000    0    0
        4     2P (3/2)   4.0000    0    0
        5     3S (1/2)   2.0000    0    0
        6     3P (1/2)   2.0000    0    0
        7     3P (3/2)   4.0000    0    0
        8     3D (3/2)   4.0000    0    0
        9     3D (5/2)   6.0000    0    0
       10     4S (1/2)   2.0000    0    0
       11     4P (1/2)   2.0000    0    0
       12     4P (3/2)   4.0000    0    0
       13     5S (1/2)   1.0000    0    1
       14     5P (1/2)   1.0000    0    2
       15     5P (3/2)   0.0000    0    2
       16     4D (3/2)   1.0000    0    3
       17     4D (5/2)   0.0000    0    3
       18     6S (1/2)   1.0000    0    4
       19     6P (1/2)   1.0000    0    5
       20     6P (3/2)   0.0000    0    5
