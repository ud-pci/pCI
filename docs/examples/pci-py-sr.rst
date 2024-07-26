Sr with pCI-py
==============

In this example, we will utilize the pCI-py scripts detailed :doc:`here <../pci-py>` to generate atomic data for neutral Sr. 

config.yml
----------
The ``config.yml`` file defines the entire calculation from beginning to end. You can find a detailed explanation of the contents of this configuration file :doc:`here <../pci-py>`. 

.. raw:: html

    <details>
    <summary>
        Click here to see the contents of the sample <a href="../src/config.yml" download>config.yml</a> we will be using.
    </summary>

.. code-block:: 

    # General parameters
    system:
        name: Sr
        isotope: 
        include_breit: True

    ## Parameters used by basis programs
    basis:
        cavity_radius: 60
        orbitals:
            core: 1s 2s 2p 3s 3p 3d 4s 4p 
            valence: 5s 5p 4d 6s 6p 5d 4f
            nmax: 35
            lmax: 6
        b_splines:
            nmax: 40
            lmax: 6
            k: 7
        val_aov:
            s: 4
            p: 4
            d: 4
            f: 3
        val_energies:
            kval: 1
            energies: 
                s: -0.28000
                p: [-0.22000, -0.22000]
                d: [-0.31000, -0.31000]
                f: [-0.13000, -0.13000]

.. raw:: html

    </details>

|

Parameters used by add program
------------------------------

.. code-block:: 

    add:
        # Lists of even and odd parity reference configurations
        ref_configs:
            odd: [5s1   5p1]
            even: [5s2]
        basis_set: 22spdfg
        orbitals:
            core: 1s 2s 2p 3s 3p 3d 4s 4p 
            active: [
                4-7p:  0  4,
                4-7d:  0  4,
                4-7f:  0  4,
                5-7g:  0  4,
                ]    
        excitations:
            single: True
            double: True
            triple: False

Parameters used by conf program
-------------------------------

.. code-block:: 

    conf:
        J: 0.0
        JM: 0.0
        J_selection: False
        num_energy_levels: 24
        num_dvdsn_iterations: 100
        include_lsj: True
        write_hij: False

Parameters used by dtm program
------------------------------

.. code-block:: 

    dtm:
        matrix_elements: E1
        include_rpa: True

Parameters used by portal script
--------------------------------

.. code-block:: 

    portal:
        ignore_g: True
    
Optional parameters
-------------------

.. code-block:: 

    optional:
        qed:
            include: False
            rotate_basis: False

        isotope_shifts: 
            include: False
            K_is: 1
            C_is: 0.01

        code_method: [ci+all-order, ci+second-order]
        generate_directories: False
        run_ao_codes: True
        run_ci: True
        pci_version: default


Running pCI-py scripts
----------------------
1. Write config.yml to root directory.
   
    * Make sure to re-define ``optional.pci_version = default`` to the correct version if parameters such as the radial grid size has to be changed. On UD computers, one can use the command ``vpkg_versions pci`` to list all version of pCI.
    * Make sure to set ``optional.run_ao_codes = True`` in order for job scripts to be automatically submitted to SLURM for all-order codes to run.

2. Run basis.py 
    
    * This script will generate basis sets for CI+all-order and CI+MBPT in their respective directories.
  
    * Input:

        * ``config.yml`` (specifically blocks ``system``, ``basis``, ``optional``)

    * Output:

        * ``/CI+all-order/basis/``
        * ``/CI+second-order/basis/``

3. Run add.py
   
    * Make sure to set ``optional.run_ci = True`` in order for job scripts to be automatically submitted to SLURM for ``conf`` to run.
  
    * Input:

        * ``config.yml`` (specifically blocks ``system``, ``basis``, ``conf``, ``optional``)
        * ``/CI+all-order/basis/BASS.INP``
        * ``/CI+second-order/basis/BASS.INP``
  
    * Output:

        * ``/CI+all-order/even/CONFFINAL.RES``
        * ``/CI+all-order/odd/CONFFINAL.RES``
        * ``/CI+second-order/even/CONFFINAL.RES``
        * ``/CI+second-order/odd/CONFFINAL.RES``

4. Run dtm.py
   
    * Make sure to set ``dtm.matrix_elements = E1`` so ``E1.RES``, which contains a table of E1 transitions, is formed.
    * Make sure to set ``dtm.include_rpa = True`` to include RPA corrections.

    * Input:
  
        * ``config.yml`` (specifically blocks ``system``, ``conf``, ``dtm``)
        * ``/CI+all-order/even/``, ``/CI+all-order/odd/``, ``/CI+second-order/even/``, ``/CI+second-order/odd/``
        
            * ``CONF.INP``
            * ``CONF.DAT``
            * ``CONF.DET``
            * ``CONF.XIJ``
            * ``CONFSTR.RES``

    * Output:
  
        * ``/CI+all-order/dtm/``, ``/CI+second-order/dtm/``
  
            * ``TM.RES``
            * ``E1.RES``
 
5. Run gen_portal_csv.py
   
    * Make sure to set ``portal.ignore_g = True`` if configurations with :math:`g` orbitals or terms with :math:`G` are to be ignored in final outputs.
    * If ``CI+second-order`` directory is not found, uncertainties will be set to 0.
    
    * Input:
  
        * ``/CI+all-order/even/``, ``/CI+all-order/odd/``, ``/CI+second-order/even/``,``/CI+second-order/odd/``

            * ``CONFFINAL.RES``
  
        * ``/CI+all-order/dtm/``, ``/CI+second-order/dtm/`` (optional)
  
            * ``E1.RES``
            * ``E1MBPT.RES``
            
    * Output:

        * ``Sr1_Energies.csv``
        * ``Sr1_Matrix_Elements.csv``
        * ``Sr1_Transition_Rates.csv``

6. Run calc_lifetimes.py

    * Input:
    
        * ``Sr1_Transition_Rates.csv``

    * Output:
    
        * ``Sr1_Lifetimes_Error_Check.csv``
        * ``Sr1_Transition_Rates_Error_Check.csv``
