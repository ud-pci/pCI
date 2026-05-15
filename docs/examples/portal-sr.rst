Generating portal data for Sr
=============================

In this example, we will utilize the pCI-py scripts detailed :doc:`here <../pci-py>` to generate atomic data for neutral Sr. The configuration is based on the UD Darwin cluster.

config.yml
----------
The ``config.yml`` file defines the entire calculation from beginning to end. You can find a detailed explanation of the contents of this configuration file :doc:`here <../pci-py>`. 

.. raw:: html

    <details>
    <summary>
        Click here to see the full contents of the sample <a href="../src/config_Sr.yml" download>config.yml</a> we will be using.
    </summary>

.. code-block:: 

    # System parameters
    system:
        bin_directory: ""
        generate_directories: False
        run_codes: True
        on_hpc: True
        pci_version: default

    # HPC parameters
    hpc:
        partition: large-mem
        nodes: 1
        tasks_per_node: 64
        submit_job = False

    # Atomic information
    atom:
        name: Sr
        isotope: 
        include_breit: True
        code_method: [ci+all-order, ci+second-order]

    # Parameters used by basis programs
    basis:
        cavity_radius: 70
        diagonalized: True
        orbitals:
            core: 1s 2s 2p 3s 3p 3d 4s 4p 
            valence: 5s 5p 4d 6s 6p 5d 7s 7p 6d
            nmax: 35
            lmax: 5
        b_splines:
            nmax: 40
            lmax: 6
            k: 7
        val_aov:
            s: 5
            p: 5
            d: 5
            f: 3
        val_energies:
            kval: 1
            energies: 
                s: -0.28000
                p: [-0.22000, -0.22000]
                d: [-0.31000, -0.31000]
                f: [-0.13000, -0.13000]

    # Parameters used by add program
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

    # Parameters used by conf program
    conf:
        odd:
            J: 0.0
            JM: 0.0
            J_selection: False
            num_energy_levels: 24
            num_dvdsn_iterations: 50
        even:
            J: 0.0
            JM: 0.0
            J_selection: False
            num_energy_levels: 24
            num_dvdsn_iterations: 50
        include_lsj: True
        write_hij: False

    # Parameters used by dtm program
    dtm:
        include_rpa: True
        DM: 
            matrix_elements: 
            level_range: 
                odd: 
                even: 
        TM:
            matrix_elements: E1
            from:
                parity: odd
                level_range: 1 24
            to:
                parity: even 
                level_range: 1 24

    # Parameters used by portal script
    portal:
        ignore_g: True

    # Optional parameters
    optional:
        qed:
            include: False

        isotope_shifts: 
            include: False
            K_is: 0
            C_is: 0

.. raw:: html

    </details>

|

General parameters
------------------

.. code-block:: 

    # System parameters
    system:
        bin_directory: ""
        generate_directories: False
        run_codes: True
        on_hpc: True
        pci_version: default

    # HPC parameters
    hpc:
        partition: large-mem
        nodes: 1
        tasks_per_node: 64
        submit_job = False

    # Atomic information
    atom:
        name: Sr
        isotope: 
        include_breit: True
        code_method: [ci+all-order, ci+second-order]

Parameters used by basis programs
---------------------------------

.. code-block:: 

    # Parameters used by basis programs
    basis:
        cavity_radius: 70
        diagonalized: True
        orbitals:
            core: 1s 2s 2p 3s 3p 3d 4s 4p 
            valence: 5s 5p 4d 6s 6p 5d 7s 7p 6d
            nmax: 35
            lmax: 5
        b_splines:
            nmax: 40
            lmax: 6
            k: 7
        val_aov:
            s: 5
            p: 5
            d: 5
            f: 3
        val_energies:
            kval: 1
            energies: 
                s: -0.28000
                p: [-0.22000, -0.22000]
                d: [-0.31000, -0.31000]
                f: [-0.13000, -0.13000]


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
        odd:
            J: 0.0
            JM: 0.0
            J_selection: False
            num_energy_levels: 24
            num_dvdsn_iterations: 50
        even:
            J: 0.0
            JM: 0.0
            J_selection: False
            num_energy_levels: 24
            num_dvdsn_iterations: 50
        include_lsj: True
        write_hij: False

Parameters used by dtm program
------------------------------

.. code-block:: 

    dtm:
        include_rpa: True
        DM: 
            matrix_elements: 
            level_range: 
                odd: 
                even: 
        TM:
            matrix_elements: E1
            from:
                parity: odd
                level_range: 1 24
            to:
                parity: even 
                level_range: 1 24

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

Running pCI-py scripts
----------------------
1. Write config.yml to root directory.
   
    * Make sure to re-define ``optional.pci_version = default`` to the correct version if parameters such as the radial grid size has to be changed. On UD computers, one can use the command ``vpkg_versions pci`` to list all version of pCI.
    * Make sure to set ``optional.run_codes = True`` in order for job scripts to be automatically submitted to SLURM.

2. Run basis.py 
    
    * This script will generate basis sets for CI+all-order and CI+MBPT in their respective directories.
  
    * Input:

        * ``config.yml`` (specifically blocks ``system``, ``hpc``, ``atom``, ``basis``, ``optional``)

    * Output:

        * ``/CI+all-order/basis/``
        * ``/CI+second-order/basis/``

3. Run ci.py
   
    * Make sure to set ``optional.run_codes = True`` in order for job scripts to be automatically submitted to SLURM for CI programs to run.
  
    * Input:

        * ``config.yml`` (specifically blocks ``system``, ``hpc``, ``atom``, ``basis``, ``conf``, ``optional``)
        * ``/CI+all-order/basis/BASS.INP``
        * ``/CI+second-order/basis/BASS.INP``
  
    * Output:

        * ``/CI+all-order/even0/CONFFINAL.RES``
        * ``/CI+all-order/odd0/CONFFINAL.RES``
        * ``/CI+second-order/even0/CONFFINAL.RES``
        * ``/CI+second-order/odd0/CONFFINAL.RES``

4. Run dtm.py
   
    * Make sure to set ``dtm.TM.matrix_elements = E1`` so ``E1.RES``, which contains a table of E1 transitions, is formed.
    * Make sure to set ``dtm.include_rpa = True`` to include RPA corrections.

    * Input:
  
        * ``config.yml`` (specifically blocks ``system``, ``hpc``, ``atom``, ``conf``, ``dtm``)
        * ``/CI+all-order/even0/``, ``/CI+all-order/odd0/``, ``/CI+second-order/even0/``, ``/CI+second-order/odd0/``
        
            * ``CONF.INP``
            * ``CONF.DAT``
            * ``CONF.DET``
            * ``CONF.XIJ``
            * ``CONFSTR.RES``

    * Output:
  
        * ``/CI+all-order/tm/``, ``/CI+second-order/tm/``
  
            * ``TM.RES``
            * ``E1.RES``
 
5. Run gen_portal_csv.py
   
    * Make sure to set ``portal.ignore_g = True`` if configurations with :math:`g` orbitals or terms with :math:`G` are to be ignored in final outputs.
    * If ``CI+second-order`` directory is not found, uncertainties will be set to 0.
    
    * Input:
  
        * ``/CI+all-order/even0/``, ``/CI+all-order/odd0/``, ``/CI+second-order/even0/``, ``/CI+second-order/odd0/``

            * ``CONFFINAL.RES``
  
        * ``/CI+all-order/tm/``, ``/CI+second-order/tm/`` (optional)
  
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
