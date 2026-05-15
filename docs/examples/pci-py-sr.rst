Atomic properties of neutral Sr
===============================

In this example, we will utilize the pCI-py scripts detailed :doc:`here <../pci-py>` to calculate energies, :math:`E1` reduced matrix elements, and dc and ac polarizabilities for the :math:`^1S_0` state of neutral Sr. 

We treat Sr as a divalent ion, constructing the basis set in the :math:`V^{N-2}` approximation, where :math:`N` is the number of electrons. The ``basis.py`` script will construct a basis set in both CI+all-order and CI+MBPT approximations, in their respective directories. Both are constructed as non-diagonalized basis sets in a spherical cavity radius of 70 a.u. The initial HFD self-consistency procedure is carried out for the core :math:`1s`, :math:`2s`, :math:`2p`, :math:`3s`, :math:`3p`, :math:`3d`, :math:`4s`, and :math:`4p` orbitals, and then the valence :math:`5s`, :math:`5p`, :math:`4d`, :math:`6s`, :math:`6p`, :math:`5d`, :math:`7s`, :math:`7p`, and :math:`6d` orbitals are constructed in the frozen-core potential. The remaining virtual orbitals are formed using 40 B-spline orbitals of order 7, including partial waves up to :math:`l=5`. The coupled-cluster equations are then solved to all-order and to second-order in the respective directories. These basis set is stored in the file ``HFD.DAT``, while corrections from all-order and MBPT are stored in the form of effective radial integrals in the files ``SGC.CON`` and ``SCRC.CON``. 

Next, the ``ci.py`` script generates lists of configurations by exciting electrons from the odd-parity :math:`5s 5p` and even-parity :math:`5s^2` configurations to all orbitals up to :math:`17spdfg`, then run the CI computations. The energies and wave functions of the lowest 12 states with :math:`J=1` are calculated for the odd-parity CI run, and the lowest 24 with :math:`J=0` for the even-parity CI run. 

The ``dtm.py`` script is then run to calculate :math:`E1` reduced matrix elements, including RPA corrections. In ``config.yml`` we specify transitions from the first 3 odd states to the first even state. This will compute reduced matrix elements for the :math:`^3P_1^o-\,^1S_0` and :math:`^1P_1^o-\,^1S_0` transitions. Running ``dtm``, we obtain :math:`|\langle ^3P_1^o || E1 || ^1S_0\rangle|=0.155` a.u. and :math:`|\langle ^1P_1^o || E1 || ^1S_0\rangle|=5.275` a.u.

Finally, the ``pol.py`` script is run to calculate dc and ac valence polarizabilities for the :math:`^1S_0` state. Doing so, we obtain dc polarizability :math:`\alpha(^1S_0)=192.78` a.u., and ac polarizability :math:`\alpha(^1S_0)=242.97` a.u. at :math:`\lambda=1000` nm. To calculate polarizabilities for the $^3P_0^o$ state, one has to set ``pol.parity`` to ``odd``, then swap the values of ``conf.odd.J`` and ``conf.odd.JM`` with ``conf.even.J`` and ``conf.even.JM``. The ``ci.py`` script is re-run to calculate wave functions in the respective projections, and then ``pol.py`` is re-run to obtain polarizabilities. 

As an extra test, one can calculate valence polarizabilities for the :math:`^3P_1^o` state by setting ``conf.odd.J`` and ``conf.odd.JM`` to ``1``, ``conf.even.J`` and ``conf.even.JM`` to ``0``, and re-running ``pol.py``. Doing so, we obtain scalar polarizability :math:`\alpha_0(^3P_1)=194.84` a.u., vector polarizability :math:`\alpha_2(^3P_1)=23.45` a.u., and total ac polarizability :math:`\alpha(^3P_1)=218.29` a.u. at :math:`\lambda=1000` nm. 

Higher precision of the polarizabilities can be attained by including core polarizabilities, which are calculated with a different program not included in this work.

config.yml
----------
The ``config.yml`` file defines the entire calculation from beginning to end using the pCI-py scripts. A detailed explanation of the contents of this configuration file can be found :doc:`here <../pci-py>`. 

.. collapse:: Click here to see config.yml

    The sample configuration file we used for this calculation can be downloaded :download:`here <../src/config_Sr.yml>`.

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
            basis_set: 17spdfg
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
                J: 1.0
                JM: 1.0
                J_selection: False
                num_energy_levels: 12
                num_dvdsn_iterations: 50
            even:
                J: 0.0
                JM: 0.0
                J_selection: False
                num_energy_levels: 24
                num_dvdsn_iterations: 50
            include_lsj: True
            write_hij: True

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
                    level_range: 1 3
                to:
                    parity: even 
                    level_range: 1 1

        # Parameters used by pol program
        pol:
            parity: even
            level: 1
            method: 1
            field_type: static, dynamic
            wavelength_range: 1000 1000
            step_size: 0

        # Optional parameters
        optional:
            qed:
                include: False

            isotope_shifts: 
                include: False
                K_is: 0
                C_is: 0

|