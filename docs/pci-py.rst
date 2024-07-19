pci-py scripts
==============

In this section, we describe supplementary python scripts used to automate the entire pCI process, from the generation of basis sets and configuration lists to the actual calculations themselves. These python scripts read a user-defined YAML-formatted configuration file ``config.yml`` to set general parameters defining the atomic system of interest, as well as parameters defining the types of the calculations.

_These python scripts are currently designed primarily for massive data generation via CI+all-order and CI+MBPT calculations on the UD clusters. They can be modified in the future for more general calculations._

config.yml
----------

The ``config.yml`` file is a YAML-formatted configuration file that defines important parameters of the atomic system of interest. This config file is divided into several sections:

* general parameters
* parameters used by basis programs (only read by basis.py)
* parameters used by add program (only read by add.py)
* parameters used by conf program (only read by add.py)
* parameters used by dtm program (only read by dtm.py)
* optional parameters

The general parameters are in the block ``system`` and is read by all python scripts:

.. code-block:: 

    system:
        name: Sr
        isotope: 
        include_breit: True

The ``system.name`` field takes the name of the atomic system of interest.  
The ``system.isotope`` field takes a specified isotope number. If left blank, the script will automatically find the atomic weight from nuclear radii table.  
The ``system.include_breit`` field takes in a boolean to set whether to include Breit corrections or not.

basis.py
--------

The ``basis.py`` script automatically generates the basis set given information about the atomic system of interest via a configuration file ``config.yml``. This script reads the ``basis`` block:

.. code-block:: 

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

* The ``basis.cavity_radius`` field takes in a value to set the size of the spherical cavity.  
* The ``basis.orbitals`` block requires specification of the core and valence orbitals to be included in the basis. It also requires a maximum principal quantum number ``nmax`` and maximum partial wave ``lmax`` for basis set orbital generation.
* The ``basis.b_splines`` block requires specification of the maximum principal quantum number ``nmax`` (number of splines), maximum partial wave  ``lmax``, and order of splines.  
* The ``basis.val_aov`` block requires the number of valence orbitals to include for each partial wave.  
* The ``basis.val_energies`` block allows specification of the energies of the valence orbitals. 

In addition, the script will read the ``optional`` block:

.. code-block:: 

    optional:
        qed:
            include: False
            rotate_basis: False

        isotope_shifts: 
            include: False
            K_is: 1
            C_is: 0.01

        code_exec: ci+all-order
        run_ao_codes: False

* The ``optional.qed`` block allows users to specify inclusion of QED corrections (this is currently not supported).  
* The ``optional.isotope_shifts`` block allows users to specify isotope shift calculations by switching the ``include`` value to ``True`` and specifying keys ``K_is`` and ``C_is``.  
* The ``optional.code_exec`` field allows users to specify automatic execution of basis set codes depending on value.  
* The ``optional.run_ao_codes`` field allows users to specify whether they would like to run the all-order set of codes after construction of the basis set.

add.py
------

The ``add.py`` script automatically generates the list of configurations given information about the atomic system of interest via a configuration file ``config.yml``. Note that if ``BASS.INP`` is not present, the order of the conigurations will be by $nl$, and not how it is defined in the basis set construction. This script reads the ``add`` block:

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

* The ``add.ref_configs`` block requires a list of reference configurations to excite electrons from to construct the list of configurations for the CI calculation. This list will not be constructed if left blank for a specified parity.  
* The ``add.basis_set`` block requires specification of the basis set designated by ``nspdfg``, where ``n`` specifies the principal quantum number, and ``spdfg`` specifies inclusion of s, p, d, f, and g orbitals. One can include higher partial waves by appending to the end of the list ``h``, ``i``, ``k``, ...  
* The ``add.orbitals`` block allows full customization of allowed orbital occupancies. For example, ``1-2s: 2 2`` defines the 1s and 2s orbitals to be closed, ``2p: 6 6`` defines the 2p orbital to be closed, ``3-7p: 0 6`` defines 3p to 7p orbitals to be completely open to allow up to 6 electrons on those orbitals.  
* The ``add.excitations`` block defines the number of allowed excitations.

This script also reads parameters for the CI execution from the ``conf`` block:

.. code-block:: 

    conf:
        J: 0.0
        JM: 0.0
        J_selection: False
        num_energy_levels: 24
        num_dvdsn_iterations: 100
        include_lsj: True
        write_hij: False

* The ``conf.J`` field defines the total angular momentum of the energy levels.  
* The ``conf.JM`` field defines the projection of the total angular momentum.  
* The ``conf.J_selection`` field defines whether the user wants energy levels of a specific J value defined by J and JM.  
* The ``conf.num_energy_levels`` field defines the number of energy levels to be calculated.  
* The ``conf.num_dvdsn_iterations`` field defines the total number of Davidson iterations allowed.  
* The ``conf.include_lsj`` field defines whether the user wants expectation values \(L^2\) and \(S^2\) to be calculated.  
* The ``conf.write_hij`` field defines whether the user wants the Hamiltonian matrix to be written to file ``CONF.HIJ``.  

dtm.py
------

The ``dtm.py`` script automatically generates the density and/or transition matrix elements given output files from ``conf`` runs. This script reads the ``dtm`` block from ``config.yml``:

.. code-block:: 

    dtm:
        matrix_elements: E1
        include_rpa: True

* The ``dtm.matrix_elements`` field defines the types of matrix elements to be calculated. This can be a single matrix element or an array of matrix elements such as ``[E1, M2]``
* The ``dtm.include_rpa`` field defines whether the user would like to include RPA corrections.

gen_portal_csv.py
-----------------

The ``gen_portal_csv.py`` script generates csv-formatted datafiles of atomic energy levels and matrix elements given output files from ``conf`` and ``dtm`` runs. This script reads the ``portal`` block from ``config.yml``:

.. code-block:: 

    portal:
        ignore_g: True
        min_uncertainty: 1.5

* The ``portal.ignore_g`` field removes all atomic properties with :math:`>g` in the configuration or :math:`>G` terms.
* The ``portal.min_uncertainty`` field sets a minimum uncertainty in percentage for matrix value uncertainties. The script has predefined minimum uncertainties set for a few systems. 

calc_lifetimes.py
-----------------

The ``calc_lifetimes.py`` script generates csv-formatted datafiles of lifetimes and transition rates given output files from ``gen_portal_csv.py``.


_More information about atomic data generation for the UD ATOM portal can be found `here <portal_codes.rst>`_._
