Isotope shift calculations
==========================

To run isotope shift calculations using the pCI code package, one simply has to add a few keys to the ``HFD.INP``, ``BASS.INP`` and ``CONF.INP`` files, then run the CI procedure as usual. The main keys used for isotope shifts are the following:  

* ``K_is`` - defines type of isotope shift to calculate  
  
    * ``K_is = 0`` - no isotope shift
    * ``K_is = 1`` - field shift
    * ``K_is = 2`` - specific mass shift
    * ``K_is = 3`` - normal mass shift
    * ``K_is = 4`` - total mass shift  
  
* ``C_is`` - prefactor for IS perturbation (C_is = dR_nucl/R_nucl)
* ``Klow`` - lower component key 

.. note::

	For isotope shifts, it is important to recompile ``basc`` and ``conf`` with arrays for 2-electron integrals set to ``double precision`` type. This can be done in the ``params.f90`` file by changing the value of ``Integer, Parameter :: type2_real = sp`` from ``sp``  to ``dp``. You can confirm this change in the title of the ``BASC.RES`` and ``CONF.RES`` output files, which should mention ``double precision for 2e integrals``.

via pCI-py scripts
------------------

Isotope shift calculations can also be automated using the pCI-py scripts. To do this, simply set the keys in the ``config.yml`` file. For example, we can set our field shift calculation with ``C_is=0.01``:

.. code-block:: 

    optional:
        isotope_shifts: 
            include: True
            K_is: 1
            C_is: 0.01

