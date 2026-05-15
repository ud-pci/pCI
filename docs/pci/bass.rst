bass - constructing the basis set
---------------------------------

The ``bass`` program constructs the basis set starting from the HFD orbitals for the core and valence orbitals, formed by the ``hfd`` program. Then virtual orbitals are added to account for correlations and are constructed from either (1) previously constructed HFD orbitals or (2) B-splines. A reasonable basis set should consist of orbitals mainly localized at the same distances from the origin as the valence orbitals. 

In the first case, virtual orbitals are formed using a recurrent procedure\ :footcite:p:`BDT77`.  The lowest virtual orbitals can be constructed from the HFD orbitals. 
The large component of the radial Dirac bispinor, :math:`f_{n'l'j'}`, is obtained from a function :math:`f_{nlj}` constructed previously by multiplying it by :math:`r^{l' - l}\, \sin(kr)`. Here :math:`l'` and :math:`l` are the orbital quantum numbers of the new and old orbitals (:math:`l' \geq l`) and the coefficient :math:`k` is determined by the properties of the radial grid. The small component :math:`g_{n'l'j'}` is found from the kinetic balance condition:

.. math::
    
    g_{n'l'j'} =\frac{\boldsymbol{\sigma} \bf p}{2mc} f_{n'l'j'} ,

where :math:`\boldsymbol{\sigma}` are the Pauli matrices, :math:`\bf p` and :math:`m` are the electron momentum and mass, and :math:`c` is the speed of light.
The newly constructed functions are then orthonormalized to the functions of the same symmetry.

Another option is to construct large components of the orbitals from B-splines. Small components are still formed with the kinetic balance method. A more detailed description of this program is given in the 2015 CI-MBPT paper\ :footcite:p:`KozPorSaf15`. 

**References**

.. footbibliography:: 