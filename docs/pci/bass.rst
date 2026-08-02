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

Input Files:

* ``BASS.INP`` - basis set parameters and orbital list
* ``HFD.DAT`` - core and valence orbitals from ``hfd``

Output Files:

* ``HFD.DAT`` - updated with newly constructed virtual orbitals appended to the existing set
* ``BASS.RES`` - program log

The ``bass`` program reads its parameters from ``BASS.INP``. The first line is a free-form label for the atom or ion. This is followed by a block of keyword parameters:

.. code-block::

     Fe 16+
      Z = 26.0
     Am = 52.0
     Nso=    4  # number of core orbitals (defines DF operator)
     Nv =  129  # number of valence & virtual orbitals
     Ksg=    1  # defines Hamiltonian: 1 - DF, 3 - DF+Breit
     Kdg=    0  # diagonalization of Hamiltonian (0 - no, 1, 2 - yes)
     orb= 4s 1  # first orbital for diagonalization
     Kkin    1  # kinetic balance (0, 1, or 2)
     orb= 5s 1  # first orbital to apply kinetic balance
     orb= 2p 3  # last frozen orbital
     orb= 0p 3  # last orbital in basis set
    kout= 0     # detail level in the output
    kbrt= 2     # 0, 1, 2 - Coulomb, Gaunt, Breit
    k_is= 0     # isotope shift (0 - off, 1 - VS, 2 - SMS, 3 - NMS, 4 - MS)
    c_is= 0.0   # isotope shift scaling factor

The value of ``Ksg`` can take the following values:

* ``Ksg = 1`` - Dirac-Fock Hamiltonian
* ``Ksg = 3`` - Dirac-Fock Hamiltonian with Breit correction

The value of ``Kdg`` can take the following values:

* ``Kdg = 0`` - no diagonalization
* ``Kdg = 1`` or ``Kdg = 2`` - diagonalize the Hamiltonian matrix

The ``orb`` keyword specifies orbitals using the format ``nl j``, where ``n`` is the principal quantum number, ``l`` is the angular momentum letter, and ``j`` is the numerator of :math:`j = j_\text{value}/2`. It is used four times:

1. First orbital for Hamiltonian diagonalization
2. First orbital to apply kinetic balance
3. Last frozen orbital (from ``HFD.DAT``, not recomputed)
4. Last orbital to retain in the final basis set

The value of ``Kkin`` can take the following values:

* ``Kkin = 0`` - :math:`E = V = 0` kinetic balance
* ``Kkin = 1`` - :math:`E = 0,\, V_c` kinetic balance
* ``Kkin = 2`` - :math:`E_\mathrm{HF},\, V_c` kinetic balance

After the keyword block, a separator line (``---...``), the list of ``Nso`` core orbital quantum numbers is given in fixed format using the pCI floating-point notation, followed by the list of virtual orbitals to construct. Each virtual orbital line is also in fixed format, specifying the pCI orbital identifier and optionally a basis type ``Kbas`` and starting orbital:

.. code-block::

    ----------------------------------------------------------
         0.1002     0.2002    -0.2102     0.2104

      1  0.3001
      2 -0.3101
      3  0.3101
      6  0.4001  3  0.4001

In the pCI orbital notation, the absolute value encodes ``n`` and ``l`` (e.g., ``0.4001`` → :math:`n=4`, :math:`l=0`), and the sign distinguishes :math:`j = l + 1/2` (positive) from :math:`j = l - 1/2` (negative). Blank lines between orbitals are permitted.

The optional ``Kbas`` field on a virtual orbital line specifies the construction method. When ``Kbas`` and the following starting orbital ``Qnl1`` are omitted, the orbital is built from the lowest available same-symmetry orbital in ``HFD.DAT`` by the standard recurrent procedure.

Running bass
~~~~~~~~~~~~

To run ``bass``, run the command:

.. code-block::

    bass

**References**

.. footbibliography::
