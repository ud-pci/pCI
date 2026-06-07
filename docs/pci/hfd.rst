hfd - hartree-fock-dirac
------------------------

The ``hfd`` program solves restricted Hartree-Fock-Dirac (HFD) equations self-consistently under the central field approximation to find four-component Dirac-Fock (DF) orbitals and eigenvalues of the HFD Hamiltonian. The program provides the initial approximation, storing both basis radial orbitals

.. math::

    \phi_{nlj}\equiv r\left(\begin{array}{c}f_{nlj}\\-g_{nlj}\end{array}\right),

as well as the radial derivatives of the orbitals :math:`\partial_r\phi_{nlj}`, to the file ``HFD.DAT``.
A more detailed description of this program is given in the 2015 CI-MBPT paper\ :footcite:p:`KozPorSaf15`.

Input Files:

* ``HFD.INP`` - atomic system definition, HFD parameters, and orbital list
* ``HFD.DAT`` (optional) - pre-existing orbitals; required when ``kl = 1`` or ``kl = 2``

Output Files:

* ``HFD.DAT`` - basis radial orbitals :math:`\phi_{nlj}` and radial derivatives :math:`\partial_r\phi_{nlj}`
* ``HFD.RES`` - program log

The ``hfd`` program reads its parameters from ``HFD.INP``. The first line is a free-form label for the atom or ion. This is followed by a block of keyword parameters, which can appear in any order and are case-insensitive:

.. code-block::

     KL =   0     # calculation mode
     NS =  18     # total number of orbitals
     NSO=   2     # number of closed core shells
     Z  =  26.0   # atomic number
     AM =  56.000 # atomic weight
     JM =  -2.0   # angular momentum projection (-2 to disable)
     R2 =  20.0   # radial grid outer boundary (a.u.)
     kbr=   2     # Breit interaction (0 - Coulomb, 1 - Gaunt, 2 - Full Breit)
    rnuc=   0.0   # nuclear rms radius (fm); 0 = derived from AM
    k_is=   0     # isotope shift (0 - off, 1 - VS, 2 - SMS, 3 - NMS, 4 - MS)
    c_is=   0.0   # isotope shift scaling factor
    n_is=   0     # number of core shells for SMS potential (default: NSO)

The value of ``KL`` can take the following values:

* ``KL = 0`` - start a new self-consistent HFD calculation from scratch
* ``KL = 1`` - start a new calculation, reading frozen orbitals (``KP = 1``) from an existing ``HFD.DAT``
* ``KL = 2`` - read all orbitals from an existing ``HFD.DAT`` without re-solving (e.g., to add Breit corrections)

The value of ``k_is`` selects the type of isotope shift perturbation added to the HFD potential:

* ``k_is = 0`` - no isotope shift
* ``k_is = 1`` - volume shift (VS)
* ``k_is = 2`` - specific mass shift (SMS)
* ``k_is = 3`` - normal mass shift (NMS)
* ``k_is = 4`` - total mass shift (SMS + NMS)

``c_is`` is a scaling factor applied to the isotope shift perturbation. ``n_is`` specifies the number of core shells included in the SMS core potential; it defaults to ``NSO`` if not set. ``rnuc`` overrides the nuclear rms radius (in fm) used in the nuclear potential; it can be omitted entirely or set to zero, in which case it is derived automatically from the atomic weight ``AM``.

After the keyword block, three header lines are followed by the orbital list. Each line specifies one orbital in fixed format:

.. code-block::

         NL   J       QQ     KP   NC

  1     1S (1/2)   2.0000    1    0
  2     2S (1/2)   2.0000    1    0
  ...
 10     4S (1/2)   1.0000    0    1

The columns are:

* **index** - orbital index (informational only)
* **NL** - principal quantum number and angular momentum letter (e.g., ``4S``, ``2P``, ``3D``)
* **J** - total angular momentum as a fraction (e.g., ``(1/2)``, ``(3/2)``)
* **QQ** - occupation number
* **KP** - ``0`` to solve the HFD equation for this orbital; ``1`` to read it from an existing ``HFD.DAT`` (frozen)
* **NC** - configuration group index (``0`` for core; ``1``, ``2``, ... for valence groups solved sequentially)

Running hfd
~~~~~~~~~~~

To run ``hfd``, run the command:

.. code-block::

    hfd

**References**

.. footbibliography::
