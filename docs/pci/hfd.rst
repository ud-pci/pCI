hfd - hartree-fock-dirac
------------------------

The ``hfd`` program solves restricted Hartree-Fock-Dirac (HFD) equations self-consistently under the central field approximation to find four-component Dirac-Fock (DF) orbitals and eigenvalues of the HFD Hamiltonian. The program provides the initial approximation, storing both basis radial orbitals

.. math::
    
    \phi_{nlj}\equiv r\left(\begin{array}{c}f_{nlj}\\-g_{nlj}\end{array}\right),

as well as the radial derivatives of the orbitals :math:`\partial_r\phi_{nlj}`, to the file ``HFD.DAT``. 