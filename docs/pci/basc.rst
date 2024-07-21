basc - calculating radial integrals
-----------------------------------

After the configuration list has been created, the next step is to calculate the radial integrals using the program ``basc``. ``basc`` calculates one-electron and two-electron radial integrals, which are used by the ``conf`` program to form the Hamiltonian in the CI space. The one-electron radial integrals correspond to the DF potential of the core, and the two-electron radial integrals account for the Coulomb and Breit interactions between the valence electrons. The matrix elements of the Coulomb interaction for the multipolarity :math:`k` can be written as

.. math::
    
    \langle c,d|V_q^k|a,b\rangle \equiv G_q^k(ca) G_q^k(bd) R_{abcd}^k,


where the angular factors :math:`G_q^k(fi)` (known as relativistic Gaunt coefficients) are given by

.. math::

    G_q^k(fi)=(-1)^{m_f+1/2}\delta_p\sqrt{(2j_i+1)(2j_f+1)}
        \begin{pmatrix} 
         j_f & j_i & k \\  
        -m_f & m_i & q
        \end{pmatrix}
        \begin{pmatrix} 
        j_f & j_i & k \\  
        1/2 & -1/2 & 0
        \end{pmatrix},
  
and :math:`R_{abcd}^k` are the relativistic Coulomb radial integrals, and :math:`\delta_p` accounts for the parity selection rule

.. math:: 

    \delta_p=\xi(l_i+l_f+k), \hspace{0.2in}\xi(n)=\Bigg\{
    \begin{matrix}
    1 & \text{if \( n \) is even,} \\ 
    0 & \text{if \( n \) is odd.}
    \end{matrix} 

The Breit interaction has the same form as the Coulomb interaction, but without the parity selection rule. 

The ``basc`` reads in the files ``HFD.DAT`` and ``CONF.INP`` to determine which radial integrals are needed. These integrals are calculated and written to the files ``CONF.INT``. The relativistic Gaunt coefficients are written to the file ``CONF.GNT``, and the file ``CONF.DAT`` is also formed, storing the basis radial orbitals :math:`\phi_{nlj}`, as well as functions :math:`\chi_{nlj} = h_\text{DF}^r\phi_{nlj}`. 
