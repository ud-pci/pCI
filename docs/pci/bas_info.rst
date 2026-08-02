bas_info - basis set diagnostics
---------------------------------

The ``bas_info`` program is an interactive diagnostic utility for inspecting a basis set stored in ``HFD.DAT`` or ``CONF.DAT``. It prints a table of orbital energies, RMS radii, and node counts, and checks partial-wave completeness of the basis. Additional numerical tests can then be run interactively.

Input Files:

* ``HFD.DAT`` or ``CONF.DAT`` - basis set from ``hfd`` or ``pbasc``

Output Files:

* ``BAS_INFO.RES`` - program output

Running bas_info
~~~~~~~~~~~~~~~~

To run ``bas_info``, run the command:

.. code-block::

    bas_info

The program first prompts for the DAT filename:

.. code-block::

    Give file name or ENTER to quit:

Enter the filename (e.g. ``HFD.DAT``), ``1`` as a shortcut for ``HFD.DAT``, or ``2`` for ``CONF.DAT``. Pressing Enter without input exits the program.

After reading the file, ``bas_info`` prints the orbital table and partial-wave completeness check automatically, then presents the optional-test menu:

.. code-block::

    Choose optional tests (0 to quit)
    1. closure r^2=r|n><n|r is checked for each
       orbital and for three sets of
       intermediate states |n> which
       meets dipole selection rules.
    2. identity int_0^infty{1/2*f*df/dr}=0 is checked
       for each orbital
    3. Checks orthogonality
    4. Checks Taylor expansion at the origin
    5. prints any orbital of your choice

Enter a test number to run it; enter ``0`` to quit. The menu repeats until ``0`` is entered.

Default output
~~~~~~~~~~~~~~

**Orbital table** — printed for every orbital in the basis:

* orbital index, quantum numbers :math:`n`, :math:`l`, :math:`j`
* energy eigenvalue :math:`-P(ii+1)`
* large component :math:`P` and small component :math:`Q` at the last grid point (a check that both decay to zero)
* RMS radius :math:`R_\text{rms} = \sqrt{\langle r^2 \rangle}`
* number of nodes in :math:`P`

A warning is printed and the program pauses if the norm :math:`\int (P^2 + Q^2)\,dr` deviates from 1 by more than :math:`10^{-3}`.

**Partial-wave completeness** — for each distinct :math:`(l, 2j)` channel present in the basis, the program accumulates :math:`\sum_n \langle n \rangle P_n(r)` and reports its average value (as a percentage of 1) over five equal segments of the radial grid. Values close to 100% indicate that the partial wave is well-represented over that segment of the radial grid.

Optional tests
~~~~~~~~~~~~~~

**Test 1 — Closure (dipole sum rule)**

Checks the sum rule

.. math::

    \langle a | r^2 | a \rangle = \sum_n |\langle a | r | n \rangle|^2

for each orbital :math:`a`, summed over all intermediate orbitals :math:`n` that satisfy the dipole selection rules (:math:`\Delta l = 0, \pm 1`; :math:`\Delta j = 0, \pm 1`). The ratio RHS/LHS is printed for each intermediate :math:`(l, 2j)` channel; a value of 1 indicates the channel is complete for that selection rule.

**Test 2 — Derivative identity**

Checks the identity

.. math::

    \int_0^\infty \left(P\frac{dP}{dr} + Q\frac{dQ}{dr}\right) dr = 0

for each orbital, evaluating it both with the numerical derivative and with the derivative stored in the DAT file. Both results are normalised by :math:`\max(P^2 + Q^2)` and printed side by side.

**Test 3 — Orthonormality**

Computes all overlap integrals :math:`\langle n_i | n_k \rangle` between pairs of orbitals with the same :math:`(l, 2j)`. Diagonal elements should be 1 (normalization) and off-diagonal elements should be 0 (orthogonality).

**Test 4 — Taylor expansion at origin**

Prompts for orbitals (1) or derivatives (2). For each orbital, compares the numerical grid values against the Taylor series

.. math::

    P_\text{Taylor}(r) = r^\gamma \sum_{m=0}^{9} c_m \left(\frac{r}{R_1}\right)^m

stored in the DAT file. The ratio (numerical)/(Taylor) is printed for several points near the origin, including a few extrapolated points before the first grid point. A ``>`` marker flags the first on-grid point.

**Test 5 — Print orbital**

Prompts for an orbital index and writes a table of :math:`R(i)`, :math:`P(i)`, :math:`Q(i)` for all radial grid points to ``BAS_INFO.RES``. Enter ``0`` to return to the main menu.
