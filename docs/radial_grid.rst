Radial grid
===========

Description
-----------

One-electron orbitals in pCI are represented on a one-dimensional radial grid. Inside the nucleus, the orbitals are described as a Taylor expansion in :math:`r/R_\text{nuc}`, where :math:`R_\text{nuc}` is the nuclear radius. Outside the nucleus, the orbitals are defined on a discrete set of radial grid points :math:`r_i`.

The radial grid uses a linear-logarithmic mapping from grid index :math:`i` to radial coordinate :math:`r_i`:

.. math::

    (i-1) H = \alpha (r_i - R_1) + \beta \ln\!\left(\frac{r_i}{R_1}\right),

where :math:`H` is a fixed step size, :math:`R_1` is the first grid point, and :math:`\alpha`, :math:`\beta` are weights controlling the balance between the linear and logarithmic terms. The logarithmic term produces a dense grid near the nucleus where orbitals vary rapidly, while the linear term extends the grid to large :math:`r` without requiring too many points.

The cavity radius :math:`R_2` (set via ``R2`` in ``HFD.INP``) fixes :math:`\alpha` by requiring the last grid point (:math:`i = N`) to land exactly at :math:`R_2`:

.. math::

    \alpha = \frac{(N-1)H - \beta\ln(R_2/R_1)}{R_2 - R_1}.

A larger cavity radius requires more grid points to maintain the same grid density, and the maximum number of grid points :math:`N` is set at compile time via the ``IP6`` parameter in ``lib/include/hfd.par``.

Changing the radial grid
------------------------

The default grid has **500 points**, suitable for cavity radii up to approximately 40 a.u. For larger cavities (e.g., 70 a.u. for neutral atoms), the grid must be expanded to accommodate more points. This requires recompiling several programs with updated array dimension parameters.

Because the grid size is encoded as compile-time Fortran parameters (``PARAMETER`` statements and ``.par`` include files), you must rebuild all affected programs after making these changes.

The following shows the changes needed for a **1500-point grid**:

1. ``bdhf``

    * Set ``NGP=1500`` in:

        * ``lib/all-order/basis/bdhf/bdhf.f``
        * ``lib/all-order/basis/bdhf/wfun.f``
        * ``lib/all-order/basis/bdhf/wint.f``
        * ``lib/all-order/basis/bdhf/yfun1.f``
        * ``lib/all-order/basis/bdhf/zfun.f``
        * ``lib/all-order/basis/bdhf/zint.f``

2. ``bspl``

    * Set ``NHF=1500`` in:

        * ``lib/all-order/basis/bspl/bspl.f``

    * Set ``NGP=1500`` in:

        * ``lib/all-order/basis/bspl/wfun.f``
        * ``lib/all-order/basis/bspl/wint.f``
        * ``lib/all-order/basis/bspl/yfun1.f``
        * ``lib/all-order/basis/bspl/yint1.f``

3. ``bas_wj``, ``bas_x``, ``hfd``, ``bass``, ``rpa``, ``rpa_dtm``

    These programs share dimension parameters via ``.par`` include files:

    * Set ``IPx6=1500`` in ``lib/all-order/include/basx.par``
    * Set ``IP6=1470`` in ``lib/include/hfd.par``

4. All-order and MBPT codes (``allcore-ci``, ``valsd-ci``, ``sdvw-ci``, ``second-ci``)

    * Set ``NHF=1500`` and ``NGP=1500`` in ``lib/all-order/include/global.par``

5. ``spl.in`` and ``bas_wj.in``

    * In the last line of each file (format: ``hp h0 nmax``), set ``nmax=1500`` to match the new grid size

After editing these files, rebuild the affected programs:

.. code-block::

    $ cd lib
    $ make install

.. note::

    The parameters ``NGP``, ``NHF``, ``IPx6``, and ``IP6`` define the maximum array sizes allocated at compile time. They must be at least as large as the number of grid points :math:`N` required for the chosen cavity radius :math:`R_2`. Setting them larger than needed wastes memory but does not cause errors; setting them too small will cause the program to abort with an error message.
