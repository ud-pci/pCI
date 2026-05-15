Changing the radial grid
========================

One feature(?) of the pCI code package stems from the original CI-MBPT code package. One-electron orbitals outside the nucleus are defined on a radial grid. Inside the nucleus, they are described in a form of the Taylor expansion over :math:`r/R`, where :math:`R` is the nuclear radius. In the codes used to construct the basis set, the parameter ``R2`` sets the size of the radial grid in a.u.  

The default version of the codes come with a grid of ``xxx`` nodes. In order to rebuild the code with a different size grid, several parameters must be changed in the basis set codes. For example, the following parameters sets the codes for a 1500-pt grid:

1. ``bdhf``
   
    * set ``NGP=1500`` in:
    
        * ``bdhf.f``
        * ``wfun.f``
        * ``wint.f``
        * ``yfun1.f``
        * ``zfun.f``
        * ``zint.f``

2. ``bspl`` 

    * set ``NHF=1500`` in:

        * ``bspl.f``

    * set ``NGP=1500`` in:

        * ``wfun.f``
        * ``wint.f``
        * ``yfun1.f``
        * ``yint1.f``
  
3. ``bas_wj``, ``bas_x``, ``hfd``, ``bass``, ``rpa``, ``rpa_dtm``

    * set ``IPx6=1500`` in:

        * ``bas_x.par``

    * set ``IP6=1470`` in:

        * ``hfd.par``

4. ``allcore-rle-ci``, ``valsd-rle-cis``, ``sdvw-rle-cis``, ``second-cis``

    * set ``NHF=1500`` in:

        * ``all.par``
        * ``second-cis.par``
        * ``archiv.par``

    * set ``NGP=1500`` in:

        * ``archiv.par``

5. ``spl.in`` and ``bas_wj.in``

    * set number of points to 0.1 and 1500