Basis for CI+all-order and CI+MBPT
==================================

*The following instructions assume familiarity with the main programs of the pCI package.*

CI+all-order basis sets
-----------------------

In this section, we describe the general method of building basis sets for the CI+all-order and CI+MBPT code packages. As the CI and all-order code packages were developed separately, they use basis sets in different formats. CI uses ``HFD.DAT`` and all-order uses ``hfspl.1`` and ``hfspl.2`` files. For all-order/CI+all-order/CI+MBPT calculations, most of the basis is constructed via the B-splines. However, one needs too many B-splines to reproduce core and lower valence stats with high enough accuracy for heavier atoms. Therefore, core and a few few valence electron wave functions are taken from Dirac-Hartree-Fock (DHF), and a combined basis with splines is built. More splines will mean larger basis for the CI as well, so we would like to avoid this. There also seems to be an issue with the basis converter codes for larger numbers of B-splines.

All-order or MBPT calculations involve sums over all possible states. To computer these accurately, one needs a large basis. Generally, we use ``lmax=6`` and ``Nmax=35`` for each of the partial waves :math:`(1-35)s, (2-35)p_{1/2}, (2-35)p_{3/2}, \dots`. CI does not require such large basis sets and it is reducing for CI computations. 

.. note::

    More splines will mean larger basis for the CI as well, so we would like to avoid this. There also seems to be an issue with the basis converter codes for larger numbers of B-splines.

.. note::

    There is a technical issue of CI codes using Taylor expansion inside the nucleus, while all-order codes use a radial grid that starts from the origin. This makes format conversion somewhat imperfect near the nucleus.

Instructions
------------

Now we will discuss the steps to build the basis set for CI+all-order and CI+MBPT calculations.

1. Produce B-splines
    
	* ``./tdhf < bas_wj.in`` - solves DHF equations (reads ``bas_wj.in`` and writes ``fort.1``)
	* ``./nspl40 < spl.in`` - produces B-spline basis (reads ``fort.1`` and writes ``hfspl.1`` and ``hfspl.2``)

2. Run HFD code
   
	* ``./hfd`` - solves DHF equations (reads ``HFD.INP`` and writes ``HFD.DAT``)  

3. Convert B-spline basis to HFD.DAT format
	
    * ``./bas_wj`` - converts ``hfspl.1`` and ``hfspl.2`` B-spline files and writes ``WJ.DAT`` and ``BASS.INP``  

	.. note::
	    
	    There are a couple of changes that have to be made to the resulting ``BASS.INP`` file:

	    * The line ``lst= 0s 1# last orbital to be kept in the basis set`` has to be deleted.  
	    * The line ``orb= 0s 1# first orbital to apply kin.bal.`` has to be changed to reflect the first orbital not from ``hfd``.  
	    * The line ``orb= 5s 1# first orbital for diagonalization`` can either be kept or changed to reflect the first orbital not from ``hfd``.  
	    * The lines ``1  0.5001  3  0.5001`` that include orbitals from ``hfd`` have to be changed to ``1  0.5001``

4. Build combined basis from ``HFD.DAT`` and ``WJ.DAT``
	
    * ``./bass`` - reads ``BASS.INP``, ``HFD.DAT`` and ``WJ.DAT`` and writes new ``HFD.DAT``  

5. Convert combined HFD.DAT to format used for all-order and second-order codes
	
    * ``rm hfspl.1 hfspl.2`` - erase the files ``hfspl.1`` and ``hfspl.2``
  
	* ``./bas_x`` - reads ``HFD.DAT`` and writes ``hfspl.1`` and ``hfspl.2``  
  
.. note::

    ``bdhf`` and ``bspl40`` are used in place of ``tdhf`` and ``nspl40`` if Breit corrections are included.

.. note::
    
	``hfd`` can work with partially opened shells, whereas ``tdhf`` and ``bdhf`` cannot.


The final ``hfspl.1`` and ``hfspl.2`` files are the basis set files that can be read from the :doc:`all-order part of the package <all-order>`.