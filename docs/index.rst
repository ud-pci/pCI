pCI documentation |release|
===========================

pCI is a high-precision relativistic atomic structure software package written in Fortran and designed for use on modern HPC systems. 
The calculations are based on the configuration interaction method, and extendable to utilize second-order perturbation theory or the all-order method to attain even higher accuracy. 

pCI is developed and maintained at the University of Delaware. All of its development is done in the GitHub repository under the development branch.

We welcome all contributions to the pCI software package, including bug fixes, feature enhancements, and documentation improvements. In order to contribute, open an `issue <https://github.com/ud-pci/pCI/issues/new/choose>`_ or a `pull request <https://github.com/ud-pci/pCI/pulls>`_ on GitHub. 

To learn how to use pCI, there are walk-through guides and smaller examples that demonstrate how to use different parts of the pCI software package. 

.. toctree::
   :maxdepth: 1
   :caption: Using pCI

   installation
   ud_instructions


.. toctree::
   :maxdepth: 1
   :caption: pCI Technical Details

   theory
   pci/index.rst
   extensions/index.rst
   upscaling
   supplementary
   pci-py
   isotope_shifts
   qed
   portal_codes
   radial_grid


.. toctree:: 
   :maxdepth: 1
   :caption: Examples

   examples/portal-sr.rst
   examples/pci-py-sr.rst
   examples/neutral.rst
   examples/hci.rst
   

.. toctree::
   :maxdepth: 1
   :caption: About

   about
   publications
   citation
   

pCI Resources
-------------
pCI development is hosted at GitHub:

* https://github.com/ud-pci/pCI

Citation
--------
If you use pCI, please acknowledge it by citing the main publication listed here:

C\. Cheung, M. G. Kozlov, S. G. Porsev, M. S. Safronova, I. I. Tupitsyn, and A. I. Bondarev, pCI: A parallel configuration interaction software package for high-precision atomic structure calculations, Computer Physics Communications 308, 109463 (2025).

.. code-block:: bibtex

   @article{2025pCI,
     title = {pCI: A parallel configuration interaction software package for high-precision atomic structure calculations},
     journal = {Computer Physics Communications},
     volume = {308},
     pages = {109463},
     year = {2025},
     author = {Charles Cheung and Mikhail G. Kozlov and Sergey G. Porsev and Marianna S. Safronova and Ilya I. Tupitsyn and Andrey I. Bondarev},
   }
