pCI documentation |release|
===========================

pCI is a high-precision relativistic atomic structure software package written in Fortran and designed for use on modern HPC systems. 
The calculations are based on the configuration interaction method, and extendable to utilize second-order perturbation theory or the all-order method to attain even higher accuracy. 

pCI is developed and maintained at the University of Delaware. All of its development is done in the GitHub repository under the development branch.

We welcome all contributions to the pCI software package, including bug fixes, feature enhancements, and docuemtnation improvements. In order to contribute, open an `issue <https://github.com/ud-pci/pCI/issues/new/choose>`_ or a `pull request <https://github.com/ud-pci/pCI/pulls>`_ on GitHub. 

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

   examples/pci-py-sr.rst
   examples/neutral.rst
   examples/hci.rst
   

.. toctree::
   :maxdepth: 1
   :caption: About

   about
   publications
   

pCI Resources
-------------
pCI development is hosted at GitHub:

* https://github.com/ud-pci/pCI