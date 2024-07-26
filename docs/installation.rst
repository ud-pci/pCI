Installation
============

Required libraries
------------------
In order to compile pCI the following software libraries and tools are required: 

* `Intel Fortran compiler 2020u4+ <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html>`_
* `CMake v3.12+ <https://cmake.org/download>`_
* `OpenMPI v4.1+ <https://www-lb.open-mpi.org/software/ompi>`_
* `Python v3.x <https://www.python.org/downloads/>`_

*Older versions of each library may be usable, but the listed versions have been tested.*

VALET on UD clusters
--------------------

If working on a UD cluster, a pre-built copy of the latest pCI package is available via VALET:

.. code-block:: 

   vpkg_require pci

This command automatically configures your environment to include all necessary libraries to run pCI: OpenMPI with an Intel Fortran compiler, as well as latest version of Python to run auxiliary scripts. If using the pre-built pCI distribution, you do not have to compile anything and can ignore the rest of this page.

If working on another HPC cluster, a similar method of setting up your environment may already be available using ``module load <package>``. Please check your cluster's documentations for more information. 

Obtaining the source code
-------------------------
Users can download the latest version of the pCI code package from our `GitHub repository <https://github.com/ud-pci/pCI>`_.

via direct download
~~~~~~~~~~~~~~~~~~~

From the GitHub repository, click on the green ``<> Code`` button towards the top, then click ``Download ZIP``. You can then unzip the downloaded file to obtain the source codes.

via git
~~~~~~~

From the command line, you can clone the latest version:

.. code-block:: 

   git clone https://github.com/ud-pci/pCI.git

Compiling with CMake
--------------------

The codes are built using CMake and the ``CMakeLists.txt`` files in the root and ``/src`` directories. The following are some example builds on the DARWIN computing cluster.

A standard build can be done:

.. code-block:: 

   $ vpkg_require cmake openmpi/4.1.0:intel-2020
   $ git clone https://github.com/ud-pci/pCI.git pci
   $ cd pci
   $ mkdir build
   $ cd build
   $ FC=mpifort cmake ..
      :
   $ make
   $ make install


A Debug build can be done:

.. code-block:: 

   $ vpkg_require cmake openmpi/4.1.0:intel-2020
   $ git clone https://github.com/ud-pci/pCI.git pci
   $ cd pci
   $ mkdir build-debug
   $ cd build-debug
   $ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200402-debug ..
      :
   $ make
   $ make install

An optimized build demands a little more:

.. code-block:: 

   $ vpkg_require cmake openmpi/4.1.0:intel-2020
   $ git clone https://github.com/ud-pci/pCI.git pci
   $ cd pci
   $ mkdir build-opt
   $ cd build-opt
   $ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200317-opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large   -xHost -m64" ..
      :
   $ make
   $ make install
