Installation
============

Required libraries
------------------
In order to compile pCI the following software libraries and tools are required:

* Fortran compiler.
* CMake build tool.
* MPI library to run on high-performance computing clusters. The codes have only been tested with OpenMPI so far.
* Python v3.xx

VALET on UD clusters
--------------------

If working on a UD cluster, a pre-built copy of the latest pCI package can be readily obtainable via VALET:

.. code-block:: 

   vpkg_require pci

This command automatically configures your environment to include all necessary libraries to run pCI: an Intel Fortran compiler and OpenMPI library, as well as latest version of Python to run auxiliary scripts. If using the pre-built pCI distribution, you do not have to compile anything and can ignore the rest of this page.

Obtaining the source code
-------------------------
Users can download the latest version of the pCI code package from

.. code-block:: 

   https://github.com/ud-pci/pCI

From the command line, you can clone the latest version:

.. code-block:: 

   git clone https://github.com/ud-pci/pCI.git pci

Compiling with CMake
--------------------

The codes are built using the ``CMakeLists.txt`` files. The following are some example builds on the DARWIN computing cluster.

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
