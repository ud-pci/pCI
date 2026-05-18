Installation
============

Required libraries
------------------
In order to compile pCI the following software libraries and tools are required: 

* `Intel Fortran compiler 2020u4+ <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html>`_
* `CMake v3.13+ <https://cmake.org/download>`_
* `OpenMPI v4.1+ <https://www-lb.open-mpi.org/software/ompi>`_
* `Python v3.x <https://www.python.org/downloads/>`_

*Older versions of each library may be usable, but the listed versions have been tested.*


On HPC clusters
---------------
If working on an HPC cluster, environment modules may be available using ``module load <package>`` for each required library. Please check your cluster's documentations for more information. 

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

A general build can be done:

.. code-block:: 

   $ cd pCI
   $ mkdir build
   $ cd build
   $ FC=ifort cmake -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
      :
   $ make
   $ make install


A Debug build can be done:

.. code-block:: 

   $ cd pCI
   $ mkdir build-debug
   $ cd build-debug
   $ FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../debug ..
      :
   $ make
   $ make install

An optimized build demands a little more:

.. code-block:: 

   $ cd pCI
   $ mkdir build-opt
   $ cd build-opt
   $ FC=ifort cmake -DCMAKE_BUILD_TYPE=Release -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..
      :
   $ make
   $ make install
