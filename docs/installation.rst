Installation
============

Requirements
------------
The following software libraries and tools are required to compile pCI:

* A Fortran compiler: GNU Fortran ``gfortran`` v12.2+, Intel Fortran Classic ``ifort`` v2020u4+, or LLVM-based Intel Fortran ``ifx``
* `CMake v3.13+ <https://cmake.org/download>`_
* `Python v3.9+ <https://www.python.org/downloads/>`_
* LAPACK and BLAS, or Intel Math Kernel Library (MKL) *(optional — required for* ``pconf``*,* ``ine``*, and* ``pol``*)*
* `OpenMPI v4.1+ <https://www-lb.open-mpi.org/software/ompi>`_ *(optional — required only for MPI-parallel programs)*

*Older versions may work but the listed versions have been tested.*

On HPC clusters
---------------
If working on an HPC cluster, environment modules may be available using ``module load <package>`` for each required library. Please check your cluster's documentation for more information.

Obtaining the source code
-------------------------
Users can download the latest version of the pCI code package from our `GitHub repository <https://github.com/ud-pci/pCI>`_.

Via direct download
~~~~~~~~~~~~~~~~~~~
From the GitHub repository, click on the green ``<> Code`` button towards the top, then click ``Download ZIP``. You can then unzip the downloaded file to obtain the source codes.

Via git
~~~~~~~

.. code-block::

   git clone https://github.com/ud-pci/pCI.git

Compiling with CMake
--------------------

Basic build
~~~~~~~~~~~

.. code-block::

   $ cd pCI
   $ mkdir build && cd build
   $ FC=ifort cmake -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
   $ make
   $ make install

For the LLVM-based Intel compiler, replace ``FC=ifort`` with ``FC=ifx``.

If LAPACK/MKL is not available, ``pconf``, ``ine``, and ``pol`` will be skipped automatically.

If OpenMPI is not available, omit ``-DMPI_HOME`` and CMake will automatically skip the MPI-parallel programs:

.. code-block::

   $ FC=ifort cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..

Build options
~~~~~~~~~~~~~

**Debug build:**

.. code-block::

   $ FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../debug ..

**Optimized build:**

.. code-block::

   $ FC=ifort cmake -DCMAKE_BUILD_TYPE=Release -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..

**Double-precision two-electron and isotope-shift integrals** (increases memory usage but improves numerical accuracy):

.. code-block::

   $ FC=ifort cmake -DUSE_DP_INTEGRALS=ON -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..

**ScaLAPACK diagonalization in pconf** (uses distributed-memory ScaLAPACK instead of serial LAPACK for large matrix diagonalization; requires OpenMPI and LAPACK/MKL):

.. code-block::

   $ FC=ifort cmake -DUSE_SCALAPACK=ON -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..

Troubleshooting
~~~~~~~~~~~~~~~

**-DMPI_HOME not working:** Specifying the executables directly is more robust:

.. code-block::

   $ FC=ifort cmake \
       -DMPI_Fortran_COMPILER=${OPENMPI_PREFIX}/bin/mpifort \
       -DMPIEXEC_EXECUTABLE=${OPENMPI_PREFIX}/bin/mpiexec \
       -DCMAKE_INSTALL_PREFIX=$(pwd)/../ \
       ..

**Intel MPI conflict:** If Intel MPI is also installed (e.g. via Intel oneAPI), CMake may detect it instead of OpenMPI. Point CMake to the correct executables explicitly:

.. code-block::

   $ FC=ifort cmake \
       -DMPI_Fortran_COMPILER=/path/to/openmpi/bin/mpifort \
       -DMPIEXEC_EXECUTABLE=/path/to/openmpi/bin/mpiexec \
       -DCMAKE_INSTALL_PREFIX=$(pwd)/../ \
       ..

If CMake still detects Intel MPI after correcting the paths, clear the CMake cache before reconfiguring:

.. code-block::

   $ rm -rf build/*
