# Installation
The pCI software package was designed to work primarily on HPC platforms running a Linux operating system.
It has also been tested on smaller personal laptops running Microsoft Windows with WSL2.

## Requirements
The following software libraries and tools are required to compile pCI:

* A Fortran compiler: GNU Fortran "gfortran" v12.2+, Intel Fortran Classic "ifort" v2020u4+, or LLVM-based Intel Fortran "ifx"
* CMake v3.13+
* Python v3.9+
* LAPACK and BLAS, or Intel Math Kernel Library (MKL) *(optional — required for* ``pconf``*,* ``ine``*, and* ``pol``*)*
* OpenMPI v4.1+ *(optional — required only for MPI-parallel programs)*

*Older versions may work but the listed versions have been tested.*

## Obtaining the source code

### Via direct download
From the GitHub repository, click the green ``<> Code`` button, then click ``Download ZIP`` and unzip the file.

### Via git
```
git clone https://github.com/ud-pci/pCI.git
```

## Building with CMake

### Basic build
```
$ cd pCI
$ mkdir build && cd build
$ FC=ifort cmake -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
$ make
$ make install
```

For the LLVM-based Intel compiler, replace ``FC=ifort`` with ``FC=ifx``.

If LAPACK/MKL is not available, ``pconf``, ``ine``, ``pol`` will be skipped automatically.

If OpenMPI is not available, omit ``-DMPI_HOME`` and CMake will automatically skip the MPI-parallel programs and build only the serial programs:
```
$ FC=ifx cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
```

### Build options

**Debug build:**
```
$ FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../debug ..
```

**Optimized build:**
```
$ FC=ifort cmake -DCMAKE_BUILD_TYPE=Release -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..
```

**Double-precision two-electron and isotope-shift integrals** (increases memory usage but improves numerical accuracy):
```
$ FC=ifort cmake -DUSE_DP_INTEGRALS=ON -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
```

**ScaLAPACK diagonalization in pconf** (uses distributed-memory ScaLAPACK instead of serial LAPACK for large matrix diagonalization; requires OpenMPI and LAPACK/MKL):
```
$ FC=ifort cmake -DUSE_SCALAPACK=ON -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
```

### Troubleshooting

**`-DMPI_HOME` not working:** ``-DMPI_HOME`` is a hint to CMake's MPI finder and is not always reliable — CMake may still pick up the wrong MPI installation. Specifying the executables directly is more robust:
```
FC=ifort cmake \
  -DMPI_Fortran_COMPILER=${OPENMPI_PREFIX}/bin/mpifort \
  -DMPIEXEC_EXECUTABLE=${OPENMPI_PREFIX}/bin/mpiexec \
  -DCMAKE_INSTALL_PREFIX=$(pwd)/../ \
  ..
```

**Incompatible compilers or modules:** Specify the compilers explicitly:
```
FC=ifort cmake -DCMAKE_Fortran_COMPILER=ifort -DMPI_Fortran_COMPILER=mpifort -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
```

**Intel MPI conflict:** If Intel MPI is also installed (e.g. via Intel oneAPI), CMake may detect it instead of OpenMPI. pCI requires OpenMPI, so point CMake to the correct executables explicitly:
```
FC=ifx cmake \
  -DMPI_Fortran_COMPILER=/path/to/openmpi/bin/mpifort \
  -DMPIEXEC_EXECUTABLE=/path/to/openmpi/bin/mpiexec \
  -DCMAKE_INSTALL_PREFIX=$(pwd)/../ \
  ..
```
If CMake still detects Intel MPI after correcting the paths, clear the CMake cache before reconfiguring:
```
rm -rf build/*
```

## Running the programs

**Serial programs:**
```
hfd
bass
add
sort
check_matrices
```

**MPI-parallel programs** (specify number of processors with ``-n``):
```
mpirun -n <nprocs> pbasc
mpirun -n <nprocs> pconf
mpirun -n <nprocs> pdtm
mpirun -n <nprocs> conf_pt
mpirun -n <nprocs> conf_lsj
mpirun -n <nprocs> formy
```

**OpenMP-parallel programs** (set number of threads before running):
```
export OMP_NUM_THREADS=<nthreads>
ine
pol
```

## HPC clusters

### General
If working on an HPC cluster, the required libraries may be available as environment modules:
```
module load <package>
```
Please check your cluster's documentation for available modules.

### UD clusters (VALET)
On UD clusters, a pre-built copy of the latest pCI package is available via VALET, which configures your environment without modifying startup files like ``.bashrc``.

To load pCI:
```
vpkg_require pci
```
This automatically loads OpenMPI, an Intel Fortran compiler, and Python. If using the pre-built distribution, no compilation is needed.
