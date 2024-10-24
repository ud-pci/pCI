# Installation
The pCI software package was designed to work primarily on 

## Required libraries
In order to compile pCI the following software libraries and tools are required: 

* Intel Fortran compiler 2020u4+
* CMake v3.12+ 
* OpenMPI v4.1+ 
* Python v3.x 

*Older versions of each library may be usable, but the listed versions have been tested.*

## On HPC clusters
If working on an HPC cluster, environment modules may be available using ``module load <package>`` for each required library. Please check your cluster's documentations for more information. 

## VALET on UD clusters
If working on a UD cluster, a pre-built copy of the latest pCI package is available via VALET. VALET is a tool developed by Dr. Jeffrey Frey that provides the ability to modify your environment without modifying your startup files like ``.bashrc`` and ``.bash_profile`` as it prevents jobs from failing by keeping a clean login environment. 

To load the pCI package to your environment, simply run the following command:
```
vpkg_require pci
```

This command automatically configures your environment to include all necessary libraries to run pCI: OpenMPI with an Intel Fortran compiler, as well as latest version of Python to run auxiliary scripts. If using the pre-built pCI distribution, you do not have to compile anything and can ignore the rest of this page.

## Obtaining the source code
Users can download the latest version of the pCI code package from our GitHub repository.

### via direct download
From the GitHub repository, click on the green ``<> Code`` button towards the top, then click ``Download ZIP``. You can then unzip the downloaded file to obtain the source codes.

### via git
From the command line, you can clone the latest version:

```
git clone https://github.com/ud-pci/pCI.git
```

## Building with CMake
The codes are built using CMake with the ``CMakeLists.txt`` files in the root and ``src`` directories. The following are some example builds on the DARWIN computing cluster.

A general build can be done:
```
$ cd pCI
$ mkdir build
$ cd build
$ FC=ifort cmake -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
   :
$ make
$ make install
```

If problems related to incompatible compilers or modules, one can specify:
```
FC=ifort cmake -DCMAKE_Fortran_COMPILER=ifort -DMPI_Fortran_COMPILER=mpifort -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
```

A Debug build can be done:
```
$ cd pCI
$ mkdir build-debug
$ cd build-debug
$ FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../debug ..
   :
$ make
$ make install
```

An optimized build demands a little more:
```
$ cd ..
$ mkdir build-opt
$ cd build-opt
$ FC=ifort cmake -DCMAKE_BUILD_TYPE=Release -DMPI_HOME=${OPENMPI_PREFIX} -DCMAKE_INSTALL_PREFIX=$(pwd)/../opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..
   :
$ make
$ make install
```

## Running the programs
To run the programs, simply run the executables for the programs:
```
hfd
bass
add
```

For the parallel codes, you will have to specify the number of processors.

```
mpirun -n <nprocs> pbasc
mpirun -n <nprocs> pconf
mpirun -n <nprocs> pdtm
mpirun -n <nprocs> conf_pt
mpirun -n <nprocs> conf_lsj
```

For the polarizability code ``pol``, you will have to set the number of threads.
```
export OMP_NUM_THREADS=<nthreads>
pol
```