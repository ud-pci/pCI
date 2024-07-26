# Installation
The pCI software package was designed to work primarily on 

## Required libraries
In order to compile pCI the following software libraries and tools are required: 

* Intel Fortran compiler 2020u4+
* CMake v3.12+ 
* OpenMPI v4.1+ 
* Python v3.x 

*Older versions of each library may be usable, but the listed versions have been tested.*

## VALET on UD clusters
If working on a UD cluster, a pre-built copy of the latest pCI package is available via VALET:

```
vpkg_require pci
```

This command automatically configures your environment to include all necessary libraries to run pCI. If using the pre-built pCI distribution, you do not have to compile anything and can ignore the rest of this page.

If working on another HPC cluster, a similar method of setting up your environment may already be available using ``module load <package>``. Please check your cluster's documentations for more information. 


## Obtaining the source code
Users can download the latest version of the pCI code package from our GitHub repository.

### via direct download
From the GitHub repository, click on the green ``<> Code`` button towards the top, then click ``Download ZIP``. You can then unzip the downloaded file to obtain the source codes.

### via git
From the command line, you can clone the latest version:

```
git clone https://github.com/ud-pci/pCI.git
```

## Installing with VALET
VALET is a tool developed by Dr. Jeffrey Frey that is available on the UD clusters. VALET provides the ability to modify your environment without modifying your startup files like ``.bashrc`` and ``.bash_profile`` as it prevents jobs from failing by keeping a clean login environment. To install the custom ``pCI`` package using VALET, simply load with:
```
vpkg_require pci
```

## Building with CMake
The codes are built using CMake with the ``CMakeLists.txt`` files in the root and ``src`` directories. 

A general build can be done:
```
$ cd pCI
$ mkdir build
$ cd build
$ FC=mpifort cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ..
   :
$ make
$ make install
```

A Debug build can be done:
```
$ cd pCI
$ mkdir build-debug
$ cd build-debug
$ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200402-debug ..
   :
$ make
$ make install
```

An optimized build demands a little more:

```
$ cd ..
$ mkdir build-opt
$ cd build-opt
$ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200317-opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..
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
mpirun -n <nprocs> basc
mpirun -n <nprocs> conf
mpirun -n <nprocs> dtm
mpirun -n <nprocs> conf_pt
mpirun -n <nprocs> conf_lsj
```

For the polarizability code ``ine``, you will have to set the number of threads.
```
export OMP_NUM_THREADS=<nthreads>
ine
```