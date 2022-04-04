# Compiling pCI

The CI/CI+MBPT/CI+all-order code package has only been tested on the Linux operating system.

In order to compile pCI the following software libraries and tools are required:
- Intel Fortran compiler.
- CMake build tool.
- (Optional) MPI library to run on high-performance computing clusters. The codes have only been tested with OpenMPI so far. 

## Obtaining the source code

The source codes can be cloned onto a local directory:
```
$ git clone https://github.com/ccheung93/pCI.git pCI
```

To switch to the devel branch:
```
$ git checkout devel
```

## Installing with VALET

VALET is a tool developed by Dr. Jeffrey Frey that is available on the UD clusters. VALET provides the ability to modify your environment without modifying your startup files like ```.bashrc``` and ```.bash_profile``` as it prevents jobs from failing by keeping a clean login environment. To install the ```pCI``` package using VALET, simply load with:
```
vpkg_require pci
```

The programs are then loaded into the environment.
To run the programs, just run the command for program:
```
add
ine
```
```
mpirun -n <nprocs> conf
mpirun -n <nprocs> basc
mpirun -n <nprocs> dtm
mpirun -n <nprocs> conf_pt
```

## Building with CMake

The codes are built using CMake with the 'CMakeLists.txt' file and OpenMPI with Intel compiler. 
To load these modules on Caviness:
```
vpkg_require cmake openmpi/4.0.2:intel
```
To load these modules on Darwin:
```
vpkg_require cmake openmpi/4.1.0:intel-2020
```

A general build can be done:
```
$ cd pCI
$ mkdir build
$ cd build
$ FC=mpifort cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../ ../src/
   :
$ make
$ make install
```

A Debug build can be done:
```
$ cd pCI
$ mkdir build-debug
$ cd build-debug
$ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200402-debug ../src/
   :
$ make
$ make install
```

An optimized build demands a little more:

```
$ cd ..
$ mkdir build-opt
$ cd build-opt
$ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200317-opt -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ../src/
   :
$ make
$ make install
```