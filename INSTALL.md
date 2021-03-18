# Compiling pCI

The CI/CI+MBPT/CI+all-order code package has only been tested on the Linux operating system.

In order to compile pCI the following software libraries and tools are required:
- Intel Fortran compiler.
- CMake build tool.
- (Optional) MPI library to run on high-performance computing clusters. The codes have only been tested with OpenMPI so far. 

## Obtaining the source code

The source codes can be cloned onto a local directory:
```
$ mkdir pCI
$ git clone https://github.com/ccheung93/pCI.git pCI
```

## Building with CMake

The codes are built using the 'CMakeLists.txt' file. 
A Debug build can be done:
```
$ cd pCI
$ mkdir build-debug
$ cd build-debug
$ vpkg_require cmake openmpi/4.1.0:intel-2020
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