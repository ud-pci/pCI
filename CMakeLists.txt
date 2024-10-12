cmake_minimum_required(VERSION 3.12)

project(pci 
        VERSION 1.0.0
        DESCRIPTION "pCI software package"
        LANGUAGES Fortran)
include(GNUInstallDirs)

# Need Intel Fortran compiler
if (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    message(FATAL_ERROR "Intel Fortran compiler (ifort) not found. Please ensure FC is set to 'ifort'.")
else()
    message(STATUS "Intel Fortran compiler found: ${CMAKE_Fortran_COMPILER}")
endif()

# Add default Intel compile/link flags:
set(MKL_PARALLEL_FLAGS -qopenmp -mkl=parallel)
set(MKL_PARALLEL_LINK_FLAGS -qopenmp -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
add_compile_options(-qopenmp -mkl=sequential)
add_link_options(-qopenmp -mkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)

# Need MPI
find_package(MPI REQUIRED)
execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE MPI_VERSION_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (MPI_VERSION_OUTPUT MATCHES "Open MPI" OR MPI_VERSION_OUTPUT MATCHES "OpenRTE")
    message(STATUS "OpenMPI detected: ${MPI_VERSION_OUTPUT}")
else()
    message(FATAL_ERROR "OpenMPI not detected. Please ensure OpenMPI is installed and loaded, and that mpifort points to the OpenMPI version.")
endif()

# The FindMPI returns linker flags as a space-delimited string, but target_link_options() needs a list:
string(REGEX REPLACE "[ \t\r\n]" ";" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_subdirectory(src)
