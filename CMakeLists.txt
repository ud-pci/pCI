cmake_minimum_required(VERSION 3.12)

project(pci 
        VERSION 1.0.5
        DESCRIPTION "pCI software package"
        LANGUAGES Fortran)
include(GNUInstallDirs)

# Need Intel Fortran compiler
if (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" AND NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    message(FATAL_ERROR "Intel Fortran compiler (ifort or ifx) not found. Please ensure FC is set to 'ifort' or 'ifx'.")
else()
    message(STATUS "Intel Fortran compiler found: ${CMAKE_Fortran_COMPILER}")
endif()

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

# Add default Intel compile/link flags:
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set(MKL_PARALLEL_FLAGS -qopenmp -mkl=parallel)
    set(MKL_PARALLEL_LINK_FLAGS -qopenmp -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
    add_compile_options(-qopenmp -mkl=sequential)
    add_link_options(-qopenmp -mkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    set(MKL_PARALLEL_FLAGS -qopenmp -qmkl=parallel)
    set(MKL_PARALLEL_LINK_FLAGS -qopenmp -qmkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
    add_compile_options(-qopenmp -qmkl=sequential)
    add_link_options(-qopenmp -qmkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
endif()
# The FindMPI returns linker flags as a space-delimited string, but target_link_options() needs a list:
string(REGEX REPLACE "[ \t\r\n]" ";" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_subdirectory(src)
