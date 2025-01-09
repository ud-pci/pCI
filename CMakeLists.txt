cmake_minimum_required(VERSION 3.13)

project(pci 
        VERSION 1.1.1
        DESCRIPTION "pCI software package"
        LANGUAGES Fortran)
include(GNUInstallDirs)

# Need a Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    message(STATUS "Intel Fortran compiler found: ${CMAKE_Fortran_COMPILER}")
    message(STATUS "Intel Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "GNU Fortran compiler found: ${CMAKE_Fortran_COMPILER}")
    message(STATUS "GNU Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
    message(WARNING "Optional programs for QED and all-order are not compilable with gfortran")
else()
    message(FATAL_ERROR "Intel or GNU Fortran compiler was not detected.")
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

# Need BLAS and LAPACK
# Find MKL
find_package(MKL QUIET)
if (MKL_FOUND)
    message(STATUS "MKL libary detected: ${MKL_LIBRARIES}")
else()
    # If find_package can't find MKL, manually check for MKL
    set(MKL_ROOT $ENV{MKLROOT})
    
    if (MKL_ROOT)
        set(MKL_INCLUDE_DIR "${MKL_ROOT}/include")
        set(MKL_LIB_DIR "${MKL_ROOT}/lib/intel64")
        message(STATUS "MKL library detected: ${MKL_LIB_DIR}")
    else()
        message(STATUS "MKL library not detected. Falling back to BLAS and LAPACK.")

        # Find BLAS and LAPACK
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        if (BLAS_FOUND AND LAPACK_FOUND)
            message(STATUS "BLAS library detected: ${BLAS_LIBRARIES}")
            message(STATUS "LAPACK library detected: ${LAPACK_LIBRARIES}")
        else()
            message(FATAL_ERROR "BLAS and LAPACK not detected. Please ensure BLAS and LAPACK are installed.")
        endif()

        # Find ScaLAPACK
        list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
        find_package(ScaLAPACK QUIET)
        if (ScaLAPACK_FOUND)
            message(STATUS "ScaLAPACK detected: ${ScaLAPACK_LIBRARIES}")
        else()
            message(FATAL_ERROR "ScaLAPACK not detected.")
        endif()
    endif()
endif()

# Add default Intel compile/link flags:
if (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
    set(MKL_PARALLEL_FLAGS -qopenmp -mkl=parallel)
    set(MKL_PARALLEL_LINK_FLAGS -qopenmp -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
    add_compile_options(-qopenmp -mkl=sequential)
    add_link_options(-qopenmp -mkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
elseif (${CMAKE_Fortran_COMPILER} MATCHES "ifx")
    set(MKL_PARALLEL_FLAGS -qopenmp -qmkl=parallel)
    set(MKL_PARALLEL_LINK_FLAGS -qopenmp -qmkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
    add_compile_options(-qopenmp -qmkl=sequential)
    add_link_options(-qopenmp -qmkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
elseif (${CMAKE_Fortran_COMPILER} MATCHES "gfortran")
    add_compile_options(-fopenmp -lscalapack -lblas -llapack)
    add_link_options(-fopenmp -lscalapack -lblas -llapack)
endif()
# The FindMPI returns linker flags as a space-delimited string, but target_link_options() needs a list:
string(REGEX REPLACE "[ \t\r\n]" ";" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_subdirectory(src)
