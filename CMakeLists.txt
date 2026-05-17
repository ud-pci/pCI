cmake_minimum_required(VERSION 3.13)

project(pci 
        VERSION 1.2.0
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

# Build options
option(USE_DP_INTEGRALS "Use double-precision integrals" OFF)
option(USE_SCALAPACK "Use ScaLAPACK for matrix diagonalization in pconf" OFF)

if (USE_DP_INTEGRALS)
    message(STATUS "Using double-precision integrals")
endif()
if (USE_SCALAPACK)
    message(STATUS "Using ScaLAPACK for pconf diagonalization")
endif()

# Need MPI (optional — serial-only programs are built without it)
find_package(MPI)
if(MPI_FOUND)
    execute_process(
        COMMAND ${MPIEXEC_EXECUTABLE} --version
        OUTPUT_VARIABLE MPI_VERSION_OUTPUT
        ERROR_VARIABLE MPI_VERSION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (MPI_VERSION_OUTPUT MATCHES "Open MPI" OR MPI_VERSION_OUTPUT MATCHES "OpenRTE")
        message(STATUS "OpenMPI detected: ${MPI_VERSION_OUTPUT}")
    else()
        message(WARNING "Non-OpenMPI MPI detected. pCI requires OpenMPI — MPI-parallel programs may not build or run correctly.")
    endif()
    string(REGEX REPLACE "[ \t\r\n]" ";" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
else()
    message(WARNING "OpenMPI not found. MPI-parallel programs (pbasc, pconf, conf_lsj, conf_pt, pdtm, formy) will not be built.")
endif()

# Need BLAS and LAPACK
# Find MKL
find_package(MKL QUIET)
if (MKL_FOUND)
    message(STATUS "MKL library detected: ${MKL_LIBRARIES}")
else()
    # If find_package can't find MKL, manually check for MKL
    set(MKL_ROOT $ENV{MKLROOT})
    
    if (MKL_ROOT)
        set(MKL_INCLUDE_DIR "${MKL_ROOT}/include")
        if (EXISTS "${MKL_ROOT}/lib/intel64")
            set(MKL_LIB_DIR "${MKL_ROOT}/lib/intel64")
        else()
            set(MKL_LIB_DIR "${MKL_ROOT}/lib")
        endif()
        message(STATUS "MKL library detected: ${MKL_LIB_DIR}")
    else()
        message(STATUS "MKL library not detected. Falling back to BLAS and LAPACK.")

        # Find BLAS and LAPACK
        find_package(BLAS)
        find_package(LAPACK)
        if (BLAS_FOUND)
            message(STATUS "BLAS library detected: ${BLAS_LIBRARIES}")
        endif()
        if (LAPACK_FOUND)
            message(STATUS "LAPACK library detected: ${LAPACK_LIBRARIES}")
        endif()

        # Find ScaLAPACK only if needed
        if (USE_SCALAPACK)
            list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
            find_package(ScaLAPACK REQUIRED)
            if (ScaLAPACK_FOUND)
                message(STATUS "ScaLAPACK detected: ${ScaLAPACK_LIBRARIES}")
            endif()
        endif()
    endif()
endif()

# Compiler-specific compile/link flags:
add_library(mkl_sequential INTERFACE)
add_library(mkl_parallel INTERFACE)

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    add_compile_options(-qopenmp -mkl=sequential)
    add_link_options(-qopenmp)
    target_link_options(mkl_sequential INTERFACE -mkl=sequential)
    target_link_options(mkl_parallel INTERFACE -mkl=parallel)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    add_compile_options(-qopenmp -qmkl=sequential)
    add_link_options(-qopenmp)
    target_link_options(mkl_sequential INTERFACE -qmkl=sequential)
    target_link_options(mkl_parallel INTERFACE -qmkl=parallel)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    add_compile_options(-fopenmp)
    add_link_options(-fopenmp)
    if (MKL_LIB_DIR)
        find_library(MKL_GF_LP64     mkl_gf_lp64     PATHS ${MKL_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
        find_library(MKL_SEQUENTIAL  mkl_sequential   PATHS ${MKL_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
        find_library(MKL_GNU_THREAD  mkl_gnu_thread   PATHS ${MKL_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
        find_library(MKL_CORE        mkl_core         PATHS ${MKL_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
        target_link_libraries(mkl_sequential INTERFACE ${MKL_GF_LP64} ${MKL_SEQUENTIAL} ${MKL_CORE} -lpthread -lm -ldl)
        target_link_libraries(mkl_parallel   INTERFACE ${MKL_GF_LP64} ${MKL_GNU_THREAD} ${MKL_CORE} -lgomp -lpthread -lm -ldl)
    elseif (MKL_LIBRARIES)
        target_link_libraries(mkl_sequential INTERFACE ${MKL_LIBRARIES})
        target_link_libraries(mkl_parallel   INTERFACE ${MKL_LIBRARIES})
    elseif (BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries(mkl_sequential INTERFACE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
        target_link_libraries(mkl_parallel   INTERFACE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    endif()
endif()

if (MKL_FOUND OR MKL_ROOT OR (BLAS_FOUND AND LAPACK_FOUND))
    set(LAPACK_AVAILABLE TRUE)
else()
    set(LAPACK_AVAILABLE FALSE)
    message(WARNING "BLAS/LAPACK/MKL not found. ine, pol, and pconf will not be built.")
endif()

if (USE_SCALAPACK)
    if (NOT MPI_FOUND)
        message(FATAL_ERROR "USE_SCALAPACK requires OpenMPI, which was not found.")
    endif()
    if (NOT LAPACK_AVAILABLE)
        message(FATAL_ERROR "USE_SCALAPACK requires BLAS/LAPACK/MKL, which was not found.")
    endif()
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_subdirectory(src)
