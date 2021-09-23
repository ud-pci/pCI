# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.4] - 2021-09-23
- conf v0.3.13 - bug fix for MPI_AllReduce to same buffer with large number of cores

## [0.7.3] - 2021-09-22
- conf_lsj v0.1.2 - bug fix for conf_lsj to read correct eigenvectors

## [0.7.2] - 2021-09-22
- conf_lsj v0.1.1 - conf_lsj now asks for 2 input parameters rec1 and rec2 to select range of desired levels

## [0.7.1] - 2021-09-21
- basc has been parallelized with MPI using simple static workload distribution
- minor refactoring and text revisions

## [0.6.12] - 2021-09-18
- removed writing/reading of processors from HIJ and JJJ files
- added more code documentation
- initialization of global variables moved to module
- revamped F_J2 function to remove redundant comparisons by Rspq

## [0.6.11] - 2021-09-16
- implemented dynamic workload scheduling in conf_lsj
- added completion time message to ine
- added several timing messages to conf_lsj

## [0.6.10] - 2021-09-15
- bug fix in conf - fixed energies for case of Kl4 = 2
- added more code documentation to conf

## [0.6.9] - 2021-09-14
- restructed conf program - removed conf_aux module and combined with main conf program
- parallelized contruction of initial approximation in conf

## [0.6.8] - 2021-09-08
- bug fix in conf.f90 - Kl4 is now broadcasted

## [0.6.7] - 2021-09-03
- global variable Scr is now allocatable

## [0.6.6] - 2021-09-03
- bug fix in conf.f90 - fixed required memory count display before Davidson procedure

## [0.6.5] - 2021-09-03
- bug fix in ine.f90 - added rounding of Tj0 to prevent errors with NaN

## [0.6.4] - 2021-09-02
- bug fix in ine.f90 - fixed ReadJJJ subroutine to read CONF.JJJ file

## [0.6.3] - 2021-09-01
- bug fixes in matrix_io.f90 - ReadMatrix subroutine now reads nprocs.conf for # processors

## [0.6.2] - 2021-09-01
- bug fixes in "conf", allowing use of Kl=1 - reading of CONF.HIJ and CONF.JJJ files in conf

## [0.6.1] - 2021-08-31
- "conf" program now writes CONFp.HIJ and CONFp.JJJ instead of CONF.HIJ and CONF.JJJ
- "sort.py" program converts CONFp.HIJ and CONFp.JJJ into serial CONF.HIJ and CONF.JJJ

## [0.6.0] - 2021-08-30
- Initial import of new modernized "ine" prgram
- "ine" program has been modernized
- "ine" program now uses intel mkl subroutine 'zspsv'
- "dtm" program subroutine OpenFS has been updated to conf_pt_breit version
- "conf" program subroutine WriteMatrix has been updated to write separate file nprocs.conf containing number of processors used in conf calculation
- "sort.py" program has been updated to read nprocs.conf

## [0.5.0] - 2021-08-09
- Initial import of new parallelized "conf_lsj" program
- "conf_lsj" program has been modernized and parallelized
- "conf_lsj" program is now included in the build process of the parallel package.
- "conf_lsj" has been tested with serial version and all results are identical in output files
- "conf_lsj" now utilizes dynamic memory allocation

## [0.4.0] - 2021-07-29
- Revamp of FormH procedure to use split comparison and calculation stages to save memory
- Many small bug fixes regarding memory
- Progress bar fixed

## [0.3.6] - 2021-07-14
- Removed usage of MPI windows - suspected memory leak due to pointers not being deallocated
- FormH bug fix - initial workload distribution was missing a chunk of determinants

## [0.3.5] - 2021-07-13
- Minor bug fixes and progress text changes

## [0.3.4] - 2021-07-09
- Updated progress report of FormJ to include memory per core

## [0.3.3] - 2021-07-08
- Bug fix for FormH and FormJ subroutines

## [0.3.2] - 2021-07-07
- Updated progress report of FormH and FormJ to include memory progress
- Removed print statements of number of elements for each core

## [0.3.1] - 2021-06-30
- Bug fixes in FormH and FormJ new work schedulers
- Implemented progress report at increments of 10% for FormH

## [0.3.0] - 2021-06-23
- Major parallelism revamp: FormH and FormJ subroutines now use a work scheduler to distribute workloads
- Householder diagonalization parameters changed (tol: -103 to -1021; eps: -24 to -53)

## [0.2.4] - 2021-05-14
- "Rspq" bug fix: added diagnostics for dets to make sure orbitals are in correct order

## [0.2.3] - 2021-05-10
- "matrix_io" bug fix: compatibility fix for sortJJJ.py and sortHIJ.py
- new code "sort.py": python code that converts parallel CONF.JJJ or CONF.HIJ to serial format

## [0.2.2] - 2021-04-23
- "sint1" in dtm_aux.f90: bug fix

## [0.2.1] - 2021-04-20
- "mpi_wins" bug fix: added lines to allocate and broadcast Iarr for zero-cores 

## [0.2.0] - 2021-04-09
- initial import of new modernized "basc" program
- "basc" program has been modernized and included in the build process of the parallel package.
- "basc" has been tested with serial version and all results are identical in output files
- "basc" now utilizes dynamic memory allocation

## [0.1.0] - 2021-03-28
- kv=3 functionality fully parallelized

## [0.0.5] - 2021-03-23
- new module "matrix_io": implements parallel reading and writing of matrices (Hamiltonian and J^2)

## [0.0.4] - 2021-03-22
- new type "Matrix": encapsulate indices and values of matrix elements (Hamil and Jsq)
- new module "mpi_wins": implements creating and closing MPI windows/shared memory for basis set
- new module "mpi_utils": used for DARWIN's no-ucx variant of intel-2020 to bypass 1GB MPI message limit
- used intrinsic function PACK to remove all zero valued matrix elements from Hamiltonian and Jsquared
    - function PACK results in seg fault if not using 'ulimit -s unlimited'
- Fixed discrepancy of NumH and NumJ between serial and parallel versions
- reorganized timing calls
- removed unused error variables
- several minor text edits

## [0.0.3] - 2021-03-21
- dtm now uses conf_init module for reading CONF.INP
- revamped dtm's Input subroutine 
- moved one-electron operator functions to a separate module amp_ops
- several minor text edits

## [0.0.2] - 2021-03-20
- dtm now uses determinants module for determinant-based subroutines
- made arrays storing configurations (iconf1, iconf2) consistent between all programs
- several minor text edits
- descriptions added for several subroutines in davidson module

## [0.0.1] - 2021-03-17
Fixed issue with table of J being unordered.

## [0.0.0] - 2021-03-17
Original source import. 