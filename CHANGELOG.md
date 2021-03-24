# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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