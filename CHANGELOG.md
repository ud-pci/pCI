# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.11.10] - 2022-11-21
- conf v5.10

## [0.11.9] - 2022-08-17
- ine v1.19 - fixed numerical errors in calculating a_1 value

## [0.11.8] - 2022-08-10
- conf v5.9 - bug fix for error occuring in allocation of arrays for LSJ
- conf.f90 - removed unnecessary allocations and broadcasts to arrays storing integrals in InitLSJ and AllocateLSJArrays
- conf.f90 - removed devel comments from previous version of conf
- conf.f90 - increased integer formatting of j for configurations to 2
- conf.f90 - increased integer formatting for energy levels to 3
- davidson.f90 - increased integer formatting for energy levels to 3

## [0.11.7] - 2022-08-05
- basc v2.2 - addition of Andrey's code for core Breit contributions
- basc.f90 - added import of breit_int subroutine in subroutine Core
- basc.f90 - new variable Ebcore for Breit contribution to core energy and its calculation
- basc.f90 - core energy output now includes Breit 
- basc.f90 - changed formatting for Rint2Table

## [0.11.6] - 2022-08-03
- conf v5.8 - fix for weight subroutines when missing lines and counters for non-rel configs in CONF.INP
- conf_init.f90 - removed counting of configurations (Nrnrc) and moved to subroutine Init
- conf.f90 - new code to count number of non-relativistic configurations

## [0.11.5] - 2022-07-29
- ine v1.18 - bug fix for integer overflow in Mxmpy
- ine.f90 - changed type for variable i8 to kind int64

## [0.11.4] - 2022-07-28
- ine v1.17 - implemented 2-step iteration method
- ine.f90 - new key value for kIters = 2 allows for 2-step iteration method
- ine.f90 - new key N_it4 defines number of 2-step iterations

## [0.11.3] - 2022-07-25
- conf v5.7 - bug fix for Kv = 3
- davidson.f90 - re-added broadcast of Z1 in initial stage of FormB0
- conf.f90, davidson.f90, formj2.f90 - reverted generic intrisic functions (abs, sqrt) back to double-type intrinsic functions (dabs, dsqrt)
- formj2.f90 - counter for J is now set by counter1 variable instead of size(Jsq%val)

## [0.11.2] - 2022-07-25
- create_add_inp.py - if no BASS.INP is present, the list of subshells and occupation numbers will now be ordered starting from basic configurations, then ascending by principal quantum number 

## [0.11.1] - 2022-07-23
- added initial developmental version of create_conf_inp.py to py-lib folder

## [0.11.0] - 2022-07-22
- added python scripts to new py-lib folder
- basis.py - python script that generates input files used to create basis sets
- create_add_yml.py - python script that creates yaml-formatted input file used to create ADD.INP via create_add_inp.py script
- create_add_inp.py - python script that creates ADD.INP using yaml-formatted input file
- orbitals.py - python library of functions related to working with orbitals
- read_bass.py - python library of functions related to reading BASS.INP
- get_excitations.py - python library of functions related to generating excited configurations
- get_atomic_data.py - python library of functions scraping atomic data tables online
- add.yml - sample yaml input file using Sr
- basis.yml - sample yaml input file using Ac
- moved sort.py to new py-lib folder

## [0.10.6] - 2022-07-21
- tm v3.2 - addition of transition rates to E1, E2, M1 tables
- dtm.f90 - removed cases E_i -> E_i from M1, E2 tables
- dtm.f90 - new variables wl for wavelengths and tr for transition rates

## [0.10.5] - 2022-07-20
- ine v1.15 - bug fix for E2 polarizabilities
- ine.f90 - added missing allocation statement in Ort subroutine
- ine.f90 - removed calculation for negative parity for E2 polarizabilities
- ine.f90 - remove xlamb_next from last iteration

## [0.10.4] - 2022-07-10
- dm v3.0 - implemented new tables for g-factors, A_hfs and B_hfs
- dtm.f90 - rewrote SetKeys subroutine to open files depending on dm/tm and which keys are used
- dtm.f90 - broadcasting of Iarr now uses mpi_utils subroutine BroadcastI
- dtm.f90 - added code to get terms from CONFSTR.RES 
- dtm.f90 - added code to write tables for g-factors, A_hfs and B_hfs in DM regime

## [0.10.3] - 2022-07-07
- conf v5.6 - bug fixes and CONFFINAL.RES table changes
- conf.f90 - CONF.ENG and CONF.LVL are now written correctly from a new file
- conf.f90 - increased allocation of weight strings by factor of 5
- conf.f90 - change greek symbol delta in tables to DEL and energies EV to E_n (a.u.) in CONFFINAL.RES
- conf_init.f90 - bug fix for reading blank lines between configurations from result of merge_ci or con_cut

## [0.10.2] - 2022-07-01
- conf v5.5 - bug fixes
- conf.f90 - removed unused variables and arrays
- conf.f90 - changed argument of WriteMatrix to be Type(Matrix) instead of 3 arrays
- conf.f90 - added EXTERNAL declarations for BLACS_GRIDEXIT, BLACS_EXIT
- conf.f90 - removed temporary Jsq fixes
- formj2.f90 - added support for when core 0 obtains no non-zero matrix elements from FormJ
- formj2.f90 - changed argument of WriteMatrix to be Type(Matrix) instead of 3 arrays
- matrix_io.f90 - changed argument of WriteMatrix to be Type(Matrix) instead of 3 arrays
- ine.f90 - removed unused variables
- add.f90 - removed unused variables
- mpi_utils.f90 - added EXTERNAL declaration for GetEnv

## [0.10.1] - 2022-06-30
- tm v3.1 - bug fixes
- dtm.f90 - added conditionals to track terms if range does not start from level 1
- dtm.f90 - removed unused variables and arrays
- dtm.f90 - increased format digits for NSP to 6
- dtm.f90 - broadcasting of Iarr now uses mpi_utils subroutine BroadcastI

## [0.10.0] - 2022-06-29
- tm v3.0 - implemented new tables for improved data analysis
- dtm.f90: new type Key to encapsulate all keys for matrix elements
- dtm.f90: new subroutines SetKeys and CloseKeys to open and close files for writing individual matrix element tables
- dtm.f90: TM part now reads in 2 ranges: nterm1 to nterm1f and nterm2 to nterm2f
- dtm.f90: dtm.in reads list of matrix elements starting from 3rd line: 'E1', 'E1_L', 'E1_V', 'E2', 'E3', 'M1', 'M2', 'M3', 'EDM', 'PNC', 'MQM', 'AM'
- dtm.f90: added missing Read(Nsu) part for CONF.DET
- dtm.f90: removed progress report
- dtm.f90: added broadcast of new global variable nterm1f
- dtm.f90: added code to read configurations and terms from CONFSTR.RES
- dtm.f90: added loop over first range of energy levels
- dtm.f90: moved reading of wavefunctions to inside loop over first range of energy levels
- dtm.f90: added code to write tables of matrix elements

## [0.9.6] - 2022-06-16
- combined dtm_aux.f90 with dtm.f90

## [0.9.5] - 2022-06-16
- conf v5.4 - implemented writing of weights of configurations when CONF.XIJ is written
- conf.f90: new subroutine PrintWeightsDvdsn writes weights of configurations for each energy level to CONF.LVL every kXIJ iterations

## [0.9.4] - 2022-05-27
- ine v1.14 - bug fix for kl>0 where required arrays were not allocated

## [0.9.3] - 2022-05-26
- conf v5.3 - implemented writing of top weights of configurations to end of CONFLEVELS.RES
- conf.f90: wrote code to count the top configurations in CONFLEVELS.RES and list them in weighted order
- conf.f90: new Type(WeightTable) stores data for top weights
- conf.f90: moved energies in EV and Delta(cm^-1) to after term if LSJ is used in CONFFINAL.RES
- conf.f90: added a significant figure to gfactors
- conf.f90: gfactors are no longer written if J=0

## [0.9.2] - 2022-05-19
- conf v5.2 - implemented writing of terms if LSJ is used
- conf.f90: separated writing of files from calculation of weights
- conf.f90: changed Wsave, Wpsave, strcsave to dimension(nconfs,Nlv)
- conf.f90: new function term to write term based on LSJ values
- conf.f90: changed formatting of CONFFINAL.RES table to include term for LSJ
- conf.f90: fixed formatting for different size configuration strings in CONFFINAL.RES table

## [0.9.1] - 2022-05-14
- conf v5.1 - implemented writing of energies every N davidson iteration
- conf_init.f90: fixed out-of-bound error for array Nrnrc in subroutine ReadConfigurations
- conf.f90: new subroutine PrintEnergiesDvdsn that writes table of energies to file CONF.ENG whenever CONF.XIJ is written
- conf.f90: updated all mpi subroutines to use mpi_f08
- conf.f90: fixed out-of-bound error for array Wsave, Wpsave, strcsave in subroutine PrintWeights
- conf.f90: updated CONFLEVELS.RES file to print energies with more digits

## [0.9.0] - 2022-05-13
- conf v5.0 - implemented ikarus and conf_lsj into conf
- new energy table with weights (and optionally, LSJ) are written to file CONFFINAL.RES
- list of top contributing configurations to each energy level are written to file CONFLEVELS.RES
- conf_variables.f90: global variable kLSJ added to allow user to decide whether L, S, J are printed
- conf_variables.f90: global variable Nnr added to count number of non-relativistic configurations required by subroutine PrintWeights
- conf_variables.f90: global array Nrnrc added to store number of relativistic configurations in each non-relativistic configuration
- determinants.f90: deallocation statements for Nq and Nip have been removed for subsequent LSJ calculations
- determinants.f90: subroutine Rdet has been revamped to take in the name of the file to read determinants Iarr from
- conf_init.f90: added code in subroutine ReadConfigurations to count number of non-rel. configurations Nnr and number of rel. configurations in each nr configuration Nrnrc
- formj2.f90: deallocation statements for Jz and Nh have been removed for subsequent LSJ calculations
- conf.f90: added code to calculate LSJ if key kLSJ = 1
- conf.f90: removed deallocation statements for Nvc and Nc0 in subroutine DeAllocateFormHArrays
- conf.f90: new subroutine AllocateLSJArrays and InitLSJ added to allocate and initialize arrays required to calculate LSJ
- conf.f90: new subroutines lsj, calcLSJ, lsj_det, plus_s, plus_l, p0_s, p0_l, and plus_j added to calculate LSJ
- conf.f90: new array W2 added to store weights of non-relativistic configurations
- conf.f90: new array Wsave and Wpsave added to store the top weights and indices of top weights
- conf.f90: new code written to write a table of main configuration, S, L, J, gfactor, energy in eV, energy from ground state in cm-1, weight% for main configuration, secondary configuration and secondary weight if the weight for main configuration is < 70%, for each energy level
- conf.f90: new code written to write list of top weights and configurations for each energy level
- conf.f90: new function g_factor to calculate g-factors given L, S, J

## [0.8.6] - 2022-05-19
- conf v4.5 - bug fix for seg faults when writing intermediate CONF.XIJ caused by deallocation of Jsq after LAPACK subroutines

## [0.8.5] - 2022-04-04
- added serial modernized add program

## [0.8.4] - 2022-04-04
- conf v4.4 - merged Kl=3 and Kl=4 routines to single block
- added exception for when no additional configurations were added with key Kl=3
- added additional timing tests to file CONF.PRG

## [0.8.3] - 2022-03-28
- conf v4.3 - bug fixes for ArrB vectorization
- added determinant number in cntarray

## [0.8.2] - 2022-03-01
- conf v4.2 - changed BroadcastD to MPI_AllReduce 

## [0.8.1] - 2022-02-04
- conf v4.1 - implementation of parameterized derived type for matrices to switch between single and double precision
- added minimum value of matrix in type Matrix

## [0.8.0] - 2022-01-21
- conf v4.0 - implementation of Kl=3 - additional configurations
- added global variables Nc_prev and Nd_prev to conf_variables.f90
- added subroutines IVAccumultorContinue and RVAccumulatorContinue to vaccumulator.f90
- added condition for Kl=3 to cycle when constructing initial approximation
- added new conditional under FormH for Kl=3
- changed arguments of ReadMatrix and WriteMatrix to take in three 1d arrays instead of a matrix
- added conditional to read and write previous Nc and Nd in subroutines ReadMatrix and WriteMatrix 

## [0.7.56] - 2022-02-08
- ine v1.13 - increased limit of NumJ from integer to integer*8

## [0.7.56] - 2022-02-04
- ine v1.12 - increased limit of NumH from integer to integer*8

## [0.7.55] - 2022-01-13
- conf v3.46 - bug fix for Jsq%n for Kv=3

## [0.7.54] - 2022-01-11
- ine v1.11 - ine now accepts a list of wavelengths to be entered with different step sizes
- formatting of ine input file is now: Kl, Kli, Klf, N0, N2, nlambda, (xlambda1, xlambda2, xlambdastep), ...

## [0.7.53] - 2021-12-27
- ine v1.10 - small components of J-decomposition are no longer put to 0

## [0.7.52] - 2021-11-30
- changed type of vLen, vSize, vGrowBy in vaccumulator.f90 to int64

## [0.7.51] - 2021-11-26
- conf v3.45 - updated Rint/RintS errors if IPx was changed
- split memCalcReqs FormH message into static arrays and FormH arrays
- removed FormH memory from comparison stage of FormH
- added memory summary after comparison stage of FormH

## [0.7.50] - 2021-11-24
- conf v3.44 - updated memory routines for FormH

## [0.7.49] - 2021-11-22
- conf v3.43 - bug fix for temporary fix for Jsq%n changing during LAPACK ZSYEV subroutine

## [0.7.48] - 2021-11-10
- conf v3.42 - fixed integer overflow in mpi_utils subroutines
- increased data size of count_remain and i to integer(kind=int64)
- mpi_utils subroutine now take in int*8 argument for count

## [0.7.47] - 2021-11-10
- conf v3.41 - bug fix for mpi_utils subroutine BroadcastI for extra large integer arrays

## [0.7.46] - 2021-11-1
- ine v1.9 - fix for overflow in final result of RdcE1

## [0.7.45] - 2021-10-31
- ine v1.8 - increased number of significant figures in final result in RdcE1

## [0.7.44] - 2021-10-31
- ine v1.7 - fixed formatting error in RdcE1

## [0.7.43] - 2021-10-31
- conf v3.40 - added if allocated statement when deallocating arrays after Diag4

## [0.7.42] - 2021-10-30
- conf v3.39 - bug fix in conf resolving slow down of convergence during Davidson procedure

## [0.7.41] - 2021-10-27
- conf v3.38 - bug fix in conf when reading old version of CONF.GNT

## [0.7.40] - 2021-10-27
- conf_lsj v2.3 - added g-factors to final table

## [0.7.39] - 2021-10-26
- conf v3.37 - re-implemented Jeff's Broadcast subroutines for MPI_Bcast calls for 2d arrays

## [0.7.38] - 2021-10-25
- conf v3.36 - added memory check after allocating arrays for Davidson procedure

## [0.7.37] - 2021-10-24
- conf v3.35 - added exception to set Ngaunt if not in CONF.GNT
- added temporary fix for value of Jsq%n changing during LAPACK ZSYEV subroutine
- all programs - removed all unused variables and arrays determined by compiler flag "--check all --warn all"

## [0.7.36] - 2021-10-20
- conf v3.34 - removed reading of CONF.XIJ from PrintEnergies subroutine
- moved the printing of weights from PrintEnergies subroutine to a new PrintWeights subroutine
- cleaned up formatting of PrintEnergies and PrintWeight subroutine

## [0.7.35] - 2021-10-19
- dtm v2.1 - bug fix for Kout /= 0

## [0.7.34] - 2021-10-13
- conf v3.33 - added calculation of J for energy levels if writing CONF.XIJ file

## [0.7.33] - 2021-10-13
- removed dead code 'testwigner' used for testing from ine

## [0.7.32] - 2021-10-13
- global change: replaced parameter IPgnt with global variable Ngaunt read in from CONF.GNT
- basc v2.1 - number of tabulated gaunts Ngaunt is now written to CONF.GNT and read in by subsequent codes
- conf v3.32 - reads in Ngaunt from CONF.GNT
- conf_lsj v2.2 - reads in Ngaunt from CONF.GNT
- conf_pt v2.1 - reads in Ngaunt from CONF.GNT

## [0.7.31] - 2021-10-12
- conf v0.3.31 - moved Bcast of ax to end of loop

## [0.7.30] - 2021-10-12
- conf v0.3.30 - added conditional to optimize construction of initial approximation depending on Kl

## [0.7.29] - 2021-10-12
- conf v0.3.29 - fixed print message for number of conf-s in starting approximation

## [0.7.28] - 2021-10-10
- conf v0.3.28 - added mpi error handling for WriteMatrix and ReadMatrix subroutines

## [0.7.27] - 2021-10-08
- cleaned up dead code in some conf subroutines

## [0.7.26] - 2021-10-07
- conf v0.3.27 - moved averaging of diagonal over configurations to its own subroutine AvgDiag

## [0.7.25] - 2021-10-07
- conf v0.3.26 - moved initial diagonalization to its own subroutine DiagInitApprox

## [0.7.24] - 2021-10-07
- conf v0.3.25 - implemented ScaLAPACK subroutine PDSYEVD to replace DSYEV in initial diagonalization

## [0.7.23] - 2021-10-07
- conf v0.3.24 - reduced size of arrays E and Iconverge from IPlv to Nlv

## [0.7.22] - 2021-10-06
- conf v0.3.23 - replaced all Hould routines with LAPACK 'DSYEV'
- removed global D, D1 work arrays and replaced with local W work array
- initialize work array for diagonalizing P before iterative procedure
- bug fix for Kl4=2 - broadcast Nlv in FormB0 after reading CONF.XIJ

## [0.7.21] - 2021-10-06
- conf v0.3.22 - optimized DSYEV subroutine in Diag4

## [0.7.20] - 2021-10-04
- conf v0.3.21 - bug fix for Kv=3 - added MPI_Bcast of Z1 after diagonalization
- more code documentation in Davidson procedure
- refactored parts of Davidson procedure

## [0.7.19] - 2021-10-04
- conf v0.3.20 - changed MPI_Reduce of Z1 in Init4 to MPI_AllReduce, removing MPI_Bcast from FormB0

## [0.7.18] - 2021-10-01
- conf v0.3.19 - subroutine ReadMatrix now distributes total number of matrix elements equally across cores

## [0.7.17] - 2021-09-29
- conf v0.3.18 - bug fix for writing/reading CONF.HIJ/CONF.JJJ for large scale runs 
- changed disp and disps in matrix_io.f90 to kind=MPI_OFFSET_KIND
- added error message if matrix file doesn't exist

## [0.7.16] - 2021-09-29
- conf v0.3.17 - changed formatting, increasing range for Nc4 and Nd0 for initial approximation

## [0.7.15] - 2021-09-29
- conf v0.3.16 - replaced Hould subroutine during initial diagonalization with LAPACK 'dsyev' subroutine 
- removed redundant MPI_Bcast call of D1=Tk in subroutine WriteFinalXIj

## [0.7.14] - 2021-09-28
- conf v0.3.15 - added timing for initial diagonalization

## [0.7.13] - 2021-09-26
- conf_lsj v0.2.1 - include comparison of rel. conf-s after comparison of non-rel. confs

## [0.7.12] - 2021-09-26
- ine v1.6 - added new subroutine SetParams to set parameters for job and arrays
- added new key kIters=(0-iterate and invert if diverged, 2-invert only)

## [0.7.11] - 2021-09-26
- ine v1.5 - removed IP1 dependency from params.f90 and allows it to be set it in subroutine Input

## [0.7.10] - 2021-09-26
- ine v1.4 - revamped input for range of wavelengths to accept single wavelength by inputting (wavelength1 wavelength1 x) as range
- added several code documentation for ine

## [0.7.9] - 2021-09-25
- minor text fix: total computation time of conf_lsj 

## [0.7.8] - 2021-09-25
- conf_lsj v0.2.0 - major revamp and optimization of lsj routine
- relevant eigenvectors are all stored in memory in array B1h(Nd,nlvs)
- one-electron matrix elements of l, s and j are stored in arrays instead of calculated every iteration
- in the inner loop of lsj, conf-s are skipped if they belong to different non-rel. conf-s (using new subroutine CompNRC) 
- selection rules are added to skip zero contributions
- many-electron matrix elements are calculated only once and used for all levels

## [0.7.7] - 2021-09-25
- conf v0.3.14 - key Kw revamped to be an input parameter along with Kl, Ksig and Kdsig

## [0.7.6] - 2021-09-25
- minor text fix for ine
- moved subroutines from ine_aux.f90 to ine.f90

## [0.7.5] - 2021-09-25
- ine now reads a range of wavelengths along with step instead of a single wavelength

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
