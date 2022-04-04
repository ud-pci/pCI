#
#
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
CMAKE_POLICY(SET CMP0074 NEW)

PROJECT(conf_new LANGUAGES Fortran)

INCLUDE(GNUInstallDirs)

SET(CMAKE_Fortran_FLAGS "-mkl -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64")

ADD_LIBRARY(basc_common params.f90 str_fmt.f90 conf_variables.f90 basc_variables.f90 conf_init.f90 readfff.f90 wigner.f90 sintg.f90 test_ori.f90 diff.f90 pi_pk.f90 breit.f90)

ADD_LIBRARY(conf_common params.f90 str_fmt.f90 conf_variables.f90 determinants.f90 integrals.f90 formj2.f90 matrix_io.f90 davidson.f90 vaccumulator.f90 mpi_utils.f90)

ADD_LIBRARY(lsj_common params.f90 str_fmt.f90 conf_variables.f90 determinants.f90 integrals.f90 formj2.f90 matrix_io.f90 davidson.f90 vaccumulator.f90)

ADD_LIBRARY(dtm_common params.f90 str_fmt.f90 conf_variables.f90 dtm_variables.f90 conf_init.f90 wigner.f90 determinants.f90 amp_ops.f90)

ADD_LIBRARY(ine_common params.f90 str_fmt.f90 ine_variables.f90 conf_variables.f90 determinants.f90 integrals.f90 conf_init.f90 wigner.f90 solvers.f90)

ADD_LIBRARY(conf_pt_common params.f90 conf_variables.f90 conf_pt_variables.f90 conf_init.f90 integrals.f90 determinants.f90)

ADD_EXECUTABLE(basc basc.f90)
TARGET_LINK_LIBRARIES(basc basc_common) 

ADD_EXECUTABLE(conf conf_init.f90 conf.f90)
TARGET_LINK_LIBRARIES(conf conf_common)

ADD_EXECUTABLE(conf_lsj conf_init.f90 conf_lsj.f90)
TARGET_LINK_LIBRARIES(conf_lsj lsj_common)

ADD_EXECUTABLE(dtm dtm_aux.f90 dtm.f90)
TARGET_LINK_LIBRARIES(dtm dtm_common)

ADD_EXECUTABLE(ine ine.f90)
TARGET_LINK_LIBRARIES(ine ine_common)

ADD_EXECUTABLE(conf_pt conf_pt.f90)
TARGET_LINK_LIBRARIES(conf_pt conf_pt_common)

ADD_EXECUTABLE(add params.f90 add.f90)

INSTALL(TARGETS add basc conf conf_lsj dtm ine conf_pt DESTINATION ${CMAKE_INSTALL_BINDIR})