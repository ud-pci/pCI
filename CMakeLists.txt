#
#
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
CMAKE_POLICY(SET CMP0074 NEW)

PROJECT(conf_new LANGUAGES Fortran)

INCLUDE(GNUInstallDirs)

ADD_LIBRARY(conf_common params.f90 str_fmt.f90 conf_variables.f90 determinants.f90 integrals.f90 formj2.f90 hamiltonian_io.f90 davidson.f90 vaccumulator.f90 progress.f90)
ADD_LIBRARY(dtm_common params.f90 str_fmt.f90 dtm_variables.f90 wigner.f90)
ADD_LIBRARY(conf_pt_common params.f90 conf_variables.f90 conf_pt_variables.f90 conf_init.f90 integrals.f90 determinants.f90)

ADD_EXECUTABLE(conf conf_aux.f90 conf_init.f90 conf.f90)
TARGET_LINK_LIBRARIES(conf conf_common)

ADD_EXECUTABLE(dtm dtm_aux.f90 dtm.f90)
TARGET_LINK_LIBRARIES(dtm dtm_common)

ADD_EXECUTABLE(conf_pt conf_pt_aux.f90 conf_pt.f90)
TARGET_LINK_LIBRARIES(conf_pt conf_pt_common)

INSTALL(TARGETS conf dtm conf_pt DESTINATION ${CMAKE_INSTALL_BINDIR})