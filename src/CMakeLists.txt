SET(CMAKE_Fortran_FLAGS "-qopenmp -mkl=sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64")

ADD_LIBRARY(pci_common params.f90 conf_variables.f90 basc_variables.f90 str_fmt.f90 readfff.f90 
        breit.f90 wigner.f90 sintg.f90 bspl.f90 conf_init.f90 test_ori.f90 pi_pk.f90 integrals.f90 determinants.f90 
        mpi_utils.f90 solvers.f90 diff.f90 utils.f90)

ADD_EXECUTABLE(hfd hfd.f90)
TARGET_LINK_LIBRARIES(hfd pci_common) 

ADD_EXECUTABLE(bass bass.f90)
TARGET_LINK_LIBRARIES(bass pci_common) 

ADD_EXECUTABLE(pbasc pbasc.f90)
TARGET_LINK_LIBRARIES(pbasc pci_common) 

ADD_EXECUTABLE(pconf env_var.f90 matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 pconf.f90)
TARGET_LINK_LIBRARIES(pconf pci_common)

ADD_EXECUTABLE(conf_lsj matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 conf_lsj.f90)
TARGET_LINK_LIBRARIES(conf_lsj pci_common)

ADD_EXECUTABLE(pdtm dtm_variables.f90 amp_ops.f90 pdtm.f90)
TARGET_LINK_LIBRARIES(pdtm pci_common)

ADD_EXECUTABLE(ine ine.f90)
TARGET_LINK_LIBRARIES(ine pci_common)

ADD_EXECUTABLE(pol pol.f90)
TARGET_LINK_LIBRARIES(pol pci_common)

ADD_EXECUTABLE(conf_pt conf_pt.f90)
TARGET_LINK_LIBRARIES(conf_pt pci_common)

ADD_EXECUTABLE(sort str_fmt.f90 sort.f90)

ADD_EXECUTABLE(add params.f90 add.f90)

INSTALL(TARGETS hfd bass add pbasc pconf conf_lsj pdtm ine pol conf_pt sort DESTINATION ${CMAKE_INSTALL_BINDIR})
