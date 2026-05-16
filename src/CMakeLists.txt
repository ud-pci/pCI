ADD_LIBRARY(pci_common params.f90 conf_variables.F90 basc_variables.f90 str_fmt.f90 readfff.f90 
         breit.f90 wigner.f90 sintg.f90 bspl.f90 conf_init.f90 test_ori.f90 pi_pk.f90 integrals.f90 determinants.f90 
         mpi_utils.f90 solvers.f90 diff.f90 utils.f90)

target_compile_definitions(pci_common PUBLIC ${MPI_Fortran_COMPILE_DEFINITIONS})
target_compile_options(pci_common PUBLIC ${MPI_Fortran_COMPILE_OPTIONS})
target_include_directories(pci_common PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
target_link_options(pci_common PUBLIC ${MPI_Fortran_LINK_FLAGS})
target_link_libraries(pci_common PUBLIC ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(hfd hfd.f90)
target_link_libraries(hfd PUBLIC pci_common)

ADD_EXECUTABLE(bass bass.f90)
target_link_libraries(bass PUBLIC pci_common)

ADD_EXECUTABLE(add add.f90)
target_link_libraries(add PUBLIC pci_common)

ADD_EXECUTABLE(pbasc pbasc.f90)
target_link_libraries(pbasc PUBLIC pci_common)

ADD_EXECUTABLE(pconf env_var.f90 matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 csf.f90 pconf.F90)
target_link_libraries(pconf PUBLIC pci_common)

ADD_EXECUTABLE(conf_lsj matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 conf_lsj.f90)
target_link_libraries(conf_lsj PUBLIC pci_common)

ADD_EXECUTABLE(conf_pt conf_pt.f90)
target_link_libraries(conf_pt PUBLIC pci_common)

ADD_EXECUTABLE(pdtm amp_ops.f90 pdtm.f90)
target_link_libraries(pdtm PUBLIC pci_common)

ADD_EXECUTABLE(sort sort.f90)
target_link_libraries(sort PUBLIC pci_common)

ADD_EXECUTABLE(check_matrices check_matrices.f90)
target_link_libraries(check_matrices PUBLIC pci_common)

ADD_EXECUTABLE(ine ine.f90)
target_link_options(ine PRIVATE ${MKL_PARALLEL_LINK_FLAGS})
target_link_libraries(ine PRIVATE ${MKL_PARALLEL_LINK_FLAGS} pci_common)

ADD_EXECUTABLE(pol pol.f90)
target_link_options(pol PRIVATE ${MKL_PARALLEL_LINK_FLAGS})
target_link_libraries(pol PRIVATE ${MKL_PARALLEL_LINK_FLAGS} pci_common)

ADD_EXECUTABLE(formy formy.f90)
target_link_libraries(formy PUBLIC pci_common)

INSTALL(TARGETS hfd bass add pbasc pconf conf_lsj pdtm ine pol conf_pt sort check_matrices formy DESTINATION ${CMAKE_INSTALL_BINDIR})
