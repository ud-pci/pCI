add_library(pci_serial params.f90 conf_variables.F90 basc_variables.f90 str_fmt.f90 readfff.f90
         breit.f90 wigner.f90 sintg.f90 bspl.f90 conf_init.f90 test_ori.f90 pi_pk.f90 integrals.f90 determinants.f90
         solvers.f90 diff.f90 utils.f90)

target_compile_definitions(pci_serial PUBLIC $<$<BOOL:${USE_DP_INTEGRALS}>:USE_DP_INTEGRALS>)

add_executable(hfd hfd.f90)
target_link_libraries(hfd PRIVATE pci_serial)

add_executable(bass bass.f90)
target_link_libraries(bass PRIVATE pci_serial)

add_executable(add add.f90)
target_link_libraries(add PRIVATE pci_serial)

add_executable(sort sort.f90)
target_link_libraries(sort PRIVATE pci_serial)

add_executable(check_matrices check_matrices.f90)
target_link_libraries(check_matrices PRIVATE pci_serial)

install(TARGETS hfd bass add sort check_matrices DESTINATION ${CMAKE_INSTALL_BINDIR})

if (LAPACK_AVAILABLE)
    add_executable(ine ine.f90)
    target_link_libraries(ine PRIVATE pci_serial mkl_parallel)

    add_executable(pol pol.f90)
    target_link_libraries(pol PRIVATE pci_serial mkl_parallel)

    install(TARGETS ine pol DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(MPI_FOUND)
    add_library(pci_common mpi_utils.f90)
    target_compile_definitions(pci_common PUBLIC ${MPI_Fortran_COMPILE_DEFINITIONS})
    target_compile_options(pci_common PUBLIC ${MPI_Fortran_COMPILE_OPTIONS})
    target_include_directories(pci_common PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
    target_link_options(pci_common PUBLIC ${MPI_Fortran_LINK_FLAGS})
    target_link_libraries(pci_common PUBLIC pci_serial ${MPI_Fortran_LIBRARIES})

    add_executable(pbasc pbasc.f90)
    target_link_libraries(pbasc PUBLIC pci_common)

    if (LAPACK_AVAILABLE)
        add_executable(pconf env_var.f90 matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 csf.f90 pconf.F90)
        target_link_libraries(pconf PUBLIC pci_common mkl_sequential)
        if (USE_SCALAPACK)
            target_compile_definitions(pconf PRIVATE USE_SCALAPACK)
            if (MKL_FOUND OR MKL_ROOT)
                target_link_options(pconf PRIVATE -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64)
            elseif (ScaLAPACK_FOUND)
                target_link_libraries(pconf PRIVATE ${ScaLAPACK_LIBRARIES})
            endif()
        endif()
        install(TARGETS pconf DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    add_executable(conf_lsj matrix_io.f90 vaccumulator.f90 formj2.f90 davidson.f90 conf_lsj.f90)
    target_link_libraries(conf_lsj PUBLIC pci_common)

    add_executable(conf_pt conf_pt.f90)
    target_link_libraries(conf_pt PUBLIC pci_common)

    add_executable(pdtm amp_ops.f90 pdtm.f90)
    target_link_libraries(pdtm PUBLIC pci_common)

    add_executable(formy formy.f90)
    target_link_libraries(formy PUBLIC pci_common)

    install(TARGETS pbasc conf_lsj pdtm conf_pt formy DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
