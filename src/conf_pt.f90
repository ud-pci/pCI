Program conf_pt   
    ! #############################################################################
    ! # This version of conf_pt forms corrections to eigenvectors in CONF_PT.XIJ           
    ! # and forms ordered CONF_new.INP
    ! #############################################################################
    ! This code forms the second order corrections to eigenvalues from CONF.XIJ. 
    ! It is supposed that vectors in CONF.XIJ correspond to the same Hamiltonian 
    ! in a smaller space.
    Use mpi
    Use conf_pt_aux
    Use integrals, Only : Rint
    Use determinants, Only : Dinit, Jterm, Wdet
    Implicit None
    Character(Len=8) cdiag(2)
    Integer  :: Nc_0, i, mpierr, mype, npes
    Real :: start_time, stop_time
    Data cdiag /'diag(H) ','diag(H0)'/
    !     - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ! |||         Maximal values for parameters            |||
    ! |||         which determine array dimensions:        |||
    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !                   Nsp  <= IPsp
    !                   Ns   <= IPs
    !                   Nc   <= Nc
    !                   Nd   <= IP4
    !                   Nlv  <= 3*Nlv
    !             Other dimensions are:
    !                   valence e              2*3*Nlv
    !                   positions for val. e   Nst
    !                   different Jtot         IPjd
    !     - - - - - - - - - - - - - - - - - - - - - - - - -
    !Init mpi
    Call MPI_Init(mpierr)
    !Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    !Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    If (mype == 0) Then
        ! CONF_PT.RES is opened only by the master process
        open(unit=11,status='UNKNOWN',file='CONF_PT.RES')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call Input
        Call Init
        Call Rint
        Call Dinit
        Call Jterm
        Call Wdet('CONF_PT.DET')
        Call NR_Init
        Call Weight_CI
    End If
    
    Call BcastParams(mype, npes)
    Call AllocatePTArrays(mype, npes)
    Call BcastPTArrays(mype, npes)

    Call cpu_time(start_time)
    Call PT_Init(npes, mype)
    Call cpu_time(stop_time)

    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    Call MPI_Bcast(En(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    If (mype == 0) Then
        print*, 'Total time for diagonal was', stop_time-start_time, 'seconds.'
    End If

    Call cpu_time(start_time)
    Call PTE(npes, mype)
    Call cpu_time(stop_time)
    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    If (mype == 0) Then
        print*, 'Total time for off-diagonal was', stop_time-start_time, 'seconds.'
        Call SaveVectors
        Call Weight_PT
        Call Cutoff(Nc_0)
        Call Sort(Nc_0)
        Call NewConfList(Nc_0)
        ! - - - - - - -  - - - - - - -  - - - - - 
        Close(11)
    End If

    Call MPI_Finalize(mpierr)
End Program conf_pt
! =============================================================================