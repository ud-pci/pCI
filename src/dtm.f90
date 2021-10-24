Program dtm 
    Use mpi_f08
    Use dtm_aux
    Use determinants, Only : Dinit, Jterm
    Implicit None
    
    Integer :: mpierr, mype, npes
    Real    :: start_time_total, end_time_total

    Call MPI_Init(mpierr)
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    If (mype==0) Then
        Call cpu_time(start_time_total) 
        Call Input   ! reads file CONF.INP
        Call Init    ! reads file CONF.DAT
        Call RintA   ! either reads or calculates radial integrals
        Call Dinit   ! forms list of determinants (CONF.DET)
        Call Jterm   ! produces table of levels and generates basis set of determinants from CONF.INP
    End If

    Call InitTDM

    Select Case(Kl1)
        Case(1) 
            Call FormDM(mype,npes)  ! calculates density matrix and expectation values
        Case(2)
            Call FormTM(mype,npes)  ! calculates transition matrix & amplitudes
    End Select
  
    If (mype==0) Then  
        Write( *,'(4X,"Threshold: ",E8.1,";",2X,A4," form of M1 operator ")') Trd,chm1(K_M1)
        Write(11,'(4X,"Threshold: ",E8.1,";",2X,A4," form of M1 operator ")') Trd,chm1(K_M1)
    
        If (Gj /= 0.d0) Then
            Write( *,'(4X,"Note that Gj=",F8.5," is assumed for ALL levels")') Gj
            Write(11,'(4X,"Note that Gj=",F8.5," is assumed for ALL levels")') Gj
        End If
        Close(11)
        Call cpu_time(end_time_total)
        Print*, 'dtm took', (end_time_total-start_time_total)/60 , 'minutes.'
    End If
    Call MPI_Finalize(mpierr)

End Program dtm
