program dtm 
    use mpi_f08
    use dtm_aux
    Use determinants, Only : Dinit, Jterm
    implicit none
    integer :: nsu2, mpierr, mype, npes
    real :: start_time_total, end_time_total

    call MPI_Init(mpierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    if (mype==0) then
      call cpu_time(start_time_total) 
      call Input  ! reads file CONF.INP
      call Init   ! reads file CONF.DAT
      call RintA   ! either reads or calculates radial integrals
      call Dinit   ! forms list of determinants (CONF.DET)
      call Jterm   ! produces table of levels and generates basis set of determinants from CONF.INP
    end if
    call InitTDM(mype,npes)
    select case(Kl1)
      case(1) 
        call FormDM(mype,npes)  ! calculates density matrix and expectation values
      case(2)
        call FormTM(mype,npes)  ! calculates transition matrix & amplitudes
    end select
  
    if (mype==0) then  
      write( *,'(4X,"Threshold: ",E8.1,";",2X,A4," form of M1 operator ")') Trd,chm1(K_M1)
      write(11,'(4X,"Threshold: ",E8.1,";",2X,A4," form of M1 operator ")') Trd,chm1(K_M1)
  
      if (Gj /= 0.d0) then
        write( *,'(4X,"Note that Gj=",F8.5," is assumed for ALL levels")') Gj
        write(11,'(4X,"Note that Gj=",F8.5," is assumed for ALL levels")') Gj
      end if
      close(11)
      call cpu_time(end_time_total)
      print*, 'dtm took', (end_time_total-start_time_total)/60 , 'minutes.'
    end if
    call MPI_Finalize(mpierr)
end program dtm
