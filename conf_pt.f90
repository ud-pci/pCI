program conf_pt   
    ! #############################################################################
    ! # This version of conf_pt forms corrections to eigenvectors in CONF_PT.XIJ           
    ! # and forms ordered CONF_new.INP
    ! #############################################################################
    ! This code forms the second order corrections to eigenvalues from CONF.XIJ. 
    ! It is supposed that vectors in CONF.XIJ correspond to the same Hamiltonian 
    ! in a smaller space.
    use mpi
    use conf_pt_aux
    use integrals, only : Rint
    use determinants, only : Dinit, Jterm, Wdet
    implicit none
    character(len=8) cdiag(2)
    integer  :: Nc_0, i
    integer  :: mpierr, mype, npes
    real :: start_time, stop_time
    data cdiag /'diag(H) ','diag(H0)'/
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
    call MPI_Init(mpierr)
    !Get process id
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    !Get number of processes
    call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)
    !This is only written by the first process

    if (mype == 0) then
        !This is written by every process
        open(unit=11,status='UNKNOWN',file='CONF_PT.RES')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        call Input
        call Init
        call Rint
        call Dinit
        call Jterm
        call Wdet('CONF_PT.DET')
        call NR_Init
        call Weight_CI
    end if
    
    call BcastParams(mype, npes)
    call AllocatePTArrays(mype, npes)
    call BcastPTArrays(mype, npes)

    call cpu_time(start_time)
    call PT_Init(npes, mype)
    call cpu_time(stop_time)

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(En(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (mype == 0) then
        print*, 'Total time for diagonal was', stop_time-start_time, 'seconds.'
    end if

    call cpu_time(start_time)
    call PTE(npes, mype)
    call cpu_time(stop_time)
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (mype == 0) then
        print*, 'Total time for off-diagonal was', stop_time-start_time, 'seconds.'
        call SaveVectors
        call Weight_PT
        call Cutoff(Nc_0)
        call Sort(Nc_0)
        call NewConfList(Nc_0)
        ! - - - - - - -  - - - - - - -  - - - - - 
        close(11)
    end if

    call MPI_Finalize(mpierr)
end program conf_pt
! =============================================================================