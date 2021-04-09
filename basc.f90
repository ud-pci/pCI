Program basc

    Use mpi
    Use basc_aux
    Use breit, Only : Gaunt
    Use str_fmt, Only : FormattedTime

    Implicit None

    Integer :: mype, npes, mpierr
    Integer :: Norb, nsx, nsx2, lsx
    Integer(kind=int64) :: clock_rate, start_time, end_time
    Real :: total_time
    Character(Len=16) :: timeStr

    MaxT=9           !### length of expansion at the origin
    fname='CONF.DAT'
    Kout=1

    ! Initialize MPI
    Call MPI_Init(mpierr)
    ! Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    ! Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    Call system_clock(count_rate=clock_rate)
    If (mype==0) Call system_clock(start_time)

    If (mype == 0) Then
        Open(unit=11,file='BASC.RES',status='UNKNOWN')
        Call recunit
        Call Input
        ! Norb = number of orbitals to be formed by Fbas:
        Call Init(Norb)
        if (Norb.NE.0) then
          write (*,*) ' Run "bass" to form ',Norb,' new orbitals.'
          stop
        end if
        Call Fill_N_l
        Call Core
        Call Rint0(nsx,nsx2,lsx)
        Call Rint(nsx,nsx2,lsx)
        write(*,*)' RINT finished'
        Call Gaunt
        write(*,*)' GAUNT finished'
    End If

    If (mype == 0) Then
        Close(11)
        Call system_clock(end_time)
        total_time=Real((end_time-start_time)/clock_rate)
        Call FormattedTime(total_time, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of basc was '// trim(timeStr)
    End If

    Call MPI_Finalize(mpierr)

End Program basc