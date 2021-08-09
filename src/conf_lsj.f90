Program conf_lsj
    use mpi
    use conf_aux
    use davidson, Only : Prj_J
    use determinants, Only : Wdet, FormD, Dinit, Jterm
    use integrals, Only : Rint, RintS
    use conf_init, Only : Init, InitFormH
    use formj2, Only : FormJ, J_av
    Use str_fmt, Only : FormattedTime
    Use lsj_aux
    !Use conffilepaths

    Implicit None

    Integer   :: n, k, i, j, ierr, mype, npes, mpierr, nnd
    Integer(kind=int64) :: clock_rate
    Integer(kind=int64) :: start_time, end_time
    Real :: total_time
    Real(dp)  :: t, xtj, xtl, xts, xj
    Character(Len=1024) :: strFromEnv
    Character(Len=255)  :: eValue, strfmt
    Character(Len=16)   :: memStr, timeStr
!   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Initialize MPI
    Call MPI_Init(mpierr)
    ! Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    ! Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    Call system_clock(count_rate=clock_rate)
    If (mype==0) Call system_clock(start_time)
    
    ! Give ConfFilePaths a chance to decide what filenames/paths to use:
    !Call ConfFileInit()

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    Call Get_Environment_Variable("CONF_MAX_BYTES_PER_CPU",eValue)
    read(eValue,'(I12)') memTotalPerCPU

    Kw=0 ! If Kw=0, CONF.HIJ files are not written
         ! If Kw=1, CONF.HIJ files are written
    kXIJ=10    ! kXIJ sets the interval in which CONF.XIJ is written
               ! e.g. kXIJ=5 => CONF.XIJ written every 5 davidson iterations
               ! If kXIJ=0, then no intermediate CONF.XIJ will be written
    ! Only the master core needs to initialize the conf program
    If (mype == 0) Then
        open(unit=11,status='UNKNOWN',file='CONF_LSJ.RES')
        strfmt = '(4X,"Program conf_lsj v0.1.0")'
        Write( 6,strfmt)
        Write(11,strfmt)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Call Input ! reads list of configurations from CONF.INP
        Call Init ! reads basis set information from CONF.DAT
        Call Rint ! reads radial integrals from CONF.INT
        Call Dinit      ! forms list of determinants
        Call Jterm      ! prints table with numbers of levels with given J
    End If

    Call AllocateFormHArrays(mype,npes)
    Call InitFormH(mype,npes) 
    Allocate(Tj(Nlv),Tl(Nlv),Ts(Nlv),B1(Nd),D1(Nlv))
    If (mype == 0) Open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
    
    Do n=1,Nlv
        If (mype == 0) Read(16) D1(n),Tj(n),nnd,(B1(i),i=1,Nd)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nnd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tj(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(D1(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(B1(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call lsj(B1,Nd,xtj,xtl,xts,mype,npes)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        xj=0.5d0*(sqrt(xtj+1)-1)
        Tl(n)=0.5d0*(sqrt(xtl+1)-1)
        Ts(n)=0.5d0*(sqrt(xts+1)-1)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    End Do

    ! Print table of final results and total computation time
    If (mype==0) Then
        Call PrintLSJ   !#   Output of the results

        Call system_clock(end_time)
        total_time=Real((end_time-start_time)/clock_rate)
        Call FormattedTime(total_time, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of conf was '// trim(timeStr)
    End If

    Call MPI_Finalize(mpierr)

End Program conf_lsj
