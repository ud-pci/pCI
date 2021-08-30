Program conf 
    ! ======== original version by I.I.Tupitsin =======
    ! ======== parallel version by C. Cheung ==========
    ! latest version of parallel code can be found at
    ! >>>   https://github.com/ccheung93/pCI.git  <<<
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    ! this code needs radial integrals in CONF.INT, which 
    ! are calculated by basc.
    ! the file CONF.DAT is still needed for INIT subroutine.
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   Kl - key:
    !            Kl=0 - everything is calculated.
    !            Kl=1 - information from CONF.* files is used.
    !            Kl=2 - same as Kl=0, but with Sigma and Screening.
    !            Kl=3 - same as Kl=1, but status of CONF.HIJ is
    !                   changed. That allows to add configurations
    !                   to the previous calculation.
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   Kv - key for the variant:
    !              Kv=1 - direct diagonalization and J selection
    !              Kv=2 - direct diagonalization
    !              Kv=3 - Davidson diagonalization and J selection
    !              Kv=4 - Davidson diagonalization
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   INPUT/OUTPUT channels:
    !   10 - 'CONF.INP'   input file.
    !   11 - 'CONF.RES'   text file with results.
    !   12 - 'CONF.DAT'   radial functions after orthogonalization
    !   13 - 'CONF.INT'   radial integrals
    !   15 - 'CONF.HIJ'   Hamiltonian matrix
    !   17 - 'CONF.XIJ'   eigenvectors in basis of determinants
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main variables:
    !   Ns      - number of orbitals (different)
    !   Nso     - number of core orbitals
    !   Nsp     - total number of shells (>=Ns, can be equal).
    !   Qnl(I)  - atomic configurations (I=1,Nsp)
    !   Jm      - projection of total momentum J
    !   Nc      - number of configurations
    !   Nd      - number of dets
    !   Ne      - number of valence electrons
    !   Nec     - number of core electrons
    !   Ecore   - core energy
    !   IPmr    - equivalent of 4 Bytes for the DIRECT files
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    use mpi
    use conf_aux
    use davidson, only : Prj_J
    use determinants, only : Wdet, FormD, Dinit, Jterm
    use integrals, only : Rint, RintS
    use conf_init, only : Init, InitFormH
    use formj2, only : FormJ, J_av
    Use str_fmt, Only : FormattedTime
    !Use conffilepaths

    Implicit None

    Integer   :: n, k, i, j, ierr, mype, npes, mpierr
    Integer(kind=int64) :: clock_rate
    Integer(kind=int64) :: start_time, end_time
    Real :: total_time
    Real(dp)  :: t
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

    Kw=1 ! If Kw=0, CONF.HIJ files are not written
         ! If Kw=1, CONF.HIJ files are written
    kXIJ=10    ! kXIJ sets the interval in which CONF.XIJ is written
               ! e.g. kXIJ=5 => CONF.XIJ written every 5 davidson iterations
               ! If kXIJ=0, then no intermediate CONF.XIJ will be written
    ! Only the master core needs to initialize the conf program
    If (mype == 0) Then
        open(unit=11,status='UNKNOWN',file='CONF.RES')
        strfmt = '(4X,"Program conf v0.3.8")'
        Write( 6,strfmt)
        Write(11,strfmt)
        Call Input ! reads list of configurations from CONF.INP
        Call Init ! reads basis set information from CONF.DAT
        Call Rint ! reads radial integrals from CONF.INT
        If (Ksig /= 0) Call RintS ! reads in radial integrals from SGC.CON and SCRC.CON
        Call Dinit      ! forms list of determinants
        Call Jterm      ! prints table with numbers of levels with given J
        Call Wdet('CONF.DET')       ! writes determinants to file CONF.DET
        If (Ksig*Kdsig /= 0) Call FormD
    End If

    Call AllocateFormHArrays(mype,npes)
    
    ! Evaluation of Hamiltonian
    !Call InitFormH(npes,mype)
    If (Kl <= 2 .or. Kl4 == 0) Call FormH(npes,mype)

    Call FormJ(mype, npes)   ! calculates matrix J^2 and writes it to CONF.JJJ
    
    Call DeAllocateFormHArrays(mype,npes)
    Call AllocateDvdsnArrays(mype,npes)

    Call Diag4(mype,npes) ! Davidson diagonalization
    
    Call WriteFinalXIJ(mype,npes)

    ! Print table of final results and total computation time
    If (mype==0) Then
        Call PrintEnergies   !#   Output of the results

        Call system_clock(end_time)
        total_time=Real((end_time-start_time)/clock_rate)
        Call FormattedTime(total_time, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of conf was '// trim(timeStr)
    End If

    Call MPI_Finalize(mpierr)

End Program conf
