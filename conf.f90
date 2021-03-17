! ======== original version by I.I.Tupitsin =======
! ======== parallel version by C. Cheung ==========
! ### last updated on 7/13/20
!   /01/17: file CONF.HIJ is cut in pieces of fixed number of 
!           records (parameter IPnr in hread.par)
! 12/11/14: new par file conf8.par & dimension of P matrix
!           is now defined by parameter IPlv
! 25/08/13: version with integer*8 NumH
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
!   12 - 'CONF.DAT'   radial functions after
!                     ortogonalization
!   13 - 'CONF.INT'   radial integrals
!   15 - 'CONF.HIJ'   Hamiltonian matrix
!   17 - 'CONF.XIJ'   eigenvectors in basis of dets
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
Program conf 
    use mpi
    use conf_aux
    use davidson, only : Prj_J
    use determinants, only : Wdet, FormD, Dinit, Jterm
    use integrals, only : Rint, RintS
    use conf_init, only : Input, Init
    use formj2, only : FormJ, J_av
    Use str_fmt, Only : FormattedTime
    !use conffilepaths

    implicit none

    integer   :: mype, npes, mpierr
    Integer(kind=int64) :: clock_rate
    Integer(kind=int64)      :: start_time, stop_time, start_time_tot, stop_time_tot, start1, end1
    Real :: ttime
    INTEGER, ALLOCATABLE :: data(:,:) 
    integer   :: n, k, ierr, i, j, nerr
    real(dp)  :: t
    character(len=1024) :: strFromEnv
    character(len=255) :: eValue
    Character(Len=16)     :: memStr, timeStr
!   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Initialize MPI
    call MPI_Init(mpierr)
    ! Get process id
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    ! Get number of processes
    call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    call system_clock(count_rate=clock_rate)
    if (mype==0) call system_clock(start_time_tot)
    
    ! Give ConfFilePaths a chance to decide what filenames/paths to use:
    !Call ConfFileInit()

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    call Get_Environment_Variable("CONF_MAX_BYTES_PER_CPU",eValue)
    read(eValue,'(I12)') memTotalPerCPU

    Kw=0 ! if Kw=0, CONF.HIJ files are not written
         ! if Kw=1, CONF.HIJ files are written
    kXIJ=10     ! kXIJ sets the interval in which CONF.XIJ is written
               ! e.g. kXIJ=5 => CONF.XIJ written every 5 davidson iterations
               ! if kXIJ=0, then no intermediate CONF.XIJ will be written
    ! Only the master core needs to initialize the conf program
    if (mype == 0) then
        open(unit=11,status='UNKNOWN',file='CONF.RES')
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call Input ! reads list of configurations from CONF.INP
        call Init ! reads basis set information from CONF.DAT
        call Rint ! reads radial integrals from CONF.INT
        if (Ksig /= 0) call RintS 
        call Dinit      ! forms list of determinants
        call Jterm      ! prints table with numbers of levels with given J
        call Wdet('CONF.DET')       ! writes determinants to file CONF.DET
        if (Ksig*Kdsig /= 0) call FormD
        call system_clock(start_time)
    end if

    call AllocateFormHArrays(mype,npes)
    
    ! Evaluation of Hamiltonian
    if (Kl <= 2 .or. Kl4 == 0) call FormH(npes,mype)
    
    if (mype == 0) then
        call system_clock(stop_time)
        ttime=Real((stop_time-start_time)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> FormH took '// trim(timeStr) // ' to complete'
        call system_clock(start_time)
    end if

    call FormJ(mype, npes)   ! calculates matrix J^2 and writes it to CONF.JJJ
    
    if (mype == 0) then
        call system_clock(stop_time)
        ttime=Real((stop_time-start_time)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> FormJ took '// trim(timeStr) // ' to complete'
    end if

    call DeAllocateFormHArrays(mype,npes)
    call AllocateDvdsnArrays(mype,npes)

    if (mype==0) call system_clock(start_time)
    if (Kv == 3 .or. Kv == 4) call Diag4(mype, npes) ! Davidson diagonalisation
    deallocate(H_n, H_k, H_t)
    
    if (mype == 0) then
        call system_clock(stop_time)
        ttime=Real((stop_time-start_time)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Davidson procedure took '// trim(timeStr) // ' to complete'

        D1(1:Nlv)=Tk(1:Nlv)
        nerr=0
        if (K_prj == 1) then
            call Prj_J(1,Nlv,Nlv+1,ierr,1.d-8)
            if (ierr == 0) nerr=nerr+1
        end if
        open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
    end if
    do n=1,Nlv
        call MPI_Bcast(ArrB(1:Nd,n), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call J_av(ArrB(1,n),Nd,Tj(n),ierr,mype,npes)  ! calculates expectation values for J^2
        if (mype==0) write (16) D1(n),Tj(n),Nd,(ArrB(i,n),i=1,Nd)
    end do
    if (mype==0) then
        close(unit=16)
        call Print   !#   Output of the results
        close(unit=6)
        close(unit=11)
        call system_clock(stop_time_tot)
        ttime=Real((stop_time_tot-start_time_tot)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of conf was '// trim(timeStr)
    end if
    deallocate(ArrB, B1, B2, Diag, Tk, Ndc)
    call MPI_Finalize(mpierr)
End Program
