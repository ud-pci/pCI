Program pconf 
    ! ======== original version by I.I.Tupitsin =======
    ! ======== parallel version by C. Cheung ==========
    ! latest version of parallel code can be found at
    ! >>>   https://github.com/ud-pci/pCI   <<<
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    ! this code needs radial integrals in CONF.INT, which are calculated by basc.
    ! the file CONF.DAT is still needed for Init subroutine.
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
    !   Ns      - number of basis orbitals
    !   Nso     - number of core orbitals
    !   Nsp     - total number of subshells with occupations for each configuration
    !   Nsu     - number of used orbitals
    !   Nst     - number of used orbital positions
    !   Qnl(i)  - atomic configurations (i=1,Nsp)
    !   Jm      - projection of total momentum J
    !   Nc      - number of relativistic configurations
    !   Nd      - number of determinants
    !   Ne      - number of valence electrons
    !   Nec     - number of core electrons
    !   Ecore   - core energy
    !   IPmr    - equivalent of 4 Bytes for the DIRECT files
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main arrays:
    !   Nn(i)   - principal quantum number n (i=1,Ns)
    !   Ll(i)   - orbital angular momentum quantum number l (i=1,Ns)
    !   Kk(i)   - relativistic quantum number k (i=1,Ns)
    !   Jj(i)   - total angular momentum quantum number j (i=1,Ns)
    !   Iarr    - basis set of determinants
    !   ArrB    - eigenvectors
    !   Nh(i)   - orbital index for total magnetic quantum number (i=1,Ns)
    !   Jz(i)   - total magnetic quantum number of orbital index Nh(i) (i=1,Ns)
    !   Nq(i)   - number of electron occupancy for atomic orbital (i,Nsp)
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    Use conf_variables
    Use mpi_f08
    use determinants, only : Wdet, Dinit, Jterm, Rdet
    Use integrals, only : Rint
    use formj2, only : FormJ
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime
    Use env_var

    Implicit None
    External :: BLACS_PINFO, BLACS_EXIT
    Integer   :: n, nnd, i, ierr, mype, npes, mpierr, threadLevel
    Integer :: ncsf, nccj, max_ndcs, nbas
    Integer(kind=int64) :: start_time
    Character(Len=255)  :: eValue, strfmt
    Character(Len=16)   :: timeStr
    Real(dp), Allocatable, Dimension(:) :: xj, xl, xs, ax_array

    ! Initialize MPI
    Call MPI_Init_thread(MPI_THREAD_SINGLE, threadLevel, mpierr)
    ! Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    ! Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)
    Call BLACS_PINFO(mype,npes)
    
    Call startTimer(start_time)
    
    ! Set MPI type for type_real
    Select Case(type_real)
    Case(sp)
        mpi_type_real = MPI_REAL
    Case(dp)
        mpi_type_real = MPI_DOUBLE_PRECISION
    End Select

    Select Case(type2_real)
    Case(sp)
        mpi_type2_real = MPI_REAL
    Case(dp)
        mpi_type2_real = MPI_DOUBLE_PRECISION
    End Select

    ! Give ConfFilePaths a chance to decide what filenames/paths to use:
    !Call ConfFileInit()

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    If (mype == 0 .and. .not. GetEnvIsPresent("CONF_MAX_BYTES_PER_CPU")) Then
        print*, "WARNING: CONF_MAX_BYTES_PER_CPU not set in environment"
    End If
    memTotalPerCPU = GetEnvInteger64("CONF_MAX_BYTES_PER_CPU", 0_int64)

    ! Initialization subroutines (core 0 will read input files and broadcast parameters and arrays to all others)
    If (mype == 0) Then
        ! Open the generic logging file:
        open(unit=11,status='UNKNOWN',file='CONF.RES')

        Call Input                  ! reads list of configurations from CONF.INP
        Call Init                   ! reads basis set information from CONF.DAT
        Call Rint                   ! reads radial integrals from CONF.INT
        Call RintS                  ! reads in radial integrals from SGC.CON and SCRC.CON
        Call Dinit                  ! forms list of determinants
        Call Jterm                  ! prints table with numbers of levels with given J
        Call Wdet('CONF.DET')       ! writes determinants to file CONF.DET
        Call FormD
    End If

    ! Allocate global arrays used in formation of Hamiltonian
    Call AllocateFormHArrays(mype)
    
    ! Formation of Hamiltonian matrix
    Call FormH(npes,mype)

    ! Formation of J^2 matrix
    Call FormJ(mype, npes)   
    
    ! Deallocate arrays that are no longer used in FormH
    Call DeAllocateFormHArrays(mype)

    ! Allocate arrays that are used in Davidson procedure
    Call AllocateDvdsnArrays(mype)
    
    ! Davidson diagonalization
    Call Diag4(mype,npes) 
    
    ! Print table of final energies
    If (mype==0) Call PrintEnergies

    ! Call LSJ routines if KLSJ=1
    If (KLSJ == 1) Then
        If (mype == 0) Call Rdet('CONF.DET')
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call AllocateLSJArrays(mype)
        Call InitLSJ(mype)
        Call lsj(ArrB,xj,xl,xs,mype,npes)
    End If
    
    ! Print table of final results and total computation time
    If (mype==0) Then
        Call PrintWeights
        
        Call stopTimer(start_time, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of conf was '// trim(timeStr)
    End If

    ! In rank 0, close-out unit 11
    If (mype == 0) Then
        Close(11)
    End If

    ! All done: 
    Call BLACS_EXIT(0)

Contains

    Subroutine Input
        ! This subroutine reads in parameters and list of configurations from CONF.INP and job configuration from ci.in
        Use conf_init, only : ReadConfInp, ReadConfigurations, ReadCiIn
        Implicit None
        Integer :: err_stat
        Character(Len=64) :: strfmt

        Select Case(type_real)
        Case(sp)
            strfmt = '(4X,"Program pconf v6.1 with single precision")'
        Case(dp)            
            Select Case(type2_real)
            Case(sp)
                strfmt = '(4X,"Program pconf v6.1")'
            Case(dp)
                strfmt = '(4X,"Program pconf v6.1 with double precision for 2e integrals")'
            End Select
        End Select
        
        Write( 6,strfmt)
        Write(11,strfmt)

        ! Read input parameters from file CONF.INP
        Call ReadConfInp

        If (dabs(C_is) < 1.d-6) K_is=0
        If (K_is == 0) C_is=0.d0

        ! Read job parameters from file ci.in
        Call ReadCiIn
        
        If (K_is == 2.OR.K_is == 4) Then
            Write( 6,'(/4X,"K_sms = (1-1el,2-2el,3-both) ",I1)') K_sms
            If ((K_sms-1)*(K_sms-2)*(K_sms-3) /= 0) Stop
        End If

        If (Kl == 2) Then
            Write(*,'(1X," Ksig = (0,1,2): ",I1)') Ksig 
            If (Ksig /= 0) Then
                Write(*,'(1X," Energy dependence of Sigma (1-Yes,0-No)? ",I1)') Kdsig
            End If
            Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
            Kecp=0
            Kl=0
        Else If (Kl == 3) Then
            Write(*,'(1X," Ksig = (0,1,2): ",I1)') Ksig 
            If (Ksig /= 0) Then
                Write(*,'(1X," Energy dependence of Sigma (1-Yes,0-No)? ",I1)') Kdsig
            End If
            Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
            Kecp=0
        Else
            Ksig=0
        End If

        ! Read configurations from file CONF.INP
        Call ReadConfigurations

        ! Set key for projection of J
        ! Kv = 1, 2 - Direct diagonaliztion (not implemented)
        ! Kv = 3, 4 - Davidson diagonalization
        ! If Kv = 3, then K_prj = 1 - use selection of states with projection of J
        ! If Kv = 4, then K_prj = 0 - no selection of states with projection of J
        If (Kv == 3) Then
            K_prj=1
            Write( *,'(4X,"Selection of states with J =",F5.1)') XJ_av
            Write(11,'(4X,"Selection of states with J =",F5.1)') XJ_av
        End If

        ! Read angular factors from file CONF.GNT
        Open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        Read(16,iostat=err_stat) Nlx, Ngaunt
        If (err_stat /= 0) Then
            Rewind(16)
            Read(16) Ngaunt
            print*, 'Nlx and Ngaunt could not be obtained from CONF.GNT, so defaulted to Nlx=5, Ngnt=', Ngaunt
            print*, 'INFO: Make sure basc and conf were compiled together in latest version!'
            Nlx = 5
            !Ngaunt = 2891
            !Rewind(16)
        End If

        Allocate(num_gaunts_per_partial_wave(Nlx+1))
        Allocate(In(Ngaunt))
        Allocate(Gnt(Ngaunt))

        If (err_stat /= 0) Then
            print*, 'num_gaunts_per_partial_wave could not be read and skipped'
        Else
            Read(16) num_gaunts_per_partial_wave
        End If
        Read(16) (In(i),i=1,Ngaunt)
        Read(16) (Gnt(i),i=1,Ngaunt)
        Close(unit=16)

    End Subroutine Input

    Subroutine Init
        Implicit None
        Integer  :: ic, n, j, imax, ni, ni2, kkj, llj, nnj, i, nj, nr, if, &
                    ii, i1, n2, n1, l, nmin, jlj, i0, nlmax, err_stat
        Real(dp) :: d, c1, c2, z1
        Real(dp), Dimension(IP6)  :: p, q, p1, q1 
        Real(dp), Dimension(4*IP6):: pq
        Integer, Dimension(0:33)  ::  nnn ,jjj ,nqq, nnn1, qqq1, nnn2, qqq2
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(0:33):: lll, lll1, lll2
        logical :: longbasis
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Character(Len=256) :: strfmt, err_msg

        Equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        Data Let/'s','p','d','f','g','h','i','k','l'/

        c1 = 0.01d0
        mj = 2*dabs(Jm)+0.01d0

        Open(12,file='CONF.DAT',status='OLD',access='DIRECT',recl=2*IP6*IPmr,iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.DAT is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If

        Read(12,rec=1) p
        Read(12,rec=2) q
        Read(12,rec=5) p1
        Read(12,rec=6) q1

        z1 = pq(1)
        If (dabs(Z-z1) > 1.d-6) Then
            strfmt = '("nuc. charge is changed: Z =",F12.6," ><",F12.6)'
            Write( 6,strfmt) Z,z1
            Write(11,strfmt) Z,z1
            Read(*,*)
        End If

        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Allocate(Nrnrc(Nc))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        Rnuc=pq(13)
        dR_N=pq(16)
        
        strfmt = '(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2,/4X, &
                    "Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3,5X,"Nc =",I7)'
        Write( 6,strfmt) Kl,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,strfmt) Kl,Z,Jm,Nsp,Ns,Nso,Nc

        longbasis=dabs(PQ(20)-0.98765d0) < 1.d-6
        If (longbasis) Then
            Write( *,*) ' Using variant for long basis '
            Write(11,*) ' Using variant for long basis '
            Do ni=1,Ns
                Nn(ni)=IQN(4*ni-3)
                Ll(ni)=IQN(4*ni-2)
                Kk(ni)=IQN(4*ni-1)
                Jj(ni)=IQN(4*ni)
            End Do
        Else
            if=20
            Do ni=1,Ns
                if=if+1
                Nn(ni)=pq(if)+c1
                if=if+1
                Ll(ni)=pq(if)+c1
                if=if+3
                c2=dsign(c1,pq(if))
                Kk(ni)=pq(if)+c2
                if=if+1
                c2=dsign(c1,pq(if))
                Jj(ni)=pq(if)+c2
            End Do
        End If

        Nsu=0
        Do nj=1,Nsp
            i=sign(1.d0,Qnl(nj))
            d=dabs(Qnl(nj))+1.d-14
            d=10.0*d
            nnj=d
            d=10.0d0*(d-nnj)
            llj=d
            jlj=2*llj+i
            kkj=-i*((jlj+1)/2)
            d=100.0d0*(d-llj)
            Nq(nj)=d+0.1d0
            Do ni=1,ns
                If (nnj == Nn(ni) .and. Kk(ni) == kkj) Then
                    Exit
                Else If (ni == ns) Then
                    strfmt = '(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)'
                    Write( 6,strfmt) nj,nnj,llj,kkj
                    Write(11,strfmt) nj,nnj,llj,kkj
                    Stop
                End If
            End Do
            Nip(nj)=ni
            If (Nsu < ni) Nsu=ni
        End Do

        Deallocate(Qnl)

        ! Calculate Nec - number of core electrons
        nec=0
        If (Nso /= 0) nec=sum(Nq(1:Nso))

        ! Calculate Nst - number of used orbital positions
        Do ni=1,Nsu
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                Nst=Nst+1
            End Do
        End Do
        strfmt = '(4X,"Number of actually used orbitals: Nsu =",I3, &
                    /4X,"Ne  =",I3,7X,"Nec =",I3,7X,"Nst =",I7)'
        Write( 6,strfmt) Nsu,Ne,nec,Nst
        Write(11,strfmt) Nsu,Ne,nec,Nst

        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        Do ni=nmin,Nsp
            i=i+1
            n=n+Nq(ni)
            If (n >= Ne) Then
                ic=ic+1
                If (n > Ne) Then
                    strfmt = '(/2X,"wrong number of electrons for configuration ICONF =",I6)'
                    Write( 6,strfmt) ic
                    Write(11,strfmt) ic
                  Stop
                End If
                Nvc(ic)=i
                Nc0(ic)=Nso+i0
                i0=i0+i
                n=0
                i=0
            End If
        End Do

        strfmt = '(1X,71("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        Do ni=1,Nso
            l =Ll(ni)+1
            lll(ni)=let(l)
        End Do

        strfmt = '(1X,"Core:", 6(I2,A1,"(",I1,"/2)",I2,";"), &
                           /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                           /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                           /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                           /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                           /6X,6(I2,A1,"(",I1,"/2)",I2,";"))'
        Write(11,strfmt) (Nn(i),lll(i),Jj(i),Nq(i),i=1,Nso)
        Write(11,'(1X,71("="))')

        ! Print list of relativistic configurations in format:
        ! nl(j) q; note that h (11/2) is currently printed as h (*/2) 
        ! and will result in conversion errors when debugging
        Nnr=1
        nr=0
        Do ic=1,Nc
            n1=Nc0(ic)+1
            n2=Nc0(ic)+Nvc(ic)
            nr=nr+1
            Do i=n1,n2
                i1=i-n1+1
                ni=Nip(i)
                l=Ll(ni)+1
                lll(i1)=let(l)
                jjj(i1)=Jj(ni)
                nnn(i1)=Nn(ni)
                nqq(i1)=Nq(i)
                If (Nq(i) > jjj(i1)+1) Then
                    strfmt = '(/2X,"wrong number of electrons for the shell:", &
                                I3,3X,I2,A1,I2,"/2"," (",I2,")")'
                    Write(11,strfmt) ni,nnn(i1),lll(i1),jjj(i1),nqq(i1)
                    Stop
                End If
            End Do
            n=n2-n1+1

            ! Convert relativistic configuration into non-relativistic configuration
            ni=1
            Do i=1,n
                If (nnn(i) == nnn(i-1) .and. lll(i) == lll(i-1)) Then
                    qqq1(ni-1) = qqq1(ni-1) + nqq(i)
                Else
                    nnn1(ni) = nnn(i)
                    lll1(ni) = lll(i)
                    qqq1(ni) = nqq(i)
                    ni=ni+1
                End If
            End Do

            If (ic > 1) Then
                If (ni2 /= ni) Then
                    Nnr=Nnr+1
                Else
                    Do i=1,ni
                        If (nnn2(i) == nnn1(i) .and. lll2(i) == lll1(i) .and. qqq2(i) == qqq1(i)) Then
                            Continue
                        Else
                            Nnr=Nnr+1
                            Exit
                        End If
                    End Do
                End If
            End If
            Nrnrc(Nnr) = nr
            ni2=ni
            nnn2=nnn1
            lll2=lll1
            qqq2=qqq1
            
            strfmt = '(1X,I6,"#",6(I2,A1,"(",I2,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I2,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I2,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I2,"/2)",I2,";"))'
            Write(11,strfmt) ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
        End Do

        Do i=Nnr,2,-1
            Nrnrc(i)=Nrnrc(i)-Nrnrc(i-1)
        End Do

        Write(11,'(1X,71("="))')

        If (Ksig > 0) Then
            Do ni=Nso+1,Nsu
               Read(12,rec=2*ni+7) p
               Eps(ni)=-p(ii+1)
            End Do
            strfmt = '(" HF energies are read from DAT file",/5(I5,F10.6))'
            Write(11,strfmt) (i,Eps(i),i=Nso+1,Nsu)
        End If
        Close(unit=12)

        ! Maximal number of eigenvectors Nlv:
        nlmax=IPlv
        If (Kv >= 3) Then
            nlmax=IPlv/3
        End If
        If (Nlv > nlmax) Nlv=nlmax
        Return

    End subroutine Init

    Subroutine RintS
        ! Reading of files SGC.CON and SCRC.CON with the self-energy and screening radial integrals.
        Implicit None
        Integer :: na, nb, ierr, k, nsh, khot, la, lb, ja, jb, na1, nb1, ind, &
                   nso1, khot1, k1, nsx1, nsx2, nav, kbox, i, ns1, idummy, err_stat, Nmax1, Lmax1
        Real :: x, y
        Character(Len=256) :: strfmt, err_msg
        Character(Len=3) :: key1, key2
        
        If (Ksig == 0) Return
        nsh=Nso+1

        ! Reading matrix element of Sigma from SGC.CON
        If (Kdsig /= 0) Then
            strfmt = '(5X,"Give valence energy: ",$)'
            Write(*,strfmt)
            Read(*,*) E_0
            Kexn=1
            strfmt = '(5X,"Extrapolation variant ","(1-normal,2-one side,3-nonlin.): ",I2)'
            Write(*,strfmt) Kexn
        End If

        Open(unit=18,file='SGC.CON',status='OLD',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file SGC.CON is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        Else
            Write(*,*)' Reading file SGC.CON ...'
        End If
        
        strfmt = '(7X,I3,5X,I1,24X,I1)'
        Read(18,strfmt) NmaxS,LmaxS,khot

        Do na=nsh,NmaxS
            la=Ll(na)
            ja=Jj(na)
            If (la > LmaxS) Cycle
            Do nb=na,NmaxS
                lb=Ll(nb)
                jb=Jj(nb)
                If (lb /= la .or. jb /= ja) Cycle
                Read(18,*) idummy,na1,nb1,x,y,z
                NhintS=NhintS+1
                If (na1 /= na .or. nb1 /= nb) Then
                    Write(*,*)' wrong indices for ',NhintS,' sigma:'
                    Write(*,*)' expected ',na,nb,' got ',na1,nb1
                    Stop
                End If
            End Do
        End Do
        Rewind(18)
        Allocate(Iint1S(NhintS),Rsig(NhintS),Dsig(NhintS),Esig(NhintS),Scr(10))
        strfmt = '(7X,I3,5X,I1,24X,I1)'
        Read(18,strfmt) NmaxS,LmaxS,khot
        i=0
        Do na=nsh,NmaxS
            la=Ll(na)
            ja=Jj(na)
            If (la > LmaxS) Cycle
            Do nb=na,NmaxS
                lb=Ll(nb)
                jb=Jj(nb)
                If (lb /= la .or. jb /= ja) Cycle
                Read(18,*) idummy,na1,nb1,x,y,z
                i=i+1
                If (na1 /= na .or. nb1 /= nb) Then
                    Write(*,*)' wrong indices for ',NhintS,' sigma:'
                    Write(*,*)' expected ',na,nb,' got ',na1,nb1
                    Stop
                End If
                ind=Nx*(na-nsh)+(nb-nsh+1)
                Iint1S(i)=ind
                Rsig(i)=x
                Dsig(i)=y
                Esig(i)=z
            End Do
        End Do
        xscr=10.d0
        If (Ksig == 2) Then
            xscr=0.d0
            strfmt = '(5X,I2,5X,F8.5)'
            Read(18,strfmt)
            Do k=1,10
                Read(18,strfmt) k1,Scr(k)
                If (k1 /= k-1) Then
                    Write(*,*)' wrong order of screening coefficients'
                    Stop
                End If
                xscr=xscr+Scr(k)
            End Do
        End If
        Write(*,*)' Khot =',khot
        Close(unit=18)

        If (Ksig < 2) Return
        ierr=0

        Open (unit=13,file='SCRC.CON',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file SCRC.CON is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        Else
            Write(*,*)' Reading file SCRC.CON ...'
        End If

        Read (13) ns1,nso1,nsh,khot1,nsx1,nsx2,Ksym
        Read (13) Nmax1,Lmax1,Kmax,nav,kbox
        strfmt = '(4X,"Parameters from the file SCRC.CON:", &
                    /4X,"Ns=",I3," Nso=",I2," Nsh=",I2," Khot=",I1, &
                    " Nsx=",2I3,/4X," Ksym=",I1,"NmaxS=",I3," LmaxS=",I2, &
                    " Kmax=",I2," Nav=",I3," Kbox=",I1)'
        Write(11,strfmt) ns1,nso1,nsh,khot1,nsx1,nsx2,Ksym,Nmax1,Lmax1,Kmax,nav,kbox
        If (nso1 /= Nso) Then
            Write(*,*)' RintS warning: Nso=',Nso,' <> ',Nso1
            ierr=ierr+1
        End If

        If (NmaxS /= nmax1 .or. LmaxS /= lmax1 .or. khot /= khot1) Then
            strfmt = '(3X,"Difference between SGC.CON and SCRC.CON:", &
                        /3X,"nmaxS =",I3," (SGC) =",I3," (SCRC)", &
                        /3X,"lmaxS =",I3," (SGC) =",I3," (SCRC)", &
                        /3X,"khot =",I3," (SGC) =",I3," (SCRC)")'
            Write( *,strfmt) NmaxS,Nmax1,LmaxS,Lmax1,khot,khot1
            If (NmaxS < Nmax1 .or. LmaxS < Lmax1) ierr=ierr+1
        End If
        Read (13) NgintS,nrd
        If (nrd /= Nx*Nx) Then
            Write(key1,'(I3)') int(sqrt(real(nrd)))
            Write(key2,'(I3)') Nx
            Write(*,*)' RintS: SCRC.CON formed with IPx=',Trim(AdjustL(key1)), ', but CONF.INT with IPx=', Trim(AdjustL(key2))
            ierr=ierr+1
        End If
        If (ierr > 0) Then
            Write(*,*) ierr,' error(s) detected'
            Stop
        End If

        Allocate(Rint2S(NgintS),Iint2S(NgintS),Iint3S(NgintS),Dint2S(NgintS),Eint2S(NgintS),IntOrdS(nrd))

        If (nsx2 == 0) Then ! Nsum is used to skip screening integrals
            Nsum=4*nav
        Else
            Nsum=4*Nmax1
        End If

        Read (13) (Rint2S(i), i=1,NgintS)
        Read (13) (Iint2S(i), i=1,NgintS)
        Read (13) (Iint3S(i), i=1,NgintS)
        Read (13) (IntOrdS(i), i=1,nrd)
        Read (13) (Dint2S(i), i=1,NgintS)
        Read (13) (Eint2S(i), i=1,NgintS)
        Close(unit=13)
        Write(*,*)' File SCRC.CON is Read'

        If (dabs(xscr-10.d0) > 1.d-5) Then
            Write(*,*) ' Note: Screening coefficients will be used'
            Write(*,*) ' above NmaxS ',Nmax1,' and LmaxS',Lmax1
            Read(*,*)
        End If
        Return

    End subroutine RintS

    Subroutine FormD
        ! This subroutine evaluates HF energies of the determinants (for HintS).
        Implicit None
        Integer :: i, ie, ke, n
        Real(dp) :: Emin, Emax, x
        Character(Len=128) :: strfmt

        If (Ksig*Kdsig == 0) Return

        Allocate(Diag(Nd))

        Emin= 1.d99
        Emax=-1.d99
        Do n=1,Nd
            x=0.d0
            Do i=1,Ne
                ie=Iarr(i,n)
                ke=Nh(ie)
                x=x+Eps(ke) ! Eps - HFD energies read from CONF.DAT
            End Do
            x=E_0-x
            Diag(n)=Real(x,kind=type_real)
            Emin=min1(Emin,x)
            Emax=max1(Emax,x)
        End Do
        strfmt = '(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)'
        write( 6,strfmt) Emin,Emax
        write(11,strfmt) Emin,Emax

    End Subroutine FormD
    
    Subroutine AllocateFormHArrays(mype)
        Use mpi_f08
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype

        vaBinSize = 10000000
        Call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IPlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngaunt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NhintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngint, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NgintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. Allocated(Jz)) Allocate(Jz(Nst))
        If (.not. Allocated(Nh)) Allocate(Nh(Nst))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(In)) Allocate(In(Ngaunt))
        If (.not. Allocated(Gnt)) Allocate(Gnt(Ngaunt))
        If (.not. Allocated(num_gaunts_per_partial_wave)) Allocate(num_gaunts_per_partial_wave(Nlx+1))
        If (.not. Allocated(Rint1)) Allocate(Rint1(Nhint))
        If (Kbrt == 0) Then
            If (.not. Allocated(Rint2)) Allocate(Rint2(1,Ngint))
        Else
            If (.not. Allocated(Rint2)) Allocate(Rint2(2,Ngint))
        End If
        If (.not. Allocated(Iint1)) Allocate(Iint1(Nhint))
        If (.not. Allocated(Iint2)) Allocate(Iint2(Ngint))
        If (.not. Allocated(Iint3)) Allocate(Iint3(Ngint))
        If (.not. Allocated(IntOrd)) Allocate(IntOrd(nrd))
        If (.not. Allocated(Iarr)) Allocate(Iarr(Ne,Nd))

        If (K_is /= 0) Then
            If (.not. Allocated(R_is)) Allocate(R_is(num_is))
            If (.not. Allocated(I_is)) Allocate(I_is(num_is))
        End If 

        If (Ksig /= 0) Then
            If (.not. Allocated(Scr)) Allocate(Scr(10))
            If (.not. Allocated(Rsig)) Allocate(Rsig(NhintS))
            If (.not. Allocated(Dsig)) Allocate(Dsig(NhintS))
            If (.not. Allocated(Esig)) Allocate(Esig(NhintS)) 
            If (.not. Allocated(Rint2S)) Allocate(Rint2S(NgintS))
            If (.not. Allocated(Dint2S)) Allocate(Dint2S(NgintS))
            If (.not. Allocated(Eint2S)) Allocate(Eint2S(NgintS))
            If (.not. Allocated(Iint1S)) Allocate(Iint1S(NhintS))
            If (.not. Allocated(Iint2S)) Allocate(Iint2S(NgintS))
            If (.not. Allocated(Iint3S)) Allocate(Iint3S(NgintS))
            If (.not. Allocated(IntOrdS)) Allocate(IntOrdS(nrd))
        End If

        If (mype==0) Then
            memFormH = 0_int64
            memFormH = sizeof(Nvc)+sizeof(Nc0) &
                + sizeof(Rint1)+sizeof(Rint2)+sizeof(Iint1)+sizeof(Iint2)+sizeof(Iint3)+sizeof(Iarr) &
                + sizeof(IntOrd)
            If (K_is /= 0) memFormH = memFormH+sizeof(R_is)+sizeof(I_is)
            If (Ksig /= 0) memFormH = memFormH+sizeof(Rint2S)+sizeof(Dint2S)+sizeof(Eint2S) &
                + sizeof(Iint1S)+sizeof(Iint2S)+sizeof(Iint3S) &
                + sizeof(Rsig)+sizeof(Dsig)+sizeof(Esig)+sizeof(IntOrdS)
        End If   

        Return
    End Subroutine AllocateFormHArrays

    Subroutine calcMemStaticArrays
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        
        memStaticArrays = memStaticArrays! + 209700000_int64 ! static "buffer" of 200 MiB 
        memStaticArrays = memStaticArrays + sizeof(Nn)+sizeof(Kk)+sizeof(Ll)+sizeof(Jj)+sizeof(Nf0) &
                        + sizeof(Jt)+sizeof(Njt)+sizeof(Eps)+sizeof(Diag)+sizeof(Ndc)+sizeof(Jz) &
                        + sizeof(Nh)+sizeof(In)+sizeof(Gnt)+sizeof(C)+(8+type_real)*vaBinSize

    End Subroutine calcMemStaticArrays

    Subroutine calcMemReqs
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer(Kind=int64) :: mem
        Integer(Kind=int64) :: bytesInteger, bytesDP, bytesReal
        Character(Len=16) :: memStr
    
        bytesInteger = 4_int64
        bytesReal = type_real
        bytesDP = 8_int64
    
        Call calcMemStaticArrays
        Call FormattedMemSize(memStaticArrays, memStr)
        Write(*,'(A,A,A)') 'calcMemReqs: Allocating static arrays will require at least ',Trim(memStr),' of memory per core. (These will not be deallocated)' 

        Call FormattedMemSize(memFormH, memStr)
        Write(*,'(A,A,A)') 'calcMemReqs: Allocating arrays for FormH will require at least ',Trim(memStr),' of memory per core' 
        
        memEstimate = memFormH + memStaticArrays
    
        mem = bytesReal * Nd * IPlv        & ! ArrB
            + bytesReal * Nlv * 2_int64    & ! Tk,Tj
            + bytesReal * IPlv * IPlv      & ! P
            + bytesReal * IPlv * 2_int64   & ! D,E
            + bytesReal * Nd * 2_int64     & ! B1,B2
            + bytesReal * Nd0 * Nd0        & ! Z1
            + bytesReal * Nd0                ! E1
    
        memDvdsn = mem
        Call FormattedMemSize(mem, memStr)
        Write(*,'(A,A,A)') 'calcMemReqs: Allocating arrays for Davidson procedure will require at least ', &
                            Trim(memStr),' of memory per core' 
    
        Call FormattedMemSize(memTotalPerCPU, memStr)
        If (memTotalPerCPU == 0) Then
            Write(*,'(A)') 'calcMemReqs: Available memory was not saved to the environment.'
        Else
            Write(*,'(A,A,A)') 'calcMemReqs: Total memory available to the job is ',Trim(memStr),' of memory per core' 
        End If
        Return
    End Subroutine calcMemReqs

    Subroutine InitFormH
        ! this subroutine initializes variables used for FormH and subsequent subroutines
        ! All necessary variables are broadcasted from root to all cores
        Use mpi_f08
        Use mpi_utils
        Implicit None
        
        Integer :: mpierr
        Integer(Kind=int64) :: count

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(KLSJ, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Crt4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Njd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kherr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kgerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kecp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_prj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nsu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Mj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NmaxS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(LmaxS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksym, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nsum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(C_is, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(XJ_av, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_sms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kdsig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kexn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Eps, IPs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(In, Ngaunt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gnt, Ngaunt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(num_gaunts_per_partial_wave, Nlx+1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint1, Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        If (Kbrt == 0) Then
            count = Ngint
        Else
            count = Ngint*2_int64
        End If
        Call BroadcastD(Rint2, count, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint1, Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastI(Iint2, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastI(Iint3, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IntOrd, Nx*Nx, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag, Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd,kind=int64)
        Call BroadcastI(Iarr, count, 0, 0, MPI_COMM_WORLD, mpierr)
        If (K_is /= 0) Then
            Call MPI_Bcast(R_is, num_is, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(I_is, num_is, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
        If (Ksig /= 0) Then
            Call MPI_Bcast(Scr, 10, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rsig, NhintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dsig, NhintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Esig, NhintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint1S, NhintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint2S, NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint3S, NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rint2S, NgintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dint2S, NgintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Eint2S, NgintS, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(IntOrdS, Nx*Nx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End subroutine InitFormH

    Subroutine FormH(npes, mype)
        Use mpi_f08
        Use str_fmt, Only : FormattedMemSize, FormattedTime
        Use determinants, Only : calcNd0, Gdet, CompCD, Rspq_phase1, Rspq_phase2
        Use matrix_io
        Use vaccumulator

        Implicit None

        Integer :: npes, mype, mpierr
        Integer :: k1, kx, n, ic, j, k, diff, icomp, ih4, ihmax
        Integer :: nn, kk, msg, sender, num_done, an_id, endnd
        Type(MPI_STATUS) :: status
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Integer(Kind=int64), Allocatable, Dimension(:) :: cntarray
        Integer(Kind=int64)     :: stot, s1, s2, numzero=0, nz0, maxme, maxNumElementsPerCore, mesplit, n8
        Real(kind=type_real)  :: t, tt
        Integer(Kind=int64) :: statmem, mem, maxmem, counter1, counter2, counter3
        Character(Len=16)     :: memStr, memStr2, memStr3, memStr4, memStr5, memTotStr, memTotStr2, counterStr, counterStr2, timeStr
        Integer :: iSign, iIndexes(3), jIndexes(3), nnd
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: vaGrowBy, ndGrowBy, ndsplit, ndcnt
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002

        Call startTimer(stot)
        Call InitFormH 
        Call calcNd0(Nc1, Nd0)

        If (mype == 0) Then
            If (Ksig == 2) Write(*,*) 'Screening is included'
            Call calcMemReqs
        End If

        ! Read number of processors
        If (Kl == 3) Then
            Open(66,file='progress.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Read(66) Nc_prev, Nd_prev
            Close(66) 
            If (Nd == Nd_prev) Then
                Kl = 1
                If (mype == 0) Then
                    print*, 'previously: Nc=',Nc_prev, ' Nd=', Nd_prev
                    print*, 'No new configurations to include'
                End If
            End If
        End If

        ! If Hamiltonian has already been fully constructed
        If (Kl == 1) Then 
            ! Read the Hamiltonian from file CONFp.HIJ
            Call ReadMatrix(Hamil%ind1,Hamil%ind2,Hamil%val,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)
            numzero = Count(Hamil%val(1:ih4) == 0)

            ! Add maximum memory per core from storing H to total memory count
            Call MPI_AllReduce(ih4, ihmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            memEstimate = memEstimate + ihmax*(8_int64+type_real)

        ! If Hamiltonian has not been fully constructed
        Else 
            ! Read the previous Hamiltonian from file CONFp.HIJ
            !!! Details of future implementation:
            !!!    the matrix elements will have to be saved to arrays iva1, iva2, rva1
            !!!    so those arrays can be extended when calculating new matrix elements 
            !Call ReadMatrix(Hamil%ind1,Hamil%ind2,Hamil%val,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)
            If (Kl == 3) Then
                Call ReadMatrix(iva1%vAccum,iva2%vAccum,rva1%vAccum,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)

                nz0 = count(rva1%vAccum==0)
                Call MPI_AllReduce(MPI_IN_PLACE, nz0, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
            
                If (mype == 0) print*, 'previously: Nc=',Nc_prev, ' Nd=', Nd_prev, ' NumH=', NumH-nz0, ' numzero=', nz0
            End If

            Allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne),cntarray(3))
            vaGrowBy = vaBinSize
            ndGrowBy = 1

            If (mype==0) Then
                Write(counterStr,fmt='(I16)') vaGrowBy
                Write(counterStr2,fmt='(I16)') ndGrowBy
                Write(*,'(A)') ' vaGrowBy = '//Trim(AdjustL(counterStr))//', ndGrowBy = '//Trim(AdjustL(counterStr2))
                print*, '========== Starting comparison stage of FormH =========='
            End If

            cntarray=0

            ! Get accumulator vectors setup (or re-setup if this is rank 0):
            If (Kl == 3) Then
                Call IVAccumulatorContinue(iva1, vaGrowBy)
                Call IVAccumulatorContinue(iva2, vaGrowBy)
                If (mype == 0) Call RVAccumulatorContinue(rva1, vaGrowBy)
            Else
                Call IVAccumulatorInit(iva1, vaGrowBy)
                Call IVAccumulatorInit(iva2, vaGrowBy)
                If (mype==0) Call RVAccumulatorInit(rva1, vaGrowBy)
                Nd_prev = 0
                ih4=0
                NumH=0_int64
                counter1=1_int64
                counter2=1_int64
                counter3=1_int64
            End If
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            If (mype == 0) Then        
                ! Distribute a portion of the workload of size ndGrowBy to each worker process
                Do an_id = 1, npes - 1
                   nnd = Nd_prev + 1 + ndGrowBy*(an_id-1) + 1
                   Call MPI_SEND( nnd, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                End Do

                n=Nd_prev+1
                Call Gdet(n,idet1)
                k=0
                Do ic=1,Nc 
                    kx=Ndc(ic)
                    If (k+kx > n) kx=n-k
                    If (kx /= 0) Then
                        Call Gdet(k+1,idet2)
                        Call CompCD(idet1,idet2,icomp)
                        If (icomp > 2) Then
                            k=k+kx
                        Else
                            Do k1=1,kx
                                k=k+1
                                Call Gdet(k,idet2)
                                Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                If (diff <= 2) Then
                                    nn=n
                                    kk=k
                                    Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                    If (Kdsig /= 0 .and. diff <= 2) E_k=Diag(kk)
                                    tt=Hmltn(idet1, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                                    If (tt /= 0) Then
                                        cntarray = cntarray + 1
                                        Call IVAccumulatorAdd(iva1, nn)
                                        Call IVAccumulatorAdd(iva2, kk)
                                        Call RVAccumulatorAdd(rva1, tt)
                                    End If
                                End If
                            End Do
                        End If
                    End If
                End Do

                Call startTimer(s1)

                NumH =  NumH + cntarray(1)
                num_done = 0
                If (Kl == 3) Then
                    ndsplit = (Nd-Nd_prev+1)/10
                    ndcnt = Nd_prev+1+ndsplit
                Else
                    ndsplit = Nd/10
                    ndcnt = ndsplit
                End If
                maxme = cntarray(2)
                j=9

                Do 
                    Call MPI_RECV( cntarray, 3, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    sender = status%MPI_SOURCE
             
                    !If (devmode == 1) print*, cntarray(3), cntarray(1), sender
                    If (nnd + ndGrowBy <= Nd) Then
                        nnd = nnd + ndGrowBy
                        Call MPI_SEND( nnd, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                    Else
                        msg = -1
                        Call MPI_SEND( msg, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If
            
                    NumH = NumH + cntarray(1)
                    maxme = max(cntarray(2),maxme)
                    mem = NumH * (8_int64+type_real)
                    maxmem = maxme * (8_int64+type_real)
                    statmem = memEstimate + memDvdsn - memFormH + maxmem
                    Call FormattedMemSize(statmem, memTotStr)
                    Call FormattedMemSize(memTotalPerCPU, memTotStr2)

                    If (nnd == ndcnt .and. j /= 0) Then
                        Call stopTimer(s1, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Call FormattedMemSize(NumH*(8+type_real), memStr3)
                        Write(counterStr,fmt='(I16)') NumH
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH comparison stage:', (10-j)*10, '% done in '// trim(timeStr)// '; '// &
                                                    Trim(AdjustL(counterStr)) // ' elements'
                        Write(*,'(4X,A)'), 'Memory: (HamiltonianTotal='// trim(memStr3)//', HamiltonianMaxMemPerCore='// trim(memStr2)//')'
                        If (memTotalPerCPU /= 0 .and. statmem > memTotalPerCPU) Then
                            Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but only ', &
                                                    Trim(memTotStr2) ,' is available.'
                            Stop
                        End If
                        j=j-1
                        ndcnt = ndcnt + ndsplit
                    End If
                    
                    If (num_done == npes-1) Then
                        Call stopTimer(s1, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Call FormattedMemSize(NumH*(8+type_real), memStr3)
                        Call FormattedMemSize(memStaticArrays, memStr4)
                        Call FormattedMemSize(memDvdsn, memStr5)
                        mem = memEstimate + memDvdsn - memFormH + maxmem
                        memEstimate = memEstimate + maxmem
                        Write(counterStr,fmt='(I16)') NumH
                        Call FormattedMemSize(mem, memStr)
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH comparison stage:', (10-j)*10, '% done in '// trim(timeStr)// '; '// &
                                                    Trim(AdjustL(counterStr)) // ' elements'
                        Write(*,'(4X,A)'), 'Memory: (HamiltonianTotal='// trim(memStr3)//', HamiltonianMaxMemPerCore='// trim(memStr2)//')'
                        Write(*,'(A)'), 'SUMMARY - (total = '// trim(memStr) // ', static = ' // trim(memStr4) &
                                            // ', davidson = ' // trim(memStr5) // ', Hamiltonian = ' // trim(memStr2) // ')'
                        If (memTotalPerCPU /= 0) Then
                            If (statmem > memTotalPerCPU) Then
                                Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but only ', &
                                                        Trim(memTotStr2) ,' is available.'
                                Stop
                            Else If (statmem < memTotalPerCPU) Then
                                Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, and ' , &
                                                        Trim(memTotStr2) ,' is available.'
                            End If
                        Else
                            Write(*,'(2X,A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, &
                                                        but available memory was not saved to environment'
                        End If
                        Exit
                    End If
                End Do
            Else
                cntarray=Int(ih4, kind=int64)
                Do 
                    Call MPI_RECV ( nnd, 1 , MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    If (nnd == -1) Then
                          Exit
                    Else
                        If (Nd - nnd < ndGrowBy) Then
                            endnd = Nd
                        Else
                            endnd = nnd+ndGrowBy-1
                        End If
                        cntarray(1)=0
                        Do n=nnd,endnd
                            Call Gdet(n,idet1)
                            k=0
                            Do ic=1,Nc 
                                kx=Ndc(ic)
                                If (k+kx > n) kx=n-k
                                If (kx /= 0) Then
                                    Call Gdet(k+1,idet2)
                                    Call CompCD(idet1,idet2,icomp)
                                    If (icomp > 2) Then
                                        k=k+kx
                                    Else
                                        Do k1=1,kx
                                            k=k+1
                                            Call Gdet(k,idet2)
                                            Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                            If (diff <= 2) Then
                                                nn=n
                                                kk=k
                                                cntarray = cntarray + 1
                                                Call IVAccumulatorAdd(iva1, nn)
                                                Call IVAccumulatorAdd(iva2, kk)
                                            End If
                                        End Do
                                    End If
                                End If
                            End Do
                        End Do 

                        cntarray(3) = nnd
                        Call MPI_SEND(cntarray, 3, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                    End if
                End Do
            End If

            Call IVAccumulatorCopy(iva1, Hamil%ind1, counter1)
            Call IVAccumulatorCopy(iva2, Hamil%ind2, counter2)
            If (mype == 0) Call RVAccumulatorCopy(rva1, Hamil%val, counter3)

            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
            If (mype == 0) Call RVAccumulatorReset(rva1)
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call startTimer(s1)

            Call MPI_AllReduce(counter1, maxNumElementsPerCore, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
            If (Kl == 3) Then
                mesplit = (maxNumElementsPerCore-ih4+1)/10
            Else
                mesplit = maxNumElementsPerCore/10
            End If
            numzero=0_int64
            j=1

            If (mype /= 0) Then
                Allocate(Hamil%val(counter1))
                Call startTimer(s2)
                If (Kl==3) Hamil%val(1:ih4) = rva1%vAccum(1:ih4)
                Do n8=ih4+1,counter1
                    nn=Hamil%ind1(n8)
                    kk=Hamil%ind2(n8)
                    Call Gdet(nn,idet1)
                    Call Gdet(kk,idet2)
                    Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                    Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                    If (Kdsig /= 0 .and. diff <= 2) E_k=Diag(kk)
                    t=Hmltn(idet1, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                    Hamil%val(n8)=t
                    If (Kl /=3 .and. t == 0) numzero=numzero+1
                    If (counter1 == maxNumElementsPerCore .and. mod(n8,mesplit)==0) Then
                        Call stopTimer(s1, timeStr)
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH calculation stage:', j*10, '% done in '// trim(timeStr)
                        j=j+1
                    End If
                End Do
                Call stopTimer(s2, timeStr)
            Else
                print*, '========== Starting calculation stage of FormH =========='
            End If
            If (Kl==3) numzero = count(Hamil%val==0)
            Deallocate(idet1, idet2, iconf1, iconf2, cntarray)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            If (mype==0) print*, '========== Formation of Hamiltonian matrix completed =========='
        End If

        ih8=size(Hamil%val, kind=int64)
        ih4=ih8
        
        Call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, numzero, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, iscr, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, xscr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        
        ! Write Hamiltonian to file CONFp.HIJ
        If (Kl /= 1 .and. Kw == 1)  Call WriteMatrix(Hamil,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)
        
        ! give all cores Hmin, the minimum matrix element value
        Call MPI_AllReduce(Hamil%val(1:ih8), Hamil%minval, 1, mpi_type_real, MPI_MIN, MPI_COMM_WORLD, mpierr)
        
        If (mype==0) Then
            ! Write number of non-zero matrix elements
            Write(counterStr,fmt='(I16)') NumH-numzero
            strfmt = '(4X,"NumH = ",A)'
            Write( 6,strfmt) Trim(AdjustL(counterStr))
            Write(11,strfmt) Trim(AdjustL(counterStr))

            ! Write number of zero-valued matrix elements
            Write(counterStr,fmt='(I16)') numzero
            strfmt = '(4X,"numzero = ",A)'
            Write( 6,strfmt) Trim(AdjustL(counterStr))
            Write(11,strfmt) Trim(AdjustL(counterStr))

            ! Write value of lowest matrix element
            Write(counterStr,fmt='(F14.8)') Hamil%minval
            strfmt = '(4X,"Hmin = ",A)'
            Write( 6,strfmt) Trim(AdjustL(counterStr))
            Write(11,strfmt) Trim(AdjustL(counterStr))

            ! Write number of screened integrals
            If (Ksig == 2 .and. iscr > 0) Then
                xscr=xscr/iscr
                Write(counterStr,fmt='(I16)') iscr
                strfmt = '(5X,"For ",A," integrals averaged screening ",F8.5)'
                Write ( 6,strfmt) Trim(AdjustL(counterStr)),xscr
                Write (11,strfmt) Trim(AdjustL(counterStr)),xscr
                If (Kherr+Kgerr > 0) Then
                    strfmt = '(4X,"Extrapolation warning: small denominators.", &
                    /4X,"HintS: ",I6,"; GintS: ",I7)'
                    Write ( 6,strfmt) Kherr,Kgerr
                    Write (11,strfmt) Kherr,Kgerr
               End If
            End If
        End If
    
        ! Stop timer and output computation time of FormH 
        If (mype == 0) Then
            Call stopTimer(stot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> FormH took '// trim(timeStr) // ' to complete'
        End If

        Return
    End Subroutine FormH

    Real(kind=type_real) function Hmltn(idet, is, nf, i2, i1, j2, j1) 
        ! This function calculates the Hamiltonian matrix element between determinants idet1 and idet2
        Use determinants, Only : Rspq
        Use integrals, Only : Gint, Hint
        Use formj2, Only : F_J0, F_J2
        Implicit None
        Integer, Allocatable, Dimension(:), Intent(InOut)   :: idet
        Integer, Intent(InOut)                              :: is, nf, i1, i2, j1, j2
        
        Integer     :: iq, jq, jq0
        Real(kind=type_real)    :: t
        
        t=0.d0
        Select Case(nf)
            Case(2) ! determinants differ by two functions
                t=t+Real(Gint(i2,j2,i1,j1),kind=type_real)*is 
                t=t+Gj*F_J2(idet,is,nf,i2,i1,j2,j1)
            Case(1) ! determinants differ by one function
                Do iq=1,Ne
                    i1=idet(iq)
                    If (i1 /= j1) t=t+Real(Gint(j2,i1,j1,i1),kind=type_real)*is
                End Do
                t=t+Real(Hint(j2,j1),kind=type_real)*is
            Case(0) ! determinants are equal
                Do iq=1,Ne
                    i1=idet(iq)
                    jq0=iq+1
                    If (jq0 <= Ne) Then
                        Do jq=jq0,Ne
                            j1=idet(jq)
                            t=t+Real(Gint(i1,j1,i1,j1),kind=type_real)*is
                        End Do
                    End If
                    t=t+Real(Hint(i1,i1),kind=type_real)*is
                End Do
                t=t+Real(Gj,kind=type_real)*F_J2(idet,is,nf,i2,i1,j2,j1)
        End Select
        Hmltn=t
        Return
    End function Hmltn

    Subroutine DeAllocateFormHArrays(mype)
        Use mpi_f08
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype
        Character(Len=16) :: memStr

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    
        If (mype==0) Then
            Call FormattedMemSize(memFormH, memStr)
            Write(*,'(A,A,A)') 'De-allocating ',Trim(memStr),' of memory per core from arrays for FormH'  
        End If
    
        !If (Allocated(Nvc)) Deallocate(Nvc)
        !If (Allocated(Nc0)) Deallocate(Nc0)
        If (Allocated(Rint1)) Deallocate(Rint1)
        If (Allocated(Rint2)) Deallocate(Rint2)
        If (Allocated(Iint1)) Deallocate(Iint1)
        If (Allocated(Iint2)) Deallocate(Iint2)
        If (Allocated(Iint3)) Deallocate(Iint3)
        If (Allocated(Rint2S)) Deallocate(Rint2S)
        If (Allocated(Dint2S)) Deallocate(Dint2S)
        If (Allocated(Eint2S)) Deallocate(Eint2S)
        If (Allocated(Iint1S)) Deallocate(Iint1S)
        If (Allocated(Iint2S)) Deallocate(Iint2S)
        If (Allocated(Iint3S)) Deallocate(Iint3S)
        If (Allocated(Rsig)) Deallocate(Rsig)
        If (Allocated(Dsig)) Deallocate(Dsig)
        If (Allocated(Esig)) Deallocate(Esig) 
        If (Allocated(R_is)) Deallocate(R_is)
        If (Allocated(I_is)) Deallocate(I_is)
        If (Allocated(IntOrd)) Deallocate(IntOrd)
        If (Allocated(IntOrdS)) Deallocate(IntOrdS)
        If (Allocated(Iarr)) Deallocate(Iarr)
        If (Allocated(Scr)) Deallocate(Scr)

    End Subroutine DeAllocateFormHArrays

    Subroutine AllocateDvdsnArrays(mype)
        Use mpi_f08
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer :: mpierr, mype
        Character(Len=16) :: memStr, memTotStr
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(ArrB)) Allocate(ArrB(Nd,IPlv))
        If (.not. Allocated(Tk)) Allocate(Tk(Nlv))
        If (.not. Allocated(Tj)) Allocate(Tj(Nlv))
        If (.not. Allocated(P)) Allocate(P(2*Nlv,2*Nlv))
        If (.not. Allocated(E)) Allocate(E(Nlv))
        If (.not. Allocated(Iconverge)) Allocate(Iconverge(Nlv))
        If (.not. Allocated(B1)) Allocate(B1(Nd))
        If (.not. Allocated(B2)) Allocate(B2(Nd))
        If (.not. Allocated(Z1)) Allocate(Z1(Nd0,Nd0))
        If (.not. Allocated(E1)) Allocate(E1(Nd0))
        If (.not. Allocated(Jt)) Allocate(Jt(Njd))
        
        Call MPI_Bcast(Jt, Njd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        If (KLSJ == 1) Then
            If (.not. Allocated(xj)) Allocate(xj(Nlv))
            If (.not. Allocated(xl)) Allocate(xl(Nlv))
            If (.not. Allocated(xs)) Allocate(xs(Nlv))
        End If

        If (mype==0) Then
            memDvdsn = 0_int64
            memDvdsn = sizeof(ArrB)+sizeof(Tk)+sizeof(Tj)+sizeof(P)+sizeof(E) &
                + sizeof(Iconverge)+sizeof(B1)+sizeof(B2)+sizeof(Z1)+sizeof(E1)
            Call FormattedMemSize(memDvdsn, memStr)
            Write(*,'(A,A,A)') 'Allocating arrays for Davidson procedure requires ',Trim(memStr),' of memory per core'  
            memEstimate = memEstimate - memFormH + memDvdsn
            Call FormattedMemSize(memEstimate, memStr)
            Write(*,'(A,A,A)') 'Total memory estimate for Davidson procedure is ',Trim(memStr),' of memory per core' 
            Call FormattedMemSize(memTotalPerCPU, memTotStr)
            If (memTotalPerCPU /= 0) Then
                If (memEstimate > memTotalPerCPU) Then
                    Write(*,'(A,A,A,A)'), 'At least '// Trim(memStr), ' is required to finish conf, but only ', &
                                            Trim(memTotStr) ,' is available.'
                    Stop
                Else If (memEstimate < memTotalPerCPU) Then
                    Write(*,'(A,A,A,A)'), 'At least '// Trim(memStr), ' is required to finish conf, and ' , &
                                            Trim(memTotStr) ,' is available.'
                End If
            Else
                Write(*,'(2X,A,A,A,A)'), 'At least '// Trim(memStr), ' is required to finish conf, &
                                            but available memory was not saved to environment'
            End If
        End If   

    End Subroutine AllocateDvdsnArrays

    Subroutine DiagInitApprox(mype, npes)
        Use mpi_f08
        Use str_fmt, Only : startTimer, stopTimer, FormattedTime
        Use env_var
        Implicit None
        External :: INFOG2L, PDELGET, SSYEV, DSYEV, BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT, NUMROC, DESCINIT, PDELSET, PDSYEVD, PSSYEVD, PSELSET

        Integer :: I, J, lwork, mype, npes, ifail, srcRow, srcCol
        Integer :: NPROW, NPCOL, CONTEXT, MYROW, MYCOL, LIWORK, NROWSA, NCOLSA, NUMROC
        Integer(Kind=int64) :: s1, scalapackThreshold
        Integer(Kind=int64), Dimension(2) :: blacsGrid
        Logical :: isVerbose, shouldForceSerial, shouldForceScalapack
        Real(type_real), Dimension(4) :: realtmp
        Integer, Dimension(4) :: itmp
        Integer, Dimension(9) :: DESCA, DESCZ
        Integer, Allocatable, Dimension(:)     :: IWORK ! integer work array
        Real(type_real), Allocatable, Dimension(:)    :: W     ! double precision work array
        Real(type_real), Allocatable, Dimension(:,:)  :: A, Z
        Character(len=48) :: blacsGridSpec

        If (mype == 0) Then
            ! At what point should we switch to the BLACS/ScaLAPACK algorithm?  The default
            ! value 10240 is just a guess and equates with an 800 MiB Z1:
            scalapackThreshold = GetEnvInteger64("CONF_DIAG_INIT_APPROX_THRESHOLD", 10240)
            isVerbose = GetEnvLogical("CONF_DIAG_INIT_APPROX_VERBOSE", .false.)
            shouldForceSerial = GetEnvLogical("CONF_DIAG_INIT_APPROX_FORCE_SERIAL", .false.)
            shouldForceScalapack = GetEnvLogical("CONF_DIAG_INIT_APPROX_FORCE_SCALAPACK", .false.)
            If (shouldForceSerial .and. shouldForceScalapack) shouldForceScalapack = .false.
        End If
        Call MPI_Bcast(scalapackThreshold, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(isVerbose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(shouldForceSerial, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(shouldForceScalapack, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        
        ! Start timer for initial diagonalization
        Call startTimer(s1)

        ! If running in serial or size of initial approximation is below threshold, 
        If (.not. shouldForceScalapack .and. (npes == 1 .or. Nd0 <= scalapackThreshold .or. shouldForceSerial)) Then
            If (mype==0) Then
                If (npes==1 .or. shouldForceSerial) Then
                    Write(*,"(A)") "DiagInitApprox: ScaLAPACK not used, running serially"
                Else
                    Write(*,"(A,I12,A,I12)") "DiagInitApprox: ScaLAPACK not used, ", Nd0, " <= ", scalapackThreshold
                End If
            End If
            Select Case(type_real)
            Case(sp)
                Call SSYEV('V','U',Nd0,Z1,Nd0,E1,realtmp,-1,ifail)
            Case(dp)
                Call DSYEV('V','U',Nd0,Z1,Nd0,E1,realtmp,-1,ifail)
            End Select
            lwork = Nint(realtmp(1))
            Allocate(W(lwork))
            If (mype==0 .and. isVerbose) Write(*,"(A,I12)") "DiagInitApprox: Work array allocated, real=", lwork
            Select Case(type_real)
            Case(sp)
                Call SSYEV('V','U',Nd0,Z1,Nd0,E1,W,lwork,ifail)
            Case(dp)
                Call DSYEV('V','U',Nd0,Z1,Nd0,E1,W,lwork,ifail)
            End Select
            If (mype==0 .and. isVerbose) Write(*,"(A)") "DiagInitApprox: Eigenvalues and -vectors calculated"
            Deallocate(W)
        Else
            If (mype==0) Write(*,"(A,I5,A)") "DiagInitApprox: ScaLAPACK used, will distribute across ", npes, " cpus"
            NPROW=1
            NPCOL=npes

            ! Initialize BLACS context
            CALL BLACS_GET(-1, 0, CONTEXT)
            CALL BLACS_GRIDINIT(CONTEXT, 'R', NPROW, NPCOL)
            CALL BLACS_GRIDINFO(CONTEXT, NPROW, NPCOL, MYROW, MYCOL)
            If (isVerbose) Write(*,"(A,I5,A,I12,A,I12)") "DiagInitApprox: BLACS grid init, rank ", mype, " -> ", MYROW,  ",", MYCOL

            ! Calculate the number of rows, NROWSA, and the number of columns, NCOLSA, for the local matrices A
            If (mype==0) Then
                blacsGrid(1) = GetEnvInteger("CONF_DIAG_INIT_APPROX_BLACS_ROWS", Nd0)
                blacsGrid(2) = GetEnvInteger("CONF_DIAG_INIT_APPROX_BLACS_COLS", 1)
                If (isVerbose) Write(*,"(A,I5,A,I5,A)") "DiagInitApprox: BLACS grid will be ", blacsGrid(1), " rows by ", blacsGrid(2), " columns"
            End If
            Call MPI_Bcast(blacsGrid, 2, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
            
            NROWSA = NUMROC(Nd0,blacsGrid(1),MYROW,0,NPROW)
            NCOLSA = NUMROC(Nd0,blacsGrid(2),MYCOL,0,NPCOL)
            ALLOCATE(A(NROWSA,NCOLSA), Z(NROWSA,NCOLSA))
            If (isVerbose) Write(*,"(A,I5,A,I12,A,I12)") "DiagInitApprox: Allocated local storage, rank ", mype, " -> ", NROWSA,  " x ", NCOLSA

            ! Initialize array descriptors DESCA and DESCZ
            CALL DESCINIT(DESCA, Nd0, Nd0, blacsGrid(1), blacsGrid(2), 0, 0, CONTEXT, NROWSA, ifail)
            CALL DESCINIT(DESCZ, Nd0, Nd0, blacsGrid(1), blacsGrid(2), 0, 0, CONTEXT, NROWSA, ifail)

            If (mype==0 .and. isVerbose) Write(*,"(A)") "DiagInitApprox: BLACS initialization completed"
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            ! Distribute the upper triangular part of the global matrix Z1 onto the local matrices A
            A=0.0d0
            Select Case(type_real)
            Case(sp)
                Do I=1,Nd0
                    Do J=1,I
                        CALL PSELSET(A, J, I, DESCA, Z1(I,J))
                    End Do
                End Do
            Case(dp)
                Do I=1,Nd0
                    Do J=1,I
                        CALL PDELSET(A, J, I, DESCA, Z1(I,J))
                    End Do
                End Do
            End Select
            If (mype==0 .and. isVerbose) Write(*,"(A)") 'DiagInitApprox: Matrix Z1 distributed'

            ! Query for the optimal work array sizes
            Select Case(type_real)
            Case(sp)
                CALL PSSYEVD('V', 'U', Nd0, A, 1, 1, DESCA, E1, Z, 1, 1, DESCZ, realtmp, -1, itmp, 1, ifail)
            Case(dp)
                CALL PDSYEVD('V', 'U', Nd0, A, 1, 1, DESCA, E1, Z, 1, 1, DESCZ, realtmp, -1, itmp, 1, ifail)
            End Select
            LWORK = realtmp(1)
            LIWORK = itmp(1)
            ALLOCATE(IWORK(LIWORK), W(LWORK))
            If (isVerbose) Write(*,"(A,I5,A,I12,A,I12)") "DiagInitApprox: Allocated work arrays, rank ", mype, " -> real=", LWORK,  ", integer=", LIWORK

            ! Compute the eigenvalues E1 and eigenvectors Z
            Select Case(type_real)
            Case(sp)
                CALL PSSYEVD('V', 'U', Nd0, A, 1, 1, DESCA, E1, Z, 1, 1, DESCZ, W, LWORK, IWORK, LIWORK, ifail)
            Case(dp)
                CALL PDSYEVD('V', 'U', Nd0, A, 1, 1, DESCA, E1, Z, 1, 1, DESCZ, W, LWORK, IWORK, LIWORK, ifail)
            End Select
            If (mype==0 .and. isVerbose) Write(*,"(A)") "DiagInitApprox: Eigenvalues and -vectors calculated"

            ! Un-distribute the eigenvector matrix back to Z1 in each worker:
            Do I=1,Nd0
                Do J=1,Nd0
                    CALL PDELGET('A', ' ', Z1(I,J), Z, I, J, DESCZ)
                End Do
            End Do
            If (mype==0 .and. isVerbose) Write(*,"(A)") 'DiagInitApprox: Matrix Z1 undistributed'
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        
            ! Exit BLACS context
            CALL BLACS_GRIDEXIT(CONTEXT)
            If (mype==0 .and. isVerbose) Write(*,"(A)") "DiagInitApprox: BLACS grid destroyed"

            ! Clean up
            Deallocate(IWORK, W, Z, A)
        End If
        If (mype==0) Then
            ! Stop timer and print time for initial diagonalization
            Call stopTimer(s1, timeStr)
            Write(*,"(2X,A,A,A)"), "TIMING >>> Initial diagonalization took ", trim(timeStr), " to complete"
        End If
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    End Subroutine DiagInitApprox

    Subroutine Diag4(mype, npes)
        ! this Subroutine executes the Davidson procedure for diagonalization
        ! the Subroutine Mxmpy is a computational bottleneck and was the only Subroutine to be parallelized
        ! all other subroutines are performed by the root core
        Use mpi_f08
        Use str_fmt, Only : startTimer, stopTimer, FormattedTime
        Use davidson
        Use formj2, Only : J_av
        Implicit None
        External :: DSYEV, SSYEV

        Integer  :: k1, k, i, n1, kx, i1, it, js, mype, npes, mpierr, ifail=0, lwork
        Real(kind=type_real), Dimension(4) :: realtmp
        Real(kind=type_real), Allocatable, Dimension(:) :: W ! work array
        Integer(Kind=int64) :: start_time, s1
        Real(type_real)  :: crit, ax, x, xx, vmax
        Real(dp) :: cnx
        Character(Len=16) :: timeStr, iStr
        Integer,  Allocatable, Dimension(:) :: Jn, Jk
        Real(kind=type_real), Allocatable, Dimension(:) :: Jv

        ! Initialize parameters and arrays
        Iconverge = 0
        crit=1.d-6
        ax=1.d0
        cnx=1.d0
        If (.not. Allocated(ax_array)) Allocate(ax_array(Nlv))

        ! Start timer for Davidson procedure
        Call startTimer(start_time)

        If (mype == 0) Then
            If (Kl4 == 0) Return ! Read CONF.XIJ and make 1 iteration
            If (Nc4 > Nc) Nc4=Nc
            Write (*,*) 'kl4=',Kl4,'  Nc4=',Nc4,'  Crt4=',Crt4
            Open(77,file='CONF.PRG',status='REPLACE',form='FORMATTED',action='WRITE')
        End If

        ! Construct initial approximation from Hamiltonian matrix
        Call Init4(mype)

        ! Diagonalize the initial approximation 
        Call DiagInitApprox(mype, npes)

        ! Write initial approximation to file CONF.XIJ
        Call FormB0(mype)
        
        ! Set up work array for diagonalization of energy matrix P during iterative procedure
        If (mype == 0) Then
            n1 = 2*Nlv
            Select Case(type_real)
            Case(sp)
                Call SSYEV('V','U',n1,P,n1,E,realtmp,-1,ifail)
            Case(dp)
                Call DSYEV('V','U',n1,P,n1,E,realtmp,-1,ifail)
            End Select
            lwork = Nint(realtmp(1))
            Allocate(W(lwork))
        End If

        ! temporary solution for value of Jsq%ind1 changing during LAPACK ZSYEV subroutine
        If (mype == 0) Then
            js = size(Jsq%ind1)
            Allocate(Jn(js),Jk(js),Jv(js))
            Jn=Jsq%ind1
            Jk=Jsq%ind2
            Jv=Jsq%val
        End If

        ! Davidson loop:
        kdavidson = 0
        If (mype==0) Write(*,*) 'Start with kdavidson =', kdavidson

        Do it=1,N_it
            Write(iStr,'(A)') it
            If (mype == 0) Write(77,'(A,I3,A)') '========== Iteration ', it , ' ==========' 
            ! Compare lowest admixture of an unconverged vector to convergence criteria
            If (ax > crit) Then
                If (mype == 0) Then
                    strfmt = '(1X,"** Davidson iteration ",I3," for ",I3," levels **")'
                    Write( 6,strfmt) it,Nlv
                    Write(11,strfmt) it,Nlv
                End if
                ! Orthonormalization of Nlv probe vectors
                If (mype==0) Then
                    Call startTimer(s1)
                    Do i=1,Nlv
                        Call Ortn(i,ifail)
                        If (ifail /= 0) Then
                            Write(*,*)' Fail of orthogonalization for ',i
                            Stop
                        End If
                    End Do
                    Call stopTimer(s1, timeStr)
                    Write(77,*) 'Orthonormalization of Nlv probe vectors took ' // timeStr
                    Call startTimer(s1)
                End If

                ! Formation of the left-upper block of the energy matrix P
                Call Mxmpy(1, mype)
                If (mype==0) Then 
                    Call stopTimer(s1, timeStr)
                    Write(77,*) 'Mxmpy1 took ' // timeStr
                    Call startTimer(s1)
                    Call FormP(1, vmax)
                    Call stopTimer(s1, timeStr)
                    Write(77,*) 'FormP1 took ' // timeStr
                    ! Average initial diagonal over relativistic configurations if projecting J
                    If (it == 1 .and. K_prj == 1) Call AvgDiag
                    
                    ! Formation of Nlv additional probe vectors
                    Call startTimer(s1)
                    cnx=0.d0
                    Do i=1,Nlv
                        i1=i+Nlv
                        Call Dvdsn(i,cnx)
                        If (Iconverge(i)==0) Then
                            Call Ortn(i1,ifail)
                            If (ifail /= 0) Then
                                Write(*,*)' Fail of orthogonalization for ',i1
                                Stop
                            End If
                        End If
                    End Do
                    Call stopTimer(s1, timeStr)
                    Write(77,*) 'Formation of Nlv additional probe vectors took' // timeStr
                End If
                If (K_prj == 1) Then
                    If (mype == 0) Then
                        Jsq%ind1=Jn
                        Jsq%ind2=Jk
                        Jsq%val=Jv
                    End If
                    Call Prj_J(Nlv+1,Nlv,2*Nlv+1,1.d-5,mype)
                    If (mype == 0) Then
                        Call startTimer(s1)
                        Do i=Nlv+1,2*Nlv
                            If (Iconverge(i-Nlv)==0) Then
                                Call Ortn(i,ifail)
                                If (ifail /= 0) Then
                                    If (mype == 0) Write(*,*)' Fail of orthogonalization 2 for ',i1
                                    If (kdavidson==1) Then
                                        Stop
                                    Else
                                        kdavidson=1
                                        If (mype == 0) Write(*,*) ' change kdavidson to ', kdavidson
                                        ifail=0
                                        Exit
                                    End If
                                End If
                            End If
                        End Do
                        Call stopTimer(s1, timeStr)
                        Write(77,*) 'Orthogonalization of Nlv additional probe vectors took ' // timeStr
                    End If
                End If
                Call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                If (cnx > Crt4) Then
                    ! Formation of other three blocks of the matrix P:
                    Call startTimer(s1)
                    Call Mxmpy(2, mype)
                    If (mype==0) Then
                        Call stopTimer(s1, timeStr)
                        Write(77,*) 'Mxmpy2 took ' // timeStr
                        Call startTimer(s1)
                        Call FormP(2, vmax)
                        Call stopTimer(s1, timeStr)
                        Write(77,*) 'FormP2 took ' // timeStr
                        ! >>>>> this block of code causes problem allocating xj, xl, xs arrays >>>>>
                        ! >>>>>>>>>>>>>>>> if allocation occurs after this block >>>>>>>>>>>>>>>>>>>
                        ! Evaluation of Nlv eigenvectors:
                        Call startTimer(s1)
                        Select Case(type_real)
                        Case(sp)
                            Call SSYEV('V','U',n1,P,n1,E,W,lwork,ifail)
                        Case(dp)
                            Call DSYEV('V','U',n1,P,n1,E,W,lwork,ifail)
                        End Select
                        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        Call stopTimer(s1, timeStr)
                        Write(77,*) 'Nlv eigenvectors evaluated in ' // timeStr
                        
                        ax=0.d0
                        vmax=-1.d10
                        Do i=1,Nlv
                            xx=0.d0
                            Do k=1,Nlv
                                k1=k+Nlv
                                x=dabs(P(k1,i))
                                If (xx <= x) Then
                                    xx=x
                                    kx=k
                                End If
                            End Do
                            If (ax < xx) ax=xx
                            If (vmax < E(i)) vmax=E(i)
                            strfmt = '(1X,"E(",I3,") =",F14.8,"; admixture of vector ",I3,": ",F10.7)'
                            Write( 6,strfmt) i,-(E(i)+Hamil%minval),kx,xx
                            Write(11,strfmt) i,-(E(i)+Hamil%minval),kx,xx
                            ax_array(i) = xx
                        End Do
                    End If
                    
                    Call MPI_Bcast(ax, 1, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
                    ! Write intermediate CONF.XIJ in frequency of kXIJ
                    Call startTimer(s1)
                    If (kXIJ > 0) Then
                        ! Only write each kXIJ iteration
                        If (mod(it, kXIJ) == 0) Then 
                            ! Calculate J for each energy level
                            Do n=1,Nlv
                                Call MPI_Bcast(ArrB(1:Nd,n), Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
                                If (mype == 0) Then
                                    Jsq%ind1=Jn
                                    Jsq%ind2=Jk
                                    Jsq%val=Jv
                                End If
                                Call J_av(ArrB(1,n),Nd,Tj(n),ierr)  ! calculates expectation values for J^2
                            End Do
                            If (mype == 0) Then
                                Call FormB
                                Call PrintEnergiesDvdsn(it)
                                Call PrintWeightsDvdsn(it)
                                print*, 'CONF.LVL updated'
                            End If
                        Else
                            If (mype == 0) Call FormBskip
                        End If
                    ! Else skip writing intermediate CONF.XIJ
                    Else
                        Call FormBskip
                    End If
                    Call stopTimer(s1, timeStr)
                    If (mype == 0) Write(77,*) 'Formation of next iteration of eigenvectors took ' // timeStr
                ! Davidson procedure is converged - write message and exit loop
                Else
                    If (mype == 0) Then
                        strfmt = '(" Davidson procedure converged")'
                        Write( 6,strfmt)
                        Write(11,strfmt)
                        ! Assign energies to array Tk for case when Davidson procedure is already converged 
                        If (it == 1) Tk=E(1:Nlv) 
                    End If
                    
                    Exit
                End If
            End If
        End Do

        ! temporary solution for value of Jsq%ind1 changing during LAPACK ZSYEV subroutine
        If (mype == 0) Then
            Jsq%ind1=Jn
            Jsq%ind2=Jk
            Jsq%val=Jv
            If (allocated(Jn)) Deallocate(Jn)
            If (allocated(Jk)) Deallocate(Jk)
            If (allocated(Jv)) Deallocate(Jv)
        End If
        ! Write final eigenvalues and eigenvectors to file CONF.XIJ
        Call WriteFinalXIJ(mype)

        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            Write(*,'(2X,A)'), 'TIMING >>> Davidson procedure took '// trim(timeStr) // ' to complete'
            Close(77)
        End If

        If (allocated(W)) Deallocate(W)
        If (allocated(Hamil%ind1)) Deallocate(Hamil%ind1)
        If (allocated(Hamil%ind2)) Deallocate(Hamil%ind2)
        If (allocated(Hamil%val)) Deallocate(Hamil%val)
        If (allocated(Jsq%ind1)) Deallocate(Jsq%ind1)
        If (allocated(Jsq%ind2)) Deallocate(Jsq%ind2)
        If (allocated(Jsq%val)) Deallocate(Jsq%val)
        If (allocated(Diag)) Deallocate(Diag)
        If (allocated(P)) Deallocate(P)
        If (allocated(B1)) Deallocate(B1)
        If (allocated(B2)) Deallocate(B2)
        If (allocated(Z1)) Deallocate(Z1)
        If (allocated(E1)) Deallocate(E1)
        If (allocated(Jt)) Deallocate(Jt)
        If (allocated(Njt)) Deallocate(Njt)

        Return
    End Subroutine Diag4

    Subroutine PrintEnergiesDvdsn(iter)
        ! This subroutine prints eigenvalues in increasing order
        Implicit None

        Integer :: j, ist, iter
        Real(kind=type_real) :: dt, del
        Real(kind=type_real), Dimension(Nlv) :: Tk0
        Character(Len=1), Dimension(10) :: stecp*7
        Character(Len=1), Dimension(4)  :: strsms*6
        Character(Len=1), Dimension(3)  :: strms*3
        Character(Len=512) :: strfmt
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/

        If (iter <= kXIJ) Then
            Open(unit=81,status='REPLACE',file='CONF.ENG',action='WRITE')
        Else
            Open(unit=81,status='UNKNOWN',POSITION='APPEND',file='CONF.ENG',action='WRITE')
            strfmt = '(A)'
            Write(81,strfmt) ''
        End If

        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7

        ! Print line of "=" for top of table
        strfmt = '(A,I3)'
        Write(81,strfmt) 'Davidson Iteration #', iter
        strfmt = '(4X,63("="))'
        Write(81,strfmt)

        ! If pure CI, print Nc, Nd, Gj
        If (Ksig*Kdsig == 0) Then
            strfmt = '(4X,"Energy levels (",A7," Nc=",I7," Nd=",I9,"); Gj =",F7.4, &
                        /4X,"N",6X,"JTOT",12X,"EV",16X,"ET",9X,"DEL(CM**-1)")'
            Write(81,strfmt) stecp(ist),Nc,Nd,Gj
        ! If CI+all-order/CI+MBPT, print E_0, Kexn, Nc, Nd, Gj
        Else
            strfmt = '(4X,"Energy levels ",A7,", Sigma(E =",F10.4,") extrapolation var.", &
                    I2,/4X,"(Nc=",I7," Nd=",I9,"); Gj =",F7.4,/4X,"N",6X,"JTOT",12X, &
                    "EV",16X,"ET",9X,"DEL(CM**-1)")'
            Write(81,strfmt) stecp(ist),E_0,Kexn,Nc,Nd,Gj
        End If

        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                strfmt = '(4X,"Volume shift: dR_N/R_N=",F9.5," Rnuc=",F10.7)'
                Write(81,strfmt) C_is,Rnuc
            Else
                strfmt = '(4X,A3,":",E9.2,"*(P_i Dot P_k) ",A6," Lower component key =",I2)'
                Write(81,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
            End If
        End If

        ! Print line of "-" to close header
        strfmt = '(4X,63("-"))'
        Write(81,strfmt)

        ! For each energy level, print the following:
        ! j = index for energy level
        ! Tj(j) = total angular momentum
        ! Tk(j) = valence energy of jth level
        ! DT = Tk(j)-Ecore = total energy
        ! DEL = Tk(1)-Tk(j) = energy difference between jth level and 1st level in cm**-1
        Do j=1,Nlv
            Tk0(j)=Tk(j)+4.d0*Gj*Tj(j)*(Tj(j)+1.d0)
            DT=Tk0(j)-Ecore
            DEL=(Tk0(1)-Tk0(j))*2*DPRy
            Write(81,'(3X,I3,F14.9,F14.8,F18.6,F15.2)') j,Tj(j),Tk0(j),DT,DEL
        End Do

        ! Print line of "=" for bottom of table
        strfmt = '(4X,63("="))'
        Write(81,strfmt)
        Close(81)
    End Subroutine PrintEnergiesDvdsn

    Subroutine PrintWeightsDvdsn(iter)
        ! This subroutine prints the weights of configurations
        Implicit None
        Integer :: j, k, ic, i, i1, l, n0, n1, n2, ni, ndk, m1, nspaces, nspacesg, nconfs, iter
        Real(dp) :: wsum
        Integer, Allocatable, Dimension(:,:) :: Wpsave
        Real(dp), Allocatable, Dimension(:)  :: C
        Real(dp), Allocatable, Dimension(:,:)  :: Weights, W2, Wsave
        Character(Len=512) :: strfmt, strconfig, strsp
        Character(Len=64), Allocatable, Dimension(:,:) :: strcsave
        Character(Len=3) :: strc, strq
        Integer, Dimension(0:33)  ::  nnn, nqq 
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(0:33):: lll
        Type WeightTable
            Character(Len=64), Allocatable, Dimension(:) :: strconfs, strconfsave
            Integer, Allocatable, Dimension(:) :: nconfs, nconfsave
            Real(dp), Allocatable, Dimension(:) :: wgt, wgtsave
        End Type WeightTable
        Type(WeightTable) :: wgtconfs
        Data Let/'s','p','d','f','g','h','i','k','l'/
        nconfs = 100
        strsp = ''
        If (.not. Allocated(C)) Allocate(C(Nd))
        If (.not. Allocated(Weights)) Allocate(Weights(Nc,Nlv))
        If (.not. Allocated(W2)) ALlocate(W2(Nnr, Nlv))
        If (.not. Allocated(Wsave)) Allocate(Wsave(nconfs,Nlv))
        If (.not. Allocated(Wpsave)) Allocate(Wpsave(nconfs,Nlv))
        If (.not. Allocated(strcsave)) Allocate(strcsave(nconfs,Nlv))

        If (iter <= kXIJ) Then
            Open(unit=98,status='REPLACE',file='CONF.LVL',action='WRITE')
        Else
            Open(unit=98,status='UNKNOWN',POSITION='APPEND',file='CONF.LVL',action='WRITE')
            strfmt = '(A)'
            Write(98,strfmt) ''
        End If

        ! Print line of "=" for top of table
        strfmt = '(A,I3)'
        Write(98,strfmt) 'Davidson Iteration #', iter
        Write(98,strfmt)

        ! Form matrix of weights of each configuration for each energy level
        Weights=0_dp
        W2=0_dp
        nspacesg = 0
        Wsave = 0_dp
        Do j=1,Nlv
            C(1:Nd)=ArrB(1:Nd,j)
            i=0
            Do ic=1,Nc
                Weights(ic,j)=0.d0
                ndk=Ndc(ic)
                Do k=1,ndk
                    i=i+1
                    Weights(ic,j)=Weights(ic,j)+C(i)**2
                End Do
            End Do
            k=0

            ! Calculate weights of non-relativistic configurations and store in array W2
            Do ic=1,Nnr
                Do i=1,Nrnrc(ic)
                    k=k+1
                    W2(ic,j) = W2(ic,j) + Weights(k,j)
                End Do
            End Do

            ! Save top nconfs weights and respective indices to Wsave and Wpsave
            Do i=1,nconfs
                Wsave(i,j) = maxval(W2(1:Nnr,j),1)
                Wpsave(i,j) = maxloc(W2(1:Nnr,j),1)
                W2(maxloc(W2(1:Nnr,j)),j) = 0
            End Do

            ! Save the configuration string to strcsave(k,j)
            Do k=1,nconfs
                strconfig = ''
                m1=sum(Nrnrc(1:Wpsave(k,j)))
                n1=Nc0(m1)+1
                n2=Nc0(m1)+Nvc(m1)
                Do i=n1,n2
                    i1=i-n1+1
                    ni=Nip(i)
                    l=Ll(ni)+1
                    lll(i1)=let(l)
                    nnn(i1)=Nn(ni)
                    nqq(i1)=Nq(i)
                End Do
                n=n2-n1+1
                n0=0
                Do i=2,n
                    If (nnn(i) == nnn(i-1) .and. lll(i) == lll(i-1)) Then
                        nqq(i-1) = nqq(i-1) + nqq(i)
                        nqq(i) = 0
                    End If
                End Do
                Do i=1,n
                    If (nqq(i) == 0) Cycle
                    Write(strc,'(I2)') nnn(i)
                    If (nqq(i) > 1) Then
                        Write(strq,'(I2)') nqq(i)
                    Else
                        strq = ''
                    End If
                    strconfig = Trim(AdjustL(strconfig)) // ' ' // Trim(AdjustL(strc)) // Trim(AdjustL(lll(i))) // Trim(AdjustL(strq))
                End Do
                strcsave(k,j) = strconfig
            End Do
        End Do

        ! Set weights back after assigning top values to Wsave and set maximum length of configuration string to nspacesg
        nspaces = 0
        Do j=1,Nlv
            Do k=nconfs,1,-1
                W2(Wpsave(k,j),j) = Wsave(k,j)
                If (nspacesg < (len(Trim(AdjustL(strcsave(k,j)))))) nspacesg = len(Trim(AdjustL(strcsave(k,j))))
            End Do
        End Do

        ! Write LEVELS.RES
        If (.not. Allocated(wgtconfs%strconfs)) Allocate(wgtconfs%strconfs(Nc))
        If (.not. Allocated(wgtconfs%nconfs)) Allocate(wgtconfs%nconfs(Nc))
        If (.not. Allocated(wgtconfs%strconfsave)) Allocate(wgtconfs%strconfsave(Nc))
        If (.not. Allocated(wgtconfs%nconfsave)) Allocate(wgtconfs%nconfsave(Nc))
        If (.not. Allocated(wgtconfs%wgt)) Allocate(wgtconfs%wgt(Nc))
        If (.not. Allocated(wgtconfs%wgtsave)) Allocate(wgtconfs%wgtsave(Nc))
        wgtconfs%strconfs = ''
        wgtconfs%strconfsave = ''
        wgtconfs%nconfs = 0
        wgtconfs%nconfsave = 0
        wgtconfs%wgt = 0_dp
        wgtconfs%wgtsave = 0_dp
        k=1
        Do j=1,Nlv
            Write(98,'(A)') '****************************************************'
            write(98,'(A)') ''
            strfmt = '("Level #",i3,2x,"J =",f6.3,5x," E =",f14.8)'
            Write(98, strfmt) j, Tj(j), Tk(j)+4.d0*Gj*Tj(j)*(Tj(j)+1.d0)
            Write(98,'(A)') ''
            Write(98,'(A)') ' Weight     Configuration'
            Write(98,'(A)') ''
            strfmt = '(F9.6,5X,A)'
            wsum=0_dp
            i=0
            Do While (wsum < 0.995 .and. i < nconfs)
                i=i+1
                Write(98, strfmt) Wsave(i,j), strcsave(i,j)
                If (Any(wgtconfs%strconfs == strcsave(i,j))) Then
                    Continue
                Else
                    wgtconfs%strconfs(k) = strcsave(i,j)
                    k = k + 1
                End If
                wgtconfs%nconfs(findloc(wgtconfs%strconfs, strcsave(i,j))) = wgtconfs%nconfs(findloc(wgtconfs%strconfs, strcsave(i,j))) + 1  
                wgtconfs%wgt(findloc(wgtconfs%strconfs, strcsave(i,j))) = wgtconfs%wgt(findloc(wgtconfs%strconfs, strcsave(i,j))) + Wsave(i,j)
                wsum = wsum + Wsave(i,j)
            End Do
            Write(98,'(A)') '_______'
            Write(98,'(F9.6)') wsum
            Write(98,'(A)') ''
        End Do

        k=k-1
        wgtconfs%wgt(1:k) = wgtconfs%wgt(1:k)/Nlv
        Do j=1,k
            i = maxloc(wgtconfs%wgt(1:k),1)
            wgtconfs%strconfsave(j) = wgtconfs%strconfs(i)
            wgtconfs%nconfsave(j) = wgtconfs%nconfs(i)
            wgtconfs%wgtsave(j) = maxval(wgtconfs%wgt(1:k),1)
            wgtconfs%nconfs(i) = 0
            wgtconfs%wgt(i) = 0
        End Do

        Write(98,'(A)') ''
        Write(98,'(A)')    ' Full list of configurations that appear in this file'
        Write(98,'(A)')    ' ===================================================='
        Write(98,'(A,I3)') ' Number of levels =  ', Nlv
        Write(98,'(A,I3)') ' Total number of configurations that appear in this file = ', k
        Write(98,'(A)')    ' The list is in ADD.INP format'
        Write(98,'(A)')    ' N is how many times this configuration appeared in this file'
        Write(98,'(A)')    ' W is its average (divided by the number of levels) weight'
        Write(98,'(A)') ''
        Write(98,'(A)') '     N       W         Configuration'
        
        Do j=1,k
            If (wgtconfs%wgtsave(j) > 0) Then
                Write(98,'(I6,1X,F10.5,10X,A)') wgtconfs%nconfsave(j), wgtconfs%wgtsave(j), wgtconfs%strconfsave(j)
            End If
        End Do

        Write(98,'(A)')    ' ===================================================='
        Close(98)

        Deallocate(C, Weights, W2, Wsave, Wpsave, strcsave)
        Deallocate(wgtconfs%strconfs, wgtconfs%nconfs)
        Deallocate(wgtconfs%strconfsave)
        Deallocate(wgtconfs%nconfsave)
        Deallocate(wgtconfs%wgt)
        Deallocate(wgtconfs%wgtsave)
        Return
    End Subroutine PrintWeightsDvdsn

    Subroutine WriteFinalXIJ(mype)
        Use mpi_f08
        Use davidson, Only : Prj_J
        Use formj2, Only : J_av
        Implicit None
        Integer :: i, n, ierr, mype, mpierr

        !If (K_prj == 1) Then
        !    Call Prj_J(1,Nlv,Nlv+1,1.d-8,mype)
        !End If
        If (mype == 0) Then
            open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        End If
        
        Do n=1,Nlv
            Call MPI_Bcast(ArrB(1:Nd,n), Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
            Call J_av(ArrB(1,n),Nd,Tj(n),ierr)  ! calculates expectation values for J^2
            If (mype==0) Then
                write(17) Tk(n),Tj(n),Nd,(ArrB(i,n),i=1,Nd)
            End If
        End Do

        If (mype==0) Then
            Print*, 'Final CONF.XIJ has been written'
            close(unit=17)
        End If

        Return
    End Subroutine WriteFinalXIJ

    Subroutine AllocateLSJArrays(mype)
        Use mpi_f08
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer, Intent(In) :: mype
        Integer :: mpierr

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. Allocated(Jz)) Allocate(Jz(Nst))
        If (.not. Allocated(Nh)) Allocate(Nh(Nst))
        If (.not. Allocated(Nh0)) Allocate(Nh0(Nst))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(Iarr)) Allocate(Iarr(Ne,Nd))
        If (.not. Allocated(Tk)) Allocate(Tk(Nlv))
        If (.not. Allocated(Tj)) Allocate(Tj(Nlv))
        If (.not. Allocated(xj)) Allocate(xj(Nlv))
        If (.not. Allocated(xl)) Allocate(xl(Nlv))
        If (.not. Allocated(xs)) Allocate(xs(Nlv))
        If (.not. Allocated(ArrB)) Allocate(ArrB(Nd,Nlv))

        If (mype == 0) print*, 'LSJ arrays allocated'

    End Subroutine AllocateLSJArrays

    Subroutine InitLSJ(mype)
        ! this subroutine initializes variables used for LSJ and subsequent subroutines
        ! All necessary variables are broadcasted from root to all cores
        Use mpi_f08
        Use mpi_utils
        Implicit None

        Integer :: mype, mpierr
        Integer(Kind=int64) :: count

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh0, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd,kind=int64)
        Call BroadcastI(Iarr, count, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tj, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tk, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        If (mype == 0) print*, 'LSJ arrays initialized'

        Return
    End subroutine InitLSJ

    Subroutine lsj(cc,xj,xl,xs,mype,npes)
        Use mpi_f08
        Use str_fmt, Only : startTimer, stopTimer
        Implicit None

        Real(dp), Allocatable, Dimension(:), Intent(InOut) :: xj, xl, xs
        Real(dp), Allocatable, Dimension(:,:), Intent(In) :: cc
        Integer(Kind=int64) :: s1, s2
        Integer :: mype, npes, mpierr, msg, ncsplit, nnc, nccnt, an_id, ncGrowBy, endnc, num_done
        Type(MPI_STATUS) :: status
        Integer :: n, j, sender
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002
        Character(Len=16) :: timeStr
        Logical :: moreTimers
        real(dp), dimension(:,:), allocatable :: plj, pls, p0s, pll, p0l

        ! Set moreTimers to .true. to display progress of each lsj iteration (default is .false.)
        moreTimers = .false.

        allocate(plj(nst,nst))
        allocate(pls(nst,nst))
        allocate(p0s(nst,nst))
        allocate(pll(nst,nst))
        allocate(p0l(nst,nst))
        
        If (mype == 0) write(*,*) ' calculating ME of l, s, & j...'
        do i=1,nst
            do j=1,nst
               plj(i,j)=plus_j(i,j)
               pls(i,j)=plus_s(i,j)
               p0s(i,j)=p0_s(i,j)
               pll(i,j)=plus_l(i,j)
               p0l(i,j)=p0_l(i,j)
            end do
        end do
        If (mype == 0) write(*,*) '              ... done'

        xj=0.d0
        xl=0.d0
        xs=0.d0

        Call startTimer(s1)

        ! If only 1 processor is available (serial job)
        If (npes == 1) Then
            Call calcLSJ(1,Nc,cc,xj,xl,xs,plj,pls,p0s,pll,p0l)
        ! If more than 1 processor is available 
        Else
            n=0 
            ncGrowBy = 1

            If (mype == 0) Then
                ! Distribute a portion of the workload of size ncGrowBy to each worker process
                Do an_id = 1, npes - 1
                   nnc = ncGrowBy*an_id + 1
                   Call MPI_SEND( nnc, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                End Do

                Call calcLSJ(1,ncGrowBy,cc,xj,xl,xs,plj,pls,p0s,pll,p0l)
                num_done = 0
                ncsplit = Nc/10
                nccnt = ncsplit
                j=9

                Do 
                    Call MPI_RECV(msg, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    sender = status%MPI_SOURCE

                    If (nnc + ncGrowBy <= Nc) Then
                        nnc = nnc + ncGrowBy
                        Call MPI_SEND( nnc, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                    Else
                        msg = -1
                        Call MPI_SEND( msg, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If

                    If (nnc == nccnt .and. nnc /= ncsplit*10) Then
                        Call stopTimer(s1, timeStr)
                        If (moreTimers == .true.) Write(*,'(2X,A,1X,I3,A)'), 'lsj:', (10-j)*10, '% done in '// trim(timeStr)
                        j=j-1
                        nccnt = nccnt + ncsplit
                    End If

                    If (num_done == npes-1) Then
                        Call stopTimer(s1, timeStr)
                        If (moreTimers == .true.) Write(*,'(2X,A,1X,I3,A)'), 'lsj:', (10-j)*10, '% done in '// trim(timeStr)
                        Exit
                    End If
                End Do
            Else
                Do 
                    Call MPI_RECV ( nnc, 1 , MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    Call startTimer(s2)
                    If (nnc == -1) Then
                        Exit
                    Else
                        If (Nc - nnc < ncGrowBy) Then
                            endnc = Nc
                        Else
                            endnc = nnc+ncGrowBy-1
                        End If

                        Call calcLSJ(nnc,endnc,cc,xj,xl,xs,plj,pls,p0s,pll,p0l)

                        Call MPI_SEND( msg, 1, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                        Call stopTimer(s2,timeStr)
                    End If
                End Do
            End If
            Call MPI_AllReduce(MPI_IN_PLACE, xj, Nlv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(MPI_IN_PLACE, xl, Nlv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(MPI_IN_PLACE, xs, Nlv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End If

        Deallocate(plj, pls, p0s, pll, p0l)        

        Xj=0.5d0*(dsqrt(Xj+1)-1)
        Xl=0.5d0*(dsqrt(Xl+1)-1)
        Xs=0.5d0*(dsqrt(Xs+1)-1)
            
        Return
    End Subroutine lsj

    Subroutine calcLSJ(startnc,endnc,cc,xj,xl,xs,plj,pls,p0s,pll,p0l)
        Use determinants, Only : Gdet, CompC, CompNRC
        Implicit None
        Integer, Intent(In) :: startnc, endnc
        Real(dp), Allocatable, Dimension(:), Intent(InOut) :: xj, xl, xs
        Real(dp), Allocatable, Dimension(:,:), Intent(In) :: cc
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Integer :: ic1, ic2, n1, ndn, n, k, k1n, ndk, k1, icomp
        Real(dp) :: ckn, tj, tl, ts
        real(dp), dimension(:,:), allocatable :: plj, pls, p0s, pll, p0l


        if (.not. allocated(idet1)) Allocate(idet1(Ne))
        if (.not. allocated(idet2)) Allocate(idet2(Ne))
        if (.not. allocated(iconf1)) Allocate(iconf1(Ne))
        if (.not. allocated(iconf2)) Allocate(iconf2(Ne))
        Do ic1=startnc, endnc
            ndn=Ndc(ic1)
            n=sum(Ndc(1:ic1-1))
            Do n1=1,ndn
                n=n+1
                call Gdet(n,idet1)
                k=n-1
                Do ic2=ic1,Nc
                    ndk=Ndc(ic2)
                    k1n=1
                    If (ic2.EQ.ic1) k1n=n1
                    call Gdet(k+1,idet2)
                    call CompNRC(idet1,idet2,icomp)
                    if (icomp.GT.0) then
                        k=k+ndk
                        Cycle
                    end if
                    call CompC(idet1,idet2,icomp)
                    if (icomp.GT.2) then
                        k=k+ndk
                        Cycle
                    end if
                    Do k1=k1n,ndk
                        k=k+1
                        call Gdet(k,idet2)
                        call lsj_det(idet1,idet2,tj,tl,ts,plj,pls,p0s,pll,p0l)
                        do m=1,Nlv
                            ckn=cc(n,m)*cc(k,m)
                            if (n.ne.k) ckn=2*ckn
                            xj(m)=xj(m)+ckn*tj
                            xl(m)=xl(m)+ckn*tl
                            xs(m)=xs(m)+ckn*ts
                        end do
                    End Do
                End Do
            End Do
        End Do
    End Subroutine calcLSJ

    Subroutine lsj_det(idet1,idet2,tj,tl,ts,plj,pls,p0s,pll,p0l)
        Use determinants, Only : Rspq
        Implicit None
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Real(dp) :: tj, tl, ts
        Integer :: ic, id, is, nf, i1, i2, j1, j2, iq, ja, la, ma, jq0, jq
        Integer :: ia, ib, ka, kb, kc, kd, na, nb, nc, nd
        Real(dp), Allocatable, Dimension(:,:) :: plj, pls, p0s, pll, p0l

        tj=0.d0
        tl=0.d0
        ts=0.d0
        call Rspq(idet1,idet2,is,nf,i1,j1,i2,j2)
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!        Determinants are equal
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.0) then
          tj=mj*mj
          Do iq=1,Ne
            ia=idet1(iq)
            na=Nh(ia)
            ja=Jj(na)
            la=Ll(na)
            ma=Jz(ia)
            tj=tj+ja*(ja+2)-ma**2
            ts=ts+0.75d0
            tl=tl+la*(la+1)
            jq0=iq+1
            If (jq0.LE.Ne) then
              Do jq=jq0,Ne
                ib=idet1(jq)
                tj=tj-plj(ia,ib)**2-plj(ib,ia)**2
                ts=ts+2*(p0s(ia,ia)*p0s(ib,ib)-p0s(ia,ib)**2)
                ts=ts-2*(pls(ia,ib)**2+pls(ib,ia)**2)
                tl=tl+2*(p0l(ia,ia)*p0l(ib,ib)-p0l(ia,ib)**2)
                tl=tl-2*(pll(ia,ib)**2+pll(ib,ia)**2)
              End Do
            End If
          End Do
        End If
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!       One function differs in the Determinants.
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.1) then 
            ia=i2
            ib=j2
            ka=Nh(ia)
            kb=Nh(ib)
            na=Nn(ka)
            nb=Nn(kb)
            if (na.EQ.nb) then
                do iq=1,ne
                    ic=idet1(iq)
                    if (ic.eq.i2) Cycle
                    ts=ts+2*(p0s(ia,ib)*p0s(ic,ic)-p0s(ia,ic)*p0s(ib,ic))
                    ts=ts-2*(pls(ia,ic)*pls(ib,ic)+pls(ic,ia)*pls(ic,ib))
                    tl=tl+2*(p0l(ia,ib)*p0l(ic,ic)-p0l(ia,ic)*p0l(ib,ic))
                    tl=tl-2*(pll(ia,ic)*pll(ib,ic)+pll(ic,ia)*pll(ic,ib))
                end do
                ts=ts*is
                tl=tl*is
            end if
        End If
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!        Determinants differ by two functions
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.2) then
            ia=i1
            ic=j1
            ib=i2
            id=j2
            ka=Nh(ia)
            kc=Nh(ic)
            kb=Nh(ib)
            kd=Nh(id)
            na=Nn(ka)
            nc=Nn(kc)
            nb=Nn(kb)
            nd=Nn(kd)
            if (na.EQ.nc.AND.nb.EQ.nd.OR.na.EQ.nd.AND.nb.EQ.nc) then
                 tj=plj(ia,ic)*plj(id,ib)+plj(ic,ia)*plj(ib,id)- &
                   plj(ia,id)*plj(ic,ib)-plj(id,ia)*plj(ib,ic)
                 ts=ts+2*(p0s(ia,ic)*p0s(id,ib)-p0s(ia,id)*p0s(ic,ib))
                 ts=ts+2*(pls(ia,ic)*pls(id,ib)+pls(ic,ia)*pls(ib,id)- &
                    pls(ia,id)*pls(ic,ib)-pls(id,ia)*pls(ib,ic))
                 tl=tl+2*(p0l(ia,ic)*p0l(id,ib)-p0l(ia,id)*p0l(ic,ib))
                 tl=tl+2*(pll(ia,ic)*pll(id,ib)+pll(ic,ia)*pll(ib,id)- &
                    pll(ia,id)*pll(ic,ib)-pll(id,ia)*pll(ib,ic))
                 tj=tj*is
                 ts=ts*is
                 tl=tl*is
            end if
        End If
        ts=ts*4
        tl=tl*4

        Return
    End Subroutine lsj_det

    Real(dp) Function plus_s(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)+2) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        msa=1
        msb=-1
        mja=jz(ia)
        mjb=jz(ib)
        mla=(mja-msa)/2
        mlb=mla
        If (iabs(mla).gt.la) goto 1000
        ta=0.d0
        If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
        If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
        tb=0.d0
        If (jb.eq.2*lb+1) tb= dsqrt((jb-mjb)/(2.d0*jb))
        If (jb.eq.2*lb-1) tb= dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
        t=-dsqrt(0.5d0)*ta*tb

1000    plus_s=t
        Return
    End Function plus_s

    Real(dp) Function plus_l(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)+2) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        mja=jz(ia)
        mjb=jz(ib)
        Do msa=-1,1,2
            msb=msa
            mla=(mja-msa)/2
            mlb=(mjb-msb)/2
            If (iabs(mla).gt.la) Cycle
            If (iabs(mlb).gt.lb) Cycle
            ta=0.d0
            tb=0.d0
            If (msa.eq.1) then
                If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            If (msa.eq.-1) then
                If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t-dsqrt(0.5d0*(lb*(lb+1.d0)-mlb*(mlb+1.d0)))*ta*tb
        End Do
1000    plus_l=t
        Return
    End Function plus_l

    Real(dp) Function p0_s(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
        ! We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
        ! but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        Do msa=-1,1,2
            msb=msa
            mja=jz(ia)
            mjb=jz(ib)
            mla=(mja-msa)/2
            mlb=mla
            If (iabs(mla).gt.la) Cycle
                ta=0.d0
                tb=0.d0
                If (msa.eq.1) Then
                    If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                    If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                    If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                    If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
                End If
                If (msa.eq.-1) then
                    If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                    If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                    If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                    If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t+0.5d0*msa*(ta*tb)
        End Do

1000    p0_s=t
        Return
    End Function p0_s

    Real(dp) Function p0_l(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.

        t=0.
        If (jz(ia).ne.jz(ib)) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        mja=jz(ia)
        mjb=jz(ib)
        Do msa=-1,1,2
            msb=msa
            mla=(mja-msa)/2
            mlb=mla
            If (iabs(mla).gt.la) Cycle
            ta=0.d0
            tb=0.d0
            If (msa.eq.1) then
                If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            If (msa.eq.-1) then
                If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t+mla*(ta*tb)
        End Do

1000    p0_l=t
        Return
    End Function p0_l

    Real(dp) Function Plus_j(ia,ib)
        Implicit None    
        Integer :: ia, ib, na, nb, ma, mb, ja
        Real(dp) :: t     

        t=0.d0
        na=Nh(ia)
        nb=Nh(ib)
        if (na.EQ.nb) then
           ma=Jz(ia)
           mb=Jz(ib)
           if (ma.EQ.mb+2) then
              ja=Jj(na)
              t=dsqrt(ja*(ja+2)-ma*mb+0.d0)
           end if
        end if
        Plus_j=t
        Return
    End Function Plus_j

    Subroutine PrintEnergies
        ! This subroutine prints eigenvalues in increasing order
        Implicit None

        Integer :: j, ist
        Real(kind=type_real) :: dt, del
        Character(Len=1), Dimension(10) :: stecp*7
        Character(Len=1), Dimension(4)  :: strsms*6
        Character(Len=1), Dimension(3)  :: strms*3
        Character(Len=512) :: strfmt
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/

        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7

        ! Print line of "=" for top of table
        strfmt = '(4X,63("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! If pure CI, print Nc, Nd, Gj
        If (Ksig*Kdsig == 0) Then
            strfmt = '(4X,"Energy levels (",A7," Nc=",I7," Nd=",I9,"); Gj =",F7.4, &
                        /4X,"N",6X,"JTOT",12X,"EV",16X,"ET",9X,"DEL(CM**-1)")'
            Write( 6,strfmt) stecp(ist),Nc,Nd,Gj
            Write(11,strfmt) stecp(ist),Nc,Nd,Gj
        ! If CI+all-order/CI+MBPT, print E_0, Kexn, Nc, Nd, Gj
        Else
            strfmt = '(4X,"Energy levels ",A7,", Sigma(E =",F10.4,") extrapolation var.", &
                    I2,/4X,"(Nc=",I7," Nd=",I9,"); Gj =",F7.4,/4X,"N",6X,"JTOT",12X, &
                    "EV",16X,"ET",9X,"DEL(CM**-1)")'
            Write( 6,strfmt) stecp(ist),E_0,Kexn,Nc,Nd,Gj
            Write(11,strfmt) stecp(ist),E_0,Kexn,Nc,Nd,Gj
        End If

        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                strfmt = '(4X,"Volume shift: dR_N/R_N=",F9.5," Rnuc=",F10.7)'
                Write( *,strfmt) C_is,Rnuc
                Write(11,strfmt) C_is,Rnuc
            Else
                strfmt = '(4X,A3,":",E9.2,"*(P_i Dot P_k) ",A6," Lower component key =",I2)'
                Write( *,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
                Write(11,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
            End If
        End If

        ! Print line of "-" to close header
        strfmt = '(4X,63("-"))'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! For each energy level, print the following:
        ! j = index for energy level
        ! Tj(j) = total angular momentum
        ! Tk(j) = valence energy of jth level
        ! DT = Tk(j)-Ecore = total energy
        ! DEL = Tk(1)-Tk(j) = energy difference between jth level and 1st level in cm**-1
        Do j=1,Nlv
            Tk(j)=Tk(j)+4.d0*Gj*Tj(j)*(Tj(j)+1.d0)
            DT=Tk(j)-Ecore
            DEL=(Tk(1)-Tk(j))*2*DPRy
            Write( 6,'(3X,I3,F14.9,F14.8,F18.6,F15.2)') j,Tj(j),Tk(j),DT,DEL
            Write(11,'(3X,I3,F14.9,F14.8,F18.6,F15.2)') j,Tj(j),Tk(j),DT,DEL
        End Do

        ! Print line of "=" for bottom of table
        strfmt = '(4X,63("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

    End Subroutine PrintEnergies

    Subroutine PrintWeights
        ! This subroutine prints the weights of configurations
        Implicit None
        Integer :: j, k, j1, j2, j3, ic, i, ii, i1, l, n0, n1, n2, ni, nk, ndk, m1, nspaces, nspacesg, nconfs, cnt, nspacesterm, maxlenconfig, num_blanks, num_blanks2
        Real(dp) :: wsum, gfactor, ax_crit, j_crit
        Integer, Allocatable, Dimension(:,:) :: Wpsave
        Real(dp), Allocatable, Dimension(:)  :: C
        Real(dp), Allocatable, Dimension(:,:)  :: W, W2, Wsave
        Character(Len=1), Dimension(11) :: st1, st2
        Character(Len=512) :: strfmt, strfmt2, strconfig, strconfig2, strsp
        Character(Len=64), Allocatable, Dimension(:,:) :: strcsave
        Character(Len=3) :: strc, strq
        Character(Len=6) :: strterm, strconverged
        Integer, Dimension(33)  ::  nnn, nqq 
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(33):: lll
        Character(Len=5) :: strgf, strconfadd2
        Character(Len=5), Allocatable, Dimension(:) :: strconfadd
        Logical :: has_secondary_config
        Logical, Allocatable, Dimension(:) :: converged
        Type WeightTable
            Character(Len=64), Allocatable, Dimension(:) :: strconfs, strconfsave
            Integer, Allocatable, Dimension(:) :: nconfs, nconfsave
            Real(dp), Allocatable, Dimension(:) :: wgt, wgtsave
        End Type WeightTable
        Type(WeightTable) :: wgtconfs
        data st1/11*'='/, st2/11*'-'/
        Data Let/'s','p','d','f','g','h','i','k','l'/

        nconfs = 100
        strsp = ''
        If (.not. Allocated(C)) Allocate(C(Nd))
        If (.not. Allocated(W)) Allocate(W(Nc,Nlv))
        If (.not. Allocated(W2)) ALlocate(W2(Nnr, Nlv))
        If (.not. Allocated(Wsave)) Allocate(Wsave(nconfs,Nlv))
        If (.not. Allocated(Wpsave)) Allocate(Wpsave(nconfs,Nlv))
        If (.not. Allocated(strcsave)) Allocate(strcsave(nconfs,Nlv))
        If (.not. Allocated(converged)) Allocate(converged(Nlv))

        If (kWeights == 1) Then
            Open(88,file='CONF.WGT',status='UNKNOWN')
            strfmt = '(I8,I6,F11.7)'
        End If

        Open(99,file='FINAL.RES',status='UNKNOWN')
        Open(98,file='LEVELS.RES',status='UNKNOWN')
        Open(97,file='CONFSTR.RES',status='UNKNOWN')

        ! Form array of booleans of converged levels
        converged = .False.
        Do i=1,Nlv
            If (Iconverge(i) == 1) converged(i) = .True.
        End Do

        ! Set levels with small enough admixtures to also be converged
        ax_crit = 1e-4
        Do i=1,Nlv
            If (ax_array(i) <= ax_crit) converged(i) = .True.
        End Do 

        ! TODO - set levels where J uncertainty is above some threshold to be not converged
        j_crit = 1e-4
        Do i=1,Nlv
            If (Tj(i)-Nint(Tj(i)) > j_crit) converged(i) = .False.
        End Do 

        ! Form matrix of weights of each configuration for each energy level
        W2=0_dp
        nspacesg = 0
        nspacesterm = 0
        If (mod(Ne,2) == 1) nspacesterm = 3
        Wsave = 0_dp
        Do j=1,Nlv
            If (kWeights == 1) Then
                Write(88,'(A,I3)') 'Level #', j
                Write(88,'(A25)') '========================='
                Write(88,'(A25)') '      ID    IC     W     '
            End If
            C(1:Nd)=ArrB(1:Nd,j)
            i=0
            Do ic=1,Nc
                W(ic,j)=0.d0
                ndk=Ndc(ic)
                Do k=1,ndk
                    i=i+1
                    W(ic,j)=W(ic,j)+C(i)**2
                    If (kWeights == 1) Write(88,strfmt) i, ic, C(i)**2
                End Do
            End Do
            If (kWeights == 1) Write(88,'(A25)') '========================='
            k=0

            ! Calculate weights of non-relativistic configurations and store in array W2
            Do ic=1,Nnr
                Do i=1,Nrnrc(ic)
                    k=k+1
                    W2(ic,j) = W2(ic,j) + W(k,j)
                End Do
            End Do

            ! Save top nconfs weights and respective indices to Wsave and Wpsave
            Do i=1,nconfs
                Wsave(i,j) = maxval(W2(1:Nnr,j),1)
                Wpsave(i,j) = maxloc(W2(1:Nnr,j),1)
                W2(maxloc(W2(1:Nnr,j)),j) = 0
            End Do

            ! Save the configuration string to strcsave(k,j)
            Do k=1,nconfs
                strconfig = ''
                m1=sum(Nrnrc(1:Wpsave(k,j)))
                n1=Nc0(m1)+1
                n2=Nc0(m1)+Nvc(m1)
                Do i=n1,n2
                    i1=i-n1+1
                    ni=Nip(i)
                    l=Ll(ni)+1
                    lll(i1)=let(l)
                    nnn(i1)=Nn(ni)
                    nqq(i1)=Nq(i)
                End Do
                n=n2-n1+1
                n0=0
                Do i=2,n
                    If (nnn(i) == nnn(i-1) .and. lll(i) == lll(i-1)) Then
                        nqq(i-1) = nqq(i-1) + nqq(i)
                        nqq(i) = 0
                    End If
                End Do
                Do i=1,n
                    If (nqq(i) == 0) Cycle
                    Write(strc,'(I2)') nnn(i)
                    If (nqq(i) > 1) Then
                        Write(strq,'(I2)') nqq(i)
                    Else
                        strq = ''
                    End If
                    If (strconfig == '') Then
                        strconfig = Trim(AdjustL(strc)) // Trim(AdjustL(lll(i))) // Trim(AdjustL(strq))
                    Else
                        strconfig = Trim(AdjustL(strconfig)) // ' ' // Trim(AdjustL(strc)) // Trim(AdjustL(lll(i))) // Trim(AdjustL(strq))
                    End If
                End Do
                strcsave(k,j) = strconfig
            End Do
        End Do

        ! Set weights back after assigning top values to Wsave and set maximum length of configuration string to nspacesg
        nspaces = 0
        Do j=1,Nlv
            Do k=nconfs,1,-1
                W2(Wpsave(k,j),j) = Wsave(k,j)
                If (nspacesg < (len(Trim(AdjustL(strcsave(k,j)))))) nspacesg = len(Trim(AdjustL(strcsave(k,j))))
            End Do
        End Do

        ! Check if a level has a secondary configuration with a weight over some threshold
        If (any(Wsave(1,:) < 0.7)) has_secondary_config = .True.

        ! Find longest configuration length
        maxlenconfig = 0
        Do i=1,Nlv
            If (len_trim(strcsave(1,i)) > maxlenconfig) Then
                maxlenconfig = len_trim(strcsave(1,i))
            End If
        End Do

        ! Write table of configurations, L, S, J, energies, and weights of top 2 configurations to FINAL.RES
        Do j=1,Nlv
            ! Calculate g-factors if including L, S, J
            If (KLSJ == 1) Then
                If (Nint(Tj(j)) == 0) Then
                    Write(strgf,'(A)') '-----'
                Else
                    gfactor = g_factor(Xl(j),Xs(j),Tj(j))
                    Write(strgf,'(F5.3)') gfactor
                End If
            End If

            ! Writes main configurations to CONFSTR.RES
            Write(97,'(A)') Trim(AdjustL(strcsave(1,j)))

            ! If LSJ is calculated, also include terms in CONFSTR.RES
            If (KLSJ == 1) Then
                strterm = term(Xl(j), Xs(j), Tj(j))
                Write(97,'(A)') Trim(AdjustL(strterm))
            End If

            ! Set string for convergence
            if (converged(j)) Then
                strconverged = '  True'
            Else
                strconverged = ' False'
            End If

            ! If L, S, J is needed
            If (KLSJ == 1) Then
                ! Write column names if first iteration
                If (j == 1) Then
                    num_blanks = max(0, maxlenconfig-3)
                    Write(99, '(A)') '  n' // '  ' // repeat(' ', num_blanks) // 'conf  term     E_n (a.u.)   DEL (cm^-1)     S     L     J     gf     conf%  converged'// repeat(' ', num_blanks) // ' conf2  conf2%'
                End If
                ! If main configuration has weight of less than 0.7, we have to include a secondary configuration
                num_blanks = max(0, maxlenconfig-len_trim(strcsave(1,j))+1)
                If (Wsave(1,j) < 0.7) Then
                    num_blanks2 = max(0, maxlenconfig-len_trim(strcsave(2,j))+1)
                    strfmt = '(I3,3X,A,A,1X,f14.8,f14.1,2X,f4.2,2x,f4.2,2x,f4.2,2x,A,4x,f4.1,"%",5X,A,2X,A,3X,f5.1,"%")'
                    Write(99,strfmt) j, repeat(' ', num_blanks) // Trim(AdjustL(strcsave(1,j))), AdjustR(strterm), Tk(j), (Tk(1)-Tk(j))*2*DPRy, Xs(j), Xl(j), Tj(j), strgf, Wsave(1,j)*100, strconverged, repeat(' ', num_blanks2) // Trim(AdjustL(strcsave(2,j))), Wsave(2,j)*100
                ! Else we include only the main configuration
                Else
                    strfmt = '(I3,3X,A,A,1X,f14.8,f14.1,2X,f4.2,2x,f4.2,2x,f4.2,2x,A,4x,f5.1,"%",5X,A)'
                    Write(99,strfmt) j, repeat(' ', num_blanks) // Trim(AdjustL(strcsave(1,j))), AdjustR(strterm), Tk(j), (Tk(1)-Tk(j))*2*DPRy, Xs(j), Xl(j), Tj(j), strgf, maxval(W2(1:Nnr,j))*100, strconverged
                End If
            ! If L, S, J is not needed
            Else
                ! Write column names if first iteration
                If (j == 1) Then
                    num_blanks = max(0, maxlenconfig-3)
                    Write(99, '(A)') '  n ' // '  ' // repeat(' ', num_blanks) // 'conf       J         E_n (a.u.)   DEL (cm^-1)    conf%  converged'// repeat(' ', num_blanks) // ' conf2  conf2%'
                End If

                ! If main configuration has weight of less than 0.7, we have to include a secondary configuration
                num_blanks = max(0, maxlenconfig-len_trim(strcsave(1,j))+1)
                If (Wsave(1,j) < 0.7) Then
                    num_blanks2 = max(0, maxlenconfig-len_trim(strcsave(2,j))+1)
                    strfmt = '(I3,4X,A,4X,f4.2,5x,f14.8,f14.1,3X,f4.1,"%",,5X,A,2X,A,3X,f5.1,"%")'
                    Write(99,strfmt) j, repeat(' ', num_blanks) // Trim(AdjustL(strcsave(1,j))), Tj(j), Tk(j), (Tk(1)-Tk(j))*2*DPRy, Wsave(1,j)*100, strconverged, repeat(' ', num_blanks2) // Trim(AdjustL(strcsave(2,j))), Wsave(2,j)*100
                ! Else we include only the main configuration
                Else
                    strfmt = '(I3,4X,A,4X,f4.2,5x,f14.8,f14.1,3X,f5.1,"%",5X,A)'
                    Write(99,strfmt) j, repeat(' ', num_blanks) // Trim(AdjustL(strcsave(1,j))), Tj(j), Tk(j), (Tk(1)-Tk(j))*2*DPRy, maxval(W2(1:Nnr,j))*100, strconverged
                End If
            End If 
        End Do

        ! Write LEVELS.RES
        If (.not. Allocated(strconfadd)) Allocate(strconfadd(nconfs))
        If (.not. Allocated(wgtconfs%strconfs)) Allocate(wgtconfs%strconfs(Nc))
        If (.not. Allocated(wgtconfs%nconfs)) Allocate(wgtconfs%nconfs(Nc))
        If (.not. Allocated(wgtconfs%strconfsave)) Allocate(wgtconfs%strconfsave(Nc))
        If (.not. Allocated(wgtconfs%nconfsave)) Allocate(wgtconfs%nconfsave(Nc))
        If (.not. Allocated(wgtconfs%wgt)) Allocate(wgtconfs%wgt(Nc))
        If (.not. Allocated(wgtconfs%wgtsave)) Allocate(wgtconfs%wgtsave(Nc))
        wgtconfs%strconfs = ''
        wgtconfs%nconfs = 0
        wgtconfs%wgt = 0_dp
        k=1
        Do j=1,Nlv
            Write(98,'(A)') '****************************************************'
            write(98,'(A)') ''
            strfmt = '("Level #",i3,2x,"J =",f6.3,5x," E =",f14.8)'
            Write(98, strfmt) j, Tj(j), Tk(j)
            Write(98,'(A)') ''
            Write(98,'(A)') ' Weight     Configuration'
            Write(98,'(A)') ''
            strfmt = '(F9.6,5X,A)'
            wsum=0_dp
            i=0
            Do While (wsum < 0.9999 .and. i < nconfs)
                i=i+1
                Write(98, strfmt) Wsave(i,j), strcsave(i,j)
                If (Any(wgtconfs%strconfs == strcsave(i,j))) Then
                    Continue
                Else
                    wgtconfs%strconfs(k) = strcsave(i,j)
                    k = k + 1
                End If
                wgtconfs%nconfs(findloc(wgtconfs%strconfs, strcsave(i,j))) = wgtconfs%nconfs(findloc(wgtconfs%strconfs, strcsave(i,j))) + 1  
                wgtconfs%wgt(findloc(wgtconfs%strconfs, strcsave(i,j))) = wgtconfs%wgt(findloc(wgtconfs%strconfs, strcsave(i,j))) + Wsave(i,j)
                wsum = wsum + Wsave(i,j)
            End Do
            Write(98,'(A)') '_______'
            Write(98,'(F9.6)') wsum
            Write(98,'(A)') ''
        End Do

        k=k-1
        wgtconfs%wgt(1:k) = wgtconfs%wgt(1:k)/Nlv
        Do j=1,k
            i = maxloc(wgtconfs%wgt(1:k),1)
            wgtconfs%strconfsave(j) = wgtconfs%strconfs(i)
            wgtconfs%nconfsave(j) = wgtconfs%nconfs(i)
            wgtconfs%wgtsave(j) = maxval(wgtconfs%wgt(1:k),1)
            wgtconfs%nconfs(i) = 0
            wgtconfs%wgt(i) = 0
        End Do

        Write(98,'(A)') ''
        Write(98,'(A)')    ' Full list of configurations that appear in this file'
        Write(98,'(A)')    ' ===================================================='
        Write(98,'(A,I3)') ' Number of levels =  ', Nlv
        Write(98,'(A,I3)') ' Total number of configurations that appear in this file = ', k
        Write(98,'(A)')    ' The list is in ADD.INP format'
        Write(98,'(A)')    ' N is how many times this configuration appeared in this file'
        Write(98,'(A)')    ' W is its average (divided by the number of levels) weight'
        Write(98,'(A)') ''
        Write(98,'(A)') '     N       W         Configuration'
        
        Do j=1,k
            If (wgtconfs%wgtsave(j) > 0) Then
                Write(98,'(I6,1X,F11.6,5X,A)') wgtconfs%nconfsave(j), wgtconfs%wgtsave(j), wgtconfs%strconfsave(j)
            End If
        End Do
        
        Write(98,'(A)')    ' ===================================================='
        Write(98,'(A)') ''
        Write(98,'(A)')    ' Full list of configurations in ADD.INP format'
        Write(98,'(A)') ''

        Do j=1,k
            cnt = 1
            strconfadd2 = ''
            If (wgtconfs%wgtsave(j) > 0) Then
                Do i=1,Len(Trim(AdjustL(wgtconfs%strconfsave(j))))+1
                    strconfig2 = Trim(AdjustL(wgtconfs%strconfsave(j)))
                    strconfadd2 = Trim(strconfadd2) // strconfig2(i:i)
                    If (strconfig2(i:i) == '') Then
                        If (any(strconfig2(i-1:i-1) == ["s","p","d","f","g","h","i"])) Then
                            strconfadd2 = Trim(strconfadd2) // "1"
                        End If
                        If (any(strconfadd2(2:2) == ["s","p","d","f","g","h","i"])) Then
                            strconfadd2 = " " // Trim(strconfadd2)
                        End If
                        strconfadd(cnt) = strconfadd2
                        strconfadd2 = ''
                        cnt = cnt + 1
                    End If
                End Do
                strconfadd(cnt) = strconfadd2
                If (cnt > 7) Then
                    Write(98,'("L:  ",6(A,1X), /"L:  ",6(A,1X))') strconfadd(1:cnt-1)
                Else
                    Write(98,'("L:  ",6(A,1X))') strconfadd(1:cnt-1)
                End If
            End If
        End Do

        Close(97)
        Close(98)
        Close(99)
        If (kWeights == 1) Close(88)

        ! Print the weights of each configuration
        n=(Nlv-1)/5+1
        j2=0
        Do k=1,n
            nk=5
            If (k == n) nk=Nlv-(n-1)*5
            j1=j2+1
            j2=j2+nk
            j3=j2+1
            strfmt = '(3X,66A1)'
            Write(11,strfmt) (st1,i=j1,j3)
            strfmt2 = '(15X,5(3X,I2,6X))'
            Write(11,strfmt2) (i,i=j1,j2)
            Write(11,strfmt) (st2,i=j1,j3)
            strfmt2 = '(4X,"ICONF",4X,5F11.5)'
            Write(11,strfmt2) (Tk(i),i=j1,j2)
            Write(11,strfmt) (st2,i=j1,j3)
            strfmt2 = '(2X,I7,6X,5(E9.3,2X))'
            Do ic=1,Nc
                Write(11,strfmt2) ic,(W(ic,j),j=j1,j2)
            End Do
            Write(11,strfmt) (st1,i=j1,j3)
        End Do
        Deallocate(C, W, W2, Wsave, Wpsave, strcsave, Ndc, Tk, Tj, ArrB)
        Close(11)
        
    End Subroutine PrintWeights

    Character(Len=6) Function term(L, S, J)
        ! This function returns the term, given L, S and J
        Implicit None
        Real(dp), Intent(In) :: L, S, J
        Integer :: is, il, ij
        Character(Len=1) :: strl

        is = Nint(2*S+1)
        il = Nint(L)

        Select Case (il)
        Case(0)
            strl = 'S' 
        Case(1)
            strl = 'P' 
        Case(2)
            strl = 'D' 
        Case(3)
            strl = 'F' 
        Case(4)
            strl = 'G' 
        Case(5)
            strl = 'H'
        Case(6)
            strl = 'I'
        Case(7)
            strl = 'K'  
        Case(8)
            strl = 'L'
        Case(9)
            strl = 'M'     
        Case Default
            strl = '>'    
        End Select

        If (mod(Ne,2) == 0) Then
          ij = Nint(J)
          Write(term, '(I1,A,I1)') is, strl, ij
        Else
          ij = Nint(2*J)  
          If (ij < 10) Then
             Write(term, '(I1,A,I1,A2)') is, strl, ij,'/2'
          Else
             Write(term, '(I1,A,I2,A2)') is, strl, ij,'/2'
          EndIf 
        End If  

    End Function term

    Real(dp) Function g_factor(L, S, J)
        ! This function returns the g-factor, given L, S and J
        Implicit None
        Real(dp), Intent(In) :: L, S, J

        g_factor = 1 + (J*(J+1)-L*(L+1)+S*(S+1))/(2*J*(J+1))

    End Function g_factor

End Program pconf
