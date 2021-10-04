Program conf 
    ! ======== original version by I.I.Tupitsin =======
    ! ======== parallel version by C. Cheung ==========
    ! latest version of parallel code can be found at
    ! >>>   https://github.com/ccheung93/pCI  <<<
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
    Use conf_variables
    use mpi
    use davidson, only : Prj_J
    use determinants, only : Wdet, Dinit, Jterm
    Use integrals, only : Rint
    use formj2, only : FormJ, J_av
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime

    Implicit None

    Integer   :: n, k, i, j, ierr, mype, npes, mpierr
    Integer(kind=int64) :: start_time
    Real :: total_time
    Real(dp)  :: t
    Character(Len=1024) :: strFromEnv
    Character(Len=255)  :: eValue, strfmt
    Character(Len=16)   :: memStr, timeStr

    ! Initialize MPI
    Call MPI_Init(mpierr)
    ! Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    ! Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    Call startTimer(start_time)
    
    ! Give ConfFilePaths a chance to decide what filenames/paths to use:
    !Call ConfFileInit()

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    Call Get_Environment_Variable("CONF_MAX_BYTES_PER_CPU", eValue)
    read(eValue,'(I12)') memTotalPerCPU  

    ! Initialization subroutines
    ! master processor will read input files and broadcast parameters and arrays to all processors
    If (mype == 0) Then
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
    Call AllocateFormHArrays(mype,npes)
    
    ! Formation of Hamiltonian matrix
    Call FormH(npes,mype)

    ! Formation of J^2 matrix
    Call FormJ(mype, npes)   
    
    ! Deallocate arrays that are no longer used in FormH
    Call DeAllocateFormHArrays(mype,npes)

    ! Allocate arrays that are used in Davidson procedure
    Call AllocateDvdsnArrays(mype,npes)
    
    ! Davidson diagonalization
    Call Diag4(mype,npes) 
    
    ! Write final eigenvectors to file CONF.XIJ
    Call WriteFinalXIJ(mype,npes)

    ! Print table of final results and total computation time
    If (mype==0) Then
        Call PrintEnergies

        Call stopTimer(start_time, timeStr)
        write(*,'(2X,A)'), 'TIMING >>> Total computation time of conf was '// trim(timeStr)
    End If

    ! Finalize MPI 
    Call MPI_Finalize(mpierr)

Contains

    Subroutine Input
        ! This subroutine reads in parameters and configurations from CONF.INP
        Use conf_init, only : inpstr, ReadConfInp, ReadConfigurations
        Implicit None
        Character(Len=32) :: strfmt

        ! Write name of program
        open(unit=11,status='UNKNOWN',file='CONF.RES')
        strfmt = '(4X,"Program conf v0.3.20")'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! Read input parameters from file CONF.INP
        Call ReadConfInp

        If (dabs(C_is) < 1.d-6) K_is=0
        If (K_is == 0) C_is=0.d0

        ! Read job parameters from file c.in
        ! Kl = 0 - new computation
        ! Kl = 1 - continuing computation with completed CONF.HIJ and CONF.JJJ files
        ! Kl = 2 - new computation with MBPT
        ! Kl = 3 - extending computation with new configurations (not implemented yet)
        Open(unit=99,file='c.in',status='OLD')
        Read(99,*) Kl, Ksig, Kdsig, Kw
        Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
        If (K_is == 2.OR.K_is == 4) Then
            Read(99,*) K_sms
            Write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): ', K_sms
            If ((K_sms-1)*(K_sms-2)*(K_sms-3) /= 0) Stop
        End If

        ! Kw determines whether CONF.HIJ will be written or not
        ! Kw=0 - CONF.HIJ will not be written
        ! Kw=1 - CONF.HIJ will be written
        Write( 6,'(/4X,"Kw = (0-write CONF.HIJ, 1-no CONF.HIJ) ",I1)') Kw
        Close(99)

        ! If starting new computation with MBPT
        ! Ksig = 0 - no MBPT included (same as Kl = 0)
        ! Ksig = 1 - include 1-electron MBPT corrections 
        ! Ksig = 2 - include 1-electron and 2-electron MBPT corrections
        ! Kdsig = 0 - automatic approximation of the energy dependence of Sigma
        ! Kdsig = 1 - manually include energy dependence of Sigma
        If (Kl == 2) Then
            Write(*,'(1X," Ksig = (0,1,2): ",I1)') Ksig 
            If (Ksig /= 0) Then
                Write(*,'(1X," Energy dependence of Sigma (1-Yes,0-No)? ",I1)') Kdsig
            End If
            Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
            Kecp=0
            Kl=0
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
        Read(16) (In(i),i=1,IPgnt)
        Read(16) (Gnt(i),i=1,IPgnt)
        Close(unit=16)

    End Subroutine Input

    Subroutine Init
        Implicit None
        Integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, If, &
                    ii, i1, n2, n1, l, nmin, jlj, i0, nlmax, err_stat
        Real(dp) :: d, c1, c2, z1
        Real(dp), Dimension(IP6)  :: p, q, p1, q1 
        Real(dp), Dimension(4*IP6):: pq
        Integer, Dimension(33)  ::  nnn ,jjj ,nqq 
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(33):: lll
        logical :: longbasis
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Character(Len=256) :: strfmt, err_msg

        Equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        Data Let/'s','p','d','f','g','h','i','k','l'/

        c1 = 0.01d0
        mj = 2*abs(Jm)+0.01d0

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
        If (abs(Z-z1) > 1.d-6) Then
            strfmt = '("nuc. charge is changed: Z =",F12.6," ><",F12.6)'
            Write( 6,strfmt) Z,z1
            Write(11,strfmt) Z,z1
            Read(*,*)
        End If

        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        Rnuc=pq(13)
        dR_N=pq(16)
        
        strfmt = '(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2,/4X, &
                    "Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3,5X,"Nc =",I6)'
        Write( 6,strfmt) Kl,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,strfmt) Kl,Z,Jm,Nsp,Ns,Nso,Nc

        longbasis=abs(PQ(20)-0.98765d0) < 1.d-6
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
            If=20
            Do ni=1,Ns
                If=If+1
                Nn(ni)=pq(If)+c1
                If=If+1
                Ll(ni)=pq(If)+c1
                If=If+3
                c2=dsign(c1,pq(If))
                Kk(ni)=pq(If)+c2
                If=If+1
                c2=dsign(c1,pq(If))
                Jj(ni)=pq(If)+c2
            End Do
        End If

        Nsu=0
        Do nj=1,Nsp
            i=sign(1.d0,Qnl(nj))
            d=abs(Qnl(nj))+1.d-14
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

        nec=0
        If (Nso /= 0) nec=sum(Nq(1:Nso))
        Do ni=1,Nsu
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                Nst=Nst+1
            End Do
        End Do
        strfmt = '(4X,"Number of actually Used orbitals: Nsu =",I3, &
                    /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)'
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

        Do ic=1,Nc
            n1=Nc0(ic)+1
            n2=Nc0(ic)+Nvc(ic)
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
            strfmt = '(1X,I6,"#",6(I2,A1,"(",I1,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                             /8X,6(I2,A1,"(",I1,"/2)",I2,";"))'
            Write(11,strfmt) ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
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
        Integer :: na, nb, ierr, k, nsh, nx, khot, la, lb, ja, jb, na1, nb1, ind, &
                   nso1, khot1, k1, nsx1, nsx2, nav, kbox, i, ns1, idummy, err_stat, Nmax1, Lmax1
        Real :: x, y
        Character(Len=256) :: strfmt, err_msg
        
        If (Ksig == 0) Return
        nsh=Nso+1

        ! parameter for indexation of integrals:
        nx = IPx

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
                ind=nx*(na-nsh)+(nb-nsh+1)
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
        Allocate(Rint2S(NgintS),Iint2S(NgintS),Iint3S(NgintS),Dint2S(NgintS),Eint2S(NgintS),IntOrdS(nrd))
        If (nrd /= nx*nx) Then
            Write(*,*)' RintS: IPx was changed from nrd=',nrd,'to nx*nx=',nx*nx
            ierr=ierr+1
        End If
        If (ierr > 0) Then
            Write(*,*) ierr,' error(s) detected'
            Stop
        End If
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
                x=x+Eps(ke)
            End Do
            x=E_0-x
            Diag(n)=x
            Emin=dmin1(Emin,x)
            Emax=dmax1(Emax,x)
        End Do
        strfmt = '(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)'
        write( 6,strfmt) Emin,Emax
        write(11,strfmt) Emin,Emax

    End Subroutine FormD
    
    Subroutine AllocateFormHArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr

        vaBinSize = 10000000
        Call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IPlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NhintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NgintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. Allocated(Jz)) Allocate(Jz(Nst))
        If (.not. Allocated(Nh)) Allocate(Nh(Nst))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(Rint1)) Allocate(Rint1(Nhint))
        If (.not. Allocated(Rint2)) Allocate(Rint2(IPbr,Ngint))
        If (.not. Allocated(Iint1)) Allocate(Iint1(Nhint))
        If (.not. Allocated(Iint2)) Allocate(Iint2(Ngint))
        If (.not. Allocated(Iint3)) Allocate(Iint3(Ngint))
        If (.not. Allocated(IntOrd)) Allocate(IntOrd(nrd))
        If (.not. Allocated(Iarr)) Allocate(Iarr(Ne,Nd))

        If (Ksig /= 0) Then
            If (.not. Allocated(Scr)) Allocate(Scr(10))
            If (.not. Allocated(Rsig)) Allocate(Rsig(NhintS))
            If (.not. Allocated(Dsig)) Allocate(Dsig(NhintS))
            If (.not. Allocated(Esig)) Allocate(Esig(NhintS)) 
            If (.not. Allocated(R_is)) Allocate(R_is(num_is))
            If (.not. Allocated(I_is)) Allocate(I_is(num_is))
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
                + sizeof(Rint1)+sizeof(Rint2)+sizeof(Iint1)+sizeof(Iint2)+sizeof(Iint3)+sizeof(Iarr)&
                + sizeof(IntOrd)
            If (Ksig /= 0) memFormH = memFormH+sizeof(Rint2S)+sizeof(Dint2S)+sizeof(Eint2S) &
                + sizeof(Iint1S)+sizeof(Iint2S)+sizeof(Iint3S) &
                + sizeof(Rsig)+sizeof(Dsig)+sizeof(Esig)+sizeof(R_is)+sizeof(I_is)+sizeof(IntOrdS)
            !Call FormattedMemSize(memFormH, memStr)
            !Write(*,'(A,A,A)') 'Allocating arrays for FormH requires ',Trim(memStr),' of memory per core' 
        End If   

        Return
    End Subroutine AllocateFormHArrays

    Subroutine calcMemStaticArrays
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Character(Len=16) :: memStr
        
        memStaticArrays = memStaticArrays! + 209700000_int64 ! static "buffer" of 200 MiB 
        memStaticArrays = memStaticArrays + sizeof(Nn)+sizeof(Kk)+sizeof(Ll)+sizeof(Jj)+sizeof(Nf0) &
                        + sizeof(Jt)+sizeof(Njt)+sizeof(Eps)+sizeof(Diag)+sizeof(Ndc)+sizeof(Jz) &
                        + sizeof(Nh)+sizeof(In)+sizeof(Gnt)+sizeof(C)+16*vaBinSize
        !Call FormattedMemSize(memStaticArrays, memStr)
        !Write(*,'(A,A,A)') 'calcMemReqs: Allocating static arrays requires ',Trim(memStr),' of memory per core' 
    End Subroutine calcMemStaticArrays

    Subroutine calcMemReqs
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer(Kind=int64) :: mem
        Integer :: bytesInteger, bytesDP, bytesReal
        Character(Len=16) :: memStr
    
        bytesInteger = 4
        bytesReal = 4
        bytesDP = 8
    
        Call calcMemStaticArrays
        memEstimate = memFormH + memStaticArrays
        Call FormattedMemSize(memEstimate, memStr)
        Write(*,'(A,A,A)') 'calcMemReqs: Allocating arrays for FormH will require at least ',Trim(memStr),' of memory per core' 
    
        mem = bytesDP * Nd * IPlv     & ! ArrB
            + bytesDP * Nlv * 2_dp    & ! Tk,Tj
            + bytesDP * IPlv ** 2_dp  & ! P
            + bytesDP * IPlv * 2_dp   & ! D,E
            + bytesDP * Nd * 2_dp     & ! B1,B2
            + bytesDP * Nd0 ** 2_dp   & ! Z1
            + bytesDP * Nd0 * 2_dp      ! D1,E1
    
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

    Subroutine InitFormH(npes,mype)
        ! this subroutine initializes variables used for FormH and subsequent subroutines
        ! All necessary variables are broadcasted from root to all cores
        Use mpi
        Implicit None
        Integer :: npes, mype, mpierr, i

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Crt4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
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
        Call MPI_Bcast(Eps(1:IPs), IPs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(In(1:IPgnt), IPgnt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gnt(1:IPgnt), IPgnt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint1(1:Nhint), Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint2(1:IPbr,1:Ngint), IPbr*Ngint, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint1(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint2(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint3(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IntOrd, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Do i=1,Ne
            Call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End do  
        If (Ksig /= 0) Then
            Call MPI_Bcast(Scr, 10, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Esig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint1S(1:NhintS), NhintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint2S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint3S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Eint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(IntOrdS, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(R_is(1:num_is), num_is, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(I_is(1:num_is), num_is, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End subroutine InitFormH

    Subroutine FormH(npes, mype)
        Use mpi
        Use str_fmt, Only : FormattedMemSize, FormattedTime
        Use determinants, Only : calcNd0, Gdet, CompCD, Rspq_phase1, Rspq_phase2
        Use matrix_io
        Use vaccumulator

        Implicit None

        Integer :: npes, mype, mpierr, interval, remainder, maxNumElementsPerCore=0
        Integer :: is, nf, i1, i2, j1, j2, k1, kx, n, ic, ic1, n1, n2, ihmax, &
                   n0, jq, jq0, iq, i, j, icomp, k, ih4, counter1, counter2, counter3, diff, k2, totsize
        Integer :: nn, kk, msg, status(MPI_STATUS_SIZE), sender, num_done, an_id, return_msg, endnd, minme, maxme
        logical :: finished
        Integer, allocatable, dimension(:) :: idet1, idet2, cntarray
        Integer(Kind=int64)     :: start_time, end_time, stot, etot, s1, e1, s2, e2, clock_rate, numzero=0
        real :: ttime, ttot
        real(dp)  :: t, tt
        Integer(Kind=int64) :: statmem, mem, memsum, maxmem, mesplit
        Character(Len=16)     :: memStr, memStr2, memTotStr, memTotStr2, npesStr, counterStr, counterStr2, filename, timeStr
        Integer :: iSign, iIndexes(3), jIndexes(3), nnd
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: growBy, vaGrowBy, ndGrowBy, ndsplit, ndcnt
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002
        Integer(Kind=MPI_OFFSET_KIND) :: disp       

        Call startTimer(stot)
        Call InitFormH(npes,mype) 
        Call calcNd0(ic1, Nd0)

        If (mype == 0) Then
            If (Ksig == 2) Write(*,*) 'Screening is included'
            Call calcMemReqs
        End If

        ! If continuing calculation or Hamiltonian has already been constructed
        If (Kl == 1) Then 
            ! Read the Hamiltonian from file CONFp.HIJ
            Call ReadMatrix(Hamil,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)
            numzero = Count(Hamil%t(1:ih4) == 0)

            ! Add maximum memory per core from storing H to total memory count
            Call MPI_AllReduce(ih4, ihmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            memEstimate = memEstimate + ihmax*16

        ! If starting a new computation and Hamiltonian has not been constructed
        Else 
            Allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne),cntarray(2))
            vaGrowBy = vaBinSize
            ndGrowBy = 1

            If (mype==0) Then
                Open (unit=18,file='CONFp.HIJ',status='UNKNOWN',form='UNFORMATTED')
                Close(unit=18,status='DELETE')
                Write(counterStr,fmt='(I16)') vaGrowBy
                Write(counterStr2,fmt='(I16)') ndGrowBy
                Write(*,'(A)') ' vaGrowBy = '//Trim(AdjustL(counterStr))//', ndGrowBy = '//Trim(AdjustL(counterStr2))
                print*, '========== Starting comparison stage of FormH =========='
            End If

            counter1=1
            counter2=1
            counter3=1
            cntarray=0

            ! Get accumulator vectors setup (or re-setup If this is rank 0):
            Call IVAccumulatorInit(iva1, vaGrowBy)
            Call IVAccumulatorInit(iva2, vaGrowBy)
            If (mype==0) Call RVAccumulatorInit(rva1, vaGrowBy)
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            If (mype == 0) Then        
                ! Distribute a portion of the workload of size ndGrowBy to each worker process
                Do an_id = 1, npes - 1
                   nnd = 1 + ndGrowBy*(an_id-1) + 1
                   Call MPI_SEND( nnd, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                End Do

                n=1
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
                                    tt=Hmltn(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
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

                NumH = cntarray(1)
                num_done = 0
                ndsplit = Nd/10
                ndcnt = ndsplit
                maxme = cntarray(2)
                j=9

                Do 
                    Call MPI_RECV( cntarray, 2, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    sender = status(MPI_SOURCE)
             
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
                    mem = NumH * 16_int64
                    maxmem = maxme * 16_int64
                    statmem = memEstimate + maxmem
                    Call FormattedMemSize(statmem, memTotStr)
                    Call FormattedMemSize(memTotalPerCPU, memTotStr2)

                    If (nnd == ndcnt .and. nnd /= ndsplit*10) Then
                        Call stopTimer(s1, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Write(counterStr,fmt='(I16)') NumH
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH comparison stage:', (10-j)*10, '% done in '// trim(timeStr)// '; '// &
                                                    Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)// &
                                                    ', MaxMemPerCore='//trim(memStr2)//')'
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
                        memEstimate = memEstimate + maxmem
                        Write(counterStr,fmt='(I16)') NumH
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH comparison stage:', (10-j)*10, '% done in '// trim(timeStr)// '; '// &
                                                    Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)// &
                                                    ', MaxMemPerCore='//trim(memStr2)//')'
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
                    
                        Call MPI_SEND(cntarray, 2, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                    End if
                End do
            End If

            Call IVAccumulatorCopy(iva1, Hamil%n, counter1)
            Call IVAccumulatorCopy(iva2, Hamil%k, counter2)
            If (mype == 0) Call RVAccumulatorCopy(rva1, Hamil%t, counter3)

            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
            If (mype == 0) Call RVAccumulatorReset(rva1)
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call startTimer(s1)

            Call MPI_AllReduce(counter1, maxNumElementsPerCore, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            mesplit = maxNumElementsPerCore/10
            numzero=0
            j=1

            If (mype /= 0) Then
                Allocate(Hamil%t(counter1))
                Call startTimer(s2)
                Do n=1,counter1
                    nn=Hamil%n(n)
                    kk=Hamil%k(n)
                    Call Gdet(nn,idet1)
                    Call Gdet(kk,idet2)
                    Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                    Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                    If (Kdsig /= 0 .and. diff <= 2) E_k=Diag(kk)
                    t=Hmltn(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                    Hamil%t(n)=t
                    If (t == 0) numzero=numzero+1

                    If (counter1 == maxNumElementsPerCore .and. mod(n,mesplit)==0) Then
                        Call stopTimer(s1, timeStr)
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH calculation stage:', j*10, '% done in '// trim(timeStr)
                        j=j+1
                    End If
                End Do
                Call stopTimer(s2, timeStr)
                !print*,mype,timeStr
            Else
                print*, '========== Starting calculation stage of FormH =========='
            End If
            Deallocate(idet1, idet2, iconf1, iconf2, cntarray)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            If (mype==0) print*, '========== Formation of Hamiltonian matrix completed =========='
        End If

        ih8=size(Hamil%t)
        ih4=ih8
        
        Call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, numzero, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, iscr, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(MPI_IN_PLACE, xscr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        
        ! Write Hamiltonian to file CONFp.HIJ
        If (Kl /= 1 .and. Kw == 1)  Call WriteMatrix(Hamil,ih4,NumH,'CONFp.HIJ',mype,npes,mpierr)
        
        ! give all cores Hmin, the minimum matrix element value
        Call MPI_AllReduce(Hamil%t(1:ih8), Hmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
        
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
            Write(counterStr,fmt='(F14.8)') Hmin
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

    real(dp) function Hmltn(idet1, idet2, is, nf, i2, i1, j2, j1) 
        ! This function calculates the Hamiltonian matrix element between determinants idet1 and idet2
        Use determinants, Only : Rspq
        Use integrals, Only : Gint, Hint
        Use formj2, Only : F_J0, F_J2
        Implicit None
        Integer, allocatable, dimension(:), intent(inout)   :: idet1, idet2
        Integer, intent(inout)                              :: is, nf, i1, i2, j1, j2
        
        Integer     :: iq, jq, jq0, k
        real(dp)    :: t
        
        t=0.d0
        Select Case(nf)
            Case(2) ! determinants dIffer by two functions
                t=t+Gint(i2,j2,i1,j1)*is 
                t=t+Gj*F_J2(idet1,idet2,is,nf,i2,i1,j2,j1)
            Case(1) ! determinants dIffer by one function
                Do iq=1,Ne
                    i1=idet1(iq)
                    If (i1 /= j1) t=t+Gint(j2,i1,j1,i1)*is
                End Do
                t=t+Hint(j2,j1)*is
            Case(0) ! determinants are equal
                Do iq=1,Ne
                    i1=idet1(iq)
                    jq0=iq+1
                    If (jq0 <= Ne) Then
                        Do jq=jq0,Ne
                            j1=idet1(jq)
                            t=t+Gint(i1,j1,i1,j1)*is
                        End Do
                    End If
                    t=t+Hint(i1,i1)*is
                End Do
                t=t+Gj*F_J2(idet1,idet2,is,nf,i2,i1,j2,j1)
        End Select
        Hmltn=t
        Return
    End function Hmltn

    Subroutine DeAllocateFormHArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    
        If (mype==0) Then
            Call FormattedMemSize(memFormH, memStr)
            Write(*,'(A,A,A)') 'De-allocating ',Trim(memStr),' of memory per core from arrays for FormH'  
        End If
    
        If (Allocated(Nvc)) Deallocate(Nvc)
        If (Allocated(Nc0)) Deallocate(Nc0)
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

    Subroutine AllocateDvdsnArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(ArrB)) Allocate(ArrB(Nd,IPlv))
        If (.not. Allocated(Tk)) Allocate(Tk(Nlv))
        If (.not. Allocated(Tj)) Allocate(Tj(Nlv))
        If (.not. Allocated(P)) Allocate(P(IPlv,IPlv))
        If (.not. Allocated(D)) Allocate(D(IPlv))
        If (.not. Allocated(E)) Allocate(E(IPlv))
        If (.not. Allocated(Iconverge)) Allocate(Iconverge(IPlv))
        If (.not. Allocated(B1)) Allocate(B1(Nd))
        If (.not. Allocated(B2)) Allocate(B2(Nd))
        If (.not. Allocated(Z1)) Allocate(Z1(Nd0,Nd0))
        If (.not. Allocated(D1)) Allocate(D1(3*Nd0-1))
        If (.not. Allocated(E1)) Allocate(E1(Nd0))
        If (mype==0) Then
            memDvdsn = 0_int64
            memDvdsn = sizeof(ArrB)+sizeof(Tk)+sizeof(Tj)+sizeof(P)+sizeof(D)+sizeof(E) &
                + sizeof(Iconverge)+sizeof(B1)+sizeof(B2)+sizeof(Z1)+sizeof(D1)+sizeof(E1)
            memDvdsn = memDvdsn
            Call FormattedMemSize(memDvdsn, memStr)
            Write(*,'(A,A,A)') 'Allocating arrays for Davidson procedure requires ',Trim(memStr),' of memory per core'  
            memEstimate = memEstimate - memFormH + memDvdsn
            Call FormattedMemSize(memEstimate, memStr)
            Write(*,'(A,A,A)') 'Total memory estimate for Davidson procedure is ',Trim(memStr),' of memory per core' 
        End If   

    End Subroutine AllocateDvdsnArrays

    Subroutine Diag4(mype, npes)
        ! this Subroutine executes the Davidson procedure for diagonalization
        ! the Subroutine Mxmpy is a computational bottleneck and was the only Subroutine to be parallelized
        ! all other Subroutines are performed by the master core
        Use mpi
        Use str_fmt, Only : startTimer, stopTimer, FormattedTime
        Use davidson
        Implicit None
        External DSYEV

        Integer  :: k1, k, i, j, n1, iyes, id, id2, id1, ic, id0, kskp, iter, &
                    kx, i1, i2, it, mype, npes, mpierr, ifail=0
        Integer(Kind=int64) :: start_time, s1
        Real :: ttot
        Real(dp)  :: crit=1.d-6, ax, x, xx, ss, cnx, vmax
        logical   :: lsym
        Character(Len=16) :: timeStr

        Call startTimer(start_time)

        If (mype == 0) Then
            If (Kl4 == 0) Return ! Read CONF.XIJ and make 1 iteration
            If (Nc4 > Nc) Nc4=Nc
            Write (*,*) 'kl4=',Kl4,'  Nc4=',Nc4,'  Crt4=',Crt4
        End If

        Call Init4(mype,npes)
        If (mype == 0) Then
            Call startTimer(s1)
            Call DSYEV('V','U',Nd0,Z1,Nd0,E1,D1,3*Nd0-1,ifail)
            !Call Hould(Nd0,D1,E1,Z1,ifail)
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A)'), 'TIMING >>> Initial diagonalization took '// trim(timeStr) // ' to complete'
        End If

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Do While (ifail /= 0)
           If (mype == 0) Write(6,'(4X,"Starting approximation of dim ",I4," failed")') Nd0
           Call Init4(mype,npes)
           !If (mype == 0) Call Hould(Nd0,D1,E1,Z1,ifail)
           If (mype == 0) Call DSYEV('V','U',Nd0,Z1,Nd0,E1,D1,3*Nd0-1,ifail)
        End Do

        Call FormB0(mype,npes)

        If (Nd0 == Nd) Return
        ! Davidson loop:
        iter = 0
        kdavidson = 0
        If (mype==0) Write(*,*) 'Start with kdavidson =', kdavidson
        ax = 1
        cnx = 1
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Do it=1,N_it
            iter=iter+1
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(ax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            If (ax > crit .and. iter <= n_it) Then
                If (mype == 0) Then
                    strfmt = '(1X,"** Davidson iteration ",I3," for ",I2," levels **")'
                    Write( 6,strfmt) iter,Nlv
                    Write(11,strfmt) iter,Nlv
                End if
                ! Orthonormalization of Nlv probe vectors:
                If (mype==0) Then
                    Do i=1,Nlv
                        Call Ortn(i,ifail)
                        If (ifail /= 0) Then
                            Write(*,*)' Fail of orthogonalization for ',i
                            Stop
                        End If
                    End Do
                End If
                ! Formation of the left-upper block of the energy matrix P:
                Call Mxmpy(1, mype, npes)
                Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
                If (mype==0) Then ! start master core only block
                  Call FormP(1, vmax)
                  lsym=K_prj == 1.OR.kskp == 1.OR.kskp == 3
                  If (iter == 1 .and. lsym) Then    !# averaging diagonal over confs
                    id0=1
                    Do ic=1,Nc
                      id1=Ndc(ic)
                      id2=id0+id1-1
                      If (id1 > 0) Then
                        ss=0.d0
                        Do id=id0,id2
                          ss=ss+Diag(id)
                        End Do
                        ss=ss/id1
                        Diag(id0:id2)=ss
                        id0=id0+id1
                      End If
                    End Do
                    If (id2 == Nd) Then
                      Write(*,*) ' Diagonal averaged over rel. configurations'
                    Else
                      Write(*,*) ' Error: id2=',id2,' Nc=',Nc
                      Stop
                    End If
                  End If 
                  ! Formation of Nlv additional probe vectors:
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
                End If

                If (K_prj == 1) Then
                    Call Prj_J(Nlv+1,Nlv,2*Nlv+1,1.d-5,mype,npes)
                    If (mype == 0) Then
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
                                exit
                              End If
                            End If
                          End If
                        End Do
                    End If
                End If
                Call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
                If (cnx > Crt4) Then
                    ! Formation of other three blocks of the matrix P:
                    Call Mxmpy(2, mype, npes)
                    If (mype==0) Then
                        Call FormP(2,vmax)
                        n1=2*Nlv
                        ! Evaluation of Nlv eigenvectors:
                        Call Hould(n1,D,E,P,ifail)
                        !Call DSYEV('V','U',n1,P,n1,E,D,3*n1-1,ifail)
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
                            If (mype == 0) Then
                                strfmt = '(1X,"E(",I2,") =",F14.8,"; admixture of vector ",I2,": ",F10.7)'
                                Write( 6,strfmt) i,-(E(i)+Hmin),kx,xx
                                Write(11,strfmt) i,-(E(i)+Hmin),kx,xx
                            End If
                        End Do
                        If (kXIJ > 0) Then ! Write intermediate CONF.XIJ
                            If (mod(iter,kXIJ) == 0) Then
                                Call FormB
                            Else
                                Call FormBskip
                            End If
                        Else
                            Call FormBskip
                        End If
                    End If
                Else
                    If (mype == 0) Then
                        strfmt = '(" Davidson procedure converged")'
                        Write( 6,strfmt)
                        Write(11,strfmt)
                        ! Assign values of Tk for case when Davidson procedure is already converged 
                        If (Kl4 == 2 .and. iter == 1) Tk=E(1:Nlv) 
                    End If
                    Exit
                End If
            End If
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End Do

        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            Write(*,'(2X,A)'), 'TIMING >>> Davidson procedure took '// trim(timeStr) // ' to complete'
        End If

        Deallocate(Hamil%n, Hamil%k, Hamil%t, Diag, P, D, E, B1, B2, Z1, E1)

        Return
    End Subroutine Diag4

    Subroutine WriteFinalXIJ(mype,npes)
        Use mpi
        Use formj2, Only : J_av
        Use davidson, Only : Prj_J, FormBskip
        Implicit None
        Integer :: i, n, ierr, mype, npes, mpierr

        If (K_prj == 1) Then
            Call Prj_J(1,Nlv,Nlv+1,1.d-8,mype,npes)
        End If
        If (mype == 0) Then
            open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        End If

        Do n=1,Nlv
            Call MPI_Bcast(ArrB(1:Nd,n), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call J_av(ArrB(1,n),Nd,Tj(n),ierr,mype,npes)  ! calculates expectation values for J^2
            If (mype==0) Then
                write(17) Tk(n),Tj(n),Nd,(ArrB(i,n),i=1,Nd)
            End If
        End Do

        Deallocate(ArrB, Tk, Tj, D1, Iconverge)

        If (mype==0) close(unit=17)

        Return
    End Subroutine WriteFinalXIJ

    Subroutine PrintEnergies
        Implicit None

        Integer :: j1, nk, k, j2, n, ndk, ic, i, j, idum, ist, jmax, imax, &
                   j3
        real(dp) :: cutoff, xj, dt, del, dummy, wmx, E, D
        real(dp), allocatable, dimension(:)  :: Cc, Dd, Er
        Character(Len=1), dimension(11) :: st1, st2 
        Character(Len=1), dimension(10)  :: stecp*7
        Character(Len=1), dimension(2)  :: st3*3
        Character(Len=1), dimension(4)  :: strsms*6
        Character(Len=1), dimension(3)  :: strms*3
        data st1/11*'='/,st2/11*'-'/,st3/' NO','YES'/
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/

        Allocate(Cc(Nd), Dd(Nd), W(Nc,IPlv), Er(Nlv))
        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is Used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7
        Open(unit=16,file='CONF.XIJ',status='UNKNOWN',form='unformatted')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! printing eigenvalues in increasing order
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        WRITE( 6,5)
        WRITE(11,5)
  5     FORMAT(4X,63('='))
        If (Ksig*Kdsig == 0) Then
           WRITE( 6,15) stecp(ist),Nc,Nd,Gj
           WRITE(11,15) stecp(ist),Nc,Nd,Gj
 15        FORMAT(4X,'Energy levels (',A7,' Nc=',I6,' Nd=',I9, &
                 ');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        Else
           WRITE( 6,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
           WRITE(11,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
 25        FORMAT(4X,'Energy levels ',A7,', Sigma(E =', &
                 F10.4,') extrapolation var.',I2,/4X,'(Nc=',I6, &
                 ' Nd=',I9,');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        End If
        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                Write( *,251) C_is,Rnuc
                Write(11,251) C_is,Rnuc
 251            format(4X,'Volume shIft: dR_N/R_N=',F9.5,' Rnuc=',F10.7)
            Else
                Write( *,253) strms(K_is-1),C_is,strsms(K_sms),Klow
                Write(11,253) strms(K_is-1),C_is,strsms(K_sms),Klow
 253            format(4X,A3,':',E9.2,'*(P_i Dot P_k) ',A6, &
                     ' Lower component key =',I2)
            End If
        End If
        Write( 6,255)
        Write(11,255)
 255    format(4X,63('-'))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        JMAX=min(Nlv,Nd)
        Do J=1,JMAX
            REWIND (16)
            IMAX=J-1
            IF (IMAX < 1) GOTO 250
            DO I=1,IMAX
                READ (16)
            END DO
 250        READ (16) ER(j),xj,idum,(CC(I),I=1,ND)
            Er(j)=Er(j)+4.d0*Gj*xj*(xj+1.d0)
            E=ER(j)
            DT=E-ECORE
            ! Rydberg constant is taken from "phys.par"
            DEL=(ER(1)-ER(J))*2*DPRy
            WRITE( 6,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
            WRITE(11,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
        End Do
        WRITE( 6,45)
        WRITE(11,45)
 45     FORMAT(4X,63('='))
        ! weights of configurations
        DO J=1,JMAX
            REWIND (16)
            IMAX=J-1
            IF (IMAX < 1) GOTO 270
            DO I=1,IMAX
                READ (16)
            END DO
 270        READ (16) D,DUMMY,idum,(CC(I),I=1,ND)
            I=0
            DO IC=1,NC
                D=0.d0
                NDK=NDC(IC)
                DO K=1,NDK
                    I=I+1
                    D=D+CC(I)**2
                END DO
                W(IC,J)=D
            END DO
        END DO
        N=(JMAX-1)/5+1
        J2=0
        DO 130 K=1,N
            NK=5
            IF (K == N) NK=JMAX-(N-1)*5
            J1=J2+1
            J2=J2+NK
            J3=J2+1
            WRITE(11,65) (ST1,I=J1,J3)
 65         FORMAT(3X,66A1)
            WRITE(11,75) (I,I=J1,J2)
 75         FORMAT(15X,5(3X,I2,6X))
            WRITE(11,65) (ST2,I=J1,J3)
            WRITE(11,85) (ER(I),I=J1,J2)
 85         FORMAT(4X,'ICONF',4X,5F11.5)
            WRITE(11,65) (ST2,I=J1,J3)
            DO 140 IC=1,NC
                I=IC
                WRITE(11,95) I,(W(I,J),J=J1,J2)
 95             FORMAT(2X,I6,'      ',5F11.6)
 140        CONTINUE
            WRITE(11,65) (ST1,I=J1,J3)
 130    CONTINUE
        CLOSE(unit=16)
        close(unit=6)
        close(unit=11)
        Deallocate(Cc, Dd, W, Ndc, Er)
        RETURN
    End Subroutine PrintEnergies
End Program conf
