Program conf_lsj
    use mpi
    use davidson, Only : Prj_J
    use determinants, Only : Dinit, Jterm
    use integrals, Only : Rint
    use conf_init, Only : InitFormH
    use formj2, Only : FormJ, J_av
    Use str_fmt, Only : FormattedTime
    Use lsj_aux

    Implicit None

    Integer   :: n, k, i, j, ierr, mype, npes, mpierr, nnd
    Integer(kind=int64) :: clock_rate
    Integer(kind=int64) :: start_time, end_time
    Real :: total_time
    Real(dp)  :: t, xtj, xtl, xts, xj
    Character(Len=1024) :: strFromEnv
    Character(Len=255)  :: eValue, strfmt
    Character(Len=16)   :: memStr, timeStr

    ! Initialize MPI
    Call MPI_Init(mpierr)
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    Call system_clock(count_rate=clock_rate)
    If (mype==0) Call system_clock(start_time)

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    Call Get_Environment_Variable("CONF_MAX_BYTES_PER_CPU",eValue)
    read(eValue,'(I12)') memTotalPerCPU

    ! Only the master core needs to initialize the conf program
    If (mype == 0) Then
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

Contains

    Subroutine Input
        ! This subroutine reads in parameters and configurations from CONF.INP
        Use conf_init, only : inpstr, ReadConfInp, ReadConfigurations
        Implicit None
        Integer  :: i, i1, i2, ic, nx, ny, nz, ne0, n, k
        Character(Len=1) :: name(16)
        Character(Len=32) :: strfmt

        ! Write name of program
        open(unit=11,status='UNKNOWN',file='CONF_LSJ.RES')
        strfmt = '(4X,"Program conf_lsj v0.1.0")'
        Write( 6,strfmt)
        Write(11,strfmt)
        ! - input from the file 'CONF.INP' - - - - - - - - - - - - - - - -
        Call ReadConfInp

        If (dabs(C_is) < 1.d-6) K_is=0
        If (K_is == 0) C_is=0.d0

        Open(unit=99,file='c.in',status='OLD')
        Read (99,*) Kl, Ksig, Kdsig
        Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
        If (K_is == 2.OR.K_is == 4) Then
        Read(99,*) K_sms
        Write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): ', K_sms
        If ((K_sms-1)*(K_sms-2)*(K_sms-3) /= 0) Stop
        End If
        Close(99)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call ReadConfigurations
        ! - - - - - - -  Case kl = 2  - - - - - - - - - - -
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
        ! - - - - - - -  Case kv = 1,3  - - - - - - - - - -
        K_prj=0                    !# this key is fixed for kv=2,4
        If (Kv == 1 .or. Kv == 3) Then
        K_prj=1
        Write( *,'(4X,"Selection of states with J =",F5.1)') XJ_av
        Write(11,'(4X,"Selection of states with J =",F5.1)') XJ_av
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        !If (Kl /= 1) Then ! If Kl=1, continue from a previous calculation
        !      Open(unit=16,status='UNKNOWN',file='CONF.JJJ')
        !      Close(unit=16,status='DELETE')
        !End If
        Open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        Read(16) (In(i),i=1,IPgnt)
        Read(16) (Gnt(i),i=1,IPgnt)
        Close(unit=16)
        Return
    End Subroutine Input

    Subroutine Init
        Implicit None
        Integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, If, &
                    ii, i1, n2, n1, l, nmin, jlj, i0, nlmax
        Real(dp) :: d, c1, c2, z1
        Real(dp), Dimension(IP6)  :: p, q, p1, q1 
        Real(dp), Dimension(4*IP6):: pq
        Integer, Dimension(33)  ::  nnn ,jjj ,nqq 
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(33):: lll
        logical :: longbasis
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        Data Let/'s','p','d','f','g','h','i','k','l'/
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        c1 = 0.01d0
        mj = 2*abs(Jm)+0.01d0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Open(12,file='CONF.DAT',status='OLD',access='DIRECT',recl=2*IP6*IPmr,err=700)
        Read(12,rec=1) p
        Read(12,rec=2) q
        Read(12,rec=5) p1
        Read(12,rec=6) q1
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        z1 = pq(1)
        If (abs(Z-z1) > 1.d-6) Then
            Write( 6,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Write(11,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Read(*,*)
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        Rnuc=pq(13)
        dR_N=pq(16)
        longbasis=abs(PQ(20)-0.98765d0) < 1.d-6
        Write( 6,'(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kl,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,'(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kl,Z,Jm,Nsp,Ns,Nso,Nc
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
                    Write( 6,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                    Write(11,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
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
        Write( 6,'(4X,"Number of actually Used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
        Write(11,'(4X,"Number of actually Used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
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
                    Write( 6,'(/2X,"wrong number of electrons", &
                         /2X,"for configuration ICONF =",I6)') ic
                    Write(11,'(/2X,"wrong number of electrons", &
                         /2X,"for configuration ICONF =",I6)') ic
                  Stop
                End If
                Nvc(ic)=i
                Nc0(ic)=Nso+i0
                i0=i0+i
                n=0
                i=0
            End If
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Write( 6,'(1X,71("="))')
        Write(11,'(1X,71("="))')
        Do ni=1,Nso
            l =Ll(ni)+1
            lll(ni)=let(l)
        End Do
        Write(11,'(1X,"Core:", 6(I2,A1,"(",I1,"/2)",I2,";"),  &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                        (Nn(i),lll(i),Jj(i),Nq(i),i=1,Nso)
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
                   Write(11,'(/2X,"wrong number of electrons"/ &
                        2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
                        " (",I2,")")') ni,nnn(i1),lll(i1),jjj(i1),nqq(i1)
                 Stop
                End If
            End Do
            n=n2-n1+1
            Write(11,'(1X,I6,"#",6(I2,A1,"(",I1,"/2)",I2,";"), &
                 /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                 /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                 /8X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                 ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
        End Do
        Write(11,'(1X,71("="))')
        If (Ksig > 0) Then
            Do ni=Nso+1,Nsu
               Read(12,rec=2*ni+7) p
               Eps(ni)=-p(ii+1)
            End Do
            Write(11,'(" HF energies are Read from DAT file", &
                   /5(I5,F10.6))') (i,Eps(i),i=Nso+1,Nsu)
        End If
        Close(unit=12)
        ! Maximal number of eigenvectors Nlv:
        nlmax=IPlv
        If (Kv >= 3) Then
            nlmax=IPlv/3
        End If
        If (Nlv > nlmax) Nlv=nlmax
        Return
        !     - - - - - - - - - - - - - - - - - - - - - - - - -
700     Write( 6,'(/2X,"file CONF.DAT is absent"/)')
        Write(11,'(/2X,"file CONF.DAT is absent"/)')
        Stop
        !     - - - - - - - - - - - - - - - - - - - - - - - - -
    End subroutine Init

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
End Program conf_lsj
