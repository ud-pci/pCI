Program conf_lsj
    use mpi_f08
    Use conf_variables
    use davidson, Only : Prj_J
    use determinants, Only : Dinit, Jterm
    use integrals, Only : Rint
    use formj2, Only : FormJ, J_av
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime

    Implicit None

    Integer   :: i, mype, npes, mpierr, nnd, rec1, rec2, nlvs
    Integer(kind=int64) :: s1, start_time
    Real(dp), Allocatable, Dimension(:) :: xj, xl, xs
    Real(dp), Allocatable, Dimension(:,:) :: B1h
    Character(Len=255)  :: eValue
    Character(Len=16)   :: timeStr

    Type(MPI_Datatype) :: mpi_type2_real

    ! Initialize MPI
    Call MPI_Init(mpierr)
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    Call startTimer(start_time)

    Select Case(type2_real)
    Case(sp)
        mpi_type2_real = MPI_REAL
    Case(dp)
        mpi_type2_real = MPI_DOUBLE_PRECISION
    End Select

    ! Read total memory per core from environment 
    ! Have to export CONF_MAX_BYTES_PER_CPU before job runs
    Call Get_Environment_Variable("CONF_MAX_BYTES_PER_CPU",eValue)
    read(eValue,'(I12)') memTotalPerCPU

    ! Only the master core needs to initialize the conf program
    If (mype == 0) Then
        Call Input          ! reads list of configurations from CONF.INP
        Call Init           ! reads basis set information from CONF.DAT
        Call Rint           ! reads radial integrals from CONF.INT
        Call Dinit          ! forms list of determinants
        Call Jterm          ! prints table with numbers of levels with given J
        Call ReadXIJ(B1h)   ! reads relevant wavefunctions from CONF.XIJ
    End If

    Call AllocateLSJArrays(mype)
    Call InitLSJ
    Call startTimer(s1)
    Call lsj(B1h,xj,xl,xs,mype,npes)

    ! Print table of final results and total computation time
    If (mype==0) Then
        Call stopTimer(s1, timeStr)
        Write(*,'(2X,A)') 'LSJ done in '// trim(timeStr)
        Call PrintLSJ   !#   Output of the results

        Call stopTimer(start_time, timeStr)
        write(*,'(2X,A)') 'TIMING >>> Total computation time of conf_lsj was '// trim(timeStr)
    End If

    Call MPI_Finalize(mpierr)

Contains

    Subroutine Input
        ! This subroutine reads in parameters and configurations from CONF.INP
        Use conf_init, only : ReadConfInp, ReadConfigurations
        Implicit None
        Character(Len=32) :: strfmt

        ! Read range of desired energy levels
        Write(*,'(/A)')'Give initial and final number of the level in CONF.RES'
        Read (*,*) rec1, rec2
        Write( 6,'(A,I2,A,I2)')' No.1=',rec1,'  No.2=',rec2

        ! Write name of program
        open(unit=11,status='UNKNOWN',file='CONF_LSJ.RES')
        Select Case(type2_real)
        Case(sp)
            strfmt = '(4X,"Program conf_lsj v2.4")'
        Case(dp)            
            strfmt = '(4X,"Program conf_lsj v2.4 with double precision for 2e integrals")'
        End Select
        Write( 6,strfmt)
        Write(11,strfmt)

        strfmt = '(" N1=",I2,"  N2=",I2)'
        Write(11,strfmt) rec1,rec2

        ! Read input parameters from file CONF.INP
        Call ReadConfInp

        ! Read configurations from file CONF.INP
        Call ReadConfigurations

        Open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        Read(16) Ngaunt
        Allocate(In(Ngaunt))
        Allocate(Gnt(Ngaunt))
        Read(16) (In(i),i=1,Ngaunt)
        Read(16) (Gnt(i),i=1,Ngaunt)
        Close(unit=16)
        Return
    End Subroutine Input

    Subroutine Init
        Use utils, Only : DetermineRecordLength
        Implicit None
        Integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, If, &
                    ii, i1, n2, n1, l, nmin, jlj, i0, nlmax, err_stat
        Real(dp) :: d, c1, c2, z1
        Real(dp), Dimension(IP6)  :: p, q, p1, q1 
        Real(dp), Dimension(4*IP6):: pq
        Integer, Dimension(33)  ::  nnn ,jjj ,nqq 
        Character(Len=1), Dimension(9) :: Let 
        Character(Len=1), Dimension(33):: lll
        logical :: longbasis, success
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Character(Len=256) :: strfmt, err_msg

        Equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        Data Let/'s','p','d','f','g','h','i','k','l'/

        c1 = 0.01d0
        mj = 2*abs(Jm)+0.01d0
        
        Call DetermineRecordLength(Mrec, success)
        If (.not. success) Then
            Write(*,*) 'ERROR: record length could not be determined'
            Stop
        End If

        Open(12,file='CONF.DAT',status='OLD',access='DIRECT',recl=2*IP6*Mrec,iostat=err_stat,iomsg=err_msg)
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

        Close(unit=12)

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
                    strfmt = '(/2X,"wrong number of electrons", &
                                /2X,"for configuration ICONF =",I6)'
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
                    strfmt = '(/2X,"wrong number of electrons",/2X,"for the shell:", &
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

        ! Maximal number of eigenvectors Nlv:
        nlmax=IPlv
        If (Kv >= 3) Then
            nlmax=IPlv/3
        End If
        If (Nlv > nlmax) Nlv=nlmax
        Return

    End subroutine Init

    Subroutine ReadXIJ(B1h)
        Implicit None
        Real(dp), Allocatable, Dimension(:,:), Intent(InOut) :: B1h

        Integer :: n, i, j, nnd
        Character(Len=16) :: strfmt, recStr1, recStr2

        nlvs = rec2-rec1+1
        If (.not. Allocated(B1h)) Allocate(B1h(Nd,nlvs))
        If (.not. Allocated(D1)) Allocate(D1(nlvs))
        If (.not. Allocated(Tj)) Allocate(Tj(nlvs))

        Open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        Write(recStr1,'(I4)') rec1
        Write(recStr2,'(I4)') rec2
        strfmt = '(A)'
        Write(*,strfmt) 'Starting lsj procedure for levels ' // Trim(AdjustL(recStr1)) // ' to ' // Trim(AdjustL(recStr2))
        n=1
        Do j=1,Nlv
            Read(16) D1(n),Tj(n),nnd,(B1h(i,n),i=1,Nd)
            If (j < rec1) Then
                Cycle
            Else If (j >= rec2) Then
                Exit
            Else
                n=n+1
            End If
        End Do
        Close(16)
    End Subroutine ReadXIJ

    Subroutine AllocateLSJArrays(mype)
        use mpi_f08
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype

        Call MPI_Bcast(nlvs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
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
        Call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. Allocated(Jz)) Allocate(Jz(Nst))
        If (.not. Allocated(Nh)) Allocate(Nh(Nst))
        If (.not. Allocated(Nh0)) Allocate(Nh0(Nst))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(In)) Allocate(In(Ngaunt))
        If (.not. Allocated(Gnt)) Allocate(Gnt(Ngaunt))
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
        If (.not. Allocated(D1)) Allocate(D1(nlvs))
        If (.not. Allocated(Tj)) Allocate(Tj(nlvs))
        If (.not. Allocated(Xj)) Allocate(Xj(nlvs))
        If (.not. Allocated(Xl)) Allocate(Xl(nlvs))
        If (.not. Allocated(Xs)) Allocate(Xs(nlvs))
        If (.not. Allocated(B1h)) Allocate(B1h(Nd,nlvs))

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
    End Subroutine AllocateLSJArrays

    Subroutine InitLSJ
        ! this subroutine initializes variables used for LSJ and subsequent subroutines
        ! All necessary variables are broadcasted from root to all cores
        Use mpi_f08
        Use mpi_utils
        Implicit None
        Integer :: mpierr, i
        Integer(Kind=int64) :: count

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(rec1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(rec2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Crt4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
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
        Call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_sms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kdsig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kexn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Eps(1:IPs), IPs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(In(1:Ngaunt), Ngaunt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gnt(1:Ngaunt), Ngaunt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh0(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint1(1:Nhint), Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        If (Kbrt == 0) Then
            count = Ngint
        Else
            count = Ngint*2_int64
        End If
        Call BroadcastD(Rint2, count, 0, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint1(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastI(Iint2, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastI(Iint3, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IntOrd, IPx*IPx, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd,kind=int64)
        Call BroadcastI(Iarr, count, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nnd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tj(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(D1(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Do i=1,nlvs
            Call MPI_Bcast(B1h(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        End Do

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End subroutine InitLSJ

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
                        do m=1,nlvs
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

    Subroutine lsj(cc,xj,xl,xs,mype,npes)
        use mpi_f08
        Use str_fmt, Only : startTimer, stopTimer
        Implicit None
        Real(dp), Allocatable, Dimension(:), Intent(InOut) :: xj, xl, xs
        Real(dp), Allocatable, Dimension(:,:), Intent(In) :: cc
        Integer(Kind=int64) :: s1, s2
        Integer :: mype, npes, mpierr, msg, ncsplit, nnc, nccnt, an_id, ncGrowBy, endnc, num_done
        Integer :: n, j, sender
        Type(MPI_STATUS) :: status
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
                   !print*,nnc,'to',an_id
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
                        !print*,nnc,'to',sender
                    Else
                        msg = -1
                        Call MPI_SEND( msg, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If

                    If (nnc == nccnt .and. nnc /= ncsplit*10) Then
                        Call stopTimer(s1, timeStr)
                        If (moreTimers) Write(*,'(2X,A,1X,I3,A)') 'lsj:', (10-j)*10, '% done in '// trim(timeStr)
                        j=j-1
                        nccnt = nccnt + ncsplit
                    End If

                    If (num_done == npes-1) Then
                        Call stopTimer(s1, timeStr)
                        If (moreTimers) Write(*,'(2X,A,1X,I3,A)') 'lsj:', (10-j)*10, '% done in '// trim(timeStr)
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
                        !Print*,mype,nnc,timeStr
                    End If
                End Do
            End If

            Call MPI_AllReduce(MPI_IN_PLACE, xj, nlvs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(MPI_IN_PLACE, xl, nlvs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(MPI_IN_PLACE, xs, nlvs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End If

        Deallocate(plj, pls, p0s, pll, p0l)        

        Xj=0.5d0*(dsqrt(Xj+1)-1)
        Xl=0.5d0*(dsqrt(Xl+1)-1)
        Xs=0.5d0*(dsqrt(Xs+1)-1)
            
        Return
    End Subroutine lsj

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

    Subroutine PrintLSJ
        Implicit None

        Integer :: j1, nk, k, j2, n, ndk, ic, i, j, idum, ist, jmax, imax, &
                   j3
        real(dp) :: xj, dt, del, dummy, E, D, gfactor
        real(dp), allocatable, dimension(:)  :: Cc, Dd, Er
        Character(Len=1), dimension(11) :: st1, st2 
        Character(Len=1), dimension(10)  :: stecp*7
        Character(Len=1), dimension(4)  :: strsms*6
        Character(Len=1), dimension(3)  :: strms*3
        data st1/11*'='/,st2/11*'-'/
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/
        Character(Len=256) :: strfmt, strfmt2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Cc(Nd), Dd(Nd), W(Nc,IPlv), Er(Nlv))
        ist=1+3*Kbrt          !### stecp(ist) is Used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7
        Open(unit=16,file='CONF.XIJ',status='UNKNOWN',form='unformatted')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! printing eigenvalues in increasing order
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        strfmt = '(1X,80("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        strfmt = '(" Energy levels (",A7," Nc=",I7," Nd=",I9,"); &
                      Gj =",F7.4, & 
                    /"  N",5X,"JTOT",5X,"L",7X,"S",4X,"G-factor", &
                    5X,"EV",15X,"ET",8X,"DEL(CM**-1)")'
        Write( 6,strfmt) stecp(ist),Nc,Nd,Gj
        Write(11,strfmt) stecp(ist),Nc,Nd,Gj

        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                strfmt = '(4X,"Volume shift: dR_N/R_N=",F9.5," Rnuc=",F10.7)'
                Write( *,strfmt) C_is,Rnuc
                Write(11,strfmt) C_is,Rnuc
            Else
                strfmt = '(4X,A3,":",E9.2,"*(P_i Dot P_k) ",A6, &
                     " Lower component key =",I2)'
                Write( *,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
                Write(11,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
            End If
        End If

        strfmt = '(1X,80("-"))'
        Write( 6,strfmt)
        Write(11,strfmt)

        jmax=min(Nlv,Nd)
        n=1
        Do j=1,jmax
            Rewind(16)
            imax=j-1
            If (imax >= 1) Then
                Do i=1,imax
                    READ(16)
                End Do
            End If
            Read(16) Er(J),xj,idum,(CC(I),i=1,Nd)
            Er(j)=Er(j)+4.d0*Gj*xj*(xj+1.d0)
            E=Er(J)
            DT=E-Ecore
            ! Rydberg constant is taken from "phys.par"
            DEL=(ER(1)-ER(J))*2*DPRy
            gfactor=(3*xj*(xj+1)-Xl(j)*(Xl(j)+1)+Xs(j)*(Xs(j)+1))/(2*xj*(xj+1))
            strfmt = '(I3,F10.6,3F8.4,F14.8,F15.6,F15.2)'
            If (j >= rec1 .and. j <= rec2) Then
                Write( 6,strfmt) j,xj,Xl(n),Xs(n),gfactor,E,DT,DEL
                Write(11,strfmt) j,xj,Xl(n),Xs(n),gfactor,E,DT,DEL
                n=n+1
            End If
        End Do

        strfmt = '(1X,80("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! weights of configurations
        Do j=1,jmax
            Rewind(16)
            imax=j-1
            If (IMAX >= 1) Then
                Do i=1,imax
                    Read(16)
                End Do
            End If
            Read(16) D,DUMMY,idum,(CC(I),i=1,Nd)
            i=0
            Do ic=1,Nc
                D=0.d0
                ndk=Ndc(ic)
                Do k=1,ndk
                    i=i+1
                    D=D+CC(i)**2
                End Do
                W(ic,j)=D
            End Do
        End Do

        n=(jmax-1)/5+1
        j2=0

        Do k=1,n
            nk=5
            If (k == n) nk=jmax-(n-1)*5
            j1=j2+1
            j2=j2+nk
            j3=j2+1
            strfmt2 = '(3X,66A1)'
            Write(11,strfmt2) (st1,i=j1,j3)
            strfmt = '(15X,5(3X,I2,6X))'
            Write(11,strfmt) (i,i=j1,j2)
            Write(11,strfmt2) (st2,i=j1,j3)
            strfmt = '(4X,"ICONF",4X,5F11.5)'
            Write(11,strfmt) (Er(i),i=j1,j2)
            Write(11,strfmt2) (st2,i=j1,j3)
            strfmt = '(2X,I6,"      ",5F11.6)'
            Do ic=1,Nc
                i=ic
                Write(11,strfmt) i,(W(i,j),j=j1,j2)
            End Do
            Write(11,strfmt2) (st1,i=j1,j3)
        End do
        Close(unit=16)
        Close(unit=6)
        Close(unit=11)
        Deallocate(Cc, Dd, W, Ndc, Xl, Xs, Er)
        Return
    End Subroutine PrintLSJ
End Program conf_lsj
