Program conf_pt   
    ! #############################################################################
    ! # This version of conf_pt forms corrections to eigenvectors in CONF_PT.XIJ           
    ! # and forms ordered CONF_new.INP
    ! #############################################################################
    ! This code forms the second order corrections to eigenvalues from CONF.XIJ. 
    ! It is supposed that vectors in CONF.XIJ correspond to the same Hamiltonian 
    ! in a smaller space.
    use mpi_f08
    Use conf_variables
    Use integrals, Only : Rint
    Use determinants, Only : Dinit, Jterm, Wdet
    Implicit None
    
    ! global variables for conf_pt
    integer, parameter :: IPPT = 5000
    integer     :: Nd1, Ncci, Nmax, Ncp0, Ncnr, Ncnrci, Ncnr0, Ncnew, &
                   NcOld, n_is, KmaxScr, KsymScr, NsumScr, MaxScr, ktf, kvar
    real(dp)    :: dvnrx, dvnrn
    logical     :: average_diag

    Integer  :: Nc_0, mpierr, mype, npes
    Real :: start_time, stop_time

    integer, allocatable, dimension(:)    :: Ndcnr, Nvcnr, NRR, NRN, Ndirc
    real(dp), allocatable, dimension(:)   :: B1h, En, Xj, EnG, Ey, DEnr, DVnr
    real(dp), allocatable, dimension(:,:) :: X0
    
    Type(MPI_Datatype) :: mpi_type2_real

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
    !     - - - - - - - - - - - - - - - - - - - - - - - - -
    !Init mpi
    Call MPI_Init(mpierr)
    !Get process id
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    !Get number of processes
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    ! Set MPI type for type_real
    Select Case(type2_real)
    Case(sp)
        mpi_type2_real = MPI_REAL
    Case(dp)
        mpi_type2_real = MPI_DOUBLE_PRECISION
    End Select

    If (mype == 0) Then
        ! CONF_PT.RES is opened only by the master process
        open(unit=11,status='UNKNOWN',file='CONF_PT.RES')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call Input
        Call Init
        Call Rint
        Call Dinit
        Call Jterm
        Call Wdet('CONF_PT.DET')
        Call NR_Init
        Call Weight_CI
    End If
    
    Call BcastParams
    Call AllocatePTArrays
    Call BcastPTArrays

    Call cpu_time(start_time)
    Call PT_Init(npes, mype)
    Call cpu_time(stop_time)

    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    Call MPI_Bcast(En(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    If (mype == 0) Then
        print*, 'Total time for diagonal was', stop_time-start_time, 'seconds.'
    End If

    Call cpu_time(start_time)
    Call PTE(npes, mype)
    Call cpu_time(stop_time)
    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    If (mype == 0) Then
        print*, 'Total time for off-diagonal was', stop_time-start_time, 'seconds.'
        Call SaveVectors
        Call Weight_PT
        Call Cutoff(Nc_0)
        Call Sort(Nc_0)
        Call NewConfList(Nc_0)
        ! - - - - - - -  - - - - - - -  - - - - - 
        Close(11)
    End If

    Call MPI_Finalize(mpierr)

Contains

    Subroutine Input
        ! This subroutine reads in parameters and configurations from CONF.INP
        Use conf_init, Only : ReadConfInp, ReadConfigurations
        Implicit None
        Character(Len=512) :: strfmt

        ! Write name of program
        Select Case(type2_real)
        Case(sp)
            strfmt = '(/4X,"Program CONF_PT v2.3", &
               /4X,"PT corrections to binding energy", & 
               /4X,"Zero approximation is taken from CONF.XIJ", &
               /4X,"New vectors are in CONF_PT.XIJ and", &
               /4X,"new input is in CONF_new.INP")'
        Case(dp) 
            strfmt = '(/4X,"Program CONF_PT v2.3 with double precision for 2e integrals", &
               /4X,"PT corrections to binding energy", & 
               /4X,"Zero approximation is taken from CONF.XIJ", &
               /4X,"New vectors are in CONF_PT.XIJ and", &
               /4X,"new input is in CONF_new.INP")'           
        End Select
        Write( 6,strfmt) 
        Write(11,strfmt)

        average_diag=.true.      ! diagonal is averaged over relat. config.

        ! Read input parameters from file CONF.INP
        Call ReadConfInp

        Ncci = Nc
        If (Ncpt < Ncci) Then
            Write (*,*) 'Nothing to do for NcPT =', Ncpt, ' & Nc = NcCI =', Ncci
            Stop
        End If
        Nc = Ncpt

        ! Read job parameters from file pt.in
        Open(unit=21, file='cpt.in')
        Read(21,*) ktf,kvar
        Close(21)
    
        If (abs(C_is) < 1.d-6) K_is = 0
        If (K_is == 0) C_is = 0.d0
        If (K_is == 2  .or.  K_is == 4) Then
            Write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): '
            Read(*,*) K_sms
            If ((K_sms-1)*(K_sms-2)*(K_sms-3) /=  0) Stop
        End If
        ! Read configurations from file CONF.INP
        Call ReadConfigurations

        Kl = 3

    End Subroutine Input

    Subroutine Init
        Use utils, Only : DetermineRecordLength
        Implicit None
        Integer                     :: ii, ni, If, nj, i, nnj, llj, jlj, kkj, err_stat, &
                                       nec, imax, j, n, ic, i0, nmin, i1, n2, n1, l
        Real(dp)                    :: c1, c2, z1, d
        Real(dp), Dimension(IP6)    :: p, q, p1, q1
        Real(dp), Dimension(4*IP6)  :: pq
        Real(dp), Dimension(IPs)    :: Qq1
        Integer, Dimension(4*IPs)   :: IQN
        Integer, Dimension(33)      :: nnn, jjj, nqq
        Character(Len=1)            :: Let(9), lll(33)
        Character(Len=512)          :: strfmt, err_msg
        Logical                     :: longbasis, success

        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), &
                    (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        data Let/'s','p','d','f','g','h','i','k','l'/

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

        z1 = pq(1)
        If (abs(Z-z1) > 1.d-6) Then
            Write( 6,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Write(11,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Read(*,*)
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        longbasis=abs(PQ(20)-0.98765d0) < 1.d-6
        strfmt='(4X,"Kl  =",I3,7X,"Kv  =",I3, &
                    7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
                   /4X,"Nsp =",I6,4X,"Ns  =",I3,7X,"Nso =",I3, &
                    7X,"Nc =",I6)'
        Write( 6,strfmt) Kl,Kv,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,strfmt) Kl,Kv,Z,Jm,Nsp,Ns,Nso,Nc
        If (longbasis) Then
            Write( *,*) ' Using variant for long basis '
            Write(11,*) ' Using variant for long basis '
            Do ni=1,Ns
                Nn(ni)=IQN(4*ni-3)
                Ll(ni)=IQN(4*ni-2)
                Kk(ni)=IQN(4*ni-1)
                Jj(ni)=IQN(4*ni)
            End Do
        else
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
                else If (ni == ns) Then
                    Write( 6,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                    Write(11,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                    Stop
                End If
            End Do
            Nip(nj)=ni
            If (Nsu < ni) Nsu=ni
        End Do
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

        Write( 6,'(1X,71("="))')
        Write(11,'(1X,71("="))')
        Do ni=1,Nso
            l =Ll(ni)+1
            lll(ni)=let(l)
        End Do
        If (Kout > 1) Then
        Write(11,'(1X,"Core:", 6(I2,A1,"(",I1,"/2)",I2,";"),  &
                    /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                    /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                    /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                    /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                    /6X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                    (Nn(i),lll(i),Jj(i),Nq(i),i=1,Nso)
        Write(11,'(1X,71("="))')
        End If
        ! Write( 6,55)
        If (Kout > 1) Then
            Write(11,'(1X,71("="))')
        End If
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
            If (Kout > 1) Then
                Write(11,'(1X,I6,"#",6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
            End If
        End Do
        If (Kout > 1) Then
            Write(11,'(1X,71("="))')
        End If
        Do ni=Nso+1,Nsu
            Read(12,rec=2*ni+7) p
            Eps(ni)=-p(ii+1)
        End Do
        Write(11,'(" HF energies are Read from DAT file", /5(I5,F10.6))') (i,Eps(i),i=Nso+1,Nsu)
        Close(unit=12)
        Open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        Read(16) Ngaunt
        Allocate(In(Ngaunt))
        Allocate(Gnt(Ngaunt))
        Read(16) (In(i),i=1,Ngaunt)
        Read(16) (Gnt(i),i=1,Ngaunt)
        Close(unit=16)

    End Subroutine Init

    Subroutine BcastParams
        use mpi_f08
        Implicit None
        Integer :: mpierr
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Am, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(XJ_av, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngaunt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(K_sms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kdsig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kout, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kecp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(C_is, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Cut0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ncpt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nsu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngint, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        Call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        Return
    End Subroutine BcastParams

    Subroutine AllocatePTArrays
        use mpi_f08
        Implicit None
        Integer :: mpierr
        
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        If (.not. allocated(Nvc)) allocate(Nvc(Nc))
        If (.not. allocated(Nc0)) allocate(Nc0(Nc))
        If (.not. allocated(Ndc)) allocate(Ndc(Nc))
        If (.not. allocated(Ndcnr)) allocate(Ndcnr(Nc))
        If (.not. allocated(Jz)) allocate(Jz(Nst))
        If (.not. allocated(Nh)) allocate(Nh(Nst))
        If (.not. Allocated(In)) Allocate(In(Ngaunt))
        If (.not. Allocated(Gnt)) Allocate(Gnt(Ngaunt))
        If (.not. allocated(Diag)) allocate(Diag(Nd))
        If (.not. allocated(Rint1)) allocate(Rint1(Nhint))
        If (.not. allocated(Iint1)) allocate(Iint1(Nhint))
        If (.not. allocated(Iint2)) allocate(Iint2(Ngint))
        If (.not. allocated(Iint3)) allocate(Iint3(Ngint))
        If (Kbrt == 0) Then
            If (.not. allocated(Rint2)) allocate(Rint2(1,Ngint))
        Else
            If (.not. allocated(Rint2)) allocate(Rint2(2,Ngint))
        End If
        If (.not. allocated(IntOrd)) allocate(IntOrd(nrd))
        If (.not. allocated(Iarr)) allocate(Iarr(Ne,Nd))
        If (.not. allocated(DVnr)) allocate(DVnr(Nc))

        If (K_is /= 0) Then
            If (.not. allocated(R_is)) allocate(R_is(Nhint))
            If (.not. allocated(I_is)) allocate(I_is(Nhint))
        End If

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End Subroutine AllocatePTArrays

    Subroutine BcastPTArrays
        Use mpi_f08
        Use mpi_utils
    	Implicit None
        Integer :: i, mpierr
        Integer(Kind=int64) :: count

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint1, Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint1, Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (Kbrt == 0) Then
            count = Ngint
        Else
            count = Ngint*2_int64
        End If
        Call BroadcastI(Iint2, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastI(Iint3, Ngint, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IntOrd, nrd, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpierr)
        Call BroadcastD(Rint2, count, 0, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd,kind=int64)
        Call BroadcastI(Iarr, count, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(In, Ngaunt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gnt, Ngaunt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(DVnr, Nc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag, Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        If (K_is /= 0) Then
            Call MPI_Bcast(R_is, Nhint, mpi_type2_real, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(I_is, Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        Return
    End Subroutine BcastPTArrays

    Subroutine NR_Init
        Implicit None
        Integer  :: inr, nrel, ndet, nshell, k, ier, knr1, knr2, ic
        Logical  :: dIf
        Integer, Dimension(20)   ::  inl1, inl2, inq1, inq2

        allocate(Ndcnr(Nc),Nvcnr(Nc),NRR(Nc),NRN(Nc))
        inr=1  ! NR config
        nrel=0
        ndet=0
        nshell=0
        Ncnr0=0
        Ncnrci=0
        Nmax=Nd
        Ndcnr=0
        Nvcnr=0
        NRR=0
        NRN=0
        Call Squeeze(1,knr1,inl1,inq1)
        NRR(1)=1
        Do ic=1,Nc
            nrel=nrel+1                   ! number of r.configs in nr.config
            Ndcnr(inr)=Ndcnr(inr)+Ndc(ic)
            ndet=ndet+Ndc(ic)
            nshell=nshell+Nvc(ic)
            dIf=.true.
            If (ic < Nc) Then
                Call Squeeze(ic+1,knr2,inl2,inq2)
                If (knr1 == knr2) Then
                    Do k=1,knr1
                        ier=abs(inl1(k)-inl2(k))+abs(inq1(k)-inq2(k))
                        dIf=ier /=  0
                        If (dIf) goto 100
                    End Do
                End If
            End If
    100     If (dif) Then
                NRN(inr)=nrel
                Nvcnr(inr)=nshell
                If (Kout >= 1) Write (11,5) inr,nrel,NDCnr(inr),Nvcnr(inr)
    5           format(4X,'NR con-n',i5,':',i3, &
                       ' R con-ns,',i5,' det-ns &',i4,' shells')
                nrel=0
                nshell=0
                If (ic == Ncci) Then
                    Ncnrci=inr
                End If
                If (ic > Ncci  .and.  Ncnrci == 0) Then
                    Write (*,*) ' NR_init error: Ncci=',Ncci
                    Write (*,*) ' NR config Ends at ic=',ic
                    Stop             
                End If 
                If (ic /=  Nc) Then
                    inr=inr+1
                    NRR(inr)=ic+1
                End If
            End If
            knr1=knr2
            Do k=1,knr1
                inl1(k)=inl2(k)
                inq1(k)=inq2(k)
            End Do
        End Do
        Ncnr=inr
        Write ( *,*) ' PT space:',Ncnr,' non-rel. configurations'
        Write (11,*) ' PT space:',Ncnr,' non-rel. configurations'
        Write ( *,*) ' CI space:',Ncnrci,' non-rel. configurations'
        Write (11,*) ' CI space:',Ncnrci,' non-rel. configurations'
        Return
    End Subroutine NR_Init

    ! =============================================================================
    Subroutine Weight_CI
        Implicit None
        Integer  :: i, ic, kcnr, ndk, k, j, kc, icnr, ncnr0, nc0, icci, nc0j, &
                   ncnr0j, nd0j, ndci
        Real(dp) :: wj, wmax, d, dvnrn, dvnrx, dummy, wj0
        Logical  :: tail
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Open(unit=17,file='CONF.XIJ', &
                 status='OLD',form='UNFORMATTED',err=900)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(DVnr(Nc),B1h(Nd))
        Nd0=0    
        nc0=0    
        ncnr0=0  
        wmax=0.d0
        Do j=1,Nlv
            Read (17) d,dummy,ndci,(B1h(i),i=1,ndci)
            wj=1.d0                 ! norm of WF
            i=0
            ic=0
            tail=.false.
            Do icnr=1,Ncnrci    !### loop over NR conf-ns
                kcnr=NRN(icnr)
                d=0.d0
                ndk=0
                Do kc=1,kcnr          !### loop over R conf-ns
                    ic=ic+1
                    ndk=ndk+Ndc(ic)     !### ndk= number of det-s in NR con-n
                End Do
                Do k=1,ndk
                    i=i+1
                    d=d+B1h(i)**2
                End Do
                If (j == 1) Then
                    DVnr(icnr)=d
                else
                    DVnr(icnr)=DVnr(icnr)+d
                End If
                wj=wj-d                ! norm of the remaining tail of WF
                If ((.not.tail)  .and.  i < IPPT) Then
                    nd0j=i
                    nc0j=ic  
                    ncnr0j=icnr
                    wj0=wj
                End If
                tail=wj < Cut0
            End Do
            Write(*,5) j,nd0j,nc0j,ncnr0j,wj0
    5       format(4x,'WF ',i2,' Nd0=',i6,' Nc0=',i5,' Ncnr0=',i5, &
                 ' wj=',f9.6)
            If (wj0 > wmax) wmax=wj0
     
            If (Nd0 < nd0j) Then
                Nd0=nd0j
                nc0=nc0j  
                ncnr0=ncnr0j
            End If            
        End Do
        Close(unit=17)
        Ncp0=nc0
        Write( *,15) Nd0,Ncp0,ncnr0,wmax,Cut0
        Write(11,15) Nd0,Ncp0,ncnr0,wmax,Cut0
    15  format(/4x,'PT block: Nd0=',i6,' Nc0=',i5,' Ncnr0=',i5, &
               /4x,'Actual Cutoff=',f9.6,' Cut0=',f9.6)
        
        If (Kout >= 1) Write(11,*) ' Weights of all NR configurations:'
        icci=0
        dvnrn=11.d1
        Do ic=1,Ncnr
            If (ic <= Ncnrci  .and.  DVnr(ic) < dvnrx) icci=icci+1
            If (ic <= Ncnrci  .and.  DVnr(ic) < dvnrn) dvnrn=DVnr(ic)
            If (Kout >= 1) Then
                Write(11,25) ic,DVnr(ic)
    25          format(4x,'NR conf ',i4,' weight ',f10.7)
                If (ic == Ncnrci) &
                  Write(11,*) '   ---- End of CI space ---------'
            End If
        End Do
        Return
    900 Write (*,*) ' CONF.XIJ file with WFs is absent'
        Stop
    End Subroutine Weight_CI
    ! =============================================================================
    Subroutine PT_Init(npes, mype)
    	use mpi_f08
        Implicit None
        Integer  :: i, ni, j, n, k, kx, ic, kd, ndci
        Integer :: npes, mype, mpierr
        Real(dp)   :: sd, x
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (mype==0) Then ! master pe Only
            Write ( *,5)
            Write (11,5)
    5       format(/4X,'PT_init')
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            Open(unit=16,file='CONF.XIJ',status='OLD', &
                 form='UNFORMATTED',err=900)
            Read(16) x,x,Nd1
            kd=0
            Nd0=0
            Do i=1,Nc
              kd=kd+Ndc(i)
              If (i == Ncci) ndci=kd
              If (i == Ncp0) Nd0=kd
            End Do
            If (ndci /=  Nd1) Then
              Write (*,*) ' Wrong length of vectors in CONF.XIJ:'
              Write (*,*) ' Nd1 =',Nd1,'; expected ',ndci
              Stop
            End If
        
            Write ( *,25) Ncci,Nd1,Ncp0,Nd0,Nc,Nd
            Write (11,25) Ncci,Nd1,Ncp0,Nd0,Nc,Nd
    25      format(4X,'CI block: Ncci=',I6,' Nd1=',I8,/4X, &
                  'PT block: (Ncp0=',I6,' Nd0=',I8,')*(Nc=',I7,' Nd=',I9,')')
            If (Nd0 > Nd1) Then
              Write (*,*) ' can not work with Nd0 > Nd1'
              Stop
            End If
            If (Nd0 > IPPT) Then
              Write (*,*) ' can not work with Nd0 > IPPT =',IPPT
              Stop
            End If
            If (Nd0 == 0) Then
              Write (*,*) ' Nothing to Do for Nd0=0 (Ncp0=',Ncp0,')'
              Stop
            End If
        
            Nlv=min(Nlv,Nd0)
        End If ! End master pe Only

        Call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        allocate(X0(Nd,Nlv), En(Nlv), EnG(Nlv), Xj(Nlv))

        Call DiagH(Nd1,npes,mype)      !# forms array Diag(i)=H(i,i)

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        If (mype==0) Then ! master pe Only
            If (average_diag) Then
              n=Nd1
              Do ic=Ncci+1,Nc      !# averagin diagonal for each configuration
                kx=Ndc(ic)
                sd=0.d0
                If (kx /=  0) Then
                  Do k=1,kx
                    sd=sd+Diag(n+k)
                  End Do
                  sd=sd/kx
                  Do k=1,kx
                    Diag(n+kx)=sd
                  End Do
                  n=n+kx
                End If
              End Do
              If (n /=  Nd) Then
                Write(*,*) 'PT_init error: n=',n,' must be equal to Nd=',Nd
                Stop
              else
                Write( *,*) ' Diagonal averaged over relat. config.'
                Write(11,*) ' Diagonal averaged over relat. config.'
              End If
            End If
            rewind(16)           !## starting from Nd1+1 to Nd
            Do n=1,Nlv
              Read(16,err=220,End=220) EnG(n),Xj(n),Nd1,(X0(j,n),j=1,Nd1)
              En(n)=EnG(n)+4.d0*Gj*Xj(n)*(Xj(n)+1.d0)
              ni=n
            End Do
    220     Close(16)
            Write (*,35) ni
    35      format (4X,i5,' vectors Read from file CONF.XIJ.')
            Nlv=ni
        End If
        Return
    900 Write(*,*)' File CONF.XIJ is absent.'
        Stop
    End Subroutine PT_Init

    Subroutine PTE(npes,mype)
        use mpi_f08
        Implicit None
        Integer  :: n, ic, l, i, k, ix
        Integer, Dimension(Nc)  :: Ndic
        Integer  :: npes, mype, mpierr, interval, remainder, start, End, size
        Integer  :: istart, iEnd, lastsize
        Real(dp)  :: ei, x, e0, del, ci, x2, des, des0, dEx, dVs
        Real(dp), Dimension(IPPT)  :: E1
        Real(dp), Dimension(Nlv,Nc) :: dE, dV, dE1, dV1
        Integer, Dimension(2,Nc) :: sizeNcNd
        Integer     :: totalNcNd
        Integer, Dimension(npes) :: sizes, idisp

        allocate(Ey(Nlv),Ndirc(Nc),DEnr(Nc))

        If (mype == 0) Then
        Write ( *,5)
        Write (11,5)
    5   format(/4X,'Second order corrections to the energies', &
                  ' & first order eigenvectors')
        
        E1=0.d0
        Do n=1,Nlv
            Ey(n)=0.d0
        End Do
        
        Do ic=Ncnrci+1,Ncnr
            DEnr(ic)=0.d0
            DVnr(ic)=0.d0
        End Do

        Ndirc(1)= 0
        Do ic= 2,Nc
            Ndirc(ic) = Ndirc(ic-1) + Ndc(ic-1)
        End Do
    
        Ndic=0
        Do ic= 1,Ncnr
            If (ic == 1) Then
                Ndic(1)= 0
            else
                Ndic(ic)= Ndic(ic-1)+Ndcnr(ic-1)
            End If
            Do n=1,Nlv
                dE(n,ic)=0.e0
                dV(n,ic)=0.e0
            End Do
        End Do
        End If

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ncnr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ncnrci, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndcnr, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ey, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(En, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndirc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndic, Ncnr, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(X0(1:Nd1,1:Nlv), Nd1*Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag, Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        totalNcNd = 0

        Do ic=Ncnrci+1,Ncnr
            Do l=1,Ndcnr(ic)
                totalNcNd = totalNcNd + 1
                sizeNcNd(1, ic) = ic
                sizeNcNd(2, ic) = Ndcnr(ic) ! size of each config
            End Do
        End Do
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Do ic=Ncnrci+1,Ncnr ! loop over configs in PT space

            interval = sizeNcNd(2,ic)/npes
            remainder = mod(sizeNcNd(2,ic),npes)
            !print*, 'starting ic=', ic, 'of ', Ncnr
            Do i=0,npes-1
                If (mype == npes-1) Then
                    start = 1+mype*interval
                    End = (mype+1)*interval+remainder
                    istart = Ndic(ic)+start
                    iEnd = Ndic(ic)+End
                    lastsize= End - start + 1
                    !print*,ic,mype,start,End,sizeNcNd(2,ic)
                else If (mype /= npes-1) Then
                    start = 1+mype*interval
                    End = (mype+1)*interval
                    istart = Ndic(ic)+start
                    iEnd = Ndic(ic)+End
                    !print*,ic,mype,start,End,sizeNcNd(2,ic)
                    !print*, 'rank ', mype, 'has istart=', istart, 'iEnd=', iEnd, 'and size = ', size
                End If
                size = End - start + 1
            End Do  

            Do i=1,npes-1
                sizes(i) = interval
            End Do
            
            sizes(npes) = interval + remainder

            idisp(1)=0
            Do i=2, npes
                idisp(i) = (i-1)*sizes(1)
            End Do
            Do l=start, End  ! loop over dets in configs in PT space
                i= Ndic(ic)+l
                ei= -Diag(i)
                Call FormH(i,Nd0,E1)    !# evalulation of i-th MEs of ALL
                Do n=1,Nlv              !## first order correction vectors
                    x=0.d0                !### and writing them to VEC.TMP,
                    e0=En(n)
                    Do k=1,Nd0
                        x=x+E1(k)*X0(k,n)
                    End Do
                    x2=x*x
                    del=e0-ei
                    ci=x/del
                    dE(n,ic)= dE(n,ic)+x2/del    ! contribution of one NR config.
                    dV(n,ic)= dV(n,ic)+(ci)**2
                    X0(i,n)=ci                   ! first order eigenvector
                End Do
            End Do
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Do n=1,Nlv
                Call MPI_Allreduce(dE(n,ic),dE1(n,ic), 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, MPI_COMM_WORLD, mpierr)
                Call MPI_Allreduce(dV(n,ic),dV1(n,ic), 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, MPI_COMM_WORLD, mpierr)
                Call MPI_Gatherv(X0(istart:iEnd,n), size, MPI_DOUBLE_PRECISION, &
                               X0(Ndic(ic)+1:Ndic(ic)+Ndcnr(ic),n), &
                               sizes, idisp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            End Do

            dEx=0.d0
            dVs=0.d0
            dE=dE1
            dV=dV1
            Do n=1,Nlv
                Ey(n)= Ey(n)+ dE(n,ic)    ! sum of all contributions to n's energy
                If (dE(n,ic) > dEx) Then
                    dEx= dE(n,ic)
                End If
                dVs= dVs+dV(n,ic)
            End Do
            DEnr(ic)= dEx
            DVnr(ic)= dVs
            If (mype==0) Then
                print*, 'finished ', ic, 'of ', Ncnr
            End If

        End Do
    
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        If (mype == 0) Then
        Write ( *,25) NcCI,Nd1,Ncp0,Nd0,Nc,Nd
        Write (11,25) NcCI,Nd1,Ncp0,Nd0,Nc,Nd
    25  format(/82('='), &
             /' CI space [',i6,'/',i8,']; PT space [',i5,'/',i6, &
             '] * [',i7,'/',i9,']' &
             /' N',6X,'J',12X,'Ev0',6X,'DEL(CM**-1)',7X,'dE(PT2)', &
             6X,'Ev(PT2)',3X,'DEL(CM**-1)',/82('-'))
    
        Do n=1,Nlv
            des0=(En(1)-En(n))*2*DPRy
            des =des0 + (Ey(1)-Ey(n))*2*DPRy
            Write ( *,55) n,Xj(n),En(n),des0,Ey(n),En(n)+Ey(n),des
            Write (11,55) n,Xj(n),En(n),des0,Ey(n),En(n)+Ey(n),des
    55      format(I2,F12.8,F14.8,F13.2,2F14.8,F13.2)
        End Do
        Write ( *,65)
        Write (11,65)
    65  format(82('='))
    
        If (Kout >= 1) Write (11,75)
    75  format(/4x,'PT contributions of non-relat. configurations', &
               /4x,'ic     max dE   Full weight',/4x,27('-'))
    
        dvnrx=0.d0
        ix=0
        Do ic=Ncnrci+1,Ncnr
            If (DVnr(ic) > dvnrx) Then
                dvnrx=DVnr(ic)
                ix=ic
            End If
            If (Kout >= 1) Write(11,85) ic,DEnr(ic),DVnr(ic)
    85      format(i6,2f12.7)
        End Do

        If (Kout >= 1) Write (11,95)
    95  format(4x,27('-'))
        Write ( *,105) dvnrx,ix
        Write (11,105) dvnrx,ix
    105 format('Max PT weight ',e11.3,' of config. ',i6)
        End If

        Return
    End Subroutine PTE
    ! =============================================================================
    Subroutine SaveVectors
        Implicit None
        Integer  :: n, i, j
        Real(dp)   :: s
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Write ( *,5)
        Write (11,5)
    5   format(/4X,'Normalizing & saving vectors')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Open(unit=16,file='CONF_PT.XIJ',status='UNKNOWN', &
             form='UNFORMATTED')
        Close(unit=16,status='DELETE')
        Open(unit=16,file='CONF_PT.XIJ',status='NEW', &
             form='UNFORMATTED')
        Do n=1,Nlv
            s=0.d0
            Do i=1,Nd
                s=s+X0(i,n)**2
            End Do
            s=1.d0/sqrt(s)
            Do i=1,Nd
                !print*,i,s,X0(i,n)
                X0(i,n)=s*X0(i,n)
            End Do
            EnG(n)=En(n)-4.d0*Gj*Xj(n)*(Xj(n)+1.d0)
            Write (16) EnG(n),Xj(n),Nd,(X0(j,n),j=1,Nd)
            Write( *,15) n,EnG(n),Xj(n),s
            Write(11,15) n,EnG(n),Xj(n),s
    15      format(4x,'vector ',i2,' E_G=',f13.8,' J=',f8.5,' S_norm=', &
                 f8.5,' saved')
        End Do
        Close(16)
        Write (*,35) 
    35  format (4X,' vectors saved to file CONF_PT.XIJ.')
        Return
    End Subroutine SaveVectors
    ! =============================================================================
    Subroutine Weight_PT
        Implicit None
        Integer  :: icci, ic, icpt
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        icci=0
        dvnrn=11.d1
        Do ic=1,Ncnrci
          If (DVnr(ic) < dvnrx) icci=icci+1
          If (DVnr(ic) < dvnrn) dvnrn=DVnr(ic)
        End Do
    
        icpt=0
        Do ic=Ncnrci+1,Ncnr
          If (DVnr(ic) > dvnrn) icpt=icpt+1
        End Do
        Write ( *,15) icci,Ncnrci,dvnrx,icpt,Ncnr,dvnrn
        Write (11,15) icci,Ncnrci,dvnrx,icpt,Ncnr,dvnrn
    15  format('Weights of ',i6,' (from ',i6, &
             ') NR config-s in CI below max PT w-t',e11.3, &
             /'Weights of ',i6,' (from ',i6, &
             ') config-s in PT space above min CI w-t',e11.3)
        Return
    End Subroutine Weight_PT
    ! =============================================================================
    Subroutine Cutoff(Nc_0)
        Implicit None
        Integer  :: izer, icc, icnr, k, ii, ic, iinr, last, i, &
                    Nc_0
        Integer, Dimension(15)  :: Nconf, Mconf
        Real(dp)  :: wx, wl, xx, dlx, dln, ctf, xcc
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(W(Nc,2))
        izer=0
        Do icc=1,15
           Nconf(icc)=0
           Mconf(icc)=0
        End Do
        Do icnr=1,Ncnr                     !### transforms weights to
          wx=DVnr(icnr)                    !###  log scale
          If (wx > 0.d0) Then
            wl=dlog(wx)*0.4342945d0        !### (factor = lg e)
          else
            izer=izer+1
            wl=-99.d0
          End If
          W(icnr,1)=wl
          Do icc=1,15
            xcc=2-icc
            If (wl < xcc) Nconf(icc)=Nconf(icc)+1
          End Do
          xx=dvnrx
          Do k=1,10
            If (wx > xx) Mconf(k)=Mconf(k)+1
            xx=0.5d0*xx
          End Do
        End Do
        Do icc=1,15
          If (Nconf(icc) /= 0) last=icc
        End Do
        If (Kout >= 1) Write(11,5) izer, &
            (Ncnr-Nconf(i),Nconf(i),2-i,i=1,last)
    5   format(4X,I6,' nonrel. conf-s with zero weight', &
              /(4X,I6,' weight above and ',I6,' below 10**(',I3,')'))
    
        Write(*,15) (Mconf(k),dvnrx,k-1,dvnrx/2**(k-1),k=1,10)
        Write(11,15) (Mconf(k),dvnrx,k-1,dvnrx/2**(k-1),k=1,10)
    15  format(/(4x,i6,' weights above ',d12.3,'/2**',i2,' = ',e11.3))

        Write(*,*) ' Give cutoff value of k (0-9):'
        dlx=dlog(dvnrx)*0.4342945d0
        dln=dlog(dvnrn)*0.4342945d0
        ctf=dlx - 0.301029996d0*ktf
    
        Write(*,*) 'Choose reordering variant:'
        Write(*,*) '(1) - all; (2) - keep PT block; (3) - keep CI space'
        Nc_0=0
        If (kvar == 2) Nc_0=Ncnr0
        If (kvar == 3) Nc_0=Ncnrci
    
        ii=0
        ic=0
        iinr=0
        Do icnr=1,Ncnr
          ic=ic+NRN(icnr)
          If (W(icnr,1) > ctf .or. icnr <= Nc_0) Then
            iinr=iinr+1
            ii=ii+NRN(icnr)
            W(icnr,2)=1.d0
          else
            W(icnr,2)=0.d0
          End If
        End Do
        If (ic /= Nc) Then
          Write(*,*) ' Cutoff: ic=',ic,' is not equal Nc=',Nc
        End If
        Ncnew=ii
        Write ( *,105) 10.d0**ctf,iinr,Ncnrci,Ncnew,Ncci
    105 format(4X,'For cutoff ',E10.3,' CI: Ncnr=',I6,' (was ',I6, &
             ') and Nc =',I6,' was(',I6,')')
        Return
    End Subroutine Cutoff
    ! =============================================================================
    Subroutine Sort(Nc_0)
        Implicit None
        Integer  :: ic, i1, irpl, mm, k, kc, Nc_0
        Real(dp)  :: x, y
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        irpl=0                    !### counts replacements
        Do ic=Nc_0+1,Ncnr-1
          x=-100.d0
          Do kc=ic,Ncnr           !### search for max from the rest
            y=W(kc,1)
            If (y > x) Then
              i1=kc
              x=y
            End If
          End Do
          If (ic /= i1) Then      !### placing max in front
            irpl=irpl+1
            Do k=1,2
              x=W(i1,k)
              W(i1,k)=W(ic,k)
              W(ic,k)=x
            End Do
            mm=NRN(i1)
            NRN(i1)=NRN(ic)
            NRN(ic)=mm
            mm=NRR(i1)
            NRR(i1)=NRR(ic)
            NRR(ic)=mm
            mm=Nvcnr(i1)
            Nvcnr(i1)=Nvcnr(ic)
            Nvcnr(ic)=mm
          End If
        End Do
        Write(*,*) irpl,' replacements first'
        Return
    End Subroutine Sort
    ! =============================================================================
    Subroutine NewConfList(Nc_0)
        Implicit None
        Character(Len=1), Dimension(16)  :: name
        Character(Len=1), Dimension(5)   :: ch1
        Character(Len=1)  :: str*70
        Integer  :: Nc_0, icnr, ic0, ic1, kd, kd4, i, num, numt, nci, nrci, &
                    icnr1, num1, n, n1, n2, kc, kcnr, kc4, ic
        Logical  :: inci
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Write( *,*) ' Forming new list of configurations...'
        Close(11)
        Open(unit=11,file='CONF_new.INP',status='UNKNOWN')
        Close(11,status='DELETE')
        Open(unit=11,file='CONF_new.INP',status='NEW')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (Nc_0 < Nc4) Then           ! we need to find new Nc4 
          kd4=0
          kc4=0
          Do icnr=1,Ncnr
            ic0=NRR(icnr)
            ic1=ic0+NRN(icnr)-1
            kd=0
            Do ic=ic0,ic1
              kd=kd+Ndc(ic)
            End Do
            If (kd4+kd > IP1) goto 100
            kd4=kd4+kd
            kc4=kc4+NRN(icnr)
          End Do
    100   If (kc4 > 0) Then 
            Nc4=kc4
            else
            Nc4=NRN(1)
          End If
          Write(*,*)' New Nc4=',Nc4,' (Nd4=',kd4,')'         
        End If
    
    !   input from the file 'CONF.INP'
        Open(unit=10,file='CONF.INP',status='OLD')
        Read (10,5) name
    5   format(1X,16A1)
        Write(11,'(1X,16A1,A12)') name,' (reordered)'
        Read (10,15) (ch1(i),i=1,5),Z
        Write(11,15) (ch1(i),i=1,5),Z
        Read (10,15) (ch1(i),i=1,5),Am
        Write(11,15) (ch1(i),i=1,5),Am
        Read (10,15) (ch1(i),i=1,5),XJ_av
        Write(11,15) (ch1(i),i=1,5),XJ_av
        Read (10,15) (ch1(i),i=1,5),Jm
        Write(11,15) (ch1(i),i=1,5),Jm
    15  format (5A1,F5.1)
        Read (10,25) (ch1(i),i=1,5),Nso
        Write(11,25) (ch1(i),i=1,5),Nso
        Read (10,25) (ch1(i),i=1,5),Nc
        Write(11,25) (ch1(i),i=1,5),Ncnew
        Read (10,25) (ch1(i),i=1,5),Kv
        Write(11,25) (ch1(i),i=1,5),Kv
        Read (10,25) (ch1(i),i=1,5),Nlv
        Write(11,25) (ch1(i),i=1,5),Nlv
        Read (10,25) (ch1(i),i=1,5),Ne
        Write(11,25) (ch1(i),i=1,5),Ne
    25  format (5A1,I6)
    
    200 Read (10,35) (ch1(i),i=1,5),str
    35  format(5A1,A)
        If (ch1(2) /= 'K'.AND.ch1(2) /= 'k') goto 300
        If (ch1(3) /= 'L'.AND.ch1(3) /= 'l') goto 300
        If (ch1(4) /= '4') goto 300
        If (Nc_0 == Ncnrci) Then
          kl4=2            ! initial approximation from disk
        else
          kl4=1
        End If
        Write(11,25) (ch1(i),i=1,5),kl4
        goto 200    
    
    300 If (ch1(2) /= 'N'.AND.ch1(2) /= 'n') goto 400
        If (ch1(3) /= 'C'.AND.ch1(3) /= 'c') goto 400
        If (ch1(4) /= '4') goto 400
        Write(11,25) (ch1(i),i=1,5),Nc4
        goto 200
        
    400 If(ch1(1) == ' '.AND.ch1(2) == ' ') goto 500
        Write(11,35) (ch1(i),i=1,5),str
        goto 200
    500 Close(10)
        Write(11,45) (Qnl(i),i=1,Nso)
    45  format (6(4X,F7.4))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        num=0
        numt=0
        nci=0
        nrci=0
        inci=.TRUE.
        Do icnr=1,Ncnr
          kcnr=NRN(icnr)
          numt=numt+kcnr
          If (W(icnr,2) /= 0.d0) Then
            nrci=nrci+1
          else
            If (inci) Then
              Write(11,55)
    55        format('<CI<')
              inci=.FALSE.
            End If
          End If
          ic=NRR(icnr)-1
          icnr1=mod(icnr,10000)
          Write(11,65) icnr1,10**W(icnr,1),Nvcnr(icnr)
    65    format(I4,62X,E9.2,i6)
          Do kc=1,kcnr
            ic=ic+1        !### old index of current configuration
            num=num+1      !### new index of current configuration
            If (inci) nci=nci+1
            num1=mod(num,10000)
            n1=Nc0(ic)+1
            n2=Nc0(ic)+Nvc(ic)
            Write(11,75) num1,(Qnl(n),n=n1,n2)
    75      format(I4,6(F7.4,4X),/6(4X,F7.4))
          End Do
        End Do
        If (nci == Ncnew) Then
          Write(11,85)
    85    format(1X,'>>>>>>>>>> END <<<<<<<<<<')
          Close(11)
          Write(*,*)' CONF_new.INP is formed.'
          Write(*,*)' New CI space:',nci,' PT space:',num,' config.'
        else
          Write(*,95) nci,Ncnew
    95    format(' Error: nci=',I6,' not equal to Ncnew=',I6)
          Read(*,*)
          Return
        End If
        Return
    End Subroutine NewConfList

    Subroutine Squeeze(ic,knr,inl,inq)
        Implicit None
        Integer  :: i1, i2, knr, i, ic, nl, mq, nl0, mq0
        Real(dp)  :: x
        Integer, Dimension(20)  :: inl, inq

        i1=Nc0(ic)
        i2=i1+Nvc(ic)
        knr=0            
        nl0=0
    
        i=i1+1
        x=100*abs(Qnl(i))
        nl=x+1.d-5
        mq=100*(x-nl)+1.d-5
    
    100 If (knr >= 20) Then
            Write (*,*) ' Squeeze error for ic=',ic
            Write (*,*) ' array overflow. knr=',knr
            Read (*,*)
        End If
        If (nl == nl0) Then
            knr=knr+1
            inl(knr)=nl0
            inq(knr)=mq0+mq
            nl=0
            If (i == i2) goto 200
        else
            If (nl0 /= 0) Then
                knr=knr+1
                inl(knr)=nl0
                inq(knr)=mq0
            End If
        End If
        If (i == i2) Then
            knr=knr+1
            inl(knr)=nl
            inq(knr)=mq
            goto 200
        End If
    
        nl0=nl
        mq0=mq
        i=i+1
        x=100*abs(Qnl(i))
        nl=x+1.d-5
        mq=100*(x-nl)+1.d-5
        goto 100
    
    200 continue
        Return
    End Subroutine Squeeze

    Subroutine FormH(n,nd0,E1)
        Use determinants, Only : CompC, Gdet
        Implicit None
        Integer  :: ic, kx, k, n, nd0, icomp
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Real(dp), Dimension(IPPT)  :: E1
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (.not. allocated(idet1)) allocate(idet1(Ne))
        If (.not. allocated(idet2)) allocate(idet2(Ne))
        If (.not. allocated(iconf1)) allocate(iconf1(Ne))
        If (.not. allocated(iconf2)) allocate(iconf2(Ne))

        E1=0.d0
        ! calculation of the matrix elements
        Call Gdet(n,idet1)
        Do ic=1,Nc
            kx= Ndc(ic)
            k = Ndirc(ic)
            If (k+kx > nd0) kx=nd0-k
            If (kx /= 0) Then
                Call Gdet(k+1,idet2)
                Call CompC(idet1,idet2,icomp)
                If (icomp <= 2) Then
                    Do k= Ndirc(ic)+1,Ndirc(ic)+kx
                        Call Gdet(k,idet2)
                        E1(k)= Hmltn(idet1,idet2)
                    End Do
                End If
            End If
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Return
    End Subroutine FormH

    Real(dp) function Hmltn(idet1,idet2)
        Use integrals, Only : Gint, Hint
        Use determinants, Only : Rspq
        Implicit None
        Integer  :: nf, iq, i1, i2, j1, j2, is, jq, jq0
        Integer, Allocatable, Dimension(:)  :: idet1, idet2
        Real(dp)  :: t
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call Rspq(idet1,idet2,is,nf,i1,i2,j1,j2)
        t=0.d0
        If (nf <= 2) Then
            If (nf == 2) Then
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants dIffer by two functions
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
                t=t+Gint(i2,j2,i1,j1)*is !### det_k goes first!
            else If (nf == 1) Then
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants dIffer by one function
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
                Do iq=1,Ne
                    i1=idet1(iq)
                    If (i1 /= j1) Then
                        t=t+Gint(j2,i1,j1,i1)*is
                    End If
                End Do
                t=t+Hint(j2,j1)*is
            else
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants are equal
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
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
            End If
        End If
        Hmltn=t
        Return
    End function Hmltn

    Subroutine DiagH(nd0,npes,mype)
    	use mpi_f08
        Use determinants, Only : Gdet
        Implicit None
        Integer   :: n, i, nd0, start, End, j
        Integer   :: npes, mype, interval, remainder, mpierr, size, count
        Integer, Allocatable, Dimension(:)  :: idet1, idet2
        Integer, Dimension(npes) :: disp, sizes, Ends
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (.not.allocated(idet1)) allocate(idet1(Ne))
        If (.not.allocated(idet2)) allocate(idet2(Ne))
        Do n=1,Nd
            Diag(n)=0.d0
        End Do
        Call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        interval = (Nd-(nd0+1))/npes
        remainder = mod(Nd-(nd0+1),npes)
        count = 0
        ! setting up sizes and displacements for MPI
        Do i=0,npes-1
            If (i == npes-1) Then
                start = nd0+1+i*interval
                End = nd0+1+(i+1)*interval+remainder
            else
                start = nd0+1+i*interval
                End = nd0+(i+1)*interval
            End If
            size = End - start + 1
            sizes(i+1) = size
            Ends(i+1) = End
        End Do
        Do i=0,npes-1
            If (mype == npes-1) Then
                start = nd0+1+mype*interval
                End = nd0+1+(mype+1)*interval+remainder
            else If (mype /= npes-1) Then
                start = nd0+1+mype*interval
                End = nd0+(mype+1)*interval
            End If
            size = End - start + 1
            sizes(mype+1) = size
            Ends(mype+1) = End
        End Do        
        disp(1)=0
        Do i=2, npes
            disp(i) = Ends(i-1) - nd0
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        !   calculation of the matrix elements
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (mype == 0) Write (*,*) ' Formation of diagonal'

        Do n=start, End
            Call Gdet(n,idet1)
            Do i=1,Ne
                idet2(i)=idet1(i)
            End Do
            Diag(n)=Hmltn(idet1,idet2)
            count=count+1
            If (mype==0) Then
                Do j=1,10
                    If (count==size/j) Then
                        print*, 'diag is', (10-j)*10, '% Done'
                    End If
                End Do
            End If
        End Do
        
        !print*,'pe=',mype,'has finished working with count=', count

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        
        Call MPI_Gatherv(Diag(start:End), size, MPI_DOUBLE_PRECISION, Diag(nd0+1:Nd), &
                           sizes, disp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        If (mype==0) Then
        Write (*,*) ' Diagonal formed'
    	End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Return
    End Subroutine DiagH
End Program conf_pt

