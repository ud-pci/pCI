program formy
    use params, IP1conf => IP1
    use determinants, only : Dinit, Jterm
    use str_fmt, only : startTimer, stopTimer, FormattedTime
    use mpi_f08
    implicit none

    integer :: mype, npes, mpierr
    integer(kind=int64) :: start_time
    character(len=16) :: timeStr

    Integer :: i, n, k, l, Nd2, Nddir, nsu2, icyc, nlamb, kIters, N_it4
    Integer  :: Kli, Klf, Ndir, Int_err, Ntr, Nint, Nmax, Nd0, ipmr, Kdiag, Nlft, IP1, IPad
    Integer  :: IPlv, IP4, Kt, Nlev, N0, N2, nrange
    Integer(kind=int64) :: NumH, NumJ
    Real(dp) :: Jm0, E0, E2, Tj0, Tj2, xlamb, xlamb1, xlamb2, xlambstep, XInuc, Crit1, W0, W00, Q, Elft, Hmin, dlamb
    Real(dp), Allocatable, Dimension(:) :: xlamb1s, xlamb2s, xlambsteps, xlamblist
    Integer, Allocatable, Dimension(:) :: Inte
    Integer, Allocatable, Dimension(:,:) :: Iarr0
    Real(dp), Allocatable, Dimension(:) :: Z1, X0, X1, X2, YY1, YY2, Rnt, Ev, Diag, E1
    Real(dp), Allocatable, Dimension(:,:) :: X1J, Y2J
    Real(dp), Dimension(2) :: s, ss, s0, s1, s2
    Character(Len=1), Dimension(9)          :: Let
    Character(Len=4), Dimension(13)         :: Alet
    Character(Len=4), Dimension(5)          :: Blet
    logical :: ok

    call MPI_Init(mpierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

    call startTimer(start_time)

    if (mype == 0) Then
        call SetParams                    ! Set job and array parameters
        call Input                        ! Read list of configurations from CONF.INP
        call Init                         ! Read basis set information from CONF.DAT
        call Rint                         ! Read radial integrals from CONF.INT
        call Dinit                        ! Form list of determinants
        call Jterm                        ! Print table with numbers of levels with given J
        call Init0                        !### Evaluation of the RHS of
    end if

    call BcastVector(mype)
    call Vector(kl,mype)                   !###  the equation and vectors Yi

    call stopTimer(start_time, timeStr)
    If (mype == 0) write(*,'(2X,A)') 'TIMING >>> Total computation time of formy was '// trim(timeStr)

    call MPI_Finalize(mpierr)

contains

    Subroutine BcastVector(mype)
        use mpi_f08
        use mpi_utils, only : BroadcastI
        Implicit None

        Integer :: mype, mpierr, i0
        Integer(kind=int64) :: count

        call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jm0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kli, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Klf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        if (mype == 0) i0 = sizeof(Nh)/4
        call MPI_Bcast(i0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        if (.not. allocated(Nh)) allocate(Nh(i0))
        if (.not. allocated(Jz)) allocate(Jz(i0))
        call MPI_Bcast(Nh, i0, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jz, i0, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nn, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jj, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ll, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        if (.not. allocated(Ndc)) allocate(Ndc(Nc))
        call MPI_Bcast(Ndc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        if (.not.allocated(Rnt)) allocate(Rnt(Nint))
        if (.not. allocated(Inte)) allocate(Inte(Nint))
        call MPI_Bcast(Rnt, Nint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Inte, Nint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        if (.not. allocated(X0)) allocate(X0(Nd0))
        if (.not. allocated(X2)) allocate(X2(Nd0))
        call MPI_Bcast(X0, Nd0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(X2, Nd0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        
        if (.not. allocated(Iarr)) allocate(Iarr(Ne,Nd))   
        if (.not. allocated(Iarr0)) allocate(Iarr0(Ne,Nd0))
        count = Ne*Int(Nd,kind=int64)
        call BroadcastI(Iarr, count, 0, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd0,kind=int64)
        call BroadcastI(Iarr0, count, 0, 0, MPI_COMM_WORLD, mpierr)

    End Subroutine BcastVector

    Subroutine Vector(kl,mype)   
        ! Read RHS vector |Y1>=A|X0> and LHS vector <Y2|=<X2|B,
        Use conf_variables, Only : iconf1, iconf2
        Use determinants, Only : Gdet, Rspq, CompCD, Rdet
        Implicit None
        Integer :: kl, i, ic, icomp, is, i0, i1, j0, j1, negl, idif, n, k, kx, err_stat, nsu1, mype
        Integer :: mdel, nskip, k1, nf, iq, start, end1
        Integer, Allocatable, Dimension(:) :: idet0, idet1
        Real(dp) :: trd, xn0, xn2, s1, s2
        Character(Len=256) :: strfmt, err_msg

        if (mype == 0) then
            strfmt = '(3X,63("="),/4X,"Vector: forming vectors Y1 and Y2"," Threshold: ",E9.1)'
            If (kl.NE.2) Write( 6,strfmt) trd
        end if

        trd=1.d-10                  !### used to drop small contributions
        
        If (.not. Allocated(YY1)) Allocate(YY1(Nd))
        If (.not. Allocated(YY2)) Allocate(YY2(Nd))
        If (.not. Allocated(idet0)) Allocate(idet0(Ne))
        If (.not. Allocated(idet1)) Allocate(idet1(Ne))
        If (.not. Allocated(iconf1)) Allocate(iconf1(Ne))
        If (.not. Allocated(iconf2)) Allocate(iconf2(Ne))

        YY1(:)=0.d0
        YY2(:)=0.d0
        negl=0              !### number of vector components below threshold
        idif= 1
        If (Kli.EQ.5) idif= 2

        If (mype == 0) Then
            idet0(1:Ne)=Iarr0(1:Ne,1)
            Call DefJz2(idet0)
            Q=Jm-Jm0
            mdel=dabs(Q)+1.d-5
            If ((Kli-1)*(Klf-1).EQ.0.AND.mdel.NE.0) Then
                Write(*,*) 'Vector: |M-M0| is not zero!'
                Write(*,*)' M=',Jm,'  M0=',Jm0
                Stop
            End If
        End If

        ! Divide workload by npes
        If (mype==0) Then
            start=1
        Else
            start=mype*(Nd0/npes)+1
        End If
        If (mype==npes-1) Then
            end1 = Nd0
        Else
            end1 = (mype+1)*(Nd0/npes)
        End If

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        Do n=start,end1
            xn0=X0(n)
            xn2=X2(n)
            idet0(1:Ne)=Iarr0(1:Ne,n)
            If ((dabs(xn0)+dabs(xn2)).GT.trd+trd) Then
                nskip=0             !### dets in skipped confs
                k=0
                Do ic=1,Nc
                    kx=Ndc(ic)
                    Call Gdet(k+1,idet1)
                    Call CompCD(idet0,idet1,icomp)
                    If (icomp.GT.1 .OR. Jdel.GT.idif) Then
                        k=k+kx
                        nskip=nskip+kx
                    Else
                        Do k1=1,kx
                            k=k+1
                            Call Gdet(k,idet1)
                            Call Rspq(idet0,idet1,is,nf,i0,i1,j0,j1)
                            If (nf.EQ.1) Then
                                If (dabs(xn0).GT.trd) Then
                                    YY1(k)=YY1(k)+is*Hint(j1,j0,kli)*xn0
                                End If
                                If (dabs(xn2).GT.trd) Then
                                    YY2(k)=YY2(k)+is*Hint(j0,j1,klf)*xn2
                                End If
                            End If
                            If (Kli.EQ.5 .AND. nf.EQ.0) Then
                                Do iq=1,Ne
                                    i1=idet1(iq)
                                    If (dabs(xn0).GT.trd) YY1(k)=YY1(k)+is*Hint(i1,i1,kli)*xn0
                                    If (dabs(xn2).GT.trd) YY2(k)=YY2(k)+is*Hint(i1,i1,klf)*xn2
                                End Do
                            End If
                        End Do
                    End If
                End Do
            Else
                negl=negl+1
            End If
        End Do
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_AllReduce(MPI_IN_PLACE, YY1, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_AllReduce(MPI_IN_PLACE, YY2, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

        if (mype == 0) then
            Open(33,file='INE.YY1',status='UNKNOWN',form='UNFORMATTED')
            Write(33) YY1(1:Nd)
            Close(33)
            print*, 'Vector YY1 written to INE.YY1'
            Open(33,file='INE.YY2',status='UNKNOWN',form='UNFORMATTED')
            Write(33) YY2(1:Nd)
            Close(33)
            print*, 'Vector YY2 written to INE.YY2'
        end if
    End Subroutine Vector

    Real(dp) Function Hint(i1,i0,kll)          !### <i1|A|i0>
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i1, i0, is, k, kll, ki, j0, l0, m0, kf, j1, l1, m1, jdel, mdel, lx
        Real(dp) :: xj0, xm0, xj1, xl0, xm1, q, xl1, h, a1, w3, w6
        ki=Nh(i0)
        j0=Jj(ki)
        l0=Ll(ki)
        m0=Jz(i0)
        kf=Nh(i1)
        j1=Jj(kf)
        l1=Ll(kf)
        m1=Jz(i1)
        xj0=j0/2.d0
        xm0=m0/2.d0
        xj1=j1/2.d0
        xl0=l0
        xm1=m1/2.d0
        q=xm1-xm0
        xl1=l1
        jdel=iabs(j1-j0)/2
        mdel=iabs(m1-m0)/2
        lx=(l1+l0+1)/2
        h=0.d0
        If (kll.EQ.1) then                         !### PNC
           If (jdel.EQ.0.AND.mdel.EQ.0) &          !### minus here accounts
                h=-Fint(5,kf,ki,-1)                !### for (-i) in H_pnc
        End If                                     !### (see AmpPNC in dtm.f)
        If (kll.LE.4) Then
          If (jdel.LE.1.AND.mdel.LE.1) then
             w3= FJ3(xj1,1.d0,xj0,-xm1,q,xm0)
             If (kll.EQ.2.OR.kll.EQ.4) then          !### E1
                k=(j1+m1+j0+1)/2+lx
                If (k.EQ.2*(k/2)) then
                   is=1
                Else
                   is=-1
                End If
                a1=is*dsqrt((j1+1)*(j0+1)*lx*1.d0)
                w6= FJ6(xl1,xj1,0.5d0,xj0,xl0,1.0d0)
                If (kll.EQ.2) then                   !### L-gauge
                   h=a1*w3*w6*Fint(3,kf,ki,+1)
                Else                                 !### V-gauge
                   h=a1*w3*w6*Fint(6,kf,ki,-1)
                End If
             End If
             If (kll.EQ.3) then                      !### Anapole Moment
                k=(j1+j0-m1+1)/2+lx
                If (k.EQ.2*(k/2)) then
                   is=1
                Else
                   is=-1
                End If
                h=is*w3*Fint(7,kf,ki,-1)
             End If
          End If
        End If

        If (kll.EQ.5) Then                         !### E2
            If( jdel.LE.2 .AND. mdel.LE.2 .AND. l0+l1.EQ.2*((l0+l1)/2) ) then
                If (j1+j0.GE.4) then
                    w3= FJ3(xj1,2.d0,xj0,-xm1,q,xm0)
                    k=(j1-m1+j0+1)/2
                    If (k.EQ.2*(k/2)) then
                        is= 1
                    Else
                        is=-1
                    End If
                    a1= is*dsqrt((j1+1)*(j0+1)*(2*l1+1)*(2*l0+1)*1.d0)
                    w3= w3*FJ3(xl1,2.d0,xl0,0.d0,0.d0,0.d0)
                    w6= FJ6(xl0,xj0,0.5d0,xj1,xl1,2.d0)
                    h = a1*w3*w6*Fint(10,kf,ki,+1)
                End If
            End If
        End If

        Hint=h
       Return
    End Function Hint

    Real(dp) Function Fint(is,nfin,nini,ic)       !### Search for radial integrals
        Implicit None           !### ic=+/-1 accounts for nfin <--> nini symmetry
        Integer :: isg, na, nb, is, nfin, nini, ic, ind, i
        Character(Len=128) :: strfmt
        isg=1
        na=nfin
        nb=nini
        If (na.GT.nb) then
            na=nini
            nb=nfin
            isg=ic
        End If
        ind=is*IPx*IPx+(na-Nso)*IPx+(nb-Nso)
        Do i=1,Nint
            If (ind.EQ.Inte(i)) goto 210
        End Do
        strfmt='(1X,"Fint: CANNOT FIND INTEGRAL ",3I4,1X,I6)'
        Write( 6,strfmt) is,nfin,nini,ind
        Write(11,strfmt) is,nfin,nini,ind
        Write(*, *) 'Try reconstructing DTM.INT'
        Stop
210     Fint=Rnt(i)*isg
        return
    End Function Fint

    Subroutine DefJz2(idet)    !### projection of angular momentum Jz for determinant idet
        Implicit None
        Integer :: m, i, k
        Integer, Allocatable, Dimension(:) :: idet
        m=0
        Do i=1,Ne
           k=idet(i)
           m=m+Jz(k)
        End Do
        Jm0=m/2.d0
        Return
    End Subroutine DefJz2

    Subroutine SetParams
        Implicit None

        ! Specify IP1 - dimension of the matrix to solve homogeneous equation
        ! Set IP1=IP1conf for same dimensionality as in conf
        IP1=15000

        ! Specify Nddir - dimension of the matrix for initial solution by SolEq1
        ! To solve homogeneous equation for the whole matrix, Nddir=IP1
        Nddir = IP1  

        ! Specify N_it4 - number of iterations in SolEq4
        N_it4 = 100

        Gj = 0.d0

    End Subroutine SetParams

    Subroutine Init_Char(Let,Alet,Blet)
        Implicit None
    
        Character(Len=1), Dimension(6) :: Let
        Character(Len=4), Dimension(13) :: Alet
        Character(Len=4), Dimension(5) :: Blet
    
        Let(1)= 's'
        Let(2)= 'p'
        Let(3)= 'd'
        Let(4)= 'f'
        Let(5)= 'g'
        Let(6)= 'h'
    
        Alet(1)= 'A_hf'
        Alet(2)= 'B_hf'
        Alet(3)= 'E1_L'
        Alet(4)= 'EDM '
        Alet(5)= 'PNC '
        Alet(6)= 'E1_V'
        Alet(7)= 'AM  '
        Alet(8)= 'MQM '
        Alet(9)= 'M1  '
        Alet(10)='E2  '
        Alet(11)='E3  '
        Alet(12)='M2  '
        Alet(13)='M3  '
    
        Blet(1)= 'Rint'
        Blet(2)= 'RPA1'
        Blet(3)= 'RPA2'
        Blet(4)= 'RPA3'
        Blet(5)= 'RPA4'

        Return
    End Subroutine Init_Char

    Subroutine recunit
        ! This Subroutine determines the record unit
        Implicit None

        Integer          :: lrec, Iflag, nbytes
        Character(Len=8) :: d1, t1, d2, t2

        t1='abcdefgh'
        d1='        '
        t2='hgfedcba'
        d2='        '
        lrec=0
        Iflag=1
200     lrec=lrec+1
        If (lrec.gt.8) Then
          Write(*,*)  'lrec > 8'
          Stop
        End If
        Open(unit=13,file='test.tmp',status='unknown',access='direct',recl=lrec)
        Write(13,rec=1,err=210) t1
        Write(13,rec=2,err=210) t2
        Read(13,rec=1,err=210) d1
        Read(13,rec=2,err=210) d2
        If (d1.ne.t1) goto 210
        If (d2.ne.t2) goto 210
        Iflag=0
210     Close(unit=13,status='delete')
        If (Iflag.ne.0) goto 200
        nbytes=8/lrec
        ipmr=4/nbytes

        Return
    End Subroutine recunit

    Subroutine Input
        Use conf_init, Only : ReadConfInp, ReadConfigurations
        Implicit None
    
        Integer :: i
        character(Len=5), Dimension(5) :: str
        Character(Len=128) :: strfmt
        data str /'H_pnc','E1(L)','H_am','E1(V)',' E2  '/

        Call recunit
        Call Init_Char(Let,Alet,Blet)
        
        Open(unit=11,status='UNKNOWN',file='FORMY.RES')
        Close(unit=11,status='DELETE')
        Open(unit=11,status='NEW',file='FORMY.RES')

        Write(*,'(A)') ' kIters= (0- invert and iterate if diverged, 1-invert only, 2-2-step iteration)'
        Read(*,*) kIters
        Write(*,'(A,I2)') ' kIters=', kIters

        Write(*,'(A)')' kl= (0-new, 1-use X1, 2-use X1,Y1,Y2 ):'
        Read (*,*) kl
        Write(*,'(A,I2)')' kl=',kl

        Write(*,'(A)')' R.H.S. of equation: H_p, E1(L), E2 (1,2,5):'
        Read (*,*) Kli
        Write(*,'(A,I2)')' Kli=',Kli

        Write(*,'(A)')' L.H.S.: H_p, E1(L), H_am, E1(V), E2 (1-5):'
        Read (*,*) Klf
        Write(*,'(A,I2)')' Klf=',Klf

        If (Klf.EQ.3) Then
          Write(*,'(A)')' Nuclear spin ='
          Read (*,*) XInuc
          Write(*,'(A,F10.6)')' XInuc=',XInuc
        End If
        Write(*,'(A)')' X0 is in file CONF0.XIJ, record number:'
        Read (*,*) N0
        Write(*,'(A,I2)')' N0=',N0
        If(N0.LE.0) Stop
        Write(*,'(A)')' X2 is in file CONF0.XIJ, record number:'
        Read (*,*) N2
        Write(*,'(A,I2)')' N2=',N2

        Read (*,*) nrange

        W0= 0.d0
        If (Kli.EQ.2 .AND. Klf.EQ.2 .OR. Kli.EQ.5) Then
            allocate(xlamb1s(nrange),xlamb2s(nrange),xlambsteps(nrange))
            Write(*,'(A)')' Give a list of ranges of wavelengths (nm) followed by step size (e.g. 535 537 0.5):'
            Write(*,'(A)')' (for static polarizability give 0 0 0)'
            Do i=1,nrange
                Read(*,*) xlamb1s(i), xlamb2s(i), xlambsteps(i)
                Write(*,'(A,F11.3,A,F11.3,A,F11.4)')' lambda=',xlamb1s(i),' nm to ',xlamb2s(i),' in steps of ',xlambsteps(i)
            End Do
        End If

        If (Kli.EQ.5 .AND. W0.EQ.0.d0) Then
          Write(*,'(/A)')' Give W00, so that H --> H + i*W00'
          Write(*,'(A)')' for solving (H - E0)|X1> = E2|X0>'
          Write(*,'(A)')' (a standard value: W00= 1.d-7)'
          Read (*,*) W00
          Write(*,'(A,1pD8.1)')' W00=',W00
        End If

        strfmt = '(1X,70("#"),/1X,"Program FormY v1.0",5X,"R.H.S.: ",A5," L.H.S.: ",A5)'
        Write( 6,strfmt) str(kli),str(klf)
        Write(11,strfmt) str(kli),str(klf)

        Call ReadConfInp
        Crit1=crt4
        IPad=3+Nlv ! set number of probe vectors to be 3 + number of levels in CONF.XIJ
        Call ReadConfigurations
        Return
    End Subroutine Input

    Subroutine Init
        Implicit None
        Integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, If, &
                    ii, nmin, jlj, i0, err_stat
        Real(dp) :: d, c1, c2, z1
        Real(dp), Dimension(IP6)  :: p, q, p1, q1 
        Real(dp), Dimension(4*IP6):: pq
        Character(Len=256) :: strfmt, err_msg
        logical :: longbasis
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        c1 = 0.01d0
        mj = 2*abs(Jm)+0.01d0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
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
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        z1 = pq(1)
        If (abs(Z-z1) > 1.d-6) Then
            Write( 6,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Write(11,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Read(*,*)
        End If

        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns  =pq(2)+c1
        Ii  =pq(3)+c1
        H   =pq(6)
        Kt  =pq(9)+c1
        Rnuc=pq(13)
        longbasis=abs(PQ(20)-0.98765d0) < 1.d-6
        Write( 6,'(4X,"Kli =",I3,7X,"Klf =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kli,Klf,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,'(4X,"Kli =",I3,7X,"Klf =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kli,Klf,Z,Jm,Nsp,Ns,Nso,Nc
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

        Write( 6,'(1X,71("="))')
        Write(11,'(1X,71("="))')

        ! Stop program if parameter Gj /= 0
        If (Gj.NE.0.d0) Then
            Write (*,*) ' Gj =',Gj
            Write (*,*) ' This code works only for Gj=0!'
            Stop
        End If

        ! Specify Nmax - the maximum dimension of vectors
        Nmax=IP4             

        ! Diagnostic for combination LHS and RHS operators
        ok = (Kli.EQ.1.AND.Klf.EQ.2).OR.(Kli.EQ.1.AND.Klf.EQ.4).OR. &
             (Kli.EQ.2.AND.Klf.LE.3).OR.(Kli.EQ.5.AND.Klf.EQ.5)
        If (.NOT.ok) Then
           Write(*,*)' Unknown combination Kli =',Kli,' Klf =',Klf
           Stop
        End If

        Return
    End subroutine Init

    Subroutine Rint      !### takes radial integrals
        Implicit None
        Integer :: i, is, ns1, nso1, err_stat
        Real(dp) :: rn1, x, z1
        Character(Len=256) :: strfmt, err_msg
        Integer, Allocatable, Dimension(:) :: l1
        Integer, Dimension(13) :: ki

        Open(13,file='DTM.INT',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file DTM.INT is absent or damaged"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If

        Read(13) ns1,nso1,z1,rn1
        Allocate(l1(ns1))
        x=iabs(ns1-Ns)+iabs(nso1-Nso)+dabs(z1-Z)+dabs(rn1-Rnuc)
        If (x.LT.1.d-6) Then
            Read(13) Nint,(l1(i),i=1,ns1)
            Allocate(Rnt(Nint),Inte(Nint))
            is=0
            Do i=1,Nsu
                is=is+iabs(Ll(i)-l1(i))
            End Do
            If (is.EQ.0) Then
                Read(13) (ki(i),i=1,13)
                Read(13) (Rnt(i),Inte(i),i=1,Nint)
                strfmt = '(/1X,"### Radial integrals from DTM.INT ("," Nint =",I6,") ###", &
                       /(4X,A4," calculated by ",A4))'
                Write(6, strfmt) Nint,(Alet(i),Blet(ki(i)),i=1,13)
                Write(11,strfmt) Nint,(Alet(i),Blet(ki(i)),i=1,13)
                Close(13)
                Deallocate(l1)
                Return
            Else
                Write (*,*) ' Rint: Wrong order of orbitals in DTM.INT'
            End If
        Else
            Write (*,*) ' Rint: Wrong basis set parameters in DTM.INT:'
            Write (*,*) ' Ns=',Ns,ns1,' Nso=',Nso,nso1
            Write (*,*) ' Z=',Z,z1
            Write (*,*) ' Rnuc=',Rnuc,rn1
        End If
        Stop
    End Subroutine Rint

    Subroutine Init0         
        !### Reads Vectors X0,X2
        Implicit None
        Integer :: nx, n, i, ij0, ij2, numj, k, l, err_stat
        Real(dp) :: xj0, xj2, t
        Character(Len=256) :: strfmt, err_msg

        ! Read number of used orbitals in X0/X2 space from file CONF0.DET
        Open(17,file='CONF0.DET',status='OLD',form='UNFORMATTED')
        Read (17) Nd2,nsu2                
        Nsu=max(nsu2,Nsu)     
        Nd0=Nd2
        If (.not. allocated(Iarr0)) Allocate(Iarr0(Ne,Nd0))
        Do i=1,Nd0
           Read(17) Iarr0(1:Ne,i)
        End Do            
        Close(17)                   

        ! Read eigenvectors X0/X2 (initial/final) from file CONF0.XIJ
        Open(unit=16,file='CONF0.XIJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF0.XIJ is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If

        Read(16) E0,Tj0,Nd0
        Allocate(X0(Nd0),X2(Nd0))
        rewind(16)
        nx=N0
        If (nx.LT.N2) nx=N2
        Do n=1,nx
            If (n.EQ.N0) Then
                Read(16) E0,Tj0,Nd0,(X0(i),i=1,Nd0)
                E0=-E0
                If (N2.EQ.N0) Then
                    E2=E0
                    Tj2=Tj0
                    Do i=1,Nd0
                       X2(i)=X0(i)
                    End Do
                End If
            Else
                If (n.EQ.N2) Then
                    Read(16) E2,Tj2,Nd0,(X2(i),i=1,Nd0)
                    E2=-E2
                Else
                    Read(16)
                End If
            End If
        End Do
        Close(16)

        Tj0=Anint(Tj0*100.0)/100.0

        Open(unit=18,file='CONF0.JJJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat == 0) Then
            ij0=2*Tj0+0.1d0            !### Rounding up momenta to
            Tj0=0.5d0*ij0              !### half-integer values
            ij2=2*Tj2+0.1d0
            Tj2=0.5d0*ij2
            
            xj0=0.d0
            xj2=0.d0
            Read(18) numj
            Do i=1,numj
                Read(18) k,n,t
                l=2
                If (n.EQ.k) l=1
                xj0=xj0+l*X0(k)*t*X0(n)
                xj2=xj2+l*X2(k)*t*X2(n)
                numj=numj+1
            End Do
            Close(18)
            xj0=0.5d0*(dsqrt(1.d0+xj0)-1.d0)
            xj2=0.5d0*(dsqrt(1.d0+xj2)-1.d0)
            strfmt='(3X,63("="),/4X,"Init0: J_0 =",F13.9," (rounded to", &
                    F4.1,"); E_0 =",F12.6,/10X," J_2 =",F13.9, &
                    " (rounded to",F4.1,"); E_2 =",F12.6)'
            Write( 6,strfmt) xj0,Tj0,E0,xj2,Tj2,E2
            Write(11,strfmt) xj0,Tj0,E0,xj2,Tj2,E2
        End If

        Return
    End Subroutine Init0

end program formy