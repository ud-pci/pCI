Program ine         
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Solution of the inhomogeneous linear equation
    !  (Ei-H)|X1> = A|X0>
    ! in CI space (Sternheimer method).
    ! A=H_pnc (Ei=E0) or A=E1(L) (Ei=E2).
    ! After that ME <X2|B|X1> is calculated,
    ! where B=H_pnc, E1(L,V) or H_am
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !   Keys Kli and Klf define RHS and LHS operators A and B:
    !   Kli=1 -  A = H_pnc            Klf=1 -  B = H_pnc
    !   Kli=2 -  A = E1(L-gauge)      Klf=2 -  B = E1(L-gauge)
    !                                 Klf=3 -  B = H_am
    !                                 Klf=4 -  B = E1(V-gauge)
    !   - - - - - - - - - - - - - - - - - - - - - - - - -
    !   Kli=5 - A = E2                Klf=5 -  B = E2
    !   - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   INPUT/OUTPUT files:
    !   CONF.INP    CI input for X1 space
    !   CONF.DAT    radial functions
    !   CONF.HIJ    Hamiltonian matrix in X1 space
    !   CONF.JJJ    J**2 in the CI space of vector X1
    !   CONF.XIJ    is used to calculate <n|X1>
    !   DTM.INT     Radial integrals
    !   INE.XIJ     solution and R.H.S of the inhom. eq.
    !   INE_J.XIJ   decomposition of the solution of the inhom. eq.
    !   CONF0.XIJ   eigenvectors X0 (initial) and X2 (final)
    !   CONF0.DET   determinants in X0/X2 space
    !   CONF0.JJJ   J**2 in X0/X2 space
    !   - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main variables:
    !   Ns      - number of orbitals (different)
    !   Ne      - number of valence electrons
    !   Nec     - number of core electrons
    !   Nso     - number of core orbitals
    !   - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main variables in X1 space:
    !   Nsp      - total number of shells (>=Ns, can be equal).
    !   Qnl(i)   - atomic configurations (i=1,Nsp)
    !   Jm       - projection of total momentum J
    !   Nc       - number of configurations
    !   Nd       - number of dets
    !   IPmr     - equivalent of 4 Bytes for DIRECT files
    !   - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main variables in X0 space:
    !   Jm0      - projection of total momentum J
    !   Nc0      - number of configurations
    !   Nd0      - number of dets
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||
    !         General convention: GLOBAL VARIABLES     |||
    !              ARE CAPITALISED                     |||
    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||
    Use params, ipmr1 => IPmr, IP1conf => IP1
    Use determinants, Only : Dinit, Jterm
    Use str_fmt, Only : startTimer, stopTimer

    Implicit None
    
    Integer, Parameter :: IP2 = 50000
    Integer, Parameter :: IPad = 8

    Integer :: i, n, k, l, Nd2, Nddir, nsu2, icyc, nlamb, kIters, N_it4
    Integer  :: Kli, Klf, Khe, Ndir, Int_err, Ntr, Nint, Nmax, Nd0, ipmr, Kdiag, Nlft, IP1
    Integer  :: IPlv, IP4, Kt, Nlev, N0, N2, nrange
    Integer(kind=int64) :: NumH, NumJ
    
    Real(dp) :: Jm0, E0, E2, Tj0, Tj2, xlamb, xlamb1, xlamb2, xlambstep, XInuc, Crit1, W0, W00, Q, Elft, Hmin, dlamb
    Real(dp), Allocatable, Dimension(:) :: xlamb1s, xlamb2s, xlambsteps, xlamblist
    Integer, Allocatable, Dimension(:) :: Int
    Real(dp), Allocatable, Dimension(:) :: Z1, X0, X1, X2, YY1, YY2, Rnt, Ev, Diag, E1
    Real(dp), Allocatable, Dimension(:,:) :: X1J, Y2J
    Real(dp), Dimension(2) :: s, ss, s0, s1, s2
    logical :: ok

    Character(Len=1), Dimension(9)          :: Let
    Character(Len=4), Dimension(13)         :: Alet
    Character(Len=4), Dimension(5)          :: Blet

    Integer(Kind=int64) :: start_time
    Character(Len=16) :: timeStr

    Type Matrix
        Integer,  Allocatable, Dimension(:) :: n, k
        Real(dp), Allocatable, Dimension(:) ::  t
    End Type Matrix
    Type(Matrix) :: Hamil, Jsq

    Call startTimer(start_time)

    Call SetParams                    ! Set job and array parameters
    Call Input                        ! Read list of configurations from CONF.INP
    Call Init                         ! Read basis set information from CONF.DAT
    Call Rint                         ! Read radial integrals from CONF.INT
    Call Dinit                        ! Form list of determinants
    Call Jterm                        ! Print table with numbers of levels with given J
    Call ReadHIJ                      ! Read Hamiltonian matrix in X1 space from file CONF.HIJ
    Call ReadJJJ                      ! Read J^2 matrix in X1 space from file CONF.JJJ
    Call Init0                        !### Evaluation of the RHS of
    Call Vector(kl)                   !###  the equation and vectors Yi

    If (W0 == 0.d0 .or. Kli == 5) Then
        icyc=1
    Else
        icyc=2
    End If
    
    Do k=1,nrange
        xlamb1 = xlamb1s(k)
        xlamb2 = xlamb2s(k)
        xlambstep = xlambsteps(k)
        If (k > 1) Then
            If (xlamb1 == xlamb2s(k-1)) xlamb1 = xlamb1+xlambstep
            If (xlamb2 < xlamb1 + xlambstep) xlamb2 = xlamb1
        End If
        If (xlamb2==xlamb1) xlambstep=1

        dlamb = (xlamb2-xlamb1)/xlambstep
        dlamb = Anint(dlamb*100.0)/100.0
        nlamb = dlamb+1
        print*, xlamb2, xlamb1, xlambstep, dlamb
        print*,'nlamb=',nlamb
        xlamb=xlamb1
        Do n=1,nlamb
            If (W0 /= 0.d0) W0 = 1.d+7/(xlamb*219474.63d0)
            print*, 'xlamb=',xlamb
            Do i=1,icyc
                Ndir=Nddir
                If (Kli.EQ.5) Then
                    Ndir= Nd    ! SolEq4 is not adopted yet for E2 polariz.
                    If (Ndir > IP2) Then
                        Print*, 'WARNING: Dimensionality of matrix exceeds', IP2
                    End If
                End If
                If (i.EQ.2) W0= -W0
                If (dabs(xlamb).GT.1.d-8) Then
                    If (i.EQ.2) xlamb= -xlamb
                    Write( *,'(/3X,34("-")/3X,"Calculation for lambda=",F11.3,/3X,34("-")/)') xlamb
                    Write(11,'(/3X,34("-")/3X,"Calculation for lambda=",F11.3,/3X,34("-")/)') xlamb
                Else
                    Write( *,'(/3X,22("-")/3X,"Calculation for W0 = 0",/3X,22("-")/)')
                    Write(11,'(/3X,22("-")/3X,"Calculation for W0 = 0",/3X,22("-")/)')
                End If
                Select Case(kIters)
                    Case(0)
                        print*, 'Ndir=',Ndir
                        Call SolEq1(kl) ! Direct solution
                        If (Ndir.LT.Nd) Then
                            Call SolEq4(ok)                 !### Iterative solution
                            If (Nd.LE.IP1 .AND. .NOT.ok) Then
                              Ndir= Nd
                              Write(*,*)
                              Call SolEq1(kl)
                              ok=.TRUE.
                            End If
                        End If
                    Case(1)
                        Ndir=Nd
                        IP1=Nd
                        Call SolEq1(kl) ! Direct solution
                    Case(2)
                        N_it = 2
                        Do l=1,N_it4
                            print*, '   2-step ITERATION #', l
                            If (l > 1) kl = 2
                            Call SolEq1(kl) ! Direct solution
                            If (Ndir.LT.Nd) Then
                                Call SolEq4(ok) ! Iterative solution
                            End If
                            If (ok) Exit
                        End Do
                End Select
                If (Kli.LE.2) Call  Prj ('  X1  ',Tj0,X1,X1J)     !### Projects X1 on J subspaces
                If (Kli.EQ.5) Call PrjE2('  X1  ',Tj0,X1,X1J)
                Call RdcX1J                                       !### Transforms and saves X1J
                If (Kli.LE.2) Call  Prj ('  Y2  ',Tj2,YY2,Y2J)    !### Projects Y2 on J subspaces
                If (Kli.EQ.5) Call PrjE2('  Y2  ',Tj2,YY2,Y2J)
                Call Prin                                  !### Output of the results
                If (Kli.EQ.2.AND.Klf.EQ.2) Then
                  Call RdcE1(i)                            !### Evaluation of E1 polarizability
                  If (W0.EQ.0.d0 .OR. i.EQ.2) Call C_3     !### C_3 coefficient for X2 state
                End If
                If (Kli.EQ.5) Call RdcE2(i)                !### Evaluation of E2 polarizabilty
                If (Kli.EQ.1.OR.Klf.EQ.1) Call RdcPNC      !### Final numbers for Q_w
                If (Klf.EQ.3) Call RdcAM                   !### Final numbers for AM
                If (Int_err.NE.0) Then
                  Write( *,'(4X,">>>> NOTE:",I7," radial integrals were absent")') Int_err
                  Write(11,'(4X,">>>> NOTE:",I7," radial integrals were absent")') Int_err
                End If
                If (.NOT.ok) Then
                  Write( *,'(4X,"Convergence was not reached.")')
                  Write(11,'(4X,"Convergence was not reached.")')
                End If
            End Do
            xlamb=abs(xlamb)+xlambstep
            print*,'xlamb_next=',xlamb
        End Do
    End Do

    Close(unit=11) ! Close INE.RES
    Close(unit=99) ! Close INEFINAL.RES
    Call stopTimer(start_time, timeStr)
    Write(*,'(2X,A)'), 'TIMING >>> Total computation time of ine was '// trim(timeStr)
    
Contains

    Subroutine SetParams
        Implicit None

        ! Choose Khe - variant of solution of homogeneous equation
        ! Khe=0 - old solution of homogeneous eq-n
        ! Khe=1 - new solution of homogeneous eq-n
        Khe= 1   

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
        
        Open(unit=11,status='UNKNOWN',file='INE.RES')
        Close(unit=11,status='DELETE')
        Open(unit=11,status='NEW',file='INE.RES')

        Open(unit=99,status='UNKNOWN',file='INEFINAL.RES')
        Close(unit=99,status='DELETE')
        Open(unit=99,status='NEW',file='INEFINAL.RES')

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
            
            xlamb1 = xlamb1s(1)
            xlamb2 = xlamb2s(1)
            xlambstep = xlambsteps(1)
            If (dabs(xlamb1).GT.1.d-8) Then
              W0 = 1.d+7/(xlamb1*219474.63d0)
            Else
              Write(*,'(A)')' W0 = 0'
            End If
        End If

        If (Kli.EQ.5 .AND. W0.EQ.0.d0) Then
          Write(*,'(/A)')' Give W00, so that H --> H + i*W00'
          Write(*,'(A)')' for solving (H - E0)|X1> = E2|X0>'
          Write(*,'(A)')' (a standard value: W00= 1.d-7)'
          Read (*,*) W00
          Write(*,'(A,1pD8.1)')' W00=',W00
        End If

        strfmt = '(1X,70("#"),/1X,"Program InhomEq. v1.22",5X,"R.H.S.: ",A5," L.H.S.: ",A5)'
        Write( 6,strfmt) str(kli),str(klf)
        Write(11,strfmt) str(kli),str(klf)

        Call ReadConfInp
        Crit1=crt4
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
            Allocate(Rnt(Nint),Int(Nint))
            is=0
            Do i=1,Nsu
                is=is+iabs(Ll(i)-l1(i))
            End Do
            If (is.EQ.0) Then
                Read(13) (ki(i),i=1,13)
                Read(13) (Rnt(i),Int(i),i=1,Nint)
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

    Subroutine Vector(kl)   
        ! Read RHS vector |Y1>=A|X0> and LHS vector <Y2|=<X2|B,
        Use conf_variables, Only : iconf1, iconf2
        Use determinants, Only : Gdet, Rspq, CompC
        Implicit None
        Integer :: kl, i, ic, icomp, is, i0, i1, j0, j1, negl, idif, n, k, kx, err_stat
        Integer :: mdel, nskip, k1, nf, iq
        Integer, Allocatable, Dimension(:) :: idet0, idet1
        Real(dp) :: trd, xn0, xn2, s1, s2
        Character(Len=256) :: strfmt, err_msg

        trd=1.d-10                  !### used to drop small contributions
        Open (unit=17,file='CONF0.DET',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF0.DET is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If
        Read(17) Nd0
        
        If (.not. Allocated(YY1)) Allocate(YY1(Nd))
        If (.not. Allocated(YY2)) Allocate(YY2(Nd))
        If (.not. Allocated(idet0)) Allocate(idet0(Ne))
        If (.not. Allocated(idet1)) Allocate(idet1(Ne))
        If (.not. Allocated(iconf1)) Allocate(iconf1(Ne))
        If (.not. Allocated(iconf2)) Allocate(iconf2(Ne))

        strfmt = '(3X,63("="),/4X,"Vector: forming vectors Y1 and Y2"," Threshold: ",E9.1)'
        If (kl.NE.2) Write( 6,strfmt) trd
        YY1(:)=0.d0
        YY2(:)=0.d0
        negl=0              !### number of vector components below threshold
        idif= 1
        If (Kli.EQ.5) idif= 2
        Do n=1,Nd0
            xn0=X0(n)
            xn2=X2(n)
            Read(17) (idet0(i),i=1,Ne)
            If (n.EQ.1) Then
                Call DefJz2(idet0)
                Q=Jm-Jm0
                mdel=dabs(Q)+1.d-5
                If ((Kli-1)*(Klf-1).EQ.0.AND.mdel.NE.0) Then
                    Write(*,*) 'Vector: |M-M0| is not zero!'
                    Write(*,*)' M=',Jm,'  M0=',Jm0
                    Stop
                End If
                If (kl.EQ.2) return
            End If
            If ((dabs(xn0)+dabs(xn2)).GT.trd+trd) Then
                nskip=0             !### dets in skipped confs
                k=0
                Do ic=1,Nc
                    kx=Ndc(ic)
                    Call Gdet(k+1,idet1)
                    Call CompC(idet0,idet1,icomp)
                    If (icomp.GT.1 .OR. Jdel.GT.idIf) Then
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
            If (1000*(n/1000).EQ.n) Then
                 strfmt='(1X,I6," from ",I6," (",I6," neglected, ",I6," skipped)")'
                 Write(*,strfmt) n,Nd0,negl,nskip
            End If
        End Do

        If (Kli.EQ.5 .AND. Nd0.EQ.Nd) Then
            s1= 0.d0
            s2= 0.d0
            Do i=1,Nd0
                s1= s1 + X2(i)*YY1(i)
                s2= s2 + X0(i)*YY2(i)
            End Do
            strfmt='(/4x,"<X2|Q_0|X0>=",F10.5,/4x"<X0|Q_0|X2>=",F10.5)'
            Write (*,strfmt) s1,s2
        End If
        Close(17)

        Return
    End Subroutine Vector

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
            If (ind.EQ.Int(i)) goto 210
        End Do
        strfmt='(1X,"Fint: CANNOT FIND INTEGRAL ",3I4,I6)'
        Write( 6,strfmt) is,nfin,nini,ind
        Write(11,strfmt) is,nfin,nini,ind
        Fint=0.d0
        Int_err=Int_err+1
        return
210     Fint=Rnt(i)*isg
        return
    End Function Fint

    Subroutine ReadHIJ
        Implicit None
        Integer :: err_stat
        Integer(kind=int64) :: i8
        Character(Len=256) :: err_msg

        Open(unit=15,file='CONF.HIJ',status='UNKNOWN',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
           Write(*,*) ' file CONF.HIJ is not found'
           Stop
        Else
            Write(*,*) ' reading CONF.HIJ...'
        End If

        Read(15) NumH
        print*,'NumH=',NumH
        If (.not. Allocated(Hamil%n)) Allocate(Hamil%n(NumH))
        If (.not. Allocated(Hamil%k)) Allocate(Hamil%k(NumH))
        If (.not. Allocated(Hamil%t)) Allocate(Hamil%t(NumH))
        Do i8=1,NumH
            Read(15), Hamil%k(i8), Hamil%n(i8), Hamil%t(i8)
        End Do
        Return
    End Subroutine

    Subroutine SolEq1(kl)
        Use solvers
        Use str_fmt, Only : FormattedTime
        Implicit None
        External ZSPSV
        Integer :: i, ic, id, info, j, k, k1, k2, kl, nd1, err_stat
        Integer(Kind=int64) :: i8, start1, end1, start2, end2, clock_rate
        Real :: ttime, ttime2
        Real(dp) :: t, e0n, e2n, tj0n, tj2n, err
        Real(dp), Allocatable, Dimension(:) :: scales
        Complex(dp), Allocatable, Dimension(:) ::Zc, XX1
        Integer, Allocatable, Dimension(:) :: ipiv, Mps
        Real(dp), Allocatable, Dimension(:,:) :: Az
        Character(Len=256) :: strfmt, err_msg, timeStr

        Call system_clock(count_rate=clock_rate)
        If (.not. Allocated(scales)) Allocate(scales(IP1))
        If (.not. Allocated(Mps)) Allocate(Mps(IP1))
        If (.not. Allocated(ipiv)) Allocate(ipiv(IP1))
        If (kl > 0) Then
            If (.not. Allocated(X1)) Allocate(X1(Nd))
            If (.not. Allocated(YY1)) Allocate(YY1(Nd))
            If (.not. Allocated(YY2)) Allocate(YY2(Nd))
        End If
        
        elft= E0+W0                          !### equation has the form:
        If (Kli.EQ.2.AND.Klf.NE.2) elft=E2   !### (elft-H)X1=Y1
        If (kl.GE.1) then
            Open (unit=16,file='INE.XIJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
            If (err_stat /= 0) Then
                Continue
            Else
                Read(16,iostat=err_stat,iomsg=err_msg) e0n,tj0n,Ntr,(X1(i),i=1,Ntr)
                If (err_stat /= 0) Then
                    Continue
                Else
                    err=dabs(e0n+E0)+dabs(tj0n-Tj0)
                    If (err.GT.1.d-1) then
                        strfmt = '(1X," Error in file INE.XIJ, vector X1:", &
                                /" E0=",2F10.6,"; J0=",2F10.6)'
                        Write( *,strfmt) -e0n,E0,tj0n,Tj0
                        Write(11,strfmt) -e0n,E0,tj0n,Tj0
                        Stop
                    End If
                    If (kl.EQ.2) then
                        Read(16,iostat=err_stat,iomsg=err_msg) e0n,tj0n,Nd1,(YY1(i),i=1,Nd1)
                        If (err_stat /= 0) Then
                            Write(*,*) ' Error in INE.XIJ'
                            Stop
                        End If
                        Read(16,iostat=err_stat,iomsg=err_msg) e2n,tj2n,Nd1,(YY2(i),i=1,Nd1)
                        If (err_stat /= 0) Then
                            Write(*,*) ' Error in INE.XIJ'
                            Stop
                        End If
                        err= err + dabs(e2n+E2)+dabs(tj2n-Tj2)+iabs(Nd-Nd1)
                        If (err.GT.1.d-1) then
                            strfmt = '(1X," Error in file INE.XIJ, vectors Y1,Y2:", &
                                    /" E0=",2F10.6,"; J0=",2F10.6,/" E2=",2F10.6,"; J2=",2F10.6)'
                            Write( *,strfmt) -e0n,E0,tj0n,Tj0,-e2n,E2,tj2n,Tj2
                            Write(11,strfmt) -e0n,E0,tj0n,Tj0,-e2n,E2,tj2n,Tj2
                            Stop
                        End If
                    Else
                        Write(16) -E0,Tj0,Nd,(YY1(i),i=1,Nd)
                        Write(16) -E2,Tj2,Nd,(YY2(i),i=1,Nd)
                    End If
                    Close(16)
                    Write( *,*)'  SolEq1: solution is taken from INE.XIJ'
                    Write(11,*)'  SolEq1: solution is taken from INE.XIJ'
                    Return
                End If
            End If
        End If
        print*, 'Nd=',Nd, 'IP1=',Ndir
        If (Nd.GT.Ndir) then
            Ntr=0                          !### Ntr = dimension of the
            Do ic=1,Nc                     !### matrix equation to solve here
               id=Ndc(ic)
               If ((Ntr+id).LE.Ndir) Ntr=Ntr+id
            End Do
        Else
            Ntr=Nd
        End If
        strfmt = '(4X,"SolEq1: matrix (",I5," X ",I5,") ( Nd =",I6,")")'
        Write(6,strfmt) Ntr,Ntr,Nd

        If (.not. Allocated(Z1)) Allocate(Z1(IP1*IP1))
        If (.not. Allocated(Az)) Allocate(Az(IP1,IP1))
        Z1=0.d0
        Az=0.d0
        Call system_clock(start1)
        Hmin=0_dp
        Do i8=1,NumH         !### construction of Z matrix                      
            i=Hamil%n(i8)    !### (note that it is stored as one-dimensional array)
            j=Hamil%k(i8)
            t=Hamil%t(i8)
            if (t < Hmin) Hmin = t
            If (i.GT.Ntr) Exit                  
            k1=i+Ntr*(j-1)
            k2=j+Ntr*(i-1)
            If (i.EQ.j) then
                Z1(k1) = elft - t
                Az(i,i)= elft - t
            Else
                Z1(k1)=-t
                Z1(k2)=-t
                Az(i,j)= -t
                Az(j,i)= -t
            End If
        End Do
        Close(unit=15)
        Call system_clock(end1)
        ttime=Real((end1-start1)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        Write(*,'(2X,A)'), 'SolEq1: Z matrix calculated in '// trim(timeStr)// '.'
        
        strfmt = '(3X," NumH =",I10,";  Hmin =",F12.6/)'
        Write(*,strfmt) NumH,Hmin

        If (.not. Allocated(Zc)) Allocate(Zc(Ntr*Ntr))
        If (.not. Allocated(XX1)) Allocate(XX1(Ntr))
        If (.not. Allocated(X1)) Allocate(X1(Nd))

        Call system_clock(start1)
        If (Khe.EQ.0) Then
            call Decomp(Z1,Ntr,scales,Mps)          !### decomposition
            Write(*,*)' decomposition finished'            !### of the matrix
            call Flsolv(Ntr,Z1,YY1,X1,Mps)          !### solution of the
            Write(*,*)' flsolv finished'                   !### inhom-s equation
        Else
            k = 1
            Do i=1,Ntr
                Do j=1,i
                    If (i.EQ.j) Then
                        Zc(k)= Az(i,i)
                        If (Kli.EQ.5 .AND. W0.EQ.0.d0) Zc(k)= Az(i,i)+(0,1)*W00
                    Else
                        Zc(k)= Az(i,j)
                    End If
                    k = k+1
                End Do
                XX1(i)= YY1(i)
            End Do
            Call system_clock(start2)
            Call Zspsv('U',Ntr,1,Zc,ipiv,XX1,Ntr,info)
            !Call Zspsv_F95(Zc,XX1,'U',ipiv,info)
            Call system_clock(end2)
            ttime2=Real((end2-start2)/clock_rate)
            Call FormattedTime(ttime2, timeStr)
            Write(*,'(2X,A)'), 'SolEq1: Zspsv in '// trim(timeStr)// '.'
            If (info.NE.0) then
                Write(*,*) 'info =',info
                Stop
            End If
            X1= 0.d0  ! Vector
            X1(1:Ntr)= Real(XX1(1:Ntr))
        End If
        Call system_clock(end1)
        ttime=Real((end1-start1)/clock_rate)
        Call FormattedTime(ttime, timeStr)
        Write(*,'(2X,A)'), 'SolEq1: decomp in '// trim(timeStr)// '.'

        If (Kli.EQ.5 .AND. dabs(Tj0).GT.0.51d0) Call Ort  !# Orthogonalization of X1 to X0 (If J0 /= 0)
        
        Open(unit=16,file='INE.XIJ',status='UNKNOWN',form='UNFORMATTED')
        Write(16) -E0,Tj0,Nd,(X1(i),i=1,Nd)
        Write(16) -E0,Tj0,Nd,(YY1(i),i=1,Nd)
        Write(16) -E2,Tj2,Nd,(YY2(i),i=1,Nd)
        Close(16)
        Return
    End Subroutine SolEq1

    Subroutine SolEq4(ok)   !### Iterative solution of the inhomogeneous eqiation
        Use solvers
        Implicit None
        Integer :: num, i, itr, k, l, n
        Real(dp) :: ynorm, crit, dnorm, s2, s3, s12, s13, dx, di, x2, x3, x, y, dx1
        Character(Len=9) :: str
        Logical :: ok
        Integer, Allocatable, Dimension(:) :: Mps
        Real(dp), Allocatable, Dimension(:) :: scales, Vl, Vr
        Real(dp), Allocatable, Dimension(:) :: XX3, X12, X13
        Character(Len=256) :: strfmt

        If (.not. Allocated(scales)) Allocate(scales(IPad))
        If (.not. Allocated(Mps)) Allocate(Mps(IPad))
        If (.not. Allocated(Vl)) Allocate(Vl(IPad))
        If (.not. Allocated(Vr)) Allocate(Vr(IPad))
        If (.not. Allocated(X1J)) Allocate(X1J(Nd,IPad))
        If (.not. Allocated(Y2J)) Allocate(Y2J(Nd,5))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(Ev)) Allocate(Ev(IPad))
        If (.not. Allocated(XX3)) Allocate(XX3(Nd))
        If (.not. Allocated(X12)) Allocate(X12(Nd))
        If (.not. Allocated(X13)) Allocate(X13(Nd))
        If (.not. Allocated(Z1)) Allocate(Z1(IP1*IP1))
        If (.not. Allocated(Vr)) Allocate(Vr(IPad))

        num=4
        call FormV(num)
        num=min(Nlft,IPad)
        strfmt = '(3X,"Number of vectors: ",I2)'
        write(*,strfmt) num
        ynorm=0.d0
        do i=1,Nd
            ynorm=ynorm+YY1(i)**2
            X1J(i,1)=X1(i)
            XX3(i)=X1(i)
        end do
        ynorm=dsqrt(ynorm)
        crit=Crit1*ynorm
        strfmt = '(3X,63("="),/4X,"SolEq4: ||Y1|| = ",E10.3," Crit = ",E10.3,/3X,63("="))'
        write( 6,strfmt) ynorm,crit
        write(11,strfmt) ynorm,crit
        Kdiag=1                        !### flags Mxmpy to fill Diag(i)
        ok=.FALSE.
        Do itr=1,N_it
            call Mxmpy(1,1,X1,X1)
            dnorm=0.d0
            s2=0.d0
            s3=0.d0
            s12=0.d0
            s13=0.d0
            do i=1,Nd
                dx=YY1(i)-Y2J(i,1)         !### residual vector
                dnorm=dnorm+dx*dx
                di=Diag(i)
                if (dabs(di).GT.1.d-6) then
                    dx1=dx/di
                else
                    dx1=0.d0
                end if
                X1J(i,2)=dx                !### first probe vector
                X1J(i,3)=dx1               !### second probe vector
                x2=X1J(i,1)
                x3=XX3(i)
                s2=s2+x2*x2
                s3=s3+x3*x3
                s12=s12+x2*dx
                s13=s13+x3*dx1
            end do
            s12=s12/s2
            s13=s13/s3
            s2=0.d0
            s3=0.d0
            do i=1,Nd                    !### orthogonalization of the
                x=X1J(i,2)-s12*X1J(i,1)    !#### probe vectors to those
                y=X1J(i,3)-s13*XX3(i)      !##### from previous iteration
                s2=s2+x*x
                s3=s3+y*y
                X1J(i,1)=x
                X1J(i,2)=x
                X1J(i,3)=y
                XX3(i)=y
            end do
            dnorm=dsqrt(dnorm)
            s2=dsqrt(s2)
            s3=dsqrt(s3)
            strfmt = '(3X,73("-"),/4X,"iter=",I3, &
                  /4X,"Norm of probe vectors:",2E12.3, &
                  /4X,"Norm of the residual vector:",E12.3)'
            write( *,strfmt) itr,s2,s3,dnorm
            write(11,strfmt) itr,s2,s3,dnorm
            if (dnorm.LE.crit) then
                str='converged'
                ok=.TRUE.
            end if
            X12=X1J(1:Nd,2)
            X13=X1J(1:Nd,3)
            call Mxmpy(2,2,X12,X13)
            ! Evaluation of the coefficients from minimum of the residue:
            do i=1,IPad
                Vr(i)=0.d0
                do k=1,IPad
                    Z1(IPad*(i-1)+k)=0.d0
                end do
            end do
            Vl=0.d0
            Mps=0
            do i=1,Nd
                y=YY1(i)
                do k=1,num
                    if (k.LE.3) then
                        Vl(k)=Y2J(i,k)
                    else
                        Vl(k)=X1J(i,k)*Ev(k)
                    end if
                    Vr(k)=Vr(k)+Vl(k)*y
                end do
                do k=1,num
                    do l=1,num
                        n=(k-1)*num+l
                        Z1(n)=Z1(n)+Vl(k)*Vl(l)
                    end do
                end do
            end do
            call Decomp(Z1,num,scales,Mps)
            call Flsolv(num,Z1,Vr,Vl,Mps)
            !     New vector X1:
            do i=1,Nd
                X1(i)=Vl(1)*X1(i)
                do k=2,num
                    X1(i)=X1(i)+Vl(k)*X1J(i,k)
                end do
            end do
            strfmt = '(4X,"mixing coefficients",/(4x,8F9.4))'
            write( 6,strfmt) (Vl(k),k=1,num)
            write(11,strfmt) (Vl(k),k=1,num)
    
            open (unit=16,file='INE.XIJ',status='UNKNOWN',form='UNFORMATTED')
            write(16) -E0,Tj0,Nd,(X1(i),i=1,Nd)
            write(16) -E0,Tj0,Nd,(YY1(i),i=1,Nd)
            write(16) -E2,Tj2,Nd,(YY2(i),i=1,Nd)
            close(16)
            If (ok) Then
                Exit
            Else
                str=' DIVERGED'
                Cycle
            End If
        End Do
        strfmt = '(4X,"Iteration process ",A9,". Vectors of length ", &
             I7," are saved",/3X,63("="))'
 200    write( 6,strfmt) str,Nd
        write(11,strfmt) str,Nd
        Return
    End Subroutine SolEq4

    Subroutine ReadJJJ
        Implicit None
        Integer :: err_stat
        Integer(kind=int64) :: j8
        Character(Len=256) :: err_msg

        open(unit=18,file='CONF.JJJ',status='OLD',form='UNFORMATTED',access='stream',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
           Write(*,*) ' file CONF.JJJ is not found'
           Stop
        Else
            Write(*,*) ' reading CONF.JJJ...'
        End If
        read(18) NumJ
        If (.not. Allocated(Jsq%n)) Allocate(Jsq%n(NumJ))
        If (.not. Allocated(Jsq%k)) Allocate(Jsq%k(NumJ))
        If (.not. Allocated(Jsq%t)) Allocate(Jsq%t(NumJ))
        Do j8=1,NumJ
            read(18) Jsq%k(j8),Jsq%n(j8),Jsq%t(j8)
        End Do
        print*, 'NumJ=',NumJ
        close(18)
    End Subroutine ReadJJJ

    Subroutine Prj(str,tj0,X1,X1J)
        Implicit None
        Integer :: i, j, j2, jj, n, k
        Real(dp) :: trsd, tj, tj0, tjm, tjp, tjj, t, sm, sj, x, y, &
                    ap, aj, sum, ds, smin, sp, sn
        Real(dp), Allocatable, Dimension(:) :: X1
        Real(dp), Allocatable, Dimension(:,:) :: X1J
        Character(Len=6) :: str
        Character(Len=256) :: strfmt

        If (.not. Allocated(X1J)) Allocate(X1J(Nd,IPad))

        trsd=1.d-3            !### Threshold used to eliminate small terms
        jj=2*Tj0+0.1d0
        j2=jj*jj

        X1J(1:Nd,1:3)=0.d0

        strfmt = '(3X,"Prj: Decomposition of the vector ",A6, &
             /7X,"  Norma ",8X,"J_0-1",8X,"J_0",8X,"J_0+1",9X,"Sum")'
        write( 6,strfmt) str
        write(11,strfmt) str
        
        Do j=1,NumJ
            n=Jsq%n(j)
            k=Jsq%k(j)
            t=Jsq%t(j)
            X1J(n,1)=X1J(n,1)+t*X1(k)           !### = J^2*X1
            if (n.NE.k) X1J(k,1)=X1J(k,1)+t*X1(n)
        End Do
        If (jj.GT.1) Then !### no vector X_(j-1)
            Do j=1,NumJ
               n=Jsq%n(j)
               k=Jsq%k(j)
               t=Jsq%t(j)
               X1J(n,2)=X1J(n,2)+t*X1J(k,1)        !### = J^4*X1
               If (n.NE.k) X1J(k,2)=X1J(k,2)+t*X1J(n,1)
            End Do
        End If
        sm=0.d0
        sj=0.d0
        sp=0.d0
        sn=0.d0
        do i=1,Nd
           x=X1(i)
           y=X1J(i,1)
           z=X1J(i,2)
           if (jj.GT.1) then
             ap= (z-2*j2*y+j2*(j2-4)*x)/(32*(jj+2)*(jj+1))
             aj=-(z-2*(j2+2*jj+4)*y+jj*(j2-4)*(jj+4)*x)/(16*jj*(jj+2))
             am= (z-2*(jj+2)**2*y+jj*(jj+2)**2*(jj+4)*x)/(32*jj*(jj+1))
           else
             ap= (y-jj*(jj+2)*x)/(4*(jj+2))
             aj=-(y-(jj+2)*(jj+4)*x)/(4*(jj+2))
             am= 0.d0
           end if
           X1J(i,1)=am
           X1J(i,2)=aj
           X1J(i,3)=ap
           sm=sm+am*am
           sj=sj+aj*aj
           sp=sp+ap*ap
           sn=sn+x*x
        end do
        sum=sm+sj+sp
        tj =0.5d0*jj
        tjm=tj-1.d0
        tjp=tj+1.d0
        tjj=(tjm*(tjm+1)*sm+tj*(tj+1)*sj+tjp*(tjp+1)*sp)/sum
        tjj=0.5d0*(dsqrt(1.d0+4*tjj)-1.d0)
        strfmt = '(2X,5E13.5,/4X,"J_av=",F11.8)'
        write( 6,strfmt) sn,sm,sj,sp,sum,tjj
        write(11,strfmt) sn,sm,sj,sp,sum,tjj
        ds=dabs(sum/sn-1.d0)
        if (ds.GT.trsd) then
           write (*,*) ' Prj warning: normalization has changed by',ds
           read(*,*)
        end if
        smin=dmin1(sm/sum,sj/sum,sp/sum)
        !if (smin.LT.trsd) then                   !### Small components of
        !   do n=1,Nd                             !#### the J-decomposition
        !     if (sm/sum.LT.trsd) X1J(n,1)=0.d0   !#### are put to zero
        !     if (sj/sum.LT.trsd) X1J(n,2)=0.d0
        !     if (sp/sum.LT.trsd) X1J(n,3)=0.d0
        !   end do
        !end if
        Return
    End Subroutine Prj

    Subroutine PrjE2(str,Tj0,X1,X1J)
        ! We decompose X1 as X1 = X1_J-2 + X1_J-1 + X1_J + X1_J+1 + X1_J+2
        Implicit None
        Integer :: j, j0, n, k, i
        Real(dp) :: t, sum, sum1, sum2, sum3, sum4, sum5, &
                    XJ3xXJ4, XJ3xXJ5, XJ4xXJ5, Tj0, f, sum_tot, ds
        Real(dp), Allocatable, Dimension(:) :: X1, TJ1, TJ2, TJ3, TJ4, TJ5
        Real(dp), Allocatable, Dimension(:,:) :: X1J
        Character(Len=6) :: str
        Character(Len=256) :: strfmt
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        If (.not. Allocated(TJ1)) Allocate(TJ1(Nd))
        If (.not. Allocated(TJ2)) Allocate(TJ2(Nd))
        If (.not. Allocated(TJ3)) Allocate(TJ3(Nd))
        If (.not. Allocated(TJ4)) Allocate(TJ4(Nd))
        If (.not. Allocated(TJ5)) Allocate(TJ5(Nd))
        If (.not. Allocated(X1J)) Allocate(X1J(Nd,IPad))
        j0= 2*Tj0+0.1d0

        TJ1=0.d0
        TJ2=0.d0
        TJ3=0.d0
        TJ4=0.d0
        TJ5=0.d0
        X1J(1:Nd,1:5)=0.d0

        strfmt = '(3X,"PrjE2: Decomposition of the vector ",A6, &
                /5X," Norma ",7X,"J0-2",9X,"J0-1",10X,"J0",10X, &
                              "J0+1",9X,"J0+2",9X"Sum")'
        write(*,*)
        write( 6,strfmt) str
        write(11,strfmt) str

        Do j=1,NumJ 
            n=Jsq%n(j)
            k=Jsq%k(j)
            t=Jsq%t(j)
            X1J(n,1)= X1J(n,1)+t*X1(k)
            if (n.NE.k) X1J(k,1)= X1J(k,1)+t*X1(n)     ! (4J^2)X1
        End Do
        Do i=1,Nd
            TJ1(i)= X1J(i,1) - (j0+4)*(j0+6)*X1(i)   ! X_J^(1)
            X1J(i,1)=0.d0
        End Do

        Do j=1,NumJ 
            n=Jsq%n(j)
            k=Jsq%k(j)
            t=Jsq%t(j)
            X1J(n,1)= X1J(n,1)+t*TJ1(k)
            if (n.NE.k) X1J(k,1)= X1J(k,1)+ t*TJ1(n)
        End Do
        Do i=1,Nd
            TJ2(i)= X1J(i,1) - (j0+2)*(j0+4)*TJ1(i)   ! X_J^(2)
            X1J(i,1)= 0.d0
        End Do

        Do j=1,NumJ 
            n=Jsq%n(j)
            k=Jsq%k(j)
            t=Jsq%t(j)
            X1J(n,1)= X1J(n,1)+t*TJ2(k)
            if (n.NE.k) X1J(k,1)= X1J(k,1)+ t*TJ2(n)
        End Do
        Do i=1,Nd
            TJ3(i)= X1J(i,1) - j0*(j0+2)*TJ2(i)       ! X_J^(3)
            X1J(i,1)= 0.d0
        End Do

        If (Tj0.GE.1.5d0) Then
            Do j=1,NumJ 
                n=Jsq%n(j)
                k=Jsq%k(j)
                t=Jsq%t(j)
                X1J(n,2)= X1J(n,2)+ t*TJ3(k)
                if (n.NE.k) X1J(k,2)= X1J(k,2)+ t*TJ3(n)
            End Do
            f= -1.d0/(1536*(j0-2)*j0*(j0+1)*(j0+2))
            Do i=1,Nd
                TJ4(i)= X1J(i,2) - (j0-4)*(j0-2)*TJ3(i)   ! X_J^(4)
                X1J(i,2)= f*TJ4(i)                        ! X_J-1
            End Do
        Else
            X1J(1:Nd,2)= 0.d0
        End If

        If (Tj0.GE.2.d0) Then
            Do j=1,NumJ 
                n=Jsq%n(j)
                k=Jsq%k(j)
                t=Jsq%t(j)
                X1J(n,1)= X1J(n,1)+t*TJ3(k)
                if (n.NE.k) X1J(k,1)= X1J(k,1)+ t*TJ3(n)
            End Do
            f= 1.d0/(6144*(j0-2)*(j0-1)*j0*(j0+1))
            Do i=1,Nd
                TJ5(i)= X1J(i,1) - (j0-2)*j0*TJ3(i)   ! X_J^(5)
                X1J(i,1)= f*TJ5(i)                    ! X_J-2
            End Do
        Else
            X1J(1:Nd,1)= 0.d0
        End If

        f= 1.d0/(2*(j0+2)*(j0+3))
        Do i=1,Nd
            X1J(i,3)= f*( TJ2(i)/16 - 6*(j0+1)*(j0+2)*X1J(i,2)  & ! X_J
                       - 12*j0*(j0+1)*X1J(i,1) )
        End Do

        f= -1.d0/(j0+4)
        Do i=1,Nd
            X1J(i,4)= f*( TJ1(i)/4 + 2*(j0+3)*X1J(i,3)  &           ! X_J+1
                       + 3*(j0+2)*X1J(i,2) + 4*(j0+1)*X1J(i,1) )
        End Do

        Do i=1,Nd
            X1J(i,5)= X1(i)- X1J(i,4)- X1J(i,3)- X1J(i,2)- X1J(i,1)  ! X_J+2
        End Do

!  Several tests:
        sum = 0.d0
        sum1= 0.d0
        sum2= 0.d0
        sum3= 0.d0
        sum4= 0.d0
        sum5= 0.d0
        XJ3xXJ4= 0.d0
        XJ3xXJ5= 0.d0
        XJ4xXJ5= 0.d0
        Do i= 1,Nd
          sum = sum  + X1(i)*X1(i)
          sum1= sum1 + X1J(i,1)**2
          sum2= sum2 + X1J(i,2)**2
          sum3= sum3 + X1J(i,3)**2
          sum4= sum4 + X1J(i,4)**2
          sum5= sum5 + X1J(i,5)**2

          XJ3xXJ4= XJ3xXJ4 + X1J(i,3)*X1J(i,4)
          XJ3xXJ5= XJ3xXJ5 + X1J(i,3)*X1J(i,5)
          XJ4xXJ5= XJ4xXJ5 + X1J(i,4)*X1J(i,5)
        End Do
        sum_tot= sum1+sum2+sum3+sum4+sum5

        strfmt = '(7E13.5)'
        write( 6,strfmt) sum,sum1,sum2,sum3,sum4,sum5,sum_tot
        write(11,strfmt) sum,sum1,sum2,sum3,sum4,sum5,sum_tot

        ds= dabs(sum/sum_tot - 1.d0)
        If (ds.GT.1.d-5) Then
            write (*,*) ' PrjE2 warning: normalization has changed by',ds
        End If

        Deallocate(TJ1,TJ2,TJ3,TJ4,TJ5)
        Return
    End Subroutine PrjE2

    Subroutine RdcX1J   !### Transforms X1J to the form which corresponds
        Use wigner      !### to the reduced ME and saves it in INE_J.XIJ.
        Implicit None   
        Integer :: i, jf, j, is
        Real(dp) :: Aj, xm, W

        open (unit=16,file='INE_J.XIJ',status='UNKNOWN',form='UNFORMATTED')
        close(16,status='DELETE')
        open (unit=16,file='INE_J.XIJ',status='NEW',form='UNFORMATTED')

        jf=3
        if (Kli.EQ.5) jf=5

        Do j=1,jf
            If (Kli.EQ.5) Then
                Aj= Tj0+(j-3)
            Else
                Aj= Tj0+(j-2)            !### Ang. momentum of this component
            End If
            Aj=Anint(Aj*100.0)/100.0
            If (Kli.NE.1) then         !### No transformation for RHS=H_p
                xm=dabs(Jm)-1.d-5
                if (Aj.GT.xm) then      !### Transformation to the reduced ME
                    is = Aj-Jm+1.d-5
                    is = 1 - 2*mod(is,2)
                    if (Kli.EQ.5) then
                      W=is/(FJ3(Aj,2.d0,Tj0,-Jm,Q,Jm0)+1.d-77)
                    else
                      W=is/(FJ3(Aj,1.d0,Tj0,-Jm,Q,Jm0)+1.d-77)
                    end if
                    if (dabs(W).GT.1.d+5) W=0.d0
                else
                    W=0.d0
                end if
                Do i=1,Nd
                    X1J(i,j)=W*X1J(i,j)
                End do
            End If
            Write(16) -E0,Aj,Nd,(X1J(i,j),i=1,Nd)
        End Do
        Close(16)
        Return
    End Subroutine RdcX1J

    Subroutine Prin     !### prints information about vectors X1, Y1 and Y2
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i, ic, j, jtj, n, ndk, k, l, ndd, is, err_stat
        Real(dp) :: xn, yn, zn, dx, dy, dz, xj, yj, zj, t, Ex, E, c1, c2, &
                    w3j, w2, al0, al2, tj
        Real(dp), Allocatable, Dimension(:) :: XXn 
        Character(Len=5), Dimension(5) :: st
        Character(Len=256) :: strfmt, err_msg
        data st /'H_pnc','E1(L)','H_am ','E1(V)',' E2  '/

        strfmt = '(3X,63("="))'
        write( 6,strfmt)
        write(11,strfmt)

        If (W0.EQ.0.d0) Then
            strfmt = '(4X,"solution of the eq-n (E-H0)|X1> = ",A5,"|X0>," &
                /4X,"where    E =",F13.6,",", &
                /10X,"|X0> = |",I2,","F13.6,",",2F4.1,">.", &
                /3X,63("-"))'
            write( 6,strfmt) st(Kli),Elft,N0,E0,Tj0,Jm0
            write(11,strfmt) st(Kli),Elft,N0,E0,Tj0,Jm0
        Else
            strfmt = '(4X,"solution of the eq-n (E+w-H0)|X1> = ",A5,"|X0>," &
                /4X,"where E+w =",F13.6," +"F10.6," =",F13.6, &
                /10X,"|X0> = |",I2,","F13.6,",",2F4.1,">.", &
                /3X,63("-"))'
            write( 6,strfmt) st(Kli),E0,W0,Elft,N0,E0,Tj0,Jm0
            write(11,strfmt) st(Kli),E0,W0,Elft,N0,E0,Tj0,Jm0
        End If

        strfmt = '(6X,"Nc",5X,"|",A5,"|X0>",8X,"|X1>",8X,"<X2|",A5,"|",/3X,63("-"))'
        write( 6,strfmt) st(kli),st(klf)
        write(11,strfmt) st(kli),st(klf)

        i=0
        xn=0.d0
        yn=0.d0
        zn=0.d0
        do ic=1,Nc
           dx=0.d0
           dy=0.d0
           dz=0.d0
           ndk=Ndc(ic)
           do k=1,ndk
              i=i+1
              dx=dx+X1(i)**2
              dy=dy+YY1(i)**2
              dz=dz+YY2(i)**2
           end do
           xn=xn+dx
           yn=yn+dy
           zn=zn+dz
        end do
        strfmt = '(3X,63("-"))'
        write(11,strfmt)

        xj=0.d0
        yj=0.d0
        zj=0.d0
        Do j=1,NumJ
            n=Jsq%n(j)
            k=Jsq%k(j)
            t=Jsq%t(j)
            If (n.EQ.k) Then 
                l=1
            Else
                l=2
            End If
            xj=xj+l*X1(k)*t*X1(n)
            yj=yj+l*YY1(k)*t*YY1(n)
            zj=zj+l*YY2(k)*t*YY2(n)
        End Do

        xj=0.5d0*(dsqrt(1.d0+xj/xn)-1.d0)
        yj=0.5d0*(dsqrt(1.d0+yj/yn)-1.d0)
        zj=0.5d0*(dsqrt(1.d0+zj/zn)-1.d0)

        strfmt = '(4X," J =",3F15.8)'
        write( 6,strfmt) yj,xj,zj
        write(11,strfmt) yj,xj,zj

        If (Kli.EQ.5) Then
            strfmt = '(3X,78("="),/4X,"n",8X,"J_n",11X,"Ev_n",5X, &
                    "<X2||",A5,"||n>",2X,"<n||",A5,"||X0>",3X, &
                    "contrib. |n>",/69X," to alpha0"/3X,78("-"))'
            write( 6,strfmt) st(klf),st(kli)
            write(11,strfmt) st(klf),st(kli)
        Else
            If (Kli.EQ.2 .AND. Klf.EQ.2) Then
                strfmt = '(3X,94("="),/4X,"n",8X,"J_n",11X,"Ev_n",5X, &
                         "<X2||",A5,"||n>",2X,"<n||",A5,"||X0>",3X, &
                         "contrib. |n>",4X,"contrib. |n>", &
                         /69X," to alpha_0",5X," to alpha_2",/3X,94("-"))'
                write( 6,strfmt) st(klf),st(kli)
                write(11,strfmt) st(klf),st(kli)
            Else
                strfmt = '(3X,63("="),/4X,"n",9X,"J_n",10X,"Ev_n",5X, &
                         "<X2#",A5,"#n>",4X,"<n#",A5,"#X0>",/3X,63("-"))'
                write( 6,strfmt) st(klf),st(kli)
                write(11,strfmt) st(klf),st(kli)
            End If
        End If

        Ex=E0
        if (Kli.EQ.2.AND.Klf.NE.2) Ex=E2

        open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.XIJ is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If
        If (.not. allocated(XXn)) Allocate(XXn(Nd))
        j=0
        jtj=0
        Do n=1,Nlv
            read(16) E,tj,ndd,(XXn(i),i=1,Nd)
            E=-E
            jtj=2*tj+1.d-3
            tj=0.5d0*jtj
            If (ndd.NE.Nd) then
               write(*,*)' CONF.XIJ: mismatch in the length of record:'
               write(*,*)' ndd =',ndd,' while Nd =',Nd
               Return
            End If
            j=j+1
            c1=0.d0                           !### c1=(e_x-e_n+W0)<n|X1>
            c2=0.d0                           !###   =(e_x-e_n+W0)<n|A|X0>/(e_x-e_n+W0)
            w3j=0.d0
            Do i=1,Nd                         !###   =<n|A|X0>.
                c1=c1+X1(i)*XXn(i)*(Ex-E+W0)  !### c2=<X2|B|n>.
                c2=c2+YY2(i)*XXn(i)
            End Do
            if (Klf.EQ.4) c2=c2/(E-E2)     !### V-gauge for E1
            if (Kli.NE.1) then             !### transformation to reduced ME
                k=tj-Jm+1.d-1
                is=1
                if (k.NE.2*(k/2)) is=-1
                If (Kli.EQ.5) Then
                    w3j=is*FJ3(tj,2.d0,Tj0,-Jm,Q,Jm0)+1.d-77
                Else
                    w3j=is*FJ3(tj,1.d0,Tj0,-Jm,Q,Jm0)+1.d-77
                End If
                c1=c1/w3j
                if (dabs(c1).GT.1.d+10) c1=0.d0
            end if
            if (Klf.NE.1) then             !### transformation to reduced ME
                k=Tj2-Jm0+1.d-1
                is=1
                if (k.NE.2*(k/2)) is=-1
                If (Klf.EQ.5) Then
                    w2=is*FJ3(Tj2,2.d0,tj,-Jm0,-Q,Jm)+1.d-77
                Else
                    w2=is*FJ3(Tj2,1.d0,tj,-Jm0,-Q,Jm)+1.d-77
                End If
                c2=c2/w2
                if (dabs(c2).GT.1.d+10) c2=0.d0
            end if
            If (Kli.EQ.5) Then
                if (dabs(E-Ex).LT.1.d-5) then
                    if (W0.EQ.0.d0) then
                        al0= 0.d0
                    else
                        al0= 2.d0/(5*(2*Tj0+1))* c1**2/(E-Ex-W0)   ! contrib-n of |n> to alphaE2_0
                    end if
                else
                    al0= 2.d0/(5*(2*Tj0+1))* c2**2/(E-Ex-W0)
                end if
                strfmt = '(3X,I2,F15.9,F13.6,3(1pE16.6))'
                write( 6,strfmt) j,tj,e,c2,c1,al0                 ! Should be compared with DTM
                write(11,strfmt) j,tj,e,c2,c1,al0
            End If
            If (Kli.EQ.2 .AND. Klf.EQ.2) Then
                al0= 2.d0/(3*(2*Tj0+1))*c2**2/(E-elft)                         ! contrib-n of |n> to alphaE1_0
                is= tj+Tj0+0.1d0
                is= 1-2*mod(is,2)
                al2= 4*dsqrt(5*Tj0*(2*Tj0-1)/(6*(2*Tj0+3)*(2*Tj0+1)*(Tj0+1))) & ! contrib-n of |n> to alphaE1_2
                     *FJ6(Tj0,1.d0,tj, 1.d0,Tj0,2.d0)* is* c2**2/(e-elft)
                strfmt = '(3X,I2,F15.9,F13.6,4(1pE16.6))'
                write( 6,strfmt) j,tj,e,c2,c1,al0,al2                              ! Should be compared with DTM
                write(11,strfmt) j,tj,e,c2,c1,al0,al2
            Else
                if (Kli.NE.5) then
                    strfmt = '(3X,I2,F15.9,F13.6,2E16.6)'
                    write( 6,strfmt) j,tj,e,c2,c1      !### Should be compared with DTM
                    write(11,strfmt) j,tj,e,c2,c1
                end if
            End If
        End Do

        If (Kli.EQ.5) Then
            strfmt = '(3X,78("="))'
            write( 6,strfmt)
            write(11,strfmt)
        End If

        IF (Kli.EQ.2 .AND. Klf.EQ.2) THEN
            strfmt = '(3X,94("="))'
            write( 6,strfmt)
            write(11,strfmt)
        ELSE
            If (Kli.NE.5) Then
                strfmt = '(3X,63("="))'
                write( 6,strfmt)
                write(11,strfmt)
            End If
        END IF

        Return
    End Subroutine Prin

    Subroutine RdcE1(n)                 !### Polarizability (Kli=Klf=2)
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: n, i, isk, id, k, kmin, kmax, is
        Real(dp) :: f0, tj, w3j, f, f1, f2, al, al0, al1, al2, alp
        Real(dp), Dimension(3) :: sk0, sk1, sk2
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        If (.not. Allocated(Y2J)) Allocate(Y2J(Nd,5))

        if (dabs(Tj2-Tj0).GT.1.d-5) then
           write(*,*) ' can not calculate polarizability for J2 >< J0'
           return
        end if

        f0=-2.d0
        ss(n)= 0.d0
        strfmt = '(3X,"alpha(J=",F4.1," M=",F4.1,")=",E12.5," (no Prj)")'
        if (dabs(Q).LT.1.d-5) then        !### summation without
          do i=1,Nd                       !#### Prj decomposition
            ss(n)= ss(n)+YY2(i)*X1(i)*f0  !#### requires JM=JM0, Q=0
          end do                          !### ss = alpha(J0,M0)
          write(*,*)
          write( 6,strfmt) Tj0,Jm0,ss(n)
          write(11,strfmt) Tj0,Jm0,ss(n)
        end if

        isk=0                            !### Prj decomposition
        s0(n)= 0.d0                      !### = alpha0(J0)
        s1(n)= 0.d0                      !### = alpha1(J0)
        s2(n)= 0.d0                      !### = alpha2(J0)
        call DefSum(id,kmin,kmax)        !### defines min/max J in sums
        do k=kmin,kmax
          tj=Tj0+(k-2)                   !### intermediate momentum
          tj=Anint(tj*100.0)/100.0
          is=dabs(tj-Jm0)+1.d-5
          is=1-2*mod(is,2)
          w3j=is*FJ3(Tj0,1.d0,tj,-Jm0,-Q,Jm)
          if (w3j.NE.0.d0) then
            f=f0/w3j/(3.d0*(2*Tj0+1))
          else
            if (Tj0+tj.GT.0.99d0) isk=isk+1     !### number of skiped terms
            f=0.d0
          end if
          sk0(k)=0.d0
          do i=1,Nd
            sk0(k)=sk0(k)+f*Y2J(i,k)*X1J(i,k)
          end do
          is=tj+Tj0+1.1d0
          is=1-2*mod(is,2)
!MK          f1=-is*1.5d0*dsqrt(3.d0)*(2*Tj0+1)
          If (W0.NE.0.d0) Then
            f1=is*1.5d0*dsqrt(6*Tj0*(2*Tj0+1)/(Tj0+1)) &
               *FJ6(Tj0,1.d0,tj, 1.d0,Tj0,1.d0)
          Else
            f1= 0.d0
          End If
          sk1(k)=f1*sk0(k)

          f2=-2*is*dsqrt(15*Tj0*(2*Tj0-1)*(2*Tj0+1) &
             /(2*(2*Tj0+3)*(Tj0+1)))*FJ6(Tj0,1.d0,tj, 1.d0,Tj0,2.d0)
          sk2(k)=f2*sk0(k)

          strfmt2 = '(1X,"J =",F4.1," to alpha_0:",E12.5, &
                        ", alpha_1:",E12.5,", alpha_2:",E12.5)'
          write( 6,strfmt2) tj,sk0(k),sk1(k),sk2(k)
          write(11,strfmt2) tj,sk0(k),sk1(k),sk2(k)

          s0(n)= s0(n) + sk0(k)
          s1(n)= s1(n) + sk1(k)
          s2(n)= s2(n) + sk2(k)
        end do
        if (Tj0.GT.0.51d0) then
          s(n)= s0(n)+(3*Jm0*Jm0-Tj0*(Tj0+1))/(Tj0*(2*Tj0-1.d0))*s2(n)
        else
          s(n)= s0(n)
        end if
        ! - - - - - - - - - - Tensor polarizability - - - - - - - - - - -
        strfmt3 = '(3X,"Polarizability Alpha( J=",F4.1, &
                 " M=",F4.1,") =",E12.5," a.u. (with Prj)", &
                 /3X,"=",E12.5," +(",E12.5, & 
                 ") (3M^2-J(J+1))/J*(2J-1) a.u.", &
                 /3X,"Alpha_1 =",E12.5," a.u.")'
        if (isk.EQ.0) then
            write( 6,strfmt3) Tj0,Jm0,s(n),s0(n),s2(n),s1(n)
            write(11,strfmt3) Tj0,Jm0,s(n),s0(n),s2(n),s1(n)
        end if

        If (n.EQ.2 .or. xlamb.EQ.0_dp) Then
            alp = (ss(1)+ss(2))/2.d0
            al  = (s(1)+s(2))/2.d0
            al0 = (s0(1)+s0(2))/2.d0
            al2 = (s2(1)+s2(2))/2.d0
            al1 = (s1(1)-s1(2))

            strfmt2 = '(3X,9("-")/,3X,"In total:",/3X,9("-"))'
            write( 6,strfmt2)
            write(11,strfmt2)

            if (dabs(Q).LT.1.d-5) then
              write( 6,strfmt) Tj0,Jm0,alp
              write(11,strfmt) Tj0,Jm0,alp
            end if

            write( 6,strfmt3) Tj0,Jm0,al,al0,al2,al1
            write(11,strfmt3) Tj0,Jm0,al,al0,al2,al1

            if (ok) Then
              strfmt2 = '(/1X,"RESULT: lambda=",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  CONVERGED")'
            else
              strfmt2 = '(/1X,"RESULT: lambda=",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  DIVERGED")'
            End If
            write( 6,strfmt2) abs(xlamb),al0,al2
            write(11,strfmt2) abs(xlamb),al0,al2
            write(99,strfmt2) abs(xlamb),al0,al2
        End If

        Return
    End Subroutine RdcE1

    Subroutine RdcE2(n)   !### E2 polarizability (Kli=Klf=5)
        Use wigner, Only : FJ3
        Implicit None
        Integer :: i, isk, k, n, is
        Real(dp) :: f0, tj, w3j, f, al0, alp
        Real(dp), Dimension(5) :: sk0
        Real(dp), Dimension(2) :: s0, ss
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        f0= -2.d0

        ss(n)= 0.d0
        if (dabs(Q).LT.1.d-5) then           !#  summation without Prj decomposition
          do i=1,Nd                          !##  requires JM=JM0, Q=0
            ss(n)= ss(n)+ YY2(i)*X1(i)*f0    !### ss = alpha_2(J0,M0)
          end do
        end if

        isk=0                                !### Prj decomposition
        s0(n)= 0.d0                          !### s0 = alpha_E2(J0)
        do k=1,5
          tj=Tj0+(k-3)                       !### intermediate momentum
          tj=Anint(tj*100.0)/100.0
          If (tj.GT.dabs(Jm)-1.d-5) Then
            is=dabs(tj-Jm0)+1.d-5
            is=1-2*mod(is,2)
            w3j=is*FJ3(Tj0,2.d0,tj,-Jm0,-Q,Jm)
            if (w3j.NE.0.d0) then
              f=f0/w3j/(5.d0*(2*Tj0+1))
            else
!              if (Tj0+tj.GT.1.99d0) isk=isk+1     !### number of skipped terms
              f=0.d0
            end if
          Else
            f= 0.d0
          End If

          sk0(k)=0.d0
          do i=1,Nd
            sk0(k)=sk0(k)+f*Y2J(i,k)*X1J(i,k)
          end do
          if (sk0(k).NE.0.d0) then
            strfmt = '(1X,"con-n of J =",F5.2," to scalar alpha_E2:",E14.7)'
            write( *,strfmt) tj,sk0(k)
            write(11,strfmt) tj,sk0(k)
          end if
          s0(n)= s0(n)+sk0(k)
        end do

        write(*,*)
        if (dabs(Q).LT.1.d-5) then
            strfmt2 = '(3X,"alpha_E2(J=",F4.1," M=",F4.1,")=", &
                        1pE15.7," (no Prj)"'
          write( 6,strfmt2) Tj0,Jm0,ss(n)
          write(11,strfmt2) Tj0,Jm0,ss(n)
        end if
        !     - - - - - - - - - - Polarizability with Prj- - - - - - -
        strfmt3 = '(3X,"Scalar part of alpha_E2( J=",F4.1, &
               " M=",F4.1,") =",1pE15.7," a.u. (with Prj)")'
        write( 6,strfmt3) Tj0,Jm0,s0(n)
        write(11,strfmt3) Tj0,Jm0,s0(n)

        If (n.EQ.2) Then
          alp = (ss(1)+ss(2))/2.d0
          al0 = (s0(1)+s0(2))/2.d0
            strfmt = '(3X,9("-")/,3X,"In total:",/3X,9("-"))'
          write( 6,strfmt)
          write(11,strfmt)

          if (dabs(Q).LT.1.d-5) then
            write( 6,strfmt2) Tj0,Jm0,alp
            write(11,strfmt2) Tj0,Jm0,alp
          end if

          write( 6,strfmt3) Tj0,Jm0,al0
          write(11,strfmt3) Tj0,Jm0,al0
        End If
        Return
    End Subroutine

    Subroutine RdcPNC
        Use wigner, Only : FJ3
        Implicit None
        Integer :: i, k, is, id, kmin, kmax, kf
        Real(dp) :: gap, nuc, a, f0, f, s, ss, tj, w3j
        Real(dp), Dimension(3) :: s1
        Character(Len=9), Dimension(4) :: str 
        Character(Len=256) :: strfmt
        data str /'| H_pnc |','||E1(L)||','||H_am ||','||E1(V)||'/

        if (min(kli,klf).NE.1) return
        ! gap=G_F*alpha/(8*pi*sqrt(2)) (taken from phys.par)
        gap=DPcw/16.d0
        nuc=Am-Z+1.d-8
        call DefSum(id,kmin,kmax)        !### defines min/max J in sums
        k=Tj2-Jm0+1.d-1
        is=1
        if (k.NE.2*(k/2)) is=-1
        a=is/(FJ3(Tj2,1.d0,Tj0,-Jm0,0.d0,Jm0)+1.d-77)
        f0=a*gap*(-nuc)                           !### PNC+E1L
        if (Kli.EQ.1.AND.Klf.EQ.4) f0=f0/(E0-E2)  !### PNC+E1V

        ss=0.d0                               !### summation without
        do i=1,Nd                             !#### Prj decomposition
           ss=ss+YY2(i)*X1(i)*f0
        end do
        s=0.d0                                !### uses Prj decomposition
        do k=kmin,kmax
           kf=k-id
           tj=Tj0+(k-2)
           tj=Anint(tj*100.0)/100.0
           if (Kli.EQ.1) f=f0                 !### E1+PNC
           if (Klf.EQ.1) then
             is=tj-Jm+1.d-5
             is=1-2*mod(is,2)
             w3j=FJ3(tj,1.d0,Tj0,-Jm,Q,Jm0)
             f=is*f0*w3j                       !### PNC+E1
           end if
           s1(k)=0.d0
           do i=1,Nd
              s1(k)=s1(k)+f*Y2J(i,kf)*X1J(i,k)
           end do
           strfmt = '(3X,"<",F5.2,A9,"R(J=",F5.2,")",A9,F5.2,"> = ",E12.5," (a.u.)")'
           write( 6,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
           write(11,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
           s=s+s1(k)
        end do
        strfmt = '(3X,63("-"))'
        write( 6,strfmt)
        write(11,strfmt)

        strfmt = '(4X,"<",F4.1," ||E_pnc ||",F4.1,"> = ", &
           /4X,"= i*Q_w/(-N_nuc) * (",E12.5,") a.u. (using Prj)", &
           /4X,"= i*Q_w/(-N_nuc) * (",E12.5,") a.u. (without Prj)", &
           /4X,"(N_nuc =",I3,")")'
        write( 6,strfmt) Tj2,Tj0,s,ss,nuc
        write(11,strfmt) Tj2,Tj0,s,ss,nuc

        strfmt = '(3X,63("="))'
        write( 6,strfmt)
        write(11,strfmt)
        Return
    End Subroutine RdcPNC

    Subroutine RdcAM                 ! Phase factor changed 20/03/11 [Porsev]
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i, id, is, ix, isk, k, kf, ks, kx, kj, kmin, kmax
        Real(dp) :: gap, f0, s, tj, w2, f, Fn0, Fx0, Fn2, Fx2, c0, c, cm, c1, &
                    F20, F21, F2, w6, a0, a
        Real(dp), Dimension(3) :: s1
        Character(Len=9), Dimension(4) :: str
        Character(Len=256) :: strfmt
        data str /'| H_pnc |','||E1(L)||','||H_am ||','||E1(V)||'/

        if (Kli.NE.2.OR.Klf.NE.3) return
        ! gap=G_F*alpha/(8*pi*sqrt(2)) (taken from phys.par)
        gap=DPcw/16.d0
        call DefSum(id,kmin,kmax)        !### defines min/max J in sums
        k=Tj2-Jm0+1.1d0
        is=1
        if (k.NE.2*(k/2)) is=-1
        f0=is*gap*2*dsqrt(6.d0)          !### E1L+AM

        s=0.d0                           !### uses Prj decomposition
        isk=0
        do k=kmin,kmax
           kf=k-id
           tj=Tj0+(k-2)
           tj=Anint(tj*100.0)/100.0
           w2=FJ3(Tj2,1.d0,tj,-Jm0,-Q,Jm)
           if (w2.NE.0.d0) then
              f=f0/w2
           else
              isk=isk+1
              f=0.d0
           end if
           s1(k)=0.d0
           do i=1,Nd
              s1(k)=s1(k)+f*Y2J(i,kf)*X1J(i,k)
           end do
           strfmt = '(3X,"<",F5.2,A9,"R(J=",F5.2,")",A9, &
                F5.2,"> = ",E12.5," (a.u.)")'
           write( 6,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
           write(11,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
           s=s+s1(k)
        end do
        strfmt = '(3X,63("-"))'
        write( 6,strfmt)
        write(11,strfmt)

        Fn0 = dabs(Tj0-XInuc)
        Fx0 = Tj0+XInuc
        ix  = Fx0-Fn0+1.1d0
        Fn2 = dabs(Tj2-XInuc)
        Fx2 = Tj2+XInuc
        c0=dsqrt((XInuc+1)*(2*XInuc+1)*XInuc)
        do i=1,ix
           F0 = Fn0+i-1
           F20= F0-1
           if (F20.LT.Fn2) F20=Fn2
           F21= F0+1
           if (F21.GT.Fx2) F21=Fx2
           kx= F21-F20+1.1d0
           if(kx.GT.3)kx=3
           do k=1,kx
              F2 = F20+(k-1)
              w6 = FJ6(Tj2,F2,XInuc,F0,Tj0,1.d0)
              if (w6.NE.0.d0) then
                 c1 = c0/w6
              else
                 write(*,*) ' w6=0'
                 c1 = 1.d99
              end if
              ks = F0+Tj2+XInuc+1.1d0              !### Factor cm is
              is = 1 +2*(2*(ks/2)-ks)              !### the same as for
              cm = is*dsqrt((2*F2+1)*(2*F0+1))*w6  !### M1 amplitude
              a = 0.d0
              do kj=kmin,kmax
                 tj=Tj0+(kj-2)
                 ks = XInuc+Tj2+F2+0.1d0
                 is = 1 +2*(2*(ks/2)-ks)
                 c = is*c1 &
                      *FJ6(tj,F2,XInuc,F0,Tj0,1.d0) &
                      *FJ6(F2,XInuc,Tj2,1.d0,tj,XInuc)
                 a0 = c*s1(kj)
                 a = a + a0
!                write(*,*) '  a0 =',a0
              end do
              strfmt = '(4X,"<",F4.1,"||E_am||",F4.1,"> = ", &
                   "i*kappa/I*(",E12.5,")*(",E12.5,") [Porsev]")'
              write( 6,strfmt) F2,F0,a,cm
              write(11,strfmt) F2,F0,a,cm
           end do
        end do
        strfmt = '(3X,63("="))'
        write( 6,strfmt)
        write(11,strfmt)
        Return
    End Subroutine RdcAM

    Subroutine DefSum(id,kmin,kmax) 
        !  defines the sum over intermediate states: J=J0+(k-2), id=J2-J0
        Implicit None     
        Integer :: id, kmin, kmax, imax, ki0, ki1, kf0, kf1
        Real(dp) :: dj
        Character(Len=256) :: strfmt
                                      
        dj=Tj2-Tj0
        id=dabs(dj)+0.1d0
        if (dj.LT.0.d0) id =-id       !### Tj2=Tj0+id
        imax=2
        if (Kli.EQ.1) then            !### initial state
           ki0=2
           ki1=2
           imax=imax-1
        else
           ki0=1
           ki1=3
        end if
        if (Tj0.LT.0.99d0.AND.ki0.LT.2) ki0=2
        if (Klf.EQ.1) then            !### final state
           kf0=2+id
           kf1=2+id
           imax=imax-1
        else
           kf0=1+id
           kf1=3+id
        end if
        if (Tj2.LT.0.99d0.AND.kf0.LT.(2+id)) kf0=2+id
        kmin=max(ki0,kf0)
        kmax=min(ki1,kf1)
        if (kmax.LT.kmin.OR.iabs(id).GT.imax) then
            strfmt = '(" For Kli = ",I1," and Klf = ",I1, &
                " there is no transition ",F5.2," -->",F5.2 &
                /" kmin = ",I2," kmax = ",I2," id = ",I2)'
            write(*,strfmt) Kli,Klf,Tj0,Tj2,kmin,kmax,id
        end if
        Return
    End Subroutine DefSum

    Subroutine C_3                     !### atom-wall coefficient C_3 for X2 state.
        Use wigner, Only : FJ3
        Integer :: n, ix, jj, i 
        Real(dp) :: c, c3, aj, w
        Real(dp), Dimension(3) :: sk2
        Character(Len=256) :: strfmt

        if (Kli.NE.2.OR.Klf.NE.2) return
        strfmt = '(/3X,"C_3 coefficient for the state |",I2, &
               ",",F12.6,",",F4.1,",",F4.1,">")'
        write( 6,strfmt) N2,E2,Tj2,Jm0
        write(11,strfmt) N2,E2,Tj2,Jm0

        jj=2*Tj2+0.1d0
        c=1.d0/(12*(jj+1))
        ix=3                            !### ix is the number of nonzero
        if (jj.LE.1) ix=2               !### vectors in Y2J
        c3=0.d0                         !### C_3 coefficient
        strfmt = '(3X,"contribution of J =",F4.1,": ",E12.5," a.u.")'
        do i=1,ix
          aj=Tj2+i+1-ix                 !### angular momentum of Y2J state
          w=(FJ3(aj,1.d0,Tj2,-Jm,Q,Jm0))**2
          if (W.GT.0.d0) w=c/w          !### ME => to reduced ME
          sk2(i)=0.d0                   !### J=aj contribution to C_3
          do n=1,Nd
            sk2(i)=sk2(i)+w*Y2J(n,i+3-ix)**2
          end do
          c3=c3+sk2(i)
          write( 6,strfmt) aj,sk2(i)
          write(11,strfmt) aj,sk2(i)
        end do
        strfmt = '(3X,"C_3 =",E12.5," a.u.")'
        write( 6,strfmt) c3
        write(11,strfmt) c3
       Return
    End Subroutine C_3

    Subroutine Mxmpy(num,j,x,y)     !# Y2J(j)  = (Elft-H)*x
        Implicit None               !# Y2J(j+1)= (Elft-H)*y (only for num=2)
        Integer :: j, j1, i8, num, n, k
        Real(dp) :: t
        Real(dp), Allocatable, Dimension(:) :: x, y
        j1=j+1
        i8= 0
        Do k=1,Nd
           Y2J(k,j)=0.d0
           If (num.EQ.2) Y2J(k,j1)=0.d0
        End Do
        Do i8=1,NumH
           n=Hamil%n(i8)
           k=Hamil%k(i8)
           t=Hamil%t(i8)
           If (n.EQ.k) Then
             t= t-Elft
             If (Kdiag.EQ.1) Diag(n)=-t
           End If
           Y2J(n,j)=Y2J(n,j)-t*x(k)
           If (num.EQ.2) Y2J(n,j1)=Y2J(n,j1)-t*y(k)
           If (n.NE.k) Then
             Y2J(k,j)=Y2J(k,j)-t*x(n)
             If (num.EQ.2) Y2J(k,j1)=Y2J(k,j1)-t*y(n)
           End If
        End Do
        Kdiag=0
        Close (unit=15)
        Return
    End Subroutine Mxmpy

    Subroutine Ort
        Implicit None
        Integer :: i, j, k, err_stat, ndd
        Real(dp) :: e, s, tres, tj, sum, sp
        Real(dp), Allocatable, Dimension(:) :: XXn
        Character(Len=256) :: strfmt, err_msg
    
        tres= 1.d-11    ! X0*X1 should be < tres
        
        If (.not. allocated(XXn)) Allocate(XXn(Nd))

        open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.XIJ is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If

        Do j=1,Nlv
            Read(16,iostat=err_stat,iomsg=err_msg) E,tj,ndd,(XXn(i),i=1,Nd)
            If (err_stat /= 0) Then
                write(*,'(/A)') 'Subr. Ort:  No vector XXn, needed for'
                write(*,'(9X,A)')'orthogonalization X0, is found in CONF.XIJ '
                Stop
            End If
            If (ndd.NE.Nd) Then
                write(*,'(A,I6,A,I6)')'Subr. Ort: ndd=',ndd,' /=  Nd=',Nd
                Stop
            End If
            E=-E
            IF (dabs(Tj0-tj).LT.0.01d0 .AND. dabs(E0-E).LT.1.d-6) THEN
                k= 1
                Do
                    s= 0.d0
                    Do i= 1,Nd
                        s= s + X1(i)*XXn(i)
                    End Do
                
                    sum= 0.d0
                    Do i= 1,Nd
                        X1(i)= X1(i) - s*XXn(i)
                        sum= sum + X1(i)*X1(i)
                    End Do
                
                    sp= 0.d0
                    Do i= 1,Nd
                        sp = sp + X1(i)*XXn(i)
                    End Do
                
                    If (sp.GT.tres) Then
                        k= k+1
                        If (k.EQ.30) then
                            strfmt = '("Ort: k=30 but X1*X0=",E14.7," is still > tres=",E14.7)'
                            write(*,strfmt) sp,tres
                            Stop
                        End If
                        Cycle
                    End If
                End Do
                close(16)
                Return
            End If
        End Do   
    End Subroutine Ort

    Subroutine FormV(n)       !### Additional vectors for minimization of residue
        Implicit None
        Integer :: i, n, ndd, err_stat
        Real(dp) :: dj, tj
        Character(Len=256) :: strfmt, err_msg

        Open(unit=16,file='CONF.XIJ',status='OLD',form='UNFORMATTED',iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.XIJ is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If
        If (Kli.EQ.1) Then         !### delta J depends on the key Kli
            dj=0.1d0
        Else
            dj=1.1d0
        End If
        Do While (n.LE.IPad)
            Read(16,iostat=err_stat,iomsg=err_msg) Ev(n),tj,ndd,(X1J(i,n),i=1,Nd)
            If (err_stat /= 0) Exit
            If (tj.GT.Tj0-dj.AND.tj.LT.Tj0+dj) Then
              Ev(n)=Elft+Ev(n)
              strfmt = '(4X,"Vector ",I3," Elft+E =",F10.5," J =",F5.2)'
              write( *,strfmt) n,Ev(n),tj
              n=n+1
            End If
        End Do
        Nlft=n-1
        close(16)
        Return
    End Subroutine FormV
End Program ine