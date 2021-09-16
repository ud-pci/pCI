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
    Use params
    Use ine_variables
    Use determinants, Only : Dinit, Jterm
    Use ine_aux
    Use str_fmt, Only : startTimer, stopTimer

    Implicit None
    Integer :: i, Nd2, Nddir, nsu2, icyc
    Integer(Kind=int64) :: start_time
    Character(Len=16) :: timeStr
    logical :: ok

    Call startTimer(start_time)

    Call Init_Char(Let,Alet,Blet)
    Khe= 1   ! 1- new solution of homogeneous eq-n,
             ! 0- old solution of homogeneous eq-n.

    Nddir= 1000  ! Dimension of the matrix for initial solution by SolEq1
                 ! To solve homogeneous eq-n for the whole matrix, put Ndir=IP1

    Open(unit=11,status='UNKNOWN',file='INE.RES')
    Close(unit=11,status='DELETE')
    Open(unit=11,status='NEW',file='INE.RES')
    Write(*,'(A)')' kl= (0-new, 1-use X1, 2-use X1,Y1,Y2 ):'
    Read (*,*) kl
    Write(*,'(A,I2)')' kl=',kl
    N_it = 20                         !### iterations in SolEq4
    Gj = 0.d0
    Call recunit
    Call Input
    Call Init
    If (Gj.NE.0.d0) Then
      Write (*,*) ' Gj =',Gj
      Write (*,*) ' This code works only for Gj=0!'
      stop
    End If
    Nmax=IP4                          !### max dimension of vectors
    ok = (Kli.EQ.1.AND.Klf.EQ.2).OR.(Kli.EQ.1.AND.Klf.EQ.4).OR. &
         (Kli.EQ.2.AND.Klf.LE.3).OR.(Kli.EQ.5.AND.Klf.EQ.5)
    If (.NOT.ok) Then
       Write(*,*)' Unknown combination Kli =',Kli,' Klf =',Klf
       Stop
    End If
    Call Rint                         !### Radial integrals
    Open(17,file='CONF0.DET',status='OLD',form='UNFORMATTED')
    Read (17) Nd2,nsu2                !### Number of used orbitals
    Nsu=max(nsu2,Nsu)                 !### can differ for two
    Close(17)                         !### spaces!
    Call Dinit                        !### Construction of the
    Call Jterm                        !### basisset of determinants

    If (Kli.EQ.5 .AND. IP1.LT.Nd) Then
      Write(*,'(A,I6,A,I6)') 'IP1=',IP1,' < Nd=',Nd
      Write(*,*)' This case is not coded for E2 yet !'
      Stop
    End If

    Int_err=0
    Call Init0                        !### Evaluation of the RHS of
    Call Vector(kl)                   !###  the equation and vectors Yi
    icyc= 1
    If (W0.NE.0.d0) icyc= 2

    Call ReadHIJ
    Call ReadJJJ

    Do i=1,icyc
        Ndir=Nddir
        If (Kli.EQ.5) Ndir= Nd    ! SolEq4 is not adopted yet for E2 polariz.
        If (i.EQ.2) W0= -W0
        If (dabs(xlamb).GT.1.d-8) Then
            If (i.EQ.2) xlamb= -xlamb
            Write( *,'(/3X,34("-")/3X,"Calculation for lambda=",F11.2,/3X,34("-")/)') xlamb
            Write(11,'(/3X,34("-")/3X,"Calculation for lambda=",F11.2,/3X,34("-")/)') xlamb
        Else
            Write( *,'(/3X,22("-")/3X,"Calculation for W0 = 0",/3X,22("-")/)')
            Write(11,'(/3X,22("-")/3X,"Calculation for W0 = 0",/3X,22("-")/)')
        End If
        Call SolEq1(kl)                   !### Direct solution
        If (Ndir.LT.Nd) Then
            Call SolEq4(ok)                 !### Iterative solution
            If (Nd.LE.IP1 .AND. .NOT.ok) Then
              Ndir= Nd
              Write(*,*)
              Call SolEq1(kl)
              ok=.TRUE.
            End If
        End If
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
    Close(unit=11)

    Call stopTimer(start_time, timeStr)
    write(*,'(2X,A)'), 'TIMING >>> Total computation time of ine was '// trim(timeStr)
    
End Program ine