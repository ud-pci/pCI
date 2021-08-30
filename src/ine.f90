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

    Implicit None
    Integer :: i, Nd2, Nddir, nsu2, icyc
    logical :: ok

    call Init_Char(Let,Alet,Blet)
    Khe= 1   ! 1- new solution of homogeneous eq-n,
             ! 0- old solution of homogeneous eq-n.

    Nddir= 1000  ! Dimension of the matrix for initial solution by SolEq1
                 ! To solve homogeneous eq-n for the whole matrix, put Ndir=IP1

    open(unit=11,status='UNKNOWN',file='INE.RES')
    close(unit=11,status='DELETE')
    open(unit=11,status='NEW',file='INE.RES')
    write(*,'(A)')' kl= (0-new, 1-use X1, 2-use X1,Y1,Y2 ):'
    read (*,*) kl
    write(*,'(A,I2)')' kl=',kl
    N_it = 20                         !### iterations in SolEq4
    Gj = 0.d0
    call recunit
    call Input
    call Init
    if (Gj.NE.0.d0) then
      write (*,*) ' Gj =',Gj
      write (*,*) ' This code works only for Gj=0!'
      stop
    end if
    Nmax=IP4                          !### max dimension of vectors
    ok = (Kli.EQ.1.AND.Klf.EQ.2).OR.(Kli.EQ.1.AND.Klf.EQ.4).OR. &
         (Kli.EQ.2.AND.Klf.LE.3).OR.(Kli.EQ.5.AND.Klf.EQ.5)
    if (.NOT.ok) then
       write(*,*)' Unknown combination Kli =',Kli,' Klf =',Klf
       Stop
    end if
    call Rint                         !### Radial integrals
    open(17,file='CONF0.DET',status='OLD',form='UNFORMATTED')
    read (17) Nd2,nsu2                !### Number of used orbitals
    Nsu=max(nsu2,Nsu)                 !### can differ for two
    close(17)                         !### spaces!
    call Dinit                        !### Construction of the
    call Jterm                        !### basisset of determinants

    IF (Kli.EQ.5 .AND. IP1.LT.Nd) THEN
      write(*,'(A,I6,A,I6)') 'IP1=',IP1,' < Nd=',Nd
      write(*,*)' This case is not coded for E2 yet !'
      Stop
    END IF

    Int_err=0
    call Init0                        !### Evaluation of the RHS of
    call Vector(kl)                   !###  the equation and vectors Yi
    icyc= 1
    If (W0.NE.0.d0) icyc= 2

    Call ReadHIJ
    Call ReadJJJ

    Do i=1,icyc
        Ndir=Nddir
        If (Kli.EQ.5) Ndir= Nd    ! SolEq4 is not adopted yet for E2 polariz.
        IF (i.EQ.2) W0= -W0
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
              write(*,*)
              Call SolEq1(kl)
              ok=.TRUE.
            End If
        End If
        if (Kli.LE.2) Call  Prj ('  X1  ',Tj0,X1,X1J)     !### Projects X1 on J subspaces
        if (Kli.EQ.5) Call PrjE2('  X1  ',Tj0,X1,X1J)
        call RdcX1J                                       !### Transforms and saves X1J
        if (Kli.LE.2) Call  Prj ('  Y2  ',Tj2,YY2,Y2J)    !### Projects Y2 on J subspaces
        if (Kli.EQ.5) Call PrjE2('  Y2  ',Tj2,YY2,Y2J)
        call Prin                                  !### Output of the results
        If (Kli.EQ.2.AND.Klf.EQ.2) Then
          call RdcE1(i)                            !### Evaluation of E1 polarizability
          if (W0.EQ.0.d0 .OR. i.EQ.2) call C_3     !### C_3 coefficient for X2 state
        End If
        if (Kli.EQ.5) call RdcE2(i)                !### Evaluation of E2 polarizabilty
        if (Kli.EQ.1.OR.Klf.EQ.1) call RdcPNC      !### Final numbers for Q_w
        if (Klf.EQ.3) call RdcAM                   !### Final numbers for AM
        if (Int_err.NE.0) then
          write( *,'(4X,">>>> NOTE:",I7," radial integrals were absent")') Int_err
          write(11,'(4X,">>>> NOTE:",I7," radial integrals were absent")') Int_err
        end if
        if (.NOT.ok) then
          write( *,'(4X,"Convergence was not reached.")')
          write(11,'(4X,"Convergence was not reached.")')
        end if
    End Do
    Close(unit=11)
End Program ine