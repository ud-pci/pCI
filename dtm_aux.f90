Module dtm_aux

  Use dtm_variables

  Implicit None

  Contains

    Subroutine Init_Char(Let,Alet,Blet,yes,chm1)
      Implicit None

      Character(Len=1), Dimension(6) :: Let
      Character(Len=4), Dimension(13) :: Alet
      Character(Len=4), Dimension(5) :: Blet
      Character(Len=4), Dimension(2) :: yes*3, chm1

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

      yes(1)= 'No '
      yes(2)= 'Yes'

      chm1(1)='NRel'
      chm1(2)=' Rel'
      Return
    End Subroutine Init_Char

    Subroutine OpenFS(inam,nam,iform,kan,ntype)
      Implicit None
      Integer :: i, inam, iform, kan, ntype, imax
      Character(Len=1), Dimension(12) :: fname, nam
      Character(Len=1) :: space
      Character(Len=12) :: fnam
      equivalence (fnam, fname)
      data space/' '/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     NTYPE = 0 - OLD; NTYPE = 1 - NEW, UNKNOWN
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      Do i=1,inam
         fname(i)=nam(i)
      End Do
      Do i=inam+1,12
         fname(i)=space
      End Do
      imax=12
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      If (ntype /= 0) Then
        If (iform == 0) Open(unit=kan,file=fnam,status='unknown',err=710)
        If (iform == 1) Open(unit=kan,file=fnam,status='unknown', &
              form='unformatted',err=710)
      else
        If (iform == 0) Open(unit=kan,file=fnam,status='old',err=700)
        If (iform == 1) Open(unit=kan,file=fnam,status='old', &
              form='unformatted',err=700)
      End If
      Return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
 700  Write( 6,'(/" NO FILE ",12A1)') (fname(i),i=1,imax)
      Write(11,'(/" NO FILE ",12A1)') (fname(i),i=1,imax)
      stop
!     - - - - - - - - - - - - - - - - - - - - - - - - -
 710  Write( 6,'(/" UNABLE TO OPEN FILE:",1X,12A1)') (fname(i),i=1,imax)
      Write(11,'(/" UNABLE TO OPEN FILE:",1X,12A1)') (fname(i),i=1,imax)
      stop
    End Subroutine OpenFS

    Subroutine Input       ! Reads file CONF.INP
      Implicit None
      Integer :: i, i1, i2, ic, istr, nx, ny, nz, ne0
      Real(dp) :: x
      Character(Len=1), Dimension(16) :: name
      Character(Len=4), Dimension(2) :: yes*3
      !     - - - - - - - - - - - - - - - - - - - - - - - - -
      Trd=1.d-10
      Kdm=0
      K_M1=2
      Open(unit=99,file='dtm.in')
        Read (99,*) Kl1
        Select Case(Kl1)
        Case(1)
          Read (99,*) nterm1, nterm2
        Case(2)
          Read (99,*) nterm1, nterm2, nterm2f
        End Select
      Close(99)
      Write (*,'(/4X,"Program DTM",/4X,"Cutoff parameter :",E8.1, &
                 /4X,"Full RES file - ",A3,/4X,"DM0.RES file - ",A3, &
                 /4X,"Do you want DM (1) OR TM (2)? ",I1)') Trd, yes(Kl+1), yes(Kdm+1), Kl1
      If ((Kl1-1)*(Kl1-2) /= 0) stop
      Call Init_Char(Let, Alet, Blet, yes, chm1)
      If (Kl1 == 1) Then   ! regime of Density matrix & expectation values
        Call OpenFS(6,'DM.RES',0,11,1)
        Iprt=+1      !### parity of the transition
      else                 ! regime of Transition matrix & amplitudes
        Call OpenFS(6,'TM.RES',0,11,1)
      End If
      !     Input from the file  'CONF.INP'
      Call OpenFS(8,'CONF.INP',0,10,0)
      Read (10,'(1X,16A1)') name
      If(Kl1 == 1) Write( 6,'(/4X,"DTM: Density matrices for ",16A1)') name
      If(Kl1 == 1) Write(11,'(/4X,"DTM: Density matrices for ",16A1)') name
      If(Kl1 == 2) Write( 6,'(/4X,"DTM: Transition matrices for ",16A1)') name
      If(Kl1 == 2) Write(11,'(/4X,"DTM: Transition matrices for ",16A1)') name
      Read (10,'(5X,F5.1)') Z
      Read (10,'(5X,F5.1)') Am
      Read (10,'(5X,F5.1)')
      Read (10,'(5X,F5.1)') Jm
      Read (10,'(5X,I6)') Nso
      Read (10,'(5X,I6)') Nc
      Read (10,'(5X,I6)') Kv
      Read (10,'(5X,I6)') Nlv
      Read (10,'(5X,I6)') Ne
!    - - - - - - - - - - - - - - - - - - - - - - - - -
      Allocate(Qnl(100000000)) ! Nso+5*Nc = upper bound
      Gnuc=1.d0                  !#
      Qnuc=1.d0                  !#  Optional
      Kout=0                     !#  parameters
      Kecp=0                     !#
      Gj  =0.d0                  !#
100   Call inpstr(istr)          !#
      If (istr /= 1) goto 100    !#
      Kl=Kout                    !# Kl is used in DTM instead of Kout
!    - - - - - - - - - - - - - - - - - - - - - - - - -
      If (Nso /= 0) Then
         Read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso)
      End If
      i1=Nso+1
      Do ic=1,Nc
         ne0=0
200      i2=i1+5
         Read (10,'(6(4X,F7.4))') (Qnl(i),i=i1,i2)
         Do i=i1,i2
            x=dabs(Qnl(i))+1.d-9
            If (x < 1.d-8) goto 210
            nx=10000*x
            ny=100*x
            nz=(nx-100*ny)
            ne0=ne0+nz
         End Do
210      i2=i-1
         i1=i2+1
         If (ne0 < Ne) goto 200
         If (ne0 > Ne) Then
            Write(6,'(" INPUT: too many electrons for ic =",I4)') ic
            stop
         End If
      End Do
      Nsp=i2
      Close(unit=10)
      If(Am < 1.d0) Then
        Write(6,*) ' Give nuclear parameter A: '
        Read(*,*) Anuc
      else
        Anuc=Am
      End If
      Write( 6,'(4X,"Anuc=",F6.1,", Gnuc =",F10.5,", Qnuc =",F10.5)') Anuc,Gnuc,Qnuc
      Write(11,'(4X,"Anuc=",F6.1,", Gnuc =",F10.5,", Qnuc =",F10.5)') Anuc,Gnuc,Qnuc
      Return
    End Subroutine Input
    
    Subroutine inpstr(istr)
      Implicit None
      Character(Len=1) :: txt(5)
      Character(Len=128) :: string
      Integer     :: i, istr
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      istr=0
      string(1:5) = '     '
      txt(1) = ' '
      txt(2) = ' '
      txt(3) = ' '
      txt(4) = ' '
      txt(5) = ' '
      Read(10,'(5a1,a)') (txt(i),i=1,5), string
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 200
      If (txt(2) /= 'b'   .and.   txt(2) /= 'B') goto 200
      If (txt(3) /= 'r'   .and.   txt(3) /= 'R') goto 200
      If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 200
      Read (string,*) kbrt
      Write( *,'(5a1,i6)') (txt(i),i=1,5),kbrt
      Write(11,'(5a1,i6)') (txt(i),i=1,5),kbrt
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
200   If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 210
      If (txt(2) /= '_') goto 210
      If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 210
      If (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 210
      Read (string,*) K_is
      Write( *,'(5a1,i6)') (txt(i),i=1,5),K_is
      Write(11,'(5a1,i6)') (txt(i),i=1,5),K_is
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
210   If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 220
      If (txt(2) /= '_') goto 220
      If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 220
      If (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 220
      Read (string,*) C_is
      Write( *,'(5a1,F6.3)') (txt(i),i=1,5),C_is
      Write(11,'(5a1,F6.3)') (txt(i),i=1,5),C_is
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
220   If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 230
      If (txt(2) /= 'l'   .and.   txt(1) /= 'L') goto 230
      If (txt(3) /= 'o'   .and.   txt(3) /= 'O') goto 230
      If (txt(4) /= 'w'   .and.   txt(4) /= 'W') goto 230
      Read (string,*) Klow
      Write( *,'(5a1,i6)') (txt(i),i=1,5),Klow
      Write(11,'(5a1,i6)') (txt(i),i=1,5),Klow
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
230   If (txt(1) /= ' ') goto 240
      If (txt(2) /= 'g'   .and.   txt(2) /= 'G') goto 240
      If (txt(3) /= 'j'   .and.   txt(3) /= 'J') goto 240
      If (txt(4) /= ' ') goto 240
      Read (string,*) Gj
      Write( *,'(5a1,F6.3)') (txt(i),i=1,5),Gj
      Write(11,'(5a1,F6.3)') (txt(i),i=1,5),Gj
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
240   If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 250
      If (txt(2) /= 'o'   .and.   txt(2) /= 'O') goto 250
      If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 250
      If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 250
      Read (string,*) kout
      Write( *,'(5a1,i6)') (txt(i),i=1,5),kout
      Write(11,'(5a1,i6)') (txt(i),i=1,5),kout
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
250   If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 260
      If (txt(2) /= 'e'   .and.   txt(2) /= 'E') goto 260
      If (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 260
      If (txt(4) /= 'p'   .and.   txt(4) /= 'P') goto 260
      Read (string,*) kecp
      Write( *,'(5a1,i6)') (txt(i),i=1,5),kecp
      Write(11,'(5a1,i6)') (txt(i),i=1,5),kecp
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
260   If (txt(1) /= 'g'   .and.   txt(1) /= 'G') goto 270
      If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 270
      If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 270
      If (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 270
      Read (string,*) gnuc
      Write( *,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
      Write(11,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
270   If (txt(1) /= 'q'   .and.   txt(1) /= 'Q') goto 280
      If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 280
      If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 280
      If (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 280
      Read (string,*) qnuc
      Write( *,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
      Write(11,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
280   If (txt(2) /= 'k'   .and.   txt(2) /= 'K') goto 290
      If (txt(3) /= 'l'   .and.   txt(3) /= 'L') goto 290
      If (txt(4) /= '4') goto 290
      Read (string,*) kl4
      Write( *,'(5a1,i6)') (txt(i),i=1,5),kl4
      Write(11,'(5a1,i6)') (txt(i),i=1,5),kl4
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
290   If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 300
      If (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 300
      If (txt(4) /= '4') goto 290
      Read (string,*) nc4
      Write( *,'(5a1,i6)') (txt(i),i=1,5),nc4
      Write(11,'(5a1,i6)') (txt(i),i=1,5),nc4
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
300   If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 310
      If (txt(2) /= 'r'   .and.   txt(2) /= 'R') goto 310
      If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 310
      If (txt(4) /= '4') goto 310
      Read (string,*) crt4
      Write( *,'(5a1,F8.5)') (txt(i),i=1,5),crt4
      Write(11,'(5a1,F8.5)') (txt(i),i=1,5),crt4
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
310   If (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 320
      If (txt(2) /= '_') goto 320
      If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 320
      If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 320
      Read (string,*) n_it
      Write( *,'(5a1,i6)') (txt(i),i=1,5),n_it
      Write(11,'(5a1,i6)') (txt(i),i=1,5),n_it
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
320   If (txt(1) /= 'a'   .and.   txt(1) /= 'A') goto 330
      If (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 330
      If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 330
      If (txt(4) /= 'o'   .and.   txt(4) /= 'O') goto 330
      Read (string,*) kautobas
      Write( *,'(5a1,i6)') (txt(i),i=1,5),kautobas
      Write(11,'(5a1,i6)') (txt(i),i=1,5),kautobas
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
330   If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 340
      If (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 340
      If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 340
      If (txt(4) /= '0') goto 340
      Read (string,*) Cut0
      Write( *,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
      Write(11,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
340   If (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 350
      If (txt(2) /= 'c'   .and.   txt(2) /= 'C') goto 350
      If (txt(3) /= 'p'   .and.   txt(3) /= 'P') goto 350
      If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 350
      Read (string,*) Ncpt
      Write( *,'(5a1,i6)') (txt(i),i=1,5),Ncpt
      Write(11,'(5a1,i6)') (txt(i),i=1,5),Ncpt
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
350   If (txt(2) /= ' '   .and.   txt(2) /= '-') goto 700
      If (txt(3) /= ' '   .and.   txt(3) /= '-') goto 700
      If (txt(4) /= ' '   .and.   txt(4) /= '-') goto 700
      If (txt(2) == ' '   .and.   txt(3) == ' ') backspace(10)
      istr=1
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
700   istr=2
      Write( *,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
      Write(11,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
      Read(*,*)
      stop
    End Subroutine inpstr

    Subroutine Init   ! Reads the head of the file CONF.DAT
      Implicit None
      Integer :: i, ni, n0, nmin, imax, kkj, jjj, llj, nnj, nj, klag, &
                 If, n1, n2, l, j, n, k, ic, i0, qi, nsu2
      Real(dp) :: C1, C2, Z1, r1, r2, rmax, Bt, Al, d
      Real(dp), Dimension(IP6)  :: p, q, p1, q1 
      Real(dp), Dimension(4*IP6):: PQ
      Character(Len=1) :: let(6)
      logical :: longbasis
      Integer, Dimension(4*IPs) :: IQN
      Real(dp), Dimension(IPs)  :: Qq1
      equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
      equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), &
              (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
      data let/'s','p','d','f','g','h'/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      C1=0.01d0
      Cl=DPcl
      Mj=2*dabs(Jm)+0.01d0
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      Open (12,file='CONF.DAT',status='OLD', &
           access='DIRECT',recl=2*IP6*IPmr,err=700)
      Call ReadF(12,1,P,Q,2)
      Call ReadF(12,3,P1,Q1,2)
      Z1  =PQ(1)
      If (dabs(Z-Z1) > 1.d-6) Then
        Write( 6,'(/2X,"Nuclear charge is changed"/2X,"Z1=",F12.6/2X,"Z2=",F12.6)') Z,Z1
        Write(11,'(/2X,"Nuclear charge is changed"/2X,"Z1=",F12.6/2X,"Z2=",F12.6)') Z,Z1
        stop
      End If
      NS  =PQ(2)+C1
      II  =PQ(3)+C1
      r1  =PQ(4)
      r2  =PQ(5)
      rmax=dabs(R2)
      H   =PQ(6)
      BT  =PQ(7)
      AL  =PQ(8)
      KT  =PQ(9)+C1
      NG  =PQ(10)+C1
      RNUC=PQ(13)
      KLAG=PQ(17)+C1
      longbasis=dabs(PQ(20)-0.98765d0) < 1.d-6
      WRITE( 6,5) KL,KV,KT,Z,JM,II,NG,KLAG,H, &
           AL,BT,NSP,NS,NSO,NC,RMAX,R1,RNUC
      WRITE(11,5) KL,KV,KT,Z,JM,II,NG,KLAG,H, &
           AL,BT,NSP,NS,NSO,NC,RMAX,R1,RNUC
5     format(4X,'KL  =',I3,7X,'KV  =',I3,7X,'KT  =',I3, &
            /4X,'Z   =',F6.2,4X,'JM  =',F6.2,4X,'II  =',I4, &
            /4X,'NG  =',I3,7X,'LAG =',I3,7X,'H   =',F7.4, &
            /4X,'AL  =',F7.4,3X,'BT  =',F6.2, &
            /4X,'NSP =',I5,5X,'NS  =',I4,6X,'NSO =',I3, &
             7X,'NC =',I6, &
            /4X,'R2 =',F6.2,2X,'R1 =',E11.4,2X,'RNUC =',E11.4)
      Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
      If (longbasis) Then
        Write( *,*) ' Using variant for long basis '
        Write(11,*) ' Using variant for long basis '
        Do ni=1,Ns
          Nn(ni)=Iqn(4*ni-3)
          Ll(ni)=Iqn(4*ni-2)
          Kk(ni)=Iqn(4*ni-1)
          Jj(ni)=Iqn(4*ni)
        End Do
      else
        If=20
        Do ni=1,Ns
           If=If+1
           Nn(ni)=Pq(If)+C1
           If=If+1
           LL(ni)=PQ(If)+C1
           If=If+3
           C2=dsign(C1,PQ(If))
           Kk(ni)=Pq(If)+C2
           If=If+1
           C2=dsign(C1,PQ(If))
           Jj(ni)=Pq(If)+C2
        End Do
      End If
      Do nj=1,Nsp
        i=dsign(1.d0,Qnl(nj))
        d=dabs(Qnl(nj))+1.d-14
        d=10.0*d
        nnj=d
        d=10.0d0*(d-nnj)
        llj=d
        jjj=2*llj+i
        kkj=-i*((jjj+1)/2)
        d=100.0d0*(d-llj)
        Nq(nj)=d+0.1d0
        Do i=1,Ns
          If (nnj == Nn(i) .and. Kk(i) == kkj) Then
            Nip(nj)=i
            exit
          else If (i == Ns) Then
            Write( 6,'(/2X,"no function for n =",I3," k =",I3)') nnj,kkj
            Write(11,'(/2X,"no function for n =",I3," k =",I3)') nnj,kkj
            stop
          End If
        End Do
      End Do
      Call ReadF(12,2,R,V,2)
      Nec=0
      If (Nso /= 0) Then
        Do ni=1,Nso
           Nec=Nec+Nq(ni)
        End Do
      End If
      n0=0
      nmin=Nso+1
      If (nmin <= Nsp) Then
        Do ni=nmin,Nsp
           n0=n0+Nq(ni)
        End Do
      End If
      Ne=n0/Nc
      Nst=0
      Do ni=1,Ns
         imax=2*Jj(ni)+1
         Do j=1,imax,2
            Nst=Nst+1
         End Do
      End Do
      Write( 6,'(4X,"NE  =",I4,6X,"NEC =",I4,6X,"NST =",I7)') NE,NEC,NST
      n=0
      ic=0
      i0=0
      i=0
      nmin=Nso+1
      Do ni=nmin,Nsp
        i=i+1
        n=n+Nq(ni)
        If (n < NE) cycle
        ic=ic+1
        If (n > NE) Then
          Write( 6,'(/2X,"Wrong number of electrons "/2X,"for ICONF =",I4/)') IC
          Write(11,'(/2X,"Wrong number of electrons "/2X,"for ICONF =",I4/)') IC
          stop
        End If
        Nvc(ic)=i
        Nc0(ic)=Nso+i0
        i0=i0+i
        n=0
        i=0
      End Do
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      Do ic=1,Nc
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)
        Do i=n1,n2
          ni=Nip(i)
          l =Ll(ni)+1
          j =Jj(ni)
          n =Nn(ni)
          qi=Nq( i)
          If (Nq(i) > j+1) Then 
            Write( 6,'(/2X,"Wrong number of electrons "/ &
               2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
               " (",F6.3,")")') ni,n,let(l),k,qi
            Write(11,'(/2X,"Wrong number of electrons "/ &
               2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
               " (",F6.3,")")') ni,n,let(l),k,qi
           stop
          End If
        End Do
      End Do
      Close(12)
      Deallocate(Qnl)
      Open(17,file='CONF.DET',status='OLD',form='UNFORMATTED')
      Read(17) Nd1
      Close(17)      
      Allocate(Iarr(Ne,Nd1))
      If (Kl1 == 2) Then   ! TM regime requires files CONF1.DET & CONF1.XIJ
                         ! in addition to CONF.DET & CONF.XIJ for DM regime
        Open(17,file='CONF1.DET',status='OLD',form='UNFORMATTED')
        Read(17) Nd2,nsu2                     !### Number of used orbitals
        Nsu=max(nsu2,Nsu)                     !### can differ for two
      Close(17)                             !### spaces!      Return
      End If
      Return
700   Write(*,*) ' No file CONF.DAT'
      stop
    End Subroutine Init

    Subroutine ReadF(kan,record,V1,V2,nrec)
      Implicit None
      Integer :: record, nrec, i, ii, nr1, nr2, kan
      Real(dp), Dimension(IP6) :: V1, V2
        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        Read(kan,rec=nr1) (V1(i),i=1,ii)
        If (nrec == 2) Read(kan,rec=nr2) (V2(i),i=1,ii)
      Return
    End Subroutine ReadF

    Subroutine RintA 
      ! this Subroutine takes radial integrals from file or calls Rint
      Implicit None
      Integer :: i, is, ns1, nso1, km1, ia, ib, ik, ix, nsu1
      Integer :: ns2, nso2, km2, Nint2
      Real(dp) :: z2, rn2
      Real(dp) :: z1, rn1, x
      Integer, Dimension(IPs) :: l1,l2
      Integer, Allocatable, Dimension(:) :: Intg2
      Real (dp), Allocatable, Dimension(:) :: Rnt2
      Integer, Dimension(13) :: ki, ki2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Nsu=0              !### Nsu is used to eliminate integrals
      Do i=Nso+1,Nsp
         is=Nip(i)
         If (is > Nsu) Nsu=is
      End Do
      If (Nsu >= IPx+Nso) Then
        Write(*,*)' Parameter IPx is too small for Nsu=',Nsu
        stop
      End If
      Open (13,file='DTM.INT',status='OLD',form='UNFORMATTED',err=200)
      Read (13,End=200,err=200) ns1,nso1,z1,rn1,km1
      x=iabs(ns1-Ns)+iabs(nso1-Nso)+iabs(km1-K_M1) &
          + dabs(z1-Z)+dabs(rn1-Rnuc)
      If (x < 1.d-6) Then
        Read (13) Nint,(l1(i),i=1,ns1)
        Allocate(Rnt(Nint), Intg(Nint))
        is=0
        Do i=1,Ns
          is=is+iabs(Ll(i)-l1(i))
        End Do
        If (is == 0) Then
          Read (13) (ki(i),i=1,13)
          Read (13) (Rnt(i),Intg(i),i=1,Nint)
          Write (6,'(/1X,"### Radial integrals from DTM.INT ("," Nint =",I5," ) ###", &
                 /(4X,A4," calculated by ",A4))') Nint,(Alet(i),Blet(ki(i)),i=1,13)
          Write (11,'(/1X,"### Radial integrals from DTM.INT ("," Nint =",I5," ) ###", &
                 /(4X,A4," calculated by ",A4))') Nint,(Alet(i),Blet(ki(i)),i=1,13)
          Close (13)
          nsu1=0
          Do i=1,Nint
            ix=Intg(i)
            ik=ix/(IPx*IPx)
            ia=(ix-IPx*IPx*ik)/IPx+Nso
            ib=(ix-IPx*IPx*ik-IPx*ia)+Nso
            nsu1=max(nsu1,ia,ib)
          End Do
          If (Nsu > nsu1) goto 200
         Return
        End If
      End If
      Open (13,file='DTM2.INT',status='OLD',form='UNFORMATTED',err=200)
      Read (13,End=200,err=200) ns2,nso2,z2,rn2,km2
      x=iabs(ns1-Ns)+iabs(nso1-Nso)+iabs(km1-K_M1) &
          + dabs(z1-Z)+dabs(rn1-Rnuc)
      If (x < 1.d-6) Then
        Read (13) Nint2,(l2(i),i=1,ns2)
        Allocate(Rnt2(Nint2), Intg2(Nint2))
        is=0
        Do i=1,Ns
          is=is+iabs(Ll(i)-l1(i))
        End Do
        If (is == 0) Then
          Read (13) (ki2(i),i=1,13)
          Read (13) (Rnt2(i),Intg2(i),i=1,Nint2)
          Close (13)
        End If
      End If
      
      Return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
200   Open(13,file='DTM.INT',status='UNKNOWN',form='UNFORMATTED')
      Close(13,status='DELETE')
      Call Rint
      Open(13,file='DTM.INT',status='NEW',form='UNFORMATTED')
      Write(13) Ns,Nso,Z,Rnuc,K_M1
      Write(13) Nint,(Ll(i),i=1,Ns)
      Write(13) (1,i=1,13)
      Write(13) (Rnt(i),Intg(i),i=1,Nint)
      Write(6, '(/5X,"### Radial integrals saved in DTM.INT ###")')
      Write(11,'(/5X,"### Radial integrals saved in DTM.INT ###")')
      Close(13)
      Return
    End Subroutine RintA

    Subroutine calcNint
      ! this Subroutine calculates the number of radial integrals for one-electron operators
      Implicit None
      Integer :: na, nb, ja, jb, la, lb, jab, ip, nmin, cnt
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      nmin=Nso+1
      cnt=0
      ! RADIAL INTEGRALS:
      Do na=nmin,Nsu                   !### - final state
        la=Ll(na)
        ja=Jj(na)
        Do nb=na,Nsu                  !### - initial state
          lb=Ll(nb)
          jb=Jj(nb)
          jab=iabs(ja-jb)
          ip=mod(iabs(la-lb),2)
          If (ip == 1) goto 220
          ! POSITIVE PARITY:
          If (jab > 2) goto 210
          ! Non-rel. M1-amplitude  (K_M1=1)
          If (K_M1 == 1 .and. la == lb) Then
            If (na == nb .or. ja /= jb) Then
              cnt=cnt+1
            End If
          End If
          ! M1-amplitude  (K_M1 /= 1)
          If (K_M1 /= 1) Then
            cnt=cnt+1
          End If
          ! DIPOLE HFS:
          cnt=cnt+1
210       If (jab > 4 .or. ja+jb < 4) goto 211
          ! QUADRUPOLE HFS:
          cnt=cnt+1
          ! E2 AMPLITUDE:
          cnt=cnt+1
211       If (jab > 6 .or. ja+jb < 6) cycle
          ! M3 AMPLITUDE:
          cnt=cnt+1
          cycle
          ! NEGATIVE PARITY:
220       If (jab > 6) cycle
          If (ja+jb < 6) goto 230
          ! E3 AMPLITUDE:
          cnt=cnt+1
230       If (jab > 4) cycle
          If (ja+jb < 4) goto 240
          ! MQM AMPLITUDE:
          cnt=cnt+1
          ! M2 AMPLITUDE:
          cnt=cnt+1    
240       If (iabs(ja-jb) > 2) cycle
          ! E1 AMPLITUDE (L GAUGE):
          cnt=cnt+1
          ! E1 AMPLITUDE (V GAUGE):
          cnt=cnt+1
          ! ANAPOLE MOMENT
          cnt=cnt+1
          ! EDM OF THE ELECTRON:
          If (ja == jb) Then
            cnt=cnt+1
            ! PNC AMPLITUDE:
            cnt=cnt+1
          End If
        End Do
      End Do
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Nint=cnt
      Return
    End Subroutine calcNint

    Subroutine Rint     
      Use wigner
      ! this Subroutine calculates radial integrals for one-electron operators
      Implicit None
      Integer :: kab, na, nb, ja, jb, la, lb, i, in, nmin, ih, k, &
                lg, lbs, is, las, ip, jab, n12
      Real(dp) :: x, c1, c2, c3, s1, s2, w1, w2, w3, ga, gb, gab, dn, tab, Alfd, qe
      Real(dp), Dimension(IP6)  :: p, q, a, b, ro
      Real(dp), Dimension(10) :: rcut
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      MaxT=9     !## Max power in Taylor expansion at the origin
      If (Rnuc == 0.d0) Rnuc = 1.2D-13/0.529D-8*Anuc**0.33333d0
      !  >>>>>>>>>>>>>> ECP of the core <<<<<<<<<<<<<<<<<<
      If (Kecp /= 0) Then
         Open(unit=10,file='CONF.ECP',status='OLD')
         Read (10,'(F7.3)') Alfd
         Read (10,'(F7.3)') (rcut(i),i=1,10)
         Write(11,'(5X,"ECP: Alfd =",F7.3,5X,"Rcut:" &
              /5(I2,":",F7.3))') Alfd,(i-1,rcut(i),i=1,10)
         Close(unit=10)
      End If
      !  >>>>>>>>>>>>>>    End of ECP   <<<<<<<<<<<<<<<<<<
      Nint=0
      Call calcNint
      dn=0.d0
      If (Kl == 1) Write(11,'(4X,"===== RADIAL INTEGRALS =====")')
      Open(unit=12,file='CONF.DAT',status='OLD',access='DIRECT',recl=2*IP6*IPmr)
      If (.not. Allocated(Rnt)) Allocate(Rnt(Nint))
      If (.not. Allocated(Intg)) Allocate(Intg(Nint))
      Nint=0
      Call ReadF(12,2,R,V,2)
      ih=2-Kt
      nmin=Nso+1
      ! CORE ELECTRON DENSITY AND CHARGE(R):
      Do i=1,Ii,ih
        C(i)=0.d0
      End Do
      gab=10.d0
      Do in=1,Nso
        qe=Jj(in)+1
        n12=in+4
        Call ReadF(12,n12,p,q,2)
        ga=(p(Ii+4))*2                ! changed on 1/11/11
        If (gab > ga) gab=ga
        Do i=1,Ii,ih
           C(i)=C(i)+((p(i))**2 +(q(i))**2)*qe
        End Do
      End Do
      ro(1)=Z-C(1)*R(1)/(gab+1)
      Do i=2,Ii,ih
        ro(i)=ro(i-1)-(C(i-1)+C(i))/2*(R(i)-R(i-1))
      End Do
      If (Kl == 1) Write(11,'(12F5.1)')(ro(i),i=1,Ii,10)
      ! RADIAL INTEGRALS:
      Do na=nmin,Nsu                   !### - final state
        la=Ll(na)
        ja=Jj(na)
        Call ReadF (12,na+4,p,q,2)
        ga=p(Ii+4)
        Do nb=na,Nsu                  !### - initial state
          lb=Ll(nb)
          jb=Jj(nb)
          Call ReadF (12,nb+4,a,b,2)
          gb=a(Ii+4)
          gab=ga+gb
          kab=iabs(Kk(na))+iabs(Kk(nb))
          jab=iabs(ja-jb)
          ip=mod(iabs(la-lb),2)
          If (ip == 1) goto 220
          ! POSITIVE PARITY:
          If (jab > 2) goto 210
          ! Non-rel. M1-amplitude  (K_M1=1)
          If (K_M1 == 1 .and. la == lb) Then
            If (na == nb .or. ja /= jb) Then
              Do i=1,Ii,ih
                C(i)=(p(i)*a(i)+q(i)*b(i))
              End Do
              C(ii+4)=gab
              Call Sint1(tab)
              Call AddRint(9,na,nb,0.d0,tab)
            End If
          End If
          ! M1-amplitude  (K_M1 /= 1)
          If (K_M1 /= 1) Then
            Do i=1,Ii,ih
              C(i)=0.5d0*DPcl*(p(i)*b(i)+q(i)*a(i))*R(i)
            End Do
            C(ii+4)=gab+1
            Call Sint1(tab)
            Call AddRint(9,na,nb,0.d0,tab)
          End If
          ! DIPOLE HFS:
          Do i=1,Ii,ih
            C(i)=-(p(i)*b(i)+q(i)*a(i))/(R(i)**2)
          End Do
          C(ii+4)=gab-2
          Call Sint1(tab)
          If (Am >= 1.d0) Then
            Call NclInt(kab-1,kab+1,-1.d0,-1.d0,p,b,a,q,dn)
            tab=tab-Dint+dn
          End If
          Call AddRint(1,na,nb,dn,tab)
210       If (jab > 4 .or. ja+jb < 4) goto 211
          ! QUADRUPOLE HFS:
          Do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))/(R(i)**3)
          End Do
          C(ii+4)=gab-3
          Call Sint1(tab)
          If (Am >= 1.d0) Then
            Call NclInt(kab-2,kab+2,1.d0,1.d0,p,a,q,b,dn)
            tab=tab-Dint+dn
          End If
          Call AddRint(2,na,nb,dn,tab)
          ! E2 AMPLITUDE:
          Do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*(R(i)**2)
          End Do
          C(ii+4)= gab+2
          Call Sint1(tab)
          Call AddRint(10,na,nb,dn,tab)
211       If (jab > 6 .or. ja+jb < 6) cycle
          ! M3 AMPLITUDE:
          Do i=1,Ii,ih
            C(i)= -0.5d0*DPcl*(KK(na)+KK(nb))* &  ! (-) stands because P= f*r, Q = -g*r
                  (p(i)*b(i) + q(i)*a(i))* R(i)**3
          End Do
          C(ii+4)=gab+3
          Call Sint1(tab)
          Call AddRint(13,na,nb,0.d0,tab)
          cycle
          ! NEGATIVE PARITY:
220       If (jab > 6) cycle
          If (ja+jb < 6) goto 230
          ! E3 AMPLITUDE:
          Do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*(R(i)**3)
          End Do
          C(ii+4)= gab+3
          Call Sint1(tab)
          Call AddRint(11,na,nb,0.d0,tab)
230       If (jab > 4) cycle
          If (ja+jb < 4) goto 240
          ! MQM AMPLITUDE:
          Do i=1,Ii,ih
            C(i)=-(p(i)*b(i)+q(i)*a(i))/(R(i))**3
          End Do
          C(ii+4)=gab-3
          Call Sint1(tab)
          If (Am >= 1.d0) Then
            Call NclInt(kab-2,kab+2,-1.d0,-1.d0,p,b,a,q,dn)
            tab=tab-Dint+dn
          End If
          Call AddRint(8,na,nb,dn,tab)
          ! M2 AMPLITUDE:
          Do i=1,Ii,ih
            C(i)= -2/3.d0*DPcl*(KK(na)+KK(nb))*  &  ! (-) stands because P= f*r, Q = -g*r
                  (p(i)*b(i) + q(i)*a(i))* R(i)**2
          End Do
          C(ii+4)=gab+2
          Call Sint1(tab)
          Call AddRint(12,na,nb,0.d0,tab)    
240       If (iabs(ja-jb) > 2) cycle
          ! E1 AMPLITUDE (L GAUGE):
          Do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*R(i)
            If (Kecp /= 0) C(i)=C(i)*(1-Alfd &
                 /sqrt((R(i)**2+Rcut(la+1)*Rcut(lb+1))**3))
          End Do
          C(ii+4)=gab+1
          Call Sint1(tab)
          Call AddRint(3,na,nb,Dint,tab)
          ! E1 AMPLITUDE (V GAUGE):
          las=ja-la
          lbs=jb-lb
          lg=max(la,lb)
          is=1
          k=(ja+jb)/2+la+lg
          If (k /= 2*(k/2)) is=-is
          w3 = Fj6(la+0.d0,ja/2.d0,0.5d0,jb/2.d0,lb+0.d0,1.d0)
          w1 = Fj6(0.5d0,ja/2.d0,la+0.d0,jb/2.d0,0.5d0,1.d0)/w3
          w2 = Fj6(0.5d0,ja/2.d0,lb+0.d0,jb/2.d0,0.5d0,1.d0)/w3
          Do i=1,Ii,ih
            C(i)=0.d0
            If (la == lbs) C(i) = C(i) - p(i)*b(i)*w1
            If (lb == las) C(i) = C(i) - q(i)*a(i)*w2
          End Do
          C(ii+4)=gab
          Call Sint1(tab)
          If (Am >= 1.d0) Then
            s1=-w1
            s2=-w2
            If (la /= lbs) s1=0.d0
            If (lb /= las) s2=0.d0
            Call NclInt(kab+1,kab,s1,s2,p,b,q,a,dn)
            tab=tab-Dint+dn
          End If
          tab=tab*cl*is*dsqrt(6.d0/lg)
          dn=dn*cl*is*dsqrt(6.d0/lg)
          Call AddRint(6,na,nb,dn,tab)
          ! ANAPOLE MOMENT
          If (gab < 2.5d0) Then
            Dint=0.d0
            If (Am >= 1.d0) Then
              If (la == 0) Then
                s1=-3.d0
                s2=-1.d0
              else
                s1= 1.d0
                s2= 3.d0
              End If
              Call NclInt(kab-2,kab,s1,s2,p,b,q,a,dn)
            else
              c3=1.d0/3.d0
              If (la == 0) x=-(p(Ii+5)*b(Ii+5)+c3*a(Ii+5)*q(Ii+5))
              If (lb == 0) x= (a(Ii+5)*q(Ii+5)+c3*p(Ii+5)*b(Ii+5))
              dn=x*Rnuc**(gab-2)
            End If
            Call AddRint(7,na,nb,dn,dn)
          else
            Call AddRint(7,na,nb,0.d0,0.d0)
          End If
          ! EDM OF THE ELECTRON:
          If (ja == jb) Then
            Do i=1,Ii,ih
              C(i)= (q(i)*b(i))/(R(i)**2)*Ro(i)
            End Do
            C(ii+4)=gab-2
            Call Sint1(tab)
            If (Am >= 1.d0) Then
              s1=0.d0
              s2=ro(1)/R(1)
              Call NclInt(kab+2,kab+1,s1,s2,p,a,q,b,dn)
              tab=tab-Dint+dn
            End If
            Call AddRint(4,na,nb,dn,tab)
            ! PNC AMPLITUDE:
            If (gab < 2.5d0) Then
              Dint=0.d0
              If (Am >= 1.d0) Then
                s1=-3.d0
                s2=3.d0
                Call NclInt(kab-2,kab,s1,s2,p,b,q,a,dn)
              else
                x=-(p(Ii+5)*b(Ii+5)-a(Ii+5)*q(Ii+5))
                dn=x*Rnuc**(gab-2)
              End If
              Call AddRint(5,na,nb,dn,dn)
            else
              Call AddRint(5,na,nb,0.d0,0.d0)
            End If
          End If
        End Do
      End Do
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Write( 6,'(4X,"Nint=",I7," Nsu=",I3)') Nint,Nsu
      Write(11,'(4X,"Nint=",I7," Nsu=",I3)') Nint,Nsu
      Close(12)
      Return
    End Subroutine Rint

    Subroutine Sint1(DS)        
      ! Simpson integration over r (with weight function HH*V(I))
      Implicit None
      Integer :: I, IH, I0, I1, I2, I3
      Real(dp) :: R1, R2, R3, T1, T2, F1, F2, F3, F21, F32, F321, &
                  C0, C1, C2, G, P1, P2, T, Q, HH, Gam
      Real(dp), intent(out) :: DS
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Gam=C(ii+4)
      IH=2-KT
      HH=H*IH/3.d0
      I0=IH+1
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      I1=1
      I2=I1+IH
      I3=I2+IH
      R1=R(I1)
      R2=R(I2)
      R3=R(I3)
      T1=0.d0
      T2=0.d0
      IF (GAM > 5.0d0) GOTO 200
      F1=C(I1)/R1**GAM
      F2=C(I2)/R2**GAM
      F3=C(I3)/R3**GAM
      F21=(F2-F1)/(R2-R1)
      F32=(F3-F2)/(R3-R2)
      F321=(F32-F21)/(R3-R1)
      C2=F321
      C1=F21-C2*(R1+R2)
      C0=F1-C1*R1-C2*R1**2
      G=GAM+1.d0
      T2=R1**G*(C0/G+C1/(G+1.d0)*R1+C2/(G+2.d0)*R1**2)
      T1=R2**G*(C0/G+C1/(G+1.d0)*R2+C2/(G+2.d0)*R2**2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
200   Dint=t2
      P1=C(I0)*V(I0)
      P2=C( 1)*V( 1)
      I1=I0+IH
      DO I=I1,II,IH
         Q=C(I)*V(I)
         T=HH*(Q+4.d0*P1+P2)
         T=T2+T
         T2=T1
         T1=T
         P2=P1
         P1=Q
      END DO
      ds=t
      RETURN
    End Subroutine SINT1

    Subroutine AddRint(ir,na,nb,dn,tab)
      ! adds radial integrals to array Rint and Writes to RES-file
      Implicit None
      Integer :: nna, na, la, ja, nnb, lb, jb, ir, nab1, nb
      Real(dp) :: dn, tab
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      nna=Nn(na)
      la=Ll(na)
      ja=Jj(na)
      nnb=Nn(nb)
      lb=Ll(nb)
      jb=Jj(nb)
      Nint=Nint+1
      nab1=ir*IPx*IPx+IPx*(na-Nso)+(nb-Nso)
      Rnt(Nint)=tab
      Intg(Nint)=nab1
      If (Kl == 1) Then
        If (Nint == 1) Write (11,'(/2X,68("="),/2X,"Nint",3X, &
              "Type",4X,"final",3X,"init.",4X,"Nuc.Int_0",5X, &
              "Nuc.Int",5X,"Tot.Int",/2x,68("-"))')
        Write(11,'(1X,I6,3X,A4,2X,I2,A1,1X,I1,"/2",1X, &
             I2,A1,1X,I1,"/2",3E13.5)') Nint,Alet(ir), &
            nna,Let(la+1),ja,nnb,Let(lb+1),jb,Dint,dn,tab
      End If
      Return
    End Subroutine AddRint

    Subroutine NclInt(n1,n2,s1,s2,p,q,a,b,dn)
      ! radial integration inside the nucleus: R(1)**n1*int_0,1 f(x) dx
      ! where x = r/R(1)  and  f(x) = (s1*p(x)*q(x) + s2*a(x)*b(x))*x**n2
      Implicit None
      Integer :: i, i1, j, j1, m2, n1, n2
      Real(dp) :: s1, s2, dn
      Real(dp), Dimension(IP6)  :: p, q, a, b
        dn=0.d0
        Do i=1,MaxT+1
           i1=Ii+4+i
           Do j=1,MaxT+1
              j1=Ii+4+j
              m2=i+j+n2-1
              dn = dn + (s1*p(i1)*q(j1) + s2*a(i1)*b(j1))/m2
           End Do
        End Do
        dn=dn*R(1)**n1
       Return
    End Subroutine NclInt

    Subroutine BcastDMArrays(mype, npes)
      Use mpi_f08
      Implicit None
      Integer :: mype, npes, mpierr, i
      Call MPI_Bcast(e2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(tj2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      If (mype/=0) Allocate(Iarr(Ne,Nd1))
      Do i=1,Ne
          Call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      End Do
      Do i=1,nlvs
          Call MPI_Bcast(ArrB2(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      End Do
      Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      Return
    End Subroutine BcastDMArrays

    Subroutine FormDM(mype,npes)           
      Use mpi_f08
      Use str_fmt, Only : FormattedMemSize
      Use determinants, Only : Gdet, CompC, Rspq
      Use conf_variables, Only : iconf1, iconf2
      ! calculates density matrix and expectation values
      Implicit None
      Integer :: ntr, lf, n, k, i, iq, j, ju, iu, nf, is, k1, icomp, &
                 ic2, kx, n1, ic1, imin, nx, Ndpt, j1, j2, i1, iab2, &
                 imax, nn, ntrm, ntrm1, ntrms, start, End, mpierr, pgs0, pgs, pct, size
      Integer :: npes,mype
      Integer*8 :: mem, memsum
      Real(dp) :: s, bn, bk, Tj, Etrm
      Character(Len=16)     :: memStr, npesStr
      Integer, Allocatable, Dimension(:) :: idet1, idet2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne))
      Npo=Nf0(Nso+1)    !### - last position of the core orbitals
      If (mype==0) Write (*,'(1X," DM for terms .. - .. : ", &
              I1,1X,I1)') nterm1, nterm2
      ntrm=nterm1
      ntrm1=nterm2
      nlvs=ntrm1-ntrm+1
      Allocate(B1(Nd),ArrB2(Nd,nlvs),e2s(nlvs),tj2s(nlvs))
      mem = sizeof(Iarr)+sizeof(B1)+sizeof(ArrB2)
      ! Sum all the mem sizes to get a total...
      Call MPI_AllReduce(mem, memsum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
      ! ...before overwriting mem with the maximum value across all workers:
      Call MPI_AllReduce(mem, mem, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
      If (mype==0) Then
        Write(npesStr,fmt='(I16)') npes
        Call FormattedMemSize(memsum, memStr)
        Write(*,'(A,A,A)') 'DM requires approximately ',Trim(memStr),' of memory'
        Call FormattedMemSize(mem, memStr)
        Write(*,'(A,A,A,A,A)') 'DM requires approximately ',Trim(memStr),' of memory per core with ',Trim(AdjustL(npesStr)),' cores'
      End If
      If (ntrm <= 0) Return
      If (mype==0) Then
        ! Read in wavefunctions of CONF.XIJ
        Call OpenFS(8,'CONF.XIJ',1,16,0)
        If (ntrm > 1) Then
          Do ntr=1,ntrm-1
            Read (16)
          End Do
        End If
        iab2=0
        Do ntr=ntrm,ntrm1
          iab2=iab2+1
          Read (16) Etrm,tj,ndpt,(B1(i),i=1,Nd)
          ArrB2(1:Nd,iab2) = B1(1:Nd)
          e2s(iab2)=Etrm+4.d0*Gj*tj*(tj+1.d0)
          tj2s(iab2)=tj
        End Do
        Close(16)
      End If

      Call BcastDMArrays(mype,npes)

      If (mype==0) Then
        start=1
      else
        start=mype*(Nc/npes)+1
      End If
      If (mype==npes-1) Then
        End = Nc
      else
        End = (mype+1)*(Nc/npes)
      End If
      size=End-start+1
      lf=0
      iab2=0
      Do ntr=ntrm,ntrm1
        iab2=iab2+1
        lf=lf+1
        Ro(:,:)=0.d0
        B1(1:Nd)=ArrB2(1:Nd,iab2)
        Etrm=e2s(iab2)
        tj=tj2s(iab2)
        s=0.d0
        Do i=1,Nd
          s=s+B1(i)**2
        End Do
        If (dabs(s-1.d0) > 1.d-5) goto 700
        imax=0
        imin=1000
        If (mype==0) Then
            n=0
        else
            n=sum(Ndc(1:start-1))
        End If
        pgs0=size/10
        pgs=start+pgs0
        pct=0
        Do ic1=start,End
          nx = Ndc(ic1)
          Do n1=1,nx
            n=n+1
            Ndr=n
            bn=B1(n)
            If (dabs(bn) > Trd) Then
              Call Gdet(n,idet1)
              k=0
              Do ic2=1,ic1
                kx=Ndc(ic2)
                If(ic1 == ic2) kx=n1
                Call Gdet(k+1,idet2)
                Call CompC(idet1,idet2,icomp)
                If(icomp > 1) Then
                  k=k+kx
                else
                  Do k1=1,kx
                    k=k+1
                    bk=B1(k)
                    If (dabs(bn*bk) > Trd) Then
                      Call Gdet(k,idet2)
                      Call Rspq(idet1,idet2,is,nf,iu,ju,i,j)
                      If (nf == 0) Then        !### DETERMINANTS ARE EQUAL
                        nn=2
                        If (n == k) nn=1
                        Do iq=1,Ne
                          i1=idet1(iq)-Npo
                          imax=max(i1,imax)
                          imin=min(i1,imin)
                          Ro(i1,i1)=Ro(i1,i1) + bn*bk*is*nn
                        End Do
                      End If
                      If (nf == 1) Then        !### DETERMINANTS DIFFER
                        j1=i-Npo               !#### BY ONE FUNCTION
                        j2=j-Npo
                        imax=max(imax,j1,j2)
                        imin=min(imin,j1,j2)
                        Ro(j1,j2)=Ro(j1,j2)+bn*bk*is
                        Ro(j2,j1)=Ro(j1,j2)
                      End If
                    End If
                  End Do
                End If
              End Do
              If (imax > IP1) goto 710
            End If
          End Do
            If (ic1 == pgs .and. pct < 90) Then
                pct=pct+10
                Write(*,'(2X,"core ",I3," is",I3,"% done")') mype,pct
                pgs=pgs+pgs0
            else If (ic1 == End) Then
                Write(*,'(2X,"core ",I3," has completed")') mype
            End If
          End Do
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        If (mype==0) Then
          Call MPI_Reduce(MPI_IN_PLACE, Ro(1:IP1,1:IP1), IP1*IP1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        else
          Call MPI_Reduce(Ro(1:IP1,1:IP1), Ro(1:IP1,1:IP1), IP1*IP1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        End If
        Call MPI_AllReduce(imax, imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
        If (mype==0) Then
          s=0.d0
          Do i=imin,imax
            s=s+Ro(i,i)
          End Do
          If (lf == 1 .and. Kl /= 0) Write( 6,35) Ne,s,imin,imax
          If (Kl /= 0)             Write(11,35) Ne,s,imin,imax
35        format(1X,'Ne = ',I2,'; Trace(Ro) = ',F12.8, &
                5X,'Imin, Imax =',2I4)
          Call RdcDM(ntr,Etrm,Tj,imin,imax,lf)
        End If
      End Do

      Deallocate(idet1,idet2)
      Return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
700   Write (*,*) 'FormDM: norma for vector ',ntr,' is ',s
      Return
710   Write (*,*) 'FormDM: index out of range (IP1<imax):'
      Write (*,*) 'IP1=',IP1,' imax=',imax,' n=',n
      Return
    End Subroutine FormDM

    Subroutine InitTDM(mype,npes)
      Use mpi_f08
      Implicit None
      Integer :: mype, npes, mpierr
      !----------------------------------------------
      Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Kl1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nd2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Trd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(nterm1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(nterm2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(nterm2f, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      If (.not. Allocated(Jz)) Allocate(Jz(Nst))
      If (.not. Allocated(Nh)) Allocate(Nh(Nst))
      If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
      If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
      If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
      If (.not. Allocated(Rnt)) Allocate(Rnt(Nint))
      If (.not. Allocated(Intg)) Allocate(Intg(Nint))   
      Call MPI_Bcast(Nn(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Kk(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Ll(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Jj(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Nf0(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Rnt(1:Nint), Nint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(Intg(1:Nint), Nint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      Return
    End Subroutine InitTDM

    Subroutine BcastTMArrays(mype, npes)
      Use mpi_f08
      Implicit None
      Integer :: mype, npes, mpierr, i
      Call MPI_Bcast(tj1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(e2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Call MPI_Bcast(tj2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      If (mype/=0) Allocate(Iarr(Ne,Nd1))
      Do i=1,Ne
          Call MPI_Bcast(Iarr(i,1:Nd1), Nd1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          Call MPI_Bcast(Iarr2(i,1:Nd2), Nd2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      End Do
      Call MPI_Bcast(B1(1:Nd1), Nd1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      Do i=1,nlvs
          Call MPI_Bcast(ArrB2(1:Nd2,i), Nd2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      End Do
      Return
    End Subroutine BcastTMArrays

    Subroutine FormTM(mype,npes)
      Use mpi_f08
      Use str_fmt, Only : FormattedMemSize
      Use determinants, Only : Gdet, CompC, Rspq
      Use conf_variables, Only : iconf1, iconf2
      ! calculates transition matrix & amplitudes
      Implicit None
      Integer :: lf, imax, k, i, n2, n21, n22, k1, icomp, ic, &
                iq, j, ju, iu, nf, is, kx, ks, mx, kxx, ixx, j1, j2, &
                imin, i1, n, n1, ndpt, n20, jt, iab2, start, End, pgs, pgs0, pct
      Integer :: mype, npes, mpierr, size
      Integer*8 :: mem, memsum
      Real(dp) :: tj2, bn, bk, rc, rxx, s, ms, e1, e2
      Real(dp), Allocatable, Dimension(:) :: ro1
      Integer, Allocatable, Dimension(:) :: idet1, idet2
      Character(Len=16)     :: memStr, npesStr
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne))
      Npo=Nf0(Nso+1)   !### - Last position of the core orbitals
      If (mype==0) Write (*,'(1X," TM from term .. to terms .. - .. : &
              ",I1,1X,I1,1X,I1)') nterm1, nterm2, nterm2f
      n1=nterm1
      n21=nterm2
      n22=nterm2f
      nlvs=n22-n21+1
      Allocate(e2s(nlvs),tj2s(nlvs))
      Allocate(ArrB2(Nd2,nlvs))
      Allocate(B1(Nd1),B2(Nd2))
      Allocate(Iarr2(Ne,Nd2))
      If (n1 <= 0) Return
      
      If (mype==0) Then
        ! Read in wavefunctions of CONF.XIJ
        Call OpenFS(8,'CONF.XIJ',1,16,0)
        Do n=1,n1
          Read (16) e1,tj1,ndpt,(B1(i),i=1,Nd1)
          e1=e1+4.d0*Gj*tj1*(tj1+1.d0)
        End Do
        Close (16)
        jt=2*tj1+0.1
        tj1=jt/2.d0
  
        ! Read in determinants of CONF1.DET
        Call OpenFS(9,'CONF1.DET',1,17,0)
        Read (17) Nd2
        Write (*,*) 'Ne, Nd1, Nd2',Ne,Nd1,Nd2
        Do n=1,Nd2
          Read(17) (idet2(i),i=1,Ne)
          Iarr2(1:Ne,n) = idet2(1:Ne)
        End Do
        Close(17)
  
        ! Read in wavefunctions of CONF1.XIJ
        Call OpenFS(9,'CONF1.XIJ',1,16,0)
        n20=n21-1
        If (n20 /= 0) Then
          Do n2=1,n20
            Read (16)
          End Do
        End If
        iab2=0
        Do n2=n21,n22
          iab2=iab2+1
          Read (16) e2,tj2,ndpt,(B2(i),i=1,Nd2)
          ArrB2(1:Nd2,iab2) = B2(1:Nd2)
          e2s(iab2)=e2+4.d0*Gj*tj2*(tj2+1.d0)
          jt=2*tj2+0.1
          tj2s(iab2)=jt/2.d0
        End Do
        Close(16)
      End If

      mem = sizeof(Iarr)+sizeof(Iarr2)+sizeof(B1)+sizeof(B2)+sizeof(ArrB2)
      ! Sum all the mem sizes to get a total...
      Call MPI_AllReduce(mem, memsum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
      ! ...before overwriting mem with the maximum value across all workers:
      Call MPI_AllReduce(mem, mem, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
      If (mype==0) Then
        Write(npesStr,fmt='(I16)') npes
        Call FormattedMemSize(memsum, memStr)
        Write(*,'(A,A,A)') 'TM requires approximately ',Trim(memStr),' of memory'
        Call FormattedMemSize(mem, memStr)
        Write(*,'(A,A,A,A,A)') 'TM requires approximately ',Trim(memStr),' of memory per core with ',Trim(AdjustL(npesStr)),' cores'
      End If
      
      Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      Call BcastTMArrays(mype, npes)

      If (mype==0) Then
        start=1
      else
        start=mype*(Nd2/npes)+1
      End If
      If (mype==npes-1) Then
        End = Nd2
      else
        End = (mype+1)*(Nd2/npes)
      End If
      size=End-start+1
      lf=0
      iab2=0
      Do n2=n21,n22
        lf=lf+1
        iab2=iab2+1
        B2(1:Nd2)=ArrB2(1:Nd2,iab2)
        e2=e2s(iab2)
        tj2=tj2s(iab2)
        If (dabs(tj2-tj1) < 3.1d0) Then
          Ro=0.d0
          imax=0
          imin=1000
          pgs0=size/10
          pgs=start+pgs0
          pct=0
          Do n=start,End   !### - final state
            Ndr=n
            bn=B2(n)
            idet2(1:Ne)=Iarr2(1:Ne,n)
            If (n == 1) Then
              ms=0     !### - Jz for the final state
              Do i=1,Ne
                ks=idet2(i)
                ms=ms+Jz(ks)
              End Do
              Tm2=ms/2.d0
            End If
            If (dabs(bn) > Trd) Then
              k=0
              Do ic=1,Nc
                kx=Ndc(ic)
                Call Gdet(k+1,idet1)
                Call CompC(idet2,idet1,icomp)
                If (icomp > 1) Then
                  k=k+kx
                else
                  Do k1=1,kx
                    k=k+1          !### - init. state
                    bk=B1(k)
                    If (dabs(bk*bn) > Trd) Then
                      Call Gdet(k,idet1)
                      Call Rspq(idet2,idet1,is,nf,iu,ju,i,j)
                      If (nf == 0) Then
                        Do iq=1,Ne
                          i1=idet1(iq)-Npo
                          imax=max(imax,i1)
                          imin=min(imin,i1)
                          Ro(i1,i1)=Ro(i1,i1) + bn*bk*is
                        End Do
                      End If
                      If (nf == 1) Then
                        j1=j-Npo          !### - init. state
                        j2=i-Npo          !### - final state
                        imax=max(imax,j1,j2)
                        imin=min(imin,j1,j2)
                        Ro(j1,j2)=Ro(j1,j2)+bn*bk*is
                      End If
                      If (imax > IP1) Then
                        Write(*,*) imax,'= imax > IP1 =',IP1
                        stop
                      End If
                    End If
                  End Do
                End If
              End Do
            End If
            If (n == pgs .and. pct < 90) Then
                pct=pct+10
                Write(*,'(2X,"core ",I3," is",I3,"% done")') mype,pct
                pgs=pgs+pgs0
            else If (n == End) Then
                Write(*,'(2X,"core ",I3," has completed")') mype
            End If
          End Do
          Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
          Call MPI_AllReduce(imax, imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)

          Allocate(ro1(imax))
          Do i = 1,npes-1
            Do j=1,imax
              If (mype == i) Call MPI_SEND(Ro(j,1:imax),imax,MPI_DOUBLE_PRECISION,0,&
                                            0,MPI_COMM_WORLD,mpierr)
              If (mype == 0) Then
                  Call MPI_RECV(ro1(1:imax),imax,MPI_DOUBLE_PRECISION,i,&
                                0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
                  Ro(j,1:imax)=Ro(j,1:imax)+ro1(1:imax)
              End If
            End Do
          End Do
          Deallocate(ro1)
          !If (mype==0) Then
          !  Call MPI_Reduce(MPI_IN_PLACE, Ro(1:imax,1:imax), imax**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
          !                      MPI_COMM_WORLD, mpierr)
          !else
          !  Call MPI_Reduce(Ro(1:imax,1:imax), Ro(1:imax,1:imax), imax**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
          !                      MPI_COMM_WORLD, mpierr)
          !End If
          !If (mype == 0 ) Then
          !  Do i=1,imax
          !    print*,i,Ro(1,i)
          !  End Do
          !End If
          
          If (mype==0) Then 
            s=0.d0
            rxx=0.d0
            Do i=imin,imax
              s=s+Ro(i,i)
              Do k=imin,imax
                rc=dabs(Ro(i,k))
                If (rc > rxx) Then
                  ixx=i
                  kxx=k
                  rxx=rc
                End If
              End Do
            End Do
            ixx=Nh(ixx+Npo)
            kxx=Nh(kxx+Npo)
            Iprt=1-2*mod(Ll(ixx)+Ll(kxx),2) !### - parity of the transition
            If (Kl /= 0) Write( 6,'(/1X,"Ne = ",I2,"; Trace(Ro) = ", &
                          F12.8,5X," PARITY = ",I2)') Ne,s,Iprt
            If (Kl /= 0) Write(11,'(/1X,"Ne = ",I2,"; Trace(Ro) = ", &
                          F12.8,5X," PARITY = ",I2)') Ne,s,Iprt
            Call RdcTM(n1,n2,e1,e2,tj1,tj2,imin,imax,lf)
          End If
        End If
      End Do
      ! - - - - - - - - - - - - - - - - - - - - - - - -
      Deallocate(idet1,idet2,iconf1,iconf2)
      Return
    End Subroutine FormTM

    Integer function Isig(n)
      Implicit None
      Integer :: n
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        If (n == 2*(n/2)) Then
           Isig=1
        else
           Isig=-1
        End If
       Return
    End function Isig

    Integer function Isgn(n)
      Implicit None
      Integer :: n
      Isgn = 1 + 2*(2*(n/2)-n) != (-1)**n
      Return
    End function Isgn

    Real(dp) function Gj1(x,l)
      ! used for G-factors
      Implicit None
      Integer :: l
      Real(dp) :: x
      Gj1 = (x+0.5d0)*dsqrt((2*x+1)*(4*x-2*l+1)/(2*l+1)) 
      Return
    End function Gj1

    Real(dp) function Gjj(kx,l)
      ! used for G-factors
      Implicit None
      Integer :: l, kx
      Gjj = Isgn((kx+1)/2+l)*dsqrt(2.d0*l*(l+1)/(2*l+1))
      Return
    End function Gjj

    Subroutine RdcTM (n1,n2,e1,e2,tj1,tj2,imin,imax,lf)
      Use wigner
      Implicit None
      ! FORMATION OF REDUCED TRANSITION MATRICES OF RANKS 0, 1, AND 2 FROM Ro
      Integer :: k, kmax, no, i, imin1, imax1, lf, n1, n2, &
                 lk, jk, ik2, ik1, imin, imax, kx, n, l, ml, il, &
                 mk, ik, mc, lll, l1, icc, nok, il1, il2, nol, jl
      Real(dp) :: ppl, delE, q12, tm1, tj2, tj1, e1, e2, &
                 s, xjk, c, tl, AM3, AM2, AE3, AE2, QM, PNC, EDM, &
                 G, A, B, AE1, x, tme, xml, xmk, xjl, xg, EDM1, QM1, &
                 AM1, Wc, PNC1, g1
      Integer, Dimension(3*IPx) :: ind
      Integer, Dimension(IPx) :: i1, i2
      ! tj1,tj2,tm1,Tm2 - TOTAL MOMENTA AND THEIR PROJECTIONS.
      imin1=imin+Npo
      imax1=imax+Npo
      tm1=Jm
      q12=Tm2-tm1
      delE=e1-e2
      Write( 6,5) n1,e1,n2,e2,tj1,tm1,tj2,Tm2
      Write(11,5) n1,e1,n2,e2,tj1,tm1,tj2,Tm2
5     format('====== E(',I2,') = ',F12.6, &
             ' ---> E(',I2,') = ',F12.6,'=======', &
            /'  J1 = ',F6.3,' M1 = ',F6.3,' J2 = ', &
             F6.3,' M2 = ',F6.3)
      ! DIVISION OF THE INTERVAL [imin1,imax1] INTO SHELLS. FOR
      ! EACH SHELL ONLY Q.N. Jz CHANGES. k IS INDEX OF A SHELL
      k=0
      i=imin1
500   k=k+1
        no=Nh(i)
        ind(k)=Nn(no)
        ind(k+IPx)=Ll(no)
        ind(k+2*IPx)=Jj(no)
        ! FIRST AND LAST SHELLS CAN BE SHORTER THEN 2J+1
        i1(k)=i
        i2(k)=i+(Jj(no)-Jz(i))/2
        If (i2(k) > imax1) i2(k)=imax1
        i=i2(k)+1
      If (i < imax1) goto 500
      kmax=k
      If (k > IPx) Then
        Write(*,*)' Dimension of reduced TM ',kmax,' > ',IPx
        stop
      End If
      If (lf <= 1 .and. Kl == 1) Then
        Write (11,'(10(3X,I2,3X))') (k,k=1,kmax)
        Write (11,'(10(3X,I2,3X))') (i2(k)-i1(k)+1,k=1,kmax)
        Write (11,'(10(1X,3I2,1X))') (ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      End If
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ppl = 0.d0 !###  ppl  = TRACE OF TM,
      AE1 = 0.d0 !###  AE1  = REDUCED E1 AMPLITUDE IN L GAUGE
      AE1V= 0.d0 !###  AE1V = REDUCED E1 AMPLITUDE IN V GAUGE
      A   = 0.d0 !###  A    = REDUCED Magnetic HFS AMPLITUDE
      G   = 0.d0 !###  G    = REDUCED M1 AMPLITUDE
      EDM = 0.d0 !###  EDM  = ELECTRON DIPOLE AMPLITUDE
      PNC = 0.d0 !###  PNC  = WEAK CHARGE AMPLITUDE
      AM  = 0.d0 !###  AM   = ANAPOLE MOMENT AMPLITUDE
      QM  = 0.d0 !###  QM   = Magnetic Quadrupole moment amplitude
      AE2 = 0.d0 !###  AE2  = REDUCED E2 AMPLITUDE
      AE3 = 0.d0 !###  AE3  = REDUCED E3 AMPLITUDE
      AM2 = 0.d0 !###  AM2  = REDUCED M2 AMPLITUDE
      AM3 = 0.d0 !###  AM3  = REDUCED M3 AMPLITUDE

      Do l1=1,4  !###  LOOP FOR RANKS 0, 1, 2, 3
        tl=l1-1
        icc=tj2-tm2+0.1d0
        c=Isgn(icc) * Fj3(tj2,tl,tj1,-tm2,q12,tm1)
        If (c /= 0.d0) Then
          If (l1 == 1) c=c*dsqrt(2*tj1+1)
          Do k=1,kmax
            ik1=i1(k)
            ik2=i2(k)
            nok=Nh(ik1)
            If (nok /= Nh(ik2)) goto 700
            jk=Jj(nok)
            xjk=JK/2.d0
            lk=Ll(nok)
            Do l=1,kmax
              s=0.d0
              il1=i1(l)
              il2=i2(l)
              nol=Nh(il1)
              jl=Jj(nol)
              xjl=jl/2.d0
              lll=Ll(nol)
              mc=Jz(ik1)-2
              Do ik=ik1,ik2     !### - loop over Jz
                mk=Jz(ik)
                If (mk-mc /= 2) goto 700
                mc=mk
                xmk=mk/2.d0
                Do il=il1,il2   !### - loop over Jz'
                  ml=Jz(il)
                  xml=ml/2.d0
                  x=Ro(ik-Npo,il-Npo)
                  If (x /= 0.d0) Then
                    s=s+Isgn((jl-ml)/2) &
                      *Fj3(xjl,tl,xjk,-xml,q12,xmk)*x
                  End If
                End Do
              End Do
              If (l1 == 1) s=s*dsqrt(2*xjk+1)
              tme=s/c
              If (mod(lk+lll,2) /= (1-Iprt)/2) Then
                If (dabs(tme) > 1.d-5) Then
                  Write(*,*)' RdcTM: big ME of wrong parity'
                  Write(*,*)' lk=',lk,' ll=',lll,' Rro=',tme
                  Read(*,*)
                End If
                tme=0.d0
              End If
              Rro(k,l)=tme            !### k = init, l = final
              If (tme /= 0.d0) Then
              ! Parity = - 1
                If (iprt == -1) Then
                  If (l1 == 1) Then
                    EDM=EDM+AmpEDM(tme, nol,xjl,lll, nok,xjk,lk)
                    PNC=PNC+AmpPNC(tme, nol,xjl,lll, nok,xjk,lk)
                  End If
                  If (l1 == 2) Then
                    AE1=AE1+AmpE1(tme, nol,xjl,lll, nok,xjk,lk)
                    AM=AM+AmpAM(tme, nol,xjl,lll, nok,xjk,lk)
                  End If
                  If (l1 == 3) Then
                    AM2 = AM2 + AmpM2(tme, nol,xjl,lll, nok,xjk,lk)
                    QM = QM + AmpMQM(tme, nol,xjl,lll, nok,xjk,lk)
                  End If
                  If (l1 == 4) Then
                    If (iabs(jl-jk) <= 6.and.(jl+jk) >= 6) Then
                      AE3 = AE3 + AmpE3(tme, nol,xjl,lll, nok,xjk,lk)
                    End If
                  End If
                End If
                ! Parity = + 1
                If (iprt == +1) Then
                  If (l1 == 2) Then
                    A=A+HfsA(tme, nol,xjl,lll, nok,xjk,lk)
                    If (K_M1 /= 1) Then           !%%% relativistic M1
                      g1=AmpM1(tme, nol,xjl,lll, nok,xjk,lk)
                    else                          !%%% non-relativistic M1
                      g1=0.d0
                      If (lk == lll) Then
                        If (nol == nok .or. jl /= jk) Then
                          g1=tme*Fint(9,nol,nok,+1)
                          If (jl == jk) g1=Gj1(xjk,lk)*g1
                          If (jl /= jk) g1=Gjj(jl,lll)*g1
                        End If
                      End If
                    End If
                    G=G+g1
                  End If
                  If (l1 == 3 .and. iabs(jl-jk) <= 4 .and. (jl+jk) >= 4) Then
                    !print*,AE2,AmpE2(tme, nol,xjl,lll, nok,xjk,lk), tme, nol,xjl,lll, nok,xjk,lk
                    AE2= AE2 + AmpE2(tme, nol,xjl,lll, nok,xjk,lk)
                  End If
                  If (l1 == 4) Then
                    If (iabs(jl-jk) <= 6.and.(jl+jk) >= 6) Then
                      !print*,AM3,AmpM3(tme, nol,xjl,lll, nok,xjk,lk), tme,nol,xjl,lll,nok,xjk,lk
                      AM3 = AM3 + AmpM3(tme, nol,xjl,lll, nok,xjk,lk)
                    End If
                  End If
                End If

              End If
            End Do
            If (l1 == 1) ppl=ppl+Rro(k,k)
          End Do
          If (Kl == 1 .and. kmax <= 10) Then
            Write(11,'(3X,15("-")," RANK ",F3.0,1X,15("-"))') tl
            Do k=1,kmax
              Write(11,'(10F8.5)') (Rro(k,l),l=1,kmax)
            End Do
          End If
        End If
      End Do

      ! ===           AMPLITUDES IN CONVENTIONAL Unlvs            ===
      ! ===  all phys. constants below are taken from "phys.par"  ===
      If (iprt == 1) Then
        A=-Gnuc/(DPcl*4*DPmp)*A*DPau*1.d-6
        B= Qnuc*B*DPau/(DPrb*DPrb)*1.d-30
        !print*, ppl, G, A, AE2, AM3
        Write( 6,55) ppl,G,A,AE2,AM3
        Write(11,55) ppl,G,A,AE2,AM3
55      format(' Trace =',F8.4,' <b||M1||a> =',F9.5,' mu_0',2x, &
               ' <b||H_hfs||a> =',E11.4,' MHz'/,16x, &
               ' <b||E2||a> =',E13.5,' a.u.'/,16x, &
               ' <b||M3||a> =',E13.5,' mu_0')
        If (dabs(tj1-tj2) < 1.d-6) Then
          xg = dsqrt(tj1*(tj1+1)*(2*tj1+1)+1.d-77)
          Write( 6,'(22X,"G_eff =",F10.5,14X,"A_eff =", &
                E11.4," MHz")') G/xg,A/xg
          Write(11,'(22X,"G_eff =",F10.5,14X,"A_eff =", &
                E11.4," MHz")') G/xg,A/xg
        End If
      else
        If (iprt /= -1) Then
          Write(*,*) 'RdcTM: wrong parity:',iprt
          Return
        End If
        AE1V= - AE1V/(delE+1.d-77)
        EDM1 = EDM * DPau/DPrb
        QM1  = QM  * DPau/(DPrb*DPrb)
        AM1  = AM  * DPcw*DPau*dsqrt(1.5d0)/4.d0
        Wc   = ((Z-Anuc) + Z*(1-4*DPsw))
        PNC1 = PNC * DPcw*DPau/16.d0 * Wc
        Write( 6,75) AE1,AE1V,EDM,EDM1,QM,QM1,AM,AM1,PNC,PNC1,Wc, &
                     AE3,AM2
        Write(11,75) AE1,AE1V,EDM,EDM1,QM,QM1,AM,AM1,PNC,PNC1,Wc, &
                     AE3,AM2
75      format(' AE1   L ',E12.5,' V ',E12.5,' A.U.',7X, &
          '(Reduced ME)', &
          /' EDM   = ',E12.5,' = ',E12.5,' Hz/e/cm ', &
          /' MQM   = ',E12.5,' = ',E12.5,' Hz/e/cm**2 (Reduced ME)', &
          /' AM    = ',E12.5,' = ',E12.5,' Hz',9X,'(Reduced ME)', &
          /' PNC   = ',E12.5,' = ',E12.5,' Hz',9X,'(Q_w=',F7.2,')' &
          /' <b||E3||a> =',E13.5,' a.u.' &
          /' <b||M2||a> =',E13.5,' mu_0'/)
      End If
      Return
700   Write(*,*)' RdcTM: WRONG DEFINITION OF THE SHELL',k
      Write(*,*) 'IND-S=',ik1,ik2,' SHELLS=',nok,Nh(ik2)
      stop
    End Subroutine RdcTM

    Subroutine RdcDM (ntrm,Etrm,Tj1,imin,imax,lf)
      Use wigner
      ! FORMATION OF REDUCED DENSITY MATRICES OF RANKS 0, 1, AND 2 FROM Ro
      Implicit None
      Integer :: n, kmax, k, no, i, imax1, imin1, lf, imax, imin, mk, ik, &
                 lll, il1, il2, nk, nol, jl, ml, lk, jk, ippx, il, ik1, ik2, l, &
                 kx, nok, l1, ipp, nk0, lk0, jt, mj, ntr, ntrm
      Real(dp) :: A, B, G, ppl, tj, tm, tj1, Etrm, s, g1, xjl, xmk, &
                  xpp, dme, c, xjk, tl, x
      Integer, Dimension(3*IPx) :: ind
      Integer, Dimension(IPx) :: i1, i2
      Real(dp), Dimension(IPx) :: pp
      Integer, Dimension(IPx) :: npp, lpp
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      imin1=imin+Npo
      imax1=imax+Npo
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      tj=Tj1        !### tj and tm -total momentum and its projection
      jt=tj+tj+0.1d0
      tj=jt/2.d0
      mj=dabs(Jm+Jm)+1.d-1
      If (Jm < 0) mj=-mj
      tm=mj/2.d0
      Write( 6,'("======= E(",I2,") = ",F12.6," Jtot = ",F8.5, & 
              " Mtot = ",F8.5," =======")') ntrm,Etrm,tj,tm
      Write(11,'("======= E(",I2,") = ",F12.6," Jtot = ",F8.5, & 
              " Mtot = ",F8.5," =======")') ntrm,Etrm,tj,tm
      ! DIVISION OF THE INTERVAL [imin1,imax1] INTO SHELLS. FOR
      ! EACH SHELL ONLY Q.N. Jz CHANGES. k IS INDEX OF A SHELL
      k=0
      i=imin1
500   k=k+1
        no=Nh(i)
        ind(k)=Nn(no)
        ind(k+IPx)=Ll(no)
        ind(k+2*IPx)=Jj(no)
        ! FIRST AND LAST SHELLS CAN BE SHORTER THEN 2J+1
        i1(k)=i
        i2(k)=i+(Jj(no)-Jz(i))/2
        If (i2(k) > imax1) i2(k)=imax1
        i=i2(k)+1
      If (i < imax1) goto 500
      kmax=k
      If (kmax > IPx) Then
        Write(*,*)' Dimension of reduced DM =',kmax,' > ',IPx
        stop
      End If
      If (lf <= 1 .and. Kl == 1) Then
        Write (11,'(10(3X,I2,3X))') (k,k=1,kmax)
        Write (11,'(10(3X,I2,3X))') (i2(k)-i1(k)+1,k=1,kmax)
        Write (11,'(10(1X,3I2,1X))') (ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      End If
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ppl=0.d0        !### = NUMBER OF ELECTRONS,
      G=0.d0          !### = G-FACTOR,
      A=0.d0          !### A,B - HFS CONSTANTS
      B=0.d0
      lk0=-1          ! lk0,nk0,ipp,xpp used to calculate
      nk0=-1          !!  occupation numbers for n,l shells
      ipp=0
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Do l1=1,3       !### loop for ranks 0,1,2:
        tl=l1-1
        c=Isgn((jt-mj)/2) * Fj3(tj,tl,tj,-tm,0.d0,tm)
        If (c /= 0.d0) Then
          If (l1 == 1) c=c*dsqrt(2*tj+1)
          Do k=1,kmax
            ik1=i1(k)
            ik2=i2(k)
            nok=Nh(ik1)
            If (nok /= Nh(ik2)) Then
              Write(*,*) ' RdcDM: WRONG DEFINITION OF THE SHELL',k
              Write(*,*) 'IND-S=',ik1,ik2,' SHELLS=',nok,Nh(ik2)
              stop
            End If
            jk=Jj(nok)
            xjk=jk/2.d0
            lk=Ll(nok)
            nk=Nn(nok) 
            Do l=1,kmax
              s=0.d0
              il1=i1(l)
              il2=i2(l)
              nol=Nh(il1)
              jl=Jj(nol)
              xjl=jl/2.d0
              lll=Ll(nol)
!     - - - - - - - - - - - - - - - - - - - - - - - -
              Do ik=ik1,ik2      !###  LOOP OVER JZ
                mk=Jz(ik)
                xmk=mk/2.d0
                Do il=il1,il2    !###  LOOP OVER JZ'
                  ml=Jz(il)
                  If (ml == mk) Then
                    x=Ro(ik-Npo,il-Npo)
                    If (x /= 0.d0) Then
                      s = s + Isgn((jl-ml)/2) &
                      * Fj3(xjl,tl,xjk,-xmk,0.d0,xmk) * x
                    End If
                  End If
                End Do
              End Do
              If (l1 == 1) s=s*dsqrt(2*xjk+1)
              dme=s/c
              If (mod(lk+lll,2) /= 0) Then
                If (dabs(dme) > 1.d-5) Then
                  Write(*,*)' RdcDM: large ME of wrong parity'
                  Write(*,*)' lk=',lk,' ll=',lll,' Rro=',dme
                  Read(*,*)
                End If
                dme=0.d0
              End If
              Rro(k,l)=dme
              ! HFS CONSTANTS AND G-FACTOR:
              If (dme /= 0.d0) Then
                If (l1 == 2) Then
                  A=A+HfsA(dme, nol,xjl,lll, nok,xjk,lk)
                  If (K_M1 /= 1) Then           !%%% relativistic M1
                    g1=AmpM1(dme, nol,xjl,lll, nok,xjk,lk)
                  else                          !%%% non-relativistic M1
                    g1=0.d0
                    If (lk == lll) Then
                      If (nol == nok .or. jl /= jk) Then
                        g1=dme*Fint(9,nol,nok,+1)
                        If (jl == jk) g1=Gj1(xjk,lk)*g1
                        If (jl /= jk) g1=Gjj(jl,lll)*g1
                      End If
                    End If
                  End If
                  G=G+g1
                End If
                If (l1 == 3) B=B+HfsB(dme, nol,xjl,lll, nok,xjk,lk)
              End If
            End Do
            If (l1 == 1) Then  ! occ. num-s of (n,l) shells
              ppl=ppl+Rro(k,k)
              If (lk0 == lk .and. nk0 == nk) Then
                xpp=xpp+Rro(k,k)
              else
                ipp=ipp+1
                lk0=lk
                nk0=nk
                xpp=Rro(k,k)
              End If
              npp(ipp)=nk
              lpp(ipp)=lk
              pp(ipp)=xpp
            End If
          End Do
          If (Kl == 1 .and. kmax <= 10) Then
            Write(11,'(3X,15("-")," RANK ",F3.0,1X,15("-"))') tl
            Do k=1,kmax
              Write(11,'(10F8.5)') (Rro(k,l),l=1,kmax)
            End Do
          End If
          If (l1 == 1 .and. Kdm > 0) &         !# DM_out opens the file
             Call DM_out(ntrm,kmax,ind,i1)   !## DM0.RES and Writes Rro
        End If
      End Do
      ! ===  all phys. constants below are taken from "phys.par"  ===
      If (G /= 0.d0) G=G/dsqrt(tj*(tj+1)*(2*tj+1))
      If (A /= 0.d0) &
         A=-Gnuc/(DPcl*4*DPmp)/dsqrt(tj*(tj+1)*(2*tj+1)) &
         *A*DPau*1.d-6
      If (B /= 0.d0) &
         B=-2*Qnuc*dsqrt(tj*(2*tj-1)/((tj+1)*(2*tj+1)*(2*tj+3))) &
         *B*DPau/(DPrb*DPrb)*1.d-30
      ippx=4
      Do i=5,ipp
        If (pp(i) > 0.5d-4) ippx=i
      End Do
      Write ( 6,'(" Num. of El.=",F8.4,"; G =",F8.5, &
             "; A =",E15.8," MHz; B =",E15.8," MHz" &
             /" oc.num.(n,l)",4(i3,i2,F8.4),/5(i3,i2,F8.4))') &
             ppl,G,A,B,(npp(i),lpp(i),pp(i),i=1,ippx)
      Write (11,'(" Num. of El.=",F8.4,"; G =",F8.5, &
             "; A =",E15.8," MHz; B =",E15.8," MHz" &
             /" oc.num.(n,l)",4(i3,i2,F8.4),/5(i3,i2,F8.4))') &
             ppl,G,A,B,(npp(i),lpp(i),pp(i),i=1,ippx)
      Return
    End Subroutine RdcDM

    Real(dp) function Fint(is,nfin,nini,ic)       
      ! this function searches for radial integrals of one-electron operators
      Implicit None
      Integer :: isg, na, nb, is, nfin, nini, ic, ind, i
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        isg=1
        na=nfin
        nb=nini
        If (na > nb) Then
          na=nini
          nb=nfin
          isg=ic
        End If
        ind=is*IPx*IPx+(na-Nso)*IPx+(nb-Nso)
        Do i=1,Nint
          If (ind == Intg(i)) Then
            Fint=Rnt(i)*isg
            Return
          End If
        End Do
        Write( 6,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        Write(11,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        Fint=0.d0
      Return
    End function Fint

    Subroutine AmpOut(i,x,nl,nk,xjl,xjk,ll,lk,y)
      Implicit None
      Integer :: i, nl, nk, ll, lk, nnl, nnk, jl, jk
      Real(dp) :: x, xjl, xjk, y
        If (dabs(x) < 1.d-3) Return
        nnl=Nf0(nl)
        nnk=Nf0(nk)
        nnl=Nh(nnl+1)
        nnk=Nh(nnk+1)
        nnl=Nn(nnl)
        nnk=Nn(nnk)
        jl=2*xjl+0.1d0
        jk=2*xjk+0.1d0
        Write(11,'(1X,A4,":",F8.4,I4,A1,I2,"/2  <<",I4,A1,I2,"/2",2E14.4)') &
         Alet(i),x,nnk,let(lk+1),jk,nnl,let(ll+1),jl,y
       Return
    End Subroutine AmpOut

    Subroutine DM_out(ntrm,kmax,ind,i1)  
      ! this Subroutine opens the file DM0.RES and Writes Rro
      Implicit None
      Integer :: ntrm, kmax, k, i
      Integer, Dimension(3*IPx) :: ind
      Integer, Dimension(IPx) :: i1, no
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      If (Kdm < 1) Return
      If (Kdm == 1) Then
        Kdm=2
        Do k=1,kmax
          i=i1(k)             !# decoding the index k:
          no(k)=Nh(i)         !## i1 -position, no - orbital
        End Do
        Open(17,file='DM0.RES',status='UNKNOWN')
        Close(17,status='DELETE')
        Open(17,file='DM0.RES',status='NEW')
        Write(17,'(4X,"Reduced Density Matrices of the rank 0", &
              /4X,"Z =",F5.1,"; Nc =",I6," Nd =",I8 &
              /4X,I3," active orbitals are:", &
              /"  No  n  l   j",/(1X,4I3,"/2"))') Z,Nc,Nd,kmax, &
              (no(k),ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      End If
      Write(17,'(2X,68("="),/2X,"Level ",I2)') ntrm
      Do i=1,kmax
        Write(17,'(4X,"line ",I3,/(8F9.5))') no(i),(Rro(i,k),k=1,kmax)
      End Do
      Return
    End Subroutine DM_out

    Real(dp) function HfsA(Ro, nk,xjk,lk, nl,xjl,ll)          
      ! this function calculates magnetic quadrupole hyperfine structure constant <k|A|l>
      Implicit None
      Integer :: is, k, nk, lk, nl, ll
      Real(dp) :: Ro, xjk, xjl, A, c, xjm
      A=Ro*Fint(1,nk,nl,+1)
      If (A /= 0.d0) Then
        is = 1
        k=xjl+ll+1.51d0
        If (k /= 2*(k/2)) is=-is
        xjm=xjl
        If (xjm > xjk+0.1d0) xjm=xjk
        c=0.d0
        If (dabs(xjl-xjk) > 0.1d0) c=dsqrt((2*xjk+1)*(2*xjl+1) &
           /(xjm+1))
        If (dabs(xjl-xjk) <= 0.1d0) c=(2*xjk+1)*dsqrt((2*xjk+1) &
           /(xjk*(xjk+1)))
        A=A*is*c
        If (Kl == 1) Call AmpOut(1,Ro,nl,nk,xjl,xjk,ll,lk,A)
      End If
      HfsA=A
      Return
    End function HfsA

    Real(dp) function HfsB(Ro, nk,xjk,lk, nl,xjl,ll)      
      Use wigner    
      ! this function calculates electric quadrupole hyperfine structure constant <k|B|l>
      Implicit None
      Integer :: nk, nl, lk, ll, is, k
      Real(dp) :: Ro, B, xjk, xjl, xlk, xll
        B=Ro*Fint(2,nk,nl,+1)
        If (B /= 0.d0) Then
          is = 1
          k=xjl+0.51d0
          If (k /= 2*(k/2)) is=-is
          xlk=lk
          xll=ll
          B=B*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
             *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
             *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
          If (Kl == 1) Call AmpOut(2,Ro,nl,nk,xjl,xjk,ll,lk,B)
        End If
        HfsB=B
       Return
    End function HfsB

    Real(dp) function AmpE1(Ro, nk,xjk,lk, nl,xjl,ll) 
      Use wigner        
      ! this function calculates the amplitude of the electric dipole
      ! transition matrix element <k|E1|l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k, lx
      Real(dp) :: Ro, xjk, xjl, AL, AV, xlk, xll, c
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AL=Ro*Fint(3,nk,nl,+1)
      AV=Ro*Fint(6,nk,nl,-1)
      If (AL /= 0.d0 .or. AV /= 0.d0) Then
        is = 1
        lx=max(lk,ll)
        k=xjl+lx+1.51d0
        If (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        c=dsqrt((2*xjk+1)*(2*xjl+1)*lx) &
         *Fj6(xlk,xjk,0.5d0,xjl,xll,1.d0)
        AV=is*c*AV
        AL=is*c*AL
        If (Kl == 1) Then
          Call AmpOut(3,Ro,nl,nk,xjl,xjk,ll,lk,AL)
          Call AmpOut(6,Ro,nl,nk,xjl,xjk,ll,lk,AV)
        End If
      End If
      AmpE1=AL
      AE1V=AE1V+AV
      Return
    End function AmpE1

    Real(dp) function AmpE2(Ro, nk,xjk,lk, nl,xjl,ll)    
      Use wigner    
      ! this function calculates the amplitude of the electric quadrupole
      ! transition matrix element <k||E2||l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k
      Real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE= Ro*Fint(10,nk,nl,+1)
      If (AE /= 0.d0) Then
        is = 1
        k=xjl+0.51d0
        If (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        AE= AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
              *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
              *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
        If (Kl == 1) Call AmpOut(10,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      End If
      AmpE2=AE
      Return
    End function AmpE2

    Real(dp) function AmpE3(Ro, nk,xjk,lk, nl,xjl,ll)        
      Use wigner
      ! this function calculates the amplitude of the electric octupole
      ! transition matrix element <k||E3||l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k
      Real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(11,nk,nl,+1)
      If (AE /= 0.d0) Then
        is = 1
        k=xjl-0.51d0
        If (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        AE=AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
           *Fj3(xlk,xll,3.d0,0.d0,0.d0,0.d0) &
           *Fj6(xjk,xjl,3.d0,xll,xlk,0.5d0)
        If (Kl == 1) Call AmpOut(11,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      End If
      AmpE3=AE
      Return
    End function AmpE3

    Real(dp) function AmpM1(Ro, nk,xjk,lk, nl,xjl,ll)          
      ! this function calculates the amplitude of the magnetic dipole
      ! transition matrix element <k|M1|l>
      Implicit None
      Integer :: nk, nl, k, is, ll, lk
      Real(dp) :: A, Ro, xjl, xjk, xjm, c
        A=Ro*Fint(9,nk,nl,+1)
        If (A /= 0.d0) Then
          is = 1
          k=xjl+ll+1.51d0
          If (k /= 2*(k/2)) is=-is
          xjm=xjl
          If (xjm > xjk+0.1d0) xjm=xjk
          c=0.d0
          If (dabs(xjl-xjk) > 0.1d0) &
            c= dsqrt((2*xjk+1)*(2*xjl+1)/(xjm+1))
          If (dabs(xjl-xjk) <= 0.1d0) &
            c= (2*xjk+1)*dsqrt((2*xjk+1)/(xjk*(xjk+1)))
          A=A*is*c
          If (Kl == 1) Call AmpOut(9,Ro,nl,nk,xjl,xjk,ll,lk,A)
        End If
        AmpM1=A
       Return
    End function AmpM1

    Real(dp) function AmpM2(Ro, nk,xjk,lk, nl,xjl,ll)     
      Use wigner   
      ! this function calculates the amplitude of the magnetic quadrupole
      ! transition matrix <element k||M2||l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k
      Real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(12,nk,nl,+1)
      If (AE /= 0.d0) Then
        is = 1
        k = xjk+0.51d0
        If (k /= 2*(k/2)) is=-is
        AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                Fj3(xjk,xjl,2.d0,-0.5d0,0.5d0,0.d0)
        If (Kl == 1) Call AmpOut(12,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      End If
      AmpM2=AE
      Return
    End function AmpM2

    Real(dp) function AmpM3(Ro, nk,xjk,lk, nl,xjl,ll)   
      Use wigner     
      ! this function calculates the amplitude of the magnetic octupole
      ! transition matrix element <k||M3||l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k
      Real(dp) :: Ro, xjk, xjl, AE
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(13,nk,nl,+1)
      If (AE /= 0.d0) Then
        is = 1
        k = xjk+0.51d0
        If (k /= 2*(k/2)) is=-is
        AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                Fj3(xjk,xjl,3.d0,-0.5d0,0.5d0,0.d0)
        If (Kl == 1) Call AmpOut(13,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      End If
      AmpM3=AE
      Return
    End function AmpM3

    Real(dp) function AmpEDM(Ro, nk,xjk,lk, nl,xjl,ll)        
      ! this function calculates the amplitude of the P, T-odd interaction 
      ! of the electron electric dipole moment <k|D|l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k, lx
      Real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=-2*Ro*Fint(4,nk,nl,+1)
      If (A /= 0.d0 .and. Kl == 1) Call AmpOut(4,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpEDM=A
      Return
    End function AmpEDM

    Real(dp) function AmpPNC(Ro, nk,xjk,lk, nl,xjl,ll)        
      ! this function calculates the nuclear spin independent 
      ! parity nonconserving (PNC) amplitude <k|W|l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k, lx
      Real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=-Ro*Fint(5,nk,nl,-1)
      If (A /= 0.d0 .and. Kl == 1) Call AmpOut(5,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpPNC=A
      Return
    End function AmpPNC

    Real(dp) function AmpAM(Ro, nk,xjk,lk, nl,xjl,ll)         
      ! this function calculates the amplitude of the electron interaction
      ! with the P-odd nuclear anapole moment <k|Am|l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k, lx
      Real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=Ro*Fint(7,nk,nl,-1)
      lx=max(ll,lk)
      k=xjl+0.51d0+lx
      is=1
      If (k /= 2*(k/2)) is=-is
      A=is*A
      If (A /= 0.d0 .and. Kl == 1) Call AmpOut(7,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpAM=A
      Return
    End function AmpAM

    Real(dp) function AmpMQM(Ro, nk,xjk,lk, nl,xjl,ll)        
      Use wigner
      ! this function calculates the amplitude of the nucleus 
      ! magnetic quadrupole moment <k|MQM|l>
      Implicit None
      Integer :: nk, lk, nl, ll, is, k, lx, is1, is2
      Real(dp) :: Ro, xjk, xjl, B, tll, xlk, xll, g
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      B=Ro*Fint(8,nk,nl,+1)
      If (B /= 0.d0) Then
        B=1.5d0*B*dsqrt((2*xjk+1)*(2*xjl+1))
        xlk=lk
        xll=ll
        g=xlk
        If (g < xll) g=xll
        tll=2*xjl-xll
        is1 = 1
        k=xlk+g+0.1d0
        If (k /= 2*(k/2)) is1=-is1
        is2 = 1
        k=xjl+0.51d0
        If (k /= 2*(k/2)) is2=-is2
        B = B &
          * ( is1*dsqrt(30*g) &
          * FJ9(xll,xlk,1.d0,xjl,xjk,2.d0,0.5d0,0.5d0,1.d0) &
          + is2*dsqrt(2.d0/3.d0*(2*xlk+1)*(2*tll+1)) &
          * Fj3(xlk,2.d0,tll,0.d0,0.d0,0.d0) &
          * Fj6(xlk,xjk,0.5d0,xjl,tll,2.d0) )
        If (Kl == 1) Call AmpOut(8,Ro,nl,nk,xjl,xjk,ll,lk,B)
        End If
      AmpMQM=B
      Return
    End function AmpMQM

End Module dtm_aux