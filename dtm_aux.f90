module dtm_aux

  use dtm_variables
  implicit none

  contains
    subroutine Init_Char(Let,Alet,Blet,yes,chm1)
      implicit none
      character(len=1), dimension(6) :: Let
      character(len=4), dimension(13) :: Alet
      character(len=4), dimension(5) :: Blet
      character(len=4), dimension(2) :: yes*3, chm1

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
      return
    end subroutine Init_Char

    subroutine OpenFS(inam,nam,iform,kan,ntype)
      implicit none
      integer :: i, inam, iform, kan, ntype, imax
      character(len=1), dimension(12) :: fname, nam
      character(len=1) :: space
      character(len=12) :: fnam
      equivalence (fnam, fname)
      data space/' '/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     NTYPE = 0 - OLD; NTYPE = 1 - NEW, UNKNOWN
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      do i=1,inam
         fname(i)=nam(i)
      end do
      do i=inam+1,12
         fname(i)=space
      end do
      imax=12
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ntype /= 0) then
        if (iform == 0) open(unit=kan,file=fnam,status='unknown',err=710)
        if (iform == 1) open(unit=kan,file=fnam,status='unknown', &
              form='unformatted',err=710)
      else
        if (iform == 0) open(unit=kan,file=fnam,status='old',err=700)
        if (iform == 1) open(unit=kan,file=fnam,status='old', &
              form='unformatted',err=700)
      end if
      return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
 700  write( 6,'(/" NO FILE ",12A1)') (fname(i),i=1,imax)
      write(11,'(/" NO FILE ",12A1)') (fname(i),i=1,imax)
      stop
!     - - - - - - - - - - - - - - - - - - - - - - - - -
 710  write( 6,'(/" UNABLE TO OPEN FILE:",1X,12A1)') (fname(i),i=1,imax)
      write(11,'(/" UNABLE TO OPEN FILE:",1X,12A1)') (fname(i),i=1,imax)
      stop
    end subroutine OpenFS

    subroutine Input       ! reads file CONF.INP
      implicit none
      integer :: i, i1, i2, ic, istr, nx, ny, nz, ne0
      real(dp) :: x
      character(len=1), dimension(16) :: name
      character(len=4), dimension(2) :: yes*3
      !     - - - - - - - - - - - - - - - - - - - - - - - - -
      Trd=1.d-10
      Kdm=0
      K_M1=2
      open(unit=99,file='dtm.in')
        read (99,*) Kl1
        select case(Kl1)
        case(1)
          read (99,*) nterm1, nterm2
        case(2)
          read (99,*) nterm1, nterm2, nterm2f
        end select
      close(99)
      write (*,'(/4X,"Program DTM",/4X,"Cutoff parameter :",E8.1, &
                 /4X,"Full RES file - ",A3,/4X,"DM0.RES file - ",A3, &
                 /4X,"Do you want DM (1) OR TM (2)? ",I1)') Trd, yes(Kl+1), yes(Kdm+1), Kl1
      if ((Kl1-1)*(Kl1-2) /= 0) stop
      call Init_Char(Let, Alet, Blet, yes, chm1)
      if (Kl1 == 1) then   ! regime of Density matrix & expectation values
        call OpenFS(6,'DM.RES',0,11,1)
        Iprt=+1      !### parity of the transition
      else                 ! regime of Transition matrix & amplitudes
        call OpenFS(6,'TM.RES',0,11,1)
      end if
      !     Input from the file  'CONF.INP'
      call OpenFS(8,'CONF.INP',0,10,0)
      read (10,'(1X,16A1)') name
      if(Kl1 == 1) write( 6,'(/4X,"DTM: Density matrices for ",16A1)') name
      if(Kl1 == 1) write(11,'(/4X,"DTM: Density matrices for ",16A1)') name
      if(Kl1 == 2) write( 6,'(/4X,"DTM: Transition matrices for ",16A1)') name
      if(Kl1 == 2) write(11,'(/4X,"DTM: Transition matrices for ",16A1)') name
      read (10,'(5X,F5.1)') Z
      read (10,'(5X,F5.1)') Am
      read (10,'(5X,F5.1)')
      read (10,'(5X,F5.1)') Jm
      read (10,'(5X,I6)') Nso
      read (10,'(5X,I6)') Nc
      read (10,'(5X,I6)') Kv
      read (10,'(5X,I6)') Nlv
      read (10,'(5X,I6)') Ne
!    - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(Qnl(100000000)) ! Nso+5*Nc = upper bound
      Gnuc=1.d0                  !#
      Qnuc=1.d0                  !#  Optional
      Kout=0                     !#  parameters
      Kecp=0                     !#
      Gj  =0.d0                  !#
100   call inpstr(istr)          !#
      if (istr /= 1) goto 100    !#
      Kl=Kout                    !# Kl is used in DTM instead of Kout
!    - - - - - - - - - - - - - - - - - - - - - - - - -
      if (Nso /= 0) then
         read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso)
      end if
      i1=Nso+1
      do ic=1,Nc
         ne0=0
200      i2=i1+5
         read (10,'(6(4X,F7.4))') (Qnl(i),i=i1,i2)
         do i=i1,i2
            x=dabs(Qnl(i))+1.d-9
            if (x < 1.d-8) goto 210
            nx=10000*x
            ny=100*x
            nz=(nx-100*ny)
            ne0=ne0+nz
         end do
210      i2=i-1
         i1=i2+1
         if (ne0 < Ne) goto 200
         if (ne0 > Ne) then
            write(6,'(" INPUT: too many electrons for ic =",I4)') ic
            stop
         end if
      end do
      Nsp=i2
      close(unit=10)
      if(Am < 1.d0) then
        write(6,*) ' Give nuclear parameter A: '
        read(*,*) Anuc
      else
        Anuc=Am
      end if
      write( 6,'(4X,"Anuc=",F6.1,", Gnuc =",F10.5,", Qnuc =",F10.5)') Anuc,Gnuc,Qnuc
      write(11,'(4X,"Anuc=",F6.1,", Gnuc =",F10.5,", Qnuc =",F10.5)') Anuc,Gnuc,Qnuc
      return
    end subroutine Input
    
    subroutine inpstr(istr)
      implicit none
      character(len=1) :: txt(5)
      character(len=128) :: string
      integer     :: i, istr
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      istr=0
      string(1:5) = '     '
      txt(1) = ' '
      txt(2) = ' '
      txt(3) = ' '
      txt(4) = ' '
      txt(5) = ' '
      read(10,'(5a1,a)') (txt(i),i=1,5), string
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 200
      if (txt(2) /= 'b'   .and.   txt(2) /= 'B') goto 200
      if (txt(3) /= 'r'   .and.   txt(3) /= 'R') goto 200
      if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 200
      read (string,*) kbrt
      write( *,'(5a1,i6)') (txt(i),i=1,5),kbrt
      write(11,'(5a1,i6)') (txt(i),i=1,5),kbrt
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
200   if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 210
      if (txt(2) /= '_') goto 210
      if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 210
      if (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 210
      read (string,*) K_is
      write( *,'(5a1,i6)') (txt(i),i=1,5),K_is
      write(11,'(5a1,i6)') (txt(i),i=1,5),K_is
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
210   if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 220
      if (txt(2) /= '_') goto 220
      if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 220
      if (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 220
      read (string,*) C_is
      write( *,'(5a1,F6.3)') (txt(i),i=1,5),C_is
      write(11,'(5a1,F6.3)') (txt(i),i=1,5),C_is
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
220   if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 230
      if (txt(2) /= 'l'   .and.   txt(1) /= 'L') goto 230
      if (txt(3) /= 'o'   .and.   txt(3) /= 'O') goto 230
      if (txt(4) /= 'w'   .and.   txt(4) /= 'W') goto 230
      read (string,*) Klow
      write( *,'(5a1,i6)') (txt(i),i=1,5),Klow
      write(11,'(5a1,i6)') (txt(i),i=1,5),Klow
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
230   if (txt(1) /= ' ') goto 240
      if (txt(2) /= 'g'   .and.   txt(2) /= 'G') goto 240
      if (txt(3) /= 'j'   .and.   txt(3) /= 'J') goto 240
      if (txt(4) /= ' ') goto 240
      read (string,*) Gj
      write( *,'(5a1,F6.3)') (txt(i),i=1,5),Gj
      write(11,'(5a1,F6.3)') (txt(i),i=1,5),Gj
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
240   if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 250
      if (txt(2) /= 'o'   .and.   txt(2) /= 'O') goto 250
      if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 250
      if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 250
      read (string,*) kout
      write( *,'(5a1,i6)') (txt(i),i=1,5),kout
      write(11,'(5a1,i6)') (txt(i),i=1,5),kout
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
250   if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 260
      if (txt(2) /= 'e'   .and.   txt(2) /= 'E') goto 260
      if (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 260
      if (txt(4) /= 'p'   .and.   txt(4) /= 'P') goto 260
      read (string,*) kecp
      write( *,'(5a1,i6)') (txt(i),i=1,5),kecp
      write(11,'(5a1,i6)') (txt(i),i=1,5),kecp
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
260   if (txt(1) /= 'g'   .and.   txt(1) /= 'G') goto 270
      if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 270
      if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 270
      if (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 270
      read (string,*) gnuc
      write( *,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
      write(11,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
270   if (txt(1) /= 'q'   .and.   txt(1) /= 'Q') goto 280
      if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 280
      if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 280
      if (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 280
      read (string,*) qnuc
      write( *,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
      write(11,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
280   if (txt(2) /= 'k'   .and.   txt(2) /= 'K') goto 290
      if (txt(3) /= 'l'   .and.   txt(3) /= 'L') goto 290
      if (txt(4) /= '4') goto 290
      read (string,*) kl4
      write( *,'(5a1,i6)') (txt(i),i=1,5),kl4
      write(11,'(5a1,i6)') (txt(i),i=1,5),kl4
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
290   if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 300
      if (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 300
      if (txt(4) /= '4') goto 290
      read (string,*) nc4
      write( *,'(5a1,i6)') (txt(i),i=1,5),nc4
      write(11,'(5a1,i6)') (txt(i),i=1,5),nc4
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
300   if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 310
      if (txt(2) /= 'r'   .and.   txt(2) /= 'R') goto 310
      if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 310
      if (txt(4) /= '4') goto 310
      read (string,*) crt4
      write( *,'(5a1,F8.5)') (txt(i),i=1,5),crt4
      write(11,'(5a1,F8.5)') (txt(i),i=1,5),crt4
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
310   if (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 320
      if (txt(2) /= '_') goto 320
      if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 320
      if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 320
      read (string,*) n_it
      write( *,'(5a1,i6)') (txt(i),i=1,5),n_it
      write(11,'(5a1,i6)') (txt(i),i=1,5),n_it
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
320   if (txt(1) /= 'a'   .and.   txt(1) /= 'A') goto 330
      if (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 330
      if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 330
      if (txt(4) /= 'o'   .and.   txt(4) /= 'O') goto 330
      read (string,*) kautobas
      write( *,'(5a1,i6)') (txt(i),i=1,5),kautobas
      write(11,'(5a1,i6)') (txt(i),i=1,5),kautobas
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
330   if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 340
      if (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 340
      if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 340
      if (txt(4) /= '0') goto 340
      read (string,*) Cut0
      write( *,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
      write(11,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
340   if (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 350
      if (txt(2) /= 'c'   .and.   txt(2) /= 'C') goto 350
      if (txt(3) /= 'p'   .and.   txt(3) /= 'P') goto 350
      if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 350
      read (string,*) Ncpt
      write( *,'(5a1,i6)') (txt(i),i=1,5),Ncpt
      write(11,'(5a1,i6)') (txt(i),i=1,5),Ncpt
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
350   if (txt(2) /= ' '   .and.   txt(2) /= '-') goto 700
      if (txt(3) /= ' '   .and.   txt(3) /= '-') goto 700
      if (txt(4) /= ' '   .and.   txt(4) /= '-') goto 700
      if (txt(2) == ' '   .and.   txt(3) == ' ') backspace(10)
      istr=1
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
700   istr=2
      write( *,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
      write(11,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
      read(*,*)
      stop
    end subroutine inpstr

    subroutine Init   ! reads the head of the file CONF.DAT
      implicit none
      integer :: i, ni, n0, nmin, imax, kkj, jjj, llj, nnj, nj, klag, &
                 if, n1, n2, l, j, n, k, ic, i0, qi, nsu2
      real(dp) :: C1, C2, Z1, r1, r2, rmax, Bt, Al, d
      real(dp), dimension(IP6)  :: p, q, p1, q1 
      real(dp), dimension(4*IP6):: PQ
      character(len=1) :: let(6)
      logical :: longbasis
      integer, dimension(4*IPs) :: IQN
      real(dp), dimension(IPs)  :: Qq1
      equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
      equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), &
              (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
      data let/'s','p','d','f','g','h'/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      C1=0.01d0
      Cl=DPcl
      Mj=2*dabs(Jm)+0.01d0
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      open (12,file='CONF.DAT',status='OLD', &
           access='DIRECT',recl=2*IP6*IPmr,err=700)
      call ReadF(12,1,P,Q,2)
      call ReadF(12,3,P1,Q1,2)
      Z1  =PQ(1)
      if (dabs(Z-Z1) > 1.d-6) then
        write( 6,'(/2X,"Nuclear charge is changed"/2X,"Z1=",F12.6/2X,"Z2=",F12.6)') Z,Z1
        write(11,'(/2X,"Nuclear charge is changed"/2X,"Z1=",F12.6/2X,"Z2=",F12.6)') Z,Z1
        stop
      end if
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
      allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
      if (longbasis) then
        write( *,*) ' Using variant for long basis '
        write(11,*) ' Using variant for long basis '
        do ni=1,Ns
          Nn(ni)=Iqn(4*ni-3)
          Ll(ni)=Iqn(4*ni-2)
          Kk(ni)=Iqn(4*ni-1)
          Jj(ni)=Iqn(4*ni)
        end do
      else
        if=20
        do ni=1,Ns
           if=if+1
           Nn(ni)=Pq(if)+C1
           if=if+1
           LL(ni)=PQ(if)+C1
           if=if+3
           C2=dsign(C1,PQ(if))
           Kk(ni)=Pq(if)+C2
           if=if+1
           C2=dsign(C1,PQ(if))
           Jj(ni)=Pq(if)+C2
        end do
      end if
      do nj=1,Nsp
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
        do i=1,Ns
          if (nnj == Nn(i) .and. Kk(i) == kkj) then
            Nip(nj)=i
            exit
          else if (i == Ns) then
            write( 6,'(/2X,"no function for n =",I3," k =",I3)') nnj,kkj
            write(11,'(/2X,"no function for n =",I3," k =",I3)') nnj,kkj
            stop
          end if
        end do
      end do
      call ReadF(12,2,R,V,2)
      Nec=0
      if (Nso /= 0) then
        do ni=1,Nso
           Nec=Nec+Nq(ni)
        end do
      end if
      n0=0
      nmin=Nso+1
      if (nmin <= Nsp) then
        do ni=nmin,Nsp
           n0=n0+Nq(ni)
        end do
      end if
      Ne=n0/Nc
      Nst=0
      do ni=1,Ns
         imax=2*Jj(ni)+1
         do j=1,imax,2
            Nst=Nst+1
         end do
      end do
      write( 6,'(4X,"NE  =",I4,6X,"NEC =",I4,6X,"NST =",I7)') NE,NEC,NST
      n=0
      ic=0
      i0=0
      i=0
      nmin=Nso+1
      do ni=nmin,Nsp
        i=i+1
        n=n+Nq(ni)
        if (n < NE) cycle
        ic=ic+1
        if (n > NE) then
          write( 6,'(/2X,"Wrong number of electrons "/2X,"for ICONF =",I4/)') IC
          write(11,'(/2X,"Wrong number of electrons "/2X,"for ICONF =",I4/)') IC
          stop
        end if
        Nvc(ic)=i
        Nc0(ic)=Nso+i0
        i0=i0+i
        n=0
        i=0
      end do
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      do ic=1,Nc
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)
        do i=n1,n2
          ni=Nip(i)
          l =Ll(ni)+1
          j =Jj(ni)
          n =Nn(ni)
          qi=Nq( i)
          if (Nq(i) > j+1) then 
            write( 6,'(/2X,"Wrong number of electrons "/ &
               2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
               " (",F6.3,")")') ni,n,let(l),k,qi
            write(11,'(/2X,"Wrong number of electrons "/ &
               2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
               " (",F6.3,")")') ni,n,let(l),k,qi
           stop
          end if
        end do
      end do
      close(12)
      deallocate(Qnl)
      open(17,file='CONF.DET',status='OLD',form='UNFORMATTED')
      read(17) Nd1
      close(17)      
      allocate(Iarr(Ne,Nd1))
      if (Kl1 == 2) then   ! TM regime requires files CONF1.DET & CONF1.XIJ
                         ! in addition to CONF.DET & CONF.XIJ for DM regime
        open(17,file='CONF1.DET',status='OLD',form='UNFORMATTED')
        read(17) Nd2,nsu2                     !### Number of used orbitals
        Nsu=max(nsu2,Nsu)                     !### can differ for two
      close(17)                             !### spaces!      return
      end if
      return
700   write(*,*) ' No file CONF.DAT'
      stop
    end subroutine Init

    subroutine ReadF(kan,record,V1,V2,nrec)
      implicit none
      integer :: record, nrec, i, ii, nr1, nr2, kan
      real(dp), dimension(IP6) :: V1, V2
        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        read(kan,rec=nr1) (V1(i),i=1,ii)
        if (nrec == 2) read(kan,rec=nr2) (V2(i),i=1,ii)
      return
    end subroutine ReadF

    subroutine RintA 
      ! this subroutine takes radial integrals from file or calls Rint
      implicit none
      integer :: i, is, ns1, nso1, km1, ia, ib, ik, ix, nsu1
      integer :: ns2, nso2, km2, Nint2
      real(dp) :: z2, rn2
      real(dp) :: z1, rn1, x
      integer, dimension(IPs) :: l1,l2
      integer, allocatable, dimension(:) :: Intg2
      real (dp), allocatable, dimension(:) :: Rnt2
      integer, dimension(13) :: ki, ki2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Nsu=0              !### Nsu is used to eliminate integrals
      do i=Nso+1,Nsp
         is=Nip(i)
         if (is > Nsu) Nsu=is
      end do
      if (Nsu >= IPx+Nso) then
        write(*,*)' Parameter IPx is too small for Nsu=',Nsu
        stop
      end if
      open (13,file='DTM.INT',status='OLD',form='UNFORMATTED',err=200)
      read (13,end=200,err=200) ns1,nso1,z1,rn1,km1
      x=iabs(ns1-Ns)+iabs(nso1-Nso)+iabs(km1-K_M1) &
          + dabs(z1-Z)+dabs(rn1-Rnuc)
      if (x < 1.d-6) then
        read (13) Nint,(l1(i),i=1,ns1)
        allocate(Rnt(Nint), Intg(Nint))
        is=0
        do i=1,Ns
          is=is+iabs(Ll(i)-l1(i))
        end do
        if (is == 0) then
          read (13) (ki(i),i=1,13)
          read (13) (Rnt(i),Intg(i),i=1,Nint)
          write (6,'(/1X,"### Radial integrals from DTM.INT ("," Nint =",I5," ) ###", &
                 /(4X,A4," calculated by ",A4))') Nint,(Alet(i),Blet(ki(i)),i=1,13)
          write (11,'(/1X,"### Radial integrals from DTM.INT ("," Nint =",I5," ) ###", &
                 /(4X,A4," calculated by ",A4))') Nint,(Alet(i),Blet(ki(i)),i=1,13)
          close (13)
          nsu1=0
          do i=1,Nint
            ix=Intg(i)
            ik=ix/(IPx*IPx)
            ia=(ix-IPx*IPx*ik)/IPx+Nso
            ib=(ix-IPx*IPx*ik-IPx*ia)+Nso
            nsu1=max(nsu1,ia,ib)
          end do
          if (Nsu > nsu1) goto 200
         return
        end if
      end if
      open (13,file='DTM2.INT',status='OLD',form='UNFORMATTED',err=200)
      read (13,end=200,err=200) ns2,nso2,z2,rn2,km2
      x=iabs(ns1-Ns)+iabs(nso1-Nso)+iabs(km1-K_M1) &
          + dabs(z1-Z)+dabs(rn1-Rnuc)
      if (x < 1.d-6) then
        read (13) Nint2,(l2(i),i=1,ns2)
        allocate(Rnt2(Nint2), Intg2(Nint2))
        is=0
        do i=1,Ns
          is=is+iabs(Ll(i)-l1(i))
        end do
        if (is == 0) then
          read (13) (ki2(i),i=1,13)
          read (13) (Rnt2(i),Intg2(i),i=1,Nint2)
          close (13)
        end if
      end if
      
      return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
200   open(13,file='DTM.INT',status='UNKNOWN',form='UNFORMATTED')
      close(13,status='DELETE')
      call Rint
      open(13,file='DTM.INT',status='NEW',form='UNFORMATTED')
      write(13) Ns,Nso,Z,Rnuc,K_M1
      write(13) Nint,(Ll(i),i=1,Ns)
      write(13) (1,i=1,13)
      write(13) (Rnt(i),Intg(i),i=1,Nint)
      write(6, '(/5X,"### Radial integrals saved in DTM.INT ###")')
      write(11,'(/5X,"### Radial integrals saved in DTM.INT ###")')
      close(13)
      return
    end subroutine RintA

    subroutine calcNint
      ! this subroutine calculates the number of radial integrals for one-electron operators
      implicit none
      integer :: na, nb, ja, jb, la, lb, jab, ip, nmin, cnt
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      nmin=Nso+1
      cnt=0
      ! RADIAL INTEGRALS:
      do na=nmin,Nsu                   !### - final state
        la=Ll(na)
        ja=Jj(na)
        do nb=na,Nsu                  !### - initial state
          lb=Ll(nb)
          jb=Jj(nb)
          jab=iabs(ja-jb)
          ip=mod(iabs(la-lb),2)
          if (ip == 1) goto 220
          ! POSITIVE PARITY:
          if (jab > 2) goto 210
          ! Non-rel. M1-amplitude  (K_M1=1)
          if (K_M1 == 1 .and. la == lb) then
            if (na == nb .or. ja /= jb) then
              cnt=cnt+1
            end if
          end if
          ! M1-amplitude  (K_M1 /= 1)
          if (K_M1 /= 1) then
            cnt=cnt+1
          end if
          ! DIPOLE HFS:
          cnt=cnt+1
210       if (jab > 4 .or. ja+jb < 4) goto 211
          ! QUADRUPOLE HFS:
          cnt=cnt+1
          ! E2 AMPLITUDE:
          cnt=cnt+1
211       if (jab > 6 .or. ja+jb < 6) cycle
          ! M3 AMPLITUDE:
          cnt=cnt+1
          cycle
          ! NEGATIVE PARITY:
220       if (jab > 6) cycle
          if (ja+jb < 6) goto 230
          ! E3 AMPLITUDE:
          cnt=cnt+1
230       if (jab > 4) cycle
          if (ja+jb < 4) goto 240
          ! MQM AMPLITUDE:
          cnt=cnt+1
          ! M2 AMPLITUDE:
          cnt=cnt+1    
240       if (iabs(ja-jb) > 2) cycle
          ! E1 AMPLITUDE (L GAUGE):
          cnt=cnt+1
          ! E1 AMPLITUDE (V GAUGE):
          cnt=cnt+1
          ! ANAPOLE MOMENT
          cnt=cnt+1
          ! EDM OF THE ELECTRON:
          if (ja == jb) then
            cnt=cnt+1
            ! PNC AMPLITUDE:
            cnt=cnt+1
          end if
        end do
      end do
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      Nint=cnt
      return
    end subroutine calcNint

    subroutine Rint     
      use wigner
      ! this subroutine calculates radial integrals for one-electron operators
      implicit none
      integer :: kab, na, nb, ja, jb, la, lb, i, in, nmin, ih, k, &
                lg, lbs, is, las, ip, jab, n12
      real(dp) :: x, c1, c2, c3, s1, s2, w1, w2, w3, ga, gb, gab, dn, tab, Alfd, qe
      real(dp), dimension(IP6)  :: p, q, a, b, ro
      real(dp), dimension(10) :: rcut
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      MaxT=9     !## Max power in Taylor expansion at the origin
      if (Rnuc == 0.d0) Rnuc = 1.2D-13/0.529D-8*Anuc**0.33333d0
      !  >>>>>>>>>>>>>> ECP of the core <<<<<<<<<<<<<<<<<<
      if (Kecp /= 0) then
         open(unit=10,file='CONF.ECP',status='OLD')
         read (10,'(F7.3)') Alfd
         read (10,'(F7.3)') (rcut(i),i=1,10)
         write(11,'(5X,"ECP: Alfd =",F7.3,5X,"Rcut:" &
              /5(I2,":",F7.3))') Alfd,(i-1,rcut(i),i=1,10)
         close(unit=10)
      end if
      !  >>>>>>>>>>>>>>    end of ECP   <<<<<<<<<<<<<<<<<<
      Nint=0
      call calcNint
      dn=0.d0
      if (Kl == 1) write(11,'(4X,"===== RADIAL INTEGRALS =====")')
      open(unit=12,file='CONF.DAT',status='OLD',access='DIRECT',recl=2*IP6*IPmr)
      if (.not. allocated(Rnt)) allocate(Rnt(Nint))
      if (.not. allocated(Intg)) allocate(Intg(Nint))
      Nint=0
      call ReadF(12,2,R,V,2)
      ih=2-Kt
      nmin=Nso+1
      ! CORE ELECTRON DENSITY AND CHARGE(R):
      do i=1,Ii,ih
        C(i)=0.d0
      end do
      gab=10.d0
      do in=1,Nso
        qe=Jj(in)+1
        n12=in+4
        call ReadF(12,n12,p,q,2)
        ga=(p(Ii+4))*2                ! changed on 1/11/11
        if (gab > ga) gab=ga
        do i=1,Ii,ih
           C(i)=C(i)+((p(i))**2 +(q(i))**2)*qe
        end do
      end do
      ro(1)=Z-C(1)*R(1)/(gab+1)
      do i=2,Ii,ih
        ro(i)=ro(i-1)-(C(i-1)+C(i))/2*(R(i)-R(i-1))
      end do
      if (Kl == 1) write(11,'(12F5.1)')(ro(i),i=1,Ii,10)
      ! RADIAL INTEGRALS:
      do na=nmin,Nsu                   !### - final state
        la=Ll(na)
        ja=Jj(na)
        call ReadF (12,na+4,p,q,2)
        ga=p(Ii+4)
        do nb=na,Nsu                  !### - initial state
          lb=Ll(nb)
          jb=Jj(nb)
          call ReadF (12,nb+4,a,b,2)
          gb=a(Ii+4)
          gab=ga+gb
          kab=iabs(Kk(na))+iabs(Kk(nb))
          jab=iabs(ja-jb)
          ip=mod(iabs(la-lb),2)
          if (ip == 1) goto 220
          ! POSITIVE PARITY:
          if (jab > 2) goto 210
          ! Non-rel. M1-amplitude  (K_M1=1)
          if (K_M1 == 1 .and. la == lb) then
            if (na == nb .or. ja /= jb) then
              do i=1,Ii,ih
                C(i)=(p(i)*a(i)+q(i)*b(i))
              end do
              C(ii+4)=gab
              call Sint1(tab)
              call AddRint(9,na,nb,0.d0,tab)
            end if
          end if
          ! M1-amplitude  (K_M1 /= 1)
          if (K_M1 /= 1) then
            do i=1,Ii,ih
              C(i)=0.5d0*DPcl*(p(i)*b(i)+q(i)*a(i))*R(i)
            end do
            C(ii+4)=gab+1
            call Sint1(tab)
            call AddRint(9,na,nb,0.d0,tab)
          end if
          ! DIPOLE HFS:
          do i=1,Ii,ih
            C(i)=-(p(i)*b(i)+q(i)*a(i))/(R(i)**2)
          end do
          C(ii+4)=gab-2
          call Sint1(tab)
          if (Am >= 1.d0) then
            call NclInt(kab-1,kab+1,-1.d0,-1.d0,p,b,a,q,dn)
            tab=tab-Dint+dn
          end if
          call AddRint(1,na,nb,dn,tab)
210       if (jab > 4 .or. ja+jb < 4) goto 211
          ! QUADRUPOLE HFS:
          do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))/(R(i)**3)
          end do
          C(ii+4)=gab-3
          call Sint1(tab)
          if (Am >= 1.d0) then
            call NclInt(kab-2,kab+2,1.d0,1.d0,p,a,q,b,dn)
            tab=tab-Dint+dn
          end if
          call AddRint(2,na,nb,dn,tab)
          ! E2 AMPLITUDE:
          do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*(R(i)**2)
          end do
          C(ii+4)= gab+2
          call Sint1(tab)
          call AddRint(10,na,nb,dn,tab)
211       if (jab > 6 .or. ja+jb < 6) cycle
          ! M3 AMPLITUDE:
          do i=1,Ii,ih
            C(i)= -0.5d0*DPcl*(KK(na)+KK(nb))* &  ! (-) stands because P= f*r, Q = -g*r
                  (p(i)*b(i) + q(i)*a(i))* R(i)**3
          end do
          C(ii+4)=gab+3
          call Sint1(tab)
          call AddRint(13,na,nb,0.d0,tab)
          cycle
          ! NEGATIVE PARITY:
220       if (jab > 6) cycle
          if (ja+jb < 6) goto 230
          ! E3 AMPLITUDE:
          do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*(R(i)**3)
          end do
          C(ii+4)= gab+3
          call Sint1(tab)
          call AddRint(11,na,nb,0.d0,tab)
230       if (jab > 4) cycle
          if (ja+jb < 4) goto 240
          ! MQM AMPLITUDE:
          do i=1,Ii,ih
            C(i)=-(p(i)*b(i)+q(i)*a(i))/(R(i))**3
          end do
          C(ii+4)=gab-3
          call Sint1(tab)
          if (Am >= 1.d0) then
            call NclInt(kab-2,kab+2,-1.d0,-1.d0,p,b,a,q,dn)
            tab=tab-Dint+dn
          end if
          call AddRint(8,na,nb,dn,tab)
          ! M2 AMPLITUDE:
          do i=1,Ii,ih
            C(i)= -2/3.d0*DPcl*(KK(na)+KK(nb))*  &  ! (-) stands because P= f*r, Q = -g*r
                  (p(i)*b(i) + q(i)*a(i))* R(i)**2
          end do
          C(ii+4)=gab+2
          call Sint1(tab)
          call AddRint(12,na,nb,0.d0,tab)    
240       if (iabs(ja-jb) > 2) cycle
          ! E1 AMPLITUDE (L GAUGE):
          do i=1,Ii,ih
            C(i)= (p(i)*a(i)+q(i)*b(i))*R(i)
            if (Kecp /= 0) C(i)=C(i)*(1-Alfd &
                 /sqrt((R(i)**2+Rcut(la+1)*Rcut(lb+1))**3))
          end do
          C(ii+4)=gab+1
          call Sint1(tab)
          call AddRint(3,na,nb,Dint,tab)
          ! E1 AMPLITUDE (V GAUGE):
          las=ja-la
          lbs=jb-lb
          lg=max(la,lb)
          is=1
          k=(ja+jb)/2+la+lg
          if (k /= 2*(k/2)) is=-is
          w3 = Fj6(la+0.d0,ja/2.d0,0.5d0,jb/2.d0,lb+0.d0,1.d0)
          w1 = Fj6(0.5d0,ja/2.d0,la+0.d0,jb/2.d0,0.5d0,1.d0)/w3
          w2 = Fj6(0.5d0,ja/2.d0,lb+0.d0,jb/2.d0,0.5d0,1.d0)/w3
          do i=1,Ii,ih
            C(i)=0.d0
            if (la == lbs) C(i) = C(i) - p(i)*b(i)*w1
            if (lb == las) C(i) = C(i) - q(i)*a(i)*w2
          end do
          C(ii+4)=gab
          call Sint1(tab)
          if (Am >= 1.d0) then
            s1=-w1
            s2=-w2
            if (la /= lbs) s1=0.d0
            if (lb /= las) s2=0.d0
            call NclInt(kab+1,kab,s1,s2,p,b,q,a,dn)
            tab=tab-Dint+dn
          end if
          tab=tab*cl*is*dsqrt(6.d0/lg)
          dn=dn*cl*is*dsqrt(6.d0/lg)
          call AddRint(6,na,nb,dn,tab)
          ! ANAPOLE MOMENT
          if (gab < 2.5d0) then
            Dint=0.d0
            if (Am >= 1.d0) then
              if (la == 0) then
                s1=-3.d0
                s2=-1.d0
              else
                s1= 1.d0
                s2= 3.d0
              end if
              call NclInt(kab-2,kab,s1,s2,p,b,q,a,dn)
            else
              c3=1.d0/3.d0
              if (la == 0) x=-(p(Ii+5)*b(Ii+5)+c3*a(Ii+5)*q(Ii+5))
              if (lb == 0) x= (a(Ii+5)*q(Ii+5)+c3*p(Ii+5)*b(Ii+5))
              dn=x*Rnuc**(gab-2)
            end if
            call AddRint(7,na,nb,dn,dn)
          else
            call AddRint(7,na,nb,0.d0,0.d0)
          end if
          ! EDM OF THE ELECTRON:
          if (ja == jb) then
            do i=1,Ii,ih
              C(i)= (q(i)*b(i))/(R(i)**2)*Ro(i)
            end do
            C(ii+4)=gab-2
            call Sint1(tab)
            if (Am >= 1.d0) then
              s1=0.d0
              s2=ro(1)/R(1)
              call NclInt(kab+2,kab+1,s1,s2,p,a,q,b,dn)
              tab=tab-Dint+dn
            end if
            call AddRint(4,na,nb,dn,tab)
            ! PNC AMPLITUDE:
            if (gab < 2.5d0) then
              Dint=0.d0
              if (Am >= 1.d0) then
                s1=-3.d0
                s2=3.d0
                call NclInt(kab-2,kab,s1,s2,p,b,q,a,dn)
              else
                x=-(p(Ii+5)*b(Ii+5)-a(Ii+5)*q(Ii+5))
                dn=x*Rnuc**(gab-2)
              end if
              call AddRint(5,na,nb,dn,dn)
            else
              call AddRint(5,na,nb,0.d0,0.d0)
            end if
          end if
        end do
      end do
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      write( 6,'(4X,"Nint=",I7," Nsu=",I3)') Nint,Nsu
      write(11,'(4X,"Nint=",I7," Nsu=",I3)') Nint,Nsu
      close(12)
      return
    end subroutine Rint

    subroutine Sint1(DS)        
      ! Simpson integration over r (with weight function HH*V(I))
      implicit none
      integer :: I, IH, I0, I1, I2, I3
      real(dp) :: R1, R2, R3, T1, T2, F1, F2, F3, F21, F32, F321, &
                  C0, C1, C2, G, P1, P2, T, Q, HH, Gam
      real(dp), intent(out) :: DS
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
    end subroutine SINT1

    subroutine AddRint(ir,na,nb,dn,tab)
      ! adds radial integrals to array Rint and writes to RES-file
      implicit none
      integer :: nna, na, la, ja, nnb, lb, jb, ir, nab1, nb
      real(dp) :: dn, tab
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
      if (Kl == 1) then
        if (Nint == 1) write (11,'(/2X,68("="),/2X,"Nint",3X, &
              "Type",4X,"final",3X,"init.",4X,"Nuc.Int_0",5X, &
              "Nuc.Int",5X,"Tot.Int",/2x,68("-"))')
        write(11,'(1X,I6,3X,A4,2X,I2,A1,1X,I1,"/2",1X, &
             I2,A1,1X,I1,"/2",3E13.5)') Nint,Alet(ir), &
            nna,Let(la+1),ja,nnb,Let(lb+1),jb,Dint,dn,tab
      end if
      return
    end subroutine AddRint

    subroutine NclInt(n1,n2,s1,s2,p,q,a,b,dn)
      ! radial integration inside the nucleus: R(1)**n1*int_0,1 f(x) dx
      ! where x = r/R(1)  and  f(x) = (s1*p(x)*q(x) + s2*a(x)*b(x))*x**n2
      implicit none
      integer :: i, i1, j, j1, m2, n1, n2
      real(dp) :: s1, s2, dn
      real(dp), dimension(IP6)  :: p, q, a, b
        dn=0.d0
        do i=1,MaxT+1
           i1=Ii+4+i
           do j=1,MaxT+1
              j1=Ii+4+j
              m2=i+j+n2-1
              dn = dn + (s1*p(i1)*q(j1) + s2*a(i1)*b(j1))/m2
           end do
        end do
        dn=dn*R(1)**n1
       return
    end subroutine NclInt

    subroutine Dinit
        implicit none
        integer :: i0, ni, j, n, ic, i, nmin, nem, imax
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Jz(Nst),Nh(Nst),Nq0(Nsp))
        i0=0
        do ni=1,Ns
            nem=Jj(ni)+1
            Nf0(ni)=i0
            imax=2*Jj(ni)+1
            do j=1,imax,2
                i0=i0+1
                Jz(i0)=j-nem
                Nh(i0)=ni
            end do
        end do
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        if (nmin < Nsp) then
            do ni=nmin,Nsp
                i=i+1
                Nq0(ni)=n
                n=n+Nq(ni)
                if (n >= Ne) then
                    ic=ic+1
                    Nvc(ic)=i
                    Nc0(ic)=Nso+i0
                    i0=i0+i
                    n=0
                    i=0
                end if
            end do
        end if
        return
    end subroutine Dinit

    subroutine Jterm
        implicit none
        integer      :: ndj, i, j, mt, n, im, ndi, imax, iconf, ndi1, ic1
        real(dp)       :: d
        integer, allocatable, dimension(:) :: idet
        integer, dimension(6) :: nmj
        logical :: fin
    !- - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(idet(Ne),Ndc(Nc),Jtc(Nc))
        Njd=0
        Nd=0
    !   list of determinants
        do ic1=1,Nc
            ndi=0
            ndi1=0
            iconf=ic1
            imax=1
            fin=.true.
    210     call Ndet(iconf,fin,idet)
            if (.not. fin) then
               if (M >= Mj) then
                  if (M == Mj) ndi=ndi+1
                  if (M == Mj+2) ndi1=ndi1+1
                  nd1=nd+ndi
                  if (M == Mj) then
                    call Pdet(nd1,idet)
                  end if
                  im=(M-Mj)/2+1
                  if (im > imax) imax=im
               end if
               goto 210
            end if
    !  - - - - - - - - - - - - - - - - - - - - - - - - -
            Ndc(iconf)=ndi
            Jtc(iconf)=ndi
            if (Kv == 1) Jtc(iconf)=ndi-ndi1
            Nd=Nd+Ndi
            if (ndi /= 0) then
                if (imax > 60) then
                    write(*,*) ' Jterm: array nmj is too small'
                    write(*,*) ' Dimension is 60, required:',imax
                    read(*,*)
                end if
                nmj(1:imax)=0
                fin= .true.
    220         call Ndet(iconf,fin,idet)
                if (.NOT.fin) then
                    if (M >= Mj) then
                        im=(M-Mj)/2+1
                        nmj(im)=nmj(im)+1
                    end if
                    goto 220
                end if
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !    list of terms
                im=imax+1
    230         im=im-1
                if (im >= 1) then
                    n=nmj(im)
                    if (n /= 0) then
                        mt=(im-1)*2+Mj
                        !do j=1,im
                        !    nmj(j)=nmj(j)-n
                        !end do
                        nmj(1:im)=nmj(1:im)-n
                        if (njd /= 0) then
                            do j=1,njd
                                if (Jt(j) == mt) then
                                    njt(j)=njt(j)+n
                                    goto 230
                                end if
                            end do
                        end if
                        njd=njd+1
                        if (njd > IPjd) then
                            write(*,*) ' Jterm: number of Js'
                            write(*,*) ' exceed IPjd=',IPjd
                            read(*,*)
                        end if
                        Jt(njd)=mt
                        Njt(njd)=n
                    end if
                    goto 230
                end if
            end if
        end do
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nd == 0) then
            write( 6,'(/2X,"term J =",F5.1,2X, &
               "is absent in all configurations"/)') Jm
            write(11,'(/2X,"term J =",F5.1,2X, &
               "is absent in all configurations"/)') Jm
            stop
        end if
        write( 6,'(4X,"Nd   =",I9/3X,21("=")/4X,"N",3X, &
            "  mult.",7X,"J"/3X,21("-"))') Nd
        write(11,'(4X,"Nd   =",I9/3X,21("=")/4X,"N",3X, &
            "  mult.",7X,"J"/3X,21("-"))') Nd
    240 i=0
        do j=2,njd
            if (Jt(j) < Jt(j-1)) then
                n=Jt(j)
                Jt(j)=Jt(j-1)
                Jt(j-1)=n
                n=Njt(j)
                Njt(j)=Njt(j-1)
                Njt(j-1)=n
                i=1
            end if
        end do
        if (i == 1) goto 240
        do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            write( 6,'(I5,3X,I8,F9.1)') n,Njt(n),d
            write(11,'(I5,3X,I8,F9.1)') n,Njt(n),d
        end do
        write( 6,'(3X,21("="))')
        write(11,'(3X,21("="))')
    ! - - - - - - - - - - - - - - - - - - - - - - -
        i=0
        do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            ndj=Njt(n)
        end do
        deallocate(idet)
        return
    end subroutine Jterm

    subroutine Ndet(ic,fin,idet)
        implicit none
        integer :: ic, jm, jf0, j, nqj, nj, nmin, jq, i2, j0, k, i1, n, is, &
                   n3, iq, im, i0, if0, n0, i, nqi, ni, n2, n1
        integer, allocatable, dimension(:) :: idet
        logical :: fin
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        !   construction of the next det from the list for configuration IC.
        !   fin=TRUE, if no dets in the list are left.
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        M=0
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)
        ! - - - - - - - - - - - - - - - -
        if (fin) then
            do ni=n1,n2
                nqi=Nq(ni)
                if (nqi /= 0) then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    if0=Nf0(i)-n0
                    i0 =n0+1
                    im =n0+nqi
                    do iq=i0,im
                       idet(iq)=if0+iq
                    end do
                end if
            end do
            fin=.false.
            else
            n3=n1+n2
            do is=n1,n2
                ni =n3-is
                nqi=Nq(ni)
                if (nqi /= 0) then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    i0 =n0+1
                    im =n0+nqi
                    n=Jj(i)+1+Nf0(i)-im
                    i1=im+i0
                    do k=i0,im
                        iq=i1-k
                        if (idet(iq) < n+iq) then
                            idet(iq)=idet(iq)+1
                            j0=iq+1
                            if (j0 <= im) then
                                i2=idet(iq)-iq
                                do jq=j0,im
                                    idet(jq)=jq+i2
                                end do
                            end if
                            nmin=ni+1
                            if (nmin <= n2) then
                                do nj=nmin,n2
                                    nqj=Nq(nj)
                                    if (nqj /= 0) then
                                        j=Nip(nj)
                                        n0=Nq0(nj)
                                        jf0=Nf0(j)-n0
                                        j0=n0+1
                                        jm=n0+nqj
                                        do jq=j0,jm
                                            idet(jq)=jq+jf0
                                        end do
                                    end if
                                end do
                            end if
                            fin=.false.
                            do iq=1,Ne
                               i=idet(iq)
                               m=m+Jz(i)
                            end do
                            return
                        end if
                    end do
                end if
            end do
            fin=.true.
            return
        end if
        !  - - - - - - - - - - - - - - - - - - - - - - - - -
        ! selection of dets with particular MJ
        do iq=1,Ne
           i=idet(iq)
           m=m+Jz(i)
        end do
        return
    end subroutine Ndet

    subroutine Pdet(n,idet)
        implicit none
        integer  :: n1, n, i
        integer, allocatable, dimension(:) :: idet
        !  - - - - - - - - - - - - - - - - - - - - - - - - -
        Iarr(1:Ne,n)=idet(1:Ne)    
        return
    end subroutine Pdet

    subroutine Gdet(n,idet)
        implicit none
        integer :: i, n
        integer, allocatable, dimension(:)  :: idet
        ! - - - - - - - - - - - - - - - - - - - - - - - -
        idet(1:Ne)=Iarr(1:Ne,n)
        return
    end subroutine Gdet

    subroutine CompC(idet1,idet2,icomp)
        implicit none
        integer  :: i, i1, i2, j1, j2, icomp, is, m
        integer, allocatable, dimension(:)  :: idet1, idet2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        iconf1(1:Ne)=Nh(idet1(1:Ne))   !### ic1(i) = No of the orbital occupied
        iconf2(1:Ne)=Nh(idet2(1:Ne))    !#### by the electron i
        call Rspq(iconf1,iconf2,is,icomp,i1,j1,i2,j2)
        return
    end subroutine CompC

    subroutine Rspq(id1,id2,is,nf,i1,j1,i2,j2)
        ! this subroutine compares determinants and counts number of differences in orbitals
        implicit none
        integer  :: n, is, ni, nj, i2, j2, nf, i, j, l0, l1, l2, nn0, nn1, &
                    ll0, ll1, jj0, jj1, ndi, k, iconf, i1, j1
        integer, allocatable, dimension(:)   :: id1, id2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        is=1
        ni=0
        nj=0
        i2=0
        j2=0
        nf=3
        i=1
        j=1

        do while (i <= Ne .and. j <= Ne)
            l1=id1(i) ! id1(i) is the i-th element of determinant 1
            l2=id2(j) ! id1(j) is the j-th element of determinant 2
            if (l1 == l2) then
                i=i+1
                j=j+1
          else if (l1 > l2) then
                nj=nj+1
                if (nj > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
                j1=j2
                j2=j
                j=j+1
            else
                ni=ni+1
                if (ni > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
                i1=i2
                i2=i
                i=i+1
          end if
        end do

        if (i > j) then
            nf=ni
            do k=j,Ne
                j1=j2
                j2=k
            end do
        else
            nf=nj
            do k=i,Ne
                i1=i2
                i2=k
            end do
        end if

        select case(nf)
            case(1) ! number of differences = 1
                k=iabs(j2-i2)
                if (k /= 2*(k/2)) is=-is
                i2=id1(i2)
                j2=id2(j2)
            case(2) ! number of differences = 2
                k=iabs(j2-i2)
                if (k /= 2*(k/2)) is=-is
                k=iabs(j1-i1)
                if (k /= 2*(k/2)) is=-is
                i2=id1(i2)
                j2=id2(j2)
                i1=id1(i1)
                j1=id2(j1)
        end select

        return
    end subroutine Rspq

    subroutine BcastDMArrays(mype, npes)
      use mpi_f08
      implicit none
      integer :: mype, npes, mpierr, i
      call MPI_Bcast(e2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(tj2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      if (mype/=0) allocate(Iarr(Ne,Nd1))
      do i=1,Ne
          call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      end do
      do i=1,nlvs
          call MPI_Bcast(ArrB2(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      end do
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      return
    end subroutine BcastDMArrays

    subroutine FormDM(mype,npes)           
      use mpi_f08
      use str_fmt, only : FormattedMemSize
      ! calculates density matrix and expectation values
      implicit none
      integer :: ntr, lf, n, k, i, iq, j, ju, iu, nf, is, k1, icomp, &
                 ic2, kx, n1, ic1, imin, nx, Ndpt, j1, j2, i1, iab2, &
                 imax, nn, ntrm, ntrm1, ntrms, start, end, mpierr, pgs0, pgs, pct, size
      integer :: npes,mype
      integer*8 :: mem, memsum
      real(dp) :: s, bn, bk, Tj, Etrm
      Character(Len=16)     :: memStr, npesStr
      integer, allocatable, dimension(:) :: idet1, idet2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne))
      Npo=Nf0(Nso+1)    !### - last position of the core orbitals
      if (mype==0) write (*,'(1X," DM for terms .. - .. : ", &
              I1,1X,I1)') nterm1, nterm2
      ntrm=nterm1
      ntrm1=nterm2
      nlvs=ntrm1-ntrm+1
      allocate(B1(Nd),ArrB2(Nd,nlvs),e2s(nlvs),tj2s(nlvs))
      mem = sizeof(Iarr)+sizeof(B1)+sizeof(ArrB2)
      ! Sum all the mem sizes to get a total...
      call MPI_AllReduce(mem, memsum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
      ! ...before overwriting mem with the maximum value across all workers:
      call MPI_AllReduce(mem, mem, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
      if (mype==0) then
        Write(npesStr,fmt='(I16)') npes
        Call FormattedMemSize(memsum, memStr)
        Write(*,'(A,A,A)') 'DM requires approximately ',Trim(memStr),' of memory'
        Call FormattedMemSize(mem, memStr)
        Write(*,'(A,A,A,A,A)') 'DM requires approximately ',Trim(memStr),' of memory per core with ',Trim(AdjustL(npesStr)),' cores'
      end if
      if (ntrm <= 0) return
      if (mype==0) then
        ! read in wavefunctions of CONF.XIJ
        call OpenFS(8,'CONF.XIJ',1,16,0)
        if (ntrm > 1) then
          do ntr=1,ntrm-1
            read (16)
          end do
        end if
        iab2=0
        do ntr=ntrm,ntrm1
          iab2=iab2+1
          read (16) Etrm,tj,ndpt,(B1(i),i=1,Nd)
          ArrB2(1:Nd,iab2) = B1(1:Nd)
          e2s(iab2)=Etrm+4.d0*Gj*tj*(tj+1.d0)
          tj2s(iab2)=tj
        end do
        close(16)
      end if

      call BcastDMArrays(mype,npes)

      if (mype==0) then
        start=1
      else
        start=mype*(Nc/npes)+1
      end if
      if (mype==npes-1) then
        end = Nc
      else
        end = (mype+1)*(Nc/npes)
      end if
      size=end-start+1
      lf=0
      iab2=0
      do ntr=ntrm,ntrm1
        iab2=iab2+1
        lf=lf+1
        Ro(:,:)=0.d0
        B1(1:Nd)=ArrB2(1:Nd,iab2)
        Etrm=e2s(iab2)
        tj=tj2s(iab2)
        s=0.d0
        do i=1,Nd
          s=s+B1(i)**2
        end do
        if (dabs(s-1.d0) > 1.d-5) goto 700
        imax=0
        imin=1000
        if (mype==0) then
            n=0
        else
            n=sum(Ndc(1:start-1))
        end if
        pgs0=size/10
        pgs=start+pgs0
        pct=0
        do ic1=start,end
          nx = Ndc(ic1)
          do n1=1,nx
            n=n+1
            Ndr=n
            bn=B1(n)
            if (dabs(bn) > Trd) then
              call Gdet(n,idet1)
              k=0
              do ic2=1,ic1
                kx=Ndc(ic2)
                if(ic1 == ic2) kx=n1
                call Gdet(k+1,idet2)
                call CompC(idet1,idet2,icomp)
                if(icomp > 1) then
                  k=k+kx
                else
                  do k1=1,kx
                    k=k+1
                    bk=B1(k)
                    if (dabs(bn*bk) > Trd) then
                      call Gdet(k,idet2)
                      call Rspq(idet1,idet2,is,nf,iu,ju,i,j)
                      if (nf == 0) then        !### DETERMINANTS ARE EQUAL
                        nn=2
                        if (n == k) nn=1
                        do iq=1,Ne
                          i1=idet1(iq)-Npo
                          imax=max(i1,imax)
                          imin=min(i1,imin)
                          Ro(i1,i1)=Ro(i1,i1) + bn*bk*is*nn
                        end do
                      end if
                      if (nf == 1) then        !### DETERMINANTS DIFFER
                        j1=i-Npo               !#### BY ONE FUNCTION
                        j2=j-Npo
                        imax=max(imax,j1,j2)
                        imin=min(imin,j1,j2)
                        Ro(j1,j2)=Ro(j1,j2)+bn*bk*is
                        Ro(j2,j1)=Ro(j1,j2)
                      end if
                    end if
                  end do
                end if
              end do
              if (imax > IP1) goto 710
            end if
          end do
            if (ic1 == pgs .and. pct < 90) then
                pct=pct+10
                write(*,'(2X,"core ",I3," is",I3,"% done")') mype,pct
                pgs=pgs+pgs0
            else if (ic1 == end) then
                write(*,'(2X,"core ",I3," has completed")') mype
            end if
          end do
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        if (mype==0) then
          call MPI_Reduce(MPI_IN_PLACE, Ro(1:IP1,1:IP1), IP1*IP1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        else
          call MPI_Reduce(Ro(1:IP1,1:IP1), Ro(1:IP1,1:IP1), IP1*IP1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        end if
        call MPI_AllReduce(imax, imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
        if (mype==0) then
          s=0.d0
          do i=imin,imax
            s=s+Ro(i,i)
          end do
          if (lf == 1 .and. Kl /= 0) write( 6,35) Ne,s,imin,imax
          if (Kl /= 0)             write(11,35) Ne,s,imin,imax
35        format(1X,'Ne = ',I2,'; Trace(Ro) = ',F12.8, &
                5X,'Imin, Imax =',2I4)
          call RdcDM(ntr,Etrm,Tj,imin,imax,lf)
        end if
      end do

      deallocate(idet1,idet2)
      return
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
700   write (*,*) 'FormDM: norma for vector ',ntr,' is ',s
      return
710   write (*,*) 'FormDM: index out of range (IP1<imax):'
      write (*,*) 'IP1=',IP1,' imax=',imax,' n=',n
      return
    end subroutine FormDM

    subroutine InitTDM(mype,npes)
      use mpi_f08
      implicit none
      integer :: mype, npes, mpierr
      !----------------------------------------------
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Kl1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nd2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Trd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(nterm1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(nterm2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(nterm2f, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      if (.not. allocated(Jz)) allocate(Jz(Nst))
      if (.not. allocated(Nh)) allocate(Nh(Nst))
      if (.not. allocated(Ndc)) allocate(Ndc(Nc))
      if (.not. allocated(Nvc)) allocate(Nvc(Nc))
      if (.not. allocated(Nc0)) allocate(Nc0(Nc))
      if (.not. allocated(Rnt)) allocate(Rnt(Nint))
      if (.not. allocated(Intg)) allocate(Intg(Nint))   
      call MPI_Bcast(Nn(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Kk(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ll(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Jj(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nf0(1:IPs), IPs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Rnt(1:Nint), Nint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Intg(1:Nint), Nint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      return
    end subroutine InitTDM

    subroutine BcastTMArrays(mype, npes)
      use mpi_f08
      implicit none
      integer :: mype, npes, mpierr, i
      call MPI_Bcast(tj1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(e2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(tj2s(1:nlvs), nlvs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      if (mype/=0) allocate(Iarr(Ne,Nd1))
      do i=1,Ne
          call MPI_Bcast(Iarr(i,1:Nd1), Nd1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Iarr2(i,1:Nd2), Nd2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      end do
      call MPI_Bcast(B1(1:Nd1), Nd1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      do i=1,nlvs
          call MPI_Bcast(ArrB2(1:Nd2,i), Nd2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      end do
      return
    end subroutine BcastTMArrays

    subroutine FormTM(mype,npes)
      use mpi_f08
      use str_fmt, only : FormattedMemSize
      ! calculates transition matrix & amplitudes
      implicit none
      integer :: lf, imax, k, i, n2, n21, n22, k1, icomp, ic, &
                iq, j, ju, iu, nf, is, kx, ks, mx, kxx, ixx, j1, j2, &
                imin, i1, n, n1, ndpt, n20, jt, iab2, start, end, pgs, pgs0, pct
      integer :: mype, npes, mpierr, size
      integer*8 :: mem, memsum
      real(dp) :: tj2, bn, bk, rc, rxx, s, ms, e1, e2
      real(dp), allocatable, dimension(:) :: ro1
      integer, allocatable, dimension(:) :: idet1, idet2
      Character(Len=16)     :: memStr, npesStr
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne))
      Npo=Nf0(Nso+1)   !### - Last position of the core orbitals
      if (mype==0) write (*,'(1X," TM from term .. to terms .. - .. : &
              ",I1,1X,I1,1X,I1)') nterm1, nterm2, nterm2f
      n1=nterm1
      n21=nterm2
      n22=nterm2f
      nlvs=n22-n21+1
      allocate(e2s(nlvs),tj2s(nlvs))
      allocate(ArrB2(Nd2,nlvs))
      allocate(B1(Nd1),B2(Nd2))
      allocate(Iarr2(Ne,Nd2))
      if (n1 <= 0) return
      
      if (mype==0) then
        ! read in wavefunctions of CONF.XIJ
        call OpenFS(8,'CONF.XIJ',1,16,0)
        do n=1,n1
          read (16) e1,tj1,ndpt,(B1(i),i=1,Nd1)
          e1=e1+4.d0*Gj*tj1*(tj1+1.d0)
        end do
        close (16)
        jt=2*tj1+0.1
        tj1=jt/2.d0
  
        ! read in determinants of CONF1.DET
        call OpenFS(9,'CONF1.DET',1,17,0)
        read (17) Nd2
        write (*,*) 'Ne, Nd1, Nd2',Ne,Nd1,Nd2
        do n=1,Nd2
          read(17) (idet2(i),i=1,Ne)
          Iarr2(1:Ne,n) = idet2(1:Ne)
        end do
        close(17)
  
        ! read in wavefunctions of CONF1.XIJ
        call OpenFS(9,'CONF1.XIJ',1,16,0)
        n20=n21-1
        if (n20 /= 0) then
          do n2=1,n20
            read (16)
          end do
        end if
        iab2=0
        do n2=n21,n22
          iab2=iab2+1
          read (16) e2,tj2,ndpt,(B2(i),i=1,Nd2)
          ArrB2(1:Nd2,iab2) = B2(1:Nd2)
          e2s(iab2)=e2+4.d0*Gj*tj2*(tj2+1.d0)
          jt=2*tj2+0.1
          tj2s(iab2)=jt/2.d0
        end do
        close(16)
      end if

      mem = sizeof(Iarr)+sizeof(Iarr2)+sizeof(B1)+sizeof(B2)+sizeof(ArrB2)
      ! Sum all the mem sizes to get a total...
      call MPI_AllReduce(mem, memsum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
      ! ...before overwriting mem with the maximum value across all workers:
      call MPI_AllReduce(mem, mem, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
      if (mype==0) then
        Write(npesStr,fmt='(I16)') npes
        Call FormattedMemSize(memsum, memStr)
        Write(*,'(A,A,A)') 'TM requires approximately ',Trim(memStr),' of memory'
        Call FormattedMemSize(mem, memStr)
        Write(*,'(A,A,A,A,A)') 'TM requires approximately ',Trim(memStr),' of memory per core with ',Trim(AdjustL(npesStr)),' cores'
      end if
      
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call BcastTMArrays(mype, npes)

      if (mype==0) then
        start=1
      else
        start=mype*(Nd2/npes)+1
      end if
      if (mype==npes-1) then
        end = Nd2
      else
        end = (mype+1)*(Nd2/npes)
      end if
      size=end-start+1
      lf=0
      iab2=0
      do n2=n21,n22
        lf=lf+1
        iab2=iab2+1
        B2(1:Nd2)=ArrB2(1:Nd2,iab2)
        e2=e2s(iab2)
        tj2=tj2s(iab2)
        if (dabs(tj2-tj1) < 3.1d0) then
          Ro=0.d0
          imax=0
          imin=1000
          pgs0=size/10
          pgs=start+pgs0
          pct=0
          do n=start,end   !### - final state
            Ndr=n
            bn=B2(n)
            idet2(1:Ne)=Iarr2(1:Ne,n)
            if (n == 1) then
              ms=0     !### - Jz for the final state
              do i=1,Ne
                ks=idet2(i)
                ms=ms+Jz(ks)
              end do
              Tm2=ms/2.d0
            end if
            if (dabs(bn) > Trd) then
              k=0
              do ic=1,Nc
                kx=Ndc(ic)
                call Gdet(k+1,idet1)
                call CompC(idet2,idet1,icomp)
                if (icomp > 1) then
                  k=k+kx
                else
                  do k1=1,kx
                    k=k+1          !### - init. state
                    bk=B1(k)
                    if (dabs(bk*bn) > Trd) then
                      call Gdet(k,idet1)
                      call Rspq(idet2,idet1,is,nf,iu,ju,i,j)
                      if (nf == 0) then
                        do iq=1,Ne
                          i1=idet1(iq)-Npo
                          imax=max(imax,i1)
                          imin=min(imin,i1)
                          Ro(i1,i1)=Ro(i1,i1) + bn*bk*is
                        end do
                      end if
                      if (nf == 1) then
                        j1=j-Npo          !### - init. state
                        j2=i-Npo          !### - final state
                        imax=max(imax,j1,j2)
                        imin=min(imin,j1,j2)
                        Ro(j1,j2)=Ro(j1,j2)+bn*bk*is
                      end if
                      if (imax > IP1) then
                        write(*,*) imax,'= imax > IP1 =',IP1
                        stop
                      end if
                    end if
                  end do
                end if
              end do
            end if
            if (n == pgs .and. pct < 90) then
                pct=pct+10
                write(*,'(2X,"core ",I3," is",I3,"% done")') mype,pct
                pgs=pgs+pgs0
            else if (n == end) then
                write(*,'(2X,"core ",I3," has completed")') mype
            end if
          end do
          call MPI_Barrier(MPI_COMM_WORLD, mpierr)
          call MPI_AllReduce(imax, imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)

          allocate(ro1(imax))
          do i = 1,npes-1
            do j=1,imax
              if (mype == i) call MPI_SEND(Ro(j,1:imax),imax,MPI_DOUBLE_PRECISION,0,&
                                            0,MPI_COMM_WORLD,mpierr)
              if (mype == 0) then
                  call MPI_RECV(ro1(1:imax),imax,MPI_DOUBLE_PRECISION,i,&
                                0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
                  Ro(j,1:imax)=Ro(j,1:imax)+ro1(1:imax)
              end if
            end do
          end do
          deallocate(ro1)
          !if (mype==0) then
          !  call MPI_Reduce(MPI_IN_PLACE, Ro(1:imax,1:imax), imax**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
          !                      MPI_COMM_WORLD, mpierr)
          !else
          !  call MPI_Reduce(Ro(1:imax,1:imax), Ro(1:imax,1:imax), imax**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
          !                      MPI_COMM_WORLD, mpierr)
          !end if
          !if (mype == 0 ) then
          !  do i=1,imax
          !    print*,i,Ro(1,i)
          !  end do
          !end if
          
          if (mype==0) then 
            s=0.d0
            rxx=0.d0
            do i=imin,imax
              s=s+Ro(i,i)
              do k=imin,imax
                rc=dabs(Ro(i,k))
                if (rc > rxx) then
                  ixx=i
                  kxx=k
                  rxx=rc
                end if
              end do
            end do
            ixx=Nh(ixx+Npo)
            kxx=Nh(kxx+Npo)
            Iprt=1-2*mod(Ll(ixx)+Ll(kxx),2) !### - parity of the transition
            if (Kl /= 0) write( 6,'(/1X,"Ne = ",I2,"; Trace(Ro) = ", &
                          F12.8,5X," PARITY = ",I2)') Ne,s,Iprt
            if (Kl /= 0) write(11,'(/1X,"Ne = ",I2,"; Trace(Ro) = ", &
                          F12.8,5X," PARITY = ",I2)') Ne,s,Iprt
            call RdcTM(n1,n2,e1,e2,tj1,tj2,imin,imax,lf)
          end if
        end if
      end do
      ! - - - - - - - - - - - - - - - - - - - - - - - -
      deallocate(idet1,idet2,iconf1,iconf2)
      return
    end subroutine FormTM

    integer function Isig(n)
      implicit none
      integer :: n
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (n == 2*(n/2)) then
           Isig=1
        else
           Isig=-1
        end if
       return
    end function Isig

    integer function Isgn(n)
      implicit none
      integer :: n
      Isgn = 1 + 2*(2*(n/2)-n) != (-1)**n
      return
    end function Isgn

    real(dp) function Gj1(x,l)
      ! used for G-factors
      implicit none
      integer :: l
      real(dp) :: x
      Gj1 = (x+0.5d0)*dsqrt((2*x+1)*(4*x-2*l+1)/(2*l+1)) 
      return
    end function Gj1

    real(dp) function Gjj(kx,l)
      ! used for G-factors
      implicit none
      integer :: l, kx
      Gjj = Isgn((kx+1)/2+l)*dsqrt(2.d0*l*(l+1)/(2*l+1))
      return
    end function Gjj

    subroutine RdcTM (n1,n2,e1,e2,tj1,tj2,imin,imax,lf)
      use wigner
      implicit none
      ! FORMATION OF REDUCED TRANSITION MATRICES OF RANKS 0, 1, AND 2 FROM Ro
      integer :: k, kmax, no, i, imin1, imax1, lf, n1, n2, &
                 lk, jk, ik2, ik1, imin, imax, kx, n, l, ml, il, &
                 mk, ik, mc, lll, l1, icc, nok, il1, il2, nol, jl
      real(dp) :: ppl, delE, q12, tm1, tj2, tj1, e1, e2, &
                 s, xjk, c, tl, AM3, AM2, AE3, AE2, QM, PNC, EDM, &
                 G, A, B, AE1, x, tme, xml, xmk, xjl, xg, EDM1, QM1, &
                 AM1, Wc, PNC1, g1
      integer, dimension(3*IPx) :: ind
      integer, dimension(IPx) :: i1, i2
      ! tj1,tj2,tm1,Tm2 - TOTAL MOMENTA AND THEIR PROJECTIONS.
      imin1=imin+Npo
      imax1=imax+Npo
      tm1=Jm
      q12=Tm2-tm1
      delE=e1-e2
      write( 6,5) n1,e1,n2,e2,tj1,tm1,tj2,Tm2
      write(11,5) n1,e1,n2,e2,tj1,tm1,tj2,Tm2
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
        if (i2(k) > imax1) i2(k)=imax1
        i=i2(k)+1
      if (i < imax1) goto 500
      kmax=k
      if (k > IPx) then
        write(*,*)' Dimension of reduced TM ',kmax,' > ',IPx
        stop
      end if
      if (lf <= 1 .and. Kl == 1) then
        write (11,'(10(3X,I2,3X))') (k,k=1,kmax)
        write (11,'(10(3X,I2,3X))') (i2(k)-i1(k)+1,k=1,kmax)
        write (11,'(10(1X,3I2,1X))') (ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      end if
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

      do l1=1,4  !###  LOOP FOR RANKS 0, 1, 2, 3
        tl=l1-1
        icc=tj2-tm2+0.1d0
        c=Isgn(icc) * Fj3(tj2,tl,tj1,-tm2,q12,tm1)
        if (c /= 0.d0) then
          if (l1 == 1) c=c*dsqrt(2*tj1+1)
          do k=1,kmax
            ik1=i1(k)
            ik2=i2(k)
            nok=Nh(ik1)
            if (nok /= Nh(ik2)) goto 700
            jk=Jj(nok)
            xjk=JK/2.d0
            lk=Ll(nok)
            do l=1,kmax
              s=0.d0
              il1=i1(l)
              il2=i2(l)
              nol=Nh(il1)
              jl=Jj(nol)
              xjl=jl/2.d0
              lll=Ll(nol)
              mc=Jz(ik1)-2
              do ik=ik1,ik2     !### - loop over Jz
                mk=Jz(ik)
                if (mk-mc /= 2) goto 700
                mc=mk
                xmk=mk/2.d0
                do il=il1,il2   !### - loop over Jz'
                  ml=Jz(il)
                  xml=ml/2.d0
                  x=Ro(ik-Npo,il-Npo)
                  if (x /= 0.d0) then
                    s=s+Isgn((jl-ml)/2) &
                      *Fj3(xjl,tl,xjk,-xml,q12,xmk)*x
                  end if
                end do
              end do
              if (l1 == 1) s=s*dsqrt(2*xjk+1)
              tme=s/c
              if (mod(lk+lll,2) /= (1-Iprt)/2) then
                if (dabs(tme) > 1.d-5) then
                  write(*,*)' RdcTM: big ME of wrong parity'
                  write(*,*)' lk=',lk,' ll=',lll,' Rro=',tme
                  read(*,*)
                end if
                tme=0.d0
              end if
              Rro(k,l)=tme            !### k = init, l = final
              if (tme /= 0.d0) then
              ! Parity = - 1
                if (iprt == -1) then
                  if (l1 == 1) then
                    EDM=EDM+AmpEDM(tme, nol,xjl,lll, nok,xjk,lk)
                    PNC=PNC+AmpPNC(tme, nol,xjl,lll, nok,xjk,lk)
                  end if
                  if (l1 == 2) then
                    AE1=AE1+AmpE1(tme, nol,xjl,lll, nok,xjk,lk)
                    AM=AM+AmpAM(tme, nol,xjl,lll, nok,xjk,lk)
                  end if
                  if (l1 == 3) then
                    AM2 = AM2 + AmpM2(tme, nol,xjl,lll, nok,xjk,lk)
                    QM = QM + AmpMQM(tme, nol,xjl,lll, nok,xjk,lk)
                  end if
                  if (l1 == 4) then
                    if (iabs(jl-jk) <= 6.and.(jl+jk) >= 6) then
                      AE3 = AE3 + AmpE3(tme, nol,xjl,lll, nok,xjk,lk)
                    end if
                  end if
                end if
                ! Parity = + 1
                if (iprt == +1) then
                  if (l1 == 2) then
                    A=A+HfsA(tme, nol,xjl,lll, nok,xjk,lk)
                    if (K_M1 /= 1) then           !%%% relativistic M1
                      g1=AmpM1(tme, nol,xjl,lll, nok,xjk,lk)
                    else                          !%%% non-relativistic M1
                      g1=0.d0
                      if (lk == lll) then
                        if (nol == nok .or. jl /= jk) then
                          g1=tme*Fint(9,nol,nok,+1)
                          if (jl == jk) g1=Gj1(xjk,lk)*g1
                          if (jl /= jk) g1=Gjj(jl,lll)*g1
                        end if
                      end if
                    end if
                    G=G+g1
                  end if
                  if (l1 == 3 .and. iabs(jl-jk) <= 4 .and. (jl+jk) >= 4) then
                    !print*,AE2,AmpE2(tme, nol,xjl,lll, nok,xjk,lk), tme, nol,xjl,lll, nok,xjk,lk
                    AE2= AE2 + AmpE2(tme, nol,xjl,lll, nok,xjk,lk)
                  end if
                  if (l1 == 4) then
                    if (iabs(jl-jk) <= 6.and.(jl+jk) >= 6) then
                      !print*,AM3,AmpM3(tme, nol,xjl,lll, nok,xjk,lk), tme,nol,xjl,lll,nok,xjk,lk
                      AM3 = AM3 + AmpM3(tme, nol,xjl,lll, nok,xjk,lk)
                    end if
                  end if
                end if

              end if
            end do
            if (l1 == 1) ppl=ppl+Rro(k,k)
          end do
          if (Kl == 1 .and. kmax <= 10) then
            write(11,'(3X,15("-")," RANK ",F3.0,1X,15("-"))') tl
            do k=1,kmax
              write(11,'(10F8.5)') (Rro(k,l),l=1,kmax)
            end do
          end if
        end if
      end do

      ! ===           AMPLITUDES IN CONVENTIONAL Unlvs            ===
      ! ===  all phys. constants below are taken from "phys.par"  ===
      if (iprt == 1) then
        A=-Gnuc/(DPcl*4*DPmp)*A*DPau*1.d-6
        B= Qnuc*B*DPau/(DPrb*DPrb)*1.d-30
        !print*, ppl, G, A, AE2, AM3
        write( 6,55) ppl,G,A,AE2,AM3
        write(11,55) ppl,G,A,AE2,AM3
55      format(' Trace =',F8.4,' <b||M1||a> =',F9.5,' mu_0',2x, &
               ' <b||H_hfs||a> =',E11.4,' MHz'/,16x, &
               ' <b||E2||a> =',E13.5,' a.u.'/,16x, &
               ' <b||M3||a> =',E13.5,' mu_0')
        if (dabs(tj1-tj2) < 1.d-6) then
          xg = dsqrt(tj1*(tj1+1)*(2*tj1+1)+1.d-77)
          write( 6,'(22X,"G_eff =",F10.5,14X,"A_eff =", &
                E11.4," MHz")') G/xg,A/xg
          write(11,'(22X,"G_eff =",F10.5,14X,"A_eff =", &
                E11.4," MHz")') G/xg,A/xg
        end if
      else
        if (iprt /= -1) then
          write(*,*) 'RdcTM: wrong parity:',iprt
          return
        end if
        AE1V= - AE1V/(delE+1.d-77)
        EDM1 = EDM * DPau/DPrb
        QM1  = QM  * DPau/(DPrb*DPrb)
        AM1  = AM  * DPcw*DPau*dsqrt(1.5d0)/4.d0
        Wc   = ((Z-Anuc) + Z*(1-4*DPsw))
        PNC1 = PNC * DPcw*DPau/16.d0 * Wc
        write( 6,75) AE1,AE1V,EDM,EDM1,QM,QM1,AM,AM1,PNC,PNC1,Wc, &
                     AE3,AM2
        write(11,75) AE1,AE1V,EDM,EDM1,QM,QM1,AM,AM1,PNC,PNC1,Wc, &
                     AE3,AM2
75      format(' AE1   L ',E12.5,' V ',E12.5,' A.U.',7X, &
          '(Reduced ME)', &
          /' EDM   = ',E12.5,' = ',E12.5,' Hz/e/cm ', &
          /' MQM   = ',E12.5,' = ',E12.5,' Hz/e/cm**2 (Reduced ME)', &
          /' AM    = ',E12.5,' = ',E12.5,' Hz',9X,'(Reduced ME)', &
          /' PNC   = ',E12.5,' = ',E12.5,' Hz',9X,'(Q_w=',F7.2,')' &
          /' <b||E3||a> =',E13.5,' a.u.' &
          /' <b||M2||a> =',E13.5,' mu_0'/)
      end if
      return
700   write(*,*)' RdcTM: WRONG DEFINITION OF THE SHELL',k
      write(*,*) 'IND-S=',ik1,ik2,' SHELLS=',nok,Nh(ik2)
      stop
    end subroutine RdcTM

    subroutine RdcDM (ntrm,Etrm,Tj1,imin,imax,lf)
      use wigner
      ! FORMATION OF REDUCED DENSITY MATRICES OF RANKS 0, 1, AND 2 FROM Ro
      implicit none
      integer :: n, kmax, k, no, i, imax1, imin1, lf, imax, imin, mk, ik, &
                 lll, il1, il2, nk, nol, jl, ml, lk, jk, ippx, il, ik1, ik2, l, &
                 kx, nok, l1, ipp, nk0, lk0, jt, mj, ntr, ntrm
      real(dp) :: A, B, G, ppl, tj, tm, tj1, Etrm, s, g1, xjl, xmk, &
                  xpp, dme, c, xjk, tl, x
      integer, dimension(3*IPx) :: ind
      integer, dimension(IPx) :: i1, i2
      real(dp), dimension(IPx) :: pp
      integer, dimension(IPx) :: npp, lpp
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      imin1=imin+Npo
      imax1=imax+Npo
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      tj=Tj1        !### tj and tm -total momentum and its projection
      jt=tj+tj+0.1d0
      tj=jt/2.d0
      mj=dabs(Jm+Jm)+1.d-1
      if (Jm < 0) mj=-mj
      tm=mj/2.d0
      write( 6,'("======= E(",I2,") = ",F12.6," Jtot = ",F8.5, & 
              " Mtot = ",F8.5," =======")') ntrm,Etrm,tj,tm
      write(11,'("======= E(",I2,") = ",F12.6," Jtot = ",F8.5, & 
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
        if (i2(k) > imax1) i2(k)=imax1
        i=i2(k)+1
      if (i < imax1) goto 500
      kmax=k
      if (kmax > IPx) then
        write(*,*)' Dimension of reduced DM =',kmax,' > ',IPx
        stop
      end if
      if (lf <= 1 .and. Kl == 1) then
        write (11,'(10(3X,I2,3X))') (k,k=1,kmax)
        write (11,'(10(3X,I2,3X))') (i2(k)-i1(k)+1,k=1,kmax)
        write (11,'(10(1X,3I2,1X))') (ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      end if
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ppl=0.d0        !### = NUMBER OF ELECTRONS,
      G=0.d0          !### = G-FACTOR,
      A=0.d0          !### A,B - HFS CONSTANTS
      B=0.d0
      lk0=-1          ! lk0,nk0,ipp,xpp used to calculate
      nk0=-1          !!  occupation numbers for n,l shells
      ipp=0
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      do l1=1,3       !### loop for ranks 0,1,2:
        tl=l1-1
        c=Isgn((jt-mj)/2) * Fj3(tj,tl,tj,-tm,0.d0,tm)
        if (c /= 0.d0) then
          if (l1 == 1) c=c*dsqrt(2*tj+1)
          do k=1,kmax
            ik1=i1(k)
            ik2=i2(k)
            nok=Nh(ik1)
            if (nok /= Nh(ik2)) then
              write(*,*) ' RdcDM: WRONG DEFINITION OF THE SHELL',k
              write(*,*) 'IND-S=',ik1,ik2,' SHELLS=',nok,Nh(ik2)
              stop
            end if
            jk=Jj(nok)
            xjk=jk/2.d0
            lk=Ll(nok)
            nk=Nn(nok) 
            do l=1,kmax
              s=0.d0
              il1=i1(l)
              il2=i2(l)
              nol=Nh(il1)
              jl=Jj(nol)
              xjl=jl/2.d0
              lll=Ll(nol)
!     - - - - - - - - - - - - - - - - - - - - - - - -
              do ik=ik1,ik2      !###  LOOP OVER JZ
                mk=Jz(ik)
                xmk=mk/2.d0
                do il=il1,il2    !###  LOOP OVER JZ'
                  ml=Jz(il)
                  if (ml == mk) then
                    x=Ro(ik-Npo,il-Npo)
                    if (x /= 0.d0) then
                      s = s + Isgn((jl-ml)/2) &
                      * Fj3(xjl,tl,xjk,-xmk,0.d0,xmk) * x
                    end if
                  end if
                end do
              end do
              if (l1 == 1) s=s*dsqrt(2*xjk+1)
              dme=s/c
              if (mod(lk+lll,2) /= 0) then
                if (dabs(dme) > 1.d-5) then
                  write(*,*)' RdcDM: large ME of wrong parity'
                  write(*,*)' lk=',lk,' ll=',lll,' Rro=',dme
                  read(*,*)
                end if
                dme=0.d0
              end if
              Rro(k,l)=dme
              ! HFS CONSTANTS AND G-FACTOR:
              if (dme /= 0.d0) then
                if (l1 == 2) then
                  A=A+HfsA(dme, nol,xjl,lll, nok,xjk,lk)
                  if (K_M1 /= 1) then           !%%% relativistic M1
                    g1=AmpM1(dme, nol,xjl,lll, nok,xjk,lk)
                  else                          !%%% non-relativistic M1
                    g1=0.d0
                    if (lk == lll) then
                      if (nol == nok .or. jl /= jk) then
                        g1=dme*Fint(9,nol,nok,+1)
                        if (jl == jk) g1=Gj1(xjk,lk)*g1
                        if (jl /= jk) g1=Gjj(jl,lll)*g1
                      end if
                    end if
                  end if
                  G=G+g1
                end if
                if (l1 == 3) B=B+HfsB(dme, nol,xjl,lll, nok,xjk,lk)
              end if
            end do
            if (l1 == 1) then  ! occ. num-s of (n,l) shells
              ppl=ppl+Rro(k,k)
              if (lk0 == lk .and. nk0 == nk) then
                xpp=xpp+Rro(k,k)
              else
                ipp=ipp+1
                lk0=lk
                nk0=nk
                xpp=Rro(k,k)
              end if
              npp(ipp)=nk
              lpp(ipp)=lk
              pp(ipp)=xpp
            end if
          end do
          if (Kl == 1 .and. kmax <= 10) then
            write(11,'(3X,15("-")," RANK ",F3.0,1X,15("-"))') tl
            do k=1,kmax
              write(11,'(10F8.5)') (Rro(k,l),l=1,kmax)
            end do
          end if
          if (l1 == 1 .and. Kdm > 0) &         !# DM_out opens the file
             call DM_out(ntrm,kmax,ind,i1)   !## DM0.RES and writes Rro
        end if
      end do
      ! ===  all phys. constants below are taken from "phys.par"  ===
      if (G /= 0.d0) G=G/dsqrt(tj*(tj+1)*(2*tj+1))
      if (A /= 0.d0) &
         A=-Gnuc/(DPcl*4*DPmp)/dsqrt(tj*(tj+1)*(2*tj+1)) &
         *A*DPau*1.d-6
      if (B /= 0.d0) &
         B=-2*Qnuc*dsqrt(tj*(2*tj-1)/((tj+1)*(2*tj+1)*(2*tj+3))) &
         *B*DPau/(DPrb*DPrb)*1.d-30
      ippx=4
      do i=5,ipp
        if (pp(i) > 0.5d-4) ippx=i
      end do
      write ( 6,'(" Num. of El.=",F8.4,"; G =",F8.5, &
             "; A =",E15.8," MHz; B =",E15.8," MHz" &
             /" oc.num.(n,l)",4(i3,i2,F8.4),/5(i3,i2,F8.4))') &
             ppl,G,A,B,(npp(i),lpp(i),pp(i),i=1,ippx)
      write (11,'(" Num. of El.=",F8.4,"; G =",F8.5, &
             "; A =",E15.8," MHz; B =",E15.8," MHz" &
             /" oc.num.(n,l)",4(i3,i2,F8.4),/5(i3,i2,F8.4))') &
             ppl,G,A,B,(npp(i),lpp(i),pp(i),i=1,ippx)
      return
    end subroutine RdcDM

    real(dp) function Fint(is,nfin,nini,ic)       
      ! this function searches for radial integrals of one-electron operators
      implicit none
      integer :: isg, na, nb, is, nfin, nini, ic, ind, i
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        isg=1
        na=nfin
        nb=nini
        if (na > nb) then
          na=nini
          nb=nfin
          isg=ic
        end if
        ind=is*IPx*IPx+(na-Nso)*IPx+(nb-Nso)
        do i=1,Nint
          if (ind == Intg(i)) then
            Fint=Rnt(i)*isg
            return
          end if
        end do
        write( 6,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        write(11,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        Fint=0.d0
      return
    end function Fint

    subroutine AmpOut(i,x,nl,nk,xjl,xjk,ll,lk,y)
      implicit none
      integer :: i, nl, nk, ll, lk, nnl, nnk, jl, jk
      real(dp) :: x, xjl, xjk, y
        if (dabs(x) < 1.d-3) return
        nnl=Nf0(nl)
        nnk=Nf0(nk)
        nnl=Nh(nnl+1)
        nnk=Nh(nnk+1)
        nnl=Nn(nnl)
        nnk=Nn(nnk)
        jl=2*xjl+0.1d0
        jk=2*xjk+0.1d0
        write(11,'(1X,A4,":",F8.4,I4,A1,I2,"/2  <<",I4,A1,I2,"/2",2E14.4)') &
         Alet(i),x,nnk,let(lk+1),jk,nnl,let(ll+1),jl,y
       return
    end subroutine AmpOut

    subroutine DM_out(ntrm,kmax,ind,i1)  
      ! this subroutine opens the file DM0.RES and writes Rro
      implicit none
      integer :: ntrm, kmax, k, i
      integer, dimension(3*IPx) :: ind
      integer, dimension(IPx) :: i1, no
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      if (Kdm < 1) return
      if (Kdm == 1) then
        Kdm=2
        do k=1,kmax
          i=i1(k)             !# decoding the index k:
          no(k)=Nh(i)         !## i1 -position, no - orbital
        end do
        open(17,file='DM0.RES',status='UNKNOWN')
        close(17,status='DELETE')
        open(17,file='DM0.RES',status='NEW')
        write(17,'(4X,"Reduced Density Matrices of the rank 0", &
              /4X,"Z =",F5.1,"; Nc =",I6," Nd =",I8 &
              /4X,I3," active orbitals are:", &
              /"  No  n  l   j",/(1X,4I3,"/2"))') Z,Nc,Nd,kmax, &
              (no(k),ind(k),ind(k+IPx),ind(k+2*IPx),k=1,kmax)
      end if
      write(17,'(2X,68("="),/2X,"Level ",I2)') ntrm
      do i=1,kmax
        write(17,'(4X,"line ",I3,/(8F9.5))') no(i),(Rro(i,k),k=1,kmax)
      end do
      return
    end subroutine DM_out

    real(dp) function HfsA(Ro, nk,xjk,lk, nl,xjl,ll)          
      ! this function calculates magnetic quadrupole hyperfine structure constant <k|A|l>
      implicit none
      integer :: is, k, nk, lk, nl, ll
      real(dp) :: Ro, xjk, xjl, A, c, xjm
      A=Ro*Fint(1,nk,nl,+1)
      if (A /= 0.d0) then
        is = 1
        k=xjl+ll+1.51d0
        if (k /= 2*(k/2)) is=-is
        xjm=xjl
        if (xjm > xjk+0.1d0) xjm=xjk
        c=0.d0
        if (dabs(xjl-xjk) > 0.1d0) c=dsqrt((2*xjk+1)*(2*xjl+1) &
           /(xjm+1))
        if (dabs(xjl-xjk) <= 0.1d0) c=(2*xjk+1)*dsqrt((2*xjk+1) &
           /(xjk*(xjk+1)))
        A=A*is*c
        if (Kl == 1) call AmpOut(1,Ro,nl,nk,xjl,xjk,ll,lk,A)
      end if
      HfsA=A
      return
    end function HfsA

    real(dp) function HfsB(Ro, nk,xjk,lk, nl,xjl,ll)      
      use wigner    
      ! this function calculates electric quadrupole hyperfine structure constant <k|B|l>
      implicit none
      integer :: nk, nl, lk, ll, is, k
      real(dp) :: Ro, B, xjk, xjl, xlk, xll
        B=Ro*Fint(2,nk,nl,+1)
        if (B /= 0.d0) then
          is = 1
          k=xjl+0.51d0
          if (k /= 2*(k/2)) is=-is
          xlk=lk
          xll=ll
          B=B*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
             *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
             *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
          if (Kl == 1) call AmpOut(2,Ro,nl,nk,xjl,xjk,ll,lk,B)
        end if
        HfsB=B
       return
    end function HfsB

    real(dp) function AmpE1(Ro, nk,xjk,lk, nl,xjl,ll) 
      use wigner        
      ! this function calculates the amplitude of the electric dipole
      ! transition matrix element <k|E1|l>
      implicit none
      integer :: nk, lk, nl, ll, is, k, lx
      real(dp) :: Ro, xjk, xjl, AL, AV, xlk, xll, c
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AL=Ro*Fint(3,nk,nl,+1)
      AV=Ro*Fint(6,nk,nl,-1)
      if (AL /= 0.d0 .or. AV /= 0.d0) then
        is = 1
        lx=max(lk,ll)
        k=xjl+lx+1.51d0
        if (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        c=dsqrt((2*xjk+1)*(2*xjl+1)*lx) &
         *Fj6(xlk,xjk,0.5d0,xjl,xll,1.d0)
        AV=is*c*AV
        AL=is*c*AL
        if (Kl == 1) then
          call AmpOut(3,Ro,nl,nk,xjl,xjk,ll,lk,AL)
          call AmpOut(6,Ro,nl,nk,xjl,xjk,ll,lk,AV)
        end if
      end if
      AmpE1=AL
      AE1V=AE1V+AV
      return
    end function AmpE1

    real(dp) function AmpE2(Ro, nk,xjk,lk, nl,xjl,ll)    
      use wigner    
      ! this function calculates the amplitude of the electric quadrupole
      ! transition matrix element <k||E2||l>
      implicit none
      integer :: nk, lk, nl, ll, is, k
      real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE= Ro*Fint(10,nk,nl,+1)
      if (AE /= 0.d0) Then
        is = 1
        k=xjl+0.51d0
        if (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        AE= AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
              *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
              *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
        if (Kl == 1) call AmpOut(10,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      end if
      AmpE2=AE
      return
    end function AmpE2

    real(dp) function AmpE3(Ro, nk,xjk,lk, nl,xjl,ll)        
      use wigner
      ! this function calculates the amplitude of the electric octupole
      ! transition matrix element <k||E3||l>
      implicit none
      integer :: nk, lk, nl, ll, is, k
      real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(11,nk,nl,+1)
      if (AE /= 0.d0) then
        is = 1
        k=xjl-0.51d0
        if (k /= 2*(k/2)) is=-is
        xlk=lk
        xll=ll
        AE=AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
           *Fj3(xlk,xll,3.d0,0.d0,0.d0,0.d0) &
           *Fj6(xjk,xjl,3.d0,xll,xlk,0.5d0)
        if (Kl == 1) call AmpOut(11,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      end if
      AmpE3=AE
      return
    end function AmpE3

    real(dp) function AmpM1(Ro, nk,xjk,lk, nl,xjl,ll)          
      ! this function calculates the amplitude of the magnetic dipole
      ! transition matrix element <k|M1|l>
      implicit none
      integer :: nk, nl, k, is, ll, lk
      real(dp) :: A, Ro, xjl, xjk, xjm, c
        A=Ro*Fint(9,nk,nl,+1)
        if (A /= 0.d0) then
          is = 1
          k=xjl+ll+1.51d0
          if (k /= 2*(k/2)) is=-is
          xjm=xjl
          if (xjm > xjk+0.1d0) xjm=xjk
          c=0.d0
          if (dabs(xjl-xjk) > 0.1d0) &
            c= dsqrt((2*xjk+1)*(2*xjl+1)/(xjm+1))
          if (dabs(xjl-xjk) <= 0.1d0) &
            c= (2*xjk+1)*dsqrt((2*xjk+1)/(xjk*(xjk+1)))
          A=A*is*c
          if (Kl == 1) call AmpOut(9,Ro,nl,nk,xjl,xjk,ll,lk,A)
        end if
        AmpM1=A
       return
    end function AmpM1

    real(dp) function AmpM2(Ro, nk,xjk,lk, nl,xjl,ll)     
      use wigner   
      ! this function calculates the amplitude of the magnetic quadrupole
      ! transition matrix <element k||M2||l>
      implicit none
      integer :: nk, lk, nl, ll, is, k
      real(dp) :: Ro, xjk, xjl, AE, xlk, xll
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(12,nk,nl,+1)
      if (AE /= 0.d0) Then
        is = 1
        k = xjk+0.51d0
        if (k /= 2*(k/2)) is=-is
        AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                Fj3(xjk,xjl,2.d0,-0.5d0,0.5d0,0.d0)
        if (Kl == 1) call AmpOut(12,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      end if
      AmpM2=AE
      return
    end function AmpM2

    real(dp) function AmpM3(Ro, nk,xjk,lk, nl,xjl,ll)   
      use wigner     
      ! this function calculates the amplitude of the magnetic octupole
      ! transition matrix element <k||M3||l>
      implicit none
      integer :: nk, lk, nl, ll, is, k
      real(dp) :: Ro, xjk, xjl, AE
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(13,nk,nl,+1)
      if (AE /= 0.d0) Then
        is = 1
        k = xjk+0.51d0
        if (k /= 2*(k/2)) is=-is
        AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                Fj3(xjk,xjl,3.d0,-0.5d0,0.5d0,0.d0)
        if (Kl == 1) call AmpOut(13,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      end if
      AmpM3=AE
      return
    end function AmpM3

    real(dp) function AmpEDM(Ro, nk,xjk,lk, nl,xjl,ll)        
      ! this function calculates the amplitude of the P, T-odd interaction 
      ! of the electron electric dipole moment <k|D|l>
      implicit none
      integer :: nk, lk, nl, ll, is, k, lx
      real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=-2*Ro*Fint(4,nk,nl,+1)
      if (A /= 0.d0 .and. Kl == 1) call AmpOut(4,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpEDM=A
      return
    end function AmpEDM

    real(dp) function AmpPNC(Ro, nk,xjk,lk, nl,xjl,ll)        
      ! this function calculates the nuclear spin independent 
      ! parity nonconserving (PNC) amplitude <k|W|l>
      implicit none
      integer :: nk, lk, nl, ll, is, k, lx
      real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=-Ro*Fint(5,nk,nl,-1)
      if (A /= 0.d0 .and. Kl == 1) call AmpOut(5,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpPNC=A
      return
    end function AmpPNC

    real(dp) function AmpAM(Ro, nk,xjk,lk, nl,xjl,ll)         
      ! this function calculates the amplitude of the electron interaction
      ! with the P-odd nuclear anapole moment <k|Am|l>
      implicit none
      integer :: nk, lk, nl, ll, is, k, lx
      real(dp) :: Ro, xjk, xjl, A
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      A=Ro*Fint(7,nk,nl,-1)
      lx=max(ll,lk)
      k=xjl+0.51d0+lx
      is=1
      if (k /= 2*(k/2)) is=-is
      A=is*A
      if (A /= 0.d0 .and. Kl == 1) call AmpOut(7,Ro,nl,nk,xjl,xjk,ll,lk,A)
      AmpAM=A
      return
    end function AmpAM

    real(dp) function AmpMQM(Ro, nk,xjk,lk, nl,xjl,ll)        
      use wigner
      ! this function calculates the amplitude of the nucleus 
      ! magnetic quadrupole moment <k|MQM|l>
      implicit none
      integer :: nk, lk, nl, ll, is, k, lx, is1, is2
      real(dp) :: Ro, xjk, xjl, B, tll, xlk, xll, g
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      B=Ro*Fint(8,nk,nl,+1)
      if (B /= 0.d0) then
        B=1.5d0*B*dsqrt((2*xjk+1)*(2*xjl+1))
        xlk=lk
        xll=ll
        g=xlk
        if (g < xll) g=xll
        tll=2*xjl-xll
        is1 = 1
        k=xlk+g+0.1d0
        if (k /= 2*(k/2)) is1=-is1
        is2 = 1
        k=xjl+0.51d0
        if (k /= 2*(k/2)) is2=-is2
        B = B &
          * ( is1*dsqrt(30*g) &
          * FJ9(xll,xlk,1.d0,xjl,xjk,2.d0,0.5d0,0.5d0,1.d0) &
          + is2*dsqrt(2.d0/3.d0*(2*xlk+1)*(2*tll+1)) &
          * Fj3(xlk,2.d0,tll,0.d0,0.d0,0.d0) &
          * Fj6(xlk,xjk,0.5d0,xjl,tll,2.d0) )
        if (Kl == 1) call AmpOut(8,Ro,nl,nk,xjl,xjk,ll,lk,B)
        end if
      AmpMQM=B
      return
    end function AmpMQM

end module dtm_aux