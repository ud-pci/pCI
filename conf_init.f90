module conf_init

    use conf_variables

    implicit none

    contains

    subroutine Input
        use hamiltonian_io, only : Hread
        implicit none
        integer  :: istr, i, i1, i2, ic, nx, ny, nz, ne0, ierr, n, k
        integer*8 :: i8
        real(dp) :: x, t
        character(len=1)  :: name(16), str2(2)*3, str4*5
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        data str2 /' NO','YES'/
        open(unit=99,file='c.in',status='OLD')
        read (99,*) Kl, Ksig, Kdsig
        write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
        close(99)
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        write( 6,'(4X,"Program Conf")')
        write(11,'(4X,"Program Conf")')
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!       input from the file 'CONF.INP'
        open(unit=10,file='CONF.INP',status='OLD')
        read (10,'(1X,16A1)') name
        write( *,'(4X,16A1)') name
        write(11,'(4X,16A1)') name
        read (10,'(5X,F5.1)') Z
        read (10,'(5X,F5.1)') Am
        read (10,'(5X,F5.1)') XJ_av
        read (10,'(5X,F5.1)') Jm
        read (10,'(5X,I6)') Nso
        read (10,'(5X,I6)') Nc
        read (10,'(5X,I6)') Kv
        read (10,'(5X,I6)') Nlv
        read (10,'(5X,I6)') Ne
!    - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Qnl(100000000)) ! 1B = upper bound
        IPlv=3*Nlv
        K_is=0                     !#
        Kbrt=0                     !#  Optional
        Kout=1                     !#  parameters
        Kecp=0                     !#
        C_is=0.d0                  !#
        Gj=0.d0                    !#
        n_it=20                    !#
100     call inpstr(istr)          !#
        if (istr /= 1) goto 100    !#
        if (dabs(C_is) < 1.d-6) K_is=0
        if (K_is == 0) C_is=0.d0
        if (K_is == 2.OR.K_is == 4) then
          open(unit=99,file='c.in',status='OLD')
          read (99,*) Kl, Ksig, Kdsig
          read(99,*) K_sms
          write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): ', K_sms
          close(99)
          if ((K_sms-1)*(K_sms-2)*(K_sms-3) /= 0) stop
        end if
!    - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nso /= 0) then
           read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso)
        end if
        i1=Nso+1
        do ic=1,Nc
           ne0=0
200        i2=i1+5
           read (10,'(6(4X,F7.4))') (Qnl(i),i=i1,i2)
           do i=i1,i2
              x=dabs(Qnl(i))+1.d-9
              if (x < 1.d-8) exit
              nx=10000*x
              ny=100*x
              nz=(nx-100*ny)
              ne0=ne0+nz
           end do
           i2=i-1
           i1=i2+1
           if (ne0 < Ne) goto 200
           if (ne0 > Ne) then
              write(6,'(" INPUT: too many electrons for ic =",I6)') ic
            Stop
           end if
        end do
        Nsp=i2
        close(unit=10)
!     - - - - - - -  case kl = 2  - - - - - - - - - - -
        if (Kl == 2) then
           write(*,'(1X," Ksig = (0,1,2): ",I1)') Ksig 
           if (Ksig /= 0) then
             write(*,'(1X," Energy dependence of Sigma (1-Yes,0-No)? ",I1)') Kdsig
           end if
           write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
           Kecp=0
           Kl=0
        else
           Ksig=0
        end if
!     - - - - - - -  case kv = 1,3  - - - - - - - - - -
        K_prj=0                    !# this key is fixed for kv=2,4
        if (Kv == 1.OR.kv == 3) then
           K_prj=1
           write( *,'(4X,"Selection of states with J =",F5.1)') XJ_av
           write(11,'(4X,"Selection of states with J =",F5.1)') XJ_av
        end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        !if (Kl /= 1) then ! If Kl=1, continue from a previous calculation
        !      open(unit=16,status='UNKNOWN',file='CONF.JJJ')
        !      close(unit=16,status='DELETE')
        !end if
        open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        read(16) (In(i),i=1,IPgnt)
        read(16) (Gnt(i),i=1,IPgnt)
        close(unit=16)
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
    200 if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 210
        if (txt(2) /= '_') goto 210
        if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 210
        if (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 210
        read (string,*) K_is
        write( *,'(5a1,i6)') (txt(i),i=1,5),K_is
        write(11,'(5a1,i6)') (txt(i),i=1,5),K_is
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    210 if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 220
        if (txt(2) /= '_') goto 220
        if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 220
        if (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 220
        read (string,*) C_is
        write( *,'(5a1,F6.3)') (txt(i),i=1,5),C_is
        write(11,'(5a1,F6.3)') (txt(i),i=1,5),C_is
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    220 if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 230
        if (txt(2) /= 'l'   .and.   txt(1) /= 'L') goto 230
        if (txt(3) /= 'o'   .and.   txt(3) /= 'O') goto 230
        if (txt(4) /= 'w'   .and.   txt(4) /= 'W') goto 230
        read (string,*) Klow
        write( *,'(5a1,i6)') (txt(i),i=1,5),Klow
        write(11,'(5a1,i6)') (txt(i),i=1,5),Klow
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    230 if (txt(1) /= ' ') goto 240
        if (txt(2) /= 'g'   .and.   txt(2) /= 'G') goto 240
        if (txt(3) /= 'j'   .and.   txt(3) /= 'J') goto 240
        if (txt(4) /= ' ') goto 240
        read (string,*) Gj
        write( *,'(5a1,F6.3)') (txt(i),i=1,5),Gj
        write(11,'(5a1,F6.3)') (txt(i),i=1,5),Gj
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    240 if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 250
        if (txt(2) /= 'o'   .and.   txt(2) /= 'O') goto 250
        if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 250
        if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 250
        read (string,*) kout
        write( *,'(5a1,i6)') (txt(i),i=1,5),kout
        write(11,'(5a1,i6)') (txt(i),i=1,5),kout
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    250 if (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 260
        if (txt(2) /= 'e'   .and.   txt(2) /= 'E') goto 260
        if (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 260
        if (txt(4) /= 'p'   .and.   txt(4) /= 'P') goto 260
        read (string,*) kecp
        write( *,'(5a1,i6)') (txt(i),i=1,5),kecp
        write(11,'(5a1,i6)') (txt(i),i=1,5),kecp
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    260 if (txt(1) /= 'g'   .and.   txt(1) /= 'G') goto 270
        if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 270
        if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 270
        if (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 270
        read (string,*) gnuc
        write( *,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
        write(11,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    270 if (txt(1) /= 'q'   .and.   txt(1) /= 'Q') goto 280
        if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 280
        if (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 280
        if (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 280
        read (string,*) qnuc
        write( *,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
        write(11,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    280 if (txt(2) /= 'k'   .and.   txt(2) /= 'K') goto 290
        if (txt(3) /= 'l'   .and.   txt(3) /= 'L') goto 290
        if (txt(4) /= '4') goto 290
        read (string,*) kl4
        write( *,'(5a1,i6)') (txt(i),i=1,5),kl4
        write(11,'(5a1,i6)') (txt(i),i=1,5),kl4
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    290 if (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 300
        if (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 300
        if (txt(4) /= '4') goto 290
        read (string,*) nc4
        write( *,'(5a1,i6)') (txt(i),i=1,5),nc4
        write(11,'(5a1,i6)') (txt(i),i=1,5),nc4
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    300 if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 310
        if (txt(2) /= 'r'   .and.   txt(2) /= 'R') goto 310
        if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 310
        if (txt(4) /= '4') goto 310
        read (string,*) crt4
        write( *,'(5a1,F8.5)') (txt(i),i=1,5),crt4
        write(11,'(5a1,F8.5)') (txt(i),i=1,5),crt4
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    310 if (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 320
        if (txt(2) /= '_') goto 320
        if (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 320
        if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 320
        read (string,*) n_it
        write( *,'(5a1,i6)') (txt(i),i=1,5),n_it
        write(11,'(5a1,i6)') (txt(i),i=1,5),n_it
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    320 if (txt(1) /= 'a'   .and.   txt(1) /= 'A') goto 330
        if (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 330
        if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 330
        if (txt(4) /= 'o'   .and.   txt(4) /= 'O') goto 330
        read (string,*) kautobas
        write( *,'(5a1,i6)') (txt(i),i=1,5),kautobas
        write(11,'(5a1,i6)') (txt(i),i=1,5),kautobas
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    330 if (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 340
        if (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 340
        if (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 340
        if (txt(4) /= '0') goto 340
        read (string,*) Cut0
        write( *,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
        write(11,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    340 if (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 350
        if (txt(2) /= 'c'   .and.   txt(2) /= 'C') goto 350
        if (txt(3) /= 'p'   .and.   txt(3) /= 'P') goto 350
        if (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 350
        read (string,*) Ncpt
        write( *,'(5a1,i6)') (txt(i),i=1,5),Ncpt
        write(11,'(5a1,i6)') (txt(i),i=1,5),Ncpt
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    350 if (txt(2) /= ' '   .and.   txt(2) /= '-') goto 700
        if (txt(3) /= ' '   .and.   txt(3) /= '-') goto 700
        if (txt(4) /= ' '   .and.   txt(4) /= '-') goto 700
        if (txt(2) == ' '   .and.   txt(3) == ' ') backspace(10)
        istr=1
        return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    700 istr=2
        write( *,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
        write(11,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
        read(*,*)
       stop
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine inpstr

    subroutine Init
        implicit none
        integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, if, &
                    ii, i1, n2, n1, l, nmin, jlj, i0, nlmax
        real(dp) :: d, c1, c2, z1
        real(dp), dimension(IP6)  :: p, q, p1, q1 
        real(dp), dimension(4*IP6):: pq
        integer, dimension(33)  ::  nnn ,jjj ,nqq 
        character(len=1), dimension(9) :: Let 
        character(len=1), dimension(33):: lll
        logical :: longbasis
        integer, dimension(4*IPs) :: IQN
        real(dp), dimension(IPs)  :: Qq1
        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), &
             (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        data Let/'s','p','d','f','g','h','i','k','l'/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        mj=2*dabs(Jm)+0.01d0
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        open (12,file='CONF.DAT',status='OLD', &
              access='DIRECT',recl=2*IP6*IPmr,err=700)
        read (12,rec=1) p
        read (12,rec=2) q
        read (12,rec=5) p1
        read (12,rec=6) q1
        z1 = pq(1)
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(Z-z1) > 1.d-6) then
           write( 6,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
           write(11,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
           read(*,*)
        end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        Rnuc=pq(13)
        dR_N=pq(16)
        longbasis=dabs(PQ(20)-0.98765d0) < 1.d-6
        write( 6,'(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kl,Z,Jm,Nsp,Ns,Nso,Nc
        write(11,'(4X,"Kl  =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3, &
                5X,"Nc =",I6)') Kl,Z,Jm,Nsp,Ns,Nso,Nc
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          write(11,*) ' Using variant for long basis '
          do ni=1,Ns
            Nn(ni)=IQN(4*ni-3)
            Ll(ni)=IQN(4*ni-2)
            Kk(ni)=IQN(4*ni-1)
            Jj(ni)=IQN(4*ni)
          end do
        else
          if=20
          do ni=1,Ns
            if=if+1
            Nn(ni)=pq(if)+c1
            if=if+1
            Ll(ni)=pq(if)+c1
            if=if+3
            c2=dsign(c1,pq(if))
            Kk(ni)=pq(if)+c2
            if=if+1
            c2=dsign(c1,pq(if))
            Jj(ni)=pq(if)+c2
          end do
        end if
        Nsu=0
        do nj=1,Nsp
           i=dsign(1.d0,Qnl(nj))
           d=dabs(Qnl(nj))+1.d-14
           d=10.0*d
           nnj=d
           d=10.0d0*(d-nnj)
           llj=d
           jlj=2*llj+i
           kkj=-i*((jlj+1)/2)
           d=100.0d0*(d-llj)
           Nq(nj)=d+0.1d0
           do ni=1,ns
              if (nnj == Nn(ni) .and. Kk(ni) == kkj) then
                exit
              else if (ni == ns) then
                write( 6,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                write(11,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                stop
              end if
           end do
           Nip(nj)=ni
           if (Nsu < ni) Nsu=ni
        end do
        nec=0
        if (Nso /= 0) nec=sum(Nq(1:Nso))
        do ni=1,Nsu
           imax=2*Jj(ni)+1
           do j=1,imax,2
              Nst=Nst+1
           end do
        end do
        write( 6,'(4X,"Number of actually used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
        write(11,'(4X,"Number of actually used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        do ni=nmin,Nsp
           i=i+1
           n=n+Nq(ni)
           if (n >= Ne) then
              ic=ic+1
              if (n > Ne) then
                 write( 6,'(/2X,"wrong number of electrons", &
                      /2X,"for configuration ICONF =",I6)') ic
                 write(11,'(/2X,"wrong number of electrons", &
                      /2X,"for configuration ICONF =",I6)') ic
                stop
              end if
              Nvc(ic)=i
              Nc0(ic)=Nso+i0
              i0=i0+i
              n=0
              i=0
           end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        write( 6,'(1X,71("="))')
        write(11,'(1X,71("="))')
        do ni=1,Nso
           l =Ll(ni)+1
           lll(ni)=let(l)
        end do
        write(11,'(1X,"Core:", 6(I2,A1,"(",I1,"/2)",I2,";"),  &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                        /6X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                        (Nn(i),lll(i),Jj(i),Nq(i),i=1,Nso)
        write(11,'(1X,71("="))')
        do ic=1,Nc
           n1=Nc0(ic)+1
           n2=Nc0(ic)+Nvc(ic)
           do i=n1,n2
              i1=i-n1+1
              ni=Nip(i)
              l=Ll(ni)+1
              lll(i1)=let(l)
              jjj(i1)=Jj(ni)
              nnn(i1)=Nn(ni)
              nqq(i1)=Nq(i)
              if (Nq(i) > jjj(i1)+1) then
                 write(11,'(/2X,"wrong number of electrons"/ &
                      2X,"for the shell:",I3,3X,I2,A1,I2,"/2", &
                      " (",I2,")")') ni,nnn(i1),lll(i1),jjj(i1),nqq(i1)
               Stop
              end if
           end do
           n=n2-n1+1
           write(11,'(1X,I6,"#",6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"), &
                /8X,6(I2,A1,"(",I1,"/2)",I2,";"))') &
                ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
        end do
        write(11,'(1X,71("="))')
        if(Ksig > 0) then
           do ni=Nso+1,Nsu
              read(12,rec=2*ni+7) p
              Eps(ni)=-p(ii+1)
           end do
           write(11,'(" HF energies are read from DAT file", &
                  /5(I5,F10.6))') (i,Eps(i),i=Nso+1,Nsu)
        end if
        close(unit=12)
        ! Maximal number of eigenvectors Nlv:
        nlmax=IPlv
        if (Kv >= 3) then
            nlmax=IPlv/3
        end if
        if (Nlv > nlmax) Nlv=nlmax
       return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( 6,'(/2X,"file CONF.DAT is absent"/)')
        write(11,'(/2X,"file CONF.DAT is absent"/)')
       Stop
!     - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine Init

    subroutine InitFormH(npes,mype)
      ! this subroutine initializes variables used for FormH and subsequent subroutines
      ! All necessary variables are broadcasted from root to all cores
        use mpi_f08
        implicit none
        integer :: npes, mype, mpierr, i
        if (mype==0) print*, 'start InitFormH'
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(N_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Crt4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kherr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kgerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kecp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nsu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Mj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(NmaxS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(LmaxS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ksym, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nsum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Gj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(C_is, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(K_sms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kdsig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kexn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Eps(1:IPs), IPs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(In(1:IPgnt), IPgnt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Gnt(1:IPgnt), IPgnt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Scr, 10, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nn(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kk(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ll(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jj(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Rint1(1:Nhint), Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Rint2(1:IPbr,1:Ngint), IPbr*Ngint, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint1(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint2(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint3(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(IntOrd, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        if (Ksig /= 0) then
          call MPI_Bcast(Rsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Dsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Esig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Iint1S(1:NhintS), NhintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Iint2S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Iint3S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Rint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Dint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(Eint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(IntOrdS, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(R_is(1:num_is), num_is, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
          call MPI_Bcast(I_is(1:num_is), num_is, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        end if
        !do i=1,Ne
        !    call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        !end do
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        if (mype==0) print*, 'InitFormH done'
        return
    end subroutine InitFormH

    subroutine InitMxmpy(npes,mype)
      ! this subroutine sets up sizes and displacements for MPI during Mxmpy subroutine
      use mpi_f08
      implicit none
      integer :: npes, mype, mpierr
      integer :: i
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(H_nsizes(npes), H_ndisps(npes))
      H_nstart=H_n(1)
      H_nend=H_n(ih8H)
      H_nsize=H_nend-H_nstart+1

      call MPI_Gather(H_nsize, 1, MPI_INTEGER, H_nsizes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

      H_ndisps(1) = 0
      do i=2,npes
        H_ndisps(i) = sum(H_nsizes(1:i-1))
      end do
      return
    end subroutine
    
end module conf_init