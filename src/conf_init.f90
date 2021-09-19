Module conf_init

    Use conf_variables

    Implicit None

    Private

    Public :: ReadConfInp, ReadConfigurations, inpstr, Init, InitFormH

  Contains

    Subroutine ReadConfInp
        ! This subroutine Reads input parameters from CONF.INP
        Implicit None
        Integer :: istr
        Character(Len=1) :: name(16)
        ! - - - - - - - - - - - - - - - - - - - - - -
        Open(unit=10,file='CONF.INP',status='OLD')
        Read(10,'(1X,16A1)') name
        Write(*,'(4X,16A1)') name
        Write(11,'(4X,16A1)') name
        Read(10,'(5X,F5.1)') Z, Am, XJ_av, Jm
        Read(10,'(5X,I6)') Nso, Nc, Kv, Nlv, Ne

        IPlv = 3*Nlv
        K_is = 0                ! 
        Kbrt = 0                !  
        Kout = 1                !   Optional
        Kecp = 0                !   parameters
        C_is = 0.d0             !      
        Gj   = 0.d0             !   
        n_it = 20               !
        Cut0 = 0.001            !       
        Ncpt = 0                !  
        Gnuc=1.d0               !
        Qnuc=1.d0               !

        istr = 0
        Do While (istr /= 1)
          Call inpstr(istr)  
        End Do    

        Return
    End subroutine ReadConfInp

    Subroutine ReadConfigurations
        ! This subroutine Reads configurations from CONF.INP
        Implicit None
        Integer  :: i, i1, i2, ic, ne0, nx, ny, nz
        Real(dp) :: x
        ! - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Qnl(100000000)) ! upper bound = 100 million rel. conf-s = 0.8 GB
        If (Nso /= 0) Then 
            Read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso) ! Read core conf-s
        End If
        ! Reading in configurations - - - - - - 
        i1=Nso+1
        Do ic=1,Nc
            ne0=0
            Do While (ne0 < Ne)
                i2=i1+5
                Read (10,'(6(4X,F7.4))') (Qnl(i),i=i1,i2)
                Do i=i1,i2
                   x=abs(Qnl(i))+1.d-9
                   If (x < 1.d-8) Exit
                   nx=10000*x
                   ny=100*x
                   nz=(nx-100*ny)
                   ne0=ne0+nz
                End Do
                i2=i-1
                i1=i2+1
            End Do
            If (ne0 > Ne) Then
                Write(6,'(" INPUT: too many electrons for ic =",I6)') ic
                Stop
            End If
        End Do
        Nsp=i2
        Close(unit=10)
        Return
    End subroutine ReadConfigurations

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
    200 If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 210
        If (txt(2) /= '_') goto 210
        If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 210
        If (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 210
        Read (string,*) K_is
        Write( *,'(5a1,i6)') (txt(i),i=1,5),K_is
        Write(11,'(5a1,i6)') (txt(i),i=1,5),K_is
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    210 If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 220
        If (txt(2) /= '_') goto 220
        If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 220
        If (txt(4) /= 's'   .and.   txt(4) /= 'S') goto 220
        Read (string,*) C_is
        Write( *,'(5a1,F6.3)') (txt(i),i=1,5),C_is
        Write(11,'(5a1,F6.3)') (txt(i),i=1,5),C_is
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    220 If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 230
        If (txt(2) /= 'l'   .and.   txt(1) /= 'L') goto 230
        If (txt(3) /= 'o'   .and.   txt(3) /= 'O') goto 230
        If (txt(4) /= 'w'   .and.   txt(4) /= 'W') goto 230
        Read (string,*) Klow
        Write( *,'(5a1,i6)') (txt(i),i=1,5),Klow
        Write(11,'(5a1,i6)') (txt(i),i=1,5),Klow
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    230 If (txt(1) /= ' ') goto 240
        If (txt(2) /= 'g'   .and.   txt(2) /= 'G') goto 240
        If (txt(3) /= 'j'   .and.   txt(3) /= 'J') goto 240
        If (txt(4) /= ' ') goto 240
        Read (string,*) Gj
        Write( *,'(5a1,F6.3)') (txt(i),i=1,5),Gj
        Write(11,'(5a1,F6.3)') (txt(i),i=1,5),Gj
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    240 If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 250
        If (txt(2) /= 'o'   .and.   txt(2) /= 'O') goto 250
        If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 250
        If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 250
        Read (string,*) kout
        Write( *,'(5a1,i6)') (txt(i),i=1,5),kout
        Write(11,'(5a1,i6)') (txt(i),i=1,5),kout
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    250 If (txt(1) /= 'k'   .and.   txt(1) /= 'K') goto 260
        If (txt(2) /= 'e'   .and.   txt(2) /= 'E') goto 260
        If (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 260
        If (txt(4) /= 'p'   .and.   txt(4) /= 'P') goto 260
        Read (string,*) kecp
        Write( *,'(5a1,i6)') (txt(i),i=1,5),kecp
        Write(11,'(5a1,i6)') (txt(i),i=1,5),kecp
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    260 If (txt(1) /= 'g'   .and.   txt(1) /= 'G') goto 270
        If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 270
        If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 270
        If (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 270
        Read (string,*) gnuc
        Write( *,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
        Write(11,'(5a1,F8.5)') (txt(i),i=1,5),gnuc
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    270 If (txt(1) /= 'q'   .and.   txt(1) /= 'Q') goto 280
        If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 280
        If (txt(3) /= 'u'   .and.   txt(3) /= 'U') goto 280
        If (txt(4) /= 'c'   .and.   txt(4) /= 'C') goto 280
        Read (string,*) qnuc
        Write( *,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
        Write(11,'(5a1,F8.5)') (txt(i),i=1,5),qnuc
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    280 If (txt(2) /= 'k'   .and.   txt(2) /= 'K') goto 290
        If (txt(3) /= 'l'   .and.   txt(3) /= 'L') goto 290
        If (txt(4) /= '4') goto 290
        Read (string,*) kl4
        Write( *,'(5a1,i6)') (txt(i),i=1,5),kl4
        Write(11,'(5a1,i6)') (txt(i),i=1,5),kl4
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    290 If (txt(2) /= 'n'   .and.   txt(2) /= 'N') goto 300
        If (txt(3) /= 'c'   .and.   txt(3) /= 'C') goto 300
        If (txt(4) /= '4') goto 290
        Read (string,*) nc4
        Write( *,'(5a1,i6)') (txt(i),i=1,5),nc4
        Write(11,'(5a1,i6)') (txt(i),i=1,5),nc4
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    300 If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 310
        If (txt(2) /= 'r'   .and.   txt(2) /= 'R') goto 310
        If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 310
        If (txt(4) /= '4') goto 310
        Read (string,*) crt4
        Write( *,'(5a1,F8.5)') (txt(i),i=1,5),crt4
        Write(11,'(5a1,F8.5)') (txt(i),i=1,5),crt4
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    310 If (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 320
        If (txt(2) /= '_') goto 320
        If (txt(3) /= 'i'   .and.   txt(3) /= 'I') goto 320
        If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 320
        Read (string,*) n_it
        Write( *,'(5a1,i6)') (txt(i),i=1,5),n_it
        Write(11,'(5a1,i6)') (txt(i),i=1,5),n_it
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    320 If (txt(1) /= 'a'   .and.   txt(1) /= 'A') goto 330
        If (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 330
        If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 330
        If (txt(4) /= 'o'   .and.   txt(4) /= 'O') goto 330
        Read (string,*) kautobas
        Write( *,'(5a1,i6)') (txt(i),i=1,5),kautobas
        Write(11,'(5a1,i6)') (txt(i),i=1,5),kautobas
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    330 If (txt(1) /= 'c'   .and.   txt(1) /= 'C') goto 340
        If (txt(2) /= 'u'   .and.   txt(2) /= 'U') goto 340
        If (txt(3) /= 't'   .and.   txt(3) /= 'T') goto 340
        If (txt(4) /= '0') goto 340
        Read (string,*) Cut0
        Write( *,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
        Write(11,'(5a1,F8.5)') (txt(i),i=1,5),Cut0
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    340 If (txt(1) /= 'n'   .and.   txt(1) /= 'N') goto 350
        If (txt(2) /= 'c'   .and.   txt(2) /= 'C') goto 350
        If (txt(3) /= 'p'   .and.   txt(3) /= 'P') goto 350
        If (txt(4) /= 't'   .and.   txt(4) /= 'T') goto 350
        Read (string,*) Ncpt
        Write( *,'(5a1,i6)') (txt(i),i=1,5),Ncpt
        Write(11,'(5a1,i6)') (txt(i),i=1,5),Ncpt
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    350 If (txt(2) /= ' '   .and.   txt(2) /= '-') goto 700
        If (txt(3) /= ' '   .and.   txt(3) /= '-') goto 700
        If (txt(4) /= ' '   .and.   txt(4) /= '-') goto 700
        If (txt(2) == ' '   .and.   txt(3) == ' ') backspace(10)
        istr=1
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    700 istr=2
        Write( *,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
        Write(11,'(/2x,"Unknown parameter in input file: ",5a1)') (txt(i),i=1,5)
        Read(*,*)
        Stop
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
    End subroutine inpstr

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

    Subroutine InitFormH(npes,mype)
        ! this subroutine initializes variables used for FormH and subsequent subroutines
        ! All necessary variables are broadcasted from root to all cores
        Use mpi
        Implicit None
        Integer :: npes, mype, mpierr, i

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Crt4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kl4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
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
        Call MPI_Bcast(In(1:IPgnt), IPgnt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Gnt(1:IPgnt), IPgnt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kk(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint1(1:Nhint), Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Rint2(1:IPbr,1:Ngint), IPbr*Ngint, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint1(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint2(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iint3(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IntOrd, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Do i=1,Ne
            Call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End do  
        If (Ksig /= 0) Then
            Call MPI_Bcast(Scr, 10, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dsig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Esig(1:NhintS), NhintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint1S(1:NhintS), NhintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint2S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Iint3S(1:NgintS), NgintS, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Rint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Dint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(Eint2S(1:NgintS), NgintS, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(IntOrdS, IPx*IPx, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(R_is(1:num_is), num_is, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(I_is(1:num_is), num_is, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End subroutine InitFormH
    
End module conf_init