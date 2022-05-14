Module conf_init

    Use conf_variables

    Implicit None

    Private

    Public :: ReadConfInp, ReadConfigurations, inpstr

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
        Integer  :: i, i1, i2, ic, ne0, nx, ny, nz, nr, cnt
        Real(dp) :: x
        ! - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Qnl(100000000)) ! upper bound = 100 million rel. conf-s = 0.8 GB
        Allocate(Nrnrc(Nc))
        If (Nso /= 0) Then 
            Read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso) ! Read core conf-s
        End If
        ! Reading in configurations - - - - - - 
        i1=Nso+1
        nr=0
        Do ic=1,Nc
            ne0=0
            Do While (ne0 < Ne)
                i2=i1+5
                Read (10,'(I4,F7.4,5(4X,F7.4))') nr, Qnl(i1), (Qnl(i),i=i1+1,i2)
                If (Qnl(i1) == 0) Then
                    Nnr = nr
                    Nrnrc(Nnr-1) = cnt
                    cnt = 0
                End If
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
            cnt = cnt + 1
            If (ne0 > Ne) Then
                Write(6,'(" INPUT: too many electrons for ic =",I6)') ic
                Stop
            End If
        End Do
        Nrnrc(Nnr) = cnt
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

End module conf_init