Module conf_init

    Use conf_variables
    Use utils, Only : ToUpperString

    Implicit None

    Private

    Public :: ReadConfInp, ReadConfigurations, inpstr, ReadCiIn

  Contains
    
    Subroutine PrintParamI(key, val)
        ! This subroutine prints a parameter from a key-value pair
        Implicit None
        character(len=4), intent(in) :: key
        integer, intent(in) :: val

        Write( *,'(a5,i6)') adjustr(key) // '=', val
        Write(11,'(a5,i6)') adjustr(key) // '=', val

    End Subroutine PrintParamI

    Subroutine PrintParamR(key, val, strfmt)
        ! This subroutine prints a parameter from a key-value pair
        Implicit None
        character(len=4), intent(in) :: key
        character(len=9), intent(in) :: strfmt
        real(dp), intent(in) :: val

        Write( *, strfmt) adjustr(key) // '=', val
        Write(11, strfmt) adjustr(key) // '=', val

    End Subroutine PrintParamR

    Subroutine ReadConfParams
        ! This subroutine reads parameters from the header of CONF.INP
        Implicit None

        integer :: index_equals, index_hashtag
        character(len=4) :: key, key_upper
        character(len=10) :: val
        character(len=80) :: line
        logical :: equals_in_str

        ! roll back to beginning of CONF.INP
        Rewind(10)

        ! read the first line
        Read(10, '(A)') 

        ! read parameters (lines with "=")
        equals_in_str = .true.
        Do While (equals_in_str)
            Read(10, '(A)') line
            If (index(string=line, substring="=") == 0) Then
                equals_in_str = .false.
            Else
                index_equals = index(string=line, substring="=")
                key = trim(adjustl(line(1:index_equals-1)))
                val = trim(adjustl(line(index_equals+1:len(line))))
                index_hashtag = index(string=val, substring="#") ! account for comments
                If (index_hashtag /= 0) val = trim(adjustl(val(1:index_hashtag-1)))

                ! normalize key to uppercase
                key_upper = ToUpperString(key)
                
                Select Case(key_upper)
                Case('Z')
                    Read(val, *) Z
                Case('AM')
                    Read(val, *) Am
                Case('J')
                    Read(val, *) XJ_av
                Case('JM')
                    Read(val, *) Jm
                Case('NSO')
                    Read(val, *) Nso
                Case('NC')
                    Read(val, *) Nc
                Case('KV')
                    Read(val, *) Kv
                Case('NLV')
                    Read(val, *) Nlv
                Case('NE')
                    Read(val, *) Ne
                Case('KL4')
                    Read(val, *) Kl4
                    Call PrintParamI(key, Kl4)
                Case('NC4')
                    Read(val, *) Nc4
                    Call PrintParamI(key, Nc4)
                Case('GJ')
                    Read(val, *) Gj
                    Call PrintParamR(key, Gj, '(A5,F6.3)')
                Case('CRT4')
                    Read(val, *) Crt4
                    Call PrintParamR(key, Crt4, '(A5,F8.5)')
                Case('KOUT')
                    Read(val, *) Kout
                    Call PrintParamI(key, Kout)
                Case('NCPT')
                    Read(val, *) Ncpt
                    Call PrintParamI(key, Ncpt)
                Case('CUT0')
                    Read(val, *) Cut0
                    Call PrintParamR(key, Cut0, '(A5,F8.5)')
                Case('N_IT')
                    Read(val, *) N_it
                    Call PrintParamI(key, N_it)
                Case('KBRT')
                    Read(val, *) Kbrt
                    Call PrintParamI(key, Kbrt)
                Case('K_IS')
                    Read(val, *) K_is
                    Call PrintParamI(key, K_is)
                Case('C_IS')
                    Read(val, *) C_is
                    Call PrintParamR(key, C_is, '(A5,F7.4)')
                Case('KLOW')
                    Read(val, *) Klow
                    Call PrintParamI(key, Klow)
                Case('GNUC')
                    Read(val, *) Gnuc
                    Call PrintParamR(key, Gnuc, '(A5,F8.5)')
                Case('QNUC')
                    Read(val, *) Qnuc
                    Call PrintParamR(key, Qnuc, '(A5,F8.5)')
                End Select
            End If
        End Do
        
        Return
    End Subroutine ReadConfParams

    Subroutine ReadCiIn
        ! This subroutine reads job parameters from file ci.in
        Implicit None

        integer :: index_equals, index_hashtag, err_stat
        character(len=10) :: key
        character(len=10) :: val
        character(len=80) :: line
        logical :: equals_in_str

        Open(unit=99, file='ci.in', status='OLD')

        ! Set default values for keys
        Kl = 0
        Ksig = 0
        Kdsig = 0
        Kw = 0
        KLSJ = 0
        K_sms = 3
        KWeights = 0
        KXIJ = 10
        MaxNd0 = 3000
        
        ! read parameters (lines with "=")
        equals_in_str = .true.
        Do While (equals_in_str)
            Read(99, '(A)', iostat=err_stat) line
            If (index(string=line, substring="=") == 0 .or. err_stat /= 0) Then
                equals_in_str = .false.
            Else
                index_equals = index(string=line, substring="=")
                key = line(1:index_equals-1)
                val = line(index_equals+1:len(line))
                
                index_hashtag = index(string=val, substring="#") ! account for comments
                If (index_hashtag /= 0) val = trim(adjustl(val(1:index_hashtag-1)))
                Select Case(key)
                Case('Kl')
                    ! Kl = 0 - new computation
                    ! Kl = 1 - continuing computation with completed CONF.HIJ and CONF.JJJ files
                    ! Kl = 2 - new computation with MBPT
                    ! Kl = 3 - extending computation with new configurations (not implemented yet)
                    ! Kl = 4 - compute minimum memory requirements for pconf
                    Read(val, *) Kl
                Case('Ksig')
                    ! If starting new computation with MBPT
                    ! Ksig = 0 - no MBPT included (same as Kl = 0)
                    ! Ksig = 1 - include 1-electron MBPT corrections 
                    ! Ksig = 2 - include 1-electron and 2-electron MBPT corrections
                    Read(val, *) Ksig
                Case('Kdsig')
                    ! Kdsig = 0 - automatic approximation of the energy dependence of Sigma
                    ! Kdsig = 1 - manually include energy dependence of Sigma
                    Read(val, *) Kdsig
                Case('Kw')
                    ! Kw determines whether CONF.HIJ will be written or not
                    ! Kw=0 - CONF.HIJ will not be written
                    ! Kw=1 - CONF.HIJ will be written
                    Read(val, *) Kw
                Case('KLSJ')
                    ! KLSJ determines whether CONF.HIJ will be written or not
                    ! KLSJ=0 - LSJ will not be written in FINAL.RES
                    ! KLSJ=1 - LSJ will be written FINAL.RES
                    Read(val, *) KLSJ
                Case('K_sms')
                    ! SMS to include 1-e (1), 2-e (2), or both (3)
                    Read(val, *) K_sms
                    ! Write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): ', K_sms
                Case('KXIJ')
                    ! KXIJ determines the interval in which CONF.XIJ will be written
                    ! e.g. KXIJ=10 - CONF.XIJ is written every 10 Davidson iterations
                    Read(val, *) KXIJ
                Case('KWeights')
                    ! KWeights determines whether CONF.WGT is written (1) or not (0)
                    Read(val, *) KWeights
                Case('MaxNd0')
                    ! MaxNd0 determines the size of the initial approximation in determinants
                    Read(val, *) MaxNd0
                End Select
            End If
        End Do

        Close(99)

        Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
        Write( 6,'(/4X,"Kw = (0-do not write CONF.HIJ, 1-write CONF.HIJ) ",I1)') Kw
        Write( 6,'(/4X,"KLSJ = (0-do not calculate LSJ, 1-calculate LSJ) ",I1)') KLSJ

    End Subroutine ReadCiIn 

    Subroutine ReadConfInp
        ! This subroutine reads input parameters from the header of CONF.INP
        Implicit None
        
        Character(Len=1) :: name(16)

        ! Optional parameters
        K_is = 0 
        C_is = 0.d0 
        Klow = 0
        Kbrt = 0 
        Kout = 1 
        Kecp = 0 
        Gj   = 0.d0 
        n_it = 20 
        Cut0 = 0.001 
        Ncpt = 0 
        Gnuc = 1.d0 
        Qnuc = 1.d0  

        Open(unit=10,file='CONF.INP',status='OLD')
        Read(10,'(1X,16A1)') name
        Write(*,'(4X,16A1)') name
        Write(11,'(4X,16A1)') name
        Call ReadConfParams

        ! Dimension of ArrB
        IPlv = 3*Nlv

        Return
    End subroutine ReadConfInp

    Subroutine GotoConfigurations
        ! This subroutine reads to the start of the list of configurations in CONF.INP
        Implicit None

        character(len=80) :: line
        logical :: equals_in_str

        ! roll back to beginning of CONF.INP
        Rewind(10)

        ! read the first line
        Read(10, '(A)') 

        ! read parameters (lines with "=")
        equals_in_str = .true.
        Do While (equals_in_str)
            Read(10, '(A)') line
            If (index(string=line, substring="=") == 0) Then
                If (line /= "") Backspace(10) ! handle blank line before list of core orbitals
                equals_in_str = .false.
            End If
        End Do

        Return
    End Subroutine GotoConfigurations

    Subroutine ReadConfigurations
        ! This subroutine Reads configurations from CONF.INP
        Implicit None
        Integer  :: i, i1, i2, ic, ne0, nx, ny, nz
        Real(dp) :: x
        
        Call calcNsp
        Call GotoConfigurations

        Allocate(Qnl(Nsp+6)) ! the '6' is a buffer for number of orbitals per line for conf-s
        
        If (Nso /= 0) Then 
            Read (10,'(6(4X,F7.4))') (Qnl(i),i=1,Nso) ! Read core conf-s
        End If
        Nlx=0
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
                   nx=Int(10000*x)
                   ny=Int(100*x)
                   nz=(nx-100*ny)
                   ne0=ne0+nz
                   If (mod(ny,10) > Nlx) Nlx=mod(ny,10) ! Set Nlx to be max partial wave
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

    Subroutine calcNsp
        ! This subroutine calculates the total number of orbitals in conf-s
        Implicit None

        Integer :: i, j, ic, ne0, nx, ny, nz, num_orbitals_per_line, counter, num_lines
        Real(dp) :: x
        Real(dp), Dimension(:), Allocatable :: line

        num_orbitals_per_line = 6
        Allocate(line(num_orbitals_per_line))

        counter = 0

        ! Read core orbitals
        If (Nso /= 0) Then 
            x = Nso/6
            num_lines = Ceiling(x)
            Do j=1,num_lines
                Read (10,'(6(4X,F7.4))') (line(i), i=1,num_orbitals_per_line) 
            End Do
            counter = counter + Nso
        End If

        ! Read configuration list
        Do ic=1,Nc
            ne0=0
            Do While (ne0 < Ne)
                Read (10,'(6(4X,F7.4))') (line(i), i=1,num_orbitals_per_line)
                Do i=1,num_orbitals_per_line
                   x=abs(line(i))+1.d-9
                   If (x < 1.d-8) Exit
                   nx=Int(10000*x)
                   ny=Int(100*x)
                   nz=(nx-100*ny)
                   ne0=ne0+nz
                   counter = counter + 1
                End Do
            End Do
        End Do

        Deallocate(line)
        Nsp = counter

        Return
    End Subroutine calcNsp

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
        Write( *,'(5a1,F7.4)') (txt(i),i=1,5),C_is
        Write(11,'(5a1,F7.4)') (txt(i),i=1,5),C_is
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