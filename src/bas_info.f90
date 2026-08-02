program bas_info         ! this code calculates rms radii and tests the sum rule <a|r|n><n|r|a>=<a|r^2|a>
    use params
    use basc_variables, r2old => r2
    use utils, only : DetermineRecordLength, CheckRecordLength

    implicit none

    integer :: iyes, num
    real(dp) :: rrc(5,IPs), r2(IPs)

    call init_basis
    call print_rms_radii             !### rms radii of orbitals
    call check_unity_expansion       !### Check sum of orbitals in p.w.
    write(*,*)

    do while (.true.)
        write(*,*) ' Choose optional tests (0 to quit)'
        write(*,*) ' 1. closure r^2=r|n><n|r is checked for each'
        write(*,*) '    orbital and for three sets of'
        write(*,*) '    intermediate states |n> which'
        write(*,*) '    meets dipole selection rules.'
        write(*,*) ' 2. identity int_0^infty{1/2*f*df/dr}=0 is checked'
        write(*,*) '    for each orbital'
        write(*,*) ' 3. Checks orthogonality'
        write(*,*) ' 4. Checks Taylor expansion at the origin'
        write(*,*) ' 5. prints any orbital of your choice'
        read(*,*) iyes

        if (iyes <= 0) stop
        if (iyes == 1) call check_closure          !### Check closure r|n><n|r=r^2
        if (iyes == 2) call check_derivative       !### Check identity int_0^infty{1/2*f*df/dr}=0
        if (iyes == 3) call check_orthogonality    !### Check orthonormality
        if (iyes == 4) call check_taylor_expansion !### Check Taylor expansion at the origin
        if (iyes == 5) then                        !### Print orbital of your choice
            do while (.true.)
                write(*,*) ' Orbital to print (0-quit) '
                read(*,*) num
                if (num <= 0 .or. num > Ns) stop
                call print_orbital(num)
            end do
        end if
    end do

contains

    subroutine init_basis
        !! Initializes the basis set by reading header records from the DAT file.
        !!
        !! Prompts the user for the DAT filename (default HFD.DAT; '1' = HFD.DAT, '2' = CONF.DAT).
        !! Opens BAS_INFO.RES for output and the DAT file for direct-access binary I/O.
        !! Reads records 1-3 into a flat PQ array mirroring the original EQUIVALENCE layout:
        !!   record 1 -> PQ(1:IP6) = P,  PQ(IP6+1:2*IP6) = Q   (header: Z, Ns, ii, h, Kt, ...)
        !!   record 2 -> R, V
        !!   record 3 -> PQ(2*IP6+1:3*IP6) = P1 (derivative P'), PQ(3*IP6+1:4*IP6) = Q1
        !! Extracts scalar parameters (Z, Ns, ii, h, Kt) and orbital quantum numbers Nn, Ll, Jj from the header.
        !! Two layouts are supported:
        !!   longbasis=.true.  -- quantum numbers packed as integers starting at PQ(21)
        !!   longbasis=.false. -- quantum numbers packed as reals, 6 words per orbital from PQ(21)
        !! Sets MaxT=9 (number of Taylor coefficients stored per orbital).
        use readfff, only : ReadF
        implicit none

        integer    :: ni, ios, idx
        real(dp)   :: PQ(4*IP6)   ! flat P|Q|P'|Q' concatenation, mirrors original EQUIVALENCE
        integer    :: IQN(4*IPs)  ! integer reinterpretation of PQ(21:...) for longbasis
        logical    :: longbasis
        real(dp), parameter :: round_offset = 0.01d0   ! added before int truncation of packed reals
        character(len=12)   :: fname

        call DetermineRecordLength(Mrec)
        MaxT = 9

        write(*,*)
        write(*,*) ' This program gives information about basis set:'
        write(*,*) ' list of orbitals, their energies and rms radii'
        write(*,*)
        write(*,*) ' Optional tests see below'

        fname = 'HFD.DAT'
        write(*,*) ' Give file name or ENTER to quit: '
        read(*,'(A12)') fname
        if (fname <= ' ')   stop
        if (fname == '1')   fname = 'HFD.DAT'
        if (fname == '2')   fname = 'CONF.DAT'

        ! open output file (delete any previous run)
        open(11, file='BAS_INFO.RES', status='UNKNOWN')
        close(11, status='DELETE')
        open(11, file='BAS_INFO.RES', status='NEW')

        open(12, file=fname, access='DIRECT', status='OLD', recl=IP6*Mrec, iostat=ios)
        if (ios /= 0) then
            write(*,*) ' No file ', fname, ' SORRY!'
            stop
        end if
        call CheckRecordLength(12, fname, IP6*Mrec)
        write(11,*) ' >>>> Basis set from ', fname, ' <<<<'
        write(*,*)  ' Output see in "BAS_INFO.RES" file'

        ! read header records into flat PQ array
        call ReadF(12, 1, PQ(1),         PQ(IP6+1),     2)  ! record 1: P header / Q header
        call ReadF(12, 2, R,              V,            2)  ! record 2: radial grid / potential
        call ReadF(12, 3, PQ(2*IP6+1),   PQ(3*IP6+1),   2)  ! record 3: P' / Q'

        ! copy into module arrays (host-associated)
        P = PQ(1:IP6)
        Q = PQ(IP6+1:2*IP6)

        ! extract scalar parameters from record 1 header
        Z  = PQ(1)
        Ns = int(PQ(2) + round_offset)
        ii = int(PQ(3) + round_offset)
        h  = PQ(6)
        Kt = int(PQ(9) + round_offset)
        longbasis = abs(PQ(20) - 0.98765d0) < 1.d-6

        ! extract orbital quantum numbers
        if (longbasis) then
            write( *,*) ' Using variant for long basis '
            write(11,*) ' Using variant for long basis '
            ! IQN overlays PQ starting at index 21 (mirrors original EQUIVALENCE (IQN(1),PQ(21)))
            IQN = transfer(PQ(21:21+size(IQN)/2), IQN)
            do ni = 1, Ns
                Nn(ni) = IQN(4*ni-3)
                Ll(ni) = IQN(4*ni-2)
                Jj(ni) = IQN(4*ni)
            end do
        else
            ! quantum numbers packed as reals: 6 words per orbital starting at PQ(21)
            idx = 20
            do ni = 1, Ns
                idx    = idx + 1;  Nn(ni) = int(PQ(idx) + round_offset)
                idx    = idx + 1;  Ll(ni) = int(PQ(idx) + round_offset)
                idx    = idx + 4;  Jj(ni) = int(PQ(idx) + round_offset)
            end do
        end if
    end subroutine init_basis


    subroutine print_rms_radii
        !! Prints root-mean-square radii and orbital energies for all basis orbitals.
        !!
        !! For each orbital ni, reads P and Q from the DAT file and computes:
        !!   - norm  = int (P^2 + Q^2) dr  (should be 1; a warning is printed if |norm-1| > 1e-3)
        !!   - <r^2> = int (P^2 + Q^2) r^2 dr / norm
        !!   - R_rms = sqrt(<r^2>)
        !!   - nnod  = number of nodes in P
        !! Also stores <r^2> in r2(ni) for later use in check_closure.
        !! Results are written to unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use sintg, only : Sint1
        implicit none

        integer :: ni, num_nodes, i
        real(dp) :: pr, s, x

        character(len=*), parameter :: fmt_hdr = &
               "(1X,'Z =',F5.1,', h0 =',F9.5,', Ns =',I4,', ii =',I4," // &
               "/1X,'Radial grid from R_1 =',F11.8,' to R_ii =',F9.5," // &
               "/1X,67('-'),/15X,'Root-mean-square radii:',/1X,67('-')," // &
               "/6X,'orbital',6X,'energy',8X,'P(ii+5)',7X,'Q(ii+5)'," // &
               "5X,'R_rms',2X,'Nodes')"
        character(len=*), parameter :: fmt_orb = &
               "(1X,I3,':',I3,A1,I2,'/2'," // &
               "F15.8,2(1X,E13.6),1X,F8.5,1X,I3)"
        character(len=*), parameter :: fmt_warn = &
               "(1X,I3,':',I3,A1,I2,'/2 E='," // &
               "F13.6,' norma=',E13.6,' <r>=',F8.5,' push..')"

        write(11,fmt_hdr) Z, h, Ns, ii, R(1), R(ii)
        do ni = 1, Ns
            call ReadF(12, ni+4, P, Q, 2)

            ! count nodes and compute norm = int (P^2 + Q^2) dr
            num_nodes = 0
            pr        = 0.d0
            do i = 1, ii
                C(i) = P(i)**2 + Q(i)**2
                if (P(i)*pr < 0.d0) num_nodes = num_nodes + 1
                if (P(i) /= 0.d0)   pr        = P(i)
            end do
            C(ii+4) = P(ii+4)*2
            call Sint1(s)                          ! s = norm

            ! compute <r^2> = int (P^2 + Q^2) r^2 dr
            do i = 1, ii
                C(i) = (P(i)**2 + Q(i)**2) * R(i)**2
            end do
            C(ii+4) = P(ii+4)*2 + 2.d0
            call Sint1(x)

            x       = x / s                        ! <r^2>
            r2(ni)  = x
            x       = sqrt(x)                      ! R_rms
            s       = sqrt(s)                      ! sqrt(norm)

            write(11,fmt_orb) ni, Nn(ni), OrbLabels(Ll(ni)+1), Jj(ni), -P(ii+1), P(ii+5), Q(ii+5), x, num_nodes
            if (abs(s - 1.d0) > 1.d-3) then
                write(*,fmt_warn) ni, Nn(ni), OrbLabels(Ll(ni)+1), Jj(ni), -P(ii+1), s, x
                read(*,*)
            end if
        end do
    end subroutine print_rms_radii


    subroutine check_closure
        !! Checks the dipole sum rule <a|r^2|a> = sum_n <a|r|n><n|r|a>.
        !!
        !! For each orbital a (ni), computes the RHS by summing |<a|r|n>|^2 over all
        !! intermediate orbitals n that satisfy the dipole selection rules
        !! (Delta_l = 0, +/-1;  Delta_j = 0, +/-1).  The sum is split into up to 5
        !! intermediate (l,2j) channels covering l-1, l, l+1.
        !! The printed ratio RHS/LHS = 1 when the basis is complete for that channel.
        !! Results are written to unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use sintg, only : Sint1
        implicit none

        integer :: ni, ln, lx, ik, il, jn, jx, ij, nk, ierr, i, li, ikx
        real(dp) :: s, y
        real(dp) :: P1(IP6), Q1(IP6)
        integer       :: jk_chan(5,IPs), numk(IPs)   ! 2j and count for each channel
        character(len=1) :: lk_chan(5,IPs)            ! l-label for each channel

        character(len=*), parameter :: fmt_hdr = &
               "(1X,65('-'),/15X,'Test of the sum rule r^2=r|n><n|r'," // &
               "/10X,'Given is the ratio RHS/LHS for each orbital'," // &
               "/10X,'and intermediate sets with L_1 = L-1,L,L+1:'," // &
               "/1X,65('-'))"
        character(len=*), parameter :: fmt_row = &
               "(I3,A1,I2,'/2 ',5(' # ',A1,I2,'/2',F6.3))"

        do ni = 1, Ns
            call ReadF(12, ni+4, P, Q, 2)
            ln = max(Ll(ni) - 1, 0)
            lx = Ll(ni) + 1
            ik = 0
            do il = ln, lx                         ! intermediate l = l-1, l, l+1
                jn = max(Jj(ni)-2, 1, 2*il-1)
                jx = min(Jj(ni)+2,   2*il+1)
                do ij = jn, jx, 2                  ! intermediate 2j
                    ik = ik + 1
                    lk_chan(ik,ni) = OrbLabels(il+1)
                    jk_chan(ik,ni) = ij
                    s = 0.d0
                    do nk = 1, Ns
                        ierr = abs(Ll(nk)-il) + abs(Jj(nk)-ij)
                        if (ierr == 0) then
                            call ReadF(12, nk+4, P1, Q1, 2)
                            do i = 1, ii
                                C(i) = (P(i)*P1(i) + Q(i)*Q1(i)) * R(i)
                            end do
                            C(ii+4) = P(ii+4) + P1(ii+4) + 1.d0
                            call Sint1(y)
                            s = s + y*y
                        end if
                    end do
                    rrc(ik,ni) = s / r2(ni)
                    numk(ni)   = ik
                end do
            end do
        end do

        write(11,fmt_hdr)
        do ni = 1, Ns
            li  = Ll(ni) + 1
            ikx = numk(ni)
            write(11,fmt_row) Nn(ni), OrbLabels(li), Jj(ni), &
                (lk_chan(ik,ni), jk_chan(ik,ni), rrc(ik,ni), ik=1,ikx)
        end do
    end subroutine check_closure


    subroutine check_derivative
        !! Checks the identity int_0^infty (P*dP/dr + Q*dQ/dr) dr = 0 for each orbital.
        !!
        !! For each orbital ni, reads P, Q and their stored derivatives CP, CQ from the DAT file (record ni+Ns+4). 
        !! Evaluates the integral two ways:
        !!   - using the numerical derivative (via Dif)
        !!   - using the derivative stored in the file
        !! Both results are normalised by max(P^2+Q^2) and printed side by side.
        !! Results are written to unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use sintg, only : Sint1
        use diff, only : Dif
        implicit none

        integer :: ni, i
        real(dp) :: fx, ci, c2, x, x1
        real(dp) :: P1(IP6), Q1(IP6), CP(IP6), CQ(IP6)

        character(len=*), parameter :: fmt_hdr = &
               "(1X,65('-'),/10X,'Test of the equation '," // &
               "'int_0^infty(P*P`+Q*Q`)dx = 0'," // &
               "/10X,'for numerical derivatives and those from the file'," // &
               "/17X,'num. dif.  file ')"
        character(len=*), parameter :: fmt_row = &
               "(1X,I3,':',I3,A1,I2,'/2:',2E12.5)"
        character(len=*), parameter :: fmt_sep = "(1X,65('-'))"

        write(11,fmt_hdr)
        do ni = 1, Ns
            call ReadF(12, ni+4,      P,  Q,  2)
            call ReadF(12, ni+Ns+4,   CP, CQ, 2)

            fx = 0.d0
            do i = 1, ii
                ci   = P(i)*CP(i) + Q(i)*CQ(i)
                C(i) = ci
                c2   = P(i)**2 + Q(i)**2
                if (fx < c2) fx = c2
            end do
            C(ii+4) = P(ii+4)*2 - 1.d0
            call Sint1(x)                          ! integral using file derivatives

            call Dif(P, P1, R, V, P(ii+4), ii, kt, MaxT, h)
            call Dif(Q, Q1, R, V, P(ii+4), ii, kt, MaxT, h)
            do i = 1, ii
                C(i) = P(i)*P1(i) + Q(i)*Q1(i)
            end do
            C(ii+4) = P(ii+4)*2 - 1.d0
            call Sint1(x1)                         ! integral using numerical derivatives

            write(11,fmt_row) ni, Nn(ni), OrbLabels(Ll(ni)+1), Jj(ni), x1/fx, x/fx
        end do
        write(11,fmt_sep)
    end subroutine check_derivative


    subroutine check_unity_expansion
        !! Checks the partial-wave expansion of unity: sum_n |n><n| = 1.
        !!
        !! Identifies all distinct (l, 2j) partial waves present in the basis,
        !! then for each partial wave k accumulates the sum
        !!   P1(r) = sum_{n in pw_k} (int P_n dr) * P_n(r)
        !! which should equal 1 everywhere if the partial wave is complete.
        !! The radial grid is split into 5 equal segments and the average value
        !! of min(P1, 1) in each segment is printed as a percentage.
        !! Results are written to both stdout and unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use sintg, only : Sint1
        implicit none

        integer :: npw, i, n, ln, jn, k, ii5, l_pw, j_pw, nf
        real(dp) :: s, s1, s2, s3, s4, s5, pm1
        real(dp) :: P1(IP6)
        integer  :: lpw(30), jpw(30), kpw(30)
        logical  :: found

        character(len=*), parameter :: fmt_pw = &
               "(/' Expansion of unity',/' Found',i3,' partial waves'," // &
               "/'  p.w.   l   2j    #',/(4I5))"
        character(len=*), parameter :: fmt_r = &
               "(/' p.w. R:',F8.5,5F11.5)"
        character(len=*), parameter :: fmt_pct = &
               "(i3,9x,'<  ',5(f5.1,'%  <  '))"

        ! --- identify distinct (l, 2j) partial waves ---
        npw  = 0
        lpw  = 0
        jpw  = 0
        kpw  = 0
        do n = 1, Ns
            ln    = Ll(n)
            jn    = Jj(n)
            found = .false.
            do k = 1, npw
                if (lpw(k) == ln .and. jpw(k) == jn) then
                    kpw(k) = kpw(k) + 1
                    found  = .true.
                    exit
                end if
            end do
            if (.not. found) then
                npw      = npw + 1
                lpw(npw) = ln
                jpw(npw) = jn
                kpw(npw) = 1
            end if
        end do

        write( *,fmt_pw) npw, (n, lpw(n), jpw(n), kpw(n), n=1,npw)
        write(11,fmt_pw) npw, (n, lpw(n), jpw(n), kpw(n), n=1,npw)

        ii5 = ii / 5

        write( *,fmt_r) R(1), R(ii5), R(2*ii5), R(3*ii5), R(4*ii5), R(ii)
        write(11,fmt_r) R(1), R(ii5), R(2*ii5), R(3*ii5), R(4*ii5), R(ii)

        ! --- for each partial wave, accumulate P1 = sum_n <n|n> P_n ---
        do k = 1, npw
            l_pw = lpw(k)
            j_pw = jpw(k)
            P1   = 0.d0
            nf   = 0
            do n = 1, Ns
                if (Ll(n) == l_pw .and. Jj(n) == j_pw) then
                    nf      = nf + 1
                    call ReadF(12, n+4, P, Q, 2)
                    C(1:ii) = P(1:ii)
                    C(ii+4) = P(ii+4)
                    call Sint1(s)
                    P1(1:ii) = P1(1:ii) + s * P(1:ii)
                end if
            end do
            if (nf /= kpw(k)) then
                write(*,*) ' Error for p.w.', k, ' found ', nf, ' expected ', kpw(k)
                read(*,*)
            end if

            ! average min(P1, 1) over each of the 5 radial segments (as %)
            s1 = 0.d0;  s2 = 0.d0;  s3 = 0.d0;  s4 = 0.d0;  s5 = 0.d0
            do i = 1, ii
                pm1 = min(1.d0, P1(i))
                if      (i <= ii5)   then;  s1 = s1 + 100*pm1/ii5
                else if (i <= 2*ii5) then;  s2 = s2 + 100*pm1/ii5
                else if (i <= 3*ii5) then;  s3 = s3 + 100*pm1/ii5
                else if (i <= 4*ii5) then;  s4 = s4 + 100*pm1/ii5
                else;                       s5 = s5 + 100*pm1/(ii - 4*ii5)
                end if
            end do

            write( *,fmt_pct) k, s1, s2, s3, s4, s5
            write(11,fmt_pct) k, s1, s2, s3, s4, s5
        end do
    end subroutine check_unity_expansion


    subroutine print_orbital(ni)
        !! Prints the radial wavefunction of orbital ni to BAS_INFO.RES.
        !!
        !! Reads the large (P) and small (Q) components from record ni+4 of
        !! the DAT file and writes R(i), P(i), Q(i) for i=1..ii.
        !!
        !! Arguments:
        !!   ni -- orbital index (1-based)
        use readfff, only : ReadF
        implicit none

        integer, intent(in) :: ni
        integer :: i

        character(len=*), parameter :: fmt_hdr = &
               "(1X,65('-'),/20X,'Orbital '," // &
               "I3,' (',I2,A1,I2,'/2):'," // &
               "/4X,'i',8X,'R',12X,'P',12X,'Q')"
        character(len=*), parameter :: fmt_row = "(2X,I3,3E13.5)"

        write(11,fmt_hdr) ni, Nn(ni), OrbLabels(Ll(ni)+1), Jj(ni)
        call ReadF(12, ni+4, P, Q, 2)
        write(11,fmt_row) (i, R(i), P(i), Q(i), i=1,ii)
    end subroutine print_orbital


    subroutine check_orthogonality
        !! Checks orthonormality of the basis set orbitals.
        !!
        !! For each orbital ni, computes the overlap integrals <ni|nk> with all orbitals nk <= ni that share the same (l, 2j) quantum numbers.
        !! Diagonal elements should be 1 (normalisation) and off-diagonal elements should be 0 (orthogonality).
        !! Results are written to unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use sintg, only : Sint1
        implicit none

        integer :: ni, ji, li, n, nk, jk, lk, i
        real(dp) :: s
        real(dp) :: P1(IP6), Q1(IP6)
        real(dp) :: sik(IPs)
        integer  :: nik(IPs)

        character(len=*), parameter :: fmt_hdr = &
               "(1X,65('-'),/15X,'Test of orthonormality'," // &
               "/1X,65('-'))"
        character(len=*), parameter :: fmt_row = &
               "(/4X,'<',I3,'| k > :',/4(I5,': ',F10.7))"

        write(11,fmt_hdr)
        do ni = 1, Ns
            call ReadF(12, ni+4, P, Q, 2)
            ji = Jj(ni)
            li = Ll(ni)
            n = 0
            do nk = 1, ni
                jk = Jj(nk)
                lk = Ll(nk)
                if (li == lk .and. ji == jk) then
                    n = n + 1
                    call ReadF(12, nk+4, P1, Q1, 2)
                    C(1:ii) = P(1:ii)*P1(1:ii) + Q(1:ii)*Q1(1:ii)
                    C(ii+4) = 2*P(ii+4)
                    call Sint1(s)
                    sik(n) = s
                    nik(n) = nk
                end if
            end do
            write(11,fmt_row) ni, (nik(i), sik(i), i=1,n)
        end do
    end subroutine check_orthogonality


    subroutine check_taylor_expansion
        !! Checks the Taylor expansion of orbitals (or their derivatives) at the origin.
        !!
        !! Prompts the user to select orbitals (1) or derivatives (2).
        !! For each orbital ni, reads the wavefunction from the DAT file, prints a header, 
        !! then calls test_taylor_expansion to compare the numerical values against the stored Taylor coefficients, 
        !! and Test_Origin to check the power-law behaviour near r=0. 
        !! Results are written to unit 11 (BAS_INFO.RES).
        use readfff, only : ReadF
        use test_ori, only : Test_Origin
        implicit none

        integer :: iyes, mi, ni, ierr
        character(len=10) :: name

        character(len=*), parameter :: fmt_orb = &
               "(/20X,A10,I4,' (',I2,A1,I2,'/2):'," // &
               "/8X,'R',13X,'PT',13X,'P/PT',10X,'QT',11X,'Q/QT')"

        write(*,*) ' Check Orbitals or derivatives (1,2):'
        read(*,*) iyes
        if ((iyes - 1)*(iyes - 2) /= 0) return

        if (iyes == 1) then
            mi   = 4
            name = ' orbital  '
        else
            mi   = 4 + Ns
            name = 'derivative'
        end if

        do ni = 1, Ns
            call ReadF(12, ni+mi, P, Q, 2)
            write(11,fmt_orb) name, ni, Nn(ni), OrbLabels(Ll(ni)+1), Jj(ni)
            call test_taylor_expansion(P, Q, P(ii+4))
            call Test_Origin(name, P, R, MaxT, ii, 1.d-3, ierr)
        end do
    end subroutine check_taylor_expansion


    subroutine test_taylor_expansion(P_orb, Q_orb, gamma)
        !! Compares numerical orbital values against their Taylor expansion at the origin.
        !!
        !! For each of 2*n1 radial points (n1 points extrapolated before the grid, n1 points on the grid), evaluates the Taylor series:
        !!   P_Taylor(r) = r^gamma * sum_{m=0}^{MaxT} P_orb(ii+5+m) * (r/R(1))^m
        !! and similarly for Q_orb, then prints the ratio numerical/Taylor.
        !! A '>' marker flags the first on-grid point.
        !! Results are written to unit 11 (BAS_INFO.RES).
        !!
        !! Arguments:
        !!   P_orb  -- large radial component; Taylor coefficients at P_orb(ii+5:ii+5+MaxT)
        !!   Q_orb  -- small radial component; Taylor coefficients at Q_orb(ii+5:ii+5+MaxT)
        !!   gamma  -- power-law exponent at the origin (stored in P_orb(ii+4))
        implicit none

        real(dp), intent(in) :: P_orb(IP6), Q_orb(IP6), gamma
        integer  :: n1, n2, n, m, im
        character(len=1) :: marker
        real(dp) :: r1, dr, rn, r_gamma
        real(dp) :: P_taylor, Q_taylor, P_ratio, Q_ratio, x_ratio, x_pow

        character(len=*), parameter :: fmt_row = "(2X,F12.8,2(2X,A1,E14.6,1X,A1,F9.5))"

        r1 = R(1)
        dr = R(2) - r1
        n1 = int(0.2d0 * r1 / dr)   ! number of extrapolated points before grid
        if (n1 < 2) n1 = 2
        n2 = 2 * n1

        do n = 1, n2
            marker = ' '
            if (n == n1 + 1) marker = '>'   ! marks first on-grid point

            if (n <= n1) then
                ! extrapolated point before the grid
                rn      = r1 - (n1 - n + 1) * dr
                P_ratio = 0.d0
                Q_ratio = 0.d0
            else
                ! on-grid point
                rn      = R(n - n1)
                P_ratio = P_orb(n - n1)
                Q_ratio = Q_orb(n - n1)
            end if

            ! evaluate Taylor series at rn
            r_gamma  = rn**gamma
            P_taylor = 0.d0
            Q_taylor = 0.d0
            x_pow    = 1.d0
            x_ratio  = rn / r1
            do m = 0, MaxT
                im       = ii + 5 + m
                P_taylor = P_taylor + r_gamma * x_pow * P_orb(im)
                Q_taylor = Q_taylor + r_gamma * x_pow * Q_orb(im)
                x_pow    = x_pow * x_ratio
            end do

            ! ratio: numerical / Taylor (guard against division by zero)
            P_ratio = P_ratio / (P_taylor + 1.d-39)
            Q_ratio = Q_ratio / (Q_taylor + 1.d-39)

            write(11,fmt_row) rn, marker, P_taylor, marker, P_ratio, marker, Q_taylor, marker, Q_ratio
        end do
    end subroutine test_taylor_expansion

end program bas_info
