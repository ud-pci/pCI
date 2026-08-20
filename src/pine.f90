Program pine
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
    !   Kli=5 - A = E2                Klf=5 -  B = E2
    !   - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   INPUT/OUTPUT files:
    !   CONF.INP    CI input for X1 space
    !   CONF.DAT    radial functions
    !   pCONF.HIJ   Hamiltonian matrix in X1 space (in CSR format, written by pconf v9.0+ with Kw=1)
    !   pCONF.JJJ   J**2 in the CI space of vector X1 (in CSR format, written by pconf v9.0+ with Kw=1)
    !   CONF.XIJ    is used to calculate <n|X1>
    !   DTM.INT     radial integrals
    !   INE.XIJ     solution and R.H.S of the inhom. eq.
    !   INE_J.XIJ   decomposition of the solution of the inhom. eq.
    !   CONF0.XIJ   eigenvectors X0 (initial) and X2 (final)
    !   CONF0.DET   determinants in X0/X2 space
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
    !   Mrec     - file storage units per Real(dp) for DIRECT files
    !   - - - - - - - - - - - - - - - - - - - - - - - - -
    !   main variables in X0 space:
    !   Jm0      - projection of total momentum J
    !   Nc0      - number of configurations
    !   Nd0      - number of dets
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    Use params
    Use conf_variables, Only : nd_start, nd_end, nd_local, HamilCSR_rowptr, JsqCSR_rowptr, Hamil, Jsq, type_real, Crt4
    Use matrix_io, Only : ReadMatrixCSR
    Use conf_init, Only : ReadConfInp, ReadConfigurations
    Use determinants, Only : Dinit, Jterm, Gdet, Rspq, CompCD, Rdet
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime, FormattedMemSize
    Use utils, Only : DetermineRecordLength, CheckRecordLength, CountSubstr, ToUpperString, WriteSkeleton
    Use wigner, Only : FJ3, FJ6
    Use mpi_f08

    Implicit None

    Integer  :: mype, npes, mpierr
    Integer  :: i, n, k, l, Nd0, Nd2, nsign, sgn, nlamb, nrange
    Integer  :: Kli, Klf, Nddir, Nint
    Integer  :: Kdiag, Nlft, Npv, Kt, N0, N2
    Integer  :: N_it4, N_it_ine, Npass
    Integer(kind=int64) :: start_time, t_step
    Real(dp) :: Jm0, E0, E2, Tj0, Tj2, xlamb, xlamb1, xlamb2, xlambstep
    Real(dp) :: XInuc, W0, W0_abs, W00, Q, Elft, Hmin, dlamb, ynorm_y1
    Real(dp), Allocatable, Dimension(:) :: xlamb1s, xlamb2s, xlambsteps
    Integer, Allocatable, Dimension(:) :: RintKey
    Integer, Allocatable, Dimension(:,:) :: Iarr0
    Real(dp), Allocatable, Dimension(:) :: X0, X1, X2, Y1, Y2, RintVal, Ev, dEv, tjv, Diag
    Real(dp), Allocatable, Dimension(:,:) :: X1J, Y2J
    Real(dp), Dimension(2) :: s, ss, s0, s1, s2
    Logical :: converged
    Character(Len=16) :: timeStr
    Character(Len=4)  :: version

    Call MPI_Init(mpierr)
    Call MPI_Comm_rank(MPI_COMM_WORLD, mype, mpierr)
    Call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)
    Call startTimer(start_time)

    ! Initialization
    If (mype == 0) Then
        Open(unit=11, status='UNKNOWN', file='PINE.RES')
        version = '1.0'
        Write( 6,'(4X,"Program pine v",A)') Trim(AdjustL(version))
        Write(11,'(4X,"Program pine v",A)') Trim(AdjustL(version))
        
        Call SetParams           ! Set job and array parameters
        Call Input               ! Read list of configurations from CONF.INP
        Call Init                ! Read basis set from CONF.DAT; reads CONF.DET and CONF0.DET
        Call Rint                ! Read radial integrals from DTM.INT; sorts RintKey/RintVal for bisect
        Call Dinit               ! Form list of determinants
        Call Jterm               ! Print table with numbers of levels with given J
        Call Init0               ! reads X0, X2 from CONF0.XIJ; sets E0/E2/Tj0/Tj2
    End If

    ! Broadcast all parameters and arrays to all ranks
    Call BcastParams(mype)

    ! Print memory requirement estimate
    If (mype == 0) Call calcMemReqs

    ! Read distributed H from pCONF.HIJ
    Call ReadMatrixCSR(Hamil%col, Hamil%val, HamilCSR_rowptr, 'pCONF.HIJ', mype, npes)

    ! Read distributed J^2 from pCONF.JJJ
    ! Pass H's nd_start/nd_end so J^2 uses the same row distribution as H
    Call ReadMatrixCSR(Jsq%col, Jsq%val, JsqCSR_rowptr, 'pCONF.JJJ', mype, npes, nd_start, nd_end)

    nd_local = nd_end - nd_start

    ! Build RHS vector |Y1> = A|X0> and LHS vector <Y2| = <X2|B.
    Call startTimer(t_step)
    Call Vector(mype)
    Call stopTimer(t_step, timeStr)
    If (mype == 0) Then
        Write(*,'(2X,A)') 'TIMING >>> Vector took ' // Trim(timeStr) // ' to complete'
    End If
    ynorm_y1 = dsqrt(sum(Y1**2))

    ! Main loop over wavelengths
    Do k = 1, nrange
        xlamb1    = xlamb1s(k)
        xlamb2    = xlamb2s(k)
        xlambstep = xlambsteps(k)
        If (k > 1) Then
            If (xlamb1 == xlamb2s(k-1)) xlamb1 = xlamb1 + xlambstep
            If (xlamb2 < xlamb1 + xlambstep) xlamb2 = xlamb1
        End If
        If (xlamb2 == xlamb1) xlambstep = 1

        dlamb = Anint((xlamb2 - xlamb1) / xlambstep * 100.0d0) / 100.0d0
        nlamb = Int(dlamb) + 1
        xlamb = xlamb1

        Do n = 1, nlamb
            If (dabs(xlamb) > 1.d-8) Then
                W0_abs = 1.d+7 / (xlamb * 219474.63d0)
                nsign = 2
            Else
                W0_abs = 0.0d0
                nsign = 1
            End If

            ! Load CONF.XIJ once per lambda; Ev(k) = Elft + dEv(k) computed per sign below.
            If (.not. Allocated(X1J)) Allocate(X1J(Nd, Npv))
            If (.not. Allocated(Ev))  Allocate(Ev(Npv))
            If (.not. Allocated(dEv)) Allocate(dEv(Npv))
            If (.not. Allocated(tjv)) Allocate(tjv(Npv))
            Call FormV
            If (mype == 0) Then
                Write(*,'(3X,"Number of vectors: ",I2)') Nlft
                Write(11,'(3X,"Number of vectors: ",I2)') Nlft
            End If

            Do i = 1, nsign
                sgn  = 3 - 2*i        ! i=1: sgn=+1, i=2: sgn=-1
                W0   = sgn * W0_abs
                Elft = E0 + W0
                If (Kli == 2 .AND. Klf /= 2) Elft = E2
                Ev(4:Nlft) = Elft + dEv(4:Nlft)

                If (mype == 0) Then
                    If (dabs(xlamb) > 1.d-8) Then
                        Write(*,'(/3X,34("-")/3X,"Calculation for lambda=",F11.3,/3X,34("-")/)') sgn * xlamb
                        Write(11,'(/3X,34("-")/3X,"Calculation for lambda=",F11.3,/3X,34("-")/)') sgn * xlamb
                    Else
                        Write(*,'(/3X,25("-")/3X,"Calculation for omega = 0",/3X,25("-")/)')
                        Write(11,'(/3X,25("-")/3X,"Calculation for omega = 0",/3X,25("-")/)')
                    End If
                End If
                Call startTimer(t_step)
                Call SolEq1(i)
                Call stopTimer(t_step, timeStr)
                If (mype == 0) Then
                    Write(*,'(4X,"SolEq1: ",A)') Trim(timeStr)
                    Write(11,'(4X,"SolEq1: ",A)') Trim(timeStr)
                End If

                Do l = 1, Npass
                    If (mype == 0) Then
                        Write(*,'(/4X,"--- Pass",I0," ---")') l
                        Write(11,'(/4X,"--- Pass",I0," ---")') l
                    End If
                    Call startTimer(t_step)
                    Call SolEq4(converged, i, Nlft, ynorm_y1)
                    Call stopTimer(t_step, timeStr)
                    If (mype == 0) Then
                        Write(*,'(4X,"SolEq4: ",A)') Trim(timeStr)
                        Write(11,'(4X,"SolEq4: ",A)') Trim(timeStr)
                    End If
                    If (converged) Exit
                End Do

                If (mype == 0 .AND. .NOT. converged) Then
                    Write(*,'(4X,"Convergence was not reached.")')
                    Write(11,'(4X,"Convergence was not reached.")')
                End If

                ! Prj and Prin need all ranks (J^2 SpMV via MxmpyJ).
                If (Kli /= 5) Then
                    Call Prj('  X1  ', Tj0, X1, X1J)
                    If (mype == 0) Call RdcX1J
                    Call Prj('  Y2  ', Tj2, Y2, Y2J)
                Else
                    Call PrjE2('  X1  ', Tj0, X1, X1J)
                    If (mype == 0) Call RdcX1J
                    Call PrjE2('  Y2  ', Tj2, Y2, Y2J)
                End If
                Call PrintDiag

                If (mype == 0) Then
                    If (Kli==2.AND.Klf==2) Call RdcE1(i)
                    If (Kli==5)            Call RdcE2(i)
                    If (Kli==1.OR.Klf==1) Call RdcPNC
                    If (Klf==3)            Call RdcAM
                    If (Kli==2.AND.Klf==2) Call C_3
                End If
            End Do

            xlamb = xlamb + xlambstep
        End Do
    End Do

    Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    Call stopTimer(start_time, timeStr)
    If (mype == 0) Write(*,'(2X,A)') 'TIMING >>> pine total: ' // Trim(timeStr)
    Call MPI_Finalize(mpierr)

Contains

    Subroutine calcMemReqs
        ! Print a pre-calculation memory estimate.
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer(kind=int64) :: nnz_H, nnz_J
        Integer(kind=int64) :: mem_iarr, mem_iarr0, mem_x0x2
        Integer(kind=int64) :: mem_yy, mem_soln, mem_x1j, mem_y2j
        Integer(kind=int64) :: mem_soleq4_ws, mem_H_csr, mem_J_csr
        Integer(kind=int64) :: mem_soleq1_az, mem_soleq1_buf, mem_soleq1_peak_r0
        Integer(kind=int64) :: mem_replicated, mem_peak
        Integer(kind=int64) :: nnz_H_rank, nnz_J_rank, jend_approx
        Integer(kind=int64) :: bI, bR, bDP
        Integer :: Ntr, ic, id
        Character(Len=16) :: ms

        Call PeekCSRNnz('pCONF.HIJ', nnz_H)
        Call PeekCSRNnz('pCONF.JJJ', nnz_J)

        bI  = 4_int64               ! int32 bytes
        bR  = int(type_real, int64) ! val element size (4 or 8)
        bDP = 8_int64               ! dp bytes

        ! CSR per rank (ceiling division, no row_ptr overhead included separately)
        nnz_H_rank = (nnz_H + npes - 1) / npes
        nnz_J_rank = (nnz_J + npes - 1) / npes
        mem_H_csr  = nnz_H_rank * (bI + bR)   ! col(int32) + val(type_real)
        mem_J_csr  = nnz_J_rank * (bI + bR)

        ! Replicated arrays (same size on every rank):
        mem_iarr      = int(Ne, int64) * int(Nd, int64) * bI      ! Iarr
        mem_iarr0     = int(Ne, int64) * int(Nd0, int64) * bI  ! Iarr0
        mem_x0x2      = 2_int64 * int(Nd0, int64) * bDP        ! X0, X2
        mem_yy        = 2_int64 * int(Nd, int64) * bDP             ! Y1, Y2
        mem_soln      = 2_int64 * int(Nd, int64) * bDP             ! X1, Diag
        mem_x1j       = int(Npv, int64) * int(Nd, int64) * bDP    ! X1J(Nd,Npv)
        mem_y2j       = 5_int64 * int(Nd, int64) * bDP             ! Y2J(Nd,5)
        mem_soleq4_ws = 3_int64 * int(Nd, int64) * bDP             ! XX3,X12,X13 (peak, freed after)

        mem_replicated = mem_iarr + mem_iarr0 + mem_x0x2 + mem_yy + mem_soln + mem_x1j + mem_y2j
        mem_peak = mem_replicated + mem_soleq4_ws + mem_H_csr + mem_J_csr

        ! SolEq1 cold-start peak (rank 0 only; freed before SolEq4 runs):
        !   Az(Ntr,Ntr) dense matrix + col_buf + val_buf for first Ntr rows of pCONF.HIJ
        Ntr = 0
        Do ic = 1, Nc
            id = Ndc(ic)
            If (Ntr + id <= Nddir) Ntr = Ntr + id
        End Do
        If (Ntr == 0) Ntr = min(Nddir, Nd)
        jend_approx    = (nnz_H * Ntr) / max(Nd, 1)        ! approximate nnz in first Ntr rows
        mem_soleq1_az  = int(Ntr, int64)**2 * bDP          ! Az dense block
        mem_soleq1_buf = jend_approx * (bI + bR)               ! col_buf + val_buf
        ! Rank 0 peak = base replicated + both CSR matrices + SolEq1 temporaries
        mem_soleq1_peak_r0 = mem_replicated + mem_H_csr + mem_J_csr &
                           + mem_soleq1_az + mem_soleq1_buf

        Call FormattedMemSize(max(mem_peak, mem_soleq1_peak_r0), ms)
        Write(*,'(A,A,A)') 'calcMemReqs: pine will require at least ',Trim(ms),' per rank'
    End Subroutine calcMemReqs

    Subroutine PeekCSRNnz(filename, nnz)
        ! Read only the total_nnz field from a global CSR file header (bytes 9-16).
        Character(Len=*),    Intent(In)  :: filename
        Integer(kind=int64), Intent(Out) :: nnz
        Integer :: err_stat
        Open(unit=17, file=filename, status='OLD', access='STREAM', form='UNFORMATTED', iostat=err_stat)
        If (err_stat /= 0) Then
            Write(*,'(2X,"PeekCSRNnz: cannot open ",A)') Trim(filename)
            nnz = -1_int64
            Return
        End If
        Read(17, pos=9) nnz
        Close(17)
    End Subroutine PeekCSRNnz

    Subroutine BcastParams(mype)
        ! Broadcast all rank-0 arrays and scalars to remaining ranks.
        ! Jm0 must be broadcast as MPI_DOUBLE_PRECISION, not MPI_INTEGER.
        Use mpi_f08
        Use mpi_utils, Only : BroadcastI
        Implicit None
        Integer, Intent(In) :: mype
        Integer(kind=int64) :: count

        Call MPI_Bcast(Nso,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jm,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jm0,     1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd0,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Kli,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Klf,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ns,      1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nint,    1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(E0,      1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(E2,      1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tj0,     1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Tj2,     1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Crt4,    1,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Npv,     1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlft,    1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nddir,   1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it4,   1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(N_it_ine,1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Npass,   1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nrange,  1,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)

        If (mype /= 0) Then
            If (.not. Allocated(Nh))      Allocate(Nh(Nst))
            If (.not. Allocated(Jz))      Allocate(Jz(Nst))
            If (.not. Allocated(Nvc))     Allocate(Nvc(Nc))
            If (.not. Allocated(Nc0))     Allocate(Nc0(Nc))
            If (.not. Allocated(Ndc))     Allocate(Ndc(Nc))
            If (.not. Allocated(RintVal)) Allocate(RintVal(Nint))
            If (.not. Allocated(RintKey)) Allocate(RintKey(Nint))
            If (.not. Allocated(X0))      Allocate(X0(Nd0))
            If (.not. Allocated(X2))      Allocate(X2(Nd0))
        End If

        Call MPI_Bcast(Nh,      Nst,  MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz,      Nst,  MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nn,      Ns,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj,      Ns,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ll,      Ns,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nvc,     Nc,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc0,     Nc,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ndc,     Nc,   MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(RintVal, Nint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(RintKey, Nint, MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(X0,      Nd0,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(X2,      Nd0,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        ! Broadcast xlamb range arrays (allocated by Input on rank 0)
        If (mype /= 0) Then
            If (.not. Allocated(xlamb1s))    Allocate(xlamb1s(nrange))
            If (.not. Allocated(xlamb2s))    Allocate(xlamb2s(nrange))
            If (.not. Allocated(xlambsteps)) Allocate(xlambsteps(nrange))
        End If
        Call MPI_Bcast(xlamb1s,    nrange, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(xlamb2s,    nrange, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(xlambsteps, nrange, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        ! Broadcast Iarr (Ne*Nd determinant list) and Iarr0 (Ne*Nd0)
        count = Int(Ne, kind=int64) * Int(Nd, kind=int64)
        If (mype /= 0 .and. .not. Allocated(Iarr)) Allocate(Iarr(Ne, Nd))
        Call BroadcastI(Iarr(1,1), count, 0, 0, MPI_COMM_WORLD, mpierr)
        count = Int(Ne, kind=int64) * Int(Nd0, kind=int64)
        If (mype /= 0 .and. .not. Allocated(Iarr0)) Allocate(Iarr0(Ne, Nd0))
        Call BroadcastI(Iarr0(1,1), count, 0, 0, MPI_COMM_WORLD, mpierr)

    End Subroutine BcastParams

    Subroutine Mxmpy(x, y)
        ! Distributed CSR SpMV: y = (Elft - H) * x
        !
        ! Each rank applies its local CSR rows (nd_start+1 ... nd_end), then
        ! MPI_Allreduce(SUM) assembles the full Nd-element y on every rank.
        ! x must be the full Nd-element vector (replicated on all ranks).
        ! The diagonal shift is: H_diag(n) stored as H(n,n), subtracted from Elft:
        !   y(n) += (Elft - H(n,n)) * x(n)  when n == k (diagonal)
        !   y(n) -= H(n,k) * x(k)            for off-diagonal n > k (lower triangle)
        !   y(k) -= H(n,k) * x(n)            symmetric contribution
        Use mpi_f08
        Implicit None
        Real(dp), Intent(In),  Dimension(Nd) :: x
        Real(dp), Intent(Out), Dimension(Nd) :: y

        Integer :: r, n, k
        Integer(kind=int64) :: j8
        Real(dp) :: t

        y = 0.0d0

        Do r = 1, nd_local
            n = nd_start + r
            Do j8 = HamilCSR_rowptr(r-1)+1_int64, HamilCSR_rowptr(r)
                k  = Hamil%col(j8)
                t  = Real(Hamil%val(j8), kind=dp)
                If (n == k) t = t - Elft       ! diagonal: t = H(n,n) - Elft
                y(n) = y(n) - t * x(k)
                If (n /= k) y(k) = y(k) - t * x(n)
                If (Kdiag == 1 .AND. n == k) Diag(n) = -t  ! Diag(n) = Elft - H(n,n)
            End Do
        End Do

        Call MPI_Allreduce(MPI_IN_PLACE, y, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        If (Kdiag == 1) Call MPI_Allreduce(MPI_IN_PLACE, Diag, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Kdiag = 0
    End Subroutine Mxmpy

    Subroutine MxmpyJ(x, y)
        ! Distributed CSR SpMV: y = J^2 * x
        ! Same pattern as Mxmpy; used by Prj, PrjE2, Prin.
        Use mpi_f08
        Implicit None
        Real(dp), Intent(In),  Dimension(Nd) :: x
        Real(dp), Intent(Out), Dimension(Nd) :: y
        Integer :: r, n, k
        Integer(kind=int64) :: j8
        Real(dp) :: t

        y = 0.0d0
        Do r = 1, nd_local
            n = nd_start + r
            Do j8 = JsqCSR_rowptr(r-1)+1_int64, JsqCSR_rowptr(r)
                k  = Jsq%col(j8)
                t  = Real(Jsq%val(j8), kind=dp)
                y(n) = y(n) + t * x(k)
                If (n /= k) y(k) = y(k) + t * x(n)
            End Do
        End Do
        Call MPI_Allreduce(MPI_IN_PLACE, y, Nd, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, mpierr)
    End Subroutine MxmpyJ

    Subroutine SolEq1(sign)
        ! Builds initial X1 for SolEq4.
        ! Reads the leading Ntr rows from pCONF.HIJ, builds the Ntr*Ntr dense matrix (Elft-H), 
        ! solves via DSYSV, and sets X1(1:Ntr).
        Use mpi_f08
        Implicit None
        External DSYSV
        Integer, Intent(In) :: sign

        Integer :: Ntr, ic, id, info, lwork, err_stat
        Integer :: i, j, j32
        Integer(kind=int64) :: Nd_g, total_nnz, jend, pos_col, pos_val, j8
        Integer(kind=int64), Allocatable, Dimension(:)   :: rp
        Integer,             Allocatable, Dimension(:)   :: col_buf, ipiv
        Real(type_real),     Allocatable, Dimension(:)   :: val_buf
        Real(dp),            Allocatable, Dimension(:)   :: rhs_v, wrk
        Real(dp),            Allocatable, Dimension(:,:) :: Az
        Real(dp), Dimension(1) :: dwork
        Character(Len=256) :: err_msg

        ! Round Ntr down to a configuration boundary.
        Ntr = 0
        Do ic = 1, Nc
            id = Ndc(ic)
            If (Ntr + id <= Nddir) Ntr = Ntr + id
        End Do
        If (Ntr == 0) Ntr = min(Nddir, Nd)

        If (.NOT. Allocated(X1)) Allocate(X1(Nd))
        X1 = 0.0d0

        ! --- Direct solve on the leading Ntr*Ntr block ---
        If (mype == 0) Then
            Write(*,'(4X,"SolEq1: matrix (",I5," x ",I5,"), Nd =",I8)') &
                Ntr, Ntr, Nd

            ! Read first Ntr rows from pCONF.HIJ (STREAM, 1-based byte positions):
            !   pos=1 : Nd_g      (int64)
            !   pos=9 : total_nnz (int64)
            !   pos=17: row_ptr(0:Nd_g) -- Ntr+1 entries are sufficient
            !   pos_col = 17 + (Nd_g+1)*8 : int32 col entries (upper-triangle: col >= row)
            !   pos_val = pos_col + total_nnz*4 : type_real val entries
            Allocate(rp(0:Ntr))
            Open(unit=17, file='pCONF.HIJ', status='OLD', access='STREAM', &
                 form='UNFORMATTED', iostat=err_stat, iomsg=err_msg)
            If (err_stat /= 0) Then
                Write(*,'(4X,"SolEq1: cannot open pCONF.HIJ - ",A)') Trim(err_msg)
                Deallocate(rp)
                Return
            End If
            Read(17, pos=1_int64)  Nd_g
            Read(17, pos=9_int64)  total_nnz
            Read(17, pos=17_int64) rp(0:Ntr)
            jend    = rp(Ntr)
            pos_col = 17_int64 + (Nd_g + 1_int64) * 8_int64
            pos_val = pos_col  + total_nnz * 4_int64

            Allocate(col_buf(jend), val_buf(jend))
            Read(17, pos=pos_col) col_buf(1:jend)
            Read(17, pos=pos_val) val_buf(1:jend)
            Close(17)

            ! Build Az = Elft*I - H(1:Ntr, 1:Ntr).
            ! pCONF.HIJ uses upper-triangle storage (col > row), so for row i <= Ntr the column j can reach up to Nd >> Ntr.
            ! Guard against out-of-block elements: only fill Az for j <= Ntr (the truncated subblock).
            Allocate(Az(Ntr, Ntr))
            Az = 0.0d0
            Do i = 1, Ntr
                Do j8 = rp(i-1) + 1_int64, rp(i)
                    j32 = col_buf(j8)
                    If (j32 > Ntr) Cycle   ! element crosses block boundary
                    If (i == j32) Then
                        Az(i,i) = Elft - Real(val_buf(j8), kind=dp)
                    Else
                        Az(i,j32) = -Real(val_buf(j8), kind=dp)
                        Az(j32,i) = -Real(val_buf(j8), kind=dp)
                    End If
                End Do
            End Do
            Deallocate(rp, col_buf, val_buf)

            ! Solve Az * rhs_v = Y1(1:Ntr) via DSYSV (symmetric indefinite)
            Allocate(rhs_v(Ntr), ipiv(Ntr))
            rhs_v(1:Ntr) = Y1(1:Ntr)
            Call DSYSV('U', Ntr, 1, Az, Ntr, ipiv, rhs_v, Ntr, dwork, -1, info)
            lwork = Int(dwork(1))
            Allocate(wrk(lwork))
            Call DSYSV('U', Ntr, 1, Az, Ntr, ipiv, rhs_v, Ntr, wrk, lwork, info)
            Deallocate(Az, wrk, ipiv)
            If (info /= 0) Then
                Write(*,'(4X,"SolEq1: DSYSV info=",I0,"; using X1=0")') info
                rhs_v = 0.0d0
            End If
            X1 = 0.0d0
            X1(1:Ntr) = rhs_v(1:Ntr)
            Deallocate(rhs_v)
        End If
        Call MPI_Bcast(X1, Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    End Subroutine SolEq1

    Subroutine SolEq4(converged, sign, num, ynorm_in)
        ! Probe-vector subspace minimisation: iterative solution of (Elft - H)|X1> = |Y1>.
        ! MPI notes:
        !   - Mxmpy (3x per iter) does internal Allreduce; Y2J fully replicated after.
        !   - Scalar accumulations (residual, Vr/Z1) distributed over nd_start...nd_end;
        !     results Allreduced to small arrays. Array updates (X1J, XX3, X1) replicated.
        !   - Decomp/Flsolv: serial on Npv-sized system, same result on all ranks.
        Use mpi_f08
        Use solvers
        Implicit None
        Integer,  Intent(In)  :: sign, num
        Real(dp), Intent(In)  :: ynorm_in
        Logical,  Intent(Out) :: converged

        Integer :: i, itr, k, l, n
        Real(dp) :: ynorm, crit, dnorm, s2, s3, s12, s13, dx, di, x, y, dx1, dnorm_check
        Real(dp), Dimension(5) :: buf
        Character(Len=9) :: str
        Character(Len=256) :: strfmt
        Real(dp), Allocatable, Dimension(:) :: XX3, X12, X13
        Real(dp), Allocatable, Dimension(:) :: Vl, Vr, scales, Z1
        Integer,  Allocatable, Dimension(:) :: Mps

        If (.not. Allocated(Diag))   Allocate(Diag(Nd))
        If (.not. Allocated(X1))     Allocate(X1(Nd))
        If (.not. Allocated(X1J))    Allocate(X1J(Nd, Npv))
        If (.not. Allocated(Y2J))    Allocate(Y2J(Nd, 5))
        If (.not. Allocated(Ev))     Allocate(Ev(Npv))
        If (.not. Allocated(XX3))    Allocate(XX3(Nd))
        If (.not. Allocated(X12))    Allocate(X12(Nd))
        If (.not. Allocated(X13))    Allocate(X13(Nd))
        If (.not. Allocated(Vl))     Allocate(Vl(Npv))
        If (.not. Allocated(Vr))     Allocate(Vr(Npv))
        If (.not. Allocated(scales)) Allocate(scales(Npv))
        If (.not. Allocated(Mps))    Allocate(Mps(Npv))
        If (.not. Allocated(Z1))     Allocate(Z1(Npv*Npv))

        ! --- Initialise ---
        ynorm = ynorm_in
        crit  = Crt4 * ynorm
        Do i = 1, Nd
            X1J(i,1) = X1(i)
            XX3(i)   = X1(i)
        End Do
        If (mype == 0) Then
            strfmt = '(3X,63("="),/4X,"SolEq4: ||Y1|| =",E10.3," Crit =",E10.3,/3X,63("="))'
            Write( *,strfmt) ynorm, crit
            Write(11,strfmt) ynorm, crit
        End If

        Diag       = 0.0d0   ! must be zero before Mxmpy fills local rows and Allreduces
        Kdiag      = 1
        dnorm_check = huge(1.0d0)
        converged   = .FALSE.
        str         = ' DIVERGED'

        Do itr = 1, N_it_ine
            ! Y2J(:,1) = (Elft - H) * X1; fills Diag on first call (Kdiag=1)
            Call Mxmpy(X1, Y2J(:,1))

            ! Residual scalars: distributed accumulation over local rows, then Allreduce.
            dnorm = 0.0d0
            s2    = 0.0d0
            s3    = 0.0d0
            s12   = 0.0d0
            s13   = 0.0d0
            Do i = nd_start+1, nd_end
                dx = Y1(i) - Y2J(i,1)
                dnorm = dnorm + dx*dx
                di = Diag(i)
                If (dabs(di) > 1.d-6) Then
                    dx1 = dx/di
                Else
                    dx1 = 0.0d0
                End If
                s2  = s2  + X1J(i,1)*X1J(i,1)
                s3  = s3  + XX3(i)*XX3(i)
                s12 = s12 + X1J(i,1)*dx
                s13 = s13 + XX3(i)*dx1
            End Do
            buf = [dnorm, s2, s3, s12, s13]
            Call MPI_Allreduce(MPI_IN_PLACE, buf, 5, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            dnorm = buf(1)
            s2    = buf(2)
            s3    = buf(3)
            s12   = buf(4)
            s13   = buf(5)
            If (s2 > 1.d-12) Then
                s12 = s12/s2
            Else
                s12 = 0.0d0
            End If
            If (s3 > 1.d-12) Then
                s13 = s13/s3
            Else
                s13 = 0.0d0
            End If
            ! Array fill and update: replicated (X1J/XX3 must be full on all ranks).
            s2 = 0.0d0
            s3 = 0.0d0
            Do i = 1, Nd
                dx = Y1(i) - Y2J(i,1)
                di = Diag(i)
                If (dabs(di) > 1.d-6) Then
                    dx1 = dx/di
                Else
                    dx1 = 0.0d0
                End If
                x = dx  - s12*X1J(i,1)
                y = dx1 - s13*XX3(i)
                s2 = s2 + x*x
                s3 = s3 + y*y
                X1J(i,1) = x
                X1J(i,2) = x
                X1J(i,3) = y
                XX3(i)   = y
            End Do
            dnorm = dsqrt(dnorm)
            s2    = dsqrt(s2)
            s3    = dsqrt(s3)

            If (mype == 0) Then
                strfmt = '(3X,73("-"),/4X,"iter=",I3, &
                      /4X,"Norm of probe vectors:",2E12.3, &
                      /4X,"Norm of the residual vector:",E12.3)'
                Write( *,strfmt) itr, s2, s3, dnorm
                Write(11,strfmt) itr, s2, s3, dnorm
            End If

            If (dnorm <= crit) Then
                str = 'converged'
                converged = .TRUE.
            End If

            ! Stall detection: every 20 iterations check improvement over that window.
            ! < 1% reduction means the subspace is stuck; exit and let the next pass
            ! warm-restart with a fresh subspace.
            If (itr == 1) dnorm_check = dnorm
            If (itr > 1 .AND. Mod(itr, 20) == 0) Then
                If (dnorm >= 0.99d0 * dnorm_check) Then
                    If (mype == 0) Then
                        strfmt = '(4X,"Stall: <1% reduction over 20 iters; restarting pass.")'
                        Write( *,strfmt)
                        Write(11,strfmt)
                    End If
                    Exit
                End If
                dnorm_check = dnorm
            End If

            ! (Elft-H) applied to probe vectors 2 and 3
            X12 = X1J(:,2)
            X13 = X1J(:,3)
            Call Mxmpy(X12, Y2J(:,2))
            Call Mxmpy(X13, Y2J(:,3))

            ! Build overlap matrix Z1 and RHS Vr: distributed over local rows.
            Vr(1:num) = 0.0d0
            Z1(1:num*num) = 0.0d0
            Do i = nd_start+1, nd_end
                y = Y1(i)
                Do k = 1, num
                    If (k <= 3) Then
                        Vl(k) = Y2J(i,k)
                    Else
                        Vl(k) = X1J(i,k) * Ev(k)   ! eigenvec approx: (Elft-H)*v ~= Ev(k)*v
                    End If
                    Vr(k) = Vr(k) + Vl(k)*y
                End Do
                Do k = 1, num
                    Do l = 1, num
                        n = (k-1)*num + l
                        Z1(n) = Z1(n) + Vl(k)*Vl(l)
                    End Do
                End Do
            End Do
            Call MPI_Allreduce(MPI_IN_PLACE, Vr, num,     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_Allreduce(MPI_IN_PLACE, Z1, num*num, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

            ! LU solve: Vl = Z1^{-1} * Vr  (serial, Npv-sized system)
            Call Decomp(Z1, num, scales, Mps)
            Call Flsolv(num, Z1, Vr, Vl, Mps)

            ! Update X1 as linear combination of probe vectors
            Do i = 1, Nd
                X1(i) = Vl(1)*X1(i)
                Do k = 2, num
                    X1(i) = X1(i) + Vl(k)*X1J(i,k)
                End Do
            End Do

            If (mype == 0) Then
                strfmt = '(4X,"mixing coefficients",/(4x,8F9.4))'
                Write( *,strfmt) (Vl(k), k=1,num)
                Write(11,strfmt) (Vl(k), k=1,num)
            End If

            If (converged) Exit
        End Do

        If (mype == 0) Then
            strfmt = '(4X,"Iteration process ",A9,/3X,63("="))'
            Write( *,strfmt) str
            Write(11,strfmt) str
        End If

        Deallocate(XX3, X12, X13, Vl, Vr, scales, Mps, Z1)
    End Subroutine SolEq4

    Subroutine FormV
        ! Read CONF.XIJ and load eigenvectors with |J - Tj0| < dj
        ! as additional probe vectors into X1J(:,4...Nlft).  Sets Nlft and dEv on exit.
        ! dEv(k) = raw eigenvalue offset from file; caller sets Ev(k) = Elft + dEv(k) per sign.
        ! All ranks read the file independently (shared filesystem; avoids broadcasting Nd*Nlv doubles).
        Implicit None
        Integer :: i, num, ndd, err_stat
        Real(dp) :: dj, tj, dE
        Character(Len=256) :: strfmt

        num = 4    ! columns 1-3 reserved as solver workspace
        Open(unit=16, file='CONF.XIJ', status='OLD', form='UNFORMATTED', iostat=err_stat)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.XIJ is absent"/)'
            Write(*,strfmt)
            Write(11,strfmt)
            Stop
        End If
        If (Kli == 1) Then
            dj = 0.1d0
        Else
            dj = 1.1d0
        End If
        Do While (num <= Npv)
            Read(16, iostat=err_stat) dE, tj, ndd, (X1J(i,num), i=1,Nd)
            If (err_stat /= 0) Exit
            If (tj > Tj0-dj .AND. tj < Tj0+dj) Then
                dEv(num) = dE              ! dE = -|eigenvalue| from file; Ev(k) = Elft + dEv(k)
                tjv(num) = tj
                If (mype == 0) Then
                    strfmt = '(4X,"Vector ",I3," E =",F10.5," J =",F5.2)'
                    Write(*,strfmt) num, dE, tj
                End If
                num = num+1
            End If
        End Do
        Nlft = num-1
        Close(16)
    End Subroutine FormV

    Subroutine Vector(mype)
        ! Build RHS vector |Y1> = A|X0> and LHS vector <Y2| = <X2|B by looping over all pairs of determinants 
        ! (n in X0/X2 space, k in X1 space) and accumulating one-particle matrix elements via Hint.
        Use conf_variables, Only : iconf1, iconf2
        Use determinants, Only : Gdet, Rspq, CompCD, Rdet
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: mype

        Integer :: n, k, ic, icomp, is, i0, i1, j0, j1, negl, idif, kx
        Integer :: mdel, nskip, k1, nf, iq, rng_start, rng_end
        Integer, Allocatable, Dimension(:) :: idet0, idet1
        Real(dp) :: trd, xn0, xn2
        Character(Len=256) :: strfmt, err_msg

        trd  = 1.d-10
        idif = 1
        If (Kli == 5) idif = 2

        If (.not. Allocated(idet0))  Allocate(idet0(Ne))
        If (.not. Allocated(idet1))  Allocate(idet1(Ne))
        If (.not. Allocated(iconf1)) Allocate(iconf1(Ne))
        If (.not. Allocated(iconf2)) Allocate(iconf2(Ne))
        If (.not. Allocated(Y1))    Allocate(Y1(Nd))
        If (.not. Allocated(Y2))    Allocate(Y2(Nd))

        Y1 = 0.0d0
        Y2 = 0.0d0
        negl = 0

        If (mype == 0) Then
            idet0(1:Ne) = Iarr0(1:Ne, 1)
            Call DefJz2(idet0)
            Q = Jm - Jm0
            mdel  = Int(dabs(Q) + 1.d-5)
            If ((Kli-1)*(Klf-1) == 0 .AND. mdel /= 0) Then
                Write(*,*) 'Vector: |M-M0| is not zero!'
                Stop
            End If
        End If
        Call MPI_Bcast(Q, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        ! Divide Nd0 outer loop across ranks
        rng_start = mype * (Nd0 / npes) + 1
        rng_end   = (mype+1) * (Nd0 / npes)
        If (mype == npes-1) rng_end = Nd0

        Do n = rng_start, rng_end
            xn0 = X0(n)
            xn2 = X2(n)
            idet0(1:Ne) = Iarr0(1:Ne, n)
            If ((dabs(xn0) + dabs(xn2)) > trd + trd) Then
                nskip = 0
                k     = 0
                Do ic = 1, Nc
                    kx = Ndc(ic)
                    Call Gdet(k+1, idet1)
                    Call CompCD(idet0, idet1, icomp)
                    If (icomp > 1 .OR. Jdel > idif) Then
                        k      = k + kx
                        nskip  = nskip + kx
                    Else
                        Do k1 = 1, kx
                            k = k + 1
                            Call Gdet(k, idet1)
                            Call Rspq(idet0, idet1, is, nf, i0, i1, j0, j1)
                            If (nf == 1) Then
                                If (dabs(xn0) > trd) Y1(k) = Y1(k) + is * Hint(j1, j0, Kli) * xn0
                                If (dabs(xn2) > trd) Y2(k) = Y2(k) + is * Hint(j0, j1, Klf) * xn2
                            End If
                            If (Kli == 5 .AND. nf == 0) Then
                                Do iq = 1, Ne
                                    i1 = idet1(iq)
                                    If (dabs(xn0) > trd) Y1(k) = Y1(k) + is * Hint(i1, i1, Kli) * xn0
                                    If (dabs(xn2) > trd) Y2(k) = Y2(k) + is * Hint(i1, i1, Klf) * xn2
                                End Do
                            End If
                        End Do
                    End If
                End Do
            Else
                negl = negl + 1
            End If
        End Do

        Call MPI_Allreduce(MPI_IN_PLACE, Y1, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_Allreduce(MPI_IN_PLACE, Y2, Nd, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

    End Subroutine Vector

    Real(dp) Function Hint(i1, i0, kll)
        ! <i1|A|i0>
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i1, i0, is, k, kll, ki, j0, l0, m0, kf, j1, l1, m1, jdel, mdel, lx
        Real(dp) :: xj0, xm0, xj1, xl0, xm1, q_loc, xl1, h, a1, w3, w6
        ki = Nh(i0)
        j0 = Jj(ki)
        l0 = Ll(ki)
        m0 = Jz(i0)
        kf = Nh(i1)
        j1 = Jj(kf)
        l1 = Ll(kf)
        m1 = Jz(i1)
        xj0 = j0/2.d0
        xm0 = m0/2.d0
        xj1 = j1/2.d0
        xl0 = l0
        xm1 = m1/2.d0
        q_loc = xm1-xm0
        xl1 = l1
        jdel = iabs(j1-j0)/2
        mdel = iabs(m1-m0)/2
        lx = (l1+l0+1)/2
        h = 0.d0
        If (kll == 1) Then
            If (jdel == 0 .AND. mdel == 0) h = -Fint(5, kf, ki, -1)
        End If
        If (kll <= 4) Then
            If (jdel <= 1 .AND. mdel <= 1) Then
                w3 = FJ3(xj1, 1.d0, xj0, -xm1, q_loc, xm0)
                If (kll == 2 .OR. kll == 4) Then
                    k = (j1+m1+j0+1)/2 + lx
                    is = 1
                    If (k /= 2*(k/2)) is = -1
                    a1 = is * dsqrt(Real((j1+1)*(j0+1)*lx, dp))
                    w6 = FJ6(xl1, xj1, 0.5d0, xj0, xl0, 1.0d0)
                    If (kll == 2) Then
                        h = a1*w3*w6*Fint(3, kf, ki, +1)
                    Else
                        h = a1*w3*w6*Fint(6, kf, ki, -1)
                    End If
                End If
                If (kll == 3) Then
                    k = (j1+j0-m1+1)/2 + lx
                    is = 1
                    If (k /= 2*(k/2)) is = -1
                    h = is*w3*Fint(7, kf, ki, -1)
                End If
            End If
        End If
        If (kll == 5) Then
            If (jdel <= 2 .AND. mdel <= 2 .AND. l0+l1 == 2*((l0+l1)/2)) Then
                If (j1+j0 >= 4) Then
                    w3 = FJ3(xj1, 2.d0, xj0, -xm1, q_loc, xm0)
                    k = (j1-m1+j0+1)/2
                    is = 1
                    If (k /= 2*(k/2)) is = -1
                    a1 = is*dsqrt(Real((j1+1)*(j0+1)*(2*l1+1)*(2*l0+1), dp))
                    w3 = w3*FJ3(xl1, 2.d0, xl0, 0.d0, 0.d0, 0.d0)
                    w6 = FJ6(xl0, xj0, 0.5d0, xj1, xl1, 2.d0)
                    h  = a1*w3*w6*Fint(10, kf, ki, +1)
                End If
            End If
        End If
        Hint = h
    End Function Hint

    Real(dp) Function Fint(is, nfin, nini, ic)
        ! Radial integral lookup via binary search in the sorted RintKey array.
        Implicit None
        Integer, Intent(In) :: is, nfin, nini, ic
        Integer :: isg, na, nb, ind, lo, hi, mid
        Character(Len=128) :: strfmt
        isg = 1
        na = nfin
        nb = nini
        If (na > nb) Then
            na = nini
            nb = nfin
            isg = ic
        End If
        ind = is*IPx*IPx + (na-Nso)*IPx + (nb-Nso)
        lo = 1
        hi = Nint
        Do While (lo < hi)
            mid = (lo + hi) / 2
            If (RintKey(mid) < ind) Then
                lo = mid + 1
            Else
                hi = mid
            End If
        End Do
        If (lo > Nint .OR. RintKey(lo) /= ind) Then
            strfmt = '(1X,"Fint: CANNOT FIND INTEGRAL ",3I4,I6)'
            Write(6,  strfmt) is, nfin, nini, ind
            Write(11, strfmt) is, nfin, nini, ind
            Stop
        End If
        Fint = RintVal(lo) * isg
    End Function Fint

    Subroutine DefJz2(idet)
        Implicit None
        Integer, Allocatable, Dimension(:), Intent(In) :: idet
        Integer :: m, i, k
        m = 0
        Do i = 1, Ne
            k = idet(i)
            m = m + Jz(k)
        End Do
        Jm0 = m / 2.d0
    End Subroutine DefJz2

    Subroutine Rint
        ! Read DTM.INT and build sorted (RintKey, RintVal) for bisect in Fint.
        Use utils, Only : DetermineRecordLength, CheckRecordLength
        Implicit None
        Integer :: i, is, ns1, nso1, err_stat, tmp_i
        Real(dp) :: rn1, x, z1, tmp_r
        Integer, Dimension(13) :: ki
        Integer, Allocatable, Dimension(:) :: l1
        Character(Len=256) :: strfmt, err_msg

        Open(13, file='DTM.INT', status='OLD', form='UNFORMATTED', iostat=err_stat, iomsg=err_msg)
        If (err_stat /= 0) Then
            Write(*,'(/2X,"Rint: DTM.INT absent"/)')
            Stop
        End If
        Read(13) ns1, nso1, z1, rn1
        Allocate(l1(ns1))
        x = iabs(ns1-Ns) + iabs(nso1-Nso) + dabs(z1-Z) + dabs(rn1-Rnuc)
        If (x < 1.d-6) Then
            Read(13) Nint, (l1(i), i=1,ns1)
            Allocate(RintVal(Nint), RintKey(Nint))
            Read(13) (ki(i), i=1,13)
            Read(13) (RintVal(i), RintKey(i), i=1,Nint)
            Close(13)
            Write(*,'(/1X,"### Radial integrals from DTM.INT ("," Nint =",I7," ) ###", &
                   /(4X,A4," calculated by ",A4))') Nint,(OpLabels(i),CorrLabels(ki(i)),i=1,13)
            Deallocate(l1)
            ! Insertion sort: sort RintKey ascending, permute RintVal correspondingly.
            Do i = 2, Nint
                tmp_i = RintKey(i)
                tmp_r = RintVal(i)
                x = Real(i-1, dp)
                Do While (x >= 1.d0 .AND. RintKey(Int(x)) > tmp_i)
                    RintKey(Int(x)+1) = RintKey(Int(x))
                    RintVal(Int(x)+1)  = RintVal(Int(x))
                    x = x - 1.0d0
                End Do
                RintKey(Int(x)+1) = tmp_i
                RintVal(Int(x)+1) = tmp_r
            End Do
            Return
        End If
        Write(*,*) 'Rint: wrong basis parameters in DTM.INT'
        Stop
    End Subroutine Rint

    Subroutine Init0
        ! Read X0/X2 eigenvectors from CONF0.XIJ.
        Implicit None
        Integer :: nx, n, i, err_stat
        Character(Len=256) :: err_msg

        Open(16, file='CONF0.XIJ', status='OLD', form='UNFORMATTED', iostat=err_stat, iomsg=err_msg)
        If (err_stat /= 0) Then
            Write(*,'(/2X,"Init0: CONF0.XIJ absent"/)')
            Stop
        End If
        Read(16) E0, Tj0, Nd2
        Allocate(X0(Nd0), X2(Nd0))
        Rewind(16)
        nx = N0
        If (N2 > nx) nx = N2
        Do n = 1, nx
            If (n == N0) Then
                Read(16) E0, Tj0, Nd2, (X0(i), i=1,Nd2)
                E0 = -E0
                If (N2 == N0) Then
                    E2 = E0
                    Tj2 = Tj0
                    X2 = X0
                End If
            Else If (n == N2) Then
                Read(16) E2, Tj2, Nd2, (X2(i), i=1,Nd2)
                E2 = -E2
            Else
                Read(16)
            End If
        End Do
        Close(16)
        Tj0 = Anint(Tj0*100.0d0)/100.0d0
        Tj2 = Anint(Tj2*100.0d0)/100.0d0
    End Subroutine Init0

    Subroutine SetParams
        Implicit None
        
        Nddir  = 15000
        N_it4    = 100
        N_it_ine = 100
        Npass    = 10
        Gj       = 0.d0
    End Subroutine SetParams

    Subroutine ReadIneIn
        ! Reads job parameters from ine.in.
        Implicit None
        Integer :: index_equals, index_hashtag, err_stat, i
        Character(len=20)  :: key, key_upper
        Character(len=256) :: val
        Character(len=256) :: line

        Open(unit=88, file='ine.in', status='OLD', iostat=err_stat)
        If (err_stat /= 0) Call WriteSkeleton('ine.in', [Character(len=90) :: &
            'RHS=2     # R.H.S. operator: 1=H_pnc, 2=E1(L), 5=E2',           &
            'LHS=2     # L.H.S. operator: 1=H_pnc, 2=E1(L), 3=H_am, 4=E1(V), 5=E2', &
            'N0=       # X0 record number in CONF0.XIJ',                       &
            'N2=       # X2 record number in CONF0.XIJ',                       &
            'Ranges=   # lambda1 lambda2 step [nm]; comma-sep for multiple (0 0 0 for static)', &
            '',                                                                 &
            '#Inuc=    # nuclear spin I (required when LHS=3)',                 &
            '#W00=     # E2 energy shift parameter (required when RHS=5)',      &
            '#Nddir=   # direct-solve block size (default 15000)',              &
            '#Passes=  # max restart passes (default 10)',                      &
            '#Iters=   # max iterations per pass (default 100)'])

        ! Defaults
        Kli    = 2
        Klf    = 2
        XInuc  = 0.0d0
        N0     = 0
        N2     = 0
        nrange = 0
        W0     = 0.0d0
        W00    = 0.0d0

        Do
            Read(88, '(A)', iostat=err_stat) line
            If (err_stat /= 0) Exit
            If (index(line, '=') == 0) Cycle
            If (index(AdjustL(line), '#') == 1) Cycle
            index_equals  = index(string=line, substring='=')
            key = Trim(AdjustL(line(1:index_equals-1)))
            val = line(index_equals+1:Len(line))
            index_hashtag = index(string=val, substring='#')
            If (index_hashtag /= 0) val = Trim(AdjustL(val(1:index_hashtag-1)))
            key_upper = ToUpperString(key)
            Select Case(key_upper)
                Case('RHS')
                    Read(val, *) Kli
                Case('LHS')
                    Read(val, *) Klf
                Case('INUC')
                    Read(val, *) XInuc
                Case('N0')
                    Read(val, *) N0
                Case('N2')
                    Read(val, *) N2
                Case('RANGES')
                    nrange = CountSubstr(val, ',') + 1
                    Allocate(xlamb1s(nrange), xlamb2s(nrange), xlambsteps(nrange))
                    Read(val, *, iostat=err_stat) (xlamb1s(i), xlamb2s(i), xlambsteps(i), i=1,nrange)
                    If (err_stat /= 0) Then
                        Write(*,*) 'ERROR: Ranges= could not be read from ine.in.'
                        Stop
                    End If
                Case('W00')
                    Read(val, *) W00
                Case('NDDIR')
                    Read(val, *) Nddir
                    Write(*,'(A,I6,A)') ' Nddir=', Nddir, ' overrides the default of 15000'
                Case('PASSES')
                    Read(val, *) Npass
                Case('ITERS')
                    Read(val, *) N_it4
                    N_it_ine = N_it4
                Case Default
                    Write(*,*) 'WARNING: unsupported key "', Trim(AdjustL(key)), '" in ine.in'
                End Select
        End Do
        Close(88)

        ! If no Ranges= key, set up a static (omega=0) single-point run
        If (.NOT. Allocated(xlamb1s)) Then
            nrange = 1
            Allocate(xlamb1s(nrange), xlamb2s(nrange), xlambsteps(nrange))
            xlamb1s = 0.0d0
            xlamb2s = 0.0d0
            xlambsteps = 1.0d0
        End If

        ! Echo
        Write(*,'(A,I1,A,I1)') ' RHS=', Kli, '  LHS=', Klf
        Write(*,'(A,I2,A,I2)') ' N0=', N0, '  N2=', N2
        Write(*,'(A,I0,A,I0,A)') ' Passes=', Npass, '  Iters=', N_it_ine, ' per pass'
        If (Klf == 3) Write(*,'(A,F6.2)') ' Inuc=', XInuc
        If (nrange > 0 .AND. xlamb1s(1) /= 0.0d0) Then
            Write(*,'(I2,A)') nrange, ' wavelength range(s):'
            Do i = 1, nrange
                Write(*,'(A,F11.3,A,F11.3,A,F11.4,A)') &
                    '  lambda=', xlamb1s(i), ' to ', xlamb2s(i), &
                    ' step ', xlambsteps(i), ' nm'
            End Do
        Else
            Write(*,'(A)') ' Static polarizability (omega=0)'
        End If
        If (N0 <= 0) Then
            Write(*,*) 'ERROR: N0 must be >= 1'
            Stop
        End If
    End Subroutine ReadIneIn

    Subroutine Input
        ! Read job parameters from ine.in and configuration from CONF.INP.
        Use conf_init, Only : ReadConfInp, ReadConfigurations
        Implicit None
        Character(Len=5), Dimension(5) :: str
        Data str / 'H_pnc','E1(L)','H_am ','E1(V)',' E2  ' /

        Call DetermineRecordLength(Mrec)
        Call ReadIneIn

        Call ReadConfInp
        Npv   = 3 + Nlv
        Nlft   = Npv
        Call ReadConfigurations
    End Subroutine Input

    Subroutine Init
        ! Read CONF.DAT; populate orbital quantum numbers, Nst, Nvc, Nc0.
        Use utils, Only : CheckRecordLength
        Use determinants, Only : Rdet
        Implicit None
        Integer  :: ic, n, j, imax, ni, kkj, llj, nnj, i, nj, if_, &
                    ii, nmin, jlj, i0, err_stat
        Real(dp) :: d, c1, c2, z1
        Logical  :: valid
        Real(dp), Dimension(IP6)   :: p, q, p1, q1
        Real(dp), Dimension(4*IP6) :: pq
        Character(Len=256) :: strfmt, err_msg
        Logical :: longbasis
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IPs)  :: Qq1
        Equivalence (IQN(1),PQ(21)), (Qq1(1),PQ(2*IPs+21))
        Equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))

        c1 = 0.01d0
        Mj = 2*abs(Jm)+0.01d0
        Open(12,file='CONF.DAT',status='OLD',access='DIRECT',recl=IP6*Mrec,iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file CONF.DAT is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If
        Call CheckRecordLength(12, 'CONF.DAT', IP6*Mrec)
        Read(12,rec=1) p
        Read(12,rec=2) q
        Read(12,rec=5) p1
        Read(12,rec=6) q1

        z1 = pq(1)
        If (abs(Z-z1) > 1.d-6) Then
            Write( 6,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Write(11,'("nuc. charge is changed: Z =",F12.6," ><",F12.6)') Z,z1
            Read(*,*)
        End If

        If (.not. Allocated(Nvc))  Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0))  Allocate(Nc0(Nc))
        If (.not. Allocated(Nq))   Allocate(Nq(Nsp))
        If (.not. Allocated(Nip))  Allocate(Nip(Nsp))
        Ns   = pq(2)+c1
        Ii   = pq(3)+c1
        H    = pq(6)
        Kt   = pq(9)+c1
        Rnuc = pq(13)
        longbasis = abs(PQ(20)-0.98765d0) < 1.d-6
        Write( 6,'(4X,"Kli =",I3,7X,"Klf =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3,5X,"Nc =",I6)') &
               Kli,Klf,Z,Jm,Nsp,Ns,Nso,Nc
        Write(11,'(4X,"Kli =",I3,7X,"Klf =",I3,7X,"Z   =",F6.2,4X,"Jm  =",F6.2, &
               /4X,"Nsp =",I7,5X,"Ns  =",I3,7X,"Nso =",I3,5X,"Nc =",I6)') &
               Kli,Klf,Z,Jm,Nsp,Ns,Nso,Nc
        If (longbasis) Then
            Write( *,*) ' Using variant for long basis '
            Write(11,*) ' Using variant for long basis '
            Do ni = 1, Ns
                Nn(ni) = IQN(4*ni-3)
                Ll(ni) = IQN(4*ni-2)
                Kk(ni) = IQN(4*ni-1)
                Jj(ni) = IQN(4*ni)
            End Do
        Else
            if_ = 20
            Do ni = 1, Ns
                if_ = if_+1
                Nn(ni) = pq(if_)+c1
                if_ = if_+1
                Ll(ni) = pq(if_)+c1
                if_ = if_+3
                c2 = dsign(c1,pq(if_))
                Kk(ni) = pq(if_)+c2
                if_ = if_+1
                c2 = dsign(c1,pq(if_))
                Jj(ni) = pq(if_)+c2
            End Do
        End If
        Nsu = 0
        Do nj = 1, Nsp
            i   = sign(1.d0, Qnl(nj))
            d   = abs(Qnl(nj))+1.d-14
            d   = 10.0*d
            nnj = d
            d   = 10.0d0*(d-nnj)
            llj = d
            jlj = 2*llj+i
            kkj = -i*((jlj+1)/2)
            d   = 100.0d0*(d-llj)
            Nq(nj) = d+0.1d0
            Do ni = 1, Ns
                If (nnj == Nn(ni) .and. Kk(ni) == kkj) Then
                    Exit
                Else If (ni == Ns) Then
                    Write( 6,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                    Write(11,'(/2X,"no orbital for shell ",I3,": n,l,k=",3I4)') nj,nnj,llj,kkj
                    Stop
                End If
            End Do
            Nip(nj) = ni
            If (Nsu < ni) Nsu = ni
        End Do
        Deallocate(Qnl)

        nec = 0
        If (Nso /= 0) nec = sum(Nq(1:Nso))
        Do ni = 1, Nsu
            imax = 2*Jj(ni)+1
            Do j = 1, imax, 2
                Nst = Nst+1
            End Do
        End Do
        Write( 6,'(4X,"Number of actually Used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
        Write(11,'(4X,"Number of actually Used orbitals: Nsu =",I3, &
             /4X,"Ne  =",I3,7X,"nec =",I3,7X,"Nst =",I7)') Nsu,Ne,nec,Nst
        n = 0
        ic = 0
        i0 = 0
        i = 0
        nmin = Nso+1
        Do ni = nmin, Nsp
            i = i+1
            n = n+Nq(ni)
            If (n >= Ne) Then
                ic = ic+1
                If (n > Ne) Then
                    Write( 6,'(/2X,"wrong number of electrons", &
                         /2X,"for configuration ICONF =",I6)') ic
                    Write(11,'(/2X,"wrong number of electrons", &
                         /2X,"for configuration ICONF =",I6)') ic
                    Stop
                End If
                Nvc(ic) = i
                Nc0(ic) = Nso+i0
                i0 = i0+i
                n = 0
                i = 0
            End If
        End Do
        Write( 6,'(1X,71("="))')
        Write(11,'(1X,71("="))')

        If (Gj /= 0.d0) Then
            Write(*,*) ' Gj =', Gj
            Write(*,*) ' This code works only for Gj=0!'
            Stop
        End If

        valid = (Kli==1.AND.Klf==2).OR.(Kli==1.AND.Klf==4).OR. &
                (Kli==2.AND.Klf<=3).OR.(Kli==5.AND.Klf==5)
        If (.NOT.valid) Then
            Write(*,*) ' Unknown combination Kli =',Kli,' Klf =',Klf
            Stop
        End If

        Call Rdet('CONF.DET')                ! sets Nd, allocates and fills Iarr
        Call Rdet('CONF0.DET', Iarr0, Nd0)   ! sets Nd0, allocates and fills Iarr0
    End Subroutine Init

    Subroutine Prj(str, tj0, X1, X1J)
        ! Decompose X1 into J0-1, J0, J0+1 components using J^2.
        Use mpi_f08
        Implicit None
        Character(Len=6), Intent(In) :: str
        Real(dp), Intent(In) :: tj0
        Real(dp), Allocatable, Dimension(:),   Intent(InOut) :: X1
        Real(dp), Allocatable, Dimension(:,:), Intent(InOut) :: X1J
        Integer :: i, jj, j2
        Real(dp) :: trsd, tj, tjm, tjp, tjj, sm, sj, sp, sn, sum, ds
        Real(dp) :: x, y, w, ap, aj, am
        Real(dp), Allocatable, Dimension(:) :: tmp
        Character(Len=256) :: strfmt

        If (.not. Allocated(X1J)) Allocate(X1J(Nd, Npv))
        If (.not. Allocated(tmp)) Allocate(tmp(Nd))
        trsd = 1.d-3
        jj   = 2*tj0 + 0.1d0
        j2   = jj*jj
        X1J(1:Nd,1:3) = 0.0d0

        If (mype == 0) Then
            strfmt = '(3X,"Prj: Decomposition of the vector ",A6, &
                /7X,"  Norma ",8X,"J_0-1",8X,"J_0",8X,"J_0+1",9X,"Sum")'
            Write( 6,strfmt) str
            Write(11,strfmt) str
        End If

        Call MxmpyJ(X1, X1J(:,1))           ! X1J(:,1) = J^2 * X1
        If (jj > 1) Then
            Call MxmpyJ(X1J(:,1), tmp)       ! tmp = J^4 * X1
            X1J(:,2) = tmp
        End If

        sm = 0.0d0
        sj = 0.0d0
        sp = 0.0d0
        sn = 0.0d0
        Do i = 1, Nd
            x = X1(i)
            y = X1J(i,1)
            w = X1J(i,2)
            If (jj > 1) Then
                ap =  (w - 2*j2*y + j2*(j2-4)*x) / (32*(jj+2)*(jj+1))
                aj = -(w - 2*(j2+2*jj+4)*y + jj*(j2-4)*(jj+4)*x) / (16*jj*(jj+2))
                am =  (w - 2*(jj+2)**2*y + jj*(jj+2)**2*(jj+4)*x) / (32*jj*(jj+1))
            Else
                ap =  (y - jj*(jj+2)*x)     / (4*(jj+2))
                aj = -(y - (jj+2)*(jj+4)*x) / (4*(jj+2))
                am = 0.0d0
            End If
            X1J(i,1) = am
            X1J(i,2) = aj
            X1J(i,3) = ap
            sm = sm + am*am
            sj = sj + aj*aj
            sp = sp + ap*ap
            sn = sn + x*x
        End Do
        sum = sm + sj + sp
        tj  = 0.5d0*jj
        tjm = tj-1.0d0
        tjp = tj+1.0d0
        tjj = (tjm*(tjm+1)*sm + tj*(tj+1)*sj + tjp*(tjp+1)*sp) / sum
        tjj = 0.5d0*(dsqrt(1.0d0 + 4*tjj) - 1.0d0)
        If (mype == 0) Then
            strfmt = '(2X,5E13.5,/4X,"J_av=",F11.8)'
            Write( 6,strfmt) sn,sm,sj,sp,sum,tjj
            Write(11,strfmt) sn,sm,sj,sp,sum,tjj
        End If
        ds = dabs(sum/sn - 1.0d0)
        If (ds > trsd .AND. mype == 0) Write(*,*) ' Prj warning: normalization changed by', ds
        Deallocate(tmp)
    End Subroutine Prj

    Subroutine PrjE2(str, tj0, X1, X1J)
        ! E2 decomposition: X1 = X_J-2 + X_J-1 + X_J + X_J+1 + X_J+2.
        Use mpi_f08
        Implicit None
        Character(Len=6), Intent(In) :: str
        Real(dp), Intent(In) :: tj0
        Real(dp), Allocatable, Dimension(:),   Intent(InOut) :: X1
        Real(dp), Allocatable, Dimension(:,:), Intent(InOut) :: X1J
        Integer :: i, j0
        Real(dp) :: sum, sum1, sum2, sum3, sum4, sum5, sum_tot, ds, f
        Real(dp) :: XJ3xXJ4, XJ3xXJ5, XJ4xXJ5
        Real(dp), Allocatable, Dimension(:) :: TJ1, TJ2, TJ3, TJ4, TJ5, tmp
        Character(Len=256) :: strfmt

        If (.not. Allocated(X1J)) Allocate(X1J(Nd, Npv))
        Allocate(TJ1(Nd), TJ2(Nd), TJ3(Nd), TJ4(Nd), TJ5(Nd), tmp(Nd))
        j0 = 2*tj0 + 0.1d0
        TJ1 = 0.0d0
        TJ2 = 0.0d0
        TJ3 = 0.0d0
        TJ4 = 0.0d0
        TJ5 = 0.0d0
        X1J(1:Nd,1:5) = 0.0d0

        If (mype == 0) Then
            strfmt = '(3X,"PrjE2: Decomposition of the vector ",A6, &
                /5X," Norma ",7X,"J0-2",9X,"J0-1",10X,"J0",10X,"J0+1",9X,"J0+2",9X,"Sum")'
            Write(*,*)
            Write( 6,strfmt) str
            Write(11,strfmt) str
        End If

        Call MxmpyJ(X1, X1J(:,1))
        Do i=1,Nd
            TJ1(i)=X1J(i,1)-(j0+4)*(j0+6)*X1(i)
            X1J(i,1)=0.0d0
        End Do

        Call MxmpyJ(TJ1, X1J(:,1))
        Do i=1,Nd
            TJ2(i)=X1J(i,1)-(j0+2)*(j0+4)*TJ1(i)
            X1J(i,1)=0.0d0
        End Do

        Call MxmpyJ(TJ2, X1J(:,1))
        Do i=1,Nd
            TJ3(i)=X1J(i,1)-j0*(j0+2)*TJ2(i)
            X1J(i,1)=0.0d0
        End Do

        If (tj0 >= 1.5d0) Then
            Call MxmpyJ(TJ3, X1J(:,2))
            f = -1.0d0/(1536*(j0-2)*j0*(j0+1)*(j0+2))
            Do i=1,Nd
                TJ4(i)=X1J(i,2)-(j0-4)*(j0-2)*TJ3(i)
                X1J(i,2)=f*TJ4(i)
            End Do
        Else
            X1J(1:Nd,2) = 0.0d0
        End If

        If (tj0 >= 2.0d0) Then
            Call MxmpyJ(TJ3, X1J(:,1))
            f = 1.0d0/(6144*(j0-2)*(j0-1)*j0*(j0+1))
            Do i=1,Nd
                TJ5(i)=X1J(i,1)-(j0-2)*j0*TJ3(i)
                X1J(i,1)=f*TJ5(i)
            End Do
        Else
            X1J(1:Nd,1) = 0.0d0
        End If

        f = 1.0d0/(2*(j0+2)*(j0+3))
        Do i=1,Nd
            X1J(i,3) = f*(TJ2(i)/16 - 6*(j0+1)*(j0+2)*X1J(i,2) - 12*j0*(j0+1)*X1J(i,1))
        End Do
        f = -1.0d0/(j0+4)
        Do i=1,Nd
            X1J(i,4) = f*(TJ1(i)/4 + 2*(j0+3)*X1J(i,3) + 3*(j0+2)*X1J(i,2) + 4*(j0+1)*X1J(i,1))
        End Do
        Do i=1,Nd
            X1J(i,5) = X1(i)-X1J(i,4)-X1J(i,3)-X1J(i,2)-X1J(i,1)
        End Do

        sum=0.0d0
        sum1=0.0d0
        sum2=0.0d0
        sum3=0.0d0
        sum4=0.0d0
        sum5=0.0d0
        XJ3xXJ4=0.0d0
        XJ3xXJ5=0.0d0
        XJ4xXJ5=0.0d0
        Do i=1,Nd
            sum=sum+X1(i)*X1(i)
            sum1=sum1+X1J(i,1)**2
            sum2=sum2+X1J(i,2)**2
            sum3=sum3+X1J(i,3)**2
            sum4=sum4+X1J(i,4)**2
            sum5=sum5+X1J(i,5)**2
            XJ3xXJ4=XJ3xXJ4+X1J(i,3)*X1J(i,4)
            XJ3xXJ5=XJ3xXJ5+X1J(i,3)*X1J(i,5)
            XJ4xXJ5=XJ4xXJ5+X1J(i,4)*X1J(i,5)
        End Do
        sum_tot = sum1+sum2+sum3+sum4+sum5
        If (mype == 0) Then
            strfmt = '(7E13.5)'
            Write( 6,strfmt) sum,sum1,sum2,sum3,sum4,sum5,sum_tot
            Write(11,strfmt) sum,sum1,sum2,sum3,sum4,sum5,sum_tot
            ds = dabs(sum/sum_tot - 1.0d0)
            If (ds > 1.d-5) Write(*,*) ' PrjE2 warning: normalization changed by', ds
        End If
        Deallocate(TJ1, TJ2, TJ3, TJ4, TJ5, tmp)
    End Subroutine PrjE2

    Subroutine RdcX1J
        ! Scale X1J columns by Wigner 3j factors to give reduced matrix elements.
        Use wigner
        Implicit None
        Integer :: i, jf, j, is
        Real(dp) :: Aj, xm, W

        jf = 3
        If (Kli == 5) jf = 5
        Do j = 1, jf
            If (Kli == 5) Then
                Aj = Tj0+(j-3)
            Else
                Aj = Tj0+(j-2)
            End If
            Aj = Anint(Aj*100.0)/100.0
            If (Kli /= 1) Then
                xm = dabs(Jm) - 1.d-5
                If (Aj > xm) Then
                    is = Aj - Jm + 1.d-5
                    is = 1 - 2*mod(is,2)
                    If (Kli == 5) Then
                        W = is/(FJ3(Aj,2.d0,Tj0,-Jm,Q,Jm0)+1.d-77)
                    Else
                        W = is/(FJ3(Aj,1.d0,Tj0,-Jm,Q,Jm0)+1.d-77)
                    End If
                    If (dabs(W) > 1.d+5) W = 0.0d0
                Else
                    W = 0.0d0
                End If
                Do i = 1, Nd
                    X1J(i,j) = W*X1J(i,j)
                End Do
            End If
        End Do
    End Subroutine RdcX1J

    Subroutine PrintDiag
        ! Print J-values of X1/Y1/Y2 (J^2 SpMV on all ranks) and per-state
        ! reduced matrix elements and alpha contributions from X1J/dEv/tjv (rank 0).
        Use mpi_f08
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i, ic, j, k, n, jtj, ndk, is
        Real(dp) :: xn, yn, zn, dx, dy, dz, xj, yj, zj, t, Ex, E, c1, c2, w3j, w2, &
                    al0, al2, tj
        Real(dp), Allocatable, Dimension(:) :: tmp
        Character(Len=5), Dimension(5) :: st
        Character(Len=256) :: strfmt
        Data st /'H_pnc','E1(L)','H_am ','E1(V)',' E2  '/

        Allocate(tmp(Nd))
        If (mype == 0) Then
            strfmt = '(3X,63("="))'
            Write( 6,strfmt)
            Write(11,strfmt)
            If (W0 == 0.0d0) Then
                strfmt = '(4X,"solution of the eq-n (E-H0)|X1> = ",A5,"|X0>," &
                    /4X,"where    E =",F13.6,",", &
                    /10X,"|X0> = |",I2,","F13.6,",",2F4.1,">.", &
                    /3X,63("-"))'
                Write( 6,strfmt) st(Kli),Elft,N0,E0,Tj0,Jm0
                Write(11,strfmt) st(Kli),Elft,N0,E0,Tj0,Jm0
            Else
                strfmt = '(4X,"solution of the eq-n (E+w-H0)|X1> = ",A5,"|X0>," &
                    /4X,"where E+w =",F13.6," +"F10.6," =",F13.6, &
                    /10X,"|X0> = |",I2,","F13.6,",",2F4.1,">.", &
                    /3X,63("-"))'
                Write( 6,strfmt) st(Kli),E0,W0,Elft,N0,E0,Tj0,Jm0
                Write(11,strfmt) st(Kli),E0,W0,Elft,N0,E0,Tj0,Jm0
            End If
            strfmt = '(6X,"Nc",5X,"|",A5,"|X0>",8X,"|X1>",8X,"<X2|",A5,"|",/3X,63("-"))'
            Write( 6,strfmt) st(Kli),st(Klf)
            Write(11,strfmt) st(Kli),st(Klf)
            i=0
            xn=0.0d0
            yn=0.0d0
            zn=0.0d0
            Do ic=1,Nc
                dx=0.0d0
                dy=0.0d0
                dz=0.0d0
                ndk=Ndc(ic)
                Do k=1,ndk
                    i=i+1
                    dx=dx+X1(i)**2
                    dy=dy+Y1(i)**2
                    dz=dz+Y2(i)**2
                End Do
                xn=xn+dx
                yn=yn+dy
                zn=zn+dz
            End Do
            Write(11,'(3X,63("-"))')
        End If

        ! J-expectation values via MxmpyJ (all ranks)
        xj=0.0d0
        yj=0.0d0
        zj=0.0d0
        Call MxmpyJ(X1, tmp)
        Do i=1,Nd
            xj=xj+X1(i)*tmp(i)
        End Do
        Call MxmpyJ(Y1, tmp)
        Do i=1,Nd
            yj=yj+Y1(i)*tmp(i)
        End Do
        Call MxmpyJ(Y2, tmp)
        Do i=1,Nd
            zj=zj+Y2(i)*tmp(i)
        End Do

        If (mype == 0) Then
            xn=0.0d0
            yn=0.0d0
            zn=0.0d0
            Do i=1,Nd
                xn=xn+X1(i)**2
                yn=yn+Y1(i)**2
                zn=zn+Y2(i)**2
            End Do
            xj=0.5d0*(dsqrt(1.0d0+xj/xn)-1.0d0)
            yj=0.5d0*(dsqrt(1.0d0+yj/yn)-1.0d0)
            zj=0.5d0*(dsqrt(1.0d0+zj/zn)-1.0d0)
            strfmt = '(4X," J =",3F15.8)'
            Write( 6,strfmt) yj,xj,zj
            Write(11,strfmt) yj,xj,zj
        End If
        Deallocate(tmp)

        If (mype /= 0) Return
        ! --- Rank 0 only: table of matrix elements vs probe vectors in X1J ---
        If (Kli == 5) Then
            strfmt = '(3X,78("="),/4X,"n",8X,"J_n",11X,"Ev_n",5X, &
                    "<X2||",A5,"||n>",2X,"<n||",A5,"||X0>",3X, &
                    "contrib. |n>",/69X," to alpha0"/3X,78("-"))'
            Write( 6,strfmt) st(Klf),st(Kli)
            Write(11,strfmt) st(Klf),st(Kli)
        Else If (Kli == 2 .AND. Klf == 2) Then
            strfmt = '(3X,94("="),/4X,"n",8X,"J_n",11X,"Ev_n",5X, &
                     "<X2||",A5,"||n>",2X,"<n||",A5,"||X0>",3X, &
                     "contrib. |n>",4X,"contrib. |n>", &
                     /69X," to alpha_0",5X," to alpha_2",/3X,94("-"))'
            Write( 6,strfmt) st(Klf),st(Kli)
            Write(11,strfmt) st(Klf),st(Kli)
        Else
            strfmt = '(3X,63("="),/4X,"n",9X,"J_n",10X,"Ev_n",5X, &
                     "<X2#",A5,"#n>",4X,"<n#",A5,"#X0>",/3X,63("-"))'
            Write( 6,strfmt) st(Klf),st(Kli)
            Write(11,strfmt) st(Klf),st(Kli)
        End If
        Ex = E0
        If (Kli==2 .AND. Klf/=2) Ex = E2
        j=0
        jtj=0
        Do n=1,Nlft-3
            E  = -dEv(n+3)
            jtj = 2*tjv(n+3)+1.d-3
            tj  = 0.5d0*jtj
            j=j+1
            c1=0.0d0
            c2=0.0d0
            w3j=0.0d0
            Do i=1,Nd
                c1=c1+X1(i)*X1J(i,n+3)*(Ex-E+W0)
                c2=c2+Y2(i)*X1J(i,n+3)
            End Do
            If (Klf==4) c2=c2/(E-E2)
            If (Kli/=1) Then
                k=tj-Jm+1.d-1
                is=1
                If (k/=2*(k/2)) is=-1
                If (Kli==5) Then
                    w3j=is*FJ3(tj,2.d0,Tj0,-Jm,Q,Jm0)+1.d-77
                Else
                    w3j=is*FJ3(tj,1.d0,Tj0,-Jm,Q,Jm0)+1.d-77
                End If
                c1=c1/w3j
                If (dabs(c1)>1.d+10) c1=0.0d0
            End If
            If (Klf/=1) Then
                k=Tj2-Jm0+1.d-1
                is=1
                If (k/=2*(k/2)) is=-1
                If (Klf==5) Then
                    w2=is*FJ3(Tj2,2.d0,tj,-Jm0,-Q,Jm)+1.d-77
                Else
                    w2=is*FJ3(Tj2,1.d0,tj,-Jm0,-Q,Jm)+1.d-77
                End If
                c2=c2/w2
                If (dabs(c2)>1.d+10) c2=0.0d0
            End If
            If (Kli==5) Then
                If (dabs(E-Ex)<1.d-5) Then
                    al0 = 0.0d0
                    If (W0/=0.0d0) al0 = 2.0d0/(5*(2*Tj0+1))*c1**2/(E-Ex-W0)
                Else
                    al0 = 2.0d0/(5*(2*Tj0+1))*c2**2/(E-Ex-W0)
                End If
                strfmt = '(3X,I2,F15.9,F13.6,3(1pE16.6))'
                Write( 6,strfmt) j,tj,E,c2,c1,al0
                Write(11,strfmt) j,tj,E,c2,c1,al0
            Else If (Kli==2 .AND. Klf==2) Then
                al0 = 2.0d0/(3*(2*Tj0+1))*c2**2/(E-Elft)
                is  = tj+Tj0+0.1d0
                is = 1-2*mod(is,2)
                al2 = 4*dsqrt(5*Tj0*(2*Tj0-1)/(6*(2*Tj0+3)*(2*Tj0+1)*(Tj0+1))) &
                      *FJ6(Tj0,1.d0,tj,1.d0,Tj0,2.d0)*is*c2**2/(E-Elft)
                strfmt = '(3X,I2,F15.9,F13.6,4(1pE16.6))'
                Write( 6,strfmt) j,tj,E,c2,c1,al0,al2
                Write(11,strfmt) j,tj,E,c2,c1,al0,al2
            Else If (Kli/=5) Then
                strfmt = '(3X,I2,F15.9,F13.6,2E16.6)'
                Write( 6,strfmt) j,tj,E,c2,c1
                Write(11,strfmt) j,tj,E,c2,c1
            End If
        End Do
        If (Kli==5 .OR. (Kli==2.AND.Klf==2)) Then
            strfmt = '(3X,94("="))'
            If (Kli==5) strfmt = '(3X,78("="))'
            Write( 6,strfmt)
            Write(11,strfmt)
        Else
            strfmt = '(3X,63("="))'
            Write( 6,strfmt)
            Write(11,strfmt)
        End If
    End Subroutine PrintDiag

    Subroutine RdcE1(n)
        ! E1 polarizability: scalar alpha_0, vector alpha_1, tensor alpha_2 (Kli=Klf=2).
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer, Intent(In) :: n
        Integer :: i, isk, id, k, kmin, kmax, is
        Real(dp) :: f0, tj, w3j, f, f1, f2, al, al0, al1, al2, alp
        Real(dp), Dimension(3) :: sk0, sk1, sk2
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        If (.not. Allocated(Y2J)) Allocate(Y2J(Nd,5))
        If (dabs(Tj2-Tj0) > 1.d-5) Then
            Write(*,*) ' can not calculate polarizability for J2 /= J0'
            Return
        End If
        f0 = -2.0d0
        ss(n) = 0.0d0
        strfmt = '(3X,"alpha(J=",F4.1," M=",F4.1,")=",E12.5," (no Prj)")'
        If (dabs(Q) < 1.d-5) Then
            Do i=1,Nd
            ss(n)=ss(n)+Y2(i)*X1(i)*f0
        End Do
            Write(*,*)
            Write( 6,strfmt) Tj0,Jm0,ss(n)
            Write(11,strfmt) Tj0,Jm0,ss(n)
        End If
        isk=0
        s0(n)=0.0d0
        s1(n)=0.0d0
        s2(n)=0.0d0
        Call DefSum(id,kmin,kmax)
        strfmt2 = '(1X,"J =",F4.1," to alpha_0:",E12.5,", alpha_1:",E12.5,", alpha_2:",E12.5)'
        Do k=kmin,kmax
            tj=Tj0+(k-2)
            tj=Anint(tj*100.0)/100.0
            is=dabs(tj-Jm0)+1.d-5
            is=1-2*mod(is,2)
            w3j=is*FJ3(Tj0,1.d0,tj,-Jm0,-Q,Jm)
            If (w3j /= 0.0d0) Then
                f=f0/w3j/(3.0d0*(2*Tj0+1))
            Else
                If (Tj0+tj > 0.99d0) isk=isk+1
                f=0.0d0
            End If
            sk0(k)=0.0d0
            sk1(k)=0.0d0
            sk2(k)=0.0d0
            Do i=1,Nd
                sk0(k)=sk0(k)+f*Y2J(i,k)*X1J(i,k)
            End Do
            is=tj+Tj0+1.1d0
            is=1-2*mod(is,2)
            If (W0 /= 0.0d0) Then
                f1=is*1.5d0*dsqrt(6.0d0*Tj0*(2*Tj0+1)/(Tj0+1)) &
                   *FJ6(Tj0,1.d0,tj,1.d0,Tj0,1.d0)
            Else
                f1=0.0d0
            End If
            sk1(k)=f1*sk0(k)
            f2=-2.0d0*is*dsqrt(15.0d0*Tj0*(2*Tj0-1)*(2*Tj0+1) &
               /(2.0d0*(2*Tj0+3)*(Tj0+1)))*FJ6(Tj0,1.d0,tj,1.d0,Tj0,2.d0)
            sk2(k)=f2*sk0(k)
            Write( 6,strfmt2) tj,sk0(k),sk1(k),sk2(k)
            Write(11,strfmt2) tj,sk0(k),sk1(k),sk2(k)
            s0(n)=s0(n)+sk0(k)
            s1(n)=s1(n)+sk1(k)
            s2(n)=s2(n)+sk2(k)
        End Do
        If (Tj0 > 0.51d0) Then
            s(n) = s0(n)+(3*Jm0*Jm0-Tj0*(Tj0+1)) &
                        /(Tj0*(2*Tj0-1.d0))*s2(n)
        Else
            s(n) = s0(n)
        End If
        strfmt3 = '(3X,"Polarizability Alpha( J=",F4.1," M=",F4.1,") =",E12.5, &
                  " a.u. (with Prj)",/3X,"=",E12.5," +(",E12.5, &
                  ") (3M^2-J(J+1))/J*(2J-1) a.u.",/3X,"Alpha_1 =",E12.5," a.u.")'
        If (isk == 0) Then
            Write( 6,strfmt3) Tj0,Jm0,s(n),s0(n),s2(n),s1(n)
            Write(11,strfmt3) Tj0,Jm0,s(n),s0(n),s2(n),s1(n)
        End If
        If (n == 2) Then
            alp=(ss(1)+ss(2))/2.0d0
            al =(s(1)+s(2))/2.0d0
            al0=(s0(1)+s0(2))/2.0d0
            al2=(s2(1)+s2(2))/2.0d0
            al1=(s1(1)-s1(2))
            strfmt = '(3X,9("-")/,3X,"In total:",/3X,9("-"))'
            Write( 6,strfmt)
            Write(11,strfmt)
            If (dabs(Q)<1.d-5) Then
                strfmt = '(3X,"alpha(J=",F4.1," M=",F4.1,")=",E12.5," (no Prj)")'
                Write( 6,strfmt) Tj0,Jm0,alp
                Write(11,strfmt) Tj0,Jm0,alp
            End If
            Write( 6,strfmt3) Tj0,Jm0,al,al0,al2,al1
            Write(11,strfmt3) Tj0,Jm0,al,al0,al2,al1
            If (converged) Then
                If (xlamb == 0.0d0) Then
                    strfmt2 = '(/1X,"RESULT: omega =",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  CONVERGED")'
                Else
                    strfmt2 = '(/1X,"RESULT: lambda=",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  CONVERGED")'
                End If
            Else
                If (xlamb == 0.0d0) Then
                    strfmt2 = '(/1X,"RESULT: omega =",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  DIVERGED")'
                Else
                    strfmt2 = '(/1X,"RESULT: lambda=",F14.6," alpha_0=",F17.7," alpha_2=",F17.7,"  DIVERGED")'
                End If
            End If
            Write( 6,strfmt2) abs(xlamb),al0,al2
            Write(11,strfmt2) abs(xlamb),al0,al2
        End If
    End Subroutine RdcE1

    Subroutine RdcE2(n)
        ! E2 scalar polarizability alpha_E2 (Kli=Klf=5).
        Use wigner, Only : FJ3
        Implicit None
        Integer, Intent(In) :: n
        Integer :: i, isk, k, is
        Real(dp) :: f0, tj, w3j, f, al0, alp
        Real(dp), Dimension(5) :: sk0
        Real(dp), Dimension(2) :: s0, ss
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        f0 = -2.0d0
        ss(n) = 0.0d0
        If (dabs(Q) < 1.d-5) Then
            Do i=1,Nd
            ss(n)=ss(n)+Y2(i)*X1(i)*f0
        End Do
        End If
        isk=0
        s0(n)=0.0d0
        Do k=1,5
            tj=Tj0+(k-3)
            tj=Anint(tj*100.0)/100.0
            If (tj > dabs(Jm)-1.d-5) Then
                is=dabs(tj-Jm0)+1.d-5
                is=1-2*mod(is,2)
                w3j=is*FJ3(Tj0,2.d0,tj,-Jm0,-Q,Jm)
                If (w3j /= 0.0d0) Then
                    f=f0/w3j/(5.0d0*(2*Tj0+1))
                Else
                    f=0.0d0
                End If
            Else
                f=0.0d0
            End If
            sk0(k)=0.0d0
            Do i=1,Nd
                sk0(k)=sk0(k)+f*Y2J(i,k)*X1J(i,k)
            End Do
            If (sk0(k)/=0.0d0) Then
                strfmt = '(1X,"con-n of J =",F5.2," to scalar alpha_E2:",E14.7)'
                Write( *,strfmt) tj,sk0(k)
                Write(11,strfmt) tj,sk0(k)
            End If
            s0(n)=s0(n)+sk0(k)
        End Do
        Write(*,*)
        strfmt2 = '(3X,"alpha_E2(J=",F4.1," M=",F4.1,")=",1pE15.7," (no Prj)")'
        strfmt3 = '(3X,"Scalar part of alpha_E2( J=",F4.1," M=",F4.1,") =",1pE15.7," a.u. (with Prj)")'
        If (dabs(Q)<1.d-5) Then
            Write( 6,strfmt2) Tj0,Jm0,ss(n)
            Write(11,strfmt2) Tj0,Jm0,ss(n)
        End If
        Write( 6,strfmt3) Tj0,Jm0,s0(n)
        Write(11,strfmt3) Tj0,Jm0,s0(n)
        If (n==2) Then
            alp = (ss(1)+ss(2))/2.0d0
            al0 = (s0(1)+s0(2))/2.0d0
            strfmt = '(3X,9("-")/,3X,"In total:",/3X,9("-"))'
            Write( 6,strfmt)
            Write(11,strfmt)
            If (dabs(Q)<1.d-5) Then
                Write( 6,strfmt2) Tj0,Jm0,alp
                Write(11,strfmt2) Tj0,Jm0,alp
            End If
            Write( 6,strfmt3) Tj0,Jm0,al0
            Write(11,strfmt3) Tj0,Jm0,al0
        End If
    End Subroutine RdcE2

    Subroutine RdcPNC
        ! Reduced PNC amplitude <J2 || E_pnc || J0> (Kli=1 or Klf=1).
        Use wigner, Only : FJ3
        Implicit None
        Integer :: i, k, is, id, kmin, kmax, kf
        Real(dp) :: gap, nuc, a, f0, f, s, ss, tj, w3j
        Real(dp), Dimension(3) :: s1
        Character(Len=9), Dimension(4) :: str
        Character(Len=256) :: strfmt
        Data str /'| H_pnc |','||E1(L)||','||H_am ||','||E1(V)||'/

        If (min(Kli,Klf)/=1) Return
        gap = DPcw/16.0d0
        nuc = Am - Z + 1.d-8
        Call DefSum(id,kmin,kmax)
        k=Tj2-Jm0+1.d-1
        is=1
        If (k/=2*(k/2)) is=-1
        a=is/(FJ3(Tj2,1.d0,Tj0,-Jm0,0.0d0,Jm0)+1.d-77)
        f0=a*gap*(-nuc)
        If (Kli==1.AND.Klf==4) f0=f0/(E0-E2)
        ss=0.0d0
        Do i=1,Nd
            ss=ss+Y2(i)*X1(i)*f0
        End Do
        s=0.0d0
        Do k=kmin,kmax
            kf=k-id
            tj=Tj0+(k-2)
            tj=Anint(tj*100.0)/100.0
            If (Kli==1) f=f0
            If (Klf==1) Then
                is=tj-Jm+1.d-5
                is=1-2*mod(is,2)
                w3j=FJ3(tj,1.d0,Tj0,-Jm,Q,Jm0)
                f=is*f0*w3j
            End If
            s1(k)=0.0d0
            Do i=1,Nd
                s1(k)=s1(k)+f*Y2J(i,kf)*X1J(i,k)
            End Do
            strfmt = '(3X,"<",F5.2,A9,"R(J=",F5.2,")",A9,F5.2,"> = ",E12.5," (a.u.)")'
            Write( 6,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
            Write(11,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
            s=s+s1(k)
        End Do
        Write( 6,'(3X,63("-"))')
        Write(11,'(3X,63("-"))')
        strfmt = '(4X,"<",F4.1," ||E_pnc ||",F4.1,"> = ", &
           /4X,"= i*Q_w/(-N_nuc) * (",E12.5,") a.u. (using Prj)", &
           /4X,"= i*Q_w/(-N_nuc) * (",E12.5,") a.u. (without Prj)", &
           /4X,"(N_nuc =",I3,")")'
        Write( 6,strfmt) Tj2,Tj0,s,ss,Int(nuc)
        Write(11,strfmt) Tj2,Tj0,s,ss,Int(nuc)
        Write( 6,'(3X,63("="))')
        Write(11,'(3X,63("="))')
    End Subroutine RdcPNC

    Subroutine RdcAM
        ! Anapole moment amplitude from E1 and H_am (Kli=2, Klf=3).
        Use wigner, Only : FJ3, FJ6
        Implicit None
        Integer :: i, id, is, ix, isk, k, kf, ks, kx, kj, kmin, kmax
        Real(dp) :: gap, f0, s, tj, w2, f, Fn0, Fx0, Fn2, Fx2, c0, c, cm, c1, &
                    F20, F21, F2, w6, a0, a
        Real(dp), Dimension(3) :: s1
        Character(Len=9), Dimension(4) :: str
        Character(Len=256) :: strfmt
        Data str /'| H_pnc |','||E1(L)||','||H_am ||','||E1(V)||'/

        If (Kli/=2.OR.Klf/=3) Return
        gap = DPcw/16.0d0
        Call DefSum(id,kmin,kmax)
        k=Tj2-Jm0+1.1d0
        is=1
        If (k/=2*(k/2)) is=-1
        f0=is*gap*2*dsqrt(6.0d0)
        s=0.0d0
        isk=0
        Do k=kmin,kmax
            kf=k-id
            tj=Tj0+(k-2)
            tj=Anint(tj*100.0)/100.0
            w2=FJ3(Tj2,1.d0,tj,-Jm0,-Q,Jm)
            If (w2/=0.0d0) Then
                f=f0/w2
            Else
                isk=isk+1
                f=0.0d0
            End If
            s1(k)=0.0d0
            Do i=1,Nd
                s1(k)=s1(k)+f*Y2J(i,kf)*X1J(i,k)
            End Do
            strfmt = '(3X,"<",F5.2,A9,"R(J=",F5.2,")",A9,F5.2,"> = ",E12.5," (a.u.)")'
            Write( 6,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
            Write(11,strfmt) Tj2,str(Klf),tj,str(Kli),Tj0,s1(k)
            s=s+s1(k)
        End Do
        Write( 6,'(3X,63("-"))')
        Write(11,'(3X,63("-"))')
        Fn0=dabs(Tj0-XInuc)
        Fx0=Tj0+XInuc
        ix=Fx0-Fn0+1.1d0
        Fn2=dabs(Tj2-XInuc)
        Fx2=Tj2+XInuc
        c0=dsqrt((XInuc+1)*(2*XInuc+1)*XInuc)
        Do i=1,ix
            F0=Fn0+i-1
            F20=F0-1
            If(F20<Fn2) F20=Fn2
            F21=F0+1
            If(F21>Fx2) F21=Fx2
            kx=F21-F20+1.1d0
            If(kx>3) kx=3
            Do k=1,kx
                F2=F20+(k-1)
                w6=FJ6(Tj2,F2,XInuc,F0,Tj0,1.d0)
                If (w6/=0.0d0) Then
                    c1=c0/w6
                Else
                    Write(*,*) ' w6=0'
                    c1=1.d99
                End If
                ks=F0+Tj2+XInuc+1.1d0
                is=1+2*(2*(ks/2)-ks)
                cm=is*dsqrt((2*F2+1)*(2*F0+1))*w6
                a=0.0d0
                Do kj=kmin,kmax
                    tj=Tj0+(kj-2)
                    ks=XInuc+Tj2+F2+0.1d0
                    is=1+2*(2*(ks/2)-ks)
                    c=is*c1*FJ6(tj,F2,XInuc,F0,Tj0,1.d0)*FJ6(F2,XInuc,Tj2,1.d0,tj,XInuc)
                    a0=c*s1(kj)
                    a=a+a0
                End Do
                strfmt = '(4X,"<",F4.1,"||E_am||",F4.1,"> = ","i*kappa/I*(",E12.5,")*(",E12.5,") [Porsev]")'
                Write( 6,strfmt) F2,F0,a,cm
                Write(11,strfmt) F2,F0,a,cm
            End Do
        End Do
        Write( 6,'(3X,63("="))')
        Write(11,'(3X,63("="))')
    End Subroutine RdcAM

    Subroutine C_3
        ! C_3 coefficient.
        Use wigner, Only : FJ3
        Implicit None
        Integer :: n, ix, jj, i
        Real(dp) :: c, c3, aj, w
        Real(dp), Dimension(3) :: sk2
        Character(Len=256) :: strfmt

        If (Kli/=2.OR.Klf/=2) Return
        strfmt = '(/3X,"C_3 coefficient for the state |",I2,",",F12.6,",",F4.1,",",F4.1,">")'
        Write( 6,strfmt) N2,E2,Tj2,Jm0
        Write(11,strfmt) N2,E2,Tj2,Jm0
        jj=2*Tj2+0.1d0
        c=1.0d0/(12*(jj+1))
        ix=3
        If (jj<=1) ix=2
        c3=0.0d0
        strfmt = '(3X,"contribution of J =",F4.1,": ",E12.5," a.u.")'
        Do i=1,ix
            aj=Tj2+i+1-ix
            w=(FJ3(aj,1.d0,Tj2,-Jm,Q,Jm0))**2
            If (w>0.0d0) w=c/w
            sk2(i)=0.0d0
            Do n=1,Nd
                sk2(i)=sk2(i)+w*Y2J(n,i+3-ix)**2
            End Do
            c3=c3+sk2(i)
            Write( 6,strfmt) aj,sk2(i)
            Write(11,strfmt) aj,sk2(i)
        End Do
        strfmt = '(3X,"C_3 =",E12.5," a.u.")'
        Write( 6,strfmt) c3
        Write(11,strfmt) c3
    End Subroutine C_3

    Subroutine DefSum(id, kmin, kmax)
        ! Define range of intermediate J values.
        Implicit None
        Integer, Intent(Out) :: id, kmin, kmax
        Integer :: imax, ki0, ki1, kf0, kf1
        Real(dp) :: dj
        dj=Tj2-Tj0
        id=dabs(dj)+0.1d0
        If (dj<0.0d0) id=-id
        imax=2
        If (Kli==1) Then
            ki0=2
            ki1=2
            imax=imax-1
        Else
            ki0=1
            ki1=3
        End If
        If (Tj0<0.99d0.AND.ki0<2) ki0=2
        If (Klf==1) Then
            kf0=2+id
            kf1=2+id
            imax=imax-1
        Else
            kf0=1+id
            kf1=3+id
        End If
        If (Tj2<0.99d0.AND.kf0<(2+id)) kf0=2+id
        kmin=max(ki0,kf0)
        kmax=min(ki1,kf1)
        If (kmax<kmin.OR.iabs(id)>imax) Then
            Write(*,'(" For Kli =",I1," and Klf =",I1," no transition",F5.2,"-->",F5.2, &
                /" kmin =",I2," kmax =",I2," id =",I2)') Kli,Klf,Tj0,Tj2,kmin,kmax,id
        End If
    End Subroutine DefSum

End Program pine
