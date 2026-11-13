Module davidson
    ! This module contains main subroutines required for the Davidson procedure.

    Use conf_variables
    Use hblocks, Only : h2d_active, BuildH2D, MatvecH2D, MvMark, MvWait, MvCount

    Implicit None

    Private

    ! Flag to ensure Diag is filled exactly once per Davidson run.
    Logical, Save :: diag_initialized = .False.

    ! Persistent work arrays reused across Davidson iterations to avoid repeated allocation.
    Real(type_real), Allocatable, Save :: col_in_cache(:), col_out_cache(:), C_dvdsn(:)
    ! b and M^-1 b on the Z1 region for the Olsen correction in Dvdsn
    Real(type_real), Allocatable, Save :: Cb_z1(:), Bb_z1(:)
    ! Olsen is used for levels with more than this fraction of b outside the Z1 region
    Real(dp), Parameter :: olsen_weight = 0.5d0
    ! Olsen is skipped when b is nearly orthogonal to M^-1 b
    Real(dp), Parameter :: olsen_min = 1.d-3

    ! H times the current eigenvectors (local rows), kept between iterations so it can be
    ! rotated along with the vectors instead of recomputed by Mxmpy(1).
    ! hb_valid is true while HB matches ArrB(:,1:Nlv); Ortn then applies its changes to HB too.
    Real(type_real), Allocatable, Save :: HB(:,:)
    Logical, Save :: hb_valid = .False.

    ! Layout of the ArrB columns during Davidson (b = vectors, c = corrections, q = previous steps):
    !   1..Nlv          b
    !   Nlv+1..2Nlv     H*b, overwritten by c in Dvdsn
    !   2Nlv+1..3Nlv    q
    !   3Nlv+1..4Nlv    H*c; also scratch for Prj_J, and for the new b in NewVectors
    !   4Nlv+1..5Nlv    H*q
    ! The subspace is [b c q], so P is 3Nlv x 3Nlv.
    ! q is the last step of each unconverged level; keeping it makes the iteration
    ! conjugate-gradient-like (LOBPCG) instead of steepest descent, which zigzags on slow levels.
    ! q_active(i) is true when level i has a previous step q in the subspace.
    Logical, Allocatable, Save :: q_active(:)

    Public :: Init4, FormB0, FormB, FormBskip, WriteXIJ, FormP, AvgDiag, Dvdsn, Mxmpy, Ortn, Prj_J
    Public :: SaveHb, LoadHb, RotateHb, InvalidateHb, FreeHb, PrintCnorm, GatherColumnToRoot
    Public :: ResetPrevStep, KeepPrevStep, OrtnPrevStep, SetupH2D

  Contains
    
    Subroutine Init4(mype, mpi_type_real)
        ! This subroutine constructs the initial approximation 
        ! by selecting configurations specified by parameter Nc4.
        ! The initial approximation Hamiltonian is stored in the 
        ! matrix Z1 and is constructed by selecting the top-left 
        ! block of the full Hamiltonian matrix H.
        Use mpi_f08
        Use determinants, Only : calcNd0
        Use str_fmt, Only : FormattedTime, startTimer, stopTimer
        Implicit None

        Type(MPI_Datatype) :: mpi_type_real
        Integer  :: k, n, mype, mpierr
        Integer(Kind=int64) :: l8, s1
        Real(type_real) :: t
        Character(Len=256) :: strfmt
        Character(Len=16) :: timeStr

        If (mype == 0) Then
            Call startTimer(s1)
            strfmt = '(3X,"Starting approx. includes ",I5," conf.,",I6," det.")'
            Write( 6,strfmt) Nc1,Nd0
            Write(11,strfmt) Nc1,Nd0
        End If

        ! Construct the Hamiltonian in the initial approximation
        Diag=0_type_real
        diag_initialized = .False.  ! Reset so Mxmpy will re-fill local Diag on first call
        Z1=0_type_real
        Do n=nd_start+1, nd_end
            If (n > Nd0) Then
                If (Kl /= 1 .And. Kl /= 3) Exit
                Cycle
            End If
            Do l8=HamilCSR_rowptr(n-nd_start-1)+1_int64, HamilCSR_rowptr(n-nd_start)
                k=Hamil%col(l8)
                t=Hamil%val(l8)
                If (n == k) t=t-Hamil%minval
                If (k <= Nd0) Then
                    Z1(n,k)=t
                    Z1(k,n)=t
                End If
            End Do
        End Do

        Call MPI_AllReduce(MPI_IN_PLACE, Z1, Nd0*Nd0, mpi_type_real, MPI_SUM, &
                                MPI_COMM_WORLD, mpierr)

        If (mype == 0) Then
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A)') 'TIMING >>> Initial approximation constructed in '// trim(timeStr)
        End If

    End Subroutine Init4

    Subroutine FormB0(mype, mpi_type_real)
        ! This subroutine constructs the initial eigenvectors for the Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are written to CONF.XIJ 
        Use mpi_f08
        Use formj2, Only : J_av
        Use str_fmt, Only : FormattedTime, startTimer, stopTimer
        Implicit None

        Type(MPI_Datatype) :: mpi_type_real
        Integer :: i, j, lr, idum, ndpt, ierr, num, nskip, mpierr, mype, i_global
        Integer(kind=int64) :: s1
        Real(dp) :: dummy
        Real(type_real) :: xj
        Real(type_real), Allocatable :: B1_tmp(:)  ! temporary full-Nd buffer for file-read path
        Character(Len=128) :: strfmt
        Character(Len=16) :: timeStr

        num=0
        nskip=0

        If (mype==0) Then
            If (kCSF == 0) Then
                Open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')
            Else
                Open(unit=17,file='CONF.WFS',status='UNKNOWN',form='UNFORMATTED')
            End If
            Call startTimer(s1)
        End If

        !Call MPI_Bcast(Z1, Nd0*Nd0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        B1=0_type_real
        If (abs(Kl4) /= 2) Then ! If not reading CONF.XIJ
            Do j=1,Nd0
                B1(1:Nd0)=Z1(1:Nd0,j)
                If (kCSF == 0) Then
                    Call J_av(B1,Nd0,xj,mpi_type_real,ierr)
                Else
                    ierr = 0
                    xj = Jm
                End If
                If (ierr == 0) Then
                    num=num+1
                    E(num)=-(E1(j)+Hamil%minval)
                    Tj(num)=xj
                    ! B1(1:Nd0) is the eigenvector; rows > Nd0 are zero (initial trial vectors).
                    Do lr = 1, nd_local
                        i_global = nd_start + lr
                        If (i_global <= Nd0) Then
                            ArrB(lr, num) = B1(i_global)
                        Else
                            ArrB(lr, num) = 0_type_real
                        End If
                    End Do
                    If (num >= Nlv) Exit
                Else
                    nskip=nskip+1
                    If (mype==0) Write(*,*) '  skip Ej =',-(E1(j)+Hamil%minval)
                End If
            End Do
        Else ! Reading CONF.XIJ
            If (mype == 0) Then
                Read(17,err=210) dummy,dummy,ndpt
 210            Write(*,*) ' Vector length = ',ndpt
                Rewind(17)
            End If
            ! B1 is now only Nd0 elements; need a full Nd buffer for the file read and Bcast.
            Allocate(B1_tmp(Nd))
            Do j=1,Nlv
                B1_tmp = 0_type_real
                If (mype == 0) Then
                    Read(17,end=220,err=220) E(j),Tj(j),idum,(B1_tmp(i),i=1,ndpt)
                    num=j
                End If
                Call MPI_Bcast(B1_tmp, Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
                ArrB(1:nd_local,j)=B1_tmp(nd_start+1:nd_end)
            End Do
            Call MPI_Bcast(num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
 220    if (mype==0) Rewind(17)
        If (Allocated(B1_tmp)) Deallocate(B1_tmp)
        Nlv=num

        ! Write initial eigenvectors to CONF.XIJ; full columns are collected on rank 0 only
        If (mype == 0) Then
            Do j=1,Nlv
                strfmt = '(1X,"E(",I3,") =",F12.6," Jtot =",F10.6)'
                Write( 6,strfmt) j,E(j),Tj(j)
                Write(11,strfmt) j,E(j),Tj(j)
            End Do
        End If
        If (KXIJ > 0) Then
            Block
                Real(type_real), Allocatable :: vec_gather(:)
                If (mype == 0) Then
                    Allocate(vec_gather(Nd))
                Else
                    Allocate(vec_gather(1))
                End If
                Do j=1,Nlv
                    Call GatherColumnToRoot(j, vec_gather, mpi_type_real)
                    If (mype==0) Write(17) E(j),Tj(j),Nd,(vec_gather(i),i=1,Nd)
                End Do
                Deallocate(vec_gather)
            End Block
        End If
        If (mype == 0) Then
            Close(unit=17)
            If (nskip > 0) Then
                Write(*,*) ' Selection is done', nskip,' vectors skiped'
                Write(11,*) ' Selection is done', nskip,' vectors skiped'
            End If
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A)') 'TIMING >>> Initial eigenvectors constructed in '// trim(timeStr) // '.'
        End If
        
        Return
    End Subroutine FormB0

    Subroutine GatherColumnToRoot(icol, buf, mpi_type_real)
        ! Collects column icol of the distributed ArrB into buf on rank 0 only.
        ! Used where only rank 0 needs the full vector, to write it to a file.
        ! buf must hold Nd values on rank 0; on other ranks it is not touched.
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: icol
        Real(type_real), Intent(InOut) :: buf(:)
        Type(MPI_Datatype), Intent(In) :: mpi_type_real
        Integer :: mpierr

        Call MPI_Gatherv(ArrB(1:nd_local, icol), nd_local, mpi_type_real, &
                         buf, nd_counts_all, nd_displs_all, mpi_type_real, &
                         0, MPI_COMM_WORLD, mpierr)
    End Subroutine GatherColumnToRoot

    Subroutine NewVectors
        ! Form the next vectors from the eigenvectors Y of P (first Nlv columns of P):
        ! new q = c*Y_c + q*Y_q goes to the q block,
        ! new b = b*Y_b + new q goes to block 4 (H*c is no longer needed after RotateHb).
        ! The caller copies block 4 into the b block.
        Implicit None

        Select Case(type_real)
        Case(sp)
            Call SGEMM('N','N', nd_local, Nlv, 2*Nlv, 1.0_sp, ArrB(1,Nlv+1), nd_local, &
                       P(Nlv+1,1), 3*Nlv, 0.0_sp, ArrB(1,3*Nlv+1), nd_local)
        Case(dp)
            Call DGEMM('N','N', nd_local, Nlv, 2*Nlv, 1.0_dp, ArrB(1,Nlv+1), nd_local, &
                       P(Nlv+1,1), 3*Nlv, 0.0_dp, ArrB(1,3*Nlv+1), nd_local)
        End Select
        ArrB(1:nd_local,2*Nlv+1:3*Nlv) = ArrB(1:nd_local,3*Nlv+1:4*Nlv)
        Select Case(type_real)
        Case(sp)
            Call SGEMM('N','N', nd_local, Nlv, Nlv, 1.0_sp, ArrB(1,1), nd_local, &
                       P(1,1), 3*Nlv, 1.0_sp, ArrB(1,3*Nlv+1), nd_local)
        Case(dp)
            Call DGEMM('N','N', nd_local, Nlv, Nlv, 1.0_dp, ArrB(1,1), nd_local, &
                       P(1,1), 3*Nlv, 1.0_dp, ArrB(1,3*Nlv+1), nd_local)
        End Select
    End Subroutine NewVectors

    Subroutine FormBskip(mype)
        ! Forms eigenvectors for next Davidson iteration; not written to CONF.XIJ.
        Implicit None
        Integer, Intent(In) :: mype
        Integer :: k

        Call NewVectors

        Do k=1,Nlv
            Tk(k) = -(E(k)+Hamil%minval)
            ArrB(1:nd_local,k) = ArrB(1:nd_local,3*Nlv+k)
        End Do

        If (mype==0) Write ( 6,*) ' FormBskip: Vectors not saved'
        Return
    End Subroutine FormBskip

    Subroutine FormB(mpi_type_real, mype)
        ! Forms eigenvectors for next Davidson iteration; written to CONF.XIJ.
        Use mpi_f08
        Implicit None
        Type(MPI_Datatype), Intent(In) :: mpi_type_real
        Integer, Intent(In) :: mype
        Integer :: i, k
        Real(type_real) :: t
        Real(type_real), Allocatable :: vec_gather(:)

        Call NewVectors

        ! Only rank 0 writes the file, so only rank 0 holds the full vector
        If (mype == 0) Then
            Allocate(vec_gather(Nd))
        Else
            Allocate(vec_gather(1))
        End If
        If (mype==0) Then
            If (kCSF == 0) Then
                Open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')
            Else
                Open(unit=17,file='CONF.WFS',status='UNKNOWN',form='UNFORMATTED')
            End If
        End If

        Do k=1,Nlv
            t = -(E(k)+Hamil%minval)
            Tk(k) = t
            Call GatherColumnToRoot(k+3*Nlv, vec_gather, mpi_type_real)
            If (mype==0) Write(17) t,Tj(k),Nd,(vec_gather(i),i=1,Nd)
        End Do

        Do k=1,Nlv
            ArrB(1:nd_local,k) = ArrB(1:nd_local,3*Nlv+k)
        End Do

        Deallocate(vec_gather)
        If (mype==0) Then
            Close(unit=17)
            Write ( 6,*) ' FormB: Vectors saved'
        End If
        Return
    End Subroutine FormB

    Subroutine WriteXIJ(mpi_type_real, mype)
        ! Writes eigenvectors to CONF.XIJ; gathers distributed ArrB before write.
        Use mpi_f08
        Implicit None
        Type(MPI_Datatype), Intent(In) :: mpi_type_real
        Integer, Intent(In) :: mype
        Integer :: i, k, mpierr
        Real(type_real), Allocatable :: vec_gather(:)

        Allocate(vec_gather(Nd))
        If (mype==0) Then
            If (kCSF == 0) Then
                Open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')
            Else
                Open(unit=17,file='CONF.WFS',status='UNKNOWN',form='UNFORMATTED')
            End If
        End If
        Do k=1,Nlv
            vec_gather(nd_start+1:nd_end) = ArrB(1:nd_local, k)
            Call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                                vec_gather, nd_counts_all, nd_displs_all, &
                                mpi_type_real, MPI_COMM_WORLD, mpierr)
            If (mype==0) Write(17) Tk(k),Tj(k),Nd,(vec_gather(i),i=1,Nd)
        End Do
        If (mype==0) Then
            Close(17)
            Write ( 6,*) ' WriteXIJ: Vectors saved'
        End If
        Deallocate(vec_gather)
    End Subroutine WriteXIJ

    Subroutine FormP(kp, vmax, mpi_type_real)
        ! Calculation of upper-left block (Nlv X Nlv) of matrix P for kp=1
        ! Calculation of the rest of P (3Nlv X 3Nlv, subspace [b c q]) for kp=2
        ! For kp=2, rows of converged c and unused q are decoupled and shifted above the wanted levels.
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: kp
        Real(type_real), Intent(In) :: vmax
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: i, k, i1, mpierr
        Real(type_real), Allocatable :: Ptmp(:,:)

        If (kp == 1) Then
            ! Ptmp(i,k) = ArrB(:,i) * ArrB(:,k+Nlv) for i,k = 1...Nlv
            Allocate(Ptmp(Nlv, Nlv))
            Ptmp = 0_type_real
            Select Case(type_real)
            Case(sp)
                Call SGEMM('T','N', Nlv, Nlv, nd_local, 1.0_sp, ArrB(1,1), nd_local, &
                           ArrB(1,Nlv+1), nd_local, 0.0_sp, Ptmp, Nlv)
            Case(dp)
                Call DGEMM('T','N', Nlv, Nlv, nd_local, 1.0_dp, ArrB(1,1), nd_local, &
                           ArrB(1,Nlv+1), nd_local, 0.0_dp, Ptmp, Nlv)
            End Select
            Call MPI_Allreduce(MPI_IN_PLACE, Ptmp, Nlv*Nlv, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Do i=1,Nlv
                Do k=1,i
                    P(i,k) = Ptmp(i,k)
                    P(k,i) = Ptmp(i,k)
                End Do
            End Do
            Deallocate(Ptmp)
        Else
            ! Ptmp(i,k) = ArrB(:,i) * ArrB(:,k+3*Nlv) for i=1...3*Nlv, k=1...2*Nlv,
            ! i.e. the subspace [b c q] times [H*c H*q]
            Allocate(Ptmp(3*Nlv, 2*Nlv))
            Ptmp = 0_type_real
            Select Case(type_real)
            Case(sp)
                Call SGEMM('T','N', 3*Nlv, 2*Nlv, nd_local, 1.0_sp, ArrB(1,1), nd_local, &
                           ArrB(1,3*Nlv+1), nd_local, 0.0_sp, Ptmp, 3*Nlv)
            Case(dp)
                Call DGEMM('T','N', 3*Nlv, 2*Nlv, nd_local, 1.0_dp, ArrB(1,1), nd_local, &
                           ArrB(1,3*Nlv+1), nd_local, 0.0_dp, Ptmp, 3*Nlv)
            End Select
            Call MPI_Allreduce(MPI_IN_PLACE, Ptmp, 6*Nlv*Nlv, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
            ! b-c block
            Do i=1,Nlv
                Do k=Nlv+1,2*Nlv
                    P(i,k) = Ptmp(i, k-Nlv)
                    P(k,i) = Ptmp(i, k-Nlv)
                End Do
            End Do
            ! c-c block
            Do i=Nlv+1,2*Nlv
                Do k=Nlv+1,i
                    P(i,k) = Ptmp(i, k-Nlv)
                    P(k,i) = Ptmp(i, k-Nlv)
                End Do
            End Do
            ! b-q and c-q blocks
            Do i=1,2*Nlv
                Do k=2*Nlv+1,3*Nlv
                    P(i,k) = Ptmp(i, k-Nlv)
                    P(k,i) = Ptmp(i, k-Nlv)
                End Do
            End Do
            ! q-q block
            Do i=2*Nlv+1,3*Nlv
                Do k=2*Nlv+1,i
                    P(i,k) = Ptmp(i, k-Nlv)
                    P(k,i) = Ptmp(i, k-Nlv)
                End Do
            End Do
            ! Converged levels have no c, and unused q columns are zero vectors.
            ! Their rows are decoupled and moved above all wanted eigenvalues.
            Do i=1,Nlv
                If (Iconverge(i) == 1) Then
                    i1=i+Nlv
                    P(i1,1:3*Nlv) = 0_type_real
                    P(1:3*Nlv,i1) = 0_type_real
                    P(i1,i1) = 2*vmax - P(i,i) + 1.d0
                End If
                If (.not. q_active(i)) Then
                    i1=i+2*Nlv
                    P(i1,1:3*Nlv) = 0_type_real
                    P(1:3*Nlv,i1) = 0_type_real
                    P(i1,i1) = 2*vmax - P(i,i) + 1.d0
                End If
            End Do
            Deallocate(Ptmp)
        End If
        Return
    End Subroutine FormP
    
    Subroutine AvgDiag(mpi_type_real, mype)
        ! Average the diagonal over relativistic configurations.
        Use mpi_f08
        Implicit None
        Type(MPI_Datatype), Intent(In) :: mpi_type_real
        Integer, Intent(In) :: mype

        Integer :: ic, id, id0, id1, id2, loc_id0, loc_id2, mpierr
        Real(type_real), Allocatable :: group_sums(:)
        Integer, Allocatable :: group_counts(:)

        Allocate(group_sums(Nc), group_counts(Nc))
        group_sums = 0_type_real

        ! Accumulate each rank's local contribution to each configuration group
        id0 = 1
        id2 = 0
        Do ic = 1, Nc
            If (kCSF > 0) Then
                id1 = ndcs(nc_neq(ic))
            Else
                id1 = Ndc(ic)
            End If
            group_counts(ic) = id1
            id2 = id0 + id1 - 1
            If (id1 > 0) Then
                loc_id0 = Max(id0, nd_start + 1)
                loc_id2 = Min(id2, nd_end)
                Do id = loc_id0, loc_id2
                    group_sums(ic) = group_sums(ic) + Diag(id - nd_start)
                End Do
                id0 = id0 + id1
            End If
        End Do

        ! Collect group sums from all ranks
        Call MPI_Allreduce(MPI_IN_PLACE, group_sums, Nc, mpi_type_real, MPI_SUM, &
                           MPI_COMM_WORLD, mpierr)

        ! Each rank writes the group average back to its local Diag slice
        id0 = 1
        Do ic = 1, Nc
            id1 = group_counts(ic)
            id2 = id0 + id1 - 1
            If (id1 > 0) Then
                loc_id0 = Max(id0, nd_start + 1)
                loc_id2 = Min(id2, nd_end)
                If (loc_id0 <= loc_id2) &
                    Diag(loc_id0 - nd_start : loc_id2 - nd_start) = group_sums(ic) / id1
                id0 = id0 + id1
            End If
        End Do

        If (mype == 0) Then
            If (id2 == Nd) Then
                Write(*,*) ' Diagonal averaged over rel. configurations'
            Else
                Write(*,*) ' Error: id2=',id2,' Nc=',Nc
                Stop
            End If
        End If

        Deallocate(group_sums, group_counts)
    End Subroutine AvgDiag

    Subroutine Dvdsn(j, cnx, cnorm, mpi_type_real, mype)
        ! One Davidson correction step for level j.
        !
        ! The Davidson method iteratively improves eigenvector approximations by adding
        ! correction vectors to the trial subspace. Each call processes one level j:
        !
        !   1. Compute the residual C = ArrB(:,j1) - val*ArrB(:,j), where ArrB(:,j)
        !      is the current eigenvector approximation and val = P(j,j) is its eigenvalue.
        !      C measures how far the pair (val, ArrB(:,j)) is from satisfying H*x = E*x;
        !      small ||C|| means the level has converged.
        !
        !   2. Build the correction vector (new trial vector) using a preconditioner M^-1.
        !      For determinants in the initial-approximation subspace (rows 1...Nd0), a
        !      projection preconditioner based on Z1 (the diagonalized initial guess) is used.
        !      For all other determinants, a simple diagonal preconditioner 1/(val - Diag(i))
        !      is applied. The correction is stored in ArrB(:,j1) for use by FormP/FormB.
        !
        !   3. Olsen correction, only for levels that lie mainly outside the Z1 region:
        !      correction = M^-1 C - eps * M^-1 b, with eps chosen so that it is orthogonal to b.
        !      There the diagonal preconditioner is nearly exact on the level's own main
        !      determinants (val close to their diagonal), so M^-1 C points almost along b;
        !      Ortn then removes that part and leaves mostly noise, and the level stalls.
        !      Inside the Z1 region this is not needed: the lowest Nlv Z1 columns, the levels
        !      themselves, are already left out of the preconditioner.
        !
        ! ArrB is distributed across ranks (each owns rows nd_start+1...nd_end).
        ! B1(1:Nd0) accumulates the Z1-region correction and is replicated on all ranks.
        ! B1_loc(1:nd_local) holds the diagonal-preconditioned correction for local rows > Nd0.
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: j, mype
        Real(dp), Intent(InOut) :: cnx
        ! Residual norm ||C|| of level j; printed by the caller with PrintCnorm
        Real(type_real), Intent(Out) :: cnorm
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: j1, i, i_global, lr, n01, nz, mpierr
        Real(type_real) :: val, t, s, s_Nd0, s_local
        Real(dp) :: w_out

        ! C_dvdsn(1:Nd0) holds the residual C for the Z1 region; C_dvdsn(Nd0+1) accumulates
        ! cnorm^2, allowing both to be Allreduced in a single call.
        If (.not. Allocated(C_dvdsn)) Allocate(C_dvdsn(Nd0+1))

        C_dvdsn = 0_type_real
        B1=0_type_real        ! Z1-region correction accumulator, replicated on all ranks
        B1_loc=0_type_real    ! diagonal-preconditioned correction for local rows > Nd0
        j1=j+Nlv
        val=P(j,j)

        ! Compute residual t = ArrB(:,j1) - val*ArrB(:,j) for each local row.
        ! Rows in the Z1 region (i_global <= Nd0) are stored in C_dvdsn for the Z1 preconditioner.
        ! Rows outside the Z1 region get the diagonal preconditioner applied immediately: t/(val-Diag).
        Do lr=1,nd_local
            i_global = nd_start + lr
            t = ArrB(lr,j1) - val*ArrB(lr,j)
            If (i_global <= Nd0) C_dvdsn(i_global) = t
            C_dvdsn(Nd0+1) = C_dvdsn(Nd0+1) + t*t
            If (i_global > Nd0) B1_loc(lr) = t/(val-Diag(lr))
        End Do

        ! Allreduce assembles the full residual C and cnorm^2 across all ranks in one call.
        Call MPI_Allreduce(MPI_IN_PLACE, C_dvdsn, Nd0+1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
        cnorm = dsqrt(dabs(C_dvdsn(Nd0+1)))

        Iconverge(j)=0
        If (cnorm < Crt4) Iconverge(j)=1
        If (cnorm > cnx) cnx=Real(cnorm,kind=dp)
        If (Iconverge(j)==1) Then
            ArrB(1:nd_local,J1)=0_type_real
            Return
        End If

        ! Z1 projection preconditioner for the initial-approximation subspace (rows 1...Nd0).
        ! Z1 contains the eigenvectors of the initial Hamiltonian approximation with eigenvalues E1.
        ! The correction is B1 += Z1 * diag(1/(val-E1)) * Z1^T * C, which shifts the residual
        ! into the Z1 eigenbasis, scales each component by 1/(val-E1), and transforms back.
        ! This is more effective than a plain diagonal preconditioner in the Z1 subspace.
        ! B1 is replicated on all ranks since Z1 and C_dvdsn(1:Nd0) are both replicated.
        n01=Nlv+1
        If (iabs(Kl4) == 2) n01=1
        nz = Nd0 - n01 + 1
        Call PrecondZ1(C_dvdsn(1:Nd0), B1)

        ! Weight of b outside the Z1 region (b is normalized); decides whether Olsen is needed
        w_out = 0.d0
        Do lr=1,nd_local
            If (nd_start+lr > Nd0) w_out = w_out + Real(ArrB(lr,j),dp)**2
        End Do
        Call MPI_Allreduce(MPI_IN_PLACE, w_out, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        If (w_out > olsen_weight) Call OlsenCorrection

        ! Normalize the correction vector and store it in ArrB(:,j1).
        ! B1(1:Nd0) is replicated (no Allreduce needed); B1_loc is local and requires Allreduce.
        s_Nd0 = 0_type_real
        Do i=1,Nd0
            s_Nd0 = s_Nd0 + B1(i)**2
        End Do
        s_local = sum(B1_loc(1:nd_local)**2)
        Call MPI_Allreduce(MPI_IN_PLACE, s_local, 1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
        s = 1.d0/dsqrt(s_Nd0 + s_local)

        ! Write normalized correction vector into ArrB(:,j1), partitioned by region.
        Do lr=1,nd_local
            i_global = nd_start + lr
            If (i_global <= Nd0) Then
                ArrB(lr,J1) = B1(i_global) * s
            Else
                ArrB(lr,J1) = B1_loc(lr) * s
            End If
        End Do

        Return

    Contains

        Subroutine PrecondZ1(vin, vout)
            ! vout += Z1 * diag(1/(val-E1)) * Z1^T * vin over columns n01..Nd0 of Z1.
            ! Terms with |val-E1| below 1e-6 are dropped.
            Real(type_real), Intent(In) :: vin(Nd0)
            Real(type_real), Intent(InOut) :: vout(Nd0)
            Real(type_real), Allocatable :: s_vec(:)
            Real(type_real) :: den
            Integer :: k

            If (nz <= 0) Return
            Allocate(s_vec(nz))
            ! Project onto the Z1 eigenbasis: s_vec = Z1^T * vin
            Select Case(type_real)
            Case(sp)
                Call SGEMV('T', Nd0, nz, 1.0_sp, Z1(1,n01), Nd0, vin, 1, 0.0_sp, s_vec, 1)
            Case(dp)
                Call DGEMV('T', Nd0, nz, 1.0_dp, Z1(1,n01), Nd0, vin, 1, 0.0_dp, s_vec, 1)
            End Select
            ! Scale each component by 1/(val - E1); skip near-zero denominators
            Do k=1,nz
                den = val - E1(n01+k-1)
                If (dabs(den) < 1.d-6) Then
                    s_vec(k) = 0_type_real
                Else
                    s_vec(k) = s_vec(k) / den
                End If
            End Do
            ! Transform back to the determinant basis and accumulate into vout
            Select Case(type_real)
            Case(sp)
                Call SGEMV('N', Nd0, nz, 1.0_sp, Z1(1,n01), Nd0, s_vec, 1, 1.0_sp, vout, 1)
            Case(dp)
                Call DGEMV('N', Nd0, nz, 1.0_dp, Z1(1,n01), Nd0, s_vec, 1, 1.0_dp, vout, 1)
            End Select
            Deallocate(s_vec)
        End Subroutine PrecondZ1

        Subroutine OlsenCorrection
            ! eps = (b . M^-1 C) / (b . M^-1 b);  B1, B1_loc -= eps * M^-1 b.
            ! M^-1 b is Bb_z1 on the Z1 region (replicated) and B2 on the local rows > Nd0.
            Real(dp) :: eps, dots(3)

            If (.not. Allocated(Cb_z1)) Allocate(Cb_z1(Nd0), Bb_z1(Nd0))
            Cb_z1 = 0_type_real
            Bb_z1 = 0_type_real
            B2 = 0_type_real
            Do lr=1,nd_local
                i_global = nd_start + lr
                If (i_global <= Nd0) Then
                    Cb_z1(i_global) = ArrB(lr,j)
                Else
                    B2(lr) = ArrB(lr,j)/(val-Diag(lr))
                End If
            End Do
            Call MPI_Allreduce(MPI_IN_PLACE, Cb_z1, Nd0, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call PrecondZ1(Cb_z1, Bb_z1)

            ! The Z1-region parts are replicated, so only the local parts need a reduction
            dots(1) = Sum(Real(B1_loc(1:nd_local),dp) * Real(ArrB(1:nd_local,j),dp))
            dots(2) = Sum(Real(B2(1:nd_local),dp) * Real(ArrB(1:nd_local,j),dp))
            dots(3) = Sum(Real(B2(1:nd_local),dp)**2)
            Call MPI_Allreduce(MPI_IN_PLACE, dots, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            dots(1) = dots(1) + Sum(Real(Cb_z1,dp) * Real(B1(1:Nd0),dp))
            dots(2) = dots(2) + Sum(Real(Cb_z1,dp) * Real(Bb_z1,dp))
            dots(3) = dots(3) + Sum(Real(Bb_z1,dp)**2)

            ! Skip when b is nearly orthogonal to M^-1 b, since eps is then ill-defined
            If (dabs(dots(2)) <= olsen_min * dsqrt(dots(3))) Return
            eps = dots(1) / dots(2)
            B1(1:Nd0) = B1(1:Nd0) - Real(eps,type_real) * Bb_z1
            B1_loc(1:nd_local) = B1_loc(1:nd_local) - Real(eps,type_real) * B2(1:nd_local)
        End Subroutine OlsenCorrection

    End Subroutine Dvdsn

    Subroutine PrintCnorm(j, cnorm, mype)
        ! Prints the residual norm of level j from Dvdsn.
        Implicit None
        Integer, Intent(In) :: j, mype
        Real(type_real), Intent(In) :: cnorm
        Character(Len=9) :: tag

        If (mype /= 0) Return
        tag='         '
        If (Iconverge(j) == 1) Then
            If (kdavidson == 1) Then
                tag='Stopped  '
            Else
                tag='converged'
            End If
        End If
        Write ( 6,'(4X,"|| C(",I3,") || = ",F10.7,2X,A9)') j,cnorm,tag
        Write (11,'(4X,"|| C(",I3,") || = ",F10.7,2X,A9)') j,cnorm,tag
    End Subroutine PrintCnorm

    Subroutine Mxmpy(ip, mpi_type_real, mype)
        ! Multiply the Hamiltonian H by a block of trial vectors: H * ArrB(:, i1min:i1min+Nlv-1)
        ! and store the result in ArrB(:, i2min:i2min+Nlv-1).
        ! ip=1: H*b, ip=2: H*c, ip=3: H*q (see the ArrB layout at the top of the module).
        ! This is the core operation of the Davidson method. 
        ! It applies H to the current trial vectors to build the projected matrix 
        ! P = ArrB^T * H * ArrB, from which improved eigenvalue estimates are obtained.
        !
        ! The Hamiltonian is stored in CSR format and distributed across ranks, so H*v requires
        ! communication. For each trial vector column:
        !   1. Each rank contributes its local rows of the input vector; Allgatherv assembles
        !      the full input column on all ranks (needed because off-diagonal H elements can
        !      couple any two determinants).
        !   2. Each rank computes H*v for its own row range using the local CSR rows of H.
        !      The upper-triangle-only storage means row n touching column k also contributes
        !      a symmetric k->n term, so both col_out(n) and col_out(k) are updated.
        !   3. Reduce_scatter sums these scatter contributions across ranks and delivers to
        !      each rank only the output values for its own rows.
        !   4. The local slice lands directly in ArrB.
        !
        ! On the first call (kd=1), the diagonal elements of H (shifted by Hamil%minval) are
        ! extracted into Diag for use as the preconditioner in Dvdsn.
        Use mpi_f08
        Implicit None

        Type(MPI_Datatype) :: mpi_type_real
        Integer, Intent(In) :: ip, mype

        Integer :: iv, k, n, mpierr, i2min, i1min, kd
        Real(type_real) :: t
        Integer(Kind=int64) :: l8
        ! Clock for the optional phase timers (PCI_MATVEC_TIMING)
        Real(dp) :: t0

        Select Case(ip)
        Case(1)
            i1min = 1
            i2min = Nlv+1
        Case(2)
            i1min = Nlv+1
            i2min = 3*Nlv+1
        Case Default
            i1min = 2*Nlv+1
            i2min = 4*Nlv+1
        End Select
        kd = 0

        If (.not. diag_initialized) kd=1

        ! The block-partitioned product (KHB=1) uses its own block-sized buffers
        If (.not. h2d_active .and. (.not. Allocated(col_in_cache) .or. size(col_in_cache) /= Nd)) Then
            If (Allocated(col_in_cache))  Deallocate(col_in_cache)
            If (Allocated(col_out_cache)) Deallocate(col_out_cache)
            Allocate(col_in_cache(Nd), col_out_cache(Nd))
        End If

        Do iv = 0, Nlv-1
            ! A converged level has a zero correction vector, so its product is zero.
            ! Iconverge is the same on all ranks, so every rank skips the same collectives.
            If (ip == 2) Then
                If (Iconverge(iv+1) == 1) Then
                    ArrB(1:nd_local, i2min+iv) = Real(0, kind=type_real)
                    Cycle
                End If
            End If
            ! Unused previous steps are zero too.
            If (ip == 3) Then
                If (.not. q_active(iv+1)) Then
                    ArrB(1:nd_local, i2min+iv) = Real(0, kind=type_real)
                    Cycle
                End If
            End If

            If (h2d_active) Then
                Call MatvecH2D(ArrB(1:nd_local, i1min+iv), ArrB(1:nd_local, i2min+iv), mpi_type_real)
                Cycle
            End If

            Call MvCount
            Call MvWait(1, t0)

            ! Step 1: assemble the full input column on all ranks.
            col_in_cache(nd_start+1:nd_end) = ArrB(1:nd_local, i1min+iv)
            Call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                                col_in_cache, nd_counts_all, nd_displs_all, &
                                mpi_type_real, MPI_COMM_WORLD, mpierr)
            Call MvMark(2, t0)

            ! Step 2: SpMV over local rows. Upper-triangle CSR storage means element (n,k)
            ! with n/=k contributes to both output row n and output row k (symmetric scatter).
            col_out_cache(:) = Real(0, kind=type_real)
            Do n = nd_start+1, nd_end
                Do l8 = HamilCSR_rowptr(n-nd_start-1)+1_int64, HamilCSR_rowptr(n-nd_start)
                    k = Hamil%col(l8)
                    t = Hamil%val(l8)
                    If (n == k) Then
                        t = t - Hamil%minval
                        If (kd == 1 .And. iv == 0) Diag(n - nd_start) = t
                    End If
                    col_out_cache(n) = col_out_cache(n) + t * col_in_cache(k)
                    If (n /= k) col_out_cache(k) = col_out_cache(k) + t * col_in_cache(n)
                End Do
            End Do
            Call MvMark(4, t0)
            Call MvWait(5, t0)

            ! Step 3-4: sum scatter contributions across ranks; each rank receives only its own
            ! row block (contiguous, in rank order) directly into its local slice of ArrB.
            Call MPI_Reduce_scatter(col_out_cache, ArrB(1:nd_local, i2min+iv), nd_counts_all, &
                                    mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MvMark(7, t0)
        End Do

        ! Each rank already wrote to its own Diag(n-nd_start); no Allreduce needed.
        If (kd==1) diag_initialized = .True.
        Return
    End Subroutine Mxmpy

    Subroutine SetupH2D(npes, mype, mpi_type_real)
        ! Switch H to the block layout of module hblocks (KHB=1); call after Init4.
        ! Diag is filled here from the row-distributed H, as the first 1D Mxmpy call would do,
        ! since that H is freed by BuildH2D.
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: npes, mype
        Type(MPI_Datatype), Intent(In) :: mpi_type_real
        Integer :: n
        Integer(Kind=int64) :: l8

        Do n = nd_start+1, nd_end
            Do l8 = HamilCSR_rowptr(n-nd_start-1)+1_int64, HamilCSR_rowptr(n-nd_start)
                If (Hamil%col(l8) == n) Diag(n - nd_start) = Hamil%val(l8) - Hamil%minval
            End Do
        End Do
        diag_initialized = .True.

        Call BuildH2D(npes, mype, mpi_type_real)
    End Subroutine SetupH2D

    Subroutine Ortn(j, ifail, mpi_type_real, mype)
        ! Orthogonalize trial vector j against the j-1 vectors already in the subspace,
        ! then normalize it. This ensures the trial vectors remain linearly independent,
        ! which is required for the Davidson subspace to be well-conditioned.
        !
        ! Uses classical Gram-Schmidt with re-orthogonalization: 
        ! in each pass, all j-1 dot products <ArrB(:,j), ArrB(:,l)> are computed in a single DGEMV call 
        ! and Allreduced once, then the projections are subtracted in a second DGEMV. 
        ! Passes repeat until the largest remaining overlap (smax1) drops below the threshold `ortho`, 
        ! or until the overlaps stop decreasing (indicating failure).
        ! Converged levels are excluded from orthogonalization, since their correction vectors
        ! are zeroed in Dvdsn and should not influence the remaining search).
        !
        ! After orthogonalization, the vector is normalized via Allreduce of its squared norm. 
        ! ifail is set to j if orthogonalization or normalization fails.
        Use mpi_f08
        Implicit None
        Integer, Intent(In)  :: j, mype
        Integer, Intent(Out) :: ifail
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: l, jm, it, mpierr
        Real(dp) :: smax2, smax1, critN, ortho
        Real(type_real) :: s
        Real(type_real), Allocatable :: s_arr(:)

        ortho=1.d-6
        critN=1.d-9
        jm=j-1
        smax2=1.d10
        ifail=0

        If (jm > 0) Allocate(s_arr(jm))

        it=0
        Do
            it=it+1
            smax1=0_type_real
            If (jm > 0) Then
                ! Compute all overlaps <ArrB(:,j), ArrB(:,l)> for l=1...jm in one DGEMV+Allreduce
                s_arr = 0_type_real
                Select Case(type_real)
                Case(sp)
                    Call SGEMV('T', nd_local, jm, 1.0_sp, ArrB(1,1), nd_local, &
                               ArrB(1,j), 1, 0.0_sp, s_arr, 1)
                Case(dp)
                    Call DGEMV('T', nd_local, jm, 1.0_dp, ArrB(1,1), nd_local, &
                               ArrB(1,j), 1, 0.0_dp, s_arr, 1)
                End Select
                Call MPI_Allreduce(MPI_IN_PLACE, s_arr, jm, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
                ! Skip converged levels; track largest remaining overlap for convergence check
                Do l=1,jm
                    If (l > Nlv .and. Iconverge(l-Nlv)==1) Then
                        s_arr(l) = 0_type_real
                    Else
                        If (dabs(s_arr(l)) > smax1) smax1 = dabs(s_arr(l))
                    End If
                End Do
                ! Subtract all projections at once: ArrB(:,j) -= ArrB(:,1:jm) * s_arr
                Select Case(type_real)
                Case(sp)
                    Call SGEMV('N', nd_local, jm, -1.0_sp, ArrB(1,1), nd_local, &
                               s_arr, 1, 1.0_sp, ArrB(1,j), 1)
                Case(dp)
                    Call DGEMV('N', nd_local, jm, -1.0_dp, ArrB(1,1), nd_local, &
                               s_arr, 1, 1.0_dp, ArrB(1,j), 1)
                End Select
                ! Keep H*b consistent: H is linear, so HB gets the same combination.
                If (hb_valid .and. j <= Nlv) Then
                    Select Case(type_real)
                    Case(sp)
                        Call SGEMV('N', nd_local, jm, -1.0_sp, HB(1,1), nd_local, &
                                   s_arr, 1, 1.0_sp, HB(1,j), 1)
                    Case(dp)
                        Call DGEMV('N', nd_local, jm, -1.0_dp, HB(1,1), nd_local, &
                                   s_arr, 1, 1.0_dp, HB(1,j), 1)
                    End Select
                End If
            End If
            If (smax1 < ortho) Exit
            If (smax1 > smax2 .or. it > 3) Then
                If (mype==0) Write (6,'(4X,"Fail of orthogonalization of vector ",I2)') j
                ifail=j
                If (jm > 0) Deallocate(s_arr)
                Return
            End If
            smax2=smax1
        End Do

        If (jm > 0) Deallocate(s_arr)

        s=0_type_real
        Do l=1,nd_local
            s=s+ArrB(l,j)*ArrB(l,j)
        End Do
        Call MPI_Allreduce(MPI_IN_PLACE, s, 1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
        If (s < critN) Then
            If (mype==0) Write (6,'(4X,"Fail of normalization of vector ",I2)') j
            ifail=j
        End If
        s=1.d0/dsqrt(s)
        ArrB(1:nd_local,j)=ArrB(1:nd_local,j)*s
        If (hb_valid .and. j <= Nlv) HB(1:nd_local,j)=HB(1:nd_local,j)*s
        Return
    End Subroutine Ortn

    Subroutine SaveHb
        ! Store H*b from Mxmpy(1), which is in ArrB(:,Nlv+1:2*Nlv), before Dvdsn overwrites it.
        Implicit None

        If (Allocated(HB)) Then
            If (size(HB,1) /= nd_local .or. size(HB,2) /= Nlv) Deallocate(HB)
        End If
        If (.not. Allocated(HB)) Allocate(HB(nd_local, Nlv))
        HB(1:nd_local,1:Nlv) = ArrB(1:nd_local,Nlv+1:2*Nlv)
        hb_valid = .True.
    End Subroutine SaveHb

    Subroutine LoadHb
        ! Put the stored H*b where Mxmpy(1) would have written it.
        Implicit None

        ArrB(1:nd_local,Nlv+1:2*Nlv) = HB(1:nd_local,1:Nlv)
    End Subroutine LoadHb

    Subroutine RotateHb
        ! Form H*b and H*q for the next iteration's vectors without multiplying by H.
        ! NewVectors makes q = c*Y_c + q*Y_q and b = b*Y_b + q,
        ! so the same combinations of H*c, H*q and H*b give H times the new vectors.
        ! H*c is in block 4 from Mxmpy(2); this must run before NewVectors overwrites it.
        ! Rows are done in chunks so that only a small scratch block is needed.
        Implicit None
        Integer, Parameter :: chunk = 4096
        Integer :: r0, nr
        Real(type_real), Allocatable :: hq(:,:), hb_new(:,:)

        Allocate(hq(Min(chunk, Max(nd_local,1)), Nlv), hb_new(Min(chunk, Max(nd_local,1)), Nlv))
        Do r0 = 1, nd_local, chunk
            nr = Min(chunk, nd_local - r0 + 1)
            ! H*q for the next iteration = [H*c H*q] * [Y_c; Y_q]
            Select Case(type_real)
            Case(sp)
                Call SGEMM('N','N', nr, Nlv, 2*Nlv, 1.0_sp, ArrB(r0,3*Nlv+1), nd_local, &
                           P(Nlv+1,1), 3*Nlv, 0.0_sp, hq, size(hq,1))
            Case(dp)
                Call DGEMM('N','N', nr, Nlv, 2*Nlv, 1.0_dp, ArrB(r0,3*Nlv+1), nd_local, &
                           P(Nlv+1,1), 3*Nlv, 0.0_dp, hq, size(hq,1))
            End Select
            ! H*b for the next iteration = H*b * Y_b + the H*q above
            hb_new(1:nr,1:Nlv) = hq(1:nr,1:Nlv)
            Select Case(type_real)
            Case(sp)
                Call SGEMM('N','N', nr, Nlv, Nlv, 1.0_sp, HB(r0,1), nd_local, &
                           P(1,1), 3*Nlv, 1.0_sp, hb_new, size(hb_new,1))
            Case(dp)
                Call DGEMM('N','N', nr, Nlv, Nlv, 1.0_dp, HB(r0,1), nd_local, &
                           P(1,1), 3*Nlv, 1.0_dp, hb_new, size(hb_new,1))
            End Select
            HB(r0:r0+nr-1,1:Nlv) = hb_new(1:nr,1:Nlv)
            ArrB(r0:r0+nr-1,4*Nlv+1:5*Nlv) = hq(1:nr,1:Nlv)
        End Do
        Deallocate(hq, hb_new)
    End Subroutine RotateHb

    Subroutine ResetPrevStep
        ! Drop all previous steps; the next step is a plain Davidson step.
        Implicit None

        If (.not. Allocated(q_active)) Allocate(q_active(Nlv))
        q_active = .False.
        ArrB(1:nd_local,2*Nlv+1:3*Nlv) = 0_type_real
        ArrB(1:nd_local,4*Nlv+1:5*Nlv) = 0_type_real
    End Subroutine ResetPrevStep

    Subroutine KeepPrevStep
        ! After the new vectors are formed, keep the previous step of unconverged levels only.
        ! Converged levels get no correction, so their step would only add noise to P.
        Implicit None
        Integer :: i

        Do i=1,Nlv
            q_active(i) = Iconverge(i) == 0
            If (.not. q_active(i)) Then
                ArrB(1:nd_local,2*Nlv+i) = 0_type_real
                ArrB(1:nd_local,4*Nlv+i) = 0_type_real
            End If
        End Do
    End Subroutine KeepPrevStep

    Subroutine OrtnPrevStep(mpi_type_real)
        ! Orthonormalize each previous step q against b, c and the earlier q,
        ! applying the same combination to H*q so that it stays H times q.
        ! Must run after Mxmpy(2), since it needs H*c, and with HB equal to H*b.
        ! A step that is almost inside the rest of the subspace is dropped:
        ! it adds nothing new, and normalizing it would magnify round-off in H*q.
        ! Classical Gram-Schmidt with two passes, as in Ortn.
        Use mpi_f08
        Implicit None
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: i, jq, jh, jm, pass, mpierr
        Real(type_real) :: s0, s1
        Real(dp) :: drop_tol
        Real(type_real), Allocatable :: s_arr(:)

        ! Smallest fraction of a step that must remain after orthogonalization
        If (type_real == sp) Then
            drop_tol = 1.d-2
        Else
            drop_tol = 1.d-4
        End If

        Allocate(s_arr(3*Nlv))
        Do i=1,Nlv
            If (.not. q_active(i)) Cycle
            jq = 2*Nlv+i
            jh = 4*Nlv+i
            jm = jq-1

            s0 = sum(ArrB(1:nd_local,jq)**2)
            Call MPI_Allreduce(MPI_IN_PLACE, s0, 1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)

            If (s0 > 0) Then
                Do pass=1,2
                    ! Overlaps with b, c and q(1:i-1) in one GEMV + Allreduce
                    Select Case(type_real)
                    Case(sp)
                        Call SGEMV('T', nd_local, jm, 1.0_sp, ArrB(1,1), nd_local, &
                                   ArrB(1,jq), 1, 0.0_sp, s_arr, 1)
                    Case(dp)
                        Call DGEMV('T', nd_local, jm, 1.0_dp, ArrB(1,1), nd_local, &
                                   ArrB(1,jq), 1, 0.0_dp, s_arr, 1)
                    End Select
                    Call MPI_Allreduce(MPI_IN_PLACE, s_arr, jm, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
                    ! q -= [b c q(1:i-1)] * s; H*q -= [H*b H*c H*q(1:i-1)] * s
                    Select Case(type_real)
                    Case(sp)
                        Call SGEMV('N', nd_local, jm, -1.0_sp, ArrB(1,1), nd_local, &
                                   s_arr, 1, 1.0_sp, ArrB(1,jq), 1)
                        Call SGEMV('N', nd_local, Nlv, -1.0_sp, HB(1,1), nd_local, &
                                   s_arr(1), 1, 1.0_sp, ArrB(1,jh), 1)
                        Call SGEMV('N', nd_local, Nlv, -1.0_sp, ArrB(1,3*Nlv+1), nd_local, &
                                   s_arr(Nlv+1), 1, 1.0_sp, ArrB(1,jh), 1)
                        If (i > 1) Call SGEMV('N', nd_local, i-1, -1.0_sp, ArrB(1,4*Nlv+1), nd_local, &
                                              s_arr(2*Nlv+1), 1, 1.0_sp, ArrB(1,jh), 1)
                    Case(dp)
                        Call DGEMV('N', nd_local, jm, -1.0_dp, ArrB(1,1), nd_local, &
                                   s_arr, 1, 1.0_dp, ArrB(1,jq), 1)
                        Call DGEMV('N', nd_local, Nlv, -1.0_dp, HB(1,1), nd_local, &
                                   s_arr(1), 1, 1.0_dp, ArrB(1,jh), 1)
                        Call DGEMV('N', nd_local, Nlv, -1.0_dp, ArrB(1,3*Nlv+1), nd_local, &
                                   s_arr(Nlv+1), 1, 1.0_dp, ArrB(1,jh), 1)
                        If (i > 1) Call DGEMV('N', nd_local, i-1, -1.0_dp, ArrB(1,4*Nlv+1), nd_local, &
                                              s_arr(2*Nlv+1), 1, 1.0_dp, ArrB(1,jh), 1)
                    End Select
                End Do
            End If

            s1 = sum(ArrB(1:nd_local,jq)**2)
            Call MPI_Allreduce(MPI_IN_PLACE, s1, 1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)

            If (s0 <= 0 .or. s1 <= (drop_tol**2)*s0) Then
                q_active(i) = .False.
                ArrB(1:nd_local,jq) = 0_type_real
                ArrB(1:nd_local,jh) = 0_type_real
            Else
                s1 = 1/Sqrt(s1)
                ArrB(1:nd_local,jq) = ArrB(1:nd_local,jq)*s1
                ArrB(1:nd_local,jh) = ArrB(1:nd_local,jh)*s1
            End If
        End Do
        Deallocate(s_arr)
    End Subroutine OrtnPrevStep

    Subroutine InvalidateHb
        ! Mark HB as out of date, so Ortn stops updating it until the next SaveHb.
        Implicit None

        hb_valid = .False.
    End Subroutine InvalidateHb

    Subroutine FreeHb
        Implicit None

        If (Allocated(HB)) Deallocate(HB)
        hb_valid = .False.
    End Subroutine FreeHb

    Subroutine Prj_J(lin, num, lout, trsd, mpi_type_real, mype)
        ! Filter `num` trial vectors (ArrB columns lin..lin+num-1) to have good J = XJ_av,
        ! storing the result in columns lout..lout+num-1.
        !
        ! The Hamiltonian commutes with J^2, so each eigenstate should have a definite J value.
        ! However, the Davidson trial vectors can mix different J values.
        ! This routine removes those unwanted J components by repeatedly applying the filter:
        !     (I - J^2/jx(jx+2)) for each unwanted J value jx in the list Jt(1:Njd).
        ! Each application attenuates components whose J^2 eigenvalue is jx(jx+2),
        ! i.e., belonging to angular momentum jx/2, while leaving the target J=XJ_av component unchanged. 
        ! After Njd+1 iterations, the vectors should lie almost entirely in the target J subspace.
        !
        ! Convergence is checked by computing <x|J^2|x> = sj and comparing
        !     aj = sqrt(1+sj) - 1 to the target jav = 2*XJ_av. 
        ! The iteration stops when the largest error over all vectors drops below the threshold `trsd`.
        !
        ! J^2 is applied as a symmetric SpMV using the redistributed CSR form (JsqCSR_rowptr + Jsq%col/val).
        Use mpi_f08
        Use determinants, Only : Gdet
        Implicit None

        Type(MPI_Datatype) :: mpi_type_real
        Integer :: ierr, i, n, k, j, i2, i1, jx, it, iter, lout, num, lin, jav, mpierr, imin, imax
        Integer, optional :: mype
        Integer(kind=int64) :: j8
        Real(dp) :: err1, err, f, trsd
        Real(type_real) :: t, sj, aj, s, s1, a
        logical :: proceed
        Character(Len=256) :: strfmt

        ierr=1
        If (mype == 0 .and. lout+num-1 > IPlv) Then
            strfmt = '(" Prj_J warning:", /," to project",I3," levels I need",I3," vectors & have only",I3)'
            Write( *,strfmt) Nlv,lout+num-1,IPlv
            Write(11,strfmt) Nlv,lout+num-1,IPlv
            Stop
        End If
        jav=Int(2*XJ_av+0.1d0)
        If (mype == 0) Then
            strfmt = '(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)'
            Write( *,strfmt) num,XJ_av
            Write(11,strfmt) num,XJ_av
        End If

        Call MPI_Bcast(Iconverge, Nlv, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        ! Gather input columns into full-size buffers; J^2 SpMV uses scatter (same as Mxmpy)
        ! so col_out is Allreduced and the local slice is written back to ArrB each iteration.
        Block
            Real(type_real), Allocatable :: col_in(:,:), col_out(:,:)
            Allocate(col_in(Nd,num), col_out(Nd,num))

            imin=lin
            imax=lin+num-1
            Do j=imin,imax
                i = j - lin + 1
                col_in(nd_start+1:nd_end, i) = ArrB(1:nd_local, j)
                Call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                                    col_in(:,i), nd_counts_all, nd_displs_all, &
                                    mpi_type_real, MPI_COMM_WORLD, mpierr)
            End Do

            Do iter=1,Njd+1
                it=Njd+1-iter
                If (it /= 0) Then
                    jx=Jt(it)
                    f=1.d0/(jx*(jx+2))
                End If
                If (jx /= jav .or. it == 0) Then
                    Do i=1,num
                        col_out(:,i)=0_type_real
                    End Do
                    Do n=1,nd_local
                        Do j8=JsqCSR_rowptr(n-1)+1_int64, JsqCSR_rowptr(n)
                            k=Jsq%col(j8)
                            t=Jsq%val(j8)
                            Do i=1,num
                                proceed=(lin-1)*Iconverge(i)==0
                                If (proceed) Then
                                    col_out(nd_start+n,i)=col_out(nd_start+n,i)+t*col_in(k,i)
                                    If (nd_start+n /= k) &
                                        col_out(k,i)=col_out(k,i)+t*col_in(nd_start+n,i)
                                End If
                            End Do
                        End Do
                    End Do

                    Do i=1,num
                        Call MPI_Allreduce(MPI_IN_PLACE, col_out(:,i), Nd, mpi_type_real, &
                                           MPI_SUM, MPI_COMM_WORLD, mpierr)
                        ArrB(1:nd_local, lout+i-1) = col_out(nd_start+1:nd_end, i)
                    End Do
                    err1=0_type_real

                    Do i=1,num
                        proceed=(lin-1)*Iconverge(i)==0
                        If (proceed) Then
                            sj=0.0
                            Do n=1,Nd
                                sj=sj+col_in(n,i)*col_out(n,i)
                            End Do
                            aj=dsqrt(1+sj) - 1
                            err=dabs(aj-jav)
                            err1=dmax1(err1,err)
                        End If
                    End Do
                    If (err1 < trsd) Then
                        ierr=0
                        Exit
                    End If
                    If (it /= 0) Then
                        Do i=1,num
                            proceed=(lin-1)*Iconverge(i)==0
                            If (proceed) Then
                                s=0_type_real
                                Do n=1,Nd
                                    a=col_in(n,i)-f*col_out(n,i)
                                    s=s+a*a
                                    col_in(n,i)=a
                                End Do
                                s1=1.d0/dsqrt(s)
                                col_in(:,i)=col_in(:,i)*s1
                                ArrB(1:nd_local, lin+i-1) = col_in(nd_start+1:nd_end, i)
                            End If
                        End Do
                    End If

                End If
            End Do

            Deallocate(col_in, col_out)
        End Block
        If (mype == 0) Then
            strfmt = '(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)'
            Write( *,strfmt) Njd+1-it,err1,ierr
            Write(11,strfmt) Njd+1-it,err1,ierr
            If (ierr /= 0) Then
                Write(*,*) ' Wrong J values for Probe vectors '
            End If
        End If
        Return
    End Subroutine Prj_J

End Module davidson
