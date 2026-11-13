Module hblocks
    ! Block partitioning of the Hamiltonian for the Davidson product H*v.
    !
    ! The vector rows are split into q blocks, about sqrt(npes/4) (H2DBlocks).
    ! Block a is made of the ArrB row slices of a contiguous group of ranks (group a).
    ! Each rank in group a also gets a second block index b, so it is labelled (a,b).
    ! H is stored as one triangle; each stored element (n,k) belongs to the block pair
    ! {block of n, block of k} and is given to one of the ranks labelled (a,b) or (b,a) for that pair.
    ! A rank therefore needs only the input blocks a and b, and adds output only into blocks a and b.
    ! The second indices are chosen from the number of elements in each block pair,
    ! so that pairs with more elements get more ranks (AssignLabels).
    ! Within a pair, each rank gets a contiguous range of whole rows with about equal element counts.
    ! Whole rows keep the runs of neighbouring columns together, as in the 1D layout;
    ! spreading a row's elements over the pair's ranks made the local product ~2x slower on Os16+.
    !
    ! One product y = H*x takes four collectives of about Nd/q values each:
    !   1. Allgatherv over group a gives x on block a.
    !   2. Bcast over the ranks with second index b gives x on block b.
    !   3. Each rank multiplies its own elements.
    !   4. Reduce over the same ranks sums the block b output on their root.
    !   5. Reduce_scatter over group a sums the block a output and returns each rank its own rows.
    ! The 1D layout instead holds and exchanges full-length vectors on every rank.
    !
    ! The root of the ranks with second index b is the rank of group b labelled (b,b).
    ! Its block a is block b, so it already has x there, and it adds the summed block b output
    ! into its own output before step 5.
    ! Other ranks labelled (b,b) skip steps 2 and 4.
    Use mpi_f08
    Use conf_variables

    Implicit None

    Private

    ! Target number of ranks per block pair, which sets the number of blocks (H2DBlocks)
    Integer, Parameter :: ranks_per_pair = 4

    ! True once BuildH2D has run; Mxmpy then uses MatvecH2D
    Logical, Save :: h2d_active = .False.

    ! Number of blocks, this rank's block indices, and the start row of each block (0-based)
    Integer, Save :: q = 0, my_a = 0, my_b = 0
    Integer, Allocatable, Save :: blk_start(:)
    ! First rank of each group
    Integer, Allocatable, Save :: g_start(:)
    ! Ranks of each block pair (lo,hi), lo <= hi, at pr_list(pr_start(lo*q+hi)+1 : ... + pr_cnt(lo*q+hi))
    Integer, Allocatable, Save :: pr_start(:), pr_cnt(:), pr_list(:)
    ! Rows of block a, rows of block b (0 when b = a), and La + Lb
    Integer, Save :: La = 0, Lb = 0, Lc = 0

    ! Group a, and the ranks that share second index b
    Type(MPI_Comm), Save :: row_comm, col_comm
    Logical, Save :: in_col = .False., is_col_root = .False.
    Integer, Save :: col_size = 1
    ! Row slice sizes and offsets of the members of group a, for Allgatherv and Reduce_scatter
    Integer, Allocatable, Save :: grp_counts(:), grp_displs(:)

    ! Local elements in CSR form; h2_rows lists the nonempty rows (combined index:
    ! block a rows, then block b rows), and row ir has elements h2_rowptr(ir-1)+1 ... h2_rowptr(ir)
    Integer, Save :: nr = 0
    Integer, Allocatable, Save :: h2_rows(:)
    Integer(Kind=Int64), Allocatable, Save :: h2_rowptr(:)
    Integer, Allocatable, Save :: h2_col(:)
    Real(type_real), Allocatable, Save :: h2_val(:)

    ! Input and output on the combined index, and the summed block b output on the root
    Real(type_real), Allocatable, Save :: xc(:), yc(:), ycol(:)

    ! Optional timing of each phase of H*v, for both layouts (environment variable PCI_MATVEC_TIMING).
    ! With it on, barriers before and after the local product separate waiting from communication,
    ! and a summary goes to CONF.PRG at the end of Davidson.
    ! Phases: 1 wait before the product, 2 Allgatherv, 3 Bcast, 4 local product,
    ! 5 wait after the local product (load imbalance), 6 Reduce, 7 Reduce_scatter.
    Integer, Parameter :: n_mv_phases = 7
    Logical, Save :: mv_timing = .False.
    Real(dp), Save :: mv_t(n_mv_phases) = 0.d0
    Integer(Kind=Int64), Save :: mv_count = 0

    Public :: h2d_active, BuildH2D, MatvecH2D, FreeH2D, H2DBlocks
    Public :: mv_timing, MvTimingInit, MvMark, MvWait, MvCount, MvTimingReport

  Contains

    Subroutine MvTimingInit
        ! Read PCI_MATVEC_TIMING and reset the phase timers
        Use env_var, Only : GetEnvLogical
        Implicit None
        Integer :: mpierr
        ! Rank 0 decides, so that every rank takes the same timed barriers
        mv_timing = GetEnvLogical('PCI_MATVEC_TIMING', .False.)
        Call MPI_Bcast(mv_timing, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        mv_t = 0.d0
        mv_count = 0
    End Subroutine MvTimingInit

    Subroutine MvMark(phase, t0)
        ! Add the time since t0 to a phase and restart the clock
        Implicit None
        Integer, Intent(In) :: phase
        Real(dp), Intent(InOut) :: t0
        Real(dp) :: t
        If (.not. mv_timing) Return
        t = MPI_Wtime()
        mv_t(phase) = mv_t(phase) + (t - t0)
        t0 = t
    End Subroutine MvMark

    Subroutine MvWait(phase, t0)
        ! Barrier over all ranks, timed as a waiting phase (only when timing is on)
        Implicit None
        Integer, Intent(In) :: phase
        Real(dp), Intent(InOut) :: t0
        Integer :: mpierr
        If (.not. mv_timing) Return
        t0 = MPI_Wtime()
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MvMark(phase, t0)
    End Subroutine MvWait

    Subroutine MvCount
        ! Count one product
        Implicit None
        mv_count = mv_count + 1
    End Subroutine MvCount

    Subroutine MvTimingReport(mype, npes)
        ! Write the largest and average time per phase over all ranks to CONF.PRG (unit 77, rank 0)
        Implicit None
        Integer, Intent(In) :: mype, npes
        Real(dp) :: tmax(n_mv_phases), tsum(n_mv_phases)
        Integer :: i, mpierr
        Character(Len=24), Parameter :: names(n_mv_phases) = [Character(Len=24) :: &
            'wait before product', 'Allgatherv', 'Bcast (2D)', 'local product', &
            'wait after product', 'Reduce (2D)', 'Reduce_scatter']

        If (.not. mv_timing) Return
        Call MPI_Reduce(mv_t, tmax, n_mv_phases, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Reduce(mv_t, tsum, n_mv_phases, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
        If (mype /= 0) Return
        If (h2d_active) Then
            Write(77,'(A,I0,A,I0,A,I0,A)') 'H*v timing: ', mv_count, ' products, 2D layout, ', q, ' x ', q, ' blocks'
        Else
            Write(77,'(A,I0,A)') 'H*v timing: ', mv_count, ' products, 1D layout'
        End If
        Write(77,'(2X,A24,2A12)') 'phase', 'max (s)', 'avg (s)'
        Do i = 1, n_mv_phases
            Write(77,'(2X,A24,2F12.2)') names(i), tmax(i), tsum(i)/npes
        End Do
    End Subroutine MvTimingReport

    Integer Function H2DBlocks(npes)
        ! Number of vector blocks for npes ranks: the largest q with ranks_per_pair*q*q <= npes.
        ! About 4 ranks per block pair leave room to give heavy pairs more ranks;
        ! with q = sqrt(npes) each group has only about one rank per pair, and the load balance
        ! measured on the small Fe16+ test was ~60% instead of ~90%.
        Integer, Intent(In) :: npes
        H2DBlocks = 1
        Do While (ranks_per_pair*(H2DBlocks+1)*(H2DBlocks+1) <= npes)
            H2DBlocks = H2DBlocks + 1
        End Do
    End Function H2DBlocks

    Integer Function BlockOf(n)
        ! Block that holds global row n (1-based), by bisection on blk_start
        Integer, Intent(In) :: n
        Integer :: lo, hi, mid
        lo = 0
        hi = q - 1
        Do While (lo < hi)
            mid = (lo + hi + 1) / 2
            If (blk_start(mid) < n) Then
                lo = mid
            Else
                hi = mid - 1
            End If
        End Do
        BlockOf = lo
    End Function BlockOf

    Integer Function OtherGroup(a, b, g)
        ! The group of the pair (a,b) that is not g
        Integer, Intent(In) :: a, b, g
        OtherGroup = a + b - g
    End Function OtherGroup

    Subroutine AssignLabels(own, label)
        ! Choose the second block index of every rank from the element counts of the block pairs.
        ! own(A,B) is the number of stored elements with row in block A and column in block B,
        ! so the rows of group A already hold them in the 1D layout.
        ! Each group first gives one rank to its own diagonal pair.
        ! Each off-diagonal pair then gets one rank from the group that holds more of its elements;
        ! a group is in q-1 such pairs and has at least q-1 ranks left, so this always succeeds.
        ! The remaining ranks go, one at a time, to the pair with the most elements per rank,
        ! again preferring the group that holds more of its elements.
        ! The result depends only on own, so every rank computes the same labels.
        Integer(Kind=Int64), Intent(In) :: own(0:q-1,0:q-1)
        Integer, Intent(Out) :: label(0:)
        Integer :: a, b, g, best_a, best_b, spare(0:q-1), c(0:q-1,0:q-1), r, nleft
        Real(dp) :: w, best_w

        c = 0
        Do a = 0, q-1
            c(a,a) = 1
            spare(a) = g_start(a+1) - g_start(a) - 1
        End Do
        Do a = 0, q-1
            Do b = a+1, q-1
                g = a
                If (own(b,a) > own(a,b)) g = b
                If (own(b,a) == own(a,b) .and. spare(b) > spare(a)) g = b
                If (spare(g) == 0) g = OtherGroup(a, b, g)
                c(g,OtherGroup(a, b, g)) = c(g,OtherGroup(a, b, g)) + 1
                spare(g) = spare(g) - 1
            End Do
        End Do

        nleft = Sum(spare)
        Do While (nleft > 0)
            best_w = -1.d0
            best_a = 0
            best_b = 0
            Do a = 0, q-1
                Do b = a, q-1
                    If (spare(a) == 0 .and. spare(b) == 0) Cycle
                    If (a == b) Then
                        w = Real(own(a,a), dp) / c(a,a)
                    Else
                        w = Real(own(a,b) + own(b,a), dp) / (c(a,b) + c(b,a))
                    End If
                    If (w > best_w) Then
                        best_w = w
                        best_a = a
                        best_b = b
                    End If
                End Do
            End Do
            g = best_a
            If (own(best_b,best_a) > own(best_a,best_b)) g = best_b
            If (spare(g) == 0) g = OtherGroup(best_a, best_b, g)
            c(g,OtherGroup(best_a, best_b, g)) = c(g,OtherGroup(best_a, best_b, g)) + 1
            spare(g) = spare(g) - 1
            nleft = nleft - 1
        End Do

        ! Ranks of group a take second index 0 first, then 1, and so on
        Do a = 0, q-1
            r = g_start(a)
            Do b = 0, q-1
                label(r:r+c(a,b)-1) = b
                r = r + c(a,b)
            End Do
        End Do
    End Subroutine AssignLabels

    Integer Function BlockOfRank(r)
        ! Group (block) of rank r
        Integer, Intent(In) :: r
        Integer :: a
        BlockOfRank = 0
        Do a = 1, q-1
            If (r >= g_start(a)) BlockOfRank = a
        End Do
    End Function BlockOfRank

    Integer Function CombinedIndex(n)
        ! Position of global row n in the combined index (block a rows, then block b rows)
        Integer, Intent(In) :: n
        If (n > blk_start(my_a) .and. n <= blk_start(my_a+1)) Then
            CombinedIndex = n - blk_start(my_a)
        Else
            CombinedIndex = La + n - blk_start(my_b)
        End If
    End Function CombinedIndex

    Subroutine BuildH2D(npes, mype, mpi_type_real)
        ! Move the row-distributed CSR H (after RedistributeHamCSR) into the 2D layout.
        ! The diagonal shift Hamil%minval is applied here, as Mxmpy does in the 1D product.
        ! Hamil%col, Hamil%val and HamilCSR_rowptr are freed.
        ! Needs nd_counts_all and nd_displs_all from AllocateDvdsnArrays.
        Use str_fmt, Only : startTimer, stopTimer
        Implicit None
        Integer, Intent(In) :: npes, mype
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: a, r, base, extra, color, key, mpierr, lr, n, k, i, ir, nd_loc, grank
        Integer, Allocatable :: label(:), dest(:)
        Integer(Kind=Int64), Allocatable :: own(:,:), pre(:), run(:)
        Integer(Kind=Int64) :: l8, ipos, nnz_old, nnz_new, nnz_total, nnz_max
        Integer, Allocatable :: send_cnts(:), recv_cnts(:), send_disp(:), recv_disp(:)
        Integer, Allocatable :: send_int(:), recv_n(:), recv_k(:)
        Real(type_real), Allocatable :: send_val(:), recv_val(:)
        Integer(Kind=Int64), Allocatable :: tmp_pos(:)
        Integer(Kind=Int64) :: s1
        Character(Len=16) :: timeStr

        If (mype == 0) Call startTimer(s1)
        nd_loc = nd_end - nd_start
        ! Local element count of the 1D H (ih8 is not reliable after reading pCONF.JJJ)
        nnz_old = HamilCSR_rowptr(nd_loc)

        ! Groups: q contiguous groups of npes/q ranks or one more
        q = H2DBlocks(npes)
        base = npes / q
        extra = Mod(npes, q)
        Allocate(g_start(0:q), blk_start(0:q))
        Do a = 0, q
            g_start(a) = a*base + Min(a, extra)
        End Do
        Do a = 0, q-1
            blk_start(a) = nd_displs_all(g_start(a))
        End Do
        blk_start(q) = Nd

        my_a = 0
        Do a = 0, q-1
            If (mype >= g_start(a)) my_a = a
        End Do

        ! Elements by row block and column block, summed over all ranks
        Allocate(own(0:q-1,0:q-1))
        own = 0
        Do lr = 1, nd_loc
            n = nd_start + lr
            a = BlockOf(n)
            Do l8 = HamilCSR_rowptr(lr-1)+1_int64, HamilCSR_rowptr(lr)
                i = BlockOf(Hamil%col(l8))
                own(a,i) = own(a,i) + 1
            End Do
        End Do
        Call MPI_Allreduce(MPI_IN_PLACE, own, q*q, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)

        ! Second index of every rank, and the ranks of each block pair
        Allocate(label(0:npes-1), pr_start(0:q*q-1), pr_cnt(0:q*q-1), pr_list(npes))
        Call AssignLabels(own, label)
        pr_cnt = 0
        Do r = 0, npes-1
            a = BlockOfRank(r)
            i = Min(a, label(r))*q + Max(a, label(r))
            pr_cnt(i) = pr_cnt(i) + 1
        End Do
        pr_start(0) = 0
        Do i = 1, q*q-1
            pr_start(i) = pr_start(i-1) + pr_cnt(i-1)
        End Do
        pr_cnt = 0
        Do r = 0, npes-1
            a = BlockOfRank(r)
            i = Min(a, label(r))*q + Max(a, label(r))
            pr_cnt(i) = pr_cnt(i) + 1
            pr_list(pr_start(i) + pr_cnt(i)) = r
        End Do

        my_b = label(mype)
        La = blk_start(my_a+1) - blk_start(my_a)
        Lb = 0
        If (my_b /= my_a) Lb = blk_start(my_b+1) - blk_start(my_b)
        Lc = La + Lb

        ! Group communicator, with the members' row slices for Allgatherv and Reduce_scatter
        Call MPI_Comm_split(MPI_COMM_WORLD, my_a, mype, row_comm, mpierr)
        Allocate(grp_counts(g_start(my_a+1)-g_start(my_a)), grp_displs(g_start(my_a+1)-g_start(my_a)))
        Do r = g_start(my_a), g_start(my_a+1)-1
            grp_counts(r-g_start(my_a)+1) = nd_counts_all(r)
            grp_displs(r-g_start(my_a)+1) = nd_displs_all(r) - blk_start(my_a)
        End Do

        ! Ranks sharing second index b; the root is the first rank labelled (b,b) and gets key 0
        is_col_root = .False.
        If (my_a == my_b) Then
            r = g_start(my_a)
            Do While (label(r) /= my_a)
                r = r + 1
            End Do
            is_col_root = (r == mype)
        End If
        Deallocate(label)
        in_col = (my_a /= my_b) .or. is_col_root
        color = MPI_UNDEFINED
        If (in_col) color = my_b
        key = mype + 1
        If (is_col_root) key = 0
        Call MPI_Comm_split(MPI_COMM_WORLD, color, key, col_comm, mpierr)
        col_size = 1
        If (in_col) Call MPI_Comm_size(col_comm, col_size, mpierr)

        ! Each pair's elements are numbered in row order: rows of its lower block first, then of its
        ! higher block, and each row in order of the 1D layout.
        ! pre(B) counts the elements with column block B in the rows of earlier ranks of this group,
        ! and run(B) those in this rank's rows already passed; together they place a row in its pair.
        Allocate(pre(0:q-1), run(0:q-1), dest(0:q-1))
        run = 0
        Do l8 = 1, nnz_old
            i = BlockOf(Hamil%col(l8))
            run(i) = run(i) + 1
        End Do
        Call MPI_Exscan(run, pre, q, MPI_INTEGER8, MPI_SUM, row_comm, mpierr)
        Call MPI_Comm_rank(row_comm, grank, mpierr)
        If (grank == 0) pre = 0

        ! Count elements per destination rank
        Allocate(send_cnts(0:npes-1), recv_cnts(0:npes-1), send_disp(0:npes-1), recv_disp(0:npes-1))
        send_cnts = 0
        run = 0
        Do lr = 1, nd_loc
            Call RowDests
            Do l8 = HamilCSR_rowptr(lr-1)+1_int64, HamilCSR_rowptr(lr)
                i = BlockOf(Hamil%col(l8))
                send_cnts(dest(i)) = send_cnts(dest(i)) + 1
                run(i) = run(i) + 1
            End Do
        End Do
        Call MPI_Alltoall(send_cnts, 1, MPI_INTEGER, recv_cnts, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        send_disp(0) = 0
        recv_disp(0) = 0
        Do r = 1, npes-1
            send_disp(r) = send_disp(r-1) + send_cnts(r-1)
            recv_disp(r) = recv_disp(r-1) + recv_cnts(r-1)
        End Do
        nnz_new = Int(recv_disp(npes-1), Int64) + recv_cnts(npes-1)
        Allocate(tmp_pos(0:npes-1))

        ! Pass 1: values, with the diagonal shift applied; Hamil%val is freed after packing
        Allocate(send_val(Max(nnz_old,1_int64)))
        tmp_pos = send_disp + 1
        run = 0
        Do lr = 1, nd_loc
            n = nd_start + lr
            Call RowDests
            Do l8 = HamilCSR_rowptr(lr-1)+1_int64, HamilCSR_rowptr(lr)
                k = Hamil%col(l8)
                i = BlockOf(k)
                r = dest(i)
                send_val(tmp_pos(r)) = Hamil%val(l8)
                If (n == k) send_val(tmp_pos(r)) = Hamil%val(l8) - Hamil%minval
                tmp_pos(r) = tmp_pos(r) + 1
                run(i) = run(i) + 1
            End Do
        End Do
        Deallocate(Hamil%val)
        Allocate(recv_val(Max(nnz_new,1_int64)))
        Call MPI_Alltoallv(send_val, send_cnts, send_disp, mpi_type_real, &
                           recv_val, recv_cnts, recv_disp, mpi_type_real, MPI_COMM_WORLD, mpierr)
        Deallocate(send_val)

        ! Pass 2: rows
        Allocate(send_int(Max(nnz_old,1_int64)))
        tmp_pos = send_disp + 1
        run = 0
        Do lr = 1, nd_loc
            n = nd_start + lr
            Call RowDests
            Do l8 = HamilCSR_rowptr(lr-1)+1_int64, HamilCSR_rowptr(lr)
                i = BlockOf(Hamil%col(l8))
                r = dest(i)
                send_int(tmp_pos(r)) = n
                tmp_pos(r) = tmp_pos(r) + 1
                run(i) = run(i) + 1
            End Do
        End Do
        Allocate(recv_n(Max(nnz_new,1_int64)))
        Call MPI_Alltoallv(send_int, send_cnts, send_disp, MPI_INTEGER, &
                           recv_n, recv_cnts, recv_disp, MPI_INTEGER, MPI_COMM_WORLD, mpierr)

        ! Pass 3: columns; the 1D H is freed after packing
        tmp_pos = send_disp + 1
        run = 0
        Do lr = 1, nd_loc
            Call RowDests
            Do l8 = HamilCSR_rowptr(lr-1)+1_int64, HamilCSR_rowptr(lr)
                i = BlockOf(Hamil%col(l8))
                r = dest(i)
                send_int(tmp_pos(r)) = Hamil%col(l8)
                tmp_pos(r) = tmp_pos(r) + 1
                run(i) = run(i) + 1
            End Do
        End Do
        Deallocate(Hamil%col, HamilCSR_rowptr)
        Allocate(recv_k(Max(nnz_new,1_int64)))
        Call MPI_Alltoallv(send_int, send_cnts, send_disp, MPI_INTEGER, &
                           recv_k, recv_cnts, recv_disp, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        Deallocate(send_int, tmp_pos, send_cnts, recv_cnts, send_disp, recv_disp)
        Deallocate(own, pre, run, dest)

        ! Counting sort into CSR over the nonempty rows of the combined index
        Allocate(tmp_pos(Max(Lc,1)))
        tmp_pos = 0
        Do l8 = 1, nnz_new
            recv_n(l8) = CombinedIndex(recv_n(l8))
            recv_k(l8) = CombinedIndex(recv_k(l8))
            tmp_pos(recv_n(l8)) = tmp_pos(recv_n(l8)) + 1
        End Do
        nr = Count(tmp_pos(1:Lc) > 0)
        Allocate(h2_rows(Max(nr,1)), h2_rowptr(0:nr))
        h2_rowptr(0) = 0
        ir = 0
        Do i = 1, Lc
            If (tmp_pos(i) > 0) Then
                ir = ir + 1
                h2_rows(ir) = i
                h2_rowptr(ir) = h2_rowptr(ir-1) + tmp_pos(i)
                ! From here tmp_pos(i) is the position of row i in h2_rows
                tmp_pos(i) = ir
            End If
        End Do
        Do l8 = 1, nnz_new
            recv_n(l8) = Int(tmp_pos(recv_n(l8)))
        End Do
        Deallocate(tmp_pos)
        ! Columns, then values, so that each input array is freed before the next output is made
        ! (peak 24 bytes per element, as in RedistributeHamCSR)
        Allocate(tmp_pos(Max(nr,1)))
        tmp_pos(1:nr) = h2_rowptr(0:nr-1) + 1
        Allocate(h2_col(Max(nnz_new,1_int64)))
        Do l8 = 1, nnz_new
            ir = recv_n(l8)
            ipos = tmp_pos(ir)
            h2_col(ipos) = recv_k(l8)
            tmp_pos(ir) = ipos + 1
        End Do
        Deallocate(recv_k)
        tmp_pos(1:nr) = h2_rowptr(0:nr-1) + 1
        Allocate(h2_val(Max(nnz_new,1_int64)))
        Do l8 = 1, nnz_new
            ir = recv_n(l8)
            ipos = tmp_pos(ir)
            h2_val(ipos) = recv_val(l8)
            tmp_pos(ir) = ipos + 1
        End Do
        Deallocate(recv_n, recv_val, tmp_pos)
        ih8 = nnz_new

        Allocate(xc(Max(Lc,1)), yc(Max(Lc,1)))
        If (is_col_root) Allocate(ycol(Max(La,1)))
        h2d_active = .True.

        ! Load balance: average over largest local element count
        nnz_total = nnz_new
        nnz_max = nnz_new
        Call MPI_Allreduce(MPI_IN_PLACE, nnz_total, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_Allreduce(MPI_IN_PLACE, nnz_max, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
        If (mype == 0) Then
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A,F5.1,A,I0,A,I0,A)') 'BuildH2D: H*v load balance: ', &
                100.0_dp * Real(nnz_total, dp) / Real(Int(npes, Int64) * Max(nnz_max,1_int64), dp), &
                '% efficiency (avg/max element count per rank), ', q, ' x ', q, ' blocks'
            Write(*,'(2X,A)') 'TIMING >>> BuildH2D done in ' // Trim(timeStr)
        End If

      Contains

        Subroutine RowDests
            ! Destination rank, for each column block B, of the elements of the current row.
            ! The row's position in its pair, out of the pair's element count, picks one of the
            ! pair's ranks, so each rank gets a contiguous range of rows with about equal counts.
            Integer :: B, lo, hi, pid, m
            Integer(Kind=Int64) :: tot, pos
            Do B = 0, q-1
                lo = Min(my_a, B)
                hi = Max(my_a, B)
                pid = lo*q + hi
                m = pr_cnt(pid)
                If (lo == hi) Then
                    tot = own(lo,lo)
                    pos = pre(B) + run(B)
                Else
                    tot = own(lo,hi) + own(hi,lo)
                    pos = pre(B) + run(B)
                    If (my_a == hi) pos = pos + own(lo,hi)
                End If
                If (tot > 0) Then
                    dest(B) = pr_list(pr_start(pid) + Int(Min(Int(m-1,Int64), pos*m/tot)) + 1)
                Else
                    dest(B) = pr_list(pr_start(pid) + 1)
                End If
            End Do
        End Subroutine RowDests
    End Subroutine BuildH2D

    Subroutine MatvecH2D(xloc, yloc, mpi_type_real)
        ! yloc = H * x, where xloc and yloc are this rank's row slices (nd_local rows)
        Implicit None
        Real(type_real), Intent(In) :: xloc(:)
        Real(type_real), Intent(InOut) :: yloc(:)
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: i, j, ir, mpierr
        Integer(Kind=Int64) :: l8
        Real(type_real) :: t
        ! Receive buffer for Reduce, not used on non-root ranks
        Real(type_real) :: unused(1)
        Real(dp) :: t0

        Call MvCount
        Call MvWait(1, t0)

        ! x on block a, from the group
        Call MPI_Allgatherv(xloc, nd_local, mpi_type_real, xc, grp_counts, grp_displs, &
                            mpi_type_real, row_comm, mpierr)
        Call MvMark(2, t0)
        ! x on block b, from the root (b,b), whose block a is block b
        If (in_col .and. col_size > 1) Then
            If (is_col_root) Then
                Call MPI_Bcast(xc, La, mpi_type_real, 0, col_comm, mpierr)
            Else
                Call MPI_Bcast(xc(La+1:Lc), Lb, mpi_type_real, 0, col_comm, mpierr)
            End If
        End If
        Call MvMark(3, t0)

        ! Local product; each stored off-diagonal element also acts as its transpose
        yc(1:Lc) = 0_type_real
        Do ir = 1, nr
            i = h2_rows(ir)
            Do l8 = h2_rowptr(ir-1)+1_int64, h2_rowptr(ir)
                j = h2_col(l8)
                t = h2_val(l8)
                yc(i) = yc(i) + t * xc(j)
                If (i /= j) yc(j) = yc(j) + t * xc(i)
            End Do
        End Do
        Call MvMark(4, t0)
        Call MvWait(5, t0)

        ! Sum the block b output on the root, which adds it to its own block
        If (in_col .and. col_size > 1) Then
            If (is_col_root) Then
                ycol(1:La) = 0_type_real
                Call MPI_Reduce(MPI_IN_PLACE, ycol, La, mpi_type_real, MPI_SUM, 0, col_comm, mpierr)
                yc(1:La) = yc(1:La) + ycol(1:La)
            Else
                Call MPI_Reduce(yc(La+1:Lc), unused, Lb, mpi_type_real, MPI_SUM, 0, col_comm, mpierr)
            End If
        End If
        Call MvMark(6, t0)

        ! Sum the block a output over the group; each rank receives its own rows
        Call MPI_Reduce_scatter(yc, yloc, grp_counts, mpi_type_real, MPI_SUM, row_comm, mpierr)
        Call MvMark(7, t0)
    End Subroutine MatvecH2D

    Subroutine FreeH2D
        ! Release the 2D H, buffers and communicators
        Implicit None
        Integer :: mpierr
        If (.not. h2d_active) Return
        Deallocate(h2_rows, h2_rowptr, h2_col, h2_val, xc, yc, grp_counts, grp_displs, g_start, blk_start)
        Deallocate(pr_start, pr_cnt, pr_list)
        If (Allocated(ycol)) Deallocate(ycol)
        Call MPI_Comm_free(row_comm, mpierr)
        If (in_col) Call MPI_Comm_free(col_comm, mpierr)
        h2d_active = .False.
    End Subroutine FreeH2D

End Module hblocks
