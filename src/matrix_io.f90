Module matrix_io
    ! 
    ! This module implements MPI functions that read and write matrix elements to disc. 
    !
    Use conf_variables
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime

    Implicit None

    Private

    Public :: WriteMatrix, ReadMatrix, RedistributeHamCSR

  Contains

    Subroutine MPIErrHandle(mpierr)
        Use mpi_f08
        Implicit None
        Integer, Intent(In) :: mpierr

        Select Case(mpierr)
        Case(MPI_ERR_BUFFER)
            Print*, 'Invalid buffer pointer'
            Stop
        Case(MPI_ERR_COUNT)
            Print*, 'Invalid count argument'
            Stop
        Case(MPI_ERR_TYPE)
            Print*, 'Invalid datatype argument'
            Stop
        Case(MPI_ERR_FILE)
            Print*, 'Invalid file handle'
            Stop
        Case(MPI_ERR_NO_SUCH_FILE)
            Print*, 'File does not exist'
            Stop
        Case(MPI_ERR_FILE_EXISTS)
            Print*, 'File exists'
            Stop
        Case(MPI_ERR_ACCESS)
            Print*, 'Permission denied'
            Stop
        Case(MPI_ERR_NO_SPACE)
            Print*, 'Not enough space on disk'
            Stop
        Case(MPI_ERR_QUOTA)
            Print*, 'Quota exceeded'
            Stop
        Case(MPI_SUCCESS)
            Continue
        End Select
    End Subroutine MPIErrHandle

    Subroutine WriteMatrix(mat,num_elements_per_core,num_elements_total,filename,mype,npes,mpierr)
        !
        ! This subroutine writes a parallel file
        ! Matrix 'mat' is written to the file 'filename'
        ! The matrix should be of Type(Matrix) and have indices n,k and value t
        ! The structure of the file is as follows:
        ! =============================================================
        ! | number of processors     |    counters for each processor | 
        ! | mat%n (core 0) | mat%n (core 1) | ... | mat%n (core npes) |
        ! | mat%k (core 0) | mat%k (core 1) | ... | mat%k (core npes) |
        ! | mat%t (core 0) | mat%t (core 1) | ... | mat%t (core npes) |
        ! =============================================================
        !
        Use mpi_f08
        Implicit None
        Type(Matrix(type_real)) :: mat
        Integer, Intent(In)             :: num_elements_per_core
        Integer(Kind=int64), Intent(In) :: num_elements_total

        Integer                       :: mype, npes, i, mpierr
        Type(MPI_File)                :: fh
        Integer(Kind=Int64)           :: start_time
        Integer, Dimension(npes)      :: sizes
        Integer(kind=MPI_OFFSET_KIND), Dimension(npes) :: disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=*)              :: filename
        Character(Len=16)             :: timeStr
    
        If (mype == 0) Then    
            If (filename == 'CONFp.HIJ') Then     
                Open(66,file='progress.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
                Write(66) Nc, Nd
                Close(66)   
            End If
            Open(66,file='nprocs.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Write(66) npes
            Close(66)
            print*, 'Writing ' // filename // '...'
            Call startTimer(start_time)
        End If

        sizes=0
        Call MPI_AllGather(num_elements_per_core, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        disps=0_MPI_OFFSET_KIND
        Do i=2,npes
            disps(i)=disps(i-1)+Int(sizes(i-1), kind=MPI_OFFSET_KIND)
        End Do
            
        ! Open file
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write number of matrix elements
        disp = mype * 4_MPI_OFFSET_KIND
        Call MPI_FILE_WRITE_AT_ALL(fh, disp, num_elements_per_core, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%row
        disp = (npes + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_WRITE_AT(fh, disp, mat%row, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%col
        disp = (npes + num_elements_total + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_WRITE_AT(fh, disp, mat%col, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%val
        Select Case(type_real)
        Case(sp)
            disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 4_MPI_OFFSET_KIND
            Call MPI_FILE_WRITE_AT(fh, disp, mat%val, num_elements_per_core, MPI_REAL, MPI_STATUS_IGNORE, mpierr) 
            Call MPIErrHandle(mpierr)
        Case(dp)
            disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 8_MPI_OFFSET_KIND
            Call MPI_FILE_WRITE_AT(fh, disp, mat%val, num_elements_per_core, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
            Call MPIErrHandle(mpierr)
        Case default 
            print*,'unrecognized real kind'
            Stop
        End Select

        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Print time to write file
        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            write(*,'(2X,A)') 'TIMING >>> Writing ' // filename // ' took ' // trim(timeStr)
        End If
    End Subroutine WriteMatrix

    Subroutine ReadMatrix(indices1,indices2,values,num_elements_per_core,num_elements_total,filename,mype,npes,mpierr)
        !
        ! This subroutine reads a parallel matrix file 
        !
        Use mpi_f08
        Implicit None
        Integer, Allocatable, Dimension(:)  :: indices1, indices2
        Real(kind=type_real), Allocatable, Dimension(:) :: values
        Integer, Intent(Out)             :: num_elements_per_core
        Integer(Kind=int64), Intent(Out) :: num_elements_total

        Integer                                         :: mype, npes, npes_read, i, mpierr
        Type(MPI_File)                                  :: fh
        Integer(kind=int64)                             :: start_time
        Integer, dimension(npes)                        :: sizes
        Integer(kind=MPI_OFFSET_KIND), Dimension(npes)  :: disps
        Integer(kind=MPI_OFFSET_KIND)                   :: disp=0_MPI_OFFSET_KIND   
        Character(Len=*)                                :: filename
        Character(Len=16)                               :: timeStr

        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, mpierr) 
        Call MPIErrHandle(mpierr)

        If (mype == 0) then
            ! Stop program if file does not exist
            If (mpierr /= 0) Then
                print*, 'File ',filename,' does not exist'
                Stop
            End If

            ! Read number of processors
            Open(66,file='nprocs.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Read(66) npes_read
            Close(66)
            If (npes /= npes_read) then
                print*, 'Number of processors inconsistent. '// filename //' was written with ', npes_read, ' processors,', &
                        ' but', npes, ' is available.'
                stop
            End If
            print*, 'Reading ' // filename // '...'
            Call startTimer(start_time)
        End If

        Call MPI_Bcast(Nd_prev, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc_prev, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        ! Read counters
        disp = mype * 4_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, num_elements_per_core, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr)
        Call MPIErrHandle(mpierr)

        sizes=0
        Call MPI_AllGather(num_elements_per_core, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)

        Do i=1,npes
            num_elements_total = num_elements_total + sizes(i)
        End Do

        ! Distribute total number of matrix elements equally across all cores
        num_elements_per_core = num_elements_total/npes
        if (mype == 0) num_elements_per_core = num_elements_per_core + mod(num_elements_total,npes)
        Call MPI_AllGather(num_elements_per_core, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)

        ! Calculate displacements
        disps=0_MPI_OFFSET_KIND
        Do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        End Do

        ! Allocate matrix if not allocated yet
        If (.not. allocated(indices1)) allocate(indices1(num_elements_per_core))
        If (.not. allocated(indices2)) allocate(indices2(num_elements_per_core))
        If (.not. allocated(values)) allocate(values(num_elements_per_core))

        ! Read first indices of matrix elements mat%n
        disp = (npes + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, indices1, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Read second indices of matrix elements mat%k
        disp = (npes + num_elements_total + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, indices2, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Read values of matrix element mat%t
        Select Case(type_real)
        Case(sp)
            disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 4_MPI_OFFSET_KIND
            Call MPI_FILE_READ_AT(fh, disp, values, num_elements_per_core, MPI_REAL, MPI_STATUS_IGNORE, mpierr) 
            Call MPIErrHandle(mpierr)
        Case(dp)
            disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 8_MPI_OFFSET_KIND
            Call MPI_FILE_READ_AT(fh, disp, values, num_elements_per_core, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
            Call MPIErrHandle(mpierr)
        Case default 
            print*,'unrecognized real kind'
            Stop
        End Select

        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Print time to write file
        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            write(*,'(2X,A)') 'TIMING >>> Reading ' // filename // ' took ' // trim(timeStr)
        End If
    
        Return
    End Subroutine ReadMatrix

    Subroutine RedistributeHamCSR(npes, mype, mpi_type_real)
        ! Redistribute the Hamiltonian from COO to CSR
        !
        ! Steps performed:
        !   1. MPI_Allreduce(SUM) to get global per-row element counts (row_cnts_g).
        !   2. Greedy scan to assign contiguous row ranges to ranks (rank_of_row);
        !      derive nd_start/nd_end/nd_loc for each rank.
        !   3. Scan local elements to build send counts; MPI_Alltoall transposes
        !      the send-count matrix to give each rank its receive counts.
        !      Compute send/recv displacements and new_ih8 (post-redistribution count).
        !   4-6. Three sequential pack+Alltoallv passes (col, val, row). Hamil%row serves
        !      as the rank-lookup key for all three passes so it is packed last.
        !      Peak: Hamil%row+val(12B×ih8) + col_new(4B×new) + send_val(8B×ih8) = 24B/element.
        !   7. Counting sort into CSR: scatter col then val into row order; row freed after val pass.
        !   8. Print load balance: ih8_max via MPI_MAX; rank 0 reports efficiency = avg/max.
        !
        ! After this call:
        !   - Each rank owns all nonzeros for rows [nd_start+1 : nd_end]
        !   - Hamil%col/val sorted in CSR row order; Hamil%row is freed
        !   - HamilCSR_rowptr(0:nd_loc): prefix-sum of per-local-row element counts
        !   - ih8 = new_ih8 (post-redistribution local element count)
        Use mpi_f08
        Implicit None

        Integer, Intent(In) :: npes, mype
        Type(MPI_Datatype), Intent(In) :: mpi_type_real

        Integer :: r, n, mpierr, nd_loc, loc_row, ipos
        Integer(Kind=Int64) :: i, new_ih8, NumH_total, cumsum, ih8_max
        Integer, Allocatable, Dimension(:) :: send_cnts, recv_cnts, send_disp, recv_disp
        Integer, Allocatable, Dimension(:) :: send_row, send_col, recv_row, recv_col
        Real(type_real), Allocatable, Dimension(:) :: send_val, recv_val
        Integer, Allocatable, Dimension(:) :: tmp_pos, row_cnt, tmp_col, row_cnts_g, rank_of_row
        Real(type_real), Allocatable, Dimension(:) :: tmp_val
        Integer(Kind=Int64) :: s1
        Character(Len=16) :: timeStr

        ! MPI_Alltoallv count and displacement arrays are Int32, so the per-rank
        ! element count must fit in Int32. Abort early with a clear message rather
        ! than silently wrapping to a negative count. 
        ! TODO - account for larger than Int32. 
        If (ih8 > Int(Huge(0), Int64)) Then
            If (mype == 0) Write(*,'(A,I0,A)') &
                'RedistributeHamCSR: local ih8 (', ih8, ') exceeds Int32 limit; ' // &
                'increase npes or use fewer matrix elements.'
            Call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        End If

        If (mype == 0) Then
            Call startTimer(s1)
            Write(*,'(2X,A)') 'RedistributeHamCSR: redistributing Hamiltonian by row...'
        End If

        ! --- Step 1: global per-row element counts ---
        ! Each rank tallies its own elements by row, then MPI_Allreduce(SUM) combines them. 
        ! MPI_IN_PLACE so every rank ends up with the same complete row_cnts_g(1:Nd).
        Allocate(row_cnts_g(Nd))
        row_cnts_g = 0
        Do i = 1, ih8
            row_cnts_g(Hamil%row(i)) = row_cnts_g(Hamil%row(i)) + 1
        End Do
        Call MPI_Allreduce(MPI_IN_PLACE, row_cnts_g, Nd, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
        NumH_total = 0_Int64
        Do n = 1, Nd
            NumH_total = NumH_total + Int(row_cnts_g(n), Int64)
        End Do

        ! --- Step 2: assign contiguous row ranges to ranks ---
        ! assign rank_of_row(n) BEFORE adding row n's count to cumsum so that row n goes to the rank that was current when we started it,
        ! not the rank we may flip to after adding it. 
        ! The flip condition
        !   cumsum * npes >= (r+1) * NumH_total
        ! is the integer-safe form of cumsum/NumH_total >= (r+1)/npes.
        Allocate(rank_of_row(Nd))
        r = 0
        cumsum = 0_Int64
        Do n = 1, Nd
            rank_of_row(n) = r
            cumsum = cumsum + Int(row_cnts_g(n), Int64)
            If (r < npes - 1) Then
                If (cumsum * npes >= Int(r + 1, Int64) * NumH_total) r = r + 1
            End If
        End Do
        Deallocate(row_cnts_g)

        ! nd_start is the global-to-local offset used throughout Davidson:
        !   loc_row = global_row - nd_start   (1-based local, range 1..nd_loc)
        ! nd_start is the last global row assigned to ranks < mype (0-based exclusive),
        ! nd_end is the last global row assigned to mype (1-based inclusive = nd_start+nd_loc).
        nd_start = 0
        nd_end   = -1
        Do n = 1, Nd
            If (rank_of_row(n) <  mype) nd_start = n
            If (rank_of_row(n) == mype) nd_end   = n
        End Do
        If (nd_end < 0) nd_end = nd_start  ! rank has no rows (extremely small problem)
        nd_loc = nd_end - nd_start

        ! --- Step 3: exchange element counts ---
        ! send_cnts(r) = how many elements this rank sends to rank r.
        ! MPI_Alltoall transposes the npes×npes count matrix: rank r's recv_cnts(mype) becomes this rank's send_cnts(r). 
        ! After this call recv_cnts(r) = how many elements this rank receives from rank r.
        Allocate(send_cnts(0:npes-1), recv_cnts(0:npes-1))
        send_cnts = 0
        Do i = 1, ih8
            r = rank_of_row(Hamil%row(i))
            send_cnts(r) = send_cnts(r) + 1
        End Do
        Call MPI_Alltoall(send_cnts, 1, MPI_INTEGER, recv_cnts, 1, MPI_INTEGER, &
                          MPI_COMM_WORLD, mpierr)

        ! Prefix-sum displacements for Alltoallv (0-based byte offsets into the flat buffer).
        ! new_ih8 is the total incoming element count = this rank's post-redistribution ih8.
        Allocate(send_disp(0:npes-1), recv_disp(0:npes-1))
        send_disp(0) = 0
        recv_disp(0) = 0
        Do r = 1, npes - 1
            send_disp(r) = send_disp(r-1) + send_cnts(r-1)
            recv_disp(r) = recv_disp(r-1) + recv_cnts(r-1)
        End Do
        new_ih8 = Int(recv_disp(npes-1), Int64) + recv_cnts(npes-1)

        ! --- Steps 4-6: staggered pack and Alltoallv (col, val, row) ---
        ! Three sequential passes: pack a send buffer, exchange it, free it, then repeat.
        ! rank_of_row + Hamil%row serve as the rank-lookup key; Hamil%row is packed last
        ! so it remains available as input for the col and val passes.
        ! Peak (val pass): Hamil%row(4B) + Hamil%val(8B) + col_new(4B×new) + send_val(8B)
        !                  = 24B/element of ih8.

        ! Pass 1: col
        Allocate(send_col(ih8), tmp_pos(0:npes-1))
        tmp_pos = send_disp + 1
        Do i = 1, ih8
            r    = rank_of_row(Hamil%row(i))
            ipos = tmp_pos(r)
            send_col(ipos) = Hamil%col(i)
            tmp_pos(r) = tmp_pos(r) + 1
        End Do
        Deallocate(Hamil%col)
        Allocate(recv_col(new_ih8))
        Call MPI_Alltoallv(send_col, send_cnts, send_disp, MPI_INTEGER, &
                           recv_col,  recv_cnts, recv_disp, MPI_INTEGER, &
                           MPI_COMM_WORLD, mpierr)
        Deallocate(send_col)
        Call Move_Alloc(recv_col, Hamil%col)

        ! Pass 2: val
        tmp_pos = send_disp + 1
        Allocate(send_val(ih8))
        Do i = 1, ih8
            r    = rank_of_row(Hamil%row(i))
            ipos = tmp_pos(r)
            send_val(ipos) = Hamil%val(i)
            tmp_pos(r) = tmp_pos(r) + 1
        End Do
        Deallocate(Hamil%val)
        Allocate(recv_val(new_ih8))
        Call MPI_Alltoallv(send_val, send_cnts, send_disp, mpi_type_real, &
                           recv_val,  recv_cnts, recv_disp, mpi_type_real, &
                           MPI_COMM_WORLD, mpierr)
        Deallocate(send_val)
        Call Move_Alloc(recv_val, Hamil%val)

        ! Pass 3: row (Hamil%row is both the rank-lookup key and the data being packed)
        tmp_pos = send_disp + 1
        Allocate(send_row(ih8))
        Do i = 1, ih8
            r    = rank_of_row(Hamil%row(i))
            ipos = tmp_pos(r)
            send_row(ipos) = Hamil%row(i)
            tmp_pos(r) = tmp_pos(r) + 1
        End Do
        Deallocate(tmp_pos, rank_of_row, Hamil%row)
        Allocate(recv_row(new_ih8))
        Call MPI_Alltoallv(send_row, send_cnts, send_disp, MPI_INTEGER, &
                           recv_row,  recv_cnts, recv_disp, MPI_INTEGER, &
                           MPI_COMM_WORLD, mpierr)
        Deallocate(send_row)
        Call Move_Alloc(recv_row, Hamil%row)

        Deallocate(send_cnts, send_disp, recv_cnts, recv_disp)
        ih8 = new_ih8

        ! --- Step 7: counting sort into CSR ---
        ! loc_row = row(i) - nd_start converts a global 1-based row index to a 1-based local index (1..nd_loc).
        ! HamilCSR_rowptr is the CSR prefix sum: elements for local row i are at positions rowptr(i-1)+1 .. rowptr(i).
        Allocate(row_cnt(nd_loc))
        row_cnt = 0
        Do i = 1, ih8
            loc_row = Hamil%row(i) - nd_start
            row_cnt(loc_row) = row_cnt(loc_row) + 1
        End Do

        If (Allocated(HamilCSR_rowptr)) Deallocate(HamilCSR_rowptr)
        Allocate(HamilCSR_rowptr(0:nd_loc))
        HamilCSR_rowptr(0) = 0
        Do r = 1, nd_loc
            HamilCSR_rowptr(r) = HamilCSR_rowptr(r-1) + row_cnt(r)
        End Do
        Deallocate(row_cnt)

        ! Scatter col/val into CSR order.
        ! Two separate passes (col then val) so only one output buffer is live at a time.
        ! Peak: Hamil(16B) + tmp_col(4B) = 20B/element for col pass;
        !        Hamil%row+col_csr+val(16B) + tmp_val(8B) = 24B/element for val pass.
        ! row is freed after the val pass: rowptr already encodes all row information.
        Allocate(tmp_pos(nd_loc), tmp_col(ih8))
        tmp_pos(1:nd_loc) = HamilCSR_rowptr(0:nd_loc-1) + 1
        Do i = 1, ih8
            loc_row          = Hamil%row(i) - nd_start
            ipos             = tmp_pos(loc_row)
            tmp_col(ipos)    = Hamil%col(i)
            tmp_pos(loc_row) = ipos + 1
        End Do
        Deallocate(Hamil%col)
        Call Move_Alloc(tmp_col, Hamil%col)

        tmp_pos(1:nd_loc) = HamilCSR_rowptr(0:nd_loc-1) + 1
        Allocate(tmp_val(ih8))
        Do i = 1, ih8
            loc_row          = Hamil%row(i) - nd_start
            ipos             = tmp_pos(loc_row)
            tmp_val(ipos)    = Hamil%val(i)
            tmp_pos(loc_row) = ipos + 1
        End Do
        Deallocate(Hamil%row, Hamil%val, tmp_pos)
        Call Move_Alloc(tmp_val, Hamil%val)

        ! --- Step 8: print load balance ---
        ih8_max = ih8
        Call MPI_AllReduce(MPI_IN_PLACE, ih8_max, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)

        If (mype == 0) Then
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A,F5.1,A)') 'RedistributeHamCSR: H*v load balance: ', &
                100.0_dp * Real(NumH_total, dp) / Real(Int(npes, Int64) * ih8_max, dp), &
                '% efficiency (avg/max element count per rank)'
            Write(*,'(2X,A)') 'TIMING >>> RedistributeHamCSR done in ' // Trim(timeStr)
        End If

    End Subroutine RedistributeHamCSR

End Module matrix_io