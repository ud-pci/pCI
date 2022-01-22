Module matrix_io
    ! 
    ! This module implements MPI functions that read and write matrix elements to disc. 
    !
    Use conf_variables
    Use str_fmt, Only : startTimer, stopTimer, FormattedTime

    Implicit None

    Private

    Public :: WriteMatrix, ReadMatrix

  Contains

    Subroutine MPIErrHandle(mpierr)
        Use mpi
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

    Subroutine WriteMatrix(matn,matk,matt,num_elements_per_core,num_elements_total,filename,mype,npes,mpierr)
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
        Use mpi
        Implicit None
        Integer, Allocatable, Dimension(:)  :: matn, matk
        Real(dp), Allocatable, Dimension(:) :: matt
        Integer, Intent(In)             :: num_elements_per_core
        Integer(Kind=int64), Intent(In) :: num_elements_total

        Integer                       :: mype, npes, i, mpierr, fh
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
            disps(i)=disps(i-1)+sizes(i-1)
        End Do
    
        ! Open file
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write number of matrix elements
        disp = mype * 4_MPI_OFFSET_KIND
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPIErrHandle(mpierr)
        Call MPI_FILE_WRITE(fh, num_elements_per_core, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%n
        disp = (npes + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPIErrHandle(mpierr)
        Call MPI_FILE_WRITE(fh, matn, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%k
        disp = (npes + num_elements_total + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPIErrHandle(mpierr)
        Call MPI_FILE_WRITE(fh, matk, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Write mat%t
        disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 8_MPI_OFFSET_KIND
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPIErrHandle(mpierr)
        Call MPI_FILE_WRITE(fh, matt, num_elements_per_core, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Print time to write file
        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> Writing ' // filename // ' took ' // trim(timeStr)
        End If
    End Subroutine WriteMatrix

    Subroutine ReadMatrix(matn,matk,matt,num_elements_per_core,num_elements_total,filename,mype,npes,mpierr)
        !
        ! This subroutine reads a parallel matrix file 
        !
        Use mpi
        Implicit None
        Integer, Allocatable, Dimension(:)  :: matn, matk
        Real(dp), Allocatable, Dimension(:) :: matt
        Integer, Intent(Out)             :: num_elements_per_core
        Integer(Kind=int64), Intent(Out) :: num_elements_total

        Integer                                         :: mype, npes, npes_read, i, mpierr, fh
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
            Open(66,file='progress.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Read(66) Nc_prev, Nd_prev
            Close(66)  
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
        If (.not. allocated(matn)) allocate(matn(num_elements_per_core))
        If (.not. allocated(matk)) allocate(matk(num_elements_per_core))
        If (.not. allocated(matt)) allocate(matt(num_elements_per_core))

        ! Read first indices of matrix elements mat%n
        disp = (npes + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, matn, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Read second indices of matrix elements mat%k
        disp = (npes + num_elements_total + disps(mype+1)) * 4_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, matk, num_elements_per_core, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Read values of matrix element mat%t
        disp = (npes + 2 * num_elements_total) * 4_MPI_OFFSET_KIND + disps(mype+1) * 8_MPI_OFFSET_KIND
        Call MPI_FILE_READ_AT(fh, disp, matt, num_elements_per_core, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Call MPIErrHandle(mpierr)

        ! Print time to write file
        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> Reading ' // filename // ' took ' // trim(timeStr)
        End If
    
        Return
    End Subroutine ReadMatrix

End Module matrix_io