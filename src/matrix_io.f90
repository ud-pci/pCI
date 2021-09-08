Module matrix_io
    ! 
    ! This module implements MPI functions that read and write matrix elements to disc. 
    !
    Use conf_variables
    Use str_fmt, Only : FormattedTime

    Implicit None

    Private

    Public :: WriteMatrix, ReadMatrix

  Contains

    Subroutine WriteMatrix(mat,num_elements,num_total,filename,mype,npes,mpierr)
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

        Integer                       :: mype, npes, num_elements, i, mpierr, fh
        Integer(Kind=Int64)           :: num_total, clock_rate, start_time, end_time
        Real                          :: total_time
        Integer, dimension(npes)      :: sizes, disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=*)              :: filename
        Character(Len=16)             :: timeStr
        Type(Matrix)                  :: mat
    
        If (mype == 0) Then            
            Open(66,file='nprocs.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Write(66) npes
            Close(66)
            print*, 'Writing ' // filename // '...'
            Call system_clock(count_rate=clock_rate)
            Call system_clock(start_time)
        End If

        sizes=0
        Call MPI_AllGather(num_elements, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        disps=0
        Do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        End Do
    
        ! Open file
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, mpierr) 

        ! Write number of processors
        disp=0_MPI_OFFSET_KIND
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        If (mype == 0) then
            Call MPI_FILE_WRITE(fh, npes, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 
        End If

        ! Write number of matrix elements
        disp = mype * 4 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, num_elements, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 

        ! Write mat%n
        disp = npes * 4 + disps(mype+1) * 4 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, mat%n, num_elements, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 

        ! Write mat%k
        disp = npes * 4 + num_total * 4 + disps(mype+1) * 4 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, mat%k, num_elements, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 

        ! Write mat%t
        disp = npes * 4 + 2 * num_total * 4 + disps(mype+1) * 8 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, mat%t, num_elements, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
        
        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 

        ! Print time to write file
        If (mype == 0) Then
            Call system_clock(end_time)
            total_time=Real((end_time-start_time)/clock_rate)
            Call FormattedTime(total_time, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> Writing ' // filename // ' took ' // trim(timeStr)
        End If

        Return
    End Subroutine WriteMatrix

    Subroutine ReadMatrix(mat,num_elements,num_total,filename,mype,npes,mpierr)
        !
        ! This subroutine reads a parallel matrix file 
        !
        Use mpi
        Implicit None
    
        Integer                       :: mype, npes, npes_read, num_elements, i, mpierr, fh
        Integer(kind=int64)           :: num_total, clock_rate, start_time, end_time
        Real                          :: total_time
        Integer, dimension(npes)      :: sizes, disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=*)              :: filename
        Character(Len=16)             :: timeStr
        Type(Matrix)                  :: mat

        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, mpierr) 

        ! Read number of processors
        disp = 0
        If (mype == 0) then
            Open(66,file='nprocs.conf',status='UNKNOWN',form='UNFORMATTED',access='stream')
            Read(66) npes_read
            Close(66)
            If (npes /= npes_read) then
                print*, 'Number of processors inconsistent. '// filename //' was written with ', npes_read, ' processors,', &
                        ' but', npes, ' is available.'
                stop
            End If
            print*, 'Reading ' // filename // '...'
            Call system_clock(count_rate=clock_rate)
            Call system_clock(start_time)
        End If

        ! Read counters
        disp = mype * 4 + 4
        Call MPI_FILE_READ_AT(fh, disp, num_elements, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr)

        ! Allocate matrix if not allocated yet
        If (.not. allocated(mat%n)) allocate(mat%n(num_elements))
        If (.not. allocated(mat%k)) allocate(mat%k(num_elements))
        If (.not. allocated(mat%t)) allocate(mat%t(num_elements))

        ! Calculate displacements
        sizes=0
        Call MPI_AllGather(num_elements, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(num_elements, num_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
        disps=0
        Do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        End Do

        ! Read first indices of matrix elements mat%n
        disp = npes * 4 + disps(mype+1) * 4 + 4
        Call MPI_FILE_READ_AT(fh, disp, mat%n, num_elements, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 

        ! Read second indices of matrix elements mat%k
        disp = npes * 4 + num_total * 4 + disps(mype+1) * 4 + 4
        Call MPI_FILE_READ_AT(fh, disp, mat%k, num_elements, MPI_INTEGER, MPI_STATUS_IGNORE, mpierr) 

        ! Read values of matrix element mat%t
        disp = npes * 4 + 2 * num_total * 4 + disps(mype+1) * 8 + 4
        Call MPI_FILE_READ_AT(fh, disp, mat%t, num_elements, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierr) 
        
        ! Close file
        Call MPI_FILE_CLOSE(fh, mpierr) 

        ! Print time to write file
        If (mype == 0) Then
            Call system_clock(end_time)
            total_time=Real((end_time-start_time)/clock_rate)
            Call FormattedTime(total_time, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> Reading ' // filename // ' took ' // trim(timeStr)
        End If
    
        Return
    End Subroutine ReadMatrix

End Module matrix_io