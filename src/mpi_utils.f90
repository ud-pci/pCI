Module mpi_utils
    !
    !   Written by Jeffrey Frey, UD IT-RCI, 2/12/21
    !   Editted by Charles Cheung, 11/10/21
    !
    ! This module implements MPI broadcast functions that break a buffer
    ! that may be too large into multiple chunks.  The <stride> passed
    ! to each of the subroutines dictates how many 8-byte elements should
    ! be passed to each underlying MPI_Bcast().
    !
    ! If the <stride> is zero, then the environment variable CONF_BROADCAST_MAX
    ! is checked for an integer value.  If no value is found, or the value
    ! exceeds the maximum threshold, the maximum threshold is used.  The
    ! maximum equates to 1 GiB worth of elements.
    !
    Implicit None
    
    Private
    
    Public :: BroadcastI, BroadcastR, BroadcastD

    Integer, Parameter :: MaxStride8Byte = 134217728
    Integer, Parameter :: MaxStride4Byte = 67108864
    
  Contains

    Subroutine BroadcastI(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Implicit None

        Type(MPI_Comm), Intent(In)                  :: comm
        Integer, Intent(In)                         :: count, stride, rank
        Integer, Dimension(*), Intent(InOut)        :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer                                     :: use_stride, count_remain, i

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride4Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride4Byte * 2) use_stride = MaxStride4Byte * 2
        count_remain = count
        i = 1
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_Bcast(buffer(i:i+use_stride), use_stride, MPI_INTEGER, rank, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting integer range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            Call MPI_Bcast(buffer(i:i+count_remain), count_remain, MPI_INTEGER, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting integer range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastI

    Subroutine BroadcastR(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Implicit None

        Type(MPI_Comm), Intent(In)                  :: comm
        Integer, Intent(In)                         :: count, stride, rank
        Real, Dimension(*), Intent(InOut)           :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer                                     :: use_stride, count_remain, i

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count

        i = 1
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_Bcast(buffer(i:i+use_stride), use_stride, MPI_REAL, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            Call MPI_Bcast(buffer(i:i+count_remain), count_remain, MPI_REAL, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastR

    Subroutine BroadcastD(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Use, intrinsic :: iso_fortran_env
        Implicit None

        Type(MPI_Comm), Intent(In)                  :: comm
        Integer, Intent(In)                         :: count, stride, rank
        Real(real64), Dimension(*), Intent(InOut)   :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer                                     :: use_stride, count_remain, i

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte) use_stride = MaxStride8Byte
        count_remain = count
        i = 1
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_Bcast(buffer(i:i+use_stride), use_stride, MPI_DOUBLE_PRECISION, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            Call MPI_Bcast(buffer(i:i+count_remain), count_remain, MPI_DOUBLE_PRECISION, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastD
  
End Module mpi_utils