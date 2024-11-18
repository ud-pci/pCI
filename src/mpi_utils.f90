Module mpi_utils    
    !
    !   Written by Jeffrey Frey, UD IT-RCI, 2/12/21
    !   Editted by Charles Cheung, 07/03/24    
    !
    ! This module implements MPI broadcast and allreduce functions that break a buffer
    ! that may be too large into multiple chunks.  The <stride> passed
    ! to each of the subroutines dictates how many 8-byte elements should
    ! be passed to each underlying MPI_Bcast() or MPI_AllReduce().
    !
    ! If the <stride> is zero, then the environment variable CONF_BROADCAST_MAX
    ! is checked for an integer value.  If no value is found, or the value
    ! exceeds the maximum threshold, the maximum threshold is used.  The
    ! maximum equates to 1 GiB worth of elements.
    !
    Use params

    Implicit None
    
    Private
    
    Public :: BroadcastI, BroadcastR, BroadcastD, AllReduceI, AllReduceR, AllReduceD

    Integer, Parameter :: MaxStride8Byte = 134217728
    
  Contains

    Subroutine BroadcastI(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                  :: comm
        Integer(Kind=int64), Intent(In)             :: count
        Integer, Intent(In)                         :: stride, rank
        Integer, Dimension(*), Intent(InOut)        :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer(Kind=int64)                         :: count_remain, i
        Integer                                     :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count
        i = 1_int64
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
            count_remain2 = count_remain
            Call MPI_Bcast(buffer(i:i+count_remain2), count_remain2, MPI_INTEGER, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting integer range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastI

    Subroutine BroadcastR(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                  :: comm
        Integer(Kind=int64), Intent(In)             :: count
        Integer, Intent(In)                         :: stride, rank
        Real, Dimension(*), Intent(InOut)           :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer(Kind=int64)                         :: count_remain, i
        Integer                                     :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count

        i = 1_int64
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
            count_remain2 = count_remain
            Call MPI_Bcast(buffer(i:i+count_remain2), count_remain2, MPI_REAL, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastR

    Subroutine BroadcastD(buffer, count, stride, rank, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                  :: comm
        Integer(Kind=int64), Intent(In)             :: count
        Integer, Intent(In)                         :: stride, rank
        Real(type2_real), Dimension(*), Intent(InOut)   :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer(Kind=int64)                         :: count_remain, i
        Integer                                     :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) Then
                use_stride = MaxStride8Byte
                If (type2_real == sp) use_stride = use_stride * 2
            End If
        Else
            use_stride = stride
        End If
        Select Case(type2_real)
        Case(sp)
            If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        Case(dp)
            If (stride .gt. MaxStride8Byte) use_stride = MaxStride8Byte
        End Select
        count_remain = count
        i = 1_int64
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_Bcast(buffer(i:i+use_stride), use_stride, mpi_type2_real, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            count_remain2 = count_remain
            Call MPI_Bcast(buffer(i:i+count_remain2), count_remain2, mpi_type2_real, rank, comm, mpierr)
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure broadcasting real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine BroadcastD

    Subroutine AllReduceI(buffer, count, stride, op, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                  :: comm
        Type(MPI_Op), Intent(In)                    :: op
        Integer(Kind=int64), Intent(In)             :: count
        Integer, Intent(In)                         :: stride
        Integer, Dimension(*), Intent(InOut)        :: buffer
        Integer, Intent(InOut)                      :: mpierr
        Character(Len=255)                          :: envvar
        Integer(Kind=int64)                         :: count_remain, i
        Integer                                     :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count
        i = 1_int64
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+use_stride), use_stride, MPI_INTEGER, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing integer range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            count_remain2 = count_remain
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+count_remain2), count_remain2, MPI_INTEGER, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing integer range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine AllReduceI

    Subroutine AllReduceR(buffer, count, stride, op, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                     :: comm
        Type(MPI_Op), Intent(In)                       :: op
        Integer(Kind=int64), Intent(In)                :: count
        Integer, Intent(In)                            :: stride
        Real(type2_real), Dimension(*), Intent(InOut)  :: buffer
        Integer, Intent(InOut)                         :: mpierr
        Character(Len=255)                             :: envvar
        Integer(Kind=int64)                            :: count_remain, i
        Integer                                        :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count
        i = 1_int64
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+use_stride), use_stride, mpi_type2_real, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing real range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            count_remain2 = count_remain
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+count_remain2), count_remain2, mpi_type2_real, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine AllReduceR

    Subroutine AllReduceD(buffer, count, stride, op, comm, mpierr)
        Use mpi_f08
        Implicit None

        EXTERNAL :: GetEnv
        Type(MPI_Comm), Intent(In)                     :: comm
        Type(MPI_Op), Intent(In)                       :: op
        Integer(Kind=int64), Intent(In)                :: count
        Integer, Intent(In)                            :: stride
        Real(dp), Dimension(*), Intent(InOut)          :: buffer
        Integer, Intent(InOut)                         :: mpierr
        Character(Len=255)                             :: envvar
        Integer(Kind=int64)                            :: count_remain, i
        Integer                                        :: count_remain2, use_stride

        If (stride .le. 0) Then
            Call GetEnv('CONF_BROADCAST_MAX', envvar)
            Read(envvar, '(i)') use_stride
            If (use_stride .le. 0) use_stride = MaxStride8Byte * 2
        Else
            use_stride = stride
        End If
        If (stride .gt. MaxStride8Byte * 2) use_stride = MaxStride8Byte * 2
        count_remain = count
        i = 1_int64
        mpierr = 0
        Do While (count_remain .gt. use_stride)
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+use_stride), use_stride, MPI_DOUBLE_PRECISION, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing real range ',i,':',i+use_stride
                Return
            End If
            count_remain = count_remain - use_stride
            i = i + use_stride
        End Do
        If (count_remain .gt. 0) Then
            count_remain2 = count_remain
            Call MPI_AllReduce(MPI_IN_PLACE, buffer(i:i+count_remain2), count_remain2, MPI_DOUBLE_PRECISION, op, comm, mpierr) 
            If (mpierr .ne. 0 ) Then
                Write(*,*) 'Failure reducing real range ',i,':',i+use_stride
                Return
            End If
        End If
    End Subroutine AllReduceD

End Module mpi_utils