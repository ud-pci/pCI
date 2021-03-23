Module mpi_wins
    
    Use conf_variables
    Use, intrinsic :: ISO_C_BINDING, Only : C_PTR, C_F_POINTER

    Implicit None
    
    Private
    
    Public :: CreateWindow, CloseWindow
    
  Contains

    Subroutine CreateWindow(win, mpierr)
        Use mpi
        Implicit None
        Integer :: split_type, key, disp_unit, win, mpierr
        Integer(kind=MPI_ADDRESS_KIND) :: size
        Integer(kind=int64) :: size8
        TYPE(C_PTR) :: baseptr
        Integer, Allocatable :: arrayshape(:)
        Integer :: shmrank, shmsize, shmcomm, zerokey, zerocomm, zerorank, zerosize
        Call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                           MPI_INFO_NULL, shmcomm, mpierr)
        Call MPI_COMM_rank(shmcomm, shmrank, mpierr)
        Call MPI_COMM_size(shmcomm, shmsize, mpierr)
        ! If core of each node, set zerokey to be true (=1)
        If (shmrank==0) Then
          zerokey=1
        Else
          zerokey=0
        End If
        ! create zerocomm for master cores to communicate
        Call MPI_COMM_split(MPI_COMM_WORLD,zerokey,0,zerocomm,mpierr)
        Call MPI_COMM_rank(zerocomm, zerorank, mpierr)
        Call MPI_COMM_size(zerocomm, zerosize, mpierr)

        Allocate(arrayshape(2))
        arrayshape=(/ Ne, Nd /)
        size = 0_MPI_ADDRESS_KIND
        If (shmrank == 0) Then
           size8 = Ne
           size8 = size8*Nd
           size = int(size8,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND !*8 for double ! Put the actual data size here                                                                                                                                                                
        Else
           size = 0_MPI_ADDRESS_KIND
        End If
        disp_unit = 1
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, shmcomm, baseptr, win, mpierr)
        ! Obtain the location of the memory segment
        If (shmrank /= 0) Then
          Call MPI_Win_shared_query(win, 0, size, disp_unit, baseptr, mpierr)
        End If
        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        Call C_F_POINTER(baseptr, Iarr_win,arrayshape)
        If (shmrank == 0) Then
           Iarr_win=Iarr
        End If
        Call MPI_Win_Fence(0, win, mpierr)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End Subroutine
  
    Subroutine CloseWindow(win, mpierr)
        Use mpi
        Implicit None
        Integer :: win, mpierr
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Win_Fence(0, win, mpierr)
       
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Win_Free(win, mpierr)
       
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    End Subroutine
End Module
        
