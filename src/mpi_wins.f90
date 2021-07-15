Module mpi_wins
    !
    ! This module implements subroutines that create and close MPI windows
    ! to allow all cores access to the basis set of determinants stored in Iarr
    ! 
    Use conf_variables
    Use, intrinsic :: ISO_C_BINDING, Only : C_PTR, C_F_POINTER

    Implicit None
    
    Private
    
    Public :: CreateIarrWindow, CloseIarrWindow
    
  Contains

    Subroutine CreateIarrWindow(win, mype, npes, shmrank, mpierr)
        Use mpi
        Implicit None
        Integer :: split_type, key, disp_unit, win, mpierr, i, mype, npes
        Integer(Kind=MPI_ADDRESS_KIND) :: size
        Integer(Kind=Int64) :: size8
        TYPE(C_PTR) :: baseptr
        Integer, Allocatable, Dimension(:) :: arrayshape
        Integer :: shmrank, shmsize, shmcomm, zerokey, zerocomm, zerorank, zerosize

        Call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shmcomm, mpierr)
        Call MPI_COMM_rank(shmcomm, shmrank, mpierr)
        Call MPI_COMM_size(shmcomm, shmsize, mpierr)
        ! If core of each node, set zerokey to be true (=1)
        If (shmrank==0) Then
            zerokey=1
        Else
            zerokey=0
        End If
        ! create zerocomm for master cores to communicate
        Call MPI_COMM_split(MPI_COMM_WORLD, zerokey, 0, zerocomm, mpierr)
        Call MPI_COMM_rank(zerocomm, zerorank, mpierr)
        Call MPI_COMM_size(zerocomm, zerosize, mpierr)
        If (shmrank == 0 .and. (.not. Allocated(Iarr))) then
            Allocate(Iarr(Ne,Nd)) 
        End If
        If (shmrank == 0) then
            Do i=1,Ne
                Call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, zerocomm, mpierr)
            End do   
        End If
        Allocate(arrayshape(2))
        arrayshape=(/ Ne, Nd /)
        size = 0_MPI_ADDRESS_KIND
        If (shmrank == 0) Then
            size8 = Ne
            size8 = size8*Nd
            size = Int(size8, MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND !*8 for double ! Put the actual data size here                                                                                                                                                                
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
        print*,'before',sizeof(Iarr_win)
        Call C_F_POINTER(baseptr, Iarr_win, arrayshape)
        If (shmrank == 0) Then
            Iarr_win=Iarr
        End If
        print*,'after',sizeof(Iarr_win)
        Call MPI_Win_Fence(0, win, mpierr)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Return
    End Subroutine CreateIarrWindow
  
    Subroutine CloseIarrWindow(win, mype, npes, shmrank, mpierr)
        Use mpi
        Implicit None
        Integer :: win, mype, npes, shmrank, mpierr

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Win_Fence(0, win, mpierr)
       
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Win_Free(win, mpierr)

        Nullify(Iarr_win)
        print*,'after2',sizeof(Iarr_win)
       
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    End Subroutine CloseIarrWindow
End Module
        
