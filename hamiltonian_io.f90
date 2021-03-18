Module hamiltonian_io
    ! 
    ! This module implements MPI functions that read and write matrix elements 
    ! to disc. 
    !
    Use conf_variables

    Implicit None

    Private

    Public :: Hwrite_s, Hwrite, Hread, Jwrite

  Contains

    Subroutine Hwrite_s(mype,npes,counter)
        Use mpi
        Implicit None
    
        Integer                  :: mype, npes, counter, i, j, k, mpierr
        Integer(kind=int64)      :: disp
        Integer, dimension(npes) :: sizes

        disp=0_int64
        ! Write counters
        if (mype==0) then
          open(unit=30,file='CONF.HIJs', &
              status='UNKNOWN',form='unformatted')
          close(30,status='delete')
          open(unit=30,file='CONF.HIJs', &
              status='NEW',form='unformatted')
          write(30) npes
        end if
        disp=disp+1
        do i=1,npes
          if (mype==i-1) then
            open(unit=30,file='CONF.HIJs', &
            status='UNKNOWN',position='append',form='unformatted')
            write(30) counter
            do j=1,counter
              write(30) H_n(j),H_k(j),H_t(j)
            end do
            close(30)
          end if
          Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        end do
        Return
    End Subroutine Hwrite_s

    Subroutine Hwrite (mype,npes,counter)
        !
        ! This subroutine writes CONF.HIJ
        ! The structure of CONF.HIJ is as follows:
        ! =======================================================
        ! | number of processors  | counters for each processor | 
        ! | H_n (core 0) | H_n (core 1) | ... | H_n (core npes) |
        ! | H_k (core 0) | H_k (core 1) | ... | H_k (core npes) |
        ! | H_t (core 0) | H_t (core 1) | ... | H_t (core npes) |
        ! =======================================================
        !
        Use mpi
        Implicit None

        Integer                       :: mype, npes, counter, i, mpierr, fh
        Integer, dimension(npes)      :: sizes, disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=16)             :: filename
    
        sizes=0
        Call MPI_AllGather(counter, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        disps=0
        do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        end do
    
        Write(filename,'(A)') 'CONF.HIJ'
        ! Write counters
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, fh, mpierr) 

        ! Write number of processors
        disp=0
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        if (mype == 0) then
            Call MPI_FILE_WRITE(fh, npes, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 
        end if
        ! Write counters
        disp = mype * 4 + 4

        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, counter, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        disp = npes * 4 + disps(mype+1) * 4 + 4
        
        ! Write H_n
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_n, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        disp = npes * 4 + NumH * 4 + disps(mype+1) * 4 + 4

        ! Write H_k
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_k, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        ! Write H_t
        disp = npes * 4 + 2 * NumH * 4 + disps(mype+1) * 8 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_t, counter, MPI_DOUBLE_PRECISION, & 
                            MPI_STATUS_IGNORE, mpierr) 
        
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Return
    End Subroutine Hwrite

    subroutine Hread (mype,npes,counter)
        !
        ! This subroutine reads CONF.HIJ
        !
        Use mpi
        Implicit None
    
        Integer                       :: mype, npes, npes_read, counter, i, mpierr, fh
        Integer(kind=int64)           :: ih8
        Integer, dimension(npes)      :: sizes, disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=16)             :: filename

        Write(filename,'(A)') 'CONF.HIJ'
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, & 
                           MPI_MODE_RDONLY, & 
                           MPI_INFO_NULL, fh, mpierr) 

        ! Read number of processors
        disp = 0
        if (mype == 0) then
            Call MPI_FILE_READ_AT(fh, disp, npes_read, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr)
            if (npes /= npes_read) then
                print*, 'Number of processors inconsistent. CONF.HIJ was written with ', npes_read, ' processors,', &
                        ' but', npes, ' is available.'
                stop
            end if
        end if

        ! Read counters
        disp = mype * 4 + 4
        Call MPI_FILE_READ_AT(fh, disp, counter, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr)

        if (.not. allocated(H_n)) allocate(H_n(counter))
        if (.not. allocated(H_k)) allocate(H_k(counter))
        if (.not. allocated(H_t)) allocate(H_t(counter))
        ! Calculate displacements
        sizes=0
        Call MPI_AllGather(counter, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        ih8=counter
        Call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)

        disps=0
        do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        end do

        ! Read H_n
        disp = npes * 4 + disps(mype+1) * 4 + 4

        Call MPI_FILE_READ_AT(fh, disp, H_n, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        ! Read H_k
        disp = npes * 4 + NumH * 4 + disps(mype+1) * 4 + 4

        Call MPI_FILE_READ_AT(fh, disp, H_k, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        ! Read H_t
        disp = npes * 4 + 2 * NumH * 4 + disps(mype+1) * 8 + 4

        Call MPI_FILE_READ_AT(fh, disp, H_t, counter, MPI_DOUBLE_PRECISION, & 
                            MPI_STATUS_IGNORE, mpierr) 
        
        Call MPI_FILE_CLOSE(fh, mpierr) 
    
        return
    end subroutine Hread

    subroutine Jwrite (mype,npes,counter)
        !
        ! This subroutine writes CONF.JJJ
        ! The structure of CONF.JJJ is as follows:
        ! =======================================================
        ! | number of processors  | counters for each processor | 
        ! | J_n (core 0) | J_n (core 1) | ... | J_n (core npes) |
        ! | J_k (core 0) | J_k (core 1) | ... | J_k (core npes) |
        ! | J_t (core 0) | J_t (core 1) | ... | J_t (core npes) |
        ! =======================================================
        !
        Use mpi
        Implicit None
    
        Integer                       :: mype, npes, counter, i, mpierr, fh
        Integer, dimension(npes)      :: sizes, disps
        Integer(kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=16)             :: filename
    
        sizes=0
        Call MPI_AllGather(counter, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
        disps=0
        do i=2,npes
          disps(i)=disps(i-1)+sizes(i-1)
        end do
    
        Write(filename,'(A)') 'CONF.HIJ'
        ! Write counters
        Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, fh, mpierr) 

        ! Write number of processors
        disp=0
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        if (mype == 0) then
            Call MPI_FILE_WRITE(fh, npes, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 
        end if
        ! Write counters
        disp = mype * 4 + 4

        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, counter, 1, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        disp = npes * 4 + disps(mype+1) * 4 + 4
        
        ! Write H_n
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_n, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        disp = npes * 4 + NumH * 4 + disps(mype+1) * 4 + 4

        ! Write H_k
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_k, counter, MPI_INTEGER, & 
                            MPI_STATUS_IGNORE, mpierr) 

        ! Write H_t
        disp = npes * 4 + 2 * NumH * 4 + disps(mype+1) * 8 + 4
        Call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, & 
                               MPI_INTEGER, 'native', & 
                               MPI_INFO_NULL, mpierr) 
        Call MPI_FILE_WRITE(fh, H_t, counter, MPI_DOUBLE_PRECISION, & 
                            MPI_STATUS_IGNORE, mpierr) 
        
        Call MPI_FILE_CLOSE(fh, mpierr) 
        Return
    End Subroutine Jwrite

End Module hamiltonian_io