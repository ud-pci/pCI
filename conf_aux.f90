module conf_aux

    use conf_variables

    implicit none

    contains

    subroutine AllocateFormHArrays(mype, npes)
      use mpi
      Use str_fmt, Only : FormattedMemSize
      implicit none
      integer :: mpierr, mype, npes
      character(len=16) :: memStr
      if (mype==0) deallocate(Qnl)
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(IPlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(NhintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ngint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(NgintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      if (.not. allocated(Nvc)) allocate(Nvc(Nc))
      if (.not. allocated(Nc0)) allocate(Nc0(Nc))
      if (.not. allocated(Ndc)) allocate(Ndc(Nc))
      if (.not. allocated(Jz)) allocate(Jz(Nst))
      if (.not. allocated(Nh)) allocate(Nh(Nst))
      if (.not. allocated(Diag)) allocate(Diag(Nd))
      if (.not. allocated(Rint1)) allocate(Rint1(Nhint))
      if (.not. allocated(Rint2)) allocate(Rint2(IPbr,Ngint))
      if (.not. allocated(Iint1)) allocate(Iint1(Nhint))
      if (.not. allocated(Iint2)) allocate(Iint2(Ngint))
      if (.not. allocated(Iint3)) allocate(Iint3(Ngint))
      if (.not. allocated(IntOrd)) allocate(IntOrd(nrd))
      if (Ksig /= 0) then
        if (.not. allocated(Rsig)) allocate(Rsig(NhintS))
        if (.not. allocated(Dsig)) allocate(Dsig(NhintS))
        if (.not. allocated(Esig)) allocate(Esig(NhintS)) 
        if (.not. allocated(R_is)) allocate(R_is(num_is))
        if (.not. allocated(I_is)) allocate(I_is(num_is))
        if (.not. allocated(Rint2S)) allocate(Rint2S(NgintS))
        if (.not. allocated(Dint2S)) allocate(Dint2S(NgintS))
        if (.not. allocated(Eint2S)) allocate(Eint2S(NgintS))
        if (.not. allocated(Iint1S)) allocate(Iint1S(NhintS))
        if (.not. allocated(Iint2S)) allocate(Iint2S(NgintS))
        if (.not. allocated(Iint3S)) allocate(Iint3S(NgintS))
        if (.not. allocated(IntOrdS)) allocate(IntOrdS(nrd))
      end if
      if (mype==0) then
        memFormH = 0_int64
        memFormH = sizeof(Nvc)+sizeof(Nc0) &
            + sizeof(Rint1)+sizeof(Rint2)+sizeof(Iint1)+sizeof(Iint2)+sizeof(Iint3)+sizeof(Iarr)&
            + sizeof(IntOrd)
        if (Ksig /= 0) memFormH = memFormH+sizeof(Rint2S)+sizeof(Dint2S)+sizeof(Eint2S) &
            + sizeof(Iint1S)+sizeof(Iint2S)+sizeof(Iint3S) &
            + sizeof(Rsig)+sizeof(Dsig)+sizeof(Esig)+sizeof(R_is)+sizeof(I_is)+sizeof(IntOrdS)
        Call FormattedMemSize(memFormH, memStr)
        Write(*,'(A,A,A)') 'Allocating arrays for FormH requires ',Trim(memStr),' of memory per core' 
      end if   
      return
    end subroutine AllocateFormHArrays

    subroutine DeAllocateFormHArrays(mype, npes)
      use mpi
      Use str_fmt, Only : FormattedMemSize
      implicit none
      integer :: mpierr, mype, npes
      character(len=16) :: memStr
      integer(kind=8) :: mem 
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)

      if (mype==0) then
        Call FormattedMemSize(memFormH, memStr)
        Write(*,'(A,A,A)') 'De-allocating ',Trim(memStr),' of memory per core from arrays for FormH'  
      end if

      if (allocated(Nvc)) deallocate(Nvc)
      if (allocated(Nc0)) deallocate(Nc0)
      if (allocated(Rint1)) deallocate(Rint1)
      if (allocated(Rint2)) deallocate(Rint2)
      if (allocated(Iint1)) deallocate(Iint1)
      if (allocated(Iint2)) deallocate(Iint2)
      if (allocated(Iint3)) deallocate(Iint3)
      if (allocated(Rint2S)) deallocate(Rint2S)
      if (allocated(Dint2S)) deallocate(Dint2S)
      if (allocated(Eint2S)) deallocate(Eint2S)
      if (allocated(Iint1S)) deallocate(Iint1S)
      if (allocated(Iint2S)) deallocate(Iint2S)
      if (allocated(Iint3S)) deallocate(Iint3S)
      if (allocated(Rsig)) deallocate(Rsig)
      if (allocated(Dsig)) deallocate(Dsig)
      if (allocated(Esig)) deallocate(Esig) 
      if (allocated(R_is)) deallocate(R_is)
      if (allocated(I_is)) deallocate(I_is)
      if (allocated(IntOrd)) deallocate(IntOrd)
      if (allocated(IntOrdS)) deallocate(IntOrdS)
      if (allocated(Iarr)) deallocate(Iarr)

      !if (mype==0) then
      !  print*,'Nvc',sizeof(Nvc)
      !  print*,'Nc0',sizeof(Nc0)
      !  print*,'Rint1',sizeof(Rint1)
      !  print*,'Rint2',sizeof(Rint1)
      !  print*,'Iint1',sizeof(Iint1)
      !  print*,'Iint2',sizeof(Iint2)
      !  print*,'Iint3',sizeof(Iint3)
      !  print*,'Rint2S',sizeof(Rint2S)
      !  print*,'Dint2S',sizeof(Dint2S)
      !  print*,'Eint2S',sizeof(Eint2S)
      !  print*,'Iint1S',sizeof(Iint1S)
      !  print*,'Iint2S',sizeof(Iint2S)
      !  print*,'Iint3S',sizeof(Iint3S)
      !  print*,'Rsig',sizeof(Rsig)
      !  print*,'Dsig',sizeof(Dsig)
      !  print*,'Esig',sizeof(Esig) 
      !  print*,'R_is',sizeof(R_is)
      !  print*,'I_is',sizeof(I_is)
      !  print*,'IntOrd',sizeof(IntOrd)
      !  print*,'IntOrdS',sizeof(IntOrdS)
      !  print*,'Iarr',sizeof(Iarr)
      !end if
     return
    end subroutine DeAllocateFormHArrays

    subroutine AllocateDvdsnArrays(mype, npes)
      use mpi
      Use str_fmt, Only : FormattedMemSize
      implicit none
      integer :: mpierr, mype, npes
      character(len=16) :: memStr
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      if (.not. allocated(ArrB)) allocate(ArrB(Nd,IPlv))
      if (.not. allocated(Tk)) allocate(Tk(Nlv))
      if (.not. allocated(Tj)) allocate(Tj(Nlv))
      if (.not. allocated(P)) allocate(P(IPlv,IPlv))
      if (.not. allocated(D)) allocate(D(IPlv))
      if (.not. allocated(E)) allocate(E(IPlv))
      if (.not. allocated(Iconverge)) allocate(Iconverge(IPlv))
      if (.not. allocated(B1)) allocate(B1(Nd))
      if (.not. allocated(B2)) allocate(B2(Nd))
      if (.not. allocated(Z1)) allocate(Z1(Nd0,Nd0))
      if (.not. allocated(D1)) allocate(D1(Nd0))
      if (.not. allocated(E1)) allocate(E1(Nd0))
      if (mype==0) then
        memDvdsn = 0_int64
        memDvdsn = sizeof(ArrB)+sizeof(Tk)+sizeof(Tj)+sizeof(P)+sizeof(D)+sizeof(E) &
            + sizeof(Iconverge)+sizeof(B1)+sizeof(B2)+sizeof(Z1)+sizeof(D1)+sizeof(E1)
        memDvdsn = memDvdsn+maxJcore*16_int64 ! 16 bytes to store J_n, J_k, J_t; maxJcore is max # of m.e. for a core
        !print*,'ArrB',sizeof(ArrB)
        !print*,'Tk',sizeof(Tk)
        !print*,'Tj',sizeof(Tj)
        !print*,'P',sizeof(P)
        !print*,'D',sizeof(D)
        !print*,'E',sizeof(E)
        !print*,'Iconverge',sizeof(Iconverge)
        !print*,'B1',sizeof(B1)
        !print*,'B2',sizeof(B2)
        !print*,'Z1',sizeof(Z1)
        !print*,'D1',sizeof(D1)
        !print*,'E1',sizeof(E1)
        !print*,'J_n',sizeof(J_n)
        !print*,'J_k',sizeof(J_k)
        !print*,'J_t',sizeof(J_t)
        Call FormattedMemSize(memDvdsn, memStr)
        Write(*,'(A,A,A)') 'Allocating arrays for Davidson procedure requires ',Trim(memStr),' of memory per core'  
        memEstimate = memEstimate - memFormH + memDvdsn
        Call FormattedMemSize(memEstimate, memStr)
        Write(*,'(A,A,A)') 'Total memory estimate for Davidson procedure is ',Trim(memStr),' of memory per core' 
      end if   
      return
    end subroutine AllocateDvdsnArrays

    subroutine calcMemStaticArrays
      Use str_fmt, Only : FormattedMemSize
      implicit none
      character(len=16) :: memStr

      memStaticArrays = 0_int64
      memStaticArrays = sizeof(Nn)+sizeof(Kk)+sizeof(Ll)+sizeof(Jj)+sizeof(Nf0)+sizeof(Jt)+sizeof(Njt) &
                      + sizeof(Eps)+sizeof(Diag)+sizeof(Ndc)+sizeof(Jz)+sizeof(Nh) &
                      + sizeof(In)+sizeof(Gnt)+sizeof(Scr)+sizeof(C)+sizeof(Er) 
      !print*,'Nn',sizeof(Nn)
      !print*,'Kk',sizeof(Kk)
      !print*,'Ll',sizeof(Ll)
      !print*,'Jj',sizeof(Jj)
      !print*,'Nf0',sizeof(Nf0)
      !print*,'Jt',sizeof(Jt)
      !print*,'Njt',sizeof(Njt)
      !print*,'Eps',sizeof(Eps)
      !print*,'Diag',sizeof(Diag)
      !print*,'Ndc',sizeof(Ndc)
      !print*,'Jz',sizeof(Jz)
      !print*,'Nh',sizeof(Nh)
      !print*,'In',sizeof(In)
      !print*,'Gnt',sizeof(Gnt)
      !print*,'Scr',sizeof(Scr)
      !print*,'C',sizeof(C)
      !print*,'Er',sizeof(Er)
      Call FormattedMemSize(memStaticArrays, memStr)
      Write(*,'(A,A,A)') 'calcMemReqs: Allocating static arrays requires ',Trim(memStr),' of memory per core' 
    end subroutine calcMemStaticArrays

    subroutine calcMemReqs
      Use str_fmt, Only : FormattedMemSize
      implicit none
      integer(kind=int64) :: mem
      integer :: bytesInteger, bytesDP, bytesReal
      character(len=16) :: memStr

      bytesInteger = 4
      bytesReal = 4
      bytesDP = 8

      call calcMemStaticArrays
      memEstimate = memFormH + memStaticArrays
      Call FormattedMemSize(memEstimate, memStr)
      Write(*,'(A,A,A)') 'calcMemReqs: Total memory estimate before FormH is ',Trim(memStr),' of memory per core' 

      mem = bytesDP * Nd * IPlv     & ! ArrB
          + bytesDP * Nlv * 2_dp    & ! Tk,Tj
          + bytesDP * IPlv ** 2_dp  & ! P
          + bytesDP * IPlv * 2_dp   & ! D,E
          + bytesDP * Nd * 2_dp     & ! B1,B2
          + bytesDP * Nd0 ** 2_dp   & ! Z1
          + bytesDP * Nd0 * 2_dp      ! D1,E1

      Call FormattedMemSize(mem, memStr)
      Write(*,'(A,A,A)') 'calcMemReqs: Allocating arrays for Davidson procedure will require ',Trim(memStr),' of memory per core' 

      Call FormattedMemSize(memTotalPerCPU, memStr)
      Write(*,'(A,A,A)') 'calcMemReqs: Total memory available to the job is ',Trim(memStr),' of memory per core' 

    end subroutine calcMemReqs

    subroutine FormH(npes, mype)
      use mpi
      Use str_fmt, Only : FormattedMemSize, FormattedTime
      use, intrinsic :: ISO_C_BINDING, Only : C_PTR, C_F_POINTER
      use conf_init, Only : InitFormH, InitMxmpy ! intialization subroutines for FormH and Mxmpy
      use determinants, Only : Gdet, Gdet_win, CompCD, CompD, Rspq, Rspq_phase1, Rspq_phase2
      use hamiltonian_io
      use vaccumulator
      implicit none
      integer :: npes, mype, mpierr, interval, remainder, Hlim, numBins, cnt, bincnt, npes_read
      integer :: is, nf, i1, i2, j1, j2, k1, kx, n, ic, ierr, n1, n2, int, split_type, key, disp_unit, win, &
                 n0, jq, jq0, iq, i, j, icomp, k, ih4, counter1, counter2, counter3, diff, k2, totsize
      integer :: nn, kk
      logical :: finished
      integer, allocatable, dimension(:) :: idet1, idet2, mepd
      integer, dimension(npes) :: avgs, start1, end1
      integer, dimension(npes) :: sizes, disps
      Integer(kind=int64)     :: start_time, end_time, s1, e1, s2, e2, clock_rate
      real :: ttime
      real(dp)  :: t
      integer(kind=int64) :: ih8, i8, l8, sumd, mem, memsum, ih, cntr, maxme, avgme, numme, size8
      integer, allocatable, dimension(:) :: nk
      Character(Len=16)     :: memStr, memStr2, npesStr, counterStr
      integer :: iSign, iIndexes(3), jIndexes(3)
      Type(IVAccumulator)   :: iva1, iva2
      Type(RVAccumulator)   :: rva1
      integer, Parameter    :: vaGrowBy = 10000000
      integer(kind=MPI_ADDRESS_KIND) :: size
      TYPE(C_PTR) :: baseptr
      integer,allocatable :: arrayshape(:)
      integer :: shmrank, shmsize, shmcomm, zerokey, zerocomm, zerorank, zerosize
      integer :: fh
      integer(kind=MPI_OFFSET_KIND) :: disp       
      Character(Len=16) :: filename
!     - - - - - - - - - - - - - - - - - - - - - - - -
      call system_clock(count_rate=clock_rate)

      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                         MPI_INFO_NULL, shmcomm, mpierr)
      call MPI_COMM_rank(shmcomm, shmrank, mpierr)
      call MPI_COMM_size(shmcomm, shmsize, mpierr)
      ! if core of each node, set zerokey to be true (=1)
      if (shmrank==0) then
        zerokey=1
      else
        zerokey=0
      end if
      ! create zerocomm for master cores to communicate
      call MPI_COMM_split(MPI_COMM_WORLD,zerokey,0,zerocomm,mpierr)
      call MPI_COMM_rank(zerocomm, zerorank, mpierr)
      call MPI_COMM_size(zerocomm, zerosize, mpierr)
      if (shmrank == 0 .and. (.not. allocated(Iarr))) then
        allocate(Iarr(Ne,Nd)) 
      end if
      if (shmrank == 0) then
        do i=1,Ne
            call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, zerocomm, mpierr)
        end do   
      end if
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      allocate(arrayshape(2))
      arrayshape=(/ Ne, Nd /)
      size = 0_MPI_ADDRESS_KIND
      if (shmrank == 0) then
         size8 = Ne
         size8 = size8*Nd
         size = int(size8,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND !*8 for double ! Put the actual data size here                                                                                                                                                                
      else
         size = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, shmcomm, baseptr, win, mpierr)
      ! Obtain the location of the memory segment
      if (shmrank /= 0) then
        call MPI_Win_shared_query(win, 0, size, disp_unit, baseptr, mpierr)
      end if
      ! baseptr can now be associated with a Fortran pointer
      ! and thus used to access the shared data
      call C_F_POINTER(baseptr, Iarr_win,arrayshape)
      if (shmrank == 0) then
         Iarr_win=Iarr
      end if
      call MPI_Win_Fence(0, win, mpierr)
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)

      i8=0_int64  ! integer*8
      ih8=0_int64
      ih8H=0_int64
      call InitFormH(npes,mype) ! initialize all variables required for constructing H_IJ
      
      NumH=0
      Kherr=0
      Kgerr=0
      n0=1
      Hmin=0.d0
      t=0.d0
      if (Ksig == 2) then
        iscr=0
        xscr=0
        if (mype == 0) write(*,*) 'Screening is included'
      end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     reading/forming of the file CONF.HIJ
!     - - - - - - - - - - - - - - - - - - - - - - - - -

      ih8=NumH
      allocate(idet1(Ne),idet2(Ne),ic1(Ne),ic2(Ne),nk(Nd))
      nk=0
      ! comparison stage - - - - - - - - - - - - - - - - - - - - - - - - -
      if (mype==0) then
        Nd0=IP1+1
        n1=0
        ic1=0
        n2=Nd0-1
        do ic=1,Nc4
            n1=n1+Ndc(ic)
            if (n1 < Nd0) then
                n2=n1
                ic1=ic
            end if
        end do
        Nd0=n2
        counter1=1
        counter2=1
        counter3=1

        ! Get accumulator vectors setup:
        Call IVAccumulatorInit(iva1, vaGrowBy)
        Call IVAccumulatorInit(iva2, vaGrowBy)
        Call RVAccumulatorInit(rva1, vaGrowBy)
        
        do n=1,Nd0+1
          call Gdet(n,idet1)
          k=0
          do ic=1,Nc 
            kx=Ndc(ic)
            if (k+kx > n) kx=n-k
            if (kx /= 0) then
              call Gdet(k+1,idet2)
              call CompCD(idet1,idet2,icomp)
              if (icomp > 2) then
                k=k+kx
              else
                do k1=1,kx
                  k=k+1
                  call Gdet(k,idet2)
                  call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                  if (diff <= 2) then
                    call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                    t=Hmltn(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                    if (t /= 0) then
                      nn=n
                      kk=k
                      Call IVAccumulatorAdd(iva1, nn)
                      Call IVAccumulatorAdd(iva2, kk)
                      Call RVAccumulatorAdd(rva1, t)
                    end if
                  end if
                end do
              end if
            end if
          end do
        end do
        Call IVAccumulatorCopy(iva1, H_n0, counter1)
        Call IVAccumulatorCopy(iva2, H_k0, counter2)
        Call RVAccumulatorCopy(rva1, H_t0, counter3)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)

      if (Kl /= 1) then ! if continuing calculation and CONF.HIJ is available, skip FormH

        if (mype==0) then
          call calcMemReqs
          print*, 'chunk_size=', vaGrowBy
          print*, '===== starting FormH comparison stage ====='
        end if
  
        counter1=1
        counter2=1
        counter3=1
  
        ! Get accumulator vectors setup (or re-setup if this is rank 0):
        Call IVAccumulatorInit(iva1, vaGrowBy)
        Call IVAccumulatorInit(iva2, vaGrowBy)
        !Call RVAccumulatorInit(rva1, vaGrowBy)
  
        call system_clock(s1)
        ! each core loops through their assigned determinants and stores the number of nonzero matrix elements per det
        do n=mype+1,Nd,npes
          call Gdet_win(n,idet1)
          k=0
          do ic=1,Nc 
            kx=Ndc(ic)
            if (k+kx > n) kx=n-k
            if (kx /= 0) then
              call Gdet_win(k+1,idet2)
              call CompCD(idet1,idet2,icomp)
              if (icomp > 2) then
                k=k+kx
              else
                do k1=1,kx
                  k=k+1
                  call Gdet_win(k,idet2)
                  call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                  if (diff <= 2) then
                      nk(n)=nk(n)+1
                      nn=n
                      kk=k
                      Call IVAccumulatorAdd(iva1, nn)
                      Call IVAccumulatorAdd(iva2, kk)
                  end if
                end do
              end if
            end if
          end do
        end do

        Call IVAccumulatorCopy(iva1, H_n, counter1)
        Call IVAccumulatorCopy(iva2, H_k, counter2)
  
        Call IVAccumulatorReset(iva1)
        Call IVAccumulatorReset(iva2)
        
        call system_clock(e1)
        ttime=Real((e1-s1)/clock_rate)

        Write(counterStr,fmt='(I12)') counter1
        Call FormattedTime(ttime, memStr)
        write(*,'(2X,A,1X,I3,1X,A)'), 'core', mype, 'took '// trim(memStr)// ' with num_me = ' // trim(AdjustL(counterStr))
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  
        ih8=counter1
        call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
  
        if (mype==0) then
          call MPI_Reduce(MPI_IN_PLACE, nk(1:Nd), Nd, MPI_INTEGER, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        else
          call MPI_Reduce(nk(1:Nd), nk(1:Nd), Nd, MPI_INTEGER, MPI_SUM, 0, &
                              MPI_COMM_WORLD, mpierr)
        end if
        deallocate(nk)
  
        mem=sizeof(H_n)*4
        ! Sum all the mem sizes to get a total...
        call MPI_AllReduce(mem, memsum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        ! ...before overwriting mem with the maximum value across all workers:
        call MPI_AllReduce(mem, mem, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)
        if (mype==0) then
          Write(npesStr,fmt='(I16)') npes
          Call FormattedMemSize(memsum, memStr)
          Write(*,'(A,A,A)') 'FormH: Hamiltonian matrix requires approximately ',Trim(memStr),' of memory'
          Call FormattedMemSize(mem, memStr)
          Write(*,'(A,A,A,A,A)') 'FormH: Hamiltonian matrix requires approximately ',Trim(memStr),' of memory per core with ',Trim(AdjustL(npesStr)),' cores'
          memEstimate = memEstimate + mem
          Call FormattedMemSize(memEstimate, memStr)
          Call FormattedMemSize(memTotalPerCPU, memStr2)
          if (memTotalPerCPU /= 0) then
            if (memEstimate > memTotalPerCPU) then
              Write(*,'(A,A,A,A)'), Trim(memStr), ' required to finish conf, but only ' , Trim(memStr2) ,' is available.'
              stop
            else if (memEstimate < memTotalPerCPU) then
              Write(*,'(A,A,A,A)'), Trim(memStr), ' required to finish conf, and ' , Trim(memStr2) ,' is available.'
            end if
          else
            Write(*,'(2X,A,A,A,A)'), Trim(memStr), ' required to finish conf, but available memory was not saved to environment'
          end if
        end if
  
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        if (mype==0) print*, '===== starting FormH calculation stage ====='
        allocate(H_t(counter1))
        call system_clock(s1)
        do n=1,counter1
          nn=H_n(n)
          kk=H_k(n)
          call Gdet_win(nn,idet1)
          call Gdet_win(kk,idet2)
          call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
          call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
          t=Hmltn(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
          H_t(n)=t
        end do
  
        call system_clock(e1)
        ttime=Real((e1-s1)/clock_rate)
        Call FormattedTime(ttime, memStr)
        write(*,'(2X,A,1X,I3,1X,A)'), 'core', mype, 'took '// trim(memStr)// ' to complete calculations'
  
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Win_Fence(0, win, mpierr)
  
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Win_Free(win, mpierr)
  
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  
        ! Compute NumH, the total number of non-zero matrix elements
        call MPI_AllReduce(iscr, iscr, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_AllReduce(xscr, xscr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        deallocate(idet1, idet2, ic1, ic2)
  
        if (mype==0) then
          print*, '===== FormH calculation stage completed ====='
        end if
        ! TEST - Write CONF.HIJs in serial
        !call dclock()
        !if (Kw == 1) call Hwrite_s(mype,npes,counter1)
        !call dclock()
        !if (mype == 0) print*, 'writing CONF.HIJ in serial took ', e1-s1, 'sec'
        if (Kl /= 1 .and. Kw == 1) then
          print*, 'Writing CONF.HIJ...'
          call system_clock(s1)
          if (Kw == 1) call Hwrite(mype,npes,counter1)
          call system_clock(e1)
          if (mype == 0) print*, 'Writing CONF.HIJ in parallel took ', Real((e1-s1)/clock_rate), 'sec'
        end if
      else
        if (mype==0) then
          print*, 'Reading CONF.HIJ...'
        end if
        call system_clock(s1)
        if (Kl == 1) then
          call Hread(mype,npes,counter1)
        end if
        call system_clock(e1)
        if (mype == 0) print*, 'Reading CONF.HIJ in parallel took ', Real((e1-s1)/clock_rate), 'sec'
      end if

      ! give all cores Hmin, the minimum matrix element value
      call MPI_AllReduce(H_t(1:ih8), Hmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)

      ih8H = counter1 ! global variable for total number of matrix elements for each core

      if (mype==0) then
        write( 6,'(4X,"NumH =",I12)') NumH
        write(11,'(4X,"NumH =",I12)') NumH
!       - - - - - - - - - - - - - - - - - - - - - - - -
        if (Ksig == 2 .and. iscr > 0) then
            xscr=xscr/iscr
            write ( 6,'(5X,"For ",I12," integrals averaged screening ",F8.5)') iscr,xscr
            write (11,'(5X,"For ",I12," integrals averaged screening ",F8.5)') iscr,xscr
            if (Kherr+Kgerr > 0) then
                write ( 6,'(4X,"Extrapolation warning: small denominators.", &
                    /4X,"HintS: ",I6,"; GintS: ",I7)') Kherr,Kgerr
                write (11,'(4X,"Extrapolation warning: small denominators.", &
                    /4X,"HintS: ",I6,"; GintS: ",I7)') Kherr,Kgerr
           end if
        end if
!       - - - - - - - - - - - - - - - - - - - - - - - -
      end if

      return
    end subroutine FormH

    real(dp) function Hmltn(idet1, idet2, is, nf, i2, i1, j2, j1) 
        use determinants, Only : Rspq
        use integrals, Only : Gint, Hint
        use formj2, Only : F_J0
        implicit none
        integer, allocatable, dimension(:), intent(inout)   :: idet1, idet2
        integer, intent(inout)                              :: is, nf, i1, i2, j1, j2
        
        integer     :: iq, jq, jq0, k
        real(dp)    :: t
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        
        t=0.d0
        if (Kdsig /= 0 .and. nf <= 2) E_k=Diag(k)
        select case(nf)
          case(2) ! determinants differ by two functions
            t=t+Gint(i2,j2,i1,j1)*is 
          case(1) ! determinants differ by one function
            do iq=1,Ne
                i1=idet1(iq)
                if (i1 /= j1) t=t+Gint(j2,i1,j1,i1)*is
            end do
            t=t+Hint(j2,j1)*is
          case(0) ! determinants are equal
            do iq=1,Ne
                i1=idet1(iq)
                jq0=iq+1
                if (jq0 <= Ne) then
                    do jq=jq0,Ne
                        j1=idet1(jq)
                        t=t+Gint(i1,j1,i1,j1)*is
                    end do
                end if
                t=t+Hint(i1,i1)*is
            end do
            t=t+Gj*F_J0(idet1)
        end select
        Hmltn=t
        return
    end function Hmltn

    subroutine Diag4(mype, npes)
      ! this subroutine executes the Davidson procedure for diagonalization
      ! the subroutine Mxmpy is a computational bottleneck and was the only subroutine to be parallelized
      ! all other subroutines are performed by the master core
      use mpi
      use davidson
      implicit none
      integer  :: k1, k, i, n1, iyes, id, id2, id1, ic, id0, kskp, iter, &
                  kx, i1, it, ierr, mype, npes, mpierr
      real(dp)     :: start_time, end_time
      real(dp) :: crit, ax, x, xx, ss
      logical :: lsym
 !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        crit=1.d-6
        if (mype == 0) then
            if (Kl4 == 0) Return
            if (Nc4 > Nc) Nc4=Nc
            write (*,*) 'kl4=',Kl4,'  Nc4=',Nc4,'  Crt4=',Crt4
            ! Starting vectors of dim. < Nd0 :
            Nd0=IP1+1
            call Init4
            call Hould(Nd0,IP1,D1,E1,Z1)
            do while (Ifail /= 0)
               Write(6,'(4X,"Starting approximation of dim ",I4," failed")') Nd0
               call Init4
               call Hould(Nd0,IP1,D1,E1,Z1)
            end do
            write( 6,'(1X,"Hmin =",F14.8)') Hmin
            write(11,'(1X,"Hmin =",F14.8)') Hmin
        end if
        select case(Kv)
          case(4)
            call FormB0(mype,npes)
          case(3)
            if (mype==0) call FormB0(mype,npes)
        end select
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        if (Nd0 == Nd) return
        ! Davidson loop:
        iter = 0
        kdavidson = 0
        if (mype==0) write(*,*) 'Start with kdavidson =', kdavidson
        ax = 1
        cnx = 1
        call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        do it=1,N_it
            iter=iter+1
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            call MPI_Bcast(ax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            if (ax > crit .and. iter <= n_it) then
              if (mype == 0) write( 6,'(1X,"** Davidson iteration ",I3," for ",I2," levels **")') iter,Nlv
              if (mype == 0) write(11,'(1X,"** Davidson iteration ",I3," for ",I2," levels **")') iter,Nlv
              ! Orthonormalization of Nlv probe vectors:
              if (mype==0) then
                do i=1,Nlv
                   call Ortn(i)
                   if (Ifail /= 0) then
                      write(*,*)' Fail of orthogonalization for ',i
                      Stop
                   end if
                end do
              end if
              ! Formation of the left-upper block of the energy matrix P:
              call Mxmpy(1, mype, npes)
              call MPI_Barrier(MPI_COMM_WORLD, mpierr)
              if (mype==0) then ! start master core only block
                call FormP(1)
                lsym=K_prj == 1.OR.kskp == 1.OR.kskp == 3
                if (iter == 1 .and. lsym) then    !# averaging diagonal over confs
                  id0=1
                  do ic=1,Nc
                    id1=Ndc(ic)
                    id2=id0+id1-1
                    if (id1 > 0) then
                      ss=0.d0
                      do id=id0,id2
                        ss=ss+Diag(id)
                      end do
                      ss=ss/id1
                      Diag(id0:id2)=ss
                      id0=id0+id1
                    end if
                  end do
                  if (id2 == Nd) then
                    write(*,*) ' Diagonal averaged over rel. configurations'
                  else
                    write(*,*) ' Error: id2=',id2,' Nc=',Nc
                    stop
                  end if
                end if 
                ! Formation of Nlv additional probe vectors:
                cnx=0.d0
                do i=1,Nlv
                  i1=i+Nlv
                  call Dvdsn(i)
                  if (Iconverge(i)==0) then
                    call Ortn(i1)
                    if (Ifail /= 0) then
                      write(*,*)' Fail of orthogonalization for ',i1
                      stop
                    end if
                  end if
                end do
                if (K_prj == 1) then
                  call Prj_J(Nlv+1,Nlv,2*Nlv+1,ierr,1.d-5)
                  if (ierr /= 0) then
                    write(*,*) ' Wrong J values for Probe vectors '
                  end if
                  do i=Nlv+1,2*Nlv
                    if (Iconverge(i-Nlv)==0) then
                      call Ortn(i)
                      if (Ifail /= 0) then
                        write(*,*)' Fail of orthogonalization 2 for ',i1
                        if (kdavidson==1) then
                          stop
                        else
                          kdavidson=1
                          write(*,*) ' change kdavidson to ', kdavidson
                          Ifail=0
                          exit
                        end if
                      end if
                    end if
                  end do
                end if
              end if ! end master core only block
              call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
              call MPI_Barrier(MPI_COMM_WORLD, mpierr)
              if (cnx > Crt4) then
                ! Formation of other three blocks of the matrix P:
                call Mxmpy(2, mype, npes)
                if (mype==0) then
                  call FormP(2)
                  n1=2*Nlv
                  ! Evaluation of Nlv eigenvectors:
                  call Hould(n1,IPlv,D,E,P)
                  ax=0.d0
                  vmax=-1.d10
                  do i=1,Nlv
                    xx=0.d0
                    do k=1,Nlv
                      k1=k+Nlv
                      x=dabs(P(k1,i))
                      if (xx <= x) then
                        xx=x
                        kx=k
                      end if
                    end do
                    if (ax < xx) ax=xx
                    if (vmax < E(i)) vmax=E(i)
                    if (mype == 0) then
                      write( 6,'(1X,"E(",I2,") =",F14.8,"; admixture of vector ", & 
                             I2,": ",F10.7)') i,-(E(i)+Hmin),kx,xx
                      write(11,'(1X,"E(",I2,") =",F14.8,"; admixture of vector ", & 
                             I2,": ",F10.7)') i,-(E(i)+Hmin),kx,xx
                    end if
                  end do
                  if (kXIJ > 0) then ! write intermediate CONF.XIJ
                    if (mod(iter,kXIJ) == 0) then
                      call FormB
                    else
                      call FormBskip
                    end if
                  else
                    call FormBskip
                  end if
                end if
              else
                if (mype == 0) then
                  Write( 6,'(1X,"Davidson procedure converged")')
                  Write(11,'(1X,"Davidson procedure converged")')
                end if
                return
              end if
            end if
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        end do
        return
    end subroutine Diag4

    subroutine Print
       implicit none
       integer :: j1, nk, k, j2, n, ndk, ic, i, j, idum, ist, jmax, imax, &
                  j3
       real(dp) :: cutoff, xj, dt, del, dummy, wmx, E, D
       real(dp), allocatable, dimension(:)  :: Cc, Dd
       character(len=1), dimension(11) :: st1, st2 
       character(len=1), dimension(7)  :: stecp*7
       character(len=1), dimension(2)  :: st3*3
       character(len=1), dimension(4)  :: strsms*6
       character(len=1), dimension(3)  :: strms*3
        data st1/11*'='/,st2/11*'-'/,st3/' NO','YES'/
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Cc(Nd), Dd(Nd), W(Nc,IPlv))
        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is used for output
        if (K_is == 3) K_sms=4       !### used for output
        if (Kecp == 1) ist=7
        open(unit=16,file='CONF.XIJ', &
             status='UNKNOWN',form='unformatted')
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!      printing eigenvalues in increasing order
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        WRITE( 6,5)
        WRITE(11,5)
  5     FORMAT(4X,63('='))
        if (Ksig*Kdsig == 0) then
           WRITE( 6,15) stecp(ist),Nc,Nd,Gj
           WRITE(11,15) stecp(ist),Nc,Nd,Gj
 15        FORMAT(4X,'Energy levels (',A7,' Nc=',I6,' Nd=',I9, &
                 ');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        else
           WRITE( 6,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
           WRITE(11,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
 25        FORMAT(4X,'Energy levels ',A7,', Sigma(E =', &
                 F10.4,') extrapolation var.',I2,/4X,'(Nc=',I6, &
                 ' Nd=',I9,');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        end if
        if (C_is /= 0.d0) then
          if (K_is == 1) then
            write( *,251) C_is,Rnuc
            write(11,251) C_is,Rnuc
 251        format(4X,'Volume shift: dR_N/R_N=',F9.5,' Rnuc=',F10.7)
          else
            write( *,253) strms(K_is-1),C_is,strsms(K_sms),Klow
            write(11,253) strms(K_is-1),C_is,strsms(K_sms),Klow
 253        format(4X,A3,':',E9.2,'*(P_i dot P_k) ',A6, &
                   ' Lower component key =',I2)
          end if
        end if
        write( 6,255)
        write(11,255)
 255    format(4X,63('-'))
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        JMAX=min(Nlv,Nd)
        do J=1,JMAX
!     - - - - - - - - - - - - - - - - - - - - - - - - -
           REWIND (16)
           IMAX=J-1
           IF (IMAX < 1) GOTO 250
           DO I=1,IMAX
              READ (16)
           END DO
 250       READ (16) ER(J),xj,idum,(CC(I),I=1,ND)
           Er(j)=Er(j)+4.d0*Gj*xj*(xj+1.d0)
           E=ER(J)
           DT=E-ECORE
!     =================================================
!     === Rydberg constant is taken from "phys.par" ===
!     =================================================
           DEL=(ER(1)-ER(J))*2*DPRy
           WRITE( 6,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
           WRITE(11,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
        end do
        WRITE( 6,45)
        WRITE(11,45)
 45     FORMAT(4X,63('='))
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!       weights of configurations
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        DO J=1,JMAX
           REWIND (16)
           IMAX=J-1
           IF (IMAX < 1) GOTO 270
           DO I=1,IMAX
              READ (16)
           END DO
 270       READ (16) D,DUMMY,idum,(CC(I),I=1,ND)
           I=0
           DO IC=1,NC
              D=0.d0
              NDK=NDC(IC)
              DO K=1,NDK
                 I=I+1
                 D=D+CC(I)**2
              END DO
              W(IC,J)=D
           END DO
        END DO
        N=(JMAX-1)/5+1
        J2=0
        DO 130 K=1,N
           NK=5
           IF (K == N) NK=JMAX-(N-1)*5
           J1=J2+1
           J2=J2+NK
           J3=J2+1
           WRITE(11,65) (ST1,I=J1,J3)
 65        FORMAT(3X,66A1)
           WRITE(11,75) (I,I=J1,J2)
 75        FORMAT(15X,5(3X,I2,6X))
           WRITE(11,65) (ST2,I=J1,J3)
           WRITE(11,85) (ER(I),I=J1,J2)
 85        FORMAT(4X,'ICONF',4X,5F11.5)
           WRITE(11,65) (ST2,I=J1,J3)
           DO 140 IC=1,NC
              I=IC
              WRITE(11,95) I,(W(I,J),J=J1,J2)
 95           FORMAT(2X,I6,'      ',5F11.6)
 140       CONTINUE
           WRITE(11,65) (ST1,I=J1,J3)
 130    CONTINUE
        CLOSE(unit=16)
!     - - - - - - - - - - - - - - - - - - - - - - - - -
       deallocate(Cc,Dd,W)
 1000  RETURN
      end subroutine Print
end module conf_aux
