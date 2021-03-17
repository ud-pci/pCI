! subroutines FormJ, F_J2, Plj, & J_av (J**2 on determinants)
module formj2
    use conf_variables
    implicit none

    contains
    real(dp) function F_J0(idet)
      implicit none
      integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is
      integer, allocatable, dimension(:) :: idet
      real(dp)  :: t
      t=0.d0
      if (nf == 0) then !determinants are equal
        t=mj*mj
        do iq=1,Ne
          ia=idet(iq)
          na=Nh(ia)
          ja=Jj(na)
          ma=Jz(ia)
          t=t+ja*(ja+2)-ma**2
          jq0=iq+1
          if (jq0 <= Ne) then
            do jq=jq0,Ne
              ib=idet(jq)
              t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
            end do
          end if
        end do
      end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
      F_J0=t
      return
    end function F_J0

    real(dp) function F_J2(idet1,idet2)
      use determinants, only : Rspq
      implicit none
      integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is
      integer, allocatable, dimension(:) :: idet1, idet2
      real(dp)  :: t
      t=0.d0
      call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
      if (nf == 0) then !determinants are equal
        t=mj*mj
        do iq=1,Ne
          ia=idet1(iq)
          na=Nh(ia)
          ja=Jj(na)
          ma=Jz(ia)
          t=t+ja*(ja+2)-ma**2
          jq0=iq+1
          if (jq0 <= Ne) then
            do jq=jq0,Ne
              ib=idet1(jq)
              t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
            end do
          end if
        end do
      end if
      if (nf == 2) then !determinants differ by two functions
         t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
           Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
      end if
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      F_J2=t
      return
    end function F_J2

    real(dp) function F_J2_new(idet1,idet2, is, nf, i2, i1, j2, j1)
      use determinants, only : Rspq
      implicit none
      integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is, i1, i2, j1, j2
      integer, allocatable, dimension(:) :: idet1, idet2
      real(dp)  :: t
      t=0.d0
      if (nf == 0) then !determinants are equal
        t=mj*mj
        do iq=1,Ne
          ia=idet1(iq)
          na=Nh(ia)
          ja=Jj(na)
          ma=Jz(ia)
          t=t+ja*(ja+2)-ma**2
          jq0=iq+1
          if (jq0 <= Ne) then
            do jq=jq0,Ne
              ib=idet1(jq)
              t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
            end do
          end if
        end do
      end if
      if (nf == 2) then !determinants differ by two functions
         t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
           Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
      end if
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      F_J2_new=t
      return
    end function F_J2_new

    real(dp) function Plj(ia,ib)
      implicit none
      integer :: ia, ib, na, nb, ma, mb, ja
      real(dp) :: t
      t=0.d0
      na=Nh(ia)
      nb=Nh(ib)
      if (na == nb) then
         ma=Jz(ia)
         mb=Jz(ib)
         if (ma == mb+2) then
            ja=Jj(na)
            t=dsqrt(ja*(ja+2)-ma*mb+0.d0)
         end if
      end if
      Plj=t
      return
    end function Plj

    recursive integer function triStart(mype,npes)
      implicit none
      integer :: mype, npes
      real :: t
      if (mype == 0) then
        t = 1.0
      else if (mype == npes) then
        t = Nc+1.0
      else
        t = sqrt((1.0/(npes-1.0))*(Nc**2-1.0)+triStart(mype-1,npes)**2)
      end if
      triStart = int(t)
      return
    end function triStart

    subroutine FormJ(mype, npes)
      use mpi
      Use str_fmt, Only : FormattedTime
      use, intrinsic :: ISO_C_BINDING, only : C_PTR, C_F_POINTER
      use vaccumulator
      use determinants, only : Gdet, Gdet_win, Rspq, Rspq_phase1, Rspq_phase2
       implicit none
       integer :: k1, n1, ic1, ndn, numjj, i, ij4, n, k, nn, kk, io, counter1, counter2, counter3, diff
       real :: ttime
       real(dp) :: t
       integer(kind=int64) :: size8, ij8, s1, e1, clock_rate, jstart, jend, memsum
       integer(kind=int64), dimension(npes) :: ij8s
       integer, allocatable, dimension(:) :: idet1, idet2, nk
       integer, allocatable, dimension(:) :: J_n0, J_k0
       real(dp), allocatable, dimension(:) :: J_t0
       integer :: npes, mype, mpierr, interval, remainder, startNc, endNc, sizeNc, counter
       integer :: split_type, key, disp_unit, win
       Type(IVAccumulator)   :: iva1, iva2
       Type(RVAccumulator)   :: rva1
       integer, Parameter    :: vaGrowBy = 1000000
       integer :: iSign, iIndexes(3), jIndexes(3)
       integer(kind=MPI_ADDRESS_KIND) :: size
       TYPE(C_PTR) :: baseptr
       integer, allocatable :: arrayshape(:)
       integer :: shmrank, shmsize, shmcomm, zerokey, zerocomm, zerorank, zerosize
       integer :: fh
       integer(kind=MPI_OFFSET_KIND) :: disp 
       Character(Len=16) :: filename, timeStr, memStr
       ! - - - - - - - - - - - - - - - - - - - - - - - - -
        call system_clock(count_rate=clock_rate)
        allocate(idet1(Ne),idet2(Ne),nk(Nc))
        
        ij8=0_int64
        NumJ=0_int64
        if (Kl == 1) then ! if continuing from previous calculation, skip forming CONF.JJJ
          open (unit=18,file='CONF.JJJ',status='OLD',form='UNFORMATTED',iostat=io)
          if (io /= 0) write( 6,'(4X,"Forming matrix J**2")')
          do
            read (18,iostat=io) numjj,k,n,t
            if (io > 0) then 
              print*, 'I/O ERROR at FormJ'
            else if (io < 0) then
              write(11,'(4X,"NumJ =",I12)') NumJ
              write( *,'(4X,"NumJ =",I12)') NumJ
              close(unit=18)
              exit
            else
              NumJ=numjj
            end if
          end do
        else ! else start new calculation
          if (mype == 0) then
            write( 6,'(4X,"Forming matrix J**2")')
            open (unit=18,file='CONF.JJJ',status='UNKNOWN',form='UNFORMATTED')
            close(unit=18,status='DELETE')
            open(unit=18,file='CONF.JJJ',status='NEW',form='UNFORMATTED')
          end if

          if (Kv==4) then
            counter1=1
            counter2=1
            counter3=1

            Call IVAccumulatorInit(iva1, vaGrowBy)
            Call IVAccumulatorInit(iva2, vaGrowBy)

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
  
            call system_clock(s1)

            interval = Nc/npes
            remainder = mod(Nc,npes)

            startNc = mype*interval+1
            endNc = (mype+1)*interval
            if (mype == npes-1) then
              endNc=Nc
            end if
            !do ic1=startNc, endNc
            do ic1=mype+1,Nc,npes
              ndn=Ndc(ic1)
              n=sum(Ndc(1:ic1-1))
              !print*,mype,n
              do n1=1,ndn
                n=n+1
                ndr=n
                call Gdet_win(n,idet1)
                k=n-n1
                do k1=1,n1
                  k=k+1
                  call Gdet_win(k,idet2)
                  call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                  if (diff == 0 .or. diff == 2) then
                    nn=n
                    kk=k
                    Call IVAccumulatorAdd(iva1, nn)
                    Call IVAccumulatorAdd(iva2, kk)
                  end if
                end do
              end do
            end do

            Call IVAccumulatorCopy(iva1, J_n, counter1)
            Call IVAccumulatorCopy(iva2, J_k, counter2)
    
            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
  
            allocate(J_t(counter1))
            do n=1,counter1
              nn=J_n(n)
              kk=J_k(n)
              call Gdet_win(nn,idet1)
              call Gdet_win(kk,idet2)
              t=F_J2(idet1, idet2)
              J_t(n)=t
            end do

            call system_clock(e1)
            ttime=Real((e1-s1)/clock_rate)
            Call FormattedTime(ttime, timeStr)
            write(*,'(2X,A,1X,I3,1X,A,I9)'), 'core', mype, 'took '// trim(timeStr)// ' for ij8=', counter1
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            ij8s=0_int64
            ij8=counter1
            call MPI_AllReduce(ij8, NumJ, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
            call MPI_AllReduce(ij8, maxJcore, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            call MPI_AllGather(ij8, 1, MPI_INTEGER8, ij8s, 1, MPI_INTEGER8, MPI_COMM_WORLD, mpierr)
            ij8J = ij8
            !allocate(J_n(NumJ),J_k(NumJ),J_t(NumJ))
            !J_n=0
            !J_k=0
            !J_t=0.d0
            !if (mype == 0) then
            !  jstart = 1
            !  jend = ij8
            !else if (mype == npes-1) then
            !  jstart = sum(ij8s(1:mype)) + 1
            !  jend = NumJ
            !else
            !  jstart = sum(ij8s(1:mype)) + 1
            !  jend = sum(ij8s(1:mype+1))
            !end if
  !
            !J_n(jstart:jend)=J_n0(1:ij8)
            !J_k(jstart:jend)=J_k0(1:ij8)
            !J_t(jstart:jend)=J_t0(1:ij8)
            !deallocate(J_n0,J_k0,J_t0)
  
            !ij4=NumJ
            !call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            !if (mype==0) then
            !  call MPI_Reduce(MPI_IN_PLACE, J_n(1:ij4), ij4, MPI_INTEGER, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !  call MPI_Reduce(MPI_IN_PLACE, J_k(1:ij4), ij4, MPI_INTEGER, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !  call MPI_Reduce(MPI_IN_PLACE, J_t(1:ij4), ij4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !else
            !  call MPI_Reduce(J_n(1:ij4), J_n(1:ij4), ij4, MPI_INTEGER, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !  call MPI_Reduce(J_k(1:ij4), J_k(1:ij4), ij4, MPI_INTEGER, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !  call MPI_Reduce(J_t(1:ij4), J_t(1:ij4), ij4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            !                    MPI_COMM_WORLD, mpierr)
            !end if

            !if (mype == 0) then
            !  call system_clock(s1)
            !  do i=1,NumJ
            !    write(18) i,J_k(i),J_n(i),J_t(i)
            !  end do
            !end if

            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            call MPI_Win_Fence(0, win, mpierr)
      
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            call MPI_Win_Free(win, mpierr)
      
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)

          else if (Kv==3) then
            if (mype == 0) then
              call system_clock(s1)
              n=0
              do ic1=1,Nc
                ndn=Ndc(ic1)
                do n1=1,ndn
                  n=n+1
                  ndr=n
                  call Gdet(n,idet1)
                  k=n-n1
                  do k1=1,n1
                    k=k+1
                    call Gdet(k,idet2)
                    t=F_J2(idet1,idet2)
                    if (t /= 0.d0) then
                      NumJ=NumJ+1
                      write(18) NumJ,k,n,t
                    end if
                  end do
                end do
              end do
            end if
          end if    

          if (mype == 0) then  
            call system_clock(e1)
            ttime=Real((e1-s1)/clock_rate)
             Call FormattedTime(ttime, timeStr) 
            write(*,'(2X,A)'), 'TIMING >>> Writing CONF.JJJ took '// trim(timeStr) // ' to complete'
            write(11,'(4X,"NumJ =",I12)') NumJ
            write( *,'(4X,"NumJ =",I12)') NumJ
            close(unit=18)
          end if
        end if
        deallocate(idet1,idet2)
       Return
    end subroutine FormJ
    
    subroutine J_av(X1,nx,xj,ierr,mype,npes)    !# <x1|J**2|x1>
      use mpi
      use determinants, only : Gdet
      implicit none
      integer :: ierr, i, k, n, nx, mpierr
      integer, optional :: mype, npes
      integer*8 :: mi, nj
      real(dp) :: r, t, xj, xj2
      real(dp), dimension(nx) :: X1
      integer, allocatable, dimension(:) :: idet1, idet2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(idet1(Ne),idet2(Ne))
      ierr=0
      xj=0.d0
      select case(Kv)
        case(4)
          nj=ij8J
        case(3)
          open(unit=18,file='CONF.JJJ',status='OLD',form='unformatted')
          nj=NumJ
      end select

      do i=1,nj
        select case(Kv)
        case(4)
          n=J_n(i)
          k=J_k(i)
          t=J_t(i)
        case(3)
          read(18) mi,k,n,t
        end select
        if (max(k,n) <= nx) then
          r=t*X1(k)*X1(n)
          if (n /= k) r=r+r
          xj=xj+r
        else
          exit
        end if
      end do
      ! MPI Reduce sum all xj to master core here 
      if (present(mype)) call MPI_AllReduce(xj, xj, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
      xj=0.5d0*(dsqrt(1.d0+xj)-1.d0)
      if (Kv == 3) close(unit=18)
      if(K_prj == 1) then
        if(dabs(xj-XJ_av) > 1.d-1) ierr=1
      end if
      deallocate(idet1,idet2)
      return
    end subroutine J_av
end module formj2