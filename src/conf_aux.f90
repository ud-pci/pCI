Module conf_aux

    Use conf_variables

    Implicit None

    Public

  Contains

    Subroutine Input
        ! This subroutine reads in parameters and configurations from CONF.INP
        Use conf_init, only : inpstr, ReadConfInp, ReadConfigurations
        Implicit None
        Integer  :: i, i1, i2, ic, nx, ny, nz, ne0, n, k
        Character(Len=1) :: name(16)
        Character(Len=32) :: strfmt
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        strfmt = '(4X,"Program Conf")'
        Write( 6,strfmt)
        Write(11,strfmt)
        ! - input from the file 'CONF.INP' - - - - - - - - - - - - - - - -
        Call ReadConfInp
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (dabs(C_is) < 1.d-6) K_is=0
        If (K_is == 0) C_is=0.d0

        Open(unit=99,file='c.in',status='OLD')
        Read (99,*) Kl, Ksig, Kdsig
        Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
        If (K_is == 2.OR.K_is == 4) Then
            Read(99,*) K_sms
            Write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): ', K_sms
            If ((K_sms-1)*(K_sms-2)*(K_sms-3) /= 0) Stop
        End If
        Close(99)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call ReadConfigurations
        ! - - - - - - -  Case kl = 2  - - - - - - - - - - -
        If (Kl == 2) Then
            Write(*,'(1X," Ksig = (0,1,2): ",I1)') Ksig 
            If (Ksig /= 0) Then
                Write(*,'(1X," Energy depEndence of Sigma (1-Yes,0-No)? ",I1)') Kdsig
            End If
            Write( 6,'(/4X,"Kl = (0-Start,1-Cont.,2-MBPT,3-Add) ",I1)') Kl
            Kecp=0
            Kl=0
        Else
            Ksig=0
        End If
        ! - - - - - - -  Case kv = 1,3  - - - - - - - - - -
        K_prj=0                    !# this key is fixed for kv=2,4
        If (Kv == 1.OR.kv == 3) Then
            K_prj=1
            Write( *,'(4X,"Selection of states with J =",F5.1)') XJ_av
            Write(11,'(4X,"Selection of states with J =",F5.1)') XJ_av
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        !If (Kl /= 1) Then ! If Kl=1, continue from a previous calculation
        !      Open(unit=16,status='UNKNOWN',file='CONF.JJJ')
        !      Close(unit=16,status='DELETE')
        !End If
        Open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        Read(16) (In(i),i=1,IPgnt)
        Read(16) (Gnt(i),i=1,IPgnt)
        Close(unit=16)
        Return
    End Subroutine Input

    Subroutine AllocateFormHArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr

        If (mype==0) Deallocate(Qnl)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(IPlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NhintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ngint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(NgintS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(Nvc)) Allocate(Nvc(Nc))
        If (.not. Allocated(Nc0)) Allocate(Nc0(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. Allocated(Jz)) Allocate(Jz(Nst))
        If (.not. Allocated(Nh)) Allocate(Nh(Nst))
        If (.not. Allocated(Diag)) Allocate(Diag(Nd))
        If (.not. Allocated(Rint1)) Allocate(Rint1(Nhint))
        If (.not. Allocated(Rint2)) Allocate(Rint2(IPbr,Ngint))
        If (.not. Allocated(Iint1)) Allocate(Iint1(Nhint))
        If (.not. Allocated(Iint2)) Allocate(Iint2(Ngint))
        If (.not. Allocated(Iint3)) Allocate(Iint3(Ngint))
        If (.not. Allocated(IntOrd)) Allocate(IntOrd(nrd))

        If (Ksig /= 0) Then
            If (.not. Allocated(Rsig)) Allocate(Rsig(NhintS))
            If (.not. Allocated(Dsig)) Allocate(Dsig(NhintS))
            If (.not. Allocated(Esig)) Allocate(Esig(NhintS)) 
            If (.not. Allocated(R_is)) Allocate(R_is(num_is))
            If (.not. Allocated(I_is)) Allocate(I_is(num_is))
            If (.not. Allocated(Rint2S)) Allocate(Rint2S(NgintS))
            If (.not. Allocated(Dint2S)) Allocate(Dint2S(NgintS))
            If (.not. Allocated(Eint2S)) Allocate(Eint2S(NgintS))
            If (.not. Allocated(Iint1S)) Allocate(Iint1S(NhintS))
            If (.not. Allocated(Iint2S)) Allocate(Iint2S(NgintS))
            If (.not. Allocated(Iint3S)) Allocate(Iint3S(NgintS))
            If (.not. Allocated(IntOrdS)) Allocate(IntOrdS(nrd))
        End If

        If (mype==0) Then
            memFormH = 0_int64
            memFormH = sizeof(Nvc)+sizeof(Nc0) &
                + sizeof(Rint1)+sizeof(Rint2)+sizeof(Iint1)+sizeof(Iint2)+sizeof(Iint3)+sizeof(Iarr)&
                + sizeof(IntOrd)
            If (Ksig /= 0) memFormH = memFormH+sizeof(Rint2S)+sizeof(Dint2S)+sizeof(Eint2S) &
                + sizeof(Iint1S)+sizeof(Iint2S)+sizeof(Iint3S) &
                + sizeof(Rsig)+sizeof(Dsig)+sizeof(Esig)+sizeof(R_is)+sizeof(I_is)+sizeof(IntOrdS)
            !Call FormattedMemSize(memFormH, memStr)
            !Write(*,'(A,A,A)') 'Allocating arrays for FormH requires ',Trim(memStr),' of memory per core' 
        End If   

        Return
    End Subroutine AllocateFormHArrays

    Subroutine DeAllocateFormHArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr
        Integer(Kind=8) :: mem 

        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    
        If (mype==0) Then
            Call FormattedMemSize(memFormH, memStr)
            Write(*,'(A,A,A)') 'De-allocating ',Trim(memStr),' of memory per core from arrays for FormH'  
        End If
    
        If (Allocated(Nvc)) Deallocate(Nvc)
        If (Allocated(Nc0)) Deallocate(Nc0)
        If (Allocated(Rint1)) Deallocate(Rint1)
        If (Allocated(Rint2)) Deallocate(Rint2)
        If (Allocated(Iint1)) Deallocate(Iint1)
        If (Allocated(Iint2)) Deallocate(Iint2)
        If (Allocated(Iint3)) Deallocate(Iint3)
        If (Allocated(Rint2S)) Deallocate(Rint2S)
        If (Allocated(Dint2S)) Deallocate(Dint2S)
        If (Allocated(Eint2S)) Deallocate(Eint2S)
        If (Allocated(Iint1S)) Deallocate(Iint1S)
        If (Allocated(Iint2S)) Deallocate(Iint2S)
        If (Allocated(Iint3S)) Deallocate(Iint3S)
        If (Allocated(Rsig)) Deallocate(Rsig)
        If (Allocated(Dsig)) Deallocate(Dsig)
        If (Allocated(Esig)) Deallocate(Esig) 
        If (Allocated(R_is)) Deallocate(R_is)
        If (Allocated(I_is)) Deallocate(I_is)
        If (Allocated(IntOrd)) Deallocate(IntOrd)
        If (Allocated(IntOrdS)) Deallocate(IntOrdS)
        If (Allocated(Iarr)) Deallocate(Iarr)
    
        !If (mype==0) Then
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
        !End If
        Return
    End Subroutine DeAllocateFormHArrays

    Subroutine AllocateDvdsnArrays(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer :: mpierr, mype, npes
        Character(Len=16) :: memStr
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        If (.not. Allocated(ArrB)) Allocate(ArrB(Nd,IPlv))
        If (.not. Allocated(Tk)) Allocate(Tk(Nlv))
        If (.not. Allocated(Tj)) Allocate(Tj(Nlv))
        If (.not. Allocated(P)) Allocate(P(IPlv,IPlv))
        If (.not. Allocated(D)) Allocate(D(IPlv))
        If (.not. Allocated(E)) Allocate(E(IPlv))
        If (.not. Allocated(Iconverge)) Allocate(Iconverge(IPlv))
        If (.not. Allocated(B1)) Allocate(B1(Nd))
        If (.not. Allocated(B2)) Allocate(B2(Nd))
        If (.not. Allocated(Z1)) Allocate(Z1(Nd0,Nd0))
        If (.not. Allocated(D1)) Allocate(D1(Nd0))
        If (.not. Allocated(E1)) Allocate(E1(Nd0))
        If (mype==0) Then
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
        End If   
        Return
    End Subroutine AllocateDvdsnArrays

    Subroutine calcMemStaticArrays
        Use str_fmt, Only : FormattedMemSize
        Implicit None

        Character(Len=16) :: memStr
    
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
        !Call FormattedMemSize(memStaticArrays, memStr)
        !Write(*,'(A,A,A)') 'calcMemReqs: Allocating static arrays requires ',Trim(memStr),' of memory per core' 
    End Subroutine calcMemStaticArrays

    Subroutine calcMemReqs
        Use str_fmt, Only : FormattedMemSize
        Implicit None
        Integer(Kind=int64) :: mem
        Integer :: bytesInteger, bytesDP, bytesReal
        Character(Len=16) :: memStr
    
        bytesInteger = 4
        bytesReal = 4
        bytesDP = 8
    
        Call calcMemStaticArrays
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
        If (memTotalPerCPU == 0) Then
            Write(*,'(A)') 'calcMemReqs: Available memory was not saved to the environment.'
        Else
            Write(*,'(A,A,A)') 'calcMemReqs: Total memory available to the job is ',Trim(memStr),' of memory per core' 
        End If
        Return
    End Subroutine calcMemReqs

    Subroutine FormH(npes, mype)
        Use mpi
        Use str_fmt, Only : FormattedMemSize, FormattedTime
        Use, intrinsic :: ISO_C_BINDING, Only : C_PTR, C_F_POINTER
        Use conf_init, Only : InitFormH ! initialization subroutines for FormH 
        Use determinants, Only : Gdet, Gdet_win, CompCD, CompD, Rspq, Rspq_phase1, Rspq_phase2
        Use matrix_io
        Use vaccumulator
        Use mpi_wins
        Implicit None

        Integer :: npes, mype, mpierr, interval, remainder, Hlim, numBins, cnt, cnt2
        Integer :: is, nf, i1, i2, j1, j2, k1, kx, n, ic, n1, n2, int, split_type, key, disp_unit, win, &
                   n0, jq, jq0, iq, i, j, icomp, k, ih4, counter1, counter2, counter3, diff, k2, totsize
        Integer :: nn, kk, msg, status(MPI_STATUS_SIZE), sender, num_done, an_id, return_msg, endnd, minme, maxme
        logical :: finished
        Integer, allocatable, dimension(:) :: idet1, idet2, mepd
        Integer, dimension(npes) :: avgs, start1, end1
        Integer, dimension(npes) :: sizes, disps
        Integer(Kind=int64)     :: start_time, end_time, stot, etot, s1, e1, s2, e2, clock_rate
        real :: ttime, ttot
        real(dp)  :: t, tt
        Integer(Kind=int64) :: ih8, i8, l8, sumd, statmem, mem, mem2, memsum, ih, cntr, avgme, numme, size8, minmem, maxmem
        Character(Len=16)     :: memStr, memStr2, memTotStr, memTotStr2, npesStr, counterStr, counterStr2
        Integer :: iSign, iIndexes(3), jIndexes(3), nnd
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: growBy, vaGrowBy, ndGrowBy, ndsplit, ndcnt, noncore, noncoreworld
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002
        Integer :: fh
        Integer(Kind=MPI_OFFSET_KIND) :: disp       
        Character(Len=16) :: filename, timeStr
!       - - - - - - - - - - - - - - - - - - - - - - - -
        Call system_clock(count_rate=clock_rate)
        If (mype==0) Call system_clock(stot)

        Call CreateIarrWindow(win, mpierr)
    
        i8=0_int64  ! Integer*8
        ih8=0_int64
        ih8H=0_int64
        Call InitFormH(npes,mype) ! initialize all variables required for constructing H_IJ
        
        NumH=0
        Kherr=0
        Kgerr=0
        n0=1
        Hmin=0.d0
        t=0.d0
        If (Ksig == 2) Then
            iscr=0
            xscr=0
            If (mype == 0) Write(*,*) 'Screening is included'
        End If
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!       Reading/forming of the file CONF.HIJ
!       - - - - - - - - - - - - - - - - - - - - - - - - - 
        ih8=NumH
        Allocate(idet1(Ne),idet2(Ne),iconf1(Ne),iconf2(Ne))

        If (Kl /= 1) Then ! If continuing calculation and CONF.HIJ is available, skip FormH
            Nd0=IP1+1
            n1=0
            n2=Nd0-1
            Do ic=1,Nc4
                n1=n1+Ndc(ic)
                If (n1 < Nd0) Then
                    n2=n1
                End If
            End Do
            Nd0=n2

            vaGrowBy = 10000000
            ndGrowBy = 100

            If (Nd0 < ndGrowBy) Then
                Continue
            Else
                ndGrowBy = Nd0+1 
            End If

            If (mype==0) Then
                Call calcMemReqs
                Write(counterStr,fmt='(I16)') vaGrowBy
                Write(counterStr2,fmt='(I16)') ndGrowBy
                Write(*,'(A)') ' vaGrowBy = '//Trim(AdjustL(counterStr))//', ndGrowBy = '//Trim(AdjustL(counterStr2))
                print*, '========== starting formation of Hamiltonian matrix =========='
            End If

            counter1=1
            counter2=1
            counter3=1
            cnt=0
            cnt2=0

            ! Get accumulator vectors setup (or re-setup If this is rank 0):
            Call IVAccumulatorInit(iva1, vaGrowBy)
            Call IVAccumulatorInit(iva2, vaGrowBy)
            Call RVAccumulatorInit(rva1, vaGrowBy)
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call system_clock(s1)

            If (mype == 0) Then        
                ! Distribute a portion of the workload of size ndGrowBy to each worker process
                Do an_id = 1, npes - 1
                   nnd = ndGrowBy*an_id + 1
                   Call MPI_SEND( nnd, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                End Do

                Do n=1,ndGrowBy
                    Call Gdet_win(n,idet1)
                    k=0
                    Do ic=1,Nc 
                        kx=Ndc(ic)
                        If (k+kx > n) kx=n-k
                        If (kx /= 0) Then
                            Call Gdet_win(k+1,idet2)
                            Call CompCD(idet1,idet2,icomp)
                            If (icomp > 2) Then
                                k=k+kx
                            Else
                                Do k1=1,kx
                                    k=k+1
                                    Call Gdet_win(k,idet2)
                                    Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                    If (diff <= 2) Then
                                        nn=n
                                        kk=k
                                        Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                        tt=Hmltn(idet1, idet2, iSign, dIff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                                        If (tt /= 0) Then
                                            cnt = cnt + 1
                                            cnt2 = cnt2 + 1
                                            Call IVAccumulatorAdd(iva1, nn)
                                            Call IVAccumulatorAdd(iva2, kk)
                                            Call RVAccumulatorAdd(rva1, tt)
                                        End If
                                    End If
                                End Do
                            End If
                        End If
                    End Do
                End Do
                
                NumH = cnt
                num_done = 0
                ndsplit = Nd/10
                ndcnt = ndsplit
                maxme = cnt2
                j=9

                Do 
                    Call MPI_RECV( cnt, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    Call MPI_RECV(cnt2, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    sender = status(MPI_SOURCE)
             
                    If (nnd + ndGrowBy <= Nd) Then
                        nnd = nnd + ndGrowBy
                        Call MPI_SEND( nnd, 1, MPI_INTEGER, &
                            sender, send_tag, MPI_COMM_WORLD, mpierr)
                    Else
                        msg = -1
                        Call MPI_SEND( msg, 1, MPI_INTEGER, &
                            sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If
            
                    NumH = NumH + cnt
                    maxme = max(cnt2,maxme)
                    mem = NumH * 16
                    maxmem = maxme * 16
                    statmem = memEstimate + mem
                    Call FormattedMemSize(statmem, memTotStr)
                    Call FormattedMemSize(memTotalPerCPU, memTotStr2)

                    If ((nnd <= ndcnt + 50 .and. nnd >= ndcnt - 50) .or. (nnd > ndcnt)) Then
                        Call system_clock(e1)
                        ttime=Real((e1-s1)/clock_rate)
                        Call FormattedTime(ttime, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Write(counterStr,fmt='(I16)') NumH
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH:', (10-j)*10, '% done in '// trim(timeStr)// ' with '//Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)//', '//trim(memStr2)//' for a single core)'
                        If (memTotalPerCPU /= 0 .and. statmem > memTotalPerCPU) Then
                            Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but only ' , Trim(memTotStr2) ,' is available.'
                            Stop
                        End If
                        j=j-1
                        ndcnt = ndcnt + ndsplit
                    End If
                    
                    If (num_done == npes-1) Then
                        Call system_clock(e1)
                        ttime=Real((e1-s1)/clock_rate)
                        Call FormattedTime(ttime, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Write(counterStr,fmt='(I16)') NumH
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH:', (10-j)*10, '% done in '// trim(timeStr)// ' with '//Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)//', '//trim(memStr2)//' for a single core)'
                        If (memTotalPerCPU /= 0) Then
                            If (statmem > memTotalPerCPU) Then
                                Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but only ' , Trim(memTotStr2) ,' is available.'
                                Stop
                            Else If (statmem < memTotalPerCPU) Then
                                Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, and ' , Trim(memTotStr2) ,' is available.'
                            End If
                        Else
                            Write(*,'(2X,A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but available memory was not saved to environment'
                        End If
                        Exit
                    End If
                End Do
            Else
                Do 
                    Call MPI_RECV ( nnd, 1 , MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    If (nnd == -1) Then
                          Exit
                    Else
                        If (Nd - nnd < ndGrowBy) Then
                            endnd = Nd
                        Else
                            endnd = nnd+ndGrowBy-1
                        End If

                        cnt=0
                        Do n=nnd,endnd
                            Call Gdet_win(n,idet1)
                            k=0
                            Do ic=1,Nc 
                                kx=Ndc(ic)
                                If (k+kx > n) kx=n-k
                                If (kx /= 0) Then
                                    Call Gdet_win(k+1,idet2)
                                    Call CompCD(idet1,idet2,icomp)
                                    If (icomp > 2) Then
                                        k=k+kx
                                    Else
                                        Do k1=1,kx
                                            k=k+1
                                            Call Gdet_win(k,idet2)
                                            Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                            If (diff <= 2) Then
                                                nn=n
                                                kk=k
                                                Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                                tt=Hmltn(idet1, idet2, iSign, dIff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                                                If (tt /= 0) Then
                                                    cnt = cnt + 1
                                                    cnt2 = cnt2 + 1
                                                    Call IVAccumulatorAdd(iva1, nn)
                                                    Call IVAccumulatorAdd(iva2, kk)
                                                    Call RVAccumulatorAdd(rva1, tt)
                                                End If
                                            End If
                                        End Do
                                    End If
                                End If
                            End Do
                        End Do
                    
                        Call MPI_SEND( cnt, 1, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                        Call MPI_SEND(cnt2, 1, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                    End if
                End do
            End If

            Call IVAccumulatorCopy(iva1, Hamil%n, counter1)
            Call IVAccumulatorCopy(iva2, Hamil%k, counter2)
            Call RVAccumulatorCopy(rva1, Hamil%t, counter3)

            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
            Call RVAccumulatorReset(rva1)
        
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            Hamil%n = PACK(Hamil%n, Hamil%t/=0)
            Hamil%k = PACK(Hamil%k, Hamil%t/=0)
            Hamil%t = PACK(Hamil%t, Hamil%t/=0)

            ih8=size(Hamil%t)
            ih4=ih8

            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)

            Call CloseIarrWindow(win, mpierr)
        
            ! Compute NumH, the total number of non-zero matrix elements
            Call MPI_AllReduce(iscr, iscr, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(xscr, xscr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Deallocate(idet1, idet2, iconf1, iconf2)
        
            If (mype==0) print*, '========== formation of Hamiltonian matrix completed =========='

            If (Kl /= 1 .and. Kw == 1)  Call WriteMatrix(Hamil,ih4,NumH,'CONF.HIJ',mype,npes,mpierr)
        Else
            If (mype==0) Then
                print*, 'Reading CONF.HIJ...'
            End If
            Call system_clock(s1)
            If (Kl == 1) Then
                Call ReadMatrix(Hamil,ih4,NumH,'CONF.HIJ',mype,npes,mpierr)
            End If
            Call system_clock(e1)
            If (mype == 0) print*, 'Reading CONF.HIJ in parallel took ', Real((e1-s1)/clock_rate), 'sec'
        End If
    
        ! give all cores Hmin, the minimum matrix element value
        Call MPI_AllReduce(Hamil%t(1:ih8), Hmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
    
        ih8H = ih8 ! global variable for total number of matrix elements for each core
    
        If (mype==0) Then
            Write( 6,'(4X,"NumH =",I12)') NumH
            Write(11,'(4X,"NumH =",I12)') NumH
            If (Ksig == 2 .and. iscr > 0) Then
                xscr=xscr/iscr
                Write ( 6,'(5X,"For ",I12," integrals averaged screening ",F8.5)') iscr,xscr
                Write (11,'(5X,"For ",I12," integrals averaged screening ",F8.5)') iscr,xscr
                If (Kherr+Kgerr > 0) Then
                    Write ( 6,'(4X,"Extrapolation warning: small denominators.", &
                        /4X,"HintS: ",I6,"; GintS: ",I7)') Kherr,Kgerr
                    Write (11,'(4X,"Extrapolation warning: small denominators.", &
                        /4X,"HintS: ",I6,"; GintS: ",I7)') Kherr,Kgerr
               End If
            End If
        End If
    
        If (mype == 0) Then
            Call system_clock(etot)
            ttot=Real((etot-stot)/clock_rate)
            Call FormattedTime(ttot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> FormH took '// trim(timeStr) // ' to complete'
        End If

        Return
    End Subroutine FormH

    real(dp) function Hmltn(idet1, idet2, is, nf, i2, i1, j2, j1) 
        Use determinants, Only : Rspq
        Use integrals, Only : Gint, Hint
        Use formj2, Only : F_J0, F_J2
        Implicit None
        Integer, allocatable, dimension(:), intent(inout)   :: idet1, idet2
        Integer, intent(inout)                              :: is, nf, i1, i2, j1, j2
        
        Integer     :: iq, jq, jq0, k
        real(dp)    :: t
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        
        t=0.d0
        If (Kdsig /= 0 .and. nf <= 2) E_k=Diag(k)
        Select Case(nf)
            Case(2) ! determinants dIffer by two functions
                t=t+Gint(i2,j2,i1,j1)*is 
                t=t+Gj*F_J2(idet1,idet2)
            Case(1) ! determinants dIffer by one function
                Do iq=1,Ne
                    i1=idet1(iq)
                    If (i1 /= j1) t=t+Gint(j2,i1,j1,i1)*is
                End Do
                t=t+Hint(j2,j1)*is
            Case(0) ! determinants are equal
                Do iq=1,Ne
                    i1=idet1(iq)
                    jq0=iq+1
                    If (jq0 <= Ne) Then
                        Do jq=jq0,Ne
                            j1=idet1(jq)
                            t=t+Gint(i1,j1,i1,j1)*is
                        End Do
                    End If
                    t=t+Hint(i1,i1)*is
                End Do
                t=t+Gj*F_J2(idet1,idet2)
        End Select
        Hmltn=t
        Return
    End function Hmltn

    Subroutine Diag4(mype, npes)
        ! this Subroutine executes the Davidson procedure for diagonalization
        ! the Subroutine Mxmpy is a computational bottleneck and was the only Subroutine to be parallelized
        ! all other Subroutines are performed by the master core
        Use mpi
        Use str_fmt, Only : FormattedTime
        Use davidson
        Implicit None
    
        Integer  :: k1, k, i, j, n1, iyes, id, id2, id1, ic, id0, kskp, iter, &
                    kx, i1, i2, it, mype, npes, mpierr
        Integer(Kind=int64) :: stot, etot, clock_rate
        Real :: ttot
        real(dp)     :: start_time, End_time
        real(dp) :: crit, ax, x, xx, ss
        logical :: lsym
        Character(Len=16) :: timeStr
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Call system_clock(count_rate=clock_rate)
        If (mype==0) Call system_clock(stot)
        crit=1.d-6
        If (mype == 0) Then
            If (Kl4 == 0) Return
            If (Nc4 > Nc) Nc4=Nc
            Write (*,*) 'kl4=',Kl4,'  Nc4=',Nc4,'  Crt4=',Crt4
            ! Starting vectors of dim. < Nd0 :
            Nd0=IP1+1
            Call Init4
            Call Hould(Nd0,IP1,D1,E1,Z1)
            Do while (Ifail /= 0)
               Write(6,'(4X,"Starting approximation of dim ",I4," failed")') Nd0
               Call Init4
               Call Hould(Nd0,IP1,D1,E1,Z1)
            End Do
            Write( 6,'(1X,"Hmin =",F14.8)') Hmin
            Write(11,'(1X,"Hmin =",F14.8)') Hmin
        End If
        Call FormB0(mype,npes)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        If (Nd0 == Nd) Return
        ! Davidson loop:
        iter = 0
        kdavidson = 0
        If (mype==0) Write(*,*) 'Start with kdavidson =', kdavidson
        ax = 1
        cnx = 1
        Call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Do it=1,N_it
            iter=iter+1
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(ax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            If (ax > crit .and. iter <= n_it) Then
                If (mype == 0) Write( 6,'(1X,"** Davidson iteration ",I3," for ",I2," levels **")') iter,Nlv
                If (mype == 0) Write(11,'(1X,"** Davidson iteration ",I3," for ",I2," levels **")') iter,Nlv
                ! Orthonormalization of Nlv probe vectors:
                If (mype==0) Then
                  Do i=1,Nlv
                     Call Ortn(i)
                     If (Ifail /= 0) Then
                        Write(*,*)' Fail of orthogonalization for ',i
                        Stop
                     End If
                  End Do
                End If
                ! Formation of the left-upper block of the energy matrix P:
                Call Mxmpy(1, mype, npes)
                Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
                If (mype==0) Then ! start master core only block
                  Call FormP(1)
                  lsym=K_prj == 1.OR.kskp == 1.OR.kskp == 3
                  If (iter == 1 .and. lsym) Then    !# averaging diagonal over confs
                    id0=1
                    Do ic=1,Nc
                      id1=Ndc(ic)
                      id2=id0+id1-1
                      If (id1 > 0) Then
                        ss=0.d0
                        Do id=id0,id2
                          ss=ss+Diag(id)
                        End Do
                        ss=ss/id1
                        Diag(id0:id2)=ss
                        id0=id0+id1
                      End If
                    End Do
                    If (id2 == Nd) Then
                      Write(*,*) ' Diagonal averaged over rel. configurations'
                    Else
                      Write(*,*) ' Error: id2=',id2,' Nc=',Nc
                      Stop
                    End If
                  End If 
                  ! Formation of Nlv additional probe vectors:
                  cnx=0.d0
                  Do i=1,Nlv
                    i1=i+Nlv
                    Call Dvdsn(i)
                    If (Iconverge(i)==0) Then
                      Call Ortn(i1)
                      If (Ifail /= 0) Then
                        Write(*,*)' Fail of orthogonalization for ',i1
                        Stop
                      End If
                    End If
                  End Do
                End If

                If (K_prj == 1) Then
                    Call Prj_J(Nlv+1,Nlv,2*Nlv+1,1.d-5,mype,npes)
                    If (mype == 0) Then
                        Do i=Nlv+1,2*Nlv
                          If (Iconverge(i-Nlv)==0) Then
                            Call Ortn(i)
                            If (Ifail /= 0) Then
                              If (mype == 0) Write(*,*)' Fail of orthogonalization 2 for ',i1
                              If (kdavidson==1) Then
                                Stop
                              Else
                                kdavidson=1
                                If (mype == 0) Write(*,*) ' change kdavidson to ', kdavidson
                                Ifail=0
                                exit
                              End If
                            End If
                          End If
                        End Do
                    End If
                End If
                Call MPI_Bcast(cnx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
                If (cnx > Crt4) Then
                    ! Formation of other three blocks of the matrix P:
                    Call Mxmpy(2, mype, npes)
                    If (mype==0) Then
                        Call FormP(2)
                        n1=2*Nlv
                        ! Evaluation of Nlv eigenvectors:
                        Call Hould(n1,IPlv,D,E,P)
                        ax=0.d0
                        vmax=-1.d10
                        Do i=1,Nlv
                            xx=0.d0
                            Do k=1,Nlv
                                k1=k+Nlv
                                x=dabs(P(k1,i))
                                If (xx <= x) Then
                                    xx=x
                                    kx=k
                                End If
                            End Do
                            If (ax < xx) ax=xx
                            If (vmax < E(i)) vmax=E(i)
                            If (mype == 0) Then
                                Write( 6,'(1X,"E(",I2,") =",F14.8,"; admixture of vector ", & 
                                        I2,": ",F10.7)') i,-(E(i)+Hmin),kx,xx
                                Write(11,'(1X,"E(",I2,") =",F14.8,"; admixture of vector ", & 
                                        I2,": ",F10.7)') i,-(E(i)+Hmin),kx,xx
                            End If
                        End Do
                        If (kXIJ > 0) Then ! Write intermediate CONF.XIJ
                            If (mod(iter,kXIJ) == 0) Then
                                Call FormB
                            Else
                                Call FormBskip
                            End If
                        Else
                            Call FormBskip
                        End If
                    End If
                Else
                    If (mype == 0) Then
                        Write( 6,'(1X,"Davidson procedure converged")')
                        Write(11,'(1X,"Davidson procedure converged")')
                    End If
                    Return
                End If
            End If
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End Do

        If (mype == 0) Then
            Call system_clock(etot)
            ttot=Real((etot-stot)/clock_rate)
            Call FormattedTime(ttot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> Davidson procedure took '// trim(timeStr) // ' to complete'
        End If

        deallocate(Hamil%n, Hamil%k, Hamil%t)
        Return
    End Subroutine Diag4

    Subroutine WriteFinalXIJ(mype,npes)
        Use mpi
        Use formj2, Only : J_av
        Use davidson, Only : Prj_J
        Implicit None
        Integer :: i, n, ierr, mype, npes, mpierr

        Call MPI_Bcast(Tk(1:Nlv), Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        D1(1:Nlv)=Tk(1:Nlv)
        If (K_prj == 1) Then
            Call Prj_J(1,Nlv,Nlv+1,1.d-8,mype,npes)
        End If
        If (mype == 0) Then
            open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        End If
    
        Do n=1,Nlv
            Call MPI_Bcast(ArrB(1:Nd,n), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            Call J_av(ArrB(1,n),Nd,Tj(n),ierr,mype,npes)  ! calculates expectation values for J^2
            If (mype==0) write(17) D1(n),Tj(n),Nd,(ArrB(i,n),i=1,Nd)
        End Do

        If (mype==0) close(unit=17)

        Return
    End Subroutine WriteFinalXIJ

    Subroutine PrintResults
        Implicit None

        Integer :: j1, nk, k, j2, n, ndk, ic, i, j, idum, ist, jmax, imax, &
                   j3
        real(dp) :: cutoff, xj, dt, del, dummy, wmx, E, D
        real(dp), allocatable, dimension(:)  :: Cc, Dd
        Character(Len=1), dimension(11) :: st1, st2 
        Character(Len=1), dimension(10)  :: stecp*7
        Character(Len=1), dimension(2)  :: st3*3
        Character(Len=1), dimension(4)  :: strsms*6
        Character(Len=1), dimension(3)  :: strms*3
        data st1/11*'='/,st2/11*'-'/,st3/' NO','YES'/
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Cc(Nd), Dd(Nd), W(Nc,IPlv))
        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is Used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7
        Open(unit=16,file='CONF.XIJ',status='UNKNOWN',form='unformatted')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! printing eigenvalues in increasing order
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        WRITE( 6,5)
        WRITE(11,5)
  5     FORMAT(4X,63('='))
        If (Ksig*Kdsig == 0) Then
           WRITE( 6,15) stecp(ist),Nc,Nd,Gj
           WRITE(11,15) stecp(ist),Nc,Nd,Gj
 15        FORMAT(4X,'Energy levels (',A7,' Nc=',I6,' Nd=',I9, &
                 ');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        Else
           WRITE( 6,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
           WRITE(11,25) stecp(ist),E_0,Kexn,Nc,Nd,Gj
 25        FORMAT(4X,'Energy levels ',A7,', Sigma(E =', &
                 F10.4,') extrapolation var.',I2,/4X,'(Nc=',I6, &
                 ' Nd=',I9,');  Gj =',F7.4,/4X,'N',6X,'JTOT',12X, &
                 'EV',16X,'ET',9X,'DEL(CM**-1)')
        End If
        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                Write( *,251) C_is,Rnuc
                Write(11,251) C_is,Rnuc
 251            format(4X,'Volume shIft: dR_N/R_N=',F9.5,' Rnuc=',F10.7)
            Else
                Write( *,253) strms(K_is-1),C_is,strsms(K_sms),Klow
                Write(11,253) strms(K_is-1),C_is,strsms(K_sms),Klow
 253            format(4X,A3,':',E9.2,'*(P_i Dot P_k) ',A6, &
                     ' Lower component key =',I2)
            End If
        End If
        Write( 6,255)
        Write(11,255)
 255    format(4X,63('-'))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        JMAX=min(Nlv,Nd)
        Do J=1,JMAX
            REWIND (16)
            IMAX=J-1
            IF (IMAX < 1) GOTO 250
            DO I=1,IMAX
                READ (16)
            END DO
 250        READ (16) ER(J),xj,idum,(CC(I),I=1,ND)
            Er(j)=Er(j)+4.d0*Gj*xj*(xj+1.d0)
            E=ER(J)
            DT=E-ECORE
            ! Rydberg constant is taken from "phys.par"
            DEL=(ER(1)-ER(J))*2*DPRy
            WRITE( 6,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
            WRITE(11,'(3X,I2,F14.9,F14.8,F18.6,F15.2)') j,xj,E,DT,DEL
        End Do
        WRITE( 6,45)
        WRITE(11,45)
 45     FORMAT(4X,63('='))
        ! weights of configurations
        DO J=1,JMAX
            REWIND (16)
            IMAX=J-1
            IF (IMAX < 1) GOTO 270
            DO I=1,IMAX
                READ (16)
            END DO
 270        READ (16) D,DUMMY,idum,(CC(I),I=1,ND)
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
 65         FORMAT(3X,66A1)
            WRITE(11,75) (I,I=J1,J2)
 75         FORMAT(15X,5(3X,I2,6X))
            WRITE(11,65) (ST2,I=J1,J3)
            WRITE(11,85) (ER(I),I=J1,J2)
 85         FORMAT(4X,'ICONF',4X,5F11.5)
            WRITE(11,65) (ST2,I=J1,J3)
            DO 140 IC=1,NC
                I=IC
                WRITE(11,95) I,(W(I,J),J=J1,J2)
 95             FORMAT(2X,I6,'      ',5F11.6)
 140        CONTINUE
            WRITE(11,65) (ST1,I=J1,J3)
 130    CONTINUE
        CLOSE(unit=16)
        close(unit=6)
        close(unit=11)
        RETURN
    End Subroutine PrintResults

End Module conf_aux
