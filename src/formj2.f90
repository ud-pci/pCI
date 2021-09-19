Module formj2
    ! subroutines FormJ, F_J0, F_J2, Plj, & J_av (J**2 on determinants)
    Use conf_variables

    Implicit None

  Contains

    Real(dp) Function F_J0(idet)
        Implicit None
    
        Integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is
        Integer, allocatable, dimension(:) :: idet
        Real(dp)  :: t
    
        t=0.d0
        If (nf == 0) Then !determinants are equal
            t=mj*mj
            Do iq=1,Ne
                ia=idet(iq)
                na=Nh(ia)
                ja=Jj(na)
                ma=Jz(ia)
                t=t+ja*(ja+2)-ma**2
                jq0=iq+1
                If (jq0 <= Ne) Then
                    Do jq=jq0,Ne
                        ib=idet(jq)
                        t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
                    End Do
                End If
            End Do
        End If
        F_J0=t
        Return
    End Function F_J0

    Real(dp) Function F_J2(idet1, idet2, is, nf, ia, ic, ib, id) 
        Implicit None
        Integer, allocatable, dimension(:), intent(InOut) :: idet1, idet2
        Integer, intent(InOut)                            :: is, nf, ia, ib, ic, id

        Integer :: na, ja, ma, jq0, iq, jq
        Real(dp)  :: t

        t=0.d0
        Select Case(nf)
            Case(0) 
                t=mj*mj
                Do iq=1,Ne
                    ia=idet1(iq)
                    na=Nh(ia)
                    ja=Jj(na)
                    ma=Jz(ia)
                    t=t+ja*(ja+2)-ma**2
                    jq0=iq+1
                    If (jq0 <= Ne) Then
                        Do jq=jq0,Ne
                            ib=idet1(jq)
                            t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
                        End Do
                    End If
                End Do
            Case(2)
                t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
                    Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
        End Select

        F_J2=t
        Return
    End Function F_J2

    Real(dp) Function Plj(ia, ib)
        Implicit None

        Integer :: ia, ib, na, nb, ma, mb, ja
        Real(dp) :: t

        t=0.d0
        na=Nh(ia)
        nb=Nh(ib)
        If (na == nb) Then
            ma=Jz(ia)
            mb=Jz(ib)
            If (ma == mb+2) Then
                ja=Jj(na)
                t=dsqrt(ja*(ja+2)-ma*mb+0.d0)
            End If
        End If
        Plj=t
        Return
    End Function Plj

    Subroutine FormJ(mype, npes)
        Use mpi
        Use str_fmt, Only : startTimer, stopTimer, FormattedMemSize, FormattedTime
        Use vaccumulator
        Use determinants, Only : Gdet, Gdet, Rspq, Rspq_phase1, Rspq_phase2
        Use matrix_io
        Implicit None

        Integer :: k1, n1, ic1, ndn, i, j, ij4, ijmax, n, k, nn, kk, counter1, counter2, counter3, diff
        Real :: ttime, ttot
        Real(dp) :: t, tt
        Integer(kind=int64) :: stot, etot, s1, e1, clock_rate, jstart, jend, memsum, mem, maxmem, statmem
        Integer, allocatable, dimension(:) :: idet1, idet2, nk, cntarray
        Integer :: npes, mype, mpierr, msg, maxme, endNc
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: vaGrowBy, ncGrowBy, nccnt, ncsplit
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002
        Integer :: iSign, iIndexes(3), jIndexes(3), an_id, nnc, num_done, return_msg, status(MPI_STATUS_SIZE)
        Integer :: fh, sender
        Integer(kind=MPI_OFFSET_KIND) :: disp 
        Character(Len=16) :: timeStr, memStr, memStr2, memTotStr, memTotStr2, counterStr, counterStr2, strfmt

        Call startTimer(stot)
        Allocate(idet1(Ne),idet2(Ne),cntarray(2))
        
        ! If continuing from previous calculation or J^2 matrix has already been constructed
        If (Kl == 1) Then 
            ! Read the matrix J^2 from file CONFp.JJJ
            Call ReadMatrix(Jsq,ij4,NumJ,'CONFp.JJJ',mype,npes,mpierr) 

            ! Add maximum memory per core from storing J^2 to total memory count
            Call MPI_AllReduce(ij4, ijmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            memEstimate = memEstimate + ijmax*16

        ! If starting a new computation and J^2 matrix has not been constructed
        Else
            vaGrowBy = 1000000
            ncGrowBy = 1
            
            If (mype == 0) Then
                Open (unit=18,file='CONFp.JJJ',status='UNKNOWN',form='UNFORMATTED')
                Close(unit=18,status='DELETE')
                Write(counterStr,fmt='(I16)') vaGrowBy
                Write(counterStr2,fmt='(I16)') ncGrowBy
                Write(*,'(A)') ' vaGrowBy = '//Trim(AdjustL(counterStr))//', ncGrowBy = '//Trim(AdjustL(counterStr2))
                print*, '========== Starting calculation stage of FormJ =========='
            End If
    
            counter1=1
            counter2=1
            counter3=1
            cntarray=0
    
            Call IVAccumulatorInit(iva1, vaGrowBy)
            Call IVAccumulatorInit(iva2, vaGrowBy)
            Call RVAccumulatorInit(rva1, vaGrowBy)
    
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)

            Call startTimer(s1)
    
            If (mype == 0) Then
                ! Distribute a portion of the workload of size ncGrowBy to each worker process
                Do an_id = 1, npes - 1
                   nnc = ncGrowBy*an_id + 1
                   Call MPI_SEND( nnc, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                End Do

                Do ic1=1, ncGrowBy
                    ndn=Ndc(ic1)
                    n=sum(Ndc(1:ic1-1))
                    Do n1=1,ndn
                        n=n+1
                        ndr=n
                        Call Gdet(n,idet1)
                        k=n-n1
                        Do k1=1,n1
                            k=k+1
                            Call Gdet(k,idet2)
                            Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                            If (diff == 0 .or. diff == 2) Then
                                nn=n
                                kk=k
                                Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                tt=F_J2(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                                If (tt /= 0) Then
                                    cntarray = cntarray + 1
                                    Call IVAccumulatorAdd(iva1, nn)
                                    Call IVAccumulatorAdd(iva2, kk)
                                    Call RVAccumulatorAdd(rva1, tt)
                                End If
                            End If
                        End Do
                    End Do
                End Do

                NumJ = cntarray(1)
                num_done = 0
                ncsplit = Nc/10
                nccnt = ncsplit
                maxme = cntarray(2)
                j=9

                Do 
                    Call MPI_RECV(cntarray, 2, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    sender = status(MPI_SOURCE)
             
                    If (nnc + ncGrowBy <= Nc) Then
                        nnc = nnc + ncGrowBy
                        Call MPI_SEND( nnc, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                    Else
                        msg = -1
                        Call MPI_SEND( msg, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If
                    
                    NumJ = NumJ + cntarray(1)
                    maxme = max(cntarray(2),maxme)
                    mem = NumJ * 16_int64
                    maxmem = maxme * 16_int64
                    statmem = memEstimate + maxmem
                    Call FormattedMemSize(statmem, memTotStr)
                    Call FormattedMemSize(memTotalPerCPU, memTotStr2)

                    If (nnc == nccnt .and. nnc /= ncsplit*10) Then
                        Call stopTimer(s1, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        Write(counterStr,fmt='(I16)') NumJ
                        Write(*,'(2X,A,1X,I3,A)'), 'FormJ:', (10-j)*10, '% done in '// trim(timeStr)// &
                                ' with '//Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)// &
                                ', '//trim(memStr2)//' for a single core)'
                        If (memTotalPerCPU /= 0 .and. memEstimate > memTotalPerCPU) Then
                            Write(*,'(A,A,A,A)'), 'At least '// Trim(memTotStr), ' is required to finish conf, but only ' , &
                            Trim(memTotStr2) ,' is available.'
                            Stop
                        End If
                        j=j-1
                        nccnt = nccnt + ncsplit
                    End If
                    
                    If (num_done == npes-1) Then
                        Call stopTimer(s1, timeStr)
                        Call FormattedMemSize(mem, memStr)
                        Call FormattedMemSize(maxmem, memStr2)
                        memEstimate = memEstimate + maxmem
                        Write(counterStr,fmt='(I16)') NumJ
                        Write(*,'(2X,A,1X,I3,A)'), 'FormJ:', (10-j)*10, '% done in '// trim(timeStr)// ' with '// &
                        Trim(AdjustL(counterStr)) // ' elements (Mem='// trim(memStr)//', '//trim(memStr2)//' for a single core)'
                        Exit
                    End If
                End Do
            Else
                Do 
                    Call MPI_RECV ( nnc, 1 , MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)

                    If (nnc == -1) Then
                        Exit
                    Else
                        If (Nc - nnc < ncGrowBy) Then
                            endnc = Nc
                        Else
                            endnc = nnc+ncGrowBy-1
                        End If

                        cntarray(1)=0
                        Do ic1=nnc,endnc
                            ndn=Ndc(ic1)
                            n=sum(Ndc(1:ic1-1))
                            Do n1=1,ndn
                              n=n+1
                              ndr=n
                              Call Gdet(n,idet1)
                              k=n-n1
                              Do k1=1,n1
                                k=k+1
                                Call Gdet(k,idet2)
                                Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                If (diff == 0 .or. diff == 2) Then
                                    nn=n
                                    kk=k
                                    Call Rspq_phase2(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                                    tt=F_J2(idet1, idet2, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
                                    If (tt /= 0) Then
                                        cntarray = cntarray + 1
                                        Call IVAccumulatorAdd(iva1, nn)
                                        Call IVAccumulatorAdd(iva2, kk)
                                        Call RVAccumulatorAdd(rva1, tt)
                                    End If
                                End If
                              End Do
                            End Do
                        End Do
                    
                        Call MPI_SEND( cntarray, 2, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                    End If
                End Do
            End If
        
            Call IVAccumulatorCopy(iva1, Jsq%n, counter1)
            Call IVAccumulatorCopy(iva2, Jsq%k, counter2)
            Call RVAccumulatorCopy(rva1, Jsq%t, counter3)
            
            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
            Call RVAccumulatorReset(rva1)
    
            Call stopTimer(s1, timeStr)
            !Write(*,'(2X,A,1X,I3,1X,A,I9)'), 'core', mype, 'took '// trim(timeStr)// ' for ij8=', counter1
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End If

        ij8=size(Jsq%t)
        ij4 = ij8

        Call MPI_AllReduce(ij8, NumJ, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)

        ! Write J^2 matrix to file CONFp.JJJ
        If (Kl /= 1) Call WriteMatrix(Jsq,ij4,NumJ,'CONFp.JJJ',mype,npes,mpierr)

        If (mype == 0) Then
            Write(counterStr,fmt='(I16)') NumJ
            strfmt = '(4X,"NumJ = ",A)'
            Write(11,strfmt) Trim(AdjustL(counterStr))
            Write( *,strfmt) Trim(AdjustL(counterStr))
            Call stopTimer(stot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> FormJ took '// trim(timeStr) // ' to complete'
        End If
        Deallocate(idet1,idet2,Jz,Nh,cntarray)

    End Subroutine FormJ
    
    Subroutine J_av(X1, nx, xj, ierr, mype, npes)    !# <x1|J**2|x1>
        Use mpi
        Use determinants, Only : Gdet
        Implicit None

        Integer :: ierr, i, k, n, nx, mpierr
        Integer, optional :: mype, npes
        Integer*8 :: mi, nj
        Real(dp) :: r, t, xj, xj2
        Real(dp), dimension(nx) :: X1
        Integer, allocatable, dimension(:) :: idet1, idet2

        ierr=0
        xj=0.d0
        Do i=1,ij8
            n=Jsq%n(i)
            k=Jsq%k(i)
            t=Jsq%t(i)
            If (max(k,n) <= nx) Then
                r=t*X1(k)*X1(n)
                If (n /= k) r=r+r
                xj=xj+r
            Else
                Exit
            End If
        End Do
        ! MPI Reduce sum all xj to master core here 
        If (Present(mype)) Call MPI_AllReduce(xj, xj, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        xj=0.5d0*(dsqrt(1.d0+xj)-1.d0)
        If (K_prj == 1) Then
            If (dabs(xj-XJ_av) > 1.d-1) ierr=1
        End If
        Return
    End Subroutine J_av

End Module formj2