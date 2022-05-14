Module formj2
    ! subroutines FormJ, F_J0, F_J2, Plj, & J_av (J**2 on determinants)
    Use conf_variables

    Implicit None

  Contains

    Real(type_real) Function F_J0(idet)
        Implicit None
    
        Integer :: ia, na, ja, ma, jq0, iq, ib, nf, jq
        Integer, allocatable, dimension(:) :: idet
        Real(type_real)  :: t
    
        t=0_type_real
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

    Real(type_real) Function F_J2(idet, is, nf, ia, ic, ib, id) 
        Implicit None
        Integer, allocatable, dimension(:), intent(InOut) :: idet
        Integer, intent(InOut)                            :: is, nf, ia, ib, ic, id

        Integer :: na, ja, ma, jq0, iq, jq
        Real(type_real)  :: t

        t=0_type_real
        Select Case(nf)
            Case(0) 
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
            Case(2)
                t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
                    Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
        End Select

        F_J2=t
        Return
    End Function F_J2

    Real(type_real) Function Plj(ia, ib)
        Implicit None

        Integer :: ia, ib, na, nb, ma, mb, ja
        Real(type_real) :: t

        t=0_type_real
        na=Nh(ia)
        nb=Nh(ib)
        If (na == nb) Then
            ma=Jz(ia)
            mb=Jz(ib)
            If (ma == mb+2) Then
                ja=Jj(na)
                t=sqrt(real(ja*(ja+2)-ma*mb,kind=type_real))
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

        Integer :: k1, n1, ic1, ndn, j, ij4, ijmax, n, k, nn, kk, counter1, counter2, counter3, diff
        Real(kind=type_real) :: tt
        Integer(kind=int64) :: stot, s1, mem, maxmem, statmem
        Integer, allocatable, dimension(:) :: idet1, idet2, cntarray
        Integer :: npes, mype, mpierr, msg, maxme, endNc
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: vaGrowBy, ncGrowBy, nccnt, ncsplit
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002
        Integer :: iSign, iIndexes(3), jIndexes(3), an_id, nnc, num_done, status(MPI_STATUS_SIZE)
        Integer :: sender
        Character(Len=16) :: timeStr, memStr, memStr2, memTotStr, memTotStr2, counterStr, counterStr2, strfmt

        Call startTimer(stot)
        Allocate(idet1(Ne),idet2(Ne),cntarray(2))
        
        ! If continuing from previous calculation or J^2 matrix has already been constructed
        If (Kl == 1) Then 
            ! Read the matrix J^2 from file CONFp.JJJ
            Call ReadMatrix(Jsq%ind1,Jsq%ind2,Jsq%val,ij4,NumJ,'CONFp.JJJ',mype,npes,mpierr) 

            ! Add maximum memory per core from storing J^2 to total memory count
            Call MPI_AllReduce(ij4, ijmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            memEstimate = memEstimate + ijmax*16

        ! If continuing calculation and Hamiltonian is to be extended with more configurations
        !Else If (Kl == 3) Then
        !    ! Read the matrix J^2 from file CONFp.JJJ
        !    Call ReadMatrix(Jsq%ind1,Jsq%ind2,Jsq%val,ij4,NumJ,'CONFp.JJJ',mype,npes,mpierr) 

        !    ! Add maximum memory per core from storing J^2 to total memory count
        !    Call MPI_AllReduce(ij4, ijmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
        !    memEstimate = memEstimate + ijmax*16

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
                                tt=F_J2(idet1, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
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
                    mem = NumJ * (8_int64+type_real)
                    maxmem = maxme * (8_int64+type_real)
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
                                    tt=F_J2(idet1, iSign, diff, jIndexes(3), iIndexes(3), jIndexes(2), iIndexes(2))
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
        
            Call IVAccumulatorCopy(iva1, Jsq%ind1, counter1)
            Call IVAccumulatorCopy(iva2, Jsq%ind2, counter2)
            Call RVAccumulatorCopy(rva1, Jsq%val, counter3)
            
            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
            Call RVAccumulatorReset(rva1)
    
            Call stopTimer(s1, timeStr)
            !Write(*,'(2X,A,1X,I3,1X,A,I9)'), 'core', mype, 'took '// trim(timeStr)// ' for ij8=', counter1
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End If

        ij8=size(Jsq%val)
        ij4 = ij8

        Call MPI_AllReduce(ij8, NumJ, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)

        ! Write J^2 matrix to file CONFp.JJJ
        If (Kl /= 1) Call WriteMatrix(Jsq%ind1,Jsq%ind2,Jsq%val,ij4,NumJ,'CONFp.JJJ',mype,npes,mpierr)

        If (mype == 0) Then
            Write(counterStr,fmt='(I16)') NumJ
            strfmt = '(4X,"NumJ = ",A)'
            Write(11,strfmt) Trim(AdjustL(counterStr))
            Write( *,strfmt) Trim(AdjustL(counterStr))
            Call stopTimer(stot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> FormJ took '// trim(timeStr) // ' to complete'
        End If
        Deallocate(idet1,idet2,cntarray)

    End Subroutine FormJ
    
    Subroutine J_av(X1, nx, xj, ierr)    !# <x1|J**2|x1>
        Use mpi_f08
        Implicit None

        Integer :: ierr, i, k, n, nx, mpierr
        Real(type_real) :: r, t, xj
        Real(type_real), dimension(nx) :: X1

        ierr=0
        xj=0_type_real
        Do i=1,ij8
            n=Jsq%ind1(i)
            k=Jsq%ind2(i)
            t=Jsq%val(i)
            If (max(k,n) <= nx) Then
                r=t*X1(k)*X1(n)
                If (n /= k) r=r+r
                xj=xj+r
            Else
                Cycle
            End If
        End Do
        ! MPI Reduce sum all xj to master core here 
        Call MPI_AllReduce(MPI_IN_PLACE, xj, 1, mpi_type_real, MPI_SUM, MPI_COMM_WORLD, mpierr)
        xj=0.5d0*(sqrt(1.d0+xj)-1.d0)

        If (K_prj == 1) Then
            If (abs(xj-XJ_av) > 1.d-1) ierr=1
        End If
        Return
    End Subroutine J_av

End Module formj2