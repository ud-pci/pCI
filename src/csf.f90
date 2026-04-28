Module csf

    Use conf_variables
    Use determinants, Only : Rspq
    Use formj2, Only : Plj
    ! Use davidson, Only : Hould

    Implicit None

    Private

    Public :: jbasis_init, jbasis, nonequiv_conf, formh_sym, unsym, reorder_det, idif

    Real(dp), Allocatable :: ccj_module(:)

Contains

    Subroutine jbasis(nconf,ncsf,nccj,max_ndcs,mype,npes)
        Use mpi_f08
        Use str_fmt, Only : startTimer, stopTimer
        Implicit None
        External :: DSYEV

        Integer :: mype, npes, ierr, lwork
        Integer :: nconf, ncsf, nccj, max_ndcs
        Integer :: iconf_neq, iconf, ndi, ncsfi, ic1, n1, n2, n, k, nf, is, ia, ib, ic, id, iq, na, ja, ma, jq0, jq, i1, i2, j, ifail, jtt
        Integer :: mesplit, iconf_local_count, min_rank, n0_tasks
        Integer, Allocatable, Dimension(:) :: idet1, idet2, first_iconf_for_neq
        Integer, Allocatable :: sort_order(:), task_rank(:)
        Integer(Kind=int64), Allocatable :: task_load(:)

        Real(dp) :: t, tj
        Real(dp), Allocatable :: zz(:,:), de(:), dd(:), zz_buf(:,:)
        Real(dp), Dimension(4) :: realtmp

        Type :: t_eigvec_buf
            Real(dp), Allocatable :: data(:)
        End Type
        Type(t_eigvec_buf), Allocatable :: eigvec_bufs(:)

        Character(Len=256) :: strfmt
        Integer(Kind=int64) :: start_time
        Character(Len=16) :: timeStr

        If (.not. allocated(mdcs)) Allocate(mdcs(nconf))
        If (.not. allocated(ndcs)) Allocate(ndcs(nconf_neq))
        If (.not. allocated(ndc_neq)) Allocate(ndc_neq(nconf_neq))
        If (.not. allocated(iplace_cj)) Allocate(iplace_cj(nconf_neq))
        If (.not. allocated(idet1)) Allocate(idet1(Ne))
        If (.not. allocated(idet2)) Allocate(idet2(Ne))
        If (.not. Allocated(nc_neq)) Allocate(nc_neq(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. allocated(Mdc)) Allocate(Mdc(Nc))
        If (.not. allocated(idt)) allocate(idt(Nd,Ne))
        If (.not. allocated(Nh)) Allocate(Nh(Nst))
        If (.not. allocated(Jz)) Allocate(Jz(Nst))

        ndcs = 0

        ! Pre-compute ndc_neq and the first representative iconf for each iconf_neq.
        ! All ranks have ndc/nc_neq from jbasis_init, so no MPI needed.
        ndc_neq = 0
        Allocate(first_iconf_for_neq(nconf_neq))
        first_iconf_for_neq = 0
        Do iconf = 1, nconf
            iconf_neq = nc_neq(iconf)
            If (first_iconf_for_neq(iconf_neq) == 0) Then
                first_iconf_for_neq(iconf_neq) = iconf
                ndc_neq(iconf_neq) = ndc(iconf)
            End If
        End Do

        If (mype == 0) Then
            Call startTimer(start_time)
            strfmt = '(/4x,"Calculating matrix J**2")'
            Write(*,strfmt)
            Write(11,strfmt)
            strfmt = '(6x,"Iconf_neq",3x,"Iconf",4x,"Ndi",5x,"Ncsfi",5x,"Rank")'
            Write(*,strfmt)
        End If

        strfmt = '(6x,i5,3x,i7,4x,i5,3(4x,i5))'
        ! Load-balance assignment: sort by descending ndi**3, then greedily assign to least-loaded rank.
        ! All ranks run this identically, so no communication is needed.
        Allocate(sort_order(nconf_neq), task_rank(nconf_neq), task_load(0:npes-1))
        task_load = 0_int64
        Do i1 = 1, nconf_neq
            sort_order(i1) = i1
        End Do
        Do i1 = 2, nconf_neq
            j = sort_order(i1)
            ndi = ndc_neq(j)
            i2 = i1 - 1
            Do While (i2 >= 1 .and. ndc_neq(sort_order(i2)) < ndi)
                sort_order(i2+1) = sort_order(i2)
                i2 = i2 - 1
            End Do
            sort_order(i2+1) = j
        End Do
        Do i1 = 1, nconf_neq
            iconf_neq = sort_order(i1)
            min_rank = 0
            Do i2 = 1, npes-1
                If (task_load(i2) < task_load(min_rank)) min_rank = i2
            End Do
            task_rank(iconf_neq) = min_rank
            task_load(min_rank) = task_load(min_rank) + int(ndc_neq(iconf_neq), int64)**3
        End Do
        n0_tasks = count(task_rank == 0)
        mesplit = max(1, (n0_tasks + 9) / 10)
        iconf_local_count = 0
        Deallocate(sort_order, task_load)

        ! All ranks compute their assigned tasks in parallel; no inter-rank communication in the loop
        Allocate(eigvec_bufs(nconf_neq))
        Do iconf_neq = 1, nconf_neq
            If (task_rank(iconf_neq) /= mype) Cycle
            iconf = first_iconf_for_neq(iconf_neq)
            ndi   = ndc_neq(iconf_neq)
            Allocate(zz(ndi,ndi), de(ndi))
            zz = 0.d0
            n1 = mdc(iconf)+1
            n2 = n1+ndi-1
            Call build_j2(n1, n2, ndi, idet1, idet2, zz)
            If (ndi > 0) Then
                Call DSYEV('V','U',ndi,zz,ndi,de,realtmp,-1,ifail)
                lwork = Nint(realtmp(1))
                Allocate(dd(lwork))
                Call DSYEV('V','U',ndi,zz,ndi,de,dd,lwork,ifail)
                If (ifail /= 0) Write(*,*) mype, ' rank: dsyev failed with ifail =', ifail
                Deallocate(dd)
            End If
            ncsfi = 0
            Allocate(zz_buf(ndi,ndi))
            Do i1 = 1, ndi
                tj  = 0.5d0*(dsqrt(1.d0+de(i1))-1.d0)
                jtt = 2*tj+0.0001
                If (dabs(2*tj-jtt) > 1.d-7) Then
                    Write(*,'(/2x,a/2x,a,f16.8,/2x,2(2x,a,i3))') &
                        '*** Value of j in jbasis is wrong ***','J=',tj,'Configuration:',iconf,'Rank',mype
                    Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                End If
                If (jtt == mj) Then
                    ncsfi = ncsfi+1
                    zz_buf(1:ndi,ncsfi) = zz(1:ndi,i1)
                End If
            End Do
            ndcs(iconf_neq) = ncsfi
            If (ncsfi > 0) Then
                Allocate(eigvec_bufs(iconf_neq)%data(ndi * ncsfi))
                eigvec_bufs(iconf_neq)%data = reshape(zz_buf(1:ndi,1:ncsfi), [ndi*ncsfi])
            End If
            Deallocate(zz_buf, de, zz)
            Write(*,strfmt) iconf_neq, iconf, ndi, ncsfi, mype
            If (mype == 0) Then
                iconf_local_count = iconf_local_count + 1
                If (mod(iconf_local_count, mesplit) == 0 .and. iconf_local_count < n0_tasks) Then
                    Call stopTimer(start_time, timeStr)
                    Write(*,'(2X,A,1X,I3,A)') 'J**2 calculation:', &
                        iconf_local_count * 100 / n0_tasks, '% done in '//trim(timeStr)
                End If
            End If
        End Do

        Call MPI_Barrier(MPI_COMM_WORLD, ierr)
        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            Write(*,'(2X,A,1X,I3,A)') 'J**2 calculation:', 100, '% done in '//trim(timeStr)
        End If

        ! Share ndcs across all ranks (each rank only set its own tasks)
        Call MPI_Allreduce(MPI_IN_PLACE, ndcs, nconf_neq, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

        nccj = 0
        Do iconf_neq = 1, nconf_neq
            nccj = nccj + ndcs(iconf_neq) * ndc_neq(iconf_neq)
        End Do

        iplace_cj(1) = 0
        Do iconf_neq = 2, nconf_neq
            ic1 = iconf_neq-1
            iplace_cj(iconf_neq) = iplace_cj(ic1) + ndcs(ic1)*ndc_neq(ic1)
        End Do

        ! Each rank fills its own portion of ccj_module; Allreduce(SUM) assembles on all ranks
        If (nccj > 0) Then
            Allocate(ccj_module(nccj))
            ccj_module = 0.d0
            Do iconf_neq = 1, nconf_neq
                If (task_rank(iconf_neq) /= mype) Cycle
                ncsfi = ndcs(iconf_neq)
                If (ncsfi > 0) Then
                    ndi = ndc_neq(iconf_neq)
                    i1 = iplace_cj(iconf_neq) + 1
                    i2 = i1 + ndi*ncsfi - 1
                    ccj_module(i1:i2) = eigvec_bufs(iconf_neq)%data
                End If
            End Do
            Call MPI_Allreduce(MPI_IN_PLACE, ccj_module, nccj, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        End If

        Do iconf_neq = 1, nconf_neq
            If (Allocated(eigvec_bufs(iconf_neq)%data)) Deallocate(eigvec_bufs(iconf_neq)%data)
        End Do
        Deallocate(eigvec_bufs, task_rank)

        ncsf = 0
        max_ndcs = 0
        mdcs(1) = 0
        Do iconf = 1, nconf
            iconf_neq = nc_neq(iconf)
            ncsfi = ndcs(iconf_neq)
            ncsf = ncsf + ncsfi
            If (ncsfi > max_ndcs) max_ndcs = ncsfi
            If (iconf > 1) Then
                ic1 = nc_neq(iconf-1)
                mdcs(iconf) = mdcs(iconf-1) + ndcs(ic1)
            End If
        End Do

        If (mype == 0) Then
            Write(*,'(4x,a,i7,2x,a,i8)') 'Nconf=',nconf,'Ncsf=',ncsf
            Write(11,'(4x,a,i7,2x,a,i8)') 'Nconf=',nconf,'Ncsf=',ncsf
            If (ncsf == 0) Then
                Write(*,'(/4x,a,f4.1,a)') 'Term ',jm,' is absent in all configurations'
                Write(11,'(/4x,a,f4.1,a)') 'Term ',jm,' is absent in all configurations'
                Stop
            End If
        End If

        Deallocate(first_iconf_for_neq, idet1, idet2)

        If (mype == 0) Then
            Call stopTimer(start_time, timeStr)
            Write(*,'(2X,A)') 'TIMING >>> CSF basis construction took '//trim(timeStr)//' to complete'
            Write(11,'(2X,A)') 'TIMING >>> CSF basis construction took '//trim(timeStr)//' to complete'
        End If

    End Subroutine jbasis

    Subroutine build_j2(n1, n2, ndi, idet1, idet2, zz)
        Implicit None
        Integer, Intent(In)    :: n1, n2, ndi
        Integer, Allocatable, Intent(InOut) :: idet1(:), idet2(:)
        Real(dp), Intent(InOut) :: zz(ndi,ndi)

        Integer :: n, k, nf, is, ia, ib, ic, id, iq, na, ja, ma, jq0, jq, i1, i2
        Real(dp) :: t

        Do n = n1, n2
            idet1(1:Ne) = idt(n,1:Ne)
            Do k = n1, n
                idet2(1:Ne) = idt(k,1:Ne)
                t = 0.d0
                If (k /= n) Then
                    Call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
                    If (nf /= 2) Cycle
                    t = plj(ia,ic)*plj(id,ib)+plj(ic,ia)*plj(ib,id)-plj(ia,id)*plj(ic,ib)-plj(id,ia)*plj(ib,ic)
                    t = t*is
                    If (t == 0.d0) Cycle
                Else
                    t = mj*mj
                    Do iq = 1, Ne
                        ia = idet1(iq)
                        na = nh(ia)
                        ja = jj(na)
                        ma = jz(ia)
                        t  = t + ja*(ja+2) - ma**2
                        jq0 = iq+1
                        If (jq0 > Ne) Cycle
                        Do jq = jq0, Ne
                            ib = idet1(jq)
                            t  = t - plj(ia,ib)**2 - plj(ib,ia)**2
                        End Do
                    End Do
                    If (t == 0.d0) Cycle
                End If
                i1 = n-n1+1
                i2 = k-n1+1
                zz(i1,i2) = t
                zz(i2,i1) = t
            End Do
        End Do
    End Subroutine build_j2

    Subroutine nonequiv_conf(nconf)
        Implicit None
        Integer :: nconf

        Integer :: iconf, iconf_neq, jconf, njj, nj, nii, ni, nsj, nj0, nsi, ni0, nii1, ni1, i
        Integer, Allocatable, Dimension(:) :: ind_nj
        Character(Len=256) :: strfmt

        write( *,'(/4x,a)') 'Creating list of nonequivalent  configurations...'

        Allocate(ni_conf(nconf), nf_conf(nconf), nc_neq(nconf), ind_nj(Ns))

        do iconf=1,nconf
          ni_conf(iconf)=nc0(iconf)+1
          nf_conf(iconf)=nc0(iconf)+Nvc(iconf)
        end do

        iconf_neq=0

        conf1: do iconf=1,nconf
            conf2: do jconf=1,iconf-1
                do njj=ni_conf(jconf),nf_conf(jconf)
                    nj=nip(njj)
                    ind_nj(nj)=0
                enddo
                inner: do nii=ni_conf(iconf),nf_conf(iconf)
                    ni=nip(nii)
                    Do njj=ni_conf(jconf),nf_conf(jconf)
                        nj=nip(njj)
                        if (ind_nj(nj).eq.1) Cycle
                        if (ll(ni).ne.ll(nj)) Cycle
                        if (jj(ni).ne.jj(nj)) Cycle
                        if (nq(nii).ne.nq(njj)) Cycle
                        ind_nj(nj)=1
                        Cycle inner
                    End Do
                    Cycle conf2
                End Do inner
                do nii=ni_conf(iconf),nf_conf(iconf)
                    ni=nip(nii)
                    ind_nj(ni)=0
                end do
                nsj=nf_conf(jconf)-ni_conf(jconf)+1
                do nj0=1,nsj
                    njj=nj0+ni_conf(jconf)-1
                    nj=nip(njj)
                    nsi=nf_conf(iconf)-ni_conf(iconf)+1
                    Do ni0=1,nsi
                        nii=ni0+ni_conf(iconf)-1
                        ni=nip(nii)
                        If (ind_nj(ni).eq.1) Cycle
                        If (ll(ni).ne.ll(nj)) Cycle
                        If (jj(ni).ne.jj(nj)) Cycle
                        If (nq(nii).ne.nq(njj)) Cycle
                        If (nj0.ne.ni0) Then
                            nii1=nj0+ni_conf(iconf)-1
                            ni1=nip(nii1)
                            i=nip(nii1)
                            nip(nii1)=nip(nii)
                            nip(nii)=i
                            i=nq(nii1)
                            nq(nii1)=nq(nii)
                            nq(nii)=i
                            i=ind_nj(ni1)
                            ind_nj(ni1)=ind_nj(ni)
                            ind_nj(ni)=i
                        End If
                    End Do
                End Do
                nc_neq(iconf)=nc_neq(jconf)
                Cycle conf1
            end do conf2
            iconf_neq=iconf_neq+1
            nc_neq(iconf)=iconf_neq
        end do conf1

        Deallocate(ind_nj)

        nconf_neq=iconf_neq
        strfmt = '(4x,"Number of nonequivalent configurations:",i7,/)'
        write( *,strfmt) nconf_neq
        write(11,strfmt) nconf_neq
    End Subroutine nonequiv_conf

    Subroutine formh_sym(nconf,ncsf,nccj,max_ndcs,mype,npes)
        Use mpi_f08
        Use vaccumulator
        Use determinants, Only : calcNd0
        Use str_fmt, Only : startTimer, stopTimer
        Implicit None

        Integer :: mype, npes, mpierr
        Integer :: nconf, ncsf, nccj, max_ndcs
        Integer :: numzero, n0, iconf, iconf_neq, nci, ndi, jconf, jconf_neq, ncj, ndj, idf
        Integer :: n1, n2, n, id, k1, k2, k, jd, jc, jc2, ic, iccj, jccj, counter1, counter2, counter3
        Integer :: j, mesplit, iconf_local_count
        Integer :: an_id, nnd, num_done, sender, iconf_task
        Integer :: cntarray(2)
        Integer(Kind=int64) :: ih8_max
        Type(MPI_STATUS) :: status
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: vaGrowBy, ndGrowBy, ndsplit, ndcnt
        Integer(Kind=int64)   :: s1, ih8_before
        Real(dp) :: Hmin, hij
        Real(dp), Allocatable, Dimension(:) :: ccj, buf
        Real(dp), Allocatable, Dimension(:,:) :: zzc
        Character(Len=256) :: strfmt
        Character(Len=16)     :: memStr, memStr2, memStr3, memStr4, memStr5, memTotStr, memTotStr2, counterStr, counterStr2, timeStr
        Integer, Parameter    :: send_tag = 2001, return_tag = 2002

        If (.not. allocated(ccj)) allocate(ccj(nccj))
        If (.not. allocated(zzc)) allocate(zzc(max_ndcs,max_ndcs))
        If (.not. allocated(buf)) allocate(buf(max_ndcs))
        If (.not. allocated(idet1)) allocate(idet1(Ne))
        If (.not. allocated(idet2)) allocate(idet2(Ne))

        vaGrowBy = vaBinSize
        ndGrowBy = 1

        If (mype==0) Then
            Write(counterStr,fmt='(I16)') vaGrowBy
            Write(counterStr2,fmt='(I16)') ndGrowBy
            Write(*,'(A)') ' vaGrowBy = '//Trim(AdjustL(counterStr))//', ndGrowBy = '//Trim(AdjustL(counterStr2))
            print*, '========== Starting calculation stage of FormH =========='
        End If

        Call IVAccumulatorInit(iva1, vaGrowBy)
        Call IVAccumulatorInit(iva2, vaGrowBy)
        Call RVAccumulatorInit(rva1, vaGrowBy)

        ! Reading the list of the symmetrized coefficients
        ccj(1:nccj) = ccj_module(1:nccj)

        mesplit = max(1, nconf / 10)
        j = 1
        iconf_local_count = 0
        Call startTimer(s1)

        NumH=0
        Kherr=0
        Kgerr=0
        numzero=0
        n0=1
        Hmin=0.d0
        Hamil%minval = 0.d0
        if (Ksig.EQ.2) then
           iscr=0
           xscr=0
           If (mype == 0) write(*,*) 'Screening is included'
        end if

        ! calculation of the matrix elements
        If (npes == 1) Then
            do iconf=1,nconf
                iconf_local_count = iconf_local_count + 1
                If (mod(iconf_local_count, mesplit)==0 .and. j < 10) Then
                    Call stopTimer(s1, timeStr)
                    Write(*,'(2X,A,1X,I3,A)'), 'FormH_sym calculation stage:', j*10, '% done in '// trim(timeStr)
                    j = j + 1
                End If
                iconf_neq=nc_neq(iconf)
                nci=ndcs(iconf_neq)
                if (nci.eq.0) cycle
                ndi=ndc(iconf)
                do jconf=1,iconf
                    jconf_neq=nc_neq(jconf)
                    ncj=ndcs(jconf_neq)
                    if (ncj.eq.0) cycle
                    ndj=ndc(jconf)
                    idf=idif(iconf,jconf)
                    if (idf.gt.2) cycle
                    zzc(1:nci,1:ncj)=0.d0
                    n1=mdc(iconf)+1
                    n2=n1+ndc(iconf)-1
                    do n=n1,n2
                        id=n-n1+1
                        idet1(1:Ne)=idt(n,1:Ne)
                        k1=mdc(jconf)+1
                        k2=k1+ndc(jconf)-1
                        buf(1:ncj)=0.d0
                        do k=k1,k2
                            idet2(1:Ne)=idt(k,1:Ne)
                            if (Kdsig.NE.0) E_k=Diag(k)
                            call Hmatrix(idf,idet1,idet2,hij)
                            if (dabs(hij).lt.1.d-20) cycle
                            jd=k-k1+1
                            do jc=1,ncj
                                jccj=jd+(jc-1)*ndj+iplace_cj(jconf_neq)
                                buf(jc)=buf(jc)+hij*ccj(jccj)
                            end do
                        end do
                        do ic=1,nci
                            iccj=id+(ic-1)*ndi+iplace_cj(iconf_neq)
                            do jc=1,ncj
                                zzc(ic,jc)=zzc(ic,jc)+buf(jc)*ccj(iccj)
                            end do
                        end do
                    end do
                    do ic=1,nci
                        jc2=ncj
                        if (iconf.eq.jconf) jc2=ic
                        do jc=1,jc2
                            hij=zzc(ic,jc)
                            n=ic+mdcs(iconf)
                            k=jc+mdcs(jconf)
                            if (hij.NE.0.d0) then
                                ih8=ih8+1
                                if (n.EQ.1.AND.k.EQ.1) Hmin=hij
                                if (n.EQ.k.AND.hij.LT.Hmin) Hmin=hij
                                Call IVAccumulatorAdd(iva1, n)
                                Call IVAccumulatorAdd(iva2, k)
                                Call RVAccumulatorAdd(rva1, hij)
                                Hamil%minval = Hmin
                            else
                                numzero=numzero+1
                            end if
                        end do
                    end do
                end do
            end do
            NumH=ih8
            Call stopTimer(s1, timeStr)
            Write(*,'(2X,A,1X,I3,A)') 'FormH_sym calculation stage:', 100, '% done in '//trim(timeStr)
        Else
            If (mype == 0) Then
                ! Master: distribute iconf rows dynamically to workers
                nnd = 1
                num_done = 0
                Do an_id = 1, npes - 1
                    If (nnd <= nconf) Then
                        Call MPI_SEND(nnd, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                        nnd = nnd + 1
                    Else
                        Call MPI_SEND(-1, 1, MPI_INTEGER, an_id, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If
                End Do

                Do While (num_done < npes - 1)
                    Call MPI_RECV(cntarray, 1, MPI_INTEGER, MPI_ANY_SOURCE, return_tag, MPI_COMM_WORLD, status, mpierr)
                    sender = status%MPI_SOURCE

                    iconf_local_count = iconf_local_count + 1
                    If (mod(iconf_local_count, mesplit) == 0 .and. j < 10) Then
                        Call stopTimer(s1, timeStr)
                        Write(*,'(2X,A,1X,I3,A)'), 'FormH_sym calculation stage:', j*10, '% done in '// trim(timeStr)
                        j = j + 1
                    End If

                    If (nnd <= nconf) Then
                        Call MPI_SEND(nnd, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        nnd = nnd + 1
                    Else
                        Call MPI_SEND(-1, 1, MPI_INTEGER, sender, send_tag, MPI_COMM_WORLD, mpierr)
                        num_done = num_done + 1
                    End If
                End Do
                Call stopTimer(s1, timeStr)
                Write(*,'(2X,A,1X,I3,A)') 'FormH_sym calculation stage:', 100, '% done in '//trim(timeStr)
            Else
                ! Workers: receive iconf rows and compute matrix elements
                Do
                    Call MPI_RECV(iconf_task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
                    If (iconf_task == -1) Exit

                    iconf = iconf_task
                    ih8_before = ih8
                    iconf_neq = nc_neq(iconf)
                    nci = ndcs(iconf_neq)
                    If (nci /= 0) Then
                        ndi = ndc(iconf)
                        do jconf=1,iconf
                            jconf_neq=nc_neq(jconf)
                            ncj=ndcs(jconf_neq)
                            if (ncj.eq.0) cycle
                            ndj=ndc(jconf)
                            idf=idif(iconf,jconf)
                            if (idf.gt.2) cycle
                            zzc(1:nci,1:ncj)=0.d0
                            n1=mdc(iconf)+1
                            n2=n1+ndc(iconf)-1
                            do n=n1,n2
                                id=n-n1+1
                                idet1(1:Ne)=idt(n,1:Ne)
                                k1=mdc(jconf)+1
                                k2=k1+ndc(jconf)-1
                                buf(1:ncj)=0.d0
                                do k=k1,k2
                                    idet2(1:Ne)=idt(k,1:Ne)
                                    if (Kdsig.NE.0) E_k=Diag(k)
                                    call Hmatrix(idf,idet1,idet2,hij)
                                    if (dabs(hij).lt.1.d-20) cycle
                                    jd=k-k1+1
                                    do jc=1,ncj
                                        jccj=jd+(jc-1)*ndj+iplace_cj(jconf_neq)
                                        buf(jc)=buf(jc)+hij*ccj(jccj)
                                    end do
                                end do
                                do ic=1,nci
                                    iccj=id+(ic-1)*ndi+iplace_cj(iconf_neq)
                                    do jc=1,ncj
                                        zzc(ic,jc)=zzc(ic,jc)+buf(jc)*ccj(iccj)
                                    end do
                                end do
                            end do
                            do ic=1,nci
                                jc2=ncj
                                if (iconf.eq.jconf) jc2=ic
                                do jc=1,jc2
                                    hij=zzc(ic,jc)
                                    n=ic+mdcs(iconf)
                                    k=jc+mdcs(jconf)
                                    if (hij.NE.0.d0) then
                                        ih8=ih8+1
                                        if (n.EQ.1.AND.k.EQ.1) Hmin=hij
                                        if (n.EQ.k.AND.hij.LT.Hmin) Hmin=hij
                                        Call IVAccumulatorAdd(iva1, n)
                                        Call IVAccumulatorAdd(iva2, k)
                                        Call RVAccumulatorAdd(rva1, hij)
                                        Hamil%minval = Hmin
                                    else
                                        numzero=numzero+1
                                    end if
                                end do
                            end do
                        end do
                    End If

                    cntarray(1) = int(ih8 - ih8_before)
                    Call MPI_SEND(cntarray, 1, MPI_INTEGER, 0, return_tag, MPI_COMM_WORLD, mpierr)
                End Do
            End If
        End If
        
        Call MPI_AllReduce(ih8, NumH, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(ih8, ih8_max, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mpierr)

        Call IVAccumulatorCopy(iva1, Hamil%ind1, counter1)
        Call IVAccumulatorCopy(iva2, Hamil%ind2, counter2)
        Call RVAccumulatorCopy(rva1, Hamil%val, counter3)

        Call IVAccumulatorReset(iva1)
        Call IVAccumulatorReset(iva2)
        Call RVAccumulatorReset(rva1)

        memEstimate = memEstimate + ih8_max * 16

        Call MPI_AllReduce(MPI_IN_PLACE, numzero, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(Hamil%minval, Hamil%minval, 1, mpi_type_real, MPI_MIN, MPI_COMM_WORLD, mpierr)

        If (mype == 0) Then
            write(*,*) 'NumH=',NumH
            write(*,*) 'numzer=',numzero
            write(*,*) 'Hmin=',Hamil%minval
        End If

        deallocate(buf)
        deallocate(zzc)
        deallocate(ccj)
    End Subroutine formh_sym

    Subroutine Hmatrix(icomp,idet1,idet2,t)
        Use integrals, Only : Hint, Gint
        Implicit None

        Integer :: icomp, is, nf, i1, i2, j1, j2, iq, jq, jq0
        Integer, Allocatable, Dimension(:) :: idet1, idet2

        Real(dp) :: t

        call Rspq(idet1,idet2,is,nf,i1,i2,j1,j2)
        t=0.d0
        if (nf.LE.2) then
            if (nf.EQ.2) then
                !    determinants differ by two functions
                t=t+Gint(i2,j2,i1,j1)*is !### det_k goes first!
                if (kCSF.eq.0.and.icomp.EQ.0) t=t+Gj*F_J2(idet1,idet2)
            else
                if (nf.EQ.1) then
                !     determinants differ by one function
                 do iq=1,Ne
                    i1=idet1(iq)
                    if (i1.NE.j1) then
                       t=t+Gint(j2,i1,j1,i1)*is
                    end if
                 end do
                 t=t+Hint(j2,j1)*is
                else
                    !     determinants are equal
                    do iq=1,Ne
                        i1=idet1(iq)
                        jq0=iq+1
                        if (jq0.LE.Ne) then
                           do jq=jq0,Ne
                              j1=idet1(jq)
                              t=t+Gint(i1,j1,i1,j1)*is
                           end do
                        end if
                        t=t+Hint(i1,i1)*is
                    end do
                    if (kCSF.eq.0) t=t+Gj*F_J2(idet1,idet2)
                end if
           end if
        end if
    End Subroutine Hmatrix

    Real(type_real) Function F_J2(idet1, idet2) 
        Implicit None

        Integer, allocatable, dimension(:), intent(InOut) :: idet1, idet2
        Integer :: na, ja, ma, jq0, iq, jq, is, nf, ia, ib, ic, id
        Real(type_real)  :: t

        t=0_type_real
        call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
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

    Subroutine unsym(ncsf,nconf,nccj)
        Implicit None

        Integer :: ncsf, nconf, nccj

        Integer :: jmax, j, idum, i, iconf_neq, nci, ndi, n1, n2, k, n, ic, iccj, iconf, id
        Real(dp), Allocatable, Dimension(:) :: ccj, ccs, cc

        allocate(ccj(nccj))
        allocate(ccs(ncsf))
        allocate(cc(nd))

        ccj(1:nccj) = ccj_module(1:nccj)

        open(unit=17,file='CONF.WFS',status='UNKNOWN',form='UNFORMATTED')
        jmax=nlv
        if (jmax.gt.ncsf) jmax=ncsf
        do j=1,jmax
            read (17) Tk(j),Tj(j),idum,(ccs(i),i=1,ncsf)
            cc(1:Nd)=0.d0
            do iconf=1,nconf
                iconf_neq=nc_neq(iconf)
                nci=ndcs(iconf_neq)
                if (nci.eq.0) cycle
                ndi=ndc(iconf)
                n1=mdc(iconf)+1
                n2=n1+ndi-1
                do n=n1,n2
                    id=n-n1+1
                    do ic=1,nci
                        k=ic+mdcs(iconf)
                        iccj=id+(ic-1)*ndi+iplace_cj(iconf_neq)
                        cc(n)=cc(n)+ccs(k)*ccj(iccj)
                    end do
                end do
            end do
            ArrB(1:Nd,j)=cc(1:Nd)
        end do

        close(unit=17)

        deallocate(cc)
        deallocate(ccs)
        deallocate(ccj)
        If (Allocated(ccj_module)) Deallocate(ccj_module)

    End Subroutine unsym

    Subroutine reorder_det(nconf,idt,idt_orig)
        Implicit None
        Integer :: nconf
        Integer, Allocatable, Dimension(:,:) :: idt, idt_orig
        
        Integer :: iconf, ndi, n1, n2, n, i, k1, k2, k, nf, is, ia, ib, ic, id, ilev
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Real(dp) :: t

        If (.not. Allocated(idet1)) Allocate(idet1(Ne))
        If (.not. Allocated(idet2)) Allocate(idet2(Ne))

        Do iconf=1,nconf
            ndi=ndc(iconf)
            If (ndi == 0) Cycle
            n1=mdc(iconf)+1
            n2=n1+ndc(iconf)-1
            Outer: Do n=n1,n2
                Do i=1,ne
                    idet1(i)=idt_orig(n,i)
                End Do
                k1=mdc(iconf)+1
                k2=k1+ndc(iconf)-1
                Inner: Do k=k1,k2
                    Do i=1,ne
                        idet2(i)=idt(k,i)
                    End Do
                    Call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
                    If (nf /= 0) Cycle Inner
                    If (is /= 1) Then
                        Write(*,*) '  Incorrect sign in subroutine reorder'
                        Write(*,*) '  Det number=',n
                        Stop
                    End If
                    If (k == n) Cycle Outer
                    Do i=1,ne
                        id=idt(n,i)
                        idt(n,i)=idt(k,i)
                        idt(k,i)=id
                    End Do
                    Do ilev=1,Nlv
                        t=ArrB(n,ilev)
                        ArrB(n,ilev)=ArrB(k,ilev)
                        ArrB(k,ilev)=t
                    End Do
                    Cycle Outer
                End Do Inner
                Stop 'Error in subroutine reorder !!!'
            End Do Outer
        End Do

    End Subroutine reorder_det

    Integer Function idif(iconf,jconf)
        Implicit None

        Integer :: iconf, jconf

        Integer :: id, nii, ni, nqi, njj, nj, nqj, i1, i2, i
        Integer, Allocatable, Dimension(:) :: iocc, jocc
        
        id=3
        idif=id

        If (.not. allocated(iocc)) Allocate(iocc(Ns))
        If (.not. allocated(jocc)) Allocate(jocc(Ns))
        iocc(1:Ns)=0
        jocc(1:Ns)=0
        do nii=ni_conf(iconf),nf_conf(iconf)
            if (nii == 0 .or. nii > Nsp) cycle
            ni=nip(nii)
            nqi=nq(nii)
            iocc(ni)=nqi
        enddo
        do njj=ni_conf(jconf),nf_conf(jconf)
            if (njj == 0 .or. njj > Nsp) cycle
            nj=nip(njj)
            nqj=nq(njj)
            jocc(nj)=nqj
        end do
        i1=0
        i2=0
        do ni=1,Ns
          i=iocc(ni)-jocc(ni)
          if (i.gt.0) i1=i1+i
          if (i.lt.0) i2=i2-i
          if (i1.gt.2) return
          if (i2.gt.2) return
        end do
        id=i1
        idif=id
    End Function idif

    Subroutine jbasis_init
        Use mpi_f08
        Use mpi_utils, Only : BroadcastI
        Implicit None

        Integer :: mpierr
        Integer(Kind=int64) :: count

        Call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(nconf_neq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Mj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        If (.not. Allocated(nc_neq)) Allocate(nc_neq(Nc))
        If (.not. Allocated(Ndc)) Allocate(Ndc(Nc))
        If (.not. allocated(Mdc)) Allocate(Mdc(Nc))
        If (.not. allocated(idt)) allocate(idt(Nd,Ne))
        If (.not. allocated(Nh)) Allocate(Nh(Nst))
        If (.not. allocated(Jz)) Allocate(Jz(Nst))

        Call MPI_Bcast(nc_neq, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)   
        Call MPI_Bcast(Ndc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Mdc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        count = Ne*Int(Nd, kind=int64)
        Call BroadcastI(idt, count, 0, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Nh, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jz, Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jj, Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
 
        Return
    End Subroutine jbasis_init

End Module