Module csf

    Use conf_variables
    Use determinants, Only : Rspq
    Use formj2, Only : Plj
    Use davidson, Only : Hould

    Implicit None

    Private

    Public :: jbasis, nonequiv_conf, formh_sym, unsym, reorder_det, read_ccj, idif


Contains

    Subroutine jbasis(nconf,ncsf,nccj,max_ndcs)
        Implicit None
        Integer :: nconf, ncsf, nccj, max_ndcs
        Integer :: iconf_neq, iconf, ndi, ncsfi, ic1, n1, n2, n, k, nf, is, ia, ib, ic, id, iq, na, ja, ma, jq0, jq, i1, i2, j, ifail, jtt
        Real(dp) :: t, tj

        Integer, Allocatable, Dimension(:) ::  ind_conf, idet1, idet2
        real(dp), dimension(:,:), allocatable :: zz
        real(dp), dimension(:), allocatable :: de
        real(dp), dimension(:), allocatable :: dd
        Character(Len=256) :: strfmt

        If (.not. allocated(iplace_cj)) Allocate(iplace_cj(Nc))
        If (.not. allocated(ind_conf)) Allocate(ind_conf(Nc))
        If (.not. allocated(mdcs)) Allocate(mdcs(Nc))
        If (.not. allocated(ndcs)) Allocate(ndcs(Nc))
        If (.not. allocated(ndc_neq)) Allocate(ndc_neq(Nc))
        If (.not. allocated(idet1)) Allocate(idet1(Ne))
        If (.not. allocated(idet2)) Allocate(idet2(Ne))

        ncsf=0
        nccj=0
        max_ndcs=0
        mdcs(1)=0
        iplace_cj(1)=0
        ind_conf(1:nconf_neq)=0

        strfmt = '(/4x,"Calculating matrix J**2")'
        write( *,strfmt)
        write(11,strfmt)

        strfmt = '(6x,"Iconf",8x,"Ndi",5x,"Ncsfi",5x,"Ncsf")'
        write( *,strfmt)
        write(11,strfmt)

        open(unit=18,file='conb.ccj',status='unknown',form='unformatted')

        do iconf=1,nconf
            ndi=ndc(iconf)
            allocate(zz(ndi,ndi))
            allocate(de(ndi))
            allocate(dd(ndi))
            ncsfi=0
            iconf_neq=nc_neq(iconf)
            if (ind_conf(iconf_neq).gt.0) then
                if (ndc(iconf).ne.ndc_neq(iconf_neq)) then
                    write( *,'(/2x,a)') 'Error in jbasis.'
                    write(11,'(/2x,a)') 'Error in jbasis.'
                    stop
                endif
                ncsf=ncsf+ndcs(iconf_neq)
                ncsfi=ndcs(iconf_neq)
                if (iconf.gt.1) then
                    ic1=nc_neq(iconf-1)
                    mdcs(iconf)=mdcs(iconf-1)+ndcs(ic1)
                endif
                if (iconf.eq.nconf) then
                    strfmt = '(i9,i12,2i9,i15)'
                    write( *,strfmt) iconf,ndi,ncsfi,ncsf
                    write(11,strfmt) iconf,ndi,ncsfi,ncsf
                end if
                deallocate(dd)
                deallocate(de)
                deallocate(zz)
                cycle
            endif
            ndc_neq(iconf_neq)=ndi
            zz(1:ndi,1:ndi)=0.d0
            n1=mdc(iconf)+1
            n2=n1+ndi-1
            do n=n1,n2
                idet1(1:ne)=idt(n,1:ne)
                do k=n1,n
                    idet2(1:ne)=idt(k,1:ne)
                    t=0.d0
                    if (k.ne.n) then
                        call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
                        if (nf.ne.2) cycle
                        ! determinants differ by two functions
                        t=plj(ia,ic)*plj(id,ib)+plj(ic,ia)*plj(ib,id)-plj(ia,id)*plj(ic,ib)-plj(id,ia)*plj(ib,ic)
                        t=t*is
                        if (t.eq.0.d0) cycle
                    else
                        t=mj*mj
                        do iq=1,ne
                            ia=idet1(iq)
                            na=nh(ia)
                            ja=jj(na)
                            ma=jz(ia)
                            t=t+ja*(ja+2)-ma**2
                            jq0=iq+1
                            if (jq0.gt.ne) cycle
                            do jq=jq0,ne
                                ib=idet1(jq)
                                t=t-plj(ia,ib)**2-plj(ib,ia)**2
                            end do
                        end do
                        if (t.eq.0.d0) cycle
                    end if
                    i1=n-n1+1
                    i2=k-n1+1
                    zz(i1,i2)=t
                    zz(i2,i1)=t
                end do
            end do
            if (ndi.gt.2500) then
                write( *,'(2x,a,i6)') 'Warning: Ndi=',ndi
                write(11,'(2x,a,i6)') 'Warning: Ndi=',ndi
            end if
            ! Diagonalization
            if (ndi.gt.0) call hould(ndi,dd,de,zz,ifail)
            ndcs(iconf_neq)=0
            do i1=1,ndi
                tj=0.5d0*(dsqrt(1.d0+de(i1))-1.d0)
                jtt=2*tj+0.0001
                if (dabs(2*tj-jtt).gt.1.d-7) then
                    write( *,'(/2x,a/2x,a,f16.8,/2x,a,i3)') '*** Value of j in jbasis is wrong ***','J=',tj,'Configuration:',iconf
                    write(11,'(/2x,a/2x,a,f16.8,/2x,a,i3)') '*** Value of j in jbasis is wrong ***','J=',tj,'Configuration:',iconf
                    stop
                end if
                if (jtt.eq.mj) then
                    ncsf=ncsf+1
                    ncsfi=ncsfi+1
                    ndcs(iconf_neq)=ndcs(iconf_neq)+1
                    nccj=nccj+ndi
                    write(18) (zz(j,i1),j=1,ndi)
                end if
            end do
            if (iconf_neq.gt.1) then
                ic1=iconf_neq-1
                iplace_cj(iconf_neq)=iplace_cj(ic1)+ndcs(ic1)*ndc_neq(ic1)
            end if
            if (iconf.gt.1) then
                ic1=nc_neq(iconf-1)
                mdcs(iconf)=mdcs(iconf-1)+ndcs(ic1)
            end if
            if (ndcs(iconf_neq).gt.max_ndcs) max_ndcs=ndcs(iconf_neq)
            ind_conf(iconf_neq)=1
            strfmt = '(i9,i12,2i9,i15)'
            write( *,strfmt) iconf,ndi,ncsfi,ncsf
            write(11,strfmt) iconf,ndi,ncsfi,ncsf
            deallocate(dd)
            deallocate(de)
            deallocate(zz)
        end do
        close(unit=18)

        write( *,'(/2x,a,i7,/2x,a,i8)') 'Nconf=',nconf,'Ncsf=',ncsf
        write(11,'(/2x,a,i7,/2x,a,i8)') 'Nconf=',nconf,'Ncsf=',ncsf

        if (ncsf.eq.0) then
            write( *,'(/4x,a,f4.1,a)') 'Term ',jm,' is absent in all configurations'
            write(11,'(/4x,a,f4.1,a)') 'Term ',jm,' is absent in all configurations'
            stop
        end if

    End Subroutine jbasis

    Subroutine nonequiv_conf(nconf)
        Implicit None
        Integer :: nconf

        Integer :: iconf, iconf_neq, jconf, njj, nj, nii, ni, nsj, nj0, nsi, ni0, nii1, ni1, i
        Integer, Allocatable, Dimension(:) :: ind_nj
        Character(Len=256) :: strfmt

        write( *,'(/4x,a)') 'Creating list of nonequivalent  configurations...'

        Allocate(ni_conf(nconf), nf_conf(nconf), nc_neq(Nc), ind_nj(Ns))

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

        nconf_neq=iconf_neq
        strfmt = '(4x,"Number of nonequivalent configurations:",i7)'
        write( *,strfmt) nconf_neq
        write(11,strfmt) nconf_neq
    End Subroutine nonequiv_conf

    Subroutine formh_sym(nconf,ncsf,nccj,max_ndcs)
        Use vaccumulator
        Use determinants, Only : calcNd0
        Implicit None

        Integer :: nconf, ncsf, nccj, max_ndcs
        Integer :: numzero, n0, iconf, iconf_neq, nci, ndi, jconf, jconf_neq, ncj, ndj, idf
        Integer :: n1, n2, n, id, k1, k2, k, jd, jc, jc2, ic, iccj, jccj, counter1, counter2, counter3
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer               :: vaGrowBy
        Real(dp) :: Hmin, hij
        Real(dp), Allocatable, Dimension(:) :: ccj, buf
        Real(dp), Allocatable, Dimension(:,:) :: zzc
        Character(Len=256) :: strfmt

        If (.not. allocated(ccj)) allocate(ccj(nccj))
        If (.not. allocated(zzc)) allocate(zzc(max_ndcs,max_ndcs))
        If (.not. allocated(buf)) allocate(buf(max_ndcs))
        If (.not. allocated(idet1)) allocate(idet1(Ne))
        If (.not. allocated(idet2)) allocate(idet2(Ne))

        Call calcNd0(Nc1, Nd0)

        vaGrowBy = 10000000
        Call IVAccumulatorInit(iva1, vaGrowBy)
        Call IVAccumulatorInit(iva2, vaGrowBy)
        Call RVAccumulatorInit(rva1, vaGrowBy)

        ! Reading the list of the symmetrized coefficients
        call read_ccj(ccj)

        close(unit=16)
        NumH=0
        Kherr=0
        Kgerr=0
        numzero=0
        n0=1
        Hmin=0.d0
        if (Ksig.EQ.2) then
           iscr=0
           xscr=0
           write(*,*) 'Screening is included'
        end if

        ! calculation of the matrix elements
        do iconf=1,nconf
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
        Call IVAccumulatorCopy(iva1, Hamil%ind1, counter1)
        Call IVAccumulatorCopy(iva2, Hamil%ind2, counter2)
        Call RVAccumulatorCopy(rva1, Hamil%val, counter3)

        Call IVAccumulatorReset(iva1)
        Call IVAccumulatorReset(iva2)
        Call RVAccumulatorReset(rva1)

        write(*,*) 'NumH=',NumH
        write(*,*) 'numzer=',numzero
        write(*,*) 'Hmin=',Hamil%minval

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
    End Subroutine hmatrix

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

        call read_ccj(ccj)

        open(unit=17,file='CONF.WFS',status='UNKNOWN',form='UNFORMATTED')
        jmax=nlv
        if (jmax.gt.ncsf) jmax=ncsf
        do j=1,jmax
            read (17) E1(j),Tj(j),idum,(ccs(i),i=1,ncsf)
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

    End Subroutine unsym

    Subroutine reorder_det(nconf,idt_orig)
        Implicit None
        Integer :: nconf
        Integer, Allocatable, Dimension(:,:) :: idt_orig
        
        Integer :: iconf, ndi, n1, n2, n, i, k1, k2, k, nf, is, ia, ib, ic, id, ilev
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Real(dp) :: t

        Do iconf=1,nconf
            ndi=ndc(iconf)
            If (ndi.eq.0) Cycle
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
                    If (nf.ne.0) Cycle
                    If (is.ne.1) Then
                        Write(*,*) '  Incorrect sign in subroutine reorder'
                        Write(*,*) '  Det number=',n
                        Stop
                    End If
                    If (k.eq.n) Cycle
                    Do i=1,ne
                        id=idt(n,i)
                        idt(n,i)=idt(k,i)
                        idt(k,i)=id
                    End Do
                    Do ilev=1,Nlv
                        t=ArrB(n,ilev)
                        Arrb(n,ilev)=ArrB(k,ilev)
                        Arrb(k,ilev)=t
                    End Do
                    Cycle Outer
                End Do Inner
                Stop 'Error in subroutine reorder !!!'
            End Do Outer
        End Do

    End Subroutine reorder_det

    Subroutine read_ccj(ccj)
        Implicit None

        Real(dp), Allocatable, Dimension(:) :: ccj

        Integer :: i1, i2, i, ndi, k, iconf_neq

        open(unit=18,file='conb.ccj',status='old',form='unformatted')
        ! Reading the list of the symmetrizied coefficients
        rewind(18)
        i1=1
        do iconf_neq=1,nconf_neq
            ndi=ndc_neq(iconf_neq)
            do k=1,ndcs(iconf_neq)
                i2=i1+ndi-1
                read (18) (ccj(i),i=i1,i2)
                i1=i2+1
            end do
        end do

        close(unit=18)
    End Subroutine read_ccj

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
            ni=nip(nii)
            nqi=nq(nii)
            iocc(ni)=nqi
        enddo
        do njj=ni_conf(jconf),nf_conf(jconf)
            nj=nip(njj)
            nqj=nq(njj)
            jocc(nj)=nqj
        end do
        i1=0
        i2=0
        do ni=1,ns
          i=iocc(ni)-jocc(ni)
          if (i.gt.0) i1=i1+i
          if (i.lt.0) i2=i2-i
          if (i1.gt.2) return
          if (i2.gt.2) return
        end do
        id=i1
        idif=id
    End Function idif

End Module