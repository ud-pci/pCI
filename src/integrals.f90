Module integrals

    Use conf_variables

    Implicit None

    Private

    Public :: Rint, Hint, HintS, Gint, GintS, Find_VS, Find_SMS, Gaunt, BuildGauntLUT, BuildHintLUT, BuildHintSLUT, BuildGintHash, BuildGintSHash, BuildISLUT

  Contains

    subroutine Rint
        Implicit None
        Integer     :: i, nsh, ns1, nso1, nsu1, nsp1, Nlist, &
                       k, mlow, m_is, err_stat, kbrt1
        Integer(Kind=int64) :: i8
        Character*7 :: str(3),str1(4)*4
        Character(Len=3) :: key1, key2

        Data str /'Coulomb','Gaunt  ','Breit  '/
        Data str1 /' VS','SMS','NMS',' MS'/

        nsh=Nso+1
        Nhint=0
        Ngint=0_int64

        Open(unit=13,file='CONF.INT',status='OLD',form='UNFORMATTED',iostat=err_stat)
        If (err_stat /= 0) Then
            Write( 6,'(2X," Can not find file CONF.INT...")')
            Write(11,'(2X," Can not find file CONF.INT...")')
            Stop
        End If
        Write(*,*)' Reading file CONF.INT...'
        Read (13) ns1,nso1,nsp1,Nsu1,Ecore
        If (ns1 /=  Ns) Then
            Write(*,*)' Rint warning: Ns=',Ns,' <> ',Ns1,' push...'
            Read(*,*)
        End If
        If (nso1 /=  Nso) Then
            Write(*,*)' Rint warning: Nso=',Nso,' <> ',Nso1
            Stop
        End If
        If (Nsu1 /=  Nsu) Then
            Write(*,*)' Rint warning: Nsu=',Nsu,' <> ',Nsu1
            If (Nsu1 < Nsu) Stop
        End If
        Write( 6,'(4X,"Total core energy:",F17.7)') Ecore
        Write(11,'(4X,"Total core energy:",F17.7)') Ecore

        Read (13) (Nn(i),Kk(i),Ll(i),Jj(i), i=1,Nsu)
        Read (13)
        Read (13) Nhint,kbrt1
        ! Handle mismatched key for inclusion of Breit corrections
        If (Kbrt /= kbrt1) Then
            Write(key1, '(I1)') kbrt1
            Write(key2, '(I1)') Kbrt
            Write(*,*) ' Rint error: Key for Breit was changed from Kbrt=', Trim(AdjustL(key1)), ' to Kbrt=', Trim(AdjustL(key2))
            Write(*,*) ' Rint error: Please re-run pbasc to update radial integrals with Kbrt=', Trim(AdjustL(key2))
            Stop
        End If
        Allocate(Rint1(Nhint),Iint1(Nhint))
        Read (13) (Rint1(i), i=1,Nhint)
        Read (13) (Iint1(i), i=1,Nhint)
        Read (13) Ngint,Nlist,nrd

        Nx = int(sqrt(real(nrd)))

        If (Kbrt == 0) Then
            Allocate(Rint2(1,Ngint))
        Else
            Allocate(Rint2(2,Ngint))
        End If
        Allocate(Iint2(Ngint),Iint3(Ngint),IntOrd(nrd))
        
        If (Kbrt == 0) Then
            Read (13) (Rint2(1,i8), i8=1,Ngint)
        Else
            Read (13) ((Rint2(k,i8), k=1,2), i8=1,Ngint)
        End If
        Read (13) (Iint2(i8), i8=1,Ngint)
        Read (13) (Iint3(i8), i8=1,Ngint)
        Read (13) (IntOrd(i), i=1,nrd)
        If (K_is >= 1) Then
            Read(13,iostat=err_stat) m_is,mlow,num_is
            If (err_stat /= 0) Then
                Write( *,*) 'No SMS integrals in CONF.INT!'
                Stop
            End If
            Allocate(R_is(num_is),I_is(num_is))
            If (m_is /=  K_is) Then
                Write(*,*) 'IS Integrals are for K_is=',m_is
                Read(*,*)
                Stop
            End If
            If (K_is >= 2 .and. mlow /=  Klow) Then
                Write(*,*) 'SMS Integrals are for Klow=',mlow
                Read(*,*)
            End If
            Read(13,iostat=err_stat) (R_is(i),i=1,num_is)
            Read(13,iostat=err_stat) (I_is(i),i=1,num_is)
            If (err_stat /= 0) Then
                Write( *,*) 'No SMS integrals in CONF.INT!'
                Stop
            End If
        End If
        Close(unit=13)
        Write( *,'(4X,A7," integrals read from CONF.INT")') str(Kbrt+1)
        Write(11,'(4X,A7," integrals read from CONF.INT")') str(Kbrt+1)
        If (K_is >= 1) Then
            Write( *,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
            Write(11,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
        End If

        If (Nsp /=  nsp1) Then
            Write ( 6,*) ' Nsp changed since integrals were calculated'
            Write (11,*) ' Nsp changed since integrals were calculated'
            End If
        Return
    End subroutine Rint

    Real(dp) Function Hint(ia,ib)
        Implicit None
        Integer  :: n, na, nb, n0, ind, ia, ib
        Real(dp) :: e
        e=0.d0

        If (Jz(ia) == Jz(ib)) Then
            na=Nh(ia)
            nb=Nh(ib)
            n0=na
            If (Kk(na) == Kk(nb)) Then
                If (na > nb) Then
                    n=na
                    na=nb
                    nb=n
                End If
                ind=Nx*(na-Nso-1)+(nb-Nso)
                e=HintLUT(ind)
                If (Ksig /= 0) e = e + HintS(na,nb,n0)
                If (C_is /= 0.d0) Then
                    If (K_is >= 2 .and. K_sms /= 2) e = e + C_is*Find_SMS(na,nb) ! K_is=3 gives NMS
                    If (K_is == 1) e = e + C_is*Find_VS(na,nb)
                End If
                Hint=e
                Return
            End If
        End If
        Hint=e
        Return
    End Function Hint

    Real(dp) Function Gint(i1,i2,i3,i4)
        Implicit None
        Integer  :: i2, ib, i1, ia, is, la, nd, nc, nb, na, md, mc, mb, ma, &
                    is_sms, ii, i4, id, i3, k1, is_br, i_br, kmax, kmin, &
                    ic, iac, ibd0, iac0, k, jd, jc, jb, ja, ibr, i, ld, &
                    lc, lb, ibd
        Integer(Kind=int64) :: mi8
        Real(dp) :: e, rabcd
        Logical  :: l_is, l_br, l_pr
        Character(Len=256) :: strfmt

        l_is= .not. (C_is == 0.d0)                ! False If C_is=0
        l_is=l_is .and. (K_is == 2 .or. K_is == 4)  ! False If K_is /= 2,4
        l_is=l_is .and. (K_sms >= 2)              ! False If K_sms<2
        l_br=Kbrt /= 0

        e=0.d0
        is=1
        ia=i1
        ib=i2
        ic=i3
        id=i4
        Do ii=1,2
            is_sms=1 !### extra phase for (p_i Dot p_k)
            ma=Jz(ia)
            mb=Jz(ib)
            mc=Jz(ic)
            md=Jz(id)
            If (ma+mb /= mc+md) Then
                Write(*,*) ' Gint error: ma+mb=',ma+mb,' mc+md=',mc+md
                Read(*,*)
                Stop
            End If
            na=Nh(ia)
            nb=Nh(ib)
            nc=Nh(ic)
            nd=Nh(id)
            la=Ll(na)
            lb=Ll(nb)
            lc=Ll(nc)
            ld=Ll(nd)
            i=la+lb+lc+ld
            ibr=0                     ! defines phase of Breit integral
            If (i == 2*(i/2)) Then
                ja=Jj(na)
                jb=Jj(nb)
                jc=Jj(nc)
                jd=Jj(nd)
                If (na > nc) Then
                    k =na
                    na=nc
                    nc=k
                    is_sms=-is_sms
                    ibr=ibr+1
                End If
                If (nb > nd) Then
                    k =nb
                    nb=nd
                    nd=k
                    is_sms=-is_sms
                    ibr=ibr+1
                End If
                If (na > nb) Then
                    k =na
                    na=nb
                    nb=k
                    k =nc
                    nc=nd
                    nd=k
                End If
                If (na == nb   .and.   nc > nd)Then
                    k =nc
                    nc=nd
                    nd=k
                End If
                If (ibr == 1) Then
                    ibr=-1
                Else
                    ibr=1
                End If
                iac0=Nx*(na-Nso-1)+(nc-Nso)
                ibd0=Nx*(nb-Nso-1)+(nd-Nso)
                kmin=iabs(ja-jc)/2+1
                k=iabs(jb-jd)/2+1
                If (kmin < k) kmin=k
                kmax=(ja+jc)/2+1
                k=(jb+jd)/2+1
                If (kmax > k) kmax=k
                i_br=la+lc+kmin-1
                If (i_br /= 2*(i_br/2)) Then
                    is_br=ibr
                Else
                    is_br=1
                End If
                Do k1=kmin,kmax
                    k=k1-1
                    i=k+la+lc
                    l_pr=i == 2*(i/2)
                    If (l_pr .or. l_br) Then
                        is_br=is_br*ibr
                        iac=Nx*Nx*k+iac0
                        ibd=ibd0
                        mi8 = FindGint(iac, ibd)
                        If (mi8 == 0) Then
                            strfmt = '(/4X,"Gint: integral is absent:"/4X,"K=",I2,2X,"na=",I3,2X,"nb=",I3,2X,"nc=",I3,2X,"nd=",I3)'
                            Write( 6,strfmt) k,na,nb,nc,nd
                            Write(11,strfmt) k,na,nb,nc,nd
                            Stop
                        End If
                        rabcd=Rint2(1,mi8)
                        If (l_br) rabcd=rabcd+is_br*Rint2(2,mi8)
                        If (Ksig >= 2) Then
                            If (max(na,nb,nc,nd) > Nd .or. max(la,lb,lc,ld) > Lmax) Then
                                If (k < 10) Then
                                    rabcd=Scr(k+1)*rabcd
                                    iscr=iscr+1
                                    xscr=xscr+Scr(k+1)
                                End If
                            End If
                        End If
                        If (k == 1 .and. l_is .and. l_pr) Then
                          rabcd=rabcd-is_sms*C_is*Find_SMS(na,nc)*Find_SMS(nb,nd)
                        End If
                        e=e+is &
                            *Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                            *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                            *rabcd
                    End If
                End Do
            End If
            k=ic
            ic=id
            id=k
            is=-is
        End Do
        If (Ksig >= 2) Then
            e = e + GintS(i1,i2,i3,i4)
        End If
        Gint=e
        Return
    End Function Gint

    Real(dp) Function Find_VS(na,nb)
        Implicit None
        Integer :: iab, na, nb

        If (na <= nb) Then
            iab=Nx*(na-Nso-1)+nb-Nso
        Else
            iab=Nx*(nb-Nso-1)+na-Nso
        End If
        Find_VS = ISLUT(iab)
    End Function Find_VS

    Real(dp) Function Find_SMS(na,nb)
        Implicit None
        Integer  :: is, iab, na, nb

        If (na <= nb) Then
            is=1
            iab=Nx*(na-Nso-1)+nb-Nso
        Else
            is=-1
            iab=Nx*(nb-Nso-1)+na-Nso
        End If
        Find_SMS = is * ISLUT(iab)
    End Function Find_SMS

    Real(dp) Function HintS(na,nb,n0)
        Implicit None
        Real(dp)  :: e, d, dr, dd, de, d1
        Integer   :: ind, i, na, nb, n0

        e=0.d0
        HintS=e

        If (Ll(na) > LmaxS .or. nb > NmaxS) Return
        ind=Nx*(na-Nso-1)+(nb-Nso)
        i = HintSPosLUT(ind)
        If (i == 0) Then
            Write( 6,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
            Write(11,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
            Stop
        End If
        e=Rsig(i)
        d=Dsig(i)

        If (Kdsig*d /= 0.d0) Then
            de=(E_k+Eps(n0)-Esig(i))
            d1=d/e
            dd=d1*de

            If (dd > 0.5d0) Then
                Kherr=Kherr+1
                dd=0.5d0
            End If

            If (d1 > 0.d0) Then
                ! NORMAL VARIANT
                e=e/(1.d0-dd)
            Else
                ! ANOMALOUS VARIANT
                Select Case(Kexn)
                    Case(1) ! two-side extrapolation
                        e=e*(1.d0+dd)
                    Case(2) ! one-side extrapolation
                        e=e*(1.d0+dmin1(0.d0,dd))
                    Case(3) ! nonlinear extrapolation
                        dr=1.d0+dd-0.1*de*de
                        e=e*dmax1(dr,0.d0)
                End Select
            End If
        End If

        HintS=e
        Return
    End Function HintS

    Real(dp) Function GintS(i1,i2,i3,i4)
        Implicit None
        Integer   :: is, ia, ib, ic, id, ma, mb, mc, md, na, nb, nc, nd, &
                     na0, nb0, la, lb, lc, ld, i, ja, jb, jc, jd, k, &
                     iac0, ibd0, kmn, kmx, k1, i1, i2, i3, i4, &
                     ii, iac, ibd
        Real(dp)  :: e, rabcd, dabcd, de, dd, dmin1, dr, dmax1, d1, eabcd

        e=0.d0
        GintS=e
        is=1
        ia=i1
        ib=i2
        ic=i3
        id=i4

        Do ii=1,2
            ma=Jz(ia)
            mb=Jz(ib)
            mc=Jz(ic)
            md=Jz(id)
            na=Nh(ia)
            nb=Nh(ib)
            nc=Nh(ic)
            nd=Nh(id)
            na0=na
            nb0=nb
            la=Ll(na)
            lb=Ll(nb)
            lc=Ll(nc)
            ld=Ll(nd)

            If (max(la,lb,lc,ld) > LmaxS .or. max(na,nb,nc,nd) > NmaxS .or. (na+nb+nc+nd) > Nsum) Return
            i=la+lb+lc+ld
            If (i /= 2*(i/2)) Return

            ja=Jj(na)
            jb=Jj(nb)
            jc=Jj(nc)
            jd=Jj(nd)

            If (na > nc .and. nb >= nc .and. nd >= nc) Then
                Call SwapValues(na, nc)
                Call SwapValues(nb, nd)
            Else If (na > nb .and. nc > nb .and. nd >= nb) Then
                Call SwapValues(na, nb)
                Call SwapValues(nc, nd)
            Else If (na > nd  .and.  nb > nd  .and.  nc > nd) Then
                Call SwapValues(na, nd)
                Call SwapValues(nb, nc)
            End If

            If (na == nb .and. nc > nd) Call SwapValues(nc, nd)
            If (na == nd  .and.  nb > nc)Then
                If (na /= nc) Then
                    Call SwapValues(nc, nb)
                Else
                    Call SwapValues(nd, nb)
                End If
            End If
            If (na == nc .and. nb > nd) Call SwapValues(nb, nd)
    
            If (Ksym == 0) Then      ! approximate symmetry which is assumed when Ksym=0
                If (nb > nd) Call SwapValues(nb, nd)
                If (na == nb .and. nc > nd) Call SwapValues(nc, nd)
            End If
    
            iac0=Nx*(na-Nso-1)+(nc-Nso)
            ibd0=Nx*(nb-Nso-1)+(nd-Nso)
            kmn=max(iabs(ja-jc)/2+1,iabs(jb-jd)/2+1)
            kmx=min((ja+jc)/2+1,(jb+jd)/2+1,Kmax+1)

            Do k1=kmn,kmx
                k=k1-1
                iac=Nx*Nx*k+iac0
                ibd=ibd0
                i = FindGintS(iac, ibd)
                If (i == 0) Then
                    Write(*,*) ' GintS: missing integral'
                    Write(*,'(4I5,I3)') na,nb,nc,nd,k
                    Stop
                End If

                rabcd=Rint2S(i)
                If (rabcd /= 0.d0) Then
                    dabcd=Dint2S(i)
                    If (Kdsig*dabcd /= 0.d0) Then
                        eabcd=Eint2S(i)
                        de=E_k+Eps(na0)+Eps(nb0)-eabcd
                        d1=dabcd/rabcd
                        dd=d1*de

                        If (dd > 0.5d0) Then
                            Kgerr=Kgerr+1
                            dd=0.5d0
                        End If

                        If (d1 > 0.d0) Then
                            ! NORMAL VARIANT
                            rabcd=rabcd/(1.d0-dd)
                        Else
                            ! ANOMALOUS VARIANT
                            Select Case(Kexn)
                                Case(1) ! two-side extrapolation
                                    rabcd=rabcd*(1.d0+dd)
                                Case(2) ! one-side extrapolation
                                    rabcd=rabcd*(1.d0+dmin1(0.d0,dd))
                                Case(3) ! nonlinear extrapolation
                                    dr=1.d0+dd-0.1*de*de
                                    rabcd=rabcd*dmax1(dr,0.d0)
                            End Select
                        End If
                    End If

                    e=e+is*Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                        *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                        *rabcd
                End If
            End Do

            Call SwapValues(ic, id)
            is=-is
        End Do

        GintS=e
        Return
    End Function GintS

    Subroutine SwapValues(a, b)
        Implicit None
        Integer, Intent(InOut) :: a, b
        Integer :: temp

        temp = a
        a = b
        b = temp
    End Subroutine SwapValues

    Real(dp) Function Gaunt(k,xj1,xm1,xj2,xm2)
        Implicit None
        Integer  :: is, ind, ib1, ib2, im, k, ij
        Real(dp)   :: g, x, xj1, xj2, xm1, xm2, j1, j2, m1, m2

        j1=xj1
        j2=xj2
        m1=xm1
        m2=xm2
        is = 1
        g = 0.d0
        im = Int(abs(m2-m1))
        If (k >= im) Then
            If (k == 0) Then
                If (j1 == j2) Then
                    g = 1.d0
                    gaunt = g*is
                    Return
                End If
            Else
                If (j2 < j1) Then
                   x = j2
                   j2 = j1
                   j1 = x
                   x = m2
                   m2 = m1
                   m1 = x
                   If (im /= 2*(im/2)) is = -is
                End If
                If (m1 <= 0.d0) Then
                   m1 = -m1
                   m2 = -m2
                   ij = Int(k+j1-j2)
                   If(ij /= 2*(ij/2)) is = -is
                End If
                ib1=2*Nlx+1
                ib2=ib1*ib1
                ind = Int(ib2*(ib2*k+2*(ib1*j1+j2))+ib1*(j1+m1)+(j2+m2))
                If (ind >= 0 .and. ind <= ubound(GauntLUT,1)) Then
                    g = GauntLUT(ind)
                    gaunt = g*is
                    Return
                End If
            End If
        End If
        If (K_gnt == 1) Then
            Write (*,'(1X,"Gaunt: Can not find Gaunt for k=",I2, &
               " j1=",F4.1," m1=",F5.1," j2=",F4.1," m2=",F5.1)') k,xj1,xm1,xj2,xm2
            Stop
        End If
        gaunt = g*is
        Return
    End Function Gaunt

    Integer(Kind=int64) Function FindGint(iac, ibd)
        ! O(1) average-case lookup replacing Gint's forward scan over Iint2/Iint3.
        ! Returns 1-based position into Iint2/Iint3/Rint2, or 0 if absent.
        Implicit None
        Integer, Intent(In) :: iac, ibd
        Integer(Kind=int64) :: h, slot

        h = Int(iac,8)*104729_8 + Int(ibd,8)
        slot = Mod(h, GintHashCap)
        If (slot < 0) slot = slot + GintHashCap
        Do While (GintHashPos(slot) /= 0)
            If (GintHashIac(slot) == iac .and. GintHashIbd(slot) == ibd) Then
                FindGint = GintHashPos(slot)
                Return
            End If
            slot = Mod(slot+1_int64, GintHashCap)
        End Do
        FindGint = 0_int64
    End Function FindGint

    Integer Function FindGintS(iac, ibd)
        ! O(1) average-case lookup replacing GintS's forward scan over Iint2S/Iint3S.
        ! Returns 1-based position into Iint2S/Iint3S/Rint2S, or 0 if absent.
        Implicit None
        Integer, Intent(In) :: iac, ibd
        Integer :: slot
        Integer(Kind=int64) :: h

        h = Int(iac,8)*104729_8 + Int(ibd,8)
        slot = Int(Mod(h, Int(GintSHashCap,8)))
        If (slot < 0) slot = slot + GintSHashCap
        Do While (GintSHashPos(slot) /= 0)
            If (GintSHashIac(slot) == iac .and. GintSHashIbd(slot) == ibd) Then
                FindGintS = GintSHashPos(slot)
                Return
            End If
            slot = Mod(slot+1, GintSHashCap)
        End Do
        FindGintS = 0
    End Function FindGintS

    Subroutine BuildGintHash
        ! Hash table over (Iint2(i), Iint3(i)) for i=1...Ngint, keyed by the (iac,ibd) pair Gint computes per k. 
        ! Array sized >= 2*Ngint so at most half the slots are filled; 
        ! collisions are resolved by stepping forward until an empty slot is found.
        ! GintHashPos(slot)==0 marks empty; valid entries are 1-based positions.
        Implicit None
        Integer :: slot32
        Integer(Kind=int64) :: i, h, slot, cap

        cap = 4_int64
        Do While (cap < 2_int64*max(Ngint,1_int64))
            cap = cap * 2_int64
        End Do
        GintHashCap = cap
        If (Allocated(GintHashPos)) Deallocate(GintHashIac, GintHashIbd, GintHashPos)
        Allocate(GintHashIac(0:cap-1), GintHashIbd(0:cap-1), GintHashPos(0:cap-1))
        GintHashPos = 0_int64
        Do i = 1, Ngint
            h = Int(Iint2(i),8)*104729_8 + Int(Iint3(i),8)
            slot = Mod(h, cap)
            If (slot < 0) slot = slot + cap
            Do While (GintHashPos(slot) /= 0)
                slot = Mod(slot+1_int64, cap)
            End Do
            GintHashIac(slot) = Iint2(i)
            GintHashIbd(slot) = Iint3(i)
            GintHashPos(slot) = i
        End Do
    End Subroutine BuildGintHash

    Subroutine BuildGintSHash
        ! Hash table over (Iint2S(i), Iint3S(i)) for i=1...NgintS, keyed by (iac,ibd).
        ! Array sized >= 2*NgintS so at most half the slots are filled; 
        ! collisions are resolved by stepping forward until an empty slot is found.
        ! GintSHashPos==0 marks empty.
        Implicit None
        Integer :: i, slot
        Integer(Kind=int64) :: h

        GintSHashCap = 4
        Do While (GintSHashCap < 2*max(NgintS,1))
            GintSHashCap = GintSHashCap * 2
        End Do
        If (Allocated(GintSHashPos)) Deallocate(GintSHashIac, GintSHashIbd, GintSHashPos)
        Allocate(GintSHashIac(0:GintSHashCap-1), GintSHashIbd(0:GintSHashCap-1), GintSHashPos(0:GintSHashCap-1))
        GintSHashPos = 0
        Do i = 1, NgintS
            h = Int(Iint2S(i),8)*104729_8 + Int(Iint3S(i),8)
            slot = Int(Mod(h, Int(GintSHashCap,8)))
            If (slot < 0) slot = slot + GintSHashCap
            Do While (GintSHashPos(slot) /= 0)
                slot = Mod(slot+1, GintSHashCap)
            End Do
            GintSHashIac(slot) = Iint2S(i)
            GintSHashIbd(slot) = Iint3S(i)
            GintSHashPos(slot) = i
        End Do
    End Subroutine BuildGintSHash

    Subroutine BuildHintLUT
        ! Build a direct-index LUT for one-electron integrals.
        ! Iint1(i) encodes the shell pair (na,nb) as ind=Nx*(na-Nso-1)+(nb-Nso);
        ! HintLUT(ind)=Rint1(i) replaces the O(Nhint) linear scan in Hint with O(1).
        Implicit None
        Integer :: i
        If (.not. Allocated(HintLUT)) Allocate(HintLUT(Nx*Nx))
        HintLUT = 0.0_dp
        Do i = 1, Nhint
            HintLUT(Iint1(i)) = Rint1(i)
        End Do
    End Subroutine BuildHintLUT

    Subroutine BuildHintSLUT
        ! Build a position LUT for one-electron sigma integrals.
        ! HintSPosLUT(ind)=i so HintS can look up Rsig(i)/Dsig(i)/Esig(i) in O(1)
        ! instead of scanning Iint1S(1:NhintS). Value 0 means absent (no entry).
        Implicit None
        Integer :: i
        If (.not. Allocated(HintSPosLUT)) Allocate(HintSPosLUT(Nx*Nx))
        HintSPosLUT = 0
        Do i = 1, NhintS
            HintSPosLUT(Iint1S(i)) = i
        End Do
    End Subroutine BuildHintSLUT

        Subroutine BuildISLUT
        ! Direct-index LUT for IS corrections (VS and SMS).
        ! ISLUT(iab)=R_is(n) for entries in I_is, replacing O(log num_is) binary search in Find_VS and Find_SMS with O(1) array lookup.
        ! Both Find_VS and Find_SMS share this table (same I_is/R_is data).
        Implicit None
        Integer :: n
        If (.not. Allocated(ISLUT)) Allocate(ISLUT(Nx*Nx))
        ISLUT = 0.0_dp
        Do n = 1, num_is
            ISLUT(I_is(n)) = Real(R_is(n), dp)
        End Do
    End Subroutine BuildISLUT

    Subroutine BuildGauntLUT
        ! Build a direct-lookup table for Gaunt coefficients. 
        ! In(ig) encodes the quantum numbers (k, j1, m1, j2, m2) as a single integer index; 
        ! GauntLUT maps that index to Gnt(ig), replacing an O(Ngaunt) linear scan in Gaunt() with an O(1) array lookup. 
        Implicit None
        Integer :: ig, max_ind_g

        max_ind_g = maxval(In(1:Ngaunt))
        If (.not. Allocated(GauntLUT)) Allocate(GauntLUT(0:max_ind_g))
        GauntLUT = 0.0_dp
        Do ig = 1, Ngaunt
            GauntLUT(In(ig)) = Gnt(ig)
        End Do
    End Subroutine BuildGauntLUT

End Module integrals