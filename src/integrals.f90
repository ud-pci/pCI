Module integrals

    Use conf_variables

    Implicit None

    Private

    Public :: Rint, RintS, Hint, HintS, Gint, GintS, Find_VS, Find_SMS, Gaunt

  Contains

    subroutine Rint
        Implicit None
        Integer     :: i, nsh, nx, ns1, nso1, nsu1, Nsul, nsp1, Nlist, &
                       k, mlow, m_is
        Character*7 :: str(3),str1(4)*4
        Data str /'Coulomb','Gaunt  ','Breit  '/
        Data str1 /' VS','SMS','NMS',' MS'/
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        nsh=Nso+1
        Nhint=0
        Ngint=0
        ! parameter for indexation of integrals:
        If (Nsu > IPx+nsh-2) Then
            Write (*,*) ' Rint: IPx is too small'
            Stop
        Else
            nx = IPx
        End If
        Open(unit=13,file='CONF.INT',status='OLD', &
             form='UNFORMATTED',err=700)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
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
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Read (13) (Nn(i),Kk(i),Ll(i),Jj(i), i=1,Nsu)
        Read (13)
        Read (13) Nhint,Kbrt
        Allocate(Rint1(Nhint),Iint1(Nhint))
        Read (13) (Rint1(i), i=1,Nhint)
        Read (13) (Iint1(i), i=1,Nhint)
        Read (13) Ngint,Nlist,nrd
        Allocate(Rint2(IPbr,Ngint),Iint2(Ngint),Iint3(Ngint),IntOrd(nrd))
        If (nrd /=  nx*nx) Then
            Write(*,*)' Rint error: IPx was changed'
            Stop
        End If
        Read (13) ((Rint2(k,i), k=1,IPbr), i=1,Ngint)
        Read (13) (Iint2(i), i=1,Ngint)
        Read (13) (Iint3(i), i=1,Ngint)
        Read (13) (IntOrd(i), i=1,nrd)
        If (K_is >= 1) Then
            Read(13,End=800) m_is,mlow,num_is
            Allocate(R_is(num_is),I_is(num_is))
            If (m_is /=  K_is) Then
                Write(*,*) 'IS Integrals are for K_is=',m_is
                Read(*,*)
                Stop
            End If
            If (K_is >= 2   .and.   mlow /=  Klow) Then
                Write(*,*) 'SMS Integrals are for Klow=',mlow
                Read(*,*)
            End If
            Read(13,End=800) (R_is(i),i=1,num_is)
            Read(13,End=800) (I_is(i),i=1,num_is)
        End If
        Close(unit=13)
        Write( *,'(4X,A7," integrals Read from CONF.INT")') str(Kbrt+1)
        Write(11,'(4X,A7," integrals Read from CONF.INT")') str(Kbrt+1)
        If (K_is >= 1) Then
            Write( *,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
            Write(11,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (Nsp /=  nsp1) Then
            Write ( 6,*) ' Nsp changed since integrals were calculated'
            Write (11,*) ' Nsp changed since integrals were calculated'
            End If
        Return
        !- - - - - - - - - - - - - - - - - - - - - - - - -
  700   Write( 6,'(2X," Can not find file CONF.INT...")')
        Write(11,'(2X," Can not find file CONF.INT...")')
        Stop
  800   Write( *,*) 'No SMS integrals in CONF.INT!'
        Stop
    End subroutine Rint

    Subroutine RintS
        ! Reading of files SGC.CON and SCRC.CON with the self-energy and
        ! screening radial integrals.
        Implicit None
        Integer :: na, nb, ierr, k, nsh, nx, khot, la, lb, ja, jb, na1, nb1, ind, &
                   nso1, khot1, k1, nsx1, nsx2, nav, kbox, i, ns1, idummy
        Real :: x, y, z
        !     - - - - - - - - - - - - - - - - - - - - - - - - -
        NhintS=0
        NgintS=0
        If (Ksig == 0) Return
        nsh=Nso+1
        ! parameter for indexation of integrals:
        nx = IPx
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Reading matrix element of Sigma from SGC.CON
        If (Kdsig /= 0) Then
            Write(*,'(5X,"Give valence energy: ",$)')
            Read(*,*) E_0
            Kexn=1
            Write(*,'(5X,"Extrapolation variant ","(1-normal,2-one side,3-nonlin.): ",I2)') Kexn
        End If
        Open(unit=18,file='SGC.CON',status='OLD',err=700)
        Write(*,*)' Reading file SGC.CON ...'
        Read(18,'(7X,I3,5X,I1,24X,I1)') NmaxS,LmaxS,khot
        Do na=nsh,NmaxS
            la=Ll(na)
            ja=Jj(na)
            If (la > LmaxS) Cycle
            Do nb=na,NmaxS
                lb=Ll(nb)
                jb=Jj(nb)
                If (lb /= la .or. jb /= ja) Cycle
                Read(18,*) idummy,na1,nb1,x,y,z
                NhintS=NhintS+1
                If (na1 /= na .or. nb1 /= nb) Then
                    Write(*,*)' wrong indices for ',NhintS,' sigma:'
                    Write(*,*)' expected ',na,nb,' got ',na1,nb1
                    Stop
                End If
            End Do
        End Do
        rewind(18)
        Allocate(Iint1S(NhintS),Rsig(NhintS),Dsig(NhintS),Esig(NhintS))
        Read(18,'(7X,I3,5X,I1,24X,I1)') NmaxS,LmaxS,khot
        i=0
        Do na=nsh,NmaxS
            la=Ll(na)
            ja=Jj(na)
            If (la > LmaxS) Cycle
            Do nb=na,NmaxS
                lb=Ll(nb)
                jb=Jj(nb)
                If (lb /= la .or. jb /= ja) Cycle
                Read(18,*) idummy,na1,nb1,x,y,z
                i=i+1
                If (na1 /= na .or. nb1 /= nb) Then
                    Write(*,*)' wrong indeces for ',NhintS,' sigma:'
                    Write(*,*)' expected ',na,nb,' got ',na1,nb1
                    Stop
                End If
                ind=nx*(na-nsh)+(nb-nsh+1)
                Iint1S(i)=ind
                Rsig(i)=x
                Dsig(i)=y
                Esig(i)=z
            End Do
        End Do
        xscr=10.d0
        If (Ksig == 2) Then
            xscr=0.d0
            Read(18,'(5X,I2,5X,F8.5)')
            Do k=1,10
                Read(18,'(5X,I2,5X,F8.5)') k1,Scr(k)
                If (k1 /= k-1) Then
                    Write(*,*)' wrong order of screening coefficients'
                    Stop
                End If
                xscr=xscr+Scr(k)
            End Do
        End If
        Write(*,*)' Khot =',khot
        Close(unit=18)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (Ksig < 2) Return
        ierr=0
        Open (unit=13,file='SCRC.CON',status='OLD',form='UNFORMATTED',err=710)
        Write(*,*)' Reading file SCRC.CON...'
        Read (13) ns1,nso1,nsh,khot1,nsx1,nsx2,Ksym
        Read (13) Nmax1,Lmax1,Kmax,nav,kbox
        Write(11,'(4X,"Parameters from the file SCRC.CON:", &
                /4X,"Ns=",I3," Nso=",I2," Nsh=",I2," Khot=",I1, &
                " Nsx=",2I3,/4X," Ksym=",I1,"NmaxS=",I3," LmaxS=",I2, &
                " Kmax=",I2," Nav=",I3," Kbox=",I1)') ns1,nso1,nsh, &
                    khot1,nsx1,nsx2,Ksym,Nmax1,Lmax1,Kmax,nav,kbox
        If (nso1 /= Nso) Then
            Write(*,*)' RintS warning: Nso=',Nso,' <> ',Nso1
            ierr=ierr+1
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (NmaxS /= nmax1 .or. LmaxS /= lmax1 .or. khot /= khot1) Then
            Write( *,'(3X,"DIfference between SGC.CON and SCRC.CON:", &
                /3X,"nmaxS =",I3," (SGC) =",I3," (SCRC)", &
                /3X,"lmaxS =",I3," (SGC) =",I3," (SCRC)", &
                /3X,"khot =",I3," (SGC) =",I3," (SCRC)")') & 
                NmaxS,Nmax1,LmaxS,Lmax1,khot,khot1
            If (NmaxS < Nmax1 .or. LmaxS < Lmax1) ierr=ierr+1
        End If
        Read (13) NgintS,nrd
        Allocate(Rint2S(NgintS),Iint2S(NgintS),Iint3S(NgintS),Dint2S(NgintS),Eint2S(NgintS),IntOrdS(nrd))
        If (nrd /= nx*nx) Then
            Write(*,*)' RintS: IPx was changed'
            ierr=ierr+1
        End If
        If (ierr > 0) Then
            Write(*,*) ierr,' error(s) detected'
            Stop
        End If
        If (nsx2 == 0) Then ! Nsum is used to skip screening integrals
            Nsum=4*nav
        Else
            Nsum=4*Nmax1
        End If
        Read (13) (Rint2S(i), i=1,NgintS)
        Read (13) (Iint2S(i), i=1,NgintS)
        Read (13) (Iint3S(i), i=1,NgintS)
        Read (13) (IntOrdS(i), i=1,nrd)
        Read (13) (Dint2S(i), i=1,NgintS)
        Read (13) (Eint2S(i), i=1,NgintS)
        Close(unit=13)
        Write(*,*)' File SCRC.CON is Read'
        If (dabs(xscr-10.d0) > 1.d-5) Then
            Write(*,*) ' Note: Screening coefficients will be used'
            Write(*,*) ' above NmaxS ',Nmax1,' and LmaxS',Lmax1
            Read(*,*)
        End If
        Return
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
  700   Write(6, '(2X," Can not find file SGC.CON...")')
        Write(11,'(2X," Can not find file SGC.CON...")')
        Stop
  710   Write(6, '(2X," Can not find file SCRC.CON...")')
        Write(11,'(2X," Can not find file SCRC.CON...")')
        Stop
    End subroutine RintS

    Real(dp) Function Hint(ia,ib)
        Implicit None
        Integer  :: nx, n, na, nb, n0, ind, i, nab, ia, ib
        Real(dp) :: e
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        e=0.d0
        nx = IPx
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
                ind=nx*(na-Nso-1)+(nb-Nso)
                Do i=1,Nhint
                    nab=Iint1(i)
                    If (nab == ind) Then
                        e=Rint1(i)
                        If (Ksig /= 0) e = e + HintS(na,nb,n0)
                        If (C_is /= 0.d0) Then
                            If (K_is >= 2 .and. K_sms /= 2) e = e + C_is*Find_SMS(na,nb) ! K_is=3 gives NMS
                            If (K_is == 1) e = e + C_is*Find_VS(na,nb)
                        End If
                        Hint=e
                        Return
                    End If
                End Do
                Write( 6,'(/4X,"one electron integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
                Write(11,'(/4X,"one electron integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
                Stop
            End If
        End If
        Hint=e
        Return
    End Function Hint

    Real(dp) Function Gint(i1,i2,i3,i4)
        Implicit None
        Integer  :: i2, ib, i1, ia, is, nx, la, nd, nc, nb, na, md, mc, mb, ma, &
                    is_sms, ii, i4, id, i3, ic, iac, k1, is_br, i_br, minint, &
                    iab, kmax, kmin, ibd0, iac0, k, jd, jc, jb, ja, ibr, i, ld, &
                    lc, lb, ibd, kac, kbd, mi
        Real(dp) :: e, rabcd
        Logical  :: l_is, l_br, l_pr
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        l_is= .not. (C_is == 0.d0)                ! False If C_is=0
        l_is=l_is .and. (K_is == 2 .or. K_is == 4)  ! False If K_is /= 2,4
        l_is=l_is .and. (K_sms >= 2)              ! False If K_sms<2
        l_br=Kbrt /= 0 .and. IPbr == 2
        e=0.d0
        nx = IPx
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
                iac0=nx*(na-Nso-1)+(nc-Nso)
                ibd0=nx*(nb-Nso-1)+(nd-Nso)
                kmin=iabs(ja-jc)/2+1
                k=iabs(jb-jd)/2+1
                If (kmin < k) kmin=k
                kmax=(ja+jc)/2+1
                k=(jb+jd)/2+1
                If (kmax > k) kmax=k
                iab = nx*(na-Nso-1)+(nb-Nso)
                minint = IntOrd(iab)
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
                        iac=nx*nx*k+iac0
                        ibd=ibd0
                        Do i=minint,ngint
                            mi=i
                            kac=Iint2(i)
                            kbd=Iint3(i)
                            If (kac == iac .and. kbd == ibd) Then
                              exit
                            Else If (i==ngint) Then
                              Write( 6,'(/4X,"integral is absent:"/4X,"K=",I2,2X,"na=",I2,2X,"nb=",I2, &
                                  2X,"nc=",I2,2X,"nd=",I2)') k,na,nb,nc,nd
                              Write(11,'(/4X,"integral is absent:"/4X,"K=",I2,2X,"na=",I2,2X,"nb=",I2, &
                                  2X,"nc=",I2,2X,"nd=",I2)') k,na,nb,nc,nd
                              Write(*,*) 'iac=',iac,' ibd=',ibd
                              Stop
                            End If
                        End Do
                        rabcd=Rint2(1,mi)
                        If (l_br) rabcd=rabcd+is_br*Rint2(IPbr,mi) ! changed 13/3/12
                        If (Ksig >= 2) Then
                            If (max(na,nb,nc,nd) > Nd .or. max(la,lb,lc,ld) > Lmax) Then
                                rabcd=Scr(k+1)*rabcd
                                iscr=iscr+1
                                xscr=xscr+Scr(k+1)
                            End If
                        End If
                        If (k == 1 .and. l_is .and. l_pr) Then
                          rabcd=rabcd-is_sms*C_is*Find_SMS(na,nc)*Find_SMS(nb,nd)
                        End If
                        e=e+is &
                            *Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                            *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                            *rabcd
                        minint = mi + 1
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
        Integer :: n, na, nb, iab, nmin, nmax, kab
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (na <= nb) Then
            iab=IPx*(na-Nso-1)+nb-Nso
        Else
            iab=IPx*(nb-Nso-1)+na-Nso
        End If
        !Write(*,*) ' na=',na,' nb=',nb,' iab=',iab
        nmin=1
        nmax=num_is
        If(I_is(nmax) == iab) Then
            n=nmax
            Find_VS=R_is(n)
            Return
        End If
    1   If (nmin == nmax) Then
            Write(*,*) ' Find_VS error: na=',na,' nb=',nb
            Write(*,*) ' num_is=',num_is,' iab=',iab
            Read(*,*)
            Stop
        End If
        n=(nmax+nmin)/2
        kab=I_is(n)
        ! Write(*,*) ' nmin=',nmin,' nmax=',nmax
        ! Write(*,*) ' n=',n,' kab=',kab
        ! Read(*,*)
        If (iab-kab) 100,300,200
  100   nmax=n
        Goto 1
  200   nmin=n
        Goto 1
  300   Find_VS=R_is(n)
        Return
    End Function Find_VS

    Real(dp) Function Find_SMS(na,nb)
        Implicit None
        Integer  :: is, iab, nmin, nmax, n, kab, na, nb
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (na <= nb) Then
            is=1
            iab=IPx*(na-Nso-1)+nb-Nso
        Else
            is=-1
            iab=IPx*(nb-Nso-1)+na-Nso
        End If
        nmin=1
        nmax=num_is
        If(I_is(nmax) == iab) Then
            n=nmax
            Find_SMS=is*R_is(n)
            Return
        End If
    1   If (nmin == nmax) Then
            Write(*,*) ' Find_SMS error: na=',na,' nb=',nb
            Write(*,*) ' num_is=',num_is,' iab=',iab
            Read(*,*)
            Stop
        End If
        n=(nmax+nmin)/2
        kab=I_is(n)
        If (iab-kab) 100,300,200
  100   nmax=n
        Goto 1
  200   nmin=n
        Goto 1
  300   Find_SMS=is*R_is(n)
        Return
    End Function Find_SMS

    Real(dp) Function HintS(na,nb,n0)
        Implicit None
        Real(dp)    :: e, d, dr, dd, de, d1
        Integer   :: nx, ind, nab, i, na, nb, n0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        e=0.d0
        HintS=e
        nx = IPx
        If (Ll(na) > LmaxS)   Return
        If (nb > NmaxS)       Return
        ind=nx*(na-Nso-1)+(nb-Nso)
        Do i=1,NhintS
            nab=Iint1S(i)
            If (nab == ind) Then
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
                    If (d1 > 0.d0) Then           ! <= NORMAL VARIANT
                        e=e/(1.d0-dd)
                    Else                           ! <= ANOMALOUS VARIANT
                        Select Case(Kexn)
                            Case(1) ! <= two side extrapolation
                                e=e*(1.d0+dd)
                            Case(2)  ! <= one side extrapolation
                                e=e*(1.d0+dmin1(0.d0,dd))
                            Case(3)   ! <= nonlenear extrapolation
                                dr=1.d0+dd-0.1*de*de
                                e=e*dmax1(dr,0.d0)
                        End Select
                    End If
                End If
                HintS=e
                Return
            End If
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Write( 6,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
        Write(11,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
        Stop
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        HintS=e
        Return
    End Function HintS

    Real(dp) Function GintS(i1,i2,i3,i4)
        Implicit None
        Integer   :: is, ia, ib, ic, id, ma, mb, mc, md, na, nb, nc, nd, &
                     na0, nb0, la, lb, lc, ld, i, ja, jb, jc, jd, k, iab, &
                     iac0, ibd0, kmn, kmx, minint, k1, mi, nx, i1, i2, i3, i4, &
                     ii, iac, ibd, kac, kbd
        Real(dp)    :: e, rabcd, dabcd, de, dd, dmin1, dr, dmax1, d1, eabcd
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        e=0.d0
        GintS=e
        nx = IPx
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
            If (max(la,lb,lc,ld) > LmaxS) Return
            If (max(na,nb,nc,nd) > NmaxS) Return
            If ((na+nb+nc+nd) > Nsum)     Return
            i=la+lb+lc+ld
            If (i /= 2*(i/2))             Return
            ja=Jj(na)
            jb=Jj(nb)
            jc=Jj(nc)
            jd=Jj(nd)
            If (na > nc  .and.  nb >= nc  .and.  nd >= nc) Then ! na>nc
                k =na                                     ! nb>=nc
                na=nc                                     ! nd>=nc
                nc=k
                k =nb
                nb=nd
                nd=k
                Goto 200
            End If
            If (na > nb  .and.  nc > nb  .and.  nd >= nb) Then ! na>nb
                k =na                                     ! nc>nb
                na=nb                                     ! nd>=nb
                nb=k
                k =nc
                nc=nd
                nd=k
                Goto 200
            End If
            If (na > nd  .and.  nb > nd  .and.  nc > nd) Then ! na>nd
                k =na                                     ! nb>nd
                na=nd                                     ! nc>nd
                nd=k
                k =nb
                nb=nc
                nc=k
                Goto 210
            End If
  200       If (na == nb  .and.  nc > nd)Then
               k =nc
               nc=nd
               nd=k
            End If
            If (na == nd  .and.  nb > nc)Then
                If (na /= nc) Then
                    k =nc
                    nc=nb
                    nb=k
                Else      ! Case: na=nc=nd<nb
                    k =nd
                    nd=nb
                    nb=k
                End If
            End If
            If (na == nc  .and.  nb > nd)Then
                k =nb
                nb=nd
                nd=k
            End If
    
  210       If (Ksym == 0) Then      ! approx. symmetry
                If (nb > nd) Then     ! which is assumed
                    k =nb               ! when Ksym=0             
                    nb=nd                         
                    nd=k
                End If
                If (na == nb  .and.  nc > nd) Then
                    k =nc
                    nc=nd
                    nd=k
                End If
            End If
    
            iac0=nx*(na-Nso-1)+(nc-Nso)
            ibd0=nx*(nb-Nso-1)+(nd-Nso)
            kmn=max(iabs(ja-jc)/2+1,iabs(jb-jd)/2+1)
            kmx=min((ja+jc)/2+1,(jb+jd)/2+1,Kmax+1)
            iab = nx*(na-Nso-1)+(nb-Nso)
            minint = IntOrdS(iab)
            Do k1=kmn,kmx
                k=k1-1
                iac=nx*nx*k+iac0
                ibd=ibd0
                mi=minint-1
                Do i=minint,NgintS
                    kac=Iint2S(i)
                    kbd=Iint3S(i)
                    If (kac == iac  .and.  kbd == ibd) Then
                        mi=i
                        Goto 230
                    End If
                End Do
  220           continue ! this is to control missing integrals
                Write(*,*) ' GintS: missing integral'
                Write(*,'(4I5,I3)') na,nb,nc,nd,k
                Read(*,*)
                Goto 240
  230           rabcd=Rint2S(i)
                If (rabcd == 0.d0) Then
                    minint = mi + 1
                    Cycle
                End If
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
                    If (d1 > 0.d0) Then           ! <= NORMAL VARIANT
                        rabcd=rabcd/(1.d0-dd)
                    Else                           ! <= ANOMALOUS VARIANT
                        If (Kexn == 1) Then          ! <= two side extr-n
                            rabcd=rabcd*(1.d0+dd)
                        End If
                        If (Kexn == 2) Then          ! <= one side extr-n
                            rabcd=rabcd*(1.d0+dmin1(0.d0,dd))
                        End If
                        If (Kexn == 3) Then          ! <= nonlinear extr-n
                            dr=1.d0+dd-0.1*de*de
                            rabcd=rabcd*dmax1(dr,0.d0)
                        End If
                    End If
                End If
                e=e+is &
                    *Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                    *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                    *rabcd
  240           minint = mi + 1
            End Do
            k=ic
            ic=id
            id=k
            is=-is
        End Do
        GintS=e
        Return
    End Function GintS

    Real(dp) Function Gaunt(k,xj1,xm1,xj2,xm2)
        Implicit None
        Integer  :: i, is, ind, ib1, ib2, im, k, ij
        Real(dp)   :: g, x, xj1, xj2, xm1, xm2, j1, j2, m1, m2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        j1=xj1
        j2=xj2
        m1=xm1
        m2=xm2
        is = 1
        g = 0.d0
        im = abs(m2-m1)
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
                   ij = k+j1-j2
                   If(ij /= 2*(ij/2)) is = -is
                End If
                ib1=2*IPlx+1
                ib2=ib1*ib1
                ind = ib2*(ib2*k+2*(ib1*j1+j2))+ib1*(j1+m1)+(j2+m2)
                Do i=1,IPgnt
                   If(In(i) == ind) Then
                      g = Gnt(i)
                      gaunt = g*is
                      Return
                   End If
                End Do
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

End Module integrals