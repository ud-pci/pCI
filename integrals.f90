module integrals

    use conf_variables

    implicit none

    contains
    
    real(dp) function Hint(ia,ib)
      implicit none
      integer  :: nx, n, na, nb, n0, ind, i, nab, ia, ib
      real(dp) :: e
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      e=0.d0
      nx = IPx
      if (Jz(ia) == Jz(ib)) then
        na=Nh(ia)
        nb=Nh(ib)
        n0=na
        if (Kk(na) == Kk(nb)) then
          if (na > nb) then
            n=na
            na=nb
            nb=n
          end if
          ind=nx*(na-Nso-1)+(nb-Nso)
          do i=1,Nhint
            nab=Iint1(i)
            if (nab == ind) then
              e=Rint1(i)
              if (Ksig /= 0) e = e + HintS(na,nb,n0)
              if (C_is /= 0.d0) then
                if (K_is >= 2 .and. K_sms /= 2) e = e + C_is*Find_SMS(na,nb) ! K_is=3 gives NMS
                if (K_is == 1) e = e + C_is*Find_VS(na,nb)
              end if
              Hint=e
              return
            end if
          end do
          write( 6,'(/4X,"one electron integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
          write(11,'(/4X,"one electron integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
          stop
        end if
      end if
      Hint=e
      return
    end function Hint

    real(dp) function Gint(i1,i2,i3,i4)
        implicit none
        logical  :: l_is, l_br, l_pr
        integer  :: i2, ib, i1, ia, is, nx, la, nd, nc, nb, na, md, mc, mb, ma, &
                    is_sms, ii, i4, id, i3, ic, iac, k1, is_br, i_br, minint, &
                    iab, kmax, kmin, ibd0, iac0, k, jd, jc, jb, ja, ibr, i, ld, &
                    lc, lb, ibd, kac, kbd, mi
        real(dp)   :: e, rabcd
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        l_is= .not. (C_is == 0.d0)                ! False if C_is=0
        l_is=l_is .and. (K_is == 2 .or. K_is == 4)  ! False if K_is /= 2,4
        l_is=l_is .and. (K_sms >= 2)              ! False if K_sms<2
        l_br=Kbrt /= 0 .and. IPbr == 2
        e=0.d0
        nx = IPx
        is=1
        ia=i1
        ib=i2
        ic=i3
        id=i4
        do ii=1,2
            is_sms=1 !### extra phase for (p_i dot p_k)
            ma=Jz(ia)
            mb=Jz(ib)
            mc=Jz(ic)
            md=Jz(id)
            if (ma+mb /= mc+md) then
                write(*,*) ' Gint error: ma+mb=',ma+mb,' mc+md=',mc+md
                read(*,*)
                stop
            end if
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
            if (i == 2*(i/2)) then
                ja=Jj(na)
                jb=Jj(nb)
                jc=Jj(nc)
                jd=Jj(nd)
                if (na > nc) then
                    k =na
                    na=nc
                    nc=k
                    is_sms=-is_sms
                    ibr=ibr+1
                end if
                if (nb > nd) then
                    k =nb
                    nb=nd
                    nd=k
                    is_sms=-is_sms
                    ibr=ibr+1
                end if
                if (na > nb) then
                    k =na
                    na=nb
                    nb=k
                    k =nc
                    nc=nd
                    nd=k
                end if
                if (na == nb   .and.   nc > nd)then
                    k =nc
                    nc=nd
                    nd=k
                end if
                if (ibr == 1) then
                    ibr=-1
                else
                    ibr=1
                end if
                iac0=nx*(na-Nso-1)+(nc-Nso)
                ibd0=nx*(nb-Nso-1)+(nd-Nso)
                kmin=iabs(ja-jc)/2+1
                k=iabs(jb-jd)/2+1
                if (kmin < k) kmin=k
                kmax=(ja+jc)/2+1
                k=(jb+jd)/2+1
                if (kmax > k) kmax=k
                iab = nx*(na-Nso-1)+(nb-Nso)
                minint = IntOrd(iab)
                i_br=la+lc+kmin-1
                if (i_br /= 2*(i_br/2)) then
                    is_br=ibr
                else
                    is_br=1
                end if 
                do k1=kmin,kmax
                    k=k1-1
                    i=k+la+lc
                    l_pr=i == 2*(i/2)
                    if (l_pr .or. l_br) then
                        is_br=is_br*ibr
                        iac=nx*nx*k+iac0
                        ibd=ibd0
                        do i=minint,ngint
                            mi=i
                            kac=Iint2(i)
                            kbd=Iint3(i)
                            if (kac == iac .and. kbd == ibd) then
                              exit
                            else if (i==ngint) then
                              write( 6,'(/4X,"integral is absent:"/4X,"K=",I2,2X,"na=",I2,2X,"nb=",I2, &
                                  2X,"nc=",I2,2X,"nd=",I2)') k,na,nb,nc,nd
                              write(11,'(/4X,"integral is absent:"/4X,"K=",I2,2X,"na=",I2,2X,"nb=",I2, &
                                  2X,"nc=",I2,2X,"nd=",I2)') k,na,nb,nc,nd
                              write(*,*) 'iac=',iac,' ibd=',ibd
                              stop
                            end if
                        end do
                        rabcd=Rint2(1,mi)
                        if (l_br) rabcd=rabcd+is_br*Rint2(IPbr,mi) ! changed 13/3/12
                        if (Ksig >= 2) then
                            if (max(na,nb,nc,nd) > Nd .or. max(la,lb,lc,ld) > Lmax) then
                                rabcd=Scr(k+1)*rabcd
                                iscr=iscr+1
                                xscr=xscr+Scr(k+1)
                            end if
                        end if
                        if (k == 1 .and. l_is .and. l_pr) then
                          rabcd=rabcd-is_sms*C_is*Find_SMS(na,nc)*Find_SMS(nb,nd)
                        end if
                        e=e+is &
                         *Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                         *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                         *rabcd
                        minint = mi + 1
                    end if
                end do
            end if
    
            k=ic
            ic=id
            id=k
            is=-is
        end do
        if (Ksig >= 2) then
            e = e + GintS(i1,i2,i3,i4)
        end if
        Gint=e
        return
    end function Gint

    real(dp) function Find_VS(na,nb)
        implicit none
        integer :: n, na, nb, iab, nmin, nmax, kab
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (na <= nb) then
            iab=IPx*(na-Nso-1)+nb-Nso
        else
            iab=IPx*(nb-Nso-1)+na-Nso
        end if
        !write(*,*) ' na=',na,' nb=',nb,' iab=',iab
        nmin=1
        nmax=num_is
        if(I_is(nmax) == iab) then
            n=nmax
            Find_VS=R_is(n)
            return
        end if
    1   if (nmin == nmax) then
            write(*,*) ' Find_VS error: na=',na,' nb=',nb
            write(*,*) ' num_is=',num_is,' iab=',iab
            read(*,*)
            stop
        end if
        n=(nmax+nmin)/2
        kab=I_is(n)
!            write(*,*) ' nmin=',nmin,' nmax=',nmax
!            write(*,*) ' n=',n,' kab=',kab
!            read(*,*)
        if (iab-kab) 100,300,200
    100 nmax=n
        goto 1
    200 nmin=n
        goto 1
    300 Find_VS=R_is(n)
        return
    end function Find_VS

    real(dp) function Find_SMS(na,nb)
       implicit none
       integer  :: is, iab, nmin, nmax, n, kab, na, nb
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (na <= nb) then
            is=1
            iab=IPx*(na-Nso-1)+nb-Nso
        else
            is=-1
            iab=IPx*(nb-Nso-1)+na-Nso
        end if
        nmin=1
        nmax=num_is
        if(I_is(nmax) == iab) then
            n=nmax
            Find_SMS=is*R_is(n)
            return
        end if
    1   if (nmin == nmax) then
            write(*,*) ' Find_SMS error: na=',na,' nb=',nb
            write(*,*) ' num_is=',num_is,' iab=',iab
            read(*,*)
            stop
        end if
        n=(nmax+nmin)/2
        kab=I_is(n)
        if (iab-kab) 100,300,200
    100 nmax=n
        goto 1
    200 nmin=n
        goto 1
    300 Find_SMS=is*R_is(n)
        return
    end function Find_SMS

    real(dp) function HintS(na,nb,n0)
        implicit none
        real(dp)    :: e, d, dr, dd, de, d1
        integer   :: nx, ind, nab, i, na, nb, n0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        e=0.d0
        HintS=e
        nx = IPx
        if (Ll(na) > LmaxS)   return
        if (nb > NmaxS)       return
        ind=nx*(na-Nso-1)+(nb-Nso)
        do i=1,NhintS
           nab=Iint1S(i)
           if (nab == ind) then
              e=Rsig(i)
              d=Dsig(i)
              if (Kdsig*d /= 0.d0) then
                de=(E_k+Eps(n0)-Esig(i))
                d1=d/e
                dd=d1*de
                if (dd > 0.5d0) then
                  Kherr=Kherr+1
                  dd=0.5d0
                end if
                if (d1 > 0.d0) then           ! <= NORMAL VARIANT
                  e=e/(1.d0-dd)
                else                           ! <= ANOMALOUS VARIANT
                  select case(Kexn)
                    case(1) ! <= two side extrapolation
                      e=e*(1.d0+dd)
                    case(2)  ! <= one side extrapolation
                      e=e*(1.d0+dmin1(0.d0,dd))
                    case(3)   ! <= nonlenear extrapolation
                      dr=1.d0+dd-0.1*de*de
                      e=e*dmax1(dr,0.d0)
                  end select
                end if
              end if
              HintS=e
              return
           end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        write( 6,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
        write(11,'(/4X,"HintS: integral is absent:"/4X,"na=",I2,2X,"nb=",I2/)') na,nb
        stop
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        HintS=e
        return
    end function HintS

    real(dp) function GintS(i1,i2,i3,i4)
        implicit none
        integer   :: is, ia, ib, ic, id, ma, mb, mc, md, na, nb, nc, nd, &
                     na0, nb0, la, lb, lc, ld, i, ja, jb, jc, jd, k, iab, &
                     iac0, ibd0, kmn, kmx, minint, k1, mi, nx, i1, i2, i3, i4, &
                     ii, iac, ibd, kac, kbd
        real(dp)    :: e, rabcd, dabcd, de, dd, dmin1, dr, dmax1, d1, eabcd
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        e=0.d0
        GintS=e
        nx = IPx
        is=1
        ia=i1
        ib=i2
        ic=i3
        id=i4
        do ii=1,2
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
            if (max(la,lb,lc,ld) > LmaxS) return
            if (max(na,nb,nc,nd) > NmaxS) return
            if ((na+nb+nc+nd) > Nsum)    return
            i=la+lb+lc+ld
            if (i /= 2*(i/2))             return
            ja=Jj(na)
            jb=Jj(nb)
            jc=Jj(nc)
            jd=Jj(nd)
            if (na > nc  .and.  nb >= nc  .and.  nd >= nc) then ! na>nc
                k =na                                     ! nb>=nc
                na=nc                                     ! nd>=nc
                nc=k
                k =nb
                nb=nd
                nd=k
                goto 200
            end if
            if (na > nb  .and.  nc > nb  .and.  nd >= nb) then ! na>nb
                k =na                                     ! nc>nb
                na=nb                                     ! nd>=nb
                nb=k
                k =nc
                nc=nd
                nd=k
                goto 200
            end if
            if (na > nd  .and.  nb > nd  .and.  nc > nd) then ! na>nd
                k =na                                     ! nb>nd
                na=nd                                     ! nc>nd
                nd=k
                k =nb
                nb=nc
                nc=k
                goto 210
            end if
    200     if (na == nb  .and.  nc > nd)then
               k =nc
               nc=nd
               nd=k
            end if
            if (na == nd  .and.  nb > nc)then
                if (na /= nc) then
                    k =nc
                    nc=nb
                    nb=k
                else      ! case: na=nc=nd<nb
                    k =nd
                    nd=nb
                    nb=k
                end if
            end if
            if (na == nc  .and.  nb > nd)then
                k =nb
                nb=nd
                nd=k
            end if
    
    210     if (Ksym == 0) then      ! approx. symmetry
                if (nb > nd) then     ! which is assumed
                    k =nb               ! when Ksym=0             
                    nb=nd                         
                    nd=k
                end if
                if (na == nb  .and.  nc > nd) then
                    k =nc
                    nc=nd
                    nd=k
                end if
            end if
    
            iac0=nx*(na-Nso-1)+(nc-Nso)
            ibd0=nx*(nb-Nso-1)+(nd-Nso)
            kmn=max(iabs(ja-jc)/2+1,iabs(jb-jd)/2+1)
            kmx=min((ja+jc)/2+1,(jb+jd)/2+1,Kmax+1)
            iab = nx*(na-Nso-1)+(nb-Nso)
            minint = IntOrdS(iab)
            do k1=kmn,kmx
                k=k1-1
                iac=nx*nx*k+iac0
                ibd=ibd0
                mi=minint-1
                do i=minint,NgintS
                   kac=Iint2S(i)
                   kbd=Iint3S(i)
                   if (kac == iac  .and.  kbd == ibd) then
                     mi=i
                     goto 230
                   end if
                end do
    220         continue ! this is to control missing integrals
                write(*,*) ' GintS: missing integral'
                write(*,'(4I5,I3)') na,nb,nc,nd,k
                read(*,*)
                goto 240
    230         rabcd=Rint2S(i)
                if (rabcd == 0.d0) then
                  minint = mi + 1
                  cycle
                end if
                dabcd=Dint2S(i)
                if (Kdsig*dabcd /= 0.d0) then
                   eabcd=Eint2S(i)
                   de=E_k+Eps(na0)+Eps(nb0)-eabcd
                   d1=dabcd/rabcd
                   dd=d1*de
                   if (dd > 0.5d0) then
                     Kgerr=Kgerr+1
                     dd=0.5d0
                   end if
                   if (d1 > 0.d0) then           ! <= NORMAL VARIANT
                     rabcd=rabcd/(1.d0-dd)
                   else                           ! <= ANOMALOUS VARIANT
                     if (Kexn == 1) then          ! <= two side extr-n
                       rabcd=rabcd*(1.d0+dd)
                     end if
                     if (Kexn == 2) then          ! <= one side extr-n
                       rabcd=rabcd*(1.d0+dmin1(0.d0,dd))
                     end if
                     if (Kexn == 3) then          ! <= nonlinear extr-n
                       dr=1.d0+dd-0.1*de*de
                       rabcd=rabcd*dmax1(dr,0.d0)
                     end if
                   end if
                end if
                e=e+is &
                    *Gaunt(k,ja*0.5d0,ma*0.5d0,jc*0.5d0,mc*0.5d0) &
                    *Gaunt(k,jd*0.5d0,md*0.5d0,jb*0.5d0,mb*0.5d0) &
                    *rabcd
    240         minint = mi + 1
            end do
            k=ic
            ic=id
            id=k
            is=-is
        end do
        GintS=e
        return
    end function GintS

    real(dp) function Gaunt(k,xj1,xm1,xj2,xm2)
        implicit none
        integer  :: i, is, ind, ib1, ib2, im, k, ij
        real(dp)   :: g, x, xj1, xj2, xm1, xm2, j1, j2, m1, m2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        j1=xj1
        j2=xj2
        m1=xm1
        m2=xm2
        is = 1
        g = 0.d0
        im = abs(m2-m1)
        if (k >= im) then
           if (k == 0) then
              if (j1 == j2) then
                 g = 1.d0
                 gaunt = g*is
                 return
              end if
           else
              if (j2 < j1) then
                 x = j2
                 j2 = j1
                 j1 = x
                 x = m2
                 m2 = m1
                 m1 = x
                 if (im /= 2*(im/2)) is = -is
              end if
              if (m1 <= 0.d0) then
                 m1 = -m1
                 m2 = -m2
                 ij = k+j1-j2
                 if(ij /= 2*(ij/2)) is = -is
              end if
              ib1=2*IPlx+1
              ib2=ib1*ib1
              ind = ib2*(ib2*k+2*(ib1*j1+j2))+ib1*(j1+m1)+(j2+m2)
              do i=1,IPgnt
                 if(In(i) == ind) then
                    g = Gnt(i)
                    gaunt = g*is
                    return
                 end if
             end do
           end if
        end if
        if (K_gnt == 1) then
           write (*,'(1X,"Gaunt: Can not find Gaunt for k=",I2, &
               " j1=",F4.1," m1=",F5.1," j2=",F4.1," m2=",F5.1)') k,xj1,xm1,xj2,xm2
           stop
        end if
        gaunt = g*is
        return
    end function Gaunt

    subroutine Rint
        implicit none
        integer     :: i, nsh, nx, ns1, nso1, nsu1, Nsul, nsp1, Nlist, &
                     k, mlow, m_is
        character*7 :: str(3),str1(4)*4
        data str /'Coulomb','Gaunt  ','Breit  '/
        data str1 /' VS','SMS','NMS',' MS'/
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        nsh=Nso+1
        Nhint=0
        Ngint=0
        ! parameter for indexation of integrals:
        if (Nsu > IPx+nsh-2) then
           write (*,*) ' Rint: IPx is too small'
           stop
        else
           nx = IPx
        end if
        open(unit=13,file='CONF.INT',status='OLD', &
             form='UNFORMATTED',err=700)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,*)' Reading file CONF.INT...'
        read (13) ns1,nso1,nsp1,Nsu1,Ecore
        if (ns1 /=  Ns) then
           write(*,*)' Rint warning: Ns=',Ns,' <> ',Ns1,' push...'
           read(*,*)
        end if
        if (nso1 /=  Nso) then
           write(*,*)' Rint warning: Nso=',Nso,' <> ',Nso1
           stop
        end if
        if (Nsu1 /=  Nsu) then
           write(*,*)' Rint warning: Nsu=',Nsu,' <> ',Nsu1
           if (Nsu1 < Nsu) stop
        end if
        write( 6,'(4X,"Total core energy:",F17.7)') Ecore
        write(11,'(4X,"Total core energy:",F17.7)') Ecore
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        read (13) (Nn(i),Kk(i),Ll(i),Jj(i), i=1,Nsu)
        read (13)
        read (13) Nhint,Kbrt
        allocate(Rint1(Nhint),Iint1(Nhint))
        read (13) (Rint1(i), i=1,Nhint)
        read (13) (Iint1(i), i=1,Nhint)
        read (13) Ngint,Nlist,nrd
        allocate(Rint2(IPbr,Ngint),Iint2(Ngint),Iint3(Ngint),IntOrd(nrd))
        if (nrd /=  nx*nx) then
           write(*,*)' Rint error: IPx was changed'
           stop
        end if
        read (13) ((Rint2(k,i), k=1,IPbr), i=1,Ngint)
        read (13) (Iint2(i), i=1,Ngint)
        read (13) (Iint3(i), i=1,Ngint)
        read (13) (IntOrd(i), i=1,nrd)
        if (K_is >= 1) then
            read(13,end=800) m_is,mlow,num_is
            allocate(R_is(num_is),I_is(num_is))
            if (m_is /=  K_is) then
                write(*,*) 'IS Integrals are for K_is=',m_is
                read(*,*)
                stop
            end if
            if (K_is >= 2   .and.   mlow /=  Klow) then
                write(*,*) 'SMS Integrals are for Klow=',mlow
                read(*,*)
            end if
            read(13,end=800) (R_is(i),i=1,num_is)
            read(13,end=800) (I_is(i),i=1,num_is)
        end if
        close(unit=13)
        write( *,'(4X,A7," integrals read from CONF.INT")') str(Kbrt+1)
        write(11,'(4X,A7," integrals read from CONF.INT")') str(Kbrt+1)
        if (K_is >= 1) then
            write( *,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
            write(11,'(4X,I7," integrals for ",A3," operator found")') num_is,str1(K_is)
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nsp /=  nsp1) then
            write ( 6,*) ' Nsp changed since integrals were calculated'
            write (11,*) ' Nsp changed since integrals were calculated'
            end if
        Return
        !- - - - - - - - - - - - - - - - - - - - - - - - -
    700 write( 6,'(2X," Can not find file CONF.INT...")')
        write(11,'(2X," Can not find file CONF.INT...")')
        Stop
    800 write( *,*) 'No SMS integrals in CONF.INT!'
        Stop
    end subroutine Rint

    Subroutine RintS
    ! Reading of files SGC.CON and SCRC.CON with the self-energy and
    ! screening radial integrals.
    implicit none
    integer :: na, nb, ierr, k, nsh, nx, khot, la, lb, ja, jb, na1, nb1, ind, &
               nso1, khot1, k1, nsx1, nsx2, nav, kbox, i, ns1, idummy
    real :: x, y, z
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        NhintS=0
        NgintS=0
        if (Ksig == 0) return
        nsh=Nso+1
!       parameter for indexation of integrals:
        nx = IPx
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     reading ME of Sigma from SGC.CON
        if (Kdsig /= 0) then
          write(*,'(5X,"Give valence energy: ",$)')
          read(*,*) E_0
          Kexn=1
          write(*,'(5X,"Extrapolation variant ","(1-normal,2-one side,3-nonlin.): ",I2)') Kexn
        end if
        open(unit=18,file='SGC.CON',status='OLD',err=700)
        write(*,*)' reading file SGC.CON ...'
        read(18,'(7X,I3,5X,I1,24X,I1)') NmaxS,LmaxS,khot
        do na=nsh,NmaxS
           la=Ll(na)
           ja=Jj(na)
           if (la.GT.LmaxS) cycle
           do nb=na,NmaxS
              lb=Ll(nb)
              jb=Jj(nb)
              if (lb.NE.la.OR.jb.NE.ja) cycle
              read(18,*) idummy,na1,nb1,x,y,z
              NhintS=NhintS+1
              if (na1.NE.na.OR.nb1.NE.nb) then
                 write(*,*)' wrong indices for ',NhintS,' sigma:'
                 write(*,*)' expected ',na,nb,' got ',na1,nb1
                 stop
              end if
           end do
        end do
        rewind(18)
        allocate(Iint1S(NhintS),Rsig(NhintS),Dsig(NhintS),Esig(NhintS))
        read(18,'(7X,I3,5X,I1,24X,I1)') NmaxS,LmaxS,khot
        i=0
        do na=nsh,NmaxS
           la=Ll(na)
           ja=Jj(na)
           if (la.GT.LmaxS) cycle
           do nb=na,NmaxS
              lb=Ll(nb)
              jb=Jj(nb)
              if (lb.NE.la.OR.jb.NE.ja) cycle
              read(18,*) idummy,na1,nb1,x,y,z
              i=i+1
              if (na1.NE.na.OR.nb1.NE.nb) then
                 write(*,*)' wrong indeces for ',NhintS,' sigma:'
                 write(*,*)' expected ',na,nb,' got ',na1,nb1
                 stop
              end if
              ind=nx*(na-nsh)+(nb-nsh+1)
              Iint1S(i)=ind
              Rsig(i)=x
              Dsig(i)=y
              Esig(i)=z
           end do
        end do
        xscr=10.d0
        if (Ksig == 2) then
           xscr=0.d0
           read(18,'(5X,I2,5X,F8.5)')
           do k=1,10
              read(18,'(5X,I2,5X,F8.5)') k1,Scr(k)
              if (k1.NE.k-1) then
                 write(*,*)' wrong order of screening coefficients'
                 stop
              end if
              xscr=xscr+Scr(k)
           end do
        end if
        write(*,*)' Khot =',khot
        close(unit=18)
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Ksig < 2) return
        ierr=0
        open (unit=13,file='SCRC.CON',status='OLD',form='UNFORMATTED',err=710)
        write(*,*)' Reading file SCRC.CON...'
        read (13) ns1,nso1,nsh,khot1,nsx1,nsx2,Ksym
        read (13) Nmax1,Lmax1,Kmax,nav,kbox
        write(11,'(4X,"Parameters from the file SCRC.CON:", &
               /4X,"Ns=",I3," Nso=",I2," Nsh=",I2," Khot=",I1, &
               " Nsx=",2I3,/4X," Ksym=",I1,"NmaxS=",I3," LmaxS=",I2, &
               " Kmax=",I2," Nav=",I3," Kbox=",I1)') ns1,nso1,nsh, &
                    khot1,nsx1,nsx2,Ksym,Nmax1,Lmax1,Kmax,nav,kbox
        if (nso1 /= Nso) then
           write(*,*)' RintS warning: Nso=',Nso,' <> ',Nso1
           ierr=ierr+1
        end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (NmaxS.NE.nmax1.OR.LmaxS.NE.lmax1.OR.khot.NE.khot1) then
           write( *,'(3X,"Difference between SGC.CON and SCRC.CON:", &
              /3X,"nmaxS =",I3," (SGC) =",I3," (SCRC)", &
              /3X,"lmaxS =",I3," (SGC) =",I3," (SCRC)", &
              /3X,"khot =",I3," (SGC) =",I3," (SCRC)")') & 
                 NmaxS,Nmax1,LmaxS,Lmax1,khot,khot1
           if (NmaxS.LT.Nmax1.OR.LmaxS.LT.Lmax1) ierr=ierr+1
        end if
        read (13) NgintS,nrd
        allocate(Rint2S(NgintS),Iint2S(NgintS),Iint3S(NgintS),Dint2S(NgintS),Eint2S(NgintS),IntOrdS(nrd))
        if (nrd.NE.nx*nx) then
           write(*,*)' RintS: IPx was changed'
           ierr=ierr+1
        end if
        if (ierr.GT.0) then
           write(*,*) ierr,' error(s) detected'
           stop
        end if
        if (nsx2.EQ.0) then ! Nsum is used to skip screening integrals
          Nsum=4*nav
        else
          Nsum=4*Nmax1
        end if
        read (13) (Rint2S(i), i=1,NgintS)
        read (13) (Iint2S(i), i=1,NgintS)
        read (13) (Iint3S(i), i=1,NgintS)
        read (13) (IntOrdS(i), i=1,nrd)
        read (13) (Dint2S(i), i=1,NgintS)
        read (13) (Eint2S(i), i=1,NgintS)
        close(unit=13)
        write(*,*)' File SCRC.CON is read'
        if (dabs(xscr-10.d0).GT.1.d-5) then
          write(*,*) ' Note: Screening coefficients will be used'
          write(*,*) ' above NmaxS ',Nmax1,' and LmaxS',Lmax1
          read(*,*)
        end if
       Return
!     - - - - - - - - - - - - - - - - - - - - - - - - -
 700    write(6, '(2X," Can not find file SGC.CON...")')
        write(11,'(2X," Can not find file SGC.CON...")')
       Stop
 710    write(6, '(2X," Can not find file SCRC.CON...")')
        write(11,'(2X," Can not find file SCRC.CON...")')
       Stop
    end subroutine RintS
end module integrals