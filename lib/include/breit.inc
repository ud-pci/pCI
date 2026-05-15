c     ============= 15/11/99 ==== last update: 26/11/99
c     subroutines for breit interaction calculation
c     =================================================
      Subroutine Mult(xk1,   !### min multipolarity
     >                xk2,   !### max ultipolarity
     >                ip)    !### parity
       implicit real*8 (a-h,o-z)
       common /kout/kout/small/small
       Common /XJ/xja,xjb,xjc,xjd !### total angular momenta of orbitals
     >        /XL/xla,xlb,xlc,xld !### orbital angular momenta of upper c-s

        ip=2
        xk1 = dmax1(dabs(xja-xjc),dabs(xjb-xjd))
        xk2 = dmin1(xja+xjc,xjb+xjd)
        is = xla+xlc+1+small
        is1=Isig(is)
        is = xlb+xld+1+small
        is2=Isig(is)
        if (is1.NE.is2) then
          write(*,5) xja,xjb,xjc,xjd,xla,xlb,xlc,xld
 5        format(' Mult error: different parity of Bra and Ket states',
     >           /' ja-jc:',4F5.1,/' la-ld:',4F5.1)
          stop
        end if
        ip=is1
        if (kout.GT.2) then
          write(11,*)' xla=',xla,' xlc=',xlc
          write(11,15) xk1,xk2,is1
 15       format(4X,'j_min =',F4.1,' j_max =',F4.1,' parity =',I2)
        end if
       return
      end
c     =================================================
      Function Br_core(na,nb, !### exchange breit core potential
     >                 kfl)   !### file with derivatives
      INCLUDE "conf.par"      !### see 12/11/99 and 8/05/00
      include "hfd.par"
c                             !#### for derivation details
c     Exchange core potential:
c     <b|V1_B+V2_B|a> = sum_c sum_k Gaunt^k_{a,c,c,b} R^k_{a,c,c,b}
c
       implicit real*8 (a-h,o-z)
       common /Nso/Nso/Ns/Ns/Kbrt/Kbrt/Kt/Kt/Kt1/Kt1
     >        /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)/Qq/Qq(IPs)
     >        /UP/Pa(IP6),Pc(IP6),Pd(IP6),Pb(IP6)
     >        /DOWN/Qa(IP6),Qc(IP6),Qd(IP6),Qb(IP6)
     >        /UP1/P1a(IP6),P1c(IP6),P1d(IP6),P1b(IP6)  !### derivatives
     >        /DOWN1/Q1a(IP6),Q1c(IP6),Q1d(IP6),Q1b(IP6)
     >        /XJ/xja,xjc,xjd,xjb !### total angular moments of orbitals
     >        /XL/xla,xlc,xld,xlb !### orbital angular moments of upper c-s
     >        /YL/yla,ylc,yld,ylb !### orbital angular moments of lower c-s
       common /kout/kout/small/small
       common /Ii/Ii
       dimension Rint1(10),      !### Radial integrals after regrouping
     >           Rint2(10)
	common /let/let(9)
       character*1 let
!       character*1 let(9)
!        data let/'s','p','d','f','g','h','i','k','l'/
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        z12=0.5d0
        z00=0.0d0
        Br_core=z00
        if (Nso.EQ.0) return
!
        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)
        if (ja.NE.jb.OR.la.NE.lb) return
!
        xja=0.5d0*ja
!
        a_ab=0.d0
        do nc=1,Nso
          call ReadF(12,nc+4,Pc,Qc,2)
          jc=Jj(nc)
          lc=Ll(nc)
          xjc=0.5d0*jc
!--
! Calculation of both Breit contributions using IIT function
!--
          kmin=abs(xja-xjc)
          kmax=xja+xjc
!
	  s = 0.d0
          do k=kmin,kmax
!***        gac=cg(1,k,na,nc)*qq(nc)
            xj=k
            gac=-(2*xjc+1)*FJ3(xja,xjc,xj,z12,-z12,z00)**2
	    ds = 0.d0
            if (dabs(gac).gt.1.d-7) then
              ds = breit_int(k,na,Pa,Qa,nc,Pc,Qc,nc,Pc,Qc,nb,Pb,Qb)
            endif
	    s = s+ds*gac
          enddo
!--
          if (kout.GT.2) then
            write( *,5) Nn(nc),let(lc+1),jc,s
            write(11,5) Nn(nc),let(lc+1),jc,s
 5          format(2X,'contribution of ',I2,A1,I1,'/2 :',E12.5)
          end if
          a_ab=a_ab+s
!
        end do
!--
        if (kout.GT.2) then
          write( *,15) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
          write(11,15) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
 15       format(/2X,'<',I2,A1,I1,'/2| V_B^core |',
     >           I2,A1,I1,'/2> = ',E12.5)
        end if
        Br_core=a_ab
       return
      end
c       =================================================
        function breit_int(l,na,pa,qa,nb,pb,qb,nc,pc,qc,
     &  nd,pd,qd)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /kbrt/kbrt
        common /Ll/Ll(IPs)
        real*8, dimension(IP6) :: pa,qa,pb,qb,pc,qc,pd,qd
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        breit_int=0.d0
        la=Ll(na)
        lb=Ll(nb)
        lc=Ll(nc)
        ld=Ll(nd)
        i=la+lc+lb+ld
        if (2*(i/2).ne.i) return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call breit_pot_magn(l,na,pa,qa,nc,pc,qc)
        ds1=brint_magn(l,nb,pb,qb,nd,pd,qd)
        ds2=0.d0
        if (kbrt.eq.2) then
          call breit_pot_ret(l,na,pa,qa,nc,pc,qc)
          ds2=brint_ret(l,nb,pb,qb,nd,pd,qd)
        endif
        breit_int=ds1+ds2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine breit_pot_magn(l,na,pa,qa,nc,pc,qc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	parameter (maxl=14,maxk=2*maxl+2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Ll/Ll(IPs)
        common /yy/yy(IP6,2*(maxk+1))
        real*8 pa(IP6),qa(IP6),pc(IP6),qc(IP6)
	real*8 ro(IP6),c(IP6)
        common /MaxT/MaxT
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        la=Ll(na)
        lc=Ll(nc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 k=l-1,l+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k.lt.0) goto 10
        if (l.eq.0.and.k.eq.0) goto 10
        i=la+lc+k
        if (2*(i/2).eq.i) goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k.gt.maxk) then
        write( k7,'(/2x,a/2x,a)')
     &  '*** k.gt.maxk in brin_m_fast ***', 'Increase maxk'
!        call exit1
         call exit(1)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g1=coefb( 1,k,l,na,nc)
        g2=coefb(-1,k,l,na,nc)
!
        call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
        call ykt(k,MaxT,ro,c)
!
        do i=1,IP6
          yy(i,k+1)=c(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine breit_pot_ret(l,na,pa,qa,nc,pc,qc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	parameter (maxl=14,maxk=2*maxl+2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Ll/Ll(IPs)
        common /yy/yy(IP6,2*(maxk+1))
        common /uu/uu(IP6,2*(maxk+1))
        common /vv/vv(IP6,2*(maxk+1))
        real*8 pa(IP6),qa(IP6),pc(IP6),qc(IP6)
	real*8 ro(IP6),c(IP6)
        common /MaxT/MaxT
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        la=Ll(na)
        lc=Ll(nc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 k1=l-1,l+1,2
        if (k1.lt.0) goto 40
        i=la+lc+k1
        if (2*(i/2).eq.i) goto 40
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g1=coefb( 1,k1,l,na,nc)
        g2=coefb(-1,k1,l,na,nc)
        call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 50 k2=l-1,l+1,2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k2.lt.0) goto 50
        i=la+lc+k2
        if (2*(i/2).eq.i) goto 50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk=k2
        if (k2.eq.l+1) then
          ck2=dsqrt(dk/(2*l+1))
        else
          ck2=-dsqrt((dk+1.d0)/(2*l+1))
        endif
        if (ck2.eq.0.d0) goto 50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g1=coefb( 1,k1,l,na,nc)
        g2=coefb(-1,k1,l,na,nc)
        call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k1.eq.k2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (k1.gt.maxk) then
            write( k7,'(/2x,a/2x,a)')
     &      '*** k1.gt.maxk in brin_r_fast ***', 'Increase maxk'
!            call exit1
            call exit(1)
          endif
!
          call ykt(k1,MaxT,ro,c)
!
          do i=1,IP6
            yy(i,k1+1)=c(i)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (k2.eq.k1+2) k=k1
          if (k1.eq.k2+2) k=k2
          if (k.gt.maxk) then
            write( k7,'(/2x,a/2x,a)')
     &      '*** k.gt.maxk in brin_r_fast ***', 'Increase maxk'
!            call exit1
            call exit(1)
          endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (k2.eq.k1+2) then
!
            call ukt(k,MaxT,ro,c)
!
            do i=1,IP6
              uu(i,k+1)=c(i)
            enddo
          endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (k1.eq.k2+2) then
!
           call vkt(k,MaxT,ro,c)
!
           do i=1,IP6
              vv(i,k+1)=c(i)
           enddo
         endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
50      continue
40      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        function brint_magn(l,nb,pb,qb,nd,pd,qd)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	parameter (maxl=14,maxk=2*maxl+2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Ll/Ll(IPs)
        common /yy/yy(IP6,2*(maxk+1))
        real*8, dimension(IP6) :: pb,qb,pd,qd
	real*8, dimension(IP6) :: c, ro
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        lb=Ll(nb)
        ld=Ll(nd)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s1=0.d0
        s2=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 k=l-1,l+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k.lt.0) goto 10
        if (l.eq.0.and.k.eq.0) goto 10
        i=lb+ld+k
        if (2*(i/2).eq.i) goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,IP6
          c(i)=yy(i,k+1)
        enddo
c
        call sbint(pb,qd,c,ds1)
        s1=s1+coefb(-1,k,l,nd,nb)*ds1
c
        call sbint(qb,pd,c,ds2)
        s2=s2+coefb( 1,k,l,nd,nb)*ds2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        brint_magn=-2*(s1+s2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        function brint_ret(l,nb,pb,qb,nd,pd,qd)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	parameter (maxl=14,maxk=2*maxl+2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Ll/Ll(IPs)
        common /yy/yy(IP6,2*(maxk+1))
        common /uu/uu(IP6,2*(maxk+1))
        common /vv/vv(IP6,2*(maxk+1))
        real*8 pb(IP6),qb(IP6),pd(IP6),qd(IP6)
	real*8 c(IP6)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        lb=Ll(nb)
        ld=Ll(nd)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s1=0.d0
        s2=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 k1=l-1,l+1,2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k1.lt.0) goto 40
        i=lb+ld+k1
        if (2*(i/2).eq.i) goto 40
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk=k1
        if (k1.eq.l+1) then
          ck1=dsqrt(dk/(2*l+1))
        else
          ck1=-dsqrt((dk+1.d0)/(2*l+1))
        endif
        if (dabs(ck1).lt.1.d-10) goto 40
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 50 k2=l-1,l+1,2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k2.lt.0) goto 50
        i=lb+ld+k2
        if (2*(i/2).eq.i) goto 50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk=k2
        if (k2.eq.l+1) then
          ck2=dsqrt(dk/(2*l+1))
        else
          ck2=-dsqrt((dk+1.d0)/(2*l+1))
        endif
        if (ck2.eq.0.d0) goto 50
c
        c12=ck1*ck2*dsqrt((2*k1+1.d0)*(2*k2+1.d0))
        if (k1.eq.k2) c12=2.d0*c12/(2*k1+1)
        if (dabs(c12).lt.1.d-10) goto 50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k1.eq.k2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          do i=1,IP6
            c(i)=yy(i,k1+1)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
         if (k2.eq.k1+2) k=k1
         if (k1.eq.k2+2) k=k2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
         if (k2.eq.k1+2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          do i=1,IP6
            c(i)=uu(i,k+1)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
         endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
         if (k1.eq.k2+2) then
             do i=1,IP6
               c(i)=vv(i,k+1)
             enddo
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call sbint(pb,qd,c,ds1)
        s1=s1+c12*coefb( 1,k2,l,nb,nd)*ds1
c
        call sbint(qb,pd,c,ds2)
        s2=s2+c12*coefb(-1,k2,l,nb,nd)*ds2
50      continue
40      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        brint_ret=-(s1+s2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine rho_pq_tot(coef1,coef2,p,q,a,b,ro)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
        include 'conf.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /KT/KT/ii/ii/h/h0
	common /r/r(IP6)/v/v(IP6) 
        common /MaxT/MaxT !/ulam/ulam
        real*8 p(IP6),q(IP6),a(IP6),b(IP6)
        real*8 ro(IP6)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IH=2-KT
        H=H0*IH
	ulam = 1
!--
!       imax=p(ii+3)+c1
!       if (imax > a(ii+3)) imax = a(ii+3)
!--
        imax=ii          ! mgk
        dh1=ulam*h*coef1
        dh2=ulam*h*coef2
        do i=1,imax
          ro(i)=(dh1*p(i)*b(i)+dh2*q(i)*a(i))*v(i)
        enddo
        do i=imax+1,ii
          ro(i)=0.d0
        enddo
        ro(ii+3)=imax
        ro(ii+4)=p(ii+4)+a(ii+4)
!
        d1=ulam*coef1
        d2=ulam*coef2
        do m=0,MaxT
          i=ii+5+m
          dr=0.d0
          do n=0,m
            j=ii+5+n
            dr=dr+(d1*p(j)*b(i-n)+d2*q(j)*a(i-n))
          enddo
          ro(i)=ulam*dr
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        function coefb(ibet,k,l,na,nb)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'conf.par'	
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Kk/Kk(IPs)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        u=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ka=Kk(na)
        kb=Kk(nb)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        alk=0.5d0*(k*(k+1)-l*(l+1))
        lk=k+l
        blk=1.d0
        if (k.eq.l) blk=1.d0/dsqrt(2.d0)
        u=dsqrt((2*l+1.d0)/(lk*(lk+1)*(lk+2)))/blk
        d=kb
        if (2*(lk/2).ne.lk) d=-d
        u=u*(ka+d+ibet*alk)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        coefb=u
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine sbint(a,b,c,s)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/h/h
     	common /r/r(IP6)/v/v(IP6)
        common /MaxT/MaxT
        real*8, dimension(IP6) :: a,b,c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        r0=r(i0)
        v0=v(i0)
	s=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
!        dmax=a(ii+3)
!        if (b(ii+3).lt.dmax) dmax=b(ii+3)
!        imax=dmax+0.01d0
!
        imax=ii   ! mgk 
        g=a(ii+4)+b(ii+4)+c(ii+4)-1.d0
        t0=0.d0
        p0=0.d0
        do 10 m=0,MaxT
        i=ii+5+m
        dw=0.d0
        do 20 n=0,m
        j=ii+5+n
        do k=0,m-n
        j1=ii+5+k
          dw=dw+a(j)*b(j1)*c(i-n-k)
        enddo
20      continue
        t0=t0+dw/(g+m+1)
        p0=p0+dw*(g+m)
10      continue
        t0=t0*r0**(g+1)
        p0=p0*r0**(g-1)
        f0=c(i0)*(a(i0)*b(i0))*v0/r(i0)
        dvr=(-3*v(1)/r(1)+4*v(2)/r(2)-v(3)/r(3))/(2*h)
        dv=r(1)/v(1)*dvr+v(1)/r(1)
        p0=p0*v0*v0+f0*dv
c       p0=p0*v0*v0+f0*bt/(al*r0+bt)**2
!	f0 = 0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     dt=0.d0
        do i=i0,imax
          f=c(i)*(a(i)*b(i))
          f=f*v(i)/r(i)
          dt=dt+f
        enddo
        s=t0+h*dt-h*(0.5d0*f0-h/12.d0*p0)-0.5d0*h*f
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        function cg(i,k,na,nb)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        INCLUDE "conf.par"
        INCLUDE "hfd.par"
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Ll/Ll(IPs)/Jj/Jj(IPs)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        la=Ll(na)
        ja=jj(na)
        lb=Ll(nb)
        jb=jj(nb)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (k.gt.0.or.i.gt.0) goto 200
        g=1.d0
        goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     if (na.ne.nb.and.i.lt.0) goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dja=0.5d0*ja
        djb=0.5d0*jb
        dk=k
!        chab=ch(dja,0.5d0,djb,-0.5d0,dk,0.d0)
!        y=2*chab**2/(2*k+1)
	y = 2.d0*Fj3(dja,djb,dk,0.5d0,-0.5d0,0.d0)**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (na.ne.nb) d=1.d0
        nd=ja+1
        if (na.eq.nb) d=nd/(nd-1.d0)
        g=-0.5d0*d*y
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    if (dabs(g).lt.1.d-7) g=0.d0
        cg=g
        return
        end
