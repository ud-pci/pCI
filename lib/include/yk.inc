c     =================================================
c                                 !# Last update 19.05.00
      subroutine Yk(k)            !# This routine goes with integration
      include "hfd.par"           !## routine Sint, which integrates over ro,
       implicit real*8 (a-h,o-z)  !### while Sint1 integrates over r.
       common /ii/ii/kt/kt/h/h0/c/c(IP6)
     1        /r/r(IP6)/v/v(IP6)
        dimension w(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of functions z,y with simpson formula
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of w(i)=r(i)**k
        j=k/2
        k1=k-j*2
        do i=1,ii,ih
          ww=1.d0
          rr=r(i)
          if (k1.EQ.1) ww=rr
          if (j.EQ.0) goto 10
          q1=rr*rr
          ww=ww*q1
          if (j.EQ.1) goto 10
          ww=ww*q1
          if (j.EQ.2) goto 10
          ww=ww*q1
          if (j.EQ.3) goto 10
          ww=ww*q1
 10       w(i)=ww
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        do i=1,i1,ih
           j=k1-i
           if (dabs(c(j)).GE.t) goto 200
           imax=imax-ih
        end do
 200    if (imax.LT.i0) imax=i0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of function z(k;r) by direct integration
        p1=c( 1)*v( 1)*w( 1)
        p2=c(i0)*v(i0)*w(i0)
        gam=c(ii+4)
        cr1=0.d0
        cr2=0.d0
        if (c(1).NE.0.d0) cr1=(c(i0)/c(1)*
     >     (r(1)/r(i0))**gam-1.d0)/((r(i0)-r(1))*(gam+k+2))
        if (c(i0).NE.0.d0) cr2=(c(1)/c(i0)*
     >     (r(i0)/r(1))**gam-1.d0)/((r(1)-r(i0))*(gam+k+2))
        c(1)=c(1)*r(1)/(gam+k+1)*(1.d0-cr1*r(1))
        c(i0)=c(i0)*r(i0)/(gam+k+1)*(1.d0-cr2*r(i0))
        t1=c(i0)*w(i0)
        t2=c(1)*w(1)
        i1=i0+ih
        do i=i1,imax,ih
          q1=c(i)*v(i)
          if (k.NE.0) q1=q1*w(i)
          t=hh*(q1+4.d0*p2+p1)
          t=t2+t
          t2=t1
          t1=t
          p1=p2
          p2=q1
          c(i)=t/w(i)
        end do
        t=c(imax)*w(imax)
        im1=imax+ih
        do i=im1,ii,ih
          p=t
          if (k.NE.0) p=p/w(i)
          c(i)=p
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of w(i)=r(i)**(k+1)
        do i=1,ii,ih
          t=r(i)
          if (k.NE.0) t=t*w(i)
          w(i)=t
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of function y(k;r) by solving dif. equation
c       for f=y(k;r)/r**(k+1)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        k1=k+k+1
        i1=imax+ih
        t1=c(i1)/w(i1)
        p1=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=i1+ih
        t2=c(i1)/w(i1)
        p2=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=imax+1
        do j=1,imax,ih
          i=i1-j
          q1=-c(i)*k1*v(i)/(r(i)*w(i))
          t=hh*(q1+4.d0*p1+p2)
          t=t2-t
          t2=t1
          t1=t
          p2=p1
          p1=q1
          c(i)=t*w(i)*v(i)*hh
        end do
        do i=im1,ii,ih
          c(i)=c(i)*v(i)*hh
        end do
        if (c(ii+4).LE.k+0.d0) then
          c(ii+4)=c(ii+4)+1
        else
          c(ii+4)=k+1
        end if
       return
      end
c     =================================================
c                                 !# Last update 18.05.00
      subroutine Yk2(k,n)         !# This routine goes with integration
      include "hfd.par"           !## routine Sint, which integrates over ro,
       implicit real*8 (a-h,o-z)  !### while Sint1 integrates over r.
       common /ii/ii/kt/kt/z/z/h/h0/c/c(IP6)
     1        /r/r(IP6)/v/v(IP6)
        dimension w(IP6),wn(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if(k.LT.-2.OR.k.GT.9.OR.n.LT.0.OR.n.GT.2) then    !### k in [-2,7]
          write(*,*) ' Y2k: invalid par-s: k=',k,' n=',n  !### n in [0,2]
          stop
        end if
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of functions z,y with simpson formula
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of w(i)=r(i)**k and wn(i)=r(i)**(k+n)
        j=k/2
        k1=k-j*2
        do i=1,ii,ih
          q1=r(i)
          if (k1.EQ.-1) w(i)=1.d0/q1
          if (k1.EQ. 0) w(i)=1.d0
          if (k1.EQ. 1) w(i)=q1
          if (j.EQ.0) goto 100
          q2=q1*q1
          if (j.EQ.-1) q2=1.d0/q2
          w(i)=w(i)*q2
          if (iabs(j).EQ.1) goto 100
          w(i)=w(i)*q2
          if (j.EQ.2) goto 100
          w(i)=w(i)*q2
          if (j.EQ.3) goto 100
          w(i)=w(i)*q2
 100      qn=1.d0
          if (n.GT.0) then
            do l=1,n
              qn=qn*q1
            end do
          end if
          wn(i)=w(i)*qn
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of imax where C(i) becomes negligible
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        do i=1,i1,ih
          j=k1-i
          if (dabs(c(j)).GE.t) goto 200
          imax=imax-ih
        end do
 200    if (imax.LT.i0) imax=i0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of function z(k,n;r) by direct integration:
C       z(k,n;r)=1/(r**k) * int_o^r {r**(k+n)c(r)}dr
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        kn=k+n
        p1=c( 1)*v( 1)*wn( 1)
        p2=c(i0)*v(i0)*wn(i0)
        gam=c(ii+4)
        cr1=0.d0
        cr2=0.d0
        if (c(1).NE.0.d0) cr1=(c(i0)/c(1)*
     >     (r(1)/r(i0))**gam-1.d0)/((r(i0)-r(1))*(gam+kn+2))
        if (c(i0).NE.0.d0) cr2=(c(1)/c(i0)*
     >     (r(i0)/r(1))**gam-1.d0)/((r(1)-r(i0))*(gam+kn+2))
        c(1)=c(1)*r(1)*wn(1)/w(1)/(gam+kn+1)*(1.d0-cr1*r(1))
        c(i0)=c(i0)*r(i0)*wn(i0)/w(i0)/(gam+kn+1)*(1.d0-cr2*r(i0))
        t1=c(i0)*w(i0)
        t2=c(1)*w(1)
        i1=i0+ih
        do i=i1,imax,ih
          q1=c(i)*v(i)
          if (kn.EQ.0) goto 210
          q1=q1*wn(i)
 210      t=hh*(q1+4.d0*p2+p1)
          t=t2+t
          t2=t1
          t1=t
          p1=p2
          p2=q1
          p=t
          c(i)=p/w(i)
        end do
        t=c(imax)*w(imax)
        im1=imax+ih
        do i=im1,ii,ih
          c(i)=t/w(i)
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of w(i)=r(i)**(k+n+1)
        do i=1,ii,ih
          w(i)=r(i)*wn(i)
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of function y(k,n;r) by solving dif. eq.:
c       for f=y(k;r)/r**(k+n+1)
c       f'= - (2k+n+1)/r**(k+n+2) * z(k,n;r)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        k1=k+kn+1
        i1=imax+ih
        t1=c(i1)/w(i1)
        p1=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=i1+ih
        t2=c(i1)/w(i1)
        p2=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=imax+1
        rx=r(imax)
        do l=2,n+1
          rx=rx*r(imax)
        end do
        do j=1,imax,ih
          i=i1-j
          q1=-c(i)*k1*v(i)/(r(i)*w(i))
          t=hh*(q1+4.d0*p1+p2)
          t=t2-t
          t2=t1
          t1=t
          p2=p1
          p1=q1
          c(i)=t*w(i)*v(i)*hh
          rr=r(i)
          do l=2,n+1
            rr=rr*r(i)
          end do
        end do
        do i=im1,ii,ih
          c(i)=c(i)*v(i)*hh
        end do
        if (c(ii+4).LE.k+0.d0) then
          c(ii+4)=c(ii+4)+n+1
        else
          c(ii+4)=kn+1
        end if
       return
      end
c     =================================================
      subroutine Sint(ds)        !# Integration over ro=r*v*hh
      include "hfd.par"          !## which is slightly faster than
       implicit real*8 (a-h,o-z) !### integration over r
       common /ii/ii/kt/kt/h/h0
     >        /c/c(IP6)/r/r(IP6)/v/v(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=c(ii+4)
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=1
        i2=i1+ih
        i3=i2+ih
        r1=r(i1)
        r2=r(i2)
        r3=r(i3)
        t1=0.d0
        t2=0.d0
        if (gam.GT.5.0d0) goto 200
        f1=c(i1)/(r1**gam*hh*v(i1))
        f2=c(i2)/(r2**gam*hh*v(i2))
        f3=c(i3)/(r3**gam*hh*v(i3))
        f21=(f2-f1)/(r2-r1)
        f32=(f3-f2)/(r3-r2)
        f321=(f32-f21)/(r3-r1)
        c2=f321
        c1=f21-c2*(r1+r2)
        c0=f1-c1*r1-c2*r1**2
        g=gam+1.d0
        t2=r1**g*(c0/g+c1/(g+1.d0)*r1+c2/(g+2.d0)*r1**2)
        t1=r2**g*(c0/g+c1/(g+1.d0)*r2+c2/(g+2.d0)*r2**2)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
 200    p1=c(i0)
        p2=c( 1)
        i1=i0+ih
        do 10 i=i1,ii,ih
           q=c(i)
           t=q+4.d0*p1+p2
           t=t2+t
           t2=t1
           t1=t
           p2=p1
           p1=q
 10     continue
        ds=t
       return
      end
c     =================================================
      subroutine Y0(Y)            !# Last update 25.11.06
      include "conf.par"          !## Nq -> Qq
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /nso/nso/ii/ii/kt/kt/h/h0/maxt/maxt/Qq/Qq(IPs)
     >        /c/c(IP6)/r/r(IP6)/v/v(IP6)/P/P(IP6)/Q/Q(IP6)
        dimension y(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     Core density:
        r1=r(1)
        do ni=1,Nso
           qa=Qq(ni)
           call ReadF (12,ni+4,P,Q,2)
           if (ni.EQ.1) gm1=2*P(ii+4)
           gm=2*P(ii+4)
           igm=gm-gm1+0.1d0
           rg=r1**igm
           do i=1,ii,ih
              C(i)=C(i)+qa*(P(i)**2+Q(i)**2)
           end do
           do m=0,MaxT-igm            !### expansion at the origin
             km=ii+5+m+igm            !#### of the electron density
             s=0.d0
             do i=0,m
               ki=ii+5+i
               ki1=ii+5+m-i
               s=s+(P(ki)*P(ki1)+Q(ki)*Q(ki1))
             end do
             C(km)=C(km)+s*qa*rg
           end do
        end do
        C(ii+4)=gm1
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        do i=1,i1,ih
           j=k1-i
           if (dabs(c(j)).GE.t) goto 200
           imax=imax-ih
        end do
 200    if (imax.LT.i0) imax=i0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of z(0;r)=int_o^r(ro*dr) by direct integration
        r1=R(1)
        ri=R(i0)
        gm=C(ii+4)
        rg1=r1**(gm+1)
        rgi=ri**(gm+1)
        ri1=ri/r1
        t2=0.d0
        t1=0.d0
        fr=1.d0
        do m=0,maxt              !### t1=int_0^r1 (ro*dr)
          im=ii+5+m              !### t2=int_0^r(i0) (ro*dr)
          C(im)=C(im)/(gm+m+1)
          t2=t2+C(im)*rg1
          t1=t1+C(im)*fr*rgi
          fr=fr*ri1
        end do

        p1=c( 1)*v( 1)
        p2=c(i0)*v(i0)
        c(1)=t2
        c(i0)=t1
        i1=i0+ih
        do i=i1,imax,ih
          q1=c(i)*v(i)
          t=hh*(q1+4.d0*p2+p1)
          t=t2+t
          t2=t1
          t1=t
          p1=p2
          p2=q1
          c(i)=t
        end do
        t=c(imax)
        im1=imax+ih
        do i=im1,ii,ih
          c(i)=t
        end do
        C(ii+4)=gm1+1
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c       evaluation of function y(0;r) by solving dif. equation
c       for f=y(0;r)/r
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=imax+ih
        t1=c(i1)/r(i1)
        p1=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=i1+ih
        t2=c(i1)/r(i1)
        p2=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=imax+1
        do j=1,imax,ih
          i=i1-j
          q1=-c(i)*v(i)/(r(i)*r(i))
          t=hh*(q1+4.d0*p1+p2)
          t=t2-t
          t2=t1
          t1=t
          p2=p1
          p1=q1
          y(i)=t
        end do
        do i=im1,ii,ih
          y(i)=c(i)/r(i)
        end do
ccc        write(*,*)' Y0: C(ii)=',C(ii)
ccc        read(*,*)

        y(ii+4)=0.d0                       !### Y0 -> const as r->0
        rg=r1**gm1
        s=0.d0
        do m=0,MaxT
          s=s+rg*C(ii+5+m)/(gm1+m)
        end do
        Y(ii+5)=Y(1)+s
        if (dabs(gm1-2.d0).LT.1.d-5) then  !### the case of finite nucleus
          Y(ii+6)=0.d0
          do m=2,MaxT
            Y(ii+5+m)=-rg*C(ii+3+m)/m
          end do
        else                               !### the case of point nucleus
          do m=1,MaxT
            Y(ii+5+m)=0.d0
          end do
        end if

       return
      end
