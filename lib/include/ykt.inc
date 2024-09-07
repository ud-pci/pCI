c       =================================================
        subroutine ykt(k,nmax,ro,c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /KT/KT/ii/ii/h/h0/r2/r2
        common /r/r(IP6)/v/v(IP6)
        real*8, dimension(IP6):: ro,c,w
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IH=2-KT
        H=H0*IH
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation functions z(k,r) and y(k,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk1=1.d0-dexp(-k*h)
        dk2=1.d0-dexp(-(k+1)*h)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (r2.lt.0.d0) goto 200
        i1=i0+1
        do i=i1,ii
          w(i)=0.d0
          if (k.ne.0) w(i)=1.d0-(r(i-1)/r(i))**k
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     g=ro(ii+4)+k
        t0=0.d0
        p0=0.d0
        do m=0,nmax
          i=ii+5+m
          t0=t0+ro(i)/(g+m+1)
          p0=p0+ro(i)*(g+m)
        enddo
        t0=t0*r0**(g+1-k)
        p0=p0*r0**(g-1-k)
        f0=ro(i0)
        dvr=(-3*v(1)/r(1)+4*v(2)/r(2)-v(3)/r(3))/(2*h)
        dv=r(1)/v(1)*dvr+v(1)/r(1)
        p0=h*p0*v0*v0+f0*dv
c       p0=h*p0*v0*v0+f0*bt/(al*r0+bt)**2
        fm=ro(imax)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the function z(k;r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        t=t0+0.5d0*f0+h/12.d0*p0
        c(i0)=t
        d=dk1
        i1=i0+1
        do i=i1,imax
          fi=ro(i)
          if (r2.gt.0.d0) d=w(i)
          if (k.ne.0) fi=fi-d*t
          t=t+fi
          c(i)=t
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the function y(k;r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        t=-0.5d0*fm
        dh=(2*k+1)*h/12.d0
        d=dk2
        im=imax-1
        i1=im+i0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 j=i0,im
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=i1-j
        ip=i+1
        s=dh
        if (r2.lt.0.d0) goto 210
        d=r(i)/r(ip)
        if (k.ne.0) d=d*(1.d0-w(ip))
        d=1.d0-d
        s=dh*v(i)/r(i)
210     fi=ro(ip)
        fi=fi-d*(t+fi)
        t=t+fi
        c(i)=c(i)+t-s*ro(i)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
40      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=dk1
        t=c(imax)-0.5d0*fm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 50 i=imax,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c(i)=t
        if (k.eq.0) goto 50
        if (r2.gt.0.d0) d=w(i)
        t=t-d*t
c       - - - - - - - - - - - - - - - - - - - - - - - - -
50      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g0=ro(ii+4)
        do m=0,nmax
          i=ii+5+m
          c(i)=0.d0
        enddo
        c(ii+3)=ii
        c(ii+4)=1.d0
        id=g0+0.5d0
        r0d=r0**id
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (g0-k.gt.1.d-4) then
          t0=0.d0
          do m=0,nmax
            i=ii+5+m
            t0=t0+(ro(i)/(g0+m-k)-ro(i)/(g0+m+k+1))
          enddo
          t0=t0*r0**g0
          if (k.lt.nmax) c(ii+5+k)=c(i0)/r0+t0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          do m=0,nmax
            i=ii+5+m
            if ((m+id).le.nmax) c(i+id)=
     1      ro(i)*(1.d0/(g0+k+m+1)-1.d0/(g0+m-k))*r0d
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
          c(ii+5+id)=c(1)/r(1)**(id+1)*r0d
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine ukt(k,nmax,ro,c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /KT/KT/ii/ii/h/h0/r2/r2
        common /r/r(IP6)/v/v(IP6)
        real*8, dimension(IP6):: ro,c,w
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IH=2-KT
        H=H0*IH
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation the function uk(k,r)=z(k+2,r)-z(k,r).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk1=1.d0-dexp(-k*h)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (r2.lt.0.d0) goto 200
        i1=i0+1
        do 10 i=i1,ii
        w(i)=0.d0
        if (k.ne.0) w(i)=1.d0-(r(i-1)/r(i))**k
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     g=ro(ii+4)+k
        t0=0.d0
        p0=0.d0
        t2=0.d0
        p2=0.d0
        do 20 m=0,nmax
        i=ii+5+m
        t0=t0+ro(i)/(g+m+1)
        t2=t2+ro(i)/(g+m+3)
        p0=p0+ro(i)*(g+m)
        p2=p2+ro(i)*(g+m+2)
20      continue
        t0=t0*r0**(g+1-k)
        t2=t2*r0**(g+1-k)
        p0=p0*r0**(g-1-k)
        p2=p2*r0**(g-1-k)
        f0=ro(i0)
        dvr=(-3*v(1)/r(1)+4*v(2)/r(2)-v(3)/r(3))/(2*h)
        dv=r(1)/v(1)*dvr+v(1)/r(1)
        p0=h*p0*v0*v0+f0*dv
        p2=h*p2*v0*v0+f0*dv
c        p0=h*p0*v0*v0+f0*bt/(al*r0+bt)**2
c        p2=h*p2*v0*v0+f0*bt/(al*r0+bt)**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the function uk(k;r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        zi=t0+0.5d0*f0+h/12.d0*p0
        z1=(t2+0.5d0*f0+h/12.d0*p2)-zi
        c(i0)=z1
        d=dk1
        i1=i0+1
        do 30 i=i1,imax
        if (r2.gt.0.d0) d=w(i)
c
        d12=r(i-1)/r(i)
        d12=d12*d12
        d2=d*d12+1.d0-d12
        z1=z1-d2*z1+(d-d2)*zi
c
        fi=ro(i)
        if (k.ne.0) fi=fi-d*zi
        zi=zi+fi
c
        c(i)=z1
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=dk1
        fm=ro(imax)
        t=c(imax)-0.5d0*fm
        do 50 i=imax,ii
        c(i)=t
        if (k.eq.0) goto 50
        if (r2.gt.0.d0) d=w(i)
        t=t-d*t
50      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Correction
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        zi=h/12.d0*p0
c        c(i0)=c(i0)-h/12.d0*p2+zi
        dh1=2*h/12.d0
        do 40 i=i0,imax
        c(i)=c(i)-dh1*ro(i)*v(i)/r(i)
40      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g0=ro(ii+4)
        id=g0+0.5d0
        r0d=r0**id
        do 70 m=0,nmax
        i=ii+5+m
        c(i)=0.d0
70      continue
        c(ii+3)=ii
        c(ii+4)=1.d0+g0-id
        do 80 m=0,nmax
        i=ii+5+m
        if (m+id.gt.nmax) goto 80
        c(i+id)=ro(i)*(1.d0/(g0+k+2+m+1)-1.d0/(g0+k+m+1))*r0d
80      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine vkt(k,nmax,ro,c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'hfd.par'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /KT/KT/ii/ii/h/h0/r2/r2
        common /r/r(IP6)/v/v(IP6)
        real*8, dimension(IP6):: ro,c,w
        real*8 v1(10),v2(10)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IH=2-KT
        H=H0*IH
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       calculation functions v(k,r)=x(k+2,r)-x(k,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk2=1.d0-dexp(-(k+1)*h)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (r2.lt.0.d0) goto 200
        i1=i0+1
        do 10 i=i1,ii
        d=r(i-1)/r(i)
        if (k.eq.0) w(i)=1.d0-d
        if (k.ne.0) w(i)=1.d0-d**(k+1)
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     g=ro(ii+4)+k
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the function v(k;r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fm=ro(imax)
        xi=-0.5d0*fm
        x1=0.d0
        d=dk2
        im=imax-1
        i1=im+i0
        do 30 j=i0,im
        i=i1-j
        ip=i+1
        if (r2.gt.0.d0) d=w(ip)
c
        d12=r(ip-1)/r(ip)
        d12=d12*d12
        d2=d*d12+1.d0-d12
        x1=x1-d2*x1+(d-d2)*(xi+ro(ip))
c
        fi=ro(ip)
        fi=fi-d*(xi+fi)
        xi=xi+fi
c
        c(i)=x1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (i-i0+1.le.10) then
        v1(i)=xi
        v2(i)=c(i)+xi
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Correction
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dh1=2*h/12.d0
        do 40 i=i0,imax
        c(i)=c(i)-dh1*ro(i)*v(i)/r(i)
40      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 50 i=imax+1,ii
        c(i)=0.d0
50      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g0=ro(ii+4)
        do m=0,nmax
        i=ii+5+m
        c(i)=0.d0
        enddo
        c(ii+3)=ii
        id=g0+0.5d0
        r0d=r0**id
        c(ii+4)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Near the origin
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (g0-k.gt.1.d-4) then
          t0=0.d0
          do m=0,nmax
          i=ii+5+m
          t0=t0+ro(i)/(g0+m-k)
          enddo
          t0=t0*r0**g0
          if (k.le.nmax) c(ii+5+k)=c(ii+5+k)-v1(i0)/r0-t0
          do m=0,nmax
          i=ii+5+m
          if ((m+id).le.nmax)
     &    c(i+id)=c(i+id)+ro(i)/(g0+m-k)*r0d
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          id=g0+0.5d0-1
          r0d=r0**id
          c1=v1(1)/r(1)**(id+1)
          c2=v1(3)/r(3)**(id+1)
          c3=v1(5)/r(5)**(id+1)
          x1=r(1)
          x2=r(3)
          x3=r(5)
          a2=(c2-c1)/(x2-x1)
          a1=c1-a2*x1
          if (id.le.nmax)   c(ii+5+id)=a1*r0d
          if (id+1.le.nmax) c(ii+6+id)=a2*r0d*r0
          a3=(x3-x2)*(c2-c1)-(x2-x1)*(c3-c2)
          a3=a3/((x2*x2-x1*x1)*(x3-x2)-(x3*x3-x2*x2)*(x2-x1))
          a2=(c2-c1)-a3*(x2*x2-x1*x1)
          a2=a2/(x2-x1)
          a1=c1-a2*x1-a3*x1*x1
          if (id.le.nmax)   c(ii+5+id)=c(ii+5+id)-a1*r0d
          if (id+1.le.nmax) c(ii+6+id)=c(ii+6+id)-a2*r0d*r0
          if (id+2.le.nmax) c(ii+7+id)=c(ii+7+id)-a3*r0d*r0**2
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (g0-k-2.gt.1.d-4) then
          t2=0.d0
          do m=0,nmax
          i=ii+5+m
          t2=t2+ro(i)/(g0+m-k-2)
          enddo
          t2=t2*r0**g0
          if (k+2.le.nmax)c(ii+5+k+2)=c(ii+5+k+2)+v2(i0)/r0+t2
          do m=0,nmax
          i=ii+5+m
          if ((m+id).le.nmax)
     &    c(i+id)=-ro(i)/(g0+m-k-2)*r0d+ro(i)/(g0+m-k)*r0d
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
ccccccccccccccccc
          id=g0+0.5d0-1
          r0d=r0**id
ccccccccccccccccc
          c1=v2(1)/r(1)**(id+1)
          c2=v2(3)/r(3)**(id+1)
          c3=v2(5)/r(5)**(id+1)
          x1=r(1)
          x2=r(3)
          x3=r(5)
          a2=(c2-c1)/(x2-x1)
          a1=c1-a2*x1
          if (id.le.nmax)   c(ii+5+id)=a1*r0d
          if (id+1.le.nmax) c(ii+6+id)=a2*r0d*r0
c          c(ii+5+id)=c1*r0d
          a3=(x3-x2)*(c2-c1)-(x2-x1)*(c3-c2)
          a3=a3/((x2*x2-x1*x1)*(x3-x2)-(x3*x3-x2*x2)*(x2-x1))
          a2=(c2-c1)-a3*(x2*x2-x1*x1)
          a2=a2/(x2-x1)
          a1=c1-a2*x1-a3*x1*x1
          if (id.le.nmax)   c(ii+5+id)=c(ii+5+id)+a1*r0d
          if (id+1.le.nmax) c(ii+6+id)=c(ii+6+id)+a2*r0d*r0
          if (id+2.le.nmax) c(ii+7+id)=c(ii+7+id)+a3*r0d*r0**2
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
