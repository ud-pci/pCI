Module sintg

    Use basc_variables

    Implicit None

    Private

    Public :: Sint, Sint1, Y0, Yk

  Contains

    Subroutine Sint(ds) 
        !# Integration over ro=r*v*hh which is slightly faster than integration over r
        Implicit None

        Integer :: i, ih, i0, i1, i2, i3
        Real(dp) :: hh, r1, r2, r3, t1, t2, f1, f2, f3, f21, f32, f321, c0, c1, c2, g, p1, p2, q, ds
        Real(dp) :: gam, t

        gam=c(ii+4)
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0

        i1=1
        i2=i1+ih
        i3=i2+ih
        r1=r(i1)
        r2=r(i2)
        r3=r(i3)
        t1=0.d0
        t2=0.d0
        If (gam.LE.5.0d0) Then
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
        End If
        p1=c(i0)
        p2=c( 1)
        i1=i0+ih
        Do i=i1,ii,ih
           q=c(i)
           t=q+4.d0*p1+p2
           t=t2+t
           t2=t1
           t1=t
           p2=p1
           p1=q
        End Do

        ds=t
        Return
    End Subroutine Sint

    Subroutine Sint1(DS)        
        ! Simpson integration over r (with weight function HH*V(I))
        Implicit None

        Integer :: I, IH, I0, I1, I2, I3
        Real(dp) :: R1, R2, R3, T1, T2, F1, F2, F3, F21, F32, F321, &
                    C0, C1, C2, G, P1, P2, T, Q, HH, Gam
        Real(dp), intent(out) :: DS

        Gam=C(ii+4)
        IH=2-KT
        HH=H*IH/3.d0
        I0=IH+1

        I1=1
        I2=I1+IH
        I3=I2+IH
        R1=R(I1)
        R2=R(I2)
        R3=R(I3)
        T1=0.d0
        T2=0.d0
        If (GAM <= 5.0d0) Then
            F1=C(I1)/R1**GAM
            F2=C(I2)/R2**GAM
            F3=C(I3)/R3**GAM
            F21=(F2-F1)/(R2-R1)
            F32=(F3-F2)/(R3-R2)
            F321=(F32-F21)/(R3-R1)
            C2=F321
            C1=F21-C2*(R1+R2)
            C0=F1-C1*R1-C2*R1**2
            G=GAM+1.d0
            T2=R1**G*(C0/G+C1/(G+1.d0)*R1+C2/(G+2.d0)*R1**2)
            T1=R2**G*(C0/G+C1/(G+1.d0)*R2+C2/(G+2.d0)*R2**2)
        End If
        Dint=t2
        P1=C(I0)*V(I0)
        P2=C( 1)*V( 1)
        I1=I0+IH
        Do I=I1,II,IH
           Q=C(I)*V(I)
           T=HH*(Q+4.d0*P1+P2)
           T=T2+T
           T2=T1
           T1=T
           P2=P1
           P1=Q
        End Do
        
        ds=t
        Return
    End Subroutine Sint1

    Subroutine Y0(Y)            !# Last update 25.11.06
        Use readfff
        Implicit None

        Integer :: i0, ih, ni, igm, i, km, ki, kil, imax, i1, k1, j, ki1, im, im1
        Real(dp) :: hh, r1, qa, gm, rg, s, gm1, t, t1, t2, fr, p1, p2, rg1, ri1, q1, ri, rgi
        Real(dp), Dimension(IP6) :: Y

        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0

        ! Core density:
        r1=r(1)
        Do ni=1,Nso
            qa=Qq(ni)
            call ReadF (12,ni+4,P,Q,2)
            If (ni.EQ.1) gm1=2*P(ii+4)
            gm=2*P(ii+4)
            igm=gm-gm1+0.1d0
            rg=r1**igm
            Do i=1,ii,ih
                C(i)=C(i)+qa*(P(i)**2+Q(i)**2)
            End Do
            Do m=0,MaxT-igm            !### expansion at the origin
                km=ii+5+m+igm            !#### of the electron density
                s=0.d0
                Do i=0,m
                    ki=ii+5+i
                    ki1=ii+5+m-i
                    s=s+(P(ki)*P(ki1)+Q(ki)*Q(ki1))
                End Do
                C(km)=C(km)+s*qa*rg
            End Do
        End Do
        C(ii+4)=gm1
        ! evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        Do i=1,i1,ih
             j=k1-i
             If (dabs(c(j)).GE.t) Exit
             imax=imax-ih
        End Do
        If (imax.LT.i0) imax=i0
        ! evaluation of z(0;r)=int_o^r(ro*dr) by direct integration
        r1=R(1)
        ri=R(i0)
        gm=C(ii+4)
        rg1=r1**(gm+1)
        rgi=ri**(gm+1)
        ri1=ri/r1
        t2=0.d0
        t1=0.d0
        fr=1.d0
        Do m=0,maxt              !### t1=int_0^r1 (ro*dr)
            im=ii+5+m              !### t2=int_0^r(i0) (ro*dr)
            C(im)=C(im)/(gm+m+1)
            t2=t2+C(im)*rg1
            t1=t1+C(im)*fr*rgi
            fr=fr*ri1
        End Do

        p1=c( 1)*v( 1)
        p2=c(i0)*v(i0)
        c(1)=t2
        c(i0)=t1
        i1=i0+ih
        Do i=i1,imax,ih
            q1=c(i)*v(i)
            t=hh*(q1+4.d0*p2+p1)
            t=t2+t
            t2=t1
            t1=t
            p1=p2
            p2=q1
            c(i)=t
        End Do
        t=c(imax)
        im1=imax+ih
        Do i=im1,ii,ih
            c(i)=t
        End Do
        C(ii+4)=gm1+1
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! evaluation of function y(0;r) by solving dif. equation
        ! for f=y(0;r)/r
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=imax+ih
        t1=c(i1)/r(i1)
        p1=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=i1+ih
        t2=c(i1)/r(i1)
        p2=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=imax+1
        Do j=1,imax,ih
            i=i1-j
            q1=-c(i)*v(i)/(r(i)*r(i))
            t=hh*(q1+4.d0*p1+p2)
            t=t2-t
            t2=t1
            t1=t
            p2=p1
            p1=q1
            y(i)=t
        End Do

        Do i=im1,ii,ih
          y(i)=c(i)/r(i)
        End Do

        y(ii+4)=0.d0                       !### Y0 -> const as r->0
        rg=r1**gm1
        s=0.d0
        Do m=0,MaxT
            s=s+rg*C(ii+5+m)/(gm1+m)
        End Do

        Y(ii+5)=Y(1)+s
        If (dabs(gm1-2.d0).LT.1.d-5) Then  !### the case of finite nucleus
            Y(ii+6)=0.d0
            Do m=2,MaxT
                Y(ii+5+m)=-rg*C(ii+3+m)/m
            End Do
        Else                               !### the case of point nucleus
            Do m=1,MaxT
                Y(ii+5+m)=0.d0
            End Do
        End If

       Return
    End Subroutine Y0

    Subroutine Yk(k)            
        !# This routine goes with integration routine Sint, which integrates over ro,
        !  while Sint1 integrates over r.
        Implicit None  

        Integer :: ih, i0, j, k, i, k1, imax, i1, im1
        Real(dp) :: hh, ww, rr, q1, t, p1, p2, gam, cr1, cr2, t1, t2, p
        Real(dp), Dimension(IP6) :: w

        ! evaluation of functions z,y with simpson formula
        ih=2-kt
        i0=ih+1
        hh=h0*ih/3.d0

        ! evaluation of w(i)=r(i)**k
        j=k/2
        k1=k-j*2
        Do i=1,ii,ih
            ww=1.d0
            rr=r(i)
            If (k1.EQ.1) ww=rr
            If (j.EQ.0) Then
                w(i)=ww
                Cycle
            End If
            q1=rr*rr
            ww=ww*q1
            If (j.EQ.1) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            If (j.EQ.2) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            If (j.EQ.3) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            w(i)=ww
        End Do

        ! evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        Do i=1,i1,ih
           j=k1-i
           If (dabs(c(j)).GE.t) Exit
           imax=imax-ih
        End Do

        If (imax.LT.i0) imax=i0

        ! evaluation of function z(k;r) by direct integration
        p1=c( 1)*v( 1)*w( 1)
        p2=c(i0)*v(i0)*w(i0)
        gam=c(ii+4)
        cr1=0.d0
        cr2=0.d0
        If (c(1).NE.0.d0) cr1=(c(i0)/c(1)*(r(1)/r(i0))**gam-1.d0)/((r(i0)-r(1))*(gam+k+2))
        If (c(i0).NE.0.d0) cr2=(c(1)/c(i0)*(r(i0)/r(1))**gam-1.d0)/((r(1)-r(i0))*(gam+k+2))
        c(1)=c(1)*r(1)/(gam+k+1)*(1.d0-cr1*r(1))
        c(i0)=c(i0)*r(i0)/(gam+k+1)*(1.d0-cr2*r(i0))
        t1=c(i0)*w(i0)
        t2=c(1)*w(1)
        i1=i0+ih
        Do i=i1,imax,ih
            q1=c(i)*v(i)
            If (k.NE.0) q1=q1*w(i)
            t=hh*(q1+4.d0*p2+p1)
            t=t2+t
            t2=t1
            t1=t
            p1=p2
            p2=q1
            c(i)=t/w(i)
        End Do
        t=c(imax)*w(imax)
        im1=imax+ih
        Do i=im1,ii,ih
            p=t
            If (k.NE.0) p=p/w(i)
            c(i)=p
        End Do

        ! evaluation of w(i)=r(i)**(k+1)
        Do i=1,ii,ih
            t=r(i)
            If (k.NE.0) t=t*w(i)
            w(i)=t
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! evaluation of function y(k;r) by solving dif. equation
        ! for f=y(k;r)/r**(k+1)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        k1=k+k+1
        i1=imax+ih
        t1=c(i1)/w(i1)
        p1=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=i1+ih
        t2=c(i1)/w(i1)
        p2=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=imax+1
        Do j=1,imax,ih
            i=i1-j
            q1=-c(i)*k1*v(i)/(r(i)*w(i))
            t=hh*(q1+4.d0*p1+p2)
            t=t2-t
            t2=t1
            t1=t
            p2=p1
            p1=q1
            c(i)=t*w(i)*v(i)*hh
        End Do

        Do i=im1,ii,ih
            c(i)=c(i)*v(i)*hh
        End Do

        If (c(ii+4).LE.k+0.d0) Then
            c(ii+4)=c(ii+4)+1
        Else
            c(ii+4)=k+1
        End If

        Return
    End Subroutine Yk

End Module sintg