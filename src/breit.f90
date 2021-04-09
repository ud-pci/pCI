Module breit

    Use basc_variables

    Implicit None

    Private

    Public :: Gaunt, Gaun, Br_core, breit_int, breit_pot_magn, breit_pot_ret, sbint, coefb, brint_magn, brint_ret

  Contains

    Subroutine Gaunt
        Use wigner
        Implicit None

        Integer :: l, lx, jx, ij1, ij2, jj2, kmin, kmax, kk, k, ind1, im1, im2, ind, mk, i
        Real(dp) :: J1,M1,J2,M2, x

        Open (unit=12,file='CONF.GNT',status='UNKNOWN',form='unformatted')

        l=0
        lx=IPlx+1            !### = max(l)
        jx=2*lx-1            !### = max(2j)
        Do ij1=1,lx
            j1=ij1-0.5
            Do ij2=ij1,lx
                j2=ij2-0.5
                jj2=2*ij2
                kmin=ij2-ij1+1
                kmax=ij1+ij2
                If (kmin.LT.2) kmin=2
                If (kmax.LT.kmin) Cycle
                Do kk=kmin,kmax
                    k=kk-1
                    ind1=jx*jx*k + 2*jx*j1 + 2*j2 + 0.1d0
                    Do im1=1,ij1
                        m1=im1-0.5
                        Do im2=1,jj2
                            m2=im2-j2-1
                            ind=jx*jx*ind1 + jx*(j1+m1) + (j2+m2) + 0.1d0
                            mk=(m2-m1)
                            If(iabs(mk).gt.k) Cycle
                            l=l+1
                            in(l)=ind
                            x=Gaun(K,J1,M1,J2,M2)
                            gnt(l)=x
                        End Do
                    End Do
                End Do
            End Do
        End Do

        If (l.NE.IPgnt) Then
           Write (*,*) ' number of gaunts ',l,' NE IPgnt'
           Stop
        End If

        Write(12) (in(i),i=1,l)
        Write(12) (gnt(i),i=1,l)

        Close(unit=12)

        Return
    End Subroutine Gaunt

    Real(dp) Function Gaun(K,J1,M1,J2,M2)
        Use wigner
        Implicit None

        Integer :: i, is,K
        Real(dp) :: J1, M1, J2, M2, xk, q

        is=dabs(m1+0.5d0)+0.1d0
        If (is.NE.(is/2)*2) Then
           is=-1
        Else
           is=1
        End If

        xk=k
        q=m1-m2
        Gaun=is*dsqrt((2*j1+1)*(2*j2+1)) &
             * Fj3(j1,j2,xk,-m1,m2,q) &
             * Fj3(j1,j2,xk,0.5d0,-0.5d0,0.d0)

        Return
    End Function Gaun

    Real(dp) Function Br_core(na,nb,kfl) 
        !### exchange breit core potential
        !### kfl - file with derivatives
        !     Exchange core potential:
        !     <b|V1_B+V2_B|a> = sum_c sum_k Gaunt^k_{a,c,c,b} R^k_{a,c,c,b}
        Use readfff
        Use wigner
        Implicit None

        Integer :: la, lb, lc, kmin, kmax, k, na, nb, nc, ja, jb, jc, kfl
        Real(dp) :: z12, z00,  xja, a_ab, xjc, s, xj, gac, ds
        Real(dp), Dimension(10) :: Rint1, Rint2 !### Radial integrals after regrouping
        Character(Len=512) :: strfmt

        Qb = Qd
        Pb = Pd

        z12=0.5d0
        z00=0.0d0
        Br_core=z00
        If (Nso.EQ.0) Return
 
        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)
        If (ja.NE.jb.OR.la.NE.lb) Return
       
        xja=0.5d0*ja
 
        a_ab=0.d0
        Do nc=1,Nso
            Call ReadF(12,nc+4,Pc,Qc,2)
            jc=Jj(nc)
            lc=Ll(nc)
            xjc=0.5d0*jc

            ! Calculation of both Breit contributions using IIT Function
            kmin=abs(xja-xjc)
            kmax=xja+xjc
            s = 0.d0
            Do k=kmin,kmax
                xj=k
                gac=-(2*xjc+1)*FJ3(xja,xjc,xj,z12,-z12,z00)**2
                ds = 0.d0
                If (dabs(gac).gt.1.d-7) &
                ds = breit_int(k,na,Pa,Qa,nc,Pc,Qc,nc,Pc,Qc,nb,Pb,Qb)
                s = s+ds*gac
            End Do
            If (kout.GT.2) Then
                strfmt='(2X,"contribution of ",I2,A1,I1,"/2 :",E12.5)'
                Write( *,strfmt) Nn(nc),let(lc+1),jc,s
                Write(11,strfmt) Nn(nc),let(lc+1),jc,s
            End If
            a_ab=a_ab+s
        End Do

        If (kout.GT.2) Then
            strfmt='(/2X,"<",I2,A1,I1,"/2| V_B^core |",I2,A1,I1,"/2> = ",E12.5)'
            Write( *,strfmt) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
            Write(11,strfmt) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
        End If
        Br_core=a_ab

        Return
    End Function Br_core

    Real(dp) Function breit_int(l,na,pa,qa,nb,pb,qb,nc,pc,qc,nd,pd,qd)
        Implicit None

        Integer :: la, lb, lc, ld, i, l, na, nb, nc, nd
        Real(dp) :: ds1, ds2
        Real(dp), Dimension(IP6) :: pa,qa,pb,qb,pc,qc,pd,qd

        breit_int=0.d0
        la=Ll(na)
        lb=Ll(nb)
        lc=Ll(nc)
        ld=Ll(nd)
        i=la+lc+lb+ld

        If (2*(i/2).ne.i) Return

        Call breit_pot_magn(l,na,pa,qa,nc,pc,qc)

        ds1=brint_magn(l,nb,pb,qb,nd,pd,qd)
        ds2=0.d0

        If (kbrt.eq.2) Then
            Call breit_pot_ret(l,na,pa,qa,nc,pc,qc)
            ds2=brint_ret(l,nb,pb,qb,nd,pd,qd)
        End If

        breit_int=ds1+ds2

        Return
    End Function breit_int

    Real(dp) Function brint_magn(l,nb,pb,qb,nd,pd,qd)
        Implicit None

        Integer :: lb, ld, k, l, i, nb, nd
        Real(dp) :: s1, s2, ds1, ds2
        Real(dp), Dimension(IP6) :: pb,qb,pd,qd
        Real(dp), Dimension(IP6) :: c, ro

        lb=Ll(nb)
        ld=Ll(nd)

        s1=0.d0
        s2=0.d0

        Do k=l-1,l+1
            If (k.lt.0) Cycle
            If (l.eq.0.and.k.eq.0) Cycle
            i=lb+ld+k
            If (2*(i/2).eq.i) Cycle
    
            Do i=1,IP6
              c(i)=yy(i,k+1)
            End Do
    
            Call sbint(pb,qd,c,ds1)
            s1=s1+coefb(-1,k,l,nd,nb)*ds1
    
            Call sbint(qb,pd,c,ds2)
            s2=s2+coefb( 1,k,l,nd,nb)*ds2
        End Do

        brint_magn=-2*(s1+s2)

        Return
    End Function brint_magn

    Subroutine sbint(a,b,c,s)
        Implicit None
        
        Integer :: i0, imax, i, j, j1, n, k
        Real(dp) :: r0, v0, g, t0, p0, dw, f0, dt, f, s
        Real(dp), Dimension(IP6) :: a,b,c

        i0=1
        r0=r(i0)
        v0=v(i0)

        s=0.d0
        imax=ii   ! mgk 
        g=a(ii+4)+b(ii+4)+c(ii+4)-1.d0
        t0=0.d0
        p0=0.d0
        Do m=0,MaxT
            i=ii+5+m
            dw=0.d0
                Do n=0,m
                j=ii+5+n
                    Do k=0,m-n
                    j1=ii+5+k
                      dw=dw+a(j)*b(j1)*c(i-n-k)
                    End Do
                End Do
            t0=t0+dw/(g+m+1)
            p0=p0+dw*(g+m)
        End Do

        t0=t0*r0**(g+1)
        p0=p0*r0**(g-1)
        f0=c(i0)*(a(i0)*b(i0))*v0/r(i0)
        p0=p0*v0*v0+f0*bt/(al*r0+bt)**2

200     dt=0.d0
        Do i=i0,imax
          f=c(i)*(a(i)*b(i))
          If (r2.gt.0.d0) f=f*v(i)/r(i)
          dt=dt+f
        End Do
        s=t0+h*dt-h*(0.5d0*f0-h/12.d0*p0)-0.5d0*h*f

        Return
    End Subroutine sbint

    Real(dp) Function coefb(ibet,k,l,na,nb)
        Implicit None

        Integer :: ka, kb, ibet, k, l, na, nb, lk
        Real(dp) :: u, alk, d, blk

        u=0.d0

        ka=Kk(na)
        kb=Kk(nb)

        alk=0.5d0*(k*(k+1)-l*(l+1))
        lk=k+l
        blk=1.d0
        If (k.eq.l) blk=1.d0/dsqrt(2.d0)
        u=dsqrt((2*l+1.d0)/(lk*(lk+1)*(lk+2)))/blk
        d=kb
        If (2*(lk/2).ne.lk) d=-d
        u=u*(ka+d+ibet*alk)

        coefb=u

        Return
    End Function coefb

    Real(dp) Function brint_ret(l,nb,pb,qb,nd,pd,qd)
        Implicit None

        Integer :: lb, ld, l, nb, nd, k1, k2, i, k
        Real(dp) :: s1, s2, dk, ck1, ck2, c12, ds1, ds2
        Real(dp), Dimension(IP6) :: pb,qb,pd,qd
        Real(dp), Dimension(IP6) :: c

        lb=Ll(nb)
        ld=Ll(nd)

        s1=0.d0
        s2=0.d0

        Do k1=l-1,l+1,2
            If (k1.lt.0) Cycle
            i=lb+ld+k1
            If (2*(i/2).eq.i) Cycle
    
            dk=k1
            If (k1.eq.l+1) Then
                ck1=dsqrt(dk/(2*l+1))
            Else
                ck1=-dsqrt((dk+1.d0)/(2*l+1))
            End If

            If (dabs(ck1).lt.1.d-10) Cycle
    
            Do k2=l-1,l+1,2
                If (k2.lt.0) Cycle
                i=lb+ld+k2
                If (2*(i/2).eq.i) Cycle
        
                dk=k2
                If (k2.eq.l+1) Then
                    ck2=dsqrt(dk/(2*l+1))
                Else
                    ck2=-dsqrt((dk+1.d0)/(2*l+1))
                End If
                If (ck2.eq.0.d0) Cycle
        
                c12=ck1*ck2*dsqrt((2*k1+1.d0)*(2*k2+1.d0))
                If (k1.eq.k2) c12=2.d0*c12/(2*k1+1)
                If (dabs(c12).lt.1.d-10) Cycle
        
                If (k1.eq.k2) Then
                    Do i=1,IP6
                        c(i)=yy(i,k1+1)
                    End Do
                Else
                    If (k2.eq.k1+2) k=k1
                    If (k1.eq.k2+2) k=k2
        
                    If (k2.eq.k1+2) Then
                        Do i=1,IP6
                          c(i)=uu(i,k+1)
                        End Do
                    End If
            
                    If (k1.eq.k2+2) Then
                        Do i=1,IP6
                            c(i)=vv(i,k+1)
                        End Do
                     End If
                End If
        
                Call sbint(pb,qd,c,ds1)
                s1=s1+c12*coefb( 1,k2,l,nb,nd)*ds1
        
                Call sbint(qb,pd,c,ds2)
                s2=s2+c12*coefb(-1,k2,l,nb,nd)*ds2
            End Do
        End Do

        brint_ret=-(s1+s2)

        Return
    End Function brint_ret

    Subroutine breit_pot_magn(l,na,pa,qa,nc,pc,qc)
        Implicit None

        Integer :: la, lc, l, na, nc, i, k, k7
        Real(dp) :: g1, g2

        Real(dp), Dimension(IP6) :: pa,qa,pc,qc
        Real(dp), Dimension(IP6) :: ro,c 

        la=Ll(na)
        lc=Ll(nc)

        Do k=l-1,l+1
            If (k.lt.0) Cycle
            If (l.eq.0.and.k.eq.0) Cycle
            i=la+lc+k
            If (2*(i/2).eq.i) Cycle
    
            If (k.gt.maxk) Then
                Write( k7,'(/2x,a/2x,a)')'*** k.gt.maxk in brin_m_fast ***', 'Increase maxk'
                Call exit(1)
            End If
    
            g1=coefb( 1,k,l,na,nc)
            g2=coefb(-1,k,l,na,nc)
            Call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
            Call ykt(k,ro,c)

            Do i=1,IP6
                yy(i,k+1)=c(i)
            End Do
        End Do

        Return
    End Subroutine breit_pot_magn


    Subroutine breit_pot_ret(l,na,pa,qa,nc,pc,qc)
        Implicit None

        Integer :: la, lc, l, na, nc, i, k, k1, k2, k7
        Real(dp) :: g1, g2, dk, ck2
        Real(dp), Dimension(IP6) :: pa,qa,pc,qc
        Real(dp), Dimension(IP6) :: ro,c

        la=Ll(na)
        lc=Ll(nc)

        Do k1=l-1,l+1,2
            If (k1.lt.0) Cycle
            i=la+lc+k1
            If (2*(i/2).eq.i) Cycle
    
            g1=coefb( 1,k1,l,na,nc)
            g2=coefb(-1,k1,l,na,nc)
            Call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
    
            Do k2=l-1,l+1,2
                If (k2.lt.0) Cycle
                i=la+lc+k2
                If (2*(i/2).eq.i) Cycle
        
                dk=k2
                If (k2.eq.l+1) Then
                    ck2=dsqrt(dk/(2*l+1))
                Else
                    ck2=-dsqrt((dk+1.d0)/(2*l+1))
                End If
                If (ck2.eq.0.d0) Cycle
        
                g1=coefb( 1,k1,l,na,nc)
                g2=coefb(-1,k1,l,na,nc)
                Call rho_pq_tot(g1,g2,pa,qa,pc,qc,ro)
        
                If (k1.eq.k2) Then
                    If (k1.gt.maxk) Then
                        Write( k7,'(/2x,a/2x,a)')'*** k1.gt.maxk in brin_r_fast ***', 'Increase maxk'
                        Call exit(1)
                    End If
            
                    Call ykt(k1,ro,c)
            
                    Do i=1,IP6
                        yy(i,k1+1)=c(i)
                    End Do
                Else
                    If (k2.eq.k1+2) k=k1
                    If (k1.eq.k2+2) k=k2
                    If (k.gt.maxk) Then
                        Write( k7,'(/2x,a/2x,a)')'*** k.gt.maxk in brin_r_fast ***', 'Increase maxk'
                        Call exit(1)
                    End If
            
                    If (k2.eq.k1+2) Then
                        Call ukt(k,ro,c)
                        Do i=1,IP6
                          uu(i,k+1)=c(i)
                        End Do
                    End If
            
                    If (k1.eq.k2+2) Then
                        Call vkt(k,ro,c)
                        Do i=1,IP6
                           vv(i,k+1)=c(i)
                        End Do
                    End If
                End If
            End Do
        End Do

        Return
    End Subroutine breit_pot_ret

    Subroutine rho_pq_tot(coef1,coef2,p,q,a,b,ro)
        Implicit None

        Integer :: i, imax, j, n, ih
        Real(dp) :: coef1, coef2, d1, d2, ulam, c1, dh1, dh2, dr
        Real(dp), Dimension(IP6) :: p,q,a,b,ro

        c1=0.01

        IH=2-KT
        H=H0*IH
        ulam = 1

        imax=ii          ! mgk
        dh1=ulam*h*coef1
        dh2=ulam*h*coef2
        Do i=1,imax
            ro(i)=(dh1*p(i)*b(i)+dh2*q(i)*a(i))*v(i)
        End Do

        ro(imax+1:ii)=0.d0
        ro(ii+3)=imax
        ro(ii+4)=p(ii+4)+a(ii+4)

        d1=ulam*coef1
        d2=ulam*coef2
        Do m=0,MaxT
            i=ii+5+m
            dr=0.d0
            Do n=0,m
                j=ii+5+n
                dr=dr+(d1*p(j)*b(i-n)+d2*q(j)*a(i-n))
            End Do
            ro(i)=ulam*dr
        End Do

        Return
    End Subroutine rho_pq_tot

    Subroutine ykt(k,ro,c)
        Implicit None

        Integer :: i0, imax, k, ih, i1, i, im, j, ip, id
        Real(dp) :: r0, v0, dk1, dk2, g, t0, p0, f0, fm, t, d, fi, dh, s, g0, r0d
        Real(dp), Dimension(IP6):: ro,c,w

        IH=2-KT
        H=H0*IH

        ! Calculation functions z(k,r) and y(k,r)
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk1=1.d0-dexp(-k*h)
        dk2=1.d0-dexp(-(k+1)*h)

        If (r2.ge.0.d0) Then
            i1=i0+1
            Do i=i1,ii
                w(i)=0.d0
                If (k.ne.0) w(i)=1.d0-(r(i-1)/r(i))**k
            End Do
        End If

        g=ro(ii+4)+k
        t0=0.d0
        p0=0.d0
        Do m=0,MaxT
            i=ii+5+m
            t0=t0+ro(i)/(g+m+1)
            p0=p0+ro(i)*(g+m)
        End Do
        t0=t0*r0**(g+1-k)
        p0=p0*r0**(g-1-k)
        f0=ro(i0)
        p0=h*p0*v0*v0+f0*bt/(al*r0+bt)**2 
        fm=ro(imax)

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation of the Function z(k;r)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        t=t0+0.5d0*f0+h/12.d0*p0
        c(i0)=t
        d=dk1
        i1=i0+1
        Do i=i1,imax
            fi=ro(i)
            If (r2.gt.0.d0) d=w(i)
            If (k.ne.0) fi=fi-d*t
            t=t+fi
            c(i)=t
        End Do

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation of the Function y(k;r)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        t=-0.5d0*fm
        dh=(2*k+1)*h/12.d0
        d=dk2
        im=imax-1
        i1=im+i0

        Do j=i0,im
            i=i1-j
            ip=i+1
            s=dh
            If (r2.ge.0.d0) Then
                d=r(i)/r(ip)
                If (k.ne.0) d=d*(1.d0-w(ip))
                d=1.d0-d
                s=dh*v(i)/r(i)
            End If
            fi=ro(ip)
            fi=fi-d*(t+fi)
            t=t+fi
            c(i)=c(i)+t-s*ro(i)
        End Do

        d=dk1
        t=c(imax)-0.5d0*fm

        Do i=imax,ii
            c(i)=t
            If (k.eq.0) Cycle
            If (r2.gt.0.d0) d=w(i)
            t=t-d*t
        End Do

        g0=ro(ii+4)
        Do m=0,MaxT
            i=ii+5+m
            c(i)=0.d0
        End Do
        c(ii+3)=ii
        c(ii+4)=1.d0
        id=g0+0.5d0
        r0d=r0**id

        If (g0-k.gt.1.d-4) Then
            t0=0.d0
            Do m=0,MaxT
                i=ii+5+m
                t0=t0+(ro(i)/(g0+m-k)-ro(i)/(g0+m+k+1))
            End Do
            t0=t0*r0**g0
            If (k.lt.MaxT) c(ii+5+k)=c(i0)/r0+t0
    
            Do m=0,MaxT
                i=ii+5+m
                If ((m+id).le.MaxT) c(i+id)= &
                ro(i)*(1.d0/(g0+k+m+1)-1.d0/(g0+m-k))*r0d
            End Do
        Else
            c(ii+5+id)=c(1)/r(1)**(id+1)*r0d
        End If

        Return
    End Subroutine ykt

    Subroutine ukt(k,ro,c)
        Implicit None

        Integer :: ih, k, i0, imax, i1, i, id
        Real(dp) :: r0, v0, dk1, g, t0, p0, t2, p2, f0, zi, d, d12, d2
        Real(dp) :: fi, fm, t, dh1, g0, z1, r0d
        Real(dp), Dimension(IP6) :: ro,c,w

        IH=2-KT
        H=H0*IH
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation the Function uk(k,r)=z(k+2,r)-z(k,r).
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk1=1.d0-dexp(-k*h)

        If (r2.ge.0.d0) Then
            i1=i0+1
            w(i1:ii)=0.d0
            Do i=i1,ii
                If (k.ne.0) w(i)=1.d0-(r(i-1)/r(i))**k
            End Do
        End If

        g=ro(ii+4)+k
        t0=0.d0
        p0=0.d0
        t2=0.d0
        p2=0.d0
        Do m=0,MaxT
            i=ii+5+m
            t0=t0+ro(i)/(g+m+1)
            t2=t2+ro(i)/(g+m+3)
            p0=p0+ro(i)*(g+m)
            p2=p2+ro(i)*(g+m+2)
        End Do

        t0=t0*r0**(g+1-k)
        t2=t2*r0**(g+1-k)
        p0=p0*r0**(g-1-k)
        p2=p2*r0**(g-1-k)
        f0=ro(i0)
        p0=h*p0*v0*v0+f0*bt/(al*r0+bt)**2
        p2=h*p2*v0*v0+f0*bt/(al*r0+bt)**2

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation of the Function uk(k;r)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        zi=t0+0.5d0*f0+h/12.d0*p0
        z1=(t2+0.5d0*f0+h/12.d0*p2)-zi
        c(i0)=z1
        d=dk1
        i1=i0+1
        Do i=i1,imax
            If (r2.gt.0.d0) d=w(i)
    
            d12=r(i-1)/r(i)
            d12=d12*d12
            d2=d*d12+1.d0-d12
            z1=z1-d2*z1+(d-d2)*zi
    
            fi=ro(i)
            If (k.ne.0) fi=fi-d*zi
            zi=zi+fi
    
            c(i)=z1
        End Do

        d=dk1
        fm=ro(imax)
        t=c(imax)-0.5d0*fm
        Do i=imax,ii
            c(i)=t
            If (k.eq.0) Cycle
            If (r2.gt.0.d0) d=w(i)
            t=t-d*t
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Correction
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! zi=h/12.d0*p0
        ! c(i0)=c(i0)-h/12.d0*p2+zi
        dh1=2*h/12.d0
        Do i=i0,imax
        c(i)=c(i)-dh1*ro(i)*v(i)/r(i)
        End Do

        g0=ro(ii+4)
        id=g0+0.5d0
        r0d=r0**id

        Do m=0,MaxT
            i=ii+5+m
            c(i)=0.d0
        End Do

        c(ii+3)=ii
        c(ii+4)=1.d0+g0-id
        Do m=0,MaxT
            i=ii+5+m
            If (m+id.gt.MaxT) Cycle
            c(i+id)=ro(i)*(1.d0/(g0+k+2+m+1)-1.d0/(g0+k+m+1))*r0d
        End Do

        Return
    End Subroutine ukt

    Subroutine vkt(k,ro,c)
        Implicit None

        Integer :: ih, i0, imax, i1, k, i, im, j, ip, id
        Real(dp) :: r0, v0, dk2, d, g, fm, xi, d12, d2, fi, dh1, g0, r0d, t0
        Real(dp) :: x1, c1, c2, c3, x2, x3, a2, a1, a3, t2
        Real(dp), Dimension(IP6):: ro,c,w
        Real(dp), Dimension(10) :: v1,v2

        IH=2-KT
        H=H0*IH
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! calculation functions v(k,r)=x(k+2,r)-x(k,r)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=ro(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
        dk2=1.d0-dexp(-(k+1)*h)

        If (r2.ge.0.d0) Then
            i1=i0+1
            Do i=i1,ii
                d=r(i-1)/r(i)
                If (k.eq.0) w(i)=1.d0-d
                If (k.ne.0) w(i)=1.d0-d**(k+1)
            End Do
        End If

        g=ro(ii+4)+k
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation of the Function v(k;r)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        fm=ro(imax)
        xi=-0.5d0*fm
        x1=0.d0
        d=dk2
        im=imax-1
        i1=im+i0
        Do j=i0,im
            i=i1-j
            ip=i+1
            If (r2.gt.0.d0) d=w(ip)
    
            d12=r(ip-1)/r(ip)
            d12=d12*d12
            d2=d*d12+1.d0-d12
            x1=x1-d2*x1+(d-d2)*(xi+ro(ip))
    
            fi=ro(ip)
            fi=fi-d*(xi+fi)
            xi=xi+fi
    
            c(i)=x1
    
            If (i-i0+1.le.10) Then
                v1(i)=xi
                v2(i)=c(i)+xi
            End If
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Correction
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        dh1=2*h/12.d0
        Do i=i0,imax
            c(i)=c(i)-dh1*ro(i)*v(i)/r(i)
        End Do

        c(imax+1:ii)=0.d0

        g0=ro(ii+4)
        Do m=0,MaxT
            i=ii+5+m
            c(i)=0.d0
        End Do
        c(ii+3)=ii
        id=g0+0.5d0
        r0d=r0**id
        c(ii+4)=1.d0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Near the origin
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (g0-k.gt.1.d-4) Then
            t0=0.d0
            Do m=0,MaxT
                i=ii+5+m
                t0=t0+ro(i)/(g0+m-k)
            End Do
            t0=t0*r0**g0
            If (k.le.MaxT) c(ii+5+k)=c(ii+5+k)-v1(i0)/r0-t0
            Do m=0,MaxT
                i=ii+5+m
                If ((m+id).le.MaxT)c(i+id)=c(i+id)+ro(i)/(g0+m-k)*r0d
            End Do
        Else
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
            If (id.le.MaxT)   c(ii+5+id)=a1*r0d
            If (id+1.le.MaxT) c(ii+6+id)=a2*r0d*r0
            a3=(x3-x2)*(c2-c1)-(x2-x1)*(c3-c2)
            a3=a3/((x2*x2-x1*x1)*(x3-x2)-(x3*x3-x2*x2)*(x2-x1))
            a2=(c2-c1)-a3*(x2*x2-x1*x1)
            a2=a2/(x2-x1)
            a1=c1-a2*x1-a3*x1*x1
            If (id.le.MaxT)   c(ii+5+id)=c(ii+5+id)-a1*r0d
            If (id+1.le.MaxT) c(ii+6+id)=c(ii+6+id)-a2*r0d*r0
            If (id+2.le.MaxT) c(ii+7+id)=c(ii+7+id)-a3*r0d*r0**2
        End If

        If (g0-k-2.gt.1.d-4) Then
            t2=0.d0
            Do m=0,MaxT
                i=ii+5+m
                t2=t2+ro(i)/(g0+m-k-2)
            End Do
            t2=t2*r0**g0
            If (k+2.le.MaxT)c(ii+5+k+2)=c(ii+5+k+2)+v2(i0)/r0+t2
            Do m=0,MaxT
                i=ii+5+m
                If ((m+id).le.MaxT)c(i+id)=-ro(i)/(g0+m-k-2)*r0d+ro(i)/(g0+m-k)*r0d
            End Do
        Else
            id=g0+0.5d0-1
            r0d=r0**id
    
            c1=v2(1)/r(1)**(id+1)
            c2=v2(3)/r(3)**(id+1)
            c3=v2(5)/r(5)**(id+1)
            x1=r(1)
            x2=r(3)
            x3=r(5)
            a2=(c2-c1)/(x2-x1)
            a1=c1-a2*x1
            If (id.le.MaxT)   c(ii+5+id)=a1*r0d
            If (id+1.le.MaxT) c(ii+6+id)=a2*r0d*r0
    
            a3=(x3-x2)*(c2-c1)-(x2-x1)*(c3-c2)
            a3=a3/((x2*x2-x1*x1)*(x3-x2)-(x3*x3-x2*x2)*(x2-x1))
            a2=(c2-c1)-a3*(x2*x2-x1*x1)
            a2=a2/(x2-x1)
            a1=c1-a2*x1-a3*x1*x1
            If (id.le.MaxT)   c(ii+5+id)=c(ii+5+id)+a1*r0d
            If (id+1.le.MaxT) c(ii+6+id)=c(ii+6+id)+a2*r0d*r0
            If (id+2.le.MaxT) c(ii+7+id)=c(ii+7+id)+a3*r0d*r0**2
        End If
        
        Return
    End Subroutine vkt

End Module breit