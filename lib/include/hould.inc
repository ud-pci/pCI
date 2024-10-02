c     =================================================
      subroutine Hould(n,ndim,E,D,Z)
       implicit real*8 (a-h,o-z)
       common  /Ifail/Ifail
       dimension E(ndim),D(ndim),Z(ndim,ndim)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c        Householder's metod of diagonalisation
c        D-array of eigenvalues
c        Z-matrix of eigenvectors
c        n-dimension  of matrix z
c        eps-criteria of diagonalisation
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        Ifail=0
CCC        tol=2.0d0**(-103)
CCC        eps=2.0d0**(-24)
c         tol=2.0d0**(-153)
c         eps=2.0d0**(-34)
        tol=2.d0**(-1021)
        eps=2.d0**(-52)
        if (n.GT.1) then
          do ii=2,n
            i=n-ii+2
            l=i-2
            f=Z(i,i-1)
            g=0.0d0
            if (l.GT.0) then
              do k=1,l
                g=g+Z(i,k)**2
              end do
            end if
            h=f**2+g
            if (g.LT.tol) then
              E(i)=f
              h=0.0d0
              goto 10
            end if
            l=l+1
            g=-dsign(dsqrt(h),f)
            E(i)=g
            h=h-f*g
            Z(i,i-1)=f-g
            f=0.0d0
            do j=1,l
              Z(j,i)=Z(i,j)/h
              g=0.0d0
              do k=1,j
                g=g+Z(j,k)*Z(i,k)
              end do
              k1=j+1
              if (k1.LE.l) then
                do k=k1,l
                  g=g+Z(k,j)*Z(i,k)
                end do
              end if
              e(j)=g/h
              f=f+g*Z(j,i)
            end do
            hh=f/(h+h)
            do j=1,l
              f=Z(i,j)
              g=E(j)-hh*f
              E(j)=g
              do k=1,j
                Z(j,k)=Z(j,k)-f*E(k)-g*Z(i,k)
              end do
            end do
 10         D(i)=h
          end do
        end if
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        D(1)=0.0d0
        E(1)=0.0d0
        do i=1,n
           l=i-1
           if (D(i).NE.0.0d0.AND.l.NE.0) then
             do j=1,l
               g=0.0d0
               do k=1,l
                 g=g+Z(i,k)*Z(k,j)
               end do
               do k=1,l
                 Z(k,j)=Z(k,j)-g*Z(k,i)
               end do
             end do
           end if
           D(i)=Z(i,i)
           Z(i,i)=1.0d0
           if (l.NE.0) then
             do j=1,l
               Z(i,j)=0.0d0
               Z(j,i)=0.0d0
             end do
           end if
         end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     p matrix is formed
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (n.NE.1) then
          do i=2,n
            E(i-1)=E(i)
          end do
        end if
        E(n)=0.0d0
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     triadigonalisation is finished
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        b=0.0d0
        f=0.0d0
        do l=1,n
          jj=0
          h=(dabs(D(l))+dabs(E(l)))*eps
          if (b.LT.h) b=h
          do m1=l,n
            m=m1
            em=dabs(E(m))
            if (em.LE.b) goto 260
          end do
 260      if (m.NE.l) then
 270        if (jj.EQ.30) then
              ifail=1
              return
            end if
            jj=jj+1
            g=D(l)
            p=(D(l+1)-g)/E(l)*0.5d0
            r=dsqrt(p**2+1.0d0)
            em=p+dsign(r,p)
            D(l)=E(l)/em
            h=g-D(l)
            l1=l+1
            do i=l1,n
              D(i)=D(i)-h
            end do
            f=f+h
            p=D(m)
            c=1.0d0
            s=0.0d0
            im1=m-1
            if (im1.GE.l) then
              do im=l,im1
                i=im1+l-im
                ei=E(i)
                g=c*ei
                h=c*p
                if (dabs(p).GE.dabs(ei)) then
                  c=ei/p
                  r=dsqrt(c**2+1.0d0)
                  E(i+1)=s*p*r
                  s=c/r
                  c=1.0d0/r
                else
                  c=p/ei
                  r=dsqrt(c**2+1.0d0)
                  E(i+1)=s*ei*r
                  s=1.0d0/r
                  c=c/r
                end if
                di=D(i)
                p=c*di-s*g
                D(i+1)=h+s*(c*g+s*di)
                do k=1,n
                  h=Z(k,i+1)
                  ei=Z(k,i)
                  Z(k,i+1)=s*ei+c*h
                  Z(k,i)=c*ei-s*h
                end do
              end do
            end if
            E(l)=s*p
            D(l)=c*p
            if (dabs(E(l)).GT.b) goto 270
          end if
          D(l)=D(l)+f
        end do
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     eigenvalues in D; eigenvectors in Z
c     ordering of D and Z
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        i=0
 320    i=i+1
        if (i.GE.n) return
        h=D(i)
        f=D(i+1)
        if (h.LE.f) goto 320
        D(i)=f
        D(i+1)=h
        do k=1,n
           h=Z(k,i)
           Z(k,i)=Z(k,i+1)
           Z(k,i+1)=h
        end do
        if (i.NE.1) i=i-2
        goto 320
      end