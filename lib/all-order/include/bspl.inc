c     =================================================
      Subroutine Bspl(ks,              !rank of splines
     >                ns,              !number of splines
     >                m1,              !grid multiplier
     >                mx,              !No of nodes
     >                nout,            !spline to form
     >                C)               !formed spline
      implicit real*8 (a-h,o-z)
      include "bspl.par"
      include "hfd.par"
       common /Ii/Ii
     >        /Bsp/B(IPisp,IPnsp)  !array of splines
     >        /Rsp/Ro(IPisp)       !Radial grid
       dimension C(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (ns.GT.IPnsp.OR.mx.GT.IPisp) then
          write(*,*) ' Bspl error:'
          write(*,*) ' No more than',IPnsp,' B-splines on',
     >               IPisp,' grid points'
          write(*,*) ' Here Nsp=',ns,' mx=',mx
          stop
        end if

        do n=1,ns
          do m=1,mx
            B(m,n)=0.d0
          end do
        end do

        do n=1,ns
          mn=m1*(n-1)
          do m=1,m1
            B(mn+m,n)=1.d0
          end do
        end do

        do k=2,ks
          mk=m1*k
          do n=1,ns
            mn=m1*(n-1)+1
            r1=Ro(mn)
            r2=Ro(mn+m1)
            r4=Ro(mn+mk)
            r3=Ro(mn+mk-m1)
            do m=1,mk
              mm=mn+m-1
              rm=Ro(mm)
              b1=0.d0
              if (m.LE.(mk-m1))
     >          b1=b1+(rm-r1)/(r3-r1)*B(mm,n)
              if (m.GT.m1.AND.n.LT.ns)
     >          b1=b1+(r4-rm)/(r4-r2)*B(mm,n+1)
              B(mm,n)=b1
            end do
          end do
        end do

        im=m1*ks
        i1=Ii-1
        ns=ns-ks           ! initial number of splines is restored
        do n=1,ns
          do m=1,mx
            if (m.LE.i1) then
              B(m,n)=B(m+im,n+1)
            else
              B(m,n)=0.d0
            end if
          end do
        end do

 200    do i=1,Ii
          C(i)=B(i,nout)
        end do
        C(Ii+4)=9.d0
       return
      end
c     =================================================
      Subroutine Grid_Bspl(ksp,   !rank of splines
     >                     nsp,   !number of splines
     >                     m1,mx) !Forms grid for B-splines
      implicit real*8 (a-h,o-z)   !m1 - grid multiplier, mx - No of nodes
      include "bspl.par"
      include "hfd.par"
       common /Ii/Ii              !Original No of grid nodes
     >        /R/R(IP6)           !Original radial grid
     >        /Rsp/Ro(IPisp)      !Radial grid for B-splines
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        small=1.d-6               ! small factor is used on grid ends
        i1=Ii-1
        m1=i1/Nsp
        i1=m1*Nsp                 ! last (Ii-i1) nodes are skipped!

        Nsp=Nsp+Ksp               ! Ksp extra splines used in construction
        dl=(R(2)-R(1))*small      ! step for extra nodes on the left
        dr=(R(i1+1)-R(i1))*small  ! step for extra nodes on the right
        im=m1*Ksp
        mx=i1+2*im                ! Number of grid nodes
        rl=R(1)
        rr=R(i1+1)

        do i=1,im
          Ro(i)=rl-dl*(im-i-1)
          Ro(im+i1+i)=rr+dr*(i-1)
        end do

        do i=1,i1
          Ro(i+im)=R(i)
        end do

        write(11,5) 1,R(1),i1,R(i1),Ii-i1-1,2*im,mx,m1
 5      format(/4X,'B-sline grid: R(',I2,')=',F9.6,' to R(',I3,')=',
     >         F9.6,/4X,I3,' last nodes are skipped, ',I3,
     >         ' nodes added, giving ',I3,' nodes',
     >         /4X,'multiplier: ',I3)

       return
      end