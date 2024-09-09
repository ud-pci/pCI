C     =========== 03/06/00 ===== last changes: 27/10/06
c     ||         numerical differentiation           ||
C     =================================================
      subroutine Dif(P,CP,gm)                       ! This Dif should be
      implicit real*8 (a-h,o-z)                     ! used for functions,
      include "hfd.par"                             ! which are not very
       common /Ii/Ii/kt/kt/H/H/MaxT/MaxT/Ierr/Ierr  ! smooth at R=R(1).
     >        /V/V(IP6)/R/R(IP6)                    !
       dimension P(IP6),CP(IP6)                     !
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (MaxT.EQ.0) MaxT=9      !### hfd uses Nmax instead of MaxT
        ih=2-kt
        h60=60*ih*h
        CP(ii+4)=gm-1.d0

c>>>>>> Expansion at the origin:
        r1=R(1)
        cp1=0.d0                   !### = left derivative (dP/dr)_{-}
        p1=0.d0                    !### = P(r1)
        do m=0,MaxT
          im=ii+5+m
          pim=P(im)
          CP(im)=(gm+m)*pim
          cp1=cp1+CP(im)
          p1=p1+pim
        end do

        p1=p1*r1**gm
        if (dabs(p1).GT.1.d-19) then
          scale=P(1)/p1
        else
          scale=1.d0
        end if

        if(dabs(scale-1.d0).GT.1.d-3) then
          write(*,*) ' Dif: expansion does not match the function'
          write(*,*) ' P(1)=',P(1),' Taylor=',p1
          read(*,*)
          stop
        end if

        cp1=scale*cp1*r1**(gm-1.d0)      !### = CP(r1)
        p1=scale*p1                      !### = P(r1)
c<<<<<< Expansion at the origin

        ix=ii                      !### P(ix) is last nonzero value
        do j=1,ii,ih
          i1=ii-j+1
          if (P(i1).NE.0.d0) goto 200
          ix=i1-ih
          CP(i1)=0.d0
        end do

 200    continue


c>>>>> Symmetric 7-node differentiation for internal nodes:

        ih2=ih +ih
        ih3=ih2+ih
        ih4=ih3+ih

        do i=ih3+1,ix-ih3,ih
          CP(i)=(P(i+ih3)-9*(P(i+ih2)-P(i-ih2))
     >         +45*(P(i+ih)-P(i-ih))-P(i-ih3))
     >          /(h60*V(i))
        end do

c>>>>> Asymmetric 7-node differentiation for first three nodes:

        ih5=ih4+ih
        ih6=ih5+ih

        CP(1)=(-147*P(1)+360*P(1+ih)
     >       -450*P(1+ih2)+400*P(1+ih3)-225*P(1+ih4)
     >       +72*P(1+ih5)-10*P(1+ih6))/(h60*V(1))

        CP(1+ih)=(-10*P(1)-77*P(1+ih)
     >       +150*P(1+ih2)-100*P(1+ih3)+50*P(1+ih4)
     >       -15*P(1+ih5)+2*P(1+ih6))/(h60*V(1+ih))

        CP(1+ih2)=(2*P(1)-24*P(1+ih)-35*P(1+ih2)
     >       +80*P(1+ih3)-30*P(1+ih4)+8*P(1+ih5)
     >       -P(1+ih6))/(h60*V(1+ih2))


c>>>>> Asymmetric 7-node differentiation for last three nodes:

        CP(ix)=-(-147*P(ix)+360*P(ix-ih)
     >       -450*P(ix-ih2)+400*P(ix-ih3)-225*P(ix-ih4)
     >       +72*P(ix-ih5)-10*P(ix-ih6))/(h60*V(ix))

        CP(ix-ih)=-(-10*P(ix)-77*P(ix-ih)
     >       +150*P(ix-ih2)-100*P(ix-ih3)+50*P(ix-ih4)
     >       -15*P(ix-ih5)+2*P(ix-ih6))/(h60*V(ix-ih))

        CP(ix-ih2)=-(2*P(ix)-24*P(ix-ih)-35*P(ix-ih2)
     >       +80*P(ix-ih3)-30*P(ix-ih4)+8*P(ix-ih5)
     >       -P(ix-ih6))/(h60*V(ix-ih2))

c                                                     ! We assume that this
c        if (cp1.EQ.0.d0.AND.P(1).EQ.0.d0) CP(1)=0.d0 ! may occur only for
c                                                     ! B-spline.
        err=cp1-CP(1)
        if(dabs(cp1).GT.1.d0) err=err/cp1
        if(dabs(err).GT.1.d-1) then
          write(*,15) cp1,CP(1)
  15      format(4X,'Dif: matching error at R1. Taylor:',E12.5,
     >           ', Num:',E12.5,/4X,'Using Origin for new expansion.')
          Ierr=Ierr+10*dabs(err)
          if (dabs(CP(ii+5)).GT.dabs(CP(ii+6))) then
            kap=1                              ! For kap>0 expansion
          else                                 ! goes over even powers
            kap=-1
          end if
          call Origin(CP,CP(ii+4),kap)
          read(*,*)
        end if
       return
      end
C    =================================================
      subroutine TBRL(L)             ! L additional grid points
       implicit real*8 (a-h,o-z)     ! inside the nucleus
       include "hfd.par"             ! for numerical differentiation
       common /R/R(IP6)/V/V(IP6)
     >        /Rdop/Rdop(IP6)
        call ReadF(12,1,Rdop,Rdop,1) ! grid: ro=al*r+bt*ln(r)
        h =Rdop(6)                   ! grid for ro is equidistant d_ro=h
        bt=Rdop(7)                   ! corresponding d_r is found
        al=Rdop(8)                   ! iteratively.
        L1=L+1

        do i=1,IP6
          Rdop(i)=0.d0
        end do

        rr=R(1)
        Rdop(1)=rr
        Rdop(L1+1)=V(1)
        ro=al*rr+bt*dlog(rr)
        do i=1,L
          i1=i+1
          roi1=ro-h              ! - new value for ro
          ri1=Rdop(i)            ! - first approximation for new r
          ro=al*ri1+bt*dlog(ri1) ! - approximate new ro
 200      ri1 = ri1 + (roi1-ro)/(al+bt/ri1)
          ro=al*ri1+bt*dlog(ri1) ! - new approximation
          if(dabs(ro-roi1).GT.5.d-9) goto 200
          Rdop(i1)=ri1
          Rdop(L1+i1)=ri1/(al*ri1+bt)
        end do
        write( *,5) L,(2-i,Rdop(i),Rdop(L1+i),i=1,L1)
        write(11,5) L,(2-i,Rdop(i),Rdop(L1+i),i=1,L1)
 5      format(4X,'TBRL:',I2,' grid points inside the nucleus added:'
     >         /(6X,'R(',I2,')=',E12.5,', V =',E12.5))
       return
      end
C     =================================================
      subroutine Origin(P,gam,kap)           ! 19/05/08
C     >> variant with four terms of Taylor expansion <<
C     =================================================
      implicit real*8 (a-h,o-z)
      include "hfd.par"                       !# we now match P(1),P(2),
       common /Ii/Ii/MaxT/MaxT/Ierr/Ierr/H/H  !# dP/dr(1), and dP/dr(2)
       common /R/R(IP6)/V/V(IP6)              !#
       dimension P(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (MaxT.EQ.0) MaxT=9        !### hfd uses Nmax instead of MaxT
        do i=ii+5,ii+5+MaxT
          P(i)=0.d0
        end do

        if (P(1).EQ.0.d0) return

        if (kap.LT.0) then           !### for positive kappa expansion
          k=0                        !#### goes over odd powers
          rn=1.d0
        else
          k=1
          rn=R(1)
        end if
        g=gam+k

        dp1 =(-147*P(1)+360*P(2)           !# 7-node right derivative
     >       -450*P(3)+400*P(4)-225*P(5)   !# to be matched with
     >       +72*P(6)-10*P(7))/(60*h*V(1)) !# Taylor expansion
        dp2 =(-10*P(1)-77*P(2)
     >       +150*P(3)-100*P(4)+50*P(5)
     >       -15*P(6)+2*P(7))/(60*h*V(2))


 1      p1=P(1)*rn/R(1)**g     !# (see 26/10/06)
        p2=P(2)*rn/R(2)**g     !#  with addition
        d1=dp1 *rn/R(1)**(g-1) !#   (19/05/08)
        d2=dp2 *rn/R(2)**(g-1) !#
        r1=R(1)*R(1)
        y=R(2)*R(2)/r1

        p21=(p2-p1)/(y-1)
        p11=(d1-g*p1)/2
        p22=(d2-g*p2)/(2*y)
        p2111=(p21-p11)/(y-1)
        p2211=(p22-p11)/(2*(y-1))

        c3=2*(p2211-p2111)/(y-1)
        c2=p2111-(y+2)*c3
        c1=p21-(y+1)*c2-(y*y+y+1)*c3
        c0=p1-c1-c2-c3

        P(ii+k+5) =c0
        P(ii+k+7) =c1
        P(ii+k+9) =c2
        P(ii+k+11)=c3

        rk=(R(2)/R(1))**k
        pr1=(c0+c1+c2+c3)*R(1)**gam / P(1)               !#
        dr1=(g*c0+(g+2)*c1+(g+4)*c2                      !#
     >     +(g+6)*c3)*R(1)**(gam-1) / dp1                !# all these parameters
        pr2=rk*(c0+c1*y+c2*y*y+c3*y**3)*R(2)**gam / P(2) !# should be equal to 1
        dr2=rk*(g*c0+(g+2)*c1*y+(g+4)*c2*y*y             !#
     >     +(g+6)*c3*y**3)*R(2)**(gam-1) / dp2           !#

        ier=1.d6*(dabs(pr1-1.d0)+dabs(dr1-1.d0)
     >     +dabs(pr2-1.d0)+dabs(dr2-1.d0))
        if (ier.GE.1) then
          write( *,5) gam,kap,pr1,dr1,pr2,dr2
          write(11,5) gam,kap,pr1,dr1,pr2,dr2
 5        format(4X,'Expansion error for gam = ',F5.1,' kap = ',I2,
     >          /4X,'pr1 = ',E15.7,' dr1 = ',E15.7
     >          /4X,'pr2 = ',E15.7,' dr2 = ',E15.7
     >          /4X,'shall we quit? ')
          read(*,*) iyes
          if (iyes.NE.0) stop
          Ierr=Ierr+ier
        end if

       return
      end
C     =================================================
      subroutine Origin3(P,gam,kap)
C     >> variant with three terms of Taylor expansion <<
C     =================================================
      implicit real*8 (a-h,o-z)
      include "hfd.par"                       !# 27/20/06:
       common /Ii/Ii/MaxT/MaxT/Ierr/Ierr/H/H  !# Instead of matching P(1),P(2),P(3)
       common /R/R(IP6)/V/V(IP6)              !# we now match P(1),P(2),dP/dr(1)
       dimension P(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (MaxT.EQ.0) MaxT=9        !### hfd uses Nmax instead of MaxT
        do i=ii+5,ii+5+MaxT
          P(i)=0.d0
        end do

        if (P(1).EQ.0.d0) return

        if (kap.LT.0) then           !### for positive kappa expansion
          k=0                        !#### goes over odd powers
          rn=1.d0
        else
          k=1
          rn=R(1)
        end if
        g=gam+k

        dpr =(-147*P(1)+360*P(2)           !# 7-node right derivative
     >       -450*P(3)+400*P(4)-225*P(5)   !# to be matched with
     >       +72*P(6)-10*P(7))/(60*h*V(1)) !# Taylor expansion


 1      f1=P(1)*rn/R(1)**g     !# Taylor expansion should match
        f2=P(2)*rn/R(2)**g     !# P(1), P(2), and dP/dr(1)
        f3=dpr *rn/R(1)**(g-1) !# (see 26/10/06)
        r1=R(1)*R(1)
        r2=R(2)*R(2)/r1

        f21=(f2-f1)/(r2-1.d0)
        f31=(f3-g*f1)/(r2-1.d0)/2
        f321=f21/(r2-1.d0)-f31
        c2=f321
        c1=f21-c2*(1.d0+r2)
        c0=f1-c1-c2

        P(ii+k+5)=c0
        P(ii+k+7)=c1
        P(ii+k+9)=c2

        pr1=(c0+c1+c2)*R(1)**gam                   !# ?=? P(1)
        dpl=(g*c0+(g+2)*c1+(g+4)*c2)*R(1)**(gam-1) !# ?=? dP/dr(1)
        rp =P(1)/pr1
        rd =dpr/dpl

        ier=1.d6*(dabs(rp-1.d0)+dabs(rd-1.d0))
        if (ier.GE.1) then
          write( *,5) gam,kap,pr1,P(1),dpl,dpr
          write(11,5) gam,kap,pr1,P(1),dpl,dpr
 5        format(4X,'Expansion error for gam = ',F5.1,' kap = ',I2,
     >          /4X,'P(1) = ',E15.7,' instead of ',E15.7
     >          /4X,'P`(1)= ',E15.7,' instead of ',E15.7
     >          /4X,'shall we quit? ')
          read(*,*) iyes
          if (iyes.NE.0) stop
          Ierr=Ierr+ier
        end if

       return
      end
C     =================================================
      subroutine Dif5(P,CP,k)             ! 5 node derivative
      implicit real*8 (a-h,o-z)           ! without Taylor expansion
      include "hfd.par"                   ! first k nodes are skipped
       common /Ii/Ii/kt/kt/H/H
     >        /V/V(IP6)/R/R(IP6)
       dimension P(IP6),CP(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        ih=2-kt
        i0=1+k*ih
        h12=12*ih*h

        ix=ii                      !### P(ix) is last nonzero value
        do j=1,ii,ih
          i1=ii-j+1
          if (P(i1).NE.0.d0) goto 200
          ix=i1-ih
          CP(i1)=0.d0
        end do

 200    continue

c>>>>> Symmetric 5-node differentiation for internal nodes:

        ih2=2*ih
        ih3=3*ih

        do i=ih2+i0,ix-ih2,ih
          CP(i)=(-P(i+ih2)+P(i-ih2)
     >         +8*(P(i+ih)-P(i-ih)))
     >          /(h12*V(i))
        end do

c>>>>> Last two nodes. Extrapolation by (ix-i)**2*(a+b*(i+3-ix)):

        CP(ix-ih)=CP(ix-ih2)/2-CP(ix-ih3)/9
        CP(ix)=0.d0

c>>>>> Short distances are not treated here:
       do i=1,i0+ih,ih
        CP(i)=0.d0
       end do
       return
      end
C     =================================================
      subroutine Cut_Short(P)         ! Smoothly brings P to
       implicit real*8 (a-h,o-z)      ! zero between r_1 and r_k.
       include "hfd.par"              ! Anzatz used:
       common /Ii/Ii/H/H/kt/kt/MaxT/MaxT ! (r-r1)**2*[a(r-rk)**2+b(r-rk)+c]
       common /R/R(IP6)/V/V(IP6)
       dimension P(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        ih=2-kt
        h12=12*ih*h
        ih2=2*ih
        ih3=3*ih

        k=15                          ! defines grid node r_k

        if (MaxT.EQ.0) MaxT=9         ! hfd uses Nmax instead of MaxT
        do i=ii+4,ii+5+MaxT
          P(i)=0.d0
        end do

        kp=k+ih
        km=k-ih
        pk=P(k)
        dpk=(-(P(k+ih2)-P(k-ih2))+8*(P(kp)-P(km))) ! first derivative at r_k
     >       /(h12*V(k))
        ddpk=(-(P(k+ih2)+P(k-ih2))+16*(P(kp)+P(km))-30*pk)
     >       /(h12*V(k))**2                        ! second derivative at r_k

        r1=R(1)
        rk=R(k)
        dr = rk-r1
        dr2=dr*dr
        c  = pk/dr2                                 ! parameters of anzatz
        b  = dpk/dr2 - 2*c/dr                       ! which preserves C_2
        a  = 0.5d0*ddpk/dr2 - c/dr2 - 2*b/dr        ! smoothness

        do i=1,k,ih
          ri=R(i)
          P(i)=(ri-r1)**2 * (a*(ri-rk)**2 + b*(ri-rk) + c)
        end do

        dpk1=(-(P(k+ih2)-P(k-ih2))+8*(P(kp)-P(km))) ! new value for P'(rk)
     >       /(h12*V(k))

       return
      end