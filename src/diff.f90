Module diff
    ! This module contains subroutines for numerical differentiation
    Use basc_variables

    Implicit None

    Private

    Public :: Dif, Dif5, Origin, Cut_Short

  Contains

    Subroutine Dif(P,CP,gm)          
        ! This subroutine should be used for functions which are not very smooth at R=R(1).
        Implicit None

        Integer :: i, ih, im, ix, j, i1, ih2, ih3, ih4, ih5, ih6, ierr, kap
        Real(dp) :: h60, gm, r1, cp1, p1, pim, scale, err
        Real(dp), Dimension(IP6) :: P, CP                 !

        If (MaxT.EQ.0) MaxT=9      !### hfd uses Nmax instead of MaxT
        ih=2-kt
        h60=60*ih*h
        CP(ii+4)=gm-1.d0

        !>>>>>> Expansion at the origin:
        r1=R(1)
        cp1=0.d0                   !### = left derivative (dP/dr)_{-}
        p1=0.d0                    !### = P(r1)
        do m=0,MaxT
            im=ii+5+m
            pim=P(im)
            CP(im)=(gm+m)*pim
            cp1=cp1+CP(im)
            p1=p1+pim
        End do

        p1=p1*r1**gm
        If (dabs(p1).GT.1.d-19) Then
            scale=P(1)/p1
        Else
            scale=1.d0
        End If

        If(dabs(scale-1.d0).GT.1.d-3) Then
            Write(*,*) ' Dif: expansion does not match the function'
            Write(*,*) ' P(1)=',P(1),' Taylor=',p1
            Read(*,*)
            Stop
        End If

        cp1=scale*cp1*r1**(gm-1.d0)      !### = CP(r1)
        p1=scale*p1                      !### = P(r1)
        !<<<<<< Expansion at the origin

        ix=ii                      !### P(ix) is last nonzero value
        do j=1,ii,ih
            i1=ii-j+1
            If (P(i1).NE.0.d0) Exit
            ix=i1-ih
            CP(i1)=0.d0
        End do

        !>>>>> Symmetric 7-node differentiation for internal nodes:
        ih2=ih +ih
        ih3=ih2+ih
        ih4=ih3+ih

        do i=ih3+1,ix-ih3,ih 
          CP(i)=(P(i+ih3)-9*(P(i+ih2)-P(i-ih2)) &
               +45*(P(i+ih)-P(i-ih))-P(i-ih3)) &
                /(h60*V(i))
        End do

        !>>>>> Asymmetric 7-node differentiation for first three nodes:
        ih5=ih4+ih
        ih6=ih5+ih

        CP(1)=(-147*P(1)+360*P(1+ih) &
             -450*P(1+ih2)+400*P(1+ih3)-225*P(1+ih4) &
             +72*P(1+ih5)-10*P(1+ih6))/(h60*V(1))

        CP(1+ih)=(-10*P(1)-77*P(1+ih) &
             +150*P(1+ih2)-100*P(1+ih3)+50*P(1+ih4) &
             -15*P(1+ih5)+2*P(1+ih6))/(h60*V(1+ih))

        CP(1+ih2)=(2*P(1)-24*P(1+ih)-35*P(1+ih2) &
             +80*P(1+ih3)-30*P(1+ih4)+8*P(1+ih5) &
             -P(1+ih6))/(h60*V(1+ih2))


        !>>>>> Asymmetric 7-node differentiation for last three nodes:
        CP(ix)=-(-147*P(ix)+360*P(ix-ih) &
             -450*P(ix-ih2)+400*P(ix-ih3)-225*P(ix-ih4) &
             +72*P(ix-ih5)-10*P(ix-ih6))/(h60*V(ix))

        CP(ix-ih)=-(-10*P(ix)-77*P(ix-ih) &
             +150*P(ix-ih2)-100*P(ix-ih3)+50*P(ix-ih4) &
             -15*P(ix-ih5)+2*P(ix-ih6))/(h60*V(ix-ih))

        CP(ix-ih2)=-(2*P(ix)-24*P(ix-ih)-35*P(ix-ih2) &
             +80*P(ix-ih3)-30*P(ix-ih4)+8*P(ix-ih5) &
             -P(ix-ih6))/(h60*V(ix-ih2))

        !                                              ! We assume that this
        ! If (cp1.EQ.0.d0.AND.P(1).EQ.0.d0) CP(1)=0.d0 ! may occur only for
        !                                              ! B-spline.
        err=cp1-CP(1)
        ierr=0
        If (dabs(cp1).GT.1.d0) err=err/cp1
        If (dabs(err).GT.1.d-1) Then
            Write(*,'(4X,"Dif: matching error at R1. Taylor:",E12.5, &
                   ", Num:",E12.5,/4X,"Using Origin for new expansion.")') cp1,CP(1)
            Ierr=Ierr+10*dabs(err)
            If (dabs(CP(ii+5)).GT.dabs(CP(ii+6))) Then
                kap=1                            ! For kap>0 expansion
            Else                                 ! goes over even powers
                kap=-1
            End If
            call Origin(CP,CP(ii+4),kap)
            Read(*,*)
        End If

        Return
    End Subroutine Dif

    Subroutine Dif5(P,CP,k)             
        ! 5 node derivative without Taylor expansion, first k nodes are skipped
        Implicit None           
        
        Integer :: ih, i0, i1, ix, k, j, ih2, ih3, i
        Real(dp) :: h12
        Real(dp), Dimension(IP6) :: P, CP

        ih=2-kt
        i0=1+k*ih
        h12=12*ih*h
  
        ix=ii                      !### P(ix) is last nonzero value
        Do j=1,ii,ih
            i1=ii-j+1
            if (P(i1).NE.0.d0) Exit
            ix=i1-ih
            CP(i1)=0.d0
        End Do
    
        ! Symmetric 5-node differentiation for internal nodes:
        ih2=2*ih
        ih3=3*ih
        Do i=ih2+i0,ix-ih2,ih
            CP(i)=(-P(i+ih2)+P(i-ih2)+8*(P(i+ih)-P(i-ih)))/(h12*V(i))
        End Do
  
        ! Last two nodes. Extrapolation by (ix-i)**2*(a+b*(i+3-ix)):
        CP(ix-ih)=CP(ix-ih2)/2-CP(ix-ih3)/9
        CP(ix)=0.d0
  
        ! Short distances are not treated here:
        Do i=1,i0+ih,ih
            CP(i)=0.d0
        End Do

        Return
    End Subroutine Dif5

    Subroutine Origin(P,gam,kap)           ! 19/05/08
        !     >> variant with four terms of Taylor expansion <<
        !# we now match P(1),P(2), dP/dr(1), and dP/dr(2)
        Implicit None
        Integer :: k, kap, ier, iyes, ierr
        Real(dp) :: gam, rn, g, dp1, dp2, p1, p2, d2, r1, y, p21, p11, d1
        Real(dp) :: p22, p2111, p2211, c2, c3, c1, c0, rk, pr1, dr1, pr2, dr2
        Real(dp), Dimension(IP6) :: P
        Character(Len=512) :: strfmt

        If (MaxT.EQ.0) MaxT=9        !### hfd uses Nmax instead of MaxT

        P(ii+5:ii+5+MaxT)=0.d0

        If (P(1).EQ.0.d0) Return

        If (kap.LT.0) Then           !### for positive kappa expansion
            k=0                      !#### goes over odd powers
            rn=1.d0
        Else
            k=1
            rn=R(1)
        End If

        g=gam+k

        dp1 =(-147*P(1)+360*P(2) &            !# 7-node right derivative
             -450*P(3)+400*P(4)-225*P(5) &    !# to be matched with
             +72*P(6)-10*P(7))/(60*h*V(1))    !# Taylor expansion
        dp2 =(-10*P(1)-77*P(2) &
             +150*P(3)-100*P(4)+50*P(5) &
             -15*P(6)+2*P(7))/(60*h*V(2))


        p1=P(1)*rn/R(1)**g     !# (see 26/10/06)
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
        dr1=(g*c0+(g+2)*c1+(g+4)*c2+(g+6)*c3)*R(1)**(gam-1) / dp1                !# all these parameters
        pr2=rk*(c0+c1*y+c2*y*y+c3*y**3)*R(2)**gam / P(2) !# should be equal to 1
        dr2=rk*(g*c0+(g+2)*c1*y+(g+4)*c2*y*y+(g+6)*c3*y**3)*R(2)**(gam-1) / dp2    

        ier=1.d6*(dabs(pr1-1.d0)+dabs(dr1-1.d0)+dabs(pr2-1.d0)+dabs(dr2-1.d0))
        If (ier.GE.1) Then
            strfmt='(4X,"Expansion error for gam = ",F5.1," kap = ",I2, &
                /4X,"pr1 = ",E15.7," dr1 = ",E15.7 &
                /4X,"pr2 = ",E15.7," dr2 = ",E15.7 &
                /4X,"shall we quit? ")'
            Write( *,strfmt) gam,kap,pr1,dr1,pr2,dr2
            Write(11,strfmt) gam,kap,pr1,dr1,pr2,dr2
            Read(*,*) iyes
            If (iyes.NE.0) Stop
            Ierr=Ierr+ier
        End If

        Return
    End Subroutine Origin

    subroutine Cut_Short(P)         
        ! Smoothly brings P to zero between r_1 and r_k.
        ! Anzatz used: (r-r1)**2*[a(r-rk)**2+b(r-rk)+c]
        Implicit None

        Integer :: ih, ih2, ih3, k, i, kp, km
        Real(dp) :: h12, pk, dpk, ddpk, r1, ri, rk, dr, dr2, a, b, c, dpk1
        Real(dp), Dimension(IP6) :: P

        ih=2-kt
        h12=12*ih*h
        ih2=2*ih
        ih3=3*ih
 
        k=15                          ! defines grid node r_k
        if (MaxT.EQ.0) MaxT=9         ! hfd uses Nmax instead of MaxT
        Do i=ii+4,ii+5+MaxT
            P(i)=0.d0
        End Do
 
        kp=k+ih
        km=k-ih
        pk=P(k)
        ! first derivative at r_k
        dpk=(-(P(k+ih2)-P(k-ih2))+8*(P(kp)-P(km)))/(h12*V(k))
        ! second derivative at r_k
        ddpk=(-(P(k+ih2)+P(k-ih2))+16*(P(kp)+P(km))-30*pk)/(h12*V(k))**2                        
 
        r1=R(1)
        rk=R(k)
        dr = rk-r1
        dr2=dr*dr
        c  = pk/dr2                                 ! parameters of anzatz
        b  = dpk/dr2 - 2*c/dr                       ! which preserves C_2
        a  = 0.5d0*ddpk/dr2 - c/dr2 - 2*b/dr        ! smoothness
 
        Do i=1,k,ih
            ri=R(i)
            P(i)=(ri-r1)**2 * (a*(ri-rk)**2 + b*(ri-rk) + c)
        End Do
        ! new value for P'(rk)
        dpk1=(-(P(k+ih2)-P(k-ih2))+8*(P(kp)-P(km)))/(h12*V(k))
 
        Return
    End Subroutine Cut_Short

End Module diff