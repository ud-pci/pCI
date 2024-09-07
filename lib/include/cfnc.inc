C     ======================================== 22/11/05 =====
c     || Functions Norb, Kph, CFnc, & D_inf for MBPT codes  ||
C     =======================================================
      Function Norb (n,l,j) !### converts (n,l,j) => Norb
       implicit real*8 (a-h,j,o-z)
      INCLUDE "conf.par"
       common /Ns/Ns
       common /Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        k=2*j+0.1d0
        do i=1,Ns
          if (n.EQ.Nn(i).AND.l.EQ.Ll(i).AND.k.EQ.Lj(i)) then
            Norb=i
            return
          end if
        end do
        Norb=0
       return
      end
C     =================================================
      Function Kph(xk,ip)     !### for Coulomb photon ip=(-1)**xk
                              !### if that is not true, xk => xk+1
       implicit real*8 (a-h,j,o-z)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        k=dabs(xk)+0.1d0
        if (ip.EQ.Isig(k)) then
          kp=k
        else
          kp=k+1
        end if
        Kph=kp
       return
      end
C     =================================================
      Subroutine CFnc(kq,ja,la,jb,lb,pa,qa,pb,qb,c1) !### evaluates Yk
      implicit real*8 (a-h,o-z)
      include "hfd.par"
       common /Ii/Ii/Kt/Kt/C_SMS/C_SMS
     >        /C/C(IP6)
       dimension c1(IP6),pa(IP6),qa(IP6),pb(IP6),qb(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        ih=2-Kt
        do i=1,Ii,ih
          C(i)=pa(i)*pb(i)+qa(i)*qb(i)
        end do
        C(Ii+4)=pa(Ii+4)+pb(Ii+4)
        call Yk(kq)
        do i=1,IP6
          c1(i)=C(i)
        end do
        if (C_SMS.NE.0.d0.AND.kq.EQ.1) then
          c1(Ii+1)=P_eff(ja,la,jb,lb,pa,qa,pb,qb)
        else
          c1(Ii+1)=0.d0
        end if
       return
      end
C     =================================================
      Function D_inf(d1, ! d1=D(Lix)   # extrapolation
     >               d2, ! d2=D(Lix-1) # to D(infty)
     >               d3) ! d3=D(Lix-2) # using (l+a)^{-3} law
      implicit real*8 (a-h,o-z)
      include "hfd.par"
       common /Lix/Lix/Kout/Kout/K_inf/K_inf
       dimension x(3)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        err=1.d5
        d0=d1
        a0=-9.d0
        imax=K_inf           ! variants: 1,3,5. imax=1,=>a=0;
        i1=(imax+1)/2        ! imax=3,5,=> a in [0.5,0,-0.5] or [1,.5,0,-.5,-1]
        do i=1,imax
          a=0.5d0*(i1-i)     ! extrapolation expression: D_l = D_inf -A/(l+a)^3
          if (Lix+a.LT.2.99d0) goto 100
          do k=1,3
            l=Lix+1-k
            x(k)=1.d0/(l+a)**3
          end do
          d12=(x(2)*d1-x(1)*d2)/(x(2)-x(1))
          d13=(x(3)*d1-x(1)*d3)/(x(3)-x(1))
          er=dabs(d12-d13)/(dabs(d12)+dabs(d13))
          if (er.LT.err) then
            err=er
            d0=d12
            a0=a
          end if
        end do
 100    D_inf=d0
        if (Kout.GE.2) write(11,5) d3,d2,d1,d0,a0
 5      format(4X,'D_inf: ',4E14.6,' a =',F4.1)
       return
      end
C     =================================================
      Function R_Lmax(str,! name of the diagram
     >                r1, ! r1=R(Lix)   # extrapolation R(Lix)
     >                r2, ! r2=R(Lix-1) # to R(Lkx) assuming that
     >                r3) ! r3=R(Lix-2) # R(L+1)-R(L)=C/(L+1)^s
      implicit real*8 (a-h,o-z)
      include "hfd.par"
       common /Lix/Lix/Lkx/Lkx/Kout/Kout/N_inf/N_inf
       character*3 str
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        small=1.d-10
        R_Lmax=r1                          !### Lkx = required L_max,
        if(Lix.GE.Lkx.OR.Lix.LE.1) return  !#### Lix = actual L_max

        if((r1-r2)*(r2-r3).LT.small) then  !### If R(L) oscillates,
          N_inf=N_inf+1                    !#### no extrapolation is done
          return
        end if

        x=Lix
        a=dlog((x-1.d0)/x)
        s=dlog((r1-r2)/(r2-r3))/a
        c=(r1-r2)*x**s

        rx=r1                           !### extrapolated value
        do l=Lix+1,Lkx
          x=l
          rx=rx+c/x**s
        end do

        R_Lmax=rx
        if (Kout.GE.2) write(*,5) str,r3,r2,r1,rx,s
 5      format(4X,'Extrapolation for ',A3,/4X,4E13.6,';  s =',F6.2)
       return
      end
