C     ======================================== 02/02/04
C     radial part of the second order diagrams
C     =================================================
      Function Rp2(ke,ja,la,jb,lb,c1,c2,pa,qa,pb,qb)
      implicit real*8 (a-h,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ii/Ii/Kt/Kt/C_SMS/C_SMS/Klow/Klow
     >        /C/C(IP6)/R/R(IP6)/Ll/Ll(IPs)/Lj/Lj(IPs)
       dimension c1(IP6),c2(IP6),pa(IP6),qa(IP6),pb(IP6),qb(IP6),
     >      pe(IP6),qe(IP6)
        ih=2-Kt
        je=Lj(ke)
        le=Ll(ke)
        call ReadFF(12,ke+4,pe,qe,2)
        do i=1,Ii,ih
          C(i)=c1(i)/R(i)*(pe(i)*pa(i)+qe(i)*qa(i))
        end do
        C(Ii+4)=c1(Ii+4)+pe(Ii+4)+pa(Ii+4)-1
        call Sint(r1)
        do i=1,Ii,ih
          C(i)=c2(i)/R(i)*(pe(i)*pb(i)+qe(i)*qb(i))
        end do
        C(Ii+4)=c2(Ii+4)+pe(Ii+4)+pb(Ii+4)-1
        call Sint(r2)

c       SMS contribution:
        s1=c1(Ii+1)
        s2=c2(Ii+1)
        if (s1.NE.0.d0) then
          r1=r1-C_SMS*s1*P_eff(ja,la,je,le,pa,qa,pe,qe)
        end if
        if (s2.NE.0.d0) then
          r2=r2-C_SMS*s2*P_eff(jb,lb,je,le,pb,qb,pe,qe)
        end if
        s=r1*r2
        Rp2=s
       return
      end