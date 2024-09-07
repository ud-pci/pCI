C     =================================================
      SUBROUTINE SINT1(DS)        !### Simpson integration over r
      include "hfd.par"           !#### (with weight function HH*V(I))
       implicit real*8 (a-h,o-z)
       COMMON /II/II/KT/KT/H/H0/Dint/Dint
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Gam=C(ii+4)
        IH=2-KT
        HH=H0*IH/3.d0
        I0=IH+1
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        I1=1
        I2=I1+IH
        I3=I2+IH
        R1=R(I1)
        R2=R(I2)
        R3=R(I3)
        T1=0.d0
        T2=0.d0
        IF (GAM.GT.5.0d0) GOTO 200
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
C     - - - - - - - - - - - - - - - - - - - - - - - - -
200     Dint=t2
        P1=C(I0)*V(I0)
        P2=C( 1)*V( 1)
        I1=I0+IH
        DO 10 I=I1,II,IH
           Q=C(I)*V(I)
           T=HH*(Q+4.d0*P1+P2)
           T=T2+T
           T2=T1
           T1=T
           P2=P1
           P1=Q
10      CONTINUE
        ds=t
       RETURN
      END