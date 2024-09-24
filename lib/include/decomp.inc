C     ==============================================================
C     Subroutines DECOMP and FLSOLV from V.Dzuba (changed to real*8)
C     ==============================================================
      SUBROUTINE DECOMP(NA,A,N,SCALES,MPS)
      implicit real*8 (a-h,o-z)
      DIMENSION A(NA)
      DIMENSION MPS(N),SCALES(N)
      DO 5 I=1,N
      MPS(I)=I
      ROWNRM=0.d0
      JK=0
      DO 2 J=1,N
      RR=dABS(A(I+JK))
      IF(ROWNRM.LT.RR)ROWNRM=RR
2     JK=JK+N
5     SCALES(I)=1.d0/ROWNRM
C            GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      KK=-N
      N1=N-1
      DO 16 K=1,N1
        KK=KK+N
        BIG=0.d0
        DO 11 I=K,N
          IP=MPS(I)
          SIZE=dABS(A(KK+IP))*SCALES(IP)
          IF(SIZE.LE.BIG)GO TO 11
          BIG=SIZE
          IDXPIV=I
11      CONTINUE
        IF(IDXPIV.EQ.K)GO TO 15
        J=MPS(K)
        MPS(K)=MPS(IDXPIV)
        MPS(IDXPIV)=J
15      KP=MPS(K)
        PIVOT=-1.d0/A(KP+KK)
        KP1=K+1
        DO 16 I=KP1,N
          IP=MPS(I)
          JK=IP+KK
          EM=A(JK)*PIVOT
          A(JK)=-EM
          JK=N+KK+KP
          IP1=IP-KP
          DO 16 J=KP1,N
            IKK=JK+IP1
            A(IKK)=A(IKK)+EM*A(JK)
16    JK=JK+N
      RETURN
      END
C     =================================================
      SUBROUTINE FLSOLV(NA,NB,A,B,X,MPS)
      implicit real*8 (a-h,o-z)
      DIMENSION A(NA),B(NB),X(NB)
      DIMENSION MPS(NB)
      NP1=NB+1
      IP=MPS(1)
      X(1)=B(IP)
      DO 2 I=2,NB
      IP=MPS(I)
      IM1=I-1
      SUM=0.d0
      JK=0
      DO 1 J=1,IM1
      SUM=SUM+A(IP+JK)*X(J)
1     JK=JK+NB
2     X(I)=B(IP)-SUM
      IP=MPS(NB)
      X(NB)=X(NB)/A(IP+(NB-1)*NB)
      DO 4 IBACK=2,NB
      I=NP1-IBACK
      IP=MPS(I)
      IP1=I+1
      SUM=0.d0
      JK=I*NB
      DO 3 J=IP1,NB
      SUM=SUM+A(IP+JK)*X(J)
3     JK=JK+NB
4     X(I)=(X(I)-SUM)/A(IP+(I-1)*NB)
      RETURN
      END