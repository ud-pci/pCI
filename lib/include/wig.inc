      FUNCTION FJ9(A11,A12,A13,A21,A22,A23,A31,A32,A33)
         implicit real*8 (a-h,o-z)
      DIMENSION UP(3),DOUN(3)
      CALL DELFJ(A11,A12,A13,N)
      GO TO(10,3),N
10    CALL DELFJ(A21,A22,A23,N)
      GO TO(11,3),N
11    CALL DELFJ(A31,A32,A33,N)
      GO TO(12,3),N
12    CALL DELFJ(A11,A21,A31,N)
      GO TO(13,3),N
13    CALL DELFJ(A12,A22,A32,N)
      GO TO(14,3),N
14    CALL DELFJ(A13,A23,A33,N)
      GO TO(15,3),N
15    UP(1)=A11+A33
      DOUN(1)=dABS(A11-A33)
      UP(2)=A32+A21
      DOUN(2)=dABS(A32-A21)
      UP(3)=A12+A23
      DOUN(3)=dABS(A12-A23)
      X=UP(1)
      DO 1 N=2,3
      IF(UP(N)-X) 6,1,1
6     X=UP(N)
1     CONTINUE
      Y=DOUN(1)
      DO 2 N=2,3
      IF(DOUN(N)-Y) 2,2,7
7     Y=DOUN(N)
2     CONTINUE
      IF(X-Y)3,8,8
8     D=X-Y+1.d0
      N=D
      IF(dABS(D-N)-0.00001d0) 9,3,3
9     FJ9=0.d0
      DO 4 K=1,N
      A=Y+K-1.d0
      M=2*A+0.00001d0
      FJ9=FJ9+(-1)**M*(M+1)*FJ6(A11,A21,A31,A32,A33,A)*
     *FJ6(A12,A22,A32,A21,A,A23)*FJ6(A13,A23,A33,A,A11,A12)
4     CONTINUE
      GO TO 5
3     FJ9=0.d0
5     RETURN
      END
c     =====================================================
      SUBROUTINE DELFJ(A,B,C,N)
         implicit real*8 (a-h,o-z)
C     N=2 MEANS DELFJ=0.d0
      X=A+B
      Y=dABS(A-B)
      IF(X-C) 1,2,2
2     IF(Y-C) 3,3,1
3     N=1
      GO TO 4
1     N=2
4     RETURN
      END

C     =================================================
      FUNCTION FJ3(A1,A2,A3,A4,A5,A6)
         implicit real*8 (a-h,o-z)
      COMMON/ZallFJ/ALL
      DIMENSION X(9),LU(9),LI(9),Y(3),MA(3),Z(3),NA(3)
C      WRITE(6,100) A1,A2,A3,A4,A5,A6
100    FORMAT(5X,6(F5.1,2X),E15.5)
      IF(dABS(A4+A5+A6)-0.0001d0) 20,2,2
20    X(1)=A1+A2-A3
      X(2)=A1-A2+A3
      X(3)=-A1+A2+A3
      X(4)=A1+A4
      X(5)=A1-A4
      X(6)=A2+A5
      X(7)=A2-A5
      X(8)=A3+A6
      X(9)=A3-A6
      EPS=0.0001d0
      DO 1 K=1,9
      IF(X(K)) 2,3,3
3     LU(K)=X(K)
      IF(dABS(LU(K)-X(K))-EPS) 1,2,2
1     CONTINUE
      GO TO  4
2     FJ3=0.d0
      GO TO 5
4     Y(1)=A1+A2-A3
      Y(2)=A1-A4
      Y(3)=A2+A5
      Z(1)=0.d0
      Z(2)=-A3+A2-A4
      Z(3)=-A3+A1+A5
      DO 6 K=1,3
      MA(K)=Y(K)
      NA(K)=Z(K)
      IF(dABS(MA(K)-Y(K))-EPS) 7,2,2
7     IF(dABS(NA(K)-Z(K))-EPS) 6,2,2
6     CONTINUE
      CALL PORAD(MA,3)
      CALL PORAD(NA,3)
      IF(MA(1)) 2,8,8
8     IF(MA(1)-NA(3)) 2,9,9
9     LI(1)=A1+A2+A3+1.d0
      DO 10 K=1,2
      LI(K+1)=MA(K+1)-MA(1)
      LI(K+3)=MA(K+1)-MA(1)
      LI(K+5)=NA(3)-NA(K)
10    LI(K+7)=NA(3)-NA(K)
      CALL PORAD(LU,9)
      CALL PORAD(LI,9)
      ALL=0.d0
      DO 11 K=1,9
      IF(LU(K)-LI(K)) 12,11,13
12    CALL NGFJ(-1,LU(K),LI(K))
      GO TO 11
13    CALL NGFJ(1,LI(K),LU(K))
11    CONTINUE
      FJ3=0.d0
      UN=ALL/2
      KM=MA(1)-NA(3)+1
      DO 14 KU=1,KM
      K=NA(3)+KU-1
      ALL=0.d0
      CALL NGFJ(-1,0,MA(1)-K)
      CALL NGFJ(-1,MA(2)-MA(1),MA(2)-K)
      CALL NGFJ(-1,MA(3)-MA(1),MA(3)-K)
      CALL NGFJ(-1,0,K-NA(3))
      CALL NGFJ(-1,NA(3)-NA(1),K-NA(1))
      CALL NGFJ(-1,NA(3)-NA(2),K-NA(2))
      L=K+A1-A2-A6
14    FJ3=FJ3+(-1)**L*dEXP(ALL)
      FJ3=dEXP(UN)*FJ3
5     CONTINUE
C     WRITE(2,100)A1,A2,A3,A4,A5,A6,FJ3
      RETURN
      END

C     =================================================
      FUNCTION FJ6(A1,A2,A3,A4,A5,A6)
         implicit real*8 (a-h,o-z)
      COMMON/ZallFJ/ALL
      DIMENSION MA(4),NA(3),X(12),LU(13),LI(13)
      X1=A1+A2+A3
      X2=A1+A5+A6
      X3=A4+A2+A6
      X4=A4+A5+A3
      Y1=A1+A2+A4+A5
      Y2=A2+A3+A5+A6
      Y3=A3+A1+A6+A4
      EPS=0.0001d0
      MA(1)=X1+EPS
      MA(2)=X2+EPS
      MA(3)=X3+EPS
      MA(4)=X4+EPS
      NA(1)=Y1+EPS
      NA(2)=Y2+EPS
      NA(3)=Y3+EPS
      IF(dABS(MA(1)-X1)-EPS) 20,12,12
20    IF(dABS(MA(2)-X2)-EPS) 21,12,12
21    IF(dABS(MA(3)-X3)-EPS) 22,12,12
22    IF(dABS(MA(4)-X4)-EPS) 23,12,12
23    IF(dABS(NA(1)-Y1)-EPS) 24,12,12
24    IF(dABS(NA(2)-Y2)-EPS) 25,12,12
25    IF(dABS(NA(3)-Y3)-EPS) 26,12,12
26    GO TO 10
12    FJ6=0.d0
      GO TO 11
10    CALL PORAD(MA,4)
      CALL PORAD(NA,3)
      IF(MA(4)-NA(1)) 27,27,12
27    X(1)=A1+A2-A3
      X(2)=A1-A2+A3
      X(3)=-A1+A2+A3
      X(4)=A1+A5-A6
      X(5)=A1-A5+A6
      X(6)=-A1+A5+A6
      X(7)=A4+A2-A6
      X(8)=A4-A2+A6
      X(9)=-A4+A2+A6
      X(10)=A4+A5-A3
      X(11)=A4-A5+A3
      X(12)=-A4+A5+A3
      ALL=0.d0
      DO 1 K=1,12
      IF(X(K)) 12,1,1
1     LU(K)=X(K)+EPS
      LU(13)=MA(4)+1
      DO 2 K=1,3
      LI(K)=MA(K)+1
      LI(3+K)=MA(4)-MA(K)
2     LI(6+K)=LI(3+K)
      DO 3 K=1,2
      LI(9+K)=NA(K+1)-NA(1)
3     LI(11+K)=LI(9+K)
      CALL PORAD(LU,13)
      CALL PORAD(LI,13)
      ALL=0.d0
      DO 4 K=1,13
      IF(LU(K)-LI(K)) 5,4,7
5     CALL NGFJ(-1,LU(K),LI(K))
      GO  TO  4
7     CALL NGFJ(1,LI(K),LU(K))
4     CONTINUE
      FJ6=0.d0
      UN=ALL/2
      KM=NA(1)-MA(4)+1
      DO 6 L=1,KM
      K=MA(4)+L-1
      ALL=0.d0
      CALL NGFJ(1,MA(4)+1,K+1)
      CALL NGFJ(-1,MA(4)-MA(1),K-MA(1))
      CALL NGFJ(-1,MA(4)-MA(2),K-MA(2))
      CALL NGFJ(-1,MA(4)-MA(3),K-MA(3))
      CALL NGFJ(-1,0,K-MA(4))
      CALL NGFJ(-1,0,NA(1)-K)
      CALL NGFJ(-1,NA(2)-NA(1),NA(2)-K)
      CALL NGFJ(-1,NA(3)-NA(1),NA(3)-K)
      FJ6=FJ6+(-1)**K*dEXP(ALL)
6     CONTINUE
      FJ6=dEXP(UN)*FJ6
11    RETURN
      END

C     =================================================
      SUBROUTINE NGFJ(J,M,N)
         implicit real*8 (a-h,o-z)
      COMMON/ZallFJ/ ALL
      IF(M-N) 4,3,3
4     K=M+1
      DO 1 L=K,N
      AL=dLOG(1.d0*L)
1     ALL=ALL+J*AL
3     RETURN
      END

C     =================================================
      SUBROUTINE OBME(M1,M2)
         implicit real*8 (a-h,o-z)
      N=M1
      M1=M2
      M2=N
      RETURN
      END

C     =================================================
      SUBROUTINE PORAD(I,N)
         implicit real*8 (a-h,o-z)
      DIMENSION I(N)
      K1=N-1
      DO 1 K=1,K1
      K2=K+1
      DO 2 L=K2,N
      IF(I(K)-I(L)) 2,2,3
3     CALL OBME(I(K),I(L))
2     CONTINUE
1     CONTINUE
      RETURN
      END
C     ======================================== 18/01/04
C     sub. CalcWig0 tabulates Wig3J(j1,j2,k,1/2,-1/2,0)
C     for use in function Wig0.
C     =================================================
      subroutine CalcWig0
      implicit real*8 (a-h,o-z)
       common /In0/In0(168)/Gnt0/Gnt0(168)
        l=0
        do ij1=1,7
          y1=ij1-0.5d0
          do ij2=ij1,7
            y2=ij2-0.5d0
            kmin=ij2-ij1+1
            kmax=ij1+ij2
            do ki=kmin,kmax
              k=ki-1
              ind=100*k + 10*ij1 + ij2
              l=l+1
              In0(l)=ind
              xk=k
              x=FJ3(y1,y2,xk,0.5d0,-0.5d0,0.d0)
              Gnt0(l)=x
            end do
          end do
        end do
        if (l.NE.168) then
           write(*,*) ' Number of 3j-symbols is <> 168:',l
           stop
        end if
        write(*,*) ' CalcWig0: array Gnt0 formed'
       return
      end
C     =================================================
      Function Wig0(k,xj1,xj2)
       implicit real*8 (a-h,o-z)
       common /In0/In0(168)/Gnt0/Gnt0(168)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        y1=xj1
        y2=xj2
        g = 0.d0
        if (y2.LT.y1) then
           x = y2
           y2 = y1
           y1 = x
        end if
        ind = 100*k+10*y1+y2+5.51
        do i=1,168
           if(In0(i).EQ.ind) then
              g = Gnt0(i)
              goto 1000
           end if
        end do
c        write(*,*)' Wig0: can not find symbol for',ind
c        write(*,*)' k,j1,j2 =',k,y1,y2
c        read(*,*)
1000    Wig0 = g
       Return
      End