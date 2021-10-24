Module solvers

    Use, Intrinsic :: iso_fortran_env

    Implicit None

    Private

    Public :: DECOMP, FLSOLV
    Integer, Parameter  :: dp   = real64 

  Contains

    Subroutine DECOMP(A,N,SCALES,MPS)
        Implicit None
        Integer :: I, IP, IP1, IKK, IDXPIV, JK, J, N1, K, KK, KP, KP1, N
        Real(dp) :: ROWNRM, RR, BIG, SIZE, PIVOT, EM
        Integer, Allocatable, Dimension(:) :: MPS
        Real(dp), Allocatable, Dimension(:) :: A, SCALES
        Do I=1,N
            MPS(I)=I
            ROWNRM=0.d0
            JK=0
            Do J=1,N
                RR=dABS(A(I+JK))
                If(ROWNRM.LT.RR)ROWNRM=RR
                JK=JK+N
            End Do
            SCALES(I)=1.d0/ROWNRM
        End Do
        ! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
        KK=-N
        N1=N-1
        Do K=1,N1
            KK=KK+N
            BIG=0.d0
            Do I=K,N
                IP=MPS(I)
                SIZE=dABS(A(KK+IP))*SCALES(IP)
                If(SIZE.LE.BIG) Cycle
                BIG=SIZE
                IDXPIV=I
            End Do
            If (IDXPIV.NE.K) Then
                J=MPS(K)
                MPS(K)=MPS(IDXPIV)
                MPS(IDXPIV)=J
            Else
                Continue
            End If
        
            KP=MPS(K)
            PIVOT=-1.d0/A(KP+KK)
            KP1=K+1
            Do I=KP1,N
              IP=MPS(I)
              JK=IP+KK
              EM=A(JK)*PIVOT
              A(JK)=-EM
              JK=N+KK+KP
              IP1=IP-KP
              Do J=KP1,N
                IKK=JK+IP1
                A(IKK)=A(IKK)+EM*A(JK)
                JK=JK+N
              End Do
              JK=JK+N
            End Do
            JK=JK+N
        End Do
        Return
    End Subroutine DECOMP
    
    Subroutine FLSOLV(NB,A,B,X,MPS)
        Implicit None
        Integer :: NB, NP1, I, IP, IP1, IM1, J, JK, IBACK
        Real(dp) :: SUM
        Integer, Allocatable, Dimension(:) :: MPS
        Real(dp), Allocatable, Dimension(:) :: A, B, X
        NP1=NB+1
        IP=MPS(1)
        X(1)=B(IP)
        Do I=2,NB
            IP=MPS(I)
            IM1=I-1
            SUM=0.d0
            JK=0
            Do J=1,IM1
                SUM=SUM+A(IP+JK)*X(J)
                JK=JK+NB
            End Do
            X(I)=B(IP)-SUM
        End Do
        IP=MPS(NB)
        X(NB)=X(NB)/A(IP+(NB-1)*NB)
        
        Do IBACK=2,NB
            I=NP1-IBACK
            IP=MPS(I)
            IP1=I+1
            SUM=0.d0
            JK=I*NB
            Do J=IP1,NB
                SUM=SUM+A(IP+JK)*X(J)
                JK=JK+NB
            End Do
            X(I)=(X(I)-SUM)/A(IP+(I-1)*NB)
        End Do
        Return
    End Subroutine FLSOLV
End Module solvers