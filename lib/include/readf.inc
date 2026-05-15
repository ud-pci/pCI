c     =================================================
c     subroutine ReadF and function Isig
c     =================================================
      Subroutine ReadF (Kan,record,V1,V2,nrec)
      include "hfd.par"
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C        FORTRAN-77, MS-DOS VERSION
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       integer record,nrec
       real*8 V1(IP6),V2(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        read(Kan,rec=nr1) (V1(I),i=1,ii)
        if (nrec.EQ.2) read(Kan,rec=nr2) (V2(I),i=1,ii)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       return
      end
C     =================================================
      Function Isig(n)
c     Isig=(-1)**n
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (n.EQ.2*(n/2)) then
           Isig=1
        else
           Isig=-1
        end if
       return
      end