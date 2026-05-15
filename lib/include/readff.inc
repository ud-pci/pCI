C     =================================================
      Subroutine ReadFF (Kan,record,V1,V2,nrec)
      include "hfd.par"
      INCLUDE "conf.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Irdff/Ird
       common /ArrP/ArrP(IP6,IPs)/ArrQ/ArrQ(IP6,IPs)
       integer record,nrec
       dimension V1(IP6),V2(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Ird.NE.1) then
           if (Ns.GT.IPs) then
              write(*,*) ' ReadFF: Basis set is too long.'
              write(*,*) ' Max length:',IPs,' Ns=',Ns
              stop
           end if
           do ni=1,Ns
              call ReadF(Kan,ni+4,V1,V2,2)
              do i=1,IP6
                 ArrP(i,ni)=V1(i)
                 ArrQ(i,ni)=V2(i)
              end do
           end do
           Ird=1
        end if
        if (record.LE.4.OR.record.GT.Ns+4) then
           call ReadF (Kan,record,V1,V2,nrec)
        else
           ni=record-4
           do i=1,IP6
              V1(i)=ArrP(i,ni)
              V2(i)=ArrQ(i,ni)
           end do
        end if
       return
      end