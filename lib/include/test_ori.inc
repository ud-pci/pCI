C     =================================================
C     Use this subroutine to test Taylor expansion for
C     orbitals and functions in HFD.DAT format
c     sub. Test_origin1 deleted on 7/11/06. See older versions!
C     =================================================
      subroutine Test_Origin(str,P,err,ierr)  ! Test of the Taylor
      implicit real*8 (a-h,o-z)               ! expansion at the origin
      include "conf.par"
      include "hfd.par"
       common /ii/ii/MaxT/MaxT/Kout/Kout
       common /R/R(IP6)
       dimension P(IP6)
       character*1 str(10),chr

        ierr=0
        gm=P(ii+4)
        r1=R(1)
        r2=R(2)
        q=r2/r1
        q0=1.d0
        pp1=0.d0                      ! Expansion value for P(1)
        pp2=0.d0                      ! Expansion value for P(2)
        do k=0,MaxT
          ik=ii+5+k
          pp1=pp1+P(ik)
          pp2=pp2+q0*P(ik)
          q0=q0*q
        end do

        chr='.'
        if (pp1*P(1).GT.0.d0) then     ! Leading power is not always
          if (dabs(gm).LT.1.d-3.AND.   ! given by P(ii+4).
     >        dabs(pp1).GT.1.d-9) then ! In this case it is determined
            gm=dlog(P(1)/pp1)/dlog(r1) ! here, which is marked by
            chr='!'                    ! exclamation mark in the output
            igm=gm+0.1d0               !### This assumes finite nucleus,
            gm=igm                     !### when gamma is always integer
          end if
        end if

        pp1=r1**gm*pp1
        pp2=r2**gm*pp2
        err1=P(1)-pp1
        err2=P(2)-pp2
        if(dabs(pp1).GT.1.d-2) then
          err1=err1/(P(1)+pp1)
          err2=err2/(P(2)+pp2)
        end if
        if (dabs(err1).GT.err) then
          if (Kout.GT.0)  write(11,5) str
          if (Kout.GT.1)  write( *,5) str
 5        format(4X,'Error in Taylor expansion at the origin for ',
     >           10A1,':')
          if (Kout.GT.1)
     >      write(11,15) P(ii+4),(i,P(i),i=1,2),(m,P(ii+5+m),m=0,MaxT)
 15       format(2X,'Input:  gam=',F8.5,
     >          2(4X,'P(',I1,')=',E13.5),/2(4X,'C_',I1,' = ',E11.4))
          if (Kout.GT.0)  write(11,25) gm,chr,err1,err2
          if (Kout.GT.1)  write( *,25) gm,chr,err1,err2
 25       format(2X,'Output:  gm=',F8.5,A1,
     >           ' Relative error at r1 and r2',2E13.5)
          if (Kout.GT.1)  read(*,*)
          ierr=dabs(err1/err)
        end if
        P(ii+4)=gm
       return
      end
