c     ========================================================
c     Parameter IPs in rpa.par must be the same as in conf.par
c     ========================================================
      Subroutine Check_Dim (itest)
      include "conf.par"
       if (itest.NE.IPs) then
         write(*,5) itest,IPs
 5       format(//5x,'Wrong value of parameter IPs in rpa.par:',
     >          ' IPs=',i4,/5x,'In conf.par IPs=',i4,
     >          /5x,'Push...')
         read(*,*)
         stop
       end if
       return
      end
c     ========================================================