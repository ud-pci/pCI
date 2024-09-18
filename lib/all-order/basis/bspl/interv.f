
      subroutine interv (xt,lxt,x,left,mflag)
      implicit doubleprecision(a-h,o-z)
      dimension xt(lxt)
      data ilo /1/
 
      ihi = ilo + 1
      if(ihi .lt. lxt)                go to 20
      if(x .ge. xt(lxt))              go to 110
      if(lxt . le. 1)                 go to 90
      ilo = lxt - 1
      ihi = lxt
 
 20   if(x .ge. xt(ihi))              go to 40
      if(x .ge. xt(ilo))              go to 100
 
      istep = 1
 31   continue
         ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)              go to 35
         if (x .ge. xt(ilo))          go to 50
         istep = istep*2
      go to 31
 35   ilo = 1
      if(x .lt. xt(1))                go to 90
                                      go to 50
 
 40   istep = 1
 41   continue
         ilo = ihi
         ihi = ilo + istep
         if(ihi .ge. lxt)             go to 45
         if(x .lt. xt(ihi))           go to 50
         istep = istep*2
      go to 41
 
 45   if (x .ge. xt(lxt))             go to 110
      ihi = lxt
 50   continue
            middle = (ilo + ihi)/2
            if(middle .eq. ilo)       go to 100
            if(x .lt. xt(middle))     go to 53
            ilo = middle
          go to 50
 53       ihi = middle
      go to 50
 
 90   mflag = -1
      left = 1
      RETURN
 
 100  mflag = 0
      left = ilo
      RETURN
 
 110  mflag = 1
      left = lxt
      RETURN
      end
