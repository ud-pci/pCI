      subroutine yfun(x,y,l,m,*)
*     this is a new version
************************************************************************
*
*  this program calculates the Hartree's y-functions
*
*  x  : input function
*  y  : output hartree's y-function : y(l,r)/r =
*            r**(-l-1) * Int {0,r} [ r**l x(r) dr ] 
*               + r**l * Int{r,inf}[ r**(-l-1) x(r) ] 
*  l  : order of the y function (must be .ge. 0)
*  m  : number of tabulation points for the input function x
*
************************************************************************
C	some added changes by APD 3/18/98
C     since it's possible that m < NGP
C     so the first call is assumed to be with m = max
C     where max is the last point on the grid 
      implicit doubleprecision(a-h,o-z)
      logical first
      parameter(NGP=500,LX=15)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      dimension x(NGP),y(NGP)
      dimension u(NGP),v(NGP),w(NGP)
      dimension xl(NGP,-LX-1:LX)
      data first/.TRUE./
*
      if(first) then
         first = .FALSE.
         do 400 i = 1,NGP
            xl(i,0) = 1.0d0
 400     continue
         do 460 ll = 1,LX
            do 450 i = 1,NGP
               xl(i,ll) = xl(i,ll-1)*r(i)
 450        continue
 460     continue
         do 480 ll = 1,LX
            xl(1,-ll) = 0.0d0
C	 Below it was  do 470 i=2,NGP
            do 470 i=2,max
               xl(i,-ll) = 1d0/xl(i,ll)
 470        continue
 480     continue 
         xl(1,-LX-1) = 0.0d0
C	 Below it was  do 490 i=2,NGP
         do 490 i=2,max
            xl(i,-LX-1) = xl(i,-LX)/r(i)
 490     continue
      endif 
      v(1)=0.0
      w(1)=0.0
      if(l.lt.0) then
         write(6,*) ' Error in yint: l < 0 '
         return 1
      elseif(l.eq.0) then 
         do 100  i= 2,m
            v(i) = x(i)*rp(i)
            w(i) = x(i)*rpor(i)
 100     continue
         call yint(v,w,y,u,m,h)
         ym = y(m)
         y(1) = u(1)
         do 120 i = 2,m
            y(i) = y(i)*xl(i,-1) + u(i)
 120     continue 
         m1 = m + 1
         if(m1.le.max) then
           do 140 i = m1,max
              y(i) = ym*xl(i,-1)
 140       continue
         endif
         return
      else
         if(l.gt.LX) then
            write(6,1000) l
 1000       format('  Error in yint:  the value of l = ',
     &      i4,' is too large')
            return 1
         endif
         do 180 i = 2,m
            v(i) = x(i)*rp(i)*xl(i,l)
            w(i) = x(i)*rp(i)*xl(i,-l-1)
 180     continue
         call yint(v,w,y,u,m,h)
         ym = y(m)
         do 200 i = 2,m
            y(i) = y(i)*xl(i,-l-1) + u(i)*xl(i,l)
 200     continue
         m1 = m + 1
         if(m1.le.max) then
            do 210  i = m1,max
               y(i) = ym*xl(i,-l-1)
 210        continue
         endif
      endif
      return
      end

