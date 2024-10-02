      subroutine zfun(x,y,l,m,*)
************************************************************************
*
*  this program calculates the Magnetic Integrals
*
*  x  : input function
*
*        -l        (l-1)             -(l+2)        (l+1)
*  y  : r   Int[ r'     x(r') dr'] - r       Int[ r'     x(r') dr']
*
*
*  l  : order of the y function (must be .ge. 0)
*  m  : number of tabulation points for the input function x
*
************************************************************************
      implicit doubleprecision(a-h,o-z)
      parameter(NGP=500,LX=15)
      logical first
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      dimension x(NGP),y(NGP)
      dimension u(NGP),v(NGP),w(NGP)
      dimension xl(NGP,-LX-2:LX+1)
      data first/.TRUE./
*
      if(first) then
         first = .FALSE.
         do 400 i = 1,NGP
            xl(i,0) = 1.0d0
 400     continue
         do 460 ll = 1,LX+1
            do 450 i = 1,NGP
               xl(i,ll) = xl(i,ll-1)*r(i)
 450        continue
 460     continue
         do 480 ll = 1,LX+1
            xl(1,-ll) = 0.0d0
            do 470 i=2,NGP
               xl(i,-ll) = 1d0/xl(i,ll)
 470        continue
 480     continue 
         xl(1,-LX-2) = 0.0d0
         do 490 i=2,NGP
            xl(i,-LX-2) = xl(i,-LX-1)/r(i)
 490     continue
      endif 
*
      v(1) = 0.0
      w(1) = 0.0
      if(l.le.0) then
         write(6,*) ' Error in zfun: l <= 0 '
         return 1
      else
         if(l.gt.LX) then
            write(6,1000) l
 1000       format('  Warning in zfun:  the value of l =',
     &      i4,' may be too large')
         endif
         do 100 i = 2,m
            v(i) = x(i)*rp(i)*xl(i,l-1)
            w(i) = x(i)*rp(i)*xl(i,l+1)
 100     continue
         call zint(v,w,y,u,m,h)
         ym = y(m)
         um = u(m)
         do 120 i = 2,m
            y(i) = y(i)*xl(i,-l) - u(i)*xl(i,-l-2)
 120     continue
         m1 = m + 1
         if(m1.le.NGP) then
            do 130 i = m1,NGP
               y(i) = ym*xl(i,-l) - um*xl(i,-l-2)
 130        continue
         endif 
      endif
      return
      end
