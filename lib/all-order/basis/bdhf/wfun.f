
      subroutine wfun(a,b,x,y,l,m,*)
************************************************************************
*
*  this program calculates the Retardation Integrals
*
*  x  : input function
*
*        -l        (l-1)             -(l+2)       (l+1)
*  x  : r   Int[ r'     a(r') dr'] - r     Int[ r'      a(r') dr']
*
*        (l-1)       -l               l+1        -(l+2)
*  y  : r     Int[ r'   b(r') dr'] - r     Int[ r'      b(r') dr']
*
*
*  l  : order of the y function (must be .ge. 0)
*  m  : number of tabulation points for the input function x
*
************************************************************************
      implicit doubleprecision(a-h,o-z)
      parameter(NGP=500)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      dimension x(NGP),y(NGP),a(NGP),b(NGP)
      dimension s(NGP),t(NGP),u(NGP),v(NGP),w(NGP)
      dimension c(NGP),d(NGP),e(NGP),f(NGP)
      data ja/16/
*
*  note that ja is determined by the order of the integration method
*  used in subroutine wint
*
      u(1) = 0.0
      v(1) = 0.0
      w(1) = 0.0
      if(l.le.0)  then
         write(6,*) ' Error in wfun: l < 0 '
         return 1
      else
         if(l.gt.ja) then
            write(6,1000) l
 1000       format('  Warning in wfun:  the value of l =',
     &      i4,' may be too large')
         endif
         do 100 i = 2,m
            w(i) = r(i)**l
            u(i) = a(i)*rpor(i)*w(i)
            v(i) = u(i)*r(i)*r(i)
            s(i) = b(i)*rp(i)/w(i)
            t(i) = s(i)/(r(i)*r(i))
 100     continue
         call wint(u,v,s,t,c,d,e,f,m,h)
         cm = c(m)
         dm = d(m)
         x(1) = 0.0
         y(1) = 0.0
         do 120 i = 2,m
            x(i) = ( c(i) - d(i)/(r(i)*r(i)) )/w(i)
            y(i) = w(i)*( e(i)/r(i) - r(i)*f(i) )
 120     continue
         m1 = m + 1
         if(m1.le.NGP) then
            do 130 i = m1,NGP
               x(i) = ( cm - dm/r(i)**2 )/r(i)**l
               y(i) = 0.0
 130        continue
         endif
      endif
      return
      end
