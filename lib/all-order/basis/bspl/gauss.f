
      subroutine gauss(n,x,w)
      implicit doubleprecision(a-h,o-z)
******************************************************************
*
*  Gaussian coordinates and weights for the interval [0..1]
*        adapted from "setgau" in Numerical Recipes
*******************************************************************
      dimension x(n),w(n)
      data eps /1d-15/
      pi = dacos(-1.0d0)
      m = (n+1)/2
      do 120 i = 1,m
         z = cos(pi*(i-0.25d0)/(n+0.5d0))
 100     continue
            p1 = 1.d0
            p2 = 0.d0
            do 110 j = 1,n
               p3 = p2
               p2 = p1
               p1 = (( 2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 110        continue
            pp = n*(z*p1-p2)/(z*z-1.d0)
            z1 = z
            z = z1-p1/pp
         if(dabs(z-z1).gt.eps) go to 100
         x(i) = 0.5d0*(1.0d0 - z)
         x(n+1-i) = 0.5d0*(1.0d0 + z)
         w(i) = 1.d0/((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
 120  continue
      return
      end
