      subroutine yint (v,w,y,z,m,h)
************************************************************************
*  this program calculates the indefinite integrals y and z using the
*  lagrange integration formula
*
*  y(r) = integral of v from 0 to r
*  z(r) = integral of w from r to infinity
*  m is the maximum tabulation point of v and w (virtual infinity)
*  h is the step size of the radial grid
************************************************************************
      implicit doubleprecision(a-h,o-z)
      include "global.par"
	parameter(NO=8,NP=NO/2)
	common/yindat/a(NO,NP),b(NP)
      dimension v(NGP),w(NGP),y(NGP),z(NGP)
      dimension dy(NGP),dz(NGP)
      
      y(1) = 0.0d0
      z(m) = 0.0d0
      do 200 i = 2,NP
         k = m - i + 1
         y(i) = y(i-1)
         z(k) = z(k+1)
         ii = i - 1
         do 190 j = 1,NO
            y(i) = y(i) + a(j,ii)*v(j)
            z(k) = z(k) + a(j,ii)*w(m-j+1)
 190     continue
 200  continue
      im = NP + 1
      in = m - NP + 1
******c$omp  parallel do private(i,k) shared(b,dy,dz,v,w,im,in,m)
      do 240 i = im,in
         k = m - i + 1
            dy(i) = b(1)*(v(i+1-1) + v(i-1))
     2            + b(2)*(v(i+1) + v(i-2))
     3            + b(3)*(v(i+2) + v(i-3))
     4            + b(4)*(v(i+3) + v(i-4))
c     5            + b(5)*(v(i+4) + v(i-5))
            dz(k) = b(1)*(w(k) + w(k+1))
     2            + b(2)*(w(k-1) + w(k+2))
     3            + b(3)*(w(k-2) + w(k+3))
     4            + b(4)*(w(k-3) + w(k+4))
c     5            + b(5)*(w(k-4) + w(k+4))
 240     continue        

*240  continue

*     do 240 i = im,in
*        k = m - i + 1
*        dy(i) = 0.0d0
*        dz(k) = 0.0d0
*        do 230 j = 1,NP
*           dy(i) = dy(i) + b(j)*(v(i+j-1) + v(i-j))
*           dz(k) = dz(k) + b(j)*(w(k-j+1) + w(k+j))
*230     continue        
*240  continue


*****c$omp  parallel do private(i,k) shared(dy,dz,y,z,im,in,m)
      do 250 i = im,in
         k = m - i + 1
         y(i) = y(i-1) + dy(i)
         z(k) = z(k+1) + dz(k)
 250  continue

      in = in + 1
      do 300 i = in,m
         k = m - i + 1
         y(i) = y(i-1)
         z(k) = z(k+1)
         do 280 j = 1,NO
            y(i) = y(i) + a(j,k)*v(m-j+1)
            z(k) = z(k) + a(j,k)*w(j)
 280     continue
 300  continue
      return
      end
