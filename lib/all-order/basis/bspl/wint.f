 
      subroutine wint (u,v,s,t,c,d,e,f,m,h)
************************************************************************
*
*  this program calculates the indefinite integrals u,v,s,t using the
*  lagrange integration formula
*
*  c(r) = integral of u from 0 to r
*  d(r) = integral of v from 0 to r
*  e(r) = integral of s from r to infinity
*  f(r) = integral of t from r to infinity
*  m is the maximum tabulation point  (virtual infinity)
*  h is the step size of the radial grid
*
************************************************************************
      implicit doubleprecision(a-h,o-z)
      logical first 
      parameter(NGP=500)
      dimension u(NGP),v(NGP),s(NGP),t(NGP)
      dimension c(NGP),d(NGP),e(NGP),f(NGP)
      dimension dc(NGP),dd(NGP),de(NGP),df(NGP)
*
*     lagrange 6 point integration formula
*     ************************************
*      parameter(NO=6,NP=NO/2)
*      dimension aa(NP,NO),a(NO,NP),b(NP)
*      data  da/1440.d0/
*      data  aa/  475.d0,   -27.d0,    11.d0,
*     2          1427.d0,   637.d0,   -93.d0,
*     3          -798.d0,  1022.d0,   802.d0,
*     4           482.d0,  -258.d0,   802.d0,
*     5          -173.d0,    77.d0,   -93.d0,
*     6            27.d0,   -11.d0,    11.d0/
*      data   b/  802.d0,   -93.d0,    11.d0/
************************************************************************
*     lagrange 8 point integration formula
*     ************************************
*      parameter(NO=8,NP=NO/2)
*      dimension aa(NP,NO),a(NO,NP),b(NP)
*      data  da/120960.d0/
*      data  aa/  36799.d0,  -1375.d0,   351.d0,  -191.d0,
*     2          139849.d0,  47799.d0, -4183.d0,  1879.d0,
*     3         -121797.d0, 101349.d0, 57627.d0, -9531.d0,
*     4          123133.d0, -44797.d0, 81693.d0, 68323.d0,
*     5          -88547.d0,  26883.d0,-20227.d0, 68323.d0,
*     6           41499.d0, -11547.d0,  7227.d0, -9531.d0,
*     7          -11351.d0,   2999.d0, -1719.d0,  1879.d0,
*     8            1375.d0,   -351.d0,   191.d0,  -191.d0/   
*      data   b/  68323.d0,  -9531.d0,  1879.d0,  -191.d0/
************************************************************************
*     lagrange 10 point integration formula
*     ************************************
      parameter(NO=10,NP=NO/2)
      dimension aa(NP,NO),a(NO,NP),b(NP)
      data da/7257600.d0/
      data aa/2082753.d0,  -57281.d0,   10625.d0,   -3969.d0,   2497.d0, 
     2        9449717.d0, 2655563.d0, -163531.d0,   50315.d0, -28939.d0,
     3      -11271304.d0, 6872072.d0, 3133688.d0, -342136.d0, 162680.d0,
     4       16002320.d0,-4397584.d0, 5597072.d0, 3609968.d0,-641776.d0,
     5      -17283646.d0, 3973310.d0,-2166334.d0, 4763582.d0,4134338.d0,
     6       13510082.d0,-2848834.d0, 1295810.d0,-1166146.d0,4134338.d0, 
     7       -7394032.d0, 1481072.d0, -617584.d0,  462320.d0,-641776.d0, 
     8        2687864.d0, -520312.d0,  206072.d0, -141304.d0, 162680.d0,
     9        -583435.d0,  110219.d0,  -42187.d0,   27467.d0, -28939.d0,
     a          57281.d0,  -10625.d0,    3969.d0,   -2497.d0,   2497.d0/
      data b/ 4134338.d0, -641776.d0,  162680.d0,  -28939.d0,   2497.d0/
*********1*********2*********3*********4*********5*********6*********7**
      data first /.true./ 
*
*  note that a different even order method can be used by replacing the
*  dimension and data statements in above block 
************************************************************************
*   initialize the intigration coefficients on first entry
      if(first) then 
         hd = h/da
         do 150 i = 1,NP
            do 100 j = 1,NO
               a(j,i) = aa(i,j)*hd
 100        continue
            b(i) = b(i)*hd   
 150     continue
         first = .false.
      endif
      c(1) = 0.0
      d(1) = 0.0
      e(m) = 0.0
      f(m) = 0.0
      do 180 i = 2,NP
         k = m-i+1
         c(i) = c(i-1)
         d(i) = d(i-1)
         e(k) = e(k+1)
         f(k) = f(k+1)
         ii = i-1
         do 170 j = 1,NO
              c(i) = c(i) + a(j,ii)*u(j)
              d(i) = d(i) + a(j,ii)*v(j)
              e(k) = e(k) + a(j,ii)*s(m-j+1)
              f(k) = f(k) + a(j,ii)*t(m-j+1)
 170     continue
 180  continue
      im = NP + 1
      in = m - NP + 1
      do 200 i = im,in
         dc(i) =  b(1)*(u(i  )+u(i-1))+b(2)*(u(i+1)+u(i-2))
     &           +b(3)*(u(i+2)+u(i-3))+b(4)*(u(i+3)+u(i-4))
     &           +b(5)*(u(i+4)+u(i-5))
         dd(i) =  b(1)*(v(i  )+v(i-1))+b(2)*(v(i+1)+v(i-2))
     &           +b(3)*(v(i+2)+v(i-3))+b(4)*(v(i+3)+v(i-4))
     &           +b(5)*(v(i+4)+v(i-5))
         de(i) =  b(1)*(s(i  )+s(i-1))+b(2)*(s(i+1)+s(i-2))
     &           +b(3)*(s(i+2)+s(i-3))+b(4)*(s(i+3)+s(i-4))
     &           +b(5)*(s(i+4)+s(i-5))
         df(i) =  b(1)*(t(i  )+t(i-1))+b(2)*(t(i+1)+t(i-2))
     &           +b(3)*(t(i+2)+t(i-3))+b(4)*(t(i+3)+t(i-4))
     &           +b(5)*(t(i+4)+t(i-5))
 200  continue
      do 250 i = im,in
         c(i) = c(i-1) + dc(i)
         d(i) = d(i-1) + dd(i)
 250  continue
      km = m - NP
      kn = NP
      do 260 k = km,kn,-1
         e(k) = e(k+1) + de(k+1)
         f(k) = f(k+1) + df(k+1)
 260  continue
      in = in + 1
      do 350 i = in,m
          k = m - i + 1
          c(i) = c(i-1)
          d(i) = d(i-1)
          e(k) = e(k+1)
          f(k) = f(k+1)
          do 300 j = 1,NO
              c(i) = c(i) + a(j,k)*u(m-j+1)
              d(i) = d(i) + a(j,k)*v(m-j+1)
              e(k) = e(k) + a(j,k)*s(j)
              f(k) = f(k) + a(j,k)*t(j)
 300      continue
 350  continue
      return
      end
