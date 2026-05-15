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
      logical first 
      parameter(NGP=500)
      dimension v(NGP),w(NGP),y(NGP),z(NGP)
      dimension dy(NGP),dz(NGP)
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

       do 240 i = im,in
          k = m - i + 1
          dy(i) = 0.0d0
          dz(k) = 0.0d0
          do 230 j = 1,NP
             dy(i) = dy(i) + b(j)*(v(i+j-1) + v(i-j))
             dz(k) = dz(k) + b(j)*(w(k-j+1) + w(k+j))
  230     continue        
  240  continue

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
