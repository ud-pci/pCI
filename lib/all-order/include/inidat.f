      subroutine inidat
      implicit DOubleprecision(a-h,o-z)
      include "global.par"
      parameter(NMAX=100,LX=15)
	parameter(NO=8,NP=NO/2)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/clebsc/p(NMAX)
      common/yfndat/xl(NGP,-LX-1:LX)      
      common/yindat/a(NO,NP),b(NP)
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
c      parameter(NO=8,NP=NO/2)
 
      dimension aa(NP,NO),bb(NP)
      data  da/120960.d0/
      data  aa/  36799.d0,  -1375.d0,   351.d0,  -191.d0,
     2          139849.d0,  47799.d0, -4183.d0,  1879.d0,
     3         -121797.d0, 101349.d0, 57627.d0, -9531.d0,
     4          123133.d0, -44797.d0, 81693.d0, 68323.d0,
     5          -88547.d0,  26883.d0,-20227.d0, 68323.d0,
     6           41499.d0, -11547.d0,  7227.d0, -9531.d0,
     7          -11351.d0,   2999.d0, -1719.d0,  1879.d0,
     8            1375.d0,   -351.d0,   191.d0,  -191.d0/   
      data   bb/  68323.d0,  -9531.d0,  1879.d0,  -191.d0/
************************************************************************
*     lagrange 10 point integration formula
*     ************************************

c

c      data da/7257600.d0/
c      data aa/2082753.d0,  -57281.d0,   10625.d0,   -3969.d0,   2497.d0, 
c     2        9449717.d0, 2655563.d0, -163531.d0,   50315.d0, -28939.d0,
c     3      -11271304.d0, 6872072.d0, 3133688.d0, -342136.d0, 162680.d0,
c     4       16002320.d0,-4397584.d0, 5597072.d0, 3609968.d0,-641776.d0,
c     5      -17283646.d0, 3973310.d0,-2166334.d0, 4763582.d0,4134338.d0,
c     6       13510082.d0,-2848834.d0, 1295810.d0,-1166146.d0,4134338.d0, 
c     7       -7394032.d0, 1481072.d0, -617584.d0,  462320.d0,-641776.d0, 
c     8        2687864.d0, -520312.d0,  206072.d0, -141304.d0, 162680.d0,
c     9        -583435.d0,  110219.d0,  -42187.d0,   27467.d0, -28939.d0,
c     a          57281.d0,  -10625.d0,    3969.d0,   -2497.d0,   2497.d0/
c      data bb/4134338.d0, -641776.d0,  162680.d0,  -28939.d0,   2497.d0/
 
****    initialize the factorial array for d6j

      p(1)=1.d0
      p(2)=1.d0

      DO i=3,NMAX
         p(i)=p(i-1)*(i-1)
      END DO

****   initialize the integration coefficients for yint
 
      hd = h/da
      DO i = 1,NP
         DO j = 1,NO
            a(j,i) = aa(i,j) * hd
         END DO
         b(i) = bb(i) * hd   
      END DO

****   initialize r^l and r^(-l-1) for yfun

      DO i = 1,max
         xl(i,0) = 1.0d0
      END DO
      DO ll = 1,LX
         DO  i = 1,max
               xl(i,ll) = xl(i,ll-1)*r(i)
         END DO
      END DO
      DO ll = 1,LX
         xl(1,-ll) = 0.0d0
         DO i=2,max
               xl(i,-ll) = 1d0/xl(i,ll)
         END DO
      END DO 
      xl(1,-LX-1) = 0.0d0
      DO i=2,max
            xl(i,-LX-1) = xl(i,-LX)/r(i)
      END DO
       
      RETURN
	END
