
      function rint (f,na,nb,nq,h)
      implicit double precision(a-h,o-z)
c
c  this program calculates the integral of the function f from point na
c  to point nb using a nq points quadrature ( nq is any integer between
c  1 and 14 ).  h is the grid size.
c                                      written by c. c. j. roothaan
c
      dimension c(105),c1(78),c2(27),d(14),f(nb)
      equivalence (c1(1),c(1)),(c2(1),c(79))
      data c1/
     a 1.d0,
     b 2.d0,1.d0,
     c 23.d0,28.d0,9.d0,
     d 25.d0,20.d0,31.d0,8.d0,
     e 1413.d0,1586.d0,1104.d0,1902.d0,475.d0,
     f 1456.d0,1333.d0,1746.d0,944.d0,1982.d0,459.d0,
     g 119585.d0,130936.d0,89437.d0,177984.d0,54851.d0,176648.d0,
     g  36799.d0,
     h 122175.d0,111080.d0,156451.d0,46912.d0,220509.d0,29336.d0,
     h 185153.d0, 35584.d0,
     i 7200319.d0, 7783754.d0,5095890.d0,12489922.d0,-1020160.d0,
     i16263486.d0,  261166.d0,11532470.d0,2082753.d0,
     j 7305728.d0,  6767167.d0, 9516362.d0, 1053138.d0,18554050.d0,
     j-7084288.d0, 20306238.d0,-1471442.d0,11965622.d0, 2034625.d0,
     k  952327935.d0, 1021256716.d0,  636547389.d0,1942518504.d0,
     k-1065220914.d0, 3897945600.d0,-2145575886.d0,3373884696.d0,
     k -454944189.d0, 1637546484.d0,  262747265.d0,
     l  963053825.d0,  896771060.d0, 1299041091.d0, -196805736.d0,
     l 3609224754.d0,-3398609664.d0, 6231334350.d0,-3812282136.d0,
     l 4207237821.d0, -732728564.d0, 1693103359.d0,  257696640.d0/
      data c2 / 5206230892907.d0,5551687979302.d0,3283609164916.d0,
     m 12465244770050.d0,-13155015007785.d0,39022895874876.d0,
     m-41078125154304.d0,53315213499588.d0,-32865015189975.d0,
     m 28323664941310.d0,-5605325192308.d0,  9535909891802.d0,
     m  1382741929621.d0,
     n  5252701747968.d0,  4920175305323.d0,  7268021504806.d0,
     n -3009613761932.d0, 28198302087170.d0,-41474518178601.d0,
     n 76782233435964.d0,-78837462715392.d0, 81634716670404.d0,
     n-48598072507095.d0, 34616887868158.d0, -7321658717812.d0,
     n  9821965479386.d0,  1360737653653.d0/
      data d/2.d0,2.d0,24.d0,24.d0,1440.d0,1440.d0,120960.d0,
     a  120960.d0,7257600.d0,7257600.d0,958003200.d0,958003200.d0,
     b  5230697472000.d0,5230697472000.d0/

      a = 0.0d0
      l = na
      m = nb
      i = nq*(nq+1)/2
      do 100 j = 1,nq
         a = a + c(i)*( f(l) + f(m) )
         l = l + 1
         m = m - 1
         i = i - 1
 100  continue
      a = a/d(nq)
      do 200 n = l,m
        a = a + f(n)
 200  continue
      rint = a*h
      return
      end
