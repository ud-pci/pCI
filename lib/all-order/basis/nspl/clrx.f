      
      function clrx(kap1,k,kap2)
      implicit doubleprecision(a-h,o-z)
c
c   the value of clrx is the 3-j symbol
c
c            j1      k        j2
c           0.5      0      -0.5
c
c   kap1=kappa value for j1
c   kap2=kappa value for j2
c
c   triangular conditions are tested in code
c
c   see : angular momentum (second edition) by brink and satchler
c              page 138
c
c  no subroutines called.
c
      kma=iabs(kap1)
      kmb=iabs(kap2)
      jp=kma+kmb-1
      jm=kma-kmb
      x=1d0/dble(kma*kmb)
      ix=jp-k
      j=1
      go to 10
    1 ix=jp+k+1
      j=1
      go to 14
    2 ix=jm+k
      j=2
      go to 10
    3 ix=k-jm
      j=3
      go to 10
    4 continue
      y=dsqrt(x)
      ip=k
      jq=jp+k
      if (mod(jq,2).eq.0) go to 5
      y=-y
      ip=k+1
    5 x=1d0
      ix=(ip+jp)/2
      j=4
      go to 10
    6 ix=(jp-ip)/2
      j=2
      go to 14
    7 ix=(jm+ip-1)/2
      j=3
      go to 14
    8 ix=(ip-1-jm)/2
      j=4
      go to 14
    9 clrx=y*x
      iphase=(jp+ip-2)/2
      if (mod(iphase,2).eq.1) clrx=-clrx
      return
c
   10 if (ix) 18,13,11
   11 do 12 i=1,ix
      x=x*dble(i)
   12 continue
   13 go to (1,3,4,6),j
c
   14 if (ix) 18,17,15
   15 do 16 i=1,ix
      x=x/dble(i)
   16 continue
   17 go to (2,7,8,9),j
c
   18 clrx=0d0
      return
      end
