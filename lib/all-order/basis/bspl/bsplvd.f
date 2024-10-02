      subroutine bsplvd(t,k,x,left,dbiatx,nderiv)
      implicit doubleprecision(a-h,o-z)
*******  changes from de Boor
      parameter(NX=40,KX=15)
      parameter (NK=NX+KX)
      dimension a(KX,KX),dbiatx(KX,KX)
*  also a missing in arg list
*****************
      dimension t(NK)
 
      mhigh = max0(min0(nderiv,k),1)
      kp1 = k + 1
      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
      if(mhigh .eq. 1) go to 99
 
      ideriv = mhigh
      do 15 m = 2,mhigh
          jp1mid = 1
          do 11 j = ideriv,k
             dbiatx(j,ideriv) = dbiatx(jp1mid,1)
             jp1mid = jp1mid + 1
 11       continue
          ideriv = ideriv - 1
          call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
 15   continue
 
      jlow = 1
      do 20 i = 1,k
         do 19 j = jlow,k
            a(j,i) = 0.0
 19      continue
         jlow = i
         a(i,i) = 1.0
 20   continue
 
      do 40 m = 2,mhigh
         kp1mm = kp1 - m
         fkp1mm = kp1mm
         il = left
         i = k
 
         do 25 ldummy = 1,kp1mm
            factor = fkp1mm/(t(il+kp1mm) - t(il))
            do 24 j = 1,i
               a(i,j) = (a(i,j) -a(i-1,j))*factor
 24         continue
            il = il - 1
            i = i - 1
 25      continue
 
         do 36 i = 1,k
              sum = 0.0
              jlow = max0(i,m)
              do 35 j = jlow,k
                 sum = sum + a(j,i)*dbiatx(j,m)
 35           continue
              dbiatx(i,m) = sum
 36      continue
 40   continue
 
 99   RETURN
      end

