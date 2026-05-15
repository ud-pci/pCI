 
      subroutine bsplvb(t,jhigh,index,x,left,biatx)
      implicit doubleprecision(a-h,o-z)
      parameter(JMAX=20,NX=50,KX=15)
      parameter (NK=NX+KX)
      dimension biatx(jhigh),t(nk),deltal(jmax),deltar(jmax)
      data j/1/
 
      go to (10,20),index
 
 10   j = 1
      biatx(1) = 1.0
      if(j .ge. jhigh)                 go to 99
 
 20   continue
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.0
         do 26 i = 1,j
             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
 26      continue
         biatx(jp1) = saved
         j = jp1
         if(j .lt. jhigh)
     *go to 20
 
 99   RETURN
      end

