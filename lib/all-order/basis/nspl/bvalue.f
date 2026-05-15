      function bvalue(t,bcoef,n,k,x,jderiv)
      implicit doubleprecision(a-h,o-z)
      parameter(KMAX=20,NX=50,KX=15)
      parameter (NK=NX+KX)
      dimension bcoef(n),t(NK),aj(KMAX),dl(KMAX),dr(KMAX)
 
      bvalue = 0.0
      if(jderiv .ge. k)                go to 99
      call interv (t,n+k,x,i,mflag)
      if(mflag .ne. 0)                 go to 99
 
      km1 = k - 1
      if(km1 . gt. 0)                  go to 1
      bvalue = bcoef(i)
                                       go to 99
 
 1    jcmin = 1
      imk = i - k
      if(imk .ge. 0)                    go to 8
      jcmin = 1 - imk
 
      do 5 j = 1,i
         dl(j) = x - t(i+1-j)
 5    continue
      do 6 j = i,km1
         aj(k-j)  = 0.0
         dl(j) = dl(i)
 6    continue
                                      go to 10
 
 8    do 9 j = 1,km1
         dl(j) = x - t(i+1-j)
 9    continue
 
 10   jcmax = k
      nmi = n - i
      if(nmi .ge. 0)                  go to 18
      jcmax = k + nmi
      do 15 j = 1,jcmax
          dr(j) = t(i+j) - x
 15   continue
      do 16 j = jcmax,km1
          aj(j+1) = 0.0
          dr(j) = dr(jcmax)
 16   continue
                                     go to 20
 
 18   do 19 j = 1,km1
        dr(j) = t(i+j) - x
 19   continue
 
 20   do 21 jc = jcmin,jcmax
          aj(jc) = bcoef(imk + jc)
 21   continue
 
      if(jderiv .eq. 0)             go to 30
      do 23 j = 1,jderiv
         kmj = k - j
         fkmj = kmj
         ilo = kmj
         do 22 jj = 1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
 22      continue
 23   continue
 
 30   if(jderiv .eq. km1)           go to 39
      do 33 j = jderiv+1,km1
         kmj = k - j
         ilo = kmj
         do 32 jj = 1,kmj
           aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
           ilo = ilo - 1
 32      continue
 33   continue
 39   bvalue = aj(1)
 
 99   RETURN
      end

