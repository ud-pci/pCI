
      subroutine st(km,nn)
      implicit real*8 (a-h,o-z)
      include "global.par"
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore       
      nn1=0
      do 1 i=1,ncore
        if (km.eq.ko(i)) then
         if (nn1.lt.no(i)) then
           nn1=no(i)
         endif
        endif
1     continue
      CALL klj(km,kapm,lm,jm,indm,n0m)
      if (nn1.eq.0) then
       nn=1
      else
       nn=nn1+1-n0m
      endif
      end
     


 
      SUBROUTINE mrr (n,k,i)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      IF (n.eq.1.and.k.eq.-1) i=1

      IF (n.eq.2.and.k.eq.-1) i=2
      IF (n.eq.2.and.k.eq. 1) i=3
      IF (n.eq.2.and.k.eq.-2) i=4

      IF (n.eq.3.and.k.eq.-1) i=5
      IF (n.eq.3.and.k.eq. 1) i=6
      IF (n.eq.3.and.k.eq.-2) i=7
      IF (n.eq.3.and.k.eq. 2) i=8
      IF (n.eq.3.and.k.eq.-3) i=9

      IF (n.eq.4.and.k.eq.-1) i=10
      IF (n.eq.4.and.k.eq. 1) i=11
      IF (n.eq.4.and.k.eq.-2) i=12
      IF (n.eq.4.and.k.eq. 2) i=13
      IF (n.eq.4.and.k.eq.-3) i=14
      IF (n.eq.4.and.k.eq. 3) i=15
      IF (n.eq.4.and.k.eq.-4) i=16

      IF (n.eq.5.and.k.eq.-1) i=17
      IF (n.eq.5.and.k.eq. 1) i=18
      IF (n.eq.5.and.k.eq.-2) i=19
      IF (n.eq.5.and.k.eq. 2) i=20
      IF (n.eq.5.and.k.eq.-3) i=21
      IF (n.eq.5.and.k.eq. 3) i=22
      IF (n.eq.5.and.k.eq.-4) i=23
      IF (n.eq.5.and.k.eq. 4) i=24
      IF (n.eq.5.and.k.eq.-5) i=25

      IF (n.eq.6.and.k.eq.-1) i=26
      IF (n.eq.6.and.k.eq. 1) i=27
      IF (n.eq.6.and.k.eq.-2) i=28
      IF (n.eq.6.and.k.eq. 2) i=29
      IF (n.eq.6.and.k.eq.-3) i=30
      IF (n.eq.6.and.k.eq. 3) i=31
      IF (n.eq.6.and.k.eq.-4) i=32
      IF (n.eq.6.and.k.eq. 4) i=33
      IF (n.eq.6.and.k.eq.-5) i=34
      IF (n.eq.6.and.k.eq. 5) i=35
      IF (n.eq.6.and.k.eq.-6) i=36

      RETURN
      END


        
      subroutine trg(g1,g2,g3,i)
      implicit real*8 (a-h,o-z)
      i=0
      do 1 g=g1,g2
        if (g.eq.g3) then
        i=1
        goto 2
        endif
1     continue 
c       write (*,*) 'Warning in trg!'
2     continue
      return
      end
c      
      function delta(a,b,c)
      implicit real*8 (a-h,o-z)    
      delta=DSQRT((f1(a+b-c)*f1(a-b+c)*f1(-a+b+c))/f1(a+b+c+1))
      if(((a+b-c).lt.0).or.((a-b+c).lt.0).or.((-a+b+c).lt.0)) delta=0
      return
      end
c      
      function f1(n)      
      real*8 f1,n
       f1=1.d0        
       if (n.eq.0) goto 6
       do  i=1,n
        f1=f1*i
       enddo 
6     return
      end 
      
      
      function s(k,k1,k2) 
      implicit doubleprecision(a-h,o-z)
      call klj(k1,kap1,l1,j1,ind,n)
      call klj(k2,kap2,l2,j2,ind,n)
      a=l1+l2+k 
      b=a/2-INT(a/2)
      if (b.ne.0.0) then
      c1=0.d0
      else
      c1=1.d0
      endif 
      coef=c3j(j1,j2,2*k,-1,1,0)
      s=2.d0*c1*coef*SQRT(kap1*kap2*1.d0)*((-1)**kap1)
      end    
c     
      function c3j(j1,j2,j3,m1,m2,m3)
      implicit real*8 (a-h,o-z)
      g1=0.5d0*j1
      g2=0.5d0*j2
      g3=0.5d0*j3
      p1=0.5d0*m1
      p2=0.5d0*m2
      p3=0.5d0*m3
      if (p1+p2+p3.ne.0.0) goto 100
      call trg(abs(g1-g2),g1+g2,g3,i1)
      call trg (-g1,g1,p1,i2)
      call trg (-g2,g2,p2,i3)
      call trg (-g3,g3,p3,i4)
      i=i1*i2*i3*i4
      if (i.eq.0) then
       c3j=0.d0
       goto 100 
      endif
      a1=f2(g1-p1)*f2(g1+p1)*f2(g2-p2)*f2(g2+p2)*f2(g3-p3)*f2(g3+p3)   
      a2=f2(g1+g2-g3)*f2(g1-g2+g3)*f2(-g1+g2+g3)
      a3=f2(g1+g2+g3+1)
      a=SQRT((a1*a2)/a3)
      imax=min(g2+p2,g1-p1,g1+g2-g3)
      imin=max(0.d0,-g3+g2-p1,-g3+g1+p2)
      sum=0.d0
      do 4 i=imin,imax
        a4=(-1)**(i+g1-g2-p3)
        a5=f2(i+0.d0)*f2(g1+g2-g3-i)*f2(g1-p1-i)*f2(g2+p2-i)
        a6=f2(g3-g2+p1+i)*f2  (g3-g1-p2+i)
        if ((a5*a6).eq.0) then
         c3j=0.d0
         goto 100
        endif 
        sum=sum+a4/(a5*a6)
4     continue
      c3j=a*sum
100   continue
      end
c
      function f2(n)
      real*8 f2,n
       f2=0.d0
       if (n.lt.0) goto 6
       f2=1.d0        
       if (n.eq.0) goto 6
       do  i=1,n
        f2=f2*i
       enddo 
6     return
      end       
c
      subroutine klj(k,kap,l,j,ind,n0)
      kap=iabs(k)
      j=2*kap-1 
      if (k.lt.0) then
      ind=-2*k-1
      n0=kap-1
      l=-k-1
      else
      n0=k
      l=k
      ind=2*k
      endif
      end    
         

      subroutine indk (ind,k,kap,l,j)
      call odd (ind,i)
      if (i.eq.0.d0) then
       k=-(ind+1)/2
       else
       k=ind/2
      endif
      kap=iabs(k)
      if (k.lt.0) then
      l=-k-1
      else
      l=k
      endif
      j=2*kap-1
      end     

      subroutine indk1 (ind,k,kap,l,j,n0)
      call odd (ind,i)
      if (i.eq.0.d0) then
       k=-(ind+1)/2
       else
       k=ind/2
      endif
      kap=iabs(k)
      if (k.lt.0) then
      l=-k-1
      n0=kap-1
      else
      l=k
      n0=k
      endif
      j=2*kap-1
      end     

c
      subroutine odd (l,i)          
      b=l/2.0-INT(l/2.0)
      if (b.ne.0.0) then
      i=0
      else
      i=1
      endif
      return
      end 



      subroutine trgi(j1,j2,j,i)
      i=0
      do 1 ii=j1,j2
        if (j.eq.ii) then
        i=1
        goto 2
        endif
1     continue
2     continue
      return
      end

      subroutine odd2 (l,i)          
      b=l/4.0-INT(l/4.0)
      if (b.ne.0.0) then
      i=0
      else
      i=1
      endif
      return
      end 
