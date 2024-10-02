      subroutine enclab (n,kap,*)
c
c  this program encodes the principal and angular quantum numbers 'n'
c  and 'kap' into an orbital label 'lab'
c
      implicit character*1(l)
      common/orblab/lab(4)
      dimension l1(10),l2(10),l3(10),l4(2)
      data l1/' ','1','2','3','4','5','6','7','8','9'/,
     1     l2/'0','1','2','3','4','5','6','7','8','9'/,
     2     l3/'s','p','d','f','g','h','i','j','k','l'/,
     3     l4/' ','*'/,lx/'&'/,nx1/10/
c
      i=1
      if(n.lt.1.or.n.gt.99)  then
          write(6,1000) n
 1000     format( ' Enclab: principal quantum number out of range',
     1            ' n =',i20)
          return 1
      endif
      n1=n/10+1
      n2=n-n1*10+11
      lab(1)=l1(n1)
      lab(2)=l2(n2)
      n3=iabs(kap)
      n4=1
      if(kap.le.0) go to 20
      n3=n3+1
      n4=2
 20   if(n3.lt.1.or.n3.gt.min0(n,nx1)) then
          write(6,1010) kap
 1010     format( ' Enclab: angular quantum number out of range',
     1            ' kap =',i20)
          return 1
      endif
      lab(3)=l3(n3)
      lab(4)=l4(n4)
      return
      end
