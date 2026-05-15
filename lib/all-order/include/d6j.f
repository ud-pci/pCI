      function d6j(ia,ib,ie,id,ic,ig)
      implicit doubleprecision(a-h,o-z)
      parameter(NMAX=100)
 
c   calculate the 6-j sysmbol:  ( ia ib ie )
c                               ( id ic ig )
c    this version is the old version used in rrpa
c    we use it here to accomodate the high values
c    of angular momentum that occur
c 
      common/clebsc/p(NMAX)
      
  
 1000 l1= iabs(ib-ie)
      if(ia-l1) 4000,12,12
   12 if(ia-ib-ie) 13,13,4000
   13 l1= iabs(ic-ig)
      if(ia-l1) 4000,14,14
   14 if(ia-ic-ig) 15,15,4000
   15 l1= iabs(id-ic)
      if(ie-l1) 4000,16,16
   16 if(ie-id-ic) 17,17,4000
   17 l1= iabs(id-ig)
      if(ib-l1) 4000,18,18
   18 if(ib-id-ig) 19,19,4000
   19 l1=ia+ib+ie-((ia+ib+ie)/2)*2
      if(l1) 4000,20,4000
   20 l1=ia+ic+ig-((ia+ic+ig)/2)*2
      if(l1) 4000,21,4000
   21 l1=ic+id+ie-((ic+id+ie)/2)*2
      if(l1) 4000,22,4000
   22 l1=ib+id+ig-((ib+id+ig)/2)*2
      if(l1) 4000,23,4000
 
   23 continue
 
c  calculate nonzero d6j using page 99 of edmonds' book
 
      i1=(ia+ib-ie)/2+1
      i2=(ia-ib+ie)/2+1
      i3=(-ia+ib+ie)/2+1
      i4=(ia+ib+ie+2)/2+1
      if(i4-NMAX) 26,26,5000
   26 continue
      d11=p(i1)
      d12=p(i2)
      d13=p(i3)/p(i4)
      d11=sqrt(d11)
      d12=sqrt(d12)
      d13=sqrt(d13)
      d1=d11*d12*d13
      i1=(ia+ic-ig)/2+1
      i2=(ia-ic+ig)/2+1
      i3=(-ia+ic+ig)/2+1
      i4=(ia+ic+ig+2)/2+1
      if(i4-NMAX) 27,27,5000
   27 continue
      d11=p(i1)
      d12=p(i2)
      d13=p(i3)/p(i4)
      d11=sqrt(d11)
      d12=sqrt(d12)
      d13=sqrt(d13)
      d2=d11*d12*d13
      i1=(id+ib-ig)/2+1
      i2=(id-ib+ig)/2+1
      i3=(-id+ib+ig)/2+1
      i4=(id+ib+ig+2)/2+1
      if(i4-NMAX) 28,28,5000
   28 continue
      d11=p(i1)
      d12=p(i2)
      d13=p(i3)/p(i4)
      d11=sqrt(d11)
      d12=sqrt(d12)
      d13=sqrt(d13)
      d3=d11*d12*d13
      i1=(id+ic-ie)/2+1
      i2=(id-ic+ie)/2+1
      i3=(-id+ic+ie)/2+1
      i4=(id+ic+ie+2)/2+1
      if(i4-NMAX) 29,29,5000
   29 continue
      d11=p(i1)
      d12=p(i2)
      d13=p(i3)/p(i4)
      d11=sqrt(d11)
      d12=sqrt(d12)
      d13=sqrt(d13)
      d4=d11*d12*d13
      ss=0.
      k1=ia+ib+ie
      k2=ia+ic+ig
      k3=id+ib+ig
      k4=id+ic+ie
      k5=ia+ib+id+ic
      k6=ib+ie+ic+ig
      k7=ie+ia+ig+id
      kmax=min0(k5,k6,k7)
      kmin=max0(k1,k2,k3,k4)
      item=max0(kmin,k5,k6,k7)
      if(item-NMAX) 30,30,5000
   30 continue
      l1=kmin/2-1
      sign1=1.0
      if(mod(l1,2).ne.0) sign1=-1.0
      k=kmin
   32 continue
      sign1=-sign1
      in=(k+2)/2+1
      i1=(k-k1)/2+1
      i2=(k-k2)/2+1
      i3=(k-k3)/2+1
      i4=(k-k4)/2+1
      i5=(k5-k)/2+1
      i6=(k6-k)/2+1
      i7=(k7-k)/2+1
      t=sign1*p(in)/p(i1)
      t=t/p(i2)
      t=t/p(i3)
      t=t/p(i4)
      t=t/p(i5)
      t=t/p(i6)
      t=t/p(i7)
      ss=ss+t
      k=k+2
      if(k-kmax) 32,32,34
   34 continue
      d6j=ss*d1*d2*d3*d4
      return
 
 4000 d6j=0.
      return
 
 5000 write(6,3000) ia,ib,ie,id,ic,ig
 3000 format('failure in calculating 6-j  {',2(i3,','),i3,'}'/
     x       '                            {',2(i3,','),i3,'}')
      d6j = 0.
      return
      end
