      implicit real*8 (a-h,o-z)
      PARAMETER (NSORB=600)
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      open (unit=14,file='SGC.CON')

      max_orb=88
      nav1=max_orb
      lvmax=4
      key_en=1
      ncore=1

      call bass_inp
      Kt=1

      Kval=key_en
      Khot=0
      C_SMS=0.00d0
      klow=1
      write(14,55) max_orb,lvmax,Kt,Kval,Khot,C_SMS,Klow
 55   format(1X,' Nmax=',I3,'lmax=',I1,'  kt=',I2,'  kval=',I2,
     >         '  Khot=',I2,' C_SMS=',F9.6,' Klow=',I1)

 95   format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))
      sigma=0.d0
      ic=0
      do 531 i=ncore+1,max_orb
      do 532 j=ncore+1,max_orb
      iio=i+j
      if (iio.le.nav1*2) then
      if (i.le.j) then
      if (kor(i).eq.kor(j)) then
      ka=kor(i)

      CALL klj(ka,kapa,la,ja,ind,n0a)
      if (la.le.lvmax) then

      ni=nor(i)-n0a
      nj=nor(j)-n0a
13    format (i7,5i5,2f14.8)
      ic=ic+1
      write (14,95) ic,i,j,sigma,sigma,sigma


****************************************************
      endif
      endif
      endif
      endif
532   continue
531   continue
************needed extra
      write (14,*) '  screening coefficients'
      do 777 k=0,9
      write (14,10) k,1.d0
10    format (5X,i2,5x,f8.5)
777   continue
      write (14,*) ' >>>>>>>>> END <<<<<<<<<'
********************************************

      write (*,*) 'Total restricted Sigma 1 count          = ', ic

*********************
      close (14)

      end


      subroutine bass_inp
      implicit real*8 (a-h,o-z)
      PARAMETER (NSORB=600)
      character*1 nsign
      character*2 line
      character*7 a
      dimension nsign(nsorb),n1(nsorb),n2(nsorb),l(nsorb)
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      open (unit=1,file='BASS.INP')
      read (1,'(//)')
      read (1,'(5x,i5)') nso
      write (*,*) nso
      do 1 i=1,50
       read (1,'(2x,a2)') line
       if (line.eq.'--') goto 2
1     continue
2     continue
3     format (6(4x,a1,i1,1x,i1,i1,2x))
5     format (a1,i1,1x,i1,i1)
****************** read in core **************************
      read (1,3) ((nsign(i),n1(i),n2(i),l(i)),i=1,nso)
****************** read set ******************************
      j=nso
      do 6 i=1,nsorb-nso
      read (1,'(4x,a7)',END=555) a
      if (a(3:3).eq.'.') then
      j=j+1
      open(unit=21,file='fort.2')
      write (21,'(a7)') a
      rewind(21)
      read (21,5) nsign(j),n1(j),n2(j),l(j)
      close(21)
      endif
6     continue
555   norb=j
      if (norb.gt.nsorb) then
      write (*,*)
     * 'NUMBER OF ORBITALS IN BASS.INP EXCEEDS NSORB = ',nsorb
      write (*,*) 'STOP'
      stop
      endif
**********************************************************
      do 4 i=1,norb
      nor(i)=10*n1(i)+n2(i)
      if (nsign(i).eq.' ') kor(i)=-(l(i)+1)
      if (nsign(i).eq.'-') kor(i)=l(i)
c      write (*,'(i4,2x,i2,1x,i3)') i,nor(i),kor(i)
4     continue
      close(1)
**********************************************************
      return
      end

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
