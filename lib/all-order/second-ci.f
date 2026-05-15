******************************************************************************
*
*      Version 5.0 Reads in all-order input if such exists and
*                replaces needed parts of the SGC.CON and SCRC.CON files
*                  May 9, 2008
*
******************************************************************************
*      Version 4.1 Dependent subroutines (rint,yfun,yint,d6j,libD)
*                  are replaced to make a consistent package, inidat was added
*      May 5, 2008
*
*      second-ci.f Version 4.0
*      THIS VERSION REQUIRES BASS.INP
*      COMPLETED: May 3, 2008
*________________________________________________
*      Fast code for second-order Sigma_1
*      Completed Version 1.0 - November 27, 2007
*      Verion 2.0 April 17, 2008
*      Derivatives for sigma 1 added, reading MK basis set is added ,
*      input changed
*      Version 3.0 Writing second-order Sigma_2
*      Version 4.0 Calculates Sigma 2 in correct order with derivatives
*      kval =1 and kval =2 options added
*      parameters are put in a separate file
*******************************************************************************
*########################
*# Explanation of input #
*########################
*
* 4  # number of core shells ncore
* 1 -1 # n and kappa of the core shells
* 2 -1
* 2  1
* 2 -2
* 25 5  # nmax and lmax in correlation diagram summations
* 93      # Nmax  (max_orb) # of orbitals (from BASS.INPT) for how many Sigmas to calculate
* 3   # Lmax (lvmax) for how many Sigmas to calculate, superceeds Nmax
* 9   # Kmax (kvmax) for higest partial wave for Sigma2 S_Kmax(ijkl)
* 93    33 # nav1 nav, maximum average of orbitals to calculate Sigma for
*       nav1=(i+j)/2 for Sigma(ij) and nav=(i+j+k+l)/4 for Sigma_K(ijkl)
* 1     # key_allorder = 1 if needs to read all-order input, 0 otherwise
* 1   # kval (key_en): key for energies, DHF if kval=0, lowest for a given kappa for kval=1
*        if kval=2, the energies are inputted in the format below
*        the program expects that all energies up to lvmax will be inputted
* Example:
* 2
* 0  -0.28000
* 1  -0.22000  -0.22000
* 2  -0.31000  -0.31000
* 3  -0.13000  -0.13000
*********************************************************************************

      implicit real*8 (a-h,o-z)
      character*2 lab1,lab2,lab3,lab4
      character*3 key_inv
      include "second.par"
      DIMENSION t(nkx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      DIMENSION g(nhf,nx),f(nhf,nx),e(nl)
      DIMENSION wco(ns),nco(ns),kco(ns),mco(ns)
      common /sigma1st/ sigma1(nxx1,nxx1,nk1),sigmad(nxx1,nxx1,nk1)
      COMMON /eninput/ ecorr(nk1)
      do 531 m=1,nxx1
        do 532 n=1,nxx1
       do 533 ind=1,nk1
            sigma1(m,n,ind)=0.d0
            sigmad(m,n,ind)=0.d0
533      continue
532     continue
531   continue



c      open (15,file='inf.vw',form='FORMATTED',status='OLD',err=710)
      OPEN (unit=1,file='hfspl.1',form='unformatted',
     >      status='OLD',err=720)
      OPEN (unit=2, file='hfspl.2', form='unformatted',
     *access='direct',recl=nrec,status='OLD',err=730)
      READ (1) nmax,k,t,lmax
      READ (1) jmax,(wco(i),nco(i),kco(i),mco(i),i=1,jmax)
      READ (1) r,rp,rpor,h,max
      nmax1=nmax
      lmax1=lmax
! MGK>
      read(1,end=137,err=137) nmaxx
      goto 138
 137  write(*,*) ' No record for nmaxx. Use autofill.'
      do i=1,2*lmax+1
        nmaxx(i)=nmax
      end do
 138  continue
! MGK<

      call inidat
*******  find last point for which r(ip) <= t(n+1) nv

      mp = 0
      t0 = (1d0 + 10*EPSILON(t0))*t(nmax+1)
      do 1009 ip = 1,max
          if( r(ip).le.t0 ) mp = ip
1009  continue
      max = mp
      write(6,1030) r(2),r(max),max
1030  format(/' r(2) =',1pe15.6,2x,'r(max) =',1pe15.6,
     &      '     max =',i6/)

      write (*,*) max
*******************READ IN SPLINES *********************************

       do 87 ind=1,2*lmax+1
         READ (2,rec=ind) g,f,e
         do 88 i=1,nmaxx(ind)
           DO 89 j=1,max
             ff(j,i,ind)=f(j,i)
             gg(j,i,ind)=g(j,i)
89          CONTINUE
             ee(i,ind)=e(i+nmax)
88         continue
87     continue

**************  READ IN CORE  ***************************************
      inmax=0
      read (5,*) ncore
      DO 1 i=1,ncore
        READ (5,*) no(i),ko(i)
        CALL klj (ko(i),k,l,jj,ind,n0o)
       if (ind.gt.inmax) then
       inmax=ind
       endif
        DO 2 j=1,max
          fo(j,i)=ff(j,no(i)-n0o,ind)
          go(j,i)=gg(j,no(i)-n0o,ind)
2      CONTINUE
        eo(i)=ee(no(i)-n0o,ind)
        write (*,'(2i3,f12.7)') no(i),ko(i),eo(i)
1     CONTINUE

      read (5,*) nmax,lmax

      nnj=2*lmax+1

      if (nnj.gt.nk) then
      write (*,*) ' Parameter NK is exceeded, STOP'
      stop
      endif

      if (nmax.gt.nxx) then
      write (*,*) ' Parameter NXX is exceeded, STOP'
      stop
      endif

      if (ncore.gt.nn) then
      write (*,*) ' Parameter NN is exceeded, STOP'
      stop
      endif

      do i=1,2*lmax+1
        if (nmaxx(i).GT.nmax) nmaxx(i)=nmax
        write(*,*) ' PW=',i,' nmax=',nmaxx(i)
      end do

c      read (5,*) nvmax,lvmax
******************************************************
*       This key determines if sigma_1 matrix should
*       be inverted. Original reason for inverting the
*       matrix was the order of x,y indexes in the original SD
*       code. The matrix was inverted to bring the upper half
*       in closer agreement with WRJ code which set energy to lowest
c       read (5,*) key_inv
      key_inv='no '
*************************************************************
*********************************************
*      This key allows to cut of lmax in k sum for comparison with WRJ
*
c       read (5,*) key_k
c       if (key_k.eq.1) then
c       write (*,*) 'Partial waves are restricted for k at lmax'
c       endif
      key_k=0
c>>>mgk
c instead of reading BASS.INP the code now reads HFD.DAT
c this allows to work with Johnson's basis sets, where BASS.INP
c is not used and HFD.DAT is formed with bas_wj from hfspl.2.
c****************** read in BASS.INP *******************
c      call bass_inp
      call init_hfd     ! reads HFD.DAT to get list of orbitals
c<<<mgk
***************************************************************
********* READ MAX NUMBER OF ORBITALS TO CALCULATE SIGMA 1 FOR
*****************************************************************
      read (5,*) max_orb
      read (5,*) lvmax
      read (5,*) kvmax
      nnj=2*lvmax+1

      if (nnj.gt.nk1) then
      write (*,*) ' Parameter NK1 is exceeded, STOP'
      stop
      endif


      read (5,*) nav1, nav

******************************************************
*   This key determines which energy to put in
*    0 - as in SD code
*    1 - lowest n for a given kappa
*    2 - reads input below
*********************************************************
      read (5,*) key_allorder
      read (5,*) key_en
      if (key_en.eq.0) then
      write (*,*) 'The HF energy is used, kval = 0'
      endif
      if (key_en.eq.1) then
       write (*,*) 'The energy is set to lowest n energy for each ind'
      endif
      if (key_en.eq.2) then
      write (*,*) 'The energy is set to the inputted data'
      read (5,*) j,ecorr(1)
      i=2
      if (lvmax.ne.0) then
       do 78 j=1,lvmax
         read  (5,*) ii,ecorr(i),ecorr(i+1)
         i=i+2
78       continue
      endif
      endif
      call orb(max_orb,lvmax)
*****************************************************************
c      goto 999
****************************
      write (*,*) nmax,lmax,nmax1,lmax1

      call count_xmnab
      call count_xmnra
      call count_sigma2code(max_orb,key4,nav,kvmax)

      t0 = mclock()
      call xmnab
      t1 = mclock()
      dt = (t1 - t0)/1000.0d0
c      write(6,1090) dT
 1090 format('Xmnab time :',f8.2,'sec')

      t0 = mclock()
      call xmnra
      t1 = mclock()
      dt = (t1 - t0)/1000.0d0
      write(6,1091) dT
 1091 format('Xmnra time :',f8.2,'sec')
*************** CALCULATION OF SIGMA 1 *************
      t0 = mclock()
      call term1_s1(key_en)
      call term2_s1(key_en,key_k)
      t1 = mclock()
      dt = (t1 - t0)/1000.0d0
c      write(6,1093) dT
 1093 format('Sigma 1 time :',f8.2,'sec')

* >>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>
      if (key_allorder.eq.1) then
        write (*,*) 'READING ALL-ORDER INPUT'
        call read_allorder(key_en)
        call merge_sigma1(max_orb)
      endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call write_out(max_orb,lvmax,key_en,nav1)

*************** CALCULATION OF SIGMA 2 *************
      write (*,*)
      key4=1
999   t0 = mclock()
      call sigma2code(max_orb,key4,key_en,nav,kvmax,key_allorder)
      t1 = mclock()
      dt = (t1 - t0)/1000.0d0
      write(6,1094) dT
 1094 format('Sigma 2 time :',f8.2,'sec')
******************************************************
      stop
 720  write(*,*) ' No file "hfspl.1" ifound '
      stop
 730  write(*,*) ' No file "hfspl.2" found '
      end
********
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>

      subroutine merge_sigma1(max_orb)
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxx/nmaxx(nx)
******** ALL_ORDER **********************************
      COMMON /val_all/ nv_all(nz),kv_all(nz),nval_all
      COMMON /energy_all/ evt_all(nz)
      COMMON /s1_all/ rsigma(nz,nx),drsigma(nz,nx)
*****************************************************
      common /sigma1st/ sigma1(nxx1,nxx1,nk1),sigmad(nxx1,nxx1,nk1)
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /energy/ en(nxx1,nxx1,nk1)
**********************************************************
 95   format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))

      ic=0
      iall=0
      do 531 i=ncore+1,max_orb
      do 532 j=ncore+1,max_orb
      if (i.le.j) then
      if (kor(i).eq.kor(j)) then
      ka=kor(i)

      CALL klj(ka,kapa,la,ja,ind,n0a)

      ni=nor(i)-n0a
      nj=nor(j)-n0a
13    format (i7,5i5,2f14.8)
      ic=ic+1
*************** find all-order sigma *******************
      do 12 ii=1,nval_all
      if (nv_all(ii).eq.nor(i).and.kv_all(ii).eq.ka) then
      r=rsigma(ii,nj)
      dr=drsigma(ii,nj)
      enr=evt_all(ii)

c   write (*,95) ic,i,j,sigma1(nj,ni,ind),sigmad(nj,ni,ind),
c     *en(nj,ni,ind)
c   write (*,95) ic,i,j,r,dr,enr
      sigma1(nj,ni,ind)=r
      sigmad(nj,ni,ind)=dr
      en(nj,ni,ind)=enr

      iall=iall+1
c      write (*,*)
      endif
12    continue
********************************************************
      endif
      endif
532   continue
531   continue
      write (*,*) 'Total unrestricted Sigma 1 count        = ', ic
      write (*,*) 'Number of Sigma 1 replaced by all-order = ', iall

      return
      end

      subroutine read_allorder(key_en)
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxx/nmaxx(nx)
      COMMON /s2_all/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ich_all
      COMMON /val_all/ nv_all(nz),kv_all(nz),nval_all
      DIMENSION nmaxx9(nx)
      COMMON /i_all/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /energy_all/ evt_all(nz)
      COMMON /ds2_all/ dxvw(nxx,nxx,0:NVWH)
      COMMON /s1_all/ rsigma(nz,nx),drsigma(nz,nx)

      OPEN (unit=17,file='sigma1')
      OPEN (unit=7,file='pair.vw',form='unformatted')
      read (7) ich_all
      read (7) nmaxx9
      read (7) ipm,ipn
      read (7) nv_all,kv_all,nval_all
      read (7) ipvw
      read (7) kval,evt_all
      if (kval.ne.key_en) then
      write (*,*) 'kval mismatch with all-order input, STOP'
      stop
      endif
      do 30 i=1,nx
      if (nmaxx(i).ne.nmaxx9(i)) then
      write (*,*) 'Mismath in basis set with all-order data, STOP'
      stop
      endif
30    continue

      do 27 i=1,ich_all
        indn=ipn(i)
        indm=ipm(i)
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)
        do 37 n=1,nmaxn
          read (7) (xvw(m,n,i),m=1,nmaxm)
37      continue

        do 327 n=1,nmaxn
          read (7) (dxvw(m,n,i),m=1,nmaxm)
327       continue

27    continue
      close (7)
      do 1 i=1,nval_all
      nvv=nv_all(i)
      kvv=kv_all(i)
      CALL klj(kvv,kapv,lv,jv,indv,n0v)
      indm=indv
      CALL indk1(indm,km,kapm,lm,jm,m0)
      call st (km,nnm)
      do 200 m=nnm,nmaxx(indm)
      read (17,9) nvv1,kvv1,mm1,rsigma(i,m),drsigma(i,m)
9     format (3i4,2e20.8)
      mmm=m+m0
      if (nvv1.ne.nvv.or.kvv1.ne.kvv.or.mmm.ne.mm1) then
      write (*,*) 'mismatch in all-order sigma read in, STOP'
      stop
      endif
200   continue
1     continue
      close(17)
      return
      end

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine orb(max_orb,lvmax1)
      implicit real*8 (a-h,o-z)
      include "second.par"
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      common /nmaxxv/ nmaxxv(nx),lvmax
      lvmax=lvmax1
************
      icount=0
      do 2 ind=1,2*lvmax+1
        nmaxxv(ind)=0
        CALL indk1(ind,kv,kapv,lv,jv,n0v)
        do 1 i=1,max_orb
           if (kor(i).eq.kv) then
           if (nor(i).gt.nmaxxv(ind)) then
            nmaxxv(ind)=nor(i)
           endif
         endif
1       continue
90      format ('ind = ',2i3,'  nmax = ',i3)
        if (nmaxxv(ind).ne.0) then
         nmaxxv(ind)=nmaxxv(ind)-n0v
c>>>>>>>>>>>>>>>>>>>>>> nxx1 limit added on valence >>>>>>>
         n99=nmaxxv(ind)
         if (n99.gt.nxx1) then 
            write (*,*) 'Replaced by nxx1'
            nmaxxv(ind)=nxx1
         endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        endif
        write (*,90) ind,kv,nmaxxv(ind)
        icount=icount+nmaxxv(ind)
2     continue
c      write (*,*) 'icount',icount
      if (icount.ne.max_orb) then
c      write (*,*) 'Counting mismatch in orb, stop'
c      stop
      endif
      return
      end

      subroutine bass_inp
      implicit real*8 (a-h,o-z)
      include "second.par"
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

c>>>mgk
C     =================================================
      subroutine init_hfd
      implicit real*8 (a-h,o-z)
      include "hfd.par"  !### check directory for hfd.par!
      include "second.par"
       common /bassinp/ nor(nsorb),kor(nsorb),norb
       common /ipmr/ipmr
       dimension P(IP6),Q(IP6),P1(IP6),Q1(IP6),PQ(4*IP6)
       logical longbasis
       dimension IQN(4000) ! dimension should be .GE. 4*IPs, IPs from conf.par
       equivalence (IQN(1),PQ(21))
       equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),
     >        (P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c small number:
        c1=0.01d0
        call recunit
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(12,file='HFD.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*IPmr,err=700)
        call ReadF (12,1,P,Q,2)
        call ReadF (12,3,P1,Q1,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        norb =PQ(2)+C1
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          do ni=1,norb
            nor(ni)=IQN(4*ni-3)
            kor(ni)=IQN(4*ni-1)
**
         write (*,'(3i5)') ni, nor(ni),kor(ni)
**
          end do
        else
          if=20
          do ni=1,norb
            if=if+1
            nor(ni)=PQ(if)+c1
            if=if+4
            c2=dsign(c1,PQ(if))
            kor(ni)=PQ(if)+c2
            if=if+1
**
         write (*,'(3i5)') ni, nor(ni),kor(ni)
**
          end do
        end if
        close(12)
       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
75      format(/2X,'file HFD.DAT is absent'/)
       stop
      end
c     =================================================
      Subroutine ReadF (Kan,record,V1,V2,nrec)
      include "hfd.par"
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C        FORTRAN-77, MS-DOS VERSION
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       integer record,nrec
       real*8 V1(IP6),V2(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        read(Kan,rec=nr1) (V1(I),i=1,ii)
        if (nrec.EQ.2) read(Kan,rec=nr2) (V2(I),i=1,ii)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       return
      end
c<<<mgk


      subroutine write_out(max_orb,lvmax,key_en,nav1)
      implicit real*8 (a-h,o-z)
      include "second.par"
      common /sigma1st/ sigma1(nxx1,nxx1,nk1),sigmad(nxx1,nxx1,nk1)
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /energy/ en(nxx1,nxx1,nk1)
      open (unit=14,file='SGC.CON')
      Kt=1
      Kval=key_en
      Khot=0
      C_SMS=0.00d0
      klow=1
      write(14,55) max_orb,lvmax,Kt,Kval,Khot,C_SMS,Klow
 55   format(1X,' Nmax=',I3,'lmax=',I1,'  kt=',I2,'  kval=',I2,
     >         '  Khot=',I2,' C_SMS=',F9.6,' Klow=',I1)

 95   format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))

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
      write (14,95) ic,i,j,sigma1(nj,ni,ind),sigmad(nj,ni,ind),
     *en(nj,ni,ind)


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
      return
      end



      subroutine llabel(k,lab)
      character*2 lab
      if (k.eq.-1) lab='s '

      if (k.eq.1)  lab='p*'
      if (k.eq.-2) lab='p '

      if (k.eq.2)  lab='d*'
      if (k.eq.-3) lab='d '

      if (k.eq.3)  lab='f*'
      if (k.eq.-4) lab='f '

      if (k.eq.4)  lab='g*'
      if (k.eq.-5) lab='g '

      if (k.eq.5)  lab='h*'
      if (k.eq.-6) lab='h '
      return
      end


*****************************************************************************************
      subroutine xmnab
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxxv/ nmaxxv(nx),lvmax
      COMMON /radmxab/ ip1(nk,nk,nn,nn,0:kk),x1(nxx,nxx,0:nch1),ichan1
      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf),gc(nhf),fc(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),uy(nhf,nxx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/ nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0
      do 531 i=0,nch1
        do 532 m=1,nxx
       do 533 nx5=1,nxx
            x1(m,nx5,i)=0.d0
533      continue
532     continue
531   continue

      do 661 i5=0,kk
       do 662 i4=1,nn
         do 663 i3=1,nn
         do 664 i2=1,nk
          do 665 i1=1,nk
                ip1(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

      do 4 j=1,ncore
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        DO 13 il=1,max
          gb(il)=go(il,j)
          fb(il)=fo(il,j)
13      CONTINUE
c       write (*,*) nb,kb,gb(100)
*****************************************************
        do 2 indx=1,2*lmax+1
          CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
          call st (kx5,nnx)
****************************************************************
*             Calculate  Xk(mxab)
****************************************************************
          kmin=iabs(kapb-kapx)
          kmax=kapx+kapb-1
          DO 11 k=kmin,kmax
            call odd (lb+lx1+k,i1)
            if (i1.eq.0) goto 555
*********************************************************************
          do 29 nx5=nnx,nmaxx(indx)
              do 19 ij=1,max
                v(ij)=gb(ij)*gg(ij,nx5,indx)+fb(ij)*ff(ij,nx5,indx)
19            continue
              call yfun(v,u,k,max,*901)
              do 122 ij=1,max
                uy(ij,nx5)=u(ij)
122           continue
29          continue
*********************************************************************
            do 14 i=1,ncore
              ka=ko(i)
              na=no(i)
              CALL klj(ka,kapa,la,ja,inda,n0a)
              DO 12 il=1,max
                ga(il)=go(il,i)
                fa(il)=fo(il,i)
12            CONTINUE
c       write (*,*) na,ka,ga(100)

              do 23 indm=1,2*lmax+1
                CALL indk1(indm,km,kapm,lm,jm,n0m)
                call st (km,nnm)
***********************************************************************
                call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                call odd (la+lm+k,i2)
                if (i1*i2.eq.0) goto 700

                index=index+1
                ip1(indm,indx,i,j,k)=index
              if (index.ge.nch1) then
               write (*,*) 'NUMBER OF Xmxab CHANNELS IS EXCEEDED'
               stop
              endif
                c=((-1)**(k))*s(k,km,ka)*s(k,kx5,kb)
                do 61 m=nnm,nmaxx(indm)
                  do 6  nx5=nnx,nmaxx(indx)
                    do 20 ij=1,max
                      vv(ij)=(ga(ij)*gg(ij,m,indm)+
     *                        fa(ij)*ff(ij,m,indm))*rp(ij)*uy(ij,nx5)
20                  continue
                    x1(m,nx5,index)=c*rint(vv,1,max,11,h)
c   write (*,'(4e12.2)') ga(100),gg(100,m,indm),rp(100),uy(100,nx5)
6                 continue
61              continue
700             continue
23            continue
14          continue
555         continue
11        continue
*****************************************************************
2       continue
4     continue

      ichan1=index

901   continue
      return
      end



*****************************************************************************************
      subroutine xmnra
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxxv/ nmaxxv(nx),lvmax
      COMMON /radmnxa/ ip2(nk,nk,nk,nn,0:kk),
     * x2(nxx,nxx,nxx,0:nch2),ichan2
      DIMENSION ga(nhf),fa(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),uy(nhf,nxx)
      DIMENSION vv1(nhf,nxx,nxx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0
      do 531 i=0,nch2
        do 532 m=1,nxx
        do 533 n=1,nxx
          do 534 nx5=1,nxx
              x2(m,n,nx5,i)=0.d0
534         continue
533        continue
532     continue
531   continue

      do 661 i5=0,kk
       do 662 i4=1,nn
         do 663 i3=1,nk
         do 664 i2=1,nk
          do 665 i1=1,nk
                ip2(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

************************************************************************

      do 4 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 13 il=1,max
          ga(il)=go(il,i)
          fa(il)=fo(il,i)
13      CONTINUE
*****************************************************
        do 23 indn=1,2*lmax+1
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (kn,nnn)
****************************************************************
*             Calculate  Xk(mnxa)
****************************************************************
          kmin=iabs(kapa-kapn)
          kmax=kapa+kapn-1
          DO 11 k=kmin,kmax
            call odd (la+ln+k,i1)
            if (i1.eq.0) goto 555
*********************************************************************
          do 29 n=nnn,nmaxx(indn)
              do 19 ij=1,max
                v(ij)=ga(ij)*gg(ij,n,indn)+fa(ij)*ff(ij,n,indn)
19            continue
              call yfun(v,u,k,max,*901)
              do 122 ij=1,max
                uy(ij,n)=u(ij)*rp(ij)
122           continue
29          continue
*********************************************************************
            do 238 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,n0m)
              call st (km,nnm)
              do 2 indx=1,2*lmax+1
                CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
                call st (kx5,nnx)

***********************************************************************
                call trgi(iabs(kapm-kapx),(kapm+kapx-1),k,i1)
                call odd (lx1+lm+k,i2)
                if (i1*i2.eq.0) goto 700

                index=index+1
                ip2(indm,indn,indx,i,k)=index
              if (index.ge.nch2) then
               write (*,*) 'NUMBER OF Xmnxa CHANNELS IS EXCEEDED, STOP'
               stop
              endif
                c=((-1)**(k))*s(k,km,kx5)*s(k,kn,ka)
**************************************************************************
                do 76  nx5=nnx,nmaxx(indx)
                do 761 m=nnm,nmaxx(indm)
                    do 720 ij=1,max
                  vv1(ij,m,nx5)=gg(ij,nx5,indx)*gg(ij,m,indm)+
     *                   ff(ij,nx5,indx)*ff(ij,m,indm)
720                  continue
761                 continue
76              continue
************************************************************************
                  do 63 n=nnn,nmaxx(indn)
                do 6  nx5=nnx,nmaxx(indx)
                do 61 m=nnm,nmaxx(indm)
                    do 20 ij=1,max
                  vv(ij)=vv1(ij,m,nx5)*uy(ij,n)
20                  continue
                    x2(m,n,nx5,index)=c*rint(vv,1,max,11,h)
61                 continue
6                 continue
63              continue
************************************************************************

700             continue
2              continue
238          continue
555         continue
11        continue
*****************************************************************
23      continue
4     continue

      ichan2=index

901   continue
      return
      end
*****************************************************************************************
      subroutine count_xmnab
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/ nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0

      do 4 j=1,ncore
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
*****************************************************
        do 2 indx=1,2*lmax+1
          CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
          call st (kx5,nnx)
****************************************************************
*             Calculate  Xk(mxab)
****************************************************************
          kmin=iabs(kapb-kapx)
          kmax=kapx+kapb-1
          DO 11 k=kmin,kmax
            call odd (lb+lx1+k,i1)
            if (i1.eq.0) goto 555
*********************************************************************
            do 14 i=1,ncore
              ka=ko(i)
              na=no(i)
              CALL klj(ka,kapa,la,ja,inda,n0a)

              do 23 indm=1,2*lmax+1
                CALL indk1(indm,km,kapm,lm,jm,n0m)
                call st (km,nnm)
***********************************************************************
                call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                call odd (la+lm+k,i2)
                if (i1*i2.eq.0) goto 700
      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded in xmnab, STOP'
      stop
      endif

                index=index+1
700             continue
23            continue
14          continue
555         continue
11        continue
*****************************************************************
2       continue
4     continue

      ichan1=index

      write (*,*)
      write (*,*) 'PRE-COUNTING ALL CHANNELS'
      write (*,*)
      if (index.le.NCH1) then
      write (*,122) index,nch1
      else
      write (*,124) index,nch1
      write (*,*) 'Parameter NCH1 is exceeded, STOP'
      stop
      endif
122   format ('Number of X_k(mnab) channels      = ',i7,
     *' < NCH1 = ',i7,' OK')

124   format ('Number of X_k(mnab) channels      = ',i7,
     *' > NCH1 = ',i7)

      return
      end



*****************************************************************************************
      subroutine count_xmnra
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0
************************************************************************

      do 4 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
*****************************************************
        do 23 indn=1,2*lmax+1
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (kn,nnn)
****************************************************************
*             Calculate  Xk(mnxa)
****************************************************************
          kmin=iabs(kapa-kapn)
          kmax=kapa+kapn-1
          DO 11 k=kmin,kmax
            call odd (la+ln+k,i1)
            if (i1.eq.0) goto 555
*********************************************************************
            do 238 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,n0m)
              call st (km,nnm)
              do 2 indx=1,2*lmax+1
                CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
                call st (kx5,nnx)

***********************************************************************
                call trgi(iabs(kapm-kapx),(kapm+kapx-1),k,i1)
                call odd (lx1+lm+k,i2)
                if (i1*i2.eq.0) goto 700
      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded xmnra, STOP'
      stop
      endif

                index=index+1
************************************************************************

700             continue
2              continue
238          continue
555         continue
11        continue
*****************************************************************
23      continue
4     continue

      ichan2=index
      if (index.le.NCH2) then
      write (*,122) index,nch2
      else
      write (*,124) index,nch2
      write (*,*) 'Parameter NCH2 is exceeded, STOP'
      stop
      endif
122   format ('Number of X_k(mnra) channels      = ',i7,
     *' < NCH2 = ',i7,' OK')

124   format ('Number of X_k(mnra) channels      = ',i7,
     *' > NCH2 = ',i7)

      return
      end

*****************************************************************************************
      subroutine term1_s1(key_en)
      implicit double precision (a-h,o-z)
      include "second.par"
      common /sigma1st/ sigma1(nxx1,nxx1,nk1),sigmad(nxx1,nxx1,nk1)
      common /nmaxxv/ nmaxxv(nx),lvmax
      COMMON /radmxab/ ip1(nk,nk,nn,nn,0:kk),x1(nxx,nxx,0:nch1),ichan1
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /eninput/ ecorr(nk1)

      do 531 m=1,nxx1
        do 532 n=1,nxx1
       do 533 ind=1,nk1
            sigma1(m,n,ind)=0.d0
            sigmad(m,n,ind)=0.d0
533      continue
532     continue
531   continue
c      write (*,*) '1 key',key_en

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
      ea=eo(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
        eb=eo(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
*****************************************************
          do 2 indx=1,2*lvmax+1
            CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
            call st (kx5,nnx)
            indy=indx
            CALL indk1(indy,ky,kapy,ly,jy,n0y)
            call st (ky,nny)
            do 24 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,n0m)
              call st (km,nnm)

              kmin=max0(iabs(kapm-kapa),iabs(kapy-kapb))
              kmax=min0((kapm+kapa-1),(kapy+kapb-1))

              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ly+lb+k,i2)
                if (i1*i2.eq.0) goto 555
      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                index1=ip1(indm,indy,i,j,k)
                c=-1.d0/((2.d0*k+1.d0)*(jy+1.d0))

                index2=ip1(indm,indx,i,j,k)
                do 61 m=nnm,nmaxx(indm)
                  do 62  nx5=nnx,nmaxxv(indx)
                   ex=ee(nx5,indx)

                    do 63  ny=nny,nmaxxv(indy)
                   ey=ee(ny,indy)
******************************************************
                  if (key_en.eq.0) then
                   e0=ee(ny,indy)
                  endif
                  if (key_en.eq.1) then
                   e0=ee(nny,indy)
                  endif
                  if (key_en.eq.2) then
                   e0=ecorr(indy)
                  endif
*****************************************************

                  del=1.d0/(ea+eb-ex-ee(m,indm)+(e0-ey))

                  tt=c*del*x1(m,ny,index1)*x1(m,nx5,index2)
                      sigma1(nx5,ny,indx)=sigma1(nx5,ny,indx)+tt
                 sigmad(nx5,ny,indx)=sigmad(nx5,ny,indx)-tt*del
c            write (*,'(4e12.3)') c,del,x1(m,ny,index1),x1(m,nx5,index2)
63                 continue
62                 continue
61              continue
**********************************************************************
              k1min=max0(iabs(kapm-kapb),iabs(kapa-kapx))
              k1max=min0((kapm+kapb-1),(kapa+kapx-1))
              DO 111 k1=k1min,k1max
                call odd (lm+lb+k1,i1)
                call odd (lx1+la+k1,i2)
                if (i1*i2.eq.0) goto 556
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                cc=c*(2.d0*k+1.d0)*d6j(jm,ja,2*k,jx,jb,2*k1)
                index2=ip1(indm,indx,j,i,k1)

                do 611 m=nnm,nmaxx(indm)
                  do 612  nx5=nnx,nmaxxv(indx)
                   ex=ee(nx5,indx)

                    do 613  ny=nny,nmaxxv(indy)
                   ey=ee(ny,indy)
******************************************************
                  if (key_en.eq.0) then
                   e0=ee(ny,indy)
                  endif
                  if (key_en.eq.1) then
                   e0=ee(nny,indy)
                  endif
                  if (key_en.eq.2) then
                   e0=ecorr(indy)
                  endif
*****************************************************

                  del=1.d0/(ea+eb-ex-ee(m,indm)+(e0-ey))



      tt=cc*del*x1(m,ny,index1)*x1(m,nx5,index2)
                      sigma1(nx5,ny,indx)= sigma1(nx5,ny,indx)+tt
                   sigmad(nx5,ny,indx)= sigmad(nx5,ny,indx)-tt*del

613                 continue
612                 continue
611              continue
556            continue
111           continue
*************************************************************
555             continue
11            continue
24          continue
2         continue
4        continue
*****************************************************************
1       continue

901   continue
      return
      end

*****************************************************************************************
      subroutine term2_s1(key_en,key_k)
      implicit double precision (a-h,o-z)
      include "second.par"
      COMMON /radmnxa/ ip2(nk,nk,nk,nn,0:kk),
     * x2(nxx,nxx,nxx,0:nch2),ichan2
      common /sigma1st/ sigma1(nxx1,nxx1,nk1),sigmad(nxx1,nxx1,nk1)
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /energy/ en(nxx1,nxx1,nk1)
      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf),gc(nhf),fc(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),uy(nhf,nxx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /eninput/ ecorr(nk1)

c      write (*,*) '2 key',key_en
      do 531 m=1,nxx1
      do 5311 n=1,nxx1

       do 533 ind=1,nk1
            en(m,n,ind)=0.d0
533      continue
5311   continue

531   continue

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
      ea=eo(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
*****************************************************
        do 2 indx=1,2*lvmax+1
          CALL indk1(indx,kx5,kapx,lx1,jx,n0x)
          call st (kx5,nnx)
          indy=indx
          CALL indk1(indy,ky,kapy,ly,jy,n0y)
          call st (ky,nny)
          do 24 indm=1,2*lmax+1

            CALL indk1(indm,km,kapm,lm,jm,n0m)
            call st (km,nnm)
            do 25 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
*************************************************************
             if (key_k.eq.0) then
             kmin=max0(iabs(kapm-kapx),iabs(kapn-kapa))
              kmax=min0((kapm+kapx-1),(kapn+kapa-1))
              endif
            if (key_k.eq.1) then
              kmin=0
              kmax=lmax
              endif
*************************************************************
              DO 11 k=kmin,kmax
                call trgi(iabs(kapm-kapx),(kapm+kapx-1),k,i3)
                call trgi(iabs(kapn-kapa),(kapn+kapa-1),k,i4)
                if (i3*i4.eq.0) goto 555

                call odd (lm+lx1+k,i1)
                call odd (ln+la+k,i2)
                if (i1*i2.eq.0) goto 555

      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif


                index1=ip2(indm,indn,indx,i,k)
                c=1.d0/((2.d0*k+1.d0)*(jy+1.d0))

                index2=ip2(indm,indn,indy,i,k)
                do 61 m=nnm,nmaxx(indm)
                  do 62 n=nnn,nmaxx(indn)
                    do 63  ny=nny,nmaxxv(indy)
******************************************************
                  if (key_en.eq.0) then
                   e0=ee(ny,indy)
                  endif
                  if (key_en.eq.1) then
                   e0=ee(nny,indy)
                  endif
                  if (key_en.eq.2) then
                   e0=ecorr(indy)
                  endif
*****************************************************

                  del=1.d0/(e0+ea-ee(m,indm)-ee(n,indn))
                    do 64  nx5=nnx,nmaxxv(indx)
      tt=c*del*x2(m,n,nx5,index1)*x2(m,n,ny,index2)
                   sigma1(nx5,ny,indx)=sigma1(nx5,ny,indx)+tt
                   sigmad(nx5,ny,indx)=sigmad(nx5,ny,indx)-tt*del

64                 continue
63                 continue
62                 continue
61              continue

              k1min=max0(iabs(kapm-kapa),iabs(kapn-kapy))
              k1max=min0((kapm+kapa-1),(kapn+kapy-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

              DO 111 k1=k1min,k1max
                call odd (lm+la+k1,i1)
                call odd (ln+ly+k1,i2)
                if (i1*i2.eq.0) goto 556
                cc=c*(2.d0*k+1.d0)*d6j(jm,jy,2*k,jn,ja,2*k1)
                index2=ip2(indn,indm,indy,i,k1)

                do 661 m=nnm,nmaxx(indm)
                  do 662 n=nnn,nmaxx(indn)
                    do 663  ny=nny,nmaxxv(indy)

******************************************************
                  if (key_en.eq.0) then
                   e0=ee(ny,indy)
                  endif
                  if (key_en.eq.1) then
                   e0=ee(nny,indy)
                  endif
                  if (key_en.eq.2) then
                   e0=ecorr(indy)
                  endif
*****************************************************

                  del=1.d0/(e0+ea-ee(m,indm)-ee(n,indn))
                    do 664  nx5=nnx,nmaxxv(indx)
      tt=cc*del*x2(m,n,nx5,index1)*x2(n,m,ny,index2)
                     sigma1(nx5,ny,indx)=sigma1(nx5,ny,indx)+tt
                     sigmad(nx5,ny,indx)=sigmad(nx5,ny,indx)-tt*del
                 en(nx5,ny,indy)=e0



664                 continue
663                 continue
662                 continue
661              continue

556            continue
111           continue
*************************************************************
555             continue
11            continue
25          continue
24         continue
2        continue
*****************************************************************
1       continue
c       write (*,*) '3s3s',sigma1(3,3,1),sigmad(3,3,1)
c       write (*,*) '3s4s',sigma1(4,3,1),sigmad(4,3,1)
c       write (*,*) '4s4s',sigma1(4,4,1),sigmad(4,4,1)
c       write (*,*) '3p3p',sigma1(2,2,2),sigmad(2,2,2)

901   continue
      return
      end

      subroutine sigma2code(max_orb,key4,key_en,nav,kvmax,key_allorder)
      implicit double precision (a-h,o-z)
c      integer*2 inv,inm,inw,inn,ink
************** NOTE : CHANGE LATER ***********************
      real*4 sigma2,der,energy
***********************************************************
      include "second.par"
c      dimension inv(NSIG),inm(NSIG),inw(NSIG),inn(NSIG),ink(NSIG)
      dimension sigma2(NSIG),der(NSIG),energy(NSIG)
      DIMENSION iint2(NSIG),iint3(NSIG)
      dimension IntOrb(ipx*ipx),tx(nxx),ty(nxx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      common /nmaxx/nmaxx(nx)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /radmnxa/ ip2(nk,nk,nk,nn,0:kk),
     * x2(nxx,nxx,nxx,0:nch2),ichan2
      COMMON /radmxab/ ip1(nk,nk,nn,nn,0:kk),x1(nxx,nxx,0:nch1),ichan1
      COMMON /eninput/ ecorr(nk1)
*     >>>>>>>>>>> Al-order >>>>>>>>>>>>>>>>
      COMMON /s2_all/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ich_all
      COMMON /val_all/ nv_all(nz),kv_all(nz),nval_all
      COMMON /i_all/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /energy_all/ evt_all(nz)
      COMMON /ds2_all/ dxvw(nxx,nxx,0:NVWH)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      iall=0
      do 1 i=1,ipx*ipx
      IntOrb(i)=0
1     continue
      index=0
      do 71 iv=ncore+1,max_orb
        kv=kor(iv)
        CALL klj(kv,kapv,lv,jv,indv,n0v)
      if (lv.gt.lvmax) goto 981
        nv=nor(iv)-n0v

      if (nv.gt.nxx1) then
      write (*,*) ' Parameter NXX1 is exceeded, STOP'
      stop
      endif


        ev=ee(nv,indv)
        do 72  iw=iv,max_orb
        kw=kor(iw)
          CALL klj(kw,kapw,lw,jw,indw,n0w)
        if (lw.gt.lvmax) goto 982
        nw=nor(iw)-n0w
          ew=ee(nw,indw)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>> ALL-order merge >>>>
*           check for v w
************************************************
      iiv=0
      iiw=0
      do 12 ii=1,nval_all
        if (nv_all(ii).eq.nor(iv).and.kv_all(ii).eq.kv) iiv=ii
12    continue

      do 13 ii=1,nval_all
        if (nv_all(ii).eq.nor(iw).and.kv_all(ii).eq.kw) iiw=ii
13    continue
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

********************************************
          iab=ipx*(iv-ncore-1)+(iw-ncore)
        IntOrb(iab)=index+1
********************************************
          DO 73  im=iv,max_orb
          km=kor(im)
            CALL klj(km,kapm,lm,jm,indm,n0m)
          if (lm.gt.lvmax) goto 983
          m=nor(im)-n0m
            em=ee(m,indm)
*************************************************
            if (im.gt.iv) then
             ihh=iv
            else
             ihh=iw
            endif
            do 74  in=ihh,max_orb
************************************************
            if (iv.eq.iw) then
               if (im.gt.in) then
              goto 800
             endif
            endif
************************************************
            if (iv.eq.in) then
               if (iw.gt.im) then
              goto 800
             endif
           endif
**************************************************
           kn=kor(in)
             CALL klj(kn,kapn,ln,jn,indn,n0n)
           if (ln.gt.lvmax) goto 984
           n=nor(in)-n0n
             en=ee(n,indn)
**************************************
           iii=iv+iw+in+im
             if (iii.gt.nav*4) goto 800
******************************************************
           if (key_en.eq.0) then
            evt=ev
            ewt=ew
           endif
           if (key_en.eq.1) then
              CALL st (kv,nnv)
              CALL st (kw,nnw)
            evt=ee(nnv,indv)
            ewt=ee(nnw,indw)
           endif
           if (key_en.eq.2) then
            evt=ecorr(indv)
            ewt=ecorr(indw)
           endif
*************************************
              kmin=max0(iabs(kapv-kapm),iabs(kapw-kapn))
              kmax=min0((kapv+kapm-1),(kapw+kapn-1),kvmax)
      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

              DO 121 k=kmin,kmax
                call odd (lv+lm,i1)
                call odd (lw+ln,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 500
                index=index+1
              if (index.gt.NSIG) then
               write (*,*) 'EXCEEDING MAX NUMBER OF SIGMA 2, STOP'
               stop
              endif
******************************************************************
c                inv(index)=iv
c                inm(index)=im
c                inw(index)=iw
c             inn(index)=in
c                ink(index)=k
***********SETTING UP ARRAYS FOR SCRC.CON ************************
              energy(index)=evt+ewt
              iint2(index)=ipx*ipx*k+ipx*(iv-ncore-1)+(im-ncore)
              iint3(index)=ipx*(iw-ncore-1)+(in-ncore)
*******************************************************************
*********************

                res=0.d0
              resd=0.d0
*******************************************************************
                res1=0.d0
                res2d=0.d0
                res2e=0.d0
                res5d=0.d0
                res5e=0.d0
                res3=0.d0
                res4=0.d0
                res5=0.d0
                res6=0.d0
                res7=0.d0
**********************************************************
*     Channel set up is complete
*     Now calculating sigma 2
*********************************************************
*********** START TERM 1 *********************************
                DO 810 ij=1,ncore
                kc=ko(ij)
                  ec=eo(ij)
                  CALL klj(kc,kapc,lc,jc,indc,n0c)
                  DO 820 jj=1,ncore
                  kd=ko(jj)
                    ed=eo(jj)
                    CALL klj(kd,kapd,ld,jd,indd,n0d)
                    cc1=(-1)**(kapm+kapn+kapc+kapd)
                  cc=cc1*(2.d0*k+1.d0)
                    k1min=max0(iabs(kapv-kapc),iabs(kapw-kapd))
                    k1max=min0((kapv+kapc-1),(kapw+kapd-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                    DO 830 l=k1min,k1max
                    call odd (lv+lc+l,i45)
                    call odd (lw+ld+l,i46)
                    if (i45*i46.eq.0) goto 840
                    index1=ip1(indv,indw,ij,jj,l)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 1 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     k2min=max0(iabs(kapm-kapc),iabs(kapn-kapd))
                     k2max=min0((kapm+kapc-1),(kapn+kapd-1))
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 850 k1=k2min,k2max
                     call odd (lm+lc+k1,i47)
                     call odd (ln+ld+k1,i48)
                     if (i47*i48.eq.0) goto 860
      zz=cc*d6j(jm,jv,2*k,2*l,2*k1,jc)*d6j(jn,jw,2*k,2*l,2*k1,jd)
                     if (zz.eq.0.d0) goto 860
                     index2=ip1(indm,indn,ij,jj,k1)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 1 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     del=1.d0/(ec+ed-em-en+evt-ev+ewt-ew)
                     tt=zz*del*x1(nv,nw,index1)*x1(m,n,index2)
                       res1=res1+tt
                     res=res+tt
                       resd=resd-tt*del
860                    CONTINUE
850                  CONTINUE
840                  CONTINUE
830                  CONTINUE
820             CONTINUE
810             CONTINUE
********** END TERM 1 ************************************
*********************************************************
* This part is common to terms 2-7
*********************************************************
                c2=(-1)**(kapw+kapn+k)/(2.d0*k+1.d0)
                c5=(-1)**(kapv+kapm+k)/(2.d0*k+1.d0)
                DO 200 indr=1,2*lmax+1
                  CALL indk1(indr,kr,kapr,lr,jr,n0r)
                  CALL st (kr,nnr)
                  DO 300 ij=1,ncore
                  kc=ko(ij)
                    ec=eo(ij)
                    CALL klj(kc,kapc,lc,jc,indc,n0c)
                  kk1=kapc+kapr-1
                  kk2=iabs(kapc-kapr)
                  if (k.le.kk1.and.k.ge.kk2) then
*****************This part is common for TERMS 2, 3, 5, and 6 ****
                   call odd (lr+lc+k,i3)
                     if (i3.eq.0) goto 110
*****************************************************************
*                Start of TERMs 2 and 5
*****************************************************************
*        Direct part so terms 2 and 5
******************************************************************
                   call odd (lv+lm+k,i4)
                     if (i4.eq.0) goto 120
                   call odd (lw+ln+k,i7)
                     if (i7.eq.0) goto 140
                   index1=ip2(indw,indr,indn,ij,k)
                   index2=ip2(indm,indr,indv,ij,k)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 2d for SIGMA 2'
                     stop
                    endif
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 2d for SIGMA 2'
                     stop
                    endif
**********************************************************************
                   index3=ip2(indv,indr,indm,ij,k)
                   index4=ip2(indn,indr,indw,ij,k)
*********************************************************************
                    if (index3.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 5d for SIGMA 2'
                     stop
                    endif
                    if (index4.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 5d for SIGMA 2'
                     stop
                    endif
**********************************************************************
                   DO 64  nr=nnr,nmaxx(indr)
                     del=1.d0/(evt+ewt-ew+ec-em-ee(nr,indr))
                     del1=1.d0/(ewt+evt-ev+ec-en-ee(nr,indr))
              tt=c2*del*x2(nw,nr,n,index1)*x2(m,nr,nv,index2)
                tt1=c5*del1*x2(nv,nr,m,index3)*x2(n,nr,nw,index4)
                      res2d=res2d+tt
                      res5d=res5d+tt1
                    res=res+tt+tt1
                      resd=resd-tt*del-tt1*del1
64                   CONTINUE
140                  CONTINUE
***************Exchange part of Term 2*********************
                     k1min=max0(iabs(kapr-kapn),iabs(kapw-kapc))
                     k1max=min0((kapr+kapn-1),(kapw+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 129 k1=k1min,k1max
                     call odd (lr+ln+k1,i5)
                     call odd (lw+lc+k1,i6)
                     if (i5*i6.eq.0) goto 130
                       zz=c2*(2.d0*k+1.d0)*d6j(jw,jn,2*k,jr,jc,2*k1)
                     index1=ip2(indr,indw,indn,ij,k1)
                     index2=ip2(indm,indr,indv,ij,k)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 2e for SIGMA 2'
                     stop
                    endif
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 2e for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 65  nr=nnr,nmaxx(indr)
                       del=1.d0/(evt+ewt-ew+ec-em-ee(nr,indr))
                 tt=zz*del*x2(nr,nw,n,index1)*x2(m,nr,nv,index2)
                         res2e=res2e+tt
                       res=res+tt
                         resd=resd-tt*del
65                     CONTINUE
130                    CONTINUE
129                  CONTINUE
120                  CONTINUE
**************** END of Term 2 *************************
***************Exchange part of Term 5*********************

                   call odd (lw+ln+k,i17)
                     if (i17.eq.0) goto 220
                   index2=ip2(indn,indr,indw,ij,k)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 5e for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     k1min=max0(iabs(kapr-kapm),iabs(kapv-kapc))
                     k1max=min0((kapr+kapm-1),(kapv+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 329 k1=k1min,k1max
                     call odd (lr+lm+k1,i18)
                     call odd (lv+lc+k1,i19)
                     if (i18*i19.eq.0) goto 230
                       zz=(2.d0*k+1.d0)*d6j(jv,jm,2*k,jr,jc,2*k1)
                     index1=ip2(indr,indv,indm,ij,k1)
*********************************************************************
                     if (index1.eq.0) then
                      write (*,*) 'index1 = 0 in TERM 5e for SIGMA 2'
                      stop
                     endif
**********************************************************************
                       DO 165  nr=nnr,nmaxx(indr)
                       del=1.d0/(ewt+evt-ev+ec-en-ee(nr,indr))
                 tt=zz*del*c5*x2(nr,nv,m,index1)*x2(n,nr,nw,index2)
                         res5e=res5e+tt
                       res=res+tt
                         resd=resd-tt*del
165                    CONTINUE
230                    CONTINUE
329                  CONTINUE
220                  continue
******************** End of exchange part of term 5 **********

*****************START OF TERM 3 *****************************
                   call odd (lw+ln+k,i8)
                     if (i8.eq.0) goto 150
                   index1=ip2(indw,indr,indn,ij,k)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 3 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     k1min=max0(iabs(kapr-kapv),iabs(kapm-kapc))
                     k1max=min0((kapr+kapv-1),(kapm+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 229 k1=k1min,k1max
                     call odd (lr+lv+k1,i9)
                     call odd (lm+lc+k1,i10)
                     if (i9*i10.eq.0) goto 160
                       zz=c2*(2.d0*k+1.d0)*d6j(jm,jv,2*k,jr,jc,2*k1)
                     index2=ip2(indr,indm,indv,ij,k1)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 3 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 66  nr=nnr,nmaxx(indr)
                       del=1.d0/(evt+ewt-ew+ec-em-ee(nr,indr))
            tt=zz*del*x2(nw,nr,n,index1)*x2(nr,m,nv,index2)
                        res3=res3+tt
                      res=res+tt
                        resd=resd-tt*del
66                     CONTINUE
160                    CONTINUE
229                  CONTINUE
150                  CONTINUE
*****************END OF TERM 3 *******************************
*****************START OF TERM 6 *****************************
                   call odd (lv+lm+k,i20)
                     if (i20.eq.0) goto 250
                   index1=ip2(indv,indr,indm,ij,k)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 6 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     k1min=max0(iabs(kapr-kapw),iabs(kapn-kapc))
                     k1max=min0((kapr+kapw-1),(kapn+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 279 k1=k1min,k1max
                     call odd (lr+lw+k1,i21)
                     call odd (ln+lc+k1,i22)
                     if (i21*i22.eq.0) goto 260
                       zz=c5*(2.d0*k+1.d0)*d6j(jn,jw,2*k,jr,jc,2*k1)
                     index2=ip2(indr,indn,indw,ij,k1)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 6 for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 661  nr=nnr,nmaxx(indr)
                       del=1.d0/(ewt+evt-ev+ec-en-ee(nr,indr))
            tt=zz*del*x2(nv,nr,m,index1)*x2(nr,n,nw,index2)
                        res6=res6+tt
                      res=res+tt
                        resd=resd-tt*del
661                     CONTINUE
260                    CONTINUE
279                  CONTINUE
250                  CONTINUE
*****************END OF TERM 6 *******************************

110                  CONTINUE
*********** Common end for terms 2, 3, 5, 6,
                    endif
*****************START OF TERM 4b MK notes version  *****************************
c          goto 9000
                     if (key4.eq.0) goto 600
                     kkmin=max0(iabs(kapv-kapn),iabs(kapr-kapc),
     *                          iabs(kapm-kapw))
                  kkmax=min0((kapv+kapn-1),(kapr+kapc-1),(kapm+kapw-1))
      if (kkmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 999 l=kkmin,kkmax
                       cc0=(2.d0*k+1.d0)*d6j(jm,jv,2*k,jn,jw,2*l)
                       c4b=cc0*(-1)**(kapv+kapn+l)*(2.d0*l+1.d0)
                       k1min=max0(iabs(kapr-kapn),iabs(kapv-kapc))
                       k1max=min0((kapr+kapn-1),(kapv+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif
                       DO 239 k1=k1min,k1max
                       call odd (lr+ln+k1,i11)
                       call odd (lv+lc+k1,i12)
                       if (i11*i12.eq.0) goto 170
                       zz2=c4b*d6j(jv,jn,2*l,jr,jc,2*k1)
                       k2min=max0(iabs(kapr-kapw),iabs(kapm-kapc))
                       k2max=min0((kapr+kapw-1),(kapm+kapc-1))
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     index1=ip2(indr,indv,indn,ij,k1)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 4b for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 249 k2=k2min,k2max
                       call odd (lr+lw+k2,i13)
                       call odd (lm+lc+k2,i14)
                       if (i13*i14.eq.0) goto 180
                         zz=zz2*d6j(jm,jw,2*l,jr,jc,2*k2)
                       index2=ip2(indr,indm,indw,ij,k2)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 4b for SIGMA 2'
                     stop
                    endif
**********************************************************************
                         DO 67  nr=nnr,nmaxx(indr)
                         del=1.d0/(ewt+evt-ev+ec-em-ee(nr,indr))
                      tt=zz*del*x2(nr,nv,n,index1)*x2(nr,m,nw,index2)
                           res4=res4+tt
                         res=res+tt
                           resd=resd-tt*del
67                       CONTINUE
180                      CONTINUE
249                    CONTINUE
170                    CONTINUE
239                  CONTINUE
999                  CONTINUE
600                  CONTINUE
9000            continue

*****************END OF TERM 4b *******************************
*****************START OF ALTERNATIVE VERSION TERM 4b MK notes version  **
******* COMMENTED OUT - TURNED OUT TO BE SLOWER ***************************
       goto 807
                     if (key4.eq.0) goto 620
                     kkmin=max0(iabs(kapv-kapn),iabs(kapr-kapc),
     *                          iabs(kapm-kapw))
                  kkmax=min0((kapv+kapn-1),(kapr+kapc-1),(kapm+kapw-1))
      if (kkmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 991 l=kkmin,kkmax
                       cc0=(2.d0*k+1.d0)*d6j(jm,jv,2*k,jn,jw,2*l)
                       c4b=cc0*(-1)**(kapv+kapn+l)/(2.d0*l+1.d0)
                       k1min=max0(iabs(kapr-kapn),iabs(kapv-kapc))
                       k1max=min0((kapr+kapn-1),(kapv+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                       dd=c4b*(2.d0*l+1.d0)*(2.d0*l+1.d0)
******************************************************************
                       DO 672  nr=nnr,nmaxx(indr)
                       tx(nr)=0.d0
                       ty(nr)=0.d0
672                    CONTINUE
******************************************************************
                       DO 2391 k1=k1min,k1max
                       call odd (lr+ln+k1,i11)
                       call odd (lv+lc+k1,i12)
                       if (i11*i12.eq.0) goto 1701
                         zz2=d6j(jv,jn,2*l,jr,jc,2*k1)
                       index1=ip2(indr,indv,indn,ij,k1)
                       if (index1.eq.0) then
                       write (*,*) 'index1=0 in TERM 4b for SIGMA 2'
                        stop
                       endif
                         DO 671  nr=nnr,nmaxx(indr)
                         tx(nr)=tx(nr)+zz2*x2(nr,nv,n,index1)
671                      CONTINUE
1701                      continue
2391                    continue
*********************************************************************

                       k2min=max0(iabs(kapr-kapw),iabs(kapm-kapc))
                       k2max=min0((kapr+kapw-1),(kapm+kapc-1))
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                       DO 2491 k2=k2min,k2max
                       call odd (lr+lw+k2,i13)
                       call odd (lm+lc+k2,i14)
                       if (i13*i14.eq.0) goto 1801
                         zz=d6j(jm,jw,2*l,jr,jc,2*k2)
                       index2=ip2(indr,indm,indw,ij,k2)
                       if (index2.eq.0) then
                        write (*,*) 'index2=0 in TERM 4b for SIGMA 2'
                        stop
                       endif
                         DO 676  nr=nnr,nmaxx(indr)
                        ty(nr)=ty(nr)+zz*x2(nr,m,nw,index2)
676                       CONTINUE
1801                      CONTINUE
2491                    CONTINUE
**********************************************************************
                         DO 677  nr=nnr,nmaxx(indr)
                         del=1.d0/(ewt+evt-ev+ec-em-ee(nr,indr))
                         tt=dd*del*tx(nr)*ty(nr)
                           res4=res4+tt
                         res=res+tt
                           resd=resd-tt*del
677                       CONTINUE
*************************************************************************
991                  CONTINUE
620                  CONTINUE

807       continue
*****************END OF TERM 4b *******************************
*****************START OF TERM 4a original all-order code version  ************

                     if (key4.eq.1) goto 601

                     k1min=max0(iabs(kapr-kapn),iabs(kapw-kapc))
                     k1max=min0((kapr+kapn-1),(kapw+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif
                     DO 259 k1=k1min,k1max
                     call odd (lr+ln+k1,i15)
                     call odd (lw+lc+k1,i16)
                     if (i15*i16.eq.0) goto 190
           zz2=c2*(2.d0*k+1.d0)*(2.d0*k+1.d0)*d6j(jw,jn,2*k,jr,jc,2*k1)
                   index1=ip2(indr,indw,indn,ij,k1)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 4a for SIGMA 2'
                     stop
                    endif
**********************************************************************
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif
                     k2min=max0(iabs(kapr-kapv),iabs(kapm-kapc))
                     k2max=min0((kapr+kapv-1),(kapm+kapc-1))
                     DO 269 k2=k2min,k2max
                     call odd (lr+lv+k2,i17)
                     call odd (lm+lc+k2,i18)
                     if (i17*i18.eq.0) goto 210
                       zz=zz2*d6j(jm,jv,2*k,jr,jc,2*k2)
                     index2=ip2(indr,indm,indv,ij,k2)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 4a for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 68  nr=nnr,nmaxx(indr)
                       del=1.d0/(evt+ewt-ew+ec-em-ee(nr,indr))
             tt=zz*del*x2(nr,nw,n,index1)*x2(nr,m,nv,index2)
                        res=res+tt
                        resd=resd-tt*del
68                     CONTINUE
210                    CONTINUE
269                  CONTINUE
190                  CONTINUE
259                  CONTINUE
601                  CONTINUE
******************************************************************

*****************START OF TERM 7b MK notes version  *****************************
                     if (key4.eq.0) goto 602
                     kkmin=max0(iabs(kapv-kapn),iabs(kapr-kapc),
     *                          iabs(kapm-kapw))
             kkmax=min0((kapv+kapn-1),(kapr+kapc-1),(kapm+kapw-1))
      if (kkmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 998 l=kkmin,kkmax
                      cc0=(2.d0*k+1.d0)*d6j(jm,jv,2*k,jn,jw,2*l)
                      c7b=cc0*(-1)**(kapw+kapm+l)/(2.d0*l+1.d0)
                       k1min=max0(iabs(kapr-kapm),iabs(kapw-kapc))
                       k1max=min0((kapr+kapm-1),(kapw+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                       DO 339 k1=k1min,k1max
                       call odd (lr+lm+k1,i24)
                       call odd (lw+lc+k1,i25)
                       if (i24*i25.eq.0) goto 370
        zz2=c7b*(2.d0*l+1.d0)*(2.d0*l+1.d0)*d6j(jw,jm,2*l,jr,jc,2*k1)
                         k2min=max0(iabs(kapr-kapv),iabs(kapn-kapc))
                         k2max=min0((kapr+kapv-1),(kapn+kapc-1))
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                       index1=ip2(indr,indw,indm,ij,k1)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 7b for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 349 k2=k2min,k2max
                       call odd (lr+lv+k2,i26)
                       call odd (ln+lc+k2,i27)
                       if (i26*i27.eq.0) goto 380
                         zz=zz2*d6j(jn,jv,2*l,jr,jc,2*k2)
                       index2=ip2(indr,indn,indv,ij,k2)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 7b for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 367  nr=nnr,nmaxx(indr)
                       del=1.d0/(evt+ewt-ew+ec-en-ee(nr,indr))
               tt=zz*del*x2(nr,nw,m,index1)*x2(nr,n,nv,index2)
                       res7=res7+tt
                      res=res+tt
                        resd=resd-tt*del
367                     CONTINUE
380                    CONTINUE
349                  CONTINUE
370                  CONTINUE
339                  CONTINUE
998                  CONTINUE
602                  CONTINUE
*****************END OF TERM 7b *******************************
*****************START OF TERM 7a original all-order code version  ************

                     if (key4.eq.1) goto 603
                     k1min=max0(iabs(kapr-kapm),iabs(kapv-kapc))
                     k1max=min0((kapr+kapm-1),(kapv+kapc-1))
      if (k1max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 459 k1=k1min,k1max
                     call odd (lr+lm+k1,i28)
                     call odd (lv+lc+k1,i29)
                     if (i28*i29.eq.0) goto 490
           zz2=c5*(2.d0*k+1.d0)*(2.d0*k+1.d0)*d6j(jv,jm,2*k,jr,jc,2*k1)
                   index1=ip2(indr,indv,indm,ij,k1)
*********************************************************************
                    if (index1.eq.0) then
                     write (*,*) 'index1 = 0 in TERM 7a for SIGMA 2'
                     stop
                    endif
**********************************************************************
                     k2min=max0(iabs(kapr-kapw),iabs(kapn-kapc))
                     k2max=min0((kapr+kapw-1),(kapn+kapc-1))
      if (k2max.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

                     DO 469 k2=k2min,k2max
                     call odd (lr+lw+k2,i30)
                     call odd (ln+lc+k2,i31)
                     if (i30*i31.eq.0) goto 410
                       zz=zz2*d6j(jn,jw,2*k,jr,jc,2*k2)
                     index2=ip2(indr,indn,indw,ij,k2)
*********************************************************************
                    if (index2.eq.0) then
                     write (*,*) 'index2 = 0 in TERM 7a for SIGMA 2'
                     stop
                    endif
**********************************************************************
                       DO 468  nr=nnr,nmaxx(indr)
                       del=1.d0/(ewt+evt-ev+ec-en-ee(nr,indr))
             tt=zz*del*x2(nr,nv,m,index1)*x2(nr,n,nw,index2)
                        res=res+tt
                        resd=resd-tt*del
468                     CONTINUE
410                    CONTINUE
469                  CONTINUE
490                  CONTINUE
459                  CONTINUE
603                  CONTINUE
******************************************************************
300               CONTINUE
200             CONTINUE
***************** Common end for terms 2-7 ***************
                ck=1.d0/(((-1)**(k))*s1(k,km,kv)*s1(k,kn,kw))
                ress=res1+res2d+res2e+res3+res4+res5d+res5e+res6+res7
              sy1=res*ck
                sy2=resd*ck
                sigma2(index)=sy1
              der(index)=sy2
       goto 877
      write (*,301) res1*ck
      write (*,302) res2d*ck
      write (*,303) res2e*ck
      write (*,304) res3*ck
      write (*,305) res4*ck
      write (*,306) res5d*ck
      write (*,307) res5e*ck
      write (*,308) res6*ck
      write (*,309) res7*ck

      write (*,'(i5,9i4,2i7,3e15.6)') index,nv+n0v,kv,nw+n0w,
     *kw,m+n0m,km,n+n0n,kn,k,iint2(index),iint3(index),
     *res*ck,resd*ck,energy(index)
877   continue

301   format ('Term 1  = ',e15.8)
302   format ('Term 2d = ',e15.8)
303   format ('Term 2e = ',e15.8)
304   format ('Term 3  = ',e15.8)
305   format ('Term 4  = ',e15.8)
306   format ('Term 5d = ',e15.8)
307   format ('Term 5e = ',e15.8)
308   format ('Term 6  = ',e15.8)
309   format ('Term 7  = ',e15.8)
* >>>>>>>>> FIND ALL_ORDER
          if (iiv.ne.0.and.iiw.ne.0) then
          iall=iall+1
          ind_all=ipvw(iiv,iiw,indm,indn,k)
          r=xvw(m,n,ind_all)
          dr=dxvw(m,n,ind_all)
          de=evt_all(iiv)+evt_all(iiw)
c      write (*,'(i5,9i4,3e15.6)') index,nv+n0v,kv,nw+n0w,
c     *kw,m+n0m,km,n+n0n,kn,k,res*ck,resd*ck,energy(index)
c      write (*,'(i5,9i4,3e15.6)') index,nv+n0v,kv,nw+n0w,
c     *kw,m+n0m,km,n+n0n,kn,k,r,dr,de
              sigma2(index)=r
              der(index)=dr
                energy(index)=de
c        write (*,*)
        endif
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
500             CONTINUE
121           CONTINUE
800           CONTINUE
984         CONTINUE
74          CONTINUE
983       CONTINUE
73        CONTINUE
982       CONTINUE
72      CONTINUE
981     CONTINUE
71    CONTINUE
       write (*,*) 'Total Sigma2 count                     = ', index
       write (*,*) 'Number of Sigma2 replaced by all order = ', iall

**************** WRITEOUT**********************************
      open(unit=20,file='SCRC.CON',status='UNKNOWN',
     *form='unformatted')
      Nso=ncore
      Nsh=ncore
      Khot=0
      nsx1=ncore
      nsx2=0
      Ksym=1
      Kbox=1
      igint=index
      write (20) norb,Nso,Nsh,Khot,nsx1,nsx2,Ksym
      write (20) max_orb,lvmax,kvmax,Nav,Kbox
      write (20) igint,ipx*ipx,igint
      write (20) (sigma2(i),i=1,igint)
      write (20) (iint2(i),i=1,igint)
      write (20) (iint3(i),i=1,igint)
      write (20) (IntOrb(i),i=1,ipx*ipx)
      write (20) (der(i),i=1,igint)
      write (20) (energy(i),i=1,igint)
c      write (20) (inv(i),inw(i),inm(i),inn(i),ink(i),i=1,igint)

      close (20)
*********************************************************
      return
      end

      subroutine count_sigma2code(max_orb,key4,nav,kvmax)
      implicit double precision (a-h,o-z)
      include "second.par"
      common /nmaxx/nmaxx(nx)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxxv/ nmaxxv(nx),lvmax
      common /bassinp/ nor(nsorb),kor(nsorb),norb
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      index=0
      do 71 iv=ncore+1,max_orb
        kv=kor(iv)
        CALL klj(kv,kapv,lv,jv,indv,n0v)
      if (lv.gt.lvmax) goto 981
        nv=nor(iv)-n0v
        do 72  iw=iv,max_orb
        kw=kor(iw)
          CALL klj(kw,kapw,lw,jw,indw,n0w)
        if (lw.gt.lvmax) goto 982
        nw=nor(iw)-n0w
********************************************
          DO 73  im=iv,max_orb
          km=kor(im)
            CALL klj(km,kapm,lm,jm,indm,n0m)
          if (lm.gt.lvmax) goto 983
          m=nor(im)-n0m
*************************************************
            if (im.gt.iv) then
             ihh=iv
            else
             ihh=iw
            endif
            do 74  in=ihh,max_orb
************************************************
            if (iv.eq.iw) then
               if (im.gt.in) then
              goto 800
             endif
            endif
************************************************
            if (iv.eq.in) then
               if (iw.gt.im) then
              goto 800
             endif
           endif
**************************************************
           kn=kor(in)
             CALL klj(kn,kapn,ln,jn,indn,n0n)
           if (ln.gt.lvmax) goto 984
           n=nor(in)-n0n
**************************************
           iii=iv+iw+in+im
             if (iii.gt.nav*4) goto 800
*************************************
              kmin=max0(iabs(kapv-kapm),iabs(kapw-kapn))
              kmax=min0((kapv+kapm-1),(kapw+kapn-1),kvmax)
      if (kmax.gt.kk) then
      write (*,*) ' Parameter KK is exceeded, STOP'
      stop
      endif

              DO 121 k=kmin,kmax
                call odd (lv+lm,i1)
                call odd (lw+ln,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 500
                index=index+1

*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
500             CONTINUE
121           CONTINUE
800           CONTINUE
984         CONTINUE
74          CONTINUE
983       CONTINUE
73        CONTINUE
982       CONTINUE
72      CONTINUE
981     CONTINUE
71    CONTINUE
      write (*,*)
      if (index.le.NSIG) then
      write (*,122) index,nsig
      else
      write (*,124) index,nsig
      write (*,*) 'Parameter NSIG is exceeded, STOP'
      stop
      endif
122   format ('Number of SIGMA2_k(vwmn) channels = ',i9,
     *' < NSIG = ',i9,' OK')

124   format ('Number of SIGMA2_k(vwmn) channels = ',i9,
     *' > NSIG = ',i9)


      return
      end


      function s1(k,k1,k2)
      implicit doubleprecision(a-h,o-z)
      call klj(k1,kap1,l1,j1,ind,n)
      call klj(k2,kap2,l2,j2,ind,n)
      a=l1+l2+k
      b=a/2-INT(a/2)
c      if (b.ne.0.0) then
c      c1=0.d0
c      else
      c1=1.d0
c      endif
      coef=c3j(j1,j2,2*k,-1,1,0)
      s1=2.d0*c1*coef*SQRT(kap1*kap2*1.d0)*((-1)**kap1)
      end


c     ========================================
      include "d6j.f"
      include "libD.f"
      include "rint.f"
      include "yfun.f"
      include "yint.f"
      include "inidat.f"
      include "rec_unit.inc"
