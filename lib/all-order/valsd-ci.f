* Version 4.0 April 20, 2012
* Stabilizer is inserted, insertions are marked with RLE
* Derivatives are not iterative and nothing needs to be done
* When sigmas are fixed, so will be derivatives
**************************************
* Version 3.1 October 17, 2008
*  Number of nd state iterations is given separately
**************************************
*     NEW utterly modified version 
*     Version 3.0 May 8, 2008
*     rho equations are replaced by sigma's equations,
*     derivatives are added, options with kval = 1 and 2 are added 
*************************************************************************     
*     Version 2.0 May 6, 2008
*     Parameters are taken out, stops are put in if parameters are 
*     exceeded, counting subroutines are added.
**************************************************************************
* In this version the term31 is excluded and energy correction 
* deltae is put to zero. That variant corresponds to SD+CI approach
**************************************************************************
*************************************************************************
*                                                                       *
*   Explanation of input  (no empty lines allowed)                      *  
*                                                                       *
*************************************************************************
* This input can be used by any of the all-order package codes          *
*                                                                       * 
* allcore-cis ignores all below max number of iterations                *
* valsd-cis   needs the entire input                                    *
* sdvw-cis    ignores kval input but needed list of valence             *
*             electrons below it                                        *
*************************************************************************
*                                                                       *
* 4	   # number of core shells ncore                                  *
* 1 -1   # n and kappa of the core shells                               *
* 2 -1                                                                  *
* 2  1                                                                  *
* 2 -2                                                                  *
* 25 5   # nmax and lmax in correlation diagram summations              *
* 0  0   # internal paremeters, do not change for valsd-cis and sdvw-cis*
*        # 1  0 option can be used for allcore-cis                      *
*               to use previously iterated file                         *
* 10  5  # max number of iterations, max number of iterations for nd states*
* 1  	   # kval (key_en): key for energies,                             *
*        # DHF if kval=0, lowest for a given kappa for kval=1           *
*        # if kval=2, the energies are inputted in the format below     *
************************************************************************* 
* Example of kval=2 input (ignored by sdvw-cis):                        *
*                                                                       *
* 2                      # kval                                         *
* 3                      # lvmax for the input to follow                *
* 0  -0.28000            # l=0 energies                                 * 
* 1  -0.22000  -0.22000  # l=1 energies p1/2 p3/2                       *
* 2  -0.31000  -0.31000                                                 *
* 3  -0.13000  -0.13000                                                 *
*                                                                       *
*************************************************************************
* 7       # nval number of valence orbitals in list to follow           *
* 3  -1   # n and kappa of the valence orbitals                         *
* 4  -1                                                                 *
* 5  -1                                                                 *
* 3   1                                                                 *
* 4   1                                                                 *
* 3   2                                                                 *
* 4   2                                                                 *
*************************************************************************
      implicit real*8 (a-h,o-z)
	include "all.par"
      DIMENSION t(nkx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      DIMENSION g(nhf,nx),f(nhf,nx),e(nl)
      DIMENSION wco(ns),nco(ns),kco(ns),mco(ns)
      COMMON /x12ini/ inmax
*>>>>>>> RLE >>>>>>>>>>>>>>>>>>>>
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nxx)
      COMMON /all/ xvo_all(nxx,nxx,nvh,nit),rvo_all(nxx,nit)
      COMMON /energysd/ delsd(nit)
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
* >>>>>>>>> SIGMA >>>>>
      DIMENSION ecorr(nk)
* >>>>>>>>>>>>>>>>>>>>>

      open (17,file='sigma')
      OPEN (unit=9,file='val2',form='unformatted')
c      open (15,file='inf.aov',form='FORMATTED',status='OLD',err=710)
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
      read (5,*) key1,key2

*>>>>>>> RLE >>>>>>>>>>>>>>>>>>>>
      read (5,*) nitmax999
      read (5,*) ikey,max_rle,i999
      read (5,*) damp
* i999 is used in code only
* nitmax is not used 
* It is added to use the same input file as for the core 
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

c      read (5,*) nits,nitd
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
      write (*,*) nmax,lmax,nmax1,lmax1
      do i=1,2*lmax+1
        if (nmaxx(i).GT.nmax) nmaxx(i)=nmax
        write(*,*) ' PW=',i,' nmax=',nmaxx(i)
      end do
* >>>>>>>> SIGMA >>>>>>>>>
      read (5,*) kval
     
      if (kval.eq.0) then 
	 write (*,*) 'The HF energy is used, kval = 0'
      endif
	if (kval.eq.1) then 
	 write (*,*) 'The energy is set to lowest n energy for each ind'
      endif
	if (kval.eq.2) then
	 read (5,*) lvmax 
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
      if (kval.lt.0.or.kval.gt.2) then 
	 write (*,*) 'WRONG INPUT IN KVAL, STOP'
	endif
* >>>>>>>>>>>>>>>>>>>>>>>>
****************************************************************
      read (5,*) nbval
      icount_val2=0

	do 1999 iil=1,nbval
        READ (5,*) nv,kv,nitmax
        CALL klj (kv,kapv,lv,jv,indv,n0v)
        DO 29 j=1,max
          fv(j)=ff(j,nv-n0v,indv)
          gv(j)=gg(j,nv-n0v,indv)
29      CONTINUE
        ev=ee(nv-n0v,indv)
        call st (kv,nnv)

*>>>> SIGMA >>>>>>>>>>>  
      if (kval.eq.0) then 
       evt=ev
	endif

      if (kval.eq.1) then 
       evt=ee(nnv,indv)
	endif

      if (kval.eq.2) then 
	 if (lv.gt.lvmax) then
      write (*,*) 'Lvmax is exceeded, no energy input for kval=2, STOP'
       stop
	endif
	 evt=ecorr(indv)
	endif


*>>>> SIGMA >>>>>>>>>>>        	   

       if (indv.gt.inmax) then
       inmax=indv
       endif

      call xcore
	call countx
	call countrho

      call xval
      call rhoval(key1,iv)

      call rhocore
      if (iv.eq.0) goto 700

17    FORMAT ('Time required = ',f9.3,'sec')

        call energy1(res1)
        call energy2(res2)
        call energy3(res3)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>
        call energy4(evt,res4)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        deltae=res1+res2+res3+res4
c        write (*,'(5f20.10)') res1,res2,res3,res4,deltae
        write (*,880) ev
        write (*,671) evt

      write (*,67) deltae
	 write (*,*) 
67    format  ('Starting correlation energy = ',f12.8) 

880    format ('DHF energy                  = ',f12.8) 
671    format ('Ev tilde energy             = ',f12.8) 
      
c	  if (lv.eq.2.or.lv.eq.3) then 
c	   nitmax=nitd
c	  else
c	   nitmax=nits
c	  endif
cRLE******************************************************************
       irle=0
       
cRLE*************************************************************

      do 1000 ik=1,nitmax

        t1=mclock()

        call outin
* >>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>
        call terma1(evt)
        call terma2(evt)
        call terma3(evt)
        call terma4(evt)

        call term1(evt)
        call term2(evt)
c        call term31
        call term32(evt)
        call term41(evt)
        call term42(evt)
        call term51(evt)
        call term52(evt)
        t2=mclock()
        td=(t2-t1)/1000
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c700     CONTINUE
        deltae0=0.d0
        call outt(deltae0)
        call energy1(res1)
        call energy2(res2)
        call energy3(res3)
* >>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>
        call energy4(evt,res4)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        deltae1=res1+res2+res3+res4
c        write (*,'(i5,5f20.10)') ik,res1,res2,res3,res4,deltae1
        or=abs((deltae-deltae1)/deltae)
        write (*,907) nv,kv,ik,deltae,deltae1,or,td
907     format (2i3,' Iter',i3,3f15.8,' Time = ',f7.1,' sec')


        if (or.lt.0.00001) goto 3000
c         if (or.lt.0.0000001) goto 3000


        deltae=deltae1
c >>>>>>>>> START RLE **********Accumulate m iterations *****
        if (ikey.ne.0) then 
        irle=irle+1
        delsd(irle)=deltae1

        do 27 i=1,ivchan
          do 37 m=1,nxx      
            do 748 n=1,nxx
              xvo_all(m,n,i,irle)=xvo(m,n,i)
748          continue
37        continue
27      continue

      do 47 i=1,ncore
        do 718 n=1,nxx
        rvo_all(n,irle)=rvo(n)
718     continue
47     continue

c      call rho_read
       if (irle.eq.max_rle) then 
c        call energy_rle

        if (ikey.eq.1) then 
         call rle(evt,ene,deltae,max_rle)
        endif

        if (ikey.eq.2) then 
         call diis(evt,ene,deltae,max_rle)
        endif

        deltae=ene
        deltae1=ene

       irle=0
      endif
      endif

*>>>>>>>>>>>>END RLE >>>>>>>>>>>>

1000  continue
3000  continue
c      write (*,'(f20.10)')  deltae1
700   continue
     
      call outrho(deltae1,icount_val2,evt,kval)
	call outsigma
1999  continue
      write (*,*) 
	write (*,*) 'IMPORTANT! TOTAL NUMBER OF VALENCE CHANNELS'
	write (*,*) 'SET NVH2 PARAMETER in SDVW-CI TO EXCEED THIS VALUE'
       write (*,*) icount_val2
122   format ('Number of rho_k(mnav) channels      = ',i9)
       write (*,*)

      close (7)
	close (17)
890   format ('RES',i5,3f20.10)
      stop
 710  write(*,*) ' No input file "inf.aov" found '
      read(*,*)
      stop
 720  write(*,*) ' No file "hfspl.1" ifound '
      stop
 730  write(*,*) ' No file "hfspl.2" found '
      end

      subroutine rle(evt,ene,deltae,max_rle)
      implicit double precision (a-h,o-z)
      integer*4 ipvt
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nxx)
      COMMON /all/ xvo_all(nxx,nxx,nvh,nit),rvo_all(nxx,nit)

      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)    
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      common /tau/ Rm(nit,nit),al(nit),bl(nit)
      dimension z(nit),ipvt(nit)
      COMMON /energysd/ delsd(nit)

c      del=deltae

      CALL klj(kv,kapv,lv,jv,indv,n0v)
      do 112 i=1,nit
        al(i)=0.d0
        do 113 j=1,nit
          Rm(i,j)=0.d0
113       continue
112    continue

      write (*,*) 'Print energies' 
      write (*,*) 
      do 117 i=1,nit
        write (*,'(i5,f13.8)') i,delsd(i) 
117    continue
      write (*,*) 'DONE'
	write (*,*) 

c ********* Accumulate  al

      nt=max_rle-1
      do 1600 i=1,nt


        do 12 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          do 13 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 14 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
   
              kmin=max0(iabs(kapm-kapv),iabs(kapn-kapb))
              kmax=min0((kapm+kapv-1),(kapn+kapb-1))
              DO 15 lk=kmin,kmax

                call odd (lm+lv+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 16
               
            index1=ipxv1(i2,indm,indn,lk)
            index2=ipxv(i2,indm,indn,lk)                 

           zz=1.d0/((2.d0*lk+1.d0)*(jv+1.d0))
           do 17 m=nnm,nmaxx(indm)
             do 18 n=nnn,nmaxx(indn)               
                 al(i)=al(i)-zz*xvo_all(m,n,index2,i)*xv1(m,n,index1)  
18            continue 
17          continue
16       continue
15            continue 
14          continue
13       continue
12          continue
1600  continue

c ********* Accumulate Rm 

  
      do 1000 i=1,nt
       do 2000 j=1,nt


         indn=indv
         CALL indk1(indn,kn,kapn,ln,jn,n0)
         call st (kn,nnn) 
         do 61 n=nnn,nmaxx(indn)	
	     if (n+n0.ne.nv) then       
c           den1=ev+delsd(j)-ee(n,indn)
c           den2=ev+delsd(j+1)-ee(n,indn)
            den1=1.d0
            den2=1.d0

           Rm(i,j)=Rm(i,j)+rvo_all(n,i)*den1*rvo_all(n,j)-
     *                     rvo_all(n,i)*den2*rvo_all(n,j+1)
	    endif
61     continue

        do 22 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          do 23 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 24 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************    
              kmin=max0(iabs(kapm-kapv),iabs(kapn-kapb))
              kmax=min0((kapm+kapv-1),(kapn+kapb-1))
              DO 25 lk=kmin,kmax

                call odd (lm+lv,ii1)
                call odd (ln+lb,ii2)
                if ((ii1.eq.0.and.ii2.ne.0).or.
     *              (ii1.ne.0.and.ii2.eq.0)) goto 26

            index2=ipxv(i2,indm,indn,lk)                 
    
           zz=1.d0/((2.d0*lk+1.d0)*(jv+1.d0))
           do 27 m=nnm,nmaxx(indm)
             do 28 n=nnn,nmaxx(indn)
c               den1=ev+eb+delsd(j)-ee(m,indm)-ee(n,indn) 
c               den2=ev+eb+delsd(j+1)-ee(m,indm)-ee(n,indn)
            den1=1.d0
            den2=1.d0


        Rm(i,j)=Rm(i,j)
     *       +zz*xvo_all(m,n,index2,i)*den1*xvo_all(m,n,index2,j)
     *       -zz*xvo_all(m,n,index2,i)*den2*xvo_all(m,n,index2,j+1)
28            continue 
27          continue			                     
26          continue 
25          continue
24          continue
23          continue 
22          continue

2000    continue
1000  continue

      do 212 i=1,nt
        do 213 j=1,nt
          Rm(i,j)=Rm(i,j)-al(i)
213       continue
212    continue


      write (*,*) 'Print alpha'
      do 801 i3=1,nt
       write (*,'(i4,f15.5)') i3,al(i3)
        bl(i3)=-al(i3)
801    continue
      write (*,*) 'Print R'
      do 802 i3=1,nt
       write (*,'(i4,11f15.5)') i3,(Rm(i3,i4),i4=1,nt)
802    continue

      call dgeco(Rm,nit,nt,ipvt,rcond,z)
      ising=0
      write (6,*) 'rcond=',rcond
      if (rcond.lt.1.d-18) ising=1
      if (ising.eq.1) then 
	 write (*,*) 'Singularity'
	 goto 901
      endif

      call dgesl(Rm,nit,nt,ipvt,bl,0)

	cc=0.d0
      write (*,*) 'Print coefficients'
      do 80 i3=1,nt
	 cc=cc+bl(i3)
         al(i3)=bl(i3)
       write (*,'(i4,f15.5)') i3,al(i3)
80    continue
      write (*,*)
90    format ('Total = ',f12.6)
      write (*,90) cc

********** Define new rho *************
      
      do 275 i=1,ivchan
        do 37 n=1,nxx
          do 78 m=1,nxx
              xvo(m,n,i)=0.d0
            do 85 i3=1,nt
              xvo(m,n,i)=xvo(m,n,i)+al(i3)*xvo_all(m,n,i,i3)
85          continue
78         continue
37       continue
275    continue


       do 718 n=1,nxx
        rvo(n)=0.d0
         do 905 i3=1,nt
           rvo(n)=rvo(n)+al(i3)*rvo_all(n,i3)
905       continue
718      continue

***************************************
        call energy1(res1)
        call energy2(res2)
        call energy3(res3)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>
        call energy4(evt,res4)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ene=res1+res2+res3+res4

      write (*,850) ene
850   format ('New energy = ',f15.8)
901   continue
      return
      end

      subroutine diis(evt,ene,deltae,max_rle)
      implicit double precision (a-h,o-z)
      integer*4 ipvt
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nxx)
      COMMON /all/ xvo_all(nxx,nxx,nvh,nit),rvo_all(nxx,nit)

      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)    
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      common /tau/ Rm(nit,nit),al(nit),bl(nit)
      dimension z(nit),ipvt(nit)
      COMMON /energysd/ delsd(nit)


      CALL klj(kv,kapv,lv,jv,indv,n0v)
      do 112 i=1,nit
        al(i)=0.d0
        bl(i)=0.d0
        do 113 j=1,nit
          Rm(i,j)=0.d0
113       continue
112    continue
       nt=max_rle-1

      write (*,*) 'Print energies' 
      write (*,*) 
      do 117 i=1,nt
        write (*,'(i5,f13.8)') i,delsd(i) 
117    continue
      write (*,*) 'DONE'
	write (*,*)

c ********* Accumulate  ak

      ak=0.d0



        do 812 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          do 813 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 814 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
   
              kmin=max0(iabs(kapm-kapv),iabs(kapn-kapb))
              kmax=min0((kapm+kapv-1),(kapn+kapb-1))
              DO 815 lk=kmin,kmax

                call odd (lm+lv+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 816
               
            index1=ipxv1(i2,indm,indn,lk)              

           zz=1.d0/((2.d0*lk+1.d0)*(jv+1.d0))
           do 817 m=nnm,nmaxx(indm)
             do 818 n=nnn,nmaxx(indn)               
                 ak=ak+zz*xv1(m,n,index1)*xv1(m,n,index1)  
818            continue 
817          continue
816       continue
815            continue 
814          continue
813       continue
812          continue
            write (*,*) 'ak=',ak

c 
c ********* Accumulate  al


      do 1600 i=1,nt


        do 12 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          do 13 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 14 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
   
              kmin=max0(iabs(kapm-kapv),iabs(kapn-kapb))
              kmax=min0((kapm+kapv-1),(kapn+kapb-1))
              DO 15 lk=kmin,kmax

                call odd (lm+lv+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 16
               
            index1=ipxv1(i2,indm,indn,lk)
            index2=ipxv(i2,indm,indn,lk)                 

           zz=1.d0/((2.d0*lk+1.d0)*(jv+1.d0))
           do 17 m=nnm,nmaxx(indm)
             do 18 n=nnn,nmaxx(indn)   
                den=1.d0
c               den=ev+eb+delsd(i)-ee(m,indm)-ee(n,indn)

                 al(i)=al(i)-zz*den*xv1(m,n,index1)*
     *       (xvo_all(m,n,index2,i)
     *       -xvo_all(m,n,index2,i+1))

18            continue 
17          continue
16       continue
15            continue 
14          continue
13       continue
12          continue
1600  continue



c ********* Accumulate Rm 

  
      do 1000 i=1,nt
       do 2000 j=1,nt


         indn=indv
         CALL indk1(indn,kn,kapn,ln,jn,n0)
         call st (kn,nnn) 
         do 61 n=nnn,nmaxx(indn)	
	     if (n+n0.ne.nv) then       
c           den1=ev+delsd(i)-ee(n,indn)
c           den2=ev+delsd(j)-ee(n,indn)
            den1=1.d0
            den2=1.d0

           Rm(i,j)=Rm(i,j)+den1*den2*
     *               (rvo_all(n,i)*rvo_all(n,j)
     *               +rvo_all(n,i+1)*rvo_all(n,j+1)
     *               -rvo_all(n,i+1)*rvo_all(n,j)
     *               -rvo_all(n,i)*rvo_all(n,j+1))
	    endif
61     continue

        do 22 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          do 23 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 24 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************    
              kmin=max0(iabs(kapm-kapv),iabs(kapn-kapb))
              kmax=min0((kapm+kapv-1),(kapn+kapb-1))
              DO 25 lk=kmin,kmax

                call odd (lm+lv,ii1)
                call odd (ln+lb,ii2)
                if ((ii1.eq.0.and.ii2.ne.0).or.
     *              (ii1.ne.0.and.ii2.eq.0)) goto 26

            index2=ipxv(i2,indm,indn,lk)                 
    
           zz=1.d0/((2.d0*lk+1.d0)*(jv+1.d0))
           do 27 m=nnm,nmaxx(indm)
             do 28 n=nnn,nmaxx(indn)

            den1=1.d0
            den2=1.d0

c               den1=ev+eb+delsd(i)-ee(m,indm)-ee(n,indn) 
c               den2=ev+eb+delsd(j)-ee(m,indm)-ee(n,indn)
        Rm(i,j)=Rm(i,j)
     *    +zz*den1*den2*
     *   (xvo_all(m,n,index2,i)*xvo_all(m,n,index2,j)
     *   +xvo_all(m,n,index2,i+1)*xvo_all(m,n,index2,j+1)
     *   -xvo_all(m,n,index2,i+1)*xvo_all(m,n,index2,j)
     *   -xvo_all(m,n,index2,i)*xvo_all(m,n,index2,j+1))

28            continue 
27          continue			                     
26          continue 
25          continue
24          continue
23          continue 
22          continue

2000    continue
1000  continue
899   format ('ak=',f20.10)
       write (*,899) ak

      do 214 i=1,nt
          al(i)=al(i)-ak
214    continue

      write (*,*) 'Print alpha'
      do 801 i3=1,nt
       write (*,'(i4,f15.7)') i3,al(i3)
        bl(i3)=-al(i3)
801    continue

      do 212 i=1,nt
        do 213 j=1,nt
          Rm(i,j)=Rm(i,j)-al(i)-al(j)-ak
213       continue
212    continue



      write (*,*) 'Print R'
      do 802 i3=1,nt
       write (*,'(i4,11f15.7)') i3,(Rm(i3,i4),i4=1,nt)
802    continue

      call dgeco(Rm,nit,nt,ipvt,rcond,z)
      ising=0
      write (6,*) 'rcond=',rcond
      if (rcond.lt.1.d-18) ising=1
      if (ising.eq.1) then 
	 write (*,*) 'Singularity'
	 goto 901
      endif

      call dgesl(Rm,nit,nt,ipvt,bl,0)

	cc=0.d0
      write (*,*) 'Print coefficients'
      do 80 i3=1,nt
	 cc=cc+bl(i3)
         al(i3)=bl(i3)
       write (*,'(i4,f15.5)') i3,al(i3)
80    continue
      write (*,*)
90    format ('Total = ',f12.6)
      write (*,90) cc

********** Define new rho *************
      
      do 275 i=1,ivchan
        do 37 n=1,nxx
          do 78 m=1,nxx
              xvo(m,n,i)=0.d0
            do 85 i3=1,nt
              xvo(m,n,i)=xvo(m,n,i)+al(i3)*xvo_all(m,n,i,i3)
85          continue
78         continue
37       continue
275    continue


       do 718 n=1,nxx
        rvo(n)=0.d0
         do 905 i3=1,nt
           rvo(n)=rvo(n)+al(i3)*rvo_all(n,i3)
905       continue
718      continue

***************************************
        call energy1(res1)
        call energy2(res2)
        call energy3(res3)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>
        call energy4(evt,res4)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ene=res1+res2+res3+res4

      write (*,850) ene
850   format ('New energy = ',f15.8)
901   continue
      return
      end



      subroutine energy1(res)
      implicit double precision (a-h,o-z)
      include "all.par" 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      res=0.d0
      CALL klj(kv,kapv,lv,jv,indv,n0v)
      do 4 j=1,ncore
        kb=ko(j)
        nb=no(j)
* >>>>>>> SIGMA >>>>>>>>>>>
        eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        do 2 indm=1,2*lmax+1
          CALL indk1(indm,km,kapm,lm,jm,m0)
          call st (km,nnm)
          if (km.eq.kb) then
           k=0
           call odd (lv+lv,i5)
           call odd (lm+lb,i6)
           if (i5*i6.ne.0) then
            z1=(-1)**(kapv+kapv+kapb+kapm)
            z2=SQRT((jb+1.d0)/(jv+1.d0))
            index1=ipxv1(j,indv,indm,k)
            do 3 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>
	        eu=1.d0/(eb-ee(m,indm))
              res=res+z1*z2*xv1(nv-n0v,m,index1)*rmi(m,j)*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>
3           continue
           endif
          endif
2       continue
4     continue
      return
      end

      subroutine energy2(res)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      res=0.d0
      CALL klj(kv,kapv,lv,jv,indv,n0v)
      do 4 j=1,ncore
        kb=ko(j)
        nb=no(j)
* >>>>>>> SIGMA >>>>>>>>>>>
        eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        do 2 indm=1,2*lmax+1
          CALL indk1(indm,km,kapm,lm,jm,m0)
          call st (km,nnm)
          if (km.eq.kb) then
           kmin=max0(iabs(kapm-kapv),iabs(kapv-kapb))
           kmax=min0((kapm+kapv-1),(kapv+kapb-1))
           DO 111 k=kmin,kmax
             call odd (lm+lv+k,i1)
             call odd (lv+lb+k,i2)
             if (i1*i2.eq.0) goto 100
             z1=(-1)**(kapv+kapv+kapb+kapm+kapv+kapb-1+k)
             z2=1.d0/(jv+1.d0)
             index1=ipxv1(j,indm,indv,k)
             do 3 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>
	        eu=1.d0/(eb-ee(m,indm))
               res=res+z1*z2*xv1(m,nv-n0v,index1)*rmi(m,j)*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>
3             continue
100          continue
111        continue
          endif
2       continue
4     continue
      return
      end

      subroutine energy3(res)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      res=0.d0
      CALL klj(kv,kapv,lv,jv,indv,n0v)
      do 3 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>> SIGMA >>>>>>>>>>>
        ea=eo(i)
        eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>
          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            kmin=max0(iabs(kapm-kapa),iabs(kapv-kapb))
            kmax=min0((kapm+kapa-1),(kapv+kapb-1))
            DO 11 k=kmin,kmax
              call odd (lm+la+k,i1)
              call odd (lv+lb+k,i2)

              if (i1*i2.eq.0) goto 555
              if (i.ge.j) then
               index1=ipx1(i,j,indm,indv,k)
               index2=ipxc(i,j,indm,indv,k)
              else
               index1=ipx1(j,i,indv,indm,k)
               index2=ipxc(j,i,indv,indm,k)
              endif
              z1=(-1)/((2*k+1.d0)*(jv+1.d0))
              do 38 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
               eu=1.d0/(ea+eb-ee(m,indm)-ev)
               if (i.ge.j ) then
                res=res+z1*xci(m,nv-n0v,index2)*x1(m,nv-n0v,index1)*eu
               else
                res=res+z1*xci(nv-n0v,m,index2)*x1(nv-n0v,m,index1)*eu
               endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
38             continue
**********************************************************
              k1min=max0(iabs(kapm-kapb),iabs(kapa-kapv))
              k1max=min0((kapm+kapb-1),(kapa+kapv-1))
              DO 21 k1=k1min,k1max

                call odd (lm+lb,i1)
                call odd (lv+la,i2)

                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 515
              if (i.ge.j) then
               index1=ipx1(i,j,indm,indv,k)
               index2=ipxc(i,j,indv,indm,k1)
              else
               index1=ipx1(j,i,indv,indm,k)
               index2=ipxc(j,i,indm,indv,k1)
              endif
                zz=(2*k+1.d0)*d6j(jm,ja,2*k,jv,jb,2*k1)
                do 31 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
               eu=1.d0/(ea+eb-ee(m,indm)-ev)
               if (i.ge.j ) then
                    res=res+zz*z1*xci(nv-n0v,m,index2)*
     *                            x1(m,nv-n0v,index1)*eu
                else
                    res=res+zz*z1*xci(m,nv-n0v,index2)*
     *                            x1(nv-n0v,m,index1)*eu
               endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
31               continue
515             continue
21            continue
555             continue
11            continue
2           continue
4         continue
3       continue
901   continue
      return
      end



      subroutine outrho(delta,icount_val2,evt,kval)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /dsigma/ drsigma(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /sigma/ rsigma(nx)

      ivv=1
      write (9) nv,kv,ev,delta
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>
	write (9) kval,evt
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      write (9) nmax,lmax,ivchan
      write (9) ivv,ivb,ivm,ivn,ivl
      write (9) ipxv

      do 27 i=1,ivchan
        indn=ivn(i)
        indm=ivm(i)
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)
        do 37 n=1,nmaxn
          write (9) (xvo(m,n,i),m=1,nmaxm)
37       continue
* >>>>> SIGMA >>>>>>>>>>
        do 317 n=1,nmaxn
          write (9) (dxvo(m,n,i),m=1,nmaxm)
317       continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>
27     continue

* >>>>> SIGMA >>>>>>>>>>
c     write (9) (rvo(n),n=1,nmax)
      write (9) (rsigma(n),n=1,nmax)
      write (9) (drsigma(n),n=1,nmax)
* >>>>>>>>>>>>>>>>>>>>>
      icount_val2=icount_val2+ivchan
      return
      end
c



      subroutine xcore
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2
      DIMENSION ga(nhf),fa(nhf),gb(nhf),fb(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /x12ini/ inmax
      lin=2*lmax+1
      if (lin.lt.inmax) then
       lin=inmax
      endif
c      write (*,*) 'lin =',lin
      index=0
      index1=0
      do 531 i=1,nxk
        do 532 m=1,nxx
          do 533 n=1,nxx
            x1(m,n,i)=0.d0
            x2(m,n,i)=0.d0

533       continue
532     continue
531   continue

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 12 il=1,max
          ga(il)=go(il,i)
          fa(il)=fo(il,i)
12      CONTINUE
        ea=eo(i)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          DO 13 il=1,max
            gb(il)=go(il,j)
            fb(il)=fo(il,j)
13        CONTINUE
          eb=eo(j)
          if (i.ge.j) then

          do 191 ij=1,max
            v2(ij)=ga(ij)*gb(ij)+fa(ij)*fb(ij)
191       continue

          do 2 indm=1,lin
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,lin
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnab) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                av=0
                index=index+1
                ipx1(i,j,indm,indn,k)=index
                c=((-1)**(k))*s(k,km,ka)*s(k,kn,kb)
                do 3 m=1,nmaxx(indm)
                  do 19 ij=1,max
                    v(ij)=ga(ij)*gg(ij,m,indm)+fa(ij)*ff(ij,m,indm)
19                continue
                  call yfun(v,u,k,max,*901)
                  do 6 n=1,nmaxx(indn)
                    do 20 ij=1,max
                      vv(ij)=(gb(ij)*gg(ij,n,indn)+
     *                        fb(ij)*ff(ij,n,indn))*rp(ij)*u(ij)
20                  continue
                    x1(m,n,index)=c*rint(vv,1,max,11,h)
                    av=av+x1(m,n,index)*x1(m,n,index)
6                 continue
3               continue
555             continue
11            continue
************************************************************
****************************************************************
*             Calculate  Xk(manb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
              kmax=min0((kapm+kapn-1),(kapa+kapb-1))
              DO 111 k=kmin,kmax
                call odd (lm+ln+k,i1)
                call odd (la+lb+k,i2)
                if (i1*i2.eq.0) goto 515
                av1=0
                index1=index1+1
                ipx2(i,j,indm,indn,k)=index1
                c=((-1)**(k))*s(k,km,kn)*s(k,ka,kb)
                call yfun(v2,u,k,max,*901)
                do 31 m=1,nmaxx(indm)
                  do 61 n=1,nmaxx(indn)
                    do 201 ij=1,max
                      vv(ij)=(gg(ij,m,indm)*gg(ij,n,indn)+
     *                        ff(ij,m,indm)*ff(ij,n,indn))*rp(ij)*u(ij)
201                  continue
                    x2(m,n,index1)=c*rint(vv,1,max,11,h)
c                    av1=av1+x2(m,n,index)*x2(m,n,index)
61                 continue
31               continue

515             continue
111            continue
************************************************************

5           continue
2         continue
187       format ('CHANNEL',i6,'  j1=',i3,'  j2=',i3,'  km=',i3,
     *    '  kn=',i3,'  k=',i3,f12.6)
c           write (*,187) index,i,j,km,kn,k,av
        endif
4       continue
1     continue

901   continue
      return
      end




      subroutine xval
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /xval2/ ipxv2(nn,nk,nk,0:kk),xv2(nxx,nxx,nxv)

      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /x12ini/ inmax
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

      lin=2*lmax+1
      if (lin.lt.inmax) then
       lin=inmax
      endif
c      write (*,*) 'lin =',lin

      index=0
      index1=0
      do 531 i=1,nxv
        do 532 m=1,nxx
          do 533 n=1,nxx
            xv1(m,n,i)=0.d0
            xv2(m,n,i)=0.d0

533       continue
532     continue
531   continue


        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 12 il=1,max
          ga(il)=gv(il)
          fa(il)=fv(il)
12      CONTINUE
        ea=ev
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          DO 13 il=1,max
            gb(il)=go(il,j)
            fb(il)=fo(il,j)
13        CONTINUE
          eb=eo(j)

          do 191 ij=1,max
            v2(ij)=ga(ij)*gb(ij)+fa(ij)*fb(ij)
191       continue

          do 2 indm=1,lin
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,lin
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnvb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in xval -1 , STOP'
	stop
      endif
		   

              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                av=0

                index=index+1
	          if (index.gt.nxv) then 
      write (*,*) 'Number of x1 channels in valsd exceeds NXV, STOP'
                 STOP
	          endif

                ipxv1(j,indm,indn,k)=index
                c=((-1)**(k))*s(k,km,ka)*s(k,kn,kb)
                do 3 m=1,nmaxx(indm)
                  do 19 ij=1,max
                    v(ij)=ga(ij)*gg(ij,m,indm)+fa(ij)*ff(ij,m,indm)
19                continue
                  call yfun(v,u,k,max,*901)
                  do 6 n=1,nmaxx(indn)
                    do 20 ij=1,max
                      vv(ij)=(gb(ij)*gg(ij,n,indn)+
     *                        fb(ij)*ff(ij,n,indn))*rp(ij)*u(ij)
20                  continue
c                 write (*,*) m,n,index
                    xv1(m,n,index)=c*rint(vv,1,max,11,h)
                    av=av+xv1(m,n,index)*xv1(m,n,index)
6                 continue
3               continue
c          write (*,187) index,j,km,kn,k,av
555             continue
11            continue
************************************************************
****************************************************************
*             Calculate  Xk(mvnb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
              kmax=min0((kapm+kapn-1),(kapa+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in xval -2 , STOP'
	stop
      endif

              DO 111 k=kmin,kmax
                call odd (lm+ln+k,i1)
                call odd (la+lb+k,i2)
                if (i1*i2.eq.0) goto 515
                av1=0
                index1=index1+1
	          if (index1.gt.nxv) then 
      write (*,*) 'Number of x2 channels in valsd exceeds NXV, STOP'
                 STOP
	          endif

                ipxv2(j,indm,indn,k)=index1
                c=((-1)**(k))*s(k,km,kn)*s(k,ka,kb)
                call yfun(v2,u,k,max,*901)
                do 31 m=1,nmaxx(indm)
                  do 61 n=1,nmaxx(indn)
                    do 201 ij=1,max
                      vv(ij)=(gg(ij,m,indm)*gg(ij,n,indn)+
     *                        ff(ij,m,indm)*ff(ij,n,indn))*rp(ij)*u(ij)
201                  continue
                    xv2(m,n,index1)=c*rint(vv,1,max,11,h)
c                    av1=av1+xv2(m,n,index)*xv2(m,n,index)
61                 continue
31               continue
515             continue


111            continue
************************************************************

5           continue
2         continue
187       format ('CHANNEL',i6,'  j2=',i3,'  km=',i3,
     *    '  kn=',i3,'  k=',i3,f12.6)
c          write (*,187) index,i,j,km,kn,k,av
4       continue

c      write (*,*) index,index1

901   continue
      return
      end

      subroutine countx
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /x12ini/ inmax
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

      lin=2*lmax+1
      if (lin.lt.inmax) then
       lin=inmax
      endif
c      write (*,*) 'lin =',lin
        index=0
	index1=0
        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)

          do 2 indm=1,lin
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,lin
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnvb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in xval -1 , STOP'
	stop
      endif
		   

              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                index=index+1
555             continue
11            continue
************************************************************
****************************************************************
*             Calculate  Xk(mvnb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
              kmax=min0((kapm+kapn-1),(kapa+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in xval -2 , STOP'
	stop
      endif
              DO 111 k=kmin,kmax
                call odd (lm+ln+k,i1)
                call odd (la+lb+k,i2)
                if (i1*i2.eq.0) goto 515
                index1=index1+1
515             continue
111            continue
************************************************************

5           continue
2         continue
4       continue
**************************************************
      write (*,*)
      write (*,899)nv,kv
899   format ('PRE-COUNTING ALL VALENCE CHANNELS FOR n = ',i2,1x,
     *'kappa = ',i2) 
	write (*,*)
	if (index.le.NXV.and.index1.le.NXV) then 
	 write (*,122) index,nxv
	 write (*,123) index1,nxv
      else
	 write (*,124) index,nxv
	 write (*,125) index1,nxv
       write (*,*) 'Parameter NXV is exceeded, STOP'
	 stop
	endif
122   format ('Number of X_k(mnvb) channels      = ',i7,
     *' < NXV = ',i7,' OK')
123   format ('Number of X_k(mvnb) channels      = ',i7,
     *' < NXV = ',i7,' OK')

124   format ('Number of X_k(mnvb) channels      = ',i7,
     *' > NXV = ',i7)
125   format ('Number of X_k(mvnb) channels      = ',i7,
     *' > NXV = ',i7)



      return
      end


      subroutine rhocore
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)

      do 531 i=1,nch
        do 532 m=1,nxx
          do 533 n=1,nxx
            xci(m,n,i)=0.d0
533       continue
532     continue
531   continue

      OPEN (unit=7,file='pair.3',form='unformatted')

      read (7) nnmax,llmax,icchan,jlo,jhi
      read (7) ipa,ipb,ipm,ipn,ipl
      read (7) ipxc

      do 27 i=1,icchan
        indn=ipn(i)
        indm=ipm(i)
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)
        do 37 n=1,nmaxn
          read (7) (xci(m,n,i),m=1,nmaxm)
37       continue
27     continue

      do 47 i=jlo,jhi
          read (7) (rmi(n,i),n=1,nmax)
47     continue
      close (7)
      if (llmax.lt.lmax.or.nnmax.lt.nmax) then
       write (*,*) 'WRONG LMAX OR NMAX IN THE CORE FILE'
       stop
      endif
**** SIGMA >>>>>>>>>>>> Make sigmas from core rhos >>>>>>>>>>>>>>>
      do 217 i=1,icchan
        indn=ipn(i)
        indm=ipm(i)
	  ia=ipa(i)
	  ib=ipb(i)
        ea=eo(ia)
	  eb=eo(ib) 
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)
        do 317 m=1,nmaxm
        do 327 n=1,nmaxn
            xci(m,n,i)=xci(m,n,i)*(ea+eb-ee(m,indm)-ee(n,indn))
327       continue
317       continue
217     continue

      do 417 i=jlo,jhi
        ea=eo(i)
	  ka=ko(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        indn=inda
	  do 500 n=1,nmax
          rmi(n,i)=rmi(n,i)*(ea-ee(n,indn))
500     continue
417     continue
*** SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

901   continue
      return
      end

      subroutine rhoval(key,iv)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /drhovin/  dxvi(nxx,nxx,nvh),drvi(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      DIMENSION ivb9(NVH),ivm9(NVH),ivn9(NVH),ivl9(NVH)
      DIMENSION ipxv9(nn,nk,nk,0:kk)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

      do 531 i=1,nvh
        do 532 m=1,nxx
          do 533 n=1,nxx
            xvo(m,n,i)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            dxvo(m,n,i)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

533       continue
532     continue
531   continue

       index=0
        ka=kv
        na=nv

        CALL klj(ka,kapa,la,ja,inda,n0a)
         indv=inda
	  ea=ev

        do 80 m=1,nmaxx(inda)
         rvo(m)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            drvo(m)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 80     continue

        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(j)

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
****************************************************************
*            Initialize  rhok(mnab)                            *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in rhoval , STOP'
	stop
      endif

              DO 11 k=kmin,kmax
c              write (*,'(11i3)') jm,ja,jn,jb,k,lm,la,ln,lb,lm+la,ln+lb
                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 551

                index=index+1
	          if (index.gt.nvh) then 
      write (*,*) 'Number of rho channels in valsd exceeds NVH, STOP'
                 STOP
	          endif

                ipxv(j,indm,indn,k)=index
                ivv=1
                ivb(index)=j
                ivm(index)=indm
                ivn(index)=indn
                ivl(index)=k

                do 732 m=1,nxx
                  do 733 n=1,nxx
                    xvo(m,n,index)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    dxvo(m,n,index)=0.d0
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

733               continue
732             continue

                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)

                if (i1*i2.eq.0) goto 555
                index1=ipxv1(j,indm,indn,k)
                av=0
                do 3 m=nnm,nmaxx(indm)
                  do 6 n=nnn,nmaxx(indn)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*                    eu=1.d0/(ea+eb-ee(m,indm)-ee(n,indn))
                  	eu=1.d0
                    dxvo(m,n,index)=0.d0
*>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    xvo(m,n,index)=xv1(m,n,index1)*eu

6                 continue

3               continue
555             continue
551             continue
11            continue
5           continue
2         continue
4       continue
1     continue
      ivchan=index
      iv=ivchan
c      write (*,*) 'Number of channels = ',ivchan

            if (key.eq.1) then
	 write (*,*) ' THIS PROGRAM CAN NOT BE RESTARTED WITH KEY=1'
	stop
	endif


901   continue
      return
      end

      subroutine countrho
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

       index=0
        ka=kv
        na=nv

        CALL klj(ka,kapa,la,ja,inda,n0a)
         indv=inda
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
****************************************************************
*            Initialize  rhok(mnab)                            *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in rhoval , STOP'
	stop
      endif

              DO 11 k=kmin,kmax
                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 551

                index=index+1
551             continue
11            continue
5           continue
2         continue
4       continue
1     continue

	if (index.le.NVH) then 
	 write (*,122) index,nvh
      else
	 write (*,124) index,nvh
       write (*,*) 'Parameter NVH is exceeded, STOP'
	 stop
	endif
122   format ('Number of rho_k(mnvb) channels    = ',i7,
     *' < NVH = ',i7,' OK')

124   format ('Number of rho_k(mnvb) channels    = ',i7,
     *' > NVH = ',i7)

      write (*,*)


      return
      end


      subroutine energy4(evt,res)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv

       res=0.d0
        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=ev
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(j)
          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*            Initialize  rhok(mnab)                            *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax

                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)

                if (i1*i2.eq.0) goto 555
                 index1=ipxv1(j,indm,indn,k)
                 index2=ipxv(j,indm,indn,k)

                zz=1.d0/((ja+1.d0)*(2*k+1.d0))
                do 3 m=nnm,nmaxx(indm)
                  do 6 n=nnn,nmaxx(indn)
*>>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    eu=1.d0/(evt+eb-ee(m,indm)-ee(n,indn))
                    res=res+zz*xvo(m,n,index2)*xv1(m,n,index1)*eu
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
6                 continue
3               continue
********************************************************8
              k1min=max0(iabs(kapm-kapb),iabs(kapa-kapn))
              k1max=min0((kapm+kapb-1),(kapa+kapn-1))
              DO 21 k1=k1min,k1max

                call odd (lm+lb,i1)
                call odd (ln+la,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 515
                 index1=ipxv1(j,indm,indn,k)
                 index2=ipxv(j,indn,indm,k1)

                z1=(2*k+1.d0)*d6j(jm,ja,2*k,jn,jb,2*k1)
                do 31 m=nnm,nmaxx(indm)
                  do 61 n=nnn,nmaxx(indn)
*>>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    eu=1.d0/(evt+eb-ee(m,indm)-ee(n,indn))
                  res=res+zz*z1*xvo(n,m,index2)*xv1(m,n,index1)*eu
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
61                 continue
31               continue
515             continue
21            continue
555             continue
11            continue
5           continue
2         continue
4       continue
901   continue
      return
      end


      subroutine outin
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /drhovin/  dxvi(nxx,nxx,nvh),drvi(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)

        do 2 m=1,nmax
          rvi(m)=rvo(m)
          rvo(m)=0.d0
* >>>>>>>>>> SIGMA >>>>>>>>>
          drvi(m)=drvo(m)
          drvo(m)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>
2       continue
      do 3 i=1,ivchan
        indn=ivn(i)
        indm=ivm(i)
        do 4 n=1,nmaxx(indn)
          do 5 m=1,nmaxx(indm)
            xvi(m,n,i)=xvo(m,n,i)
            xvo(m,n,i)=0.d0
* >>>>>>>>>> SIGMA >>>>>>>>>>
            dxvi(m,n,i)=dxvo(m,n,i)
            dxvo(m,n,i)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>
5         continue
4       continue

3     continue
      return
      end

      subroutine term1(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do 30 index=1,ivchan
        j=ivb(index)
        indm=ivm(index)
        indn=ivn(index)
        k=ivl(index)
        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        do 3 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 4 jj=1,ncore
            kd=ko(jj)
            nd=no(jj)
            CALL klj(kd,kapd,ld,jd,indd,n0d)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>
            ec=eo(ii)
            ed=eo(jj)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            l1min=max0(iabs(kapc-kapa),iabs(kapd-kapb))
            l1max=min0((kapc+kapa-1),(kapd+kapb-1))
            DO 5 l=l1min,l1max
              call odd (la+lc+l,i1)
              call odd (lb+ld+l,i2)
              if (i1*i2.eq.0) goto 515
              index1=ipxv1(j,indc,indd,l)
              k1min=max0(iabs(kapm-kapc),iabs(kapn-kapd),iabs(k-l))
              k1max=min0((kapm+kapc-1),(kapn+kapd-1),(k+l))
              DO 9 k1=k1min,k1max
                call odd (lm+lc,i1)
                call odd (ln+ld,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 535
                z1=(-1)**(kapm+kapn+kapa+kapb)
                z2=d6j(jm,ja,2*k,2*l,2*k1,jc)
                z3=z1*z2*(2*k+1)*d6j(jn,jb,2*k,2*l,2*k1,jd)
                if (ii.ge.jj) then
                 index2=ipxc(ii,jj,indm,indn,k1)
                else
                 index2=ipxc(jj,ii,indn,indm,k1)
                endif
                do 10 n=nnn,nmaxx(indn)
                  do 11 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         eu=1.d0/(ec+ed-ee(m,indm)-ee(n,indn)+evt-ev)
                    if (ii.ge.jj) then
         tt=z3*xv1(nc-n0c,nd-n0d,index1)*xci(m,n,index2)*eu
         xvo(m,n,index)=xvo(m,n,index)+tt
         dxvo(m,n,index)=dxvo(m,n,index)-tt*eu
                    else
	   tt=z3*xv1(nc-n0c,nd-n0d,index1)*xci(n,m,index2)*eu
         xvo(m,n,index)=xvo(m,n,index)+tt
         dxvo(m,n,index)=dxvo(m,n,index)-tt*eu
                    endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

11                continue
10              continue
535             continue
9             continue
515           continue
5           continue
4         continue
3       continue
30    continue
      return
      end


      subroutine outsigma
      include "all.par"
* >>>>>> SIGMA >>>>>>>>>>>>>
      COMMON /dsigma/ drsigma(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /sigma/ rsigma(nx)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      common /nmaxx/nmaxx(nx)

        CALL klj(kv,kapv,lv,jv,indv,n0v)
        indm=indv
        CALL indk1(indm,km,kapm,lm,jm,m0)
        call st (km,nnm)
        do 200 m=nnm,nmaxx(indm)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	write (17,'(3i4,2e20.8)') nv,kv, m+m0,rsigma(m),drsigma(m) 
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
200     continue
      return
	 end


      subroutine outt(deltae)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      COMMON /sigma/ rsigma(nx)
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /dsigma/ drsigma(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      del=deltae
c      write (*,*) del

        CALL klj(kv,kapv,lv,jv,indv,n0v)
        indm=indv
c        write (*,*) 'rho =', rvo(6)
        CALL indk1(indm,km,kapm,lm,jm,m0)
        call st (km,nnm)
        do 200 m=nnm,nmaxx(indm)
	  rsigma(m)=rvo(m)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>
	  drsigma(m)=drvo(m)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>

          if (m+m0.ne.nv) then
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c              eu=1.d0/(del+ev-ee(m,indm))
	       eu=1.d0
             rvo(m)=rvo(m)*eu
             drvo(m)=drvo(m)*eu
          else
            rvo(m)=0.d0
            drvo(m)=0.d0

          endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

200     continue

      do 3 iu=1,ivchan
        j=ivb(iu)
        indm=ivm(iu)
        indn=ivn(iu)
        k=ivl(iu)
        kb=ko(j)
        nb=no(j)

        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        CALL klj(kv,kapv,lv,jv,indv,n0v)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        call st (km,nnm)
        call st (kn,nnn)
        eb=eo(j)
        index=ipxv1(j,indm,indn,k)
        av=0

        av=0.d0
        call odd (lv+lm+k,i1)
        call odd (lb+ln+k,i2)
        if (i1*i2.eq.0) then
         do 10 n=nnn,nmaxx(indn)
           do 11 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c             eu=1.d0/(del+ev+eb-ee(m,indm)-ee(n,indn))
	        eu=1.d0
             xvo(m,n,iu)=xvo(m,n,iu)*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
11         continue
10       continue
        else
         do 101 n=nnn,nmaxx(indn)
           do 111 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c             eu=1.d0/(del+ev+eb-ee(m,indm)-ee(n,indn))
	        eu=1.d0
             xvo(m,n,iu)=(xv1(m,n,index)+xvo(m,n,iu))*eu
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
111        continue
101      continue
        endif
3     continue
      return
      end

      subroutine term51(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do 30 index=1,ivchan
        j=ivb(index)
        indm=ivm(index)
        indn=ivn(index)
        k=ivl(index)
        CALL klj(kv,kapv,lv,jv,indv,n0v)
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        do 40 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
* >>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>
          ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          call odd (lv+lc+k,i1)
          call odd (lb+ln+k,i2)
          if (kc.eq.km.and.(i1*i2).ne.0) then
           iu=ipxv1(j,indc,indn,k)
           do 4 n=nnn,nmaxx(indn)
             do 41 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(ec-ee(m,indm)+evt-ev)
	          tt=-xv1(nc-n0c,n,iu)*rmi(m,ii)*eu
                xvo(m,n,index)=xvo(m,n,index)+tt
                dxvo(m,n,index)=dxvo(m,n,index)-tt*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
41           continue
4          continue
          endif
40      continue
30    continue
901   continue
      return
      end

      subroutine term52(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 30 index=1,ivchan
        j=ivb(index)
        indm=ivm(index)
        indn=ivn(index)
        k=ivl(index)
        CALL klj(kv,kapv,lv,jv,indv,n0v)
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        do 40 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
* >>>>>>> SIGMA >>>>>>>>>>
          ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          call odd (lb+lc+k,i1)
          call odd (lm+lv+k,i2)
          if (kc.eq.kn.and.(i1*i2).ne.0) then
           iu=ipxv1(j,indm,indc,k)
           do 4 n=nnn,nmaxx(indn)
             do 41 m=nnm,nmaxx(indm)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(ec-ee(n,indn)+evt-ev)
                tt=-xv1(m,nc-n0c,iu)*rmi(n,ii)*eu
                xvo(m,n,index)=xvo(m,n,index)+tt
                dxvo(m,n,index)=dxvo(m,n,index)-tt*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
41           continue
4          continue
          endif
40      continue
30    continue
901   continue
      return
      end

      subroutine terma1(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)

        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>
          eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
           if (indn.eq.indb.and.indm.eq.inda) then
          CALL indk1(indm,km,kapm,lm,jm,m0)
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (km,nnm)
          call st (kn,nnn)
          k=0
          call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
          call trgi(iabs(kapn-kapb),(kapn+kapb-1),k,i2)
          call odd (lb+ln,i5)
          call odd (lm+la,i6)
          if (i1*i2*i5*i6.ne.0) then
           index1=ipxv1(j,indm,indn,k)

          z1=((-1)**(kapb-kapn))*SQRT((jb+1.d0)/(ja+1.d0))
          do 4 n=nnn,nmaxx(indn)
            do 15 m=nnm,nmaxx(indm)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(eb-ee(n,indn)+evt-ev)
	          tt=z1*xv1(m,n,index1)*rmi(n,j)*eu
                rvo(m)=rvo(m)+tt
                drvo(m)=drvo(m)-tt*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
15          continue
4         continue
700        continue
         endif
        endif
5       continue

2      continue
40      continue
      return
      end


      subroutine terma2(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /xval2/ ipxv2(nn,nk,nk,0:kk),xv2(nxx,nxx,nxv)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>
          eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
           if (indn.eq.indb.and.indm.eq.inda) then

          CALL indk1(indm,km,kapm,lm,jm,m0)
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (km,nnm)
          call st (kn,nnn)
          kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
          kmax=min0((kapm+kapn-1),(kapa+kapb-1))
          DO 111 k=kmin,kmax
            call odd (lm+ln+k,i1)
            call odd (la+lb+k,i2)
            if (i1*i2.eq.0) goto 100

             index1=ipxv2(j,indm,indn,k)

            z1=((-1)**(kapa+kapb+k-1+kapb-kapa))/(ja+1.d0)
            do 4 n=nnn,nmaxx(indn)
              do 15 m=nnm,nmaxx(indm)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(eb-ee(n,indn)+evt-ev)
                tt=z1*xv2(m,n,index1)*rmi(n,j)*eu
                 rvo(m)=rvo(m)+tt
                 drvo(m)=drvo(m)-tt*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

15            continue
4           continue
100         continue
111       continue
        endif
5       continue

2      continue

40      continue
c     write (*,'(i5,e23.10)') i,rvo(6,i)

      return
      end

      subroutine terma3(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan2
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)
* >>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
          eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 25 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            if (indm.eq.inda) then
          do 3 ii=1,ncore
            kc=ko(ii)
            nc=no(ii)
* >>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
            ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            CALL klj(kc,kapc,lc,jc,indc,n0c)
            do 2 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0n)
              call st (kr,nnr)
              call st (kn,nnn)
              kmin=max0(iabs(kapc-kapn),iabs(kapa-kapb))
              kmax=min0((kapc+kapn-1),(kapa+kapb-1))
              DO 111 k=kmin,kmax

                call odd (lc+ln+k,i1)
                call odd (la+lb+k,i2)
                if (i1*i2.eq.0) goto 100

                if (ii.ge.j) then
                 index1=ipx1(ii,j,indn,inda,k)
                 index2=ipxc(ii,j,indn,indm,k)
                else
                 index1=ipx1(j,ii,inda,indn,k)
                 index2=ipxc(j,ii,indm,indn,k)
                endif

                z1=(-1)/((ja+1.d0)*(2*k+1.d0))
                do 4 n=nnn,nmaxx(indn)
                  do 15 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
                eu=1.d0/(ec+eb-ee(m,indm)-ee(n,indn)+evt-ev)
                    if (ii.ge.j) then
	               tt=z1*x1(n,na-n0a,index1)*
     *                                   xci(n,m,index2)*eu
                     rvo(m)=rvo(m)+tt
                     drvo(m)=drvo(m)-tt*eu

                    else
	               tt=z1*x1(na-n0a,n,index1)*
     *                                   xci(m,n,index2)*eu
                     rvo(m)=rvo(m)+tt
                     drvo(m)=drvo(m)-tt*eu
                    endif
15                continue
4               continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                k1min=max0(iabs(kapb-kapn),iabs(kapm-kapc))
                k1max=min0((kapb+kapn-1),(kapm+kapc-1))
                DO 211 k1=k1min,k1max
                  call odd (lb+ln,i5)
                  call odd (lm+lc,i6)
                  if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 900

                  if (ii.ge.j) then
                   index2=ipxc(ii,j,indm,indn,k1)
                  else
                   index2=ipxc(j,ii,indn,indm,k1)
                 endif

                zz=z1*(2*k+1.d0)*d6j(jn,jc,2*k,jm,jb,2*k1)
                do 41 n=nnn,nmaxx(indn)
                  do 151 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(ec+eb-ee(m,indm)-ee(n,indn)+evt-ev)

                    if (ii.ge.j) then
	               tt=zz*x1(n,na-n0a,index1)*
     *                                   xci(m,n,index2)*eu
                     rvo(m)=rvo(m)+tt
                     drvo(m)=drvo(m)-tt*eu

                    else
	               tt=zz*x1(na-n0a,n,index1)*
     *                                   xci(n,m,index2)*eu
                     rvo(m)=rvo(m)+tt

                     drvo(m)=drvo(m)-tt*eu
                    endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
151                continue
41               continue
900             continue
211             continue
100             continue
111           continue
2           continue
3         continue
       endif
25      continue
40      continue
      return
      end


      subroutine term41(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      DIMENSION x(nx,nx),y(nx,nx)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovin/  dxvi(nxx,nxx,nvh),drvi(nx)
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
	DIMENSION dy(nx,nx) 
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
* >>>>>> SIGMA >>>>>>>>>>>
          ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 22 indr=1,2*lmax+1
            CALL indk1(indr,kr,kapr,lr,jr,n0r)
            call st (kr,nnr)
            do 21 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,m0)
              call st (km,nnm)
              kmin=max0(iabs(kapr-kapc),iabs(kapm-kapa))
              kmax=min0((kapr+kapc-1),(kapm+kapa-1))
              DO 500 k=kmin,kmax
                call odd (lc+lr,i5)
                call odd (lm+la,i6)
                if ((i5.eq.0.and.i6.ne.0).or.
     *              (i5.ne.0.and.i6.eq.0)) goto 900

***********************************************

              zf=((-1)**(k+kapr+kapc))/(2.d0*k+1)
              do 171 nr=nnr,nmaxx(indr)
                do 172 m=nnm,nmaxx(indm)
                  y(m,nr)=0.d0
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>
                  dy(m,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue
              index2=ipxv(ii,indm,indr,k)
              do 12 nr=nnr,nmaxx(indr)
                do 10 m=nnm,nmaxx(indm)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>
	             eu=1.d0/(evt+ec-ee(m,indm)-ee(nr,indr))
	             tt=xvi(m,nr,index2)*eu
	             tt1=dxvi(m,nr,index2)*eu
                   y(m,nr)=y(m,nr)+tt
                   dy(m,nr)=dy(m,nr)-tt*eu+tt1
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
10              continue
12            continue

              k1min=max0(iabs(kapm-kapc),iabs(kapr-kapa))
              k1max=min0((kapm+kapc-1),(kapr+kapa-1))
              DO 9 k1=k1min,k1max
                z2=(2*k+1)*d6j(jm,ja,2*k,jr,jc,2*k1)

                 index2=ipxv(ii,indr,indm,k1)

                do 121 nr=nnr,nmaxx(indr)
                  do 101 m=nnm,nmaxx(indm)
* >>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	             eu=1.d0/(evt+ec-ee(m,indm)-ee(nr,indr))
                   tt =z2*xvi(nr,m,index2)*eu
	             tt1=z2*dxvi(nr,m,index2)*eu
                   y(m,nr)=y(m,nr)+tt
                   dy(m,nr)=dy(m,nr)-tt*eu+tt1
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
101                continue
121              continue
9              continue
               do 33 j=1,ncore
                 kb=ko(j)
                 nb=no(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)
                 do 23 indn=1,2*lmax+1
                   CALL indk1(indn,kn,kapn,ln,jn,n0)
                   call st (kn,nnn)
                   call trgi(iabs(kapn-kapb),(kapn+kapb-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (ln+lb,i5)
                   call odd (lm+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700


                   do 165 n=nnn,nmaxx(indn)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,n)=0.d0
164                  continue
165                continue


                   call odd (lc+lr+k,i5)
                   call odd (ln+lb+k,i6)
                   if (i5*i6.eq.0) goto 800

                   if (ii.ge.j) then
                    index1=ipx1(ii,j,indr,indn,k)
                   else
                    index1=ipx1(j,ii,indn,indr,k)
                   endif
                   z2=zf*((-1)**(kapc-kapr))
                   do 123 n=nnn,nmaxx(indn)
                     do 103 nr=nnr,nmaxx(indr)
                       if (ii.ge.j) then
                        x(nr,n)=x(nr,n)+z2*x1(nr,n,index1)
                       else
                        x(nr,n)=x(nr,n)+z2*x1(n,nr,index1)
                       endif
103                  continue
123                continue
800                continue

                   zz=(-1)**(kapb+kapc+kapr+kapn)
                   k2min=max0(iabs(kapc-kapb),iabs(kapn-kapr))
                   k2max=min0((kapc+kapb-1),(kapn+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+lb+k2,i1)
                     call odd (ln+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155

                     if (ii.ge.j) then
                      index1=ipx2(ii,j,indn,indr,k2)
                     else
                      index1=ipx2(j,ii,indr,indn,k2)
                     endif
                     z1=zf*(2*k+1.d0)*d6j(jb,jn,2*k,jr,jc,2*k2)
                     do 124 n=nnn,nmaxx(indn)
                       do 104 nr=nnr,nmaxx(indr)
                         if (ii.ge.j) then
                          x(nr,n)=x(nr,n)+z1*x2(n,nr,index1)
                         else
                          x(nr,n)=x(nr,n)+z1*zz*x2(nr,n,index1)
                         endif
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipxv(j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
                        do 107 m=nnm,nmaxx(indm)
                         xvo(m,n,index)=xvo(m,n,index)+
     *                                  x(nr,n)*y(m,nr)

* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                         dxvo(m,n,index)=dxvo(m,n,index)+
     *                                  x(nr,n)*dy(m,nr)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
107                     continue
117                   continue
127                 continue
700              continue
23               continue
33             continue
900           continue
500           continue
21          continue
22        continue
31      continue
      return
      end

      subroutine term42(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /xval1/ ipxv1(nn,nk,nk,0:kk),xv1(nxx,nxx,nxv)
      COMMON /xval2/ ipxv2(nn,nk,nk,0:kk),xv2(nxx,nxx,nxv)
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      DIMENSION x(nx,nx),y(nx,nx)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
	DIMENSION dy(nx,nx) 
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 32 j=1,ncore
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
* >>>>>> SIGMA >>>>>>>>
         eb=eo(j)
	   ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 22 indr=1,2*lmax+1
            CALL indk1(indr,kr,kapr,lr,jr,n0r)
            call st (kr,nnr)
            do 21 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
              kmin=max0(iabs(kapr-kapc),iabs(kapn-kapb))
              kmax=min0((kapr+kapc-1),(kapn+kapb-1))
              DO 500 k=kmin,kmax
                call odd (lc+lr,i5)
                call odd (ln+lb,i6)
                if ((i5.eq.0.and.i6.ne.0).or.
     *              (i5.ne.0.and.i6.eq.0)) goto 900

***********************************************

              zf=((-1)**(k+kapr+kapc))/(2.d0*k+1)
              do 171 nr=nnr,nmaxx(indr)
                do 172 n=nnn,nmaxx(indn)
                  y(n,nr)=0.d0
* >>>>>>>>>> SIGMA >>>>>>>
                  dy(n,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue

              if (j.ge.ii) then
               index2=ipxc(j,ii,indn,indr,k)
              else
               index2=ipxc(ii,j,indr,indn,k)
              endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 n=nnn,nmaxx(indn)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>
                  eu=1.d0/(ec+eb-ee(n,indn)-ee(nr,indr)+evt-ev)
                  if (j.ge.ii) then
	             tt=xci(n,nr,index2)*eu
                   y(n,nr)=y(n,nr)+tt
                   dy(n,nr)=dy(n,nr)-tt*eu
                  else
	             tt=xci(nr,n,index2)*eu
                   y(n,nr)=y(n,nr)+tt
                   dy(n,nr)=dy(n,nr)-tt*eu
                 endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
10              continue
12            continue

              k1min=max0(iabs(kapn-kapc),iabs(kapr-kapb))
              k1max=min0((kapn+kapc-1),(kapr+kapb-1))
              DO 9 k1=k1min,k1max
                z2=(2*k+1)*d6j(jn,jb,2*k,jr,jc,2*k1)
                if (j.ge.ii) then
                 index2=ipxc(j,ii,indr,indn,k1)
                else
                 index2=ipxc(ii,j,indn,indr,k1)
                endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 n=nnn,nmaxx(indn)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>
                  eu=1.d0/(ec+eb-ee(n,indn)-ee(nr,indr)+evt-ev)
                    if (j.ge.ii) then
	               tt=z2*xci(nr,n,index2)*eu
                     y(n,nr)=y(n,nr)+tt
                     dy(n,nr)=dy(n,nr)-tt*eu
                    else
	               tt=z2*xci(n,nr,index2)*eu
                     y(n,nr)=y(n,nr)+tt
                     dy(n,nr)=dy(n,nr)-tt*eu
                    endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
101                continue
121              continue
9              continue
                 ka=kv
                 na=nv

                 CALL klj(ka,kapa,la,ja,inda,n0a)

                 do 23 indm=1,2*lmax+1
                   CALL indk1(indm,km,kapm,lm,jm,m0)
                   call st (km,nnm)
                   call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (ln+lb,i5)
                   call odd (lm+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700


                   do 165 m=nnm,nmaxx(indm)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,m)=0.d0
164                  continue
165                continue


                   call odd (lc+lr+k,i5)
                   call odd (lm+la+k,i6)
                   if (i5*i6.eq.0) goto 800
                   index1=ipxv1(ii,indm,indr,k)

                   z2=zf*((-1)**(kapc-kapr))
                   do 123 m=nnm,nmaxx(indm)
                     do 103 nr=nnr,nmaxx(indr)

                        x(nr,m)=x(nr,m)+z2*xv1(m,nr,index1)

103                  continue
123                continue
800                continue

                   zz=(-1)**(kapc+kapa+kapr+kapm)
                   k2min=max0(iabs(kapc-kapa),iabs(kapm-kapr))
                   k2max=min0((kapc+kapa-1),(kapm+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+la+k2,i1)
                     call odd (lm+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155
                      index1=ipxv2(ii,indr,indm,k2)
                     z1=zf*(2*k+1.d0)*d6j(ja,jm,2*k,jr,jc,2*k2)
                     do 124 m=nnm,nmaxx(indm)
                       do 104 nr=nnr,nmaxx(indr)

                          x(nr,m)=x(nr,m)+z1*zz*xv2(nr,m,index1)

104                    continue
124                  continue
155                  continue
119                continue
                   index=ipxv(j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
                        do 107 m=nnm,nmaxx(indm)
                         xvo(m,n,index)=xvo(m,n,index)+
     *                                  x(nr,m)*y(n,nr)
* >>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                         dxvo(m,n,index)=dxvo(m,n,index)+
     *                                  x(nr,m)*dy(n,nr)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
107                     continue
117                   continue
127                 continue
700              continue
23               continue
900           continue
500           continue
21          continue
22        continue
31      continue
32    continue
      return
      end

      subroutine term31
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION vv(nhf),u(nhf),v(nhf),fa(nhf),fb(nhf),ga(nhf),gb(nhf)
      DIMENSION v1(nhf),v2(nhf),uu(nhf,nx)
      do 32 j=1,ncore
        kb=ko(j)
        nb=no(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        DO 13 il=1,max
          gb(il)=go(il,j)
          fb(il)=fo(il,j)
13      CONTINUE
        do 23 indn=1,2*lmax+1
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (kn,nnn)
          kmin=iabs(kapn-kapb)
          kmax=kapn+kapb-1
          DO 500 k=kmin,kmax
            call odd (ln+lb+k,i5)
            if (i5.eq.0) goto 900
            do 4 n=nnn,nmaxx(indn)
              do 20 ij=1,max
                v(ij)=gb(ij)*gg(ij,n,indn)+fb(ij)*ff(ij,n,indn)
20            continue
              call yfun(v,u,k,max,*901)
              do 19 ij=1,max
                uu(ij,n)=u(ij)
19            continue
4           continue
              ka=kv
              na=nv
              CALL klj(ka,kapa,la,ja,inda,n0a)
              DO 15 il=1,max
                ga(il)=gv(il)
                fa(il)=fv(il)
15            CONTINUE
              do 26 indr=1,2*lmax+1

              if (indr.eq.inda) then
              CALL indk1(indr,kr,kapr,lr,jr,n0r)
              call st (kr,nnr)
                do 33 ij=1,max
                  v1(ij)=0.d0
                  v2(ij)=0.d0
33              continue
              do 10 nr=nnr,nmaxx(indr)
                do 11 ij=1,max
                  v1(ij)=v1(ij)+gg(ij,nr,indr)*rvi(nr)
                  v2(ij)=v2(ij)+ff(ij,nr,indr)*rvi(nr)
11              continue
10            continue
              do 25 indm=1,2*lmax+1
                CALL indk1(indm,km,kapm,lm,jm,m0)
                call st (km,nnm)
                call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                call odd (la+lm+k,i2)
                if (i1*i2.eq.0) goto 800
                c=((-1)**(k))*s(k,km,kr)*s(k,kn,kb)
                index=ipxv(j,indm,indn,k)
                do 165 n=nnn,nmaxx(indn)
                  do 160 m=nnm,nmaxx(indm)
                    do 161 ij=1,max
                      vv(ij)=(v1(ij)*gg(ij,m,indm)+
     *                        v2(ij)*ff(ij,m,indm))*uu(ij,n)*rp(ij)
161                 continue
             xvo(m,n,index)=xvo(m,n,index)+c*rint(vv,1,max,11,h)
160               continue
165             continue
800             continue
25            continue
             endif
26          continue
900         continue
500       continue
23      continue
32    continue
901   continue
      return
      end

      subroutine term32(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      DIMENSION dv1(nhf),dv2(nhf),dvv(nhf)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION vv(nhf),u(nhf),v(nhf),fa(nhf),fb(nhf),ga(nhf),gb(nhf)
      DIMENSION v1(nhf),v2(nhf),uu(nhf,nx)
        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 13 il=1,max
          ga(il)=gv(il)
          fa(il)=fv(il)
13      CONTINUE
        do 23 indm=1,2*lmax+1
          CALL indk1(indm,km,kapm,lm,jm,m0)
          call st (km,nnm)
          kmin=iabs(kapm-kapa)
          kmax=kapm+kapa-1
          DO 500 k=kmin,kmax
            call odd (lm+la+k,i5)
            if (i5.eq.0) goto 900
            do 4 m=nnm,nmaxx(indm)
              do 20 ij=1,max
                v(ij)=ga(ij)*gg(ij,m,indm)+fa(ij)*ff(ij,m,indm)
20            continue
              call yfun(v,u,k,max,*901)
              do 19 ij=1,max
                uu(ij,m)=u(ij)
19            continue
4           continue
            do 30 j=1,ncore
              kb=ko(j)
              nb=no(j)
              CALL klj(kb,kapb,lb,jb,indb,n0b)
              DO 15 il=1,max
                gb(il)=go(il,j)
                fb(il)=fo(il,j)
15            CONTINUE
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>
               eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          do 26 indr=1,2*lmax+1

              if (indr.eq.indb) then
              CALL indk1(indr,kr,kapr,lr,jr,n0r)
              call st (kr,nnr)

                  do 33 ij=1,max
                  v1(ij)=0.d0
                  v2(ij)=0.d0
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  dv1(ij)=0.d0
                  dv2(ij)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
33              continue
              do 10 nr=nnr,nmaxx(indr)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                 eu=1.d0/(eb-ee(nr,indr)+evt-ev)

                do 11 ij=1,max

	            tt1=gg(ij,nr,indr)*rmi(nr,j)*eu
	            tt2=ff(ij,nr,indr)*rmi(nr,j)*eu

                  v1(ij)=v1(ij)+tt1
                  v2(ij)=v2(ij)+tt2

                  dv1(ij)=dv1(ij)-tt1*eu
                  dv2(ij)=dv2(ij)-tt2*eu

11              continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

10            continue
              do 25 indn=1,2*lmax+1
                CALL indk1(indn,kn,kapn,ln,jn,n0)
                call st (kn,nnn)
                call trgi(iabs(kapn-kapr),(kapn+kapr-1),k,i1)
                call odd (lr+ln+k,i2)
                if (i1*i2.eq.0) goto 800
                c=((-1)**(k))*s(k,kn,kr)*s(k,km,ka)
                index=ipxv(j,indm,indn,k)
                do 165 n=nnn,nmaxx(indn)
                  do 160 m=nnm,nmaxx(indm)
                    do 161 ij=1,max
                      vv(ij)=(v1(ij)*gg(ij,n,indn)+
     *                        v2(ij)*ff(ij,n,indn))*uu(ij,m)*rp(ij)
* >>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                      dvv(ij)=(dv1(ij)*gg(ij,n,indn)+
     *                         dv2(ij)*ff(ij,n,indn))*uu(ij,m)*rp(ij)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
161                 continue
             xvo(m,n,index)=xvo(m,n,index)+c*rint(vv,1,max,11,h)
* >>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
             dxvo(m,n,index)=dxvo(m,n,index)+c*rint(dvv,1,max,11,h)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
160               continue
165             continue
800             continue
25            continue
             endif
26          continue
30          continue
900         continue
500       continue
23      continue
901   continue
      return
      end

      subroutine terma4(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION vv(nhf),u(nhf),v(nhf),fa(nhf),fb(nhf),ga(nhf),gb(nhf)
      DIMENSION v1(nhf),v2(nhf),uu(nhf,nx),x(nx,nx)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DIMENSION dx(nx,nx),dv1(nhf),dv2(nhf),dvv(nhf)
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /drhovin/  dxvi(nxx,nxx,nvh),drvi(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 32 j=1,ncore
        kb=ko(j)
        nb=no(j)
* >>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        DO 13 il=1,max
          gb(il)=go(il,j)
          fb(il)=fo(il,j)
13      CONTINUE
        do 23 indn=1,2*lmax+1
          CALL indk1(indn,kn,kapn,ln,jn,n0)
          call st (kn,nnn)
          kmin=iabs(kapn-kapb)
          kmax=kapn+kapb-1
          DO 500 k=kmin,kmax
            call odd (ln+lb+k,i5)
            if (i5.eq.0) goto 900
            do 4 n=nnn,nmaxx(indn)
              do 20 ij=1,max
                v(ij)=gb(ij)*gg(ij,n,indn)+fb(ij)*ff(ij,n,indn)
20            continue
              call yfun(v,u,k,max,*901)
              do 19 ij=1,max
                uu(ij,n)=u(ij)
19            continue
4           continue
            ka=kv
            na=nv
            CALL klj(ka,kapa,la,ja,inda,n0a)
            DO 15 il=1,max
              ga(il)=gv(il)
              fa(il)=fv(il)
15          CONTINUE
            do 25 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,m0)
              call st (km,nnm)

              if (indm.eq.inda) then
              CALL indk1(indm,km,kapm,lm,jm,n0m)
              call st (km,nnm)
              do 24 indr=1,2*lmax+1
                CALL indk1(indr,kr,kapr,lr,jr,n0r)
                call st (kr,nnr)
                call trgi(iabs(kapm-kapr),(kapm+kapr-1),k,i1)
                call odd (lr+lm+k,i2)
                if (i1*i2.eq.0) goto 800
                 index2=ipxv(j,indr,indn,k)

            zz=((-1)**(kapa+kapr+kapb+kapn))/((ja+1.d0)*(2*k+1.d0))
                do 165 n=nnn,nmaxx(indn)
                  do 160 nr=nnr,nmaxx(indr)

* >>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>
                    eu=1.d0/(evt+eb-ee(nr,indr)-ee(n,indn))
                    x(nr,n)=xvi(nr,n,index2)*eu
                    dx(nr,n)=-x(nr,n)*eu+dxvi(nr,n,index2)*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
160               continue
165             continue
                k2min=max0(iabs(kapr-kapb),iabs(kapn-kapa))
                k2max=min0((kapr+kapb-1),(kapn+kapa-1))
                DO 119 k2=k2min,k2max
                  call odd (lr+lb,i5)
                  call odd (ln+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 155


                   index2=ipxv(j,indn,indr,k2)

                  z1=(2*k+1.d0)*d6j(jr,ja,2*k,jn,jb,2*k2)
                  do 185 n=nnn,nmaxx(indn)
                    do 180 nr=nnr,nmaxx(indr)
* >>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>
                      eu=1.d0/(evt+eb-ee(nr,indr)-ee(n,indn))
                      tt=z1*xvi(n,nr,index2)*eu
                      x(nr,n)=x(nr,n)+tt
                      dx(nr,n)=dx(nr,n)-tt*eu+z1*dxvi(n,nr,index2)*eu
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
180                 continue
185               continue
155               continue
119               continue
                  c=((-1)**(k))*s(k,km,kr)*s(k,kb,kn)

                  do 33 ij=1,max
                    v1(ij)=0.d0
                    v2(ij)=0.d0
* >>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>
                    dv1(ij)=0.d0
                    dv2(ij)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
33                continue
                 do 109 n=nnn,nmaxx(indn)
                   do 10 nr=nnr,nmaxx(indr)
                     do 11 ij=1,max
                     v1(ij)=v1(ij)+gg(ij,nr,indr)*x(nr,n)*uu(ij,n)
                     v2(ij)=v2(ij)+ff(ij,nr,indr)*x(nr,n)*uu(ij,n)
* >>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                     dv1(ij)=dv1(ij)+gg(ij,nr,indr)*dx(nr,n)*uu(ij,n)
                     dv2(ij)=dv2(ij)+ff(ij,nr,indr)*dx(nr,n)*uu(ij,n)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
11                   continue
10                 continue
109              continue
                 do 169 m=nnm,nmaxx(indm)
                   do 161 ij=1,max
                     vv(ij)=(v1(ij)*gg(ij,m,indm)+
     *                       v2(ij)*ff(ij,m,indm))*rp(ij)
* >>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                     dvv(ij)=(dv1(ij)*gg(ij,m,indm)+
     *                        dv2(ij)*ff(ij,m,indm))*rp(ij)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
161                continue
                   rvo(m)=rvo(m)+c*zz*rint(vv,1,max,11,h)
* >>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                   drvo(m)=drvo(m)+c*zz*rint(dvv,1,max,11,h)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

169              continue
800              continue
24             continue
           endif
25         continue
900          continue
500        continue
23      continue
32    continue
901   continue
      return
      end

      subroutine term2(evt)
      implicit double precision (a-h,o-z)
      include "all.par"
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhov/ ipxv(nn,nk,nk,0:kk),ivchan
      COMMON /rhovin/  xvi(nxx,nxx,nvh),rvi(nx)
      COMMON /rhovout/ xvo(nxx,nxx,nvh),rvo(nx)
      common /ipv/ ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH)
      COMMON /val/ gv(nhf),fv(nhf),ev,nv,kv
      DIMENSION v(nhf),u(nhf),uu(nhf,nx,nx),v1(nhf,nx),v2(nhf,nx),
     *vv(nhf),uv1(nhf,nx),uv2(nhf,nx)
      DIMENSION uu1(nhf,nx,0:kk,nk),uu2(nhf,nx,0:kk,nk)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhovout/ dxvo(nxx,nxx,nvh),drvo(nx)
      COMMON /drhovin/  dxvi(nxx,nxx,nvh),drvi(nx)
      DIMENSION dv1(nhf,nx),dv2(nhf,nx),duv1(nhf,nx),duv2(nhf,nx)
      DIMENSION duu1(nhf,nx,0:kk,nk),duu2(nhf,nx,0:kk,nk),dvv(nhf)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ka=kv
        na=nv
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 6 j=1,ncore
          kb=ko(j)
          nb=no(j)
* >>>>>>>>> SIGMA >>>>>>>>>>>>>>
          eb=eo(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kb,kapb,lb,jb,indb,n0b)
c          t1=mclock()
          do 1 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 91 indn=1,nk
              DO 92 k=0,kk
                do 93 m=1,nx
                  do 94 ij=1,nhf
                    uu1(ij,m,k,indn)=0.d0
                    uu2(ij,m,k,indn)=0.d0
* >>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>
                    duu1(ij,m,k,indn)=0.d0
                    duu2(ij,m,k,indn)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
94                continue
93              continue
92            continue
91          continue

            do 2 indr=1,2*lmax+1
              CALL indk1(indr,kr,kapr,lr,jr,n0r)
              call st (kr,nnr)
              l1min=iabs(kapm-kapr)
              l1max=kapm+kapr-1
              DO 500 l=l1min,l1max
                call odd (lm+lr+l,i5)
                if (i5.eq.0) goto 900
                 do 3 m=nnm,nmaxx(indm)
                   do 4 nr=nnr,nmaxx(indr)
                     do 20 ij=1,max
                       v(ij)=gg(ij,nr,indr)*gg(ij,m,indm)+
     *                       ff(ij,nr,indr)*ff(ij,m,indm)
20                   continue
                     call yfun(v,u,l,max,*901)
                     do 19 ij=1,max
                       uu(ij,nr,m)=u(ij)
19                     continue
4                  continue
3                continue
                 do 7 inds=1,2*lmax+1
                   CALL indk1(inds,ks,kaps,ls,js,n0s)
                   call st (ks,nns)
                   call odd (lr+la,i5)
                   call odd (ls+lb,i6)
                    if ((i5.eq.0.and.i6.ne.0).or.
     *                  (i5.ne.0.and.i6.eq.0)) goto 155
                    k1min=max0(iabs(kapr-kapa),iabs(kaps-kapb))
                    k1max=min0((kapr+kapa-1),(kaps+kapb-1))
                    DO 119 k1=k1min,k1max
                      index1=ipxv(j,indr,inds,k1)
                      do 31 nr=nnr,nmaxx(indr)
                        do 203 ij=1,max
                          v1(ij,nr)=0.d0
                          v2(ij,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>
                          dv1(ij,nr)=0.d0
                          dv2(ij,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
203                     continue
                        do 41 nss=nns,nmaxx(inds)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          eu=1.d0/(evt+eb-ee(nr,indr)-ee(nss,inds))

                          do 201 ij=1,max

	    tt1=gg(ij,nss,inds)*xvi(nr,nss,index1)*eu
	    tt2=ff(ij,nss,inds)*xvi(nr,nss,index1)*eu

	    dtt1=gg(ij,nss,inds)*dxvi(nr,nss,index1)*eu
	    dtt2=ff(ij,nss,inds)*dxvi(nr,nss,index1)*eu

          v1(ij,nr)=v1(ij,nr)+tt1
          v2(ij,nr)=v2(ij,nr)+tt2

          dv1(ij,nr)=dv1(ij,nr)-tt1*eu+dtt1
          dv2(ij,nr)=dv2(ij,nr)-tt2*eu+dtt2

201                       continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

41                      continue
31                    continue

                      do 51 m=nnm,nmaxx(indm)
                        do 206 ij=1,max
                          uv1(ij,m)=0.d0
                          uv2(ij,m)=0.d0
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                          duv1(ij,m)=0.d0
                          duv2(ij,m)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
206                     continue
                        do 57 nr=nnr,nmaxx(indr)
                          do 202 ij=1,max
                            uv1(ij,m)=uv1(ij,m)+v1(ij,nr)*uu(ij,nr,m)
                            uv2(ij,m)=uv2(ij,m)+v2(ij,nr)*uu(ij,nr,m)
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                          duv1(ij,m)=duv1(ij,m)+dv1(ij,nr)*uu(ij,nr,m)
                          duv2(ij,m)=duv2(ij,m)+dv2(ij,nr)*uu(ij,nr,m)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
202                       continue
57                      continue
51                    continue

                      do 188 indn=1,2*lmax+1
                        CALL indk1(indn,kn,kapn,ln,jn,n0)
                        call st (kn,nnn)
                        call odd (la+lm,i5)
                        call odd (ln+lb,i6)
                        if ((i5.eq.0.and.i6.ne.0).or.
     *                      (i5.ne.0.and.i6.eq.0)) goto 800

                        call trgi(iabs(kapn-kaps),(kapn+kaps-1),l,i1)
                        call odd (ln+ls+l,i2)
                        if (i1*i2.eq.0) goto 800
                        c=((-1)**(l))*s(l,km,kr)*s(l,kn,ks)

                        z1=((-1)**(kapa+kapb+kapn+kapm))
               kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb),iabs(k1-l))
               kmax=min0((kapm+kapa-1),(kapn+kapb-1),(k1+l))
      if (kmax.gt.kk) then 
	  write (*,*) ' Parameter KK is exceeded in term 2, STOP'
	stop
      endif
		   

                        DO 219 k=kmin,kmax
                          c1=(2*k+1.d0)*d6j(jm,ja,2*k,2*k1,2*l,jr)
                          cfin=c*c1*z1*d6j(jn,jb,2*k,2*k1,2*l,js)
                          do 59 m=nnm,nmaxx(indm)
                            do 280 ij=1,max
               uu1(ij,m,k,indn)=uu1(ij,m,k,indn)+uv1(ij,m)*cfin
               uu2(ij,m,k,indn)=uu2(ij,m,k,indn)+uv2(ij,m)*cfin
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
               duu1(ij,m,k,indn)=duu1(ij,m,k,indn)+duv1(ij,m)*cfin
               duu2(ij,m,k,indn)=duu2(ij,m,k,indn)+duv2(ij,m)*cfin
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
280                         continue
59                        continue
219                     continue
800                     continue
188                   continue
119                 continue
155                 continue
7                 continue
900               continue
500             continue
2             continue
              do 189 indn=1,2*lmax+1
                CALL indk1(indn,kn,kapn,ln,jn,n0)
                call st (kn,nnn)
                call odd (la+lm,i5)
                call odd (ln+lb,i6)
                if ((i5.eq.0.and.i6.ne.0).or.
     *              (i5.ne.0.and.i6.eq.0)) goto 810

                kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
                kmax=min0((kapm+kapa-1),(kapn+kapb-1))
                DO 319 k=kmin,kmax
                  index=ipxv(j,indm,indn,k)
                  do 52 n=nnn,nmaxx(indn)
                    do 53 m=nnm,nmaxx(indm)
                      do 205 ij=1,max
                        vv(ij)=(uu1(ij,m,k,indn)*gg(ij,n,indn)+
     *                          uu2(ij,m,k,indn)*ff(ij,n,indn))*rp(ij)
* >>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        dvv(ij)=(duu1(ij,m,k,indn)*gg(ij,n,indn)+
     *                           duu2(ij,m,k,indn)*ff(ij,n,indn))*rp(ij)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
205                   continue
                      xvo(m,n,index)=xvo(m,n,index)+
     *                rint(vv,1,max,11,h)
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                      dxvo(m,n,index)=dxvo(m,n,index)+
     *                rint(dvv,1,max,11,h)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
53                  continue
52                continue
319             continue
810             continue
189           continue
1           continue
6        continue
901   continue
      return
      end
c     ========================================
      include "inidat.f"
      include "d6j.f"
      include "libD.f"
      include "rint.f"
      include "yfun.f"
      include "yint.f"
      include "linpakd.f"

