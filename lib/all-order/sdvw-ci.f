* * April 20, 2012 Read in is modified to match core and valence files
*     NEW utterly modified version 
*     Version 3.0 May 7, 2008
*     rho equations are replaced by sigma's equations,
*     derivatives are added 
**********************************************************
*
*      Version 2.0 May 6, 2008
*      Renamed and taken parameters out, added counting 
*      soubroutines and parameter checks
*                                         
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
* 10     # max number of iterations                                     *
* 1	   # kval (key_en): key for energies,                             *
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
	character*2 lab1,lab2,lab3,lab4
	include "allvw.par"
      PARAMETER (n9=1000)
      DIMENSION t(nkx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      DIMENSION g(nhf,nx),f(nhf,nx),e(nl)
      DIMENSION wco(ns),nco(ns),kco(ns),mco(ns)
      COMMON /x12ini/ inmax
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      dimension n19(n9),k19(n9),n29(n9),k29(n9),n39(n9),k39(n9),n49(n9)
      dimension k49(n9),k9(n9)
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
      
	call inidat 
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
      read (5,*) nitmax
*>>>>>>> RLE >>>>>>>>>>>>>>>>>>>>
      read (5,*) ikey,max_rle,i999
      read (5,*) damp
* It is added to use the same input file as for the core and valence
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

c>>>ms: added 10/04/09 
      do i=1,2*lmax+1
        if (nmaxx(i).GT.nmax) nmaxx(i)=nmax
        write(*,*) ' PW=',i,' nmax=',nmaxx(i)
      end do
c<<<ms

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
*      The kval input below is ignored by the code and is added
*      only for consistency
* >>>>>>>> SIGMA >>>>>>>>>
      read (5,*) kval
     

	if (kval.eq.2) then
	 read (5,*) lvmax 
	 read (5,*) j,ecorr
	 i=2
	  if (lvmax.ne.0) then  
	   do 78 j=1,lvmax
	     read  (5,*) ii,ecorr,ecorr1
	     i=i+2
78       continue	    
	  endif
      endif
      if (kval.lt.0.or.kval.gt.2) then 
	 write (*,*) 'WRONG INPUT IN KVAL, STOP'
	endif
* >>>>>>>>>>>>>>>>>>>>>>>>
****************************************************************
        read (5,*) nval
      if (nval.gt.nz) then 
	  write (*,*) ' Parameter NZ is exceeded, STOP'
	stop
      endif


	do 10 i=1,nval
        READ (5,*) nv(i),kv(i),i99
c        write (*,*) nv(i),kv(i)
        kvv=kv(i)
	  nvv=nv(i)
        CALL klj (kvv,kapv,lv,jv,indv,n0v)
        DO 29 j=1,max
          fv(j,i)=ff(j,nvv-n0v,indv)
          gv(j,i)=gg(j,nvv-n0v,indv)
29      CONTINUE
        ev(i)=ee(nvv-n0v,indv)
        write (*,'(2i3,f12.7)') nv(i),kv(i),ev(i)
10    continue
c       if (indv.gt.inmax) then
c       inmax=indv
c       endif
******************************************************
      write (*,*) nmax,lmax,nmax1,lmax1
c	open (unit=8,file='in.txt')
* If calculating radial integrals, uncomment next line
c	call xvalw
      goto 9888
      read (8,*) num
	do 1387 i=1,num
	 read (8,*) n19(i),k19(i),n29(i),k29(i),n39(i),k39(i),n49(i),
     * k49(i),k9(i)
       n1=n19(i)
       k1=k19(i)
       n2=n29(i)
       k2=k29(i)
       n3=n39(i)
       k3=k39(i)
       n4=n49(i)
       k4=k49(i)
       k=k9(i)
	call llabel(k1,lab1)
	call llabel(k2,lab2)
	call llabel(k3,lab3)
	call llabel(k4,lab4)

c      call search_rad(n1,k1,n2,k2,n3,k3,n4,k4,k,r1)
c      call search_rad(n2,k2,n1,k1,n4,k4,n3,k3,k,r2)
c      call search_rad(n3,k3,n4,k4,n1,k1,n2,k2,k,r3)

c	write (*,9) n1,lab1,n2,lab2,n3,lab3,n4,lab4,k,r1,r2,r3
1387  continue
c	 call outrad
c      stop
*     Counting all channels
9888   continue
      call COUNT_VW
	call COUNT1
	call COUNT_xfor4

	call rhocore 

      call xcnvw
	call rhoval(kval)
	call xfor4
      call ini
*********************************************

	call term1
      call term51
	call term52

      call term41a
c      call term41b
      call term41c
      call term41d

      call term42a
c      call term42b
      call term42c
      call term42d
9     format (4(i3,a2),i5,3f15.10)
      goto 877
	do 17 i=1,num

       n1=n19(i)
       k1=k19(i)
       n2=n29(i)
       k2=k29(i)
       n3=n39(i)
       k3=k39(i)
       n4=n49(i)
       k4=k49(i)
       k=k9(i)


	call llabel(k1,lab1)
	call llabel(k2,lab2)
	call llabel(k3,lab3)
	call llabel(k4,lab4)

      call search_rad(n1,k1,n2,k2,n3,k3,n4,k4,k,r1,dr1)
      call search_rad(n2,k2,n1,k1,n4,k4,n3,k3,k,r2,dr2)
      call search_rad(n3,k3,n4,k4,n1,k1,n2,k2,k,r3,dr3)

	write (*,9) n1,lab1,n2,lab2,n3,lab3,n4,lab4,k,r1,r2,r3
	write (*,9) n1,lab1,n2,lab2,n3,lab3,n4,lab4,k,dr1,dr2,dr3
      write (*,*) 
17     continue      
877   continue
	 call outrad(kval)
*******************************************************	  
      stop
 710  write(*,*) ' No input file "inf.aov" found '
      read(*,*)
      stop
 720  write(*,*) ' No file "hfspl.1" ifound '
      stop
 730  write(*,*) ' No file "hfspl.2" found '
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

      subroutine COUNT1
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
**************************************
* COUNTING X_k(civw) channels, NXCN  *
**************************************
      index=0
      do 14 i=1,nval
        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
*****************************************************
        do 3 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 2 indn=1,2*lmax+1
            CALL indk1(indn,kn,kapn,ln,jn,n0)
            call st (kn,nnn)
****************************************************************
*             COUNT  Xk(cnvw) m can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapc-kapa),iabs(kapn-kapb))
              kmax=min0((kapc+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lc+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                index=index+1
	         if (index.ge.nxcn) then
c	          write (*,*) 'NUMBER OF CNVW CHANNELS IS EXCEEDED, STOP'
c	          stop
	         endif
555             continue
11            continue
*****************************************************************

2           continue
3         continue
4       continue
14    continue
***********************************************************************      
	if (index.le.NXCN) then 
	 write (*,122) index,nxcn
      else
	 write (*,123) index,nxcn
       write (*,*) 'Parameter NXCN is exceeded, STOP'
	 stop
	endif
122   format ('Number of X_k(civw) channels      = ',i9,' < NXCN = ',i9,
     *' OK')
123   format ('Number of X_k(civw) channels      = ',i9,' > NXCN = ',i9)
   
********************************************************************
      return
      end



****************************************************************************************
      subroutine outrad(kval)
      implicit double precision (a-h,o-z)
      include "allvw.par"
     
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      common /nmaxx/nmaxx(nx)
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      OPEN (unit=17,file='sigma1')

      OPEN (unit=7,file='pair.vw',form='unformatted')
      write (7) ichan
      write (7) nmaxx
	write (7) ipm,ipn
      write (7) nv,kv,nval
      write (7) ipvw
* >>>> SIGMA >>>>>>>>>>
      write (7) kval,evt
* >>>>>>>>>>>>>>>>>>>>>
      do 27 i=1,ichan
        indn=ipn(i)
        indm=ipm(i)
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)
        do 37 n=1,nmaxn
          write (7) (xvw(m,n,i),m=1,nmaxm)
37       continue
* >>>> SIGMA >>>>>>>>>>
        do 327 n=1,nmaxn
          write (7) (dxvw(m,n,i),m=1,nmaxm)
327       continue
* >>>>>>>>>>>>>>>>>>>>>

27     continue
       close (7)
      do 1 i=1,nval

	 nvv=nv(i)
	 kvv=kv(i)
         CALL klj(kvv,kapv,lv,jv,indv,n0v)
        indm=indv
        CALL indk1(indm,km,kapm,lm,jm,m0)
        call st (km,nnm)
        do 200 m=nnm,nmaxx(indm)
	write (17,'(3i4,2e20.8)') nvv,kvv, m+m0,rsigma(i,m),drsigma(i,m)  
200     continue
1     continue
      close(17)
      return
      end
c


      subroutine search_rad(m,km,n,kn,nvv,kvv,nw,kw,k,r,dr)
      implicit double precision (a-h,o-z)
      include "allvw.par"

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
* >>>>>>>SIGMA 
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
* >>>>>>>
      call sv(nvv,kvv,iv)
      call sv(nw,kw,iw)
 
      CALL klj(kvv,kapv,lv,jv,indv,n0v)
      CALL klj(kw,kapw,lw,jw,indw,n0w)
      CALL klj(km,kapm,lm,jm,indm,m0)
      CALL klj(kn,kapn,ln,jn,indn,n0)

      index=ipvw(iv,iw,indm,indn,k)

***************************************************************
	if (index.eq.0) then 
	 write (*,*) 'NO SUCH INTERGAL',nvv,kvv,nw,kw,m,km,n,kn,k
	 stop
      endif
***************************************************************

	m9=m-m0
	n9=n-n0
	r=xvw(m9,n9,index)
	dr=dxvw(m9,n9,index)

      return
	end

*****************************************************************************************
      subroutine sv(nvv,kvv,ival)
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      do 1 i=1,nval
	 if (nv(i).eq.nvv.and.kv(i).eq.kvv) then 
	  ival=i
	goto 800
	endif
1     continue
800   continue
	return
	end

      subroutine COUNT_VW
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /x12ini/ inmax
      index=0
      do 14 i=1,nval
        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
*****************************************************

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnvw) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 559
                 
			 index=index+1
559    continue
11            continue
*****************************************************************

5           continue
2         continue
4       continue
14    continue
***********************************************************************      
      write  (*,*)
	write (*,*) 'PRE-COUNTING ALL CHANNELS'
	write (*,*)
	if (index.le.NVWH) then 
	 write (*,122) index,nvwh
      else
	 write (*,123) index,nvwh
       write (*,*) 'Parameter NVWH is exceeded, STOP'
	 stop
	endif
122   format ('Number of Sigma2_k(mnvw) channels = ',i9,' < NVWH = ',
     *i9,' OK')
123   format ('Number of Sigma2_k(mnvw) channels = ',i9,' > NVWH = ',
     *i9)

   
********************************************************************

      return
      end
*****************************************************************************************
      subroutine xvalw
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /x12ini/ inmax
      icount=0
      index=0
      do 531 i=0,NVWH
        do 532 m=1,nxx
          do 533 n=1,nxx
            xvw(m,n,i)=0.d0
533       continue
532     continue
531   continue

      do 661 i5=0,10
       do 662 i4=1,nk
         do 663 i3=1,nk
         do 664 i2=1,nz
          do 665 i1=1,nz
                ipvw(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

      do 14 i=1,nval
      DO 12 il=1,max
          ga(il)=gv(il,i)
          fa(il)=fv(il,i)
12      CONTINUE
        ea=ev(i)
        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          DO 13 il=1,max
            gb(il)=gv(il,j)
            fb(il)=fv(il,j)
13        CONTINUE
          eb=ev(j)
*****************************************************

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
		   
c            do 5 indn=indm,2*lmax+1
* Changed for consistency with ini
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnvw) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 559
                 
			 index=index+1
			 if (index.ge.NVWH) then
	          write (*,*) 'NUMBER OF CHANNELS IS EXCEEDED, STOP'
	          stop
	         endif
                ipvw(i,j,indm,indn,k)=index
	 ipm(index)=indm
	 ipn(index)=indn
	 ippv(index)=i
	 ippw(index)=j
	 ipk(index)=k	     

                do 66 m=1,nmaxx(indm)
                  do 677 n=1,nmaxx(indn)
                    xvw(m,n,index)=0.d0
677                 continue
66               continue

                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555

c                c=((-1)**(k))*s(k,km,ka)*s(k,kn,kb)
                 c=1.d0
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
                    icount=icount+1
                    xvw(m,n,index)=c*rint(vv,1,max,11,h)
6                 continue
3               continue
555             continue
559    continue
11            continue
*****************************************************************

5           continue
c         endif
2         continue
187       format ('CHANNEL',i6,'  j2=',i3,'  km=',i3,
     *    '  kn=',i3,'  k=',i3,f12.6)
4       continue
14    continue
c      write (*,*) index
      ichan=index
901   continue
      return
      end
*****************************************************************************************
      subroutine ini
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
* >>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /x12ini/ inmax
      icount=0
      index=0
      do 531 i=0,NVWH
        do 532 m=1,nxx
          do 533 n=1,nxx
            xvw(m,n,i)=0.d0
* >>>>>>>>> SIGMA >>>>>>>>>>
            dxvw(m,n,i)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>
533       continue
532     continue
531   continue

      do 661 i5=0,10
       do 662 i4=1,nk
         do 663 i3=1,nk
         do 664 i2=1,nz
          do 665 i1=1,nz
                ipvw(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

      do 14 i=1,nval
        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
*****************************************************

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
c            if (inda.ge.indm) then 
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnvw) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax

                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 555

                index=index+1
                ipvw(i,j,indm,indn,k)=index
	          ipm(index)=indm
	          ipn(index)=indn
	          ippv(index)=i
	          ippw(index)=j
	          ipk(index)=k

	         if (index.ge.NVWH) then
	          write (*,*) 'NUMBER OF CHANNELS IS EXCEEDED, STOP'
	          stop
	         endif
                do 3 m=1,nmaxx(indm)
                  do 6 n=1,nmaxx(indn)
                    icount=icount+1
                    xvw(m,n,index)=0.d0
* >>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>
                    dxvw(m,n,index)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
6                 continue
3               continue
555             continue
11            continue
*****************************************************************

5           continue
c         endif
2         continue
4       continue
14    continue

      ichan=index
901   continue
      return
      end


*****************************************************************************************
      subroutine xcnvw
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radcbvw/ ipcn(nz,nz,nn,nk,0:kk),xcn(nxx,0:nxcn),ichancn
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      DIMENSION gb(nhf),fb(nhf),ga(nhf),fa(nhf),gc(nhf),fc(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0
      do 531 i=0,nxcn
        do 532 m=1,nxx
            xcn(m,i)=0.d0
532     continue
531   continue

      do 661 i5=0,kk
       do 662 i4=1,nk
         do 663 i3=1,nn
         do 664 i2=1,nz
          do 665 i1=1,nz
                ipcn(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

      do 14 i=1,nval
      DO 12 il=1,max
          ga(il)=gv(il,i)
          fa(il)=fv(il,i)
12      CONTINUE
        ka=kv(i)
        na=nv(i)
c	  write (*,*) i,index
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          DO 13 il=1,max
            gb(il)=gv(il,j)
            fb(il)=fv(il,j)
13        CONTINUE
*****************************************************
        do 3 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          DO 133 il=1,max
            gc(il)=go(il,ii)
            fc(il)=fo(il,ii)
133        CONTINUE
          

          do 2 indn=1,2*lmax+1
            CALL indk1(indn,kn,kapn,ln,jn,n0)
            call st (kn,nnn)

****************************************************************
*             Calculate  Xk(cnvw) m can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapc-kapa),iabs(kapn-kapb))
              kmax=min0((kapc+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lc+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                index=index+1
                ipcn(i,j,ii,indn,k)=index
	         if (index.ge.nxcn) then
	          write (*,*) 'NUMBER OF CNVW CHANNELS IS EXCEEDED, STOP'
	          stop
	         endif
                c=((-1)**(k))*s(k,kc,ka)*s(k,kn,kb)
                  do 19 ij=1,max
                    v(ij)=ga(ij)*gc(ij)+fa(ij)*fc(ij)
19                continue
                  call yfun(v,u,k,max,*901)
                  do 6 n=1,nmaxx(indn)
                    do 20 ij=1,max
                      vv(ij)=(gb(ij)*gg(ij,n,indn)+
     *                        fb(ij)*ff(ij,n,indn))*rp(ij)*u(ij)
20                  continue
                    xcn(n,index)=c*rint(vv,1,max,11,h)
6                 continue
555             continue
11            continue
*****************************************************************

2           continue
3         continue
4       continue
14    continue
c      write (*,*) 'indexcnvw=',index
      ichancn=index
c	write (*,*) icount
901   continue
      return
      end



      subroutine rhocore
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
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


c      if (llmax.lt.lmax.or.nnmax.lt.nmax) then
c       write (*,*) 'WRONG LMAX OR NMAX IN THE CORE FILE'
c       stop
c      endif
901   continue
      return
      end

      subroutine term1
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radcbvw/ ipcn(nz,nz,nn,nk,0:kk),xcn(nxx,0:nxcn),ichancn
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 30 index=1,ichan
        i=ippv(index)
        j=ippw(index)

* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        indm=ipm(index)
        indn=ipn(index)
        k=ipk(index)

        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        kb=kv(j)
        nb=nv(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)

        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0

        do 3 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 4 jj=1,ncore
            kd=ko(jj)
            nd=no(jj)
            CALL klj(kd,kapd,ld,jd,indd,n0d)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            ec=eo(ii)
		  ed=eo(jj)    
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            l1min=max0(iabs(kapc-kapa),iabs(kapd-kapb))
            l1max=min0((kapc+kapa-1),(kapd+kapb-1))
            DO 5 l=l1min,l1max
              call odd (la+lc+l,i1)
              call odd (lb+ld+l,i2)
              if (i1*i2.eq.0) goto 515
              index1=ipcn(i,j,ii,indd,l)
	
	if (index1.eq.0) then
	 write (*,*) 'Term 1 - missing Xcnvw channel, STOP'
	 stop
	endif
              k1min=max0(iabs(kapm-kapc),iabs(kapn-kapd),iabs(k-l))
              k1max=min0((kapm+kapc-1),(kapn+kapd-1),(k+l))
              DO 9 k1=k1min,k1max
                call odd (lm+lc,i1)
                call odd (ln+ld,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 535
                z1=(-1)**(kapm+kapn+kapa+kapb)
                z2=d6j(jm,ja,2*k,2*l,2*k1,jc)
                z3=cq*z1*z2*(2.d0*k+1.d0)*d6j(jn,jb,2*k,2*l,2*k1,jd)
                if (ii.ge.jj) then
                 index2=ipxc(ii,jj,indm,indn,k1)
                else
                 index2=ipxc(jj,ii,indn,indm,k1)
                endif
                do 10 n=nnn,nmaxx(indn)
                  do 11 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
           eu=1.d0/(ec+ed-ee(m,indm)-ee(n,indn)+eevt-eev+eewt-eew)
                    if (ii.ge.jj) then
	  tt=z3*xcn(nd-n0d,index1)*xci(m,n,index2)*eu

         xvw(m,n,index)=xvw(m,n,index)+tt
         dxvw(m,n,index)=dxvw(m,n,index)-tt*eu
                    else
         tt=z3*xcn(nd-n0d,index1)*xci(n,m,index2)*eu

         xvw(m,n,index)=xvw(m,n,index)+tt
         dxvw(m,n,index)=dxvw(m,n,index)-tt*eu

                    endif
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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



      subroutine term51
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radcbvw/ ipcn(nz,nz,nn,nk,0:kk),xcn(nxx,0:nxcn),ichancn
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 30 index=1,ichan
        i=ippv(index)
        j=ippw(index)
        indm=ipm(index)
        indn=ipn(index)
        k=ipk(index)

        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        kb=kv(j)
        nb=nv(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0

        do 40 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>
          ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          call odd (la+lc+k,i1)
          call odd (lb+ln+k,i2)
          if (kc.eq.km.and.(i1*i2).ne.0) then
           iu=ipcn(i,j,ii,indn,k)
	if (iu.eq.0) then
	 write (*,*) 'Term 51 - missing Xcnvw channel, STOP'
	 stop
	endif

           do 4 n=nnn,nmaxx(indn)
             do 41 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                eu=1.d0/(ec-ee(m,indm)+eevt-eev+eewt-eew)
                tt = -cq*xcn(n,iu)*rmi(m,ii)*eu

                xvw(m,n,index)=xvw(m,n,index)+tt
                dxvw(m,n,index)=dxvw(m,n,index)-tt*eu

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
41           continue
4          continue
          endif
40      continue
30    continue
901   continue
      return
      end

      subroutine term52
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /radcbvw/ ipcn(nz,nz,nn,nk,0:kk),xcn(nxx,0:nxcn),ichancn
      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
	COMMON /add/ ipm(NVWH),ipn(NVWH),ippv(NVWH),ippw(NVWH),ipk(NVWH)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nx,nn)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do 30 index=1,ichan
        i=ippv(index)
        j=ippw(index)
        indm=ipm(index)
        indn=ipn(index)
        k=ipk(index)

        ka=kv(i)
        na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        kb=kv(j)
        nb=nv(j)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        call st (km,nnm)
        call st (kn,nnn)
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0
        do 40 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>
          ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

          call odd (lb+lc+k,i1)
          call odd (lm+la+k,i2)
          if (kc.eq.kn.and.(i1*i2).ne.0) then
           iu=ipcn(j,i,ii,indm,k)
	if (iu.eq.0) then
	 write (*,*) 'Term 52 - missing Xcnvw channel, STOP'
	 stop
	endif
           do 4 n=nnn,nmaxx(indn)
             do 41 m=nnm,nmaxx(indm)
* >>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            eu=1.d0/(ec-ee(n,indn)+eevt-eev+eewt-eew)
            tt=-cq*xcn(m,iu)*rmi(n,ii)*eu

            xvw(m,n,index)=xvw(m,n,index)+tt
            dxvw(m,n,index)=dxvw(m,n,index)-tt*eu
* >>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

41           continue
4          continue
          endif
40      continue
30    continue
901   continue
      return
      end

      subroutine rhoval(kval)
      implicit double precision (a-h,o-z)
      include "allvw.par"
      common /nmaxx/nmaxx(nx)
****************************** COMPLETE v ARRAY ************************
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      DIMENSION dxvi9(nxx,nxx,nvh),drvi(nx)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DIMENSION ivb(NVH),ivm(NVH),ivn(NVH),ivl(NVH),rt(nx)
      DIMENSION ipxv9(nn,nk,nk,0:kk),xvi9(nxx,nxx,nvh),rvi(nx)
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      index=0

      do 531 i=0,nvh2
        do 532 m=1,nxx
          do 533 n=1,nxx
            xvi(m,n,i)=0.d0
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            dxvi(m,n,i)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

533       continue
532     continue
531   continue

      do 661 i5=0,10
       do 662 i4=1,nk
         do 663 i3=1,nk
         do 664 i2=1,nn
          do 665 i1=1,nz
                ipxv(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue
      write (*,*) 'READING VALENCE VAL2 FILE'
	 write (*,*)
      OPEN (unit=9,file='val2',form='unformatted')
      do 1 ii=1,nval
       nvv=nv(ii)
	 kvv=kv(ii)
        CALL klj(kvv,kapv,lv,jv,indv,n0v)
       read (9) nnval,kval,ehf,dehf
	write (*,'(i3,1x,i2,2f13.7)') nnval,kval,ehf,dehf
       if (nnval.ne.nvv.or.kval.ne.kvv) then 
	  write (*,*) 'Mismatch in valence SD file, STOP'
	  stop
	 endif
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       read (9) kval,evt9
        evt(ii) = evt9
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       read (9) nmax9,lmax9,ivchan9
       read (9) ivv,ivb,ivm,ivn,ivl
       read (9) ipxv9

       do 27 i=1,ivchan9
         j=ivb(i)
         indm=ivm(i)
         indn=ivn(i)
         k=ivl(i)
         index=index+1
	if (index.gt.NVH2) then 
       write (*,*) 'Parameter NVH2 is exceeded, STOP'
       write (*,*) 'See valsd-ci output for total number of valence'
	 write (*,*) 'channels, reset  NVH2'
	 STOP
	endif
	   ipxv(ii,j,indm,indn,k)=index

         do 37 n=1,nmaxx(indn)
           read (9) (xvi9(m,n,i),m=1,nmaxx(indm))
           do 77 m=1,nmaxx(indm)
             xvi(m,n,index)=xvi9(m,n,i)
77        continue
37      continue
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         do 317 n=1,nmaxx(indn)
           read (9) (dxvi9(m,n,i),m=1,nmaxx(indm))
           do 717 m=1,nmaxx(indm)
             dxvi(m,n,index)=dxvi9(m,n,i)
717        continue
317      continue
27      continue
*********** MAY NEED TO BE CORRECTED ***********************************************
      read (9) (rvi(m),m=1,nmax9)
      read (9) (drvi(m),m=1,nmax9)
c      read (9) (rt(m),m=1,nmax9)
c      read (9) (rvi(m),m=1,nmaxx(indv))
**************NEXT LINE - CORRECTLY TRUNCATED ARRAY *********************************
      do 79 m=1,nmaxx(indv)
        rsigma(ii,m)=rvi(m)
        drsigma(ii,m)=drvi(m)
79    continue
*************************************************
1     continue
      close (9)
      ivchan=index
	 write (*,*)
	 write (*,122) index,nvh2
122   format ('Number of rho_k(mnva) channels    = ',i9,
     *' < NVH2 = ',i9,' OK')
      write (*,*)

      return
      end

	subroutine xfor4
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      DIMENSION ga(nhf),fa(nhf),gb(nhf),fb(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval

c      write (*,*) 'lin =',lin
      index=0
      index1=0
      do 531 i=0,NXV2
        do 532 m=1,nxx
          do 533 n=1,nxx
            x1(m,n,i)=0.d0
            x2(m,n,i)=0.d0
533       continue
532     continue
531   continue
      do 661 i5=0,10
       do 662 i4=1,nk
         do 663 i3=1,nk
         do 664 i2=1,nz
          do 665 i1=1,nn
                ipx1(i1,i2,i3,i4,i5)=0
                ipx2(i1,i2,i3,i4,i5)=0
665           continue
664          continue
663       continue
662     continue
661   continue

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 12 il=1,max
          ga(il)=go(il,i)
          fa(il)=fo(il,i)
12      CONTINUE
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          DO 13 il=1,max
            gb(il)=gv(il,j)
            fb(il)=fv(il,j)
13        CONTINUE

          do 191 ij=1,max
            v2(ij)=ga(ij)*gb(ij)+fa(ij)*fb(ij)
191       continue

          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnav) m,n can not be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                index=index+1
	
	if (index.gt.NXV2) then
	 write (*,*) 'EXCEED MAX CHANNEL NUMBER in XFOR4-1, STOP'
	stop
      endif
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
                index1=index1+1

	if (index1.gt.NXV2) then
	 write (*,*) 'EXCEED MAX CHANNEL NUMBER in XFOR4, STOP'
	 stop
      endif

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
61                 continue
31               continue

515             continue
111            continue
************************************************************

5           continue
2         continue
4       continue
1     continue
c      write (*,*) 'XFOR4number of channels',index,index1
901   continue
      return
      end

      subroutine count_xfor4
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval

      index=0
      index1=0

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,nval
          kb=kv(j)
          nb=nv(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 2 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            do 5 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)

****************************************************************
*             Calculate  Xk(mnav) m,n can not be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 11 k=kmin,kmax
                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                if (i1*i2.eq.0) goto 555
                index=index+1
555             continue
11            continue
****************************************************************
*             Calculate  Xk(manb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
              kmax=min0((kapm+kapn-1),(kapa+kapb-1))
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
1     continue
	if (index.le.NXV2.and.index1.le.NXV2) then 
	 write (*,122) index,nxv2
	 write (*,123) index1,nxv2
      else
	 write (*,124) index,nxv2
	 write (*,125) index1,nxv2
       write (*,*) 'Parameter NVH2 is exceeded, STOP'
	 stop
	endif
122   format ('Number of X_k(mnav) channels      = ',i9,
     *' < NXV2 = ',i9,' OK')
123   format ('Number of X_k(manv) channels      = ',i9,
     *' < NXV2 = ',i9,' OK')

124   format ('Number of X_k(mnav) channels      = ',i9,
     *' > NXV2 = ',i9)
125   format ('Number of X_k(manv) channels      = ',i9,
     *' > NXV2 = ',i9)

      write (*,*)
   
********************************************************************

      return
      end

***********************
      subroutine term41a
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
	DIMENSION  dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      DIMENSION x(nx,nx),y(nx,nx)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
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
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
                  dy(m,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue
              index2=ipxv(i,ii,indm,indr,k)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 m=nnm,nmaxx(indm)
************************* DIRECT 1 ONLY ************************************
                   y(m,nr)=y(m,nr)+xvi(m,nr,index2)
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
                  dy(m,nr)=dy(m,nr)+dxvi(m,nr,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*************************************************************
10              continue
12            continue

              k1min=max0(iabs(kapm-kapc),iabs(kapr-kapa))
              k1max=min0((kapm+kapc-1),(kapr+kapa-1))
              DO 9 k1=k1min,k1max
                z2=(2.d0*k+1.d0)*d6j(jm,ja,2*k,jr,jc,2*k1)

         index2=ipxv(i,ii,indr,indm,k1)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 m=nnm,nmaxx(indm)
************************** EXCHANGE 1 OUT ************************************
                     y(m,nr)=y(m,nr)+z2*xvi(nr,m,index2)
* >>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>
                  dy(m,nr)=dy(m,nr)+z2*dxvi(nr,m,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

***************************************************************
101                continue
121              continue
9              continue
               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                 do 23 indn=1,2*lmax+1
                   CALL indk1(indn,kn,kapn,ln,jn,n0)
                   call st (kn,nnn)
                   call trgi(iabs(kapn-kapb),(kapn+kapb-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (ln+lb,i5)
                   call odd (lm+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0


                   do 165 n=nnn,nmaxx(indn)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,n)=0.d0
164                  continue
165                continue


                   call odd (lc+lr+k,i5)
                   call odd (ln+lb+k,i6)
                   if (i5*i6.eq.0) goto 800

                    index1=ipx1(ii,j,indr,indn,k)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif

                   z2=zf*((-1)**(kapc-kapr))
                   do 123 n=nnn,nmaxx(indn)
                     do 103 nr=nnr,nmaxx(indr)
************************* DIRECT 2 ONLY ************************************
                        x(nr,n)=x(nr,n)+z2*x1(nr,n,index1)
****************************************************************************
103                  continue
123                continue
800                continue

                   index=ipvw(i,j,indm,indn,k)
* >>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                   do 107 m=nnm,nmaxx(indm)
                    do 127 nr=nnr,nmaxx(indr)

	eu=1.d0/(eevt+ec-ee(m,indm)-ee(nr,indr)+eewt-eew)

                      do 117 n=nnn,nmaxx(indn)

	                  tt=cq*x(nr,n)*y(m,nr)*eu
	                  tt1=cq*x(nr,n)*dy(m,nr)*eu

                        xvw(m,n,index)=xvw(m,n,index)+tt
                        dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1

117                     continue
127                   continue
107                 continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
700              continue
23               continue
33             continue
900           continue
500           continue
21          continue
22        continue
31      continue
336    continue
      return
      end
**************

***********************
      subroutine term41c
      implicit double precision (a-h,o-z)
      include "allvw.par"

      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
      DIMENSION dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      DIMENSION x(nx,nx),y(nx,nx)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
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
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(m,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue
              index2=ipxv(i,ii,indm,indr,k)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 m=nnm,nmaxx(indm)
************************* DIRECT 1 IN ************************************
                   y(m,nr)=y(m,nr)+xvi(m,nr,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(m,nr)=dy(m,nr)+dxvi(m,nr,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*************************************************************
10              continue
12            continue

               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                 do 23 indn=1,2*lmax+1
                   CALL indk1(indn,kn,kapn,ln,jn,n0)
                   call st (kn,nnn)
                   call trgi(iabs(kapn-kapb),(kapn+kapb-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (ln+lb,i5)
                   call odd (lm+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0


                   do 165 n=nnn,nmaxx(indn)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,n)=0.d0
164                  continue
165                continue

                   z2=zf*((-1)**(kapc-kapr))

                   zz=(-1)**(kapb+kapc+kapr+kapn)
                   k2min=max0(iabs(kapc-kapb),iabs(kapn-kapr))
                   k2max=min0((kapc+kapb-1),(kapn+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+lb+k2,i1)
                     call odd (ln+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155
                      index1=ipx2(ii,j,indn,indr,k2)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif


                     z1=zf*(2*k+1.d0)*d6j(jb,jn,2*k,jr,jc,2*k2)
                     do 124 n=nnn,nmaxx(indn)
                       do 104 nr=nnr,nmaxx(indr)
************************* EXCHANGE 2 OUT ************************************
                          x(nr,n)=x(nr,n)+z1*x2(n,nr,index1)
*****************************************************************************
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipvw(i,j,indm,indn,k)

* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    do 127 nr=nnr,nmaxx(indr)
                        do 107 m=nnm,nmaxx(indm)
	eu=1.d0/(eevt+ec-ee(m,indm)-ee(nr,indr)+eewt-eew)

                      do 117 n=nnn,nmaxx(indn)
	                 tt=cq*x(nr,n)*y(m,nr)*eu
	                 tt1=cq*x(nr,n)*dy(m,nr)*eu

                       xvw(m,n,index)=xvw(m,n,index)+tt
                       dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1

117                     continue
107                   continue
127                 continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
700              continue
23               continue
33             continue
900           continue
500           continue
21          continue
22        continue
31      continue
336    continue
      return
      end
**************
***********************
      subroutine term41d
      implicit double precision (a-h,o-z)
      include "allvw.par"

      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
	DIMENSION dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      DIMENSION x(nx,nx),y(nx,nx)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 22 indr=1,2*lmax+1
            CALL indk1(indr,kr,kapr,lr,jr,n0r)
            call st (kr,nnr)
            do 21 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn)
              kmin=max0(iabs(kapr-kapc),iabs(kapn-kapa))
              kmax=min0((kapr+kapc-1),(kapn+kapa-1))

              DO 500 l=kmin,kmax
                call odd (lc+lr,i5)
                call odd (ln+la,i6)
                if ((i5.eq.0.and.i6.ne.0).or.
     *              (i5.ne.0.and.i6.eq.0)) goto 900

***********************************************

              zf=((-1)**(l+kapr+kapc))/(2.d0*l+1.d0)
              do 171 nr=nnr,nmaxx(indr)
                do 172 n=nnn,nmaxx(indn)
                  y(n,nr)=0.d0
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

172             continue
171           continue

              k1min=max0(iabs(kapn-kapc),iabs(kapr-kapa))
              k1max=min0((kapn+kapc-1),(kapr+kapa-1))
              DO 9 k1=k1min,k1max
                z2=(2.d0*l+1.d0)*d6j(jn,ja,2*l,jr,jc,2*k1)

         index2=ipxv(i,ii,indr,indn,k1)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 n=nnn,nmaxx(indn)
************************** EXCHANGE 1 OUT ************************************
                     y(n,nr)=y(n,nr)+z2*xvi(nr,n,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  dy(n,nr)=dy(n,nr)+z2*dxvi(nr,n,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

***************************************************************
101                continue
121              continue
9              continue
               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                 do 23 indm=1,2*lmax+1
                   CALL indk1(indm,km,kapm,lm,jm,m0)
                   call st (km,nnm)
                   call trgi(iabs(kapm-kapb),(kapm+kapb-1),l,i1)
                   if (i1.eq.0) goto 700
                   call odd (ln+lb,i5)
                   call odd (lm+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700

****************8 Insert k cycle here ***********************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 505 k=kmin,kmax
*************************************************************
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0
          brr=(2.d0*k+1.d0)*d6j(jm,ja,2*k,jn,jb,2*l)


                   do 165 m=nnm,nmaxx(indm)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,m)=0.d0
164                  continue
165                continue

                   k2min=max0(iabs(kapc-kapb),iabs(kapm-kapr))
                   k2max=min0((kapc+kapb-1),(kapm+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+lb+k2,i1)
                     call odd (lm+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155
                      index1=ipx2(ii,j,indm,indr,k2)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif


                     z1=zf*(2.d0*l+1.d0)*d6j(jb,jm,2*l,jr,jc,2*k2)
                     do 124 m=nnm,nmaxx(indm)
                       do 104 nr=nnr,nmaxx(indr)
************************* EXCHANGE 2 in ************************************
                          x(nr,m)=x(nr,m)+z1*x2(m,nr,index1)
*****************************************************************************
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipvw(i,j,indm,indn,k)

* >>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
	eu=1.d0/(eevt+ec-ee(n,indn)-ee(nr,indr)+eewt-eew)

                        do 107 m=nnm,nmaxx(indm)
	                    tt=brr*cq*x(nr,m)*y(n,nr)*eu
	                    tt1=brr*cq*x(nr,m)*dy(n,nr)*eu

                         xvw(m,n,index)=xvw(m,n,index)+tt
                         dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1

107                     continue
117                   continue
127                 continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
505              continue
700              continue
23               continue
33             continue
900           continue
500           continue
21          continue
22        continue
31      continue
336    continue
      return
      end
**************
***********************************************************************************
      subroutine term42a
      implicit double precision (a-h,o-z)
      include "allvw.par"
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
	DIMENSION dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      DIMENSION x(nx,nx),y(nx,nx)
               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)

        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
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
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue
              index2=ipxv(j,ii,indn,indr,k)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 n=nnn,nmaxx(indn)
************************* DIRECT 1 IN ************************************
                   y(n,nr)=y(n,nr)+xvi(n,nr,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=dy(n,nr)+dxvi(n,nr,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
10              continue
12            continue

              k1min=max0(iabs(kapn-kapc),iabs(kapr-kapb))
              k1max=min0((kapn+kapc-1),(kapr+kapb-1))
              DO 9 k1=k1min,k1max
                z2=(2.d0*k+1.d0)*d6j(jn,jb,2*k,jr,jc,2*k1)

         index2=ipxv(j,ii,indr,indn,k1)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 n=nnn,nmaxx(indn)
************************* DIRECT 1 OUT ************************************
                     y(n,nr)=y(n,nr)+z2*xvi(nr,n,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=dy(n,nr)+z2*dxvi(nr,n,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

***************************************************************************
101                continue
121              continue
9              continue

                 do 23 indm=1,2*lmax+1
                   CALL indk1(indm,km,kapm,lm,jm,m0)
                   call st (km,nnm)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*******************************
                   call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (lm+la,i5)
                   call odd (ln+lb,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0


                   do 165 m=nnm,nmaxx(indm)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,m)=0.d0
164                  continue
165                continue


                   call odd (lc+lr+k,i5)
                   call odd (lm+la+k,i6)
                   if (i5*i6.eq.0) goto 800

                    index1=ipx1(ii,i,indr,indm,k)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif

                   z2=zf*((-1)**(kapc-kapr))
                   do 123 m=nnm,nmaxx(indm)
                     do 103 nr=nnr,nmaxx(indr)
************************* DIRECT 2 IN ************************************
                        x(nr,m)=x(nr,m)+z2*x1(nr,m,index1)
***************************************************************************
103                  continue
123                continue
800                continue

                   index=ipvw(i,j,indm,indn,k)
                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
* >>>>>>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	eu=1.d0/(eewt+ec-ee(n,indn)-ee(nr,indr)+eevt-eev)
                        do 107 m=nnm,nmaxx(indm)
	             tt=cq*x(nr,m)*y(n,nr)*eu
	             tt1=cq*x(nr,m)*dy(n,nr)*eu

                   xvw(m,n,index)=xvw(m,n,index)+tt
                   dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1
                                   
107                     continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
117                   continue
127                 continue
700              continue
336             continue
23               continue
900           continue
500           continue
21          continue
22        continue
31      continue
33    continue
      return
      end
**********************************

***************************
      subroutine term42c
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
	DIMENSION dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      DIMENSION x(nx,nx),y(nx,nx)
               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)

        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
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
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue
              index2=ipxv(j,ii,indn,indr,k)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 n=nnn,nmaxx(indn)
************************* DIRECT 1 IN ************************************
                   y(n,nr)=y(n,nr)+xvi(n,nr,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(n,nr)=dy(n,nr)+dxvi(n,nr,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
10              continue
12            continue

                 do 23 indm=1,2*lmax+1
                   CALL indk1(indm,km,kapm,lm,jm,m0)
                   call st (km,nnm)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*******************************
                   call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1)
                   if (i1.eq.0) goto 700
                   call odd (lm+la,i5)
                   call odd (ln+lb,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700
        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0


                   do 165 m=nnm,nmaxx(indm)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,m)=0.d0
164                  continue
165                continue


                   z2=zf*((-1)**(kapc-kapr))

                   zz=(-1)**(kapa+kapc+kapr+kapm)
                   k2min=max0(iabs(kapc-kapa),iabs(kapm-kapr))
                   k2max=min0((kapc+kapa-1),(kapm+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+la+k2,i1)
                     call odd (lm+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155
                      index1=ipx2(ii,i,indm,indr,k2)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif


                     z1=zf*(2.d0*k+1.d0)*d6j(ja,jm,2*k,jr,jc,2*k2)
                     do 124 m=nnm,nmaxx(indm)

                       do 104 nr=nnr,nmaxx(indr)
*************************  EXCHANGE 2 OUT ************************************
                          x(nr,m)=x(nr,m)+z1*x2(m,nr,index1)
******************************************************************************
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipvw(i,j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
* >>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          eu=1.d0/(eewt+ec-ee(n,indn)-ee(nr,indr)+eevt-eev)

			            do 107 m=nnm,nmaxx(indm)
	           tt=cq*x(nr,m)*y(n,nr)*eu
	           tt1=cq*x(nr,m)*dy(n,nr)*eu

                 xvw(m,n,index)=xvw(m,n,index)+tt
                 dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1

107                     continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
117                   continue
127                 continue
700              continue
336             continue
23               continue
900           continue
500           continue
21          continue
22        continue
31      continue
33    continue
      return
      end
***************************
      subroutine term42d
      implicit double precision (a-h,o-z)
      include "allvw.par"
      COMMON /xcore1/ ipx1(nn,nz,nk,nk,0:kk),x1(nxx,nxx,0:NXV2),ichan1
      COMMON /xcore2/ ipx2(nn,nz,nk,nk,0:kk),x2(nxx,nxx,0:NXV2),ichan2
      COMMON /irhov/ ipxv(nz,nn,nk,nk,0:kk),ivchan
      COMMON /rhov/ xvi(nxx,nxx,0:nvh2),rsigma(nz,nx)

      COMMON /radvw/ ipvw(nz,nz,nk,nk,0:kk),xvw(nxx,nxx,0:NVWH),ichan
      COMMON /val/ gv(nhf,nz),fv(nhf,nz),ev(nz),nv(nz),kv(nz),nval
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore
* >>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON /drhov/ dxvi(nxx,nxx,0:nvh2),drsigma(nz,nx)
	COMMON /denergy/ evt(nz)
      COMMON /dradvw/ dxvw(nxx,nxx,0:NVWH)
	DIMENSION dy(nx,nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      DIMENSION x(nx,nx),y(nx,nx)
               do 33 j=1,nval
                 kb=kv(j)
                 nb=nv(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)

        do 31 ii=1,ncore
          kc=ko(ii)
          nc=no(ii)
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          do 22 indr=1,2*lmax+1
            CALL indk1(indr,kr,kapr,lr,jr,n0r)
            call st (kr,nnr)
            do 21 indm=1,2*lmax+1
              CALL indk1(indm,km,kapm,lm,jm,m0)
              call st (km,nnm)
              kmin=max0(iabs(kapr-kapc),iabs(kapm-kapb))
              kmax=min0((kapr+kapc-1),(kapm+kapb-1))
              DO 500 l=kmin,kmax
                call odd (lc+lr,i5)
                call odd (lm+lb,i6)
                if ((i5.eq.0.and.i6.ne.0).or.
     *              (i5.ne.0.and.i6.eq.0)) goto 900

***********************************************

              zf=((-1)**(l+kapr+kapc))/(2.d0*l+1.d0)
              do 171 nr=nnr,nmaxx(indr)
                do 172 m=nnm,nmaxx(indm)
                  y(m,nr)=0.d0
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(m,nr)=0.d0
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
172             continue
171           continue

              k1min=max0(iabs(kapm-kapc),iabs(kapr-kapb))
              k1max=min0((kapm+kapc-1),(kapr+kapb-1))
              DO 9 k1=k1min,k1max
                z2=(2.d0*l+1.d0)*d6j(jm,jb,2*l,jr,jc,2*k1)

         index2=ipxv(j,ii,indr,indm,k1)
          if (index2.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH rho_mnva, STOP'
	     stop
	    endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 m=nnm,nmaxx(indm)
************************* DIRECT 1 in************************************
                     y(m,nr)=y(m,nr)+z2*xvi(nr,m,index2)
* >>>>>>>>>>>>>>>>SIGMA >>>>>>>
                  dy(m,nr)=dy(m,nr)+z2*dxvi(nr,m,index2)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
***************************************************************************
101                continue
121              continue
9              continue

                 do 23 indn=1,2*lmax+1
                   CALL indk1(indn,kn,kapn,ln,jn,n0)
                   call st (kn,nnn)
      do 336 i=1,nval
          ka=kv(i)
          na=nv(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
* >>>>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>
	  eev=ev(i)
	  eew=ev(j)
	  eevt=evt(i)
	  eewt=evt(j)
	  ec=eo(ii)
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


*******************************
                   call trgi(iabs(kapn-kapa),(kapn+kapa-1),l,i1)
                   if (i1.eq.0) goto 700
                   call odd (lm+la,i5)
                   call odd (ln+lb,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 700
****************8 Insert k cycle here ***********************
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 505 k=kmin,kmax
*************************************************************
          brr=(2.d0*k+1.d0)*d6j(jm,ja,2*k,jn,jb,2*l)

        c0=((-1)**(k))*s1(k,km,ka)*s1(k,kn,kb)
        cq=1.d0/c0


                   do 165 n=nnn,nmaxx(indn)
                     do 164 nr=nnr,nmaxx(indr)
                       x(nr,n)=0.d0
164                  continue
165                continue

                   k2min=max0(iabs(kapc-kapa),iabs(kapn-kapr))
                   k2max=min0((kapc+kapa-1),(kapn+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+la+k2,i1)
                     call odd (ln+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155
                      index1=ipx2(ii,i,indn,indr,k2)
          if (index1.eq.0) then
	     write (*,*) 'TERM41 - NO SUCH g_mncv, STOP'
	     stop
	    endif


                     z1=zf*(2.d0*l+1.d0)*d6j(ja,jn,2*l,jr,jc,2*k2)
                     do 124 n=nnn,nmaxx(indn)

                       do 104 nr=nnr,nmaxx(indr)
*************************  EXCHANGE 2 in ************************************
                          x(nr,n)=x(nr,n)+z1*x2(n,nr,index1)
******************************************************************************
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipvw(i,j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                        do 107 m=nnm,nmaxx(indm)
* >>>>>>>>>>>> SIGMA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	eu=1.d0/(eewt+ec-ee(m,indm)-ee(nr,indr)+eevt-eev)

                      do 117 n=nnn,nmaxx(indn)

	            tt= brr*cq*x(nr,n)*y(m,nr)*eu
	            tt1= brr*cq*x(nr,n)*dy(m,nr)*eu

                  xvw(m,n,index)=xvw(m,n,index)+tt
                  dxvw(m,n,index)=dxvw(m,n,index)-tt*eu+tt1
                                
117                     continue
107                   continue
* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

127                 continue
505           continue
700          continue
336         continue
23               continue
900           continue
500           continue
21          continue
22        continue
31      continue
33    continue
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
c     

c     ========================================
      include "d6j.f"
      include "libD.f"
      include "rint.f"
      include "yfun.f"
      include "yint.f"
	include "inidat.f"
