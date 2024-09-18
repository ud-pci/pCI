* this is allcore-rlen.f code renamed just for convinience
****** revised version of April 20, 2012, extra parameters kk1, nhf1, nx are removed
** lowdin version - option 2 (lapack) commented out   
** NEW VERSION - January 15 modified to run with CI+all-order basis,
* i.e. different number of orbitals for different partial waves
* This is fully modified version *
************************************
      implicit real*8 (a-h,o-z)    
      include "all.par"
      DIMENSION t(nkx)
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)    
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore 
      DIMENSION g(nhf,nx),f(nhf,nx),e(nl)
      DIMENSION wco(ns),nco(ns),kco(ns),mco(ns)
      COMMON /x12ini/ inmax
      common /nmaxx/nmaxx(nx)

      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn) 
      COMMON /rho_all/ xco_all(nxx,nxx,nch,nit),rmo_all(nxx,nn,nit)  

      OPEN (unit=1,file='hfspl.1',form='unformatted')
      OPEN (unit=2, file='hfspl.2', form='unformatted',   
     *access='direct',recl=nrec) 
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
****************************************************************
      read (5,*) nmax,lmax
      read (5,*) key1,key2
      read (5,*) nitmax
	read (5,*) ikey,max_rle,keym
      read (5,*) damp
      

c********** NEW 2012 ********
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
c      read(*,*)
      call countx
      call countrho

*********END NEW 2012************



      call rhocore
      call rhoini(key1)
17    FORMAT ('Time required = ',f9.3,'sec')
      call energy(deltae)
      write (*,'(f15.8)') deltae

cRLE******************************************************************
       irle=0
        
cRLE*************************************************************
      do 1000 ik=1,nitmax
        
        call outin
        call terma1
        call terma2
        call terma3
        call terma4
  

        call term1

        t1=mclock()
        call term21
        t2=mclock()
        td=(t2-t1)/1000
c        WRITE (*,17) td 


        call term31
        call term32
        call term41
        call term42
        call term5(1)
        call term5(2)


        call outt(deltae,damp,ik,key2)

        call energy(deltae1)
        or=abs((deltae-deltae1)/deltae)
        write (*,907) ik,deltae,deltae1,or
907     format ('Iter',i5,3f15.8)
        call outrho
c        if (or.lt.0.0000001) goto 3000
        if (or.lt.0.00001) goto 3000
        deltae=deltae1

c RLE **********Accumulate m iterations ****************
       if (ikey.ne.0) then 

        irle=irle+1

        do 27 i=1,icchan
          do 37 m=1,nxx      
            do 78 n=1,nxx
              xco_all(m,n,i,irle)=xco(m,n,i)
78          continue
37        continue
27      continue

      do 47 i=1,ncore
        do 718 n=1,nxx
        rmo_all(n,i,irle)=rmo(n,i)
718     continue
47     continue

***********************************************************
777    continue
c      call rho_read
       if (irle.eq.max_rle) then 
        call energy_rle(max_rle)

        if (ikey.eq.1) then 
         write (*,*)
         write (*,*) ' RLE CONVERGENCE SCHEME IS USED'
         if (keym.eq.1) then 
         write (*,*) ' Using linpack codes'
         endif
         if (keym.eq.2) then 
         write (*,*) ' Using lapack codes'
         endif
         write (*,*) ' Number of iterations stored = ',max_rle
         write (*,*) ' Number of coefficents tau = ',max_rle-1
         write (*,*)
        endif

        if (ikey.eq.2) then 
         write (*,*)
         write (*,*) ' DIIS CONVERGENCE SCHEME IS USED'
         if (keym.eq.1) then 
         write (*,*) ' Using linpack codes'
         endif
         if (keym.eq.2) then 
         write (*,*) ' Using lapack codes'
         endif
         write (*,*) ' Number of iterations stored = ',max_rle
         write (*,*) ' Number of coefficents tau = ',max_rle-1
         write (*,*)
        endif

        if (ikey.eq.1) then 
         call rle(ene,max_rle,keym)
        endif

        if (ikey.eq.2) then 
         call diis1(ene,max_rle,keym)
        endif


        deltae=ene
        deltae1=ene
        irle=0
       endif
      endif
*********************************************
1000  continue
3000  continue
      write (*,'(f20.10)')  deltae1
      call outrho

      stop
      end

      subroutine rle(ene,max_rle,keym)
      implicit double precision (a-h,o-z)
      integer*4 ipvt
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      COMMON /rho_all/ xco_all(nxx,nxx,nch,nit),rmo_all(nxx,nn,nit)  

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      dimension Rm(nit,nit),al(nit),bl(nit,nyy),cl(nit)
      dimension z(nit),IPIV(nit),ipvt(nit)

      do 112 i=1,nit
        al(i)=0.d0
        cl(i)=0.d0
        bl(i,1)=0.d0
        do 113 j=1,nit
          Rm(i,j)=0.d0
113       continue
112    continue


c ********* Accumulate  al

      nt=max_rle-1
      do 1600 i=1,nt

      do 11 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 12 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
             if (i1.ge.i2) then 
          do 13 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 14 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
*            Initialize  rhok(mnab)                            * 
****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 15 lk=kmin,kmax

                call odd (lm+la+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 16
               
          if (i1.ge.i2) then
            index1=ipx1(i1,i2,indm,indn,lk)
            index2=ipxc(i1,i2,indm,indn,lk)                 
          else
            index1=ipx1(i2,i1,indn,indm,lk)
            index2=ipxc(i2,i1,indn,indm,lk)                 
           endif            
c            zz=1.d0    
           zz=1.d0/(2.d0*lk+1.d0)
           do 17 m=nnm,nmaxx(indm)
             do 18 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
               if (i1.ge.i2) then 
                
                 al(i)=al(i)-zz*xco_all(m,n,index2,i)*x1(m,n,index1)  

                else

               al(i)=al(i)-zz*xco_all(n,m,index2,i)*x1(n,m,index1)  
              
                endif                     
18            continue 
17          continue
16       continue
15            continue 
14          continue
13       continue
          endif
12          continue
11       continue
1600  continue

c ********* Accumulate Rm 

  
      do 1000 i=1,nt
       do 2000 j=1,nt

       do 1 i1=1,ncore 
         ka=ko(i1)
         na=no(i1)
         CALL klj(ka,kapa,la,ja,inda,n0a)
         ea=eo(i1)

         indn=inda
         CALL indk1(indn,kn,kapn,ln,jn,n0)
         call st (kn,nnn) 
         do 61 n=nnn,nmaxx(indn)	      
           den=ea-ee(n,indn)
           Rm(i,j)=Rm(i,j)+rmo_all(n,i1,i)*den*rmo_all(n,i1,j)-
     *                     rmo_all(n,i1,i)*den*rmo_all(n,i1,j+1)
61     continue

1      continue

      do 21 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 22 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)      
          if (i1.ge.i2) then 

          do 23 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 24 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 25 lk=kmin,kmax

                call odd (lm+la,ii1)
                call odd (ln+lb,ii2)
                if ((ii1.eq.0.and.ii2.ne.0).or.
     *              (ii1.ne.0.and.ii2.eq.0)) goto 26

          if (i1.ge.i2) then
            index2=ipxc(i1,i2,indm,indn,lk)                 
          else
            index2=ipxc(i2,i1,indn,indm,lk)                 
           endif            
    
           zz=1.d0/(2.d0*lk+1.d0)
c           zz=1.d0
           do 27 m=nnm,nmaxx(indm)
             do 28 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
               if (i1.ge.i2) then 
                
      Rm(i,j)=Rm(i,j)+zz*xco_all(m,n,index2,i)*den*xco_all(m,n,index2,j)
     *             -zz*xco_all(m,n,index2,i)*den*xco_all(m,n,index2,j+1)

                else

      Rm(i,j)=Rm(i,j)+zz*xco_all(n,m,index2,i)*den*xco_all(n,m,index2,j)
     *             -zz*xco_all(n,m,index2,i)*den*xco_all(n,m,index2,j+1)
              
                endif  
28            continue 
27          continue			                     
26          continue 
25          continue
24          continue
23          continue 
            endif
22          continue
21       continue

2000    continue
1000  continue

      do 212 i=1,nt
        do 213 j=1,nt
          Rm(i,j)=Rm(i,j)-al(i)
213       continue
212    continue
***************************************


      write (*,*) '*********************************************'
      write (*,*) ' Print RLE output'
      write (*,*) '*********************************************'
      write (*,*)
      write (*,*) 'Print alpha, alpha + R tau = 0'
      write (*,*)

      do 801 i3=1,nt
       al(i3)=al(i3)-ak
       write (*,'(i4,f12.6)') i3,al(i3)
        bl(i3,1)=-al(i3)
        cl(i3)=-al(i3)

801    continue
      write (*,*)
      write (*,*) 'Print R'
      write (*,*)

      do 802 i3=1,nt
       write (*,'(i4,8f12.6)') i3,(Rm(i3,i4),i4=1,nt)
802    continue

****************************** LINPACK VERSION ***************
      if (keym.eq.1) then 
      call dgeco(Rm,nit,nt,ipvt,rcond,z)
      ising=0
      write (6,*) 'rcond=',rcond
      if (rcond.lt.1.d-18) ising=1
      if (ising.eq.1) then 
	 write (*,*) 'Singularity'
	 goto 901
      endif
      call dgesl(Rm,nit,nt,ipvt,cl,0)
	cc=0.d0
      write (*,*)
      write (*,*) 'Print coefficients tau'
      write (*,*)
      do 80 i3=1,nt
	 cc=cc+cl(i3)
         al(i3)=cl(i3)
       write (*,'(i4,f12.6)') i3,al(i3)
80    continue
      write (*,*)
90    format ('Total = ',f12.6)
      write (*,90) cc
      write (*,*) 
      endif
****************************************************************

****************************** LAPACK VERSION ***************
      if (keym.eq.2) then 

      NRHS=1
c      call DGESV(nit,NRHS,Rm,nit,ipiv,bl,nit,INFO)
      write (*,*) 'info',info
      if (INFO.ne.0) then 
       write (*,*) 'INFO is not ZERO, stop'
       stop
      endif

	cc=0.d0
      write (*,*)
      write (*,*) 'Print coefficients tau'
      write (*,*)

      do 808 i3=1,nt
	 cc=cc+bl(i3,1)
         al(i3)=bl(i3,1)
       write (*,'(i4,f12.6)') i3,al(i3)
808    continue
      write (*,*)
      write (*,90) cc
      endif
********** Define new rho *************
      
      do 275 i=1,icchan
        do 37 n=1,nxx
          do 78 m=1,nxx
              xco(m,n,i)=0.d0
            do 85 i3=1,nt
              xco(m,n,i)=xco(m,n,i)+al(i3)*xco_all(m,n,i,i3)
85          continue
78         continue
37       continue
275    continue


      do 47 i=1,ncore
       do 718 n=1,nxx
        rmo(n,i)=0.d0
         do 905 i3=1,nt
           rmo(n,i)=rmo(n,i)+al(i3)*rmo_all(n,i,i3)
905       continue
718      continue
47     continue

***************************************

      call energy(ene)
      write (*,850) ene
      write (*,*)
850   format ('New energy = ',e15.8)
      write (*,*) 'END OF RLE output'
***************************************
901   continue
      return
      end


      subroutine diis(ene,max_rle,keym)
      implicit double precision (a-h,o-z)
      integer*4 ipvt
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      COMMON /rho_all/ xco_all(nxx,nxx,nch,nit),rmo_all(nxx,nn,nit)  

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      dimension Rm(nit,nit),al(nit),bl(nit),cl(nit)
      dimension z(nit),IPIV(nit),ipvt(nit)

      do 112 i=1,nit
        al(i)=0.d0
        bl(i)=0.d0
        do 113 j=1,nit
          Rm(i,j)=0.d0
113       continue
112    continue
********************** accumulate ak
      nt=max_rle-1
      ak=0.d0

      do 511 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 512 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
             if (i1.ge.i2) then 
          do 513 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 514 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
*            Initialize  rhok(mnab)                            * 
****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 515 lk=kmin,kmax

                call odd (lm+la+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 516
               
            index1=ipx1(i1,i2,indm,indn,lk)

           zz=1.d0/(2.d0*lk+1.d0)
           do 517 m=nnm,nmaxx(indm)
             do 518 n=nnn,nmaxx(indn)
                 ak=ak+zz*x1(m,n,index1)*x1(m,n,index1)  

518            continue 
517          continue
516       continue
515            continue 
514          continue
513       continue
          endif
512          continue
511       continue
      write (*,*) 'ak=', ak
c ********* Accumulate  al

      nt=max_rle-1
      do 1600 i=1,nt

      do 11 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 12 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
             if (i1.ge.i2) then 
          do 13 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 14 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
*            Initialize  rhok(mnab)                            * 
****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 15 lk=kmin,kmax

                call odd (lm+la+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 16
               
            index1=ipx1(i1,i2,indm,indn,lk)
            index2=ipxc(i1,i2,indm,indn,lk)                 

           zz=1.d0/(2.d0*lk+1.d0)
           do 17 m=nnm,nmaxx(indm)
             do 18 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
                
         al(i)=al(i)-zz*den*x1(m,n,index1)*
     *(xco_all(m,n,index2,i)-xco_all(m,n,index2,i+1))

                                  
18            continue 
17          continue
16       continue
15            continue 
14          continue
13       continue
          endif
12          continue
11       continue
1600  continue

c ********* Accumulate Rm 

  
      do 1000 i=1,nt
       do 2000 j=1,nt

       do 1 i1=1,ncore 
         ka=ko(i1)
         na=no(i1)
         CALL klj(ka,kapa,la,ja,inda,n0a)
         ea=eo(i1)

         indn=inda
         CALL indk1(indn,kn,kapn,ln,jn,n0)
         call st (kn,nnn) 
         do 61 n=nnn,nmaxx(indn)	      
           den=ea-ee(n,indn)
           Rm(i,j)=Rm(i,j)+den*den*
     *      (rmo_all(n,i1,i)*rmo_all(n,i1,j)
     *      +rmo_all(n,i1,i+1)*rmo_all(n,i1,j+1)
     *      -rmo_all(n,i1,i+1)*rmo_all(n,i1,j)
     *      -rmo_all(n,i1,i)*rmo_all(n,i1,j+1))
61     continue

1      continue

      do 21 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 22 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)      
          if (i1.ge.i2) then 

          do 23 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 24 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 25 lk=kmin,kmax

                call odd (lm+la,ii1)
                call odd (ln+lb,ii2)
                if ((ii1.eq.0.and.ii2.ne.0).or.
     *              (ii1.ne.0.and.ii2.eq.0)) goto 26

            index2=ipxc(i1,i2,indm,indn,lk)                 
    
           zz=1.d0/(2.d0*lk+1.d0)
           do 27 m=nnm,nmaxx(indm)
             do 28 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
                
      Rm(i,j)=Rm(i,j)+zz*den*den*(
     *       xco_all(m,n,index2,i)*xco_all(m,n,index2,j)
     *      +xco_all(m,n,index2,i+1)*xco_all(m,n,index2,j+1)
     *      -xco_all(m,n,index2,i+1)*xco_all(m,n,index2,j)
     *      -xco_all(m,n,index2,i)*xco_all(m,n,index2,j+1))
28            continue 
27          continue			                     
26          continue 
25          continue
24          continue
23          continue 
            endif
22          continue
21       continue

2000    continue
1000  continue

      do 912 i=1,nt
          al(i)=al(i)-ak
912    continue

      do 212 i=1,nt
        do 213 j=1,nt
          Rm(i,j)=Rm(i,j)-al(i)-al(j)-ak
213       continue
212    continue
***************************************


      write (*,*) '*********************************************'
      write (*,*) ' Print DIIS output'
      write (*,*) '*********************************************'
      write (*,*)
      write (*,*) 'Print alpha, alpha + R tau = 0'
      write (*,*)

      do 801 i3=1,nt
       write (*,'(i4,f12.6)') i3,al(i3)
        bl(i3)=-al(i3)

801    continue
      write (*,*)
      write (*,*) 'Print R'
      write (*,*)

      do 802 i3=1,nt
       write (*,'(i4,8f12.6)') i3,(Rm(i3,i4),i4=1,nt)
802    continue

****************************** LINPACK VERSION ***************
      if (keym.eq.1) then 
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
      write (*,*)
      write (*,*) 'Print coefficients tau'
      write (*,*)
      do 80 i3=1,nt
	 cc=cc+bl(i3)
         al(i3)=bl(i3)
       write (*,'(i4,f12.6)') i3,al(i3)
80    continue
      write (*,*)
90    format ('Total = ',f12.6)
      write (*,90) cc
      write (*,*) 
      endif
****************************************************************

********** Define new rho *************
      
      do 275 i=1,icchan
        do 37 n=1,nxx
          do 78 m=1,nxx
              xco(m,n,i)=0.d0
            do 85 i3=1,nt
              xco(m,n,i)=xco(m,n,i)+al(i3)*xco_all(m,n,i,i3)
85          continue
78         continue
37       continue
275    continue


      do 47 i=1,ncore
       do 718 n=1,nxx
        rmo(n,i)=0.d0
         do 905 i3=1,nt
           rmo(n,i)=rmo(n,i)+al(i3)*rmo_all(n,i,i3)
905       continue
718      continue
47     continue

***************************************

      call energy(ene)
      write (*,850) ene
      write (*,*)
850   format ('New energy = ',e15.8)
      write (*,*) 'END OF DIIS output'
***************************************
901   continue
      return
      end



      subroutine diis1(ene,max_rle,keym)
      implicit double precision (a-h,o-z)
      integer*4 ipvt
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      COMMON /rho_all/ xco_all(nxx,nxx,nch,nit),rmo_all(nxx,nn,nit)  

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      dimension Rm(nit,nit),al(nit),bl(nit,nyy),cl(nit)
      dimension z(nit),ipvt(nit),ipiv(nit)

      do 112 i=1,nit
        al(i)=0.d0
        cl(i)=0.d0
        bl(i,1)=0.d0
        do 113 j=1,nit
          Rm(i,j)=0.d0
113       continue
112    continue

********* Accumulate ak^2 ******
      ak=0.d0
      do 811 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 812 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)          
          if (i1.ge.i2) then 
          do 813 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 814 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
*            Initialize  rhok(mnab)                            * 
****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 815 lk=kmin,kmax

                call odd (lm+la+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 816
               
          if (i1.ge.i2) then
            index1=ipx1(i1,i2,indm,indn,lk)                
          else
            index1=ipx1(i2,i1,indn,indm,lk)                
           endif            
           zz=1.d0/(2.d0*lk+1.d0)
           do 817 m=nnm,nmaxx(indm)
             do 818 n=nnn,nmaxx(indn)
               if (i1.ge.i2) then 
                
           ak=ak+zz*x1(m,n,index1)*x1(m,n,index1)

                else

           ak=ak+zz*x1(n,m,index1)*x1(n,m,index1)
              
                endif                     
818            continue 
817          continue
816       continue
815            continue 
814          continue
813       continue
        endif

812          continue
811       continue
877     format ('ak = ',f24.7)
        write (*,877) ak
c ********* Accumulate  al

      nt=max_rle-1
      do 1600 i=1,nt

      do 11 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 12 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)       
             if (i1.ge.i2) then 

          do 13 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 14 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************
*            Initialize  rhok(mnab)                            * 
****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 15 lk=kmin,kmax

                call odd (lm+la+lk,ii1)
                call odd (ln+lb+lk,ii2)
                 
                if (ii1*ii2.eq.0) goto 16
               
          if (i1.ge.i2) then
            index1=ipx1(i1,i2,indm,indn,lk)
            index2=ipxc(i1,i2,indm,indn,lk)                 
          else
            index1=ipx1(i2,i1,indn,indm,lk)
            index2=ipxc(i2,i1,indn,indm,lk)                 
           endif            
           zz=1.d0/(2.d0*lk+1.d0)
           do 17 m=nnm,nmaxx(indm)
             do 18 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
               if (i1.ge.i2) then 
*****************
           al(i)=al(i)-
     *        zz*den*x1(m,n,index1)*(xco_all(m,n,index2,i)
     *                              -xco_all(m,n,index2,i+1))

                else

           al(i)=al(i)-
     *        zz*den*x1(n,m,index1)*(xco_all(n,m,index2,i)
     *                              -xco_all(n,m,index2,i+1))              
                endif                     
18            continue 
17          continue
16       continue
15            continue 
14          continue
13       continue
           endif
12          continue
11       continue
1600  continue

c ********* Accumulate Rm 

  
      do 1000 i=1,nt
       do 2000 j=1,nt

       do 1 i1=1,ncore 
         ka=ko(i1)
         na=no(i1)
         CALL klj(ka,kapa,la,ja,inda,n0a)
         ea=eo(i1)

         indn=inda
         CALL indk1(indn,kn,kapn,ln,jn,n0)
         call st (kn,nnn) 
         do 61 n=nnn,nmaxx(indn)	      
           den=ea-ee(n,indn)
      Rm(i,j)=Rm(i,j)+
     *        den*den*(rmo_all(n,i1,i)*rmo_all(n,i1,j)
     *                +rmo_all(n,i1,i+1)*rmo_all(n,i1,j+1)
     *                -rmo_all(n,i1,i+1)*rmo_all(n,i1,j)
     *                -rmo_all(n,i1,i)*rmo_all(n,i1,j+1))
61     continue

1      continue

      do 21 i1=1,ncore 
        ka=ko(i1)
        na=no(i1)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i1)
        do 22 i2=1,ncore  
          kb=ko(i2)
          nb=no(i2)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(i2)     
               if (i1.ge.i2) then 

          do 23 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)            
            do 24 indn=1,2*lmax+1
              CALL indk1(indn,kn,kapn,ln,jn,n0)
              call st (kn,nnn) 

****************************************************************    
              kmin=max0(iabs(kapm-kapa),iabs(kapn-kapb))
              kmax=min0((kapm+kapa-1),(kapn+kapb-1))
              DO 25 lk=kmin,kmax

                call odd (lm+la,ii1)
                call odd (ln+lb,ii2)
                if ((ii1.eq.0.and.ii2.ne.0).or.
     *              (ii1.ne.0.and.ii2.eq.0)) goto 26

          if (i1.ge.i2) then
            index2=ipxc(i1,i2,indm,indn,lk)                 
          else
            index2=ipxc(i2,i1,indn,indm,lk)                 
           endif            
    
           zz=1.d0/(2.d0*lk+1.d0)

           do 27 m=nnm,nmaxx(indm)
             do 28 n=nnn,nmaxx(indn)
               den=ea+eb-ee(m,indm)-ee(n,indn) 
               if (i1.ge.i2) then 
                
      Rm(i,j)=Rm(i,j)+zz*den*den*
     *            (xco_all(m,n,index2,i)*xco_all(m,n,index2,j)
     *            +xco_all(m,n,index2,i+1)*xco_all(m,n,index2,j+1)
     *            -xco_all(m,n,index2,i+1)*xco_all(m,n,index2,j)
     *            -xco_all(m,n,index2,i)*xco_all(m,n,index2,j+1))


                else

      Rm(i,j)=Rm(i,j)+zz*den*den*
     *            (xco_all(n,m,index2,i)*xco_all(n,m,index2,j)
     *            +xco_all(n,m,index2,i+1)*xco_all(n,m,index2,j+1)
     *            -xco_all(n,m,index2,i+1)*xco_all(n,m,index2,j)
     *            -xco_all(n,m,index2,i)*xco_all(n,m,index2,j+1))

c12345678901234567892123456789312345678901234567892123456789612345678
              
                endif  
28            continue 
27          continue			                     
26          continue 
25          continue
24          continue
23          continue 
         endif
22          continue
21       continue

2000    continue
1000  continue

      do 809 i3=1,nt
       al(i3)=al(i3)-ak
809    continue


      do 212 i=1,nt
        do 213 j=1,nt
          Rm(i,j)=Rm(i,j)-al(i)-al(j)-ak
213       continue
212    continue


      write (*,*) '*********************************************'
      write (*,*) ' Print DIIS output'
      write (*,*) '*********************************************'
      write (*,*)
      write (*,*) 'Print alpha, alpha + R tau = 0'
      write (*,*)

      do 801 i3=1,nt
       write (*,'(i4,f17.12)') i3,al(i3)
        bl(i3,1)=-al(i3)
        cl(i3)=-al(i3)

801    continue
      write (*,*)
      write (*,*) 'Print R'
      write (*,*)

      do 802 i3=1,nt
       write (*,'(i4,8f17.12)') i3,(Rm(i3,i4),i4=1,nt)
802    continue

****************************** LINPACK VERSION ***************
      write (*,*) 'keym',keym  
      if (keym.eq.1) then 
      write (*,*) 'linpack'

      call dgeco(Rm,nit,nt,ipvt,rcond,z)
      ising=0
      write (6,*) 'rcond=',rcond
      if (rcond.lt.1.d-18) ising=1
      if (ising.eq.1) then 
	 write (*,*) 'Singularity'
	 goto 901
      endif
      call dgesl(Rm,nit,nt,ipvt,cl,0)
	cc=0.d0
      write (*,*)
      write (*,*) 'Print coefficients tau'
      write (*,*)
      do 807 i3=1,nt
	 cc=cc+cl(i3)
         al(i3)=cl(i3)
       write (*,'(i4,f14.8)') i3,al(i3)
807    continue
      write (*,*)
      write (*,90) cc
      write (*,*) 
      endif
****************************************************************

****************************** LAPACK VERSION ***************
      if (keym.eq.2) then 
      write (*,*) 'lapack'
      NRHS=1
c      call DGESV(nit,NRHS,Rm,nit,ipiv,bl,nit,INFO)
      write (*,*) 'info',info
      if (INFO.ne.0) then 
       write (*,*) 'INFO is not ZERO, stop'
       stop
      endif

	cc=0.d0
      write (*,*)
      write (*,*) 'Print coefficients tau'
      write (*,*)

      do 80 i3=1,nt
	 cc=cc+bl(i3,1)
         al(i3)=bl(i3,1)
       write (*,'(i4,f14.8)') i3,al(i3)
80    continue
      write (*,*)
90    format ('Total = ',f14.8)
      write (*,90) cc
       endif
********** Define new rho *************
      
      do 275 i=1,icchan
        do 37 n=1,nxx
          do 78 m=1,nxx
              xco(m,n,i)=0.d0
            do 85 i3=1,nt
              xco(m,n,i)=xco(m,n,i)+al(i3)*xco_all(m,n,i,i3)
85          continue
78         continue
37       continue
275    continue


      do 47 i=1,ncore
       do 718 n=1,nxx
        rmo(n,i)=0.d0
         do 905 i3=1,nt
           rmo(n,i)=rmo(n,i)+al(i3)*rmo_all(n,i,i3)
905       continue
718      continue
47     continue

***************************************

      call energy(ene)
      write (*,850) ene
      write (*,*)
850   format ('New energy = ',e15.8)
      write (*,*) 'END OF DIIS output'
901   continue
      return
      end



      subroutine energy_rle(max_rle)
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      COMMON /rho_all/ xco_all(nxx,nxx,nch,nit),rmo_all(nxx,nn,nit)  
      write (*,*) 
      write (*,*) 'Check energies'
      write (*,*) 

      do 1000 iii=1,max_rle
      res=0.d0
      do 1 i=1,ncore 
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i)
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
                if (i.ge.j) then
                 index1=ipx1(i,j,indm,indn,k)
                 index2=ipxc(i,j,indm,indn,k)                 
                else
                 index1=ipx1(j,i,indn,indm,k)
                 index2=ipxc(j,i,indn,indm,k)                 
                endif                
                zz=0.5/(2*k+1.d0)
                do 3 m=nnm,nmaxx(indm)
                  do 6 n=nnn,nmaxx(indn)
                     if (i.ge.j) then
                      res=res+zz*xco_all(m,n,index2,iii)*x1(m,n,index1)        
                     else
                      res=res+zz*xco_all(n,m,index2,iii)*x1(n,m,index1)       
                     endif                     
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
                if (i.ge.j) then
                 index1=ipx1(i,j,indm,indn,k)                
                 index2=ipxc(i,j,indn,indm,k1)                 
                else                          
                 index1=ipx1(j,i,indn,indm,k)                
                 index2=ipxc(j,i,indm,indn,k1)                 
                endif                
                zz=0.5*d6j(jm,ja,2*k,jn,jb,2*k1)
                do 31 m=nnm,nmaxx(indm)
                  do 61 n=nnn,nmaxx(indn)
                     if (i.ge.j) then
                      res=res+zz*xco_all(n,m,index2,iii)*x1(m,n,index1)
                     else
                      res=res+zz*xco_all(m,n,index2,iii)*x1(n,m,index1)
                     endif                     
61                 continue 
31               continue
515             continue
21            continue
555             continue
11            continue
5           continue
2         continue
4       continue
1     continue
      write (*,'(i4,e16.8)') iii,res
1000  continue
      write (*,*) 'DONE'
901   continue
      return
      end


      subroutine outrho 
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)  

      OPEN (unit=7,file='pair.3',form='unformatted')
      jlo=1
      jhi=ncore

      write (7) nmax,lmax,icchan,jlo,jhi
      write (7) ipa,ipb,ipm,ipn,ipl
      write (7) ipxc
      
      do 27 i=1,icchan
        indn=ipn(i)
        indm=ipm(i)
        nmaxn=nmaxx(indn)
        nmaxm=nmaxx(indm)


        do 37 n=1,nmaxn
          write (7) (xco(m,n,i),m=1,nmaxm)
37       continue
27     continue
       do 47 i=jlo,jhi
          write (7) (rmo(n,i),n=1,nxx)
47     continue
       close (7)

        return
        end
c


      subroutine rhocore
      implicit double precision (a-h,o-z)
      
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2   
      DIMENSION ga(nhf),fa(nhf),gb(nhf),fb(nhf)
      DIMENSION u(nhf),v(nhf),vv(nhf),v2(nhf)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /x12ini/ inmax
      common /nmaxx/nmaxx(nx)

      lin=2*lmax+1
      if (lin.lt.inmax) then
       lin=inmax
      endif
      write (*,*) 'lin =',lin    

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
c                    if (m.eq.7.and.n.eq.7) then 
c              write (*,'(4i3,10e13.3)') i,j,ka,kb,v(100),u(100),ga(100),
c     *gg(100,m,indm),rp(100),v2(100),vv(100),c,x1(m,n,index),av
c                 endif
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
 
        endif
4       continue
         write (*,187) index,i,j,km,kn,k,av
1     continue
      icmax=index
      write (*,*) index

      jlo=1
      jhi=ncore
      ichan1=index
      ichan2=index1
      write (*,*) index,index1

901   continue
      return
      end



      subroutine rhoini(key)
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2  
 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)      
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      DIMENSION ipa9(NCH),ipb9(NCH),ipm9(NCH),ipn9(NCH),ipl9(NCH)
      DIMENSION ipxc9(nn,nn,nk,nk,0:kk)
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)

      do 531 i=1,nch
        do 532 m=1,nxx
          do 533 n=1,nxx
            xco(m,n,i)=0.d0
533       continue
532     continue
531   continue
      index=0
      do 1 i=1,ncore 
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i)
        do 80 m=1,nxx
         rmo(m,i)=0.d0
 80     continue        
        do 4 j=1,ncore  
          kb=ko(j)
          nb=no(j)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          eb=eo(j)          
          if (i.ge.j) then
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
                call odd (lm+la,i1)
                call odd (ln+lb,i2)
                if ((i1.eq.0.and.i2.ne.0).or.
     *              (i1.ne.0.and.i2.eq.0)) goto 551
     
                index=index+1
                ipxc(i,j,indm,indn,k)=index 
                ipa(index)=i
                ipb(index)=j
                ipm(index)=indm
                ipn(index)=indn
                ipl(index)=k

                do 732 m=1,nxx
                  do 733 n=1,nxx
                    xco(m,n,index)=0.d0
733               continue
732             continue

                call odd (lm+la+k,i1)
                call odd (ln+lb+k,i2)
                 
                if (i1*i2.eq.0) goto 555
                index1=ipx1(i,j,indm,indn,k)
                av=0
                do 3 m=nnm,nmaxx(indm)
                  do 6 n=nnn,nmaxx(indn)
                    eu=1.d0/(ea+eb-ee(m,indm)-ee(n,indn))
                    xco(m,n,index)=x1(m,n,index1)*eu
                    av=av+xco(m,n,index)*xco(m,n,index)
6                 continue 
3               continue
555             continue
551             continue
11            continue
5           continue
2         continue
187       format ('CHANNEL',i6,'  j1=',i3,'  j2=',i3,'  km=',i3,
     *    '  kn=',i3,'  k=',i3,f12.6) 
          write (*,187) index,i,j,km,kn,k,av
        endif
4       continue
1     continue
      icchan=index
      write (*,*) 'Number of channels = ',icchan
      if (key.eq.1) then
      OPEN (unit=7,file='pair.3',form='unformatted')
      read (7) nmax9,lmax9,icchan9,jlo9,jhi9
      read (7) ipa9,ipb9,ipm9,ipn9,ipl9
      read (7) ipxc9

      do 27 i=1,icchan9
          ii=ipa9(i)
          j=ipb9(i)
          indm=ipm9(i)
          indn=ipn9(i)
          k=ipl9(i)
          nmaxm=nmaxx(indm)
          nmaxn=nmaxx(indn)
          index=ipxc(ii,j,indm,indn,k)
        do 37 n=1,nmaxn
          read (7) (xci(m,n,i),m=1,nmaxm)
          do 77 m=1,nmaxm
            xco(m,n,index)=xci(m,n,i)
77        continue
37       continue
27     continue
      do 47 i=jlo9,jhi9
          read (7) (rmo(n,i),n=1,nxx)
47     continue
      close (7)

      do 39 i=1,icchan
        do 49 n=1,nxx
          do 59 m=1,nxx
            xci(m,n,i)=0.d0
59         continue
49       continue
39     continue
       endif

901   continue
      return
      end

      subroutine energy(res)
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)  
          
      res=0.d0
      do 1 i=1,ncore 
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        ea=eo(i)
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
                if (i.ge.j) then
                 index1=ipx1(i,j,indm,indn,k)
                 index2=ipxc(i,j,indm,indn,k)                 
                else
                 index1=ipx1(j,i,indn,indm,k)
                 index2=ipxc(j,i,indn,indm,k)                 
                endif                
                zz=0.5/(2*k+1.d0)
                do 3 m=nnm,nmaxx(indm)
                  do 6 n=nnn,nmaxx(indn)
                     if (i.ge.j) then
                      res=res+zz*xco(m,n,index2)*x1(m,n,index1)        
                     else
                      res=res+zz*xco(n,m,index2)*x1(n,m,index1)       
                     endif 

c                  if (m.eq.7.and.n.eq.7) then 
c               write (*,*) '1'
c                write (*,*) zz,xco(m,n,index2),x1(m,n,index1)        
c                  endif

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
                if (i.ge.j) then
                 index1=ipx1(i,j,indm,indn,k)                
                 index2=ipxc(i,j,indn,indm,k1)                 
                else                          
                 index1=ipx1(j,i,indn,indm,k)                
                 index2=ipxc(j,i,indm,indn,k1)                 
                endif                
                zz=0.5*d6j(jm,ja,2*k,jn,jb,2*k1)
                do 31 m=nnm,nmaxx(indm)
                  do 61 n=nnn,nmaxx(indn)
                     if (i.ge.j) then
                      res=res+zz*xco(n,m,index2)*x1(m,n,index1)
                     else
                      res=res+zz*xco(m,n,index2)*x1(n,m,index1)
                     endif         
c                   if (m.eq.7.and.n.eq.7) then 
c               write (*,*) '2'
c                write (*,*) zz,xco(n,m,index2),x1(m,n,index1)  
c                  endif
           
61                 continue 
31               continue

515             continue
21            continue
555             continue
11            continue
5           continue
2         continue
4       continue
1     continue
901   continue
      return
      end


      subroutine outin
      implicit double precision (a-h,o-z)
      
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      common /nmaxx/nmaxx(nx)

      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)  
187   format ('CHANNEL',i6,2e18.6) 

      do 1 i=1,ncore 
        do 2 m=1,nxx
          rmi(m,i)=rmo(m,i)
          rmo(m,i)=0.d0
2       continue
1     continue
      do 3 i=1,icchan
c        write (*,187) i,xci(6,7,i),xco(6,7,i)
        do 4 n=1,nxx
          do 5 m=1,nxx
            xci(m,n,i)=xco(m,n,i)
            xco(m,n,i)=0.d0
5         continue
4       continue
c        write (*,187) i,xci(6,7,i),xco(6,7,i)
c      write (*,*) 
3     continue
      return
      end
      
      subroutine term1
      implicit double precision (a-h,o-z)
      
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      common /nmaxx/nmaxx(nx)

      DIMENSION x(nx,nx)
      do 30 index=1,icchan
        i=ipa(index)
        j=ipb(index)
        indm=ipm(index)
        indn=ipn(index)
        k=ipl(index)
        ka=ko(i)
        na=no(i)
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
***************************************************
            l1min=max0(iabs(kapc-kapa),iabs(kapd-kapb))
            l1max=min0((kapc+kapa-1),(kapd+kapb-1))
            DO 5 l=l1min,l1max
              call odd (la+lc+l,i1)
              call odd (lb+ld+l,i2)
              if (i1*i2.eq.0) goto 515
              index1=ipx1(i,j,indc,indd,l)
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
                    if (ii.ge.jj) then
         xco(m,n,index)=xco(m,n,index)+z3*x1(nc-n0c,nd-n0d,index1)*
     *                  xci(m,n,index2)
                    else
         xco(m,n,index)=xco(m,n,index)+z3*x1(nc-n0c,nd-n0d,index1)*
     *                  xci(n,m,index2)
                    endif
11                continue
10              continue
535             continue
9             continue
515           continue
5           continue
4         continue
3       continue
187       format ('CHANNEL',i6,'  i=',i3,'  j=',i3,'  km=',i3,
     *    '  kn=',i3,'  k=',i3,2e18.6) 
c          write (*,187) index,i,j,km,kn,k,xco(6,7,index),xco(7,6,index)
30    continue
      return
      end


      subroutine outt(deltae,damp,iter,key)
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn) 
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
	common /nmaxx/nmaxx(nx)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1) 
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1  
c      write (*,*) icchan
      if (key.eq.0) then
       del=0.d0
      else
       del=deltae
      endif
      do 100 i=1,ncore
        ka=ko(i)
        ea=eo(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        indm=inda
        CALL indk1(indm,km,kapm,lm,jm,m0)
        call st (km,nnm)
        do 200 m=nnm,nmaxx(indm)
             eu=1.d0/(del+ea-ee(m,indm))
             rmo(m,i)=rmo(m,i)*eu
c          write (*,'(i3,e18.8)') i,rmo(6,i)
200     continue
100   continue
********** DAMP *************************
       if (iter.ge.3) then 
      do 110 i=1,ncore
        ka=ko(i)
        ea=eo(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        indm=inda
        CALL indk1(indm,km,kapm,lm,jm,m0)
        call st (km,nnm)

       do 220 m=nnm,nmaxx(indm)
             dd1=damp/(1.d0+damp)
             dd2=1.d0/(1.d0+damp)
             rmo(m,i)=dd1*rmi(m,i)+dd2*rmo(m,i)
220      continue
110   continue

       endif
******************************************
    
      do 3 iu=1,icchan
        i=ipa(iu)
        j=ipb(iu)
        indm=ipm(iu)
        indn=ipn(iu)
        k=ipl(iu)
        ka=ko(i)
        kb=ko(j)
        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        CALL klj(kb,kapb,lb,jb,indb,n0b)
        call st (km,nnm)            
        call st (kn,nnn) 
        ea=eo(i)
        eb=eo(j)
        index=ipx1(i,j,indm,indn,k)
c        write (*,187) iu,i,j,km,kn,k,xco(6,7,iu),xco(7,6,iu)
        call odd (la+lm+k,i1)
        call odd (lb+ln+k,i2)
        if (i1*i2.eq.0) then
         do 10 n=nnn,nmaxx(indn)
           do 11 m=nnm,nmaxx(indm)
             eu=1.d0/(del+ea+eb-ee(m,indm)-ee(n,indn))
             xco(m,n,iu)=xco(m,n,iu)*eu
11         continue
10       continue
        else
         do 101 n=nnn,nmaxx(indn)
           do 111 m=nnm,nmaxx(indm)
             eu=1.d0/(del+ea+eb-ee(m,indm)-ee(n,indn))
             xco(m,n,iu)=(x1(m,n,index)+xco(m,n,iu))*eu
111        continue
101      continue
        endif
3     continue
187   format ('CHANNEL',i6,'  i=',i3,'  j=',i3,'  km=',i3,
     *'  kn=',i3,'  k=',i3,2e18.6) 
c      write (*,*) index

********** DAMP *************************
       if (iter.ge.3) then 

      do 33 iu=1,icchan

        indm=ipm(iu)
        indn=ipn(iu)
    
        CALL indk1(indm,km,kapm,lm,jm,m0)
        CALL indk1(indn,kn,kapn,ln,jn,n0)
       
        call st(km,nnm)            
        call st(kn,nnn) 

       do 230 m=nnm,nmaxx(indm)
          do 240 n=nnn,nmaxx(indn)
 
             dd1=damp/(1.d0+damp)
             dd2=1.d0/(1.d0+damp)
        
         xco(m,n,iu)=dd1*xci(m,n,iu)+dd2*xco(m,n,iu)
240      continue
230      continue
33     continue
       endif

******************************************
      return
      end
      

      subroutine term5(key)
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      do 30 index=1,icchan
        if (key.eq.1) then
         i=ipa(index)
         j=ipb(index)
         indm=ipm(index)
         indn=ipn(index)
        endif
        if (key.eq.2) then
         j=ipa(index)
         i=ipb(index)
         indn=ipm(index)
         indm=ipn(index)
        endif
        k=ipl(index)
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
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
          CALL klj(kc,kapc,lc,jc,indc,n0c)
          call odd (la+lc+k,i1)
          call odd (lb+ln+k,i2)
          if (kc.eq.km.and.(i1*i2).ne.0) then
           if (i.ge.j) then
            iu=ipx1(i,j,indc,indn,k)
           else
            iu=ipx1(j,i,indn,indc,k)
           endif
c           write (*,*) i,j,indc,indn,k,iu,index
           if (key.eq.1) then
           do 4 n=nnn,nmaxx(indn)
             do 41 m=nnm,nmaxx(indm)
               if (i.ge.j) then
                xco(m,n,index)=xco(m,n,index)-x1(nc-n0c,n,iu)*rmi(m,ii)        
               else
                xco(m,n,index)=xco(m,n,index)-x1(n,nc-n0c,iu)*rmi(m,ii)
               endif
41           continue
4          continue
           endif
           if (key.eq.2) then
           do 49 m=nnm,nmaxx(indm)
             do 419 n=nnn,nmaxx(indn)
               if (i.ge.j) then
                xco(n,m,index)=xco(n,m,index)-x1(nc-n0c,n,iu)*rmi(m,ii)
               else
                xco(n,m,index)=xco(n,m,index)-x1(n,nc-n0c,iu)*rmi(m,ii)
               endif
419           continue
49          continue
           endif
          endif
40      continue 
30    continue
901   continue
      return
      end

      subroutine terma1
      implicit double precision (a-h,o-z)
      
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1   
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1     
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      do 30 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)   
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
 
          if (i.ge.j) then
           index1=ipx1(i,j,indm,indn,k)                            
          else                          
           index1=ipx1(j,i,indn,indm,k)                             
          endif  
          z1=((-1)**(kapb-kapn))*SQRT((jb+1.d0)/(ja+1.d0))
c          write (*,*) i,j,km,kn,km,ka,kindex1
          do 4 n=nnn,nmaxx(indn)
            do 15 m=nnm,nmaxx(indm)
c        write (*,*) m,n
              if (i.ge.j) then   
                rmo(m,i)=rmo(m,i)+z1*x1(m,n,index1)*rmi(n,j)
               else
                rmo(m,i)=rmo(m,i)+z1*x1(n,m,index1)*rmi(n,j)
              endif
15          continue
4         continue 
700        continue
         endif
        endif
5       continue

2      continue
40      continue
c      write (*,'(i5,e23.10)') i,rmo(6,i)
30    continue
      return
      end


      subroutine terma2
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      do 30 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)   
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
          zz=(-1)**(kapn+kapm+kapa+kapb)
          DO 111 k=kmin,kmax 
            call odd (lm+ln+k,i1)
            call odd (la+lb+k,i2)
            if (i1*i2.eq.0) goto 100
            if (j.ge.i) then
             index1=ipx2(j,i,indm,indn,k)                            
            else                          
             index1=ipx2(i,j,indn,indm,k)                             
            endif  
            z1=((-1)**(kapa+kapb+k-1))/(ja+1.d0)
            do 4 n=nnn,nmaxx(indn)
              do 15 m=nnm,nmaxx(indm)
                if (j.ge.i) then   
                 rmo(m,i)=rmo(m,i)+z1*x2(m,n,index1)*rmi(n,j)
                else
                 rmo(m,i)=rmo(m,i)+z1*zz*x2(n,m,index1)*rmi(n,j)
                endif
15            continue
4           continue 
100         continue
111       continue
        endif
5       continue

2      continue

40      continue
c     write (*,'(i5,e23.10)') i,rmo(6,i)
30    continue
      return
      end

      subroutine terma3
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      do 30 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 40 j=1,ncore
          kb=ko(j)
          nb=no(j)   
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 25 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)
            if (indm.eq.inda) then
          do 3 ii=1,ncore 
            kc=ko(ii)
            nc=no(ii)
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
                    if (ii.ge.j) then   
                     rmo(m,i)=rmo(m,i)+z1*x1(n,na-n0a,index1)*
     *                                   xci(n,m,index2)
                    else
                     rmo(m,i)=rmo(m,i)+z1*x1(na-n0a,n,index1)*
     *                                   xci(m,n,index2)
                    endif
15                continue
4               continue 

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
                    if (ii.ge.j) then   
                     rmo(m,i)=rmo(m,i)+zz*x1(n,na-n0a,index1)*
     *                                   xci(m,n,index2)
                    else
                     rmo(m,i)=rmo(m,i)+zz*x1(na-n0a,n,index1)*
     *                                   xci(n,m,index2)
                    endif
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
30    continue
      return
      end


      subroutine term41
      implicit double precision (a-h,o-z)
      
      include "all.par"
      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2   
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      common /nmaxx/nmaxx(nx)

      DIMENSION x(nx,nx),y(nx,nx),z(nx,nx)
      do 32 i=1,ncore 
        ka=ko(i)
        na=no(i)
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
172             continue
171           continue

              if (i.ge.ii) then
               index2=ipxc(i,ii,indm,indr,k)
              else
               index2=ipxc(ii,i,indr,indm,k)
              endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 m=nnm,nmaxx(indm)
                  if (i.ge.ii) then
                   y(m,nr)=y(m,nr)+xci(m,nr,index2)
                  else
                   y(m,nr)=y(m,nr)+xci(nr,m,index2)
                 endif
10              continue
12            continue

              k1min=max0(iabs(kapm-kapc),iabs(kapr-kapa))
              k1max=min0((kapm+kapc-1),(kapr+kapa-1))
              DO 9 k1=k1min,k1max
                z2=(2*k+1)*d6j(jm,ja,2*k,jr,jc,2*k1)
                if (i.ge.ii) then
                 index2=ipxc(i,ii,indr,indm,k1)
                else
                 index2=ipxc(ii,i,indm,indr,k1)
                endif

                do 121 nr=nnr,nmaxx(indr)
                  do 101 m=nnm,nmaxx(indm)
                    if (i.ge.ii) then
                     y(m,nr)=y(m,nr)+z2*xci(nr,m,index2)
                    else
                     y(m,nr)=y(m,nr)+z2*xci(m,nr,index2)
                    endif
101                continue
121              continue
9              continue
               do 33 j=1,ncore 
                 kb=ko(j)
                 nb=no(j)
                 CALL klj(kb,kapb,lb,jb,indb,n0b)
                 if (i.ge.j) then
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
                   index=ipxc(i,j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
                        do 107 m=nnm,nmaxx(indm)
                         xco(m,n,index)=xco(m,n,index)+
     *                                  x(nr,n)*y(m,nr)

107                     continue
117                   continue
127                 continue
700              continue
23               continue
                endif
33             continue
900           continue
500           continue
21          continue
22        continue
31      continue
32    continue
      return
      end

      subroutine term42
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /xcore1/ ipx1(nn,nn,nk,nk,0:kk),x1(nxx,nxx,nxk),ichan1 
      COMMON /xcore2/ ipx2(nn,nn,nk,nk,0:kk),x2(nxx,nxx,nxk),ichan2   
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      DIMENSION x(nx,nx),y(nx,nx),z(nx,nx)
      do 32 j=1,ncore 
        kb=ko(j)
        nb=no(j)
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
172             continue
171           continue

              if (j.ge.ii) then
               index2=ipxc(j,ii,indn,indr,k)
              else
               index2=ipxc(ii,j,indr,indn,k)
              endif

              do 12 nr=nnr,nmaxx(indr)
                do 10 n=nnn,nmaxx(indn)
                  if (j.ge.ii) then
                   y(n,nr)=y(n,nr)+xci(n,nr,index2)
                  else
                   y(n,nr)=y(n,nr)+xci(nr,n,index2)
                 endif
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
                    if (j.ge.ii) then
                     y(n,nr)=y(n,nr)+z2*xci(nr,n,index2)
                    else
                     y(n,nr)=y(n,nr)+z2*xci(n,nr,index2)
                    endif
101                continue
121              continue
9              continue
               do 33 i=1,ncore 
                 ka=ko(i)
                 na=no(i)
                 CALL klj(ka,kapa,la,ja,inda,n0a)
                 if (i.ge.j) then
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

                   if (ii.ge.i) then
                    index1=ipx1(ii,i,indr,indm,k)
                   else
                    index1=ipx1(i,ii,indm,indr,k)
                   endif
                   z2=zf*((-1)**(kapc-kapr))
                   do 123 m=nnm,nmaxx(indm)
                     do 103 nr=nnr,nmaxx(indr)
                       if (ii.ge.i) then
                        x(nr,m)=x(nr,m)+z2*x1(nr,m,index1)
                       else
                        x(nr,m)=x(nr,m)+z2*x1(m,nr,index1)
                       endif
103                  continue
123                continue
800                continue

                   zz=(-1)**(kapa+kapc+kapr+kapm)
                   k2min=max0(iabs(kapc-kapa),iabs(kapm-kapr))
                   k2max=min0((kapc+kapa-1),(kapm+kapr-1))
                   DO 119 k2=k2min,k2max
                     call odd (lc+la+k2,i1)
                     call odd (lm+lr+k2,i2)
                     if (i1*i2.eq.0) goto 155

                     if (ii.ge.i) then
                      index1=ipx2(ii,i,indm,indr,k2)
                     else
                      index1=ipx2(i,ii,indr,indm,k2)
                     endif
                     z1=zf*(2*k+1.d0)*d6j(ja,jm,2*k,jr,jc,2*k2)
                     do 124 m=nnm,nmaxx(indm)
                       do 104 nr=nnr,nmaxx(indr)
                         if (ii.ge.i) then
                          x(nr,m)=x(nr,m)+z1*x2(m,nr,index1)
                         else
                          x(nr,m)=x(nr,m)+z1*zz*x2(nr,m,index1)
                         endif
104                    continue
124                  continue
155                  continue
119                continue
                   index=ipxc(i,j,indm,indn,k)


                    do 127 nr=nnr,nmaxx(indr)
                      do 117 n=nnn,nmaxx(indn)
                        do 107 m=nnm,nmaxx(indm)
                         xco(m,n,index)=xco(m,n,index)+
     *                                  x(nr,m)*y(n,nr)
107                     continue
117                   continue
127                 continue
700              continue
23               continue
                endif
33             continue
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
      common /nmaxx/nmaxx(nx)

      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
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
            do 30 i=1,ncore
              ka=ko(i)
              na=no(i)
              CALL klj(ka,kapa,la,ja,inda,n0a)
              DO 15 il=1,max
                ga(il)=go(il,i)
                fa(il)=fo(il,i)
15            CONTINUE    
              if (i.ge.j) then
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
                  v1(ij)=v1(ij)+gg(ij,nr,indr)*rmi(nr,i)
                  v2(ij)=v2(ij)+ff(ij,nr,indr)*rmi(nr,i) 
11              continue
10            continue 
              do 25 indm=1,2*lmax+1
                CALL indk1(indm,km,kapm,lm,jm,m0)
                call st (km,nnm)    
                call trgi(iabs(kapm-kapa),(kapm+kapa-1),k,i1) 
                call odd (la+lm+k,i2)
                if (i1*i2.eq.0) goto 800
                c=((-1)**(k))*s(k,km,kr)*s(k,kn,kb)
                index=ipxc(i,j,indm,indn,k)
                do 165 n=nnn,nmaxx(indn)
                  do 160 m=nnm,nmaxx(indm)
                    do 161 ij=1,max 
                      vv(ij)=(v1(ij)*gg(ij,m,indm)+
     *                        v2(ij)*ff(ij,m,indm))*uu(ij,n)*rp(ij) 
161                 continue
             xco(m,n,index)=xco(m,n,index)+c*rint(vv,1,max,11,h)
160               continue
165             continue 
800             continue
25            continue
             endif
26          continue
           endif
30          continue 
900         continue
500       continue
23      continue
32    continue 
901   continue
      return
      end

      subroutine term32
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION vv(nhf),u(nhf),v(nhf),fa(nhf),fb(nhf),ga(nhf),gb(nhf)
      DIMENSION v1(nhf),v2(nhf),uu(nhf,nx)
      do 32 i=1,ncore 
        ka=ko(i)
        na=no(i)   
        CALL klj(ka,kapa,la,ja,inda,n0a)
        DO 13 il=1,max
          ga(il)=go(il,i)
          fa(il)=fo(il,i)
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
              if (i.ge.j) then
          do 26 indr=1,2*lmax+1

              if (indr.eq.indb) then
              CALL indk1(indr,kr,kapr,lr,jr,n0r)
              call st (kr,nnr)

                  do 33 ij=1,max 
                  v1(ij)=0.d0
                  v2(ij)=0.d0
33              continue
              do 10 nr=nnr,nmaxx(indr)
                do 11 ij=1,max 
                  v1(ij)=v1(ij)+gg(ij,nr,indr)*rmi(nr,j)
                  v2(ij)=v2(ij)+ff(ij,nr,indr)*rmi(nr,j) 
11              continue
10            continue 
              do 25 indn=1,2*lmax+1
                CALL indk1(indn,kn,kapn,ln,jn,n0)
                call st (kn,nnn)    
                call trgi(iabs(kapn-kapr),(kapn+kapr-1),k,i1) 
                call odd (lr+ln+k,i2)
                if (i1*i2.eq.0) goto 800
                c=((-1)**(k))*s(k,kn,kr)*s(k,km,ka)
                index=ipxc(i,j,indm,indn,k)
                do 165 n=nnn,nmaxx(indn)
                  do 160 m=nnm,nmaxx(indm)
                    do 161 ij=1,max 
                      vv(ij)=(v1(ij)*gg(ij,n,indn)+
     *                        v2(ij)*ff(ij,n,indn))*uu(ij,m)*rp(ij) 
161                 continue
             xco(m,n,index)=xco(m,n,index)+c*rint(vv,1,max,11,h)
160               continue
165             continue 
800             continue
25            continue
             endif
26          continue
           endif
30          continue 
900         continue
500       continue
23      continue
32    continue 
901   continue
      return
      end

      subroutine terma4
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      common /ips/ ipa(NCH),ipb(NCH),ipm(NCH),ipn(NCH),ipl(NCH)
      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)
      DIMENSION vv(nhf),u(nhf),v(nhf),fa(nhf),fb(nhf),ga(nhf),gb(nhf)
      DIMENSION v1(nhf),v2(nhf),uu(nhf,nx),x(nx,nx)
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
            do 30 i=1,ncore
              ka=ko(i)
              na=no(i)
              CALL klj(ka,kapa,la,ja,inda,n0a)
              DO 15 il=1,max
                ga(il)=go(il,i)
                fa(il)=fo(il,i)
15            CONTINUE    
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

                if (i.ge.j) then 
                 index2=ipxc(i,j,indr,indn,k)
                else
                 index2=ipxc(j,i,indn,indr,k)
                endif
            zz=((-1)**(kapa+kapr+kapb+kapn))/((ja+1.d0)*(2*k+1.d0))
                do 165 n=nnn,nmaxx(indn)
                  do 160 nr=nnr,nmaxx(indr)

                  if (i.ge.j) then
                    x(nr,n)=xci(nr,n,index2)
                   else
                    x(nr,n)=xci(n,nr,index2)
                  endif
160               continue
165             continue
                k2min=max0(iabs(kapr-kapb),iabs(kapn-kapa))
                k2max=min0((kapr+kapb-1),(kapn+kapa-1))
                DO 119 k2=k2min,k2max
                  call odd (lr+lb,i5)
                  call odd (ln+la,i6)
                   if ((i5.eq.0.and.i6.ne.0).or.
     *                (i5.ne.0.and.i6.eq.0)) goto 155


                  if (i.ge.j) then 
                   index2=ipxc(i,j,indn,indr,k2)
                  else
                   index2=ipxc(j,i,indr,indn,k2)
                  endif
                  z1=(2*k+1.d0)*d6j(jr,ja,2*k,jn,jb,2*k2)
                  do 185 n=nnn,nmaxx(indn)
                    do 180 nr=nnr,nmaxx(indr)
                      if (i.ge.j) then
                       x(nr,n)=x(nr,n)+z1*xci(n,nr,index2)
                      else
                       x(nr,n)=x(nr,n)+z1*xci(nr,n,index2)
                      endif
180                 continue
185               continue
155               continue
119               continue
                  c=((-1)**(k))*s(k,km,kr)*s(k,kb,kn)

                  do 33 ij=1,max 
                    v1(ij)=0.d0
                    v2(ij)=0.d0
33                continue
                 do 109 n=nnn,nmaxx(indn)
                   do 10 nr=nnr,nmaxx(indr)
                     do 11 ij=1,max 
                       v1(ij)=v1(ij)+gg(ij,nr,indr)*x(nr,n)*uu(ij,n)
                       v2(ij)=v2(ij)+ff(ij,nr,indr)*x(nr,n)*uu(ij,n) 
11                   continue
10                 continue 
109              continue
                 do 169 m=nnm,nmaxx(indm)
                   do 161 ij=1,max 
                     vv(ij)=(v1(ij)*gg(ij,m,indm)+
     *                       v2(ij)*ff(ij,m,indm))*rp(ij) 
161                continue
                   rmo(m,i)=rmo(m,i)+c*zz*rint(vv,1,max,11,h)
169              continue
800              continue 
24             continue
             endif
25         continue
30           continue
900          continue
500        continue
23      continue
32    continue 
901   continue
      return
      end

      subroutine term2
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)     
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      DIMENSION v(nhf),u(nhf),uu(nhf,nx,nx),v1(nhf,nx),v2(nhf,nx),
     *vv(nhf),uv1(nhf,nx),uv2(nhf,nx)
      DIMENSION uu1(nhf,nxx,0:kk,kk),uu2(nhf,nxx,0:kk,kk)
      do 5 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 6 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          if (i.ge.j) then
          write (*,'(2i5)') i,j
          t1=mclock()
          do 1 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 91 indn=1,kk
              DO 92 k=0,kk
                do 93 m=1,nxx
                  do 94 ij=1,nhf
                    uu1(ij,m,k,indn)=0.d0
                    uu2(ij,m,k,indn)=0.d0
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
                      index1=ipxc(i,j,indr,inds,k1)
                      do 31 nr=nnr,nmaxx(indr)
                        do 203 ij=1,max
                          v1(ij,nr)=0.d0
                          v2(ij,nr)=0.d0
203                     continue
                        do 41 nss=nns,nmaxx(inds)
                          do 201 ij=1,max
                v1(ij,nr)=v1(ij,nr)+gg(ij,nss,inds)*xci(nr,nss,index1)
                v2(ij,nr)=v2(ij,nr)+ff(ij,nss,inds)*xci(nr,nss,index1)
201                       continue
41                      continue
31                    continue

                      do 51 m=nnm,nmaxx(indm)
                        do 206 ij=1,max
                          uv1(ij,m)=0.d0
                          uv2(ij,m)=0.d0
206                     continue
                        do 57 nr=nnr,nmaxx(indr)
                          do 202 ij=1,max
                            uv1(ij,m)=uv1(ij,m)+v1(ij,nr)*uu(ij,nr,m)
                            uv2(ij,m)=uv2(ij,m)+v2(ij,nr)*uu(ij,nr,m)
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
                        DO 219 k=kmin,kmax
                          c1=(2*k+1.d0)*d6j(jm,ja,2*k,2*k1,2*l,jr)  
                          cfin=c*c1*z1*d6j(jn,jb,2*k,2*k1,2*l,js)
                          do 59 m=nnm,nmaxx(indm)
                            do 280 ij=1,max
               uu1(ij,m,k,indn)=uu1(ij,m,k,indn)+uv1(ij,m)*cfin
               uu2(ij,m,k,indn)=uu2(ij,m,k,indn)+uv2(ij,m)*cfin
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
                  index=ipxc(i,j,indm,indn,k)
                  do 52 n=nnn,nmaxx(indn)
                    do 53 m=nnm,nmaxx(indm)
                      do 205 ij=1,max
                        vv(ij)=(uu1(ij,m,k,indn)*gg(ij,n,indn)+
     *                          uu2(ij,m,k,indn)*ff(ij,n,indn))*rp(ij)
205                   continue
                      xco(m,n,index)=xco(m,n,index)+
     *                rint(vv,1,max,11,h)
53                  continue
52                continue
319             continue
810             continue
189           continue
1           continue
            t2=mclock()
            dt=(t2-t1)/100.d0
c            write (*,'(i3,f10.2)') indm,dt
           endif
6        continue
5      continue
901   continue
      return
      end

 

      subroutine term21
      implicit double precision (a-h,o-z)
      
      include "all.par"
      common /nmaxx/nmaxx(nx)

      common /spl/ gg(nhf,nx,kx1),ff(nhf,nx,kx1),ee(nx,kx1)     
      COMMON /radial/ r(nhf),rp(nhf),rpor(nhf),h,max
      COMMON /ldata/ nmax,lmax,nmax1,lmax1       
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore  
      COMMON /irhoc/ ipxc(nn,nn,nk,nk,0:kk),icchan
      COMMON /rhocin/ xci(nxx,nxx,nch),rmi(nxx,nn)
      COMMON /rhocout/ xco(nxx,nxx,nch),rmo(nxx,nn)
      DIMENSION v(nhf),u(nhf),uu(nhf,nx,nx),v1(nhf,nx),v2(nhf,nx),
     *vv(nhf),uv1(nhf,nx),uv2(nhf,nx)
      DIMENSION uu1(nhf,nxx,0:kk,kk),uu2(nhf,nxx,0:kk,kk)
      common /iqqq/ iqa(nnh),iqb(nnh)
      iiq=0
      do 109 i1=1,ncore
        do 110 i2=1,ncore
          if (i1.ge.i2) then
          iiq=iiq+1
          iqa(iiq)=i1
          iqb(iiq)=i2
c          write (*,'(3i5)') i1,i2,iiq,nnh
         endif
110     continue
109   continue
      nnmq=iiq

      do 5 ittt=1,nnmq
      i=iqa(ittt)
      j=iqb(ittt)
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          do 1 indm=1,2*lmax+1
            CALL indk1(indm,km,kapm,lm,jm,m0)
            call st (km,nnm)

            do 91 indn=1,kk
              DO 92 k=0,kk
                do 93 m=1,nxx
                  do 94 ij=1,nhf
                    uu1(ij,m,k,indn)=0.d0
                    uu2(ij,m,k,indn)=0.d0
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
                      index1=ipxc(i,j,indr,inds,k1)
                      do 31 nr=nnr,nmaxx(indr)
                        do 203 ij=1,max
                          v1(ij,nr)=0.d0
                          v2(ij,nr)=0.d0
203                     continue
                        do 41 nss=nns,nmaxx(inds)
                          do 201 ij=1,max
                v1(ij,nr)=v1(ij,nr)+gg(ij,nss,inds)*xci(nr,nss,index1)
                v2(ij,nr)=v2(ij,nr)+ff(ij,nss,inds)*xci(nr,nss,index1)
201                       continue
41                      continue
31                    continue

                      do 51 m=nnm,nmaxx(indm)
                        do 206 ij=1,max
                          uv1(ij,m)=0.d0
                          uv2(ij,m)=0.d0
206                     continue
                        do 57 nr=nnr,nmaxx(indr)
                          do 202 ij=1,max
                            uv1(ij,m)=uv1(ij,m)+v1(ij,nr)*uu(ij,nr,m)
                            uv2(ij,m)=uv2(ij,m)+v2(ij,nr)*uu(ij,nr,m)
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
                        DO 219 k=kmin,kmax
                          c1=(2*k+1.d0)*d6j(jm,ja,2*k,2*k1,2*l,jr)  
                          cfin=c*c1*z1*d6j(jn,jb,2*k,2*k1,2*l,js)
                          do 59 m=nnm,nmaxx(indm)
                            do 280 ij=1,max
               uu1(ij,m,k,indn)=uu1(ij,m,k,indn)+uv1(ij,m)*cfin
               uu2(ij,m,k,indn)=uu2(ij,m,k,indn)+uv2(ij,m)*cfin
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
                  index=ipxc(i,j,indm,indn,k)
                  do 52 n=nnn,nmaxx(indn)
                    do 53 m=nnm,nmaxx(indm)
                      do 205 ij=1,max
                        vv(ij)=(uu1(ij,m,k,indn)*gg(ij,n,indn)+
     *                          uu2(ij,m,k,indn)*ff(ij,n,indn))*rp(ij)
205                   continue
                      xco(m,n,index)=xco(m,n,index)+
     *                rint(vv,1,max,11,h)
53                  continue
52                continue
319             continue
810             continue
189           continue
1           continue
5      continue
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
      lin=2*lmax+1
      if (lin.lt.inmax) then
       lin=inmax
      endif
      index=0
      index1=0

      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          if (i.ge.j) then

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
      if (kmax.gt.kk) then 
      write (*,*) ' Parameter KK is exceeded in rhocore -1 , STOP'
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
*             Calculate  Xk(manb) m,n can be core this time!!! *
****************************************************************
              kmin=max0(iabs(kapm-kapn),iabs(kapa-kapb))
              kmax=min0((kapm+kapn-1),(kapa+kapb-1))
      if (kmax.gt.kk) then 
	write (*,*) ' Parameter KK is exceeded in rhocore -2, STOP'
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
        endif
4       continue
1     continue
**************************************************
      write (*,*)
      write (*,*) 'PRE-COUNTING ALL CORE CHANNELS'
	write (*,*)
	if (index.le.NXK.and.index1.le.NXK) then 
	 write (*,122) index,nxk
	 write (*,123) index1,nxk
      else
	 write (*,124) index,nxk
	 write (*,125) index1,nxk
       write (*,*) 'Parameter NXK is exceeded, STOP'
	 stop
	endif
122   format ('Number of X_k(mnab) channels      = ',i7,
     *' < NXK = ',i7,' OK')
123   format ('Number of X_k(manb) channels      = ',i7,
     *' < NXK = ',i7,' OK')

124   format ('Number of X_k(mnab) channels      = ',i7,
     *' > NXK = ',i7)
125   format ('Number of X_k(manb) channels      = ',i7,
     *' > NXK = ',i7)



      return
      end



      subroutine countrho
      implicit double precision (a-h,o-z)
	    
      include "all.par"
      COMMON /ldata/ nmax,lmax,nmax1,lmax1
      common /nmaxx/nmaxx(nx)
      COMMON /core/ go(nhf,ns),fo(nhf,ns),eo(ns),no(ns),ko(ns),ncore

      index=0
      do 1 i=1,ncore
        ka=ko(i)
        na=no(i)
        CALL klj(ka,kapa,la,ja,inda,n0a)
        do 4 j=1,ncore
          kb=ko(j)
          nb=no(j)
          CALL klj(kb,kapb,lb,jb,indb,n0b)
          if (i.ge.j) then
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
	  write (*,*) ' Parameter KK is exceeded in rhoini , STOP'
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
        endif
4       continue
1     continue
	if (index.le.Nch) then 
	 write (*,122) index,nch
      else
	 write (*,124) index,nch
       write (*,*) 'Parameter NCH is exceeded, STOP'
	 stop
	endif
122   format ('Number of rho_k(mnab) channels    = ',i7,
     *' < NCH = ',i7,' OK')

124   format ('Number of rho_k(mnab) channels    = ',i7,
     *' > NCH = ',i7)

      write (*,*)



      return
      end


c     ========================================
      include "d6j.f"
      include "libD.f"
      include "rint.f"
      include "yfun.f"
      include "yint.f"
      include "inidat.f"
      include "linpakd.f"