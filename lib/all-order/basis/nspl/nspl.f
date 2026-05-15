      PROGRAM DHFORB
************************************************
*     Driving program for creating HF orbitals
*         to be used in calculating second
*           order correlation corrections
*
*         multigrid version  8/24/87
*         input new grid parameters: 3/25/97
*
************************************************
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      parameter(NK=KX+NX)
      common/newgri/t1(NK),t2(NK),t3(NK)
 
      open(unit=1,form="unformatted",file="fort.1",status="old")
      read(5,*) lmax
      write(6,1001) lmax
 1001 format(/' lmax =',i3)
 5    read(5,*) rmax
      write(6,1002) rmax
 1002 format(' rmax = ',f16.5)
 
***********************************************
*       enter number of collocation points : n
*                       order of spline    : k
***********************************************
      read(5,*) n,k
      write(6,1003) n,k
 1003 format('  n = ',i4,10x,'k = ',i3)
      if(n.gt.NX .or. n.lt.(k+2)) go to 999
      if(k.gt.KX .or. k.lt.2) go to 999
 
***************************************************
*    find          knots              t(i)  i=1,n+k
****************************************************
 
      dr0 = 0.0001
      call mkgrid(n,k,rmax,dr0,t1,*5)
      write(6,1010) (t1(i),i=1,n+k)
 1010 format(/' knots : t1'/(1p,5d14.6))
      dr0 = 0.001
      call mkgrid(n,k,rmax,dr0,t2,*5)
      write(6,1011) (t2(i),i=1,n+k)
 1011 format(/' knots : t2'/(1p,5d14.6))
      dr0 = 0.01
      call mkgrid(n,k,rmax,dr0,t3,*5)
      write(6,1012) (t3(i),i=1,n+k)
 1012 format(/' knots : t3'/(1p,5d14.6))
 
********  call in core HF functions **********newgr
 
      call infun
 
*********  set up matrix for metric and derivs *****
 
      call putout(n,k,lmax)
 
 999  stop
      end
 
      subroutine putout(n,k,lmax)
      implicit doubleprecision(a-h,o-z)
******************************************************************
*      Routine to write spline data on a direct access file
*      for future use in calculating correlation energy
*      3/27/97  added call to newgrid for better norms
******************************************************************
      external dsygv
      include "global.par"
      parameter (NL=2*NX,NREC=8*NL*(NHF+1),NK=KX+NX)
      dimension cd(NL,NL),ev(NL)
      dimension t(NK)
      dimension g(NHF,NX),f(NHF,NX)
      common/radial/r(NHF),rp(NHF),rpor(NHF),h,max
      common/shells/wco(NS),nco(NS),kco(NS),mco(NS),jmax
      common/newgri/t1(NK),t2(NK),t3(NK)
      dimension rn(NHF),rpn(NHF),rporn(NHF)
      data index /0/
      common /rho0/ en,cn,tn
****  read new grid parameters, e.g.:  (0.2, 0.03125, 500) or (0,0,0)

      read(5,*) hp,h0,nmax
      IF(nmax.gt.NHF) THEN
         write(6,1200) 
 1200    format(' nmax =',i10,' > ',i10,' is too large')
         stop
      ELSEIF(nmax.eq.0) THEN
****   default nmax = max
         nmax = max
      END IF
            
      IF(h0.eq.0d0) THEN
****   default  h0 = h
         h0 = h
      END IF
      
      IF(hp.eq.0d0) THEN
****  default if hp = 0 use hf grid
         DO i = 1,max
           rn(i) = r(i)
           rpn(i) = rp(i)
           rporn(i) = rpor(i)
         END DO
         h0 = h
         nmax = max
         imax = 0        
         DO i = 1,max
            IF(r(i).le.t1(n+k)) THEN
               imax = i
            END IF
         END DO
         write(6,1210) imax,r(imax),h,rp(1)
 1210    format(' Grid: r_n = r_0 [ exp( h(n-1) ) - 1 ]'/
     a       '                    max =',i5/
     1       '                 r(max) =',1p,e15.7/
     2       '                      h =',e15.7/
     3       '                     r0 =',e15.7/)

      ELSE
      
****  if hp != 0 use modified version of Vladimir's grid

         rmax = t1(n+k)
         call newgrid(h0,hp,rmax,nmax,r0,rn,rpn,rporn)
         IF((r0. le. 0d0) .or. (r0 .ge. 0.1d0)) THEN
           write(6,1220) r0
 1220      format(' r0 =',1p,e12.4,'  improper; try again!')
           stop
         END IF
         a = h0/hp
         write(6,1230) nmax,rmax,h0,r0,hp,a
 1230    format(' Grid: log(r_n/r_0+1) + a r_n = h(n-1); a = h/hp' /
     a       '                    max =',i5/
     1       '                 r(max) =',1p,e15.7/
     2       '                      h =',e15.7/
     3       '                     r0 =',e15.7/
     4       '                     hp =',e15.7/
     5       '                     a  =',e15.7/)

      END IF

      if(index.eq.0) then
         open(unit=3,file='hfspl.1',form='unformatted')
         write(3) n,k,t1,lmax
         write(3) jmax,(wco(j),nco(j),kco(j),mco(j),j=1,jmax)
         write(3) rn,rpn,rporn,h0,nmax
         write(3) en,cn,tn
         close(unit=3)
         open(unit=2,file='hfspl.2',form='unformatted',
     &        access='direct',recl=NREC)
 
         nf = 2*n
      endif
 
****  loop over angular momentum values
 
      do 400 lang = 0,lmax
 
      if(lang.eq.0) then
         do 10 i=1,NK
            t(i) = t1(i)
 10      continue
         call setmat(n,k,t)
      elseif(lang.eq.2) then
         do 20 i=1,NK
            t(i) = t2(i)
 20      continue
         call setmat(n,k,t)
      elseif(lang.eq.4) then
         do 30 i=1,NK
            t(i) = t3(i)
 30      continue
         call setmat(n,k,t)
      endif
 
         do  300 ikap = 1,2
           kap = 3*lang + 1 -(2*lang+1)*ikap
           if(kap.eq.0) go to 300
           index = index + 1
 
**********  solve dirac equation *******
 
           call soleqn(k,t,n,kap,ev,cd,iflag)
c	write (6,*) iflag
           if(iflag.ne.0) go to 901
 
**********  calculate corrections
 
           do 100 n0 = n+1,n+5
           write(6,1001) cd(1,n0),cd(n+1,n0),cd(n,n0),cd(nf,n0),n0-n,
     &     kap,ev(n0)
 1001      format(/' g0,f0 =',1p,2d20.10/
     &        ' gn,fn =',1p,2d20.10/
     &        '  w(',i3,i3,') =',1pd20.10)
 100       continue
 
**********  interpolate onto new-grid *********************
          n0 = n+1
          do 200 nn = 1,n
               mm = nn+n
               do 150 i = 1,nmax
                  g(i,nn) = bvalue(t,cd(1,mm),n,k,rn(i),0)
                  f(i,nn) = bvalue(t,cd(n0,mm),n,k,rn(i),0)
 150           continue
 200       continue
 
           write(2,rec=index) g,f,ev
 300     continue
 400  continue
      close(unit=2)
      open (unit=15,file='r')
      do 1670 i=1,max
        write (15,*) i,rn(i)
1670  continue
      close (15)

      return
 901  stop
      end
 
      subroutine hfexc(k,t,n,kang)
      implicit doubleprecision(a-h,o-z)
***********************************************************
*
*    set up hf exchange term for 1s core
*
***********************************************************
      include "global.par"
      PARAMETER (NGF=KX*(NX+1-KX),NK=KX+NX)
      common/excmat/rgg(NX,NX),rgf(NX,NX),rff(NX,NX),vdr(NX,NX)
      common/radial/r(NHF),rp(NHF),rpor(NHF),h,max
      common/shells/wco(NS),nco(NS),kco(NS),mco(NS),jmax
      common/shield/znc(NHF),zhf(NHF),ghf(NHF,NS),fhf(NHF,NS)
      common/tramat/bi(KX*KX,NX),di(KX*KX,NX),xg(KX),wg(KX)
      dimension t(NK),gx(NGF,NS),fx(NGF,NS),zx(NGF),uxx(NGF),vxx(NGF)
      dimension u(NHF),v(NHF),ux(NHF),vx(NHF),biatx(KX)
      data nint /8/
 
         n0 = n+1
         do 40 i = 1,NHF
           u(i) = 0.0
 40      continue
         mss = 0
         do 55 js = 1,jmax
            ms = mco(js)
            deg = 2*iabs(kco(js))
            if(ms.gt.mss) mss = ms
            do 50 i = 1,ms
                u(i) = u(i) + deg*( ghf(i,js)**2 + fhf(i,js)**2 )
 50         continue
 55      continue
         call yfun(u,zhf,0,ms,*901)
 
         do 60 i = 1,max
             zhf(i) = r(i)*zhf(i)
 60      continue
 
****  interpolate zhf, ghf and fhf onto spline grid
 
         do 100 l = k,n
            dl = (t(l+1)-t(l))
            x0 = t(l)
            ind = k*(l-k)
 
****  loop over gaussian points
 
            do 90 m = 1,k
              xm = dl*xg(m) + x0
              do 85 js = 1,jmax
                ms = mco(js)
                rmax = r(ms)
                if(xm.le.rmax) then
                   call interp (ghf(1,js),r,ms,gx(ind+m,js),xm,nint)
                   call interp (fhf(1,js),r,ms,fx(ind+m,js),xm,nint)
                else
                   gx(ind+m,js) = 0.0
                   fx(ind+m,js) = 0.0
                endif
 85           continue
              call interp (zhf,r,max,zx(ind+m),xm,nint)
 90         continue
 100     continue
 
          mm = k*(n+1-k)
 
*      write(6,1900) zx(mm)
 
****   set up direct hf potential  matrix
      do 440 i = 1,n
         do 440 j = 1,n
            vdr(i,j) = 0.0
 440  continue
 
      do 500 i = 1,n
         jhi = min0(i+k-1,n)
         do 490 j = i,jhi
*        sum over intermediate segments
            low = max0(i,k)
            lhi = min0(i+k-1,n)
            do 480 l = low,lhi
               sev = 0.0
               dl = (t(l+1)-t(l))
               x0 = t(l)
*           sum over gaussian weights
               ini = k*(l-i)
               inj = k*(l-j)
               ind = k*(l-k)
               IF(inj.ge.0) THEN
               do 450 m = 1,k
                  xm = dl*xg(m) + x0
                  gm = bi(ini+m,i)*bi(inj+m,j)
                  sev = sev + wg(m)*gm*zx(ind+m)/xm
 450           continue
               END IF
               vdr(i,j) = vdr(i,j) + dl*sev
 480        continue
 490     continue
 500  continue
*   store remainder
      do 600 i = 2,n
         do 550 j = 1,i-1
            vdr(i,j) =  vdr(j,i)
 550     continue
 600  continue
 
****  setup HF exchange
 
*****   set up vl(Bi*g,r) and vl(Bi*f,r)
      do 555 ix=1,NX
        do 555 jx = 1,NX
          rgg(ix,jx) = 0.0
          rgf(ix,jx) = 0.0
          rff(ix,jx) = 0.0
 555  continue
 
      do 315 js = 1,jmax
      ks = iabs(kco(js))
      ls = kco(js)
      if(ls.lt.0) ls = -kco(js)-1
      kin = iabs(kang)
      lin = kang
      if(lin.lt.0) lin = -kang-1
      lmin = iabs(ks-kin)
      lmax = ks + kin -1
 
 
      do 310 lang = lmin,lmax
      lt = ls+lin+lang
      if(mod(lt,2).ne.0) go to 310
      angfac = cre(ks,lang,kin)**2/(2*kin)
      do 300 i = 1,n
          rimin = t(i)
          rimax = t(i+k)
 
          mig = 1
          do 210 ig = 1,max
             if(r(ig).ge.rimin.and.r(ig).le.rimax)  then
                if(r(ig).gt.r(mig)) mig = ig
                 call interv (t,n+k,r(ig),left,mflag)
                 call bsplvb(t,k,1,r(ig),left,biatx)
                 u(ig) = ghf(ig,js)*biatx(i-left+k)
                 v(ig) = fhf(ig,js)*biatx(i-left+k)
              else
                 u(ig) = 0.0
                 v(ig) = 0.0
              endif
 210      continue

           if (mig.lt.10) then
           write (*,109) mig
109        format (' mig in hfexc is equal to ',i5,
     *   ' It was changed to 10')
           mig=10
           endif


          call yfun(u,ux,lang,mig,*901)
          call yfun(v,vx,lang,mig,*901)
 
****  interpolate ux, vx onto spline grid
 
         if(lang.eq.0) go to 2091
          do 209 ig = 1,max
              ux(ig) = ux(ig)*r(ig)**lang
              vx(ig) = vx(ig)*r(ig)**lang
 209      continue
 2091     continue
 
          do 215 l = k,n
              dl = (t(l+1)-t(l))
              x0 = t(l)
              ind = k*(l-k)
             do 212 m = 1,k
                 xm = dl*xg(m) + x0
                    call interp (ux,r,max,uxx(ind+m),xm,nint)
                    call interp (vx,r,max,vxx(ind+m),xm,nint)
 212          continue
 215      continue
 
         do 250 j = 1,n
*        sum over intermediate segments
            low = max0(j,k)
            lhi = min0(j+k-1,n)
            do 230 l = low,lhi
               sgg = 0.0
               sgf = 0.0
               sff = 0.0
               dl = (t(l+1)-t(l))
               x0 = t(l)
*           sum over gaussian weights
               inj = k*(l-j)
               ind = k*(l-k)
               IF(inj.ge.0) THEN
               do 220 m = 1,k
                  xm = dl*xg(m) + x0
                  gm = wg(m)*bi(inj+m,j)/xm**lang
                  sgg = sgg + gm*uxx(ind+m)*gx(ind+m,js)
                  sgf = sgf + gm*uxx(ind+m)*fx(ind+m,js)
                  sff = sff + gm*vxx(ind+m)*fx(ind+m,js)
 220           continue
               END IF
               rgg(i,j) = rgg(i,j) + dl*sgg*angfac
               rgf(i,j) = rgf(i,j) + dl*sgf*angfac
               rff(i,j) = rff(i,j) + dl*sff*angfac
 
 230        continue
 250     continue
 300  continue
 310  continue
 315  continue
      return
 
 901  stop ' yfun error in hfexc'
      end
 
      subroutine soleqn(k,t,n,kap,ev,ad,iflag)
      implicit doubleprecision(a-h,o-z)
      external dsygv
      character*1 jobz,uplo
***********************************************************
*
*    set up matrices for 1-dim dirac equation
*
***********************************************************
      include "global.par"
      parameter (NL=2*NX,LWORK=4*NL,NK=KX+NX)
      dimension t(NK)
      dimension ad(NL,NL),bd(NL,NL)
      dimension ev(NL),work(LWORK)
      common/excmat/rgg(NX,NX),rgf(NX,NX),rff(NX,NX),vdr(NX,NX)
      common/matrix/a(NX,NX),b(NX,NX),c(NX,NX),v(NX,NX)
      data cl/137.0359895d0/, itype/1/, jobz/'V'/, uplo/'U'/
 
 
      call hfexc(k,t,n,kap)
 
      cl2 = 2*cl**2
      nf = 2*n
      np1 = n+1
      do 50 i = 1,n
         do 30 j = 1,n
            ad(i,j) = v(i,j) + vdr(i,j) - rgg(i,j)
            bd(i,j) = b(i,j)
 30      continue
         do 40 j = np1,nf
            jj = j-n
            ad(i,j) = cl*(a(i,jj)-kap*c(i,jj)) - rgf(i,jj)
            bd(i,j) = 0.0
 40      continue
 50   continue
 
      do 80 i = np1,nf
         ii = i-n
         do 60 j = 1,n
            ad(i,j) = ad(j,i)
            bd(i,j) = bd(j,i)
 60      continue
         do 70 j = np1,nf
            jj = j-n
            ad(i,j) = -cl2*b(ii,jj) + v(ii,jj) + vdr(ii,jj)
     &                - rff(ii,jj)
            bd(i,j) = b(ii,jj)
 70      continue
 80   continue
 
      if(kap.lt.0) then
         ad(1,1) = ad(1,1) + cl
      else
         ad(1,1) = ad(1,1) + cl2
      endif
      ad(1,n+1) = ad(1,n+1) - 0.5*cl
      ad(n+1,1) = ad(n+1,1) - 0.5*cl
 
      ad(n,n) = ad(n,n) + cl/2.0
      ad(nf,nf) = ad(nf,nf) - cl/2.0
 
***********************************************************
*     call eispack to solve general eigenvalue problem
*     A c = eps B c
***********************************************************
*     call rsg(NL,nf,ad,bd,ev,matz,cd,fv1,fv2,ierr)
      call dsygv(itype,jobz,uplo,nf,ad,NL,bd,NL,ev,
     1           work,LWORK,info) 
c          write(6,*) ' info =',info
      iflag = info
      if(info.ne.0) then
         write(6,*) ' LAPACK WARNING info =',info
	   stop
         go to 999
      endif
      return
 
 999  iflag = info
      end
 
      subroutine setmat(n,k,t)
      implicit doubleprecision(a-h,o-z)
**********************************************
*
*      b (i,j) = Int[ Bi(x) Bj(x) dx ]
*
*      b(j,i) = b(i,j)
*
*      for i<=j  :
*                       t(i+k)
*               b(i,j) =  Int [Bi(x) Bj(x) dx]  j <i+k
*                         t(j)
*               b(i,j) =           0            j>=i+k
*
*      therefore
*                        l=i+k-1  t(l+1)
*               b(i,j)  =  Sum   { Int [ Bi(x) Bj(x) dx ] }
*                          l=j     t(l)
*
*       where t(l), l = 1,n+k = knot sequence
*
*     t(l+1)                             k
*      Int [ f(x) dx ]  = [t(l+1)-t(l)] Sum [ w(m) f(x(m)) ]
*     t(l)                              m=1
*
*      x(m) = [t(l+1)-t(l)] z(m) + t(l)
*      z(m) = gaussian k-point coordinates for [0,1]
*      w(m) = gaussian k-point weights     for [0,1]
*
******************************************************************
 
      include "global.par"
      parameter(NK=KX+NX)
      common/matrix/a(NX,NX),b(NX,NX),c(NX,NX),v(NX,NX)
      common/tramat/bi(KX*KX,NX),di(KX*KX,NX),xg(KX),wg(KX)
      dimension t(NK)
      dimension db(KX,KX)
      data nderiv/2/
 
*  initialize  a and b
      do 20 j = 1,n
        do 10  i = 1,n
          a(i,j) = 0.0
          b(i,j) = 0.0
          c(i,j) = 0.0
          v(i,j) = 0.0
 10     continue
 20   continue
 
      call gauss(k,xg,wg)
 
*****  store matrix B(jg,i) i = 1,...,n  jg = 1,....,k*k
*
*    values of B(x,i) for gaussian points on support interval
*  x  = t(i) + xg(1..k) , t(i+1) + xg(1..k),.. t(i+k-1) + xg(1..k)
*       also   B(x,i)' on same interval
*
*************
 
      do 50 l = k,n
          dl = (t(l+1)-t(l))
          x0 = t(l)
*     sum over gaussian weights
          imin = l-k+1
          imax = l
          do 40 m = 1,k
              xm = dl*xg(m) + x0
              call bsplvd(t,k,xm,l,db,nderiv)
 
              do 30 i = imin,imax
                    ind = k*(l-i)+m
                    bi(ind,i) = db(i-l+k,1)
                    di(ind,i) = db(i-l+k,2)
 30           continue
 40       continue
 50   continue
 
*  set up non zero loops  i = 1 ... n  ,j = i,i+k-1
      do 500 i = 1,n
         jhi = min0(i+k-1,n)
         do 400 j = i,jhi
*        sum over intermediate segments
            low = max0(i,k)
            lhi = min0(i+k-1,n)
            do 300 l = low,lhi
               sea = 0.0
               seb = 0.0
               sec = 0.0
               sev = 0.0
               dl = (t(l+1)-t(l))
               x0 = t(l)
*           sum over gaussian weights
               ini = k*(l-i)
               inj = k*(l-j)
               IF(inj.ge.0) THEN
               do 200 m = 1,k
                  xm = dl*xg(m) + x0
                  gm = bi(ini+m,i)*bi(inj+m,j)
                  hm = bi(ini+m,i)*di(inj+m,j)
                  sea = sea + wg(m)*hm
                  seb = seb + wg(m)*gm
                  sec = sec + wg(m)*gm/xm
                  sev = sev + wg(m)*gm*vf(xm)
 200           continue
               END IF
               if(i.ne.j) a(i,j) = a(i,j) + dl*sea
               b(i,j) = b(i,j) + dl*seb
               c(i,j) = c(i,j) + dl*sec
               v(i,j) = v(i,j) + dl*sev
 300        continue
 400     continue
 500  continue
 
*   store remainder
      do 600 i = 2,n
         do 550 j = 1,i-1
            a(i,j) = -a(j,i)
            b(i,j) =  b(j,i)
            c(i,j) =  c(j,i)
            v(i,j) =  v(j,i)
 550     continue
 600  continue
 
      return
      end
 
      subroutine mkgrid(n,k,rmax,dr0,t,*)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
********************************************************
*
*  input :
*         n = number of collacation points
*         k = order of spline > 2
*        a,b = 1st and last points on interval
*  output :
*         t(1) = t(2) = ... t(k) = a
*         t(n+1) = t(n+2) = ... t(n+k)
*         t(i) = a + (i-k)*h   i = k+1,n
*           with  h = (b-a)/(n-k+1)
*
***********************************************************
      parameter (NK=KX+NX)
      dimension t(NK)
 
      do 50 i = 1,k
         t(i) = 0.0
 50   continue
*   ipow = no of internal points + 2 endpoints
      ipow = n-k+1
      call step(ipow,rmax,dr0,t0,h,*901)
 
      write(6,1000) ' t0 =',t0,'  h =',h
 1000 format(/a,1pd14.6,a,1pd14.6/)
      do 100 i = k+1,n
         t(i) =  t0*(exp(h*(i-k))-1.0)
 100  continue
 
      do 150 i = n+1,n+k
         t(i) = rmax
 150  continue
 
      return
 901  return 1
      end
 
      subroutine step(ipow,rmax,dr,t0,h,*)
      implicit doubleprecision(a-h,o-z)
      data small /1.0d-10/,nmax/54/
 
         f(x) = 1.0 + con*x - (1.0+x)**ipow
 
       con = rmax/dr
       fa = f(small)
       a = small
       if(fa.le.0.0) then
          write(6,*) 'con=',con,'  ipow=',ipow
          write(6,*) ' a =',a,  '    fa=',fa
          write(6,*) ' problem with fa !!!'
          return 1
       endif
       b = 1.0
       fb = f(b)
 
       if(fb.ge.0.0) then
            write(6,*) ' b =',b,'   fb=',fb
            write(6,*) ' problem with fb !!!'
            return 1
       endif
       dx = 1.0
       i = 0
 200   i = i+1
       dx = 0.5*dx
       c = 0.5*(a+b)
       fc = f(c)
       if(fb*fc.gt.0.0) then
          b = c
          fb = fc
       else
          a = c
          fa = fc
       endif
c      write(6,*) ' i=',i,' c=',c,'  f=',fc
       if(fc.ne.0.0.and.i.lt.nmax) go to 200
 
       h =  log(1.0+c)
       t0 = dr/c
 
       return
       end
 
      function vf(x)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      common/shield/znc(NHF),zhf(NHF),ghf(NHF,NS),fhf(NHF,NS)
      common/radial/r(NHF),rp(NHF),rpor(NHF),h,max
 
      call interp (znc,r,max,y,x,6)
      vf = -y/x
 
      return
      end
 
      subroutine infun
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      character*4 ident
      common/radial/r(NHF),rp(NHF),rpor(NHF),h,max
      common/shells/wco(NS),nco(NS),kco(NS),mco(NS),jmax
      common/shield/znc(NHF),zhf(NHF),ghf(NHF,NS),fhf(NHF,NS)
      common /rho0/ en,cn,tn
      read(1) ident,jmax,jz,nat,inuc,ion,io,in,inf
      write(6,1000) ident
 1000 format(/' Core input for ',a)
      mss = 0
      do 10 js = 1,jmax
         read(1) n,kap,k,l,jdg,ms,iof,wf
         write(6,1010) n,kap,wf/2
         wco(js) = wf/2
         nco(js) = n
         kco(js) = kap
         mco(js) = ms
         if(ms.gt.mss) mss = ms
 1010    format(/'   n=',i3,5x,' kap=',i3,5x,' w =',1pd20.12)
 10   continue
      read(1) r,rp,rpor,h,max
      do 20 js = 1,jmax
          read(1) (ghf(i,js),i=1,max),(fhf(i,js),i=1,max)
 20   continue
 
      read(1) znc,zhf,cn,tn,en
 
      if(in.ne.0) jmax=in-1
      write(6,1020) r(mco(jmax))
 1020 format(/'   R(core) =',f10.4,'au'/)
 
      return
      end
 
      subroutine interp (ya,xa,n,y,x,np)
      implicit doubleprecision(a-h,o-z)
************************************************************************
*  this program uses aitken's method to do interpolations
 
*  ya    : input array
*  xa    : original grid where ya is defined
*  n     : length of arrays xa and ya
*  y     : output interpolated value
*  x     : point at which y is defined
*  np    : order of the interpolation method (NMAX >= np >= 2)
 
*  n.b.    xa must be monotonically increasing
 
************************************************************************
 
      parameter(NMAX=20)
      dimension ya(n),xa(n)
      dimension ty(NMAX),tx(NMAX),sy(NMAX)
      data ilo /0/
 
*    find ilo such that xa(ilo) =< x < xa(ilo+1)
*    assume last value of ilo is satisfactory as initial guess
 
      if(ilo.eq.0.or.ilo.gt.n) then
         ilo = 0
         ihi = n + 1
         go to 300
      endif
 
      inc = 1
      if(x.gt.xa(ilo)) then
 100     ihi = ilo+inc
         if(ihi.gt.n) then
             ihi = n + 1
         else if(x.ge.xa(ihi)) then
             ilo = ihi
             inc = inc + inc
             go to 100
         endif
      else
         ihi = ilo
 200     ilo = ihi - inc
         if(ilo.lt.1) then
            ilo = 0
         else if(x.lt.xa(ilo)) then
            ihi = ilo
            inc = inc + inc
            go to 200
         endif
      endif
 
 300  if(ihi-ilo.eq.1) go to 400
      im = (ihi + ilo)/2
      if(x.gt.xa(im)) then
         ilo = im
      else
         ihi = im
      endif
      go to 300
 400  continue
 
*  aitken's np point method of interpolation
 
      if(np.gt.NMAX) np = NMAX
 
      nh = np/2
      ia = max0(ilo-nh,0)
      ia = min0(ia,n-np)
 
      do 450 k = 1,np
         ty(k)=ya(ia+k)
         tx(k)=xa(ia+k)
 450  continue
      do 600 j=2,np
         do 500 k=j,np
             dx=(x-tx(j-1))/(tx(k)-tx(j-1))
             sy(k)=(1.0-dx)*ty(j-1)+dx*ty(k)
 500     continue
         do 550 k=j,np
             ty(k)=sy(k)
 550     continue
 600  continue
      y=ty(np)
      return
      end
      
      subroutine newgrid(h0,hp,rmax,max,r0,r,rp,rpor)
      implicit double precision(a-h,o-z)
c
c   Set up new radial grid:  log(r/r0+1) + a*r = (n-1)*h
c   here  a = h/hp, where hp is the effective step size of the radial
c   grid in the asymptotic region.  When a=0, r = r0*(exp(n-1)*h)-1)
c   Ripped off from the RRPA code, May 21, 1997
c
c   INPUT: h0 
c        : hp
c        : rmax
c        : max
c  OUTPUT: r0
c        : r
c        : rp
c        : rpor
c
      include "global.par"
      dimension r(NGP),rp(NGP),rpor(NGP)
      data itmax/100/

      r(1) = 0.d0
      rpor(1) = 0.d0
      a = h0/hp
      arg1 = (max-1) * hp
      arg2 = rmax

      IF(arg1 .le. arg2) THEN
        write(6,1000) arg1,arg2
 1000   format(' (max-1) * hp',1p,d12.4,'  must be > Rmax',d14.4)
        stop
      END IF
      r0 = rmax/(exp((max-1)*h0-a*rmax)-1.d0)
      rp(1) = r0 / (a*r0+1.d0)

      
      eps = 100 * EPSILON(r0)

       DO k = 2,max
         xi = r(k-1)
         xpi = rp(k-1)
         t = (k-1)*h0
         DO it = 1,itmax
            del = (t-a*xi-log(xi/r0+1.d0))*xpi
            xj = xi + del
            xpj = (xj+r0)/(a*(xj+r0)+1.d0)
            IF(abs(del/xj).lt.eps) go to 24
            xi=xj
            xpi=xpj
         END DO
         write(6,*) ' Too many iterations in newgrid'
         stop
         
 24      continue
         r(k) = xj
         rp(k) = xpj
         rpor(k) = xpj/xj
      END DO

      return
      end

c     ========================================
      include "bsplvb.f"
      include "bsplvd.f"
      include "bvalue.f"
      include "clrx.f"
      include "cre.f"
      include "gauss.f"
      include "interv.f"
      include "rint.f"
      include "rsg.f"
      include "yfun1.f"
      include "yint1.f"