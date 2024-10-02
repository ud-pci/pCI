*******BREIT FROM NEW TDHF , BREIT AND QED NOT REMOVED

      program tdhf
      implicit doubleprecision(a-h,o-z)
***********************************************************************
*
*        RS/6000 version 
*     
*    addition of excited states     14/5/83
*    500 point grid version         3/23/85
*    formats modified for PC       12/12/85
*    modified for CRAY              7/21/86
*    Breit                          3/20/87
*    QED + distorted nucleus        4/15/89
*    clean up messer mess          22/10/91
*    replace dminv by lapack       11/10/94
***********************************************************************
      call datain(xa,ex,r0,hh,mm,*901)
      call setgri(r0,hh,mm,*901)
      call hart(ex)
      call pick
      open(unit=1,form='unformatted',file='fort.1')
      rewind 1
      call output(xa)
      rewind 1
      close(unit=1)
 901  stop
      end

      block data
      implicit doubleprecision(a-h,o-z)
      common/phycon/alpha,pi,ev,bohr,ryd
      data alpha /   137.035 989 5 d0/ 
      data pi    / 3.141 592 653 589 793 d0/
      data bohr  / 0.529 177 249 d0/
      data ev    /13.605 698 1 d0/
      data ryd   /   219 474.631 42 d0/
      end

      subroutine datain(xa,ex,r0,hh,mm,*)
      implicit doubleprecision(a-h,o-z)
**********************************************************************
*
*    Read in data and calculate DHF exchange factors
*
*********************************************************************
      include "global.par"
      character*4 ident,lab,label
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/lam/clam(5,10,5)
     &      /phycon/alpha,pi,ev,bohr,ryd
      common/orblab/label
      common/nucdat/rnuc,cnuc,tnuc,anuc,b2,b4,en,iparm
**********************************************************************

      read(5,1000) ident,jmax,jz,nat,nuc,ion,io,in,inf
 1000 format(a,8i4)
      write(6,1010)
 1010 format('   *****  input data stream image  *****'/)
      write(6,1000) ident,jmax,jz,nat,nuc,ion,io,in,inf
      if(jmax.gt.NS) then
        write(6,1011) ' Too many input orbitals !!  shell dim = ',NS
 1011   format(a,i4)
        go to 901
      endif
      j1=jz
      kmax=0
      do 100 j=1,jmax
           read(5,1020) n(j),kap(j),iof(j),wh(j)
 1020      format(3i4,f12.4)
           k(j)=iabs(kap(j))
           if(k(j).gt.kmax) kmax=k(j)
           jdg(j)=2*k(j)
           if(inf.eq.1.and.j.ge.in) jdg(j)=1
           if(inf.eq.2.and.j.eq.jmax) jdg(j)=1
           l(j)=kap(j)
           if(l(j).lt.0) l(j)=-l(j)-1
           call enclab(n(j),kap(j),*901)
           lab(j)=label
           an=n(j)
           ak=k(j)
           az=j1/alpha
           gam=dsqrt(ak**2-az**2)
           en2=an**2-2*(an-ak)*(ak-gam)
           en1=dsqrt(en2)
*****   wh(j) is input in a.u. but stored internally in Ry.
           wh(j) = 2*wh(j)              
           if(wh(j).ge.0.0) then
              wh(j)=-2*j1**2/(en2+en1*(an-ak+gam))
              if(inf.eq.1.and.j.ge.in)
     &           wh(j)=-(ion+1)**2/(n(j)-n(in-1)-.325)**2
           endif                                           
           write(6,1030) n(j),kap(j),iof(j),lab(j),wh(j)
 1030      format(3i4,4x,a,4x,f16.4)
 100  continue

*************  inf = 1  ==> frozen core + valence

      if(inf.eq.1) then
           read(5,*) xa
           write(6,1060) xa
 1060      format(' xalpha =',f12.4)
           ion=ion+1
           itemp=in
           in=jmax
           jmax=itemp-1
      endif  

      read(5,*) r0,hh,mm

      read(5,*) iparm
      if(iparm.eq.1) then
         read(5,*) rnuc,cnuc,tnuc
      else
         read(5,*) cnuc,anuc,b2,b4
      endif
 
      read(5,*) ex

********* calculate DHF exchange angular momentum coupling factors

      max=2*kmax
      do 400  j = 1,max
           do 300 ka = 1,kmax
                do 200  kb = ka,kmax
                     as=cre(ka,j-1,kb)
                     as=as**2/(4*ka*kb)
                     clam(ka,j,kb)=as
                     clam(kb,j,ka)=as
 200            continue
 300       continue
 400  continue         

      return
 901  return 1
      end

      subroutine setgri(r0,hh,mm,*)
      implicit doubleprecision(a-h,o-z)
***********************************************************************
*
*  routine to set up grid and nuclear potential
*  revised 28/4/83/
*  revised for PC 13/12/85
*
***********************************************************************
      include "global.par"
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/charge/znuc(NGP),z(NGP),zn(NGP),rho(NGP),y(NGP),q(NGP)
     &      /nucdat/rnuc,cnuc,tnuc,anuc,b2,b4,en,iparm
     &      /phycon/alpha,pi,ev,bohr,ryd      
      data nmax /NGP/,hdef/0.03125/,rdef/5e-4/

      if(r0.le.0.0) r0=rdef
      h = hh
      if(h.le.0.0) h=hdef
      max = mm
      if(max.le.0.or.max.gt.nmax) max=nmax
      write(6,2400) r0,h,max
 2400 format(/' grid parameters:'
     &       /' r0=',f13.5,'  h=',f9.5,'   max=',i6)

**********    set up grid

      r(1)=0.0
      rp(1)=r0
      rpor(1)=0.0
      do 200 i=2,max
           rp(i)=dexp((i-1)*h)
           r(i)=r0*(rp(i)-1.0)
           rp(i)=r0*rp(i)
           rpor(i)=rp(i)/r(i)
 200  continue
************  calculate nuclear potential
      if(iparm.eq.1) then
         call potent(jz,nat,rnuc,cnuc,tnuc,en,z,inuc,*901)
         write(6,1900) rnuc,cnuc,tnuc,inuc
 1900    format(/' r(rms) =',f7.5,'  c =',f7.5,'   t =',f4.2,
     &          /' i(nuc) =',i5/)
      else
         call potnew(jz,cnuc,anuc,b2,b4,rnuc,en,z,inuc)    
         write(6,1950) rnuc,cnuc,anuc,b2,b4,inuc
 1950    format(/' r(rms) =',f7.5,'  c =',f7.5,'   a =',f7.5/
     &           '     b2 =',f7.5,' b4 =',f7.5/
     &           ' i(nuc) =',i5/)
      endif

      do 100 i = 1,max
         znuc(i)=z(i)
         zn(i)=z(i)
 100  continue
      return
 901  return 1
      end

      subroutine output(xalpha)
      implicit doubleprecision(a-h,o-z)
***********************************************************************
*
*    1)  Calculate total energy
*    2)  decide on valence action
*    3)  prepare output list
*    4)  determine scalar products
*    5)  determine Breit corrections
*    6)  write wave functions on tape
*
***********************************************************************
      include "global.par"
      PARAMETER(NSCAL=10)
      character*4 ident,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/shell/rk(NGP),rl(NGP),dumm(16)
      common/charge/znuc(NGP),z(NGP),yhf(NGP),y(NGP),u(NGP),ysum(NGP)
      common/nucdat/rnuc,cnuc,tnuc,anuc,b2,b4,en,iparm
      common/lam/clam(5,10,5)
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension ebr(NS)
      dimension over(NSCAL),icnt(NSCAL),js(NSCAL,NSCAL)
      dimension vpint(NS),seint(NS),dfs(NS)
      dimension vk(NGP),vl(NGP)
*********  calculate total energy

      etot=0.0
      do 70 jk = 1,jmax
         mk=ms(jk)
         kk=k(jk)

         do 20 i = 1,max
             rk(i)=0.0
             rl(i)=0.0
 20      continue 

         do 60 jl = 1,jmax
            if(jk.eq.jmax.and.jl.eq.jmax.and.inf.eq.2) go to 60
            ml=ms(jl)
            kl=k(jl)
            do 25 i = 1,ml
               u(i) = g(i,jl)**2 + f(i,jl)**2
 25         continue
            call yfun(u,y,0,ml,*901)
            adg = jdg(jl)
            do 30 i = 1,mk
               rk(i) = rk(i) - adg*y(i)*g(i,jk)
               rl(i) = rl(i) - adg*y(i)*f(i,jk)
 30         continue
            mlow=min0(mk,ml)
            do 35 i=1,mlow
               u(i) = g(i,jl)*g(i,jk) + f(i,jl)*f(i,jk)
 35         continue
            low = iabs(kk-kl)
            lhi = kk + kl - 1
            do 50 ll=low,lhi
               ltest = l(jk) + l(jl) + ll
               if(mod(ltest,2).ne.0) go to 50
               call yfun(u,y,ll,mlow,*901)
               cl = clam(kk,ll+1,kl)
               adg = cl*jdg(jl)
               do 40 i = 1,ml
                  rk(i) = rk(i) + adg*y(i)*g(i,jl)
                  rl(i) = rl(i) + adg*y(i)*f(i,jl)
 40            continue 
 50         continue
 60      continue
              call brpot(jk,vk,vl)
              do 85 i = 1,mk
                 rk(i) = rk(i) - vk(i)
                 rl(i) = rl(i) - vl(i)
85            continue


         do 65 i=1,mk
            u(i) = (g(i,jk)*rk(i) + f(i,jk)*rl(i))*rp(i)
 65      continue
         term = rint(u,1,mk,7,h)
         etot = etot + jdg(jk)*(wf(jk) + term)
 70   continue

***********  calculate valence orbitals if inf = 1

      if(inf.eq.1) then
            call valenc(xalpha)
            ion=ion-1
            itemp=jmax
            jmax=in
            in=itemp+1
      endif

************   calculate and print scalar products

      write(6,1000)
 1000 format(/' *****  scalar products   *****')

      ic = 0
      do 120 kin = 1,4
         do 110 jin = 1,2
             ic = ic + 1
             kaps = -kin
             if(jin.eq.2)  kaps = kin
             ins = 0
             do 100 j = 1,jmax
                 if(kap(j).ne.kaps) go to 100
                 ins = ins + 1
                 js(ins,ic) = j
 100         continue
             icnt(ic) = ins

 110     continue
 120  continue

      ic = 0
      do 160 kin = 1,4
         do 150 jin = 1,2
            ic = ic + 1
            kaps = -kin
            if(jin.eq.2)  kaps = kin
            jns = 0
            do 140 j = 1,jmax
               if(kap(j).ne.kaps) go to 140
               jns = jns + 1
               ins = 0
               do 130 jl = 1,j
                  if(kap(jl).ne.kaps) go to 130
                  ins = ins + 1
                  mm=min0(ms(j),ms(jl))
                  do 125 i=1,mm
                      u(i)=(g(i,j)*g(i,jl)+f(i,j)*f(i,jl))*rp(i)
 125              continue
                  over(ins) = rint(u,1,mm,7,h)
                  if(jl.eq.j) over(ins) = over(ins) - 1.0
                  if(jns.eq.1) then
                     write(6,1010) (lab(js(ii,ic)),ii = 1,icnt(ic))
 1010                format(/10x,8('    ',a,'> '))
                  endif
 130           continue
               write(6,1020) lab(j),(over(ii),ii=1,ins)
 1020          format(' <',a,'|   ',1p,8e10.1)
 140        continue
 150     continue
 160  continue

**********  evaluate breit interaction
*    ebr(ja) = correction to ja'th eigenvalue
*    brc     = breit correction to core (for inf = 0 or 1)
*    brc     = breit correction to atom (for inf = 2)
**********

      call breit(ebr,brc) 
      call qed(vpint,seint)
      qtot = 0.0
      if(inf.eq.1) then
        jup = in-1
      else
        jup = jmax
      endif
      do 200 j = 1,jup
         qtot = qtot + jdg(j)*(seint(j)+vpint(j))
 200  continue

******   estimate of finite nuclear size correction
      do 210 i = 2,max
         y(i) = (dfloat(jz)-znuc(i))/r(i)
 210  continue
         y(1) = 0.0
      do 230 j = 1,jmax
          mj = ms(j)
          do 220 i = 1,mj
             u(i) = (g(i,j)**2+f(i,j)**2)*y(i)*rp(i)
 220     continue
         dfs(j) = rint(u,1,mj,7,h)
 230  continue

***********  write output summary

      write(6,1030) ident
 1030 format(/6x,' Dirac Hartree Fock Levels for ',a/)
      write(6,1040)
 1040 format(' shell',4x,' m',3x,'atomic units','       breit',
     &                        6x,'      qed',6x,'      f.s.'/)

      do 300 j = 1,jmax
          wpt=wf(j)
          wau=.5*wpt
          qqed = seint(j)+vpint(j)
          write(6,1050) lab(j),ms(j),wau,ebr(j),qqed,dfs(j)
 1050     format(a,i8,f17.9,3f15.7)
 300  continue
      efin=etot
      wau=.5*efin
      write(6,1060) wau,brc,qtot
 1060 format(//4x,' etotal:',3f15.7)
      if(inf.eq.1) then
          etot = etot+wf(in)
          wau = 0.5*etot
          brc = brc + ebr(in)
          qtot = qtot+seint(in)+vpint(in)
          write(6,1070) wau,brc,qtot
 1070     format(/'  with val :',3f15.7)
      endif

************  write answers on file 1
      open(unit=1,file='fort.1',form='unformatted')
      rewind 1
      write(1) ident,jmax,jz,nat,nuc,ion,io,in,inf
      do 325 j=1,jmax
          write(1) n(j),kap(j),k(j),l(j),jdg(j),ms(j),iof(j),wf(j)
 325  continue
      write(1) r,rp,rpor,h,max
      do 335 j=1,jmax
         write(1) (g(i,j),i=1,max),(f(i,j),i=1,max)
 335   continue
      do 340 i = 1,max
         ysum(i) = ysum(i)*r(i)
 340   continue

      write(1) znuc,ysum,cnuc,tnuc,en

      return
 901  stop
      end

      subroutine rside(jk)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/shell/rk(NGP),rl(NGP),dumm(16)
      common/soln/x(NGP),w(NGP)
      common/charge/znuc(NGP),z(NGP),yhf(NGP),y(NGP),u(NGP),ysum(NGP)
      common/lam/clam(5,10,5)
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension tk(NGP),tl(NGP)
      dimension vk(NGP),vl(NGP)
      data ii /32/

      if(ii.ge.32) then
         do 10 i = 1,max
             u(i) = 0.0
 10      continue
         mx = 0
         do 30 j = 1,jmax
            mj = ms(j)
            if(mj.gt.mx) mx = mj
            do 20 i = 1,mj
               u(i) = u(i) + jdg(j)*(g(i,j)**2 + f(i,j)**2)
 20         continue
 30      continue
         call yfun(u,ysum,0,mx,*901)
         ii = 0
      endif
      ii = ii + 1
      mk = ms(jk)
      kk = k(jk)
      do 40 i = 1,mk
         u(i) = g(i,jk)**2 + f(i,jk)**2
 40   continue    
      call yfun(u,x,0,mk,*901)
      lhi = 2*kk - 1
      if(lhi.ge.2) then
         do 60 ll = 2,lhi,2
            call yfun(u,y,ll,mk,*901)
            cl = clam(kk,ll+1,kk)
            adg = (jdg(jk)-1)*cl*(2.0d0*kk)/(2.0d0*kk-1.0d0)
            do 50 i = 1,mk
               x(i) = x(i) + adg*y(i)
 50         continue 
 60      continue
      endif
      do 70 i = 1,mk
         zz = yhf(i) + (x(i)-ysum(i))*r(i)
         rk(i) = zz*g(i,jk)
         rl(i) = zz*f(i,jk)
 70   continue
      mk1 = mk + 1
      if(mk1.le.max) then
         do 75 i = mk1,max
            rk(i) = 0.0
            rl(i) = 0.0
 75      continue
      endif
      do 110 jl = 1,jmax
         if(jl.eq.jk) go to 110
         kl = k(jl)
         ml = ms(jl)
         mlow = min0(mk,ml)
         do 80 i = 1,mlow
            u(i) = g(i,jl)*g(i,jk) + f(i,jl)*f(i,jk)
 80      continue
         low = iabs(kk-kl)
         lhi = kk + kl - 1
         do 100 ll=low,lhi
            ltest = l(jk) + l(jl) + ll
            if(2*(ltest/2).ne.ltest) go to 100
            call yfun(u,y,ll,mlow,*901)
            cl = clam(kk,ll+1,kl)
            adg = cl*jdg(jl)
            do 90 i = 1,ml
               rk(i)=rk(i) + adg*y(i)*r(i)*g(i,jl)
               rl(i)=rl(i) + adg*y(i)*r(i)*f(i,jl)
 90         continue
 100     continue
 110  continue
      call brpot(jk,vk,vl)
c       write (*,*) jk,mk
      do 271 i = 1,mk
         rk(i) = rk(i) - vk(i)*r(i)
         rl(i) = rl(i) - vl(i)*r(i)
 271  continue

******* lagrange multipliers  ********
      if(inf.ne.2) return
      
      if(jk.lt.jmax) then
      if(kap(jk).ne.kap(jmax)) return

***  lagrange multiplier for core state with proper kap

         jv = jmax
         mv = ms(jv)
         kv = k(jv)
         do 120 i = 1,mv
            u(i) = g(i,jv)**2 + f(i,jv)**2
 120     continue
         call yfun(u,y,0,mv,*901)
         adg = -1.0d0/(2.0d0*kv)
         do 130 i = 1,mv
            tk(i) = adg*y(i)*g(i,jv)
            tl(i) = adg*y(i)*f(i,jv)
 130     continue
         lhi = 2*kv - 1
         if(lhi.ge.2) then
            do 150 ll = 2,lhi,2
                call yfun(u,y,ll,mv,*901)
                adg = clam(kv,ll+1,kv)/(2.0d0*kv-1.0d0)
                do 140 i = 1,mv
                   tk(i) = tk(i) + adg*y(i)*g(i,jv)
                   tl(i) = tl(i) + adg*y(i)*f(i,jv)
 140            continue
 150       continue
         endif
         mm = min0(mk,mv)
         do 160 i = 1,mm
            u(i) = (tk(i)*g(i,jk) + tl(i)*f(i,jk))*rp(i)
 160     continue

         alam = rint(u,1,mm,7,h)
         write(6,1000) jk,jv,alam
 1000    format(10x,'lm(',i2,',',i2,')=',1pd14.7)
         do 200 i = 1,mv
            rk(i) = rk(i) + alam*r(i)*g(i,jv)
            rl(i) = rl(i) + alam*r(i)*f(i,jv)
 200     continue
         return
      else

*** lagrange multipliers for valence state
         jv = jk
         mv = mk
         kv = kk
         do 210 i = 1,mv
            u(i) = g(i,jv)**2 + f(i,jv)**2
 210     continue

         call yfun(u,y,0,mk,*901)
         adg = -1.0d0/(2.0d0*kv)
         do 220 i = 1,mv
            tk(i) = adg*y(i)*g(i,jv)
            tl(i) = adg*y(i)*f(i,jv)
 220     continue
         lhi = 2*kv - 1
         if(lhi.ge.2) then
            do 250 ll = 2,lhi,2
                call yfun(u,y,ll,mv,*901)
                adg = clam(kv,ll+1,kv)/(2.0d0*kv-1.0d0)
                do 240 i = 1,mv
                   tk(i) = tk(i) + adg*y(i)*g(i,jv)
                   tl(i) = tl(i) + adg*y(i)*f(i,jv)
 240            continue
********  add in breit  ******
  250   continue

         endif

         do 300 jl = 1,jmax-1
            if(kap(jl).ne.kap(jv)) go to 300
            ml = ms(jl)
            mm = min0(mk,ml)
            do 260 i = 1,mm
               u(i) = (g(i,jl)*tk(i) + f(i,jl)*tl(i))*rp(i)
 260        continue
            alam = rint(u,1,mm,7,h)
            write(6,1000) jk,jl,alam
            alam = jdg(jl)*alam
            do 270 i = 1,ml
               rk(i) = rk(i) + alam*g(i,jl)*r(i)
               rl(i) = rl(i) + alam*f(i,jl)*r(i)
 270        continue
 300     continue
      endif 
      return
 901  stop
      end

      subroutine brpot(jj,vl,vs)
***********************************************************************
*
*     purpose : Evaluate Single Particle Potential of Breit Interaction
*               V(r)*F(jj,r)
*
*     input: jj  =  index of HF orbital upon which V acts
*
*     output:  vl,vs large and small components of V*F
*
************************************************************************
      include "global.par"
      implicit doubleprecision(a-h,o-z)
      character*4 ident,lbl
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lbl(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1       nx1(NS),kp(NS),kx1(NS),lx(NS),jx(NS),mx(NS),ny(NS),ky(NS),
     2       iy(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/phycon/alpha,pi,ev,bohr,ryd
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
     
      dimension vl(NGP),vs(NGP)
      dimension paj(NGP),qaj(NGP),uaj(NGP),vaj(NGP)
      dimension x1(NGP),x2(NGP),z1(NGP),z2(NGP),y(NGP)

      if(inf.eq.0.or.inf.eq.1) then
         ncore = jmax
      else
         ncore = in - 1
      endif

      do 100 i=1,max
         vl(i) = 0.0
         vs(i) = 0.0
 100  continue
      
      mj = mx(jj)
      kj = kx1(jj)
      lj = lx(jj)
      kpj = kp(jj)

      do 500 ja = 1,ncore
         ma = mx(ja)
         ka = kx1(ja)
         kpa = kp(ja)
         la = lx(ja)
         maj = min(ma,mj)
         do 120 i=1,maj
            uaj(i) = g(i,ja)*f(i,jj) - f(i,ja)*g(i,jj)
            vaj(i) = g(i,ja)*f(i,jj) + f(i,ja)*g(i,jj)
 120     continue

         llo = abs(kj-ka)
         lhi = kj+ka-1
         do 450 ll = llo,lhi
            angfac = cre(kj,ll,ka)**2/(2*kj)
            lt = la+lj+ll
            if(mod(lt,2).eq.0) then

*****  even parity branch
*****  n.b.   ll = 0  cancels between  magnetic and retardation

               if(ll.gt.0) then
                  cp = -angfac*ll*(ll+1)/dfloat((2*ll+1)*(2*ll-1))
                  cq = -angfac*ll*(ll+1)/dfloat((2*ll+1)*(2*ll+3))
                  cx = angfac*ll*(ll+1)/dfloat(2*(2*ll+1))
                  fp = (kpj-kpa)/dfloat(ll)
                  fq = (kpj-kpa)/dfloat(ll+1)
                  do 150 i = 1,maj
                      paj(i) =  uaj(i) + fp*vaj(i)
                      qaj(i) = -uaj(i) + fq*vaj(i)
 150              continue
                  call yfun(paj,x1,ll-1,maj,*901)
                  call yfun(qaj,x2,ll+1,maj,*901)
                  call wfun(paj,qaj,z1,z2,ll,maj,*901)

                  do 160 i = 1,ma
                     vl(i) = vl(i)
     &                     + ( cp*( 1.0 - fp)*x1(i)
     &                     +   cq*(-1.0 - fq)*x2(i)
     &                     +   cx*(-1.0 - fq)*z1(i)
     &                     +   cx*( 1.0 - fp)*z2(i) )*f(i,ja)
                     vs(i) = vs(i)
     &                     + ( cp*(-1.0 - fp)*x1(i)
     &                     +   cq*( 1.0 - fq)*x2(i)
     &                     +   cx*( 1.0 - fq)*z1(i)
     &                     +   cx*(-1.0 - fp)*z2(i) )*g(i,ja)

 160              continue
               endif
            else

*****  odd parity branch

               if(kpa+kpj.ne.0) then
                  call yfun(vaj,y,ll,maj,*901)
                  cy = angfac*(kpa+kpj)**2/dfloat(ll*(ll+1))
                  do 170 i = 1,ma
                     vl(i) = vl(i) + cy*y(i)*f(i,ja)
                     vs(i) = vs(i) + cy*y(i)*g(i,ja)
 170              continue
               endif
            endif
 450     continue
 500  continue
      return

 901  write(6,*) ' Error in yfun called in brpot'
      return
      end



      subroutine breit(ebr,brc)
***********************************************************************
*
*     Calculate breit correction 
*
***********************************************************************
      include "global.par"
      implicit doubleprecision(a-h,o-z)
      character*4 ident,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/lam/clam(5,10,5)
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension pba(NGP),qba(NGP),uba(NGP),vba(NGP)
      dimension x1(NGP),x2(NGP),y(NGP),z1(NGP),ebr(NS)

*********  calculate breit  energy

      if(inf.eq.0)  then 
         jcore = jmax
      elseif(inf.eq.2) then
         jcore = jmax-1
      else
         jcore = in-1
      endif

      brc = 0.0
      do 900 ja = 1,jmax
         ma = ms(ja)
         ka = k(ja)
         kpa = kap(ja)
         la = l(ja)
         kpa = kap(ja)
         ebr(ja) = 0.0

         do 800 jb = 1,jcore
            mb = ms(jb)
            kb = k(jb)
            kpb = kap(jb)
            lb = l(jb)
            mba = min(ma,mb)
            do 100 i=1,mba
               uba(i) = g(i,jb)*f(i,ja)-f(i,jb)*g(i,ja)
               vba(i) = g(i,jb)*f(i,ja)+f(i,jb)*g(i,ja)
 100        continue
            llo = iabs(ka - kb)
            lhi = kb + ka - 1
            tba = 0.0
            rba = 0.0
            do 700 ll = llo,lhi
               lt = la + lb + ll
               if(mod(lt,2).eq.0) then
                  if(ll.eq.0) then
                     fq = kpa-kpb
                     do 200 i = 1,mba 
                         qba(i) = -uba(i) + fq*vba(i)
 200                 continue
                     call yfun(qba,x2,1,mba,*901)
                     do 210 i = 1,mba
                        x2(i) = -qba(i)*x2(i)*rp(i)
 210                 continue
                     xl = rint(x2,1,mba,7,h)/3d0
                     zl = xl
                     tba = tba - xl*clam(ka,ll+1,kb)
                     rba = rba + zl*clam(ka,ll+1,kb)
                  else
                     fp = (kpa-kpb)/dfloat(ll)
                     fq = (kpa-kpb)/dfloat(ll+1)
                     do 250 i = 1,mba
                         pba(i) =  uba(i) + fp*vba(i)
                         qba(i) = -uba(i) + fq*vba(i)
 250                 continue
                     call yfun(pba,x1,ll-1,mba,*901)
                     call yfun(qba,x2,ll+1,mba,*901)
                     call zfun(pba,z1,ll,mba,*901)
                     do 260 i = 1,mba
                         x1(i) = -pba(i)*x1(i)*rp(i)
                         x2(i) = -qba(i)*x2(i)*rp(i)
                         z1(i) = -qba(i)*z1(i)*rp(i)
 260                 continue
                     xl1 = rint(x1,1,mba,7,h)
                     xl2 = rint(x2,1,mba,7,h)
                     zl1 = rint(z1,1,mba,7,h)
                     xl = ll*xl1/dfloat(2*ll-1)
     &               + (ll+1)*xl2/dfloat(2*ll+3)
                     zl = ( ll**2*xl1/dfloat(2*ll-1)
     &                +     (ll+1)**2*xl2/dfloat(2*ll+3)
     &                +     ll*(ll+1)*zl1 )/dfloat(2*ll+1)
                     tba = tba - xl*clam(ka,ll+1,kb)
                     rba = rba + zl*clam(ka,ll+1,kb)
                  endif
               else
                  if(kpb+kpa.ne.0) then
                    call yfun(vba,y,ll,mba,*901)
                    do 270 i = 1,mba
                       y(i) = y(i)*vba(i)*rp(i)
 270                continue
                    yl = -(kpb+kpa)**2*rint(y,1,mba,7,h)/
     &                         dfloat(ll*(ll+1))
                    tba = tba - yl*clam(ka,ll+1,kb)
                  endif
               endif
 700        continue
c           write(6,*) 'T(',lab(ja),lab(jb),')=',tba,rba
            ebr(ja) = ebr(ja) + jdg(jb)*(tba + rba)
 800     continue
         if(ja.le.jcore) brc = brc + 0.5*jdg(ja)*ebr(ja)
         if(inf.eq.2.and.ja.eq.jmax) brc = brc + ebr(ja)
 900  continue
      return
 901  write(6,*) ' Error in yfun called in Breit '
      return
      end

      subroutine valenc(xalpha)
      implicit doubleprecision(a-h,o-z)
****************************************************************
*
*     determine valence orbitals
*
****************************************************************
      include "global.par"
      character*4 ident,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/shell/rk(NGP),rl(NGP),dumm(16)
      common/charge/znuc(NGP),z(NGP),yhf(NGP),y(NGP),u(NGP),ysum(NGP)
      common/lam/clam(5,10,5)
     &      /phycon/alpha,pi,ev,bohr,ryd
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension zex(NGP),vk(NGP),vl(NGP)
      data eps/1d-7/,ieps/100/,th/0.333333333333333d0/

******   set up direct potential and starting df functions ******

      cex=xalpha*(3.0/(32.0*pi**2))**th
      do 10  i = 1,max 
  10      u(i)=0.0
      mex=0
      do 20  j = 1,jmax
           mm=ms(j)
           if(mm.gt.mex) mex=mm
           do 15  i = 1,mm
  15          u(i)=u(i)+jdg(j)*(g(i,j)**2+f(i,j)**2)
  20  continue
      call yfun(u,y,0,mex,*901)
      do 100 i = 1,max
           zex(i)=0.0
           if(i.gt.mex) go to 100
              zex(i)=cex*(r(i)*u(i))**th
  100      z(i)=znuc(i) - ysum(i)*r(i) + zex(i)

******** loop over all valence states

      j1=jmax+1
      do 500 jv = j1,in
         call master(jv)
         write(6,2050) -0.5*wh(jv),-0.5*wf(jv),kc(jv)
 2050    format(/5x,'valence electron iteration'/
     1    5x,'wh=',f14.7,4x,'wf=',f14.7,4x,'# iter=',i3)

******* exchange iteration loop   ********

      wh(jv)=wf(jv)
      mv=ms(jv)
      kv=k(jv)
      ddw=0.0

      do 120 i = 1,mv
          g(i,jv)=gh(i,jv)
 120      f(i,jv)=fh(i,jv)
      if(mv.lt.max) then
            mv1 = mv+1
            do 130 i = mv1,max
                  g(i,jv)=0.0
 130              f(i,jv)=0.0
      endif

      it=0
 200  continue

*******  calculate inner shell overlap integrals  *****
      noin=2
      do 230 j = 1,jmax
          if(kap(j).ne.kap(jv)) go to 230
          mm=min0(ms(j),mv)
          do 210 i = 1,mm
 210          u(i)=(g(i,j)*g(i,jv)+f(i,j)*f(i,jv))*rp(i)
          over=rint(u,1,mm,7,h)
          write(6,2100) lab(j),lab(jv),over
 2100     format('     <',a,'|',a,'> =', 1pe10.2)
          noin=noin+1
 230  continue

******  calculate valence norm   *****

      do 212 i = 1,mv
 212      u(i)=(g(i,jv)**2+f(i,jv)**2)*rp(i)
      over=1.0-rint(u,1,mv,7,h)
      write(6,2101) lab(jv),lab(jv),over
 2101 format(' 1.0-<',a,'|',a,'> =', 1pe10.2/)

********** setup rk and rl   ***********

      do 310 i = 1,mv
          rk(i) = -zex(i)*g(i,jv)
 310      rl(i) = -zex(i)*f(i,jv)
      do 320 i = mv1,max
          rk(i)=0.0
 320      rl(i)=0.0
      do 400 jl = 1,jmax
         kl=k(jl)
         ml=ms(jl)
         mlow=min0(ml,mv)
         do 340 i = 1,mlow
 340        u(i)=g(i,jl)*g(i,jv)+f(i,jl)*f(i,jv)
         llow=iabs(kv-kl)+1
         lhi=kv+kl
         do 360 ll = llow,lhi
            lt=l(jl)+l(jv)+ll+1
            if(mod(lt,2).ne.0) go to 360
            call yfun(u,y,ll-1,mlow,*901)
            cl=clam(kl,ll,kv)
            adg=cl*jdg(jl)
            do 350 i = 1,ml
                rk(i)=rk(i)+adg*y(i)*g(i,jl)*r(i)
                rl(i)=rl(i)+adg*y(i)*f(i,jl)*r(i)
 350       continue
 360     continue
 400  continue

*** add in breit

      call brpot(jv,vk,vl)
      do 365 i = 1,mv
           rk(i) = rk(i) - vk(i)*r(i)
           rl(i) = rl(i) - vl(i)*r(i)
 365  continue

      call solve(jv)
      del=dw(jv)-ddw
      ddw=dw(jv)
      it=it+1
      if(it.gt.ieps) go to 450
      write(6,2060) it,-0.5*wf(jv),0.5*del
 2060 format(4x,'iteration',i3,4x,' w(au) =',1pe14.6,4x,' dw =',e10.2)

      if(dabs(del).gt.eps) go to 200
 450  continue
 500  continue
      return
 901  stop
      end

      subroutine hart(exfac)
      implicit doubleprecision(a-h,o-z)
**********************************************************************
*
*   Routine to obtain approximate starting potential
*
**********************************************************************
      include "global.par"
      character*4 ident,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/charge/znuc(NGP),z(NGP),zn(NGP),rho(NGP),y(NGP),q(NGP)
      common/shell/g(NGP),f(NGP),gf(16)
     &      /phycon/alpha,pi,ev,bohr,ryd
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),gt(NGP,NS),ft(NGP,NS)
      
      eps=1.
      delt=.001
      icr=max
      sum=0.
      do 1 j=1,jmax
    1 sum=sum+jdg(j)
       
      exfac = exfac/sum
      if(exfac.eq.0.0) 
     *   exfac = 1.0/sum

	 ttf=exfac*sum
      write(6,1199) ttf
 1199 format(' exfac = ',f8.4)

      write(6,1000) ident
 1000 format(/'          Hartree iteration for ',a/)
      delw=.001
      h1=.5
      h2=.5
      mtest=-1
      ntest=0
   20 ntest=ntest+1
      epss=eps

      do 8 i=1,max
    8 q(i)=0.0
      do 2 j=1,jmax
      dgn=jdg(j)
      call master(j)
      m=ms(j)
      do 22 i=2,m
   22 q(i)=q(i)+dgn*y(i)
      if(m.eq.max) go to 2
      m1=m+1
      do 12 i=m1,max
   12 q(i)=q(i)+dgn
    2 continue
      do 3 i=1,max
      z1=znuc(i)-(1.0-exfac)*q(i)
      zn(i)=h1*z1+h2*z(i)
    3 continue
      eps=0.
      if(mtest.lt.1) go to 45
      do 40 i=1,max
   40 zn(i)=znuc(i)-(1.0-exfac)*q(i)
   45 continue
      do 7 j=1,jmax
      mj=ms(j)
      do 77 i=2,mj
      rrh=gh(i,j)**2+fh(i,j)**2
   77 y(i)=rrh*(z(i)-zn(i))*rpor(i)
      y(1)=0.
      dw(j)=2*rint(y,1,mj,7,h)
      eps=dmax1(eps,dabs(dw(j)/wf(j)))
      write(6,102) lab(j),-0.5*wh(j),-0.5*wf(j),kc(j)
  102 format('  w(',a,') =',2f15.5,6x,'k =',i3)
 7    wh(j)=wf(j)+dw(j)
      do 9 i=1,max
    9 z(i)=zn(i)
      write(6,103) ntest,eps
  103 format('  loop ',i3,'   rel-err =',1pe9.1/)
      if(ntest.gt.30) go to 10
      if(eps.gt.delt.and.eps.lt.epss) go to 20
   10 if(mtest) 23,24,25
   23 mtest=0
      delw=1.e-5
      ntest=0
      go to 20
   24 mtest=1
      delw=(10.0)**(-io)
      ntest=0
      go to 20
 25   continue
      write(6,105) ident
  105 format(10x,'Relativistic Hartree Energy Levels for ',a/)
      write(6,106)
  106 format(4x,5h jmax,5x,2h z,5x,4h atw,4x,4h nuc,4x,4h ion,4x,3h io,
     $5x,3h in,5x,4h inf)
      write(6,107) jmax,jz,nat,nuc,ion,io,in,inf
  107 format(8i8/)
      write(6,108)
  108 format(4x,6h shell,4x,2h n,5x,4h kap,4x,4h deg,4x,
     $ 13h atomic units/)
      do 30 j=1,jmax
      call master(j)
      wh(j)=wf(j)
      wpt=-wf(j)
      wau=.5*wpt
      write(6,109) lab(j),n(j),kap(j),jdg(j),wau
  109 format(4x,a,3i8,f16.6)
      do 29 i=1,max
      gt(i,j)=gh(i,j)
29    ft(i,j)=fh(i,j)
   30 continue
      do 31 i=1,max
   31 zn(i)=(1-exfac)*q(i)
      return
      end

      subroutine solve(j)
      implicit doubleprecision(a-h,o-z)
******************************************************************
*
*      Solve DHF equations for orbital jv using 
*           Fock's Green's function
*
******************************************************************
      include "global.par"
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/soln/gp(NGP),fp(NGP)
      common/drive/gk(NGP),fk(NGP),gfk(16)
      common/shell/rk(NGP),rl(NGP),dumm(16)
      common/charge/znuc(NGP),z(NGP),yhf(NGP),w(NGP),u(NGP),ysum(NGP)
      common/trash/s1(NGP),t1(NGP),s2(NGP),t2(NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension x(NGP),v(NGP),rkold(NGP),rlold(NGP)
      data jsave /0/

      m1=ms(j)
      w1=wh(j)
      k1=kap(j)
      if(j.ne.jsave) go to 60
      if(iof(j).eq.0) go to 70
      h2=iof(j)/10d0
      h1=1d0-h2
      do 75 i=1,m1
      rk(i)=h1*rk(i)+h2*rkold(i)
75    rl(i)=h1*rl(i)+h2*rlold(i)
      go to 70
60    continue
      do 1 i=1,max
      gk(i)=gh(i,j)
1     fk(i)=fh(i,j)
      do 2 i=1,16
2     gfk(i)=gfh(i,j)
      call outi(m1,k1,w1)
      do 3 i=1,max
      s1(i)=gp(i)
3     t1(i)=fp(i)
      call ini(m1,k1,w1)
      do 4 i=1,max
      s2(i)=gp(i)
4     t2(i)=fp(i)
70    jsave=j
      do 5 i=2,m1
      rkold(i)=rk(i)
      rlold(i)=rl(i)
      u(i)=(rk(i)*gh(i,j)+rl(i)*fh(i,j))*rpor(i)
5     w(i)=(g(i,j)*gh(i,j)+f(i,j)*fh(i,j))*rp(i)
      u(1)=0
      w(1)=0
      rkold(1)=0.0
      rlold(1)=0.0
      top=rint(u,1,m1,7,h)
      bot=rint(w,1,m1,7,h)
      dw(j)=-2.*top/bot
      wf(j)=wh(j)+dw(j)
      do 6 i=2,m1
      rk(i)=rk(i)/r(i)+.5*dw(j)*g(i,j)
6     rl(i)=rl(i)/r(i)+.5*dw(j)*f(i,j)
      rk(1)=0
      rl(1)=0
      do 7 i=1,m1
7     u(i)=(rk(i)*gh(i,j)+rl(i)*fh(i,j))*rp(i)
      call yint(u,u,v,w,m1,h)
      do 8 i=1,m1
      gk(i)=s1(i)*w(i)+s2(i)*v(i)
8     fk(i)=t1(i)*w(i)+t2(i)*v(i)
      do 9 i=1,m1
      u(i)=(s1(i)*rk(i)+t1(i)*rl(i))*rp(i)
9     w(i)=(s2(i)*rk(i)+t2(i)*rl(i))*rp(i)
      call yint(u,w,v,x,m1,h)
      do 10 i=1,m1
      yy=v(i)+x(i)
      gk(i)=(gk(i)+yy*gh(i,j))/alpha
10    fk(i)=(fk(i)+yy*fh(i,j))/alpha
      do 11 i=1,m1
      w(i)=(gk(i)**2+fk(i)**2)*rp(i)
11    u(i)=(gk(i)*gh(i,j)+fk(i)*fh(i,j))*rp(i)
      sp=rint(u,1,m1,7,h)
      tint=rint(w,1,m1,7,h)
      c=dsqrt(1d0-tint)
      do 12 i=1,m1
      u(i)=g(i,j)**2+f(i,j)**2
      g(i,j)=c*gh(i,j)-gk(i)
      f(i,j)=c*fh(i,j)-fk(i)
12    w(i)=g(i,j)**2+f(i,j)**2
      if(j.gt.jmax) return
      call yfun(u,v,0,m1,*901)
      call yfun(w,x,0,m1,*901)
      do 13 i=1,max
13    ysum(i)=ysum(i)+jdg(j)*(x(i)-v(i))
      return
 901  stop
      end

      subroutine pick
      implicit doubleprecision(a-h,o-z)
****************************************************************
*
*     Select orbital with largest error and solve DHF equations
* 
****************************************************************
      include "global.par"
      character*4 lab,ident
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      dimension ddw(NS),erorb(NS)

******   initialization   ***
      mpr = 0
      jup=100*jmax
      jtot=0
      emax=1.0d30
      write(6,100) ident
100   format(/4x,'Dirac-Fock iteration loops for ',a/)
      do 10 j=1,jmax
      jtot=jtot+1
      if(mpr.eq.0) 
     * write(6,101) lab(j),-0.5*wf(j),emax,jtot
101   format(2x,a,f16.8,4x,f10.1,4x,'initial ',i4)
      jtot=jtot+1
      call rside(j)
      call solve(j)
      ddw(j)=dw(j)
      write(6,101) lab(j),-0.5*wf(j),emax,jtot
      call rside(j)
      call solve(j)
10    continue

******   set up convergence criterion   ***

      eps=1d0-7
      if(io.eq.0) go to 20
      nx1=io
      if(io.lt.6) nx1=6
      if(io.gt.9) nx1=9
      eps=10.d0**(-nx1)
20    jtot=jtot+1

******   choose orbital with largest error   ***

      do 30 j=1,jmax
30    erorb(j)=dabs((ddw(j)-dw(j))/wf(j))
      j=1
      emax=erorb(1)
      do 40 jj=2,jmax
      if(erorb(jj).le.emax) go to 40
      j=jj
      emax=erorb(jj)
40    continue
      if(emax.le.eps) go to 70
      write(6,102) lab(j),-0.5*wf(j),emax,jtot
102   format(2x,a,f16.8,4x,1pe10.1,4x,'iterate ',i4)
      ddw(j)=dw(j)
      call rside(j)
      call solve(j)
      if(jtot.ge.jup) return
      go to 20

******   final convergence test   ***

70    jtot=jtot-1
      elim=emax
      do 60 jj=1,jmax
      jtot=jtot+1
      ddw(jj)=dw(jj)
      call rside(jj)
      call solve(jj)
      emax=dabs((ddw(jj)-dw(jj))/wf(jj))
      if(emax.gt.elim) elim=emax
      write(6,103) lab(jj),-0.5*wf(jj),emax,jtot
103   format(2x,a,f16.8,4x,1pe10.1,4x,'final   ',i4)
60    continue
      if(elim.gt.eps) go to 20
      return
      end

      subroutine potnew(jz,c,a,b2,b4,rms,en,z,inuc)
      include "global.par"
      parameter(NR=5,NT=10)
      implicit doubleprecision(a-h,o-z)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
     &      /phycon/alpha,pi,ev,bohr,ryd
     &      /nucpar/aa,bb2,bb4,cc
      dimension xt(NT),wt(NT),xr(NR),wr(NR)
      dimension g(NGP),f(NGP),u(NGP),z(NGP)
      conv = 1.0d5*bohr
      aa = a/conv
      bn2 = 0.25d0*dsqrt(5.0d0/pi)
      bb2 = bn2*b2
      bn4 = 3.d0/(16.d0*dsqrt(pi))
      bb4 = bn4*b4
      cc = c/conv
      cau = cc
      do 100 i = 1,max
         g(i) = 0.0
         f(i) = 0.0
         u(i) = 0.0
         if(r(i).lt.cau) inuc = i
         if(r(i).lt.4*cau) im = i
 100  continue
      im = im + 1

      call setgau(xt,wt,NT)
      call setgau(xr,wr,NR)

      do 300 ir = 1,im-1

         r1 = r(ir)
         r2 = r(ir+1)
         dr = 0.5*(r2-r1)
         rav = 0.5*(r2+r1)
         do 200 irr = 1,NR
            rr = dr*xr(irr) + rav

            rhoav = 0.0
            do 150 it = 1,NT
               rhoav = rhoav + 0.5*wt(it)*rho(rr,0.5*(xt(it)+1.0d0))
 150        continue

            del = dr*wr(irr)*rr*rhoav
            g(ir+1) = g(ir+1) + rr*del
            u(ir+1) = u(ir+1) + rr*rr*rr*del
            f(ir) = f(ir) + del
 200     continue
 300  continue
      do 350 i = 2,max
          g(i) = g(i) + g(i-1)
          u(i) = u(i) + u(i-1)
 350  continue
      do 360 i = max,2,-1
          f(i-1) = f(i-1) + f(i)
 360  continue

      gnorm = g(max)
      unorm = u(max)
      znorm = dfloat(jz)/gnorm
      c3 = cc**3/3.0d0
      c5 = cc**5/5.0d0
      en = gnorm/c3
      em = unorm/c5
      rms = c*dsqrt(0.6d0*em/en)

      z(1) = 0.0
      do 400 i = 2,max
         z(i) = znorm*( g(i)+r(i)*f(i) )
 400  continue 
      return
      end

      function rho(r,x)
      implicit doubleprecision(a-h,o-z)
      common/nucpar/a,b2,b4,c

      rr = c*( 1d0 + b2*( 3.0d0*x*x - 1.0d0 ) +
     &               b4*( ( 35d0*x*x-30.0d0 )*x*x + 3.0d0 ) )
      arg = (r-rr)/a
      if(arg.le.80.d0) then
           rho = 1.0d0/( 1.0d0 + dexp(arg) )
      else
           rho  = 0.0d0
      endif

      return
      end

      subroutine setgau(x,w,n)
      implicit doubleprecision(a-h,o-z)
      dimension x(n),w(n)
      data eps /1d-15/
      pi = dacos(-1.0d0)
      m = (n+1)/2
      do 120 i = 1,m
         z = cos(pi*(i-0.25d0)/(n+0.5d0))
 100     continue
            p1 = 1.d0
            p2 = 0.d0
            do 110 j = 1,n
               p3 = p2
               p2 = p1
               p1 = (( 2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 110        continue
            pp = n*(z*p1-p2)/(z*z-1.d0)
            z1 = z
            z = z1-p1/pp
         if(dabs(z-z1).gt.eps) go to 100
         x(i) = -z
         x(n+1-i) = z
         w(i) = 2.d0/((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
 120  continue
      return
      end

      subroutine potent(jz,ia,rnuc,cnuc,tnuc,en,z,inuc,*)
      implicit doubleprecision(a-h,o-z)
c ---------------------------------------------------------------------
c  subroutine to determine the nuclear potential
c
c  rnuc  =  rms radius (f) (input or output)
c  cnuc  =  radial parameter in fermi distribution (f) (input or output)
c  tnuc  =  nuclear skin thickness parameter (f) (input only)
c  z(i),i=1,max :  output charge distribution
c  inuc  =       # grid points with r(i) .le. cnuc
c          ----------------------------------------------
c  ia=0; rnuc=cnuc=tnuc=0    : coulomb field
c  rnuc = cnuc = tnuc = 0    : empirical uniform charge distribution
c  rnuc .ne. 0               : uniform distribution with rms rad = rnuc
c  cnuc .ne. 0 , tnuc = 0    : uniforn distribution with radius  = cnuc
c  cnuc .ne. 0 , tnuc .ne. 0 : fermi distribution
c
c ---------------------------------------------------------------------
      include "global.par"
      common/phycon/alpha,pi,ev,bohr,ryd
      dimension z(NGP)
c
      if(rnuc.ne.0.0) go to 300
      if(cnuc.ne.0.0) go to 200
      if(tnuc.eq.0.0) go to 100

      write(6,*) ' **** Error in potent: cnuc.eq.0 and tnuc.ne.0'
      go to 901
c
c  empirical default :  rnuc=cnuc=tnuc=0
c
 100  continue
      if(ia.ne.0) go to 110
      rnuc = 0.0
      go to 150

*  Johnson and Soff ADNDT 33, 410 (1985)
 110  rnuc = 0.836d0*dfloat(ia)**(1.d0/3.d0) + 0.570d0

 150  cnuc = dsqrt(5.d0/3.d0)*rnuc
      c = 1.d-5*cnuc/bohr
      t = 0.0
      call zfill(jz,c,t,en,rms,z,inuc)
      return
c
c   rnuc = 0 , cnuc.ne.0  :  fermi distribution
c
 200  continue
      c = 1.d-5*cnuc/bohr
      t = 1.d-5*tnuc/(4.d0*dlog(3.d0)*bohr)
      call zfill(jz,c,t,en,rms,z,inuc)
      rnuc = 1.d5*rms*bohr
      return
c
c  cnuc = 0  rnuc.ne.0  uniform distribution
c
 300  continue
      if((cnuc.ne.0.0).or.(tnuc.ne.0.0)) then
         write(6,*) 
     &   ' **** Error in potent: rnuc.ne.0 and tnuc or cnuc given'
      endif
      go to 150
 901  return 1
      end

      subroutine aitsum(sum,fac,n)
      implicit doubleprecision(a-h,o-z)
c ---------------------------------------------------------------------
c    routine to calculate sum from 1 to inf. of the series
c              (-1)**(l-1)*fac**l/l**n
c           using aitken's acceleration scheme
c    restrictions  :  n.ge.1   fac.le.1
c ---------------------------------------------------------------------

      dimension x(3,4),w(3,4),xx(4),ww(4)
      data xx /0.322547689619392d0,1.745761101158346d0,
     &         4.536620296921128d0,9.395070912301134d0/
      data ww /0.603154104341634d0,0.357418692437800d0,
     &         0.038887908515005d0,0.000539294705561d0/
      data x/.935822227524088d0,3.305407289332278d0,7.758770483143634d0,
     &     1.517387080677412d0,4.311583133719520d0, 9.171029785603068d0,
     &     2.141216276717724d0,5.315517126176787d0,10.54326659710549d0,
     &     2.796496075793495d0,6.318244090839835d0,11.88525983336667d0/
      data w/.588681481039660d0,0.391216059222310d0,0.020102459738030d0,
     &      1.037494961490425d0,0.905750004703065d0,0.056755033806510d0,
     &      2.836328204670725d0,2.951294312715705d0,0.212377482613570d0,
     &     10.55459043765716d0,12.45156474508734d0, 0.99384481725550d0/
      data nmax/99/,eps/1d-15/
c
      s1=1.0
      s2=s1-fac/dfloat(2)**n
      sac=fac**2
      s3=s2+sac/dfloat(3)**n
      if(dabs(s2-s3).gt.eps*dabs(s3)) go to 20
      sum=fac*s3
      return
 20   den=s1+s3-2*s2
      s2a=s3-(s3-s2)**2/den
      do 100 l=4,nmax
      s1=s2
      s2=s3
      sac=-fac*sac
      s3=s2+sac/dfloat(l)**n
      if(dabs(s3-s2).gt.eps*dabs(s3)) go to 50
      sum=s3*fac
      return
 50   sa1=sa2
      den=s1+s3-2*s2
      if(dabs(den).gt.eps*dabs(s3)) go to 90
      write(6,*) ' ***** Error in aitken: problems with convergence'
      write(6,1010) n,fac,l,del,sa2
 1010 format(' n=',i4,4x,'fac=',d13.6,4x,'l=',i4,4x,'del=',d10.2,4x,
     & 'sa2=',d10.2)
      go to 120
 90   sa2=s3-(s3-s2)**2/den
      del=dabs(sa2-sa1)
      if(del.lt.eps*dabs(sa2)) go to 120
 100  continue
      m=nmax+1
      k=n-1
      if(k.lt.0.or.k.gt.4) go to 105
      fw=0.0
      sw=0.0
      if(k.ne.0) go to 102
      do 101 j=1,4
      fw=fw+ww(j)/(1.d0+fac*dexp(-xx(j)/m))
 101  sw=sw+ww(j)
      go to 104
 102  do 103 j=1,3
      fw=fw+w(j,k)/(1.0+fac*dexp(-x(j,k)/m))
 103  sw=sw+w(j,k)
 104  rem=(-fac)**nmax*fw/(sw*dfloat(m)**n)
      sum=fac*(s3+rem)
      return
 105  write(6,*)
     & ' **** Error in aitken: series not converged after 100 steps'
      sum=fac*sa2
      write(6,1000) fac,n,del,sum
 1000 format(' fac=',d14.8,4x,'n=',i6,4x,'del=',d10.2,4x,'sum=',d20.13)
      return
 120  sum=fac*sa2
      return
      end

      subroutine zfill(jz,c,t,en,rms,z,inuc)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
c ---------------------------------------------------------------------
c    routine to fill in the potential for a uniform or fermi charge
c                     distribution
c ---------------------------------------------------------------------
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
     &      /phycon/alpha,pi,ev,bohr,ryd
      dimension z(NGP)
      data cut/36./
      if(t.eq.0.0) go to 200
      a=c/t
      if(a.lt.cut) go to 50
      fac=0.0
      s3=0.0
      s5=0.0
      go to 60
 50   fac=dexp(-a)
      call aitsum(s3,fac,3)
      call aitsum(s5,fac,5)
 60   en=1.0+(pi/a)**2+6.0*s3/a**3
      em=1.0+10.0*(pi/a)**2/3.0+7.0*(pi/a)**4/3.0+120.0*s5/a**5
      rms=c*dsqrt(3.d0*em/(5.d0*en))
      z(1)=0.0
      inuc=1
      do 100 i=2,max
      b=dabs(r(i)-c)/t
      if(b.lt.cut) go to 70
      sac=0.0
      p2=0.0
      p3=0.0
      go to 80
 70   sac=dexp(-b)
      call aitsum(p2,sac,2)
      call aitsum(p3,sac,3)
 80   if(r(i).lt.c) go to 90
      z(i)=jz*(1.0+(pi/a)**2+6.0*(s3-p3)/a**3
     & -3.0*r(i)*p2/(c*a**2))/en
      go to 100
 90   z(i)=jz*(r(i)*(1.5-0.5*(r(i)/c)**2+0.5*(pi/a)**2+
     & 3.0*p2/a**2)+6.0*t*(s3-p3)/a**2)/(en*c)
      inuc=i
 100  continue
      return
 200  inuc=1
      do 300 i=1,max
      if(r(i).ge.c) go to 250
      z(i)=jz*r(i)*(1.5-0.5*(r(i)/c)**2)/c
      inuc=i
      go to 300
 250  z(i)=jz
 300  continue
      en = 1.0
      rms=dsqrt(0.6d0)*c
      return
      end

      subroutine master(j)
      implicit doubleprecision(a-h,o-z)
      character*4 lab,ident
      include "global.par"
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common/shell/g(NGP),f(NGP),a0(7),g0,b0(7),f0
      common/charge/znuc(NGP),z(NGP),zn(NGP),rho(NGP),y(NGP),dumm(NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),gt(NGP,NS),ft(NGP,NS)
      dimension u(NGP)
      data delw/1.e-9/
      
      kp=kap(j)
      ka=k(j)
      na=n(j)
      wmore=0.
      wless=0.
      more=0
      less=0
      ntf=na-ka
      ntg=ntf
      if(kp.gt.0) ntg=ntg-1
      w=wh(j)
      kc(j)=0
    1 kc(j)=kc(j)+1
      if(kc(j).lt.30)  go to 101
      del=0.0
      if(kc(j).gt.50) go to 25
      write(6,8000) na,kp,kc(j)
 8000 format(/4x,'n=',i4,4x,'kap=',i4,4x,'kc=',i4)
      write(6,8001) w,wmore,wless
 8001 format(4x,'w=',f14.8,4x,'wmore=',f14.8,4x,'wless=',f14.8)
  101 continue
      m1=max+1
    2 m1=m1-1
      if((w*r(m1)+2.*z(m1))*r(m1)+1600.) 2,2,3
    3 m=m1
      ms(j)=m
      m1=m1+1
    4 m1=m1-1
      if(w*r(m1)+2.*z(m1)) 4,5,5
    5 nc=m1
      if(kc(j).lt.30)  go to 501
      write(6,8002) nc,m
 8002 format(3x,'nc=',i4,4x,'  m=',i4)
 501  continue
      call inh(nc,m,kp,w)
      fn=f(nc)
      gn=g(nc)
      call outh(nc,kp,w)
      rat=g(nc)/gn
      fn=rat*fn
      m1=nc+1
      do 6 i=m1,m
      g(i)=rat*g(i)
    6 f(i)=rat*f(i)
      nzf=0
      nzg=0
      sg=sign(1.0d0,g(2))
      sf=sign(1.0d0,f(2))
      do 7 i=3,m
      if(sign(1.0d0,g(i)).eq.sg) go to 8
      sg=-sg
      nzg=nzg+1
    8 if(sign(1.0d0,f(i)).eq.sf) go to 7
      sf=-sf
      nzf=nzf+1
    7 continue
      if(kc(j).le.30) go to 701
      write(6,8003) ntf,nzf,ntg,nzg
 8003 format(3x,'ntf,nzf,ntg,nzg=',4i6)
 701  continue
      if(nzf-ntf) 9,10,11
   10 if(nzg-ntg) 9,20,11
    9 less=less+1
      if((less.eq.1).or.(w.gt.wless)) wless=w
      if(more.eq.0) go to 12
   13 w=.5*(wmore+wless)
      go to 1
   12 w=.90*wless
      go to 1
   11 more=more+1
      if((more.eq.1).or.(w.lt.wmore)) wmore=w
      if(less.ne.0) go to 13
      w=1.10*wmore
      go to 1
   20 do 19 i=1,m
      rho(i) = g(i)**2+f(i)**2
   19 u(i)=rho(i)*rp(i)
      s=0.5*rint(u,1,m,7,h)/alpha
      del=g(nc)*(fn-f(nc))/s
      w=w+del
      if((less.ne.0).and.(w.lt.wless)) go to 22
      if((more.ne.0).and.(w.gt.wmore)) go to 23
      if(w.lt.0) go to 24
      w=1.5*(w-del)
      go to 1
   22 w=.5*(w-del+wless)
      go to 1
   23 w=.5*(w-del+wmore)
      go to 1
   24 if(dabs(del/w)-delw) 25,25,1
   25 wf(j)=w-del
      nct(j)=nc
      call yfun(rho,y,0,m,*901)
      an2=1.0/rint(u,1,m,7,h)
      an1=dsqrt(an2)
      do 30 i=1,m
      fh(i,j)=an1*f(i)
      gh(i,j)=an1*g(i)
   30 y(i) = an2*y(i)*r(i)
      do 31 i=1,7
      gfh(i,j)=a0(i)
31    gfh(i+8,j)=b0(i)
      gfh(8,j)=an1*g0
      gfh(16,j)=an1*f0
      return
 901  stop
      end

      subroutine inh(ncl,max,kap,w)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,lim
      common/shell/g(NGP),f(NGP),a0(8),b0(8)
      common/charge/znuc(NGP),z(NGP),dumm(4*NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd
      dimension a(8),gp(7),fp(7)
      data a,dd/  36799.d0, 139849.d0,-121797.d0, 123133.d0,
     & -88547.d0, 41499.d0, -11351.d0,   1375.d0, 120960.d0/

      con = 2*alpha
      lhi=max+6
      low=max
      if(lhi.le.lim) go to 202
      write(6,2020) max
 2020 format(/' **** warning  max =',i5,' in inh *****')
      lhi=lim
      low=lim-6
 202  continue
      f1=dsqrt(1.0+w/con**2)
      f2=dsqrt(-w)
      alam=f1*f2
      eps=1.+2.*w/con**2
      z0=z(max+3)
      sig=eps*z0/alam
      a0(1)=.5*(alam*kap+z0)
      g2=kap**2-(z0/alpha)**2
      g1=kap-z0*alam/alpha**2
      do 40 k=1,7
      if(k.eq.7) go to 40
      a0(k+1)=.5*(g2-(k-sig)**2)*a0(k)/k
   40 b0(k)=a0(k)*(g1+(k-sig)*eps)/(k*alam)
      do 10 i=low,lhi
      k=lhi+1-i
      rho=alam*r(i)
      fac=rho**sig*dexp(-rho)
      af=0.d0
      bf=1.d0
      rhn=1.d0
      do 11 kn=1,7
      rhn=rhn*rho
      af=af+a0(kn)/rhn
   11 bf=bf+b0(kn)/rhn
      g(i)=fac*(f1*bf+2*f2*af/con**2)
      f(i)=fac*(.5*f2*bf-f1*af)/alpha
      x=(r(i)*w+2.d0*z(i))/con
      gp(k)=-rpor(i)*(kap*g(i)+(con*r(i)+x)*f(i))
   10 fp(k)=rpor(i)*(kap*f(i)+x*g(i))
      u=g(low)
      v=f(low)
      lop=low-ncl
      d8 = dd/h
      do 20 k=1,lop
      i=low-k
      all=rpor(i)*kap
      akk=-all
      x=(r(i)*w+2.*z(i))/con
      amm=-rpor(i)*(con*r(i)+x)
      ann=rpor(i)*x
      brk=1.d0+a(1)*akk/d8
      brl=1.d0+a(1)*all/d8
      brn=-a(1)*ann/d8
      brm=-a(1)*amm/d8
      del=brk*brl-brm*brn
      rk=v-(a(2)*fp(7)+a(3)*fp(6)+a(4)*fp(5)+
     1a(5)*fp(4)+a(6)*fp(3)+a(7)*fp(2)+a(8)*fp(1))/d8
      sk=u-(a(2)*gp(7)+a(3)*gp(6)+a(4)*gp(5)+
     1a(5)*gp(4)+a(6)*gp(3)+a(7)*gp(2)+a(8)*gp(1))/d8
      u=(brl*sk+brm*rk)/del
      v=(brk*rk+brn*sk)/del
      g(i)=u
      f(i)=v
      do 30 j=1,6
      gp(j)=gp(j+1)
   30 fp(j)=fp(j+1)
      gp(7)=akk*g(i)+amm*f(i)
   20 fp(7)=all*f(i)+ann*g(i)
      return
      end

      subroutine outh(max,kap,w)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      parameter(NDI=5,LWK=15)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,lim
      common/shell/g(NGP),f(NGP),a0(7),g0,b0(7),f0
      common/charge/znuc(NGP),z(NGP),dumm(4*NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd
      dimension aa(NDI,NDI),bb(NDI),d(NDI,NDI),e(NDI,NDI),ak(NDI),
     1 al(NDI),am(NDI),an(NDI),tf(NDI),tg(NDI),p(NDI),q(NDI),s(NDI),
     2 a(8),gp(7),fp(7)
      dimension ipiv(NDI),work(LWK)
      data a,dd/  36799.d0, 139849.d0,-121797.d0, 123133.d0,
     & -88547.d0, 41499.d0, -11351.d0,   1375.d0, 120960.d0/
      data bb,aa/ -24.d0,   6.d0,  -4.d0,   6.d0, -24.d0,
     1           -130.d0, -60.d0,  30.d0, -40.d0, 150.d0,
     2            240.d0, -40.d0,-120.d0, 120.d0,-400.d0,
     3           -120.d0, 120.d0,  40.d0,-240.d0, 600.d0,
     4             40.d0, -30.d0,  60.d0, 130.d0,-600.d0,
     5             -6.d0,   4.d0,  -6.d0,  24.d0, 274.d0/
      data nl/4/
      con = 2*alpha
      g(1)=0.
      f(1)=0.
      az=2.*z(1)/con
      gam=dsqrt(kap**2-az**2)
      if(kap.gt.0) go to 25
      g0=1.
      f0=az/(gam-kap)
      go to 40
   25 f0=-1.
      g0=az/(gam+kap)
   40 gg=g0
      ff=f0
      do 10 k=1,nl
      it=5*(k-1)+1
      do 1 i=1,5
      fac=120.*h*rpor(i+it)
      ak(i)=-fac*(kap+gam)
      al(i)=fac*(kap-gam)
      x=(r(i+it)*w+2.*z(i+it))/con
      am(i)=-fac*(x+con*r(i+it))
      an(i)=fac*x
      p(i)=-bb(i)*gg
    1 q(i)=-bb(i)*ff
      do 2 i=1,5
      do 2 j=1,5
      d(i,j)=aa(i,j)
      e(i,j)=aa(i,j)
      if(j.ne.i) go to 2
      d(i,i)=d(i,i)-ak(i)
      e(i,i)=e(i,i)-al(i)
    2 continue
      call dgetrf(NDI,NDI,d,NDI,ipiv,info)
      call dgetri(NDI,d,NDI,ipiv,work,LWK,info)
*     call dminv(d,5,det,l,m)
      do 3 i=1,5
      s(i)=0.
      do 3 j=1,5
      e(i,j)=e(i,j)-an(i)*d(i,j)*am(j)
    3 s(i)=s(i)+d(i,j)*p(j)
      call dgetrf(NDI,NDI,e,NDI,ipiv,info)
      call dgetri(NDI,e,NDI,ipiv,work,LWK,info)
*     call dminv(e,5,det,l,m)
      do 4 i=1,5
      tf(i)=0.
      do  4 j=1,5
    4 tf(i)=tf(i)+e(i,j)*(an(j)*s(j)+q(j))
      do 5 i=1,5
      tg(i)=s(i)
      do 5 j=1,5
    5 tg(i)=tg(i)+d(i,j)*am(j)*tf(j)
      do 6 i=1,5
      rsig=r(i+it)**gam
      g(i+it)=rsig*tg(i)
    6 f(i+it)=rsig*tf(i)
      gg=tg(5)
   10 ff=tf(5)
      d8 = dd/h
      loop=5*nl+1
      imin=loop-6
      do 80 i=imin,loop
      k=i+1-imin
      x=(r(i)*w+2.*z(i))/con
      gp(k)=-rpor(i)*(kap*g(i)+(con*r(i)+x)*f(i))
   80 fp(k)=rpor(i)*(kap*f(i)+x*g(i))
      u=g(loop)
      v=f(loop)
      min=loop+1
      do 81 i=min,max
      all=rpor(i)*kap
      akk=-all
      x=(r(i)*w+2.*z(i))/con
      amm=-rpor(i)*(con*r(i)+x)
      ann=rpor(i)*x
      brk=1.0-a(1)*akk/d8
      brl=1.0-a(1)*all/d8
      brn=a(1)*ann/d8
      brm=a(1)*amm/d8
      del=brk*brl-brm*brn
      rk=v+(a(2)*fp(7)+a(3)*fp(6)+a(4)*fp(5)+
     1a(5)*fp(4)+a(6)*fp(3)+a(7)*fp(2)+a(8)*fp(1))/d8
      sk=u+(a(2)*gp(7)+a(3)*gp(6)+a(4)*gp(5)+
     1a(5)*gp(4)+a(6)*gp(3)+a(7)*gp(2)+a(8)*gp(1))/d8
      u=(brl*sk+brm*rk)/del
      v=(brk*rk+brn*sk)/del
      g(i)=u
      f(i)=v
      do 85 j=1,6
      gp(j)=gp(j+1)
   85 fp(j)=fp(j+1)
      gp(7)=akk*g(i)+amm*f(i)
   81 fp(7)=all*f(i)+ann*g(i)
      return
      end

      subroutine ini(max,kap,w)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
c   ***   routine for fock green  function   ***
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,lim
      common/soln/g(NGP),f(NGP)
      common/drive/rg(NGP),rf(NGP),a0(7),g0,b0(7),f0
      common/charge/znuc(NGP),z(NGP),dumm(4*NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd      
      dimension a(8),gp(7),fp(7)
      dimension e0(7),d0(7)
      data a,dd/  36799.d0, 139849.d0,-121797.d0, 123133.d0,
     & -88547.d0, 41499.d0, -11351.d0,   1375.d0, 120960.d0/

      con = 2*alpha
      lhi=max+6
      low=max
      lave=max+3
      if(lhi.le.lim)go to 46
      lhi=lim
      low=lim-6
      lave=lim-3
      write(6,*)
     & ' **** Error in ini: scale too small to accomodate orbital'
 46   continue
      f1=dsqrt(1.0+w/con**2)
      f2=dsqrt(-w)
      alam=f1*f2
      eps=1+2*w/con**2
      z0=z(lave)
      sig=eps*z0/alam
      alc=z0/alam**2
      g1=kap-z0*alam/alpha**2
      g2=alam/alpha**2
      g3=alam*kap+z0
      g4=sig+1-eps*kap
      d0(1)=0.0
      e0(1)=.5*(sig+1+eps*kap)
      ami=0.0
      do 40 k=1,7
      if(k.eq.1) go to 41
      d0(k)=((g1+(k-1-sig)*eps)*e0(k)-eps*b0(k)-g2*a0(k)
     1 +alc*(alam*b0(k-1)-eps*a0(k-1)))/((k-1)*alam)
   41 if(k.eq.7) go to 40
      e0(k+1)=.5*(g3*d0(k)+(g4-k)*e0(k)+b0(k)+ami)
   40 ami=alc*a0(k)
      d8 = dd/h 
      do 10 i=low,lhi
      k=lhi+1-i
      rho=alam*r(i)
      fac=rho**sig*dexp(-rho)
      asum=0.
      bsum=1.
      dsum=eps/alam
      esum=0.
      rlog=dlog(rho)*alc
      rhn=1.
      do 11 kn=1,7
      rhn=rhn*rho
      asum=asum+a0(kn)/rhn
      bsum=bsum+b0(kn)/rhn
      dsum=dsum+d0(kn)/rhn
   11 esum=esum+e0(kn)/rhn
      si=rho*fac*(f1*dsum+2*f2*esum/con**2)
      ti=rho*fac*(.5*f2*dsum-f1*esum)/alpha
      gi=fac*(f1*bsum+2*f2*asum/con**2)
      fi=fac*(.5*f2*bsum-f1*asum)/alpha
      if(i.eq.low) eta=rg(i)/gi
      anorm=alpha*eta/alam
      g(i)=anorm*(rlog*gi+si)
      f(i)=anorm*(rlog*fi+ti)
      x=(r(i)*w+2.*z(i))/con
      gp(k)=-rpor(i)*(kap*g(i)+(con*r(i)+x)*f(i)+r(i)*eta*fi)
   10 fp(k)=rpor(i)*(kap*f(i)+x*g(i)+r(i)*eta*gi)
      u=g(low)
      v=f(low)
      lop=low-2
      g(1)=0.
      f(1)=0.
      do 20 k=1,lop
      i=low-k
      all=rpor(i)*kap
      akk=-all
      x=(r(i)*w+2.*z(i))/con
      amm=-rpor(i)*(con*r(i)+x)
      ann=rpor(i)*x
      ax=-rp(i)*rf(i)
      ay=rp(i)*rg(i)
      brk=1.0+a(1)*akk/d8
      brl=1.0+a(1)*all/d8
      brn=-a(1)*ann/d8
      brm=-a(1)*amm/d8
      del=brk*brl-brm*brn
      rk=a(1)*ay+a(2)*fp(7)+a(3)*fp(6)+a(4)*fp(5)+
     1a(5)*fp(4)+a(6)*fp(3)+a(7)*fp(2)+a(8)*fp(1)
      sk=a(1)*ax+a(2)*gp(7)+a(3)*gp(6)+a(4)*gp(5)+
     1a(5)*gp(4)+a(6)*gp(3)+a(7)*gp(2)+a(8)*gp(1)
      rk=v-rk/d8
      sk=u-sk/d8
      u=(brl*sk+brm*rk)/del
      v=(brk*rk+brn*sk)/del
      g(i)=u
      f(i)=v
      do 30 j=1,6
      gp(j)=gp(j+1)
   30 fp(j)=fp(j+1)
      gp(7)=akk*g(i)+amm*f(i)+ax
   20 fp(7)=all*f(i)+ann*g(i)+ay
      return
      end

      subroutine outi(max,kap,w)
      implicit doubleprecision(a-h,o-z)
      include "global.par"
      parameter(NDI=5,LWK=15)
c   ***   routine for fock green  function   ***
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,lim
      common/soln/g(NGP),f(NGP)
      common/drive/rg(NGP),rf(NGP),a0(7),g0,b0(7),f0
      common/charge/znuc(NGP),z(NGP),dumm(4*NGP)
     &      /phycon/alpha,pi,ev,bohr,ryd
      dimension aa(NDI,NDI),bb(NDI),d(NDI,NDI),e(NDI,NDI),ak(NDI),
     1 al(NDI),am(NDI),an(NDI),tf(NDI),tg(NDI),p(NDI),q(NDI),s(NDI),
     2 a(8),gp(7),fp(7)
      dimension ipiv(NDI),work(LWK)
      data a,dd/  36799.d0, 139849.d0,-121797.d0, 123133.d0,
     & -88547.d0, 41499.d0, -11351.d0,   1375.d0, 120960.d0/
      data bb,aa/ -24.d0,   6.d0,  -4.d0,   6.d0, -24.d0,
     1           -130.d0, -60.d0,  30.d0, -40.d0, 150.d0,
     2            240.d0, -40.d0,-120.d0, 120.d0,-400.d0,
     3           -120.d0, 120.d0,  40.d0,-240.d0, 600.d0,
     4             40.d0, -30.d0,  60.d0, 130.d0,-600.d0,
     5             -6.d0,   4.d0,  -6.d0,  24.d0, 274.d0/
      data nl/4/
      con = 2*alpha
      g(1)=0.
      f(1)=0.
      az=z(1)/alpha
      gam=dsqrt(kap**2-az**2)
      gam=gam+1.
      dn=2*gam-1
      gg=-(dn-2*kap)*f0/dn
      ff=(dn+2*kap)*g0/dn
      do 10 k=1,nl
      it=5*(k-1)+1
      do 1 i=1,5
      fac=120.*h*rpor(i+it)
      ak(i)=-fac*(kap+gam)
      al(i)=fac*(kap-gam)
      x=(r(i+it)*w+2.*z(i+it))/con
      am(i)=-fac*(x+con*r(i+it))
      an(i)=fac*x
      fac=fac/r(i+it)**(gam-1)
      p(i)=-fac*rf(i+it)-bb(i)*gg
    1 q(i)=fac*rg(i+it)-bb(i)*ff
      do 2 i=1,5
      do 2 j=1,5
      d(i,j)=aa(i,j)
      e(i,j)=aa(i,j)
      if(j.ne.i) go to 2
      d(i,i)=d(i,i)-ak(i)
      e(i,i)=e(i,i)-al(i)
    2 continue
      call dgetrf(NDI,NDI,d,NDI,ipiv,info)
      call dgetri(NDI,d,NDI,ipiv,work,LWK,info)
*     call dminv(d,5,det,l,m)
      do 3 i=1,5
      s(i)=0.
      do 3 j=1,5
      e(i,j)=e(i,j)-an(i)*d(i,j)*am(j)
    3 s(i)=s(i)+d(i,j)*p(j)
      call dgetrf(NDI,NDI,e,NDI,ipiv,info)
      call dgetri(NDI,e,NDI,ipiv,work,LWK,info)
*     call dminv(e,5,det,l,m)
      do 4 i=1,5
      tf(i)=0.
      do  4 j=1,5
    4 tf(i)=tf(i)+e(i,j)*(an(j)*s(j)+q(j))
      do 5 i=1,5
      tg(i)=s(i)
      do 5 j=1,5
    5 tg(i)=tg(i)+d(i,j)*am(j)*tf(j)
      do 6 i=1,5
      rsig=r(i+it)**gam
      g(i+it)=rsig*tg(i)
    6 f(i+it)=rsig*tf(i)
      gg=tg(5)
   10 ff=tf(5)
      loop=5*nl+1
      imin=loop-6
      d8=dd/h
      do 80 i=imin,loop
      k=i+1-imin
      x=(r(i)*w+2.*z(i))/con
      gp(k)=-rpor(i)*(kap*g(i)+(con*r(i)+x)*f(i)+r(i)*rf(i))
   80 fp(k)=rpor(i)*(kap*f(i)+x*g(i)+r(i)*rg(i))
      u=g(loop)
      v=f(loop)
      min=loop+1
      do 81 i=min,max
      all=rpor(i)*kap
      akk=-all
      x=(r(i)*w+2.*z(i))/con
      amm=-rpor(i)*(con*r(i)+x)
      ann=rpor(i)*x
      ax=-rp(i)*rf(i)
      ay=rp(i)*rg(i)
      brk=1.0-a(1)*akk/d8
      brl=1.0-a(1)*all/d8
      brn=a(1)*ann/d8
      brm=a(1)*amm/d8
      del=brk*brl-brm*brn
      rk=a(1)*ay+a(2)*fp(7)+a(3)*fp(6)+a(4)*fp(5)+
     1a(5)*fp(4)+a(6)*fp(3)+a(7)*fp(2)+a(8)*fp(1)
      sk=a(1)*ax+a(2)*gp(7)+a(3)*gp(6)+a(4)*gp(5)+
     1a(5)*gp(4)+a(6)*gp(3)+a(7)*gp(2)+a(8)*gp(1)
      rk=v+rk/d8
      sk=u+sk/d8
      u=(brl*sk+brm*rk)/del
      v=(brk*rk+brn*sk)/del
      g(i)=u
      f(i)=v
      do 85 j=1,6
      gp(j)=gp(j+1)
   85 fp(j)=fp(j+1)
      gp(7)=akk*g(i)+amm*f(i)+ax
   81 fp(7)=all*f(i)+ann*g(i)+ay
      return
      end

      subroutine qed(vpint,sfint)
      include "global.par"
      implicit doubleprecision(a-h,o-z)
c
c  this routine evaluates the qed corrections to the energy levels
c  due to vacuum polarisation (correct to first order) and a crude
c  approximation to the self energy.
c  the v.p. contribution is calculated using the results of fullerton
c  and rinker phys. rev. a  vol 13  p 1283 (1976) while the s.e.
c  contribution is estimated for s, p- and p orbitals by interpolating
c  among the values given by p. mohr for coulomb type wavefunctions
c  after an effective nuclear charge, zeff, is obtained by finding the
c  zeff required to give a coulomb orbital with the same average r
c  as the mcdf orbital.
c
c
      character*4 ident,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      common gh(NGP,NS),fh(NGP,NS),gfh(16,NS),g(NGP,NS),f(NGP,NS)
      dimension vac(NGP),u(NGP)
      dimension vpint(NS),sfint(NS)
c
      clocal=137.0359895d0
      pi=dacos(-1.0d0)

*  calculates the vacuum polarisation potential at each of the
*  grid points.

      zz = jz
      call vacpol(zz,vac)

****  obtain contribution from each orbital
      write(6,1010)
 1010 format(/12x,' <r>',10x,' Zeff')

      do 300 j = 1,jmax

****  obtain v.p. contribution for orbital j

         mj = ms(j)
         do 100 i=1,mj
              u(i)=(g(i,j)**2+f(i,j)**2)*vac(i)*rp(i)
 100     continue
         vpint(j) = rint(u,1,mj,7,h)

****  obtain s.e. contribution for orbital i.

         if (kap(j).eq.-1.or.kap(j).eq.1.or.kap(j).eq.-2) then

****  find average r for mcdf orbital
            

            do 200 i = 1,mj
               u(i) = (g(i,j)**2+f(i,j)**2)*r(i)*rp(i)
 200        continue
            rav=rint(u,1,mj,7,h)

****  find effective z of coulomb orbital with same average r

            call zefr(j,rav,zeff)
            write(6,1000) lab(j),rav,zeff
 1000       format(4x,a,1p,e12.3,0p,f12.3)

****  interpolate among p. mrs data

            call fzalf(zeff,j,valu)
****  scale as required

            valu=( zeff**4/clocal**3)*valu/(pi*( dfloat(n(j))**3))
         else
            valu = 0.0
         endif
         sfint(j) = valu
 300  continue
      return
      end

      subroutine vacpol(z,vac)
      include "global.par"
      implicit doubleprecision(a-h,o-z)
c
c  this routine sets up the vacuum polarization potential for
c  a point charge z at each grid point using the analytic
c  functions defined by l. wayne fullerton and g. a. rinker jr.
c  in phys. rev. a   vol 13, page 1283, (1976).
c
c  no subroutines called.
c
      common/radial/r(NGP),rp(NGP),rpor(NGP),h,max
      dimension vac(NGP)
c
c  the following are the analytic functions needed:
c
      p(x) = -0.71740181754d 00 + x*( 1.1780972274d 00 +
     x  x*( -0.37499963087d 00 + x*( 0.13089675530d 00 +
     x  x*( -0.038258286439d 00 + x*( -0.0000242972873d 00 +
     x  x*( -0.0003592014867d 00 + x* (-1.71700907d-05)))))))
c
      b(x) = -64.0514843293d 00 + x*( 0.711722714285d 00 + x)
c
      cf(x) = 64.0514843287d 00 + x*( -0.711722686403d 00 +
     x  x* 8.042207748d-04)
c
      d(x) = 217.2386409d 00 + x*( 1643.364528d 00 +
     x  x*( 2122.244512d 00 + x*( -45.12004044d 00 + x )))
c
      e(x) = 115.5589983d 00 + x*( 1292.191441d 00 +
     x  x*( 3831.198012d 00 + x* 2904.410075d 00 ))
c
      pi=dacos(-1.0d0)
      clocal=137.0359895d0
      factor = -z/(1.5*pi*clocal)
      do 100 i=2,max
         x = 2.0*r(i)*clocal
         if (x .le. 1.0) then
            y=x*x
            vac(i) = factor*( p(x) + dlog(x)*b(y)/cf(y))/r(i)
         elseif (x.gt. 163.0) then
            vac(i)=0.0
         else
            y=1.0/x
            vac(i)=factor*dexp(-x)*d(y) / e(y) /(r(i)*x**(1.5d0))
         endif
 100  continue
      vac(1) = 0.0
      return
      end

      subroutine fzalf(zeff,j,slfe)
      include "global.par"
      parameter(NUMVAL=11)
      implicit doubleprecision(a-h,o-z)
      character*4 ident,lab
c
c  this routine obtains an estimate of the self energy contribution
c  to the energy resulting from either 1s, 2s, 2p- or 2p orbital in
c  the field of a point nucleus with effective charge zeff.
c  the values are interpolated among the values supplied by p.j. mohr
c
c  subroutine called: interp
c
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)

      dimension val1s(NUMVAL),val2s(NUMVAL),val2p1(NUMVAL),
     & val2p3(NUMVAL),arg(NUMVAL)

c  1s data :
      data val1s/ 4.654d0,  3.246d0,  2.5519d0,  2.1351d0,
     x            1.8644d0, 1.6838d0, 1.5675d0,  1.5032d0,
     x            1.4880d0, 1.5317d0, 1.6614d0/
c  2s data :
      data val2s/ 4.8930d0, 3.5063d0, 2.8391d0,  2.4550d0,
     x            2.2244d0, 2.0948d0, 2.0435d0,  2.0650d0,
     x            2.1690d0, 2.3870d0, 2.7980d0/
c  2p- data :
      data val2p1/ -0.1145d0,  -0.0922d0,  -0.0641d0,  -0.0308d0,
     x              0.0082d0,   0.0549d0,   0.1129d0,   0.1884d0,
     x              0.2934d0,   0.4530d0,   0.7250d0/
c  2p data
      data val2p3/ 0.1303d0, 0.1436d0, 0.1604d0, 0.1794d0,
     x             0.1999d0, 0.2215d0, 0.2440d0, 0.2671d0,
     x             0.2906d0, 0.3141d0, 0.3367d0/
c  z data values :
      data arg/ 10.0d0,  20.0d0, 30.0d0,  40.0d0,  50.0d0,
     x          60.0d0,  70.0d0, 80.0d0,  90.0d0, 100.0d0,
     x         110.0d0/
      
      if (n(j).eq.1) then
c     ...1s case...
         call interp(arg,val1s,zeff,numval,slfe)
      return
      endif
      slfe = 0.0
c     ....ns case...
      if (kap(j).eq.-1) call interp(arg,val2s,zeff,numval,slfe)
c     ....np- case....
      if (kap(j).eq.1) call interp(arg,val2p1,zeff,numval,slfe)
c     ....np case....
      if (kap(j).eq.-2) call interp(arg,val2p3,zeff,numval,slfe)
      return
      end

      subroutine interp(arg,val,x,n,y)
      implicit doubleprecision(a-h,o-z)
c
c  uses lagrange interpolation formula to obtain value of val(x).
c  arg(i),val(i) ,i=1,n contain the data values to interpolate
c  among.
c
c  no subroutines called.
c
      dimension arg(n),val(n)
c
      y=0.0
      do 3 l=1,n
      pl=1.0
      do 2 j=1,n
      if (l-j)1,2,1
    1 pl = (x-arg(j)) * pl / (arg(l)-arg(j))
    2 continue
      y=y+pl*val(l)
    3 continue
      return
      end

      subroutine zefr(j,r,ze)
      include "global.par"
      implicit doubleprecision(a-h,o-z)
c
      dimension a(3),b(3)
c
c   this subroutine calculates an effective charge ze such
c   that in a coulomb field of charge ze the mean radius
c   <r> of orbital j is equal to the input value r .
c      find ze such that   <r>-r=0
c   the newton method is used with accuracy acu and maximum
c   number of iterations itr .
c
c   no subroutines called
c
      character*4 ident,lab
      common/ell/ident,lab(NS)
      common/elm/jmax,jz,nat,nuc,ion,io,in,inf,
     1n(NS),kap(NS),k(NS),l(NS),jdg(NS),ms(NS),nct(NS),kc(NS),
     2iof(NS),wh(NS),wf(NS),dw(NS),gs(NS),fs(NS)
      data c /137.0359895d0/

      nn=n(j)
      kk=kap(j)
      ll=l(j)
      eps=1.0d-05
      acu=1.0d-06
      itr=30
      mk=iabs(kk)
      nns=nn*nn
c
c   initial estimate using non-relativistic formula
c
      a(3)=dfloat(nns+nns+nns-ll*(ll+1))/(r+r)
      if(a(3).gt.dfloat(jz)) a(3) = dfloat(jz)
      ze = a(3)
      fk=dfloat(kk)
      fks=fk*fk
      fnr=dfloat(nn-mk)
c
c   begin iterations on the relativistic formula
c
      i=1
    1 a(1)=a(3)-eps
      a(2)=a(3)+eps
      do 2 ii=1,3
      wa=a(ii)/c
      wa=wa*wa
      wb=dsqrt(fks-wa)+fnr
      wa=wb*wb+wa
      bn=dsqrt(wa)
      wa=wb*(wa+wa+wa-fks)-bn*fk
      wb=a(ii)*bn
    2 b(ii)=wa/(wb+wb)-r
      wb=a(3)
      wa=eps*b(3)/(b(2)-b(1))
      a(3)=a(3)-(wa+wa)
      if (dabs(a(3)-wb)-acu) 5,5,3
    3 i=i+1
      if (i-itr) 1,1,4
    4 write (6,300) lab(j)
      return
c
    5 ze=a(3)
      return
c
  300 format (' iteration limit exceeded in zefr for ',a)
      end

c     ========================================
      include "clrx.f"
      include "cre.f"
      include "enclab.f"
      include "rint.f"
      include "wfun.f"
      include "wint.f"
      include "yfun1.f"
      include "yint1.f"
      include "zfun.f"
      include "zint.f"