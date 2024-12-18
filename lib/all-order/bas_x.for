C     ======================================== 09/01/07
      program Bas_x           !### last update 13/09/09
      include "conf.par"
      include "hfd.par"
      include "basx.par"
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c       Bas_x reads basis set from HFD.DAT and rewrites
c       it in a form consistent with WJ and MS codes
c       That requires:
c       1. Change of the grid so that R(1)=0 and,
c          therefore, no Taylor expansion is necessary
C       2. Orbitals are stored in partial wave order.
c       Output files:
c           hfspl.1  - radial grid
c           hfspl.2  - orbitals and energies
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsp/Nsp/Nso/Nso/Nsv/Nsv/Ierr/Ierr
     1        /MaxT/MaxT/Ii/Ii/Nc/Nc/Ne/Ne/Vfst/n4,j4,ch4
     3        /FNAME/FNAME/Kkin/Kkin/Kbas/Kbas(IPsp)
       common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /Qnl/Qnl(IPsp)/Qnl1/Qnl1(IPsp)
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /P/P(IP6)/Q/Q(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
       common /kout/kout/small/small/let/let(6)
       common /ipmr/ipmr
       character*1 let,ch4,FNAME*12
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        MaxT=9           !### length of expansion at the origin
        kout=1           !### output details
        call recunit
        irec=2*IP6*IPmr  !### record length in DAT files
        let(1)='s'
        let(2)='p'
        let(3)='d'
        let(4)='f'
        let(5)='g'
        let(6)='h'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c        open(unit=11,file='BAS_X.RES',status='UNKNOWN')
        call Init  ! reads the grid and QN of orbitals from HFD.DAT
        call SortPW
        call NewGrid
        call PWorbitals
       stop
      end
C     =================================================
      include "readf.inc"
      include "rintms.f"
      include "rec_unit.inc"
C     =================================================
      subroutine Init
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "phys.par"
      include "hfd.par"
      include "basx.par"
       common /Ns/Ns/II/II/Kt/Kt/Ng/Ng/Rmax/Rmax
     1        /Z/Z/Cl/Cl/H/H/Rnuc/Rnuc
       common /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /R/R(IP6)/V/V(IP6)
       common /ipmr/ipmr
       dimension P(IP6),Q(IP6),P1(IP6),Q1(IP6),PQ(4*IP6)
       logical longbasis
       dimension IQN(4*IPs),Qq1(IPs)
       equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
       Equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),
     >        (P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c small number:
        c1=0.01d0
c speed of light is taken from "phys.par":
        Cl=DPcl
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(12,file='HFD.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*IPmr,err=700)
        call ReadF (12,1,P,Q,2)
        call ReadF (12,2,R,V,2)
        call ReadF (12,3,P1,Q1,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Z   =PQ(1)
        Ns  =PQ(2)+C1
        II  =PQ(3)+C1
        R1  =PQ(4)
        Rmax=dabs(PQ(5))
        H   =PQ(6)
        Bt  =PQ(7)
        Al  =PQ(8)
        Kt  =PQ(9)+C1
        Ng  =PQ(10)+C1
        Rnuc=PQ(13)
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
        write( *,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Ns
c        write(11,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Ns
 5      format (/4X,'Z   = ',F6.2,5X,'Kt  =',I3,  7X,'II =',I4,
     >          /4X,'H   =',F7.4, 5X,'R2  =',F6.2,4X,'R1 =',E11.4,
     >          /4X,'Rnuc=',E11.4,1X,'Al  =',F7.4,3X,'Bt =',F5.2,
     >          /4X,'Ns  =',I5)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
c          write(11,*) ' Using variant for long basis '
          do ni=1,Ns
            Nn(ni)=IQN(4*ni-3)
            Ll(ni)=IQN(4*ni-2)
            Kk(ni)=IQN(4*ni-1)
            Jj(ni)=IQN(4*ni)
            Qq(ni)=Qq1(ni)
          end do
        else
          if=20
          do ni=1,Ns
            if=if+1
            Nn(ni)=PQ(if)+c1
            if=if+1
            Ll(ni)=PQ(if)+c1
            if=if+1
            Qq(ni)=PQ(if)
            if=if+2
            c2=dsign(c1,PQ(if))
            Kk(ni)=PQ(if)+c2
            if=if+1
            c2=dsign(c1,PQ(if))
            Jj(ni)=PQ(if)+c2
          end do
        end if
        close(12)
       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
        write(11,75)
75      format(/2X,'file HFD.DAT is absent'/)
       stop
      end
C     =================================================
      subroutine SortPW
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/lmax/lmax/npwx/npwx
     >        /Npw/Npw(IPxo)/Norb/Norb(IPxo,IPxo)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       write(*,*) ' Sort partial waves from HFD.DAT'
       imax=0
       npwx=0
       lmax=0
       do i=1,IPxo
         Npw(i)=0
         do k=1,IPxo
           Norb(i,k)=0
         end do
       end do

       do n=1,Ns
         ln=Ll(n)
         jn=Jj(n)
         kn=Kk(n)
         ind=ln+iabs(kn)
         if (ind.GT.imax) imax=ind
         if (ln.GT.lmax) lmax=ln
         Npw(ind)=Npw(ind)+1
         if (Npw(ind).GT.npwx) npwx=Npw(ind)
         kpw=Npw(ind)
         Norb(ind,kpw)=n
       end do
       write(*,*) imax,' partial waves found; lmax=',lmax
       if (imax.GT.IPxo) then
         write(*,*) ' IPxo=',IPxo,' is too small!'
         read(*,*)
         stop
       end if
       if (imax.NE.2*lmax+1) then
         write(*,*) ' Error: PW number differs from 2*lmax+1'
         read(*,*)
         stop
       end if
       do i=1,imax
         nx=Npw(i)
         write(*,*) ' PW ',i,' includes ',nx,' orbitals'
c         do n=1,nx
c           no=Norb(i,n)
c           write(*,5) no, Nn(no),LL(no),Jj(no)
c 5         format(4X,I3,' n=',I2,' l=',I1,' j='I2,'/2')
c         end do
c         read(*,*)
       end do
c       read(*,*)
       Return
      end
C     =================================================
      subroutine NewGrid
      implicit real*8 (a-h,o-z)
      include "basx.par"
      include "hfd.par"
       common /ii/ii/iix/iix/lmax/lmax/Kt/Kt/npwx/npwx
     >        /H/H/Rnuc/Rnuc/Rmax/Rmax
       common /R/R(IP6)/V/V(IP6)/Rn/Rn(IPx6)/Vn/Vn(IPx6)
     >        /Rpor/Rpor(IPx6)/Npw/Npw(IPxo)
       dimension tx(IPxo+IPkx)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=R(1)
        r2=R(2)
        v1=V(1)*h
        xk0=r1/v1
        k=2*xk0
        iter=1
 1      write(*,*) ' k0=',xk0,' k=',k,' ii=',ii

        if (ii+k.GT.IPx6) then
          maxk=IPx6-ii
          write(*,*) ' for IPx6=',IPx6,' max k =',maxk
          if (maxk.GT.xk0.AND.iter.EQ.1) then
            k=maxk
            write(*,*) ' trying to use k=',k
          end if
        end if

        if (iter.GT.1) then
          write(*,*) ' choose new k: '
          read(*,*) k
        end if

        rox=k*h
        cx=(v1*k-r1)/rox**2
        write(*,*) ' cx=',cx
        kx=k+1
        ierr=0
        do i=1,kx+1
          ix=kx-i
          Rn(i)=r1-v1*ix+cx*(h*ix)**2
          Vn(i)=V(1)-2*cx*h*ix
          if (i.LT.kx) then
            write(*,5) i,Rn(i),Vn(i)
 5          format(4X,I3,' R   =',E12.5,' V   =',E12.5)
          else
            drn=(Rn(i)-R(i-k))/R(i-k)
            dvn=(Vn(i)-V(i-k))/V(i-k)
            if (i.EQ.kx.AND.(dabs(drn)+dabs(0.03d0*dvn)).GT.1.d-6)
     >        ierr=ierr+1
            if (i.EQ.kx+1.AND.(dabs(drn)+dabs(0.03d0*dvn)).GT.5.d-3)
     >        ierr=ierr+1
            write(*,15) i,Rn(i),Vn(i),R(i-k),V(i-k)
 15         format(4X,I3,' Rnew=',E12.5,' Vnew=',E12.5
     >             /7X,' Rold=',E12.5,' Rold=',E12.5)
          end if
        end do
        if (ierr.GT.0) then
          write(*,*) ' Not a perfect match at point ',kx
          write(*,*) ' Is it OK? [0,1] '
          read(*,*) iyes
        else
          iyes=1
        end if
        iter=iter+1
        if (iyes.NE.1) goto 1
        do i=1,ii
         ix=i+kx-1
         Rn(ix)=R(i)
         Vn(ix)=V(i)
        end do
        iix=ii+k
        do i=1,iix
          Rpor(i)=Vn(i)/Rn(i)
        end do
        write(*,*) ' New grid formed, iix=',iix
        do i=1,IPxo+IPkx
          tx(i)=Rmax
        end do
        write(*,*) ' tx =',Rmax,' npwx=',npwx
        write(*,*) ' Forming file hfspl.1...'
        open(12,file='hfspl.1',form='unformatted')
        write(12) npwx,kt,tx,lmax
        write(12) 0
        write(12) Rn,Vn,Rpor,h,iix
        write(12) Npw
        close(12)
        write(*,*) ' File hfspl.1 formed'
c        read(*,*)
       Return
      end
C     =================================================
      subroutine PWorbitals
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/lmax/lmax/MaxT/MaxT/npwx/npwx
     >        /Npw/Npw(IPxo)/Norb/Norb(IPxo,IPxo)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
       common /ii/ii/iix/iix/Kt/Kt/H/H/Rnuc/Rnuc
       common /P/P(IP6)/Q/Q(IP6)/Pn/Pn(IPx6)/Qn/Qn(IPx6)
     >        /En/En(2*IPxo)/Rn/Rn(IPx6)/Vn/Vn(IPx6)/C/C(IPx6)
       common /ipmr/ipmr
       dimension g(IPx6,IPxo),f(IPx6,IPxo)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       write(*,*) ' Rewriting orbitals on new grid'
       open(12,file='HFD.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*IPmr,err=700)
       open(13,file='hfspl.2',access='DIRECT',
     >       status='UNKNOWN',recl=4*4*IPxo*(IPx6+1)) ! note that 4 stands for IPmr, so
                                                      ! record may be longer than required
       tlr=1.d-6
       write(*,*) ' Orthonormality tolerance:',tlr
       imax=2*lmax+1
       kx=iix-ii
       rx1=Rn(kx+1)
       if (dabs(rx1-rnuc).GT.1.d-10) then
         write(*,*) ' Error: rx1=',rx1,' .NE. rnuc=',rnuc
         read(*,*)
         stop
       end if
       do ind=1,imax
         nx=Npw(ind)
         write(*,*) ' PW ',ind,' nx=',nx
         ibad=0
         do n=1,nx
           nold=Norb(ind,n)
           call ReadF (12,nold+4,P,Q,2)
           gam=P(ii+4)

           do j=1,kx+1
             r=Rn(j)
             rg=r**gam
             pp=0.d0
             qq=0.d0
             rk=1.d0
             do k=0,MaxT
               ik=ii+5+k
               pp=pp+P(ik)*rk
               qq=qq+Q(ik)*rk
               rk=rk*r/Rnuc
             end do
             Pn(j)=pp*rg
             Qn(j)=qq*rg
           end do
           pnt=Pn(kx+1)  ! Taylor value for rnuc
           qnt=Qn(kx+1)  !

           do j=1,ii
             Pn(j+kx)=P(j)
             Qn(j+kx)=Q(j)
           end do
           perr=dabs(pnt/Pn(kx+1)-1.d0)
           qerr=dabs(qnt/Qn(kx+1)-1.d0)
           if(perr+qerr.GT.1.d-6) then
             write(*,5) n,nold,perr,qerr
 5           format(4X,'Error for n=',I2,' orb ',I3,':',2E15.5,
     >             /4X,'-Taylor expansion does not match the function')
             read(*,*)
           end if

           En(npwx+n)=-P(ii+1)

           do j=1,iix
             g(j,n)=Pn(j)
             f(j,n)=Qn(j)
             C(j)=(Pn(j)**2+Qn(j)**2)*Vn(j)
           end do
           sn=rintms(c,1,iix,7,h)
           write(*,15) n,nold,sn,En(npwx+n)
 15        format(4X,I2,' n_old=',I3,' norma=',F10.6,' E=',F15.8)
           if (dabs(sn-1.d0).GT.tlr) then
             ibad=ibad+1
           end if
         end do
         if (ibad.GT.0) then
           write(*,*) ibad,' orbitals with bad norma from ',nx
           write(*,*) ' Start orthogonalization: '
           call ort(g,f,nx,tlr)
         end if
         write(13,rec=ind) g,f,En
         write(*,*) ' PW ',ind,' writen to file; npwx=',npwx
c         read(*,*)
       end do
       close(13)
       write(*,*) ' File hfspl.2 formed.'
c       read(*,*)
       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
        write(11,75)
75      format(/2X,'file HFD.DAT is absent'/)
       stop
      end
c     =====================================================
      subroutine ort(g,f,nx,tlr)
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /iix/iix/H/H
       common /C/C(IPx6)/Vn/Vn(IPx6)
       dimension g(IPx6,IPxo),f(IPx6,IPxo)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       iter=1
 1     err=0.d0
       do n=1,nx
         s=0.d0
         scal=0.d0
         do m=1,n
           do j=1,iix
             C(j)=(g(j,n)*g(j,m)+f(j,n)*f(j,m))*Vn(j)
           end do
           scal=rintms(c,1,iix,7,h)
c           write(*,*) ' n=',n,' m=',m,' scal=',scal
           if (m.LT.n) then
             err0=0.d0
             s=s+scal**2
             do j=1,iix
               g(j,n)=g(j,n)-scal*g(j,m)
               f(j,n)=f(j,n)-scal*f(j,m)
             end do
           else
             err0=1.d0
             xnorm=1.d0/dsqrt(scal)
             do j=1,iix
               g(j,n)=g(j,n)*xnorm
               f(j,n)=f(j,n)*xnorm
             end do
           end if
           err1=dabs(scal-err0)
           if (err.LT.err1) err=err1
         end do
c         read(*,*)
       end do
       write(*,*) ' iter=',iter,' err=',err
       iter=iter+1
       if (iter.LE.3.AND.err.GT.tlr) goto 1
c       read(*,*)
       return
      end
