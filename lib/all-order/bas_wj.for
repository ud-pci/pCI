C     ======================================== 04/03/09
      program Bas_WJ          !### last update 02/10/09
      include "conf.par"
      include "hfd.par"
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c       Bas_WJ reads basis set from files hfspl.1 and hfspl.2
c       and rewrites it to file FNAME1 in HFD.DAT format
c       02/10/09: also forms file FNAME2=BASS.INP
c
c       INPUT: bas_wj.in (=atom.in - console input used by tdhf code)
c              spl.in - the file nspl code uses
c              hfspl.1 & hfspl.2 - the output of nspl code
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsp/Nsp/Nso/Nso/Nsv/Nsv/Ierr/Ierr
     1        /MaxT/MaxT/Ii/Ii/Nc/Nc/Ne/Ne/Vfst/n4,j4,ch4
     3        /FNAME1/FNAME1/FNAME2/FNAME2/Kkin/Kkin/Kbas/Kbas(IPsp)
       common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /Qnl/Qnl(IPsp)/Qnl1/Qnl1(IPsp)
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /P/P(IP6)/Q/Q(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
       common /kout/kout/small/small/let/let(6)
       character*1 let,ch4,FNAME1*12,FNAME2*12
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        MaxT=9           !### length of expansion at the origin
        kout=1           !### output details
        let(1)='s'
        let(2)='p'
        let(3)='d'
        let(4)='f'
        let(5)='g'
        let(6)='h'
        FNAME1='WJ.DAT'
        FNAME2='BASS.INP'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        call Input ! reads data from bas_wj.in & spl.in
        call Init  ! reads the grid from hfspl.1 & orbitals from hfspl.2
        call SortPW  ! select & reorder orbitals
        read(*,*)
        call DefRmax ! defines Rmax for new grid
        call NewGrid ! formes new grid
        write(*,*) ' Push...'
        read(*,*)
        call NewBasis
        write(*,*) ' Want to form BASS.INP file (0,1)?'
        read(*,*) iyes
        if (iyes.EQ.1) call BassInp
       stop
      end
C     =================================================
      include "readf.inc"
      include "rintms.f"
C     =================================================
      subroutine Input
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "phys.par"
      include "hfd.par"
       common /name/name/Ns/Ns/Nso/Nso/II/II/Rmax/Rmax/Lmax/Lmax
     1        /Z/Z/A/A/Rnuc/Rnuc/cnuc/cnuc/tnuc/tnuc
       common /Kk/Kk(IPs)/Nn/Nn(IPs)
       Character name*4
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=12,file='bas_wj.in',status='OLD',err=700)
        read(12,*) name, Ns, iz, ia, id, id, id, Nv, ival
        Z=iz
        A=ia
        if (ival.EQ.0) then
          Nso=Ns                ! number of core shells
        else
          Nso=Nv-1
        end if
        do i=1,Ns
          read(12,*) Nn(i),Kk(i)
        end do

        read(12,*)
        read(12,*)
        read(12,*)

        read(12,*) dummy,cnuc,tnuc
        close(12)


        open(unit=12,file='spl.in',status='OLD',err=710)

        read(12,*) Lmax
        read(12,*) Rmax
        read(12,*)
        read(12,*) r1,r2,Ii
        write(*,5) name,Z,A,Nso,Ns,Rmax,Ii,cnuc,tnuc,
     >             (i,NN(i),KK(i),i=1,Ns)
 5      format(//'Basis set conversion for ',a4,' Z=',F5.1,' A=',F5.1,
     >        /'Nso=',I3,' Ns=',I4,
     >        /'Rmax=',F5.1,' II=',I4,
     >        /'c_nuc=',F8.4,' t_nuc=',F8.4,
     >        /'Shell, n, kappa',/(3x,I2,I3,2x,I4))
        close(12)

       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
75      format(/2X,'file bas_wj.in is absent'/)
       stop
710     write( *,85)
85      format(/2X,'file spl.in is absent'/)
       stop
      end
C     =================================================
      subroutine Init
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /II/II/Rmax/Rmax/Lmax/Lmax
     >        /iix/iix/npwx/npwx/Kt/Kt/H/H
       common /Rn/Rn(IPx6)/Vn/Vn(IPx6)/Rpor/Rpor(IPx6)
     >        /g0/g0(IPx6,IPxo,IPpw)/f0/f0(IPx6,IPxo,IPpw)
     >        /En0/En0(2*IPxo,IPpw)
       dimension tx(IPxo+IPkx)
C     - - - - - - - - - - - - - - - - - - - - - - - - -

        write(*,*) ' Reading file hfspl.1...'
        open(12,file='hfspl.1',form='unformatted',status='OLD',err=700)

        read(12,err=710,end=710) npwx,kt,tx,lmax1

        kt=1  ! I do not know what kt means in hfspl.1!

        read(12,err=710,end=710)

        read(12,err=710,end=710) Rn,Vn,Rpor,h,iix

        close(12)

        ierr=0
        if (lmax1.NE.Lmax) then
          write(*,*) 'In hfspl.1 Lmax=',lmax1,
     >               ', while spl.in gives ',Lmax
          Lmax=min(Lmax,lmax1)
          write(*,*) 'We will use Lmax=',Lmax
          ierr=ierr+1
        end if

        if (iix.NE.Ii) then
          write(*,*) 'In hfspl.1 Iix=',Iix,
     >               ', while spl.in gives Ii=',Ii
          Ii=min(iix,Ii)
          write(*,*) 'We will use Ii=',Ii
          ierr=ierr+1
        end if

        iix=0
        do i=1,Ii
          if (Rn(i).LE.Rmax) iix=i
        end do

        write(*,5) h,iix,Rn(iix),npwx
 5      format(/'File hfspl.1 read.',/' h=',F8.4,' Rn(',I4,')=',F10.2,
     >         ' npwx=',I2)

        write(*,15) Lmax
 15     format(/' Lmax=',i2,' give Lmax for HFD basis set: '$)
        read(*,*) Lmax

        write(*,25)
 25     format(/' Reading file hfspl.2...')


        open(13,file='hfspl.2',access='DIRECT',            !### note, that 4 stands instead
     >       status='OLD',recl=4*4*IPxo*(IPx6+1),err=720)  !### of IPmr. Thus, record may be
                                                           !### longer, than necessary!
        if (IPpw.LT.2*Lmax+1) then
          write(*,*) ' IPpw=',IPpw,' is too small for Lmax=',Lmax
          stop
        end if

        do k=1,2*Lmax+1
          read(13,rec=k,err=730)
     >                   ((g0(m,n,k),m=1,IPx6),n=1,IPxo),
     >                   ((f0(m,n,k),m=1,IPx6),n=1,IPxo),
     >                   (En0(m,k),m=1,2*IPxo)
        end do
        close(13)
        write(*,*) ' File read OK'

       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
75      format(/2X,'file hfspl.1 is absent'/)
       stop
710     write( *,85)
85      format(/2X,'unexpected end of file hfspl.1'/)
       stop
720     write( *,95)
95      format(/2X,'file hfspl.2 is absent'/)
       stop
730     write( *,105)
105     format(/2X,'unexpected end of file hfspl.2'/)
       stop
      end
C     =================================================
      subroutine SortPW
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/lmax/lmax/npwx/npwx
     >        /Npw/Npw(IPpw)/Norb/Norb(IPxo,IPpw)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       write(*,*)
       write(*,*) ' Sort partial waves from hfspl.2'

       write(*,*) ' npwx=',npwx,'. Give new n_max for each l:'
       npwx1=0

       do ind=1,2*lmax+1
         li=ind/2             ! = l
         ki=ind-li
         if (ki.GT.li) ki=-ki ! = kappa
         if (ki.GT.0.OR.li.EQ.0) then
           write(*,5) li
 5         format(' n_max for l=',i2,': ',$)
           read(*,*) n_max
           Npw(ind)=n_max-li
           if (Npw(ind).GT.npwx1) npwx1=Npw(ind)
         else
           Npw(ind)=Npw(ind-1)
         end if
       end do

       if (npwx1.GT.IPxo) then
         write(*,*) 'IPxo=',IPxo,' is too small for npwx=',npwx1
         stop
       end if

       nhfd=0

       do ind=1,2*lmax+1      ! finding orbitals from the list
                              ! in bas_wj.in
         li=ind/2             ! = l
         ki=ind-li
         if (ki.GT.li) ki=-ki ! = kappa
c         write(*,*) 'ind=',ind,' l=',li,' kappa=',ki,' Npw=',Npw(ind)

         mx=Npw(ind)
         do m=1,mx
           nm=li+m
           i1=0
           do i=1,Ns
             if(Nn(i).EQ.nm.AND.Kk(i).EQ.ki) then
               i1=i
               goto 100
             end if
           end do
 100       Norb(m,ind)=i1
           if (i1.NE.0) then
             nhfd=nhfd+1
             Ll(i1)=li
             Jj(i1)=2*iabs(ki)-1
             write(*,*) 'orb ',i1,' n=',nm,' kappa=',ki,' jj=',Jj(i1)
           end if
         end do
       end do

       write(*,*) nhfd,' orbitals found, expected ',Ns
       if (Ns.NE.nhfd) stop
       write(*,*) ' Adding more valence orbitals in shell order'

       do n=1,npwx1+lmax   ! principle QN
         indx=min(2*n-1,2*lmax+1)
         i1=0
         do ind=1,indx   ! PW for given n
           li=ind/2             ! = l
           ki=ind-li
           if (ki.GT.li) ki=-ki ! = kappa
           m=n-li
           if (m.LE.0.OR.m.GT.Npw(ind)) goto 200 ! orbital is out of range
           if (Norb(m,ind).NE.0)        goto 200 ! orbital is already included
           Ns=Ns+1
           if (Ns.GT.IPs) then
             write (*,*) ' Error: number of orbitals exceeds IPs=',IPs
             read (*,*)
             stop
           end if
           i1=i1+1
           Norb(m,ind)=Ns
           Nn(Ns)=n
           Ll(Ns)=li
           Kk(Ns)=ki
           Jj(Ns)=2*iabs(ki)-1
           write(*,*) 'orb ',Ns,' n=',n,' kappa=',ki,' jj=',Jj(Ns)
 200     end do
c         if (i1.GT.0) read(*,*)
       end do
       nsmax=(4*IP6-21)/3
       write(*,15) Ns,nsmax,IPs
 15    format(i5,' orbitals found',
     >        /' (should be less than (4*IP6-21)/3=',
     >        i3,' and less than IPs=',i3,')')
c       read(*,*)
       if (nsmax.LT.Ns.OR.IPs.LT.Ns) stop
       Return
      end
C     =================================================
      subroutine DefRmax
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/iix/iix/ii/ii/Rmax/Rmax/h/h
     >        /Kk/Kk(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
       common /g0/g0(IPx6,IPxo,IPpw)/f0/f0(IPx6,IPxo,IPpw)
     >        /Rn/Rn(IPx6)/Vn/Vn(IPx6)/Pn/Pn(IPx6)/Qn/Qn(IPx6)
     >        /Norb/Norb(IPxo,IPpw)
       dimension C(IPx6),P(IPx6),Q(IPx6)
       logical norm,      ! orbitals are normalized if norm=TRUE
     >         ortho,     ! orbitals are orthogonalized if ortho=TRUE
     >         bad
C     - - - - - - - - - - - - - - - - - - - - - - - - -

        norm =.TRUE.        ! change this to switch normalization on/off
        ortho=.TRUE.        ! change this to switch orthogonalization on/off
        norm=norm.OR.ortho  ! orthogonalization followed by normalization

        i4=0
        i6=0
        i8=0
        ibad=0
        iort=0
        tlr=1.d-3

        do m=1,Ns
          nm=Nn(m)
          lm=Ll(m)
          km=Kk(m)
          ind=lm+iabs(km)
          npw=nm-lm

          if (Norb(npw,ind).NE.m) then
            write(*,*) ' Wrong QN for orbital ',m
            write(*,*) ' n=',nm,' l=',lm,' kappa=',km
            stop
          end if

          do i=1,ii
            Pn(i)=g0(i,npw,ind)
            Qn(i)=f0(i,npw,ind)
          end do

          if (ortho) then
            do j=1,npw-1
              do i=1,ii
                C(i)=(Pn(i)*g0(i,j,ind)
     >               +Qn(i)*f0(i,j,ind))*Vn(i)
              end do
              sn=rintms(C,1,iix,7,h)

              if (dabs(sn).GT.tlr) then
                iort=iort+1
                write(*,5) nm,km,j+lm,km,sn
 5              format(4X,'<',2I2,'|',2I2,'> =',F10.6)
              end if

              do i=1,ii
                Pn(i)=Pn(i)-sn*g0(i,j,ind)
                Qn(i)=Qn(i)-sn*f0(i,j,ind)
              end do
            end do
c            if (iort.GT.0) read(*,*)
          end if

          do i=1,ii
            s=Pn(i)**2+Qn(i)**2
            C(i)=s*Vn(i)
            if (s.GT.1.d-8.AND.i.GT.i8) i8=i
            if (s.GT.1.d-6.AND.i.GT.i6) i6=i
            if (s.GT.1.d-4.AND.i.GT.i4) i4=i
          end do
          sn=rintms(C,1,iix,7,h)
          bad=dabs(sn-1.d0).GT.tlr
          if (norm) then
            sn1=1.d0/dsqrt(sn)
            do i=1,ii
              Pn(i)=Pn(i)*sn1
              Qn(i)=Qn(i)*sn1
              s=Pn(i)**2+Qn(i)**2
              C(i)=s*Vn(i)
              g0(i,npw,ind)=Pn(i)
              f0(i,npw,ind)=Qn(i)
            end do
            sn2=rintms(C,1,iix,7,h)
          else
            sn2=sn
          end if
          if (bad) then
            write(*,15) m,nm,km,sn,sn2
 15          format(2X,I4,' n=',I2,' k=',I2,' Old & new norma=',2F10.6)
            ibad=ibad+1
          end if
        end do
        if (ibad.GT.0) then
          write(*,*) ibad,
     >       ' orbitals with bad normalization found'
        end if
        if (ortho) then
          write(*,*) ' #### All orbitals are orthonormalized!'
        else
          if (norm) then
            write(*,*) ' #### All orbitals are normalized!'
          end if
        end if

        if (i8.GT.iix) then
          write(*,25) iix,Rn(iix),i4,Rn(i4),i6,Rn(i6),i8,Rn(i8)
 25       format(4x,'Orbitals stretch out of the box ',i3,2x,f6.2,':',
     >         /(4x,i3,2x,f6.2))
          iix=i8
          Rmax=Rn(iix)
          write(*,35) iix,Rmax
 35       format(4x,'New values: iix=',i4,' Rmax=',f7.3)
          read(*,*)
        end if
       Return
      end
C     =================================================
      subroutine NewGrid
      implicit real*8 (a-h,o-z)
      include "basx.par"
      include "hfd.par"
      include "phys.par"
       common /A/A/Rnuc/Rnuc/cnuc/cnuc/tnuc/tnuc/Rmax/Rmax/MaxT/MaxT
     >        /ii/ii/iin/iin/iix/iix/lmax/lmax/Kt/Kt/H/H
       common /R/R(IP6)/V/V(IP6)/Rn/Rn(IPx6)/Vn/Vn(IPx6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        fermi=1.d-13/DPrb
        a3=A**0.3333333d0

c Parameter tnuc is 0.9-0.1 skin width. It is linked to Woods-Sakson
c parameter t:
        t=tnuc/(4.d0*dlog(3.d0))

c rnuc1 gives approx. same rms raduis as Fermi distribution
c (see ./Work/Maple/fermi1.pdf) :

        rnuc1=dsqrt(cnuc**2+7*3.14159**2/3*t**2)*fermi

c rnuc2 is what is used in hfd:

        rnuc2=(1.115d0*a3+2.151d0/a3-1.742d0/A)*fermi

        write(*,*) 'Rnuc here ',rnuc1,', expected ',rnuc2
        write(*,*) '(cnuc =',cnuc*fermi,')'

        do i=1,ii
          i1=i
          if (Rn(i1).GT.rnuc1) goto 100
        end do

 100    d1=Rn(i1)-rnuc1
        d2=rnuc1-Rn(i1-1)
        if (d1.GT.d2) i1=i1-1
        Rnuc=Rn(i1)

        k=0
        do i=i1,iix
          k=k+1
          R(k)=Rn(i)
          V(k)=Vn(i)
        end do
        iin=i1
        iix=k

        write(*,5) iix, R(1),R(iix),iin-1,ii-iin-iix+1
 5      format(/' New grid has ',I3,' nodes; R1=',F8.6,' Rmax=',F8.3,
     >         /i3,' nodes inside nucleus are skipped'
     >         /i3,' nodes above Rmax are skipped')

        if (iix.GT.IP6-MaxT-5) then
          write(*,*) ' IP6=',IP6,' is too small'
          stop
        end if

       Return
      end
C     =================================================
      subroutine NewBasis
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/Nso/Nso/Z/Z/A/A/let/let(6)/FNAME1/FNAME1
       common /iix/iix/Kt/Kt/H/H/Rnuc/Rnuc/npwx/npwx
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /R/R(IP6)/V/V(IP6)/Qq/Qq(IPs)

       common /g0/g0(IPx6,IPxo,IPpw)/f0/f0(IPx6,IPxo,IPpw)
     >        /En0/En0(2*IPxo,IPpw)/Pn/Pn(IPx6)/Qn/Qn(IPx6)
     >        /Norb/Norb(IPxo,IPpw)

       dimension P(IP6),Q(IP6),P1(IP6),Q1(IP6)

       character*1 let, FNAME1*12
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(12,file=FNAME1,access='DIRECT',
     >       status='UNKNOWN',recl=2*IP6*IPmr)
        close(12,status='DELETE')
        open(12,file=FNAME1,access='DIRECT',
     >       status='NEW',recl=2*IP6*IPmr)

C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,*)
        write(*,*) ' Forming file ',FNAME1,'...'
        call HeadRec

        write(*,*)

        do m=1,Ns
          nm=Nn(m)
          lm=Ll(m)
          km=Kk(m)
          jm=Jj(m)
          ind=lm+iabs(km)
          npw=nm-lm

          if (Norb(npw,ind).NE.m) then
            write(*,*) ' Wrong QN for orbital ',m
            write(*,*) ' n=',nm,' l=',lm,' kappa=',km
            stop
          end if

          do i=1,IPx6
            Pn(i)=g0(i,npw,ind)
            Qn(i)=f0(i,npw,ind)
          end do

          call Convert(P,Q,P1,Q1,km) ! converts orbital to HFD format

          P(iix+1)=-En0(npwx+nm-lm,ind)
          call WriteF (12,m+4,P,Q,2)
          call WriteF (12,m+Ns+4,P1,Q1,2)

          write(*,15) m,nm,let(lm+1),jm,-P(iix+1),P(iix+4),Qq(m)
 15        format(i4,2x,i2,a1,i2,'/2','   E=',f16.6,
     >            '  gamma=',f4.1,'  qq=',f5.1)

        end do

       close(12)
       write(*,*) ' File ',FNAME1,' formed.'
c       read(*,*)
       Return
      end
C     =================================================
      subroutine HeadRec
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns/Nso/Nso/Z/Z/A/A/npwx/npwx/FNAME1/FNAME1
       common /iix/iix/Kt/Kt/H/H/Rnuc/Rnuc/Rmax/Rmax
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /R/R(IP6)/V/V(IP6)/Qq/Qq(IPs)


       dimension IQN(4*IPs),Qq1(IPs),P(IP6),Q(IP6),
     >           P1(IP6),Q1(IP6),PQ(4*IP6)

       equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
       equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),
     >        (P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))
       character*12 FNAME1

        do i=1,4*IP6
          PQ(i)=0.d0
        end do

        PQ(1) = Z
        PQ(2) = Ns
        PQ(3) = iix
        PQ(4) = R(1)
        PQ(5) = Rmax
        PQ(6) = h
        PQ(9) = Kt
        PQ(12)= A
        PQ(13)= Rnuc
        PQ(14)=-1.0d0  ! = average over relativistic configuration
        PQ(20)= 0.98765d0 ! = long basis variant

        write(*,5) Z,A,Ns,iix,Rnuc,R(1),Rmax,h,kt
 5      format(/' Z=',f5.1,' A=',f5.1,' Ns=',i4,' ii=',i4,/' Rnuc=',
     >        f8.6,' R1=',f8.6,' R2=',f8.3,' h=',f6.4,' kt=',i1)

        do ni=1,Ns
          IQN(4*ni-3)=Nn(ni)
          IQN(4*ni-2)=Ll(ni)
          IQN(4*ni-1)=Kk(ni)
          IQN(4*ni)  =Jj(ni)
          if (ni.LE.Nso) then
            Qq1(ni) = Jj(ni)+1
          else
            Qq1(ni) = 0
          end if
          Qq(ni)=Qq1(ni)
        end do

        call WriteF (12,1,P,Q,2)
        call WriteF (12,2,R,V,2)
        call WriteF (12,3,P1,Q1,2)

       write(*,*) ' Head records of ',FNAME1,' formed.'
       Return
      end
C     =================================================
      subroutine WriteF (kan,record,v1,v2,nrec)
      include "hfd.par"
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     fortran-77, ms-dos version
c     - - - - - - - - - - - - - - - - - - - - - - - - -
       integer record,nrec
       real*8 v1(IP6),v2(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        write(kan,rec=nr1) (v1(i),i=1,ii)
        if (nrec.eq.2) write(kan,rec=nr2) (v2(i),i=1,ii)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
       return
      end
C     =================================================
      subroutine Convert(P,Q,P1,Q1,kap)
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /ii/ii/iin/iin/iix/iix/Kt/Kt/H/H
     >        /Pn/Pn(IPx6)/Qn/Qn(IPx6)/Rn/Rn(IPx6)/V/V(IP6)
       dimension P(IP6),Q(IP6),P1(IP6),Q1(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        h60=60*h

        if (iin.LE.3) then
          write(*,*) ' Short tail. Can not calculate derivative:'
          write(*,*) ' iin=',iin,' need ',3
          stop
        end if

        irt=0
c        if (ii.LT.iin+iix+3-1) irt=3   ! asymmetric derivative

        do i=1,IP6
          P(i)=0.d0
          Q(i)=0.d0
          P1(i)=0.d0
          Q1(i)=0.d0
        end do

        do i=1,iix-irt
          i1=i+iin-1
          P(i)=Pn(i1)
          Q(i)=Qn(i1)

          P1(i)=(Pn(i1+3)-9*(Pn(i1+2)-Pn(i1-2)) ! 7 node symmetric derivative
     >         +45*(Pn(i1+1)-Pn(i1-1))-Pn(i1-3))  !# (borrowed from dif.inc)
     >          /(h60*V(i))
          Q1(i)=(Qn(i1+3)-9*(Qn(i1+2)-Qn(i1-2))
     >         +45*(Qn(i1+1)-Qn(i1-1))-Qn(i1-3))
     >          /(h60*V(i))
        end do

c>>>>> Asymmetric 7-node differentiation for last three nodes:
        if (irt.GT.0) then
          P1(iix)=-(-147*Pn(iix)+360*Pn(iix-1)
     >       -450*Pn(iix-2)+400*Pn(iix-3)-225*Pn(iix-4)
     >       +72*Pn(iix-5)-10*Pn(iix-6))/(h60*V(iix))

          P1(iix-1)=-(-10*Pn(iix)-77*Pn(iix-1)
     >       +150*Pn(iix-2)-100*Pn(iix-3)+50*Pn(iix-4)
     >       -15*Pn(iix-5)+2*Pn(iix-6))/(h60*V(iix-1))

          P1(iix-2)=-(2*Pn(iix)-24*Pn(iix-1)-35*Pn(iix-2)
     >       +80*Pn(iix-3)-30*Pn(iix-4)+8*Pn(iix-5)
     >       -Pn(iix-6))/(h60*V(iix-2))
          Q1(iix)=-(-147*Qn(iix)+360*Qn(iix-1)
     >       -450*Qn(iix-2)+400*Qn(iix-3)-225*Qn(iix-4)
     >       +72*Qn(iix-5)-10*Qn(iix-6))/(h60*V(iix))

          Q1(iix-1)=-(-10*Qn(iix)-77*Qn(iix-1)
     >       +150*Qn(iix-2)-100*Qn(iix-3)+50*Qn(iix-4)
     >       -15*Qn(iix-5)+2*Qn(iix-6))/(h60*V(iix-1))

          Q1(iix-2)=-(2*Qn(iix)-24*Qn(iix-1)-35*Qn(iix-2)
     >       +80*Qn(iix-3)-30*Qn(iix-4)+8*Qn(iix-5)
     >       -Qn(iix-6))/(h60*V(iix-2))

          write(*,775) (P1(i),i=iix-4,iix),(Q1(i),i=iix-4,iix)
 775      format(5f12.6)
          read(*,*)
        end if

        P(iix+4)=iabs(kap)
        P1(iix+4)=iabs(kap)-1
        call Origin(P,P1,kap)
        call Origin(Q,Q1,-kap)

        tmax=0.d0
        do i=iix+iin,ii
          t=(Pn(i))**2+(Qn(i))**2
          if (t.GT.tmax) then
            tmax=t
            i1=i
          end if
        end do

        if (tmax.GT.1.d-6) then
          write(*,*) ' Tail outside the box: ',tmax,' i=',i1
          read(*,*)
        end if

      Return
      end
C     =================================================
      subroutine Origin(P,CP,kap)               ! based on dif.inc
      implicit real*8 (a-h,o-z)
      include "hfd.par"                         !# we match Taylor expansion to
       common /iix/iix/MaxT/MaxT/Ierr/Ierr/H/H  !# P(1),P(2), CP(1)=dP/dr(1),
       common /R/R(IP6)/V/V(IP6)                !# and CP(2)=dP/dr(2)
       dimension P(IP6),CP(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=iix+5,iix+5+MaxT
          P(i)=0.d0
          CP(i)=0.d0
        end do

        if (kap.LT.0) then           !### for positive kappa expansion
          k=0                        !#### goes over odd powers
          rn=1.d0
        else
          k=1
          rn=R(1)
        end if

        igam=iabs(kap)
        ig=igam+k

 1      p1=P(1)*rn/R(1)**ig       !# (see 26/10/06)
        p2=P(2)*rn/R(2)**ig       !#  with addition
        d1=CP(1)*rn/R(1)**(ig-1)  !#   (19/05/08)
        d2=CP(2)*rn/R(2)**(ig-1)  !#
        r1=R(1)*R(1)
        y=R(2)*R(2)/r1

        p21=(p2-p1)/(y-1)
        p11=(d1-ig*p1)/2
        p22=(d2-ig*p2)/(2*y)
        p2111=(p21-p11)/(y-1)
        p2211=(p22-p11)/(2*(y-1))

        c3=2*(p2211-p2111)/(y-1)
        c2=p2111-(y+2)*c3
        c1=p21-(y+1)*c2-(y*y+y+1)*c3
        c0=p1-c1-c2-c3

        P(iix+k+5) =c0
        P(iix+k+7) =c1
        P(iix+k+9) =c2
        P(iix+k+11)=c3

        CP(iix+k+5) =ig*c0
        CP(iix+k+7) =(ig+2)*c1
        CP(iix+k+9) =(ig+4)*c2
        CP(iix+k+11)=(ig+6)*c3

        rk=(R(2)/R(1))**k
        pr1=(c0+c1+c2+c3)*R(1)**igam / P(1)               !#
        dr1=(ig*c0+(ig+2)*c1+(ig+4)*c2                    !#
     >     +(ig+6)*c3)*R(1)**(igam-1) / CP(1)             !# all these parameters
        pr2=rk*(c0+c1*y+c2*y*y+c3*y**3)*R(2)**igam / P(2) !# should be equal to 1
        dr2=rk*(ig*c0+(ig+2)*c1*y+(ig+4)*c2*y*y           !#
     >     +(ig+6)*c3*y**3)*R(2)**(igam-1) / CP(2)        !#


        ier=1.d6*(dabs(pr1-1.d0)+dabs(dr1-1.d0)
     >     +dabs(pr2-1.d0)+dabs(dr2-1.d0))


        if (ier.GE.1) then
          write( *,5) igam,kap,pr1,dr1,pr2,dr2
 5        format(4X,'Expansion error for gam = ',I2,' kap = ',I2,
     >          /4X,'pr1 = ',E15.7,' dr1 = ',E15.7
     >          /4X,'pr2 = ',E15.7,' dr2 = ',E15.7
     >          /4X,'shall we quit? ')
          read(*,*) iyes
          if (iyes.NE.0) stop
          Ierr=Ierr+ier
        end if

       return
      end
C     =================================================
      subroutine TestQN
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /Ns/Ns
       common /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
        write(*,*)' Check QN:'
        do i=1,Ns
          write(*,5) i,Nn(i),Ll(i),Jj(i),Kk(i)
5         format(5i4)
        end do
        read(*,*)
       Return
      end
C     =================================================
      subroutine BassInp       ! forms file BASS.INP
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "basx.par"
      include "hfd.par"
       common /name/name/Ns/Ns/Nso/Nso/Z/Z/A/A/npwx/npwx
     >        /FNAME2/FNAME2
       common /iix/iix/Kt/Kt/H/H/Rnuc/Rnuc/Rmax/Rmax/let/let(6)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /R/R(IP6)/V/V(IP6)/Qq/Qq(IPs)
       character*12 FNAME2

       Character name*4,let*1
       dimension IQN(4*IPs),Qq1(IPs)


        open (11,file=FNAME2,status='UNKNOWN')
        close(11,status='DELETE')
        open (11,file=FNAME2,status='NEW')

        Nv=Ns-Nso
        Ksg=1
        Kdg=1
        Kkin=1

        i1=Nso+1  ! first orbital for diagonalization
        n1=Nn(i1)
        l1=Ll(i1)
        j1=Jj(i1)

        i2=Nso    ! last frozen orbital
        n2=Nn(i2)
        l2=Ll(i2)
        j2=Jj(i2)

        kout=0
        kbrt=0

        write(11,'(1X,A4)') name
        write(11,5) Z,A,Nso,Nv
 5      format('  Z =',F5.1,/' Am =',F5.1,/' Nso=',I5,
     >   '# number of core orbitals (defines DF operator)',
     >   /' Nv =',I5,'# number of valence & virtual orbitals')

        write(11,15) Ksg,Kdg
 15     format(' Ksg=',I5,
     >   '# defines Hamiltonian: 1-DF, 3-DF+Breit',
     >   /' Kdg=',I5,'# diagonalization of Hamiltonian (0=no,1,2=yes)')

        write(11,25) n1,let(l1+1),j1,n2,let(l2+1),j2
 25     format(' orb=',I2,A1,I2,'# first orbital for diagonalization',
     >   /' Kkin    1# kinetic balance (0,1,or 2)',
     >   /' orb= 0s 1# first orbital to apply kin.bal.',
     >   /' orb=',I2,A1,I2,'# last frozen orbital')

        write(11,35) Kout,Kbrt
 35     format(' Merr    0# allowed number of warnings.',
     >   /' lst= 0s 1# last orbital to be kept in the basis set',
     >   /'kout=',I2,'   # detail rate in the output',
     >   /'kbrt=',I2,'   # 0,1,2 - Coulomb, Gaunt, Breit',
     >   /'----------------------------------------------------------')

        do i=1,Ns
          if (i.LE.Nso) then
            nq=Jj(i)+1
          else
            nq=1
          end if
          x=0.1d0*Nn(i)+0.01d0*Ll(i)+0.0001d0*nq
          if (Jj(i).LT.2*Ll(i)) x=-x
          Qq1(i)=x
        end do

        write(11,45) (Qq1(i),i=1,Nso)
 45     format(6(4X,F7.4))

        do i=1,Nv
          if (i.EQ.1.OR.Ll(i+Nso).EQ.0) write(11,*)
          write(11,55) i,Qq1(i+Nso),Qq1(i+Nso)
 55       format(I3,1X,F7.4,'  3 ',F7.4)
        end do

        write(11,65)
 65     format('>>>>>>>>>> END <<<<<<<<<<')
        write(*,*) ' File ',FNAME2,' is formed'

        close(11)
       Return
      end
