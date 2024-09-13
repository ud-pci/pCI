c    =====================   04/12/95   ============================
                               ! last update  7-Feb-2014 (Newark)
C  SP:  calculation of M3 amplitudes is added
                               !      update 13-Dec-2013 (Newark)
C  SP:  "Gaun" is replaced by relativistic M1
                               !      update 10-Dec-2013 (Newark)
C  SP:  calculation of E3 and M2 amplitudes is added
                               !      update 17-Feb-2012 (Newark)
C  SP:  calculation of E2 amplitudes is added
C    ===============================================================
      Program RPA              !### previous update on 01/11/11
c#    Solution of the RPA equation in the matrix form
c#    for 10 operators including:
c#       E1 amplitudes in L and V form, E2 amplitudes,
c#       hfs constants A and B,
c#       EDM, PNC, AM and MQM,
c#    Program can be easily adopted for new operators.
c#    Breit corrections can be added to RPA equations (Kbrt=1).

c  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c  || RPA equations are solved SELF-CONSISTENTLY for the ||
c  || following Nhf shells:                              ||
c  ||    Nso-Nhf+1, Nso-Nhf+2, ..., Nso                  ||
c  || Note that sums over virtual states include         ||
c  ||    Nsh+1, Nsh+2, ...,Ns                            ||
c  || If Nsv <= Nsh, then Core-Core MEs for shells       ||
c  ||    Nsv... Nsh are found.                           ||
c  || Then, RPA corrections are found for valence shells ||
c  ||    Nsv, Nsv+1, ..., Nmax                           ||
c  || with additional restriction: max(li,lk) <= Lmax    ||
c  || Valence Eqs have the same form as core Eqs         ||
c  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ns/Ns/Nsv/Nsv/Nsh/Nsh/Nr/Nr/Nr1/Nr1/Kl/Kl/Jmax/Jmax
     >        /Kval/Kval/Nmax/Nmax/Lmax/Lmax/Kout/Kout/Small/Small
     >        /k_denom/k_denom
       common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Lj/Lj(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)
     >        /Kl1/Kl1(IPcs)
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /P/P(IP6)/Q/Q(IP6)/A/A(IP6)/B/B(IP6)
       common /ipmr/ipmr
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        call recunit
        itest=IPs
        call Check_Dim(itest) ! IPs must be the same as in conf.par
        Jmax=10 !### number of terms of the Taylor expansion at the origin
        Small=1.d-8
        write(*,*)' RPA for core (1) or for valence (2):'
        read(*,*) kv
        if ((kv-1)*(kv-2).NE.0) stop
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=10,file='MBPT.INP',status='OLD')
        open(unit=11,file='RPA.RES',status='UNKNOWN')
        open(unit=12,file='HFD.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*IPmr,err=700)
        if (kv.EQ.2) open(unit=13,file='CONF.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*IPmr,err=710)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Kout=2
        k_denom=1       !### two variants of valence denominators
                        !#### are coded (see SolVal).
        call Input
        call Init
        call CalcWig0
        if (kv.EQ.2.AND.k_denom.EQ.2) call DefEv(13,Lmax)
        close(unit=10)
        Kout=0
        do i=1,IPcs
          if (Kl1(i).EQ.1) then
            Kl=i
            call Rinit
            call RHS
            if (kv.EQ.2) then
              call Rinput(ierr)      !# ierr=0 means that solution is
            else                     !## known from the file RPA_N.INT
              ierr=1
            end if
            if (ierr.NE.0) then
              call LHS
              call Iteration
              call Check
              call Rsave(1)
            end if
            if (kv.EQ.1) then
              call Prin
            else
              if (Nsv.LE.Nsh) then   !# Core-Core MEs are used in SolVal
                call CoreCore        !## these MEs are added to Core-Valence
                Nr=Nr1               !### MEs in arrays Fp and Fm, Nr being
              end if                 !#### the total number of MEs.
              call Vinit
              call RHSv
              call SolVal
              call Rsave(kv)
            end if
          end if
        end do
       stop
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
        write(11,75)
75      format(/2X,'file HFD.DAT is absent'/)
       stop
C     - - - - - - - - - - - - - - - - - - - - - - - - -
710     write( *,85)
        write(11,85)
85      format(/2X,'file CONF.DAT is absent'/)
       stop
      end
C     =================================================
      include "wig.inc"
      include "readf.inc"
      include "readff.inc"
      include "decomp.inc"
      include "sint1.inc"
      include "yk.inc"
      include "ykt.inc"
      include "eval.inc"
      include "breit.inc"
      include "check_dim.inc"
      include "rec_unit.inc"
c     =================================================
      Subroutine Input
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Nss/Nss/Nsh/Nsh/Kt2/Kt2/Nhf/Nhf/Lmax/Lmax/Nso/Nso
     1        /Nsv/Nsv/Nmax/Nmax/Kmg/Kmg/Omega/Omega/Kex/Kex/Kbrt/Kbrt
     2        /Alet/Alet(IPcs)/Kl1/Kl1(IPcs)
       character*4 Alet,lbl
       character*7 blet(2)
       character*7 blet_br(3)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        data blet/'WITHOUT','  with '/
        data blet_br/'WITHOUT','  with ','  with '/
        Alet(1) ='A_hf'
        Alet(2) ='B_hf'
        Alet(3) ='E1_L'
        Alet(4) ='EDM '
        Alet(5) ='PNC '
        Alet(6) ='E1_V'
        Alet(7) ='ANM '
        Alet(8) ='MQM '
        Alet(9) =' M1 '
        Alet(10)=' E2 '
        Alet(11)=' E3 '
        Alet(12)=' M2 '
        Alet(13)=' M3 '
C       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    read (10,5,end=1000) lbl
 5      format(A4)
        if (lbl.NE.'MBPT') goto 200
        read (10,15) Nso     !# - CI core
        read (10,15) Nsh     !# - defines SCF field
        read (10,15) Nss     !# - Last virtual shell for RPA
        read (10,15) Nsv     !# - First valence shell
                             !##   (normaly Nso+1, but can be arbitrary)
        read (10,15) Nmax    !# - Last valence shell for effective
                             !##   radial integrals
        read (10,15) Lmax    !# - max L for radial integrals
                             !##   in the valence space
        read (10,15)         !# - max multipolarity (not used)
        read (10,15) Kt2     !# - accuracy for Coulomb radial integrals
        read (10,15) Kbrt    !# - Breit corrections
 15     format(5X,I3)
        n=0
        m=0
        rewind(10)
 210    read (10,5,end=1000) lbl
        if (lbl.NE.'RPA ') goto 210
        do i=1,IPcs
           read (10,15) Kl1(i)
           if (Kl1(i).EQ.1) then
              n=n+1
           else
              m=m+iabs(Kl1(i))
           end if
        end do
        if (n.EQ.0.OR.m.NE.0) then
           write(*,*) 'Ambigious input file:'
           write(*,*) 'RHS operators not defined'
           stop
        end if
        read (10,15)
        read (10,15) Nhf        !# - SCF procedure includes Nhf shells
        read(10,25) Kmg,Omega   !# - frequency of the perturbation
                                !##   (ignored if Kmg=0)
 25     format(5X,I3,/6X,F9.6)
        if (Kmg.EQ.0) Omega=0.d0
        read (10,15) Kex        !# - skip(0) or include(1) exchange
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*, 35) blet(Kex+1),Omega,blet_br(Kbrt+1)
        write(11,35) blet(Kex+1),Omega,blet_br(Kbrt+1)
 35     format(1X,56('#'),/' ##### RPA ',A7,' exchange for frequency ',
     >       F9.6,' #######',/' ####### and ',A7,' Breit corrections ',
     >       'to MEs #########',/1X,56('#'))
       return
 1000   write (*,*) ' Can not read file MBPT.INP'
       stop
      end
C     =================================================
      Subroutine Init
      implicit real*8 (a-h,o-z)
      INCLUDE "rpa.par"
      INCLUDE "phys.par"
      INCLUDE "hfd.par"
       common /Ns/Ns/Ii/Ii/Kt1/Kt1/Kt2/Kt2/Z/Z/Nhf/Nhf
     1        /Cl/Cl/H/H/Rnuc/Rnuc
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Eps/Eps(IPs)
     >        /R/R(IP6)/V/V(IP6)
       common /AL/AL/BT/BT/R2/R2	! cAB
       dimension p(IP6),q(IP6),p1(IP6),q1(IP6),pq(4*IP6)
       logical longbasis
       dimension IQN(4*IPs),Qq1(IPs)
       equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
       equivalence (p(1),pq(1)), (q(1),pq(IP6+1)),
     >      (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0         !# - small number
        Cl=DPcl           !# - speed of light is taken from "phys.par"
        call ReadF (12,1,p,q,2)
        call ReadF (12,2,R,V,2)
        call ReadF (12,3,p1,q1,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Z   =pq(1)
        Ns  =pq(2)+c1
        Ii  =pq(3)+c1
        H   =pq(6)
        Kt1 =pq(9)+c1
        Rnuc=pq(13)
        al = pq(8) ! cAB
        bt = pq(7)
        r2 = pq(5)
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
        write( *,25) Z,Ii,Ns,Nhf,Rnuc,Kt1,Kt2
        write(11,25) Z,Ii,Ns,Nhf,Rnuc,Kt1,Kt2
 25     format (4X,'Z =',F6.2,4X,'II =',I4,4X,'NS =',I3,4X,'Nhf =',I3,
     1          /4X,'RNUC=',E11.4,2X,'KT1=',I2,2X,'KT2=',I2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          write(11,*) ' Using variant for long basis '
          do ni=1,Ns
            Nn(ni)=IQN(4*ni-3)
            Ll(ni)=IQN(4*ni-2)
            Kk(ni)=IQN(4*ni-1)
            Lj(ni)=IQN(4*ni)
          end do
        else
          if=20
          do ni=1,Ns
            if=if+1
            Nn(ni)=pq(if)+c1
            if=if+1
            Ll(ni)=pq(if)+c1
            if=if+3
            c2=dsign(c1,pq(if))
            Kk(ni)=pq(if)+c2
            if=if+1
            c2=dsign(c1,pq(if))
            Lj(ni)=pq(if)+c2
          end do
        end if
        do ni=1,Ns
           call ReadF (12,ni+4,p,q,1)
           Eps(ni)=-p(Ii+1)
        end do
        write(11,35) (i,Eps(i),i=1,Ns)
 35     format(5X,'Hartree-Fock energies of the orbitals:',
     >       /6(I4,F8.2))
       Return
      end
C     =================================================
      Subroutine Rinit     !### definition of RPA space
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nss/Nss/Nsh/Nsh/Nhf/Nhf/Nr/Nr/Nr1/Nr1
     1        /Nso/Nso/Nsv/Nsv/Ip/Ip/Ij/Ij/Kl/Kl/Ita/Ita
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /Alet/Alet(IPcs)
       dimension ip1(IPcs),ij1(IPcs),ita1(IPcs)
       character*4 Alet
       character*1 let(9)
        data let /'s','p','d','f','g','h','i','k','l'/
        data ip1 /+1,+1,-1,-1,-1,-1,-1,-1,+1,+1,-1,-1,+1/
        data ij1 / 2, 4, 2, 0, 0, 2, 2, 4, 2, 4, 6, 4, 6/
        data ita1/+1,+1,+1,+1,-1,-1,-1,+1,+1,+1,+1,+1,+1/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     selection rules for given Kl:
        ip = ip1(Kl)            !# - parity
        ij = ij1(Kl)            !# - 2*rank
        ita= ita1(Kl)           !# - real(+1)/imaginary(-1)
        write( *,5) Alet(Kl),ip,ij/2,(1-Ita)/2
        write(11,5) Alet(Kl),ip,ij/2,(1-Ita)/2
 5      format (///4X,'RHS operator ',A4,': parity =',I2,'; rank =',
     >       I2,' phase =',I2,'*pi/2',/4X,40('-'))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Nss=min(Nss,Ns)
        Nr=0                    !# - Dimension of RPA space
        do i1=1,Nhf             !# Note that core loop ALWAYS starts with
           i=Nso-i1+1           !## Nso and goes down. That means that
           li=Ll(i)             !### shells Nso+1,..Nsh are not included!
           ji=Lj(i)
           num=0
           do k=Nsh+1,Nss       !# Virtual loop starts with Nsh+1,
              lk=Ll(k)          !## so shells Nso+1,..Nsh are skipped in
              jk=Lj(k)          !### both loops
              ipp=Isig(li)*Isig(lk)*ip
              jm=iabs(ji-jk)
              js=ji+jk
              if (ipp.EQ.1.AND.jm.LE.ij.AND.js.GE.ij) then
                 Nr=Nr+1
                 num=num+1
                 Iri(Nr)=i
                 Irk(Nr)=k
              end if
           end do
           write(*, 15) Nn(i),let(li+1),Lj(i),num
           write(11,15) Nn(i),let(li+1),Lj(i),num
 15        format(4X,I3,A1,I2,'/2  num =',I3)
           if (Nr.GT.IPr) then
              write(*,*) 'Dimension of RPA space is already ',Nr
              write(*,*) 'current core shell is',i1,' of total',Nhf
              stop
           end if
        end do
        write(*, 25) Nr,Nss
        write(11,25) Nr,Nss
 25     format(4X,40('-'),/4X,'Dimension of RPA space =',
     >       I5,' Nss =',I3)
        if (Nr.EQ.0) stop
        Nr1=Nr
        nmin=min(Nso,Nsv-1)
        if (Nsv.LE.Nsh) then         !#  These MEs do not enter RPA equations
           ix=Nhf-Nso+nmin           !##  for core, but are used to calculate
           do i1=1,ix                !###  RPA MEs for valence shells.
              i=nmin+1-i1            !####  Note, that Nsv.LE.Nso is allowed
              li=Ll(i)
              ji=Lj(i)
              do k=Nsv,Nsh
                 lk=Ll(k)
                 jk=Lj(k)
                 ipp=Isig(li)*Isig(lk)*ip
                 jm=iabs(ji-jk)
                 js=ji+jk
                 if (ipp.EQ.1.AND.jm.LE.ij.AND.js.GE.ij) then
                    Nr1=Nr1+1
                    Iri(Nr1)=i
                    Irk(Nr1)=k
                 end if
              end do
           end do
           write(*, 35) Nsv,Nsh,Nr1-Nr
           write(11,35) Nsv,Nsh,Nr1-Nr
 35        format(' First valence shell',I4,' is LE Nsh=',I4,
     >        ', so',I4,' core-core MEs required')
           if (Nr1.GT.IPv) then
             write (*,*) ' IPv=',IPv,' < Nr1=',Nr1
             stop
           end if
        end if
       Return
      end
C     =================================================
      Subroutine RHS   !### Reduced ME for RHS operator
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Nso/Nso/Nsh/Nsh/Nr1/Nr1/Kl/Kl/Ii/Ii
     1        /Z/Z/Ita/Ita/Kt/Kt/Kt1/Kt1
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Jj/Jj(IPs)
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /F0/F0(IPv)
     >        /Ro/Ro(IP6)/C/C(IP6)/R/R(IP6)
     >        /UP/Pa(IP6),Pc(IP6),Pd(IP6),Pb(IP6)
     >        /DOWN/Qa(IP6),Qc(IP6),Qd(IP6),Qb(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     CORE ELECTRON DENSITY AND Z_eff(R):
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Kt=Kt1
        ih=2-Kt
        do i=1,IPs     !### file breit.inc uses Jj instead of Lj
          Jj(i)=Lj(i)
        end do
        if (Kl.EQ.4) then
           do i=1,Ii,ih
              C(i)=0.d0
           end do
           gab=10.d0
           do in=1,Nso    ! Nsh is changed to Nso to agree with dtm (1.11.11)
              qe=Lj(in)+1
              call ReadFF(12,in+4,pa,qa,2)
              ga=(pa(Ii+4))*2
              if (gab.GT.ga) gab=ga
              do i=1,Ii,ih
                 C(i)=C(i)+((pa(i))**2 +(qa(i))**2)*qe
              end do
           end do
           Ro(1)=Z-C(1)*R(1)/(gab+1)
           do i=1+ih,Ii,ih
              Ro(i)=Ro(i-ih)-(C(i-ih)+C(i))/2*(R(i)-R(i-ih))
           end do
           Write(11,15) (Ro(i),i=1,Ii,10*ih)
 15        format(5X,'Ro(r):',/(12F5.1))
        end if
C - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ia0=0
        do i=1,Nr1
           ia=Iri(i)
           na=Nn(ia)
           la=Ll(ia)
           ja=Lj(ia)
           ka=Kk(ia)
           if (ia0.NE.ia) then
              call ReadFF(12,ia+4,pa,qa,2)
              ia0=ia
           end if
           ib=Irk(i)
           nb=Nn(ib)
           lb=Ll(ib)
           jb=Lj(ib)
           kb=Kk(ib)
           call ReadFF(12,ib+4,pb,qb,2)
           F0(i)=Fnc(la,ja,ka,na,lb,jb,kb,nb,pa,qa,pb,qb)
        end do
       Return
      end
C     =================================================
      Function Fnc(la,ja,ka,na,lb,jb,kb,nb,pa,qa,pb,qb)
c#    Radial integral for ME <b|F|a> of the RHS operator.
c#    Note the inverse order of indeces in comparisson with DTM.
c#    That changes sign of PNC, E1V and AM operators.
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Kl/Kl/Ii/Ii/Kt/Kt/Dint/Dint/Cl/Cl/Z/Z
     >        /R/R(IP6)/V/V(IP6)/C/C(IP6)/Ro/Ro(IP6)
       dimension pa(IP6),qa(IP6),pb(IP6),qb(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        tab=0.d0
        dn=0.d0
        ih=2-Kt
        gab=pa(Ii+4)+pb(Ii+4)
        kab=iabs(ka)+iabs(kb)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     1                  DIPOLE HFS:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.1) then
           do i=1,Ii,ih
              C(i)=-(pa(i)*qb(i)+qa(i)*pb(i))/(R(i)**2)
           end do
           C(Ii+4)=gab-2
           call Sint1(tab)
           call NclInt(kab-1,kab+1,-1.d0,-1.d0,pa,qb,qa,pb,dn)
           tab=tab-Dint+dn
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     2                  QUADUROPE HFS:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.2) then
           do i=1,Ii,ih
              C(i)= (pa(i)*pb(i)+qa(i)*qb(i))/(R(i)**3)
           end do
           C(Ii+4)=gab-3
           call Sint1(tab)
           call NclInt(kab-2,kab+2,1.d0,1.d0,pa,pb,qa,qb,dn)
           tab=tab-Dint+dn
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     3             E1 AMPLITUDE (L GAUGE):
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.3) then
           do i=1,Ii,ih
              C(i)= (pa(i)*pb(i)+qa(i)*qb(i))*R(i)
           end do
           C(Ii+4)=gab+1
           call Sint1(tab)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     4             EDM OF THE ELECTRON:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.4) then
           do i=1,Ii,ih
              C(i)= (qa(i)*qb(i)) * Ro(i)/(R(i)**2)
           end do
           C(Ii+4)=gab-2
           call Sint1(tab)
           call NclInt(kab+2,kab+1,0.d0,Ro(1)/R(1),pa,pb,qa,qb,dn)
           tab=tab-Dint+dn
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     5                PNC AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.5.AND.gab.LT.2.5d0) then
           call NclInt(kab-2,kab,-3.d0,+3.d0,pa,qb,qa,pb,dn)
           tab=-dn                   !### minus accounts for i <-> f ###
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     6             E1 AMPLITUDE (V GAUGE):
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.6) then
           lg=max(la,lb)
           is=-Isig((ja+jb)/2+la+lg) !### minus accounts for i <-> f ###
           w3 = FJ6(lb+0.d0,jb/2.d0,0.5d0,ja/2.d0,la+0.d0,1.d0)
           w1 = FJ6(0.5d0,jb/2.d0,lb+0.d0,ja/2.d0,0.5d0,1.d0)/w3
           w2 = FJ6(0.5d0,jb/2.d0,la+0.d0,ja/2.d0,0.5d0,1.d0)/w3
           do i=1,Ii,ih
              C(i)=0.d0
              if(lb.EQ.ja-la) C(i) = C(i) - qa(i)*pb(i)*w1
              if(la.EQ.jb-lb) C(i) = C(i) - pa(i)*qb(i)*w2
           end do
           C(Ii+4)=gab
           call Sint1(tab)
           s1=0.d0
           s2=0.d0
           IF(lb.EQ.ja-la) s1=-w1
           if(la.EQ.jb-lb) s2=-w2
           call NclInt(kab+1,kab,s1,s2,pb,qa,qb,pa,dn)
           tab=tab-Dint+dn
           tab=tab*Cl*is*dsqrt(6.d0/lg)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     7                  ANAPOLE MOMENT
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.7.AND.gab.LT.2.5d0) then
           if (la.EQ.0) then
              s1=-3.d0
              s2=-1.d0
           else
              s1= 1.d0
              s2= 3.d0
           end if
           call NclInt(kab-2,kab,s1,s2,pa,qb,qa,pb,dn)
           tab=-dn                   !### minus accounts for i <-> f ###
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     8                MQM AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.8) then
           do i=1,Ii,ih
              C(i)=-(pa(i)*qb(i)+qa(i)*pb(i))/(R(i))**3
           end do
           C(Ii+4)=gab-3
           call Sint1(tab)
           call NclInt(kab-2,kab+2,-1.d0,-1.d0,pa,qb,qa,pb,dn)
           tab=tab-Dint+dn
        end if
!   Added by SP: 13-Dec-2103
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     9         rel. M1-amplitude  (K_M1=2 in "dtm")
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.9) then
          do i=1,Ii,ih
            C(i)=0.5d0*Cl*(pa(i)*qb(i)+qa(i)*pb(i))*R(i)
          end do
          C(ii+4)=gab+1
          call Sint1(tab)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     10               E2 AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.10) then
          do i=1,Ii,ih
             C(i)= (pa(i)*pb(i)+qa(i)*qb(i))*(R(i)**2)
          end do
          C(Ii+4)=gab+2
          call Sint1(tab)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C     11               E3 AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.11) then
          do i=1,Ii,ih
             C(i)= (pa(i)*pb(i)+qa(i)*qb(i))*(R(i)**3)
          end do
          C(Ii+4)=gab+3
          call Sint1(tab)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
*     12            M2 AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (Kl.EQ.12) THEN
          do i=1,Ii,ih
            C(i)= -2/3.d0*(ka+kb)*Cl*
     >            (pa(i)*qb(i) + qa(i)*pb(i))* R(i)**2
          end do
          C(ii+4)= gab+2
          call Sint1(tab)
        END IF
C     - - - - - - - - - - - - - - - - - - - - - - - - -
*     13            M3 AMPLITUDE:
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (Kl.EQ.13) THEN
          do i=1,Ii,ih
            C(i)= -0.5d0*(ka+kb)*Cl*    ! (-) stands because P= f*r,  Q = -g*r
     >                   (pa(i)*qb(i) + qa(i)*pb(i))* R(i)**3
          end do
          C(ii+4)= gab+3
          call Sint1(tab)
        END IF
C     -------------------------------------------------
        Fnc=tab
       Return
      end
C     =================================================
      Function Factor(la,ja,lb,jb)
c#    Angular factor: <b||F||a> = Factor * Radial_Integral_ba
       implicit real*8 (a-h,o-z)
       common /Kl/Kl/Ij/Ij
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xja=ja*0.5d0
        xla=la
        xjb=jb*0.5d0
        xlb=lb
        lx=max(la,lb)
        xlx=lx
        jmin=min(ja,jb)
        xjn=jmin*0.5d0
        xab=(ja+1)*(jb+1)
        xlt=ja-la
        xij=Ij/2.d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kl.EQ.1) then
           f=Isig((ja+1)/2+la)*dsqrt(xab/(xjn+1))
           if (ja.EQ.jb) f=f*dsqrt((ja+1)/xja)
        end if
        if (Kl.EQ.2 .OR. Kl.EQ.10)
     >     f=Isig((ja+1)/2)*dsqrt(xab*(2*la+1)*(2*lb+1))
     >       *FJ3(xlb,xla,xij,0.d0,0.d0,0.d0)
     >       *FJ6(xjb,xja,xij,xla,xlb,0.5d0)
        if (Kl.EQ.3.OR.Kl.EQ.6) f=Isig((ja-1)/2+lx)*dsqrt(xab*lx)
     >       *FJ6(xlb,xjb,0.5d0,xja,xla,xij)
        if (Kl.EQ.4.OR.Kl.EQ.5) f=dsqrt((2*xja+1))
        if (Kl.EQ.7) f=Isig((ja+1)/2+lx)
        if (Kl.EQ.8) f=dsqrt(xab)*(Isig(lb+lx)*dsqrt(30*xlx)
     >       *FJ9(xla,xlb,1.d0,xja,xjb,xij,0.5d0,0.5d0,1.d0)
     >       +Isig((ja+1)/2)*dsqrt(2.d0/3*(2*lb+1)*(2*xlt+1))
     >       *FJ3(xlt,xlb,xij,0.d0,0.d0,0.d0)
     >       *FJ6(xlb,xjb,0.5d0,xja,xlt,xij))

        if (Kl.EQ.9) then
          is = 1
          k= xjb+lb+1.51d0
          if (k.NE.2*(k/2)) is=-is
          
          if (ja.NE.jb) f= is* dsqrt(xab/(xjn+1))
          if (ja.EQ.jb)
     >      f= is* (2*xja+1)*dsqrt((2*xja+1)/(xja*(xja+1)))
        end if
        
        if (Kl.EQ.11)
     >    f= Isig((ja-1)/2)*dsqrt(xab*(2*la+1)*(2*lb+1))
     >       *FJ3(xlb,xla,xij,0.d0,0.d0,0.d0)
     >       *FJ6(xjb,xja,xij,xla,xlb,0.5d0)

        if (Kl.EQ.12)
     >    f= Isig((jb+1)/2)*dsqrt(xab)
     >       *FJ3(xjb,xja,2.d0,-0.5d0,0.5d0,0.d0)

        if (Kl.EQ.13)
     >    f= Isig((jb+1)/2)*dsqrt(xab)
     >       *FJ3(xjb,xja,3.d0,-0.5d0,0.5d0,0.d0)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Factor = f
        if (f.EQ.0.d0) then
           write(*,*) 'Factor: zero factor for'
           write(*,*) 'la,ja=',la,xja
           write(*,*) 'lb,jb=',lb,xjb
           stop
        end if
       Return
      end
C     =================================================
      Subroutine LHS
c#    Matrices Tp,Tm,Dp,Dm for the left hand side of RPA equations
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Nr/Nr/Omega/Omega/Ita/Ita/Kbrt/Kbrt
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /Tp/Tp(IPr,IPr)/Tm/Tm(IPr,IPr)
     >        /Dp/Dp(IPr)/Dm/Dm(IPr)
     >        /Eps/Eps(IPs)
     >        /UP/P1(IP6),A1(IP6),P2(IP6),A2(IP6)
     >        /DOWN/Q1(IP6),B1(IP6),Q2(IP6),B2(IP6)
       character*13 clmb_br(3)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        clmb_br(1)='Pure Coulomb '
        clmb_br(2)='Coulomb-Breit'
        clmb_br(3)='Coulomb-Breit'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        idel=IPdl
        nu=0
        nu1=0
        write(*, 5) clmb_br(Kbrt+1),Nr*(Nr+1)/2
 5      format(5X,'Formation of ',A13,' matrices Tp and Tm ('
     >       ,I8,' MEs):',/)
        do n=1,Nr
           ni=Iri(n)
           nk=Irk(n)
           call ReadFF (12,ni+4,P1,Q1,2)
           call ReadFF (12,nk+4,P2,Q2,2)
           do m=1,n
              mi=Iri(m)
              mk=Irk(m)
              call ReadFF (12,mi+4,A1,B1,2)
              call ReadFF (12,mk+4,A2,B2,2)
              call Clmb(ni,nk,mi,mk,t1,t2,t3,ff)
              if (Kbrt.GE.1) call Brt(ni,nk,mi,mk,t1,t2,t3)
              Tp(n,m)=((1+Ita)*t1+t2+Ita*t3)*ff
              Tp(m,n)=((1+Ita)*t1+t2+Ita*t3)/ff
              Tm(n,m)=((1-Ita)*t1+t2-Ita*t3)*ff
              Tm(m,n)=((1-Ita)*t1+t2-Ita*t3)/ff
              nu=nu+1
              if (dabs(t1)+dabs(t2)+dabs(t3).NE.0.d0) nu1=nu1+1
              if (nu-idel*(nu/idel).EQ.0) then
                 write(*,25) nu
 25              Format('+',3X,'nu =',I8)
              end if
           end do
        end do
        write(*, 15) clmb_br(Kbrt+1),nu,nu-nu1
 15     format(5X,A13,' matrices Tp and Tm formed (',
     >       I8,' MEs, ', I5,'  zeros).')
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        do n=1,Nr  !### - formation of the diagonal matrices Dp and Dm
           i=Iri(n)
           k=Irk(n)
           s=Eps(i)-Eps(k)
           Dp(n)=s/(s*s-Omega*Omega)
           Dm(n)=Omega/(s*s-Omega*Omega)
        end do
       return
      end
C     =================================================
      Subroutine Clmb(i,k,n,m,t1,t2,t3,ff) !### Reduced Coulomb MEs
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ip/Ip/Ij/Ij/Kex/Kex
     >        /Ll/Ll(IPs)/Lj/Lj(IPs)
     >        /UP/P1(IP6),A1(IP6),P2(IP6),A2(IP6)
     >        /DOWN/Q1(IP6),B1(IP6),Q2(IP6),B2(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        t1=0.d0
        t2=0.d0
        t3=0.d0
        ji=Lj(i)
        jk=Lj(k)
        jn=Lj(n)
        jm=Lj(m)
        kr=Ij/2
        h=0.5d0
c factor ff accounts for the transition: reduced ME --> rad. int.
        ff=Factor(Ll(n),jn,Ll(m),jm)/Factor(Ll(i),ji,Ll(k),jk)
C t1:
c parity selection rule for direct term:
        if (Ip.EQ.Isig(kr)) then
          s1=Isig((ji+jn)/2+1)*Cwig(jk,ji,jm,jn,kr)/(2*kr+1)
          if (s1.NE.0.d0) then
             call CFnc(kr,P2,Q2,P1,Q1)
             t1=s1*Rint(A1,B1,A2,B2)
          end if
        end if
        if (Kex.EQ.0) return
C t2:
        mn=max(jk-jm,jm-jk,ji-jn,jn-ji)/2
        mx=min(jk+jm,ji+jn)/2
        ipp=Isig(Ll(k)+Ll(m))
c parity selection rule for exchange term:
        if (Isig(mn).NE.ipp) mn=mn+1
        do kq=mn,mx,2
           s2=Isig((ji+jn)/2+1+kr+kq)*Cwig(jk,jm,ji,jn,kq)
     >       *FJ6(jk*h,ji*h,kr+0.d0,jn*h,jm*h,kq+0.d0)
           if (s2.NE.0.d0) then
              call CFnc(kq,A1,B1,P1,Q1)
              t2=t2+s2*Rint(P2,Q2,A2,B2)
           end if
        end do
C t3:
        mn=max(jk-jn,jn-jk,ji-jm,jm-ji)/2
        mx=min(jk+jn,ji+jm)/2
        ipp=Isig(Ll(k)+Ll(n))
        if (Isig(mn).NE.ipp) mn=mn+1
        do kq=mn,mx,2
           s3=Isig((ji+jn)/2+1+kr+kq)*Cwig(jk,jn,ji,jm,kq)
     >       *FJ6(jk*h,ji*h,kr+0.d0,jm*h,jn*h,kq+0.d0)
           if (s3.NE.0.d0) then
              call CFnc(kq,P2,Q2,A1,B1)
              t3=t3+s3*Rint(A2,B2,P1,Q1)
           end if
        end do
       return
      end
C     =================================================
      Subroutine Brt(i,k,n,m,t1,t2,t3) !### Reduced Breit MEs. Note that
      INCLUDE "rpa.par"                !#### here there is no parity
      implicit real*8 (a-h,o-z)       !##### selection rule: p=(-1)^k
      INCLUDE "hfd.par"
       common /Ip/Ip/Ij/Ij/Kex/Kex/Kbrt/Kbrt/small/small
     >        /Ll/Ll(IPs)/Lj/Lj(IPs)
     >        /UP/P1(IP6),A1(IP6),P2(IP6),A2(IP6)
     >        /DOWN/Q1(IP6),B1(IP6),Q2(IP6),B2(IP6)
     >        /XJ/xja,xjb,xjc,xjd !### total angular momenta of orbitals
     >        /XL/xla,xlb,xlc,xld !### orbital angular momenta of upper c-s
     >        /YL/yla,ylb,ylc,yld !### orbital angular momenta of lower c-s
       dimension Rint2(10)        !### Radial integrals after regrouping
       real*8, dimension(IP6) :: Pi,Qi,Pk,Qk,Pn,Qn,Pm,Qm
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Kbrt.LT.1) return
!
        Pi = P1; Qi = Q1
        Pk = P2; Qk = Q2
        Pn = A1; Qn = B1
        Pm = A2; Qm = B2
!
        ji=Lj(i)
        jk=Lj(k)
        jn=Lj(n)
        jm=Lj(m)

        li=Ll(i)
        lk=Ll(k)
        ln=Ll(n)
        lm=Ll(m)

        kr=Ij/2
        h=0.5d0

        xja=h*ji
        xjb=h*jn
        xjc=h*jk
        xjd=h*jm

        xla=li
        xlb=ln
        xlc=lk
        xld=lm

        yla=2*xja-xla
        ylb=2*xjb-xlb
        ylc=2*xjc-xlc
        yld=2*xjd-xld

C t1:
        s1=Isig((ji+jn)/2+1)*Cwig(jk,ji,jm,jn,kr)/(2*kr+1)
        if (s1.NE.0.d0) then

          call Mult(xk1,xk2,ipb)
          xin=(kr-xk1-small)*(xk2+small-kr)      !### xk1.LE.kr.LE.xk2
          if (iabs(ipb).EQ.1.AND.xin.GT.0.d0) then
            ds=breit_int(kr,i,Pi,Qi,n,Pn,Qn,k,Pk,Qk,m,Pm,Qm)
            t1=t1+s1*ds
          end if
        end if

        if (Kex.EQ.0) return
C t2:
        xja=h*jm
        xjd=h*ji

        xla=lm
        xld=li

        yla=2*xja-xla
        yld=2*xjd-xld

        do i1=1,IP6
          xi=A2(i1)
          A2(i1)=P1(i1)
          P1(i1)=xi
          xi=B2(i1)
          B2(i1)=Q1(i1)
          Q1(i1)=xi
        end do
        call Mult(xk1,xk2,ipb)
        if (iabs(ipb).EQ.1) then

          mn=max(jk-jm,jm-jk,ji-jn,jn-ji)/2
          mx=min(jk+jm,ji+jn)/2

          do kq=mn,mx,1
            s2=Isig((ji+jn)/2+1+kr+kq)*Cwig(jk,jm,ji,jn,kq)
     >        *FJ6(jk*h,ji*h,kr+0.d0,jn*h,jm*h,kq+0.d0)
            ds=breit_int(kq,m,Pm,Qm,n,Pn,Qn,k,Pk,Qk,i,Pi,Qi)
            if (s2.NE.0.d0) t2=t2+s2*ds
          end do
        end if

C t3:
        xja=h*jn
        xjb=h*jm

        xla=ln
        xlb=lm

        yla=2*xja-xla
        ylb=2*xjb-xlb

        do i1=1,IP6
          xi=A1(i1)
          A1(i1)=P1(i1)
          P1(i1)=xi
          xi=B1(i1)
          B1(i1)=Q1(i1)
          Q1(i1)=xi
        end do
        call Mult(xk1,xk2,ipb)
        if (iabs(ipb).EQ.1) then
          mn=max(jk-jn,jn-jk,ji-jm,jm-ji)/2
          mx=min(jk+jn,ji+jm)/2
          do kq=mn,mx,1
             s3=Isig((ji+jn)/2+1+kr+kq)*Cwig(jk,jn,ji,jm,kq)
     >         *FJ6(jk*h,ji*h,kr+0.d0,jm*h,jn*h,kq+0.d0)
             ds=breit_int(kq,m,Pm,Qm,i,Pi,Qi,k,Pk,Qk,n,Qn,Qn)
            if (s3.NE.0.d0) t3=t3+s3*ds
          end do
        end if

        do i1=1,IP6            !### (P1,Q1) and (P2,Q2) should be
          P1(i1)=A2(i1)        !#### at their places for future use
          Q1(i1)=B2(i1)        !##### while (A1,B1) and (A2,B2) are
        end do                 !###### read from the file each time
       return
      end
C     =================================================
      Subroutine CFnc(kq,pa,qa,pb,qb) !### - function Yk
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ii/Ii/Kt/Kt/Kt2/Kt2
     >        /Kk/Kk(IPs)
     >        /C/C(IP6)
       dimension pa(IP6),qa(IP6),pb(IP6),qb(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Kt=Kt2
        ih=2-Kt
        do i=1,Ii,ih
           C(i)=pa(i)*pb(i)+qa(i)*qb(i)
        end do
        C(Ii+4)=pa(Ii+4)+pb(Ii+4)
        call Yk(kq)
       return
      end
C     =================================================
      Function Rint(pa,qa,pb,qb)  !#  radial integral for the
      implicit real*8 (a-h,o-z)   !##  Coulomb interaction
      INCLUDE "hfd.par"
       common /Ii/Ii/Kt/Kt
     >        /C/C(IP6)/R/R(IP6)
       dimension pa(IP6),qa(IP6),pb(IP6),qb(IP6)
        ih=2-Kt
        do i=1,Ii,ih
           C(i)=C(i)/R(i)*(pa(i)*pb(i)+qa(i)*qb(i))
        end do
        C(Ii+4)=C(Ii+4)+pa(Ii+4)+pb(Ii+4)-1
        call Sint(s)
        Rint=s
       return
      end
C     =================================================
      Function Cwig(j1,j2,j3,j4,jq) !# Angular factor for the RPA matrices.
       implicit real*8 (a-h,o-z)    !## Note that parameters ji are INTEGER!
       h=0.5d0
       Cwig = dsqrt((j1+1)*(j2+1)*(j3+1)*(j4+1)*1.d0)
     >      * Wig0(jq,j1*h,j2*h) * Wig0(jq,j3*h,j4*h)
       return
      end
C     =================================================
      Subroutine Iteration      !#  Iterative solution of the system of two
      INCLUDE "rpa.par"         !##  coupled RPA equations. If omega=0, first
       implicit real*8 (a-h,o-z)!###  of the two equations is trivial, and
                                !####  one iteration solves the problem.
       common /Nr/Nr/Kmg/Kmg/Kl/Kl/Omega/Omega
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0/F0(IPv)
     >        /Fp0/Fp0(IPv)   /Fm0/Fm0(IPv)
     >        /Dp/Dp(IPr)     /Dm/Dm(IPr)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        acc=0.0002       !# - convergence parameter
        mx=100           !# - max number of iterations
        err1=0.d0
        write(*, 5)
        write(11,5)
 5      format (1X,15('#'),
     >       ' Iterative solution of the RPA system ',15('#'))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,Nr
           Fp0(i)=F0(i)
           Fm0(i)=0.d0
           Fm(i)=0.d0
           Fp(i)=0.d0
        end do
        do iter=1,mx
           write(*, 15) iter
           write(11,15) iter
 15        format (1X,'iteration',I3)
           if (kmg.EQ.1) call SolvEq(1,kmg,err1)
           call SolvEq(2,kmg,err2)
           do i=1,Nr
              Fp0(i)=Fp(i)
              Fm0(i)=Fm(i)
           end do
           if (err1+err2.LE.acc.OR.kmg.EQ.0) return
        end do
       return
      end
C     =================================================
      Subroutine SolvEq(keq,kmg,cng) !#  Solution of the RPA equations.
      INCLUDE "rpa.par"              !##  Keys: keq = 1,2 correspond to the
       implicit real*8 (a-h,o-z)     !###  first/second subsystem of coupled
       common /Nr/Nr                 !####  equations. kmg = 0,1 corresponds
                                     !#####  to omega=0 or omega <> 0.
     >        /Tp/Tp(IPr,IPr) /Tm/Tm(IPr,IPr) /T/T(IPr*IPr)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0/F0(IPv)
     >        /Fp0/Fp0(IPv)   /Fm0/Fm0(IPv)
     >        /Dp/Dp(IPr)     /Dm/Dm(IPr)
       dimension Idec(IPr),Scales(IPr),X(IPr),Y(IPr)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     formation of the matrix:
        do i=1,Nr
           do k=1,Nr
              if (keq.EQ.1) then
                 sl=-Tm(i,k)*Dp(k)
              else
                 sl=-Tp(i,k)*Dp(k)
              end if
              if (i.EQ.k) sl=1.d0+sl
              T(Nr*(k-1)+i)=sl
           end do
        end do
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     formation of the RHS:
        do i=1,Nr
           if (keq.EQ.1) then
              sr=0.d0
           else
              sr=F0(i)
           end if
           if (kmg.EQ.1) then                  !# Fp0 and Fm0 are known
              do k=1,Nr                        !## from the previous
                 if (keq.EQ.1) then            !### iteration
                    sr=sr+Tm(i,k)*Dm(k)*Fp0(k)
                 else
                    sr=sr+Tp(i,k)*Dm(k)*Fm0(k)
                 end if
              end do
           end if
           Y(i)=sr
        end do
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     solution of the equation:
        call Decomp(IPr*IPr,T,Nr,Scales,Idec)
        call Flsolv(IPr*IPr,Nr,T,Y,X,Idec)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     check of the correctness:
        do i=1,Nr
           do k=1,Nr
              if (keq.EQ.1) then
                 sl=-Tm(i,k)*Dp(k)
              else
                 sl=-Tp(i,k)*Dp(k)
              end if
              if (i.EQ.k) sl=1.d0+sl
              T(Nr*(k-1)+i)=sl
           end do
        end do
        s1=0.d0
        s2=0.d0
        do i=1,Nr
           s3=0.d0
           do k=1,Nr
              s3=s3+T(Nr*(k-1)+i)*X(k)
           end do
           s1=s1+(s3-Y(i))**2
           s2=s2+Y(i)**2
        end do
        err=dsqrt(s1/s2)
        if (err.GT.1.d-6) then
           write(*, 5) dsqrt(s2),err
           write(11,5) dsqrt(s2),err
 5         format(3X,'Correctness test: |Y| =',E12.5,
     >          ' ERROR = |TX-Y|/|Y| =',F9.6,' push..')
           read (*,*)
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (keq.EQ.1) then
           do i=1,Nr
              Fm(i)=X(i)
              Y(i)=Fm0(i)
           end do
        else
           do i=1,Nr
              Fp(i)=X(i)
              Y(i)=Fp0(i)
           end do
        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     check of the convergence:
        s1=0.d0
        s2=0.d0
        s3=0.d0
        do i=1,Nr
           s1=s1+Y(i)*X(i)
           s2=s2+X(i)**2
           s3=s3+(X(i)-Y(i))**2
        end do
        cng=dsqrt(s3/s2)
        write(*, 15) keq,s1/s2,cng
        write(11,15) keq,s1/s2,cng
 15     format(3X,'Eq.',I1,': <X_old|X_new>/|X_new|**2 =',F8.4,
     >       ' |X_new-X_old|/|X_new|=',F8.4)
       return
      end
C     =================================================
      Subroutine Check
c#    Check if the vector (Fm,Fp) is the solution
c#    of the RPA equations:
c#    |A11 A12| |Fm| = | 0|
c#    |A21 A22| |Fp|   |F0|
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Nr/Nr
     >        /Tp/Tp(IPr,IPr) /Tm/Tm(IPr,IPr)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0/F0(IPv)
     >        /Dp/Dp(IPr)     /Dm/Dm(IPr)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        sm1=0.d0
        sm2=0.d0
        x=0.d0
        do i=1,Nr
           sum1=Fm(i)
           sum2=Fp(i)
           do k=1,Nr
              sum1=sum1-Tm(i,k)*(Dp(k)*Fm(k)+Dm(k)*Fp(k))
              sum2=sum2-Tp(i,k)*(Dp(k)*Fp(k)+Dm(k)*Fm(k))
           end do
           sm1=sm1 + sum1**2
           sm2=sm2 + (sum2-F0(i))**2
           x=x + F0(i)**2
        end do
        err=dsqrt((sm1+sm2)/x)
        write(*, 5) dsqrt(sm1),dsqrt(sm2),dsqrt(x),err
        write(11,5) dsqrt(sm1),dsqrt(sm2),dsqrt(x),err
 5      format(1X,50('-'),
     >       /3X,'|A_11*Fm + A_12*Fp|       =',E10.3,
     >       /3X,'|A_21*Fm + A_22*Fp - F0|  =',E10.3,
     >       /3X,'|F0|                      =',E10.3,
     >       /3X,'|A(Fm,Fp)`- (0,F0)`|/|F0| =',E10.3)
       return
      end
C     =================================================
      Subroutine Prin !### - Radial MEs of the RPA operator
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nr/Nr/Kl/Kl/Omega/Omega/Ita/Ita
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /F0/F0(IPv)/Fp/Fp(IPv)/Fm/Fm(IPv)
     >        /Alet/Alet(IPcs)
       character*1 let(9)
       character*4 Alet
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        data let/'s','p','d','f','g','h','i','k','l'/
        write(*, 5) Alet(Kl),Omega
        write(11,5) Alet(Kl),Omega
 5      format(1X,20('#'),' RPA for ',A4,'(w =',F9.6,') ',
     >       20('#'),/4X,'i',4X,'init',5X,'final',
     >       7X,'F0',11X,'F^+',10X,'F')
        s1=0.d0
        s2=0.d0
        s3=0.d0
        do i=1,Nr
           ia=Iri(i)
           na=Nn(ia)
           la=Ll(ia)
           ja=Lj(ia)
           ib=Irk(i)
           nb=Nn(ib)
           lb=Ll(ib)
           jb=Lj(ib)
           s1=s1+F0(i)**2
           s2=s2+(Fp(i)+Fm(i)-F0(i))**2
           s3=s3+(Fp(i)-Fm(i)-F0(i))**2
           write(11,15) i,na,let(la+1),ja,nb,let(lb+1),jb,
     >          F0(i),(Fp(i)+Fm(i)),(Fp(i)-Fm(i))
 15        format(2I5,A1,I1,'/2',I5,A1,I1,'/2',3E14.5)
        end do
        s2=dsqrt(s2/s1)
        s3=dsqrt(s3/s1)
        write( *, 25) s2,s3
        write( 11,25) s2,s3
 25     format(1X,63('-'),/3X,'|F^+ -F0|/|F0| =',F8.3,
     >       3X,'|F-F0|/|F0| =',F8.3)
       Return
      end
C     =================================================
      Subroutine Vinit !### - Definition of the valence space
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ns/Ns/Nss/Nss/Nsh/Nsh/Nv/Nv/Lmax/Lmax
     1        /Ip/Ip/Ij/Ij/Kl/Kl/Ita/Ita/Nsv/Nsv/Nmax/Nmax
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Ivi/Ivi(IPv)/Ivk/Ivk(IPv)
       dimension p(IP6),q(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        call ReadF (13,1,p,q,2)
        ns1=p(2)+0.1d0
        Nv=0                       !# - Dimension of the valence space
        Nmax=min(Nmax,Nss,ns1)
        do i=Nsv,Nmax
           li=Ll(i)
           ji=Lj(i)
           num=0
           do k=i,Nmax
              lk=Ll(k)
              jk=Lj(k)
              ipp=Isig(li)*Isig(lk)*ip
              jm=iabs(ji-jk)
              js=ji+jk
              if (ipp.EQ.1.AND.jm.LE.ij.AND.js.GE.ij
     >            .AND.max(li,lk).LE.Lmax) then
                 Nv=Nv+1
                 num=num+1
                 Ivi(Nv)=i
                 Ivk(Nv)=k
                 if (Nv.GT.IPv) then
                    write(*,*) 'Nv=',Nv,' > IPv=',IPv
                    write(*,*) 'current shell is',i
                    stop
                 end if
              end if
           end do
        end do
        write(*, 15) Nv,Lmax
        write(11,15) Nv,Lmax
 15     format(4X,'Valence space dimension =',I4,
     >        ' (max(li,lk) <=',I2,')')
       Return
      end
C     =================================================
      Subroutine RHSv      !# - Reduced MEs for the RHS
                           !## operator in the valence subspace
      INCLUDE "rpa.par"
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Nv/Nv/Kl/Kl/Ii/Ii
     1        /Z/Z/Ita/Ita/Kt/Kt/Kt1/Kt1
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Jj/Jj(IPs)
     >        /Ivi/Ivi(IPv)/Ivk/Ivk(IPv)
     >        /F0v/F0v(IPv)
     >        /UP/Pa(IP6),Pc(IP6),Pd(IP6),Pb(IP6)
     >        /DOWN/Qa(IP6),Qc(IP6),Qd(IP6),Qb(IP6)
        Kt=Kt1
        ia0=0
        do i=1,Nv
           ia=Ivi(i)
           na=Nn(ia)
           la=Ll(ia)
           ja=Lj(ia)
           ka=Kk(ia)
           if (ia0.NE.ia) then
              call ReadF(13,ia+4,pa,qa,2)
              ia0=ia
           end if
           ib=Ivk(i)
           nb=Nn(ib)
           lb=Ll(ib)
           jb=Lj(ib)
           kb=Kk(ib)
           call ReadF(13,ib+4,pb,qb,2)
           F0v(i)=Fnc(la,ja,ka,na,lb,jb,kb,nb,pa,qa,pb,qb)
        end do
       Return
      end
C     =================================================
      Subroutine CoreCore        !# Calculation of the effective reduced
      INCLUDE "rpa.par"          !## Core-Core MEs. They are used in SolVal
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Nr/Nr/Nr1/Nr1/Kl/Kl/Nsh/Nsh
     >        /Omega/Omega/Ita/Ita/Kbrt/Kbrt
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /Fp0/Fp0(IPv)  /Fm0/Fm0(IPv)  /F0/F0(IPv)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)
     >        /UP/P1(IP6),A1(IP6),P2(IP6),A2(IP6)
     >        /DOWN/Q1(IP6),B1(IP6),Q2(IP6),B2(IP6)
     >        /Alet/Alet(IPcs)
       character*1 let(9)
       character*4 Alet
       character*13 clmb_br(3)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        data let/'s','p','d','f','g','h','i','k','l'/
        clmb_br(1)='Pure Coulomb '
        clmb_br(2)='Coulomb-Breit'
        clmb_br(3)='Coulomb-Breit'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nr1.LE.Nr) return
        write(*, 5) Alet(Kl),Omega,clmb_br(Kbrt+1)
        write(11,5) Alet(Kl),Omega,clmb_br(Kbrt+1)
 5      format(1X,61('#'),/1X,4('#'),' Effective Core-Core integrals ',
     >       'for ',A4,'(w =',F9.6,') ',3('#'),/1X,11('#'),2X,A13,
     >       ' approximation is used',2X,11('#'),
     >       /4X,'i',6X,'n',8X,'k',9X,'F0',11X,'F^+',10X,'F')
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     formation of the Coulomb (Coulomb-Breit) matrices Tp and Tm
c     and solution of the equations for Fp and Fm:
        write(*, 15) Nr*(Nr1-Nr)
 15     format(25X,'Number of Coulomb MEs <=',I8)
        s1=0.d0
        s2=0.d0
        s3=0.d0
        do i=Nr+1,Nr1
           ia=Iri(i)
           na=Nn(ia)
           la=Ll(ia)
           ja=Lj(ia)
           ib=Irk(i)
           nb=Nn(ib)
           lb=Ll(ib)
           jb=Lj(ib)
           sp=F0(i)
           sm=0.d0
           call ReadFF (12,ia+4,P1,Q1,2)
           call ReadFF (12,ib+4,P2,Q2,2)
           do m=1,Nr
              mi=Iri(m)
              mk=Irk(m)
              call ReadFF (12,mi+4,A1,B1,2)
              call ReadFF (12,mk+4,A2,B2,2)
              call Clmb(ia,ib,mi,mk,t1,t2,t3,ff)
              if (Kbrt.GE.1) call Brt(ia,ib,mi,mk,t1,t2,t3)
              u1=(t1+t2)*ff
              u2=(t1+t3)*ff
              s=Eps(mi)-Eps(mk)
              Dp=s/(s*s-Omega*Omega)
              Dm=Omega/(s*s-Omega*Omega)
              sp = sp + (u1+Ita*u2)*Dp*Fp0(m) + (u1+Ita*u2)*Dm*Fm0(m)
              sm = sm + (u1-Ita*u2)*Dp*Fm0(m) + (u1-Ita*u2)*Dm*Fp0(m)
           end do
           Fp0(i)=sp
           Fm0(i)=sm
           s1=s1+F0(i)**2
           s2=s2+(sp-sm-F0(i))**2
           s3=s3+(sp+sm-F0(i))**2
           write ( *,25) i,na,let(la+1),ja,nb,let(lb+1),jb,
     >                   F0(i),sp+sm,sp-sm
           write (11,25) i,na,let(la+1),ja,nb,let(lb+1),jb,
     >                   F0(i),sp+sm,sp-sm
 25        format (2I5,A1,I1,'/2',I5,A1,I1,'/2',3E14.5)
        end do
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        s2=dsqrt(s2/s1)
        s3=dsqrt(s3/s1)
        write( *,45) s2,s3
        write(11,45) s2,s3
 45     format(1X,78('-'),/3X,'|F^+ -F0|/|F0| =',E10.3,
     >       3X,'|F -F0|/|F0| =',E10.3)
       return
      end
C     =================================================
      Subroutine SolVal            !# Calculation of the effective reduced
      INCLUDE "rpa.par"            !## MEs for valence space
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Nr/Nr/Nv/Nv/Kl/Kl/Nsv/Nsv/Nsh/Nsh/k_denom/k_denom
     >        /Omega/Omega/Ita/Ita/Kbrt/Kbrt/Small/Small
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Iri/Iri(IPv)/Irk/Irk(IPv)
     >        /Ivi/Ivi(IPv)/Ivk/Ivk(IPv)
     >        /Fp/Fp(IPv)    /Fm/Fm(IPv)   /F0v/F0v(IPv)
     >        /Fp0/Fp0(IPv)  /Fm0/Fm0(IPv)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /UP/P1(IP6),A1(IP6),P2(IP6),A2(IP6)
     >        /DOWN/Q1(IP6),B1(IP6),Q2(IP6),B2(IP6)
     >        /Alet/Alet(IPcs)
       character*1 let(9)
       character*4 Alet
       character*13 clmb_br(3)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        data let/'s','p','d','f','g','h','i','k','l'/
        clmb_br(1)='Pure Coulomb '
        clmb_br(2)='Coulomb-Breit'
        clmb_br(3)='Coulomb-Breit'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*, 5) Alet(Kl),Omega,clmb_br(Kbrt+1)
        write(11,5) Alet(Kl),Omega,clmb_br(Kbrt+1)
 5      format(1X,61('#'),/1X,5('#'),' Effective radial integrals ',
     >       'for ',A4,'(w =',F9.6,') ',5('#'),/1X,11('#'),2X,A13,
     >       ' approximation is used',2X,11('#'),/1X,61('#'))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     formation of the Coulomb (Coulomb-Breit) matrices Tp and Tm
c     and solution of the equations for Fp and Fm:
        write(11,15)
 15     format(4X,'i',6X,'n',8X,'k',9X,'F0',11X,'F^+',10X,'F')
        write(*, 25) Nv*Nr
 25     format(25X,'Number of Coulomb MEs <=',I8)

c        write(*,*)' SolVal: there is unresolved issue with denominators'
c        write(*,*)' You have to look at the code again before using it'
c        write(*,*)' The two possible variants are marked below'

        idel=IPdl
        nu=0
        nu1=0
        s1=0.d0
        s2=0.d0
        s3=0.d0
        do i=1,Nv
          ia=Ivi(i)
          na=Nn(ia)
          la=Ll(ia)
          ja=Lj(ia)
          ib=Ivk(i)
          nb=Nn(ib)
          lb=Ll(ib)
          jb=Lj(ib)
          sp=F0v(i)
          sm=0.d0
          call ReadF (13,ia+4,P1,Q1,2)
          call ReadF (13,ib+4,P2,Q2,2)
          do m=1,Nr
            mi=Iri(m)
            mk=Irk(m)
            if (mi.LT.Nsv) then                 !# note that valence MEs
              call ReadFF (12,mi+4,A1,B1,2)     !## do not depend on
              call ReadFF (12,mk+4,A2,B2,2)     !### shells Nsv - Nsh
              call Clmb(ia,ib,mi,mk,t1,t2,t3,ff)
              if (Kbrt.GE.1) call Brt(ia,ib,mi,mk,t1,t2,t3)
              u1=(t1+t2)*ff
              u2=(t1+t3)*ff
              s=Eps(mi)-Eps(mk)

c>>>>>>>>>>>>>> Two variants for denominators >>>>>>>>>>>>>>>>>>>>>>>>>
c Variant of denominators that provides gauge invariance for E1 (22/02/01)
              if (k_denom.EQ.1) then
                Dp=s/(s*s-Omega*Omega)
                Dm=Omega/(s*s-Omega*Omega)
                sp = sp + (u1+Ita*u2)*Dp*Fp0(m)
     >                  + (u1+Ita*u2)*Dm*Fm0(m)
                sm = sm + (u1-Ita*u2)*Dp*Fm0(m)
     >                  + (u1-Ita*u2)*Dm*Fp0(m)
C Previous variant that allows to change valence energies
              else
                emin=dmin1(Eps1(ia),Eps1(ib))
                if (ia.GE.Nsv) then               !# denominators in P space
                  Da=1.d0/(s+emin-Eps2(ia))       !## correspond to MBPT
                  Db=1.d0/(s+emin-Eps2(ib))       !### rules, while in Q
                else                              !#### space --- to RPA
                  Da=1.d0/s                       !##### rules (with zero
                  Db=1.d0/s                       !###### frequency).
                end if
                sp = sp + (u1*Da + Ita*u2*Db)*Fp0(m)
                sm = sm + (u1*Da - Ita*u2*Db)*Fm0(m)
              end if
c <<<<<<<<<<<<< end of variants of denominator <<<<<<<<<<<<<<<<<<<<<<

              nu=nu+1
              if (dabs(t1)+dabs(t2)+dabs(t3).NE.0.d0) then
                nu1=nu1+1
              end if
              if (nu-idel*(nu/idel).EQ.0) then
                write(*,35) nu
 35             Format('+',3X,'nu =',I8)
              end if
            end if
          end do
          Fp(i)=sp
          Fm(i)=sm
          s1=s1+F0v(i)**2
          s2=s2+(sp-sm-F0v(i))**2
          s3=s3+(sp+sm-F0v(i))**2
          write(11,45) i,na,let(la+1),ja,nb,let(lb+1),jb,
     >                 F0v(i),sp+sm,sp-sm
 45       format(2I5,A1,I1,'/2',I5,A1,I1,'/2',3E14.5)
        end do
        write(*, 55) nu,nu-nu1
 55     format(5X,I8,' Coulomb MEs ( ', I6,' zeros) calculated.')
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        s2=dsqrt(s2/s1)
        s3=dsqrt(s3/s1)
        write( *,65) s2,s3,k_denom
        write(11,65) s2,s3,k_denom
 65     format(1X,68('-'),/3X,'|F^+ -F0|/|F0| =',E10.3,
     >       3X,'|F-F0|/|F0| =',E10.3,
     >       /3X,'Denominators correspond to var. ',I1)
        if (Nsv.LE.Nsh) then
           write( *,75) Nsv,Nsh
           write(11,75) Nsv,Nsh
 75        format(3X,'Note: valence MEs do not depend on core shells ',
     >        I3,'...,',I3,/)
        end if
       return
      end
C     =================================================
      Subroutine Rsave(kv)  !#  Formation of file RPA_n.INT, where n=Kl
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsh/Nsh/Nhf/Nhf/Nr/Nr/Nv/Nv/Kl/Kl
     1        /Nsv/Nsv/Nmax/Nmax/Lmax/Lmax
       common /Iri/Iri(IPv)   /Irk/Irk(IPv)
     >        /Ivi/Ivi(IPv)   /Ivk/Ivk(IPv)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0v/F0v(IPv)
     >        /Fp0/Fp0(IPv)   /Fm0/Fm0(IPv)   /F0/F0(IPv)
     >        /Eps1/Eps1(IPs)
       character*1 cn(9),cf(9)
       character*9 FNAME
       equivalence (FNAME,cf(1))
        data cn /'1','2','3','4','5','6','7','8','9'/
        data cf /'R','P','A','_','n','.','I','N','T'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
!old        cf(5)=cn(Kl)
        if (Kl.LT.10) then
          cf(4)='_'
          cf(5)=cn(Kl)
        else
          cf(4)='1'
          if (Kl.EQ.10) then
            cf(5)='0'
          else
            cf(5)=cn(Kl-10)
          end if
        end if
        if (kv.NE.2) Nv=0
        write(*,5) FNAME
 5      format(5X,' Forming file ',A9)
        open(15,file=FNAME,status='UNKNOWN',form='UNFORMATTED')
        write(15) Ns,Nsh,Nhf,Nr,Nv,Nsv,Nmax,Lmax,(Eps1(i),i=Nsv,Nmax)
        write(15) (Iri(i),Irk(i),i=1,Nr)
        write(15) (F0(i),Fp0(i),Fm0(i),i=1,Nr)
        if (Nv.GE.1) then
          write(15) (Ivi(i),Ivk(i),i=1,Nv)
          write(15) (F0v(i),Fp(i),Fm(i),i=1,Nv)
          write(15) 0            !### number of subtraction integrals
        end if
        close(15)
        write(*,15) FNAME
 15     format(5X,' File ',A9,' formed')
       return
      end
C     =================================================
      Subroutine Rinput(ierr)
c#    Input from the file RPA_n.INT, where n=Kl
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsh/Nsh/Nhf/Nhf/Nr/Nr/Nr1/Nr1/Nv/Nv/Kl/Kl
       common /Iri/Iri(IPv)   /Irk/Irk(IPv)
     >        /Ivi/Ivi(IPv)   /Ivk/Ivk(IPv)
     >        /Fp0/Fp0(IPv)   /Fm0/Fm0(IPv)   /F0/F0(IPv)
       dimension ff0(IPv)
       character*1 cn(9),cf(9)
       character*9 FNAME
       equivalence (FNAME,cf(1))
        data cn /'1','2','3','4','5','6','7','8','9'/
        data cf /'R','P','A','_','n','.','I','N','T'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        err=1.d-6
        if (Kl.LT.10) then
          cf(4)='_'
          cf(5)=cn(Kl)
        else
          cf(4)='1'
          if (Kl.EQ.10) then
            cf(5)='0'
          else
            cf(5)=cn(Kl-10)
          end if
        end if

        ierr=1
        open(15,file=FNAME,status='OLD',form='UNFORMATTED',err=1000)
        ierr=0
        write(*,5) FNAME
 5      format(5X,' Reading file ',A9)
        read(15) ns0,nsh0,nhf0,nr0
        ierr=iabs(Ns-ns0)+iabs(Nsh-nsh0)+iabs(Nhf-nhf0)
        ierr=ierr+min(iabs(Nr1-nr0),iabs(Nr-nr0))
        if (ierr.EQ.0) then
          read(15) (Iri(i),Irk(i),i=1,Nr)
          read(15) (ff0(i),Fp0(i),Fm0(i),i=1,Nr)
          s=0.d0
          s1=0.d0
          do i=1,Nr
            s = s + dabs(F0(i)-ff0(i))
            s1 = s1 + dabs(F0(i))
          end do
          s=s/s1
          if (s.GT.err) then
            ierr=1
            write( *,15) s
            write(11,15) s
 15         format(5X,' RHS differs from that in the file: s =',E10.3)
          end if
        end if
        close(15)
        if (ierr.EQ.0) then
          write( *,25) FNAME
          write(11,25) FNAME
 25       format(4X,'Solution of RPA equations is taken from file ',A9)
        end if
 1000  return
      end
C     =================================================
      subroutine NclInt(n1,n2,s1,s2,p,q,a,b,dn)
c radial integration inside the nucleus: R(1)**n1*int_0,1 f(x) dx
c where x = r/R(1)
c and   f(x) = (s1*p(x)*q(x) + s2*a(x)*b(x))*x**n2
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ii/Ii/Jmax/Jmax/R/R(IP6)
       dimension p(IP6),q(IP6),a(IP6),b(IP6)
        dn=0.d0
        do i=1,Jmax
           i1=Ii+4+i
           do j=1,Jmax
              j1=Ii+4+j
              m2=i+j+n2-1
              dn = dn + (s1*p(i1)*q(j1) + s2*a(i1)*b(j1))/m2
           end do
        end do
        dn=dn*R(1)**n1
       return
      end
