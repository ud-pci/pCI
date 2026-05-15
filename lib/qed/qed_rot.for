C     =================================================
      program QED_ROT         !### last update 30/05/15
      include "conf.par"
      include "hfd.par"
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c       Based on Bas_B
c       Rotates orbitals to diagonalize the one electron
c       Hamiltonian with QED corrections.
c       Two step diagonalization of the basis set is
c       made to allow higher accuracy.
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Nsv/Nsv/Kt/Kt/Kt1/Kt1/Sf/Sf
     1        /MaxT/MaxT/Ierr/Ierr/Merr/Merr/Ii/Ii/Nc/Nc/Ne/Ne
     3        /Nhint/Nhint/Nsint/Nsint/FNAME/FNAME/C_is/C_is
     4        /Ksg/Ksg/Kdg/Kdg/Kkin/Kkin/Kbrt/Kbrt/Khf/Khf/K_is/K_is
     5        /Kfst/n3,j3,ch3/Vfst/n4,j4,ch4/Blst/n5,j5,ch5
       common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /let/let(9)/Nf0/Nf0(IPs)/Qq/Qq(IPs)
     >        /Nq/Nq(IPsp)/Nip/Nip(IPsp)
     >        /QNL/QNL(IPsp)
     >        /Ndc/Ndc(IPc)/Nc0/Nc0(IPc)/Nvc/Nvc(IPc)
     >        /rint1/rint1(10*IPh),iint1(10*IPh)
     >        /rint2/rint2(10*IPh),iint2(10*IPh)
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)/Y/Y(IP6)
     >        /P/P(IP6)/Q/Q(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
       common /kout/kout/small/small
       common /ipmr/ipmr
       character*1 let,ch3,ch4,ch5,FNAME*12
        data let /'s','p','d','f','g','h','i','k','l'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        call recunit
        MaxT=9           !### length of expansion at the origin
        kout=1           !### output details
        irec=2*IP6*IPmr  !### record length in DAT files
        iter=1
        ier=0
        small=1.d-8
        Kt1=0            !### used as Kt for Breit integrals
        Kbrt=0           !### 0 - Coulomb, 1 - Gaunt, 2 - Full Breit
        K_is=0           !### 0 - no IS, 1 - Volume shift, 2 - SMS
        C_is=0.d0        !### scaling factor for IS
        Nsint=0          !### number of MEs of Sigma
        Khf=0            !### =1 if self-consistency is reached
        FNAME='QED_ROT.DAT'

 1      Ierr=0           !### number of or warnings
        open(unit=11,file='QED_ROT.RES',status='UNKNOWN')
        call Init0
        call Input

        write( *,5) iter,ier
        write(11,5) iter,ier
 5      format(1X,70('=')/4X,'Iteration ',I3,' (had ',I3,' warnings)',
     >        /1X,70('='))

C        Nsv=N_orb(n4,ch4,j4) !### last "frozen" orbital
        Nsv=0
        if (Khf.EQ.1) Nsv=max(Nsv,Nso)
        n1=Nsv+1
        if (Nsv.GT.Nso) then      !### dumping factor for i_diag
          Sf=0.75d0
        else
          Sf=1.00d0
        end if
        write( *,*) ' Last frosen orbital is',Nsv
        write(11,*) ' Last frosen orbital is',Nsv

        open(12,file=FNAME,status='OLD',access='DIRECT',recl=irec)
        open(13,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)

        if (Ns.GE.n1) then

          call V_core                         !### Y(r)=V_core-Z/r

C          nmin=N_orb(n3,ch3,j3)               !### first orbital to apply
C          if (nmin.EQ.0) nmin=Ns+1            !#### kinetic balance condition
          nmin=Ns+1
          write( *,*) ' Kinetic ballance is not used'
          write(11,*) ' Kinetic ballance is not used'

          do ni=n1,Ns
            if (ni.GE.nmin) call Change_Q(ni) !### kinetic balance for Q
            call Ort(ni,Ns,iort)
            if (iort.GE.1) then
              if (ni.GE.nmin) then
                write( *,15) ni,iort
                write(11,15) ni,iort
 15             format(4X,'Big change in Q for orbital',I4,' iort =',I4)
              else
                write( *,25) ni,iort
                write(11,25) ni,iort
 25             format(4X,'No orthogonality for',I4,' iort =',I4)
              end if
              read(*,*)
            end if
          end do
        end if

        do ni=1,2*Ns+4
          call ReadF (12,ni,P,Q,2)
          call WriteF(13,ni,P,Q,2)
        end do
        close(13)                  !### sub-ne Core changes records
        close(12)                  !#### Ns+5,...2Ns+4 in the file FNAME

        call Core                  !### calculates Dirac-Fock operator
        Call Rint                  !### Radial integrals depend on Ksg,Kbrt
        idg=Kdg
        Call Rot(idg)              !### idg=0,1,2 - diagonalization status

        nr=0
        if (idg.NE.0) then
          do ni=n1,Ns
            call Ort(ni,Ns,iort)
            if (iort.GE.1) then
              write( *,35) ni,iort
              write(11,35) ni,iort
 35           format(4X,I3,': orthogonatily is lost. iort=',I4)
              nr=nr+1
            end if
          end do
        end if
        close (12)
        close (13,status='DELETE')
        if (nr.GT.0) then
          read(*,*)
          stop
        end if

        if (idg.NE.Kdg) then
          write( *,45) idg
          write(11,45) idg
 45       format(/' Diagonalization key was changed to',I2)
          Ierr=Ierr+100
        end if

        nlst=Ns
        write( *,*) ' All orbitals are kept in the basis set'
        write(11,*) ' All orbitals are kept in the basis set'
        if (nlst.GT.0.AND.nlst.LT.Ns) then
          open(12,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)
          open(13,file=FNAME,status='NEW',access='DIRECT',recl=irec)
          call ReadF (12,1,P,Q,2)
          P(2)=nlst
          call WriteF (13,1,P,Q,2)
          do n=2,nlst+4
            call ReadF (12,n,P,Q,2)
            call WriteF (13,n,P,Q,2)
          end do
          do n=1,nlst
            call ReadF (12,n+Ns+4,P,Q,2)
            call WriteF (13,n+nlst+4,P,Q,2)
          end do
          close(12)
          close(13)
          write( *,65) nlst,FNAME
          write(11,65) nlst,FNAME
 65       format(/' First ',I3,' orbitals saved in ',A12)
        end if

        if (Ierr.GT.0) then
          write( *,75) Ierr
          write(11,75) Ierr
 75       format(/' There were',I4,' warnings.')
        end if

        if (Nsint.GT.0.AND.idg.GT.0) then
          write( *,85)
          write(11,85)
 85       format(' You may need to recalculate Sigma now.')
        end if

        write( *,95) Nso,Kdg,Kbrt,Ksg
        write(11,95) Nso,Kdg,Kbrt,Ksg
 95     format(/' *****   Summary   *****',
     >         /' SCF field for',i3,' closed shells.',
     >         /' Diagonalization:',i2,' Breit:',i2,' QED:',i2,
     >         /' *****   *******   *****')
        
       stop
      end
c     =================================================
      include "bas1.inc"
      include "breit.inc"
      include "hould.inc"
      include "i_diag.inc"
C    include "inpstr.inc"
      include "yk.inc"
      include "ykt.inc"
      include "sint1.inc"
      include "readf.inc"
      include "dif.inc"
      include "wig.inc"
      include "test_ori.inc"
      include "rec_unit.inc"
C     =================================================
      subroutine Init0
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "phys.par"
      include "hfd.par"
       common /Ns/Ns/Nso/Nso/Ng/Ng
     1        /II/II/Kt/Kt/Z/Z/Cl/Cl
     2        /H/H/Rnuc/Rnuc/FNAME/FNAME
       common /Qnl/Qnl(IPsp)/Nq/Nq(IPsp)/Nip/Nip(IPsp)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Qq/Qq(IPs)
     >        /R/R(IP6)/V/V(IP6)
       common /ipmr/ipmr/AL/AL/BT/BT/R2/R2
       dimension P(IP6),Q(IP6),P1(IP6),Q1(IP6),PQ(4*IP6),Qw(IPs)
       real*8 JM
       character*12 FNAME
       logical longbasis
       dimension IQN(4*IPs),Qq1(IPs)
       equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
       equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),
     >        (P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c small number:
        c1=0.01d0
c speed of light is taken from "phys.par":
        Cl=DPcl
        Mj=2*dabs(Jm)+0.01d0
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
        R2  =dabs(PQ(5))
        H   =PQ(6)
        Bt  =PQ(7)
        Al  =PQ(8)
        Kt  =PQ(9)+C1
        Ng  =PQ(10)+C1
        Rnuc=PQ(13)
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          write(11,*) ' Using variant for long basis '
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
            Qw(ni)=0.d0
          end do
        end if

        Nso=0
        do ni=1,Ns
          xq=Jj(ni)+1.d0
          if (xq.GT.QQ(ni)+c1) then
            goto 210
          else 
            Nso=Nso+1
          end if
        end do

 210    write( *,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Ns,Nso
        write(11,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Ns,Nso
 5      format (/4X,'Program QED_ROT',
     >          //4X,'Z   = ',F6.2,5X,'Kt  =',I3,  7X,'II =',I4,
     >          /4X,'H   =',F7.4, 5X,'R2  =',F6.2,4X,'R1 =',E11.4,
     >          /4X,'Rnuc=',E11.4,1X,'Al  =',F7.4,3X,'Bt =',F5.2,
     >          /4X,'Ns  =',I4,   8X,'Nso =',I3)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(13,file=FNAME,status='UNKNOWN',
     >       access='DIRECT',recl=2*IP6*IPmr)
        do ni=1,4
          call ReadF (12,ni,P,Q,2)
          call WriteF(13,ni,P,Q,2)
        end do
        do ni=1,Ns
          call ReadF (12,ni+4,P,Q,2)
          call WriteF(13,ni+4,P,Q,2)
          call ReadF (12,ni+4+Ns,P,Q,2)
          call WriteF(13,ni+4+Ns,P,Q,2)
        end do
        close(13)
        close(12)
       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
        write(11,75)
75      format(/2X,'file HFD.DAT is absent'/)
       STOP
      end
C     =================================================
      subroutine Input
      include "conf.par"
       implicit real*8 (a-h,o-z)
       common /Nso/Nso/Nsv/Nsv/Z/Z/let/let(9)
     >        /Kkin/Kkin/Ksg/Ksg/Kdg/Kdg/Merr/Merr/Kbrt/Kbrt
       common /Qnl/Qnl(IPsp)/Qnl1/Qnl1(IPsp)
     >        /Dfst/n1,j1,ch1
     >        /Kfst/n3,j3,ch3/Vfst/n4,j4,ch4/Blst/n5,j5,ch5
       character*1 name(16),Let,str1(4)*8,str2(2)*9,
     >           ch1,ch3,ch4,ch5
       data str1,str2 /'  DF-C  ','DF-C-QED','  DF-B  ','DF-B-QED',
     >      ' General ','1st order'/
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,*)' Diagonalization: (1) general, (2) First order?'
        read (*,*) Kdg
        if (Kdg.LE.0) stop
        write(*,*)' No QED (1), or with QED (2)?'
        read (*,*) Ksg
        write(*,*)' Coulomb (0), Gaunt (1), or full Breit (2)?'
        read (*,*) Kbrt

        i1=(Ksg-1)*(Ksg-2)
        i2=Kdg*(Kdg-1)*(Kdg-2)
        i3=Kbrt*(Kbrt-1)*(Kbrt-2)
        if (i1.NE.0.OR.i2.NE.0.OR.i3.NE.0) then
          write(*,*)' Wrong keys: Ksg=',Ksg,' Kdg=',Kdg,' Kbrt=',Kbrt
          stop
        end if

        n1 =1
        ch1='s'
        j1 =1
        Kkin = 1
        n3=0
        ch3='s'
        n4=1
        ch4='s'
        n5=0
        ch5='s'
C    - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,65) Z,Nso,str1(Ksg+2*Kbrt),str2(Kdg)
        write(11,65) Z,Nso,str1(Ksg+2*Kbrt),str2(Kdg)
 65     format (4X,'Atom: Z =',F5.1,'  Nso =',I3,
     >         /4X,'Hamiltonian:',A8,
     >         /4X,'Diagonalization: ',A9)
       return
      end
C     ========================================================
      subroutine Rint
      INCLUDE "conf.par"
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Nsv/Nsv/Ii/Ii/Kt/Kt
     1        /Nhint/Nhint/MaxT/MaxT/Ecore/Ecore/FNAME/FNAME
     2        /Ksg/Ksg/Kbrt/Kbrt/Dint/Dint
       common /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /Kk/Kk(IPs)
     >        /rint1/rint1(10*IPh),iint1(10*IPh)
     >        /Nq/Nq(IPsp)/Nip/Nip(IPsp)
     >        /UP/Pa(IP6),Pc(IP6),Pd(IP6),Pb(IP6)
     >        /DOWN/Qa(IP6),Qc(IP6),Qd(IP6),Qb(IP6)
     >        /UP1/P1a(IP6),P1c(IP6),P1d(IP6),P1b(IP6)  !### derivatives
     >        /DOWN1/Q1a(IP6),Q1c(IP6),Q1d(IP6),Q1b(IP6)
       common /ipmr/ipmr
       dimension CPa(IP6),CQa(IP6)
       character*12 FNAME
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Nhint=0
        if (Ns.LE.Nsv) return
        fis=0.d0              !### Field IS
        sg=0.d0               !### ME of Sigma
        br=0.d0               !### Breit ME
        nx = ns - Nsv + 1     !### parameter for indexation of integrals
        irec=2*IP6*IPmr
        open(12,file=FNAME,status='OLD',access='DIRECT',recl=irec)
        open(13,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)
        call ReadF (12,1,Pa,Qa,2)
        call WriteF(13,1,Pa,Qa,2)   !### now Ecore is written to HFD.DAT
        ih=2-Kt
        nmin=Nsv+1
C     - - - - - - - - - - - - - - - - - - - - - - - - -
C       evaluation of one-electron integrals
        if (Kbrt.GE.1) write(*,*) ' Calculating Breit integrals...'
        write(11,5) 
 5      format(4X,'   na   nb',8X,'DF',13X,'QED',12X,'Br')
        do na=nmin,Ns
          call ReadF (12,na+4,Pa,Qa,2)
          call ReadF (12,na+Ns+4,CPa,CQa,2)
          if (Kbrt.EQ.2) call ReadF (13,na+Ns+4,P1a,Q1a,2)
          do nb=na,Ns
            if (Kk(na).NE.Kk(nb)) goto 20
            call ReadF (12,nb+4,Pb,Qb,2)
            if (Kbrt.EQ.2) call ReadF (13,nb+Ns+4,P1b,Q1b,2)
            do i=1,Ii,ih
              C(i)=CPa(i)*Pb(i)+CQa(i)*Qb(i)
            end do
            C(ii+4)=Pb(ii+4)+Pb(ii+4)-1.d0
            call Sint1(ds)
c            if (iabs(Kk(na)).EQ.1) then         !## Correction of the
c              call NclInt(CPa,CQa,Pb,Qb,dint1)  !## integral inside the
c              ds=ds+dint1-Dint                  !## nucleus for j=1/2
c            end if
            ierr = 0 ! cAB
            if (Ksg.EQ.2)  sg = Sigma(na,nb,ierr)  !### QED potential
            if (Kbrt.GE.1) br = Br_Core(na,nb,13)  !### Breit
            Nhint=Nhint+1
            nab=nx*(na-Nsv-1)+(nb-Nsv)
            rint1(nhint)=ds+sg+br
            iint1(nhint)=nab
            write(11,15) na,nb,ds,sg,br
 15         format(4X,2I5,4F15.8)
 20       end do
        end do
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        close(12)
        close(13)
 1000   write( *,25) Nhint,ierr
        write(11,25) Nhint,ierr
 25     format(4X,'Nhint=',I5,' Number of zeroes for Sigma:',I5)
       return
      end
C     ========================================================
      function Sigma(na,nb,ierr)
      INCLUDE "conf.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsv/Nsv/Nsint/Nsint
     >        /rint2/rint2(10*IPh),iint2(10*IPh)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c       parameter for indexation of integrals:
        nx = ns - Nsv + 1
        Sigma=0.d0

c       first time:

        if (Nsint.EQ.0) then
          ierr=0
          write(*,*) ' Reading QED matrix elements from qedpot.dat:'
          open(unit=15,file='qedpot.dat',status='OLD',err=700)
          do ir=1,10*IPh
            read(15,*,end=10,err=10) ia,ib,s
            write(*,15) ir,ia,ib,s
 15         format(1X,I4,5X,I3,5X,I3,5X,E12.5)
            nab=nx*(ia-Nsv-1)+(ib-Nsv)
            iint2(ir)=nab
            rint2(ir)=s
            Nsint=ir
          end do
 10       write( *,25) Nsint
          write(11,25) Nsint
 25       format(5X,I5,' matrix elements of QED potential read')
          close (15)
          if (Nsint.EQ.0) Stop
        end if

c       search for the radial integral:

        nab0=nx*(na-Nsv-1)+(nb-Nsv)
        do i=1,Nsint
          nab=iint2(i)
          if (nab.EQ.nab0) then
            Sigma=rint2(i)
            return
          end if
        end do
        ierr=ierr+1
       return
 700   write(*,*) ' Need file qedpot.dat for QED corrections!'
       stop
      end
C     ========================================================
      subroutine Rot(idg)
      include "conf.par"
      include "hfd.par"
      implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Nsv/Nsv/Emax/Emin,Emax/Khf/Khf/Kout/Kout
     1        /let/let(9)/Kkin/Kkin/Ierr/Ierr/Cl/Cl/FNAME/FNAME
       common /Kk/Kk(IPs)/Jj/Jj(IPs)/Ll/Ll(IPs)
     >        /P/P(IP6)/Q/Q(IP6)
     >        /Dfst/n1,j1,ch1
     >        /PW/kw(IPs),numw(IPs),nb(IPs),jjw(IPs),ltr(IPs)
       common /ipmr/ipmr
       dimension zz(50,50),hh(50,50),ez(50),dz(50)
       character*1 let,ltr,FNAME*12,ch1
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        nd1=50                    !### array dimension
        Emin= 1.d6                !### min(abs(e_i))
        Emax=-1.d6                !### max(abs(e_i))
        dhf=0.d0                  !### change of DF operator
        zmax=1.d-2                !### max |H_ik/(E_i-E_k)| gor idg=2
        trd=1.d-8                 !### convergence threshold
        irec=2*IP6*IPmr
        nmin=Nsv+1
        if (idg.NE.0) then        !### define diagonalization domain
          nmin=max(nmin,N_orb(n1,ch1,j1))
        end if

        call Count_PW(nmin,nw)    !### counting partial waves

        open(12,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)
        open(13,file=FNAME,status='OLD',access='DIRECT',recl=irec)
        do ni=1,2*Ns+4
          call ReadF (12,ni,p,q,2)
          call WriteF (13,ni,p,q,2)
        end do

        do 10 n=1,nw                    !### rotation of orbitals of
           nd=numw(n)                   !#### partial wave n

           do i=1,nd                    !### zz is transformation matrix
              do k=1,nd
                 zz(i,k)=0.d0
              end do
           end do

           i=0
           do ni=nmin,Ns                !### identification of orbitals
              if (Kk(ni).EQ.kw(n)) then !#### for p. wave n
                 i=i+1
                 nb(i)=ni
              end if
           end do
           if (i.NE.nd) then
              write(*, 15) ltr(n),jjw(n),nd,i
 15           format(4X,A1,I2,'/2 partial wave. Dimension =',I2,
     >             ' orbitals found ',I2)
              stop
           end if

           call Form_zz(nd,nd1,idg,zz)  !### forms energy matrix for p.w.

           if (idg.EQ.0) then           !### no diagonalization
             do i=1,nd
               ez(i)=zz(i,i)
               zz(i,i)=1.d0             !### zz_ik=delta_ik
             end do
           end if

           if (idg.EQ.2) then            !### this key works ONLY when
             z_nd=0.d0                   !#### starting approximation
             do i=1,nd                   !##### is good
               do k=i+1,nd
                 zik=dabs(zz(i,k)/(zz(i,i)-zz(k,k)))
                 if (z_nd.LT.zik) z_nd=zik
               end do
             end do
             if (z_nd.GT.zmax) then      !### switch to normal
               idg=1                     !#### diagonalization
               write(*,*)' z_nd=',z_nd,' idg changed to 1. Push...'
               read(*,*)
             else
               do i=1,nd
                 do k=1,nd
                   hh(i,k)=zz(i,k)
                   zz(i,k)=0.d0
                 end do
                 ez(i)=hh(i,i)
                 zz(i,i)=1.d0
               end do
               call I_Diag(nd,nd1,hh,zz,ez,dz,erx)
               if (erx.LT.trd) erx=trd
               Ierr=Ierr+dlog(erx/trd)
               write( *,*)' Iter_Diag: Ierr=',Ierr,' erx=',erx
               write(11,*)' Iter_Diag: Ierr=',Ierr,' erx=',erx
             end if
           end if

           if (idg.EQ.1) then            !### use idg=2 for more accurate
             call Hould(nd,nd1,dz,ez,zz) !#### diagonalization
           end if

           if (Kout.GT.0)
     >       write(*, 25) ltr(n),jjw(n),(ez(i),i=1,nd)
           write(11,25) ltr(n),jjw(n),(ez(i),i=1,nd)
 25        format(1X,'Eigen values for wave: (',A1,I2,'/2)',
     >            /(4E15.8))
           if (idg.NE.0) then
             if (Kout.GT.0)
     >         write(*, 35) (dabs(zz(i,i)),i=1,nd)
             write(11,35) (dabs(zz(i,i)),i=1,nd)
 35          format(1X,'Diagonal coefficients:',/(7F10.5))
             do j=2,nd                         !### test of unitarity
               do k=1,j-1
                 djk=0.d0
                 do i=1,nd
                   djk=djk+zz(i,j)*zz(i,k)
                 end do
                 if (dabs(djk).GT.1.d-8) then
                   write( *,45) j,k,djk
                   write(11,45) j,k,djk
 45                format(1X,'Diagonalization error: <',
     >                    I3,'|',I3,'> =',E10.3)
                   read(*,*)
                 end if
               end do
             end do
           end if

           dhf0=dhf
           call New_Orb(nd,nd1,zz,ez,dhf)    !# formation of new orbitals

           if (dhf-dhf0.GT.1.d-6) then
             write(*,*) ' PW ',n,': dif =',dhf-dhf0
             write(*,'(10i5)') (nb(i),i=1,nd)
                        
             do i=1,nd
               if (nb(i).LE.Nso) then
                 write(*,355) (zz(j,i),j=1,nd)
 355             format(/8f10.6,/(8f10.6))
               end if
             end do
           end if

 10     end do

        write( *,55) Emin,Emax,Kkin
        write(11,55) Emin,Emax,Kkin
 55     format(4X,'E_min = ',E10.3,'; E_max = ',E10.3,'; Kkin=',I2)

        if (nmin.LE.Nso.AND.idg.GT.0) then
          write( *,65) dhf,Ierr
          write(11,65) dhf,Ierr
 65       format(4X,'Core orbitals have changed by ',E10.3,' Ierr=',I4)
          if (dhf.GT.1.d-6) then
            Ierr=Ierr+20*dlog(dhf/1.d-6)
          else
            Khf=1
            write( *,75)
            write(11,75)
 75         format(4X,'Self-consistency is reached.')
          end if
        end if
       return
      end
C     ========================================================
      subroutine Count_PW(nmin,nw)    !### counts partial waves
      include "conf.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/let/let(9)
       common /Kk/Kk(IPs)/Jj/Jj(IPs)/Ll/Ll(IPs)
     >        /PW/kw(IPs),numw(IPs),nb(IPs),jjw(IPs),ltr(IPs)
       character*1 let,ltr
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,5) nmin
        write(11,5) nmin
 5      format(4X,'Counting partial waves starting from orbital ',I3)

        nw=1                      !### p.w. number
        kw(1)=Kk(nmin)
        jjw(1)=jj(nmin)
        ltr(1)=let(Ll(nmin)+1)
        numw(1)=1
        do ni=nmin+1,Ns
           n=0
           ki=Kk(ni)
           do i=1,nw
              if (ki.EQ.kw(i)) then
                 n=n+1
                 i1=i
              end if
           end do
           if (n.EQ.0) then
              nw=nw+1
              kw(nw)=ki
                if(ki.LT.0) then
                 jjw(nw)= -2*ki-1
                 ltr(nw)= let(iabs(ki))
                else
                 jjw(nw)= 2*ki-1
                 ltr(nw)= let(ki+1)
                end if
              numw(nw)=1
           else
              numw(i1)=numw(i1)+1
           end if
        end do
        write (*, 15) (ltr(n),jjw(n),numw(n),n=1,nw)
        write (11,15) (ltr(n),jjw(n),numw(n),n=1,nw)
 15     format(3X,A1,I2,'/2','-wave: ',I2,' orbitals')
       return
      end
C     ========================================================
      subroutine Form_zz(nd,nd1,idg,zz)
      include "conf.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nsv/Nsv/Nhint/Nhint
     >        /rint1/rint1(10*IPh),iint1(10*IPh)
     >        /PW/kw(IPs),numw(IPs),nb(IPs),jjw(IPs),ltr(IPs)
       character*1 ltr
       dimension zz(nd1,nd1)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        nx = Ns - Nsv + 1         !### used for indexation of integrals
        do i=1,nd
           ni=nb(i)
           if (idg.EQ.0) then
             kx=i
           else
             kx=Nd
           end if
           do k=i,kx
              nk=nb(k)
              nik=nx*(ni-Nsv-1)+(nk-Nsv)
              do ind=1,Nhint
                 if (nik.EQ.Iint1(ind)) then
                    t=Rint1(ind)
                    zz(i,k)=t
                    zz(k,i)=t
                    goto 200
                 end if
              end do
              write(*,*) 'no integral for i,k=',i,k
              stop
 200          continue
           end do

c         write(11,'(/I3,/(5E15.7))') i,(zz(i,k),k=1,nd)

        end do
       return
      end
C     ========================================================
      subroutine New_Orb(nd,nd1,zz,ez,dhf)   !### Orbitals after
      include "conf.par"                     !#### diagonalization
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Ii/Ii/Emax/Emin,Emax
       common /P/P(IP6)/Q/Q(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
     >        /PW/kw(IPs),numw(IPs),nb(IPs),jjw(IPs),ltr(IPs)
     >        /UP/A(IP6),B(IP6),CA(IP6),CB(IP6)
       dimension zz(nd1,nd1),ez(nd1)
       character*1 ltr
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        small=1.d-15                     !### cutoff parameter
        crit=0.9d0                       !### how to calculate derivatives
        n_num=0
        if (nb(1).LE.Nso)
     >    call Check_core(nd,nd1,zz,ez)  !### Search for ghosts in the core
        do i=1,nd
          ni=nb(i)
          do j=1,IP6
            a(j) =0.d0
            b(j) =0.d0
            ca(j) =0.d0
            cb(j) =0.d0
          end do

          do k=1,nd
            c=zz(k,i)
            if (dabs(c).LT.small) goto 200
            nk=nb(k)
            call ReadF (13,nk+4,p,q,2)
            if (ni.LE.Nso.AND.nk.GT.Nso) dhf=dhf+dabs(c)
            if (ni.EQ.nk) sorb=dabs(c)
            do j=1,IP6
              if (j.LE.ii.OR.j.GE.ii+5) then
                a(j) = a(j) + c*p(j)
                b(j) = b(j) + c*q(j)
              else
                a(j) = p(j)
                b(j) = q(j)
              end if
            end do
 200      end do

          a(ii+1)=-ez(i)                   !### new energy assigned
          if ( Emax.LT.ez(i) ) emax=ez(i)
          if ( Emin.GT.ez(i) ) emin=ez(i)

          if (sorb.LT.crit) then           !### direct differentiation
            call Dif(a,ca,a(ii+4))         !#### is used only for large
            call Dif(b,cb,a(ii+4))         !##### rotations
            n_num=n_num+1
          else

            do k=1,nd
              c=zz(k,i)
              if (dabs(c).LT.small) goto 210
              nk=nb(k)

              call ReadF (13,nk+Ns+4,cp,cq,2)
              do j=1,IP6
                if (j.LE.ii.OR.j.GE.ii+5) then
                  ca(j) = ca(j) + c*cp(j)
                  cb(j) = cb(j) + c*cq(j)
                else
                  ca(j) = cp(j)
                  cb(j) = cq(j)
                end if
              end do
 210        end do

          end if

          if (a(1).LT.0.d0) then           !### convension: a(0)>0
            do j=1,IP6
              if (j.LE.ii.OR.j.GE.ii+5) then
                a(j) = -a(j)
                b(j) = -b(j)
                ca(j) = -ca(j)
                cb(j) = -cb(j)
              end if
            end do
          end if

          call WriteF (12,ni+4,a,b,2)
          call WriteF (12,ni+Ns+4,ca,cb,2)
        end do

        if (n_num.GT.0) then
          write( *,5) n_num
          write(11,5) n_num
 5        format(4X,'New_Orb: num dif-n used for',I3,' orbitals')
        end if
       return
      end
C     ========================================================
      subroutine Check_core(nd,nd1,zz,ez)    !### Search for ghosts
      include "conf.par"                     !#### in the core
       implicit real*8 (a-h,o-z)
       common /Nso/Nso
       common /PW/kw(IPs),numw(IPs),nb(IPs),jjw(IPs),ltr(IPs)
       dimension zz(nd1,nd1),ez(nd1)
       character*1 ltr
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        itr=0
        i1=0
        do i=1,nd                            !### find core orbitals
          if (nb(i).LE.Nso) i1=i1+1
        end do
        write(11,*) i1,' core orbitals found'

        do i=1,i1
          m=0
          zx=0.d0
          do k=1,nd                          !### find eigenvector with max
            if (dabs(zz(i,k)).GT.zx) then    !#### weight of orbital nb(i)
              zx=dabs(zz(i,k))
              m=k
            end if
          end do

          if (m.NE.i) then                   !### transpose eigenvectors
            itr=itr+1
            eghost=ez(i)
            ez(i)=ez(m)
            ez(m)=eghost
            write( *,5) i,m,eghost
            write(11,5) i,m,eghost
 5          format(4X,'Transposing vectors',I3,' and',I3,
     >             ', e_ghost=',E13.5)
            do k=1,nd
              v=zz(k,i)
              zz(k,i)=zz(k,m)
              zz(k,m)=v
            end do
          end if
        end do
        if (itr.GT.0) then
          write( *,15) itr
          write(11,15) itr
 15       format(4X,'Ghosts detected! ',I3,' transpositions made.')
          read(*,*)
        end if
       return
      end
