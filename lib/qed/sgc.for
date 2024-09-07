c     ==================   10/08/95   =====================
      Program Sigma_C             !### last update 14/04/08
c     Evaluation of the SET of second order diagrams for Sigma.
c     Diagrams are reduced to the radial integrals analytically.
c     Results are WRITTEN TO SG.RES in the radial integral format.
c     Input file: MBPT.INP
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      INCLUDE "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Kt/Kt/au/au/Nsh/Nsh/Nsx/Nsx/Nmax/Nmax
     1        /Lix/Lix/Dlix/Dlx1,Dlx2/Kout/Kout/Kval/Kval/Nsv/Nsv
     2        /QNn/na,nb,nc,nd
     3        /QNl/la,lb,lc,ld
     4        /QNj/ja,jb,jc,jd
     5        /Orbit/ka,kb,kc,kd
     6        /C_SMS/C_SMS/Klow/Klow/K_inf/K_inf
       Common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Lj/Lj(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
     >        /Khot/Khot/Chot/Chot(10)
       common /ipmr/ipmr
       character*1 let(9), lbl *4
       character*128 string
       integer Mrec
        data let/'s','p','d','f','g','h','i','k','l'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c |||         Maximal values for parameters            |||
c |||         which determine array dimensions:        |||
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c                   Ns   <=   IPs
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c Mrec determine the length of the records in the DIRECT files:
        call recunit
        Mrec=IPmr
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=10,file='MBPT.INP',status='OLD')
        open(unit=11,file='SGC.CON',status='UNKNOWN')
        open(unit=12,file='HFD.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*Mrec,err=700)
        open(unit=13,file='CONF.DAT',access='DIRECT',
     >       status='OLD',recl=2*IP6*Mrec,err=700)
        write( 6,5)
 5      format(4X,'PROGRAM Sigma_C')
        Kout=0                                   !### use sg.for to see details
C     - - - - - - - - - - - - - - - - - - - - - - - - -
 200    read (10,15) lbl                         !### search for MBPT
 15     format(A4)                               !#### part of the input
        if (lbl.NE.'MBPT') goto 200
        call reads(10,string)
        read (string,*) Nso                      !### CI core
        call reads(10,string)
        read (string,*) Nsh                      !### defines SCF field
        call reads(10,string)
        read (string,*) nss                      !### last virtual shell
        call reads(10,string)
        read (string,*) nsv                      !### first vacant shell
        call reads(10,string)
        read (string,*) Nmax
        call reads(10,string)
        read (string,*) lmax
        read (10,*)
        call reads(10,string)
        read (string,*) kt1
        read (10,*)                          !### Breit is not included
c 25     format(5X,I3)
        rewind(10)
 210    read (10,15) lbl                         !### search for Sigma
        if (lbl.NE.'Sig:') goto 210              !#### part of the input
        call reads(10,string)
        read (string,*) Khot                     !### HOT key
        call reads(10,string)
        read (string,*) Nav                      !### switches nsx1 to nsx2
        call reads(10,string)
        read (string,*) nsx1                     !### active core shells
        call reads(10,string)
        read (string,*) nsx2                     !###   "      "     "
        call reads(10,string)
        read (string,*) K_inf                    !### extrapolation to L_max=infty
        if (Khot.EQ.1) then
           rewind(10)
 220       read (10,15) lbl                      !### search for HOT
           if (lbl.NE.'HOT:') goto 220           !#### part of the input
           do k=1,10
             call reads(10,string) 
             read (string,*) Chot(k)             !### Screening coefficients
c 35          format(5X,F5.3)                    !#### for multipolarity k-1
           end do
        end if
        if (Nsv.LE.Nsh) then
           rewind(10)
 230       read (10,15) lbl                      !### subtraction
           if (lbl.NE.'SBT ') goto 230           !#### part of the input
           do k=Nsv,Nsh
              read (10,45) k1,Qq(k1)             !### ocupation number
 45           format(5X,I3,10X,F8.4)             !#### of the shell k
              if (k1.NE.k) then
                write(*,*)' Error in SBT part of input'
                write(*,*)' expected shell ',k,' received ',k1
                stop
              end if
           end do
        end if
        rewind(10)
        C_SMS=0.d0
 240    read (10,15,end=250,err=250) lbl         !### search for SMS
        if (lbl.NE.'SMS:') goto 240              !#### part of the input
        read (10,'(6X,F9.6)') C_SMS              !### scaling parameter
        call reads(10,string)
        read (string,*) Klow                        !### low component (0=NO)
c        if (Nso.NE.Nsh.AND.C_SMS.NE.0.d0) then
c          write(*,*) ' SMS correction defined only for Nso=Nsh'
c          read(*,*)
c          stop
c        end if
C     - - - - - - - - - - - - - - - - - - - - - - - - -
 250    Nsx=nsx1
        call Init(Nmax)
        call CalcWig0
        Kt=min(kt1,Kt)
        Ns=min(Ns,nss)
        call DefEv(13,lmax)
        close(unit=10)
        write(6, 55) Nmax,lmax,Kt,Kval,Khot,C_SMS,Klow
        write(11,55) Nmax,lmax,Kt,Kval,Khot,C_SMS,Klow
 55     format(1X,' Nmax=',I3,'lmax=',I1,'  kt=',I2,'  kval=',I2,
     >         '  Khot=',I2,' C_SMS=',F9.6,' Klow=',I1)
        do n=Nso+1,Nmax
          ln=Ll(n)
          if (ln.LE.lmax) then
            write(6, 65) n,Nn(n),let(ln+1),Lj(n),Eps1(n),Eps2(n),Eps(n)
 65         format(1X,I3,3X,I2,A1,I1,'/2',' E =',F10.5,'(val P) ',
     >             ' E =',F10.5,'(val Q) ',F10.5,'(virt)')
          end if
        end do
        write(6, 75)
 75     format(1X,30('-'))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        irpt=0
        do ka1=Nso+1,Nmax
           na =Nn(ka1)
           la =Ll(ka1)
           lja=Lj(ka1)
           ja =lja*0.5d0
           if (la.LE.lmax) then
              do kb1=ka1,Nmax
                 ka =ka1         !### this assignment is necessary because
                 kb =kb1         !### of the ka<-->kb transposition below
                 if (ka+kb.LE.Nav*2) then
                    Nsx=nsx1
                 else
                    Nsx=nsx2
                 end if
                 nb =Nn(kb)
                 lb =Ll(kb)
                 ljb=Lj(kb)
                 jb =ljb*0.5d0
                 if (lb.EQ.la.AND.ljb.EQ.lja) then
                    irpt=irpt+1
                    dsg =0.d0                         !### d/dE Sigma
                    Dlx1=0.d0                         !### Sigma(Lix-1)
                    Dlx2=0.d0                         !### Sigma(Lix-2)
                    x = 0.d0                          !### x   = Sigma
                    if (Nsx.NE.0) then
                      call ReadF (13,ka+4,Pa,Qa,2)
                      call ReadF (13,kb+4,Pb,Qb,2)
                      x = 0.d0
     >                  + D11(dsg)
     >                  + D12(dsg)
     >                  + D13(dsg)
     >                  + D14(dsg)
c SUBTRACTION diagrams:
                      if (Nsh.GE.Nsv) then
                        x = x
     >                  + D1ab(dsg)
     >                  + D1cd(dsg)
     >                  + D1e(dsg)
     >                  + D1f(dsg)
     >                  + D1g(dsg)
     >                  + D1h(dsg)
     >                  + D1a1(dsg)
     >                  + D1b1(dsg)
     >                  + D1c1(dsg)
     >                  + D1d1(dsg)
                      end if
                    end if
                    x_inf=D_inf(x,Dlx1,Dlx2)
                    write(6, 85) na,let(la+1),lja,nb,let(lb+1),ljb,
     >                           Eps1(ka),x_inf,x_inf*au,dsg,
     >                           Lix-2,Dlx2,Lix-1,Dlx1,Lix,x,x_inf
 85                 format(1X,'>>>> Sigma(',I2,A1,I2,'/2,',I2,A1,I2,
     >                '/2, E = ',F10.6,') =',/1X,'=',E12.5,' au =',
     >                F9.1,' (1/cm) d/dE Sig =',E11.4,/3(' Sig(',
     >                I1,')=',E11.4),' Sig_inf=',E11.4,/1X,73('='))
                    write(11,95) irpt,ka1,kb1,x_inf,dsg,Eps1(ka1)
 95                 format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))
                    close(11)
                    open(unit=11,file='SGC.CON',status='OLD')
                    do ir=1,irpt+1
                       read(11,95)
                    end do
                 end if
              end do
           end if
        end do
        write(11,*)' screening coefficients:'
        do k=1,10
           if (Khot.EQ.1) then
              write(11,105) (k-1),Chot(k)
           else
              write(11,105) (k-1),1.d0
           end if
 105       format(5X,I2,5X,F8.5)
        end do
        write(11,*)' >>>>>>>>> END <<<<<<<<<'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=11)
        close(unit=12)
        close(unit=13)
        if (Nsv.NE.Nso+1) then
           write(*,*) 'Note that Nsv=',Nsv,' <> Nso+1 =',Nso+1
        end if
       stop
C     - - - - - - - - - - - - - - - - - - - - - - - - -
 700    write( 6,115)
115     format(/2X,'file HFD.DAT/CONF.DAT is absent'/)
       stop
      end
C     =================================================
      include "readf.inc"
      include "readff.inc"
      include "eval.inc"
      include "wig.inc"
      include "sg1.inc"
      include "cfnc.inc"
      include "rp2.inc"
      include "pi_pk.inc"
      include "dif.inc"
      include "sint1.inc"
      include "rec_unit.inc"
C     =================================================
      Subroutine reads(kan,string)
      character*128 string
      integer kan
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      read(kan,'(a)') string
      string(1:5)='     '
      return
      end
C     =================================================
      Subroutine Init(nmax)
      implicit real*8 (a-h,j,m,o-z)
      INCLUDE "phys.par"
      INCLUDE "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Ii/Ii/Kt/Kt/Z/Z
     1        /Lix/Lix/Cl/Cl/H/H/Rnuc/Rnuc/au/au/Jlst/Jlst/Kval/Kval
       common /Kk/Kk(IPs)/Lj/Lj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /R/R(IP6)/V/V(IP6)
     >        /Eps/Eps(IPs)
       dimension p(IP6),q(IP6),p1(IP6),q1(IP6),pq(4*IP6)
       logical longbasis
       dimension IQN(4*IPs),Qq1(IPs)
       equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
       equivalence (p(1),pq(1)), (q(1),pq(IP6+1)),
     >      (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c small number:
        c1=0.01d0
c     ===============================================
c     ===  values for speed of light and Rydberg  ===
c     ===  constant are taken from "phys.par"     ===
c     ===============================================
        Cl=DPcl
        au=2*DPRy
        call ReadF (12,1,p,q,2)
        call ReadF (12,2,R,V,2)
        call ReadF (12,3,p1,q1,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        Z   =pq(1)
        Ns  =pq(2)+c1
        Ii  =pq(3)+c1
        H   =pq(6)
        Bt  =pq(7)
        Al  =pq(8)
        Kt  =pq(9)+c1
        Rnuc=pq(13)
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
        write( 6,25) Z,Ii,Ns,Nsx,H,Rnuc,Al,Bt,Kt,Kval
c       write(11,25) Z,Ii,Ns,Nsx,H,Rnuc,Al,Bt,Kt,Kval
 25     format (4X,'Z =',F6.2,4X,'II =',I4,4X,'NS =',I3,4X,'NSX =',I3,
     1          /4X,'H =',F7.4,2X,'RNUC=',E11.4,2X,'AL  =',
     2          F7.4,2X,'BT =',F5.2,/4X,'KT=',I2,4X,'Kval=',I2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        llst=0
        Lix=0
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          do ni=1,Ns
            Nn(ni)=IQN(4*ni-3)
            Ll(ni)=IQN(4*ni-2)
            Lix=max(Lix,Ll(ni))
            Kk(ni)=IQN(4*ni-1)
            Lj(ni)=IQN(4*ni)
            if (llst.LT.Lj(ni).AND.ni.LE.Nso) llst=Lj(ni)
          end do
        else
          if=20
          do ni=1,Ns
             if=if+1
             Nn(ni)=pq(if)+c1
             if=if+1
             Ll(ni)=pq(if)+c1
             Lix=max(Lix,Ll(ni))
             if=if+3
             c2=dsign(c1,pq(if))
             Kk(ni)=pq(if)+c2
             if=if+1
             c2=dsign(c1,pq(if))
             Lj(ni)=pq(if)+c2
             if (llst.LT.Lj(ni).AND.ni.LE.Nso) llst=Lj(ni)
          end do
        end if
c>>>> check of the basis set CONF.DAT: <<<<
        ier=0
        call ReadF (13,1,p,q,2)
        call ReadF (13,3,p1,q1,2)
        Ns1 =pq(2)+c1
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
        write( 6,*) ' File CONF.DAT: Ns =',Ns1
        nmax=min(Ns1,nmax)
        if (longbasis) then
          write( *,*) ' Using variant for long basis '
          do ni=1,Ns1
            nn1=IQN(4*ni-3)
            ll1=IQN(4*ni-2)
            kk1=IQN(4*ni-1)
            if (nn1+ll1+kk1.NE.Nn(ni)+Ll(ni)+Kk(ni)) then
              write(*,*) ' CONF.DAT: unexpected orbital ',ni
c              write(*,*) ' got: n=',nn1,' l=',ll1,' k=',kk1
c              write(*,*) ' exp: n=',nn(ni),' l=',ll(ni),' k=',kk(ni)
c              read(*,*)
              ier=ier+1
            end if
          end do
        else
          if=20
          do ni=1,Ns1
            if=if+1
            nn1=pq(if)+c1
            if=if+1
            ll1=pq(if)+c1
            if=if+3
            c2=dsign(c1,pq(if))
            kk1=pq(if)+c2
            if=if+1
            if (nn1+ll1+kk1.NE.Nn(ni)+Ll(ni)+Kk(ni)) then
              write(*,*) ' CONF.DAT: unexpected orbital ',ni
c              write(*,*) ' got: n=',nn1,' l=',ll1,' k=',kk1
c              write(*,*) ' exp: n=',nn(ni),' l=',ll(ni),' k=',kk(ni)
c              read(*,*)
              ier=ier+1
            end if
          end do
        end if
        if (ier.NE.0) then
           write(*,*) ier,' errors detected in CONF.INP'
           stop
        end if
c>>>>>>>>>>>>>> end of check <<<<<<<<<<<<<<
        Jlst=llst/2.d0
        do ni=1,Ns
           call ReadF (12,ni+4,p,q,1)
           Eps(ni)=-p(Ii+1)
        end do
       Return
      end
