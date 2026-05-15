C              ### file bas1.f ### last update 13/04/08
C     13/04/08: subr. Init changed to allow long basis sets
C     25/11/06: Array Qq for fractional occupancies restored
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
      subroutine Init(Norb)
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "phys.par"
      include "hfd.par"
       common /Ns/Ns/Nsp/Nsp/Nso/Nso/Nc/Nc/Ne/Ne/Nst/Nst
     1        /II/II/Kt/Kt/Kv/Kv/Ng/Ng/Mj/Mj/Z/Z/Cl/Cl
     2        /Jm/Jm/H/H/Rnuc/Rnuc/FNAME/FNAME/Kfile/Kfile
       common /AL/AL/BT/BT/R2/R2	! cAB
       common /Qnl/Qnl(IPsp)/Nq/Nq(IPsp)/Nip/Nip(IPsp)
     >        /Kk/Kk(IPs)/Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)
     >        /Qq/Qq(IPs)
     >        /R/R(IP6)/V/V(IP6)
     >        /Nvc/Nvc(IPc)/Nc0/Nc0(IPc)/Kbas/Kbas(IPsp)
       common /ipmr/ipmr
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
        z1  =PQ(1)
        if (dABS(Z-z1).GT.1.d-6) goto 730
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
        write( *,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc
        write(11,5) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc
 5      format (/4X,'Z   = ',F6.2,5X,'Kt  =',I3,  7X,'II =',I4,
     >          /4X,'H   =',F7.4, 5X,'R2  =',F6.2,4X,'R1 =',E11.4,
     >          /4X,'Rnuc=',E11.4,1X,'Al  =',F7.4,3X,'Bt =',F5.2,
     >          /4X,'Nsp =',I5,   7X,'Ns  =',I3,  7X,'Nso=',I3,
     >           4X,'Nc =',I4)
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
        Norb=0
        Kfile=0
        nsb=Ns
        do nj=1,Nsp
          i=dsign(1.d0,Qnl(nj))
          d=dabs(Qnl(nj))+1.d-14
          d=10.0*d
          nnj=d
          d=10.0d0*(d-nnj)
          llj=d
          jjj=2*llj+i
          kkj=-i*((jjj+1)/2)
          d=100.0d0*(d-llj)
          Nq(nj)=d+0.1d0
          if (nj.LE.Nso) then
            do ni = 1,Nso
              if (Nn(ni).EQ.Nn(nj).AND.Ll(ni).EQ.Ll(nj))
     >          Qw(ni) = Qw(ni) + Nq(nj)   ! number of e on NR shell
            end do
          end if
          do i=1,nsb
            ni=i
            if (nnj.EQ.Nn(ni).AND.Kk(ni).EQ.kkj) goto 210
          end do
          Norb=Norb+1
          ni=Ns+Norb
          nsb=Ns+Norb
          Nn(ni)=nnj
          Kk(ni)=kkj
          Ll(ni)=llj
          Jj(ni)=jjj
          if (Kbas(ni).EQ.3) Kfile=1
          if (Kbas(ni).EQ.4) kfile=1
 210      Nip(nj)=ni
        end do
        Write(*,*) Norb, ' new orbitals to be formed'
        nsmax=(4*IP6-20)/6
        nsmax1=(4*IP6-20)/3
        if (nsb.GT.nsmax.AND..NOT.longbasis) then
          write(*,*) ' Ns =',nsb,' > ',nsmax
          write(*,*) ' switch to long basis variant '
c          read(*,*)
        end if
        if (nsb.GT.nsmax1) then
          write(*,*) ' For IP6 =',IP6
          write(*,*) ' maximum length of basis set is ',nsmax1
          stop
        end if
        if (Nso.EQ.0) goto 220
        do ni=1,Nso
          Qq(ni) = Qw(ni)*(Jj(ni)+1.d0)/(4*Ll(ni)+2)    !# N_e on a shell
        end do
 220    n0=0
        ne=0
        nmin=Nso+1
        if (nmin.GT.Nsp) goto 230
        do ni=nmin,Nsp
          n0=n0+Nq(ni)
        end do
        ne=n0/Nc
 230    Nst=0
        do ni=1,Ns
          imax=2*Jj(ni)+1
          do j=1,imax,2
            Nst=Nst+1
          end do
        end do
        write( *,35) ne,Nst
 35     format(4X,'Ne  =',I3,7X,'Nst =',I4)
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        do ni=nmin,Nsp
          i=i+1
          n=n+Nq(ni)
          if (n.LT.ne) goto 240
          ic=ic+1
          if (n.GT.ne) goto 740
          Nvc(ic)=i
          Nc0(ic)=Nso+i0
          i0=i0+i
          n=0
          i=0
 240      continue
        end do
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
          call WriteF(13,ni+4+Ns+Norb,P,Q,2)
        end do
        close(13)
        close(12)
       Return
C     - - - - - - - - - - - - - - - - - - - - - - - - -
700     write( *,75)
        write(11,75)
75      format(/2X,'file HFD.DAT is absent'/)
       STOP
C     - - - - - - - - - - - - - - - - - - - - - - - - -
730     write( *,105) Z,Z1
        write(11,105) Z,Z1
105     format(/2X,'nuc. charge is changed'/2X,'Z1=',F12.6/2X,
     1         'Z2=',F12.6)
       stop
C     - - - - - - - - - - - - - - - - - - - - - - - - -
740     write( *,115) IC
        write(11,115) IC
115     format(/2X,'wrong number of electrons'/
     1         2X,'for configuration ICONF =',I4/)
       stop
      end
C     =================================================
      subroutine Orbit2 (P,Q,A,B,kkj,llj)
c lower component of the new orbital is expressed as follows:
c    kkin=0:     from Dirac equation with  E=0 and V=0
c    kkin=1,2:   from Dirac equation with  E=0 and V=Z/r
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /Ii/Ii/Z/Z/Cl/Cl/Kkin/Kkin/MaxT/MaxT
       common /R/R(IP6)/Y/Y(IP6)
       dimension P(IP6),Q(IP6),A(IP6),B(IP6)
       character*7 str2(3)
        data str2 /'   0   ','V_c-Z/r','V_c-Z/r'/
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        gj=P(ii+4)
!        write(*,*) ' Orbit2 calls Dif(P)'
        call Dif(P,A,gj)              !### first derivative
        if (Kkin.EQ.2) then
          call NonRelE(P,A,llj,e_nr)  !### Non-relativistic energy
        else
          e_nr=0.d0
        end if
        P(ii+1)=-e_nr

        write( *,5) str2(Kkin+1),e_nr
        write(11,5) str2(Kkin+1),e_nr
5       format(4X,'Orbit2: V =',A7,' E =',E12.5)

        do i=1,ii                     !### small component Q
          ee=2*Cl
          if (Kkin.NE.0) then
            ee=ee+(e_nr-Y(i))/Cl
          end if
          Q(i)=-1.d0/ee*(A(i)+kkj/R(i)*P(i))
        end do
        Q(ii+4)=gj
        call Origin(Q,gj,-kkj)
        write(*,*) ' Orbit2 calls Dif(Q)'
        call Dif(Q,B,gj)              !### first derivative

        call Test_Origin('Orbit_2: P',P,1.d-6,ir1)
        call Test_Origin('Orbit_2: A',A,2.d-3,ir2)
        call Test_Origin('Orbit_2: Q',Q,1.d-5,ir3)
        call Test_Origin('Orbit_2: B',B,2.d-2,ir4)
        ir=ir1+ir2+ir3+ir4
        if(ir.GT.0) then
          write(11,15) ir1,ir2,ir3,ir4
          write( *,15) ir1,ir2,ir3,ir4
 15       format(4X,'Orbit2: Bad matching at the origin.',
     >           //4X,'Errors for P, P`, Q, Q` =',4I4)
          if (ir.GT.2) read(*,*)
          Ierr=Ierr+ir
        end if
       return
      end
C     =================================================
      subroutine NonRelE(P,CP,l,e_nr)      !### Non-relativistic energy
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /Ii/Ii
       common /R/R(IP6)/C/C(IP6)/Y/Y(IP6)
       dimension P(IP6),CP(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        xl=l*(l+1)
        do i=1,ii
          C(i)=P(i)**2
        end do
        C(ii+4)=2*P(ii+4)
        call Sint1(s_nr)            !### non-relativistic normalization
        s_nr=1.d0/s_nr
        do i=1,ii
          s_kin=0.5d0*(CP(i)**2+xl*(P(i)/R(i))**2)
          s_pot=P(i)**2*Y(i)
          C(i)=s_nr*(s_kin+s_pot)
        end do
        C(ii+4)=2*P(ii+4)-2
        call Sint1(e_nr)
       return
      end
C     =================================================
      subroutine V_core      !### calculates V_core-Z/R
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "hfd.par"
       common /Ii/Ii/Z/Z
     >        /R/R(IP6)/C/C(IP6)/Y/Y(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,IP6
           C(i)=0.d0                !C is used in Y0
        end do
        call Y0(Y)                  !Y=V_core
        do i=1,ii
          Y(i)=Y(i)-Z/R(i)
        end do
        Y(ii+5)=Y(ii+5)-1.5d0*Z/R(1)
        Y(ii+7)=Y(ii+7)+0.5d0*Z/R(1)
       return
      end
C     =================================================
      function dV_nuc(P,Q,A,B)  !### calculates volume shift
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "hfd.par"
       common /Ii/Ii/Z/Z/MaxT/MaxT/R/R(IP6)
       dimension P(IP6),Q(IP6),A(IP6),B(IP6)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=R(1)
        gm=2*P(ii+4)
        V0=3*Z*r1**gm

        de=0.d0
        do m=0,MaxT
          im=ii+m+5
          s=0.d0
          do k=0,m
            ik=ii+k+5
            s=s+P(ik)*A(im-k)+Q(ik)*B(im-k)
          end do
          de=de+s/((gm+m+1)*(gm+m+3))
        end do
        dV_nuc=V0*de
       return
      end
C     =================================================
      subroutine Change_Q(ni)
c lower component of the orbital 'ni' is constructed from
c Dirac equation and compared with that from the file.
c For Kkin>0 Dirac equation is solved for V=-Z(r)/r, otherwise V=0
c For Kkin>1 E=-P(Ii+1), otherwise E=0
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "hfd.par"
       common /Ns/Ns/Ii/Ii/Z/Z/Cl/Cl/let/let(9)/Ierr/Ierr/Kkin/Kkin
     >        /R/R(IP6)/P/P(IP6)/Q/Q(IP6)/Y/Y(IP6)
     >        /C/C(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
     >        /Jj/Jj(IPs)/Nn/Nn(IPs)/Ll/Ll(IPs)/Kk/Kk(IPs)
       character*1 let
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        small=1.d-6
        kn=Kk(ni)
        n =Nn(ni)
        j =Jj(ni)
        l1=Ll(ni)+1
        call ReadF(12,ni+4,P,Q,2)

c        write(*,5) n,let(l1),j
c  5     format(' Change_Q for ',I2,A1,I1,'/2')

        gj=P(Ii+4)
        if (Kkin.EQ.2) then
          ep=-P(ii+1)
        else
          ep=0.d0
        end if
        do i=1,Ii
          C(i)=Q(i)**2
        end do
        C(Ii+4)=2*P(Ii+4)
        call Sint1(s)

c        write(*,*) ' gj=',gj,' ep=',ep,' s=',s

        call ReadF(12,ni+ns+4,CP,CP,1)
        d=2*Cl+ep/Cl
        do i=1,Ii
          ee=d
          if (Kkin.GE.1) ee=d-Y(i)/Cl
          q1=-1.d0/ee * (CP(i) + kn/R(i)*P(i))
          dq=q1-Q(i)
          C(i)=dq**2
c          if (mod(i,10).EQ.0) write(*,*) R(i),Q(i),q1
          Q(i)=q1
        end do
        C(Ii+4)=2*P(Ii+4)
        call Sint1(s1)
        s1=dsqrt(s1/s)

        if (s1.GT.small) then      !### Q1 substitutes Q
          call Origin(Q,gj,-kn)    !### Expansion at the origin
          call Dif(Q,CQ,gj)        !### derivative of the lower component
          call Test_Origin('Change_Q:Q',Q,1.d-5,ir1)
          call Test_Origin('Change_Q:B',CQ,2.d-2,ir2)
          call WriteF(12,ni+4,P,Q,2)
          call WriteF(12,ni+ns+4,CP,CQ,2)
          Ierr=Ierr+dlog(s1/small)+ir1+ir2

          write( *,25) n,let(l1),j,ep,dsqrt(s),s1
          write(11,25) n,let(l1),j,ep,dsqrt(s),s1
 25       format(4X,I2,A1,I2,'/2 (E=',F12.6,
     >          '): ||Q|| =',F10.6,' ||dQ||/||Q||=',F10.6)
        end if
c        read(*,*)
       return
      end
C     =================================================
      subroutine Ort(ni,ns1,ir)  !### ir - linear dependence parameter
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "hfd.par"
       common /Ns/Ns/II/II/Kt/Kt/Ierr/Ierr/MaxT/MaxT
       common /Kk/Kk(IPs)/C/C(IP6)/R/R(IP6)/P/P(IP6)/Q/Q(IP6)
       dimension S(IPs),A(IP6),B(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        small=1.d-6
        ir=0
c        write(*,*)' Orthogonalizing orbital ',ni
        call ReadF (12,ni+4,P,Q,2)           !### -  normalization
        ih=2-Kt
        do i=1,ii,ih
           C(i)=P(i)**2+Q(i)**2
        end do
        C(ii+4)=2*P(ii+4)
        call Sint1(ds)
        t=1.d0/dsqrt(ds)
        if (dabs(t-1.d0).GT.0.1d0) then
          ir=ir+1
          write ( *,5) ni,ni,ds
          write (11,5) ni,ni,ds
 5        format(4X,'<',I3,'|',I3,'> =',F15.7)
        end if
        call ReadF(12,ni+ns1+4,A,B,2)
        call Test_Origin('Ort0:dP/dr',A,2.d-3,ir1)
        call Test_Origin('Ort0:dQ/dr',B,2.d-2,ir2)
        ir=ir+ir1+ir2
        do i=1,IP6
          if (i.LE.ii.OR.i.GE.ii+5) then
            P(i)=P(i)*t
            Q(i)=Q(i)*t
            A(i)=A(i)*t
            B(i)=B(i)*t
          end if
        end do

        call WriteF(12,ni+4,P,Q,2)
        call WriteF(12,ni+ns1+4,A,B,2)
        call Test_Origin('Ort1:    P',P,1.d-6,ir1)
        call Test_Origin('Ort1:    Q',Q,1.d-5,ir2)
        call Test_Origin('Ort1:dP/dr',A,2.d-3,ir3)
        call Test_Origin('Ort1:dQ/dr',B,2.d-2,ir4)
        ir=ir+ir1+ir2+ir3+ir4
        if (ni.EQ.1) return

        sum=0.d0
        do 20 nj=1,ni-1                      !### - orthogonalization
          if (Kk(ni).NE.Kk(nj)) goto 20
          call ReadF (12,nj+4,A,B,2)
          do i=1,ii,ih
           C(i)=P(i)*A(i)+Q(i)*B(i)
          end do
          C(ii+4)=P(ii+4)+A(ii+4)
          call Sint1(ds)
          S(nj)=ds
          sum=sum+ds**2
          if (ds**2.GT.0.9999d0) then
            write(*,15) ni,nj,ds
 15         format(4X,'Fatal: <',I3,'|',I3,'> =',F9.6)
            stop
          end if

          d=1.d0/dsqrt(1.d0-ds*ds)
          do i=1,ii,ih
            P(i)=(P(i)-ds*A(i))*d
            Q(i)=(Q(i)-ds*B(i))*d
          end do

          pp=0.d0                            !### pp and qq are used to
          qq=0.d0                            !#### check the expansion
          do i=ii+5,ii+5+MaxT                !##### at the origin
            P(i)=(P(i)-ds*A(i))*d
            Q(i)=(Q(i)-ds*B(i))*d
            pp=pp+P(i)
            qq=qq+Q(i)
          end do
          rr=R(1)
          g=P(ii+4)
          rg=rr**g
          if (pp.NE.0.d0) then
            pp=P(1) / (rg*pp)
          else
            pp=1+P(1)
          end if
          if (qq.NE.0.d0) then
            qq=Q(1) / (rg*qq)
          else
            qq=1+Q(1)
          end if
          do i=ii+5,ii+5+MaxT  !### Rescaling of
            P(i)=pp*P(i)       !#### Taylor expansion
            Q(i)=qq*Q(i)
          end do

          erp=dabs(pp-1.d0)/small   !### - Expansion error for P
          erq=dabs(qq-1.d0)/small   !### - Expansion error for Q
          err=0.2d0/(1.d0-dabs(ds)) !### - Linear dependence
          ers=erp+erq+err

          if (ers.GT.1.d0) then     !### Error printout
            ir=ir+err
            write(11,25) ni,nj,ds,pp,qq
 25         format(4x,'<ni=',I3,'|nj=',I3,'> =',F9.6,
     >             ' Origin: P1/pp =',F9.6,' Q1/qq =',F9.6)
            if (ers.GT.5.d0) write(*,25) ni,nj,ds,pp,qq
            if(ers.GT.20.d0) then
              write (*,35) ir,ni,ers
 35           format(4X,'Severe: linear dependence ',I4,
     >               ' for orbital ',I3,' err=',F8.1)
            end if

            if(erp+erq.GT.1) then
              write ( *,45) ni
              write (11,45) ni
 45           format(4X,'Taylor expansion changed for ',I3)
              call Test_Origin('Ort2:    P',P,1.d-6,ir1)
              call Test_Origin('Ort2:    Q',Q,1.d-5,ir2)
              ir=ir+ir1+ir2
            end if
          end if
 20     continue

        if (P(1).LT.0.d0) then     !### >>> Phase convention: P(0)>0 <<<<
          iphase=-1
          do i=1,IP6
            if (i.LE.ii.OR.i.GE.ii+5) then
              P(i) = -P(i)
              Q(i) = -Q(i)
            end if
          end do
        else
          iphase=1
        end if
        call WriteF(12,ni+4,P,Q,2)

c Derivatives of P & Q:

        if (sum.GT.1.d-1) then               !### direct differentiation
          write(*,*) ' Dif(P):'
          call Dif(P,A,P(ii+4))
          write(*,*) ' Dif(Q):'
          call Dif(Q,B,P(ii+4))
        else                                 !### derivatives calculated
          call ReadF (12,ni+4+ns1,A,B,2)     !#### using matrix S
          snorm=1.d0                         !### snorm accounts for
          do nj=1,ni-1                       !#### the change in
            if (Kk(ni).EQ.Kk(nj)) then       !##### normalization during
              call ReadF (12,nj+4+ns1,P,Q,2) !###### consequent orthog-n
              ds=snorm*S(nj)
              d=1.d0/dsqrt(1.d0-ds*ds)
              snorm=snorm*d
              do i=1,IP6
                if (i.LE.ii.OR.i.GE.ii+5) then
                  A(i)=(A(i)-ds*P(i))*d
                  B(i)=(B(i)-ds*Q(i))*d
                end if
              end do
            end if
          end do
          do j=1,IP6
            if (j.LE.ii.OR.j.GE.ii+5) then
              A(j) =iphase*A(j)
              B(j) =iphase*B(j)
            end if
          end do
        end if
        call Test_Origin('Ort2:dP/dr',A,2.d-3,ir1)
        call Test_Origin('Ort2:dQ/dr',B,2.d-2,ir2)

        call WriteF(12,ni+4+ns1,A,B,2)
        ir=ir+ir1+ir2
        if (ir.GT.0) write(11,55) ir,ni
 55     format(4X,'Ort: linear dependence ',I4,' for orbital ',I3)
        Ierr=Ierr+ir
       return
      end
C     =================================================
      subroutine Core
      implicit real*8 (a-h,o-z)
      include "conf.par"
      include "hfd.par"
       common /Ns/Ns/Nso/Nso/ii/ii/Kt/Kt/H/H0/Kout/Kout
     1        /Z/Z/Cl/Cl/Ecore/Ecore/FNAME/FNAME/MaxT/MaxT
     2        /Nn/Nn(IPs)/Ll/Ll(IPs)/Kk/Kk(IPs)/Jj/Jj(IPs)
     3        /Qq/Qq(IPs)
       common /C/C(IP6)/R/R(IP6)/V/V(IP6)
     >        /P/P(IP6)/Q/Q(IP6)
       common /ipmr/ipmr
       dimension A(IP6),B(IP6),CP(IP6),CQ(IP6),CC(IP6),
     1           Y(IP6),dd(20),ss(20)
       character*12 FNAME
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        open(12,file=FNAME,status='OLD',
     >       access='DIRECT',recl=2*IP6*IPmr)
        Ecore=0.d0
        Hcore=0.d0
        ih=2-kt
        hh=h0*ih/3.d0
        do i=1,IP6
           C(i)=0.d0
           Y(i)=0.d0
        end do
        if (Nso.EQ.0) goto 200
        call Y0(Y)               !### Coulomb potential of the core
        err1=2.d-4
        err2=2.d-3
        nigd=0
        irr=0
 200    do ni=1,Ns               !### one-electron part of the Hamiltonian
          na=Nn(ni)
          la=Ll(ni)
          ja=Jj(ni)
          qa=Qq(ni)
          t=Cl*Kk(ni)
          call ReadF (12,ni+4,P,Q,2)
          call ReadF (12,ni+Ns+4,A,B,2)
          if (Kout.GT.0) then
             write(11,*) ' Core: forming CP and CQ for orbital ',ni
             write(*,*) ' Core: Testing P,Q,A,B for orb. ',ni,'..'
          end if
          call Test_Origin('Core:    P',P,1.d-6,ir1)
          call Test_Origin('Core:    Q',Q,1.d-6,ir2)
          call Test_Origin('Core:    A',A,err1,ir3)
          call Test_Origin('Core:    B',B,err1,ir4)
          ir=ir1+ir2+ir3+ir4
          do i=1,ii,ih
            d= Cl*B(i)-(Z*P(i)+t*Q(i))/R(i)
            s=-Cl*A(i)-(Z*Q(i)+t*P(i))/R(i)-2*Cl*Cl*Q(i)
            C(i)=d*P(i)+s*Q(i)
            CP(i)=d+Y(i)*P(i)
            CQ(i)=s+Y(i)*Q(i)
          end do
          r1=R(1)
          gm=P(ii+4)
          CP(ii+4)=gm-1.d0
          CQ(ii+4)=gm-1.d0
          igm=gm+1.d-5
          dgm=gm-igm                   !### dgm=0 for finite nucleus
          do m=0,MaxT
            m1=m+1
            im=ii+5+m
            dd(m1)=Cl*B(im)-t*Q(im)
            ss(m1)=-Cl*A(im)-t*P(im)
            if (m.GE.1) ss(m1)=ss(m1)-2*Cl*Cl*r1*Q(im-1)
            if (dgm.GT.1.d-6) then     !### case of the pointlike nucleus
              dd(m1)=dd(m1)-Z*P(im)
              ss(m1)=ss(m1)-Z*Q(im)
            else                       !### case of the finite nucleus
              if (m.GE.1) then
                dd(m1)=dd(m1)-1.5d0*Z*P(im-1)
                ss(m1)=ss(m1)-1.5d0*Z*Q(im-1)
              end if
              if (m.GE.3) then
                dd(m1)=dd(m1)+0.5d0*Z*P(im-3)
                ss(m1)=ss(m1)+0.5d0*Z*Q(im-3)
              end if
            end if
            cpp=0.d0
            cqq=0.d0
            if (m.GE.1) then
              do k=0,m-1
                ik=ii+5+k
                cpp=cpp+P(ik)*Y(im-k-1)
                cqq=cqq+Q(ik)*Y(im-k-1)
              end do
            end if
            CP(im)=dd(m1)+r1*cpp
            CQ(im)=ss(m1)+r1*cqq
          end do
c>>>> test of expansion at the origin:
c          if (Kout.GT.0)
c     >      write(*,*) ' Core: Testing CP & CQ for orb. ',ni,'...'
c          call Test_Origin('Core:   CP',CP,err2,ir1)
c          call Test_Origin('Core:   CQ',CQ,err2,ir2)
c<<<< test of expansion at the origin
          if (Nso.EQ.0) goto 240
          if (ni.LE.Nso) then
            C(ii+4)=2*P(ii+4)-1.d0
            call Sint1(ds)
            Hcore=Hcore+qa*ds
            Ecore=Ecore+qa*ds
            do i=1,ii,ih
              C(i)=Y(i)*(P(i)**2+Q(i)**2)
            end do
            C(ii+4)=2*P(ii+4)
            call Sint1(ds)
            Ecore=Ecore+0.5d0*qa*ds
          end if
          do nj=1,Nso
            call ReadF (12,nj+4,A,B,2)
            nb=Nn(nj)
            lb=Ll(nj)
            jb=Jj(nj)
            qb=Qq(nj)
!!            qb0=qb                          !# used for k=0
!!            if (na.EQ.nb.AND.la.EQ.lb) then !# ni and nj belong to
!!              qb0=(jb+1)                    !# one NR shell and
!!            end if                          !# qb0 is adjusted.
            kmin=iabs(ja-jb)/2+1
            kmax=(ja+jb)/2+1
            do k1=kmin,kmax
              k=k1-1
              ip=la+lb+k
              if (ip.EQ.2*(ip/2)) then
                do i=1,ii,ih
                  C(i)=P(i)*A(i)+Q(i)*B(i)
                end do
                C(ii+4)=P(ii+4)+A(ii+4)
                call Yk(k)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
                xja=ja/2.d0
                xjb=jb/2.d0
                xk=k
                d=qb*(FJ3(xk,xja,xjb,0.d0,-0.5d0,0.5d0))**2
                do i=1,ii,ih
                  yy=d*C(i)/(V(i)*hh*R(i))
                  CP(i)=CP(i)-yy*A(i)
                  CQ(i)=CQ(i)-yy*B(i)
c                  CC(i)=yy
                end do
c>>>> expansion at the origin (finite nucleus assumed)
c                CC(ii+4)=k
c                call Origin(CC,CC(ii+4),-1)
c                mg=A(ii+4)+k-gm+0.1d0
c                if (mg.LT.0) then
c                  write(*,*) ' Core: mg<0 for ni,nj=',ni,nj,' mg=',mg
c                  read(8,*)
c                  stop
c                end if
c                rmg=r1**(mg+1) !## Note that CP(ii+4)=gm-1, not gm!
c                do m=mg,MaxT-1
c                  im=ii+5+m
c                  cpp=0.d0
c                  cqq=0.d0
c                  do mk=0,m-mg
c                    ik=ii+5+mk
c                    cpp=cpp+CC(ik)*A(im-mk-mg)
c                    cqq=cqq+CC(ik)*B(im-mk-mg)
c                  end do
c                  CP(im+1)=CP(im+1)-rmg*cpp
c                  CQ(im+1)=CQ(im+1)-rmg*cqq
c                end do
c<<<< expansion at the origin

                if (nj.LE.ni.AND.ni.LE.Nso) then
 230              t=qa
                  if (ni.EQ.nj.AND.k.NE.0)
     >               t=0.5d0*(qa-1.d0)*(ja+1.d0)/ja
                  if (ni.EQ.nj.AND.k.EQ.0) t=0.5d0*qa
                  d=t*d
                  do i=1,ii,ih
                    C(i)=C(i)*(P(i)*A(i)+Q(i)*B(i))/R(i)
                  end do
                  C(ii+4)=P(ii+4)+A(ii+4)+k
                  call Sint(ds)             !# integration over ro!
                  Ecore=Ecore-ds*d
                end if
              end if
            end do
          end do
c>>>> test of expansion at the origin:
c          if (Kout.GT.0)
c     >      write(*,*) ' Core: Testing CP & CQ for orb. ',ni,'....'
c          call Test_Origin('Core:   CP',CP,err2,ir1)
c          call Test_Origin('Core:   CQ',CQ,err2,ir2)
c          ir=ir+ir1+ir2
          if (ir.GT.0) then
            irr=irr+1
            if (Kout.GT.0)
     >        write( *,5) ni,ir
            write(11,5) ni,ir
 5          format(' Orbital',I4,': Expansion at the origin badness ',
     >             I6)
          end if
          if(irr.EQ.0) nigd=ni
c<<<< test of expansion at the origin
 240      call WriteF(12,ni+4+Ns,CP,CQ,2)
        end do
        call ReadF (12,1,P,Q,2)
        P(18)=Ecore
        P(19)=Hcore
        call WriteF(12,1,P,Q,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,15) Hcore,Ecore
        write(11,15) Hcore,Ecore
 15     format(/4X,'one-el. and total core energy:',F17.7,
     1       4X,F17.7)
        if(irr.GT.0) then
          write( *,25) nigd,irr
          write(11,25) nigd,irr
 25       format(4X,'Only',I4,' first orbitals are good at the origin',
     >          /4X,'for other orbitals there were',I6,' errors.',
     >          /4X,'See RES file and use Kout=2 for more details.')
          if (Kout.GT.0) then
            write(*,*) '  Push...'
            read(*,*)
          end if
        end if
        close(12)
       return
      end
C     =================================================
      function N_orb(n,ch,j)
      include "conf.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Let/let(9)
       common /Nn/Nn(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
       character*1 Let,ch
        N_orb=0
        do i=1,9
          l=i-1
          if(ch.EQ.Let(i)) goto 200
        end do
        write( *,5) n,ch,j,ch,Let
        write(11,5) n,ch,j,ch,Let
 5      format(/4X,'N_orb error for input: n=',i3,', ch= ',
     >         a1,', j=',i3/4X,'no l associated with letter ',A1,
     >         /4X,'Known letters are:',9A2)
        stop
 200    do i=1,Ns
          if (n.EQ.Nn(i).AND.l.EQ.Ll(i).AND.j.EQ.Jj(i)) then
            N_orb=i
            return
          end if
        end do
      end
C     =================================================
      subroutine Tail(ni,     !### forces Large component to zero at R(ii)
     >                ns,     !### No of orbitals
     >                itail)  !### =0 no changes & =1 when changes made
      INCLUDE "conf.par"
      include "hfd.par"
       implicit real*8 (a-h,o-z)
       common /II/II/Ierr/Ierr
     >        /let/let(9)/JJ/JJ(IPs)/NN/NN(IPs)/LL/LL(IPs)
     >        /R/R(IP6)/P/P(IP6)/Q/Q(IP6)/CP/CP(IP6)/CQ/CQ(IP6)
       character*1 let
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        itail=0
        ix= 7                                   !### length of the tail
        small=1.d-6

        call ReadF (12,ni+4,P,Q,2)
        if (dabs(P(ii)).LT.small) return

        call ReadF (12,ni+ns+4,CP,CQ,2)         !### At last ix nodes
        itail=1                                 !### P(i) -> P(i)*cos(y_i)
                                                !### y_{ii-ix}=0, y_ii=pi/2
        r0=R(ii-ix)
        rx=R(ii)-r0
        c=3.141592653589793d0/(2.d0*rx)
        do i1=1,ix
          i=ii-ix+i1
          x=c*(R(i)-r0)
          pi=P(i)
          cpi=CP(i)
          dx=dcos(x)
          P(i)=pi*dx
          CP(i)=cpi*dx-c*pi*dsin(x)
        end do

        write( *,5) nn(ni),let(ll(ni)+1),jj(ni),pi,P(ii)
        write(11,5) nn(ni),let(ll(ni)+1),jj(ni),pi,P(ii)
 5      format(1X,'Tail(',I3,A1,I2,'/2): P(ii) old',E12.5,' new',E12.5)
        ierr=ierr+dlog(1.d0+dabs(pi)/small)

        call WriteF (12,ni+4,P,Q,2)
        call WriteF (12,ni+ns+4,CP,CQ,2)

       return
      end
C     =================================================
      subroutine NclInt(p,q,a,b,dn)
c radial integration inside the nucleus: int_0,1 f(x) dx
c where x = r/R(1)  and  f(x) = (p(x)*q(x) + a(x)*b(x))
      implicit real*8 (a-h,o-z)
      INCLUDE "hfd.par"
       common /Ii/Ii/MaxT/MaxT/R/R(IP6)
       dimension p(IP6),q(IP6),a(IP6),b(IP6)
        dn=0.d0
        g12=p(ii+4)+a(ii+4)
        do i=0,MaxT
           i1=Ii+5+i
           do j=0,MaxT-i
              j1=Ii+5+j
              dn = dn + (p(i1)*a(j1)+q(i1)*b(j1))
     >            /(g12+i+j+1)
           end do
        end do
        dn=dn*R(1)**(g12+1)
       return
      end
