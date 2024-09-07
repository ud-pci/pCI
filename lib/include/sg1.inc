c### 28/12/96: d_Sigma/d_E is changed.
c### 14/12/98: in subtraction diagrams (2j_i+1) => Qq(i).
c### 08/02/99: Eps2: valence energies in Q subspace.
c### 31/01/01: dSigma/dE is printed for each diagram.
c### 02/02/01: input/output changed (/kan/ deleted).
c### 13/09/02: diagrams D1b, D1d & D1g are added.
c### 18/01/04: function Wig0 moved to the file wig0.inc
c### 27/01/04: functions Norb,Kph,CFnc => to cfnc.inc
c### 17/05/05: Lix = max_l for intermediate states and xs1,xs2
c###           to extrapolate Sigma to Lix->infty.
c### 08/11/06: subtraction diagrams changed to account for SMS
c### 14/04/08: Parameter IPs is introduced
C     =================================================
      include "yk.inc"
C     =================================================
      Function D11(dsg)
c     Evaluation of the first diagram for Sigma:
c
c       a _________ e _______ b
c             (  __ d __(
c              )/       \)
c               \__ c __/
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc1=1,Nsx
          kc=Nso+1-kc1
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon and 'e' state:
            kpp=Isig(lc-ld)
            kpe=Isig(la)*kpp
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(jc-jd,kpp)
            kmax=jc+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s1=(2*jc+1)*(2*jd+1)/(2*k1+1)*Wig0(k1,jc,jd)**2
                call CFnc(k1,ljc,lc,ljd,ld,Pc,Qc,Pd,Qd,c1)
                do 20 ke=Nsv,Ns
                  le=Ll(ke)
                  lje=Lj(ke)
                  je=lje/2.d0
c     restrictions on the last orbital:
                  if (Isig(le).NE.kpe)          goto 20
                  if (dabs(ja-je).GT.(k1+0.1d0))goto 20
                  if ((ja+je).LT.(k1-0.1d0))    goto 20
                  lmx=max(ld,le)                  !### max l of virt.states
                  s=s1*(2*je+1)*Wig0(k1,ja,je)**2
                  e=Eps1(ka)+Eps(kc)-Eps(kd)-Eps(ke)
                  if (Khot.EQ.1) then             !### HOT correction
                    s=s*Chot(k1+1)
                  end if
                  if (s.NE.0.d0) then
                    dt=s*Rp2(ke,lja,la,ljb,lb,c1,c1,Pa,Qa,Pb,Qb)/e
                    xs=xs+dt
                    if (lmx.LE.Lix-1) xs1=xs1+dt
                    if (lmx.LE.Lix-2) xs2=xs2+dt
                    ds=ds-dt/e
                  end if
 20             end do
              end do
            end if
          end do
        end do
        D11=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D11 (L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D11 (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D12(dsg)
c     Evaluation of the second diagram for Sigma:
c
c       a _________ e __
c              )___ c _/)k1
c           k2(/___ d _(_____ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc1=1,Nsx
          kc=Nso+1-kc1
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photons and 'e' state:
            kp1=Isig(lb-ld)
            kp2=Isig(ld-lc)
            kpe=Isig(la)*kp2
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for the first photon momentum:
            kmi1=Kph(jb-jd,kp1)
            kma1=jb+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for the second photon momentum:
            kmi2=Kph(jd-jc,kp2)
            kma2=jd+jc+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmi1.LE.kma1.AND.kmi2.LE.kma2) then
              do k1=kmi1,kma1,2
                call CFnc(k1,ljb,lb,ljd,ld,Pb,Qb,Pd,Qd,c1)
                do k2=kmi2,kma2,2
                  s1=Isig(k1+k2)*(2*jc+1)*(2*jd+1)
     >                 *Wig0(k1,ja,jd)*Wig0(k2,jc,jd)
                  call CFnc(k2,ljc,lc,ljd,ld,Pc,Qc,Pd,Qd,c2)
                  do 20 ke=Nsv,Ns
                    le=Ll(ke)
                    lje=Lj(ke)
                    je=lje/2.d0
c     restrictions on the last orbital:
                    if (Isig(le).NE.kpe)          goto 20
                    if (dabs(ja-je).GT.(k2+0.1d0))goto 20
                    if ((ja+je).LT.(k2-0.1d0))    goto 20
                    if (dabs(jc-je).GT.(k1+0.1d0))goto 20
                    if ((jc+je).LT.(k1-0.1d0))    goto 20
                    lmx=max(ld,le)                !### max l of virt.states
                    s=s1*FJ6(ja,jd,k1+0.d0,jc,je,k2+0.d0)
                    if (s.NE.0.d0) then
                      s=s*(2*je+1)*Wig0(k1,jc,je)*Wig0(k2,ja,je)
                      e=Eps1(ka)+Eps(kc)-Eps(kd)-Eps(ke)
                      if (Khot.EQ.1) then         !### HOT correction
                        s=s*Chot(k1+1)*Chot(k2+1)
                      end if
                      dt=s*Rp2(ke,ljc,lc,lja,la,c1,c2,Pc,Qc,Pa,Qa)/e
                      xs=xs+dt
                      if (lmx.LE.Lix-1) xs1=xs1+dt
                      if (lmx.LE.Lix-2) xs2=xs2+dt
                      ds=ds-dt/e
                    end if
 20               end do
                end do
              end do
            end if
          end do
        end do
        D12=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D12 (L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D12 (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D13(dsg)
c     Evaluation of the third diagram for Sigma:
c
c       a _______________
c               ___ e __/
c              /________ )___ b
c             (  __ d __(
c              )/       \)
c               \__ c __/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc1=1,Nsx
          kc=Nso+1-kc1
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon and 'e' state:
            kpp=Isig(lc-ld)
            kpe=Isig(la)*kpp
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(jc-jd,kpp)
            kmax=jc+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s1=-(2*jc+1)*(2*jd+1)/(2*k1+1)*Wig0(k1,jc,jd)**2
                call CFnc(k1,ljd,ld,ljc,lc,Pd,Qd,Pc,Qc,c1)
                do 20 ke1=1,Nsx
                  ke=Nso+1-ke1
                  le=Ll(ke)
                  lje=Lj(ke)
                  je=lje/2.d0
c     restrictions on the last orbital:
                  if (Isig(le).NE.kpe)          goto 20
                  if (dabs(ja-je).GT.(k1+0.1d0))goto 20
                  if ((ja+je).LT.(k1-0.1d0))    goto 20
                  s=s1*(2*je+1)*Wig0(k1,ja,je)**2
                  delE=Eps1(ka)-Eps2(ka)
                  e=delE+Eps(ke)+Eps(kc)-Eps(kd)-Eps2(kb)
                  if (Khot.EQ.1) then               !### HOT correction
                    s=s*Chot(k1+1)
                  end if
                  if (s.NE.0.d0) then
                    dt=s*Rp2(ke,lja,la,ljb,lb,c1,c1,Pa,Qa,Pb,Qb)/e
                    xs=xs+dt
                    if (ld.LE.Lix-1) xs1=xs1+dt
                    if (ld.LE.Lix-2) xs2=xs2+dt
                    ds=ds-dt/e               !### changed 28/12/96
                  end if
 20             end do
              end do
            end if
          end do
        end do
        D13=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D13 (L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D13 (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D14(dsg)
c     Evaluation of the fourth diagram for Sigma:
c
c       a _______________
c               ___ e __/)
c              /___ d __(
c             (  __ c __/
c              )/___________ b
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      include "hfd.par"
      INCLUDE "conf.par"
      implicit real*8 (a-h,j,m,o-z)
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc1=1,Nsx
          kc=Nso+1-kc1
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photons and 'e' state:
            kp1=Isig(lb-lc)
            kp2=Isig(ld-lc)
            kpe=Isig(la)*kp2
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for the first photon momentum:
            kmi1=Kph(jb-jc,kp1)
            kma1=jb+jc+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for the second photon momentum:
            kmi2=Kph(jd-jc,kp2)
            kma2=jd+jc+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmi1.LE.kma1.AND.kmi2.LE.kma2) then
              do k1=kmi1,kma1,2
                call CFnc(k1,ljb,lb,ljc,lc,Pb,Qb,Pc,Qc,c1)
                do k2=kmi2,kma2,2
                  s1=Isig(k1+k2+1)*(2*jc+1)*(2*jd+1)
     >                 *Wig0(k1,ja,jc)*Wig0(k2,jc,jd)
                  call CFnc(k2,ljd,ld,ljc,lc,Pd,Qd,Pc,Qc,c2)
                  do 20 ke1=1,Nsx
                    ke=Nso+1-ke1
                    le=Ll(ke)
                    lje=Lj(ke)
                    je=lje/2.d0
c     restrictions on the last orbital:
                    if (Isig(le).NE.kpe)          goto 20
                    if (dabs(ja-je).GT.(k2+0.1d0))goto 20
                    if ((ja+je).LT.(k2-0.1d0))    goto 20
                    if (dabs(jd-je).GT.(k1+0.1d0))goto 20
                    if ((jd+je).LT.(k1-0.1d0))    goto 20
                    s=s1*FJ6(ja,jc,k1+0.d0,jd,je,k2+0.d0)
                    if (s.NE.0.d0) then
                      s=s*(2*je+1)*Wig0(k1,jd,je)*Wig0(k2,ja,je)
                      delE=Eps1(ka)-Eps2(ka)
                      e=delE+Eps(ke)+Eps(kc)-Eps(kd)-Eps2(kb)
                      if (Khot.EQ.1) then     !### HOT correction
                        s=s*Chot(k1+1)*Chot(k2+1)
                      end if
                      dt=s*Rp2(ke,ljd,ld,lja,la,c1,c2,Pd,Qd,Pa,Qa)/e
                      xs=xs+dt
                      if (ld.LE.Lix-1) xs1=xs1+dt
                      if (ld.LE.Lix-2) xs2=xs2+dt
                      ds=ds-dt/e            !### changed 28/12/96
                    end if
 20               end do
                end do
              end do
            end if
          end do
        end do
        D14=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D14 (L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D14 (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1ab(dsg)                            !### changed on 07/04/97
c     Evaluation of the subtraction diagrams 'a' and 'b' for Sigma:
c                                                         ###################
c       a ___________________ b   a ___________________ b #both photons have#
c                __ e __(                 )__ e __        #k=0, so there is #
c               /       \)               (/       \       #no SMS correction#
c              (\__ d __/      +          \__ d __/)      ###################
c               )                                 (
c              O c                                 O c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        call CFnc(0,lja,la,ljb,lb,Pa,Qa,Pb,Qb,c1)
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c2)
          do kd1=1,Nsx
            kd=Nso+1-kd1
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
            s=-Qq(kc)*(2*jd+1)      !### changed 14/12/98
            if (Khot.EQ.1) then     !### HOT correction
              s=s*Chot(1)
            end if
            do 20 ke=Nsv,Ns
              le=Ll(ke)
              lje=Lj(ke)
c     restrictions on the last orbital:
              if (le.NE.ld.OR.lje.NE.ljd) goto 20
              e =Eps1(ka)-Eps2(ka)+Eps(kd)-Eps(ke)
              e1=Eps1(ka)+Eps(kd)-Eps(ke)-Eps2(kb)
              dt=s*Rp2(ke,ljd,ld,ljd,ld,c1,c2,Pd,Qd,Pd,Qd)*(1/e+1/e1)
              xs=xs+dt
              if (le.LE.Lix-1) xs1=xs1+dt
              if (le.LE.Lix-2) xs2=xs2+dt
              ds=ds-dt*(1/e+1/e1-2/(e+e1))
 20         end do
          end do
        end do
        D1ab=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1ab(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1ab(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1cd(dsg)                            !### changed on 08/11/06
c     Evaluation of the subtraction diagrams 'c' and 'd' for Sigma:
c
c       a _________________ b  a ________________ b ###################
c              ___ e _(              )_ e ___       #first photon have#
c             (\c     \)     +      (/     c/)      #k=0, so SMS cor- #
c              )/_ d _/              \_ d _\(       #rections sum up  #
c                                                   ###################
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        call CFnc(0,lja,la,ljb,lb,Pa,Qa,Pb,Qb,c1)
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd1=1,Nsx
            kd=Nso+1-kd1
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon:
            kpp=Isig(lc-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(jc-jd,kpp)
            kmax=jc+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s=Qq(kc)*(2*jd+1)*Wig0(k1,jd,jc)**2!### changed 14/12/98
                if (Khot.EQ.1) then     !### HOT correction
                  s=s*Chot(k1+1)
                end if
                call CFnc(k1,ljd,ld,ljc,lc,Pd,Qd,Pc,Qc,c2) !# changed 08/11/06
                do 20 ke=Nsv,Ns
                  le=Ll(ke)
                  lje=Lj(ke)
c     restrictions on the last orbital:
                  if (le.NE.ld.OR.lje.NE.ljd) goto 20
                  e =Eps1(ka)-Eps2(ka)+Eps(kd)-Eps(ke)
                  e1=Eps1(ka)+Eps(kd)-Eps(ke)-Eps2(kb)
                  dt=s*Rp2(ke,ljd,ld,ljc,lc,c1,c2,Pd,Qd,Pc,Qc)
     >              *(1/e+1/e1)
                  xs=xs+dt
                  if (le.LE.Lix-1) xs1=xs1+dt
                  if (le.LE.Lix-2) xs2=xs2+dt
                  ds=ds-dt*(1/e+1/e1-2/(e+e1))
 20             end do
              end do
            end if
          end do
        end do
        D1cd=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1cd(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1cd(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1a1(dsg)                        !# changed 08/11/06
c     Evaluation of the subtraction diagram 'a1' for Sigma:
c
c       a _____________
c              __ e _/(
c             /__ d ___)___ b
c            (
c             )
c            O c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon:
            kpp=Isig(lb-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(jb-jd,kpp)
            kmax=jb+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s=Qq(kc)*(2*jd+1)*Wig0(k1,ja,jd)**2!### changed 14/12/98
                if (Khot.EQ.1) then     !### HOT correction
                  s=s*Chot(k1+1)
                end if
                call CFnc(k1,ljd,ld,ljb,lb,Pd,Qd,Pb,Qb,c1)
                do 20 ke1=1,Nsx
                  ke=Nso+1-ke1
                  le=Ll(ke)
                  lje=Lj(ke)
c     restrictions on the last orbital:
                  if (le.NE.ld.OR.lje.NE.ljd) goto 20
                  e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps(kd)
                  dt=s*Rp2(ke,lja,la,ljd,ld,c1,c2,Pa,Qa,Pd,Qd)/e
                  xs=xs+dt
                  if (ld.LE.Lix-1) xs1=xs1+dt
                  if (ld.LE.Lix-2) xs2=xs2+dt
                  ds=ds-dt/e                    !### changed 28/12/96
 20             end do
              end do
            end if
          end do
        end do
        D1a1=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1a1(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1a1(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1b1(dsg)                        !# changed 08/11/06
c     Evaluation of the subtraction diagram 'b1' for Sigma:
c                     O c
c                    (
c       a _______ d __)
c             )__ e _/
c            (/____________ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon:
            kpp=Isig(la-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(ja-jd,kpp)
            kmax=ja+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s=Qq(kc)*(2*jd+1)*Wig0(k1,ja,jd)**2!### changed 14/12/98
                if (Khot.EQ.1) then     !### HOT correction
                  s=s*Chot(k1+1)
                end if
                call CFnc(k1,ljd,ld,lja,la,Pd,Qd,Pa,Qa,c1)
                do 20 ke1=1,Nsx
                  ke=Nso+1-ke1
                  le=Ll(ke)
                  lje=Lj(ke)
c     restrictions on the last orbital:
                  if (le.NE.ld.OR.lje.NE.ljd) goto 20
                  e=Eps1(ka)-Eps2(kb)+Eps(ke)-Eps(kd)
                  dt=s*Rp2(ke,ljb,lb,ljd,ld,c1,c2,Pb,Qb,Pd,Qd)/e
                  xs=xs+dt
                  if (ld.LE.Lix-1) xs1=xs1+dt
                  if (ld.LE.Lix-2) xs2=xs2+dt
                  ds=ds-dt/e                    !### changed 28/12/96
 20             end do
              end do
            end if
          end do
        end do
        D1b1=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1b1(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1b1(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1c1(dsg)                        !# changed 08/11/06
c     Evaluation of the subtraction diagram 'c' for Sigma:
c
c       a _______________
c                __ e _/(
c              c/)_ d ___)___ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photons:
            kp1=Isig(la-ld)
            kp2=Isig(lc-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmi1=Kph(ja-jd,kp1)
            kma1=ja+jd+0.1d0
            kmi2=Kph(jc-jd,kp2)
            kma2=jc+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmi1.LE.kma1.AND.kmi2.LE.kma2) then
              do k1=kmi1,kma1,2
                s=-Qq(kc)*(2*jd+1)*Wig0(k1,ja,jd)**2!### changed 14/12/98
                if (Khot.EQ.1) then     !### HOT correction
                   s=s*Chot(k1+1)
                end if
                call CFnc(k1,ljd,ld,ljb,lb,Pd,Qd,Pb,Qb,c1)
                do k2=kmi2,kma2,2
                  s1=s*Wig0(k2,jc,jd)**2
                  call CFnc(k2,ljd,ld,ljc,lc,Pd,Qd,Pc,Qc,c2)
                  do 20 ke1=1,Nsx
                    ke=Nso+1-ke1
                    le=Ll(ke)
                    lje=Lj(ke)
c     restrictions on the last orbital:
                    if (le.NE.ld.OR.lje.NE.ljd) goto 20
                    e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps(kd)
                    dt=s1*Rp2(ke,lja,la,ljc,lc,c1,c2,Pa,Qa,Pc,Qc)/e
                    xs=xs+dt
                    if (ld.LE.Lix-1) xs1=xs1+dt
                    if (ld.LE.Lix-2) xs2=xs2+dt
                    ds=ds-dt/e                  !### changed 28/12/96
 20               end do
                end do
              end do
            end if
          end do
        end do
        D1c1=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1c1(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1c1(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1d1(dsg)                        !# changed 08/11/06
c     Evaluation of the subtraction diagram 'd' for Sigma:
c
c       a _________ d ___ c
c               )__ e _(/
c              (/____________ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /Khot/Khot/Chot/Chot(10)/Lix/Lix/Dlix/Dlx1,Dlx2
     2        /QNl/la,lb,lc,ld
     3        /QNj/ja,jb,jc,jd
     4        /Orbit/ka,kb,kc,kd
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        xs1=0.d0
        xs2=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Ns
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photons:
            kp1=Isig(la-ld)
            kp2=Isig(lc-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmi1=Kph(ja-jd,kp1)
            kma1=ja+jd+0.1d0
            kmi2=Kph(jc-jd,kp2)
            kma2=jc+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmi1.LE.kma1.AND.kmi2.LE.kma2) then
              do k1=kmi1,kma1,2
                s=-Qq(kc)*(2*jd+1)*Wig0(k1,ja,jd)**2!### changed 14/12/98
                if (Khot.EQ.1) then     !### HOT correction
                  s=s*Chot(k1+1)
                end if
                call CFnc(k1,ljd,ld,lja,la,Pd,Qd,Pa,Qa,c1)
                do k2=kmi2,kma2,2
                  s1=s*Wig0(k2,jc,jd)**2
                  call CFnc(k2,ljd,ld,ljc,lc,Pd,Qd,Pc,Qc,c2)
                  do 20 ke1=1,Nsx
                    ke=Nso+1-ke1
                    le=Ll(ke)
                    lje=Lj(ke)
c     restrictions on the last orbital:
                    if (le.NE.ld.OR.lje.NE.ljd) goto 20
                    e=Eps1(ka)-Eps2(kb)+Eps(ke)-Eps(kd)
                    dt=s1*Rp2(ke,ljb,lb,ljc,lc,c1,c2,Pb,Qb,Pc,Qc)/e
                    xs=xs+dt
                    if (ld.LE.Lix-1) xs1=xs1+dt
                    if (ld.LE.Lix-2) xs2=xs2+dt
                    ds=ds-dt/e
 20               end do
                end do
              end do
            end if
          end do
        end do
        D1d1=xs
        Dlx1=Dlx1+xs1
        Dlx2=Dlx2+xs2
        dsg=dsg+ds
        if (Kout.GT.1) then
          write(6, 25) Lix-2,xs2,Lix-1,xs1,Lix,xs
          write(11,25) Lix-2,xs2,Lix-1,xs1,Lix,xs
        end if
 25     format(4X,'D1d1(L_max=',I2,')  =',E12.5)
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1d1(E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1e(dsg)
c     Evaluation of the subtraction diagram 'e' for Sigma:
c
c                        O c
c                       (
c       a _______________)
c                __ e __/
c               /____________ b
c              (
c               )
c              O d
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /QNl/la,lb,lc,ld
     2        /QNj/ja,jb,jc,jd
     3        /Orbit/ka,kb,kc,kd/Dlix/Dlx1,Dlx2
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c1)
          do kd=Nsv,Nsh
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd/2.d0
            call ReadFF (12,kd+4,Pd,Qd,2)
            call CFnc(0,ljd,ld,ljd,ld,Pd,Qd,Pd,Qd,c2)
            s=-Qq(kc)*Qq(kd)   !### changed 14/12/98
            do 20 ke1=1,Nsx
              ke=Nso+1-ke1
              le=Ll(ke)
              lje=Lj(ke)
c     restrictions on the last orbital:
              if (le.NE.la.OR.lje.NE.lja) goto 20
              e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps2(kb)
              dt=s*Rp2(ke,lja,la,ljb,lb,c1,c2,Pa,Qa,Pb,Qb)/e
              xs=xs+dt
              ds=ds-dt/e                        !### changed 28/12/96
 20         end do
          end do
        end do
        D1e=xs
        Dlx1=Dlx1+xs                            !### This diagram does
        Dlx2=Dlx2+xs                            !#### not depend on Lix
        dsg=dsg+ds
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1e (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1f(dsg)                         !# changed 08/11/06
c     Evaluation of the subtraction diagram 'f' for Sigma:
c
c                        O c
c                       (
c       a _______________)
c                __ e __/
c             d /)___________ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /QNl/la,lb,lc,ld
     2        /QNj/ja,jb,jc,jd
     3        /Orbit/ka,kb,kc,kd/Dlix/Dlx1,Dlx2
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c1)
          do kd=Nsv,Nsh
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd*0.5d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon:
            kpp=Isig(la-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(ja-jd,kpp)
            kmax=ja+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
              do k1=kmin,kmax,2
                s=Qq(kc)*Qq(kd)*Wig0(k1,jd,ja)**2  !### changed 14/12/98
                call CFnc(k1,ljb,lb,ljd,ld,Pb,Qb,Pd,Qd,c2)
                do 20 ke1=1,Nsx
                  ke=Nso+1-ke1
                  le=Ll(ke)
                  lje=Lj(ke)
c     restrictions on the last orbital:
                  if (le.NE.la.OR.lje.NE.lja) goto 20
                  e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps2(kb)
                  dt=s*Rp2(ke,lja,la,ljd,ld,c1,c2,Pa,Qa,Pd,Qd)/e
                  xs=xs+dt
                  ds=ds-dt/e                    !### changed 28/12/96
 20             end do
              end do
            end if
          end do
        end do
        D1f=xs
        Dlx1=Dlx1+xs                            !### This diagram does
        Dlx2=Dlx2+xs                            !#### not depend on Lix
        dsg=dsg+ds
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1f (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1g(dsg)                         !# changed 08/11/06
c     Evaluation of the subtraction diagram 'f' for Sigma:
c
c       a _______________ d
c                __ e _(/
c               /____________ b
c              (
c               )
c              O c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /QNl/la,lb,lc,ld
     2        /QNj/ja,jb,jc,jd
     3        /Orbit/ka,kb,kc,kd/Dlix/Dlx1,Dlx2
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc/2.d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          call CFnc(0,ljc,lc,ljc,lc,Pc,Qc,Pc,Qc,c1)
          do kd=Nsv,Nsh
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd*0.5d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photon:
            kpp=Isig(la-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmin=Kph(ja-jd,kpp)
            kmax=ja+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmin.LE.kmax) then
               do k1=kmin,kmax,2
                 s=Qq(kc)*Qq(kd)*Wig0(k1,jd,ja)**2  !### changed 14/12/98
                 call CFnc(k1,lja,la,ljd,ld,Pa,Qa,Pd,Qd,c2)
                 do 20 ke1=1,Nsx
                   ke=Nso+1-ke1
                   le=Ll(ke)
                   lje=Lj(ke)
c     restrictions on the last orbital:
                   if (le.NE.la.OR.lje.NE.lja) goto 20
                   e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps2(kb)
                   dt=s*Rp2(ke,ljb,lb,ljd,ld,c1,c2,Pb,Qb,Pd,Qd)/e
                   xs=xs+dt
                   ds=ds-dt/e                   !### changed 28/12/96
 20              end do
               end do
             end if
           end do
        end do
        D1g=xs
        Dlx1=Dlx1+xs                            !### This diagram does
        Dlx2=Dlx2+xs                            !#### not depend on Lix
        dsg=dsg+ds
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1g (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
C     =================================================
      Function D1h(dsg)
c     Evaluation of the subtraction diagram 'f' for Sigma:
c
c       a _______________
c                __ e _(/ c
c             d /)___________ b
c
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit real*8 (a-h,j,m,o-z)
      include "hfd.par"
      INCLUDE "conf.par"
       common /Ns/Ns/Nso/Nso/Nsh/Nsh/Nsx/Nsx/Kout/Kout/Nsv/Nsv
     1        /QNl/la,lb,lc,ld
     2        /QNj/ja,jb,jc,jd
     3        /Orbit/ka,kb,kc,kd/Dlix/Dlx1,Dlx2
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Qq/Qq(IPs)
     >        /Eps/Eps(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
     >        /C1/C1(IP6)/C2/C2(IP6)
     >        /Porb/Pa(IP6),Pb(IP6),Pc(IP6),Pd(IP6)
     >        /Qorb/Qa(IP6),Qb(IP6),Qc(IP6),Qd(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        xs=0.d0
        ds=0.d0
        lja=ja+ja+0.1d0
        ljb=jb+jb+0.1d0
        do kc=Nsv,Nsh
          lc=Ll(kc)
          ljc=Lj(kc)
          jc=ljc*0.5d0
          call ReadFF (12,kc+4,Pc,Qc,2)
          do kd=Nsv,Nsh
            ld=Ll(kd)
            ljd=Lj(kd)
            jd=ljd*0.5d0
            call ReadFF (12,kd+4,Pd,Qd,2)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     parity of photons:
            kp1=Isig(la-lc)
            kp2=Isig(la-ld)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     limits for photons momenta:
            kmi1=Kph(ja-jc,kp1)
            kma1=ja+jc+0.1d0
            kmi2=Kph(ja-jd,kp2)
            kma2=ja+jd+0.1d0
C     - - - - - - - - - - - - - - - - - - - - - - - - -
c     sums over photon momenta:
            if (kmi1.LE.kma1.AND.kmi2.LE.kma2) then
              do k1=kmi1,kma1,2
                s=-Qq(kc)*Qq(kd)*Wig0(k1,ja,jc)**2  !### changed 14/12/98
                call CFnc(k1,lja,la,ljc,lc,Pa,Qa,Pc,Qc,c1)
                do k2=kmi2,kma2,2
                  s1=s*Wig0(k2,ja,jd)**2
                  call CFnc(k2,ljb,lb,ljd,ld,Pb,Qb,Pd,Qd,c2)
                  do 20 ke1=1,Nsx
                    ke=Nso+1-ke1
                    le=Ll(ke)
                    lje=Lj(ke)
c     restrictions on the last orbital:
                    if (le.NE.la.OR.lje.NE.lja) goto 20
                    e=Eps1(ka)-Eps2(ka)+Eps(ke)-Eps2(kb)
                    dt=s1*Rp2(ke,ljc,lc,ljd,ld,c1,c2,Pc,Qc,Pd,Qd)/e
                    xs=xs+dt
                    ds=ds-dt/e                  !### changed 28/12/96
 20               end do
                end do
              end do
            end if
          end do
        end do
        D1h=xs
        Dlx1=Dlx1+xs                            !### This diagram does
        Dlx2=Dlx2+xs                            !#### not depend on Lix
        dsg=dsg+ds
        if (Kout.GT.0) then
          write(6, 15) Eps1(ka),xs,ds
          write(11,15) Eps1(ka),xs,ds
        end if
 15     format(1X,'>> D1h (E=',F7.3,') =',E12.5,'a.u. d_e =',E12.5)
       return
      end
