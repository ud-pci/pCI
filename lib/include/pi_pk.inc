c     ============= 19/09/03 ==== last update: 28/10/05
c     Effecitve radial integrals for SMS (pi dot pk)
c     =================================================
      Function P_eff (ja,la,jc,lc, !### effective vertex for (pi dot pk)
     >                Pa,Qa,Pc,Qc) !### see 19/09/03 and 16/09/05
      implicit real*8 (a-h,o-z)    !### for derivation details
      INCLUDE "hfd.par"            
      INCLUDE "phys.par"
       common /Kt/Kt/Ii/Ii/klow/klow/Z/Z
       common /C/C(IP6)/R/R(IP6)
       dimension Pa(IP6),Pc(IP6),Qa(IP6),Qc(IP6),Da(IP6)

c        write(*,*) ' kt=',kt,' klow=',klow,' ii=',ii
c        read(*,*) 

        P_eff=0.d0
        if (iabs(la-lc).NE.1) return

        ila=ja-la ! l for low component
        ilc=jc-lc ! "      "      "
        lmax=max(la,lc)
        ilmax=max(ila,ilc)

        ih=2-Kt
        do i=1,Ii,ih
          C(i)=0.d0
        end do

        ga=Pa(ii+4)               !### upper component
        C(ii+4)=Pc(ii+4)+ga-1
        call Dif(Pa,Da,ga)

        do i=1,Ii,ih
          C(i)=Pc(i)*(Da(i)+(la-lc)*lmax/R(i)*Pa(i))
        end do

        if (klow.GT.0) then       !### lower component
          call Dif(Qa,Da,ga)
          do i=1,Ii,ih
            C(i)=C(i)+Qc(i)*(Da(i)+(ila-ilc)*ilmax/R(i)*Qa(i))
          end do
        end if

        if (klow.GT.1) then       !### relativistic correction
          half=0.5d0
          xjc=half*jc
          xja=half*ja
          xjx=dmax1(xja,xjc)
          xjn=dmin1(xja,xjc)
          xla=la
          xlc=lc
          yla=ila
          ylc=ilc
          s0=dsqrt((xjx-xjn+1.d0)/(2*xja+1.d0))
     >      /FJ3(xjc,xja,1.d0,half,-half,0.d0)

          if (dabs(s0).GT.1.d9) then
            write(*,*) ' s0 = infty'
            write(*,*) ja,jc
            write(*,*) la,lc
            read(*,*)
          end if

          if (mod(la+1,2).NE.0) s0=-s0
          s1=0.d0
          if (lc.EQ.ila) s1=dsqrt((xjc+3*xja-2*ila+1)/(2*ila+1))
          ss=0.5d0*Z/DPcl*(s0*s1+1.0d0)
          do i=1,Ii,ih
            C(i)=C(i)+ss/R(i)*Pc(i)*Qa(i)
          end do

          s2=0.d0
          if (ilc.EQ.la) s2=dsqrt((xjc+3*xja-2*la+1)/(2*la+1))
          ss=0.5d0*Z/DPcl*(s0*s2-1.0d0)
          do i=1,Ii,ih
            C(i)=C(i)+ss/R(i)*Qc(i)*Pa(i)
          end do
        end if

        call Sint1(ds)
        P_eff=ds
       return
      end
c     =================================================
      Function SMS_core(na,nb, !### exchange core potential for (pi_pk)
     >                  kan)   !### - file with orbitals and derivatives
      INCLUDE "conf.par"       !### <b|V_SMS|a> =
      include "hfd.par"        !### sum_c Gaunt^1_{a,c,c,b} R^1_{a,c,c,b}
       implicit real*8 (a-h,o-z)
       common /Nso/Nso
     >        /Nn/Nn(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /UP/Pa(IP6),Pc(IP6),Qa(IP6),Qc(IP6)
       common /kout/kout
        character*1 let(9)
        data let /'s','p','d','f','g','h','i','k','l'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -

        a_ab=0.d0
        if (Nso.EQ.0) goto 1000

        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)


        if (ja.NE.jb.OR.la.NE.lb) goto 1000

        z12=0.5d0
        z00=0.0d0
        xja=z12*ja
        do 10 nc=1,Nso
          jc=Jj(nc)
          lc=Ll(nc)
          if (iabs(ja-jc).GT.2.OR.iabs(la-lc).NE.1) goto 10
          xjc=z12*jc
          f0=(jc+1)*FJ3(xja,xjc,1.d0,z12,-z12,z00)**2
          call ReadF(kan,nc+4,Pc,Qc,2)
          call ReadF(kan,nb+4,Pa,Qa,2)
          s1=+f0*P_eff(jc,lc,jb,lb,Pc,Qc,Pa,Qa)
          call ReadF(kan,na+4,Pa,Qa,2)
          s1=-s1*P_eff(jc,lc,ja,la,Pc,Qc,Pa,Qa)
          if (kout.GE.2) write(*,5) Nn(nc),let(lc+1),jc,s1
 5        format(4X,'contribution of ',I2,A1,I1,'/2 =',E12.5)
          a_ab=a_ab+s1
 10     end do

        if (kout.GE.1) then
          write( *,15) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
          write(11,15) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
 15       format(/2X,'<',I2,A1,I1,'/2| V_SMS^core |',
     >           I2,A1,I1,'/2> = ',E12.5)
        end if

 1000   SMS_core = a_ab
       return
      end
c     =================================================
      Function V_nms (na,nb,   !### Normal mass shift
     >                  kan)   !### - file with orbitals and derivatives
      implicit real*8 (a-h,o-z)
      INCLUDE "conf.par"
      INCLUDE "phys.par"
      include "hfd.par"
       common /kt/kt/ii/ii/k_is/k_is/klow/klow/Z/Z
       common /kout/kout
     >        /Nn/Nn(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /C/C(IP6)/R/R(IP6)
     >        /UP/Pa(IP6),Pb(IP6),P1a(IP6),P1b(IP6)
     >        /DOWN/Qa(IP6),Qb(IP6),Q1a(IP6),Q1b(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        a_ab=0.d0
        if (k_is.LE.2) goto 1000 ! NMS is calculated for k_is=3,4

        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)

        if (ja.NE.jb.OR.la.NE.lb) goto 1000

        ih=2-kt
        call ReadF(kan,na+4,Pa,Qa,2)
        call ReadF(kan,nb+4,Pb,Qb,2)
        ga=Pa(ii+4)
        call Dif(Pa,P1a,ga)
        call Dif(Pb,P1b,ga)
        ip=la*(la+1)
        cis2=0.5d0                !  No c_is here
        if (klow.GT.0) then
          call Dif(Qa,Q1a,ga)
          call Dif(Qb,Q1b,ga)
          ila=ja-la
          iq=ila*(ila+1)
          xja=0.5d0*ja
          clj=-xja*(xja+1)+la*(la+1)+0.75d0 - 2.d0
          aZ=Z/DPcl
        end if

        do i=1,ii,ih
          ri=R(i)
          ri2=ri*ri
          C(i)=P1b(i)*P1a(i)+ip*Pb(i)*Pa(i)/ri2
          if (klow.GT.0) then
            C(i)=C(i)+Q1b(i)*Q1a(i)+iq*Qb(i)*Qa(i)/ri2
          end if
          if (klow.GE.2) then
            xi = 2/ri*(P1b(i)*Qa(i)+Qb(i)*P1a(i))
     >         + clj/ri2*(Pb(i)*Qa(i)+ Qb(i)*Pa(i))
            C(i)=C(i)+aZ*xi
          end if
          C(i)=cis2*C(i)
        end do
        C(ii+4)=2*Pa(ii+4)-2
        call Sint1(a_ab)

 1000   V_nms = a_ab
       return
      end
