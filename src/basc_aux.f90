Module basc_aux

    Use basc_variables

    Implicit None

    Public

  Contains

    Subroutine recunit
        ! This Subroutine determines the record unit
        Implicit None

        Integer          :: lrec, iflag, nbytes
        Character(Len=8) :: d1, t1, d2, t2

        t1='abcdefgh'
        d1='        '
        t2='hgfedcba'
        d2='        '
        lrec=0
        iflag=1
200     lrec=lrec+1
        If (lrec.gt.8) Then
          Write(*,*)  'lrec > 8'
          Stop
        End If
        Open(unit=13,file='test.tmp',status='unknown',access='direct',recl=lrec)
        Write(13,rec=1,err=210) t1
        Write(13,rec=2,err=210) t2
        Read(13,rec=1,err=210) d1
        Read(13,rec=2,err=210) d2
        If (d1.ne.t1) goto 210
        If (d2.ne.t2) goto 210
        iflag=0
210     Close(unit=13,status='delete')
        If (iflag.ne.0) goto 200
        nbytes=8/lrec
        ipmr=4/nbytes

        Return
    End Subroutine recunit

    Subroutine Input
        Use conf_init, Only : ReadConfInp, ReadConfigurations
        Implicit None

        Integer           :: Ncci, key
        Character(Len=1)  :: NAME(16),str2(2)*3,str4*5
        Character(Len=32) :: strfmt
        data str2 /' NO','YES'/

        strfmt = '(4X,"PROGRAM BasC")'
        Write( *,strfmt)
        Write(11,strfmt)
        ! - input from the file 'CONF.INP' - - - - - - - - - - - - - - - -
        Call ReadConfInp

        If (Ncpt.GT.Nc) Then
            Write (*,*) ' NcPT =',Ncpt, ' & Nc = ',Nc
            Write (*,*) ' which to use for integrals (1-NcPT,2-Nc)?'
            Read (*,*) key
            If (key.EQ.1) Then
                Write( *,*) 'NcPT =',Ncpt,' used for integrals'
                Write(11,*) 'NcPT =',Ncpt,' used for integrals'
                Ncci=Nc
                Nc=Ncpt
            End If
        End If

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Call ReadConfigurations
        Return
       
    End Subroutine Input

    Subroutine Init(Norb)
        Use readfff
        Implicit None

        Integer :: i, i0, ic, imax, j, n, n0, ni, nmin, nsmax, nsmax1, If
        Integer :: Norb, nsb, kkj, llj, jjj, nj, ng, nnj, err_stat
        Real(dp) :: c1, c2, r1, a2, z1, d
        Integer, Dimension(4*IPs) :: IQN
        Real(dp), Dimension(IP6) ::  P, Q, P1, Q1
        Real(dp), Dimension(4*IP6) :: PQ
        Real(dp), Dimension(IPs)    :: Qw, Qq1
        logical :: longbasis
        Character(Len=512) :: strfmt
        Character(Len=256) :: err_msg


        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),(P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))

        ! small number:
        c1=0.01d0
        ! speed of light is taken from "phys.par":
        Cl=DPcl
        Mj=2*dabs(Jm)+0.01d0

        Open(12,file='HFD.DAT',access='DIRECT',status='OLD',recl=2*IP6*IPmr,iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            strfmt='(/2X,"file HFD.DAT is absent"/)'
            Write( *,strfmt)
            Write(11,strfmt)
            Stop
        End If

        Call ReadF (12,1,P,Q,2)
        Call ReadF (12,2,R,V,2)
        Call ReadF (12,3,P1,Q1,2)

        z1  =PQ(1)
        If (dABS(Z-z1).GT.1.d-6) Then
            strfmt='(/2X,"nuc. charge is changed"/2X,"Z1=",F12.6/2X,"Z2=",F12.6)'
            Write( *,strfmt) Z,Z1
            Write(11,strfmt) Z,Z1
            Stop
        End If
        Ns  =PQ(2)+C1
        II  =PQ(3)+C1
        R1  =PQ(4)
        R2  =dabs(PQ(5))
        H   =PQ(6)
        H0  =H
        Bt  =PQ(7)
        Al  =PQ(8)
        Kt  =PQ(9)+C1
        Ng  =PQ(10)+C1
        Rnuc=PQ(13)

        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))

        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
        strfmt = '(/4X,"Z   = ",F6.2,5X,"Kt  =",I3,  7X,"II =",I4,    &
                /4X,"H   =",F7.4, 5X,"R2  =",F6.2,4X,"R1 =",E11.4, &
                /4X,"Rnuc=",E11.4,1X,"Al  =",F7.4,3X,"Bt =",F5.2,  &
                /4X,"Nsp =",I5,   7X,"Ns  =",I3,  7X,"Nso=",I3,4X,"Nc =",I4)'
        Write( *,strfmt) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc
        Write(11,strfmt) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc

        If (longbasis) Then
            Write( *,*) ' Using variant for long basis '
            Write(11,*) ' Using variant for long basis '
            Do ni=1,Ns
                Nn(ni)=IQN(4*ni-3)
                Ll(ni)=IQN(4*ni-2)
                Kk(ni)=IQN(4*ni-1)
                Jj(ni)=IQN(4*ni)
                Qq(ni)=Qq1(ni)
            End Do
        Else
            If=20
            Do ni=1,Ns
                If=If+1
                Nn(ni)=PQ(If)+c1
                If=If+1
                Ll(ni)=PQ(If)+c1
                If=If+1
                Qq(ni)=PQ(If)
                If=If+2
                c2=dsign(c1,PQ(If))
                Kk(ni)=PQ(If)+c2
                If=If+1
                c2=dsign(c1,PQ(If))
                Jj(ni)=PQ(If)+c2
                Qw(ni)=0.d0
            End Do
        End If

        Norb=0
        Kfile=0
        nsb=Ns
        Do nj=1,Nsp
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
            If (nj.LE.Nso) Then
                Do ni = 1,Nso
                    If (Nn(ni).EQ.Nn(nj).AND.Ll(ni).EQ.Ll(nj)) Qw(ni) = Qw(ni) + Nq(nj)   ! number of e on NR shell
                End Do
            End If
            Do i=1,nsb
                ni=i
                If (nnj.EQ.Nn(ni).AND.Kk(ni).EQ.kkj) goto 210
            End Do
            Norb=Norb+1
            ni=Ns+Norb
            nsb=Ns+Norb 
            Nn(ni)=nnj
            Kk(ni)=kkj
            Ll(ni)=llj
            Jj(ni)=jjj
            If (Kbas(ni).EQ.3) Kfile=1
            If (Kbas(ni).EQ.4) kfile=1
 210        Nip(nj)=ni
        End Do
        Write(*,*) Norb, ' new orbitals to be formed'
        nsmax=(4*IP6-20)/6
        nsmax1=(4*IP6-20)/3
        If (nsb.GT.nsmax.AND..NOT.longbasis) Then
            Write(*,*) ' Ns =',nsb,' > ',nsmax
            Write(*,*) ' switch to long basis variant '
        End If
        If (nsb.GT.nsmax1) Then
            Write(*,*) ' For IP6 =',IP6
            Write(*,*) ' maximum length of basis set is ',nsmax1
            Stop
        End If
        If (Nso.NE.0) Then
            Do ni=1,Nso
                Qq(ni) = Qw(ni)*(Jj(ni)+1.d0)/(4*Ll(ni)+2)    !# N_e on a shell
            End Do
        End If
        n0=0
        ne=0
        nmin=Nso+1
        If (nmin.LE.Nsp) Then
            Do ni=nmin,Nsp
                n0=n0+Nq(ni)
            End Do
            ne=n0/Nc
        End If
        Nst=0
        Do ni=1,Ns
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                Nst=Nst+1
            End Do
        End Do
        strfmt = '(4X,"Ne  =",I3,7X,"Nst =",I4)'
        Write( *,strfmt) ne,Nst
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        Do ni=nmin,Nsp
            i=i+1
            n=n+Nq(ni)
            If (n.LT.ne) Cycle
            ic=ic+1
            If (n.GT.ne) Then
                strfmt='(/2X,"wrong number of electrons"/2X,"for configuration ICONF =",I4/)'
                Write( *,strfmt) IC
                Write(11,strfmt) IC
                Stop
            End If
            Nvc(ic)=i
            Nc0(ic)=Nso+i0
            i0=i0+i
            n=0
            i=0
        End Do
        Open(13,file=FNAME,status='UNKNOWN',access='DIRECT',recl=2*IP6*IPmr)
        Do ni=1,4
            Call ReadF (12,ni,P,Q,2)
            Call WriteF(13,ni,P,Q,2)
        End Do
        Do ni=1,Ns
            Call ReadF (12,ni+4,P,Q,2)
            Call WriteF(13,ni+4,P,Q,2)
            Call ReadF (12,ni+4+Ns,P,Q,2)
            Call WriteF(13,ni+4+Ns+Norb,P,Q,2)
        End Do
        Close(13)
        Close(12)
        Return

    End Subroutine Init

    Subroutine Fill_N_l         
        ! This subroutine finds n_max for each l
        Implicit None
    
        Integer :: i, l_max, n_i, n2, l_i, l1, l
        Real(dp) :: x
        Character(Len=512) :: strfmt

        N_l=0
        l_max=0

        Do i=1,Nsp
            x=dabs(QNL(i))+1.d-9
            n_i=10*x               ! principle QN
            n2=100*x               
            l_i=n2-10*n_i          ! orbital QN  
            l1=l_i+1
            If (n_i.GT.N_l(l1)) N_l(l1)=n_i          ! - max n for given l
            If (l_i.GT.l_max) l_max=l_i
        End Do

        strfmt='(/4x,"Max principle QN for given l:")'
        Write ( *,strfmt)
        Write (11,strfmt)

        Do l=0,l_max
            l1=l+1
            strfmt='(4x,"n_max(",i2,") = ",i2)'
            Write ( *,strfmt) l,N_l(l1)
            Write (11,strfmt) l,N_l(l1)
        End Do

        Return  
    End Subroutine Fill_N_l

    Subroutine Core
        Use wigner
        Use readfff
        Use sintg
        Use test_ori
        Implicit None

        Integer :: ih, i, nigd, irr, ni, na, la, ja, ir, ir1, ir2, ir3, ir4, igm, m, m1, im
        Integer :: ik, nj, k, nb, lb, jb, k1, ip, kmin, kmax
        Real(dp) :: hh, Hcore, err1, err2, qa, t, r1, gm, dgm, cpp, cqq, qb, s, d, xk, xja, xjb, ds, yy
        Real(dp), Dimension(IP6) :: A, B, CP, CQ, CC, Y
        Real(dp), Dimension(20) :: dd, ss
        Character(Len=512) :: strfmt

        Open(12,file=FNAME,status='OLD',access='DIRECT',recl=2*IP6*IPmr)
        Ecore=0.d0
        Hcore=0.d0
        ih=2-kt
        hh=h0*ih/3.d0
        C=0.d0
        Y=0.d0

        !### Coulomb potential of the core
        If (Nso.NE.0) Then
            Call Y0(Y)               
            err1=2.d-4
            err2=2.d-3
            nigd=0
            irr=0
        End If

        Do ni=1,Ns               !### one-electron part of the Hamiltonian
            na=Nn(ni)
            la=Ll(ni)
            ja=Jj(ni)
            qa=Qq(ni)
            t=Cl*Kk(ni)
            Call ReadF (12,ni+4,P,Q,2)
            Call ReadF (12,ni+Ns+4,A,B,2)
            If (Kout.GT.0) Then
                Write(11,*) ' Core: forming CP and CQ for orbital ',ni
                Write(*,*) ' Core: Testing P,Q,A,B for orb. ',ni,'..'
            End If
            Call Test_Origin('Core:    P',P,1.d-6,ir1)
            Call Test_Origin('Core:    Q',Q,1.d-6,ir2)
            Call Test_Origin('Core:    A',A,err1,ir3)
            Call Test_Origin('Core:    B',B,err1,ir4)
            ir=ir1+ir2+ir3+ir4
            Do i=1,ii,ih
                d= Cl*B(i)-(Z*P(i)+t*Q(i))/R(i)
                s=-Cl*A(i)-(Z*Q(i)+t*P(i))/R(i)-2*Cl*Cl*Q(i)
                C(i)=d*P(i)+s*Q(i)
                CP(i)=d+Y(i)*P(i)
                CQ(i)=s+Y(i)*Q(i)
            End Do
            r1=R(1)
            gm=P(ii+4)
            CP(ii+4)=gm-1.d0
            CQ(ii+4)=gm-1.d0
            igm=gm+1.d-5
            dgm=gm-igm                   !### dgm=0 for finite nucleus
            Do m=0,MaxT
                m1=m+1
                im=ii+5+m
                dd(m1)=Cl*B(im)-t*Q(im)
                ss(m1)=-Cl*A(im)-t*P(im)
                If (m.GE.1) ss(m1)=ss(m1)-2*Cl*Cl*r1*Q(im-1)
                If (dgm.GT.1.d-6) Then     !### case of the pointlike nucleus
                    dd(m1)=dd(m1)-Z*P(im)
                    ss(m1)=ss(m1)-Z*Q(im)
                Else                       !### case of the finite nucleus
                    If (m.GE.1) Then
                        dd(m1)=dd(m1)-1.5d0*Z*P(im-1)
                        ss(m1)=ss(m1)-1.5d0*Z*Q(im-1)
                    End If
                    If (m.GE.3) Then
                        dd(m1)=dd(m1)+0.5d0*Z*P(im-3)
                        ss(m1)=ss(m1)+0.5d0*Z*Q(im-3)
                    End If
                End If
                cpp=0.d0
                cqq=0.d0
                If (m.GE.1) Then
                    Do k=0,m-1
                        ik=ii+5+k
                        cpp=cpp+P(ik)*Y(im-k-1)
                        cqq=cqq+Q(ik)*Y(im-k-1)
                    End Do
                End If
                CP(im)=dd(m1)+r1*cpp
                CQ(im)=ss(m1)+r1*cqq
            End Do
    
            If (Nso.EQ.0) Then
                ! test of expansion at the origin
                Call WriteF(12,ni+4+Ns,CP,CQ,2)
                Cycle
            End If

            If (ni.LE.Nso) Then
                C(ii+4)=2*P(ii+4)-1.d0
                Call Sint1(ds)
                Hcore=Hcore+qa*ds
                Ecore=Ecore+qa*ds
                Do i=1,ii,ih
                    C(i)=Y(i)*(P(i)**2+Q(i)**2)
                End Do
                C(ii+4)=2*P(ii+4)
                Call Sint1(ds)
                Ecore=Ecore+0.5d0*qa*ds
            End If

            Do nj=1,Nso
                Call ReadF (12,nj+4,A,B,2)
                nb=Nn(nj)
                lb=Ll(nj)
                jb=Jj(nj)
                qb=Qq(nj)
                kmin=iabs(ja-jb)/2+1
                kmax=(ja+jb)/2+1
                Do k1=kmin,kmax
                    k=k1-1
                    ip=la+lb+k
                    If (ip.EQ.2*(ip/2)) Then
                        Do i=1,ii,ih
                            C(i)=P(i)*A(i)+Q(i)*B(i)
                        End Do
                        C(ii+4)=P(ii+4)+A(ii+4)
                        Call Yk(k)
                
                        xja=ja/2.d0
                        xjb=jb/2.d0
                        xk=k
                        d=qb*(FJ3(xk,xja,xjb,0.d0,-0.5d0,0.5d0))**2
                        Do i=1,ii,ih
                            yy=d*C(i)/(V(i)*hh*R(i))
                            CP(i)=CP(i)-yy*A(i)
                            CQ(i)=CQ(i)-yy*B(i)
                        End Do
                
                        If (nj.LE.ni.AND.ni.LE.Nso) Then
                            t=qa
                            If (ni.EQ.nj.AND.k.NE.0) t=0.5d0*(qa-1.d0)*(ja+1.d0)/ja
                            If (ni.EQ.nj.AND.k.EQ.0) t=0.5d0*qa
                            d=t*d
                            Do i=1,ii,ih
                                C(i)=C(i)*(P(i)*A(i)+Q(i)*B(i))/R(i)
                            End Do
                            C(ii+4)=P(ii+4)+A(ii+4)+k
                            Call Sint(ds)  !# integration over ro!
                            Ecore=Ecore-ds*d
                        End If
                    End If
                End Do
            End Do
    
            If (ir.GT.0) Then
                irr=irr+1
                strfmt='(" Orbital",I4,": Expansion at the origin badness ",I6)'
                If (Kout.GT.0) Write( *,strfmt) ni,ir
                Write(11,strfmt) ni,ir
            End If

            If(irr.EQ.0) nigd=ni

            ! test of expansion at the origin
            Call WriteF(12,ni+4+Ns,CP,CQ,2)
        End Do

        Call ReadF (12,1,P,Q,2)
        P(18)=Ecore
        P(19)=Hcore
        Call WriteF(12,1,P,Q,2)

        strfmt='(/4X,"one-el. and total core energy:",F17.7,4X,F17.7)'
        Write( *,strfmt) Hcore,Ecore
        Write(11,strfmt) Hcore,Ecore

        If(irr.GT.0) Then
            strfmt='(4X,"Only",I4," first orbitals are good at the origin", &
                  /4X,"for other orbitals there were",I6," errors.", &
                  /4X,"See RES file and use Kout=2 for more details.")'
            Write( *,strfmt) nigd,irr
            Write(11,strfmt) nigd,irr
            If (Kout.GT.0) Then
                Write(*,*) '  Push...'
                Read(*,*)
            End If
        End If
        Close(12)

        Return
    End Subroutine Core

    Subroutine Rint0(nsx,nsx2,lsx) 
        ! This subroutine counts the number of radial integrals
        Implicit None
        Integer :: nsx, lsx, i, is, ngint2, nmin, n0, na, nb, nc, nd, la, lb, lc, ld
        Integer :: ja, jb, jc, jd, n1, n2, n3, n4, iis, kmin, kmax, k1, k, ierr, nsx2
        Character(Len=512) :: strfmt

        ! nsx and lsx are used to eliminate integrals:
        nsx=0
        lsx=0
        Do i=Nso+1,Nsp
           is=Nip(i)
           If (is.GT.nsx) nsx=is
           If (Ll(is).GT.lsx) lsx=Ll(is)
        End Do

        Call Nxt_to_max(nsx,nsx2)

        nhint=0
        ngint=0
        ngint2=0

        ! parameter for indexation of integrals:
        If (nsx.GT.IPx+Nso-1) Then
           Write (*,*) ' Rint0: IPx is too small'
           Stop
        End If
        If (nsx.NE.Nso) Then
            Write(*,*) 'Counting of radial integrals'
            nmin=Nso+1
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            ! counting of one-electron integrals
            Do na=nmin,Nsx
                If (Ll(na).GT.lsx) Cycle
                Do nb=na,Nsx
                    If (Kk(na).NE.Kk(nb)) Cycle
                    nhint=nhint+1
                End Do
            End Do
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            ! counting of integrals R(K;NA,NB,NC,ND)
            If (Ne.GT.1) Then
             n0=Nso+1
             Do n1=n0,Nsx
              If (Ll(n1).GT.lsx)    Cycle
              Do n2=n1,Nsx
               If (Ll(n2).GT.lsx)   Cycle
               Do n3=n1,Nsx
                If (Ll(n3).GT.lsx)  Cycle
                Do n4=n3,Nsx
                 If (Ll(n4).GT.lsx) Cycle
                 na=n1
                 nb=n2
                 nc=n3
                 nd=n4
                 Do iis=1,2
                  If (nb.LE.nd) Then
                   la=Ll(na)
                   If (Nn(na).GT.N_l(la+1)) Continue
                   lb=Ll(nb)
                   If (Nn(nb).GT.N_l(lb+1)) Continue
                   lc=Ll(nc)
                   If (Nn(nc).GT.N_l(lc+1)) Continue
                   ld=Ll(nd)
                   If (Nn(nd).GT.N_l(ld+1)) Continue
                   i=la+lb+lc+ld
                   If (i.NE.2*(i/2)) Exit
                   ja=Jj(na)
                   jb=Jj(nb)
                   jc=Jj(nc)
                   jd=Jj(nd)
                   kmin=max(iabs(ja-jc)/2+1,iabs(jb-jd)/2+1)
                   kmax=min((ja+jc)/2+1,(jb+jd)/2+1)
                   Do k1=kmin,kmax
                     k=k1-1
                     i=k+la+lc
                     If (i.NE.2*(i/2).AND.kbrt.EQ.0) Cycle
                     If (n1.LE.nsx2) Then
                       ngint=ngint+1
                     Else
                       ngint2=ngint2+1
                     End If
                   End Do
                   If (n1.EQ.n2.OR.n3.EQ.n4) Exit
                   nc=n4
                   nd=n3
                   Continue
                  End If
                 End Do
                End Do
               End Do
              End Do
             End Do
            End If
        End If

        Write  (*,*) ' Nmax=',nsx,', Lmax=',lsx
        Write (11,*) ' Nmax=',nsx,', Lmax=',lsx

        strfmt='(4X,"NHINT=",I5)  '
        Write( *,strfmt) nhint
        Write(11,strfmt) nhint

        nsx2=nsx
        ngint=ngint+ngint2

        strfmt='(4X,"NGINT=",I9," including",I9," unused")'
        Write( *,strfmt) ngint,ngint2
        Write(11,strfmt) ngint,ngint2
        
        Return
    End Subroutine Rint0

    Subroutine Nxt_to_max(nsx,nsx2) 
        ! next to max orbital in each configuration
        Implicit None             

        Integer :: nsx2, i, ic, kx, ncx, kcx, nek, k, nsx

        nsx2=0
        i=Nso
        Do ic=1,Nc
            kx=i+Nvc(ic)
            ncx=0                 !### ncx= max orbital in conf ic
            kcx=0                 !### index of max orbital
            Do k=i+1,kx
                If (Nip(k).GT.ncx) Then
                    ncx=Nip(k)
                    kcx=k
                    nek=Nq(k)
                End If
            End Do
            If (nek.GE.2) Then    !### next to max orbital
                nsx2=max(nsx2,ncx)
            Else
                Do k=i+1,kx
                    If (Nip(k).GT.nsx2.AND.k.NE.kcx) nsx2=Nip(k)
                End Do
            End If
            i=kx
        End Do

        If (i.NE.Nsp) Then
            Write(*,*) 'Nxt_to_max error: i=',i,' <> ',Nsp
            Read(*,*)
        End If

        If (nsx2.GT.nsx) Then
            Write(*,*) 'Nxt_to_max error: nsx2=',nsx2,' > ',nsx
            Read(*,*)
            nsx2=nsx
        End If

        Return
    End Subroutine Nxt_to_max

    Subroutine Rint(nsx,nsx2,lsx)
        Use breit
        Use readfff
        Use sintg
        Use str_fmt, Only : FormattedTime
        Implicit None

        Integer :: lnx, num_is, nx, i, Idel, idel1, nsx, ih, nmin, na, lsx
        Integer :: nna, lla, jjd, lb, la, nnc, nnb, iis, nb, n4, n3, n2, n1
        Integer :: nsx2, kmin, kmax, k, k1, lc, ld, ja, jb, jc, jd, Iab, nm_br
        Integer :: ln, kac, kbd, n0, nad, jja, nnd, lld
        Real(dp) :: dint1, r_e1, r_e2, fis, tab_br, tad_br, r_br, rabcd, r_c, tad
        Real(dp), Dimension(IP6) :: p, q, a, b, cp, cq, ca, cb
        Integer, Dimension(11) :: ilogr
        Real(dp), Dimension(20) :: Rint_br
        Integer(kind=int64) :: clock_rate, start_time, end_time
        Real :: total_time
        Character(Len=16) :: timeStr
        Character(Len=512) :: strfmt
        Character(Len=1) :: let(9)
        Logical*1  ::  l_br
        data let/'s','p','d','f','g','h','i','k','l'/

        Call system_clock(count_rate=clock_rate)

        ! nsx and lsx are used to eliminate integrals
        allocate(Rint1(nhint),Iint1(nhint),I_is(nhint),R_is(nhint))
        allocate(Rint2(IPbr,ngint),Iint2(ngint),Iint3(ngint))
        allocate(IntOrd(IPx*IPx))
        small=1.d-8
        nhint=0
        ngint=0
        lnx=0
        kt1=kt
        fis=0.d0                   !## rad. int. for Volume shift (VS)
        num_is=0                   !## number of IS integrals (VS or SMS)
        !       parameter for indexation of integrals:
        nx = IPx
        Do i=1,11                  ! this array is used to count
            ilogr(i)=0             !# two-electron integrals of
        End Do                     !## particular abs value
        If (Kout.EQ.0) Idel=50000
        If (Kout.EQ.1) Idel=5000
        If (Kout.GE.2) Idel=1
        idel1=max(idel,10000)
        !     >>>>>>>>>>>>>> ECP of the core >>>>>>>>>>>>>>>>>>
        If (Kecp.NE.0) Then
            Open(unit=10,file='CONF.ECP',status='OLD')
            strfmt='(F7.3)'
            Read (10,strfmt) Alfd
            Read (10,strfmt) (rcut(i),i=1,10)
            strfmt='(5X,"ECP: Alfd =",F7.3,5X,"Rcut:"/5(I2,":",F7.3))'
            Write(11,strfmt) Alfd,(i-1,rcut(i),i=1,10)
            Close(unit=10)
        End If
        !     <<<<<<<<<<<<<<    End of ECP   <<<<<<<<<<<<<<<<<<
        strfmt='(14X,"Radial integrals",2X,/4x,67("="),  &
               /5X,"N",12X,"NA",8X,"NB",8X,"H0(NA,NB)",&
               6X,"Br(Na,Nb)",6X,"VS(NA,NB)",/4x,67("-"))'
        Write(11,strfmt)

        If (nsx.EQ.Nso) Then
            Continue
        Else
            Open(12,file=FNAME,status='OLD',access='DIRECT',recl=2*IP6*IPmr)
            ih=2-Kt
            nmin=Nso+1
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            !   evaluation of one-electron integrals
            Call system_clock(start_time)
            Do na=nmin,Nsx
                If (Ll(na).GT.lsx) Cycle
                Call ReadFF (12,na+4,Pa,Qa,2)      ! p,q
                Call Readf (12,na+Ns+4,cp,cq,2)
                nna=nn(na)
                lla=ll(na)+1
                If (nna.GT.N_l(lla)) Cycle
                jja=jj(na)
                Do nd=na,Nsx
                    If (Kk(na).NE.Kk(nd)) Cycle
                    tab_br=0.d0
                    nnd=nn(nd)
                    lld=ll(nd)+1
                    If (nnd.GT.N_l(lld)) Cycle
                    jjd=jj(nd)
                    Call ReadFF (12,nd+4,Pd,Qd,2)    ! a,b
                    Call readf (12,nd+Ns+4,ca,cb,2)
                    Do i=1,Ii,ih
                        C(i)=cp(i)*Pd(i)+cq(i)*Qd(i)
                        ! >>>>>>>>>>>>>>> EFFECTIVE CORE POTENTIAL >>>>>>>>>>>>>>>>>
                        If(Kecp.NE.0) C(i)=C(i)- &
                             0.5d0*(Pa(i)*Pd(i)+Qa(i)*Qd(i)) &
                             *Alfd/(R(i)**2+Rcut(lla)**2)**2
                    End Do
                    C(ii+4)=Pa(ii+4)+Pa(ii+4)-1.d0
                    Call Sint1(tad)
                    If (jja.EQ.1) Then                    !## Correction of the
                        Call NclInt(cp,cq,Pd,Qd,dint1)    !## integral inside the
                        tad=tad+dint1-Dint                !## nucleus for j=1/2
                    End If
                    If (K_is.EQ.1) fis = dV_nuc(Pa,Qa,Pd,Qd)     !### Volume IS
                    If (Kbrt.GE.1) tad_br=Br_core(na,nd,12)      !### Core Breit
                    nhint=nhint+1
                    nad=nx*(na-nso-1)+(nd-nso)
                    Write(11,'(I6,8X,I3,A1,1X,I1,"/2",2X,I3,A1,1X,I1,"/2",3F15.8)') &
                        nhint,nna,let(lla),jja,nnd,let(lld),jjd,tad,tad_br,fis
                    Rint1(nhint)=tad+tad_br
                    Iint1(nhint)=nad
                    R_is(nhint)=fis
                    I_is(nhint)=nad
                End Do
            End Do
            
            Call system_clock(end_time)
            total_time=Real((end_time-start_time)/clock_rate)
            Call FormattedTime(total_time, timeStr)
            Write(*,'(2X,A)'), 'TIMING >>> Evaluation of one-electron integrals took '// trim(timeStr)
            If (K_is.EQ.1) num_is=nhint
    
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            !   evaluation of integrals R(K;NA,NB,NC,ND)
            Call system_clock(start_time)
            If (Ne.NE.1) Then
                kt1=Kt
                ih=2-Kt
    
                strfmt = '(3x,69("=")/5X,"N",4X,"K",5X,"NA",8X,"NB",8X,"NC", &
                            7X,"ND",3X,"R(K;NA,NB,NC,ND)  R_br"/3x,69("-"))'
                Write( *,strfmt)
                Write(11,strfmt)
    
                n0=Nso+1
                r_br=0.d0
                l_br=Kbrt.NE.0.AND.IPbr.EQ.2
                nm_br=0
                Do n1=n0,Nsx2
                    If (Ll(n1).GT.lsx) Cycle
                    Do n2=n1,Nsx
                        If (Ll(n2).GT.lsx) Cycle
                        Iab= nx*(n1-Nso-1)+(n2-nso)
                        IntOrd(Iab)= ngint+1
                        Do n3=n1,Nsx
                            If (Ll(n3).GT.lsx) Cycle
                            Do n4=n3,Nsx
                                If (Ll(n4).GT.lsx) Cycle
                                na=n1
                                nb=n2
                                nc=n3
                                nd=n4
                                Do iis=1,2
                                    If (nb.LE.nd) Then
                                        nna=Nn(na)
                                        nnb=Nn(nb)
                                        nnc=Nn(nc)
                                        nnd=Nn(nd)
                                        la=Ll(na)
                                        If (nna.GT.N_l(la+1)) Continue
                                        lb=Ll(nb)
                                        If (nnb.GT.N_l(lb+1)) Continue
                                        lc=Ll(nc)
                                        If (nnc.GT.N_l(lc+1)) Continue
                                        ld=Ll(nd)
                                        If (nnd.GT.N_l(ld+1)) Continue
                                        i=la+lb+lc+ld
                                        If (i.NE.2*(i/2)) Exit
                                        ja=Jj(na)
                                        jb=Jj(nb)
                                        jc=Jj(nc)
                                        jd=Jj(nd)
                                        kmin=max(iabs(ja-jc)/2+1,iabs(jb-jd)/2+1)
                                        kmax=min((ja+jc)/2+1,(jb+jd)/2+1)
                                        Call ReadFF (12,na+4,Pa,Qa,2)
                                        Call ReadFF (12,nc+4,Pc,Qc,2)
                                        Call ReadFF (12,nb+4,Pb,Qb,2)
                                        Call ReadFF (12,nd+4,Pd,Qd,2)
                                        If (l_br) Then  ! valence Gaunt included
                                            xja=0.5d0*ja
                                            xla=la
                                            yla=ja-la
                                            xjb=0.5d0*jb
                                            xlb=lb
                                            ylb=jb-lb
                                            xjc=0.5d0*jc
                                            xlc=lc
                                            ylc=jc-lc
                                            xjd=0.5d0*jd
                                            xld=ld
                                            yld=jd-ld
                                            Do k = kmin, kmax
                                                Rint_br(k) = breit_int(k-1,na,Pa,Qa,nb,Pb,Qb,nc,Pc,Qc,nd,Pd,Qd)
                                            enddo
                                        End If
                                        Do k1=kmin,kmax
                                            k=k1-1
                                            i=k+la+lc
                                            If (i.NE.2*(i/2).AND..NOT.l_br) Cycle
                                            If (i.NE.2*(i/2)) Then
                                                r_c=0.d0
                                            Else 
                                                Do i=1,Ii,ih
                                                    C(i)=Pa(i)*Pc(i)+Qa(i)*Qc(i)
                                                End Do
                                                C(ii+4)=Pa(ii+4)+Pc(ii+4)
                                                Call Yk(k)
                                                Do i=1,Ii,ih
                                                    C(i)=(Pb(i)*Pd(i)+Qb(i)*Qd(i))*C(i)/R(i)
                                                End Do
                                                C(ii+4)=Pd(ii+4)+Pb(ii+4)+k
                                                Call Sint(r_c)
                                                ! >>>>>>>>>>>>>>> EFFECTIVE CORE POTENTIAL >>>>>>>>>>>>>>>>>
                                                If (k.EQ.1.AND.Kecp.NE.0) Then
                                                    Do i=1,ii,ih
                                                        C(i)=(Pa(i)*Pc(i)+Qa(i)*Qc(i))/(R(i)**2+Rcut(la+1)*Rcut(lc+1))
                                                    End Do
                                                    C(ii+4)=Pa(ii+4)+Pc(ii+4)
                                                    Call Sint1(r_e1)
                                                    Do i=1,ii,ih
                                                        C(i)=(Pb(i)*Pd(i)+Qb(i)*Qd(i))/(R(i)**2+Rcut(lb+1)*Rcut(ld+1))
                                                    End Do
                                                    C(ii+4)=Pb(ii+4)+Pd(ii+4)
                                                    Call Sint1(r_e2)
                                                    r_c=r_c-r_e1*r_e2*Alfd
                                                End If
                                            End If
                                            rabcd=r_c
                                            If (l_br) r_br=Rint_br(k1)  ! - Breit/Gaunt(+Breit/retardation for kbrt=2)
                                            ngint=ngint+1
                                            kac=nx*nx*K+nx*(na-Nso-1)+(nc-Nso)
                                            kbd=nx*(nb-Nso-1)+(nd-Nso)
                                            strfmt='(1X,I8,1X,I2,2X,I3,A1,1X,I1,"/2",1X,I3, &
                                                 A1,1X,I1,"/2",1X,I3,A1,1X,I1,"/2", &
                                                 1X,I3,A1,1X,I1,"/2",2F13.7)'
                                            If (ngint.EQ.(ngint/idel1)*idel1) &
                                                 Write (*,strfmt) ngint,k,nna,let(la+1),ja, &
                                                                       nnb,let(lb+1),jb, &
                                                                       nnc,let(lc+1),jc, &
                                                                       nnd,let(ld+1),jd,r_c,r_br
                                            If (ngint.EQ.(ngint/Idel)*Idel) &
                                                 Write (11,strfmt) ngint,k,nna,let(la+1),ja, &
                                                                       nnb,let(lb+1),jb, & 
                                                                       nnc,let(lc+1),jc, &
                                                                       nnd,let(ld+1),jd,r_c,r_br
                                            rint2(1,ngint)=rabcd
                                            If (l_br) rint2(IPbr,ngint)=r_br
                                            iint2(ngint)=kac
                                            iint3(ngint)=kbd
                                            ln=1.d0-dlog10(dabs(rabcd)+dabs(r_br)+1.d-9)
                                            ln=min(ln,11)
                                            ln=max(ln,1)
                                            lnx=max(lnx,ln)
                                            ilogr(ln)=ilogr(ln)+1
                                            If (l_br.AND.dabs(r_br).GT.1.d-9) nm_br=nm_br+1
                                        End Do
                                        If (n1.EQ.n2.OR.n3.EQ.n4) Exit
                                        nc=n4
                                        nd=n3
                                    End If
                                End Do
                            End Do
                        End Do
                    End Do
                End Do
                Kt=kt1
            Else
                Continue
            End If
        End If

        Call system_clock(end_time)
        total_time=Real((end_time-start_time)/clock_rate)
        Call FormattedTime(total_time, timeStr)
        Write(*,'(2X,A)'), 'TIMING >>> Evaluation of two-electron integrals took '// trim(timeStr)

        Close(unit=12)

        strfmt = '(2X,"NHINT=",I5," NGINT=", &
               I8, &
               /(4X,"10**(-",I2,") < |R| < 10**(-",I2,") num =",I8))'
        Write( *,strfmt) nhint,ngint,(i,i-1,ilogr(i),i=1,lnx)
        Write(11,strfmt) nhint,ngint,(i,i-1,ilogr(i),i=1,lnx)

        ! >>>>>>>>>>>>>> MS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        If (K_is.GE.2) Then
            Open(13,file='HFD.DAT',status='OLD',access='DIRECT',recl=2*IP6*IPmr)
            Call Rint_MS(nsx,lsx,num_is)
            Close(13)
        End If
        ! <<<<<<<<<<<<<< MS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        Write(*,*)' Formation of the file CONF.INT...'
        Open (unit=13,file='CONF.INT',status='UNKNOWN',form='UNFORMATTED')
        Write (13) Ns,Nso,Nsp,Nsx,Ecore
        Write (13) (Nn(i),Kk(i),Ll(i),Jj(i), i=1,Nsx)
        Write (13) (Nq(i),Nip(i), i=1,Nsp)
        Write (13) nhint,kbrt
        Write (13) (Rint1(i), i=1,nhint)
        Write (13) (Iint1(i), i=1,nhint)
        Write (13) ngint,0,nx*nx
        Write (13) ((Rint2(k,i),k=1,IPbr), i=1,ngint)
        Write (13) (Iint2(i), i=1,ngint)
        Write (13) (Iint3(i), i=1,ngint)
        Write (13) (IntOrd(i), i=1,nx*nx)

        If (K_is.GE.1) Then
            Write(13) K_is,Klow,num_is
            Write(13) (R_is(i),i=1,num_is)
            Write(13) (I_is(i),i=1,num_is)
        End If

        Close(unit=13)

        Write(*,*)' File CONF.INT formed'
        Write(*,*)' nsx = ',nsx,' lsx =',lsx
        Write(*,*)' nhint=',nhint,' ngint=',ngint,' num_is=',num_is

        If (l_br) Then
            If (kbrt.eq.1) Write(*,*)' Valence Magnetic Breit is included:'
            If (kbrt.eq.2) Write(*,*)' Valence Full Breit is included:'
            Write (*,*) nm_br,' integrals exceed 1.d-9'
        End If

       Return

    End Subroutine Rint

    Real(dp) Function dV_nuc(P,Q,A,B)  !### calculates volume shift
        Implicit none
        Integer :: m, im, ik, k
        Real(dp) :: r1, gm, V0, de, s
        Real(dp), Dimension(IP6) :: P, Q, A, B

        r1=R(1)
        gm=2*P(ii+4)
        V0=3*Z*r1**gm

        de=0.d0
        Do m=0,MaxT
            im=ii+m+5
            s=0.d0
            Do k=0,m
                ik=ii+k+5
                s=s+P(ik)*A(im-k)+Q(ik)*B(im-k)
            End Do
            de=de+s/((gm+m+1)*(gm+m+3))
        End Do

        dV_nuc=V0*de

        Return
    End Function dV_nuc

    Subroutine NclInt(p,q,a,b,dn)
        ! radial integration inside the nucleus: int_0,1 f(x) dx
        ! where x = r/R(1)  and  f(x) = (p(x)*q(x) + a(x)*b(x))
        Implicit none
        Integer :: i, i1, j, j1
        Real(dp) :: dn, g12
        Real(dp), dimension(IP6) :: p, q, a, b

        dn=0.d0
        g12=p(ii+4)+a(ii+4)
        Do i=0,MaxT
            i1=Ii+5+i
            Do j=0,MaxT-i
                j1=Ii+5+j
                dn = dn + (p(i1)*a(j1)+q(i1)*b(j1))/(g12+i+j+1)
            End Do
        End Do
        dn=dn*R(1)**(g12+1)

        Return
    End Subroutine NclInt

    Subroutine Rint_MS(nsx,lsx,num_is)
        Use pi_pk
        Use readfff
        Implicit None

        Integer :: nsx, lsx, num_is, nmin, na, nna, lla, jja, nb, nnb, llb, jjb
        Real(dp) :: tab
        character(Len=1) :: let(9), case*5,chms(3)*3
        logical*1 one_e,two_e
        Character(Len=512) :: strfmt
        data let/'s','p','d','f','g','h','i','k','l'/
        data chms/'SMS','NMS','MS '/

        If (nsx.EQ.Nso.OR.K_is.LE.1) Return

        num_is=0
        strfmt='(14X,A3," verteces",2X,/4x,48("="),/5X,"N",12X, &
                    "NA",8X,"NB",5X,"P(NA,NB)",/4x,48("-"))'
        Write( *,strfmt) chms(K_is-1)
        Write(11,strfmt) chms(K_is-1)

        nmin=Nso+1

        ! evaluation of one-electron SMS and NMS verteces
        Do na=nmin,nsx
            If (Ll(na).GT.lsx) Cycle
            nna=nn(na)
            lla=ll(na)
            jja=jj(na)
            Do nb=na,nsx
                If (Ll(nb).GT.lsx) Cycle
                nnb=nn(nb)
                llb=ll(nb)
                jjb=jj(nb)
        
                one_e=lla.EQ.llb.AND.jja.EQ.jjb
                two_e=iabs(lla-llb).EQ.1.AND.iabs(jja-jjb).LE.2
                two_e=two_e.AND.K_is.NE.3             ! No two-e integrals for NMS
                If (.NOT.one_e.AND..NOT.two_e) Cycle
        
                tab=0.d0
                If(one_e) Then
                    case='one_e'
                    If (K_is.NE.3) tab=SMS_core(na,nb,13)  ! SMS included for K_is=2,4
                    If (K_is.NE.2) tab=tab+V_nms(na,nb,13) ! NMS included for K_is=3,4
                End If
                If (two_e) Then
                    case='two_e'
                    Call Readf (13,na+4,Pa,Qa,2) !### Do not move to the outer loop!
                    Call Readf (13,nb+4,Pb,Qb,2)
                    tab=P_eff(jja,lla,jjb,llb,Pa,Qa,Pb,Qb)
                End If

                num_is=num_is+1
                
                strfmt='(I6,2X,A5,1X,I3,A1,1X,I1,"/2",2X,I3,A1,1X,I1,"/2",F17.7)'
                Write( *,strfmt) num_is,case,nna,let(lla+1),jja,nnb,let(llb+1),jjb,tab
                Write(11,strfmt) num_is,case,nna,let(lla+1),jja,nnb,let(llb+1),jjb,tab
                R_is(num_is)=tab
                I_is(num_is)=IPx*(na-Nso-1)+(nb-Nso)
            End Do
        End Do

        Write(*,*)' R_is finished'

        Return
    End Subroutine Rint_MS

End Module basc_aux