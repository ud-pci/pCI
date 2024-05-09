Program bass
    Use basc_variables, Kbasparam => Kbas, Qqparam => Qq, letparam => let, Rint2param => Rint2
    Use readfff
    Implicit None
    
    Integer :: irec, Nsint, Kdg, Ksg, Khf, Kkin, Ierr, Norb, Nsv, ni, iort, i, n,  &
                nlst, nr, nmin, idg, ier_A, n11, n1, j1, n3, j3, n4, j4, n5, j5, l, ii3, Ns3
    Real(dp) :: Sf, Emin, Emax

    Logical :: intrpl
    Integer, Dimension(IPs) :: Jj3, Nn3, Ll3
    Integer, Dimension(IPs) :: kw, numw, nb, jjw ! PW
    Integer, Allocatable, Dimension(:) :: Kbas, Rint2
    Integer, Dimension(IP6) :: Nintrpl
    Real(dp), Dimension(IP6) :: Y, CP, CQ, CA, CB, A, B, R3, V3
    Real(dp), Allocatable, Dimension(:) :: Qnl1, Qq
    Character(Len=1) :: let(9), ch1, ch3, ch4, ch5, FNAME*12, str_is(4)*3, ltr(IPs), FNAME3*12
    Character(Len=256) :: strfmt

    data let /'s','p','d','f','g','h','i','k','l'/
    data str_is /' VS','SMS','NMS',' MS'/

    Y=0_dp
    CP=0_dp
    CQ=0_dp
    CA=0_dp
    CB=0_dp
    A=0_dp
    B=0_dp
    R3=0_dp
    V3=0_dp

    call recunit
    IPmr=1           ! word length in direct access files
    MaxT=9           !### length of expansion at the origin
    kout=1           !### output details
    irec=2*IP6*IPmr  !### record length in DAT files
    small=1.d-8
    FNAME='BAS_A.DAT'
    Kt1=0            !### used as Kt for Breit integrals
    Kbrt=0           !### 0 - Coulomb, 1 - Gaunt, 2 - Full Breit
    K_is=0           !### 0 - no IS, 1 - Volume shift, 2 - SMS
    C_is=0.d0        !### scaling factor for IS
    Nsint=0          !### number of MEs of Sigma
    Khf=0            !### =1 if self-consistency is reached
    kautobas=0       !### =1 if automatic list of orbitals is used

    ! First part of the code (bas_a)
    Ierr=0           !### number of nonfatal errors or warnings
    open(unit=11,file='BASS.RES',status='UNKNOWN')
    call Input
    call Init(Norb)      !### Norb = number of orbitals to be formed

    Nsv=N_orb(n4,ch4,j4) !### last "frozen" orbital
    if (Ns.GT.Nsv) then  !### Orthogonalization of existing orbitals
        open(12,file=FNAME,status='OLD',access='DIRECT',recl=irec)
        do ni=Nsv+1,Ns
            call Ort(ni,Ns+Norb,iort)
            if (iort.GE.10) then
                strfmt = '(4X,"Error: orbital ",I3," is close to linear dependence")'
                write( *,strfmt) ni
                close(12,status='DELETE')
                stop
            end if
        end do
        close(12)
    end if

    open(12,file=FNAME,status='OLD',access='DIRECT',recl=irec)
    if (Norb.GT.0) call Fbas(Norb)  !### forms new orbitals

    open(13,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)
    do ni=1,2*Ns+4
        call ReadF (12,ni,P,Q,2)
        call WriteF(13,ni,P,Q,2)
    end do
    strfmt = '(/4X,"Basis set of ",I3," orbitals is formed.",/4X,"Level of linear dependence of initial orbitals: ",I6)'
    write( *,strfmt) Ns,Ierr
    write(11,strfmt) Ns,Ierr

    ! Second part of the code (bas_b)
    ier_A=Ierr
    Ierr=0           !### number of or warnings

    strfmt = '(1X,70("=")/4X,"Rotation of the basis set",/1X,70("="))'
    write( *,strfmt)
    write(11,strfmt)

    Nsv=N_orb(n4,ch4,j4) !### last "frozen" orbital
    if (Khf.EQ.1) Nsv=max(Nsv,Nso)
    n11=Nsv+1
    if (Nsv.GT.Nso) then      !### dumping factor for i_diag
        Sf=0.75d0
    else
        Sf=1.00d0
    end if

    if (Ns.GE.n11) then
        call V_core                         !### Y(r)=V_core-Z/r

        nmin=N_orb(n3,ch3,j3)               !### first orbital to apply
        if (nmin.EQ.0) nmin=Ns+1            !#### kinetic balance condition

        do ni=n11,Ns
            if (ni.GE.nmin) call Change_Q(ni) !### kinetic balance for Q
            call Ort(ni,Ns,iort)
            if (iort.GE.1) then
                if (ni.GE.nmin) then
                    strfmt = '(4X,"Big change in Q for orbital",I4," iort =",I4)'
                    write( *,strfmt) ni,iort
                    write(11,strfmt) ni,iort
                else
                    strfmt = '(4X,"No orthogonality for",I4," iort =",I4)'
                    write( *,strfmt) ni,iort
                    write(11,strfmt) ni,iort
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
        do ni=n11,Ns
            call Ort(ni,Ns,iort)
            if (iort.GE.1) then
                strfmt = '(4X,I3,": orthogonality is lost. iort=",I4)'
                write( *,strfmt) ni,iort
                write(11,strfmt) ni,iort
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
        strfmt = '(/" Diagonalization key was changed to",I2)'
        write( *,strfmt) idg
        write(11,strfmt) idg
        Ierr=Ierr+100
    end if

    if (n5.GT.0) then
        nlst=N_orb(n5,ch5,j5) !### last orbital to be kept in the set
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
            Ns=nlst
            strfmt = '(/" First ",I3," orbitals saved in ",A12)'
            write( *,strfmt) Ns,FNAME
            write(11,strfmt) Ns,FNAME
        end if
    end if

    strfmt = '(/25x,"Basis set is formed:",/5(i6,i3,a1,i2,"/2"))'
    write( *,strfmt) (i,Nn(i),let(LL(i)+1),JJ(i),i=1,Ns)
    write(11,strfmt) (i,Nn(i),let(LL(i)+1),JJ(i),i=1,Ns)

    if (ier_A+Ierr.GT.0) then
        strfmt = '(/" Constuction of orbitals:",I4," warnings;"/" Rotation of orbitals:   ",I4," warnings.")'
        write( *,strfmt) ier_A,Ierr
        write(11,strfmt) ier_A,Ierr
    end if

    if (Nsint.GT.0.AND.idg.GT.0) then
        strfmt = '(" You may need to recalculate Sigma now.")'
        write( *,strfmt)
        write(11,strfmt)
    end if

    if (dabs(C_is).GT.1.d-5.AND.K_is.GT.0) then
        strfmt = '(" Note, that ",A3," added with C_is = ",F8.4)'
        write( *,strfmt) str_is(K_is),C_is
        write(11,strfmt) str_is(K_is),C_is
    end if

Contains

    Integer Function N_orb(n,ch,j)
        Implicit None

        Integer :: n, j
        Character(Len=1) :: ch
        Character(Len=256) :: strfmt

        N_orb=0
        do i=1,9
            l=i-1
            if(ch.EQ.Let(i)) goto 200
        end do
        strfmt = '(/4X,"N_orb error for input: n=",i3,", ch= ",a1,", j=",i3/ &
              4X,"no l associated with letter ",A1,/4X,"Known letters are:",9A2)'
        write( *,strfmt) n,ch,j,ch,Let
        write(11,strfmt) n,ch,j,ch,Let
        stop

200     do i=1,Ns
            if (n.EQ.Nn(i).AND.l.EQ.Ll(i).AND.j.EQ.Jj(i)) then
                N_orb=i
                return
            end if
        end do
    End Function N_orb

    Subroutine recunit
        ! Determination of the record unit.
        Implicit None
        
        Integer :: lrec, iflag, ipmr, nbytes
        Character(Len=8) :: d1,t1,d2,t2

        t1='abcdefgh'
        d1='        '
        t2='hgfedcba'
        d2='        '
        lrec=0
        iflag=1
200     lrec=lrec+1
        if (lrec.gt.8) then
          write(*,*)  'lrec > 8'
          stop
        end if
        open(unit=13,file='test.tmp',status='unknown',access='direct',recl=lrec)
        write(13,rec=1,err=210) t1
        write(13,rec=2,err=210) t2
        read(13,rec=1,err=210) d1
        read(13,rec=2,err=210) d2
        if (d1.ne.t1) goto 210
        if (d2.ne.t2) goto 210
        iflag=0
210     close(unit=13,status='delete')
        if (iflag.ne.0) goto 200
        nbytes=8/lrec
        ipmr=4/nbytes

        Return

    End Subroutine recunit

    Subroutine Input
        Use conf_init, Only : inpstr
        Implicit None

        Integer :: i1, i2, i3, istr, ic
        Character(Len=1) :: name(16), str1(4)*5, str2(3)*3, str3(3)*9
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        data str1,str2,str3 /'   DF','Brcnr','Breit','Br+Bt',' NO','YES','YES','E = V = 0','E=0, V_c ','E_hf, V_c'/

        strfmt = '(/4X,"Program Bass v1.1")'
        write( *,strfmt)
        write(11,strfmt)

        open(unit=10,file='BASS.INP',status='OLD')
        strfmt = '(1X,16A1)'
        read (10,strfmt) name

        strfmt = '(5X,F5.1)'
        strfmt2= '(5X,I5)'
        strfmt3= '(5X,I2,A1,I2)'
        read (10,strfmt) Z
        read (10,strfmt) Am
        read (10,strfmt2) Nso
        read (10,strfmt2) Nc
        read (10,strfmt2) Ksg
        read (10,strfmt2) Kdg
        read (10,strfmt3) n1,ch1,j1
        read (10,strfmt2) Kkin
        read (10,strfmt3) n3,ch3,j3
        read (10,strfmt3) n4,ch4,j4
        read (10,strfmt3) n5,ch5,j5
  
        i1=(Ksg-1)*(Ksg-2)
        i2=Kdg*(Kdg-1)*(Kdg-2)
        i3=(Kkin+1)*(3-Kkin)
        if (i1.NE.0.OR.i2.NE.0.OR.i3.LE.0) then
            write(*,*)' Wrong keys: Ksg=',Ksg,' Kdg=',Kdg
            write(*,*)' Kkin=',Kkin
            stop
        end if

        istr = 0
        Do While (istr /= 1)
            Call inpstr(istr)
        End Do
        strfmt = '(4X,"Atom:",16A1," Z =",F5.1,"  Am =",F5.1,/4X,"Hamiltonian:",A5," for ",I3," shells. "/4X, &
                    "Diagonalization: ",A3," from",I3,A1,I1,"/2",/4X,"Last frozen orbital",I4,A1,I1,"/2"/4X, &
                    "Kinetic balance: ",A9," for",I4,A1,I1,"/2 and up.")'
        write( *,strfmt) name,Z,Am,str1(Ksg),Nso,str2(Kdg+1),n1,ch1,j1,n4,ch4,j4,str3(Kkin+1),n3,ch3,j3
        write(11,strfmt) name,Z,Am,str1(Ksg),Nso,str2(Kdg+1),n1,ch1,j1,n4,ch4,j4,str3(Kkin+1),n3,ch3,j3

        Allocate(Qnl(1000000))
        if (Nso.eq.0) goto 200
        strfmt = '(6(4X,F7.4))'
        read (10,strfmt) (Qnl(i),i=1,Nso)
200     Nsp=Nso
        Allocate(Kbas(1000*Nsp), Qnl1(1000*Nsp))
        if (kautobas.NE.0) then
            call GenOrbList
        else
            do ic=1,Nc
                Nsp=Nsp+1
210             strfmt = '(4X,F7.4,2X,I1,1X,F7.4)'
                read (10,strfmt) Qnl(Nsp),Kbas(Nsp),Qnl1(Nsp)
                if (dabs(Qnl(Nsp)).LT.1.d-2) goto 210
            end do
        end if
        close(unit=10)
        strfmt = '(4X,"Basis set consists of ",I3," orbitals")'
        write( *,strfmt) Nsp
        write(11,strfmt) Nsp

        Return
    End Subroutine Input

    Subroutine GenOrbList
        Implicit None

        Integer :: lmax, ngen, nnn, nnx, id, k, nl, nk, lk, n1, l, n
        Real(dp) :: x, y
        Integer, Dimension(0:8) :: nmin, nmax
        Character(Len=1) :: str1(5)*1, str2*20
        Character(Len=256) :: strfmt

        strfmt = '(/4X,"Generating list of orbitals")'
        write( *,strfmt)
        write(11,strfmt)

100     read (10,'(5a1,a)') (str1(i),i=1,5),str2
        if (str1(1).NE.'l'.AND.str1(1).NE.'L') goto 100
        if (str1(2).NE.'m'.AND.str1(2).NE.'M') goto 100
        if (str1(3).NE.'a'.AND.str1(3).NE.'A') goto 100
        if (str1(4).NE.'x'.AND.str1(4).NE.'X') goto 100
  
        read (str2,*) lmax
        do l=0,lmax
            read (10,'(5a1,a)') (str1(i),i=1,5),str2
            read (str2,*) nmax(l)
            strfmt = '(4x,a1,"-wave; n_max = ",i2)'
            write ( *,strfmt) let(l+1),nmax(l)
            write (11,strfmt) let(l+1),nmax(l)
        end do
  
        ngen=0
        nnn=100
        nnx=0
        do l=0,lmax
            nl=0
            id=1
            if (l.GT.0) id=2
            do k=1,Nsp
                n1=100*dabs(Qnl(k))+small
                nk=10*dabs(Qnl(k))+small
                lk=n1-10*nk
                if (lk.EQ.l) nl=max(nl,nk)
            end do
            nmin(l)=max(l+1,nl+1)
            nnn=min(nmin(l),nnn)
            nnx=max(nmax(l),nnx)
            ngen=ngen+(nmax(l)-nmin(l)+1)*id
        end do
        
        strfmt = '(/4x,"Need ",i3," orbitals for shells ",i2," - ",i2)'
        write(*,strfmt) ngen,nnn,nnx
  
        Nc=0
        strfmt = '(i6,i5,a1,i2,"/2:",f9.4,i3,f9.4)'
        do n=nnn,nnx
            write(*,*)
            do l=0,lmax
                if(n.GE.nmin(l).AND.n.LE.nmax(l)) then
                    Nc=Nc+1
                    Nsp=Nsp+1
                    x=(1000*n+100*l+1)*1.d-4
                    y=(1000*(n-1)+100*l+1)*1.d-4
                    if (n.EQ.l+1) y=(1000*(n-1)+100*(l-1)+1)*1.d-4
                    if (l.GT.0) then
                        Qnl(Nsp)=-x
                        Qnl1(Nsp)=-y
                        Kbas(Nsp)=0
                        write (*,strfmt) Nsp,n,let(l+1),2*l-1,Qnl(Nsp),Kbas(Nsp),Qnl1(Nsp)
                        Nc=Nc+1
                        Nsp=Nsp+1
                    end if
                    Qnl(Nsp)=x
                    Qnl1(Nsp)=y
                    Kbas(Nsp)=0
                  
                    write (*,strfmt) Nsp,n,let(l+1),2*l+1,Qnl(Nsp),Kbas(Nsp),Qnl1(Nsp)
                end if
            end do
        end do
  
        if (Nc.NE.ngen) then
            write(*,*) ' Number of orbitals:',Nc,', expected:',ngen
            stop
        end if

        return
    End Subroutine GenOrbList

    Subroutine Init(Norb)
        Use readfff
        Implicit None
        Integer, Intent(Out) :: Norb

        Integer :: i, i0, ic, imax, j, n, n0, ni, nmin, nsmax, nsmax1, If
        Integer :: nsb, kkj, llj, jjj, nj, ng, nnj, err_stat
        Real(dp) :: c1, c2, r1, z1, d, H0
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
        Qw=0_dp

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
        If (dABS(Z-z1) > 1.d-6) Then
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
        Allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp),Qq(IPs))

        strfmt = '(/4X,"Z   = ",F6.2,5X,"Kt  =",I3,  7X,"II =",I4,    &
                /4X,"H   =",F7.4, 5X,"R2  =",F6.2,4X,"R1 =",E11.4, &
                /4X,"Rnuc=",E11.4,1X,"Al  =",F7.4,3X,"Bt =",F5.2,  &
                /4X,"Nsp =",I7,   7X,"Ns  =",I3,  7X,"Nso=",I3,4X,"Nc =",I7)'
        Write( *,strfmt) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc
        Write(11,strfmt) Z,Kt,II,H,R2,R1,Rnuc,Al,Bt,Nsp,Ns,Nso,Nc

        longbasis=dabs(PQ(20)-0.98765d0) < 1.d-6
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
            If (nj <= Nso) Then
                Do ni=1,Nso
                    If (Nn(ni) == Nn(nj) .and. Ll(ni) == Ll(nj)) Qw(ni) = Qw(ni) + Nq(nj)   ! number of e on NR shell
                End Do
            End If
            Do i=1,nsb
                ni=i
                If (nnj == Nn(ni) .and. Kk(ni) == kkj) goto 210
            End Do
            Norb=Norb+1
            ni=Ns+Norb
            nsb=Ns+Norb 
            Nn(ni)=nnj
            Kk(ni)=kkj
            Ll(ni)=llj
            Jj(ni)=jjj
            If (Kbas(ni) == 3) Kfile=1
            If (Kbas(ni) == 4) kfile=1
 210        Nip(nj)=ni
        End Do
        Write(*,*) Norb, ' new orbitals to be formed'
        nsmax=(4*IP6-20)/6
        nsmax1=(4*IP6-20)/3
        If (nsb > nsmax .and. .not.longbasis) Then
            Write(*,*) ' Ns =',nsb,' > ',nsmax
            Write(*,*) ' switch to long basis variant '
        End If
        If (nsb > nsmax1) Then
            Write(*,*) ' For IP6 =',IP6
            Write(*,*) ' maximum length of basis set is ',nsmax1
            Stop
        End If
        If (Nso /= 0) Then
            Do ni=1,Nso
                Qq(ni) = Qw(ni)*(Jj(ni)+1.d0)/(4*Ll(ni)+2)    !# N_e on a shell
            End Do
        End If
        n0=0
        ne=0
        nmin=Nso+1
        If (nmin <= Nsp) Then
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
            If (n < ne) Cycle
            ic=ic+1
            If (n > ne) Then
                strfmt='(/2X,"wrong number of electrons"/2X,"for configuration ICONF =",I4/)'
                Write( *,strfmt) ic
                Write(11,strfmt) ic
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

    Subroutine Core
        Use wigner
        Use readfff
        Use test_ori

        Implicit None

        Integer :: ih, i, nigd, irr, ni, na, la, ja, ir, ir1, ir2, ir3, ir4, igm, m, m1, im
        Integer :: ik, nj, k, nb, lb, jb, k1, ip, kmin, kmax
        Real(dp) :: hh, Hcore, err1, err2, qa, t, r1, gm, dgm, cpp, cqq, qb, s, d, xk, xja, xjb, ds, yy
        Real(dp) :: hgc, dsb, Ebcore, H0
        Real(dp), Dimension(IP6) :: A, B, CP, CQ, Y
        Real(dp), Dimension(20) :: dd, ss
        Character(Len=512) :: strfmt

        open(12,file=FNAME,status='OLD',access='DIRECT',recl=2*IP6*IPmr)
        Ecore=0.d0
        Hcore=0.d0
        ih=2-kt
        h0=h
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
            call Test_Origin('Core:    P',P,R,MaxT,ii,1.d-6,ir1)
            call Test_Origin('Core:    Q',Q,R,MaxT,ii,1.d-6,ir2)
            call Test_Origin('Core:    A',A,R,MaxT,ii,err1,ir3)
            call Test_Origin('Core:    B',B,R,MaxT,ii,err1,ir4)
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
                      
                        xja=ja/2.d0
                        xjb=jb/2.d0
                        xk=k
                        d=qb*(FJ3(xk,xja,xjb,0.d0,-0.5d0,0.5d0))**2
                        do i=1,ii,ih
                            yy=d*C(i)/(V(i)*hh*R(i))
                            CP(i)=CP(i)-yy*A(i)
                            CQ(i)=CQ(i)-yy*B(i)
                        end do
                      
                        if (nj.LE.ni.AND.ni.LE.Nso) then
                            t=qa
                            if (ni.EQ.nj.AND.k.NE.0) t=0.5d0*(qa-1.d0)*(ja+1.d0)/ja
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

            if (ir.GT.0) then
                irr=irr+1
                strfmt = '(" Orbital",I4,": Expansion at the origin badness ",I6)'
                if (Kout.GT.0) write( *,strfmt) ni,ir
                write(11,strfmt) ni,ir
            end if
            if(irr.EQ.0) nigd=ni
            ! test of expansion at the origin
 240        call WriteF(12,ni+4+Ns,CP,CQ,2)
        end do

        call ReadF (12,1,P,Q,2)
        P(18)=Ecore
        P(19)=Hcore
        call WriteF(12,1,P,Q,2)

        strfmt = '(/4X,"one-el. and total core energy:",F17.7,4X,F17.7)'
        write( *,strfmt) Hcore,Ecore
        write(11,strfmt) Hcore,Ecore

        if(irr.GT.0) then
            strfmt = '(4X,"Only",I4," first orbitals are good at the origin", &
                        /4X,"for other orbitals there were",I6," errors.", &
                        /4X,"See RES file and use Kout=2 for more details.")'
            write( *,strfmt) nigd,irr
            write(11,strfmt) nigd,irr

            if (Kout.GT.0) then
                write(*,*) '  Push...'
                read(*,*)
            end if
        end if
        close(12)

        Return
    End Subroutine Core

    Subroutine Yk(k)            
        !# This routine goes with integration routine Sint, which integrates over ro,
        !  while Sint1 integrates over r.
        Implicit None  

        Integer :: ih, i0, j, k, i, k1, imax, i1, im1
        Real(dp) :: hh, ww, rr, q1, t, p1, p2, gam, cr1, cr2, t1, t2, p
        Real(dp), Dimension(IP6) :: w

        ! evaluation of functions z,y with simpson formula
        ih=2-kt
        i0=ih+1
        h0=h
        hh=h0*ih/3.d0

        ! evaluation of w(i)=r(i)**k
        j=k/2
        k1=k-j*2
        Do i=1,ii,ih
            ww=1.d0
            rr=r(i)
            If (k1.EQ.1) ww=rr
            If (j.EQ.0) Then
                w(i)=ww
                Cycle
            End If
            q1=rr*rr
            ww=ww*q1
            If (j.EQ.1) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            If (j.EQ.2) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            If (j.EQ.3) Then
                w(i)=ww
                Cycle
            End If
            ww=ww*q1
            w(i)=ww
        End Do

        ! evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        Do i=1,i1,ih
           j=k1-i
           If (dabs(c(j)).GE.t) Exit
           imax=imax-ih
        End Do

        If (imax.LT.i0) imax=i0

        ! evaluation of function z(k;r) by direct integration
        p1=c( 1)*v( 1)*w( 1)
        p2=c(i0)*v(i0)*w(i0)
        gam=c(ii+4)
        cr1=0.d0
        cr2=0.d0
        If (c(1).NE.0.d0) cr1=(c(i0)/c(1)*(r(1)/r(i0))**gam-1.d0)/((r(i0)-r(1))*(gam+k+2))
        If (c(i0).NE.0.d0) cr2=(c(1)/c(i0)*(r(i0)/r(1))**gam-1.d0)/((r(1)-r(i0))*(gam+k+2))
        c(1)=c(1)*r(1)/(gam+k+1)*(1.d0-cr1*r(1))
        c(i0)=c(i0)*r(i0)/(gam+k+1)*(1.d0-cr2*r(i0))
        t1=c(i0)*w(i0)
        t2=c(1)*w(1)
        i1=i0+ih
        Do i=i1,imax,ih
            q1=c(i)*v(i)
            If (k.NE.0) q1=q1*w(i)
            t=hh*(q1+4.d0*p2+p1)
            t=t2+t
            t2=t1
            t1=t
            p1=p2
            p2=q1
            c(i)=t/w(i)
        End Do
        t=c(imax)*w(imax)
        im1=imax+ih
        Do i=im1,ii,ih
            p=t
            If (k.NE.0) p=p/w(i)
            c(i)=p
        End Do

        ! evaluation of w(i)=r(i)**(k+1)
        Do i=1,ii,ih
            t=r(i)
            If (k.NE.0) t=t*w(i)
            w(i)=t
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! evaluation of function y(k;r) by solving dif. equation
        ! for f=y(k;r)/r**(k+1)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        k1=k+k+1
        i1=imax+ih
        t1=c(i1)/w(i1)
        p1=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=i1+ih
        t2=c(i1)/w(i1)
        p2=-c(i1)*k1*v(i1)/(r(i1)*w(i1))
        i1=imax+1
        Do j=1,imax,ih
            i=i1-j
            q1=-c(i)*k1*v(i)/(r(i)*w(i))
            t=hh*(q1+4.d0*p1+p2)
            t=t2-t
            t2=t1
            t1=t
            p2=p1
            p1=q1
            c(i)=t*w(i)*v(i)*hh
        End Do

        Do i=im1,ii,ih
            c(i)=c(i)*v(i)*hh
        End Do

        If (c(ii+4).LE.k+0.d0) Then
            c(ii+4)=c(ii+4)+1
        Else
            c(ii+4)=k+1
        End If

        Return
    End Subroutine Yk

    Subroutine Sint(ds)
        Implicit None
        Integer :: ih, i0, i1, i2, i3, i
        Real(dp) :: h1, v0, r0, gam, r1, rr2, r3, t0, p0, f0, g, f1, f2, f3, f21, f32, f321, c0, c1, c2, dt, f, s, t, ds, q, hh, p1, p2, t1, t2

        gam=c(ii+4)
        ih=2-kt
        i0=ih+1
        h0=h
        hh=h0*ih/3.d0

        i1=1
        i2=i1+ih
        i3=i2+ih
        r1=r(i1)
        r2=r(i2)
        r3=r(i3)
        t1=0.d0
        t2=0.d0
        if (gam.GT.5.0d0) goto 200
        f1=c(i1)/(r1**gam*hh*v(i1))
        f2=c(i2)/(r2**gam*hh*v(i2))
        f3=c(i3)/(r3**gam*hh*v(i3))
        f21=(f2-f1)/(r2-r1)
        f32=(f3-f2)/(r3-r2)
        f321=(f32-f21)/(r3-r1)
        c2=f321
        c1=f21-c2*(r1+r2)
        c0=f1-c1*r1-c2*r1**2
        g=gam+1.d0
        t2=r1**g*(c0/g+c1/(g+1.d0)*r1+c2/(g+2.d0)*r1**2)
        t1=r2**g*(c0/g+c1/(g+1.d0)*r2+c2/(g+2.d0)*r2**2)

200     p1=c(i0)
        p2=c( 1)
        i1=i0+ih
        do i=i1,ii,ih
            q=c(i)
            t=q+4.d0*p1+p2
            t=t2+t
            t2=t1
            t1=t
            p2=p1
            p1=q
        End Do
        ds=t

        Return
    End Subroutine Sint

    Subroutine Ort(ni,ns1,ir)  !### ir - linear dependence parameter
        Use diff, Only: Dif
        Use test_ori
        Implicit None

        Integer :: ni, ns1, ir, ih, ir1, ir2, ir3, ir4, nj, iphase, j, i
        Real(dp) :: small, ds, t, sum, d, pp, rr, g, rg, qq, erp, erq, err, ers, snorm

        Real(dp), Dimension(IPs) :: S
        Real(dp), Dimension(IP6) :: A, B

        Character(Len=256) :: strfmt

        small=1.d-6
        ir=0
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
            strfmt = '(4X,"<",I3,"|",I3,"> =",F15.7)'
            write ( *,strfmt) ni,ni,ds
            write (11,strfmt) ni,ni,ds
        end if
        call ReadF(12,ni+ns1+4,A,B,2)
        call Test_Origin('Ort0:dP/dr',A,R,MaxT,ii,2.d-3,ir1)
        call Test_Origin('Ort0:dQ/dr',B,R,MaxT,ii,2.d-2,ir2)
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
        call Test_Origin('Ort1:    P',P,R,MaxT,ii,1.d-6,ir1)
        call Test_Origin('Ort1:    Q',Q,R,MaxT,ii,1.d-5,ir2)
        call Test_Origin('Ort1:dP/dr',A,R,MaxT,ii,2.d-3,ir3)
        call Test_Origin('Ort1:dQ/dr',B,R,MaxT,ii,2.d-2,ir4)
        ir=ir+ir1+ir2+ir3+ir4
        if (ni.EQ.1) return
  
        sum=0.d0
        Do nj=1,ni-1                      !### - orthogonalization
            if (Kk(ni).NE.Kk(nj)) Cycle
            call ReadF (12,nj+4,A,B,2)
            do i=1,ii,ih
                C(i)=P(i)*A(i)+Q(i)*B(i)
            end do
            C(ii+4)=P(ii+4)+A(ii+4)
            call Sint1(ds)
            S(nj)=ds
            sum=sum+ds**2
            if (ds**2.GT.0.9999d0) then
                strfmt = '(4X,"Fatal: <",I3,"|",I3,"> =",F9.6)'
                write(*,strfmt) ni,nj,ds
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
                strfmt = '(4x,"<ni=",I3,"|nj=",I3,"> =",F9.6," Origin: P1/pp =",F9.6," Q1/qq =",F9.6)'
                write(11,strfmt) ni,nj,ds,pp,qq
                if (ers.GT.5.d0) write(*,strfmt) ni,nj,ds,pp,qq
                if(ers.GT.20.d0) then
                    strfmt = '(4X,"Severe: linear dependence ",I4," for orbital ",I3," err=",F8.1)'
                    write (*,strfmt) ir,ni,ers
                end if
  
                if(erp+erq.GT.1) then
                    strfmt = '(4X,"Taylor expansion changed for ",I3)'
                    write ( *,strfmt) ni
                    write (11,strfmt) ni
                    call Test_Origin('Ort2:    P',P,R,MaxT,ii,1.d-6,ir1)
                    call Test_Origin('Ort2:    Q',Q,R,MaxT,ii,1.d-5,ir2)
                    ir=ir+ir1+ir2
                end if
            end if
        End Do
  
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
  
        ! Derivatives of P & Q:
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
        call Test_Origin('Ort2:dP/dr',A,R,MaxT,ii,2.d-3,ir1)
        call Test_Origin('Ort2:dQ/dr',B,R,MaxT,ii,2.d-2,ir2)
  
        call WriteF(12,ni+4+ns1,A,B,2)
        ir=ir+ir1+ir2
        if (ir.GT.0) Then
            strfmt = '(4X,"Ort: linear dependence ",I4," for orbital ",I3)'
            write(11,strfmt) ir,ni
        End If

        Ierr=Ierr+ir

        Return
    End Subroutine Ort

    Subroutine Sint1(DS)        
        ! Simpson integration over r (with weight function HH*V(I))
        Implicit None

        Integer :: I, IH, I0, I1, I2, I3
        Real(dp) :: R1, R2, R3, T1, T2, F1, F2, F3, F21, F32, F321, &
                    C0, C1, C2, G, P1, P2, T, Q, HH, Gam
        Real(dp), intent(out) :: DS

        Gam=C(ii+4)
        IH=2-KT
        HH=H*IH/3.d0
        I0=IH+1
        I1=1
        I2=I1+IH
        I3=I2+IH
        R1=R(I1)
        R2=R(I2)
        R3=R(I3)
        T1=0.d0
        T2=0.d0

        If (GAM <= 5.0d0) Then
            F1=C(I1)/R1**GAM
            F2=C(I2)/R2**GAM
            F3=C(I3)/R3**GAM
            F21=(F2-F1)/(R2-R1)
            F32=(F3-F2)/(R3-R2)
            F321=(F32-F21)/(R3-R1)
            C2=F321
            C1=F21-C2*(R1+R2)
            C0=F1-C1*R1-C2*R1**2
            G=GAM+1.d0
            T2=R1**G*(C0/G+C1/(G+1.d0)*R1+C2/(G+2.d0)*R1**2)
            T1=R2**G*(C0/G+C1/(G+1.d0)*R2+C2/(G+2.d0)*R2**2)
        End If
        P1=C(I0)*V(I0)
        P2=C( 1)*V( 1)
        I1=I0+IH
        Do I=I1,II,IH
           Q=C(I)*V(I)
           T=HH*(Q+4.d0*P1+P2)
           T=T2+T
           T2=T1
           T1=T
           P2=P1
           P1=Q
        End Do
        
        ds=t
        Return
    End Subroutine Sint1

    Subroutine V_core      !### calculates V_core-Z/R
        Implicit None

        Integer :: i

        do i=1,IP6
            C(i)=0.d0                !C is used in Y0
        end do

        call Y0(Y)                  !Y=V_core

        do i=1,ii
          Y(i)=Y(i)-Z/R(i)
        end do
        Y(ii+5)=Y(ii+5)-1.5d0*Z/R(1)
        Y(ii+7)=Y(ii+7)+0.5d0*Z/R(1)

        Return
    End Subroutine V_core
    
    Subroutine Change_Q(ni)
        Use test_ori
        Use diff, Only : Dif, Origin
        ! lower component of the orbital 'ni' is constructed from
        ! Dirac equation and compared with that from the file.
        ! For Kkin>0 Dirac equation is solved for V=-Z(r)/r, otherwise V=0
        ! For Kkin>1 E=-P(Ii+1), otherwise E=0
        Implicit None

        Integer :: ni, kn, n, j, l1, ir1, ir2
        Real(dp) :: small, gj, ep, d, ee, q1, dq, s1, s

        small=1.d-6
        kn=Kk(ni)
        n =Nn(ni)
        j =Jj(ni)
        l1=Ll(ni)+1
        call ReadF(12,ni+4,P,Q,2)
        
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
    
    
        call ReadF(12,ni+ns+4,CP,CP,1)
        d=2*Cl+ep/Cl
        do i=1,Ii
            ee=d
            if (Kkin.GE.1) ee=d-Y(i)/Cl
            q1=-1.d0/ee * (CP(i) + kn/R(i)*P(i))
            dq=q1-Q(i)
            C(i)=dq**2
            Q(i)=q1
        end do
        C(Ii+4)=2*P(Ii+4)
        call Sint1(s1)
        s1=dsqrt(s1/s)
    
        if (s1.GT.small) then      !### Q1 substitutes Q
            call Origin(Q,gj,-kn)    !### Expansion at the origin
            call Dif(Q,CQ,gj)        !### derivative of the lower component
            call Test_Origin('Change_Q:Q',Q,R,MaxT,ii,1.d-5,ir1)
            call Test_Origin('Change_Q:B',CQ,R,MaxT,ii,2.d-2,ir2)
            call WriteF(12,ni+4,P,Q,2)
            call WriteF(12,ni+ns+4,CP,CQ,2)
            Ierr=Ierr+dlog(s1/small)+ir1+ir2
            
            strfmt = '(4X,I2,A1,I2,"/2 (E=",F12.6,"): ||Q|| =",F10.6," ||dQ||/||Q||=",F10.6)'
            write( *,strfmt) n,let(l1),j,ep,dsqrt(s),s1
            write(11,strfmt) n,let(l1),j,ep,dsqrt(s),s1
        end if
        Return
    End Subroutine Change_Q

    Subroutine Y0(Y)      
        Implicit None

        Integer :: ih, i0, igm, km, ki, ki1, imax, i1, k1, j, im1, im
        Real(dp) :: hh, r1, qa, gm1, gm, rg, s, t, ri, rg1, rgi, ri1, t1, t2, fr, p1, p2, q1
        Real(dp), Dimension(IP6) :: Y

        ih=2-kt
        i0=ih+1
        hh=h*ih/3.d0

        ! Core density:
        r1=r(1)
        do ni=1,Nso
            qa=Qq(ni)
            call ReadF (12,ni+4,P,Q,2)
            if (ni.EQ.1) gm1=2*P(ii+4)
            gm=2*P(ii+4)
            igm=gm-gm1+0.1d0
            rg=r1**igm
            do i=1,ii,ih
                C(i)=C(i)+qa*(P(i)**2+Q(i)**2)
            end do
            do m=0,MaxT-igm            !### expansion at the origin
                km=ii+5+m+igm            !#### of the electron density
                s=0.d0
                do i=0,m
                    ki=ii+5+i
                    ki1=ii+5+m-i
                    s=s+(P(ki)*P(ki1)+Q(ki)*Q(ki1))
                end do
                C(km)=C(km)+s*qa*rg
            end do
        end do
        C(ii+4)=gm1

        ! evaluation of imax
        t=1.d-11
        imax=ii-2*ih
        i1=imax
        k1=i1+1
        do i=1,i1,ih
            j=k1-i
            if (dabs(c(j)).GE.t) goto 200
            imax=imax-ih
        end do
200     if (imax.LT.i0) imax=i0

        ! evaluation of z(0;r)=int_o^r(ro*dr) by direct integration
        r1=R(1)
        ri=R(i0)
        gm=C(ii+4)
        rg1=r1**(gm+1)
        rgi=ri**(gm+1)
        ri1=ri/r1
        t2=0.d0
        t1=0.d0
        fr=1.d0
        do m=0,maxt              !### t1=int_0^r1 (ro*dr)
            im=ii+5+m              !### t2=int_0^r(i0) (ro*dr)
            C(im)=C(im)/(gm+m+1)
            t2=t2+C(im)*rg1
            t1=t1+C(im)*fr*rgi
            fr=fr*ri1
        end do
  
        p1=c( 1)*v( 1)
        p2=c(i0)*v(i0)
        c(1)=t2
        c(i0)=t1
        i1=i0+ih
        do i=i1,imax,ih
            q1=c(i)*v(i)
            t=hh*(q1+4.d0*p2+p1)
            t=t2+t
            t2=t1
            t1=t
            p1=p2
            p2=q1
            c(i)=t
        end do
        t=c(imax)
        im1=imax+ih
        do i=im1,ii,ih
            c(i)=t
        end do
        C(ii+4)=gm1+1

        ! evaluation of function y(0;r) by solving dif. equation for f=y(0;r)/r
        i1=imax+ih
        t1=c(i1)/r(i1)
        p1=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=i1+ih
        t2=c(i1)/r(i1)
        p2=-c(i1)*v(i1)/(r(i1)*r(i1))
        i1=imax+1
        do j=1,imax,ih
            i=i1-j
            q1=-c(i)*v(i)/(r(i)*r(i))
            t=hh*(q1+4.d0*p1+p2)
            t=t2-t
            t2=t1
            t1=t
            p2=p1
            p1=q1
            y(i)=t
        end do
        do i=im1,ii,ih
            y(i)=c(i)/r(i)
        end do
  
        y(ii+4)=0.d0                       !### Y0 -> const as r->0
        rg=r1**gm1
        s=0.d0
        do m=0,MaxT
            s=s+rg*C(ii+5+m)/(gm1+m)
        end do
        Y(ii+5)=Y(1)+s
        if (dabs(gm1-2.d0).LT.1.d-5) then  !### the case of finite nucleus
            Y(ii+6)=0.d0
            do m=2,MaxT
                Y(ii+5+m)=-rg*C(ii+3+m)/m
            end do
        else                               !### the case of point nucleus
            do m=1,MaxT
                Y(ii+5+m)=0.d0
            end do
        end if
  
        Return
    End Subroutine Y0

    Subroutine Rint
        Use breit, Only : Br_Core
        Use pi_pk, Only : SMS_Core, V_nms
        Implicit None

        Integer :: ih, nx, na, nb, nab, i
        Real(dp) :: fis, sg, br, ds, sms
        Real(dp), Dimension(IP6) :: CPa, CQa

        character*12 chis*6

        Allocate(Rint1(1000000), Iint1(1000000))

        chis='      '
        if (dabs(C_is).GT.1.d-5.AND.K_is.EQ.1) chis=' VS'
        if (dabs(C_is).GT.1.d-5.AND.K_is.EQ.2) chis='SMS'
        if (dabs(C_is).GT.1.d-5.AND.K_is.EQ.3) chis='NMS'
        if (dabs(C_is).GT.1.d-5.AND.K_is.EQ.4) chis=' MS'
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

        ! evaluation of one-electron integrals
        if (Kbrt.GE.1) write(*,*) ' Calculating Breit integrals...'
        strfmt = '(4X,"   na   nb",8X,"DF",13X,"Sg",13X,"Br",12X,A6)'
        write(11,strfmt) chis
        do na=nmin,Ns
            call ReadF (12,na+4,Pa,Qa,2)
            call ReadF (12,na+Ns+4,CPa,CQa,2)
            if (Kbrt.EQ.2) call ReadF (13,na+Ns+4,P1a,Q1a,2)
            do nb=na,Ns
                if (Kk(na).NE.Kk(nb)) Cycle
                call ReadF (12,nb+4,Pb,Qb,2)
                if (Kbrt.EQ.2) call ReadF (13,nb+Ns+4,P1b,Q1b,2)
                do i=1,Ii,ih
                    C(i)=CPa(i)*Pb(i)+CQa(i)*Qb(i)
                end do
                C(ii+4)=Pb(ii+4)+Pb(ii+4)-1.d0
                call Sint1(ds)
                ierr = 0 ! cAB
                if (Ksg.EQ.2)  sg = Sigma(na,nb,ierr)  !### Brueckner
                if (Kbrt.GE.1) br = Br_Core(na,nb,c,r,v,ii,kt)  !### Breit
                sms=0.d0                               !### SMS, NMS, or both
                if (dabs(C_is).GT.1.d-5) then
                    if (K_is.EQ.1) fis = C_is*dV_nuc(Pa,Qa,Pb,Qb)    !### Volume IS
                    if (K_is.EQ.2.OR.K_is.EQ.4) sms = C_is*SMS_core(na,nb,13)     !### SMS
                    if (K_is.GE.3) sms = sms+C_is*V_nms(na,nb,13)    !### NMS
                end if
                Nhint=Nhint+1
                nab=nx*(na-Nsv-1)+(nb-Nsv)
                rint1(nhint)=ds+sg+br+sms+fis
                iint1(nhint)=nab
                strfmt = '(4X,2I5,4F15.8)'
                write(*,strfmt) na,nb,ds,sg,br,sms+fis
                write(11,strfmt) na,nb,ds,sg,br,sms+fis
            end do
        end do

        close(12)
        close(13)
          
1000    strfmt = '(4X,"Nhint=",I5," Number of zeroes for Sigma:",I5)'
        write( *,strfmt) Nhint,ierr
        write(11,strfmt) Nhint,ierr

        Return
    End Subroutine Rint

    Real(dp) Function Sigma(na,nb,ierr)
        Implicit None

        Integer :: nx, na, nb, ierr, ir, nab0, nab, ia, ib, khot, kval, nmax, IPgS, lmax
        Real(dp) :: s
        ! parameter for indexation of integrals:
        nx = ns - Nsv + 1
        Sigma=0.d0

        IPgS = 1000000
        ! first time:
        if (Nsint.EQ.0) then
            ierr=0
            open(unit=15,file='SGC.CON',status='OLD')
            strfmt = '(7X,I3,5X,I1,5X,I2,7X,I2,7X,I2)'
            read(15,strfmt,err=10) nmax,lmax,kt,kval,khot
            write(*,*) ' nmax=',nmax,' lmax=',lmax,' khot=',khot
            do ir=1,IPgS
                read(15,*,err=10) n,ia,ib,s
                strfmt = '(1X,I4,5X,I3,5X,I3,5X,E12.5)'
                write(*,strfmt) n,ia,ib,s
                nab=nx*(ia-Nsv-1)+(ib-Nsv)
                iint2(n)=nab
                rint2(n)=s
                Nsint=n
            end do
            strfmt = '(5X,I5," matrix elements of Sigma read")'
10          write( *,strfmt) Nsint
            write(11,strfmt) Nsint
            close (15)
            if (Nsint.EQ.0) Stop
        end if

        ! search for the radial integral:
        nab0=nx*(na-Nsv-1)+(nb-Nsv)
        do i=1,Nsint
            nab=iint2(i)
            if (nab.EQ.nab0) then
                Sigma=rint2(i)
                return
            end if
        end do
        ierr=ierr+1

        Return
    End Function Sigma

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

    Subroutine Fbas(Norb)     !### formation of Norb new orbitals
        Use bspl, Only : FormBspl, Grid_Bspl
        Use diff, Only : Origin
        Implicit None

        Integer :: Norb, iw, kcg, ksin, nj, nnj, llj, kkj, i1, n1, l1, kbk, nbk, &
                    nnk, llk, jjk, kkk, nk, m1, mx, if, longbasis, itail, lort, nim, &
                    iwx, jjj, j1
        Real(dp) :: gj, qn
        Integer, Dimension(4*IPs) :: Iqn
        Real(dp), Dimension(IPs) :: Qq1
        Real(dp), Dimension(4*IP6) :: PQ
        Real(dp), Dimension(IP6) :: A, B, P, Q

        Character(Len=256) :: strfmt
        
        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),(A(1),PQ(2*IP6+1)),(B(1),PQ(3*IP6+1))
        
        A=0_dp
        B=0_dp
        P=0_dp
        Q=0_dp

        iw=40
        strfmt = '(4X,"Fbas: tolerance: ",I2)'
        write(11,strfmt) iw
  
        kcg=0                 !### number of changes of key Kbas
        small=1.d-8
        ksin=2                !### ksin=1 => r**n; ksin=2 => sin(r)*r**(n-1)
  
        call V_core           !### Y(r)=V_core-Z/r
  
        if (Kfile.GT.0) call Init_File
  
        do nj=1,Norb          !### Norb new orbitals to be constructed
            ni=Ns+nj
            nnj=Nn(ni)
            llj=Ll(ni)
            jjj=Jj(ni)
            kkj=Kk(ni)
            gj=iabs(kkj)
            if(Am.EQ.0.d0) gj=dsqrt(kkj**2-(Z/cl)**2)
            
            do i=Nso+1,Nsp      !### search for this orbital
                i1=i              !### in the Qnl array
                qn=Qnl(i)
                n1=10*dabs(qn)
                l1=100*dabs(qn)-10*n1+small
                if (qn.LT.small) then
                    j1=2*l1-1
                else
                    j1=2*l1+1
                end if
                if (nnj.EQ.n1.AND.llj.EQ.l1.AND.jjj.EQ.j1) goto 110
            end do
            write(*,*)' Fbas error. In Qnl no orbit:',nnj,llj,jjj
            stop
110         qn=Qnl1(i1)
  
            if (qn.EQ.0.d0) then
                kbk=0
                nnk=nnj-1
                llk=llj
                jjk=jjj
            else
                kbk=Kbas(i1)
                nnk=10*dabs(qn)+small
                llk=100*dabs(qn)-10*nnk+small
                nbk=10000*dabs(qn)-1000*nnk-100*llk+small
                if (qn.LT.small) then
                    jjk=2*llk-1
                else
                    jjk=2*llk+1
                end if
            end if

            if (kbk.EQ.0) then                     !# parent orbital is used #
                iwx=iw
                kkk=-(jjk-2*llk)*((jjk+1)/2)
                strfmt = '(/" >> ",I3,": Forming",I3,A1,I2,"/2 from",I3,A1,I2,"/2")'
                write( *,strfmt) ni,nnj,let(llj+1),jjj,nnk,let(llk+1),jjk
                write(11,strfmt) ni,nnj,let(llj+1),jjj,nnk,let(llk+1),jjk
                nim=ni-1
                do i=1,nim
                    nk=i
                    if (nnk.EQ.Nn(nk).AND.Kk(nk).EQ.kkk) goto 200
                end do
                strfmt = '(2X,"Error: no function n =",I3," k =",I3)'
                write( *,strfmt) nnk,kkk
                write(11,strfmt) nnk,kkk
                stop
200             call ReadF (12,nk+4,P,Q,2)
                call Orbit1(P,kkj,kkk,ksin)
            end if

            if (kbk.EQ.1.OR.kbk.EQ.2) then        ! B-spline is used 
                strfmt = '(/" >> ",I3,": Forming",I3,A1,I2,"/2 from "," B-spline ",I3,"[",I2,"/",I1,"]")'
                write( *,strfmt) ni,nnj,let(llj+1),jjj,nnk,nbk,llk
                write(11,strfmt) ni,nnj,let(llj+1),jjj,nnk,nbk,llk
                iwx=iw+5
                call Grid_Bspl(llk,nbk,m1,mx,R,ii)    !Forms mx grid nodes Ro from R and defines multiplier m1
                call FormBspl(llk,nbk,m1,mx,nnk,P,ii)
                
                P(ii+4)=gj
                call Origin(P,gj,kkj)             !### expansion at the origin
            end if

            if(kbk.EQ.3) then  
                iwx=iw-5   !# Extra file is used for P
                call OrbitF(P,nnk,llk,jjk,gj)  
            end if

            if (kbk.EQ.4) then      !# Extra file is used for P & Q
                iwx=iw-5                      
                call OrbitFG(P,Q,A,B,nnk,llk,jjk)
            else
                call Orbit2(P,Q,A,B,kkj,llj)         !### forming Q, P`, and Q`
            end if

            call WriteF(12,ni+4,P,Q,2)
            call WriteF(12,ni+4+Ns+Norb,A,B,2)
        
            call Ort(ni,Ns+Norb,iort)
            if (iort.GT.iwx) then
                write( *,*) ' Orthogonalization failed!'
                write(11,*) ' Orthogonalization failed!'
                if (kbk.EQ.0) then
                    stop
                else
                    write( *,*) ' Using kbk=0 to form this orbital...'
                    write(11,*) ' Using kbk=0 to form this orbital...'
                    Qnl1(i1)=0.d0
                    kcg=kcg+1
                    Ierr=Ierr-iort+1
                    goto 110
                end if
            end if
  
            lort=0
   210      if (iort.GE.1) then
                write( *,*) ' Reforming orbital ',ni,' iort=',iort
                write(11,*) ' Reforming orbital ',ni,' iort=',iort
                call Tail(ni,Ns+Norb,itail)      !### enforces that P(ii)=0
                call ReadF(12,ni+4,P,Q,2)
                call Orbit2(P,Q,A,B,kkj,llj)
                call WriteF(12,ni+4,P,Q,2)
                call WriteF(12,ni+4+Ns+Norb,A,B,2)
                call Ort(ni,Ns+Norb,iort)
                lort=lort+1
                if (lort.LT.3) goto 210
                strfmt = '("  Can not orthogonalize orbital ",I3,"  after ",I2," attempts")'
                write ( *,strfmt) ni,lort
                write (11,strfmt) ni,lort

                if (kbk.EQ.0) then
                    stop
                else
                    write( *,*) ' Using kbk=0 to form this orbital...'
                    write(11,*) ' Using kbk=0 to form this orbital...'
                    Qnl1(i1)=0.d0
                    kcg=kcg+1
                    Ierr=Ierr-iort+1
                    goto 110
                end if
            end if
        end do
        Ns=Ns+Norb
  
        if (3*IPs+22.GT.4*IP6) then
            write(*,*) ' Too many orbitals even for long basis variant!'
            write(*,*) ' Maximum number is ',(4*IP6-22)/3
            stop
        end if
  
        call ReadF (12,1,P,Q,2)           !### Saving Ns and quantum numbers
        call ReadF (12,3,A,B,2)
        PQ(2)=Ns
        longbasis=Ns.GT.(4*IP6-20)/6
  
        if (longbasis) then
            if (dabs(PQ(20)-0.98765d0).GT.1.d-6) then  ! PQ uses short basis format
                if=20                                    !# so we need to rewrite occ. numbers
                do ni=1,Ns                               ! reading occ. numbers from PQ
                    if=if+3
                    if (ni.LE.Ns-Norb) then
                        Qq1(ni)=PQ(if)
                    else
                        Qq1(ni)=0.d0
                    end if
                    if=if+3
                end do
            end if
  
            PQ(20)=0.98765d0                           ! - long basis mark
  
            do ni=1,Ns
              IQN(4*ni-3)=Nn(ni)
              IQN(4*ni-2)=Ll(ni)
              IQN(4*ni-1)=Kk(ni)
              IQN(4*ni)=Jj(ni)
            end do
        else                                         ! short basis format
            if=20
            do ni=1,Ns
                if=if+1
                PQ(if)=Nn(ni)
                if=if+1
                PQ(if)=Ll(ni)
                if=if+3
                PQ(if)=Kk(ni)
                if=if+1
                PQ(if)=Jj(ni)
            end do
        end if

        call WriteF(12,1,P,Q,2)
        call WriteF(12,3,A,B,2)
  
        if (kfile.GT.0) close(16)
  
        if (kcg.GT.0) then
            strfmt = '(/4X,"For ",I3," orbitals Kbas was changed to zero")'
            write( *,strfmt) kcg
            write(11,strfmt) kcg
        end if

        Return
    End Subroutine Fbas

    Subroutine Rot(idg)
        Implicit None

        Integer :: idg, nd1, nw, k, j
        Real(dp) :: dhf, zmax, trd, z_nd, zik, erx, djk
        Real(dp), Dimension(50) :: ez, dz
        Real(dp), Dimension(50,50) :: zz, hh
        Character(Len=256) :: strfmt
        
        nd1=50                    !### array dimension
        Emin= 1.d6                !### min(abs(e_i))
        Emax=-1.d6                !### max(abs(e_i))
        dhf=0.d0                  !### change of DF operator
        zmax=1.d-2                !### max |H_ik/(E_i-E_k)| gor idg=2
        trd=1.d-8                 !### convergence threshold
        nw=0
        irec=2*IP6*IPmr
        nmin=Nsv+1
        nmin=max(nmin,N_orb(n1,ch1,j1)) !### define diagonalization domain
                                        !### if kdg=0, then energies are
                                        !###    assigned as expectation values

        call Count_PW(nmin,nw)    !### counting partial waves
  
        open(12,file='HFD.DAT',status='OLD',access='DIRECT',recl=irec)
        open(13,file=FNAME,status='OLD',access='DIRECT',recl=irec)
        do ni=1,2*Ns+4
            call ReadF (12,ni,p,q,2)
            call WriteF (13,ni,p,q,2)
        end do
  
        do n=1,nw                    !### rotation of orbitals of
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
                strfmt = '(4X,A1,I2,"/2 partial wave. Dimension =",I2," orbitals found ",I2)'
                write(*, strfmt) ltr(n),jjw(n),nd,i
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

            strfmt = '(1X,"Eigen values for wave: (",A1,I2,"/2)",/(4E15.8))'
            if (Kout.GT.0) write(*, strfmt) ltr(n),jjw(n),(ez(i),i=1,nd)
            write(11,strfmt) ltr(n),jjw(n),(ez(i),i=1,nd)
            if (idg.NE.0) then
                strfmt = '(1X,"Diagonal coefficients:",/(7F10.5))'
                if (Kout.GT.0) write(*, strfmt) (dabs(zz(i,i)),i=1,nd)
                write(11,strfmt) (dabs(zz(i,i)),i=1,nd)
                do j=2,nd                         !### test of unitarity
                    do k=1,j-1
                        djk=0.d0
                        do i=1,nd
                            djk=djk+zz(i,j)*zz(i,k)
                        end do
                        if (dabs(djk).GT.1.d-8) then
                            strfmt = '(1X,"Diagonalization error: <",I3,"|",I3,"> =",E10.3)'
                            write( *,strfmt) j,k,djk
                            write(11,strfmt) j,k,djk
                            read(*,*)
                        end if
                    end do
                end do
            end if
  
            ! formation of new orbitals is done even for idg=0 to assign energies
            call New_Orb(nd,nd1,zz,ez,dhf)  
        end do                               
  
        strfmt = '(4X,"E_min = ",E10.3,"; E_max = ",E10.3,"; Kkin=",I2)'
        write( *,strfmt) Emin,Emax,Kkin
        write(11,strfmt) Emin,Emax,Kkin
  
        if (nmin.LE.Nso.AND.idg.GT.0) then
            strfmt = '(4X,"Core orbitals have changed by ",E10.3," Ierr=",I4)'
            write( *,strfmt) dhf,Ierr
            write(11,strfmt) dhf,Ierr

            if (dhf.GT.1.d-6) then
                Ierr=Ierr+20*dlog(dhf/1.d-6)
            else
              Khf=1
              strfmt = '(4X,"Self-consistency is reached. Push...")'
              write( *,strfmt)
              write(11,strfmt)
            end if
        end if

        Return
    End Subroutine Rot

    Subroutine Count_PW(nmin,nw)    !### counts partial waves
        Implicit None

        Integer :: nmin, nw, ki, i1
        Character(Len=256) :: strfmt

        strfmt = '(4X,"Counting partial waves starting from orbital ",I3)'
        write( *,strfmt) nmin
        write(11,strfmt) nmin
        
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

        strfmt = '(3X,A1,I2,"/2","-wave: ",I2," orbitals")'
        write (*, strfmt) (ltr(n),jjw(n),numw(n),n=1,nw)
        write (11,strfmt) (ltr(n),jjw(n),numw(n),n=1,nw)

        Return
    End Subroutine Count_PW

    Subroutine Form_zz(nd,nd1,idg,zz)
        Implicit None

        Integer :: nd, nd1, idg, nx, ni, kx, nik, ind, k, nk
        Real(dp) :: t
        Real(dp), Dimension(nd1, nd1) :: zz

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
200         continue
            end do
        end do

        Return
    End Subroutine Form_zz

    Subroutine New_Orb(nd,nd1,zz,ez,dhf)   !### Orbitals after diagonalization
        Use diff, Only : Dif
        Implicit None

        Integer :: nd, nd1, n_num, ni, j, k, i0, i1, nk, i
        Real(dp) :: small, crit, dhf, c, sorb
        Real(dp), Dimension(nd1) :: ez
        Real(dp), Dimension(nd1, nd1) :: zz

        small=1.d-15                     !### cutoff parameter
        crit=0.9d0                       !### how to calculate derivatives
        n_num=0
        if (nb(1).LE.Nso) call Check_core(nd,nd1,zz,ez)  !### Search for ghosts in the core
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
                if (dabs(c).LT.small) cycle
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
            end do
            
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
                    if (dabs(c).LT.small) cycle
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
                end do 
            end if
          
            i0=1
            do i1=1,ii
                if (dabs(a(i1)).ne.0.d0) then
                    i0=i1
                    exit
                endif
            enddo

            if (a(i0).LT.0.d0) then !### convension: a(0)>0
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
            strfmt = '(4X,"New_Orb: num dif-n used for",I3," orbitals")'
            write( *,strfmt) n_num
            write(11,strfmt) n_num
        end if

        Return
    End Subroutine New_Orb

    Subroutine Check_core(nd,nd1,zz,ez)    !### Search for ghosts in the core
        Implicit None

        Integer :: nd, nd1, itr, i1, i, m, k
        Real(dp) :: zx, eghost, v
        Real(dp), Dimension(nd1) :: ez
        Real(dp), Dimension(nd1, nd1) :: zz

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
                strfmt = '(4X,"Transposing vectors",I3," and",I3,", e_ghost=",E13.5)'
                write( *,strfmt) i,m,eghost
                write(11,strfmt) i,m,eghost
                do k=1,nd
                    v=zz(k,i)
                    zz(k,i)=zz(k,m)
                    zz(k,m)=v
                end do
            end if
        end do

        if (itr.GT.0) then
            strfmt = '(4X,"Ghosts detected! ",I3," transpositions made.")'
            write( *,strfmt) itr
            write(11,strfmt) itr
            read(*,*)
        end if

        Return
    End Subroutine Check_core

    subroutine Hould(n,ndim,Ee,Dd,Zz)
        implicit none
        integer :: n, ii, k, j, l, i, jj, m1, ndim, im, im1, k1, l1, ifail
        real(dp) :: tol, f, g, b, r, em, h, hh, c, ei, s, p, di, eps
        real(dp), dimension(ndim) :: Ee, Dd
        real(dp), dimension(ndim,ndim) :: Zz
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!        Householder's method of diagonalization
!        D-array of eigenvalues
!        Z-matrix of eigenvectors
!        n-dimension  of matrix z
!        eps-criteria of diagonalization
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        Ifail=0
        tol=2.0d0**(-103)
        eps=2.0d0**(-24)
        if (n > 1) then
          do ii=2,n
            i=n-ii+2
            l=i-2
            f=Zz(i,i-1)
            g=0.0d0
            if (l > 0) then
              do k=1,l
                g=g+Zz(i,k)**2
              end do
            end if
            h=f**2+g
            if (g < tol) then
              Ee(i)=f
              h=0.0d0
              Dd(i)=h
              cycle
            end if
            l=l+1
            g=-dsign(dsqrt(h),f)
            Ee(i)=g
            h=h-f*g
            Zz(i,i-1)=f-g
            f=0.0d0
            do j=1,l
              Zz(j,i)=Zz(i,j)/h
              g=0.0d0
              do k=1,j
                g=g+Zz(j,k)*Zz(i,k)
              end do
              k1=j+1
              if (k1 <= l) then
                do k=k1,l
                  g=g+Zz(k,j)*Zz(i,k)
                end do
              end if
              ee(j)=g/h
              f=f+g*Zz(j,i)
            end do
            hh=f/(h+h)
            do j=1,l
              f=Zz(i,j)
              g=Ee(j)-hh*f
              Ee(j)=g
              Zz(j,1:j)=Zz(j,1:j)-f*Ee(1:j)-g*Zz(i,1:j)
            end do
            Dd(i)=h
          end do
        end if
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        Dd(1)=0.0d0
        Ee(1)=0.0d0
        do i=1,n
           l=i-1
           if (Dd(i) /= 0.0d0 .and. l /= 0) then
             do j=1,l
               g=0.0d0
               do k=1,l
                 g=g+Zz(i,k)*Zz(k,j)
               end do
               Zz(1:l,j)=Zz(1:l,j)-g*Zz(1:l,i)
             end do
           end if
           Dd(i)=Zz(i,i)
           Zz(i,i)=1.0d0
           if (l /= 0) then
             Zz(i,1:l)=0.0d0
             Zz(1:l,i)=0.0d0
           end if
         end do
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     p matrix is formed
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (n /= 1) then
          Ee(1:n-1)=Ee(2:n)
        end if
        Ee(n)=0.0d0
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     triadigonalisation is finished
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        b=0.0d0
        f=0.0d0
        do l=1,n
          jj=0
          h=(dabs(Dd(l))+dabs(Ee(l)))*eps
          if (b < h) b=h
          do m1=l,n
            m=m1
            em=dabs(Ee(m))
            if (em <= b) exit
          end do
          if (m /= l) then
            do while (dabs(Ee(l)) > b)
              if (jj == 30) then
                ifail=1
                return
              end if
              jj=jj+1
              g=Dd(l)
              p=(Dd(l+1)-g)/Ee(l)*0.5d0
              r=dsqrt(p**2+1.0d0)
              em=p+dsign(r,p)
              Dd(l)=Ee(l)/em
              h=g-Dd(l)
              l1=l+1
              Dd(l1:n)=Dd(l1:n)-h
              f=f+h
              p=Dd(m)
              c=1.0d0
              s=0.0d0
              im1=m-1
              if (im1 >= l) then
                do im=l,im1
                  i=im1+l-im
                  ei=Ee(i)
                  g=c*ei
                  h=c*p
                  if (dabs(p) >= dabs(ei)) then
                    c=ei/p
                    r=dsqrt(c**2+1.0d0)
                    Ee(i+1)=s*p*r
                    s=c/r
                    c=1.0d0/r
                  else
                    c=p/ei
                    r=dsqrt(c**2+1.0d0)
                    Ee(i+1)=s*ei*r
                    s=1.0d0/r
                    c=c/r
                  end if
                  di=Dd(i)
                  p=c*di-s*g
                  Dd(i+1)=h+s*(c*g+s*di)
                  do k=1,n
                    h=Zz(k,i+1)
                    ei=Zz(k,i)
                    Zz(k,i+1)=s*ei+c*h
                    Zz(k,i)=c*ei-s*h
                  end do
                end do
              end if
              Ee(l)=s*p
              Dd(l)=c*p
            end do
          end if
          Dd(l)=Dd(l)+f
        end do
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     eigenvalues in D; eigenvectors in Z
!     ordering of D and Z
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        i=0
        do while (i < n)
          i=i+1
          if (i >= n) return
          h=Dd(i)
          f=Dd(i+1)
          if (h <= f) cycle
          Dd(i)=f
          Dd(i+1)=h
          do k=1,n
             h=Zz(k,i)
             Zz(k,i)=Zz(k,i+1)
             Zz(k,i+1)=h
          end do
          if (i /= 1) i=i-2
        end do
        return
    end subroutine Hould

    subroutine I_Diag(nd,nd1,H,Z,D,X,res)
        ! Iterative diagonalization of diagonally-dominated matrix
        ! nd - real dimension
        ! nd1 - declared dimension
        ! H - initial matrix (not changed)
        ! Z - matrix of eigenvectors in columns
        ! D - eigenvalues
        ! X - reserved space
        ! res - max residue
        Implicit None

        Integer :: nd, nd1, nzx, k
        Real(dp) :: res, res1, s, e, zx, xk, a, a2, xd
        Real(dp), Dimension(nd1) :: D, X
        Real(dp), Dimension(nd1,nd1) :: H, Z

        Character(Len=256) :: strfmt

        if (Sf.LT.1.d-3) Sf=1.d0                 ! Sf - damping factor; Sf = 0 - no damping
   
        call Test_Z(nd,nd1,H,Z,D,res1)
   
        do n=1,nd
            s=0.d0                                 !### new eigenvalue:
            e=0.d0                                 !#### e_n=<X|Z_n>/<Z_n|Z_n>
            zx=0.d0
            nzx=0
            do k=1,nd
                xk=0.d0
                a=Z(k,n)
                a2=a*a
                if (a2.GT.zx) then
                    zx=a2
                    nzx=k                              !nzx - dominant projection
                end if
                s=s+a2
                do m=1,nd
                    if(m.NE.k) then
                        xk=xk+H(k,m)*Z(m,n)
                    else
                        xd=H(k,k)*a
                    end if
                end do
                X(k)=xk                              !### X=(H-Diag)*Z_n
                e=e+a*(xk+xd)
            end do
            e=e/s                                  !### new eigenvalue

            D(n)=e
          
            do k=1,nd                              !### new eigenvector
                if (k.EQ.nzx) then
                    X(k)=Z(k,n)
                else
                    X(k)=Sf*X(k)/(e-H(k,k)) + (1.d0-Sf)*Z(k,n)
                end if
            end do
          
            do m=1,n-1                             !### orthogonalization
                s=0.d0
                do k=1,nd
                    s=s+Z(k,m)*X(k)
                end do

                do k=1,nd
                    X(k)=X(k)-s*Z(k,m)
                end do
            end do
          
            s=0.d0                                 !### normalization
            do k=1,nd
                a=X(k)
                s=s+a*a
            end do
            s=1.d0/dsqrt(s)
          
            do k=1,nd
              Z(k,n)=s*X(k)
            end do
        end do
   
        call Test_Z(nd,nd1,H,Z,D,res)

        strfmt = '(4X,"I_Diag: residue was",E10.3," and now is",E10.3)'
        write( *,strfmt) res1,res
        write(11,strfmt) res1,res

        if (res/res1.GT.0.9d0) then
            read(*,*)
            stop
        end if

        Return
    End Subroutine I_Diag

    Subroutine Test_Z(nd,nd1,H,Z,D,x_max) 
        ! x_max - max residue for eigenvalue eq-n
        Implicit None

        Integer :: nd, nd1, n, k
        Real(dp) :: x, x_max, xk, dxk
        Real(dp), Dimension(nd1) :: D
        Real(dp), Dimension(nd1,nd1) :: H, Z

        Character(Len=256) :: strfmt

        x_max=0.d0
        do n=1,nd
            x=0.d0
            strfmt = '(4X,"Eivenvector ",I3," E=",E15.8)'
            write(11,strfmt) n,D(n)
            do k=1,nd
                xk=0.d0
                do m=1,nd
                  xk=xk+H(k,m)*Z(m,n)
                end do
                dxk=(xk-D(n)*Z(k,n))
                x=x+dxk*dxk
            end do
            x=dsqrt(x)
            if (x_max.LT.x) x_max=x
            strfmt = '(4X,"Residue for n=",I3," is ",E12.5)'
            write(11,strfmt) n,x
        end do

        Return
    End Subroutine Test_Z

    Subroutine Init_File
        Implicit None

        Integer :: err_stat, kt3, i0, k0, k, k1, ia
        Real(dp) :: c1, z3, rnu3
        Integer, Dimension(4*IPs) :: Iqn
        Real(dp), Dimension(IP6) :: P, Q, P1, Q1
        Real(dp), Dimension(4*IP6) :: PQ
        Real(dp), Dimension(IPs) :: Qq1
        Character(Len=256) :: err_msg, strfmt
        Logical :: longbasis

        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+22))
        equivalence (P(1),PQ(1)),(Q(1),PQ(IP6+1)),(P1(1),PQ(2*IP6+1)),(Q1(1),PQ(3*IP6+1))
  
        write(*,*) ' Give file name for additional orbitals: '
        read(*,'(A12)') FNAME3
        open(16,file=FNAME3,access='DIRECT',status='OLD',recl=2*IP6*IPmr,iostat=err_stat,iomsg=err_msg)
        If (err_stat /= 0) Then
            write( *,*) ' No file ',FNAME3
            write(11,*) ' No file ',FNAME3
            stop
        End If
  
        c1=0.01d0
  
        call ReadF (16,1,P,Q,2)
        call ReadF (16,2,R3,V3,2)
        call ReadF (16,3,P1,Q1,2)
  
        z3  =PQ(1)
        Ns3 =PQ(2)+C1
        ii3 =PQ(3)+C1
        kt3 =PQ(9)+C1
        rnu3=PQ(13)
        strfmt = '(/4X,A12," opened. Z = ",F6.2,5X,"Kt =",I3,7X,"II =",I4,/4X,"Rnuc =",E11.4," Ns =",I3)'
        write( *,strfmt) FNAME3,z3,kt3,ii3,rnu3,Ns3
        write(11,strfmt) FNAME3,z3,kt3,ii3,rnu3,Ns3
  
        ierr=iabs(ii3-ii)+iabs(kt3-kt)
        intrpl=ierr.NE.0.OR.Rnuc.NE.rnu3
        if (intrpl) then
            strfmt = '(4X,"Different grid:",/4X,"ii =",I4,", ",I4, &
                       ";   kt =",I2,", ",I2,/4X,"Rnuc =",2E12.4, &
                       /4X,"R2   =",2E12.4,//4X,"Start interpolation:")'
            write( *,strfmt) ii,ii3,kt,kt3,Rnuc,rnu3,R(ii),R3(ii3)
            write(11,strfmt) ii,ii3,kt,kt3,Rnuc,rnu3,R(ii),R3(ii3)

            if (R(1).LT.R3(1)) then
                do i=1,ii
                    i0=i
                    if (R(i).GE.R3(1)) then
                        goto 100
                    else
                        Nintrpl(i)=0
                    end if
                end do
            else
                i0=1
            end if
          
100         write ( *,*) ' First ',i0-1,' nodes are skipped'
            write (11,*) ' First ',i0-1,' nodes are skipped'
            k0=1
            do i=i0,ii
                do k=k0,ii3
                    k1=k
                    if (R(i).LT.R3(k)) goto 200
                end do
                k0=ii3
                Nintrpl(i)=9999
                goto 300
200             k0=k1
                Nintrpl(i)=k1
300             write ( *,*) 'node ',i,' <-> ',Nintrpl(i)
                write (11,*) 'node ',i,' <-> ',Nintrpl(i)
            end do
        end if
  
        longbasis=dabs(PQ(20)-0.98765d0).LT.1.d-6
  
        if (longbasis) then
            do ni=1,Ns3
              Nn3(ni)=IQN(4*ni-3)
              Ll3(ni)=IQN(4*ni-2)
              Jj3(ni)=IQN(4*ni)
            end do
        else
            ia=20
            do ni=1,Ns3
              ia=ia+1
              Nn3(ni)=PQ(ia)+c1
              ia=ia+1
              Ll3(ni)=PQ(ia)+c1
              ia=ia+4
              Jj3(ni)=PQ(ia)+dsign(c1,PQ(ia))
            end do
        end if

        Return
    End Subroutine Init_File

    Subroutine Orbit1 (P,kkj,kkk,ksin)
        Use diff, Only : Origin
        ! upper component of the new orbital is made from the old one in
        ! the following way. ksin=1: P~=r^n * P
        !                    ksin=2: P~=r^(n-1) * r_max/pi*sin(pi*r/r_max) * P
        Implicit None

        Integer :: kkj, kkk, ksin, icase, np
        Real(dp) :: g, gj, ct, si, rrn
        Real(dp), Dimension(IP6) :: P
        Character(Len=3) :: str1(2)

        data str1 /' r ','sin'/

        g=P(ii+4)
        gj=iabs(kkj)
        if(Am.EQ.0.d0) gj=dsqrt(kkj**2-(Z/cl)**2)
        P(ii+4)=gj
        icase=gj-g+1.d-5
        if (kkj*kkk.LT.0) icase=icase+1
        if (icase.LT.0) Then
            write(*,*) ' Orbit: this case is not coded...'
            stop
        End If
        np=icase
        ct=3.1415926526d0/R(Ii)
        do i=1,ii                     !### large component P
            if (ksin.EQ.1) then
                si=R(i)
            else
                si=dsin(ct*R(i))/ct
            end if
            rrn=R(i)**np
            P(i)=rrn*si*P(i)
        end do
        call Origin(P,gj,kkj)         !### expansion at the origin
        strfmt = '(4X,"(r**",I2,"*",A3,")")'
        write(11,strfmt) icase,str1(ksin)

        Return
    End Subroutine Orbit1

    Subroutine Orbit2 (P,Q,A,B,kkj,llj)
        Use diff, Only : Origin, Dif
        Use test_ori
        ! lower component of the new orbital is expressed as follows:
        !    kkin=0:     from Dirac equation with  E=0 and V=0
        !    kkin=1,2:   from Dirac equation with  E=0 and V=Z/r
        Implicit None

        Integer :: kkj, llj, i, ir, ir1, ir2, ir3, ir4, ierr
        Real(dp) :: gj, ee, e_nr
        Real(dp), Dimension(IP6) :: P, Q, A, B
        Character(Len=7) :: str2(3)
        Character(Len=256) :: strfmt

        Data str2 /'   0   ','V_c-Z/r','V_c-Z/r'/

        gj=P(ii+4)
        call Dif(P,A,gj)              !### first derivative
        if (Kkin.EQ.2) then
            call NonRelE(P,A,llj,e_nr)  !### Non-relativistic energy
        else
            e_nr=0.d0
        end if
        P(ii+1)=-e_nr
        
        strfmt = '(4X,"Orbit2: V =",A7," E =",E12.5)'
        write( *,strfmt) str2(Kkin+1),e_nr
        write(11,strfmt) str2(Kkin+1),e_nr

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
    
        call Test_Origin('Orbit_2: P',P,R,MaxT,ii,1.d-6,ir1)
        call Test_Origin('Orbit_2: A',A,R,MaxT,ii,2.d-3,ir2)
        call Test_Origin('Orbit_2: Q',Q,R,MaxT,ii,1.d-5,ir3)
        call Test_Origin('Orbit_2: B',B,R,MaxT,ii,2.d-2,ir4)
        ir=ir1+ir2+ir3+ir4
        if(ir.GT.0) then
            strfmt = '(4X,"Orbit2: Bad matching at the origin.",//4X,"Errors for P, P`, Q, Q` =",4I4)'
            write(11,strfmt) ir1,ir2,ir3,ir4
            write( *,strfmt) ir1,ir2,ir3,ir4
            if (ir.GT.2) read(*,*)
            Ierr=Ierr+ir
        end if

        Return
    End Subroutine Orbit2

    Subroutine OrbitF(P,na,la,ja,gj)
        Implicit None

        Integer :: na, la, ja, ka, ind, ier, ni
        Real(dp) :: gj, gj1, scale, p1, rg
        Real(dp), Dimension(IP6) :: P, A
        Character(Len=13) :: lcase(3)
  
        data lcase /'interpolation','scaling','unchanged'/
  
        scale = 0_dp
        do ni=1,Ns3
            n=ni
            ier=iabs(na-Nn3(n))+iabs(la-Ll3(n))+iabs(ja-Jj3(n))
            if (ier.EQ.0) goto 200
        end do
        
        strfmt = '(4X,"Can not find orbital ",I3,A1,I2,"/2 in ",A12)'
        write( *,strfmt) na,let(la+1),ja,FNAME3
        write(11,strfmt) na,let(la+1),ja,FNAME3
        stop
  
200     call ReadF(16,n+4,P,P,1)
  
        if (intrpl) then
            call ReadF(16,n+Ns3+4,A,A,1)
            gj1=P(ii3+4)
            ka=((2*la-ja)*(ja+1))/2
            if (dabs(gj-gj1).GT.1.d-6) then
                write (*,*) ' Error: mismatch'
                write (*,*) ' gj=',gj,gj1
                write (*,*) ' kappa=',ka
                read(*,*)
                stop
            end if
            call Interpolation (P,A,gj,ka)
            ind=1
        else
            rg=R(1)**gj        !### scaling of the expansion at the origin
            p1=0.d0
            do i=ii+5,ii+5+MaxT
                p1=p1+P(i)
            end do
          
            if (dabs(scale-1.d0).GT.1.d-9) then
                ind=2
            else
                ind=3
            end if
          
            if (P(1)*p1.NE.0.d0) then
                scale=P(1)/(p1*rg)
                do i=ii+5,ii+5+MaxT
                    P(i)=P(i)*scale
                end do
            else
                if (dabs(p1-P(1)).GT.1.d-9) then
                    write(*,*) ' OrbitF warning for orbital ',na,' l=',la
                    write(*,*) ' P(1)=',P(1)
                    write(*,*) ' expansion at the origin gives ',p1
                    read(*,*)
                end if
            end if
        end if

        strfmt = '(/4X,"Orbital ",I3,A1,I2,"/2 is taken from ",A12," changes: ",A13)'
        write( *,strfmt) na,let(la+1),ja,FNAME3,lcase(ind)
        write(11,strfmt) na,let(la+1),ja,FNAME3,lcase(ind)

        Return
    End Subroutine OrbitF

    Subroutine OrbitFG(P,Q,A,B,na,la,ja)
        Use diff, Only : Dif
        Implicit None

        Integer :: na, la, ja, n, ier, ka
        Real(dp) :: gj
        Real(dp), Dimension(IP6) :: P, Q, A, B
  
        do ni=1,Ns3
            n=ni
            ier=iabs(na-Nn3(n))+iabs(la-Ll3(n))+iabs(ja-Jj3(n))
            if (ier.EQ.0) goto 200
        end do
        
        strfmt = '(4X,"Can not find orbital ",I3,A1,I2,"/2 in ",A12)'
        write( *,strfmt) na,let(la+1),ja,FNAME3
        write(11,strfmt) na,let(la+1),ja,FNAME3
        stop
  
200     call ReadF(16,n+4,P,Q,2)
        call ReadF(16,n+Ns3+4,A,B,2)
  
        if (intrpl) then
            gj=P(ii3+4)
            ka=((2*la-ja)*(ja+1))/2
            call Interpolation (P,A,gj,ka)
            call Dif(P,A,gj)
            call Interpolation (Q,B,gj,-ka)
            call Dif(Q,B,gj)
        end if
        strfmt = '(/4X,"Orbital ",I3,A1,I2,"/2 is taken from ",A12," (both P & Q)")'
        write( *,strfmt) na,let(la+1),ja,FNAME3
        write(11,strfmt) na,let(la+1),ja,FNAME3

        Return
    End Subroutine OrbitFG

    Subroutine interpolation(P,A,gj,kkj)
        Use diff, Only : Origin
        ! Interpolation formula between points a and b:
        !   F(x)=F(a) + d(x-a) + s(x-a)(x-b) + t(x-a)^2(x-b),
        ! where d=(F(b)-F(a))/(b-a) and
        ! parameters s and t are chosen to fit derivatives F'(a) and F'(b):
        ! F'(a)=d+s(a-b); F'(b)=d+s(b-a)+t(b-a)^2
        Implicit None

        Integer :: kkj, i, k1, k2
        Real(dp) :: gj, x, pa, pb, ra, rb, dr, d, da, db, s, t, rx
        Real(dp), Dimension(IP6) :: P, P1, A
        
        do i=1,ii
            k2=Nintrpl(i)
            k1=k2-1
            x=R(i)
            if (k2.NE.0.AND.k2.NE.9999) then
                pa=P(k1)
                pb=P(k2)
                ra=R3(k1)
                rb=R3(k2)
                dr=rb-ra
                d=(pb-pa)/dr
                da=A(k1)          ! = derivative d/dr
                db=A(k2)
                s=(d-da)/dr
                t=(db-d-s*dr)/dr**2
                P1(i)=P(k1) + d*(x-ra) + s*(x-ra)*(x-rb) + t*(x-ra)**2*(x-rb)
            else
                if (k2.EQ.9999) P1(i)=0.d0
                if (k2.EQ.0) then
                    da=A(1)
                    db=A(2)
                    rx=x-R3(1)
                    d=da+0.5d0*(db-da)/(R3(2)-R3(1))*rx  ! linear extrapolation
                    P1(i)=P(1)+d*rx                      ! for derivative
                end if
            end if
        end do
    
        do i=1,4
            P1(ii+i)=P(ii3+i)
        end do
    
        do i=1,ii+4
            P(i)=P1(i)
        end do
    
        call Origin(P,gj,kkj)
        
        Return
    End Subroutine Interpolation

    Subroutine Tail(ni, ns, itail)  
        ! forces Large component to zero at R(ii)
        ! ns - No of orbitals
        ! itail - =0 no changes & =1 when changes made
        Implicit None

        Integer :: ni, ns, itail, ix, ierr, i1
        Real(dp) :: small, c, x, pi, cpi, dx, rx, r0

        itail=0
        ix= 7                                   !### length of the tail
        small=1.d-6
        ierr=0
   
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
   
        strfmt = '(1X,"Tail(",I3,A1,I2,"/2): P(ii) old",E12.5," new",E12.5)'
        write( *,strfmt) nn(ni),let(ll(ni)+1),jj(ni),pi,P(ii)
        write(11,strfmt) nn(ni),let(ll(ni)+1),jj(ni),pi,P(ii)
        ierr=ierr+dlog(1.d0+dabs(pi)/small)
   
        call WriteF (12,ni+4,P,Q,2)
        call WriteF (12,ni+ns+4,CP,CQ,2)
   
        Return
    End Subroutine Tail

    Subroutine NonRelE(P,CP,l,e_nr)      !### Non-relativistic energy
        Implicit None

        Integer :: l, i
        Real(dp) :: e_nr, xl, s_kin, s_pot, s_nr
        Real(dp), Dimension(IP6) :: P, CP

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

        Return
    End Subroutine NonRelE

End Program bass