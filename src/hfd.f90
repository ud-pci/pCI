Program hfd
    Use params, Ncc => Nc, Nqparams => Nq, Nnparams => Nn, Llparams => Ll, Kkparams => Kk, Jjparams => Jj

    Implicit None

    Integer :: Mrec, iconf, MaxNc, is, ni, kt, ki, kf, n12, l, i, ii, nii, nit, iwr, iis, MaxT, &
                m1, m2, m3, ifin, kbr, n_is, Nmax, Imax, Nsmin, Nsmax, Ni0, Niter, ng, Lag, Klag
    Real(dp) :: eps0, eps1, e, d, del, r2, Rnucl, Gs, Es, Al, Bt, Cl, G0, R1, xja, xjc, Hs, Bts, Als
    Real(dp), Dimension(10) :: Vnuc
    Real(dp), Dimension(IP6) :: P, Q, A, B, C, CC, R, V, Y, CP, CQ, RO, ROO, W, ROIJ, RS, FS
    Real(dp), Dimension(IP6) :: Pa, Pc, Qa, Qc, P1a, P1c, Q1a, Q1c, Py, Qy
    Real(dp), Dimension(IPns) :: eadd
    Integer, Dimension(64) :: IGAM
    Real(dp), Dimension(64) :: GAM
    Integer, Dimension(IPns) :: Nc, Kp, Nf, Nq, Nn, Ll, Kk, Jj
    Real(dp), Dimension(IPns) :: Qq, Qw, Zat
    Character(Len=1) :: let(5), str(4)*8, str1(2)*5, name(16)
    Character(Len=256) :: strfmt
    
    let(1)='S'
    let(2)='P'
    let(3)='D'
    let(4)='F'
    let(5)='G'
    str(1)=' Volume '
    str(2)='Specific'
    str(3)=' Normal '
    str(4)='  Mass  '
    str1(1)='Gaunt'
    str1(2)='Breit'
    
    ! Initialize global arrays with zeros
    MaxT=0
    P=0_dp
    Q=0_dp
    A=0_dp
    B=0_dp
    C=0_dp
    CC=0_dp
    R=0_dp
    V=0_dp
    Y=0_dp
    CP=0_dp
    CQ=0_dp
    RO=0_dp
    ROO=0_dp
    ROIJ=0_dp
    RO=0_dp
     C=0_dp
    W=0_dp
    Pa=0_dp
    Pc=0_dp
    Qa=0_dp
    Qc=0_dp
    P1a=0_dp
    P1c=0_dp
    Q1a=0_dp
    Q1c=0_dp
    Py=0_dp
    Qy=0_dp
    eadd=0_dp
    igam=0
    gam=0_dp
    Nc=0
    Kp=0
    Nf=0
    Nq=0
    Nn=0
    Ll=0
    Kk=0
    Jj=0
    Qq=0_dp
    Qw=0_dp
    Zat=0_dp
    M1=0
    M2=0
    M3=0
    
    Call OpenFS('HFD.RES',11,1)
    Call recunit  ! determines word length
    Mrec=ipmr
    Call INPUT
    eps0=1.d-7         !### eps0 defines covergence criterion
    eps1=1.d-3 
    If (kl.EQ.2) eps0=0.1d0*eps0
    Call NUCL(0)
    Call INIT1

    iconf=0
    Do While (iconf <= maxnc)
        Call INIT2(iconf)
        Call INIT3(iconf)
        Do is=1,Ns
            ni=is
            If (NC(ni).NE.iconf) Cycle
            If (KP(ni).EQ.1) Cycle
            If (KL.EQ.1.AND.KP(ni).GE.0) Cycle
            If (KL.EQ.2.AND.KP(ni).EQ.0) Cycle
            KT=0
            KP(ni)=0
            KF=0
            N12=ni+4
            Call READF(12,N12,P,Q,2)
            Call NUCL(ni)
            If (KL.NE.3) Then
                Call FED(ni)
                Call NORM
            End If
            L=LL(ni)+1
            E=P(II+1)
            D=P(II+5)**2+Q(II+5)**2
            strfmt = '(I3,2X,I2,A1,I2,"/2 (",F6.3,") ",F12.6,3X,E13.6,4(2X,I3))'
            Write( *,strfmt) ni,NN(ni),LET(L),JJ(ni),QQ(ni),E,D,NIT,M1,M2,M3
            Write(11,strfmt) ni,NN(ni),LET(L),JJ(ni),QQ(ni),E,D,NIT,M1,M2,M3
            N12=ni+4
            Call WRITEF(12,N12,P,Q,2)
            Call READF(12,1,P,P,1)
            I=20+6*NS+2*(ni-1)+1
            P(I)=KP(ni)
            Call WRITEF(12,1,P,P,1)
        End Do
        
        ! THE BEGINING OF THE ITERATION PROCEDURE
        kout=0  !### output of Breit_Core
        ifin=0
        NITER=0
400     NITER=NITER+1
        Call ITER(KI,IFIN,iconf,DEL)
        Call DENS(0,iconf)
        If (ifin.EQ.1) GOTO 500
        Do is=1,Ns
            ni=is
            If (NC(ni).NE.iconf) Cycle
            If (KP(ni).EQ.1) Cycle
            KF=0
            n12=ni+4
            Call READF(12,n12,P,Q,2)
            If (NF(ni).EQ.1) GOTO 200
            Call Y0(ni)
            Call PLEX(ni)
            If (k_is.GE.2) Call sms_core(ni)
            If (kbr.GE.1)  Call Breit_Core(ni)
            Call NUCL(ni)
            Call FED(ni)
            Call NORM
            If (KF.EQ.0) GOTO 200
            Call Y0(ni)
            Call PLEX(ni)
            KF=2
            Call NUCL(ni)
            Call FED(ni)
            Call NORM
200         Call INTENS(ni,KI,DEL)
        End Do
        GOTO 400

        ! CALCULATION OF FIRST DERIVATIVES
500     If (iconf.EQ.0) eps0=0.1d0*eps0
        
        Do nii=1,Ns
            ni=nii
            If (Nc(ni).NE.iconf) Cycle
            N12=ni+4
            Call READF(12,N12,P,Q,2)
            If (kp(ni).NE.1) then
                Call Y0(ni)
                Call PLEX(ni)
                If (k_is.GE.2) Call sms_core(ni)
                If (kbr.GE.1)  Call Breit_Core(ni)
            end If
            Call PQDIF(ni)
        End Do

        If (iconf.EQ.0) Call Y0(0)
        iconf=iconf+1
    End Do

    If (k_is.GE.1) then
        strfmt = '(3X,A8," isotope shift * ",F7.4," added")'
        Write( *,strfmt) str(k_is),c_is
        Write(11,strfmt) str(k_is),c_is
    end If

    If (kbr.GE.1) then
        strfmt = '(3X,A5," exchange potential of the core (Nso=",I3,") added")'
        Write( *,strfmt) str1(kbr),nso
        Write(11,strfmt) str1(kbr),nso
    End If

    Call CLOSEF(11)
    Call CLOSEF(12)

Contains

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

    Subroutine OpenFS(fnam,kan,ntype)
        Implicit None

        Integer :: err_stat, kan, ntype
        Character*(*) :: fnam
        Character(Len=256) :: strfmt, err_msg
        
        ! ntype = 0 - old, ntype = 1 - new

        If (ntype.NE.0) Then
            Open(unit=kan,file=fnam,status='UNKNOWN',iostat=err_stat,iomsg=err_msg)
            If (err_stat /= 0) Then
                strfmt = '(/" NO FILE ",A)'
                Write( *,strfmt) fnam
                Write(11,strfmt) fnam
                Stop
            End If
        Else
            Open(unit=kan,file=fnam,status='OLD',iostat=err_stat,iomsg=err_msg)
            If (err_stat /= 0) Then
                strfmt = '(/" UNABLE TO Open FILE:",1X,A)'
                Write( *,strfmt) fnam
                Write(11,strfmt) fnam
                Stop
            End If
        End If

        Return
    End Subroutine OpenFS

    Subroutine CloseF(kan)
        Implicit None

        Integer :: kan

        Close(unit=kan)

        Return
    End Subroutine CloseF

    Subroutine WriteF(kan,record,V1,V2,nrec)
        Implicit None

        Integer :: kan, record, nrec, i, ii, nr1, nr2
        Real(dp) :: V1(IP6), V2(IP6)

        ii=IP6
        nr1=2*record-1
        nr2=nr1+1

        Write(kan,rec=nr1) (V1(i),i=1,ii)
        If (nrec.EQ.2) Write(kan, rec=nr2) (V2(i),i=1,ii)

        Return
    End Subroutine WriteF

    Subroutine ReadF(kan,record,V1,V2,nrec)
        Implicit None

        Integer :: kan, record, nrec, i, ii, nr1, nr2
        Real(dp) :: V1(IP6), V2(IP6)

        ii=IP6
        nr1=2*record-1
        nr2=nr1+1

        Read(kan,rec=nr1) (V1(i),i=1,ii)
        If (nrec.EQ.2) Read(kan, rec=nr2) (V2(i),i=1,ii)

        Return
    End Subroutine ReadF

    Subroutine OpenFD(fnam,kan,lr,nrec)
        Implicit None

        Integer :: err_stat, kan, lr, nrec, lrec, nr, nr1, nr2, nblo
        Character*(*) :: fnam
        Character(Len=256) :: strfmt, err_msg

        lrec=2*lr*Mrec
        nblo=0

        IF (nrec.NE.0) Then
            nr=IABS(NREC)
            nr1=nr/128
            nr2=nr-nr1*128
            nblo=lrec*nr1+(lrec*nr2-1)/128+1
            Open(unit=kan,file=fnam,status='UNKNOWN')
            Close(unit=kan,status='DELETE')
            Open(unit=kan,file=fnam,status='NEW',access='DIRECT',recl=lrec,iostat=err_stat,iomsg=err_msg)
            If (err_stat /= 0) Then
                strfmt = '(/" UNABLE TO ALLOCATE",I5,1X,"BLOCKS FOR FILE:",A)'
                Write( *,strfmt) nblo,fnam
                Write(11,strfmt) nblo,fnam
                Stop
            End If
        Else
            Open(unit=kan,file=fnam,status='OLD',access='DIRECT',recl=lrec,iostat=err_stat,iomsg=err_msg)
            If (err_stat /= 0) Then
                strfmt = '(/" NO FILE ",A)'
                Write( *,strfmt) nblo,fnam
                Write(11,strfmt) nblo,fnam
                Stop
            End If
        End If

        Return
    End Subroutine OpenFD

    Subroutine inpstr(istr)
        Implicit None

        Integer :: istr, i
        character*1 txt(5)
        character*128 string
        Character(Len=256) :: strfmt, strfmt2, strfmt3

        istr=0
        string(1:5)='     '
        txt(1)=' '
        txt(2)=' '
        txt(3)=' '
        txt(4)=' '
        txt(5)=' '
        n_is = 0

        strfmt = '(5a1,a)'
        read (10,strfmt) (txt(i),i=1,5),string

        strfmt = '(5a1,i4)'
        strfmt2 = '(5a1,f7.2)'
        strfmt3 = '(5a1,f7.4)'

        if (txt(2).ne.'k'.and.txt(2).ne.'K') goto 200
        if (txt(3).ne.'l'.and.txt(3).ne.'L') goto 200
        if (txt(4).ne.' ') goto 200
        read (string,*) kl
        write( *,strfmt) (txt(i),i=1,5),kl
        write(11,strfmt) (txt(i),i=1,5),kl
        return

 200    if (txt(2).ne.'n'.and.txt(2).ne.'N') goto 210
        if (txt(3).ne.'s'.and.txt(3).ne.'S') goto 210
        if (txt(4).ne.' ') goto 210
        read (string,*) ns
        write( *,strfmt) (txt(i),i=1,5),ns
        write(11,strfmt) (txt(i),i=1,5),ns
        return

 210    if (txt(2).ne.'n'.and.txt(2).ne.'N') goto 220
        if (txt(3).ne.'s'.and.txt(3).ne.'S') goto 220
        if (txt(4).ne.'o'.and.txt(4).ne.'O') goto 220
        read (string,*) nso
        write( *,strfmt) (txt(i),i=1,5),nso
        write(11,strfmt) (txt(i),i=1,5),nso
        return
        
 220    if (txt(2).ne.'z'.and.txt(2).ne.'Z') goto 230
        if (txt(3).ne.' '.and.txt(3).ne.' ') goto 230
        if (txt(4).ne.' '.and.txt(4).ne.' ') goto 230
        read (string,*) z
        write( *,strfmt2) (txt(i),i=1,5),z
        write(11,strfmt2) (txt(i),i=1,5),z
        return
        
 230    if (txt(2).ne.'a'.and.txt(2).ne.'A') goto 240
        if (txt(3).ne.'m'.and.txt(3).ne.'M') goto 240
        if (txt(4).ne.' ') goto 240
        read (string,*) am
        write( *,strfmt2) (txt(i),i=1,5),am
        write(11,strfmt2) (txt(i),i=1,5),am
        return
        
 240    if (txt(2).ne.'j'.and.txt(2).ne.'J') goto 250
        if (txt(3).ne.'m'.and.txt(3).ne.'M') goto 250
        if (txt(4).ne.' ') goto 250
        read (string,*) jm
        write( *,strfmt2) (txt(i),i=1,5),jm
        write(11,strfmt2) (txt(i),i=1,5),jm
        return
        
 250    if (txt(2).ne.'r'.and.txt(2).ne.'R') goto 260
        if (txt(3).ne.'2') goto 260
        if (txt(4).ne.' ') goto 260
        read (string,*) r2
        write( *,strfmt2) (txt(i),i=1,5),r2
        write(11,strfmt2) (txt(i),i=1,5),r2
        return

 260    if (txt(1).ne.'k'.and.txt(1).ne.'K') goto 270
        if (txt(3).ne.'i'.and.txt(3).ne.'I') goto 270
        if (txt(4).ne.'s'.and.txt(4).ne.'S') goto 270
        read (string,*) k_is
        write( *,strfmt) (txt(i),i=1,5),k_is
        write(11,strfmt) (txt(i),i=1,5),k_is
        return
        
 270    if (txt(1).ne.'c'.and.txt(1).ne.'C') goto 280
        if (txt(3).ne.'i'.and.txt(3).ne.'I') goto 280
        if (txt(4).ne.'s'.and.txt(4).ne.'S') goto 280
        read (string,*) c_is
        write( *,strfmt3) (txt(i),i=1,5),c_is
        write(11,strfmt3) (txt(i),i=1,5),c_is
        return
        
 280    if (txt(1).ne.'n'.and.txt(1).ne.'N') goto 290
        if (txt(3).ne.'i'.and.txt(3).ne.'I') goto 290
        if (txt(4).ne.'s'.and.txt(4).ne.'S') goto 290
        read (string,*) n_is
        write( *,strfmt) (txt(i),i=1,5),n_is
        write(11,strfmt) (txt(i),i=1,5),n_is
        return
        
 290    if (txt(2).ne.'l'.and.txt(2).ne.'L') goto 300
        if (txt(3).ne.'o'.and.txt(3).ne.'O') goto 300
        if (txt(4).ne.'w'.and.txt(4).ne.'W') goto 300
        read (string,*) klow
        write( *,strfmt) (txt(i),i=1,5),klow
        write(11,strfmt) (txt(i),i=1,5),klow
        return
        
 300    if (txt(2).ne.'k'.and.txt(2).ne.'K') goto 310
        if (txt(3).ne.'b'.and.txt(3).ne.'B') goto 310
        if (txt(4).ne.'r'.and.txt(4).ne.'R') goto 310
        read (string,*) kbr
        write( *,strfmt) (txt(i),i=1,5),kbr
        write(11,strfmt) (txt(i),i=1,5),kbr
        return
        
 310    if (txt(1).ne.'r'.and.txt(1).ne.'R') goto 320
        if (txt(2).ne.'n'.and.txt(2).ne.'N') goto 320
        if (txt(3).ne.'u'.and.txt(3).ne.'U') goto 320
        read (string,*) rnuc
        write( *,strfmt3) (txt(i),i=1,5),rnuc
        write(11,strfmt3) (txt(i),i=1,5),rnuc
        return
        
 320    if (txt(2).ne.' '.and.txt(2).ne.'-') goto 700
        if (txt(3).ne.' '.and.txt(3).ne.'-') goto 700
        if (txt(4).ne.' '.and.txt(4).ne.'-') goto 700
        backspace 10
        istr=1
        return
         
 700    istr=2
        strfmt = '(/2x,"Wrong string in input file: ",5a1)'
        write( *,strfmt) (txt(i),i=1,5)
        write(11,strfmt) (txt(i),i=1,5)
        Stop
        
    End Subroutine inpstr
    
    Subroutine Input
        Implicit None
        
        Integer :: istr, li, ji, Ns0, nmin, kpa, nca, na, la, ja, k, nj, ni, nii, i, iconf
        Real(dp) :: qa
        Integer, Dimension(IPns) :: Nn1, Ll1, Kk1, Jj1, KP1, Nc1
        Real(dp), Dimension(IPns) :: Qq1, Dqa
        Character*1 :: txt(36), lt
        Character(Len=256) :: strfmt

        ! INPUT DATA FROM FILE 'HFD.INP'
        Call OpenFS('HFD.INP',10,0)

        strfmt = '(1X,16A1)'
        Read(10,strfmt) NAME

        strfmt = '(/2X,"PROGRAM HFD (version with IS and Breit) v1.1",4X,16A1)'
        Write( *,strfmt) NAME
        Write(11,strfmt) NAME


        ! Default values:
        Kl=0
        Ns=IPns
        Nso=0
        Jm=-2_dp
        R2=50.d0
        k_is=0         !### 0 - No IS, 1 - VS, 2,3,4 - SMS,NMS,MS
        c_is=0.d0      !### prefactor for IS perturbation
        kbr=0          !### 0 - Coulomb, 1 - Gaunt, 2 - Breit
        rnuc=0.d0      !### 0 - Rnucl is defined by AM;
                       !### >0 - Rnucl=sqrt(5/3)*rnuc*fermi

        istr = 0
        Do While (istr /= 1)
            Call inpstr(istr)
        End Do
        
        If (n_is.EQ.0) n_is=Nso  !### defines SMS core potential

        strfmt = '(36A1)'
        READ (10,strfmt) TXT
        WRITE( *,strfmt) TXT
        WRITE(11,strfmt) TXT
        READ (10,strfmt) TXT
        WRITE( *,strfmt) TXT
        WRITE(11,strfmt) TXT
        READ (10,strfmt) TXT
        WRITE( *,strfmt) TXT
        WRITE(11,strfmt) TXT

        MAXNC=0
        Do NII=1,NS
            NI=NII
            KP1(NI)=0
            NC1(NI)=0
            strfmt = '(1X,6A1,I2,A1,1X,A1,I1,3A1,2X,F7.4,4X,I1,3X,I2)'
            READ (10,strfmt) (TXT(I),I=1,6),NN1(NI),LT,TXT(7),JJ1(NI),TXT(8),TXT(9),TXT(10),QQ1(NI),KP1(NI),NC1(NI)
            Do I=1,5
                LI=I-1
                IF (LT.EQ.LET(I)) GOTO 200
            End Do

            ns=nii-1
            goto 201
  200       WRITE( *,strfmt) (TXT(I),I=1,6),NN1(NI),LT,TXT(7),JJ1(NI),TXT(8),TXT(9),TXT(10),QQ1(NI),KP1(NI),NC1(NI)
            WRITE(11,strfmt) (TXT(I),I=1,6),NN1(NI),LT,TXT(7),JJ1(NI),TXT(8),TXT(9),TXT(10),QQ1(NI),KP1(NI),NC1(NI)

            LL1(NI)=LI
            JI=JJ1(NI)
            IF (JI.NE.2*LI+1.AND.JI.NE.2*LI-1) Then
                strfmt = '(/1X,"WRONG CONFIGURATION NI =",I3)'
                WRITE( *,strfmt) NI
                WRITE(11,strfmt) NI
                STOP
            End If 
            IF (JI.EQ.2*LI+1) KK1(NI)=-(LI+1)
            IF (JI.EQ.2*LI-1) KK1(NI)=LI
            IF (NC1(NI).GT.MAXNC) MAXNC=NC1(NI)
        End Do

  201   Do ICONF=0,MAXNC
            DQA(ICONF+1)=0.D0
        End Do
        ! END OF INPUT

        ! UNPACKING OF ARRAY QNL(NI)
        NS0=NSO
        IF (JM.GT.-1.1D0) NS0=NS

        IF (NS0.EQ.0) GOTO 210
        Do NI=1,NS0
            NN(NI)=NN1(NI)
            LL(NI)=LL1(NI)
            KK(NI)=KK1(NI)
            JJ(NI)=JJ1(NI)
            KP(NI)=KP1(NI)
            NC(NI)=NC1(NI)
            QQ(NI)=QQ1(NI)
            NQ(NI)=QQ(NI)+0.01D0
            ICONF=NC(NI)
            DQA(ICONF+1)=DQA(ICONF+1)+QQ(NI)
            QW(NI)=1.D0
        End Do

210     If (NS0.EQ.NS) GOTO 1000
        NS0=NS
        NS=NSO
        NMIN=NSO+1
        Do NI=NMIN,NS0
            KPA=KP1(NI)
            NCA=NC1(NI)
            NA=NN1(NI)
            LA=LL1(NI)
            JA=JJ1(NI)
            QA=QQ1(NI)
            DQA(NCA+1)=DQA(NCA+1)+QA
            K=0
            IF (NMIN.GT.NS) GOTO 220
            Do NJ=NMIN,NS
                IF (NN(NJ).NE.NA) Cycle
                IF (LL(NJ).NE.LA) Cycle
                IF (NC(NJ).NE.NCA) Cycle
                QW(NJ)=QW(NJ)+QA
                IF (JJ(NJ).EQ.JA) KP(NJ)=KPA
                K=1
            End Do
    220     IF (K.EQ.1) Cycle
            NS=NS+1
            NN(NS)=NA
            LL(NS)=LA
            KP(NS)=KPA
            NC(NS)=NCA
            QW(NS)=QA
            IF (LA.GT.0) GOTO 230
            JJ(NS)=1
            KK(NS)=-1
            Cycle
    230     JJ(NS)=LA+LA-1
            KK(NS)=LA
            NS=NS+1
            NN(NS)=NA
            LL(NS)=LA
            KP(NS)=KPA
            NC(NS)=NCA
            QW(NS)=QA
            JJ(NS)=LA+LA+1
            KK(NS)=-LA-1
        End Do

        Do NI=NMIN,NS
            QQ(NI)=QW(NI)*(JJ(NI)+1.D0)/(4*LL(NI)+2)
            QW(NI)=QW(NI)/(4*LL(NI)+2)
            NQ(NI)=QQ(NI)+0.01D0
        End Do

1000    Do ICONF=0,MAXNC
            I=ICONF+1
            ZAT(I)=Z-DQA(I)
            IF (ICONF.GT.0) ZAT(I)=ZAT(I)-DQA(1)
        End Do
        CALL CLOSEF(10)
        RETURN

    End Subroutine Input


    Subroutine Nucl(ni)
        Implicit None

        Integer :: ni, k, n, j, m
        Real(dp) :: fermi, r0, gam, t, t1, t2, s, s1, s2, e, d

        Real(dp), Dimension(11) :: Yy, P1, Q1, Vp, Vq, Xp, Xq

        Yy=0_dp
        P1=0_dp
        Q1=0_dp
        Vp=0_dp
        Vq=0_dp
        Xp=0_dp
        Xq=0_dp

        CL=DPcl
        NMAX=9

        ! CALCULATION OF THE NUCLEAR RADIUS
        RNUCL=0.d0

        If (AM.LT.1.d0) Then
            Continue
        Else
            FERMI=1.d-13/DPrb
            D=AM**0.333333333d0
            RNUCL=(1.115d0*D+2.151d0/D-1.742d0/AM)*FERMI
        End If

        If (rnuc.GT.0.d0) Then
            rnucl=dsqrt(5.d0/3.d0)*rnuc*fermi
        End If

        R1=DEXP(-4.d0)/Z
        IF (RNUCL.GT.0.d0) R1=RNUCL

        Vnuc(1:Nmax+1) = 0.d0

        If (AM.LT.1.d0) Then
            VNUC(1)=-Z
        Else
            VNUC(2)=-1.5d0*Z
            VNUC(4)=0.5d0*Z
            If (k_is.EQ.1) then
                ! adding isotope shift perturbation: !### These lines account for
                VNUC(2)=VNUC(2)*(1.d0-c_is)          !### volume shift. Note, that
                VNUC(4)=VNUC(4)*(1.d0-3*c_is)        !### c_is=dRnucl/Rnucl
            End If
        End If

        IF (NI.GT.0) Then
            Continue
        Else
            G0=DSQRT(1.d0-(VNUC(1)/CL)**2)
            Return
        End If

        R0=R1
        K=KK(NI)
        E=P(II+1)
        GAM=DSQRT(K**2-(VNUC(1)/CL)**2)

        Do N=0,NMAX
            I=N+1
            YY(I+1)= Y(II+I)*R0
            XP(I+1)=CP(II+I)*R0
            XQ(I+1)=CQ(II+I)*R0
        End Do
        YY(1)=0.d0
        XP(1)=0.d0
        XQ(1)=0.d0
        Do N=0,NMAX
            I=N+1
            VP(I)=VNUC(I)+YY(I)
            VQ(I)=VNUC(I)+YY(I)
        End Do
        VP(2)=VP(2)+(E-2*CL*CL)*R0
        VQ(2)=VQ(2)+E*R0

        T=IABS(K)
        S=GAM+T
        D=DSQRT((P(II+5)**2+Q(II+5)**2)*0.5d0*S/(T*P(II+2)))
        T=D*DABS(VNUC(1))/(S*CL)

        If (K.GT.0) Then
            P1(1)=T
            Q1(1)=-D
        Else
            P1(1)=D
            Q1(1)=T
        End If
        Do N=1,NMAX
            I=N+1
            D=N*(N+2*GAM)*CL*CL
            T1=(N+GAM-K)*CL
            T2=(N+GAM+K)*CL
            S1=XQ(I)
            S2=XP(I)
            Do M=1,N
                J=M+1
                S1=S1+P1(I-M)*VQ(J)
                S2=S2+Q1(I-M)*VP(J)
            End Do
            P1(I)=(-VP(1)*S1+T1*S2)/D
            Q1(I)=(-VQ(1)*S2-T2*S1)/D
        End Do
        P(II+4)=GAM
        Q(II+4)=GAM 	

        Do N=0,NMAX
            I=N+1
            P(II+5+N)=P1(I)
            Q(II+5+N)=Q1(I)
        End Do

        Return
    End Subroutine Nucl

    Subroutine Init1
        Implicit None

        Integer :: lrec, kpn, nrec, if, ih, ms, max, i, ni, nii, n12
        Real(dp) :: c1, rm, d0, rmax, d
        Character(Len=256) :: strfmt

        c1=0.01D0

        KLAG=1
        NG=0
        KT=0
        RMAX=DABS(R2)
        IMAX=IP6-17              !### = 239 for IP6=256
        IMAX=((IMAX+1)/2)*2-1
        IF (R2.LE.0.D0) Then
            ! LOGARITHMIC GRID (FISHER)
            H=1.D0/32.D0      ! this parameter is linked to number of grid points
            AL=0.D0
            BT=1.D0
            II=DLOG(RMAX/R1)/H+1
            II=((II+1)/2)*2-1
            If (II.LE.IMAX) Then
                Continue
            Else
                strfmt = '(/"  II > IMAX"/"  IMAX=",I3/"  II  =",I3/"  H   =",F6.3, &
                            "  BT  =",F6.3/"  R1  =",E10.3,"  R2  =",F7.3)'
                WRITE( *,strfmt) IMAX,II,H,BT,R1,RMAX
                STOP
            End If
        Else
            ! GRID RO=AL*R+BT*LN(R)
            II=IMAX
            H=1.D0/32.D0       ! this parameter is linked to number of grid points
            BT=1.D0
            AL=(H*(IMAX-1)-BT*DLOG(R2/R1))/(R2-R1)
            IF (AL.GE.0.D0) Then
                Continue
            Else
                II=(BT*DLOG(R2/R1))/H+1
                strfmt = '(/"  II > IMAX"/"  IMAX=",I3/"  II  =",I3/"  H   =",F6.3, &
                            "  BT  =",F6.3/"  R1  =",E10.3,"  R2  =",F7.3)'
                WRITE( *,strfmt) IMAX,II,H,BT,R1,R2
                STOP
            End If
        End If

        strfmt = '(/"  NSP =",I4/"  II  =",I4/"  R1  =",E12.5/"  H   =",F7.4/"  AL  =",F7.4,"  BT =",F5.2)'
        WRITE( *,strfmt) NS,II,R1,H,AL,BT
        WRITE(11,strfmt) NS,II,R1,H,AL,BT

        IF (JM.LT.0.D0.AND.JM.GT.-1.1D0) Then
            strfmt = '(/2X,"RELATIVISTIC CONFIGURATION AVERAGE")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
        End If

        IF (JM.LE.-1.1D0) Then
            strfmt = '(/2X,"NONRELATIVISTIC CONFIGURATION AVERAGE")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
        End If

        IF (JM.GE.0.D0) Then
            strfmt = 'FORMAT(/2X,"TERM   J =",F5.1)'
            WRITE( *,strfmt) JM
            WRITE(11,strfmt) JM
        End If

        ! FILE HFD.DAT
        ! key KL defines initial approximation
        ! KL = 0 - calculation starts from beginning using standard initial approximation
        ! KL = 1 - continuation of the previous calculation (used when self-consistency has not been reached)
        ! KL = 2 - initial approximation is taken from pre-existing HFD.DAT
        ! KL = 3 - all orbitals in the initial approximation are set to 0
        IF (KL.EQ.1.OR.KL.EQ.2) Then
            LREC=IP6 ! record length
            CALL OPENFD('HFD.DAT',12,LREC,0)
            CALL TEST
        Else
            ! key KP defines where the orbital is found
            ! KP = 0 means that orbital is found from HFD equations and written to HFD.DAT
            ! KP = 1 means that the orbital is taken from pre-existing HFD.DAT and kept frozen
            
            ! Check if any orbitals have to be calculated from HFD equations
            KPN=0
            Do NI=1,NS
                If (KP(NI).EQ.1) Then
                    KPN=1
                Else
                    KP(NI)=-1
                End If
            End Do

            ! If any orbitals are taken from pre-existing HFD.DAT, test them
            IF ((KPN.EQ.1).OR.(JM.GE.0)) Then
                LREC=IP6
                CALL OPENFD('HFD.DAT',12,LREC,0)
                CALL TEST
            Else
                LREC=IP6
                NREC=NS+NS+5
                CALL OPENFD('HFD.DAT',12,LREC,NREC)
            End If
        End If

        IF (KL.NE.1) Then
            Continue
        Else
            ! KEY KL=1
            CALL READF (12,2,R,V,2)
            Return
        End If

        IF (KL.NE.2.AND.KPN.EQ.0) Then
            Continue
        Else
            ! OLD GRID
            CALL READF (12,1,CQ,CQ,1)
            CALL READF (12,2,CP,CP,1)
        End If

        NI0=1
        CALL TBR
        P(1) =Z
        P(2) =NS
        P(3) =II
        P(4) =R1
        P(5) =R2
        P(6) =H
        P(7) =BT
        P(8) =AL
        P(9) =KT
        P(10)=NG
        P(11)=NSO
        P(12)=AM
        P(13)=RNUCL
        P(14)=JM
        P(15)=NI0
        P(17)=MAXNC
        IF=20
        Do NI=1,NS
            IF=IF+1
            P(IF)=NN(NI)
            IF=IF+1
            P(IF)=LL(NI)
            IF=IF+1
            P(IF)=QQ(NI)
            IF=IF+1
            P(IF)=QW(NI)
            IF=IF+1
            P(IF)=KK(NI)
            IF=IF+1
            P(IF)=JJ(NI)
        End Do
        Do NI=1,NS
            IF=IF+1
            P(IF)=KP(NI)
            IF=IF+1
            P(IF)=NC(NI)
        End Do
        CALL WRITEF(12,1,P,P,1)
        CALL WRITEF(12,2,R,V,2)
        IF (KL.EQ.1) Return

        ! INTERPOLATION
        IH=2-KT
        IIS=CQ(3)+C1
        HS =CQ(6)
        BTS=CQ(7)
        ALS=CQ(8)
        Do NII=1,NS
            NI=NII
            IF (KP(NI).EQ.-1) Cycle
            IF (KL.EQ.0.AND.KP(NI).NE.1) Cycle
            IF (KL.EQ.3.AND.KP(NI).NE.1) Cycle
            N12=NI+4
            CALL READF (12,N12,P,Q,2)
            GS  =P(IIS+4)
            ES  =P(IIS+1)
            IMAX=P(IIS+3)+C1
            MS  =P(IIS+15)+C1
            Do I=1,IP6
                CQ(I)=P(I)
            End Do
            M1=MS
            RM=CP(MS)
            D0=DABS(RM-R(M1))
            Do I=1,II,IH
                P(I)=FLAGR(R(I))
                IF (I.EQ.2*(I/2)) Cycle
                D=DABS(RM-R(I))
                IF (D.GE.D0) Cycle
                M1=I
                D0=D
            End Do
            MAX=8+NMAX
            Do I=1,MAX
                P(II+I)=CQ(IIS+I)
            End Do
            P(II+3)=II
            P(II+15)=M1
            Do I=1,IP6
                CQ(I)=Q(I)
            End Do
            Do I=1,II,IH
                Q(I)=FLAGR(R(I))
            End Do
            MAX=8+NMAX
            Do I=1,MAX
                Q(II+I)=CQ(IIS+I)
            End Do
            CALL WRITEF(12,N12,P,Q,2)
        End Do

        Return
    End Subroutine Init1

    Real(dp) Function FLAGR(X)
        Implicit None

        Integer :: ih, ih2, i1, i2, i3, i4, i5, i6
        Real(dp) :: x, ro, pm1, pm2, pm3, pp1, pp2, fi, p, h, d

        IF (X.GT.CP(IMAX)) Then
            D=DSQRT(ES)*(X-CP(IMAX))
            FLAGR=0.D0
            IF (D.LT.30.D0) FLAGR=CQ(IMAX)*DEXP(-D)
            RETURN
        End If

        IF (X.LT.CP(   1)) Then
            FLAGR=X**(GS)*(CQ(II+5)+CQ(II+6)+CQ(II+7)+CQ(II+8))
            RETURN
        End If

        IH=2-KT
        H=IH*HS
        IH2=IH+IH
        RO=AL*(X-CP(1))+BT*DLOG(X/CP(1))
        I3=1+IH2

        If (X.LE.CP(I3)) Then
            Continue
        Else
            I3=IMAX-(IH2+IH)
            If (X.GE.CP(I3)) Then
                Continue
            Else
                I3=RO/HS+1
                If (IH.EQ.1) Then
                    Continue
                Else
                    I3=((I3-1)/2)*2+1
                End If
            End If
        End If
        
        P=(RO-(I3-1)*HS)/H
        I2=I3-IH
        I1=I2-IH
        I4=I3+IH
        I5=I4+IH
        I6=I5+IH
        PM1=P-1.D0
        PM2=P-2.D0
        PM3=P-3.D0
        PP1=P+1.D0
        PP2=P+2.D0
        FI=0.1D0*PM2*PM1*P*PP1*(PP2*CQ(I6)-PM3*CQ(I1))
        FI=FI+0.5D0*PM3*PM1*P*PP2*(PM2*CQ(I2)-PP1*CQ(I5))
        FI=FI+PM3*PM2*PP1*PP2*(P*CQ(I4)-PM1*CQ(I3))
        FLAGR=FI/12.D0
        
        Return
    End Function FlagR

    Subroutine Init2(iconf)
        Implicit None

        Integer :: iconf, ni, j, n, nii, n12, l
        Real(dp) :: c1, d0, qi, dn, e
        Character*1 LM(2)
        Character(Len=256) :: strfmt
        DATA LM/' ','*'/

        strfmt = '(/2X,"ICONF =",I2,4X,"ZAT =",F7.3)'
        WRITE( *,strfmt) ICONF,ZAT(ICONF+1)
        WRITE(11,strfmt) ICONF,ZAT(ICONF+1)

        NSMIN=NS
        Do NI=1,NS
            IF (NC(NI).NE.ICONF) Cycle
            IF (NI.LT.NSMIN) NSMIN=NI
            NSMAX=NI
        End Do

        ! KEY KL=1
        Do NI=1,NS
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).NE.-1) GOTO 200
        End Do
        Return

200     C1=0.01D0

        strfmt = '(/2X,"INITIAL APPROXIMATION FROM DISK"/)'
        WRITE( *,strfmt)
        WRITE(11,strfmt)

        strfmt = '(10X,"JJ",4X,"QQ",4X,"KP",7X,"EE",9X,"D(0)",11X,"NORM",10X,"M1")'
        WRITE( *,strfmt)

        strfmt = '(10X,"JJ",4X,"QQ",4X,"KP",7X,"EE",9X,"D(0)",11X,"NORM")'
        WRITE(11,strfmt)

        Do NII=1,NS
            NI=NII
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).EQ.-1) Cycle
            N12=NI+4
            CALL READF (12,N12,P,Q,2)
            M1=P(II+15)+C1
            D0=P(II+5)**2+Q(II+5)**2
            L =LL(NI)+1
            J =JJ(NI)
            N =NN(NI)
            QI=QQ(NI)
            E =P(II+1)
            DN=P(II+2)
            M=1
            IF (NQ(NI).NE.J+1) M=2
            strfmt = '(I4,1X,I2,A1,I2,"/2"," (",F6.3,")",A1,1X,I1,F13.6,2X,E13.6,2X,E13.6,I6)'
            WRITE( *,strfmt) NI,N,LET(1),J,QI,LM(M),KP(NI),E,D0,DN,M1
            WRITE(11,strfmt) NI,N,LET(L),J,QI,LM(M),KP(NI),E,D0,DN
        End Do

        Return
    End Subroutine Init2

    Subroutine Init3(iconf)
        Implicit None

        Integer :: iconf, n, k, nj, nii, ni, i, njj, ih, lrec, n12, l
        Real(dp) :: dl, da, dz0, dz, qi, qj, qt, zef, rm1, a0, e, d
        Real(dp), Dimension(30) :: FAC
        Character(Len=256) :: strfmt

        FAC=0_dp

        Do NI=1,NS
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).EQ.-1) GOTO 200
        End Do

        Return

200     IF (KL.NE.3) Then
            strfmt = '(/"  GASHPAR INITIAL APPROXIMATION")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
        Else
            strfmt = '(/2X,"INITIAL APPROXIMATION:P(I)=0;Q(I)=0")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
        End If

        strfmt = '(6X,"NL",2X,"JJ",4X,"QQ",10X,"EE",10X,"D(0)",8X,"NIT",2X,"M1",3X,"M2",3X,"M3")'
        WRITE( *,strfmt)
        WRITE(11,strfmt)

        FAC(1)=0.d0
        Do I=2,30
            D=I-1.d0
            FAC(I)=FAC(I-1)+DLOG(D)
        End Do

        D=0.d0
        Do NI=1,NS
            IF (NC(NI).GT.0.AND.NC(NI).NE.ICONF) Cycle
            D=D+QQ(NI)
        End Do

        D=D-1.0D0
        DL=Z**0.333333D0
        DA=DL*1.19D0
        DL=DL*0.2075D0
        DZ0=0.D0

        Do NII=1,NS
            NI=NII
            IF (NC(NI).GT.0.AND.NC(NI).NE.ICONF) Cycle
            N=NN(NI)
            L=LL(NI)
            K=KK(NI)
            QI=QQ(NI)
            QJ=0.D0
            NJ=0
            Do NJJ=1,NS
                IF (NC(NJJ).GT.0.AND.NC(NJJ).NE.ICONF) Cycle
                IF (NI.EQ.NJJ) Cycle
                IF (NN(NI).NE.NN(NJJ)) Cycle
                IF (LL(NI).NE.LL(NJJ)) Cycle
                NJ=NJJ
                QJ=QQ(NJ)
            End Do
            QT=QI+QJ
            DZ=DZ0+0.5D0*(QT-1.0D0)
            IF (NJ.LT.NI) DZ0=DZ0+QT
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).NE.-1) Cycle
            ZEF=Z-DZ
            IF (ZEF.LE.0.01D0) ZEF=0.01D0
            Do I=1,IP6
                P(I)=0.d0
                Q(I)=0.d0
                Y(I)=0.d0
                CP(I)=0.d0
                CQ(I)=0.d0
            End Do
            E=0.5d0*(ZEF/N)**2
            P(II+1)=E
            P(II+2)=1.d0
            RM1=(L+0.5d0)**2/(Z+DSQRT(Z*Z-2*E*(L+0.5d0)**2))
            M1=(AL*(RM1-R(1))+BT*DLOG(RM1/R(1)))/H+1.D0
            IH=2-KT
            I=1+2*IH
            IF (M1.LT.I) M1=I
            M1=(M1/2)*2+1
            P(II+15)=M1
            P(II+3)=II
            P(II+5)=0.D0
            Q(II+5)=0.D0
            A0=DEXP(DLOG(2.D0*ZEF/N)+(L+0.5D0)*DLOG(2.D0*Z/N)-FAC(L+L+2)+0.5D0*(FAC(N+L+1)-FAC(N-L)))/DSQRT(2.d0*N)
            IF (K.GT.0) Q(II+5)=-A0*(LL(NI)+K+1)/(2*CL)
            IF (K.LT.0) P(II+5)=A0
            IF (KL.NE.3) Y(II+1)=D*(DA+DL)
            IF (KL.EQ.3) Then
                Continue
            Else
                Do I=1,II
                    Y(I)=D*(1.0d0-DEXP(-DL*R(I))/(1.0d0+DA*R(I)))-Z
                End Do
            End If

            P(II+16)=1.D0
            P(II+17)=1.D0
            Q(II+15)=M2
            N12=NI+4
            CALL WRITEF(12,N12,P,Q,2)
        End Do

        ! SAVING OF THE INITIAL APPROXIMATION
        N12=NS+NS+5
        CALL WRITEF(12,N12,P,Q,2)
        CALL CLOSEF(12)

        ! FILE HFD.DAT
        LREC=IP6
        CALL OPENFD('HFD.DAT',12,LREC,0)

        Return
    End Subroutine Init3

    Subroutine Fed(ni)
        Implicit None

        Integer :: ni, ih, k, n, mm, im, j, i0, i1, i2, j1, j2, j3, j4, n1, n3, nt1, nt3, nj, i, m, l
        Real(dp) :: eps, h1, hh, gam, e, e0, st, t, t1, s, hk, t2, t3, tp, tq, qm2, rj, pi, qi, pj, qj, d, dw, dc, da, db, dt, dq, pm2, e1, e3, dt1, dt3, dpp
        Real(dp), Dimension(IP6) :: W

        W=0_dp
  300   EPS=EPS0*0.1D0
        IF (KT.EQ.0) EPS=EPS1*0.1D0
        IH=2-KT
        H1=H*IH
        HH=0.5D0*H1
        K=KK(NI)
        GAM=P(II+4)
        L=LL(NI)
        N=NN(NI)
        E=P(II+1)
        E0=E
        M1=P(II+15)+0.01D0
        MM=II
        IMAX=II
        
        ! CALCULATION OF THE INHOMOGENEOUS PART
        ST=DSQRT(1.D0/P(II+2))
        IF (KF.EQ.1) ST=0.D0
        T=HH/CL
        T1=T*ST
        S=E
        HK=HH*K
        Do I=1,II,IH
            C(I)=HK
            A(I)=T*(S*R(I)+Y(I))
            CP(I)= CP(I)*T1
            CQ(I)=-CQ(I)*T1
            IF (R2.LT.0.D0) Cycle
            D=V(I)/R(I)
            C(I)=C(I)*D
            A(I)=A(I)*D
            CP(I)=CP(I)*D
            CQ(I)=CQ(I)*D
        End Do
        
        IM=II-IH
        IF (KT.EQ.0) GOTO 320

        ! CALCULATION OF THE POINT M3 (KT=1)
        IM=II+1
        M3=II
        S=1.D0/(H1*H1)
        Do I=1,II,IH
            J=IM-I
            IF (0.5d0*E*V(J)*V(J).GE.S) Cycle
            M3=J
            GOTO 310
        End Do

  310   MM=P(II+3)+0.01D0
        T3=DSQRT(0.5d0*E)/(2*CL)
        IM=MM-3*IH
        T=CQ(IM)+T3*CP(IM)
        T1=0.D0
        IF (T.NE.0.D0) T1=DLOG(DABS(T)/V(IM))
        I2=IM
        IF (M3.GT.IM) GOTO 320
        J1=IM+M3
        T2=T1
        Do I1=M3,IM,IH
            T1=T2
            I=J1-I1
            J=I-IH
            T=CQ(J)+T3*CP(J)
            T2=0.D0
            IF (T.NE.0.D0) T2=DLOG(DABS(T)/V(J))
            D=0.D0
            IF (T2.NE.0.D0.AND.T1.NE.0.D0) D=(T1-T2)/(R(I)-R(J))
            I2=I
            IF (D.EQ.0.D0) Cycle
            IF (E*R(I)*R(I).GT.500.D0) Cycle
            D=D*D
            IF (D*V(I)*V(I).LT.S) GOTO 320
        End Do

  320   IM=MM-3*IH
        Do I=1,IM,IH
            J=I+IH
            CP(I)=CP(J)+CP(I)
            CQ(I)=CQ(J)+CQ(I)
        End Do

        IF (KT.EQ.0) GOTO 330
        MM=I2+3*IH
        IM=MM-2*IH
        Do I=IM,IMAX,IH
            J=I-IH
            T=CQ(I)+T3*CP(I)
            T2=0.D0
            IF (T.NE.0.D0) T2=DLOG(DABS(T)/V(I))
            D=0.D0
            IF (T1.NE.0.D0.AND.T2.NE.0.D0) D=(T2-T1)/(R(I)-R(J))
            T1=T2
            CP(I)=-D*V(I)*HH
            CQ(I)=-T
        End Do

  330   IM=MM-3*IH
        S=ST/3.D0
        T=ST/8.D0
        D=ST/120.D0
        I0=1+IH+IH
        Do I=I0,IM,IH
            J =I+IH
            J1=J+IH
            J2=J1+IH
            I1=I-IH
            I2=I1-IH
            TP=S*(P(J)-P(I))-T*(P(J1)-P(I1))+D*(P(J2)-P(I2))
            TQ=S*(Q(J)-Q(I))-T*(Q(J1)-Q(I1))+D*(Q(J2)-Q(I2))
            CP(I)=CP(I)+TP
            CQ(I)=CQ(I)+TQ
        End Do

        I=1+IH
        J =I+IH
        J1=J+IH
        J2=J1+IH
        J3=J2+IH
        I1=I-IH
        TP=D*(9*P(I1)-25*P(I)+20*P(J)-5*P(J2)+P(J3))
        TQ=D*(9*Q(I1)-25*Q(I)+20*Q(J)-5*Q(J2)+Q(J3))
        CP(I)=CP(I)+TP
        CQ(I)=CQ(I)+TQ

        I=1
        J =I+IH
        J1=J+IH
        J2=J1+IH
        J3=J2+IH
        J4=J3+IH
        I1=I-IH
        TP=D*(29*P(I)-115*P(J)+180*P(J1)-140*P(J2)+55*P(J3)-9*P(J4))
        TQ=D*(29*Q(I)-115*Q(J)+180*Q(J1)-140*Q(J2)+55*Q(J3)-9*Q(J4))
        CP(I)=CP(I)+TP
        CQ(I)=CQ(I)+TQ

        N1=0
        N3=0
        NT1=0
        NIT=0
        NT3=0

        ! THE BEGINING OF THE ITERATION PROCEDURE
 600    IF (NIT.LE.40) GOTO 500
        strfmt = '(/"  NIT=",I3," E=",E15.6,"  Q(M2)=",E15.6," QM2=",E15.6/"  M1=",I5,"  M2=",I5,"  M3=",I5,"  NJ=",I3," N1=",I3)'
        WRITE( *,strfmt) NIT,E0,Q(M2),QM2,M1,M2,M3,NJ,N1
        WRITE(11,strfmt) NIT,E0,Q(M2),QM2,M1,M2,M3,NJ,N1
        IF (NIT.LE.50) GOTO 500
        strfmt = '(/"  NIT>50","  KF=",I1)'
        WRITE( *,strfmt) KF
        WRITE(11,strfmt) KF

        IF (KF.NE.0) CALL CLOSEF(11)
        IF (KF.NE.0) STOP
        KF=1
        
        GOTO 300
  500   NIT=NIT+1
        D=(E-E0)*HH/CL
        S=CL*H1
        ST=0.D0
        Do I=1,II,IH
            T=A(I)+D*V(I)
            A(I)=T
            B(I)=S*V(I)-A(I)
            W(I)=C(I)*C(I)+A(I)*B(I)
            IF (ST.GT.T) ST=T
        End Do
        IF (ST.LT.0.D0) GOTO 340
        N1=N1+1
        E0=E
        E=E*(1.D0-0.5D0*N1/(N1+1.D0))
        GOTO 600
  340   N1=0
        IF (KT.EQ.1) GOTO 350

        ! CALCULATION OF THE POINT M3 (KT=0)
        IM=II+1
        M3=II
        S=0.8D0
        Do I=1,II,IH
            J=IM-I
            IF (W(J).GE.S) Cycle
            M3=J
            GOTO 350
        End Do

  350   IM=M1+M3
        S=0.5D0*K*(K+1)*HH/CL
        Do I=M1,M3,IH
            J=IM-I
            RJ=R(J)
            IF (A(J)+S*V(J)/(RJ*RJ).GT.0.D0.AND.W(J).GT.0.D0) Cycle
            M2=J
            GOTO 360
        End Do

        E0=E
        E=E*0.8D0
        GOTO 600
        
  360   IF (M2.LT.M3-2*IH) GOTO 370
        E0=E
        E=E*1.2D0
        GOTO 600

  370   I0=1
        PI=0.D0
        QI=0.D0
        Do M=0,NMAX
            J=II+5+M
            PI=PI+P(J)
            QI=QI+Q(J)
        End Do
        ST=R(I0)**GAM
        PI=PI*ST
        QI=QI*ST
        P(I0)=(1.D0+C(I0))*PI+B(I0)*QI
        Q(I0)=(1.D0-C(I0))*QI+A(I0)*PI

        ! create functions P(I),Q(I) (I=1,M2)
        NJ=N-L-1
        I1=I0+IH
        PJ=P(I0)
        QJ=Q(I0)

        Do I=I1,M2,IH
            J=I-IH
            DW=W(J)
            DC=C(J)
            DA=A(J)
            DB=B(J)
            DT=0.5D0*(1.D0-DW)
            DQ=((DW+DC)*QJ-DA*PJ)/DT+CQ(J)
            DPP=((DW-DC)*PJ-DB*QJ)/DT+CP(J)
            J1=1
            IF (PJ.LE.0.D0) J1=-1
            PJ=PJ+DPP
            QJ=QJ+DQ
            P(I)=PJ
            Q(I)=QJ
            J2=1
            IF (PJ.LE.0.D0) J2=-1
            IF (J1.NE.J2) NJ=NJ-1
        End Do
        PM2=PJ
        QM2=QJ

        IF (NJ.EQ.0) GOTO 380
        I=1
        IF (NJ.LT.0) I=-1
        E0=E
        E=E*(1.D0-0.3D0*I)
        GOTO 600

  380   IM=MM-2*IH
        I1=IM+IMAX
        Do I=IM,IMAX,IH
            J=I1-I
            DW=W(J)
            DW=DSQRT(DW)
            DB=B(J)
            DC=C(J)
            PJ=(DW-DC)/DB
            QJ=CQ(J)/(DW+CP(J))
            P(J)=PJ
            Q(J)=QJ
        End Do

        ! wy~islenie Q(I)/P(I); (I=M2,IM)
        QJ=QJ*(1.D0-DC-PJ*DB)
        P(IM)=PJ
        Q(IM)=QJ
        J=IM-IH
  390   DC=C(J)
        DA=A(J)
        DB=B(J)
        DW=W(J)
        DQ=DC+DC
        DT=1.D0+DW+DQ+2*DB*PJ
        IF (J.LE.M3) DPP=(DA+DA+(1.D0+DW-DQ)*PJ)/DT
        IF (J.GT.M3) DPP=(DSQRT(DW)-DC)/DB
        QJ=(QJ+PJ*CP(J)-CQ(J))/DT*(1.D0-DW)
        PJ=DPP
        IF (J.EQ.M2) GOTO 400
        P(J)=PJ
        Q(J)=QJ
        J=J-IH
        GOTO 390

        !s{iwanie w to~ke M2
  400   DQ=QM2
        QM2=PJ*PM2+QJ
        DPP=QM2-DQ
        DQ=DQ+K/R(M2)*PM2/(2*CL)
        T=DABS(DPP)/(DABS(DQ)+EPS/CL)
        IF (T.LT.5*EPS.AND.NIT.GT.2) GOTO 440
        IF (DPP*B(M2).LE.0.D0) GOTO 410
        J=1
        E3=E
        DT3=DPP
        NT3=1
        GOTO 420
  410   E1=E
        DT1=DPP
        NT1=1
        J=-1
  420   IF (NT1.EQ.NT3) GOTO 430
        I=N-L-1
        N3=N3+1
        D=I-2*(I/2)-0.5D0
        T=0.01D0
        IF (KT.EQ.0) T=0.1D0
        IF (KT.EQ.0.OR.NIT.GT.10) T=0.3D0
        E0=E
        E=E*(1.D0+T*(J*N3)*D/(N3+3.D0))
        GOTO 600
  430   N3=0
        E0=E
        E=(DT3*E1-DT1*E3)/(DT3-DT1)
        GOTO 600

        ! wy~islenie funkcij P(I),Q(I); (I=M2,II)
  440   I1=M2+IH
        I2=M2
        IM=MM-2*IH
        Do I=I1,IM,IH
            J=I-IH
            D=C(J)+C(J)
            S=B(J)+B(J)
            T=((1.D0+W(J)-D)*P(J)-S*Q(J))/(1.D0-W(J))
            PJ=T+CP(J)
            QJ=P(I)*PJ+Q(I)
            IF (DABS(PJ).LT.1.D-12) GOTO 450
            I2=I
            P(I)=PJ
            Q(I)=QJ
        End Do

        ! postroenie nowyh funkcij P(I),Q(I)
  450   Do J=1,I2,IH
            PJ=P(J)
            QJ=Q(J)
            T=1.D0-W(J)
            P(J)=((1.D0-C(J))*PJ-B(J)*QJ)/T
            Q(J)=((1.D0+C(J))*QJ-A(J)*PJ)/T
        End Do

        IF (I2.LT.IM) GOTO 460

        I1=IM+IH
        I2=IM
        Do I=I1,II,IH
            J=I-IH
            D=(R(I)-R(J))/(HH*V(I))
            S=DSQRT(W(I))
            T=CP(I)-S
            T1=T*D
            IF (DABS(T1).LT.0.01D0) T2=-D
            IF (DABS(T1).GE.0.01D0) T2=(1.D0-DEXP(T1))/T
            T1=P(J)*DEXP(-S*D)
            T2=T2*B(I)*Q(I)
            PJ=T1+T2
            IF (DABS(PJ).LT.1.D-12) GOTO 460
            I2=I
            T2=T1
            Q(I)=PJ*P(I)+Q(I)
            P(I)=PJ
        End Do
        GOTO 1000

  460   IMAX=I2
        IF (IMAX.EQ.II) GOTO 1000
        IM=IMAX+IH
        Do I=IM,II,IH
            P(I)=0.D0
            Q(I)=0.D0
        End Do

  1000  P(II+1)=E
        P(II+3)=IMAX
        Q(II+15)=M2
        Q(II+3)=IMAX 

        Return
    End Subroutine Fed

    Subroutine NORM
        Implicit None
        
        Integer :: ih, i0, imax, n, m, i, j
        Real(dp) :: r0, v0, gam, g, t0, p0, d, f0, dt, f, h1

        IH=2-KT
        H1=H*IH
        I0=1
        IMAX=P(II+3)+0.01D0
        R0=R(I0)
        V0=V(I0)

        GAM=P(II+4)
        G=GAM+GAM
        T0=0.D0
        P0=0.D0
        Do M=0,NMAX
            I=II+5+M
            D=0.D0
            Do N=0,M
                J=II+5+N
                D=D+P(J)*P(I-N)+Q(J)*Q(I-N)
            End Do
            T0=T0+D/(G+M+1)
            P0=P0+D*(G+M)
        End Do
        T0=T0*R0**(G+1)
        P0=P0*R0**(G-1)
        F0=(P(I0)*P(I0)+Q(I0)*Q(I0))*V0
        P0=P0*V0*V0+F0*BT/(AL*R0+BT)**2

        DT=0.D0
        Do I=I0,IMAX,IH
            F=(P(I)*P(I)+Q(I)*Q(I))*V(I)
            DT=DT+F
        End Do
        D=T0+H1*DT-H1*(0.5D0*F0-H1/12.D0*P0)-0.5D0*H1*F
        Q(II+2)=D
        D=1.D0/DSQRT(D)
        Do I=1,IMAX,IH
            P(I)=P(I)*D
            Q(I)=Q(I)*D
        End Do

        P(II+2)=1.D0
        Do M=0,NMAX
            I=II+5+M
            P(I)=P(I)*D
            Q(I)=Q(I)*D
        End Do

        Return
    End Subroutine Norm

    Subroutine Iter(KI,IFIN,ICONF,DEL)
        Implicit None

        Integer :: nmax, ifin, ni, ki, iconf, nt, kpn, n12
        Real(dp) :: eps, del
        Character(Len=256) :: strfmt

        NMAX=30
        IFIN=0
        KPN=0

        EPS=EPS0
        IF (KT.EQ.0) EPS=EPS1

        Do NI=1,NS
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).EQ.1) Cycle
            KPN=0
        End Do

        IF (KPN.EQ.1) IFIN=1

        LAG=KLAG
        IF (NITER.EQ.1) LAG=0
        IF (NITER.GT.NMAX) Then
            strfmt = '(2X,"DIVERGENCE OF SELF-CONSISTENCY")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
            IFIN=1
            Return
        End If

        IF (NITER.GT.1) GOTO 200
        Do NI=1,NS
            NF(NI)=0
            IF (KL.EQ.1.AND.NI.LT.NI0) NF(NI)=1
        End Do
        KI=2
        IF (KL.EQ.0.OR.KL.EQ.3) KI=3
        IF (KL.EQ.1.OR.KL.EQ.2) LAG=KLAG
        GOTO 210

  200   KI=IABS(1-KI)
        IF (DEL.GE.EPS) GOTO 210
        If (KT.EQ.1) Then
            strfmt = '(2X,76("-")/2X,"END OF SELF-CONSISTENCY")'
            WRITE( *,strfmt)
            strfmt = '(2X,67("-")/2X,"END OF SELF-CONSISTENCY")'
            WRITE(11,strfmt)
            IFIN=1
            RETURN
        End If
        KT=1
        CALL PQINT(ICONF)
        KI=2
        CALL READF (12,1,P,Q,2)
        P(9)=KT
        CALL WRITEF(12,1,P,Q,2)
        Do NI=1,NS
            NF(NI)=0
        End Do

  210   NIT=0
        M2=0
        M3=0
        NT=II/(2-KT)
        NT=(NT/2)*2+1
                       
        If (NITER.EQ.1) Then
            strfmt = '(/2X,76("="))' 
            WRITE( *,strfmt)
            strfmt = '(/2X,67("="))'
            WRITE(11,strfmt)
        End If

        If (NITER.GT.1) Then
            strfmt = '(2X,76("="))'
            WRITE( *,strfmt)
            strfmt = '(2X,67("="))'
            WRITE(11,strfmt)
        End If

        IF (KI.EQ.0) Then
            strfmt = '(2X,"CONF=",I3,3X,"ITER = ",I3,3X,"NT = ",I3,3X,"AITKEN ITERATION")'
            WRITE( *,strfmt) ICONF,NITER,NT
            WRITE(11,strfmt) ICONF,NITER,NT
            Do NI=1,NS
                IF (NC(NI).NE.ICONF) Cycle
                IF (KP(NI).EQ.1) Cycle
                IF (NF(NI).EQ.1) Cycle
                N12=NI+NS+4
                CALL READF (12,N12,P,Q,2)
                N12=NI+4
                CALL WRITEF(12,N12,P,Q,2)
            End Do
        Else
            strfmt = '(2X,"CONF=",I3,3X,"ITER = ",I3,3X,"NT = ",I3,3X,"DIRECT ITERATION")'
            WRITE( *,strfmt) ICONF,NITER,NT
            WRITE(11,strfmt) ICONF,NITER,NT
        End If        

        strfmt = '(2X,67("-")/5X,"JJ",4X,"QQ",9X,"EE",8X,"dE",8X,"D(0)",9X,"DE",6X,"DN")'
        Write(11,strfmt)
        DEL=0.D0

        Return
    End Subroutine Iter

    Subroutine Dens(NI,ICONF)
        Implicit None

        Integer :: ih, i0, i, j, n, ni, nj, id, ig, ig0, ik, iconf, m, n12
        Real(dp) :: dr, r0d, s, dh, t, th, h1, d
        Real(dp), Dimension(IP6) :: RO

        RO=0_dp
        IH=2-KT
        H1=H*IH
        I0=1
        ! CALCULATION OF THE TOTAL DENSITY (NI=0)
        If (NI.EQ.0) Then
            Do I=1,IP6
                RO(I)=0.D0
            End Do
            RO(II+4)=G0
        Else
            CALL READF(12,4,RO,RO,1)
            T=QQ(NI)/P(II+2)
            TH=T*H1
        End If

        Do NJ=1,NS
            IF (NI.GT.0.AND.NI.NE.NJ) Cycle
            IF (NI.EQ.0.AND.NC(NJ).GT.0.AND.NC(NJ).NE.ICONF) Cycle
            N12=NJ+4
            CALL READF(12,N12,A,B,2)
            D=QQ(NJ)/A(II+2)
            DH=D*H1
            IMAX=A(II+3)+0.01D0
            Do I=I0,IMAX,IH
                S=DH*(A(I)**2+B(I)**2)
                IF (NI.GT.0) S=TH*(P(I)**2+Q(I)**2)-S
                RO(I)=RO(I)+S*V(I)
            End Do
            IK=IABS(KK(NJ))
            IG=IK+IK
            IG0=G0+0.5D0
            ID=IG-IG0
            R0D=R(1)**ID
            Do M=0,NMAX
                I=II+5+M
                DR=0.D0
                Do N=0,M
                    J=II+5+N
                    DR=DR+D*(A(J)*A(I-N)+B(J)*B(I-N))
                    IF (NI.GT.0) DR=-DR+T*(P(J)*P(I-N)+Q(J)*Q(I-N))
                End Do
                IF (M+ID.LE.NMAX) RO(I+ID)=RO(I+ID)+DR*R0D
            End Do
        End Do
        CALL WRITEF(12,4,RO,RO,1)
        Return
    End Subroutine Dens

    Subroutine Y0(NI)             ! changed 30/11/04
        ! CALCULATION OF THE COULOMB POTENTIAL FOR K=0 (Y(I),I=1,II)
        Implicit None

        Integer :: ih, i0, i1, j, n, id, ig0, ig, ik, ni, nj, ip, im, i, imax, m
        Real(dp) :: r0, v0, dt, fm, p0, t0, g, g0, dr, r0d, dh, ww, di, dk, s, f, f0, t, const, h1, d
        Real(dp), Dimension(IP6) :: RO, C

        RO=0_dp
        C=0_dp

        IH=2-KT
        H1=H*IH
        CONST=0.1D0
        I0=1
        R0=R(I0)
        V0=V(I0)
        DK=1.D0-DEXP(-H1)
        Iconf=0
        IF (ni.gt.0) ICONF=NC(NI)
        DI=ZAT(ICONF+1)

        CALL READF(12,4,RO,RO,1)

        If (NI.EQ.0) Then
            Continue
        Else
            WW=(QW(NI)-1)/(4*LL(NI)+1)
            Do NJ=1,NS
                IF (NC(NJ).NE.ICONF) Cycle
                IF (NN(NI).NE.NN(NJ)) Cycle
                IF (LL(NI).NE.LL(NJ)) Cycle
                IF (JJ(NI).EQ.JJ(NJ).AND.NI.NE.NJ) Cycle
                IF (NJ.NE.NI) Then
                    CALL READF (12,NJ+4,A,B,2)
                    D=(JJ(NJ)+1)*WW/A(II+2)
                Else
                    D=(JJ(NI)*WW-QW(NI))/P(II+2)
                    Do I=1,IP6
                        A(I)=P(I)
                        B(I)=Q(I)
                    End Do
                End If
                DH=D*H1
                IMAX=A(II+3)+CONST
                Do I=I0,IMAX,IH
                    RO(I)=RO(I)+DH*(A(I)**2+B(I)**2)*V(I)
                End Do
                IK=IABS(KK(NJ))
                IG=IK+IK
                G0=RO(II+4)
                IG0=G0+0.5D0
                ID=IG-IG0
                R0D=R(1)**ID
                Do M=0,NMAX-ID   ! Taylor expansion for electron density
                    I=II+5+M
                    DR=0.D0
                    Do N=0,M
                        J=II+5+N
                        DR=DR+D*(A(J)*A(I-N)+B(J)*B(I-N))
                End Do
                RO(I+ID)=RO(I+ID)+DR*R0D
                End Do
            End Do
        End If

        G=RO(II+4)
        T0=0.D0
        P0=0.D0
        Do M=0,NMAX
          I=II+5+M
          T0=T0+RO(I)/(G+M+1)
          P0=P0+RO(I)*(G+M)
        End Do
        T0=T0*R0**(G+1)
        P0=P0*R0**(G-1)
        F0=RO(I0)
        P0=H1*P0*V0*V0+F0*BT/(AL*R0+BT)**2

        ! CALCULATION OF THE IMAX
        IMAX=II+IH
  240   IMAX=IMAX-IH
        IF (RO(IMAX).EQ.0.D0.AND.IMAX.GT.I0+IH) GOTO 240
        FM=RO(IMAX)

        ! CALCULATION OF THE FUNCTION Z(0;R)
        DT=T0+0.5D0*F0+H1/12.D0*P0
        C(I0)=DT
        I1=I0+IH
        Do I=I1,IMAX,IH
            F=RO(I)
            DT=DT+F
            C(I)=DT
        End Do

        ! CALCULATION OF THE FUNCTION Y(0;R)
        T=-0.5D0*FM
        DT=T
        DH=H1/12.D0
        D=DK
        IM=IMAX-IH
        I1=IM+I0
        Do J=I0,IM,IH
            I=I1-J
            IP=I+IH
            S=DH
            If (R2.LT.0.D0) Then
                Continue
            Else
                D=R(I)/R(IP)
                D=1.D0-D
                S=DH*V(I)/R(I)
            End If
            F=RO(IP)
            F=F-D*(T+F)
            DT=DT+F
            T=DT
            C(I)=C(I)+DT-(S*RO(I))
        End Do
        DT=C(IMAX)-0.5D0*FM
        Do I=IMAX,II,IH
            C(I)=DT
        End Do

        ! CALCULATION OF THE Y(I)-Z AND CORRECTION
        DT=-1.D0-DI
        IF (NI.EQ.0) DT=-DI
        DT=DT-C(II)
        IF (NITER.EQ.1.AND.KL.EQ.3) DT=-Z
        I1=II+1
        Do I=1,II,IH
            J=I1-I
            Y(J)=C(J)+DT
        End Do

        G0=RO(II+4)
        ID=G0+0.5D0
        R0D=R0**ID
        T0=0.D0
        Do M=0,NMAX
            I=II+5+M
            T0=T0+(RO(I)/(G0+M)-RO(I)/(G0+M+1))
        End Do
        T0=T0*R0**G0
        DO M=0,NMAX
            I=II+1+M
            Y(I)=0.D0
        End Do
        Y(II+1)=C(I0)/R0+T0
        Do M=0,NMAX
            I=II+5+M
            IF ((M+ID).LE.NMAX) Y(II+1+M+ID)=RO(I)*(1.D0/(G0+M+1)-1.D0/(G0+M))*R0D
        End Do

        IF (NI.EQ.0) CALL WRITEF (12,4,RO,Y,2)

        Return
    End Subroutine Y0

    Subroutine Plex(NI)
        Implicit None
        
        Integer :: ni, ih, na, ja, ka, nj, nb, jb, kb, kmin, kmax, k1, k, id, n, i, j, nk, ilag, kc, n12, nnc, m
        Real(dp) :: h1, c1, c2, qa, qb, gab, s, eab, dr, r0d, dpp, dq, dab, qc, gac, dac, d
        Real(dp), Dimension(10) :: ELAG
        
        iwr=0

        ! CALCULATION OF THE COULOMB (Y(I),I=1,II) AND EXCHANGE (U(I),I=1,II) POTENTIALS
        IH=2-KT
        H1=H*IH
        C1=0.01D0

        IF (NG.EQ.0) Then
            Continue
        Else
            CALL READF(12,3,ROIJ,ROIJ,1)
            Do I=1,NG
                C2=DSIGN(C1,ROIJ(2*I-1))
                IGAM(I)=ROIJ(2*I-1)+C2
                GAM (I)=ROIJ(2*I)
            End Do
        End If
        
        NA=NI
        JA=JJ(NA)
        QA=QQ(NA)
        KA=KK(NA)

        Do I=1,10
            ELAG(I)=0.D0
        End Do
        Do I=1,IP6
            CP(I)=0.D0
            CQ(I)=0.D0
        End Do

        Do NJ=1,NS
            NB=NJ
            IF (NC(NB).GT.0.AND.NC(NB).NE.NC(NA)) Cycle  !### configuration
            JB=JJ(NB)
            QB=QQ(NB)
            KB=KK(NB)
            if (NA.EQ.NB) then
                Do I=1,IP6
                    A(I)=P(I)
                    B(I)=Q(I)
                End Do
            else
                N12=NB+4
                CALL READF(12,N12,A,B,2)
            end if

            ! COULOMB (K>0)
            KMIN=3
            kmax=min(ja,jb)+1
            If (KMAX.LT.KMIN) Then            !### skipped for j=1/2 waves
                Continue
            Else
                Do K1=KMIN,KMAX                    !###
                    K=K1-1                              !### K is multipolarity
                    GAB=COEF(-1,K,NA,NB)
                    IF (DABS(GAB).LT.1.D-7) Cycle
                    IMAX=A(II+3)+C1                     !### Last used grid node
                    Do I=1,IMAX,IH
                        ROIJ(I)=H1*V(I)*(A(I)**2+B(I)**2)    !### density
                    End Do
                    ROIJ(II+3)=IMAX
                    ROIJ(II+4)=2*A(II+4)                  !### gamma for density
                    Do M=0,NMAX                      !### expansion at the origin
                        I=II+5+M
                        DR=0.D0
                        Do N=0,M
                            J=II+5+N
                            DR=DR+(A(J)*A(I-N)+B(J)*B(I-N))
                        End Do
                        ROIJ(I)=DR
                    End Do

                    CALL YK(K)                          !### Yk: uses ROIJ; forms C
                    D=GAB/A(II+2)                       !### A(ii+2)=norma
                    IF (NA.EQ.NB) D=D*2
                    Do I=1,II,IH
                        Y(I)=Y(I)+D*C(I)                  !### = Coulomb potential
                    End Do                           !#### Note that symmetric
                    Do M=0,NMAX                     !##### part is calculated by
                        I=II+1+M                          !###### Y0 separately
                        Y(I)=Y(I)+D*C(I)
                    End Do
                End Do
            End If
        
            ! EXCHANGE
            IF (NA.EQ.NB) Cycle
            imax=dmax1(p(ii+3),a(ii+3))+c1
            Do I=1,IMAX,IH
                ROIJ(I)=H1*V(I)*(P(I)*A(I)+Q(I)*B(I))  !### exchange density
            End Do
            ROIJ(II+3)=IMAX
            ROIJ(II+4)=P(II+4)+A(II+4)              !### gamma
            Do M=0,NMAX                       !### expansion at the origin
                I=II+5+M
                DR=0.D0
                Do N=0,M
                  J=II+5+N
                  DR=DR+P(J)*A(I-N)+Q(J)*B(I-N)
                End Do
                ROIJ(I)=DR
            End Do
            KMIN=IABS(JA-JB)/2+1
            KMAX=(JA+JB)/2+1
            EAB=0.D0
            
            Do K1=KMIN,KMAX
                K=K1-1
                I=LL(NA)+LL(NB)+K
                IF (I.NE.2*(I/2)) Cycle
                GAB=COEF(1,K,NA,NB)
                D=GAB/A(II+2)                       !### A(ii+2)=norma
                CALL YK(K)                          !### Yk: uses ROIJ, forms C
                Do I=1,II,IH
                    S=D*C(I)
                    CP(I)=CP(I)+S*B(I)
                    CQ(I)=CQ(I)+S*A(I)
                End Do
                ID=IABS(KB)-IABS(KA)
                R0D=R(1)**ID
                Do M=0,NMAX
                    I=II+1+M
                    DPP=0.D0
                    DQ=0.D0
                    Do N=0,M
                        J=II+5+N
                        DPP=DPP+B(J)*C(I-N)
                        DQ=DQ+A(J)*C(I-N)
                    End Do
                    IF (M+ID.LT.0) Cycle
                    IF (M+ID.GT.NMAX) Cycle
                    CP(I+ID)=CP(I+ID)+D*DPP*R0D
                    CQ(I+ID)=CQ(I+ID)+D*DQ*R0D
                End Do

                IF (LAG.EQ.0) Cycle
                ! CALCULATION OF THE LAGRANGE MULTIPLIER
                IF (KA.NE.KB) GOTO 360
                IF (NN(NA).EQ.NN(NB)) GOTO 360
                IF (DABS(QA-QB).LT.1.D-3) GOTO 360
                DAB=2*QA*QB/(QA-QB)
                GAB=2*COEF(-1,K,NA,NA)-COEF(-1,K,NB,NA)-COEF(1,K,NB,NA)
                IF (DABS(GAB).LT.1.D-7) GOTO 340
                Do I=1,II,IH
                    W(I)=C(I)*(P(I)**2+Q(I)**2)
                End Do
                W(II+4)=2*P(II+4)+K
                CALL SINT(S)                         !### Sint integrates W
                EAB=EAB+DAB*GAB*S/P(II+2)
      340       GAB=2*COEF(-1,K,NB,NB)-COEF(-1,K,NA,NB)-COEF(1,K,NA,NB)
                IF (DABS(GAB).LT.1.D-7) GOTO 350
                Do I=1,II,IH
                    W(I)=C(I)*(A(I)**2+B(I)**2)
                End Do
                W(II+4)=2*A(II+4)+K
                CALL SINT(S)
                EAB=EAB-DAB*GAB*S/A(II+2)
      350       Do NK=1,NS
                    NNC=NK
                    IF (NC(NNC).NE.NC(NA)) Cycle
                    QC=QQ(NNC)
                    IF (NNC.EQ.NA.OR.NNC.EQ.NB) Cycle
                    GAB=COEF(-1,K,NA,NNC)-COEF(-1,K,NB,NNC)
                    IF (DABS(GAB).LT.1.D-7) Cycle
                    N12=NNC+4
                    CALL READF (12,N12,W,V,2)
                    Do I=1,II,IH
                        W(I)=C(I)*(W(I)**2+V(I)**2)
                    End Do
                    CALL READF (12,2,R,V,2)
                    W(II+4)=2*W(II+4)+K
                    D=W(II+2)
                    CALL SINT(S)
                    EAB=EAB+DAB*GAB*S/D
                End Do

      360       ILAG=0
                Do NK=1,NS
                    NNC=NK
                    IF (NC(NNC).NE.NC(NA)) Cycle
                    KC=KK(NNC)
                    QC=QQ(NNC)
                    IF (NA.EQ.NNC) Cycle
                    IF (KA.NE.KC) Cycle
                    IF (NN(NNC).EQ.NN(NA)) Cycle
                    ILAG=ILAG+1
                    IF (NNC.EQ.NB) Cycle
                    IF (DABS(QA-QC).LT.1.D-3) Cycle
                    DAC=2*QA*QC/(QA-QC)
                    GAC=COEF(1,K,NA,NB)-COEF(1,K,NNC,NB)
                    IF (DABS(GAC).LT.1.D-7) Cycle
                    N12=NNC+4
                    CALL READF (12,N12,W,V,2)
                    Do I=1,II,IH
                        W(I)=C(I)*(W(I)*A(I)+V(I)*B(I))
                    End Do
                    CALL READF (12,2,R,V,2)
                    W(II+4)=W(II+4)+A(II+4)+K
                    D=DSQRT(W(II+2))*A(II+2)
                    CALL SINT(S)
                    ELAG(ILAG)=ELAG(ILAG)+DAC*GAC*S/D
                End Do
            End Do
            
            IF (LAG.EQ.0) Cycle
            IF (DABS(EAB).LT.1.D-7) Cycle
        
            iwr=iwr+1  !### warning: lagrange factors present!
        
            D=0.5D0*EAB/(A(II+2)*QA)
            Do I=1,II,IH
                CP(I)=CP(I)+D*B(I)*R(I)
                CQ(I)=CQ(I)+D*A(I)*R(I)
            End Do
            Do M=0,NMAX
                I=II+M+1
                CP(I)=CP(I)+D*B(I+4)
                CQ(I)=CQ(I)+D*A(I+4)
            End Do
        End Do

        IF (LAG.EQ.0) Return

        ILAG=0
        Do NJ=1,NS
            NB=NJ
            IF (NC(NB).NE.NC(NA)) Cycle
            KB=KK(NB)
            IF (NB.EQ.NA) Cycle
            IF (KB.NE.KA) Cycle
            IF (NN(NB).EQ.NN(NA)) Cycle
            ILAG=ILAG+1
            EAB=ELAG(ILAG)
            IF (DABS(EAB).LT.1.D-7) Cycle
            N12=NB+4
            CALL READF (12,N12,A,B,2)
            D=0.5D0*EAB/(DSQRT(A(II+2))*QA)
            Do I=1,II,IH
                CP(I)=CP(I)+D*B(I)*R(I)
                CQ(I)=CQ(I)+D*A(I)*R(I)
            End Do
            Do M=0,NMAX
                I=II+M+1
                CP(I)=CP(I)+D*B(I+4)
                CQ(I)=CQ(I)+D*A(I+4)
            End Do
        End Do

        Return
    End Subroutine Plex

    Subroutine Sint(s)
        Implicit None
        Integer :: ih, i0, i1, i2, i3, i
        Real(dp) :: h1, v0, r0, gam, r1, rr2, r3, t0, p0, f0, g, f1, f2, f3, f21, f32, f321, c0, c1, c2, dt, f, s, t

        IH=2-KT
        h1=H*IH
        I0=1
        V0=V(I0)
        R0=R(I0)
        GAM=W(II+4)

        I1=I0
        I2=I1+IH
        I3=I2+IH
        R1=R(I1)
        rr2=R(I2)
        R3=R(I3)
        T0=0.D0
        P0=0.D0
        F0=0.D0
        IF (GAM.GT.5.D0) GOTO 200
        G=GAM+1.D0
        F1=W(I1)/R1**G
        F2=W(I2)/rr2**G
        F3=W(I3)/R3**G
        F21=(F2-F1)/(rr2-R1)
        F32=(F3-F2)/(R3-rr2)
        F321=(F32-F21)/(R3-R1)
        C2=F321
        C1=F21-C2*(R1+rr2)
        C0=F1-C1*R1-C2*R1**2
        T0=R0**G*(C0/G+C1/(G+1.D0)*R0+C2/(G+2.D0)*R0**2)
        F0=W(I0)*V0/R0
        G=G-1.D0
        P0=R0**(G-1.D0)*(C0*G+C1*(G+1.D0)*R0+C2*(G+2.D0)*R0**2)
        P0=P0*V0*V0+F0*BT/(AL*R0+BT)**2

  200   DT=0.D0
        Do I=I0,II,IH
        F=W(I)
        IF (R2.GT.0.D0) F=F*V(I)/R(I)
        DT=DT+F
        End Do
        T=T0+H1*DT-0.5D0*H1*F0+H1*H1/12.D0*P0
        S=T

        Return
    End Subroutine

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

    Real(dp) Function Coef(I,K,NA,NB)
        Implicit None
        
        Integer :: j, ja, jb, n, n1, n2, la, lb, na, nb, i, k
        Real(dp) :: qa, g, ww, qb, wa, d, y
        Integer, Dimension(70) :: IG
        Real(dp), Dimension(70) :: GG

        ! GG=2*CH(J1,0.5,J2,-0.5,K,0)**2/(2*K+1)
        DATA GG/ &
            1*0.10000000000000D+01,2*0.33333333333333D+00, &
            2*0.20000000000000D+00,2*0.14285714285714D+00, &
            2*0.11111111111111D+00,1*0.90909090909091D-01, &
            1*0.50000000000000D+00,1*0.33333333333333D-01, &
            1*0.10000000000000D+00,1*0.12857142857143D+00, &
            1*0.20000000000000D+00,1*0.28571428571429D-01, &
            1*0.57142857142857D-01,1*0.95238095238095D-01, &
            1*0.12857142857143D+00,1*0.23809523809524D-01, &
            1*0.39682539682540D-01,1*0.75757575757576D-01, &
            1*0.95238095238095D-01,1*0.20202020202020D-01, &
            1*0.30303030303030D-01,1*0.62937062937063D-01, &
            1*0.33333333333333D+00,1*0.95238095238095D-02, &
            1*0.76190476190476D-01,1*0.25396825396825D-01, &
            1*0.31746031746032D-01,1*0.72150072150072D-01, &
            1*0.14285714285714D+00,1*0.95238095238095D-02, &
            1*0.47619047619048D-01,2*0.21645021645022D-01, &
            1*0.58275058275058D-01,1*0.95238095238095D-01, &
            1*0.86580086580087D-02,1*0.34632034632035D-01, &
            1*0.18648018648019D-01,1*0.16317016317016D-01, &
            1*0.48951048951049D-01,1*0.25000000000000D+00, &
            1*0.39682539682540D-02,1*0.59523809523810D-01, &
            1*0.97402597402597D-02,1*0.29220779220779D-01, &
            1*0.18731268731269D-01,1*0.14568764568765D-01, &
            1*0.47591297591298D-01,1*0.11111111111111D+00, &
            1*0.43290043290043D-02,1*0.38961038961039D-01, &
            1*0.89910089910090D-02,1*0.20979020979021D-01, &
            1*0.16317016317016D-01,1*0.10878010878011D-01, &
            1*0.40312628547923D-01,1*0.20000000000000D+00, &
            1*0.20202020202020D-02,1*0.48484848484848D-01, &
            1*0.47952047952048D-02,1*0.25174825174825D-01, &
            1*0.83916083916084D-02,1*0.14918414918415D-01, &
            1*0.14333379039261D-01,1*0.80625257095845D-02, &
            1*0.34371820130334D-01/

        DATA  IG/ &
            011,111,113,213,215,315,317,417,419,519, &
            033,133,233,333,135,235,335,435,237,337, &
            437,537,339,439,539,639, 55,155,255,355, &
            455,555,157,257,357,457,557,657,259,359, &
            459,559,659,759, 77,177,277,377,477,577, &
            677,777,179,279,379,479,579,679,779,879, &
            099,199,299,399,499,599,699,799,899,999/

        IMAX=70
        G=0.D0

        LA=LL(NA)
        JA=JJ(NA)
        QA=QQ(NA)
        WA=QW(NA)
        LB=LL(NB)
        JB=JJ(NB)
        QB=QQ(NB)
        WW=(WA-1)/(4*LA+1)

        IF (K.GT.0.OR.I.GT.0) GOTO 200
        IF (NA.EQ.NB) G=0.5D0*((QA-WA)+JA*WW)
        IF (NA.EQ.NB) GOTO 1000
        IF (NA.NE.NB) G=QB
        IF (NC(NA).EQ.NC(NB).AND.NN(NA).EQ.NN(NB).AND.LA.EQ.LB) G=(QB+(JB+1)*WW)
        GOTO 1000

  200   IF (NA.NE.NB.AND.I.LT.0) GOTO 220

        ! K+LA+LB - EVEN
        J=K+LA+LB
        IF (J.NE.2*(J/2)) GOTO 1000
        N2=K*100+JB*10+JA
        N1=K*100+JA*10+JB
        Y=0.D0
        Do J=1,IMAX
            N=IG(J)
            IF (N1.NE.N.AND.N2.NE.N) Cycle
            Y=GG(J)
            Exit
        End Do
        IF (NA.NE.NB) D=QB
        IF (NA.EQ.NB) D=0.5D0*((QA-WA)*(JA+1.D0)/JA+WW*(JA+1))
        IF (NA.NE.NB.AND.NC(NA).EQ.NC(NB).AND.NN(NA).EQ.NN(NB).AND.LA.EQ.LB)D=(QB+WW*(JB+1))
        G=-0.5D0*D*Y
  220   IF (NG.EQ.0) GOTO 1000
        IF (DABS(QA).LT.1.D-7) GOTO 1000
        N1=I*(4096*K+NA*64+NB)
        N2=I*(4096*K+NB*64+NA)
        Do J=1,NG
            N=IGAM(J)
            IF (N.NE.N1.AND.N.NE.N2) Cycle
            G=G+GAM(J)
            GOTO 230
        End Do
  230   G=G/QA
  1000  IF (DABS(G).LT.1.D-7) G=0.D0
        COEF=G
        Return
    End Function Coef

    Subroutine sms_core(na)  
        Use diff, Only : Dif
        Use wigner, Only : FJ3
        !### Exchange core potential due to sms
        !### see 12/11/99 for derivation details
        Implicit None

        Integer :: na, ih, ja, la, icore, nc, jc, lc, ip, ila, ilc, iq, lx, ly, i
        Real(dp) :: small, sn, err, f0, s1, cp1, cq1, e_sms

        Character :: let(5), ch1
        Character(Len=256) :: strfmt
        if (k_is.GE.3) then      ! normal mass shift is included for k_is=3,4
            call V_nms(na)         
        else
            eadd(na)=0.d0
        end if
        if (k_is.EQ.3) return    ! specific mass shift is included for k_is=2,4
        small=0.1d-8
        ih=2-kt

        sn=1.d0/P(ii+2)   !### normalization factor

        do i=1,IP6
            Pa(i)=P(i)
            Qa(i)=Q(i)
        end do
        call Dif(Pa,P1a,R,V,Pa(ii+4),ii,kt,MaxT,h)
        if (klow.GE.1) call Dif(Qa,Q1a,R,V,Pa(ii+4),ii,kt,MaxT,h)

        ja=Jj(na)
        la=Ll(na)
        xja=0.5d0*ja

        do i=1,IP6        !### (Py,Qy) = V^sms_core * (Pa,Qa)
            Py(i)=0.d0
            Qy(i)=0.d0
        end do
        Py(ii+4)=5.d0     !### this will be reduced later
        icore=0
        do nc=1,n_is
            jc=Jj(nc)
            lc=Ll(nc)
            ip=lc-la
            if (iabs(ja-jc).GT.2.OR.iabs(ip).NE.1) Cycle
            xjc=0.5d0*jc
            ila=ja-la   !
            ilc=jc-lc   ! these are for lower component
            iq=ilc-ila  !
            call ReadF(12,nc+4,Pc,Qc,2)
            err=dabs(1.d0-Pc(ii+2))
            if (err.GT.small) then
                strfmt = '(4X,"sms_core warning:",/4X,"norm of orbital ",I2," differs from unity by",E10.2)'
                write(*,strfmt) nc,err
                read(*,*)
            end if
            call Dif(Pc,P1c,R,V,Pc(ii+4),ii,kt,MaxT,h)
            if (klow.GE.1) call Dif(Qc,Q1c,R,V,Pc(ii+4),ii,kt,MaxT,h)
            f0=(jc+1)*FJ3(xja,xjc,1.d0,0.5d0,-0.5d0,0.d0)**2 ! we assume closed
            s1=+c_is*f0*P_eff(na,ja,la,nc,jc,lc)             !! shells for core
            if (dabs(s1).LT.small) Cycle
            icore=icore+1
            lx=ip*(la+lc+1)/2
            ly=iq*(ila+ilc+1)/2
            do i=1,ii,ih
                Py(i)=Py(i) + s1*(P1c(i)+lx*Pc(i)/R(i))
                if (klow.GE.1) Qy(i)=Qy(i) + s1*(Q1c(i)+ly*Qc(i)/R(i))
            end do

            if (Py(ii+4).GT.P1c(ii+4)) Py(ii+4)=P1c(ii+4)
        end do
        if (icore.EQ.0) return

        cp1=CP(1)
        cq1=CQ(1)
        do i=1,ii,ih
            CP(i)=CP(i)+Qy(i)*R(i)
            CQ(i)=CQ(i)+Py(i)*R(i)
            C(i)=(Pa(i)*Py(i)+Qa(i)*Qy(i))*sn
        end do
        C(ii+4)=Pa(ii+4)+Py(ii+4)
        call Sint1(e_sms)
        eadd(na)=eadd(na)+e_sms

        cp1=CP(1)/cp1
        cq1=CQ(1)/cq1
        do i=ii+1,ii+10
            CP(i)=CP(i)*cp1
            CQ(i)=CQ(i)*cq1
        end do

        ch1=let(la+1)
        strfmt = '(2X,"E_sms(",I2,A1,I1,") = ",E14.7)'
        if (kout.GT.1) write( *,strfmt) Nn(na),ch1,ja,e_sms
        if (kout.GT.0) write(11,strfmt) Nn(na),ch1,ja,e_sms
        Return
    End Subroutine SMS_Core

    Real(dp) Function P_eff (na,ja,la,nc,jc,lc ) 
        Use wigner, Only : FJ3
        !### effective vertex for (pi dot pk)
        !### see 19/09/03 for derivation details
        Implicit None
   
        Integer :: na, ja, la, nc, jc, lc, ila, ilc, lmax, ilmax, ih, i
        Real(dp) :: ga, half, xjc, xja, xjx, xjn, xla, xlc, yla, ylc, s0, ds, s1, s2, ss

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
   
        do i=1,Ii,ih
            C(i)=Pc(i)*(P1a(i)+(la-lc)*lmax/R(i)*Pa(i))
        end do
   
        if (klow.GT.0) then       !### lower component
            do i=1,Ii,ih
                C(i)=C(i)+Qc(i)*(Q1a(i)+(ila-ilc)*ilmax/R(i)*Qa(i))
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
            s0=dsqrt((xjx-xjn+1.d0)/(2*xja+1.d0))/FJ3(xjc,xja,1.d0,half,-half,0.d0)
            if (mod(la+1,2).NE.0) s0=-s0
            s1=0.d0
            if (lc.EQ.ila) s1=dsqrt((xjc+3*xja-2*ila+1)/(2*ila+1))
            ss=Z/DPcl*(s0*s1+1.0d0)

            do i=1,Ii,ih
                C(i)=C(i)+ss/R(i)*Pc(i)*Qa(i)
            end do
          
            s2=0.d0
            if (ilc.EQ.la) s2=dsqrt((xjc+3*xja-2*la+1)/(2*la+1))
            ss=Z/DPcl*(s0*s2-1.0d0)

            do i=1,Ii,ih
                C(i)=C(i)+ss/R(i)*Qc(i)*Pa(i)
            end do
        end if
   
        call Sint1(ds)
        P_eff=ds

        Return
    End Function P_eff

    Subroutine Breit_Core(na)  
        !### Exchange core potential due to Breit
        !### see 12/11/99 for derivation details
        !     Exchange core potential:
        !     <b|V_B|a> = sum_c sum_k Gaunt^k_a,c,c,b R^k_a,c,c,b
        Implicit None

        Integer :: na, ih, ja, la, nc, jc, i
        Real(dp) :: small, sn, err, cp1, cq1, e_breit
        Character :: ch1
        Character(Len=256) :: strfmt

        small=0.1d-8
        ih=2-kt

        sn=1.d0/P(ii+2)   !### normalization factor

        do i=1,IP6
            Pa(i)=P(i)
            Qa(i)=Q(i)
        end do
        ja=Jj(na)
        la=Ll(na)
        xja=0.5d0*ja

        do i=1,IP6        !### (Py,Qy) = V_core_Breit * (Pa,Qa)
            Py(i)=0.d0
            Qy(i)=0.d0
        end do
        Py(ii+4)=5.d0     !### this will be reduced later
        do nc=1,Nso
            call ReadF(12,nc+4,Pc,Qc,2)
            err=dabs(1.d0-Pc(ii+2))
            if (err.GT.small) then
                strfmt = '(4X,"Breit_Core warning:",/4X,"norm of orbital ",I2," differs from unity by",E10.2)'
                write(*,strfmt) nc,err
                read(*,*)
            end if
            jc=Jj(nc)
            xjc=0.5d0*jc

            ! Using subroutines of IIT to calculate Breit [magnetic and retardation]
            call breit_mag(na,nc)
            if (kbr == 2) call breit_ret(na,nc)
        end do
        cp1=CP(1)
        cq1=CQ(1)
        do i=1,ii,ih
            CP(i)=CP(i)+Qy(i)
            CQ(i)=CQ(i)+Py(i)
            C(i)=(Pa(i)*Py(i)+Qa(i)*Qy(i))*sn/R(i)
        end do
            C(ii+4)=Pa(ii+4)+Py(ii+4)-1
        call Sint1(e_breit)
        eadd(na)=e_breit
        
        cp1=CP(1)/cp1
        cq1=CQ(1)/cq1
        
        do i=ii+1,ii+10
            CP(i)=CP(i)*cp1
            CQ(i)=CQ(i)*cq1
        end do

        ch1=let(la+1)
        strfmt = '(2X,"E_breit(",I2,A1,"_",I1,"/2 ) = ",E14.7)'
        if (kout.GT.1) write( *,strfmt) Nn(na),ch1,ja,e_breit
        if (kout.GT.0) write(11,strfmt) Nn(na),ch1,ja,e_breit

        Return
    End Subroutine Breit_Core

    Subroutine Intens(NI,KI,DEL)
        Implicit None
        Integer :: ni, ki, ih, nfi, nj, ndi, i, nii, n12, l, iconf
        Real(dp) :: del, eps, e0, d0, d1, e1, de, dn, deli, p1, p2, p3, a1, a2, a3, b1, b2, b3, s1, s2, x, y
        Character(Len=256) :: strfmt

        EPS=EPS1
        IF (KT.EQ.1) EPS=EPS0
        ICONF=NC(NI)
        IH=2-KT
        N12=NI+4
        CALL READF (12,N12,A,B,2)
        E0=A(II+1)
        D0=(A(II+5)**2+B(II+5)**2)/A(II+2)
        E1=E0
        D1=D0
        DE=A(II+16)
        DN=A(II+17)
        DELI=DE
        IF (DN.GT.DE) DELI=DN
        IF (NF(NI).EQ.1) GOTO 200
        E0=P(II+1)
        DE=DABS((E0-E1)/E0)
        D0=(P(II+5)**2+Q(II+5)**2)/P(II+2)
        DN=DABS((D0-D1)/D0)
        DELI=DE
        IF (DN.GT.DE) DELI=DN
        P(II+16)=DE
        P(II+17)=DN
        IF (DELI.GT.DEL) DEL=DELI
  200   NFI=0
        IF (DELI.GE.0.1D0*EPS) GOTO 210
        NFI=1
        GOTO 230
  210   IF (KI.NE.1) GOTO 220

        ! uskorenie shodimosti
        N12=NI+NS+4
        CALL READF (12,N12,CP,CQ,2)
        P1=1.D0/DSQRT(CP(II+5)**2+CQ(II+5)**2)
        P2=1.D0/DSQRT( A(II+5)**2+ B(II+5)**2)
        P3=1.D0/DSQRT( P(II+5)**2+ Q(II+5)**2)
        Do I=1,II,IH
            A1=CP(I)*P1
            A2= A(I)*P2
            A3= P(I)*P3
            B1=CQ(I)*P1
            B2= B(I)*P2
            B3= Q(I)*P3
            S1=A1-2*A2+A3
            S2=B1-2*B2+B3
            X=0.5D0
            Y=0.5D0
            IF (DABS(S1).GT.1.D-10*DABS(A3)) X=(A3-A2)/S1
            IF (DABS(S2).GT.1.D-11*DABS(B3)) Y=(B3-B2)/S2
            IF (X.LE.0.D0) X=0.D0
            IF (Y.LE.0.D0) Y=0.D0
            IF (X.GE.0.7D0) X=0.7D0
            IF (Y.GE.0.7D0) Y=0.7D0
            A(I)=(X*A2+(1.D0-X)*A3)/P3
            B(I)=(Y*B2+(1.D0-Y)*B3)/P3
        End Do
        Do I=1,IP6
            CP(I)=P(I)
            CQ(I)=Q(I)
        End Do
        Do I=1,II,IH
            P(I)=A(I)
            Q(I)=B(I)
        End Do
        CALL NORM
        Do I=1,IP6
            A(I)=P(I)
            B(I)=Q(I)
            P(I)=CP(I)
            Q(I)=CQ(I)
        End Do
  220   N12=NI+NS+4
        CALL WRITEF(12,N12,A,B,2)  
  230   IF (KL.EQ.3) GOTO 240
        NJ=NI+1
        IF (NJ.GT.NS) GOTO 240
        IF (NC(NJ).NE.NC(NI)) GOTO 240
        IF (KP(NJ).EQ.1) GOTO 240
        IF (NN(NJ).NE.NN(NI)) GOTO 240
        IF (LL(NJ).NE.LL(NI)) GOTO 240
        N12=NS+NS+5
        CALL WRITEF(12,N12,P,Q,2)
        GOTO 260
  240   CALL DENS(NI,ICONF)
        N12=NI+4
        CALL WRITEF(12,N12,P,Q,2)
        ND=NS
        Do NII=1,NS
            IF (NC(NII).NE.NC(NI)) Cycle
            IF (NII.LE.NI) Cycle
            NDI=NII-NI
            IF (NDI.LT.ND) ND=NDI
        End Do
        IF (NI0.GT.NSMAX) GOTO 250
        IF (NI.GE.NI0) NI0=NI+ND
        IF (NI.EQ.NSMAX) NI0=NSMIN
  250   IF (KL.EQ.3) GOTO 260
        NJ=NI-1
        IF (NJ.LT.1) GOTO 260
        IF (NC(NJ).NE.NC(NI)) GOTO 260
        IF (KP(NJ).EQ.1) GOTO 260
        IF (NN(NJ).NE.NN(NI)) GOTO 260
        IF (LL(NJ).NE.LL(NI)) GOTO 260
        N12=NS+NS+5
        CALL READF (12,N12,P,Q,2)
        CALL DENS(NJ,ICONF)
        N12=NJ+4
        CALL WRITEF(12,N12,P,Q,2)  
  260   CALL READF (12,1,A,B,2)
        A(15)=NI0
        CALL WRITEF(12,1,A,B,2)
        NF(NI)=NFI
        L=LL(NI)+1
        strfmt = '(1X,I2,A1,I2,"/2 ",F6.3,F13.6,F10.6,2X,E13.6,1X,2E8.1,1X,I2,2(2X,I3))'
        WRITE(11,strfmt) NN(NI),LET(L),JJ(NI),QQ(NI),E0,-eadd(ni),D0,DE,DN
      
        NIT=0
        M2=0
        M3=0
        Return
    End Subroutine Intens

    Subroutine PQDif(NI)
        Implicit None

        Integer :: ni, j, n, im, kpi, k, nj, ih, i, m, l, iconf, m2
        Real(dp) :: c1, de, rm2, d0, sa, sb, r1gm, gm, xq, yi, qi, pi, rk, cr, t, ep, xp, eq, e, d

        Character(Len=256) :: strfmt
        Character*1 str_is(5)*3,str_br(3)*7

        data str_is /' NO',' VS','SMS','NMS',' MS'/
        data str_br /'Coulomb','Gaunt  ','Breit  '/

        ! wy~islenie perwyh proizwodnyh
        C1=0.01D0
        IH=2-KT
        Iconf=0
        If (ni.gt.0) ICONF=NC(NI)
        Do NJ=1,NS
            If (NC(NJ).NE.ICONF) Cycle
            If (NI.GT.NJ) Exit
            strfmt = '(/2X,76("=")/3X,16A1,2X,"ICONF =",I2,/3X,"Interaction: ",A7,"; R_N=",E11.4, &
                    "; Isotope shift: ",A3,"; c_is=",F7.4,/3X,"N",6X,"JJ",3X,"QQ",2X,"KP",4X, &
                    "EE (a.u.)",8X,"dE",8X,"D(0)",8X,"DELT",5X,"R(M2)"/2X,76("-"))'
            WRITE( *,strfmt) NAME,ICONF,str_br(kbr+1),RNUCL,str_is(k_is+1),c_is
            WRITE(11,strfmt) NAME,ICONF,str_br(kbr+1),RNUCL,str_is(k_is+1),c_is
            Exit
        End Do

        K=KK(NI)
        KPI=KP(NI)
        IF (KPI.EQ.1) Then
            Continue
        Else
            EQ=P(II+1)/CL
            EP=CL+CL-EQ
            D=1.D0/DSQRT(P(II+2))
            T=1.D0/DSQRT(Q(II+2))
            Do I=1,II,IH
                CR=1.D0/(CL*R(I))
                RK=K/R(I)
                PI=P(I)*D
                QI=Q(I)*D
                YI=Y(I)*CR
                XP=CP(I)*CR*T
                XQ=CQ(I)*CR*T
                A(I)=-( (EP-YI)*QI+RK*PI-XP)
                B(I)= (-(EQ+YI)*PI+RK*QI-XQ)
            End Do

            !>>>> origin (added 20/11/04):
            gm=P(ii+4)
            A(ii+4)=gm-1.d0
            B(ii+4)=gm-1.d0 ! cAB
            sa=0.d0
            sb=0.d0
            do m=0,Nmax
                im=ii+5+m
                A(im)=(gm+m)*P(im)
                sa=sa+A(im)
                B(im)=(gm+m)*Q(im)
                sb=sb+B(im)
            end do
            r1gm=R(1)**(gm-1.d0)
            sa=sa*r1gm
            sb=sb*r1gm
            sa=(A(1)-sa)/sa
            sb=(B(1)-sb)/sb
            if (dabs(sa)+dabs(sb).GT.3.d-3) then
                write(*,*) ' PQdif: matching error for ni=',ni
                write(*,*) ' mismatch for A:',sa,' for B:',sb
                read(*,*)
            end if
            !<<<< origin
            CALL WRITEF(12,NI+4+NS,A,B,2)
        End If
        N =NN(NI)
        L =LL(NI)
        J =JJ(NI)
        QI=QQ(NI)
        !!!!!!!!!!!!!!! E is in  a.u.  !!!!!!!!!!!!!!!!!!!
        E =P(II+1)
        D0=P(II+5)**2+Q(II+5)**2
        ! M2=Q(II+3)+C1 ! cAB
        M2=Q(II+15)+C1
        RM2=R(M2)
        DE=P(II+16)
        IF (DE.LT.P(II+17)) DE=P(II+17)
        if (k_is.EQ.1) call add_vs(ni)
        strfmt = '(I4,1X,I2,A1,I2,"/2 ",F6.3,1X,I1,F15.8,F12.8,1X,E12.5,E10.2,F7.3)'
        WRITE( *,strfmt) NI,N,LET(L+1),J,QI,KPI,E,-eadd(ni),D0,DE,RM2
        WRITE(11,strfmt) NI,N,LET(L+1),J,QI,KPI,E,-eadd(ni),D0,DE,RM2

        IF (NI.EQ.NSMAX) Then
            strfmt = '(2X,76("="))'
            WRITE( *,strfmt)
            WRITE(11,strfmt)
        End If

        If (NI.LT.NS) Then
            Continue
        Else
            CALL READF (12,1,A,B,2)
            NI0=0
            A(15)=NI0
            CALL WRITEF(12,1,A,B,2)
        End If

        Return
    End Subroutine PQDif

    Subroutine Test
        Implicit None

        Integer :: itest, ns1, ng, nso1, if, ni, nj, njj, nsm, nj1, li, n1, k1, n12
        Real(dp) :: c1, z1, r21, am1, c2, q1
        Integer, Dimension(IPns) :: Nn1, Kk1, Kp1, Nc1
        Real(dp), Dimension(IPns) :: Qq1
        Real(dp) :: Jm1
        Character(Len=256) :: strfmt

        C1=0.01D0
        Nn1=0
        Kk1=0
        Kp1=0
        Nc1=0
        Qq1=0_dp

        ITEST=0
        CALL READF (12,1,P,P,1)
        Z1  =P(1)
        NS1  =P(2)+C1
        R1   =P(4)
        R21  =P(5)
        H    =P(6)
        BT   =P(7)
        AL   =P(8)
        KT   =P(9)+C1
        NG   =P(10)+C1
        NSO1 =P(11)+C1
        AM1  =P(12)
        RNUCL=P(13)
        JM1  =P(14)
        NI0  =P(15)+C1

        IF (DABS(Z  -  Z1       ).GT.EPS0) GOTO 300
        IF (DABS(NS - NS1 + 0.d0).GT.EPS0) GOTO 300
        IF (DABS(R2 - R21 + 0.d0).GT.EPS0) GOTO 300
        IF (DABS(NSO-NSO1 + 0.d0).GT.EPS0) GOTO 300
        IF (DABS(AM - AM1       ).GT.EPS0) GOTO 300
        IF (DABS(JM - JM1 + 0.d0).GT.EPS0) GOTO 300
        GOTO 360

  300   ITEST=1
        strfmt = '(/2X,"INPUT DATA WERE CHANGED"//10X,"HFD.INP",26X,"HFD.DAT"/)'
        WRITE( *,strfmt)
        WRITE(11,strfmt)

        IF (DABS(NS-NS1+0.d0).LT.EPS0) GOTO 310
        strfmt = '(2X,"  NS  =",I5,25X,I5)'
        WRITE( *,strfmt) NS,NS1
        WRITE(11,strfmt) NS,NS1

  310   IF (DABS(NSO-NSO1+0.d0).LT.EPS0) GOTO 320
        strfmt = '(2X," NSO  =",I5,25X,I5)'
        WRITE( *,strfmt) NSO,NSO1
        WRITE(11,strfmt) NSO,NSO1

  320   IF (DABS(Z-Z1).LT.EPS0) GOTO 330
        strfmt = '(2X,"   Z  =",F12.6,18X,F12.6)'
        WRITE( *,strfmt) Z,Z1
        WRITE(11,strfmt) Z,Z1
  
  330   IF (DABS(R2-R21).LT.EPS0) GOTO 340
        strfmt = '(2X,"  R2  =",F12.6,18X,F12.6)'
        WRITE( *,strfmt) R2,R21
        WRITE(11,strfmt) R2,R21

  340   IF (DABS(JM-JM1+0.d0).LT.EPS0) GOTO 350
        strfmt = '(2X,"  JM  =",F12.6,18X,F12.6)'
        WRITE( *,strfmt) JM,JM1
        WRITE(11,strfmt) JM,JM1

  350   IF (DABS(AM-AM1).LT.EPS0) GOTO 360
        strfmt = '(2X,"  AM  =",F12.6,18X,F12.6)'
        WRITE( *,strfmt) AM,AM1
        WRITE(11,strfmt) AM,AM1

        If (KL.EQ.1) Then
            IF (ITEST.NE.0.AND.KL.EQ.1) STOP
            Return
        End If

  360   IF=20
        Do NJ=1,NS1
            IF=IF+1
            NN1(NJ)=P(IF)+C1
            IF=IF+1
            IF=IF+1
            QQ1(NJ)=P(IF)
            IF=IF+1
            IF=IF+1
            C2=DSIGN(C1,P(IF))
            KK1(NJ)=P(IF)+C2
            IF=IF+1
        End Do
        Do NJ=1,NS1
            IF=IF+1
            C2=DSIGN(C1,P(IF))
            KP1(NJ)=kp(nj) !P(IF)+C2 - changed 7/12/98 ### these changes ####
            IF=IF+1        !                           ### are made to be ###
            NC1(NJ)=nc(nj) !P(IF)+C1 - changed 7/12/98 ### able to use ######
        End Do             !                           ### bass orbitals ####

        Do NI=1,NS
            NJ=0
            Do NJJ=1,NS1
                IF (NC1(NJJ).NE.NC(NI)) Cycle
                IF (NN1(NJJ).NE.NN(NI)) Cycle
                IF (KK1(NJJ).NE.KK(NI)) Cycle
                NJ=NJJ
                IF (KP1(NJ).LT.0) KP(NI)=-1
                Exit
            End Do
              
            IF (NI.EQ.NJ) Cycle
            IF (NJ.GT.0) then        !### this allows not to call QQ1(0)
              if (DABS(QQ(NI)-QQ1(NJ)).LT.EPS0) GOTO 390
            end if

            IF (ITEST.EQ.2) GOTO 380
            strfmt = '(/2X,"CONFIGURATION WAS CHANGED"//10X,"HFD.INP",26X,"HFD.DAT"//8X,"NL",2X,"JJ",4X,"QQ",2X,"NC",16X,"NL",2X,"JJ",4X,"QQ",2X,"NC")'
            WRITE( *,strfmt)
            WRITE(11,strfmt)

            ITEST=2
            IF (KL.EQ.1) Cycle

            NSM=NS1
            IF (NS1.GE.NS) GOTO 380
            NSM=NS
            Do NJ1=1,NS1
                NJJ=NS1+1-NJ1
                N12=NJJ+NS1+4
                CALL READF (12,N12,A,B,2)
                N12=NJJ+NS+4
                CALL WRITEF (12,N12,A,B,2)
            End Do

      380   IF (NJ.EQ.0) GOTO 390
            LI=LL(NI)
            strfmt = '(I5,2X,I2,A1,I2,"/2"," (",F6.3,")",3X,I1,10X,I4,2X,I2,A1,I2,"/2"," (",F6.3,")",3X,I1)'
            WRITE( *,strfmt) NI,NN(NI),LET(LI+1),JJ(NI),QQ(NI),NC(NI),NI,NN(NI),LET(LI+1),JJ(NI),QQ1(NJ),NC(NI)
            WRITE(11,strfmt) NI,NN(NI),LET(LI+1),JJ(NI),QQ(NI),NC(NI),NI,NN(NI),LET(LI+1),JJ(NI),QQ1(NJ),NC(NI)
            GOTO 400

      390   strfmt = '(I5,2X,I2,A1,I2,"/2"," (",F6.3,")",3X,I1,20X,"----")'
            WRITE( *,strfmt) NI,NN(NI),LET(LL(NI)+1),JJ(NI),QQ(NI),NC(NI)
            WRITE(11,strfmt) NI,NN(NI),LET(LL(NI)+1),JJ(NI),QQ(NI),NC(NI)
            KP(NI)=-1
            Cycle

      400   IF (NI.GT.NS1) GOTO 410
            NJJ=NI
            N1=NN1(NJ)
            NN1(NJ)=NN1(NI)
            NN1(NI)=N1
            K1=KK1(NJ)
            KK1(NJ)=KK1(NI)
            KK1(NI)=K1
            Q1=QQ1(NJ)
            QQ1(NJ)=QQ1(NI)
            QQ1(NI)=Q1
            N1=NC1(NJ)
            NC1(NJ)=NC1(NI)
            NC1(NI)=N1

            N12=NJ+4
            CALL READF (12,N12,P,Q,2)
            N12=NI+4
            CALL READF (12,N12,A,B,2)
            N12=NJ+4
            CALL WRITEF(12,N12,A,B,2)
            N12=NI+4
            CALL WRITEF(12,N12,P,Q,2)

            N12=NJ+NSM+4
            CALL READF (12,N12,P,Q,2)
            N12=NI+NSM+4
            CALL READF (12,N12,A,B,2)
            N12=NJ+NSM+4
            CALL WRITEF(12,N12,A,B,2)
            N12=NI+NSM+4
            CALL WRITEF(12,N12,P,Q,2)
            Cycle

      410   N12=NJ+4
            CALL READF (12,N12,P,Q,2)
            N12=NI+4
            CALL WRITEF(12,N12,P,Q,2)
            N12=NJ+NSM+4
            CALL READF (12,N12,P,Q,2)
            N12=NI+NSM+4
            CALL WRITEF(12,N12,P,Q,2)
        End Do

        IF (NS.GE.NS1) Then
            Continue
        Else
            Do NJJ=1,NS1
                N12=NJJ+NS1+4
                CALL READF (12,N12,A,B,2)
                N12=NJJ+NS+4
                CALL WRITEF (12,N12,A,B,2)
            End Do
        End If

        IF (ITEST.NE.0.AND.KL.EQ.1) STOP

        RETURN
    End Subroutine Test

    Subroutine Tbr
        Implicit None

        Integer :: i
        Real(dp) :: al2, bt2, h2, d2, t2, p2

        AL2=AL
        BT2=BT
        H2=H
        R(1)=R1
        D2=R1
        T2=AL2*D2+BT2
        V(1)=D2/T2
        P2=AL2*D2+BT2*DLOG(D2)

        Do I=2,II
            P2=P2+H2
  200       T2=AL2*D2+BT2*DLOG(D2)
            T2=(P2-T2)/(AL2*D2+BT2)
            D2=D2*(1.D0+T2)
            IF (DABS(T2).GT.0.0005D-5) GOTO 200
            T2=AL2*D2+BT2
            V(I)=D2/T2
            R(I)=D2
        End Do

        Return
    End Subroutine Tbr

    Subroutine Yk(K)
        Implicit None

        Integer :: k, ih, i0, i1, im, j, ip, id, i, m
        Real(dp) :: h1, r0, v0, dk1, dk2, g, g0, t0, p0, f0, fm, t, f, dh, s, r0d, dt, d

        ! CALCULATION FUNCTIONS Z(K,R) AND Y(K,R)
        IH=2-KT
        H1=H*IH
        I0=1
        IMAX=ROIJ(II+3)+0.01D0
        R0=R(I0)
        V0=V(I0)
        DK1=1.D0-DEXP(-K*H1)
        DK2=1.D0-DEXP(-(K+1)*H1)
        IF (R2.LT.0.D0) Then
            Continue
        Else
            I1=I0+IH
            Do I=I1,II,IH
                W(I)=0.D0
                IF (K.NE.0) W(I)=1.D0-(R(I-IH)/R(I))**K
            End Do
        End If

        G=ROIJ(II+4)+K
        T0=0.D0
        P0=0.D0
        Do M=0,NMAX
            I=II+5+M
            T0=T0+ROIJ(I)/(G+M+1)
            P0=P0+ROIJ(I)*(G+M)
        End Do
        T0=T0*R0**(G+1-K)
        P0=P0*R0**(G-1-K)
        F0=ROIJ(I0)
        P0=H1*P0*V0*V0+F0*BT/(AL*R0+BT)**2
        FM=ROIJ(IMAX)

        ! FUNCTION Z(K,R)
        T=T0+0.5D0*F0+H1/12.D0*P0
        DT=T
        C(I0)=T
        D=DK1
        I1=I0+IH
        Do I=I1,IMAX,IH
            F=ROIJ(I)
            IF (R2.GT.0.D0) D=W(I)
            IF (K.NE.0) F=F-D*T
            DT=DT+F
            T=DT
            C(I)=T
        End Do

        ! FUNCTION Y(K,R)
        T=-0.5D0*FM
        DT=T
        DH=(2*K+1)*H1/12.D0
        D=DK2
        IM=IMAX-IH
        I1=IM+I0
        Do J=I0,IM,IH
            I=I1-J
            IP=I+IH
            S=DH
            IF (R2.LT.0.D0) Then
                Continue
            Else
                D=R(I)/R(IP)
                IF (K.NE.0) D=D*(1.D0-W(IP))
                D=1.D0-D
                S=DH*V(I)/R(I)
            End If
            F=ROIJ(IP)
            F=F-D*(T+F)
            DT=DT+F
            T=DT
            C(I)=C(I)+T-S*ROIJ(I)
        End Do
        D=DK1
        T=C(IMAX)-0.5D0*FM
        DT=T
        Do I=IMAX,II,IH
            C(I)=T
            IF (K.EQ.0) Cycle
            IF (R2.GT.0.D0) D=W(I)
            DT=DT-D*T
            T=DT
        End Do

        G0=ROIJ(II+4)
        ID=G0+0.5D0
        R0D=R0**ID
        T0=0.D0
        Do M=0,NMAX
            I=II+5+M
            T0=T0+(ROIJ(I)/(G0+M-K)-ROIJ(I)/(G0+M+K+1))
        End Do

        T0=T0*R0**G0
        Do M=0,NMAX
            I=II+1+M
            C(I)=0.D0
        End Do

        IF (K.LT.NMAX) C(II+1+K)=C(I0)/R0+T0
        Do M=0,NMAX
            I=II+5+M
            IF ((M+ID).LE.NMAX) C(II+1+M+ID)=ROIJ(I)*(1.D0/(G0+K+M+1)-1.D0/(G0+M-K))*R0D
        End Do

        Return
    End Subroutine Yk

    Subroutine PQInt(ICONF)
        Implicit None

        Integer :: is, i1, i2, iconf, j, ip, im1, im2, im3, ip1, ip2, jp2, jp4, jp6, jm2, jm4, i, n12, ni
        Real(dp) :: pi, qi

        ! 6 POINTS INTERPOLATION
        Do IS=1,NS
            NI=IS
            IF (NC(NI).NE.ICONF) Cycle
            IF (KP(NI).EQ.1) Cycle
            N12=NI+4
            CALL READF(12,N12,P,Q,2)
            I1=II-1
            Do I=2,I1,2
            P(I)=0.D0
            Q(I)=0.D0
            End Do
            IMAX=P(II+3)+0.01D0
            I1=IMAX-5
            Do I=2,I1,2
                J=I-1
                IF (J.LT.5) J=5
                IP=I-J
                IM1=IP-2
                IM2=IP-4
                IM3=IP-6
                IP1=IP+2
                IP2=IP+4
                JP2=J+2
                JP4=J+4
                JP6=J+6
                JM2=J-2
                JM4=J-4
                PI=0.1D0*(IM2*IM1*IP*IP1)*(IP2*P(JP6)-IM3*P(JM4))
                PI=PI+0.5D0*(IM3*IM1*IP*IP2)*(IM2*P(JM2)-IP1*P(JP4))
                PI=PI+(IM3*IM2*IP1*IP2)*(IP*P(JP2)-IM1*P(J))
                QI=0.1D0*(IM2*IM1*IP*IP1)*(IP2*Q(JP6)-IM3*Q(JM4))
                QI=QI+0.5D0*(IM3*IM1*IP*IP2)*(IM2*Q(JM2)-IP1*Q(JP4))
                QI=QI+(IM3*IM2*IP1*IP2)*(IP*Q(JP2)-IM1*Q(J))
                P(I)=PI/384.D0
                Q(I)=QI/384.D0
            End Do

            I1=IMAX-3
            I2=IMAX-1
            Do I=I1,I2,2
                J=I-1
                JP2=J+2
                JM2=J-2
                P(I)=(3*P(JP2)+6*P(J)-P(JM2))/8.D0
                Q(I)=(3*Q(JP2)+6*Q(J)-Q(JM2))/8.D0
            End Do
            N12=NI+4
            CALL WRITEF(12,N12,P,Q,2)
        End Do

        Return
    End Subroutine PQInt

    Subroutine add_vs(ni)         
        !### Calculates volume shift for orbital ni
        Implicit None   

        Integer :: ni, m, im, k, ik
        Real(dp) :: r1, gm, v0, de, s
 
        if (k_is.NE.1) return
 
        r1=r(1)
        gm=2*p(ii+4)
        V0=3*Z*r1**gm*c_is
 
        de=0.d0
        do m=0,nmax
          im=ii+m+5
          s=0.d0
          do k=0,m
            ik=ii+k+5
            s=s+p(ik)*p(im-k)+q(ik)*q(im-k)
          end do
          de=de+s/((gm+m+1)*(gm+m+3))
        end do
        eadd(ni)=V0*de

        Return
    End Subroutine add_vs

    Subroutine V_nms(na)  
        Use diff, Only : Dif5, Cut_Short
        !### Normal mass shift potential
        !### see 30/09/05 for derivation details

        Integer :: na, ih, ja, la, ila, ip, iq, i
        Real(dp) :: small, sn, gm, xja, cis2, clj, az, ri, pi, qi, e_sms
        Character :: ch1

        small=0.1d-8
        ih=2-kt

        sn=1.d0/P(ii+2)   !### normalization factor
        gm=P(ii+4)

        do i=1,IP6
            Pa(i)=P(i)
            Qa(i)=Q(i)
        end do

        call Dif5(Pa,Py,0,V,ii,kt,h)     ! first derivative
        call Dif5(Py,P1a,2,V,ii,kt,h)    ! second derivative
        if (klow.GT.0) then
            call Dif5(Qa,Qy,0,V,ii,kt,h)   ! first derivative
            call Dif5(Qy,Q1a,2,V,ii,kt,h)  ! second derivative
        end if

        ja=Jj(na)
        xja=0.5d0*ja
        la=Ll(na)
        ila=ja-la

        ip=la*(la+1)
        iq=ila*(ila+1)
        cis2=0.5d0*c_is
        clj=-xja*(xja+1)+la*(la+1)+0.75d0  !!!(4*la-2*xja+1)/(2*la+1)!!!
        aZ=Z/DPcl
        do i=1,ii,ih               ! (Py,Qy) = V_nms * (Pa,Qa)
            ri=R(i)
            if (klow.GE.2) then
                pi=aZ*(clj/ri**2*Q(i)-2/ri*Qy(i))
                qi=aZ*((clj-2)/ri**2*P(i)+2/ri*Py(i))
            else
                pi=0.d0
                qi=0.d0
            end if
            Py(i)=cis2*(-P1a(i)+ip*Pa(i)/ri**2+pi)
            if (klow.GT.0) then
                Qy(i)=cis2*(-Q1a(i)+iq*Qa(i)/ri**2+qi)
            else
                Qy(i)=0.d0
            end if
        end do

        call Cut_Short(Py,R,V,ii,kt,MaxT,h)
        if (klow.GT.0) call Cut_Short(Qy,R,V,ii,kt,MaxT,h)

        do i=1,ii,ih
            CP(i)=CP(i)+Qy(i)*R(i)
            CQ(i)=CQ(i)+Py(i)*R(i)
            C(i)=(P(i)*Py(i)+Q(i)*Qy(i))*sn
        end do
        C(ii+4)=P(ii+4)+Py(ii+4)
        call Sint1(e_sms)
        eadd(na)=e_sms

        ch1=let(la+1)
        strfmt = '(2X,"E_nms(",I2,A1,I1,") = ",E14.7)'
        if (kout.GT.1) write( *,strfmt) Nn(na),ch1,ja,e_sms
        if (kout.GT.0) write(11,strfmt) Nn(na),ch1,ja,e_sms

        Return
    End Subroutine V_nms

    Subroutine breit_mag(na,nc)
        Implicit None

        Integer :: l, lmin, lmax, na, nc, i, m
        Real(dp) :: gac
        Real(dp), Dimension(IP6) :: Pw, Qw

        ! Direct
        lmin=0
        lmax=int(2.d0*dmin1(xja,xjc))
        Do l=lmin,lmax
            gac=coeft(-1,l,na,nc) 
            Pw = 0.d0; Qw = 0.d0
            if (dabs(gac).gt.1.d-7) call mag_coulomb(l,na,Pa,Qa,nc,Pc,Qc,gac,Pw,Qw)

            Do i=1,ii
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do
            
            Do m=0,nmax
                i=ii+5+m
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do     
        End Do


        ! Exchange
        if (na == nc) return
        !RO=0_dp
        ! C=0_dp
        lmin=int(dabs(xja-xjc))
        lmax=int(xja+xjc)

        Do l=lmin,lmax
            gac=coeft(1,l,na,nc)
            Pw = 0.d0
            Qw = 0.d0
            
            if (dabs(gac).gt.1.d-7) call mag_exch(l,na,Pa,Qa,nc,Pc,Qc,gac,Pw,Qw)

            Do i=1,ii
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do

            Do m=0,nmax
                i=ii+5+m
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do
        End Do

        Return
    End Subroutine breit_mag

    Real(dp) Function coeft(i,k,na,nb)
        Implicit None

        Integer :: ni, i, k, na, nb, la, ja, lb, jb, n1, n2, j, n, maxns, k7, k11
        Real(dp) :: eps, g, wa, qb, ww, d, qa, y

        Integer, Dimension(70) :: ig
        Real(dp), Dimension(70) :: gg
        Integer, Dimension(IPns) :: nq


        Do ni = 1, IPns
            nq(ni) = qq(ni)+0.01d0
        End Do

        ! gg=2*ch(j1,0.5,j2,-0.5,k,0)**2/(2*k+1)

        gg( 1)=0.10000000000000d+01
        gg( 2)=0.33333333333333d+00
        gg( 3)=0.33333333333333d+00
        gg( 4)=0.20000000000000d+00
        gg( 5)=0.20000000000000d+00
        gg( 6)=0.14285714285714d+00
        gg( 7)=0.14285714285714d+00
        gg( 8)=0.11111111111111d+00
        gg( 9)=0.11111111111111d+00
        gg(10)=0.90909090909091d-01
        gg(11)=0.50000000000000d+00
        gg(12)=0.33333333333333d-01
        gg(13)=0.10000000000000d+00
        gg(14)=0.12857142857143d+00
        gg(15)=0.20000000000000d+00
        gg(16)=0.28571428571429d-01
        gg(17)=0.57142857142857d-01
        gg(18)=0.95238095238095d-01
        gg(19)=0.12857142857143d+00
        gg(20)=0.23809523809524d-01
        gg(21)=0.39682539682540d-01
        gg(22)=0.75757575757576d-01
        gg(23)=0.95238095238095d-01
        gg(24)=0.20202020202020d-01
        gg(25)=0.30303030303030d-01
        gg(26)=0.62937062937063d-01
        gg(27)=0.33333333333333d+00
        gg(28)=0.95238095238095d-02
        gg(29)=0.76190476190476d-01
        gg(30)=0.25396825396825d-01
        gg(31)=0.31746031746032d-01
        gg(32)=0.72150072150072d-01
        gg(33)=0.14285714285714d+00
        gg(34)=0.95238095238095d-02
        gg(35)=0.47619047619048d-01
        gg(36)=0.21645021645022d-01
        gg(37)=0.21645021645022d-01
        gg(38)=0.58275058275058d-01
        gg(39)=0.95238095238095d-01
        gg(40)=0.86580086580087d-02
        gg(41)=0.34632034632035d-01
        gg(42)=0.18648018648019d-01
        gg(43)=0.16317016317016d-01
        gg(44)=0.48951048951049d-01
        gg(45)=0.25000000000000d+00
        gg(46)=0.39682539682540d-02
        gg(47)=0.59523809523810d-01
        gg(48)=0.97402597402597d-02
        gg(49)=0.29220779220779d-01
        gg(50)=0.18731268731269d-01
        gg(51)=0.14568764568765d-01
        gg(52)=0.47591297591298d-01
        gg(53)=0.11111111111111d+00
        gg(54)=0.43290043290043d-02
        gg(55)=0.38961038961039d-01
        gg(56)=0.89910089910090d-02
        gg(57)=0.20979020979021d-01
        gg(58)=0.16317016317016d-01
        gg(59)=0.10878010878011d-01
        gg(60)=0.40312628547923d-01
        gg(61)=0.20000000000000d+00
        gg(62)=0.20202020202020d-02
        gg(63)=0.48484848484848d-01
        gg(64)=0.47952047952048d-02
        gg(65)=0.25174825174825d-01
        gg(66)=0.83916083916084d-02
        gg(67)=0.14918414918415d-01
        gg(68)=0.14333379039261d-01
        gg(69)=0.80625257095845d-02
        gg(70)=0.34371820130334d-01
        ig( 1)= 11
        ig( 2)=111
        ig( 3)=113
        ig( 4)=213
        ig( 5)=215
        ig( 6)=315
        ig( 7)=317
        ig( 8)=417
        ig( 9)=419
        ig(10)=519
        ig(11)= 33
        ig(12)=133
        ig(13)=233
        ig(14)=333
        ig(15)=135
        ig(16)=235
        ig(17)=335
        ig(18)=435
        ig(19)=237
        ig(20)=337
        ig(21)=437
        ig(22)=537
        ig(23)=339
        ig(24)=439
        ig(25)=539
        ig(26)=639
        ig(27)= 55
        ig(28)=155
        ig(29)=255
        ig(30)=355
        ig(31)=455
        ig(32)=555
        ig(33)=157
        ig(34)=257
        ig(35)=357
        ig(36)=457
        ig(37)=557
        ig(38)=657
        ig(39)=259
        ig(40)=359
        ig(41)=459
        ig(42)=559
        ig(43)=659
        ig(44)=759
        ig(45)= 77
        ig(46)=177
        ig(47)=277
        ig(48)=377
        ig(49)=477
        ig(50)=577
        ig(51)=677
        ig(52)=777
        ig(53)=179
        ig(54)=279
        ig(55)=379
        ig(56)=479
        ig(57)=579
        ig(58)=679
        ig(59)=779
        ig(60)=879
        ig(61)= 99
        ig(62)=199
        ig(63)=299
        ig(64)=399
        ig(65)=499
        ig(66)=599
        ig(67)=699
        ig(68)=799
        ig(69)=899
        ig(70)=999

        imax=70
        eps=1.d-4
        g=0.d0

        la=ll(na)
        ja=jj(na)
        qa=qq(na)
        wa=qw(na)
        lb=ll(nb)
        jb=jj(nb)
        qb=qq(nb)
        ww=(wa-1)/(4*la+1)

        if (k.gt.0.or.i.gt.0) Then
            Continue
        Else
            ! Coulomb
            if (qa.ge.eps) then
                if (na.eq.nb) then
                    g=0.5d0*((qa-wa)+ja*ww)
                else
                    g=qb
                    if (nc(na).eq.nc(nb).and.nn(na).eq.nn(nb).and.la.eq.lb) g=(qb+(jb+1)*ww)
                endif
            else
                if (na.eq.nb) g=0.5d0*qa
                if (na.ne.nb) g=qb          

            endif
            goto 1000
        End If

        ! Exchange
        If (na.ne.nb.and.i.lt.0) then
            Continue
        Else
            n2=k*100+jb*10+ja
            n1=k*100+ja*10+jb
            y=0.d0
            Do j=1,imax
                n=ig(j)
                if (n1.ne.n.and.n2.ne.n) Cycle
                y=gg(j)
                Exit
            End Do

            if (qa.ge.eps) then
                if (na.eq.nb) then
                    d=0.5d0*((qa-wa)*(ja+1.d0)/ja+ww*(ja+1))
                else
                    d=qb
                    if (nc(na).eq.nc(nb).and. nn(na).eq.nn(nb).and.la.eq.lb) d=(qb+ww*(jb+1))
                endif
            else
                d=qb
            endif
            g=-0.5d0*d*y
        End If

        if (ng.eq.0) goto 1000
        if (dabs(qa).lt.eps.or.dabs(qb).lt.eps) goto 1000
        if (nq(na).eq.ja+1.or.nq(nb).eq.jb+1) goto 1000

        n1=i*(maxns*maxns*k+na*maxns+nb)
        n2=i*(maxns*maxns*k+nb*maxns+na)

        Do j=1,ng
            n=igam(j)
            if (n.ne.n1.and.n.ne.n2) Cycle
            g=g+gam(j)/qa
            goto 1000
        End Do

        write( k7,'(/2x,a)') 'Gaunt coefficient did not found in coef.' 
        write(k11,'(/2x,a)') 'Gaunt coefficient did not found in coef.' 
        write( k7,'(2x,a,i4,2x,a,i4,2x,a,i2,2x,a,i2)') 'na=',na,'nb=',nb,'k=',k,'i=',i
        write(k11,'(2x,a,i4,2x,a,i4,2x,a,i2,2x,a,i2)') 'na=',na,'nb=',nb,'k=',k,'i=',i

        call exit(1)

1000    if (dabs(g).lt.1.d-7) g=0.d0
        coeft=g
        Return
    End Function coeft

    Subroutine breit_ret(na,nc)
        Implicit None

        Integer :: na, nc, lmin, lmax, l, i, m
        Real(dp) :: gac
        Real(dp), dimension(IP6) :: Pw,Qw
        
        ! Direct
        lmin=0
        lmax=int(2.d0*dmin1(xja,xjc))
        Do l=lmin,lmax
            gac=coeft(-1,l,na,nc) 
            Pw = 0.d0; Qw = 0.d0
            if (dabs(gac).gt.1.d-7) call ret_coulomb(l,na,Pa,Qa,nc,Pc,Qc,gac,Pw,Qw)

            Do i=1,ii
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do

            Do m=0,nmax
                i=ii+5+m
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do
        End Do 

        ! Exchange
        if (na == nc) return
        !RO=0_dp
        ! C=0_dp
        lmin=int(dabs(xja-xjc))
        lmax=int(xja+xjc)

        Do l=lmin,lmax
            gac=coeft(1,l,na,nc)
            Pw = 0.d0
            Qw = 0.d0
            if (dabs(gac).gt.1.d-7) call ret_exch(l,na,Pa,Qa,nc,Pc,Qc,gac,Pw,Qw)
            Do i=1,ii
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do

            Do m=0,nmax
                i=ii+5+m
                Py(i)=Py(i)+Qw(i)
                Qy(i)=Qy(i)+Pw(i)
            End Do
        End Do 

        Return
    End Subroutine breit_ret

    Subroutine mag_coulomb(l,na,pa,qa,nb,pb,qb,gab,wp,wq)
        Use breit, Only : Ykt
        Implicit None

        Integer :: l, na, nb, la, lb, ka, kb, k, n, j, id, i, m
        Real(dp) :: gab, gab1, gab2, c1, s, r0d, dpp, dq, d
        Real(dp), Dimension(IP6) :: pa, qa, pb, qb, wp, wq, ro=0_dp, c=0_dp

        c1=0.01

        la=ll(na)
        ka=kk(na)
        lb=ll(nb)
        kb=kk(nb)

        ! Calculation cross density
        call rho_pq(pb,qb,ro)
        
        Do k=l-1,l+1
            if (k.lt.0) Cycle
            if (l.eq.0.and.k.eq.0) Cycle
            i=la+la+k
            if (2*(i/2).eq.i) Cycle
            gab1=coefb(1,k,l,na,na)*coefb(1,k,l,nb,nb)
            gab2=coefb(-1,k,l,na,na)*coefb(1,k,l,nb,nb)
            if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle
            call ykt(k,nmax,ro,c,r,v,r2,ii,kt)

            d=-2*gab*(gab1+gab2)
            if (na.eq.nb) d=2*d
            do i=1,ii
              s=d*c(i)
              wp(i)=wp(i)+s*pa(i)
              wq(i)=wq(i)+s*qa(i)
            enddo

            id=iabs(kb)-iabs(ka)
            r0d=r(1)**id

            Do m=0,nmax
                i=ii+5+m
                dpp=0.d0
                dq=0.d0
                do n=0,m
                j=ii+5+n
                  dpp=dpp+pa(j)*c(i-n)
                  dq=dq+qa(j)*c(i-n)
                end do
                if (m+id.lt.0) Cycle
                if (m+id.gt.nmax) Cycle
                wp(i+id)=wp(i+id)+d*dpp*r0d
                wq(i+id)=wq(i+id)+d*dq*r0d
            End Do
        End Do

        Return
    End Subroutine mag_coulomb

    Real(dp) Function coefb(ibet,k,l,na,nb)
        Implicit None

        Integer :: ibet, k, l, na, nb, ka, kb, lk
        Real(dp) :: u, alk, blk, d

        u=0.d0

        ka=kk(na)
        kb=kk(nb)

        alk=0.5d0*(k*(k+1)-l*(l+1))
        lk=k+l
        blk=1.d0
        if (k.eq.l) blk=1.d0/dsqrt(2.d0)
        u=dsqrt((2*l+1.d0)/(lk*(lk+1)*(lk+2)))/blk
        d=kb
        if (2*(lk/2).ne.lk) d=-d
        u=u*(ka+d+ibet*alk)

        coefb=u

        Return
    End Function coefb

    Subroutine mag_exch(l,na,pa,qa,nb,pb,qb,gab,wp,wq)
        Use breit, Only : Ykt
        Implicit None

        Integer :: l, na, nb, la, lb, ka, kb, k, id, n, j, i, m
        Real(dp) :: gab, gab1, gab2, d1, d2, r0d, dpp, dq
        Real(dp), dimension(IP6) :: pa,qa,pb,qb,wp,wq,ro=0_dp,c=0_dp

        la=ll(na)
        ka=kk(na)
        lb=ll(nb)
        kb=kk(nb)

        ! Calculation cross density
        call rho_pq(pb,qa,ro)

        Do k=l-1,l+1
            if (k.lt.0) Cycle
            if (l.eq.0.and.k.eq.0) Cycle
            i=la+lb+k
            if (2*(i/2).eq.i) Cycle
            gab1=coefb(1,k,l,na,nb)*coefb(-1,k,l,na,nb)
            gab2=coefb(-1,k,l,na,nb)*coefb(-1,k,l,na,nb)
            if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle

            call ykt(k,nmax,ro,c,r,v,r2,ii,kt)

            d1=-2*gab*gab1
            d2=-2*gab*gab2
            do i=1,ii
              wq(i)=wq(i)+d1*c(i)*qb(i)
              wp(i)=wp(i)+d2*c(i)*pb(i)
            enddo

            id=iabs(kb)-iabs(ka)
            r0d=r(1)**id

            Do m=0,nmax
                i=ii+5+m
                dpp=0.d0
                dq=0.d0
                do n=0,m
                    j=ii+5+n
                    dpp=dpp+pb(j)*c(i-n)
                    dq=dq+qb(j)*c(i-n)
                enddo
                if (m+id.lt.0) Cycle
                if (m+id.gt.nmax) Cycle
                wq(i+id)=wq(i+id)+d1*dq*r0d
                wp(i+id)=wp(i+id)+d2*dpp*r0d
            End Do
        End Do

        ! Calculation cross density
        call rho_pq(pa,qb,ro)

        Do k=l-1,l+1
            if (k.lt.0) Cycle
            if (l.eq.0.and.k.eq.0) Cycle
            i=la+lb+k
            if (2*(i/2).eq.i) Cycle
            gab1=coefb(1,k,l,na,nb)*coefb(1,k,l,na,nb)
            gab2=coefb(-1,k,l,na,nb)*coefb(1,k,l,na,nb)
            if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle

            call ykt(k,nmax,ro,c,r,v,r2,ii,kt)

            d1=-2*gab*gab1
            d2=-2*gab*gab2
            do i=1,ii
                wq(i)=wq(i)+d1*c(i)*qb(i)
                wp(i)=wp(i)+d2*c(i)*pb(i)
            enddo

            id=iabs(kb)-iabs(ka)
            r0d=r(1)**id

            Do m=0,nmax
                i=ii+5+m
                dpp=0.d0
                dq=0.d0
                do n=0,m
                    j=ii+5+n
                    dpp=dpp+pb(j)*c(i-n)
                    dq=dq+qb(j)*c(i-n)
                end do
                if (m+id.lt.0) Cycle
                if (m+id.gt.nmax) Cycle
                wq(i+id)=wq(i+id)+d1*dq*r0d
                wp(i+id)=wp(i+id)+d2*dpp*r0d
            End Do
        End Do

        Return
    End Subroutine mag_exch

    Subroutine ret_coulomb(l,na,pa,qa,nb,pb,qb,gab,wp,wq)
        Use breit, Only : ykt, ukt, vkt
        Implicit None
        
        Integer :: l, na, nb, la, lb, ka, kb, k1, k2, k, id, n, j, i, m
        Real(dp) :: gab, c1, dk, ck1, ck2, gab1, gab2, c12, s, r0d, dpp, dq, d
        Real(dp), dimension(IP6) :: pa,qa,pb,qb,wp,wq,ro=0_dp,c=0_dp
        
        c1=0.01

        la=ll(na)
        ka=kk(na)
        lb=ll(nb)
        kb=kk(nb)

        ! Calculation cross density
        call rho_pq(pb,qb,RO)
        Do k1=l-1,l+1,2
            if (k1.lt.0) Cycle
            i=la+la+k1
            if (2*(i/2).eq.i) Cycle

            dk=k1
            if (k1.eq.l+1) then
                ck1=dsqrt(dk/(2*l+1))
            else
                ck1=-dsqrt((dk+1.d0)/(2*l+1))
            end if

            Do k2=l-1,l+1,2
                if (k2.lt.0) Cycle
                i=lb+lb+k2
                if (2*(i/2).eq.i) Cycle

                dk=k2
                if (k2.eq.l+1) then
                    ck2=dsqrt(dk/(2*l+1))
                else
                    ck2=-dsqrt((dk+1.d0)/(2*l+1))
                endif
                if (ck2.eq.0.d0) Cycle

                gab1=coefb( 1,k1,l,na,na)*coefb(-1,k2,l,nb,nb)
                gab2=coefb(-1,k1,l,na,na)*coefb(-1,k2,l,nb,nb)
                if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle

                If (k1.eq.k2) then
                    call ykt(k1,nmax,RO,C,r,v,r2,ii,kt)
                Else
                    if (k2.eq.k1+2) k=k1
                    if (k1.eq.k2+2) k=k2

                    if (k1.eq.k2+2) then
                        call ukt(k,nmax,RO,C,r,v,r2,ii,kt)
                    endif

                    if (k2.eq.k1+2) then
                        call vkt(k,nmax,RO,C,r,v,r2,ii,kt)
                    endif
                End If

                c12=ck1*ck2*dsqrt((2*k1+1.d0)*(2*k2+1.d0))
                if (k1.eq.k2) c12=2.d0*c12/(2*k1+1)

                d=-gab*(gab1+gab2)*c12
                if (na.eq.nb) d=2*d
                do i=1,ii
                    s=d* C(i)
                    wp(i)=wp(i)+s*pa(i)
                    wq(i)=wq(i)+s*qa(i)
                end do

                id=iabs(kb)-iabs(ka)
                r0d=r(1)**id

                Do m=0,nmax
                    i=ii+5+m
                    dpp=0.d0
                    dq=0.d0
                    Do n=0,m
                        j=ii+5+n
                        dpp=dpp+pa(j)*C(i-n)
                        dq=dq+qa(j)*C(i-n)
                    End Do
                    if (m+id.lt.0) Cycle
                    if (m+id.gt.nmax) Cycle
                    wp(i+id)=wp(i+id)+d*dpp*r0d
                    wq(i+id)=wq(i+id)+d*dq*r0d
                End Do
            End Do
        End Do

        Return
    End Subroutine ret_coulomb

    Subroutine ret_exch(l,na,pa,qa,nb,pb,qb,gab,wp,wq)
        Use breit, Only : ykt, ukt, vkt
        Implicit None
        
        Integer :: l, na, nb, la, lb, ka, kb, k1, k2, k, id, n, j, i, m
        Real(dp) :: gab, dk, ck1, ck2, gab1, gab2, c12, d1, d2, r0d, dpp, dq
        Real(dp), Dimension(IP6) :: pa,qa,pb,qb,wp,wq,ro=0_dp,c=0_dp

        la=ll(na)
        ka=kk(na)
        lb=ll(nb)
        kb=kk(nb)

        ! Calculation cross density
        call rho_pq(pb,qa,RO)
        Do k1=l-1,l+1,2
            if (k1.lt.0) Cycle
            i=la+lb+k1
            if (2*(i/2).eq.i) Cycle

            dk=k1
            if (k1.eq.l+1) then
                ck1=dsqrt(dk/(2*l+1))
            else
                ck1=-dsqrt((dk+1.d0)/(2*l+1))
            end if

            Do k2=l-1,l+1,2
                if (k2.lt.0) Cycle
                i=la+lb+k2
                if (2*(i/2).eq.i) Cycle

                dk=k2
                if (k2.eq.l+1) then
                    ck2=dsqrt(dk/(2*l+1))
                else
                    ck2=-dsqrt((dk+1.d0)/(2*l+1))
                endif
                if (ck2.eq.0.d0) Cycle

                gab1=coefb( 1,k1,l,na,nb)*coefb(1,k2,l,nb,na)
                gab2=coefb(-1,k1,l,na,nb)*coefb(1,k2,l,nb,na)
                if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle
                if (k1.eq.k2) then
                    call ykt(k1,nmax,RO, C,r,v,r2,ii,kt)
                else
                    if (k2.eq.k1+2) k=k1
                    if (k1.eq.k2+2) then
                        k=k2
                        call ukt(k,nmax,RO, C,r,v,r2,ii,kt)
                    end if

                    if (k2.eq.k1+2) then
                       call vkt(k,nmax,RO, C,r,v,r2,ii,kt)
                    endif
                endif
                

                c12=ck1*ck2*dsqrt((2*k1+1.d0)*(2*k2+1.d0))
                if (k1.eq.k2) c12=2.d0*c12/(2*k1+1)

                d1=-gab*gab1*c12
                d2=-gab*gab2*c12
                do i=1,ii
                    wq(i)=wq(i)+d1* C(i)*qb(i)
                    wp(i)=wp(i)+d2* C(i)*pb(i)
                enddo

                id=iabs(kb)-iabs(ka)
                r0d=r(1)**id

                Do m=0,nmax
                    i=ii+5+m
                    dpp=0.d0
                    dq=0.d0
                    do n=0,m
                        j=ii+5+n
                        dpp=dpp+pb(j)* C(i-n)
                        dq=dq+qb(j)* C(i-n)
                    end do
                    if (m+id.lt.0) Cycle
                    if (m+id.gt.nmax) Cycle
                    wq(i+id)=wq(i+id)+d1*dq*r0d
                    wp(i+id)=wp(i+id)+d2*dpp*r0d
                End Do
            End Do
        End Do
        ! Calculation cross density
        call rho_pq(pa,qb,RO)
        Do k1=l-1,l+1,2
            if (k1.lt.0) Cycle
            i=la+lb+k1
            if (2*(i/2).eq.i) Cycle

            dk=k1
            if (k1.eq.l+1) then
                ck1=dsqrt(dk/(2*l+1))
            else
                ck1=-dsqrt((dk+1.d0)/(2*l+1))
            endif

            Do k2=l-1,l+1,2
                if (k2.lt.0) Cycle
                i=la+lb+k2
                if (2*(i/2).eq.i) Cycle

                dk=k2
                if (k2.eq.l+1) then
                    ck2=dsqrt(dk/(2*l+1))
                else
                    ck2=-dsqrt((dk+1.d0)/(2*l+1))
                endif
                if (ck2.eq.0.d0) Cycle

                gab1=coefb( 1,k1,l,na,nb)*coefb(-1,k2,l,nb,na)
                gab2=coefb(-1,k1,l,na,nb)*coefb(-1,k2,l,nb,na)
                if (dabs(gab1).lt.1.d-7.and.dabs(gab2).lt.1.d-7) Cycle
                if (k1.eq.k2) then
                    call ykt(k1,nmax,RO, C,r,v,r2,ii,kt)
                else
                    if (k2.eq.k1+2) k=k1
                    if (k1.eq.k2+2) k=k2
                    if (k1.eq.k2+2) then
                        call ukt(k,nmax,RO, C,r,v,r2,ii,kt)
                    endif
                    if (k2.eq.k1+2) then
                        call vkt(k,nmax,RO, C,r,v,r2,ii,kt)
                    endif
                endif
                

                c12=ck1*ck2*dsqrt((2*k1+1.d0)*(2*k2+1.d0))
                if (k1.eq.k2) c12=2.d0*c12/(2*k1+1)

                d1=-gab*gab1*c12
                d2=-gab*gab2*c12
                do i=1,ii
                    wq(i)=wq(i)+d1* C(i)*qb(i)
                    wp(i)=wp(i)+d2* C(i)*pb(i)
                enddo

                id=iabs(kb)-iabs(ka)
                r0d=r(1)**id

                Do m=0,nmax
                    i=ii+5+m
                    dpp=0.d0
                    dq=0.d0
                    Do n=0,m
                        j=ii+5+n
                        dpp=dpp+pb(j)* C(i-n)
                        dq=dq+qb(j)* C(i-n)
                    End Do
                    if (m+id.lt.0) Cycle
                    if (m+id.gt.nmax) Cycle
                    wq(i+id)=wq(i+id)+d1*dq*r0d
                    wp(i+id)=wp(i+id)+d2*dpp*r0d
                End Do
            End Do
        End Do

        Return
    End Subroutine ret_exch

    Subroutine rho_pq(a,b,ro)
        Implicit None

        Integer :: imax, ih, n, j, i, m
        Real(dp) :: c1, h1, ulam, dh, dr
        Real(dp), Dimension(IP6) :: a, b, ro

        c1=0.01

        imax=a(ii+3)+c1
        if (b(ii+3) < imax) imax = b(ii+3)+c1

        ulam=1
        
        ih=2-kt
        h1=h*ih
        dh=ulam*h1

        Do i=1,imax
            ro(i)=a(i)*b(i)*dh*v(i)
        End Do

        Do i=imax+1,ii
          ro(i)=0.d0
        End Do

        ro(ii+3)=imax
        ro(ii+4)=a(ii+4)+b(ii+4)

        Do m=0,nmax
            i=ii+5+m
            dr=0.d0
            Do n=0,m
                j=ii+5+n
                dr=dr+(a(j)*b(i-n))
            End Do
            ro(i)=ulam*dr
        End Do

        Return
    End Subroutine rho_pq

End Program hfd