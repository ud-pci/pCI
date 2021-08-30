Module lsj_aux
    ! subroutines lsj, lsj_det, pls, pll, p0s, p0l
    Use conf_variables

    Implicit None

  Contains

    Subroutine lsj(cc,nmax,xj,xl,xs,mype,npes)
        Use mpi
        Use determinants, Only : Gdet
        Implicit None
        Real(dp), Allocatable, Dimension(:) :: cc
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Integer :: mype, npes, mpierr
        Integer :: n, n1, ic1, ic2, k, k1, k1n, ndn, ndk, nmax
        Real(dp) :: xj, xl, xs, ckn, tj, tl, ts

        xj=0.d0
        xl=0.d0
        xs=0.d0
        n=0

        Allocate(idet1(Ne),idet2(Ne))
        Do ic1=mype+1,Nc,npes
            ndn=Ndc(ic1)
            n=sum(Ndc(1:ic1-1))
            Do n1=1,ndn
                n=n+1
                call Gdet(n,idet1)
                k=n-1
                Do ic2=ic1,Nc
                    ndk=Ndc(ic2)
                    k1n=1
                    If (ic2.EQ.ic1) k1n=n1
                    Do k1=k1n,ndk
                        k=k+1
                        call Gdet(k,idet2)
                        call lsj_det(idet1,idet2,tj,tl,ts)
                        ckn=cc(n)*cc(k)
                        If (n.ne.k) ckn=2*ckn
                        xj=xj+ckn*tj
                        xl=xl+ckn*tl
                        xs=xs+ckn*ts
                    End Do
                End Do
            End Do
        End Do
        Call MPI_AllReduce(xj, xj, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(xl, xl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_AllReduce(xs, xs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        Deallocate(idet1,idet2)
        Return
    End Subroutine lsj


    Subroutine lsj_det(idet1,idet2,tj,tl,ts)
        Use formj2, Only : Plj
        Use determinants, Only : Rspq
        Implicit None
        Integer, Allocatable, Dimension(:) :: idet1, idet2
        Real(dp) :: tj, tl, ts
        Integer :: ic, id, is, nf, i1, i2, j1, j2, iq, ia, ja, la, na, ma, j, jq0, jq, ib

        tj=0.d0
        tl=0.d0
        ts=0.d0
        call Rspq(idet1,idet2,is,nf,i1,j1,i2,j2)
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!        Determinants are equal
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.0) then
          tj=mj*mj
          Do iq=1,Ne
            ia=idet1(iq)
            na=Nh(ia)
            ja=Jj(na)
            la=Ll(na)
            ma=Jz(ia)
            tj=tj+ja*(ja+2)-ma**2
            ts=ts+0.75d0
            tl=tl+la*(la+1)
            jq0=iq+1
            If (jq0.LE.Ne) then
              Do jq=jq0,Ne
                ib=idet1(jq)
                tj=tj-Plj(ia,ib)**2-Plj(ib,ia)**2
                ts=ts+2*(p0s(ia,ia)*p0s(ib,ib)-p0s(ia,ib)**2)
                ts=ts-2*(pls(ia,ib)**2+pls(ib,ia)**2)
                tl=tl+2*(p0l(ia,ia)*p0l(ib,ib)-p0l(ia,ib)**2)
                tl=tl-2*(pll(ia,ib)**2+pll(ib,ia)**2)
              End Do
            End If
          End Do
        End If
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!       One function differs in the Determinants.
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.1) then
          ia=i2
          ib=j2
          Do iq=1,ne
            If (iq.eq.j) Cycle
            ic=idet1(iq)
            ts=ts+2*(p0s(ia,ib)*p0s(ic,ic)-p0s(ia,ic)*p0s(ib,ic))
            ts=ts-2*(pls(ia,ic)*pls(ib,ic)+pls(ic,ia)*pls(ic,ib))
            tl=tl+2*(p0l(ia,ib)*p0l(ic,ic)-p0l(ia,ic)*p0l(ib,ic))
            tl=tl-2*(pll(ia,ic)*pll(ib,ic)+pll(ic,ia)*pll(ic,ib))
          End Do
          ts=ts*is
          tl=tl*is
        End If
!       - - - - - - - - - - - - - - - - - - - - - - - - -
!        Determinants differ by two functions
!       - - - - - - - - - - - - - - - - - - - - - - - - -
        If (nf.EQ.2) then
           ia=i1
           ic=j1
           ib=i2
           id=j2
           tj=Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
             Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic)
           ts=ts+2*(p0s(ia,ic)*p0s(id,ib)-p0s(ia,id)*p0s(ic,ib))
           ts=ts+2*(pls(ia,ic)*pls(id,ib)+pls(ic,ia)*pls(ib,id)- &
              pls(ia,id)*pls(ic,ib)-pls(id,ia)*pls(ib,ic))
           tl=tl+2*(p0l(ia,ic)*p0l(id,ib)-p0l(ia,id)*p0l(ic,ib))
           tl=tl+2*(pll(ia,ic)*pll(id,ib)+pll(ic,ia)*pll(ib,id)- &
              pll(ia,id)*pll(ic,ib)-pll(id,ia)*pll(ib,ic))
           tj=tj*is
           ts=ts*is
           tl=tl*is
        End If
        ts=ts*4
        tl=tl*4

        Return
    End Subroutine lsj_det

    Real(dp) Function pls(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)+2) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        msa=1
        msb=-1
        mja=jz(ia)
        mjb=jz(ib)
        mla=(mja-msa)/2
        mlb=mla
        If (iabs(mla).gt.la) goto 1000
        ta=0.d0
        If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
        If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
        tb=0.d0
        If (jb.eq.2*lb+1) tb= dsqrt((jb-mjb)/(2.d0*jb))
        If (jb.eq.2*lb-1) tb= dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
        t=-dsqrt(0.5d0)*ta*tb

1000    pls=t
        Return
    End Function pls

    Real(dp) Function pll(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)+2) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        mja=jz(ia)
        mjb=jz(ib)
        Do msa=-1,1,2
            msb=msa
            mla=(mja-msa)/2
            mlb=(mjb-msb)/2
            If (iabs(mla).gt.la) Cycle
            If (iabs(mlb).gt.lb) Cycle
            ta=0.d0
            tb=0.d0
            If (msa.eq.1) then
                If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            If (msa.eq.-1) then
                If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t-dsqrt(0.5d0*(lb*(lb+1.d0)-mlb*(mlb+1.d0)))*ta*tb
        End Do
1000    pll=t
        Return
    End Function pll

    Real(dp) Function p0s(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
        ! We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
        ! but maybe different ja and jb.
        t=0.
        If (jz(ia).ne.jz(ib)) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        Do msa=-1,1,2
            msb=msa
            mja=jz(ia)
            mjb=jz(ib)
            mla=(mja-msa)/2
            mlb=mla
            If (iabs(mla).gt.la) Cycle
                ta=0.d0
                tb=0.d0
                If (msa.eq.1) Then
                    If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                    If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                    If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                    If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
                End If
                If (msa.eq.-1) then
                    If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                    If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                    If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                    If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t+0.5d0*msa*(ta*tb)
        End Do

1000    p0s=t
        Return
    End Function p0s

    Real(dp) Function p0l(ia,ib)
        Implicit None
        Integer :: ia, ib, na, nb, ja, jb, la, lb, msa, msb, mja, mjb, mla, mlb
        Real(dp) :: t, ta, tb
!       We suggest that Int(Pa*Pb)=1.d0 for na=nb and la=lb
!       but maybe different ja and jb.

        t=0.
        If (jz(ia).ne.jz(ib)) goto 1000
        na=nh(ia)
        nb=nh(ib)
        If (nn(na).ne.nn(nb)) goto 1000
        If (ll(na).ne.ll(nb)) goto 1000
        ja=jj(na)
        jb=jj(nb)
        la=ll(na)
        lb=ll(nb)
        mja=jz(ia)
        mjb=jz(ib)
        Do msa=-1,1,2
            msb=msa
            mla=(mja-msa)/2
            mlb=mla
            If (iabs(mla).gt.la) Cycle
            ta=0.d0
            tb=0.d0
            If (msa.eq.1) then
                If (ja.eq.2*la+1) ta= dsqrt((ja+mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=-dsqrt((ja-mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb= dsqrt((jb+mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=-dsqrt((jb-mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            If (msa.eq.-1) then
                If (ja.eq.2*la+1) ta=dsqrt((ja-mja)/(2.d0*ja))
                If (ja.eq.2*la-1) ta=dsqrt((ja+mja+2.d0)/(2.d0*ja+4.d0))
                If (jb.eq.2*lb+1) tb=dsqrt((jb-mjb)/(2.d0*jb))
                If (jb.eq.2*lb-1) tb=dsqrt((jb+mjb+2.d0)/(2.d0*jb+4.d0))
            End If
            t=t+mla*(ta*tb)
        End Do

1000    p0l=t
        Return
    End Function p0l

    Subroutine PrintLSJ
        Implicit None

        Integer :: j1, nk, k, j2, n, ndk, ic, i, j, idum, ist, jmax, imax, &
                   j3
        real(dp) :: cutoff, xj, dt, del, dummy, wmx, E, D
        real(dp), allocatable, dimension(:)  :: Cc, Dd
        Character(Len=1), dimension(11) :: st1, st2 
        Character(Len=1), dimension(10)  :: stecp*7
        Character(Len=1), dimension(2)  :: st3*3
        Character(Len=1), dimension(4)  :: strsms*6
        Character(Len=1), dimension(3)  :: strms*3
        data st1/11*'='/,st2/11*'-'/,st3/' NO','YES'/
        data stecp/'COULOMB','C+MBPT1','C+MBPT2', &
                   'GAUNT  ','G+MBPT1','G+MBPT2', &
                   'BREIT  ','B+MBPT1','B+MBPT2','ECP    '/
        data strsms/'(1-e) ','(2-e) ','(full)','      '/
        data strms/'SMS','NMS',' MS'/
        Character(Len=256) :: strfmt, strfmt2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Cc(Nd), Dd(Nd), W(Nc,IPlv))
        ist=(Ksig+1)+3*Kbrt          !### stecp(ist) is Used for output
        If (K_is == 3) K_sms=4       !### Used for output
        If (Kecp == 1) ist=7
        Open(unit=16,file='CONF.XIJ',status='UNKNOWN',form='unformatted')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! printing eigenvalues in increasing order
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        strfmt = '(4X,82("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        If (Ksig*Kdsig == 0) Then
            strfmt = '(4X," Energy levels (",A7," Nc=",I7," Nd=",I9,"); &
                    Gj =",F7.4,/,4X,"  N",6X,"JTOT",9X,"L",9X,"S",12X, &
                    "EV",14X,"ET",7X,"DEL(CM**-1)")'
            Write( 6,strfmt) stecp(ist),Nc,Nd,Gj
            Write(11,strfmt) stecp(ist),Nc,Nd,Gj
        Else
            strfmt = '(4X,"Energy levels ",A7,", Sigma(E =", &
                 F10.4,") extrapolation var.",I2,/4X,"(Nc=",I6, &
                 " Nd=",I9,");  Gj =",F7.4,/4X,"N",6X,"JTOT",12X, &
                 "EV",16X,"ET",9X,"DEL(CM**-1)")'
            Write( 6,strfmt) stecp(ist),E_0,Kexn,Nc,Nd,Gj
            Write(11,strfmt) stecp(ist),E_0,Kexn,Nc,Nd,Gj
        End If

        If (C_is /= 0.d0) Then
            If (K_is == 1) Then
                strfmt = '(4X,"Volume shift: dR_N/R_N=",F9.5," Rnuc=",F10.7)'
                Write( *,strfmt) C_is,Rnuc
                Write(11,strfmt) C_is,Rnuc
            Else
                strfmt = '(4X,A3,":",E9.2,"*(P_i Dot P_k) ",A6, &
                     " Lower component key =",I2)'
                Write( *,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
                Write(11,strfmt) strms(K_is-1),C_is,strsms(K_sms),Klow
            End If
        End If

        strfmt = '(4X,82("-"))'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        jmax=min(Nlv,Nd)
        Do j=1,jmax
            Rewind(16)
            imax=j-1
            If (imax >= 1) Then
                Do i=1,imax
                    READ(16)
                End Do
            End If
            Read(16) Er(J),xj,idum,(CC(I),i=1,Nd)
            Er(j)=Er(j)+4.d0*Gj*xj*(xj+1.d0)
            E=Er(J)
            DT=E-Ecore
            ! Rydberg constant is taken from "phys.par"
            DEL=(ER(1)-ER(J))*2*DPRy
            strfmt = '(4X,I3,F14.9,2F10.5,F14.8,F15.6,F15.2)'
            Write( 6,strfmt) j,xj,Tl(j),Ts(j),E,DT,DEL
            Write(11,strfmt) j,xj,Tl(j),Ts(j),E,DT,DEL
        End Do

        strfmt = '(4X,82("="))'
        Write( 6,strfmt)
        Write(11,strfmt)

        ! weights of configurations
        Do j=1,jmax
            Rewind(16)
            imax=j-1
            If (IMAX >= 1) Then
                Do i=1,imax
                    Read(16)
                End Do
            End If
            Read(16) D,DUMMY,idum,(CC(I),i=1,Nd)
            i=0
            Do ic=1,Nc
                D=0.d0
                ndk=Ndc(ic)
                Do k=1,ndk
                    i=i+1
                    D=D+CC(i)**2
                End Do
                W(ic,j)=D
            End Do
        End Do

        n=(jmax-1)/5+1
        j2=0

        Do k=1,n
            nk=5
            If (k == n) nk=jmax-(n-1)*5
            j1=j2+1
            j2=j2+nk
            j3=j2+1
            strfmt2 = '(3X,66A1)'
            Write(11,strfmt2) (st1,i=j1,j3)
            strfmt = '(15X,5(3X,I2,6X))'
            Write(11,strfmt) (i,i=j1,j2)
            Write(11,strfmt2) (st2,i=j1,j3)
            strfmt = '(4X,"ICONF",4X,5F11.5)'
            Write(11,strfmt) (Er(i),i=j1,j2)
            Write(11,strfmt2) (st2,i=j1,j3)
            strfmt = '(2X,I6,"      ",5F11.6)'
            Do ic=1,Nc
                i=ic
                Write(11,strfmt) i,(W(i,j),j=j1,j2)
            End Do
            Write(11,strfmt2) (st1,i=j1,j3)
        End do
        Close(unit=16)
        Close(unit=6)
        Close(unit=11)
        Deallocate(Cc, Dd, W, Ndc, Tl, Ts)
        Return
    End Subroutine PrintLSJ
    
End Module lsj_aux