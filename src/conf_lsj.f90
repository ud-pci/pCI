Module conf_lsj
    ! subroutines lsj, lsj_det, pls, pll, p0s, p0l
    Use conf_variables

    Implicit None

  Contains

    Subroutine lsj(cc,nmax,xj,xl,xs)
        Implicit None
        Real(dp), Dimension(nmax) :: cc(nmax)
        Integer, Allocatable, Dimension(:) :: idet1, idet2

        xj=0.d0
        xl=0.d0
        xs=0.d0
        n=0

        Do ic1=1,Nc
            ndn=Ndc(ic1)
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

        Return
    End Subroutine lsj


    Subroutine lsj_det(idet1,idet2,tj,tl,ts)
        Implicit None
        Integer, Allocatable, Dimension(:) :: idet1, idet2

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

    Function pls(ia,ib)
        Implicit None
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

    Function pll(ia,ib)
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

    Function p0s(ia,ib)
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
        Do 10 msa=-1,1,2
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

    Function p0l(ia,ib)
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

End Module conf_lsj