Module pi_pk

    Use basc_variables
    
    Implicit None

    Private

    Public :: SMS_core, V_nms, P_eff

  Contains

    Real(dp) Function SMS_core(na,nb,kan)  
        !### exchange core potential for (pi_pk)
        Use wigner
        Use readfff
        Implicit None

        Integer :: ja, jb, la, lb, na, nb, kan, jc, lc
        Real(dp) :: a_ab, z12, z00, f0, s1
        Character(Len=1), Dimension(9) :: let 
        Character(Len=512) :: strfmt
        data let /'s','p','d','f','g','h','i','k','l'/

        a_ab=0.d0
        If (Nso.EQ.0) Then
            SMS_core = a_ab
            Return
        End If

        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)


        If (ja.NE.jb.OR.la.NE.lb) Then
            SMS_core = a_ab
            Return
        End If

        z12=0.5d0
        z00=0.0d0
        xja=z12*ja
        Do nc=1,Nso
            jc=Jj(nc)
            lc=Ll(nc)
            If (iabs(ja-jc).GT.2.OR.iabs(la-lc).NE.1) Cycle
            xjc=z12*jc
            f0=(jc+1)*FJ3(xja,xjc,1.d0,z12,-z12,z00)**2
            Call ReadF(kan,nc+4,Pc,Qc,2)
            Call ReadF(kan,nb+4,Pa,Qa,2)
            s1=+f0*P_eff(jc,lc,jb,lb,Pc,Qc,Pa,Qa)
            Call ReadF(kan,na+4,Pa,Qa,2)
            s1=-s1*P_eff(jc,lc,ja,la,Pc,Qc,Pa,Qa)
            If (kout.GE.2) write(*,'(4X,"contribution of ",I2,A1,I1,"/2 =",E12.5)') &
                              Nn(nc),let(lc+1),jc,s1
            a_ab=a_ab+s1
        End Do

        If (kout.GE.1) Then
            strfmt='(/2X,"<",I2,A1,I2,"/2| V_SMS^core |",I2,A1,I2,"/2> = ",E12.5)'
            write( *,strfmt) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
            write(11,strfmt) Nn(na),let(la+1),ja,Nn(nb),let(lb+1),jb,a_ab
        End If

        SMS_core = a_ab

        Return
    End Function SMS_core

    Real(dp) Function P_eff(ja,la,jc,lc,Pa,Qa,Pc,Qc) 
        !### effective vertex for (pi dot pk)
        Use diff
        Use wigner
        Use sintg
        Implicit None

        Integer :: la, lc, ila, ja, ilc, jc, lmax, ilmax, ih, i
        Real(dp) :: ga, half, xjx, xjn, s0, s1, ss, s2, ds
        Real(dp), Dimension(IP6) :: Pa,Pc,Qa,Qc,Da

        P_eff=0.d0
        If (iabs(la-lc).NE.1) Return

        ila=ja-la ! l for low component
        ilc=jc-lc ! "      "      "
        lmax=max(la,lc)
        ilmax=max(ila,ilc)

        ih=2-Kt
        C(ii:ih)=0.d0

        ga=Pa(ii+4)               !### upper component
        C(ii+4)=Pc(ii+4)+ga-1
        Call Dif(Pa,Da,R,V,ga,ii,kt,MaxT,h)

        Do i=1,Ii,ih
            C(i)=Pc(i)*(Da(i)+(la-lc)*lmax/R(i)*Pa(i))
        End Do

        If (klow.GT.0) Then       !### lower component
            Call Dif(Qa,Da,R,V,ga,ii,kt,MaxT,h)
            Do i=1,Ii,ih
                C(i)=C(i)+Qc(i)*(Da(i)+(ila-ilc)*ilmax/R(i)*Qa(i))
            End Do
        End If

        If (klow.GT.1) Then       !### relativistic correction
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
    
            If (dabs(s0).GT.1.d9) Then
                write(*,*) ' s0 = infty'
                write(*,*) ja,jc
                write(*,*) la,lc
                read(*,*)
            End If
    
            If (mod(la+1,2).NE.0) s0=-s0
            s1=0.d0
            If (lc.EQ.ila) s1=dsqrt((xjc+3*xja-2*ila+1)/(2*ila+1))
            ss=0.5d0*Z/DPcl*(s0*s1+1.0d0)
            Do i=1,Ii,ih
                C(i)=C(i)+ss/R(i)*Pc(i)*Qa(i)
            End Do
    
            s2=0.d0
            If (ilc.EQ.la) s2=dsqrt((xjc+3*xja-2*la+1)/(2*la+1))
            ss=0.5d0*Z/DPcl*(s0*s2-1.0d0)
            Do i=1,Ii,ih
                C(i)=C(i)+ss/R(i)*Qc(i)*Pa(i)
            End Do
        End If

        Call Sint1(ds)
        P_eff=ds

        Return
    End Function P_eff

    Real(dp) Function V_nms(na,nb,kan)  
        !### Normal mass shift 
        Use diff
        Use readfff
        Use sintg
        Implicit None

        Integer :: ja, jb, la, lb, na, nb, kan, ih, ip, ila, iq, i
        Real(dp) :: ga, cis2, clj, az, ri, ri2, xi, a_ab

        If (k_is.LE.2) Then ! NMS is calculated for k_is=3,4
            V_nms = a_ab
            Return
        End If

        ja=Jj(na)
        jb=Jj(nb)
        la=Ll(na)
        lb=Ll(nb)

        If (ja.NE.jb.OR.la.NE.lb) Then
            V_nms = a_ab
            Return
        End If

        ih=2-kt
        Call ReadF(kan,na+4,Pa,Qa,2)
        Call ReadF(kan,nb+4,Pb,Qb,2)
        ga=Pa(ii+4)
        Call Dif(Pa,P1a,R,V,ga,ii,kt,MaxT,h)
        Call Dif(Pb,P1b,R,V,ga,ii,kt,MaxT,h)
        ip=la*(la+1)
        cis2=0.5d0                !  No c_is here
        If (klow.GT.0) Then
            Call Dif(Qa,Q1a,R,V,ga,ii,kt,MaxT,h)
            Call Dif(Qb,Q1b,R,V,ga,ii,kt,MaxT,h)
            ila=ja-la
            iq=ila*(ila+1)
            xja=0.5d0*ja
            clj=-xja*(xja+1)+la*(la+1)+0.75d0 - 2.d0
            aZ=Z/DPcl
        End If

        Do i=1,ii,ih
            ri=R(i)
            ri2=ri*ri
            C(i)=P1b(i)*P1a(i)+ip*Pb(i)*Pa(i)/ri2
            If (klow.GT.0) Then
                C(i)=C(i)+Q1b(i)*Q1a(i)+iq*Qb(i)*Qa(i)/ri2
            End If
            If (klow.GE.2) Then
                xi = 2/ri*(P1b(i)*Qa(i)+Qb(i)*P1a(i)) &
                    + clj/ri2*(Pb(i)*Qa(i)+ Qb(i)*Pa(i))
                C(i)=C(i)+aZ*xi
            End If
            C(i)=cis2*C(i)
        End Do
        C(ii+4)=2*Pa(ii+4)-2
        Call Sint1(a_ab)

        V_nms = a_ab

        Return
    End Function V_nms
    
End Module pi_pk