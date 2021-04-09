Module test_ori
    ! This module contains a subroutine to test Taylor expansion for
    ! orbitals and functions in HFD.DAT format
    Use basc_variables

    Implicit None
    
    Private

    Public :: Test_Origin

  Contains

    Subroutine Test_Origin(str,P,err,ierr)  
        ! Test of the Taylor expansion at the origin
        Implicit None

        Integer :: ierr, i, k, ik, igm
        Real(dp) :: gm, r1, r2, q, q0, pp1, pp2, err1, err2, err
        Real(dp), Dimension (IP6) :: P
        Character(Len=1) :: str(10), chr
        Character(Len=512) :: strfmt

        ierr=0
        gm=P(ii+4)
        r1=R(1)
        r2=R(2)
        q=r2/r1
        q0=1.d0
        pp1=0.d0                      ! Expansion value for P(1)
        pp2=0.d0                      ! Expansion value for P(2)
        Do k=0,MaxT
            ik=ii+5+k
            pp1=pp1+P(ik)
            pp2=pp2+q0*P(ik)
            q0=q0*q
        End Do

        chr='.'
        If (pp1*P(1).GT.0.d0) Then     ! Leading power is not always given by P(ii+4).
            If (dabs(gm).LT.1.d-3.AND.dabs(pp1).GT.1.d-9) Then ! In this case it is determined
                gm=dlog(P(1)/pp1)/dlog(r1)                     ! here, which is marked by
                chr='!'                                        ! exclamation mark in the output
                igm=gm+0.1d0                                   !### This assumes finite nucleus,
                gm=igm                                         !### when gamma is always integer
            End If
        End If

        pp1=r1**gm*pp1
        pp2=r2**gm*pp2
        err1=P(1)-pp1
        err2=P(2)-pp2
        If(dabs(pp1).GT.1.d-2) Then
            err1=err1/(P(1)+pp1)
            err2=err2/(P(2)+pp2)
        End If

        If (dabs(err1).GT.err) Then
            strfmt='(4X,"Error in Taylor expansion at the origin for ",10A1,":")'
            If (Kout.GT.0)  Write(11,strfmt) str
            If (Kout.GT.1)  Write( *,strfmt) str
            strfmt='(2X,"Input:  gam=",F8.5,2(4X,"P(",I1,")=",E13.5),/2(4X,"C_",I1," = ",E11.4))'
            If (Kout.GT.1)  Write(11,strfmt) P(ii+4),(i,P(i),i=1,2),(m,P(ii+5+m),m=0,MaxT)
            strfmt='(2X,"Output:  gm=",F8.5,A1," Relative error at r1 and r2",2E13.5)'
            If (Kout.GT.0)  Write(11,strfmt) gm,chr,err1,err2
            If (Kout.GT.1)  Write( *,strfmt) gm,chr,err1,err2
            If (Kout.GT.1)  read(*,*)
            ierr=dabs(err1/err)
        End If

        P(ii+4)=gm

        Return
    End subroutine Test_Origin

End Module test_ori