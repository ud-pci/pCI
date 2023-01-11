Module amp_ops
    ! This module contains functions to calculate amplitudes of one-electron operators
    ! One-electron operators include:
    !       AE1  = REDUCED E1 AMPLITUDE IN L GAUGE
    !       AE1V = REDUCED E1 AMPLITUDE IN V GAUGE
    !       A    = REDUCED Magnetic HFS AMPLITUDE
    !       G    = REDUCED M1 AMPLITUDE
    !       EDM  = ELECTRON DIPOLE AMPLITUDE
    !       PNC  = WEAK CHARGE AMPLITUDE
    !       AM   = ANAPOLE MOMENT AMPLITUDE
    !       QM   = Magnetic Quadrupole moment amplitude
    !       AE2  = REDUCED E2 AMPLITUDE
    !       AE3  = REDUCED E3 AMPLITUDE
    !       AM2  = REDUCED M2 AMPLITUDE
    !       AM3  = REDUCED M3 AMPLITUDE
    Use dtm_variables
    Implicit None

    Private

    Public :: AmpE1, AmpE2, AmpE3, AmpM1, AmpM2, AmpM3, AmpEDM, AmpPNC
    Public :: AmpAM, AmpMQM, HfsA, HfsB, Fint, AmpOut

  Contains 

    Real(dp) function Fint(is,nfin,nini,ic)       
        ! this function searches for radial integrals of one-electron operators
        Implicit None
        Integer :: isg, na, nb, is, nfin, nini, ic, ind, i
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        isg=1
        na=nfin
        nb=nini
        If (na > nb) Then
            na=nini
            nb=nfin
            isg=ic
        End If
        ind=is*IPx*IPx+(na-Nso)*IPx+(nb-Nso)
        Do i=1,Nint
            If (ind == Intg(i)) Then
                Fint=Rnt(i)*isg
                Return
            End If
        End Do
        Write( 6,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        Write(11,'(1X,"Fint: NO INTEGRAL ",A4,2I4,I8)') Alet(is),nfin,nini,ind
        Fint=0.d0
        Return
    End function Fint
    
    Subroutine AmpOut(i,x,nl,nk,xjl,xjk,ll,lk,y)
        Implicit None
        Integer :: i, nl, nk, ll, lk, nnl, nnk, jl, jk
        Real(dp) :: x, xjl, xjk, y
            If (dabs(x) < 1.d-3) Return
            nnl=Nf0(nl)
            nnk=Nf0(nk)
            nnl=Nh(nnl+1)
            nnk=Nh(nnk+1)
            nnl=Nn(nnl)
            nnk=Nn(nnk)
            jl=2*xjl+0.1d0
            jk=2*xjk+0.1d0
            Write(11,'(1X,A4,":",F8.4,I4,A1,I2,"/2  <<",I4,A1,I2,"/2",2E14.4)') &
             Alet(i),x,nnk,let(lk+1),jk,nnl,let(ll+1),jl,y
        Return
    End Subroutine AmpOut

    Real(dp) function AmpE1(Ro, nk,xjk,lk, nl,xjl,ll) 
        ! this function calculates the amplitude of the electric dipole
        ! transition matrix element <k|E1|l>
        Use wigner        
        Implicit None
        Integer :: nk, lk, nl, ll, is, k, lx
        Real(dp) :: Ro, xjk, xjl, AL, AV, xlk, xll, c
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        AL=Ro*Fint(3,nk,nl,+1)
        AV=Ro*Fint(6,nk,nl,-1)
        If (AL /= 0.d0 .or. AV /= 0.d0) Then
            is = 1
            lx=max(lk,ll)
            k=xjl+lx+1.51d0
            If (k /= 2*(k/2)) is=-is
            xlk=lk
            xll=ll
            c=dsqrt((2*xjk+1)*(2*xjl+1)*lx)*Fj6(xlk,xjk,0.5d0,xjl,xll,1.d0)
            AV=is*c*AV
            AL=is*c*AL
            If (Kl == 1) Then
                Call AmpOut(3,Ro,nl,nk,xjl,xjk,ll,lk,AL)
                Call AmpOut(6,Ro,nl,nk,xjl,xjk,ll,lk,AV)
            End If
        End If
        AmpE1=AL
        AE1V=AE1V+AV
        Return
    End function AmpE1

    Real(dp) function AmpE2(Ro, nk,xjk,lk, nl,xjl,ll)  
        ! this function calculates the amplitude of the electric quadrupole
        ! transition matrix element <k||E2||l>  
        Use wigner    
        Implicit None
        Integer :: nk, lk, nl, ll, is, k
        Real(dp) :: Ro, xjk, xjl, AE, xlk, xll
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        AE= Ro*Fint(10,nk,nl,+1)
        If (AE /= 0.d0) Then
            is = 1
            k=xjl+0.51d0
            If (k /= 2*(k/2)) is=-is
            xlk=lk
            xll=ll
            AE= AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
                  *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
                  *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
            If (Kl == 1) Call AmpOut(10,Ro,nl,nk,xjl,xjk,ll,lk,AE)
        End If
        AmpE2=AE
        Return
    End function AmpE2

    Real(dp) function AmpE3(Ro, nk,xjk,lk, nl,xjl,ll)    
        ! this function calculates the amplitude of the electric octupole
        ! transition matrix element <k||E3||l>    
        Use wigner
        Implicit None
        Integer :: nk, lk, nl, ll, is, k
        Real(dp) :: Ro, xjk, xjl, AE, xlk, xll
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        AE=Ro*Fint(11,nk,nl,+1)
        If (AE /= 0.d0) Then
            is = 1
            k=xjl+0.51d0
            If (k /= 2*(k/2)) is=-is
            xlk=lk
            xll=ll
            AE=AE*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
                *Fj3(xlk,xll,3.d0,0.d0,0.d0,0.d0) &
                *Fj6(xjk,xjl,3.d0,xll,xlk,0.5d0)
            If (Kl == 1) Call AmpOut(11,Ro,nl,nk,xjl,xjk,ll,lk,AE)
        End If
        AmpE3=AE
        Return
    End function AmpE3

    Real(dp) function AmpM1(Ro, nk,xjk,lk, nl,xjl,ll)          
        ! this function calculates the amplitude of the magnetic dipole
        ! transition matrix element <k|M1|l>
        Implicit None
        Integer :: nk, nl, k, is, ll, lk
        Real(dp) :: A, Ro, xjl, xjk, xjm, c
        A=Ro*Fint(9,nk,nl,+1)
        If (A /= 0.d0) Then
            is = 1
            k=xjl+ll+1.51d0
            If (k /= 2*(k/2)) is=-is
            xjm=xjl
            If (xjm > xjk+0.1d0) xjm=xjk
            c=0.d0
            If (dabs(xjl-xjk) > 0.1d0) &
                c= dsqrt((2*xjk+1)*(2*xjl+1)/(xjm+1))
            If (dabs(xjl-xjk) <= 0.1d0) &
                c= (2*xjk+1)*dsqrt((2*xjk+1)/(xjk*(xjk+1)))
            A=A*is*c
            If (Kl == 1) Call AmpOut(9,Ro,nl,nk,xjl,xjk,ll,lk,A)
        End If
        AmpM1=A
        Return
    End function AmpM1

    Real(dp) function AmpM2(Ro, nk,xjk,lk, nl,xjl,ll)     
      ! this function calculates the amplitude of the magnetic quadrupole
      ! transition matrix <element k||M2||l>
      Use wigner   
      Implicit None
      Integer :: nk, lk, nl, ll, is, k
      Real(dp) :: Ro, xjk, xjl, AE
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      AE=Ro*Fint(12,nk,nl,+1)
      If (AE /= 0.d0) Then
            is = 1
            k = xjk+0.51d0
            If (k /= 2*(k/2)) is=-is
            AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                    Fj3(xjk,xjl,2.d0,-0.5d0,0.5d0,0.d0)
            If (Kl == 1) Call AmpOut(12,Ro,nl,nk,xjl,xjk,ll,lk,AE)
      End If
      AmpM2=AE
      Return
    End function AmpM2

    Real(dp) function AmpM3(Ro, nk,xjk,lk, nl,xjl,ll)   
        ! this function calculates the amplitude of the magnetic octupole
        ! transition matrix element <k||M3||l>
        Use wigner     
        Implicit None
        Integer :: nk, lk, nl, ll, is, k
        Real(dp) :: Ro, xjk, xjl, AE
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        AE=Ro*Fint(13,nk,nl,+1)
        If (AE /= 0.d0) Then
            is = 1
            k = xjk+0.51d0
            If (k /= 2*(k/2)) is=-is
            AE= AE* is* dsqrt((2*xjk+1)*(2*xjl+1))* &
                    Fj3(xjk,xjl,3.d0,-0.5d0,0.5d0,0.d0)
            If (Kl == 1) Call AmpOut(13,Ro,nl,nk,xjl,xjk,ll,lk,AE)
        End If
        AmpM3=AE
        Return
    End function AmpM3

    Real(dp) function AmpEDM(Ro, nk,xjk,lk, nl,xjl,ll)        
        ! this function calculates the amplitude of the P, T-odd interaction 
        ! of the electron electric dipole moment <k|D|l>
        Implicit None
        Integer :: nk, lk, nl, ll
        Real(dp) :: Ro, xjk, xjl, A
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        A=-2*Ro*Fint(4,nk,nl,+1)
        If (A /= 0.d0 .and. Kl == 1) Call AmpOut(4,Ro,nl,nk,xjl,xjk,ll,lk,A)
        AmpEDM=A
        Return
    End function AmpEDM

    Real(dp) function AmpPNC(Ro, nk,xjk,lk, nl,xjl,ll)        
        ! this function calculates the nuclear spin independent 
        ! parity nonconserving (PNC) amplitude <k|W|l>
        Implicit None
        Integer :: nk, lk, nl, ll
        Real(dp) :: Ro, xjk, xjl, A
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        A=-Ro*Fint(5,nk,nl,-1)
        If (A /= 0.d0 .and. Kl == 1) Call AmpOut(5,Ro,nl,nk,xjl,xjk,ll,lk,A)
        AmpPNC=A
        Return
    End function AmpPNC

    Real(dp) function AmpAM(Ro, nk,xjk,lk, nl,xjl,ll)         
        ! this function calculates the amplitude of the electron interaction
        ! with the P-odd nuclear anapole moment <k|Am|l>
        Implicit None
        Integer :: nk, lk, nl, ll, is, k, lx
        Real(dp) :: Ro, xjk, xjl, A
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        A=Ro*Fint(7,nk,nl,-1)
        lx=max(ll,lk)
        k=xjl+0.51d0+lx
        is=1
        If (k /= 2*(k/2)) is=-is
        A=is*A
        If (A /= 0.d0 .and. Kl == 1) Call AmpOut(7,Ro,nl,nk,xjl,xjk,ll,lk,A)
        AmpAM=A
        Return
    End function AmpAM

    Real(dp) function AmpMQM(Ro, nk,xjk,lk, nl,xjl,ll)        
        ! this function calculates the amplitude of the nucleus 
        ! magnetic quadrupole moment <k|MQM|l>
        Use wigner
        Implicit None
        Integer :: nk, lk, nl, ll, k, is1, is2
        Real(dp) :: Ro, xjk, xjl, B, tll, xlk, xll, g
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        B=Ro*Fint(8,nk,nl,+1)
        If (B /= 0.d0) Then
            B=1.5d0*B*dsqrt((2*xjk+1)*(2*xjl+1))
            xlk=lk
            xll=ll
            g=xlk
            If (g < xll) g=xll
            tll=2*xjl-xll
            is1 = 1
            k=xlk+g+0.1d0
            If (k /= 2*(k/2)) is1=-is1
            is2 = 1
            k=xjl+0.51d0
            If (k /= 2*(k/2)) is2=-is2
            B = B &
              * ( is1*dsqrt(30*g) &
              * FJ9(xll,xlk,1.d0,xjl,xjk,2.d0,0.5d0,0.5d0,1.d0) &
              + is2*dsqrt(2.d0/3.d0*(2*xlk+1)*(2*tll+1)) &
              * Fj3(xlk,2.d0,tll,0.d0,0.d0,0.d0) &
              * Fj6(xlk,xjk,0.5d0,xjl,tll,2.d0) )
            If (Kl == 1) Call AmpOut(8,Ro,nl,nk,xjl,xjk,ll,lk,B)
        End If
        AmpMQM=B
        Return
    End function AmpMQM

    Real(dp) function HfsA(Ro, nk,xjk,lk, nl,xjl,ll)          
        ! this function calculates magnetic quadrupole hyperfine structure constant <k|A|l>
        Implicit None
        Integer :: is, k, nk, lk, nl, ll
        Real(dp) :: Ro, xjk, xjl, A, c, xjm
        A=Ro*Fint(1,nk,nl,+1)
        If (A /= 0.d0) Then
            is = 1
            k=xjl+ll+1.51d0
            If (k /= 2*(k/2)) is=-is
            xjm=xjl
            If (xjm > xjk+0.1d0) xjm=xjk
            c=0.d0
            If (dabs(xjl-xjk) > 0.1d0) c=dsqrt((2*xjk+1)*(2*xjl+1)/(xjm+1))
            If (dabs(xjl-xjk) <= 0.1d0) c=(2*xjk+1)*dsqrt((2*xjk+1)/(xjk*(xjk+1)))
            A=A*is*c
            If (Kl == 1) Call AmpOut(1,Ro,nl,nk,xjl,xjk,ll,lk,A)
        End If
        HfsA=A
        Return
    End function HfsA

    Real(dp) function HfsB(Ro, nk,xjk,lk, nl,xjl,ll)   
        ! this function calculates electric quadrupole hyperfine structure constant <k|B|l>
        Use wigner    
        Implicit None
        Integer :: nk, nl, lk, ll, is, k
        Real(dp) :: Ro, B, xjk, xjl, xlk, xll
            B=Ro*Fint(2,nk,nl,+1)
            If (B /= 0.d0) Then
              is = 1
              k=xjl+0.51d0
              If (k /= 2*(k/2)) is=-is
              xlk=lk
              xll=ll
              B=B*is*dsqrt((2*xjk+1)*(2*xjl+1)*(2*xlk+1)*(2*xll+1)) &
                 *Fj3(xlk,xll,2.d0,0.d0,0.d0,0.d0) &
                 *Fj6(xjk,xjl,2.d0,xll,xlk,0.5d0)
              If (Kl == 1) Call AmpOut(2,Ro,nl,nk,xjl,xjk,ll,lk,B)
            End If
            HfsB=B
        Return
    End function HfsB

End Module amp_ops