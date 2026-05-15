Module bspl

    Use params

    Implicit None

    Private

    Public :: FormBspl, Grid_Bspl

    Integer, Parameter :: IPisp = 1200          ! number of nodes on the grid
    Integer, Parameter :: IPmsp = 70            ! number of B-splines
    Real(dp), Dimension(IPisp, IPmsp) :: Bsp      ! array of splines
    Real(dp), Dimension(IPisp) :: Rsp            ! radial grid for B-splines

Contains

    Subroutine FormBspl(ks, ms, m1, mx, nount, C, ii)
        ! ks - rank of splines
        ! ms - number of splines
        ! m1 - grid multiplier
        ! mx - No of nodes
        ! nout - spline to form
        ! C - formed spline
        Implicit None

        Integer :: ks, ms, m1, mx, nount, i, i1, ii, n, k, mk, mm, im, nout, mn
        Real(dp) :: r1, r2, r3, r4, rm, b1
        Real(dp), Dimension(IP6) :: C

        if (ms.GT.IPmsp.OR.mx.GT.IPisp) then
            write(*,*) ' Bspl error:'
            write(*,*) ' No more than',IPmsp,' B-splines on',IPisp,' grid points'
            write(*,*) ' Here Nsp=',ms,' mx=',mx
            stop
        end if

        do n=1,ms
            do m=1,mx
                Bsp(m,n)=0.d0
            end do
        end do

        do n=1,ms
            mn=m1*(n-1)
            do m=1,m1
                Bsp(mn+m,n)=1.d0
            end do
        end do

        do k=2,ks
            mk=m1*k
            do n=1,ms
                mn=m1*(n-1)+1
                r1=Rsp(mn)
                r2=Rsp(mn+m1)
                r4=Rsp(mn+mk)
                r3=Rsp(mn+mk-m1)
                do m=1,mk
                    mm=mn+m-1
                    rm=Rsp(mm)
                    b1=0.d0
                    if (m.LE.(mk-m1)) b1=b1+(rm-r1)/(r3-r1)*Bsp(mm,n)
                    if (m.GT.m1.AND.n.LT.ms) b1=b1+(r4-rm)/(r4-r2)*Bsp(mm,n+1)
                    Bsp(mm,n)=b1
                end do
            end do
        end do

        im=m1*ks
        i1=Ii-1
        ms=ms-ks           ! initial number of splines is restored
        do n=1,ms
            do m=1,mx
                if (m.LE.i1) then
                    Bsp(m,n)=Bsp(m+im,n+1)
                else
                    Bsp(m,n)=0.d0
                end if
            end do
        end do

        do i=1,Ii
            C(i)=Bsp(i,nout)
        end do
        C(Ii+4)=9.d0

        Return
    End Subroutine FormBspl

    Subroutine Grid_Bspl(ksp, msp, m1, mx, R, ii)
        ! ksp - rank of splines
        ! msp - number of splines
        ! m1,mx - Forms grid for B-splines, m1 - grid multiplier, mx - No of nodes
        ! R - original radial grid
        ! ii - original number of grid points
        Implicit None

        Integer :: ksp, i1, m1, msp, im, mx, i, ii
        Real(dp) :: small, dl, dr, rl, rr
        Real(dp), Dimension(IP6) :: R

        small=1.d-6               ! small factor is used on grid ends
        i1=Ii-1
        m1=i1/msp
        i1=m1*msp                 ! last (Ii-i1) nodes are skipped!
   
        msp=msp+Ksp               ! Ksp extra splines used in construction
        dl=(R(2)-R(1))*small      ! step for extra nodes on the left
        dr=(R(i1+1)-R(i1))*small  ! step for extra nodes on the right
        im=m1*Ksp
        mx=i1+2*im                ! Number of grid nodes
        rl=R(1)
        rr=R(i1+1)
   
        do i=1,im
            Rsp(i)=rl-dl*(im-i-1)
            Rsp(im+i1+i)=rr+dr*(i-1)
        end do
   
        do i=1,i1
          Rsp(i+im)=R(i)
        end do
   
        Return
    End Subroutine Grid_Bspl

End Module bspl