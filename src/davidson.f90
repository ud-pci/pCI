Module davidson
    ! This module contains main subroutines required for the Davidson procedure.

    Use conf_variables

    Implicit None

    Private

    Public :: Init4, Hould, FormB0, FormB, FormBskip, FormP, Dvdsn, Mxmpy, Ortn, Prj_J, Vread, VWrite

  Contains
    
    Subroutine Init4(mype,npes)
        ! This subroutine constructs the initial approximation 
        ! by selecting configurations specified by parameter Nc4.
        ! The initial approximation Hamil is stored in the 
        ! matrix Z1 and is constructed by selecting the top-left 
        ! block of the full Hamil matrix H.
        ! The diagonal elements of H are stored in the array Diag.
        Use mpi
        Use determinants, Only : calcNd0
        Implicit None
        Integer  :: k, ierr, n, n1, n2, ic, ic1, mype, npes, mpierr
        Integer(Kind=int64) :: l8
        Real(dp) :: t

        ! Dimension of the approximation to start with:
        Call calcNd0(ic1, Nd0)
        If (mype == 0) Then
            Write( 6,'(3X,"Starting approx. includes ",I3," conf.,",I4," det.")') ic1,Nd0
            Write(11,'(3X,"Starting approx. includes ",I3," conf.,",I4," det.")') ic1,Nd0
        End If

        ! Construct the Hamil in the initial approximation:
        Diag=0.d0
        Z1=0.d0
        Do l8=1,ih8H
            n=Hamil%n(l8)
            k=Hamil%k(l8)
            t=Hamil%t(l8)
            If (n == k) t=t-Hmin
            If (n <= Nd0) Then
                Z1(n,k)=t
                Z1(k,n)=t
            Else
                Exit
            End If
        End Do
        If (mype==0) Then
            Call MPI_Reduce(MPI_IN_PLACE, Z1, Nd0*Nd0, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                MPI_COMM_WORLD, mpierr)
          Else
            Call MPI_Reduce(Z1, Z1, Nd0*Nd0, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                MPI_COMM_WORLD, mpierr)
        End If
    End Subroutine Init4

    Subroutine Hould(n,ndim,Ee,Dd,Zz)
        ! This subroutine diagonalizes the matrix Zz using 
        ! Householder's method of diagonalization.
        ! Dd-array of eigenvalues
        ! Zz-matrix of eigenvectors
        ! n-dimension of matrix z
        ! eps-criteria of diagonalization
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Implicit None

        Integer :: n, ii, k, j, l, i, jj, m1, ndim, im, im1, li, k1, l1
        Real(dp) :: tol, f, g, b, r, em, h, hh, c, ei, s, p, di, eps
        Real(dp), dimension(:), allocatable :: Ee, Dd
        Real(dp), dimension(:,:), allocatable :: Zz

        Ifail=0
        tol=2.0d0**(-1021) !-103 or -1021
        eps=2.0d0**(-53) ! -24 or -53
        If (n > 1) Then
            Do ii=2,n
                i=n-ii+2
                l=i-2
                f=Zz(i,i-1)
                g=0.0d0
                If (l > 0) Then
                    Do k=1,l
                        g=g+Zz(i,k)**2
                    End Do
                End If
                h=f**2+g
                If (g < tol) Then
                    Ee(i)=f
                    h=0.0d0
                    Dd(i)=h
                    Cycle
                End If
                l=l+1
                g=-dsign(dsqrt(h),f)
                Ee(i)=g
                h=h-f*g
                Zz(i,i-1)=f-g
                f=0.0d0
                Do j=1,l
                    Zz(j,i)=Zz(i,j)/h
                    g=0.0d0
                    Do k=1,j
                          g=g+Zz(j,k)*Zz(i,k)
                        End Do
                    k1=j+1
                    If (k1 <= l) Then
                        Do k=k1,l
                            g=g+Zz(k,j)*Zz(i,k)
                        End Do
                    End If
                    ee(j)=g/h
                    f=f+g*Zz(j,i)
                End Do
                hh=f/(h+h)
                Do j=1,l
                    f=Zz(i,j)
                    g=Ee(j)-hh*f
                    Ee(j)=g
                    Zz(j,1:j)=Zz(j,1:j)-f*Ee(1:j)-g*Zz(i,1:j)
                End Do
                Dd(i)=h
            End Do
        End If

        Dd(1)=0.0d0
        Ee(1)=0.0d0
        Do i=1,n
            l=i-1
            If (Dd(i) /= 0.0d0 .and. l /= 0) Then
                 Do j=1,l
                     g=0.0d0
                     Do k=1,l
                         g=g+Zz(i,k)*Zz(k,j)
                     End Do
                     Zz(1:l,j)=Zz(1:l,j)-g*Zz(1:l,i)
                 End Do
            End If
            Dd(i)=Zz(i,i)
            Zz(i,i)=1.0d0
            If (l /= 0) Then
                Zz(i,1:l)=0.0d0
                Zz(1:l,i)=0.0d0
            End If
         End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! p matrix is formed
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (n /= 1) Then
            Ee(1:n-1)=Ee(2:n)
        End If
        Ee(n)=0.0d0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! triadigonalisation is finished
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        b=0.0d0
        f=0.0d0
        Do l=1,n
            jj=0
            h=(dabs(Dd(l))+dabs(Ee(l)))*eps
            If (b < h) b=h
            Do m1=l,n
                m=m1
                em=dabs(Ee(m))
                If (em <= b) Exit
            End Do
            If (m /= l) Then
                Do while (dabs(Ee(l)) > b)
                    If (jj == 30) Then
                        ifail=1
                        Return
                    End If
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
                    If (im1 >= l) Then
                        Do im=l,im1
                            i=im1+l-im
                            ei=Ee(i)
                            g=c*ei
                            h=c*p
                            If (dabs(p) >= dabs(ei)) Then
                                c=ei/p
                                r=dsqrt(c**2+1.0d0)
                                Ee(i+1)=s*p*r
                                s=c/r
                                c=1.0d0/r
                            Else
                                c=p/ei
                                r=dsqrt(c**2+1.0d0)
                                Ee(i+1)=s*ei*r
                                s=1.0d0/r
                                c=c/r
                            End If
                            di=Dd(i)
                            p=c*di-s*g
                            Dd(i+1)=h+s*(c*g+s*di)
                            Do k=1,n
                                h=Zz(k,i+1)
                                ei=Zz(k,i)
                                Zz(k,i+1)=s*ei+c*h
                                Zz(k,i)=c*ei-s*h
                            End Do
                        End Do
                    End If
                    Ee(l)=s*p
                    Dd(l)=c*p
                End Do
            End If
            Dd(l)=Dd(l)+f
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! eigenvalues in D; eigenvectors in Z
        ! ordering of D and Z
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i=0
        Do while (i < n)
            i=i+1
            If (i >= n) Return
            h=Dd(i)
            f=Dd(i+1)
            If (h <= f) Cycle
            Dd(i)=f
            Dd(i+1)=h
            Do k=1,n
                h=Zz(k,i)
                Zz(k,i)=Zz(k,i+1)
                Zz(k,i+1)=h
            End Do
            If (i /= 1) i=i-2
        End Do
        Return
    End Subroutine Hould

    Subroutine FormB0(mype,npes)
        ! This subroutine constructs the initial eigenvectors for the Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are written to CONF.XIJ 
        Use mpi
        Use formj2, Only : J_av
        Implicit None
        Integer :: i, j, idum, ndpt, ierr, num, nskip, mpierr, mype, npes
        Real(dp) :: dummy, xj
        Character(Len=128) :: strfmt

        num=0
        nskip=0
        
        If (mype==0) Open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')

        Call MPI_Bcast(Z1,Nd0*Nd0,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)

        If (abs(Kl4) /= 2) Then ! If not reading CONF.XIJ
            Do j=1,Nd0
                Do i=1,Nd
                    If (j == 1)   B1(i)=0.d0
                    If (i <= Nd0) B1(i)=Z1(i,j)
                End Do
                Call J_av(B1,Nd0,xj,ierr,mype,npes)
                If (ierr == 0) Then
                    num=num+1
                    E(num)=-(E1(j)+Hmin)
                    Tj(num)=xj
                    ArrB(1:Nd,num)=B1(1:Nd)
                    If (num >= Nlv) goto 220
                Else
                    nskip=nskip+1
                    If (mype==0) Write(*,*) '  skip Ej =',-(E1(j)+Hmin)
                End If
            End Do
        Else ! Reading CONF.XIJ
            If (mype == 0) Then
                B1=0.d0
                Read(17,err=210) dummy,dummy,ndpt
 210            Write(*,*) ' Vector length = ',ndpt
                Rewind(17)
                Do J=1,Nlv
                    Read(17,end=220,err=220) E(J),Tj(J),idum,(B1(i),i=1,ndpt)
                    num=j
                    ArrB(1:Nd,j)=B1(1:Nd)
                End Do
            End If
        End If
 220    if (mype==0) Rewind(17)
        Nlv=num
        If (mype == 0) Then
            Do j=1,Nlv
                strfmt = '(1X,"E(",I2,") =",F12.6," Jtot =",F10.6)'
                Write( 6,strfmt) j,E(j),Tj(j)
                Write(11,strfmt) j,E(j),Tj(j)
                If (kXIJ > 0) Write(17) E(j),Tj(j),Nd,(ArrB(i,j),i=1,Nd)
            End Do
            close (unit=17)
            If (nskip > 0) Then
                Write(*,*) ' Selection is done', nskip,' vectors skiped'
                Write(11,*) ' Selection is done', nskip,' vectors skiped'
            End If
        End If
        Return
    End Subroutine FormB0

    Subroutine FormBskip
        ! This subroutine forms eigenvectors for next Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are not written to CONF.XIJ.
        Implicit None
        Integer :: i, idum, k, l
        Real(dp) :: t, tjk, ek

        ! Eigenvectors obtained by Davidson procedure are not written to CONF.XIJ file
        Do k=1,Nlv
            B2=0.d0
            Do l=1,2*Nlv
                t=P(l,k)
                B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
            End Do
            t=-(E(k)+Hmin)
            Tk(k)=t
            ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        End Do
        ArrB(1:Nd,1:Nlv)=ArrB(1:Nd,2*Nlv+1:3*Nlv)
        Write ( 6,*) '(1X,"FormBskip: Vectors not saved")'
        Return
    End Subroutine FormBskip

    Subroutine FormB
        ! This subroutine forms eigenvectors for next Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are written to CONF.XIJ.
        Implicit None
        Integer :: i, idum, k, l
        Real(dp) :: t, tjk, ek
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Eigenvectors obtained by Davidson procedure are written to CONF.XIJ file
        Open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        Do k=1,Nlv
            B2=0.d0
            Do l=1,2*Nlv
                t=P(l,k)
                B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
            End Do
            t=-(E(k)+Hmin)
            Tk(k)=t
            Write(17) t,Tj(k),Nd,(B2(i),i=1,Nd)
            ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        End Do
        ArrB(1:Nd,1:Nlv)=ArrB(1:Nd,2*Nlv+1:3*Nlv)
        close (unit=17)
        Write ( 6,*) '(1X,"FormB: Vectors saved")'
        Return
    End Subroutine FormB

    Subroutine FormP (kp)
        ! Calculation of upper-left block (Nlv X Nlv) of matrix P for kp=1 
        ! Calculation of other three blocks (2Nlv X 2Nlv) of matrix P for kp=2
        Implicit None
        Integer :: ix, i, kx, kp, k1, k2, j, k, i1
        Real(dp) :: x

        If (kp == 1) Then
            Do i=1,Nlv
              Do k=1,i
                k1=k+Nlv
                x=0.d0
                Do j=1,Nd
                    x=x+ArrB(j,i)*ArrB(j,k1)
                End Do
                P(i,k)=x
                P(k,i)=x
              End Do
            End Do
        Else                  ! kp=2
            Do i=1,Nlv
                kx=2*Nlv
                Do k=Nlv+1,kx
                    k1=k+Nlv
                    x=0.d0
                    Do j=1,Nd
                        x=x+ArrB(j,i)*ArrB(j,k1)
                    End Do
                    P(i,k)=x
                    P(k,i)=x
                End Do
            End Do
            Do i=Nlv+1,2*Nlv
                Do k=Nlv+1,i
                    k1=k+Nlv
                    x=0.d0
                    Do j=1,Nd
                        x=x+ArrB(j,i)*ArrB(j,k1)
                    End Do
                    P(i,k)=x
                    P(k,i)=x
                End Do
            End Do
            Do i=1,Nlv                   ! Cleanup for converged levels
              If (Iconverge(i) == 1) Then
                i1=i+Nlv
                Do k=1,2*Nlv
                    P(i1,k)=0.d0
                    P(k,i1)=0.d0
                End Do
                P(i1,i1)=2*vmax-P(i,i)+1.d0
              End If
            End Do
        End If
        Return
    End Subroutine FormP
    
    Subroutine Dvdsn (j)
        ! Main part of the Davidson algorithm. 
        ! Formation of the J-th residual vector and the corresponding new probe vector.
        ! New probe vectors are stored in the second block of ArrB.
        Implicit None
        Integer :: j1, i, j, l, k, n01
        Real(dp) :: val, t, cnorm, s
        character(len=9) :: char
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        cnorm=0.d0
        j1=j+Nlv
        val=P(j,j)
        Do i=1,Nd
            t=ArrB(i,j1)-val*ArrB(i,j)
            If (i <= Nd0) C(i)=t
            cnorm=cnorm+t*t
            If (i > Nd0) Then
                B1(i)=t/(val-Diag(i))
                Cycle
            End If
            B1(i)=0.d0
            Cycle
        End Do
        cnorm=dsqrt(dabs(cnorm))
        Iconverge(j)=0
        char='         '
        If (cnorm < Crt4) Then
            If (kdavidson == 1) Then
                char='Stopped  '
                Iconverge(j)=1
            Else
                char='converged'
            End If
        End If
        Write ( 6,'(4X,"|| C(",I2,") || = ",F10.7,2X,A9)') j,cnorm,char
        Write (11,'(4X,"|| C(",I2,") || = ",F10.7,2X,A9)') j,cnorm,char
        If (cnorm > cnx) cnx=cnorm
        If (Iconverge(j)==1) Then
            B1(1:Nd)=0.d0
            Return
        End If
        n01=Nlv+1
        If (iabs(Kl4) == 2) n01=1
        If (n01 > Nd0) Then
            continue
        Else
            Do i=n01,Nd0
                s=0.d0
                Do k=1,Nd0
                    s=s+Z1(k,i)*C(k)
                End Do
                t=val-E1(i)
                If (dabs(t) < 1.d-6) Cycle
                B1(1:Nd0)=B1(1:Nd0)+(s/t)*Z1(1:Nd0,i)
            End Do
        End If
        ! Normalization of vector B:
        s=0.d0
        Do i=1,Nd
            s=s+B1(i)**2
        End Do
        s=1.d0/dsqrt(s)
        ArrB(1:Nd,J1)=B1(1:Nd)*s
        Return
    End Subroutine Dvdsn

    Subroutine Mxmpy(ip,mype,npes)
        ! This subroutine evaluates the products Q_i = H_ij * B_j
        ! if ip=1, products are stored in the second block of ArrB
        ! if ip=2, products are stored in the third block of ArrB
        Use mpi
        Implicit None
        Integer :: i, i1, i2, k, ierr, ip, nlp, n, mype, npes, mpierr, &
                    i2min, i2max, i1min, i1max, kd
        Real :: start_time, end_time
        Real(dp) :: t
        Integer*8 :: l8
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
        nlp=ip*Nlv

        i2min=1+nlp
        i2max=Nlv+nlp
        i1min=i2min-Nlv
        i1max=i2max-Nlv
        kd=0
        Call MPI_Bcast(Diag(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        If (Diag(1)==0.d0) Then
            Diag(2:Nd)=0.d0
            kd=1
        End If

        Do i=1,nlp ! 
            Call MPI_Bcast(ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        End Do
        ArrB(1:Nd,i2min:i2max)=0.d0
        Do l8=1, ih8H
            n=Hamil%n(l8)
            k=Hamil%k(l8)
            t=Hamil%t(l8)
            If (n == k) Then
                t=t-Hmin
                If (kd == 1) Diag(n)=t
            End If
            ArrB(n,i2min:i2max)=ArrB(n,i2min:i2max)+t*ArrB(k,i1min:i1max)
            If (n /= k) ArrB(k,i2min:i2max)=ArrB(k,i2min:i2max)+t*ArrB(n,i1min:i1max)
        End Do
        ! Only master core needs ArrB for other calculations
        Do i=i2min,i2max
            If (mype==0) Then
                Call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            Else
                Call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            End If
        End Do
        If (kd==1) Then
            If (mype==0) Then
                Call MPI_Reduce(MPI_IN_PLACE, Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, & 
                                MPI_COMM_WORLD, mpierr)
            Else
                Call MPI_Reduce(Diag(1:Nd), Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, & 
                                MPI_COMM_WORLD, mpierr)
            End If   
        End If
        Return
    End Subroutine Mxmpy

    Subroutine Ortn (j)
        ! Orthogonalization of J-th vector to J-1 previous ones.
        ! Normalization of the J-th vector.
        Implicit None
        Integer :: i, l, j, jm, it
        Real(dp) :: s, smax2, smax1, critN, ortho
        ! Orthogonality criteria:
        ortho=1.d-6
        ! Critical value for norma of new vector:
        critN=1.d-9
        jm=j-1
        smax2=1.d10

        it=0
        do
            it=it+1
            smax1=0.d0
            Do l=1,jm
                If (l > Nlv .and. Iconverge(l-Nlv)==1) Cycle
                s=0.d0
                Do i=1,Nd
                    s=s+ArrB(i,j)*ArrB(i,l)
                End Do
                If (dabs(s) > smax1) smax1=dabs(s)
                ArrB(1:Nd,j)=ArrB(1:Nd,j)-s*ArrB(1:Nd,l)
            End Do
            If (smax1 < ortho) Exit
            If (smax1 > smax2 .or. it > 3) Then
                Write (6,'(4X,"Fail of orthogonalization of vector ",I2)') j
                Ifail=j
                Return
            End If
            smax2=smax1
        End Do

        ! Normalization of J-th vector:
        s=0.d0
        Do i=1,Nd
            s=s+ArrB(i,j)*ArrB(i,j)
        End Do
        If (s < critN) Then
            Write (6,'(4X,"Fail of normalization of vector ",I2)') j
            Ifail=j
        End If
        s=1.d0/dsqrt(s)
        ArrB(1:Nd,j)=ArrB(1:Nd,j)*s
        Return
    End Subroutine Ortn

    Subroutine Prj_J (lin, num, lout, trsd, mype, npes) 
        Use mpi
        Use formj2, Only : F_J2
        Use determinants, Only : Gdet
        Implicit None
        Integer :: ierr, i, n, k, j, i2, i1, jx, it, iter, nsj, lout, num, lin, jav, j_av, mpierr, imin, imax
        Integer, optional :: mype, npes
        Integer*8 :: nj
        Real(dp) :: err1, err, sj, t, f, s, s1, a, aj, trsd, start, end
        logical :: proceed
        Integer, allocatable, dimension(:) :: idet1, idet2
        !     - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(idet1(Ne),idet2(Ne))
        ierr=1
        If (mype == 0 .and. lout+num-1 > IPlv) Then
            Write( *,'(" Prj_J warning:", /," to project",I3," levels I need",I3, &
                " vectors & have Only",I3)') Nlv,lout+num-1,IPlv
            Write(11,'(" Prj_J warning:", /," to project",I3," levels I need",I3, &
                " vectors & have Only",I3)') Nlv,lout+num-1,IPlv
            Stop
        End If
        jav=2*XJ_av+0.1d0
        If (mype == 0) Then
            Write( *,'(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)') num,XJ_av
            Write(11,'(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)') num,XJ_av
        End If

        Call MPI_Bcast(Njd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jt, IPjd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iconverge, IPlv, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)        
        imin=lin
        imax=lout-1
        Do j=imin,imax
            Call MPI_Bcast(ArrB(1:Nd,j), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        End Do
        Do iter=1,Njd+1          !# each iteration eliminates one J
            it=Njd+1-iter        !## starting from J_max down to J_min
            If (it /= 0) Then
                jx=Jt(it)            !# - angular momentum (doubled)
                f=1.d0/(jx*(jx+2))   !# normalization factor for subtraction
            End If
            If (jx /= jav .or. it == 0) Then
                Do i=1,num
                    i2=lout+i-1
                    ArrB(1:Nd,i2)=0.d0
                End Do
                Do j=1,ij8J
                    n=Jsq%n(j)
                    k=Jsq%k(j)
                    t=Jsq%t(j)
                    Do i=1,num
                        proceed=(lin-1)*Iconverge(i)==0
                        If (proceed) Then
                            i1=lin+i-1
                            i2=lout+i-1
                            ArrB(n,i2)=ArrB(n,i2)+t*ArrB(k,i1)
                            If (n /= k) ArrB(k,i2)=ArrB(k,i2)+t*ArrB(n,i1)
                        End If
                    End Do
                End Do

                Do i=lout,lout+num-1
                    If (mype==0) Then
                        Call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                          MPI_COMM_WORLD, mpierr)
                    Else
                        Call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                          MPI_COMM_WORLD, mpierr)
                    End If
                    Call MPI_Bcast(ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                End Do
                err1=0.d0
            
                Do i=1,num
                    proceed=(lin-1)*Iconverge(i)==0
                    If (proceed) Then
                        i1=lin+i-1
                        i2=lout+i-1
                        sj=0.0             !# <X1|J**2|X1> is Used to check convergence
                        Do n=1,Nd                       
                            sj=sj+ArrB(n,i1)*ArrB(n,i2) 
                        End Do
                        aj=dsqrt(1+sj) - 1               !# eigenvalue for 2*J
                        err=dabs(aj-jav)                 !# difference from target
                        err1=dmax1(err1,err)
                    End If
                End Do
                If (err1 < trsd) Then             !# convergency check
                    ierr=0
                    Exit
                End If
                If (it /= 0) Then !# eliminating jx component
                    Do i=1,num
                        proceed=(lin-1)*Iconverge(i)==0
                        If (proceed) Then
                            i1=lin+i-1
                            i2=lout+i-1
                            s=0.d0
                            Do n=1,Nd                        !# X1 = X1 - f*Xj
                                a=ArrB(n,i1)-f*ArrB(n,i2)      !## f=1/(2j*(2j+2))
                                s=s+a*a
                                ArrB(n,i1)=a
                            End Do
                            s1=1.d0/dsqrt(s)
                            ArrB(1:Nd,i1)=ArrB(1:Nd,i1)*s1  !# normalization of X1
                        End If
                    End Do
                End If

            End If
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End Do
        If (mype == 0) Then
            Write( *,'(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)') Njd+1-it,err1,ierr
            Write(11,'(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)') Njd+1-it,err1,ierr
            If (ierr /= 0) Then
                Write(*,*) ' Wrong J values for Probe vectors '
            End If
        End If
        Return
    End Subroutine Prj_J

    Subroutine Vread (b,j)
        ! This subroutine reads eigenvector j from ArrB(:,j)
        Implicit None
        Integer :: i, j
        Real(dp), dimension(Nd) :: b
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        b(1:Nd)=ArrB(1:Nd,j)
        Return
    End Subroutine Vread

    Subroutine VWrite (b,j)
        ! This subroutine writes eigenvector j to ArrB(:,j)
        Implicit None
        Integer :: i, j
        Real(dp), dimension(Nd) :: b
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ArrB(1:Nd,j)=b(1:Nd)
        Return
    End Subroutine VWrite
    
End Module davidson
