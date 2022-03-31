Module davidson
    ! This module contains main subroutines required for the Davidson procedure.

    Use conf_variables

    Implicit None

    Private

    Public :: Init4, FormB0, FormB, FormBskip, WriteXIJ, FormP, AvgDiag, Dvdsn, Mxmpy, Ortn, Prj_J

  Contains
    
    Subroutine Init4(mype)
        ! This subroutine constructs the initial approximation 
        ! by selecting configurations specified by parameter Nc4.
        ! The initial approximation Hamilonian is stored in the 
        ! matrix Z1 and is constructed by selecting the top-left 
        ! block of the full Hamiltonian matrix H.
        ! The diagonal elements of H are stored in the array Diag.
        Use mpi_f08
        Use determinants, Only : calcNd0
        Implicit None
        Integer  :: k, n, mype, mpierr
        Integer(Kind=int64) :: l8
        Real(type_real) :: t
        Character(Len=256) :: strfmt

        If (mype == 0) Then
            strfmt = '(3X,"Starting approx. includes ",I5," conf.,",I6," det.")'
            Write( 6,strfmt) Nc1,Nd0
            Write(11,strfmt) Nc1,Nd0
        End If

        ! Construct the Hamilonian in the initial approximation
        Diag=0_type_real
        Z1=0_type_real
        Do l8=1,ih8
            n=Hamil%ind1(l8)
            k=Hamil%ind2(l8)
            t=Hamil%val(l8)
            If (n == k) t=t-Hamil%minval
            If (n <= Nd0) Then
                Z1(n,k)=t
                Z1(k,n)=t
            Else
                ! If Kl=1, matrix elements will be spread among num_cores instead of num_cores-1,
                ! so we must ensure that all determinants below Nd0 are accounted for
                If (Kl == 1) Then
                    Cycle
                Else If (Kl == 3) Then
                    Cycle
                Else
                    Exit
                End If
            End If
        End Do

        Call MPI_AllReduce(MPI_IN_PLACE, Z1, Nd0*Nd0, mpi_type_real, MPI_SUM, &
                                MPI_COMM_WORLD, mpierr)

    End Subroutine Init4

    Subroutine FormB0(mype)
        ! This subroutine constructs the initial eigenvectors for the Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are written to CONF.XIJ 
        Use mpi_f08
        Use formj2, Only : J_av
        Implicit None
        Integer :: i, j, idum, ndpt, ierr, num, nskip, mpierr, mype
        Real(dp) :: dummy
        Real(type_real) :: xj
        Character(Len=128) :: strfmt

        num=0
        nskip=0

        If (mype==0) Open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')

        Call MPI_Bcast(Z1, Nd0*Nd0, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)

        B1=0_type_real
        If (abs(Kl4) /= 2) Then ! If not reading CONF.XIJ
            Do j=1,Nd0
                B1(1:Nd0)=Z1(1:Nd0,j)
                Call J_av(B1,Nd0,xj,ierr)
                If (ierr == 0) Then
                    num=num+1
                    E(num)=-(E1(j)+Hamil%minval)
                    Tj(num)=xj
                    ArrB(1:Nd,num)=B1(1:Nd)
                    If (num >= Nlv) Exit
                Else
                    nskip=nskip+1
                    If (mype==0) Write(*,*) '  skip Ej =',-(E1(j)+Hamil%minval)
                End If
            End Do
        Else ! Reading CONF.XIJ
            If (mype == 0) Then
                Read(17,err=210) dummy,dummy,ndpt
 210            Write(*,*) ' Vector length = ',ndpt
                Rewind(17)
                Do j=1,Nlv
                    Read(17,end=220,err=220) E(j),Tj(j),idum,(B1(i),i=1,ndpt)
                    num=j
                    ArrB(1:Nd,j)=B1(1:Nd)
                End Do
            End If
            Call MPI_Bcast(num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        End If
 220    if (mype==0) Rewind(17)
        Nlv=num

        ! Write eigenvalues and eigenvectors of initial approximation to CONF.XIJ
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
        Integer :: k, l
        Real(type_real) :: t

        ! Eigenvectors obtained by Davidson procedure are not written to CONF.XIJ file
        Do k=1,Nlv
            B2=0_type_real
            Do l=1,2*Nlv
                t=P(l,k)
                B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
            End Do
            t=-(E(k)+Hamil%minval)
            Tk(k)=t
            ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        End Do

        Do k=1,Nlv
            ArrB(1:Nd,k)=ArrB(1:Nd,2*Nlv+k)
        End Do
        
        Write ( 6,*) ' FormBskip: Vectors not saved'
        Return
    End Subroutine FormBskip

    Subroutine FormB
        ! This subroutine forms eigenvectors for next Davidson iteration
        ! and stores them in the first block of ArrB.
        ! Eigenvectors are written to CONF.XIJ.
        Implicit None
        Integer :: i, k, l
        Real(type_real) :: t

        ! Eigenvectors obtained by Davidson procedure are written to CONF.XIJ file
        Open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        Do k=1,Nlv
            B2=0_type_real
            Do l=1,2*Nlv
                t=P(l,k)
                B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
            End Do
            t=-(E(k)+Hamil%minval)
            Tk(k)=t
            Write(17) t,Tj(k),Nd,(B2(i),i=1,Nd)
            ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        End Do
        Do k=1,Nlv
            ArrB(1:Nd,k)=ArrB(1:Nd,2*Nlv+k)
        End Do
        close (unit=17)
        Write ( 6,*) ' FormB: Vectors saved'
        Return
    End Subroutine FormB

    Subroutine WriteXIJ
        Implicit None
        Integer :: i, k
        ! Eigenvectors obtained by Davidson procedure are written to CONF.XIJ file
        Open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        Do k=1,Nlv
            Write(17) Tk(k),Tj(k),Nd,(ArrB(i,k),i=1,Nd)
        End Do
        Close(17)
        Write ( 6,*) ' WriteXIJ: Vectors saved'
    End Subroutine WriteXIJ

    Subroutine FormP(kp, vmax)
        ! Calculation of upper-left block (Nlv X Nlv) of matrix P for kp=1 
        ! Calculation of other three blocks (2Nlv X 2Nlv) of matrix P for kp=2
        Implicit None
        Integer, Intent(In) :: kp
        Real(type_real), Intent(In) :: vmax

        Integer :: i, kx, k1, j, k, i1
        Real(type_real) :: x

        If (kp == 1) Then
            Do i=1,Nlv
              Do k=1,i
                k1=k+Nlv
                x=0_type_real
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
                    x=0_type_real
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
                    x=0_type_real
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
                    P(i1,k)=0_type_real
                    P(k,i1)=0_type_real
                End Do
                P(i1,i1)=2*vmax-P(i,i)+1.d0
              End If
            End Do
        End If
        Return
    End Subroutine FormP
    
    Subroutine AvgDiag
        ! This subroutine averages the diagonal over configurations
        Implicit None

        Integer :: ic, id, id0, id1, id2
        Real(type_real) :: ss

        id0=1
        Do ic=1,Nc
            id1=Ndc(ic)
            id2=id0+id1-1
            If (id1 > 0) Then
                ss=0_type_real
                Do id=id0,id2
                    ss=ss+Diag(id)
                End Do
                ss=ss/id1
                Diag(id0:id2)=ss
                id0=id0+id1
            End If
        End Do

        If (id2 == Nd) Then
            Write(*,*) ' Diagonal averaged over rel. configurations'
        Else
            Write(*,*) ' Error: id2=',id2,' Nc=',Nc
            Stop
        End If

    End Subroutine AvgDiag

    Subroutine Dvdsn (j, cnx)
        ! Main part of the Davidson algorithm. 
        ! Formation of the J-th residual vector and the corresponding new probe vector.
        ! New probe vectors are stored in the second block of ArrB.
        Implicit None
        Integer, Intent(In) :: j
        Real(dp), Intent(InOut) :: cnx

        Integer :: j1, i, k, n01
        Real(type_real) :: val, t, cnorm, s
        character(len=9) :: char

        cnorm=0_type_real
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
            B1(i)=0_type_real
            Cycle
        End Do
        cnorm=sqrt(abs(cnorm))
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
        If (cnorm > cnx) cnx=Real(cnorm,kind=dp)
        If (Iconverge(j)==1) Then
            B1(1:Nd)=0_type_real
            Return
        End If
        n01=Nlv+1
        If (iabs(Kl4) == 2) n01=1
        If (n01 > Nd0) Then
            continue
        Else
            Do i=n01,Nd0
                s=0_type_real
                Do k=1,Nd0
                    s=s+Z1(k,i)*C(k)
                End Do
                t=val-E1(i)
                If (abs(t) < 1.d-6) Cycle
                B1(1:Nd0)=B1(1:Nd0)+(s/t)*Z1(1:Nd0,i)
            End Do
        End If
        ! Normalization of vector B:
        s=0_type_real
        Do i=1,Nd
            s=s+B1(i)**2
        End Do
        s=1.d0/sqrt(s)
        ArrB(1:Nd,J1)=B1(1:Nd)*s
        Return
    End Subroutine Dvdsn

    Subroutine Mxmpy(ip, mype)
        ! This subroutine evaluates the products Q_i = H_ij * B_j
        ! if ip=1, products are stored in the second block of ArrB
        ! if ip=2, products are stored in the third block of ArrB
        Use mpi_f08
        Use mpi_utils, Only : BroadcastD
        Implicit None
        Integer, Intent(In) :: ip, mype

        Integer :: i, k, nlp, n, mpierr, i2min, i2max, i1min, i1max, kd
        Real(type_real) :: t
        Integer(Kind=int64) :: l8, count

        nlp=ip*Nlv

        i2min=1+nlp
        i2max=Nlv+nlp
        i1min=i2min-Nlv
        i1max=i2max-Nlv
        kd=0
        Call MPI_Bcast(Diag(1), 1, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
        If (Diag(1)==0_type_real) Then
            Diag(2:Nd)=0_type_real
            kd=1
        End If

        count=Nd*nlp
        !Call BroadcastD(ArrB, count, 0, 0, MPI_COMM_WORLD, mpierr)
        If (mype /= 0) ArrB = Real(0, kind=type_real)
        Do i=1,nlp
            Call MPI_AllReduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, mpi_type_real, MPI_SUM, &
                                  MPI_COMM_WORLD, mpierr)
        End Do

        ArrB(1:Nd,i2min:i2max)=0_type_real
        Do l8=1, ih8
            n=Hamil%ind1(l8)
            k=Hamil%ind2(l8)
            t=Hamil%val(l8)
            If (n == k) Then
                t=t-Hamil%minval
                If (kd == 1) Diag(n)=t
            End If
            ArrB(n,i2min:i2max)=ArrB(n,i2min:i2max)+t*ArrB(k,i1min:i1max)
            If (n /= k) ArrB(k,i2min:i2max)=ArrB(k,i2min:i2max)+t*ArrB(n,i1min:i1max)
        End Do

        ! Only master core needs ArrB for other calculations
        Do i=i2min,i2max
            If (mype==0) Then
                Call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, mpi_type_real, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            Else
                Call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, mpi_type_real, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            End If
        End Do
        If (kd==1) Then
            If (mype==0) Then
                Call MPI_Reduce(MPI_IN_PLACE, Diag(1:Nd), Nd, mpi_type_real, MPI_SUM, 0, & 
                                MPI_COMM_WORLD, mpierr)
            Else
                Call MPI_Reduce(Diag(1:Nd), Diag(1:Nd), Nd, mpi_type_real, MPI_SUM, 0, & 
                                MPI_COMM_WORLD, mpierr)
            End If 
        End If
        Return
    End Subroutine Mxmpy

    Subroutine Ortn(j, ifail)
        ! Orthogonalization of J-th vector to J-1 previous ones.
        ! Normalization of the J-th vector.
        Implicit None
        Integer, Intent(In)  :: j
        Integer, Intent(Out) :: ifail

        Integer :: i, l, jm, it
        Real(dp) :: smax2, smax1, critN, ortho
        Real(type_real) :: s
        
        ! Orthogonality criteria:
        ortho=1.d-6

        ! Critical value for normalization of new vector:
        critN=1.d-9
        jm=j-1
        smax2=1.d10

        ! Orthogonalization of J-th vector
        it=0
        Do
            it=it+1
            smax1=0_type_real
            Do l=1,jm
                If (l > Nlv) Then
                    If (Iconverge(l-Nlv)==1) Cycle
                End If
                s=0_type_real
                Do i=1,Nd
                    s=s+ArrB(i,j)*ArrB(i,l)
                End Do
                If (abs(s) > smax1) smax1=abs(s)
                ArrB(1:Nd,j)=ArrB(1:Nd,j)-s*ArrB(1:Nd,l)
            End Do
            If (smax1 < ortho) Exit
            If (smax1 > smax2 .or. it > 3) Then
                Write (6,'(4X,"Fail of orthogonalization of vector ",I2)') j
                ifail=j
                Return
            End If
            smax2=smax1
        End Do

        ! Normalization of J-th vector:
        s=0_type_real
        Do i=1,Nd
            s=s+ArrB(i,j)*ArrB(i,j)
        End Do
        If (s < critN) Then
            Write (6,'(4X,"Fail of normalization of vector ",I2)') j
            ifail=j
        End If
        s=1.d0/sqrt(s)
        ArrB(1:Nd,j)=ArrB(1:Nd,j)*s
        Return
    End Subroutine Ortn

    Subroutine Prj_J(lin, num, lout, trsd, mype) 
        Use mpi_f08
        Use determinants, Only : Gdet
        Implicit None
        Integer :: ierr, i, n, k, j, i2, i1, jx, it, iter, lout, num, lin, jav, mpierr, imin, imax
        Integer, optional :: mype
        Real(dp) :: err1, err, f, trsd
        Real(type_real) :: t, sj, aj, s, s1, a
        logical :: proceed
        Character(Len=256) :: strfmt

        ierr=1
        If (mype == 0 .and. lout+num-1 > IPlv) Then
            strfmt = '(" Prj_J warning:", /," to project",I3," levels I need",I3," vectors & have only",I3)'
            Write( *,strfmt) Nlv,lout+num-1,IPlv
            Write(11,strfmt) Nlv,lout+num-1,IPlv
            Stop
        End If
        jav=2*XJ_av+0.1d0
        If (mype == 0) Then
            strfmt = '(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)'
            Write( *,strfmt) num,XJ_av
            Write(11,strfmt) num,XJ_av
        End If

        Call MPI_Bcast(Njd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Jt, IPjd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        Call MPI_Bcast(Iconverge, Nlv, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        imin=lin
        imax=lout-1
        Do j=imin,imax
            Call MPI_Bcast(ArrB(1:Nd,j), Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
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
                    ArrB(1:Nd,i2)=0_type_real
                End Do
                Do j=1,ij8
                    n=Jsq%ind1(j)
                    k=Jsq%ind2(j)
                    t=Jsq%val(j)
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
                        Call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, mpi_type_real, MPI_SUM, 0, &
                                          MPI_COMM_WORLD, mpierr)
                    Else
                        Call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, mpi_type_real, MPI_SUM, 0, &
                                          MPI_COMM_WORLD, mpierr)
                    End If
                    Call MPI_Bcast(ArrB(1:Nd,i), Nd, mpi_type_real, 0, MPI_COMM_WORLD, mpierr)
                End Do
                err1=0_type_real
            
                Do i=1,num
                    proceed=(lin-1)*Iconverge(i)==0
                    If (proceed) Then
                        i1=lin+i-1
                        i2=lout+i-1
                        sj=0.0             !# <X1|J**2|X1> is Used to check convergence
                        Do n=1,Nd                       
                            sj=sj+ArrB(n,i1)*ArrB(n,i2) 
                        End Do
                        aj=sqrt(1+sj) - 1               !# eigenvalue for 2*J
                        err=abs(aj-jav)                 !# difference from target
                        err1=max1(err1,err)
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
                            s=0_type_real
                            Do n=1,Nd                        !# X1 = X1 - f*Xj
                                a=ArrB(n,i1)-f*ArrB(n,i2)      !## f=1/(2j*(2j+2))
                                s=s+a*a
                                ArrB(n,i1)=a
                            End Do
                            s1=1.d0/sqrt(s)
                            ArrB(1:Nd,i1)=ArrB(1:Nd,i1)*s1  !# normalization of X1
                        End If
                    End Do
                End If

            End If
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        End Do
        If (mype == 0) Then
            strfmt = '(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)'
            Write( *,strfmt) Njd+1-it,err1,ierr
            Write(11,strfmt) Njd+1-it,err1,ierr
            If (ierr /= 0) Then
                Write(*,*) ' Wrong J values for Probe vectors '
            End If
        End If
        Return
    End Subroutine Prj_J
    
End Module davidson
