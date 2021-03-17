module davidson

    use conf_variables

    implicit none

    contains
    subroutine Init4
        implicit none
        integer  :: k, ierr, n, n1, n2, ic, ic1
        real(dp) :: t
        integer*8 :: num8, Ibuf0, ibuf, l8
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Dimension of the approximation to start with:
        n1=0
        ic1=0
        n2=Nd0-1
        do ic=1,Nc4
            n1=n1+Ndc(ic)
            if (n1 < Nd0) then
                n2=n1
                ic1=ic
            end if
        end do
        Nd0=n2
        write( 6,'(3X,"Starting approx. includes ",I3," conf.,", &
             I4," det.")') ic1,Nd0
        write(11,'(3X,"Starting approx. includes ",I3," conf.,", &
             I4," det.")') ic1,Nd0
        ! Construct the Hamiltonian in the initial approximation:
        Diag=0.d0
        Z1=0.d0
        do l8=1,NumH
            n=H_n0(l8)
            k=H_k0(l8)
            t=H_t0(l8)
            if (n == k) t=t-Hmin
            if (n <= Nd0) then
                Z1(n,k)=t
                Z1(k,n)=t
            else
                exit
            end if
        end do
        return
    end subroutine Init4

!     =================================================
    subroutine Hould(n,ndim,Ee,Dd,Zz)
        implicit none
        integer :: n, ii, k, j, l, i, jj, m1, ndim, im, im1, li, k1, l1
        real(dp) :: tol, f, g, b, r, em, h, hh, c, ei, s, p, di, eps
        real(dp), dimension(:), allocatable :: Ee, Dd
        real(dp), dimension(:,:), allocatable :: Zz
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

    subroutine FormB0(mype,npes)
       use mpi
       use formj2, only : J_av
       implicit none
       integer :: i, j, idum, ndpt, ierr, num, nskip, mpierr, mype, npes
       real(dp) :: dummy, xj
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        num=0
        nskip=0
        
        if (mype==0) open(unit=17,file='CONF.XIJ',status='UNKNOWN',form='UNFORMATTED')

        if (Kv == 4) call MPI_Bcast(Z1,Nd0*Nd0,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)

        if (iabs(Kl4) /= 2) then
           do j=1,Nd0
              do i=1,Nd
                 if (j == 1)   B1(i)=0.d0
                 if (i <= Nd0) B1(i)=Z1(i,j)
              end do
              select case(Kv)
                case(4)
                  call J_av(B1,Nd0,xj,ierr,mype,npes)
                case(3)
                  call J_av(B1,Nd0,xj,ierr)
              end select
              if(ierr == 0) then
                num=num+1
                E(num)=-(E1(j)+Hmin)
                Tj(num)=xj
                ArrB(1:Nd,num)=B1(1:Nd)
                if(num >= Nlv) goto 220
              else
                nskip=nskip+1
                if (mype==0) write(*,*) '  skip Ej =',-(E1(j)+Hmin)
              end if
           end do
        else
           B1=0.d0
           read(17,err=210) dummy,dummy,ndpt
 210       write(*,*) ' Vector length = ',ndpt
           rewind(17)
           do J=1,Nlv
              read(17,end=220,err=220) E(J),Tj(J),idum,(B1(i),i=1,ndpt)
              num=j
              ArrB(1:Nd,j)=B1(1:Nd)
           end do
        end if
 220    rewind(17)
        Nlv=num
        if (mype == 0) then
          do j=1,Nlv
             write( 6,'(1X,"E(",I2,") =",F12.6," Jtot =",F10.6)') j,E(j),Tj(j)
             write(11,'(1X,"E(",I2,") =",F12.6," Jtot =",F10.6)') j,E(j),Tj(j)
             if (kXIJ > 0) write(17) E(j),Tj(j),Nd,(ArrB(i,j),i=1,Nd)
          end do
          close (unit=17)
          if (nskip > 0) then
             write(*,*) ' Selection is done', nskip,' vectors skiped'
             write(11,*) ' Selection is done', nskip,' vectors skiped'
          end if
        end if
       return
    end subroutine FormB0

    subroutine FormBskip
       implicit none
       integer :: i, idum, k, l
       real(dp) :: t, tjk, ek
       !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! Eigenvectors obtained by Davidson procedure are not written to CONF.XIJ file
        do k=1,Nlv
          B2=0.d0
          do l=1,2*Nlv
              t=P(l,k)
              B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
          end do
          t=-(E(k)+Hmin)
          Tk(k)=t
          ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        end do
        ArrB(1:Nd,1:Nlv)=ArrB(1:Nd,2*Nlv+1:3*Nlv)
        write ( 6,'(1X,"FormBskip: Vectors not saved")')
       return
    end subroutine FormBskip

    subroutine FormB
       implicit none
       integer :: i, idum, k, l
       real(dp) :: t, tjk, ek
       !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! Eigenvectors obtained by Davidson procedure are written to CONF.XIJ file
        open(unit=17,file='CONF.XIJ',status='OLD',form='UNFORMATTED')
        do k=1,Nlv
          B2=0.d0
          do l=1,2*Nlv
              t=P(l,k)
              B2(1:Nd)=B2(1:Nd)+t*ArrB(1:Nd,l)
          end do
          t=-(E(k)+Hmin)
          Tk(k)=t
          write(17) t,Tj(k),Nd,(B2(i),i=1,Nd)
          ArrB(1:Nd,k+2*Nlv)=B2(1:Nd)
        end do
        ArrB(1:Nd,1:Nlv)=ArrB(1:Nd,2*Nlv+1:3*Nlv)
        close (unit=17)
        write ( 6,'(1X,"FormB: Vectors saved")')
       return
    end subroutine FormB

    subroutine FormP (kp)
      ! Calculation of upper-left block (Nlv X Nlv) of matrix P for kp=1 
      ! Calculation of other three blocks (2Nlv X 2Nlv) of matrix P for kp=2
      implicit none
      integer :: ix, i, kx, kp, k1, k2, j, k, i1
      real(dp) :: x
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (kp == 1) then
        do i=1,Nlv
          do k=1,i
            k1=k+Nlv
            x=0.d0
            do j=1,Nd
              x=x+ArrB(j,i)*ArrB(j,k1)
            end do
            P(i,k)=x
            P(k,i)=x
          end do
        end do
      else                  ! kp=2
        do i=1,Nlv
          kx=2*Nlv
          do k=Nlv+1,kx
            k1=k+Nlv
            x=0.d0
            do j=1,Nd
              x=x+ArrB(j,i)*ArrB(j,k1)
            end do
            P(i,k)=x
            P(k,i)=x
          end do
        end do
        do i=Nlv+1,2*Nlv
          do k=Nlv+1,i
            k1=k+Nlv
            x=0.d0
            do j=1,Nd
              x=x+ArrB(j,i)*ArrB(j,k1)
            end do
            P(i,k)=x
            P(k,i)=x
          end do
        end do
        do i=1,Nlv                   ! Cleanup for converged levels
          if (Iconverge(i) == 1) then
            i1=i+Nlv
            do k=1,2*Nlv
              P(i1,k)=0.d0
              P(k,i1)=0.d0
            end do
            P(i1,i1)=2*vmax-P(i,i)+1.d0
          end if
        end do
      end if
      return
    end subroutine FormP
    
    subroutine Dvdsn (j)
    ! Main part of the Davidson algorithm. 
    ! Formation of the J-th residual vector and the corresponding new probe vector
      implicit none
      integer :: j1, i, j, l, k, n01
      real(dp) :: val, t, cnorm, s
      character(len=9) :: char
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      cnorm=0.d0
      j1=j+Nlv
      val=P(j,j)
      do i=1,Nd
        t=ArrB(i,j1)-val*ArrB(i,j)
        if (i <= Nd0) C(i)=t
        cnorm=cnorm+t*t
        if (i > Nd0) then
          B1(i)=t/(val-Diag(i))
          cycle
        end if
        B1(i)=0.d0
        cycle
      end do
      cnorm=dsqrt(dabs(cnorm))
      Iconverge(j)=0
      char='         '
      if(cnorm < Crt4) then
        if (kdavidson == 1) then
          char='stopped  '
          Iconverge(j)=1
        else
          char='converged'
        end if
      end if
      write ( 6,'(4X,"|| C(",I2,") || = ",F10.7,2X,A9)') j,cnorm,char
      write (11,'(4X,"|| C(",I2,") || = ",F10.7,2X,A9)') j,cnorm,char
      if (cnorm > cnx) cnx=cnorm
      if (Iconverge(j)==1) then
        B1(1:Nd)=0.d0
        return
      end if
      n01=Nlv+1
      if (iabs(Kl4) == 2) n01=1
      if (n01 > Nd0) then
        continue
      else
        do i=n01,Nd0
          s=0.d0
          do k=1,Nd0
             s=s+Z1(k,i)*C(k)
          end do
          t=val-E1(i)
          if (dabs(t) < 1.d-6) cycle
          B1(1:Nd0)=B1(1:Nd0)+(s/t)*Z1(1:Nd0,i)
        end do
      end if
      ! Normalization of vector B:
      s=0.d0
      do i=1,Nd
         s=s+B1(i)**2
      end do
      s=1.d0/dsqrt(s)
      ArrB(1:Nd,J1)=B1(1:Nd)*s
      return
    end subroutine Dvdsn

    subroutine Mxmpy(ip,mype,npes)
        ! This subroutine evaluates the products Q_i = H_ij * B_j
        use mpi
        implicit none
        integer :: i, i1, i2, k, ierr, ip, nlp, n, mype, npes, mpierr, &
                    i2min, i2max, i1min, i1max, kd
        real :: start_time, end_time
        real(dp) :: t
        integer*8 :: l8
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
        nlp=ip*Nlv

        i2min=1+nlp
        i2max=Nlv+nlp
        i1min=i2min-Nlv
        i1max=i2max-Nlv
        kd=0
        call MPI_Bcast(Diag(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        if (Diag(1)==0.d0) then
          Diag(2:Nd)=0.d0
          kd=1
        end if

        do i=1,nlp ! 
            call MPI_Bcast(ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        end do
        ArrB(1:Nd,i2min:i2max)=0.d0
        do l8=1, ih8H
            n=H_n(l8)
            k=H_k(l8)
            t=H_t(l8)
            if (n == k) then
                t=t-Hmin
                if (kd == 1) Diag(n)=t
            end if
            ArrB(n,i2min:i2max)=ArrB(n,i2min:i2max)+t*ArrB(k,i1min:i1max)
            if (n /= k) ArrB(k,i2min:i2max)=ArrB(k,i2min:i2max)+t*ArrB(n,i1min:i1max)
        end do
        ! only master core needs ArrB for other calculations
        do i=i2min,i2max
            if (mype==0) then
              call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            else
              call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                  MPI_COMM_WORLD, mpierr)
            end if
        end do
        if (kd==1) then
          if (mype==0) then
            call MPI_Reduce(MPI_IN_PLACE, Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, & 
                              MPI_COMM_WORLD, mpierr)
          else
            call MPI_Reduce(Diag(1:Nd), Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, & 
                              MPI_COMM_WORLD, mpierr)
          end if   
        end if
        return
    end subroutine Mxmpy

      Subroutine Ortn (j)
      ! Orthogonalization of J-th vector to J-1 previous ones.
      ! Normalization of the J-th vector.
       implicit none
       integer :: i, l, j, jm, it
       real(dp) :: s, smax2, smax1, critN, ortho
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
          do l=1,jm
             if (l > Nlv .and. Iconverge(l-Nlv)==1) cycle
             s=0.d0
             do i=1,Nd
                s=s+ArrB(i,j)*ArrB(i,l)
             end do
             if (dabs(s) > smax1) smax1=dabs(s)
             ArrB(1:Nd,j)=ArrB(1:Nd,j)-s*ArrB(1:Nd,l)
          end do
          if (smax1 < ortho) exit
          if (smax1 > smax2 .or. it > 3) then
            write (6,'(4X,"Fail of orthogonalization of vector ",I2)') j
            Ifail=j
            return
          end if
          smax2=smax1
        end do

        ! Normalization of J-th vector:
        s=0.d0
        do i=1,Nd
           s=s+ArrB(i,j)*ArrB(i,j)
        end do
        if (s < critN) then
          write (6,'(4X,"Fail of normalization of vector ",I2)') j
          Ifail=j
        end if
        s=1.d0/dsqrt(s)
        ArrB(1:Nd,j)=ArrB(1:Nd,j)*s
       return
    end subroutine Ortn

    subroutine Prj_J (lin, num, lout, ierr, trsd, mype, npes) 
       use mpi
       use formj2, only : F_J2
       use determinants, only : Gdet
       implicit none
       integer :: ierr, i, n, k, j, i2, i1, jx, it, iter, nsj, lout, num, lin, jav, j_av, mpierr
       integer, optional :: mype, npes
       integer*8 :: nj
       real(dp) :: err1, err, sj, t, f, s, s1, a, aj, trsd, start, end
       logical :: proceed
       integer, allocatable, dimension(:) :: idet1, idet2
!     - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(idet1(Ne),idet2(Ne))
        ierr=1
        if (lout+num-1 > IPlv) then
          write( *,'(" Prj_J warning:", /," to project",I3," levels I need",I3, &
            " vectors & have only",I3)') Nlv,lout+num-1,IPlv
          write(11,'(" Prj_J warning:", /," to project",I3," levels I need",I3, &
            " vectors & have only",I3)') Nlv,lout+num-1,IPlv
          stop
        end if
        jav=2*XJ_av+0.1d0
        write( *,'(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)') num,XJ_av
        write(11,'(4X,"Prj_J: projecting",I3," vectors on subspace J=",F5.1)') num,XJ_av
        if (Kv == 3) open(unit=18,file='CONF.JJJ',status='OLD',form='UNFORMATTED')
        do iter=1,Njd+1          !# each iteration eliminates one J
          it=Njd+1-iter        !## starting from J_max down to J_min
          if (it /= 0) then
            jx=Jt(it)            !# - angular momentum (doubled)
            f=1.d0/(jx*(jx+2))   !# normalization factor for subtraction
          end if
          if (jx /= jav .or. it == 0) then
            do i=1,num
              i2=lout+i-1
              ArrB(1:Nd,i2)=0.d0
            end do
            if (Kv == 3) rewind(18)
            do j=1,NumJ
              select case(Kv)
              case(4)
                n=J_n(j)
                k=J_k(j)
                t=J_t(j)
              case(3)
                read(18) nj,k,n,t
              end select
              do i=1,num
                proceed=(lin-1)*Iconverge(i)==0
                if (proceed) then
                  i1=lin+i-1
                  i2=lout+i-1
                  ArrB(n,i2)=ArrB(n,i2)+t*ArrB(k,i1)
                  if (n /= k) ArrB(k,i2)=ArrB(k,i2)+t*ArrB(n,i1)
                end if
              end do
            end do
            if (Kv == 4) then
              do i=lout,lout+num-1
                if (mype==0) then
                  call MPI_Reduce(MPI_IN_PLACE, ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                      MPI_COMM_WORLD, mpierr)
                else
                  call MPI_Reduce(ArrB(1:Nd,i), ArrB(1:Nd,i), Nd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                      MPI_COMM_WORLD, mpierr)
                end if
              end do
            end if
            err1=0.d0
            do i=1,num
              proceed=(lin-1)*Iconverge(i)==0
              if (proceed) then
                i1=lin+i-1
                i2=lout+i-1
                sj=0.0             !# <X1|J**2|X1> is used to check convergence
                do n=1,Nd                       
                  sj=sj+ArrB(n,i1)*ArrB(n,i2) 
                end do
                aj=dsqrt(1+sj) - 1               !# eigenvalue for 2*J
                err=dabs(aj-jav)                 !# difference from target
                err1=dmax1(err1,err)
              end if
            end do
            if (err1 < trsd) then             !# convergency check
              ierr=0
              exit
            end if
            if (it /= 0) then !# eliminating jx component
              do i=1,num
                proceed=(lin-1)*Iconverge(i)==0
                if (proceed) then
                  i1=lin+i-1
                  i2=lout+i-1
                  s=0.d0
                  do n=1,Nd                        !# X1 = X1 - f*Xj
                    a=ArrB(n,i1)-f*ArrB(n,i2)      !## f=1/(2j*(2j+2))
                    s=s+a*a
                    ArrB(n,i1)=a
                  end do
                  s1=1.d0/dsqrt(s)
                  ArrB(1:Nd,i1)=ArrB(1:Nd,i1)*s1  !# normalization of X1
                end if
              end do
            end if
          end if
        end do
        if (Kv == 3) close(18)
        write( *,'(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)') Njd+1-it,err1,ierr
        write(11,'(4X,"iter =",I3,"; max error =",E12.3,"; ierr =",I2)') Njd+1-it,err1,ierr
       return
    end subroutine Prj_J

    subroutine Vread (b,j)
      implicit none
      integer :: i, j
      real(dp), dimension(Nd) :: b
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      b(1:Nd)=ArrB(1:Nd,j)
      return
    end subroutine Vread

    subroutine Vwrite (b,j)
      implicit none
      integer :: i, j
      real(dp), dimension(Nd) :: b
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ArrB(1:Nd,j)=b(1:Nd)
      return
    end subroutine Vwrite
    

end module davidson
