Module determinants
    !
    ! This module implements subroutines that 
    !
    Use conf_variables

    Implicit None

    Private

    Public :: FormD, Dinit, Jterm, Ndet, Pdet, Wdet, Rdet, Rspq, Rspq_phase1, Rspq_phase2
    Public :: Gdet, Gdet_win, CompC, CompD, CompCD

  Contains
    
    subroutine FormD
      ! evaluates HF energies of the determinants (for HintS)
      implicit none
      integer :: i, ie, ke, n
      real(dp) :: Emin, Emax, x
      !- - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(Diag(Nd))
      Emin= 1.d99
      Emax=-1.d99
      do n=1,Nd
        x=0.d0
        do i=1,Ne
          ie=Iarr(i,n)
          ke=Nh(ie)
          x=x+Eps(ke)
        end do
        x=E_0-x
        Diag(n)=x
        Emin=dmin1(Emin,x)
        Emax=dmax1(Emax,x)
      end do
      write( 6,'(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)') Emin,Emax
      write(11,'(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)') Emin,Emax
      return
    end subroutine FormD

    subroutine Dinit
      implicit none
      integer :: i0, ni, j, n, ic, i, nmin, nem, imax
      !- - - - - - - - - - - - - - - - - - - - - - - - -
      i0=0
      do ni=1,Ns
        nem=Jj(ni)+1
        imax=2*Jj(ni)+1
        do j=1,imax,2
          i0=i0+1
        end do
      end do
      allocate(Jz(i0),Nh(i0),Nq0(Nsp))    
      i0=0
      do ni=1,Ns
        nem=Jj(ni)+1
        Nf0(ni)=i0
        imax=2*Jj(ni)+1
        do j=1,imax,2
          i0=i0+1
          Jz(i0)=j-nem
          Nh(i0)=ni
        end do
      end do
      n=0
      ic=0
      i0=0
      i=0
      nmin=Nso+1
      if (nmin < Nsp) then
        do ni=nmin,Nsp
          i=i+1
          Nq0(ni)=n
          n=n+Nq(ni)
          if (n >= Ne) then
            ic=ic+1
            Nvc(ic)=i
            Nc0(ic)=Nso+i0
            i0=i0+i
            n=0
            i=0
          end if
        end do
      end if
      return
    end subroutine Dinit

    subroutine Jterm
        implicit none
        integer      :: ndj, i, j, mt, n, im, ndi, nd1, imax, iconf, ndi1, ic1
        real(dp)       :: d
        integer, allocatable, dimension(:) :: idet
        integer, dimension(6) :: nmj
        logical :: fin
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(idet(Ne),Ndc(Nc),Jtc(Nc))
        Njd=0
        Nd=0
        ! Calculate Nd
        do ic1=1,Nc
            ndi=0
            iconf=ic1
            fin=.true.
    260     call Ndet(iconf,fin,idet)
            if (.not. fin) then
               if (M == Mj) ndi=ndi+1
               goto 260
            end if
            Nd=Nd+Ndi
        end do

        allocate(Iarr(Ne,Nd))

        Njd=0
        Nd=0
         do ic1=1,Nc
            ndi=0
            ndi1=0
            iconf=ic1
            imax=1
            fin=.true.
    210     call Ndet(iconf,fin,idet)
            if (.not. fin) then
               if (M >= Mj) then
                  if (M == Mj) ndi=ndi+1
                  if (M == Mj+2) ndi1=ndi1+1
                  nd1=nd+ndi
                  if (M == Mj) then
                    call Pdet (nd1,idet)
                  end if
                  im=(M-Mj)/2+1
                  if (im > imax) imax=im
               end if
               goto 210
            end if
    !  - - - - - - - - - - - - - - - - - - - - - - - - -
            Ndc(iconf)=ndi
            Jtc(iconf)=ndi
            if (kv == 1) Jtc(iconf)=ndi-ndi1
            Nd=Nd+Ndi
            if (ndi /= 0) then
                if (imax > 60) then
                    write(*,*) ' Jterm: array nmj is too small'
                    write(*,*) ' Dimension is 60, required:',imax
                    read(*,*)
                end if
                do im=1,imax
                    nmj(im)=0
                end do
                fin= .true.
    220         call Ndet(iconf,fin,idet)
                if (.NOT.fin) then
                    if (M >= Mj) then
                        im=(M-Mj)/2+1
                        nmj(im)=nmj(im)+1
                    end if
                    goto 220
                end if
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    !    list of terms
                im=imax+1
    230         im=im-1
                if (im >= 1) then
                    n=nmj(im)
                    if (n /= 0) then
                        mt=(im-1)*2+Mj
                        do j=1,im
                            nmj(j)=nmj(j)-n
                        end do
                        if (njd /= 0) then
                            do j=1,njd
                                if (Jt(j) == mt) then
                                    njt(j)=njt(j)+n
                                    goto 230
                                end if
                            end do
                        end if
                        njd=njd+1
                        if (njd > IPjd) then
                            write(*,*) ' Jterm: number of Js'
                            write(*,*) ' exceed IPjd=',IPjd
                            read(*,*)
                        end if
                        Jt(njd)=mt
                        Njt(njd)=n
                    end if
                    goto 230
                end if
            end if
        end do
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nd == 0) then
            write( 6,5) Jm
            write(11,5) Jm
    5       format(/2X,'term J =',F5.1,2X, &
               'is absent in all configurations'/)
            stop
        end if
        write( 6,15) Nd
        write(11,15) Nd
    15  format(4X,'Nd   =',I9/3X,21('=')/4X,'N',3X, &
            '  mult.',7X,'J'/3X,21('-'))
    240 i=0
        do j=2,njd
            if (Jt(j) < Jt(j-1)) then
                n=Jt(j)
                Jt(j)=Jt(j-1)
                Jt(j-1)=n
                n=Njt(j)
                Njt(j)=Njt(j-1)
                Njt(j-1)=n
                i=1
            end if
        end do
        if (i == 1) goto 240
        do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            write( 6,25) n,Njt(n),d
            write(11,25) n,Njt(n),d
    25      format(I5,3X,I8,F9.1)
        end do
        write( 6,35)
        write(11,35)
    35  format(3X,21('='))
    ! - - - - - - - - - - - - - - - - - - - - - - -
        i=0
        do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            ndj=Njt(n)
        end do
    end subroutine Jterm

    subroutine Ndet(ic,fin,idet)
        implicit none
        integer :: ic, jm, jf0, j, nqj, nj, nmin, jq, i2, j0, k, i1, n, is, &
                   n3, iq, im, i0, if0, n0, i, nqi, ni, n2, n1
        integer, allocatable, dimension(:) :: idet
        logical :: fin
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        !   construction of the next det from the list for configuration IC.
        !   fin=TRUE, if no dets in the list are left.
        !- - - - - - - - - - - - - - - - - - - - - - - - -
        M=0
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)
        ! - - - - - - - - - - - - - - - -
        if (fin) then
            do ni=n1,n2
                nqi=Nq(ni)
                if (nqi /= 0) then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    if0=Nf0(i)-n0
                    i0 =n0+1
                    im =n0+nqi
                    do iq=i0,im
                       idet(iq)=if0+iq
                    end do
                end if
            end do
            fin=.false.
            else
            n3=n1+n2
            do is=n1,n2
                ni =n3-is
                nqi=Nq(ni)
                if (nqi /= 0) then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    i0 =n0+1
                    im =n0+nqi
                    n=Jj(i)+1+Nf0(i)-im
                    i1=im+i0
                    do k=i0,im
                        iq=i1-k
                        if (idet(iq) < n+iq) then
                            idet(iq)=idet(iq)+1
                            j0=iq+1
                            if (j0 <= im) then
                                i2=idet(iq)-iq
                                do jq=j0,im
                                    idet(jq)=jq+i2
                                end do
                            end if
                            nmin=ni+1
                            if (nmin <= n2) then
                                do nj=nmin,n2
                                    nqj=Nq(nj)
                                    if (nqj /= 0) then
                                        j=Nip(nj)
                                        n0=Nq0(nj)
                                        jf0=Nf0(j)-n0
                                        j0=n0+1
                                        jm=n0+nqj
                                        do jq=j0,jm
                                            idet(jq)=jq+jf0
                                        end do
                                    end if
                                end do
                            end if
                            fin=.false.
                            do iq=1,Ne
                               i=idet(iq)
                               m=m+Jz(i)
                            end do
                            return
                        end if
                    end do
                end if
            end do
            fin=.true.
            return
        end if
        !  - - - - - - - - - - - - - - - - - - - - - - - - -
        ! selection of dets with particular MJ
        do iq=1,Ne
           i=idet(iq)
           m=m+Jz(i)
        end do
        return
    end subroutine Ndet

    subroutine Pdet(n,idet)
      implicit none
      integer  ::  n
      integer, allocatable, dimension(:) :: idet
      !  - - - - - - - - - - - - - - - - - - - - - - - - -
      Iarr(1:Ne,n)=idet(1:Ne) 
      return
    end subroutine Pdet

    subroutine Wdet(str)
      implicit none
      integer  :: i, n
      character(len=*) :: str
      open (16,file=str,status='UNKNOWN',form='UNFORMATTED')
      Write(16) Nd,Nsu
      do i=1,Nd
         write(16) Iarr(1:Ne,i)
      end do
      close(16)
      return
    end subroutine Wdet

    subroutine Rdet
      implicit none
      integer :: i, n
      open (16,file='CONF.DET',status='OLD', &
            form='UNFORMATTED')
      read(16) Nd,Nsu
      do i=1,Nd
         read(16,end=710) Iarr(1:Ne,i)
      end do
      close(16)
      return
710   write(*,*)' Rdet: end of CONF.DET for idet=',i
      stop
    end subroutine Rdet

    subroutine Rspq(id1,id2,is,nf,i1,j1,i2,j2)
      ! this subroutine compares determinants and counts number of differences in orbitals
      implicit none
      integer  :: n, is, ni, nj, i2, j2, nf, i, j, l0, l1, l2, nn0, nn1, &
                  ll0, ll1, jj0, jj1, ndi, k, iconf, i1, j1
      integer, allocatable, dimension(:)   :: id1, id2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      is=1
      ni=0
      nj=0
      i2=0
      j2=0
      nf=3
      i=1
      j=1

      do while (i <= Ne .and. j <= Ne)
        l1=id1(i) ! id1(i) is the i-th element of determinant 1
        l2=id2(j) ! id1(j) is the j-th element of determinant 2
        if (l1 == l2) then
          i=i+1
          j=j+1
        else if (l1 > l2) then
          nj=nj+1
          if (nj > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          j1=j2
          j2=j
          j=j+1
        else
          ni=ni+1
          if (ni > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          i1=i2
          i2=i
          i=i+1
        end if
      end do

      if (i > j) then
        nf=ni
        do k=j,Ne
          j1=j2
          j2=k
        end do
      else
        nf=nj
        do k=i,Ne
          i1=i2
          i2=k
        end do
      end if

      select case(nf)
        case(1) ! number of differences = 1
          k=iabs(j2-i2)
          if (k /= 2*(k/2)) is=-is
          i2=id1(i2)
          j2=id2(j2)
        case(2) ! number of differences = 2
          k=iabs(j2-i2)
          if (k /= 2*(k/2)) is=-is
          k=iabs(j1-i1)
          if (k /= 2*(k/2)) is=-is
          i2=id1(i2)
          j2=id2(j2)
          i1=id1(i1)
          j1=id2(j1)
      end select

      return
    end subroutine Rspq
    
    subroutine Rspq_phase1(id1,id2,is,nf,i,j)
      ! this subroutine compares determinants and counts number of differences in orbitals
      ! it does NOT post-process the located indices to determine scaling factor (is) or
      ! correct differing indices
      !
      ! The two integer index vectors, "i" and "j", are three elements in size and are
      ! used as:
      !
      !     i(1)        running index over comparison
      !     i(2)        first differing index (a.k.a. i2)
      !     i(3)        second differing index (a.k.a. i1)
      !
      implicit none
      integer, allocatable, dimension(:), intent(InOut) :: id1, id2
      integer, intent(InOut)                            :: is, nf
      integer, intent(InOut)                            :: i(3), j(3)
      
      integer                                           :: l1, l2, ni, nj
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      is=1
      nf=3
      i(1) = 1
      j(1) = 1
      i(2) = 0
      j(2) = 0
      i(3) = 0
      j(3) = 0
      
      ni=0
      nj=0

      do while (i(1) <= Ne .and. j(1) <= Ne)
        l1 = id1(i(1)) ! id1(i) is the i-th element of determinant 1
        l2 = id2(j(1)) ! id2(j) is the j-th element of determinant 2
        if (l1 == l2) then
          i(1) = i(1) + 1
          j(1) = j(1) + 1
      else if (l1 > l2) then
          nj = nj + 1
          if (nj > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          j(3) = j(2)
          j(2) = j(1)
          j(1) = j(1) + 1
        else
          ni = ni + 1
          if (ni > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          i(3) = i(2)
          i(2) = i(1)
          i(1) = i(1) + 1
      end if
      end do

      if (i(1) > j(1)) then
        nf = ni
        ! Correct j1/j2 to be the tail end indices; j is the index beyond the last compared:
        select case(Ne - j(1))
          case (0)
            j(3) = j(2)
            j(2) = Ne
          case (1)
            j(3) = Ne - 1
            j(2) = Ne
        end select
      else
        nf = nj
        ! Correct i1/i2 to be the tail end indices; i is the index beyond the last compared:
        select case(Ne - i(1))
          case (0)
            i(3) = i(2)
            i(2) = Ne
          case (1)
            i(3) = Ne - 1
            i(2) = Ne
        end select
      end if
    end subroutine Rspq_phase1

    subroutine Rspq_phase2(id1,id2,is,nf,i,j)
      ! this subroutine follows after Rspq_phase1() to determine the correct
      ! differing indices for the two determinants
      implicit none
      integer, allocatable, dimension(:), intent(InOut) :: id1, id2
      integer, intent(InOut)                            :: is, nf
      integer, intent(InOut)                            :: i(3), j(3)
      
      integer                                           :: k
      ! - - - - - - - - - - - - - - - - - - - - - - - - -

      select case(nf)
        case (1) ! number of differences = 1
          k = iabs(j(2) - i(2))
          if (mod(k,2) == 1) is = -is
          i(2) = id1(i(2))
          j(2) = id2(j(2))
        case (2) ! number of differences = 2
          k = iabs(j(2) - i(2))
          if (mod(k,2) == 1) is = -is
          k = iabs(j(3) - i(3))
          if (mod(k,2) == 1) is = -is
          i(2) = id1(i(2))
          j(2) = id2(j(2))
          i(3) = id1(i(3))
          j(3) = id2(j(3))
      end select
    end subroutine Rspq_phase2

    subroutine Gdet(n,idet)
      ! this subroutine generates the determinant with index n
      implicit none
      integer :: i, n
      integer, allocatable, dimension(:)  :: idet
      ! - - - - - - - - - - - - - - - - - - - - - - - -
      idet(1:Ne)=Iarr(1:Ne,n)
      return
    end subroutine Gdet

    subroutine Gdet_win(n,idet)
      ! this subroutine generates the determinant with index n
      implicit none
      integer :: i, n
      integer, allocatable, dimension(:)  :: idet
      ! - - - - - - - - - - - - - - - - - - - - - - - -
      idet(1:Ne)=Iarr_win(1:Ne,n)
      return
    end subroutine Gdet_win

    subroutine CompCD(idet1,idet2,icomp)
      implicit none
      integer  :: i, i1, i2, icomp, m
      integer, allocatable, dimension(:)  :: idet1, idet2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ic1(1:Ne)=Nh(idet1(1:Ne))   ! ic1(i) = No of the orbital occupied by the electron i
      ic2(1:Ne)=Nh(idet2(1:Ne))    
      call CompD(ic1,ic2,icomp)
      return
    end subroutine CompCD

    subroutine CompD(id1,id2,nf)
      ! this subroutine compares determinants and counts number of differences in orbitals
      ! return the number of differences nf
      implicit none
      integer  :: ni, nj, nf, i, j, l1, l2, k
      integer, allocatable, dimension(:)   :: id1, id2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ni=0
      nj=0
      nf=3
      i=1
      j=1

      do while (i <= Ne .and. j <= Ne)
        l1=id1(i) ! id1(i) is the i-th element of determinant 1 (number of orbital occupied)
        l2=id2(j) ! id1(j) is the j-th element of determinant 2 (number of orbital occupied)
        if (l1 == l2) then
          i=i+1
          j=j+1
        else if (l1 > l2) then
          nj=nj+1
          if (nj > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          j=j+1
        else
          ni=ni+1
          if (ni > 2) return ! if difference > 2 then matrix element between id1 and id2 will be 0
          i=i+1
        end if
      end do
      nf = max(ni,nj)
      return
    end subroutine CompD

    subroutine CompC(idet1,idet2,icomp)
      implicit none
      integer  :: i, i1, i2, j1, j2, icomp, is, m
      integer, allocatable, dimension(:)  :: idet1, idet2
      ! - - - - - - - - - - - - - - - - - - - - - - - - -
      ic1(1:Ne)=Nh(idet1(1:Ne))   !### ic1(i) = No of the orbital occupied
      ic2(1:Ne)=Nh(idet2(1:Ne))    !#### by the electron i
      call Rspq(ic1,ic2,is,icomp,i1,j1,i2,j2)
      return
    end subroutine CompC

end module determinants