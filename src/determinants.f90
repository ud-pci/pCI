Module determinants
    !
    ! This module implements subroutines related to determinants.
    !
    Use conf_variables

    Implicit None

    Private

    Public :: calcNd0, Dinit, Jterm, Ndet, Pdet, Wdet, Rdet, Rspq, Rspq_phase1, Rspq_phase2
    Public :: Gdet, CompC, CompD, CompD2, CompCD, CompNRC, FormBarr
    Public :: print_bits, convert_bit_rep_to_int_rep, convert_int_rep_to_bit_rep, compare_bit_dets, get_det_indexes

    Integer, Parameter, Public :: bits_per_int = 32
    Integer, Public :: num_ints_bit_rep
    Integer, Dimension(:), Allocatable, Public :: bdet, bdet1, bdet2
    Integer, Dimension(:,:), Allocatable, Public :: Barr

  Contains

    Subroutine calcNd0(ic1, n2)
        Implicit None
        Integer, Intent(out) :: ic1, n2
        Integer :: ic, n1, n0

        ic1=0
        n0=MaxNd0+1
        n1=0

        Do ic=1,Nc4
            if (ic > Nc) Exit
            n1=n1+Ndc(ic)
            If (n1 < n0) Then
                n2=n1
                ic1=ic
            End If
        End Do    
        
    End Subroutine calcNd0
    
    Subroutine Dinit
        Implicit None
        Integer :: i0, ni, j, n, ic, i, nmin, nem, imax
        Logical :: NRorb

        i0=0
        Do ni=1,Ns
            nem=Jj(ni)+1
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                i0=i0+1
            End Do
        End Do

        Allocate(Jz(i0),Nh(i0),Nh0(i0),Nq0(Nsp))
  
        i0=0
        Do ni=1,Ns
            if (ni.GT.1) then
                NRorb=Nn(ni-1).EQ.Nn(ni).AND.Ll(ni-1).EQ.Ll(ni)
            else
                NRorb=.FALSE.
            end if
            nem=Jj(ni)+1 
            Nf0(ni)=i0
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                i0=i0+1
                Jz(i0)=j-nem ! total magnetic quantum number m_j (unit of 1/2)
                Nh(i0)=ni    ! orbital index for Jz 
                if (NRorb) then
                    Nh0(i0)=ni-1
                else
                    Nh0(i0)=ni
                end if
            End Do
        End Do
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        If (nmin < Nsp) Then
            Do ni=nmin,Nsp
                i=i+1
                Nq0(ni)=n
                n=n+Nq(ni)
                If (n >= Ne) Then
                    ic=ic+1
                    Nvc(ic)=i
                    Nc0(ic)=Nso+i0
                    i0=i0+i
                    n=0
                    i=0
                End If
            End Do
        End If
        Return
    End Subroutine Dinit

    Subroutine Jterm
        Implicit None

        Integer         :: j, mt, n, im, ndi, nd1, imax, iconf, ndi1, iconf1, jmax
        Real(dp)        :: d
        Integer, allocatable, dimension(:) :: idet, nmj
        logical :: fin, swapped

        allocate(idet(Ne),Ndc(Nc),Jtc(Nc))

        Njd=0
        Nd=0
        jmax=0

        ! Calculate total number of determinants Nd
        Do iconf1=1,Nc
            ndi=0
            iconf=iconf1
            fin=.true.
            Call Ndet(iconf,fin,idet)
            Do While (.not. fin)
                im = (M-Mj)/2+1
                If (im > jmax) jmax = im
                If (M == Mj) Then
                    ndi = ndi + 1
                End If
                Call Ndet(iconf,fin,idet)
            End Do
            Nd=Nd+Ndi
        End Do

        ! Allocate global arrays:
        ! Iarr - 2d array of determinants, dimension(Ne, Nd)
        ! Jt - array of J term
        ! Njt - array of number of determinants with J term
        if (.not. allocated(Iarr)) allocate(Iarr(Ne,Nd))
        if (.not. allocated(Jt)) allocate(Jt(jmax))
        if (.not. allocated(Njt)) allocate(Njt(jmax))
        if (.not. allocated(nmj)) allocate(nmj(jmax))

        Njd=0
        Nd=0

        ! Loop over configurations
        Do iconf1=1,Nc
            ndi=0
            ndi1=0
            iconf=iconf1
            imax=1
            fin=.true.
            Call Ndet(iconf,fin,idet)

            ! Construct next determinant in configuration iconf
            Do While (.not. fin)
                If (M >= Mj) Then
                    If (M == Mj) ndi=ndi+1
                    If (M == Mj+2) ndi1=ndi1+1

                    nd1=nd+ndi

                    If (M == Mj) Call Pdet(nd1,idet) ! Save determinant to Iarr(1:Ne,nd1)

                    im=(M-Mj)/2+1
                    If (im > imax) imax=im
                End If
                Call Ndet(iconf,fin,idet)
            End Do

            Ndc(iconf)=ndi
            Jtc(iconf)=ndi
            If (kv == 1) Jtc(iconf)=ndi-ndi1
            Nd=Nd+Ndi

            If (ndi /= 0) Then
                nmj(1:imax)=0
                fin= .true.
                Call Ndet(iconf,fin,idet)

                Do While (.NOT.fin) 
                    If (M >= Mj) Then
                        im=(M-Mj)/2+1
                        nmj(im)=nmj(im)+1
                    End If
                    Call Ndet(iconf,fin,idet)
                End Do

                ! list of terms
                Do im=imax, 1, -1
                    n=nmj(im)
                    If (n /= 0) Then
                        mt=(im-1)*2+Mj
                        Do j=1,im
                            nmj(j)=nmj(j)-n
                        End Do

                        ! Check if J term is already in Jt array
                        Do j=1,njd
                            If (Jt(j) == mt) Then
                                njt(j)=njt(j)+n
                                Exit
                            End If
                        End Do
                        
                        ! Add a new J term if not already present
                        If (j > Njd) Then
                            Njd=Njd+1
                            Jt(njd)=mt
                            Njt(njd)=n
                        End If
                    End If
                End Do
            End If
        End Do

        ! Handle case when no determinants were constructed for J terms
        If (Nd == 0) Then
            write( 6,'(/2X,"term J =",F5.1,2X,"is absent in all configurations"/)') Jm
            write(11,'(/2X,"term J =",F5.1,2X,"is absent in all configurations"/)') Jm
            stop
        End If

        ! Write table of number of determinants with respective J terms
        write( 6,'(4X,"Nd   =",I9/3X,23("=")/4X,"N",3X,"  mult.",7X,"J"/3X,23("-"))') Nd
        write(11,'(4X,"Nd   =",I9/3X,23("=")/4X,"N",3X,"  mult.",7X,"J"/3X,23("-"))') Nd

        ! Sort Jt and Njt
        Do
            swapped = .false.
            Do j=2,njd
                If (Jt(j) < Jt(j-1)) Then
                    ! swap Jt values
                    n=Jt(j)
                    Jt(j)=Jt(j-1)
                    Jt(j-1)=n

                    ! swap Njt values
                    n=Njt(j)
                    Njt(j)=Njt(j-1)
                    Njt(j-1)=n
                    
                    swapped = .true.
                End If
            End Do

            ! Exit loop if no elements were swapped
            If (.not. swapped) Exit
        End Do

        ! Print number of determinants for each J term
        Do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            write( 6,'(I5,3X,I8,F9.1)') n,Njt(n),d
            write(11,'(I5,3X,I8,F9.1)') n,Njt(n),d
        End Do
        write( 6,'(3X,23("="))')
        write(11,'(3X,23("="))')

        Deallocate(idet, Jtc, Nq0, nmj)

        Return
    End Subroutine Jterm

    Subroutine Ndet(ic,fin,idet)
        ! This subroutine constructs the next determinant from the list for configuration ic.
        ! fin=TRUE If no determinants in the list are left.
        !
        Implicit None
        Integer :: ic, jm, jf0, j, nqj, nj, nmin, jq, i2, j0, k, i1, n, is, &
                   n3, iq, im, i0, if0, n0, i, nqi, ni, n2, n1
        Integer, allocatable, dimension(:) :: idet
        logical :: fin

        M=0
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)

        If (fin) Then
            Do ni=n1,n2
                nqi=Nq(ni)
                If (nqi /= 0) Then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    if0=Nf0(i)-n0
                    i0 =n0+1
                    im =n0+nqi
                    Do iq=i0,im
                       idet(iq)=if0+iq
                    End Do
                End If
            End Do
            fin=.false.
            Else
            n3=n1+n2
            Do is=n1,n2
                ni =n3-is
                nqi=Nq(ni)
                If (nqi /= 0) Then
                    i  =Nip(ni)
                    n0 =Nq0(ni)
                    i0 =n0+1
                    im =n0+nqi
                    n=Jj(i)+1+Nf0(i)-im
                    i1=im+i0
                    Do k=i0,im
                        iq=i1-k
                        If (idet(iq) < n+iq) Then
                            idet(iq)=idet(iq)+1
                            j0=iq+1
                            If (j0 <= im) Then
                                i2=idet(iq)-iq
                                Do jq=j0,im
                                    idet(jq)=jq+i2
                                End Do
                            End If
                            nmin=ni+1
                            If (nmin <= n2) Then
                                Do nj=nmin,n2
                                    nqj=Nq(nj)
                                    If (nqj /= 0) Then
                                        j=Nip(nj)
                                        n0=Nq0(nj)
                                        jf0=Nf0(j)-n0
                                        j0=n0+1
                                        jm=n0+nqj
                                        Do jq=j0,jm
                                            idet(jq)=jq+jf0
                                        End Do
                                    End If
                                End Do
                            End If
                            fin=.false.
                            Do iq=1,Ne
                               i=idet(iq)
                               m=m+Jz(i)
                            End Do
                            Return
                        End If
                    End Do
                End If
            End Do
            fin=.true.
            Return
        End If
        ! selection of determinants with particular MJ
        Do iq=1,Ne
            i=idet(iq)
            m=m+Jz(i)
        End Do
        Return
    End Subroutine Ndet

    Subroutine Pdet(n,idet)
        Implicit None
        Integer  ::  n
        Integer, allocatable, dimension(:) :: idet
        !  - - - - - - - - - - - - - - - - - - - - - - - - -
        Iarr(1:Ne,n)=idet(1:Ne) 
        Return
    End Subroutine Pdet

    Subroutine Wdet(str)
        ! This subroutine writes the basis set of determinants to the file 'str'.
        !
        Implicit None
        Integer  :: i
        Character(Len=*) :: str

        Open (16,file=str,status='UNKNOWN',form='UNFORMATTED')
        Write(16) Nd,Nsu
        Do i=1,Nd
           write(16) Iarr(1:Ne,i)
        End Do
        close(16)
        Return
    End Subroutine Wdet

    Subroutine Rdet(str)
        ! This subroutine reads the basis set of determinants from the file CONF.DET.
        !
        Implicit None
        Integer :: i
        Character(Len=*) :: str

        If (.not. Allocated(Iarr)) Allocate(Iarr(Ne,Nd))
        Open (16,file=str,status='OLD',form='UNFORMATTED')
        Read(16) Nd,Nsu
        Do i=1,Nd
           Read(16,end=710) Iarr(1:Ne,i)
        End Do
        close(16)
        Return
710     write(*,*)' Rdet: end of CONF.DET for idet=',i
        stop
    End Subroutine Rdet

    Subroutine Rspq(id1,id2,is,nf,i1,j1,i2,j2)
        ! this subroutine compares determinants and counts number of differences in orbitals
        !
        Implicit None
        Integer  :: i, n, ic, is, ni, nj, i2, j2, nf, j, l0, l1, l2, nn0, nn1, &
                    ll0, ll1, jj0, jj1, ndi, k, iconf, i1, j1
        Integer, Allocatable, dimension(:)   :: id1, id2
        Character(Len=1), Dimension(5) :: let
        data let/'s','p','d','f','g'/

        is=1
        ni=0
        nj=0
        i2=0
        j2=0
        nf=3
        i=1
        j=1
        l0=0
  
        Do While (i <= Ne .and. j <= Ne)
            l1=id1(i) ! id1(i) is the i-th element of determinant 1
            l2=id2(j) ! id1(j) is the j-th element of determinant 2
            If (l1 < l0) Then   !### diagnostics for det1:
                l1=Nh(l1)
                l0=Nh(l0)
                nn1=Nn(l1)
                nn0=Nn(l0)
                ll1=Ll(l1)
                ll0=Ll(l0)
                jj1=Jj(l1)
                jj0=Jj(l0)
                ndi=0
                Do ic=1,Nc
                    ndi=ndi+Ndc(ic)
                    iconf=ic
                    If (ndi >= Ndr) Then
                        write (*,'(1X,"RSPQ: Wrong order of shells ",I2,A1,I2,"/2", &
                                " and ",I2,A1,I2,"/2 in configuration ",I5)') &
                                nn0,let(ll0+1),jj0,nn1,let(ll1+1),jj1,iconf
                        write (*,'(4X,"Det",I1,": ",15I4)') 1,(Id1(n),n=1,Ne)
                        write (*,'(4X,"Det",I1,": ",15I4)') 2,(Id2(n),n=1,Ne)
                        Stop
                    End If
                End Do
            End If
            l0=l1
            If (l1 == l2) Then
                i=i+1
                j=j+1
            Else If (l1 > l2) Then
                nj=nj+1
                If (nj > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                j1=j2
                j2=j
                j=j+1
            Else
                ni=ni+1
                If (ni > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                i1=i2
                i2=i
                i=i+1
            End If
        End Do
  
        If (i > j) Then
            nf=ni
            Do k=j,Ne
                j1=j2
                j2=k
            End Do
        Else
            nf=nj
            Do k=i,Ne
                i1=i2
                i2=k
            End Do
        End If
  
        Select Case(nf)
            Case(1) ! number of differences = 1
                k=iabs(j2-i2)
                If (k /= 2*(k/2)) is=-is
                i2=id1(i2)
                j2=id2(j2)
            Case(2) ! number of differences = 2
                k=iabs(j2-i2)
                If (k /= 2*(k/2)) is=-is
                k=iabs(j1-i1)
                If (k /= 2*(k/2)) is=-is
                i2=id1(i2)
                j2=id2(j2)
                i1=id1(i1)
                j1=id2(j1)
        End Select
        Return
    End Subroutine Rspq
    
    Subroutine Rspq_phase1(id1,id2,is,nf,i,j)
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
        Implicit None
        Integer, allocatable, dimension(:), intent(InOut) :: id1, id2
        Integer, intent(InOut)                            :: is, nf
        Integer, intent(InOut)                            :: i(3), j(3)
        
        Integer                                           :: l0, l1, l2, ni, nj, n, ic, nn0, &
                                                             nn1, ll0, ll1, jj0, jj1, ndi, iconf
        Character(Len=1), Dimension(5) :: let
        data let/'s','p','d','f','g'/

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

        l0=0
  
        Do While (i(1) <= Ne .and. j(1) <= Ne)
            l1 = id1(i(1)) ! id1(i) is the i-th element of determinant 1
            l2 = id2(j(1)) ! id2(j) is the j-th element of determinant 2
            If (l1 < l0) Then   !### diagnostics for det1:
                l1=Nh(l1)
                l0=Nh(l0)
                nn1=Nn(l1)
                nn0=Nn(l0)
                ll1=Ll(l1)
                ll0=Ll(l0)
                jj1=Jj(l1)
                jj0=Jj(l0)
                ndi=0
                Do ic=1,Nc
                    ndi=ndi+Ndc(ic)
                    iconf=ic
                    If (ndi >= Ndr) Then
                        write (*,'(1X,"RSPQ: Wrong order of shells ",I2,A1,I2,"/2", &
                                " and ",I2,A1,I2,"/2 in configuration ",I5)') &
                                nn0,let(ll0+1),jj0,nn1,let(ll1+1),jj1,iconf
                        write (*,'(4X,"Det",I1,": ",15I4)') 1,(Id1(n),n=1,Ne)
                        write (*,'(4X,"Det",I1,": ",15I4)') 2,(Id2(n),n=1,Ne)
                        Stop
                    End If
                End Do
            End If
            l0=l1
            If (l1 == l2) Then
                i(1) = i(1) + 1
                j(1) = j(1) + 1
            Else If (l1 > l2) Then
                nj = nj + 1
                If (nj > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                j(3) = j(2)
                j(2) = j(1)
                j(1) = j(1) + 1
            Else
                ni = ni + 1
                If (ni > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                i(3) = i(2)
                i(2) = i(1)
                i(1) = i(1) + 1
            End If
        End Do
  
        If (i(1) > j(1)) Then
            nf = ni
            ! Correct j1/j2 to be the tail end indices; j is the index beyond the last compared:
            Select Case(Ne - j(1))
                Case (0)
                    j(3) = j(2)
                    j(2) = Ne
                Case (1)
                    j(3) = Ne - 1
                    j(2) = Ne
            End Select
        Else
            nf = nj
            ! Correct i1/i2 to be the tail end indices; i is the index beyond the last compared:
            Select Case(Ne - i(1))
                Case (0)
                    i(3) = i(2)
                    i(2) = Ne
                Case (1)
                    i(3) = Ne - 1
                    i(2) = Ne
            End Select
        End If
    End Subroutine Rspq_phase1

    Subroutine Rspq_phase2(id1,id2,is,nf,i,j)
        ! this subroutine follows after Rspq_phase1() to determine the correct
        ! differing indices for the two determinants
        Implicit None
        Integer, allocatable, dimension(:), intent(InOut) :: id1, id2
        Integer, intent(InOut)                            :: is, nf
        Integer, intent(InOut)                            :: i(3), j(3)
        
        Integer                                           :: k
  
        Select Case(nf)
            Case (1) ! number of differences = 1
                k = iabs(j(2) - i(2))
                If (mod(k,2) == 1) is = -is
                i(2) = id1(i(2))
                j(2) = id2(j(2))
            Case (2) ! number of differences = 2
                k = iabs(j(2) - i(2))
                If (mod(k,2) == 1) is = -is
                k = iabs(j(3) - i(3))
                If (mod(k,2) == 1) is = -is
                i(2) = id1(i(2))
                j(2) = id2(j(2))
                i(3) = id1(i(3))
                j(3) = id2(j(3))
        End Select
    End Subroutine Rspq_phase2

    Subroutine Gdet(n,idet)
        ! this subroutine generates the determinant with index n
        Implicit None
        Integer :: n
        Integer, allocatable, dimension(:)  :: idet

        idet(1:Ne)=Iarr(1:Ne,n)
        Return
    End Subroutine Gdet

    Subroutine CompCD(idet1,idet2,icomp)
        Implicit None
        Integer, Intent(Out)  :: icomp
        Integer, allocatable, dimension(:), Intent(In)  :: idet1, idet2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        iconf1(1:Ne)=Nh(idet1(1:Ne))   ! iconf1(i) = No of the orbital occupied by the electron i
        iconf2(1:Ne)=Nh(idet2(1:Ne))    
        Call CompD(iconf1,iconf2,icomp)
        Return
    End Subroutine CompCD

    Subroutine CompD2(id1,id2,nf)
        ! this subroutine compares determinants and counts number of differences in orbitals
        ! return the number of differences nf
        Implicit None
        Integer  :: nf, i, imax
        Integer, allocatable, dimension(:)   :: id1, id2
        Integer, dimension(Nsu) :: det1, det2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        det1=0
        det2=0
        det1(id1(1:Ne))=1
        det2(id2(1:Ne))=1
        nf = 0
        imax=max(maxval(id1),maxval(id2))
        !print*,imax
        Do i=1,imax
            nf = nf + popcnt(xor(det1(i),det1(i))) 
        End Do
        !Do i=1,Ne
        !    nf = nf + popcnt(xor(id1(i),id2(i)))
        !End Do
        !print*,'before',id1(1:Ne),id2(1:Ne),nf
        nf = ishft(nf,-1)
        !print*,'after',id1(1:Ne),id2(1:Ne),nf
        Return
    End Subroutine CompD2


    Subroutine CompD(id1,id2,nf)
        ! this subroutine compares determinants and counts number of differences in orbitals
        ! return the number of differences nf
        Implicit None
        Integer  :: ni, nj, nf, i, j, l1, l2
        Integer, allocatable, dimension(:)   :: id1, id2

        ni=0
        nj=0
        nf=3
        i=1
        j=1
  
        Do While (i <= Ne .and. j <= Ne)
            l1=id1(i) ! id1(i) is the i-th element of determinant 1 (number of orbital occupied)
            l2=id2(j) ! id1(j) is the j-th element of determinant 2 (number of orbital occupied)
            If (l1 == l2) Then
                i=i+1
                j=j+1
            Else If (l1 > l2) Then
                nj=nj+1
                If (nj > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                j=j+1
            Else
                ni=ni+1
                If (ni > 2) Return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                i=i+1
            End If
        End Do
        nf = max(ni,nj)
        Return
    End Subroutine CompD

    Subroutine CompC(idet1,idet2,icomp)
        ! this subroutine compares configurations and counts number of differences in orbitals
        ! return the number of differences nf
        Implicit None
        Integer  :: i1, i2, j1, j2, icomp, is
        Integer, allocatable, dimension(:)  :: idet1, idet2

        Jdel=0
        iconf1(1:Ne)=Nh(idet1(1:Ne))   !### iconf1(i) = No of the orbital occupied by the electron i
        iconf2(1:Ne)=Nh(idet2(1:Ne))    
        Call Rspq(iconf1,iconf2,is,icomp,i1,j1,i2,j2)
        If (icomp == 1) then
            Jdel=iabs(Jj(i2)-Jj(j2))/2
        End If
        Return
    End Subroutine CompC

    Subroutine CompNRC (idet1,idet2,icomp) 
        ! compares two determinants and determines the difference between
        ! corresponding non-relativistic configurations
        Implicit None
        Integer, Intent(Out) :: icomp
        Integer, Allocatable, Dimension(:) :: idet1, idet2

        iconf1(1:Ne)=Nh0(idet1(1:Ne))   !### iconf1(i) = No of the NR orbital occupied by the electron i
        iconf2(1:Ne)=Nh0(idet2(1:Ne))    
        Call CompD(iconf1,iconf2,icomp)
        Return
    End Subroutine CompNRC

    Subroutine FormBarr
        Implicit None
        
        Integer :: i
        Integer, Dimension(:), Allocatable :: idet, bdet

        print*, 'Forming Barr...'

        ! Calculate number of integers needed to store basis orbitals in bit representation
        num_ints_bit_rep = (Nst + bits_per_int - 1) / bits_per_int

        ! Allocate bit_rep equivalent of Iarr
        Allocate(Barr(num_ints_bit_rep, Nd), bdet(num_ints_bit_rep), idet(Ne))

        ! convert integer representation to bit representation
        Do i=1,Nd
            idet = Iarr(1:Ne, i)
            Call convert_int_rep_to_bit_rep(idet, bdet, Ne)
            Barr(1:num_ints_bit_rep, i) = bdet
        End Do

        Deallocate(bdet, idet)

    End Subroutine FormBarr

    Function print_bits(num, n_bits) Result(bitstr)
        ! print bit string from bit representation of determinant
        Implicit None
        Integer, Intent(In) :: num, n_bits
        Character(Len=n_bits) :: bitstr
        Integer :: i

        Do i=1,n_bits
            If (btest(num, i-1)) Then
                bitstr(i:i) = '1'
            Else
                bitstr(i:i) = '0'
            End If
        End Do

    End Function print_bits

    Subroutine convert_int_rep_to_bit_rep(int_rep, bit_rep, size_int_rep)
        Implicit none
        Integer, Dimension(:), Allocatable, Intent(In) :: int_rep ! integer representation
        Integer, Dimension(:), Allocatable, Intent(InOut) :: bit_rep ! bit representation
        Integer, Intent(In) :: size_int_rep
        Integer :: i, index, bit_position ! loop index, integer index, bit position

        ! initialize bit representation array to 0 (all orbitals unoccupied)
        bit_rep = 0

        ! loop over occupied orbitals and set the corresponding bits
        Do i=1, size_int_rep
            ! find which integer the orbital falls in
            index = (int_rep(i) - 1) / bits_per_int + 1
            
            ! find which bit position within the integer the orbital falls in
            bit_position = mod(int_rep(i) - 1, bits_per_int)
            
            ! set the bit corresponding to the orbital
            bit_rep(index) = ibset(bit_rep(index), bit_position)
        End Do

    End Subroutine convert_int_rep_to_bit_rep

    Subroutine convert_bit_rep_to_int_rep(bit_rep, int_rep)
        Implicit None
        Integer, Dimension(:), Allocatable, Intent(In) :: bit_rep ! bit representation
        Integer, Dimension(:), Allocatable :: int_rep             ! integer representation
        Integer :: i, j, cnt                                      ! loop index, integer index, bit position
        Character(len=bits_per_int) :: bit_int

        ! initialize int representation array to 0
        int_rep = 0

        ! loop over occupied orbitals and set the corresponding integer positions
        cnt = 1
        Do i=1, size(bit_rep)
            ! we can skip '0' elements of the bit_rep since they correspond to no occupancy 
            If (bit_rep(i) == 0) Cycle
            bit_int = print_bits(bit_rep(i), bits_per_int)
            Do j = 1, bits_per_int
                If (bit_int(j:j) == '1') Then
                    int_rep(cnt) = j + (i-1)*bits_per_int
                    cnt = cnt + 1
                End If
            End Do
        End Do

    End Subroutine convert_bit_rep_to_int_rep

    Function compare_bit_dets(bdet1, bdet2, n_ints) Result(ndiffs)
        ! print bit string from bit representation of determinant
        Implicit None
        Integer, Dimension(:), Allocatable, Intent(In) :: bdet1, bdet2
        Integer :: i, n_ints, nf, ndiffs

        nf=0

        Do i=1,n_ints
            ! we can skip the i-th element of bdet1 and bdet2 if there are no occupancies
            If (bdet1(i) /= 0 .or. bdet2(i) /= 0) Then

                ! count differing bits
                nf = nf + popcnt(ieor(bdet1(i), bdet2(i)))

                ! we only care if the number of differences between determinants is < 3
                If (nf >= 6) Exit
            End If
        End Do

        ndiffs = (nf + 1) / 2

    End Function compare_bit_dets

    Subroutine get_det_indexes(bdet1, bdet2, n_ints, ndiffs, is, det1_indexes, det2_indexes)
        Implicit None

        Integer, Dimension(:), Allocatable, Intent(In) :: bdet1, bdet2
        Integer, Intent(In) :: n_ints, ndiffs
        Integer, Intent(Out) :: is
        Integer, Dimension(3), Intent(InOut) :: det1_indexes, det2_indexes
        Integer, Dimension(2) :: l1, l2
        Integer :: i, j, i1, j1, bit_position, temp, num_zero_bits, l
        Logical :: first_found, bt1, bt2

        i1=1
        j1=1
        det1_indexes = 0
        det2_indexes = 0

        num_zero_bits = 0
        first_found = .false.

        ! loop through integers representing bitstring determinants
        Do i=1,n_ints
            ! skip the i-th element of bdet1 and bdet2 if there are no occupancies
            If (bdet1(i) == 0 .and. bdet2(i) == 0) Cycle

            Do j=0,bits_per_int-1
                bt1 = btest(bdet1(i), j)
                bt2 = btest(bdet2(i), j)

                ! mark first occupancy
                If (bt1 .or. bt2) first_found = .true.

                ! count zero bits after first occupancy
                If (.not. bt1 .and. .not. bt2) Then
                    if (first_found) num_zero_bits = num_zero_bits + 1
                Else If (bt1 .neqv. bt2) Then
                    bit_position = (32 * (i-1)) + j + 1

                    If (bt1) Then
                        i1 = i1 + 1
                        det1_indexes(i1) = bit_position
                        l1(i1 - 1) = bit_position - num_zero_bits
                    end If

                    If (bt2) Then
                        j1 = j1 + 1
                        det2_indexes(j1) = bit_position
                        l2(j1 - 1) = bit_position - num_zero_bits
                    End If

                    if ((ndiffs == 1 .and. i1 == 2 .and. j1 == 2) .or. (ndiffs == 2 .and. i1 == 3 .and. j1 == 3)) exit
                End If
            End Do
        End Do

        is = 1
        Select Case(ndiffs)
            Case(1)
                If (num_zero_bits > 0) then
                    l = iabs(l2(1) - l1(1) - 1)
                    If (mod(l,2) == 1) is = -is 
                End If
            Case(2)
                l = iabs(l2(1) - l1(1) - mod(num_zero_bits,2))
                If (mod(l,2) == 1) is = -is 
                l = iabs(l2(2) - l1(2) - mod(num_zero_bits,2))
                If (mod(l,2) == 1) is = -is
                temp = det2_indexes(2)
                det2_indexes(2) = det2_indexes(3)
                det2_indexes(3) = temp
                temp = det1_indexes(2)
                det1_indexes(2) = det1_indexes(3)
                det1_indexes(3) = temp
        End Select

    End Subroutine get_det_indexes

End Module determinants