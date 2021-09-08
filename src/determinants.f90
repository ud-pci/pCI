Module determinants
    !
    ! This module implements Subroutines related to determinants.
    !
    Use conf_variables

    Implicit None

    Private

    Public :: FormD, Dinit, Jterm, Ndet, Pdet, Wdet, Rdet, Rspq, Rspq_phase1, Rspq_phase2
    Public :: Gdet, CompC, CompD, CompD2, CompCD

  Contains
    
    Subroutine FormD
        ! This subroutine evaluates HF energies of the determinants (for HintS).
        Implicit None
        Integer :: i, ie, ke, n
        Real(dp) :: Emin, Emax, x
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Allocate(Diag(Nd))

        Emin= 1.d99
        Emax=-1.d99
        Do n=1,Nd
            x=0.d0
            Do i=1,Ne
                ie=Iarr(i,n)
                ke=Nh(ie)
                x=x+Eps(ke)
            End Do
            x=E_0-x
            Diag(n)=x
            Emin=dmin1(Emin,x)
            Emax=dmax1(Emax,x)
        End Do
        write( 6,'(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)') Emin,Emax
        write(11,'(4X,"FormD: min(E0-Ek)=",F12.6," max(E0-Ek)=",F12.6)') Emin,Emax
        Return
    End Subroutine FormD

    Subroutine Dinit
        Implicit None
        Integer :: i0, ni, j, n, ic, i, nmin, nem, imax
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=0
        Do ni=1,Ns
            nem=Jj(ni)+1
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                i0=i0+1
            End Do
        End Do

        Allocate(Jz(i0),Nh(i0),Nq0(Nsp))
  
        i0=0
        Do ni=1,Ns
            nem=Jj(ni)+1
            Nf0(ni)=i0
            imax=2*Jj(ni)+1
            Do j=1,imax,2
                i0=i0+1
                Jz(i0)=j-nem
                Nh(i0)=ni
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

        Integer         :: ndj, i, j, mt, n, im, ndi, nd1, imax, iconf, ndi1, iconf1
        Real(dp)        :: d
        Integer, allocatable, dimension(:) :: idet
        Integer, dimension(60) :: nmj
        logical :: fin
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(idet(Ne),Ndc(Nc),Jtc(Nc))
        Njd=0
        Nd=0
        ! Calculate Nd
        Do iconf1=1,Nc
            ndi=0
            iconf=iconf1
            fin=.true.
            Call Ndet(iconf,fin,idet)
            Do While (.not. fin)
                If (M == Mj) ndi = ndi + 1
                Call Ndet(iconf,fin,idet)
            End Do
            Nd=Nd+Ndi
        End Do

        if (.not. allocated(Iarr)) allocate(Iarr(Ne,Nd))

        Njd=0
        Nd=0
        Do iconf1=1,Nc
            ndi=0
            ndi1=0
            iconf=iconf1
            imax=1
            fin=.true.
            Call Ndet(iconf,fin,idet)
            Do While (.not. fin)
                If (M >= Mj) Then
                    If (M == Mj) ndi=ndi+1
                    If (M == Mj+2) ndi1=ndi1+1
                    nd1=nd+ndi
                    If (M == Mj) Then
                      Call Pdet (nd1,idet)
                    End If
                    im=(M-Mj)/2+1
                    If (im > imax) imax=im
                End If
                Call Ndet(iconf,fin,idet)
            End Do
    !  -- -  - - - - - - - - - - - - - - - - - - - - - -
            Ndc(iconf)=ndi
            Jtc(iconf)=ndi
            If (kv == 1) Jtc(iconf)=ndi-ndi1
            Nd=Nd+Ndi
            If (ndi /= 0) Then
                If (imax > 60) Then
                    write(*,*) ' Jterm: array nmj is too small'
                    write(*,*) ' Dimension is 60, required:',imax
                    Read(*,*)
                End If
                Do im=1,imax
                    nmj(im)=0
                End Do
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
                im=imax+1
    230         im=im-1
                If (im >= 1) Then
                    n=nmj(im)
                    If (n /= 0) Then
                        mt=(im-1)*2+Mj
                        Do j=1,im
                            nmj(j)=nmj(j)-n
                        End Do
                        If (njd /= 0) Then
                            Do j=1,njd
                                If (Jt(j) == mt) Then
                                    njt(j)=njt(j)+n
                                    goto 230
                                End If
                            End Do
                        End If
                        njd=njd+1
                        If (njd > IPjd) Then
                            write(*,*) ' Jterm: number of Js'
                            write(*,*) ' exceed IPjd=',IPjd
                            Read(*,*)
                        End If
                        Jt(njd)=mt
                        Njt(njd)=n
                    End If
                    goto 230
                End If
            End If
        End Do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        If (Nd == 0) Then
            write( 6,'(/2X,"term J =",F5.1,2X,"is absent in all configurations"/)') Jm
            write(11,'(/2X,"term J =",F5.1,2X,"is absent in all configurations"/)') Jm
            stop
        End If

        ! Writing table of J
        write( 6,'(4X,"Nd   =",I9/3X,23("=")/4X,"N",3X,"  mult.",7X,"J"/3X,23("-"))') Nd
        write(11,'(4X,"Nd   =",I9/3X,23("=")/4X,"N",3X,"  mult.",7X,"J"/3X,23("-"))') Nd
        i = 1
        Do While (i == 1)
            i=0
            Do j=2,njd
                If (Jt(j) < Jt(j-1)) Then
                    n=Jt(j)
                    Jt(j)=Jt(j-1)
                    Jt(j-1)=n
                    n=Njt(j)
                    Njt(j)=Njt(j-1)
                    Njt(j-1)=n
                    i=1
                End If
            End Do
        End Do
        Do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            write( 6,'(I5,3X,I8,F9.1)') n,Njt(n),d
            write(11,'(I5,3X,I8,F9.1)') n,Njt(n),d
        End Do
        write( 6,'(3X,23("="))')
        write(11,'(3X,23("="))')
     
        i=0
        Do j=1,njd
            d=Jt(j)*0.5d0
            n=j
            ndj=Njt(n)
        End Do

        Deallocate(idet, Nq, Nip, Jtc, Nq0)

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
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        M=0
        n1=Nc0(ic)+1
        n2=Nc0(ic)+Nvc(ic)
        ! - - - - - - - - - - - - - - - -
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
        Integer  :: i, n
        Character(Len=*) :: str
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Open (16,file=str,status='UNKNOWN',form='UNFORMATTED')
        Write(16) Nd,Nsu
        Do i=1,Nd
           write(16) Iarr(1:Ne,i)
        End Do
        close(16)
        Return
    End Subroutine Wdet

    Subroutine Rdet
        ! This subroutine reads the basis set of determinants from the file CONF.DET.
        !
        Implicit None
        Integer :: i, n
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Open (16,file='CONF.DET',status='OLD',form='UNFORMATTED')
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
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
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
        ! The two Integer index vectors, "i" and "j", are three elements in size and are
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
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
  
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
        Integer :: i, n
        Integer, allocatable, dimension(:)  :: idet
        ! - - - - - - - - - - - - - - - - - - - - - - - -
        idet(1:Ne)=Iarr(1:Ne,n)
        Return
    End Subroutine Gdet

    Subroutine CompCD(idet1,idet2,icomp)
        Implicit None
        Integer  :: i, i1, i2, icomp, m
        Integer, allocatable, dimension(:)  :: idet1, idet2
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
        Integer  :: ni, nj, nf, i, j, l1, l2, k, imax
        Integer, allocatable, dimension(:)   :: id1, id2
        Integer, dimension(IP6) :: det1, det2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        det1=0
        det2=0
        Do i=1,Ne
            det1(id1(i))=1
            det2(id2(i))=1
        End Do
        nf = 0
        imax=max(maxval(id1),maxval(id2))
        print*,imax
        Do i=1,imax
            nf = nf + popcnt(xor(det1(i),det1(i))) 
        End Do
        !Do i=1,Ne
        !    nf = nf + popcnt(xor(id1(i),id2(i)))
        !End Do
        print*,'before',id1(1:Ne),id2(1:Ne),nf
        nf = ishft(nf,-1)
        print*,'after',id1(1:Ne),id2(1:Ne),nf
        Return
    End Subroutine CompD2


    Subroutine CompD(id1,id2,nf)
        ! this subroutine compares determinants and counts number of differences in orbitals
        ! return the number of differences nf
        Implicit None
        Integer  :: ni, nj, nf, i, j, l1, l2, k
        Integer, allocatable, dimension(:)   :: id1, id2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
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
        Integer  :: i, i1, i2, j1, j2, icomp, is, m
        Integer, allocatable, dimension(:)  :: idet1, idet2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        Jdel=0
        iconf1(1:Ne)=Nh(idet1(1:Ne))   !### iconf1(i) = No of the orbital occupied
        iconf2(1:Ne)=Nh(idet2(1:Ne))    !#### by the electron i
        Call Rspq(iconf1,iconf2,is,icomp,i1,j1,i2,j2)
        If (icomp == 1) then
          Jdel=iabs(Jj(i2)-Jj(j2))/2
        End If
        Return
    End Subroutine CompC

end module determinants