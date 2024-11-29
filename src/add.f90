Program add
    Use params, nlvl => Nlv, nnparam => Nn
    Implicit None

    ! set initial size and growth size for storing max number of configurations in Ac, NOz
    Integer, Parameter :: growth_size = 1000

    Integer :: Ncor, NsvNR, Mult, l, keyPT, max_num_shells
    Character(Len=32) :: strfmt
    Real(dp), Allocatable, Dimension(:, :) :: Ac
    Integer, Allocatable, Dimension(:) :: NOz, Knr
    Real(dp), Allocatable, Dimension(:) :: V, Nqmin, Nqmax

    ! Print program header
    strfmt = '(/4X,"PROGRAM add v3.0"/)'
    Write(6, strfmt)

    ! Read inputs from ADD.INP
    Call Input 

    ! Construct new non-relativistic configurations from excitations
    Call Constr

    ! Test parity
    Call Parity
    
    ! Expand non-relativistic configurations into relativistic configurations
    Call Expand

    ! Choose CI or CI+PT
    Call CI_or_PT(Ncpt,keyPT)

    ! Form CONF.INP
    Call PrintConfINP(Ncpt,keyPT)     ! Forms CONF.INP

Contains

    Subroutine Input
        Implicit None
        Integer :: nsmc, ic, ne0, i, i1, i2, iz, j, nx, ny, nz, ncheck, ndum, nlvi, nl1, &
                    iskip, ics, ji, k, jk, idif, qqmax, nozmax
        Real(dp) :: x
        Integer, Allocatable, Dimension(:) :: nyi, myi, nyk, myk
        Integer, Allocatable, Dimension(:) :: Nq1, Nq2, Nq3, Nn
        Real(dp), Allocatable, Dimension(:) :: Qnl
        Real(dp), Allocatable, Dimension(:,:) :: Ac0
        Character(Len=1), Allocatable, Dimension(:) :: let, chr
        Character(Len=1), Dimension(4) :: txt
        Character(Len=128) :: string, strfmt
        Logical :: new

        string(1:5) = '     '
        txt(1)=' '
        txt(2)=' '
        txt(3)=' '
        txt(4)=' '

        ! Read in the header of ADD.INP
        Open(unit=10, file='ADD.INP')

        strfmt = '(2(5X,I3))'
        Read(10, strfmt) Ncor    ! number of basic configurations
        Read(10, strfmt) NsvNR   ! number of active non-relativistic shells
        Read(10, strfmt) Mult    ! multiplicity of excitations
        Read(10, strfmt) Ne      ! number of valence electrons

        Write(*,*) ' Ncor=',Ncor,' NsvNR=',NsvNR
        Write(*,*) ' Mult=',Mult,' Ne   =',Ne

        nsmc=0
        iskip=0     ! - skipped rel. configurations
        ics=0       ! - basic non-rel. configurations 

        Allocate(Qnl(40)) ! max number of shells in a configuration
        Allocate(Nn(NsvNR), Nq1(NsvNR), Nq2(NsvNR), Nq3(NsvNR), let(NsvNR), chr(NsvNR))
        Allocate(NOz(growth_size))

        ! Run through each core configuration
        Do ic=1,Ncor
            ne0=0
            i1=1
            iz=0 ! number of shells
            
            Do While (ne0 < Ne)
                i2=i1+5

                ! Read core configuration
                strfmt = '(4A1,A)'
                Read(10, strfmt) (txt(i),i=1,4), string

                ! If configuration is written in readable form
                If (txt(1) == 'L') Then 
                    strfmt = '(6(I2,A1,I2,1X))'
                    Read(string, strfmt) (Nq1(i),chr(i),Nq2(i),i=i1,i2)
                    Do i=i1,i2
                        Call ConvertChar(chr(i),Nq3(i),0)
                        Qnl(i)=(1000*Nq1(i)+100*Nq3(i)+Nq2(i))/10000.d0
                    End Do
                ! If configuration is written in digital form
                Else
                    strfmt = '(F7.4,5(4X,F7.4))'
                    Read(string, strfmt) (Qnl(i),i=i1,i2)
                End If

                ! Loop through each shell in configuration
                Do i=i1,i2
                    x=abs(Qnl(i))+1.d-7
                    If (x.LT.1.d-6) Exit
                    nx=10000*x
                    ny=100*x
                    nz=(nx-100*ny)
                    If (nz.EQ.0) Then
                        strfmt = '(" no electrons in ",I2," shell"," of configuration ",I5)'
                        Write(6,strfmt) iz,ic
                        Stop
                    End If
                    ne0=ne0+nz
                    iz=iz+1
                End Do

                If (ne0 > Ne) Then
                    Write(*,*) ' Too many electrons for ic =', ic
                    Stop
                End If

                i1=iz+1
            End Do

            If (.not. Allocated(Ac0)) Allocate(Ac0(Ncor, iz+Mult))
            Noz(ic)=iz          !### number of shells in config-n
            do j=1,iz
                Ac0(ic,j)=Qnl(j)
            end do
        End Do

        Deallocate(Qnl)

        nozmax = maxval(NOz)+1
        If (.not. Allocated(nyi)) Allocate(nyi(nozmax))
        If (.not. Allocated(myi)) Allocate(myi(nozmax))
        If (.not. Allocated(nyk)) Allocate(nyk(nozmax))
        If (.not. Allocated(myk)) Allocate(myk(nozmax))

        Do ic=1,Ncor
            call Squeeze(ic,Ac0,ji,nyi,myi) 
            new=.TRUE.
            outer: do k=1,ics
                call Squeeze(k,Ac0,jk,nyk,myk)
                if (ji.EQ.jk) then   ! compare two non.rel. config-s
                    inner: do j=1,ji
                        idif = iabs(nyi(j)-nyk(j)) + iabs(myi(j)-myk(j))
                        if (idif.GT.0) Cycle outer
                    end do inner
                    new=.FALSE.
                    Exit outer
                end if
            end do outer
        
            if (new) then
                ics=ics+1
                Noz(ics)=ji
                nsmc=max(nsmc,ji) !### max number of shells in config-n
                do j=1,ji
                    Ac0(ics,j)=(100.d0*nyi(j)+myi(j))/10000.d0
                end do
                strfmt = '(I3,1X,F7.4,5(4X,F7.4),/6(4X,F7.4))'
                write(*,strfmt) ics,(Ac0(ics,j),j=1,ji)
            else
                iskip=iskip+1
            end if         
        End Do

        Deallocate(nyi, myi, nyk, myk)

        max_num_shells = nsmc+2*Mult
        Allocate(Ac(growth_size, max_num_shells))
        Do ic=1,Ncor
            Do j=1,nsmc
                Ac(ic,j) = Ac0(ic,j)
            End Do
        End Do
        Deallocate(Ac0)

        If (iskip.GT.0) Then
            write(*,*) iskip,' config-s skipped from ',Ncor
            Ncor=Ncor-iskip
        End If

        ! Read any empty lines before list of active non-relativistic shells
        strfmt = '(6(2x,I2,A1,1x,I2,1x,I2))'
        ndum=0
        Do While (ndum < 1)
            Read(10, strfmt) ndum
        End Do
        Backspace(10)

        Read(10, strfmt) (Nn(i),let(i),Nq1(i),Nq2(i), i=1,NsvNR)
        Write(*, strfmt) (Nn(i),let(i),Nq1(i),Nq2(i), i=1,NsvNR)

        ! Read list of active non-relativistic shells
        Do i=1,NsvNR
            Call ConvertChar(let(i),Ll(i),-1)

            If (Ll(i).LT.0) Then
                Write(*,*) 'Orbital "',let(i),'" for i =',i,' is undefined'
                Read(*,*)
                Stop
            End If
        End Do

        qqmax = 10*Nn(NsvNR)+Ll(NsvNR)
        allocate(Nqmin(qqmax), Nqmax(qqmax), V(NsvNR))

        Do i=1,NsvNR
            nlvi=10*Nn(i)+Ll(i)
            Nqmin(nlvi)=Nq1(i)
            Nqmax(nlvi)=Nq2(i)

            V(i) = (10000*Nn(i) + 1000*Ll(i))/100000. + 5.0E-8
        End Do

        Do j=1,NsvNR
            nl1=100*(abs(V(j))+0.001)
            Nn(j)=nl1/10
            Ll(j)=nl1-10*Nn(j)
        End Do

        Nc=Ncor

        Return
    End Subroutine Input

    Subroutine ConvertChar(char,ichar,init) 
        ! Convert s,p,d -> 0,1,2
        Implicit None
        Integer :: init, ichar
        Character(Len=1) :: char
        ichar=init
        If (char.EQ.'s'.OR.char.EQ.'S') ichar = 0
        If (char.EQ.'p'.OR.char.EQ.'P') ichar = 1
        If (char.EQ.'d'.OR.char.EQ.'D') ichar = 2
        If (char.EQ.'f'.OR.char.EQ.'F') ichar = 3
        If (char.EQ.'g'.OR.char.EQ.'G') ichar = 4
        If (char.EQ.'h'.OR.char.EQ.'H') ichar = 5
        If (char.EQ.'i'.OR.char.EQ.'I') ichar = 6
        If (char.EQ.'k'.OR.char.EQ.'K') ichar = 7
        If (char.EQ.'l'.OR.char.EQ.'L') ichar = 8
        If (char.EQ.'m'.OR.char.EQ.'M') ichar = 9
        Return
    End Subroutine ConvertChar

    Subroutine Expand
        ! Expand the list of non-relativistic configurations into a list of relativistic configurations
        Implicit None
        Integer :: ic, icnr, kc, ji, k, n, jk, i, nx, mx, mmax1, mmax2, mmin1, m1, k1, ivar, irmax, k1max
        Integer :: j, ir, j1, l, iq1, iq2, idif, ncr, jkmax, nozmax, num_rel_confs, num_rel_confs0
        Real(dp) :: qnl
        Integer, Allocatable, Dimension(:)  :: nyi, myi, nyk, myk, ivc, ni, li, iv, NozN, temp
        Integer, Allocatable, Dimension(:, :)   :: mx1, mx2, ivv
        Real(dp), Allocatable, Dimension(:, :) :: AcN, temp2

        Write(*,*) ' forming relativistic conf-s...'
        Write(*,*) ' Expanding ',Nc,' configurations'

        ic=0
        icnr=0
        kc=0
        ji=0

        nozmax = maxval(NOz)+1

        Allocate(nyi(nozmax), myi(nozmax))
        Allocate(nyk(nozmax), myk(nozmax))
        Allocate(ivc(nozmax), ni(nozmax), li(nozmax), iv(nozmax))

        ! Pre-calculate array sizes
        num_rel_confs=0 ! number of relativistic configurations
        irmax=0         ! max number of rel. configurations in a non-relativistic configuration
        k1max=0         ! max number of j=l-1/2 to j=l+1/2 shells

        Do ic=1,Nc
            ! Process the current configuration
            Call Squeeze(ic,Ac,jk,nyk,myk)

            ! Number of relativistic cofiguration in the current non-relativistic configuration
            ir=1 
            Do i=1,jk
                nx=nyk(i)
                ni(i)=nx/10             ! quantum number n
                li(i)=nx-10*ni(i)       ! quantum number l
              
                If (li(i).EQ.0) Then
                    iv(i)=1             ! variants of occupations
                Else
                    mx=myk(i)
                    mmax1=2*li(i)
                    mmax2=mmax1+2
                    mmin1=max(0,mx-mmax2)
                    m1=min(mx,mmax1)    ! max occupation for j=l-1/2 shell
                    iv(i)=m1-mmin1+1    ! variants of occupations
                    Do k=m1,mmin1,-1
                        k1=m1-k+1
                    End Do
                    if (k1 > k1max) k1max = k1
                End If
                ir=ir*iv(i)
            End Do
            if (ir > irmax) irmax = ir
            num_rel_confs = num_rel_confs + ir
        End Do

        If (.not. Allocated(Knr)) Allocate(Knr(num_rel_confs))
        If (.not. Allocated(NozN)) Allocate(NozN(num_rel_confs))
        If (.not. Allocated(AcN)) Allocate(AcN(num_rel_confs, max_num_shells))
        If (.not. Allocated(ivv)) Allocate(ivv(irmax, nozmax))
        If (.not. Allocated(mx1)) Allocate(mx1(nozmax, k1max))
        If (.not. Allocated(mx2)) Allocate(mx2(nozmax, k1max))

        ! Reallocate Ac, NOz
        num_rel_confs0 = size(NOz)
        Call move_alloc(Ac, temp2)
        Call move_alloc(NOz, temp)
        Allocate(NOz(num_rel_confs), Ac(num_rel_confs, max_num_shells))
        NOz(1:num_rel_confs0) = temp
        Ac(1:num_rel_confs0,1:max_num_shells) = temp2
        Deallocate(temp, temp2)

        ! Loop over non-relativistic configurations
        Do ic=1,Nc
            ! Print message when all core configurations are expanded
            If (ic == Ncor+1) then
                Write (*,*) Ncor,' initial config-s expanded in ',kc,' relat. config-s'
                ncr=kc
            End If

            ! Process the current configuration
            Call Squeeze(ic,Ac,jk,nyk,myk)

            ! Compare two non-relativistic configurations
            If (ji.EQ.jk) Then
                ! Compare element by element and skip the configuration if a difference is found
                Do j=1,ji
                  idif = iabs(nyi(j)-nyk(j)) + iabs(myi(j)-myk(j))
                  If (idif > 0) Exit
                End Do

                ! If all elements match, update the configuration and continue
                ji=jk
                Do i=1,ji
                    nyi(i)=nyk(i)
                    myi(i)=myk(i)
                End Do
                Cycle
            End If

            ! New non-relativistic configuration
            icnr=icnr+1
            Knr(icnr)=kc+1

            ! Number of relativistic cofiguration in the current non-relativistic configuration
            ir=1 
            Do i=1,jk
                nx=nyk(i)
                ni(i)=nx/10             ! quantum number n
                li(i)=nx-10*ni(i)       ! quantum number l
              
                If (li(i).EQ.0) Then
                    iv(i)=1             ! variants of occupations
                    mx1(i,1)=0          ! occupation of j=l-1/2 shell 
                    mx2(i,1)=myk(i)     ! occupation of j=l+1/2 shell 
                Else
                    mx=myk(i)
                    mmax1=2*li(i)
                    mmax2=mmax1+2
                    mmin1=max(0,mx-mmax2)
                    m1=min(mx,mmax1)    ! max occupation for j=l-1/2 shell
                    iv(i)=m1-mmin1+1    ! variants of occupations
                    Do k=m1,mmin1,-1
                        k1=m1-k+1
                        mx1(i,k1)=k     ! occupation of j=l-1/2 shell
                        mx2(i,k1)=mx-k  ! occupation of j=l+1/2 shell
                    End Do
                End If
                ir=ir*iv(i)
            End Do

            ! Check if number of relativistic configurations exceeds ivv array limit
            If (ir.GT.500) Then
                write(*,*) ' number of relativistic configurations ',ir
                write(*,*) ' in one NR configuration exceed array size 500'
                stop
            End If

            ivc(jk)=1
            Do i=jk-1,1,-1
                ivc(i)=ivc(i+1)*iv(i+1)
            End Do
            
            Do i=1,ir
                Do j=1,jk
                    ivar=(i-1)/ivc(j)+1
                    Do While (ivar > iv(j)) 
                        ivar=ivar-iv(j)
                    End Do
                    ivv(i,j)=ivar
                End do
            End do 

            Do i=1,ir
                kc=kc+1
                j1=0
                Do j=1,jk
                    n=ni(j)
                    l=li(j)
                    ivar=ivv(i,j)
                    iq1=mx1(j,ivar)
                    iq2=mx2(j,ivar)
                    If (iq1.GT.0) Then
                        j1=j1+1
                        qnl=-(1000*n+100*l+iq1)/10000.d0
                        AcN(kc,j1)=qnl
                    End If
                    If (iq2.GT.0) Then
                        j1=j1+1
                        qnl= (1000*n+100*l+iq2)/10000.d0
                        AcN(kc,j1)=qnl
                    End If
                End Do
                NozN(kc)=j1
            End Do       
        End Do

        ! Print message when all configurations of proper parity are expanded
        Write (*,*) Nc,' config-s expanded in ',kc,' relat. config-s'
        Nc=kc
        Ncor=ncr
        Do k=1,kc
            Noz(k)=NozN(k)
            Do n=1,Noz(k)
                Ac(k,n)=AcN(k,n)
            End Do
        End Do

        Write(*,*) ' relativistic conf-s: ',Nc

        Deallocate(NozN, AcN, nyi, myi, nyk, myk, ivc, ni, li, iv, mx1, mx2)

        Return
    End Subroutine Expand


    Subroutine Squeeze(ic, Ac, k, nyi, myi)
        ! This subroutine reads config ic from Ac and returns arrays nl and occ. num. for each shell
        !
        ! ic - current config in Ac
        ! k - number of non-rel shells in config
        ! nyi - array of QN n & l
        ! myi - array of occupation numbers
        Implicit None
        Integer :: k, j, j1, ic
        Integer, Allocatable, Dimension(:), Intent(Out) :: nyi, myi
        Integer, Allocatable, Dimension(:) :: nzi, mzi
        Real(dp), Allocatable, Dimension(:,:) :: Ac

        If (.not. Allocated(nzi)) Allocate(nzi(NOz(ic)+1))
        If (.not. Allocated(mzi)) Allocate(mzi(NOz(ic)+1))
        If (.not. Allocated(nyi)) Allocate(nyi(NOz(ic)+1))
        If (.not. Allocated(myi)) Allocate(myi(NOz(ic)+1))

        Do j=1,NOz(ic)
            j1=j
            nzi(j) = 100*abs(Ac(ic,j)) + 0.1
            mzi(j) = 10000*abs(Ac(ic,j)) - 100*nzi(j) + 0.1
        End Do
        nzi(j1+1)=0

        k=0
        j=0

        Do while (j.LT.NOz(ic))
            k=k+1
            j=j+1
            nyi(k)=nzi(j)
            ! If two shells have same quantum numbers n & l
            If (nzi(j).EQ.nzi(j+1)) Then  
                myi(k)=mzi(j)+mzi(j+1)
                j=j+1
            Else  
                myi(k)=mzi(j)
            End If
        End Do

        Deallocate(nzi, mzi)

        Return
    End Subroutine Squeeze

    Subroutine Constr
        Implicit None
        Integer :: Nc0, i, k, j, jj, l, ii, mlast
        Real(dp) :: B
        Integer, Allocatable, Dimension(:) :: temp
        Real(dp), Allocatable, Dimension(:) :: Ac1
        Real(dp), Allocatable, Dimension(:,:) :: temp2
        Character(Len=128) :: strfmt

        If (.not. Allocated(Ac1)) Allocate(Ac1(max_num_shells))
        
        ! Loop over number of excitations
        Do l=1,Mult
            Write(*,*) ' multiplicity of excitation:    ',l

            Nc0=Nc ! number of configurations to allow excitations froms

            ! Loop over currently constructed non-rel. configurations
            Do i=1,Nc0
                ! Loop over number of shells in non-rel. configurations
                Do k=1,NOz(i)
                    ! Loop over number of active non-rel. shells
                    Do jj=1,NsvNR
                        If (k.EQ.1) Then
                          Continue
                        Else
                            Do ii=1,k-1
                              Ac1(ii) = Ac(i,ii)
                            End Do
                        End If
                        If (Ac(i,k).LT.0.) Ac1(k)= Ac(i,k) + 0.0001
                        If (Ac(i,k).GT.0.) Ac1(k)= Ac(i,k) - 0.0001
                        If (NOz(i).GT.k) Then
                            Do j=k+1,NOz(i)
                                Ac1(j)= Ac(i,j)
                            End Do
                        End if
                        mlast= NOz(i)+1
                        B = abs(V(jj)) + 0.0001
                        Ac1(mlast) = sign(B,V(jj))
                        Call Shell(Ac1,mlast)
                        Call Reorder(Ac1,mlast)
                        If (New(Ac1,mlast)) Then
                            Nc=Nc+1
                            If (Nc > size(NOz)) Then
                                Call move_alloc(NOz, temp)
                                Call move_alloc(Ac, temp2)

                                Allocate(NOz(size(temp)+growth_size), Ac(size(temp2,1)+growth_size, max_num_shells))
                                NOz(1:size(temp)) = temp
                                Ac(1:size(temp2,1),1:max_num_shells) = temp2
                                Deallocate(temp, temp2)
                            End If
                            NOz(Nc) = mlast
                            Do j=1,mlast
                                Ac(Nc,j)=Ac1(j)
                            End Do
                        End If
                    End Do
                End Do
            End Do
            Write(*,*) ' number of NR conf-s after Constr: ',Nc
        End Do

        Deallocate(Ac1)

        Return
    End Subroutine Constr

    Subroutine Shell(Ac1, nozi)
        Implicit None
        Integer :: j, nozi, k, nx, ny, iz, nac, nxy
        Real(dp) :: B, R
        Real(dp), Allocatable, Dimension(:) :: Ac1

        Do j=1,nozi-1
            k=j+1
   100      B= Ac1(j) - Ac1(k)
            If (abs(B).GE.0.005) goto 110     ! shells are different
            nx=10000*abs(Ac1(k)) + 0.1
            ny=100*abs(Ac1(k)) + 0.1
            Nxy=nx-100*ny                     ! electrons on shell k
            Ac1(j)= abs(Ac1(j)) + Nxy*0.0001  ! summing electrons
            Ac1(j)= sign(Ac1(j),Ac1(k))       ! restoring the sign
            R= 0.01*ny
            Ac1(k)= sign(R,Ac1(k))            ! now shell k has 0 electrons
   110      If (k.NE.nozi) Then
                k=k+1
                goto 100
            End if
        End Do

        iz=0
        Do j=1,nozi                         ! excluding shells with 0 electrons
            nx=10000*abs(Ac1(j))+0.1
            ny=100*abs(Ac1(j))+0.1
            nac=nx-100*ny
            If (nac.NE.0) Then
                Ac1(j-iz) = Ac1(j)
            Else
                iz=iz+1
            End If
        End Do
        nozi=nozi-iz

        Return
    End Subroutine Shell

    Subroutine Reorder(Ac1, nozi)
        Implicit None
        Integer :: n1, j, m, k, k1, k2, nozi
        Real(dp) :: P1, P2, A
        Real(dp), Allocatable, Dimension(:) :: Ac1
        Character(Len=128) :: strfmt

        strfmt = '(2x,"Ac1(",I4,")=",F10.4," is not in V:",/(6F10.4))'
        n1=NsvNR
        Do j=1,nozi-1
            m=j
            k1=0
            Do k=1,n1
                P1= abs(Ac1(m)-V(k))
                If (P1.LT.0.002) k1=k
                If (k1.NE.0) goto 100
            End Do
            Write(*, strfmt) m,Ac1(j),(V(k),k=1,n1)
            Stop
 100        m=m+1
            k2=0
            Do k=1,n1
                P2= abs(Ac1(m)-V(k))
                If (P2.LT.0.002) k2=k
                If (k2.NE.0) goto 110
            End Do
            Write(*, strfmt) m,Ac1(m),(V(k),k=1,n1)
            Stop
 110        If (k1.LE.k2) goto 120
            A = Ac1(j)
            Ac1(j)=Ac1(m)
            Ac1(m) = A
            k1=k2
 120        If (m.NE.nozi) goto 100
        End Do

        Return
    End Subroutine Reorder

    Logical Function New(Ac1, noz1)
        ! Returns TRUE if configuration is new and valid
        Implicit None
        Integer :: noz1, i, j, i1, nnj, llj, jlj, nq, j1, i2, nnj2, llj2, jlj2, nq2, ierr, nlj, nozi
        Real(dp) :: d, d2, del
        Real(dp), Allocatable, Dimension(:) :: Ac1

        New=.FALSE.
  
        Do j=1,noz1                  !### check of the Pauli principle:
            i1=sign(1.01_dp,Ac1(j))     !#### nq.LE.(4*l+2)
            d=abs(Ac1(j))+1.E-6
            d=10*d
            nnj=d
            d=10*(d-nnj)
            llj=d
            jlj=2*llj+i1
            d=100*(d-llj)
            nq=d
            If (nq.GT.4*llj+2) Return
        End Do

        Do j=1,noz1                  ! does this config.
            i1=sign(1.01_dp,Ac1(j))     !  belong to RAS?
            d=abs(Ac1(j))+1.E-6
            d=10*d
            nnj=d
            d=10*(d-nnj)
            llj=d
            jlj=2*llj+i1
            d=100*(d-llj)
            nq=d

            if (i1.LT.0) return

            nlj=10*nnj+llj
            ! check that occupation numbers are within allowed limits
            If (nq.LT.Nqmin(nlj).OR.nq.GT.Nqmax(nlj)) Return 
        End Do

        ! Check and see if configurations are new
        Outer: Do i=1,Nc            
            nozi=NOz(i)
            If (nozi.EQ.noz1) Then
                Inner: Do j=1,nozi
                    del=abs(Ac(i,j)-Ac1(j))
                    If (del.GT.5.E-5) Cycle Outer
                End Do Inner
                Return
            End If
        End Do Outer
  
        New=.TRUE.

        Return
    End Function New 

    Subroutine Parity
        ! Check parity of configurations
        Implicit None
        Integer :: l1, j, n1j, l1j, n1q, is1, mz, n, ln1, nnj, lnj, nnq, isn
        Real(dp) :: d

        Write(*,*) " test of parity..."

        l1=0
        Do j=1,NOz(1)
            d=abs(Ac(1,j))+1.E-6
            d=10.*d
            n1j=d
            d=10.*(d-n1j)
            l1j=d
            d=100.*(d-l1j)
            n1q=d
            If (n1q.NE.(n1q/2)*2) l1= l1 + l1j
        End Do
        
        is1=0
        If (l1.NE.(l1/2)*2) is1=1
        mz=0

        Do n=2,Nc
            ln1=0
            Do j=1,NOz(n)
                d=ABS(Ac(n,j))+1.E-6
                d=10.*d
                nnj=d
                d=10.*(d-nnj)
                lnj=d
                d=100.*(d-lnj)
                nnq=d
                If (nnq.NE.(nnq/2)*2) ln1= ln1 + lnj
            End Do
            isn=0
            If (ln1.NE.(ln1/2)*2) isn=1
            If (is1.NE.isn) Then
                mz=mz+1
                Nc=Nc-1
            Else
                Do j=1,NOz(n)
                    Ac(n-mz,j)=Ac(n,j)
                End Do
                NOz(n-mz)=NOz(n)
            End If
        End Do

        Write(*,*) " NR conf-s of proper parity: ", Nc

        Return
    End Subroutine Parity

    Subroutine CI_or_PT(Ncpt, keyPT)
        Implicit None
        Integer :: Ncpt, keyPT, k
        Character(Len=1), Dimension(5) :: com
        Character(Len=70) :: string
        Character(Len=128) :: strfmt

        ! last line before head of CONF.INP must start with ">"
        strfmt = '(A)'
        Do While (com(1).NE.">")
            Read (10, strfmt) com(1)      
        End Do

        Ncpt=0
        Cut0=0.d0

 100    Read (10,15) (com(k),k=1,5),string
 15     Format(5A1,A)
        If (com(2).ne." ".and.com(2).ne."-") goto 110
        If (com(3).ne." ".and.com(3).ne."-") goto 110
        If (com(4).ne." ".and.com(4).ne."-") goto 110
        goto 200

 110    If (com(1).ne."n".and.com(1).ne."N") goto 120
        If (com(2).ne."c".and.com(2).ne."C") goto 120
        If (com(3).ne."p".and.com(3).ne."P") goto 120
        If (com(4).ne."t".and.com(4).ne."T") goto 120
        Read (string,*) Ncpt
        goto 100

 120    If (com(1).ne."c".and.com(1).ne."C") goto 100
        If (com(2).ne."u".and.com(2).ne."U") goto 100
        If (com(3).ne."t".and.com(3).ne."T") goto 100
        If (com(4).ne."0") goto 100
        Read (string,*) Cut0
        goto 100

 200    If (Ncpt+Cut0.GT.0) then
            Write(*,*) " Choose variant: CI (0), or PT [Nc -> Ncpt] (1)"
            Read(*,*) keyPT
        Else
            keyPT=0
        End If
        If (keyPT.eq.1) Then
            Ncpt=Nc
            Nc=Ncor
        Else
            Ncpt=0
        End If  

        Return
    End Subroutine CI_or_PT

    Subroutine PrintConfINP(Ncpt, keyPT)
        Implicit None
        Integer :: i, k, Ncpt, keyPT, Ncl, ic, ic1, j, icnr, ne1
        Real(kind=dp) :: Jm
        Real(dp), Allocatable, Dimension(:) :: Q
        Character(Len=1), Dimension(5) :: com
        Character(Len=70) :: string
        Character(Len=128) :: strfmt

        Open(unit=11, file="CONF.INP")
        Open(unit=12, file="CONF_.INP")

        ! last line before head of CONF.INP must start with ">"
        Rewind(10)
        strfmt = '(A)'
        Do While (com(1).NE.">")
            Read (10, strfmt) com(1)      
        End Do

        read (10, strfmt) string
        write(11, strfmt) string
        write(12, strfmt) string
        
        do i=1,4
          read (10,15) (com(k),k=1,5),string
          read (string,*) z
          write(11,"(5A1,F5.1)") (com(k),k=1,5),z
          write(12,"(5A1,F5.1)") (com(k),k=1,5),z
          If (i==4) Jm= z 
        end do
        Mj=2*dabs(Jm)+0.01d0 
 15     format(5A1,A)
        read (10,15) (com(k),k=1,5),string
        if (com(2).ne."n".and.com(2).ne."N") goto 200
        if (com(3).ne."s".and.com(3).ne."S") goto 200
        if (com(4).ne."o".and.com(4).ne."O") goto 200
        read (string,*) Nso
        write(11,25) (com(k),k=1,5),Nso
        write(12,25) (com(k),k=1,5),Nso
 25     format(5A1,I4)
        read (10,15) (com(k),k=1,5)
        if (com(2).ne."n".and.com(2).ne."N") goto 300
        if (com(3).ne."c".and.com(3).ne."C") goto 300
        if (com(4).ne." ") goto 300
        write(11,35) (com(k),k=1,5),Nc
        write(12,35) (com(k),k=1,5),Nc
 35     format(5A1,I7)

 100    read (10,15) (com(k),k=1,5),string
        if (com(2).ne." ".and.com(2).ne."-") goto 110
        if (com(3).ne." ".and.com(3).ne."-") goto 110
        if (com(4).ne." ".and.com(4).ne."-") goto 110
        if (com(2).eq." ".and.com(3).eq." ") then
          backspace(10)
        else
          write(11,15) (com(k),k=1,5),string
          write(12,15) (com(k),k=1,5),string
        end if
        goto 140
 110    if (com(2).ne."n".and.com(2).ne."N") goto 120
        if (com(3).ne."e".and.com(3).ne."E") goto 120
        read (string,*) ne1
        if (Ne.NE.ne1) then
          write(*,*) " Ne=",ne1,"; expected: ",Ne
          read(*,*)
          stop
        end if
 120    if (com(1).ne."n".and.com(1).ne."N") goto 130
        if (com(2).ne."c".and.com(2).ne."C") goto 130
        if (com(3).ne."p".and.com(3).ne."P") goto 130
        if (com(4).eq."t".or.com(4).eq."T") write(11,45) (com(k),k=1,5),Ncpt
 45     format(5A1,I6,"     #")
        goto 100

 130    write(11,15) (com(k),k=1,5),string
        write(12,15) (com(k),k=1,5),string
        goto 100

 140    Allocate(Q(Nso))
        if (Nso.NE.0) then
          read (10,55) (Q(i),i=1,Nso)
        end if
        write(11,55) (Q(i),i=1,Nso) ! for Nso=0 empty line written
        write(12,55) (Q(i),i=1,Nso)
        Deallocate(Q)

 55     format (6(4X,F7.4))
 65     format (I4,F7.4,5(4X,F7.4))

        icnr=1
        kc=Knr(icnr)

        if (keyPT.eq.1) then 
          Ncl=Ncpt
        else
          Ncl=Nc
        end if

        do ic=1,Ncl
          ic1=ic-(ic/10000)*10000
          if (ic.EQ.kc) then    ! New non.rel. config-n
            write(11,65) icnr
            write(12,65) icnr
            icnr=icnr+1
            kc=Knr(icnr)
          end if   
          write (12,'(I4,F7.4,*(4X,F7.4))') ic1,(Ac(ic,j),j=1,NOz(ic))
          if (NOz(ic).LE.6) then
            write (11,65) ic1,(Ac(ic,j),j=1,NOz(ic))
          else
            if (NOz(ic).GT.6.AND.NOz(ic).LE.12) then
              write (11,65) ic1,(Ac(ic,j),j=1,6)
              write (11,55) (Ac(ic,j),j=7,NOz(ic))
            else
              write (*,75) ic,NOz(ic)
 75           format ("NOz(",I5,")=",I3," is greater than 12")
              Stop
            end if
          end if
        end do

        ! Close ADD.INP and CONF.INP
        Close(12)
        Close(11)
        Close(10)
        Write(*,*) " File CONF.INP is formed"

        ! Dellocate Arrays
        Deallocate(NOz, Knr, Ac)

        Return
 200    write (*,*) " Expected Nso=, got: ",com
        read(*,*)
        stop
 300    write (*,*) " Expected Nc =, got: ",com
        read(*,*)
        stop
    End Subroutine PrintConfINP

End Program add