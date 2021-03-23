Module formj2
    ! subroutines FormJ, F_J2, Plj, & J_av (J**2 on determinants)
    Use conf_variables

    Implicit None

  Contains

    Real(dp) Function F_J0(idet)
        Implicit None
    
        Integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is
        Integer, allocatable, dimension(:) :: idet
        Real(dp)  :: t
    
        t=0.d0
        If (nf == 0) Then !determinants are equal
            t=mj*mj
            Do iq=1,Ne
                ia=idet(iq)
                na=Nh(ia)
                ja=Jj(na)
                ma=Jz(ia)
                t=t+ja*(ja+2)-ma**2
                jq0=iq+1
                If (jq0 <= Ne) Then
                    Do jq=jq0,Ne
                        ib=idet(jq)
                        t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
                    End Do
                End If
            End Do
        End If
        F_J0=t
        Return
    End Function F_J0

    Real(dp) Function F_J2(idet1,idet2)
        Use determinants, Only : Rspq
        Implicit None

        Integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is
        Integer, allocatable, dimension(:) :: idet1, idet2
        Real(dp)  :: t

        t=0.d0
        Call Rspq(idet1,idet2,is,nf,ia,ic,ib,id)
        If (nf == 0) Then !determinants are equal
            t=mj*mj
            Do iq=1,Ne
                ia=idet1(iq)
                na=Nh(ia)
                ja=Jj(na)
                ma=Jz(ia)
                t=t+ja*(ja+2)-ma**2
                jq0=iq+1
                If (jq0 <= Ne) Then
                    Do jq=jq0,Ne
                        ib=idet1(jq)
                        t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
                    End Do
                End If
            End Do
        End If
        If (nf == 2) Then !determinants differ by two functions
            t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
                Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
        End If
        F_J2=t
        Return
    End Function F_J2

    Real(dp) Function F_J2_new(idet1,idet2, is, nf, i2, i1, j2, j1)
        Use determinants, Only : Rspq
        Implicit None

        Integer :: ia, na, ja, ma, jq0, iq, ib, ic, id, nf, jq, is, i1, i2, j1, j2
        Integer, allocatable, dimension(:) :: idet1, idet2
        Real(dp)  :: t

        t=0.d0
        If (nf == 0) Then !determinants are equal
            t=mj*mj
            Do iq=1,Ne
                ia=idet1(iq)
                na=Nh(ia)
                ja=Jj(na)
                ma=Jz(ia)
                t=t+ja*(ja+2)-ma**2
                jq0=iq+1
                If (jq0 <= Ne) Then
                    Do jq=jq0,Ne
                        ib=idet1(jq)
                        t=t-Plj(ia,ib)**2-Plj(ib,ia)**2
                    End Do
                End If
            End Do
        End If
        If (nf == 2) Then !determinants differ by two functions
            t=is*(Plj(ia,ic)*Plj(id,ib)+Plj(ic,ia)*Plj(ib,id)- &
                Plj(ia,id)*Plj(ic,ib)-Plj(id,ia)*Plj(ib,ic))
        End If
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        F_J2_new=t
        Return
    End Function F_J2_new

    Real(dp) Function Plj(ia,ib)
        Implicit None

        Integer :: ia, ib, na, nb, ma, mb, ja
        Real(dp) :: t

        t=0.d0
        na=Nh(ia)
        nb=Nh(ib)
        If (na == nb) Then
            ma=Jz(ia)
            mb=Jz(ib)
            If (ma == mb+2) Then
                ja=Jj(na)
                t=dsqrt(ja*(ja+2)-ma*mb+0.d0)
            End If
        End If
        Plj=t
        Return
    End Function Plj

    Recursive Integer Function triStart(mype,npes)
        ! This function divides Nc configurations into equal workloads
        Implicit None
    
        Integer :: mype, npes
        Real :: t
    
        If (mype == 0) Then
            t = 1.0
        Else If (mype == npes) Then
             t = Nc+1.0
        Else
            t = sqrt((1.0/(npes-1.0))*(Nc**2-1.0)+triStart(mype-1,npes)**2)
        End If
        triStart = int(t)
        Return
    End Function triStart

    Subroutine FormJ(mype, npes)
        Use mpi
        Use str_fmt, Only : FormattedTime
        Use vaccumulator
        Use determinants, Only : Gdet, Gdet_win, Rspq, Rspq_phase1, Rspq_phase2
        Use mpi_wins
        Implicit None

        Integer :: k1, n1, ic1, ndn, numjj, i, ij4, n, k, nn, kk, io, counter1, counter2, counter3, diff
        Real :: ttime, ttot
        Real(dp) :: t
        Integer(kind=int64) :: size8, ij8, stot, etot, s1, e1, clock_rate, jstart, jend, memsum
        Integer(kind=int64), dimension(npes) :: ij8s
        Integer, allocatable, dimension(:) :: idet1, idet2, nk
        Integer :: npes, mype, mpierr, interval, remainder, startNc, endNc, sizeNc, counter, win
        Type(IVAccumulator)   :: iva1, iva2
        Type(RVAccumulator)   :: rva1
        Integer, Parameter    :: vaGrowBy = 1000000
        Integer :: iSign, iIndexes(3), jIndexes(3)
        Integer :: fh
        Integer(kind=MPI_OFFSET_KIND) :: disp 
        Character(Len=16) :: filename, timeStr, memStr

        Call system_clock(count_rate=clock_rate)
        If (mype==0) Call system_clock(stot)
        Allocate(idet1(Ne),idet2(Ne),nk(Nc))
        
        ij8=0_int64
        NumJ=0_int64
        If (Kl == 1) Then ! If continuing from previous calculation, skip forming CONF.JJJ
          Open (unit=18,file='CONF.JJJ',status='OLD',form='UNFORMATTED',iostat=io)
          If (io /= 0) Write( 6,'(4X,"Forming matrix J**2")')
          Do
            Read (18,iostat=io) numjj,k,n,t
            If (io > 0) Then 
              print*, 'I/O ERROR at FormJ'
            Else If (io < 0) Then
              Write(11,'(4X,"NumJ =",I12)') NumJ
              Write( *,'(4X,"NumJ =",I12)') NumJ
              Close(unit=18)
              Exit
            Else
              NumJ=numjj
            End If
          End Do
        Else ! Else start new calculation
          If (mype == 0) Then
            Write( 6,'(4X,"Forming matrix J**2")')
            Open (unit=18,file='CONF.JJJ',status='UNKNOWN',form='UNFORMATTED')
            Close(unit=18,status='DELETE')
            Open(unit=18,file='CONF.JJJ',status='NEW',form='UNFORMATTED')
          End If
   
          If (Kv==4) Then
            counter1=1
            counter2=1
            counter3=1
   
            Call IVAccumulatorInit(iva1, vaGrowBy)
            Call IVAccumulatorInit(iva2, vaGrowBy)
   
            Call CreateWindow(win, mpierr)
   
            Call system_clock(s1)
   
            interval = Nc/npes
            remainder = mod(Nc,npes)
   
            startNc = mype*interval+1
            endNc = (mype+1)*interval
            If (mype == npes-1) Then
              endNc=Nc
            End If
            Do ic1=mype+1,Nc,npes
              ndn=Ndc(ic1)
              n=sum(Ndc(1:ic1-1))
              Do n1=1,ndn
                n=n+1
                ndr=n
                Call Gdet_win(n,idet1)
                k=n-n1
                Do k1=1,n1
                  k=k+1
                  Call Gdet_win(k,idet2)
                  Call Rspq_phase1(idet1, idet2, iSign, diff, iIndexes, jIndexes)
                  If (diff == 0 .or. diff == 2) Then
                    nn=n
                    kk=k
                    Call IVAccumulatorAdd(iva1, nn)
                    Call IVAccumulatorAdd(iva2, kk)
                  End If
                End Do
              End Do
            End Do
   
            Call IVAccumulatorCopy(iva1, Jsq%n, counter1)
            Call IVAccumulatorCopy(iva2, Jsq%k, counter2)
       
            Call IVAccumulatorReset(iva1)
            Call IVAccumulatorReset(iva2)
   
            Allocate(Jsq%t(counter1))
            Do n=1,counter1
                nn=Jsq%n(n)
                kk=Jsq%k(n)
                Call Gdet_win(nn,idet1)
                Call Gdet_win(kk,idet2)
                t=F_J2(idet1, idet2)
                Jsq%t(n)=t
            End Do

            Jsq%n = PACK(Jsq%n, Jsq%t/=0)
            Jsq%k = PACK(Jsq%k, Jsq%t/=0)
            Jsq%t = PACK(Jsq%t, Jsq%t/=0)

            Call system_clock(e1)
            ttime=Real((e1-s1)/clock_rate)
            Call FormattedTime(ttime, timeStr)
            Write(*,'(2X,A,1X,I3,1X,A,I9)'), 'core', mype, 'took '// trim(timeStr)// ' for ij8=', counter1
            Call MPI_Barrier(MPI_COMM_WORLD, mpierr)
   
            ij8s=0_int64
            ij8=size(Jsq%t)
            Call MPI_AllReduce(ij8, NumJ, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpierr)
            Call MPI_AllReduce(ij8, maxJcore, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
            Call MPI_AllGather(ij8, 1, MPI_INTEGER8, ij8s, 1, MPI_INTEGER8, MPI_COMM_WORLD, mpierr)
            ij8J = ij8

            Call CloseWindow(win, mpierr)
   
          Else If (Kv==3) Then
            If (mype == 0) Then
              Call system_clock(s1)
              n=0
              Do ic1=1,Nc
                ndn=Ndc(ic1)
                Do n1=1,ndn
                  n=n+1
                  ndr=n
                  Call Gdet(n,idet1)
                  k=n-n1
                  Do k1=1,n1
                    k=k+1
                    Call Gdet(k,idet2)
                    t=F_J2(idet1,idet2)
                    If (t /= 0.d0) Then
                      NumJ=NumJ+1
                      Write(18) NumJ,k,n,t
                    End If
                  End Do
                End Do
              End Do
            End If
          End If    
   
          If (mype == 0) Then  
            Call system_clock(e1)
            ttime=Real((e1-s1)/clock_rate)
             Call FormattedTime(ttime, timeStr) 
            Write(*,'(2X,A)'), 'TIMING >>> Writing CONF.JJJ took '// trim(timeStr) // ' to complete'
            Write(11,'(4X,"NumJ =",I12)') NumJ
            Write( *,'(4X,"NumJ =",I12)') NumJ
            Close(unit=18)
          End If
        End If

        If (mype == 0) Then
            Call system_clock(etot)
            ttot=Real((etot-stot)/clock_rate)
            Call FormattedTime(ttot, timeStr)
            write(*,'(2X,A)'), 'TIMING >>> FormJ took '// trim(timeStr) // ' to complete'
        End If
        Deallocate(idet1,idet2)
        Return
    End Subroutine FormJ
    
    Subroutine J_av(X1,nx,xj,ierr,mype,npes)    !# <x1|J**2|x1>
        Use mpi
        Use determinants, Only : Gdet
        Implicit None

        Integer :: ierr, i, k, n, nx, mpierr
        Integer, optional :: mype, npes
        Integer*8 :: mi, nj
        Real(dp) :: r, t, xj, xj2
        Real(dp), dimension(nx) :: X1
        Integer, allocatable, dimension(:) :: idet1, idet2

        Allocate(idet1(Ne),idet2(Ne))
        ierr=0
        xj=0.d0
        Select Case(Kv)
            Case(4)
                nj=ij8J
            Case(3)
                Open(unit=18,file='CONF.JJJ',status='OLD',form='unformatted')
                nj=NumJ
        End Select
    
        Do i=1,nj
            Select Case(Kv)
            Case(4)
                n=Jsq%n(i)
                k=Jsq%k(i)
                t=Jsq%t(i)
            Case(3)
                Read(18) mi,k,n,t
            End Select
            If (max(k,n) <= nx) Then
                r=t*X1(k)*X1(n)
                If (n /= k) r=r+r
                xj=xj+r
            Else
                Exit
            End If
        End Do
        ! MPI Reduce sum all xj to master core here 
        If (Present(mype)) Call MPI_AllReduce(xj, xj, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        xj=0.5d0*(dsqrt(1.d0+xj)-1.d0)
        If (Kv == 3) Close(unit=18)
        If (K_prj == 1) Then
            If (dabs(xj-XJ_av) > 1.d-1) ierr=1
        End If
        Deallocate(idet1,idet2)
        Return
    End Subroutine J_av

End Module formj2