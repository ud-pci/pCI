module conf_pt_aux

    use conf_pt_variables
    implicit none

    contains

    subroutine Input
        use conf_init, only : inpstr
        ! This subroutine reads in from CONF.INP the following variables
        ! name, Z, Am, XJ_av, Jm, Nso, Nc, Kv, Nlv
        ! K_is, Kbrt, Kout, Kecp, C_is, Gj, Cut0, Ncpt
        ! Qnl - atomic configurations
        ! Nsp, Kl
        implicit none
        character(len=1) :: name(16)
        integer :: i, j, k, i1, i2, ic, icc, ne0, nx, ny, nz, istr
        real(dp)  :: x
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        write( 6,5) 
        write(11,5)
    5   format(/4X,'Program CONF_PT', &
               /4X,'PT corrections to binding energy', & 
               /4X,'Zero approximation is taken from CONF.XIJ', &
               /4X,'New vectors are in CONF_PT.XIJ and', &
               /4X,'new input is in CONF_new.INP')
    ! - input from the file 'CONF.INP' - - - - - - - - - - - - - - - -
        open(unit=10, status='OLD',file='CONF.INP')
        read(10,'(1X,16A1)') name
        write(*,'(4X,16A1)') name
        write(11,'(4X,16A1)') name
        read(10,'(5X,F5.1)') Z, Am, XJ_av, Jm
        read(10,'(5X,I6)') Nso, Nc, Kv, Nlv, Ne
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Qnl(100000000)) ! 1B = upper bound
        IPlv=3*Nlv
        K_is = 0                ! 
        Kbrt = 0                !  
        Kout = 1                !   Optional
        Kecp = 0                !   parameters
        C_is = 0.d0             !      
        Gj   = 0.d0             !      
        Cut0 = 0.001            !       
        Ncpt = 0                !   
        Ksig=0
        Kdsig=0
        average_diag=.true.      ! diagonal is averaged over relat. config.
    100 call inpstr(istr)       !
        if (istr /=  1) goto 100 !
        Ncci = Nc
        if (Ncpt < Ncci) then
            write (*,*) 'Nothing to do for NcPT =', Ncpt, ' & Nc = NcCI =', Ncci
            stop
        end if
        Nc = Ncpt
    
        if (abs(C_is) < 1.d-6) K_is = 0
        if (K_is == 0) C_is = 0.d0
        if (K_is == 2  .or.  K_is == 4) then
            write(*,*) ' SMS to include 1-e (1), 2-e (2), both (3): '
            read(*,*) K_sms
            if ((K_sms-1)*(K_sms-2)*(K_sms-3) /=  0) stop
        end if
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nso /=  0) then
            read(10,75) (Qnl(i),i=1,Nso)
    75      format(6(4X,F7.4))
        end if
    ! - - - reading in configurations - - - - - - - - - - - - - - - - -
        i1 = Nso + 1
        do ic = 1, Nc
            icc=ic
            ne0=0
    200     i2=i1+5
            read(10,75,err=900,end=900) (Qnl(i),i=i1,i2)
            do i = i1,i2
                x = abs(Qnl(i))+1.d-9
                if (x < 1.d-8) goto 210
                nx=10000*x
                ny=100*x
                nz=(nx-100*ny)
                ne0=ne0+nz
            end do
    210     i2=i-1
            i1=i2+1
            if (ne0 < Ne) goto 200
            if (ne0 > Ne) then
                write(6,'(I6)') 'INPUT: too many electrons for ic =', ic
                Stop
            end if
        end do
        Nsp=i2
        close(unit=10)
    ! - - - - - - - - - - - - - - - - - - - - - - - - -  
        Kl = 3
        return
    900 write (*,*) 'Error while reading configuration ', icc
        write (*,*) 'Ncpt =', Ncpt, ' Ncci =', Ncci
        Stop
    end subroutine Input

    subroutine Init
        use conf_init, only : inpstr
        implicit none
        real(dp)                     :: c1, c2, z1, d
        real(dp), dimension(IP6)     :: p, q, p1, q1
        real(dp), dimension(4*IP6)   :: pq
        real(dp), dimension(IPs)     :: Qq1
        integer                    :: ii, ni, if, nj, i, nnj, llj, jlj, kkj, &
                                     nec, imax, j, n, ic, i0, nmin, i1, n2, n1, & l
        integer, dimension(4*IPs)  :: IQN
        integer, dimension(33)     :: nnn, jjj, nqq
        character(len=1)           :: Let(9), lll(33)
        logical                    :: longbasis
        equivalence (IQN(1),PQ(21)),(Qq1(1),PQ(2*IPs+21))
        equivalence (p(1),pq(1)), (q(1),pq(IP6+1)), &
                    (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        data Let/'s','p','d','f','g','h','i','k','l'/
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        c1 = 0.01d0
        mj = 2*abs(Jm)+0.01d0
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        open(12,file='CONF.DAT',status='OLD', &
              access='DIRECT',recl=2*IP6*IPmr,err=700)
        read(12,rec=1) p
        read(12,rec=2) q
        read(12,rec=5) p1
        read(12,rec=6) q1
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        z1 = pq(1)
        if (abs(Z-z1) > 1.d-6) then
           write( 6,5) Z,z1
           write(11,5) Z,z1
    5      format('nuc. charge is changed: Z =',F12.6,' ><',F12.6)
           read(*,*)
        end if
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Nvc(Nc),Nc0(Nc),Nq(Nsp),Nip(Nsp))
        Ns = pq(2)+c1
        ii = pq(3)+c1
        longbasis=abs(PQ(20)-0.98765d0) < 1.d-6
        write( 6,15) Kl,Kv,Z,Jm,Nsp,Ns,Nso,Nc
        write(11,15) Kl,Kv,Z,Jm,Nsp,Ns,Nso,Nc
    15  format (4X,'Kl  =',I3,7X,'Kv  =',I3, &
                    7X,'Z   =',F6.2,4X,'Jm  =',F6.2, &
                   /4X,'Nsp =',I6,4X,'Ns  =',I3,7X,'Nso =',I3, &
                    7X,'Nc =',I6)
        if (longbasis) then
            write( *,*) ' Using variant for long basis '
            write(11,*) ' Using variant for long basis '
            do ni=1,Ns
                Nn(ni)=IQN(4*ni-3)
                Ll(ni)=IQN(4*ni-2)
                Kk(ni)=IQN(4*ni-1)
                Jj(ni)=IQN(4*ni)
            end do
        else
            if=20
            do ni=1,Ns
                if=if+1
                Nn(ni)=pq(if)+c1
                if=if+1
                Ll(ni)=pq(if)+c1
                if=if+3
                c2=dsign(c1,pq(if))
                Kk(ni)=pq(if)+c2
                if=if+1
                c2=dsign(c1,pq(if))
                Jj(ni)=pq(if)+c2
            end do
        end if
        Nsu=0
        do nj=1,Nsp
            i=sign(1.d0,Qnl(nj))
            d=abs(Qnl(nj))+1.d-14
            d=10.0*d
            nnj=d
            d=10.0d0*(d-nnj)
            llj=d
            jlj=2*llj+i
            kkj=-i*((jlj+1)/2)
            d=100.0d0*(d-llj)
            Nq(nj)=d+0.1d0
            do ni=1,ns
                if (nnj == Nn(ni) .and. Kk(ni) == kkj) goto 200
            end do
            write( 6,25) nj,nnj,llj,kkj
            write(11,25) nj,nnj,llj,kkj
    25      format(/2X,'no orbital for shell ',I3,': n,l,k=',3I4)
            stop
    200     Nip(nj)=ni
            if (Nsu < ni) then
                Nsu=ni
            end if
        end do
        nec=0
        if (Nso /=  0) then
           do ni=1,Nso
              nec=nec+Nq(ni)
            end do
        end if
        do ni=1,Nsu
           imax=2*Jj(ni)+1
           do j=1,imax,2
              Nst=Nst+1
          end do
        end do
        write( 6,35) Nsu,Ne,nec,Nst
    35  format(4X,'Number of actually used orbitals: Nsu =',I3, &
               /4X,'Ne  =',I3,7X,'nec =',I3,7X,'Nst =',I7)
        n=0
        ic=0
        i0=0
        i=0
        nmin=Nso+1
        do ni=nmin,Nsp
           i=i+1
           n=n+Nq(ni)
           if (n >= Ne) then
                ic=ic+1
                if (n > Ne) then
                    write( 6,45) ic
                    write(11,45) ic
    45              format(/2X,'wrong number of electrons'/ &
                            2X,'for configuration ICONF =',I4/)
                    stop
                end if
                Nvc(ic)=i
                Nc0(ic)=Nso+i0
                i0=i0+i
                n=0
                i=0
           end if
        end do
    !     - - - - - - - - - - - - - - - - - - - - - - - - -
        write( 6,55)
        write(11,55)
    55  format(1X,71('='))
        do ni=1,Nso
           l =Ll(ni)+1
           lll(ni)=let(l)
        end do
        if (Kout > 1) then
            write(11,65) (Nn(i),lll(i),Jj(i),Nq(i),i=1,Nso)
    65      format (1X,'Core:',6(I2,A1,'(',I1,'/2)',I2,';'),&
                           /6X,6(I2,A1,'(',I1,'/2)',I2,';'),&
                           /6X,6(I2,A1,'(',I1,'/2)',I2,';'),&
                           /6X,6(I2,A1,'(',I1,'/2)',I2,';'),&
                           /6X,6(I2,A1,'(',I1,'/2)',I2,';'),&
                           /6X,6(I2,A1,'(',I1,'/2)',I2,';'))
        end if
    !   write( 6,55)
        if (Kout > 1) then
            write(11,55)
        end if
        do ic=1,Nc
            n1=Nc0(ic)+1
            n2=Nc0(ic)+Nvc(ic)
            do i=n1,n2
                i1=i-n1+1
                ni=Nip(i)
                l=Ll(ni)+1
                lll(i1)=let(l)
                jjj(i1)=Jj(ni)
                nnn(i1)=Nn(ni)
                nqq(i1)=Nq(i)
                if (Nq(i) > jjj(i1)+1) then
                    write(11,75) ni,nnn(i1),lll(i1),jjj(i1),nqq(i1)
     75             format(/2X,'wrong number of electrons'/ &
                            2X,'for the shell:',I3,3X,I2,A1,I2,'/2', &
                            ' (',F6.3,')')
                    Stop
                end if
            end do
            n=n2-n1+1
            if (Kout > 1) then
                write(11,85) ic,(nnn(i),lll(i),jjj(i),nqq(i),i=1,n)
     85         format(I5,'#',6(I2,A1,'(',I1,'/2)',I2,';'), &
                    /6X,6(I2,A1,'(',I1,'/2)',I2,';'), &
                    /6X,6(I2,A1,'(',I1,'/2)',I2,';'))
            end if
        end do
        if (Kout > 1) then
            write(11,55)
        end if
        do ni=Nso+1,Nsu
            read(12,rec=2*ni+7) p
            Eps(ni)=-p(ii+1)
        end do
        write(11,95) (i,Eps(i),i=Nso+1,Nsu)
     95 format(' HF energies are read from DAT file', /5(I5,F10.6))
        close(unit=12)
        open(unit=16,file='CONF.GNT',status='OLD',form='UNFORMATTED')
        read(16) (In(i),i=1,IPgnt)
        read(16) (Gnt(i),i=1,IPgnt)
        close(unit=16)
        Return
    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    700 write( 6,105)
        write(11,105)
    105 format(/2X,'file CONF.DAT is absent'/)
        stop
    end subroutine Init

    subroutine BcastParams(mype,npes)
        use mpi
        implicit none
        integer :: mype, npes, mpierr
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Am, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(XJ_av, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nso, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(nrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nlv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(K_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(K_sms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ksig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kdsig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kbrt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kout, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kecp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(C_is, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Gj, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Cut0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ncpt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Lmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nsu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(num_is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ngint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nhint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        return
    end subroutine BcastParams

    subroutine AllocatePTArrays(mype,npes)
        use mpi
        implicit none
        integer :: i, mpierr, mype, npes
        
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        if (.not. allocated(Nvc)) allocate(Nvc(Nc))
        if (.not. allocated(Nc0)) allocate(Nc0(Nc))
        if (.not. allocated(Ndc)) allocate(Ndc(Nc))
        if (.not. allocated(Ndcnr)) allocate(Ndcnr(Nc))
        if (.not. allocated(Jz)) allocate(Jz(Nst))
        if (.not. allocated(Nh)) allocate(Nh(Nst))
        if (.not. allocated(Diag)) allocate(Diag(Nd))
        if (.not. allocated(Rint1)) allocate(Rint1(Nhint))
        if (.not. allocated(Iint1)) allocate(Iint1(Nhint))
        if (.not. allocated(Iint2)) allocate(Iint2(Ngint))
        if (.not. allocated(Iint3)) allocate(Iint3(Ngint))
        if (.not. allocated(Rint2)) allocate(Rint2(IPbr,Ngint))
        if (.not. allocated(IntOrd)) allocate(IntOrd(nrd))
        if (.not. allocated(Iarr)) allocate(Iarr(Ne,Nd))
        if (.not. allocated(DVnr)) allocate(DVnr(Nc))

        if (Ksig /= 0) then
            if (.not. allocated(R_is)) allocate(R_is(Nhint))
            if (.not. allocated(I_is)) allocate(I_is(Nhint))
        end if

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        return
    end subroutine AllocatePTArrays

    subroutine BcastPTArrays(mype,npes)
        use mpi
    	implicit none
        integer :: i, mype, npes, mpierr

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nn(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Kk(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ll(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jj(1:Ns), Ns, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nh(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Jz(1:Nst), Nst, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndc(1:Nc), Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Rint1(1:Nhint), Nhint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint1(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint2(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Iint3(1:Ngint), Ngint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(IntOrd(1:nrd), nrd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Rint2(1:IPbr,1:Ngint), IPbr*Ngint, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
        do i=1,Ne
            call MPI_Bcast(Iarr(i,1:Nd), Nd, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        end do
        call MPI_Bcast(In, IPgnt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Gnt, IPgnt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(DVnr(1:Nc), Nc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Diag(1:Nd), Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        if (Ksig /= 0) then
            call MPI_Bcast(R_is(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            call MPI_Bcast(I_is(1:Nhint), Nhint, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        return
    end subroutine BcastPTArrays

    subroutine NR_Init
        implicit none
        integer  :: inr, nrel, ndet, nshell, k, ier, knr1, knr2, ic, inq
        real(dp)     :: x
        logical  :: dif
        integer, dimension(20)   ::  inl1, inl2, inq1, inq2
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Ndcnr(Nc),Nvcnr(Nc),NRR(Nc),NRN(Nc))
        inr=1                          ! NR config
        nrel=0
        ndet=0
        nshell=0
        Ncnr0=0
        Ncnrci=0
        Nmax=Nd
        Ndcnr=0
        Nvcnr=0
        NRR=0
        NRN=0
        call Squeeze(1,knr1,inl1,inq1)
        NRR(1)=1
        do ic=1,Nc
            nrel=nrel+1                   ! number of r.configs in nr.config
            Ndcnr(inr)=Ndcnr(inr)+Ndc(ic)
            ndet=ndet+Ndc(ic)
            nshell=nshell+Nvc(ic)
            dif=.true.
            if (ic < Nc) then
                call Squeeze(ic+1,knr2,inl2,inq2)
                if (knr1 == knr2) then
                    do k=1,knr1
                        ier=abs(inl1(k)-inl2(k))+abs(inq1(k)-inq2(k))
                        dif=ier /=  0
                        if (dif) goto 100
                    end do
                end if
            end if
    100     if (dif) then
                NRN(inr)=nrel
                Nvcnr(inr)=nshell
                if (Kout >= 1) write (11,5) inr,nrel,NDCnr(inr),Nvcnr(inr)
    5           format(4X,'NR con-n',i5,':',i3, &
                       ' R con-ns,',i5,' det-ns &',i4,' shells')
                nrel=0
                nshell=0
                if (ic == Ncci) then
                    Ncnrci=inr
                end if
                if (ic > Ncci  .and.  Ncnrci == 0) then
                    write (*,*) ' NR_init error: Ncci=',Ncci
                    write (*,*) ' NR config ends at ic=',ic
                    stop             
                end if 
                if (ic /=  Nc) then
                    inr=inr+1
                    NRR(inr)=ic+1
                end if
            end if
            knr1=knr2
            do k=1,knr1
                inl1(k)=inl2(k)
                inq1(k)=inq2(k)
            end do
        end do
        Ncnr=inr
        write ( *,*) ' PT space:',Ncnr,' non-rel. configurations'
        write (11,*) ' PT space:',Ncnr,' non-rel. configurations'
        write ( *,*) ' CI space:',Ncnrci,' non-rel. configurations'
        write (11,*) ' CI space:',Ncnrci,' non-rel. configurations'
        return
    end subroutine NR_Init

    ! =============================================================================
    subroutine Weight_CI
        implicit none
        logical :: tail
        integer :: i, ic, kcnr, ndk, k, j, kc, icnr, ncnr0, nc0, icci, nc0j, &
                   ncnr0j, nd0j, ndci
        real(dp)  :: wj, wmax, d, dvnrn, dvnrx, dummy, wj0
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=17,file='CONF.XIJ', &
                 status='OLD',form='UNFORMATTED',err=900)
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(DVnr(Nc),B1h(Nd))
        Nd0=0    
        nc0=0    
        ncnr0=0  
        wmax=0.d0
        do j=1,Nlv
            read (17) d,dummy,ndci,(B1h(i),i=1,ndci)
            wj=1.d0                 ! norm of WF
            i=0
            ic=0
            tail=.false.
            do icnr=1,Ncnrci    !### loop over NR conf-ns
                kcnr=NRN(icnr)
                d=0.d0
                ndk=0
                do kc=1,kcnr          !### loop over R conf-ns
                    ic=ic+1
                    ndk=ndk+Ndc(ic)     !### ndk= number of det-s in NR con-n
                end do
                do k=1,ndk
                    i=i+1
                    d=d+B1h(i)**2
                end do
                if (j == 1) then
                    DVnr(icnr)=d
                else
                    DVnr(icnr)=DVnr(icnr)+d
                end if
                wj=wj-d                ! norm of the remaining tail of WF
                if ((.not.tail)  .and.  i < IPPT) then
                    nd0j=i
                    nc0j=ic  
                    ncnr0j=icnr
                    wj0=wj
                end if
                tail=wj < Cut0
            end do
            write(*,5) j,nd0j,nc0j,ncnr0j,wj0
    5       format(4x,'WF ',i2,' Nd0=',i6,' Nc0=',i5,' Ncnr0=',i5, &
                 ' wj=',f9.6)
            if (wj0 > wmax) wmax=wj0
     
            if (Nd0 < nd0j) then
                Nd0=nd0j
                nc0=nc0j  
                ncnr0=ncnr0j
            end if            
        end do
        close(unit=17)
        Ncp0=nc0
        write( *,15) Nd0,Ncp0,ncnr0,wmax,Cut0
        write(11,15) Nd0,Ncp0,ncnr0,wmax,Cut0
    15  format(/4x,'PT block: Nd0=',i6,' Nc0=',i5,' Ncnr0=',i5, &
               /4x,'Actual Cutoff=',f9.6,' Cut0=',f9.6)
        
        if (Kout >= 1) write(11,*) ' Weights of all NR configurations:'
        icci=0
        dvnrn=11.d1
        do ic=1,Ncnr
            if (ic <= Ncnrci  .and.  DVnr(ic) < dvnrx) icci=icci+1
            if (ic <= Ncnrci  .and.  DVnr(ic) < dvnrn) dvnrn=DVnr(ic)
            if (Kout >= 1) then
                write(11,25) ic,DVnr(ic)
    25          format(4x,'NR conf ',i4,' weight ',f10.7)
                if (ic == Ncnrci) &
                  write(11,*) '   ---- End of CI space ---------'
            end if
        end do
        return
    900 write (*,*) ' CONF.XIJ file with WFs is absent'
        stop
    end subroutine Weight_CI
    ! =============================================================================
    subroutine PT_Init(npes, mype)
    	use mpi
        implicit none
        integer  :: i, ni, j, n, k, kx, ic, kd, ndci
        integer :: npes, mype, mpierr
        real(dp)   :: sd, x
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (mype==0) then ! master pe only
            write ( *,5)
            write (11,5)
    5       format(/4X,'PT_init')
            ! - - - - - - - - - - - - - - - - - - - - - - - - -
            open(unit=16,file='CONF.XIJ',status='OLD', &
                 form='UNFORMATTED',err=900)
            read(16) x,x,Nd1
            kd=0
            Nd0=0
            do i=1,Nc
              kd=kd+Ndc(i)
              if (i == Ncci) ndci=kd
              if (i == Ncp0) Nd0=kd
            end do
            if (ndci /=  Nd1) then
              write (*,*) ' Wrong length of vectors in CONF.XIJ:'
              write (*,*) ' Nd1 =',Nd1,'; expected ',ndci
              stop
            end if
        
            write ( *,25) Ncci,Nd1,Ncp0,Nd0,Nc,Nd
            write (11,25) Ncci,Nd1,Ncp0,Nd0,Nc,Nd
    25      format(4X,'CI block: Ncci=',I6,' Nd1=',I8,/4X, &
                  'PT block: (Ncp0=',I6,' Nd0=',I8,')*(Nc=',I7,' Nd=',I9,')')
            if (Nd0 > Nd1) then
              write (*,*) ' can not work with Nd0 > Nd1'
              Stop
            end if
            if (Nd0 > IPPT) then
              write (*,*) ' can not work with Nd0 > IPPT =',IPPT
              Stop
            end if
            if (Nd0 == 0) then
              write (*,*) ' Nothing to do for Nd0=0 (Ncp0=',Ncp0,')'
              Stop
            end if
        
            Nlv=min(Nlv,Nd0)
        end if ! end master pe only

        call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        allocate(X0(Nd,Nlv), En(Nlv), EnG(Nlv), Xj(Nlv))

        call DiagH(Nd1,npes,mype)      !# forms array Diag(i)=H(i,i)

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        if (mype==0) then ! master pe only
            if (average_diag) then
              n=Nd1
              do ic=Ncci+1,Nc      !# averagin diagonal for each configuration
                kx=Ndc(ic)
                sd=0.d0
                if (kx /=  0) then
                  do k=1,kx
                    sd=sd+Diag(n+k)
                  end do
                  sd=sd/kx
                  do k=1,kx
                    Diag(n+kx)=sd
                  end do
                  n=n+kx
                end if
              end do
              if (n /=  Nd) then
                write(*,*) 'PT_init error: n=',n,' must be equal to Nd=',Nd
                stop
              else
                write( *,*) ' Diagonal averaged over relat. config.'
                write(11,*) ' Diagonal averaged over relat. config.'
              end if
            end if
            rewind(16)           !## starting from Nd1+1 to Nd
            do n=1,Nlv
              read(16,err=220,end=220) EnG(n),Xj(n),Nd1,(X0(j,n),j=1,Nd1)
              En(n)=EnG(n)+4.d0*Gj*Xj(n)*(Xj(n)+1.d0)
              ni=n
            end do
    220     close(16)
            write (*,35) ni
    35      format (4X,i5,' vectors read from file CONF.XIJ.')
            Nlv=ni
        end if
        return
    900 write(*,*)' File CONF.XIJ is absent.'
        stop
    end subroutine PT_Init
    ! =============================================================================
    subroutine PTE(npes,mype)
        use mpi
        implicit none
        integer  :: n, ic, l, i, k, ncoef, ifrac, ix, j
        integer, dimension(Nc)  :: Ndic
        integer  :: npes, mype, mpierr, interval, remainder, start, end, size, xstart
        integer  :: xend, xsize, Ndcount, istart, iend, disp, lastsize
        real  :: start_time, stop_time
        real(dp)  :: ei, x, e0, del, ci, x2, des, des0, dEx, dVs
        real(dp), dimension(IPPT)  :: E1
        real(dp), dimension(Nlv,Nc) :: dE, dV, dE1, dV1
        integer, dimension(2,Nc) :: sizeNcNd
        integer     :: totalNcNd, pecounter, counter
        integer, dimension(npes) :: sizes, sizes1, idisp
        logical     :: first
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(Ey(Nlv),Ndirc(Nc),DEnr(Nc))

        if (mype == 0) then
        write ( *,5)
        write (11,5)
    5   format(/4X,'Second order corrections to the energies', &
                  ' & first order eigenvectors')
        
        E1=0.d0
        do n=1,Nlv
            Ey(n)=0.d0
        end do
        
        do ic=Ncnrci+1,Ncnr
            DEnr(ic)=0.d0
            DVnr(ic)=0.d0
        end do

        Ndirc(1)= 0
        do ic= 2,Nc
            Ndirc(ic) = Ndirc(ic-1) + Ndc(ic-1)
        end do
    
        Ndic=0
        do ic= 1,Ncnr
            if (ic == 1) Then
                Ndic(1)= 0
            else
                Ndic(ic)= Ndic(ic-1)+Ndcnr(ic-1)
            end If
            do n=1,Nlv
                dE(n,ic)=0.e0
                dV(n,ic)=0.e0
            end do
        end do
        end if

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Nd1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ncnr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ncnrci, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndcnr, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ey, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(En, Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndirc, Nc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Ndic, Ncnr, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(X0(1:Nd1,1:Nlv), Nd1*Nlv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(Diag, Nd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        totalNcNd = 0

        do ic=Ncnrci+1,Ncnr
            do l=1,Ndcnr(ic)
                totalNcNd = totalNcNd + 1
                sizeNcNd(1, ic) = ic
                sizeNcNd(2, ic) = Ndcnr(ic) ! size of each config
            end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        do ic=Ncnrci+1,Ncnr ! loop over configs in PT space

            interval = sizeNcNd(2,ic)/npes
            remainder = mod(sizeNcNd(2,ic),npes)
            !print*, 'starting ic=', ic, 'of ', Ncnr
            do i=0,npes-1
                if (mype == npes-1) then
                    start = 1+mype*interval
                    end = (mype+1)*interval+remainder
                    istart = Ndic(ic)+start
                    iend = Ndic(ic)+end
                    lastsize= end - start + 1
                    !print*,ic,mype,start,end,sizeNcNd(2,ic)
                else if (mype /= npes-1) then
                    start = 1+mype*interval
                    end = (mype+1)*interval
                    istart = Ndic(ic)+start
                    iend = Ndic(ic)+end
                    !print*,ic,mype,start,end,sizeNcNd(2,ic)
                    !print*, 'rank ', mype, 'has istart=', istart, 'iend=', iend, 'and size = ', size
                end if
                size = end - start + 1
            end do  

            do i=1,npes-1
                sizes(i) = interval
            end do
            
            sizes(npes) = interval + remainder

            idisp(1)=0
            do i=2, npes
                idisp(i) = (i-1)*sizes(1)
            end do
            do l=start, end  ! loop over dets in configs in PT space
                i= Ndic(ic)+l
                ei= -Diag(i)
                call FormH(i,Nd0,E1)    !# evalulation of i-th MEs of ALL
                do n=1,Nlv              !## first order correction vectors
                    x=0.d0                !### and writing them to VEC.TMP,
                    e0=En(n)
                    do k=1,Nd0
                        x=x+E1(k)*X0(k,n)
                    end do
                    x2=x*x
                    del=e0-ei
                    ci=x/del
                    dE(n,ic)= dE(n,ic)+x2/del    ! contribution of one NR config.
                    dV(n,ic)= dV(n,ic)+(ci)**2
                    X0(i,n)=ci                   ! first order eigenvector
                end do
            end do
            call MPI_Barrier(MPI_COMM_WORLD, mpierr)
            do n=1,Nlv
                call MPI_Allreduce(dE(n,ic),dE1(n,ic), 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, MPI_COMM_WORLD, mpierr)
                call MPI_Allreduce(dV(n,ic),dV1(n,ic), 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, MPI_COMM_WORLD, mpierr)
                call MPI_Gatherv(X0(istart:iend,n), size, MPI_DOUBLE_PRECISION, &
                               X0(Ndic(ic)+1:Ndic(ic)+Ndcnr(ic),n), &
                               sizes, idisp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            end do

            dEx=0.d0
            dVs=0.d0
            dE=dE1
            dV=dV1
            do n=1,Nlv
                Ey(n)= Ey(n)+ dE(n,ic)    ! sum of all contributions to n's energy
                if (dE(n,ic) > dEx) then
                    dEx= dE(n,ic)
                end if
                dVs= dVs+dV(n,ic)
            end do
            DEnr(ic)= dEx
            DVnr(ic)= dVs
            if (mype==0) then
                print*, 'finished ', ic, 'of ', Ncnr
            end if

        end do
    
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        if (mype == 0) then
        write ( *,25) NcCI,Nd1,Ncp0,Nd0,Nc,Nd
        write (11,25) NcCI,Nd1,Ncp0,Nd0,Nc,Nd
    25  format(/82('='), &
             /' CI space [',i6,'/',i8,']; PT space [',i5,'/',i6, &
             '] * [',i7,'/',i9,']' &
             /' N',6X,'J',12X,'Ev0',6X,'DEL(CM**-1)',7X,'dE(PT2)', &
             6X,'Ev(PT2)',3X,'DEL(CM**-1)',/82('-'))
    
        do n=1,Nlv
            des0=(En(1)-En(n))*2*DPRy
            des =des0 + (Ey(1)-Ey(n))*2*DPRy
            write ( *,55) n,Xj(n),En(n),des0,Ey(n),En(n)+Ey(n),des
            write (11,55) n,Xj(n),En(n),des0,Ey(n),En(n)+Ey(n),des
    55      format(I2,F12.8,F14.8,F13.2,2F14.8,F13.2)
        end do
        write ( *,65)
        write (11,65)
    65  format(82('='))
    
        if (Kout >= 1) write (11,75)
    75  format(/4x,'PT contributions of non-relat. configurations', &
               /4x,'ic     max dE   Full weight',/4x,27('-'))
    
        dvnrx=0.d0
        ix=0
        do ic=Ncnrci+1,Ncnr
            if (DVnr(ic) > dvnrx) then
                dvnrx=DVnr(ic)
                ix=ic
            end if
            if (Kout >= 1) write(11,85) ic,DEnr(ic),DVnr(ic)
    85      format(i6,2f12.7)
        end do

        if (Kout >= 1) write (11,95)
    95  format(4x,27('-'))
        write ( *,105) dvnrx,ix
        write (11,105) dvnrx,ix
    105 format('Max PT weight ',e11.3,' of config. ',i6)
        end if

        return
    end subroutine PTE
    ! =============================================================================
    subroutine SaveVectors
        implicit none
        integer  :: n, i, j
        real(dp)   :: s
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        write ( *,5)
        write (11,5)
    5   format(/4X,'Normalizing & saving vectors')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=16,file='CONF_PT.XIJ',status='UNKNOWN', &
             form='UNFORMATTED')
        close(unit=16,status='DELETE')
        open(unit=16,file='CONF_PT.XIJ',status='NEW', &
             form='UNFORMATTED')
        do n=1,Nlv
            s=0.d0
            do i=1,Nd
                s=s+X0(i,n)**2
            end do
            s=1.d0/sqrt(s)
            do i=1,Nd
                !print*,i,s,X0(i,n)
                X0(i,n)=s*X0(i,n)
            end do
            EnG(n)=En(n)-4.d0*Gj*Xj(n)*(Xj(n)+1.d0)
            write (16) EnG(n),Xj(n),Nd,(X0(j,n),j=1,Nd)
            write( *,15) n,EnG(n),Xj(n),s
            write(11,15) n,EnG(n),Xj(n),s
    15      format(4x,'vector ',i2,' E_G=',f13.8,' J=',f8.5,' S_norm=', &
                 f8.5,' saved')
        end do
        close(16)
        write (*,35) 
    35  format (4X,' vectors saved to file CONF_PT.XIJ.')
        return
    end subroutine SaveVectors
    ! =============================================================================
    subroutine Weight_PT
        implicit none
        integer  :: icci, ic, icpt
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        icci=0
        dvnrn=11.d1
        do ic=1,Ncnrci
          if (DVnr(ic) < dvnrx) icci=icci+1
          if (DVnr(ic) < dvnrn) dvnrn=DVnr(ic)
        end do
    
        icpt=0
        do ic=Ncnrci+1,Ncnr
          if (DVnr(ic) > dvnrn) icpt=icpt+1
        end do
        write ( *,15) icci,Ncnrci,dvnrx,icpt,Ncnr,dvnrn
        write (11,15) icci,Ncnrci,dvnrx,icpt,Ncnr,dvnrn
    15  format('Weights of ',i6,' (from ',i6, &
             ') NR config-s in CI below max PT w-t',e11.3, &
             /'Weights of ',i6,' (from ',i6, &
             ') config-s in PT space above min CI w-t',e11.3)
        return
    end subroutine Weight_PT
    ! =============================================================================
    subroutine Cutoff(Nc_0)
        implicit none
        integer  :: izer, icc, icnr, k, ii, ic, iinr, kvar, ktf, last, i, &
                    Nc_0
        integer, dimension(15)  :: Nconf, Mconf
        real(dp)  :: wx, wl, xx, dlx, dln, ctf, xcc
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(W(Nc,2))
        izer=0
        do icc=1,15
           Nconf(icc)=0
           Mconf(icc)=0
        end do
        do icnr=1,Ncnr                     !### transforms weights to
          wx=DVnr(icnr)                    !###  log scale
          if (wx > 0.d0) then
            wl=dlog(wx)*0.4342945d0        !### (factor = lg e)
          else
            izer=izer+1
            wl=-99.d0
          end if
          W(icnr,1)=wl
          do icc=1,15
            xcc=2-icc
            if (wl < xcc) Nconf(icc)=Nconf(icc)+1
          end do
          xx=dvnrx
          do k=1,10
            if (wx > xx) Mconf(k)=Mconf(k)+1
            xx=0.5d0*xx
          end do
        end do
        do icc=1,15
          if (Nconf(icc) /= 0) last=icc
        end do
        if (Kout >= 1) write(11,5) izer, &
            (Ncnr-Nconf(i),Nconf(i),2-i,i=1,last)
    5   format(4X,I6,' nonrel. conf-s with zero weight', &
              /(4X,I6,' weight above and ',I6,' below 10**(',I3,')'))
    
        write(*,15) (Mconf(k),dvnrx,k-1,dvnrx/2**(k-1),k=1,10)
        write(11,15) (Mconf(k),dvnrx,k-1,dvnrx/2**(k-1),k=1,10)
    15  format(/(4x,i6,' weights above ',d12.3,'/2**',i2,' = ',e11.3))

        open(unit=21, file='cpt.in')
        read(21,*) ktf,kvar
        close (21)

        write(*,*) ' Give cutoff value of k (0-9):'
        dlx=dlog(dvnrx)*0.4342945d0
        dln=dlog(dvnrn)*0.4342945d0
        ctf=dlx - 0.301029996d0*ktf
    
        write(*,*) 'Choose reordering variant:'
        write(*,*) '(1) - all; (2) - keep PT block; (3) - keep CI space'
        Nc_0=0
        if (kvar == 2) Nc_0=Ncnr0
        if (kvar == 3) Nc_0=Ncnrci
    
        ii=0
        ic=0
        iinr=0
        do icnr=1,Ncnr
          ic=ic+NRN(icnr)
          if (W(icnr,1) > ctf .or. icnr <= Nc_0) then
            iinr=iinr+1
            ii=ii+NRN(icnr)
            W(icnr,2)=1.d0
          else
            W(icnr,2)=0.d0
          end if
        end do
        if (ic /= Nc) then
          write(*,*) ' Cutoff: ic=',ic,' is not equal Nc=',Nc
        end if
        Ncnew=ii
        write ( *,105) 10.d0**ctf,iinr,Ncnrci,Ncnew,Ncci
    105 format(4X,'For cutoff ',E10.3,' CI: Ncnr=',I6,' (was ',I6, &
             ') and Nc =',I6,' was(',I6,')')
        return
    end subroutine Cutoff
    ! =============================================================================
    subroutine Sort(Nc_0)
        implicit none
        integer  :: ic, i1, irpl, mm, k, kc, Nc_0
        real(dp)  :: x, y
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        irpl=0                    !### counts replacements
        do ic=Nc_0+1,Ncnr-1
          x=-100.d0
          do kc=ic,Ncnr           !### search for max from the rest
            y=W(kc,1)
            if (y > x) then
              i1=kc
              x=y
            end if
          end do
          if (ic /= i1) then      !### placing max in front
            irpl=irpl+1
            do k=1,2
              x=W(i1,k)
              W(i1,k)=W(ic,k)
              W(ic,k)=x
            end do
            mm=NRN(i1)
            NRN(i1)=NRN(ic)
            NRN(ic)=mm
            mm=NRR(i1)
            NRR(i1)=NRR(ic)
            NRR(ic)=mm
            mm=Nvcnr(i1)
            Nvcnr(i1)=Nvcnr(ic)
            Nvcnr(ic)=mm
          end if
        end do
        write(*,*) irpl,' replacements first'
        return
    end subroutine Sort
    ! =============================================================================
    subroutine NewConfList(Nc_0)
        implicit none
        character(len=1), dimension(16)  :: name
        character(len=1), dimension(5)   :: ch1
        character(len=1)  :: str*70
        integer  :: Nc_0, icnr, ic0, ic1, kd, kd4, i, num, numt, nci, nrci, &
                    icnr1, num1, n, n1, n2, kc, kcnr, kc4, ic
        logical  :: inci
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,*) ' Forming new list of configurations...'
        close(11)
        open(unit=11,file='CONF_new.INP',status='UNKNOWN')
        close(11,status='DELETE')
        open(unit=11,file='CONF_new.INP',status='NEW')
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Nc_0 < Nc4) then           ! we need to find new Nc4 
          kd4=0
          kc4=0
          do icnr=1,Ncnr
            ic0=NRR(icnr)
            ic1=ic0+NRN(icnr)-1
            kd=0
            do ic=ic0,ic1
              kd=kd+Ndc(ic)
            end do
            if (kd4+kd > IP1) goto 100
            kd4=kd4+kd
            kc4=kc4+NRN(icnr)
          end do
    100   if (kc4 > 0) then 
            Nc4=kc4
            else
            Nc4=NRN(1)
          end if
          write(*,*)' New Nc4=',Nc4,' (Nd4=',kd4,')'         
        end if
    
    !   input from the file 'CONF.INP'
        open(unit=10,file='CONF.INP',status='OLD')
        read (10,5) name
    5   format(1X,16A1)
        write(11,'(1X,16A1,A12)') name,' (reordered)'
        read (10,15) (ch1(i),i=1,5),Z
        write(11,15) (ch1(i),i=1,5),Z
        read (10,15) (ch1(i),i=1,5),Am
        write(11,15) (ch1(i),i=1,5),Am
        read (10,15) (ch1(i),i=1,5),XJ_av
        write(11,15) (ch1(i),i=1,5),XJ_av
        read (10,15) (ch1(i),i=1,5),Jm
        write(11,15) (ch1(i),i=1,5),Jm
    15  format (5A1,F5.1)
        read (10,25) (ch1(i),i=1,5),Nso
        write(11,25) (ch1(i),i=1,5),Nso
        read (10,25) (ch1(i),i=1,5),Nc
        write(11,25) (ch1(i),i=1,5),Ncnew
        read (10,25) (ch1(i),i=1,5),Kv
        write(11,25) (ch1(i),i=1,5),Kv
        read (10,25) (ch1(i),i=1,5),Nlv
        write(11,25) (ch1(i),i=1,5),Nlv
        read (10,25) (ch1(i),i=1,5),Ne
        write(11,25) (ch1(i),i=1,5),Ne
    25  format (5A1,I6)
    
    200 read (10,35) (ch1(i),i=1,5),str
    35  format(5A1,A)
        if (ch1(2) /= 'K'.AND.ch1(2) /= 'k') goto 300
        if (ch1(3) /= 'L'.AND.ch1(3) /= 'l') goto 300
        if (ch1(4) /= '4') goto 300
        if (Nc_0 == Ncnrci) then
          kl4=2            ! initial approximation from disk
        else
          kl4=1
        end if
        write(11,25) (ch1(i),i=1,5),kl4
        goto 200    
    
    300 if (ch1(2) /= 'N'.AND.ch1(2) /= 'n') goto 400
        if (ch1(3) /= 'C'.AND.ch1(3) /= 'c') goto 400
        if (ch1(4) /= '4') goto 400
        write(11,25) (ch1(i),i=1,5),Nc4
        goto 200
        
    400 if(ch1(1) == ' '.AND.ch1(2) == ' ') goto 500
        write(11,35) (ch1(i),i=1,5),str
        goto 200
    500 close(10)
        write(11,45) (Qnl(i),i=1,Nso)
    45  format (6(4X,F7.4))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        num=0
        numt=0
        nci=0
        nrci=0
        inci=.TRUE.
        do icnr=1,Ncnr
          kcnr=NRN(icnr)
          numt=numt+kcnr
          if (W(icnr,2) /= 0.d0) then
            nrci=nrci+1
          else
            if (inci) then
              write(11,55)
    55        format('<CI<')
              inci=.FALSE.
            end if
          end if
          ic=NRR(icnr)-1
          icnr1=mod(icnr,10000)
          write(11,65) icnr1,10**W(icnr,1),Nvcnr(icnr)
    65    format(I4,62X,E9.2,i6)
          do kc=1,kcnr
            ic=ic+1        !### old index of current configuration
            num=num+1      !### new index of current configuration
            if (inci) nci=nci+1
            num1=mod(num,10000)
            n1=Nc0(ic)+1
            n2=Nc0(ic)+Nvc(ic)
            write(11,75) num1,(Qnl(n),n=n1,n2)
    75      format(I4,6(F7.4,4X),/6(4X,F7.4))
          end do
        end do
        if (nci == Ncnew) then
          write(11,85)
    85    format(1X,'>>>>>>>>>> END <<<<<<<<<<')
          close(11)
          write(*,*)' CONF_new.INP is formed.'
          write(*,*)' New CI space:',nci,' PT space:',num,' config.'
        else
          write(*,95) nci,Ncnew
    95    format(' Error: nci=',I6,' not equal to Ncnew=',I6)
          read(*,*)
          return
        end if
        return
    end subroutine NewConfList

    subroutine Squeeze(ic,knr,inl,inq)
        implicit none
        integer  :: i1, i2, knr, i, ic, nl, mq, nl0, mq0
        logical :: dif
        real(dp)  :: x
        integer, dimension(20)  :: inl, inq
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=Nc0(ic)
        i2=i1+Nvc(ic)
        knr=0            
        nl0=0
    
        i=i1+1
        x=100*abs(Qnl(i))
        nl=x+1.d-5
        mq=100*(x-nl)+1.d-5
    
    100 if (knr >= 20) then
            write (*,*) ' Squeeze error for ic=',ic
            write (*,*) ' array overflow. knr=',knr
            read (*,*)
        end if
        if (nl == nl0) then
            knr=knr+1
            inl(knr)=nl0
            inq(knr)=mq0+mq
            nl=0
            if (i == i2) goto 200
        else
            if (nl0 /= 0) then
                knr=knr+1
                inl(knr)=nl0
                inq(knr)=mq0
            end if
        end if
        if (i == i2) then
            knr=knr+1
            inl(knr)=nl
            inq(knr)=mq
            goto 200
        end if
    
        nl0=nl
        mq0=mq
        i=i+1
        x=100*abs(Qnl(i))
        nl=x+1.d-5
        mq=100*(x-nl)+1.d-5
        goto 100
    
    200 continue
        return
    end subroutine Squeeze

    subroutine FormH(n,nd0,E1)
        use determinants, only : CompC, Gdet
        implicit none
        integer  :: ic, kx, k, n, nd0, icomp
        integer, allocatable, dimension(:) :: idet1, idet2
        real(dp), dimension(IPPT)  :: E1
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (.not. allocated(idet1)) allocate(idet1(Ne))
        if (.not. allocated(idet2)) allocate(idet2(Ne))
        if (.not. allocated(ic1)) allocate(ic1(Ne))
        if (.not. allocated(ic2)) allocate(ic2(Ne))

        E1=0.d0
        ! calculation of the matrix elements
        call Gdet(n,idet1)
        do ic=1,Nc
            kx= Ndc(ic)
            k = Ndirc(ic)
            if (k+kx > nd0) kx=nd0-k
            if (kx /= 0) then
                call Gdet(k+1,idet2)
                call CompC(idet1,idet2,icomp)
                if (icomp <= 2) then
                    do k= Ndirc(ic)+1,Ndirc(ic)+kx
                        call Gdet(k,idet2)
                        E1(k)= Hmltn(idet1,idet2)
                    end do
                end if
            end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        return
    end subroutine FormH

    real(dp) function Hmltn(idet1,idet2)
        use integrals, only : Gint, Hint
        use determinants, only : Rspq
        implicit none
        integer  :: nf, iq, i1, i2, j1, j2, is, jq, jq0
        integer, allocatable, dimension(:)  :: idet1, idet2
        real(dp)  :: t
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        call Rspq(idet1,idet2,is,nf,i1,i2,j1,j2)
        t=0.d0
        if (nf <= 2) then
            if (nf == 2) then
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants differ by two functions
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
                t=t+Gint(i2,j2,i1,j1)*is !### det_k goes first!
            else if (nf == 1) then
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants differ by one function
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
                do iq=1,Ne
                    i1=idet1(iq)
                    if (i1 /= j1) then
                        t=t+Gint(j2,i1,j1,i1)*is
                    end if
                end do
                t=t+Hint(j2,j1)*is
            else
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! determinants are equal
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
                do iq=1,Ne
                    i1=idet1(iq)
                    jq0=iq+1
                    if (jq0 <= Ne) then
                        do jq=jq0,Ne
                            j1=idet1(jq)
                            t=t+Gint(i1,j1,i1,j1)*is
                        end do
                    end if
                    t=t+Hint(i1,i1)*is
                end do
            end if
        end if
        Hmltn=t
        return
    end function Hmltn

    subroutine DiagH(nd0,npes,mype)
    	use mpi
        use determinants, only : Gdet
        implicit none
        integer   :: n, i, nd0, start, end, pe, j
        integer   :: npes, mype, interval, remainder, mpierr, size, count
        integer, allocatable, dimension(:)  :: idet1, idet2
        integer, dimension(npes) :: disp, sizes, ends
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (.not.allocated(idet1)) allocate(idet1(Ne))
        if (.not.allocated(idet2)) allocate(idet2(Ne))
        do n=1,Nd
            Diag(n)=0.d0
        end do
        call MPI_Bcast(nd0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        interval = (Nd-(nd0+1))/npes
        remainder = mod(Nd-(nd0+1),npes)
        count = 0
        ! setting up sizes and displacements for MPI
        do i=0,npes-1
            if (i == npes-1) then
                start = nd0+1+i*interval
                end = nd0+1+(i+1)*interval+remainder
            else
                start = nd0+1+i*interval
                end = nd0+(i+1)*interval
            end if
            size = end - start + 1
            sizes(i+1) = size
            ends(i+1) = end
        end do
        do i=0,npes-1
            if (mype == npes-1) then
                start = nd0+1+mype*interval
                end = nd0+1+(mype+1)*interval+remainder
            else if (mype /= npes-1) then
                start = nd0+1+mype*interval
                end = nd0+(mype+1)*interval
            end if
            size = end - start + 1
            sizes(mype+1) = size
            ends(mype+1) = end
        end do        
        disp(1)=0
        do i=2, npes
            disp(i) = ends(i-1) - nd0
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        !   calculation of the matrix elements
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        if (mype == 0) write (*,*) ' Formation of diagonal'

        do n=start, end
            call Gdet(n,idet1)
            do i=1,Ne
                idet2(i)=idet1(i)
            end do
            Diag(n)=Hmltn(idet1,idet2)
            count=count+1
            if (mype==0) then
                do j=1,10
                    if (count==size/j) then
                        print*, 'diag is', (10-j)*10, '% done'
                    end if
                end do
            end if
        end do
        
        !print*,'pe=',mype,'has finished working with count=', count

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        
        call MPI_Gatherv(Diag(start:end), size, MPI_DOUBLE_PRECISION, Diag(nd0+1:Nd), &
                           sizes, disp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

        call MPI_Barrier(MPI_COMM_WORLD, mpierr)

        if (mype==0) then
        write (*,*) ' Diagonal formed'
    	end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        return
    end subroutine DiagH

end module conf_pt_aux