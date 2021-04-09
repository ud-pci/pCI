module conf_pt_variables

    use conf_variables

    implicit none
    
    public

    integer, parameter :: IPPT = 5000
    ! global variables for conf_pt
    integer     :: Nd1, Ncci, Nmax, Ncp0, Ncnr, Ncnrci, Ncnr0, Ncnew, &
                   NcOld, n_is, KmaxScr, KsymScr, NsumScr, MaxScr, ktf, kvar
    real(dp)    :: dvnrx, dvnrn
    logical     :: average_diag

    integer, allocatable, dimension(:)    :: Ndcnr, Nvcnr, NRR, NRN, Ndirc
    real(dp), allocatable, dimension(:)   :: B1h, En, Xj, EnG, Ey, DEnr, DVnr
    real(dp), allocatable, dimension(:,:) :: X0

    save

end module conf_pt_variables