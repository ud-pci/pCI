Module conf_variables

    Use params

    Implicit None
    
    Public

    integer             :: Kexn, Kp, Ksig, Kdsig, K_prj, Ngint, Nhint, Ifail, Kw
    integer             :: NmaxS, LmaxS, Nmax1, Lmax1, NmaxScr, LmaxScr, kchunks
    integer             :: IPlv, nrd, kdavidson, K_sms, kXIJ, kJJJ, num_is
    integer             :: Kherr, Kgerr, Idel, Lmax, NhintS, NgintS, Kmax, Ksym, Nsum
    integer             :: IPdel, Jdel, ifile, ihf, Nd0
    integer             :: startH, endH, sizeH, startJ, endJ, sizeJ, maxJcore
    integer             :: H_nstart, H_nend, H_nsize, H_nmax, bin
    integer(kind=int64) :: iscr, NumH, NumJ, Ibuf0, ih8H, counter1, ij8J
    integer(kind=int64) :: memTotalPerCPU, memDvdsn, memFormH, memStaticArrays, memEstimate
    real(dp)            :: Ecore, Hmin, XJ_av, dR_N, E_k, xscr, K_gnt, Cnx, vmax, E_0

    integer, POINTER :: Iarr_win(:,:)

    integer, allocatable, dimension(:)     :: IntOrd, IntOrdS, Iconverge, iconf1, iconf2
    integer, allocatable, dimension(:)     :: Iint1, Iint2, Iint3, Iint1S, Iint2S, Iint3S, I_is
    integer, allocatable, dimension(:)     :: H_nsizes, H_ndisps
    real, allocatable, dimension(:)        :: Rsig, Dsig, Esig, Rint2S, Dint2S, Eint2S, R_is
    real, allocatable, dimension(:,:)      :: Rint2
    real(dp), allocatable, dimension(:)    :: Rint1, Tk, Tj, Diag, D, E
    real(dp), allocatable, dimension(:,:)  :: Z1, P, W

    integer, dimension(IPgnt)  :: In
    real(dp), dimension(IPs)   :: EPs
    real(dp), dimension(IPgnt) :: Gnt
    real(dp), dimension(10)    :: Scr
    real(dp), dimension(IP1)   :: C, Er

    Type Matrix
        Integer, Allocatable, Dimension(:)  :: n, k
        Real(dp), Allocatable, Dimension(:) ::  t
    End Type Matrix

    Type(Matrix) :: Hamil, Hamil0, Jsq
    save
    
End Module