Module conf_variables

    Use params

    Implicit None
    
    Public

    ! Set Kw to determine whether CONF.HIJ will be written or not
    ! Kw=0 - CONF.HIJ will not be written
    ! Kw=1 - CONF.HIJ will be written
    Integer, Parameter :: Kw=1 

    ! Set kXIJ to determine the interval in which CONF.XIJ will be written
    ! e.g. kXIJ=10 - CONF.XIJ is written every 10 Davidson iterations
    Integer, Parameter :: kXIJ=10

    Integer             :: Kexn, Ksig, Kdsig, K_prj, Ngint, Nhint, Ifail 
    Integer             :: NmaxS, LmaxS, Nmax1, Lmax1, NmaxScr, LmaxScr
    Integer             :: IPlv, nrd, kdavidson, K_sms, num_is
    Integer             :: Kherr, Kgerr, Lmax, NhintS, NgintS, Kmax, Ksym, Nsum
    Integer             :: Nd0, maxJcore, vaBinSize
    Integer(Kind=Int64) :: iscr, NumH, NumJ, counter1, ij8J, ih8H
    Integer(Kind=Int64) :: memTotalPerCPU, memDvdsn, memFormH, memStaticArrays, memEstimate
    Real(dp)            :: Ecore, Hmin, XJ_av, dR_N, E_k, xscr, K_gnt, Cnx, vmax, E_0

    Integer,  Allocatable, Dimension(:)    :: IntOrd, IntOrdS, Iconverge, iconf1, iconf2
    Integer,  Allocatable, Dimension(:)    :: Iint1, Iint2, Iint3, Iint1S, Iint2S, Iint3S, I_is
    Real,     Allocatable, Dimension(:)    :: Rsig, Dsig, Esig, Rint2S, Dint2S, Eint2S, R_is, Scr
    Real,     Allocatable, Dimension(:,:)  :: Rint2
    Real(dp), Allocatable, Dimension(:)    :: Rint1, Tk, Tj, Tl, Ts, Diag, D, E
    Real(dp), Allocatable, Dimension(:,:)  :: Z1, P, W

    Integer,  Dimension(IPgnt) :: In
    Real(dp), Dimension(IPs)   :: Eps
    Real(dp), Dimension(IPgnt) :: Gnt
    Real(dp), Dimension(IP1)   :: C

    Type Matrix
        Integer,  Allocatable, Dimension(:) :: n, k
        Real(dp), Allocatable, Dimension(:) :: t
    End Type Matrix

    Type(Matrix) :: Hamil, Jsq

    Save
    
End Module