Module conf_variables

    Use params

    Implicit None
    
    Public

    ! Set kXIJ to determine the interval in which CONF.XIJ will be written
    ! e.g. kXIJ=10 - CONF.XIJ is written every 10 Davidson iterations
    Integer, Parameter :: kXIJ=10

    Integer             :: Kexn=0, Ksig=0, Kdsig=0, K_prj=0, K_sms=0, Kw=0
    Integer             :: NmaxS=0, LmaxS=0, Ngint=0, Nhint=0, NhintS=0, NgintS=0
    Integer             :: IPlv, nrd, kdavidson, num_is, Nc1
    Integer             :: Kherr=0, Kgerr=0, Lmax, Kmax, Ksym, Nsum
    Integer             :: Nd0=0, vaBinSize
    Integer(Kind=Int64) :: iscr=0_int64, NumH=0_int64, NumJ=0_int64, counter1=0_int64, ij8=0_int64, ih8=0_int64
    Integer(Kind=Int64) :: memTotalPerCPU=0_int64, memDvdsn=0_int64, memFormH=0_int64, memStaticArrays=0_int64, memEstimate=0_int64
    Real(dp)            :: Ecore=0_dp, Hmin=0_dp, XJ_av, dR_N, E_k, xscr=0_dp, K_gnt, E_0

    Integer,  Allocatable, Dimension(:)    :: IntOrd, IntOrdS, Iconverge, iconf1, iconf2, In
    Integer,  Allocatable, Dimension(:)    :: Iint1, Iint2, Iint3, Iint1S, Iint2S, Iint3S, I_is
    Real,     Allocatable, Dimension(:)    :: Rsig, Dsig, Esig, Rint2S, Dint2S, Eint2S, R_is, Scr
    Real,     Allocatable, Dimension(:,:)  :: Rint2
    Real(dp), Allocatable, Dimension(:)    :: Gnt, Rint1, Tk, Tj, Tl, Ts, Diag, D, E
    Real(dp), Allocatable, Dimension(:,:)  :: Z1, P, W

    Real(dp), Dimension(IPs)   :: Eps
    Real(dp), Dimension(IP1)   :: C

    Type Matrix
        Integer,  Allocatable, Dimension(:) :: n, k
        Real(dp), Allocatable, Dimension(:) :: t
    End Type Matrix

    Type(Matrix) :: Hamil, Jsq

    Save
    
End Module