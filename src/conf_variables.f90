Module conf_variables

    Use params

    Implicit None
    
    Public

    Integer             :: Kexn=0, Ksig=0, Kdsig=0, K_prj=0, K_sms=0, Kw=0, kLSJ=0, KXIJ, KWeights, MaxNd0
    Integer             :: NmaxS=0, LmaxS=0, Nhint=0, NhintS=0, NgintS=0
    Integer             :: IPlv, nrd, kdavidson, num_is, Nc1, Nc_prev, Nd_prev
    Integer             :: Kherr=0, Kgerr=0, Lmax, Kmax, Ksym, Nsum, Nnr
    Integer             :: Nd0=0, vaBinSize, Nlx
    Integer(Kind=Int64) :: iscr=0_int64, NumH=0_int64, NumJ=0_int64, counter1=0_int64, ij8=0_int64, ih8=0_int64, Ngint=0_int64
    Integer(Kind=Int64) :: memTotalPerCPU=0_int64, memDvdsn=0_int64, memFormH=0_int64, memStaticArrays=0_int64, memEstimate=0_int64
    Real(dp)            :: Ecore=0_dp, XJ_av, dR_N, E_k, xscr=0_dp, K_gnt, E_0

    Integer(Kind=Int64), Allocatable, Dimension(:) :: IntOrd
    Integer,  Allocatable, Dimension(:)    :: IntOrdS, Iconverge, iconf1, iconf2, In, Nrnrc
    Integer,  Allocatable, Dimension(:)    :: Iint1, Iint2, Iint3, Iint1S, Iint2S, Iint3S, I_is
    Integer,  Allocatable, Dimension(:)    :: num_gaunts_per_partial_wave ! counts number of gaunt factors calculated in each partial wave
    Real(dp), Allocatable, Dimension(:)    :: Gnt, Rint1, Tl, Ts, D
    Real(dp), Allocatable, Dimension(:,:)  :: W

    ! Set type_real to determine whether to use single precision (sp) or double precision (dp) for Hamiltonian
    Integer, Parameter :: type_real=dp

    ! Set type2_real to determine whether to use single precision (sp) or double precision (dp) for two-electron and IS integrals
    Integer, Parameter :: type2_real=sp
    
    Real(type2_real), Allocatable, Dimension(:)    :: Rsig, Dsig, Esig, Rint2S, Dint2S, Eint2S, R_is, Scr
    Real(type2_real), Allocatable, Dimension(:,:)  :: Rint2
    
    Real(type_real),  Allocatable, Dimension(:,:) :: ArrB, P, Z1
    Real(type_real),  Allocatable, Dimension(:) :: Diag, Tj, Tk, E, E1
    Real(type_real), Allocatable, Dimension(:) :: B1, B2

    Type :: Matrix(knd)
        Integer, kind :: knd
        Real(kind=knd) :: minval
        Integer,  Allocatable, Dimension(:) :: ind1, ind2
        Real(kind=knd), Allocatable, Dimension(:) :: val
    End Type Matrix

    Type(Matrix(type_real)) :: Hamil, Jsq

    Save
    
End Module