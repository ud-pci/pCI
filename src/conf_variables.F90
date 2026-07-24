Module conf_variables

    Use params

    Implicit None
    
    Public

    Integer             :: Kexn=0, Ksig=0, Kdsig=0, K_prj=0, K_sms=0, Kw=0, kLSJ=0, KXIJ, KWeights, MaxNd0, kCSF=0
    Integer             :: NmaxS=0, LmaxS=0, Nhint=0, NhintS=0, NgintS=0
    Integer             :: IPlv, nrd, kdavidson, num_is, Nc1, Nc_prev, Nd_prev
    Integer             :: Kherr=0, Kgerr=0, Lmax, Kmax, Ksym, Nsum, Nnr
    Integer             :: Nd0=0, vaBinSize, Nlx
    Integer(Kind=Int64) :: iscr=0_int64, NumH=0_int64, NumJ=0_int64, counter1=0_int64, ij8=0_int64, ih8=0_int64, Ngint=0_int64
    Integer(Kind=Int64) :: memTotalPerCPU=0_int64, memDvdsn=0_int64, memFormH=0_int64, memStaticArrays=0_int64, memEstimate=0_int64
    Real(dp)            :: Ecore=0_dp, XJ_av, dR_N, E_k, xscr=0_dp, K_gnt, E_0

    ! CSF variables
    Integer             :: nconf_neq=0
    Character(Len=1), Dimension(16) :: name
    Integer, Allocatable, Dimension(:)   :: ndcs, mdc, mdcs, iplace_cj, nc_neq, ndc_neq, ni_conf, nf_conf

    Integer(Kind=Int64), Allocatable, Dimension(:) :: IntOrd
    Integer,  Allocatable, Dimension(:)    :: IntOrdS, Iconverge, iconf1, iconf2, In, Nrnrc
    Integer,  Allocatable, Dimension(:)    :: Iint1, Iint2, Iint3, Iint1S, Iint2S, Iint3S, I_is
    Integer,  Allocatable, Dimension(:)    :: num_gaunts_per_partial_wave ! counts number of gaunt factors calculated in each partial wave
    Real(dp), Allocatable, Dimension(:)    :: Gnt, Rint1, Tl, Ts, D, GauntLUT
    Real(dp), Allocatable, Dimension(:,:)  :: W

    ! Integral lookup tables and hash tables
    Integer(Kind=Int64) :: GintHashCap=0_int64
    Integer             :: GintSHashCap=0
    Integer(Kind=Int64), Allocatable, Dimension(:) :: GintHashPos
    Integer,  Allocatable, Dimension(:)    :: HintSPosLUT
    Integer,  Allocatable, Dimension(:)    :: GintHashIac, GintHashIbd
    Integer,  Allocatable, Dimension(:)    :: GintSHashPos, GintSHashIac, GintSHashIbd
    Real(dp), Allocatable, Dimension(:)    :: HintLUT, ISLUT

    ! Set type_real to determine whether to use single precision (sp) or double precision (dp) for Hamiltonian
    Integer, Parameter :: type_real=dp

    ! Set type2_real to determine whether to use single precision (sp) or double precision (dp) for two-electron and IS integrals
#ifdef USE_DP_INTEGRALS
    Integer, Parameter :: type2_real=dp
#else
    Integer, Parameter :: type2_real=sp
#endif
    
    Real(type2_real), Allocatable, Dimension(:)    :: Rsig, Dsig, Esig, Rint2S, Dint2S, Eint2S, R_is, Scr
    Real(type2_real), Allocatable, Dimension(:,:)  :: Rint2
    
    Real(type_real),  Allocatable, Dimension(:,:) :: ArrB, P, Z1
    Real(type_real),  Allocatable, Dimension(:) :: Diag, Tj, Tk, E, E1
    Real(type_real), Allocatable, Dimension(:) :: B1, B1_loc, B2

    Type :: Matrix(knd)
        Integer, kind :: knd
        Real(kind=knd) :: minval
        Integer,  Allocatable, Dimension(:) :: row, col
        Real(kind=knd), Allocatable, Dimension(:) :: val
    End Type Matrix

    Type(Matrix(type_real)) :: Hamil, Jsq

    ! CSR storage for Hamil after RedistributeHamCSR.
    ! nd_start/nd_end define a rank's row range [nd_start+1 : nd_end].
    ! HamilCSR_rowptr(r) = cumulative nonzero count through local row r (0-based), so nonzeros for local row r occupy positions
    ! HamilCSR_rowptr(r-1)+1 ... HamilCSR_rowptr(r).
    Integer :: nd_start=0, nd_end=0, nd_local=0
    Integer, Allocatable, Dimension(:) :: HamilCSR_rowptr, JsqCSR_rowptr
    
    ! MPI_Allgatherv metadata for distributed ArrB column operations (set in AllocateDvdsnArrays)
    Integer, Allocatable, Dimension(:) :: nd_counts_all, nd_displs_all

    Save
    
End Module