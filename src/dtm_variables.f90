Module dtm_variables
    
    Use params

    Implicit None
    
    Public
    
    ! Global variables for dtm program
    Integer  :: Kt, Kdm, Ii, Ng, Iprt, K_M1, Kl1, Njf, Nmax, Nd1, Nd2, nterm1, nterm1f, nterm2, nterm2f
    Integer  :: Nint, Npo, nlvs
    Real(dp) :: Cl, Trd, H0, Anuc, Dint, MaxT, all, Tm2, AE1V, tj1

    Integer, Allocatable, Dimension(:)      :: Intg, iconf1, iconf2
    Integer, Allocatable, Dimension(:,:)    :: Iarr2
    Real(dp), Allocatable, Dimension(:)     :: Rnt, e2s, tj2s
    Real(dp), Allocatable, Dimension(:,:)   :: ArrB2
    Real(dp), Dimension(IP1,IP1)            :: Ro
    Real(dp), Dimension(IPx,IPx)            :: Rro
    Real(dp), Dimension(IP6)                :: C, R, V
    Character(Len=4), Dimension(2)          :: chm1
    Character(Len=1), Dimension(6)          :: Let
    Character(Len=4), Dimension(14)         :: Alet
    Character(Len=4), Dimension(5)          :: Blet

    Save

End Module