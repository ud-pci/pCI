Module ine_variables

    Use params, ipmr1 => IPmr, IP1conf => IP1

    Implicit None
    
    Public

    Integer  :: Kli, Klf, Khe, Ndir, Int_err, Ntr, Nint, Nmax, Nd0, ipmr, Kdiag, Nlft, IP1
    Integer  :: IPlv, IP4, Kt, NumJ, Nlev, N0, N2, nrange
    Integer(kind=int64) :: NumH
    Real(dp) :: Jm0, E0, E2, Tj0, Tj2, xlamb, xlamb1, xlamb2, xlambstep, XInuc, Crit1, W0, W00, Q, Elft, Hmin
    Real(dp), Allocatable, Dimension(:) :: xlamb1s, xlamb2s, xlambsteps, xlamblist

    Integer, Parameter :: IPad = 8

    Integer, Allocatable, Dimension(:) :: Int
    Real(dp), Allocatable, Dimension(:) :: Z1, X0, X1, X2, YY1, YY2, Rnt, Ev, Diag
    Real(dp), Allocatable, Dimension(:,:) :: X1J, Y2J

    Character(Len=1), Dimension(9)          :: Let
    Character(Len=4), Dimension(13)         :: Alet
    Character(Len=4), Dimension(5)          :: Blet

    Type Matrix
        Integer,  Allocatable, Dimension(:) :: n, k
        Real(dp), Allocatable, Dimension(:) ::  t
    End Type Matrix

    Type(Matrix) :: Hamil, Jsq

    Save
    
End Module