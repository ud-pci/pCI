Module basc_variables

    Use params, ipmr1 => IPmr
    Use conf_variables, Only : Ecore, Rint1, Iint1, IntOrd, Rint2, Iint2, Iint3, R_is, I_is, In, Gnt

    Implicit None

    Integer :: MaxT, Kfile, kt, kt1, ii, ipmr, nhint, ngint
    Integer, Allocatable, Dimension(:) :: Kbas
    Integer, Dimension(10) :: N_l
    Real(dp), Dimension(10) :: Rcut
    Real(dp) :: H0, Dint, Cl, Alfd, small
    Real(dp), Dimension(IPs) :: Qq
    Real(dp), Dimension(IP6) :: C, R, V, P, Q
    Real(dp), Dimension(IP6) :: Pa, Pb, Pc, Pd, Qa, Qb, Qc, Qd
    Real(dp), Dimension(IP6) :: P1a, P1b, P1c, P1d, Q1a, Q1b, Q1c, Q1d

    Integer, Parameter :: maxl = 14
    Integer, Parameter :: maxk = 2*maxl+2
    Real(dp), Dimension(IP6,2*maxk+1) :: yy, uu, vv
    Real(dp) :: r2, al, bt, rm


    Real(dp) :: xja, xjb, xjc, xjd !### total angular moments of orbitals
    Real(dp) :: xla,xlb,xlc,xld !### orbital angular moments of upper c-s
    Real(dp) :: yla,ylb,ylc,yld !### orbital angular moments of lower c-s
    Character(len=1), Dimension(9) :: let

    Public 
    
    Save

End Module basc_variables