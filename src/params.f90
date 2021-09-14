Module params
    !
    ! This module contains numerical values for physical constants,
    !                      parameters determining dimensions of arrays,
    !                      and global variables and arrays used in parallel programs
    !
    Use, Intrinsic :: iso_fortran_env

    Implicit None
    
    Public

    ! double precision
    Integer, Parameter  :: dp   = real64   

    ! defining numerical values for physical constants
    Real(dp), Parameter :: DPcl = 137.0359922594e0_dp ! speed of light in a.u.
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*1.264911064067 ! =sqrt(8/5)
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*1.069044967650 ! =sqrt(8/7)
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*0.942809041582 ! =sqrt(8/9)
    Real(dp), Parameter :: DPRy = 109737.315709e0_dp  ! Rydberg constant in 1/cm
    Real(dp), Parameter :: DPau = 6.5796839e15_dp     ! au = a.u. for energy is Hz
    Real(dp), Parameter :: DPmp = 1836.153e0_dp       ! mp = proton mass in a.u.
    Real(dp), Parameter :: DPrb = 0.529177249e-8_dp   ! rb = Bohr radius in cm
    Real(dp), Parameter :: DPsw = 0.2319_dp           ! sin**2(Qeta_W)
    Real(dp), Parameter :: DPcw = 1.00050e-14_dp      ! sqrt(2)*G_F*alpha/pi (a.u.)
    Real(dp), Parameter :: DPef = 5.14225e9_dp        ! au for electric field in V/cm
    
    ! defining Parameters which determine dimensions of main arrays
    !                     Array dimension          Associated variable
    Integer, Parameter :: IPs   =    600         ! Ns - number of orbitals
    Integer, Parameter :: IPjd  =     33         ! Njd - number of possible J's   
    Integer, Parameter :: IPx   =    400         ! Nx - used for indexation of integrals   
    Integer, Parameter :: IP1   =  15000         ! Nd1 - number of determinants for 
                                                 !       direct diagonalization   
    Integer, Parameter :: IPgnt =   2891         ! Ngaunt - number of tabulated Gaunts
    Integer, Parameter :: IPbr  =      2         ! 1 - no valence Breit
                                                 ! 2 - valence Breit/Gaunt  

    ! defining parameters which determine dimensions of main arrays in hfd
    Integer, Parameter :: IP6   =    470         ! record length for DAT files
    Integer, Parameter :: IPmr  =      1         ! =4 if word=1B (Pentium)
                                                 ! =1 if word=4B (Alpha-processor)   
    
    ! Global variables 
    Integer :: Ns, Nsp, Nso, Nsu, Ne, Nec, Nc, Nc4, Nd, Nlv, Ndr, Njd, Nst, Ncpt, N_it, M, Mj
    Integer :: Kl, Kl4, Klow, Kc, Kv, Kbrt, Kout, Kecp, K_is, Kautobas, Jdel, testwigner
    Real(dp) :: Z, H, Jm, Gj, gnuc, Rnuc, Qnuc, Cut0, Crt4, C_is, Am

    ! Global arrays
    Integer, Dimension(IPjd)              :: Jt, Njt
    Integer, Dimension(IPs)               :: Nn, Kk, Ll, Jj, Nf0
    Integer, Allocatable, Dimension(:)    :: Jz, Nh, Nip, Nq, Nq0, Nc0, Ndc, Nvc, Jtc
    Real(dp), Allocatable, Dimension(:)   :: D1, E1, B1, B2, Qnl
    Integer, Allocatable, Dimension(:,:)  :: Iarr
    Real(dp), Allocatable, Dimension(:,:) :: ArrB

    Save
                                       
End Module