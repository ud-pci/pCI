Module params
    !
    ! This module contains numerical values for physical constants,
    !                      parameters determining dimensions of arrays,
    !                      and global variables and arrays used in parallel programs
    !
    Use, Intrinsic :: iso_fortran_env, Only : sp => real32, dp => real64, int64
    Use mpi_f08, Only : MPI_Datatype

    Implicit None
    
    Public

    ! defining numerical values for physical constants
    Real(dp), Parameter :: DPcl = 137.0359922594_dp    ! speed of light in a.u.
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*1.264911064067 ! =sqrt(8/5)
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*1.069044967650 ! =sqrt(8/7)
    !Real(dp), Parameter :: DPcl = 137.0359922594d0*0.942809041582 ! =sqrt(8/9)
    Real(dp), Parameter :: DPRy = 109737.315709_dp     ! Rydberg constant in 1/cm
    Real(dp), Parameter :: DPau = 6.5796839e15_dp     ! au = a.u. for energy is Hz
    Real(dp), Parameter :: DPmp = 1836.153e0_dp       ! mp = proton mass in a.u.
    Real(dp), Parameter :: DPrb = 0.529177249e-8_dp   ! rb = Bohr radius in cm
    Real(dp), Parameter :: DPsw = 0.2319_dp           ! sin**2(Qeta_W)
    Real(dp), Parameter :: DPcw = 1.00050e-14_dp      ! sqrt(2)*G_F*alpha/pi (a.u.)
    Real(dp), Parameter :: DPef = 5.14225e9_dp        ! au for electric field in V/cm
    Real(dp), Parameter :: DPma = 1.15965218046e-3_dp ! electron magnetic moment anomaly
    
    ! defining Parameters which determine dimensions of main arrays
    !                     Array dimension          Associated variable
    Integer, Parameter :: IPs   =    600         ! Ns - number of orbitals
    Integer, Parameter :: IPx   =    440         ! Nx - used for indexation of integrals   
    Integer, Parameter :: IP1   =   3000         ! Nd1 - number of determinants for direct diagonalization   
    Integer, Parameter :: IPbr  =      2         ! 1 - no valence Breit, 2 - valence Breit/Gaunt  

    ! defining parameters which determine dimensions of main arrays in hfd
    Integer, Parameter :: IP6   =    470         ! record length for DAT files
    Integer, Parameter :: IPmr  =      1         ! word length in direct access files
                                                 ! =4 if word=1B (Pentium)
                                                 ! =1 if word=4B (Alpha-processor)   

    ! defining parameters which determine dimensions of main arrays in add
    Integer, Parameter :: IPad1 = 1000000        ! max number of configurations

    ! Set type_real to determine whether to use single precision (sp) or double precision (dp) for Hamiltonian
    Integer, Parameter :: type_real=dp
    Type(MPI_Datatype) :: mpi_type_real

    ! Set type2_real to determine whether to use single precision (sp) or double precision (dp) for two-electron and IS integrals
    Integer, Parameter :: type2_real=sp
    Type(MPI_Datatype) :: mpi_type2_real

    ! Set key to turn on extra outputs for developing
    Integer, Parameter :: devmode=1
    
    ! Global variables 
    Integer :: Ns, Nsp, Nso, Nsu, Ne, Nec, Nc, Nc4, Nd, Nlv, Ndr, Njd, Nst, Ncpt, N_it, Ngaunt, M, Mj
    Integer :: Kl, Kl4, Klow, Kc, Kv, Kbrt, Kout, Kecp, K_is, Kautobas, Jdel
    Real(dp) :: Z, H, Gj, gnuc, Rnuc, Qnuc, Cut0, C_is, Am, Jm, Crt4

    ! Global arrays
    Integer, Allocatable, Dimension(:)    :: Jt, Njt
    Integer, Dimension(IPs)               :: Nn, Kk, Ll, Jj, Nf0
    Integer, Allocatable, Dimension(:)    :: Jz, Nh, Nh0, Nip, Nq, Nq0, Nc0, Ndc, Nvc, Jtc
    Real(dp), Allocatable, Dimension(:)   :: D1, Qnl
    Integer, Allocatable, Dimension(:,:)  :: Iarr

    Real(type_real), Allocatable, Dimension(:) :: B1, B2

    Save
                                       
End Module