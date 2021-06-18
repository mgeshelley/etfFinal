!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Module to hold global parameters and variables
module parameters

    implicit none
    integer, parameter :: dp = selected_real_kind(15,300) !< Precision
    
    
    !> @defgroup CONSTANTS Global parameters
    !> @{
    real(kind=dp), parameter :: pi = 3.1415926535897932_dp !< \f$\pi\f$
    real(kind=dp), parameter :: hbar2_2m = 20.73553_dp !< \f$\frac{\hbar^2}{2m}\f$ for SLy forces
    real(kind=dp), parameter :: rho0 = 0.16_dp !< \f$\rho_0\f$, nuclear saturation density \f$[fm^{-3}]\f$
    real(kind=dp), parameter :: e2 = 1.439978408596513_dp !< \f$e^2\f$, electric charge squared \f$[MeV\cdot fm]\f$
    real(kind=dp), parameter :: mec2 = 0.51099895_dp !< Electron rest mass \f$[MeV/c^2]\f$
    real(kind=dp), parameter :: mnc2 = 939.56542052_dp ! < Neutron rest mass \f$[MeV/c^2]\f$
    real(kind=dp), parameter :: mpc2 = 938.27208816_dp ! < Proton rest mass \f$[MeV/c^2]\f$
    real(kind=dp), parameter :: hbarc = 197.3269804_dp !< \f$\hbar*c\f$
    !> @}
    
    !> @defgroup HBAR2_2M Extra constants for calculating \f$\frac{\hbar^2}{2m}\f$ for BSk forces
    !> @{
    real(kind=dp), parameter :: mamuc2 = 931.49386_dp !< Atomic mass unit \f$u\f$
    real(kind=dp), parameter :: xmh = 7.28896940_dp
    real(kind=dp), parameter :: rydb = 13.6056981e-6_dp !< Rydberg constant
    real(kind=dp), parameter :: xmp = xmh - mec2 + rydb
    real(kind=dp), parameter :: xmn = 8.07132281_dp
    !> @}
    
    !> @defgroup FRACTIONS Numerical fractions
    !> @{
    real(kind=dp), parameter :: fr12 = 1._dp/2, fr14 = 1._dp/4, fr18 = 1._dp/8, fr1_16 = 1._dp/16, fr1_32 = 1._dp/32
    real(kind=dp), parameter :: fr13 = 1._dp/3, fr16 = 1._dp/6, fr1_12 = 1._dp/12, fr1_24 = 1._dp/24, fr1_36 = 1._dp/36
    real(kind=dp), parameter :: fr23 = 2._dp/3, fr29 = 2._dp/9
    real(kind=dp), parameter :: fr32 = 3._dp/2, fr34 = 3._dp/4, fr35 = 3._dp/5, fr38 = 3._dp/8, fr3_10 = 3._dp/10
    real(kind=dp), parameter :: fr43 = 4._dp/3, fr53 = 5._dp/3, fr54 = 5._dp/4, fr83 = 8._dp/3, fr260_3 = 260._dp/3
    !> @}
    
    
    !> @defgroup SKYRME_PARS Skyrme parameters
    !> @{
    real(kind=dp)                   :: W0 !< Spin-orbit strength \f$W_0\f$
    real(kind=dp)                   :: t0,x0,t1,x1,t2,x2,t3,x3,sigma
    real(kind=dp),  dimension(0:1)  :: hbar2_2m_q !< \f$\frac{\hbar^2}{2m}\f$ with isospin dependence, as required by BSk forces
    logical                         :: J2_terms !< Whether functional uses \f$\textbf{J}^2\f$ terms
    !> @}
    
    !> @defgroup SKYRME_COEFFS Skyrme coefficients
    !> @{
    real(kind=dp) :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13
    !> @}
    
    !> @defgroup BSK_EXTRA_PARS Extra parameters for BSk forces
    !> @{
    real(kind=dp) :: alpha,beta,gamma
    real(kind=dp) :: t4,x4,t5,x5
    real(kind=dp) :: fn_pos,fn_neg,fp_pos,fp_neg,epsilon_Lambda
    !> @}
    
    
    !> @defgroup MESH Mesh variable
    !> @{
    real(kind=dp),  dimension(:),   allocatable :: r !< Mesh for r
    !> @}
    
    !> @defgroup DENSITIES Global density arrays
    !> @{
    real(kind=dp),  dimension(:,:), allocatable :: rho_q, del_rho_q, del2_rho_q
    real(kind=dp),  dimension(:),   allocatable :: rho_t, del_rho_t, del2_rho_t
    !> @}
    
    !> @defgroup EFF_MASSES Global effective mass arrays
    !> @{
    real(kind=dp),  dimension(:,:), allocatable :: f_q, del_f_q, del2_f_q
    !> @}
    
    !> @defgroup OTHER_DENSITIES Global arrays for other densities and fields
    !> @{
    real(kind=dp),  dimension(:,:),     allocatable :: W_q, div_W_q, U_q, J_2_q, J_q, del_J_q, div_J_2_q, div_J_q
    real(kind=dp),  dimension(:,:),     allocatable :: tau_TF_q, tau_2_L_q, tau_2_NL_q, tau_ETF_q
    real(kind=dp),  dimension(:),       allocatable :: J_t, del_J_t, div_J_t, tau_ETF_t
    real(kind=dp),  dimension(:),       allocatable :: e_density_field, e_density_sky
    real(kind=dp),  dimension(:),       allocatable :: v_c_di, v_c_ex, e_density_c_di, e_density_c_ex, v_c_pe
    real(kind=dp),  dimension(:,:,:),   allocatable :: tau_2_cont_q
    real(kind=dp),  dimension(:),       allocatable :: rho_ch !< Charge density
    !> @}
    
    !> @defgroup FOURTH_ORDER Global arrays for densities and fields for 4th-order terms
    real(kind=dp),  dimension(:,:),     allocatable :: d2_rho_q, d3_rho_q, d4_rho_q, d2_f_q, d3_f_q, d4_f_q
    real(kind=dp),  dimension(:,:),     allocatable :: A_q, d1_A_q, d2_A_q, d3_A_q, d4_A_q
    real(kind=dp),  dimension(:,:),     allocatable :: tau_4_no_spin_q, tau_4_so_q, J_4_q, div_J_4_q
    real(kind=dp),  dimension(:,:,:),   allocatable :: tau_4_cont_q
    !> @}
    
    !> @defgroup WS_PROPERTIES Variables for storing properties of Wigner-Seitz cell
    !> @{
    real(kind=dp),  dimension(0:1)  :: N_q !< Number of neutrons and protons
    real(kind=dp)                   :: N_t !< Total number of particles
    real(kind=dp)                   :: e_field !< Field energy
    real(kind=dp)                   :: e_skyrme !< Skyrme energy
    real(kind=dp),  dimension(0:1)  :: e_kinetic_q !< Kinetic energy for neutron and protons
    real(kind=dp)                   :: e_kinetic_t !< Total kinetic energy
    real(kind=dp)                   :: e_so !< Spin-orbit energy
    real(kind=dp)                   :: e_coulomb_di !< Direct Coulomb energy
    real(kind=dp)                   :: e_coulomb_ex !< Exhange Coulomb energy
    real(kind=dp)                   :: e_coulomb !< Total Coulomb energy
    real(kind=dp)                   :: e_total !< Total energy
    real(kind=dp),  dimension(0:1)  :: mu_q !< Neutron and proton chemical potential
    real(kind=dp)                   :: mu_e !< Electron chemical potential
    real(kind=dp)                   :: mu_c !< Coulomb interaction contribution to electron chemical potential
    real(kind=dp)                   :: delta_mu !< Overall chemical potential (beta-equilibrium condition)
    real(kind=dp)                   :: pressure_nucl !< Nuclear pressure
    real(kind=dp)                   :: pressure_e !< Electron pressure
    real(kind=dp)                   :: pressure_ex !< Coulomb exchange pressure
    real(kind=dp)                   :: pressure_t !< Total pressure
    !> @}
    
    
    !> @defgroup INPUT_PARS Input parameters
    !> @{
    integer                         :: run_mode !< Mode to run code in (normal, test)
    logical                         :: verbose !< Whether to print all extra messages about run
    real(kind=dp)                   :: dr !< Mesh spacing of r \f$(fm)\f$
    logical                         :: profile_r_max !< Whether to take "r_max" from "r_ws"
    real(kind=dp)                   :: r_max  !< Max value of r \f$(fm)\f$
    logical                         :: specify_n !< Whether to specify "n" instead of "dr"
    integer                         :: n !< Number of mesh points
    integer                         :: force !< Force to use
    integer                         :: etf_order !< Order at which to calculate kinetic energy densities
    logical                         :: coulomb_on !< Use Coulomb interaction or not
    logical                         :: electrons_on !< Add electrons
    integer                         :: Nmaxstate !< Number of states that can be stored
    integer                         :: Lmax !< Maximum angular momentum of states to find
    logical                         :: calc_chem_pots !< Calculate chemical potentials
    logical,        dimension(0:1)  :: strutinsky_on !< Use Strutinsky Integral (SI) correction for neutrons and protons
    real(kind=dp)                   :: strut_r_max !< Max value of r to use for box for single particle states \f$(fm)\f$
    logical                         :: emax0_mod !< Whether to modify "emax0" to always = 0
    integer                         :: neutron_pairing !< Type of pairing calculation to perform for neutrons
    logical                         :: proton_bcs !< Whether to do BCS for protons
    real(kind=dp)                   :: pair_qp_cut !< Cut-off for quasiparticle energy for solving gap equation
    logical                         :: smooth_qp_cut !< Whether to use smooth cut-off on quasiparticle energy
    real(kind=dp)                   :: pair_k_max !< Maximum momentum for integral in gap equation
    real(kind=dp)                   :: pair_v0 !< Interaction strength parameter
    real(kind=dp)                   :: pair_eta !< Interaction parameter eta
    real(kind=dp)                   :: pair_alpha !< Interaction parameter alpha
    real(kind=dp)                   :: pair_tol !< Gap equation self-consistency tolerance
    real(kind=dp)                   :: pair_mix !< Mix of previous iteration in gap equation self-consistency loop
    real(kind=dp)                   :: pair_dk !< Step size in integral in gap equation
    real(kind=dp)                   :: pair_del_init !< Initial guess for pairing gap
    integer                         :: pair_max_iters !< Maximum iterations for gap equation self-consistency loop 
    real(kind=dp)                   :: pair_num_tol !< BCS number equation self-consistency tolerance
    !> @}
    
    !> @defgroup EXT_PROFILES External profiles
    !> @{
    real(kind=dp),  dimension(5)    :: n_profile, p_profile
    real(kind=dp)                   :: num_n
    integer                         :: num_p
    real(kind=dp)                   :: r_ws
    !> @}
    
    
    !> @defgroup FORMATS Format statements
    !> @{
    character(18),  dimension(10)   :: n_floats = (/'( 1(es24.16e3,1x))','( 2(es24.16e3,1x))','( 3(es24.16e3,1x))', &
                                        & '( 4(es24.16e3,1x))','( 5(es24.16e3,1x))','( 6(es24.16e3,1x))','( 7(es24.16e3,1x))', &
                                        & '( 8(es24.16e3,1x))','( 9(es24.16e3,1x))','(10(es24.16e3,1x))'/)
    character(27)                   :: info_format = '(a27,1x,a1,1x,es24.16e3)'
    !> @}
    
    !> @defgroup STRINGS Global strings for printing output
    !> @{
    character(8)                                :: force_string
    character(80),                  parameter   :: line_break = '########################################&
                                                                &########################################' !< Separating output
    character(62),                  parameter   :: table_string = '----------------------------+-------------------------'
    character(1),   dimension(0:1), parameter   :: iso_string = (/'n','p'/) !< Isospin labels
    character(4),                   parameter   :: space4 = '    ' !< Whitespace of length 4
    !> @}
    
    
    !> @defgroup STRUTINSKY Variables for use in strutinsky module
    !> @{
    ! Mesh points for use in "boundary"
    integer                                     :: strut_n
    ! Inputs and fields
    real(kind=dp),  dimension(:,:), allocatable :: dhmen,d2hmen,hb2m,vpot,vso
    real(kind=dp),  dimension(0:1)              :: Ecut !< Cutoff energies for neutrons and protons
    ! Outputs
    real(kind=dp),  dimension(:,:), allocatable :: PS !< Storage for single particle wavefunctions
    integer                                     :: NNST !< Number of states found within cutoff energy
    integer,        dimension(:),   allocatable :: JJP !< Storage for J of states found
    integer,        dimension(:),   allocatable :: LP !< Storage for L of states found
    real(kind=dp),  dimension(:),   allocatable :: EP !< Storage for energy of states found
    ! Shell correction
    real(kind=dp),  dimension(0:1)              :: e_sc_q !< Storage for shell correction energies
    real(kind=dp)                               :: e_sc_t !< Storage for total shell correction energy
    !> @}
    
    !> @defgroup PAIRING Variables for use in pairing module
    !> @{
    ! Neutrons
    real(kind=dp),  dimension(:),   allocatable :: pair_gap_n !< Neutron pairing gap
    real(kind=dp),  dimension(:),   allocatable :: delta_n !< Neutron pairing field
    ! Proton BCS
    real(kind=dp),  dimension(:),   allocatable :: pair_gap_p !< Proton pairing gap
    real(kind=dp),  dimension(:),   allocatable :: delta_p !< Proton pairing field
    real(kind=dp),  dimension(:),   allocatable :: rho_p_bcs !< Proton density (from wavefunctions)
    real(kind=dp),  dimension(:),   allocatable :: rho_anom_p_bcs !< Anomalous density
    real(kind=dp),  dimension(:),   allocatable :: strength_bcs !< Interaction strength for protons
    real(kind=dp),  dimension(:),   allocatable :: occ_v !< Particle occupation probabilities
    real(kind=dp),  dimension(:),   allocatable :: occ_u !< Hole occupation probabilities
    real(kind=dp),  dimension(:),   allocatable :: e_qp_p !< Quasiparticle energies
    real(kind=dp)                               :: num_p_bcs !< Number of protons calculated from occupations
    real(kind=dp)                               :: e_pair_bcs !< BCS pairing energy
    ! Energies
    real(kind=dp),  dimension(0:1)              :: e_pair_q !< Storage for pairing condensation energies
    real(kind=dp)                               :: e_pair_t !< Storage for total pairing condensation energy
    !> @}
    
    !> @defgroup ELECTRONS Variables for electron contributions
    !> @{
    real(kind=dp) :: rho_e !< Electron density
    real(kind=dp) :: e_kinetic_e !< Total electron kinetic energy
    real(kind=dp) :: e_coulomb_e !< Total electron-electron potential energy from Coulomb interaction
    real(kind=dp) :: e_coulomb_pe !< Total proton-electron potential energy from Coulomb interaction
    !> @}
    
    
    ! Namelists
    namelist /params/ run_mode,verbose,dr,profile_r_max,r_max,specify_n,n, &
                    & force,etf_order,coulomb_on,electrons_on, &
                    & Nmaxstate,Lmax, &
                    & calc_chem_pots, &
                    & strutinsky_on,strut_r_max,emax0_mod, &
                    & neutron_pairing,proton_bcs, &
                    & pair_qp_cut,smooth_qp_cut,pair_k_max,pair_v0,pair_eta,pair_alpha, &
                    & pair_tol,pair_mix,pair_dk,pair_del_init,pair_max_iters,pair_num_tol !< Input parameter namelist
    
    namelist /profiles/ n_profile, p_profile, num_n, num_p, r_ws !< External density profiles namelist
    
    
    contains
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to initialise Skyrme force parameters with values for a given
    !> force specified in 'input.dat'.
    subroutine force_initialise()
        implicit none
        
        select case(force)
            ! SLy4
            case(1)
                force_string = 'SLy4'
                t0 = -2488.913_dp
                t1 = 486.818_dp
                t2 = -546.395_dp
                t3 = 13777._dp
                x0 = 0.834_dp
                x1 = -0.3438_dp
                x2 = -1._dp
                x3 = 1.354_dp
                sigma = fr16
                W0 = 123._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SLy7
            case(2)
                force_string = 'SLy7'
                t0 = -2482.41_dp
                t1 = 457.97_dp
                t2 = -419.85_dp
                t3 = 13677._dp
                x0 = 0.846_dp
                x1 = -0.511_dp
                x2 = -1._dp
                x3 = 1.391_dp
                sigma = fr16
                W0 = 126._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SkM*
            case(3)
                force_string = 'SkM*'
                t0 = -2645._dp
                t1 = 410._dp
                t2 = -135_dp
                t3 = 15595._dp
                x0 = 0.09_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 0._dp
                sigma = fr16
                W0 = 120._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SkM*, with W0 = 130.0
            case(4)
                force_string = 'SkM*_130'
                t0 = -2645._dp
                t1 = 410._dp
                t2 = -135_dp
                t3 = 15595._dp
                x0 = 0.09_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 0._dp
                sigma = fr16
                W0 = 130._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! BSk14
            case(5)
                force_string = 'BSk14'
                t0 = -1822.67_dp
                t1 = 377.470_dp
                t2 = -2.41056_dp
                t3 = 11406.3_dp
                x0 = 0.302096_dp
                x1 = -0.823575_dp
                x2 = 61.9411_dp
                x3 = 0.473460_dp
                sigma = fr3_10 ! Written as "gamma" in original paper
                W0 = 135.565_dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk17
            case(6)
                force_string = 'BSk17'
                t0 = -1837.33_dp
                t1 = 389.102_dp
                t2 = -3.1742_dp
                t3 = 11523.8_dp
                x0 = 0.411377_dp
                x1 = -0.832102_dp
                x2 = 49.4875_dp
                x3 = 0.654962_dp
                sigma = fr3_10 ! Written as "gamma" in original paper
                W0 = 145.885_dp
                fn_pos = 1.00_dp
                fn_neg = 1.04_dp
                fp_pos = 1.05_dp
                fp_neg = 1.05_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk19
            case(7)
                force_string = 'BSk19'
                t0 = -4115.21_dp
                t1 = 403.072_dp
                t2 = 1.e-3_dp
                t3 = 23670.4_dp
                t4 = -60._dp
                t5 = -90._dp
                x0 = 0.398848_dp
                x1 = -0.137960_dp
                x2 = -1055550_dp
                x3 = 0.375201_dp
                x4 = -6.0_dp
                x5 = -13.0_dp
                W0 = 110.802_dp
                alpha = fr1_12
                beta = fr13
                gamma = fr1_12
                fn_pos = 1.00_dp
                fn_neg = 1.05_dp
                fp_pos = 1.10_dp
                fp_neg = 1.17_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk20
            case(8)
                force_string = 'BSk20'
                t0 = -4056.04_dp
                t1 = 438.219_dp
                t2 = 1.e-3_dp
                t3 = 23256.6_dp
                t4 = -100.000_dp
                t5 = -120.000_dp
                x0 = 0.56913_dp
                x1 = -0.392047_dp
                x2 = -1147640_dp
                x3 = 0.614276_dp
                x4 = -3.00000_dp
                x5 = -11.0000_dp
                W0 = 110.228_dp
                alpha = fr1_12
                beta = fr16
                gamma = fr1_12
                fn_pos = 1.00_dp
                fn_neg = 1.06_dp
                fp_pos = 1.09_dp
                fp_neg = 1.16_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk21
            case(9)
                force_string = 'BSk21'
                t0 = -3961.39_dp
                t1 = 396.131_dp
                t2 = 1.e-3_dp
                t3 = 22588.2_dp
                t4 = -100.000_dp
                t5 = -150.000_dp
                x0 = 0.885231_dp
                x1 = 0.0648452_dp
                x2 = -1390380_dp
                x3 = 1.03928_dp
                x4 = 2.00000_dp
                x5 = -11.0000_dp
                W0 = 109.622_dp
                alpha = fr1_12
                beta = fr12
                gamma = fr1_12
                fn_pos = 1.00_dp
                fn_neg = 1.05_dp
                fp_pos = 1.07_dp
                fp_neg = 1.13_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk22
            case(10)
                force_string = 'BSk22'
                t0 = -3978.97_dp
                t1 = 404.461_dp
                t2 = 1.e-3_dp
                t3 = 22704.7_dp
                t4 = -100.000_dp
                t5 = -150.000_dp
                x0 = 0.472558_dp
                x1 = 0.0627540_dp
                x2 = -1396130_dp
                x3 = 0.514386_dp
                x4 = 2.00000_dp
                x5 = -11.0000_dp
                W0 = 111.109_dp
                alpha = fr1_12
                beta = fr12
                gamma = fr1_12
                fn_pos = 1.00_dp
                fn_neg = 1.05_dp
                fp_pos = 1.07_dp
                fp_neg = 1.13_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! BSk24
            case(11)
                force_string = 'BSk24'
                t0 = -3970.28718809175_dp
                t1 = 395.7662092819614_dp
                t2 = 1.e-5_dp
                t3 = 22648.57756432157_dp
                t4 = -100.000_dp
                t5 = -150.000_dp
                x0 = 0.8943714063887188_dp
                x1 = 0.05635350781437905_dp
                x2 = -138961080.022762_dp
                x3 = 1.051191623994705_dp
                x4 = 2.00000_dp
                x5 = -11.0000_dp
                W0 = 108.4052877131591_dp
                alpha = fr1_12
                beta = fr12
                gamma = fr1_12
                fn_pos = 1.00_dp
                fn_neg = 1.06_dp
                fp_pos = 1.09_dp
                fp_neg = 1.16_dp
                epsilon_Lambda = 16._dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            ! LNS
            case(12)
                force_string = 'LNS'
                t0 = -2484.97_dp
                t1 = 266.735_dp
                t2 = -337.135_dp
                t3 = 14588.2_dp
                x0 = 0.06277_dp
                x1 = 0.65845_dp
                x2 = -0.95382_dp
                x3 = -0.03413_dp
                sigma = fr16
                W0 = 96._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! FPLyon
            case(13)
                force_string = 'FPLyon'
                t0 = -2498.9_dp
                t1 = 382.19_dp
                t2 = -336.96_dp
                t3 = 15230.5_dp
                x0 = 0.5469_dp
                x1 = -0.7624_dp
                x2 = -0.6813_dp
                x3 = 0.8094_dp
                sigma = 0.18832_dp
                W0 = 119.58_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SQMC650
            case(14)
                force_string = 'QMC650'
                t0 = -2462.6784_dp
                t1 = 436.0993_dp
                t2 = -151.9441_dp
                t3 = 14154.4809_dp
                x0 = 0.13_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 0._dp
                sigma = fr16
                W0 = 110.5048_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SQMC700
            case(15)
                force_string = 'QMC700'
                t0 = -2429.1323_dp
                t1 = 370.9804_dp
                t2 = -96.6917_dp
                t3 = 13773.6340_dp
                x0 = 0.1_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 0._dp
                sigma = fr16
                W0 = 104.585_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SII
            case(16)
                force_string = 'SII'
                t0 = -1169.9_dp
                t1 = 586.6_dp
                t2 = -27.1_dp
                t3 = 9331.1_dp
                x0 = 0.34_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 1._dp
                sigma = 1._dp
                W0 = 105._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SIV
            case(17)
                force_string = 'SIV'
                t0 = -1205.6_dp
                t1 = 765.0_dp
                t2 = 35.0_dp
                t3 = 5000.0_dp
                x0 = 0.05_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = 1._dp
                sigma = 1._dp
                W0 = 150._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SKa
            case(18)
                force_string = 'SKa'
                t0 = -1602.78_dp
                t1 = 570.88_dp
                t2 = -67.7_dp
                t3 = 8000.0_dp
                x0 = -0.02_dp
                x1 = 0._dp
                x2 = 0._dp
                x3 = -0.286_dp
                sigma = fr13
                W0 = 125._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SKz-1
            case(19)
                force_string = 'SKz-1'
                t0 = -2471.1_dp
                t1 = 439.85_dp
                t2 = -299.14_dp
                t3 = 13732.8_dp
                x0 = -0.2665_dp
                x1 = 1.2968_dp
                x2 = -0.8899_dp
                x3 = -0.7282_dp
                sigma = 0.1694_dp
                W0 = 120._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .false.
                
            ! SLy5
            case(20)
                force_string = 'SLy5'
                t0 = -2483.45_dp
                t1 = 484.23_dp
                t2 = -556.69_dp
                t3 = 13757.0_dp
                x0 = 0.776_dp
                x1 = -0.317_dp
                x2 = -1._dp
                x3 = 1.263_dp
                sigma = fr16
                W0 = 125._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! KDE
            case(21)
                force_string = 'KDE'
                t0 = -2532.8842_dp
                t1 = 403.7285_dp
                t2 = -394.5578_dp
                t3 = 14575.0234_dp
                x0 = 0.7707_dp
                x1 = -0.5229_dp
                x2 = -0.8956_dp
                x3 = 1.1716_dp
                sigma = 0.169_dp
                W0 = 128.0572_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! KDE0v1
            case(22)
                force_string = 'KDE0v1'
                t0 = -2553.0843_dp
                t1 = 411.6963_dp
                t2 = -419.8712_dp
                t3 = 14603.6069_dp
                x0 = 0.6483_dp
                x1 = -0.3472_dp
                x2 = -0.9268_dp
                x3 = 0.9475_dp
                sigma = 0.1673_dp
                W0 = 124.42_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! NRAPR
            case(23)
                force_string = 'NRAPR'
                t0 = -2719.7_dp
                t1 = 417.64_dp
                t2 = -66.687_dp
                t3 = 15042.0_dp
                x0 = 0.16154_dp
                x1 = -0.047986_dp
                x2 = 0.02717_dp
                x3 = 0.13611_dp
                sigma = 0.14416_dp
                W0 = 83.916_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! ska25s20
            case(24)
                force_string = 'ska25s20'
                t0 = -2180.48_dp
                t1 = 281.49_dp
                t2 = -160.44_dp
                t3 = 14577.8_dp
                x0 = 0.13712_dp
                x1 = -0.8_dp
                x2 = 0._dp
                x3 = 0.06208_dp
                sigma = fr14
                W0 = 206.4_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! ska45s20
            case(25)
                force_string = 'ska45s20'
                t0 = -1537.89_dp
                t1 = 245.47_dp
                t2 = -157.2_dp
                t3 = 12461.12_dp
                x0 = 0.12318_dp
                x1 = -0.8_dp
                x2 = 0._dp
                x3 = -0.05608_dp
                sigma = 0.45_dp
                W0 = 196.8_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! SKRA
            case(26)
                force_string = 'SKRA'
                t0 = -2895.4_dp
                t1 = 405.5_dp
                t2 = -89.1_dp
                t3 = 16660.0_dp
                x0 = 0.08_dp
                x1 = 0._dp
                x2 = 0.2_dp
                x3 = 0._dp
                sigma = 0.1422_dp
                W0 = 129._dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! NRAPRii
            case(27)
                force_string = 'NRAPRii'
                t0 = -2719.7_dp
                t1 = 417.64_dp
                t2 = -66.687_dp
                t3 = 15042.0_dp
                x0 = 0.16154_dp
                x1 = -0.047986_dp
                x2 = 0.02717_dp
                x3 = 0.13611_dp
                sigma = 0.14416_dp
                W0 = 167.832_dp
                hbar2_2m_q(:) = 20.73553_dp
                J2_terms = .true.
                
            ! BSk13
            case(28)
                force_string = 'BSk13'
                t0 = -1773.86_dp
                t1 = 307.766_dp
                t2 = -73.428_dp
                t3 = 12254.3_dp
                x0 = 0.253076_dp
                x1 = -0.471466_dp
                x2 = 1.27406_dp
                x3 = 0.215541_dp
                sigma = fr13
                W0 = 119.44_dp
                hbar2_2m_q(0) = hbarc**2 / (2._dp * (mamuc2 + xmn))
                hbar2_2m_q(1) = hbarc**2 / (2._dp * (mamuc2 + xmp))
                J2_terms = .false.
                
            case default
                print '(/a)', 'STOP: no force selected'
                stop
            
        end select
        
        if (verbose .eqv. .true.) then
            print '(/a)', 'Force = '//force_string
        end if
        
        ! Skyrme coefficients as defined in Bonche et al., Nucl. Phys. A (1987)
        B1 = fr12*t0 * (1._dp + fr12*x0)
        B2 = -fr12*t0 * (fr12 + x0)
        B3 = fr14 * ( t1*(1._dp + fr12*x1) + t2*(1._dp + fr12*x2) )
        B4 = -fr14 * ( t1*(fr12 + x1) - t2*(fr12 + x2) )
        B5 = -fr1_16 * ( 3._dp*t1*(1._dp + fr12*x1) - t2*(1._dp + fr12*x2) )
        B6 = fr1_16 * ( 3._dp*t1*(fr12 + x1) + t2*(fr12 + x2) )
        B7 = fr1_12 * t3 * (1._dp + fr12*x3)
        B8 = -fr1_12 * t3 * (fr12 + x3)
        B9 = -fr12*W0
        B10 = fr14*t0*x0
        B11 = -fr14*t0
        B12 = fr1_24*t3*x3
        B13 = -fr1_24*t3
        
    end subroutine force_initialise
    
    
end module parameters