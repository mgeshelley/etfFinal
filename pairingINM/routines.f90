module routines

    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    ! Constants
    real(kind=dp), parameter :: pi = 3.1415926535897932_dp
    real(kind=dp), parameter :: hbar_2m = 20.75_dp
    real(kind=dp), parameter :: rho0 = 0.16_dp

    contains
    
    function eff_mass(rho)
        implicit none
        
        real(kind=dp), intent(in) :: rho
        real(kind=dp) :: eff_mass
        real(kind=dp) :: t1,x1,t2,x2
!         real(kind=dp) :: tw23
        
!         tw23 = 2._dp**(2.0_dp/3)
        
!         ! SLy4 / SLy230b
!         t1 = 486.82_dp
!         x1 = -0.3438_dp
!         t2 = -546.39_dp
!         x2 = -1._dp
        
!         ! LNS
!         t1 = 266.735
!         x1 = 0.65845
!         t2 = -337.135
!         x2 = -0.95382
        
!         ! NRAPR(ii)
!         t1 = 417.64_dp
!         x1 = -0.047986_dp
!         t2 = -66.687_dp
!         x2 = 0.02717_dp
        
!         ! SKRA
!         t1 = 405.5_dp
!         x1 = 0._dp
!         t2 = -89.1_dp
!         x2 = 0.2_dp
        
!         ! SQMC700
!         t1 = 370.9804_dp
!         x1 = 0._dp
!         t2 = -96.6917_dp
!         x2 = 0._dp
        
!         ! KDE0v1
!         t1 = 403.7285_dp
!         x1 = -0.5229_dp
!         t2 = -394.5578_dp
!         x2 = -0.8956_dp
        
!         ! SII
!         t1 = 586.6_dp
!         x1 = 0._dp
!         t2 = -27.1_dp
!         x2 = 0._dp

!         ! SIV
!         t1 = 765.0_dp
!         x1 = 0._dp
!         t2 = 35.0_dp
!         x2 = 0._dp

!         ! SKa
!         t1 = 570.88_dp
!         x1 = 0._dp
!         t2 = -67.7_dp
!         x2 = 0._dp

!         ! SQMC650
!         t1 = 436.0993_dp
!         x1 = 0._dp
!         t2 = -151.9441_dp
!         x2 = 0._dp

!         ! SkM*
!         t1 = 410._dp
!         x1 = 0._dp
!         t2 = -135_dp
!         x2 = 0._dp

!         ! Skz-1
!         t1 = 439.85_dp
!         x1 = 1.2968_dp
!         t2 = -299.14_dp
!         x2 = -0.8899_dp

        ! KDE
        t1 = 403.7285_dp
        x1 = -0.5229_dp
        t2 = -394.5578_dp
        x2 = -0.8956_dp
        
!         ! Dutra et al.
!         eff_mass = (tw23 &
!                     + 0.125_dp/hbar_2m*rho*((t1*(x1+2._dp)+t2*(x2+2._dp))*tw23 &
!                     + (0.5_dp*(t2*(2._dp*x2+1._dp)-t1*(2._dp*x1+1._dp)))*tw23*2._dp))**(-1)
        
!         ! Constant effective mass = 1
!         eff_mass = 1._dp
        
!         ! LYVA1
!         t1 = 207.300_dp
!         x1 = -0.1688_dp
!         t2 = 527.930_dp
!         x2 = -1.0131_dp
        
        ! Davesne et al.
        eff_mass = (1._dp + 0.0625_dp/hbar_2m* &
                                    (3._dp*t1 + (5._dp+4._dp*x2)*t2 - (2._dp*x1+1._dp)*t1 + (2._dp*x2+1._dp)*t2)*rho)**(-1)
!         eff_mass = (1._dp + 3.9466018098455597_dp*rho)**(-1)
        
    end function eff_mass
    
    
    subroutine cutoff_inter_pars(cutoff_set,v0,e_cut,eta_n,alpha_n)
        implicit none
        
        integer,       intent(in)       :: cutoff_set
        real(kind=dp), intent(inout)    :: v0,e_cut,eta_n,alpha_n
        
        select case(cutoff_set)
            case(1)
                v0 = -448._dp
                e_cut = 120._dp
                eta_n = 0.947_dp
                alpha_n = 0.554_dp
            case(2)
                v0 = -542._dp
                e_cut = 80._dp
                eta_n = 1.01_dp
                alpha_n = 0.525_dp
            case(3)
                v0 = -746._dp
                e_cut = 40._dp
                eta_n = 1.10_dp
                alpha_n = 0.485_dp
            case(4)
                v0 = -1024._dp
                e_cut = 20._dp
                eta_n = 0.931_dp
                alpha_n = 0.378_dp
            case(5) ! From Ale paper (2012)
                v0 = -458.0_dp
                e_cut = 120._dp
                eta_n = 0.71_dp
                alpha_n = 0.51_dp
            case(6) ! From Ale HFB_INM code, DDDI sharp cutoff
                v0 = -430._dp
                e_cut = 60._dp
                eta_n = 0.75_dp
                alpha_n = 0.45_dp
            case(7) ! From Pastore et al. (2013)
                V0 = -530._dp
                e_cut = 20._dp
                eta_n = 0.7_dp
                alpha_n = 0.45_dp
            case(8) ! Modified version of case(7) to have 1/3 max value
                V0 = -355._dp
                e_cut = 20._dp
                eta_n = 0.67_dp
                alpha_n = 0.45_dp
            case(9) ! Modified DDDI sharp cutoff to have max of ~1MeV
                V0 = -362._dp
                e_cut = 60._dp
                eta_n = 0.74_dp
                alpha_n = 0.34_dp
            case(10) ! LNS fitted to have max of ~3MeV
                V0 = -446._dp
                e_cut = 60._dp
                eta_n = 0.9_dp
                alpha_n = 0.43_dp
            case(11) ! Other functionals with strength adjusted to have max of ~3MeV
                V0 = -402._dp
                e_cut = 60._dp
                eta_n = 0.7_dp
                alpha_n = 0.45_dp
        end select
    
    end subroutine cutoff_inter_pars

end module routines
