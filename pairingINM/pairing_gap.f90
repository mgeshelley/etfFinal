program pairing_gap
    
    use routines
    implicit none
    
    
    ! Errors
    integer :: ierr
    
    
    ! Quantities for pairing gap equation
    real(kind=dp) :: energy,nu_n,m_star_m,k,dk,kf,k_low,k_high,rho,deln_guess,deln,integral
    integer :: n_grid, i
    ! Interaction
    real(kind=dp) :: v0, g, eta_n, alpha_n
    ! Ranges
    real(kind=dp) :: e_cut, kf_max, dkf
    integer :: cutoff_set
    ! Self-consistence
    real(kind=dp) :: tol, mix
    integer :: iters, tot_iters
    ! Smooth cutoff
    logical :: smooth_cutoff
    
    ! Parameter I/O
    namelist /par/ cutoff_set,e_cut,v0,eta_n,alpha_n,dk,tol,mix
    
    
    !!! PROGRAM START !!!
    
    ! Energy cutoff
!     e_cut = 120._dp ! Table 1
!     alpha = -0.55_dp ! Table 1
    
!     v0 = -448._dp ! Table 1
!     v0 = hbar_2m * 4.0_dp*pi**2*alpha
!     g = 0.5_dp ! Interaction
    
    cutoff_set = 11
    smooth_cutoff = .false.
    
    ! Self-consistence parameters
    tol = 1e-6_dp ! Tolerance
    mix = 0.5_dp ! Fraction of old Delta_n to use for next guess
    
    dk = 1e-4_dp ! Step in k in integral
    
    ! Interaction strength
    call cutoff_inter_pars(cutoff_set,v0,e_cut,eta_n,alpha_n)

    write(*,nml=par)
    
    
    open(unit=10,file='gap.dat',status='unknown',iostat=ierr)
    
    ! Fermi mtm to calculate pairing gap for
    kf = 0.0_dp ! Start density
    kf_max = 1.70_dp ! End density
    dkf = 1e-2_dp ! Density step
    
    tot_iters = 0
    
    densities: do while (kf .lt. kf_max)
        iters = 0
        
        ! Density
        rho = kf**3 / (3._dp*pi**2)
        
        ! Effective mass
        m_star_m = eff_mass(rho)
        
        ! Effective neutron chemical potential
        nu_n = hbar_2m*kf**2/m_star_m
        
        ! Interaction strength
        g = 1._dp - eta_n*((rho/rho0)**alpha_n)
        
        deln = 1._dp ! Start guess for pairing gap
        
!         print *, nu_n, g, m_star_m, rho
        
        ! First prescription for cutoff momenta
        k_low = 0._dp
        k_high = 3.0_dp ! sqrt(e_cut / (2._dp*hbar_2m/m_star_m))
        
        ! Third prescription for cutoff momenta
!         if (nu_n .lt. sqrt((e_cut/2)**2 - deln**2)) then
!         if (e_cut .gt. deln) then
!             k_low = 0._dp
!         else
!             k_low = sqrt((nu_n - sqrt((e_cut/2.)**2 - deln**2)) / (hbar_2m/m_star_m))
!         end if
!         k_high = sqrt((nu_n + sqrt((e_cut/2.)**2 - deln**2)) / (hbar_2m/m_star_m))
        
        ! Number of integration points
        n_grid = nint((k_high-k_low) / dk)
        
        self_con: do while (abs(deln-deln_guess) .gt. tol)
            iters = iters + 1 ! Track number of iterations in self-consistent loop
            tot_iters = tot_iters + 1
            
            deln_guess = deln
            
            ! Calculate integral
            integral = 0._dp
            
            ! Integration
            integration: do i=0, n_grid
                k = real(i, kind=dp) * dk
                energy = sqrt((hbar_2m*k**2/m_star_m - nu_n)**2 + deln_guess**2)
                
                if (energy .le. e_cut) then
                    integral = integral + ( dk*k**2*deln_guess / energy )
                    
                else
                    if (smooth_cutoff .eqv. .true.) then
                        integral = integral + ( dk*k**2*deln_guess / energy ) * exp(-(energy-20._dp)**2 / 100._dp)
                    end if
                    
                end if
            end do integration
            
            ! Self-consistent step
            deln = mix*deln + (1._dp-mix)* (-v0*g/(4._dp*pi**2))*integral 
            
        end do self_con
        
        write(unit=10,fmt=*) rho, kf, deln, nu_n, iters, m_star_m, k_low, k_high
        
        kf = kf + dkf
    end do densities
    
    close(unit=10,iostat=ierr)
    
    print *, tot_iters
    
end program pairing_gap
