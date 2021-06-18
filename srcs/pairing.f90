!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Module to hold routines for carrying out Strutinsky correction
module pairing
    
    use parameters
    use routines
    use test
    use strutinsky
    use omp_lib
    implicit none
    
    
    contains
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the pairing energy for the WS cell, using the
    !> method specified in "input.in"
    subroutine calc_e_pair()
        implicit none
        real(kind=dp) :: mu !< Effective chemical potential
        integer :: i !< Loop index
        
        select case(neutron_pairing)
            ! Local density approximation (LDA); pairing value for neutron gas
            ! used for entire WS cell
            case(1)
                ! Calculate neutron chemical potential using density and mass
                ! values at edge of cell (index "n") (just equals Fermi energy)
                mu = calc_e_f_q(rho_q(n,0),f_q(n,0),0)
                
                ! Calculate gap for infinite nuclear matter
                pair_gap_n(n) = calc_gap(r(n),rho_q(n,0),rho_q(n,1),f_q(n,0),mu)
                
                ! Calculate LDA approximation for pairing field
                delta_n(1:n) = -fr38 * pair_gap_n(n)**2 / mu * rho_q(1:n,0)
                
                ! Integrate LDA pairing field across WS cell to get total pairing energy
                e_pair_q(0) = WS_integral(delta_n(:),n)
                
            ! Local density approximation (LDA); pairing field calculated for
            ! neutron density at each point in WS cell
            case(2)
                densities: do i=1, n
                    ! If "exact same" density and effective mass as previous
                    ! mesh point, use same value for pairing field and gap
                        ! Tolerance of epsilon*10^4 for "exact same" preserves
                        ! "correct" value of neutron pairing correction to 1keV
                    if ((i .gt. 1) .and. &
                            & (abs(rho_q(i,0) - rho_q(i-1,0)) .le. epsilon(rho_q(i,0))*1.e4_dp) .and. & 
                            & (abs(f_q(i,0) - f_q(i-1,0)) .le. epsilon(f_q(i,0))*1.e4_dp)) then
                        ! Values from previous density
                        pair_gap_n(i) = pair_gap_n(i-1)
                        delta_n(i) = delta_n(i-1)
                        
                        ! Skip to next iteration
                        cycle densities
                    end if
                    
                    ! Calculate neutron chemical potential
                    mu = calc_e_f_q(rho_q(i,0),f_q(i,0),0)
                
                    ! Calculate gap for infinite nuclear matter (just equals
                    ! Fermi energy)
                    pair_gap_n(i) = calc_gap(r(i),rho_q(i,0),rho_q(i,1),f_q(i,0),mu)
                    
                    ! Calculate LDA approximation for pairing field
                    delta_n(i) = -fr38 * pair_gap_n(i)**2 / mu * rho_q(i,0)
                end do densities
                
                ! Integrate LDA pairing field across WS cell to get total pairing energy 
                e_pair_q(0) = WS_integral(delta_n(:),n)
                
            ! If nothing selected, exit
            case default
                print '(a)', 'No pairing mode selected'
                stop
                
        end select
        
    end subroutine calc_e_pair
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the pairing gap in infinite neutron matter, for
    !> specified density and effective mass
    !> @param[in] r     Radius
    !> @param[in] rhon  Neutron density
    !> @param[in] rhop  Proton density
    !> @param[in] fq    Effective mass
    !> @param[in] mu    Effective chemical potential
    function calc_gap(r,rhon,rhop,fq,mu)
        implicit none
        real(kind=dp) :: calc_gap
        real(kind=dp), intent(in) :: r,rhon,rhop,fq,mu
        real(kind=dp) :: g !< Interaction strength (density dependence)
        integer :: k_n !< Number of grid points for integral in gap equation
        integer :: iters !< Number of iterations for self-consistence loop
        integer :: i !< Loop index
        real(kind=dp) :: del !< INM pairing gap
        real(kind=dp) :: del_old !< Previous iteration value of gap
        real(kind=dp) :: k !< Momentum
        real(kind=dp) :: dk !< Integration step
        real(kind=dp) :: qp_energy !< Quasiparticle energy
        real(kind=dp) :: integral !< Integral in gap equation
        real(kind=dp) :: hbms !< hbar2_2m_n * fq
        real(kind=dp) :: k_f !< Fermi momentum
        real(kind=dp) :: e_f !< Fermi energy
        
        if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .le. 15) then
            print '(a)', 'Pairing not implemented for BSk15 and earlier'
            stop
            
        ! BSk16 and onwards: Analytical expression for pairing gap in INM (PRC 82, 014313 (2010))
        else if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 16) then
            ! Fermi momentum
            k_f = (3._dp*pi**2 * rhon)**fr13
            
            select case(trim(force_string))
                ! Bare interaction (BSk16 uses effective mass)
                case("BSk16")
                    hbms = fq * hbar2_2m_q(0)
                    
                ! BSk17-29
                case default
                    hbms = hbar2_2m_q(0)
                    
            end select
            
            ! Interaction strength
            g = InterazChamelAnalyt(rhon,rhop,hbms,0)
            
            ! Fermi energy
            e_f = hbms * k_f**2
            
            ! Pairing gap
            del = 2._dp*e_f &
                & * exp(4._dp*pi**2*hbms / (k_f * InterazChamelAnalyt(rhon,rhop,hbms,0)) + fr12*LLambda(epsilon_Lambda/e_f))
            
        else
            ! Standard Skyrme forces
            g = pair_v0 * (1._dp - pair_eta * ((rhon+rhop)/rho0)**pair_alpha)
            
            ! Reset dk from value of pair_dk
            dk = pair_dk
            
            ! Number of integration points
            k_n = nint(pair_k_max/dk)
            
            ! Start guess for pairing gap
            del = pair_del_init
            
            ! Initialise del_old to large number
            del_old = 1.e9_dp
            
            ! Iterations
            iters = 0
            
            ! Self-consistence loop
            self_con: do while (abs(del-del_old) .gt. pair_tol)
                ! Count iterations
                iters = iters + 1
                
                ! If too many iterations have passed, convergence has not happened
                if (iters .gt. pair_max_iters) then
                    if (verbose .eqv. .true.) then
                        print '(/a)', 'WARNING: max iterations reached in solving gap equation, halving dk'
                        print '(a,1x,es24.16e3)', 'Radius =', r
                        print '(a,1x,es24.16e3)', 'Density =', rhon
                        print '(a,1x,es24.16e3)', 'Effective mass =', fq
                        print '(a,1x,es24.16e3)', 'Effective chemical potential =', mu
                        print '(a,1x,es24.16e3)', 'Interaction strength =', g
                        print '(a,1x,es24.16e3)', 'Current precision =', abs(del-del_old)
                    end if
                    
                    ! Halve dk before reattempting to solve equation
                    dk = dk / 2._dp
                    
                    if (verbose .eqv. .true.) then
                        print '(a,1x,es24.16e3)', 'New dk =',dk
                    end if
                    
                    ! Recalculate number of grid points
                    k_n = nint(pair_k_max/dk)
                    
                    ! Reset number of iterations
                    iters = 1
                    
                    ! New start guess for pairing gap (average of last 2 non-converged values)
                    del = fr12*(del+del_old)
                end if
                
                ! Update del_old to have del value from last iteration
                del_old = del
                
                ! Initialise integral sum to zero
                integral = 0._dp
                
                ! Integration over momenta
                !$OMP parallel do private(i,k,qp_energy) shared(k_n,dk,fq,mu,del,pair_qp_cut) &
                !$OMP & reduction(+:integral) num_threads(omp_get_max_threads())
                integration: do i=0, k_n
                    ! Momentum
                    k = real(i, kind=dp) * dk
                    
                    ! Quasiparticle energy
                        ! NB: Neutron mass hard-coded
                    qp_energy = sqrt((hbar2_2m_q(0)*fq*k**2 - mu)**2 + del**2)
                    
                    ! Only contribute to integral if less than cut-off energy
                    if (qp_energy .le. pair_qp_cut) then
                        ! Integral in gap equation
                        integral = integral + ( dk*k**2*del / qp_energy )
                    
                    else
                        if (smooth_qp_cut .eqv. .true.) then
                            ! Implement smooth cut-off with Gaussian form
                            integral = integral + ( dk*k**2*del / qp_energy ) * exp(-(qp_energy-pair_qp_cut)**2 / 100._dp)
                        end if
                        
                    end if
                end do integration
                !$OMP end parallel do
                
                ! Self-consistent step
                del = pair_mix * del_old + (1._dp - pair_mix) * (-g / (4._dp * pi**2)) * integral 
            end do self_con
            
            ! Update start guess for next call of function
    !         pair_del_init = del
            
        end if
            
        ! Result
        calc_gap = del
        
    end function calc_gap
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to perform BCS for protons
    subroutine bcs_protons()
        implicit none
        integer :: i !< Loop index
        real(kind=dp) :: mu_lower !< Lower bound on chemical potential
        real(kind=dp) :: mu_upper !< Upper bound on chemical potential
        real(kind=dp) :: num_p_bcs_old !< Previous iteration value of num_p_bcs
        integer :: iters !< Number of iterations for self-consistence loop
        
        
!         real(kind=dp) :: dens
!         dens=0.001_dp
!         do i=1,120
!             write(99,*) dens, 2._dp*sqrt(calc_e_f_q(dens,1._dp,0)*epsilon_Lambda) &
!                 & * exp(4._dp*pi**2*hbar2_2m_q(0)/((3._dp*pi**2*dens)**fr13 * InterazChamelAnalyt(dens,0._dp,hbar2_2m_q(0),0))) &
!                 & * exp(fr12*LLambda(epsilon_Lambda/calc_e_f_q(dens,1._dp,0))) / sqrt(epsilon_Lambda/calc_e_f_q(dens,1._dp,0))
!             dens = dens + 0.001_dp
!         end do
        
        
        !!! Calculate sp states !!!
        ! Call boundary routine
        call boundary(dr,Lmax,Ecut(1),strut_n,PS,NNST,JJP,LP,EP,1)
        
        ! Sort states
        call state_sort(JJP(1:NNST),LP(1:NNST),EP(1:NNST),PS(1:NNST,1:n))
        
        ! Normalise and square wavefunctions
        all_states_0: do i=1, NNST
            PS(i,1:n) = PS(i,1:n)**2 / (4._dp*pi*r(1:n)**2)
        end do all_states_0
        
        ! Initial bounds on guess for chemical potential
        mu_lower = minval(U_q(:,1))
        
        mu_upper = 0._dp
        
        ! Initialise num_p_bcs
        num_p_bcs = 0._dp
        
        ! Initialise num_p_bcs_old to large number
        num_p_bcs_old = 1.e9_dp
        
        ! Iterations
        iters = 0
        
        ! Use bisection to method to find solution to BCS number equation
        number: do while (abs(num_p_bcs-N_q(1)) .gt. pair_num_tol)
            ! Update num_p_bcs_old to have num_p_bcs value from last iteration
            num_p_bcs_old = num_p_bcs
            
            ! Count iterations
            iters = iters + 1
            
            ! Trial chemical potential in middle of upper and lower bounds
            mu_q(1) = fr12 * (mu_lower + mu_upper)
            
            ! New value of (BCS) number of protons
            num_p_bcs = calc_bcs_num_p(mu_q(1))
            
            ! Convergence reached, if many iterations and no change in num_p_bcs
            if ((iters .gt. pair_max_iters) .and. (abs(num_p_bcs-num_p_bcs_old) .lt. pair_num_tol)) then
                exit number
            end if
            
            ! Update bounds
            if (num_p_bcs .gt. N_q(1)) then
                ! Decrease upper bound if num_p_bcs too high
                mu_upper = mu_q(1)
            
            else
                ! Increase lower bound if num_p_bcs too low
                mu_lower = mu_q(1)
                
            end if
        end do number
        
        ! Make chemical potential 0 if not needed (for consistency)
!         if (calc_chem_pots .eqv. .false.) then
!             mu_q(1) = 0._dp
!         end if
        
    end subroutine bcs_protons
    
    
    !> @author M. Shelley, A. Pastore
    !> @brief
    !> Function to solve gap equation at given chemical potential "mu_p",
    !> returning the (BCS) number of protons
    !> @param[in]   mu_p    Proton chemical potential
    function calc_bcs_num_p(mu_p)
        implicit none
        real(kind=dp) :: calc_bcs_num_p
        real(kind=dp), intent(in) :: mu_p
        integer :: i !< Loop index
        real(kind=dp) :: e_pair_bcs_old !< Previous iteration value of e_pair_bcs
        real(kind=dp) :: xxx !< Factor for modifying U and V for BSk forces
        real(kind=dp) :: yyy !< Factor for determining when V is pinned to 0 for BSk forces
        
        ! Initial guess for gap for each sp state
        pair_gap_p(1:NNST) = -pair_del_init ! Using convention that gap is negative
        
        ! Initialise e_pair_bcs
        e_pair_bcs = 0._dp
        
        ! Initialise e_pair_bcs_old to large number
        e_pair_bcs_old = 1.e9_dp
        
        self_con: do while (abs(e_pair_bcs-e_pair_bcs_old) .gt. pair_tol)
!             ! Initialise pairing condensation energy
!             e_pair_q(1) = 0._dp
            
            ! Update e_pair_bcs_old to have e_pair_bcs value from last iteration
            e_pair_bcs_old = e_pair_bcs
            
            ! Quasiparticle energies
            e_qp_p(1:NNST) = sqrt((EP(1:NNST)-mu_p)**2 + pair_gap_p(1:NNST)**2)
            
            ! Occupation probabilities
            all_states_1: do i=1, NNST
                occ_v(i) = sqrt(fr12 * ( 1._dp - (EP(i)-mu_p) / e_qp_p(i) ))
                occ_u(i) = sqrt(fr12 * ( 1._dp + (EP(i)-mu_p) / e_qp_p(i) ))
            end do all_states_1
            
            ! Reset densities to 0
            rho_p_bcs(:) = 0._dp
            rho_anom_p_bcs(:) = 0._dp
            
            
            !!! U, V, DENSITIES, STRENGTH !!!
            ! BSk16 and onwards
            if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 16) then
                ! Loop over all states, sum up densities
                all_states_3: do i=1, NNST
                    ! Special modifications to V and U factors
                    xxx = (EP(i) - mu_p - epsilon_Lambda) / 0.25_dp
                    yyy = EP(i) - mu_p - epsilon_Lambda + 0.5_dp*log(0.005_dp)
                    
                    if (xxx .gt. 0._dp) then
                        occ_u(i) = occ_u(i) / sqrt(1._dp + exp(xxx))
                        occ_v(i) = occ_v(i) / sqrt(1._dp + exp(xxx))
                    end if
                    
                    if (yyy .gt. 0._dp) then
                        occ_v(i) = 0._dp
                    end if
                    
                    ! Proton density from wavefunctions
                    rho_p_bcs(1:n) = rho_p_bcs(1:n) &
                                        & + (JJP(i)+1) * occ_v(i)**2 * PS(i,1:n)
                    
                    ! Anomalous density
                    rho_anom_p_bcs(1:n) = rho_anom_p_bcs(1:n) &
                                        & + (JJP(i)+1) * occ_v(i) * occ_u(i) * PS(i,1:n)
                end do all_states_3
                
                ! Interaction strength
                radii: do i=1, n
                    strength_bcs(i) = InterazChamelAnalyt(rho_q(i,0),rho_p_bcs(i),hbar2_2m_q(1),1)
                end do radii
                
                ! Pairing renormalising
                if (mod(num_p,2) .eq. 0) then
                    strength_bcs(1:n) = strength_bcs(1:n) * fp_pos
                
                else
                    strength_bcs(1:n) = strength_bcs(1:n) * fp_neg
                
                end if
            
            ! Standard Skyrme forces, with DDDI pairing
            else
                ! Loop over all states, sum up densities
                all_states_2: do i=1, NNST
                    ! Quasiparticle cutoff
                    if (e_qp_p(i) .lt. pair_qp_cut) then
                        ! Proton density from wavefunctions
                        rho_p_bcs(1:n) = rho_p_bcs(1:n) &
                                            & + (JJP(i)+1) * occ_v(i)**2 * PS(i,1:n)
                        
                        ! Anomalous density
                        rho_anom_p_bcs(1:n) = rho_anom_p_bcs(1:n) &
                                            & + (JJP(i)+1) * occ_v(i) * occ_u(i) * PS(i,1:n)
                        
                    else
                        ! Smooth cutoff above pair_qp_cut
                        if (smooth_qp_cut .eqv. .true.) then
                            ! Proton density from wavefunctions
                            rho_p_bcs(1:n) = rho_p_bcs(1:n) &
                                    & + (JJP(i)+1) * occ_v(i)**2 * PS(i,1:n) * exp(-(e_qp_p(i)-pair_qp_cut)**2 / 100._dp)
                            
                            ! Anomalous density
                            rho_anom_p_bcs(1:n) = rho_anom_p_bcs(1:n) &
                                    & + (JJP(i)+1) * occ_v(i) * occ_u(i) * PS(i,1:n) * exp(-(e_qp_p(i)-pair_qp_cut)**2 / 100._dp)
                        end if
                        
                    end if
                end do all_states_2
                
                ! Interaction strength (from 0 to 1, has density dependence) * strength
                strength_bcs(1:n) = pair_v0 * (1._dp - pair_eta * (((rho_q(1:n,0) + rho_p_bcs(1:n)) / rho0)**pair_alpha))
                
            end if
            
            
            !!! FIELD, GAPS, ENERGY !!!
            ! Pairing field
            delta_p(1:n) = fr12 * strength_bcs(1:n) * rho_anom_p_bcs(1:n)
            
            ! Pairing gaps
            all_states_4: do i=1, NNST
                pair_gap_p(i) = WS_integral(delta_p(1:n) * PS(i,1:n), n)
!                 pair_gap_p(i) = pair_mix * pair_gap_p(i) + (1._dp-pair_mix) * WS_integral(delta_p(1:n) * PS(i,1:n), n)
            end do all_states_4
            
            ! Number of protons
            num_p_bcs = sum((JJP(1:NNST) + 1) * (occ_v(1:NNST)**2))
            
            ! BCS pairing energy
            e_pair_bcs = fr12 * WS_integral(delta_p(1:n) * rho_anom_p_bcs(1:n), n)
            
!             ! Pairing condensation energy
!             all_states_5: do i=1, NNST
!                 ! Quasiparticle cutoff
!                 if (e_qp_p(i) .lt. pair_qp_cut) then
!                     e_pair_q(1) = e_pair_q(1) - fr14 * ((JJP(i)+1)*pair_gap_p(i)**2 / e_qp_p(i))
!                 end if
!             end do all_states_5
        end do self_con
        
        ! Return converged value of num_p_bcs
        calc_bcs_num_p = num_p_bcs
        
    end function calc_bcs_num_p
    
    
    
    !> @author A. Pastore, M. Shelley
    !> @brief
    !> Function to calculate the interaction strength for the effective contact
    !> pairing force. Modified by M. Shelley.
    !> @param[in]   rhon    Neutron density
    !> @param[in]   rhop    Proton density
    !> @param[in]   hbm     \f$\frac{\hbar^2}{2m^*}\f$
    !> @param[in]   itz     Isospin
    function InterazChamelAnalyt(rhon,rhop,hbm,itz)
        implicit none
        real(kind=dp) :: InterazChamelAnalyt
        real(kind=dp),  intent(in) :: rhon,rhop,hbm
        integer,        intent(in) :: itz
        real(kind=dp) :: ef !< Fermi energy
        real(kind=dp) :: YY !< Neutron-proton composition
        real(kind=dp) :: xk0 !< "Average" Fermi momentum
        real(kind=dp) :: integral !< Integral needed for interaction
        real(kind=dp) :: kf !< Fermi momentum
        
        ! Fermi momentum for neutrons or protons
        if (itz .eq. 0) then
            kf = (3._dp * pi**2 * abs(rhon))**fr13
            
        else
            kf = (3._dp * pi**2 * abs(rhop))**fr13
            
        endif
        
        ! Fermi energy
        ef = kf**2 * hbm
        
        ! Neutron-proton composition
        YY = (rhon - rhop) / (rhon + rhop)
        
        ! "Average" Fermi momentum
        xk0 = (fr32 * pi**2 * abs(rhon + rhop))**fr13
        
        if (kf .lt. 1.e-10_dp) then
            ! Low density limit of integral
            integral = 2._dp * sqrt(epsilon_Lambda)
        
        else
            ! Weak-coupling (all terms) analytical expression for integral
            integral = sqrt(ef) * (2._dp * log(2._dp * ef / Delta_parametric(kf,YY,itz,xk0)) + LLambda(epsilon_Lambda / ef))
        
        end if
        
        ! Final expression for interaction
        InterazChamelAnalyt = -8._dp * pi**2 * (hbm**fr32) / integral
        
    end function InterazChamelAnalyt
    
    
    !> @author A. Pastore
    !> @brief
    !> Function to calculate the pairing cutoff for the effective contact
    !> pairing force
    !> @param[in]   x   (DDCI cutoff) / (Fermi energy)
    function LLambda(x)
        implicit none
        real(kind=dp) :: LLambda
        real(kind=dp), intent(in) :: x
        
        LLambda = log(16._dp * x) + 2._dp * sqrt(1._dp + x) - 2._dp * log(1._dp + sqrt(1._dp + x)) - 4._dp
        
    end function LLambda
    
    
    !> @author A. Pastore, M. Shelley
    !> @brief
    !> Function to calculate the pairing gap using the analytical BSk
    !> expression. Modified by M. Shelley.
    !> @param[in]   kf  Fermi momentum
    !> @param[in]   YY  Neutron-proton composition
    !> @param[in]   itz Isospin
    !> @param[in]   xk0 "Average" fermi momentum
    function Delta_parametric(kf,YY,itz,xk0)
        implicit none
        real(kind=dp) :: Delta_parametric
        real(kind=dp),  intent(in) :: kf,xk0,YY
        integer,        intent(in) :: itz
        real(kind=dp) :: d0,k1,k2,k3,b0,b1,b2,b3 !< Interaction parameters
        real(kind=dp) :: rho0,rn,rp !< Total, neutron and proton densities
        real(kind=dp) :: aux,dn,ds !< Auxiliary quantities for calculations
        
        rho0 = xk0**3 / (fr32 * pi**2)
        rp = rho0 * (1._dp - YY) / 2._dp
        rn = rho0 - rp
        
        select case(trim(force_string))
            ! Bare interaction (BSk16 uses effective mass)
            case("BSk16")
                d0 = sqrt(910.603_dp)
                k1 = 1.38297_dp
                k2 = 1.57068_dp
                k3 = 0.905237_dp
                
                if (kf .lt. 1.57_dp)then
                    aux = (kf - k2)**2 / ((kf - k2)**2  +  k3**2)
                    dn = d0 * kf**2 / (kf**2  +  k1**2) * aux
                    
                else
                    dn = 0._dp
                    
                endif
                
                Delta_parametric = dn
                
            ! Cao et al. (free spectrum and screening)
                ! Forces BSk17+ use effective mass of 1 in pairing interaction
            case default
                d0 = 14.9003_dp
                k1 = 1.18847_dp
                k2 = 1.51854_dp
                k3 = 0.639489_dp
                
                b0 = 204.396_dp
                b1 = 0.943146_dp
                b2 = 1.52786_dp
                b3 = 2.11577_dp
                
                if (kf .lt. 1.528_dp)then
                    aux = (kf - k2)**2 / ((kf - k2)**2 + k3**2)
                    dn = d0 * kf**2 / (kf**2 + k1**2) * aux
                    
                else
                    dn = 0._dp
                    
                endif
                
                if (xk0 .lt. 1.5184_dp)then
                    aux = (xk0 - b2)**2 / ((xk0**2 + b1**2) * ((xk0 - b2)**2 + b3**2))
                    ds = b0 * xk0**3 * aux / b2
                    ds = (1._dp - abs(YY)) * ds
                    
                    if (ds .lt. 0) then
                        ds = 0._dp
                    end if
                    
                else
                    ds = 0._dp
                    
                endif
                
                if (itz .eq. 0) then
                    Delta_parametric = ds + YY * (YY + 1._dp) / 2._dp * dn
                end if
                
                if (itz .eq. 1) then
                    Delta_parametric = ds + YY * (YY - 1._dp) / 2._dp * dn
                end if
                
                if (Delta_parametric .lt. 0._dp) then
                    Delta_parametric = 0._dp
                end if
                
        end select
        
    end function Delta_parametric


end module pairing