!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Module to hold routines for carrying out Strutinsky correction
module strutinsky
    
    use parameters
    use routines
    use rspace
    implicit none
    
    
    contains
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to allocate + initialise variables needed in 'boundary'
    !> routine, and allocate arrays needed for BCS
    subroutine sp_states_setup()
        implicit none
        integer :: allo_stat !< Error status for array allocation
        integer :: q !< isospin
        
        ! Allocate fields
        allocate(dhmen(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating dhmen'
        allocate(d2hmen(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d2hmen'
        allocate(hb2m(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating hb2m'
        allocate(vpot(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating vpot'
        allocate(vso(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating vso'
                
        ! Allocate output arrays
        allocate(PS(1:Nmaxstate,1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating PS'
        allocate(JJP(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating JJP'
        allocate(LP(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating LP'
        allocate(EP(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating EP'
        
        
        !!! Initialise field arrays !!!
        ! Central fields
        vpot(1:n,0:1) = U_q(1:n,0:1)
        
        isospin: do q=0, 1
            ! Change effective mass ratios 'f_q' into \f$\frac{\hbar^2}{2m^*}\f$
            dhmen(1:n,q) = del_f_q(1:n,q) * hbar2_2m_q(q)
            d2hmen(1:n,q) = del2_f_q(1:n,q) * hbar2_2m_q(q)
            hb2m(1:n,q) = f_q(1:n,q) * hbar2_2m_q(q)
            
            ! Spin-orbit fields
                ! Multiply by 2/r to be compatible with routines in 'rspace'
            vso(1:n,q) = W_q(1:n,q) * 2._dp / r(1:n)
        end do isospin
        
        
        !!! Arrays for BCS for protons !!!
        allocate(pair_gap_p(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating pair_gap_p'
        allocate(delta_p(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating delta_p'
        allocate(rho_p_bcs(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating rho_p_bcs'
        allocate(rho_anom_p_bcs(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating rho_anom_p_bcs'
        allocate(strength_bcs(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating strength_bcs'
        allocate(occ_v(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating occ_v'
        allocate(occ_u(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating occ_u'
        allocate(e_qp_p(1:Nmaxstate),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating e_qp_p'
        
        
    end subroutine sp_states_setup
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to deallocate arrays needed in 'boundary' routine, and
    !> deallocate arrays needed for BCS
    subroutine sp_states_deallocate()
        implicit none
        integer :: deallo_stat !< Error status for array deallocation
        
        ! Deallocate field arrays
        deallocate(dhmen,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating dhmen'
        deallocate(d2hmen,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d2hmen'
        deallocate(hb2m,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating hb2m'
        deallocate(vpot,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating vpot'
        deallocate(vso,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating vso'
        
        ! Deallocate output arrays
        deallocate(PS,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating PS'
        deallocate(JJP,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating JJP'
        deallocate(LP,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating LP'
        deallocate(EP,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating EP'
        
        
        !!! Arrays for BCS for protons !!!
        deallocate(pair_gap_p,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating pair_gap_p'
        deallocate(delta_p,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating delta_p'
        deallocate(rho_p_bcs,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating rho_p_bcs'
        deallocate(rho_anom_p_bcs,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating rho_anom_p_bcs'
        deallocate(strength_bcs,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating strength_bcs'
        deallocate(occ_v,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating occ_v'
        deallocate(occ_u,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating occ_u'
        deallocate(e_qp_p,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating e_qp_p'
        
    end subroutine sp_states_deallocate
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to sort single particle states by their energies. Use bubble
    !> sort.
    !> @param[in,out]   sp_J    Array for (2)J's of single particle states
    !> @param[in,out]   sp_L    Array for L's of single particle states
    !> @param[in,out]   sp_E    Array for energies of single particle states
    !> @param[in,out]   sp_wfns Array for wavefunctions of single particle
    !>                          states
    subroutine state_sort(sp_J,sp_L,sp_E,sp_wfns)
        implicit none
        integer,        dimension(:),   intent(inout)   :: sp_J,sp_L
        real(kind=dp),  dimension(:),   intent(inout)   :: sp_E
        real(kind=dp),  dimension(:,:), intent(inout)   :: sp_wfns
        integer :: i !< Loop index
        logical :: swapped !< For determining if more swaps are needed
        integer :: int_swap !< Store for J and L values during swaps
        real(kind=dp) :: real_swap !< Store for E values during swaps
        real(kind=dp),  dimension(1:n) :: real_array_swap !< Store for wfns values during swaps
        
        swapped = .true.
        
        ! Keep doing passes over states until no swaps happened in previous
        ! pass
        pass: do while (swapped .eqv. .true.)
            swapped = .false.
            
            ! Do a pass over all consecutive pairs of states
            states: do i=2, NNST
                ! Compare pair of states: swap pair if energies are out of
                ! order
                if (sp_E(i-1) .gt. sp_E(i)) then
                    ! J
                    int_swap = sp_J(i)
                    sp_J(i) = sp_J(i-1)
                    sp_J(i-1) = int_swap
                    ! L
                    int_swap = sp_L(i)
                    sp_L(i) = sp_L(i-1)
                    sp_L(i-1) = int_swap
                    ! E
                    real_swap = sp_E(i)
                    sp_E(i) = sp_E(i-1)
                    sp_E(i-1) = real_swap
                    ! Wfns
                    real_array_swap(:) = sp_wfns(i,:)
                    sp_wfns(i,:) = sp_wfns(i-1,:)
                    sp_wfns(i-1,:) = real_array_swap(:)
                    
                    ! Remember that swap has taken place, so that sorting
                    ! continues
                    swapped = .true.
                end if
            end do states
        end do pass
        
    end subroutine state_sort
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the Strutinsky shell correction energy for
    !> isospin \f$q\f$. First, sum over occupied states of s.p. energies.
    !> Second, integrate over W-S cell of smoothed ETF densities and fields.
    !> @param[in] q Isospin
    subroutine calc_e_sc(q)
        implicit none
        integer, intent(in) :: q !< Isospin
        integer :: tot_parts !< Number of particles to count up to
        integer :: state !< Index for state
        integer :: num_parts_counted !< Counter for number of particles
        integer :: occupation !< \f$2J+1\f$ for each state
        real(kind=dp) :: e_sp_sum !< Sum of energies of single particle states
        real(kind=dp) :: smoothed_integral !< Integral over smoothed ETF densities and fields
!         real(kind=dp) :: test_1,test_2,test_3
!         integer :: ierr !< I/O error storage
        
        !!! Sum over occupied s.p. states !!!
        ! Call boundary routine
        call boundary(dr,Lmax,Ecut(q),strut_n,PS,NNST,JJP,LP,EP,q)
        
        ! Sort states
        call state_sort(JJP(1:NNST),LP(1:NNST),EP(1:NNST),PS(1:NNST,:))
        
        ! Number of neutrons or protons for occupying states
        tot_parts = nint(N_q(q))
        
        ! Counter for states
        state = 0
        
        ! Counter for particles
        num_parts_counted = 0
        
        ! For summing energies of single particle states
        e_sp_sum = 0._dp
        
        ! Loop over calculated single-particle states
        states: do while (num_parts_counted .lt. tot_parts)
            ! Increment for storage of next state
            state = state + 1
            
            occupation = JJP(state) + 1 ! 2J+1
            
            if (num_parts_counted + occupation .le. tot_parts) then
                ! If not all particles counted, add (2J+1) times energy of state
                e_sp_sum = e_sp_sum + EP(state) * occupation
                
                ! Add (2J+1) particles to counter
                num_parts_counted = num_parts_counted + occupation
!                 mu_q(1) = fr12*(EP(state)+EP(state+1))
!                 mu_q(1) = EP(state+1)
                
            else
                ! Add energy of state >= 1 time so that all particles are
                ! accounted for
                e_sp_sum = e_sp_sum + EP(state) * (tot_parts-num_parts_counted)
                
                ! All particles counted
                num_parts_counted = tot_parts
!                 mu_q(1) = EP(state)
                
            end if
        end do states
        
        !!! Strutinsky correction calculation
        ! Integral of smoothed quantities
        smoothed_integral = WS_integral( hb2m(1:n,q)*(tau_TF_q(1:n,q)+tau_2_L_q(1:n,q)+tau_2_NL_q(1:n,q))&!tau_ETF_q(1:n,q) &
                                            & + rho_q(1:n,q)*U_q(1:n,q) &
                                            & + J_q(1:n,q)*W_q(1:n,q) , n )
        
        
!                                             & + J_2_q(1:n,q)*W_q(1:n,q) , n )
!         smoothed_integral = smoothed_integral + e_kinetic_t
        
        
!                 ! Details of correction calculation
!                 if (q .eq. 1) then
!                     ! Components of integral
!                     test_1 = WS_integral(hb2m(1:n,q)*tau_ETF_q(1:n,q),n)
!                     test_2 = WS_integral(rho_q(1:n,q)*U_q(1:n,q),n)
!                     test_3 = WS_integral(J_q(1:n,q)*W_q(1:n,q),n)
!                     
!                     open(unit=300,file='strut_correction_details.dat',status='unknown',iostat=ierr)
!                     if (ierr.ne.0) stop
!                     
!                     ! Write details to file
!                     write(300,n_floats(5)) e_sp_sum,smoothed_integral,test_1,test_2,test_3
!                     
!                     close(unit=300,iostat=ierr,status='keep')
!                     if (ierr.ne.0) stop
!                 end if
        
        ! Strutinsky correction is (sum of s.p. energies) - (integral of smoothed quantities)
        e_sc_q(q) = e_sp_sum - smoothed_integral
        
    end subroutine calc_e_sc
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the Strutinsky shell correction energy for
    !> protons, using BCS occupation probabilities
    subroutine calc_e_sc_pair()
        implicit none
        real(kind=dp) :: e_sp_sum !< Sum of energies of single particle states
        integer :: i !< Loop index
        real(kind=dp) :: smoothed_integral !< Integral over smoothed ETF densities and fields
!         real(kind=dp) :: test_1,test_2,test_3
!         integer :: ierr !< I/O error storage
        
        ! For summing energies of single particle states
        e_sp_sum = 0._dp
        
        ! Loop over all states
        all_states: do i=1, NNST
            ! Multiple each energy by (2J+1) and by the particle occupation probability (V^2)
            e_sp_sum = e_sp_sum + occ_v(i)**2 * (JJP(i)+1) * EP(i)
        end do all_states
        
        ! Integral of smoothed quantities
        smoothed_integral = WS_integral( hb2m(1:n,1)*tau_ETF_q(1:n,1) &
                                            & + rho_q(1:n,1)*U_q(1:n,1) &
                                            & + J_q(1:n,1)*W_q(1:n,1) , n )
        
!         ! Details of correction calculation
!         ! Components of integral
!         test_1 = WS_integral(hb2m(1:n,1)*tau_ETF_q(1:n,1),n)
!         test_2 = WS_integral(rho_q(1:n,1)*U_q(1:n,1),n)
!         test_3 = WS_integral(J_q(1:n,1)*W_q(1:n,1),n)
!         
!         open(unit=300,file='strut_correction_details.dat',status='unknown',iostat=ierr)
!         if (ierr.ne.0) stop
!         
!         ! Write details to file
!         write(300,n_floats(5)) e_sp_sum,smoothed_integral,test_1,test_2,test_3
!         
!         close(unit=300,iostat=ierr,status='keep')
!         if (ierr.ne.0) stop
        
        ! Strutinsky correction is (sum of s.p. energies) - (integral of smoothed quantities)
        e_sc_q(1) = e_sp_sum - smoothed_integral
        
    end subroutine calc_e_sc_pair
    
    
end module strutinsky