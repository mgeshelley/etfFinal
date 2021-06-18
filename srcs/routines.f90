!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Module to hold main routines: densities, energies
module routines

    use parameters
    implicit none

    integer :: file_unit_0, file_unit_1 !< Unit numbers for opening files
    
    
    contains
    
    
    function str2int(string)
        implicit none
        integer :: str2int
        character(len=*), intent(in) :: string
        
        read(string,*) str2int
        
    end function str2int
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to allocate arrays for all densities and effective masses,
    !> and for their derivatives
    subroutine allocate_arrays()
        implicit none
        integer :: allo_stat !< Error status for array allocation
        
        ! Allocate densities and density derivatives
            ! Add 3 points for extrapolations for 5- or 7-point derivatives
        allocate(rho_q(-2:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating rho_q'
        allocate(del_rho_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del_rho_q'
        allocate(del2_rho_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del2_rho_q'
        allocate(rho_t(-2:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating rho_t'
        allocate(del_rho_t(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del_rho_t'
        allocate(del2_rho_t(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del2_rho_t'
        
        ! Allocate effective mass arrays
            ! Add 3 points for extrapolations for 5- or 7-point derivatives
        allocate(f_q(-2:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating f_q'
        allocate(del_f_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del_f_q'
        allocate(del2_f_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del2_f_q'
        
        ! Allocate arrays for spin-orbit fields (and their divergences) and
        ! central nuclear potentials
        allocate(W_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating W_q'
        allocate(div_W_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating div_W_q'
        allocate(U_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating U_q'
        
        ! Allocate arrays for spin current densities and their derivatives
        allocate(J_2_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating J_2_q'
        allocate(J_q(-2:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating J_q'
        allocate(del_J_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del_J_q'
        allocate(div_J_2_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating div_J_2_q'
        allocate(div_J_4_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating div_J_4_q'
        allocate(div_J_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating div_J_q'
        allocate(J_t(-2:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating J_t'
        allocate(del_J_t(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating del_J_t'
        allocate(div_J_t(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating div_J_t'
        
        ! Allocate arrays for kinetic energy density and its contributions
        allocate(tau_TF_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_TF_q'
        allocate(tau_2_L_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_2_L_q'
        allocate(tau_2_NL_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_2_NL_q'
        allocate(tau_ETF_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_ETF_q'
        allocate(tau_ETF_t(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_ETF_t'
        allocate(tau_2_cont_q(1:n,1:4,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_2_cont_q'
        
        ! Allocate arrays for densities, fields, and derivatives, for fourth
        ! order corrections
        allocate(d2_rho_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d2_rho_q'
        allocate(d3_rho_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d3_rho_q'
        allocate(d4_rho_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d4_rho_q'
        allocate(d2_f_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d2_f_q'
        allocate(d3_f_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d3_f_q'
        allocate(d4_f_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d4_f_q'
        allocate(A_q(-2:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating A_q'
        allocate(d1_A_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d1_A_q'
        allocate(d2_A_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d2_A_q'
        allocate(d3_A_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d3_A_q'
        allocate(d4_A_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating d4_A_q'
        allocate(tau_4_so_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_4_so_q'
        allocate(tau_4_no_spin_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_4_no_spin_q'
        allocate(J_4_q(1:n,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating J_4_q'
        allocate(tau_4_cont_q(1:n,1:4,0:1),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating tau_4_cont_q'
        
        ! Allocate arrays for Skyrme and Coulomb energy densities, Coulomb
        ! potentials, and charge density
        allocate(e_density_sky(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating e_density_sky'
        allocate(e_density_field(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating e_density_field'
        allocate(v_c_di(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating v_c_di'
        allocate(v_c_ex(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating v_c_ex'
        allocate(e_density_c_di(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating e_density_c_di'
        allocate(e_density_c_ex(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating e_density_c_ex'
        allocate(rho_ch(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating rho_ch'
        
        
        ! Allocate arrays for neutron pairing gap and pairing field
        allocate(pair_gap_n(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating pair_gap_n'
        allocate(delta_n(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating delta_n'
        
        ! Allocate arrays for electron contributions
        allocate(v_c_pe(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating v_c_pe'
        
    end subroutine allocate_arrays
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate number of mesh points for a uniform grid,
    !> allocate array for r, and then populate r with the mesh points
    subroutine initialise_uniform_mesh()
        implicit none
        integer :: allo_stat !< Error status for array allocation
        integer :: i !< Loop index
        
        ! If reading r_ws from profile, set r_max = r_ws
        if (profile_r_max .eqv. .true.) then
            r_max = r_ws
        end if
        
        ! If specify_n is true, set dr according to r_max and n
        if (specify_n .eqv. .true.) then
            dr = r_max/real(n,kind=dp)
            
        ! If false, set n according to r_max and dr
        else
            n = nint(r_max/dr)
            
            ! dr probably does not divide perfectly into r_max, so adjust it
            ! using calculated value of n
            dr = r_max/real(n,kind=dp)
            
        end if
        
        ! Allocate r
        allocate(r(1:n),stat=allo_stat)
        if (allo_stat .ne. 0) stop 'Error in allocating r'
        
        ! Populate r with mesh points
        mesh: do i=1, n
            r(i) = dr * real(i,kind=dp)
        end do mesh
        
        ! If reading r_ws from profile, and specified strut_r_max > r_ws, set
        ! strut_r_max = r_ws
        if (profile_r_max .eqv. .true. .and. strut_r_max .gt. r_ws) then
            strut_r_max = r_ws
        end if
        
        ! Number of mesh points for grid for "boundary" for Strutinsky Integral
        ! correction
        strut_n = nint(strut_r_max/dr)
        
    end subroutine initialise_uniform_mesh
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to deallocate array for r
    subroutine deallocate_mesh()
        implicit none
        integer :: deallo_stat !< Error status for array deallocation
        
        ! Deallocate r
        deallocate(r,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating r'
        
    end subroutine deallocate_mesh
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to deallocate arrays for all densities and effective masses,
    !> and for their derivatives
    subroutine deallocate_arrays()
        implicit none
        integer :: deallo_stat !< Error status for array deallocation
        
        ! Deallocate densities and derivatives
        deallocate(rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating rho_q'
        deallocate(del_rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del_rho_q'
        deallocate(del2_rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del2_rho_q'
        deallocate(rho_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating rho_t'
        deallocate(del_rho_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del_rho_t'
        deallocate(del2_rho_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del2_rho_t'
        
        ! Deallocate effective mass arrays
        deallocate(f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating f_q'
        deallocate(del_f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del_f_q'
        deallocate(del2_f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del2_f_q'
        
        ! Allocate arrays for spin-orbit fields (and their divergences) and
        ! central nuclear potentials
        deallocate(W_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating W_q'
        deallocate(div_W_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating div_W_q'
        deallocate(U_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating U_q'
        
        ! Allocate arrays for spin current densities and their derivatives
        deallocate(J_2_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating J_2_q'
        deallocate(J_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating J_q'
        deallocate(del_J_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del_J_q'
        deallocate(div_J_2_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating div_J_2_q'
        deallocate(div_J_4_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating div_J_4_q'
        deallocate(div_J_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating div_J_q'
        deallocate(J_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating J_t'
        deallocate(del_J_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating del_J_t'
        deallocate(div_J_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating div_J_t'
        
        ! Deallocate arrays for kinetic energy density and its contributions
        deallocate(tau_TF_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_TF_q'
        deallocate(tau_2_L_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_2_L_q'
        deallocate(tau_2_NL_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_2_NL_q'
        deallocate(tau_ETF_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_ETF_q'
        deallocate(tau_ETF_t,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_ETF_t'
        deallocate(tau_2_cont_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_2_cont_q'
        
        ! Deallocate arrays for densities, fields, and derivatives, for fourth
        ! order corrections
        deallocate(d2_rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d2_rho_q'
        deallocate(d3_rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d3_rho_q'
        deallocate(d4_rho_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d4_rho_q'
        deallocate(d2_f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d2_f_q'
        deallocate(d3_f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d3_f_q'
        deallocate(d4_f_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d4_f_q'
        deallocate(A_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating A_q'
        deallocate(d1_A_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d1_A_q'
        deallocate(d2_A_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d2_A_q'
        deallocate(d3_A_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d3_A_q'
        deallocate(d4_A_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating d4_A_q'
        deallocate(tau_4_so_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_4_so_q'
        deallocate(tau_4_no_spin_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_4_no_spin_q'
        deallocate(J_4_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating J_4_q'
        deallocate(tau_4_cont_q,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating tau_4_cont_q'
        
        ! Deallocate arrays for Skyrme and Coulomb energy densities, Coulomb
        ! potentials, and charge density
        deallocate(e_density_sky,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating e_density_sky'
        deallocate(e_density_field,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating e_density_field'
        deallocate(v_c_di,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating v_c_di'
        deallocate(v_c_ex,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating v_c_ex'
        deallocate(e_density_c_di,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating e_density_c_di'
        deallocate(e_density_c_ex,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating e_density_c_ex'
        deallocate(rho_ch,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating rho_ch'
        
        ! Deallocate arrays for neutron pairing gap and pairing field
        deallocate(pair_gap_n,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating pair_gap_n'
        deallocate(delta_n,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating delta_n'
        
        ! Deallocate arrays for electron contributions
        deallocate(v_c_pe,stat=deallo_stat)
        if (deallo_stat .ne. 0) stop 'Error in deallocating v_c_pe'
        
    end subroutine deallocate_arrays
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the neutron or proton density at a given radius r,
    !> using supplied input parameters
    !> @param[in] rhoqgas   Asymptotic density in gas far from surface
    !> @param[in] rhoqliq   Asymptotic density in cluster far from surface
    !> @param[in] rq        Cluster radius
    !> @param[in] aq        Surface diffuseness
    !> @param[in] gammaq    Parameter allowing for "asymmetric surface"
    !> @param[in] r         Radius
    function calc_rho_q(rhoqgas,rhoqliq,rq,aq,gammaq,r)
        implicit none
        real(kind=dp) :: calc_rho_q
        real(kind=dp), intent(in) :: rhoqgas,rhoqliq,rq,aq,gammaq,r
        
        ! Modified Fermi-Dirac profile
        calc_rho_q = rhoqgas + ( rhoqliq - rhoqgas ) / ( 1._dp + exp((r-rq)/aq) )**gammaq
!         if (((rq-r_max)/(r-r_max))**2 .gt. 100._dp) then
!             calc_rho_q = rhoqgas
!             
!         else
!             calc_rho_q = rhoqgas + ( rhoqliq - rhoqgas ) / ( 1._dp + exp(((rq-r_max)/(r-r_max))**2 - 1._dp) * exp((r-rq)/aq) )
!             
!         end if
        
        ! Peg all densities less than 10^{-80} (i.e. protons far from cluster)
        ! to 10^{-80}
        if (calc_rho_q .lt. 1.e-80_dp) then
            calc_rho_q = 1.e-80_dp
        end if
        
    end function calc_rho_q
    
    
!     !> @author M. Shelley
!     !> @brief
!     !> Subroutine to calculate analytical derivatives of a Fermi-Dirac density
!     !> profile
!     !> @param[in]       order   Order of derivative
!     !> @param[in]       pars    Profile parameters
!     !> @param[in,out]   deriv   Derivative array
!     subroutine fd_profile_derivative(order,pars,deriv)
!         implicit none
!         integer,                        intent(in)      :: order
!         real(kind=dp),  dimension(1:5), intent(in)      :: pars
!         real(kind=dp),  dimension(1:n), intent(inout)   :: deriv
!         real(kind=dp),  dimension(1:n)  :: numer_fac !< Factor appearing in numerators
!         real(kind=dp),  dimension(1:n)  :: denom_fac !< Factor appearing in denominators
!         
!         numer_fac = (pars(2)-pars(1)) * exp((r(:)+pars(3))/pars(4))
!         denom_fac = exp(pars(3)/pars(4)) + exp(r(:)/pars(4))
!         
!         select case(order)
!             ! First order
!             case(1)
!                 deriv = -1._dp * numer_fac / (pars(4)*denom_fac**2)
!                 
!             case default
!                 print '(a)', 'STOP: no order selected'
!                 stop
!         end select
!         
!     end subroutine fd_profile_derivative
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the effective mass ratio
    !> \f$f_q=\frac{m}{m^*_q}\f$ for neutrons or protons, for standard Skyrme
    !> or for BSk forces
    !> @param[in]       rhot    Total density
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       q       Isospin
    !> @param[in,out]   fq      Neutron or proton effective mass
    subroutine calc_f_q(rhot,rhoq,q,fq)
        implicit none
        real(kind=dp),  dimension(:),   intent(in)      :: rhot,rhoq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(:),   intent(inout)   :: fq
        
        
        ! BSk forces with extra terms
        if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 18) then
            fq = 1._dp + fr14/hbar2_2m_q(q) * &
                & ( t1 * ((1._dp+fr12*x1)*rhot - (fr12+x1)*rhoq) &
                &   + t2 * ((1._dp+fr12*x2)*rhot + (fr12+x2)*rhoq) &
                &   + t4 * ((1._dp+fr12*x4)*rhot - (fr12+x4)*rhoq) * rhot**beta &
                &   + t5 * ((1._dp+fr12*x5)*rhot + (fr12+x5)*rhoq) * rhot**gamma )
            
        ! Standard Skyrme forces
        else
            fq = 1._dp + (B3*rhot + B4*rhoq) / hbar2_2m_q(q)
            
        end if
        
    end subroutine calc_f_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the spin-orbit field \f$\textbf{W}_q\f$ for
    !> neutrons or protons
    !> @param[in]       delrho  Gradient of total density
    !> @param[in]       delrhoq Gradient of neutron or proton density
    !> @param[in,out]   Wq      Array for spin-orbit field
    subroutine calc_W_q(delrho,delrhoq,Wq)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: delrho,delrhoq
        real(kind=dp),  dimension(1:n), intent(inout)   :: Wq
        
        ! All arrays are 1:n
        Wq(:) = -B9 * (delrho+delrhoq)
        
    end subroutine calc_W_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the zeroth-order contribution to the kinetic
    !> energy density for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in,out]   tauTF   Array for zeroth-order contributions
    subroutine calc_tau_TF(rhoq,tauTF)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq
        real(kind=dp),  dimension(1:n), intent(inout)   :: tauTF
        
        ! All arrays are 1:n
        tauTF(:) = fr35 * (3._dp*pi**2)**fr23 * rhoq**fr53
        
    end subroutine calc_tau_TF
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the local second-order contribution to the kinetic
    !> energy density for neutrons or protons
    !> @param[in]       rhoq        Neutron or proton density
    !> @param[in]       delrhoq     Gradient of proton or neutron density
    !> @param[in]       del2rhoq    Laplacian of proton or neutron density
    !> @param[in,out]   tau2Lq      Array for local second-order contributions
    !> @param[in,out]   t2contq     Array for individual second-order
    !>                              contributions
    subroutine calc_tau_2_L(rhoq,delrhoq,del2rhoq,tau2Lq,t2contq)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,delrhoq,del2rhoq
        real(kind=dp),  dimension(1:n), intent(inout)   :: tau2Lq,t2contq
        
        ! All arrays are 1:n
        tau2Lq = fr1_36*delrhoq**2/rhoq + fr13*del2rhoq
        
        ! Store second-order contribution (gradients of density)
        t2contq = tau2Lq
        
    end subroutine calc_tau_2_L
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the non-local second-order contribution to the
    !> kinetic energy density for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       delrhoq Gradient of neutron or proton density
    !> @param[in]       fq      Effective mass for neutrons or protons
    !> @param[in]       delfq   Gradient of effective mass for neutrons or
    !>                          protons
    !> @param[in]       del2fq  Laplacian of effective mass for neutrons or
    !>                          protons
    !> @param[in]       Wq      Spin-orbit field
    !> @param[in]       q       Isospin
    !> @param[in,out]   tau2NLq Array for non-local second-order contributions
    !> @param[in,out]   t2contq Array for individual second-order contributions
    subroutine calc_tau_2_NL(rhoq,delrhoq,fq,delfq,del2fq,Wq,q,tau2NLq,t2contq)
        implicit none
        real(kind=dp),  dimension(1:n),     intent(in)      :: rhoq,delrhoq,fq,delfq,del2fq,Wq
        integer,                            intent(in)      :: q
        real(kind=dp),  dimension(1:n),     intent(inout)   :: tau2NLq
        real(kind=dp),  dimension(1:n,2:4), intent(inout)   :: t2contq
        
        ! All arrays are 1:n
        
        ! Terms containing gradients of effective mass
        t2contq(:,2) = fr16 * rhoq * (del2fq/fq) - fr1_12 * rhoq * (delfq/fq)**2
        
        ! Terms containing gradients of density and of effective mass
        t2contq(:,3) = fr16 * (delrhoq*delfq) / fq
        
        ! Terms containing spin-orbit contirbutions
        t2contq(:,4) = fr12/hbar2_2m_q(q)**2 * rhoq * (Wq/fq)**2
        
        ! Store total second-order contribution (non-local)
        tau2NLq(:) = t2contq(:,2) + t2contq(:,3) + t2contq(:,4)
        
    end subroutine calc_tau_2_NL
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the fourth-order contributions to the kinetic
    !> energy density (without spin-orbit contributions) \f$\textbf{tau}_q\f$
    !> for neutrons or protons
    !> @param[in]       rhoq        Neutron or proton density
    !> @param[in]       d1rhoq      First derivative of neutron or proton
    !>                              density
    !> @param[in]       d2rhoq      Second derivative of neutron or proton
    !>                              density
    !> @param[in]       d3rhoq      Third derivative of neutron or proton
    !>                              density
    !> @param[in]       d4rhoq      Fourth derivative of neutron or proton
    !>                              density
    !> @param[in]       fq          Effective mass for neutrons or protons
    !> @param[in]       d1fq        First derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in]       d2fq        Second derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in]       d3fq        Third derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in]       d4fq        Fourth derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in,out]   tau4nospinq Array for fourth-order contributions (no
    !>                              spin-orbit)
    !> @param[in,out]   t4contq     Array for individual fourth-order
    !>                              contributions
    subroutine calc_tau_4_no_spin(rhoq,d1rhoq,d2rhoq,d3rhoq,d4rhoq,fq,d1fq,d2fq,d3fq,d4fq,tau4nospinq,t4contq)
        implicit none
        real(kind=dp),  dimension(1:n),     intent(in)      :: rhoq,d1rhoq,d2rhoq,d3rhoq,d4rhoq,fq,d1fq,d2fq,d3fq,d4fq
        real(kind=dp),  dimension(1:n),     intent(inout)   :: tau4nospinq
        real(kind=dp),  dimension(1:n,1:3), intent(inout)   :: t4contq
        
        ! All arrays (including r) are 1:n
        
                ! Terms containing gradients of density
        t4contq(:,1) = (3._dp*pi**2)**(-fr23) * rhoq**fr13 / 4320._dp * &
                & ( 24._dp/rhoq * (d4rhoq + 4._dp*d3rhoq/r) &
                &   - 8._dp/rhoq**2 * (11._dp*d3rhoq*d1rhoq + 7._dp*d2rhoq**2 + 36._dp*d2rhoq*d1rhoq/r - (d1rhoq/r)**2) &
                &   + fr83/rhoq**3 * (81._dp*d2rhoq*d1rhoq**2 + 70._dp*d1rhoq**3/r) &
                &   - 96._dp*(d1rhoq/rhoq)**4 )
        
                ! Terms containing gradients of effective mass
        t4contq(:,2) = (3._dp*pi**2)**(-fr23) * rhoq**fr13 / 4320._dp * &
                & ( - 36._dp/fq * (d4fq + 4._dp*d3fq/r) &
                &   + 18._dp/fq**2 * (4._dp*d3fq*d1fq + 3._dp*d2fq**2 + 4._dp*d2fq*d1fq/r - 4._dp*(d1fq/r)**2) &
                &   - 144._dp/fq**3 * d2fq*d1fq**2 &
                &   + 54._dp*(d1fq/fq)**4 )
        
                ! Terms containing gradients of density and of effective mass
        t4contq(:,3) = (3._dp*pi**2)**(-fr23) * rhoq**fr13 / 4320._dp * &
                & ( 12._dp/(fq*rhoq) * ( 3._dp*d1fq*d3rhoq + 2._dp*d2fq*d2rhoq - 2._dp*d3fq*d1rhoq &
                &                           + 6._dp*d1fq*d2rhoq/r - 4._dp*d2fq*d1rhoq/r + 2._dp*d1fq*d1rhoq/r**2 ) &
                &   + 12._dp/(rhoq*fq**2) * (d1fq*d2fq*d1rhoq - 2._dp*d2rhoq*d1fq**2 - 2._dp*d1rhoq*d1fq**2/r) &
                &   - 4._dp/(fq*rhoq**2) * (11._dp*d2fq*d1rhoq**2 + 32._dp*d1fq*d1rhoq*d2rhoq + 26._dp*d1fq*d1rhoq**2/r) &
                &   + 30._dp*((d1fq*d1rhoq)/(fq*rhoq))**2 &
                &   + fr260_3*d1fq*(d1rhoq/rhoq)**3/fq )
        
        ! Store total fourth-order contribution (no spin-orbit)
        tau4nospinq(:) = t4contq(:,1) + t4contq(:,2) + t4contq(:,3)
        
    end subroutine calc_tau_4_no_spin
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the fourth-order contributions to the kinetic
    !> energy density (spin-orbit contributions) \f$\textbf{tau}_q\f$ for
    !> neutrons or protons
    !> @param[in]       rhoq        Neutron or proton density
    !> @param[in]       d1rhoq      First derivative of neutron or proton
    !>                              density
    !> @param[in]       fq          Effective mass for protons or neutrons
    !> @param[in]       d1fq        First derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in]       d2fq        Second derivative of effective mass for
    !>                              neutrons or protons
    !> @param[in]       d1Aq        First derivative of "composite" neutron or
    !>                              proton density \f$A_q\f$
    !> @param[in]       d2Aq        Second derivative of "composite" neutron or
    !>                              proton density \f$A_q\f$
    !> @param[in]       d3Aq        Third derivative of "composite" neutron or
    !>                              proton density \f$A_q\f$
    !> @param[in]       q           Isospin
    !> @param[in,out]   tau_4_so_q  Array for fourth-order contributions (spin-
    !>                              orbit)
    !> @param[in,out]   t4contq     Array for individual fourth-order
    !>                              contributions
    subroutine calc_tau_4_so(rhoq,d1rhoq,fq,d1fq,d2fq,d1Aq,d2Aq,d3Aq,q,tau_4_so_q,t4contq)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,d1rhoq,fq,d1fq,d2fq,d1Aq,d2Aq,d3Aq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(1:n), intent(inout)   :: tau_4_so_q,t4contq
        
        ! All arrays (including r) are 1:n
        tau_4_so_q = (3._dp*pi**2)**(-fr23) * fr1_16*(W0/hbar2_2m_q(q))**2 * rhoq**fr13 / fq**2 * &
                & ( fr12 * (2._dp*d1Aq*d3Aq + d2Aq**2 + 6._dp*d1Aq*d2Aq/r - (d1Aq/r)**2) &
                &   - d1Aq/fq * (d2fq*d1Aq + 2._dp*d1fq*d2Aq + 2._dp*d1fq*d1Aq/r) &
                &   + fr34*(d1Aq/fq)**2 * (2._dp*d1fq**2 - fr14*(W0*d1Aq/hbar2_2m_q(q))**2) &
                &   + fr13*d1rhoq*d1Aq/(fq*rhoq) * (fq*(d2Aq + d1Aq/r) - d1fq*d1Aq) )
        
        ! Store fourth-order contribution (spin-orbit)
        t4contq = tau_4_so_q
        
    end subroutine calc_tau_4_so
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the second-order contributions to the spin
    !> current density \f$\textbf{J}_q\f$ for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       Wq      Spin-orbit field
    !> @param[in]       fq      Effective mass for neutrons or protons
    !> @param[in]       q       Isospin
    !> @paran[in,out]   J2q     Array for second-order contributions
    subroutine calc_J_2_q(rhoq,Wq,fq,q,J2q)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,Wq,fq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(1:n), intent(inout)   :: J2q
        
        ! All arrays are 1:n
        J2q = -1._dp/hbar2_2m_q(q) * rhoq*Wq/fq

    end subroutine calc_J_2_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the fourth-order contributions to the spin
    !> current density \f$\textbf{J}_q\f$ for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       d1rhoq  First derivative of neutron or proton density
    !> @param[in]       fq      Effective mass for protons or neutrons
    !> @param[in]       d1fq    First derivative of effective mass for neutrons
    !>                          or protons
    !> @param[in]       d2fq    Second derivative of effective mass for
    !>                          neutrons or protons
    !> @param[in]       d1Aq    First derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d2Aq    Second derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d3Aq    Third derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       q       Isospin
    !> @param[in,out]   J4q     Array for fourth-order contributions
    subroutine calc_J_4_q(rhoq,d1rhoq,fq,d1fq,d2fq,d1Aq,d2Aq,d3Aq,q,J4q)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,d1rhoq,fq,d1fq,d2fq,d1Aq,d2Aq,d3Aq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(1:n), intent(inout)   :: J4q
        
        ! All arrays (including r) are 1:n
        J4q = (3._dp*pi**2)**(-fr23) * fr18*W0/hbar2_2m_q(q) * rhoq**fr13 / fq * &
            & ( - d3Aq - 2._dp*d2Aq/r + 2._dp*d1Aq/r**2 &
            &   + (d2fq*d1Aq + d1fq*d2Aq + d1fq*d1Aq/r) / fq &
            &   - d1Aq/fq**2 * (d1fq**2 - (fr12*W0/hbar2_2m_q(q))**2 * d1Aq**2) &
            &   - fr13*d1rhoq/rhoq * (d2Aq + d1Aq/r - d1fq*d1Aq/fq) )
        
    end subroutine calc_J_4_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the second-order contributions to the divergence
    !> of the current density \f$\textbf{J}_q\f$ for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       d1rhoq  First derivative of neutron or proton density
    !> @param[in]       fq      Effective mass for protons or neutrons
    !> @param[in]       d1fq    First derivative of effective mass for neutrons
    !>                          or protons
    !> @param[in]       d1Aq    First derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d2Aq    Second derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       q       Isospin
    !> @param[in,out]   divJ2q  Array for second-order contributions
    subroutine calc_div_J_2_q(rhoq,d1rhoq,fq,d1fq,d1Aq,d2Aq,q,divJ2q)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,d1rhoq,fq,d1fq,d1Aq,d2Aq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(1:n), intent(inout)   :: divJ2q
        
        ! All arrays (including r) are 1:n
        divJ2q = -fr12*W0/(hbar2_2m_q(q)*fq) * (rhoq*(2._dp*d1Aq/r + d2Aq) + d1rhoq*d1Aq - rhoq*d1fq*d1Aq/fq)
        
    end subroutine calc_div_J_2_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the fourth-order contributions to the divergence
    !> of the current density \f$\textbf{J}_q\f$ for neutrons or protons
    !> @param[in]       rhoq    Neutron or proton density
    !> @param[in]       d1rhoq  First derivative of neutron or proton density
    !> @param[in]       d2rhoq  Second derivative of neutron or proton density
    !> @param[in]       fq      Effective mass for protons or neutrons
    !> @param[in]       d1fq    First derivative of effective mass for neutrons
    !>                          or protons
    !> @param[in]       d2fq    Second derivative of effective mass for neutrons
    !>                          or protons
    !> @param[in]       d3fq    Third derivative of effective mass for neutrons
    !>                          or protons
    !> @param[in]       d1Aq    First derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d2Aq    Second derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d3Aq    Third derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       d4Aq    Fourth derivative of "composite" neutron or
    !>                          proton density \f$A_q\f$
    !> @param[in]       q       Isospin
    !> @param[in,out]   divJ4q  Array for fourth-order contributions
    subroutine calc_div_J_4_q(rhoq,d1rhoq,d2rhoq,fq,d1fq,d2fq,d3fq,d1Aq,d2Aq,d3Aq,d4Aq,q,divJ4q)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)      :: rhoq,d1rhoq,d2rhoq,fq,d1fq,d2fq,d3fq,d1Aq,d2Aq,d3Aq,d4Aq
        integer,                        intent(in)      :: q
        real(kind=dp),  dimension(1:n), intent(inout)   :: divJ4q
        real(kind=dp) :: W0fac !< Factor containing \f$W_0\f$ that appears several times
        
        W0fac = fr12*W0/hbar2_2m_q(q)
        
        ! All arrays (including r) are 1:n
        divJ4q = (3._dp*pi**2)**(-fr23) * fr14*W0fac * rhoq**fr13 / fq * &
            & ( - (d4Aq + 4._dp*d3Aq/r) &
            &   + ( d3fq*d1Aq + 2._dp*d2fq*d2Aq + 2._dp*d1fq*d3Aq &
                    & + 5._dp*d1fq*d2Aq/r + 3._dp*d2fq*d1Aq/r - d1fq*d1Aq/r ) / fq &
            &   - ( 3._dp*d2Aq*d1fq**2 + 4._dp*d1fq*d2fq*d1Aq + 4._dp*d1Aq*d1fq**2/r &
                    & - (W0fac*d1Aq)**2 * (3._dp*d2Aq + 2._dp*d1Aq/r) ) / fq**2 &
            &   + (3._dp*d1fq*d1Aq/fq**3) * (d1fq**2 - (W0fac*d1Aq)**2) &
            &   - fr13/rhoq * ( 2._dp*d1rhoq*d3Aq + d2rhoq*d2Aq &
                    & + d2rhoq*d1Aq/r + 5._dp*d1rhoq*d2Aq/r + d1rhoq*d1Aq/r**2 ) &
            &   + fr13/(fq*rhoq) * ( 3._dp*d1fq*d1rhoq*d2Aq + d1fq*d2rhoq*d1Aq &
                    & + 2._dp*d2fq*d1rhoq*d1Aq + 4._dp*d1fq*d1rhoq*d1Aq/r ) &
            &   - fr13*d1rhoq*d1Aq/(rhoq*fq**2) * (3._dp*d1fq**2 - (W0fac*d1Aq)**2) &
            &   + fr29*(d1rhoq/rhoq)**2 * (d2Aq + d1Aq/r - d1fq*d1Aq/fq) )
        
    end subroutine calc_div_J_4_q
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to calculate the central potential \f$U_q\f$ for neutrons or
    !> protons, for standard Skyrme or for BSk forces with extra terms
    !> @param[in] rhot      Total density
    !> @param[in] rhon      Neutron density
    !> @param[in] rhop      Proton density
    !> @param[in] delrhot   First derivative of total density
    !> @param[in] delrhon   First derivative of neutron density
    !> @param[in] delrhop   First derivative of proton density
    !> @param[in] del2rhot  Laplacian of total density
    !> @param[in] del2rhon  Laplacian of neutron density
    !> @param[in] del2rhop  Laplacian of proton density
    !> @param[in] taut      Total kinetic density
    !> @param[in] taun      Neutron kinetic density
    !> @param[in] taup      Proton kinetic density
    !> @param[in] divJt     Divergence of total spin current density
    !> @param[in] divJq     Divergence of proton or neutron spin current
    !>                      density
    !> @param[in] q         Isospin
    function calc_U_q(rhot,rhon,rhop,delrhot,delrhon,delrhop,del2rhot,del2rhon,del2rhop,taut,taun,taup,divJt,divJq,q)
        implicit none
        real(kind=dp) :: calc_U_q
        real(kind=dp),  intent(in)  :: rhot,rhon,rhop,delrhot,delrhon,delrhop,del2rhot,del2rhon,del2rhop
        real(kind=dp),  intent(in)  :: taut,taun,taup,divJt,divJq
        integer,        intent(in)  :: q
        real(kind=dp),  dimension(0:1)  :: rhonp !< Storage for neutron and proton density together
        real(kind=dp),  dimension(0:1)  :: delrhonp !< Storage for neutron and proton density derivatives together
        real(kind=dp),  dimension(0:1)  :: del2rhonp !< Storage for neutron and proton density second derivatives together
        real(kind=dp),  dimension(0:1)  :: taunp !< Storage for neutron and proton kinetic density together
        
        rhonp(0) = rhon
        rhonp(1) = rhop
        delrhonp(0) = delrhon
        delrhonp(1) = delrhop
        del2rhonp(0) = del2rhon
        del2rhonp(1) = del2rhop
        taunp(0) = taun
        taunp(1) = taup
        
        ! BSk forces with extra terms
        if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 18) then
            calc_U_q = t0 * ( (1._dp+fr12*x0)*rhot - (fr12+x0)*rhonp(q) ) &
                & + fr14*t1 * ( (1._dp+fr12*x1)*(taut-fr32*del2rhot) - (fr12+x1)*(taunp(q)-fr32*del2rhonp(q)) ) &
                & + fr14*t2 * ( (1._dp+fr12*x2)*(taut+fr12*del2rhot) + (fr12+x2)*(taunp(q)+fr12*del2rhonp(q)) ) &
                & + fr1_12*t3 * ( ( 1._dp+fr12*x3)*(2._dp+alpha)*rhot**(alpha+1._dp) &
                    & - (fr12+x3)*(2._dp*rhot**alpha*rhonp(q) + alpha*rhot**(alpha-1._dp)*(rhon**2+rhop**2)) ) &
                & + fr18*t4 * &
                    & ( (1._dp+fr12*x4) * rhot**(beta-1._dp) * &
                        & ( 2._dp*(1._dp+beta)*rhot*taut &
                        & - (2._dp*beta+3._dp)*(fr12*beta*delrhot**2 + rhot*del2rhot) ) &
                    & + (fr12+x4) * rhot**(beta-2._dp) * &
                        & ( 3._dp*beta*rhot*delrhot*delrhonp(q) &
                        & + 3._dp*rhot**2*del2rhonp(q) &
                        & - 2._dp*rhot**2*taunp(q) &
                        & + beta*(beta-1._dp)*rhonp(q)*delrhot**2 &
                        & + beta*rhot*rhonp(q)*del2rhot &
                        & - fr12*beta*rhot * &
                            & ( delrhonp(0)**2 + 4._dp*rhonp(0)*taunp(0) - 2._dp*rhonp(0)*del2rhonp(0) &
                            & + delrhonp(1)**2 + 4._dp*rhonp(1)*taunp(1) - 2._dp*rhonp(1)*del2rhonp(1) ) ) ) &
                & + fr14*t5 * rhot**(gamma-1._dp) * &
                    & ( (1._dp+fr12*x5) * ( (1._dp+gamma)*rhot*taut + fr14*gamma*delrhot**2 + fr12*rhot*del2rhot ) &
                    & + (fr12+x5) * &
                        & ( rhot*taunp(q) &
                        & + fr12*rhot*del2rhonp(q) &
                        & + gamma * (rhonp(0)*taunp(0) - fr14*delrhonp(0)**2 + rhonp(1)*taunp(1) - fr14*delrhonp(1)**2) &
                        & + fr12*gamma*delrhot*delrhonp(q) ) ) &
                & - fr12*W0 * (divJt + divJq)
        
        ! Standard Skyrme forces
        else
            calc_U_q = 2._dp * (B1*rhot + B2*rhonp(q)) &
                & + B3*taut + B4*taunp(q) &
                & + 2._dp * (B5*del2rhot + B6*del2rhonp(q)) &
                & + (2._dp+sigma) * B7 * rhot**(sigma+1._dp) &
                & + B8 * ( sigma*rhot**(sigma-1._dp)*(rhon**2+rhop**2) + 2._dp*rhot**sigma*rhonp(q) ) &
                & + B9 * (divJt + divJq)
            
        end if
        
    end function calc_U_q
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate matter density derivatives, and all fields: 
    !> effective masses (and derivatives), spin-orbit fields, spin current
    !> densities (and derivatives)
    subroutine eval_fields()
        implicit none
        integer :: q !< Loop index
        
        ! Derivatives for individual densities
        isospin_1: do q=0, 1
            ! Extrapolate back to enable derivatives to be taken
            call extrapolate_back_3(rho_q(:,q))
            
            ! 1st derivative
            call d1(n,dr,rho_q(-1:n,q),del_rho_q(:,q))

            ! Laplacian in spherical symmetry
            call lap_sphe_symm(rho_q(:,q),del_rho_q(:,q),del2_rho_q(:,q))
            
            ! Evaluate effective mass using densities
            call calc_f_q(rho_t(1:n),rho_q(1:n,q),q,f_q(1:n,q))
            
            ! Extrapolate back to enable derivatives to be taken
            call extrapolate_back_3(f_q(:,q))
            
            ! 1st derivative of effective mass
            call d1(n,dr,f_q(-1:n,q),del_f_q(:,q))
            
            ! Laplacian in spherical symmetry
            call lap_sphe_symm(f_q(:,q),del_f_q(:,q),del2_f_q(:,q))
            
        end do isospin_1
        
        ! Extrapolate back to enable derivatives of total density to be taken
        call extrapolate_back_3(rho_t(:))
        
        ! 1st derivative of total density in spherical symmetry
        del_rho_t(:) = del_rho_q(:,0) + del_rho_q(:,1)
        
        ! Laplacian of total density in spherical symmetry
        call lap_sphe_symm(rho_t(:),del_rho_t(:),del2_rho_t(:))
        
        ! 2nd order contributions
        if (etf_order .eq. 2 .or. etf_order .eq. 4) then
            isospin_2: do q=0, 1
                ! Evaluate spin-orbit field, and spin current density (2nd
                ! order contribution)
                call calc_W_q(del_rho_t(:),del_rho_q(:,q),W_q(:,q))
                call calc_J_2_q(rho_q(1:n,q),W_q(:,q),f_q(1:n,q),q,J_2_q(1:n,q))
                
                ! Store contribution of 2nd order term for spin current density
                J_q(1:n,q) = J_2_q(1:n,q)
                
                ! New composite density A_q
                A_q(1:n,q) = rho_t(1:n) + rho_q(1:n,q)
                
                ! Extrapolate back to enable derivatives to be taken
                call extrapolate_back_3(A_q(:,q))
                
                ! Derivatives of A_q
                call d1(n,dr,A_q(-1:n,q),d1_A_q(:,q))
                call d2(n,dr,A_q(-1:n,q),d2_A_q(:,q))
                
                ! 2nd order contribution to divergence of spin current density
                call calc_div_J_2_q(rho_q(1:n,q),del_rho_q(:,q), &
                                    & f_q(1:n,q),del_f_q(:,q), &
                                    & d1_A_q(:,q),d2_A_q(:,q),q,div_J_2_q(:,q))
                
                ! Add contribution of 2nd order terms
                div_J_q(1:n,q) = div_J_q(1:n,q) + div_J_2_q(1:n,q)
                
                if (etf_order .eq. 4) then
                    ! Extra derivatives of densities
                    call d2(n,dr,rho_q(-1:n,q),d2_rho_q(:,q)) ! Different from Laplacian
                    call d3(rho_q(:,q),d3_rho_q(:,q))
                    call d4(rho_q(:,q),d4_rho_q(:,q))
                    
                    ! Extra derivatives of effective masses
                    call d2(n,dr,f_q(-1:n,q),d2_f_q(:,q)) ! Different from Laplacian
                    call d3(f_q(:,q),d3_f_q(:,q))
                    call d4(f_q(:,q),d4_f_q(:,q))
                    
                    ! Derivatives of A_q
                    call d3(A_q(:,q),d3_A_q(:,q))
                    call d4(A_q(:,q),d4_A_q(:,q))
                    
                    ! 4th order contribution to spin current density
                    call calc_J_4_q(rho_q(1:n,q),del_rho_q(:,q),f_q(1:n,q),del_f_q(:,q),d2_f_q(:,q), &
                                    & d1_A_q(:,q),d2_A_q(:,q),d3_A_q(:,q),q,J_4_q(:,q))
                    
                    ! Add contribtion of 4th order terms
                    J_q(1:n,q) = J_q(1:n,q) + J_4_q(1:n,q)
                    
                    ! 4th order contribution to divergence of spin current density
                    call calc_div_J_4_q(rho_q(1:n,q),del_rho_q(:,q),d2_rho_q(:,q), &
                                        & f_q(1:n,q),del_f_q(:,q),d2_f_q(:,q),d3_f_q(:,q), &
                                        & d1_A_q(:,q),d2_A_q(:,q),d3_A_q(:,q),d4_A_q(:,q),q,div_J_4_q(:,q))
                    
                    ! Add contribution of 4th order terms
                    div_J_q(1:n,q) = div_J_q(1:n,q) + div_J_4_q(1:n,q)
                end if
                
                ! Extrapolate back to enable derivatives to be taken
!                 call extrapolate_back_3(J_q(:,q))
                
                ! 1st derivative of spin current density
!                 call d1(n,dr,J_q(-1:n,q),del_J_q(:,q))
                
                ! Divergence of spin current density
!                 call div_sphe_symm(J_q(:,q),del_J_q(:,q),div_J_q(:,q))
            end do isospin_2
            
            ! Total spin current density
            J_t(-2:n) = J_q(-2:n,0) + J_q(-2:n,1)
            
            ! Spin-orbit energy
            e_so = WS_integral(-B9 * (J_q(1:n,0)*del_rho_q(1:n,0) + J_q(1:n,1)*del_rho_q(1:n,1) + J_t(1:n)*del_rho_t(1:n)),n)
            
            ! J^2 terms contribution to spin-orbit energy
            if (J2_terms .eqv. .true.) then
                e_so = e_so &
                    & + WS_integral(-fr1_16*(t1*x1 + t2*x2)*J_t(1:n)**2 + fr1_16*(t1-t2)*(J_q(1:n,0)**2 + J_q(1:n,1)**2),n)
            end if
            
            ! Divergence of total spin current density
            div_J_t(1:n) = div_J_q(1:n,0) + div_J_q(1:n,1)
            
!             ! Extrapolate back to enable derivatives to be taken
!             call extrapolate_back_3(J_t(:))
!             
!             ! 1st derivative of total spin current density
!             call d1(n,dr,J_t(-1:n),del_J_t(:))
!             
!             ! Divergence of total spin current density
!             call div_sphe_symm(J_t(:),del_J_t(:),div_J_t(:))
            
        end if
        
    end subroutine eval_fields
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the kinetic energy density \f$\tau_{ETF}\f$ with
    !> the (extended) Thomas-Fermi approximation at order specified by
    !> "etf_order"
    subroutine eval_tau_ETF()
        implicit none
        integer :: q !< Loop index
        
        ! Calculate kinetic energy density at different radii, for neutrons and
        ! protons
        isospin: do q=0, 1
            ! 0th order contribution
            call calc_tau_TF(rho_q(1:n,q),tau_TF_q(:,q))
            
            ! Store contribtion of 0th order term
            tau_ETF_q(1:n,q) = tau_TF_q(1:n,q)
            
            ! 2nd order contributions
            if (etf_order .eq. 2 .or. etf_order .eq. 4) then
                ! Local contribution
                call calc_tau_2_L(rho_q(1:n,q),del_rho_q(:,q),del2_rho_q(:,q),tau_2_L_q(:,q),tau_2_cont_q(:,1,q))
                
                ! Non-local contribution
                call calc_tau_2_NL(rho_q(1:n,q),del_rho_q(:,q),f_q(1:n,q),del_f_q(:,q),del2_f_q(:,q),W_q(:,q), &
                                    & q,tau_2_NL_q(:,q),tau_2_cont_q(:,2:4,q))
                
                ! Add contribution of 2nd order terms
                tau_ETF_q(1:n,q) = tau_ETF_q(1:n,q) + tau_2_L_q(1:n,q) + tau_2_NL_q(1:n,q)
            end if
            
            ! 4th order contributions
            if (etf_order .eq. 4) then
                ! Spin-independent part
                call calc_tau_4_no_spin(rho_q(1:n,q),del_rho_q(:,q),d2_rho_q(:,q),d3_rho_q(:,q),d4_rho_q(:,q), &
                                        & f_q(1:n,q),del_f_q(:,q),d2_f_q(:,q),d3_f_q(:,q),d4_f_q(:,q), &
                                        & tau_4_no_spin_q(:,q),tau_4_cont_q(:,1:3,q))
                
                ! Spin-dependent part
                call calc_tau_4_so(rho_q(1:n,q),del_rho_q(:,q),f_q(1:n,q),del_f_q(:,q),d2_f_q(:,q), &
                                    & d1_A_q(:,q),d2_A_q(:,q),d3_A_q(:,q),q,tau_4_so_q(:,q),tau_4_cont_q(:,4,q))
                
                ! Add contribution of 4th order terms
                tau_ETF_q(1:n,q) = tau_ETF_q(1:n,q) + tau_4_no_spin_q(1:n,q) + tau_4_so_q(1:n,q)
            end if
            
            ! Kinetic energies
            e_kinetic_q(q) = WS_integral(hbar2_2m_q(q)*tau_ETF_q(1:n,q),n)
        end do isospin
        
        ! Total kinetic energy density
        tau_ETF_t(1:n) = tau_ETF_q(1:n,0) + tau_ETF_q(1:n,1)
        
        ! Total kinetic energy
        e_kinetic_t = sum(e_kinetic_q)
        
    end subroutine eval_tau_ETF
    
    
    
    !> @author A. Pastore, M. Shelley
    !> @brief
    !> Subroutine to calculate the charge density from the matter densities.
    !> Modified by M. Shelley.
    !> @param[in]       rho     Neutron and proton densities
    !> @param[in]       Neutr   Number of neutrons
    !> @param[in]       Nprot   Number of protons
    !> @param[in]       Ngrid1  Number of mesh points
    !> @param[in]       del1    Step size
    !> @param[in,out]   rhoch   Array for charge density
    subroutine charge(rho,Neutr,Nprot,Ngrid1,del1,rhoch)
        implicit none
        integer,        intent(in)      :: Ngrid1
        real(kind=dp),  intent(in)      :: rho(ngrid1,0:1),Neutr,Nprot
        real(kind=dp),  intent(in)      :: del1
        real(kind=dp),  intent(inout)   :: rhoch(Ngrid1)
        integer :: i,j
        real(kind=dp) :: amass,a,aux,bta,hbmedio,ct0,r,rp,hb(0:1)
        real(kind=dp) :: apc,spi,bn1,bn2,an1,an2,an1c,an2c
        real(kind=dp) :: aux2,exn1,exn2,ap1,ap2,ap3,ap1c,ap2c,ap3c
        
        bn1=0.387_dp
        bn2=0.467_dp
        
        amass=Nprot+Neutr
!         hb(:)=1._dp/dmshb0(:)
        hb(:) = hbar2_2m_q(:)
        hbmedio=2._dp/(1._dp/hb(0)+1._dp/hb(1))
        ct0=2*hbmedio
        bta=41._dp/(ct0*amass**(1._dp/3._dp))
        
        spi=sqrt(pi)
        rhoch=0._dp
        
        a=0.53401666_dp!!sqrt(2._dp/3._dp)*ravp
        a=a-1._dp/(amass*bta)
        ap1=0.431566_dp -1._dp/(amass*bta)
        ap2=0.13914_dp -1._dp/(amass*bta)
        ap3=1.52554_dp-1._dp/(amass*bta)  
        
        ap1c=1._dp/(sqrt(ap1)*spi)
        ap2c=1._dp/(sqrt(ap2)*spi)
        ap3c=1._dp/(sqrt(ap3)*spi)
        
        an1=bn1-1._dp/(amass*bta)
        an2=bn2-1._dp/(amass*bta)
        
        apc=1._dp/(sqrt(a)*spi)
        an1c=1._dp/(sqrt(an1)*spi)
        an2c=1._dp/(sqrt(an2)*spi)
            
!         write(*,*)a,apc,an1c,an2c,an1,an2
        
        do i=1,Ngrid1
            r=i*del1
            do j=1,Ngrid1
                rp=j*del1
                exn1=exp(-(r-rp)**2/an1)-exp(-(r+rp)**2/an1)
                exn2=exp(-(r-rp)**2/an2)-exp(-(r+rp)**2/an2)
                aux2=(exn1*an1c-exn2*an2c)*rho(j,0)
                
                aux=rho(j,1)*apc*(exp(-(r-rp)**2/a)-exp(-(r+rp)**2/a))
                rhoch(i)= rhoch(i)+ rp*(aux+aux2)*del1/r
            enddo
!             write(*,*)i,rhoch(i)
            
!             if(rhoch(i).lt.0._dp)rhoch(i)=0._dp
            if (rhoch(i) .lt. 0._dp) then
                rhoch(i:n) = 0._dp
                exit
            end if
        enddo
        
    end subroutine charge
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the Coulomb potentials and energy densities,
    !> using either the proton or charge density
    !> @param[in]   prot_dens       Density to use (proton or charge)
    !> @param[in]   calc_exchange   Whether to calculate exchange potential
    subroutine eval_coulomb(prot_dens,calc_exchange)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in) :: prot_dens
        logical,                        intent(in) :: calc_exchange
        integer :: i !< Loop index
        real(kind=dp) :: int1,int2 !< Integrals for direct potential
        
        !!! Direct Coulomb potential !!!
        int1 = 0._dp
        int2 = 0._dp
        
        ! Evaluate on the mesh the three integrals needed for the direct
        ! potential
            ! Euler method is sufficiently accurate (~10^-3 keV/A) for Coulomb
        integrals: do i=1, n
            int1 = int1 + r(i)**2 * prot_dens(i) * 4._dp*pi * dr
            int2 = int2 + r(i) * prot_dens(i) * 4._dp*pi * dr
            v_c_di(i) = (int1 / r(i)) - int2
        end do integrals
        
        ! Combine the three integrals to get the total direct potential
            ! Using last value of "int2" (i.e. the 0 to infinity integral)
        v_c_di(1:n) = e2 * (v_c_di(1:n) + int2)
        
        
        !!! Exchange Coulomb potential !!!
        if (calc_exchange .eqv. .true.) then
            v_c_ex(1:n) = -e2 * (3._dp/pi * prot_dens(1:n))**fr13
        
        else
            ! Set to 0 if not calculated
            v_c_ex(1:n) = 0._dp
            
        end if
        
        
        !!! Energy densities and energies
        ! Direct Coulomb energy density
        e_density_c_di(1:n) = prot_dens(1:n) * v_c_di(1:n) / 2._dp
        
        ! Direct Coulomb energy
        e_coulomb_di = WS_integral(e_density_c_di(1:n),n)
        
        ! Exchange Coulomb energy density
        e_density_c_ex(1:n) = fr34 * prot_dens(1:n) * v_c_ex(1:n)
        
        ! Exchange Coulomb energy
        e_coulomb_ex = WS_integral(e_density_c_ex(1:n),n)
        
        ! Coulomb energy
        e_coulomb = WS_integral(e_density_c_di(1:n) + e_density_c_ex(1:n),n)
        
    end subroutine eval_coulomb
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the proton-electron potential which contributes
    !> to \f$U_q\f$, all energy contributions, and the pressure, coming from a
    !> homogeneous electron gas. Also calculates the electron chemical
    !> potential.
    subroutine eval_electrons()
        implicit none
        real(kind=dp) :: np !< Number of protons
        real(kind=dp) :: k_fe !< Electron fermi momentum
        real(kind=dp) :: x !< Relativity parameter
        real(kind=dp) :: phi !< Relativity function for electron-electron potential energy
        real(kind=dp) :: pe_term2 !< Term appearing in proton-electron energy/ chemical potential
        
        !!! Proton-electron potential !!!
        ! Number of protons
        np = N_q(1)
        
        ! Potential
        v_c_pe(1:n) = fr12 * np * e2 / r_max * ((r(1:n)/r_max)**2 - 3._dp)
        
        
        !!! Electron kinetic energy !!!
        ! Electron density
        rho_e = np / (fr43 * pi * r_max**3)
        
        ! Electron fermi momentum
        k_fe = (3._dp*pi**2 * rho_e)**fr13
        
        ! Relativity parameter (x-->0 as rho-->0)
        x = (hbarc / mec2) * k_fe
        
        ! Kinetic energy for relativistic electrons (excluding electron rest mass)
        e_kinetic_e = np * mec2 * (fr38/x**3 * (x*(1._dp + 2._dp*x**2)*sqrt(1._dp+x**2) - asinh(x)) - 1._dp)
        
        
        !!! Pressure !!!
        pressure_e = rho_e * mec2 * fr38/x**3 * (x*(fr23*x**2 - 1._dp)*sqrt(1._dp+x**2) + asinh(x))
        
        ! Coulomb exchange pressure
        pressure_ex = fr18 * e2 * (3._dp/pi)**fr13 * rho_e**fr43
        
        
        !!! Electron-electron potential energy !!!
        ! Relativity function \f$\Phi(x)\f$
        phi = -fr12/x**4 * (3._dp*x**2 + x**4 - 6._dp*x*sqrt(1._dp+x**2)*asinh(x) + 3._dp*asinh(x)**2)
        
        ! Potential energy for relativistic electrons
        e_coulomb_e = fr35 * e2 * np**2 / r_max * (1._dp - fr54 * (fr32/(pi*np))**fr23 * phi)
        
        
        !!! Proton-electron potential energy !!!
        ! Calculate quantity that also appears in proton-electron Coulomb chemical potential
        pe_term2 = (fr12 * e2 / r_max**3) * WS_integral(rho_q(1:n,1)*r(1:n)**2,n)
        
        ! Potential energy
        e_coulomb_pe = (-fr32 * np**2 * e2 / r_max) + np * pe_term2
        
        
        !!! Electron chemical potential !!!
        if (calc_chem_pots .eqv. .true.) then
            ! Proton-electron Coulomb chemical potential
            mu_c = -fr3_10 * e2 * np / r_max + pe_term2
            
            ! Electron chemical potential
            mu_e = sqrt((hbarc*k_fe)**2 + mec2**2) - mec2! + mu_c
!             mu_e = sqrt((hbarc*k_fe)**2 + mec2**2) - mec2 + fr12*e2*(3._dp*rho_e/pi)**fr13
        end if
        
    end subroutine eval_electrons
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the central fields \f$U_q\f$
    subroutine eval_U_q()
        implicit none
        integer :: i,q !< Loop indices
        
        isospin: do q=0, 1
            ! Central fields
            radii: do i=1, n
                U_q(i,q) = calc_U_q( rho_t(i),rho_q(i,0),rho_q(i,1), &
                                    & del_rho_t(i),del_rho_q(i,0),del_rho_q(i,1), &
                                    & del2_rho_t(i),del2_rho_q(i,0),del2_rho_q(i,1), &
                                    & tau_ETF_t(i),tau_ETF_q(i,0),tau_ETF_q(i,1), &
                                    & div_J_t(i),div_J_q(i,q),q )
            end do radii
        end do isospin
        
        ! Add direct and exchange Coulomb potentials to proton central potential
        if (coulomb_on .eqv. .true.) then
            U_q(1:n,1) = U_q(1:n,1) + v_c_di(1:n) + v_c_ex(1:n)
        end if
        
        ! Add proton-electron Coulomb potential
        if (electrons_on .eqv. .true.) then
            U_q(1:n,1) = U_q(1:n,1) + v_c_pe(1:n)
        end if
        
    end subroutine eval_U_q
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the Skyrme energy density, after first
    !> separately evaluating the field energy density
    subroutine eval_skyrme_energy_density()
        implicit none
        integer :: i !< Loop indices
        real(kind=dp) :: rhot,rhon,rhop,taut,taun,taup !< Dummy variables for clarity
        real(kind=dp) :: drhot,drhon,drhop!,Jt,Jn,Jp !< Dummy variables for clarity
!         integer :: q
!         real(kind=dp), dimension(1:n,0:1) :: cont1, cont2
        
        !!! Calculate field energy density at different radii !!!
        ! BSk forces
        if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 18) then
            radii_1: do i=1, n
                rhot = rho_t(i)
                rhon = rho_q(i,0)
                rhop = rho_q(i,1)
                taut = tau_ETF_t(i)
                taun = tau_ETF_q(i,0)
                taup = tau_ETF_q(i,1)
                
                
!                 taut = 0._dp
!                 taun = 0._dp
!                 taup = 0._dp
                
                
                drhot = del_rho_t(i)
                drhon = del_rho_q(i,0)
                drhop = del_rho_q(i,1)
!                 Jt = J_t(i)
!                 Jp = J_q(i,0)
!                 Jn = J_q(i,1)
                
                e_density_field(i) = &
                        ! Standard Skyrme terms (with sigma --> alpha)
                        & + fr14*t0 * ( (2._dp+x0)*rhot**2 - (2._dp*x0+1._dp)*(rhon**2+rhop**2) ) &
                        & + fr1_24*t3*rhot**alpha * ((2._dp+x3)*rhot**2 - (2._dp*x3+1._dp)*(rhon**2+rhop**2) ) &
                        & + fr18 * (t1*(2._dp+x1) + t2*(2._dp+x2)) * rhot*taut &
                        & + fr18 * (t2*(2._dp*x2+1._dp) - t1*(2._dp*x1+1._dp)) * (rhon*taun + rhop*taup) &
                        & + fr1_32 * (3._dp*t1*(2._dp+x1) - t2*(2._dp+x2)) * drhot**2 &
                        & - fr1_32 * (3._dp*t1*(2._dp*x1+1._dp) + t2*(2._dp*x2+1._dp)) * (drhon**2 + drhop**2) &
                        ! Extra BSk terms
                        & + fr14*t4 * ( (1._dp + fr12*x4) * (rhot*taut + fr34*drhot**2) &
                            & - (fr12+x4) * (rhon*taun+rhop*taup + fr34*(drhon**2+drhop**2)) ) * rhot**beta &
                        & + beta*fr18*t4 * ( (1._dp + fr12*x4) * rhot * drhot**2 &
                            & - (fr12+x4) * drhot * (rhon*drhon+rhop*drhop) ) * rhot**(beta-1._dp) &
                        & + fr14*t5 * ( (1._dp + fr12*x5) * (rhot*taut - fr14*drhot**2) &
                            & + (fr12+x5) * (rhon*taun+rhop*taup - fr14*(drhon**2+drhop**2)) ) * rhot**gamma
            end do radii_1
            
        else
            ! Standard Skyrme forces
            radii_2: do i=1, n
                rhot = rho_t(i)
                rhon = rho_q(i,0)
                rhop = rho_q(i,1)
                taut = tau_ETF_t(i)
                taun = tau_ETF_q(i,0)
                taup = tau_ETF_q(i,1)
                drhot = del_rho_t(i)
                drhon = del_rho_q(i,0)
                drhop = del_rho_q(i,1)
!                 Jt = J_t(i)
!                 Jp = J_q(i,0)
!                 Jn = J_q(i,1)
                
                ! Eq. (5) from Bartel and Bencheikh, EPJA (2002)
                e_density_field(i) = &
                        & + B1*rhot**2 + B2 * (rhon**2 + rhop**2) &
                        & + B3*rhot*taut + B4 * (rhon*taun + rhop*taup) &
                        & - B5*drhot**2 - B6 * (drhon**2 + drhop**2) &
                        & + rhot**sigma * (B7*rhot**2 + B8*(rhon**2+rhop**2) )
            end do radii_2
            
        end if
        
        ! Field energy (excluding kinetic and Coulomb)
        e_field = WS_integral(e_density_field(1:n),n)
        
        ! Add field and kinetic energy densities together to calculate Skyrme energy density
        e_density_sky(1:n) = e_density_field(1:n) + hbar2_2m_q(0)*tau_ETF_q(1:n,0) + hbar2_2m_q(1)*tau_ETF_q(1:n,1)
        
        ! Add contribution of spin-orbit
        e_density_sky(1:n) = e_density_sky(1:n) &
                            & - B9 * (J_q(1:n,0)*del_rho_q(1:n,0) + J_q(1:n,1)*del_rho_q(1:n,1) + J_t(1:n)*del_rho_t(1:n))
        
        ! J^2 terms
        if (J2_terms .eqv. .true.) then
            e_density_sky(1:n) = e_density_sky(1:n) &
                            & - fr1_16*(t1*x1 + t2*x2)*J_t(1:n)**2 + fr1_16*(t1-t2)*(J_q(1:n,0)**2 + J_q(1:n,1)**2)
        end if
        
!         e_density_sky(1:n) = e_density_field(1:n)
!         do q=0, 1
!             e_density_sky(1:n) = e_density_sky(1:n) &
!                                 & + f_q(1:n,q)*hbar2_2m_q(q) * (tau_TF_q(1:n,q)+tau_2_L_q(1:n,q)+tau_2_NL_q(1:n,q))
!             e_density_sky(1:n) = e_density_sky(1:n) - B9 * J_2_q(1:n,q)*del_rho_q(1:n,q)
!         end do
!         e_density_sky(1:n) = e_density_sky(1:n) - B9 * (J_2_q(1:n,0)+J_2_q(1:n,1)) * del_rho_t(1:n)
        
        
        ! Add Coulomb energy densities
        if (coulomb_on .eqv. .true.) then
            e_density_sky(1:n) = e_density_sky(1:n) + e_density_c_di(1:n) + e_density_c_ex(1:n)
        end if
        
        ! Skyrme energy (including Coulomb) of WS cell
        e_skyrme = WS_integral(e_density_sky(1:n),n)
        
        
!         e_kinetic_t = 0._dp
!         do q=0, 1
!             call lap_sphe_symm(A_q(-2:n,q),d1_A_q(1:n,q),div_W_q(1:n,q))
!             div_W_q(1:n,q) = div_W_q(1:n,q) * (-B9)
!             
!             cont1(1:n,q) = rho_q(1:n,q)**fr13 * ( &
!                         & 1._dp/270._dp*f_q(1:n,q)*(del2_rho_q(1:n,q)/rho_q(1:n,q))**2 &
!                         & -1._dp/240._dp*f_q(1:n,q)*del2_rho_q(1:n,q)/rho_q(1:n,q)*(del_rho_q(1:n,q)/rho_q(1:n,q))**2 &
!                         & +1._dp/810._dp*f_q(1:n,q)*(del_rho_q(1:n,q)/rho_q(1:n,q))**2 &
!                         & -1._dp/240._dp*del2_f_q(1:n,q)**2/f_q(1:n,q) &
!                         & +1._dp/120._dp*del2_f_q(1:n,q)*del_f_q(1:n,q)**2/f_q(1:n,q)**2 &
!                         & -1._dp/240._dp*del_f_q(1:n,q)**4/f_q(1:n,q)**3 &
!                         & +1._dp/360._dp*del2_f_q(1:n,q)*del_f_q(1:n,q)*del_rho_q(1:n,q)/(f_q(1:n,q)*rho_q(1:n,q)) &
!                         & -1._dp/360._dp*del2_rho_q(1:n,q)*del_f_q(1:n,q)**2/(f_q(1:n,q)*rho_q(1:n,q)) &
!                         & -7._dp/2160._dp*del2_f_q(1:n,q)*(del_rho_q(1:n,q)/rho_q(1:n,q))**2 &
!                         & +1._dp/540._dp*(del_rho_q(1:n,q)/rho_q(1:n,q))**2*del_f_q(1:n,q)**2/f_q(1:n,q) &
!                         & +7._dp/2160._dp*(del_f_q(1:n,q)*del_rho_q(1:n,q))**2/(f_q(1:n,q)*rho_q(1:n,q)**2) &
!                         & -11._dp/3240._dp*(del_rho_q(1:n,q)/rho_q(1:n,q))**2*del_f_q(1:n,q)*del_rho_q(1:n,q)/rho_q(1:n,q) &
!                         & +7._dp/1080._dp*del2_rho_q(1:n,q)*del_f_q(1:n,q)*del_rho_q(1:n,q)/rho_q(1:n,q)**2 &
!                         & +1._dp/180._dp*del2_f_q(1:n,q)*del2_rho_q(1:n,q)/rho_q(1:n,q) )
!             e_kinetic_t = e_kinetic_t + (3._dp*pi**2)**(-fr23) * WS_integral(cont1(1:n,q),n)
!             
!             cont2(1:n,q) = rho_q(1:n,q)**fr13/f_q(1:n,q) * ( &
!                         & fr14*div_W_q(1:n,q)**2 &
!                         & -fr38*div_W_q(1:n,q)*del_f_q(1:n,q)*W_q(1:n,q)/f_q(1:n,q) &
!                         & +fr1_16*W_q(1:n,q)**2*del2_f_q(1:n,q)/f_q(1:n,q) &
!                         & +fr18*(del_f_q(1:n,q)*W_q(1:n,q))**2/f_q(1:n,q)**2 &
!                         & +fr1_24*div_W_q(1:n,q)*del_rho_q(1:n,q)*W_q(1:n,q)/rho_q(1:n,q) &
!                         & +1._dp/48._dp*W_q(1:n,q)**2*del2_rho_q(1:n,q)/rho_q(1:n,q) &
!                         & -1._dp/72._dp*W_q(1:n,q)**2*(del_rho_q(1:n,q)/rho_q(1:n,q))**2 &
!                         & +fr18/hbar2_2m_q(q)**2*W_q(1:n,q)**4/f_q(1:n,q)**2 )
!             e_kinetic_t = e_kinetic_t + (3._dp*pi**2)**(-fr23) * fr12 / hbar2_2m_q(q) * WS_integral(cont2(1:n,q),n)
!         end do
!         e_skyrme = e_skyrme + e_kinetic_t
        
    end subroutine eval_skyrme_energy_density
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the nuclear contribution to the pressure, and
    !> then the total pressure.
    subroutine pressure()
        implicit none
        real(kind=dp) :: tau_edge !< Kinetic density of neutrons at edge of box
        real(kind=dp) :: rho_n_edge !< Density of neutrons at edge of box
        real(kind=dp) :: C0n, C1n, C0tau, C1tau, dC0n_drho, dC1n_drho, dC0tau_drho, dC1tau_drho !< Constants
        
        ! Value of neutron gas at edge of box
        rho_n_edge = n_profile(1)
        
        ! Kinetic density of neutrons at edge of box
        tau_edge = fr35 * (3._dp*pi**2)**fr23 * rho_n_edge**fr53
        
        ! Constants
        if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 18) then
            C0n = fr38*t0 + fr1_16*t3 * rho_n_edge**alpha
            C1n = -fr14*t0*(fr12+x0) - fr1_24*t3*(fr12+x3) * rho_n_edge**alpha
            
            C0tau = 3._dp/16._dp*t1 + fr14*t2*(fr54+x2) &
                & + 3._dp/16._dp*t4 * rho_n_edge**beta + fr14*t5*(fr54+x5) * rho_n_edge**gamma
            C1tau = -fr18*t1*(fr12+x1) + fr18*t2*(fr12+x2) &
                & - fr18*t4*(fr12+x4) * rho_n_edge**beta + fr18*t5*(fr12+x5) * rho_n_edge**gamma
            
            dC0n_drho = fr1_16*t3*alpha * rho_n_edge**(alpha-1._dp)
            dC1n_drho = -fr1_24*t3*(fr12+x3)*alpha * rho_n_edge**(alpha-1._dp)
            
            dC0tau_drho = 3._dp/16._dp*t4*beta * rho_n_edge**(beta-1._dp) &
                        & + fr14*t5*(fr54+x5)*gamma * rho_n_edge**(gamma-1._dp)
            dC1tau_drho = -fr18*t4*(fr12+x4)*beta * rho_n_edge**(beta-1._dp) &
                        & + fr18*t5*(fr12+x5)*gamma * rho_n_edge**(gamma-1._dp)
            
        else
            C0n = fr38*t0 + fr1_16*t3 * rho_n_edge**sigma
            C1n = -fr14*t0*(fr12+x0) - fr1_24*t3*(fr12+x3) * rho_n_edge**sigma
            C0tau = 3._dp/16._dp*t1 + fr14*t2*(fr54+x2)
            C1tau = -fr18*t1*(fr12+x1) + fr18*t2*(fr12+x2)
            dC0n_drho = fr1_16*t3*sigma * rho_n_edge**(sigma-1._dp)
            dC1n_drho = -fr1_24*t3*(fr12+x3)*sigma * rho_n_edge**(sigma-1._dp)
            dC0tau_drho = 0._dp
            dC1tau_drho = 0._dp
            
        end if
        
        ! Nuclear pressure (Eq. B28 from PRC 85, 065803 (2012))
            ! NB A mistake is printed for coefficient C1n (Eq. B30b)
            ! (1 + x3) should be (1/2 + x3)
        pressure_nucl = fr13*sum(hbar2_2m_q(:))*tau_edge &
                    & + rho_n_edge**2 * (C0n + C1n) &
                    & + fr53 * rho_n_edge * tau_edge * (C0tau + C1tau) &
                    & + rho_n_edge**3 * (dC0n_drho + dC1n_drho) &
                    & + rho_n_edge**2 * tau_edge * (dC0tau_drho + dC1tau_drho)
        
        ! Total pressure
        pressure_t = pressure_nucl + pressure_e + pressure_ex
        
    end subroutine
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the number of particles in the Wigner-Seitz
    !> cell
    subroutine calc_particle_number()
        implicit none
        integer :: q !< Loop index
        
        isospin: do q=0, 1
            N_q(q) = WS_integral(rho_q(1:n,q),n)
        end do isospin
        
        N_t = sum(N_q)
        
    end subroutine calc_particle_number
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the various densities, fields, derivatives,
    !> energies, and particle numbers in the WS cell with densities
    !> \f$rho_q\f$
    subroutine eval_WS_quantities()
        implicit none
        integer :: i !< Loop index
        
        ! Evaluate densities using external profile parameters
        densities: do i=1, n
            rho_q(i,0) = calc_rho_q(n_profile(1),n_profile(2),n_profile(3),n_profile(4),n_profile(5),r(i))
            rho_q(i,1) = calc_rho_q(p_profile(1),p_profile(2),p_profile(3),p_profile(4),p_profile(5),r(i))
            
            ! Total density
            rho_t(i) = rho_q(i,0) + rho_q(i,1)
        end do densities
        
        ! Initialise force
        call force_initialise()
        
        ! Calculate particle numbers
        call calc_particle_number()
        
        ! Evaluate all fields using matter densities
        call eval_fields()
        
        ! Calculate kinetic energy densities
        call eval_tau_ETF()
        
        ! Coulomb fields, needed for central fields
        if (coulomb_on .eqv. .true.) then
            ! BSk15 and onwards
            if ((force_string(1:3) .eq. 'BSk') .and. str2int(force_string(4:5)) .ge. 15) then
                ! Calculate charge density
                call charge(rho_q(1:n,:),N_q(0),N_q(1),n,dr,rho_ch)
                
                ! Use charge density, set Coulomb exchange to 0
                call eval_coulomb(rho_ch,.false.)
                
            ! Standard Skyrme forces
            else
                ! Use proton density, calculate Coulomb exchange
                call eval_coulomb(rho_q(1:n,1),.true.)
                
            end if
        end if
        
        ! Energy contributions of homogeneous electron gas
        if (electrons_on .eqv. .true.) then
            call eval_electrons()
        end if
        
        ! Evaluate central fields, adding Coulomb contribution for protons
        call eval_U_q()
        
        ! Evalute energy density and calculate Skyrme energy
        call eval_skyrme_energy_density()
        
        ! Total energy is Skyrme energy + electron contributions
        e_total = e_skyrme + e_kinetic_e + e_coulomb_e + e_coulomb_pe
        
        ! Pressure
        call pressure()
        
    end subroutine eval_WS_quantities
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate the Fermi energy at a given density and effective
    !> mass, for neutrons or protons
    !> @param[in]   rhoq    Neutron or proton density
    !> @param[in]   fn      Neutron or proton effective mass
    function calc_e_f_q(rhoq,fq,q)
        implicit none
        real(kind=dp) :: calc_e_f_q
        real(kind=dp),  intent(in)  :: rhoq,fq
        integer,        intent(in)  :: q
        
        calc_e_f_q = hbar2_2m_q(q) * fq * (3._dp * pi**2 * rhoq)**fr23
        
    end function calc_e_f_q
    
    
    
    !> @author A. Pastore, M. Shelley
    !> @brief
    !> This subroutine computes the first derivative of function evaluated on
    !> the meshpoints 1,...,npt. The input is the function f with extrapolated
    !> values in -1, 0. Modified by M. Shelley.
    !> @param[in]       n   Number of meshpoints
    !> @param[in]       f   Input function (1:n) with extrapolated values in
    !>                      -1,0
    !> @param[in]       h   Step size
    !> @param[in,out]   df  Array for first derivatives
    subroutine d1(n,h,f,df)
        implicit none
        integer,        intent(in)      :: n
        real(kind=dp),  intent(in)      :: f(-1:n),h
        real(kind=dp),  intent(inout)   :: df(1:n)
        real(kind=dp) :: h_12
        integer :: i
        
        df=0._dp
        
        h_12 = 12._dp * h
        
        do i = 1, n - 2
        df(i) = ( 8._dp * ( f(i+1) - f(i-1) ) - f(i+2) + f(i-2) ) / h_12
        end do
        
        df(n-1) = ( - f(n-4) + 6._dp * f(n-3) - 18._dp * f(n-2) &
            + 10._dp * f(n-1) + 3._dp * f(n) ) / h_12
        df(n) = ( 3._dp * f(n-4) - 16._dp * f(n-3) + 36._dp * f(n-2) &
            - 48._dp * f(n-1) + 25._dp * f(n) ) / h_12
        
        return ! Exit
        
    end subroutine d1
    
    
    !> @author A. Pastore, M. Shelley
    !> @brief
    !> This subroutine computes the second derivative of function evaluated on
    !> the meshpoints 1,...,npt. The input is the function f with extrapolated
    !> values in -1, 0. Modified by M. Shelley.
    !> @param[in]       n   Number of meshpoints
    !> @param[in]       f   Input function (1:n) with extrapolated values in
    !>                      -1,0
    !> @param[in]       h   Step size
    !> @param[in,out]   d2f Array for second derivatives
    subroutine d2(n,h,f,d2f)
        implicit none
        integer,        intent(in)      :: n
        real(kind=dp),  intent(in)      :: f(-1:n),h
        real(kind=dp),  intent(inout)   :: d2f(1:n)
        real(kind=dp) :: hh_12
        integer :: i
        
        d2f=0._dp
        
        hh_12 = 12._dp * h * h
        
        do i = 1, n - 2
        d2f(i) = ( - 30._dp * f(i) + 16._dp * ( f(i+1) + f(i-1) ) &
                - f(i+2) - f(i-2) ) / hh_12
        end do
        
        d2f(n-1) = ( - f(n-4) + 4._dp * f(n-3) + 6._dp * f(n-2) &
            - 20._dp * f(n-1) + 11._dp * f(n) ) / hh_12
        d2f(n) = ( 11._dp * f(n-4) - 56._dp * f(n-3) + 114._dp * f(n-2) &
            - 104._dp * f(n-1) + 35._dp * f(n) ) / hh_12
        
        return ! Exit
        
    end subroutine d2
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to carry out Laplacian in spherical symmetry.
    !> @param[in]       f       Input function (1:n) with extrapolated values
    !>                          in -2,-1,0
    !> @param[in]       df_dr   First derivative array
    !> @param[in,out]   d2f_dr2 Array for Laplacian
    subroutine lap_sphe_symm(f,df_dr,d2f_dr2)
        implicit none
        real(kind=dp),  dimension(-2:n), intent(in)     :: f
        real(kind=dp),  dimension(1:n),  intent(in)     :: df_dr
        real(kind=dp),  dimension(1:n),  intent(inout)  :: d2f_dr2
!         real(kind=dp),  dimension(1:n)                  :: deriv_store
        
        ! Calculate second derivative of f
        call d2(n,dr,f(-1:n),d2f_dr2(:))
        
        ! Calculate (d2f/dr2) + 2/r*(df/dr)
        d2f_dr2(1:n) = d2f_dr2(1:n) + 2._dp/r(1:n)*df_dr(1:n)
        
        
!         call d1(n,dr,(r(:)**2) * df_dr,deriv_store)
!         d2f_dr2(:) = deriv_store / r(:)**2
        
        
!         call d2(n,dr,r(:)*f,d2f_dr2)
!         d2f_dr2(:) = d2f_dr2(:) / r(:)
        
        
!         call extrapolate_back_3(d2f_dr2(:))
        
    end subroutine lap_sphe_symm
    
    
!     !> @author M. Shelley
!     !> @brief
!     !> Subroutine to calculate divergence in spherical symmetry.
!     !> @param[in]       f       Input function (1:n) with extrapolated values
!     !>                          in -2,-1,0
!     !> @param[in]       df_dr   First derivative array
!     !> @param[in,out]   divf    Array for divergence
!     subroutine div_sphe_symm(f,df_dr,divf)
!         implicit none
!         real(kind=dp),  dimension(-2:n), intent(in)     :: f
!         real(kind=dp),  dimension(1:n),  intent(in)     :: df_dr
!         real(kind=dp),  dimension(1:n),  intent(inout)  :: divf
!         real(kind=dp),  dimension(-2:n)                 :: deriv_store
!         
!         ! Calculate 2*f/r + (df/dr)
! !         divf(1:n) = 2._dp*f(1:n)/r(1:n) + df_dr(1:n)
!         
!         ! Calculate (1/r**2) * d/dr(f*r**2)
!         deriv_store(1:n) = f(1:n) * r(1:n)**2
!         call extrapolate_back_3(deriv_store)
!         call d1(n,dr,deriv_store(-1:n),divf(:))
!         divf(:) = divf(:) / r(1:n)**2
!         
!     end subroutine div_sphe_symm
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate third derivative of function f evaluated on
    !> evenly-spaced meshpoints from 1 to n. Works with extrapolated values in
    !> -2,-1,0. Uses 7-point stencil.
    !> @param[in]       f       Input function (1:n) with extrapolated values
    !>                          in -2,-1,0
    !> @param[in,out]   d3f_dr3 Array for third derivative
    subroutine d3(f,d3f_dr3)
        implicit none
        real(kind=dp),  dimension(-2:n), intent(in)     :: f
        real(kind=dp),  dimension(1:n),  intent(inout)  :: d3f_dr3
        real(kind=dp) :: h3m8 !< Denominator for derivative evaluations
        integer :: i !< Loop index
        
        h3m8 = 8._dp * dr**3
        
        ! Central difference for i=1 to n-3
        radii: do i=1, n-3
            d3f_dr3(i) = ( f(i-3) - f(i+3) &
                        & - 8._dp * (f(i-2) - f(i+2)) &
                        & + 13._dp * (f(i-1) - f(i+1)) ) / h3m8
        end do radii
        
        ! Non-central differences for cell boundary
        d3f_dr3(n-2) = ( - f(i-4) + f(i+2) &
                        & + 8._dp * (f(i-3) + f(i+1)) &
                        & - 29._dp*f(i-2) + 48._dp*f(i-1) - 35._dp*f(i) ) / h3m8
        
        d3f_dr3(n-1) = ( f(i-5) - 8._dp*f(i-4) + 29._dp*f(i-3) - 64._dp*f(i-2) &
                    & + 83._dp*f(i-1) - 56._dp*f(i) + 15._dp*f(i+1) ) / h3m8
        
        d3f_dr3(n) = ( 15._dp*f(i-6) - 104._dp*f(i-5) + 307._dp*f(i-4) - 496._dp*f(i-3) &
                    & + 461._dp*f(i-2) - 232._dp*f(i-1) + 49._dp*f(i) ) / h3m8
        
    end subroutine d3
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate fourth derivative of function f evaluated on
    !> evenly-spaced meshpoints from 1 to n. Works with extrapolated values in
    !> -2,-1,0. Uses 7-point stencil.
    !> @param[in]       f       Input function (1:n) with extrapolated values
    !>                          in -2,-1,0
    !> @param[in,out]   d4f_dr4 Array for fourth derivative
    subroutine d4(f,d4f_dr4)
        implicit none
        real(kind=dp),  dimension(-2:n), intent(in)     :: f
        real(kind=dp),  dimension(1:n),  intent(inout)  :: d4f_dr4
        real(kind=dp) :: h4m6 !< Denominator for derivative evaluations
        integer :: i !< Loop index
        
        h4m6 = 6._dp * dr**4
        
        ! Central difference for i=1 to n-3
        radii: do i=1, n-3
            d4f_dr4(i) = ( - f(i-3) - f(i+3) &
                        & + 12._dp * (f(i-2) + f(i+2)) &
                        & - 39._dp * (f(i-1) + f(i+1)) &
                        & + 56._dp * f(i) ) / h4m6
        end do radii
        
        ! Non-central differences for cell boundary
        d4f_dr4(n-2) = ( - f(i-4) + 6._dp*f(i-3) - 9._dp*f(i-2) - 4._dp*f(i-1) &
                    & + 21._dp*f(i) - 18._dp*f(i+1) + 5._dp*f(i+2) ) / h4m6
        
        d4f_dr4(n-1) = ( 5._dp*f(i-5) - 36._dp*f(i-4) + 111._dp*f(i-3) - 184._dp*f(i-2) &
                    & + 171._dp*f(i-1) - 84._dp*f(i) + 17._dp*f(i+1) ) / h4m6
        
        d4f_dr4(n) = ( 17._dp*f(i-6) - 114._dp*f(i-5) + 321._dp*f(i-4) - 484._dp*f(i-3) &
                    & + 411._dp*f(i-2) - 186._dp*f(i-1) + 35._dp*f(i) ) / h4m6
        
    end subroutine d4
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to extrapolate back 3 points on the r mesh, to faciliate the
    !> calculation of first and second derivatives using a 5-point stencil.
    !> Values in grid start at 4, extrapolated values are put in elements
    !> 3,2,1. Coefficients come from solving 4th-order polynomial, assuming
    !> that f'(0) = 0, and f(-h) = f(h).
    !> @param[in,out]   grid    Array of quantity evaluated on r mesh, with 3
    !>                          elements before r = dr
    subroutine extrapolate_back_3(grid)
        implicit none
        real(kind=dp), dimension(:), intent(inout) :: grid
        
!         grid(1) = grid(5)
!         grid(2) = grid(4)
!         
!         grid(3) = (15._dp*grid(4) - 6._dp*grid(5) + grid(6)) / 10._dp
        
        
!         grid(3) = 2._dp*grid(4) - grid(5)
!         grid(2) = 2._dp*grid(3) - grid(4)
!         grid(1) = 2._dp*grid(2) - grid(3)
        
        
!         grid(3) = 3._dp*grid(4) - 3._dp*grid(5) + grid(6)
!         grid(2) = 3._dp*grid(3) - 3._dp*grid(4) + grid(5)
!         grid(1) = 3._dp*grid(2) - 3._dp*grid(3) + grid(4)
        
        
        grid(3) = 5._dp*grid(4) - 10._dp*grid(5) + 10._dp*grid(6) - 5._dp*grid(7) + grid(8)
        grid(2) = 5._dp*grid(3) - 10._dp*grid(4) + 10._dp*grid(5) - 5._dp*grid(6) + grid(7)
        grid(1) = 5._dp*grid(2) - 10._dp*grid(3) + 10._dp*grid(4) - 5._dp*grid(5) + grid(6)
        
    end subroutine extrapolate_back_3
    
    
    !> @author M. Shelley
    !> @brief
    !> Function to integrate a density or field over the whole W-S cell, using
    !> Simpson's rule (+ Simpson's 3/8 rule if odd number of points)
    !> @param[in] quantity  Array of quantity, evaluated on r mesh, to be
    !>                      integrated
    !> @param[in] n_max     Maximum mesh point in array for integration
    function WS_integral(quantity,n_max)
        implicit none
        real(kind=dp) :: WS_integral
        real(kind=dp), dimension(:),    intent(in)  :: quantity
        integer,                        intent(in)  :: n_max
        integer :: i !< Loop index
        real(kind=dp) :: sum_1 !< First sum
        real(kind=dp) :: sum_2 !< Second sum
        
        WS_integral = 0._dp
        
!         ! Euler method
!         radii: do i=1, size(quantity)
!             WS_integral = WS_integral + quantity((i)) * r(i)**2
!         end do radii
!         
!         WS_integral = WS_integral * 4._dp*pi * dr
        
        
        !!! Simpson's rule !!!
        ! Value at first point (r=0) always 0, so not included in sums
        
        ! Composite Simpson's rule requires even number of partitions
        if (mod(n_max,2) .eq. 0) then
            ! First sum
            sum_1 = 0._dp
            loop_1: do i=1, n_max-1, 2
                sum_1 = sum_1 + quantity(i) * r(i)**2
            end do loop_1
            
            ! Second sum
            sum_2 = 0._dp
            loop_2: do i=2, n_max-2, 2
                sum_2 = sum_2 + quantity(i) * r(i)**2
            end do loop_2
            
            ! h/3 factor for Simpson's rule, 4*pi for spherical integral
            WS_integral = dr/3._dp * 4._dp*pi * (4._dp*sum_1 + 2._dp*sum_2 + quantity(n_max) * r(n_max)**2)
            
        ! If n_max is odd, use Simpson's rule up to n_max-3, then Simpson's 3/8 rule
        ! for final 3 partitions (4 mesh points)
        else
            ! First sum
            sum_1 = 0._dp
            loop_1a: do i=1, n_max-3, 2
                sum_1 = sum_1 + quantity(i) * r(i)**2
            end do loop_1a
            
            ! Second sum
            sum_2 = 0._dp
            loop_2a: do i=2, n_max-4, 2
                sum_2 = sum_2 + quantity(i) * r(i)**2
            end do loop_2a
            
            ! h/3 factor for Simpson's rule, 4*pi for spherical integral
            WS_integral = dr/3._dp * 4._dp*pi * (4._dp*sum_1 + 2._dp*sum_2 + quantity(n_max) * r(n_max)**2)
            
            ! Final 3 partitions with Simpson's rule
            ! 3h/8 factor for 3/8 rule, 4*pi for spherical integral
            WS_integral = WS_integral + 3._dp*dr/8._dp * 4._dp*pi * &
                                        & ( quantity(n_max-3)*r(n_max-3)**2 + 3._dp*quantity(n_max-2)*r(n_max-2)**2 &
                                        & + 3._dp*quantity(n_max-1)*r(n_max-1)**2 + quantity(n_max)*r(n_max)**2 )
            
        end if
        
    end function WS_integral
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write neutron and proton densities to files
    subroutine write_densities()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        file_unit_0 = 200
        file_unit_1 = 201
        
        open(unit=file_unit_0,file='densities_n.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        open(unit=file_unit_1,file='densities_p.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        densities: do i=1, n
            write(unit=file_unit_0,fmt=n_floats(4)) r(i),rho_q(i,0),tau_ETF_q(i,0),J_q(i,0)
            write(unit=file_unit_1,fmt=n_floats(4)) r(i),rho_q(i,1),tau_ETF_q(i,1),J_q(i,1)
        end do densities
        
        close(unit=file_unit_0,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        close(unit=file_unit_1,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
    end subroutine write_densities
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write neutron and proton fields to files
    subroutine write_fields()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        file_unit_0 = 202
        file_unit_1 = 203
        
        open(unit=file_unit_0,file='fields_n.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        open(unit=file_unit_1,file='fields_p.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        densities: do i=1, n
            write(unit=file_unit_0,fmt=n_floats(6)) r(i),U_q(i,0),f_q(i,0)*hbar2_2m_q(0),W_q(i,0), &
                                                    & pair_gap_n(i),delta_n(i)
            write(unit=file_unit_1,fmt=n_floats(9)) r(i),U_q(i,1),f_q(i,1)*hbar2_2m_q(1),W_q(i,1), &
                                                    & rho_p_bcs(i),rho_anom_p_bcs(i),strength_bcs(i),delta_p(i)
        end do densities
        
        close(unit=file_unit_0,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        close(unit=file_unit_1,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
    end subroutine write_fields
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write proton single particle energies to file, and extra
    !> details if BCS has been performed
    subroutine write_sp_states_p()
        implicit none
        character(1), dimension(0:10) :: l !< Orbital labels
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index

        !< Orbital labels
        l = (/'s','p','d','f','g','h','i','j','k','l','m'/)
        
        open(unit=204,file='sp_states_p.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! If no BCS carried out, only write L, J, and E
        if ((proton_bcs .eqv. .false.) .and. (calc_chem_pots .eqv. .false.)) then
            ! Column headings
            write(unit=204,fmt='(1x,a5,2x,a12)') 'State','Energy (MeV)'
            
            states_1: do i=1, NNST
                ! Write sorted states to file
                write(unit=204,fmt='(a1,1x,i2,a2,2x,f11.6)') l(LP(i)),JJP(i),'/2',EP(i)
            end do states_1
            
        ! Write all BCS details
        else
            ! Column headings
            write(unit=204,fmt='(1x,a5,1x,a10,1x,a10,1x,a8,1x,a8,1x,a9)') &
                    & 'State','E (MeV)','E_QP (MeV)','V^2','U^2','Gap (MeV)'
            
            states_2: do i=1, NNST
                ! Write all details to file
                write(unit=204,fmt='(a1,1x,i2,a2,1x,f10.6,1x,f10.6,1x,f8.6,1x,f8.6,1x,f9.6)') &
                    & l(LP(i)),JJP(i),'/2',EP(i),e_qp_p(i),occ_v(i)**2,occ_u(i)**2,pair_gap_p(i)
            end do states_2
            
        end if
        
        close(unit=204,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
    end subroutine write_sp_states_p
    
    
end module routines