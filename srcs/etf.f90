!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Program to calculate the equation of state in the inner crust with the 
!> extended Thomas-Fermi (ETF) approach.
program etf
    
    use parameters
    use routines
    use rspace
    use strutinsky
    use pairing
    use test
    implicit none
    
    integer :: i !< Loop index
    integer :: ierr !< I/O error storage
    
    character(len=24),  dimension(11)       :: args !< Store for command line arguments read in as strings
    real(kind=dp),      dimension(11)       :: args_real !< Store for real versions of cmd-line arguments
    
    real(kind=dp),      dimension(200,14)   :: old_mins_pars !< Store for reading in minima parameters from old inner crust runs
    
    ! Read input parameters from 'input.in'
    open(unit=10,file='input.in',status='old',iostat=ierr)
    if (ierr.ne.0) stop
    read(10,nml=params)
    close(unit=10,iostat=ierr,status='keep')
    if (ierr.ne.0) stop
    
    if (run_mode .ne. 1) then
        print '(a,2(/a))', line_break,'ETF program start',line_break
    end if
    
    select case(run_mode)
        ! Minimisation mode: minimal output
        case(1)
            ! Read all 10 profile parameters, and r_ws, from command line
            do i = 1, 11
                ! Store each parameter in "args"
                call get_command_argument(i,args(i))
                
                ! Each parameter is padded with whitespace up to 24 characters
                read(args(i),'(f24.21)') args_real(i)
            end do
            
            ! Store command line arguments as profile parameters
            do i = 1, 5
                n_profile(i) = args_real(i)
                p_profile(i) = args_real(i+5)
            end do
            
            ! Radius of WS cell
            r_ws = args_real(11)
            
            ! "verbose" must be false, so that only energy and delta_mu are printed
            verbose = .false.
            
            ! Allocate memory
            call initialise_uniform_mesh()
            call allocate_arrays()
            
            ! All densities, fields, energies, and quantum corrections
            call total_energy_with_corrections(.false.)
            
            ! Print outputs for minimisation
            print '(es24.16e3,1x,es24.16e3)', e_total/N_t, delta_mu
            
            
        ! Test mode
        case(2)
            ! Print input parameters
            write (*,nml=params)
            
            print '(3(/a))', line_break,'Running in test mode',line_break
            
            ! Read parameters for testing 
            call read_test_params()
            
            ! Don't read WS cell radius from profile.in
            profile_r_max = .false.
            ! Use full ETF for test
            etf_order = 4
            ! Coulomb is on for test
            coulomb_on = .true.
            ! Electrons are off for test
            electrons_on = .false.
            ! Use complete arrays for Strutinsky test
            strut_r_max = r_max
            ! Do not modify emax0 in boundary
            emax0_mod = .false.
            
            ! Initialise force
            call force_initialise()
            
            ! Allocate memory
            call initialise_uniform_mesh()
            call allocate_arrays()
            
            ! Calculate densities and write to 'densities.dat'
            call test_densities()
            
            ! Calculate particle numbers
            call calc_particle_number()
            
            ! Evaluate all fields using matter densities
            call eval_fields()
            
            ! Calculate kinetic energy densities
            call eval_tau_ETF()
            
            ! Coulomb fields, needed for central fields
            call eval_coulomb(rho_q(1:n,1),.true.)
            
            ! Evaluate central fields, adding Coulomb contribution for protons
            call eval_U_q()
            
            ! Evalute energy density and calculate Skyrme energy
            call eval_skyrme_energy_density()
            
            ! Test particle number routines
            call test_calc_particle_number()
            
            ! Calculate effective mass at test values and write to 'eff_mass.dat'
            call test_eff_mass()
            
            ! Write density derivatives to 'derivatives.dat'
            call test_density_derivs()
            
            ! Write kinetic energy densities to file
            call test_kinetic_densities()
            
            ! Write spin current densities to file
            call test_spin_current_densities()
            
            ! Write order-by-order kinetic densities
            call orders_kinetic_densities()
            
            ! Write central potentials to file
            call test_central_potentials()
            
            ! Test Skyrme energy (including Coulomb)
            call test_calc_skyrme_energy()
            
            ! Test Coulomb energy
            call test_calc_coulomb_energy()
            
            ! Test Strutinsky integral correction routines for Pb208
            if (sample_profile .eq. 2 .or. sample_profile .eq. 3 .or. sample_profile .eq. 4 .or. sample_profile .eq. 5) then
                Nmaxstate = 200
                call sp_states_setup()
                
                Ecut(0) = -0._dp !< Cutoff energy for neutrons
                Ecut(1) = maxval(U_q(:,1)) ! Maximum of proton potential (including Coulomb)
                
                Lmax = 9 !< Maximum L for Pb208
                call test_sp_energies()
                
                ! Correct particle numbers for test
                N_q(0) = 126._dp
                N_q(1) = 82._dp
                
                ! Strutinsky on for neutrons and protons
                strutinsky_on(:) = .true.
                
                call test_si_correction()
                
                call sp_states_deallocate()
            end if
            
            ! Benchmarks
            call bartel_bencheikh_benchmark()
            
            
        ! Read external density profiles and particle numbers; calculate WS cell densities, fields, energies
        case(3)
            ! Print input parameters
            write (*,nml=params)
            
            print '(3(/a))', line_break,'Running in mode to read external profiles',line_break
            
            ! Read profile parameters from 'profile.in'
            open(unit=11,file='profile.in',status='old',iostat=ierr)
            if (ierr.ne.0) stop
            read(11,nml=profiles)
            close(unit=11,iostat=ierr,status='keep')
            if (ierr.ne.0) stop
            
            ! Print parameters to screen
            write (*,nml=profiles)
            
            ! Allocate memory
            call initialise_uniform_mesh()
            call allocate_arrays()
            
            ! All densities, fields, energies, and quantum corrections
            call total_energy_with_corrections(.true.)
            
            
        ! Read profile parameters minimum configuration file; calculate pressure
        case(4)
            ! Allocate memory
            call initialise_uniform_mesh()
            call allocate_arrays()
            
            ! Initialise force
            call force_initialise()
            
            ! File containing parameters of minima
            open(unit=20,file='mins_pars.dat',status='old',iostat=ierr)
            if (ierr.ne.0) stop
            
            ! File for writing pressure data
            open(unit=21,file='pressure.dat',status='unknown',iostat=ierr)
            if (ierr.ne.0) stop
            
            do i=1, 200
                ! Read parameters from file
                read(20,*) old_mins_pars(i,:)
                
                ! Z, R_WS, rho_n_gas
                N_q(1) = old_mins_pars(i,2)
                r_max = old_mins_pars(i,4)
                n_profile(1) = old_mins_pars(i,5)
                
                ! "Normal" and Coulomb exchange pressure
                call eval_electrons()
                
                ! Nuclear and total pressure
                call pressure()
                
                ! Total pressure
                write(unit=21,fmt=n_floats(5)) old_mins_pars(i,1), pressure_nucl, pressure_e, pressure_ex, pressure_t
                
            end do
            
            ! Close files
            close(unit=20,iostat=ierr,status='keep')
            if (ierr.ne.0) stop
            close(unit=21,iostat=ierr,status='keep')
            if (ierr.ne.0) stop
            
            
        ! If no mode selected, exit
        case default
            print '(a)', 'STOP: no running mode selected'
            stop
    end select
    
    ! Deallocate memory
    call deallocate_mesh()
    call deallocate_arrays()
    
    if (run_mode .ne. 1) then
        print '(3(/a))', line_break,'Program end',line_break
    end if
    
    
    contains
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to evaluate all WS quantities using supplied profile
    !> parameters, and calculate total energy, with pairing and shell
    !> corrections included. Write densities and fields to files, and all info
    !> to screen, if specified.
    !> @param[in]   write_files_print_info  Whether to write to files and screen
    subroutine total_energy_with_corrections(write_files_print_info)
        implicit none
        logical, intent(in) :: write_files_print_info
        
        ! All calculations of densities, fields, and energies
        call eval_WS_quantities()
        
        ! Pairing correction for neutrons
        if (neutron_pairing .ne. 0) then
            call calc_e_pair()
            
            ! Add pairing condensation energy to total energy
            e_total = e_total + e_pair_q(0)
        end if
        
        ! Allocate arrays for Strutinsky and pairing
        call sp_states_setup()
        
        ! "boundary" parameters
!         Ecut(1) = maxval(U_q(:,1)) ! Maximum of proton potential (including Coulomb)
        Ecut(1) = 0._dp
        
        ! Pairing correction for protons
        if ((proton_bcs .eqv. .true.) .or. (calc_chem_pots .eqv. .true.)) then
            ! Perform BCS calculation
            call bcs_protons()
            
            ! Add pairing energy to total energy
            e_total = e_total + e_pair_bcs
            
            if (write_files_print_info .eqv. .true.) then
                ! Write proton states to file
                call write_sp_states_p()
            end if
        end if
        
        ! Strutinsky correction for protons
        if (strutinsky_on(1) .eqv. .true.) then
            if (proton_bcs .eqv. .false.) then
                ! Standard Strutinsky energy correction
                call calc_e_sc(1)
                
                ! Write proton states to file
                call write_sp_states_p()
                
            else
                ! Strutinsky energy correction using BCS particle occupations
                call calc_e_sc_pair()
                
            end if
            
            ! Total correction energy
            e_sc_t = sum(e_sc_q)
            
            ! Total WS cell energy
            e_total = e_total + e_sc_t
        end if
        
        if (write_files_print_info .eqv. .true.) then
            ! Write all densities and fields to files
            call write_densities()
            call write_fields()
        end if
        
        ! Deallocate arrays for Strutinsky and pairing
        call sp_states_deallocate()
        
        if (calc_chem_pots .eqv. .true.) then
            ! Calculate neutron chemical potential, using values at edge of
            ! cell, for density, mass, and central potential
            mu_q(0) = calc_e_f_q(rho_q(n,0),f_q(n,0),0) + U_q(n,0)

            ! Overall chemical potential difference
            delta_mu = mec2 + mu_e + mpc2 + mu_q(1) - mnc2 - mu_q(0)
        end if
        
        ! Incorporate neutron-proton mass difference into total energy
        e_total = e_total - N_q(1)*0.782343_dp
        
        if (write_files_print_info .eqv. .true.) then
            ! Print all energies
            print '(3(/a))', line_break,'Final energies',line_break
            print '(a54)', '============================+========================='
            print '(8x,a12,8x,a1,6x,a12,7x)', 'Contribution', '|', 'Energy [MeV]'
            print '(a54)', '============================+========================='
            print info_format, 'Neutron kinetic', '|', e_kinetic_q(0)
            print info_format, 'Proton kinetic', '|', e_kinetic_q(1)
            print info_format, 'Total kinetic', '|', e_kinetic_t
            print '(a54)', table_string
            print info_format, 'Field', '|', e_field
            print '(a54)', table_string
            print info_format, 'Spin-orbit', '|', e_so
            print '(a54)', table_string
            print info_format, 'Direct Coulomb', '|', e_coulomb_di
            print info_format, 'Exchange Coulomb', '|', e_coulomb_ex
            print info_format, 'Total Coulomb', '|', e_coulomb
            print '(a54)', table_string
            print info_format, 'Total Skyrme', '|', e_skyrme
            print '(a54)', table_string
            print info_format, 'Proton shell-correction', '|', e_sc_q(1)
            print '(a54)', table_string
            print info_format, 'Neutron pairing-effect', '|', e_pair_q(0)
            print info_format, 'Proton pairing energy', '|', e_pair_bcs
            print '(a54)', table_string
            print info_format, 'Electron kinetic', '|', e_kinetic_e
            print info_format, 'Electron-electron Coulomb', '|', e_coulomb_e
            print info_format, 'Proton-electron Coulomb', '|', e_coulomb_pe
            print '(a54)', table_string
            print info_format, 'Neutron chemical potential', '|', mu_q(0)
            print info_format, 'Proton chemical potential', '|', mu_q(1)
            print info_format, 'Electron chemical potential', '|', mu_e
            print info_format, 'Overall chemical potential', '|', delta_mu
            print '(a54)', table_string
            print info_format, 'Corrected total energy', '|', e_total
            print info_format, 'Total E/A', '|', e_total/N_t
            print '(a54)', table_string
            print info_format, 'Nuclear pressure', '|', pressure_nucl
            print info_format, 'Electron pressure', '|', pressure_e
            print info_format, 'Exchange Coulomb pressure', '|', pressure_ex
            print info_format, 'Total pressure', '|', pressure_t
            print '(a54)', table_string
        end if
        
    end subroutine total_energy_with_corrections
    
    
end program etf