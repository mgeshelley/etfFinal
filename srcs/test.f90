!> @author
!> M. Shelley
!
! DESCRIPTION:
!> Module to hold test routines
module test

    use parameters
    use routines
    use rspace
    use strutinsky
    implicit none
    
    integer :: file_unit, file_unit_2 !< Unit numbers for opening files
    integer :: sample_profile !< Test density profile to use
    
    
    contains
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to read from "test_input.in" parameters used for test
    !> routines
    subroutine read_test_params()
        implicit none
        integer :: ierr
        
        namelist /test_params/ sample_profile
        
        open(unit=file_unit,file='test_input.in',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Read parameters into namelist
        read(unit=file_unit,nml=test_params)
        
        close(unit=file_unit,status='keep',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Print input parameters
        write (*,nml=test_params)
        
    end subroutine read_test_params
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to calculate the density profiles using sample parameters
    subroutine test_densities()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Densities test !!!'
        
        file_unit = 100
        open(unit=file_unit,file='test_densities.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        densities: do i=1, n
            select case(sample_profile)
                case(1)
                    ! Values from fit to profile for inner crust
                    ! (SLy4, Z=40, N=40000, Rbox=80fm)
                    rho_q(i,0) = calc_rho_q(0.01855_dp,0.10437_dp,6.923_dp,0.93_dp,1._dp,r(i))
                    rho_q(i,1) = calc_rho_q(0._dp,0.0328_dp,6.43_dp,0.64_dp,1._dp,r(i))
                    
                case(2)
                    ! Values from fit to profile for Pb208 (SLy4);
                    ! maximum radius = 25fm
                    rho_q(i,0) = calc_rho_q(0._dp,0.0915_dp,6.75_dp,0.58_dp,1._dp,r(i))
                    rho_q(i,1) = calc_rho_q(0._dp,0.0651_dp,6.60_dp,0.54_dp,1._dp,r(i))
                    
                case(3)
                    ! Values from fit to profile for Pb208 (SLy4);
                    ! maximum radius = 25fm; including gamma/=1
                    rho_q(i,0) = calc_rho_q(0._dp,0.0916_dp,6.9_dp,0.62_dp,1.2_dp,r(i))
                    rho_q(i,1) = calc_rho_q(0._dp,0.0654_dp,7.1_dp,0.65_dp,1.9_dp,r(i))
                    
                case(4)
                    ! Values from fit to profile for Pb208 (SLy4);
                    ! maximum radius = 25fm; fitted for correct particle numbers
                    rho_q(i,0) = calc_rho_q(0._dp,0.091617566995_dp,6.9121873559_dp,0.61842936209_dp,1.1944159607_dp,r(i))
                    rho_q(i,1) = calc_rho_q(0._dp,0.065396291244_dp,7.2164739134_dp,0.65994797291_dp,2.0261839357_dp,r(i))
                    
                case(5)
                    ! Values from fit to profile for Sn176 (SLy4);
                    ! maximum radius = 60fm; tuned for particle numbers
!                     rho_q(i,0) = calc_rho_q(0._dp,0.09831199_dp,6.54533682_dp,0.62897015_dp,1._dp,r(i))
!                     rho_q(i,1) = calc_rho_q(0._dp,0.04804052_dp,6.17824503_dp,0.45532648_dp,1._dp,r(i))
                    
                    ! Values from fit to profile for Sn300 (SLy4);
                    ! maximum radius = 60fm; tuned for particle numbers
!                     rho_q(i,0) = calc_rho_q(0.00013887_dp,0.09801905_dp,6.53141337_dp,0.61879205_dp,1._dp,r(i))
!                     rho_q(i,1) = calc_rho_q(0._dp,0.04795780_dp,6.18172782_dp,0.45572951_dp,1._dp,r(i))
                    
                    ! Values from fit to profile for Sn11100 (SLy4);
                    ! maximum radius = 60fm; tuned for particle numbers
!                     rho_q(i,0) = calc_rho_q(0.01204156_dp,0.10467315_dp,7.09336363_dp,0.78098322_dp,1._dp,r(i))
!                     rho_q(i,1) = calc_rho_q(0._dp,0.03664600_dp,6.75433567_dp,0.51364787_dp,1._dp,r(i))
                    
                    ! Values from fit to profile for Sn21900 (SLy4);
                    ! maximum radius = 60fm; tuned for particle numbers
!                     rho_q(i,0) = calc_rho_q(0.02394670_dp,0.10030844_dp,7.99890799_dp,0.88631534_dp,1._dp,r(i))
!                     rho_q(i,1) = calc_rho_q(0._dp,0.026225538_dp,7.55425390_dp,0.56802125_dp,1._dp,r(i))
                    
                case(6)
                    ! Values from fit to profile for Pb208 (SkM*);
                    ! maximum radius = 12fm; including gamma/=1; tuned for particle numbers
                    rho_q(i,0) = calc_rho_q(0._dp,0.09091936_dp,7.18597577_dp,0.65216328_dp,1.53824468_dp,r(i))
                    rho_q(i,1) = calc_rho_q(0._dp,0.06223230_dp,6.95562565_dp,0.53101747_dp,1.38297704_dp,r(i))
                    
                case default
                    print '(a)', 'STOP: no running mode selected'
                    stop
            end select
            
            rho_t(i) = rho_q(i,0) + rho_q(i,1)
            write(unit=file_unit,fmt=n_floats(4)) r(i), rho_t(i), rho_q(i,0), rho_q(i,1)
        end do densities
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_densities.dat"'
        
    end subroutine test_densities
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test some derivatives of the total and individual
    !> densities
    subroutine test_density_derivs()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Density derivatives test !!!'
        
        ! Write derivatives to file
        file_unit = 101
        open(unit=file_unit,file='test_derivatives.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        densities: do i=1, n
            write(unit=file_unit,fmt=n_floats(7)) r(i), del_rho_t(i), del2_rho_t(i), del_rho_q(i,0), &
                                                    & del2_rho_q(i,0), del_rho_q(i,1), del2_rho_q(i,1)
        end do densities
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_derivatives.dat"'
        
    end subroutine test_density_derivs
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test effective masses
    !> \f$fm^{-3}\f$.
    subroutine test_eff_mass()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Effective masses test !!!'
        
        file_unit = 102
        open(unit=file_unit,file='test_eff_mass.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write effective masses to file 'test_eff_mass.dat'
        densities: do i=1, n
            write(unit=file_unit,fmt=n_floats(3)) r(i), f_q(i,0), f_q(i,1)
        end do densities
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_eff_mass.dat"'
        
    end subroutine test_eff_mass
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test kinetic energy densities
    subroutine test_kinetic_densities()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Kinetic energy densities test !!!'
        
        file_unit = 103
        open(unit=file_unit,file='test_kinetic_densities.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write kinetic energy densities to file 'test_kinetic_densities.dat'
        radii: do i=1, n
            write(unit=file_unit,fmt=n_floats(4)) r(i), tau_ETF_t(i), tau_ETF_q(i,0), tau_ETF_q(i,1)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_kinetic_densities.dat"'
        
    end subroutine test_kinetic_densities
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write spin current densities and spin-orbit fields to file
    subroutine test_spin_current_densities()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Spin current densities test !!!'
        
        file_unit = 104
        open(unit=file_unit,file='test_spin_current_densities.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write spin current densities and spin-orbit fields to file
        ! 'test_spin_current_densities.dat'
        radii: do i=1, n
            write(unit=file_unit,fmt=n_floats(5)) r(i), J_q(i,0), J_q(i,1), W_q(i,0), W_q(i,1)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_spin_current_densities.dat"'
        
    end subroutine test_spin_current_densities
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write order-by-order comparison of kinetic densities to
    !> file
    subroutine orders_kinetic_densities()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Order-by-order kinetic energy densities test !!!'
        
        file_unit = 105
        open(unit=file_unit,file='order_by_order_kinetic_densities.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write kinetic energy densities to file
        ! 'order_by_order_kinetic_densities.dat'
        radii: do i=1, n
            write(unit=file_unit,fmt=n_floats(7)) r(i), & 
                                    & tau_TF_q(i,0), tau_TF_q(i,0)+tau_2_L_q(i,0)+tau_2_NL_q(i,0), tau_ETF_q(i,0), &
                                    & tau_TF_q(i,1), tau_TF_q(i,1)+tau_2_L_q(i,1)+tau_2_NL_q(i,1), tau_ETF_q(i,1)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "order_by_order_kinetic_densities.dat"'
        
    end subroutine orders_kinetic_densities 
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to write spin central potentials to file
    subroutine test_central_potentials()
        implicit none
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        
        print '(/a)', '!!! Central potentials test !!!'
        
        file_unit = 106
        open(unit=file_unit,file='test_central_potentials.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write central potentials to file 'test_central_potentials.dat'
        radii: do i=1, n
            write(unit=file_unit,fmt=n_floats(4)) r(i), U_q(i,0), U_q(i,1), U_q(i,1)-v_c_di(i)-v_c_ex(i)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Densities written to "test_central_potentials.dat"'
        
    end subroutine test_central_potentials
    
    
    !> @author M. Shelley
    !> @brief
    !> Utility subroutine that can be called anywhere at any time, for writing
    !> a quantity (evaluated over WS cell) to file
    !> @param[in] quantity  Array to be written to file (1:n)
    subroutine write_ws_cell_array_quantity(quantity)
        implicit none
        real(kind=dp),  dimension(1:n), intent(in)  :: quantity
        integer :: ierr !< I/O error storage
        integer :: i !< Loop index
        integer, save  :: test_counter = 0 !< Counts invocations of subroutine, for writing to different file name
        character(len=1024) :: filename !< Name for file
        
        ! Invocation number
        test_counter = test_counter + 1
        
        ! Create filename
        write (filename,"(a20,i1,a4)") "test_array_quantity_", test_counter, ".dat" 
        
        file_unit = 1000
        open(unit=file_unit,file=trim(filename),status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write quantity to file "test_array_quantity_#.dat"
        radii: do i=1, n
            write(unit=file_unit,fmt=n_floats(4)) r(i), quantity(i)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
    
    end subroutine write_ws_cell_array_quantity
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test calculation of neutron and proton particle numbers
    subroutine test_calc_particle_number()
        implicit none
        integer :: q !< isospin
        
        print '(/a)', '!!! Particle number test !!!'
        
        do q=0, 1
            print '(1x,a5,2(1x,a1),1x,f9.3)', 'Total',iso_string(q),'=',N_q(q)
        end do
        
        print '(1x,a9,1x,f9.3)', 'N total =', N_t
        
    end subroutine test_calc_particle_number
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test calculation of Coulomb energy
    subroutine test_calc_coulomb_energy()
        implicit none
        
        print '(/a)', '!!! Coulomb energy test !!!'
        
        print '(a17,1x,f9.3,1x,a3)', 'Coulomb energy =', e_coulomb, 'MeV'
        
    end subroutine test_calc_coulomb_energy
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test calculation of Skyrme energy
    subroutine test_calc_skyrme_energy()
        implicit none
        
        print '(/a)', '!!! Skyrme energy test !!!'
        
        print '(a16,1x,f10.3,1x,a3)', 'Skyrme energy =', e_skyrme, 'MeV'
        
    end subroutine test_calc_skyrme_energy
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test calculation of single particle energies
    subroutine test_sp_energies()
        implicit none
        character(1), dimension(0:10) :: l !< Orbital labels
        integer, dimension(0:1) :: file_unit_q
        integer :: ierr !< I/O error storage
        integer :: q !< Isospin
        integer :: i !< Loop index
        
        print '(/a)', '!!! Single particle energies test !!!'
        
        l = (/'s','p','d','f','g','h','i','j','k','l','m'/)
        
        file_unit_q(0) = 107
        open(unit=file_unit_q(0),file='test_sp_levels_n.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        file_unit_q(1) = 108
        open(unit=file_unit_q(1),file='test_sp_levels_p.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        isospin: do q=0, 1
            ! Call boundary routine for Pb208 with SLy4
            call boundary(dr,Lmax,Ecut(q),n,PS,NNST,JJP,LP,EP,q)
            
            ! Print states if found, otherwise print error message
            if (NNST .gt. 0) then
                ! Number of states found
                write(unit=file_unit_q(q),fmt='(a16,1x,2(a1),1x,i3)') 'States found for',iso_string(q),':',NNST
                
                ! Sort states
                call state_sort(JJP(1:NNST),LP(1:NNST),EP(1:NNST),PS(1:NNST,:))
                
                ! Print J, L, energies of single particle states
                write(unit=file_unit_q(q),fmt='(a5,2x,a12)') 'State','Energy (MeV)'
                states: do i=1, Nmaxstate
                    if (JJP(i) .ne. 0) then
                        write(unit=file_unit_q(q),fmt='(a1,i2,a2,1x,f12.6)') l(LP(i)),JJP(i),'/2',EP(i)
                    end if
                end do states
            else
                print '(2a)', space4,'- !!! WARNING: no states found !!!'
            end if
        end do isospin
        
        close(unit=file_unit_q(0),iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        close(unit=file_unit_q(1),iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Neutron single-particle levels written to "sp_levels_n.dat"'
        print '(2a)', space4,'- Proton single-particle levels written to "sp_levels_p.dat"'
        
    end subroutine test_sp_energies
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to test the Strutinsky integral correction
    subroutine test_si_correction()
        implicit none
        integer :: q !< Isospin
        
        print '(/a)', '!!! Strutinsky integral correction test !!!'
        
        isospin: do q=0, 1
            call calc_E_sc(q)
            
            print '(1x,a27,2(1x,a1),1x,f10.3,1x,a3)', 'Shell correction energy for',iso_string(q),'=', E_sc_q(q), 'MeV'
        end do isospin
        
        print '(1x,a31,1x,f10.3,1x,a3)', 'Total shell correction energy =', E_sc_t, 'MeV'
        print '(1x,a31,1x,f10.3,1x,a3)', 'Corrected total energy =', e_skyrme + E_sc_t, 'MeV'
        
    end subroutine test_si_correction
    
    
    
    !> @author M. Shelley
    !> @brief
    !> Subroutine to check different contributions to \f$\tau_{ETF}\f$ at 2nd
    !> order, for comparison with results in Bartel and Bencheikh (2002)
    subroutine bartel_bencheikh_benchmark()
        implicit none
        integer :: i !< Loop index
        integer :: ierr !< I/O error storage
        
        print '(/a)', 'Bartel and Bencheikh benchmark:'
        
        file_unit = 109
        open(unit=file_unit,file='bartel_bencheikh_benchmark.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        file_unit_2 = 110
        open(unit=file_unit_2,file='bartel_bencheikh_extra.dat',status='unknown',iostat=ierr)
        if (ierr.ne.0) stop
        
        ! Write neutron kinetic energy density contributions to file
        ! 'bartel_bencheikh_benchmark.dat'
        radii: do i=1, n
                write(unit=file_unit,fmt=n_floats(9)) r(i), tau_TF_q(i,0), tau_2_L_q(i,0)+tau_2_NL_q(i,0), &
                                            & tau_2_cont_q(i,1,0), tau_2_cont_q(i,2,0), tau_2_cont_q(i,3,0), tau_2_cont_q(i,4,0), &
                                            & J_2_q(i,0), J_4_q(i,0)
                
                write(unit=file_unit_2,fmt=n_floats(10)) r(i), rho_q(i,0), rho_q(i,1), &
                                            & U_q(i,0), U_q(i,1)-v_c_di(i)-v_c_ex(i), &
                                            & tau_4_no_spin_q(i,0)+tau_4_so_q(i,0), &
                                            & tau_4_cont_q(i,1,0), tau_4_cont_q(i,2,0), tau_4_cont_q(i,3,0), tau_4_cont_q(i,4,0)
        end do radii
        
        close(unit=file_unit,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        close(unit=file_unit_2,iostat=ierr,status='keep')
        if (ierr.ne.0) stop
        
        print '(2a)', space4,'- Neutron densities written to "bartel_bencheikh_benchmark.dat"'
        print '(2a)', space4,'- Other neutron densities and fields written to "bartel_bencheikh_extra.dat"'
        
    end subroutine bartel_bencheikh_benchmark
    
    
end module test