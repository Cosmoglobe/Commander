!===============================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
submodule (comm_zodi_mod) zodi_solcentric_smod
contains

    module subroutine update_zodi_lookups(tod, zodi_model, start_time, end_time)
    !Once a set of zodi_model parameters is generated, fills in all lookup 
    !tables (tod%pixcache%map_zodi) with observed zodi predicted by the new model 
    !
    ! Parameters:
    ! tod : class(comm_tod)
    !     The TOD object holding the spline objects to update.
    ! zodi_model : type(ZodiModel)
    !     The zodiacal emission model.
    ! start_time, end_time : real(sp)
    !     Start and end time of the mission in MJD.
    !
        implicit none
        class(comm_tod),        intent(inout)     :: tod
        type(ZodiModel),        intent(in)        :: zodi_model
        real(sp),               intent(in)        :: start_time, end_time 

        real(dp)                                  :: vec_sol(3), vec_ecl(3), obs_pos(3), earth_pos(3)
        integer(i4b)                              :: d, t, p, j, pix_sol, pix_ecl, nslices
        real(sp)                                  :: current_time, signal_zodi

        do d = 1, tod%ndet !d is detector index from now on
            nslices = ((end_time - start_time)/zodi_model%zodi_timeslice) + 2
            current_time = start_time !MJD of the timeslice under consideration
        
            do t = 1, nslices
                do p = 1, tod%pixcache%nobs_zodi !p is pixel solar centric index from now on

                    !observer and Earth position in Heliocentric Ecliptic coordinates given the observation time
                    do j = 1, 3
                        earth_pos(j) = splint_simple(tod%x_earth_spline(j), real(current_time, dp))
                        obs_pos(j) = splint_simple(tod%x_obs_spline(j), real(current_time, dp))
                    end do

                    !take a pixel and find the correspondent pointing in Ecliptic coords
                    pix_sol = tod%pixcache%ind2pix_zodi(p) !solar-centric
                    call pix2vec_ring(tod%pixcache%nside_zodi, pix_sol, vec_sol)
                    call solar_centered_to_ecliptic(current_time, obs_pos, vec_sol, vec_ecl)

                    !Evaluates zodi at current_time as measured in pix_ecl, then uses that to set zodi in pix_sol
                    call compute_zodi_model(tod, d, vec_ecl, obs_pos, earth_pos, zodi_model, signal_zodi)
                    tod%pixcache%map_zodi(p, t, d) = signal_zodi
                end do

                !Next timeslice for detector d
                if (t == (nslices-1)) then
                    !end-start might not be exactly divisible per nslices so I force the last slice to be end_time no matter what
                    current_time = end_time 
                else 
                    current_time = current_time + zodi_model%zodi_timeslice
                end if
            end do
        end do
    end subroutine update_zodi_lookups

    module subroutine get_zodi_from_lookups(tod, pixels, scan, det, zodi_model, final_zodi)
    !uses the lookup maps to evaluate zodi for a series of pixels,
    !for now only takes a linear interpolation between the two closest
    !timeslices to the observation time, might implement a nicer thing later
    !
    ! Parameters:
    ! tod : class(comm_tod)
    !     The TOD object holding the spline objects to update.
    ! pixels : integer, dimension(ntod)
    !     The pixel indices of each time-ordered observation in Ecliptic coordinates.
    ! scan : integer
    !     The scan number.
    ! det : integer
    !     Index of the detector under consideration.
    ! zodi_model : type(ZodiModel)
    !     The zodiacal emission model.
    ! final_zodi : integer, dimension(ntod)
    !     Vector where I store the actual zodi signal values (in Ecliptic coords).
        implicit none
        class(comm_tod),               intent(inout)    :: tod
        integer(i4b),    dimension(:), intent(in)       :: pixels
        integer(i4b),                  intent(in)       :: scan, det
        type(ZodiModel),               intent(in)       :: zodi_model
        real(sp),        dimension(:), intent(out)      :: final_zodi

        real(sp)                                        :: current_time, start_time, end_time, val1, val2, time1, time2
        integer(i4b)                                    :: ntod, n, pix_sol, pix_ecl, p, s1, s2
        real(dp)                                        :: vec_sol(3), vec_ecl(3), obs_pos(3)

        ntod = size(pixels)
        do n = 1, ntod
            !Find the pixel
            current_time = tod%scans(scan)%d(det)%downsamp_obs_time(n) !Angela: is this the right time??
            obs_pos = tod%scans(scan)%x0_obs !Ecliptic observer
            pix_ecl = pixels(n)
            call pix2vec_ring(tod%nside, pix_ecl, vec_ecl) 
            call ecliptic_to_solar_centered(current_time, obs_pos, vec_ecl, vec_sol)
            call pix2vec_ring(tod%pixcache%nside_zodi, pix_sol, vec_sol)
            p = tod%pixcache%pix2ind_zodi(pix_sol)

            !Find the timeslices that are relevant for interpolation
            start_time = tod%global_min_mjd
            end_time = tod%global_max_mjd
            s1 = ((current_time - start_time)/zodi_model%zodi_timeslice) + 1  !index of timeslice right before current_time
            s2 = s1 + 1 !index of timeslice right after current_time
            time1 = start_time + ((s1-1)*current_time)
            time2 = start_time + ((s2-1)*current_time)
            if (time2 > end_time) then !this might happen if s2 is the index of the last slice which is forced to be end_time
                time2 = end_time  
            end if 

            !Perform simple linear interpolation
            val1 = tod%pixcache%map_zodi(p, s1, det)
            val2 = tod%pixcache%map_zodi(p, s2, det)
            final_zodi(n) = ((val1*(time2-current_time))+(val2*(current_time-time1)))/(time2-time1)            
        end do
    end subroutine get_zodi_from_lookups

    module subroutine solar_centered_to_ecliptic(time, obs_pos, vec_sol, vec_ecl)
    !Takes a vector in solar-centric coordinates, computes the corrispondent 
    !vector in Ecliptic coordinates using Eulero's angles matrix.
    !
    ! Parameters:
    ! time : real(sp)
    !     Time of the observation in MJD.
    ! obs_pos(3) : real(dp)
    !     Position of the observer in Heliocentric Ecliptic coordinates.
    ! vec_sol(3) : real(dp)
    !     Observed point in solar-centric coordinates.
    ! vec_ecl(3) : real(dp)
    !     Observed point in Heliocentric Ecliptic coordinates.
    !
        implicit none
        real(sp),               intent(in)        :: time
        real(dp),               intent(in)        :: obs_pos(3), vec_sol(3) 
        real(dp),               intent(inout)     :: vec_ecl(3)

        real(dp)                                  :: x_sun(3), lon_sun, colat_sun, M_e2s(3,3), phi!, theta, sinth, costh, sinph, cosph 

        x_sun = - obs_pos !Sun's position at time t in Heliocentric Ecliptic coordinates
        call vec2ang(x_sun, colat_sun, lon_sun) !Sun's colatitude and longitude in Heliocentric Ecliptic coords

        !Eulero's angles
        phi = -lon_sun !the matrix is defined with -phi
        !theta = (PI/2.0)-colat_sun !for the rotation I want theta is latitude
        call compute_euler_matrix_zyz(phi, 0.d0, 0.d0, M_e2s) !latitude correction is a small effect so we can neglect it

        vec_ecl = matmul(transpose(M_e2s), vec_sol) !transpose matrix=inverse matrix for rotations

    end subroutine solar_centered_to_ecliptic

    module subroutine ecliptic_to_solar_centered(time, obs_pos, vec_ecl, vec_sol)
    !Takes a vector in Ecliptic coordinates, computes the corrispondent 
    !vector in solar-centric coordinates using Eulero's angles matrix.
    !
    ! Parameters:
    ! time : real(sp)
    !     Time of the observation in MJD.
    ! obs_pos(3) : real(dp)
    !     Position of the observer in Heliocentric Ecliptic coordinates.
    ! vec_sol(3) : real(dp)
    !     Observed point in solar-centric coordinates.
    ! vec_ecl(3) : real(dp)
    !     Observed point in Heliocentric Ecliptic coordinates.
    !
        implicit none
        real(sp),               intent(in)        :: time
        real(dp),               intent(in)        :: obs_pos(3), vec_ecl(3) 
        real(dp),               intent(inout)     :: vec_sol(3)

        real(dp)                                  :: x_sun(3), lon_sun, colat_sun, M_e2s(3,3), phi!, theta, sinth, costh, sinph, cosph 

        x_sun = - obs_pos !Sun's position at time t in Heliocentric Ecliptic coordinates
        call vec2ang(x_sun, colat_sun, lon_sun) !Sun's colatitude and longitude in Heliocentric Ecliptic coords

        !Eulero's angles
        phi = -lon_sun !the matrix is defined with -phi
        !theta = (PI/2.0)-colat_sun !for the rotation I want theta is latitude
        call compute_euler_matrix_zyz(phi, 0.d0, 0.d0, M_e2s) !latitude correction is a small effect so we can neglect it

        vec_sol = matmul(M_e2s, vec_ecl) 

    end subroutine ecliptic_to_solar_centered

    module subroutine compute_zodi_model(tod, det, vec_ecl, obs_pos, earth_pos, zodi_model, signal_zodi)
    !Computes zodi_model along los defined between obs_pos and pix_ecl 
    !using the bandpass of detector det, stores the value in signal_zodi.
    !
    ! Parameters:
    ! tod : class(comm_tod)
    !     The TOD object holding the spline objects to update.
    ! det : integer
    !     Index of the detector under consideration.
    ! vec_ecl(3) : real(dp)
    !     Vector of observed pixel in Heliocentric Ecliptic coordinates.
    ! obs_pos(3), earth_pos(3) : real(dp)
    !     Observer's and Earth's position in Ecliptic coordinates.
    ! zodi_model : type(ZodiModel)
    !     The zodiacal emission model.
    ! signal_zodi : real(dp)
    !     Computed zodi emission.
    !
        implicit none
        class(comm_tod),        intent(inout)     :: tod
        integer(i4b),           intent(in)        :: det
        real(dp),               intent(in)        :: vec_ecl(3), obs_pos(3), earth_pos(3)
        type(ZodiModel),        intent(in)        :: zodi_model
        real(sp),               intent(inout)     :: signal_zodi

        real(dp)                                  :: s_tot, s_scat, s_therm, earth_lon, R_obs, R_max, unit_vector(3), al, em, zodi_sed_ampl
        integer(i4b)                              :: k, l
        logical(lgt)                              :: scattering, thermal
        real(dp), allocatable, dimension(:)       :: b_nu, mbb_integrand

        R_obs = norm2(obs_pos)
        earth_lon = atan(earth_pos(2), earth_pos(1))

        !Angela: check this
        scattering = tod%central_freq >= zodi_model%nu_min_scatter
        thermal    = tod%central_freq <= zodi_model%nu_max_thermal

        !Angela: do I need cache?

        s_tot = 0.
        do k = 1, zodi_model%n_comps

            !Get line of sight integration range
            call get_sphere_intersection(vec_ecl, obs_pos, R_obs, comp_LOS(k)%R_min, R_min)
            call get_sphere_intersection(vec_ecl, obs_pos, R_obs, comp_LOS(k)%R_max, R_max)

            do l = 1, 3
                ! Convert quadrature range from [-1, 1] to [R_min, R_max]
                comp_LOS(k)%X_unit(l,:) = (0.5*(R_max-R_MIN))*comp_LOS(k)%gauss_nodes+(0.5*(R_max+R_MIN))
                comp_LOS(k)%X_unit(l,:) = comp_LOS(k)%X_unit(l, :)*vec_ecl(l)
                comp_LOS(k)%X(l,:) = comp_LOS(k)%X_unit(l, :)+obs_pos(l)
            end do
            comp_LOS(k)%R = norm2(comp_LOS(k)%X, dim=1)

            !Compute phase function if scattering is included
            if (scattering) then
                if (trim(zodi_model%phasefunc_type) == 'Wright') then
                    allocate(b_nu(size(tod%bp(0)%p%nu)))
                    call get_blackbody_emission(tod%bp(0)%p%nu, 5772.d0, b_nu) 
                    comp_LOS(k)%F_sol = tsum(tod%bp(0)%p%nu, tod%bp(0)%p%tau*b_nu)/comp_LOS(k)%R**2 !Angela: use detector bandpass for something?
                    deallocate(b_nu)
                else
                    comp_LOS(k)%F_sol = zodi_model%F_sun(tod%zodiband)/comp_LOS(k)%R**2 !Angela:what is tod%zodiband??
               end if
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call zodi_model%get_phase_function(comp_LOS(k)%Theta, tod%zodiband, comp_LOS(k)%Phi)
            end if

            ! Get dust grain temperature along the line of sight, and compute splined blackbody emission
            if (trim(zodi_model%phasefunc_type) == 'Wright' .and. k > 1) then
               comp_LOS(k)%T = exp(5.5301d0) / sqrt(comp_LOS(k)%R)
            else
               comp_LOS(k)%T = zodi_model%T_0 * exp(-zodi_model%delta * log(comp_LOS(k)%R))
            end if

            !Compute number density
            call zodi_model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)

            !Prepare scattering and emission
            s_scat = 0.; s_therm = 0.
            !Angela: scattering always computed in the old way (for now) 
            if (scattering) s_scat  = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5d0*(R_max-R_MIN)*1d20 !Angela: R_MIN or R_min???
            al = zodi_model%comps(k)%c%albedo(tod%id)   !Angela: id???

            ! BB+MBB thing only for cloud at this stage
            if ((trim(zodi_model%comp_labels(k)) == 'cloud') .and. (trim(zodi_model%phasefunc_type) /= 'Wright')) then 
            
               !Thermal emission=   BlackBody     for sed_cutoff<nu
               !                =   Modified BB   for nu<sed_cutoff
                if (thermal) then
                    allocate(mbb_integrand(size(comp_LOS(k)%B_nu)))

                    !call get_blackbody_emission(tod%bp(0)%p%nu, comp_LOS(k)%T, comp_LOS(k)%B_nu)  !Compute exact BlackBody, check frequency
                    call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

                    zodi_sed_ampl = zodi_model%comps(k)%c%sed_cutoff**(-zodi_model%comps(k)%c%sed_b) !amplitude of MBB such that MBB+BB function is continuous in nu_cutoff

                    if (zodi_model%comps(k)%c%sed_cutoff <= tod%central_freq) then !THIS SHOULD BE INTEGRATION OF BANDPASS, HERE ONLY SIMPLE VERSION
                        s_therm = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights)*0.5d0*(R_max-R_MIN)*1d20 !BB integrated along line of sight (I think)
                    else 
                        mbb_integrand = zodi_sed_ampl*tod%central_freq**zodi_model%comps(k)%c%sed_b*comp_LOS(k)%B_nu !HERE ONLY SIMPLE VERSION
                        s_therm = sum(comp_LOS(k)%n*mbb_integrand*comp_LOS(k)%gauss_weights)*0.5d0*(R_max-R_MIN)*1d20 !MBB integrated along line of sight (I think)
                    end if

                    deallocate(mbb_integrand)
                end if

                   !Total intensity 
                   s_tot = s_tot + al*s_scat + s_therm                
            else 

                !Compute thermal emission in the old way
                call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu) 
                if (thermal) s_therm = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights)* 0.5d0*(R_max-R_MIN)*1d20
                em = zodi_model%comps(k)%c%emissivity(tod%id)

                !Compute total intensity in the old way
                if (trim(zodi_model%phasefunc_type) == 'Wright') then
                    s_tot = s_tot + al*s_scat + em*s_therm
                else
                    s_tot = s_tot + al*s_scat + em*(1.d0-al)*s_therm
                end if
            end if
        end do

        ! Return total zodi signal
        signal_zodi = s_tot

    end subroutine compute_zodi_model

end submodule zodi_solcentric_smod

