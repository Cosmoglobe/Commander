!================================================================================
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

module comm_zodi_mod
    !   """
    !   The zodi module handles simulating and fitting zodiacal emission at tod level.
    !
    !   Methods
    !   -------
    !   initialize_zodi_mod
    !       Initializes the zodi_mod.

    !   get_zodi_emission
    !       Method that simulates the zodiacal emission observed by the instrument for
    !       a given pixel.
    !   """

    use comm_utils
    use comm_param_mod
    use comm_bp_mod
    use spline_1D_mod
    implicit none

    private
    public :: initialize_zodi_mod, get_zodi_emission
    integer(i4b) :: GAUSS_QUAD_ORDER
    real(dp) :: T_0, DELTA, LOS_CUT, EPS, PLANCK_TERM1_DELTA, PLANCK_TERM2_DELTA, DELTA_T_ZODI, previous_chunk_obs_time
    real(dp), dimension(:), allocatable :: UNIQUE_NSIDES, PLANCK_TERM1_BP, PLANCK_TERM2_BP
    real(dp), allocatable, dimension(:,:) :: tabulated_earth_pos
    real(dp), allocatable, dimension(:) :: tabulated_earth_time
    real(sp), allocatable, dimension(:) :: tabulated_zodi

    character(len=32) :: freq_correction_type
    type(spline_type) :: spline_x, spline_y, spline_z


    type, abstract :: ZodiComponent
        ! Abstract base Zodical component class.

        ! Pointers to the next/prev links in the linked list
        class(ZodiComponent), pointer :: next_link => null()
        class(ZodiComponent), pointer :: prev_link => null()

        ! Shared component variables
        real(dp) :: emissivity
        real(dp) :: x_0, y_0, z_0
        real(dp) :: incl, Omega
        real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl

        contains
            ! Shared component procedures
            procedure(initialize_interface), deferred :: initialize
            procedure(density_interface), deferred :: get_density

            ! Linked list procedures
            procedure :: next
            procedure :: add
    end type ZodiComponent

    abstract interface
        subroutine initialize_interface(self)
            import dp, ZodiComponent
            class(ZodiComponent)  :: self
        end subroutine initialize_interface

        subroutine density_interface(self, x, y, z, theta, n)
            ! Returns the dust density (n) of the component at heliocentric 
            ! coordinates (x, y, z) and the earths longitude (theta).

            import i4b, dp, ZodiComponent
            class(ZodiComponent) :: self
            real(dp), intent(in), dimension(:) :: x, y, z
            real(dp), intent(in) :: theta
            real(dp), intent(out), dimension(:) :: n
            real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane
        end subroutine density_interface
    end interface

    ! Individual class components
    ! -------------------------------------------------------------------------
    type, extends(ZodiComponent) :: Cloud
        real(dp)  :: n_0, alpha, beta, gamma, mu
        contains
            procedure :: initialize => initialize_cloud
            procedure :: get_density => get_density_cloud
    end type Cloud

    type, extends(ZodiComponent) :: Band
        real(dp)                :: n_0, delta_zeta, delta_r, R_0, v, p
        contains
            procedure :: initialize => initialize_band
            procedure :: get_density => get_density_band
    end type Band

    type, extends(ZodiComponent) :: Ring
        real(dp)                :: n_0, R_0, sigma_r, sigma_z
        contains
            procedure :: initialize => initialize_ring
            procedure :: get_density => get_density_ring
    end type Ring

    type, extends(ZodiComponent) :: Feature
        real(dp)                :: n_0, R_0, sigma_r, sigma_z, theta_0, sigma_theta
        contains
            procedure :: initialize => initialize_feature
            procedure :: get_density => get_density_feature
    end type Feature

    ! Initializing global ZodiComponent list and instances
    class(ZodiComponent), pointer :: comp_list => null()
    type(Cloud),          target  :: cloud_comp
    type(Band),           target  :: band1_comp, band2_comp, band3_comp
    type(Ring),           target  :: ring_comp
    type(Feature),        target  :: feature_comp

    ! Initializing container for precomputed pixel to unitvectors for each unique 
    ! nside in the Commander run.
    type UnitVector
        real(dp), dimension(:,:), allocatable :: elements
    end type UnitVector

    type UnitVectorList
        type(UnitVector), dimension(:), allocatable :: vectors
    end type UnitVectorList
    type(UnitVectorList) :: unit_vectors


contains
    subroutine initialize_zodi_mod(cpar)
        !   """
        !   Initialize the zodi module by allocating global variables, instanciating
        !   zodi components, and precomputing the pixel to unitvector maps.
        !
        !   Parameters
        !   ----------
        !   cpar:
        !      Parameter file variables.
        !
        !   """
        implicit none

        type(comm_params), intent(in) :: cpar
        class(ZodiComponent), pointer :: comp

        integer(i4b) :: i, j, npix, nside, unit, n_earthpos
        logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature, apply_color_correction
        logical(lgt) :: use_unit_emissivity
        character(len=1024) :: tabulated_earth_pos_filename
        real(dp) :: emissivity
        real(dp), dimension(3) :: vec
        real(dp), dimension(3,3) :: ecliptic_to_galactic_matrix
        integer(i4b), dimension(:), allocatable :: sorted_unique_nsides
        real(dp), dimension(:,:), allocatable :: galactic_vec
        real(dp), dimension(6) :: EMISSIVITY_PLANCK_30, EMISSIVITY_PLANCK_44, EMISSIVITY_PLANCK_70, &
                                  EMISSIVITY_PLANCK_100, EMISSIVITY_PLANCK_143, EMISSIVITY_PLANCK_217, &
                                  EMISSIVITY_PLANCK_353, EMISSIVITY_PLANCK_545, EMISSIVITY_PLANCK_857
        real(dp), dimension(6) :: EMISSIVITY_DIRBE_01, EMISSIVITY_DIRBE_02, EMISSIVITY_DIRBE_03, &
                                  EMISSIVITY_DIRBE_04, EMISSIVITY_DIRBE_05, EMISSIVITY_DIRBE_06,  &
                                  EMISSIVITY_DIRBE_07, EMISSIVITY_DIRBE_08, EMISSIVITY_DIRBE_09, &
                                  EMISSIVITY_DIRBE_10

        use_cloud = .true.
        use_band1 = .true.
        use_band2 = .true.
        use_band3 = .true.
        use_ring = .true.
        use_feature = .true.

        use_unit_emissivity = .false.

        ! freq_correction_type = "delta"
        ! freq_correction_type = "bandpass"
        freq_correction_type = "color"

        T_0 = 286.d0 ! temperature at 1 AU
        DELTA = 0.46686260d0 ! rate at which temperature falls with radius
        LOS_CUT = 5.2d0
        GAUSS_QUAD_ORDER = 100
        EPS = 3.d-14
        DELTA_T_ZODI = 0.5 ! clear zodi cache after 1 day

        ! Planck emissivities (cloud, band1, band2, band3, ring, feature)
        EMISSIVITY_PLANCK_857 = (/0.301d0, 1.777d0, 0.716d0, 2.870d0, 0.578d0, 0.423d0/)
        EMISSIVITY_PLANCK_545 = (/0.223d0, 2.235d0, 0.718d0 , 3.193d0, 0.591d0, -0.182d0/)
        EMISSIVITY_PLANCK_353 = (/0.168d0, 2.035d0, 0.436d0, 2.400d0, -0.211d0, 0.676d0/)
        EMISSIVITY_PLANCK_217 = (/0.031d0, 2.024d0, 0.338d0, 2.507d0, -0.185d0, 0.243d0/)
        EMISSIVITY_PLANCK_143 = (/-0.014d0, 1.463d0, 0.530d0, 1.794d0, -0.252d0, -0.002d0/)
        EMISSIVITY_PLANCK_100 = (/0.003d0, 1.129d0, 0.674d0, 1.106d0, 0.163d0, 0.252d0/)

        ! DIRBE emissivities (cloud, band1, band2, band3, ring, feature)
        EMISSIVITY_DIRBE_01 = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
        EMISSIVITY_DIRBE_02 = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
        EMISSIVITY_DIRBE_03 = (/1.6598924040649741d0, 1.6598924040649741d0, 1.6598924040649741d0, &
                                1.6598924040649741d0, 1.6598924040649741d0, 1.6598924040649741d0/)
        EMISSIVITY_DIRBE_04 = (/0.99740908486652979d0, 0.35926451958350442d0, 0.35926451958350442d0, &
                                0.35926451958350442d0, 1.0675116768340536d0, 1.0675116768340536d0/)
        EMISSIVITY_DIRBE_05 = (/0.95766914805948866d0, 1.0127926948497732d0, 1.0127926948497732d0, &
                                0.35926451958350442d0, 1.0608768682182081d0, 1.0608768682182081d0/)
        EMISSIVITY_DIRBE_06 = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
        EMISSIVITY_DIRBE_07 = (/0.73338832616768868d0, 1.2539242027824944d0, 1.2539242027824944d0, &
                                1.2539242027824944d0, 0.87266361378785184d0, 0.87266361378785184d0/)
        EMISSIVITY_DIRBE_08 = (/0.64789881802224070d0, 1.5167023376593836d0, 1.5167023376593836d0, &
                                1.5167023376593836d0, 1.0985346556794289d0, 1.0985346556794289d0/)
        EMISSIVITY_DIRBE_09 = (/0.67694205881047387d0, 1.1317240279481993d0, 1.1317240279481993d0, &
                                1.1317240279481993d0, 1.1515825707787077d0, 1.1515825707787077d0/)
        EMISSIVITY_DIRBE_10 = (/0.51912085401950736d0, 1.3996145963796358d0, 1.3996145963796358d0, &
                                1.3996145963796358d0, 0.85763800994217443d0, 0.85763800994217443d0/)

        ! Initialize Zodi components
        if (use_unit_emissivity) emissivity = 1.d0

        if (use_cloud) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(1)
            cloud_comp = Cloud(emissivity=emissivity, x_0=0.011887800744346281d0, y_0=0.0054765064662263777d0, &
                               z_0=-0.0021530908020710744d0, incl=2.0335188072390769d0, Omega=77.657955554097114d0, &
                               n_0=1.1344373881427960d-7, alpha=1.3370696705930281d0, beta=4.1415004157586637d0, &
                               gamma=0.94206179393358036d0, mu=0.18873176489090190d0)
            comp => cloud_comp
            call add_component_to_list(comp)
        end if

        if (use_band1) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(2)
            band1_comp = Band(emissivity=emissivity, x_0=0.0, y_0=0.0, z_0=0.0, &
                              incl=0.56438265154389733d0, Omega=80.d0, n_0=5.5890290403228370d-10, &
                              delta_zeta=8.7850534408713035d0, delta_r=1.5d0, R_0=3.d0, v=0.10000000149011612d0, p=4.d0)
            comp => band1_comp
            call add_component_to_list(comp)
        end if

        if (use_band2) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(3)
            band2_comp = Band(emissivity=emissivity, x_0=0.d0, y_0=0.d0, z_0=0.d0, &
                              incl=1.2000000476837158d0, Omega=30.347475578624532d0, n_0=1.9877609422590801d-09, &
                              delta_zeta=1.9917032425777641d0, delta_r=0.94121881201651147d0, R_0=3.d0, &
                              v=0.89999997615814209d0, p=4.d0)
            comp => band2_comp
            call add_component_to_list(comp)
        end if

        if (use_band3) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(4)
            band3_comp = Band(emissivity=emissivity, x_0=0.d0, y_0=0.d0, z_0=0.d0, &
                              incl=0.80000001192092896d0, Omega=80.d0, n_0=1.4369827283512384d-10, &
                              delta_zeta=15.d0, delta_r=1.5d0, R_0=3.d0, v=0.050000000745058060d0, p=4.d0)
            comp => band3_comp
            call add_component_to_list(comp)
        end if

        if (use_ring) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(5)
            ring_comp = Ring(emissivity=emissivity, x_0=0.d0, y_0=0.d0, z_0=0.d0, &
                             incl=0.48707166006819241d0, Omega=22.278979678854448d0, &
                             n_0=1.8260527826501675d-8, R_0=1.0281924326308751d0, &
                             sigma_r=0.025000000372529030d0, sigma_z=0.054068037356978099d0)
            comp => ring_comp
            call add_component_to_list(comp)
        end if

        if (use_feature) then
            if (.not. use_unit_emissivity) emissivity = EMISSIVITY_DIRBE_06(6)
            feature_comp = Feature(emissivity=emissivity, x_0=0.d0, y_0=0.d0, z_0=0.d0, &
                                   incl=0.48707166006819241d0, Omega=22.278979678854448d0, &
                                   n_0=2.0094267183590947d-8, R_0=1.0579182694524214d0, &
                                   sigma_r=0.10287314662396611d0, sigma_z=0.091442963768716023d0, &
                                   theta_0=-10.d0, sigma_theta=12.115210933938741d0)
            comp => feature_comp
            call add_component_to_list(comp)
        end if

        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
        end do

        ! Precompute unit vectors in ecliptic for all galactic pixel indices
        ! per unique data nside.
        call get_gal_to_ecl_conversion_matrix(ecliptic_to_galactic_matrix)

        sorted_unique_nsides = unique_sort(pack(cpar%ds_nside, cpar%ds_nside /= 0))
        allocate(unit_vectors%vectors(size(sorted_unique_nsides)))

        do i = 1, size(sorted_unique_nsides)
            nside = sorted_unique_nsides(i)
            npix = nside2npix(nside)
            allocate(unit_vectors%vectors(i)%elements(npix, 3))
            allocate(galactic_vec(npix, 3))
        
            do j = 0, npix - 1
                call pix2vec_ring(nside, j, vec)
                galactic_vec(j + 1, 1) = vec(1)
                galactic_vec(j + 1, 2) = vec(2)
                galactic_vec(j + 1, 3) = vec(3)
            end do
        
            unit_vectors%vectors(i)%elements = matmul(galactic_vec, ecliptic_to_galactic_matrix)
            deallocate(galactic_vec)
        end do

    ! Reading in tabulated earth position
    unit = getlun()
    tabulated_earth_pos_filename = trim(cpar%datadir)//'/'//trim("earth_pos_1980-2050_ephem_de432s.txt")
    open(unit, file=trim(tabulated_earth_pos_filename))
    read(unit, *) n_earthpos
    read(unit, *) ! skip header
    allocate(tabulated_earth_pos(3, n_earthpos))
    allocate(tabulated_earth_time(n_earthpos))
    do i = 1, n_earthpos
      read(unit,*) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
    end do
    close(unit)

    ! Create spline objects for the Earths position given an observatino time
    call spline_simple(spline_x, tabulated_earth_time, tabulated_earth_pos(1, :), regular=.true.)
    call spline_simple(spline_y, tabulated_earth_time, tabulated_earth_pos(2, :), regular=.true.)
    call spline_simple(spline_z, tabulated_earth_time, tabulated_earth_pos(3, :), regular=.true.)

    previous_chunk_obs_time = 0 ! initialize previous chunk observation time

    end subroutine initialize_zodi_mod

    subroutine get_zodi_emission(nside, pix, sat_pos, obs_time, bandpass, s_zodi)
        !   """
        !   Routine which computes the zodiacal light emission at a given nside
        !   resolution for a chunk of obs_time-ordered data.
        !
        !   Arguments:
        !   ----------
        !   nside: int
        !       Map grid resolution.
        !   pix: array
        !       Pixel array containing the pixels from which to compute the
        !       zodi signal with dimensions (n_tod, n_det)
        !   sat_pos: real
        !       Satellite longitude for given time-order data chunk
        !   obs_time: real
        !       Time of observation in MJD.
        !   bandpass: bandpass object
        !       bandpass object containing the updates bandpass for each detector.
        !
        !   Returns:
        !   --------
        !   s_zodi: array
        !       Zodiacal emission for current time-ordered data chunk
        !       with dimensions (n_tod, n_det)
        !
        !   """
        implicit none
        class(ZodiComponent), pointer :: comp

        integer(i4b), intent(in) :: nside
        integer(i4b), dimension(1:,1:), intent(in) :: pix
        real(dp), dimension(3), intent(in) :: sat_pos
        real(dp), intent(in) :: obs_time
        class(comm_bp_ptr), dimension(:), intent(in) :: bandpass
        real(sp), dimension(1:,1:), intent(out) :: s_zodi

        integer(i4b) :: i, j, k, n_detectors, n_tods, pixel_index, los_step
        real(dp) :: u_x, u_y, u_z, x1, y1, z1, dx, dy, dz, x_obs, y_obs, z_obs
        real(dp) :: lon_earth, R_obs, R_max, nu_det
        real(dp), dimension(3) :: earth_pos
        real(dp), dimension(:), allocatable :: blackbody_emission_delta, blackbody_emission_c, nu_ratio
        real(dp), dimension(:,:), allocatable :: unit_vector_map, blackbody_emission_bp, b_nu_ratio
        real(dp), dimension(GAUSS_QUAD_ORDER) :: x_helio, y_helio, z_helio, R_los, gauss_weights, &
                                                 R_helio, dust_grain_temperature, &
                                                 los_density, comp_emission, bp_integrated_blackbody_emission, b_nu_colorcorr

        ! Allocate tabulated zodi array if first chunk
        if (.not. allocated(tabulated_zodi)) allocate(tabulated_zodi(nside2npix(nside)))

        ! Reset quantites from previous chunk
        comp_emission = 0.d0
        R_los = 0.d0
        gauss_weights = 0.d0
        s_zodi = 0.d0

        print *, obs_time - previous_chunk_obs_time, DELTA_T_ZODI
        ! if observer has mobed by more than DELTA_T_ZODI days, reset the cached zodi values
        if (obs_time - previous_chunk_obs_time > DELTA_T_ZODI) tabulated_zodi = 0.d0

        ! Interpolate earths position given the obs_time and tabulated earth position
        earth_pos(1) = splint_simple(spline_x, obs_time)
        earth_pos(2) = splint_simple(spline_y, obs_time)
        earth_pos(3) = splint_simple(spline_z, obs_time)

        ! Extracting n time-orderd data and n detectors for current chunk
        n_tods = size(pix,1)
        n_detectors = size(pix,2)

        ! Get precomputed pixel to unit vector values given the nside
        unit_vector_map = get_unit_vector_map(nside)

        x_obs = sat_pos(1)
        y_obs = sat_pos(2)
        z_obs = sat_pos(3)
        R_obs = sqrt(x_obs**2 + y_obs**2 + z_obs**2)
        lon_earth = atan(earth_pos(2), earth_pos(1))

        select case (trim(freq_correction_type))
        case ("delta")
            allocate(blackbody_emission_delta(GAUSS_QUAD_ORDER))
            do j = 1, n_detectors
                PLANCK_TERM1_DELTA = (2 * h * bandpass(j)%p%nu_c**3) / (c*c)
                PLANCK_TERM2_DELTA = (h * bandpass(j)%p%nu_c)/ k_B
                do i = 1, n_tods
                    pixel_index = pix(i, j) + 1
                    if (tabulated_zodi(pixel_index) == 0.d0) then
                        u_x = unit_vector_map(pixel_index, 1)
                        u_y = unit_vector_map(pixel_index, 2)
                        u_z = unit_vector_map(pixel_index, 3)

                        call get_R_max(u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs, R_max)
                        call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=GAUSS_QUAD_ORDER, x=R_los, w=gauss_weights)

                        x_helio = R_los * u_x + x_obs
                        y_helio = R_los * u_y + y_obs
                        z_helio = R_los * u_z + z_obs
                        R_helio = sqrt(x_helio**2 + y_helio**2 + z_helio**2)

                        call get_dust_grain_temperature(R=R_helio, T=dust_grain_temperature)
                        call get_blackbody_emission_delta(T=dust_grain_temperature, b_nu=blackbody_emission_delta)

                        comp => comp_list
                        k = 1
                        do while (associated(comp))
                            call comp%get_density(x=x_helio, y=y_helio, z=z_helio, theta=lon_earth, n=los_density)
                            comp_emission = comp%emissivity * los_density * blackbody_emission_delta
                            s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission * gauss_weights)
                            comp => comp%next()
                            k = k + 1
                        end do
                        tabulated_zodi(pixel_index) = s_zodi(i, j)

                    else
                        ! Looking up tabulated emission
                        s_zodi(i, j) = tabulated_zodi(pixel_index)
                    end if

                end do
            end do
            deallocate(blackbody_emission_delta)

        case ("bandpass")
            do j = 1, n_detectors
                allocate(PLANCK_TERM1_BP(bandpass(j)%p%n))
                allocate(PLANCK_TERM2_BP(bandpass(j)%p%n))
                allocate(blackbody_emission_bp(GAUSS_QUAD_ORDER, bandpass(j)%p%n))
                PLANCK_TERM1_BP = (2 * h * bandpass(j)%p%nu**3) / (c*c)
                PLANCK_TERM2_BP = (h * bandpass(j)%p%nu)/ k_B
                do i = 1, n_tods
                    pixel_index = pix(i, j) + 1
                    if (tabulated_zodi(pixel_index) == 0.d0) then
                        u_x = unit_vector_map(pixel_index, 1)
                        u_y = unit_vector_map(pixel_index, 2)
                        u_z = unit_vector_map(pixel_index, 3)

                        call get_R_max(u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs, R_max)
                        call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=GAUSS_QUAD_ORDER, x=R_los, w=gauss_weights)

                        x_helio = R_los * u_x + x_obs
                        y_helio = R_los * u_y + y_obs
                        z_helio = R_los * u_z + z_obs
                        R_helio = sqrt(x_helio**2 + y_helio**2 + z_helio**2)

                        call get_dust_grain_temperature(R=R_helio, T=dust_grain_temperature)
                        call get_blackbody_emission_bp(T=dust_grain_temperature, b_nu=blackbody_emission_bp)
                        ! Bandpass integrate blackbody emission at each step along the line-of-sight
                        do los_step = 1, GAUSS_QUAD_ORDER
                            bp_integrated_blackbody_emission(los_step) = bandpass(j)%p%SED2F(blackbody_emission_bp(los_step, :))
                        end do

                        comp => comp_list
                        k = 1
                        do while (associated(comp))
                            call comp%get_density(x=x_helio, y=y_helio, z=z_helio, theta=lon_earth, n=los_density)
                            comp_emission = comp%emissivity * los_density * bp_integrated_blackbody_emission
                            s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission * gauss_weights)
                            comp => comp%next()
                            k = k + 1
                        end do
                        tabulated_zodi(pixel_index) = s_zodi(i, j)

                    else
                        ! Looking up tabulated emission
                        s_zodi(i, j) = tabulated_zodi(pixel_index)
                    end if

                end do
                deallocate(PLANCK_TERM1_BP, PLANCK_TERM2_BP, blackbody_emission_bp)
            end do

        case ("color")
            do j = 1, n_detectors
                allocate(PLANCK_TERM1_BP(bandpass(j)%p%n))
                allocate(PLANCK_TERM2_BP(bandpass(j)%p%n))
                allocate(b_nu_ratio(GAUSS_QUAD_ORDER, bandpass(j)%p%n))
                allocate(nu_ratio(bandpass(j)%p%n))
                allocate(blackbody_emission_c(GAUSS_QUAD_ORDER))
                allocate(blackbody_emission_bp(GAUSS_QUAD_ORDER, bandpass(j)%p%n))
                PLANCK_TERM1_BP = (2 * h * bandpass(j)%p%nu**3) / (c*c)
                PLANCK_TERM2_BP = (h * bandpass(j)%p%nu)/ k_B
                PLANCK_TERM1_DELTA = (2 * h * bandpass(j)%p%nu_c**3) / (c*c)
                PLANCK_TERM2_DELTA = (h * bandpass(j)%p%nu_c)/ k_B
                do i = 1, n_tods
                    pixel_index = pix(i, j) + 1
                    if (tabulated_zodi(pixel_index) == 0.d0) then
                        u_x = unit_vector_map(pixel_index, 1)
                        u_y = unit_vector_map(pixel_index, 2)
                        u_z = unit_vector_map(pixel_index, 3)

                        call get_R_max(u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs, R_max)
                        call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=GAUSS_QUAD_ORDER, x=R_los, w=gauss_weights)

                        x_helio = R_los * u_x + x_obs
                        y_helio = R_los * u_y + y_obs
                        z_helio = R_los * u_z + z_obs
                        R_helio = sqrt(x_helio**2 + y_helio**2 + z_helio**2)

                        call get_dust_grain_temperature(R=R_helio, T=dust_grain_temperature)
                        call get_blackbody_emission_delta(T=dust_grain_temperature, b_nu=blackbody_emission_c)
                        call get_blackbody_emission_bp(T=dust_grain_temperature, b_nu=blackbody_emission_bp)

                        do los_step = 1, GAUSS_QUAD_ORDER
                            b_nu_ratio(los_step, :) = blackbody_emission_bp(los_step, :) / blackbody_emission_c(los_step)
                            b_nu_colorcorr(los_step) = tsum(bandpass(j)%p%nu, b_nu_ratio(los_step, :) * bandpass(j)%p%tau)
                        end do
                        nu_ratio = bandpass(j)%p%nu_c / bandpass(j)%p%nu
                        b_nu_colorcorr = b_nu_colorcorr / tsum(bandpass(j)%p%nu, nu_ratio * bandpass(j)%p%tau)

                        comp => comp_list
                        k = 1
                        do while (associated(comp))
                            call comp%get_density(x=x_helio, y=y_helio, z=z_helio, theta=lon_earth, n=los_density)
                            comp_emission = comp%emissivity * los_density * blackbody_emission_c * b_nu_colorcorr
                            s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission * gauss_weights)
                            comp => comp%next()
                            k = k + 1
                        end do
                        tabulated_zodi(pixel_index) = s_zodi(i, j)

                    else
                        ! Looking up tabulated emission
                        s_zodi(i, j) = tabulated_zodi(pixel_index)
                    end if

                end do
                deallocate(PLANCK_TERM1_BP, PLANCK_TERM2_BP, blackbody_emission_bp, blackbody_emission_c, b_nu_ratio, nu_ratio)
            end do
        end select

        previous_chunk_obs_time = obs_time ! Store prevous chunks obs time

    end subroutine get_zodi_emission


    subroutine get_gal_to_ecl_conversion_matrix(matrix)
        ! Ecliptic to galactic rotation matrix
        implicit none
        real(dp), dimension(3,3) :: matrix

        matrix(1,1) =  -0.054882486d0
        matrix(1,2) =  -0.993821033d0
        matrix(1,3) =  -0.096476249d0
        matrix(2,1) =   0.494116468d0
        matrix(2,2) =  -0.110993846d0
        matrix(2,3) =   0.862281440d0
        matrix(3,1) =  -0.867661702d0
        matrix(3,2) =  -0.000346354d0
        matrix(3,3) =   0.497154957d0

        ! call invert_matrix_dp(matrix)
    end subroutine get_gal_to_ecl_conversion_matrix


    function get_unit_vector_map(nside) result(unit_vector_map)
        ! Routine which selects coordinate transformation map based on resolution
        implicit none
        integer(i4b), intent(in) :: nside
        integer(i4b) :: i, npix, nside_idx
        real(dp), dimension(:,:), allocatable :: unit_vector_map

        npix = nside2npix(nside)
        do i = 1, size(unit_vectors%vectors)
            if (size(unit_vectors%vectors(i)%elements(:,1)) == npix) then
                unit_vector_map = unit_vectors%vectors(i)%elements
            end if
        end do
    end function get_unit_vector_map


    subroutine get_R_max(u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs, R_max)
        ! Computes the length of the line of sight such that it stops exactly at LOS_CUT.
        implicit none
        real(dp), intent(in) :: u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs
        real(dp), intent(out) :: R_max
        real(dp) :: lon, lat, cos_lat, b, d, q

        lon = atan(u_y, u_x)
        lat = asin(u_z)
        cos_lat = cos(lat)

        b = 2.d0 * (x_obs * cos_lat * cos(lon) + y_obs * cos_lat * sin(lon))
        d = R_obs**2 - LOS_CUT**2
        q = -0.5d0 * b * (1.d0 + sqrt(b**2 - (4.d0 * d)) / abs(b))
        R_max = max(q, d / q)

    end subroutine get_R_max


    subroutine get_dust_grain_temperature(R, T)
        implicit none
        real(dp), dimension(:), intent(in) :: R
        real(dp), dimension(:), intent(out) :: T
        T = T_0 * R ** (-DELTA)
    end subroutine get_dust_grain_temperature


    subroutine get_blackbody_emission_bp(T, b_nu)
        implicit none
        real(dp), dimension(:), intent(in) :: T
        real(dp), dimension(:, :), intent(out) :: b_nu
        integer(i4b) :: i
        do i = 1, GAUSS_QUAD_ORDER
            b_nu(i, :) = PLANCK_TERM1_BP/(exp(PLANCK_TERM2_BP/T(i)) - 1.d0)
        end do
        b_nu = b_nu * 1d20 !Convert from W/s/m^2/sr to MJy/sr
    end subroutine get_blackbody_emission_bp

    subroutine get_blackbody_emission_delta(T, b_nu)
        implicit none
        real(dp), dimension(:), intent(in) :: T
        real(dp), dimension(:), intent(out) :: b_nu
        b_nu = (PLANCK_TERM1_DELTA/(exp(PLANCK_TERM2_DELTA/T) - 1.d0)) * 1d20 !Convert from W/s/m^2/sr to MJy/sr
    end subroutine get_blackbody_emission_delta

    ! Deferred ZodiComponent procedures
    subroutine initialize_cloud(self)
        implicit none
        class(Cloud) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_cloud

    subroutine get_density_cloud(self, x, y, z, theta, n)
        implicit none
        class(Cloud) :: self
        real(dp), dimension(:), intent(in) :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n
        integer(i4b) :: i
        real(dp) :: R, Z_midplane, zeta, g, x_prime, y_prime, z_prime

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x_0
            y_prime = y(i) - self%y_0
            z_prime = z(i) - self%z_0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
            zeta = abs(Z_midplane/R)

            if (zeta < self%mu) then
                g = (zeta * zeta) / (2.d0 * self%mu)
            else
                g = zeta - (0.5d0 * self%mu)
            end if

            n(i) = self%n_0 * R**(-self%alpha) * exp(-self%beta * g**self%gamma)
        end do
    end subroutine get_density_cloud

    subroutine initialize_band(self)
        implicit none
        class(Band) :: self
        self%delta_zeta = self%delta_zeta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_band

    subroutine get_density_band(self, x, y, z, theta, n)
        implicit none
        class(Band) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n
        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, zeta, zeta_over_delta_zeta, term1, term2, term3, term4

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x_0
            y_prime = y(i) - self%y_0
            z_prime = z(i) - self%z_0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
            zeta = abs(Z_midplane/R)

            zeta_over_delta_zeta = zeta / self%delta_zeta
            term1 = self%R_0 * self%n_0 / R
            term2 = exp(-(zeta_over_delta_zeta**6))

            ! Differs from eq 8 in K98 by a factor of 1/self.v. See Planck XIV
            ! section 4.1.2.
            term3 = 1.d0 + (zeta_over_delta_zeta**self%p) / self%v
            term4 = 1.d0 - exp(-((R / self%delta_r) ** 20))

            n(i) = term1 * term2 * term3 * term4
        end do
    end subroutine get_density_band

    subroutine initialize_ring(self)
        implicit none
        class(Ring) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_ring

    subroutine get_density_ring(self, x, y, z, theta, n)
        implicit none
        class(Ring) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n
        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, term1, term2

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x_0
            y_prime = y(i) - self%y_0
            z_prime = z(i) - self%z_0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

            term1 = -((R - self%R_0) ** 2) / self.sigma_r**2
            term2 = abs(Z_midplane/self.sigma_z)

            n(i) = self%n_0 * exp(term1 - term2)
        end do
    end subroutine get_density_ring

    subroutine initialize_feature(self)
        implicit none
        class(Feature) :: self
        self%theta_0 = self%theta_0 * deg2rad
        self%sigma_theta = self%sigma_theta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_feature

    subroutine get_density_feature(self, x, y, z, theta, n)
        implicit none
        class(Feature) :: self
        real(dp), dimension(:), intent(in) :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n
        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, theta_prime, exp_term

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x_0
            y_prime = y(i) - self%y_0
            z_prime = z(i) - self%z_0
            theta_prime = atan2(y(i), x(i)) - (theta + self%theta_0)

            ! Constraining the angle to the limit [-pi, pi]
            do while (theta_prime < -pi)
                theta_prime = theta_prime + 2.d0*pi
            end do
            do while (theta_prime > pi)
                theta_prime = theta_prime - 2.d0*pi
            end do

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

            exp_term = ((R - self%R_0) ** 2 / self%sigma_r**2) + (abs(Z_midplane) / self%sigma_z) + (theta_prime**2 / self%sigma_theta**2)

            n(i) = self%n_0 * exp(-exp_term)
        end do
    end subroutine get_density_feature

    function unique_sort(array) result(unique_sorted_array)
        implicit none
        integer :: idx, min_val, max_val
        integer, intent(in), dimension(:) :: array
        integer, dimension(size(array)) :: unique
        integer, dimension(:), allocatable :: unique_sorted_array

        idx = 0
        min_val = minval(array) - 1
        max_val = maxval(array)
        do while (min_val < max_val)
            idx = idx + 1
            min_val = minval(array, mask=array > min_val)
            unique(idx) = min_val
        enddo

        allocate(unique_sorted_array(idx), source=unique(1:idx))
    end function unique_sort


    ! ZodiComponent linked list methods
    function next(self)
        ! Routine which selects the next link in the linked list
        class(ZodiComponent) :: self
        class(ZodiComponent), pointer :: next
        next => self%next_link
    end function next

    subroutine add(self,link)
        ! Routine which add a new object and link to the linked list
        class(ZodiComponent), target  :: self
        class(ZodiComponent), pointer :: link
        class(ZodiComponent), pointer :: comp

        comp => self
        do while (associated(comp%next_link))
            comp => comp%next_link
        end do
        link%prev_link => comp
        comp%next_link    => link
    end subroutine add

    subroutine add_component_to_list(comp)
        implicit none
        class(ZodiComponent), pointer :: comp

        if (.not. associated(comp_list)) then
            comp_list => comp
        else
            call comp_list%add(comp)
        end if
    end subroutine add_component_to_list

end module comm_zodi_mod
