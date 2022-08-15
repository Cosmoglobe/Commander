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
    implicit none

    private
    public :: initialize_zodi_mod, get_zodi_emission
    integer(i4b) :: GAUSS_QUAD_ORDER, n_comps
    real(dp) :: T_0, delta, LOS_CUT, EPS
    real(dp), dimension(:), allocatable :: unique_nsides, x_helio, y_helio, z_helio, R_los, gauss_weights
    real(dp), dimension(:, :), allocatable :: zodi_emission

    type, abstract :: ZodiComponent
        ! Abstract base Zodical component class.

        ! Pointers to the next/prev links in the linked list
        class(ZodiComponent), pointer :: next_link => null()
        class(ZodiComponent), pointer :: prev_link => null()

        ! Shared component variables
        real(dp) :: emissivity
        real(dp) :: x0, y0, z0
        real(dp) :: Incl, Omega
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
            real(dp) :: x_prime, y_prime, z_prime, R, Z_c
        end subroutine density_interface
    end interface

    ! Individual class components
    ! -------------------------------------------------------------------------
    type, extends(ZodiComponent) :: Cloud
        real(dp)  :: n0, alpha, beta, gamma, mu

        contains
            procedure :: initialize => initialize_cloud
            procedure :: get_density => get_density_cloud
    end type Cloud

    type, extends(ZodiComponent) :: Band
        real(dp)                :: n0, Dz, Dr, R0, Vi, Vr, P_i, P_r
        real(dp), allocatable   :: ViInv, DrInv, DzRinv

        contains
            procedure :: initialize => initialize_band
            procedure :: get_density => get_density_band
    end type Band

    type, extends(ZodiComponent) :: Ring
        real(dp)                :: nsr, Rsr, sigmaRsr, sigmaZsr
        real(dp), allocatable   :: sigmaRsr2Inv, sigmaZsrInv

        contains
            procedure :: initialize => initialize_ring
            procedure :: get_density => get_density_ring
    end type Ring

    type, extends(ZodiComponent) :: Feature
        real(dp)                :: ntf, Rtf, sigmaRtf, sigmaZtf, thetatf, sigmaThetatf
        real(dp), allocatable   :: thetatfR,  sigmaRtfInv, sigmaZtfInv, sigmaThetatfRinv
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

        integer(i4b) :: i, j, npix, nside
        logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature, use_unit_emissivity
        real(dp) :: emissivity
        real(dp), dimension(3) :: vec
        real(dp), dimension(3,3) :: ecliptic_to_galactic_matrix
        integer(i4b), dimension(:), allocatable :: sorted_unique_nsides
        real(dp), dimension(:,:), allocatable :: galactic_vec
        real(dp), dimension(6) :: em30, em44, em70, em100, em143, em217, em353, em545, em857

        T_0 = 286.d0 ! temperature at 1 AU
        delta = 0.46686260 ! rate at which temperature falls with radius
        LOS_CUT = 5.2
        GAUSS_QUAD_ORDER = 100
        EPS = 3.d-14

        use_cloud = .true.
        use_band1 = .false.
        use_band2 = .false.
        use_band3 = .false.
        use_ring = .false.
        use_feature = .false.

        use_unit_emissivity = .true.

        !       (cloud, band1, band2, band3, ring, feature)
        em857 = (/0.301, 1.777, 0.716, 2.870, 0.578, 0.423/)
        em545 = (/0.223, 2.235, 0.718 , 3.193, 0.591, -0.182/)
        em353 = (/0.168, 2.035, 0.436, 2.400, -0.211, 0.676/)
        em217 = (/0.031, 2.024, 0.338, 2.507, -0.185, 0.243/)
        em143 = (/-0.014, 1.463, 0.530, 1.794, -0.252, -0.002/)
        em100 = (/0.003, 1.129, 0.674, 1.106, 0.163, 0.252/)


        ! Initialize Zodi components
        if (use_unit_emissivity) then
            emissivity = 1.d0
        end if

        if (use_cloud) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(1)
            end if
            cloud_comp = Cloud(emissivity=emissivity, x0=0.011887801, y0=0.0054765065, &
                               z0=-0.0021530908, Incl=2.0335188, Omega=77.657956, &
                               n0=1.1344374e-7, alpha=1.3370697, beta=4.1415004, &
                               gamma=0.94206179, mu=0.18873176)
            comp => cloud_comp
            call add_component_to_list(comp)
        end if

        if (use_band1) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(2)
            end if
            band1_comp = Band(emissivity=emissivity, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=0.56438265, Omega=80d0, n0=5.5890290d-10, &
                              Dz=8.7850534, Dr=1.5, R0=3.d0, Vi=0.1, Vr=0.05, &
                              P_i=4.d0, P_r=1.d0)
            comp => band1_comp
            call add_component_to_list(comp)
        end if

        if (use_band2) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(3)
            end if
            band2_comp = Band(emissivity=emissivity, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=1.2, Omega=30.347476, n0=1.9877609d-09, &
                              Dz=1.9917032, Dr=0.94121881, R0=3.d0, &
                              Vi=0.89999998, Vr=0.15, P_i=4.d0, P_r=1.d0)
            comp => band2_comp
            call add_component_to_list(comp)
        end if

        if (use_band3) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(4)
            end if
            band3_comp = Band(emissivity=emissivity, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=0.8, Omega=80.0, n0=1.4369827d-10,     &
                              Dz=15.0, Dr=1.5, R0=3.d0, Vi=0.05, Vr=-1.0, &
                              P_i=4.d0, P_r=1.d0)
            comp => band3_comp
            call add_component_to_list(comp)
        end if

        if (use_ring) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(5)
            end if
            ring_comp = Ring(emissivity=emissivity, x0=0.d0, y0=0.d0, z0=0.d0, &
                             Incl=0.48707166d0, Omega=22.27898d0, &
                             nsr=1.8260528d-8, Rsr=1.0281924d0, &
                             sigmaRsr=0.025d0, sigmaZsr=0.054068037d0)
            comp => ring_comp
            call add_component_to_list(comp)
        end if

        if (use_feature) then
            if (.not. use_unit_emissivity) then
                emissivity = em857(6)
            end if
            feature_comp = Feature(emissivity=emissivity, x0=0.d0, y0=0.d0, z0=0.d0, &
                                   Incl=0.48707166d0, Omega=22.27898d0, &
                                   ntf=2.0094267d-8, Rtf=1.0579183d0, &
                                   sigmaRtf=0.10287315d0, sigmaZtf=0.091442964d0, &
                                   thetatf=-10.d0, sigmaThetatf=12.115211d0)
            comp => feature_comp
            call add_component_to_list(comp)
        end if

        ! Executes initialization routines for all activated components
        n_comps = 0
        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
            n_comps = n_comps + 1
        end do

        ! Allocating line-of-sight related arrays
        allocate(x_helio(GAUSS_QUAD_ORDER))
        allocate(y_helio(GAUSS_QUAD_ORDER))
        allocate(z_helio(GAUSS_QUAD_ORDER))
        allocate(R_los(GAUSS_QUAD_ORDER))
        allocate(zodi_emission(n_comps, GAUSS_QUAD_ORDER))
        allocate(gauss_weights(GAUSS_QUAD_ORDER))

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
            ! unit_vectors%vectors(i)%elements = galactic_vec
            deallocate(galactic_vec)
        end do

    end subroutine initialize_zodi_mod

    subroutine get_zodi_emission(nside, pix, sat_pos, nu, s_zodi)
        !   """
        !   Routine which computes the zodiacal light emission at a given nside
        !   resolution for a chunk of time-ordered data.
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
        !   nu: array
        !       Array containing the all detector frequencies with
        !       dimensions (n_det)
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
        real(dp), dimension(1:), intent(in) :: nu
        real(sp), dimension(1:,1:), intent(out) :: s_zodi

        integer(i4b) :: i, j, k, n_det, n_tod, pixel_index
        real(dp) :: u_x, u_y, u_z
        real(dp) :: x1, y1, z1
        real(dp) :: dx, dy, dz
        real(dp) :: x_obs, y_obs, z_obs, lon_earth, R_obs, R_max
        real(dp), dimension(:), allocatable :: tabulated_zodi
        real(dp), dimension(:,:), allocatable :: unit_vector_map

        allocate(tabulated_zodi(nside2npix(nside)))
        tabulated_zodi = 0.d0
        zodi_emission = 0.d0
        R_los = 0.d0
        gauss_weights = 0.d0
        s_zodi = 0.d0
        ! Extracting n time-orderd data and n detectors for current chunk
        n_tod = size(pix,1)
        n_det = size(pix,2)

        unit_vector_map = get_unit_vector_map(nside) ! TODO: This should give me ecliptic unit vectors from galactic pixel index

        x_obs = sat_pos(1)
        y_obs = sat_pos(2)
        z_obs = sat_pos(3)
        R_obs = sqrt(x_obs**2 + y_obs**2 + z_obs**2)
        lon_earth = atan(y_obs, x_obs) ! TODO: this currently returns sat lon and not earth lon (unless this is basicaly always the same)

        do j = 1, n_det
            do i = 1, n_tod
                pixel_index = pix(i, j) + 1! TODO: make sure the +1 is correct. Does healpix indices in fortran also start at 0?
                if (tabulated_zodi(pixel_index) == 0.d0) then
                    u_x = unit_vector_map(pixel_index, 1)
                    u_y = unit_vector_map(pixel_index, 2)
                    u_z = unit_vector_map(pixel_index, 3)

                    call get_R_max(u_x, u_y, u_z, x_obs, y_obs, z_obs, R_obs, R_max)
                    call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=GAUSS_QUAD_ORDER, x=R_los, w=gauss_weights)

                    x_helio = R_los * u_x + x_obs
                    y_helio = R_los * u_y + y_obs
                    z_helio = R_los * u_z + z_obs

                    call get_zodi_emission_los(x=x_helio, y=y_helio, z=z_helio, nu=nu(j), theta=lon_earth, s_zodi=zodi_emission)
                    do k = 1, n_comps
                        s_zodi(i, j) = s_zodi(i, j) + sum(zodi_emission(k, :) * gauss_weights)
                    end do

                    tabulated_zodi(pixel_index) = s_zodi(i, j)

                else
                    ! Looking up tabulated emission
                    s_zodi(i, j) = tabulated_zodi(pixel_index)
                end if

            end do
        end do

        ! Converting to MJy/sr
        s_zodi = s_zodi * 1d20

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

        T = T_0 * R ** (-delta)
    end subroutine get_dust_grain_temperature


    subroutine get_blackbody_emission(nu, T, b_nu)
        implicit none
        real(dp), intent(in) :: nu
        real(dp), dimension(:), intent(in) :: T
        real(dp), dimension(:), intent(out) :: b_nu
        real(dp) :: term1, term2

        term1 = (2.d0 * h * nu**3) / (c*c)
        term2 = (h * nu) / k_B
        b_nu = term1/(exp(term2/T) - 1.d0)
    end subroutine get_blackbody_emission

    subroutine get_zodi_emission_los(x, y, z, nu, theta, s_zodi)
        implicit none
        real(dp), dimension(:), intent(in) :: x, y, z
        real(dp), intent(in) :: nu, theta
        real(dp), dimension(:, :), intent(out) :: s_zodi

        integer(i4b) :: i
        real(dp), dimension(GAUSS_QUAD_ORDER) :: R_helio, dust_grain_temperature, blackbody_emission, los_density

        class(ZodiComponent), pointer :: comp

        R_helio = sqrt(x**2 + y**2 + z**2)

        call get_dust_grain_temperature(R=R_helio, T=dust_grain_temperature)
        call get_blackbody_emission(nu=nu, T=dust_grain_temperature, b_nu=blackbody_emission)


        comp => comp_list
        i = 1
        do while (associated(comp))
            call comp%get_density(x=x, y=y, z=z, theta=theta, n=los_density)
            s_zodi(i, :) = comp%emissivity * los_density * blackbody_emission
            comp => comp%next()
            i = i + 1
        end do

    end subroutine get_zodi_emission_los


    ! Deferred ZodiComponent procedures
    subroutine initialize_cloud(self)
        implicit none
        class(Cloud) :: self

        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%Incl * deg2rad)
        self%cos_incl = cos(self%Incl * deg2rad)
    end subroutine initialize_cloud

    subroutine get_density_cloud(self, x, y, z, theta, n)
        implicit none

        class(Cloud) :: self
        real(dp), dimension(:), intent(in) :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n

        integer(i4b) :: i
        real(dp) :: R, Z_c, zeta, g, x_prime, y_prime, z_prime

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x0
            y_prime = y(i) - self%y0
            z_prime = z(i) - self%z0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_c = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
            zeta = abs(Z_c)/R

            if (zeta < self%mu) then
                g = 0.5 * zeta * zeta / self%mu
            else
                g = zeta - 0.5 * self%mu
            end if

            n(i) = self%n0 * R**(-self%alpha) * exp(-self%beta * g**self%gamma)
        end do
    end subroutine get_density_cloud

    subroutine initialize_band(self)
        implicit none
        class(Band) :: self

        self%ViInv = 1.d0/self%Vi
        self%DrInv = 1.d0/self%Dr
        self%DzRInv = 1.d0/(self%Dz * deg2rad)
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%Incl * deg2rad)
        self%cos_incl = cos(self%Incl * deg2rad)
    end subroutine initialize_band

    subroutine get_density_band(self, x, y, z, theta, n)
        implicit none

        class(Band) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n

        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_c, zeta, ZDz, ZDz2, ZDz4, ZDz6, ViTerm, WtTerm

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x0
            y_prime = y(i) - self%y0
            z_prime = z(i) - self%z0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_c = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
            zeta = abs(Z_c)/R
            ZDz = zeta * self%DzRInv
            ZDz2 = ZDz * ZDz
            ZDz4 = ZDz2 * ZDz2
            ZDz6 = ZDz4 * ZDz2
            ViTerm = 1.d0 + ZDz4 * self%ViInv
            WtTerm = 1.d0 - exp(-(R*self%DrInv)**20)

            n(i) = self%n0 * exp(-ZDz6) * ViTerm * WtTerm * self%R0/R
        end do
    end subroutine get_density_band

    subroutine initialize_ring(self)
        implicit none
        class(Ring) :: self

        self%sigmaRsr2Inv = 1.d0 / (self%sigmaRsr * self%sigmaRsr)
        self%sigmaZsrInv = 1.d0 / self%sigmaZsr
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%Incl * deg2rad)
        self%cos_incl = cos(self%Incl * deg2rad)
    end subroutine initialize_ring

    subroutine get_density_ring(self, x, y, z, theta, n)
        implicit none

        class(Ring) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n

        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_c

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x0
            y_prime = y(i) - self%y0
            z_prime = z(i) - self%z0

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_c = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

            n(i) = self%nsr * exp(-((R - self%Rsr)/self%sigmaRsr)**2 - abs(Z_c)*self%sigmaZsrInv)
        end do
    end subroutine get_density_ring

    subroutine initialize_feature(self)
        implicit none
        class(Feature) :: self

        self%thetatfR = self%thetatf * deg2rad
        self%sigmaRtfInv = 1.d0 / self%sigmaRtf
        self%sigmaZtfInv = 1.d0 / self%sigmaZtf
        self%sigmaThetatfRinv = 1.d0  /(self%sigmaThetatf * deg2rad)
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%Incl * deg2rad)
        self%cos_incl = cos(self%Incl * deg2rad)
    end subroutine initialize_feature

    subroutine get_density_feature(self, x, y, z, theta, n)
        implicit none

        class(Feature) :: self
        real(dp), dimension(:), intent(in) :: x, y, z
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n

        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_c, theta_prime

        do i = 1, GAUSS_QUAD_ORDER
            x_prime = x(i) - self%x0
            y_prime = y(i) - self%y0
            z_prime = z(i) - self%z0
            theta_prime = atan2(y(i), x(i)) - (theta + self%thetatfR)

            ! Constraining the angle to the limit [-pi, pi]
            do while (theta_prime < -pi)
                theta_prime = theta_prime + 2*pi
            end do
            do while (theta_prime > pi)
                theta_prime = theta_prime - 2*pi
            end do

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_c = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

            n(i) = self%ntf * exp(-((R-self%Rtf)*self%sigmaRtfInv)**2 - abs(Z_c)*self%sigmaZtfInv - (theta*self%sigmaThetatfRinv)**2)
        end do
    end subroutine get_density_feature

    function unique_sort(array) result(unique_sorted_array)
        implicit none
        integer :: idx, min_x, max_x
        integer, intent(in), dimension(:) :: array
        integer, dimension(size(array)) :: unique
        integer, dimension(:), allocatable :: unique_sorted_array

        idx = 0
        min_x = minval(array) - 1
        max_x = maxval(array)
        do while (min_x<max_x)
            idx = idx + 1
            min_x = minval(array, mask=array > min_x)
            unique(idx) = min_x
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
