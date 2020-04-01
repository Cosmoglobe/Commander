module comm_zodi_mod
    !   """
    !   Module which returns the zodiacal light emission computed for a given
    !   chunk of time-ordered data.
    !
    !   Main Methods
    !   ------------
    !   initialize_zodi_mod(cpar)
    !       Initializes the zodi_mod. Pre-computes galactic to ecliptic
    !       pixel coordinates, and initializes the different ZodiComponents.
    !   compute_zodi_template(nside, pix, nu, s_zodi)
    !       Routine which calculates and returns the zodiacal emission over
    !       a line-of-sight for a chunck of time-ordered data at a given
    !       frequency nu.
    !
    !   """
    use comm_utils
    use comm_param_mod

    implicit none

    private
    public :: initialize_zodi_mod, compute_zodi_template

    integer(i4b) :: n_LOS
    real(dp)     :: T0, const1, const2, delta
    real(dp)     :: R_max, R_sat
    real(dp)     :: t1, t2
    real(dp), dimension(:), allocatable :: x, y, z
    real(dp), dimension(:), allocatable :: blackbody, emission
    real(dp), dimension(:), allocatable :: density, tot_density

    ! =========================================================================
    !                      ZodiComponent Class Definition
    ! =========================================================================
    type, abstract :: ZodiComponent
        ! Pointers to the next/prev links in the linked list
        class(ZodiComponent), pointer :: nextLink => null()
        class(ZodiComponent), pointer :: prevLink => null()

        ! Shared component variables
        real(dp)                :: emissivity
        real(dp)                :: x0, y0, z0
        real(dp)                :: Incl, Omega
        real(dp), allocatable   :: sinOmega, cosOmega, sinIncl, cosIncl

        contains
            ! Shared component procedures
            procedure(init),     deferred :: initialize
            procedure(getDens), deferred :: getDensity

            ! Linked list procedures
            procedure :: next
            procedure :: add
    end type ZodiComponent

    abstract interface
        subroutine init(self)
            ! Routine which initializes and precomputes frequently used values
            ! for the a ZodiComponent.
            import dp, ZodiComponent
            class(ZodiComponent)  :: self
        end subroutine init

        subroutine getDens(self, x, y, z, density, lon)
            ! Routine which computes the density of a ZodiComponent at a
            ! given x, y, z coordinate.
            import i4b, dp, ZodiComponent
            class(ZodiComponent)                :: self
            real(dp), intent(in),  dimension(:) :: x, y, z
            real(dp), intent(out), dimension(:) :: density
            real(dp), intent(in),  optional     :: lon
            real(dp)                            :: xprime, yprime, zprime
            real(dp)                            :: R, Z_c
        end subroutine getDens
    end interface

    ! Individual class components
    ! -------------------------------------------------------------------------
    type, extends(ZodiComponent) :: Cloud
        real(dp)  :: n0, alpha, beta, gamma, mu

        contains
            procedure :: initialize => initializeCloud
            procedure :: getDensity => getDensityCloud
    end type Cloud

    type, extends(ZodiComponent) :: Band
        real(dp)                :: n0, Dz, Dr, R0, Vi, Vr, P_i, P_r
        real(dp), allocatable   :: ViInv, DrInv, DzRinv

        contains
            procedure :: initialize => initializeBand
            procedure :: getDensity => getDensityBand
    end type Band

    type, extends(ZodiComponent) :: Ring
        real(dp)                :: nsr, Rsr, sigmaRsr, sigmaZsr
        real(dp), allocatable   :: sigmaRsr2Inv, sigmaZsrInv

        contains
            procedure :: initialize => initializeRing
            procedure :: getDensity => getDensityRing
    end type Ring

    type, extends(ZodiComponent) :: Feature
        real(dp)                :: ntf, Rtf, sigmaRtf, sigmaZtf, thetatf, &
                                   sigmaThetatf
        real(dp), allocatable   :: thetatfR,  sigmaRtfInv, sigmaZtfInv, &
                                   sigmaThetatfRinv
        contains
            procedure :: initialize => initializeFeature
            procedure :: getDensity => getDensityFeature
    end type Feature

    ! Initializing global ZodiComponent list and objects
    class(ZodiComponent), pointer :: comp_list => null()
    type(Cloud),          target  :: cloud_comp
    type(Band),           target  :: band1_comp, band2_comp, band3_comp
    type(Ring),           target  :: ring_comp
    type(Feature),        target  :: feature_comp

    ! Derived type which stores dynamically sized pixel coordinate maps
    ! -------------------------------------------------------------------------
    type Vector3D
        real(dp), dimension(:,:), allocatable :: elements
    end type Vector3D

    type RaggedArray
        type(Vector3D), dimension(:), allocatable :: vectors
    end type RaggedArray
    type(RaggedArray) :: coord_maps

contains
    ! =========================================================================
    !                          Main zodi_mod Routines
    ! =========================================================================
    subroutine initialize_zodi_mod(cpar)
        !   """
        !   Routine which initializes the zodi_mod. Galactic to ecliptic
        !   x,y,z coordinates are precomputed to reduce computation at later
        !   stages.
        !
        !   Arguments:
        !   ----------
        !   cpar: derived type
        !       Object containing parameters from the parameterfile.
        !
        !   """
        implicit none

        type(comm_params), intent(in)  :: cpar
        class(ZodiComponent), pointer  :: comp

        integer(i4b) :: i, j, npix, nside
        logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, &
                        use_ring, use_feature

        real(dp),     dimension(3)                :: vec
        real(dp),     dimension(3,3)              :: gal2ecl_matrix
        integer(i4b), dimension(:),   allocatable :: ds_nside_unique, &
                                                     nside_unique
        real(dp),     dimension(:,:), allocatable :: ecliptic_vec

        ! Model parameters
        ! ---------------------------------------------------------------------
        ! Parameters for the temperature at a radial distance from the Sun
        T0 = 286.d0         ! temperature at 1 AU
        delta = 0.46686260  ! rate at which temperature falls with radius

        ! Line-of-sight integration parameters
        R_max = 6            ! max radial integration distance from the Sun [AU]
        R_sat = 1.1          ! satellite radial distance from the Sun [AU]
        n_LOS = 50           ! n integration steps

        ! Zodi component selection
        use_cloud = .true.
        use_band1 = .false.
        use_band2 = .false.
        use_band3 = .false.
        use_ring = .false.
        use_feature = .false.

        ! Initializing zodi components
        ! ---------------------------------------------------------------------
        if (use_cloud == .true.) then
            cloud_comp = Cloud(emissivity=1.d0, x0=0.011887801, y0=0.0054765065, &
                               z0=-0.0021530908, Incl=2.0335188, Omega=77.657956, &
                               n0=1.1344374e-7, alpha=1.3370697, beta=4.1415004, &
                               gamma=0.94206179, mu=0.18873176)
            comp => cloud_comp
            call add2Complist(comp)
        end if

        if (use_band1 == .true.) then
            band1_comp = Band(emissivity=1.d0, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=0.56438265, Omega=80d0, n0=5.5890290d-10, &
                              Dz=8.7850534, Dr=1.5, R0=3.d0, Vi=0.1, Vr=0.05, &
                              P_i=4.d0, P_r=1.d0)
            comp => band1_comp
            call add2Complist(comp)
        end if

        if (use_band2 == .true.) then
            band2_comp = Band(emissivity=1.d0, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=1.2, Omega=30.347476, n0=1.9877609d-09, &
                              Dz=1.9917032, Dr=0.94121881, R0=3.d0, &
                              Vi=0.89999998, Vr=0.15, P_i=4.d0, P_r=1.d0)
            comp => band2_comp
            call add2Complist(comp)
        end if

        if (use_band3 == .true.) then
            band3_comp = Band(emissivity=1.d0, x0=0.d0, y0=0.d0, z0=0.d0, &
                              Incl=0.8, Omega=80.0, n0=1.4369827d-10,     &
                              Dz=15.0, Dr=1.5, R0=3.d0, Vi=0.05, Vr=-1.0, &
                              P_i=4.d0, P_r=1.d0)
            comp => band3_comp
            call add2Complist(comp)
        end if

        if (use_ring == .true.) then
            ring_comp = Ring(emissivity=1.d0, x0=0.d0, y0=0.d0, z0=0.d0, &
                             Incl=0.48707166d0, Omega=22.27898d0, &
                             nsr=1.8260528d-8, Rsr=1.0281924d0, &
                             sigmaRsr=0.025d0, sigmaZsr=0.054068037d0)
            comp => ring_comp
            call add2Complist(comp)
        end if

        if (use_feature == .true.) then
            feature_comp = Feature(emissivity=1.d0, x0=0.d0, y0=0.d0, z0=0.d0, &
                                   Incl=0.48707166d0, Omega=22.27898d0, &
                                   ntf=2.0094267d-8, Rtf=1.0579183d0, &
                                   sigmaRtf=0.10287315d0, sigmaZtf=0.091442964d0, &
                                   thetatf=-10.d0, sigmaThetatf=12.115211d0)
            comp => feature_comp
            call add2Complist(comp)
        end if

        ! Executes initialization routines for all activated components
        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
        end do

        ! Allocating line-of-sight related arrays
        ! ---------------------------------------------------------------------
        allocate(x(n_LOS))
        allocate(y(n_LOS))
        allocate(z(n_LOS))
        allocate(density(n_LOS))
        allocate(tot_density(n_LOS))
        allocate(blackbody(n_LOS))
        allocate(emission(n_LOS))

        ! Precomputing ecliptic to galactic coordinates per pixels for all
        ! relevant nsides
        ! ---------------------------------------------------------------------

        ! Get all unique LFI nsides from cpar (This part should be updated when
        ! BeyondPlanck moves on to using HFI data)
        allocate(ds_nside_unique(cpar%numband))
        j = 0
        do i = 1, cpar%numband
            if (trim(cpar%ds_tod_type(i)) /= 'none') then
                if (any(ds_nside_unique /= cpar%ds_nside(i))) then
                    j = j + 1
                    ds_nside_unique(j) = cpar%ds_nside(i)
                end if
            end if
        end do

        allocate(coord_maps%vectors(j))
        allocate(nside_unique(j))
        do i = 1, j
            nside_unique(i) = ds_nside_unique(i)
        end do
        deallocate(ds_nside_unique)

        ! Getting pixel coordinates through HEALPix pix2vec_ring for each
        ! relevant nside
        call getEcl2GalMatrix(gal2ecl_matrix)
        do i = 1, size(nside_unique)
            nside = nside_unique(i)
            npix = nside2npix(nside)

            allocate(coord_maps%vectors(i)%elements(npix,3))
            allocate(ecliptic_vec(npix, 3))
            do j = 0, npix-1
                call pix2vec_ring(nside, j, vec)
                ecliptic_vec(j+1,1) = vec(1)
                ecliptic_vec(j+1,2) = vec(2)
                ecliptic_vec(j+1,3) = vec(3)
            end do

            ! Transforming to ecliptic coordinates
            coord_maps%vectors(i)%elements = matmul(ecliptic_vec, gal2ecl_matrix)
            deallocate(ecliptic_vec)
        end do

    end subroutine initialize_zodi_mod

    subroutine compute_zodi_template(nside, pix, sat_pos, nu, s_zodi)
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

        integer(i4b),                   intent(in)  :: nside
        integer(i4b), dimension(1:,1:), intent(in)  :: pix
        real(dp),     dimension(2),     intent(in)  :: sat_pos
        real(dp),     dimension(1:),    intent(in)  :: nu
        real(sp),     dimension(1:,1:), intent(out) :: s_zodi

        integer(i4b) :: i, j, k, n_det, n_tod
        real(dp)     :: x0, y0, z0
        real(dp)     :: u_x, u_y, u_z
        real(dp)     :: x1, y1, z1
        real(dp)     :: dx, dy, dz
        real(dp)     :: s, ds
        real(dp)     :: R_squared, R_cos_theta, R_LOS
        real(dp)     :: integral
        real(dp)     :: longitude_sat
        real(dp), dimension(:,:), allocatable :: coord_map

        ! Resetting arrays for each call
        s_zodi(:,:) = 0.0
        density(:) = 0.0
        tot_density(:) = 0.0
        blackbody(:) = 0.0
        emission(:) = 0.0

        n_tod = size(pix,1)  ! n time-ordered data for current chunk
        n_det = size(pix,2)  ! n detectors

        ! Selecting coordinate map containing pixel coordinates
        coord_map = getCoordMap(nside)

        ! Satellite longitude in radians
        longitude_sat = sat_pos(1)*deg2rad

        ! Computing the zodiacal emission for current time-ordered data chunk
        do j = 1, n_det
            ! Computing the consant terms in Planck's law B(T)
            const1 = (2.d0*h*217.d9**3)/(c*c)
            ! const1 = (2.d0*h*nu(j)**3)/(c*c)
            const2 = (h*217.d9)/k_B
            ! const2 = (h*nu(j))/k_B

            do i = 1, n_tod
                ! Computing heliocentric satellite coordinates (x0, y0, z0)
                x0 = R_sat*cos(longitude_sat)
                y0 = R_sat*sin(longitude_sat)
                z0 = 0.d0

                ! Looking up pre-computed HEALPix unit vector pointing to
                ! pixel at infinity
                u_x = coord_map(pix(i,j),1)
                u_y = coord_map(pix(i,j),2)
                u_z = coord_map(pix(i,j),3)

                ! Finding the coordinates at the end of the line-of-sight
                ! (x1, y1, z1) used to construct the line-of-sight array.
                ! The proof for these calculations are provided in my thesis
                R_squared = x0**2 + y0**2 + z0**2
                R_cos_theta = x0*u_x + y0*u_y + z0*u_z
                R_LOS = - R_cos_theta + sqrt(R_cos_theta**2 - R_squared + R_max**2)

                x1 = x0 + R_LOS*u_x
                y1 = y0 + R_LOS*u_y
                z1 = z0 + R_LOS*u_z

                ! Constructing line-of-sight array
                dx = (x1-x0)/(n_LOS-1)
                dy = (y1-y0)/(n_LOS-1)
                dz = (z1-z0)/(n_LOS-1)
                ds = sqrt(dx**2 + dy**2 + dz**2)

                do k = 1, n_LOS
                    x(k) = x0 + (k-1)*dx
                    y(k) = y0 + (k-1)*dy
                    z(k) = z0 + (k-1)*dz

                    ! Computing radial distance from the Sun
                    s = sqrt(x(k)**2 + y(k)**2 + z(k)**2)

                    ! Computing blackbody emission B(T) at radial distance s
                    ! from the Sun through the Planck function
                    blackbody(k) = const1/(exp(const2/(T0*s**(-delta))) - 1.d0)
                end do

                ! Looping over all active zodi components
                comp => comp_list
                do while (associated(comp))

                    ! Get density at given x, y, z coordinate along the line-of-sight
                    call comp%getDensity(x, y, z, density, longitude_sat)

                    ! Adding up density along the line-of-sight
                    tot_density = tot_density + density

                    comp => comp%next()
                end do

                ! Emission produced by the density
                emission = tot_density*blackbody

                ! Integrating up the emission along the line-of-sight
                ! with the trapezoidal method
                call trapezoidal(emission, ds, n_LOS, integral)

                ! Saving emission
                s_zodi(i,j) = s_zodi(i,j) + integral
            end do
        end do

        ! Converting to MJs/sr before returning signal
        ! s_zodi = s_zodi*1d20
    end subroutine compute_zodi_template

    ! =========================================================================
    !                         Functions and subroutines
    ! =========================================================================
    subroutine getEcl2GalMatrix(matrix)
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
    end subroutine getEcl2GalMatrix

    function getCoordMap(nside) result(coord_map)
        ! Routine which selects coordinate transformation map based on resolution
        implicit none

        integer(i4b), intent(in) :: nside

        integer(i4b)                           :: i, npix
        real(dp), dimension(:,:), allocatable  :: coord_map

        npix = nside2npix(nside)
        do i = 1, size(coord_maps%vectors)
            if (size(coord_maps%vectors(i)%elements(:,1)) == npix) then
                coord_map = coord_maps%vectors(i)%elements
            end if
        end do
    end function getCoordMap

    subroutine trapezoidal(f, ds, n, result)
        ! Trapezoidal integration method
        implicit none

        real(dp), dimension(:), intent(in) :: f
        integer(i4b), intent(in)           :: n
        real(dp), intent(in)               :: ds
        real(dp), intent(out)              :: result

        result = 0.5*f(1) + 0.5*f(n)
        result = result + sum(f(2:n-1))
        result = result*ds
    end subroutine trapezoidal

    ! =========================================================================
    !                         Zodi Components Routines
    ! =========================================================================
    subroutine initializeCloud(self)
        implicit none
        class(Cloud) :: self

        self%sinOmega = sin(self%Omega * deg2rad)
        self%cosOmega = cos(self%Omega * deg2rad)
        self%sinIncl = sin(self%Incl * deg2rad)
        self%cosIncl = cos(self%Incl * deg2rad)
    end subroutine initializeCloud

    subroutine getDensityCloud(self, x, y, z, density, lon)
        implicit none

        class(Cloud) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), dimension(:), intent(out) :: density
        real(dp), intent(in), optional      :: lon

        integer(i4b) :: i
        real(dp)     :: R, Z_c, zeta, g
        real(dp)     :: xprime, yprime, zprime

        do i = 1, n_LOS
            xprime = x(i) - self%x0
            yprime = y(i) - self%y0
            Zprime = Z(i) - self%z0

            R = sqrt(xprime*xprime + yprime*yprime + zprime*zprime)
            Z_c = (xprime*self%sinOmega - yprime*self%cosOmega)*self%sinIncl &
                + zprime*self%cosIncl
            zeta = abs(Z_c)/R

            if (zeta < self%mu) then
                g = 0.5 * zeta * zeta / self%mu
            else
                g = zeta - 0.5 * self%mu
            end if

            density(i) = self%n0 * R**(-self%alpha) * exp(-self%beta &
                    * g**self%gamma) * self%emissivity
        end do

    end subroutine getDensityCloud

    subroutine initializeBand(self)
        implicit none
        class(Band) :: self

        self%ViInv = 1.d0/self%Vi
        self%DrInv = 1.d0/self%Dr
        self%DzRInv = 1.d0/(self%Dz * deg2rad)
        self%sinOmega = sin(self%Omega * deg2rad)
        self%cosOmega = cos(self%Omega * deg2rad)
        self%sinIncl = sin(self%Incl * deg2rad)
        self%cosIncl = cos(self%Incl * deg2rad)
    end subroutine initializeBand

    subroutine getDensityBand(self, x, y, z, density, lon)
        implicit none

        class(Band) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), dimension(:), intent(out) :: density
        real(dp), intent(in), optional      :: lon

        integer(i4b) :: i
        real(dp) :: xprime, yprime, zprime
        real(dp) :: R, Z_c
        real(dp) :: zeta, ZDz, ZDz2, ZDz4, ZDz6, ViTerm, WtTerm

        do i = 1, n_LOS
            xprime = x(i) - self%x0
            yprime = y(i) - self%y0
            zprime = z(i) - self%z0

            R = sqrt(xprime*xprime + yprime*yprime + zprime*zprime)
            Z_c = (xprime*self%sinOmega - yprime*self%cosOmega)*self%sinIncl &
                + zprime*self%cosIncl
            zeta = abs(Z_c)/R
            ZDz = zeta * self%DzRInv
            ZDz2 = ZDz * ZDz
            ZDz4 = ZDz2 * ZDz2
            ZDz6 = ZDz4 * ZDz2
            ViTerm = 1.d0 + ZDz4 * self%ViInv
            WtTerm = 1.d0 - exp(-(R*self%DrInv)**20)

            density(i) = self%n0 * exp(-ZDz6) * ViTerm * WtTerm * self%R0/R &
                    * self%emissivity
        end do
    end subroutine getDensityBand

    subroutine initializeRing(self)
        implicit none
        class(Ring) :: self

        self%sigmaRsr2Inv = 1.d0 / (self%sigmaRsr * self%sigmaRsr)
        self%sigmaZsrInv = 1.d0 / self%sigmaZsr
        self%sinOmega = sin(self%Omega * deg2rad)
        self%cosOmega = cos(self%Omega * deg2rad)
        self%sinIncl = sin(self%Incl * deg2rad)
        self%cosIncl = cos(self%Incl * deg2rad)
    end subroutine initializeRing

    subroutine getDensityRing(self, x, y, z, density, lon)
        implicit none

        class(Ring) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), dimension(:), intent(out) :: density
        real(dp), intent(in), optional      :: lon

        integer(i4b) :: i
        real(dp) :: xprime, yprime, zprime
        real(dp) :: R, Z_c

        do i = 1, n_LOS
            xprime = x(i) - self%x0
            yprime = y(i) - self%y0
            zprime = z(i) - self%z0

            R = sqrt(xprime*xprime + yprime*yprime + zprime*zprime)
            Z_c = (xprime*self%sinOmega - yprime*self%cosOmega)*self%sinIncl &
                + zprime*self%cosIncl

            density(i) = self%nsr * exp(-((R-self%Rsr)/self%sigmaRsr)**2 &
                       - abs(Z_c)*self%sigmaZsrInv) * self%emissivity
        end do
    end subroutine getDensityRing

    subroutine initializeFeature(self)
        implicit none
        class(Feature) :: self

        self%thetatfR = self%thetatf * deg2rad
        self%sigmaRtfInv = 1.d0 / self%sigmaRtf
        self%sigmaZtfInv = 1.d0 / self%sigmaZtf
        self%sigmaThetatfRinv = 1.d0  /(self%sigmaThetatf * deg2rad)
        self%sinOmega = sin(self%Omega * deg2rad)
        self%cosOmega = cos(self%Omega * deg2rad)
        self%sinIncl = sin(self%Incl * deg2rad)
        self%cosIncl = cos(self%Incl * deg2rad)
    end subroutine initializeFeature

    subroutine getDensityFeature(self, x, y, z, density, lon)
        implicit none

        class(Feature) :: self
        real(dp), dimension(:), intent(in)  :: x, y, z
        real(dp), dimension(:), intent(out) :: density
        real(dp), intent(in), optional      :: lon

        integer(i4b) :: i
        real(dp) :: xprime, yprime, zprime
        real(dp) :: R, Z_c
        real(dp) :: theta

        do i = 1, n_LOS
            xprime = x(i) - self%x0
            yprime = y(i) - self%y0
            zprime = z(i) - self%z0
            theta = atan2(y(i), x(i)) - (lon + self%thetatfR)

            ! Constraining the angle to the limit [-pi, pi]
            do while (theta < -pi)
                theta = theta + 2*pi
            end do
            do while (theta > pi)
                theta = theta - 2*pi
            end do

            R = sqrt(xprime*xprime + yprime*yprime + zprime*zprime)
            Z_c = (xprime*self%sinOmega - yprime*self%cosOmega)*self%sinIncl &
                + zprime*self%cosIncl

            density(i) = self%ntf * exp(-((R-self%Rtf)*self%sigmaRtfInv)**2 &
                    - abs(Z_c)*self%sigmaZtfInv &
                    - (theta*self%sigmaThetatfRinv)**2) * self%emissivity
        end do
    end subroutine getDensityFeature

    ! =========================================================================
    !                   Linked list routines (ZodiComponent)
    ! =========================================================================
    function next(self)
        ! Routine which selects the next link in the linked list
        class(ZodiComponent) :: self
        class(ZodiComponent), pointer :: next
        next => self%nextLink
    end function next

    subroutine add(self,link)
        ! Routine which add a new object and link to the linked list
        class(ZodiComponent), target  :: self
        class(ZodiComponent), pointer :: link
        class(ZodiComponent), pointer :: comp

        comp => self
        do while (associated(comp%nextLink))
            comp => comp%nextLink
        end do
        link%prevLink => comp
        comp%nextLink    => link
    end subroutine add

    subroutine add2Complist(comp)
        implicit none
        class(ZodiComponent), pointer :: comp

        if (.not. associated(comp_list)) then
            comp_list => comp
        else
            call comp_list%add(comp)
        end if
    end subroutine add2Complist
end module comm_zodi_mod
