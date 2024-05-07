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
module comm_cmb_relquad_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_cmb_relquad_comp

  !**************************************************
  !           CMB relativistic quadrupole component
  !**************************************************
  type, extends (comm_template_comp) :: comm_cmb_relquad_comp
     integer(i4b) :: nside, lmax, ndet
     real(dp)     :: T_cmb
     real(dp)     :: v_sun, lon_sun, lat_sun
     class(comm_mapinfo),              pointer        :: info
     logical(lgt),        allocatable, dimension(:)   :: band_active
     real(dp),            allocatable, dimension(:,:) :: q            ! Bandpass integrated SED
     type(F_int_ptr),     allocatable, dimension(:,:) :: F_int        ! SED integrator
   contains
     procedure :: S             => evalSED_relquad
     procedure :: getBand       => evalRelquadBand
     procedure :: projectBand   => projectRelquadBand
     procedure :: update_template
     procedure :: read_definition_file
     procedure :: update_F_int  => updateIntF_relquad
     procedure :: initHDFComp   => initTemplateRelQuadHDF
  end type comm_cmb_relquad_comp

  interface comm_cmb_relquad_comp
     procedure constructor_relquad
  end interface comm_cmb_relquad_comp

contains

  function constructor_relquad(cpar, id, id_abs) result(c)
    !****************************************************************************************
    !> Purpose: Initializes and returns a comm_cmb_relquad_comp object
    ! 
    ! Inputs:
    !! @param[in]    cpar   -- parameter file structure, defined in comm_param_mod
    !! @param[in]    id     -- signal component counter, counting only active components
    !! @param[in]    id_abs -- absolute signal component counter, with order matching the parameter file
    ! 
    ! Outputs:
    ! @return        constructor -- allocated comm_cmb_relquad_comp object
    !****************************************************************************************
    implicit none
    class(comm_params),           intent(in) :: cpar
    integer(i4b),                 intent(in) :: id, id_abs
    class(comm_cmb_relquad_comp), pointer    :: c

    integer(i4b) :: i, j, k, n
    character(len=16), dimension(1000) :: comp_label

    ! Initialize general parameters
    allocate(c)

    c%class         = cpar%cs_class(id_abs)
    c%type          = cpar%cs_type(id_abs)
    c%label         = cpar%cs_label(id_abs)
    c%id            = id
    c%nmaps         = 1    ! Only used for CR book-keeping; must be 1 for templates
    c%outprefix     = trim(cpar%cs_label(id_abs))
    c%cg_scale      = 1.d0
    c%band          = 0
    c%myid          = cpar%myid_chain
    c%comm          = cpar%comm_chain
    c%numprocs      = cpar%numprocs_chain
    c%P             = [cpar%cs_P_gauss(id_abs,1,1), cpar%cs_P_gauss(id_abs,2,1)]
    c%init_from_HDF = cpar%cs_initHDF(id_abs)

    if (c%myid == 0) then
       c%ncr = 1
       c%x   = cpar%cs_theta_def(1,id_abs)
    else
       c%ncr = 0
    end if

    ! Component specific parameters
    c%T_cmb   = cpar%T_cmb
    c%v_sun   = 370.082d3   ! m/s;  Planck 2015
    c%lon_sun = 264.000d0   ! deg;  Planck 2015
    c%lat_sun = 48.24d0     ! deg;  Planck 2015
    c%nside   = 128
    c%lmax    = 2
    c%nmaps   = 1
    c%ndet    = maxval(data%ndet)
    c%info    => comm_mapinfo(c%comm, &
         & c%nside, c%lmax, c%nmaps, .false.)
    c%T       => comm_map(c%info)

    ! Construct spatial template
    call c%update_template

    ! Read band definition file
    allocate(c%band_active(numband))
    call c%read_definition_file(trim(cpar%cs_SED_template(1,id_abs)))

    ! Precompute mixmat for each band; temperature only
    allocate(c%F_int(numband,0:c%ndet))
    allocate(c%q(numband,0:c%ndet))

    ! Set up CG sampling groups
    allocate(c%active_samp_group(cpar%cg_num_samp_groups))
    c%active_samp_group = .false.
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(c%label) == trim(comp_label(j))) then
             c%active_samp_group(i) = .true.
             if (n == 1) c%cg_unique_sampgroup = i ! Dedicated sampling group for this component
             exit
          end if
       end do
    end do

  end function constructor_relquad
  
  !****************************************************************************************
  !> Purpose: Return SED of relativistic quadrupole correction; see Notari and Quartin (2015)
  ! 
  ! Inputs:
  !> @param[in]    self     -- current object
  !> @param[in]    nu       -- frequency in Hz
  !> @param[in]    band     -- not used; included for interface consistency
  !> @param[in]    pol      -- not used; included for interface consistency
  !> @param[in]    theta    -- not used; included for interface consistency
  ! 
  ! Outputs:
  !> @return       evalSED  -- SED in uK_RJ, measured in absolute units
  !****************************************************************************************
  function evalSED_relquad(self, nu, band, pol, theta)
    class(comm_cmb_relquad_comp), intent(in)           :: self
    real(dp),                     intent(in), optional :: nu
    integer(i4b),                 intent(in), optional :: band
    integer(i4b),                 intent(in), optional :: pol
    real(dp), dimension(1:),      intent(in), optional :: theta
    real(dp)                                           :: evalSED_relquad

    real(dp) :: x, cmb2RJ, beta

    beta    = self%v_sun/c
    x       = h*nu / (k_B*self%T_cmb)
    cmb2RJ  = (x**2 * exp(x)) / (exp(x)-1.d0)**2
    evalSED_relquad = cmb2RJ * beta**2 * 1.d6*self%T_cmb * (x/2.d0)*(exp(x)+1.d0)/(exp(x)-1.d0)

  end function evalSED_relquad


  !****************************************************************************************
  !> Purpose: Return component map as observed by a given detector, m_out = B*M*alm_in,
  !>          where B denotes beam convolution and M denotes the mixing matrix
  ! 
  ! Inputs:
  !> @param[in]     self     -- current object
  !> @param[in]     band     -- data set ID, counting only active sets
  !> @param[in]     amp_in   -- input amplitude alms, with size (1:nalm, nmaps). If not 
  !                             present, use the parameters stored in self%T%alm
  !> @param[in]     pix      -- not used here, included only for interface consistency
  !> @param[in]     alm_out  -- output array contains alms instead of pixelized sky map
  !> @param[in]     det      -- detector index; if not present, assume band-averaged output
  ! 
  ! Outputs:
  !> @return        evalRelquadBand -- map (or alms) of current component as seen by given detector
  !****************************************************************************************
  function evalRelquadBand(self, band, amp_in, pix, alm_out, det)
    implicit none
    class(comm_cmb_relquad_comp),                 intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: evalRelquadBand

    integer(i4b) :: i, j, np, nmaps, lmax, nmaps_comp, d
    logical(lgt) :: alm_out_
    real(dp)     :: t1, t2
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: m

    alm_out_ = .false.; if (present(alm_out)) alm_out_ = alm_out
    d = 0; if (present(det)) d = det

    if (.not. self%band_active(band)) then
       if (alm_out_) then
          if (.not. allocated(evalRelquadBand)) allocate(evalRelquadBand(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       else
          if (.not. allocated(evalRelquadBand)) allocate(evalRelquadBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       end if
       evalRelquadBand = 0.d0
       return
    end if

    ! Initialize amplitude map
    nmaps =  min(data(band)%info%nmaps, self%nmaps)
    info  => comm_mapinfo(data(band)%info%comm, data(band)%info%nside, data(band)%info%lmax, nmaps, nmaps==3)
    m     => comm_map(info)
    if (present(amp_in)) then
       m%alm(:,1:nmaps) = amp_in(:,1:nmaps)
    else
       call self%T%alm_equal(m)
    end if

    ! Scale to correct frequency through multiplication with mixing matrix
    m%alm = m%alm * self%q(band,d)
       
    ! Convolve with band-specific beam
    call data(band)%B(d)%p%conv(trans=.false., map=m)
    if (.not. alm_out_) call m%Y()

    ! Return correct data product
    if (alm_out_) then
       if (.not. allocated(evalRelquadBand)) allocate(evalRelquadBand(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) evalRelquadBand = 0.d0
       evalRelquadBand(:,1:nmaps) = m%alm(:,1:nmaps)
    else
       if (.not. allocated(evalRelquadBand)) allocate(evalRelquadBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) evalRelquadBand = 0.d0
       evalRelquadBand(:,1:nmaps) = m%map(:,1:nmaps)
    end if
       
    ! Clean up
    call m%dealloc(); deallocate(m)
    nullify(info)
    
  end function evalRelquadBand
  

  !****************************************************************************************
  !> Purpose: Return projected frequency map with respect to current component, 
  !>          alm_out = M^t*B^t*map_in, where B denotes beam convolution and M 
  !>          denotes the mixing matrix. This is the transpose of the above function
  ! 
  ! Inputs:
  !> @param[in]     self     -- current object
  !> @param[in]     band     -- data set ID, counting only active sets
  !> @param[in]     map      -- input map, with size (1:npix, nmaps). If not present, use
  !                             the parameters stored in self%T%alm
  !> @param[in]     alm_in   -- if true, the 'map' variable actually contains alms
  !> @param[in]     det      -- detector index; if not present, assume band-averaged output
  ! 
  ! Outputs:
  !> @return        projectRelquadBand -- alms = M^t*B^t*map_in
  !****************************************************************************************
  function projectRelquadBand(self, band, map, alm_in, det)
    implicit none
    class(comm_cmb_relquad_comp),                 intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: projectRelquadBand

    integer(i4b) :: i, nmaps, d
    logical(lgt) :: alm_in_
    class(comm_mapinfo), pointer :: info_in => null(), info_out => null()
    class(comm_map),     pointer :: m       => null(), m_out    => null()

    if (.not. self%band_active(band)) then
       if (.not. allocated(projectRelquadBand)) allocate(projectRelquadBand(0:self%info%nalm-1,self%info%nmaps))
       projectRelquadBand = 0.d0
       return
    end if

    d = 0; if (present(det)) d = det

    nmaps     =  min(self%nmaps, map%info%nmaps)
    info_in   => comm_mapinfo(self%info%comm, map%info%nside, map%info%lmax, nmaps, nmaps==3)
    info_out  => comm_mapinfo(self%info%comm, map%info%nside, self%lmax, nmaps, nmaps==3)
    m         => comm_map(info_in)
    m_out     => comm_map(info_out)
    alm_in_ = .false.; if (present(alm_in)) alm_in_ = alm_in

    ! Convolve with band-specific beam
    if (alm_in) then
       m%alm = map%alm(:,1:nmaps)
    else
       m%map = map%map(:,1:nmaps)
       call m%Yt()
    end if
    call data(band)%B(d)%p%conv(trans=.true., map=m)
    
    ! Apply mixing matrix
    m%alm(:,1) = m%alm(:,1) * self%q(band,d)

    ! Output map
    if (.not. allocated(projectRelquadBand)) allocate(projectRelquadBand(0:self%info%nalm-1,self%info%nmaps))
    call m%alm_equal(m_out)
    projectRelquadBand = m_out%alm

    call m%dealloc(); deallocate(m)
    call m_out%dealloc(); deallocate(m_out)

  end function projectRelquadBand


  !****************************************************************************************
  !> Purpose: Update spatial relativistic quadrupole template; if input variables are present,
  !>          set object variables equal to input. If not, use internal variables
  ! 
  ! Inputs:
  !> @param[in]     self     -- current object
  !> @param[in]     v        -- solar dipole velocity in m/s; optional
  !> @param[in]     lon      -- solar dipole longitude in degrees; optional
  !> @param[in]     lat      -- solar dipole latitude in degrees; optional
  ! 
  ! Outputs:
  !                 Internal object variables, self%T%{map,v,lon,lat}
  !****************************************************************************************
  subroutine update_template(self, v, lon, lat)
    implicit none
    class(comm_cmb_relquad_comp),    intent(inout)          :: self
    real(dp),                        intent(in),   optional :: v, lon, lat

    integer(i4b) :: i, j
    real(dp)     :: vec0(3), vec(3), theta, phi

    if (present(v))   self%v_sun   = v
    if (present(lon)) self%lon_sun = lon
    if (present(lat)) self%lat_sun = lat
    
    ! Set up velocity vector, beta, in units of v/c
    theta = (90.d0-self%lat_sun) * pi/180.d0
    phi   = self%lon_sun         * pi/180.d0
    call ang2vec(theta, phi, vec0)

    ! Construct template
    do i = 0, self%info%np-1
       call pix2vec_ring(self%info%nside, self%info%pix(i+1), vec)
       self%T%map(i,1) = dot_product(vec0, vec)**2 - 1.d0/3.d0
    end do
    
    ! Compute spherical harmonics up to lmax=2
    call self%T%YtW()

  end subroutine update_template

  !****************************************************************************************
  !> Purpose: Update band integration lookup table (self%F_int) and mixing matrix (self%q).
  !>          Done during initialization or when bandpasses are modified.
  ! 
  ! Inputs:
  !> @param[in]     self     -- current object
  !> @param[in]     band     -- frequency band index; optional
  !
  ! Outputs:
  !                 Internal object variables, self%{F_int,q}
  !****************************************************************************************
  subroutine updateIntF_relquad(self, band)
    implicit none
    class(comm_cmb_relquad_comp),            intent(inout)        :: self
    integer(i4b),                            intent(in), optional :: band

    integer(i4b) :: i, j, first_band, last_band
    real(dp)     :: f

    if (present(band)) then
       first_band = band
       last_band  = band
    else
       first_band = 1
       last_band  = numband
    end if

    do i = first_band, last_band
       if (.not. self%band_active(i)) then
          self%q(i,:) = 0.d0
          cycle
       end if
       do j = 0, data(i)%ndet
          self%F_int(i,j)%p => comm_F_int_0D(self, data(i)%bp(j)%p, 1)
          self%q(i,j)       = self%F_int(i,j)%p%eval([0.d0]) * data(i)%gain * self%cg_scale
       end do
    end do
    
  end subroutine updateIntF_relquad


  !****************************************************************************************
  !> Purpose: Read component definition file. 
  ! 
  !> File format: Each line contains {frequency label, .true./.false.}, e.g.,
  !>          030   .true.    # Quadrupole component is active
  !>          044   .true.    # Quadrupole component is active
  !>          Ka    .false.   # Quadrupole component is inactive, ie., assumed already subtracted
  ! 
  ! Inputs:
  !> @param[in]     self     -- current object
  !> @param[in]     filename -- input filename with path relative to workdir
  !
  ! Outputs:
  !     Internal object variables, self%band_active
  !****************************************************************************************
  subroutine read_definition_file(self, filename)
    implicit none
    class(comm_cmb_relquad_comp),  intent(inout)  :: self
    character(len=*),              intent(in)     :: filename

    integer(i4b)      :: i, unit
    logical(lgt)      :: active
    character(len=80) :: line, label
    logical(lgt), allocatable, dimension(:) :: ok

    allocate(ok(numband))
    ok = .false.

    unit  = getlun()
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       read(line,*) label, active
       do i = 1, numband
          if (trim(label) == trim(data(i)%label)) then
             self%band_active(i) = active
             ok(i)               = .true.
             exit
          end if
       end do
    end do
1   close(unit)
    
    if (any(.not. ok)) then
       write(*,*) 'ERROR -- all bands must be defined in CMB relativistic quadrupole file = ', trim(filename)
       stop
    end if

    deallocate(ok)

  end subroutine read_definition_file

  ! Dump current sample to HEALPix FITS file
  subroutine initTemplateRelQuadHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_cmb_relquad_comp), intent(inout) :: self
    type(comm_params),            intent(in)    :: cpar    
    type(hdf_file),               intent(in)    :: hdffile
    character(len=*),             intent(in)    :: hdfpath
  end subroutine initTemplateRelQuadHDF


end module comm_cmb_relquad_comp_mod
