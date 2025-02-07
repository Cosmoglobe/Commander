module comm_zodi_mod
   use comm_tod_mod
   use comm_zodi_comp_mod
   use comm_bp_mod
   implicit none

   private
   public initialize_zodi_mod, get_s_zodi, zodi_model, get_zodi_emission, update_zodi_splines, output_tod_params_to_hd5, read_tod_zodi_params, get_zodi_emissivity_albedo
   public get_s_tot_zodi, ZodiModel, zodi_model_to_ascii, ascii_to_zodi_model, print_zodi_model
   public band_monopole, band_update_monopole

   type :: ZodiCompLOS
      real(dp) :: R_min, R_max
      real(dp), allocatable, dimension(:) :: gauss_nodes, gauss_weights
      real(dp), allocatable, dimension(:) :: R, T, n, F_sol, Theta, Phi, B_nu
      real(dp), allocatable, dimension(:, :) :: X, X_unit
   end type

   type :: ZodiModel
      class(ZodiComponentContainer), allocatable :: comps(:)
      character(len=24), allocatable :: comp_labels(:), general_labels(:), par_labels(:), par_labels_full(:)
      integer(i4b) :: n_comps, n_params, n_common_params, n_general_params
      logical(lgt) :: joint_mono
      real(dp)     :: min_solar_elong, max_solar_elong
      real(dp)     :: T_0, delta
      real(dp), dimension(10) :: F_sun = [2.3405606d8, 1.2309874d8, 64292872d0, 35733824d0, 5763843d0, 1327989.4d0, 230553.73d0, 82999.336d0, 42346.605d0, 14409.608d0] * 1d-20 ! convert to specific intensity units
      real(dp), dimension(10) :: C0 = [-0.94209999, -0.52670002, -0.4312, -0.4312, 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C1 = [0.1214, 0.18719999, 0.1715, 0.1715, 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C2 = [-0.1648, -0.59829998, -0.63330001, -0.63330001, 0., 0., 0., 0., 0., 0.]

      integer(i4b) :: npar_tot
      integer(i4b), allocatable, dimension(:,:) :: theta_stat
      logical(lgt), allocatable, dimension(:,:) :: sampgroup_active_band
      integer(i4b), allocatable, dimension(:)   :: theta2band
      real(dp),     allocatable, dimension(:,:) :: theta_prior
      real(dp),     allocatable, dimension(:,:) :: theta_scale

      ! Stationary model
!      real(dp), allocatable, dimension(:)   :: amp_static
!      real(dp), allocatable, dimension(:,:) :: map_static
    contains
      procedure :: init_general_params, model_to_chain, params_to_model2, model_to_params2, comp_from_chain, get_par_ind, init_general_priors_and_scales
   end type ZodiModel

   type(ZodiModel), target :: zodi_model
   integer(i4b)            :: numband
   character(len=128)      :: zodi_refband
   character(len=128), allocatable, dimension(:) :: band_labels
   character(len=128), allocatable, dimension(:) :: band_instlabels
   character(len=128), allocatable, dimension(:) :: band_todtype
   real(dp),           allocatable, dimension(:) :: band_nu_c
   real(dp),           allocatable, dimension(:) :: band_monopole
   logical(lgt),       allocatable, dimension(:,:) :: band_update_monopole  ! (numband,0:n_sampgroup)
   
   real(dp) :: R_MIN = 3.d-14, R_CUTOFF = 5.2, EPS = TINY(1.0_dp), delta_t_reset
   real(dp), allocatable :: T_grid(:), B_nu_integrals(:)
   type(ZodiCompLOS), allocatable, dimension(:) :: comp_LOS
   
contains
  subroutine initialize_zodi_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar
    integer(i4b) :: i, j, k, ierr, ind, npar, ntok, gauss_degree
    real(dp) :: min_temp = 40.0, max_temp = 550.0
    integer(i4b) :: n_interp_points = 100
    character(len=256) :: file_path
    real(dp), allocatable :: comp_params(:, :)
    character(len=128), allocatable :: comp_labels(:)
    character(len=128) :: tokens(100)
    
    

    end subroutine initialize_zodi_mod

    
   subroutine init_general_params(self, general_params)
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: general_params(:)
      self%T_0 = general_params(1)
      self%delta = general_params(2)
   end subroutine

    subroutine init_general_priors_and_scales(self, prior, scale)
      implicit none
      class(ZodiModel),           intent(in)    :: self
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      prior(:,1) = [250.d0, 300.d0, 286.d0, 5.d0] ! T_0
      scale(1,:) = [286.d0, 3.d0]
      prior(:,2) = [0.4d0, 0.5d0, 0.467d0, 0.004d0] ! delta
      scale(2,:) = [0.4d0, 0.01d0]      
    end subroutine init_general_priors_and_scales






    subroutine model_to_params2(self, x, labels)
      ! Dumps a zodi model to a parameter vector `x`. If `labels` is present, it is populated with
      ! the corresponding parameter labels.
      class(ZodiModel), intent(in) :: self
      real(dp), intent(out) :: x(:)
      character(len=*), allocatable, optional, intent(inout) :: labels(:)
      character(len=128), allocatable :: labels_copy(:), comp_label_upper(:)
      integer(i4b) :: i, j, running_idx

      if (size(x) /= self%n_params) stop "Error: argument 'x' has the wrong size. must be `size(zodi_model%n_params)`"
      if (present(labels)) then
         if (allocated(labels)) stop "`labels` must not be allocated at the time of passing it in to `model_to_params`"
         allocate(comp_label_upper(self%n_comps))
      end if

      running_idx = 0
      do i = 1, self%n_comps
         x(running_idx + 1) = self%comps(i)%c%n_0
         x(running_idx + 2) = self%comps(i)%c%incl
         x(running_idx + 3) = self%comps(i)%c%Omega
         x(running_idx + 4) = self%comps(i)%c%x_0
         x(running_idx + 5) = self%comps(i)%c%y_0
         x(running_idx + 6) = self%comps(i)%c%z_0
         running_idx = running_idx + self%n_common_params
         select type (comp => self%comps(i)%c)
         class is (ZodiCloud)
            x(running_idx + 1) = comp%alpha
            x(running_idx + 2) = comp%beta
            x(running_idx + 3) = comp%gamma
            x(running_idx + 4) = comp%mu
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiBand)
            x(running_idx + 1) = comp%delta_zeta
            x(running_idx + 2) = comp%delta_r
            x(running_idx + 3) = comp%v
            x(running_idx + 4) = comp%p
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiRing)
            x(running_idx + 1) = comp%R_0
            x(running_idx + 2) = comp%sigma_r
            x(running_idx + 3) = comp%sigma_z
            x(running_idx + 4) = comp%theta_0
            x(running_idx + 5) = comp%sigma_theta
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            x(running_idx + 1) = comp%R_0
            x(running_idx + 2) = comp%sigma_r
            x(running_idx + 3) = comp%sigma_z
            x(running_idx + 4) = comp%theta_0
            x(running_idx + 5) = comp%sigma_theta
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiInterstellar)
            x(running_idx + 1) = comp%R
            x(running_idx + 2) = comp%alpha
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFan)
            x(running_idx + 1) = comp%Q
            x(running_idx + 2) = comp%P
            x(running_idx + 3) = comp%gamma
            x(running_idx + 4) = comp%Z_midplane_0
            x(running_idx + 5) = comp%R_outer
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiComet)
            x(running_idx + 1) = comp%P
            x(running_idx + 2) = comp%Z_midplane_0
            x(running_idx + 3) = comp%R_inner
            x(running_idx + 4) = comp%R_outer
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiWrightCloudRing)
            x(running_idx + 1)  = comp%p1
            x(running_idx + 2)  = comp%p3
            x(running_idx + 3)  = comp%p4
            x(running_idx + 4)  = comp%p5
            x(running_idx + 5)  = comp%p6
            x(running_idx + 6)  = comp%p7
            x(running_idx + 7)  = comp%p8
            x(running_idx + 8)  = comp%p9
            x(running_idx + 9)  = comp%p10
            x(running_idx + 10) = comp%p13
            x(running_idx + 11) = comp%p14
            x(running_idx + 12) = comp%p15
            x(running_idx + 13) = comp%p16
            x(running_idx + 14) = comp%p17
            x(running_idx + 15) = comp%p18
            x(running_idx + 16) = comp%p19
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiWrightBand)
            x(running_idx + 1)  = comp%q1
            x(running_idx + 2)  = comp%q5
            x(running_idx + 3)  = comp%q6
            x(running_idx + 4)  = comp%q7
            x(running_idx + 5)  = comp%q8
            x(running_idx + 6)  = comp%R_1
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         if (present(labels)) then
            labels_copy = self%comps(i)%labels
            comp_label_upper(i) = self%comp_labels(i)
            call toupper(comp_label_upper(i))
            do j = 1, size(labels_copy)
               labels_copy(j) = trim(adjustl(comp_label_upper(i)))//'_'//trim(adjustl(labels_copy(j))) 
            end do
               labels = [labels, labels_copy]
         end if
      end do
      x(running_idx + 1) = self%T_0
      x(running_idx + 2) = self%delta
      if (present(labels)) then
         labels = [labels, self%general_labels]
      end if
    end subroutine model_to_params2

   subroutine params_to_model2(self, x)
      ! Dumps a zodi model to a parameter vector `x`.
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      integer(i4b) :: i, running_idx

      if (size(x) /= self%n_params) then
         write(*,*) "Error: argument 'x' has the wrong size. must be `size(zodi_model%n_params)`", size(x), self%n_params
         stop
      end if

      running_idx = 0
      do i = 1, self%n_comps
         self%comps(i)%c%n_0 = x(running_idx + 1)
         self%comps(i)%c%incl = mod(x(running_idx + 2), 360.) ! degree prior
         self%comps(i)%c%Omega = mod(x(running_idx + 3), 360.) ! degree prior
         self%comps(i)%c%x_0 = x(running_idx + 4)
         self%comps(i)%c%y_0 = x(running_idx + 5)
         self%comps(i)%c%z_0 = x(running_idx + 6)
         running_idx = running_idx + self%n_common_params
         
         ! The order of these operations much match the order tabulated in the labels in `InterplanetaryDustParamLabels`
         select type (comp => self%comps(i)%c)
         class is (ZodiCloud)
            comp%alpha = x(running_idx + 1)
            comp%beta = x(running_idx + 2)
            comp%gamma = x(running_idx + 3)
            comp%mu = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiBand)
            comp%delta_zeta = mod(x(running_idx + 1), 360.) ! degree prior
            comp%delta_r = x(running_idx + 2)
            comp%v = x(running_idx + 3)
            comp%p = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiRing)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            comp%theta_0 = mod(x(running_idx + 4), 360.) ! degree prior
            comp%sigma_theta = mod(x(running_idx + 5), 360.) ! degree prior
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            comp%theta_0 = mod(x(running_idx + 4), 360.) ! degree prior
            comp%sigma_theta = mod(x(running_idx + 5), 360.) ! degree prior
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiInterstellar)
            comp%R = x(running_idx + 1)
            comp%alpha = x(running_idx + 2)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFan)
            comp%Q = x(running_idx + 1)
            comp%P = x(running_idx + 2)
            comp%gamma = x(running_idx + 3)
            comp%Z_midplane_0 = x(running_idx + 4)
            comp%R_outer = x(running_idx + 5)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiComet)
            comp%P = x(running_idx + 1)
            comp%Z_midplane_0 = x(running_idx + 2)
            comp%R_inner = x(running_idx + 3)
            comp%R_outer = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiWrightCloudRing)
            comp%p1  = x(running_idx + 1)
            comp%p3  = x(running_idx + 2)
            comp%p4  = x(running_idx + 3)
            comp%p5  = x(running_idx + 4)
            comp%p6  = x(running_idx + 5)
            comp%p7  = x(running_idx + 6)
            comp%p8  = x(running_idx + 7)
            comp%p9  = x(running_idx + 8)
            comp%p10 = x(running_idx + 9)
            comp%p13 = x(running_idx + 10)
            comp%p14 = x(running_idx + 11)
            comp%p15 = x(running_idx + 12)
            comp%p16 = x(running_idx + 13)
            comp%p17 = x(running_idx + 14)
            comp%p18 = x(running_idx + 15)
            comp%p19 = x(running_idx + 16)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiWrightBand)
            comp%q1  = x(running_idx + 1)
            comp%q5  = x(running_idx + 2)
            comp%q6  = x(running_idx + 3)
            comp%q7  = x(running_idx + 4)
            comp%q8  = x(running_idx + 5)
            comp%R_1 = x(running_idx + 6)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         call self%comps(i)%c%init()
      end do
      self%T_0 = x(running_idx + 1)
      self%delta = x(running_idx + 2)
    end subroutine params_to_model2

   subroutine model_to_chain(self, cpar, iter)
      ! Dumps the zodi model to the chain file
      class(ZodiModel), intent(in) :: self
      type(comm_params), intent(in) :: cpar
      integer(i4b), intent(in) :: iter

      integer(i4b) :: i, j, hdferr, ierr, unit, n_comp_params, param_idx
      logical(lgt) :: exist, init, new_header
      character(len=6) :: itext
      character(len=4) :: ctext
      character(len=512) :: zodi_path, comp_path, comp_group_path, param_path, chainfile, hdfpath, param_label, general_group_path, path
      character(len=32), allocatable :: labels(:)
      real(dp), allocatable :: params(:)
      type(hdf_file) :: file

      if (.not. cpar%myid_chain == 0) return

      call int2string(cpar%mychain, ctext)
      chainfile = trim(adjustl(cpar%outdir))//'/chain'// &
          & '_c'//trim(adjustl(ctext))//'.h5'

      inquire (file=trim(chainfile), exist=exist)
      call open_hdf_file(chainfile, file, 'b')

      call int2string(iter, itext)
      zodi_path = trim(adjustl(itext))//'/zodi'
      call create_hdf_group(file, trim(adjustl(zodi_path)))


      general_group_path = trim(adjustl(zodi_path))//'/general'
      call create_hdf_group(file, trim(adjustl(general_group_path)))

      comp_group_path = trim(adjustl(zodi_path))//'/comps'
      call create_hdf_group(file, trim(adjustl(comp_group_path)))

      allocate(params(self%n_params))
      call self%model_to_params2(params, labels)
      
      param_idx = 0 
      do i = 1, self%n_comps
         comp_path = trim(adjustl(comp_group_path))//'/'//trim(adjustl(self%comp_labels(i)))//'/'
         call create_hdf_group(file, trim(adjustl(comp_path)))
         n_comp_params = size(self%comps(i)%labels)
         do j = 1, n_comp_params
            param_label = trim(adjustl(comp_path))//'/'//trim(adjustl(self%comps(i)%labels(j)))
            call write_hdf(file, trim(adjustl(param_label)), params(param_idx + j))
         end do
         param_idx = param_idx + n_comp_params
         call write_hdf(file, trim(adjustl(comp_path))//'/emissivity', self%comps(i)%c%emissivity)
         call write_hdf(file, trim(adjustl(comp_path))//'/albedo', self%comps(i)%c%albedo)
      end do
      if (param_idx + self%n_general_params /= self%n_params) stop "Error: param_idx + self%n_general_params /= self%n_params"
      do i = 1, self%n_general_params
         param_label = trim(adjustl(general_group_path))//'/'//trim(adjustl(self%general_labels(i)))
         call write_hdf(file, trim(adjustl(param_label)), params(param_idx + i))
      end do

      ! Static component
      !path = trim(adjustl(zodi_path))//'/static'
      !call create_hdf_group(file, trim(adjustl(path)))
      !call write_hdf(file, trim(adjustl(path))//'/map', self%map_static)
      !call write_hdf(file, trim(adjustl(path))//'/amp', self%amp_static)

      call close_hdf_file(file)
      deallocate(params,labels)
   end subroutine

   subroutine comp_from_chain(self, cpar, params, comp_idx)
      ! Initialize a component from a chain
      class(ZodiModel), target, intent(inout) :: self
      type(comm_params), intent(in) :: cpar
      real(dp), intent(inout) :: params(:, :)
      integer(i4b), intent(in) :: comp_idx

   end subroutine

   subroutine get_s_zodi(band, s_therm, s_scat, s_zodi, comp_id)
      ! Evaluates the zodiacal signal (eq. 20 in ZodiPy paper [k98 model]) given
      ! integrated thermal zodiacal emission and scattered zodiacal light.
      !
      ! Parameters:
      ! -----------
      ! s_scat :
      !     Integrated contribution from scattered sunlight light.
      ! s_therm :
      !     Integrated contribution from thermal interplanetary dust emission.
      ! s_zodi :
      !     Zodiacal signal.
      ! emissivity :
      !     Emissivity of the zodiacal components.
      ! albedo :
     !     Albedo of the zodiacal components.
     implicit none
     integer(i4b),                  intent(in)  :: band
     real(sp),     dimension(:, :), intent(in)  :: s_scat, s_therm
     real(sp),     dimension(:),    intent(out) :: s_zodi
     integer(i4b),                  intent(in), optional :: comp_id
     
     integer(i4b) :: i, first
     real(dp)     :: al, em

     first = 1; if (present(comp_id)) first = comp_id
     
     s_zodi = 0.
     do i = first, first+size(s_therm,2)-1
        al     = zodi_model%comps(i)%c%albedo(band)
        em     = zodi_model%comps(i)%c%emissivity(band)
        !write(*,*) i, em, al, any(s_scat(:,i)/=s_scat(:,i)), any(s_therm(:,i)/=s_therm(:,i))
        s_zodi = s_zodi + ((s_scat(:,i-first+1) * al) + (1. - al) * em * s_therm(:,i-first+1))
     end do
   end subroutine get_s_zodi

   
   function get_par_ind(self, comp, comp_str, param, em_band, al_band, em_string, al_string, mono_band, mono_string)
     implicit none
     class(ZodiModel),     intent(in)           :: self
     class(ZodiComponentContainer), intent(in), target, optional :: comp
     character(len=*),     intent(in), optional :: comp_str
     character(len=*),     intent(in), optional :: param
     integer(i4b),         intent(in), optional :: em_band, al_band, mono_band
     character(len=*),     intent(in), optional :: em_string, al_string, mono_string
     integer(i4b)                               :: get_par_ind
     
     integer(i4b) :: i, band
     class(ZodiComponentContainer), pointer :: c
     character(len=128) :: str1, str2

     get_par_ind = 0
     
     
   end function get_par_ind

   subroutine samp_group2stat(cpar, samp_group, active, stat)
     implicit none     
     type(comm_params), intent(in) :: cpar
     integer(i4b),                   intent(in)    :: samp_group
     logical(lgt),     dimension(:), intent(in)    :: active
     integer(i4b),     dimension(:), intent(inout) :: stat

     integer(i4b) :: i, c, j, k, first, last, n_params, n, m, ind, em_global, al_global, c_to, c_from, band
     character(len=128) :: tokens(100), comp_param(2), wire_from(2), wire_to(2), label, param_label_tokens(10), em_from(2), em_to(2)
     character(len=2048) :: str
     

   end subroutine samp_group2stat



   subroutine get_zodi_emission(tod, pix, scan, det, s_zodi_scat, s_zodi_therm, model, always_scattering, use_lowres_pointing, comp)
      ! Returns the predicted zodiacal emission for a scan (chunk of time-ordered data).
      !
      ! Parameters
      ! ----------
      ! tod : class(comm_tod)
      !     The TOD object holding the spline objects to update.
      ! pix : integer(i4b), dimension(ntod)
      !     The pixel indices of each time-ordered observation.
      ! det : integer(i4b)
      !     The detector index.
      ! scan : integer(i4b)
      !     The scan number.
      ! s_zodi_scat : real(sp), dimension(ntod, ncomps)
      !     Contribution from scattered sunlight light.
      ! s_zodi_therm : real(sp), dimension(ntod, ncomps)
      !     Contribution from thermal interplanetary dust emission.
      ! model : type(ZodiModel)
      !     The zodiacal emission model.
      ! always_scattering : logical(lgt), optional
      !     If present, this overrides the default behavior of only including scattering when the albedo is non-zero.
      ! use_lowres_pointing : logical(lgt), optional
      !     If present, the input pixels are converted to low resolution pixels before evaluating the zodiacal emission.
      ! comp : integer(i4b), optional
      !     If present, only evaluate the zodiacal emission for this component.
      !
      ! Returns
      ! -------
      ! s_zodi_scat : real(sp), dimension(ntod, ncomps, ndet)
      !     Contribution from scattered sunlight light.
      ! s_zodi_therm : real(sp), dimension(ntod, ncomps, ndet)
      !     Contribution from thermal interplanetary dust emission.
      !
      class(comm_tod), intent(inout) :: tod
      integer(i4b), intent(in) :: pix(:), scan, det
      real(sp), dimension(:, :), intent(inout) :: s_zodi_scat, s_zodi_therm
      type(ZodiModel), intent(in) :: model
      logical(lgt), intent(in), optional :: always_scattering, use_lowres_pointing
      integer(i4b), intent(in), optional :: comp

      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr, cache_hits
      logical(lgt) :: scattering, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, obs_time, phase_normalization, C0, C1, C2, lat, lon
      real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3)
      !real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS

      s_zodi_scat = 0.
      s_zodi_therm = 0.
      n_tod = size(pix, dim=1)

      dt_tod = (1./tod%samprate)
      obs_pos = tod%scans(scan)%x0_obs
      earth_pos = tod%scans(scan)%x0_earth
      R_obs = norm2(obs_pos)
      obs_time = tod%scans(scan)%t0(1)
      earth_lon = atan(earth_pos(2), earth_pos(1))

      C0 = zodi_model%C0(tod%zodiband)
      C1 = zodi_model%C1(tod%zodiband)
      C2 = zodi_model%C2(tod%zodiband)
      phase_normalization = get_phase_normalization(C0, C1, C2)
      if (present(always_scattering)) then
         scattering = always_scattering
      else
         scattering = .false.
         do i = 1, zodi_model%n_comps
            if (zodi_model%comps(i)%c%albedo(tod%id) > EPS) then
               scattering = .true.
               exit
            end if
         end do
      end if
      ! select the correct cache
      !use_lowres = .false.

      cache_hits = 0
!!$      open(58,file="zodi.dat",recl=1024)
      do i = 1, n_tod
         ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
         ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
         if (use_lowres) then
            obs_time = tod%scans(scan)%downsamp_obs_time(i)
         else
            obs_time = obs_time + dt_tod ! assumes a time continuous TOD
         end if
         if ((obs_time - tod%zodi_cache_time) >= delta_t_reset) then
            R_obs = norm2(obs_pos)
            earth_lon = atan(earth_pos(2), earth_pos(1))
            call tod%clear_zodi_cache(obs_time)
         end if

         ! Get lookup index for cache. If the pixel is already cached, used that value.
         if (use_lowres) then
            lookup_idx = tod%pix2ind_lowres(tod%udgrade_pix_zodi(pix(i)))
            unit_vector = tod%ind2vec_ecl_lowres(:, lookup_idx)
         else
            lookup_idx = tod%pix2ind(pix(i))
            !write(*,*) 'q1', tod%scanid(scan), lookup_idx
            unit_vector = tod%ind2vec_ecl(:, lookup_idx)
         end if

         do k = 1, model%n_comps
            ! If comp is present we only evaluate the zodi emission for that component.
            ! If comp == 0 then we evaluate the zodi emission for all components.
            if (present(comp)) then
               if (k /= comp .and. comp /= 0) cycle
            end if

            ! Get line of sight integration range
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_min, R_min)
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_max, R_max)

            do l = 1, 3
               ! Convert quadrature range from [-1, 1] to [R_min, R_max]
               comp_LOS(k)%X_unit(l, :) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
               comp_LOS(k)%X_unit(l, :) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
               comp_LOS(k)%X(l, :) = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
            end do
            comp_LOS(k)%R = norm2(comp_LOS(k)%X, dim=1)
!!$            do l = 1, size(comp_LOS(k)%R)
!!$               if (sum(comp_LOS(k)%X_unit(:,l)**2) == 0.d0) then
!!$                  write(*,*) 'comp', k, l
!!$                  write(*,*) 'R_MIN', R_max, R_min
!!$                  write(*,*) 'gauss', comp_LOS(k)%gauss_nodes(l)
!!$                  write(*,*) 'unit', unit_vector
!!$                  write(*,*) 'X_unit', comp_LOS(k)%X_unit(:,l)
!!$               end if
!!$            end do
            
            if (scattering) then
               comp_LOS(k)%F_sol = model%F_sun(tod%zodiband)/comp_LOS(k)%R**2
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call get_phase_function(comp_LOS(k)%Theta, C0, C1, C2, phase_normalization, comp_LOS(k)%Phi)
            end if

            call get_dust_grain_temperature(comp_LOS(k)%R, comp_LOS(k)%T, model%T_0, model%delta)
!!$            write(*,*) tod%info%myid, k, size(comp_LOS(k)%T), size(tod%zodi_B_nu_spl_obj(det)%x), size(tod%zodi_B_nu_spl_obj(det)%y)
!!$            write(*,*) tod%info%myid, k, comp_LOS(k)%T
!!$            write(*,*) tod%info%myid, k, tod%zodi_B_nu_spl_obj(det)%x
!!$            write(*,*) tod%info%myid, k, tod%zodi_B_nu_spl_obj(det)%y
!!$            write(*,*) tod%info%myid, k, tod%zodi_B_nu_spl_obj(det)%y2

            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            if (scattering) then
               s_zodi_scat(i, k) = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5*(R_max - R_MIN) * 1d20
            end if
            s_zodi_therm(i, k) = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights) * 0.5 * (R_max - R_MIN) * 1d20
         end do
!!$         call vec2ang(unit_vector, lat, lon)
!!$         write(58,*) i, lon*180.d0/pi, 90.d0-180.d0/pi*lat, 0.958*sum(s_zodi_therm(i,:)), sum(s_zodi_scat(i,:)), sum(s_zodi_therm(i,:))+sum(s_zodi_scat(i,:))
!!$
!!$         write(*,*) "X", comp_LOS(1)%X(1,:)
!!$         write(*,*) "Y", comp_LOS(1)%X(2,:)
!!$         write(*,*) "Z", comp_LOS(1)%X(3,:)
!!$         write(*,*) "R", comp_LOS(1)%R
!!$         !write(*,*) "scat", comp_LOS(1)%F_sol*comp_LOS(1)%Phi
!!$         write(*,*) "n", comp_LOS(1)%n
!!$         !write(*,*) "F_sun", model%F_sun(tod%zodiband)
!!$         !write(*,*) "F", comp_LOS(1)%F_sol*1.d20
!!$         !write(*,*) "Phi", comp_LOS(1)%Phi
!!$         !write(*,*) "s", comp_LOS(1)%F_sol*comp_LOS(1)%Phi*1d20 * 0.255d0 + (1.d0-0.255d0) * 1.d0 * comp_LOS(1)%B_nu* 1.d0
!!$         write(*,*) "T_0, delta", model%T_0, model%delta
!!$         write(*,*) "T", comp_LOS(1)%T
!!$         write(*,*) "s", comp_LOS(1)%B_nu*0.958
      end do

!!$      close(58)
!!$      call mpi_finalize(i)
!!$      stop

    end subroutine get_zodi_emission

   ! Functions for evaluating the zodiacal emission
   ! -----------------------------------------------------------------------------------
   subroutine get_sphere_intersection(unit_vector, obs_pos, R_obs, R_cutoff, R_intersection)
      ! Computes R_max (the length of the LOS such that it stops exactly at los_cutoff_radius).

      real(dp), intent(in), dimension(:) :: unit_vector, obs_pos
      real(dp), intent(in) :: R_obs, R_cutoff
      real(dp), intent(out) :: R_intersection
      real(dp) :: lon, lat, cos_lat, b, d, q

      if (R_obs > R_cutoff) then
         R_intersection = EPS
         return
      end if

      lon = atan(unit_vector(2), unit_vector(1))
      lat = asin(unit_vector(3))
      cos_lat = cos(lat)
      b = 2.*(obs_pos(1)*cos_lat*cos(lon) + obs_pos(2)*cos_lat*sin(lon))
      d = R_obs**2 - R_cutoff**2
      q = -0.5*b*(1.+sqrt(b**2 - (4.*d))/abs(b))
      R_intersection = max(q, d/q)
   end subroutine get_sphere_intersection

   subroutine get_dust_grain_temperature(R, T_out, T_0, delta)
      real(dp), dimension(:), intent(in) :: R
      real(dp), dimension(:), intent(out) :: T_out
      real(dp), intent(in) :: T_0, delta
      T_out = T_0*R**(-delta)
   end subroutine get_dust_grain_temperature

   subroutine get_blackbody_emission(nus, T, b_nu)
      real(dp), intent(in) :: nus(:), T
      real(dp), dimension(:), intent(out) :: b_nu
      integer(i4b) :: i
      real(dp) :: x
      do i = 1, size(nus)
         if (x < 0.001d0) then
            ! Use RJ approximation
            b_nu(i) = 2.d0
         else if (x > 50.d0) then
            ! Use Wien approximation
            b_nu(i) = 2.d0
         else
            ! Use exact expression
            b_nu(i) = 2.d0
         end if
      end do
      !b_nu = b_nu * 1d20 ! Convert from W/(m^2*sr*Hz) to MJy/sr
   end subroutine get_blackbody_emission

   subroutine get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, scattering_angle)
      real(dp), intent(in) :: X_helio_vec_LOS(:, :), X_vec_LOS(:, :), R_helio_LOS(:)
      real(dp), dimension(:), intent(out) :: scattering_angle
      real(dp), allocatable, dimension(:) :: cos_theta, R_LOS

      allocate(cos_theta(size(X_vec_LOS, dim=1)))
      allocate(R_LOS(size(X_vec_LOS, dim=1)))

      R_LOS = norm2(X_vec_LOS, dim=1)
      if (any(abs(R_LOS*R_helio_LOS) < 1e-6)) then
         write(*,*) 'Error in get_scattering_angle'
         write(*,*) 'X_vec_LOS = ', X_vec_LOS
         write(*,*) 'R_LOS = ', R_LOS
         write(*,*) 'helio = ', R_helio_LOS
      end if
      cos_theta = sum(X_helio_vec_LOS*X_vec_LOS, dim=1)/(R_LOS*R_helio_LOS)
      ! clip cos(theta) to [-1, 1]
      where (cos_theta > 1)
         cos_theta = 1
      elsewhere(cos_theta < -1)
         cos_theta = -1
      end where
      scattering_angle = acos(-cos_theta)
   end subroutine get_scattering_angle

   subroutine get_phase_function(Theta, C0 , C1 , C2, N, phase_function)
      real(dp), intent(in) :: Theta(:), C0, C1, C2, N
      real(dp), intent(out) :: phase_function(:)
      phase_function = N *  (C0 + (C1 * Theta) + exp(C2 * Theta))
   end subroutine

   function get_phase_normalization(C0, C1, C2) result(N)
      real(dp), intent(in) :: C0, C1, C2
      real(dp) :: term1, term2, term3, term4, N

      term1 = 2.*pi
      term2 = 2.*C0
      term3 = pi*C1
      term4 = (exp(C2 * pi) + 1.)/(C2**2 + 1.)
      N = 1. / (term1 * (term2 + term3 + term4))
   end function

   subroutine initialize_earth_pos_spline(cpar)
      ! Returns the spline object which is used to evaluate the earth position

      type(comm_params), intent(in) :: cpar

      integer :: i, n_earthpos, unit
      real(dp), allocatable :: tabulated_earth_time(:), tabulated_earth_pos(:, :)
      unit = getlun()
      open (unit, file=trim(trim(cpar%datadir)//'/'//trim("earth_pos_1980-2050_ephem_de432s.txt")))
      read (unit, *) n_earthpos
      read (unit, *) ! skip header

      allocate (tabulated_earth_pos(3, n_earthpos))
      allocate (tabulated_earth_time(n_earthpos))
      do i = 1, n_earthpos
         read (unit, *) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
      end do
      close (unit)
   end subroutine initialize_earth_pos_spline

   subroutine update_zodi_splines(tod, bandpass, det, model)
      ! Updates the spectral spline objects in the TOD object.
      !
      ! In the K98 model, several spectral parameters are tabulated at individual frequencies,
      ! which we need to evaluate over the bandpass. In a future version, we may want to fit
      ! a modified blackbody which would allow us to drop using some of these spline objects.
      !
      !  -----------------------------------------------------------------------------------------
      ! | The difficulty with this functino is that it needs the access to the bandpass, so is is |
      ! | very limited in where it can be excecuted in commander.                                 |
      !  -----------------------------------------------------------------------------------------
      !
      ! Parameters
      ! ----------
      ! tod : class(comm_tod)
      !     The TOD object holding the spline objects to update.
      ! bandpass : class(comm_bp_ptr)
      !   The bandpass object holding the bandpass frequencies, and the SED2F function
      !   (bandpass integration).
      ! det : integer(i4b)
      !   The detector to update the spline objects for.

      class(comm_tod), intent(inout) :: tod
      class(comm_bp_ptr), intent(in) :: bandpass
      integer(i4b), intent(in) :: det
      type(ZodiModel), intent(inout) :: model

      real(dp), allocatable :: b_nu(:)
      integer(i4b) :: i, j
      real(dp)     :: K, Inu0(1)

   end subroutine update_zodi_splines

   subroutine output_tod_params_to_hd5(cpar, model, iter)
     implicit none
      ! Writes the zodi model to an hdf file
      type(comm_params), intent(in) :: cpar
      !class(comm_tod), intent(inout) :: tod
      type(ZodiModel),   intent(in) :: model
      integer(i4b),      intent(in) :: iter

      integer(i4b) :: i, j, hdferr, ierr, unit
      logical(lgt) :: exist, init, new_header
      character(len=6) :: itext
      character(len=4) :: ctext
      character(len=512) :: zodi_path, tod_path, band_path, det_path, comp_path, chainfile, hdfpath, path
      character(len=10), allocatable :: param_names(:)
      real(dp), allocatable :: param_values(:)
      type(hdf_file) :: file
      type(h5o_info_t) :: object_info

      if (.not. cpar%myid_chain == 0) return

      call int2string(cpar%mychain, ctext)
      chainfile = trim(adjustl(cpar%outdir))//'/chain'// &
      & '_c'//trim(adjustl(ctext))//'.h5'

      inquire (file=trim(chainfile), exist=exist)
      call open_hdf_file(chainfile, file, 'b')

      call int2string(iter, itext)
      zodi_path = trim(adjustl(itext))//'/zodi'
      !tod_path = trim(adjustl(zodi_path))//'/tod'

      ! Dynamic components
      !call create_hdf_group(file, trim(adjustl(tod_path)))
      !band_path = trim(adjustl(tod_path))//'/'//trim(adjustl(tod%freq))
      !call create_hdf_group(file, trim(adjustl(band_path)))
      do i = 1, model%n_comps
         comp_path = trim(adjustl(zodi_path))//'/'//trim(adjustl(model%comp_labels(i)))//'/'
         call create_hdf_group(file, trim(adjustl(comp_path)))
         call write_hdf(file, trim(adjustl(comp_path))//'/emissivity', model%comps(i)%c%emissivity)
         call write_hdf(file, trim(adjustl(comp_path))//'/albedo', model%comps(i)%c%albedo)
      end do

      ! Static component
      !path = trim(adjustl(zodi_path))//'/static'
      !call create_hdf_group(file, trim(adjustl(path)))
      !call write_hdf(file, trim(adjustl(path))//'/map', model%map_static)
      !call write_hdf(file, trim(adjustl(path))//'/amp', model%amp_static)

      call close_hdf_file(file)
   end subroutine

   subroutine read_tod_zodi_params(cpar, model, tod)
     implicit none
      type(comm_params), intent(in)  :: cpar
      type(ZodiModel), intent(inout) :: model
      class(comm_tod), intent(in)    :: tod

      logical(lgt) :: exist
      integer(i4b) :: i, l, comp, unit, ierr, initsamp
      character(len=6) :: itext, itext2

      type(hdf_file) :: file
      real(dp) :: lambda, lambda_min, lambda_max
      real(dp), allocatable :: emissivity(:), albedo(:)
      character(len=512) :: chainfile, emissivity_path, albedo_path, band_path, comp_path, tod_path, group_name

   end subroutine read_tod_zodi_params


   subroutine get_s_tot_zodi(zodi_model, tod, det, scan, s, pix_dynamic, pix_static, s_therm, s_scat)
      implicit none
      class(ZodiModel),                 intent(in)               :: zodi_model
      class(comm_tod),                  intent(inout)            :: tod
      integer(i4b),                     intent(in)               :: det, scan
      real(sp),         dimension(:),   intent(out)              :: s
      integer(i4b),     dimension(:,:), intent(in),     optional :: pix_dynamic,  pix_static
      real(sp),         dimension(:,:), intent(out),    optional :: s_therm, s_scat
      
      integer(i4b) :: i, j, h, ntod, nhorn, ncomp, band
      real(sp)     :: w
      real(sp),     allocatable, dimension(:)   :: s_zodi
      real(sp),     allocatable, dimension(:,:) :: s_scat_, s_therm_

      s = 0.


    end subroutine get_s_tot_zodi

    subroutine get_zodi_emissivity_albedo(zodi, band, emissivity, albedo)
     implicit none
     class(ZodiModel),               intent(in)            :: zodi
     integer(i4b),                   intent(in)            :: band
     real(dp),         dimension(:), intent(out), optional :: emissivity, albedo

     integer(i4b) :: i

     do i = 1, zodi%n_comps
        if (present(emissivity)) emissivity(i) = zodi%comps(i)%c%emissivity(band)
        if (present(albedo))     albedo(i)     = zodi%comps(i)%c%albedo(band)
     end do

   end subroutine get_zodi_emissivity_albedo

   subroutine zodi_model_to_ascii(cpar, model, filename, overwrite)
      ! Dumps the zodi model to an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(in) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      logical(lgt), intent(in), optional :: overwrite

   end subroutine

   subroutine ascii_to_zodi_model(cpar, model, filename)
      ! Reads in and updates the zodi model from an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(inout) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      type(hash_tbl_sll) :: htbl

      integer(i4b) :: i, j, io, io_status, ierr, n_comps
      logical(lgt) :: exists
      character(len=512) :: key, val, line
      character(len=128), allocatable :: labels(:)
      characteR(len=128) :: toks(100)
      characteR(len=512) :: concatenated_string
      real(dp), allocatable :: params(:)

   end subroutine


   subroutine print_zodi_model(theta, samp_group)
     implicit none
     real(dp),     allocatable, intent(in) :: theta(:)
     integer(i4b),              intent(in) :: samp_group
      
      integer(i4b) :: i, j, k, l, n, idx, n_cols, col, comp, n_params
      logical(lgt) :: newline
      
      
   end subroutine

   
   
 end module comm_zodi_mod
