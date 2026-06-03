module comm_zodi_mod
   use comm_tod_mod
   use comm_zodi_comp_mod
   use comm_bp_mod
   implicit none

   private
   public initialize_zodi_mod, zodi_model, get_zodi_emission, update_zodi_splines, output_tod_params_to_hd5, read_tod_zodi_params, get_zodi_emission_mbb
   public ZodiModel, zodi_model_to_ascii, ascii_to_zodi_model, print_zodi_model, compute_zodi
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
      character(len=24) :: phasefunc_type, bandpass_type
      integer(i4b) :: n_comps, n_common_params, n_general_params, n_phase_params, numband
      logical(lgt) :: joint_mono
      real(dp)     :: min_solar_elong, max_solar_elong
      real(dp)     :: nu_min_scatter, nu_max_thermal
      real(dp)     :: T_0, delta
      real(dp), dimension(10) :: F_sun = [2.3405606d8, 1.2309874d8, 64292872d0, 35733824d0, 5763843d0, 1327989.4d0, 230553.73d0, 82999.336d0, 42346.605d0, 14409.608d0] * 1d-20 ! convert to specific intensity units
      real(dp), allocatable, dimension(:) :: par_phase

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
      procedure :: init_comps, init_general_params, model_to_chain, params_to_model, model_to_params, comp_from_chain, get_par_ind, init_general_priors_and_scales, get_phase_function
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
   type(spline_type) :: earth_pos_spl_obj(3)

   real(dp), dimension(10) :: C0_K98 = [-0.94209999, -0.52670002, -0.4312, -0.4312, 0., 0., 0., 0., 0., 0.]
   real(dp), dimension(10) :: C1_K98 = [0.1214, 0.18719999, 0.1715, 0.1715, 0., 0., 0., 0., 0., 0.]
   real(dp), dimension(10) :: C2_K98 = [-0.1648, -0.59829998, -0.63330001, -0.63330001, 0., 0., 0., 0., 0., 0.]

   
contains
  subroutine initialize_zodi_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar
    integer(i4b) :: i, j, k, ierr, ind, npar, ntok, gauss_degree, band
    real(dp) :: min_temp = 40.0, max_temp = 550.0
    integer(i4b) :: n_interp_points = 100
    character(len=256) :: file_path
    real(dp), allocatable :: comp_params(:, :)
    character(len=128), allocatable :: comp_labels(:)
    character(len=128) :: tokens(100)
    
    
      ! Find number of bands and labels
    numband         = count(cpar%ds_active)
    zodi_model%numband = numband
    allocate(band_labels(numband),band_instlabels(numband),band_todtype(numband),band_nu_c(numband),band_monopole(numband),band_update_monopole(numband,0:cpar%zs_num_samp_groups))
    band_labels     = pack(cpar%ds_label, cpar%ds_active)
    band_instlabels = pack(cpar%ds_instlabel, cpar%ds_active)
    band_todtype    = pack(cpar%ds_tod_type, cpar%ds_active)
    band_nu_c       = pack(cpar%ds_nu_c, cpar%ds_active)
    zodi_refband    = cpar%zs_refband
      
    ! Set model and zodi_mod parameters from cpar
    zodi_model%nu_min_scatter   = cpar%zs_nu_min_scatter * 1d9
    zodi_model%nu_max_thermal   = cpar%zs_nu_max_thermal * 1d9
    zodi_model%phasefunc_type   = cpar%zs_phasefunc
    zodi_model%bandpass_type    = cpar%zs_bandpass
    zodi_model%n_comps          = cpar%zs_ncomps
    allocate(comp_params(zodi_model%n_comps, size(cpar%zs_comp_params, dim=1)))
    zodi_model%general_labels   = cpar%zodi_param_labels%general
    zodi_model%comp_labels      = cpar%zs_comp_labels(1:zodi_model%n_comps)
    zodi_model%n_common_params  = size(cpar%zodi_param_labels%common)
    zodi_model%joint_mono       = cpar%zs_joint_mono
    
    if (trim(zodi_model%phasefunc_type) == 'K98') then
       zodi_model%n_phase_params   = 3*numband
       zodi_model%n_general_params = 2 + zodi_model%n_phase_params
       allocate(zodi_model%par_phase(zodi_model%n_phase_params))
       do i  = 1, numband
          read(band_instlabels(i),*) band
          if (band >= 1 .and. band <= 10) then
             zodi_model%par_phase(band+0*numband) = C0_K98(band)
             zodi_model%par_phase(band+1*numband) = C1_K98(band)
             zodi_model%par_phase(band+2*numband) = C2_K98(band)
          else
             write(*,*) 'Invalid band label for K98 phase function = ', band_instlabels(i)
          end if
       end do
    else if (trim(zodi_model%phasefunc_type) == 'Wright') then
       zodi_model%n_phase_params   = 2
       zodi_model%n_general_params = 2 + zodi_model%n_phase_params
       allocate(zodi_model%par_phase(zodi_model%n_phase_params))
       zodi_model%par_phase(1) = -0.3133d0 ! p20
       zodi_model%par_phase(2) =  0.5749d0 ! p21
    else if (trim(zodi_model%phasefunc_type) == 'Hong') then
       zodi_model%n_general_params = size(cpar%zodi_param_labels%general)
       zodi_model%n_phase_params   = 5
       allocate(zodi_model%par_phase(zodi_model%n_phase_params))
       zodi_model%par_phase(1) =  0.700d0 ! g_1
       zodi_model%par_phase(2) = -0.200d0 ! g_2       
       zodi_model%par_phase(3) = -0.810d0 ! g_3
       zodi_model%par_phase(4) =  0.330d0 ! w_2
       zodi_model%par_phase(5) =  0.005d0 ! w_3       
    else
       write(*,*) 'Unsupported zodi phase function type:', trim(zodi_model%phasefunc_type)
       stop
    end if
      
    comp_params = cpar%zs_comp_params(:, :, 1)
    do i = 1, zodi_model%n_comps
       if (trim(adjustl(cpar%zs_init_hdf(i))) /= 'none') then 
          call zodi_model%comp_from_chain(cpar, comp_params, i)
       end if
    end do
    call zodi_model%init_comps(comp_params, cpar%zs_comp_types, cpar%zodi_param_labels)
    call zodi_model%init_general_params(cpar%zs_general_params(:, 1))

    ! Find total number of free parameters, and set up parameter mapping
    zodi_model%npar_tot = zodi_model%n_general_params + 2*zodi_model%n_comps*numband + numband
    zodi_model%comps(1)%start_ind = zodi_model%n_general_params + 1
    do i = 1, zodi_model%n_comps
       zodi_model%npar_tot = zodi_model%npar_tot + zodi_model%comps(i)%npar
       if (i < zodi_model%n_comps) then
          zodi_model%comps(i+1)%start_ind = zodi_model%comps(i)%start_ind + &
               & zodi_model%comps(i)%npar + 2*numband
       end if
    end do
      
    ! stat =  0  -> sample freely
    ! stat = -1  -> fix to input
    ! stat = -2  -> fix to zero
    ! stat = -3  -> fix to unity
    ! stat >  0  -> set equal to parameter stat
    allocate(zodi_model%theta_stat(zodi_model%npar_tot,0:cpar%zs_num_samp_groups))
    allocate(zodi_model%theta2band(zodi_model%npar_tot))
    allocate(zodi_model%theta_prior(4,zodi_model%npar_tot)) ! [min,max,mean,rms]
    allocate(zodi_model%theta_scale(zodi_model%npar_tot,2))
    allocate(zodi_model%par_labels(zodi_model%npar_tot))
    allocate(zodi_model%par_labels_full(zodi_model%npar_tot))
      
    ! Set up sampling groups
    allocate(zodi_model%sampgroup_active_band(numband,cpar%zs_num_samp_groups))
    zodi_model%sampgroup_active_band = .false.
    do i = 1, cpar%zs_num_samp_groups
       call get_tokens(cpar%zs_samp_group_bands(i), ',', tokens, ntok) 
       do j = 1, ntok
          k = get_string_index(band_labels, tokens(j))
          zodi_model%sampgroup_active_band(k,i) = .true. 
       end do
       call samp_group2stat(cpar, i, zodi_model%sampgroup_active_band(:,i), zodi_model%theta_stat(:,i))
    end do
    do i = 1, zodi_model%npar_tot
       zodi_model%theta_stat(i,0) = maxval(zodi_model%theta_stat(i,1:cpar%zs_num_samp_groups))
    end do
    do i = 1, numband
       band_update_monopole(i,0) = any(band_update_monopole(i,1:cpar%zs_num_samp_groups))
    end do

      ! Initialize parameter-band mapping
      zodi_model%theta2band(1:zodi_model%n_general_params) = 0 ! General parameters affect all bands
      do i = 1, zodi_model%n_comps
         ind  = zodi_model%comps(i)%start_ind
         npar = zodi_model%comps(i)%npar
         zodi_model%theta2band(ind:ind+npar-1) = 0 ! Shape parameters affect all bands
         do j = 1, numband
            zodi_model%theta2band(ind+npar+j-1) = j   ! Emissivity only affect band j
            zodi_model%theta2band(ind+npar+numband+j-1) = j ! The same for albedo
         end do
      end do

      ! Initialize priors and scale factors
      call zodi_model%init_general_priors_and_scales(zodi_model%theta_prior, &
           & zodi_model%theta_scale)

      if (trim(zodi_model%phasefunc_type) == 'Hong') then
         zodi_model%par_labels(1:zodi_model%n_general_params) = &
              & zodi_model%general_labels
         zodi_model%par_labels_full(1:zodi_model%n_general_params) = &
              & zodi_model%par_labels(1:zodi_model%n_general_params)
      else
         zodi_model%par_labels(1:2) = zodi_model%general_labels(1:2)
         zodi_model%par_labels_full(1:2) = zodi_model%par_labels(1:2)
         do i = 3, zodi_model%n_general_params
            zodi_model%par_labels(i)      = 'skip'
            zodi_model%par_labels_full(i) = 'skip'
         end do
      end if
      do i = 1, zodi_model%n_comps
         ! Shape parameters
         ind = zodi_model%comps(i)%start_ind
         call zodi_model%comps(i)%c%init_priors_and_scales(ind, &
              & zodi_model%theta_prior, zodi_model%theta_scale)
         zodi_model%par_labels(ind:ind+zodi_model%comps(i)%npar-1) = &
              & zodi_model%comps(i)%labels
         do j = ind, ind+zodi_model%comps(i)%npar-1
            zodi_model%par_labels_full(j) = &
              & trim(zodi_model%comp_labels(i))//':'//trim(zodi_model%par_labels(j))
         end do
            
         ! Emissivity and albedo
         ind = zodi_model%comps(i)%start_ind + zodi_model%comps(i)%npar-1
         do j = 1, numband
            zodi_model%theta_prior(:,ind+j) = [0.d0, 10.d0, 1.d0, -1.d0] ! Emissivity
            zodi_model%theta_scale(ind+j,:)   = [1.d0,0.1d0]
            zodi_model%par_labels(ind+j)    = 'em@'//trim(band_labels(j))
            zodi_model%par_labels_full(ind+j)  = trim(zodi_model%comp_labels(i))//':em@'//trim(band_labels(j))
            
            zodi_model%theta_prior(:,ind+numband+j) = [0.d0, 1.d0, 0.3d0, -1.d0] ! Albedo
            zodi_model%theta_scale(ind+numband+j,:)   = [1.d0, 0.01d0]
            zodi_model%par_labels(ind+numband+j)    = 'al@'//trim(band_labels(j))
            zodi_model%par_labels_full(ind+numband+j) = trim(zodi_model%comp_labels(i))//':al@'//trim(band_labels(j))

            ! Set DIRBE emissivity rms by hand
            if (trim(band_instlabels(j)) == '06') zodi_model%theta_scale(ind+j,:)   = [1.d0,0.01d0]
            if (trim(band_instlabels(j)) == '07') zodi_model%theta_scale(ind+j,:)   = [1.d0,0.02d0]
            if (trim(band_instlabels(j)) == '08') zodi_model%theta_scale(ind+j,:)   = [1.d0,0.03d0]
         end do
      end do
      ! Monopoles
      do j = 1, numband
         ind = zodi_model%npar_tot - numband + j
         zodi_model%theta_prior(:,ind) = [0.d0, 1d30, 0.d0, -1.d0] ! Priors
         zodi_model%theta_scale(ind,:) = [1.d0,0.01d0]
         zodi_model%par_labels(ind)    = 'm@'//trim(band_labels(j))
         zodi_model%par_labels_full(ind) = 'm@'//trim(band_labels(j))
      end do
      
      if (cpar%myid_chain == 0) then
         write(*,*) ' Total number of free zodi parameters = ', count(zodi_model%theta_stat(:,0)==0)
      end if

      allocate(comp_LOS(zodi_model%n_comps))
      allocate(B_nu_integrals(n_interp_points))

      ! hyper parameters
      delta_t_reset = cpar%zs_delta_t_reset

      ! Set up Gauss-Legendre quadrature
      do i = 1, zodi_model%n_comps
         gauss_degree = cpar%zs_los_steps(i)
         allocate (comp_LOS(i)%gauss_nodes(gauss_degree), comp_LOS(i)%gauss_weights(gauss_degree))
         if (cpar%zs_r_min(i) == 0.) then
            comp_LOS(i)%R_min = R_MIN
         else 
               comp_LOS(i)%R_min = cpar%zs_r_min(i)
         end if
         comp_LOS(i)%R_max = cpar%zs_r_max(i)
         call leggaus(gauss_degree, comp_LOS(i)%gauss_nodes, comp_LOS(i)%gauss_weights)
         allocate(comp_LOS(i)%X_unit(3, gauss_degree))
         allocate(comp_LOS(i)%X(3, gauss_degree))
         allocate(comp_LOS(i)%R(gauss_degree))
         allocate(comp_LOS(i)%T(gauss_degree))
         allocate(comp_LOS(i)%n(gauss_degree))
         allocate(comp_LOS(i)%F_sol(gauss_degree))
         allocate(comp_LOS(i)%Theta(gauss_degree))
         allocate(comp_LOS(i)%Phi(gauss_degree))
         allocate(comp_LOS(i)%B_nu(gauss_degree))
      end do

      ! Set up interpolation grid for evaluating temperature
      allocate (T_grid(n_interp_points))
      call linspace(min_temp, max_temp, T_grid)

      ! Read earth position from file and set up spline object
      call initialize_earth_pos_spline(cpar)

    end subroutine initialize_zodi_mod

    
   subroutine init_general_params(self, general_params)
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: general_params(:)
      self%T_0   = general_params(1)
      self%delta = general_params(2)
      if (self%phasefunc_type == 'Hong') then
         self%par_phase(1:5)   = general_params(3:7)
      end if
   end subroutine

    subroutine init_general_priors_and_scales(self, prior, scale)
      implicit none
      class(ZodiModel),           intent(in)    :: self
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      prior(:,1) = [250.d0, 330.d0, 286.d0, 5.d0] ! T_0
      scale(1,:) = [286.d0, 3.d0]
      prior(:,2) = [0.3d0, 0.5d0, 0.467d0, 0.004d0] ! delta
      scale(2,:) = [0.4d0, 0.01d0]
      prior(:,3) = [-1d0, 1d0, 0.70d0, 0.1d0] ! g_1
      scale(3,:) = [1d0, 0.01d0]
      prior(:,4) = [-1d0, 1d0, -0.20d0, 0.1d0] ! g_2
      scale(4,:) = [1d0, 0.01d0]
      prior(:,5) = [-1d0, 1d0, -0.81d0, 0.1d0] ! g_3
      scale(5,:) = [1d0, 0.01d0]
      prior(:,6) = [0d0, 1d0, 0.330d0, 0.1d0] ! w_2
      scale(6,:) = [1d0, 0.01d0]
      prior(:,7) = [0d0, 1d0, 0.005d0, 0.0001d0] ! w_3
      scale(7,:) = [1d0, 0.0005d0]

    end subroutine init_general_priors_and_scales



    subroutine init_comps(self, params, comp_types, param_labels)
      implicit none
      ! Initializes the components in the zodi model and computes the number of parameters in the model.
      class(ZodiModel), target, intent(inout) :: self
      real(dp), intent(in) :: params(:, :)
      character(len=*), intent(in) :: comp_types(:)
      class(InterplanetaryDustParamLabels), intent(in) :: param_labels
      integer(i4b) :: i, ierr
      allocate (self%comps(self%n_comps))
      ! NOTE: The order of the parameters in the `params` array must match the below order of readin
      do i = 1, self%n_comps
         !self%comps(i)%labels = [param_labels%common]
         select case (trim(adjustl(comp_types(i))))
         case ('cloud')
            allocate (ZodiCloud::self%comps(i)%c)
            self%comps(i)%c = ZodiCloud(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & sed_ampl=params(i, 7), &
                & sed_cutoff=params(i, 8), &
                & sed_b=params(i, 9), &
                & alpha=params(i, 10), &
                & beta=params(i, 11), &
                & gamma=params(i, 12), &
                & mu=params(i, 13) &
                &)
            allocate(self%comps(i)%labels(13))
            self%comps(i)%labels = [param_labels%common, param_labels%cloud]
         case ('band')
            allocate (ZodiBand::self%comps(i)%c)
            self%comps(i)%c = ZodiBand(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & sed_ampl=params(i, 7), &
                & sed_cutoff=params(i, 8), &
                & sed_b=params(i, 9), &
                & delta_zeta=params(i, 10), &
                & delta_r=params(i, 11), &
                & v=params(i, 12), &
                & p=params(i, 13) &
            &)
            allocate(self%comps(i)%labels(13))
            self%comps(i)%labels = [param_labels%common, param_labels%band]
         case ('ring')
            allocate (ZodiRing::self%comps(i)%c)
            self%comps(i)%c = ZodiRing(&
               & n_0=params(i, 1), &
               & incl=params(i, 2), &
               & Omega=params(i, 3), &
               & x_0=params(i, 4), &
               & y_0=params(i, 5), &
               & z_0=params(i, 6), &
               & sed_ampl=params(i, 7), &
               & sed_cutoff=params(i, 8), &
               & sed_b=params(i, 9), &
               & R_0=params(i, 10), &
               & sigma_r=params(i, 11), &
               & sigma_z=params(i, 12), &
               & theta_0=params(i, 13), &
               & sigma_theta=params(i, 14) &
            &)
            allocate(self%comps(i)%labels(14))
            self%comps(i)%labels = [param_labels%common, param_labels%ring]
         case ('feature')
            allocate (ZodiFeature::self%comps(i)%c)
            self%comps(i)%c = ZodiFeature(&
               & n_0=params(i, 1), &
               & incl=params(i, 2), &
               & Omega=params(i, 3), &
               & x_0=params(i, 4), &
               & y_0=params(i, 5), &
               & z_0=params(i, 6), &
               & sed_ampl=params(i, 7), &
               & sed_cutoff=params(i, 8), &
               & sed_b=params(i, 9), &
               & R_0=params(i, 10), &
               & sigma_r=params(i, 11), &
               & sigma_z=params(i, 12), &
               & theta_0=params(i, 13), &
               & sigma_theta=params(i, 14) &
               &)
            allocate(self%comps(i)%labels(14))
            self%comps(i)%labels = [param_labels%common, param_labels%feature]
         case ('interstellar')
            allocate (ZodiInterstellar::self%comps(i)%c)
            self%comps(i)%c = ZodiInterstellar(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & sed_ampl=params(i, 7), &
                 & sed_cutoff=params(i, 8), &
                 & sed_b=params(i, 9), &
                 & R=params(i, 10), &
                 & alpha=params(i, 11) &
                 &)
            allocate(self%comps(i)%labels(11))
            self%comps(i)%labels = [param_labels%common, param_labels%interstellar]
         case ('fan')
            allocate (ZodiFan::self%comps(i)%c)
            self%comps(i)%c = ZodiFan(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & sed_ampl=params(i, 7), &
                 & sed_cutoff=params(i, 8), &
                 & sed_b=params(i, 9), &
                 & Q=params(i, 10), &
                 & P=params(i, 11), &
                 & gamma=params(i, 12), &
                 & Z_midplane_0=params(i, 13), &
                 & R_outer=params(i, 14) &
            &)
            allocate(self%comps(i)%labels(14))
            self%comps(i)%labels = [param_labels%common, param_labels%fan]
         case ('comet')
            allocate (ZodiFan::self%comps(i)%c)
            self%comps(i)%c = ZodiComet(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & sed_ampl=params(i, 7), &
                 & sed_cutoff=params(i, 8), &
                 & sed_b=params(i, 9), &
                 & P=params(i, 10), &
                 & Z_midplane_0=params(i, 11), &
                 & R_inner=params(i, 12), &
                 & R_outer=params(i, 13) &
                 &)
            allocate(self%comps(i)%labels(13))
            self%comps(i)%labels = [param_labels%common, param_labels%comet]
         case ('wrightcloudring')
            allocate (ZodiWrightCloudRing::self%comps(i)%c)
            self%comps(i)%c = ZodiWrightCloudRing(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & sed_ampl=params(i, 7), &
                & sed_cutoff=params(i, 8), &
                & sed_b=params(i, 9), &
                & p1=params(i, 10), &
                & p3=params(i, 11), &
                & p4=params(i, 12), &
                & p5=params(i, 13), &
                & p6=params(i, 14), &
                & p7=params(i, 15), &
                & p8=params(i, 16), &
                & p9=params(i, 17), &
                & p10=params(i, 18), &
                & p13=params(i, 19), &
                & p14=params(i, 20), &
                & p15=params(i, 21), &
                & p16=params(i, 22), &
                & p17=params(i, 23), &
                & p18=params(i, 24), &
                & p19=params(i, 25) &
                &)
            allocate(self%comps(i)%labels(25))
            self%comps(i)%labels = [param_labels%common, param_labels%wrightcloudring]
         case ('wrightband')
            allocate (ZodiWrightBand::self%comps(i)%c)
            self%comps(i)%c = ZodiWrightBand(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & sed_ampl=params(i, 7), &
                & sed_cutoff=params(i, 8), &
                & sed_b=params(i, 9), &
                & q1=params(i, 10), &
                & q5=params(i, 11), &
                & q6=params(i, 12), &
                & q7=params(i, 13), &
                & q8=params(i, 14), &
                & R_1=params(i, 15) &
                &)
            allocate(self%comps(i)%labels(15))
            self%comps(i)%labels = [param_labels%common, param_labels%wrightband]
         case default
            print *, 'Invalid zodi component type in zodi `init_from_params`:', trim(adjustl(comp_types(i)))
            stop
         end select
         self%comps(i)%npar = size(self%comps(i)%labels)
         call self%comps(i)%c%init()

         ! Initialize emissivity and albedo
         allocate(self%comps(i)%c%emissivity(numband),self%comps(i)%c%albedo(numband))
         self%comps(i)%c%emissivity(numband) = 1.d0
         self%comps(i)%c%albedo(numband)     = 0.d0
      end do

    end subroutine init_comps


   subroutine model_to_params(self, x, samp_group, labels)
     implicit none
     class(ZodiModel),               intent(in)           :: self
     real(dp),         dimension(:), intent(out)          :: x
     integer(i4b),                   intent(in), optional :: samp_group
     character(len=*), allocatable, optional, intent(inout) :: labels(:)

     integer(i4b) :: i, j, idx
     real(dp), allocatable, dimension(:) :: z
     character(len=128), allocatable :: labels_copy(:), comp_label_upper(:)
     
     allocate(z(self%npar_tot))

     if (present(labels)) then
        if (allocated(labels)) stop "`labels` must not be allocated at the time of passing it in to `model_to_params`"
        allocate(comp_label_upper(self%n_comps))
     end if
     
     ! General parameters
     z(1) = self%T_0
     z(2) = self%delta
     if (trim(self%phasefunc_type) == 'Hong') then
        z(3:7) = self%par_phase(1:5)
     end if
     if (present(labels)) then
        if (trim(self%phasefunc_type) == 'Hong') then
           labels = self%general_labels
        else
           labels = self%general_labels(1:2)
           do i = 1, self%n_phase_params
              labels = [labels, "SKIP"]
           end do
        end if
     end if
     idx = self%n_general_params
     
     ! Component parameters
     do i = 1, self%n_comps
        ! Shape parameters
        call self%comps(i)%c%model2param(z(idx+1:idx+self%comps(i)%npar))
        idx = idx + self%comps(i)%npar

        ! Emissivity and albedo
        do j = 1, numband
           if (band_todtype(j) == "none") then
              z(idx        +j) = 0.d0
              z(idx+numband+j) = 0.d0
           else
              z(idx        +j) = self%comps(i)%c%emissivity(j)
              z(idx+numband+j) = self%comps(i)%c%albedo(j)
           end if
        end do
        idx = idx + 2*numband

         if (present(labels)) then
            labels_copy = self%comps(i)%labels
            comp_label_upper(i) = self%comp_labels(i)
            call toupper(comp_label_upper(i))
            do j = 1, size(labels_copy)
               labels_copy(j) = trim(adjustl(comp_label_upper(i)))//'_'//trim(adjustl(labels_copy(j))) 
            end do
            labels = [labels, labels_copy]
            ! Add placeholders for emissivity and albedo
            do j = 1, numband
               labels = [labels, ["skip","skip"]]
            end do
         end if
     end do

     ! Monopoles
     do i = 1, numband
        idx = idx+1
        if (band_todtype(i) /= "none") then
           z(idx) = band_monopole(i) !get_monopole_amp(data(i)%label)
           !write(*,*) 'get', i, trim(data(i)%label), z(idx)
        else
           z(idx) = 0.d0
        end if
        if (present(labels)) labels = [labels, "skip"]
     end do
     
     if (present(samp_group)) then
        x = pack(z, self%theta_stat(:,samp_group)==0)
        if (present(labels)) labels = pack(labels, self%theta_stat(:,samp_group)==0)
     else
        x = z
     end if
     deallocate(z)     
     
   end subroutine model_to_params

   subroutine params_to_model(self, x, samp_group)
     implicit none
     class(ZodiModel),               intent(inout)        :: self
     real(dp),         dimension(:), intent(in)           :: x
     integer(i4b),                   intent(in), optional :: samp_group

     integer(i4b) :: i, j, idx
     real(dp), allocatable, dimension(:) :: z, z_prev
     
     allocate(z(self%npar_tot))

     ! Initialize full parameter vector
     if (present(samp_group)) then
        allocate(z_prev(self%npar_tot))
        call self%model_to_params(z_prev)
        idx = 1
        do i = 1, self%npar_tot
           if (self%theta_stat(i,samp_group) == 0) then
              z(i) = x(idx)
              idx = idx+1
           else if (self%theta_stat(i,samp_group) == -1) then
              z(i) = z_prev(i)
           else if (self%theta_stat(i,samp_group) == -2) then
              z(i) = 0.d0
           else if (self%theta_stat(i,samp_group) == -3) then
              z(i) = 1.d0
           end if
        end do
        do i = 1, self%npar_tot
           if (self%theta_stat(i,samp_group) > 0) then
              z(i) = z(self%theta_stat(i,samp_group))
           end if
        end do
        deallocate(z_prev)
     else
        z = x
     end if
     
     ! General parameters
     self%T_0   = z(1) 
     self%delta = z(2)
     if (trim(self%phasefunc_type) == 'Hong') then
        self%par_phase(1:5) = z(3:7) 
     end if
     
     ! Component parameters
     idx = self%n_general_params
     do i = 1, self%n_comps
        ! Shape parameters
        call self%comps(i)%c%param2model(z(idx+1:idx+self%comps(i)%npar))
        call self%comps(i)%c%init()

        ! Emissivity and albedo
        idx = idx + self%comps(i)%npar
        do j = 1, numband
           if (band_todtype(j) /= "none") then
              self%comps(i)%c%emissivity(j) = z(idx+j)
              self%comps(i)%c%albedo(j)     = z(idx+numband+j)
           end if
        end do
        idx = idx + 2*numband
     end do

     ! Monopoles
     do i = 1, numband
        idx = idx+1
        if (band_todtype(i) /= "none") then
           band_monopole(i) = z(idx) 
        end if
     end do

     deallocate(z)
     
   end subroutine params_to_model




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

      allocate(params(self%npar_tot))
      call self%model_to_params(params, labels=labels)

      do i = 1, self%n_general_params
         if (trim(labels(i)) == 'skip' .or. trim(labels(i)) == 'SKIP') cycle
         param_label = trim(adjustl(general_group_path))//'/'//trim(adjustl(self%general_labels(i)))
         call write_hdf(file, trim(adjustl(param_label)), params(i))
      end do

      param_idx = self%n_general_params      
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

      call close_hdf_file(file)
      deallocate(params,labels)
   end subroutine
   
   subroutine comp_from_chain(self, cpar, params, comp_idx)
     ! Initialize a component from a chain
     class(ZodiModel), target, intent(inout) :: self
     type(comm_params), intent(in) :: cpar
     real(dp), intent(inout) :: params(:, :)
     integer(i4b), intent(in) :: comp_idx
     
     logical(lgt) :: exist
     integer(i4b) :: i, j, l, ierr, initsamp
     character(len=6) :: itext
     
     type(hdf_file) :: file
     
     character(len=32), allocatable :: param_labels(:)
     character(len=2048) :: chainfile, group_name
     
     if (cpar%myid == cpar%root) then
        if (trim(cpar%zs_init_hdf(comp_idx)) == 'default') then
           call get_chainfile_and_samp(trim(cpar%init_chain_prefixes(1)), chainfile, initsamp)
        else
           call get_chainfile_and_samp(trim(cpar%zs_init_hdf(comp_idx)), chainfile, initsamp)
        end if
        inquire (file=trim(chainfile), exist=exist)
        if (.not. exist) call report_error('Zodi init chain does not exist = '//trim(chainfile))
        l = len(trim(chainfile))
        if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
        
        call open_hdf_file(trim(chainfile), file, "r")
        
        call int2string(initsamp, itext)
        group_name = trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(comp_idx))))
        if (.not. hdf_group_exists(file, group_name)) then
           print *, "zodi component: ", trim(adjustl(cpar%zs_comp_labels(comp_idx))), "not found in chain:", trim(chainfile)
           stop
        end if
        
        param_labels = cpar%zodi_param_labels%get_labels(trim(adjustl(cpar%zs_comp_types(comp_idx))), add_common=.true.)
        do j = 1, size(param_labels)
           call read_hdf(file, trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(comp_idx)))// &
                & '/'//trim(adjustl(param_labels(j)))), params(comp_idx, j))
        end do
        deallocate(param_labels)
     end if
     ! call mpi_bcast(params, sum(shape(params)), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
     call mpi_bcast(params, size(params, dim=1) * size(params, dim=2), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
     
   end subroutine comp_from_chain
   
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
     
     if (present(comp) .or. present(comp_str)) then
        if (present(comp_str)) then
           i =  get_string_index(self%comp_labels, comp_str)
           c => zodi_model%comps(i)
        else
           c => comp
        end if

        if (present(param)) then
           get_par_ind = c%start_ind + get_string_index(c%labels, param)-1
        else if (present(em_band)) then
           get_par_ind = c%start_ind + c%npar + em_band - 1
        else if (present(al_band)) then
           get_par_ind = c%start_ind + c%npar + numband + al_band - 1
        else if (present(em_string)) then
           get_par_ind = c%start_ind + c%npar + get_string_index(band_labels, em_string) - 1
        else if (present(al_string)) then
           get_par_ind = c%start_ind + c%npar + numband + get_string_index(band_labels, al_string) - 1
        else
           write(*,*) 'get_par_ind error: Need parameter specification'
           stop
        end if
     else
        if (present(param)) then
           ! General parameter
           get_par_ind = get_string_index(self%general_labels, param)
        else if (present(mono_string) .or. present(mono_band)) then
           ! Monopoles are always last
           if (present(mono_string)) then
              band = get_string_index(band_labels, mono_string)
           else
              band = mono_band
           end if
           get_par_ind = self%npar_tot - numband + band
        else
           write(*,*) 'get_par_ind error: Need parameter specification'
           stop
        end if
     end if
     
   end function get_par_ind

   subroutine samp_group2stat(cpar, samp_group, active, stat)
     implicit none     
     type(comm_params), intent(in) :: cpar
     integer(i4b),                   intent(in)    :: samp_group
     logical(lgt),     dimension(:), intent(in)    :: active
     integer(i4b),     dimension(:), intent(inout) :: stat

     integer(i4b) :: i, c, j, k, first, last, n_params, n, m, ind, em_global, al_global, c_to, c_from, band
     character(len=128) :: tokens(100), comp_param(2), wire_from(2), wire_to(2), label, param_label_tokens(10), em_from(2), em_to(2), al_from(2), al_to(2)
     character(len=2048) :: str
     
     ! Default: Fix everything at input
     str  = cpar%zs_samp_groups(samp_group)
     stat = -1
     
     ! Parse user directives
     call get_tokens(str, ',', tokens, n_params) 

     em_global = 0; al_global = 0
     if (cpar%zs_em_global /= 'none') em_global = get_string_index(zodi_model%comp_labels, cpar%zs_em_global)
     if (cpar%zs_al_global /= 'none') al_global = get_string_index(zodi_model%comp_labels, cpar%zs_al_global)
     band_update_monopole(:,samp_group) = .false.
     do i = 1, n_params
        call get_tokens(tokens(i), ':', comp_param, num=n)
        if (n == 1) then
           label = comp_param(1)
           call get_tokens(label, '@', comp_param, num=n)
           if (n == 1) then
              ! General 
              ind = zodi_model%get_par_ind(param=comp_param(1))
              stat(ind) = 0
           else if (n == 2) then
              if (trim(comp_param(1)) == 'm') then
                 ! Monopole
                 band = get_string_index(band_labels, comp_param(2))
                 if (active(band) .and. band_todtype(band) /= 'none') then
                    !ind       = zodi_model%get_par_ind(mono_band=band)
                    band_update_monopole(band,samp_group) = .true.
                    ind = zodi_model%get_par_ind(mono_band=band)
                    stat(ind) = 0
                 end if
                 cycle
              else
                 write(*,*) 'Unsupported zodi parameter in samp_group2stat: ', trim(comp_param(1))
                 stop
              end if
           end if
        else if (n == 2) then
           if (trim(comp_param(1)) == 'em') then
              band = get_string_index(band_labels, comp_param(2))
              if (.not. active(band)) cycle
              do j = 1, zodi_model%n_comps
                 ind       = zodi_model%get_par_ind(comp=zodi_model%comps(j), em_string=comp_param(2))
                 stat(ind) = 0
              end do
              cycle
           else if (trim(comp_param(1)) == 'al') then
              band = get_string_index(band_labels, comp_param(2))
              if (.not. active(band)) cycle
              do j = 1, zodi_model%n_comps
                 ind       = zodi_model%get_par_ind(comp=zodi_model%comps(j), al_string=comp_param(2))
                 stat(ind) = 0
              end do
              cycle
           end if

           ! Component parameters
           label = comp_param(2)
           c     = get_string_index(zodi_model%comp_labels, comp_param(1))
           first = zodi_model%comps(c)%start_ind
           if (trim(label) == 'all') then
              last  = first + zodi_model%comps(c)%npar - 1
              stat(first:last) = 0 ! Activate all
              do j = 1, numband
                 if (.not. active(j)) cycle
                 stat(last+j)         = 0
                 stat(last+numband+j) = 0
              end do
              last = last + 2*numband
           else if (trim(label(1:2)) == 'em') then
              ! Emissivity
              call get_tokens(label, '@', comp_param, num=n)
              if (n == 1) then
                 first = first + zodi_model%comps(c)%npar
                 do j = 1, numband
                    if (active(j)) stat(first+j-1) = 0 
                 end do
              else
                 band = get_string_index(band_labels, comp_param(2))
                 if (.not. active(band)) cycle
                 ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_string=comp_param(2))
                 stat(ind) = 0
              end if
           else if (trim(label(1:2)) == 'al' .and. trim(label) /= 'alpha') then
              ! Albedo
              call get_tokens(label, '@', comp_param, num=n)
              if (n == 1) then
                 first = first + zodi_model%comps(c)%npar + numband
                 do j = 1, numband
                    if (active(j)) stat(first+j-1) = 0 
                 end do
              else
                 band = get_string_index(band_labels, comp_param(2))
                 if (.not. active(band)) cycle
                 ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_string=comp_param(2))
                 stat(ind) = 0
              end if
           else
              ! Shape parameter
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), param=comp_param(2))
              stat(ind) = 0
           end if
        else
           write(*,*) 'Invalid zodi samp group element = ', trim(tokens(i))
           stop
        end if
     end do

     ! Set up monopoles
!!$     do i = 1, numband
!!$        if (active(i) .and. band_todtype(i) /= 'none') then
!!$           ind = zodi_model%get_par_ind(mono_band=i)
!!$           stat(ind) = 0
!!$        end if
!!$     end do

     ! Apply explicit parameter wiring
     call get_tokens(cpar%zs_wiring, ',', tokens, n_params) 
     do i = 1, n_params
        if (trim(tokens(1)) == 'none') exit
        
        call get_tokens(tokens(i), '>', comp_param, num=n)
        if (n /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', tokens(i)
           stop
        end if

        call get_tokens(comp_param(1), ':', wire_from, num=n)
        call get_tokens(comp_param(2), ':', wire_to,   num=m)
        if (n /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(tokens(i))!, trim(comp_param(1)), trim(wire_from(1)), trim(wire_from(2))
           stop
        end if
        if (m /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(tokens(i))!, trim(comp_param(2)), trim(wire_to(1)), trim(wire_to(2))
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(comp_param(2))!, trim(wire_to(1)), trim(wire_to(2))
           stop
        end if

        if (trim(wire_from(2)(1:2)) == 'em') then
           call get_tokens(wire_from(2), '@', em_from, num=n)
           call get_tokens(wire_to(2), '@', em_to, num=n)
           if (n == 1) then
              ! Attach all bands
              do c = 1, numband
                 if (.not. active(c)) cycle
                 c_from = zodi_model%get_par_ind(comp_str=wire_from(1), em_band=c)
                 c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   em_band=c)
                 stat(c_from) = c_to
              end do
           else
              ! Attach specified band
              j = get_string_index(band_labels, em_from(2))
              k = get_string_index(band_labels, em_to(2))
              if (.not. active(j) .or. .not. active(k)) cycle
              c_from = zodi_model%get_par_ind(comp_str=wire_from(1), em_string=em_from(2))
              c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   em_string=em_to(2))
              stat(c_from) = c_to
           end if
        else if (trim(wire_from(2)(1:2)) == 'al') then
           call get_tokens(wire_from(2), '@', al_from, num=n)
           call get_tokens(wire_to(2), '@', al_to, num=n)
           if (n == 1) then
              ! Attach all bands
              do c = 1, numband
                 if (.not. active(c)) cycle
                 c_from = zodi_model%get_par_ind(comp_str=wire_from(1), al_band=c)
                 c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   al_band=c)
                 stat(c_from) = c_to
              end do
           else
              ! Attach specified band
              j = get_string_index(band_labels, al_from(2))
              k = get_string_index(band_labels, al_to(2))
              if (.not. active(j) .or. .not. active(k)) cycle
              c_from = zodi_model%get_par_ind(comp_str=wire_from(1), al_string=al_from(2))
              c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   al_string=al_to(2))
              stat(c_from) = c_to
           end if
        else
           write(*,*) 'Unsupported zodi wire = ', trim(tokens(i))
           stop
        end if
     end do

     ! Apply global directives
     if (em_global > 0) then
        do i = 1, zodi_model%n_comps
           if (i == em_global) cycle
           do j = 1, numband
              if (.not. active(j)) cycle
              ind   = zodi_model%get_par_ind(comp=zodi_model%comps(em_global), em_band=j)
              first = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(first) = ind
           end do
        end do
     end if

     if (al_global > 0) then
        do i = 1, zodi_model%n_comps
           if (i == al_global) cycle
           do j = 1, numband
              if (.not. active(j)) cycle
              ind   = zodi_model%get_par_ind(comp=zodi_model%comps(al_global), al_band=j)
              first = zodi_model%get_par_ind(comp=zodi_model%comps(i), al_band=j)
              stat(first) = ind
           end do
        end do
     end if

     ! Match emissivity and albedo for bands with identical instruments
     do i = 1, numband
        if (.not. active(i)) cycle
        do j = i+1, numband
           if (.not. active(j)) cycle
           if (trim(band_instlabels(j)) == trim(band_instlabels(i))) then
              do c = 1, zodi_model%n_comps
                 ind   = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_band=i)
                 first = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_band=j)
                 stat(first) = ind
                 
                 ind   = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_band=i)
                 first = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_band=j)
                 stat(first) = ind
              end do
           end if
        end do
     end do

     ! Apply absolute constraints
     do j = 1, numband
        if (band_todtype(j) == 'none') then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -2 ! Fix emissivity to zero
              stat(ind+numband) = -2 ! Fix albedo to zero
           end do
        end if

        if (band_todtype(j) /= 'none' .and. band_nu_c(j) > 1d14) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -3 ! Fix emissivity to unity
           end do
        end if

        if (band_todtype(j) /= 'none' .and. band_nu_c(j) < 50000d9) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), al_band=j)
              stat(ind)         = -2 ! Fix albedo to zero
           end do
        end if

        if (trim(zodi_refband) == trim(band_labels(j))) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -3 ! Fix emissivity to unity
           end do
        end if
     end do

     ! Short-cut multi-leg wires
     do i = 1, size(stat)
        if (stat(i) > 0) then
           j = stat(i)
           do while (stat(j) > 0)
              stat(i) = stat(j)
              j       = stat(i)
              if (j <= 0) exit
           end do
        end if
     end do

   end subroutine samp_group2stat



   subroutine get_zodi_emission(tod, pix, scan, det, model, s_zodi, use_lowres_pointing, comp)
     implicit none
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
      class(comm_tod),               intent(inout)        :: tod
      integer(i4b),    dimension(:), intent(in)           :: pix
      integer(i4b),                  intent(in)           :: scan, det
      type(ZodiModel),               intent(in)           :: model
      real(sp),        dimension(:), intent(out)          :: s_zodi
      logical(lgt),                  intent(in), optional :: use_lowres_pointing
      integer(i4b),                  intent(in), optional :: comp

      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr
      logical(lgt) :: scattering, thermal, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, obs_time, lat, lon, s_tot, s_scat, s_therm, al, em
      real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3)
      real(dp), allocatable, dimension(:)   :: b_nu

      integer(i4b) :: cache_hits
      real(sp), allocatable, dimension(:,:)     :: cache
      real(dp)                                  :: cache_time, init_cache_time
      
      use_lowres = .false.; if (present(use_lowres_pointing)) use_lowres = use_lowres_pointing
      n_tod      = size(pix)
      dt_tod     = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_pos    = tod%scans(scan)%x0_obs
      earth_pos  = tod%scans(scan)%x0_earth
      R_obs      = norm2(obs_pos)
      obs_time   = tod%scans(scan)%t0(1)
      earth_lon  = atan(earth_pos(2), earth_pos(1))
      scattering = tod%central_freq >= model%nu_min_scatter
      thermal    = tod%central_freq <= model%nu_max_thermal

      ! Initialize cache
      if (use_lowres) then
         !allocate(cache(tod%zodi_cache_nobs_lowres, tod%ndet))
      else
         allocate(cache(tod%pixcache%nobs,                   tod%ndet))
      end if
      cache_time = tod%scans(1)%t0(1)
      cache      = -1
      
      cache_hits = 0
!!$      open(58,file="zodi.dat",recl=1024)
      do i = 1, n_tod
         ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
         ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
         if (use_lowres) then
            obs_time = tod%scans(scan)%d(det)%downsamp_obs_time_full(i)
         else
            obs_time = obs_time + dt_tod ! assumes a time continuous TOD
         end if

         ! Reset cache
         if ((obs_time - cache_time) >= delta_t_reset) then
            do j = 1, 3
               earth_pos(j) = splint_simple(tod%x_earth_spline(j), obs_time)
               obs_pos(j)   = splint_simple(tod%x_obs_spline(j),   obs_time)
            end do
            R_obs      = norm2(obs_pos)
            earth_lon  = atan(earth_pos(2), earth_pos(1))
            cache      = -1
            cache_time =  obs_time
         end if

         ! Get lookup index for cache. If the pixel is already cached, used that value.
         if (use_lowres) then
            !lookup_idx  = tod%pix2ind_lowres(tod%udgrade_pix_zodi(pix(i))) 
            !unit_vector = tod%ind2vec_ecl_lowres(:, lookup_idx)
         else
            lookup_idx  = tod%pixcache%pix2ind(pix(i))
            if (lookup_idx == -2) write(*,*) 'oor', tod%scanid(scan), det, i, pix(i)
            unit_vector = tod%pixcache%ind2vec_ecl(:, lookup_idx)
         end if
         
         if (cache(lookup_idx, det) > 0.d0) then
            s_zodi(i)  = cache(lookup_idx, det)
            cache_hits = cache_hits + 1
            cycle
         end if

         ! If not present in cache, compute signal for each zodi component from scratch, and add signals together
         s_tot = 0.
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
               comp_LOS(k)%X_unit(l,:) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
               comp_LOS(k)%X_unit(l,:) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
               comp_LOS(k)%X(l,:)      = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
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

            ! Compute phase function if scattering is included
            if (scattering) then
               if (trim(model%phasefunc_type) == 'Wright') then
                  allocate(b_nu(size(tod%bp(0)%p%nu)))
                  call get_blackbody_emission(tod%bp(0)%p%nu, 5772.d0, b_nu) 
                  comp_LOS(k)%F_sol = tsum(tod%bp(0)%p%nu, tod%bp(0)%p%tau*b_nu)/comp_LOS(k)%R**2
                  deallocate(b_nu)
               else
                  comp_LOS(k)%F_sol = model%F_sun(tod%zodiband)/comp_LOS(k)%R**2
               end if
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call model%get_phase_function(comp_LOS(k)%Theta, tod%zodiband, comp_LOS(k)%Phi)
            end if

            ! Get dust grain temperature along the line of sight, and compute splined blackbody emission
            if (trim(model%phasefunc_type) == 'Wright' .and. k > 1) then
               !comp_LOS(k)%T = exp(5.5301d0) * comp_LOS(k)%R**(-0.5d0)
               comp_LOS(k)%T = exp(5.5301d0) / sqrt(comp_LOS(k)%R)
            else
               !comp_LOS(k)%T = model%T_0 * comp_LOS(k)%R**(-model%delta)
               comp_LOS(k)%T = model%T_0 * exp(-model%delta * log(comp_LOS(k)%R))
            end if
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

!!$            ! Compute density along line of sight
!!$            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            
            ! Compute predicted signal in MJy/sr; store computed signal in cache
            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            s_scat = 0.; s_therm = 0.
            if (scattering) s_scat  = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5d0*(R_max-R_MIN)*1d20
            if (thermal)    s_therm = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights)                  * 0.5d0*(R_max-R_MIN)*1d20

            al = zodi_model%comps(k)%c%albedo(tod%id)
            em = zodi_model%comps(k)%c%emissivity(tod%id)
            if (trim(zodi_model%phasefunc_type) == 'Wright') then
               s_tot = s_tot + al*s_scat + em          *s_therm
            else
               s_tot = s_tot + al*s_scat + em*(1.d0-al)*s_therm
            end if
         end do

         ! Return total zodi signal, and store final value in cache for later use
         s_zodi(i)              = s_tot
         cache(lookup_idx, det) = s_tot

         
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

      deallocate(cache)
      
    end subroutine get_zodi_emission


    subroutine get_zodi_emission_adaptive(tod, pix, scan, det, model, s_zodi, accuracy, samprate_min, comp)
      implicit none
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
      class(comm_tod),               intent(in)           :: tod
      integer(i4b),    dimension(:), intent(in)           :: pix
      integer(i4b),                  intent(in)           :: scan, det
      type(ZodiModel),               intent(in)           :: model
      real(sp),        dimension(:), intent(out)          :: s_zodi
      real(sp),                      intent(in), optional :: accuracy
      real(sp),                      intent(in), optional :: samprate_min
      integer(i4b),                  intent(in), optional :: comp

      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr, nsamp, dn, iter, maxiter, ncurr, n_new, offset, offset2
      logical(lgt) :: scattering, thermal, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, obs_time, lat, lon, s_tot, s_scat, s_therm, al, em, acc
      real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3), samprate
      type(spline_type) :: s_int
      integer(i4b), allocatable, dimension(:) :: ind
      real(dp),     allocatable, dimension(:) :: b_nu
      real(dp),     allocatable, dimension(:) :: s, s_pred
            
      acc        = 1d-3; if (present(accuracy))     acc      = accuracy  ! Minimum relative accuracy
      samprate   = 1.d0; if (present(samprate_min)) samprate = samprate_min
      dn         = max(1,int(tod%samprate/samprate)) ! Default sampling distance
      maxiter    = 10
      
      n_tod      = size(pix)
      dt_tod     = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_pos    = tod%scans(scan)%x0_obs
      earth_pos  = tod%scans(scan)%x0_earth
      R_obs      = norm2(obs_pos)
      obs_time   = tod%scans(scan)%t0(1)
      earth_lon  = atan(earth_pos(2), earth_pos(1))
      scattering = tod%central_freq >= model%nu_min_scatter
      thermal    = tod%central_freq <= model%nu_max_thermal

      allocate(s(n_tod), ind(n_tod), s_pred(n_tod))

      ! Initialize adaptive grid
      ncurr = 0
      do i = 1, n_tod, dn
         ncurr      = ncurr +1
         ind(ncurr) = i
         s(ncurr)   = evaluate_zodi(i) ! Get zodi at i
      end do
      if (ind(ncurr) < n_tod) then
         ! Add endpoint manually
         ncurr      = ncurr + 1
         ind(ncurr) = n_tod
         s(ncurr)   = evaluate_zodi(n_tod) ! Get zodi at n_tot
      end if
      ind(1)     = -ind(1)     ! Fix first point
      ind(ncurr) = -ind(ncurr) ! Fix last point

!!$      open(58,file="spline1.dat")
!!$      do i = 1, ncurr
!!$         write(58,*) i, ind(i), abs(ind(i)), s(i)
!!$      end do
!!$      close(58)
      
      ! Compute spline
      call spline(s_int, real(abs(ind(1:ncurr)),dp), s(1:ncurr))

      
      ! Refine grid until all segments have reached sufficient accuracy
      do iter  = 1, maxiter

         ! Count new points to be evaluated; skip those with both negative indices
         offset = 0
         do i = ncurr, 2, -1
            if (abs(ind(i))-abs(ind(i-1)) > 1 .and. (ind(i) >= 0 .or. ind(i-1) >= 0)) then
               offset = offset + 1
            end if
         end do
         offset2 = offset
         
         ! Make space for new points; fill with new values
         !write(*,*) "offset = ", offset
         s_pred(1:ncurr+offset) = 0.
         do i = ncurr, 2, -1
!            write(*,*) i, i+offset, offset, ind(i), s(i)
            ind(i+offset) = ind(i)
            s(i+offset)   = s(i)
            if (abs(ind(i))-abs(ind(i-1)) > 1 .and. (ind(i) >= 0 .or. ind(i-1) >= 0)) then
               ! Evaluate and insert new midpoint
               offset           = offset-1               
               j                = (abs(ind(i-1))+abs(ind(i)))/2
               ind(i+offset)    = j
               s(i+offset)      = evaluate_zodi(j)
               s_pred(i+offset) = splint(s_int, real(j,dp))
               !write(*,*) 'adding', j, s(i+offset), s_pred(i+offset), abs((s(i+offset)-s_pred(i+offset))/s(i+offset))
            end if
         end do
         ncurr = ncurr + offset2

         ! Check accuracy
         do i = 2, ncurr-1
            if (abs(ind(i)) == abs(ind(i-1))+1) ind(i-1:i) = -abs(ind(i-1:i))
            if (s_pred(i) /= 0.) then
               ! Recent point
               if (abs((s(i)-s_pred(i))/s(i)) < acc) ind(i-1:i+1) = -abs(ind(i-1:i+1))
            end if
         end do
         
!!$         open(58,file="spline2.dat")
!!$      do i = 1, ncurr
!!$         write(58,*) i, ind(i), abs(ind(i)), s(i), abs((s(i)-s_pred(i))/s(i))
!!$      end do
!!$      close(58)

         ! Compute spline
         call free_spline(s_int)
         call spline(s_int, real(abs(ind(1:ncurr)),dp), s(1:ncurr))

!!$      open(58,file="spline3.dat")
!!$      do i = 1, n_tod
!!$         write(58,*) i, splint(s_int, real(i,dp))
!!$      end do
!!$      close(58)

         
         if (all(ind < 0)) exit
      end do

      ! Interpolate to full scan
      do i = 1, n_tod
         s_zodi(i) = splint(s_int, real(i,dp))
      end do

      !write(*,fmt='(a,i8,i4,i12,i12)') '  Zodi spline = ', tod%scanid(scan), det, ncurr, n_tod
      call free_spline(s_int)
         
    contains

      function evaluate_zodi(ind) result(s)
        implicit none
        integer(i4b), intent(in) :: ind
        real(dp)                 :: s

        integer(i4b) :: j, k
        real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3), samprate
        real(dp) :: earth_lon, R_obs, R_min, R_max, lat, lon, s_tot, s_scat, s_therm, al, em

        unit_vector = tod%pixcache%ind2vec_ecl(:, tod%pixcache%pix2ind(pix(ind)))
        obs_time    = tod%scans(scan)%t0(1) + (ind-1)*dt_tod 
        do j = 1, 3
           earth_pos(j) = splint_simple(tod%x_earth_spline(j), obs_time)
           obs_pos(j)   = splint_simple(tod%x_obs_spline(j),   obs_time)
        end do
        R_obs      = norm2(obs_pos)
        earth_lon  = atan(earth_pos(2), earth_pos(1))

        s = 0.d0
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
              comp_LOS(k)%X_unit(l,:) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
              comp_LOS(k)%X_unit(l,:) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
              comp_LOS(k)%X(l,:)      = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
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

           ! Get dust grain temperature along the line of sight, and compute splined blackbody emission
           if (trim(model%phasefunc_type) == 'Wright' .and. k > 1) then
               !comp_LOS(k)%T = exp(5.5301d0) * comp_LOS(k)%R**(-0.5d0)
               comp_LOS(k)%T = exp(5.5301d0) / sqrt(comp_LOS(k)%R)
           else
               !comp_LOS(k)%T = model%T_0 * comp_LOS(k)%R**(-model%delta)
               comp_LOS(k)%T = model%T_0 * exp(-model%delta * log(comp_LOS(k)%R))
           end if
           call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

           
           ! Compute density
           call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
           
           if (scattering) then
              ! Compute phase function if scattering is included
              if (trim(model%phasefunc_type) == 'Wright') then
                 allocate(b_nu(size(tod%bp(0)%p%nu)))
                 call get_blackbody_emission(tod%bp(0)%p%nu, 5772.d0, b_nu) 
                 comp_LOS(k)%F_sol = tsum(tod%bp(0)%p%nu, tod%bp(0)%p%tau*b_nu)/comp_LOS(k)%R**2
                 deallocate(b_nu)
              else
                 comp_LOS(k)%F_sol = model%F_sun(tod%zodiband)/comp_LOS(k)%R**2
              end if
              call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
              call model%get_phase_function(comp_LOS(k)%Theta, tod%zodiband, comp_LOS(k)%Phi)
              s_scat = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5d0*(R_max-R_MIN)*1d20
           else
              s_scat = 0.d0
           end if

           if (thermal) then
              s_therm = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights) &
                   & * 0.5d0*(R_max-R_MIN)*1d20
           else
              s_therm = 0.d0
           end if

           ! Get albedo and emussivity
           al = zodi_model%comps(k)%c%albedo(tod%id)
           em = zodi_model%comps(k)%c%emissivity(tod%id)
           if (trim(zodi_model%phasefunc_type) == 'Wright') then
              s = s + al*s_scat + em          *s_therm
           else
              s = s + al*s_scat + em*(1.d0-al)*s_therm
           end if
        end do

      end function evaluate_zodi

    end subroutine get_zodi_emission_adaptive


   subroutine get_zodi_emission_mbb(tod, pix, scan, det, model, s_zodi, use_lowres_pointing, comp)
     implicit none
      ! Returns the predicted zodiacal emission for a scan (chunk of time-ordered data).
      ! Compute thermal term fitting a BB+MBB model
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
      class(comm_tod),               intent(inout)        :: tod
      integer(i4b),    dimension(:), intent(in)           :: pix
      integer(i4b),                  intent(in)           :: scan, det
      type(ZodiModel),               intent(in)           :: model
      real(sp),        dimension(:), intent(out)          :: s_zodi
      logical(lgt),                  intent(in), optional :: use_lowres_pointing
      integer(i4b),                  intent(in), optional :: comp

      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr
      logical(lgt) :: scattering, thermal, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, obs_time, lat, lon, s_tot, s_scat, s_therm, al, em
      real(dp) :: mbb_cutoff, mbb_ampl, mbb_b
      real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3)
      real(dp), allocatable, dimension(:)   :: b_nu

      integer(i4b) :: cache_hits
      real(sp), allocatable, dimension(:,:)     :: cache
      real(dp)                                  :: cache_time, init_cache_time
      
      use_lowres = .false.; if (present(use_lowres_pointing)) use_lowres = use_lowres_pointing
      n_tod      = size(pix)
      dt_tod     = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_pos    = tod%scans(scan)%x0_obs
      earth_pos  = tod%scans(scan)%x0_earth
      R_obs      = norm2(obs_pos)
      obs_time   = tod%scans(scan)%t0(1)
      earth_lon  = atan(earth_pos(2), earth_pos(1))
      scattering = tod%central_freq >= model%nu_min_scatter
      thermal    = tod%central_freq <= model%nu_max_thermal

      ! Initialize cache
      if (use_lowres) then
         !allocate(cache(tod%zodi_cache_nobs_lowres, tod%ndet))
      else
         allocate(cache(tod%pixcache%nobs,                   tod%ndet))
      end if
      cache_time = tod%scans(1)%t0(1)
      cache      = -1
      
      cache_hits = 0
!!$      open(58,file="zodi.dat",recl=1024)
      do i = 1, n_tod
         ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
         ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
         if (use_lowres) then
            obs_time = tod%scans(scan)%d(det)%downsamp_obs_time_full(i)
         else
            obs_time = obs_time + dt_tod ! assumes a time continuous TOD
         end if

         ! Reset cache
         if ((obs_time - cache_time) >= delta_t_reset) then
            do j = 1, 3
               earth_pos(j) = splint_simple(tod%x_earth_spline(j), obs_time)
               obs_pos(j)   = splint_simple(tod%x_obs_spline(j),   obs_time)
            end do
            R_obs      = norm2(obs_pos)
            earth_lon  = atan(earth_pos(2), earth_pos(1))
            cache      = -1
            cache_time =  obs_time
         end if

         ! Get lookup index for cache. If the pixel is already cached, used that value.
         if (use_lowres) then
            !lookup_idx  = tod%pix2ind_lowres(tod%udgrade_pix_zodi(pix(i))) 
            !unit_vector = tod%ind2vec_ecl_lowres(:, lookup_idx)
         else
            lookup_idx  = tod%pixcache%pix2ind(pix(i))
            if (lookup_idx == -2) write(*,*) 'oor', tod%scanid(scan), det, i, pix(i)
            unit_vector = tod%pixcache%ind2vec_ecl(:, lookup_idx)
         end if
         
         if (cache(lookup_idx, det) > 0.d0) then
            s_zodi(i)  = cache(lookup_idx, det)
            cache_hits = cache_hits + 1
            cycle
         end if

         ! If not present in cache, compute signal for each zodi component from scratch, and add signals together
         s_tot = 0.
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
               comp_LOS(k)%X_unit(l,:) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
               comp_LOS(k)%X_unit(l,:) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
               comp_LOS(k)%X(l,:)      = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
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

            ! Compute phase function if scattering is included
            if (scattering) then
               if (trim(model%phasefunc_type) == 'Wright') then
                  allocate(b_nu(size(tod%bp(0)%p%nu)))
                  call get_blackbody_emission(tod%bp(0)%p%nu, 5772.d0, b_nu) 
                  comp_LOS(k)%F_sol = tsum(tod%bp(0)%p%nu, tod%bp(0)%p%tau*b_nu)/comp_LOS(k)%R**2
                  deallocate(b_nu)
               else
                  comp_LOS(k)%F_sol = model%F_sun(tod%zodiband)/comp_LOS(k)%R**2
               end if
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call model%get_phase_function(comp_LOS(k)%Theta, tod%zodiband, comp_LOS(k)%Phi)
            end if

            ! Get dust grain temperature along the line of sight, and compute splined blackbody emission
            if (trim(model%phasefunc_type) == 'Wright' .and. k > 1) then
               !comp_LOS(k)%T = exp(5.5301d0) * comp_LOS(k)%R**(-0.5d0)
               comp_LOS(k)%T = exp(5.5301d0) / sqrt(comp_LOS(k)%R)
            else
               !comp_LOS(k)%T = model%T_0 * comp_LOS(k)%R**(-model%delta)
               comp_LOS(k)%T = model%T_0 * exp(-model%delta * log(comp_LOS(k)%R))
            end if
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

!!$            ! Compute density along line of sight
!!$            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            
            ! Compute predicted signal in MJy/sr; store computed signal in cache
            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            s_scat = 0.; s_therm = 0.
            if (scattering) s_scat  = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5d0*(R_max-R_MIN)*1d20
            if (thermal)    s_therm = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights)                  * 0.5d0*(R_max-R_MIN)*1d20

            ! BB+MBB thing only for cloud at this stage
            if (trim(model%comp_labels(k)) == 'cloud') then 
               al = zodi_model%comps(k)%c%albedo(tod%id)
               em = zodi_model%comps(k)%c%emissivity(tod%id)
               if (trim(zodi_model%phasefunc_type) == 'Wright') then
                  s_tot = s_tot + al*s_scat + em*s_therm
               else
                  mbb_cutoff = model%comps(k)%c%sed_cutoff
                  if (tod%central_freq < mbb_cutoff) then
                     mbb_ampl = model%comps(k)%c%sed_ampl
                     mbb_b = zodi_model%comps(k)%c%sed_b
                     s_tot = s_tot + al*s_scat + mbb_ampl*(tod%central_freq**mbb_b)*s_therm
                     !how to write this such that BB and MBB form a continuous function?
                  else 
                     s_tot = s_tot + al*s_scat + s_therm
                  end if
               end if 
            else 
               al = zodi_model%comps(k)%c%albedo(tod%id)
               em = zodi_model%comps(k)%c%emissivity(tod%id)
               if (trim(zodi_model%phasefunc_type) == 'Wright') then
                  s_tot = s_tot + al*s_scat + em          *s_therm
               else
                  s_tot = s_tot + al*s_scat + em*(1.d0-al)*s_therm
               end if
            end if
         end do

         ! Return total zodi signal, and store final value in cache for later use
         s_zodi(i)              = s_tot
         cache(lookup_idx, det) = s_tot

         
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

      deallocate(cache)
      
    end subroutine get_zodi_emission_mbb

    
   subroutine compute_zodi(zodi_model, tod, sd, det, comp)
      implicit none
      class(ZodiModel),     intent(in)             :: zodi_model
      class(comm_tod),      intent(in)             :: tod
      class(comm_scandata), intent(inout)          :: sd
      integer(i4b),         intent(in),   optional :: det
      integer(i4b),         intent(in),   optional :: comp
      
      integer(i4b) :: i, j, d, h, hp, ntod, nhorn, ncomp, ndet
      real(sp)     :: w
      real(dp)     :: t1, t2

      ntod  = sd%ntod
      nhorn = tod%nhorn
      ndet  = tod%ndet;           if (present(det)) ndet = 1
      ncomp = zodi_model%n_comps; if (present(comp)) ncomp = 1

      ! Compute non-stationary zodi TOD through line-of-sight integration for each horn, and add together
      do j = 1, ndet
         d = j; if (present(det)) d = det
         do h = 1, nhorn
            hp = h; if (nhorn == 1) hp = 0
            ! Fast, but approximate
            call get_zodi_emission_adaptive(tod, sd%pix(:,j,h), sd%scan, d, zodi_model, sd%s_zodi(:,j,hp), accuracy=1e-3, samprate_min=1.0)
            ! Evaluate each sample
            !call get_zodi_emission(tod, pix_dynamic(:,h), scan, det, zodi_model, s_zodi, comp=comp)
         end do
      end do
      
    end subroutine compute_zodi

    
   ! Functions for evaluating the zodiacal emission
   ! -----------------------------------------------------------------------------------
    ! Computes R_max (the length of the LOS such that it stops exactly at los_cutoff_radius).
   subroutine get_sphere_intersection(unit_vector, obs_pos, R_obs, R_cutoff, R_intersection)
     implicit none
     real(dp), intent(in), dimension(:) :: unit_vector, obs_pos
     real(dp), intent(in) :: R_obs, R_cutoff
     real(dp), intent(out) :: R_intersection
     
     real(dp) :: lon, lat, cos_lat, b, d, q
     
     if (R_obs > R_cutoff) then
        R_intersection = EPS
        return
     end if
     
     lon     = atan(unit_vector(2), unit_vector(1))
     lat     = asin(unit_vector(3))
     cos_lat = cos(lat)
     b       = 2.d0*(obs_pos(1)*cos_lat*cos(lon) + obs_pos(2)*cos_lat*sin(lon))
     d       = R_obs**2 - R_cutoff**2
     q       = -0.5d0*b*(1.d0+sqrt(b**2 - (4.d0*d))/abs(b))
     R_intersection = max(q, d/q)
   end subroutine get_sphere_intersection

   subroutine get_blackbody_emission(nus, T, b_nu)
     implicit none
     real(dp),               intent(in)  :: nus(:), T
     real(dp), dimension(:), intent(out) :: b_nu
     
     integer(i4b) :: i
     real(dp) :: x
     do i = 1, size(nus)
        x = h*nus(i)/(k_B*T)
        if (x < 0.001d0) then
           ! Use RJ approximation
           b_nu(i) = 2.d0*nus(i)**2*k_B*T/c**2
        else if (x > 50.d0) then
           ! Use Wien approximation
           b_nu(i) = 2.d0*h*nus(i)**3/c**2 * exp(-x)
        else
           ! Use exact expression
           b_nu(i) = 2.d0*h*nus(i)**3/c**2 / (exp(x) - 1.d0)
        end if
     end do
     !b_nu = b_nu * 1d20 ! Convert from W/(m^2*sr*Hz) to MJy/sr
   end subroutine get_blackbody_emission

   subroutine get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, scattering_angle)
     implicit none
     real(dp),               intent(in) :: X_helio_vec_LOS(:, :), X_vec_LOS(:, :), R_helio_LOS(:)
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
     where (cos_theta > 1.d0)
        scattering_angle = 0.d0
     elsewhere(cos_theta < -1.d0)
        scattering_angle = pi
     elsewhere
        scattering_angle = acos(-cos_theta)
     end where
     deallocate(cos_theta, R_LOS)
   end subroutine get_scattering_angle

   subroutine get_phase_function(self, Theta, band, Phi)
     implicit none
     class(ZodiModel), intent(in)           :: self
     real(dp),         intent(in)           :: Theta(:)
     integer(i4b),     intent(in)           :: band
     real(dp),         intent(out)          :: Phi(:)

     real(dp) :: C0, C1, C2, p20, p21, g1, g2, g3, w1, w2, w3, norm
     integer(i4b) :: n
     real(dp), allocatable :: den1(:), den2(:), den3(:), cosTheta(:)
     
     if (trim(self%phasefunc_type) == 'K98') then
        n    = self%numband
        C0   = self%par_phase(band+0*n)
        C1   = self%par_phase(band+1*n)
        C2   = self%par_phase(band+2*n)
        norm =  1.d0 / (2.d0*pi * (2.d0*C0 + pi*C1 + (exp(C2 * pi) + 1.d0)/(C2**2 + 1.d0)))
        Phi = norm * (C0 + (C1 * Theta) + exp(C2 * Theta))
     else if (trim(self%phasefunc_type) == 'Wright') then
        cosTheta = cos(Theta)
        p20  = self%par_phase(1)
        p21  = self%par_phase(2)
        Phi = exp(-p20*cosTheta + p21*cosTheta**2)
     else if (trim(self%phasefunc_type) == 'Hong') then
        cosTheta = cos(Theta)
        g1  = self%par_phase(1)
        g2  = self%par_phase(2)
        g3  = self%par_phase(3)
        w1  = 1.d0-sum(self%par_phase(4:5))
        w2  = self%par_phase(4)
        w3  = self%par_phase(5)

        den1 = (1.d0+g1**2-2d0*g1*cosTheta)
        den2 = (1.d0+g2**2-2d0*g2*cosTheta)
        den3 = (1.d0+g3**2-2d0*g3*cosTheta)

        Phi =       w1 * (1d0-g1**2)/(den1*sqrt(den1))
        Phi = Phi + w2 * (1d0-g2**2)/(den2*sqrt(den2))
        Phi = Phi + w3 * (1d0-g3**2)/(den3*sqrt(den3)) 
        Phi = Phi / (4.d0*pi)
     else
        write(*,*) 'Unsupported zodi phase function type:', trim(self%phasefunc_type)
        stop
     end if

   end subroutine get_phase_function

   subroutine initialize_earth_pos_spline(cpar)
      ! Returns the spline object which is used to evaluate the earth position

      type(comm_params), intent(in) :: cpar

      integer :: i, n_earthpos, unit
      real(dp), allocatable :: tabulated_earth_time(:), tabulated_earth_pos(:, :)
      unit = getlun()
      !open (unit, file=trim(trim(cpar%datadir)//'/'//trim(cpar%ephemerides_file)))
      open (unit, file=trim(cpar%ephemerides_file))
      read (unit, *) n_earthpos
      read (unit, *) ! skip header

      allocate (tabulated_earth_pos(3, n_earthpos))
      allocate (tabulated_earth_time(n_earthpos))
      do i = 1, n_earthpos
         read (unit, *) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
      end do
      close (unit)
      do i = 1, 3
         call spline_simple(earth_pos_spl_obj(i), tabulated_earth_time, tabulated_earth_pos(i, :), regular=.true.)
      end do
   end subroutine initialize_earth_pos_spline

   subroutine update_zodi_splines(tod, bandpass, det, model)
     implicit none
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
      class(comm_tod),    intent(inout) :: tod
      class(comm_bp_ptr), intent(in)    :: bandpass
      integer(i4b),       intent(in)    :: det
      type(ZodiModel),    intent(inout) :: model

      real(dp), allocatable :: b_nu(:)
      integer(i4b) :: i, j
      real(dp)     :: K, Inu0(1)

      allocate (b_nu(bandpass%p%n))
      do i = 1, size(B_nu_integrals)
         call get_blackbody_emission( bandpass%p%nu,    T_grid(i), b_nu) ! MJy/sr
         if (trim(model%bandpass_type) == 'DIRBE') then
            call get_blackbody_emission([bandpass%p%nu_c], T_grid(i), Inu0) ! Center frequency for color correction
            K     = tsum(bandpass%p%nu, bandpass%p%tau * b_nu/Inu0(1)) / tsum(bandpass%p%nu, bandpass%p%tau * bandpass%p%nu_c/bandpass%p%nu) ! Color correction
            B_nu_integrals(i) = K * Inu0(1)
         else if (trim(model%bandpass_type) == 'Wright') then
            B_nu_integrals(i) = tsum(bandpass%p%nu, bandpass%p%tau*b_nu)
         else if (trim(model%bandpass_type) == 'default') then
            b_nu = 1.d0/(2.d0*bandpass%p%nu**2*k_b/c**2 * 1d14) * b_nu ! uK_RJ
            B_nu_integrals(i) = tod%bp(det)%p%SED2F(b_nu)              ! Converts to data units
         else
            write(*,*) 'Update_zodi_splines -- unknown bandpass type = ', trim(model%bandpass_type)
         end if
      end do
      call spline_simple(tod%zodi_b_nu_spl_obj(det), T_grid, B_nu_integrals, regular=.true.)
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
      character(len=2048) :: chainfile, emissivity_path, albedo_path, band_path, comp_path, tod_path, group_name

      !allocate(tod%zodi_emissivity(tod%zodi_n_comps))
      !allocate(tod%zodi_albedo(tod%zodi_n_comps))
      lambda_min = 0.1
      lambda_max = 4.
      do i = 1, cpar%zs_ncomps
         if (trim(adjustl(cpar%zs_init_hdf(i))) == 'none') then
            model%comps(i)%c%emissivity(tod%id) = 1.
            lambda = (c / tod%nu_c(1)) * 1e6 ! in microns
            if ((lambda_min < lambda) .and. (lambda_max > lambda)) then
               model%comps(i)%c%albedo(tod%id) = 0.5
            else 
               model%comps(i)%c%albedo(tod%id) = 0.
            end if
            cycle
         end if
         
         if (cpar%myid == cpar%root) then
            unit = getlun()
            if (trim(cpar%zs_init_hdf(i)) == 'default') then
               call get_chainfile_and_samp(trim(cpar%init_chain_prefixes(1)), chainfile, initsamp)
            else
               call get_chainfile_and_samp(trim(cpar%zs_init_hdf(i)), chainfile, initsamp)
            end if
            inquire(file=trim(chainfile), exist=exist)
            if (.not. exist) call report_error('Zodi init chain does not exist = ' // trim(chainfile))
            l = len(trim(chainfile))
            if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
            call open_hdf_file(trim(chainfile), file, "r")
            
            call int2string(initsamp, itext)
            
            tod_path = trim(adjustl(itext//"/zodi/tod/"))
            band_path = trim(adjustl(tod_path))//trim(adjustl(tod%freq))

            if (.not. hdf_group_exists(file, band_path)) then
               print *, "Zodi init chain does not contain emissivities or albedos for band: " // trim(adjustl(tod%freq))
               stop
            end if 

            comp_path = trim(adjustl(band_path))//'/'//trim(adjustl(model%comp_labels(i)))//'/'
            if (hdf_group_exists(file, comp_path)) then
               call read_hdf(file, trim(adjustl(comp_path))//'/emissivity', model%comps(i)%c%emissivity(tod%id))
               call read_hdf(file, trim(adjustl(comp_path))//'/albedo', model%comps(i)%c%albedo(tod%id))
            else 
               model%comps(i)%c%emissivity(tod%id) = 1.
               model%comps(i)%c%albedo(tod%id)     = 0.
            end if 
         end if
      end do
      !call mpi_bcast(tod%zodi_emissivity, size(tod%zodi_emissivity), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      !call mpi_bcast(tod%zodi_albedo, size(tod%zodi_albedo), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
   end subroutine read_tod_zodi_params

   ! Get total zodi for a given TOD, including (optionally) both dynamic and static components, and accounting for horns

    subroutine zodi_model_to_ascii(cpar, model, filename, overwrite)
      ! Dumps the zodi model to an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(in) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      logical(lgt), intent(in), optional :: overwrite

      integer(i4b) :: io, i, j, running_idx
      logical(lgt) :: exists, overwrite_
      real(dp), allocatable :: params(:)
      integer(i4b), allocatable :: comp_switch_indices(:)
      character(len=128), allocatable :: labels(:)
      character(len=512) :: concatenated_string, val

      if (present(overwrite)) then
         overwrite_ = overwrite
      else
         overwrite_ = .false.
      end if

      if (cpar%myid_chain /= cpar%root) return
      inquire(file=trim(adjustl(filename)), exist=exists)
      if (exists .and. (.not. overwrite_)) then
         print *, "zodi asciifile: " // trim(adjustl(filename)) // " exists and overwrite = .false."
         stop
      end if

      open(newunit=io, file=trim(adjustl(filename)), action="write")
      allocate(params(model%npar_tot))
      call model%model_to_params(params, labels=labels)

      allocate(comp_switch_indices(0:model%n_comps))

      ! Newline after general parameters
      if (trim(model%phasefunc_type) == 'Hong') then
         comp_switch_indices(0) = model%n_general_params
      else
         comp_switch_indices(0) = 2
      end if

      ! Newline after each component
      running_idx = model%n_general_params
      do i = 1, model%n_comps
         running_idx = running_idx + size(model%comps(i)%labels)
         comp_switch_indices(i) = running_idx
      end do

      do i = 1, model%npar_tot
         if (trim(labels(i)) == 'skip' .or. trim(labels(i)) == 'SKIP') cycle
         if (any(comp_switch_indices == i)) then
            write(io, fmt='(a, T25, a, ES12.5, a)') trim(adjustl(labels(i))), "= ", params(i)
            write(io, *)
         else
            write(io, fmt='(a, T25, a, ES12.5)') trim(adjustl(labels(i))), "= ", params(i)
         end if
      end do

      write(io, fmt='(a)') ''
      do i = 1, numband
         if (trim(band_todtype(i)) == 'none') cycle
         !if (.not. data(i)%tod%subtract_zodi) cycle
         concatenated_string = ""
         do j = 1, model%n_comps
            write(val, fmt='(ES12.5)') model%comps(j)%c%emissivity(i) 
            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
         end do
         write(io, fmt='(a, T25, a, a)') trim(adjustl("EMISSIVITY_"//trim(adjustl(band_labels(i))))), "= ", trim(adjustl(concatenated_string(2:)))
      end do

      write(io, fmt='(a)') ''
      do i = 1, numband
         if (trim(band_todtype(i)) == 'none') cycle
         !if (.not. data(i)%tod%subtract_zodi) cycle
         concatenated_string = ""
         do j = 1, model%n_comps
            write(val, fmt='(ES12.5)') model%comps(j)%c%albedo(i)
            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
         end do
         write(io, fmt='(a, T25, a, a)') trim(adjustl("ALBEDO_"//trim(adjustl(band_labels(i))))), "= ", trim(adjustl(concatenated_string(2:)))
      end do

      ! Output static zodi amplitude
!      write(io, fmt='(a)') ''
!      do i = 1, numband
!         if (trim(band_todtype(i)) == 'none') cycle
!         !if (.not. data(i)%tod%subtract_zodi) cycle
!         write(val, fmt='(ES12.5)') model%amp_static(i)
!         write(io, fmt='(a, T25, a, a)') trim(adjustl("AMP_STATIC_"//trim(adjustl(band_labels(i))))), "= ", trim(adjustl(val))
!      end do
      
      close(io)
   end subroutine

   subroutine ascii_to_zodi_model(cpar, model, filename)
     implicit none
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
     
     allocate(params(model%npar_tot))
     !if (cpar%myid_chain == cpar%root) then
     inquire(file=trim(adjustl(filename)), exist=exists)
     if (.not. exists) then
        print *, "zodi asciifile: " // trim(adjustl(filename)) // " does not exist"
        stop
     end if
     
     call init_hash_tbl_sll(htbl, tbl_len=500)
     
     open(newunit=io, file=trim(adjustl(filename)), action="read")
     io_status = 0
     do while (io_status == 0)
        read(io, "(a)", iostat=io_status) line
        if (io_status == 0 .and. line /= "") then
           j = index(line, "=")
           if (j == 0) then
              print *, "Error: invalid line in ascii file: ", trim(adjustl(line))
              close(io)
              stop
           end if
           
           key = trim(adjustl(line(:j-1)))
           val = trim(adjustl(line(j+1:)))
           call tolower(key)
           call put_hash_tbl_sll(htbl, trim(adjustl(key)), trim(adjustl(val)))
        end if
     end do
     close(io)
     
     call model%model_to_params(params, labels=labels)
     params = 0.
     if (size(labels) /= size(params)) then
        write(*,*) labels
        write(*,*) params
        write(*,*) size(labels), size(params)
        stop "Error: size of labels and params do not match"
     end if
     do i = 1, size(labels)
        if (trim(labels(i)) == 'skip' .or. trim(labels(i)) == 'SKIP') cycle
        call get_parameter_hashtable(htbl, labels(i), par_dp=params(i))
     end do
     !end if

     !call mpi_bcast(params, size(params), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
     call model%params_to_model(params)
     
     do i = 1, numband
        if (trim(band_todtype(i)) == 'none') cycle
        !if (.not. data(i)%tod%subtract_zodi) cycle
        !if (cpar%myid == 0) then
        call get_parameter_hashtable(htbl, trim(adjustl("EMISSIVITY_"//trim(adjustl(band_labels(i))))), par_string=concatenated_string)
        call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
        if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model emissivity"
        do j = 1, n_comps
           read(toks(j), *) model%comps(j)%c%emissivity(i)
        end do
        
        call get_parameter_hashtable(htbl, trim(adjustl("ALBEDO_"//trim(adjustl(band_labels(i))))), par_string=concatenated_string)
        call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
        if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model albedo"
        do j = 1, n_comps
           read(toks(j), *) model%comps(j)%c%albedo(i)
        end do
        
        !            call get_parameter_hashtable(htbl, trim(adjustl("AMP_STATIC_"//trim(adjustl(band_labels(i))))), par_dp=zodi_model%amp_static(i))
        !end if
        !call mpi_bcast(data(i)%tod%zodi_emissivity, size(data(i)%tod%zodi_emissivity), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
        !call mpi_bcast(data(i)%tod%zodi_albedo, size(data(i)%tod%zodi_albedo), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
     end do
   end subroutine ascii_to_zodi_model
   
   subroutine print_zodi_model(theta, samp_group)
     implicit none
     real(dp),     allocatable, intent(in) :: theta(:)
     integer(i4b),              intent(in) :: samp_group
      
      integer(i4b) :: i, j, k, l, n, idx, n_cols, col, comp, n_params
      logical(lgt) :: newline
      
      n_cols = 5
      n_params = size(theta)

      ! General parameters
      k = 1
      if (any(zodi_model%theta_stat(1:zodi_model%n_general_params,samp_group)==0)) then
         write(*, "(a)", advance="no") 'General: '
         col = 1
         do i = 1, zodi_model%n_general_params
            if (zodi_model%theta_stat(i,samp_group)==0) then
               if (col > n_cols .or. i == zodi_model%n_general_params) then
                  write(*, "(a,a,g0.4,a)") trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
                  col = 1
               else
                  write(*, "(a,a,g0.4,a)", advance="no") trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
                  col = col+1
               end if
               k = k+1
            end if
         end do
      end if

      ! Component parameters
      do j = 1, zodi_model%n_comps
         idx = zodi_model%comps(j)%start_ind
         n   = zodi_model%comps(j)%npar + 2*numband
         if (all(zodi_model%theta_stat(idx:idx+n-1,samp_group)/=0)) cycle
         write(*, "(a,a)", advance="no") trim(adjustl(zodi_model%comp_labels(j))),': '
         col = 1
         do i = idx, idx+n-1
            newline = (i==idx+n-1) .or. col == n_cols-1
            if (zodi_model%theta_stat(i,samp_group)==0) then
               write(*, "(a,a,a,g0.4,a)", advance="no") "  ", trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
               k   = k+1
               col = col+1
            end if
            if (newline .and. col>1) then
               write(*,*)
               col = 1
            end if
         end do
      end do

      ! Monopoles
      col = 1
      write(*, "(a)", advance="no") 'Mono : '
      do j = 1, numband
         idx = zodi_model%npar_tot - numband + j
         if (zodi_model%theta_stat(idx,samp_group)==0) then
            write(*, "(a,a,a,g0.4,a)", advance="no") "  ", trim(adjustl(band_labels(j))), "=", theta(k), ",  "
            k   = k+1
            col = col+1
         end if
         if ((col == n_cols-1 .and. col>1) .or. j == numband) then
            write(*,*)
            col = 1
         end if
      end do
      
   end subroutine

   
   
 end module comm_zodi_mod
