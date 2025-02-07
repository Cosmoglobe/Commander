module comm_zodi_samp_mod
   use comm_zodi_mod
   use comm_chisq_mod
   use powell_mod
   implicit none

   private
   public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky, compute_downsamp_zodi, sample_zodi_group, sample_linear_zodi
   public minimize_zodi_with_powell, get_chisq_priors, precompute_lowres_zodi_lookups, create_zodi_glitch_mask
   public apply_zodi_glitch_mask, sample_static_zodi_map!, sample_static_zodi_amps
   
   real(dp), allocatable :: chisq_previous, step_size, prior_vec(:, :), prior_vec_powell(:, :), step_sizes_emissivity(:, :), step_sizes_albedo(:, :), step_sizes_ipd(:), step_sizes_n0(:), theta_0(:)
   real(dp), allocatable :: powell_emissivity(:, :), powell_albedo(:, :)
   integer(i4b) :: n_samp_bands, reference_emissivity_band
   real(dp) :: EPS = TINY(1.0_dp)
   real(dp), dimension(2) :: emissivity_prior, albedo_prior
   logical(lgt), allocatable :: powell_included_params(:), ref_band(:)
   character(len=128), allocatable, dimension(:) :: implemented_sampling_algorithms
   
contains
  subroutine initialize_zodi_samp_mod(cpar)
    implicit none
      ! Initialize the zodi sampling module.
      !
      ! Parameters
      ! ----------
      ! cpar: comm_params
      !    Parameter file variables.
      type(comm_params), intent(inout) :: cpar

      
    end subroutine initialize_zodi_samp_mod

   function get_boxwidth(samprate_lowres, samprate) result(box_width)
      ! Returns the boxcar width for downsampling the zodi tods
      real(dp), intent(in) :: samprate_lowres, samprate
      real(dp) :: box_width
      box_width = max(real(nint(samprate/samprate_lowres),dp),1.d0)
   end function


   subroutine sample_linear_zodi(cpar, handle, gibbs_iter, model, verbose)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model
      logical(lgt), intent(in), optional :: verbose

      !call compute_downsamp_zodi(cpar, model)
      !if (cpar%myid == cpar%root) print *, new_line('A'), "sampling and subtracting monopole"
      !call sample_and_subtract_monopole(cpar, handle)
      !if (cpar%myid == cpar%root) print *, new_line('A'), "sampling n0"
      call sample_n0_emissivity_and_albedo(cpar, handle, gibbs_iter, model)
   end subroutine

   subroutine sample_and_subtract_monopole(cpar, handle)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      
   end subroutine

   subroutine sample_n0_emissivity_and_albedo(cpar, handle, gibbs_iter, model)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model

   end subroutine

   subroutine sample_zodi_group(cpar, handle, gibbs_iter, model, verbose)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model
      logical(lgt), intent(in), optional :: verbose

   end subroutine

   function get_chisq_priors(params, samp_group) result(chisq_prior)
     implicit none 
     real(dp), intent(in) :: params(:)
     integer(i4b), intent(in) :: samp_group
      
      real(dp) :: chisq_prior
      integer(i4b) :: i, j
      logical(lgt) :: prior_is_violated

      chisq_prior = 0.d0
      j = 1
      do i = 1, zodi_model%npar_tot
         if (zodi_model%theta_stat(i,samp_group) /= 0) cycle

         if (params(j) < zodi_model%theta_prior(1,i) .or. &
              & params(j) > zodi_model%theta_prior(2,i)) then
            write(*,fmt='(a,a, 3e16.8)') 'Out of bounds -- ', trim(zodi_model%par_labels(i))//':', zodi_model%theta_prior(1,i), params(j), zodi_model%theta_prior(2,i)
            !write(*,*) 'Parameter out of bounds', i, zodi_model%theta_prior(1,i), params(j), zodi_model%theta_prior(2,i)
            chisq_prior = 1.d30
            return
         end if
         
         if (zodi_model%theta_prior(4,i) > 0.d0) then
            chisq_prior = chisq_prior + (params(j) - zodi_model%theta_prior(3,i))**2 / zodi_model%theta_prior(4,i)**2
         end if

         j = j+1
      end do
    end function get_chisq_priors

   subroutine downsamp_invariant_structs(cpar)
      ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
      type(comm_params), intent(in) :: cpar

   end subroutine
   

   subroutine precompute_lowres_zodi_lookups(cpar)
      ! Loop over each band with zodi and precompute lookup tables for lowres zodi caching
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, l, scan, n_lowres_obs, nobs_downsamp, nobs_lowres, pix_high, pix_low, ierr
      integer(i4b), allocatable :: pix2ind_highres(:), ind2pix_highres(:)
   end subroutine

   subroutine project_and_downsamp_sky(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ext(2), upper_bound, padding, ierr, ntod, nhorn, npix, ndet, nmaps, ndelta
   end subroutine

   subroutine compute_downsamp_zodi(cpar, model)
      type(comm_params), intent(in) :: cpar
      type(ZodiModel), intent(inout) :: model
      integer(i4b) :: i, j, scan, ierr, ndet

   end subroutine

   subroutine create_zodi_glitch_mask(cpar, handle)
     type(comm_params), intent(in) :: cpar
     type(planck_rng), intent(inout) :: handle
   end subroutine

   subroutine apply_zodi_glitch_mask(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ierr, non_glitch_size
      real(sp), allocatable :: res(:)
      real(sp), allocatable :: downsamp_scat_comp(:, :), downsamp_therm_comp(:, :)
      integer(i4b), allocatable :: downsamp_pix(:)

   end subroutine

   subroutine minimize_zodi_with_powell(cpar, iter, handle, samp_group)
      implicit none 
      type(comm_params), intent(in)      :: cpar
      integer(i4b),      intent(in)      :: iter
      type(planck_rng),  intent(inout)   :: handle
      integer(i4b),      intent(in)      :: samp_group


 end subroutine minimize_zodi_with_powell


   subroutine get_s_zodi_with_n0(s_therm, s_scat, s_zodi, emissivity, albedo, n_0_comp_ratio, comp)
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(inout) :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo, n_0_comp_ratio
      integer(i4b) :: comp
      integer(i4b) :: i

   end subroutine get_s_zodi_with_n0

   subroutine parse_samp_group_strings(samp_group_str, param_labels, param_indices)
      character(len=*), intent(in) :: samp_group_str
      character(len=*), intent(in) :: param_labels(:)
      integer(i4b), allocatable, intent(inout) :: param_indices(:)
      character(len=128) :: tokens(100), comp_param(2), label, param_label_tokens(10)
      character(len=128), allocatable :: tokens_trunc(:)
      integer(i4b) :: i, j, n_params
      logical(lgt) :: found


   end subroutine

 
   subroutine min_max_param_vec_with_priors(params, priors)
      real(dp), intent(inout) :: params(:)
      real(dp), intent(in) :: priors(:, :)
      integer(i4b) :: i

   end subroutine

   subroutine randomize_zodi_init(x, samp_group, cpar, handle, rms)
     implicit none
     real(dp),          dimension(:), intent(inout) :: x
     integer(i4b),                    intent(in)    :: samp_group
     type(comm_params),               intent(in)    :: cpar
     type(planck_rng),                intent(inout) :: handle
     real(dp),                        intent(in), optional :: rms
     
     integer(i4b) :: i, j, ierr
     real(dp)     :: eps


   end subroutine randomize_zodi_init


   subroutine model_to_params(zodi, x, samp_group)
     implicit none
     class(ZodiModel),               intent(in)           :: zodi
     real(dp),         dimension(:), intent(out)          :: x
     integer(i4b),                   intent(in), optional :: samp_group

     integer(i4b) :: i, j, idx
     real(dp), allocatable, dimension(:) :: z
     
     allocate(z(zodi%npar_tot))

     x = z

     
   end subroutine model_to_params

   subroutine params_to_model(zodi, x, samp_group)
     implicit none
     class(ZodiModel),               intent(inout)        :: zodi
     real(dp),         dimension(:), intent(in)           :: x
     integer(i4b),                   intent(in), optional :: samp_group

     integer(i4b) :: i, j, idx
     real(dp), allocatable, dimension(:) :: z, z_prev
     
     
   end subroutine params_to_model
   
   subroutine sample_static_zodi_map(cpar, handle)
     implicit none
      type(comm_params), intent(inout) :: cpar
      type(planck_rng), intent(inout) :: handle

      
    end subroutine sample_static_zodi_map


    
end module
