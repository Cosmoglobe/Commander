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
module comm_nonlin_mod
  use comm_param_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_gain_mod
  use comm_line_comp_mod
  use comm_diffuse_comp_mod
  use comm_signal_mod
  use comm_utils
  use InvSamp_mod
  use powell_mod
  use comm_output_mod
  implicit none

interface


  module subroutine sample_nonlin_params(cpar, iter, handle, handle_noise)
    !
    ! Routine that loops through all components and samples the spectral parameters that are defined with
    ! non-zero RMS values
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! handle_noise: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handles are updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle, handle_noise    

    integer(i4b) :: i, j, p
    real(dp)     :: t1, t2
    logical(lgt) :: samp_cg
    class(comm_comp),    pointer :: c    => null()
    character(len=512), allocatable, dimension(:) :: comp_labels

    
  end subroutine sample_nonlin_params


  module subroutine sample_specind_alm(cpar, iter, handle, comp_id, par_id)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, r, q, p, pl, np, nlm, l_, m_, idx, delta, burnin
    integer(i4b) :: nsamp, out_every, check_every, num_accepted, smooth_scale, id_native, ierr, ind, nalm_tot_reg
    integer(i4b) :: p_min, p_max, nalm_tot, pix, region
    real(dp)     :: t1, t2, ts, dalm, thresh, steplen
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior, alms_mean, alms_var, chisq_jeffreys, chisq_temp
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext
    character(len=3) :: tag
    character(len=15) :: regfmt
    character(len=9) :: ar_tag
    character(len=1000) :: outmessage
    character(len=512) :: filename

    logical :: accepted, exist, doexit, optimize, apply_prior
    class(comm_mapinfo), pointer :: info => null()
    class(comm_mapinfo), pointer :: info_theta => null()
    class(comm_map),     pointer :: theta => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_map),     pointer :: theta_smooth => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_comp),    pointer :: c    => null()
    type(map_ptr),     allocatable, dimension(:) :: df

    real(dp),          allocatable, dimension(:,:,:)  :: alms, regs, buffer3
    real(dp),          allocatable, dimension(:,:)    :: m
    real(dp),          allocatable, dimension(:)      :: buffer, buffer2, rgs, chisq, theta_pixreg_prop, theta_delta_prop
    integer(c_int),    allocatable, dimension(:)      :: maxit
    real(dp)     :: theta_max, theta_min
    logical      :: outside_limit


  end subroutine sample_specind_alm

  module subroutine sample_specind_local(cpar, iter, handle, comp_id, par_id)
    !
    ! Routine that sets up the sampling using the local sampling routine  for the spectral parameter given by
    ! par_id for the component given by the comp_id parameter. 
    ! Then it calls on the specific sampling routine and finally updates the components spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx, p_ind, p_min, p_max
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind, ind_pol
    real(dp)     :: t1, t2, ts, dalm, fwhm_prior, temp_theta
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext
    logical :: accepted, exist, doexit, skip
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_map),     pointer :: temp_res  => null()
    class(comm_map),     pointer :: temp_map  => null()
    class(comm_map),     pointer :: temp_map2  => null()
    class(comm_comp),    pointer :: c    => null()
    class(comm_map),     pointer :: theta_single_pol => null() ! Spectral parameter of one poltype index (too be smoothed)
    real(dp),          allocatable, dimension(:,:,:)   :: alms
    real(dp),          allocatable, dimension(:,:) :: m
    real(dp),          allocatable, dimension(:) :: buffer, rgs, chisq

    integer(c_int),    allocatable, dimension(:,:) :: lm
    integer(i4b), dimension(MPI_STATUS_SIZE) :: mpistat

   
  end subroutine sample_specind_local


  module subroutine sample_specind_multi(cpar, iter, handle, comp_labels)
    !
    ! Routine that sets up the sampling using the local sampling routine  for the spectral parameter given by
    ! par_id for the component given by the comp_id parameter. 
    ! Then it calls on the specific sampling routine and finally updates the components spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    character(len=512), dimension(:), intent(in)    :: comp_labels

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx, p_ind, p_min, p_max, id, ncomp, npar, pix_out
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind, ind_pol, pol, pix
    real(dp)     :: t1, t2, ts, dalm, fwhm_prior, temp_theta
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    real(dp)     :: s, lnL_old, lnL_new, eps
    real(dp), allocatable, dimension(:) :: theta, theta_old
    character(len=2) :: itext
    logical :: accepted, exist, doexit, skip
    class(comm_mapinfo), pointer :: info => null()
    class(comm_mapinfo), pointer :: info_lowres => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_map),     pointer :: raw  => null()
    class(comm_map),     pointer :: raw2  => null()
    class(comm_map),     pointer :: raw_lowres  => null()
    class(comm_map),     pointer :: temp_res  => null()
    class(comm_map),     pointer :: temp_map  => null()
    class(comm_map),     pointer :: temp_map2  => null()
    class(comm_map),     pointer :: chisq_lowres  => null()
    class(comm_map),     pointer :: chisq_old  => null()
    class(comm_comp),    pointer :: c_in    => null()
    class(comm_diffuse_comp),    pointer :: c    => null()
    type(diff_ptr),    allocatable, dimension(:) :: comps
    real(dp),          allocatable, dimension(:,:,:)   :: alms
    real(dp),          allocatable, dimension(:,:) :: m
    real(dp),          allocatable, dimension(:) :: buffer, rgs, chisq

    integer(c_int),    allocatable, dimension(:,:) :: lm
    integer(i4b), dimension(MPI_STATUS_SIZE) :: mpistat



  end subroutine sample_specind_multi
  

  module subroutine gather_alms(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(in)    :: alm
    integer(c_int), dimension(1:,0:), intent(in) :: lm
    real(dp), dimension(0:,0:,1:), intent(inout) :: alms
    integer(i4b),                intent(in)    :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind


  end subroutine gather_alms

  module subroutine distribute_alms(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(inout)    :: alm
    integer(c_int), dimension(1:,0:), intent(in)   :: lm
    real(dp), dimension(0:,0:,1:),  intent(in)       :: alms
    integer(i4b),                intent(in)       :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind
    

  end subroutine distribute_alms

  module subroutine compute_corrlen(x, fix, n, maxit, corrlen)
    implicit none

    real(dp), dimension(:,:),    intent(in)    :: x
    logical(lgt), dimension(:),  intent(in)      :: fix        
    integer(i4b),                  intent(in)    :: n
    integer(i4b),                  intent(in)    :: maxit
    integer(i4b),                  intent(out)   :: corrlen

    real(dp),          allocatable, dimension(:) :: C_
    integer(c_int),    allocatable, dimension(:) :: N_      
    real(dp)     :: x_mean, x_var
    integer(i4b) :: pl, p, q, k, corrlen_init, delta

  end subroutine compute_corrlen


  ! Sample spectral parameters using straight Powell or inversion sampler
  module subroutine sampleDiffuseSpecInd_simple(cpar, handle, comp_id, par_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: par_id, comp_id, iter

    integer(i4b) :: pix, ierr, pol, status
    real(dp)     :: s, x(3), lnL_old, lnL_new, eps, theta_old
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()

    
  end subroutine sampleDiffuseSpecInd_simple


  ! Sample spectral parameters using straight Powell or inversion sampler
  module subroutine sampleDiffuseSpecInd_powell(cpar, handle, comp_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, iter

    integer(i4b) :: i, pix, ierr, pol, status
    real(dp)     :: s, x(3), lnL_old, lnL_new, eps
    real(dp), allocatable, dimension(:) :: theta, theta_old
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()

    
  end subroutine sampleDiffuseSpecInd_powell


  !Here comes all subroutines for sampling diffuse components locally
  ! Sample spectral parameters
  module subroutine sampleDiffuseSpecInd_nonlin(cpar, handle, comp_id, par_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: par_id, comp_id, iter

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_min, i_max, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer_lnL, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, accept_rate
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample
    class(comm_mapinfo),            pointer :: info => null()


  end subroutine sampleDiffuseSpecInd_nonlin



  module subroutine sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
    !
    ! Routine that samples diffuse type component spectral parameters in pixel regions
    ! 
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       A parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       Integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       Integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! buffer_lnL: double precision array
    !       An array copy of the current spectral parameter map, truncated by the absolute parameter limits.
    !       Array with dimension (0:npix-1,nmaps), where npix is the number of pixels given by the components 
    !       resolution parameter (Nside, see HEALPix), and nmaps is the number of polarizations (1 if only Temperature;
    !       3 if polarization is included, TQU)
    !
    ! p: integer
    !       Index counter for the polarization type that is to be sampled. Sets what map polarization(s) to be sampled.
    !
    ! Returns:
    !       No explicit parameter is returned, except for the sampled spectral parameter through the 
    !       'buffer_lnL' parameter.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    real(dp),               dimension(0:,:), intent(inout)        :: buffer_lnL
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, par_id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, m, n, q, pr, max_pr, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_mono,l_mono,m_mono, N_theta_MC, l_count
    integer(i4b) :: i_min, i_max, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps, nsamp
    real(dp)     :: a, b, a_tot, b_tot, s, t0, t1, t2, t3, t4, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: buff_init, theta_min, theta_max, running_accept, running_dlnL
    real(dp)     :: running_correlation, correlation_limit
    logical(lgt) :: ok, outside, sample_mono
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_comp),         pointer :: c2 => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    class(comm_md_comp),      pointer :: c_mono => null()
    !Following is for the local sampler
    real(dp)     :: mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, lnL_prior, lnL_init
    real(dp)     :: accept_rate, accept_scale, lnL_sum, proplen, chisq_jeffreys, avg_dlnL, lnL_total_init
    real(dp)     :: old_theta, new_theta, prior_rms_scaling
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, pix_count, buff1_i(1), buff2_i(1), burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    integer(i4b) :: npixreg, smooth_scale, arr_ind, np_lr, np_fr, myid_pix, unit
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in, sampled_nprop, sampled_proplen, first_nprop, sample_accepted
    character(len=512) :: filename, postfix, fmt_pix, npixreg_txt, monocorr_type
    character(len=6) :: itext
    character(len=4) :: ctext
    character(len=2) :: pind_txt
    character(len=512), dimension(1000) :: tokens
    real(dp),      allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),      allocatable, dimension(:) :: old_thetas, new_thetas, init_thetas, sum_theta
    real(dp),      allocatable, dimension(:) :: theta_corr_arr
    real(dp),      allocatable, dimension(:) :: old_theta_smooth, new_theta_smooth, dlnL_arr
    real(dp),  allocatable, dimension(:,:,:) :: theta_MC_arr
    integer(i4b),  allocatable, dimension(:) :: band_i, pol_j, accept_arr
    class(comm_mapinfo),             pointer :: info_fr => null() !full resolution
    class(comm_mapinfo),             pointer :: info_fr_single => null() !full resolution, nmaps=1
    class(comm_mapinfo),             pointer :: info_lr => null() !low resolution
    class(comm_mapinfo),             pointer :: info_lr_single => null() !low resolution, nmaps=1
    class(comm_mapinfo),             pointer :: info_data => null() !data band map info
    class(comm_map),                 pointer :: theta_single_fr => null() ! Spectral parameter of polt_id index
    class(comm_map),                 pointer :: theta_single_lr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_lr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_lr_hole => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_fr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: mask_lr => null() ! lowres mask
    class(comm_map),                 pointer :: mask_mono => null() ! lowres mask
    class(comm_map),                 pointer :: res_map => null() ! lowres  residual map
    class(comm_map),                 pointer :: ones_map => null() ! Constant sky
    class(comm_map),                 pointer :: Ninv_map => null() ! Inverse covariance diagonal
    class(comm_map),                 pointer :: temp_res => null() 
    class(comm_map),                 pointer :: temp_map => null() 
    class(comm_map),                 pointer :: temp_noise => null() 
    class(comm_map),                 pointer :: temp_mixing => null() 
    class(comm_map),                 pointer :: temp_chisq => null() 
    type(map_ptr), allocatable, dimension(:) :: df
    real(dp),      allocatable, dimension(:) :: monopole_val, monopole_rms, monopole_mu, monopole_mixing
    real(dp),      allocatable, dimension(:) :: old_mono, new_mono
    real(dp), allocatable, dimension(:,:,:,:) :: multipoles_trace !trace for proposed and accepted multipoles. (2,Nsamp,numband,0:3), the final dimension can be adjusted if only monopoles are estimated
    logical(lgt),  allocatable, dimension(:) :: monopole_active
    real(dp),    allocatable, dimension(:,:) :: reduced_data
    real(dp),    allocatable, dimension(:,:) :: harmonics, harmonics2
    real(dp),                 dimension(0:3) :: multipoles, md_b   
    real(dp),             dimension(0:3,0:3) :: md_A
    real(dp),                 dimension(3)   :: vector
    integer(i4b) :: i_md, j_md, k_md, max_prop


  end subroutine sampleDiffuseSpecIndPixReg_nonlin

  module function comp_lnL_marginal_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_marginal_diagonal

    integer(i4b) :: i, j, mat_len
    real(dp)     :: MNd,MNM,invMNM
    real(dp), dimension(:), allocatable :: MN
    !  Function to evaluate the log-likelihood for a pixel across all bands in one polarization,
    !  returning the log-likelihood value, equivalent to the highest likelihood chisq for the given theta/mixing matrix.
    !
    !  It does so by substituting the equation for the highest likelihood amplitude into the chisq equation and
    !  making a few assumptions on the behaviour of the theta values.
    !
    !  This evaluation is equivalen to the maximum chisq likelihood evaluation, 
    !  with the exception of a constant difference term.
    !
    !  This version of the evaluation is a little faster than the maximum chisq evaluation, but it does not allow
    !  for combined monopole sampling, as one can show that it will be biased compared to the true value of the
    !  sampled parameter.
    !
    !  Note: This evaluation assumes no correlation between pixels and are only evaluating on a pixel-by-pixel basis.
    !        Also, The length of each input array must be equal, and be larger than or equal to 'arr_len' if the latter
    !        is present.
    !
    !  Arguments:
    !  ------------------
    !  mixing: double precision array, unknown length
    !     An arraycontaining the pixel specific mixing matrix values for the different bands in the evaluation.
    !  invN_arr: double precision array, unknown length
    !     An array with the pixel specific inverse noise variance values for the different bands in the evaluation.
    !  data: double precision array, unknown length
    !     An array with the pixel specific (reduced) data values for the different bands in the evaluation.
    !  use_det: logical
    !     A logical flag specifying to add the logarithm of the determinant of the inverse MNM^T matrix to the
    !     log-likelihood. This is the case for the marginal likelihood evaluation.
    !  arr_len: integer(i4b), optional
    !     If present, one will only evaluate the first arr_len elements in the input arrays.
    !
    !  Returns:
    !  -----------------
    !  comp_lnL_marginal_diagonal: double precision real number
    !     The evaluated log-likelihood value for the pixel, the ridge/marginal likelihood.
    !

  end function comp_lnL_marginal_diagonal

  module function comp_lnL_max_chisq_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_max_chisq_diagonal

    integer(i4b) :: i, j, mat_len
    real(dp)     :: MNd,MNM,invMNM, amp, chisq
    real(dp), dimension(:), allocatable :: MN
    !  
    !  Evaluates "ridge likelihood" from BP13, arXiv:2201.08188. Take equation
    !  23 from this, but replace the integral with the the maximum likelihood
    !  solution for a given beta.
    !
    !  Implicitly assumes that we are evaluating for a single component.
    !
    !  Function to evaluate the log-likelihood for a pixel across all bands in one polarization,
    !  returning the highest likelihood chisq value of the log-likelihood for the given theta/mixing matrix.
    !
    !  It does so by solving the equation for the highest likelihood amplitude and using this amplitude
    !  to evaluate the chisq.
    !
    !  This evaluation is equivalen to the ridge/marginal evaluation, 
    !  with the exception of a constant difference term.
    !
    !  The benefit of this evaluation is that it allows for combined monopole sampling, 
    !  where one can show that the ridge/marginal evaluation does not return (or are able to find) 
    !  the true value of the sampled parameter, but will be biased. 
    !
    !  Note: This evaluation assumes no correlation between pixels and are only evaluating on a pixel-by-pixel basis.
    !        Also, The length of each input array must be equal, and be larger than or equal to 'arr_len' if the latter
    !        is present.
    !
    !  Arguments:
    !  ------------------
    !  mixing: double precision array, unknown length
    !     An array containing the pixel specific mixing matrix values for the different bands in the evaluation.
    !  invN_arr: double precision array, unknown length
    !     An array with the pixel specific inverse noise variance values for the different bands in the evaluation.
    !  data: double precision array, unknown length
    !     An array with the pixel specific (reduced) data values for the different bands in the evaluation.
    !  use_det: logical
    !     A logicalflag specifying to add the logarithm of the determinant of the inverse MNM^T matrix to the
    !     log-likelihood. This is the case for the marginal likelihood equivalent.
    !  arr_len: integer(i4b), optional
    !     If present, one will only evaluate the first arr_len elements in the input arrays.
    !
    !  Returns:
    !  -----------------
    !  comp_lnL_max_chisq_diagonal: double precision real number
    !     The evaluated log-likelihood value for the pixel, assuming the maximum likelihood chisq given the mixing matrix.
    !
    !  

  end function comp_lnL_max_chisq_diagonal


end interface

end module comm_nonlin_mod
