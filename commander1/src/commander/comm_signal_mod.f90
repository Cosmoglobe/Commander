module comm_signal_mod
  use comm_mp_mod
  use comm_direct_sol_mod
  use comm_cgd_mod
  use comm_Cl_sampling_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************


  character(len=128),                            private :: sampling_method
  character(len=128),                            private :: solution_method


contains

  subroutine initialize_signal_mod(paramfile)
    implicit none

    character(len=128),                   intent(in) :: paramfile

    call get_parameter(paramfile, 'SAMPLING_METHOD',    par_string=sampling_method)
    call get_parameter(paramfile, 'SOLUTION_METHOD',    par_string=solution_method)

    if (trim(sampling_method) == 'spline_sampler') then
!      call initialize_cl_spline_mod(paramfile)
   else
      call initialize_Cl_sampling_mod(paramfile)
   end if


  end subroutine initialize_signal_mod

  subroutine cleanup_signal_mod
    implicit none

    if (trim(sampling_method) == 'spline_sampler') then
!       call cleanup_cl_spline_mod
    else
      call cleanup_Cl_sampling_mod
   end if

  end subroutine cleanup_signal_mod

  

  ! Routine for sampling the signal component
  ! Inputs:  Residual maps, current power spectrum
  ! Outputs: Sampled a_lms and properly smoothed signal components for each band
  subroutine sample_signal_component(iter, precond_type, s_coeff, stat)
    implicit none

    integer(i4b), intent(in)    :: iter, precond_type
    type(genvec), intent(inout) :: s_coeff
    integer(i4b), intent(inout) :: stat

    integer(i4b) :: ncomp, nmaps
    type(genvec) :: rhs
    real(dp), allocatable, dimension(:,:) :: alm
    
    call allocate_genvec(rhs)
    call compute_rhs_product(.true., .true., .true., rhs)
!    write(*,*) 'sum = ', sum(abs(rhs%fg_amp))
!    stop
    if (trim(solution_method) == 'CG') then
       call draw_constrained_realization_by_CG(precond_type, rhs, s_coeff, stat)
    else if (trim(solution_method) == 'brute_force') then
       call solve_fundamental_system_direct(rhs, s_coeff, stat)
    else
       write(*,*) 'Error -- unknown solution method = ', trim(solution_method)
       stop
    end if
    call deallocate_genvec(rhs)
    
  end subroutine sample_signal_component



  ! Routine for sampling a power spectrum given a signal map
  subroutine sample_Cls(alms, rng_handle, cls)
    implicit none

    type(planck_rng),              intent(inout) :: rng_handle
    real(dp), dimension(1:,1:),    intent(in)    :: alms
    real(dp), dimension(0:,1:),    intent(inout) :: cls

    integer(i4b) :: i, j, l, m

    if (trim(sampling_method) == 'joint_sampler') then
       call sample_inv_wishart_spectrum(alms, rng_handle, cls)
    else if (trim(sampling_method) == 'spline_sampler') then
!       do l = 0, lmax
!          sigma_l(:,:,l) = sigma_l(:,:,l) / real(2*l+1,dp)
!       end do
!       call sample_splined_spectrum(sigma_l, rng_handle, cls)
    else
       write(*,*) 'comm_signal_mod: Unknown Cl sampler = ', trim(sampling_method)
       stop
    end if
    
  end subroutine sample_Cls

    
end module comm_signal_mod
