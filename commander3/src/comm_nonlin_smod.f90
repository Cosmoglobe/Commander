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
submodule (comm_nonlin_mod) comm_nonlin_smod

contains

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



  end subroutine sample_specind_multi
  

  module subroutine gather_alms(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(in)    :: alm
    integer(c_int), dimension(1:,0:), intent(in) :: lm
    real(dp), dimension(0:,0:,1:), intent(inout) :: alms
    integer(i4b),                intent(in)    :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind

    do k = 0, nalm-1
       ! Gather all alms
       l = lm(1,k)
       m = lm(2,k)
       ind = l**2 + l + m
       alms(i,ind,pl_tar) = alm(k,pl)
    end do

  end subroutine gather_alms

  module subroutine distribute_alms_nonlin(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(inout)    :: alm
    integer(c_int), dimension(1:,0:), intent(in)   :: lm
    real(dp), dimension(0:,0:,1:),  intent(in)       :: alms
    integer(i4b),                intent(in)       :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind
    
    do k = 0, nalm-1
       ! Distribute alms
       l = lm(1,k)
       m = lm(2,k)
       ind = l**2 + l + m
       alm(k,pl_tar) = alms(i,ind,pl)
    end do

  end subroutine distribute_alms_nonlin

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

    ! Calculate Correlation length
    delta = 100
    corrlen_init = 1
    corrlen = corrlen_init
    allocate(C_(delta))
    allocate(N_(delta))
    
    !open(58, file='correlation_function.dat', recl=10000)
          
    ! Calculate correlation function per parameter
    do p = 1, n
       if (fix(p)) cycle ! Skip fixed regions

       x_mean = mean(x(1:maxit,p))
       x_var = variance(x(1:maxit,p))
       
       N_ = 0
       C_ = 0.d0
       do q = 1, maxit
          do k = 1, delta
             if (q+k > maxit) cycle
             C_(k) = C_(k) + (x(q,p)-x_mean)*(x(q+k,p)-x_mean)
             N_(k) = N_(k) + 1 ! Less samples every q
          end do
       end do
       
       where (N_>0) C_ = C_/N_
       if ( x_var > 0 ) C_ = C_/x_var

      ! write(58,*) p, C_ ! Write to file

       ! Find correlation length
       do k = corrlen_init, delta
          if (C_(k) > 0.1) then
             if (k > corrlen) corrlen = k
          else
             exit
          end if
       end do
    end do
    deallocate(C_, N_)
    !close(58)
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
    real(dp)     :: x(3), lnL_old, lnL_new, eps, theta_old
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()

    !eps = 1d-9

    !c           => compList     
    !do while (comp_id /= c%id)
    !   c => c%nextComp()
    !end do
    !select type (c)
    !class is (comm_diffuse_comp)
    !   c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    !end select


!   ! do pol = 1, c_lnL%theta_smooth(par_id)%p%info%nmaps
    !pol = 1
    !do pix = 0, c_lnL%theta_smooth(par_id)%p%info%np-1
    !   if (c_lnL%theta_smooth(par_id)%p%map(pix,pol) <= c_lnL%p_uni(1,par_id)+eps) then
    !      x(1) = c_lnL%p_uni(1,par_id)+eps
    !      x(3) = min(x(1) + c_lnL%p_gauss(2,par_id), c_lnL%p_uni(2,par_id)-eps)
    !      x(2) = 0.5d0*(x(1)+x(3))
    !   else if (c_lnL%theta_smooth(par_id)%p%map(pix,pol) >= c_lnL%p_uni(2,par_id)-eps) then
    !      x(3) = c_lnL%p_uni(2,par_id)-eps
    !      x(1) = max(x(3) - c_lnL%p_gauss(2,par_id), c_lnL%p_uni(1,par_id)+eps)
    !      x(2) = 0.5d0*(x(1)+x(3))
    !   else
    !      x(2) = c_lnL%theta_smooth(par_id)%p%map(pix,pol)
    !      x(1) = max(x(2) - c_lnL%p_gauss(2,par_id), 0.5d0*(x(2) + c_lnL%p_uni(1,par_id)))
    !      x(3) = min(x(2) + c_lnL%p_gauss(2,par_id), 0.5d0*(x(2) + c_lnL%p_uni(2,par_id)))
    !      if (x(2) == x(1) .or. x(2) == x(3)) x(2) = 0.5d0*(x(1)+x(3))
    !   end if
    !   theta_old = x(2)

    !   if (mod(pix,1000) == 0) lnL_old = lnL_simple(c_lnL%theta_smooth(par_id)%p%map(pix,pol))
    !   c_lnL%theta_smooth(par_id)%p%map(pix,pol) = &
    !        & sample_InvSamp(handle, x, lnL_simple, prior=c_lnL%p_uni(:,par_id), &
    !        & optimize=(trim(cpar%operation)=='optimize'), status=status)
    !   if (mod(pix,1000) == 0) then
    !      lnL_new = lnL_simple(c_lnL%theta_smooth(par_id)%p%map(pix,pol))
    !      write(*,fmt='(i8,2f8.3,2e16.8)') c_lnL%theta_smooth(par_id)%p%info%pix(pix+1), theta_old, c_lnL%theta_smooth(par_id)%p%map(pix,pol), lnL_old, lnL_new
    !   end if
    !end do
    !!end do
!!$c!all mpi_finalize(ierr)
!!$ !   stop
    !call c_lnL%theta_smooth(par_id)%p%udgrade(c_lnL%theta(par_id)%p)

  contains

    function lnL_simple(x)
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_simple
      
      integer(i4b) :: i, j, n
      real(dp)     :: s, res, theta(6)
      real(dp), allocatable, dimension(:,:) :: f_precomp

      lnL_simple = 0d0
      

    end function lnL_simple
    
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
    real(dp)     :: lnL_old, lnL_new, eps
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

    integer(i4b) :: i, j, k, n, p, q, pix, ierr, ind(1), n_ok, id
    integer(i4b) :: npar, np, nmaps
    real(dp)     :: t1, t2, delta_lnL_threshold
    real(dp)     :: mu
    real(dp)     :: theta_min, theta_max
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    real(dp),     allocatable, dimension(:,:) :: buffer_lnL
    !Following is for the local sampler
    integer(i4b) :: p_min, p_max
    class(comm_mapinfo),            pointer :: info => null()

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    first_call          = .false.


    id = par_id
    c           => compList     
    do while (comp_id /= c%id)
       c => c%nextComp()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
       info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
            & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)


       npar      = c_lnL%npar
       np        = c_lnL%x_smooth%info%np
       nmaps     = c_lnL%x_smooth%info%nmaps

       theta_min = c_lnL%p_uni(1,id)
       theta_max = c_lnL%p_uni(2,id)


       allocate(buffer_lnL(0:c_lnL%theta(id)%p%info%np-1,c_lnL%theta(id)%p%info%nmaps))
       buffer_lnL=max(min(c_lnL%theta(id)%p%map,theta_max),theta_min) 
       
       
       do p = 1,c_lnL%poltype(id)
          if (c_lnL%lmax_ind_pol(p,id) >= 0) cycle !this set of polarizations are not to be local sampled (is checked before this point)
          if (c_lnL%poltype(id) > 1 .and. cpar%only_pol .and. p == 1) cycle !only polarization (poltype > 1)
          if (p > c_lnL%nmaps) cycle ! poltype > number of maps
          ! Return if all prior RMS's are zero

          call wall_time(t1)
          if (c_lnL%pol_pixreg_type(p,id) /= 0) then
             if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) '| Sampling poltype index', p, &
                  & 'of ', c_lnL%poltype(id) !Needed?
             call sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
             call wall_time(t2)
             if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) '| poltype:',c_lnL%poltype(id),' pol:', &
                  & p,'CPU time specind = ', real(t2-t1,sp)
          else
             write(*,*) '| Undefined spectral index sample region'
             write(*,*) '| Component: ',trim(c_lnL%label),', ind: ',trim(c_lnL%indlabel(id))
             stop
          end if



       end do

       !after sampling is done we assign the spectral index its new value(s)
       c_lnL%theta(id)%p%map = buffer_lnL !unsmoothed map (not post_proc smoothed)

       ! Ask for CG preconditioner update
       if (c_lnL%cg_unique_sampgroup > 0) recompute_diffuse_precond = .true.

       ! deallocate

       deallocate(buffer_lnL)

    end select

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


  end subroutine sampleDiffuseSpecIndPixReg_nonlin

  module function comp_lnL_marginal_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_marginal_diagonal

    integer(i4b) :: i, j
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

    if (present(arr_len)) then
       allocate(MN(arr_len))
       MN=mixing(1:arr_len)*invN_arr(1:arr_len)
       MNd=sum(MN*data(1:arr_len))
       MNM=sum(MN*mixing(1:arr_len))
    else
       allocate(MN(size(mixing)))
       MN=mixing*invN_arr
       MNd=sum(MN*data)
       MNM=sum(MN*mixing)
    end if

    comp_lnL_marginal_diagonal = 0.d0

    if (MNM /= 0.d0) then 
       invMNM=1.d0/MNM !invert 1x1 matrix
    else
       comp_lnL_marginal_diagonal=-1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    comp_lnL_marginal_diagonal = 0.5d0*MNd*invMNM*MNd

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_marginal_diagonal = comp_lnL_marginal_diagonal - 0.5d0*log(invMNM) 

    deallocate(MN)
  end function comp_lnL_marginal_diagonal

  module function comp_lnL_max_chisq_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_max_chisq_diagonal

    integer(i4b) :: i, j
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

    if (present(arr_len)) then
       allocate(MN(arr_len))
       MN=mixing(1:arr_len)*invN_arr(1:arr_len)
       MNd=sum(MN*data(1:arr_len))
       MNM=sum(MN*mixing(1:arr_len))
    else
       allocate(MN(size(mixing)))
       MN=mixing*invN_arr
       MNd=sum(MN*data)
       MNM=sum(MN*mixing)
    end if

    comp_lnL_max_chisq_diagonal = 0.d0

    if (MNM /= 0.d0) then 
       invMNM=1.d0/MNM !invert 1x1 matrix
       amp = MNd*invMNM
    else
       comp_lnL_max_chisq_diagonal=-1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    if (present(arr_len)) then
       chisq=sum((data(1:arr_len)-mixing(1:arr_len)*amp)**2 * invN_arr(1:arr_len))
    else
       chisq=sum((data-mixing*amp)**2 * invN_arr)
    end if

    comp_lnL_max_chisq_diagonal = -0.5d0*chisq

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_max_chisq_diagonal = comp_lnL_max_chisq_diagonal - 0.5d0*log(invMNM) 

    deallocate(MN)
  end function comp_lnL_max_chisq_diagonal


end submodule comm_nonlin_smod
