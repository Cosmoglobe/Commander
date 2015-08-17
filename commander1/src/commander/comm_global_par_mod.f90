module comm_global_par_mod
  use comm_mp_mod
  use sort_utils
  implicit none

  ! Internal data for global RMS sampling
  integer(i4b),                                private :: comp_rms, par_rms
  real(dp),     allocatable, dimension(:,:,:), private :: index_map_rms
  type(genvec),                                private :: s_rms
  real(dp),     allocatable, dimension(:,:,:), private :: x_init

  character(len=512), private :: chaindir

contains

  subroutine initialize_global_par_mod(paramfile)
    implicit none

    character(len=*), intent(in) :: paramfile

    call get_parameter(paramfile, 'CHAIN_DIRECTORY', par_string=chaindir)

  end subroutine initialize_global_par_mod

  subroutine sample_global_fg_par(handle, s, residuals, inv_N, fg_param_map, stat)
    implicit none

    type(planck_rng),                 intent(inout) :: handle
    type(genvec),                     intent(inout) :: s
    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N
    real(dp), dimension(0:,1:,1:),    intent(inout) :: fg_param_map
    integer(i4b),                     intent(inout) :: stat

    ! Sample parameter RMS's
    call sample_par_rms(handle, s, fg_param_map, stat)

    ! Optimize priors
    call optimize_priors(handle, s, residuals, inv_N, fg_param_map, stat)
    
  end subroutine sample_global_fg_par


  subroutine optimize_priors(handle, s, residuals, inv_N, fg_param_map, stat)
    implicit none

    type(planck_rng),                 intent(inout) :: handle
    type(genvec),                     intent(inout) :: s
    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N
    real(dp), dimension(0:,1:,1:),    intent(inout) :: fg_param_map
    integer(i4b),                     intent(inout) :: stat

    integer(i4b) :: i, j, k, l, q, nstep, nsub, cmb_comp
    real(dp)     :: p_old, p_prop, eta, corr_old, corr_prop
    
    nstep = 3
    nsub  = 3

    ! Find CMB component
    do i = 1, num_fg_comp
       if (trim(fg_components(i)%type) == 'cmb') cmb_comp = i
    end do

    ! Perform a Monte Carlo search for minimum cross-correlations
    do i = 1, num_fg_comp
       do j = 1, fg_components(i)%npar
          if (.not. fg_components(i)%optimize_prior(j) .or. fg_components(i)%gauss_prior(j,2) <= 0.d0) cycle
          
          corr_old = compute_crosscorr(s%fg_amp(:,:,i), s%fg_amp(:,:,cmb_comp), &
               & fg_components(i)%priormask)
          
          do l = 1, nstep
             ! Propose new prior mean
             p_old  = fg_components(i)%gauss_prior(j,1)
             p_prop = max(min(p_old + 1.d0 * rand_gauss(handle) * fg_components(i)%gauss_prior(j,2), &
                  & fg_components(i)%priors(j,2)), fg_components(i)%priors(j,1))
             fg_components(i)%gauss_prior(j,1) = p_prop

             do q = 1, nsub
                ! Sample spectral indices with new prior
                call sample_spectral_param_map(s, residuals, inv_N, s%fg_amp, fg_param_map, stat)
                call update_fg_pix_response_maps(fg_param_map)
                
                ! Sample amplitudes
                call enforce_pos_amps(chaindir, residuals, inv_N, s%fg_amp, fg_param_map, .false.) 
             end do

             ! Compute new cross-correlation coefficient
             corr_prop = compute_crosscorr(s%fg_amp(:,:,i), s%fg_amp(:,:,cmb_comp), &
                  & fg_components(i)%priormask)

             write(*,*)
             write(*,*) i, j, real(p_old,sp), real(corr_old,sp)
             write(*,*) i, j, real(p_prop,sp), real(corr_prop,sp)

             if (corr_prop < corr_old) then
                ! Accept, and leave prior at proposed value
             else
                ! Reject, and set prior back to old value
                fg_components(i)%gauss_prior(j,1) = p_old
             end if

          end do
          
       end do
    end do

  end subroutine optimize_priors

  function compute_crosscorr(map1, map2, mask)
    implicit none

    real(dp), dimension(0:,1:), intent(in) :: map1, map2, mask
    real(dp)                               :: compute_crosscorr

    integer(i4b) :: i, j, k, nmaps, npix
    real(dp) :: mu1, mu2, rms1, rms2, corr, corr0

    nmaps = size(map1,2)
    corr = 0.d0
    do i = 1, nmaps
       npix = sum(mask(:,i))
       if (npix == 0) cycle
       mu1   = sum(mask(:,i)*map1(:,i)) / npix
       mu2   = sum(mask(:,i)*map2(:,i)) / npix
       rms1  = sqrt(sum((mask(:,i)*(map1(:,i)-mu1))**2)/npix)
       rms2  = sqrt(sum((mask(:,i)*(map2(:,i)-mu2))**2)/npix)
       corr0 = abs(sum(mask(:,i)*(map1(:,i)-mu1)*(map2(:,i)-mu2)) / (npix-1) / (rms1*rms2))
       corr  = max(corr, corr0)
    end do

!!$    open(58,file='corr.dat')
!!$    do i = 0, size(map1,1)-1
!!$       if (mask(i,1) > 0.5d0) then
!!$          write(58,*) map1(i,1), map2(i,1)
!!$       end if
!!$    end do
!!$    close(58)
!    call mpi_finalize(i)
!    stop

    compute_crosscorr = corr

  end function compute_crosscorr

  subroutine sample_par_rms(handle, s, fg_param_map, stat)
    implicit none

    type(planck_rng),                 intent(inout) :: handle
    type(genvec),                     intent(in)    :: s
    real(dp), dimension(0:,1:,1:),    intent(in)    :: fg_param_map
    integer(i4b),                     intent(inout) :: stat

    integer(i4b) :: i, j, k, n, status
    real(dp)     :: par, par_old, lnL, chisq_fullsky, chisq_highlat, x_i(3), lnL_bf, par_bf, p(1)
   
    if (.not. allocated(x_init)) then
       allocate(x_init(3,num_fg_comp,2))
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             if (fg_components(i)%p_rms_gauss(j,2) <= 0.d0 .or. &
                  & fg_components(i)%p_rms_uni(j,2) <= fg_components(i)%p_rms_uni(j,1)) cycle
             do k = 1, 3
                x_init(k,i,j) = (fg_components(i)%p_rms_uni(j,2)-fg_components(i)%p_rms_uni(j,1))/4.d0 * k
             end do
          end do
       end do
    end if

    allocate(index_map_rms(0:npix-1, nmaps, num_fg_par))
    index_map_rms = fg_param_map
    call allocate_genvec(s_rms)
    call genvec_set_equal(s, s_rms)

    ! Sample spectral index RMS parameters
    do i = 1, num_fg_comp
       do j = 1, fg_components(i)%npar
          if (fg_components(i)%p_rms_gauss(j,2) <= 0.d0 .or. &
               & fg_components(i)%p_rms_uni(j,2) <= fg_components(i)%p_rms_uni(j,1)) cycle
          comp_rms = i
          par_rms  = j

          par_old = fg_components(i)%p_rms(j)
          x_i = x_init(:,i,j)
          call QuickSort_real(x_i)
!          par = sample_ARS(handle, lnL_par_rms_ARS, x_i, &
!               & fg_components(i)%p_rms_uni(j,1:2), status=status, bounds_ok=.true.)
!          par = sample_InvSamp(handle, x_i, lnL_par_rms_ARS, fg_components(i)%p_rms_uni(j,1:2), status, tolerance_= 1.d-1)

          p(1) = par_old
          !call powell(p, lnL_par_rms_powell, status)
          par = p(1)

!!$          par_bf = par_old
!!$          lnL_bf = -1.d30
!!$          open(58,file='p.dat')
!!$          do k = 0, 300
!!$             par = k*0.001d0
!!$             lnL = lnL_par_rms_ARS(par)
!!$             if (lnL > lnL_bf) then
!!$                par_bf = par
!!$                lnL_bf = lnL
!!$             else
!!$                exit
!!$             end if
!!$             write(58,*) par, lnL
!!$             write(*,*) par, lnL
!!$          end do
!!$          close(58)
!!$
!!$          par = par_bf
!!$          write(*,*) par_old, par
!!$
          write(*,*) 'foer', real(x_init(:,i,j),sp), status
          if (status == 0) then
             call update_par_rms(comp_rms, par_rms, par)
             x_init(1:2,i,j) = x_init(2:3,i,j)
             x_init(3,  i,j) = par
          else
             call update_par_rms(comp_rms, par_rms, par_old)
          end if
          call update_all_eff_fg_spectrum(comp_rms)
          call update_fg_pix_response_maps(index_map_rms)
          write(*,*) 'etter', real(x_init(:,i,j),sp)
       end do
    end do
   
    call deallocate_genvec(s_rms)
    deallocate(index_map_rms)
    
  end subroutine sample_par_rms


  function lnL_par_rms_powell(x)
    use healpix_types
    implicit none
    real(dp), dimension(1:), intent(in), optional :: x
    real(dp)                                      :: lnL_par_rms_powell

    lnL_par_rms_powell = -lnL_par_rms_ARS(x(1))

  end function lnL_par_rms_powell



  function lnL_par_rms_ARS(x)
    use healpix_types
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: lnL_par_rms_ARS

    real(dp)     :: chisq, chisq_prior, chisq_fullsky, chisq_highlat

    integer(i4b) :: k
    real(dp)     :: s1
    real(dp), allocatable, dimension(:,:) :: map

    ! Compute chi-square
    if (x >= fg_components(comp_rms)%p_rms_uni(par_rms,1) .and. &
         & x <= fg_components(comp_rms)%p_rms_uni(par_rms,2)) then
       call update_par_rms(comp_rms, par_rms, x)
    else
       lnL_par_rms_ARS = -1.d30
       return
    end if
    call update_all_eff_fg_spectrum(comp_rms)
    call update_fg_pix_response_maps(index_map_rms)
    call compute_chisq(map_id, s_rms, chisq_highlat=chisq)    

    ! Compute prior term
    chisq_prior = (x-fg_components(comp_rms)%p_rms_gauss(par_rms,1))**2 / &
         & fg_components(comp_rms)%p_rms_gauss(par_rms,2)**2

    ! Return result of chi-square and prior
    lnL_par_rms_ARS = -0.5d0 * (chisq + chisq_prior)

    write(*,*) real(x,sp), real(chisq,sp), real(chisq_prior,sp)
    if (abs(x-0.238235354d0) < 0.001d0 .and. chisq > 49556070.d0) then
       allocate(map(0:npix-1,1))
       call compute_chisq(map_id, s_rms, chisq_highlat=chisq, chisq_map=map)    
       call write_map('chisq.fits', map)

       open(58,file='s.dat')
       do k = 1, 19
          s1 = get_effective_fg_spectrum(fg_components(5), k, [26.5d0])
          if (trim(bp(k)%unit) == 'uK_cmb') then
             write(*,*) real(bp(k)%nu_c,sp), s1 / bp(k)%a2t
             write(58,*) bp(k)%nu_c/1.d9, s1 / bp(k)%a2t
          else
             write(*,*) real(bp(k)%nu_c,sp), s1 * bp(k)%f2t / bp(k)%a2t
             write(58,*) bp(k)%nu_c/1.d9, s1 * bp(k)%f2t / bp(k)%a2t
          end if
       end do
       close(58)
       stop
    end if

  end function lnL_par_rms_ARS


end module comm_global_par_mod
