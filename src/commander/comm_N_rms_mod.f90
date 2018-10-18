module comm_N_rms_mod
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N_rms, comm_N_rms_ptr
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: siN
   contains
     ! Data procedures
     procedure :: invN     => matmulInvN_1map
     procedure :: N        => matmulN_1map
     procedure :: sqrtInvN => matmulSqrtInvN_1map
     procedure :: rms      => returnRMS
     procedure :: rms_pix  => returnRMSpix
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

  type comm_N_rms_ptr
     type(comm_N_rms), pointer :: p
  end type comm_N_rms_ptr



!!$  interface matmulInvN
!!$     module procedure matmulInvN_1map, matmulInvN_2map
!!$  end interface matmulInvN
!!$  
!!$  interface matmulSqrtInvN
!!$     module procedure matmulSqrtInvN_1map, matmulSqrtInvN_2map
!!$  end interface matmulSqrtInvN
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask)
    implicit none
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out)   :: regnoise
    class(comm_N_rms),                  pointer       :: constructor
    class(comm_map),                    pointer, intent(in), optional :: procmask

    integer(i4b)       :: i, ierr, tmp, nside_smooth
    real(dp)           :: sum_tau, sum_tau2, sum_noise, npix, t1, t2
    character(len=512) :: dir, cache
    character(len=4)   :: itext
    type(comm_mapinfo), pointer :: info_smooth
    class(comm_map),    pointer :: invW_tau
    
    ! General parameters
    allocate(constructor)
    call int2string(info%myid, itext)
    dir   = trim(cpar%datadir) // '/'
    cache = trim(dir) // 'invNlm_' // trim(cpar%ds_label(id)) // '_proc' // itext // '.unf'

    ! Component specific parameters
    constructor%type    = cpar%ds_noise_format(id_abs)
    constructor%nmaps   = info%nmaps
    constructor%pol     = info%nmaps == 3
!    call wall_time(t1)
!    constructor%siN     => comm_map(info, trim(dir)//trim(cpar%ds_noise_rms(id_abs)))
!    call wall_time(t2)
!    if (info%myid == 0) write(*,*) 'read = ', t2-t1
    if (id_smooth == 0) then
       constructor%nside   = info%nside
       constructor%np      = info%np
       constructor%siN     => comm_map(info, trim(dir)//trim(cpar%ds_noise_rms(id_abs)))
       call wall_time(t1)
       call uniformize_rms(handle, constructor%siN, cpar%ds_noise_uni_fsky(id_abs), regnoise)
       call wall_time(t2)
       if (info%myid == 0) write(*,*) 'uniformize = ', t2-t1
       call wall_time(t1)
       constructor%siN%map = constructor%siN%map * mask%map ! Apply mask
       call wall_time(t2)
       if (info%myid == 0) write(*,*) 'apply_mask = ', t2-t1
       if (present(procmask)) then
          where (procmask%map < 0.5d0)
             constructor%siN%map = constructor%siN%map * 20.d0 ! Boost noise by 20 in processing mask
          end where
       end if
    else
       tmp         =  getsize_fits(trim(dir)//trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)), nside=nside_smooth)
       info_smooth => comm_mapinfo(info%comm, nside_smooth, cpar%lmax_smooth(id_smooth), &
            & constructor%nmaps, constructor%pol)
       constructor%nside   = info_smooth%nside
       constructor%np      = info_smooth%np
       call wall_time(t1)
       constructor%siN     => comm_map(info_smooth, trim(dir)//trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)))
       call wall_time(t2)
       if (info%myid == 0) write(*,*) 'read = ', t2-t1
    end if
    call wall_time(t1)
    where (constructor%siN%map > 0.d0) 
       constructor%siN%map = 1.d0 / constructor%siN%map
    elsewhere
       constructor%siN%map = 0.d0
    end where
    call wall_time(t2)
    if (info%myid == 0) write(*,*) 'siN = ', t2-t1

    ! Set siN to its mean; useful for debugging purposes
    if (cpar%set_noise_to_mean) then
       do i = 1, constructor%nmaps
          sum_noise = sum(constructor%siN%map(:,i))
          npix      = size(constructor%siN%map(:,i))
          call mpi_allreduce(MPI_IN_PLACE, sum_noise,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, npix,       1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          constructor%siN%map(:,i) = sum_noise/npix
       end do
    end if

    ! Set up diagonal covariance matrix
    if (id_smooth == 0 .and. trim(cpar%cg_precond) == 'diagonal') then
       constructor%invN_diag     => comm_map(info)
       constructor%invN_diag%map = constructor%siN%map**2
       call compute_invN_lm(cache, constructor%invN_diag)
    end if

    ! Compute alpha_nu for pseudo-inverse preconditioner
    call wall_time(t1)
    if (id_smooth == 0 .and. trim(cpar%cg_precond) == 'pseudoinv') then
       allocate(constructor%alpha_nu(constructor%nmaps))
       invW_tau     => comm_map(constructor%siN)
       invW_tau%map =  invW_tau%map**2
       call invW_tau%Yt()
       call invW_tau%Y()
       ! Temperature
       sum_tau  = sum(invW_tau%map(:,1))
       sum_tau2 = sum(invW_tau%map(:,1)**2)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       if (sum_tau > 0.d0) then
          constructor%alpha_nu(1) = sqrt(sum_tau2/sum_tau)
       else
          constructor%alpha_nu(1) = 0.d0
       end if

       if (constructor%nmaps == 3) then
          sum_tau  = sum(invW_tau%map(:,2:3))
          sum_tau2 = sum(invW_tau%map(:,2:3)**2)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          if (sum_tau > 0.d0) then
             constructor%alpha_nu(2:3) = sqrt(sum_tau2/sum_tau)
          else
             constructor%alpha_nu(2:3) = 0.d0
          end if
       end if
    end if
    call wall_time(t2)
    if (info%myid == 0) write(*,*) 'precond = ', t2-t1

  end function constructor

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    map%map = (self%siN%map)**2 * map%map
  end subroutine matmulInvN_1map

  ! Return map_out = N * map
  subroutine matmulN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    where (self%siN%map > 0.d0)
       map%map = map%map / (self%siN%map)**2 
    elsewhere
       map%map = 0.d0
    end where
  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    map%map = self%siN%map * map%map
  end subroutine matmulSqrtInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_2map(self, map, res)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(in)              :: map
    class(comm_map),   intent(inout)           :: res
    res%map = (self%siN%map)**2 * map%map
  end subroutine matmulInvN_2map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_2map(self, map, res)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(in)              :: map
    class(comm_map),   intent(inout)           :: res
    res%map = self%siN%map * map%map
  end subroutine matmulSqrtInvN_2map

  ! Return RMS map
  subroutine returnRMS(self, res)
    implicit none
    class(comm_N_rms), intent(in)    :: self
    class(comm_map),   intent(inout) :: res
    where (self%siN%map > 0.d0)
       res%map = 1.d0/self%siN%map
    elsewhere
       res%map = infinity
    end where
  end subroutine returnRMS
  
  ! Return rms for single pixel
  function returnRMSpix(self, pix, pol)
    implicit none
    class(comm_N_rms),   intent(in)    :: self
    integer(i4b),        intent(in)    :: pix, pol
    real(dp)                           :: returnRMSpix
    if (self%siN%map(pix,pol) > 0.d0) then
       returnRMSpix = 1.d0/self%siN%map(pix,pol)
    else
       returnRMSpix = infinity
    end if
  end function returnRMSpix

  subroutine uniformize_rms(handle, rms, fsky, regnoise)
    implicit none
    type(planck_rng),                   intent(inout) :: handle
    class(comm_map),                    intent(inout) :: rms
    real(dp),                           intent(in)    :: fsky
    real(dp),         dimension(0:,1:), intent(out)   :: regnoise

    integer(i4b) :: i, j, nbin=1000, ierr, b
    real(dp)     :: limits(2), dx, threshold, sigma
    real(dp), allocatable, dimension(:) :: F

    if (fsky <= 0.d0) then
       regnoise = 0.d0
       return
    end if

    allocate(F(nbin))
    do j = 1, rms%info%nmaps
       ! Find pixel histogram across cores
       limits(1) = minval(rms%map(:,j))
       limits(2) = maxval(rms%map(:,j))
       call mpi_allreduce(MPI_IN_PLACE, limits(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, rms%info%comm, ierr)       
       call mpi_allreduce(MPI_IN_PLACE, limits(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, rms%info%comm, ierr)       
       dx = (limits(2)-limits(1))/nbin
       if (dx == 0.d0) cycle
       F = 0.d0
       do i = 0, rms%info%np-1
          b    = max(min(int((rms%map(i,j)-limits(1))/dx),nbin),1)
          F(b) = F(b) + 1.d0
       end do
       call mpi_allreduce(MPI_IN_PLACE, F, nbin, MPI_DOUBLE_PRECISION, MPI_SUM, rms%info%comm, ierr)

       ! Compute cumulative distribution
       do i = 2, nbin
          F(i) = F(i-1) + F(i)
       end do
       F = F / maxval(F)

       ! Find threshold
       i = 1
       do while (F(i) < fsky)
          i = i+1
       end do
       threshold = limits(1) + dx*(i-1)

       ! Update RMS map, and draw corresponding noise realization
       do i = 0, rms%info%np-1
          if (rms%map(i,j) < threshold) then
             sigma         = sqrt(threshold**2 - rms%map(i,j)**2)
             rms%map(i,j)  = threshold                  ! Update RMS map to requested limit
             regnoise(i,j) = sigma * rand_gauss(handle) ! Draw corresponding noise realization
          else
             regnoise(i,j) = 0.d0
          end if
       end do
    end do
    deallocate(F)

  end subroutine uniformize_rms

end module comm_N_rms_mod
