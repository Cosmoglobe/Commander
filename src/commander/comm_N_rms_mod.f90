module comm_N_rms_mod
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N_rms, comm_N_rms_ptr
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: siN
     class(comm_map), pointer :: siN_lowres
     class(comm_map), pointer :: rms0
   contains
     ! Data procedures

     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map_lowres
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS
     procedure :: rms_pix     => returnRMSpix
     procedure :: update_N    => update_N_rms
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

  type comm_N_rms_ptr
     type(comm_N_rms), pointer :: p
  end type comm_N_rms_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask)
    implicit none
    class(comm_N_rms),                  pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask

    integer(i4b)       :: i, ierr, tmp, nside_smooth
    real(dp)           :: sum_tau, sum_tau2, sum_noise, npix, t1, t2
    character(len=512) :: dir
    type(comm_mapinfo), pointer :: info_smooth

    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type              = cpar%ds_noise_format(id_abs)
    constructor%nmaps             = info%nmaps
    constructor%pol               = info%nmaps == 3
    constructor%uni_fsky          = cpar%ds_noise_uni_fsky(id_abs)
    constructor%set_noise_to_mean = cpar%set_noise_to_mean
    constructor%cg_precond        = cpar%cg_precond
    constructor%info              => info
    if (id_smooth == 0) then
       constructor%nside        = info%nside
       constructor%nside_lowres = min(info%nside,128)
       constructor%np           = info%np
       if (cpar%ds_regnoise(id_abs) /= 'none') then
          constructor%rms_reg => comm_map(constructor%info, trim(dir)//'/'//trim(cpar%ds_regnoise(id_abs)))
       end if
       if (present(procmask)) then
          call constructor%update_N(info, handle, mask, regnoise, procmask=procmask, &
               & noisefile=trim(dir)//trim(cpar%ds_noisefile(id_abs)))
       else
          call constructor%update_N(info, handle, mask, regnoise, &
               & noisefile=trim(dir)//trim(cpar%ds_noisefile(id_abs)))
       end if
    else
       tmp         =  getsize_fits(trim(dir)//trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)), nside=nside_smooth)
       info_smooth => comm_mapinfo(info%comm, nside_smooth, cpar%lmax_smooth(id_smooth), &
            & constructor%nmaps, constructor%pol)
       constructor%nside   = info_smooth%nside
       constructor%np      = info_smooth%np
       constructor%siN     => comm_map(info_smooth, trim(dir)//trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)))

       where (constructor%siN%map > 0.d0) 
          constructor%siN%map = 1.d0 / constructor%siN%map
       elsewhere
          constructor%siN%map = 0.d0
       end where

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

    end if

  end function constructor


  subroutine update_N_rms(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_rms),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, ierr
    real(dp)     :: sum_tau, sum_tau2, sum_noise, npix, t1, t2
    class(comm_map),     pointer :: invW_tau, iN
    class(comm_mapinfo), pointer :: info_lowres

    if (present(noisefile)) then
       self%rms0     => comm_map(mask%info, noisefile)
    else
       self%rms0%map = map%map
    end if
    if (associated(self%rms_reg)) then
       self%rms0%map = sqrt(self%rms0%map**2 + self%rms_reg%map**2) 
    end if
    self%siN     => comm_map(self%rms0)
    call uniformize_rms(handle, self%siN, self%uni_fsky, mask, regnoise)
    self%siN%map = self%siN%map * mask%map ! Apply mask
    if (present(procmask)) then
       where (procmask%map < 0.5d0)
          self%siN%map = self%siN%map * 20.d0 ! Boost noise by 20 in processing mask
       end where
    end if

    where (self%siN%map > 0.d0) 
       self%siN%map = 1.d0 / self%siN%map
    elsewhere
       self%siN%map = 0.d0
    end where

    ! Set siN to its mean; useful for debugging purposes
    if (self%set_noise_to_mean) then
       do i = 1, self%nmaps
          sum_noise = sum(self%siN%map(:,i))
          npix      = size(self%siN%map(:,i))
          call mpi_allreduce(MPI_IN_PLACE, sum_noise,  1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, npix,       1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
          self%siN%map(:,i) = sum_noise/npix
       end do
    end if

    if (trim(self%cg_precond) == 'diagonal') then
       ! Set up diagonal covariance matrix
       if (.not. associated(self%invN_diag)) self%invN_diag => comm_map(mask%info)
       self%invN_diag%map = self%siN%map**2
       call compute_invN_lm(self%invN_diag)
    else if (trim(self%cg_precond) == 'pseudoinv') then
       ! Compute alpha_nu for pseudo-inverse preconditioner
       if (.not. allocated(self%alpha_nu)) allocate(self%alpha_nu(self%nmaps))
       invW_tau     => comm_map(self%siN)
       invW_tau%map =  invW_tau%map**2
       call invW_tau%Yt()
       call invW_tau%Y()
       ! Temperature
       sum_tau  = sum(invW_tau%map(:,1))
       sum_tau2 = sum(invW_tau%map(:,1)**2)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
       if (sum_tau > 0.d0) then
          self%alpha_nu(1) = sqrt(sum_tau2/sum_tau)
       else
          self%alpha_nu(1) = 0.d0
       end if

       if (self%nmaps == 3) then
          sum_tau  = sum(invW_tau%map(:,2:3))
          sum_tau2 = sum(invW_tau%map(:,2:3)**2)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mask%info%comm, ierr)
          if (sum_tau > 0.d0) then
             self%alpha_nu(2:3) = sqrt(sum_tau2/sum_tau)
          else
             self%alpha_nu(2:3) = 0.d0
          end if
          call invW_tau%dealloc()
       end if
    end if

    ! Set up lowres map
    if (.not.associated(self%siN_lowres)) then
       info_lowres => comm_mapinfo(self%info%comm, self%nside_lowres, 0, self%nmaps, self%pol)
       self%siN_lowres => comm_map(info_lowres)
    end if
    iN => comm_map(self%siN)
    iN%map = iN%map**2
    call iN%udgrade(self%siN_lowres)
    call iN%dealloc()
    self%siN_lowres%map = sqrt(self%siN_lowres%map) * (self%nside/self%nside_lowres)

  end subroutine update_N_rms

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    map%map = (self%siN%map)**2 * map%map
  end subroutine matmulInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_1map_lowres(self, map)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    map%map = (self%siN_lowres%map)**2 * map%map
  end subroutine matmulInvN_1map_lowres

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

  subroutine uniformize_rms(handle, rms, fsky, mask, regnoise)
    implicit none
    type(planck_rng),                   intent(inout) :: handle
    class(comm_map),                    intent(inout) :: rms
    real(dp),                           intent(in)    :: fsky
    class(comm_map),                    intent(in)    :: mask
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
       if (all(mask%map(:,j) < 0.5d0)) cycle
       ! Find pixel histogram across cores
       limits(1) = minval(rms%map(:,j), mask%map(:,j) > 0.5d0)
       limits(2) = maxval(rms%map(:,j), mask%map(:,j) > 0.5d0)
       call mpi_allreduce(MPI_IN_PLACE, limits(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, rms%info%comm, ierr)       
       call mpi_allreduce(MPI_IN_PLACE, limits(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, rms%info%comm, ierr)       
       dx = (limits(2)-limits(1))/nbin
       if (dx == 0.d0) then
          regnoise(:,j) = 0.d0
          cycle
       end if
       F = 0.d0
       do i = 0, rms%info%np-1
          if (mask%map(i,j) <= 0.5d0) cycle
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
          if (rms%map(i,j) < threshold .and. mask%map(i,j) > 0.5d0) then
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
