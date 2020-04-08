module comm_N_QUcov_mod
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N_QUcov, comm_N_QUcov_ptr
  
  type, extends (comm_N) :: comm_N_QUcov
     real(dp),        allocatable, dimension(:,:) :: Ncov
     real(dp),        allocatable, dimension(:,:) :: siN
     real(dp),        allocatable, dimension(:,:) :: iN
     class(comm_map), pointer                     :: siN_diag => null()
   contains
     ! Data procedures
     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS
     procedure :: rms_pix     => returnRMSpix
     procedure :: update_N    => update_N_QUcov
  end type comm_N_QUcov

  interface comm_N_QUcov
     procedure constructor
  end interface comm_N_QUcov

  type comm_N_QUcov_ptr
     type(comm_N_QUcov), pointer :: p => null()
  end type comm_N_QUcov_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask)
    implicit none
    class(comm_N_QUcov),                pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask

    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type              = cpar%ds_noise_format(id_abs)
    constructor%nmaps             = 3
    constructor%pol               = .true.
    constructor%uni_fsky          = 1.d0 
    constructor%set_noise_to_mean = .false.
    constructor%cg_precond        = cpar%cg_precond
    constructor%nside             = info%nside
    constructor%nside_lowres      = info%nside
    constructor%npix              = 12*info%nside**2
    constructor%np                = info%np
    constructor%myid              = info%myid
    constructor%comm              = info%comm
    constructor%nprocs            = info%nprocs
    constructor%info              => info
    constructor%pol_only          = .true.
    call constructor%update_N(info, handle, mask=mask, noisefile=trim(dir)//trim(cpar%ds_noisefile(id_abs)))

  end function constructor

  subroutine update_N_QUcov(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_QUcov),                 intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, ierr, status, unit
    logical(lgt) :: exist
    character(len=512) :: filename
    real(dp)     :: sum_tau, sum_tau2, val
    real(sp),        dimension(:,:), pointer     :: Ninv_sp => null()
    real(dp),        dimension(:,:), allocatable :: Ninv, Ncov, sNinv, buffer, mask_fullsky
    class(comm_map),                 pointer     :: invW_tau => null()

    if (allocated(self%iN)) return

    ! Set up dense QU invN covariance matrix
    allocate(self%iN(2*self%npix,2*self%np), self%siN(2*self%npix,2*self%np), self%Ncov(2*self%npix,2*self%np))
    allocate(buffer(2*self%npix,2*self%npix))

    filename = trim(noisefile)//'_precomp.unf'
    inquire(file=trim(filename), exist=exist)

    if (.not. exist) then
       allocate(mask_fullsky(0:self%npix-1,self%nmaps))
       call mask%bcast_fullsky_map(mask_fullsky)
    end if

    if (self%myid == 0) then
       allocate(Ninv(2*self%npix,2*self%npix))
       allocate(sNinv(2*self%npix,2*self%npix))
       allocate(Ncov(2*self%npix,2*self%npix))

       unit = getlun()
       if (exist) then
          write(*,*) '   Reading precomputed matrices from ', trim(filename)
          open(unit,file=trim(filename),form='unformatted')
          read(unit) Ninv
          read(unit) sNinv
          read(unit) Ncov
          close(unit)
       else
          write(*,*) '   Eigen-decomposing ', trim(noisefile)
          allocate(Ninv_sp(2*self%npix,2*self%npix))
          call WMAP_Read_NInv(noisefile, status, Ninv_sp)

          ! Convert from nest to ring format
          do i = 1, 2*self%npix ! Rows
             call convert_nest2ring(self%nside, Ninv_sp(          1:  self%npix,i))
             call convert_nest2ring(self%nside, Ninv_sp(self%npix+1:2*self%npix,i))
          end do
          do i = 1, 2*self%npix ! Columns
             call convert_nest2ring(self%nside, Ninv_sp(i,          1:  self%npix))
             call convert_nest2ring(self%nside, Ninv_sp(i,self%npix+1:2*self%npix))
          end do
          Ncov = Ninv_sp
          deallocate(Ninv_sp)

          ! Compute N mask and apply Q mask to both fields
          call invert_matrix_with_mask(Ncov)
          do i = 1, self%npix
             if (mask_fullsky(i-1,2) == 0.d0) then
                Ncov(i,:)           = 0.d0
                Ncov(:,i)           = 0.d0
                Ncov(i+self%npix,:) = 0.d0
                Ncov(:,i+self%npix) = 0.d0
             end if
          end do

          Ninv = Ncov
          call compute_hermitian_root_with_mask(Ninv, -1.d0, sNinv, -0.5d0)
          open(unit,file=trim(filename),form='unformatted')
          write(unit) Ninv
          write(unit) sNinv
          write(unit) Ncov
          close(unit)
       end if
    end if
    if (.not. exist) deallocate(mask_fullsky)

    if (self%myid == 0) buffer = Ninv
    call mpi_bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    do i = 1, 2
       self%iN(:,(i-1)*self%np+1:i*self%np) = buffer(:,(i-1)*self%npix+info%pix+1)
    end do

    if (self%myid == 0) buffer = sNinv
    call mpi_bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    do i = 1, 2
       self%siN(:,(i-1)*self%np+1:i*self%np) = buffer(:,(i-1)*self%npix+info%pix+1)
    end do

    if (self%myid == 0) buffer = Ncov
    call mpi_bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    do i = 1, 2
       self%Ncov(:,(i-1)*self%np+1:i*self%np) = buffer(:,(i-1)*self%npix+info%pix+1)
    end do
    if (self%myid == 0) deallocate(Ninv, Ncov, sNinv)
    deallocate(buffer)


    ! Set up diagonal invN
    self%siN_diag => comm_map(info)
    self%siN_diag%map(:,1) = 0.d0 ! Set temperature to zero
    do i = 1, info%np
       val = self%Ncov(info%pix(i)+1,i) ! N(Q,Q)
       if (val > 0.d0) then
          self%siN_diag%map(i-1,2) = 1.d0/sqrt(val)
       else
          self%siN_diag%map(i-1,2) = 0.d0
       end if
       val = self%Ncov(info%npix+info%pix(i)+1,info%np+i) ! N(U,U)
       if (val > 0.d0) then
          self%siN_diag%map(i-1,3) = 1.d0/sqrt(val)
       else
          self%siN_diag%map(i-1,3) = 0.d0
       end if
    end do

!!$    call self%siN_diag%writeFITS('test.fits')
!!$    stop

    if (trim(self%cg_precond) == 'diagonal') then
       ! Set up diagonal covariance matrix
       if (.not. associated(self%invN_diag)) self%invN_diag => comm_map(info)
       self%invN_diag%map = self%siN_diag%map**2
       call compute_invN_lm(self%invN_diag)
    else if (trim(self%cg_precond) == 'pseudoinv') then
       ! Compute alpha_nu for pseudo-inverse preconditioner
       if (.not. allocated(self%alpha_nu)) allocate(self%alpha_nu(self%nmaps))
       invW_tau     => comm_map(self%siN_diag)
       invW_tau%map =  invW_tau%map**2
       call invW_tau%Yt()
       call invW_tau%Y()
       ! Temperature
       sum_tau  = sum(invW_tau%map(:,1))
       sum_tau2 = sum(invW_tau%map(:,1)**2)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       if (sum_tau > 0.d0) then
          self%alpha_nu(1) = sqrt(sum_tau2/sum_tau)
       else
          self%alpha_nu(1) = 0.d0
       end if

       ! Polarization
       sum_tau  = sum(invW_tau%map(:,2:3))
       sum_tau2 = sum(invW_tau%map(:,2:3)**2)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       if (sum_tau > 0.d0) then
          self%alpha_nu(2:3) = sqrt(sum_tau2/sum_tau)
       else
          self%alpha_nu(2:3) = 0.d0
       end if
       call invW_tau%dealloc()
    end if

  end subroutine update_N_QUcov



  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: map
    
    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m

    allocate(m(2*self%np), invN_m(2*self%npix))
    m(1:self%np)           = map%map(:,2)
    m(self%np+1:2*self%np) = map%map(:,3)
    invN_m                 = matmul(self%iN, m)
    call mpi_allreduce(MPI_IN_PLACE, invN_m, 2*self%npix, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)       
    map%map(:,1) = 0.d0
    map%map(:,2) = invN_m(self%info%pix+1)
    map%map(:,3) = invN_m(self%npix+self%info%pix+1)
    deallocate(m, invN_m)

  end subroutine matmulInvN_1map


  ! Return map_out = N * map
  subroutine matmulN_1map(self, map)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: map

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m

    allocate(m(2*self%np), invN_m(2*self%npix))
    m(1:self%np)           = map%map(:,2)
    m(self%np+1:2*self%np) = map%map(:,3)
    invN_m                 = matmul(self%Ncov, m)
    call mpi_allreduce(MPI_IN_PLACE, invN_m, 2*self%npix, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)       
    map%map(:,1) = 0.d0
    map%map(:,2) = invN_m(self%info%pix+1)
    map%map(:,3) = invN_m(self%npix+self%info%pix+1)
    deallocate(m, invN_m)

  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m

    allocate(m(2*self%np), invN_m(2*self%npix))
    m(1:self%np)           = map%map(:,2)
    m(self%np+1:2*self%np) = map%map(:,3)
    invN_m                 = matmul(self%siN, m)
    call mpi_allreduce(MPI_IN_PLACE, invN_m, 2*self%npix, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)       
    map%map(:,1) = 0.d0
    map%map(:,2) = invN_m(self%info%pix+1)
    map%map(:,3) = invN_m(self%npix+self%info%pix+1)
    deallocate(m, invN_m)

  end subroutine matmulSqrtInvN_1map

!!$  ! Return map_out = invN * map
!!$  subroutine matmulInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$
!!$    integer(i4b) :: ierr
!!$    real(dp), allocatable, dimension(:) :: m, invN_m
!!$
!!$    allocate(m(2*self%np), invN_m(2*self%npix))
!!$    m(1:self%np)           = map%map(:,2)
!!$    m(self%np+1:2*self%np) = map%map(:,3)
!!$    invN_m                 = matmul(self%iN, m)
!!$    call mpi_allreduce(MPI_IN_PLACE, invN_m, 2*self%npix, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)       
!!$    res%map(:,1) = 0.d0
!!$    res%map(:,2) = invN_m(self%info%pix+1)
!!$    res%map(:,3) = invN_m(self%npix+self%info%pix+1)
!!$    deallocate(m, invN_m)
!!$
!!$  end subroutine matmulInvN_2map
!!$  
!!$  ! Return map_out = sqrtInvN * map
!!$  subroutine matmulSqrtInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$
!!$    integer(i4b) :: ierr
!!$    real(dp), allocatable, dimension(:) :: m, invN_m
!!$
!!$    allocate(m(2*self%np), invN_m(2*self%npix))
!!$    m(1:self%np)           = map%map(:,2)
!!$    m(self%np+1:2*self%np) = map%map(:,3)
!!$    invN_m                 = matmul(self%siN, m)
!!$    call mpi_allreduce(MPI_IN_PLACE, invN_m, 2*self%npix, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr) 
!!$
!!$    res%map(:,1) = 0.d0      
!!$    res%map(:,2) = invN_m(self%info%pix+1)
!!$    res%map(:,3) = invN_m(self%npix+self%info%pix+1)
!!$    deallocate(m, invN_m)
!!$
!!$  end subroutine matmulSqrtInvN_2map

  ! Return RMS map
  subroutine returnRMS(self, res)
    implicit none
    class(comm_N_QUcov), intent(in)    :: self
    class(comm_map),     intent(inout) :: res
    where (self%siN_diag%map > 0.d0)
       res%map = 1.d0/self%siN_diag%map
    elsewhere
       res%map = infinity
    end where
  end subroutine returnRMS
  
  ! Return rms for single pixel
  function returnRMSpix(self, pix, pol)
    implicit none
    class(comm_N_QUcov),   intent(in)    :: self
    integer(i4b),          intent(in)    :: pix, pol
    real(dp)                             :: returnRMSpix
    if (self%siN_diag%map(pix,pol) > 0.d0) then
       returnRMSpix = 1.d0/self%siN_diag%map(pix,pol)
    else
       returnRMSpix = infinity
    end if
  end function returnRMSpix

end module comm_N_QUcov_mod
