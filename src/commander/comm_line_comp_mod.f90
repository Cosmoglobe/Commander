module comm_line_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_line_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_line_comp

  !**************************************************
  !           Line emission component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_line_comp
     integer(i4b)                            :: ref_band
     real(dp)                                :: line2RJ_ref
     integer(i4b), allocatable, dimension(:) :: ind2band
     real(dp),     allocatable, dimension(:) :: line2RJ
   contains
     procedure :: S    => evalSED
     procedure :: sampleSpecInd => sampleLineRatios
  end type comm_line_comp

  interface comm_line_comp
     procedure constructor
  end interface comm_line_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_line_comp), pointer   :: constructor

    integer(i4b) :: i, j, k, l, nline, b, n, ierr
    real(dp)     :: f
    logical(lgt) :: ref_exist
    character(len=512), allocatable, dimension(:) :: label
    real(dp),           allocatable, dimension(:) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:) :: poltype
    type(comm_mapinfo), pointer :: info
    
    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Read line template file
    call read_line_template(trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), &
         & nline, label, mu, sigma, line2RJ, poltype)

    ! Check how many lines are included in current run
    n         = 0
    ref_exist = .false.
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) n = n+1
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) ref_exist = .true.
    end do
    if (.not. ref_exist) call report_error("Line component reference band does not exist")

    allocate(constructor%ind2band(n))
    constructor%npar = n
    allocate(constructor%theta_def(n), constructor%p_gauss(2,n), constructor%p_uni(2,n))
    allocate(constructor%poltype(n), constructor%indlabel(n), constructor%line2RJ(n))
    n         = 0
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) then
             n = n+1
             constructor%ind2band(n)  = i
             constructor%theta_def(n) = mu(j)
             constructor%p_gauss(1,n) = mu(j)
             constructor%p_gauss(2,n) = sigma(j)
             constructor%p_uni(1,n)   = -100.d0   !mu(j)-5*sigma(j)
             constructor%p_uni(2,n)   =  100.d0   !mu(j)+5*sigma(j)
             constructor%poltype(n)   = poltype(j)
             constructor%indlabel(n)  = label(j)
             constructor%line2RJ(n)   = line2RJ(j)
             exit
          end if
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) then
          constructor%ref_band    = i
          constructor%line2RJ_ref = constructor%line2RJ(n)
       end if
    end do

    ! Initialize spectral index maps
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(n))
    do i = 1, n
       constructor%theta(i)%p     => comm_map(info)
       constructor%theta(i)%p%map = constructor%theta_def(i)
       if (constructor%lmax_ind >= 0) call constructor%theta(i)%p%YtW_scalar
    end do


    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    j = 1
    do l = 1, 3
       do i = 1, numband
          if (l > 1) then
             constructor%F_int(l,i,j)%p => constructor%F_int(l-1,i,j)%p
             cycle
          end if
          if (any(constructor%ind2band == i)) then
             do k = 0, data(i)%ndet
                constructor%F_int(l,i,k)%p => comm_F_line(constructor, data(i)%bp(k)%p, .true., &
                     & constructor%line2RJ(j) / constructor%line2RJ_ref * data(i)%RJ2data(k), j)
             end do
             j = j+1
          else
             do k = 0, data(i)%ndet
                constructor%F_int(l,i,j)%p => comm_F_line(constructor, data(i)%bp(k)%p, .false., 0.d0, j)
             end do
          end if
       end do
    end do
    
    ! Initialize mixing matrix
    if (cpar%init_samp < 0 .or. trim(cpar%init_chain_prefix) == 'none' &
         & .or. .not. constructor%init_from_HDF) &
         & call constructor%updateMixmat

    deallocate(label, mu, sigma, line2RJ, poltype)

  end function constructor

  ! Definition:
  !    SED  = delta_{band,
  function evalSED(self, nu, band, pol, theta)
    class(comm_line_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, ind

    if (band == self%ref_band) then
       evalSED = 1.d0
    else 
       do i = 1, self%npar
          if (band == self%ind2band(i)) exit
       end do
       if (i > self%npar) then
          evalSED = 0.d0
       else
          evalSED = theta(i) * self%line2RJ(i) / self%line2RJ_ref
       end if
    end if

  end function evalSED

  subroutine read_line_template(filename, nline, label, mu, sigma, line2RJ, poltype)
    implicit none
    character(len=*),                              intent(in)  :: filename
    integer(i4b),                                  intent(out) :: nline
    character(len=512), allocatable, dimension(:), intent(out) :: label
    real(dp),           allocatable, dimension(:), intent(out) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:), intent(out) :: poltype

    integer(i4b)        :: i, unit
    character(len=1024) :: line

    unit  = getlun()

    ! Find number of lines
    nline = 0
    open(unit, file=trim(filename), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
    end do
1   close(unit)
    
    allocate(label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline))
    open(unit, file=trim(filename), recl=1024)
    nline = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
       read(line,*) label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline)
    end do
2   close(unit)

  end subroutine read_line_template

  ! Sample line ratios
  subroutine sampleLineRatios(self, handle, id)
    implicit none
    class(comm_line_comp),                   intent(inout)        :: self
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id

    integer(i4b)    :: i, j, l, n, m, band, ierr
    real(dp)        :: A, b, mu, sigma, par, sigma_p, scale, w
    class(comm_map), pointer :: invN_amp, amp, mask
    character(len=2) :: id_text
    
    band = self%ind2band(id)
    !if (band == self%ref_band) return

    ! Construct mask
    if (associated(self%indmask)) then
       if (data(band)%info%nside /= self%indmask%info%nside) then
          call report_error("Mask udgrade in line_comp not yet supported.")
       else
          mask => self%indmask
       end if
    end if
    
    ! Compute likelihood term
    w            = self%theta(id)%p%map(1,1)
    amp          => comm_map(data(band)%info)
    invN_amp     => comm_map(data(band)%info)
    amp%map      =  self%getBand(band)/w
    invN_amp%map = amp%map
    call data(band)%N%invN(invN_amp)     ! Inverse noise variance weighted amplitude map
    
!!$    call int2string(id, id_text)
!!$    call mask%writeFITS('co_mask'//id_text//'.fits')
!!$    call amp%writeFITS('co_amp'//id_text//'.fits')
!!$    call data(band)%res%writeFITS('co_res'//id_text//'.fits')
!!$    call data(band)%N%invN_diag%writeFITS('co_invN'//id_text//'.fits')

    ! Reduce across processors
    if (associated(self%indmask)) then
       A = sum(invN_amp%map * mask%map * amp%map)
       b = sum(invN_amp%map * mask%map * data(band)%res%map)
    else
       A = sum(invN_amp%map * amp%map)
       b = sum(invN_amp%map * data(band)%res%map)
    end if
    call mpi_allreduce(MPI_IN_PLACE, A, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    
    call amp%dealloc()
    call invN_amp%dealloc()
    
    ! Compute new line ratio; just root processor
    if (self%x%info%myid == 0) then

!       write(*,*) 'A,b = ', A, b

       if (A > 0.d0) then
          mu    = b / A
          sigma = sqrt(1.d0 / A)
       else if (self%p_gauss(2,id) > 0.d0) then
          mu    = 0.d0
          sigma = 1.d30
       else
          mu    = self%p_uni(1,id) + (self%p_uni(2,id)-self%p_uni(1,id))*rand_uni(handle)
          sigma = 0.d0
       end if

!       write(*,*) '  mu, sigma = ', mu, sigma
       
       ! Add prior
       if (self%p_gauss(2,id) > 0.d0) then
          sigma_p = self%p_gauss(2,id) !/ sqrt(real(npix_reg,dp))
          mu      = (mu*sigma_p**2 + self%p_gauss(1,id) * sigma**2) / (sigma_p**2 + sigma**2)
          sigma   = sqrt(sigma**2 * sigma_p**2 / (sigma**2 + sigma_p**2))
       end if

!       write(*,*) '  mu_prior, sigma_prior = ', mu, sigma
       
       ! Draw sample
       par = -1.d30
       if (trim(self%operation) == 'optimize') then
          if (mu < self%p_uni(1,id)) then
             par = self%p_uni(1,id)
          else if (mu > self%p_uni(2,id)) then
             par = self%p_uni(2,id)
          else
             par = mu
          end if
       else
          do while (par < self%p_uni(1,id) .or. par > self%p_uni(2,id))
             if (mu < self%p_uni(1,id)) then
                par = rand_trunc_gauss(handle, mu, self%p_uni(1,id), sigma)
             else if (mu > self%p_uni(2,id)) then
                par = 2.d0*mu-rand_trunc_gauss(handle, mu, 2.d0*mu-self%p_uni(2,id), sigma)
             else
                par = mu + sigma * rand_gauss(handle)
             end if
          end do
       end if
       
       write(*,*) '  Line ratio i = ', id, ' = ', par
    end if
    
    ! Distribute new relative line ratio, and update
    call mpi_bcast(par, 1, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)

    if (band == self%ref_band) then
       self%x%map = self%x%map * par  ! Rescale amplitude map, but leave mixing matrix
       self%x%alm = self%x%alm * par  
       do i = 1, self%npar            ! Rescale line ratios at other frequencies
          if (self%ind2band(i) == self%ref_band) cycle
          self%theta(i)%p%map = self%theta(i)%p%map / par
          if (self%lmax_ind >= 0) then
             self%theta(i)%p%alm(:,1) = self%theta(i)%p%alm(:,1) / par
          end if
       end do
    else
       self%theta(id)%p%map = par
       if (self%lmax_ind >= 0) then
          self%theta(id)%p%alm(:,1) = par * sqrt(4.d0*pi)
       end if
    end if
    call self%updateMixmat()

  end subroutine sampleLineRatios
  
end module comm_line_comp_mod
