module comm_Cl_mod
  use comm_param_mod
  use comm_map_mod
  use math_tools
  use comm_hdf_mod
  use comm_bp_utils
  implicit none

  private
  public comm_Cl, write_sigma_l

  integer(i4b), parameter :: TT = 1
  integer(i4b), parameter :: TE = 2
  integer(i4b), parameter :: TB = 3
  integer(i4b), parameter :: EE = 4
  integer(i4b), parameter :: EB = 5
  integer(i4b), parameter :: BB = 6
  logical(lgt), private   :: only_pol
  
  type :: comm_Cl
     ! General parameters
     class(comm_mapinfo), pointer :: info
     character(len=512)           :: type  ! {none, binned, power_law, exp}
     character(len=512)           :: label ! {none, binned, power_law, exp}
     character(len=512)           :: unit
     character(len=512)           :: outdir
     integer(i4b)                 :: lmax, nmaps, nspec, l_apod
     integer(i4b)                 :: poltype  ! {1 = {T+E+B}, 2 = {T,E+B}, 3 = {T,E,B}}
     real(dp)                     :: nu_ref(3), RJ2unit(3)
     real(dp),         allocatable, dimension(:,:)   :: Dl
     real(dp),         allocatable, dimension(:,:,:) :: sqrtS_mat, S_mat, sqrtInvS_mat

     ! Bin parameters
     integer(i4b) :: nbin
     character(len=1), allocatable, dimension(:,:) :: stat
     integer(i4b), allocatable, dimension(:,:)     :: bins
     
     ! Power law/exponential parameters
     character(len=512) :: plfile
     integer(i4b)       :: iter = 0
     real(dp)           :: prior(2), lpiv
     real(dp), allocatable, dimension(:) :: amp, beta
   contains
     ! Data procedures
     procedure :: S        => matmulS
     procedure :: sqrtS    => matmulSqrtS
     procedure :: sqrtInvS => matmulSqrtInvS
     procedure :: sampleCls
     procedure :: read_binfile
     procedure :: read_Cl_file
     procedure :: binCls
     procedure :: writeFITS
     procedure :: updatePowlaw
     procedure :: updateExponential
     procedure :: updateGaussian
     procedure :: updateS
     procedure :: getCl
  end type comm_Cl

  interface comm_Cl
     procedure constructor
  end interface comm_Cl


contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    type(comm_mapinfo), target,         intent(in) :: info
    integer(i4b),                       intent(in) :: id, id_abs
    class(comm_Cl),                     pointer    :: constructor

    integer(i4b)       :: l, nmaps
    character(len=512) :: datadir, binfile
    logical(lgt)       :: pol

    allocate(constructor)
    
    ! General parameters
    constructor%type   = cpar%cs_cltype(id_abs)
    if (trim(constructor%type) == 'none') return
    
    constructor%info   => info
    constructor%label  = cpar%cs_label(id_abs)
    constructor%lmax   = cpar%cs_lmax_amp(id_abs)
    constructor%l_apod = cpar%cs_l_apod(id_abs)
    constructor%unit   = cpar%cs_unit(id_abs)
    constructor%nu_ref = cpar%cs_nu_ref(id_abs,:)
    constructor%nmaps  = 1; if (cpar%cs_polarization(id_abs)) constructor%nmaps = 3
    constructor%nspec  = 1; if (cpar%cs_polarization(id_abs)) constructor%nspec = 6
    constructor%outdir = cpar%outdir
    datadir            = cpar%datadir
    nmaps              = constructor%nmaps
    only_pol           = cpar%only_pol

    ! Set up conversion factor between RJ and native component unit
    ! D_l is defined in component units, while S, invS etc are defined in RJ
    do l = 1, 3
       select case (trim(constructor%unit))
       case ('uK_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l))
       case ('mK_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l)) * 1d-3
       case ('K_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l)) * 1d-6
       case ('MJy/sr') 
          constructor%RJ2unit(l) = comp_bnu_prime_RJ(constructor%nu_ref(l)) * 1e14
       case ('K km/s') 
          constructor%RJ2unit(l) = 1.d0
       case ('y_SZ') 
          constructor%RJ2unit(l) = 2.d0*constructor%nu_ref(l)**2*k_b/c**2 / &
               & (comp_bnu_prime(constructor%nu_ref(l)) * comp_sz_thermo(constructor%nu_ref(l)))
       case ('uK_RJ') 
          constructor%RJ2unit(l) = 1.d0
       case default
          call report_error('Unsupported unit: ' // trim(constructor%unit))
       end select
    end do

    allocate(constructor%Dl(0:constructor%lmax,constructor%nspec))
    allocate(constructor%sqrtS_mat(nmaps,nmaps,0:constructor%lmax))
    allocate(constructor%sqrtInvS_mat(nmaps,nmaps,0:constructor%lmax))
    allocate(constructor%S_mat(nmaps,nmaps,0:constructor%lmax))

    if (trim(constructor%type) == 'binned') then
       call constructor%read_binfile(trim(datadir) // '/' // trim(cpar%cs_binfile(id_abs)))
       call constructor%read_Cl_file(trim(datadir) // '/' // trim(cpar%cs_clfile(id_abs)))
       call constructor%binCls
       if (cpar%only_pol) then
          constructor%stat(:,1:3) = '0'
          constructor%Dl(:,TT)     = 0.d0
          if (constructor%nspec == 6) constructor%Dl(:,[TE,TB]) = 0.d0
       end if
    else if (trim(constructor%type) == 'power_law') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updatePowlaw(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps), cpar%only_pol)
    else if (trim(constructor%type) == 'exp') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updateExponential(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps), cpar%only_pol)
    else if (trim(constructor%type) == 'gauss') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updateGaussian(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps), cpar%only_pol)
    else
       call report_error("Unknown Cl type: " // trim(constructor%type))
    end if

    !constructor%Dl = constructor%Dl / c%RJ2unit_    ! Define prior in output units
    call constructor%updateS
    
  end function constructor


  subroutine updatePowlaw(self, amp, beta, only_pol)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in)    :: amp, beta
    logical(lgt),                  intent(in)    :: only_pol

    integer(i4b) :: i, j, l, i_min
    
    self%amp   = amp
    self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = 1, self%lmax
          self%Dl(l,j) = amp(i) * (real(l,dp)/real(self%lpiv,dp))**beta(i) 
       end do
       self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updatePowlaw

  subroutine updateExponential(self, amp, beta, only_pol)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in)    :: amp, beta
    logical(lgt),                  intent(in)    :: only_pol

    integer(i4b) :: i, j, l, i_min
    
    self%amp   = amp
    self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = 1, self%lmax
          self%Dl(l,j) = amp(i) * exp(-beta(i)*(real(l,dp)/real(self%lpiv,dp))) 
       end do
       self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updateExponential

  subroutine updateGaussian(self, amp, beta, only_pol)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in)    :: amp, beta
    logical(lgt),                  intent(in)    :: only_pol

    integer(i4b) :: i, j, l, i_min
    
    self%amp   = amp
    self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = 0, self%lmax
          self%Dl(l,j) = amp(i) * exp(-l*(l+1)*(beta(i)*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
!          if (self%info%myid == 0) then
!             write(*,*) l, j, amp(i), beta(i), exp(-l*(l+1)*(beta(i)*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2), self%Dl(l,j)
!          end if
       end do
    end do

  end subroutine updateGaussian

  subroutine updateS(self)
    implicit none
    class(comm_Cl),                intent(inout) :: self

    integer(i4b) :: i, j, k, l
    logical(lgt) :: ok(self%nmaps)

    self%S_mat        = 0.d0
    self%sqrtS_mat    = 0.d0
    self%sqrtInvS_mat = 0.d0
    do l = 0, self%lmax
       k = 1
       do i = 1, self%nmaps
          do j = i, self%nmaps
             !k = i*(1-i)/2 + (i-1)*self%nmaps + j
             if (l == 0) then
                self%sqrtS_mat(i,j,l) = self%Dl(l,k)
             else
                self%sqrtS_mat(i,j,l) = self%Dl(l,k) / (l*(l+1)/(2.d0*pi))
             end if
             if (self%RJ2unit(i)*self%RJ2unit(j) == 0.d0) write(*,*) trim(self%label), self%info%myid, l, k, i, j, real(self%RJ2unit(i),sp), real(self%RJ2unit(j),sp)
             self%sqrtS_mat(i,j,l) = self%sqrtS_mat(i,j,l) / &
                  & (self%RJ2unit(i)*self%RJ2unit(j))
             self%sqrtS_mat(j,i,l) = self%sqrtS_mat(i,j,l) 
             if (i == j) ok(i) = self%Dl(l,k) > 0.d0
             k = k+1
          end do
       end do
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtS_mat(i,:,l) = 0.d0
             self%sqrtS_mat(:,i,l) = 0.d0
             self%sqrtS_mat(i,i,l) = 1.d0
          end if
       end do

!       if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(:,:,l)
       !if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(1,1,l)*self%sqrtS_mat(2,2,l) - self%sqrtS_mat(2,1,l)*self%sqrtS_mat(1,2,l)

       ! Change to RJ units
       !self%sqrtS_mat(:,:,l) = self%sqrtS_mat(:,:,l) / self%RJ2unit**2

       self%sqrtInvS_mat(:,:,l) = self%sqrtS_mat(:,:,l)
       call compute_hermitian_root(self%sqrtS_mat(:,:,l), 0.5d0)
!       if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(:,:,l), self%sqrtInvS_mat(:,:,l)
       
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtS_mat(i,:,l) = 0.d0
             self%sqrtS_mat(:,i,l) = 0.d0
          end if
       end do
       self%S_mat(:,:,l) = matmul(self%sqrtS_mat(:,:,l), self%sqrtS_mat(:,:,l))
       
       call compute_hermitian_root(self%sqrtInvS_mat(:,:,l), -0.5d0)
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtInvS_mat(i,:,l) = 0.d0
             self%sqrtInvS_mat(:,i,l) = 0.d0
          end if
       end do

    end do

!!$    call mpi_finalize(i)
!!$    stop

  end subroutine updateS

  subroutine read_binfile(self, binfile)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)  :: binfile

    integer(i4b)       :: unit, l1, l2, n
    character(len=512) :: line
    character(len=1), dimension(6) :: stat

    unit = getlun()

    ! Find number of bins
    n = 0
    open(unit,file=trim(binfile))
    do while (.true.)
       read(unit,'(a)',end=1) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       read(line,*) l1, l2
       if (l1 <= self%lmax .and. l2 <= self%lmax .and. l1 >= 0 .and. l2 >= 0) n = n+1
    end do
1   close(unit)

    if (n == 0) call report_error("Error: Cl bin file " // trim(binfile) &
         & // " does not contain any valid entries")

    self%nbin = n
    allocate(self%stat(n,6), self%bins(n,2))
    open(unit,file=trim(binfile))
    n = 0
    do while (.true.)
       read(unit,'(a)',end=2) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       read(line,*) l1, l2
       if (l1 <= self%lmax .and. l2 <= self%lmax .and. l1 >= 0 .and. l2 >= 0) then
          n = n+1
          read(line,*) l1, l2, stat   !self%bins(n,1), self%bins(n,2), self%stat(n,:)
          self%bins(n,1) = l1
          self%bins(n,2) = l2
          self%stat(n,:) = stat
       end if
    end do
2   close(unit)

  end subroutine read_binfile

  subroutine matmulS(self, map, alm, info)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    integer(i4b) :: i, l
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if

    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .true.)
          map%alm(i,:) = f_apod**2 * matmul(self%S_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .true.)
          if (l <= self%lmax) then
             alm(i,:) = f_apod**2 * matmul(self%S_mat(:,:,l), alm(i,:))
          else
             alm(i,:) = 0.d0
          end if
       end do
    end if
  end subroutine matmulS

  subroutine matmulSqrtS(self, map, alm, info, diag)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    logical(lgt),                      intent(in),    optional :: diag
    integer(i4b) :: i, l, j
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if
    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .true.)
          if (present(diag)) then
             do j = 1, self%info%nmaps
                map%alm(i,j) = f_apod * sqrt(self%S_mat(j,j,l)) * map%alm(i,j)
             end do
          else
             map%alm(i,:) = f_apod * matmul(self%sqrtS_mat(:,:,l), map%alm(i,:))
          end if
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .true.)
          if (l <= self%lmax) then
             if (present(diag)) then
                do j = 1, self%info%nmaps
                   alm(i,j) = f_apod * sqrt(self%S_mat(j,j,l)) * alm(i,j)
                end do
             else
                alm(i,:) = f_apod * matmul(self%sqrtS_mat(:,:,l), alm(i,:))
             end if
          else
             alm(i,:) = 0.d0
          end if
          !if (self%info%myid == 0 .and. i == 4) write(*,*) l, real(f_apod,sp), real(self%sqrtS_mat(:,:,l),sp),  real(alm(i,:),sp), self%Dl(l,:), trim(self%label)
       end do
    end if
  end subroutine matmulSqrtS

  subroutine matmulSqrtInvS(self, map, alm, info)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    integer(i4b) :: i, l
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if
    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .false.)
          map%alm(i,:) = f_apod * matmul(self%sqrtInvS_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, .false.)
          if (l <= self%lmax) then
             alm(i,:) = f_apod * matmul(self%sqrtInvS_mat(:,:,l), alm(i,:))
          else
             alm(i,:) = 0.d0
          end if
       end do
    end if
  end subroutine matmulSqrtInvS

  function get_Cl_apod(l, l_apod, lmax, positive)
    implicit none
    integer(i4b), intent(in) :: l, l_apod, lmax
    logical(lgt), intent(in) :: positive
    real(dp)                 :: get_Cl_apod
    real(dp), parameter :: alpha = log(1d3)
    if (l_apod > 0) then
       if (l <= l_apod) then
          get_Cl_apod = 1.d0
       else if (l > lmax) then
          get_Cl_apod = 0.d0
       else
          get_Cl_apod = exp(-alpha * (l-l_apod)**2 / real(lmax-l_apod+1,dp)**2)
       end if
    else
       if (l >= abs(l_apod)) then
          get_Cl_apod = 1.d0
       else if (l == 0 .or. l > lmax) then
          get_Cl_apod = 0.d0
       else
          get_Cl_apod = exp(-alpha * (abs(l_apod)-l)**2 / real(abs(l_apod)-1,dp)**2)
       end if
    end if
    if (.not. positive .and. get_Cl_apod /= 0.d0) get_Cl_apod = 1.d0 / get_Cl_apod
  end function get_Cl_apod

  subroutine read_Cl_file(self, clfile)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)    :: clfile

    integer(i4b)       :: unit, l, n
    character(len=512) :: line
    real(dp)           :: Dl_TT, Dl_TE, Dl_EE, Dl_BB
    real(dp),          allocatable, dimension(:,:) :: clin
    character(len=80),              dimension(120) :: header
    n    = len(trim(clfile))
    unit = getlun()

    self%Dl = 0.d0
    if (clfile(n-3:n) == 'fits') then
       ! Assume input file is in standard Healpix FITS format
       allocate(clin(0:self%lmax,4))
       call fits2cl(clfile, clin, self%lmax, 4, header)
       do l = 0, self%lmax
          self%Dl(l,TT) = clin(l,1) * (l*(l+1)/(2.d0*pi))
          if (self%nmaps == 3) then
             self%Dl(l,EE) = clin(l,2)
             self%Dl(l,BB) = clin(l,3)
             self%Dl(l,TE) = clin(l,4)
          end if
       end do
       deallocate(clin)
    else
       ! Assume input file is in ASCII format with {l, TT, TE, EE, BB} ordering, and
       ! normalization is C_l * l(l+1)/2pi
       open(unit,file=trim(clfile))
       do while (.true.)
          read(unit,'(a)',end=3) line
          line = trim(line)
          if (line(1:1) == '#' .or. trim(line) == '') cycle
          if (self%nmaps == 1) then
             read(line,*) l, Dl_TT
             if (l <= self%lmax) then
                self%Dl(l,TT) = Dl_TT 
             end if
          else
             !write(*,*) trim(line)
             read(line,*) l, Dl_TT, Dl_TE, Dl_EE, Dl_BB
             if (l <= self%lmax) then
                self%Dl(l,TT) = Dl_TT 
                self%Dl(l,TE) = Dl_TE 
                self%Dl(l,EE) = Dl_EE 
                self%Dl(l,BB) = Dl_BB 
             end if
          end if
       end do
3      close(unit)
    end if

    if (self%Dl(0,TT) == 0.d0) then
       if (self%info%myid == 0) write(*,*) 'Warning: Input spectrum file ' // trim(clfile) // ' has vanishing monopole'
    end if
    
  end subroutine read_Cl_file

  subroutine binCls(self)
    implicit none
    class(comm_Cl), intent(inout) :: self

    real(dp)     :: val, n
    integer(i4b) :: j, k, l

    do k = 1, self%nbin
       do j = 1, self%nspec
          if (self%stat(k,j) == '0') then
             self%Dl(self%bins(k,1):self%bins(k,2),j) = 0.d0
             cycle
          end if
          val = 0.d0
          n   = 0.d0
          do l = self%bins(k,1), self%bins(k,2)
             val = val + (2*l+1) * self%Dl(l,j)
             n   = n   + 2*l+1
          end do
          self%Dl(self%bins(k,1):self%bins(k,2),j) = val/n
       end do
    end do
    
  end subroutine binCls

  subroutine sampleCls(self, map, handle)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    type(planck_rng), intent(inout) :: handle

    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       call sample_Cls_inverse_wishart(self, map, handle)
    case ('power_law')
       call sample_Cls_powlaw(self, map)
    case ('exp')
       call sample_Cls_powlaw(self, map)
    end select

    call self%updateS
    
  end subroutine sampleCls

  subroutine sample_Cls_inverse_wishart(self, map, handle)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    class(comm_map),  intent(in)    :: map
    type(planck_rng), intent(inout) :: handle

    integer(i4b) :: bin, b, i, j, k, l, m, n, p, ind, b1, b2, col, ierr
    logical(lgt), allocatable, dimension(:,:) :: pattern
    integer(i4b), allocatable, dimension(:) :: i2p
    real(dp), allocatable, dimension(:,:)   :: y, y_t
    real(dp), allocatable, dimension(:,:)   :: C_b
    real(dp), allocatable, dimension(:,:)   :: sigma, s, sigma_l

    allocate(sigma_l(0:self%lmax,self%nspec))

    call map%getSigmaL(sigma_l)

    if (self%info%myid == 0) then
       allocate(sigma(self%nmaps,self%nmaps), pattern(self%nmaps,self%nmaps))

       do bin = 1, self%nbin
          if (.not. any(self%stat(bin,:) == 'S') .and. .not. any(self%stat(bin,:) == 'M')) cycle
          
          sigma = 0.d0
          n     = 0.d0
          do l = self%bins(bin,1), self%bins(bin,2)
             n   = n + 2*l+1
             k = 1
             do i = 1, self%nmaps
                do j = i, self%nmaps
                   sigma(i,j) = sigma(i,j) + (2*l+1) * sigma_l(l,k) * l*(l+1)/(2.d0*pi)
                   k          = k+1
                end do
             end do
          end do

          do i = 1, self%nmaps
             do j = i, self%nmaps
                sigma(j,i) = sigma(i,j)
             end do
          end do

          ! Set up sampling pattern
          k = 1
          do i = 1, self%nmaps
             do j = i, self%nmaps
                pattern(i,j) = self%stat(bin,k) == 'S' .or. self%stat(bin,k) == 'M'
                pattern(j,i) = pattern(i,j)
                k            = k+1
             end do
          end do

          ! Keep sampling until all elements have been treated
          do while (any(pattern))
             ! Find which elements to include this time
             do col = 1, self%nmaps
                if (any(pattern(:,col))) exit
             end do

             ! Extract the appropriate segment
             p = count(pattern(:,col))
             allocate(s(p,p), y(p,1), y_t(1,p), i2p(p), C_b(p,p))
             j = 1
             do i = 1, self%nmaps
                if (pattern(i,col)) then
                   i2p(j) = i
                   j      = j+1
                end if
             end do
             s = sigma(i2p,i2p)
             call invert_matrix(s)
             call cholesky_decompose_single(s)
             
             ! Draw sample
             C_b = 0.d0
             do i = 1, n - p - 1
                do j = 1, p
                   y(j,1) = rand_gauss(handle)
                end do
                y_t(1,:) = matmul(s, y(:,1))
                y(:,1)   = y_t(1,:)
                C_b      = C_b + matmul(y(:,1:1), y_t(1:1,:))
             end do
             call invert_matrix(C_b)

             ! Copy information over to output Dl array
             do i = 1, p
                do j = i, p
                   ind = i2p(i)*(1-i2p(i))/2 + (i2p(i)-1)*self%nmaps + i2p(j)
                   self%Dl(self%bins(bin,1):self%bins(bin,2),ind) = C_b(i,j) 
                   !write(*,*) bin, i, j, C_b(i,j), sigma_l(self%bins(bin,1),ind) * self%bins(bin,1)*(self%bins(bin,1)+1)/2/pi
                end do
             end do
          
             ! Remove current elements from pattern, and prepare for next round
             pattern(i2p,i2p) = .false.
             deallocate(s, y, y_t, i2p, C_b)
             
          end do
       end do
       deallocate(sigma, pattern)
    end if

    call mpi_bcast(self%Dl, size(self%Dl), MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)

    deallocate(sigma_l)

  end subroutine sample_Cls_inverse_wishart

  subroutine sample_Cls_powlaw(self, map)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    
    
  end subroutine sample_Cls_powlaw

  subroutine writeFITS(self, chain, iter, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    integer(i4b),     intent(in)    :: chain, iter
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath

    character(len=4) :: ctext
    character(len=6) :: itext

    if (trim(self%type) == 'none') return
    if (self%info%myid /= 0) return

    call int2string(chain, ctext)
    call int2string(iter,  itext)
    
    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext, hdffile=hdffile, hdfpath=hdfpath)
    case ('power_law')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    case ('exp')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    end select

  end subroutine writeFITS

  subroutine write_Dl_to_FITS(self, postfix, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(in) :: self
    character(len=*), intent(in) :: postfix
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath
    
    integer(i4b) :: i, l, unit
    character(len=512) :: filename
    real(dp), allocatable, dimension(:,:) :: sigmal

    if (trim(self%outdir) == 'none') return

    unit = getlun()
    filename = trim(self%outdir) // '/cls_' // trim(self%label) // '_' // trim(postfix) // '.dat'
    open(unit,file=trim(filename), recl=1024)
    if (self%nspec == 1) then
       write(unit,*) '# Columns are {l, Dl_TT}'
    else
       write(unit,*) '# Columns are {l, Dl_TT, Dl_TE, Dl_TB, Dl_EE, Dl_EB, Dl_BB}'
    end if
    do l = 0, self%lmax
       if (self%nspec == 1) then
          write(unit,fmt='(i6,e16.8)') l, self%Dl(l,:)
       else
          write(unit,fmt='(i6,6e16.8)') l, self%Dl(l,:)
       end if
    end do
    close(unit)

    if (present(hdffile)) call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl', self%Dl)       
    
  end subroutine write_Dl_to_FITS

  subroutine write_powlaw_to_FITS(self, postfix, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)    :: postfix
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath

    integer(i4b) :: unit
    character(len=512) :: filename

    if (trim(self%outdir) == 'none') return

    unit    = getlun()
    filename = trim(self%outdir) // '/cls_' // trim(self%label) // '_powlaw_' // &
         & trim(postfix) // '.dat' 

    if (self%iter == 0) then
       open(unit,file=trim(filename), recl=1024)
       if (self%nspec == 1) then
          write(unit,*) '# Columns are {iteration, A_TT, beta_TT}'
       else
          write(unit,*) '# Columns are {iteration, A_TT, A_EE, A_BB, beta_TT, beta_EE, beta_BB}'
       end if
       close(unit)
    end if

    self%iter = self%iter + 1
    open(unit,file=trim(filename), recl=1024, position='append')
    if (self%nspec == 1) then
       write(unit,fmt='(i8,2e16.8)') self%iter, self%amp, self%beta
    else
       write(unit,fmt='(i8,6e16.8)') self%iter, self%amp, self%beta
    end if
    close(unit)

    if (present(hdffile)) then
       call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
    end if
    
  end subroutine write_powlaw_to_FITS

  subroutine write_sigma_l(filename, sigma_l)
    implicit none
    character(len=*),                   intent(in) :: filename
    real(dp),         dimension(0:,1:), intent(in) :: sigma_l

    integer(i4b) :: unit, l, lmax, nspec

    lmax  = size(sigma_l,1)-1
    nspec = size(sigma_l,2)

    unit = getlun()
    open(unit,file=trim(filename), recl=1024)
    if (nspec == 1) then
       write(unit,*) '# Columns are {l, Dl_TT}'
    else
       write(unit,*) '# Columns are {l, Dl_TT, Dl_TE, Dl_TB, Dl_EE, Dl_EB, Dl_BB}'
    end if
    do l = 0, lmax
       if (nspec == 1) then
          write(unit,fmt='(i6,e16.8)') l, sigma_l(l,:) * l*(l+1)/(2.d0*pi)
       else
          write(unit,fmt='(i6,6e16.8)') l, sigma_l(l,:) * l*(l+1)/(2.d0*pi)
       end if
    end do
    close(unit)

  end subroutine write_sigma_l

  function getCl(self, l, p)
    implicit none
    class(comm_Cl),  intent(in) :: self    
    integer(i4b),    intent(in) :: l, p
    real(dp)                    :: getCl

    integer(i4b) :: j

    j = p*(1-p)/2 + (p-1)*self%nmaps + p
    if (l == 0) then
       getCl = self%Dl(l,j)
    else
       getCl = self%Dl(l,j) / (l*(l+1)/2.d0/pi)
    end if
    getCl = getCl * get_Cl_apod(l, self%l_apod, self%lmax, .true.)**2

  end function getCl
  
end module comm_Cl_mod
