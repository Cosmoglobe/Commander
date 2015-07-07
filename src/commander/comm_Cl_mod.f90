module comm_Cl_mod
  use comm_param_mod
  use comm_map_mod
  use math_tools
  implicit none

  private
  public comm_Cl

  integer(i4b), parameter :: TT = 1
  integer(i4b), parameter :: TE = 2
  integer(i4b), parameter :: TB = 3
  integer(i4b), parameter :: EE = 4
  integer(i4b), parameter :: EB = 5
  integer(i4b), parameter :: BB = 6
  
  type :: comm_Cl
     ! General parameters
     class(comm_mapinfo), pointer :: info
     character(len=512)           :: type  ! {none, binned, power_law}
     character(len=512)           :: label ! {none, binned, power_law}
     character(len=512)           :: outdir
     integer(i4b)                 :: lmax, nmaps, nspec
     integer(i4b)                 :: poltype  ! {1 = {T+E+B}, 2 = {T,E+B}, 3 = {T,E,B}}
     real(dp),         allocatable, dimension(:,:)   :: Dl
     real(dp),         allocatable, dimension(:,:,:) :: sqrtS_mat, S_mat

     ! Bin parameters
     integer(i4b) :: nbin
     character(len=1), allocatable, dimension(:,:) :: stat
     integer(i4b), allocatable, dimension(:,:)     :: bins
     
     ! Power law parameters
     character(len=512) :: plfile
     integer(i4b)       :: iter = 0
     real(dp)           :: prior(2), lpiv
     real(dp), allocatable, dimension(:) :: amp, beta
   contains
     ! Data procedures
     procedure :: S        => matmulS
     procedure :: sqrtS    => matmulSqrtS
     procedure :: sampleCls
     procedure :: read_binfile
     procedure :: read_Cl_file
     procedure :: binCls
     procedure :: writeFITS
     procedure :: updatePowlaw
     procedure :: updateS
  end type comm_Cl

  interface comm_Cl
     procedure constructor
  end interface comm_Cl
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    type(comm_mapinfo), target,         intent(in) :: info
    integer(i4b),                       intent(in) :: id
    class(comm_Cl),                     pointer    :: constructor

    integer(i4b)       :: l, nmaps
    character(len=512) :: datadir, binfile
    logical(lgt)       :: pol

    allocate(constructor)
    
    ! General parameters
    constructor%type   = cpar%cs_cltype(id)
    if (trim(constructor%type) == 'none') return
    
    constructor%info   => info
    constructor%label  = cpar%cs_label(id)
    constructor%lmax   = cpar%cs_lmax_amp(id)
    constructor%nmaps  = 1; if (cpar%cs_polarization(id)) constructor%nmaps = 3
    constructor%nspec  = 1; if (cpar%cs_polarization(id)) constructor%nspec = 6
    constructor%outdir = cpar%outdir
    datadir            = cpar%datadir
    nmaps              = constructor%nmaps

    allocate(constructor%Dl(0:constructor%lmax,constructor%nspec))
    allocate(constructor%sqrtS_mat(nmaps,nmaps,0:constructor%lmax))
    allocate(constructor%S_mat(nmaps,nmaps,0:constructor%lmax))

    if (trim(constructor%type) == 'binned') then
       call constructor%read_binfile(trim(datadir) // '/' // trim(cpar%cs_binfile(id)))
       call constructor%read_Cl_file(trim(datadir) // '/' // trim(cpar%cs_clfile(id)))
       call constructor%binCls
    else if (trim(constructor%type) == 'power_law') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id)
       constructor%prior         = cpar%cs_cl_prior(id,:)
       constructor%poltype       = cpar%cs_cl_poltype(id)
       call constructor%updatePowlaw(cpar%cs_cl_amp_def(id,1:nmaps), cpar%cs_cl_beta_def(id,1:nmaps))
    else
       call report_error("Unknown Cl type: " // trim(constructor%type))
    end if

    call constructor%updateS
    
  end function constructor


  subroutine updatePowlaw(self, amp, beta)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in)    :: amp, beta

    integer(i4b) :: i, j, l
    
    self%amp   = amp
    self%beta  = beta
    self%Dl    = 0.d0
    do i = 1, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = 1, self%lmax
          self%Dl(l,j) = amp(i) * (l/self%lpiv)**beta(i)
       end do
       self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updatePowlaw

  subroutine updateS(self)
    implicit none
    class(comm_Cl),                intent(inout) :: self

    integer(i4b) :: i, j, k, l
    logical(lgt) :: ok(self%nmaps)

    self%sqrtS_mat = 0.d0
    self%S_mat     = 0.d0
    do l = 1, self%lmax
       do i = 1, self%nmaps
          do j = i, self%nmaps
             k = i*(1-i)/2 + (i-1)*self%nmaps + j
             self%sqrtS_mat(i,j,l) = self%Dl(l,k) / (l*(l+1)/(2.d0*pi))
             if (i == j) ok(i) = self%Dl(l,k) > 0.d0
             if (.not. ok(i)) self%sqrtS_mat(i,j,l) = 1.d0
          end do
       end do
       call compute_hermitian_root(self%sqrtS_mat(:,:,l), 0.5d0)
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtS_mat(i,:,l) = 0.d0
             self%sqrtS_mat(:,i,l) = 0.d0
          end if
       end do
       self%S_mat(:,:,l) = matmul(self%sqrtS_mat(:,:,l), self%sqrtS_mat(:,:,l))
       
    end do

  end subroutine updateS

  subroutine read_binfile(self, binfile)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)  :: binfile

    integer(i4b)       :: unit, l1, l2, n
    character(len=512) :: line

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
          read(line,*) self%bins(n,1), self%bins(n,2), self%stat(n,:)
       end if
    end do
2   close(unit)

  end subroutine read_binfile

  subroutine matmulS(self, map, alm)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    integer(i4b) :: i, l
    if (trim(self%type) == 'none') return

    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          map%alm(i,:) = matmul(self%S_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          alm(i,:) = matmul(self%S_mat(:,:,l), alm(i,:))
       end do
    end if
  end subroutine matmulS

  subroutine matmulSqrtS(self, map, alm)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    integer(i4b) :: i, l
    if (trim(self%type) == 'none') return
    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          map%alm(i,:) = matmul(self%sqrtS_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          alm(i,:) = matmul(self%sqrtS_mat(:,:,l), alm(i,:))
       end do
    end if
  end subroutine matmulSqrtS

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
       ! Assume input file is in ASCII format with {l, TT, EE, BB, TE} ordering, and
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
             read(line,*) l, Dl_TT, Dl_EE, Dl_BB, Dl_TE
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

  end subroutine read_Cl_file

  subroutine binCls(self)
    implicit none
    class(comm_Cl), intent(inout) :: self

    real(dp)     :: val, n
    integer(i4b) :: j, k, l

    do k = 1, self%nbin
       do j = 1, self%nspec
          if (self%stat(k,j) == '0' .or. self%stat(k,j) == 'C') cycle
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

  subroutine sampleCls(self, map)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map

    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       call sample_Cls_inverse_wishart(self, map)
    case ('power_law')
       call sample_Cls_powlaw(self, map)
    end select
    
  end subroutine sampleCls

  subroutine sample_Cls_inverse_wishart(self, map)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map


  end subroutine sample_Cls_inverse_wishart

  subroutine sample_Cls_powlaw(self, map)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    
    
  end subroutine sample_Cls_powlaw

  subroutine writeFITS(self, chain, iter)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    integer(i4b),     intent(in)    :: chain, iter

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
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
    case ('power_law')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext)
    end select

  end subroutine writeFITS

  subroutine write_Dl_to_FITS(self, postfix)
    implicit none
    class(comm_Cl),   intent(in) :: self
    character(len=*), intent(in) :: postfix

    integer(i4b) :: i, l, unit
    character(len=512) :: filename

    unit = getlun()
    filename = trim(self%outdir) // '/cls_' // trim(self%label) // '_' // trim(postfix) // '.dat'
    open(unit,file=trim(filename), recl=1024)
    if (self%nspec == 1) then
       write(unit,*) '# Columns are {l, Dl_TT}'
    else
       write(unit,*) '# Columns are {l, Dl_TT, Dl_EE, Dl_BB, Dl_TE, Dl_TB, Dl_EB}'
    end if
    do l = 0, self%lmax
       if (self%nspec == 1) then
          write(unit,fmt='(i6,f16.3)') l, self%Dl(l,:)
       else
          write(unit,fmt='(i6,6f16.3)') l, self%Dl(l,:)
       end if
    end do
    close(unit)

  end subroutine write_Dl_to_FITS

  subroutine write_powlaw_to_FITS(self, postfix)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)    :: postfix

    integer(i4b) :: unit
    character(len=512) :: filename

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
       write(unit,fmt='(i8,2f16.3)') self%iter, self%amp, self%beta
    else
       write(unit,fmt='(i8,6f16.3)') self%iter, self%amp, self%beta
    end if
    close(unit)
    
  end subroutine write_powlaw_to_FITS
  
  
end module comm_Cl_mod
