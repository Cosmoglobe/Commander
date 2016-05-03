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

    integer(i4b) :: i, j, nline, b, n, ierr
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
             constructor%p_uni(1,n)   = mu(j)-5*sigma(j)
             constructor%p_uni(2,n)   = mu(j)+5*sigma(j)
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
       call constructor%theta(i)%p%YtW
    end do


    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    j = 1
    do i = 1, numband
       if (any(constructor%ind2band == i)) then
          constructor%F_int(i)%p => comm_F_line(constructor, data(i)%bp, .true., &
               & constructor%line2RJ(j) / constructor%line2RJ_ref * data(i)%RJ2data(), j)
          j = j+1
       else
          constructor%F_int(i)%p => comm_F_line(constructor, data(i)%bp, .false., 0.d0, j)
       end if
    end do
    
    ! Initialize mixing matrix
    call constructor%updateMixmat

    deallocate(label, mu, sigma, line2RJ, poltype)

  end function constructor

  ! Definition:
  !    SED  = delta_{band,
  function evalSED(self, nu, band, theta)
    class(comm_line_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
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

  
end module comm_line_comp_mod
