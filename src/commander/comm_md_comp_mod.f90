module comm_md_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_line_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_md_comp, initialize_md_comps

  !**************************************************
  !           Monopole/dipole component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_md_comp
     !integer(i4b)                            :: ref_band
   contains
     procedure :: S    => evalSED
  end type comm_md_comp

  interface comm_md_comp
     procedure constructor
  end interface comm_md_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs, band, label, mu, rms, def)
    implicit none
    type(comm_params),              intent(in) :: cpar
    integer(i4b),                   intent(in) :: id, id_abs, band
    character(len=*),               intent(in) :: label
    real(dp),         dimension(4), intent(in) :: mu, def
    real(dp),         dimension(2), intent(in) :: rms
    class(comm_md_comp), pointer               :: constructor

    integer(i4b) :: i, j, l, m
    type(comm_mapinfo), pointer :: info

!    write(*,*) 'mu', trim(label), real(mu,sp)
!    write(*,*) 'rms', trim(label), real(rms,sp)

    ! General parameters
    allocate(constructor)

    ! Initialize comm_comp_mod parameters
    constructor%id              = id
    constructor%active          = cpar%cs_include(id_abs)
    constructor%label           = data(band)%label
    constructor%type            = cpar%cs_type(id_abs)
    constructor%class           = cpar%cs_class(id_abs)
    constructor%cg_samp_group   = cpar%cs_cg_samp_group(id_abs)
    constructor%unit            = data(band)%unit
    constructor%nu_ref          = data(band)%bp%nu_c
    constructor%cg_scale        = 1.d0
    constructor%myid            = cpar%myid_chain
    constructor%comm            = cpar%comm_chain
    constructor%numprocs        = cpar%numprocs_chain
    !constructor%ref_band = band

    ! Set up conversion factor between RJ and native component unit
    select case (trim(constructor%unit))
    case ('uK_cmb')
       constructor%RJ2unit_ = data(band)%bp%a2t
    case ('MJy/sr') 
       constructor%RJ2unit_ = data(band)%bp%a2t / data(band)%bp%f2t
    case ('uK_RJ') 
       constructor%RJ2unit_ = 1.d0
    case default
       call report_error('Unsupported unit: ' // trim(constructor%unit))
    end select

    ! Initialize comm_diffuse_comp_mod parameters
    constructor%pol      = .false.
    constructor%nside    = data(band)%map%info%nside
    constructor%lmax_amp = 1
    constructor%l_apod   = 2
    constructor%lmax_ind = 0
    constructor%cltype   = 'binned'
    constructor%nmaps    = 1
    !info          => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_amp, &
    !     & constructor%nmaps, constructor%pol)
    info          => comm_mapinfo(cpar%comm_chain, 128, constructor%lmax_amp, &
         & constructor%nmaps, constructor%pol)

    ! Diffuse preconditioner variables
    call add_to_npre(1,constructor%nside,1,1)

    ! Initialize amplitude and prior maps
    constructor%x   => comm_map(info)
    constructor%mu  => comm_map(info)
    constructor%ncr =  size(constructor%x%alm)
    do i = 0, constructor%x%info%nalm-1
       call constructor%x%info%i2lm(i,l,m)
       if (l == 0) then ! Monopole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi) * def(1) / constructor%RJ2unit_
          constructor%mu%alm(i,1) = sqrt(4.d0*pi) * mu(1)  / constructor%RJ2unit_
       end if
       if (l == 1 .and. m == -1) then ! Y dipole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(3) / constructor%RJ2unit_
          constructor%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(3)  / constructor%RJ2unit_
       end if
       if (l == 1 .and. m ==  0) then ! Z dipole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(4) / constructor%RJ2unit_
          constructor%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(4)  / constructor%RJ2unit_
       end if
       if (l == 1 .and. m ==  1) then ! X dipole
          constructor%x%alm(i,1)  = -sqrt(4.d0*pi/3.d0) * def(2) / constructor%RJ2unit_
          constructor%mu%alm(i,1) = -sqrt(4.d0*pi/3.d0) * mu(2)  / constructor%RJ2unit_
       end if
    end do


    ! Allocate mixing matrix
    allocate(constructor%F(numband), constructor%F_mean(numband,constructor%nmaps), constructor%F_null(numband))
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & constructor%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
       if (i == band) then
          constructor%F(i)%p      => comm_map(info)
          constructor%F_null(i)   = .false.
          constructor%F(i)%p%map  = constructor%RJ2unit_
          constructor%F(i)%p%alm  = constructor%RJ2unit_ * sqrt(4.d0*pi)
          constructor%F_mean(i,:) = constructor%RJ2unit_
       else
          constructor%F_null(i)   = .true.
          !constructor%F(i)%p%map  = 0.d0
          !constructor%F(i)%p%alm  = 0.d0
          constructor%F_mean(i,:) = 0.d0
       end if
    end do

    ! Initialize output beam
    constructor%B_out => comm_B_bl(cpar, constructor%x%info, 0, 0, fwhm=0.d0, init_realspace=.false.)

    ! Initialize power spectrum
    allocate(constructor%Cl)
    constructor%Cl%type   = 'binned'
    constructor%Cl%info   => constructor%x%info
    constructor%Cl%label  = 'md_'//trim(data(band)%label)
    constructor%Cl%lmax   = 1
    constructor%Cl%nmaps  = 1
    constructor%Cl%nspec  = 1
    constructor%Cl%outdir = 'none'
    allocate(constructor%Cl%Dl(0:1,1), constructor%Cl%sqrtS_mat(1,1,0:1))
    allocate(constructor%Cl%sqrtInvS_mat(1,1,0:1), constructor%Cl%S_mat(1,1,0:1))
    constructor%Cl%Dl(0,1) = 4.d0*pi      * rms(1)**2 / constructor%RJ2unit_**2
    constructor%Cl%Dl(1,1) = 4.d0*pi/3.d0 * rms(2)**2 / constructor%RJ2unit_**2
    constructor%Cl%S_mat(1,1,0:1)        = rms**2     / constructor%RJ2unit_**2
    constructor%Cl%sqrtS_mat(1,1,0:1)    = rms        / constructor%RJ2unit_
    if (rms(1) > 0.d0) constructor%Cl%sqrtInvS_mat(1,1,0) = 1.d0/rms(1) * constructor%RJ2unit_
    if (rms(2) > 0.d0) constructor%Cl%sqrtInvS_mat(1,1,1) = 1.d0/rms(2) * constructor%RJ2unit_

    ! Initialize md_mod specific parameters
    constructor%npar = 0

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    do i = 1, numband
       if (i == band) then
          constructor%F_int(i)%p => comm_F_line(constructor, data(i)%bp, .true., 1.d0, -1)
       else
          constructor%F_int(i)%p => comm_F_line(constructor, data(i)%bp, .true., 0.d0, -1)
       end if
    end do

  end function constructor

  ! Definition:
  !    SED  = delta_{band,ref_band}
  function evalSED(self, nu, band, theta)
    class(comm_md_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, ind

    evalSED = 1.d0

  end function evalSED

  function initialize_md_comps(cpar, id, id_abs, n)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_md_comp), pointer   :: initialize_md_comps

    integer(i4b)        :: i, unit
    real(dp)            :: mu(4), rms(2), def(4)
    character(len=1024) :: line, label
    class(comm_comp), pointer :: c

    unit  = getlun()

    ! Find number of lines
    n = 0
    open(unit, file=trim(cpar%cs_SED_template(1,id_abs)), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       read(line,*) label, def, mu, rms
       do i = 1, numband
          if (trim(label) == trim(data(i)%label)) exit
       end do
       if (i > numband) cycle

       if (n == 0) then
          initialize_md_comps => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
       else
          c => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
          call initialize_md_comps%add(c)
       end if
       n = n+1
    end do
1   close(unit)
  
    if (n < numband .and. cpar%myid == 0) then
       write(*,'(a,i6)') '  Warning: Number of channels without a monopole/dipole definition = ', numband-n
    end if
  
  end function initialize_md_comps

  
end module comm_md_comp_mod
