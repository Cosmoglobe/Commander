module comm_B_bl_mod
  use comm_param_mod
  use comm_map_mod
  use comm_B_mod
  use comm_utils
  use spline_1D_mod
  implicit none

  private
  public comm_B_bl

  type, extends (comm_B) :: comm_B_bl
   contains
     ! Data procedures
     procedure :: conv => matmulB
  end type comm_B_bl

  interface comm_B_bl
     procedure constructor
  end interface comm_B_bl

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, fwhm, init_realspace)
    implicit none
    type(comm_params),                  intent(in)           :: cpar
    type(comm_mapinfo), target,         intent(in)           :: info
    integer(i4b),                       intent(in)           :: id, id_abs
    real(dp),                           intent(in), optional :: fwhm
    logical(lgt),                       intent(in), optional :: init_realspace
    class(comm_B_bl),   pointer                              :: constructor

    integer(i4b)       :: i
    logical(lgt)       :: init_real
    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type  =  'b_l'    
    constructor%info  => info
    if (present(fwhm)) then
       call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm)
    else
       call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, &
            & beamfile=trim(dir)//trim(cpar%ds_blfile(id_abs)), &
            & pixwin=trim(dir)//trim(cpar%ds_pixwin(id_abs)))
    end if

    ! Initialize real-space profile
    init_real = .true.; if (present(init_realspace)) init_real = init_realspace
    if (init_real) then
       call compute_radial_beam(constructor%info%nmaps, constructor%b_l, constructor%b_theta)
       constructor%r_max = maxval(constructor%b_theta(1)%x)
    end if
    
  end function constructor
  
  subroutine matmulB(self, trans, map)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map

    integer(i4b) :: i, l

    do i = 0, map%info%nalm-1
       l = map%info%lm(1,i)
       if (l <= self%info%lmax) then
          map%alm(i,:) = map%alm(i,:) * self%b_l(l,:)
       else
          map%alm(i,:) = 0.d0
       end if
    end do

  end subroutine matmulB
  
end module comm_B_bl_mod
