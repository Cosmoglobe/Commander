module comm_ptsrc_comp_mod
  use comm_param_mod
  use comm_template_comp_mod
  use comm_F_int_mod
  use comm_F_int_0D_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  use pix_tools
  implicit none

  private
  public comm_ptsrc_comp

  !**************************************************
  !            Template class
  !**************************************************
  type, extends (comm_template_comp) :: comm_ptsrc_comp
     character(len=512) :: id_src
     real(dp)           :: glon, glat, f_beam
     type(F_int_ptr), dimension(:), allocatable :: F_int  ! SED integrator
   contains
     procedure          :: updateF
     procedure          :: S => evalSED
  end type comm_ptsrc_comp

  interface comm_ptsrc_comp
     procedure constructor
  end interface comm_ptsrc_comp
  
contains

  function constructor(cpar, id, id_abs, glon, glat, amp, beta, id_src, P_ref)
    implicit none
    class(comm_params),       intent(in) :: cpar
    integer(i4b),             intent(in) :: id, id_abs
    real(dp),                 intent(in) :: glon, glat
    real(dp), dimension(:),   intent(in) :: amp
    real(dp), dimension(:,:), intent(in) :: beta
    character(len=*),         intent(in) :: id_src
    class(comm_ptsrc_comp),   pointer    :: P_ref
    class(comm_ptsrc_comp),   pointer    :: constructor

    integer(i4b) :: i, j, k, nlist, npix, listpix(0:10000-1), hits(10000)
    real(dp)     :: vec0(3), vec(3), r
    
    ! General parameters
    allocate(constructor)

    ! Initialize general parameters
    constructor%class     = cpar%cs_class(id_abs)
    constructor%type      = cpar%cs_type(id_abs)
    constructor%label     = cpar%cs_label(id_abs)
    constructor%nmaps     = 1; if (cpar%cs_polarization(id_abs)) constructor%nmaps = 3
    constructor%glon      = glon * DEG2RAD
    constructor%glat      = glat * DEG2RAD
    constructor%nu_ref    = cpar%cs_nu_ref(id_abs)
    constructor%outprefix = trim(cpar%cs_label(id_abs))
    constructor%id_src    = id_src
    allocate(constructor%poltype(1))
    constructor%poltype   = cpar%cs_poltype(1,id_abs)

    ! Allocate data structures
    allocate(constructor%x(constructor%nmaps))
    allocate(constructor%T(numband))

    ! Initialize beam templates and frequency scaling for each band
    do i = 1, numband
       constructor%T(i)%band    = i
       constructor%T(i)%fullsky = .false.
       constructor%T(i)%nside   = data(i)%info%nside
       constructor%T(i)%nmaps   = data(i)%info%nmaps

       ! Search for pixels controlled by current processor
       npix = 12*constructor%T(i)%nside**2
       call ang2vec(0.5d0*pi-constructor%glat, constructor%glon, vec0)
       call query_disc(constructor%T(i)%nside, vec0, data(i)%B%r_max, listpix, nlist)
       constructor%T(i)%np = 0
       k                   = 1   ! map counter
       j                   = 0   ! source counter
       do while (k <= data(i)%info%np .and. j <= nlist-1)
          if (listpix(j) == data(i)%info%pix(k)) then
             constructor%T(i)%np       = constructor%T(i)%np + 1
             hits(constructor%T(i)%np) = k
             j                         = j+1
             k                         = k+1
          else if (listpix(j) < data(i)%info%pix(k)) then
             j = j+1
          else
             k = k+1
          end if
       end do

       allocate(constructor%T(i)%pix(constructor%T(i)%np,2))
       allocate(constructor%T(i)%map(constructor%T(i)%np,constructor%T(i)%nmaps))
       constructor%T(i)%pix(:,2) = hits(1:constructor%T(i)%np)
       do j = 1, constructor%T(i)%np
          constructor%T(i)%pix(j,1) = data(i)%info%pix(constructor%T(i)%pix(j,2))
          call pix2vec_ring(constructor%T(i)%nside, constructor%T(i)%pix(j,1), vec)
          call angdist(vec0, vec, r)
          constructor%T(i)%map(j,:) = r
       end do

       ! Initialize frequency scaling structure
       allocate(constructor%T(i)%f(constructor%T(i)%nmaps))
       
    end do

    ! Initialize frequency scaling parameters
    allocate(constructor%F_int(numband))
    constructor%x = amp
    select case (trim(constructor%type))
    case ("radio")
       constructor%npar = 2   ! (alpha, beta)
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
       constructor%p_uni   = cpar%cs_p_uni(id_abs,:,:)
       constructor%p_gauss = cpar%cs_p_gauss(id_abs,:,:)

       allocate(constructor%theta(constructor%npar,constructor%nmaps))
       do i = 1, numband
          if (associated(P_ref)) then
             constructor%F_int(i)%p => P_ref%F_int(i)%p
          else
             constructor%F_int(i)%p => comm_F_int_2D(constructor, data(i)%bp)
          end if
       end do
    case ("fir")
       constructor%npar = 2   ! (beta, T_d)
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
       constructor%p_uni   = cpar%cs_p_uni(id_abs,:,:)
       constructor%p_gauss = cpar%cs_p_gauss(id_abs,:,:)

       allocate(constructor%theta(constructor%npar,constructor%nmaps))
       do i = 1, numband
          if (associated(P_ref)) then
             constructor%F_int(i)%p => P_ref%F_int(i)%p
          else
             constructor%F_int(i)%p => comm_F_int_2D(constructor, data(i)%bp)
          end if
       end do
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
       constructor%p_uni   = cpar%cs_p_uni(id_abs,:,:)
       constructor%p_gauss = cpar%cs_p_gauss(id_abs,:,:)
    case ("sz")
       constructor%npar = 0   ! (none)
       do i = 1, numband
          if (associated(P_ref)) then
             constructor%F_int(i)%p => P_ref%F_int(i)%p
          else
             constructor%F_int(i)%p => comm_F_int_0D(constructor, data(i)%bp)
          end if
       end do
    case default
       call report_error("Unknown point source model: " // trim(constructor%type))
    end select

    call constructor%updateF(beta)
    
  end function constructor

  subroutine updateF(self, beta)
    implicit none
    class(comm_ptsrc_comp),                 intent(inout) :: self
    real(dp),               dimension(:,:), intent(in)    :: beta

    integer(i4b) :: i
    
    self%theta = beta
    do i = 1, numband

       ! Temperature
       self%T(i)%f(1) = self%F_int(i)%p%eval(self%theta(:,1)) * data(i)%RJ2data()

       ! Polarization
       if (self%T(i)%nmaps == 3) then
          ! Stokes Q
          if (self%poltype(1) < 2) then
             self%T(i)%f(2) = self%T(i)%f(1)
          else
             self%T(i)%f(2) = self%F_int(i)%p%eval(self%theta(:,2)) * self%RJ2unit()
          end if
          
          ! Stokes U
          if (self%poltype(1) < 3) then
             self%T(i)%f(3) = self%T(i)%f(2)
          else
             self%T(i)%f(3) = self%F_int(i)%p%eval(self%theta(:,3)) * self%RJ2unit()
          end if
       end if
       
    end do
    
  end subroutine updateF

  function evalSED(self, nu, band, theta)
    class(comm_ptsrc_comp),    intent(in)           :: self
    real(dp),                  intent(in), optional :: nu
    integer(i4b),              intent(in), optional :: band
    real(dp), dimension(1:),   intent(in), optional :: theta
    real(dp)                                        :: evalSED

    real(dp) :: x
    
    select case (trim(self%type))
    case ("radio")
       evalSED = exp(theta(1) * (nu/self%nu_ref) + theta(2) * (log(nu/self%nu_ref))**2)
    case ("fir")
       x = h/(k_B*theta(2))
       evalSED = (exp(x*self%nu_ref)-1.d0)/(exp(x*nu)-1.d0) * (nu/self%nu_ref)**(theta(1)+1.d0)
    case ("sz")
       evalSED = 0.d0
       call report_error('SZ not implemented yet')
    end select
    
  end function evalSED
  
end module comm_ptsrc_comp_mod
