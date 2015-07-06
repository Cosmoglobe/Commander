module comm_template_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  implicit none

  private
  public comm_template_comp

  type Tnu
     real(dp), allocatable, dimension(:,:) :: map
  end type Tnu
  
  !**************************************************
  !            Template class
  !**************************************************
  type, extends (comm_comp) :: comm_template_comp
     real(dp)                                :: x     ! Amplitude
     real(dp),     allocatable, dimension(:) :: theta ! Spectral parameters
     integer(i4b), allocatable, dimension(:) :: ind   ! Frequency index per channel
     type(Tnu),    allocatable, dimension(:) :: T     ! Spatial template, one per channel
     real(dp),     allocatable, dimension(:) :: f     ! Frequency scaling, one per channel
   contains
     procedure :: initTemplate
     procedure :: S => evalSED
     procedure :: dumpFITS => dumpTemplateToFITS
     procedure :: getBand  => evalTemplateBand
     procedure :: projectBand  => projectTemplateBand
  end type comm_template_comp

contains

  subroutine initTemplate(self, cpar, id)
    implicit none
    class(comm_template_comp)             :: self
    type(comm_params),         intent(in) :: cpar
    integer(i4b),              intent(in) :: id

    call self%initComp(cpar, id)

    ! Initialize variables specific to template source type
    
  end subroutine initTemplate

  ! Currently returning SED in data units, not brightness temperature
  function evalSED(self, nu, band, theta)
    class(comm_template_comp), intent(in)           :: self
    real(dp),                  intent(in), optional :: nu
    integer(i4b),              intent(in), optional :: band
    real(dp), dimension(1:),   intent(in), optional :: theta
    real(dp)                                        :: evalSED

    if (self%ind(band) > 0) then
       evalSED = self%f(self%ind(band))
    else
       evalSED = 0.d0
    end if
    
  end function evalSED

  function evalTemplateBand(self, band, amp_in, pix)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    real(dp),        dimension(:,:), allocatable                        :: evalTemplateBand

  end function evalTemplateBand
  
  ! Return component projected from map
  function projectTemplateBand(self, band, map)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    real(dp),        dimension(:,:), allocatable                        :: projectTemplateBand
  end function projectTemplateBand
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpTemplateToFITS(self, postfix, dir)
    implicit none
    character(len=*),                        intent(in)           :: postfix
    class(comm_template_comp),               intent(in)           :: self
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i
    character(len=512) :: filename
    
  end subroutine dumpTemplateToFITS
  
end module comm_template_comp_mod
