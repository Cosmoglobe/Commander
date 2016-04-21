module comm_template_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  use comm_hdf_mod
  implicit none

  private
  public comm_template_comp

  type Tnu
     logical(lgt) :: fullsky
     integer(i4b) :: nside, np, nmaps, band
     integer(i4b), allocatable, dimension(:,:) :: pix   ! Pixel list, both in map and absolute order
     real(dp),     allocatable, dimension(:,:) :: map   ! (0:np-1,nmaps)
     real(dp),     allocatable, dimension(:)   :: f     ! Frequency scaling (nmaps)
  end type Tnu
  
  !**************************************************
  !            Template class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_template_comp
     character(len=512) :: outprefix
     real(dp),     allocatable, dimension(:)   :: x     ! Amplitude (nmaps)
     real(dp),     allocatable, dimension(:,:) :: theta ! Spectral parameters (npar,nmaps)
     type(Tnu),    allocatable, dimension(:)   :: T     ! Spatial template (nband)
   contains
     procedure :: dumpFITS => dumpTemplateToFITS
     procedure :: getBand  => evalTemplateBand
     procedure :: projectBand  => projectTemplateBand
  end type comm_template_comp

contains

  function evalTemplateBand(self, band, amp_in, pix, alm_out)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    real(dp),        dimension(:,:), allocatable                        :: evalTemplateBand

  end function evalTemplateBand
  
  ! Return component projected from map
  function projectTemplateBand(self, band, map, alm_in)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    real(dp),        dimension(:,:), allocatable                        :: projectTemplateBand
  end function projectTemplateBand
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpTemplateToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    implicit none
    class(comm_template_comp),               intent(in)           :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i
    character(len=512) :: filename
    
  end subroutine dumpTemplateToFITS
  
end module comm_template_comp_mod
