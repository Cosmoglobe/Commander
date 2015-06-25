module comm_diffuse_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  implicit none

  private
  public comm_diffuse_comp
  
  !**************************************************
  !            Diffuse component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nmaps, nx, x0
     logical(lgt)       :: pol
     integer(i4b)       :: lmax_amp, lmax_ind, lpiv
     real(dp), allocatable, dimension(:,:) :: cls

     class(comm_map),               pointer :: mask
     class(comm_map),               pointer :: x      ! Spatial parameters
     class(comm_map), dimension(:), pointer :: theta  ! Spectral parameters
   contains
     procedure :: initDiffuse
!!$     procedure :: a        => evalDiffuseAmp
!!$     procedure :: F        => evalDiffuseMixmat
!!$     procedure :: sim      => simDiffuseComp
!!$     procedure :: dumpHDF  => dumpDiffuseToHDF
!!$     procedure :: dumpFITS => dumpDiffuseToFITS
  end type comm_diffuse_comp

contains

  subroutine initDiffuse(self, cpar, id)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id

    type(comm_mapinfo), pointer :: info
    
    call self%initComp(cpar, id)

    ! Initialize variables specific to diffuse source type
    self%pol      = cpar%cs_polarization(id)
    self%nside    = cpar%cs_nside(id)
    self%lmax_amp = cpar%cs_lmax_amp(id)
    self%lmax_ind = cpar%cs_lmax_ind(id)
    self%nmaps    = 1; if (self%pol) self%nmaps = 3
    info          => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id)) == 'zero' .or. trim(cpar%cs_input_amp(id)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       self%x => comm_map(info, cpar%cs_input_amp(id))
       call self%x%YtW(cleanup=.true.)
    end if
    self%ncr = size(self%x%alm)
       

  end subroutine initDiffuse

!!$  ! Evaluate amplitude map in brightness temperature at reference frequency
!!$  function evalDiffuseAmp(self, nside, nmaps, pix, x_1D, x_2D)
!!$    class(comm_diffuse_comp),                  intent(in)           :: self
!!$    integer(i4b),                              intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:),   intent(in), optional :: pix
!!$    real(dp),                  dimension(:),   intent(in), optional :: x_1D
!!$    real(dp),                  dimension(:,:), intent(in), optional :: x_2D
!!$    real(dp),     allocatable, dimension(:,:)                       :: evalDiffuseAmp
!!$
!!$    if (present(x_1D)) then
!!$       ! Input array is complete packed data vector
!!$       
!!$    else if (present(x_2D)) then
!!$       ! Input array is alms(:,:)
!!$       
!!$    else
!!$       ! Use internal array
!!$       
!!$    end if
!!$    
!!$  end function evalDiffuseAmp
!!$
!!$  ! Evaluate amplitude map in brightness temperature at reference frequency
!!$  function evalDiffuseMixmat(self, nside, nmaps, pix)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    integer(i4b),                            intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:), intent(in), optional :: pix
!!$    real(dp),     allocatable, dimension(:,:)                     :: evalDiffuseMixmat
!!$  end function evalDiffuseMixmat
!!$
!!$  ! Generate simulated component
!!$  function simDiffuseComp(self, handle, nside, nmaps, pix)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    type(planck_rng),                        intent(inout)        :: handle
!!$    integer(i4b),                            intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:), intent(in), optional :: pix
!!$    real(dp),     allocatable, dimension(:,:)                     :: simDiffuseComp
!!$  end function simDiffuseComp
!!$  
!!$  ! Dump current sample to HDF chain files
!!$  subroutine dumpDiffuseToHDF(self, filename)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    character(len=*),                        intent(in)           :: filename
!!$  end subroutine dumpDiffuseToHDF
!!$  
!!$  ! Dump current sample to HEALPix FITS file
!!$  subroutine dumpDiffuseToFITS(self, dir)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    character(len=*),                        intent(in)           :: dir
!!$  end subroutine dumpDiffuseToFITS
  
end module comm_diffuse_comp_mod
