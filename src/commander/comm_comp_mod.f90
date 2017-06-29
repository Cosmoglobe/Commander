module comm_comp_mod
  use comm_param_mod
  use comm_bp_utils
  use comm_bp_mod
  use comm_map_mod
  use comm_cr_utils
  use comm_data_mod
  use comm_hdf_mod
  implicit none

  private
  public  :: comm_comp, ncomp, compList!, dumpCompMaps
  
  !**************************************************
  !        Generic component class definition
  !**************************************************
  type, abstract :: comm_comp
     ! Linked list variables
     class(comm_comp), pointer :: nextLink => null()
     class(comm_comp), pointer :: prevLink => null()

     ! Data variables
     logical(lgt)       :: active, init_from_HDF
     integer(i4b)       :: npar, ncr, id, nmaps, myid, comm, numprocs, cg_samp_group
     character(len=512) :: label, class, type, unit, operation
     real(dp)           :: nu_ref, RJ2unit_
     character(len=512), allocatable, dimension(:)   :: indlabel
     integer(i4b),       allocatable, dimension(:)   :: poltype
     real(dp),           allocatable, dimension(:)   :: theta_def
     real(dp),           allocatable, dimension(:,:) :: p_gauss
     real(dp),           allocatable, dimension(:,:) :: p_uni
     integer(i4b),       allocatable, dimension(:)   :: smooth_scale
     real(dp),           allocatable, dimension(:)   :: nu_min_ind
     real(dp),           allocatable, dimension(:)   :: nu_max_ind

   contains
     ! Linked list procedures
     procedure :: next    ! get the link after this link
     procedure :: prev    ! get the link before this link
     procedure :: setNext ! set the link after this link
     procedure :: add     ! add new link at the end

     ! Data procedures
     procedure                          :: initComp
     procedure(evalSED),       deferred :: S
     procedure(evalBand),      deferred :: getBand
     procedure(projectBand),   deferred :: projectBand
     procedure                          :: dumpSED
     procedure(dumpFITS),      deferred :: dumpFITS
     procedure(initHDF),       deferred :: initHDF
     procedure                          :: RJ2unit
     procedure(sampleSpecInd), deferred :: sampleSpecInd
     procedure                          :: CG_mask
     procedure(updateMixmat),  deferred :: updateMixmat
  end type comm_comp

  abstract interface
     ! Evaluate SED in brightness temperature normalized to reference frequency
     function evalSED(self, nu, band, theta)
       import comm_comp, dp, i4b
       class(comm_comp),        intent(in)           :: self
       real(dp),                intent(in), optional :: nu
       integer(i4b),            intent(in), optional :: band
       real(dp), dimension(1:), intent(in), optional :: theta
       real(dp)                                      :: evalSED
     end function evalSED

     ! Return effective signal at given frequency band
     function evalBand(self, band, amp_in, pix, alm_out)
       import i4b, dp, comm_comp, lgt
       class(comm_comp),                             intent(in)            :: self
       integer(i4b),                                 intent(in)            :: band
       integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
       real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
       logical(lgt),                                 intent(in),  optional :: alm_out
       real(dp),        dimension(:,:), allocatable                        :: evalBand
     end function evalBand

     ! Return component projected from map
     function projectBand(self, band, map, alm_in)
       import i4b, dp, comm_comp, comm_map, lgt
       class(comm_comp),                             intent(in)            :: self
       integer(i4b),                                 intent(in)            :: band
       class(comm_map),                              intent(in)            :: map
       logical(lgt),                                 intent(in), optional  :: alm_in
       real(dp),        dimension(:,:), allocatable                        :: projectBand
     end function projectBand

!!$     ! Evaluate amplitude map in brightness temperature at reference frequency
!!$     function evalAmp(self, nside, nmaps, pix, x_1D, x_2D)
!!$       import comm_comp, dp, i4b
!!$       class(comm_comp),                          intent(in)           :: self
!!$       integer(i4b),                              intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:),   intent(in), optional :: pix
!!$       real(dp),                  dimension(:),   intent(in), optional :: x_1D
!!$       real(dp),                  dimension(:,:), intent(in), optional :: x_2D
!!$       real(dp),     allocatable, dimension(:,:)                       :: evalAmp
!!$     end function evalAmp
!!$
!!$     ! Evaluate mixing matrix in brightness temperature at reference frequency
!!$     function evalMixmat(self, nside, nmaps, pix)
!!$       import comm_comp, dp, i4b
!!$       class(comm_comp),                        intent(in)           :: self
!!$       integer(i4b),                            intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:), intent(in), optional :: pix
!!$       real(dp),     allocatable, dimension(:,:)                     :: evalMixmat
!!$     end function evalMixmat
!!$
!!$     ! Generate simulated component
!!$     function simComp(self, handle, nside, nmaps, pix)
!!$       import planck_rng, comm_comp, dp, i4b
!!$       class(comm_comp),                        intent(in)           :: self
!!$       type(planck_rng),                        intent(inout)        :: handle
!!$       integer(i4b),                            intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:), intent(in), optional :: pix
!!$       real(dp),     allocatable, dimension(:,:)                     :: simComp
!!$     end function simComp
!!$
!!$     ! Dump current sample to HDF chain files
!!$     subroutine dumpHDF(self, filename)
!!$       import comm_comp
!!$       class(comm_comp),                        intent(in)           :: self
!!$       character(len=*),                        intent(in)           :: filename
!!$     end subroutine dumpHDF
!!$
     ! Dump current sample to HEALPix FITS file
     subroutine dumpFITS(self, iter, chainfile, output_hdf, postfix, dir)
       import comm_comp, i4b, hdf_file, lgt
       class(comm_comp),                        intent(in)           :: self
       integer(i4b),                            intent(in)           :: iter
       type(hdf_file),                          intent(in)           :: chainfile
       logical(lgt),                            intent(in)           :: output_hdf
       character(len=*),                        intent(in)           :: postfix
       character(len=*),                        intent(in)           :: dir
     end subroutine dumpFITS

     ! Initialize from HDF chain file
     subroutine initHDF(self, cpar, hdffile, hdfpath)
       import comm_comp, comm_params, hdf_file
       class(comm_comp),                        intent(inout)        :: self
       type(comm_params),                       intent(in)           :: cpar
       type(hdf_file),                          intent(in)           :: hdffile
       character(len=*),                        intent(in)           :: hdfpath
     end subroutine initHDF

     ! Sample spectral parameters
     subroutine sampleSpecInd(self, handle, id)
       import comm_comp, planck_rng, i4b
       class(comm_comp),                        intent(inout)        :: self
       type(planck_rng),                        intent(inout)        :: handle
       integer(i4b),                            intent(in)           :: id
     end subroutine sampleSpecInd

     ! Update mixing matrices
     subroutine updateMixmat(self, theta, beta)
       import comm_comp, comm_map, dp
       class(comm_comp),                        intent(inout)        :: self
       class(comm_map), dimension(:),           intent(in), optional :: theta
       real(dp),        dimension(:,:,:),       intent(in), optional :: beta
     end subroutine updateMixmat
       
  end interface

  !**************************************************
  !             Auxiliary variables
  !**************************************************

  integer(i4b)              :: n_dump = 1000
  real(dp),    dimension(2) :: nu_dump = [0.1d9, 3000.d9] 

  !**************************************************
  !             Internal module variables
  !**************************************************
  integer(i4b)              :: ncomp
  class(comm_comp), pointer :: compList, currComp => null()
  
contains

  subroutine initComp(self, cpar, id, id_abs)
    implicit none
    class(comm_comp)               :: self
    type(comm_params),  intent(in) :: cpar
    integer(i4b),       intent(in) :: id, id_abs

    self%id              = id
    self%active          = cpar%cs_include(id_abs)
    self%label           = cpar%cs_label(id_abs)
    self%type            = cpar%cs_type(id_abs)
    self%class           = cpar%cs_class(id_abs)
    self%unit            = cpar%cs_unit(id_abs)
    self%nu_ref          = cpar%cs_nu_ref(id_abs)
    self%myid            = cpar%myid_chain    
    self%numprocs        = cpar%numprocs_chain
    self%cg_samp_group   = cpar%cs_cg_samp_group(id_abs)
    self%init_from_HDF   = cpar%cs_initHDF(id_abs)
    self%operation       = cpar%operation

    ! Set up conversion factor between RJ and native component unit
    select case (trim(self%unit))
    case ('uK_cmb')
       self%RJ2unit_ = comp_a2t(self%nu_ref)
    case ('MJy/sr') 
       self%RJ2unit_ = comp_bnu_prime_RJ(self%nu_ref) * 1e14
    case ('K km/s') 
       self%RJ2unit_ = -1.d30
    case ('y_SZ') 
       self%RJ2unit_ = 2.d0*self%nu_ref**2*k_b/c**2 / &
               & (comp_bnu_prime(self%nu_ref) * comp_sz_thermo(self%nu_ref))
    case ('uK_RJ') 
       self%RJ2unit_ = 1.d0
    case default
       call report_error('Unsupported unit: ' // trim(self%unit))
    end select

  end subroutine initComp

  subroutine dumpSED(self, unit)
    implicit none
    class(comm_comp), intent(in) :: self
    integer(i4b),     intent(in) :: unit

    integer(i4b) :: i
    real(dp)     :: nu, dnu, S

    if (trim(self%type) == 'line') then
       do i = 1, numband
          S = self%S(band=i, theta=self%theta_def(1:self%npar))
          if (S /= 0.d0) write(unit,*) data(i)%bp%nu_c*1d-9, S
       end do
    else
       nu = nu_dump(1)
       dnu = (nu_dump(2)/nu_dump(1))**(1.d0/(n_dump-1))
       do i = 1, n_dump
          if (self%npar > 0) then
             S = self%S(nu, theta=self%theta_def(1:self%npar))
          else
             S = self%S(nu)
          end if
          if (S /= 0.d0) write(unit,*) nu*1d-9, S
          nu = nu*dnu
       end do
    end if
    
  end subroutine dumpSED

!!$  subroutine dumpCompMaps(postfix, dir)
!!$    implicit none
!!$    character(len=*), intent(in) :: postfix, dir
!!$
!!$    class(comm_comp), pointer :: c
!!$
!!$    c => compList
!!$    do while (associated(c))
!!$       call c%dumpFITS(postfix, dir)
!!$       c => c%next()
!!$    end do
!!$    
!!$  end subroutine dumpCompMaps

  function RJ2unit(self, bp)
    implicit none

    class(comm_comp), intent(in)           :: self
    class(comm_bp),   intent(in), optional :: bp
    real(dp)                               :: RJ2unit

    if (present(bp)) then
       if (trim(self%unit) == 'K km/s') then
          RJ2unit = 1.d0 / bp%lineAmp_RJ(self%nu_ref)
       else
          RJ2unit = self%RJ2unit_
       end if
    else
       RJ2unit = self%RJ2unit_
    end if
    
  end function RJ2unit

  subroutine CG_mask(self, samp_group, mask)
    implicit none

    class(comm_comp),               intent(in)    :: self
    integer(i4b),                   intent(in)    :: samp_group
    real(dp),         dimension(:), intent(inout) :: mask

    call cr_mask(self%id, samp_group==self%cg_samp_group, mask)

  end subroutine CG_mask
  
  function next(self)
    class(comm_comp) :: self
    class(comm_comp), pointer :: next
    next => self%nextLink
  end function next

  function prev(self)
    class(comm_comp) :: self
    class(comm_comp), pointer :: prev
    prev => self%prevLink
  end function prev
  
  subroutine setNext(self,next)
    class(comm_comp) :: self
    class(comm_comp), pointer :: next
    self%nextLink => next
  end subroutine setNext

  subroutine add(self,link)
    class(comm_comp), target  :: self
    class(comm_comp), pointer :: link

    class(comm_comp), pointer :: c
    
    c => self
    do while (associated(c%nextLink))
       c => c%nextLink
    end do
    link%prevLink => c
    c%nextLink    => link
  end subroutine add
  
end module comm_comp_mod
