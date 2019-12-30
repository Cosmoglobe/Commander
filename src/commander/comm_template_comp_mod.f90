module comm_template_comp_mod
  use math_tools
  use comm_param_mod
  use comm_comp_mod
  use comm_F_int_mod
  use comm_F_int_0D_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  use pix_tools
  use comm_hdf_mod
  use comm_cr_utils
  use comm_cr_precond_mod
  use locate_mod
  use spline_1D_mod
  implicit none

  private
  public comm_template_comp, initTemplatePrecond, updateTemplatePrecond, applyTemplatePrecond, &
       & initialize_template_comps

  !**************************************************
  !            Compact object class
  !**************************************************
  type template

  end type template
  
  type, extends (comm_comp) :: comm_template_comp
     character(len=512)         :: outprefix
     real(dp)                   :: cg_scale
     integer(i4b)               :: band
     real(dp),  dimension(1,1)  :: x         ! Amplitude
     real(dp),  dimension(2)    :: P         ! Gaussian prior on amplitude (mean,sigma)
     real(dp),  dimension(2)    :: P_cg      ! Gaussian prior on amplitude for CG (mean,sigma)
     class(comm_map), pointer   :: T    => null()   ! Template
     class(comm_map), pointer   :: mask => null()   ! Template mask
   contains
     procedure :: dumpFITS      => dumpTemplateToFITS
     procedure :: getBand       => evalTemplateBand
     procedure :: projectBand   => projectTemplateBand
     procedure :: S             => evalSED
     procedure :: initHDF       => initTemplateHDF
     procedure :: sampleSpecInd => sampleTempSpecInd
     procedure :: updateMixmat  => updateTempMixmat
     procedure :: update_F_int  => updateTempFInt
  end type comm_template_comp

  interface comm_template_comp
     procedure constructor
  end interface comm_template_comp

  type template_ptr
     class(comm_template_comp), pointer :: p => null()
  end type template_ptr
  
  integer(i4b) :: npre         =   0
  integer(i4b) :: comm_pre     =  -1
  integer(i4b) :: myid_pre     =  -1
  integer(i4b) :: numprocs_pre =  -1
  logical(lgt) :: recompute_template_precond = .true.
  
contains

  function constructor(cpar, id, id_abs, band, label, mapfile, maskfile, mu, rms, def)
    implicit none
    class(comm_params),       intent(in) :: cpar
    integer(i4b),             intent(in) :: id, id_abs, band
    character(len=*),         intent(in) :: label, mapfile, maskfile
    real(dp),                 intent(in) :: mu, rms, def
    class(comm_template_comp),   pointer    :: constructor

    character(len=512) :: dir

    ! General parameters
    allocate(constructor)

    ! Initialize general parameters
    constructor%class     = cpar%cs_class(id_abs)
    constructor%type      = cpar%cs_type(id_abs)
    constructor%label     = cpar%cs_label(id_abs)
    constructor%id        = id
    constructor%nmaps     = 1    ! Only used for CR book-keeping; must be 1 for templates
    constructor%outprefix = trim(cpar%cs_label(id_abs))
    constructor%cg_scale  = 1.d0
    constructor%cg_samp_group  = cpar%cs_cg_samp_group(id_abs)
    constructor%myid      = cpar%myid
    constructor%comm      = cpar%comm_chain
    constructor%numprocs  = cpar%numprocs_chain
    constructor%P         = [mu,rms]
    constructor%band      = band
    npre                  = npre + 1
    comm_pre              = cpar%comm_chain
    myid_pre              = cpar%myid
    numprocs_pre          = cpar%numprocs_chain

    if (constructor%myid == 0) then
       constructor%ncr = 1
       constructor%x   = def
    else
       constructor%ncr = 0
    end if

    ! Read template and mask
    dir = trim(cpar%datadir) // '/'
    constructor%T => comm_map(data(band)%info, trim(dir)//trim(mapfile))
    if (trim(maskfile) /= 'fullsky') then
       constructor%mask  => comm_map(data(band)%info, trim(dir)//trim(maskfile))
       constructor%P_cg  =  [mu,1.d-6]
    else
       constructor%P_cg  =  constructor%P      
    end if

  end function constructor


  function initialize_template_comps(cpar, id, id_abs, n)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_template_comp), pointer   :: initialize_template_comps

    integer(i4b)        :: i, unit
    real(dp)            :: mu, rms, def
    character(len=1024) :: line, label, mapfile, maskfile
    class(comm_comp), pointer :: c => null()

    unit  = getlun()

    ! Find number of lines
    n = 0
    open(unit, file=trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       read(line,*) label, mapfile, maskfile, mu, rms, def
       do i = 1, numband
          if (trim(label) == trim(data(i)%label)) exit
       end do
       if (i > numband) cycle

       if (n == 0) then
          initialize_template_comps => comm_template_comp(cpar, id+n, id_abs, i, label, mapfile, maskfile, &
               & mu, rms, def)
       else
          c => comm_template_comp(cpar, id+n, id_abs, i, label, mapfile, maskfile, mu, rms, def)
          call initialize_template_comps%add(c)
       end if
       n = n+1
    end do
1   close(unit)
  
  end function initialize_template_comps

  function evalSED(self, nu, band, pol, theta)
    class(comm_template_comp),  intent(in)           :: self
    real(dp),                   intent(in), optional :: nu
    integer(i4b),               intent(in), optional :: band
    integer(i4b),               intent(in), optional :: pol
    real(dp), dimension(1:),    intent(in), optional :: theta
    real(dp)                                        :: evalSED

    if (band == self%band) then
       evalSED = 1.d0
    else
       evalSED = 0.d0
    end if
    
  end function evalSED

  function evalTemplateBand(self, band, amp_in, pix, alm_out, det)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: evalTemplateBand

    integer(i4b) :: i, j, p, q, ierr
    real(dp)     :: a

    if (.not. allocated(evalTemplateBand)) &
         & allocate(evalTemplateBand(0:data(band)%info%np-1,data(band)%info%nmaps))

    if (band /= self%band) then
       evalTemplateBand = 0.d0
       return
    end if

    if (self%myid == 0) then
       if (present(amp_in)) then
          a = amp_in(1,1)
       else
          a = self%x(1,1)
       end if
    end if
    call mpi_bcast(a, 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

    evalTemplateBand = a * self%T%map
    
  end function evalTemplateBand
  
  ! Return component projected from map
  function projectTemplateBand(self, band, map, alm_in, det)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: projectTemplateBand

    integer(i4b) :: i, j, q, p, ierr, d
    real(dp)     :: val, val2
    real(dp), allocatable, dimension(:,:) :: amp, amp2
    
    if (.not. allocated(projectTemplateBand)) allocate(projectTemplateBand(1,1))

    if (band /= self%band) then
       projectTemplateBand = 0.d0
       return
    end if

    val = sum(map%map * self%T%map)
    call mpi_reduce(val, val2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    if (self%myid == 0) projectTemplateBand = val2

  end function projectTemplateBand
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpTemplateToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    class(comm_template_comp),               intent(in)           :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    if (self%myid == 0) write(*,*) '     Template amplitude = ', self%x

    return
    
  end subroutine dumpTemplateToFITS

  ! Dump current sample to HEALPix FITS file
  subroutine initTemplateHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_template_comp),    intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath

    integer(i4b)       :: i, j
    real(dp)           :: md(4)
    character(len=512) :: path
    real(dp), allocatable, dimension(:,:,:) :: theta

!!$    path = trim(adjustl(hdfpath))//trim(adjustl(self%label)) // '/'
!!$    if (self%myid == 0) then
!!$       call read_hdf(hdffile, trim(adjustl(path))//'/amp', self%x)
!!$       self%x = self%x/self%cg_scale
!!$       
!!$       allocate(theta(self%nsrc,self%nmaps,self%npar))
!!$       call read_hdf(hdffile, trim(adjustl(path))//'/specind', theta)
!!$       do i = 1, self%nsrc
!!$          do j = 1, self%nmaps
!!$             self%src(i)%theta(:,j) = theta(i,j,:) 
!!$          end do
!!$       end do
!!$       deallocate(theta)
!!$    end if
  end subroutine initTemplateHDF


  subroutine initTemplatePrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, j1, j2, k1, k2, q, l, m, n, p, p1, p2, n1, n2, myid, ierr, cnt
    real(dp)     :: t1, t2
    logical(lgt) :: skip
    class(comm_comp),          pointer :: c => null(), c1 => null(), c2 => null()
    class(comm_template_comp), pointer :: pt1 => null(), pt2 => null()
    class(comm_map),           pointer :: invN_T => null()
    real(dp),     allocatable, dimension(:,:) :: mat, mat2

    if (npre == 0) return
    if (allocated(P_cr%invM_temp)) return

    call mpi_comm_rank(comm, myid, ierr)
        
    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(mat(npre,npre), mat2(npre,npre))
    mat = 0.d0
    i1  = 0
    c1 => compList
    do while (associated(c1))
       skip = .true.
       select type (c1)
       class is (comm_template_comp)
          pt1  => c1
          skip = .false.
       end select
       if (skip) then
          c1 => c1%next()
          cycle
       end if
       i1 = i1+1

       invN_T => comm_map(pt1%T)
       call data(pt1%band)%N%invN(invN_T) ! Multiply with invN

       i2 = 0
       c2 => compList
       do while (associated(c2))
          skip = .true.
          select type (c2)
          class is (comm_template_comp)
             pt2 => c2
             skip = .false.
          end select
          if (skip) then
             c2 => c2%next()
             cycle
          end if
          i2 = i2+1
          !if (i2 < i1) cycle

          mat(i1,i2) = sum(invN_T%map * pt2%T%map)
          !mat(i2,i1) = mat(i1,i2)
          c2 => c2%next()
       end do
       c1 => c1%next()
    end do

    ! Collect contributions from all cores
    call mpi_reduce(mat, mat2, size(mat2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    ! Invert matrix to finalize preconditioner
    allocate(P_cr%invM_temp(1,1))
    if (myid == 0) then
       ! Multiply with sqrtS from left side
       i1  = 0
       c1 => compList
       do while (associated(c1))
          skip = .true.
          select type (c1)
          class is (comm_template_comp)
             pt1  => c1
             skip = .false.
          end select
          if (skip) then
             c1 => c1%next()
             cycle
          end if
          i1         = i1+1
          mat2(i1,:) = mat2(i1,:) * pt1%P_cg(2)
          c1 => c1%next()
       end do
       ! Multiply with sqrtS from right side
       i1  = 0
       c1 => compList
       do while (associated(c1))
          skip = .true.
          select type (c1)
          class is (comm_template_comp)
             pt1  => c1
             skip = .false.
          end select
          if (skip) then
             c1 => c1%next()
             cycle
          end if
          i1         = i1+1
          mat2(:,i1) = mat2(:,i1) * pt1%P_cg(2)
          c1 => c1%next()
       end do
       ! Add unity
       do i1 = 1, npre
          mat2(i1,i1) = mat2(i1,i1) + 1.d0
       end do
       ! Invert matrix to finalize preconditioner
       call invert_matrix_with_mask(mat2)
       allocate(P_cr%invM_temp(1,1)%M(npre,npre))
       P_cr%invM_temp(1,1)%M = mat2
    end if

    deallocate(mat,mat2)
    
  end subroutine initTemplatePrecond

  subroutine updateTemplatePrecond(samp_group)
    implicit none
    integer(i4b),          intent(in) :: samp_group

    ! Placeholder for now; already fully initialized
    if (npre == 0) return
       
  end subroutine updateTemplatePrecond


  subroutine applyTemplatePrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    logical(lgt)              :: skip
    real(dp),                  allocatable, dimension(:,:) :: amp
    real(dp),                  allocatable, dimension(:)   :: y
    class(comm_comp),          pointer                     :: c => null() 
    class(comm_template_comp), pointer                     :: pt => null()

    if (npre == 0 .or. myid_pre /= 0) return
    
    ! Reformat linear array into y(npre) structure
    allocate(y(npre))
    y = 0.d0
    l = 1
    c => compList
    do while (associated(c))
       skip = .true.
       select type (c)
       class is (comm_template_comp)
          pt => c
          skip = .false.
       end select
       if (skip) then
          c => c%next()
          cycle
       end if
       call cr_extract_comp(pt%id, x, amp)
       y(l) = amp(0,1)
       l  = l + 1
       c => c%next()
       deallocate(amp)
    end do

    ! Multiply with preconditioner
    y = matmul(P_cr%invM_temp(1,1)%M, y)

    ! Reformat y(npre) structure back into linear array
    l = 1
    c => compList
    do while (associated(c))
       skip = .true.
       select type (c)
       class is (comm_template_comp)
          pt => c
          skip = .false.
       end select
       if (skip) then
          c => c%next()
          cycle
       end if
       allocate(amp(0:0,1:1))
       amp(0,1) = y(l)
       call cr_insert_comp(pt%id, .false., amp, x)
       l = l + 1
       c => c%next()
       deallocate(amp)
    end do
    
    deallocate(y)

  end subroutine applyTemplatePrecond

  ! Sample spectral parameters
  subroutine sampleTempSpecInd(self, handle, id)
    implicit none
    class(comm_template_comp),               intent(inout)        :: self
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
  end subroutine sampleTempSpecInd


  subroutine updateTempMixmat(self, theta, beta, band)
    implicit none
    class(comm_template_comp),         intent(inout)        :: self
    class(comm_map), dimension(:),     intent(in), optional :: theta
    real(dp),        dimension(:,:,:), intent(in), optional :: beta  ! (npar,nmaps,nsrc)
    integer(i4b),                      intent(in), optional :: band

  end subroutine updateTempMixmat

  subroutine updateTempFInt(self, band)
    implicit none
    class(comm_template_comp), intent(inout)          :: self
    integer(i4b),              intent(in),   optional :: band

  end subroutine updateTempFInt

  
end module comm_template_comp_mod
