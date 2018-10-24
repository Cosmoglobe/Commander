module comm_map_mod
  use sharp
  use healpix_types
  use fitstools
  use pix_tools
  use udgrade_nr
  use iso_c_binding, only : c_ptr, c_double
  use head_fits
  use comm_hdf_mod
  use extension
  use comm_param_mod
  implicit none

!  include "mpif.h"
      
  private
  public comm_map, comm_mapinfo, map_ptr

  type :: comm_mapinfo
     ! Linked list variables
     class(comm_mapinfo), pointer :: nextLink => null()
     class(comm_mapinfo), pointer :: prevLink => null()

     ! Data variables
     type(sharp_alm_info)  :: alm_info
     type(sharp_geom_info) :: geom_info_T, geom_info_P
     logical(lgt) :: pol
     integer(i4b) :: comm, myid, nprocs
     integer(i4b) :: nside, npix, nmaps, nspec, nring, np, lmax, nm, nalm
     integer(c_int), allocatable, dimension(:)   :: rings
     integer(c_int), allocatable, dimension(:)   :: ms
     integer(c_int), allocatable, dimension(:)   :: mind
     integer(c_int), allocatable, dimension(:,:) :: lm
     integer(c_int), allocatable, dimension(:)   :: pix
     real(c_double), allocatable, dimension(:,:) :: W
   contains
     procedure     :: lm2i
     procedure     :: i2lm
     procedure     :: dealloc => comm_mapinfo_finalize
  end type comm_mapinfo

  type :: comm_map
     ! Linked list variables
     class(comm_map), pointer :: nextLink => null()
     class(comm_map), pointer :: prevLink => null()

     ! Data variables
     class(comm_mapinfo), pointer :: info
     real(c_double), allocatable, dimension(:,:) :: map
     real(c_double), allocatable, dimension(:,:) :: alm
   contains
     ! Data routines
     procedure     :: Y    => exec_sharp_Y
     procedure     :: Yt   => exec_sharp_Yt
     procedure     :: YtW  => exec_sharp_YtW
     procedure     :: WY   => exec_sharp_WY
     procedure     :: Y_EB => exec_sharp_Y_EB
     procedure     :: Y_scalar    => exec_sharp_Y_scalar
     procedure     :: Yt_scalar   => exec_sharp_Yt_scalar
     procedure     :: YtW_scalar  => exec_sharp_YtW_scalar
     procedure     :: writeFITS
     procedure     :: readFITS
     procedure     :: readHDF
     procedure     :: dealloc => deallocate_comm_map
     procedure     :: alm_equal
     procedure     :: add_alm
     procedure     :: udgrade
     procedure     :: getSigmaL
     procedure     :: getCrossSigmaL
     procedure     :: smooth

     ! Linked list procedures
     procedure :: next    ! get the link after this link
     procedure :: prev    ! get the link before this link
     procedure :: setNext ! set the link after this link
     procedure :: add     ! add new link at the end
  end type comm_map

  type map_ptr
     class(comm_map), pointer :: p
  end type map_ptr
  
  
  interface comm_mapinfo
     procedure constructor_mapinfo
  end interface comm_mapinfo

  interface comm_map
     procedure constructor_map, constructor_clone
  end interface comm_map

  ! Library of mapinfos; resident in memory during the analysis
  class(comm_mapinfo), pointer, private :: mapinfos

contains

  !**************************************************
  !             Constructors
  !**************************************************
  function constructor_mapinfo(comm, nside, lmax, nmaps, pol)
    implicit none
    integer(i4b),                 intent(in) :: comm, nside, lmax, nmaps
    logical(lgt),                 intent(in) :: pol
    class(comm_mapinfo), pointer             :: constructor_mapinfo

    integer(i4b) :: myid, nprocs, ierr
    integer(i4b) :: l, m, i, j, k, iring, np, ind 
    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b), allocatable, dimension(:) :: pixlist
    real(dp),     allocatable, dimension(:,:) :: weights
    character(len=5)   :: nside_text
    character(len=512) :: weight_file, healpixdir
    class(comm_mapinfo), pointer :: p, p_prev, p_new

    ! Check if requested mapinfo already exists in library; if so return a link to that object
    p => mapinfos
    do while (associated(p)) 
       if (p%nside == nside .and. p%lmax == lmax .and. p%nmaps == nmaps .and. p%pol == pol) then
          constructor_mapinfo => p
          return
       else
          p_prev => p
          p      => p%nextLink
       end if
    end do

    ! Set up new mapinfo object
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    allocate(p_new)
    p_new%comm   = comm
    p_new%myid   = myid
    p_new%nprocs = nprocs
    p_new%nside  = nside
    p_new%nmaps  = nmaps
    p_new%nspec  = nmaps*(nmaps+1)/2
    p_new%lmax   = lmax
    p_new%pol    = pol
    p_new%npix   = 12*nside**2

    ! Select rings and pixels
    allocate(pixlist(0:4*nside-1))
    p_new%nring = 0
    p_new%np    = 0
    do i = 1+myid, 2*nside, nprocs
       call in_ring(nside, i, 0.d0, pi, pixlist, np)
       p_new%nring = p_new%nring + 1
       p_new%np    = p_new%np    + np
       if (i < 2*nside) then ! Symmetric about equator
          p_new%nring = p_new%nring + 1
          p_new%np    = p_new%np    + np
       end if
    end do
    allocate(p_new%rings(p_new%nring))
    allocate(p_new%pix(p_new%np))
    j = 1
    k = 1
    do i = 1+myid, 2*nside, nprocs
       call in_ring(nside, i, 0.d0, pi, pixlist, np)
       p_new%rings(k) = i
       p_new%pix(j:j+np-1) = pixlist(0:np-1)
       k = k + 1
       j = j + np
       if (i < 2*nside) then ! Symmetric about equator
          call in_ring(nside, 4*nside-i, 0.d0, pi, pixlist, np)
          p_new%rings(k) = 4*nside-i
          p_new%pix(j:j+np-1) = pixlist(0:np-1)
          k = k + 1
          j = j + np
       end if
    end do
    deallocate(pixlist)
    call QuickSort_int(p_new%rings)
    call QuickSort_int(p_new%pix)

    ! Select m's
    p_new%nm = 0
    do m = myid, lmax, nprocs
       p_new%nm   = p_new%nm   + 1
       if (m == 0) then
          p_new%nalm = p_new%nalm + lmax+1
       else
          p_new%nalm = p_new%nalm + 2*(lmax-m+1)
       end if
    end do
    allocate(p_new%ms(p_new%nm))
    allocate(p_new%mind(0:p_new%lmax))
    allocate(p_new%lm(2,0:p_new%nalm-1))
    ind = 0
    p_new%mind = -1
    do i = 1, p_new%nm
       m                           = myid + (i-1)*nprocs
       p_new%ms(i)   = m
       p_new%mind(m) = ind
       if (m == 0) then
          do l = m, lmax
             p_new%lm(:,ind) = [l,m]
             ind                           = ind+1
          end do
       else
          do l = m, lmax
             p_new%lm(:,ind) = [l,+m]
             ind                           = ind+1
             p_new%lm(:,ind) = [l,-m]
             ind                           = ind+1
          end do
       end if
    end do
    
    ! Read ring weights, and create SHARP info structures
    call sharp_make_mmajor_real_packed_alm_info(lmax, ms=p_new%ms, &
         & alm_info=p_new%alm_info)
    call getEnvironment('HEALPIX', healpixdir) 
    call int2string(nside, nside_text)
    weight_file = trim(healpixdir)//'/data/weight_ring_n' // nside_text // '.fits'
    if (nmaps == 1) then
       allocate(p_new%W(2*nside,1))
       call read_dbintab(weight_file, p_new%W, 2*nside, 1, nullval, anynull)
       p_new%W = p_new%W + 1.d0
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,1), geom_info=p_new%geom_info_T)
    else
       allocate(p_new%W(2*nside,2))
       call read_dbintab(weight_file, p_new%W, 2*nside, 2, nullval, anynull)
       p_new%W(:,:) = p_new%W + 1.d0
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,1), geom_info=p_new%geom_info_T)
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,2), geom_info=p_new%geom_info_P)
    end if

    ! Set up new object, and add to list
    p => mapinfos
    if (.not. associated(p)) then
       mapinfos => p_new
    else
       do while (associated(p%nextLink)) 
          p => p%nextLink
       end do
       p%nextLink => p_new
    end if
    constructor_mapinfo => p_new

  end function constructor_mapinfo

  function constructor_map(info, filename, mask_misspix, udgrade)
    implicit none
    class(comm_mapinfo),                              intent(in),  target   :: info
    character(len=*),                                 intent(in),  optional :: filename
    logical(lgt),                                     intent(in),  optional :: udgrade
    real(dp),            allocatable, dimension(:,:), intent(out), optional :: mask_misspix
    class(comm_map),     pointer                                            :: constructor_map

    allocate(constructor_map)
    constructor_map%info => info
    allocate(constructor_map%map(0:info%np-1,info%nmaps))
    allocate(constructor_map%alm(0:info%nalm-1,info%nmaps))

    if (present(filename)) then
       if (present(mask_misspix)) then
          allocate(mask_misspix(0:info%np-1,info%nmaps))
          call constructor_map%readFITS(filename, mask=mask_misspix, udgrade=udgrade)
       else
          call constructor_map%readFITS(filename, udgrade=udgrade)
       end if
    else
       constructor_map%map = 0.d0
    end if
    constructor_map%alm = 0.d0
    
  end function constructor_map
  
  function constructor_clone(map)
    implicit none
    class(comm_map),         intent(in)  :: map
    class(comm_map), pointer             :: constructor_clone

    allocate(constructor_clone)
    constructor_clone%info => map%info
    allocate(constructor_clone%map(0:map%info%np-1,map%info%nmaps))
    allocate(constructor_clone%alm(0:map%info%nalm-1,map%info%nmaps))
    constructor_clone%map = map%map
    constructor_clone%alm = map%alm
    
  end function constructor_clone

  subroutine deallocate_comm_map(self, clean_info)
    implicit none

    class(comm_map), intent(inout)          :: self
    logical(lgt),    intent(in),   optional :: clean_info
    class(comm_map), pointer :: link

    logical(lgt) :: clean_info_

    clean_info_ = .false. !; if (present(clean_info)) clean_info_ = clean_info ! Never clean info with new library
    if (allocated(self%map)) deallocate(self%map)
    if (allocated(self%alm)) deallocate(self%alm)
    if (clean_info_ .and. associated(self%info)) call self%info%dealloc()
    nullify(self%info)

    if (associated(self%nextLink)) then
       ! Deallocate all links
       link => self%nextLink
       do while (associated(link))
          if (allocated(link%map)) deallocate(link%map)
          if (allocated(link%alm)) deallocate(link%alm)
          if (clean_info_ .and. associated(link%info)) call link%info%dealloc()
          nullify(link%info)
          link => link%nextLink
       end do
       nullify(self%nextLink)
    end if

  end subroutine deallocate_comm_map

  subroutine comm_mapinfo_finalize(self)
    implicit none

    class(comm_mapinfo) :: self
    
    if (allocated(self%rings)) then
       deallocate(self%rings, self%ms, self%mind, self%lm, self%pix, self%W)
       call sharp_destroy_alm_info(self%alm_info)
       call sharp_destroy_geom_info(self%geom_info_T)
       if (self%nmaps == 3) call sharp_destroy_geom_info(self%geom_info_P)
    end if

  end subroutine comm_mapinfo_finalize
  
  !**************************************************
  !             Spherical harmonic transforms
  !**************************************************

  subroutine exec_sharp_Y(self)
    implicit none

    class(comm_map), intent(inout)          :: self

    if (.not. allocated(self%map)) allocate(self%map(0:self%info%np-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_Y, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info_T, comm=self%info%comm)
       if (self%info%nmaps == 3) then
          call sharp_execute(SHARP_Y, 2, 2, self%alm(:,2:3), self%info%alm_info, &
               & self%map(:,2:3), self%info%geom_info_P, comm=self%info%comm)
       end if
    else
       call sharp_execute(SHARP_Y, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info_T, comm=self%info%comm)       
    end if
    
  end subroutine exec_sharp_Y

  subroutine exec_sharp_WY(self)
    implicit none

    class(comm_map), intent(inout)          :: self

    if (.not. allocated(self%map)) allocate(self%map(0:self%info%np-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_WY, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info_T, comm=self%info%comm)
       if (self%info%nmaps == 3) then
          call sharp_execute(SHARP_WY, 2, 2, self%alm(:,2:3), self%info%alm_info, &
               & self%map(:,2:3), self%info%geom_info_P, comm=self%info%comm)
       end if
    else
       call sharp_execute(SHARP_WY, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info_T, comm=self%info%comm)       
    end if
    
  end subroutine exec_sharp_WY

  subroutine exec_sharp_Y_scalar(self)
    implicit none

    class(comm_map), intent(inout)          :: self
    integer(i4b) :: i

    if (.not. allocated(self%map)) allocate(self%map(0:self%info%np-1,self%info%nmaps))
    do i = 1, self%info%nmaps
       call sharp_execute(SHARP_Y, 0, 1, self%alm(:,i:i), self%info%alm_info, &
            & self%map(:,i:i), self%info%geom_info_T, comm=self%info%comm)
    end do
    
  end subroutine exec_sharp_Y_scalar
  
  subroutine exec_sharp_Y_EB(self)
    implicit none

    class(comm_map), intent(inout)          :: self

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info
    
    info => comm_mapinfo(self%info%comm, self%info%nside, self%info%lmax, self%info%nmaps, .false.)

    if (.not. allocated(self%map)) allocate(self%map(0:self%info%np-1,self%info%nmaps))
    do i = 1, self%info%nmaps
       call sharp_execute(SHARP_Y, 0, 1, self%alm(:,i:i), info%alm_info, &
            & self%map(:,i:i), info%geom_info_T, comm=info%comm)
    end do

    deallocate(info)
    
  end subroutine exec_sharp_Y_EB

  subroutine exec_sharp_Yt(self)
    implicit none

    class(comm_map), intent(inout) :: self

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_Yt, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info_T, comm=self%info%comm)
       if (self%info%nmaps == 3) then
          call sharp_execute(SHARP_Yt, 2, 2, self%alm(:,2:3), self%info%alm_info, &
               & self%map(:,2:3), self%info%geom_info_P, comm=self%info%comm)
       end if
       
    else
       call sharp_execute(SHARP_Yt, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info_T, comm=self%info%comm)       
    end if
    
  end subroutine exec_sharp_Yt

  subroutine exec_sharp_Yt_scalar(self)
    implicit none

    class(comm_map), intent(inout) :: self
    integer(i4b) :: i

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    do i = 1, self%info%nmaps
       call sharp_execute(SHARP_Yt, 0, 1, self%alm(:,i:i), self%info%alm_info, &
            & self%map(:,i:i), self%info%geom_info_T, comm=self%info%comm)
    end do
    
  end subroutine exec_sharp_Yt_scalar

  subroutine exec_sharp_YtW(self)
    implicit none

    class(comm_map), intent(inout) :: self

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_YtW, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info_T, comm=self%info%comm)
       if (self%info%nmaps == 3) then
          call sharp_execute(SHARP_YtW, 2, 2, self%alm(:,2:3), self%info%alm_info, &
               & self%map(:,2:3), self%info%geom_info_P, comm=self%info%comm)
       end if
    else
       call sharp_execute(SHARP_YtW, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info_T, comm=self%info%comm)
    end if
    
  end subroutine exec_sharp_YtW


  subroutine exec_sharp_YtW_scalar(self)
    implicit none

    class(comm_map), intent(inout) :: self
    integer(i4b) :: i

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    do i = 1, self%info%nmaps
       call sharp_execute(SHARP_YtW, 0, 1, self%alm(:,i:i), self%info%alm_info, &
            & self%map(:,i:i), self%info%geom_info_T, comm=self%info%comm)
    end do
    
  end subroutine exec_sharp_YtW_scalar
  
  !**************************************************
  !                   IO routines
  !**************************************************

  subroutine writeFITS(self, filename, comptype, nu_ref, unit, ttype, spectrumfile, &
       & hdffile, hdfpath)
    implicit none

    class(comm_map),  intent(in) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),         intent(in), optional :: nu_ref
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath

    integer(i4b) :: i, j, l, m, ind, nmaps, npix, np, nlm, ierr
    real(dp),     allocatable, dimension(:,:) :: map, alm, buffer
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), allocatable, dimension(:,:) :: lm
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat
    
    ! Only the root actually writes to disk; data are distributed via MPI
    npix  = self%info%npix
    nmaps = self%info%nmaps
    if (self%info%myid == 0) then

       ! Distribute to other nodes
       call update_status(status, "fits1")
       allocate(p(npix), map(0:npix-1,nmaps))
       map(self%info%pix,:) = self%map
       do i = 1, self%info%nprocs-1
          call mpi_recv(np,       1, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          call mpi_recv(p(1:np), np, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          allocate(buffer(np,nmaps))
          call mpi_recv(buffer, np*nmaps, &
               & MPI_DOUBLE_PRECISION, i, 98, self%info%comm, mpistat, ierr)
          map(p(1:np),:) = buffer(1:np,:)
          deallocate(buffer)
       end do
       call update_status(status, "fits2")
       call write_map(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
              call update_status(status, "fits3")

       if (present(hdffile)) then
          allocate(alm(0:(self%info%lmax+1)**2-1,self%info%nmaps))
          do j = 0, self%info%nalm-1
             l   = self%info%lm(1,j)
             m   = self%info%lm(2,j)
             ind = l**2 + l + m
             alm(ind,:) = self%alm(j,:)
          end do
          do i = 1, self%info%nprocs-1
             call mpi_recv(nlm,      1, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
             allocate(lm(2,0:nlm-1))
             call mpi_recv(lm, size(lm), MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
             allocate(buffer(0:nlm-1,nmaps))
             call mpi_recv(buffer, size(buffer), &
                  & MPI_DOUBLE_PRECISION, i, 98, self%info%comm, mpistat, ierr)
             do j = 0, nlm-1
                l   = lm(1,j)
                m   = lm(2,j)
                ind = l**2 + l + m                
                alm(ind,:) = buffer(j,:)
             end do
             deallocate(lm, buffer)
          end do
          call update_status(status, "fits4")
          call write_hdf(hdffile, trim(adjustl(hdfpath)//'alm'),   alm)
          call write_hdf(hdffile, trim(adjustl(hdfpath)//'map'),   map)
          call write_hdf(hdffile, trim(adjustl(hdfpath)//'lmax'),  self%info%lmax)
          call write_hdf(hdffile, trim(adjustl(hdfpath)//'nmaps'), self%info%nmaps)
          call update_status(status, "fits5")
          deallocate(alm)
       end if

       deallocate(p, map)

    else
       call mpi_send(self%info%np,  1,              MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%info%pix, self%info%np,   MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%map,      size(self%map), MPI_DOUBLE_PRECISION, 0, 98, &
            & self%info%comm, ierr)

       if (present(hdffile)) then
          call mpi_send(self%info%nalm, 1,                  MPI_INTEGER, 0, 98, self%info%comm, ierr)
          call mpi_send(self%info%lm,   size(self%info%lm), MPI_INTEGER, 0, 98, self%info%comm, ierr)
          call mpi_send(self%alm,       size(self%alm),     MPI_DOUBLE_PRECISION, 0, 98, &
               & self%info%comm, ierr)
       end if
    end if
    
  end subroutine writeFITS
  
  subroutine readFITS(self, filename, mask, udgrade)
    implicit none

    class(comm_map),                    intent(inout) :: self
    real(dp),         dimension(0:,1:), intent(out), optional :: mask
    logical(lgt),                       intent(in),  optional :: udgrade
    character(len=*),                   intent(in)    :: filename

    integer(i4b) :: i, j, np, npix, ordering, nside, nmaps, ierr, badpix
    real(dp),     allocatable, dimension(:,:) :: map, buffer, map_in
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat

    ! Check file consistency 
    npix = getsize_fits(trim(filename), ordering=ordering, nside=nside, nmaps=nmaps)
    if (nmaps < self%info%nmaps) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nmaps in ' // trim(filename)
       call mpi_finalize(ierr)
       stop
    end if
    if (nside /= self%info%nside .and. .not. (present(udgrade))) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nside in ' // trim(filename)
       call mpi_finalize(ierr)
       stop
    end if

    ! Only the root actually reads from disk; data are distributed via MPI
    if (self%info%myid == 0) then

       ! Read map and convert to RING format if necessary
       allocate(map(0:self%info%npix-1,self%info%nmaps))
       if (present(udgrade)) then
          allocate(map_in(0:npix-1,nmaps))
          call input_map(filename, map_in, npix, nmaps)
          if (ordering == 1) then
             call udgrade_ring(map_in, nside, map, nside_out=self%info%nside)
          else
             call udgrade_nest(map_in, nside, map, nside_out=self%info%nside)
          end if
          deallocate(map_in)
       else
          call input_map(filename, map, self%info%npix, self%info%nmaps)
       end if

       if (present(mask)) then
          badpix = count(abs(map+1.6375d30) < 1d25 .or. map > 1.d30)
          if (badpix > 0) then
             write(*,fmt='(a,i8)') '   Number of bad pixels = ', badpix
          end if
       end if
       
       if (ordering == 2) then
          do i = 1, self%info%nmaps
             call convert_nest2ring(self%info%nside, map(:,i))
          end do
       end if

       ! Distribute to other nodes
       allocate(p(self%info%npix))
       self%map = map(self%info%pix,:)
       do i = 1, self%info%nprocs-1
          call mpi_recv(np,       1, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          call mpi_recv(p(1:np), np, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          allocate(buffer(np,self%info%nmaps))
          buffer = map(p(1:np),:)
          call mpi_send(buffer, np*self%info%nmaps, MPI_DOUBLE_PRECISION, i, 98, &
               & self%info%comm, ierr)
          map(p(1:np),:) = buffer
          deallocate(buffer)
       end do
       deallocate(p, map)
    else
       call mpi_send(self%info%np,               1, MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%info%pix, self%info%np,   MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_recv(self%map,      size(self%map), MPI_DOUBLE_PRECISION, 0, 98, &
            &  self%info%comm, mpistat, ierr)
    end if

    if (present(mask)) then
       do j = 1, self%info%nmaps
          do i = 0, self%info%np-1
             if (abs(self%map(i,j)+1.6375d30) < 1d25 .or. self%map(i,j) == 0.d0 .or. self%map(i,j) > 1.d30) then
                mask(i,j) = 0.d0
             else
                mask(i,j) = 1.d0
             end if
          end do
       end do
    end if
    
  end subroutine readFITS

  subroutine readHDF(self, hdffile, hdfpath)
    implicit none

    class(comm_map),  intent(inout) :: self
    type(hdf_file),   intent(in)    :: hdffile
    character(len=*), intent(in)    :: hdfpath

    integer(i4b) :: i, l, m, j, lmax, nmaps, ierr, nalm, npix
    real(dp),     allocatable, dimension(:,:) :: alms, map
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat

    lmax  = self%info%lmax
    npix  = self%info%npix
    nmaps = self%info%nmaps
    nalm  = (lmax+1)**2

    ! Only the root actually reads from disk; data are distributed via MPI
    allocate(alms(0:nalm-1,nmaps))
    if (self%info%myid == 0) call read_hdf(hdffile, trim(adjustl(hdfpath))//'alm', alms)
    call mpi_bcast(alms, size(alms),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
    do i = 0, self%info%nalm-1
       call self%info%i2lm(i, l, m)
       j = l**2 + l + m
       self%alm(i,:) = alms(j,:)
    end do
    deallocate(alms)

    ! Only the root actually reads from disk; data are distributed via MPI
    allocate(map(0:npix-1,nmaps))
    if (self%info%myid == 0) call read_hdf(hdffile, trim(adjustl(hdfpath))//'map', map)
    call mpi_bcast(map, size(map),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
    do i = 1, nmaps
       self%map(:,i) = map(self%info%pix,i)
    end do
    deallocate(map)
    
  end subroutine readHDF

  subroutine write_map(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map
    character(len=*),                   intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),                           intent(in), optional :: nu_ref

    integer(i4b)   :: npix, nlheader, nmaps, i, nside
    logical(lgt)   :: exist, polarization

    character(len=80), dimension(1:120)    :: header
    character(len=16) :: unit_, ttype_

    npix         = size(map(:,1))
    nside        = nint(sqrt(real(npix,sp)/12.))
    nmaps        = size(map(0,:))
    polarization = (nmaps == 3)
    unit_        = '';       if (present(unit)) unit_  = unit
    ttype_       = 'Stokes'; if (present(unit)) ttype_ = ttype


    !-----------------------------------------------------------------------
    !                      write the map to FITS file
    !  This is copied from the synfast.f90 file in the Healpix package
    !-----------------------------------------------------------------------
    
    nlheader = SIZE(header)
    do i=1,nlheader
       header(i) = ""
    enddo

    ! start putting information relative to this code and run
    call add_card(header)
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
    call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
    call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
    call add_card(header) ! blank line
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"POLCCONV","COSMO"," Coord. convention for polarisation (COSMO/IAU)")
    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
    call add_card(header) ! blank line
    call add_card(header,"POLAR",polarization," Polarisation included (True/False)")

    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "I_"//ttype_,"Stokes I")
    call add_card(header,"TUNIT1", unit_,"Map unit")
    call add_card(header)

    if (polarization) then
       call add_card(header,"TTYPE2", "Q_"//ttype_,"Stokes Q")
       call add_card(header,"TUNIT2", unit_,"Map unit")
       call add_card(header)
       
       call add_card(header,"TTYPE3", "U_"//ttype_,"Stokes U")
       call add_card(header,"TUNIT3", unit_,"Map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Commander Keywords                        ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","Commander is a code for global CMB analysis    ")
    call add_card(header,"COMMENT","developed in collaboration between the University")
    call add_card(header,"COMMENT","of Oslo and Jet Propulsion Laboratory (NASA).  ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    if (present(comptype)) call add_card(header,"COMPTYPE",trim(comptype), "Component type")
    if (present(nu_ref))   call add_card(header,"NU_REF",  nu_ref,         "Reference frequency")
    if (present(spectrumfile)) call add_card(header,"SPECFILE",  trim(spectrumfile), &
         & "Reference spectrum")
    call add_card(header,"COMMENT","-----------------------------------------------")

    call output_map(map, header, "!"//trim(filename))

  end subroutine write_map

  ! Note: Allocates a full-size internal map. Should not be used extensively.
  subroutine udgrade(self, map_out)
    implicit none
    class(comm_map), intent(in)    :: self
    class(comm_map), intent(inout) :: map_out

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:,:) :: m_in, m_out, buffer

    if (self%info%nside == map_out%info%nside) then
       map_out%map = self%map
       return
    end if

    allocate(m_in(0:self%info%npix-1,self%info%nmaps))
    allocate(m_out(0:map_out%info%npix-1,map_out%info%nmaps))
    allocate(buffer(0:map_out%info%npix-1,map_out%info%nmaps))
    m_in                  = 0.d0
    m_in(self%info%pix,:) = self%map
    call udgrade_ring(m_in, self%info%nside, m_out, map_out%info%nside)
    call mpi_allreduce(m_out, buffer, size(m_out), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
    map_out%map = buffer(map_out%info%pix,:)
    deallocate(m_in, m_out, buffer)

  end subroutine udgrade


  subroutine smooth(self, fwhm)
    implicit none
    class(comm_map), intent(inout)    :: self
    real(dp),        intent(in)       :: fwhm

    integer(i4b) :: i, j, l, m
    real(dp)     :: sigma, bl, fact_pol
    real(dp), allocatable, dimension(:,:) :: m_in, m_out

    if (fwhm <= 0.d0) return

    call self%YtW
    sigma    = fwhm * pi/180.d0/60.d0 / sqrt(8.d0*log(2.d0))
    fact_pol = exp(2.d0*sigma**2)
    do i = 0, self%info%nalm-1
       call self%info%i2lm(i,l,m)
       do j = 1, self%info%nmaps
          bl = exp(-0.5*l*(l+1)*sigma**2)
          if (j > 1) bl = bl * fact_pol
          self%alm(i,j) = self%alm(i,j) * bl
       end do
    end do
    call self%Y    

  end subroutine smooth

  !**************************************************
  !                   Utility routines
  !**************************************************

  subroutine alm_equal(self, map)
    implicit none
    class(comm_map), intent(in)    :: self
    class(comm_map), intent(inout) :: map

    integer(i4b) :: i, j, l, m, q

    map%alm = 0.d0
    do i = 0, map%info%nalm-1
       call map%info%i2lm(i,l,m)
       call self%info%lm2i(l,m,j)
       if (j == -1) cycle
       do q = 1, min(self%info%nmaps, map%info%nmaps)
          map%alm(i,q) = self%alm(j,q)
       end do
    end do

  end subroutine alm_equal

  subroutine add_alm(self, alm, info)
    implicit none
    class(comm_map),                       intent(inout) :: self
    real(dp),            dimension(0:,1:), intent(in)    :: alm
    class(comm_mapinfo),                   intent(in)    :: info

    integer(i4b) :: i, j, l, m, q

    do i = 0, info%nalm-1
       call info%i2lm(i,l,m)
       call self%info%lm2i(l,m,j)
       if (j == -1) cycle
       do q = 1, min(self%info%nmaps, info%nmaps)
          self%alm(j,q) = self%alm(j,q) + alm(i,q)
       end do
    end do

  end subroutine add_alm
  
  subroutine lm2i(self, l, m, i)
    implicit none
    class(comm_mapinfo)             :: self
    integer(i4b),       intent(in)  :: l, m
    integer(i4b),       intent(out) :: i

    if (l > self%lmax) then
       i = -1
       return
    end if
    
    if (m == 0) then
       i = self%mind(m) + l
    else
       i = self%mind(abs(m)) + 2*(l-abs(m))
       if (m < 0) i = i+1
    end if
    
  end subroutine lm2i

  subroutine i2lm(self, i, l, m)
    implicit none
    class(comm_mapinfo)             :: self
    integer(i4b),       intent(in)  :: i
    integer(i4b),       intent(out) :: l, m

    if (i > self%nalm) then
       l = -1
       m = -1
       return
    end if
    
    l = self%lm(1,i)
    m = self%lm(2,i)    
    
  end subroutine i2lm
  
  function next(self)
    class(comm_map) :: self
    class(comm_map), pointer :: next
    next => self%nextLink
  end function next

  function prev(self)
    class(comm_map) :: self
    class(comm_map), pointer :: prev
    prev => self%prevLink
  end function prev
  
  subroutine setNext(self,next)
    class(comm_map) :: self
    class(comm_map), pointer :: next
    self%nextLink => next
  end subroutine setNext

  subroutine add(self,link)
    class(comm_map)         :: self
    class(comm_map), target :: link

    class(comm_map), pointer :: c
    
    if (associated(self%nextLink)) then
       c => self%nextLink
       do while (associated(c%nextLink))
          c => c%nextLink
       end do
       !link%prevLink => c
       c%nextLink    => link
    else
       !link%prevLink => self
       self%nextLink => link
    end if

  end subroutine add

  subroutine getSigmaL(self, sigma_l)
    implicit none
    class(comm_map),                   intent(in)  :: self
    real(dp),        dimension(0:,1:), intent(out) :: sigma_l

    integer(i4b) :: l, m, i, j, k, ind, nspec, lmax, nmaps, ierr

    lmax  = self%info%lmax
    nmaps = self%info%nmaps
    nspec = nmaps*(nmaps+1)/2
    sigma_l = 0.d0
    do ind = 0, self%info%nalm-1
       call self%info%i2lm(ind,l,m)
       k   = 1
       do i = 1, nmaps
          do j = i, nmaps
             sigma_l(l,k) = sigma_l(l,k) + self%alm(ind,i)*self%alm(ind,j)
             k            = k+1
          end do
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE, sigma_l, size(sigma_l), MPI_DOUBLE_PRECISION, &
         & MPI_SUM, self%info%comm, ierr)

    do l = 0, lmax
       sigma_l(l,:) = sigma_l(l,:) / real(2*l+1,dp)
    end do

  end subroutine getSigmaL

  subroutine getCrossSigmaL(self, map2, sigma_l)
    implicit none
    class(comm_map),                   intent(in)  :: self, map2
    real(dp),        dimension(0:,1:), intent(out) :: sigma_l

    integer(i4b) :: l, m, i, j, k, ind, ind2, nspec, lmax, nmaps, ierr

    lmax  = self%info%lmax
    nmaps = self%info%nmaps
    nspec = nmaps*(nmaps+1)/2
    sigma_l = 0.d0
    do ind = 0, self%info%nalm-1
       call self%info%i2lm(ind,l,m)
       call map2%info%lm2i(l,m,ind2)
       k   = 1
       do i = 1, nmaps
          do j = i, nmaps
             sigma_l(l,k) = sigma_l(l,k) + self%alm(ind,i)*map2%alm(ind2,j)
             k            = k+1
          end do
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE, sigma_l, size(sigma_l), MPI_DOUBLE_PRECISION, &
         & MPI_SUM, self%info%comm, ierr)

    do l = 0, lmax
       sigma_l(l,:) = sigma_l(l,:) / real(2*l+1,dp)
    end do

  end subroutine getCrossSigmaL
  
end module comm_map_mod
