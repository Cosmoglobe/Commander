module comm_map_mod
  use sharp
  use healpix_types
  use fitstools
  use pix_tools
  use udgrade_nr
  use iso_c_binding, only : c_ptr, c_double
  use head_fits
  implicit none

  include "mpif.h"
      
  private
  public comm_map, comm_mapinfo, map_ptr

  type :: comm_mapinfo
     type(sharp_alm_info)  :: alm_info
     type(sharp_geom_info) :: geom_info
     logical(lgt) :: pol
     integer(i4b) :: comm, myid, nprocs
     integer(i4b) :: nside, npix, nmaps, nring, np, lmax, nm, nalm
     integer(c_int), allocatable, dimension(:)   :: rings
     integer(c_int), allocatable, dimension(:)   :: ms
     integer(c_int), allocatable, dimension(:)   :: mind
     integer(c_int), allocatable, dimension(:,:) :: lm
     integer(c_int), allocatable, dimension(:)   :: pix
     real(c_double), allocatable, dimension(:)   :: W
   contains
     procedure     :: lm2i
     procedure     :: i2lm
     final :: comm_mapinfo_finalize
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
     procedure     :: writeFITS
     procedure     :: readFITS
     procedure     :: dealloc => deallocate_comm_map
     procedure     :: alm_equal
     procedure     :: add_alm
     procedure     :: udgrade

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
    integer(i4b) :: l, m, i, j, iring, np, ind
    integer(i4b), allocatable, dimension(:) :: pixlist

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    allocate(constructor_mapinfo)
    constructor_mapinfo%comm   = comm
    constructor_mapinfo%myid   = myid
    constructor_mapinfo%nprocs = nprocs
    constructor_mapinfo%nside  = nside
    constructor_mapinfo%nmaps  = nmaps
    constructor_mapinfo%lmax   = lmax
    constructor_mapinfo%pol    = pol
    constructor_mapinfo%npix   = 12*nside**2

    ! Select rings and pixels
    allocate(pixlist(0:4*nside-1))
    constructor_mapinfo%nring = 0
    constructor_mapinfo%np    = 0
    do i = 1+myid, 4*nside-1, nprocs
       call in_ring(nside, i, 0.d0, pi, pixlist, np)
       constructor_mapinfo%nring = constructor_mapinfo%nring + 1
       constructor_mapinfo%np    = constructor_mapinfo%np    + np
    end do
    allocate(constructor_mapinfo%rings(constructor_mapinfo%nring))
    allocate(constructor_mapinfo%pix(constructor_mapinfo%np))
    j = 1
    do i = 1, constructor_mapinfo%nring
       constructor_mapinfo%rings(i) = 1 + myid + (i-1)*nprocs
       call in_ring(nside, constructor_mapinfo%rings(i), 0.d0, pi, pixlist, np)
       constructor_mapinfo%pix(j:j+np-1) = pixlist(0:np-1)
       j = j + np
    end do
    deallocate(pixlist)

    ! Select m's
    constructor_mapinfo%nm = 0
    do m = myid, lmax, nprocs
       constructor_mapinfo%nm   = constructor_mapinfo%nm   + 1
       if (m == 0) then
          constructor_mapinfo%nalm = constructor_mapinfo%nalm + lmax+1
       else
          constructor_mapinfo%nalm = constructor_mapinfo%nalm + 2*(lmax-m+1)
       end if
    end do
    allocate(constructor_mapinfo%ms(constructor_mapinfo%nm))
    allocate(constructor_mapinfo%mind(0:constructor_mapinfo%lmax))
    allocate(constructor_mapinfo%lm(2,0:constructor_mapinfo%nalm-1))
    ind = 0
    constructor_mapinfo%mind = -1
    do i = 1, constructor_mapinfo%nm
       m                           = myid + (i-1)*nprocs
       constructor_mapinfo%ms(i)   = m
       constructor_mapinfo%mind(m) = ind
       if (m == 0) then
          do l = m, lmax
             constructor_mapinfo%lm(:,ind) = [l,m]
             ind                           = ind+1
          end do
       else
          do l = m, lmax
             constructor_mapinfo%lm(:,ind) = [l,+m]
             ind                           = ind+1
             constructor_mapinfo%lm(:,ind) = [l,-m]
             ind                           = ind+1
          end do
       end if
    end do
    
    ! Read ring weights
    allocate(constructor_mapinfo%W(constructor_mapinfo%nring))
    constructor_mapinfo%W = 1.d0

    ! Create SHARP info structures
    call sharp_make_mmajor_real_packed_alm_info(lmax, ms=constructor_mapinfo%ms, &
         & alm_info=constructor_mapinfo%alm_info)
    call sharp_make_healpix_geom_info(nside, rings=constructor_mapinfo%rings, &
         & geom_info=constructor_mapinfo%geom_info)
    
  end function constructor_mapinfo

  function constructor_map(info, filename)
    implicit none
    class(comm_mapinfo),           intent(in), target   :: info
    character(len=*),              intent(in), optional :: filename
    class(comm_map),     pointer                        :: constructor_map

    allocate(constructor_map)
    constructor_map%info => info
    allocate(constructor_map%map(0:info%np-1,info%nmaps))
    allocate(constructor_map%alm(0:info%nalm-1,info%nmaps))

    if (present(filename)) then
       call constructor_map%readFITS(filename)
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

  subroutine deallocate_comm_map(self)
    implicit none

    class(comm_map), intent(inout) :: self

    if (allocated(self%map)) deallocate(self%map)
    if (allocated(self%alm)) deallocate(self%alm)
    nullify(self%info)

  end subroutine deallocate_comm_map

  subroutine comm_mapinfo_finalize(self)
    implicit none

    type(comm_mapinfo) :: self
    
    deallocate(self%rings, self%ms, self%mind, self%lm, self%pix, self%W)
    call sharp_destroy_alm_info(self%alm_info)
    call sharp_destroy_geom_info(self%geom_info)
    
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
            & self%map(:,1:1), self%info%geom_info, comm=self%info%comm)
       if (self%info%nmaps == 3) then
!          call sharp_execute(SHARP_Y, 2, 1, self%alm(:,2:3), self%info%alm_info, &
!               & self%map(:,2:3), self%info%geom_info, comm=self%info%comm)
          self%map(:,2:3) = 0.d0
       end if
    else
       call sharp_execute(SHARP_Y, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info, comm=self%info%comm)       
    end if
    
  end subroutine exec_sharp_Y

  subroutine exec_sharp_Yt(self)
    implicit none

    class(comm_map), intent(inout) :: self

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_Yt, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info, comm=self%info%comm)
       if (self%info%nmaps == 3) then
          !call sharp_execute(SHARP_Yt, 2, 1, self%alm(:,2:3), self%info%alm_info, &
          !     & self%map(:,2:3), self%info%geom_info, comm=self%info%comm)
          self%alm(:,2:3) = 0.d0
       end if
       
    else
       call sharp_execute(SHARP_Yt, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info, comm=self%info%comm)       
    end if
    
  end subroutine exec_sharp_Yt

  subroutine exec_sharp_YtW(self)
    implicit none

    class(comm_map), intent(inout) :: self

    if (.not. allocated(self%alm)) allocate(self%alm(0:self%info%nalm-1,self%info%nmaps))
    if (self%info%pol) then
       call sharp_execute(SHARP_YtW, 0, 1, self%alm(:,1:1), self%info%alm_info, &
            & self%map(:,1:1), self%info%geom_info, comm=self%info%comm)
       if (self%info%nmaps == 3) then
!          call sharp_execute(SHARP_YtW, 2, 1, self%alm(:,2:3), self%info%alm_info, &
!               & self%map(:,2:3), self%info%geom_info, comm=self%info%comm)
          self%alm(:,2:3) = 0.d0
       end if
    else
       call sharp_execute(SHARP_YtW, 0, self%info%nmaps, self%alm, self%info%alm_info, &
            & self%map, self%info%geom_info, comm=self%info%comm)
    end if
    
  end subroutine exec_sharp_YtW
  
  !**************************************************
  !                   IO routines
  !**************************************************

  subroutine writeFITS(self, filename, comptype, nu_ref, unit, ttype, spectrumfile)
    implicit none

    class(comm_map),  intent(in) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),         intent(in), optional :: nu_ref

    integer(i4b) :: i, nmaps, npix, np, ierr
    real(dp),     allocatable, dimension(:,:) :: map, buffer
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat
    
    ! Only the root actually writes to disk; data are distributed via MPI
    npix  = self%info%npix
    nmaps = self%info%nmaps
    if (self%info%myid == 0) then

       ! Distribute to other nodes
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
       call write_map(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
       deallocate(p, map)

    else
       call mpi_send(self%info%np,  1,              MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%info%pix, self%info%np,   MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%map,      size(self%map), MPI_DOUBLE_PRECISION, 0, 98, &
            & self%info%comm, ierr)
    end if
    
  end subroutine writeFITS
  
  subroutine readFITS(self, filename)
    implicit none

    class(comm_map),  intent(inout) :: self
    character(len=*), intent(in)    :: filename

    integer(i4b) :: i, np, npix, ordering, nside, nmaps, ierr
    real(dp),     allocatable, dimension(:,:) :: map
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat

    ! Check file consistency 
    npix = getsize_fits(trim(filename), ordering=ordering, nside=nside, nmaps=nmaps)
    if (nmaps < self%info%nmaps) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nmaps in ' // trim(filename)
       call mpi_finalize(ierr)
       stop
    end if
    if (nside /= self%info%nside) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nside in ' // trim(filename)
       call mpi_finalize(ierr)
       stop
    end if

    ! Only the root actually reads from disk; data are distributed via MPI
    if (self%info%myid == 0) then

       ! Read map and convert to RING format if necessary
       allocate(map(0:self%info%npix-1,self%info%nmaps))
       call input_map(filename, map, self%info%npix, self%info%nmaps)
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
          call mpi_send(map(p(1:np),:), np*self%info%nmaps, MPI_DOUBLE_PRECISION, i, 98, &
               & self%info%comm, ierr)
       end do
       deallocate(p, map)
    else
       call mpi_send(self%info%np,               1, MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%info%pix, self%info%np,   MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_recv(self%map,      size(self%map), MPI_DOUBLE_PRECISION, 0, 98, &
            &  self%info%comm, mpistat, ierr)
    end if
    
  end subroutine readFITS

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
    real(dp), allocatable, dimension(:,:) :: m_in, m_out

!    if (self%info%nside == map_out%info%nside) then
!       map_out%map = self%map
!       return
!    end if

    allocate(m_in(0:self%info%npix-1,self%info%nmaps))
    allocate(m_out(0:map_out%info%npix-1,map_out%info%nmaps))
    m_in                  = 0.d0
    m_in(self%info%pix,:) = self%map
    call udgrade_ring(m_in, self%info%nside, m_out, map_out%info%nside)
    call mpi_allreduce(MPI_IN_PLACE, m_out, size(m_out), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
    map_out%map = m_out(map_out%info%pix,:)
    deallocate(m_in, m_out)

  end subroutine udgrade

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
  
end module comm_map_mod
