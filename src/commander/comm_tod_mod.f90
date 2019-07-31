module comm_tod_mod
  use comm_utils
  use comm_param_mod
  use comm_hdf_mod
  use comm_map_mod
  use comm_fft_mod
  use comm_huffman_mod
  use comm_conviqt_mod
  USE ISO_C_BINDING
  implicit none

  private
  public comm_tod, initialize_tod_mod

  type :: comm_detscan
     character(len=10) :: label                           ! Detector label
     real(dp)          :: gain, gain_sigma                ! Gain; assumed constant over scan
     real(dp)          :: sigma0, alpha, fknee            ! Noise parameters
     real(dp)          :: samprate                        ! Sample rate in Hz
     real(dp)          :: polang                          ! Detector polarization angle
     real(dp)          :: mbang                           ! Main beams angle
     real(dp)          :: mono                            ! Monopole
     real(dp)          :: chisq
     real(dp)          :: chisq_prop
     real(dp)          :: chisq_masked
     logical(lgt)      :: accept
     integer(i4b)      :: nside                           ! Nside for pixelized pointing
     real(sp),     allocatable, dimension(:)   :: tod     ! Detector values in time domain, (ntod)     
     byte,         allocatable, dimension(:)   :: flag    ! Compressed detector flag; 0 is accepted, /= 0 is rejected
     byte,         allocatable, dimension(:)   :: pix     ! Compressed pixelized pointing, (ntod,nhorn)
     byte,         allocatable, dimension(:)   :: psi     ! Compressed polarization angle, (ntod,nhorn)
     real(dp),     allocatable, dimension(:)   :: N_psd   ! Noise power spectrum density; in uncalibrated units
     real(dp),     allocatable, dimension(:)   :: nu_psd  ! Noise power spectrum bins; in Hz
  end type comm_detscan

  type :: comm_scan
     integer(i4b)   :: ntod                                        ! Number of time samples
     real(dp)       :: proctime    = 0.d0                          ! Processing time in seconds
     real(dp)       :: n_proctime  = 0                             ! Number of completed loops
     real(dp)       :: v_sun(3)                                    ! Observatory velocity relative to Sun in km/s
     real(dp)       :: t0                                          ! MJD for first sample; not currently used, so set to zero for now
     type(huffcode) :: hkey                                        ! Huffman decompression key
     class(comm_detscan), allocatable, dimension(:)     :: d       ! Array of all detectors
  end type comm_scan

  type, abstract :: comm_tod 
     character(len=512) :: freq
     character(len=512) :: filelist
     character(len=512) :: procmaskf1
     character(len=512) :: procmaskf2
     character(len=512) :: instfile
     integer(i4b) :: comm, myid, numprocs                         ! MPI parameters
     integer(i4b) :: comm_shared, myid_shared                     ! MPI parameters
     integer(i4b) :: comm_inter, myid_inter                     ! MPI parameters
     integer(i4b) :: nmaps                                        ! Number of Stokes parameters
     integer(i4b) :: ndet                                         ! Number of active detectors
     integer(i4b) :: nhorn                                        ! Number of horns
     integer(i4b) :: nscan, nscan_tot                              ! Number of scans
     integer(i4b) :: npsi                                         ! Number of discretized psi steps
     real(dp)     :: samprate                                     ! Sample rate in Hz
     logical(lgt) :: output_all                                   ! Oitput all samples
     logical(lgt) :: init_from_HDF                                   ! Oitput all samples
     integer(i4b),       allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),           allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     real(sp),           allocatable, dimension(:)     :: sin2psi  ! Lookup table of sin(2psi) 
     real(sp),           allocatable, dimension(:)     :: cos2psi  ! Lookup table of cos(2psi) 
     real(sp),           allocatable, dimension(:)     :: psi     ! Lookup table of psi
     type(comm_scan),    allocatable, dimension(:)     :: scans    ! Array of all scans
     integer(i4b),       allocatable, dimension(:)     :: scanid   ! List of scan IDs
     integer(i4b),       allocatable, dimension(:)     :: nscanprproc   ! List of scan IDs
     character(len=512), allocatable, dimension(:)     :: hdfname  ! List of HDF filenames for each ID
     character(len=512), allocatable, dimension(:)     :: label    ! Detector labels
     class(comm_map), pointer                          :: procmask ! Mask for gain and n_corr
     class(comm_map), pointer                          :: procmask2 ! Mask for gain and n_corr
     class(comm_mapinfo), pointer                      :: info     ! Map definition
     class(comm_mapinfo), pointer                      :: slinfo   ! Sidelobe map info
     class(map_ptr),     allocatable, dimension(:)     :: slbeam   ! Sidelobe beam data (ndet)
     class(conviqt_ptr), allocatable, dimension(:)     :: slconv   ! SL-convolved maps (ndet)
     real(dp),           allocatable, dimension(:,:)   :: bp_delta  ! Bandpass parameters (ndet, ndelta)
   contains
     procedure                        :: read_tod
     procedure                        :: get_scan_ids
     procedure                        :: dumpToHDF
     procedure                        :: initHDF
     procedure(process_tod), deferred :: process_tod
  end type comm_tod

  abstract interface
     subroutine process_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
       import i4b, comm_tod, comm_map, map_ptr, dp, planck_rng
       implicit none
       class(comm_tod),                     intent(inout) :: self
       character(len=*),                    intent(in)    :: chaindir
       integer(i4b),                        intent(in)    :: chain, iter
       type(planck_rng),                    intent(inout) :: handle
       type(map_ptr),     dimension(:,:),   intent(inout) :: map_in            
       real(dp),          dimension(:,:,:), intent(inout) :: delta
       class(comm_map),                     intent(inout) :: map_out, rms_out  
     end subroutine process_tod
  end interface


contains

  subroutine initialize_tod_mod(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar

    logical(lgt), save :: initialized = .false.

    if (initialized) return

    call initialize_fft_mod(cpar)

  end subroutine initialize_tod_mod


  !**************************************************
  !             Utility routines
  !**************************************************

  subroutine read_tod(self, detlabels)
    implicit none
    class(comm_tod),                intent(inout)  :: self
    character(len=*), dimension(:), intent(in)     :: detlabels

    integer(i4b) :: i, j, det, ierr, npsi, nside
    real(dp)     :: t1, t2, psi, fsamp
    real(dp), allocatable, dimension(:) :: mbang, polang

    if (self%nscan ==0) then
       call mpi_barrier(self%comm, ierr)
       return
    end if

    call wall_time(t1)
    allocate(self%scans(self%nscan), mbang(self%ndet), polang(self%ndet))
    do i = 1, self%nscan
       call read_hdf_scan(self%scans(i), self%myid, self%hdfname(i), self%scanid(i), self%ndet, &
            & detlabels, i>1, npsi, nside, fsamp, mbang, polang)
       if (i == 1) then
          self%npsi     = npsi
          self%samprate = fsamp
       end if
       do det = 1, self%ndet
          self%scans(i)%d(det)%accept = all(self%scans(i)%d(det)%tod==self%scans(i)%d(det)%tod)
          if (.not. self%scans(i)%d(det)%accept) then
             write(*,fmt='(a,i8,a,i3, i10)') 'Input TOD contain NaN -- scan =', &
                  & self%scanid(i), ', det =', det, count(self%scans(i)%d(det)%tod/=self%scans(i)%d(det)%tod)
             write(*,fmt='(a,a)') '    filename = ', &
                  & trim(self%hdfname(i))
          end if
          self%scans(i)%d(1)%mbang = -22.46 * DEG2RAD 
          self%scans(i)%d(2)%mbang = -22.46 * DEG2RAD 
          self%scans(i)%d(3)%mbang =  22.45 * DEG2RAD 
          self%scans(i)%d(4)%mbang =  22.45 * DEG2RAD 
          !self%scans(i)%d(:)%mbang =   0.d0
       end do
    end do

    ! Precompute trigonometric functions
    allocate(self%sin2psi(0:self%npsi-1), self%cos2psi(0:self%npsi-1))
    allocate(self%psi(0:self%npsi-1))
    do i = 0, self%npsi-1
       psi             = (i+0.5d0)*2.d0*pi/real(self%npsi,dp)
       self%psi(i)     = psi
       self%sin2psi(i) = sin(2.d0*psi)
       self%cos2psi(i) = cos(2.d0*psi)
    end do

    deallocate(polang, mbang)
    call mpi_barrier(self%comm, ierr)
    call wall_time(t2)
    if (self%myid == 0) write(*,fmt='(a,i4,a,i6,a,f8.1,a)') &
         & '    Myid = ', self%myid, ' -- nscan = ', self%nscan, &
         & ', TOD IO time = ', t2-t1, ' sec'

  end subroutine read_tod

  subroutine read_hdf_scan(self, myid, filename, scan, ndet, detlabels, use_input, npsi, nside, fsamp, mbang, polang)
    implicit none
    integer(i4b),                   intent(in)    :: myid
    character(len=*),               intent(in)    :: filename
    integer(i4b),                   intent(in)    :: scan, ndet
    class(comm_scan),               intent(inout) :: self
    logical(lgt),                   intent(in)    :: use_input
    integer(i4b),                   intent(inout) :: npsi, nside
    real(dp),                       intent(inout) :: fsamp
    real(dp), dimension(:),         intent(inout) :: mbang, polang

    character(len=*), dimension(:), intent(in)     :: detlabels

    integer(i4b)       :: i, j, k, n, m, nhorn, ext(1), ierr
    real(dp)           :: psi, t1, t2, t3, t4, t_tot(6)
    character(len=6)   :: slabel
    character(len=32)   :: out
    character(len=8)   :: out2
    character(len=128) :: field
    integer(hid_t)     :: nfield, err, obj_type
    type(hdf_file)     :: file
    integer(i4b), allocatable, dimension(:) :: hsymb
    real(dp),     allocatable, dimension(:) :: buffer_sp
    integer(i4b), allocatable, dimension(:) :: htree

    call wall_time(t3)
    t_tot = 0.d0

    call wall_time(t1)
    call int2string(scan, slabel)

    call open_hdf_file(filename, file, "r")

    ! Find array sizes
    call get_size_hdf(file, slabel // "/" // trim(detlabels(1)) // "/tod", ext)
    !nhorn     = ext(1)
    n         = ext(1)
    m         = get_closest_fft_magic_number(n)
    self%ntod = m

    ! Read common scan data
    call read_hdf(file, slabel // "/common/vsun",  self%v_sun)
    !call read_hdf(file, slabel // "/common/time",  self%t0)
    self%t0 = 0.d0 ! Not currently used, set to zero for now
    if (.not. use_input) call read_hdf(file, slabel // "/common/npsi",  npsi)
    call wall_time(t2)
    t_tot(1) = t2-t1

    ! Read detector scans
    allocate(self%d(ndet), buffer_sp(n))
    do i = 1, ndet
       call wall_time(t1)
       field = detlabels(i)
       call wall_time(t2)
       t_tot(2) = t_tot(2) + t2-t1
       call wall_time(t1)
       allocate(self%d(i)%tod(m))
       self%d(i)%label = trim(field)
       call read_hdf(file, slabel // "/" // trim(field) // "/gain",   self%d(i)%gain)
       call read_hdf(file, slabel // "/" // trim(field) // "/sigma0", self%d(i)%sigma0)
       call read_hdf(file, slabel // "/" // trim(field) // "/alpha",  self%d(i)%alpha)
       call read_hdf(file, slabel // "/" // trim(field) // "/fknee",  self%d(i)%fknee)
       if (.not. use_input) then
          call read_hdf(file, slabel // "/" // trim(field) // "/polang", polang(i))
!          call read_hdf(file, slabel // "/" // trim(field) // "/mbang",  mbang(i))
          call read_hdf(file, slabel // "/common/nside",                 nside)
          call read_hdf(file, slabel // "/common/fsamp",                 fsamp)
       end if
       self%d(i)%polang   = polang(i)
       self%d(i)%mbang    = mbang(i) * DEG2RAD 
       self%d(i)%nside    = nside
       self%d(i)%samprate = fsamp
       call wall_time(t2)
       t_tot(3) = t_tot(3) + t2-t1
       call wall_time(t1)
       call read_hdf(file, slabel // "/" // trim(field) // "/tod",    buffer_sp)
       self%d(i)%tod = buffer_sp(1:m)
       call wall_time(t2)
       t_tot(4) = t_tot(4) + t2-t1
   
       ! Read Huffman coded data arrays
       call wall_time(t1)
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/pix",  self%d(i)%pix)
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/psi",  self%d(i)%psi)
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/flag", self%d(i)%flag)
       call wall_time(t2)
       t_tot(5) = t_tot(5) + t2-t1
    end do
    deallocate(buffer_sp)

    ! Initialize Huffman key
    call wall_time(t1)
    call read_alloc_hdf(file, slabel // "/common/huffsymb", hsymb)
    call read_alloc_hdf(file, slabel // "/common/hufftree", htree)
    call hufmak_precomp(hsymb,htree,self%hkey)
    !call hufmak(hsymb,hfreq,self%hkey)
    deallocate(hsymb, htree)
    call wall_time(t2)
    t_tot(6) = t_tot(6) + t2-t1

    call close_hdf_file(file)

    call wall_time(t4)

    if (myid == 0) then
       write(*,*)
       write(*,*) '  IO init   = ', t_tot(1)
       write(*,*) '  IO field  = ', t_tot(2)
       write(*,*) '  IO scalar = ', t_tot(3)
       write(*,*) '  IO tod    = ', t_tot(4)
       write(*,*) '  IO comp   = ', t_tot(5)
       write(*,*) '  IO huff   = ', t_tot(6)
       write(*,*) '  IO total  = ', t4-t3
    end if

  end subroutine read_hdf_scan


  subroutine get_scan_ids(self, filelist)
    implicit none
    class(comm_tod),   intent(inout) :: self    
    character(len=*),  intent(in)    :: filelist

    integer(i4b)       :: unit, j, k, np, ind(1), i, n, n_tot, ierr
    real(dp),           allocatable, dimension(:) :: weight
    integer(i4b),       allocatable, dimension(:) :: scanid, id
    integer(i4b),       allocatable, dimension(:) :: proc
    real(dp),           allocatable, dimension(:) :: pweight
    character(len=512), allocatable, dimension(:) :: filename

    np = self%numprocs
    if (self%myid == 0) then
       unit = getlun()
       open(unit, file=trim(filelist))
       read(unit,*) n_tot
       allocate(id(n_tot), filename(n_tot), scanid(n_tot), weight(n_tot), proc(n_tot), pweight(0:np-1))
       do i = 1, n_tot
          id(i) = i
          read(unit,*) scanid(i), filename(i), weight(i)
       end do

       ! Sort according to weight
       pweight = 0.d0
       call QuickSort(id, weight)
       do i = n_tot, 1, -1
          ind             = minloc(pweight)-1
          proc(id(i))     = ind(1)
          pweight(ind(1)) = pweight(ind(1)) + weight(i)
       end do
       deallocate(id, pweight, weight)
    end if

    call mpi_bcast(n_tot, 1,  MPI_INTEGER, 0, self%comm, ierr)
    if (self%myid /= 0) then
       allocate(filename(n_tot), scanid(n_tot), proc(n_tot))
    end if
    call mpi_bcast(filename, 512*n_tot,  MPI_CHARACTER, 0, self%comm, ierr)
    call mpi_bcast(scanid,       n_tot,  MPI_INTEGER,   0, self%comm, ierr)
    call mpi_bcast(proc,         n_tot,  MPI_INTEGER,   0, self%comm, ierr)

    self%nscan     = count(proc == self%myid)
    self%nscan_tot = maxval(scanid)
    allocate(self%scanid(self%nscan), self%hdfname(self%nscan))
    j = 1
    do i = 1, n_tot
       if (proc(i) == self%myid) then
          self%scanid(j)  = scanid(i)
          self%hdfname(j) = filename(i)
          j               = j+1
       end if
    end do

    if (self%myid == 0) then
       allocate(self%nscanprproc(0:self%numprocs-1))
       do i = 0, self%numprocs-1
          self%nscanprproc(i) = count(proc == i)
       end do
    end if

    deallocate(filename, scanid, proc) 

  end subroutine get_scan_ids


  subroutine dumpToHDF(self, chainfile, iter)
    implicit none
    class(comm_tod),                   intent(in)    :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile

    integer(i4b)       :: i, j, k, npar, ierr
    character(len=6)   :: itext
    character(len=512) :: path
    real(dp), allocatable, dimension(:,:,:) :: output

    npar = 7
    allocate(output(self%nscan_tot,self%ndet,npar))

    ! Collect all parameters 
    output = 0.d0
    do j = 1, self%ndet
       do i = 1, self%nscan
          k             = self%scanid(i)
          output(k,j,1) = self%scans(i)%d(j)%gain
          output(k,j,2) = self%scans(i)%d(j)%sigma0
          output(k,j,3) = self%scans(i)%d(j)%alpha
          output(k,j,4) = self%scans(i)%d(j)%fknee
          output(k,j,5) = self%scans(i)%d(j)%polang
          output(k,j,6) = self%scans(i)%d(j)%mono
          output(k,j,7) = merge(1.d0,0.d0,self%scans(i)%d(j)%accept)
       end do
    end do

    if (self%myid == 0) then
       call mpi_reduce(mpi_in_place, output, size(output), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    else
       call mpi_reduce(output,       output, size(output), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    end if

    if (self%myid == 0) then
       call int2string(iter, itext)
       path = trim(adjustl(itext))//'/tod/'//trim(adjustl(self%freq))//'/'
       call create_hdf_group(chainfile, trim(adjustl(path)))
       call write_hdf(chainfile, trim(adjustl(path))//'gain',   output(:,:,1))
       call write_hdf(chainfile, trim(adjustl(path))//'sigma0', output(:,:,2))
       call write_hdf(chainfile, trim(adjustl(path))//'alpha',  output(:,:,3))
       call write_hdf(chainfile, trim(adjustl(path))//'fknee',  output(:,:,4))
       call write_hdf(chainfile, trim(adjustl(path))//'polang', output(:,:,5))
       call write_hdf(chainfile, trim(adjustl(path))//'mono',   output(:,:,6))
       call write_hdf(chainfile, trim(adjustl(path))//'accept', output(:,:,7))
       call write_hdf(chainfile, trim(adjustl(path))//'bp_delta', self%bp_delta)
    end if

    deallocate(output)
 
  end subroutine dumpToHDF

  subroutine initHDF(self, chainfile, iter)
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile

    integer(i4b)       :: i, j, k, npar, ierr
    character(len=6)   :: itext
    character(len=512) :: path
    real(dp), allocatable, dimension(:,:,:) :: output

    npar = 7
    allocate(output(self%nscan_tot,self%ndet,npar))

    if (self%myid == 0) then
       call int2string(iter, itext)
       path = trim(adjustl(itext))//'/tod/'//trim(adjustl(self%freq))//'/'
       call read_hdf(chainfile, trim(adjustl(path))//'gain',     output(:,:,1))
       call read_hdf(chainfile, trim(adjustl(path))//'sigma0',   output(:,:,2))
       call read_hdf(chainfile, trim(adjustl(path))//'alpha',    output(:,:,3))
       call read_hdf(chainfile, trim(adjustl(path))//'fknee',    output(:,:,4))
       call read_hdf(chainfile, trim(adjustl(path))//'polang',   output(:,:,5))
       call read_hdf(chainfile, trim(adjustl(path))//'mono',     output(:,:,6))
       call read_hdf(chainfile, trim(adjustl(path))//'accept',   output(:,:,7))
       call read_hdf(chainfile, trim(adjustl(path))//'bp_delta', self%bp_delta)
    end if

    call mpi_bcast(output, size(output), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%bp_delta, size(self%bp_delta), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)

    do j = 1, self%ndet
       do i = 1, self%nscan
          k             = self%scanid(i)
          self%scans(i)%d(j)%gain   = output(k,j,1)
          self%scans(i)%d(j)%sigma0 = output(k,j,2)
          self%scans(i)%d(j)%alpha  = output(k,j,3)
          self%scans(i)%d(j)%fknee  = output(k,j,4)
          self%scans(i)%d(j)%polang = output(k,j,5)
          self%scans(i)%d(j)%mono   = output(k,j,6)
          self%scans(i)%d(j)%accept = output(k,j,7) == 1.d0
       end do
    end do

    deallocate(output)
 
  end subroutine initHDF


end module comm_tod_mod
