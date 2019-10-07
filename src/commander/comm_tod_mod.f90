module comm_tod_mod
  use comm_utils
  use comm_param_mod
  use comm_hdf_mod
  use comm_map_mod
  use comm_fft_mod
  use comm_huffman_mod
  use comm_conviqt_mod
  use comm_zodi_mod
  USE ISO_C_BINDING
  implicit none

  private
  public comm_tod, initialize_tod_mod

  type :: comm_detscan
     character(len=10) :: label                           ! Detector label
     real(dp)          :: gain, gain_sigma                ! Gain; assumed constant over scan
     real(dp)          :: sigma0, alpha, fknee            ! Noise parameters
     real(dp)          :: chisq
     real(dp)          :: chisq_prop
     real(dp)          :: chisq_masked
     logical(lgt)      :: accept
     real(sp),     allocatable, dimension(:)  :: tod        ! Detector values in time domain, (ntod)     
     byte,         allocatable, dimension(:)  :: flag       ! Compressed detector flag; 0 is accepted, /= 0 is rejected
     byte,         allocatable, dimension(:)  :: pix        ! Compressed pixelized pointing, (ntod,nhorn)
     byte,         allocatable, dimension(:)  :: psi        ! Compressed polarization angle, (ntod,nhorn)
     real(dp),     allocatable, dimension(:)  :: log_n_psd  ! Noise power spectrum density; in uncalibrated units
     real(dp),     allocatable, dimension(:)  :: log_n_psd2 ! Second derivative (for spline)
     real(dp),     allocatable, dimension(:)  :: log_nu     ! Noise power spectrum bins; in Hz
  end type comm_detscan

  type :: comm_scan
     integer(i4b)   :: ntod                                        ! Number of time samples
     real(dp)       :: proctime    = 0.d0                          ! Processing time in seconds
     real(dp)       :: n_proctime  = 0                             ! Number of completed loops
     real(dp)       :: v_sun(3)                                    ! Observatory velocity relative to Sun in km/s
     real(dp)       :: t0                                          ! MJD for first sample
     type(huffcode) :: hkey                                        ! Huffman decompression key
     class(comm_detscan), allocatable, dimension(:)     :: d       ! Array of all detectors
  end type comm_scan

  type, abstract :: comm_tod 
     character(len=512) :: freq
     character(len=512) :: filelist
     character(len=512) :: procmaskf1
     character(len=512) :: procmaskf2
     character(len=512) :: instfile
     character(len=512) :: operation
     character(len=512) :: outdir
     logical(lgt) :: first_call
     integer(i4b) :: comm, myid, numprocs                         ! MPI parameters
     integer(i4b) :: comm_shared, myid_shared, numprocs_shared    ! MPI parameters
     integer(i4b) :: comm_inter, myid_inter                       ! MPI parameters
     integer(i4b) :: nmaps                                        ! Number of Stokes parameters
     integer(i4b) :: ndet                                         ! Number of active detectors
     integer(i4b) :: nhorn                                        ! Number of horns
     integer(i4b) :: nscan, nscan_tot                              ! Number of scans
     integer(i4b) :: first_scan, last_scan
     integer(i4b) :: npsi                                         ! Number of discretized psi steps
     integer(i4b) :: flag0

     real(dp)     :: samprate                                     ! Sample rate in Hz
     real(dp), allocatable, dimension(:)     :: polang                                      ! Detector polarization angle
     real(dp), allocatable, dimension(:)     :: mbang                                       ! Main beams angle
     real(dp), allocatable, dimension(:)     :: mono                                        ! Monopole
     real(dp), allocatable, dimension(:)     :: fwhm, elip, psi_ell                         ! Beam parameter
     real(dp), allocatable, dimension(:,:,:) :: prop_bp         ! proposal matrix, L(ndet,ndet,ndelta),  for bandpass sampler
     real(dp), allocatable, dimension(:)     :: prop_bp_mean    ! proposal matrix, sigma(ndelta), for mean
     integer(i4b)      :: nside                           ! Nside for pixelized pointing
     integer(i4b)      :: nobs                            ! Number of observed pixeld for this core
     integer(i4b) :: output_n_maps                                ! Output n_maps
     logical(lgt) :: init_from_HDF                                   ! Read from HDF file
     integer(i4b) :: output_4D_map                                ! Output 4D maps
     logical(lgt) :: subtract_zodi                                ! Subtract zodical light
     integer(i4b),       allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),           allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     real(dp),           allocatable, dimension(:)     :: sin2psi  ! Lookup table of sin(2psi) 
     real(dp),           allocatable, dimension(:)     :: cos2psi  ! Lookup table of cos(2psi) 
     real(sp),           allocatable, dimension(:)     :: psi      ! Lookup table of psi
     real(dp),           allocatable, dimension(:,:)   :: pix2vec  ! Lookup table of pix2vec
     real(dp),           allocatable, dimension(:,:)   :: L_prop_mono  ! Proposal matrix for monopole sampling
     type(comm_scan),    allocatable, dimension(:)     :: scans    ! Array of all scans
     integer(i4b),       allocatable, dimension(:)     :: scanid   ! List of scan IDs
     integer(i4b),       allocatable, dimension(:)     :: nscanprproc   ! List of scan IDs
     integer(i4b),       allocatable, dimension(:)     :: partner ! Partner detector; for symmetrizing flags     
     integer(i4b),       allocatable, dimension(:)     :: horn_id  ! Internal horn number per detector
     character(len=512), allocatable, dimension(:)     :: hdfname  ! List of HDF filenames for each ID
     character(len=512), allocatable, dimension(:)     :: label    ! Detector labels
     class(comm_map), pointer                          :: procmask ! Mask for gain and n_corr
     class(comm_map), pointer                          :: procmask2 ! Mask for gain and n_corr
     class(comm_mapinfo), pointer                      :: info     ! Map definition
     class(comm_mapinfo), pointer                      :: slinfo   ! Sidelobe map info
     class(map_ptr),     allocatable, dimension(:)     :: slbeam   ! Sidelobe beam data (ndet)
     class(conviqt_ptr), allocatable, dimension(:)     :: slconv   ! SL-convolved maps (ndet)
     real(dp),           allocatable, dimension(:,:)   :: bp_delta  ! Bandpass parameters (0:ndet, npar)
     integer(i4b),       allocatable, dimension(:)     :: pix2ind, ind2pix
   contains
     procedure                        :: read_tod
     procedure                        :: get_scan_ids
     procedure                        :: dumpToHDF
     procedure                        :: initHDF
     procedure                        :: get_det_id
     procedure                        :: initialize_bp_covar
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
       class(comm_map),                     intent(inout) :: map_out
       class(comm_map),                     intent(inout) :: rms_out  
     end subroutine process_tod
  end interface


contains

  subroutine initialize_tod_mod(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar

    logical(lgt), save :: initialized = .false.

    if (initialized) return

    call initialize_fft_mod(cpar)
    if (cpar%include_tod_zodi) call initialize_zodi_mod(cpar)

  end subroutine initialize_tod_mod


  !**************************************************
  !             Utility routines
  !**************************************************

  subroutine read_tod(self, detlabels)
    implicit none
    class(comm_tod),                intent(inout)  :: self
    character(len=*), dimension(:), intent(in)     :: detlabels

    integer(i4b) :: i, j, n, det, ierr, npsi, nside, ndet_tot
    real(dp)     :: t1, t2, psi, fsamp
    type(hdf_file)     :: file

    real(dp), dimension(:), allocatable           :: mbang_buf, polang_buf
    character(len=1024)                           :: det_buf
    character(len=128), dimension(:), allocatable :: dets


    ! Read common fields
    allocate(self%polang(self%ndet), self%mbang(self%ndet), self%mono(self%ndet))
    self%mono = 0.d0
    if (self%myid == 0) then
       call open_hdf_file(self%hdfname(1), file, "r")


       !TODO: figure out how to make this work
       call read_hdf_string2(file, "/common/det",    det_buf, n)
       !call read_hdf(file, "/common/det",    det_buf)
       !write(det_buf, *) "27M, 27S, 28M, 28S"
       !write(det_buf, *) "18M, 18S, 19M, 19S, 20M, 20S, 21M, 21S, 22M, 22S, 23M, 23S"
       ndet_tot = num_tokens(det_buf(1:n), ",")
       allocate(polang_buf(ndet_tot), mbang_buf(ndet_tot), dets(ndet_tot))
       call get_tokens(trim(adjustl(det_buf(1:n))), ',', dets)
!!$       do i = 1, ndet_tot
!!$          write(*,*) i, trim(adjustl(dets(i)))
!!$       end do
       !write(*,*) ndet_tot
       call read_hdf(file, "common/nside",  self%nside)
       call read_hdf(file, "common/npsi",   self%npsi)
       call read_hdf(file, "common/fsamp",  self%samprate)
       call read_hdf(file, "common/polang", polang_buf)
       call read_hdf(file, "common/mbang",  mbang_buf)      


!!$          do j = 1, ndet_tot
!!$             write(*,*) j, trim(dets(j))
!!$          end do

       do i = 1, self%ndet
          do j = 1, ndet_tot
             if(trim(adjustl(detlabels(i))) == trim(adjustl(dets(j)))) then
                exit
             end if
          end do
          if (j > ndet_tot) then
             write(*,*) ' Error -- detector not found in HDF file: ', trim(adjustl(detlabels(i)))
             stop
          end if
          self%polang(i) = polang_buf(j)
          self%mbang(i) = mbang_buf(j)
       end do
       deallocate(polang_buf, mbang_buf, dets)
       call close_hdf_file(file)
    end if
    call mpi_bcast(self%nside,    1,     MPI_INTEGER,          0, self%comm, ierr)
    call mpi_bcast(self%npsi,     1,     MPI_INTEGER,          0, self%comm, ierr)
    call mpi_bcast(self%samprate, 1,     MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(self%polang,   self%ndet,  MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(self%mbang,    self%ndet,  MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

    call wall_time(t1)
    allocate(self%scans(self%nscan))
    do i = 1, self%nscan
       call read_hdf_scan(self%scans(i), self%myid, self%hdfname(i), self%scanid(i), self%ndet, &
            & detlabels)
       do det = 1, self%ndet
          self%scans(i)%d(det)%accept = all(self%scans(i)%d(det)%tod==self%scans(i)%d(det)%tod)
          if (.not. self%scans(i)%d(det)%accept) then
             write(*,fmt='(a,i8,a,i3, i10)') 'Input TOD contain NaN -- scan =', &
                  & self%scanid(i), ', det =', det, count(self%scans(i)%d(det)%tod/=self%scans(i)%d(det)%tod)
             write(*,fmt='(a,a)') '    filename = ', &
                  & trim(self%hdfname(i))
          end if
       end do
    end do

    ! Precompute trigonometric functions
    allocate(self%sin2psi(self%npsi), self%cos2psi(self%npsi))
    allocate(self%psi(self%npsi))
    do i = 1, self%npsi
       psi             = (i-0.5d0)*2.d0*pi/real(self%npsi,dp)
       self%psi(i)     = psi
       self%sin2psi(i) = sin(2.d0*psi)
       self%cos2psi(i) = cos(2.d0*psi)
    end do

    call mpi_barrier(self%comm, ierr)
    call wall_time(t2)
    if (self%myid == 0) write(*,fmt='(a,i4,a,i6,a,f8.1,a)') &
         & '    Myid = ', self%myid, ' -- nscan = ', self%nscan, &
         & ', TOD IO time = ', t2-t1, ' sec'

  end subroutine read_tod

  subroutine read_hdf_scan(self, myid, filename, scan, ndet, detlabels)
    implicit none
    integer(i4b),                   intent(in)    :: myid
    character(len=*),               intent(in)    :: filename
    integer(i4b),                   intent(in)    :: scan, ndet
    class(comm_scan),               intent(inout) :: self
    character(len=*), dimension(:), intent(in)     :: detlabels

    integer(i4b)       :: i, j, k, n, m, nhorn, ext(1), ierr
    real(dp)           :: psi, t1, t2, t3, t4, t_tot(6), scalars(4)
    character(len=6)   :: slabel
    character(len=32)   :: out
    character(len=8)   :: out2
    character(len=128) :: field
    integer(hid_t)     :: nfield, err, obj_type
    type(hdf_file)     :: file
    integer(i4b), allocatable, dimension(:)       :: hsymb
    real(dp),     allocatable, dimension(:)       :: buffer_sp
    integer(i4b), allocatable, dimension(:)       :: htree

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
    call read_hdf(file, slabel // "/common/time",  self%t0)
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
       call read_hdf(file, slabel // "/" // trim(field) // "/scalars",   scalars)
       self%d(i)%gain = scalars(1)
       self%d(i)%sigma0 = scalars(2)
       self%d(i)%fknee = scalars(3)
       self%d(i)%alpha = scalars(4)
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

!!$    if (myid == 0) then
!!$       write(*,*)
!!$       write(*,*) '  IO init   = ', t_tot(1)
!!$       write(*,*) '  IO field  = ', t_tot(2)
!!$       write(*,*) '  IO scalar = ', t_tot(3)
!!$       write(*,*) '  IO tod    = ', t_tot(4)
!!$       write(*,*) '  IO comp   = ', t_tot(5)
!!$       write(*,*) '  IO huff   = ', t_tot(6)
!!$       write(*,*) '  IO total  = ', t4-t3
!!$    end if

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

       n_tot = 0
       open(unit, file=trim(filelist))
       read(unit,*) n
       do i = 1, n
          read(unit,*) j
          if (j >= self%first_scan .and. j <= self%last_scan) n_tot = n_tot+1
       end do
       close(unit)

       if (n_tot == 0) then
          write(*,*) 'Error: No accepted scans in filelist: ', trim(filelist)
          stop
       end if

       open(unit, file=trim(filelist))
       read(unit,*) n
       allocate(id(n_tot), filename(n_tot), scanid(n_tot), weight(n_tot), proc(n_tot), pweight(0:np-1))
       j = 1
       do i = 1, n
          read(unit,*) scanid(j), filename(j), weight(j)
          if (scanid(j) < self%first_scan .or. scanid(j) > self%last_scan) cycle
          id(j) = j
          j     = j+1
          if (j > n_tot) exit
       end do
       close(unit)

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

    npar = 6
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
          output(k,j,5) = merge(1.d0,0.d0,self%scans(i)%d(j)%accept)
          output(k,j,6) = self%scans(i)%d(j)%chisq
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
       call write_hdf(chainfile, trim(adjustl(path))//'accept', output(:,:,5))
       call write_hdf(chainfile, trim(adjustl(path))//'chisq',  output(:,:,6))
       call write_hdf(chainfile, trim(adjustl(path))//'polang', self%polang)
       call write_hdf(chainfile, trim(adjustl(path))//'mono',   self%mono)
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

    npar = 5
    allocate(output(self%nscan_tot,self%ndet,npar))

    if (self%myid == 0) then
       call int2string(iter, itext)
       path = trim(adjustl(itext))//'/tod/'//trim(adjustl(self%freq))//'/'
       call read_hdf(chainfile, trim(adjustl(path))//'gain',     output(:,:,1))
       call read_hdf(chainfile, trim(adjustl(path))//'sigma0',   output(:,:,2))
       call read_hdf(chainfile, trim(adjustl(path))//'alpha',    output(:,:,3))
       call read_hdf(chainfile, trim(adjustl(path))//'fknee',    output(:,:,4))
       call read_hdf(chainfile, trim(adjustl(path))//'accept',   output(:,:,5))
       call read_hdf(chainfile, trim(adjustl(path))//'polang',   self%polang)
       call read_hdf(chainfile, trim(adjustl(path))//'mono',     self%mono)
       call read_hdf(chainfile, trim(adjustl(path))//'bp_delta', self%bp_delta)
    end if

    call mpi_bcast(output, size(output), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%bp_delta, size(self%bp_delta), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%polang, size(self%polang), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%mono, size(self%mono), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)

    do j = 1, self%ndet
       do i = 1, self%nscan
          k             = self%scanid(i)
          self%scans(i)%d(j)%gain   = output(k,j,1)
          self%scans(i)%d(j)%sigma0 = output(k,j,2)
          self%scans(i)%d(j)%alpha  = output(k,j,3)
          self%scans(i)%d(j)%fknee  = output(k,j,4)
          self%scans(i)%d(j)%accept = output(k,j,5) == 1.d0
       end do
    end do

    deallocate(output)
 
  end subroutine initHDF

  function get_det_id(self, label)
    implicit none
    class(comm_tod),                intent(inout) :: self
    character(len=*),               intent(in)    :: label
    integer(i4b)                                  :: get_det_id

    integer(i4b) :: i

    do i = 1, self%ndet
       if (trim(adjustl(label)) == trim(adjustl(self%label(i)))) then
          get_det_id = i
          return
       end if
    end do

    write(*,*) 'Error: Requested detector ', trim(label), ' not found'
    stop

  end function get_det_id

  subroutine initialize_bp_covar(self, filename)
    implicit none
    class(comm_tod),   intent(inout) :: self
    character(len=*),  intent(in)    :: filename

    integer(i4b) :: i, j, k, ndet, npar, n, unit, par
    real(dp)     :: val
    character(len=16)   :: label, det1, det2
    character(len=1024) :: line

    unit = getlun()
    ndet = self%ndet
    npar = size(self%bp_delta,2)

    allocate(self%prop_bp(ndet,ndet,npar))
    allocate(self%prop_bp_mean(npar))
    self%prop_bp      = 0.d0
    self%prop_bp_mean = 0.d0

    open(unit,file=trim(filename))
    do while (.true.)
       read(unit,'(a)',end=34) line
       line = trim(adjustl(line))
       if (line(1:1) == ' ' .or. line(1:1) == '#') then
          cycle
       else if (line(1:4) == 'INIT') then
          read(line,*) label, det1, par, val
          if (trim(adjustl(det1)) == 'MEAN') then
             self%bp_delta(0,par) = val
          else
             j = self%get_det_id(det1)
             self%bp_delta(j,par) = val
          end if
       else if (line(1:4) == 'PROP') then
          read(line,*) label, det1, det2, par, val
          if (trim(adjustl(det1)) == 'MEAN') then
             self%prop_bp_mean(par) = sqrt(val)
          else
             j = self%get_det_id(det1)
             if (trim(adjustl(det2)) == trim(adjustl(det1))) then
                k = j
             else
                k = self%get_det_id(det2)
             end if
             self%prop_bp(j,k,par) = val
             self%prop_bp(k,j,par) = val
          end if
       else
          write(*,*) 'Unsupported entry in ', trim(filename)
          stop
       end if
       
    end do
34  close(unit)
       
    ! Compute square root; mean will be projected out after proposal generation
    do par = 1, npar
       call compute_hermitian_root(self%prop_bp(:,:,par), 0.5d0)
    end do

  end subroutine initialize_bp_covar


end module comm_tod_mod
