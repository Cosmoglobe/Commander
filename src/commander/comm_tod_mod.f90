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
     real(dp)          :: chisq
     integer(i4b)      :: nside                           ! Nside for pixelized pointing
     integer(i4b)      :: npsi                            ! Number of discrete polarization angles
     real(sp),     allocatable, dimension(:)   :: tod     ! Detector values in time domain, (ntod)     
     byte,         allocatable, dimension(:)   :: flag    ! Compressed detector flag; 0 is accepted, /= 0 is rejected
     byte,         allocatable, dimension(:)   :: pix     ! Compressed pixelized pointing, (ntod,nhorn)
     byte,         allocatable, dimension(:)   :: psi     ! Compressed polarization angle, (ntod,nhorn)
     real(dp),     allocatable, dimension(:)   :: N_psd   ! Noise power spectrum density; in uncalibrated units
     real(dp),     allocatable, dimension(:)   :: nu_psd  ! Noise power spectrum bins; in Hz
  end type comm_detscan

  type :: comm_scan
     integer(i4b)   :: ntod                                        ! Number of time samples
     real(dp)       :: v_sun(3)                                    ! Observatory velocity relative to Sun in km/s
     real(dp)       :: t0                                          ! MJD for first sample
     type(huffcode) :: hkey                                        ! Huffman decompression key
     class(comm_detscan), allocatable, dimension(:)     :: d       ! Array of all detectors
  end type comm_scan

  type, abstract :: comm_tod 
     integer(i4b) :: myid, numprocs                               ! MPI parameters
     integer(i4b) :: nmaps                                        ! Number of Stokes parameters
     integer(i4b) :: ndet                                         ! Number of active detectors
     integer(i4b) :: nhorn                                        ! Number of horns
     integer(i4b) :: nscan                                        ! Number of scans
     integer(i4b) :: npsi                                         ! Number of discretized psi steps
     real(dp)     :: samprate                                     ! Sample rate in Hz
     integer(i4b),       allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),           allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     real(sp),           allocatable, dimension(:,:)   :: sin2psi  ! Lookup table of sin(2psi) for each det
     real(sp),           allocatable, dimension(:,:)   :: cos2psi  ! Lookup table of cos(2psi) for each det
     type(comm_scan),    allocatable, dimension(:)     :: scans    ! Array of all scans
     integer(i4b),       allocatable, dimension(:)     :: scanid   ! List of scan IDs
     character(len=512), allocatable, dimension(:)     :: hdfname  ! List of HDF filenames for each ID
     class(comm_map), pointer                          :: procmask ! Mask for gain and n_corr
     class(comm_mapinfo), pointer                      :: info     ! Map definition
     class(comm_mapinfo), pointer                      :: slinfo   ! Sidelobe map info
     class(map_ptr),     allocatable, dimension(:)     :: slbeam   ! Sidelobe beam data (ndet)
     class(conviqt_ptr), allocatable, dimension(:)     :: slconv   ! SL-convolved maps (ndet)
   contains
     procedure                        :: read_tod
     procedure                        :: get_scan_ids
     procedure(process_tod), deferred :: process_tod
  end type comm_tod

  abstract interface
     subroutine process_tod(self, map_in, map_out, rms_out)
       import comm_tod, comm_map, map_ptr
       implicit none
       class(comm_tod),                 intent(inout) :: self
       type(map_ptr),     dimension(:), intent(inout) :: map_in            
       class(comm_map),                 intent(inout) :: map_out, rms_out  
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

  subroutine read_tod(self)
    implicit none
    class(comm_tod),                intent(inout)  :: self

    integer(i4b) :: i, j
    real(dp)     :: t1, t2, psi

    if (self%nscan ==0) return

    call wall_time(t1)
    allocate(self%scans(self%nscan))
    do i = 1, self%nscan
       call read_hdf_scan(self%scans(i), self%hdfname(i), self%scanid(i), self%ndet, &
            & self%npsi)
    end do

    ! Precompute trigonometric functions
    allocate(self%sin2psi(0:self%npsi-1,self%ndet), self%cos2psi(0:self%npsi-1,self%ndet))
    do j = 1, self%ndet
       do i = 0, self%npsi-1
          psi          = (i+0.5d0)*2.d0*pi/real(self%npsi,dp)+self%scans(1)%d(j)%polang
          self%sin2psi(i,j) = sin(2.d0*psi)
          self%cos2psi(i,j) = cos(2.d0*psi)
       end do
    end do


    call wall_time(t2)
    if (self%myid == 0) write(*,fmt='(a,i4,a,i6,a,f8.1,a)') &
         & '    Myid = ', self%myid, ' -- nscan = ', self%nscan, &
         & ', TOD IO time = ', t2-t1, ' sec'

  end subroutine read_tod

  subroutine read_hdf_scan(self, filename, scan, ndet, npsi)
    implicit none
    character(len=*),                              intent(in)    :: filename
    integer(i4b),                                  intent(in)    :: scan, ndet
    integer(i4b),                                  intent(out)   :: npsi
    class(comm_scan),                              intent(inout) :: self

    integer(i4b)       :: i, j, k, n, m, nhorn, ext(1), ierr
    real(dp)           :: psi
    character(len=6)   :: slabel
    character(len=32)   :: out
    character(len=8)   :: out2
    character(len=128) :: field
    integer(hid_t)     :: nfield, err, obj_type
    type(hdf_file)     :: file
    integer(i4b), allocatable, dimension(:) :: hsymb, hfreq
    real(dp),     allocatable, dimension(:) :: buffer_sp

    call int2string(scan, slabel)

    call open_hdf_file(filename, file, "r")

    ! Find array sizes
    call h5gget_obj_info_idx_f(file%filehandle, slabel, 1, field, obj_type, err)
    if (trim(field) == 'common') then
       call h5gget_obj_info_idx_f(file%filehandle, slabel, 2, field, obj_type, err)
    end if
    call get_size_hdf(file, slabel // "/" // trim(field) // "/tod", ext)
    !nhorn     = ext(1)
    n         = ext(1)
    m         = get_closest_fft_magic_number(n)
    self%ntod = m

    ! Read common scan data
    call read_hdf(file, slabel // "/common/vsun", self%v_sun)
    call read_hdf(file, slabel // "/common/time", self%t0)

    ! Read detector scans
    allocate(self%d(ndet), buffer_sp(n))
    i = 0
    do j = 0, ndet
       call h5gget_obj_info_idx_f(file%filehandle, slabel, j, field, obj_type, err)
       if (trim(field) == 'common') cycle
       i = i+1
!write(*,*) i, allocated(self%d(i)%tod)
       allocate(self%d(i)%tod(m))
       self%d(i)%label = trim(field)
       call read_hdf(file, slabel // "/" // trim(field) // "/gain",   self%d(i)%gain)
       call read_hdf(file, slabel // "/" // trim(field) // "/sigma0", self%d(i)%sigma0)
       call read_hdf(file, slabel // "/" // trim(field) // "/alpha",  self%d(i)%alpha)
       call read_hdf(file, slabel // "/" // trim(field) // "/fknee",  self%d(i)%fknee)
       call read_hdf(file, slabel // "/" // trim(field) // "/tod",    buffer_sp)
       self%d(i)%tod = buffer_sp(1:m)
       call read_hdf(file, slabel // "/" // trim(field) // "/nside",  self%d(i)%nside)
       call read_hdf(file, slabel // "/" // trim(field) // "/polang", self%d(i)%polang)
       call read_hdf(file, slabel // "/common/fsamp",                 self%d(i)%samprate)
       call read_hdf(file, slabel // "/" // trim(field) // "/npsi",  npsi)
   
       ! Read Huffman coded data arrays
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/pix",  self%d(i)%pix)
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/psi",  self%d(i)%psi)
       call read_hdf_opaque(file, slabel // "/" // trim(field) // "/flag", self%d(i)%flag)
!!$       write(*,'(10z2)') self%d(i)%pix(1:10)
!!$
!!$    do k = 1, 8
!!$       if (btest(self%d(i)%pix(1),k-1)) then
!!$          write(*,*) k,'1'
!!$       else
!!$          write(*,*) k,'0'
!!$       end if
!!$    end do
!!$
!!$       write(*,*) 
!!$       stop

!!$       if (scan==1) then
!!$          do k=0, 31
!!$             if (btest(self%d(i)%flag(1),k)) then
!!$                out(k+1:k+1) = '1'
!!$             else
!!$                out(k+1:k+1) = '0'
!!$             end if
!!$          end do
!!$          write(*,*)  self%d(i)%flag(1), out
!!$
!!$          do k=0, 7
!!$             if (btest(self%d(i)%flag2(2),k)) then
!!$                out2(k+1:k+1) = '1'
!!$             else
!!$                out2(k+1:k+1) = '0'
!!$             end if
!!$          end do
!!$          write(*,*)  self%d(i)%flag2(2), out2
!!$
!!$       end if
    end do
    deallocate(buffer_sp)

    if (i /= ndet) then
       write(*,*) 'ERROR: Too few detectors in TOD file = ', trim(filename)
       stop
    end if

    ! Initialize Huffman key
    call read_alloc_hdf(file, slabel // "/common/huffsymb", hsymb)
    call read_alloc_hdf(file, slabel // "/common/hufffreq", hfreq)
!!$    if (scan==1) then
!!$    do i =1, size(hsymb)
!!$       write(*,*) i, hsymb(i), hfreq(i)
!!$    end do
!!$ end if
    call hufmak(hsymb,hfreq,self%hkey)
    deallocate(hsymb, hfreq)

    call close_hdf_file(file)

  end subroutine read_hdf_scan

  subroutine get_scan_ids(self, filelist)
    implicit none
    class(comm_tod),   intent(inout) :: self    
    character(len=*),  intent(in)    :: filelist

    integer(i4b)       :: unit, j, k
    integer(hid_t)     :: i, n, n_tot, err, obj_type
    character(len=128) :: label
    character(len=512) :: filename
    type(hdf_file)     :: file

    ! Find total number of scans for this core
    unit = getlun()
    n_tot = 0
    j     = 0
    open(unit, file=trim(filelist))
    do while (.true.)
       read(unit,'(a)',end=58) filename
       j = j+1
       if (mod(j-1,self%numprocs) /= self%myid) cycle
       call open_hdf_file(filename, file, "r")
       call h5gn_members_f(file%filehandle, "/", n, err)
       call close_hdf_file(file)
       n_tot = n_tot + n
    end do
58  close(unit)

    self%nscan = n_tot
    if (n_tot == 0) return

    ! Read in all filenames and associated scan ids
    allocate(self%scanid(n_tot), self%hdfname(n_tot))
    j     = 0
    k     = 0
    open(unit, file=trim(filelist))
    do while (.true.)
       read(unit,'(a)',end=59) filename
       j = j+1
       if (mod(j-1,self%numprocs) /= self%myid) cycle
       call open_hdf_file(filename, file, "r")
       call h5gn_members_f(file%filehandle, "/", n, err)
       do i = 0, n-1
          call h5gget_obj_info_idx_f(file%filehandle, "/", i, label, obj_type, err)
          k            = k+1
          self%hdfname(k) = filename
          read(label,*) self%scanid(k)
       end do
       call close_hdf_file(file)
    end do
59  close(unit)

  end subroutine get_scan_ids

end module comm_tod_mod
