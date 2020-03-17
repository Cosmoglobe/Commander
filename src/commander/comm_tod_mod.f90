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
  public comm_tod, initialize_tod_mod, fill_all_masked

  type :: comm_detscan
     character(len=10) :: label                           ! Detector label
     real(dp)          :: gain, dgain, gain_sigma         ! Gain; assumed constant over scan
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
     integer(i4b)      :: ext_lowres(2)             ! Shape of downgraded TOD including padding
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
     character(len=512) :: initfile
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

     real(dp)     :: central_freq                                 !Central frequency
     real(dp)     :: samprate, samprate_lowres                      ! Sample rate in Hz
     logical(lgt) :: orb_abscal
     real(dp), allocatable, dimension(:)     :: gain0                                      ! Mean gain
     real(dp), allocatable, dimension(:)     :: polang                                      ! Detector polarization angle
     real(dp), allocatable, dimension(:)     :: mbang                                       ! Main beams angle
     real(dp), allocatable, dimension(:)     :: mono                                        ! Monopole
     real(dp), allocatable, dimension(:)     :: fwhm, elip, psi_ell, mb_eff                         ! Beam parameter
     real(dp), allocatable, dimension(:)     :: nu_c                                        ! Center frequency
     real(dp), allocatable, dimension(:,:,:) :: prop_bp         ! proposal matrix, L(ndet,ndet,ndelta),  for bandpass sampler
     real(dp), allocatable, dimension(:)     :: prop_bp_mean    ! proposal matrix, sigma(ndelta), for mean
     integer(i4b)      :: nside                           ! Nside for pixelized pointing
     integer(i4b)      :: nobs                            ! Number of observed pixeld for this core
     integer(i4b) :: output_n_maps                                ! Output n_maps
     character(len=512) :: init_from_HDF                          ! Read from HDF file
     integer(i4b) :: output_4D_map                                ! Output 4D maps
     logical(lgt) :: subtract_zodi                                ! Subtract zodical light
     integer(i4b),       allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),           allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     real(sp),           allocatable, dimension(:)     :: sin2psi  ! Lookup table of sin(2psi) 
     real(sp),           allocatable, dimension(:)     :: cos2psi  ! Lookup table of cos(2psi) 
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
     class(comm_map), pointer                          :: procmask => null() ! Mask for gain and n_corr
     class(comm_map), pointer                          :: procmask2 => null() ! Mask for gain and n_corr
     class(comm_mapinfo), pointer                      :: info => null()    ! Map definition
     class(comm_mapinfo), pointer                      :: slinfo => null()  ! Sidelobe map info
     class(map_ptr),     allocatable, dimension(:)     :: slbeam   ! Sidelobe beam data (ndet)
     class(conviqt_ptr), allocatable, dimension(:)     :: slconv   ! SL-convolved maps (ndet)
     real(dp),           allocatable, dimension(:,:)   :: bp_delta  ! Bandpass parameters (0:ndet, npar)
     real(dp),           allocatable, dimension(:,:)   :: spinaxis ! For load balancing
     integer(i4b),       allocatable, dimension(:)     :: pix2ind, ind2pix, ind2sl
     real(sp),           allocatable, dimension(:,:)   :: ind2ang
     real(dp),           allocatable, dimension(:, :) :: orb_dp_s !precomputed s integrals for orbital dipole sidelobe term 
   contains
     procedure                        :: read_tod
     procedure                        :: get_scan_ids
     procedure                        :: dumpToHDF
     procedure                        :: initHDF
     procedure                        :: get_det_id
     procedure                        :: initialize_bp_covar
     procedure(process_tod), deferred :: process_tod
     procedure                        :: construct_sl_template
     procedure                        :: output_scan_list
     procedure                        :: downsample_tod
     procedure                        :: compute_chisq
     procedure                        :: get_total_chisq
     procedure                        :: symmetrize_flags
     procedure                        :: decompress_pointing_and_flags
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

    integer(i4b) :: i, j, n, det, ierr, ndet_tot
    real(dp)     :: t1, t2
    real(sp)     :: psi
    type(hdf_file)     :: file

    integer(i4b), dimension(:), allocatable       :: ns   
    real(dp), dimension(:), allocatable           :: mbang_buf, polang_buf
    character(len=100000)                         :: det_buf
    character(len=128), dimension(:), allocatable :: dets


    ! Read common fields
    allocate(self%polang(self%ndet), self%mbang(self%ndet), self%mono(self%ndet), self%gain0(0:self%ndet))
    self%mono = 0.d0
    if (self%myid == 0) then
       call open_hdf_file(self%initfile, file, "r")


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
       call read_hdf_scan(self%scans(i), self, self%hdfname(i), self%scanid(i), self%ndet, &
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

    ! Initialize mean gain
    allocate(ns(0:self%ndet))
    self%gain0 = 0.d0
    ns         = 0
    do i = 1, self%nscan
       do j = 1, self%ndet 
          if (.not. self%scans(i)%d(j)%accept) cycle
          self%gain0(j) = self%gain0(j) + self%scans(i)%d(j)%gain
          ns(j)         = ns(j) + 1
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, self%gain0, self%ndet+1, &
         & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, ns,         self%ndet+1, &
         & MPI_INTEGER,          MPI_SUM, self%comm, ierr)
    self%gain0(0) = sum(self%gain0)/sum(ns)
    where (ns > 0)
       self%gain0 = self%gain0 / ns - self%gain0(0)
    end where

    do i = 1, self%nscan
       do j = 1, self%ndet 
          self%scans(i)%d(j)%dgain = self%scans(i)%d(j)%gain - self%gain0(0) - self%gain0(j)
!          self%scans(i)%d(j)%dgain = 0.d0
!          self%scans(i)%d(j)%gain  = self%gain0(0) + self%gain0(j)
       end do
    end do

    ! Precompute trigonometric functions
    allocate(self%sin2psi(self%npsi), self%cos2psi(self%npsi))
    allocate(self%psi(self%npsi))
    do i = 1, self%npsi
       psi             = (i-0.5)*2.0*pi/real(self%npsi,sp)
       self%psi(i)     = psi
       self%sin2psi(i) = sin(2.0*psi)
       self%cos2psi(i) = cos(2.0*psi)
    end do

    call mpi_barrier(self%comm, ierr)
    call wall_time(t2)
    if (self%myid == 0) write(*,fmt='(a,i4,a,i6,a,f8.1,a)') &
         & '    Myid = ', self%myid, ' -- nscan = ', self%nscan, &
         & ', TOD IO time = ', t2-t1, ' sec'

  end subroutine read_tod

  subroutine read_hdf_scan(self, tod, filename, scan, ndet, detlabels)
    implicit none
    class(comm_scan),               intent(inout) :: self
    class(comm_tod),                intent(in)    :: tod
    character(len=*),               intent(in)    :: filename
    integer(i4b),                   intent(in)    :: scan, ndet
    character(len=*), dimension(:), intent(in)     :: detlabels

    integer(i4b)       :: i, n, m, ext(1)
    real(dp)           :: t1, t2, t3, t4, t_tot(6), scalars(4)
    character(len=6)   :: slabel
    character(len=128) :: field
    type(hdf_file)     :: file
    integer(i4b), allocatable, dimension(:)       :: hsymb
    real(sp),     allocatable, dimension(:)       :: buffer_sp
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
    !m = n
    m         = get_closest_fft_magic_number(n)
!!$    m         = get_closest_fft_magic_number(2*n)
!!$    do while (mod(m,2) == 1)
!!$       m = get_closest_fft_magic_number(m-1)
!!$    end do
!!$    m = m/2
    self%ntod = m
    self%ext_lowres(1)   = -5    ! Lowres padding
    self%ext_lowres(2)   = int(self%ntod/int(tod%samprate/tod%samprate_lowres)) + 1 + self%ext_lowres(1)

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

    integer(i4b)       :: unit, j, k, np, ind(1), i, n, m, n_tot, ierr
    real(dp)           :: w_tot, w, v0(3), v(3)
    real(dp),           allocatable, dimension(:)   :: weight, sid
    real(dp),           allocatable, dimension(:,:) :: spinpos, spinaxis
    integer(i4b),       allocatable, dimension(:)   :: scanid, id
    integer(i4b),       allocatable, dimension(:)   :: proc
    real(dp),           allocatable, dimension(:)   :: pweight
    character(len=512), allocatable, dimension(:)   :: filename

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
       else if (n_tot==1) then
         self%nscan = n_tot
         open(unit, file=trim(filelist))
         read(unit,*) n
         allocate(filename(n_tot), scanid(n_tot), proc(n_tot), spinpos(2,n_tot), weight(n_tot))
         j = 1
         do i = 1, n
            read(unit,*) scanid(j), filename(j), weight(j), spinpos(1:2,j)
            if (scanid(j) < self%first_scan .or. scanid(j) > self%last_scan) exit
         end do
         proc(1) = 0
         self%initfile = filename(1)
         close(unit)

       else

       
         open(unit, file=trim(filelist))
         read(unit,*) n
         allocate(id(n_tot), filename(n_tot), scanid(n_tot), weight(n_tot), proc(n_tot), pweight(0:np-1), sid(n_tot), spinaxis(n_tot,3), spinpos(2,n_tot))
         j = 1
         do i = 1, n
            read(unit,*) scanid(j), filename(j), weight(j), spinpos(1:2,j)
            if (scanid(j) < self%first_scan .or. scanid(j) > self%last_scan) cycle
            id(j)  = j
            sid(j) = scanid(j)
            call ang2vec(spinpos(1,j), spinpos(2,j), spinaxis(j,1:3))
            if (j == 1) self%initfile = filename(j)
            j      = j+1
            if (j > n_tot) exit
         end do
         close(unit)

         ! Compute symmetry axis
         v0 = 0.d0
         do i = 2, n_tot
            v(1) = spinaxis(1,2)*spinaxis(i,3)-spinaxis(1,3)*spinaxis(i,2)
            v(2) = spinaxis(1,3)*spinaxis(i,1)-spinaxis(1,1)*spinaxis(i,3)
            v(3) = spinaxis(1,1)*spinaxis(i,2)-spinaxis(1,2)*spinaxis(i,1)
            if (v(3) < 0.d0) v  = -v
            if (sum(v*v) > 0.d0)  v0 = v0 + v / sqrt(sum(v*v))
         end do
         v0 = v0 / sqrt(v0*v0)
!        v0(1) = 1
       

!!$
!!$       ! Compute angle between i'th and first vector
!!$       do i = 1, n
!!$          v(1) = spinaxis(1,2)*spinaxis(i,3)-spinaxis(1,3)*spinaxis(i,2)
!!$          v(2) = spinaxis(1,3)*spinaxis(i,1)-spinaxis(1,1)*spinaxis(i,3)
!!$          v(3) = spinaxis(1,1)*spinaxis(i,2)-spinaxis(1,2)*spinaxis(i,1)          
!!$       end do
         do i = n_tot, 1, -1
            v(1) = spinaxis(1,2)*spinaxis(i,3)-spinaxis(1,3)*spinaxis(i,2)
            v(2) = spinaxis(1,3)*spinaxis(i,1)-spinaxis(1,1)*spinaxis(i,3)
            v(3) = spinaxis(1,1)*spinaxis(i,2)-spinaxis(1,2)*spinaxis(i,1)
            sid(i) = acos(max(min(sum(spinaxis(i,:)*spinaxis(1,:)),1.d0),-1.d0))
            if (sum(v*v0) < 0.d0) sid(i) = -sid(i) ! Flip sign 
         end do

!!$       ! Sort according to weight
!!$       pweight = 0.d0
!!$       call QuickSort(id, weight)
!!$       do i = n_tot, 1, -1
!!$          ind             = minloc(pweight)-1
!!$          proc(id(i))     = ind(1)
!!$          pweight(ind(1)) = pweight(ind(1)) + weight(i)
!!$       end do
!!$       deallocate(id, pweight, weight)

       ! Sort according to scan id
         proc    = -1
         call QuickSort(id, sid)
         w_tot = sum(weight)
         j     = 1
         do i = np-1, 1, -1
            w = 0.d0
            do k = 1, n_tot
               if (proc(k) == i) w = w + weight(k) 
            end do
            do while (w < w_tot/np .and. j <= n_tot)
               proc(id(j)) = i
               w           = w + weight(id(j))
               if (w > 1.2d0*w_tot/np) then
                  ! Assign large scans to next core
                  proc(id(j)) = i-1
                  w           = w - weight(id(j))
               end if
               j           = j+1
            end do
         end do
         do while (j <= n_tot)
            proc(id(j)) = 0
            j = j+1
         end do
         pweight = 0.d0
         do k = 1, n_tot
            pweight(proc(id(k))) = pweight(proc(id(k))) + weight(id(k))
         end do
         write(*,*) '  Min/Max core weight = ', minval(pweight)/w_tot*np, maxval(pweight)/w_tot*np
         deallocate(id, pweight, weight, sid, spinaxis)

       ! Distribute according to consecutive PID
!!$       do i = 1, n_tot
!!$          proc(i) = max(min(int(real(i-1,sp)/real(n_tot-1,sp) * np),np-1),0)
!!$       end do

      end if
   end if

    call mpi_bcast(n_tot, 1,  MPI_INTEGER, 0, self%comm, ierr)
    if (self%myid /= 0) then
       allocate(filename(n_tot), scanid(n_tot), proc(n_tot), spinpos(2,n_tot))
    end if
    call mpi_bcast(filename, 512*n_tot,  MPI_CHARACTER, 0, self%comm, ierr)
    call mpi_bcast(scanid,       n_tot,  MPI_INTEGER,   0, self%comm, ierr)
    call mpi_bcast(proc,         n_tot,  MPI_INTEGER,   0, self%comm, ierr)
    call mpi_bcast(spinpos,    2*n_tot,  MPI_DOUBLE_PRECISION,   0, self%comm, ierr)

    self%nscan     = count(proc == self%myid)
    allocate(self%scanid(self%nscan), self%hdfname(self%nscan), self%spinaxis(self%nscan,2))
    j = 1
    do i = 1, n_tot
       if (proc(i) == self%myid) then
          self%scanid(j)     = scanid(i)
          self%hdfname(j)    = filename(i)
          self%spinaxis(j,:) = spinpos(:,i)
          j                  = j+1
       end if
    end do

    if (self%myid == 0) then
       allocate(self%nscanprproc(0:self%numprocs-1))
       do i = 0, self%numprocs-1
          self%nscanprproc(i) = count(proc == i)
       end do
    end if

    deallocate(filename, scanid, proc, spinpos) 

  end subroutine get_scan_ids


  subroutine dumpToHDF(self, chainfile, iter, map, rms)
    implicit none
    class(comm_tod),                   intent(in)    :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile
    class(comm_map),                   intent(in)    :: map, rms

    integer(i4b)       :: i, j, k, l, npar, ierr
    real(dp)           :: mu
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
       ! Fill in defaults (closest previous)
       do j = 1, self%ndet
          do i = 1, 4
             do k = 1, self%nscan_tot
                if (output(k,j,i) == 0.d0) then
                   l = k
                   if (k == 1) then
                      do while (output(l,j,i) == 0.d0 .and. l < self%nscan) 
                         l = l+1
                      end do
                   else
                      do while (output(l,j,i) == 0.d0 .and. l > 1) 
                         l = l-1
                      end do
                   end if
                   output(k,j,i) = output(l,j,i)
                end if
             end do
!!$             if (output(
!!$             mu = sum(output(:,j,i)) / count(output(:,j,i) /= 0.d0)
!!$             where (output(:,j,i) == 0.d0) 
!!$                output(:,j,i) = mu
!!$             end where
          end do
       end do

!!$       do j = 1, self%ndet
!!$          do i = 1, 4
!!$             mu = sum(output(:,j,i)) / count(output(:,j,i) /= 0.d0)
!!$             where (output(:,j,i) == 0.d0) 
!!$                output(:,j,i) = mu
!!$             end where
!!$          end do
!!$       end do

       call int2string(iter, itext)
       path = trim(adjustl(itext))//'/tod/'//trim(adjustl(self%freq))//'/'
       !write(*,*) 'path', trim(path)
       call create_hdf_group(chainfile, trim(adjustl(path)))
       call write_hdf(chainfile, trim(adjustl(path))//'gain',   output(:,:,1))
       call write_hdf(chainfile, trim(adjustl(path))//'sigma0', output(:,:,2))
       call write_hdf(chainfile, trim(adjustl(path))//'alpha',  output(:,:,3))
       call write_hdf(chainfile, trim(adjustl(path))//'fknee',  output(:,:,4))
       call write_hdf(chainfile, trim(adjustl(path))//'accept', output(:,:,5))
       call write_hdf(chainfile, trim(adjustl(path))//'chisq',  output(:,:,6))
       call write_hdf(chainfile, trim(adjustl(path))//'polang', self%polang)
       call write_hdf(chainfile, trim(adjustl(path))//'gain0',  self%gain0)
       call write_hdf(chainfile, trim(adjustl(path))//'mono',   self%mono)
       call write_hdf(chainfile, trim(adjustl(path))//'bp_delta', self%bp_delta)
    end if

    call map%writeMapToHDF(chainfile, path, 'map')
    call rms%writeMapToHDF(chainfile, path, 'rms')

    deallocate(output)
 
  end subroutine dumpToHDF

  subroutine initHDF(self, chainfile, iter, map, rms)
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile
    class(comm_map),                   intent(inout) :: map, rms

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
       call read_hdf(chainfile, trim(adjustl(path))//'gain0',    self%gain0)

       ! Redefine gains; should be removed when proper initfiles are available
       self%gain0(0) = sum(output(:,:,1))/count(output(:,:,1)>0.d0)
       !write(*,*) self%gain0(0), minval(output(:,:,1)), maxval(output(:,:,1))
       !stop
       do i = 1, self%ndet
          self%gain0(i) = sum(output(:,i,1))/count(output(:,i,1)>0.d0) - self%gain0(0)
       end do
    end if

    call mpi_bcast(output, size(output), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%bp_delta, size(self%bp_delta), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%polang, size(self%polang), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%mono, size(self%mono), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)
    call mpi_bcast(self%gain0, size(self%gain0), MPI_DOUBLE_PRECISION, 0, &
         & self%comm, ierr)

    do j = 1, self%ndet
       do i = 1, self%nscan
          k             = self%scanid(i)
          self%scans(i)%d(j)%gain   = output(k,j,1)
          self%scans(i)%d(j)%dgain  = output(k,j,1)-self%gain0(0)-self%gain0(j)
          self%scans(i)%d(j)%sigma0 = output(k,j,2)
          self%scans(i)%d(j)%alpha  = output(k,j,3)
          self%scans(i)%d(j)%fknee  = output(k,j,4)
          self%scans(i)%d(j)%accept = output(k,j,5) == 1.d0
       end do
    end do

    call map%readMapFromHDF(chainfile, trim(adjustl(path))//'map')
    call rms%readMapFromHDF(chainfile, trim(adjustl(path))//'rms')

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

    get_det_id = -1

  end function get_det_id

  subroutine initialize_bp_covar(self, filename)
    implicit none
    class(comm_tod),   intent(inout) :: self
    character(len=*),  intent(in)    :: filename

    integer(i4b) :: j, k, ndet, npar, unit, par
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
             if (j > 0) self%bp_delta(j,par) = val
          end if
       else if (line(1:4) == 'PROP') then
          read(line,*) label, det1, det2, par, val
          if (trim(adjustl(det1)) == 'MEAN') then
             self%prop_bp_mean(par) = sqrt(val)
          else
             j = self%get_det_id(det1)
             if (j < 0) cycle
             if (trim(adjustl(det2)) == trim(adjustl(det1))) then
                k = j
             else
                k = self%get_det_id(det2)
             end if
             if (k < 0) cycle
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


  !construct a sidelobe template in the time domain
  subroutine construct_sl_template(self, slconv, pix, psi, s_sl, polangle)
    implicit none
    class(comm_tod),                     intent(in)    :: self
    class(comm_conviqt),                 intent(in)    :: slconv
    integer(i4b),        dimension(:),   intent(in)    :: pix, psi
    real(dp),                            intent(in)   :: polangle
    real(sp),            dimension(:),   intent(out)   :: s_sl
    
    integer(i4b) :: j, pix_, pix_prev, psi_prev
    real(dp)     :: psi_

    pix_prev = -1; psi_prev = -1
    do j=1, size(pix)
       pix_    = self%ind2sl(self%pix2ind(pix(j)))
       if (pix_prev == pix_ .and. psi(j) == psi_prev) then
          s_sl(j) = s_sl(j-1)
       else
          psi_    = self%psi(psi(j))-polangle 
          s_sl(j) = slconv%interp(pix_, psi_)
          pix_prev = pix_; psi_prev = psi(j)
       end if
    end do

  end subroutine construct_sl_template

  subroutine output_scan_list(self, slist)
    implicit none
    class(comm_tod),                               intent(in)    :: self
    character(len=512), allocatable, dimension(:), intent(inout) :: slist

    integer(i4b)     :: i, j, mpistat(MPI_STATUS_SIZE), unit, ns, ierr
    character(len=4) :: pid

    if (self%myid == 0) then
       call int2string(self%myid, pid)
       unit = getlun()
       open(unit,file=trim(self%outdir)//'/filelist_'//trim(self%freq)//'.txt', recl=512)
       write(unit,*) sum(self%nscanprproc)
       do i = 1, self%nscan
          if (trim(slist(i)) == '') cycle
          write(unit,*) trim(slist(i))
       end do
       deallocate(slist)
       do j = 1, self%numprocs-1
          ns = self%nscanprproc(j)
          allocate(slist(ns))
          call mpi_recv(slist, 512*ns, MPI_CHARACTER, j, 98, self%comm, mpistat, ierr)
          do i = 1, ns
             if (trim(slist(i)) == '') cycle
             write(unit,*) trim(slist(i))
          end do
          deallocate(slist)
       end do
       close(unit)
    else
       call mpi_send(slist, 512*self%nscan, MPI_CHARACTER, 0, 98, self%comm, ierr)
       deallocate(slist)
    end if
  end subroutine output_scan_list



  subroutine downsample_tod(self, tod_in, ext, tod_out, mask)
    implicit none
    class(comm_tod),                    intent(in)     :: self
    real(sp), dimension(:),             intent(in)     :: tod_in
    integer(i4b),                       intent(inout)  :: ext(2)
    real(sp), dimension(ext(1):ext(2)), intent(out), optional :: tod_out
    real(sp), dimension(:),             intent(in),  optional :: mask
 
    integer(i4b) :: i, j, k, n, step, ntod, w, npad

    ntod = size(tod_in)
    npad = 5
    step = int(self%samprate/self%samprate_lowres)
    w    = 2*step    ! Boxcar window width
    n    = int(ntod / step) + 1
    if (.not. present(tod_out)) then
       ext = [-npad, n+npad]
       return
    end if

    do i = -npad, n+npad
       j = max(i*step - w, 1)
       k = min(i*step + w, ntod)

       if (j > k) then
          tod_out(i) = 0.
       else
          !write(*,*) i, shape(tod_in), j, k 
          if (present(mask)) then
             tod_out(i) = sum(tod_in(j:k)*mask(j:k)) / (2*w+1)
          else
             tod_out(i) = sum(tod_in(j:k)) / (2*w+1)
          end if
       end if
       !write(*,*) i, tod_out(i), sum(mask(j:k)), sum(tod_in(j:k))
    end do

!!$    if (self%myid == 0) then
!!$       open(58, file='filter.dat')
!!$!       do i = 1, ntod
!!$!          write(58,*) i, tod_in(i)
!!$!       end do
!!$!       write(58,*)
!!$       do i = -npad, n+npad
!!$          write(58,*) i*step, tod_out(i)
!!$       end do
!!$       close(58)
!!$    end if
!!$    call mpi_finalize(i)
!!$    stop

  end subroutine downsample_tod



! Fills masked region with linear function between the mean of 20 points at each end
  subroutine fill_masked_region(d_p, mask, i_start, i_end, ntod)
    implicit none
    real(sp), intent(inout)  :: d_p(:)
    real(sp), intent(in)     :: mask(:)
    integer(i4b), intent(in) :: i_end, i_start, ntod 
    real(dp)     :: mu1, mu2
    integer(i4b) :: i, n_mean, earliest, latest
    n_mean = 20
    earliest = max(i_start - (n_mean + 1), 1)
    latest = min(i_end + (n_mean + 1), ntod)
    if (i_start == 1) then  ! masked region at start for scan
       if (i_end < ntod) then
          mu2 = sum(d_p(i_end:latest) * mask(i_end:latest)) / sum(mask(i_end:latest))
          d_p(i_start:i_end) = mu2
       else
          write(*,*) "Entire scan masked, this should not happen (in comm_tod_mod.fill_masked_region)"
       end if
    else if (i_end == ntod) then  ! masked region at end of scan
       mu1 = sum(d_p(earliest:i_start) * mask(earliest:i_start)) / sum(mask(earliest:i_start))
       d_p(i_start:i_end) = mu1
    else   ! masked region in middle of scan
       mu1 = sum(d_p(earliest:i_start) * mask(earliest:i_start)) / sum(mask(earliest:i_start))
       mu2 = sum(d_p(i_end:latest) * mask(i_end:latest)) / sum(mask(i_end:latest))
       do i = i_start, i_end
          d_p(i) = mu1 + (mu2 - mu1) * (i - i_start + 1) / (i_end - i_start + 2)
       end do
    end if
  end subroutine fill_masked_region


! Identifies and fills masked region
  subroutine fill_all_masked(d_p, mask, ntod, sample, sigma_0, handle)
    implicit none
    real(sp),         intent(inout)  :: d_p(:)
    real(sp),         intent(in)     :: mask(:)
    real(sp),         intent(in), optional     :: sigma_0
    type(planck_rng), intent(inout), optional  :: handle
    integer(i4b),     intent(in) :: ntod
    logical(lgt),     intent(in) :: sample
    integer(i4b) :: j_end, j_start, j, k
    logical(lgt) :: init_masked_region, end_masked_region
    
    ! Fill gaps in data 
    init_masked_region = .true.
    end_masked_region = .false.
    do j = 1,ntod
       if (mask(j) == 1.) then
          if (end_masked_region) then
             j_end = j - 1
             call fill_masked_region(d_p, mask, j_start, j_end, ntod)
             ! Add noise to masked region
             if (sample) then
                do k = j_start, j_end
                   d_p(k) = d_p(k) + sigma_0 * rand_gauss(handle)
                end do
             end if
             end_masked_region = .false.
             init_masked_region = .true.
          end if
       else
          if (init_masked_region) then
             init_masked_region = .false.
             end_masked_region = .true.
             j_start = j
          end if
       end if
    end do
    ! if the data ends with a masked region
    if (end_masked_region) then
       j_end = ntod
       call fill_masked_region(d_p, mask, j_start, j_end, ntod)
       if (sample) then
          do k = j_start, j_end
             d_p(k) = d_p(k) + sigma_0 * rand_gauss(handle)
          end do
       end if
    end if

  end subroutine fill_all_masked


  ! Compute chisquare
  subroutine compute_chisq(self, scan, det, mask, s_sky, s_spur, &
       & n_corr, absbp)
    implicit none
    class(comm_tod),                 intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sky, s_spur
    real(sp),          dimension(:), intent(in)     :: n_corr
    logical(lgt),                    intent(in), optional :: absbp

    real(dp)     :: chisq, d0, g
    integer(i4b) :: i, n

    chisq       = 0.d0
    n           = 0
    g           = self%scans(scan)%d(det)%gain 
    do i = 1, self%scans(scan)%ntod
       if (mask(i) < 0.5) cycle 
       n     = n+1
       d0    = self%scans(scan)%d(det)%tod(i) - &
            & (g * s_spur(i) + n_corr(i))
       chisq = chisq + (d0 - g * s_sky(i))**2 
    end do

    if (self%scans(scan)%d(det)%sigma0 <= 0.d0) then
       if (present(absbp)) then
          self%scans(scan)%d(det)%chisq_prop   = 0.d0
       else
          self%scans(scan)%d(det)%chisq        = 0.d0
       end if
    else
       chisq      = chisq      / self%scans(scan)%d(det)%sigma0**2
       if (present(absbp)) then
          self%scans(scan)%d(det)%chisq_prop   = chisq
       else
          self%scans(scan)%d(det)%chisq        = (chisq - n) / sqrt(2.d0*n)
       end if
    end if
    ! write(*,*) "chi2 :  ", scan, det, self%scanid(scan), &
    !      & self%scans(scan)%d(det)%chisq, self%scans(scan)%d(det)%sigma0
    !if(self%scans(scan)%d(det)%chisq > 2000.d0 .or. isNaN(self%scans(scan)%d(det)%chisq)) then
      !write(*,*) "chisq", scan, det, sum(mask), sum(s_sky), sum(s_sl), sum(s_orb), sum(n_corr)
    !end if
    
  end subroutine compute_chisq

  function get_total_chisq(self, det)
    implicit none
    class(comm_tod),     intent(in)  :: self
    integer(i4b),        intent(in)  :: det
    real(dp)                         :: get_total_chisq

    integer(i4b) :: i, ierr
    real(dp)     :: chisq, buffer

    chisq = 0.d0
    do i = 1, self%nscan ! Sum chisq for all scans
       if (.not. self%scans(i)%d(det)%accept) cycle
       chisq = chisq + self%scans(i)%d(det)%chisq 
    end do
    call mpi_reduce(chisq, buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         & 0, self%info%comm, ierr)

    get_total_chisq = buffer

  end function get_total_chisq

  subroutine symmetrize_flags(self, flag)
    implicit none
    class(comm_tod),                          intent(inout) :: self
    integer(i4b),         dimension(1:,1:),   intent(inout) :: flag

    integer(i4b) :: i, det

    do det = 1, self%ndet
       do i = 1, size(flag,2)
          if (iand(flag(i,det),self%flag0) .ne. 0) then
             flag(i,self%partner(det)) = flag(i,det)
          end if
       end do
    end do

  end subroutine symmetrize_flags

  subroutine decompress_pointing_and_flags(self, scan, det, pix, psi, flag)
    implicit none
    class(comm_tod),                    intent(in)  :: self
    integer(i4b),                       intent(in)  :: scan, det
    integer(i4b),        dimension(:),  intent(out) :: pix, psi, flag

    call huffman_decode2(self%scans(scan)%hkey, self%scans(scan)%d(det)%pix,  pix)
    call huffman_decode2(self%scans(scan)%hkey, self%scans(scan)%d(det)%psi,  psi, imod=self%npsi-1)
    call huffman_decode2(self%scans(scan)%hkey, self%scans(scan)%d(det)%flag, flag)

!!$    if (det == 1) psi = modulo(psi + 30,self%npsi)
!!$    if (det == 2) psi = modulo(psi + 20,self%npsi)
!!$    if (det == 3) psi = modulo(psi - 10,self%npsi)
!!$    if (det == 4) psi = modulo(psi - 15,self%npsi)

!!$    do j = 2, self%scans(scan)%ntod
!!$       pix(j)  = pix(j-1)  + pix(j)
!!$       psi(j)  = psi(j-1)  + psi(j)
!!$       flag(j) = flag(j-1) + flag(j)
!!$    end do
!!$    psi = modulo(psi,4096)

!!$    call int2string(scan,stext)
!!$    call int2string(det,dtext)
!!$    open(58,file='psi'//stext//'_'//dtext//'.dat')
!!$    do j = 1, self%scans(scan)%ntod
!!$       if (pix(j) == 6285034) then
!!$          write(58,*) scan, psi(j), j
!!$       end if
!!$    end do
!!$    close(58)

  end subroutine decompress_pointing_and_flags



end module comm_tod_mod
