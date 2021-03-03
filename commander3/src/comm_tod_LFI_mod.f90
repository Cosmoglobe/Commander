!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types
  use comm_huffman_mod
  use comm_hdf_mod
  use comm_fft_mod
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_4D_map_mod
  use comm_tod_driver_mod
  use comm_utils
  implicit none

  private
  public comm_LFI_tod

  type, extends(comm_tod) :: comm_LFI_tod
   contains
     procedure     :: process_tod        => process_LFI_tod
     procedure     :: simulate_LFI_tod
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info, tod_type)
    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_LFI_tod),     pointer    :: constructor

    integer(i4b) :: i, j, k, nside_beam, lmax_beam, nmaps_beam, ndelta, np_vec, ierr, offset
    real(dp)     :: f_fill, f_fill_lim(3), theta, phi, pixVal
    character(len=512) :: datadir
    logical(lgt) :: pol_beam
    type(hdf_file) :: h5_file
    integer(i4b), allocatable, dimension(:) :: pix
    real(dp), dimension(3) :: v
    character(len=10) :: fieldname
    class(tod_pointer), pointer :: selfptr

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%tod_type        = tod_type
    constructor%myid            = cpar%myid_chain
    constructor%comm            = cpar%comm_chain
    constructor%numprocs        = cpar%numprocs_chain
    constructor%myid_shared     = cpar%myid_shared
    constructor%comm_shared     = cpar%comm_shared
    constructor%myid_inter      = cpar%myid_inter
    constructor%comm_inter      = cpar%comm_inter
    constructor%info            => info
    constructor%init_from_HDF   = cpar%ds_tod_initHDF(id_abs)
    constructor%freq            = cpar%ds_label(id_abs)
    constructor%operation       = cpar%operation
    constructor%outdir          = cpar%outdir
    constructor%first_call      = .true.
    constructor%first_scan      = cpar%ds_tod_scanrange(id_abs,1)
    constructor%last_scan       = cpar%ds_tod_scanrange(id_abs,2)
    constructor%flag0           = cpar%ds_tod_flag(id_abs)
    constructor%orb_abscal      = cpar%ds_tod_orb_abscal(id_abs)
    constructor%nscan_tot       = cpar%ds_tod_tot_numscan(id_abs)
    constructor%output_4D_map   = cpar%output_4D_map_nth_iter
    constructor%output_aux_maps = cpar%output_aux_maps
    constructor%subtract_zodi   = cpar%include_TOD_zodi
    constructor%central_freq    = cpar%ds_nu_c(id_abs)
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%halfring_split  = cpar%ds_tod_halfring(id_abs)
    constructor%nside_param     = cpar%ds_nside(id_abs)
    constructor%compressed_tod  = .false.
    constructor%correct_sl      = .true.
    constructor%verbosity       = cpar%verbosity
    constructor%orb_4pi_beam    = .true.
    constructor%symm_flags      = .true.
    constructor%chisq_threshold = 7.d0

    !----------------------------------------------------------------------------------
    ! Simulation Routine
    !----------------------------------------------------------------------------------
    constructor%sims_output_dir = cpar%sims_output_dir
    constructor%enable_tod_simulations = cpar%enable_tod_simulations
    !----------------------------------------------------------------------------------

    call mpi_comm_size(cpar%comm_shared, constructor%numprocs_shared, ierr)

    if (constructor%first_scan > constructor%last_scan) then
       write(*,*) 'Error: First scan larger than last scan for ', trim(constructor%freq)
       call mpi_finalize(ierr)
       stop
    end if

    datadir = trim(cpar%datadir)//'/'
    constructor%filelist    = trim(datadir)//trim(cpar%ds_tod_filelist(id_abs))
    constructor%procmaskf1  = trim(datadir)//trim(cpar%ds_tod_procmask1(id_abs))
    constructor%procmaskf2  = trim(datadir)//trim(cpar%ds_tod_procmask2(id_abs))
    constructor%instfile    = trim(datadir)//trim(cpar%ds_tod_instfile(id_abs))

    call constructor%get_scan_ids(constructor%filelist)

    constructor%procmask => comm_map(constructor%info, constructor%procmaskf1)
    constructor%procmask2 => comm_map(constructor%info, constructor%procmaskf2)
    do i = 0, constructor%info%np-1
       if (any(constructor%procmask%map(i,:) < 0.5d0)) then
          constructor%procmask%map(i,:) = 0.d0
       else
          constructor%procmask%map(i,:) = 1.d0
       end if
       if (any(constructor%procmask2%map(i,:) < 0.5d0)) then
          constructor%procmask2%map(i,:) = 0.d0
       else
          constructor%procmask2%map(i,:) = 1.d0
       end if
    end do


    constructor%nmaps    = info%nmaps
    constructor%ndet     = num_tokens(cpar%ds_tod_dets(id_abs), ",")
    constructor%nhorn    = 1
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    allocate(constructor%label(constructor%ndet))
    allocate(constructor%partner(constructor%ndet))
    allocate(constructor%horn_id(constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do

    call constructor%read_jumplist(datadir, cpar%ds_tod_jumplist(id_abs))

    if (trim(cpar%ds_bpmodel(id_abs)) == 'additive_shift') then
       ndelta = 1
    else if (trim(cpar%ds_bpmodel(id_abs)) == 'powlaw_tilt') then
       ndelta = 1
    else
       write(*,*) 'Unknown bandpass model'
       stop
    end if
    allocate(constructor%bp_delta(0:constructor%ndet,ndelta))
    constructor%bp_delta = 0.d0

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

    ! Initialize beams
    allocate(constructor%fwhm(constructor%ndet))
    allocate(constructor%elip(constructor%ndet))
    allocate(constructor%psi_ell(constructor%ndet))
    allocate(constructor%mb_eff(constructor%ndet))
    allocate(constructor%nu_c(constructor%ndet))

    call open_hdf_file(constructor%instfile, h5_file, 'r')
    nside_beam = 512
    nmaps_beam = 3
    pol_beam   = .true.
    call read_hdf(h5_file, trim(adjustl(constructor%label(1)))//'/'//'sllmax', lmax_beam)
    constructor%slinfo => comm_mapinfo(cpar%comm_chain, nside_beam, lmax_beam, &
         & nmaps_beam, pol_beam)

    allocate(constructor%slbeam(constructor%ndet), constructor%slconv(constructor%ndet))
    allocate(constructor%mbeam(constructor%ndet))

    do i = 1, constructor%ndet
       constructor%slbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., "sl", &
            & trim(constructor%label(i)))
       constructor%mbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., "beam", trim(constructor%label(i)))
       call constructor%mbeam(i)%p%Y()
    end do

    do i = 1, constructor%ndet
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'fwhm', constructor%fwhm(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'elip', constructor%elip(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'psi_ell', constructor%psi_ell(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'mbeam_eff', constructor%mb_eff(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'centFreq', constructor%nu_c(i))
    end do

    constructor%mb_eff = 1.d0
    constructor%mb_eff = constructor%mb_eff / mean(constructor%mb_eff)
    constructor%nu_c   = constructor%nu_c * 1d9

    call close_hdf_file(h5_file)

    ! Lastly, create a vector pointing table for fast look-up for orbital dipole
    np_vec = 12*constructor%nside**2 !npix
    allocate(constructor%pix2vec(3,0:np_vec-1))
    do i = 0, np_vec-1
       call pix2vec_ring(constructor%nside, i, constructor%pix2vec(:,i))
    end do

    ! Construct observed pixel array
    allocate(constructor%pix2ind(0:12*constructor%nside**2-1))
    constructor%pix2ind = -1
    do i = 1, constructor%nscan
       allocate(pix(constructor%scans(i)%ntod))
       offset = 0
       if(constructor%halfring_split == 2) then
         offset = constructor%scans(i)%ntod
       end if
       do j = 1, constructor%ndet
          call huffman_decode(constructor%scans(i)%hkey, &
               constructor%scans(i)%d(j)%pix(1)%p, pix)
          constructor%pix2ind(pix(1)) = 1
          do k = 2, constructor%scans(i)%ntod
             pix(k)  = pix(k-1)  + pix(k)
             constructor%pix2ind(pix(k)) = 1
          end do
       end do
       deallocate(pix)
    end do
    constructor%nobs = count(constructor%pix2ind == 1)
    allocate(constructor%ind2pix(constructor%nobs))
    allocate(constructor%ind2sl(constructor%nobs))
    allocate(constructor%ind2ang(2,constructor%nobs))
    j = 1
    do i = 0, 12*constructor%nside**2-1
       if (constructor%pix2ind(i) == 1) then
          constructor%ind2pix(j) = i
          constructor%pix2ind(i) = j
          call pix2ang_ring(constructor%nside, i, theta, phi)
          call ang2pix_ring(nside_beam, theta, phi, constructor%ind2sl(j))
          constructor%ind2ang(:,j) = [theta,phi]
          j = j+1
       end if
    end do
    f_fill = constructor%nobs/(12.*constructor%nside**2)
    call mpi_reduce(f_fill, f_fill_lim(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, constructor%info%comm, ierr)
    if (constructor%myid == 0) then
       write(*,*) '  Min/mean/max TOD-map f_sky = ', real(100*f_fill_lim(1),sp), real(100*f_fill_lim(3)/constructor%info%nprocs,sp), real(100*f_fill_lim(2),sp)
    end if

    ! Set all baseline corrections to zero
    do i = 1, constructor%nscan
       constructor%scans(i)%d%baseline = 0.d0
    end do

    ! Allocate orbital dipole object
    allocate(constructor%orb_dp)
    constructor%orb_dp => comm_orbdipole(constructor%mbeam)

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
    logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp), allocatable, dimension(:)       :: procmask, procmask2
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    ! Toggle optional operations
    sample_rel_bandpass   = .true.                 
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses...
    select_data           = self%first_call        ! only perform data selection the first time
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration

    ! Initialize local variables
    ndelta          = size(delta,3)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    self%output_n_maps = 3
    if (self%output_aux_maps > 0) then
       if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
    end if

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Distribute maps
    allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
    call distribute_sky_maps(self, map_in, 1.e-6, map_sky) ! uK to K

    ! Distribute processing masks
    allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
    call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
    call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    deallocate(m_buf)

    ! Precompute far sidelobe Conviqt structures
    if (.not. self%correct_sl) then
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
    call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
    call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)

    ! Prepare intermediate data structures
    call binmap%init(self, .true., sample_rel_bandpass)
    if (sample_abs_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    ! Perform loop over scans
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_sky_prop=.true.)
       else
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
       end if
       allocate(s_buf(sd%ntod,sd%ndet))

       ! Calling Simulation Routine
       if (self%enable_tod_simulations) then
          call self%simulate_LFI_tod(i, sd%s_tot, handle)
          call sd%dealloc
          cycle
       end if

       ! Sample correlated noise
       call sample_n_corr(self, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)

       ! Compute noise spectrum parameters
       call sample_noise_psd(self, handle, i, sd%mask, sd%s_tot, sd%n_corr)

       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j))
       end do

       ! Select data
       if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute binned map
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)
       
       ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
       if (self%output_4D_map > 0) then
          if (mod(iter-1,self%output_4D_map) == 0) then
             call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
                  & samptext, self%scanid(i), self%nside, self%npsi, &
                  & self%label, self%horn_id, real(self%polang*180/pi,sp), &
                  & real(self%scans(i)%d%sigma0/self%scans(i)%d%gain,sp), &
                  & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
                  & self%scans(i)%d(:)%accept)
          end if
       end if

       ! Bin TOD
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

       ! Update scan list
       call wall_time(t1)
       self%scans(i)%proctime   = self%scans(i)%proctime   + t2-t1
       self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
       if (output_scanlist) then
          write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
               & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),&
               & real(self%spinaxis(i,:),sp)
       end if

       ! Clean up
       call sd%dealloc
       deallocate(s_buf, d_calib)

    end do

    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    call syncronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d6, chisq_S=chisq_S, mask=procmask2)
    else
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d6)
    end if
    map_out%map = binmap%outmaps(1)%p%map

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
       self%bp_delta = delta(:,:,1)
    end if
   
    ! Output maps to disk
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
    if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 7) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

    ! Clean up
    call binmap%dealloc()
    if (allocated(slist)) deallocate(slist)
    deallocate(map_sky, procmask, procmask2)
    if (self%correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
       end do
    end if

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

    call update_status(status, "tod_end"//ctext)

  end subroutine process_LFI_tod

  ! ************************************************
  !
  !> @brief Commander3 native simulation module. It
  !! simulates correlated noise and then rewrites
  !! the original timestreams inside the files.
  !
  !> @author Maksym Brilenkov
  !
  !> @param[in]
  !> @param[out]
  !
  ! ************************************************
  subroutine simulate_LFI_tod(self, scan_id, s_tot, handle)
    use omp_lib
    implicit none
    class(comm_LFI_tod), intent(inout) :: self
    ! Parameter file variables
    !type(comm_params),                     intent(in)    :: cpar
    ! Other input/output variables
    real(sp), allocatable, dimension(:,:), intent(in)    :: s_tot   !< total sky signal
    integer(i4b),                          intent(in)    :: scan_id !< current PID
    type(planck_rng),                      intent(inout) :: handle
    ! Simulation variables
    real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
    real(sp)                              :: gain   !< detector's gain value
    real(sp)                              :: sigma0
    real(sp) :: nu_knee
    real(sp) :: alpha
    real(sp) :: samprate
    real(sp) :: fft_norm
    integer(i4b)                          :: ntod !< total amount of ODs
    integer(i4b)                          :: ndet !< total amount of detectors
    ! HDF5 variables
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=3)   :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err !< MPI error status
    integer(i4b)       :: nomp !< Number of threads available
    integer(i4b)       :: omp_err !< OpenMP error status
    integer(i4b) :: omp_get_max_threads
    integer(i4b) :: n, nfft
    integer*8    :: plan_back
    real(sp) :: nu
    real(sp), allocatable, dimension(:,:) :: n_corr
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    character(len=10) :: processor_label   !< to have a nice output to screen

    ! shortcuts
    ntod = self%scans(scan_id)%ntod
    ndet = self%ndet

    ! Simulating 1/f noise
    !write(*,*) "Simulating correlated noise"
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i, j, k, dt, dv, sigma0, nu, nu_knee, alpha)
    allocate(dt(nfft), dv(0:n-1), n_corr(ntod, ndet))
    !$OMP DO SCHEDULE(guided)
    do j = 1, ndet
      ! skipping iteration if scan was not accepted
      if (.not. self%scans(scan_id)%d(j)%accept) cycle
      ! getting gain for each detector (units, V / K)
      ! (gain is assumed to be CONSTANT for EACH SCAN)
      gain   = self%scans(scan_id)%d(j)%gain
      sigma0 = self%scans(scan_id)%d(j)%sigma0
      samprate = self%samprate
      alpha    = self%scans(scan_id)%d(j)%alpha
      ! knee frequency
      nu_knee  = self%scans(scan_id)%d(j)%fknee
      ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
      fft_norm = sqrt(1.d0 * nfft)
      !
      !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      dv(0) = fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      do k = 1, (n - 1)
        nu = k * (samprate / 2) / (n - 1)
        !dv(k) = sigma0 * cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(1 + (nu / nu_knee)**alpha) /sqrt(2          .0)
        dv(k) = sigma0 * cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt((nu / nu_knee)**alpha) /sqrt(2.0)
      end do
      ! Executing Backward FFT
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      n_corr(:, j) = dt(1:ntod)
      !write(*,*) "n_corr ", n_corr(:, j)
    end do
    !$OMP END DO
    deallocate(dt, dv)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_back)

    ! Allocating main simulations' array
    allocate(tod_per_detector(ntod, ndet))       ! Simulated tod

    ! Main simulation loop
    do i = 1, ntod
      do j = 1, ndet
        ! skipping iteration if scan was not accepted
        if (.not. self%scans(scan_id)%d(j)%accept) cycle
        ! getting gain for each detector (units, V / K)
        ! (gain is assumed to be CONSTANT for EACH SCAN)
        gain   = self%scans(scan_id)%d(j)%gain
        !write(*,*) "gain ", gain
        sigma0 = self%scans(scan_id)%d(j)%sigma0
        !write(*,*) "sigma0 ", sigma0
        ! Simulating tods
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = 0
      end do
    end do

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = 'LFI_0'
    myindex = index(trim(mystring), trim(mysubstring))
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(mystring(myindex:))
    !write(*,*) "hdf5name "//trim(self%hdfname(scan_id))
    !write(*,*) "currentHDFFile "//trim(currentHDFFile)
    ! Converting PID number into string value
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    write(*,*) "Process: "//trim(processor_label)//" started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) trim(currentHDFFile)
    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop

    dims(1) = ntod
    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_error)
    ! Open an existing file - returns hdf5_file_id
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    do j = 1, ndet
      detectorLabel = self%label(j)
      ! Open an existing dataset.
      call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
      ! Write tod data to a dataset
      call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      ! Close the dataset.
      call h5dclose_f(dset_id, hdf5_error)
    end do
    ! Close the file.
    call h5fclose_f(hdf5_file_id, hdf5_error)
    ! Close FORTRAN interface.
    call h5close_f(hdf5_error)

    !write(*,*) "hdf5_error",  hdf5_error
    ! freeing memory up
    deallocate(n_corr, tod_per_detector)
    write(*,*) "Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

    ! lastly, we need to copy an existing filelist.txt into simulation folder
    ! and change the pointers to new files
    !if (self%myid == 0) then
    !  call system("cp "//trim(filelist)//" "//trim(simsdir))
    !  !mystring = filelist
    !  !mysubstring = ".txt"
    !  !myindex = index(trim(mystring), trim(mysubstring))
    !end if

    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop
  end subroutine simulate_LFI_tod

end module comm_tod_LFI_mod
