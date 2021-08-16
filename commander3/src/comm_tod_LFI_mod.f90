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
  !   Module which contains all the LFI time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates
  !       all data needed for TOD processing
  !   process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !
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
  use comm_tod_adc_mod
  implicit none

  private
  public comm_LFI_tod

  type, extends(comm_tod) :: comm_LFI_tod
     integer(i4b) :: nbin_spike
     integer(i4b) :: nbin_adc
     logical(lgt) :: use_dpc_adc
     real(dp),          allocatable, dimension(:)       :: mb_eff
     real(dp),          allocatable, dimension(:,:)     :: diode_weights
     type(spline_type), allocatable, dimension(:,:)     :: ref_splint ! ndet
     type(adc_pointer), allocatable, dimension(:,:)     :: adc_corrections ! ndet, n_diode
     real(dp),          allocatable, dimension(:,:)     :: spike_templates ! nbin, ndet
     real(dp),          allocatable, dimension(:,:)     :: spike_amplitude ! nscan, ndet
     real(dp),          allocatable, dimension(:,:,:)   :: R               ! nscan, ndet, ndiode/2
     character(len=10)                                  :: adc_mode        ! gauss, dpc, none
   contains
     procedure     :: process_tod             => process_LFI_tod
     procedure     :: diode2tod_inst          => diode2tod_LFI
     procedure     :: load_instrument_inst    => load_instrument_LFI
     procedure     :: dumpToHDF_inst          => dumpToHDF_LFI
     procedure     :: construct_corrtemp_inst => construct_corrtemp_LFI
     procedure     :: filter_reference_load
     procedure     :: compute_ref_load_filter
     procedure     :: get_nsmooth
     procedure     :: preprocess_L1_to_L2
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(handle, cpar, id_abs, info, tod_type)
    !
    ! Constructor function that gathers all the instrument parameters in a pointer
    ! and constructs the objects
    !
    ! Arguments:
    ! ----------
    ! handle:   type(planck_rng)
    !           Healpix random number type
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    ! id_abs:   integer
    !           The index of the current band within the parameters, related to cpar
    ! info:     map_info structure
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data

    implicit none
    type(planck_rng),          intent(inout) :: handle
    type(comm_params),         intent(in)    :: cpar
    integer(i4b),              intent(in)    :: id_abs
    class(comm_mapinfo),       target        :: info
    character(len=128),        intent(in)    :: tod_type
    class(comm_LFI_tod),       pointer       :: constructor

    real(sp), dimension(:,:),    allocatable :: diode_data, corrected_data
    integer(i4b), dimension(:),  allocatable :: flag

    integer(i4b) :: i, j, k, nside_beam, lmax_beam, nmaps_beam, ierr, filter_count, nsmooth
    logical(lgt) :: pol_beam
    character(len=50) :: name

    real(dp), dimension(:),   allocatable :: nus
    real(sp), dimension(:,:), allocatable :: filtered
    real(dp), dimension(:),   allocatable :: nu_saved
    real(dp), dimension(:,:), allocatable :: filter_sum

    ! Allocate object
    allocate(constructor)

    ! Set up noise PSD type and priors
    constructor%freq            = cpar%ds_label(id_abs)    
    if (trim(constructor%freq) == '030') then
       constructor%n_xi            = 6
       constructor%noise_psd_model = 'oof_gauss'    
       allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
       allocate(constructor%xi_n_P_rms(constructor%n_xi))
       constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0, 1.d6, 0.d0, 0.d0] ! [sigma0, fknee, alpha, g_amp, g_loc, g_sig]; sigma0 is not used

       constructor%xi_n_nu_fit     = [0.d0, 3*1.225d0]    ! More than max(7*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.010d0, 0.45d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-2.5d0, -0.4d0]   ! alpha
       constructor%xi_n_P_uni(4,:) = [0.0d0,   1d0]     ! g_amp
       constructor%xi_n_P_uni(5,:) = [1.35d0,  1.35d0 ] ! g_loc
       constructor%xi_n_P_uni(6,:) = [0.4d0,   0.4d0]   ! g_sig
    else if (trim(constructor%freq) == '044') then
       constructor%n_xi            = 6
       constructor%noise_psd_model = 'oof_gauss'
       allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
       allocate(constructor%xi_n_P_rms(constructor%n_xi))
       constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0, 1.d6, 0.d0, 0.d0] ! [sigma0, fknee, alpha, g_amp, g_loc, g_sig]; sigma0 is not used

       constructor%xi_n_nu_fit     = [0.d0, 3*1.00d0]    ! More than max(2*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.002d0, 0.40d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-2.5d0, -0.4d0]   ! alpha
       constructor%xi_n_P_uni(4,:) = [0.0d0,   1d0]     ! g_amp
       constructor%xi_n_P_uni(5,:) = [1.35d0,  1.35d0 ] ! g_loc
       constructor%xi_n_P_uni(6,:) = [0.4d0,   0.4d0]   ! g_sig
    else if (trim(constructor%freq) == '070') then
       constructor%n_xi            = 3
       constructor%noise_psd_model = 'oof'
       allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
       allocate(constructor%xi_n_P_rms(constructor%n_xi))
       constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0] ! [sigma0, fknee, alpha]; sigma0 is not used
       constructor%xi_n_nu_fit     = [0.d0, 0.140d0]    ! More than max(2*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.001d0, 0.25d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-3.0d0, -0.4d0]   ! alpha
    else
       write(*,*) 'Invalid LFI frequency label = ', trim(constructor%freq)
       stop
    end if

    ! Initialize instrument-specific parameters
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn           = 1
    constructor%ndiode          = 4
    constructor%sample_L1_par   = .true. !.false.
    constructor%compressed_tod  = .true.
    constructor%correct_sl      = .true.
    constructor%orb_4pi_beam    = .true.
    constructor%use_dpc_adc     = .false. !.true.
    constructor%symm_flags      = .true.
    constructor%chisq_threshold = 20.d6 ! 9.d0
    constructor%nmaps           = info%nmaps
    constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")

    nside_beam                  = 512
    nmaps_beam                  = 3
    pol_beam                    = .true.
    constructor%nside_beam      = nside_beam

    ! Initialize common parameters
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)

    ! Get detector labels
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
    
    ! Define detector partners
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do

    ! Define diode labels
    do i = 1, constructor%ndet
       if (index(constructor%label(i), 'M') /= 0) then
          constructor%diode_names(i,1) = 'ref00'
          constructor%diode_names(i,2) = 'sky00'
          constructor%diode_names(i,3) = 'ref01'
          constructor%diode_names(i,4) = 'sky01'
       else
          constructor%diode_names(i,1) = 'ref10'
          constructor%diode_names(i,2) = 'sky10'
          constructor%diode_names(i,3) = 'ref11'
          constructor%diode_names(i,4) = 'sky11'
       end if
    end do

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Setting polarization angles to DPC post-analysis values
    if (trim(constructor%freq) == '030') then
       constructor%polang = -[-3.428, -3.428, 2.643, 2.643]*pi/180.
    else if (trim(constructor%freq) == '044') then
       constructor%polang = -[-2.180, -2.180,  7.976, 7.976, -4.024, -4.024]*pi/180.
    else if (trim(constructor%freq) == '070') then
       constructor%polang = -[ 0.543, 0.543,  1.366, 1.366,  -1.811, -1.811, -1.045, -1.045,  -2.152, -2.152,  -0.960, -0.960]*pi/180.
    end if

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    call constructor%precompute_lookups()

    ! allocate LFI specific instrument file data
    constructor%nbin_spike      = nint(constructor%samprate*sqrt(3.d0))
    allocate(constructor%mb_eff(constructor%ndet))
    allocate(constructor%diode_weights(constructor%ndet, 2))
    allocate(constructor%spike_templates(0:constructor%nbin_spike-1, constructor%ndet))
    allocate(constructor%spike_amplitude(constructor%nscan,constructor%ndet))
    allocate(constructor%adc_corrections(constructor%ndet, constructor%ndiode))
    allocate(constructor%ref_splint(constructor%ndet,constructor%ndiode/2))
    allocate(constructor%R(constructor%nscan,constructor%ndet,constructor%ndiode/2))

    ! Declare adc_mode 
    constructor%adc_mode = 'gauss'
    constructor%nbin_adc = 500

    ! Load the instrument file
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)
    constructor%spike_amplitude = 0.d0

    ! Compute ADC correction tables for each diode
    if (.not. constructor%L2_exist) then

       if (.not. constructor%use_dpc_adc) then
          if (constructor%myid == 0) write(*,*) '   Building ADC correction tables'

          ! Determine v_min and v_max for each diode
          call update_status(status, "ADC_start")
          do i = 1, constructor%ndet

             do j=1, constructor%ndiode ! init the adc correction structures
                constructor%adc_corrections(i,j)%p => comm_adc(cpar,info,constructor%nbin_adc)
             end do

             do k = 1, constructor%nscan ! determine vmin and vmax for each diode
                if (.not. constructor%scans(k)%d(i)%accept) cycle
                allocate(diode_data(constructor%scans(k)%ntod, constructor%ndiode))
                allocate(flag(constructor%scans(k)%ntod))
                call constructor%decompress_diodes(k, i, diode_data, flag=flag)
                do j = 1, constructor%ndiode
                   call constructor%adc_corrections(i,j)%p%find_horn_min_max(diode_data(:,j), flag,constructor%flag0)
                end do
                deallocate(diode_data, flag)

             end do ! end loop over scans

             do j = 1, constructor%ndiode ! allreduce vmin and vmax
                ! All reduce min and max
                call mpi_allreduce(mpi_in_place,constructor%adc_corrections(i,j)%p%v_min,1,MPI_REAL,MPI_MIN,constructor%comm,ierr)
                call mpi_allreduce(mpi_in_place,constructor%adc_corrections(i,j)%p%v_max,1,MPI_REAL,MPI_MAX,constructor%comm,ierr)
                call constructor%adc_corrections(i,j)%p%construct_voltage_bins
             end do
          end do
          call update_status(status, "ADC_range")

          ! Now bin rms for all scans and compute the correction table
          if (constructor%myid == 0) write(*,*) '    Bin RMS for ADC corrections'
          do k = 1, constructor%nscan ! compute and bin the rms as a function of voltage for each scan
             allocate(diode_data(constructor%scans(k)%ntod, constructor%ndiode))
             allocate(flag(constructor%scans(k)%ntod))
             do i = 1, constructor%ndet
                if (.not. constructor%scans(k)%d(i)%accept) cycle
                call constructor%decompress_diodes(k, i, diode_data, flag)
                do j = 1, constructor%ndiode
                   call constructor%adc_corrections(i,j)%p%bin_scan_rms(diode_data(:,j), flag,constructor%flag0) 
                end do
             end do
             deallocate(diode_data, flag)
          end do
          call update_status(status, "ADC_bin")

          if (constructor%myid == 0) write(*,*) '    Generate ADC correction tables'
          do i = 1, constructor%ndet
             do j = 1, constructor%ndiode
                ! Build the actual adc correction tables (adc_in, adc_out)
                name = trim(constructor%label(i))//'_'//trim(constructor%diode_names(i,j))
                if (constructor%myid == 0) write(*,*) '    Building table for '// trim(name)
                call constructor%adc_corrections(i,j)%p%build_table(handle, name)
             end do
          end do
          call update_status(status, "ADC_table")
       end if

       ! Compute reference load filter spline
       if (constructor%myid == 0) write(*,*) '   Build reference load filter'
       nsmooth = constructor%get_nsmooth()
       allocate(filter_sum(constructor%ndiode/2,nsmooth))
       allocate(nu_saved(nsmooth))
       nu_saved = 0.d0
       do i=1, constructor%ndet
          filter_count = 0
          filter_sum   = 0.d0
          do k = 1, constructor%nscan
             if (.not. constructor%scans(k)%d(i)%accept) cycle

             allocate(diode_data(constructor%scans(k)%ntod, constructor%ndiode), corrected_data(constructor%scans(k)%ntod, constructor%ndiode))
             call constructor%decompress_diodes(k, i, diode_data)

             if (any(diode_data == 0.)) then
                write(*,*) 'Contains zeros', constructor%scanid(k), i
                constructor%scans(k)%d(i)%accept = .false.
                deallocate(diode_data, corrected_data)
                cycle
             end if

             ! corrected_data = diode_data
             do j = 1, constructor%ndiode
                call constructor%adc_corrections(i,j)%p%adc_correct(diode_data(:,j), corrected_data(:,j), constructor%scanid(k),i,j)
             end do
             if (any(abs(corrected_data(:,[1,3])) > 10)) then
                constructor%scans(k)%d(i)%accept = .false.
                deallocate(diode_data, corrected_data)
                cycle
             end if

             ! compute the ref load transfer function
             call constructor%compute_ref_load_filter(corrected_data, filter_sum, nu_saved, ierr)
             if (ierr == 0) filter_count = filter_count + 1
             
             deallocate(diode_data, corrected_data)




          end do

             ! Mpi average the load filter over all cores, save as a spline
       call mpi_allreduce(MPI_IN_PLACE, filter_count, 1, MPI_INTEGER, MPI_SUM, constructor%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, filter_sum, size(filter_sum), MPI_DOUBLE_PRECISION, MPI_SUM, constructor%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, nu_saved,   size(nu_saved), MPI_DOUBLE_PRECISION, MPI_MAX, constructor%info%comm, ierr)

       ! Remove empty bins
       j = 1
       do while (j <= nsmooth)
          if (filter_sum(1,j) == 0.d0) then
             filter_sum(:,j+1:nsmooth) = filter_sum(:,j:nsmooth-1)
             nu_saved(j+1:nsmooth)     = nu_saved(j:nsmooth-1)
             nsmooth                   = nsmooth-1
          else
             j = j+1
          end if
       end do

       ! HKE: Should probably manually add a regularization bin at the end, so that the spline doesn't go crazy for frequencies between the last bin center and fsamp/2
       if (filter_count > 0) then
          filter_sum = filter_sum/filter_count
       else
          filter_sum = 1.d0
          do j = 1, size(nu_saved)
             nu_saved(j) = constructor%samprate/2 * (j-1)/real(size(nu_saved)-1,dp)
          end do
       end if
       !if (constructor%myid == 0) write(*,*) nu_saved
       do j=1, constructor%ndiode/2
          call spline_simple(constructor%ref_splint(i,j), nu_saved, filter_sum(j,:))
         !if (constructor%myid == 0) then
        !  open(100, file=trim(constructor%outdir)//'/load_filter_'//trim(constructor%label(i))//'_'//trim(constructor%diode_names(i,j+2))//'.dat')
        !  do k = 1, size(nu_saved)
        !    write(100, fmt='(f30.8,f30.8)') nu_saved(k), filter_sum(j,k)
        !  end do
        !  close(100)
        !  write(*,*) 'Writing file ', trim(constructor%outdir)//'/load_filter_'//trim(constructor%label(i))//'_'//trim(constructor%diode_names(i,j))//'.dat'
        !end if
       end do

       end do
       call update_status(status, "ADC_table")

    end if


    ! Allocate sidelobe convolution data structures
    allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
    constructor%orb_dp => comm_orbdipole(constructor%mbeam)

    ! Initialize all baseline corrections to zero
    do i = 1, constructor%nscan
       constructor%scans(i)%d%baseline = 0.d0
    end do

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    !
    ! Routine that processes the LFI time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_LFI_tod class
    !           Points to output of the constructor
    ! chaindir: string
    !           Directory for output files
    ! chain:    integer
    !           Index number of the chain being run
    ! iter:     integer
    !           Gibbs iteration number
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same point
    ! map_in:   array
    !           Array of dimension (ndet,ndelta) with pointer to maps,
    !           with both access to maps and changing them.
    !           ndet is the number of detectors and
    !           ndelta is the number of bandpass deltas being considered
    ! delta:    array
    !           Array of bandpass corrections with dimensions (0:ndet,npar,ndelta)
    !           where ndet is number of detectors, npar is number of parameters
    !           and ndelta is the number of bandpass deltas being considered
    !
    ! Returns:
    ! ----------
    ! map_out: comm_map class
    !          Final output map after TOD processing combined for all detectors
    ! rms_out: comm_map class
    !          Final output rms map after TOD processing combined for all detectors

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
    character(len=512)  :: prefix, postfix, prefix4D, filename, Sfilename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, sigma0
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1 .and. mod(iter,2) == 0     ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = size(delta,3) > 1 .and. mod(iter,2) == 1     ! don't sample absolute bandpasses
    select_data           = self%first_call        ! only perform data selection the first time
    output_scanlist       = mod(iter-1,1) == 0    ! only output scanlist every 10th iteration

    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta-1
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
    if (self%correct_sl) then
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if
!    write(*,*) 'qqq', self%myid
!    if (.true. .or. self%myid == 78) write(*,*) 'a', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
!!$    call mpi_finalize(ierr)
!!$    stop

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Pre-process L1 data into L2 data if requested, and set ndiode = 1 to skip directly to L2 later on
    if (.not. self%sample_L1_par .and. self%ndiode > 1) then
       call self%preprocess_L1_to_L2(map_sky, procmask)
       self%ndiode = 1
       self%compressed_tod = .false.
    end if
    call update_status(status, "L1_to_L2")


    ! Sample 1Hz spikes
    call sample_1Hz_spikes(self, handle, map_sky, procmask, procmask2); call update_status(status, "tod_1Hz")

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2); call update_status(status, "tod_gain1")
    call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2); call update_status(status, "tod_gain2")
    call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2); call update_status(status, "tod_gain3")
    call sample_gain_psd(self, handle)

    ! Prepare intermediate data structures
    call binmap%init(self, .true., sample_rel_bandpass)
    if (sample_abs_bandpass .or. sample_rel_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    ! Perform loop over scans
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
!          if (.true. .or. self%myid == 78) write(*,*) 'b', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_sky_prop=.true.)
       else
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
       end if
       allocate(s_buf(sd%ntod,sd%ndet))

       ! Calling Simulation Routine
       if (self%enable_tod_simulations) then
          call simulate_tod(self, i, sd%s_tot, handle)
          call sd%dealloc
          cycle
       end if

       ! Sample correlated noise
       call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)

       ! Compute noise spectrum parameters
       call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
       end do

       ! Select data
       if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute binned map
       allocate(d_calib(binmap%nout,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)
       
       ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
       if (self%output_4D_map > 0) then
          if (mod(iter-1,self%output_4D_map) == 0) then
             allocate(sigma0(sd%ndet))
             do j = 1, sd%ndet
                sigma0(j) = self%scans(i)%d(j)%N_psd%sigma0/self%scans(i)%d(j)%gain
             end do
             call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
                  & samptext, self%scanid(i), self%nside, self%npsi, &
                  & self%label, self%horn_id, real(self%polang*180/pi,sp), sigma0, &
                  & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
                  & self%scans(i)%d(:)%accept)
             deallocate(sigma0)
          end if
       end if

       ! Bin TOD
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

       ! Update scan list
       call wall_time(t2)
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

    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    call synchronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       Sfilename = trim(prefix) // 'Smap'// trim(postfix)
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d6, chisq_S=chisq_S, &
            & Sfilename=Sfilename, mask=procmask2)
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
    if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

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
  
  
  subroutine load_instrument_LFI(self, instfile, band)
    !
    ! Reads the LFI specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_LFI_tod
    !    the LFI tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_LFI_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band

    integer(i4b) :: i, j
    real(dp) :: weight
    character(len=1) :: id
    character(len=512) :: path

    if(index(self%label(band), 'M') /= 0) then
       self%diode_names(band,:) = ['ref00','sky00','ref01','sky01']
       id = '0'
    else
       self%diode_names(band,:) = ['ref10','sky10','ref11','sky11']
       id = '1'
    end if

    ! Read in mainbeam_eff
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'mbeam_eff', self%mb_eff(band))

    ! read in the diode weights
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'diodeWeight', weight)
    self%diode_weights(band,1) = weight
    self%diode_weights(band,2) = 1.d0 - weight

    if (self%use_dpc_adc .and. .not. self%L2_exist) then
       ! Read ADC corrections
       path = trim(adjustl(self%label(band)))//'/'//'adc91-'//id//'0'
       self%adc_corrections(band,1)%p => comm_adc(instfile, path, .true.) !first half, load
       self%adc_corrections(band,2)%p => comm_adc(instfile, path, .false.) !first half, load
       path = trim(adjustl(self%label(band)))//'/'//'adc91-'//id//'1'
       self%adc_corrections(band,3)%p => comm_adc(instfile, path, .true.) !first half, load
       self%adc_corrections(band,4)%p => comm_adc(instfile, path, .false.) !first half, load
    end if


  end subroutine load_instrument_LFI
  
  subroutine diode2tod_LFI(self, scan, map_sky, procmask, tod)
    ! 
    ! Generates detector-coadded TOD from low-level diode data
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! scan:     int
    !           Scan ID number
    ! procmask: array of sp
    !           processing mask that cuts out the galaxy
    !
    ! Returns
    ! ----------
    ! tod:      ntod x ndet sp array
    !           Output detector TOD generated from raw diode data
    !
    implicit none
    class(comm_LFI_tod),                       intent(inout) :: self
    integer(i4b),                              intent(in)    :: scan
    real(sp),          dimension(0:,1:,1:,1:), intent(in)    :: map_sky
    real(sp),          dimension(0:),          intent(in)    :: procmask
    real(sp),          dimension(:,:),         intent(out)   :: tod

    integer(i4b) :: i,j,k,half,horn,n_mask, n_unmask
    real(sp), allocatable, dimension(:,:) :: diode_data, corrected_data, s_sky, mask
    integer(i4b), allocatable, dimension(:,:,:) :: pix, psi
    integer(i4b), allocatable, dimension(:,:)   :: flag
    real(dp) :: r1, r2, sum1, sum2, A(3,3,2), b(3,2), x(3,2)

    allocate(diode_data(self%scans(scan)%ntod, self%ndiode))
    allocate(corrected_data(self%scans(scan)%ntod, self%ndiode))
    allocate(pix(self%scans(scan)%ntod,self%ndet,self%nhorn), psi(self%scans(scan)%ntod,self%ndet,self%nhorn), flag(self%scans(scan)%ntod,self%ndet))
    allocate(s_sky(self%scans(scan)%ntod, self%ndet), mask(self%scans(scan)%ntod, self%ndet))

    do i=1, self%ndet
      call self%decompress_pointing_and_flags(scan, i, pix(:,i,:), psi(:,i,:), flag(:,i))
    end do
    call project_sky(self, map_sky(:,:,:,1), pix(:,:,1), psi(:,:,1), flag, &
         & procmask, scan, s_sky, mask)

    do i=1, self%ndet

       if (.not. self%scans(scan)%d(i)%accept) cycle

        ! Decompress diode TOD for current scan
        call self%decompress_diodes(scan, i, diode_data)

        ! Apply ADC corrections
        do j=1, self%ndiode
          call self%adc_corrections(i,j)%p%adc_correct(diode_data(:,j), corrected_data(:,j))
          !corrected_data(:,j) = diode_data(:,j)
          !do k = 1, 10
          !   write(*,*) diode_data(k,j), corrected_data(k,j)
          !end do
          !stop

          !corrected_data(:,j) = diode_data(:,j)
        end do

        ! Wiener-filter load data         
        call self%filter_reference_load(i, corrected_data)

        ! Compute the gain modulation factors

!!$        r1 = 0.d0
!!$        r2 = 0.d0
!!$        sum1 = 0.d0
!!$        sum2 = 0.d0
!!$        n_mask = 0
!!$        n_unmask = 0
!!$
!!$        do k = 1, size(corrected_data(:,1))
!!$          if (iand(flag(k), self%flag0) .ne. 0) cycle
!!$
!!$          sum1 = sum1 + corrected_data(k,1)
!!$          sum2 = sum2 + corrected_data(k,3)
!!$          n_unmask = n_unmask + 1
!!$
!!$          if (procmask(pix(k)) .ne. 0) then 
!!$            r1 = r1 + corrected_data(k,2)
!!$            r2 = r2 + corrected_data(k,4)
!!$            n_mask = n_mask + 1
!!$          end if
!!$
!!$        end do
!!$
!!$        if (r1 == 0.d0 .or. r2 == 0.d0 .or. sum1 == 0.d0 .or. sum2 == 0.d0) then
!!$           self%scans(scan)%d(i)%accept = .false.
!!$           cycle
!!$        end if

        A = 0.d0
        b = 0.d0
        do k = 1, size(corrected_data(:,1))
          if (mask(k,i) == 0.) cycle

          A(1,1,1) = A(1,1,1) + 1.d0
          A(1,2,1) = A(1,2,1) + s_sky(k,i)
          A(1,3,1) = A(1,3,1) + corrected_data(k,1)
          A(2,2,1) = A(2,2,1) + s_sky(k,i)          * s_sky(k,i)
          A(2,3,1) = A(2,3,1) + s_sky(k,i)          * corrected_data(k,1)
          A(3,3,1) = A(3,3,1) + corrected_data(k,1) * corrected_data(k,1)
          b(1,1)   = b(1,1)   +                       corrected_data(k,2)
          b(2,1)   = b(2,1)   + s_sky(k,i)          * corrected_data(k,2)
          b(3,1)   = b(3,1)   + corrected_data(k,1) * corrected_data(k,2)

          A(1,1,2) = A(1,1,2) + 1.d0
          A(1,2,2) = A(1,2,2) + s_sky(k,i)
          A(1,3,2) = A(1,3,2) + corrected_data(k,3)
          A(2,2,2) = A(2,2,2) + s_sky(k,i)          * s_sky(k,i)
          A(2,3,2) = A(2,3,2) + s_sky(k,i)          * corrected_data(k,3)
          A(3,3,2) = A(3,3,2) + corrected_data(k,3) * corrected_data(k,3)
          b(1,2)   = b(1,2)   +                       corrected_data(k,4)
          b(2,2)   = b(2,2)   + s_sky(k,i)          * corrected_data(k,4)
          b(3,2)   = b(3,2)   + corrected_data(k,3) * corrected_data(k,4)
       end do
       if (A(1,1,1) == 0.d0) then
           self%scans(scan)%d(i)%accept = .false.
           cycle
        end if        
        do j = 1, 3
           do k = j+1, 3
              A(k,j,:) = A(j,k,:)
           end do
        end do
        call solve_system_real(A(:,:,1), x(:,1), b(:,1))
        call solve_system_real(A(:,:,2), x(:,2), b(:,2))
        
        ! average sky value/average load value
        self%R(scan,i,1) = x(3,1) !(r1/n_mask)/(sum1/n_unmask)
        self%R(scan,i,2) = x(3,2) !(r2/n_mask)/(sum2/n_unmask)
        !write(*,*) self%scanid(scan), i, real(self%R(scan,i,:),sp)

        ! Compute output differenced TOD

        !w1(sky00 - R*ref00) + w2(sky01 - R*ref01)
        !if(self%myid == 0) write(*,*) r1, r2, n_mask, size(diode_data(:,1))
        !tod(:,i) = self%diode_weights(i,1) * (diode_data(:,2) - self%R(scan,i,1) * diode_data(:,1)) + self%diode_weights(i,2)*(diode_data(:,4) - self%R(scan,i,2) * diode_data(:,3))
!        do k = 0, 100
!           self%R(scan,i,2) = 0.70 + 0.003*k
           tod(:,i) = self%diode_weights(i,1) * (corrected_data(:,2) - self%R(scan,i,1) * corrected_data(:,1)) + self%diode_weights(i,2)*( corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3))
!           tod(:,i) = corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3)
!           write(*,*) self%R(scan,i,2), sum(tod(:,i)**2), variance(1.d0*tod(:,i))
!        end do
        !stop

!!$    if (self%scanid(scan) == 3 .and. i == 1) then
!!$        open(58,file='res_L2fromL1_030_pid3.dat', recl=1024)
!!$        do j = 1, size(tod,1)
!!$           write(58,*) j, corrected_data(j,2) - x(1,1) - s_sky(j,1)*x(2,1) - corrected_data(j,1)*x(3,1)!, corrected_data(j,:)
!!$        end do
!!$        close(58)
!!$     end if


        !tod(:,i) = self%diode_weights(i,1) * (corrected_data(:,2) - filtered_data(:,1)) + self%diode_weights(i,2)*( corrected_data(:,4) - filtered_data(:,3))
        !tod(:,i) = (corrected_data(:,1) - corrected_data(:,3)) + (corrected_data(:,2) - corrected_data(:,4))
        

        !stop
        
    end do
!    stop

    if (self%scanid(scan) == 3) then
        open(58,file='comm3_L2fromL1_030_pid3.dat', recl=1024)
        do j = 1, size(tod,1)
           write(58,*) j, tod(j,:), diode_data(j,:)!, corrected_data(j,:)
        end do
        close(58)
     end if

    deallocate(diode_data, corrected_data, pix, psi, flag, s_sky, mask)

!stop

!call mpi_finalize(i)
!stop

  end subroutine diode2tod_LFI

  function get_nsmooth(self)
    implicit none
    class(comm_LFI_tod),  intent(in)   :: self
    integer(i4b)                       :: get_nsmooth  
    integer(i4b) :: j
    real(sp)     :: fbin, nu

    fbin         = 1.2 ! multiplicative bin scaling factor
    get_nsmooth  = 2
    nu           = 0.01
    do while (nu <= self%samprate/2)
       get_nsmooth = get_nsmooth + 1
       nu          = nu * fbin
    end do
  end function get_nsmooth

  subroutine compute_ref_load_filter(self, data_in, binned_out, nu_out, err)
    ! 
    ! Computes the binned weiner filter for the reference load
    !
    ! Arguments:
    ! ----------
    ! 
    ! self:     comm_tod_LFI object
    !           TOD processing class
    ! data_in:  float array (ntod, ndiode)
    !           input diode timestreams
    !
    ! Returns:
    ! --------
    !
    ! binned_out : float array
    !              array of filter transfer function for ref load
    ! nu_out     : float_array
    !              frequencies that index binned_out
    ! err        : error flag; 0 if OK, 1 if no data
    implicit none
    class(comm_LFI_tod),          intent(in)    :: self
    real(sp),     dimension(:,:), intent(in)    :: data_in
    real(dp),     dimension(:,:), intent(inout) :: binned_out
    real(dp),     dimension(:),   intent(inout) :: nu_out
    integer(i4b),                 intent(out)   :: err

    integer(i4b) :: i, j, k, nfft, n, nsmooth
    real(dp)     :: num, denom, fsamp, fbin, nu, upper, subsum, nu_low, delta_nu
    integer*8    :: plan_fwd

    real(sp),     allocatable, dimension(:) :: dt_sky, dt_ref
    real(dp),     allocatable, dimension(:) :: filter
    complex(spc), allocatable, dimension(:) :: dv_sky, dv_ref

    ! This test should be replaced with something more fine-tuned
    if (all(data_in == 0.)) then
       err = 1
       return
    else
       err = 0
    end if

    n       = size(data_in(:,1))
    nfft    = n/2+1
    fsamp   = self%samprate
    nsmooth = self%get_nsmooth()

    allocate(dt_sky(n), dt_ref(n), dv_sky(0:nfft-1), dv_ref(0:nfft-1), filter(0:nfft-1))
    
    call sfftw_plan_dft_r2c_1d(plan_fwd, n, dt_ref, dv_ref, fftw_estimate + fftw_unaligned)

    do i = 1, self%ndiode/2

      ! Check if data is all zeros
      dt_ref = data_in(:, 2*i -1) 
      dt_sky = data_in(:, 2*i)
      
      ! FFT of ref signal
      call sfftw_execute_dft_r2c(plan_fwd, dt_ref, dv_ref)

      ! FFT of sky signal
      call sfftw_execute_dft_r2c(plan_fwd, dt_sky, dv_sky)     

      ! Compute cross correlation
      do j = 0, nfft-1
         num   =      real(dv_sky(j)*conjg(dv_ref(j)) + dv_ref(j)*conjg(dv_sky(j)),dp)
         denom = sqrt(real(dv_sky(j)*conjg(dv_sky(j)) * dv_ref(j)*conjg(dv_ref(j)),dp))
         !write(*,*) j, num, denom
         if (denom < 1d-100) then
            filter(j) = 0.
         else 
            filter(j) = 0.5*num/denom
         end if
      end do

      ! Set first bin to unity
      nu_low          = 0.01d0 ! Lower limit in Hz; crosscorr is defined to be 1 below this
      binned_out(i,1) = binned_out(i,1) + 1.d0
      delta_nu        = ind2freq(2, fsamp, nfft)
      nu_out(1)       = sqrt(delta_nu * nu_low)

      ! Bin with logarithmic bin width above nu_low
      fbin         = 1.2                      ! Bin scaling factor
      j            = nint(0.01d0/delta_nu)    ! First frequency to consider
      nu           = ind2freq(j, fsamp, nfft) ! Current frequency
      upper        = nu                       ! Lower limit of first bin; init
      nsmooth      = 2                        ! Bin counter
      do while (nu < fsamp/2)
         upper           = min(upper*fbin, fsamp/2)
         nu_out(nsmooth) = nu    ! Start of bin
         subsum          = 0.d0  ! Summing variable
         k               = 0     ! Number of frequencies in current bin
         do while (nu <= upper)
            subsum = subsum + filter(j-1)
            k      = k+1
            j      = j+1
            nu     = nu + delta_nu
         end do
         !write(*,*) nsmooth, nu_out(nsmooth), nu, fsamp/2, sqrt(nu_out(nsmooth)*nu)
         nu_out(nsmooth) = sqrt(nu_out(nsmooth) * nu)
         if (k > 0) binned_out(i, nsmooth) = binned_out(i, nsmooth) + subsum/k
         nsmooth = nsmooth+1
      end do

    end do

    !if(self%myid == 0) write(*,*) binned_out, size(binned_out), size(nu_out)

    call sfftw_destroy_plan(plan_fwd)

    deallocate(dt_sky, dt_ref, dv_sky, dv_ref, filter)

  end subroutine compute_ref_load_filter


  subroutine filter_reference_load(self, det, data)
    class(comm_LFI_tod),               intent(in)      :: self
    integer(i4b),                      intent(in)      :: det
    real(sp), dimension(:,:),          intent(inout)   :: data

    integer(i4b) :: i, j, nfft, n
    integer*8    :: plan_fwd, plan_back

    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv

    n       = size(data(:,1))
    nfft    = n/2+1

    allocate(dt(n), dv(0:nfft-1))

    call sfftw_plan_dft_r2c_1d(plan_fwd,  n, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, n, dv, dt, fftw_estimate + fftw_unaligned)

!!$    open(58,file='raw.dat')
!!$    do i = 1, n
!!$       write(58,*) data(i,:)
!!$    end do
!!$    close(58)

    do i = 1, self%ndiode/2

      ! Check if data is all zeros
      dt = data(:, 2*i -1)
      if(all(dt == 0)) cycle

      ! FFT of ref signal
      call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

      ! Filter ref with cross correlation transfer function
!      open(58,file='filter.dat')
      do j=1, size(dv) -1
        dv(j) = dv(j) * splint(self%ref_splint(det,i), ind2freq(j, self%samprate, nfft))
!        write(58,*) ind2freq(j, self%samprate, nfft), splint(self%ref_splint(i), ind2freq(j, self%samprate, nfft))
      end do
!     close(58)

      ! IFFT ref signal
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      
      ! Normalize
      data(:, 2*i-1) = dt/n

    end do

    call sfftw_destroy_plan(plan_fwd)
    call sfftw_destroy_plan(plan_back)

    deallocate(dt, dv)

!!$    open(58,file='filtered.dat')
!!$    do i = 1, n
!!$       write(58,*) data(i,:)
!!$    end do
!!$    close(58)
!!$    stop

  end subroutine filter_reference_load

  subroutine dumpToHDF_LFI(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_LFI_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    character(len=10) :: diode_name
    integer(i4b) :: ierr, i, j
    real(dp), allocatable, dimension(:,:)   :: amp, amp_tot
    real(dp), allocatable, dimension(:,:,:) :: R, R_tot
    real(dp), allocatable, dimension(:,:,:,:) :: ref_filter, adc_corr

    allocate(amp(self%nscan_tot,self%ndet), amp_tot(self%nscan_tot,self%ndet))
    amp = 0.d0
    amp(self%scanid,:) = self%spike_amplitude
    call mpi_reduce(amp, amp_tot, size(amp), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

    allocate(R(self%nscan_tot,self%ndet,size(self%R,3)),R_tot(self%nscan_tot,self%ndet,size(self%R,3)))
    R = 0.d0
    R(self%scanid,:,:) = self%R
    call mpi_reduce(R, R_tot, size(R), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

    if (self%myid == 0) then
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_temp', self%spike_templates)
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_ampl', amp_tot)
       call write_hdf(chainfile, trim(adjustl(path))//'R_factor', R_tot)
       call write_hdf(chainfile, trim(adjustl(path))//'w_diode', self%diode_weights)

       if (allocated(self%ref_splint(1,1)%x)) then
          allocate(ref_filter(size(self%ref_splint(1,1)%y),2,self%ndet,size(self%ref_splint(1,:))))
          do i = 1, self%ndet
             do j = 1, size(self%ref_splint(1,:))
                ref_filter(:,1,i,j) = self%ref_splint(i,j)%x
                ref_filter(:,2,i,j) = self%ref_splint(i,j)%y
             end do
          end do
          call write_hdf(chainfile, trim(adjustl(path))//'ref_filter', ref_filter)
          deallocate(ref_filter)
       end if

       if (associated(self%adc_corrections(1,1)%p)) then
          allocate(adc_corr(size(self%adc_corrections(1,1)%p%adc_in),2,self%ndet,size(self%adc_corrections(1,:))))
          do i = 1, self%ndet
             do j = 1, size(self%adc_corrections(1,:))
                adc_corr(:,1,i,j) = self%adc_corrections(i,j)%p%adc_in
                adc_corr(:,2,i,j) = self%adc_corrections(i,j)%p%adc_out
             end do
          end do
          call write_hdf(chainfile, trim(adjustl(path))//'adc_corr', adc_corr)
          deallocate(adc_corr)
       end if
    end if

    deallocate(amp, amp_tot, R, R_tot)

  end subroutine dumpToHDF_LFI

  subroutine sample_1Hz_spikes(tod, handle, map_sky, procmask, procmask2)
    !   Sample LFI specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    implicit none
    class(comm_LFI_tod),                          intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, k, bin, ierr, nbin
    real(dp)     :: dt, t_tot, t, A, b, mval
    real(dp)     :: t1, t2
    character(len=6) :: scantext
    real(dp), allocatable, dimension(:)     :: nval
    real(sp), allocatable, dimension(:)     :: res
    real(dp), allocatable, dimension(:,:)   :: s_sum
    real(dp), allocatable, dimension(:,:,:) :: s_bin
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '   --> Sampling 1Hz spikes'

    dt    = 1.d0/tod%samprate   ! Sample time
    t_tot = 1.d0                ! Time range in sec
    nbin  = tod%nbin_spike      ! Number of bins 

    allocate(s_bin(0:nbin-1,tod%ndet,tod%nscan), s_sum(0:nbin-1,tod%ndet), nval(0:nbin-1))

    ! Compute template per scan
    s_bin = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       tod%apply_inst_corr = .false. ! Disable 1Hz correction for just this call
       call sd%init_singlehorn(tod, i, map_sky, procmask, procmask2)
       tod%apply_inst_corr = .true.  ! Enable 1Hz correction again

       allocate(res(tod%scans(i)%ntod))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle

          !write(*,*) tod%scanid(i), j, maxval(abs(sd%tod(:,j))), tod%scans(i)%d(j)%gain, maxval(abs(sd%s_sky(:,j))), maxval(abs(sd%s_sl(:,j))), maxval(abs(sd%s_orb(:,j)))
          !res = sd%tod(:,j)/tod%scans(i)%d(j)%gain - (sd%s_sky(:,j) + &
          !     & sd%s_sl(:,j) + sd%s_orb(:,j))
          do k = 1, tod%scans(i)%ntod
             if (sd%tod(k,j) /= sd%tod(k,j)) then
                write(*,*) tod%scanid(i), j, k, sd%tod(k,j), tod%scans(i)%d(j)%gain, sd%s_sky(k,j), sd%s_sl(k,j), sd%s_orb(k,j)
             end if
             !res(k) = 1/tod%scans(i)%d(j)%gain - (sd%s_sky(k,j) + &
             !     & sd%s_sl(k,j) + sd%s_orb(k,j))
             res(k) = sd%tod(k,j)/tod%scans(i)%d(j)%gain - (sd%s_sky(k,j) + &
                  & sd%s_sl(k,j) + sd%s_orb(k,j))
          end do

          nval = 0.d0
          do k = 1, tod%scans(i)%ntod
             if (sd%mask(k,j) == 0.) cycle
             t = modulo(tod%scans(i)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
             bin = min(int(t*nbin),nbin-1)
             s_bin(bin,j,i) = s_bin(bin,j,i)  + res(k)
             nval(bin)      = nval(bin)       + 1.d0
          end do
          if (all(nval > 0)) then
             s_bin(:,j,i) = s_bin(:,j,i) / nval
             s_bin(:,j,i) = s_bin(:,j,i) - mean(s_bin(1:nbin/3,j,i))
          else
             s_bin(:,j,i) = 0.d0
             tod%scans(i)%d(j)%accept = .false.
          end if
       end do

!!$       if (trim(tod%freq) == '070') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_22S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,10,i)
!!$          end do
!!$          close(58)
!!$       end if
!!$
!!$       if (trim(tod%freq) == '044') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_26S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,6,i)
!!$          end do
!!$          close(58)
!!$       end if

       ! Clean up
        call sd%dealloc
       deallocate(res)
    end do

    ! Compute smoothed templates
    s_sum = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       do j = 1, tod%ndet
          s_sum(:,j) = s_sum(:,j) + s_bin(:,j,i)
       end do
    end do
    call mpi_allreduce(mpi_in_place, s_sum,  size(s_sum),  &
         & MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)

    ! Normalize to maximum of unity
    do j = 1, tod%ndet
       !s_sum(:,j) = s_sum(:,j) - median(s_sum(:,j)) 
       mval = maxval(abs(s_sum(:,j))) 
       do k = 0, nbin-1
          s_sum(k,j) = s_sum(k,j) / mval
       end do
       tod%spike_templates(:,j) = s_sum(:,j) 
    end do

    ! Compute amplitudes per scan and detector
    tod%spike_amplitude = 0.
    do i = 1, tod%nscan
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          b = sum(s_sum(:,j)*s_bin(:,j,i)) / tod%scans(i)%d(j)%N_psd%sigma0**2
          A = sum(s_sum(:,j)**2)           / tod%scans(i)%d(j)%N_psd%sigma0**2
          if (A == 0.d0) then
             tod%spike_amplitude(i,j) = 0.d0
          else
             tod%spike_amplitude(i,j) = b/A
             if (trim(tod%operation) == 'sample') &
                  & tod%spike_amplitude(i,j) = tod%spike_amplitude(i,j) + &
                  &                            rand_gauss(handle) / sqrt(A)
          end if
       end do
    end do

    ! Clean up
    deallocate(s_bin, s_sum, nval)

  end subroutine sample_1Hz_spikes

  subroutine construct_corrtemp_LFI(self, scan, pix, psi, s)
    !  Construct an LFI instrument-specific correction template; for now contains 1Hz template only
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_LFI_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    s = 0.
    return

    dt    = 1.d0/self%samprate   ! Sample time
    t_tot = 1.d0                ! Time range in sec
    nbin  = self%nbin_spike      ! Number of bins 

    do j = 1, self%ndet
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t = modulo(self%scans(scan)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
          b = min(int(t*nbin),nbin-1)
          s(k,j) = self%spike_amplitude(scan,j) * self%spike_templates(b,j)
       end do
    end do

  end subroutine construct_corrtemp_LFI


  subroutine preprocess_L1_to_L2(self, map_sky, procmask)
    implicit none
    class(comm_LFI_tod),                          intent(inout) :: self
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask

    integer(i4b) :: i, j, k, m, n, npix, unit, barrier, mpistat(MPI_STATUS_SIZE), ierr
    character(len=4)   :: id
    character(len=6)   :: scantext
    character(len=512) :: filename
    type(hdf_file) :: h5_file
    real(dp), allocatable, dimension(:,:)     :: m_buf
 !   real(sp), allocatable, dimension(:)       :: procmask
    real(sp), allocatable, dimension(:,:)     :: tod

!    npix = 12*self%nside**2

    ! Distribute processing masks
!!$    allocate(m_buf(0:npix-1,self%nmaps), procmask(0:npix-1))
!!$    call self%procmask%bcast_fullsky_map(m_buf); procmask  = m_buf(:,1)
!!$    deallocate(m_buf)
    
    !unit = getlun()
    !open(unit, file=trim(self%L2file), form='unformatted')

    if (self%L2_exist) then
       if (self%myid == 0) write(*,*) "Reading L2 from ", trim(self%L2file)
       call open_hdf_file(self%L2file, h5_file, 'r')
    end if
    
    ! Reduce all scans
    do i = 1, self%nscan
       !call update_status(status, "L1_to_L2")

       ! Generate detector TOD
       n = self%scans(i)%ntod
       allocate(tod(n, self%ndet))
       if (self%L2_exist) then
          call int2string(self%scanid(i), scantext)
          call read_hdf(h5_file, scantext, tod)
       else
          if (self%myid == 0 .and. i == 1) write(*,*) "Converting L1 to L2"
          call self%diode2tod_inst(i, map_sky, procmask, tod)
       end if

       ! Store relevant data
       do j = 1, self%ndet
          if (any(isnan(tod(:,j)))) self%scans(i)%d(j)%accept = .false.
          if (.not. self%scans(i)%d(j)%accept) cycle 
          allocate(self%scans(i)%d(j)%tod(n))
          self%scans(i)%d(j)%tod = tod(:,j)
       end do

!!$       ! Find effective TOD length
!!$       if (self%halfring_split == 0) then
!!$          m = get_closest_fft_magic_number(n)
!!$       else if (self%halfring_split == 1 .or. self%halfring_split == 2) then
!!$          m = get_closest_fft_magic_number(n/2)
!!$       else 
!!$          write(*,*) "Unknown halfring_split value in read_hdf_scan"
!!$          stop
!!$       end if
!!$       if (real(m-n,dp)/real(n,dp) > 0.001d0) then
!!$          write(*,*) 'Warning: More than 0.1% of scan', self%scanid(i), ' removed by FFTW cut'
!!$       end if
!!$
!!$       if (m /= n) write(*,*) 'ERROR', self%scanid(i), m, n
       
       ! Free up old arrays
       do j = 1, self%ndet
!!$          if (.not. self%scans(i)%d(j)%accept) cycle 
!!$
!!$          allocate(self%scans(i)%d(j)%tod(m))
!!$          if (self%halfring_split == 2) then
!!$             self%scans(i)%d(j)%tod = tod(m+1:2*m,j)
!!$          else
!!$             self%scans(i)%d(j)%tod = tod(1:m,j)
!!$          end if
          if (allocated(self%scans(i)%d(j)%ztod))   deallocate(self%scans(i)%d(j)%ztod)
          if (allocated(self%scans(i)%d(j)%diode))  deallocate(self%scans(i)%d(j)%diode)
          if (allocated(self%scans(i)%d(j)%zdiode)) then
             do k = self%ndiode, 1, -1
                deallocate(self%scans(i)%d(j)%zdiode(k)%p) 
             end do
             deallocate(self%scans(i)%d(j)%zdiode)
          end if
        end do
        call huff_deallocate(self%scans(i)%todkey)
        deallocate(tod)
     end do
    if (self%L2_exist) call close_hdf_file(h5_file)

    if (.not. self%L2_exist) then
       ! Write preprocessed data to file in round-robin manner
       barrier = self%nscan
       if (self%myid == 0) then
          write(*,*) "Writing L2 to ", trim(self%L2file)
          call open_hdf_file(self%L2file, h5_file, 'w')
          call write_hdf(h5_file, "freq", trim(self%freq))
       end if
       if (self%myid > 0) then
          call mpi_recv(barrier, 1, MPI_INTEGER, self%myid-1, 98, self%comm, mpistat, ierr)
          call open_hdf_file(self%L2file, h5_file, 'b')
       end if
       do i = 1, self%nscan
          call int2string(self%scanid(i), scantext)
          n = self%scans(i)%ntod
          allocate(tod(n,self%ndet))
          do j = 1, self%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle 
             tod(:,j) = self%scans(i)%d(j)%tod
          end do
          call write_hdf(h5_file, scantext, tod)
          deallocate(tod)
       end do
       call close_hdf_file(h5_file)
       if (self%myid < self%numprocs-1) then
          call mpi_send(barrier, 1, MPI_INTEGER, self%myid+1, 98, self%comm, ierr)      
       end if
       call mpi_barrier(self%comm, ierr)
    end if

     !deallocate(procmask)

  end subroutine preprocess_L1_to_L2

end module comm_tod_LFI_mod

