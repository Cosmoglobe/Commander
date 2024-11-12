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
module comm_tod_HFI_mod
  !   Module which contains all the HFI time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates
  !       all data needed for TOD processing
  !   process_HFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
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
  implicit none

  private
  public comm_HFI_tod

  type, extends(comm_tod) :: comm_HFI_tod
     integer(i4b), allocatable, dimension(:,:)   :: mod_phase
   contains
     procedure     :: process_tod             => process_HFI_tod
     procedure     :: read_tod_inst           => read_tod_inst_HFI
     procedure     :: read_scan_inst          => read_scan_inst_HFI
     procedure     :: initHDF_inst            => initHDF_HFI
     procedure     :: dumpToHDF_inst          => dumpToHDF_HFI
     procedure     :: construct_corrtemp_inst => construct_corrtemp_hfi
     procedure     :: apply_nonlin_corr_inst  => apply_nonlin_corr_hfi
  end type comm_HFI_tod

  interface comm_HFI_tod
     procedure constructor
  end interface comm_HFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info, tod_type)
    !
    ! Constructor function that gathers all the instrument parameters in a pointer
    ! and constructs the objects
    !
    ! Arguments:
    ! ----------
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
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_HFI_tod),     pointer    :: constructor

    integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ierr
    logical(lgt) :: pol_beam

    ! Allocate object
    allocate(constructor)

    ! Set up noise PSD type and priors
    constructor%freq            = cpar%ds_label(id_abs)
    constructor%n_xi            = 3
    constructor%noise_psd_model = 'oof'
    allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
    allocate(constructor%xi_n_P_rms(constructor%n_xi))

    ! just so that it actually runs
    constructor%xi_n_P_uni(2,:) = [0.10d0, 1.0d0]  ! fknee
    constructor%xi_n_P_uni(3,:) = [-2.5d0, -0.4d0]   ! alpha
    constructor%xi_n_nu_fit     = [0.d0, 1.225d0] ! I took it from freq=30 for LFI, so not true

    constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0] ! [sigma0, fknee, alpha]; sigma0 is not used

    constructor%n_cray_temps    = 3

    ! Initialize common parameters
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)

    ! Initialize instrument-specific parameters
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn           = 1
    constructor%ndiode          = 1
    constructor%compressed_tod  = .true.
    constructor%correct_sl      = .false.
    constructor%orb_4pi_beam    = .false.
    constructor%symm_flags      = .false.
    constructor%chisq_threshold = 1500.d0 !20.d0 ! 9.d0
    constructor%nmaps           = info%nmaps
    constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), "," )
    constructor%ntime           = 1
    constructor%partner         = -1
    constructor%n_cray_temps    = 3

    nside_beam                  = 512
    nmaps_beam                  = 3
    pol_beam                    = .true.
    constructor%nside_beam      = nside_beam
    
    ! Get detector labels
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
 
    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Initialize bandpass mean and proposal matrix
    !call constructor%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    call constructor%precompute_lookups()

    ! Load the instrument file
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Allocate sidelobe convolution data structures
    !allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
    
    allocate(constructor%orb_dp)
    constructor%orb_dp => comm_orbdipole(constructor%mbeam)
    
    ! Initialize all baseline corrections to zero
    do i = 1, constructor%nscan
       constructor%scans(i)%d%baseline = 0.d0
    end do

    ! Allocate modulation phase
    allocate(constructor%mod_phase(constructor%ndet,constructor%nscan))
    constructor%mod_phase = 1.d0
    
  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_HFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    !
    ! Routine that processes the HFI time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_HFI_tod class
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
    class(comm_HFI_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
    type(map_ptr),       dimension(1:,1:),    intent(inout), optional :: map_gain       ! (ndet,1)

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
    logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, sigma0
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    select_data           = self%first_call        ! only perform data selection the first time
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration

    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    self%output_n_maps = 3
    if (self%output_aux_maps > 0 .or. .true.) then
       if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
    end if

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Distribute maps
    allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
    call distribute_sky_maps(self, map_in, 1.e0, map_sky) ! already in K_cmb

    ! Distribute processing masks
    allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
    call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
    call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    deallocate(m_buf)

    call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
    !call self%procmask%writeFITS("mask.fits")
    
!    call mpi_finalize(ierr)
!    stop

    ! Precompute far sidelobe Conviqt structures
    !if (self%correct_sl) then
    !   if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
    !   do i = 1, self%ndet
          !TODO: figure out why this is rotated
    !      call map_in(i,1)%p%YtW()  ! Compute sky a_lms
    !      self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
    !           & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
    !           & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
    !   end do
    !end if
!    write(*,*) 'qqq', self%myid
!    if (.true. .or. self%myid == 78) write(*,*) 'a', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
!!$    call mpi_finalize(ierr)
!!$    stop

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Fit low-level non-linearity parameters
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, skip_nonlin=.true.)

       ! Estimate modulation baselines; separate for odd and even samples
       if (self%first_call) then
          call sample_hfi_baselines(sd, self, i, handle, subtract_s_tot=.false.)
       else
          call sample_hfi_baselines(sd, self, i, handle)
       end if

       ! Fix modulation phase
       if (self%first_call) call set_modulation_phase(sd, self, i)

       ! Estimate ADC corrections
       !    Not implemented yet
       
       ! Clean up
       call dealloc_scandata(sd)
    end do
    
    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    call sample_calibration(self, 'abscal', handle, map_sky, procmask2, procmask2)
    !call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
    !call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)

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

    ! Fit higher-level corrections
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
!          if (.true. .or. self%myid == 78) write(*,*) 'b', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
          call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
          call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_sky_prop=.true.)
       else
          call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
       end if

!!$       ! Calling Simulation Routine
!!$       if (self%enable_tod_simulations) then
!!$          call simulate_tod(self, i, sd%s_tot, sd%n_corr, handle)
!!$          call sd%dealloc
!!$          cycle
!!$       end if

       ! ************************************************************************************
       !       Perform low-level TOD processing; raw modulated TOD -> clean demodulated TOD
       ! ************************************************************************************

       if (self%first_call) then
          ! Search for jumps
          !    Not implemented yet
       end if

       ! Apply jump corrections
       !    Not implemented yet       

       ! Fit and subtract dark bolometer signal; should this come before or after cosmic rays..?
       !    Not implemented yet       

       ! remove cosmic rays
!!$       do j=1, sd%ndet
!!$        call self%cray(j)%p%build_cray_templates()
!!$        call self%cray(j)%p%fit_cray_amplitudes()
!!$       end do

       ! Fit and subtract 4K lines
       !    Not implemented yet       

       ! Fit and deconvolve bolometer transfer function
       !    Not implemented yet       



       ! ************************************************************************************
       !       Perform high-level TOD processing; clean demodulated TOD -> maps
       ! ************************************************************************************

       ! Sample correlated noise
       !call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)
       sd%n_corr = 0.

       ! Compute noise spectrum parameters
       !call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), tod=sd%tod(:, j))
       end do

       ! Select data
       !if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute calibrated TOD for mapmaking
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)
       

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
       call dealloc_scandata(sd)
       deallocate(d_calib)

    end do


    ! calculate gain as a median of all scans 
    !do j = 1, self%ndet
    !   allocate(temp(self%nscan))

    !   do i = 1, self%nscan
    !      temp(i) = self%scans(i)%d(j)%gain
    !   end do

    !   gain_median = median(temp)

    !   do i = 1, self%nscan
    !      self%scans(i)%d(j)%gain = gain_median
    !   end do

    !   deallocate(temp)
    !end do


    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    call synchronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       if (self%nmaps > 1) then
         call finalize_binned_map(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2)
       else
         call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0, chisq_S=chisq_s, mask=procmask2)
       end if
    else
       if(self%nmaps > 1) then
         call finalize_binned_map(self, binmap, rms_out, 1.d0)
       else 
         call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0)
       end if
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

  end subroutine process_HFI_tod


  subroutine sample_hfi_baselines(self, tod, scan, handle, subtract_s_tot)
    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !                     tod%scans(scan)%d(:)%gain (temporary solution)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
    type(planck_rng),                     intent(inout) :: handle
    logical(lgt),                         intent(in), optional :: subtract_s_tot
    
    real(dp) :: eta, A1, A2, x, b1, b2, sgn,gal_mean
    integer(i4b) :: i, j, n
    logical(lgt) :: sub_s

    sub_s = .true.; if (present(subtract_s_tot)) sub_s = subtract_s_tot
    
    
    ! tod%scans(scan)%d(i)%gain - the gain constant over a scan [real number]
    ! sd = self --- self%s_tot - sky signal model
    ! self%s_tot(self%ntod, self%ndet) - how s_tot structured

    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       sgn = tod%mod_phase(i,scan)
       
       ! Odd samples
       A1 = 0.d0; b1 = 0; gal_mean = 0.d0; n = 0
       do j = 1, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A1 = A1 + 1.d0
          b1 = b1 + self%tod(j,i)
          if (sub_s) b1 = b1 - sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
       end do
       A1 = A1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b1 = b1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       tod%scans(scan)%d(i)%baseline  = b1/A1 + rand_gauss(handle)/sqrt(A1)

       ! Even samples
       A2 = 0.d0; b2 = 0.d0
       do j = 2, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A2 = A2 + 1.d0
          b2 = b2 + self%tod(j,i)
          if (sub_s) b1 = b1 + sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
       end do
       A2 = A2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b2 = b2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       tod%scans(scan)%d(i)%baseline2  = b2/A2 + rand_gauss(handle)/sqrt(A2)

       !write(*,*) "baseline=", tod%scans(scan)%d(i)%baseline, tod%scans(scan)%d(i)%baseline2, tod%mod_phase(i,scan)

    end do

  end subroutine sample_hfi_baselines

  subroutine set_modulation_phase(self, tod, scan)
    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !                     tod%scans(scan)%d(:)%gain (temporary solution)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
    
    real(dp) :: mu, n
    integer(i4b) :: i, j
    
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle       
       
       mu = 0.d0; n = 0.d0
       do j = 1, self%ntod, 2
          if (self%pix(j,i,1) > 0.48*tod%info%npix .and. self%pix(j,i,1) < 0.52*tod%info%npix) then
             if (mod(j,2) == 1) then
                mu = mu + self%tod(j,1)-tod%scans(scan)%d(i)%baseline
             else
                mu = mu - self%tod(j,1)-tod%scans(scan)%d(i)%baseline
             end if
             n  = n  + 1.d0
          end if
       end do

       if (n == 0) then
          write(*,*) "Scan disabled in set_modulation_phase"
          tod%scans(scan)%d(i)%accept = .false.
       else
          mu = mu/n

          ! saving the phase to the tod object
          if (mu < 0.d0) then
             tod%mod_phase(i,scan) = -1
          end if
       end if
    end do

  end subroutine set_modulation_phase
  

  subroutine demodulate_tod(self, tod, scan)
    ! 
    ! Demodulate HFI TOD
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     Scan ID
    !
    ! Returns
    ! ----------
    !   None, but updates self%tod
    !
    implicit none
    class(comm_scandata),                         intent(inout) :: self
    class(comm_hfi_tod),                          intent(in)    :: tod
    integer(i4b),                                 intent(in)    :: scan

    integer(i4b) :: i, j
    real(sp)     :: sgn
    logical :: exists

    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle       
       sgn = tod%mod_phase(i,scan)
       
       ! Subtract baselines and flip sign of every other sample
       do j = 1, self%ntod
           if (mod(j,2) == 1) then
               self%tod(j,i) =  sgn*(self%tod(j,i) - tod%scans(scan)%d(i)%baseline)
           else
               self%tod(j,i) = -sgn*(self%tod(j,i) - tod%scans(scan)%d(i)%baseline2)
           end if
       end do

    end do

!!$    open(58,file='tod.dat')
!!$    do j = 1, self%ntod
!!$       write(58,*) j, self%tod(j,1), self%mask(j,1)
!!$    end do
!!$    close(58)
    
  end subroutine demodulate_tod

  
  subroutine read_tod_inst_HFI(self, file)
    ! 
    ! Reads HFI-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_HFI_tod),                 intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_HFI
  
  subroutine read_scan_inst_HFI(self, file, slabel, detlabels, scan)
    ! 
    ! Reads HFI-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle
    ! slabel:   string
    !           Scan label, e.g., "000001/"
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! scan:     derived class (comm_scan)
    !           Scan object
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    implicit none
    class(comm_HFI_tod),                 intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_HFI

  subroutine initHDF_HFI(self, chainfile, path)
    ! 
    ! Initializes HFI-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
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
    class(comm_HFI_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_HFI
  
  subroutine dumpToHDF_HFI(self, chainfile, path)
    ! 
    ! Writes HFI-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
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
    class(comm_HFI_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_HFI

  subroutine construct_corrtemp_hfi(self, scan, pix, psi, s)
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
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    s = 0.

  end subroutine construct_corrtemp_hfi

  subroutine apply_nonlin_corr_hfi(self, scan, sd)
    !  Construct and apply HFI instrument-specific non-linear corrections
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
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd

    ! Apply ADC corrections to raw self%tod
    !    Not implemented yet

    ! Demodulate TOD
    call demodulate_tod(sd, self, scan)
    
  end subroutine apply_nonlin_corr_hfi

end module comm_tod_HFI_mod
