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
module comm_tod_chipass_mod
   !   Module which contains all the chipass time ordered data processing and routines
   !   for a given frequency band
   !
   !   Main Methods
   !   ------------
   !   constructor(cpar, id_abs, info, tod_type)
   !       Initialization routine that reads in, allocates and associates 
   !       all data needed for TOD processing
   !   process_chipass_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
   !       Routine which processes the time ordered data
   !
  use comm_tod_mapmaking_mod
  use comm_tod_driver_mod
  use comm_tod_pixhist_mod
  implicit none

   private
   public comm_chipass_tod

   type, extends(comm_tod) :: comm_chipass_tod
      integer(i4b)  :: tsys_order
      real(dp)      :: tsys_eta0
      real(dp), allocatable, dimension(:)      :: el_min, el_max
      real(dp), allocatable, dimension(:,:)    :: tsys_fit ! ndet, tsys_order+1
      type(spline_type), allocatable, dimension(:) :: tsys_spline
      class(comm_dynmask), pointer :: dynmask
    contains
      procedure     :: process_tod          => process_chipass_tod
      procedure     :: read_scan_inst          => read_scan_inst_chipass
      procedure     :: construct_corrtemp_inst => construct_corrtemp_chipass
      procedure     :: initHDF_inst            => initHDF_chipass
      procedure     :: dumpToHDF_inst          => dumpToHDF_chipass
   end type comm_chipass_tod

   interface comm_chipass_tod
      procedure constructor_chipass
   end interface comm_chipass_tod

contains
   !**************************************************
   !             Constructor
   !**************************************************
   function constructor_chipass(cpar, id, id_abs, info, tod_type) result(c)
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
      type(comm_params),       intent(in) :: cpar          !comm_param structure, list of all the input parameters
      integer(i4b),            intent(in) :: id, id_abs        !index of the current band within the parameters 
      class(comm_mapinfo),     target     :: info
      character(len=128),      intent(in) :: tod_type      !
      class(comm_chipass_tod), pointer    :: c

      integer(i4b) :: i, j, nside_beam, lmax_beam, nmaps_beam, ierr
      logical(lgt) :: pol_beam
      real(dp)     :: samprate

      call timer%start(TOD_INIT, id_abs)

      ! Allocate object
      allocate(c)

      ! Set up noise PSD type and priors
      samprate          = 0.2d0 ! Hz
      c%freq            = cpar%ds_label(id_abs)
      c%n_xi            = 3
      c%noise_psd_model = 'oof'
      allocate(c%xi_n_nu_fit(c%n_xi,2))
      allocate(c%xi_n_P_uni(c%n_xi,2))
      allocate(c%xi_n_P_rms(c%n_xi))
      c%xi_n_P_rms      = [1.d0, 0.01d0, 0.2d0] 
      ! [sigma0, fknee, alpha]; sigma0 is not used
      c%xi_n_nu_fit(1,:) = [0.001d0, 0.99*samprate/2]        ! sigma0
      c%xi_n_nu_fit(2,:) = [0.001d0, 0.99*samprate/2]         ! fknee
      c%xi_n_nu_fit(3,:) = [0.001d0, 0.99*samprate/2]         ! alpha
      c%xi_n_P_uni(1,:)  = [0.001d0, 10.d0]       ! sigma0
      !c%xi_n_P_uni(2,:)  = [0.00001d0, 0.1d0]  ! fknee
      !c%xi_n_P_uni(3,:)  = [-3.0d0,   -0.4d0]  ! alpha
      c%xi_n_P_uni(2,:)  = [0.002d0, 0.0021d0]  ! fknee
      !
      c%xi_n_P_uni(3,:)  = [-1.801d0, -1.799d0]  ! alpha

      ! Initialize common parameters
      call c%tod_constructor(cpar, id, id_abs, info, tod_type)
      
      ! Initialize instrument-specific parameters
      !read(c%freq(1:2),*) c%zodiband
      c%enable_fft_magic = .false.
      c%samprate_lowres = 0.2d0  ! Lowres samprate in Hz
      c%nhorn           = 1
      c%ndiode          = 1
      c%baseline_order  = 1
      c%apply_inst_corr = .true.
      c%compressed_tod  = .false.
      c%correct_sl      = .false.
      c%correct_orb     = .false.
      c%orb_4pi_beam    = .false.
      c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
      c%symm_flags      = .false.
      c%read_elev       = .true.
      c%read_az         = .false.
      c%per_slew_baseline = .true.
      c%max_nslew       = 10
      ! c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      c%chisq_threshold = 30
      c%nmaps           = info%nmaps
      if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
         c%ndet         = count_detectors(cpar%ds_tod_dets(id_abs)) !, cpar%datadir)
      else
         c%ndet         = num_tokens(cpar%ds_tod_dets(id_abs), ",")
      end if
            
      nside_beam                  = 1024
      nmaps_beam                  = 1
      pol_beam                    = .false.
      c%nside_beam      = nside_beam
      c%nside_pixhist   = 64
      
      ! allocate CHIPASS-specific instrument file data
      c%tsys_order      = 7
      c%tsys_eta0       = 58.0d0
      allocate(c%tsys_fit(c%ndet, 0:c%tsys_order))
      c%tsys_fit = 0.d0

      ! Get detector labels
      if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
         call get_detectors(cpar%ds_tod_dets(id_abs), c%label)
      else
         call get_tokens(trim(adjustl(cpar%ds_tod_dets(id_abs))), ",", c%label)
      end if
      
      ! Read the actual TOD
      call c%read_tod(c%label)
      
      ! Initialize bandpass mean and proposal matrix
      !call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
      c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false.)
      call c%precompute_lookups()

      ! Load the instrument file
      call c%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

   ! commenting this out for now
      !do i=1, c%ndet
      !  call init_noise_model(c, i)
      !end do

      ! Allocate sidelobe convolution data structures
      !allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)

      !c%orb_dp => comm_orbdipole(comm=c%comm)

      ! Initialize dynamic mask
      c%dynmask => comm_dynmask(c, cpar)
      c%dynmask%output_scan             = 500
      c%dynmask%output_det              = 1
      c%dynmask%threshold_singlesamp    = 10  ! Exclude 10 sigma outliers
      c%dynmask%window_excessRMS(1)     = 7   ! Search for windows of 7 samples
      c%dynmask%threshold_excessRMS(1)  = 10  ! Exclude 10 sigma outliers; planets?

      c%dynmask%apply_pixhist           = .false. !.true.
      !c%dynmask%remove_isolated_samples = .false.
      !c%dynmask%threshold_longchunks    = -0.3
      !c%dynmask%window_longchunks       = -2000
      !c%dynmask%threshold_cr            = -6.   ! sigma
      !c%dynmask%width_cr_mask           = -5
      !allocate(c%dynmask%templ_cr(-5:5))
      !c%dynmask%templ_cr    = 1.
      !c%dynmask%templ_cr(0) = 10.
      !c%dynmask%templ_cr(1) = 5.
      !c%dynmask%templ_cr(2) = 2.5
      !c%dynmask%templ_cr(3) = 1.25
      !c%dynmask%templ_cr(4) = 1.
      !c%dynmask%templ_cr    = c%dynmask%templ_cr - mean(real(c%dynmask%templ_cr,dp))

      ! Reject obviously bad scans
      do i = 1, c%nscan
         do j = 1, c%ndet
            if (c%scans(i)%d(j)%N_psd%sigma0 == 0.d0) c%scans(i)%d(j)%accept = .false.
         end do 
      end do      

      ! Compute elevation range
      allocate(c%el_min(c%ndet), c%el_max(c%ndet))
      c%el_min = 90.d0
      c%el_max =  0.d0
      do j = 1, c%ndet
         do i = 1, c%nscan
            if (c%scans(i)%d(j)%accept) then
               c%el_min = min(c%el_min, minval(c%scans(i)%d(j)%elev,c%scans(i)%d(j)%elev>0.d0))
               c%el_max = max(c%el_max, maxval(c%scans(i)%d(j)%elev,c%scans(i)%d(j)%elev>0.d0))
            end if
         end do 
      end do
      call mpi_allreduce(mpi_in_place, c%el_min, size(c%el_min), &
           & MPI_DOUBLE_PRECISION, MPI_MIN,  c%comm, ierr)
      call mpi_allreduce(mpi_in_place, c%el_max, size(c%el_max), &
           & MPI_DOUBLE_PRECISION, MPI_MAX,  c%comm, ierr)
      if (c%myid == 0) write(*,fmt='(a,2f8.3)') '  Elevation range, det 1 = ', c%el_min(1), c%el_max(1)

      
      call timer%stop(TOD_INIT, id_abs)
    end function constructor_chipass

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_chipass_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      ! 
      ! Routine that processes the chipass Calibrated Individual Observations. 
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
      ! Writes maps to disc in fits format
      ! 
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_chipass_tod class
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
      class(comm_chipass_tod),                  intent(inout) :: self
      character(len=*),                         intent(in)    :: chaindir
      integer(i4b),                             intent(in)    :: chain, iter
      type(planck_rng),                         intent(inout) :: handle
      type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
      real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map),                          intent(inout) :: map_out      ! Combined output map
      class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr),       dimension(1:),       intent(inout), optional :: map_gain       ! (ndet,1)
      real(dp)            :: t1, t2
      integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps, tod_start_idx, n_tod_tot, n_comps_to_fit, oper_default
      logical(lgt)        :: select_data, sample_gain, output_scanlist, sample_ncorr, sample_baseline, sample_tsys
      type(comm_binmap)   :: binmap, binmap2
      type(comm_scandata) :: sd
      character(len=4)    :: ctext, myid_text
      !character(len=2)    :: zodi_param_text
      !character(len=1)    :: up_down_text
      character(len=6)    :: samptext, scantext
      character(len=512)  :: prefix, postfix!, prefix4D, prefix_atlas, postfix_atlas
      character(len=512), allocatable, dimension(:) :: slist
      !real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi
      real(sp), allocatable, dimension(:,:,:)   :: d_calib, d_calib2
      !real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
      !real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf
      !real(dp), allocatable, dimension(:, :)    :: A_T_A, A_T_A_reduced
      !real(dp), allocatable, dimension(:)       :: AY, AY_reduced, X
      !real(dp), allocatable, dimension(:, :, :) :: s_therm_tot, s_scat_tot ! (n_tod_tot, ncomps, ndet)
      !real(dp), allocatable, dimension(:, :)    :: res_tot ! (n_tod_tot, ndet)
      !real(dp), allocatable, dimension(:)       :: mask_tot ! (n_tod_tot)
      type(hdf_file) :: tod_file

      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)
      
      call timer%start(TOD_ALLOC, self%band)
      
      ! Toggle optional operations
      !sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
      !output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
      !use_k98_samp_groups   = .false.                          ! fits one overall albedo and episolon for the dust bands, and one for ring + feature
      !if (self%myid == 0) write(*, *) '[comm_tod_chipass_mod.f90] size(delta,3) > 1:', size(delta,3) > 1, size(delta,1), size(delta,2), size(delta,3)
      !sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      !sample_abs_bandpass   = .false.                         ! don't sample absolute bandpasses

      if (.false.) then
         ! Debug
         select_data     = iter == 2
         sample_baseline = iter > 2
         sample_tsys     = iter > 2
         sample_gain     = iter > 2                         ! Gain sampling
         sample_ncorr    = iter > 2
      else if (trim(self%init_from_HDF) == 'none') then ! OBS FIXME bug when BAND_TOD_INI_:FROM_HDF=default and INIT_CHAIN=none
         select_data     = .false. !iter == 10                   ! in param file. Takes you to 'else' below
         sample_baseline = iter > 1
         sample_tsys     = iter > 1
         sample_gain     = iter > 2                        ! Gain sampling
         sample_ncorr    = iter > 3
      else
         select_data     = self%first_call  !iter == 1                    
         sample_baseline = iter > 1
         sample_tsys     = iter > 1
         sample_gain     = iter > 1                         ! Gain sampling
         sample_ncorr    = iter > 0
      end if
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
       
      ! Initialize local variables
      ndelta          = size(delta,3)
      self%n_bp_prop  = 0 !ndelta-1
      nside           = map_out%info%nside
      nmaps           = map_out%info%nmaps
      npix            = 12*nside**2
      self%output_n_maps = 8
      !if (output_zodi_comps) self%output_n_maps = self%output_n_maps + zodi_model%n_comps

      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'

      ! Define useful sd operation codes
      oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
           & SD_SKY,SD_INST,SD_NCORR])

      ! Initialize index-based sky map and mask
      ! NOTE: CHIPASS TOD given in Jy beam^-1 but parameter file assumes MJy/sr
      ! for a 14.3 arcmin beam the conversion factor is roughly 0.057793
      call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
      call timer%stop(TOD_ALLOC, self%band)
      call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
      call update_status(status, "tod_init")


      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      !if (self%myid == 0) then
      !   write(*,*) 'gain init scan 1 det 1:', self%scans(1)%d(1)%gain
      !   write(*,*) 'baseline init scan 1 det 1:', self%scans(1)%d(1)%baseline
      !   if (self%per_slew_baseline) write(*,*) 'baseline_slew init scan 1 det 1:', self%scans(1)%d(1)%baseline_slew
      !   write(*,*) 'tsys_fit init det 1:', self%tsys_fit(1, :)
      !end if

      ! Sample gain components in separate TOD loops; marginal with respect to n_corr
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_gain:', sample_gain
      if (sample_gain) then
         ! 'abscal': the global constant gain factor
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sampling calibration'
         call sample_calibration(self, 'abscal', oper_default, handle)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sampling calibration done'
         ! 'relcal': the gain factor that is constant in time but varying between detectors
         call sample_calibration(self, 'relcal', oper_default, handle)         
         ! 'deltaG': the time-variable and detector-variable gain
         !call sample_calibration(self, 'deltaG', oper_default, handle)
      end if

      ! TODO: sample tsys_eta0

      ! sample Tsys; marginal over baseline -- must be directly followd by baseline sampling
      if (sample_tsys) then
         if (self%myid == 0) then
            write(*,*) '|    --> Sampling Tsys'
         end if
         call update_status(status, "Tsys")
         call sample_chipass_Tsys(self, oper_default, handle)
      end if

      ! sample baseline
      if (sample_baseline) then
         if (self%myid == 0) then
            write(*,*) '|    --> Sampling baseline'
         end if
         call update_status(status, "baseline")
         do i = 1, self%nscan
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, spur_level=0)
            call timer%start(TOD_BASELINE, self%band)
            if (self%per_slew_baseline) then
               call sample_chipass_baseline_per_slew(self, i, sd%tod, sd%s_tot(:,:,0,1), sd%mask)
            else
               call sample_chipass_baseline(self, i, sd%tod, sd%s_tot(:,:,0,1), sd%mask)
            end if
            call timer%stop(TOD_BASELINE, self%band)
            call dealloc_scan_data(sd)
         end do
         call timer%start(TOD_WAIT, self%band)
         call mpi_barrier(self%comm, ierr)
         call timer%stop(TOD_WAIT, self%band)
      end if

      !if (self%myid == 0) then
      !   write(*,*) 'gain new scan 1 det 1:', self%scans(1)%d(1)%gain
      !   write(*,*) 'baseline new scan 1 det 1:', self%scans(1)%d(1)%baseline
      !   if (self%per_slew_baseline) write(*,*) 'baseline_slew new scan 1 det 1:', self%scans(1)%d(1)%baseline_slew
      !   write(*,*) 'tsys_fit new det 1:', self%tsys_fit(1, :)
      !end if

      ! Create pixel histograms
      if (self%first_call) call compute_tod_pixhist(self, spur_level=100)
      call update_status(status, "tod_hist"//ctext)
      
       ! Create dynamic mask
      if (select_data) then
         if (self%myid == 0) write(*,*) '   --> Creating dynamic mask'
         if (self%myid == 0) write(*,*) 'scan 1 det 1 N_psd%sigma0 1:', self%scans(1)%d(1)%N_psd%sigma0
         do i = 1, self%nscan
            ! Skip scan if no accepted data
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd) 
            do j = 1, sd%ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               call self%dynmask%create(sd, j)
            end do
            !if ((self%myid == 0) .and. (i == 1)) then
            !   write(*,*) '    [chipass] tod:', sd%tod(1:10,1)
            !   write(*,*) '    [chipass] s_tot adj:', self%scans(i)%d(1)%gain*sd%s_tot(1:10,1,0,1)
            !   write(*,*) '    [chipass] sig res:', &
            !           & (sd%tod(1:10,1) - self%scans(i)%d(1)%gain*sd%s_tot(1:10,1,0,1))/self%scans(i)%d(1)%N_psd%sigma0
            !end if
            call dealloc_scan_data(sd)
         end do
         ! Synchronize and output flagging statistics in first iteration
         call self%dynmask%report
      end if

      
      ! Prepare intermediate data structures
      call binmap%init(self, .true., .false.)
      call binmap2%init(self, .true., .false.)

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
         call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)

         ! Sample correlated noise
         if (sample_ncorr) then
            call sample_n_corr(self, sd, handle)
            call sample_noise_psd(self, sd, handle, chaindir)
         else
            sd%n_corr = 0.d0
            call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true.)
         end if

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_tod_chisq(sd, j)
         end do

         ! Select data
         !if (select_data) call remove_bad_data(self, i, sd%flag)
         if (iter==3) call remove_bad_data(self, i, sd%flag)
         
         ! Compute binned map
         allocate(d_calib(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib = 0.d0
         call compute_calibrated_data(self, i, sd, d_calib)
         allocate(d_calib2(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib2 = 0.d0 ! baseline, Tsys
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            !do k = 1, self%scans(i)%ntod
            !   do l = 0, self%tsys_order
            !      d_calib2(2, :, j) = d_calib2(2, :, j) + self%tsys_fit(j,l) * (self%scans(i)%d(j)%elev(k) - self%tsys_eta0)**l
            !   end do
            !end do
            if (allocated(self%tsys_spline)) then
               do k = 1, self%scans(i)%ntod
                  d_calib2(2, k, j) = d_calib2(2, k, j) + splint(self%tsys_spline(j), self%scans(i)%d(j)%elev(k))
               end do
               !d_calib2(2, :, j) = d_calib2(2, :, j) + splint(self%tsys_spline(j), self%scans(i)%d(j)%elev)
               !call splint_simple_multi(self%tsys_spline(j), self%scans(i)%d(j)%elev, d_calib2(2, :, j))
            else
               do l = 0, self%tsys_order
                  d_calib2(2, :, j) = d_calib2(2, :, j) + self%tsys_fit(j,l) * (self%scans(i)%d(j)%elev - self%tsys_eta0)**l
               end do
            end if
            d_calib2(1, :, j) = sd%s_inst(:, j) - d_calib2(2, :, j)
         end do

         ! For debugging: write TOD to hdf
         if (.true.) then
            if (mod(self%scanid(i), 500) == 0) then
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod/self%scans(i)%d(1)%gain)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call write_hdf(tod_file, '/flag', sd%flag)
               call write_hdf(tod_file, '/calib', d_calib(1, :, :))
               call write_hdf(tod_file, '/s_sky', sd%s_sky)
               call write_hdf(tod_file, '/n_corr', sd%n_corr/self%scans(i)%d(1)%gain)
               call write_hdf(tod_file, '/s_inst', sd%s_inst/self%scans(i)%d(1)%gain)
               !call write_hdf(tod_file, '/s_sl', sd%s_sl)
               !call write_hdf(tod_file, '/s_orb', sd%s_orb)
               call write_hdf(tod_file, '/res', d_calib(2, :, :))
               call write_hdf(tod_file, '/mask', sd%mask)
               call write_hdf(tod_file, '/sigma0', self%scans(i)%d(1)%N_psd%sigma0)
               call close_hdf_file(tod_file)
            end if
         end if

         ! Bin TOD
         call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)
         call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib2, binmap2)

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
         call dealloc_scan_data(sd)
         deallocate(d_calib)
         deallocate(d_calib2)
      end do

      if (self%myid == 0) write(*,*) 'scan 1 det 1 N_psd%sigma0 2:', self%scans(1)%d(1)%N_psd%sigma0

      if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

      ! Output latest scan list with new timing information
      if (output_scanlist) call self%output_scan_list(slist)

      ! Solve for maps
      call synchronize_binmap(binmap2, self)
      call finalize_binned_map_unpol(self, binmap2, rms_out, 1.d0)
      call synchronize_binmap(binmap, self)
      call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0)
      map_out%map = binmap%outmaps(1)%p%map

      ! Inpaint missing pixels
      call map_out%inpaint_misspix(rms_out, map_in(1,1)%p, 30.d0, handle)
      
      ! Output maps to disk
      call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
      call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
      if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
      if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
      ! added
      !if (self%myid == 0) write(*,*) '   [comm_tod_chipass_mod] output_n_maps:', self%output_n_maps
      if (self%output_n_maps > 7) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'inst_corr'//trim(postfix))
      call binmap2%outmaps(1)%p%writeFITS(trim(prefix)//'baseline'//trim(postfix))
      call binmap2%outmaps(2)%p%writeFITS(trim(prefix)//'Tsys_eta'//trim(postfix))

      ! Clean up
      call binmap%dealloc()
      call binmap2%dealloc()
      if (allocated(slist)) deallocate(slist)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

      call update_status(status, "tod_end"//ctext)
      
      call timer%stop(TOD_TOT, self%band)
   end subroutine process_chipass_tod


   subroutine sample_chipass_baseline(tod, scan, raw, s_tot, mask)
    !
    !   Sample CHIPASS baseline
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   scan:     local scan ID
    !   raw:      raw tod in du
    !   s_tot:    total signal model in mK
    !   mask:     list of accepted samples
    !
    implicit none
    class(comm_chipass_tod),                intent(inout) :: tod
    integer(i4b),                           intent(in)    :: scan
    real(sp), dimension(1:,1:),             intent(in)    :: raw, s_tot, mask

    integer(i4b) :: i, j, k, n
    real(dp)     :: t, dt
    real(dp), allocatable, dimension(:) :: x, y

    allocate(x(tod%scans(scan)%ntod), y(tod%scans(scan)%ntod))
    dt = 1.d0 / tod%scans(scan)%ntod

    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       t = 0.d0
       n = 0
       do k = 1, tod%scans(scan)%ntod
          t      = t + dt
          if (mask(k,j) > 0.5) then
             n    = n + 1
             x(n) = t
             y(n) = raw(k,j) - tod%scans(scan)%d(j)%gain * s_tot(k,j)
             if (allocated(tod%tsys_spline)) then
                y(n) = y(n) - splint(tod%tsys_spline(j),tod%scans(scan)%d(j)%elev(k))
             else
                do i = 0, tod%tsys_order
                   y(n) = y(n) - tod%tsys_fit(j, i) * (tod%scans(scan)%d(j)%elev(k) - tod%tsys_eta0)**i
                end do
             end if
          end if
       end do

       if (n > tod%baseline_order+1) call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline)
    end do

    deallocate(x, y)

  end subroutine sample_chipass_baseline

  subroutine sample_chipass_baseline_per_slew(tod, scan, raw, s_tot, mask)
    !
    !   Sample CHIPASS baseline
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   scan:     local scan ID
    !   raw:      raw tod in du
    !   s_tot:    total signal model in mK
    !   mask:     list of accepted samples
    !
    implicit none
    class(comm_chipass_tod),                intent(inout) :: tod
    integer(i4b),                           intent(in)    :: scan
    real(sp), dimension(1:,1:),             intent(in)    :: raw, s_tot, mask

    integer(i4b) :: i, j, k, n, s, ntod_slew
    real(dp)     :: t, dt
    real(dp), allocatable, dimension(:) :: x, y

    do s = 1, tod%scans(scan)%nslew

       ntod_slew = tod%scans(scan)%slew_inds(s, 2) - tod%scans(scan)%slew_inds(s, 1) + 1
       allocate(x(ntod_slew), y(ntod_slew))
       dt = 1.d0 / ntod_slew

       !if ((tod%myid == 0) .and. (scan == 1)) write(*,*) 'slew, inds, n:', s, tod%scans(scan)%slew_inds(s, :), ntod_slew

       do j = 1, tod%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          t = 0.d0
          n = 0
          do k = tod%scans(scan)%slew_inds(s, 1), tod%scans(scan)%slew_inds(s, 2)
             t      = t + dt
             if (mask(k,j) > 0.5) then
                n    = n + 1
                x(n) = t
                y(n) = raw(k,j) - tod%scans(scan)%d(j)%gain * s_tot(k,j)
                if (allocated(tod%tsys_spline)) then
                   y(n) = y(n) - splint(tod%tsys_spline(j),tod%scans(scan)%d(j)%elev(k))
                else
                   do i = 0, tod%tsys_order
                      y(n) = y(n) - tod%tsys_fit(j, i) * (tod%scans(scan)%d(j)%elev(k) - tod%tsys_eta0)**i
                   end do
                end if
             end if
          end do
          if (n > tod%baseline_order+1) call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline_slew(s, :))
       end do

       deallocate(x, y)

    end do

  end subroutine sample_chipass_baseline_per_slew

  subroutine read_scan_inst_chipass(self, file, slabel, detlabels, scan)
    ! 
    ! Reads CHIPASS-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_chipass_tod)
    !           CHIPASS-specific TOD object
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
    class(comm_chipass_tod),             intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan

    integer(i4b) :: j

    if (self%read_elev) then
       do j = 1, self%ndet
          call read_hdf(file, slabel//'/'//trim(detlabels(j))//'/elev', scan%d(j)%elev)
       end do
    end if

    if (self%read_az) then
       do j = 1, self%ndet
          call read_hdf(file, slabel//'/'//trim(detlabels(j))//'/az', scan%d(j)%az)
       end do
    end if

    if (self%per_slew_baseline) then
       call read_hdf(file, slabel//'/common/nslew', scan%nslew)
       allocate(scan%slew_inds(scan%nslew,2))
       call read_hdf(file, slabel//'/common/slew_inds', scan%slew_inds)
    end if

  end subroutine read_scan_inst_chipass

  subroutine sample_chipass_Tsys(self, oper_default, handle)
    !
    !  Sample CHIPASS Tsys(eta)
    !
    !  Arguments:
    !  ----------
    !  self:     derived class (comm_chipass_tod)
    !            CHIPASS-specific TOD object
    !
    !  Returns:
    !  --------
    !  None, but updates TOD object
    !
    implicit none
    class(comm_chipass_tod), intent(inout) :: self
    integer(i4b),            intent(in)    :: oper_default
    type(planck_rng),        intent(inout) :: handle

    character(len=2) :: dtext
    type(comm_scandata) :: sd
    integer(i4b)        :: i, j, k, l, m, nout, ierr, bin, nbin
    real(dp)            :: t, dt, x, y, eta, del, el_min, el_max, tsys, alpha, var_scale
    real(dp), allocatable, dimension(:)   :: b, el, mu, sigma
    real(dp), allocatable, dimension(:,:) :: A

    allocate(A(0:self%tsys_order,0:self%tsys_order), b(0:self%tsys_order))
    if (.not. allocated(self%tsys_spline)) allocate(self%tsys_spline(self%ndet))

    del       = 1.d0  ! Delta elevation; bin width
    alpha     = 3.d4  ! Spline stiffness
    var_scale = 10.d0 ! Scaling factor between real and white noise variance for binned Tsys

    do j = 1, self%ndet
       el_min = self%el_min(j)
       el_max = self%el_max(j)
       nbin   = int((el_max-el_min)/del)+1
       allocate(mu(nbin), sigma(nbin), el(nbin))
       mu    = 0.d0
       sigma = 0.d0
       do i = 1, nbin
          el(i) = el_min + real(i-0.5,dp)*del ! Elevation in degrees
       end do

       ! Compute binned Tsys
       do i = 1, self%nscan
          if (.not. any(self%scans(i)%d%accept)) cycle
          call init_scan_data(self, i, oper_default, TODMASK_PROC, sd, spur_level=0)
          do k = 1, self%scans(i)%ntod
             if (sd%mask(k,j) == 1) then
                ! Prepare data 
                x   = self%scans(i)%d(j)%elev(k) - self%tsys_eta0
                y   = sd%tod(k,j) - self%scans(i)%d(j)%gain * sd%s_tot(k,j,0,1)
                bin = min(max(int((self%scans(i)%d(j)%elev(k)-el_min)/del),1),nbin)
                mu(bin)    = mu(bin)    +    y/self%scans(i)%d(j)%N_psd%sigma0**2
                sigma(bin) = sigma(bin) + 1.d0/self%scans(i)%d(j)%N_psd%sigma0**2
!!$                do l = 0, self%baseline_order
!!$                   y = y - self%scans(i)%d(j)%baseline(l) * t**l
!!$                end do
             end if
          end do
          call dealloc_scan_data(sd)
       end do
    
       ! Reduce binned Tsys
       if (self%myid == 0) then
          call mpi_reduce(mpi_in_place, mu, size(mu), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
          call mpi_reduce(mpi_in_place, sigma, size(sigma), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       else
          call mpi_reduce(mu, mu, size(mu), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
          call mpi_reduce(sigma, sigma, size(sigma), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       end if
    
       ! Sample splined Tsys
       if (self%myid == 0) then

          where(sigma > 0.d0)
             mu    = mu / sigma
             sigma = sqrt(1./sigma) ! Sigma
          elsewhere
             mu    = 0.d0
             sigma = 1.d10
          end where

          ! Method 1: Polynomial fit to binned Tsys
          A = 0.d0
          b = 0.d0
          do i = 1, nbin
             x = el(i) - self%tsys_eta0
             do l = 0, self%tsys_order
                b(l) = b(l) + mu(i) * x**l / sigma(i)**2
                do m = 0, self%tsys_order
                   A(l,m) = A(l,m) + x**l * x**m / sigma(i)**2
                end do
             end do
          end do
          
          ! Compute best-fit solution
          call invert_matrix(A, cholesky=.true.)
          self%tsys_fit(j,:) = matmul(A, b)

          ! Add random fluctuation
          call compute_hermitian_root(A, 0.5d0)
          do l = 0, self%tsys_order
             b(l) = rand_gauss(handle)
          end do
          self%tsys_fit(j,:) = self%tsys_fit(j,:) + matmul(A, b)

          ! Method 2: Smoothed spline with random variations
          do i = 1, nbin
             mu(i) = mu(i) + rand_gauss(handle) * sigma(i)
          end do
          call smooth_spline(el, mu, self%tsys_spline(j), "inv_var", alpha, var_scale*sigma**2)
          
          ! Output to ascii
          call int2string(j, dtext)
          open(58, file=trim(self%outdir)//'/Tsys'//dtext//'.dat', recl=128)
          do i = 1, nbin
             tsys = 0.
             do l = 0, self%tsys_order
                tsys = tsys + self%tsys_fit(j,l) * (el(i) - self%tsys_eta0)**l
             end do
             write(58,fmt='(f8.3,4e16.8)') el(i), mu(i), sigma(i), tsys, splint(self%tsys_spline(j), el(i))
          end do
          close(58)
       else
          if (.not. allocated(self%tsys_spline(j)%x)) then
             ! Initialize spline structures
             allocate(self%tsys_spline(j)%x(nbin),self%tsys_spline(j)%y(nbin),self%tsys_spline(j)%y2(nbin))
             self%tsys_spline(j)%boundary = 1d30
             self%tsys_spline(j)%regular  = .false.
             self%tsys_spline(j)%linear   = .false.
             self%tsys_spline(j)%verbose  = .false.
          end if
       end if

       ! Distribute new sample
       call mpi_bcast(self%tsys_fit(j,:), size(self%tsys_fit(j,:)), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
       call mpi_bcast(self%tsys_spline(j)%x,  nbin, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
       call mpi_bcast(self%tsys_spline(j)%y,  nbin, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
       call mpi_bcast(self%tsys_spline(j)%y2, nbin, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

       ! Clean up
       deallocate(el, mu, sigma)
    end do
    
    deallocate(A, b)

  end subroutine sample_chipass_Tsys

  subroutine construct_corrtemp_chipass(self, sd, det)
    !
    !  Construct a CHIPASS instrument-specific correction template; for now contains tsys
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !        CHIPASS-specific TOD object
    !  sd:   comm_scandata object
    !  det:  int
    !        detector number
    !
    !  Returns:
    !  --------
    !  None, but modifies sd
    !
    implicit none
    class(comm_chipass_tod), intent(in)          :: self
    class(comm_scandata), intent(inout)          :: sd
    integer(i4b),         intent(in),   optional :: det

    type(hdf_file)   :: file
    !character(len=6) :: slabel
    integer(i4b)     :: i, j, k, d, scan, s
    real(dp)         :: t, dt, t2, dt2
    !real(dp), allocatable, dimension(:) :: elev

    scan      = sd%scan
    dt        = 1.d0 / self%scans(scan)%ntod
    sd%s_inst = 0.d0
    do j = 1, sd%ndet
       d = j; if (present(det)) d = det
       t = 0.d0
       if (.not. self%scans(scan)%d(d)%accept) cycle
       s = 1
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          if (sd%mask(k,j) .le. 0.5) cycle
          if (self%per_slew_baseline) then
             if (k > self%scans(scan)%slew_inds(s, 2)) s = s + 1
                !if(s < self%scans(scan)%nslew) if (k == self%scans(scan)%slew_inds(s+1, 1)) s = s + 1
             dt2 = 1.d0 / (self%scans(scan)%slew_inds(s, 2) - self%scans(scan)%slew_inds(s, 1) + 1)
             if (k .ge. self%scans(scan)%slew_inds(s, 1)) then
                t2 = (k - self%scans(scan)%slew_inds(s, 1))*dt2
                if ((t2 < 0.d0) .or. (t2 > 1.d0)) write(*,*) 'bad t2; t2, k, s', t2, k, s
                do i = 0, self%baseline_order
                   sd%s_inst(k,j) = sd%s_inst(k,j) + self%scans(scan)%d(d)%baseline_slew(s, i) * t2**i
                end do
             end if
          else
             do i = 0, self%baseline_order
                sd%s_inst(k,j) = sd%s_inst(k,j) + self%scans(scan)%d(d)%baseline(i) * t**i
             end do
          end if
          if (allocated(self%tsys_spline)) then
             sd%s_inst(k,j) = sd%s_inst(k,j) + splint(self%tsys_spline(d),self%scans(scan)%d(d)%elev(k))
          else
             do i = 0, self%tsys_order
                sd%s_inst(k,j) = sd%s_inst(k,j) + self%tsys_fit(d,i) * (self%scans(scan)%d(d)%elev(k) - self%tsys_eta0)**i
             end do
          end if
       end do
    end do
  end subroutine construct_corrtemp_chipass

  subroutine initHDF_chipass(self, chainfile, path)
    ! 
    ! Initializes CHIPASS-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:      derived class (comm_chipass_tod)
    !            CHIPASS-specific TOD object
    ! chainfile: derived type (hdf_file)
    !            Already open HDF file handle to existing chainfile
    ! path:      string
    !            HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_chipass_tod),             intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    integer(i4b) :: ierr, i, j, k, ext(3), nbin
    real(dp), allocatable, dimension(:,:,:) :: baseline, buffer
    real(dp), allocatable, dimension(:,:,:,:) :: baseline_slew
    real(dp), allocatable, dimension(:,:)   :: tsys_fit
    real(dp) :: tsys_eta0

    allocate(baseline(self%nscan_tot, self%ndet, 0:self%baseline_order))
    allocate(baseline_slew(self%nscan_tot, self%ndet, self%max_nslew, 0:self%baseline_order))
    allocate(tsys_fit(self%ndet, 0:self%tsys_order))

    if (self%myid == 0) then
       call read_hdf(chainfile, trim(adjustl(path))//'baseline', baseline)
       if (hdf_group_exists(chainfile, trim(adjustl(path))//'baseline_per_slew')) then
          call read_hdf(chainfile, trim(adjustl(path))//'baseline_per_slew', baseline_slew)
       else
          baseline_slew = 0.d0
       end if
       call read_hdf(chainfile, trim(adjustl(path))//'tsys_fit', tsys_fit)
       call read_hdf(chainfile, trim(adjustl(path))//'tsys_eta0', tsys_eta0)

       call get_size_hdf(chainfile, trim(adjustl(path))//'tsys_spline', ext)
       allocate(buffer(ext(1), ext(2), ext(3)))
       call read_hdf(chainfile, trim(adjustl(path))//'tsys_spline', buffer)
    end if

    call mpi_bcast(ext, 3, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    if (self%myid /= 0) allocate(buffer(ext(1), ext(2), ext(3)))
    call mpi_bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    
    call mpi_bcast(baseline, size(baseline), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(baseline_slew, size(baseline_slew), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(tsys_fit, size(tsys_fit), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(tsys_eta0, 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

    self%tsys_eta0 = tsys_eta0
    do j = 1, self%ndet
       self%tsys_fit(j, :) = tsys_fit(j, :)
       do i = 1, self%nscan
          k = self%scanid(i)
          if (.not. self%scans(i)%d(j)%accept) cycle
          self%scans(i)%d(j)%baseline = baseline(k, j, :)
          self%scans(i)%d(j)%baseline_slew = baseline_slew(k, j, :, :)
       end do
    end do

    if (.not. allocated(self%tsys_spline)) allocate(self%tsys_spline(self%ndet))
    do j = 1, self%ndet
       nbin = count(buffer(:,1,j) /= -1.d30)
       allocate(self%tsys_spline(j)%x(nbin),self%tsys_spline(j)%y(nbin),self%tsys_spline(j)%y2(nbin))       
       self%tsys_spline(j)%x(1:nbin)  = buffer(1:nbin,1,j)
       self%tsys_spline(j)%y(1:nbin)  = buffer(1:nbin,2,j)
       self%tsys_spline(j)%y2(1:nbin) = buffer(1:nbin,3,j)
       self%tsys_spline(j)%boundary   = 1d30
       self%tsys_spline(j)%regular    = .false.
       self%tsys_spline(j)%linear     = .false.
       self%tsys_spline(j)%verbose    = .false.
    end do
    
    deallocate(baseline, baseline_slew, tsys_fit, buffer)

  end subroutine initHDF_chipass

  subroutine dumpToHDF_chipass(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:      derived class (comm_tod)
    !            TOD object
    ! chainfile: derived type (hdf_file)
    !            Already open HDF file handle to existing chainfile
    ! path:      string
    !            HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_chipass_tod), intent(in)     :: self
    type(hdf_file),          intent(in)     :: chainfile
    character(len=*),        intent(in)     :: path

    integer(i4b) :: i, j, k, nbin, ierr
    real(dp), allocatable, dimension(:,:,:) :: buffer
    real(dp), allocatable, dimension(:,:,:,:) :: baseline_slew

    if (self%per_slew_baseline) then
       allocate(baseline_slew(self%nscan_tot, self%ndet, self%max_nslew, 0:self%baseline_order))
       do j = 1, self%ndet
          do i = 1, self%nscan
             k = self%scanid(i)
             if (.not. self%scans(i)%d(j)%accept) cycle
             baseline_slew(k, j, :, :) = self%scans(i)%d(j)%baseline_slew
          end do
       end do
       if (self%myid == 0) then
          call mpi_reduce(mpi_in_place, baseline_slew, size(baseline_slew), &
                & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       else
          call mpi_reduce(baseline_slew, baseline_slew, size(baseline_slew), &
                & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       end if
    end if
    
    if (self%myid == 0 .and. trim(self%level) == 'L1') then
       call write_hdf(chainfile, trim(adjustl(path))//'tsys_eta0', self%tsys_eta0)
       call write_hdf(chainfile, trim(adjustl(path))//'tsys_fit', self%tsys_fit)

       if (allocated(self%tsys_spline)) then
          nbin = 0
          do i = 1, self%ndet
             nbin = max(nbin,size(self%tsys_spline(i)%x))
          end do
          allocate(buffer(nbin,3,self%ndet))
          buffer = -1.d30
          do i = 1, self%ndet
             nbin = size(self%tsys_spline(i)%x)
             buffer(1:nbin,1,i) = self%tsys_spline(i)%x
             buffer(1:nbin,2,i) = self%tsys_spline(i)%y
             buffer(1:nbin,3,i) = self%tsys_spline(i)%y2
          end do
          call write_hdf(chainfile, trim(adjustl(path))//'tsys_spline', buffer)
          deallocate(buffer)
       end if

       if (self%per_slew_baseline) then
          call write_hdf(chainfile, trim(adjustl(path))//'baseline_per_slew', baseline_slew)
          deallocate(baseline_slew)
       end if
    end if

  end subroutine dumpToHDF_chipass

end module comm_tod_chipass_mod
