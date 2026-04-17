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
submodule (comm_tod_akari_mod) comm_tod_akari_smod
contains

  !**************************************************
   !             Constructor
   !**************************************************
  module function constructor_akari(cpar, id, id_abs, info, tod_type) result(c)
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
      class(comm_akari_tod),   pointer    :: c

      integer(i4b) :: i, j, nside_beam, lmax_beam, nmaps_beam, ierr
      logical(lgt) :: pol_beam
      real(dp)     :: band_samprate

      call timer%start(TOD_INIT, id_abs)

      ! Allocate object
      allocate(c)

      ! Set up noise PSD type and priors
      c%freq            = cpar%ds_label(id_abs)
      c%noise_psd_model = '2oof'

      if (trim(c%noise_psd_model) == 'oof') then
         c%n_xi            = 3
      else if (trim(c%noise_psd_model) == '2oof') then
         c%n_xi            = 5
      end if 

      allocate(c%xi_n_nu_fit(c%n_xi,2))
      allocate(c%xi_n_P_uni(c%n_xi,2))
      allocate(c%xi_n_P_rms(c%n_xi))
      
      ! For 2oof, the max psd freq. (samprate/2) is needed for freq. ranges and priors.
      ! samprate is defined in c%read_tod(c%label), but psd params must be defined before this. 
      if (trim(c%freq) == 'AKARI_N160') then
         band_samprate = 16.86 ! Hz
      else if (trim(c%freq) == 'AKARI_WIDE-L') then
         band_samprate = 16.86 ! Hz 
      else if (trim(c%freq) == 'AKARI_WIDE-S') then
         band_samprate = 25.28 ! Hz
      else if (trim(c%freq) == 'AKARI_N60') then
         band_samprate = 25.28 ! Hz
      end if

      ! Correlated noise parameters
      c%xi_n_nu_fit(1,:) = [0.05d0, 0.99*band_samprate/2]   ! Freq range for sigma0
      c%xi_n_nu_fit(2,:) = [0.d0, 0.5d0]    ! Freq range for fknee
      c%xi_n_nu_fit(3,:) = [0.d0, 0.5d0]    ! Freq range for alpha
      
      c%xi_n_P_uni(1,:)  = [1d-6, 100.d0]   ! Uniform prior for sigma0
      c%xi_n_P_uni(2,:)  = [6d-5, 1.d0]     ! Uniform prior for fknee
      c%xi_n_P_uni(3,:)  = [-4d0, -0.5d0]   ! Uniform prior for alpha
      
      ! Set rms of all parameters to 0.05 for initial test phase. 
      c%xi_n_P_rms(1)    = 1.00d0           ! Prior rms for sigma0
      c%xi_n_P_rms(2)    = 0.05d0           ! Prior rms for fknee
      c%xi_n_P_rms(3)    = 0.05d0           ! Prior rms for alpha

      if (trim(c%noise_psd_model) == '2oof') then
         ! Extend correlated noise parameters for fknee2 and alpha2
         c%xi_n_nu_fit(4,:) = [4.d0, 0.99*band_samprate/2]  ! Freq range for fknee2
         c%xi_n_nu_fit(5,:) = [4.d0, 0.99*band_samprate/2]  ! Freq range for alpha2
         
         c%xi_n_P_uni(4,:)  = [4.0d0, 0.99*band_samprate/2]  ! Uniform prior for fknee2
         c%xi_n_P_uni(5,:)  = [0.5d0, 3.d0]   ! Uniform prior for alpha2
         
         c%xi_n_P_rms(4)    = 0.05d0         ! Prior rms for fknee2
         c%xi_n_P_rms(5)    = 0.05d0         ! Prior rms for alpha2
      end if

      ! Data selection parameters
      c%chisq_threshold  = 1000d0       ! Cut scans with higher chisq
      c%sigma0_threshold = 1.d0        ! Cut scans with higher sigma0

      ! Initialize common parameters
      call c%tod_constructor(cpar, id, id_abs, info, tod_type)

      ! Initialize instrument-specific parameters
      !read(c%freq(1:2),*) c%zodiband
      c%samprate_lowres = 8.  ! Lowres samprate in Hz
      c%nhorn           = 1
      c%ndiode          = 1
      if (trim(c%level) == 'L1') then
          c%compressed_tod  = .true.
      else
          c%compressed_tod  = .false.
      end if
      c%correct_sl      = .false.                ! OBS FIXME these flags are obsolete and not used
      c%correct_orb     = .false.                ! replaced by oper_default below
      c%orb_4pi_beam    = .false.
      c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
      c%apply_inst_corr = .false.   
      c%symm_flags      = .false.
      c%use_moon_point  = cpar%sample_moon_maps
      c%use_earth_elon  = cpar%sample_earth_maps
      ! c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      ! c%chisq_threshold = 50000.
      c%nside_pixhist    = 64

      c%nmaps           = info%nmaps
      if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
         c%ndet         = count_detectors(cpar%ds_tod_dets(id_abs)) !, cpar%datadir)
      else
         c%ndet         = num_tokens(cpar%ds_tod_dets(id_abs), ",")
      end if
            
      nside_beam                  = 128
      nmaps_beam                  = 1
      pol_beam                    = .false.
      c%nside_beam      = nside_beam

      ! Get detector labels
      if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
         call get_detectors(cpar%ds_tod_dets(id_abs), c%label)
      else
         call get_tokens(trim(adjustl(cpar%ds_tod_dets(id_abs))), ",", c%label)
      end if
      
      ! Read the actual TOD
      call c%read_tod(c%label)

      ! Initialize bandpass mean and proposal matrix
      call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
      c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false.)
      call c%precompute_lookups()

      ! Load the instrument file
      call c%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

   ! commenting this out for now
      !do i=1, c%ndet
      !  call init_noise_model(c, i)
      !end do

      ! Collect Sun velocities from all scans
      call c%collect_v_sun
      
      ! Allocate sidelobe convolution data structures
      allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)

      c%orb_dp => comm_orbdipole(comm=c%comm)

      ! Initialize all baseline corrections to zero
      !do i = 1, c%nscan
      !   c%scans(i)%d%baseline = 0.d0
      !end do

      ! Initialize dynamic mask
      c%dynmask => comm_dynmask(c, cpar)
      c%dynmask%output_scan             = 1000
      c%dynmask%output_det              = 1
      c%dynmask%apply_pixhist           = .true.
      c%dynmask%remove_isolated_samples = .true.
      c%dynmask%threshold_longchunks    = 0.3
      c%dynmask%window_longchunks       = 2000
      c%dynmask%threshold_cr            = 6.   ! sigma
      c%dynmask%width_cr_mask           = 5
      allocate(c%dynmask%templ_cr(-5:5))
      c%dynmask%templ_cr    = 1.
      c%dynmask%templ_cr(0) = 10.
      c%dynmask%templ_cr(1) = 5.
      c%dynmask%templ_cr(2) = 2.5
      c%dynmask%templ_cr(3) = 1.25
      c%dynmask%templ_cr(4) = 1.
      c%dynmask%templ_cr    = c%dynmask%templ_cr - mean(real(c%dynmask%templ_cr,dp))
      

!!$      if (trim(self%freq) == 'AKARI_N160') then
!!$         flag_threshold        = [1.0,    -20., -5.0,     -3.0, -2.0, -1.5,     1.0,      0.3,     -1.0]         
!!$      else if (trim(self%freq) == 'AKARI_WIDE-L') then
!!$         flag_threshold        = [1.0,    -20., -5.0,     -3.0, -2.0, -1.5,     1.0,      0.3,     -1.0]
!!$      else if (trim(self%freq) == 'AKARI_WIDE-S') then
!!$         flag_threshold        = [1.0,    -20., -5.0,     -3.0, -2.0, -1.5,     1.0,      0.3,     -1.0]
!!$      else if (trim(self%freq) == 'AKARI_N60') then
!!$         flag_threshold        = [1.0,    -20., -5.0,     -3.0, -2.0, -1.5,     1.0,      0.3,     -1.0]
!!$      end if

      ! initialize ramp correction templates
      c%ntempl = 1     ! OBS may change
      allocate(c%nsamp_templ(c%ntempl))
      if (trim(c%freq) == 'AKARI_N160' .or. trim(c%freq) == 'AKARI_WIDE-L') then
         c%nsamp_templ = [28]
      else if (trim(c%freq) == 'AKARI_N60' .or. trim(c%freq) == 'AKARI_WIDE-S') then
         c%nsamp_templ = [44]
      end if
      allocate(c%tod_correction_templ(maxval(c%nsamp_templ),c%ntempl,c%ndet,c%nscan))
      c%tod_correction_templ = 0.d0

      
      call timer%stop(TOD_INIT, id_abs)
    end function constructor_akari

   !**************************************************
   !             Driver routine
   !**************************************************
   module subroutine process_akari_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      ! 
      ! Routine that processes the akari Calibrated Individual Observations. 
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
      ! Writes maps to disc in fits format
      ! 
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_akari_tod class
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
      class(comm_akari_tod),                    intent(inout) :: self
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
      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist, sample_zodi, use_k98_samp_groups, output_zodi_comps, sample_ncorr, only_solar_mask, sample_xi_n, sample_ramp
      type(comm_binmap)   :: binmap
      type(comm_scandata) :: sd
      character(len=4)    :: ctext, myid_text
      character(len=2)    :: zodi_param_text
      character(len=1)    :: up_down_text
      character(len=6)    :: samptext, scantext
      character(len=512)  :: prefix, postfix, prefix4D, prefix_atlas, postfix_atlas
      character(len=512), allocatable, dimension(:) :: slist
      real(sp),              dimension(9)       :: flag_threshold
      real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi
      real(sp), allocatable, dimension(:,:,:)   :: d_calib
      real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
      real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)    :: A_T_A, A_T_A_reduced
      real(dp), allocatable, dimension(:)       :: AY, AY_reduced, X
      real(dp), allocatable, dimension(:, :, :) :: s_therm_tot, s_scat_tot ! (n_tod_tot, ncomps, ndet)
      real(dp), allocatable, dimension(:, :)    :: res_tot ! (n_tod_tot, ndet)
      real(dp), allocatable, dimension(:)       :: mask_tot ! (n_tod_tot)
      type(hdf_file) :: tod_file

      integer(i4b), dimension(1) :: col_def = [1]  ! Only T
      class(comm_cgmap), pointer :: cgmap

      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)

      !call timer%start(TOD_ALLOC, self%band)

      ! Output input sky model
      call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
      
      ! Toggle optional operations
      sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
      output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
      use_k98_samp_groups   = .false.                          ! fits one overall albedo and episolon for the dust bands, and one for ring + feature
      sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                         ! don't sample absolute bandpasses
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
      sample_gain           = iter > 1                        ! Gain sampling, LB TOD sims have perfect gain
      only_solar_mask       = .false.                        ! Only apply solar mask

      if (.false.) then
         ! Debug
         select_data     = iter == 2
         sample_ncorr    = iter == 2
         sample_xi_n     = iter == 2
         sample_ramp     = iter == 2
      else if (trim(self%init_from_HDF) == 'none') then ! OBS FIXME bug when BAND_TOD_INI_:FROM_HDF=default and INIT_CHAIN=none
         select_data     = iter == 5                    ! in param file. Takes you to 'else' below
         sample_ncorr    = iter > 2
         sample_xi_n     = iter > 2
         sample_ramp   = iter > 3
      else
         select_data     = iter == 3 
         sample_ncorr    = iter > 1
         sample_xi_n     = iter > 1
         sample_ramp   = iter > 2
      end if

      oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
           & SD_SKY,SD_INST,SD_NCORR,SD_BP,SD_ORB,SD_ZODI])
      
      ! Initialize local variables
      ndelta          = size(delta,3)
      self%n_bp_prop  = ndelta-1
      nside           = map_out%info%nside
      nmaps           = map_out%info%nmaps
      npix            = 12*nside**2
      self%output_n_maps = 8
      if (output_zodi_comps) self%output_n_maps = self%output_n_maps + zodi_model%n_comps

      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'
      prefix_atlas = trim(chaindir) // '/atlas_' // trim(self%freq) // '_' // trim(zodi_param_text) // '_' // trim(up_down_text) // '_'
      postfix_atlas = '.fits'

      ! Initialize index-based sky map and mask
      call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
      call update_status(status, "tod_cache"//ctext)


      ! Write mask for debugging
      if (.false. .and. self%myid == 0) then
         print *, "writing masks"
         call open_hdf_file(trim(chaindir)//'/mask.h5', tod_file, 'w')
         call write_hdf(tod_file, '/procmask', procmask)
         call write_hdf(tod_file, '/procmask2', procmask2)
         call write_hdf(tod_file, '/procmask_zodi', procmask_zodi)
         call close_hdf_file(tod_file)
         stop
      end if


      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      ! Create pixel histograms
      if (self%first_call) call compute_tod_pixhist(self)
      
      ! Sample gain components in separate TOD loops; marginal with respect to n_corr
      if (sample_gain) then
         ! 'abscal': the global constant gain factor
         call sample_calibration(self, 'abscal', oper_default, handle)
         ! 'relcal': the gain factor that is constant in time but varying between detectors
         call sample_calibration(self, 'relcal', oper_default, handle)
         ! 'deltaG': the time-variable and detector-variable gain
         !call sample_calibration(self, 'deltaG', oper_default, handle)
      end if

      call timer%start(TOD_INSTCORR, self%band)
!!$      ! initialize ramp sampling
!!$      if (self%first_call) then
!!$         ! set number of correction templates
!!$         self%ntempl = 1     ! OBS may change
!!$         ! allocate self%nsamp_templ
!!$         allocate(self%nsamp_templ(self%ntempl, self%ndet, self%nscan))
!!$         do i = 1, self%nscan
!!$            ! skip scan if no accepted data
!!$            if (.not. any(self%scans(i)%d%accept)) cycle
!!$            ! decompress data for this scan into sd
!!$            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, spur_level=0) ! OBS may change to decompress_flags
!!$            ! find length of ramp correction template(s)
!!$            call self%init_sample_ramp(sd)
!!$            ! clean up
!!$            call dealloc_scan_data(sd)
!!$         end do
!!$         call timer%stop(TOD_INSTCORR, self%band)
!!$
!!$         if (self%myid==0) then
!!$            write(*,*) 'max templ length', maxval(self%nsamp_templ)
!!$            open(19, file="maxlength_scan.dat", recl=1024)
!!$            do i = 1, self%nscan
!!$               !do j = 1, self%ndet
!!$               if (maxval(self%nsamp_templ(1, :, i)) /= 28) then
!!$                  write(19,*) self%scanid(i), maxval(self%nsamp_templ(1, :, i))
!!$                  ! write(19,*) self%scanid(i), self%nsamp_templ(1, j, i)
!!$               end if
!!$               !end do
!!$            end do
!!$            close(19)
!!$            write(*,*) 'written to file'
!!$         end if
!!$         call timer%start(TOD_INSTCORR, self%band)
!!$         ! allocate and initialize self%tod_correction_templ
!!$         allocate(self%tod_correction_templ(maxval(self%nsamp_templ),self%ntempl,self%ndet,self%nscan))
!!$         self%tod_correction_templ = 0.d0
!!$      end if
      ! ramp sampling
      if (sample_ramp) then
          do i = 1, self%nscan
            ! skip scan if no accepted data
            if (.not. any(self%scans(i)%d%accept)) cycle
            ! decompress data for this scan into sd
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, spur_level=2)  ! OBS may change mask
            ! sample tod correction template for ramp
            call self%sample_ramp(sd)
            ! clean up
            call dealloc_scan_data(sd)
         end do
      end if
      call timer%stop(TOD_INSTCORR, self%band)

      ! Create dynamic mask
      if (select_data) then
         if (self%myid == 0) write(*,*) '   --> Creating dynamic mask'
         do i = 1, self%nscan
            ! Skip scan if no accepted data
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_GAIN, sd) ! Should be TODMASK_FLAG
            do j = 1, sd%ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               call self%dynmask%create(sd, j)
            end do
            call dealloc_scan_data(sd)
         end do
         ! Synchronize and output flagging statistics in first iteration
         call self%dynmask%report
      end if
      
      ! Prepare mapmaking data structures
      call binmap%init(self, .true., sample_rel_bandpass)
      if (sample_abs_bandpass .or. sample_rel_bandpass) then
         allocate(chisq_S(self%ndet,size(delta,3)))
         chisq_S = 0.d0
      end if
      if (output_scanlist) then
         allocate(slist(self%nscan))
         slist   = ''
      end if

      ! Initialize CG mapmaker, maptype = 1 = T-only
      cgmap => comm_cgmap(self, 1)
      
      ! Perform loop over scans to prepare data to make maps
      if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
      do i = 1, self%nscan

         ! Skip scan if no accepted data
         !write(*,*) i, self%scans(i)%d%accept
         if (.not. any(self%scans(i)%d%accept)) cycle
         call wall_time(t1)
         ! FIXME!! NCORR MASK CHANGED TEMPORARILY TO GAIN MASK. WILL LEAK SIGNAL INTO NCORR!!!
         call init_scan_data(self, i, oper_default, TODMASK_GAIN, sd)

          !if (self%myid == 0 .and. i == 1) then
         !    open(58,file='tod.dat', recl=1024)
         !    write(*,*)  "debug", shape(sd%tod), "debug"
         !    do j = 1, sd%ntod
         !       write(58,*) j, sd%tod(j,1), sd%mask(j,1), sd%flag(j,1), sd%tod(j,2), sd%mask(j,2), sd%flag(j,2), sd%tod(j,3), sd%mask(j,3), sd%flag(j,3), sd%tod(j,4), sd%mask(j,4), sd%flag(j,4), sd%tod(j,5), sd%mask(j,5), sd%flag(j,5)
         !    end do
         !    close(58)
         ! end if
         
         ! Apply fast flags
         call update_status(status, "quick_tod_flag_"//ctext)
         
         ! Sample correlated noise
         !call project_mask(self, TODMASK_ZODI, sd)
         if (sample_ncorr) then
            call sample_n_corr(self, sd, handle)
            if (sample_xi_n) then
               call sample_noise_psd(self, sd, handle, chaindir)
            else
               call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true.)
            end if
         else
            call sample_n_corr(self, sd, handle, onlymono=.true.)
            call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true.)
         end if

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_tod_chisq(sd, j)
         end do

         ! Select data
         if (select_data) call remove_bad_data(self, i, sd%flag)

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
         
         ! Compute binned map
         allocate(d_calib(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib = 0.d0
         call compute_calibrated_data(self, i, sd, d_calib)    

         ! Feed CG mapmaker calibrated and cleaned data
         call cgmap%load_data(i, self, sd, d_calib(1,:,:))
                  
         !write(*,*) "Scan = ", self%scanid(i), ', num moon = ', count(iand(sd%flag,2)==2)
         
         ! For debugging: write TOD to hdf
         if (.true.) then
            ! scan id appears to be the worst chi2
            if (self%scanid(i) == 915) then 
               !print *, self%scanid(i)
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//'_'//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call write_hdf(tod_file, '/flag', sd%flag)
               call write_hdf(tod_file, '/todz', d_calib(1, :, :))
               call write_hdf(tod_file, '/s_sky', sd%s_sky)
               call write_hdf(tod_file, '/s_inst', sd%s_inst)
               call write_hdf(tod_file, '/n_corr', sd%n_corr)
               call write_hdf(tod_file, '/s_sl', sd%s_sl)
               call write_hdf(tod_file, '/s_orb', sd%s_orb)
               call write_hdf(tod_file, '/res', d_calib(2, :, :))
               call write_hdf(tod_file, '/zodi', d_calib(7, :, :))
               call write_hdf(tod_file, '/mask', sd%mask)
!               call write_hdf(tod_file, '/accept', real(self%scans(i)%d(:)%accept,dp))
               call write_hdf(tod_file, '/sigma0', self%scans(i)%d(1)%N_psd%sigma0)
               call write_hdf(tod_file, '/gain', self%scans(i)%d%gain)
               call close_hdf_file(tod_file)
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
         call dealloc_scan_data(sd)
         deallocate(d_calib)
      end do

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
         call sample_bp(self, handle, chisq_S, delta)
         self%bp_delta = delta(:,:,1)
      end if


      ! Output maps to disk
      if (.false.) then
         ! call map_out%writeFITS(trim(prefix_atlas)//'map'//trim(postfix_atlas))
         ! call rms_out%writeFITS(trim(prefix_atlas)//'rms'//trim(postfix_atlas))
         ! if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix_atlas)//'res'//trim(postfix_atlas))
         ! if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix_atlas)//'ncorr'//trim(postfix_atlas))
         if (self%output_n_maps > 6 .and. self%subtract_zodi) call binmap%outmaps(7)%p%writeFITS(trim(prefix_atlas)//'zodi'//trim(postfix_atlas))
         if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
            do i = 1, zodi_model%n_comps
               call binmap%outmaps(8+i)%p%writeFITS(trim(prefix_atlas)//'zodi_'//trim(zodi_model%comp_labels(i))//trim(postfix_atlas))
            end do
         endif
      else
         call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
         call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
         if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
         if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
         if (self%output_n_maps > 6 .and. self%subtract_zodi) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
         if (self%output_n_maps > 7 ) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'inst'//trim(postfix))
         if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
            do i = 1, zodi_model%n_comps
               call binmap%outmaps(8+i)%p%writeFITS(trim(prefix)//'zodi_'//trim(zodi_model%comp_labels(i))//trim(postfix))
            end do
         endif
      end if
      ! if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
      !    do i = 1, zodi%n_comps
      !       call binmap%outmaps(8+i)%p%writeFITS(trim(prefix)//'zodi_'//trim(zodi_comp_names(i))//trim(postfix))
      !    end do
      ! endif

      ! Solve for CG map
      call cgmap%solve(map_out)
      call dealloc_cgmap(cgmap)
      call map_out%writeFITS(trim(prefix)//'cgmap'//trim(postfix))
      !call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

      ! Clean up
      call binmap%dealloc()
      if (allocated(slist)) deallocate(slist)
      !  if (self%correct_sl) then
      !     do i = 1, self%ndet
      !        call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
      !     end do
      !  end if

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

      call update_status(status, "tod_end"//ctext)
      
      call timer%stop(TOD_TOT, self%band)

   end subroutine process_akari_tod   


   module subroutine apply_fast_flags_akari(self, sd)
     !  Apply fast flags to sd%flag; should only depend on time, pix or flag arrays, not TOD
     !  Expensive operations should instead be added to the dynamic mask
     !
     !  Arguments:
     !  ----------
     !  self: comm_tod object
     !
     implicit none
     class(comm_akari_tod),                 intent(inout)    :: self
     class(comm_scandata),                  intent(inout)    :: sd


     integer(i4b) :: i, j, k


     ! Exclude a sample by setting bit 29 in sd%flag to 1,
     ! ie., sd%flag = sd%flag + 536870912

     ! Nils and Katrine -- Add all fast flagging operations here (like getting rid of the 2 second ramp reset etc.)
     ! If things need to be precomputed, and require more expensive analyses, it's better to add it to the dynamic mask,
     ! which is only computed once. However, that's stored in an (nmask,2) array format, giving the samples of the start and
     ! end of a rejected TOD segment. We clearly do not want the 2 second thing in there, as it would require a lot of memory
     ! just to store it. So, use this routine for things like the 2 second mask, which may be computed quickly on the fly,
     ! but the dynamic mask for things that are expensive.
     !
     ! Also note that the actual flag array should not change from Gibbs sample to Gibbs sample (at least not after some burn-in), and
     ! so the flags in this array must be deterministic between iterations
     
     ! Masking the 2 second ramp reset after each calibration lamp flash


      do i = 1, sd%ndet
         do j = 1, sd%ntod
            if (btest(sd%flag(j,i), TOD_CALLAMP1)) then
               ! Mask the next 150 samples (~4 seconds at 24 Hz) after cal lamp 1
               do k = j, min(j+150, sd%ntod)
                  sd%flag(k,i) = ibset(sd%flag(k,i), 29)
               end do
            end if

            if (btest(sd%flag(j,i), TOD_CALLAMP2)) then
               ! Mask the next 150 samples (~4 seconds at 24 Hz) after cal lamp 2
               do k = j, min(j+150, sd%ntod)
                  sd%flag(k,i) = ibset(sd%flag(k,i), 29)
               end do
            end if

            ! if (btest(sd%flag(j,i), TOD_RAMP_RESET)) then
            !    ! Mask the next 2 samples after ramp reset
            !    do k = j, min(j+30, sd%ntod)
            !       sd%flag(k,i) = ibset(sd%flag(k,i), 29)
            !    end do   
            ! end if
            
         end do
      end do

   end subroutine apply_fast_flags_akari

   module subroutine construct_corrtemp_akari(self, sd, det)
     ! construct AKARI instrument-specific correction template from ramp template
     ! put into sd%s_inst
     ! ramp template is calculated in sample_ramp routine below and stored in self%tod_correction_templ
     implicit none
     class(comm_akari_tod), intent(in)             :: self ! full data object with compressed tods
     class(comm_scandata),  intent(inout)          :: sd   ! decompressed self%scans(scan)%d(det)%tod
     integer(i4b),          intent(in),   optional :: det

     integer(i4b) :: i, j, k, l, m, ndet, scan, ntod, nramps, nsamp, ntempl, beg_det, end_det
     integer(i4b) :: flag_id(5), ndata, nflagged, start
     logical(lgt) :: flagged
     real(dp),       allocatable, dimension(:)     :: templ!, res

     ! flag bit 5 = ramp reset, 15 = cal lamp on for long wavelengts, 16 = cal lamp for short wavelength
     flag_id(1) =  5
     flag_id(2) =  5
     flag_id(3) =  5
     flag_id(4) = 15
     flag_id(5) = 16
    
     scan   = sd%scan              ! scan number
     ndet   = sd%ndet              ! number of detectors                = self%ndet
     ntod   = sd%ntod              ! number of tod samples for my scan  = self%scans(scan)%ntod
     ntempl = self%ntempl          ! number of tod correction templates to include

     !allocate(res(1:ntod))

     ! start by putting instrument signal tod to zero   ! OBS assuming s_inst not used for other systematics!
     sd%s_inst = 0.d0

     ! check whitch detectors to include
     if (present(det)) then
        beg_det = det
        end_det = det
     else
        beg_det = 1
        end_det = ndet
     end if

     ! loop over detetctors
     do m = beg_det, end_det
        !if (self%myid==0) write (*,*) 'detector nr', m, ' of ', ndet
        ! skip detector if full scan is un-accepted
        if (.not. self%scans(scan)%d(m)%accept) cycle

        ! loop through samples and find unfagged segments of length nsamp, following ramp flags
        do k = 1, ntempl
           !if (self%myid==0) write (*,*) 'templ nr', k, ' of ', ntempl
           ! Put correction template with length nsamp into templ
           nsamp     = self%nsamp_templ(k)                         ! length of give correction template (ramp event)
           allocate(templ(nsamp))
           templ     = self%tod_correction_templ(1:nsamp,k,m,scan) ! tod correction template
           ! skip unwrapping if  correction template is zero
           if (maxval(templ)==0.d0 .and. minval(templ)==0.d0) then
              deallocate(templ)
              cycle
           end if
           nramps    = 0                                           ! counting number of ramp events 
           start     = 0                                           ! start of new ramp event in tod
           ndata     = 0                                           ! counting unflagged samples
           nflagged  = 0                                           ! counting flagged samples
           flagged   = .false.
           do i = 1, ntod
              !if (self%myid==0) write (*,*) 'sample nr', i, start, nflagged, nramps
              if ( btest(sd%flag(i,m), flag_id(k) ) ) then       ! if tod sample has ramp flag
                 if (.not. flagged) then                         ! .. and previous sample was unflagged
                    ! if this is not start of first event
                    if (start /= 0) then
                       ! check that previous event had correct length
                       if (ndata == nsamp) then
                          !if (self%myid==0) write (*,*) 'B start', start+nflagged, start, nflagged, nsamp
                          ! put template into s_inst tod
                          sd%s_inst(start+nflagged:start+nflagged+nsamp-1, m) = templ
                          nramps = nramps + 1
                       end if
                    end if
                    start = i                                     ! .. then start new event
                    ndata = 0                                     ! .. and restart unflagged counter
                    nflagged = 0                                  ! .. andand restart flagged counter
                 end if
                 flagged = .true.                                 ! irrespective of whether last sample was flagged
                 nflagged = nflagged +1                           ! this one is, and flag counter increases.   
              else                                                ! else if tod sample is unflagged,
                 flagged = .false.                                ! unflagge data counter increases.
                 ndata = ndata + 1
              end if
           end do
           ! output for plotting
           !if (self%myid==0) write(*,*) 'B num ramps for template:', nramps, k, scan
           deallocate(templ)
        end do

        ! get tod data for residual = full tod - signal_tod ! OBS ignoring ncorr which is not set yet
        !res = sd%tod(:,m) - sd%s_tot(:,m,0,1) * self%scans(scan)%d(m)%gain
        !if (self%myid==0) then
        !   open(19, file="res_corr.dat", recl=1024)
        !   do l = 1, ntod
        !      write(19,*) l, res(l), sd%s_inst(l,1), res(l) - sd%s_inst(l,1)
        !   end do
        !   close(19)
        !end if

     end do
     !deallocate(res)
   end subroutine construct_corrtemp_akari

   module subroutine sample_ramp(self, sd)
     ! find template of baseline beween to ramp events, by averaging over all events in a scan
     ! Put this into self%tod_correction_templ
     implicit none
     class(comm_akari_tod), intent(inout)          :: self  ! full data object with compressed tods
     class(comm_scandata),  intent(in)             :: sd    ! decompressed self%scans(scan)%d(det)%tod

     integer(i4b) :: i, j, k, l, m, ndet, scan, ntod, nramps, nsamp, ntempl
     integer(i4b) :: flag_id(5), ndata, nflagged, start
     real(dp)     :: flag_val(3) 
     logical(lgt) :: flagged
     real(dp),       allocatable, dimension(:)     :: res, total, buffer, mask

     ! flag bit 5 = ramp reset, 15 = cal lamp on for long wavelengts, 16 = cal lamp for short wavelength
     flag_id(1) =  5
     flag_id(2) =  5
     flag_id(3) =  5
     flag_id(4) = 15
     flag_id(5) = 16
     ! output value for plotting
     flag_val(1) = 20.d0
     flag_val(2) = 30.d0
     flag_val(3) = 40.d0
    
     scan   = sd%scan              ! scan number
     ndet   = sd%ndet              ! number of detectors                = self%ndet
     ntod   = sd%ntod              ! number of tod samples for my scan  = self%scans(scan)%ntod
     ntempl = self%ntempl          ! number of tod correction templates we will generate

     allocate(res(1:ntod))

     ! loop over detectors
     do m = 1, ndet
        ! skip detector if full scan is un-accepted
        if (.not. self%scans(scan)%d(m)%accept) cycle
        ! get tod data for residual = full tod - signal_tod ! OBS ignoring ncorr which is not set yet
        res = sd%tod(:,m) - sd%s_tot(:,m,0,1) * self%scans(scan)%d(m)%gain

        ! loop through tod samples and find unflagged segments of length nsamp, following ramp flags
        do k = 1, ntempl
           nsamp     = self%nsamp_templ(k)                       ! length of give correction template (ramp event)
           if (nsamp == 0) cycle
           allocate(buffer(nsamp), total(nsamp), mask(nsamp))
           total     = 0.d0                                      ! will collect all ramp segements here
           nramps    = 0                                         ! counting number of ramp events 
           start     = 0                                         ! start of new ramp event in tod
           ndata     = 0                                         ! counting unflagged samples
           nflagged  = 0                                         ! counting flagged samples
           flagged   = .false.
           do i = 1, ntod
              !if (self%myid==0) write (*,*) 'sample nr', i, start, nflagged, nramps
              if ( btest(sd%flag(i,m), flag_id(k) ) ) then       ! if tod sample has ramp flag
                 if (.not. flagged) then                         ! .. and previous sample was unflagged
                    ! if this is not start of first event
                    if (start /= 0) then
                       ! check that previous event had correct length
                       if (ndata == nsamp) then                    
                          ! put previous event into buffer
                          buffer =  sd%tod(start+nflagged:start+nflagged+nsamp-1, m)
                          ! and get corresponding mask
                          mask   = sd%mask(start+nflagged:start+nflagged+nsamp-1, m)
                          ! if enough unmasked sample, add buffer to total sum, and count +1 ramp event
                          if ( sum(mask) > nsamp*0.5 ) then
                             total = total + buffer * mask
                             nramps = nramps + 1
                          end if
                       end if
                    end if
                    start = i                                     ! .. then start new event
                    ndata = 0                                     ! .. and restart unflagged counter
                    nflagged = 0                                  ! .. andand restart flagged counter
                 end if
                 flagged = .true.                                 ! irrespective of whether last sample was flagged
                 nflagged = nflagged +1                           ! this one is, and flag counter increases.   
              else                                                ! else if tod sample is unflagged,
                 flagged = .false.                                ! unflagge data counter increases.
                 ndata = ndata + 1
              end if
           end do
           if (nramps /= 0) then
              ! average over all ramp events to create template
              total = total/real(nramps,dp)
              ! make sure template has zero mean
              total = total - mean(total)
           end if
           ! put template into comm_akari_tod data structure
           self%tod_correction_templ(1:nsamp,k,m,scan) = total

           ! output for plotting
           !if (self%myid==0) then
           !   write(*,*) 'A num ramps for template:', nramps
           !   open(19, file="templcorr2.dat", recl=1024)
           !   do l = 1, nsamp
           !      write(19,*) l, buffer(l), total(l), buffer(l)-total(l)
           !   end do
           !   close(19)
           !end if
           
           deallocate(buffer, total, mask)
        end do
        
        !if (self%myid==0) then
        !   open(19, file="templ2.dat", recl=1024)
        !   do l = 1, self%nsamp_templ(1)
        !      write(19,*) l, self%tod_correction_templ(l,1,m,scan)
        !   end do
        !   close(19)
        !end if

     end do
     deallocate(res)
     
   end subroutine sample_ramp

   module subroutine init_sample_ramp(self, sd)
     ! scans tod flags to find length of tod segment between ramp flags (= length of correction template)
     ! puts this into self%nsamp_templ
     ! used to initialize correction template structure
     implicit none
     class(comm_akari_tod), intent(inout)          :: self  ! full data object with compressed tods
     class(comm_scandata),  intent(in)             :: sd    ! decompressed self%scans(scan)%d(det)%tod

     integer(i4b) :: i, j, k, l, nbin, b, ndet, scan, ntod, num_used, nsamp, nflags
     integer(i4b) :: num_events(3), flag_id(3), max_num, maxlength, num
     real(dp)     :: dt, invgain, flag(3), flag_val(3), length
     logical(lgt) :: flagged(3)
     !real(dp),       allocatable, dimension(:)     :: res
     integer(i4b),   allocatable, dimension(:,:)   :: nflagged, ndata, events

     ! flag bit 5 = ramp reset, 15 = cal lamp on for long wavelengts, 16 = cal lamp for short wavelength
     flag_id(1) =  5
     flag_id(2) = 15
     flag_id(3) = 16
     ! output value for plotting
     flag_val(1) = 20.d0
     flag_val(2) = 30.d0
     flag_val(3) = 40.d0

     nflags = 1                    ! OBS set here
     scan   = sd%scan              ! scan number
     ndet   = sd%ndet              ! number of detectors                = self%ndet
     ntod   = sd%ntod              ! number of tod samples for my scan  = self%scans(scan)%ntod
     dt     = 1.d0/self%samprate   ! time per/between sample in seconds
     !ntempl = self%ntempl          ! number of tod correction templates we will generate

     
!     if (self%myid==0) write (*,*) 'ntod', ntod
!     if (self%myid==0) write (*,*) 'ndet', ndet
!     if (self%myid==0) write (*,*) 'dt  ', dt
!     if (self%myid==0) write (*,*) 'length of tod (s)', dt*ntod
!     if (self%myid==0) write (*,*) 'length of ramp (s)', dt*32
!     if (self%myid==0) write (*,*) 'samprate', self%samprate
!     if (self%myid==0) write (*,*) ''

     max_num     = 10000 ! OBS
     maxlength   =   100 ! max length of ramp segment  OBS hack until we get new hdf files with all flags
     allocate(nflagged(max_num, nflags))
     allocate(ndata(0:max_num, nflags))
     allocate(events(max_num, nflags))
!     allocate(res(1:ntod))
     
     do j = 1, ndet
        ! if (self%myid==0) write (*,*) 'det scan samprate ', trim(self%label(j)), self%scanid(scan), self%samprate
        if (.not. self%scans(scan)%d(j)%accept) cycle
        ! get tod data for residual = full tod - signal_tod ! OBS ignoring ncorr which is not set yet
        !res = sd%tod(:,j) - sd%s_tot(:,j,0,1) * self%scans(scan)%d(j)%gain

        ! loop through samples and check for flags
        flagged(:)    = .false.
        num_events(:) = 0
        nflagged      = 0
        ndata         = 0
!        if (self%myid==0) open(19, file="res.dat", recl=1024)
        do i = 1, ntod
           flag(:)    = 0.d0      
           do k = 1, nflags
              if ( btest(sd%flag(i,j), flag_id(k) ) ) then
!                 flag(k) = flag_val(k)
                 if (.not. flagged(k)) num_events(k) = num_events(k) + 1
                 flagged(k) = .true.
                 nflagged(num_events(k),k) = nflagged(num_events(k),k) + 1
                 events(num_events(k),k) = i
              else
                 flagged(k) = .false.
                 ndata(num_events(k),k) = ndata(num_events(k),k) + 1
              end if
           end do
!           if (self%myid==0) write (19,*) i, res(i), flag(1:nflags)
        end do
!        if (self%myid==0) close(19)
        
!        if (self%myid==0) open(19, file="ndata.dat", recl=1024)
        do k = 1, nflags
           do i = 1, num_events(k)
!              if (self%myid==0) write(19,*) ndata(i,k)
              if (ndata(i,k) > maxlength) ndata(i,k) = 0
           end do
           length = maxval(ndata(1:num_events(k),k))
           do l = length, 1, -1
              num = 0
              do i = 1, num_events(k)
                 if (ndata(i,k) == length) num = num + 1
              end do
              if (num > 1) exit
           end do
!           self%nsamp_templ(1, j, scan) = l
           self%nsamp_templ(1) = l
        end do
 !       if (self%myid==0) close(19)
        
     end do

!     if (self%myid==0) then
!        open(19, file="length_det.dat", recl=1024)
!        do j = 1, ndet
!           write(19,*) j, self%nsamp_templ(1, j, scan)
!        end do
!        close(19)
!     end if
     
     deallocate(ndata, nflagged, events)!, res)
     
   end subroutine init_sample_ramp
   
end submodule comm_tod_akari_smod
