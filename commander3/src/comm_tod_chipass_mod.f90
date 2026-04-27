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
   implicit none

   private
   public comm_chipass_tod

   type, extends(comm_tod) :: comm_chipass_tod
   contains
      procedure     :: process_tod          => process_chipass_tod
      procedure     :: construct_corrtemp_inst => construct_corrtemp_chipass
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

      call timer%start(TOD_INIT, id_abs)

      ! Allocate object
      allocate(c)

      ! Set up noise PSD type and priors
      c%freq            = cpar%ds_label(id_abs)
      c%n_xi            = 3
      c%noise_psd_model = 'oof'
      allocate(c%xi_n_nu_fit(c%n_xi,2))
      allocate(c%xi_n_P_uni(c%n_xi,2))
      allocate(c%xi_n_P_rms(c%n_xi))

      c%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0] 
      ! [sigma0, fknee, alpha]; sigma0 is not used
      do i = 1, c%n_xi 
         c%xi_n_nu_fit(i,:) = [0.d0, 0.01d0] 
      end do
      c%xi_n_P_uni(1,:) = [0.d0, 0.d0]
      c%xi_n_P_uni(2,:) = [0.00001d0, 0.3d0]  ! fknee
      c%xi_n_P_uni(3,:) = [-3.0d0, -0.4d0]   ! alpha

      ! Initialize common parameters
      call c%tod_constructor(cpar, id, id_abs, info, tod_type)

      ! Initialize instrument-specific parameters
      !read(c%freq(1:2),*) c%zodiband
      c%samprate_lowres = 8.  ! Lowres samprate in Hz
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
      ! c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      c%chisq_threshold = 50000.
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

      ! Allocate sidelobe convolution data structures
      allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)

      c%orb_dp => comm_orbdipole(comm=c%comm)


      ! Initialize all baseline corrections to zero
      do i = 1, c%nscan
         do j = 1, c%ndet
            c%scans(i)%d(j)%baseline = 0.d0
         end do 
      end do
      
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
      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist, sample_zodi, use_k98_samp_groups, output_zodi_comps, sample_ncorr
      type(comm_binmap)   :: binmap
      type(comm_scandata) :: sd
      character(len=4)    :: ctext, myid_text
      character(len=2)    :: zodi_param_text
      character(len=1)    :: up_down_text
      character(len=6)    :: samptext, scantext
      character(len=512)  :: prefix, postfix, prefix4D, prefix_atlas, postfix_atlas
      character(len=512), allocatable, dimension(:) :: slist
      !real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi
      real(sp), allocatable, dimension(:,:,:)   :: d_calib
      !real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
      real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)    :: A_T_A, A_T_A_reduced
      real(dp), allocatable, dimension(:)       :: AY, AY_reduced, X
      real(dp), allocatable, dimension(:, :, :) :: s_therm_tot, s_scat_tot ! (n_tod_tot, ncomps, ndet)
      real(dp), allocatable, dimension(:, :)    :: res_tot ! (n_tod_tot, ndet)
      real(dp), allocatable, dimension(:)       :: mask_tot ! (n_tod_tot)
      type(hdf_file) :: tod_file

      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)

      !call timer%start(TOD_ALLOC, self%band)

      
      ! Toggle optional operations
      sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
      output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
      use_k98_samp_groups   = .false.                          ! fits one overall albedo and episolon for the dust bands, and one for ring + feature
      !if (self%myid == 0) write(*, *) '[comm_tod_chipass_mod.f90] size(delta,3) > 1:', size(delta,3) > 1, size(delta,1), size(delta,2), size(delta,3)
      sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                         ! don't sample absolute bandpasses
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
      sample_gain           = .true.                         ! Gain sampling
      sample_ncorr = .true.
         
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

      ! Define useful sd operation codes
      oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
           & SD_SKY,SD_BP,SD_INST])

      ! Initialize index-based sky map and mask
      ! NOTE: CHIPASS TOD given in Jy beam^-1 but parameter file assumes MJy/sr
            ! for a 14.3 arcmin beam the conversion factor is roughly 0.057793   
      call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
      call timer%stop(TOD_ALLOC, self%band)
      call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
      call update_status(status, "tod_init")

      ! Write mask for debugging
!!$      if (.false. .and. self%myid == 0) then
!!$         print *, "writing masks"
!!$         call open_hdf_file(trim(chaindir)//'/mask.h5', tod_file, 'w')
!!$         call write_hdf(tod_file, '/procmask', procmask)
!!$         call write_hdf(tod_file, '/procmask2', procmask2)
!!$         call write_hdf(tod_file, '/procmask_zodi', procmask_zodi)
!!$         call close_hdf_file(tod_file)
!!$         stop
!!$      end if


      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      ! sample baseline
      ! Sample baseline for current scan
      if (.true.) then
         if (self%myid == 0) then
            write(*,*) '|    --> Sampling baseline'
         end if
         self%apply_inst_corr = .false. ! Disable baseline correction for just this call
         call update_status(status, "baseline")
         do i = 1, self%nscan
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
            ! check there are valid points
            l = 0
            do j = 1, self%ndet
               do k = 1, self%scans(i)%ntod
                  if (sd%mask(k,j) > 0.5) then
                     l    = l + 1
                  end if
               end do
            end do
            if (.not. l > self%baseline_order+1) cycle
            !call sd%init_differential(self, i, map_sky, m_gain, procmask, procmask2)
            call timer%start(TOD_BASELINE, self%band)
            !if (i == 1 .and. self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] 1 sample_chipass_baseline'
            !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_baseline; self%apply_inst_corr:', self%apply_inst_corr
            !if (self%myid == 0) write(*,*) '    size(sd%mask):', size(sd%mask)
            !if (self%myid == 0) write(*,*) '    sd%mask(1,1):', sd%mask(1,1)
            call sample_chipass_baseline(self, i, sd%tod, sd%s_tot(:,:,0,1), sd%mask, handle) ! changed ind from 1 to 0
            !if (i == 1 .and. self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] 2 sample_chipass_baseline done'
            !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_baseline done'
            call timer%stop(TOD_BASELINE, self%band)
            ! apply inst corr
            call dealloc_scan_data(sd)
         end do
         self%apply_inst_corr = .true.

         call timer%start(TOD_WAIT, self%band)
         call mpi_barrier(self%comm, ierr)
         call timer%stop(TOD_WAIT, self%band)
      end if

      ! apply inst corr
      !if (self%apply_inst_corr) then
      !   do i = 1, self%nscan
      !      call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
      !      call construct_corrtemp_chipass(self, sd) ! added
      !      call timer%start(TOD_INSTCORR, self%band)
      !      do j = 1, self%ndet
      !         if (.not. self%scans(i)%d(j)%accept) cycle
      !         if (self%myid == 0) write(*, *) 'j:', j
      !         if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%s_inst(:,j):', sd%s_inst(:,j)
      !         if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%tod(:,j):', sd%tod(:,j)
      !         sd%tod(:,j) = sd%tod(:,j) - sd%s_inst(:,j)
      !         if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%tod(:,j) new:', sd%tod(:,j)
      !      end do
      !      call dealloc_scan_data(sd)
      !   end do
      !   call timer%stop(TOD_INSTCORR, self%band)
      !end if

      ! Sample gain components in separate TOD loops; marginal with respect to n_corr
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_gain:', sample_gain
      if (sample_gain) then
         ! 'abscal': the global constant gain factor
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] handle:', handle
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sampling calibration'
         call sample_calibration(self, 'abscal', oper_default, handle)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sampling calibration done'
         ! 'relcal': the gain factor that is constant in time but varying between detectors
         call sample_calibration(self, 'relcal', oper_default, handle)         
         ! 'deltaG': the time-variable and detector-variable gain
         !call sample_calibration(self, 'deltaG', handle, map_sky, m_gain, procmask, procmask2)
         !call sample_calibration(self, 'deltaG', oper_default, handle)
      end if

      ! Prepare intermediate data structures
      call binmap%init(self, .true., sample_rel_bandpass)
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] samp. abs. or rel. bp?', (sample_abs_bandpass .or. sample_rel_bandpass)
      if (sample_abs_bandpass .or. sample_rel_bandpass) then
         allocate(chisq_S(self%ndet,size(delta,3)))
         chisq_S = 0.d0
      end if
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] output_scanlist?', output_scanlist
      if (output_scanlist) then
         allocate(slist(self%nscan))
         slist   = ''
      end if

      ! Perform loop over scans
      if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] self%nscan:', self%nscan
      do i = 1, self%nscan

         ! Skip scan if no accepted data
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] skip?', (.not. any(self%scans(i)%d%accept))
         if (.not. any(self%scans(i)%d%accept)) cycle
         call wall_time(t1)
         call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
         allocate(sd%n_corr(sd%ntod, sd%ndet))
         sd%n_corr = 0d0
         !if (self%myid == 0) write(*,*) '   [comm_tod_chipass_mod] size(sd%n_corr, 1), size(sd%n_corr, 2):', size(sd%n_corr, 1), size(sd%n_corr, 2)

!!$         if (.false. .and. self%myid == 0 .and. i == 1) then
!!$            open(58,file='tod.dat', recl=1024)
!!$            do j = 1, sd%ntod
!!$               write(58,*) j, sd%tod(j,1), sd%mask(j,1), sd%flag(j,1), sd%n_corr(j,1), sd%s_tot(j,1), sd%s_sky(j,1), sd%s_bp(j,1)
!!$            end do
!!$            close(58)
!!$         end if

         ! apply inst corr
         if (self%apply_inst_corr) then
            call construct_corrtemp_chipass(self, sd) ! added
            call timer%start(TOD_INSTCORR, self%band)
            do j = 1, self%ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               !if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] j:', j
               !if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%s_inst(:,j):', sd%s_inst(:,j)
               !if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%tod(:,j):', sd%tod(:,j)
               sd%tod(:,j) = sd%tod(:,j) - sd%s_inst(:,j)
               !if (self%myid == 0) write(*, *) '    [comm_tod_chipass_mod] sd%tod(:,j) new:', sd%tod(:,j)
            end do
            call timer%stop(TOD_INSTCORR, self%band)
         end if

         ! Sample correlated noise
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_ncorr?', sample_ncorr
         if (sample_ncorr) then
            call sample_n_corr(self, sd, handle)
            call sample_noise_psd(self, sd, handle)
         else
            sd%n_corr = 0.d0
            call sample_noise_psd(self, sd, handle, only_sigma0=.true.)
         end if
         

         ! Compute chisquare for bandpass fit
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] compute_chisq_abs_bp'
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] compute_chisq_abs_bp done'

         ! Compute binned map
         allocate(d_calib(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib = 0.d0
         !write(*,*) 'a', self%scanid(i), any(sd%s_zodi_scat/=sd%s_zodi_scat), any(sd%s_zodi_therm/=sd%s_zodi_therm)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] compute_calibrated_data'
         call compute_calibrated_data(self, i, sd, d_calib)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] compute_calibrated_data done'

         ! For debugging: write TOD to hdf
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] self%scanid(i):', self%scanid(i)
         if (.true.) then
            ! scan id appears to be the worst chi2
            !if (self%scanid(i) < 10000) then
            !if (self%scanid(i) < 10) then
            if (mod(self%scanid(i), 1690) == 0) then
               !print *, self%scanid(i)
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call write_hdf(tod_file, '/todz', d_calib(1, :, :))
               call write_hdf(tod_file, '/s_sky', sd%s_sky)
               call write_hdf(tod_file, '/n_corr', sd%n_corr)
               call write_hdf(tod_file, '/s_sl', sd%s_sl)
               call write_hdf(tod_file, '/s_orb', sd%s_orb)
               call write_hdf(tod_file, '/res', d_calib(2, :, :))
               !call write_hdf(tod_file, '/zodi', d_calib(7, :, :))
               call write_hdf(tod_file, '/mask', sd%mask)
               call write_hdf(tod_file, '/sigma0', self%scans(i)%d(1)%N_psd%sigma0)
!!$               do k = 1, size(sd%s_zodi_therm, dim=2)
!!$                  call int2string(k, scantext)
!!$                  call write_hdf(tod_file , '/zodi'//scantext, d_calib(8 + k, :, :))
!!$               end do
               call close_hdf_file(tod_file)
            end if
         end if

         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] ', trim(chaindir)//'/map_in.fits'
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] dim pix: ', size(sd%pix, 1), size(sd%pix, 2), size(sd%pix, 3)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] dim tod: ', size(sd%tod, 1), size(sd%tod, 2)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] self%label(1):', trim(self%label(1))
         if (.false.) then
            ! scan id appears to be the worst chi2
            if (self%scanid(i) == 10000) then
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/tod_point_scan_'//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call close_hdf_file(tod_file)
            end if
         end if

         ! Bin TOD
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] bin_TOD'
         call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] bin_TOD done'

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
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] synchronize_binmap'
      call synchronize_binmap(binmap, self)
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] synchronize_binmap done'
!!$      if (sample_rel_bandpass) then
!!$         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_rel_bandpass = T; self%nmaps:', self%nmaps
!!$         if (self%nmaps > 1) then
!!$            call finalize_binned_map(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2)
!!$         else
!!$            call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0, chisq_S=chisq_s, mask=procmask2)
!!$         end if
!!$      else
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_rel_bandpass = F; self%nmaps:', self%nmaps
         if(self%nmaps > 1) then
            call finalize_binned_map(self, binmap, rms_out, 1.d0)
         else 
            call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0)
         end if
!!$      end if
      map_out%map = binmap%outmaps(1)%p%map

      ! Sample bandpass parameters
      !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample bp?', (sample_rel_bandpass .or. sample_abs_bandpass)
      if (sample_rel_bandpass .or. sample_abs_bandpass) then
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_bp'
         call sample_bp(self, handle, chisq_S, delta)
         !if (self%myid == 0) write(*,*) '[comm_tod_chipass_mod] sample_bp done'
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
         ! added
         if (self%myid == 0) write(*,*) '   [comm_tod_chipass_mod] output_n_maps:', self%output_n_maps
         if (self%output_n_maps > 7) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'baseline'//trim(postfix))
         if (self%output_n_maps > 6 .and. self%subtract_zodi) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
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

      ! Clean up
      call binmap%dealloc()
      if (allocated(slist)) deallocate(slist)
      !deallocate(map_sky, procmask, procmask2)
      !  if (self%correct_sl) then
      !     do i = 1, self%ndet
      !        call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
      !     end do
      !  end if

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

      call update_status(status, "tod_end"//ctext)
      
      call timer%stop(TOD_TOT, self%band)
   end subroutine process_chipass_tod


   subroutine sample_chipass_baseline(tod, scan, raw, s_tot, mask, handle)
    !   Sample LFI specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   scan:     local scan ID
    !   raw:      raw tod in du
    !   s_tot:    total signal model in mK
    !   mask:     list of accepted samples
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    implicit none
    class(comm_chipass_tod),                   intent(inout) :: tod
    integer(i4b),                           intent(in)    :: scan
    real(sp),            dimension(1:,1:),  intent(in)    :: raw, s_tot, mask
    type(planck_rng),                       intent(inout) :: handle

    integer(i4b) :: i, j, k, n
    real(dp)     :: dt, t_tot, t, A, b, mval, eta
    real(dp), allocatable, dimension(:) :: x, y
    !real(dp), dimension(0:4) :: b2
    !real(dp), dimension(0:4, 0:4) :: C

    allocate(x(tod%scans(scan)%ntod), y(tod%scans(scan)%ntod))
    dt = 1.d0 / tod%scans(scan)%ntod

    !if (tod%myid == 0) write(*,*) ''
    !if (tod%myid == 0) write(*,*) '    [comm_tod_chipass_mod] scan', scan


    do j = 1, tod%ndet
       !if (tod%myid == 0) write(*,*) '    [comm_tod_chipass_mod]    det:', j
       !if (tod%myid == 0) write(*,*) '    [comm_tod_chipass_mod]    gain:', tod%scans(scan)%d(j)%gain
       if (.not. tod%scans(scan)%d(j)%accept) cycle

       t = 0.d0
       n = 0
       !if (tod%myid == 0) write(*,*) '    ktod ', 'raw ', 'gain*s_tot ', 'y'
       do k = 1, tod%scans(scan)%ntod
          t      = t + dt
          if (mask(k,j) > 0.5) then
             n    = n + 1
             x(n) = t
             y(n) = raw(k,j) - tod%scans(scan)%d(j)%gain * s_tot(k,j)
             !if (tod%myid == 0) write(*,*) k, raw(k,j), tod%scans(scan)%d(j)%gain * s_tot(k,j), y(n)
          end if
       end do

       if (n > tod%baseline_order+1) call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline)
       !if (tod%myid == 0 .and. mod(scan, 1690) == 0) write(*,*) '    [comm_tod_chipass_mod] n, ntod, baseline:', &
       !        & n, tod%scans(scan)%ntod, tod%scans(scan)%d(j)%baseline

       !do i = 0, 4
       !   b2(i) = sum(y(1:n) * x(1:n)**i)
       !   do k = i, 4
       !      C(i,k) = sum(x(1:n)**(i+k))
       !      C(k,i) = C(i,k)
       !   end do
       !end do
       !call solve_system_real(C, tod%scans(scan)%d(j)%baseline, b2)
       !write(*,*) '    sample_baseline fit_polynomial done'
    end do
    !if (tod%myid == 0) then
    !   write(*,*) '    [comm_tod_chipass_mod] scan', scan
    !   write(*,*) '        ntod:', tod%scans(scan)%ntod
    !   !write(*,*) '        x:', x(1:n)
    !   write(*,*) '        y:', y(1:n)
    !   write(*,*) '        raw:', raw(:, 13)
    !   write(*,*) '        s_tot:', s_tot(:, 13)
    !   write(*,*) '        mask:', mask(:,tod%ndet)
    !   write(*,*) '        det 13 gain:', tod%scans(scan)%d(13)%gain
    !   write(*,*) '        baseline det 13 out:', tod%scans(scan)%d(13)%baseline
    !end if

    deallocate(x, y)

  end subroutine sample_chipass_baseline

  subroutine construct_corrtemp_chipass(self, sd, det)
    !  Construct a CHIPASS instrument-specific correction template; for now contains baseline
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
    class(comm_chipass_tod), intent(in)             :: self
    class(comm_scandata), intent(inout)          :: sd
    integer(i4b),         intent(in),   optional :: det

    integer(i4b) :: i, j, k, d, nbin, b, scan
    real(dp)     :: dt, t

    !if (self%myid == 0) write(*,*) '    [comm_tod_chipass_mod] construct_corrtemp_chipass 0'
    scan      = sd%scan
    dt        = 1.d0 / self%scans(scan)%ntod
    sd%s_inst = 0.
    do j = 1, self%ndet
       d = j; if (present(det)) d = det
       t = 0.d0
       if (.not. self%scans(scan)%d(d)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          do i = 0, self%baseline_order
             sd%s_inst(k,j) = sd%s_inst(k,j) + self%scans(scan)%d(d)%baseline(i) * t**i
          end do
       end do
    end do
    !if (self%myid == 0) write(*,*) '    [comm_tod_chipass_mod] construct_corrtemp_chipass 1'

  end subroutine construct_corrtemp_chipass

end module comm_tod_chipass_mod
