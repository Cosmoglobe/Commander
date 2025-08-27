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
module comm_tod_DIRBE_mod
   !   Module which contains all the LiteBIRD time ordered data processing and routines
   !   for a given frequency band
   !
   !   Main Methods
   !   ------------
   !   constructor(cpar, id_abs, info, tod_type)
   !       Initialization routine that reads in, allocates and associates 
   !       all data needed for TOD processing
   !   process_DIRBE_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
   !       Routine which processes the time ordered data
   !
  use comm_tod_driver_mod
  use comm_tod_pixhist_mod
  use comm_tod_mapmaking_mod
   implicit none

   private
   public comm_dirbe_tod

   type, extends(comm_tod) :: comm_dirbe_tod
   contains
      procedure     :: process_tod          => process_DIRBE_tod
      procedure     :: dumpToHDF_inst       => dumpToHDF_DIRBE
   end type comm_dirbe_tod

   interface comm_dirbe_tod
      procedure constructor_dirbe
   end interface comm_dirbe_tod

contains
   !**************************************************
   !             Constructor
   !**************************************************
   function constructor_dirbe(cpar, id, id_abs, info, tod_type) result(c)
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
      class(comm_dirbe_tod),      pointer    :: c

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
      read(c%freq(1:2),*) c%zodiband
      c%samprate_lowres = 8.  ! Lowres samprate in Hz
      c%nhorn           = 1
      c%ndiode          = 1
      c%compressed_tod  = .false.
      c%correct_sl      = .false.
      c%correct_orb     = .false.
      c%orb_4pi_beam    = .false.
      c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
      c%use_moon_point  = cpar%sample_moon_maps
      c%use_earth_elon  = cpar%sample_earth_maps
      c%symm_flags      = .false.
      ! c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      c%chisq_threshold = 50000.
      c%nmaps           = info%nmaps
      c%ndet            = num_tokens(trim(cpar%ds_tod_dets(id_abs)), ",")
      c%sol_elong_range = cpar%zs_sol_elong
      c%nside_pixhist   = 64
      
      c%gain_tune_sigma0 = .false.
      c%gain_samprate    = 1.d0 / 3600.d0 / 24.d0
      c%gain_sigma_0     = 0.005
      
      nside_beam                  = 128
      nmaps_beam                  = 1
      pol_beam                    = .false.
      c%nside_beam      = nside_beam

      ! Get detector labels
      call get_tokens(trim(cpar%ds_tod_dets(id_abs)), ",", c%label)

      ! Read the actual TOD
      call c%read_tod(c%label)

      ! Initialize bandpass mean and proposal matrix
      call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
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
      !do i = 1, c%nscan
      !   c%scans(i)%d%baseline = 0.d0
      !end do

      call timer%stop(TOD_INIT, id_abs)
    end function constructor_dirbe

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_DIRBE_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      ! 
      ! Routine that processes the DIRBE Calibrated Individual Observations. 
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
      ! Writes maps to disc in fits format
      ! 
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_dirbe_tod class
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
      class(comm_dirbe_tod),                    intent(inout) :: self
      character(len=*),                         intent(in)    :: chaindir
      integer(i4b),                             intent(in)    :: chain, iter
      type(planck_rng),                         intent(inout) :: handle
      type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
      real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map),                          intent(inout) :: map_out      ! Combined output map
      class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr),       dimension(1:),       intent(inout), optional :: map_gain       ! (ndet)

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
      real(sp),              dimension(9)       :: flag_threshold
      real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi
      real(sp), allocatable, dimension(:,:,:)   :: d_calib
      real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
      real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf, freqmap
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
      sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                         ! don't sample absolute bandpasses
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
      !sample_gain           = .true. !iter > 1                         ! Gain sampling, LB TOD sims have perfect gain
!!$      if (trim(self%freq) == '01' .or. trim(self%freq) == '02' .or. &
!!$        & trim(self%freq) == '03' .or. &
!!$        & trim(self%freq) == '09' .or. trim(self%freq) == '10') then
      !if (trim(self%freq(1:2)) == '09' .or. trim(self%freq(1:2)) == '10') then
      if (trim(self%freq(1:2)) == '09' .or. trim(self%freq(1:2)) == '10') then
         sample_ncorr = .true.
         !sample_ncorr = .false.
      else
         sample_ncorr = .false.
      end if

      if (trim(self%freq(1:2)) == '05' .or. trim(self%freq(1:2)) == '06') then
         sample_gain        =  iter > 1
      else
         sample_gain = .false.
      end if
      
      if (trim(self%freq(1:2)) == '05' .or. trim(self%freq(1:2)) == '06' .or. &
        & trim(self%freq(1:2)) == '07' .or. trim(self%freq(1:2)) == '08' .or. &
        & trim(self%freq(1:2)) == '09' .or. trim(self%freq(1:2)) == '10') then
         !                     Pixhist  Extreme           RMS ranges     Single     Ranges   Pointing
         flag_threshold     = [1.0,    -20., -5.0,       -3.0, -2.0, -1.5,     1.0,      -0.30,     1.0]
      else
         !                     Pixhist    Extreme            RMS ranges     Single     Ranges   Pointing
         flag_threshold     = [1.0,    -50., -1.0,      -1.0, -1.0, -1.0,   1.0,        -1.0,      1.0]
         !flag_threshold     = [-1.0, -1.0,      -1.0, -1.0, -1.0,   1.0,        -1.0,      1.0]
      end if

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
      if (sample_ncorr) then
         oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,&
              & SD_TOD,SD_SKY,SD_SL,SD_ORB,SD_INST,SD_ZODI, SD_NCORR])
      else
         oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,&
              & SD_TOD,SD_SKY,SD_SL,SD_ORB,SD_INST,SD_ZODI])
      end if

      ! Initialize index-based sky map and mask
      call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain)

      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      ! Create pixel histograms
      if (self%first_call) call compute_tod_pixhist(self)
      
      ! Sample gain components in separate TOD loops; marginal with respect to n_corr
      if (sample_gain) then
         ! 'abscal': the global constant gain factor
         !call sample_calibration(self, 'abscal', handle, map_sky, m_gain, procmask, procmask2)
         ! 'relcal': the gain factor that is constant in time but varying between detectors
         ! call sample_calibration(self, 'relcal', handle, map_sky, m_gain, procmask, procmask2)
         ! 'deltaG': the time-variable and detector-variable gain
         call sample_calibration(self, 'deltaG', oper_default, handle, smooth=.true.)
      end if

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
         call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)

         ! Create dynamic mask
         if (self%first_call) then
            do j = 1, sd%ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               call self%create_dynamic_mask(sd, j, flag_threshold)
            end do
            call dealloc_scan_data(sd)
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
         end if

         ! Sample correlated noise
         if (sample_ncorr) then
            call sample_n_corr(self, sd, handle, nomono=.true.)
            call sample_noise_psd(self, sd, handle)
         else
            call sample_noise_psd(self, sd, handle, only_sigma0=.true.)
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
         call compute_calibrated_data(self, i, sd, d_calib)    

         ! For debugging: write TOD to hdf
         if (.false.) then
            ! scan id appears to be the worst chi2
            if (self%scanid(i) < 10000) then 
               !print *, self%scanid(i)
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call write_hdf(tod_file, '/pix_sol', self%scans(i)%d(1)%pix_sol)
               call write_hdf(tod_file, '/todz', d_calib(1, :, :))
               call write_hdf(tod_file, '/s_sky', sd%s_sky)
               !call write_hdf(tod_file, '/n_corr', sd%n_corr)
               !call write_hdf(tod_file, '/s_sl', sd%s_sl)
               !call write_hdf(tod_file, '/s_orb', sd%s_orb)
               call write_hdf(tod_file, '/res', d_calib(2, :, :))
               call write_hdf(tod_file, '/zodi', d_calib(7, :, :))
               call write_hdf(tod_file, '/mask', sd%mask)
               call write_hdf(tod_file, '/flag', sd%flag)
               call write_hdf(tod_file, '/sigma0', self%scans(i)%d(1)%N_psd%sigma0)
               !do k = 1, size(sd%s_zodi_therm, dim=2)
               !   call int2string(k, scantext)
               !   call write_hdf(tod_file , '/zodi'//scantext, d_calib(8 + k, :, :))
               !end do
               call close_hdf_file(tod_file)
            end if
         end if

         ! Bin TOD
         call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

!!$         do j = 1, binmap%nobs
!!$            if (binmap%A_map(1,j) > 0.) write(*,*) j, real(binmap%b_map(1,1,j),sp), real(binmap%A_map(1,j),sp), real(binmap%b_map(1,1,j)/binmap%A_map(1,j),sp)
!!$         end do

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

      ! Synchronize and output flagging statistics in first iteration
      if (self%first_call) call self%report_dynamic_mask_stats
      
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
      deallocate(map_sky, procmask, procmask2)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

      call update_status(status, "tod_end"//ctext)
      
      call timer%stop(TOD_TOT, self%band)
   end subroutine process_DIRBE_tod   

       subroutine dumpToHDF_DIRBE(self, chainfile, path)
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
      class(comm_dirbe_tod),               intent(in)     :: self
      type(hdf_file),                      intent(in)     :: chainfile
      character(len=*),                    intent(in)     :: path


      if (self%myid == 0) then
         if (self%map_solar_allocated == .true.) then
           call write_hdf(chainfile, trim(adjustl(path))//'map_solar',  self%map_solar)
         end if
      end if

    end subroutine dumpToHDF_DIRBE



end module comm_tod_DIRBE_mod
