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
      class(comm_akari_tod),      pointer    :: c

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
      if (trim(c%level) == 'L1') then
          c%compressed_tod  = .true.
      else
          c%compressed_tod  = .false.
      end if
      c%correct_sl      = .false.
      c%correct_orb     = .false.
      c%orb_4pi_beam    = .false.
      c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
      c%symm_flags      = .false.
      c%use_moon_point  = cpar%sample_moon_maps
      c%use_earth_elon  = cpar%sample_earth_maps
      ! c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      ! c%chisq_threshold = 50000.
      c%chisq_threshold = huge(0.0d0)
      c%nside_pixhist   = 64

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
      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist, sample_zodi, use_k98_samp_groups, output_zodi_comps, sample_ncorr, only_solar_mask, sample_xi_n
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
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
      sample_gain           = iter > 1                        ! Gain sampling, LB TOD sims have perfect gain
      only_solar_mask       = .false.                        ! Only apply solar mask
      flag_threshold        = [1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
!!$      if (trim(self%freq) == '01' .or. trim(self%freq) == '02' .or. &
!!$        & trim(self%freq) == '03' .or. &
!!$        & trim(self%freq) == '09' .or. trim(self%freq) == '10') then
      !if (trim(self%freq(1:2)) == '09' .or. trim(self%freq(1:2)) == '10') then
      ! if (trim(self%freq(1:2)) == '10') then
      !    sample_ncorr = .true.
      ! else
      !    sample_ncorr = .false.
      ! end if
      sample_ncorr = .false. !iter > 2
      sample_xi_n  = .false.

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
         !call sample_calibration(self, 'deltaG', handle, map_sky, m_gain, procmask, procmask2)
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
         !write(*,*) i, self%scans(i)%d%accept
         if (.not. any(self%scans(i)%d%accept)) cycle
         call wall_time(t1)
         call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)

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
         
         ! Create dynamic mask
         if (self%first_call) then
            do j = 1, sd%ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               if (self%scans(i)%d(j)%N_psd%sigma0 .eq. 0.d0) write(*,*) 'debug sigma0 = 0.0'
               call self%create_dynamic_mask(sd, j, flag_threshold)
            end do
            call dealloc_scan_data(sd)
            if (.not. any(self%scans(i)%d%accept)) cycle
            call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
         end if

         ! Sample correlated noise
         if (sample_ncorr) then
            call sample_n_corr(self, sd, handle)
            if (sample_xi_n) then
               call sample_noise_psd(self, sd, handle)
            else
               call sample_noise_psd(self, sd, handle, only_sigma0=.true.)
            end if
         else
            call sample_noise_psd(self, sd, handle, only_sigma0=.true.)
         end if

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_tod_chisq(sd, j)
         end do

         ! Select data
!!$         if (select_data) call remove_bad_data(self, i, sd%flag)

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
         
         ! Compute binned map
         allocate(d_calib(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib = 0.d0
         call compute_calibrated_data(self, i, sd, d_calib)    

         !write(*,*) "Scan = ", self%scanid(i), ', num moon = ', count(iand(sd%flag,2)==2)
         
         ! For debugging: write TOD to hdf
         if (.false.) then
            ! scan id appears to be the worst chi2
            write(*,*) "DEBUG Writing scan ", self%scanid(i)
            if (self%scanid(i) > 1075 .and. self%scanid(i) < 1085) then 
               !print *, self%scanid(i)
               call int2string(self%scanid(i), scantext)
               call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//scantext//'.h5', tod_file, 'w')
               call write_hdf(tod_file, '/tod', sd%tod)
               call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
               call write_hdf(tod_file, '/flag', sd%flag)
               call write_hdf(tod_file, '/todz', d_calib(1, :, :))
               call write_hdf(tod_file, '/s_sky', sd%s_sky)
               call write_hdf(tod_file, '/n_corr', sd%n_corr)
               call write_hdf(tod_file, '/s_sl', sd%s_sl)
               call write_hdf(tod_file, '/s_orb', sd%s_orb)
               call write_hdf(tod_file, '/res', d_calib(2, :, :))
               call write_hdf(tod_file, '/zodi', d_calib(7, :, :))
               call write_hdf(tod_file, '/mask', sd%mask)
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

      ! Synchronize and output flagging statistics in first iteration
      if (self%first_call) call self%report_dynamic_mask_stats
      
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

   

end submodule comm_tod_akari_smod
