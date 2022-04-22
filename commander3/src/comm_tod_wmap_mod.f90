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
module comm_tod_WMAP_mod
  !   Module which contains all the WMAP time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
   use comm_tod_mod
   use comm_param_mod
   use comm_map_mod
   use comm_conviqt_mod
   use pix_tools
   use healpix_types
   use comm_huffman_mod
   use comm_hdf_mod
   use comm_fft_mod
   use spline_1D_mod
   use comm_4D_map_mod
   use comm_zodi_mod
   use comm_tod_mapmaking_mod
   use comm_tod_pointing_mod
   use comm_tod_gain_mod
   use comm_tod_bandpass_mod
   use comm_tod_orbdipole_mod
   use comm_tod_driver_mod
   use comm_utils
   implicit none

   private
   public comm_WMAP_tod

   type, extends(comm_tod) :: comm_WMAP_tod
      integer(i4b) :: nside_M_lowres, nmaps_M_lowres
      logical(lgt) :: comp_S
      character(len=20), allocatable, dimension(:) :: labels ! names of fields
      real(dp), allocatable, dimension(:,:)        :: M_lowres, M_diag
   contains
      procedure     :: process_tod             => process_WMAP_tod
      procedure     :: precompute_M_lowres
      procedure     :: apply_map_precond       => apply_wmap_precond
      procedure     :: construct_corrtemp_inst => construct_corrtemp_wmap
   end type comm_WMAP_tod

   interface comm_WMAP_tod
      procedure constructor
   end interface comm_WMAP_tod

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
      !           Information about the maps for this band, e.g.,
      !           how the maps are distributed in memory
      ! tod_type: string
      !           Instrument specific tod type
      !
      ! Returns
      ! ----------
      ! constructor: pointer
      !              Pointer that contains all instrument data

      implicit none
      type(comm_params),      intent(in) :: cpar
      integer(i4b),           intent(in) :: id_abs
      class(comm_mapinfo),    target     :: info
      character(len=128),     intent(in) :: tod_type
      class(comm_WMAP_tod),   pointer    :: constructor

      integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam
      logical(lgt) :: pol_beam

      call timer%start(TOD_INIT, id_abs)

      ! Initialize common parameters
      allocate (constructor)

      ! Set up noise PSD type and priors
      constructor%freq            = cpar%ds_label(id_abs)
      !constructor%n_xi            = 3
      !constructor%noise_psd_model = 'oof'

      constructor%n_xi            = 5
      constructor%noise_psd_model = 'oof_f'
      constructor%comp_S          = .false.

      allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
      allocate(constructor%xi_n_nu_fit(constructor%n_xi,2))
      allocate(constructor%xi_n_P_rms(constructor%n_xi))
  
      ! Jarosik 2003 Table 2 gives knee frequencies between 0.09 mHz and 
      ! 46.5 mHz. 
      !constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2]   ! [sigma0, fknee, alpha]; sigma0 is not used
      constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2, -1.0, -1.0]   ! [sigma0, fknee, alpha, slope, intercept]; sigma0 is not used
      constructor%xi_n_P_uni(4,:) = [0.0, 0.1]            ! slope
      constructor%xi_n_nu_fit(4,:) = [0.1, 1.0]       ! slope nu_fit
      constructor%xi_n_P_uni(5,:) = [-1,1]             ! intercept
      constructor%xi_n_nu_fit(5,:) = [0.1, 1.0]       ! intercept nu_fit
      if (trim(constructor%freq) == '023-WMAP_K') then
         constructor%xi_n_nu_fit(2,:) = [0.0, 0.200]    
         constructor%xi_n_nu_fit(3,:) = [0.0, 0.200]    
         constructor%xi_n_P_uni(2,:) = [0.00001, 0.005]  ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '030-WMAP_Ka') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.200]    
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.200]    
         constructor%xi_n_P_uni(2,:) = [0.0001, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '040-WMAP_Q1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.200]    
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.200]    
         constructor%xi_n_P_uni(2,:) = [0.0001, 0.02]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '040-WMAP_Q2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.200]   
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.200]   
         constructor%xi_n_P_uni(2,:) = [0.0003, 0.02]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '060-WMAP_V1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.200]  
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.200]  
         constructor%xi_n_P_uni(2,:) = [0.0005, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '060-WMAP_V2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.200] 
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.200] 
         constructor%xi_n_P_uni(2,:) = [0.0005, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '090-WMAP_W1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 2.000]
         constructor%xi_n_nu_fit(3,:)     = [0.0, 2.000]
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.00]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '090-WMAP_W2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 2.000]
         constructor%xi_n_nu_fit(3,:)     = [0.0, 2.000]
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '090-WMAP_W3') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 2.000] 
         constructor%xi_n_nu_fit(3,:)     = [0.0, 2.000] 
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else if (trim(constructor%freq) == '090-WMAP_W4') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 2.000]  
         constructor%xi_n_nu_fit(3,:)     = [0.0, 2.000]  
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
      else
         write(*,*) 'Invalid WMAP frequency label = ', trim(constructor%freq)
         stop
      end if

      call constructor%tod_constructor(cpar, id_abs, info, tod_type)
      if (constructor%enable_tod_simulations) constructor%chisq_threshold = 1d6

      ! Set up WMAP specific parameters
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn           = 2
      constructor%ndiode          = 1
      constructor%baseline_order  = 1
      constructor%apply_inst_corr = .true.
      if (trim(constructor%level) == 'L1') then
          constructor%compressed_tod  = .true.
      else
          constructor%compressed_tod  = .false.
      end if
      constructor%correct_sl      = .true.
      constructor%orb_4pi_beam    = .true.
      constructor%symm_flags      = .false.
      constructor%chisq_threshold = 50.d0 ! 9.d0
      constructor%nmaps           = info%nmaps
      constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")
      constructor%verbosity       = cpar%verbosity
      constructor%x_im            = 0d0

      if (constructor%myid == 0) then
         allocate(constructor%M_diag(0:info%npix-1,info%nmaps+1))
      end if

      ! Iniitialize TOD labels
      allocate (constructor%labels(7))
      constructor%labels(1) = 'map'
      constructor%labels(2) = 'res'
      constructor%labels(3) = 'ncorr'
      constructor%labels(4) = 'bpcorr'
      constructor%labels(5) = 'orb'
      constructor%labels(6) = 'sl'
      constructor%labels(7) = 'zodi'

      ! Initialize beams
      nside_beam                  = 512
      nmaps_beam                  = 3
      pol_beam                    = .true.
      constructor%nside_beam      = nside_beam


      ! Get detector labels
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)

      ! Define detector partners
      ! I don't think this is necessary for WMAP...
      do i = 1, constructor%ndet
         if (mod(i,2) == 1) then
            constructor%partner(i) = i+1
         else
            constructor%partner(i) = i-1
         end if
         constructor%horn_id(i) = (i+1)/2
      end do

      ! Read the actual TOD
      call constructor%read_tod(constructor%label)

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(trim(cpar%datadir)//cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
      call constructor%precompute_lookups()

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

      ! Precompute low-resolution preconditioner
      call constructor%precompute_M_lowres

      ! Need precompute the main beam precomputation for both the A-horn and
      ! B-horn.
      ! Allocate sidelobe convolution data structures
      allocate(constructor%slconvA(constructor%ndet), constructor%slconvB(constructor%ndet))
      allocate(constructor%orb_dp)
      constructor%orb_dp => comm_orbdipole(constructor%mbeam)


      call timer%stop(TOD_INIT, id_abs)

   end function constructor

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      !
      ! Routine that processes the WMAP time ordered data.
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
      ! Writes maps to disc in fits format
      !
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_WMAP_tod class
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
      class(comm_WMAP_tod),             intent(inout) :: self
      character(len=*),                    intent(in) :: chaindir
      integer(i4b),                        intent(in) :: chain, iter
      type(planck_rng),                 intent(inout) :: handle
      type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in    ! (ndet,ndelta)
      real(dp),  dimension(0:, 1:, 1:), intent(inout) :: delta     ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map), intent(inout) :: map_out      ! Combined output map
      class(comm_map), intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr), dimension(1:, 1:), intent(inout), optional :: map_gain    ! (ndet,1)

      real(dp)     :: t1, t2
      integer(i4b) :: i, j, k, l, n
      integer(i4b) :: nside, npix, nmaps 
      integer(i4b) :: ierr, ndelta
      real(sp), allocatable, dimension(:, :)          :: s_buf
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: M_diag, buffer1
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono, buffer2
      character(len=512) :: prefix, postfix
      character(len=2048) :: Sfilename

      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, bp_corr, output_scanlist
      type(comm_scandata) :: sd

      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      real(sp),       allocatable, dimension(:)     :: procmask, procmask2, sigma0
      real(sp),  allocatable, dimension(:, :, :, :) :: map_sky
      class(map_ptr),     allocatable, dimension(:) :: outmaps

      ! biconjugate gradient-stab parameters
      integer(i4b) :: num_cg_iters
      real(dp) ::  epsil
      real(dp) ::  nullval
      real(dp), allocatable, dimension(:, :)    :: bicg_sol
      real(dp), allocatable, dimension(:)       :: map_full
      class(comm_map), pointer :: wmap_guess

      character(len=80), dimension(180) :: header


      type(hdf_file) :: tod_file



      ! Parameters used for testing
      real(dp) :: polang
      !polang = mod(2*PI*iter/12, 2*PI)
      polang = 0d0


!!$      do i = 1, self%ndet
!!$         call int2string(i, ctext)
!!$         call map_in(i,1)%p%writeFITS('map'//ctext//'.fits')
!!$      end do
!!$      call mpi_finalize(ierr)
!!$      stop


      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)

      ! Toggle optional operations
      sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
      bp_corr               = .true.                 ! by default, take into account differences in bandpasses. (WMAP does not do this in default analysis)
      bp_corr               = (bp_corr .or. sample_rel_bandpass) ! Bandpass is necessary to include if bandpass sampling is happening.
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration


      sample_rel_bandpass   = sample_rel_bandpass .and. .not. self%enable_tod_simulations
      sample_abs_bandpass   = sample_abs_bandpass .and. .not. self%enable_tod_simulations

      ! Initialize local variables
      ndelta          = size(delta,3)
      self%n_bp_prop  = ndelta-1
      nside           = map_out%info%nside
      nmaps           = map_out%info%nmaps
      npix            = 12*nside**2
      self%output_n_maps = 1
      if (self%output_aux_maps > 0) then
         if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 6
      end if

      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'

      ! Distribute maps
      ! Allocate total map (for monopole sampling)
      allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
      allocate(map_full(0:npix-1))
      map_full = 0.d0
      !call distribute_sky_maps(self, map_in, 1.e-3, map_sky) ! uK to mK
      call distribute_sky_maps(self, map_in, 1., map_sky, map_full) ! K to K?


      ! Distribute processing masks
      allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
      call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
      call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
      deallocate(m_buf)

      ! Precompute far sidelobe Conviqt structures
      if (self%correct_sl) then
         call timer%start(TOD_SL_PRE, self%band)
         do i = 1, self%ndet
            call map_in(i,1)%p%YtW()  ! Compute sky a_lms
            self%slconvA(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
                 & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
                 & 100, 3, 100, self%slbeam(1)%p, map_in(i,1)%p, 2)
            self%slconvB(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
                 & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
                 & 100, 3, 100, self%slbeam(3)%p, map_in(i,1)%p, 2)
                 ! lmax, nmaps, bmax, beam, map, optim
         end do
         call timer%stop(TOD_SL_PRE, self%band)
      end if
      ! In order to not completely break the rest of Commander, I am making the
      ! notation a little different for the beam indexing. Essentially, when it
      ! indexes 113, 114, 123, 124, the actual output will be A/A/B/B

      call update_status(status, "tod_init")

      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      ! Sample baseline for curren scan
      self%apply_inst_corr = .false. ! Disable baseline correction for just this call
      do i = 1, self%nscan
         if (.not. any(self%scans(i)%d%accept)) cycle
         call sd%init_differential(self, i, map_sky, procmask, procmask2, polang=polang)
         call sample_baseline(self, i, sd%tod, sd%s_tot, sd%mask, handle)
         call sd%dealloc
      end do
      self%apply_inst_corr = .true.

      ! Sample calibration
      if (.not. self%enable_tod_simulations) then
          if (trim(self%level) == 'L1') then
              call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2, polang)
              call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2, polang)
              call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2, polang)
          end if
          call sample_calibration(self, 'imbal',  handle, map_sky, procmask, procmask2, polang)
      end if


      ! Prepare intermediate data structures
      if (sample_abs_bandpass .or. sample_rel_bandpass) then
         allocate(chisq_S(self%ndet,size(delta,3)))
         chisq_S = 0.d0
      end if
      if (output_scanlist) then
         allocate(slist(self%nscan))
         slist   = ''
      end if

      allocate (M_diag(0:npix-1, nmaps+1))
      if (self%comp_S) then
         allocate ( b_map(0:npix-1, nmaps+1, self%output_n_maps))
      else
         allocate ( b_map(0:npix-1, nmaps,   self%output_n_maps))
      end if
      M_diag = 0d0
      b_map = 0d0

      ! Perform loop over scans
      if (self%myid == 0) then
          if (self%enable_tod_simulations) then
            write(*,*) '|    --> Simulating TODs'
          else
            write(*,*) '|    --> Sampling ncorr, xi_n, maps'
          end if
      endif
      do i = 1, self%nscan
         ! Skip scan if no accepted data

         if (.not. any(self%scans(i)%d%accept)) cycle
         call wall_time(t1)

         ! Prepare data
         if (sample_rel_bandpass) then
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=.true., init_s_bp_prop=.true., polang=polang)
         else if (sample_abs_bandpass) then
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=.true., init_s_sky_prop=.true., polang=polang)
         else
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=bp_corr, polang=polang)
         end if
         allocate(s_buf(sd%ntod,sd%ndet))

         ! Make simulations or Sample correlated noise
         if (self%enable_tod_simulations) then
            call simulate_tod(self, i, sd%s_tot, sd%n_corr, handle)
         else
            call timer%start(TOD_NCORR, self%band)
            call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,1,:), dospike=.false.)
            !sd%n_corr = 0.
            call timer%stop(TOD_NCORR, self%band)
         end if


         ! Compute noise spectrum parameters
         call timer%start(TOD_XI_N, self%band)
         call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)
         call timer%stop(TOD_XI_N, self%band)

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), &
              & sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
         end do

         ! Select data
         if (select_data) call remove_bad_data(self, i, sd%flag)

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

         ! Compute binned map
         allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
         call compute_calibrated_data(self, i, sd, d_calib)

!!$         do j = 1, self%ndet
!!$            write(*,*) i, j, minval(d_calib(4,:,j)), maxval(d_calib(4,:,j)), sqrt(variance(1.d0*d_calib(4,:,j)))
!!$         end do

         !if (mod(iter-1,self%output_aux_maps*10) == 0 .and. .not. self%enable_tod_simulations .and. iter .ne. 1) then
         if (.false.) then
            call int2string(self%scanid(i), scantext)
            if (self%myid == 0 .and. i == 1) write(*,*) '| Writing tod to hdf'
            call open_hdf_file(trim(chaindir)//'/tod_'//scantext//'_samp'//samptext//'.h5', tod_file, 'w')
            call write_hdf(tod_file, '/sl', sd%s_sl)
            call write_hdf(tod_file, '/s_orb', sd%s_orb)
            call write_hdf(tod_file, '/n_corr', sd%n_corr)
            call write_hdf(tod_file, '/bpcorr', sd%s_bp)
            call write_hdf(tod_file, '/s_tot', sd%s_tot)
            call write_hdf(tod_file, '/tod',   sd%tod)
            call write_hdf(tod_file, '/flag', sd%flag)
            call write_hdf(tod_file, '/pixA', sd%pix(:,1,1))
            call write_hdf(tod_file, '/pixB', sd%pix(:,1,2))
            call write_hdf(tod_file, '/psiA', sd%psi(:,1,1))
            call write_hdf(tod_file, '/psiB', sd%psi(:,1,2))
            call write_hdf(tod_file, '/x_im', self%x_im)

            do k = 1, self%ndet
              call int2string(k, scantext)
              call write_hdf(tod_file, '/xi_n_'//scantext, self%scans(i)%d(k)%N_psd%xi_n)
              call write_hdf(tod_file, '/gain_'//scantext, self%scans(i)%d(k)%gain)
            end do

            call close_hdf_file(tod_file)


         end if
         
         ! Bin TOD
         call timer%start(TOD_MAPBIN, self%band)
         call bin_differential_TOD(self, d_calib, sd%pix(:,1,:),  &
           & sd%psi(:,1,:), sd%flag(:,1), self%x_im, procmask, b_map, M_diag, i, &
           & self%comp_S)
         call timer%stop(TOD_MAPBIN, self%band)

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

      call mpi_barrier(self%comm, ierr)
      if (.not. self%enable_tod_simulations) then
        if (self%myid == 0) write(*,*) '|    --> Finalizing binned maps'

        ! Output latest scan list with new timing information
        if (output_scanlist) call self%output_scan_list(slist)


        call update_status(status, "Running allreduce on M_diag")
        call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
             & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
        call update_status(status, "Ran allreduce on M_diag")

        call update_status(status, "Running allreduce on b_map")
        call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
             & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
        call update_status(status, "Ran allreduce on b_map")


        where (M_diag == 0d0)
           M_diag = 1d0
        end where
        if (.not. self%comp_S) then
           ! If we want to not do the "better preconditioning"
           M_diag(:,4) = 0d0
        end if
        if (self%myid == 0) self%M_diag = M_diag

        call timer%start(TOD_MAPSOLVE, self%band)

        allocate(outmaps(1))
        outmaps(1)%p => comm_map(self%info)

        if (self%comp_S) then
          allocate (bicg_sol(0:npix-1, nmaps+1))
        else
          allocate (bicg_sol(0:npix-1, nmaps  ))
        end if


        ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
        do l = self%output_n_maps, 1, -1
          !if (l .ne. 6) b_map(:,:,l) = 0d0
          !b_map = 0d0
          bicg_sol = 0.0d0

          if (l == 1) then
            epsil = 1d-10
          else
            epsil = 1d-6
          end if
          num_cg_iters = 0

          ! Doing this now because it's still burning in...
          !if (mod(iter,self%output_aux_maps) == 0) then
            ! Solve for maps
            if (self%verbosity > 0 .and. self%myid == 0) then
              write(*,*) '|      Solving for ', trim(adjustl(self%labels(l)))
            end if
            call run_bicgstab(self, handle, bicg_sol, npix, nmaps, num_cg_iters, &
                           & epsil, procmask, map_full, M_diag, b_map, l, &
                           & prefix, postfix, self%comp_S)
          !end if

          call mpi_bcast(bicg_sol, size(bicg_sol),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
          call mpi_bcast(num_cg_iters, 1,  MPI_INTEGER, 0, self%info%comm, ierr)
          if (self%comp_S) then
             outmaps(1)%p%map(:,1) = bicg_sol(self%info%pix, nmaps+1)
             map_out%map = outmaps(1)%p%map
             call map_out%writeFITS(trim(prefix)//'S_'//trim(adjustl(self%labels(l)))//trim(postfix))
          end if
          do j = 1, nmaps
             outmaps(1)%p%map(:, j) = bicg_sol(self%info%pix, j)
          end do

          if (l == 1) then
             map_out%map = outmaps(1)%p%map
             rms_out%map = 1/sqrt(M_diag(self%info%pix, 1:nmaps))
             call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
             call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
          else
             call outmaps(1)%p%writeFITS(trim(prefix)//trim(adjustl(self%labels(l)))//trim(postfix))
          end if
        end do
        deallocate(bicg_sol)


        call timer%stop(TOD_MAPSOLVE, self%band)

        ! Sample bandpass parameters
        if (sample_rel_bandpass .or. sample_abs_bandpass) then
           call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
           self%bp_delta = delta(:,:,1)
        end if
      end if

      ! Clean up temporary arrays
      deallocate(procmask, procmask2)
      deallocate(b_map, M_diag)
      deallocate(map_full)
      if (allocated(chisq_S)) deallocate (chisq_S)
      if (allocated(b_mono)) deallocate (b_mono)
      if (allocated(sys_mono)) deallocate (sys_mono)
      if (allocated(slist)) deallocate (slist)

      if (allocated(outmaps)) then
         call outmaps(1)%p%dealloc
         deallocate (outmaps)
      end if

      deallocate(map_sky)

      if (self%correct_sl) then
         do i = 1, self%ndet
            call self%slconvA(i)%p%dealloc(); deallocate(self%slconvA(i)%p)
            call self%slconvB(i)%p%dealloc(); deallocate(self%slconvB(i)%p)
         end do
      end if

      call int2string(iter, ctext)
      call update_status(status, "tod_end"//ctext)
    call timer%stop(TOD_TOT, self%band)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

   end subroutine process_WMAP_tod


   subroutine precompute_M_lowres(self)
      !
      ! Routine that precomputes the low-resolution preconditioner, M_lowres = (P^t invN_w P)^{-1}
      !
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_WMAP_tod class
      !           Points to output of the constructor
      !
      implicit none
      class(comm_WMAP_tod),             intent(inout) :: self

      integer(i4b) :: i, j, k, t, p1, p2, k1, k2, ntot, npix, ierr, ntod, lpix, rpix, q, nhorn
      real(dp)     :: var, inv_sigma, lcos2psi, lsin2psi, rcos2psi, rsin2psi
      real(dp), allocatable, dimension(:)   :: wl, wr
      real(dp), allocatable, dimension(:,:) :: M
      integer(i4b), allocatable, dimension(:)         :: flag, dgrade
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi

      call update_status(status, "M_lowres")

      self%nmaps_M_lowres = 3; if (self%comp_S) self%nmaps_M_lowres = 4
      self%nside_M_lowres = 8
      npix                = 12  *self%nside_M_lowres**2
      ntot                = npix*self%nmaps_M_lowres
      nhorn               = self%nhorn

      ! Precompute udgrade lookup table
      allocate(dgrade(0:12*self%info%nside**2-1))
      q = (self%info%nside / self%nside_M_lowres)**2
      do i = 0, 12*self%info%nside**2-1
         call ring2nest(self%info%nside, i, j)
         j = j/q
         call ring2nest(self%nside_M_lowres, j, dgrade(i))
      end do

      ! Allocate local coupling matrix
      allocate(M(0:ntot-1,0:ntot-1), wl(self%nmaps_M_lowres), wr(self%nmaps_M_lowres))
      M = 0.d0

      ! Loop over scans
      do i = 1, self%nscan
         ! Skip scan if no accepted data
         if (.not. self%scans(i)%d(1)%accept) cycle

         ntod = self%scans(i)%ntod
         allocate(pix(ntod, nhorn))             ! Decompressed pointing
         allocate(psi(ntod, nhorn))             ! Decompressed pol angle
         allocate(flag(ntod))                   ! Decompressed flags
         call self%decompress_pointing_and_flags(i, 1, pix, psi, flag)

         var = 0.d0
         do k = 1, 4
            var = var + (self%scans(i)%d(k)%N_psd%sigma0/self%scans(i)%d(k)%gain)**2/4
         end do
         inv_sigma = sqrt(1.d0/var)

         do t = 1, ntod
            if (iand(flag(t),self%flag0) .ne. 0) cycle
            lpix = dgrade(pix(t, 1))
            rpix = dgrade(pix(t, 2))

            wl(1) = inv_sigma
            wl(2) = inv_sigma * self%cos2psi(psi(t,1))
            wl(3) = inv_sigma * self%sin2psi(psi(t,1))

            wr(1) = -inv_sigma
            wr(2) = -inv_sigma * self%cos2psi(psi(t,2))
            wr(3) = -inv_sigma * self%sin2psi(psi(t,2))

            do k1 = 1, self%nmaps_M_lowres
               p1 = (k1-1)*npix + lpix
               do k2 = 1, self%nmaps_M_lowres
                  p2 = (k2-1)*npix + rpix
                  M(p1,p1) = M(p1,p1) + wl(k1) * wl(k1)
                  M(p1,p2) = M(p1,p2) + wl(k1) * wr(k2)
                  M(p2,p1) = M(p2,p1) + wr(k2) * wl(k1)
                  M(p2,p2) = M(p2,p2) + wr(k2) * wr(k2)
               end do
            end do

         end do

         deallocate(pix, psi, flag)
      end do

      ! Collect contributions from all cores 
      if (self%myid == 0) then
         allocate(self%M_lowres(ntot,ntot))
         call mpi_reduce(M, self%M_lowres, size(M),  MPI_DOUBLE_PRECISION,  MPI_SUM,  0, self%comm, ierr)
         call invert_matrix(self%M_lowres)
      else
         call mpi_reduce(M, M,             size(M),  MPI_DOUBLE_PRECISION,  MPI_SUM,  0, self%comm, ierr)
      end if

      deallocate(M, dgrade, wl, wr)

      call update_status(status, "M_lowres_done")

    end subroutine precompute_M_lowres


  subroutine apply_wmap_precond(self, map, map_out)
    implicit none
    class(comm_WMAP_tod),              intent(in)    :: self
    real(dp),        dimension(0:,1:), intent(in)    :: map
    real(dp),        dimension(0:,1:), intent(out)   :: map_out

    integer(i4b) :: i, npix_lowres, n_lowres, nmaps
    real(dp) :: determ
    real(dp), allocatable, dimension(:)   :: m_lin, m

    if (self%comp_S) then
       map_out =  map/self%M_diag
    else
!!$
!!$       npix_lowres = 12*self%nside_M_lowres**2
!!$       nmaps       = self%nmaps_M_lowres
!!$
!!$       ! Apply lowres preconditioner
!!$       allocate(m_lin(0:npix_lowres*nmaps-1), m(0:size(map,1)-1))
!!$       do i = 1, nmaps
!!$          m = map(:,i)
!!$          call udgrade_ring(m, self%info%nside, m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres)
!!$          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, map_out(:,i), self%info%nside)
!!$       end do
!!$       m_lin = matmul(self%M_lowres, m_lin)
       
       ! Apply highres preconditioner to residual
       map_out = map !- map_out
       do i = 0, size(map,1)-1
          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
          map_out(i,1) =  map_out(i,1)/self%M_diag(i,1)
          map_out(i,2) = (map_out(i,2)*self%M_diag(i,3) - map_out(i,2)*self%M_diag(i,4))/determ
          map_out(i,3) = (map_out(i,3)*self%M_diag(i,2) - map_out(i,3)*self%M_diag(i,4))/determ
       end do

!!$       do i = 1, nmaps
!!$          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, m, self%info%nside)
!!$          map_out(:,i) = map_out(:,i) + m
!!$       end do
!!$
!!$       deallocate(m, m_lin)
       

!!$       do i = 0, size(map,1)-1
!!$          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
!!$          map_out(i,1) =  map(i,1)/self%M_diag(i,1)
!!$          map_out(i,2) = (map(i,2)*self%M_diag(i,3) - map(i,2)*self%M_diag(i,4))/determ
!!$          map_out(i,3) = (map(i,3)*self%M_diag(i,2) - map(i,3)*self%M_diag(i,4))/determ
!!$       end do
    end if

  end subroutine apply_wmap_precond

  subroutine sample_baseline(tod, scan, raw, s_tot, mask, handle)
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
    class(comm_wmap_tod),                   intent(inout) :: tod
    integer(i4b),                           intent(in)    :: scan
    real(sp),            dimension(1:,1:),  intent(in)    :: raw, s_tot, mask
    type(planck_rng),                       intent(inout) :: handle

    integer(i4b) :: i, j, k, n
    real(dp)     :: dt, t_tot, t, A, b, mval, eta
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
          end if
       end do

       call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline)
    end do

    deallocate(x, y)

  end subroutine sample_baseline

  subroutine construct_corrtemp_wmap(self, scan, pix, psi, s)
    !  Construct an WMAP instrument-specific correction template; for now contains 1Hz tem
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
    class(comm_wmap_tod),                  intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t

    dt = 1.d0 / self%scans(scan)%ntod
    t = 0.d0
    do j = 1, self%ndet
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          s(k,j) = 0.
          do i = 0, self%baseline_order
             s(k,j) = s(k,j) + self%scans(scan)%d(j)%baseline(i) * t**i
          end do
       end do
    end do

  end subroutine construct_corrtemp_wmap


end module comm_tod_WMAP_mod
