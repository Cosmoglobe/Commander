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
      character(len=512) :: noise_format
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

      real(dp),  allocatable, dimension(:, :)      :: m_buf



      call timer%start(TOD_INIT, id_abs)

      ! Initialize common parameters
      allocate (constructor)

      ! Set up noise PSD type and priors
      constructor%freq            = cpar%ds_label(id_abs)
      !constructor%n_xi            = 3
      !constructor%noise_psd_model = 'oof'

      constructor%n_xi            = 5
      constructor%noise_psd_model = 'oof_quad'
      constructor%comp_S          = .false.

      allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
      allocate(constructor%xi_n_nu_fit(constructor%n_xi,2))
      allocate(constructor%xi_n_P_rms(constructor%n_xi))
 
      ! Using Bennett 2013 x_im as initial guess for precomputing preconditioner 
      ! Jarosik 2003 Table 2 gives knee frequencies between 0.09 mHz and 
      ! 46.5 mHz. 
      !constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2]   ! [sigma0, fknee, alpha]; sigma0 is not used
      constructor%xi_n_P_rms      = [-1.0, 0.5, 0.5, -1.0, -1.0]   ! [sigma0, fknee, alpha, slope, intercept]; sigma0 is not used
      constructor%xi_n_P_uni(4,:) = [-0.5, 0.5]            ! slope
      constructor%xi_n_nu_fit(4,:) = [0.1, 1.0]       ! slope nu_fit
      constructor%xi_n_P_uni(5,:) = [-0.5, 0.5]             ! intercept
      constructor%xi_n_nu_fit(5,:) = [0.1, 1.0]       ! intercept nu_fit
      if (trim(constructor%freq) == '023-WMAP_K') then
         constructor%xi_n_nu_fit(2,:) = [0.0, 0.005]    
         constructor%xi_n_nu_fit(3,:) = [0.0, 0.005]    
         constructor%xi_n_P_uni(2,:) = [0.00001, 0.005]  ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [1.00e-05, 1.00e-05, 4.58e-03, 4.58e-03]
      else if (trim(constructor%freq) == '030-WMAP_Ka') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.005]    
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.005]    
         constructor%xi_n_P_uni(2,:) = [0.0001, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.00416, 0.00416, 0.00191, 0.00191]
      else if (trim(constructor%freq) == '040-WMAP_Q1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.010]    
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.010]    
         constructor%xi_n_P_uni(2,:) = [0.0001, 0.02]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.00139, 0.00139, 0.00581, 0.00581]
      else if (trim(constructor%freq) == '040-WMAP_Q2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.010]   
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.010]   
         constructor%xi_n_P_uni(2,:) = [0.0003, 0.02]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.00894, 0.00894, 0.01137, 0.01137]
      else if (trim(constructor%freq) == '060-WMAP_V1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.020]  
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.020]  
         constructor%xi_n_P_uni(2,:) = [0.0005, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.00153, 0.00153, 0.00415, 0.00415]
      else if (trim(constructor%freq) == '060-WMAP_V2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.020] 
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.020] 
         constructor%xi_n_P_uni(2,:) = [0.0005, 0.01]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.00233, 0.00233, 0.00243, 0.00243]
      else if (trim(constructor%freq) == '090-WMAP_W1') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.020]
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.020]
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.00]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.01211, 0.01211, 0.00414, 0.00414]
      else if (trim(constructor%freq) == '090-WMAP_W2') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.040]
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.040]
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.01255, 0.01255, 0.01736, 0.01736]
      else if (trim(constructor%freq) == '090-WMAP_W3') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.020] 
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.020] 
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [-0.00183, -0.00183,  0.00416,  0.00416]
      else if (trim(constructor%freq) == '090-WMAP_W4') then
         constructor%xi_n_nu_fit(2,:)     = [0.0, 0.080]  
         constructor%xi_n_nu_fit(3,:)     = [0.0, 0.080]  
         constructor%xi_n_P_uni(2,:) = [0.0005, 1.0]    ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.01]     ! alpha
         constructor%x_im = [0.02285, 0.02285, 0.02025, 0.02025]
      else
         write(*,*) 'Invalid WMAP frequency label = ', trim(constructor%freq)
         stop
      end if

      call constructor%tod_constructor(cpar, id_abs, info, tod_type)
      if (constructor%enable_tod_simulations) constructor%chisq_threshold = 1d6

      ! Set up WMAP specific parameters
      constructor%accept_threshold = 0.1d0 ! more stringent than default,
                                           ! cutting scans with more than 10% flagged
      constructor%samprate_lowres = 1.d0   ! Lowres samprate in Hz
      constructor%nhorn           = 2
      constructor%ndiode          = 1
      constructor%baseline_order  = 1
      ! Jarosik et al. uses a third-order baseline. How much of a difference
      ! would this make?
      ! constructor%baseline_order  = 3
      ! It turns out that the noise parameters get weird very quickly, with some
      ! bands immediately going to the boundaries.
      constructor%apply_inst_corr = .true.
      if (trim(constructor%level) == 'L1') then
          constructor%compressed_tod  = .true.
      else
          constructor%compressed_tod  = .false.
      end if
      constructor%correct_sl      = .true.
      constructor%orb_4pi_beam    = .true.
      constructor%symm_flags      = .false.
      constructor%chisq_threshold = 1000
      constructor%nmaps           = info%nmaps
      constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")
      constructor%verbosity       = cpar%verbosity

      ! Gain PSD Wiener filter parameters; determined by trial-and-error
      constructor%gain_tune_sigma0 = .false.
      constructor%gain_samprate    = 1.d0 / (24.d0*60.d0 * 60.d0)
      constructor%gain_sigma_0     = 3d-3 !3d-4                           ! Default from LFI
      constructor%gain_fknee       = constructor%gain_samprate      ! Default from LFI
      constructor%gain_alpha       = -1.d0                          ! Default from LFI

      if (constructor%myid == 0) then
         allocate(constructor%M_diag(0:info%npix-1,info%nmaps+1))
      end if

      ! Iniitialize TOD labels
      allocate (constructor%labels(8))
      constructor%labels(1) = 'map'
      constructor%labels(2) = 'res'
      constructor%labels(3) = 'ncorr'
      constructor%labels(4) = 'bpcorr'
      constructor%labels(5) = 'orb'
      constructor%labels(6) = 'sl'
      constructor%labels(7) = 'zodi'
      constructor%labels(8) = 'baseline'

      ! Initialize beams
      nside_beam                  = 512
      nmaps_beam                  = 3
      pol_beam                    = .true.
      constructor%nside_beam      = nside_beam


      ! Get detector labels
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)

      ! Get noise format
      constructor%noise_format   = cpar%ds_noise_format(id_abs)

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
      call constructor%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
      call constructor%precompute_lookups()

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

      ! Collect Sun velocities from all scals
      call constructor%collect_v_sun


      ! MJDs corresponding to August 10, 2001--2010
      constructor%split = (/52131,52496,52861,53227,53592,&
                          & 53957,54322,54688,55053,55418/)


      ! Need precompute the main beam precomputation for both the A-horn and
      ! B-horn.
      ! Allocate sidelobe convolution data structures
      allocate(constructor%slconvA(constructor%ndet), constructor%slconvB(constructor%ndet))
      allocate(constructor%orb_dp)
      if (constructor%orb_4pi_beam) then
         constructor%orb_dp => comm_orbdipole(beam=constructor%mbeam)
      else
         constructor%orb_dp => comm_orbdipole(comm=constructor%info%comm)
      end if

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

      real(dp)     :: t1, t2, monopole, sigma_mono
      integer(i4b) :: i, j, k, l, n
      integer(i4b) :: nside, npix, nmaps 
      integer(i4b) :: ierr, ndelta, t_mid=53765
      real(sp), allocatable, dimension(:, :)          :: s_buf
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: M_diag, buffer1
      real(dp), allocatable, dimension(:, :, :)       :: M_diag_1
      real(dp), allocatable, dimension(:)             :: II_inv, QQ_inv, UU_inv, QU_inv, det
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono, buffer2
      real(dp), allocatable, dimension(:,:,:,:)       :: b_map_1, b_map_2
      character(len=512) :: prefix, postfix
      character(len=2048) :: Sfilename

      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, bp_corr, output_scanlist, split
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
      real(dp), allocatable, dimension(:, :)    :: bicg_sol, bicg_sol_1, bicg_sol_2
      real(dp), allocatable, dimension(:, :)    :: map_full
      class(comm_map), pointer :: wmap_guess

      ! Counting data kept, lost
      integer(i8b) :: n_tot, n_flag, n_discard

      character(len=80), dimension(180) :: header


      type(hdf_file) :: tod_file
      real(sp) :: pow2  ! log2(ntod)



      ! Parameters used for testing
      real(dp) :: polang = 0d0



      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)


      call timer%start(TOD_ALLOC, self%band)

      ! Toggle optional operations
      sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
      bp_corr               = .true.                 ! by default, take into account differences in bandpasses. (WMAP does not do this in default analysis)
      bp_corr               = (bp_corr .or. sample_rel_bandpass) ! Bandpass is necessary to include if bandpass sampling is happening.
      select_data           = .false.         ! no data selection
      !select_data           = self%first_call ! only perform data selection the first time
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
      split = .false.
      if (self%output_aux_maps > 0) then
         if (self%first_call) then
           self%output_n_maps = 3
           split = .false.
         else
           if (mod(iter-1,10) == 0)  self%output_n_maps = 3
           if (mod(iter-1,20) == 0)  self%output_n_maps = 8
           if (mod(iter-1,100) == 0) split = .true.
         end if
      end if
      self%output_n_maps = 0



      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'

      ! Distribute maps
      ! Allocate total map (for monopole sampling)
      allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
      allocate(map_full(nmaps, 0:npix-1))
      !call distribute_sky_maps(self, map_in, 1.e-3, map_sky) ! uK to mK
      call distribute_sky_maps(self, map_in, 1., map_sky, map_full) ! K to K?



      ! Distribute processing masks
      allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
      call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
      call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
      deallocate(m_buf)

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
      if (split) then
          allocate (M_diag_1(0:npix-1, nmaps+1, 9))
          allocate ( b_map_1(0:npix-1, nmaps,   9, 1))
          M_diag_1 = 0d0
          b_map_1 = 0d0
      end if

      allocate(outmaps(1))
      outmaps(1)%p => comm_map(self%info)

      if (self%comp_S) then
        allocate (bicg_sol(0:npix-1, nmaps+1))
      else
        allocate (bicg_sol(0:npix-1, nmaps  ))
      end if
      if (split) then
          allocate (bicg_sol_1(0:npix-1, 9))
      end if
      call timer%stop(TOD_ALLOC, self%band)

      if (mod(iter-1, 10) == 0 .or. self%first_call) call self%precompute_M_lowres

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

      if (select_data) then 
         do i = 1, self%nscan
            self%scans(i)%d%accept = .true.
         end do
      end if


      ! Sample calibration
      if (.not. self%enable_tod_simulations) then
          if (trim(self%level) == 'L1') then
              ! Sample baseline for current scan
              if (self%myid == 0) then
                    write(*,*) '|    --> Sampling baseline'
              end if
              self%apply_inst_corr = .false. ! Disable baseline correction for just this call
              call update_status(status, "baseline")
              do i = 1, self%nscan
                 if (.not. any(self%scans(i)%d%accept)) cycle
                 call sd%init_differential(self, i, map_sky, procmask, procmask2, polang=polang)
                 call timer%start(TOD_BASELINE, self%band)
                 call sample_baseline(self, i, sd%tod, sd%s_tot, sd%mask, handle)
                 call timer%stop(TOD_BASELINE, self%band)
                 call sd%dealloc
              end do
              self%apply_inst_corr = .true.

              call timer%start(TOD_WAIT, self%band)
              call mpi_barrier(self%comm, ierr)
              call timer%stop(TOD_WAIT, self%band)

              call update_status(status, "abscal")
              call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2, polang)
              call update_status(status, "relcal")
              call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2, polang)
              call update_status(status, "deltaG")
              call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2, polang, smooth=.true.)
           else
              self%correct_sl      = .false.
              do j = 1, self%nscan
                 do i = 1, self%ndet
                    self%scans(j)%d(i)%gain = 1.d0
                 end do
              end do
              self%gain0(0) = 1
              self%gain0(1:) = 0
              self%x_im = 0
           end if
           call sample_calibration(self, 'imbal',  handle, map_sky, procmask, procmask2, polang)
      end if



      ! Perform loop over scans
      if (self%myid == 0) then
          if (self%enable_tod_simulations) then
            write(*,*) '|    --> Simulating TODs'
          else
            write(*,*) '|    --> Sampling ncorr, xi_n, maps'
          end if
      endif
      n_tot = 0
      n_discard = 0
      n_flag = 0
      do i = 1, self%nscan
         ! Skip scan if no accepted data

         if (.not. any(self%scans(i)%d%accept)) then
            if (output_scanlist) then
               write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                    & '"', 0.0, &
                    & real(self%spinaxis(i,:),sp)
            end if
            if (select_data) then 
               call sd%init_differential(self, i, map_sky, procmask, procmask2, &
                 & init_s_bp=bp_corr, polang=polang)
               n_tot = n_tot + sd%ntod/1000
               n_flag = n_flag + sd%ntod/1000
               call sd%dealloc
            end if
            cycle
         end if
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

         call timer%start(TOD_ALLOC, self%band)
         allocate(s_buf(sd%ntod,sd%ndet))
         call timer%stop(TOD_ALLOC, self%band)

         ! Make simulations or Sample correlated noise
         if (self%enable_tod_simulations) then
            call simulate_tod(self, i, sd%s_tot, sd%n_corr, handle)
         else
            call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,1,:), dospike=.false.)
            !sd%n_corr = 0.
         end if


         ! Compute noise spectrum parameters
         call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), &
              & sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
         end do

         ! Select data
         if (select_data) then 
            ! Count how many good data points are thrown out from this
            ! procedure.
            n_tot = n_tot + sd%ntod/1000
            n_flag = n_flag + count(iand(sd%flag(:,1),self%flag0) .ne. 0)/1000
            call remove_bad_data(self, i, sd%flag)
            if (.not. self%scans(i)%d(1)%accept) then
                n_discard = n_discard + (sd%ntod -count(iand(sd%flag(:,1),self%flag0) .ne. 0))/1000
            end if
            n = len(trim(self%freq)) - 1
            if ((self%freq(n:n) == 'W') .or. (self%freq(n:n) == 'V')) then
                if (sd%ntod < 2**22) then
                    write(*,*) '| Reject scan =', self%scanid(i), ' length ', sd%ntod
                    self%scans(i)%d%accept = .false.
                    n_discard = n_discard + (sd%ntod - count(iand(sd%flag(:,1),self%flag0) .ne. 0))/1000
                end if
            else
                if (sd%ntod < 2**21) then
                    write(*,*) '| Reject scan =', self%scanid(i), ' length ', sd%ntod
                    self%scans(i)%d%accept = .false.
                    n_discard = n_discard + (sd%ntod - count(iand(sd%flag(:,1),self%flag0) .ne. 0))/1000
                end if
            end if
         end if

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

         ! Compute binned map
         call timer%start(TOD_ALLOC, self%band)
         allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
         call timer%stop(TOD_ALLOC, self%band)

         call compute_calibrated_data(self, i, sd, d_calib)


         if (.false.) then
            call int2string(self%scanid(i), scantext)
            if (self%myid == 0 .and. i == 1) write(*,*) '| Writing tod to hdf'
            call open_hdf_file(trim(chaindir)//'/tod_'//scantext//'_samp'//samptext//'.h5', tod_file, 'w')
            !call write_hdf(tod_file, '/sl', sd%s_sl)
            !call write_hdf(tod_file, '/s_orb', sd%s_orb)
            call write_hdf(tod_file, '/n_corr', sd%n_corr)
            !call write_hdf(tod_file, '/bpcorr', sd%s_bp)
            call write_hdf(tod_file, '/s_tot', sd%s_tot)
            !call write_hdf(tod_file, '/s_sky', sd%s_sky)
            call write_hdf(tod_file, '/tod',   sd%tod)
            !call write_hdf(tod_file, '/flag', sd%flag)
            call write_hdf(tod_file, '/mask', sd%mask)
            !call write_hdf(tod_file, '/pixA', sd%pix(:,1,1))
            !call write_hdf(tod_file, '/pixB', sd%pix(:,1,2))
            !call write_hdf(tod_file, '/psiA', sd%psi(:,1,1))
            !call write_hdf(tod_file, '/psiB', sd%psi(:,1,2))
            !call write_hdf(tod_file, '/x_im', self%x_im)

            do k = 1, self%ndet
              call int2string(k, scantext)
              call write_hdf(tod_file, '/xi_n_'//scantext, self%scans(i)%d(k)%N_psd%xi_n)
              call write_hdf(tod_file, '/gain_'//scantext, self%scans(i)%d(k)%gain)
            end do

            call close_hdf_file(tod_file)


         end if
         
         ! Bin TOD
         ! Using procmask2 to only mask out the very brightest pixels.
         call bin_differential_TOD(self, d_calib, sd%pix(:,1,:),  &
           & sd%psi(:,1,:), sd%flag(:,1), self%x_im, procmask2, b_map, M_diag, i, &
           & self%comp_S)

         ! Temporal splits
         if (split) then

            do j = 1, 9
               if (self%scans(i)%t0(1) < self%split(j) .or. self%scans(i)%t0(1) > self%split(j+1)) cycle
               call bin_differential_TOD(self, d_calib(1:1,:,:), sd%pix(:,1,:),  &
                 & sd%psi(:,1,:), sd%flag(:,1), self%x_im, procmask2, b_map_1(:,:,j,:), M_diag_1(:,:,j), i, &
                 & self%comp_S)
            end do
         end if

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
         call timer%start(TOD_ALLOC, self%band)
         deallocate(s_buf, d_calib)
         call timer%stop(TOD_ALLOC, self%band)

      end do

      call timer%start(TOD_WAIT, self%band)
      call mpi_barrier(self%comm, ierr)
      call timer%stop(TOD_WAIT, self%band)
      if (.not. self%enable_tod_simulations) then
        if (self%myid == 0) write(*,*) '|    --> Finalizing binned maps'

        ! Output latest scan list with new timing information
        if (output_scanlist) call self%output_scan_list(slist)


        call timer%start(TOD_MPI, self%band)
        call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
             & MPI_DOUBLE_PRECISION, MPI_SUM,  self%info%comm, ierr)
        call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
             & MPI_DOUBLE_PRECISION, MPI_SUM,  self%info%comm, ierr)


        if (select_data) then
          call mpi_allreduce(mpi_in_place, n_tot, 1, &
               & MPI_INTEGER,  MPI_SUM, self%info%comm, ierr)
          call mpi_allreduce(mpi_in_place, n_flag, 1, &
               & MPI_INTEGER,  MPI_SUM, self%info%comm, ierr)
          call mpi_allreduce(mpi_in_place, n_discard, 1, &
               & MPI_INTEGER,  MPI_SUM, self%info%comm, ierr)
        end if
        call timer%stop(TOD_MPI, self%band)
        if (split) then
           call timer%start(TOD_MPI, self%band)
           call mpi_allreduce(mpi_in_place, b_map_1, size(b_map_1), &
                & MPI_DOUBLE_PRECISION, MPI_SUM,  self%info%comm, ierr)
           call mpi_allreduce(mpi_in_place, M_diag_1, size(M_diag_1), &
                & MPI_DOUBLE_PRECISION, MPI_SUM,  self%info%comm, ierr)
           call timer%stop(TOD_MPI, self%band)
        end if


        where (M_diag == 0d0)
           M_diag = 1d0
        end where
        if (self%myid == 0) self%M_diag = M_diag

        ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
        do l = self%output_n_maps, 1, -1
          !if (l .ne. 6) b_map(:,:,l) = 0d0
          if (l == 1) then
            bicg_sol = transpose(map_full)
            epsil = 1d-10
          else
            bicg_sol = 0d0
            epsil = 1d-6
          end if

          num_cg_iters = 0

          ! Solve for maps
          if (self%verbosity > 0 .and. self%myid == 0) then
            write(*,*) '|      Solving for ', trim(adjustl(self%labels(l)))
          end if
          call run_bicgstab(self, handle, bicg_sol, npix, nmaps, num_cg_iters, &
                         & epsil, procmask2, map_full, M_diag, b_map, l, &
                         & prefix, postfix, self%comp_S, 0)

          if (l == 1 .and. self%myid == 0) then
             ! Maximum likelihood monopole
             monopole = sum((bicg_sol(:,1)-map_full(1,:))*M_diag(:,1)*procmask) &
                    & / sum(M_diag(:,1)*procmask)
             if (trim(self%operation) == 'sample') then
                ! Add fluctuation term if requested
                sigma_mono = sum(M_diag(:,1) * procmask)
                if (sigma_mono > 0.d0) sigma_mono = 1.d0 / sqrt(sigma_mono)
                if (self%verbosity > 1) then
                  write(*,*) '|  monopole, fluctuation sigma'
                  write(*,*) '|  ', monopole, sigma_mono
                end if
                monopole = monopole + sigma_mono * rand_gauss(handle)
             end if
             bicg_sol(:,1) = bicg_sol(:,1) - monopole
          end if
          call mpi_bcast(bicg_sol, size(bicg_sol),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)

          if (split .and. l == 1) then
             do k = 1, 9
               if (self%verbosity > 0 .and. self%myid == 0) then
                 write(*,*) '|      Solving for map year ', k
               end if
               bicg_sol_1 = bicg_sol
               call run_bicgstab(self, handle, bicg_sol_1, npix, nmaps, num_cg_iters, &
                              & epsil, procmask2, map_full, M_diag_1(:,:,k), b_map_1(:,:,k,:), l, &
                              & prefix, postfix, self%comp_S, k)
               monopole = sum((bicg_sol_1(:,1)-map_full(1,:))*M_diag_1(:,1,k)*procmask) &
                      & / sum(M_diag_1(:,1,k)*procmask)

               call timer%start(TOD_WRITE) 
               call mpi_bcast(bicg_sol_1, size(bicg_sol_1),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
               do j = 1, nmaps
                  outmaps(1)%p%map(:, j) = bicg_sol_1(self%info%pix, j)
               end do
               map_out%map = outmaps(1)%p%map
               call int2string(k, ctext)
               call map_out%writeFITS(trim(prefix)//'map_yr'//ctext//trim(postfix))
               call timer%stop(TOD_WRITE) 
             end do
          end if
          
          call mpi_bcast(num_cg_iters, 1,  MPI_INTEGER, 0, self%info%comm, ierr)


          if (self%comp_S) then
             outmaps(1)%p%map(:,1) = bicg_sol(self%info%pix, nmaps+1)
             map_out%map = outmaps(1)%p%map
             call map_out%writeFITS(trim(prefix)//'S_'//trim(adjustl(self%labels(l)))//trim(postfix))
          end if


          do j = 1, nmaps
             outmaps(1)%p%map(:, j) = bicg_sol(self%info%pix, j)
          end do

          call timer%start(TOD_WRITE) 
          if (l == 1) then
             map_out%map = outmaps(1)%p%map
             call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))

             ! Recall:
             ! A = [[a, b],
             !      [c, d]]
             ! has the inverse
             ! A-1 = [[d, -b],
             !        [-c, a]]/(ad - bc)
             ! Note that if bc = 0, then A-1 is just the inverse of the
             ! diagonals.

             if (trim(self%noise_format) == 'rms_qucov') then
               allocate(det(0:npix-1))
               if (split) then
                 do k = 1, 9
                   det = M_diag_1(:,2,k)*M_diag_1(:,3,k) - M_diag_1(:,4,k)**2
                   rms_out%map(:,1:nmaps) = 1/sqrt(M_diag_1(self%info%pix, 1:nmaps,k))
                   rms_out%map(:,nmaps+1) = M_diag_1(self%info%pix, nmaps+1,k)
                   call int2string(k, ctext)
                   call rms_out%writeFITS(trim(prefix)//'rms_yr'//ctext//trim(postfix))
                 end do
               end if
               det = M_diag(:,2)*M_diag(:,3) - M_diag(:,4)**2
               rms_out%map(:,1) = 1/M_diag(self%info%pix, 1)
               rms_out%map(:,2) = M_diag(self%info%pix,3)/det(self%info%pix)
               rms_out%map(:,3) = M_diag(self%info%pix,2)/det(self%info%pix)
               rms_out%map(:,4) = -M_diag(self%info%pix,4)/det(self%info%pix)
               deallocate(det)
               call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
             else if (trim(self%noise_format) == 'rms') then
               rms_out%map(:,1:nmaps) = 1/sqrt(M_diag(self%info%pix, 1:nmaps))
               call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
             else
               if (self%myid == 0) write(*,*) 'unexpected noise format'
             end if
          else
             call outmaps(1)%p%writeFITS(trim(prefix)//trim(adjustl(self%labels(l)))//trim(postfix))
          end if
          call timer%stop(TOD_WRITE) 
        end do



        ! Sample bandpass parameters
        !if (sample_rel_bandpass .or. sample_abs_bandpass) then
        !   call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
        !   self%bp_delta = delta(:,:,1)
        !end if
      end if

      if (self%myid == 0  .and. select_data) then
        write(*,*) '| Data rejection statistics'
        write(*,*) '| n_tot     ', n_tot
        write(*,*) '| flagged   ', n_flag
        write(*,*) '| discarded ', n_discard
      end if

      ! Clean up temporary arrays

      call timer%start(TOD_ALLOC, self%band)
      deallocate(procmask, procmask2)
      deallocate(b_map, M_diag)
      deallocate(map_full)
      deallocate(bicg_sol)
      if (allocated(chisq_S)) deallocate (chisq_S)
      if (allocated(b_mono)) deallocate (b_mono)
      if (allocated(sys_mono)) deallocate (sys_mono)
      if (allocated(slist)) deallocate (slist)
      if (allocated(b_map_1)) deallocate(b_map_1,M_diag_1,bicg_sol_1)

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
      call timer%stop(TOD_ALLOC, self%band)

      call int2string(iter, ctext)
      call update_status(status, "tod_end"//ctext)

    call timer%stop(TOD_TOT, self%band)

    ! Parameter to check if this is first time routine has been called
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

      integer(i4b) :: i, j, k, t, p1, p2, p1_l,p1_r,p2_l,p2_r,k1, k2, ntot, npix, npix_hi, ierr, ntod, lpix, rpix, q, nhorn, lpsi, rpsi
      real(dp)     :: var, inv_sigma, lcos2psi, lsin2psi, rcos2psi, rsin2psi
      real(dp)     :: dx, xbar, f_l, f_r, mA, mB
      real(dp), allocatable, dimension(:)   :: dl, dr, pl, pr
      real(dp), allocatable, dimension(:,:) :: M
      integer(i4b), allocatable, dimension(:)         :: flag, dgrade
      real(sp),  allocatable, dimension(:)         :: procmask
      real(dp),  allocatable, dimension(:, :)      :: m_buf
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi
      type(hdf_file) :: precond_file



      call update_status(status, "M_lowres")
      if (self%myid == 0) write(*,*) '|    Computing preconditioner'

      self%nmaps_M_lowres = 3; if (self%comp_S) self%nmaps_M_lowres = 4
      self%nside_M_lowres = 8
      npix                = 12  *self%nside_M_lowres**2
      ntot                = npix*self%nmaps_M_lowres
      nhorn               = self%nhorn
      npix_hi             = 12  * self%info%nside**2
    
      allocate(m_buf(0:npix_hi-1,3), procmask(0:npix_hi-1))
      call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
      deallocate(m_buf)

      ! Computing the factors involving imbalance parameters
      dx   = (self%x_im(1) - self%x_im(3))*0.5
      xbar = (self%x_im(1) + self%x_im(3))*0.5
      !dx = 0
      !xbar = 0

      ! Precompute udgrade lookup table
      allocate(dgrade(0:12*self%info%nside**2-1))
      q = (self%info%nside / self%nside_M_lowres)**2
      do i = 0, 12*self%info%nside**2-1
         call ring2nest(self%info%nside, i, j)
         j = j/q
         call nest2ring(self%nside_M_lowres, j, dgrade(i))
      end do

      ! Allocate local coupling matrix
      allocate(M(0:ntot-1,0:ntot-1), dl(self%nmaps_M_lowres), dr(self%nmaps_M_lowres))
      allocate(pl(self%nmaps_M_lowres), pr(self%nmaps_M_lowres))
      M = 0.d0

      inv_sigma = 1
      f_l = 1
      f_r = 1
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
         ! 16 because each variable is divided by 4, variance goes as Var(aX) = a^2 Var(X)
         do k = 1, 4
            var = var + (self%scans(i)%d(k)%N_psd%sigma0/self%scans(i)%d(k)%gain)**2/16
         end do
         ! TODO
         inv_sigma = sqrt(1.d0/var)

         do t = 1, ntod
            if (iand(flag(t),self%flag0) .ne. 0) cycle
            lpix = dgrade(pix(t, 1))
            rpix = dgrade(pix(t, 2))
            lpsi = psi(t,1)
            rpsi = psi(t,2)

            f_l = procmask(pix(t,2))
            f_r = procmask(pix(t,1))

            dl(1) = 1+xbar
            dl(2) = dx * self%cos2psi(lpsi)
            dl(3) = dx * self%sin2psi(lpsi)
            dl    = dl * inv_sigma * f_l

            dr(1) = -(1-xbar)
            dr(2) = dx * self%cos2psi(rpsi)
            dr(3) = dx * self%sin2psi(rpsi)
            dr    = dr * inv_sigma * f_r


            pl(1) = dx
            pl(2) = (1+xbar) * self%cos2psi(lpsi)
            pl(3) = (1+xbar) * self%sin2psi(lpsi)
            pl    = pl * inv_sigma * f_l

            pr(1) = dx
            pr(2) = -(1-xbar) * self%cos2psi(rpsi)
            pr(3) = -(1-xbar) * self%sin2psi(rpsi)
            pr    = pr * inv_sigma * f_r

            do k1 = 1, self%nmaps_M_lowres
               p1_l = (k1-1)*npix + lpix
               p1_r = (k1-1)*npix + rpix
               do k2 = 1, self%nmaps_M_lowres
                  p2_l = (k2-1)*npix + lpix
                  p2_r = (k2-1)*npix + rpix
                  !write(*,*) p1_l, p1_r, k1, p2_l, p2_r, k2

                  if ((k1 .eq. 1 .and. k2 .eq. 1) .or. (k1 > 1 .and. k2 > 1)) then
                      ! Intensity
                      ! P_A N^-1 P_A
                      M(p1_l,p2_l) = M(p1_l,p2_l) + dl(k1) * dl(k2) 
                      ! P_B N^-1 P_B
                      M(p1_r,p2_r) = M(p1_r,p2_r) + dr(k1) * dr(k2) 
                      ! P_A N^-1 P_B
                      M(p1_l,p2_r) = M(p1_l,p2_r) + dl(k1) * dr(k2) 
                      ! P_B N^-1 P_A
                      M(p1_r,p2_l) = M(p1_r,p2_l) + dr(k1) * dl(k2) 

                      ! Polarization
                      M(p1_l,p2_l) = M(p1_l,p2_l) + pl(k1) * pl(k2) 
                      M(p1_r,p2_r) = M(p1_r,p2_r) + pr(k1) * pr(k2) 
                      M(p1_l,p2_r) = M(p1_l,p2_r) + pl(k1) * pr(k2) 
                      M(p1_r,p2_l) = M(p1_r,p2_l) + pr(k1) * pl(k2) 
                 end if
               end do
            end do

         end do

         deallocate(pix, psi, flag)
      end do

      call timer%start(TOD_WAIT, self%band)
      call mpi_barrier(self%comm, ierr)
      call timer%stop(TOD_WAIT, self%band)

      ! Collect contributions from all cores 
      if (self%myid == 0) then
         write(*,*) '|    Inverting preconditioner'
         if (.not. allocated(self%M_lowres)) allocate(self%M_lowres(ntot,ntot))
         call mpi_reduce(M, self%M_lowres, size(M),  MPI_DOUBLE_PRECISION,  MPI_SUM,  0, self%comm, ierr)

         if (.false.) then
             call open_hdf_file('precond_'//trim(self%freq)//'.h5', precond_file, 'w')
             call write_hdf(precond_file, '/M', self%M_lowres)
             call close_hdf_file(precond_file)
         end if

         call invert_matrix(self%M_lowres)
      else
         call mpi_reduce(M, M,             size(M),  MPI_DOUBLE_PRECISION,  MPI_SUM,  0, self%comm, ierr)
      end if

      deallocate(M, dgrade, dl, dr, pl, pr)


      call update_status(status, "M_lowres_done")

    end subroutine precompute_M_lowres


  subroutine apply_wmap_precond(self, map, map_out)
    implicit none
    class(comm_WMAP_tod),              intent(in)    :: self
    real(dp),        dimension(0:,1:), intent(in)    :: map
    real(dp),        dimension(0:,1:), intent(out)   :: map_out

    integer(i4b) :: i, npix_lowres, n_lowres, nmaps, n, p, q, j
    real(dp) :: determ
    real(dp), allocatable, dimension(:)   :: m_lin, m
    real(dp), allocatable, dimension(:,:) :: m_low
    !
    !   Routine follows Section 3.4.7 of Jarosik et al. 2007
    !
    !   Arguments: 
    !   ----------
    !   self:     pointer of comm_WMAP_tod class
    !             Points to output of the constructor
    !   map:      Map to be preconditioned
    !
    !   map_out:  Map after being preconditioned
    !
    !   Intermediate steps:
    !   -------------------
    !   m_lin - a linearized map, length nmaps * npix

    if (self%comp_S) then
       map_out =  map/self%M_diag
    else

       map_out = 0d0

       npix_lowres = 12*self%nside_M_lowres**2
       nmaps       = self%nmaps_M_lowres

       ! Apply lowres preconditioner
       allocate(m_lin(0:npix_lowres*nmaps-1), m(0:size(map,1)-1))
       allocate(m_low(0:size(map,1)-1, nmaps))
       do i = 1, nmaps
          m = map(:,i)
          call udgrade_ring(m, self%info%nside, m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres)
       end do
       ! m_lin is now the low resolution linearized version of the map
       m_lin = matmul(self%M_lowres, m_lin)
       ! m_lin has now been preconditioned.

       do i = 1, nmaps
          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, m_low(:,i), self%info%nside)
       end do
       
       ! Apply highres preconditioner to residual
       map_out = map - m_low
       do i = 0, size(map,1)-1
          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
          map_out(i,1) =  map_out(i,1)/self%M_diag(i,1)
          map_out(i,2) = (map_out(i,2)*self%M_diag(i,3) - map_out(i,2)*self%M_diag(i,4))/determ
          map_out(i,3) = (map_out(i,3)*self%M_diag(i,2) - map_out(i,3)*self%M_diag(i,4))/determ
       end do

       do i = 1, nmaps
          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, m, self%info%nside)
          map_out(:,i) = map_out(:,i) + m_low(:,i)
       end do

       deallocate(m, m_lin, m_low)
       

!       do i = 0, size(map,1)-1
!          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
!          map_out(i,1) =  map(i,1)/self%M_diag(i,1)
!          map_out(i,2) = (map(i,2)*self%M_diag(i,3) - map(i,2)*self%M_diag(i,4))/determ
!          map_out(i,3) = (map(i,3)*self%M_diag(i,2) - map(i,3)*self%M_diag(i,4))/determ
!       end do


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
    !  Construct an WMAP instrument-specific correction template; for now contains baseline
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
    s = 0.
    do j = 1, self%ndet
       t = 0.d0
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          do i = 0, self%baseline_order
             s(k,j) = s(k,j) + self%scans(scan)%d(j)%baseline(i) * t**i
          end do
       end do
    end do


  end subroutine construct_corrtemp_wmap


end module comm_tod_WMAP_mod
