!===============================================================================
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
module comm_tod_mod
  use comm_fft_mod
  use comm_huffman_mod
  use comm_conviqt_mod
  use comm_tod_cray_mod
  use comm_tod_orbdipole_mod
  use comm_tod_noise_psd_mod
  use comm_shared_arr_mod
  use comm_utils
  USE ISO_C_BINDING
  implicit none

  private
  public comm_tod, comm_scan, initialize_tod_mod, fill_masked_region, fill_all_masked, tod_pointer, distribute_sky_maps

  ! Structure for individual detectors
  type :: comm_detscan
     character(len=30) :: label                             ! Detector label
     real(dp)          :: gain, dgain, gain_invsigma        ! Gain; assumed constant over scan
     real(dp)          :: gain_def                          ! Default parameters
     real(dp)          :: chisq
     real(dp)          :: chisq_prop
     real(dp)          :: chisq_masked
     real(dp)          :: baseline1, baseline2
     logical(lgt)      :: accept
     class(comm_noise_psd), pointer :: N_psd                            ! Noise PSD object
     real(sp),           allocatable, dimension(:)     :: tod            ! Detector values in time domain, (ntod)
     byte,               allocatable, dimension(:)     :: ztod           ! compressed values in time domain, (ntod)
     real(sp),           allocatable, dimension(:,:)   :: diode          ! (ndiode, ntod) array of undifferenced data
     type(byte_pointer), allocatable, dimension(:)     :: zdiode         ! pointers to the compressed undifferenced diode data, len (ndiode)
     byte,               allocatable, dimension(:)     :: flag           ! Compressed detector flag; 0 is accepted, /= 0 is rejected
     integer(i4b),       allocatable, dimension(:,:)   :: mask_dyn       ! Dynamic online-generated mask, (2,ntod), each row gives a range of masked samples
     type(byte_pointer), allocatable, dimension(:)     :: pix            ! pointer array of pixels length nhorn
     type(byte_pointer), allocatable, dimension(:)     :: psi            ! pointer array of psi, length nhorn
     integer(i4b),       allocatable, dimension(:,:)   :: offset_range   ! Beginning and end tod index of every offset region
     real(sp),           allocatable, dimension(:)     :: offset_level   ! Amplitude of every offset region(step)
     integer(i4b),       allocatable, dimension(:,:)   :: jumpflag_range ! Beginning and end tod index of regions where jumps occur
     real(dp),           allocatable, dimension(:)     :: baseline       ! Polynomial coefficients for baseline function
     integer(i4b),       allocatable, dimension(:,:)   :: pix_sol        ! Discretized pointing in solar centric coordinates, for zodi and sidelobe mapping
     integer(i4b),       allocatable, dimension(:,:)   :: psi_sol        ! Discretized polarization angle in solar centric coordinates, for zodi and sidelobe mapping

     ! Zodi sampling structures (downsampled and precomputed quantities. only allocated if zodi sampling is true)
     logical(lgt),       allocatable, dimension(:)    :: zodi_glitch_mask
     integer(i4b),       allocatable, dimension(:)    :: downsamp_pix
     real(sp),           allocatable, dimension(:)    :: downsamp_tod
     real(sp),           allocatable, dimension(:)    :: downsamp_sky
     real(sp),           allocatable, dimension(:)    :: downsamp_zodi
     real(sp),           allocatable, dimension(:, :) :: downsamp_scat
     real(sp),           allocatable, dimension(:, :) :: downsamp_therm
     real(sp),           allocatable, dimension(:, :) :: downsamp_point  ! (ntod,{lat_gal, lon_gal, lat_ecl, lon_ecl, solar elongation}

     real(sp),           allocatable, dimension(:, :) :: s_scat_lowres
     real(sp),           allocatable, dimension(:, :) :: s_therm_lowres
  end type comm_detscan

  ! Stores information about all detectors at once 
  type :: comm_scan
     integer(i4b)   :: ntod                                        ! Number of time samples
     integer(i4b)   :: ext_lowres(2)                            ! Shape of downgraded TOD including padding
     real(dp)       :: proctime    = 0.d0                          ! Processing time in seconds
     real(dp)       :: n_proctime  = 0                             ! Number of completed loops
     real(dp)       :: v_sun(3)                                    ! Observatory velocity relative to Sun in km/s
     real(dp)       :: t0(3)                                       ! MJD, OBT, SCET for start of chunk
     real(dp)       :: t1(3)                                       ! MJD, OBT, SCET for end of chunk
     real(dp)       :: x0_obs(3)                                   ! Observatory position (x,y,z) for start of chunk
     real(dp)       :: x1_obs(3)                                   ! Observatory position (x,y,z) for end of chunk
     real(dp)       :: x0_earth(3)                                 ! Observatory position (x,y,z) for start of chunk
     real(dp)       :: x1_earth(3)                                 ! Observatory position (x,y,z) for end of chunk

     type(huffcode) :: hkey                                        ! Huffman decompression key
     type(huffcode) :: todkey                                      ! Huffman decompression key
     integer(i4b)   :: chunk_num                                   ! Absolute number of chunk in the data files
     integer(i4b),        allocatable, dimension(:,:)   :: zext    ! Extension of compressed diode arrays
     real(sp),            allocatable, dimension(:)     :: downsamp_obs_time ! downsampled_obs_time used for zodi sampling
     class(comm_detscan), allocatable, dimension(:)     :: d       ! Array of all detectors
  end type comm_scan


  type, abstract :: comm_tod
     character(len=512) :: freq
     character(len=512) :: instlabel
     character(len=512) :: filelist
     character(len=512) :: procmaskf1
     character(len=512) :: procmaskf2
     character(len=512) :: procmaskfzodi
     character(len=512) :: initfile
     character(len=512) :: instfile
     character(len=512) :: operation
     character(len=512) :: outdir
     character(len=512) :: sims_output_dir  !< simulation folder
     character(len=512) :: noise_psd_model  
     character(len=512) :: level !which level of tod we want, L1 or L2
     logical(lgt) :: enable_tod_simulations !< simulation parameter to run commander3 in different regime
     logical(lgt) :: first_call
     logical(lgt) :: sample_L1_par                                ! If false, reduce L1 (diode) to L2 (detector) in precomputations
     logical(lgt) :: L2_exist
     character(len=512) :: L2file
     integer(i4b) :: comm, myid, numprocs                         ! MPI parameters
     integer(i4b) :: comm_shared, myid_shared, numprocs_shared    ! MPI parameters
     integer(i4b) :: comm_inter, myid_inter                       ! MPI parameters
     integer(i4b) :: band                                        ! Absolute band ID
     integer(i4b) :: id                                          ! Relative band ID
     integer(i4b) :: zodiband                                        ! Band ID for zodi
     integer(i4b) :: nmaps                                        ! Number of Stokes parameters
     integer(i4b) :: ndet                                         ! Number of active detectors
     integer(i4b) :: nhorn                                        ! Number of horns
     integer(i4b) :: ndiode                                      ! Number of diodes that makeup each detector
     character(len=10), allocatable, dimension(:,:)  :: diode_names  ! Names of each diode, (ndet, ndiode)
     integer(i4b) :: nscan, nscan_tot                             ! Number of scans
     integer(i4b) :: first_scan, last_scan
     integer(i4b) :: npsi                                         ! Number of discretized psi steps
     integer(i4b) :: flag0
     integer(i4b) :: n_xi                                         ! Number of noise parameters
     integer(i4b) :: ntime                                        ! Number of time values
     integer(i4b) :: n_cray_temps = 0                             ! number of classes of cosmic rays we have
     integer(i4b) :: baseline_order                               ! Polynomial order for baseline
     real(dp)     :: central_freq                                 !Central frequency
     real(dp)     :: samprate, samprate_lowres                    ! Sample rate in Hz
     real(dp)     :: chisq_threshold                              ! Quality threshold in sigma
     character(len=512) :: abscal_comps            ! List of components to calibrate against
     logical(lgt) :: compressed_tod               
     logical(lgt) :: apply_inst_corr               
     logical(lgt) :: sample_abs_bp
     logical(lgt) :: symm_flags               
     class(comm_orbdipole), pointer :: orb_dp
     real(dp), allocatable, dimension(:)     :: gain0                                      ! Mean gain
     real(dp), allocatable, dimension(:)     :: polang                                      ! Detector polarization angle
     real(dp), allocatable, dimension(:,:)   :: polang_prior
        ! Detector polarization angle prior [ndet,mean/rms]
     real(dp), allocatable, dimension(:)     :: mbang                                       ! Main beams angle
     real(dp), allocatable, dimension(:)     :: mono                                        ! Monopole
     real(dp), allocatable, dimension(:)     :: fwhm, elip, psi_ell                         ! Beam parameter
     real(dp), allocatable, dimension(:)     :: nu_c                                        ! Center frequency
     real(dp), allocatable, dimension(:,:,:) :: prop_bp         ! proposal matrix, L(ndet,ndet,ndelta),  for bandpass sampler
     real(dp), allocatable, dimension(:)     :: prop_bp_mean    ! proposal matrix, sigma(ndelta), for mean
     real(sp), allocatable, dimension(:,:)   :: xi_n_P_uni      ! Uniform prior for noise PSD parameters
     real(sp), allocatable, dimension(:)     :: xi_n_P_rms      ! RMS for active noise PSD prior
     real(sp), allocatable, dimension(:,:)   :: xi_n_nu_fit     ! Frequency range used to fit noise PSD parameters, (xi_n, 2)
     integer(i4b)      :: nside, nside_param                    ! Nside for pixelized pointing
     integer(i4b)      :: nobs, nobs_lowres                     ! Number of observed pixels for this core
     integer(i4b)      :: n_bp_prop                       ! Number of consecutive bandpass proposals in each main iteration; should be 2 for MH
     integer(i4b) :: output_n_maps                                ! Output n_maps
     character(len=512) :: init_from_HDF                          ! Read from HDF file
     character(len=512) :: datadir
     integer(i4b) :: output_4D_map                                ! Output 4D maps
     integer(i4b) :: output_aux_maps                              ! Output auxiliary maps
     integer(i4b) :: halfring_split                               ! Type of halfring split 0=None, 1=HR1, 2=HR2
     logical(lgt) :: subtract_zodi                                ! Subtract zodical light (defined in the parameter file)
     logical(lgt) :: sample_zodi                                  ! Sample zodi model parameters (defined in the parameter file)
     logical(lgt) :: output_zodi_comps                            ! Output zodi components
     logical(lgt) :: use_solar_point                              ! Compute solar centric pointing, for zodi or sidelobe mapping
     real(sp)     :: sol_elong_range(2)                           ! Acceptable solar elongation range
     logical(lgt) :: correct_sl                                   ! Subtract sidelobes
     logical(lgt) :: correct_orb                                  ! Subtract CMB dipole
     logical(lgt) :: sample_mono                                  ! Subtract detector-specific monopoles
     logical(lgt) :: orb_4pi_beam                                 ! Perform 4pi beam convolution for orbital CMB dipole 
     integer(i4b),       allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),           allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     real(sp),           allocatable, dimension(:)     :: sin2psi  ! Lookup table of sin(2psi)
     real(sp),           allocatable, dimension(:)     :: cos2psi  ! Lookup table of cos(2psi)
     real(sp),           allocatable, dimension(:)     :: psi      ! Lookup table of psi
     real(dp),           allocatable, dimension(:,:)   :: L_prop_mono  ! Proposal matrix for monopole sampling
     real(dp),           allocatable, dimension(:,:)   :: satpos   ! Satellite position for all scans
     real(dp),           allocatable, dimension(:)     :: mjds     ! MJDs for all scans(nscan_tot)
     real(dp),           allocatable, dimension(:,:)   :: v_sun    ! Sun velocities for all scans (3, nscan_tot)
     type(spline_type)                                 :: x_obs_spline(3), x_earth_spline(3) ! splines to compute observer and earth positions
     type(comm_scan),    allocatable, dimension(:)     :: scans    ! Array of all scans
     integer(i4b),       allocatable, dimension(:)     :: scanid   ! List of scan IDs
     integer(i4b),       allocatable, dimension(:)     :: nscanprproc   ! List of scan IDs
     integer(i4b),       allocatable, dimension(:)     :: partner ! Partner detector; for symmetrizing flags
     integer(i4b),       allocatable, dimension(:)     :: horn_id  ! Internal horn number per detector
     real(dp),           dimension(4)                  :: x_im    ! feedhorn imbalance parameters, with duplicates
     character(len=512), allocatable, dimension(:)     :: hdfname  ! List of HDF filenames for each ID
     character(len=512), allocatable, dimension(:)     :: label    ! Detector labels
     class(comm_map), pointer                          :: procmask => null() ! Mask for gain and n_corr
     class(comm_map), pointer                          :: procmask2 => null() ! Mask for gain and n_corr
     class(comm_map), pointer                          :: procmask_zodi => null() ! Mask for sampling zodi
     !class(comm_map), pointer                          :: mask_solar => null() ! Solar centric/sidelobe mask
     real(dp),           allocatable, dimension(:,:)   :: mask_solar           ! Solar centric/sidelobe mask
     logical(lgt)                                      :: map_solar_allocated
     real(dp),           pointer,     dimension(:,:)   :: map_solar           ! Full-sky solar centric/sidelobe model
!     class(comm_map), pointer                          :: map_solar => null() ! Solar centric/sidelobe model
     class(comm_mapinfo), pointer                      :: info => null()    ! Map definition
     class(comm_mapinfo), pointer                      :: slinfo => null()  ! Sidelobe map info
     class(comm_mapinfo), pointer                      :: mbinfo => null()  ! Main beam map info
     class(map_ptr),     allocatable, dimension(:)     :: slbeam, mbeam   ! Sidelobe beam data (ndet)
     class(conviqt_ptr), allocatable, dimension(:)     :: slconv   ! SL-convolved maps (ndet)
     class(cray_ptr),    allocatable, dimension(:)     :: cray ! cosmic ray templates
     class(conviqt_ptr), allocatable, dimension(:)     :: slconvA, slconvB ! SL-convolved maps (ndet)
     real(dp),           allocatable, dimension(:,:)   :: bp_delta  ! Bandpass parameters (0:ndet, npar)
     real(dp),           allocatable, dimension(:,:)   :: spinaxis ! For load balancing
     integer(i4b),       allocatable, dimension(:)     :: pix2ind ! Array mapping all npix pixels to the uniquely observed pixels in the tod object for saving memory
     integer(i4b),       allocatable, dimension(:)     :: ind2pix, ind2sl ! Lookup tables used with pix2ind 
     real(sp),           allocatable, dimension(:,:)   :: ind2ang ! Lookup tables used with pix2ind for pixel angles
     real(dp),           allocatable, dimension(:,:)   :: ind2vec ! Lookup tables used with pix2ind for pixel unit vectors
     real(dp),           allocatable, dimension(:,:)   :: ind2vec_ecl ! Lookuptable for lowres ind to ecliptic unit vector
     real(dp),           allocatable, dimension(:,:)   :: ind2vec_ecl_lowres ! Lookuptable for lowres ind to ecliptic unit vector
     integer(i4b),       allocatable, dimension(:)     :: udgrade_pix_zodi !Lookuptable for highres pix to lowres pix
     integer(i4b),       allocatable, dimension(:)     :: pix2ind_lowres !Lookuptable for lowres zodi pixels

     character(len=128)                                :: tod_type
     integer(i4b)                                      :: nside_beam
     integer(i4b)                                      :: verbosity ! verbosity of output
     integer(i4b),       allocatable, dimension(:,:)   :: jumplist  ! List of stationary periods (ndet,njump+2)
     real(dp)                                          :: accept_threshold ! Required fraction of unflagged data in a detscan in order to be accepted 
     logical(lgt)                                      :: orbital ! flag for whether the orbital or solar dipole is used as the template in construct_dipole_template()
     ! Gain parameters
     logical(lgt)                            :: gain_tune_sigma0
     real(dp)                                :: gain_samprate
     real(dp), allocatable, dimension(:)     :: gain_sigma_0  ! size(ndet), the estimated white noise level of that scan. Not truly a white noise since our model is sigma_0**2 * (f/fknee)**alpha instead of sigma_0 ** 2 (1 + f/fknee ** alpha)
     real(dp), allocatable, dimension(:)    :: gain_fknee ! size(ndet)
     real(dp), allocatable, dimension(:)    :: gain_alpha ! size(ndet)
     real(dp) :: gain_sigma0_std ! std for metropolis-hastings sampling
     real(dp) :: gain_fknee_std ! std for metropolis-hastings sampling
     real(dp) :: gain_alpha_std ! std for metropolis-hastings sampling
     integer(i4b), allocatable, dimension(:) :: split

     ! Zodi parameters and spline objects
     integer(i4b) :: zodi_n_comps
   !   real(sp), allocatable, dimension(:, :, :) :: zodi_scat_cache, zodi_therm_cache ! Cached s_zodi array for a given processor
     real(sp), allocatable, dimension(:, :, :) :: zodi_scat_cache, zodi_therm_cache ! Cache for zodi
     real(sp), allocatable, dimension(:, :, :) :: zodi_scat_cache_lowres, zodi_therm_cache_lowres ! Cache for lowresolution zodi (used for sampling)
     real(dp)                                  :: zodi_cache_time, zodi_init_cache_time! Time of cached zodi array
     !real(dp), allocatable, dimension(:)       :: zodi_emissivity, zodi_albedo ! sampled parameters
     real(dp), allocatable, dimension(:, :)    :: zodi_spl_phase_coeffs
     real(dp), allocatable, dimension(:)       :: zodi_spl_solar_irradiance, zodi_phase_func_normalization
     type(spline_type), allocatable            :: zodi_b_nu_spl_obj(:)
     logical(lgt)                              :: zodi_tod_params_are_initialized, zodi_scattering, udgrade_zodi
   contains
     procedure                           :: read_tod
     procedure                           :: diode2tod_inst
     procedure                           :: read_tod_inst
     procedure                           :: read_scan_inst
     procedure                           :: get_scan_ids
     procedure                           :: dumpToHDF
     procedure                           :: dumpToHDF_inst
     procedure                           :: initHDF
     procedure                           :: initHDF_inst
     procedure                           :: get_det_id
     procedure                           :: initialize_bp_covar
     procedure(process_tod), deferred    :: process_tod
     procedure                           :: construct_sl_template
     procedure                           :: construct_corrtemp_inst
     procedure                           :: construct_dipole_template
     procedure                           :: construct_dipole_template_diff
     procedure                           :: output_scan_list
     procedure                           :: downsample_tod
     procedure                           :: compute_tod_chisq
     procedure                           :: get_total_chisq
     procedure                           :: symmetrize_flags
     procedure                           :: decompress_pointing_and_flags
     procedure                           :: decompress_tod
     procedure                           :: decompress_diodes
     procedure                           :: tod_constructor
     procedure                           :: load_instrument_file
     procedure                           :: load_instrument_inst
     procedure                           :: precompute_lookups
     procedure                           :: read_jumplist
     procedure                           :: remove_fixed_scans
     procedure                           :: apply_map_precond
     procedure                           :: collect_v_sun
     procedure                           :: precompute_zodi_lookups
     procedure                           :: clear_zodi_cache
     procedure                           :: create_dynamic_mask
     procedure                           :: get_s_static
  end type comm_tod
  
  abstract interface
     subroutine process_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
       import i4b, comm_tod, comm_map, map_ptr, dp, planck_rng
       implicit none
       class(comm_tod),                     intent(inout) :: self
       character(len=*),                    intent(in)    :: chaindir
       integer(i4b),                        intent(in)    :: chain, iter
       type(planck_rng),                    intent(inout) :: handle
       type(map_ptr),     dimension(:,:),   intent(inout) :: map_in
       real(dp),          dimension(:,:,:), intent(inout) :: delta
       class(comm_map),                     intent(inout) :: map_out
       class(comm_map),                     intent(inout) :: rms_out
       type(map_ptr),     dimension(:,:),   intent(inout), optional :: map_gain
     end subroutine process_tod
  end interface

  type tod_pointer
    class(comm_tod), pointer :: p => null()
  end type tod_pointer
  
contains

  subroutine initialize_tod_mod(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar

    logical(lgt), save :: initialized = .false.

    if (initialized) return

    call initialize_fft_mod(cpar)
  end subroutine initialize_tod_mod

  subroutine tod_constructor(self, cpar, id, id_abs, info, tod_type)
    ! 
    ! Common constructor function for all TOD objects; allocatates and initializes general
    ! data structures. This routine is typically called from within an instrument-specific 
    ! initialization routine, *after* defining fields such that nhorn, ndet etc.
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object to be initialized
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    ! id:       integer
    !           The index of the current band ampng accepted bands
    ! id_abs:   integer
    !           The index of the current band within the parameters, related to cpar

    ! info:     map_info structure
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_tod),                intent(inout)  :: self
    integer(i4b),                   intent(in)     :: id, id_abs
    type(comm_params),              intent(in)     :: cpar
    class(comm_mapinfo),            target         :: info
    character(len=128),             intent(in)     :: tod_type

    integer(i4b) :: i, ndelta, ierr, unit
    character(len=512) :: datadir, solar_init


  end subroutine tod_constructor

  subroutine precompute_lookups(self)
    ! 
    ! Routine that precomputes static look-up tables in a given TOD object (pix2ind, ind2pix, ind2sl, ind2ang). 
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object to be initialized
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_tod),                intent(inout)  :: self

    real(dp)     :: f_fill, f_fill_lim(3), theta, phi
    integer(i4b) :: i, j, k, l, ierr
    integer(i4b), allocatable, dimension(:) :: pix


  end subroutine precompute_lookups

  !**************************************************
  !             Utility routines
  !**************************************************

  subroutine load_instrument_file(self, nside_beam, nmaps_beam, pol_beam, comm_chain)
    implicit none
    class(comm_tod),   intent(inout) :: self
    integer(i4b),      intent(in)    :: nside_beam
    integer(i4b),      intent(in)    :: nmaps_beam
    logical(lgt),      intent(in)    :: pol_beam
    integer(i4b),      intent(in)    :: comm_chain 

    type(hdf_file) :: h5_file
    integer(i4b) :: lmax_beam, lmax_sl, i
    type(comm_mapinfo), pointer :: info_beam


  end subroutine load_instrument_file


  subroutine read_tod(self, detlabels, datadir)
    ! 
    ! Reads common TOD fields into existing TOD object
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object to be initialized
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! Returns
    ! ----------
    ! None, but updates self
    !
   implicit none
   class(comm_tod),                intent(inout)  :: self
   character(len=*), dimension(:), intent(in)     :: detlabels
   character(len=512), optional                   :: datadir
   
   


  end subroutine read_tod

  subroutine read_hdf_scan(self, tod, filename, scan, ndet, detlabels, nhorn, ndiode, diode_names)
    ! 
    ! Reads common scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scan)
    !           Scan object
    ! tod:      derived class (comm_tod)
    !           Main TOD object to which current scan belongs
    ! filename: character
    !           TOD filename
    ! scan:     int
    !           Scan ID
    ! ndet:     int
    !           Number of detectors
    ! nhorn:    int
    !           Number of horns
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! ndiode:   int 
    !           Number of diodes per combined tod
    ! diode_names : string (array (ndet, ndiode)
    !           Array of diode labels, eg. [['sky00', 'sky01', 'load00',
    !           'load01'], ['sky10', 'sky11', 'load10', 'load11'], ...]
    ! scan:     derived class (comm_scan)
    !           
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    ! TODO
    ! ----
    ! - ndet, nhorn and detlabels should be taken from tod, not inserted as separate parameters?
    ! 
    implicit none
    class(comm_scan),               intent(inout) :: self
    class(comm_tod),                intent(in)    :: tod
    character(len=*),               intent(in)    :: filename
    integer(i4b),                   intent(in)    :: scan, ndet, nhorn, ndiode
    character(len=*), dimension(:), intent(in)    :: detlabels
    character(len=*), dimension(:,:), intent(in)  :: diode_names

    integer(i4b)       :: i,j,k,l, n, m, ext(1), setsize(1)
    real(sp)           :: nu
    real(dp)           :: scalars(4), time
    character(len=6)   :: slabel
    character(len=128) :: field
    type(hdf_file)     :: file
    integer(i4b), allocatable, dimension(:)       :: hsymb
    real(sp),     allocatable, dimension(:)       :: buffer_sp, xi_n, hsymb_sp
    integer(i4b), allocatable, dimension(:)       :: htree


  end subroutine read_hdf_scan

  subroutine read_hdf_scan_data(self, tod, filename, scan, ndet, detlabels, nhorn, ndiode, diode_names)
    ! 
    ! Reads common scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scan)
    !           Scan object
    ! tod:      derived class (comm_tod)
    !           Main TOD object to which current scan belongs
    ! filename: character
    !           TOD filename
    ! scan:     int
    !           Scan ID
    ! ndet:     int
    !           Number of detectors
    ! nhorn:    int
    !           Number of horns
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! ndiode:   int 
    !           Number of diodes per combined tod
    ! diode_names : string (array (ndet, ndiode)
    !           Array of diode labels, eg. [['sky00', 'sky01', 'load00',
    !           'load01'], ['sky10', 'sky11', 'load10', 'load11'], ...]
    ! scan:     derived class (comm_scan)
    !           
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    ! TODO
    ! ----
    ! - ndet, nhorn and detlabels should be taken from tod, not inserted as separate parameters?
    ! 
    implicit none
    class(comm_scan),               intent(inout) :: self
    class(comm_tod),                intent(in)    :: tod
    character(len=*),               intent(in)    :: filename
    integer(i4b),                   intent(in)    :: scan, ndet, nhorn, ndiode
    character(len=*), dimension(:), intent(in)    :: detlabels
    character(len=*), dimension(:,:), intent(in)  :: diode_names

    integer(i4b)       :: i,j,k,l, n, m, ext(1)
    real(sp)           :: nu
    real(dp)           :: scalars(4)
    character(len=6)   :: slabel
    character(len=128) :: field   
    type(hdf_file)     :: file
    real(sp),     allocatable, dimension(:)       :: buffer_sp

    call int2string(scan, slabel)
    call open_hdf_file(filename, file, "r")

    call read_hdf(file, slabel // "/" // "common/ntod",   n)


  end subroutine read_hdf_scan_data


  subroutine read_jumplist(self, jumplist)
    implicit none
    class(comm_tod),   intent(inout) :: self
    character(len=*),  intent(in)    :: jumplist

    integer(i4b) :: i, j, n_jumps, unit


  end subroutine read_jumplist


  subroutine get_scan_ids(self, filelist)
    implicit none
    class(comm_tod),   intent(inout) :: self
    character(len=*),  intent(in)    :: filelist

    integer(i4b)       :: unit, j, k, np, ind(1), i, n, m, n_tot, ierr, p, q, flen, c
    real(dp)           :: w_tot, w_curr, w, v0(3), v(3), spin(2)
    character(len=6)   :: fileid
    character(len=512) :: infile
    real(dp),           allocatable, dimension(:)   :: weight, sid
    real(dp),           allocatable, dimension(:,:) :: spinpos, spinaxis
    integer(i4b),       allocatable, dimension(:)   :: scanid, id, filenum
    integer(i4b),       allocatable, dimension(:)   :: proc
    real(dp),           allocatable, dimension(:)   :: pweight
    character(len=512), allocatable, dimension(:)   :: filename


  end subroutine get_scan_ids


  subroutine dumpToHDF(self, chainfile, iter, map, rms)
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile
    class(comm_map),                   intent(in)    :: map, rms

    integer(i4b)       :: i, j, k, l, npar, ierr
    real(dp)           :: mu
    character(len=6)   :: itext
    character(len=512) :: path
    real(dp), allocatable, dimension(:,:,:) :: output
    real(dp), allocatable, dimension(:)     :: mjds


  end subroutine dumpToHDF

  subroutine initHDF(self, chainfile, iter, map, rms)
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: iter
    type(hdf_file),                    intent(in)    :: chainfile
    class(comm_map),                   intent(inout) :: map, rms

    integer(i4b)       :: i, j, k, npar, ierr, ext(3)
    character(len=6)   :: itext
    character(len=512) :: path
    real(dp), allocatable, dimension(:,:,:) :: output


  end subroutine initHDF

  function get_det_id(self, label)
    implicit none
    class(comm_tod),                intent(inout) :: self
    character(len=*),               intent(in)    :: label
    integer(i4b)                                  :: get_det_id

    integer(i4b) :: i

    do i = 1, self%ndet
       !write(*,*) "From get_det_id:", trim(adjustl(label)), trim(adjustl(self%label(i)))
       !write(*,*) self%ndet
       if (trim(adjustl(label)) == trim(adjustl(self%label(i)))) then
          get_det_id = i
          return
       end if
    end do

    get_det_id = -1

  end function get_det_id

  subroutine initialize_bp_covar(self, filename)
    implicit none
    class(comm_tod),   intent(inout) :: self
    character(len=*),  intent(in)    :: filename

    integer(i4b) :: j, k, ndet, npar, unit, par, ios
    real(dp)     :: val
    character(len=25)   :: label, det1, det2
    character(len=1024) :: line


  end subroutine initialize_bp_covar


  !construct a sidelobe template in the time domain
  subroutine construct_sl_template(self, slconv, pix, psi, s_sl, polangle)
    implicit none
    class(comm_tod),                     intent(in)    :: self
    class(comm_conviqt),                 intent(in)    :: slconv
    integer(i4b),        dimension(:),   intent(in)    :: pix, psi
    real(dp),                            intent(in)    :: polangle
    real(sp),            dimension(:),   intent(out)   :: s_sl

    integer(i4b) :: j,k, pix_, bpsi, psii, psiu, subsamp
    real(dp)     :: psi_, unwrap, x0, x1

    real(dp), dimension(:), allocatable :: sub_sl, x_sl
    type(spline_type) :: spline


    s_sl = 0

  end subroutine construct_sl_template




  subroutine construct_corrtemp_inst(self, scan, pix, psi, s)
    !  Construct an instrument-specific correction template
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
    class(comm_tod),                       intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    s = 0.d0

  end subroutine construct_corrtemp_inst

  
  subroutine construct_dipole_template(self, scan, pix, psi, s_dip)
    !  construct a CMB dipole template in the time domain
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
    !  s_dip: real (sp)
    !       output dipole template timestream
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: scan
    integer(i4b),    dimension(:,:),   intent(in)    :: pix, psi
    real(sp),        dimension(:,:),   intent(out)   :: s_dip

    integer(i4b) :: i, j, ntod
    real(dp)     :: v_ref(3)
    real(dp), allocatable, dimension(:,:) :: P
    logical(lgt)  :: relativistic

    s_dip = 0


  end subroutine construct_dipole_template

  subroutine construct_dipole_template_diff(self, scan, pix, psi, orbital, horn_ind, s_dip, &
                                         &  factor)
    !construct a CMB dipole template in the time domain for differential data
    !
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
    !  orbital: logical
    !       flag for whether the orbital or solar dipole is used as the template
    !  horn: integer
    !       corresponds to either horn = 1 (A) or horn = 2 (B)
    !  Returns:
    !  --------
    !  s_dip: real (sp)
    !       output dipole template timestream
    implicit none
    class(comm_tod),                   intent(inout) :: self
    integer(i4b),                      intent(in)    :: scan
    integer(i4b),    dimension(:,:),   intent(in)    :: pix, psi
    logical(lgt),                      intent(in)    :: orbital
    integer(i4b),                      intent(in)    :: horn_ind
    real(sp),        dimension(:,:),   intent(out)   :: s_dip
    real(dp),               intent(in), optional     :: factor

    integer(i4b) :: i, j, ntod
    real(dp)     :: v_ref(3), v_ref_next(3), f
    real(dp), allocatable, dimension(:,:) :: P
    logical(lgt)  :: relativistic

    s_dip = 0


  end subroutine construct_dipole_template_diff

  subroutine output_scan_list(self, slist)
    implicit none
    class(comm_tod),                               intent(in)    :: self
    character(len=512), allocatable, dimension(:), intent(inout) :: slist

    integer(i4b)     :: i, j, mpistat(MPI_STATUS_SIZE), unit, ns, ierr, num_scan, n_buff
    character(len=4) :: pid


  end subroutine output_scan_list



   subroutine downsample_tod(self, tod_in, ext, tod_out, mask, threshold, width)
      ! Downsamples a time-ordered signal by a moving average filter.
      !
      ! This function is used by calling it twice. In the first call, we providing it with only the 
      ! `tod_in` and `ext` arguments:  `call tod%downsample_tod(tod_in, ext)`. This populates the `ext` 
      ! list with the upper and lower bounds of the downsampled TOD. We then allocate the downsampled 
      ! `tod_out` as `allocate(downsampled_array(ext(1):ext(2)))`. Finally, we call the function again
      ! to get the downsampled tods: `tod%downsample_tod(tod_in, ext, downsampled_array)`.
      ! 
      ! Remember to set `tod_out` to zero before passing it to this function to avoid floating invalid.
      !
      ! Parameters
      ! ----------
      ! tod_in:
      !     The input TOD to be downsampled.
      ! ext:
      !     An integer array of length 2. The first element is the lower bound of the downsampled TOD,
      !     and the second element is the upper bound of the downsampled TOD.
      ! tod_out: optional
      !     The downsampled TOD. If not provided, the function will only populate the `ext` array.
      ! mask: optional
      !     A mask to apply to the TOD before downsampling. If not provided, no mask is applied.
      ! threshold: optional
      !     Sets all values in the downsampled tod to under the threshold to zero.
      ! step: optional
      !     The step size of the downsampled TOD. If not provided, the step size is determined by the
      !     ratio between the input TOD samplerate and the instrument lower resolution samplerate.
      !
      implicit none
      class(comm_tod),                    intent(in)     :: self
      real(sp), dimension(:),             intent(in)     :: tod_in
      integer(i4b),                       intent(inout)  :: ext(2)
      real(sp), dimension(ext(1):ext(2)), intent(out), optional :: tod_out
      real(sp), dimension(:),             intent(in),  optional :: mask
      real(sp),                           intent(in),  optional :: threshold
      real(dp),                           intent(in),  optional :: width

      integer(i4b) :: i, j, k, m, n, ntod, w, npad

   end subroutine downsample_tod



! Fills masked region with linear function between the mean of 20 points at each end
  subroutine fill_masked_region(d_p, mask, i_start, i_end, ntod, chunk)
    implicit none
    real(sp), intent(inout)  :: d_p(:)
    real(sp), intent(in)     :: mask(:)
    integer(i4b), intent(in) :: i_end, i_start, ntod, chunk
    real(dp)     :: mu1, mu2
    integer(i4b) :: i, n_mean, earliest, latest
  end subroutine fill_masked_region


! Identifies and fills masked region
  subroutine fill_all_masked(d_p, mask, ntod, sample, sigma_0, handle, chunk)
    implicit none
    real(sp),         intent(inout)  :: d_p(:)
    real(sp),         intent(in)     :: mask(:)
    real(sp),         intent(in), optional     :: sigma_0
    type(planck_rng), intent(inout), optional  :: handle
    integer(i4b),     intent(in) :: ntod
    logical(lgt),     intent(in) :: sample
    integer(i4b),     intent(in) :: chunk
    integer(i4b) :: j_end, j_start, j, k
    logical(lgt) :: init_masked_region, end_masked_region

    ! Fill gaps in data
    init_masked_region = .true.
    end_masked_region = .false.
    do j = 1,ntod
       if (mask(j) == 1.) then
          if (end_masked_region) then
             j_end = j - 1
             call fill_masked_region(d_p, mask, j_start, j_end, ntod, chunk)
             ! Add noise to masked region
             if (sample) then
                do k = j_start, j_end
                   d_p(k) = d_p(k) + sigma_0 * rand_gauss(handle)
                end do
             end if
             end_masked_region = .false.
             init_masked_region = .true.
          end if
       else
          if (init_masked_region) then
             init_masked_region = .false.
             end_masked_region = .true.
             j_start = j
          end if
       end if
    end do
    ! if the data ends with a masked region
    if (end_masked_region) then
       j_end = ntod
       call fill_masked_region(d_p, mask, j_start, j_end, ntod, chunk)
       if (sample) then
          do k = j_start, j_end
             d_p(k) = d_p(k) + sigma_0 * rand_gauss(handle)
          end do
       end if
    end if

  end subroutine fill_all_masked


  ! Compute chisquare
  subroutine compute_tod_chisq(self, scan, det, mask, s_sky, s_spur, &
       & n_corr, tod, s_jump, absbp, verbose)
    implicit none
    class(comm_tod),                 intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sky, s_spur
    real(sp),          dimension(:), intent(in)     :: n_corr
    real(sp),          dimension(:), intent(in)     :: tod
    real(sp),          dimension(:), intent(in), optional :: s_jump
    logical(lgt),                    intent(in), optional :: absbp, verbose

    
    real(dp)     :: chisq, d0, g
    integer(i4b) :: i, n

    call timer%start(TOD_CHISQ, self%band)

    chisq       = 0.d0
    n           = 0
    g           = self%scans(scan)%d(det)%gain
    do i = 1, self%scans(scan)%ntod
       if (mask(i) < 0.5) cycle
       n     = n+1
       d0    = tod(i) - (g * s_spur(i) + n_corr(i))
       if (present(s_jump)) d0 = d0 - s_jump(i)
       chisq = chisq + (d0 - g * s_sky(i))**2
    end do

    if (self%scans(scan)%d(det)%N_psd%sigma0 <= 0.d0) then
       if (present(absbp)) then
          self%scans(scan)%d(det)%chisq_prop   = 0.d0
       else
          self%scans(scan)%d(det)%chisq        = 0.d0
       end if
    else
       chisq      = chisq      / self%scans(scan)%d(det)%N_psd%sigma0**2
       if (present(absbp)) then
          self%scans(scan)%d(det)%chisq_prop   = chisq
       else
          self%scans(scan)%d(det)%chisq        = (chisq - n) / sqrt(2.d0*n)
       end if
    end if
    if (present(verbose)) then
      if (verbose) write(*,fmt='(a,i8,i3,e16.8,i10,f16.3)') "chi2 :  ", self%scanid(scan), det, &
         &self%scans(scan)%d(det)%N_psd%sigma0, n, self%scans(scan)%d(det)%chisq
    end if

    call timer%stop(TOD_CHISQ, self%band)

  end subroutine compute_tod_chisq

  function get_total_chisq(self, det)
    implicit none
    class(comm_tod),     intent(in)  :: self
    integer(i4b),        intent(in)  :: det
    real(dp)                         :: get_total_chisq

    integer(i4b) :: i, ierr
    real(dp)     :: chisq, buffer

    chisq = 0.d0
    do i = 1, self%nscan ! Sum chisq for all scans
       if (.not. self%scans(i)%d(det)%accept) cycle
       chisq = chisq + self%scans(i)%d(det)%chisq
    end do
    call mpi_reduce(chisq, buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         & 0, self%info%comm, ierr)

    get_total_chisq = buffer

  end function get_total_chisq

  subroutine symmetrize_flags(self, flag)
    implicit none
    class(comm_tod),                          intent(inout) :: self
    integer(i4b),         dimension(1:,1:),   intent(inout) :: flag

    integer(i4b) :: i, det

    do det = 1, self%ndet
       do i = 1, size(flag,1)
          if (iand(flag(i,det),self%flag0) .ne. 0) then
             flag(i,self%partner(det)) = flag(i,det)
          end if
       end do
    end do

  end subroutine symmetrize_flags

  subroutine decompress_pointing_and_flags(self, scan, det, pix, psi, flag)
    implicit none
    class(comm_tod),                    intent(in)            :: self
    integer(i4b),                       intent(in)            :: scan, det
    integer(i4b),        dimension(:),  intent(out), optional :: flag
    integer(i4b),        dimension(:,:),intent(out), optional :: psi, pix
    integer(i4b) :: i, j, k
    do i=1, self%nhorn
       if (present(pix)) then
          call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%pix(i)%p,  pix(:,i))
       end if
       if (present(psi)) then
          call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%psi(i)%p,  psi(:,i))
          if (minval(psi) < 1) then
            write(*,*) 'Psi bin ranges from ', minval(psi), maxval(psi), ', should be 1-indexed'
            stop
          end if
          if (maxval(psi) > self%npsi) then
            write(*,*) 'Psi bin ranges from ', minval(psi), maxval(psi), ', greater than npsi,', self%npsi
            stop
          end if
          if (self%polang(det) /= 0.) then
             do j = 1, size(psi,1)
                psi(j,i) = psi(j,i) + nint(self%polang(det)/(2.d0*pi)*self%npsi)
             end do
          end if
          do j = 1, size(psi,1)
             if (psi(j,i) < 1) then
                psi(j,i) = psi(j,i) + self%npsi
             else if (psi(j,i) > self%npsi) then
                psi(j,i) = psi(j,i) - self%npsi
             end if
          end do
       end if
    end do
    if (present(flag)) then
       call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%flag, flag)

       ! Apply dynamic mask if it exists
       if (allocated(self%scans(scan)%d(det)%mask_dyn)) then
          do i = 1, size(self%scans(scan)%d(det)%mask_dyn,2)
             j = self%scans(scan)%d(det)%mask_dyn(1,i)
             k = self%scans(scan)%d(det)%mask_dyn(2,i)
             flag(j:k) = huge(flag(j))
          end do
       end if
    end if

  end subroutine decompress_pointing_and_flags

  subroutine decompress_diodes(self, scan, det, diodes, flag, pix, psi)
    ! Decompress per-diode tod information
    ! 
    ! Inputs:
    ! ----------
    ! self: comm_tod
    !
    ! scan: integer
    !     scan integer label
    ! det: integer
    !     detector number
    !
    ! Returns:
    ! --------
    ! diodes : real(sp) (ntod, ndiode)
    !    full raw diode values

    implicit none
    class(comm_tod),                    intent(in)  :: self
    integer(i4b),                       intent(in)  :: scan, det
    real(sp),          dimension(:,:),  intent(out) :: diodes
    integer(i4b),      dimension(:),    intent(out), optional :: flag
    integer(i4b),      dimension(:),    intent(out), optional :: pix
    integer(i4b),      dimension(:),    intent(out), optional :: psi

    integer(i4b) :: i, j
    real(sp)     :: tot
    integer(i4b), allocatable, dimension(:) :: buff

!    allocate(buff(size(diodes,1)))
    do i = 1, self%ndiode
!HKEHKE
        call huffman_decode2_sp(self%scans(scan)%todkey, self%scans(scan)%d(det)%zdiode(i)%p, diodes(:,i))
        !tot = sum(diodes(:,i))
        !call huffman_decode3(self%scans(scan)%todkey, self%scans(scan)%d(det)%zdiode(i)%p, buff)
!        write(*,*) sum(abs(diodes(:,i)-buff)), maxval(abs(diodes(:,i)-buff))
        !diodes(:,i) = buff
        !write(*,*) diodes(1:2,i), buff(1:2)
!!$        if (self%myid ==0) then
!!$           open(58,file='test2.dat')
!!$           do j = 1, size(buff)
!!$              write(58,*) diodes(j,i), buff(j)
!!$           end do
!!$           close(58)
!!$        end if
!!$        call mpi_finalize(j)
!!$        stop

    end do
!    deallocate(buff)

    if (present(flag)) then
       call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%flag, flag)
    end if

    if (present(pix)) then ! this assumes nhorn = 1, sorry future person
      call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%pix(1)%p, pix)
    end if

    if (present(psi)) then ! this assumes nhorn = 1, sorry future person
      call huffman_decode2_int(self%scans(scan)%hkey, self%scans(scan)%d(det)%psi(1)%p,  psi, imod=self%npsi-1)
      if (self%polang(det) /= 0.) then
         do j = 1, size(psi)
            psi(j) = psi(j) + nint(self%polang(det)/(2.d0*pi)*self%npsi)
            if (psi(j) < 1) then
               psi(j) = psi(j) + self%npsi
            else if (psi(j) > self%npsi) then
               psi(j) = psi(j) - self%npsi
            end if
         end do
      end if
   end if

  end subroutine decompress_diodes


  subroutine decompress_tod(self, scan, det, tod)
    !
    ! decompresses huffman coded tods 
    !
    ! Argumnets:
    ! ----------
    ! self: comm_tod
    !
    ! scan: integer
    !     scan integer label
    ! det: integer
    !     detector number
    !
    ! Returns:
    ! --------
    ! tod: real(sp)
    !    full raw TOD values
    !
    implicit none
    class(comm_tod),                    intent(in)  :: self
    integer(i4b),                       intent(in)  :: scan, det
    real(sp),            dimension(:),  intent(out) :: tod

    byte,  allocatable, dimension(:)  :: test
    integer(i4b), allocatable, dimension(:) :: tod_int

    if( index(self%instlabel, 'LFI') /= 0) then ! raw data is naturally a float
          call huffman_decode2_sp(self%scans(scan)%todkey, self%scans(scan)%d(det)%ztod, tod)
    else ! raw data is an int
      allocate(tod_int(size(tod)))
      !write(*,*) self%scans(scan)%d(det)%label, self%scans(scan)%chunk_num, size(self%scans(scan)%d(det)%ztod)
      call huffman_decode2_int(self%scans(scan)%todkey, self%scans(scan)%d(det)%ztod, tod_int)
      tod = real(tod_int, sp)
      deallocate(tod_int)
    endif

  end subroutine decompress_tod
  
  subroutine compute_solar_centered_pointing(tod, scan, det, pix, psi, pix_sol, psi_sol)
    implicit none
    class(comm_tod),                  intent(in)  :: tod
    integer(i4b),                     intent(in)  :: scan, det
    integer(i4b),        dimension(:),intent(in)  :: pix, psi
    integer(i4b),        dimension(:),intent(out) :: pix_sol, psi_sol

    integer(i4b) :: i, j
    real(dp)     :: alpha, lat, lon, vec(3), vec0(3), M_sun(3,3), x_sun(3), theta_sun, phi_sun, M_ecl2gal(3,3)

    call ecl_to_gal_rot_mat(M_ecl2gal)
    
    do j = 1, tod%scans(scan)%ntod
       ! Set up transformation from Ecliptic to Solar centered coordinates
       alpha = real(j-1,dp) / real(tod%scans(scan)%ntod-1,dp)
       x_sun = - (1.d0-alpha) * tod%scans(scan)%x0_obs - alpha * tod%scans(scan)%x1_obs ! Sun position at time t in Ecliptic coordinates
       call vec2ang(x_sun, theta_sun, phi_sun)
       call compute_euler_matrix_zyz(-phi_sun, 0.d0, 0.d0, M_sun)
       
       ! Compute pointing in solar centered coordinates
       call pix2vec_ring(tod%nside, pix(j), vec0) ! Galactic coordinates
       vec0 = matmul(transpose(M_ecl2gal), vec0)    ! Ecliptic coordinates
       vec0 = matmul(M_sun, vec0)                   ! Solar-centered coordinates
       call vec2pix_ring(tod%nside, vec0, pix_sol(j))
    end do
    psi_sol = 0 ! Not computed yet
    
  end subroutine compute_solar_centered_pointing


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to save time-ordered-data chunk
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_tod_chunk(filename, tod)
    implicit none
    character(len=*),                   intent(in) :: filename
    real(sp),         dimension(:),     intent(in) :: tod
    ! Expects one-dimensional TOD chunk

    integer(i4b) :: unit, n_tod, t

    n_tod = size(tod)

    unit = getlun()
    open(unit,file=trim(filename), recl=1024)
    write(unit,*) '# TOD value in mK'
    do t = 1, n_tod
       write(unit,fmt='(e16.8)') tod(t)
    end do
    close(unit)
  end subroutine write_tod_chunk


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to save map array to fits file 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_fits_file(filename, array, outmaps)
    implicit none
    character(len=*),                   intent(in) :: filename
    real(dp),         dimension(0:),    intent(in) :: array
    class(map_ptr),   dimension(:),     intent(inout) :: outmaps

    integer(i4b) :: np0, m

    do m = 0, size(array) - 1
       outmaps(1)%p%map(m, 1) = array(m)
    end do

    call outmaps(1)%p%writeFITS(filename)

  end subroutine write_fits_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to save map array to fits file 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_fits_file_iqu(filename, array, outmaps)
    implicit none
    character(len=*),                    intent(in) :: filename
    real(dp),         dimension(0:, 1:), intent(in) :: array
    class(map_ptr),   dimension(:),      intent(inout) :: outmaps

    outmaps(1)%p%map = array

    call outmaps(1)%p%writeFITS(filename)

  end subroutine write_fits_file_iqu


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generic deferred routines that do not do anything
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diode2tod_inst(self, scan, map_sky, procmask, tod)
    ! 
    ! Generates detector-coadded TOD from low-level diode data
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! scan:     int
    !           Scan ID number
    !
    ! Returns
    ! ----------
    ! tod:      ntod x ndet sp array
    !           Output detector TOD generated from raw diode data
    !
    implicit none
    class(comm_tod),                           intent(inout) :: self
    integer(i4b),                              intent(in)    :: scan
    real(sp),          dimension(0:,1:,1:,1:), intent(in)    :: map_sky
    real(sp),          dimension(:),           intent(in)    :: procmask
    real(sp),          dimension(:,:),         intent(out)   :: tod
    tod = 0.
  end subroutine diode2tod_inst

  subroutine read_tod_inst(self, file)
    ! 
    ! Reads instrument-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_LB_tod)
    !           LB-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_tod),                     intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst
  
  subroutine read_scan_inst(self, file, slabel, detlabels, scan)
    ! 
    ! Reads instrument-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
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
    class(comm_tod),                     intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst

  subroutine initHDF_inst(self, chainfile, path)
    ! 
    ! Initializes instrument-specific TOD parameters from existing chain file
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
    class(comm_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_inst

  subroutine load_instrument_inst(self, instfile, band)
    !
    ! Reads the instrument specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_tod
    !    the tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_tod),                     intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band
  end subroutine load_instrument_inst

  
  subroutine dumpToHDF_inst(self, chainfile, path)
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
    class(comm_tod),                     intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_inst

  subroutine remove_fixed_scans(self)
    ! 
    ! Sets accept = .false. for known bad scans
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_tod),                     intent(inout)  :: self
  end subroutine remove_fixed_scans
  
  subroutine apply_map_precond(self, map, map_out)
    implicit none
    class(comm_tod),                   intent(in)    :: self
    real(dp),        dimension(0:,1:), intent(in)    :: map
    real(dp),        dimension(0:,1:), intent(out)   :: map_out

    map_out = map

  end subroutine apply_map_precond

  subroutine collect_v_sun(self)
    implicit none
    class(comm_tod),   intent(inout) :: self

    integer(i4b) :: i, j, ierr

    allocate(self%v_sun(3,self%nscan_tot))
    self%v_sun = 0.d0
    do i = 1, self%nscan
       self%v_sun(:,self%scanid(i)) = self%scans(i)%v_sun
    end do

    call mpi_allreduce(MPI_IN_PLACE, self%v_sun, size(self%v_sun), &
         & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)

  end subroutine collect_v_sun

   subroutine precompute_zodi_lookups(self, cpar)
      class(comm_tod),   intent(inout) :: self
      type(comm_params),       intent(in) :: cpar

      integer(i4b) :: i, j, ierr, pix, pix_high, pix_low, nest_pix, n_subpix, nobs_lowres, npix_lowres, npix_highres
      real(dp), allocatable :: x0_obs(:, :), x1_obs(:, :), x0_earth(:, :), x1_earth(:, :), t0(:), t1(:), ind2vec_zodi_temp(:, :)
      real(dp), allocatable :: x0_obs_packed(:, :), x1_obs_packed(:, :), x0_earth_packed(:, :), x1_earth_packed(:, :), t0_packed(:), t1_packed(:)
      real(dp) :: r, obs_time_end, dt_tod, SECOND_TO_DAY, rotation_matrix(3, 3)
      real(dp), allocatable :: time(:), x_obs(:, :), x_earth(:, :)

      allocate(x0_obs(3, self%nscan_tot), x1_obs(3, self%nscan_tot))
      allocate(x0_earth(3, self%nscan_tot), x1_earth(3, self%nscan_tot))
      allocate(t0(self%nscan_tot), t1(self%nscan_tot))

      x0_obs = 0.
      x1_obs = 0.
      x0_earth = 0.
      x1_earth = 0.
      t0 = 0.
      t1 = 0.

      do i = 1, self%nscan
         t0(self%scanid(i)) = self%scans(i)%t0(1)
         t1(self%scanid(i)) = self%scans(i)%t1(1)
         x0_obs(:, self%scanid(i)) = self%scans(i)%x0_obs
         x1_obs(:, self%scanid(i)) = self%scans(i)%x1_obs
         x0_earth(:, self%scanid(i)) = self%scans(i)%x0_earth
         x1_earth(:, self%scanid(i)) = self%scans(i)%x1_earth
      end do

      
      call mpi_allreduce(MPI_IN_PLACE, t0, size(t0), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call mpi_allreduce(MPI_IN_PLACE, t1, size(t1), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call mpi_allreduce(MPI_IN_PLACE, x0_obs, size(x0_obs), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call mpi_allreduce(MPI_IN_PLACE, x1_obs, size(x1_obs), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call mpi_allreduce(MPI_IN_PLACE, x0_earth, size(x0_earth), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call mpi_allreduce(MPI_IN_PLACE, x1_earth, size(x1_earth), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)

      ! filter out non zero values
      t0_packed = pack(t0, t0 /= 0.)
      t1_packed = pack(t1, t1 /= 0.)
      if (size(t0_packed) /= size(t1_packed)) then
          write(*,*) "Irregularity in number of unique start/end-times of tods, ", size(t0_packed), ' versus ', size(t1_packed), ' needed for zodi/tod interpolation'
          stop
      end if
  

      allocate(x0_obs_packed(3, size(t0_packed)), x1_obs_packed(3, size(t0_packed)))
      allocate(x0_earth_packed(3, size(t0_packed)), x1_earth_packed(3, size(t0_packed)))
      do i = 1, 3
         x0_obs_packed(i, :) = pack(x0_obs(i, :), x0_obs(i, :) /= 0.)
         x1_obs_packed(i, :) = pack(x1_obs(i, :), x1_obs(i, :) /= 0.)
         x0_earth_packed(i, :) = pack(x0_earth(i, :), x0_earth(i, :) /= 0.)
         x1_earth_packed(i, :) = pack(x1_earth(i, :), x1_earth(i, :) /= 0.)
      end do

      ! make new time, obs_pos and earth_pos arrays containing both chunk start and chunk end values
      allocate(time(size(t0_packed) * 2))
      allocate(x_obs(3, size(t0_packed) * 2))
      allocate(x_earth(3, size(t0_packed) * 2))
      do i = 1, size(t0_packed) * 2 , 2
         j = (i - 1) / 2 + 1 ! index from 1, to size(t0)
         time(i) = t0_packed(j)
         time(i + 1) = t1_packed(j)
         x_obs(:, i) = x0_obs_packed(:, j)
         x_obs(:, i + 1) = x1_obs_packed(:, j)
         x_earth(:, i) = x0_earth_packed(:, j)
         x_earth(:, i + 1) = x1_earth_packed(:, j)
      end do

      do i = 2, size(time)
         if (.not. time(i) > time(i - 1)) stop "precomputed MJD time array must be strictly increasing"
      end do

      do i = 1, 3
         call spline_simple(self%x_obs_spline(i), time, x_obs(i, :))
         call spline_simple(self%x_earth_spline(i), time, x_earth(i, :))
      end do

      self%zodi_init_cache_time = self%scans(1)%t0(1)
      call self%clear_zodi_cache()
      
      !allocate spectral quantities
      allocate(self%zodi_b_nu_spl_obj(self%ndet))

      ! allocate cache files and precompute ecliptic unit vectors
      call ecl_to_gal_rot_mat(rotation_matrix)
      allocate(self%zodi_scat_cache(self%nobs, self%zodi_n_comps, self%ndet))
      allocate(self%zodi_therm_cache(self%nobs, self%zodi_n_comps, self%ndet))
      self%zodi_scat_cache = -1.d0
      self%zodi_therm_cache = -1.d0
      allocate(self%ind2vec_ecl(3, self%nobs))
      do i = 1, self%nobs
         self%ind2vec_ecl(:,i) = matmul(self%ind2vec(:,i), rotation_matrix)
      end do
      
      ! If zodi sampling is turned on we precompute lowres zodi lookups
      if (.not. cpar%sample_zodi) return
      ! Skip if zodi nside = tod nside
      if (self%nside == zodi_nside) return
      n_subpix = (self%nside / zodi_nside)**2

      npix_lowres = 12*zodi_nside**2
      npix_highres = 12*self%nside**2

      ! Make lookup table for highres pixels to lowres pixels
      allocate(self%udgrade_pix_zodi(0:npix_highres - 1))
      do i = 0, npix_highres - 1
         call ring2nest(self%nside, i, nest_pix)
         nest_pix = nest_pix / n_subpix
         call nest2ring(zodi_nside, nest_pix, self%udgrade_pix_zodi(i))
      end do
   
   end subroutine precompute_zodi_lookups

   subroutine clear_zodi_cache(self, obs_time)
      ! Clears the zodi tod cache used to look up already computed zodi values. 
      !
      ! This cache has an associate time of observation since the cache is only valid
      ! if the time between the observations are small enough for the observer to not 
      ! have moved significantly.
      class(comm_tod),   intent(inout) :: self
      real(dp), intent(in), optional :: obs_time

      self%zodi_scat_cache = -1.d0
      self%zodi_therm_cache = -1.d0
      if (present(obs_time)) then
         self%zodi_cache_time = obs_time
      else 
         self%zodi_cache_time = self%zodi_init_cache_time
      end if
      if (allocated(self%zodi_therm_cache_lowres)) then
         self%zodi_scat_cache_lowres = -1.d0
         self%zodi_therm_cache_lowres = -1.d0
      end if
   end subroutine clear_zodi_cache

   subroutine create_dynamic_mask(self, scan, det, res, rms_range, mask)
     implicit none
     class(comm_tod),                   intent(inout) :: self
     integer(i4b),                      intent(in)    :: scan, det
     real(sp),            dimension(:), intent(in)    :: res
     real(sp),            dimension(2), intent(in)    :: rms_range
     real(sp),            dimension(:,:), intent(inout) :: mask
     
     logical(lgt) :: cut
     integer(i4b) :: i, n, ntod, nmax
     real(dp) :: box_width, rms, vec(3), elon
     integer(i4b), allocatable, dimension(:,:) :: bad, buffer
     !real(dp),     allocatable, dimension(:,:) :: mask_solar
     
     ntod = size(res)
     nmax = 1000
     
      ! Compute rms
      rms = 0.d0
      n   = 0
      do i = 1, ntod
         if (mask(i,det) /= 1.) cycle
         rms = rms + res(i)**2
         n   = n   + 1
      end do
      if (n <= 1) then
        rms = 0
      else
        rms = sqrt(rms/(n-1))
      end if

!      write(*,*) 'a'
      ! Get full-sky mask
!      if (associated(self%mask_solar)) then
!         allocate(mask_solar(0:12*self%nside**2-1,1))
!         call self%mask_solar%bcast_fullsky_map(mask_solar)
!      end if

      
!      write(*,*) 'b'
      ! Look for strong outliers and masked samples, save bad ranges
      allocate(bad(2,nmax))
      bad = -1
      n   = 0
      do i = 1, ntod

         ! Apply RMS selection criterium
         if (mask(i,det) == 1) then
            cut = res(i) < rms_range(1)*rms .or. res(i) > rms_range(2)*rms
         else
            cut = .false.
         end if
         
         ! Apply solar mask selection criterium
!         call pix2vec_ring(self%nside, self%scans(scan)%d(det)%pix_sol(i,1), vec)
!         elon = acos(min(max(vec(1),-1.d0),1.d0)) * 180.d0/pi                       ! The Sun is at (1,0,0)
         !         cut = cut .or. elon < self%sol_elong_range(1) .or. elon > self%sol_elong_range(2)
         if (allocated(self%mask_solar) .and. self%use_solar_point) then
            cut = cut .or. (self%mask_solar(self%scans(scan)%d(det)%pix_sol(i,1),1) < 0.5)
         end if

         if (cut) then
            ! Start new range if not already active
            if (bad(1,n+1) == -1) bad(1,n+1) = i
            mask(i,det) = 0.
         else
            ! Close active range
            if (bad(1,n+1) /= -1 .and. bad(2,n+1) == -1) then
               bad(2,n+1) = i-1
               n          = n+1
            end if
         end if
         
         ! Increase array size if needed
         if (n == nmax) then
            nmax = 2*nmax
            allocate(buffer(2,nmax))
            buffer = -1
            buffer(:,1:nmax/2) = bad
            deallocate(bad)
            allocate(bad(2,nmax))
            bad = buffer
            deallocate(buffer)
         end if
      end do

      ! Close open range if needed at the end
      if (bad(1,n+1) /= -1 .and. bad(2,n+1) == -1) then
         bad(2,n+1) = ntod
         n          = n+1
      end if
      
      ! Store final array
      if (n > 0) then
         allocate(self%scans(scan)%d(det)%mask_dyn(2,n))
         self%scans(scan)%d(det)%mask_dyn = bad(:,1:n)
!!$         do i = 1, n
!!$            write(*,*) i, bad(:,i)
!!$         end do
         !write(*,fmt='(a,i6,a,i6,i4)') ' Removing ', n, ' ranges in dynamic mask for scan, det', self%scanid(scan), det
      end if

      
!      write(*,*) 'c'
      deallocate(bad)
!      if (allocated(mask_solar)) deallocate(mask_solar)
   end subroutine

   subroutine distribute_sky_maps(tod, map_in, scale, map_out, map_full)
    implicit none
    class(comm_tod),                       intent(in)     :: tod
    type(map_ptr), dimension(1:,1:),       intent(inout)  :: map_in       ! (ndet,ndelta)    
    real(sp),                              intent(in)     :: scale
    real(sp),      dimension(1:,1:,0:,1:), intent(out)    :: map_out
    real(dp),      dimension(1:, 0:), intent(out), optional   :: map_full

    integer(i4b) :: i, j, k, l, npix, nmaps
    real(dp),     allocatable, dimension(:,:) :: m_buf

    npix  = map_in(1,1)%p%info%npix
    nmaps = map_in(1,1)%p%info%nmaps
    allocate(m_buf(0:npix-1,nmaps))
    if (present(map_full)) map_full = 0
    do j = 1, size(map_in,2)       ! ndelta
       do i = 1, size(map_in,1)    ! ndet
          map_in(i,j)%p%map = scale * map_in(i,j)%p%map ! unit conversion
          call map_in(i,j)%p%bcast_fullsky_map(m_buf)
          do k = 1, tod%nobs
             map_out(:,k,i,j) = m_buf(tod%ind2pix(k),:)
          end do
          if (present(map_full) .and. j .eq. 1) then
            do k = 1, nmaps
              do l = 0, npix-1
                map_full(k, l) = map_full(k, l) + m_buf(l, k)
              end do
            end do
          end if
       end do
       do k = 1, tod%nobs
          do l = 1, tod%nmaps
             map_out(l,k,0,j) = sum(map_out(l,k,1:tod%ndet,j))/tod%ndet
          end do
       end do
    end do

    if (present(map_full)) map_full = map_full/tod%ndet

    deallocate(m_buf)

  end subroutine distribute_sky_maps


  subroutine get_s_static(self, band, point, s)
     ! Evaluates the solar centric TOD
     !
     ! Parameters:
     ! -----------
     ! self  : ZodiModel object
     ! band  : Band number id
     ! point : pointing array
     ! s     : Output TOD of static zodi signal
     implicit none
     class(comm_tod),            intent(in)  :: self
     integer(i4b),               intent(in)  :: band
     integer(i4b), dimension(:), intent(in)  :: point
     real(sp),     dimension(:), intent(out) :: s

     integer(i4b) :: i

     if (.not. associated(self%map_solar)) then
        s = 0.d0
        return
     end if
     
     do i = 1, size(s)
        s(i) = self%map_solar(point(i),1)
     end do

   end subroutine get_s_static

  
   
end module comm_tod_mod
