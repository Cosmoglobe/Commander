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
module comm_tod_LB_mod
  !   Module which contains all the LiteBIRD time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_LB_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !
  use comm_tod_mod
  use comm_tod_driver_mod
  use comm_tod_simulations_mod
  use comm_conviqt_mod
  use comm_tod_mapmaking_mod
  use comm_hdf_mod
  implicit none

  !private
  !public comm_LB_tod

  type, extends(comm_tod) :: comm_LB_tod
   contains
     procedure     :: process_tod          => process_LB_tod
  end type comm_LB_tod

  interface comm_LB_tod
     procedure constructor_lb
  end interface comm_LB_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor_lb(cpar, id, id_abs, info, tod_type) result(c)
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
    class(comm_LB_tod),      pointer    :: c

    integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ierr
    logical(lgt) :: pol_beam

    ! Allocate object
    allocate(c)

    ! Set up noise PSD type and priors
    c%freq            = cpar%ds_label(id_abs)
    c%n_xi            = 3
    c%noise_psd_model = 'oof'
    allocate(c%xi_n_P_uni(c%n_xi,2))
    allocate(c%xi_n_nu_fit(c%n_xi,2))
    allocate(c%xi_n_P_rms(c%n_xi))

    ! Correlated noise parameters
    c%xi_n_nu_fit(1,:) = [1.d0, 8.5d0]   ! Freq range for sigma0 in Hz
    c%xi_n_nu_fit(2,:) = [1d-3, 1.d0]    ! Freq range for fknee in Hz
    c%xi_n_nu_fit(3,:) = [1d-3, 1.d0]    ! Freq range for alpha in Hz
      
    c%xi_n_P_uni(1,:)  = [1d-6, 1d-3]     ! Uniform prior for sigma0
    c%xi_n_P_uni(2,:)  = [1d-6, 1.d0]     ! Uniform prior for fknee
    c%xi_n_P_uni(3,:)  = [-4d0, -0.5d0]   ! Uniform prior for alpha
      
    ! Set rms of all parameters to 0.05 for initial test phase. 
    c%xi_n_P_rms(1)    = 1.00d0           ! Prior rms for sigma0
    c%xi_n_P_rms(2)    = 0.1d0            ! Prior rms for fknee
    c%xi_n_P_rms(3)    = 0.2d0            ! Prior rms for alpha

    ! Initialize common parameters
    call c%tod_constructor(cpar, id, id_abs, info, tod_type)

    ! Initialize instrument-specific parameters
    c%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    c%nhorn           = 1
    c%ndiode          = 1
    c%compressed_tod  = .false.
    c%correct_sl      = .false.
    c%correct_orb     = .true.
    c%orb_4pi_beam    = .false.
    c%symm_flags      = .true.
    c%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
    c%nmaps           = info%nmaps
    if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
       c%ndet         = count_detectors(cpar%ds_tod_dets(id_abs)) !, cpar%datadir)
    else
       c%ndet         = num_tokens(cpar%ds_tod_dets(id_abs), ",")
    end if
    nside_beam                  = 512
    nmaps_beam                  = 3
    pol_beam                    = .true.
    c%nside_beam      = nside_beam

    ! Get detector labels
    if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
        call get_detectors(cpar%ds_tod_dets(id_abs), c%label)
    else
        call get_tokens(trim(adjustl(cpar%ds_tod_dets(id_abs))), ",", c%label)
    end if
        
    ! Define detector partners
    do i = 1, c%ndet
       if (mod(i,2) == 1) then
          c%partner(i) = i+1
       else
          c%partner(i) = i-1
       end if
       c%horn_id(i) = (i+1)/2
    end do

    ! Read the actual TOD
    call c%read_tod(c%label)
    ! option to set noise parameters by hand for on-the-fly litebird sims
    if (c%on_the_fly_tod_sim .and. cpar%sim_noisepar) call overwrite_noisepar(c, cpar)
    
    ! Initialize bandpass mean and proposal matrix
    call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup table
    c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false.)
    call c%precompute_lookups()
    
    ! Load the instrument file
    call c%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Allocate sidelobe convolution data structures
    allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)
    !c%orb_dp => comm_orbdipole(c%mbeam)
    c%orb_dp => comm_orbdipole(comm=info%comm)

  end function constructor_lb

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LB_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    ! 
    ! Routine that processes the LiteBIRD time ordered data. 
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
    ! Writes maps to disc in fits format
    ! 
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_LB_tod class
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
    class(comm_LB_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
    type(map_ptr),       dimension(1:),       intent(inout), optional :: map_gain       ! (ndet)
    
    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, h, ierr, ndelta, nside, npix, nmaps, oper_default
    logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist, sample_ncorr, sample_xi_n
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist

    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    type(hdf_file) :: tod_file
    
    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)
    call timer%start(TOD_TOT, self%band) 

    ! Output input sky model
    call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
    
    ! Toggle optional operations
    sample_ncorr          = .true. !.true. OBS
    sample_xi_n           = .true. !.false. ! OBS
    sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    select_data           = self%first_call        ! only perform data selection the first time
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
    sample_gain           = .false.                ! Gain sampling, LB TOD sims have perfect gain

    ! Define useful sd operation codes
    if (sample_ncorr) then
       oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,&
            & SD_TOD,SD_SKY, SD_NCORR])
            !& SD_TOD,SD_SKY,SD_SL,SD_ORB,SD_INST,SD_ZODI, SD_NCORR])
            ! OBS FIXME include V_SUN, sidelobes and zodi in litebird files!
    else
       oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,&
            & SD_TOD,SD_SKY])
            !& SD_TOD,SD_SKY,SD_SL,SD_ORB,SD_INST,SD_ZODI])
            ! OBS FIXME include V_SUN, sidelobes and zodi in litebird files!
    end if
    
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

    ! Initialize index-based sky map and mask
    ! OBS: input sky maps are in uK as usual, while LiteBIRD tods are in K
    ! Therefore we are scaling sky maps to K, using the optional scale parameter
    call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain, scale=1e-6)  

    ! Precompute far sidelobe Conviqt structures
    if (self%correct_sl) then
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i,1)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if
    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
     if (sample_gain) then
       ! 'abscal': the global constant gain factor
       call sample_calibration(self, 'abscal', oper_default, handle)
       ! 'relcal': the gain factor that is constant in time but varying between detectors
       call sample_calibration(self, 'relcal', oper_default, handle)
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
       !write(*,*) "Scan number", i
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
       allocate(s_buf(sd%ntod,sd%ndet))

       ! call simulate_tod if we want to make and output tod sims
       if (self%enable_tod_simulations) then
          call simulate_tod(self, i, sd%s_tot(:,:,0,1), sd%n_corr, handle)
       end if

       ! call simulate_tod_on_the_fly if we want to make tod sim in memory and continue analyzing it
       ! but only if this is the first sample
       if (self%on_the_fly_tod_sim .and. self%first_call) then
          call simulate_tod_on_the_fly(self, sd, handle)
          ! we have now overwritten tod in self, and also have to update uncompresses data in sd
          call dealloc_scan_data(sd)
          call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
       end if

       ! Sample correlated noise
       if (sample_ncorr) then
          call sample_n_corr(self, sd, handle)
          if (sample_xi_n) then
             call sample_noise_psd(self, sd, handle, chaindir)
          else
             call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true.)
          end if
       else
          call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true.)
       end if
       
       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_tod_chisq(sd, j)
       end do

       ! Select data
       !if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
       
       ! Compute binned map
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)    

       ! For debugging: write TOD to hdf
       !if (.true.) then
       if (self%on_the_fly_tod_sim .and. self%first_call) then
          ! scan id appears to be the worst chi2
          if (self%scanid(i) == 1) then 
             !print *, self%scanid(i)
             call int2string(self%scanid(i), scantext)
             call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//'_'//scantext//'.h5', tod_file, 'w')
             call write_hdf(tod_file, '/tod', sd%tod)
             call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
             call write_hdf(tod_file, '/flag', sd%flag)
             call write_hdf(tod_file, '/caltod', d_calib(1, :, :))
             call write_hdf(tod_file, '/s_sky', sd%s_sky)
             !call write_hdf(tod_file, '/s_inst', sd%s_inst)
             call write_hdf(tod_file, '/n_corr', sd%n_corr)
             !call write_hdf(tod_file, '/s_sl', sd%s_sl)
             !call write_hdf(tod_file, '/s_orb', sd%s_orb)
             call write_hdf(tod_file, '/res', d_calib(2, :, :))
             !call write_hdf(tod_file, '/zodi', d_calib(7, :, :))
             call write_hdf(tod_file, '/mask', sd%mask)
             !call write_hdf(tod_file, '/accept', real(self%scans(i)%d(:)%accept,dp))
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
       deallocate(s_buf, d_calib)

    end do

    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

    call update_status(status, "finalizing maps, BP")

    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    call synchronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       call finalize_binned_map(self, binmap, rms_out, 1.d6, chisq_S=chisq_S)
    else
       call finalize_binned_map(self, binmap, rms_out, 1.d6)
    end if
    map_out%map = binmap%outmaps(1)%p%map

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, handle, chisq_S, delta)
       self%bp_delta = delta(:,:,1)
    end if

    ! Output maps to disk
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
    !if (self%output_n_maps > 3) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'hitmap'//trim(postfix))
    if (self%output_n_maps > 4) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 7) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

    ! Clean up
    call binmap%dealloc()
    if (allocated(slist)) deallocate(slist)
    if (self%correct_sl) then
       do i = 1, self%ndet
          do h = 1, self%nhorn
             call self%slconv(i,h)%p%dealloc(); deallocate(self%slconv(i,h)%p)
          end do
       end do
    end if

    ! Parameter to check if this is first time routine has been called
    self%first_call = .false.

    call update_status(status, "tod_end"//ctext)
    call timer%stop(TOD_TOT, self%band) 
  end subroutine process_LB_tod   

  subroutine overwrite_noisepar(self, cpar)
    ! setting noise parameters sigma0, fknee and alpha by hand for on-the-fly litebird sims
    ! hardcoded sigma0 values for 11-13 specific litebird channels
    ! overwriting existing values in memory for self%scans(k)%d(j)%N_psd%xi_n
    implicit none
    class(comm_LB_tod),   intent(inout) :: self
    type(comm_params),    intent(in)    :: cpar
    real(sp), allocatable, dimension(:) :: sigma0
    real(dp)                            :: root_nsamp_per_arcmin
    integer(i4b)                        :: i, j, k, numband, nband

    if (self%nscan == 0) return

    !if (self%myid==0) write(*,*) 'band ', self%band, trim(self%freq) 
    !if (self%myid==0) write(*,*) 'sigma0 before ', self%scans(1)%d(1)%N_psd%sigma0 
    !if (self%myid==0) write(*,*) 'fknee  before ', self%scans(1)%d(1)%N_psd%xi_n(2) 
    !if (self%myid==0) write(*,*) 'alpha  before ', self%scans(1)%d(1)%N_psd%xi_n(3)

    numband = count(cpar%ds_active)
    allocate(sigma0(numband))

    ! table of sigma0 values to pick from for different litebird configurations
    ! unit uK_cmb*root(arcmin ^2)  (not uK_cmb*sqrt(s))
    ! sensitivity for Q and U separately
    !         L1-40  L2-50  L1-61  L2-77  M1-94  M2-118 M1-145 M2-182 H1-217 H2-280 H1-334 H2-402 H3-570
    !sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59, 11.27, 14.93, 33.62,   0.0] !postKDP2
    !sigma0 = [35.96, 25.24, 18.90, 12.87,  7.20,  4.71,  4.41,  3.40,  8.59, 11.27, 14.93, 33.62, 233.1] !lessLF2
    !sigma0 = [35.96, 20.61, 18.90, 10.51,  8.05,  4.71,  4.94,  3.40,  8.59, 11.27, 14.93, 33.62, 233.1] !lessMF1 
    !sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  5.15,  4.41,  3.73,  8.59, 11.27, 14.93, 33.62, 233.1] !lessMF2
    !sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59,  0.00, 14.93,  0.00, 233.1] !noHF2
    !sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59,  0.00, 14.93, 29.13, 233.1] !newHF2

    if (trim(cpar%noisepar_ver) == 'postKDP2') then
       nband = 12
       sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59, 11.27, 14.93, 33.62]
    else if (trim(cpar%noisepar_ver) == 'lessLF2') then
       nband = 13
       sigma0 = [35.96, 25.24, 18.90, 12.87,  7.20,  4.71,  4.41,  3.40,  8.59, 11.27, 14.93, 33.62, 233.1]
    else if (trim(cpar%noisepar_ver) == 'lessMF1') then
       nband = 13
       sigma0 = [35.96, 20.61, 18.90, 10.51,  8.05,  4.71,  4.94,  3.40,  8.59, 11.27, 14.93, 33.62, 233.1]
    else if (trim(cpar%noisepar_ver) == 'lessMF2') then
       nband = 13
       sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  5.15,  4.41,  3.73,  8.59, 11.27, 14.93, 33.62, 233.1]
    else if (trim(cpar%noisepar_ver) == 'noHF2') then
       nband = 11
       sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59, 14.93, 233.1]
    else if (trim(cpar%noisepar_ver) == 'newHF2') then
       nband = 12
       sigma0 = [35.96, 20.61, 18.90, 10.51,  7.20,  4.71,  4.41,  3.40,  8.59, 14.93, 29.13, 233.1]
    else if (trim(cpar%noisepar_ver) == 'debug570') then
       nband = 1
       sigma0 = [233.1]
    end if

    if (numband /= nband) then
       write(*,*) trim(cpar%noisepar_ver), ' assumes', nband ,' bands, not', numband
    end if

    ! num arcmin_square on a sphere: 41252.96 degree_square * 3600 arcmin_square/degree_square
    ! numsamp for a frequency band: 65536 samples per scan * 9142 number of scans * numdets per band
    ! (all litebird sim scans has the same number of samples (ntod))
    ! nsamp_per_arcmin_square = 65536*9142*4/41252.96/3600.d0
    root_nsamp_per_arcmin = sqrt( self%scans(1)%ntod * real(self%nscan_tot * self%ndet,dp)/41252.96/3600.d0 )
 
    !write(*,*) 'root_nsamp_per_arcmin', root_nsamp_per_arcmin 
    !write(*,*) 'root_nsamp_per_arcmin', sqrt(65536*9142/41252.96/3600.d0*4)
    !write(*,*) 'ntod, nscan_tot, ndet', self%scans(1)%ntod, self%nscan_tot, self%ndet

    
    ! loop over scans and set new noise parameters
    do k = 1, self%nscan
       do j = 1, self%ndet
          ! want sigm0 (aka xi_n(1)) in K (litebird tods are in K), while table above is in uK*arcmin
          ! given sigma0 is for Q and U so the total sensitivity is sqrt(2) higher
          self%scans(k)%d(j)%N_psd%sigma0  = sigma0(self%band) * root_nsamp_per_arcmin * 1e-6 /sqrt(2.d0)
          self%scans(k)%d(j)%N_psd%xi_n(2) = 0.05                 ! fknee = 50 mHz
          self%scans(k)%d(j)%N_psd%xi_n(3) = -1                   ! alpha
       end do
    end do
    deallocate(sigma0)

    if (self%myid==0) then
       write(*,*) '|> Put sigma0 =', self%scans(1)%d(1)%N_psd%sigma0
       write(*,*) '|       fknee =', self%scans(1)%d(1)%N_psd%xi_n(2)
    end if
       
  end subroutine overwrite_noisepar
    
end module comm_tod_LB_mod
