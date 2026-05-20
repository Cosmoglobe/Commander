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
      integer(i4b)  :: tsys_order
      real(dp)      :: tsys_eta0
      real(dp), allocatable, dimension(:,:)    :: tsys_fit ! ndet, tsys_order+1
   contains
      procedure     :: process_tod          => process_chipass_tod
      procedure     :: read_scan_inst          => read_scan_inst_chipass
      procedure     :: construct_corrtemp_inst => construct_corrtemp_chipass
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
      c%xi_n_P_uni(2,:)  = [0.001d0, 0.0011d0]  ! fknee
      c%xi_n_P_uni(3,:)  = [-1.5d0,   -1.4d0]  ! alpha

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

      ! allocate CHIPASS-specific instrument file data
      c%tsys_order      = 3
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


      ! Initialize all baseline corrections to zero
      do i = 1, c%nscan
         do j = 1, c%ndet
            c%scans(i)%d(j)%baseline = 0.d0
         end do 
      end do

      ! Reject obviously bad scans
      do i = 1, c%nscan
         do j = 1, c%ndet
            if (c%scans(i)%d(j)%N_psd%sigma0 == 0.d0) c%scans(i)%d(j)%accept = .false.
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
      logical(lgt)        :: select_data, sample_gain, output_scanlist, sample_ncorr, sample_baseline, sample_tsys
      type(comm_binmap)   :: binmap
      type(comm_scandata) :: sd
      character(len=4)    :: ctext, myid_text
      !character(len=2)    :: zodi_param_text
      !character(len=1)    :: up_down_text
      character(len=6)    :: samptext, scantext
      character(len=512)  :: prefix, postfix!, prefix4D, prefix_atlas, postfix_atlas
      character(len=512), allocatable, dimension(:) :: slist
      !real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi
      real(sp), allocatable, dimension(:,:,:)   :: d_calib
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
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0             ! only output scanlist every 10th iteration
      sample_baseline       = iter > 1
      sample_tsys           = iter > 1
      sample_gain           = iter > 1                         ! Gain sampling
      sample_ncorr          = iter > 1
         
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

      if (self%myid == 0) write(*,*) 'baseline init scan 1 det 1:', self%scans(1)%d(1)%baseline
      if (self%myid == 0) write(*,*) 'tsys_fit init det 1:', self%tsys_fit(1, :)

      ! TODO: sample tsys_eta0

      ! sample Tsys
      if (sample_tsys) then
         if (self%myid == 0) then
            write(*,*) '|    --> Sampling Tsys'
         end if
         call update_status(status, "Tsys")
         call sample_chipass_Tsys(self, oper_default)
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
            call sample_chipass_baseline(self, i, sd%tod, sd%s_tot(:,:,0,1), sd%mask)
            call timer%stop(TOD_BASELINE, self%band)
            call dealloc_scan_data(sd)
         end do
         call timer%start(TOD_WAIT, self%band)
         call mpi_barrier(self%comm, ierr)
         call timer%stop(TOD_WAIT, self%band)
      end if

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

      ! Prepare intermediate data structures
      call binmap%init(self, .true., .false.)

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
         
         ! Compute binned map
         allocate(d_calib(self%output_n_maps, sd%ntod, sd%ndet))
         d_calib = 0.d0
         call compute_calibrated_data(self, i, sd, d_calib)

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
      call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0)
      map_out%map = binmap%outmaps(1)%p%map

      ! Output maps to disk
      call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
      call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
      if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
      if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
      ! added
      !if (self%myid == 0) write(*,*) '   [comm_tod_chipass_mod] output_n_maps:', self%output_n_maps
      if (self%output_n_maps > 7) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'baseline'//trim(postfix))

      ! Clean up
      call binmap%dealloc()
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
             do i = 0, tod%tsys_order
                y(n) = y(n) - tod%tsys_fit(j, i) * (tod%scans(scan)%d(j)%elev(k) - tod%tsys_eta0)**i
             end do
          end if
       end do

       if (n > tod%baseline_order+1) call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline)
    end do

    deallocate(x, y)

  end subroutine sample_chipass_baseline

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

  end subroutine read_scan_inst_chipass

  subroutine sample_chipass_Tsys(self, oper_default)
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
    integer(i4b), intent(in)               :: oper_default

    type(hdf_file)      :: file
    !character(len=6)    :: slabel
    type(comm_scandata) :: sd
    integer(i4b)        :: i, j, k, l, n
    real(dp)            :: t, dt
    real(dp), allocatable, dimension(:) :: x, y !, elev

    !allocate(x(50690*150), y(50690*150))
    n = 0
    do i = 1, self%nscan
        n = n + self%scans(i)%ntod
    end do
    allocate(x(n), y(n))
    if (self%myid == 0) write(*,*) 'ntod tot:', n

    do j = 1, self%ndet
       n = 0
       do i = 1, self%nscan
          if (.not. any(self%scans(i)%d%accept)) cycle
          call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, spur_level=0)
          dt = 1.d0 / self%scans(i)%ntod
          t  = 0.d0
          do k = 1, self%scans(i)%ntod
             t = t + dt
             if (sd%flag(k, j) < 1) then
                n    = n + 1
                x(n) = self%scans(i)%d(j)%elev(k) - self%tsys_eta0
                y(n) = sd%tod(k,j) - self%scans(i)%d(j)%gain * sd%s_tot(k,j,0,1)
                do l = 0, self%baseline_order
                   y(n) = y(n) - self%scans(i)%d(j)%baseline(l) * t**l
                end do
             end if
          end do
          call dealloc_scan_data(sd)
       end do
       call fit_polynomial(x(1:n), y(1:n), self%tsys_fit(j, 0:self%tsys_order))
    end do

    deallocate(x, y)

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
    integer(i4b)     :: i, j, k, d, scan
    real(dp)         :: t, dt
    !real(dp), allocatable, dimension(:) :: elev

    scan      = sd%scan
    dt        = 1.d0 / self%scans(scan)%ntod
    sd%s_inst = 0.d0
    do j = 1, self%ndet
       d = j; if (present(det)) d = det
       t = 0.d0
       if (.not. self%scans(scan)%d(d)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          do i = 0, self%baseline_order
             sd%s_inst(k,j) = sd%s_inst(k,j) + self%scans(scan)%d(j)%baseline(i) * t**i
          end do
          do i = 0, self%tsys_order
             sd%s_inst(k,j) = sd%s_inst(k,j) + self%tsys_fit(j, i) * (self%scans(scan)%d(j)%elev(k) - self%tsys_eta0)**i
          end do
       end do
    end do
  end subroutine construct_corrtemp_chipass

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

    if (self%myid == 0 .and. trim(self%level) == 'L1') then
       call write_hdf(chainfile, trim(adjustl(path))//'tsys_eta0', self%tsys_eta0)
       call write_hdf(chainfile, trim(adjustl(path))//'tsys_fit', self%tsys_fit)
    end if

  end subroutine dumpToHDF_chipass

end module comm_tod_chipass_mod
