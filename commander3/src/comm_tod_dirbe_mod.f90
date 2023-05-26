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

   use comm_tod_mod
   use comm_param_mod
   use comm_map_mod
   use comm_conviqt_mod
   use pix_tools
   use healpix_types
   use comm_huffman_mod
   use comm_hdf_mod
   use comm_fft_mod
   use comm_shared_arr_mod
   use spline_1D_mod
   use comm_4D_map_mod
   use comm_tod_driver_mod
   use comm_utils
   use comm_bp_mod
   use comm_tod_zodi_mod

   implicit none

   private
   public comm_dirbe_tod

   type, extends(comm_tod) :: comm_dirbe_tod
   contains
      procedure     :: process_tod          => process_DIRBE_tod
   end type comm_dirbe_tod

   interface comm_dirbe_tod
      procedure constructor
   end interface comm_dirbe_tod

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
      integer(i4b),            intent(in) :: id_abs        !index of the current band within the parameters 
      class(comm_mapinfo),     target     :: info
      character(len=128),      intent(in) :: tod_type      !
      class(comm_dirbe_tod),      pointer    :: constructor

      integer(i4b) :: i, j, nside_beam, lmax_beam, nmaps_beam, ierr, SMALL_CHUNK_THRESHOLD
      logical(lgt) :: pol_beam

      call timer%start(TOD_INIT, id_abs)

      ! Allocate object
      allocate(constructor)

      ! Set up noise PSD type and priors
      constructor%freq            = cpar%ds_label(id_abs)
      constructor%n_xi            = 3
      constructor%noise_psd_model = 'oof'
      allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
      allocate(constructor%xi_n_P_rms(constructor%n_xi))
      
      constructor%xi_n_P_rms      = [-1.0, 0.1, -0.2] ! [sigma0, fknee, alpha]; sigma0 is not used
      if (.true.) then
         ! constructor%xi_n_nu_fit     = [0.,    0.200] ! More than max(2*fknee_default)
         constructor%xi_n_P_uni(2,:) = [0.001, 0.45]  ! fknee
         constructor%xi_n_P_uni(3,:) = [-2.5, -0.4]   ! alpha
         else
         write(*,*) 'Invalid DIRBE frequency label = ', trim(constructor%freq)
         stop
      end if

      ! Initialize common parameters
      call constructor%tod_constructor(cpar, id_abs, info, tod_type)

      ! Initialize instrument-specific parameters
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn           = 1
      constructor%ndiode          = 1
      constructor%compressed_tod  = .false.
      constructor%correct_sl      = .false.
      constructor%orb_4pi_beam    = .false.
      constructor%symm_flags      = .false.
      constructor%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
      constructor%nmaps           = info%nmaps
      constructor%ndet            = num_tokens(trim(cpar%ds_tod_dets(id_abs)), ",")

      nside_beam                  = 128
      nmaps_beam                  = 1
      pol_beam                    = .false.
      constructor%nside_beam      = nside_beam
      SMALL_CHUNK_THRESHOLD = 1000

      ! Get detector labels
      call get_tokens(trim(cpar%ds_tod_dets(id_abs)), ",", constructor%label)

      ! Read the actual TOD
      call constructor%read_tod(constructor%label)

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

      if (constructor%subtract_zodi) call initialize_tod_zodi_mod(cpar, constructor)

      ! Construct lookup tables
      call constructor%precompute_lookups()

      ! Load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

   ! commenting this out for now
      !do i=1, constructor%ndet
      !  call init_noise_model(constructor, i)
      !end do

      ! Allocate sidelobe convolution data structures
      allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
      constructor%orb_dp => comm_orbdipole(constructor%mbeam)
      ! Initialize all baseline corrections to zero
      !do i = 1, constructor%nscan
      !   constructor%scans(i)%d%baseline = 0.d0
      !end do

      ! do i = 1, constructor%nscan
      !    do j = 1, constructor%ndet
      !       if (constructor%scans(i)%ntod <= SMALL_CHUNK_THRESHOLD) constructor%scans(i)%d(j)%accept = .false.
      !    end do
      ! end do
      call timer%stop(TOD_INIT, id_abs)
   end function constructor

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
      type(map_ptr),       dimension(1:,1:),    intent(inout), optional :: map_gain       ! (ndet,1)
      real(dp)            :: t1, t2
      integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps, N_LINEAR_SAMP_COMPS
      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist, sample_zodi, output_zodi_comps
      type(comm_binmap)   :: binmap
      type(comm_scandata) :: sd
      character(len=4)    :: ctext, myid_text
      character(len=6)    :: samptext, scantext
      character(len=512)  :: prefix, postfix, prefix4D
      character(len=512), allocatable, dimension(:) :: slist
      real(sp), allocatable, dimension(:)       :: procmask, procmask2
      real(sp), allocatable, dimension(:,:)     :: s_buf
      real(sp), allocatable, dimension(:,:,:)   :: d_calib
      real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
      real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)    :: A_T_A_emiss, A_T_A_albedo, A_T_A_reduced
      real(dp), allocatable, dimension(:)       :: AY_emiss , AY_albedo, AY_reduced, X
      type(hdf_file) :: tod_file
      character(len=10), allocatable, dimension(:) :: zodi_comp_names

      allocate(zodi_comp_names(zodi%n_comps))
      zodi_comp_names = [character(len=10) :: "cloud", "band1", "band2", "band3", "ring", "feature"]

      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)
      call timer%start(TOD_TOT, self%band)

      ! Toggle optional operations
      sample_zodi           = .true. ! Sample zodi parameters
      output_zodi_comps     = .true. ! Output zodi components
      sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
      select_data           = .false. !self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
      sample_gain           = .false.                ! Gain sampling, LB TOD sims have perfect gain

      ! Initialize local variables
      ndelta          = size(delta,3)
      self%n_bp_prop  = ndelta-1
      nside           = map_out%info%nside
      nmaps           = map_out%info%nmaps
      npix            = 12*nside**2
      self%output_n_maps = 8
      if (output_zodi_comps) self%output_n_maps = self%output_n_maps + zodi%n_comps
      if (self%output_aux_maps > 0) then
         if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 8
      end if

      ! Sampling parameters
      N_LINEAR_SAMP_COMPS = zodi%n_comps
      if (sample_zodi) then
         allocate(A_T_A_emiss(N_LINEAR_SAMP_COMPS, N_LINEAR_SAMP_COMPS))
         allocate(A_T_A_albedo, mold=A_T_A_emiss)
         allocate(AY_emiss(N_LINEAR_SAMP_COMPS))
         allocate(AY_albedo, mold=AY_emiss)
         allocate(A_T_A_reduced, mold=A_T_A_emiss)
         allocate(AY_reduced, mold=AY_emiss)
         allocate(X, mold=AY_emiss)
         A_T_A_emiss = 0.d0
         A_T_A_albedo = 0.d0
         A_T_A_reduced = 0.d0
         AY_emiss = 0.d0
         AY_albedo = 0.d0
         AY_reduced = 0.d0
         X = 0.d0
      end if

      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'

      ! Distribute maps
      allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
      call distribute_sky_maps(self, map_in, 1.e0, map_sky) ! uK to K
      allocate(m_gain(nmaps,self%nobs,0:self%ndet,1))
      call distribute_sky_maps(self, map_gain, 1.e0, m_gain) ! uK to K

      allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
      call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
      call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
      deallocate(m_buf)

      call update_status(status, "tod_init")

      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------

      ! Sample gain components in separate TOD loops; marginal with respect to n_corr
      if (sample_gain) then
         ! 'abscal': the global constant gain factor
         call sample_calibration(self, 'abscal', handle, map_sky, m_gain, procmask, procmask2)
         ! 'relcal': the gain factor that is constant in time but varying between detectors
         ! call sample_calibration(self, 'relcal', handle, map_sky, m_gain, procmask, procmask2)
         ! 'deltaG': the time-variable and detector-variable gain
         call sample_calibration(self, 'deltaG', handle, map_sky, m_gain, procmask, procmask2)
      end if

      ! Sample dynamic (before computing zodi)
      if (sample_zodi) call sample_dynamic_zodi_parameters(self)

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

         ! Prepare data and generate zodi timestreams
         call sd%init_singlehorn(self, i, map_sky, m_gain, procmask, procmask2, init_s_bp=.true.)

         !  sd%tod = sd%tod
         allocate(s_buf(sd%ntod,sd%ndet))

         ! Sample correlated noise
         !  call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)
         sd%n_corr = 0.d0

         ! Compute noise spectrum parameters
         ! call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
         end do

         ! Select data
         if (select_data) call remove_bad_data(self, i, sd%flag)

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
         
         ! Sample linear zodi parameters (after computing zodi)
         if (sample_zodi) then
            call accumulate_zodi_emissivities(self, sd, handle, A_T_A_emiss, AY_emiss)
            call accumulate_zodi_albedos(self, sd, handle, A_T_A_albedo, AY_albedo)
         end if

         ! Compute binned map
         allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
         d_calib = 0.d0
         call compute_calibrated_data(self, i, sd, d_calib)    

         
         if (.false.) then
            call int2string(self%scanid(i), scantext)
            if (self%myid == 0 .and. i == 1) write(*,*) '| Writing tod to hdf'
            call open_hdf_file(trim(chaindir)//'/tod_'//scantext//'_samp'//samptext//'.h5', tod_file, 'w')
            ! call write_hdf(tod_file, '/n_corr', sd%n_corr)
            !call write_hdf(tod_file, '/bpcorr', sd%s_bp)
            !call write_hdf(tod_file, '/s_tot', sd%s_tot)
            !call write_hdf(tod_file, '/s_sky', sd%s_sky)
            ! call write_hdf(tod_file, '/tod',   sd%tod)
            ! call write_hdf(tod_file, '/flag', sd%flag)
            ! call write_hdf(tod_file, '/pix', sd%pix)
            call write_hdf(tod_file, '/mask1', procmask)
            call write_hdf(tod_file, '/mask2', procmask2)
            !call write_hdf(tod_file, '/pixA', sd%pix(:,1,1))
            !call write_hdf(tod_file, '/pixB', sd%pix(:,1,2))
            !call write_hdf(tod_file, '/psiA', sd%psi(:,1,1))
            !call write_hdf(tod_file, '/psiB', sd%psi(:,1,2))
            !call write_hdf(tod_file, '/x_im', self%x_im)

            !do k = 1, self%ndet
            !  call int2string(k, scantext)
            !  call write_hdf(tod_file, '/xi_n_'//scantext, self%scans(i)%d(k)%N_psd%xi_n)
            !  call write_hdf(tod_file, '/gain_'//scantext, self%scans(i)%d(k)%gain)
            !end do

            call close_hdf_file(tod_file)
            print *, "printed mask.", trim(chaindir)//'/tod_'//scantext//'_samp'//samptext//'.h5'
            stop
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
         call sd%dealloc
         deallocate(s_buf, d_calib)

      end do

      ! Reduce the A_T_A and AY matrices from each tod object and solve the normal equations to find emissivity
      if (sample_zodi) then
         ! Sample emissivities
         call mpi_reduce(A_T_A_emiss, A_T_A_reduced, size(A_T_A_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call mpi_reduce(AY_emiss, AY_reduced, size(AY_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         if (self%myid == 0) then
            call sample_linear_parameter(A_T_A_reduced, AY_reduced, handle, X)
            print *, "Sampled emissivity: ", X
            self%zodi_emissivity = X
         end if
         call mpi_bcast(self%zodi_emissivity, size(self%zodi_emissivity), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

         ! Sample albedo
         A_T_A_reduced = 0.d0
         AY_reduced = 0.d0
         X = 0.d0
         call mpi_reduce(A_T_A_albedo, A_T_A_reduced, size(A_T_A_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call mpi_reduce(AY_albedo, AY_reduced, size(AY_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         if (self%myid == 0) then
            call sample_linear_parameter(A_T_A_reduced, AY_reduced, handle, X)
            print *, "Sampled albedo: ", X
            self%zodi_albedo = X
         end if
         call mpi_bcast(self%zodi_albedo, size(self%zodi_albedo), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

      end if
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
         call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
         self%bp_delta = delta(:,:,1)
      end if


      ! Output maps to disk
      call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
      call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
      if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
      if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
      
      !  if (self%output_n_maps > 2) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'hitmap'//trim(postfix))
      !  if (self%output_n_maps > 4) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
      !  if (self%output_n_maps > 5) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
      !  if (self%output_n_maps > 6) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
      if (self%output_n_maps > 7 .and. self%subtract_zodi) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
      if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
         do i = 1, zodi%n_comps
            call binmap%outmaps(7+i)%p%writeFITS(trim(prefix)//'zodi_'//trim(zodi_comp_names(i))//trim(postfix))
         end do
      endif

      ! Clean up
      call binmap%dealloc()
      if (allocated(slist)) deallocate(slist)
      deallocate(map_sky, procmask, procmask2)
      !  if (self%correct_sl) then
      !     do i = 1, self%ndet
      !        call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
      !     end do
      !  end if

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

      call update_status(status, "tod_end"//ctext)
      
      call timer%stop(TOD_TOT, self%band)
   end subroutine process_DIRBE_tod   


end module comm_tod_DIRBE_mod
