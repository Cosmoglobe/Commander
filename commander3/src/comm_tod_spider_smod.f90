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
submodule (comm_tod_spider_mod) comm_tod_SPIDER_mod
contains
 
   !**************************************************
   !             Constructor
   !**************************************************
   module function constructor(cpar, id_abs, info, tod_type)
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
     type(comm_params),       intent(in) :: cpar
     integer(i4b),            intent(in) :: id_abs
     class(comm_mapinfo),     target     :: info
     character(len=128),      intent(in) :: tod_type
     class(comm_SPIDER_tod),  pointer    :: constructor
 
     integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ierr
     logical(lgt) :: pol_beam
     character(len=8)     :: det, det_partner
     character(len=2)     :: row_str, row_partner_str
     integer(i4b)         :: row_int, row_partner_int
     integer(i4b)         :: det_index(1)


     call timer%start(TOD_INIT, id_abs)

 
     ! Allocate object
     allocate(constructor)
 
     ! Set up noise PSD type and priors
     constructor%freq            = cpar%ds_label(id_abs)
     constructor%n_xi            = 3
     constructor%noise_psd_model = 'oof'
     allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
     allocate(constructor%xi_n_P_rms(constructor%n_xi))
     
     constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0] ! [sigma0, fknee, alpha]; sigma0 is not used
     if (trim(constructor%freq) == 'SPIDER_150') then
        constructor%xi_n_nu_fit = [0.d0, 0.400d0]    ! More than max(2*fknee_DPC) | 350d0
        constructor%xi_n_P_uni(2,:)  = [0.0010d0, 0.45d0] ! fknee
        constructor%xi_n_P_uni(3,:)  = [-2.8d0, -0.4d0]   ! alpha
     else if (trim(constructor%freq) == 'SPIDER_90') then
        constructor%xi_n_nu_fit = [0.d0, 0.400d0]    ! More than max(2*fknee_DPC) | 0.200d0
        constructor%xi_n_P_uni(2,:)  = [0.002d0, 0.40d0]  ! fknee
        constructor%xi_n_P_uni(3,:)  = [-2.8d0, -0.4d0]   ! alpha
     else
        write(*,*) 'Invalid SPIDER frequency label = ', trim(constructor%freq)
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

     if (trim(constructor%freq) == 'SPIDER_150') then
        constructor%chisq_threshold = 100.d0 !20.d0 ! 9.d0
     else if (trim(constructor%freq) == 'SPIDER_90') then
        constructor%chisq_threshold = 20.d0 !20.d0 ! 9.d0
     end if 

     constructor%nmaps           = info%nmaps
   !   constructor%ndet            = num_tokens(trim(cpar%datadir)//'/'//trim(adjustl(cpar%ds_tod_dets(id_abs))), ",")
     
     nside_beam                  = 512
     nmaps_beam                  = 3
     pol_beam                    = .true.
     constructor%nside_beam      = nside_beam
 
     ! Get detector labels
     if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
        call get_detectors(cpar%ds_tod_dets(id_abs), cpar%datadir, constructor%label)
     else
        call get_tokens(trim(adjustl(cpar%ds_tod_dets(id_abs))), ",", constructor%label)
     end if

     ! Define detector partners
   !   do i = 1, constructor%ndet
   !      if (mod(i,2) == 1) then
   !         constructor%partner(i) = i+1
   !      else
   !         constructor%partner(i) = i-1
   !      end if
   !      constructor%horn_id(i) = (i+1)/2
   !   end do


     ! Define detector partners for SPIDER mux layout
     if (trim(constructor%freq) == 'SPIDER_150') then
        ! A/B - pairs are located in alternating rows. E.g., x1r01c09 and x1r02c09 are a pair.
        do i=1, constructor%ndet
           det = constructor%label(i)
           row_str = det(4:5)
           read(row_str, *) row_int
           if (mod(row_int,2)==0) then
              row_partner_int = row_int - 1
           else
              row_partner_int = row_int + 1
           end if
           call int2string(row_partner_int, row_partner_str)
           det_partner = det
           det_partner(4:5) = row_partner_str

           det_index = findloc(constructor%label, det_partner)
           ! If there is no partner, assign it the index of the partnerless detector itself
           if (det_index(1)==0) then
              det_index(1) = i
           end if
           constructor%partner(i) = det_index(1)
        end do

     else if (trim(constructor%freq) == 'SPIDER_90') then
      do i=1, constructor%ndet
         det = constructor%label(i)
         row_str = det(4:5)
         read(row_str, *) row_int
         if (mod(row_int,2)==0) then
            row_partner_int = row_int - 1
         else
            row_partner_int = row_int + 1
         end if
         call int2string(row_partner_int, row_partner_str)
         det_partner = det
         det_partner(4:5) = row_partner_str

         det_index = findloc(constructor%label, det_partner)
         ! If there is no partner, assign it the index of the partnerless detector itself
         if (det_index(1)==0) then
            det_index(1) = i
         end if
         constructor%partner(i) = det_index(1)
      end do
     end if


       
     ! Read the actual TOD
     call constructor%read_tod(constructor%label, cpar%datadir)
 
     ! Initialize bandpass mean and proposal matrix
     call constructor%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))
 
     ! Construct lookup tables
     call constructor%precompute_lookups()
 
     ! Load the instrument file
     call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)
 
     ! Allocate sidelobe convolution data structures
     allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
     if (constructor%orb_4pi_beam) constructor%orb_dp => comm_orbdipole(constructor%mbeam)
 
     ! Initialize all baseline corrections to zero
     do i = 1, constructor%nscan
        constructor%scans(i)%d%baseline = 0.d0
     end do

     call timer%stop(TOD_INIT, id_abs)

 
   end function constructor
 
   !**************************************************
   !             Driver routine
   !**************************************************
   module subroutine process_SPIDER_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
     !
     ! Routine that processes the SPIDER time ordered data.
     ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
     ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
     ! Writes maps to disc in fits format
     !
     ! Arguments:
     ! ----------
     ! self:     pointer of comm_SPIDER_tod class
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
     class(comm_SPIDER_tod),                   intent(inout) :: self
     character(len=*),                         intent(in)    :: chaindir
     integer(i4b),                             intent(in)    :: chain, iter
     type(planck_rng),                         intent(inout) :: handle
     type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
     real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
     class(comm_map),                          intent(inout) :: map_out      ! Combined output map
     class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
     type(map_ptr),     dimension(:,:),   intent(inout), optional :: map_gain
 

     real(dp)            :: t1, t2
     integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
     logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
     type(comm_binmap)   :: binmap
     type(comm_scandata) :: sd
     character(len=4)    :: ctext, myid_text
     character(len=6)    :: samptext, scantext
     character(len=512)  :: prefix, postfix, prefix4D, filename
     character(len=512), allocatable, dimension(:) :: slist
     real(sp), allocatable, dimension(:)       :: procmask, procmask2, sigma0
     real(sp), allocatable, dimension(:,:)     :: s_buf
     real(sp), allocatable, dimension(:,:,:)   :: d_calib
     real(sp), allocatable, dimension(:,:,:,:) :: map_sky
     real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

     real(sp),     allocatable, dimension(:,:)   :: s_jump, tod_gapfill
     real(sp),     allocatable, dimension(:,:,:) :: jump_calib
     integer(i4b), allocatable, dimension(:,:)   :: jumps, offset_range, jumpflag_range
     real(sp),     allocatable, dimension(:)     :: offset_level
     type(comm_binmap)                           :: jump_map
     character(len=4)                            :: it_label
     logical(lgt)                                :: debug
     real(sp),    allocatable, dimension(:)      :: test_array




     call int2string(iter, ctext)
     call int2string(iter, it_label)
     call update_status(status, "tod_start"//ctext)
     call timer%start(TOD_TOT, self%band)


     ! Toggle optional operations
     sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
     sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
     select_data           = (iter>4) !self%first_call        ! only perform data selection the first time
     output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
 
     ! Initialize local variables
     ndelta          = size(delta,3)
     self%n_bp_prop  = ndelta
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

     ! Distribute maps
     allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
     call distribute_sky_maps(self, map_in, 1.e-6, map_sky) ! uK to K

     ! Distribute processing masks
     allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
     call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
     call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
     deallocate(m_buf)

     ! Precompute far sidelobe Conviqt structures
     if (self%correct_sl) then
        call timer%start(TOD_SL_PRE, self%band)

        if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
        do i = 1, self%ndet
           !TODO: figure out why this is rotated
           call map_in(i,1)%p%YtW()  ! Compute sky a_lms
           self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
                & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
                & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
        end do
        call timer%stop(TOD_SL_PRE, self%band)
     end if
 !    write(*,*) 'qqq', self%myid
 !    if (.true. .or. self%myid == 78) write(*,*) 'a', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
 !!$    call mpi_finalize(ierr)
 !!$    stop
 
     call update_status(status, "tod_init")

     !------------------------------------
     ! Perform main sampling steps
     !------------------------------------
     ! Sample gain components in separate TOD loops; marginal with respect to n_corr
   !   call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
   !   call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
   !   call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)
 


     ! Prepare intermediate data structures
     call binmap%init(self, .true., sample_rel_bandpass)
     call jump_map%init(self, .true., sample_rel_bandpass)  
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

        ! Prepare data
        if (sample_rel_bandpass) then
 !          if (.true. .or. self%myid == 78) write(*,*) 'b', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
           call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
        else if (sample_abs_bandpass) then
           call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_sky_prop=.true.)
        else
           call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
        end if
        allocate(s_buf(sd%ntod,sd%ndet))

        ! Calling Simulation Routine
        if (self%enable_tod_simulations) then
           call simulate_tod(self, i, sd%s_tot, sd%n_corr, handle)
           call dealloc_scandata(sd)
           cycle
        end if









        
















        ! REMOVE JUMPS ----------------------------------
        allocate(s_jump(sd%ntod, sd%ndet))
        allocate(jumps(sd%ntod, sd%ndet))
        allocate(tod_gapfill(sd%ntod, sd%ndet))
        allocate(test_array(sd%ntod))

        debug = .false.



        do j=1, sd%ndet
         ! Throw away detectors that are more than 80% flagged. Also throw away detectors that don't have a partner. 
         if ((sum(sd%flag(:,j)) > 0.8*sd%ntod) .or. (.not. self%scans(i)%d(j)%accept) .or. (j==self%partner(j))) then
            self%scans(i)%d(j)%accept = .false.
            cycle
         end if

         
           
           ! Retrieve offsets from previous run, if they exist
           if (allocated(self%scans(i)%d(j)%offset_range)) then
              call expand_offset_list(              &
                 & self%scans(i)%d(j)%offset_range, &
                 & self%scans(i)%d(j)%offset_level, &
                 & s_jump(:,j))
           else
              s_jump(:,j) = 0
           end if

           
           ! Retrieve jump flags from previous run, if they exist
           if (allocated(self%scans(i)%d(j)%jumpflag_range)) then
            call add_jumpflags(                     &
            & self%scans(i)%d(j)%jumpflag_range, &
            & sd%flag(:,j))
         end if

         
         ! Scanning for jumps


         if (.true.) then
            call jump_scan(                                 &
            & sd%tod(:,j) - sd%s_sky(:,j) - s_jump(:,j), &
            & sd%flag(:,j),                              &
            & jumps(:,j),                                &
            & offset_range,                              &
            & offset_level,                              &
            & handle,                                    &
            & jumpflag_range,                            &
            & it_label,                                  &
            & chaindir,                                  &
            & debug)
            


              ! Add offsets to persistent list
              if (.not. allocated(self%scans(i)%d(j)%offset_range)) then
                 allocate(self%scans(i)%d(j)%offset_range(size(offset_level),2))
                 allocate(self%scans(i)%d(j)%offset_level(size(offset_level)))

                 self%scans(i)%d(j)%offset_range = offset_range
                 self%scans(i)%d(j)%offset_level = offset_level
              else
                 call update_offset_list(              &
                    & offset_range,                    &
                    & offset_level,                    &
                    & self%scans(i)%d(j)%offset_range, &
                    & self%scans(i)%d(j)%offset_level)
              end if

              ! Add jump flags to persistent list
              if (allocated(jumpflag_range)) then
                 if (.not. allocated(self%scans(i)%d(j)%jumpflag_range)) then
                    allocate(self%scans(i)%d(j)%jumpflag_range(size(jumpflag_range)/2,2))
                    self%scans(i)%d(j)%jumpflag_range = jumpflag_range
                 else
                    call update_jumpflag(jumpflag_range, self%scans(i)%d(j)%jumpflag_range)
                 end if
              end if

              call expand_offset_list(                &
                  & self%scans(i)%d(j)%offset_range,  &
                  & self%scans(i)%d(j)%offset_level,  & 
                  & s_jump(:,j))
           end if


           call gap_fill_linear(           &
              & sd%tod(:,j) - s_jump(:,j), &
              & sd%flag(:,j),              &
              & tod_gapfill(:,j),          &
              & handle,                    &
              & .true.)


           if (allocated(offset_range))   deallocate(offset_range)
           if (allocated(offset_level))   deallocate(offset_level)
           if (allocated(jumpflag_range)) deallocate(jumpflag_range)

           if (debug) then
              call tod2file(trim(adjustl(chaindir))//'/tod_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt',         sd%tod(:,j))
              call tod2file(trim(adjustl(chaindir))//'/flag_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt',        sd%flag(:,j))
              call tod2file(trim(adjustl(chaindir))//'/s_sky_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt',       sd%s_sky(:,j))
              call tod2file(trim(adjustl(chaindir))//'/s_jump_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt',      s_jump(:,j))
              call tod2file(trim(adjustl(chaindir))//'/tod_gapfill_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt', tod_gapfill(:,j))
              call tod2file(trim(adjustl(chaindir))//'/s_tot_'//trim(adjustl(self%label(j)))//'_'//trim(adjustl(it_label))//'.txt',       sd%s_tot(:,j))
           end if
        end do











































        ! Sample correlated noise
        if (self%first_call) then
           do j=1, sd%ndet
              self%scans(i)%d(j)%N_psd%xi_n(1) = 0.0018
           end do
        end if
        do j=1, sd%ndet
           self%scans(i)%d(j)%gain = 1.d0
        end do

      !   call sample_n_corr(self, tod_gapfill, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)
        call sample_n_corr(self, tod_gapfill, handle, i, 1.0-sd%flag, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.) 
        

        if (debug) then
           do j=1, sd%ndet
              call tod2file(trim(adjustl(chaindir))//'/n_corr_'//trim(adjustl(it_label))//'.txt', sd%n_corr(:,j))
              call tod2file(trim(adjustl(chaindir))//'/mask_'//trim(adjustl(it_label))//'.txt', sd%mask(:,j))
           end do
        end if

        ! Compute noise spectrum parameters
      !   call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)
      !   call sample_noise_psd(self, tod_gapfill, handle, i, sd%mask, sd%s_sky, sd%n_corr)
        call sample_noise_psd(self, tod_gapfill, handle, i, 1.0-sd%flag, sd%s_sky, sd%n_corr)


        ! Compute chisquare
        call timer%start(TOD_CHISQ, self%band)
        do j = 1, sd%ndet
           if (.not. self%scans(i)%d(j)%accept) cycle
         !   call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j))
           call self%compute_chisq(i, j, 1.0-sd%flag(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j), s_jump=s_jump(:,j))
        end do
        call timer%stop(TOD_CHISQ, self%band)



        ! Select data
        if (select_data) call remove_bad_data(self, i, sd%flag)
 
        ! Compute chisquare for bandpass fit
        if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)
 
        ! Compute binned map
        allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
        call compute_calibrated_data(self, i, sd, d_calib, jump_template=s_jump)


        allocate(jump_calib(1, sd%ntod, sd%ndet))
        jump_calib(1,:,:) = s_jump 

        ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
        if (self%output_4D_map > 0) then
           if (mod(iter-1,self%output_4D_map) == 0) then
              call timer%start(TOD_4D, self%band)
              allocate(sigma0(sd%ndet))
              do j = 1, sd%ndet
                 sigma0(j) = self%scans(i)%d(j)%N_psd%sigma0/self%scans(i)%d(j)%gain
              end do
              call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
                   & samptext, self%scanid(i), self%nside, self%npsi, &
                   & self%label, self%horn_id, real(self%polang*180/pi,sp), sigma0, &
                   & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
                   & self%scans(i)%d(:)%accept)
              deallocate(sigma0)
              call timer%stop(TOD_4D, self%band)
           end if
        end if

        ! Bin TOD
        call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib,    binmap)
        call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, jump_calib, jump_map) 

         


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
        call dealloc_scandata(sd)
        deallocate(s_buf, d_calib)
        deallocate(s_jump, jumps, tod_gapfill, jump_calib, test_array)

     end do

     if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'
 
     ! Output latest scan list with new timing information
     if (output_scanlist) call self%output_scan_list(slist)
 
     ! Solve for maps
     call synchronize_binmap(binmap, self)
     call synchronize_binmap(jump_map, self) 

     if (sample_rel_bandpass) then
        call finalize_binned_map(self, binmap, rms_out, 1.d6, chisq_S=chisq_S, mask=procmask2)
     else
        call finalize_binned_map(self, binmap, rms_out, 1.d6)
     end if
     map_out%map = binmap%outmaps(1)%p%map

     call finalize_binned_map(self, jump_map, rms_out, 1.d6) 

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
     if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
     if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
     if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
     if (self%output_n_maps > 6) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
     call jump_map%outmaps(1)%p%writeFITS(trim(prefix)//'jumps'//trim(postfix)) 

     ! Clean up
     call binmap%dealloc()
     call jump_map%dealloc() 
     if (allocated(slist)) deallocate(slist)
     deallocate(map_sky, procmask, procmask2)
     if (self%correct_sl) then
        do i = 1, self%ndet
           call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
        end do
     end if
 
     ! Parameter to check if this is first time routine has been called
     self%first_call = .false.
 
     call update_status(status, "tod_end"//ctext)
     call timer%stop(TOD_TOT, self%band)
 
   end subroutine process_SPIDER_tod
 
   
   module subroutine read_tod_inst_SPIDER(self, file)
     ! 
     ! Reads SPIDER-specific common fields from TOD fileset
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
     ! file:     derived type (hdf_file)
     !           Already open HDF file handle; only root includes this
     !
     ! Returns
     ! ----------
     ! None, but updates self
     !
     implicit none
     class(comm_SPIDER_tod),                 intent(inout)          :: self
     type(hdf_file),                      intent(in),   optional :: file
   end subroutine read_tod_inst_SPIDER
   
   module subroutine read_scan_inst_SPIDER(self, file, slabel, detlabels, scan)
     ! 
     ! Reads SPIDER-specific scan information from TOD fileset
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(in)    :: self
     type(hdf_file),                      intent(in)    :: file
     character(len=*),                    intent(in)    :: slabel
     character(len=*), dimension(:),      intent(in)    :: detlabels
     class(comm_scan),                    intent(inout) :: scan
   end subroutine read_scan_inst_SPIDER
 
   module subroutine initHDF_SPIDER(self, chainfile, path)
     ! 
     ! Initializes SPIDER-specific TOD parameters from existing chain file
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(inout)  :: self
     type(hdf_file),                      intent(in)     :: chainfile
     character(len=*),                    intent(in)     :: path
   end subroutine initHDF_SPIDER
   
   module subroutine dumpToHDF_SPIDER(self, chainfile, path)
     ! 
     ! Writes SPIDER-specific TOD parameters to existing chain file
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(in)     :: self
     type(hdf_file),                      intent(in)     :: chainfile
     character(len=*),                    intent(in)     :: path
   end subroutine dumpToHDF_SPIDER

   module subroutine write2file(filename, iter, param)
      implicit none
   
      character(len=*), intent(in)         :: filename
      real(dp), intent(in)                 :: param
      integer(i4b), intent(in)             :: iter
   
      integer(i4b)                           :: unit, io_error
      logical                                :: existing
   
      unit = 22
   
      inquire(file=trim(filename),exist=existing)
      if (existing) then
         open(unit,file=trim(filename),status='old',position='append',action='write',iostat=io_error)
      else
         open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
      end if
   
      write(unit,*), iter, param
   
      close(unit)
    end subroutine write2file
 
 end submodule comm_tod_SPIDER_mod
 
