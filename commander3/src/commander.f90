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
program commander
  use comm_param_mod
  use comm_data_mod
  use comm_signal_mod
  use comm_cr_mod
  use comm_chisq_mod
  use comm_output_mod
  use comm_comp_mod
  use comm_nonlin_mod
  use comm_tod_simulations_mod
  use comm_tod_gain_mod
  use comm_zodi_mod
  use comm_zodi_samp_mod
  use comm_tod_zodi_mod
  implicit none

  integer(i4b)        :: i, j, iargc, ierr, iter, stat, first_sample, samp_group, curr_samp, tod_freq
  real(dp)            :: t0, t1, t2, t3, dbp
  logical(lgt)        :: ok, first
  type(comm_params)   :: cpar
  type(planck_rng)    :: handle, handle_noise

  ! param_vec for zodi (REMOVE THIS AFTER ATLAS)
   ! type(zodi_model) :: current_model, previous_model
   ! real(dp), allocatable :: param_vec(:)
   
  type(comm_mapinfo), pointer :: info => null()
  type(comm_map),     pointer :: m    => null()
  class(comm_comp),   pointer :: c1   => null()

  !----------------------------------------------------------------------------------
  ! Command line arguments
  character(len=*), parameter :: version = '1.0.0'
  character(len=32)           :: arg
  integer                     :: arg_indx

  ! Giving the simple command line arguments for user to chose from.
  comm3_args: do arg_indx = 1, command_argument_count()
    call get_command_argument(arg_indx, arg)

    select case (arg)
      case ('-v', '--version')
        print '(2a)', 'Commander3 version: '//trim(version)
        print '(2a)', "Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo."
        print '(2a)', "This is free software; see the source for copying conditions. There is NO warranty;"
        print '(2a)', "not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
        call exit(0)
      case ('-h', '--help')
        call print_help()
        call exit(0)
      case default
        !print '(2a, /)', 'Unrecognised command-line option: ', arg
        !call print_help()
        !call exit(0)
        exit comm3_args
    end select
  end do comm3_args
  !----------------------------------------------------------------------------------

  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************
  call wall_time(t0)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, cpar%myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, cpar%numprocs, ierr)
  
  cpar%root = 0
  
  if (cpar%myid == cpar%root) call wall_time(t1)
  call read_comm_params(cpar)
  if (cpar%myid == cpar%root) call wall_time(t3)
  
  call initialize_mpi_struct(cpar, handle, handle_noise)
  call validate_params(cpar)  
  call init_status(status, trim(cpar%outdir)//'/comm_status.txt', cpar%numband, cpar%comm_chain)
  status%active = cpar%myid_chain == 0 !.false.
  call timer%start(TOT_RUNTIME); call timer%start(TOT_INIT)

!!$  n = 100000
!!$  q = 100000
!!$  allocate(arr(n))
!!$  do i = 1, n
!!$     allocate(arr(i)%p(q))
!!$     arr(i)%p = i
!!$     if (mod(i,1000) == 0) then
!!$        write(*,*) 'up', arr(i)%p(6)
!!$        call update_status(status, "debug")
!!$     end if
!!$  end do
!!$
!!$  do i = 1, n
!!$     deallocate(arr(i)%p)
!!$     if (mod(i,1000) == 0) then
!!$        write(*,*) 'down', i
!!$        call update_status(status, "debug2")
!!$     end if
!!$  end do
!!$  deallocate(arr)
!!$  stop
  
  if (iargc() == 0) then
     if (cpar%myid == cpar%root) write(*,*) 'Usage: commander [parfile] {sample restart}'
     call MPI_Finalize(ierr)
     stop
  end if
  if (cpar%myid == cpar%root) call wall_time(t2)

  ! Output a little information to notify the user that something is happening
  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) then
     write(*,fmt='(a)') ' ---------------------------------------------------------------------'
     write(*,fmt='(a)') ' |                             Commander3                            |'
     write(*,fmt='(a)') ' ---------------------------------------------------------------------'
     if (cpar%enable_tod_simulations) then
       write(*,fmt='(a,t70,a)')       ' |  Regime:                            TOD Simulations', '|'
     else
       write(*,fmt='(a,t70,a)')       ' |  Regime:                            Data Processing', '|'
     endif
     write(*,fmt='(a,i12,t70,a)') ' |  Number of chains                       = ', cpar%numchain, '|'
     write(*,fmt='(a,i12,t70,a)') ' |  Number of processors in first chain    = ', cpar%numprocs_chain, '|'
     write(*,fmt='(a,t70,a)')         ' |', '|'
!     write(*,fmt='(a,f12.3,a,t70,a)') ' |  Time to initialize run                 = ', t2-t0, ' sec', '|'
!     write(*,fmt='(a,f12.3,a,t70,a)') ' |  Time to read in parameters             = ', t3-t1, ' sec', '|'
     write(*,fmt='(a)') ' ---------------------------------------------------------------------'
  end if

  ! ************************************************
  ! *               Initialize modules             *
  ! ************************************************

  if (cpar%myid == cpar%root) call wall_time(t1)

  call update_status(status, "init")

  ! Initialize tod modules
  if (cpar%enable_tod_analysis) then 
    call initialize_tod_mod(cpar)
    if (cpar%include_tod_zodi) then 
      call initialize_zodi_mod(cpar)
      call initialize_tod_zodi_mod(cpar)
    end if
  end if

  call define_cg_samp_groups(cpar)
  ! Initialising Bandpass
  ! TODO: Add QUIET stuff into bandpass module
  call initialize_bp_mod(cpar);             call update_status(status, "init_bp")
  ! Initialising Data -- load it into memory?
  call initialize_data_mod(cpar, handle);   call update_status(status, "init_data")
  ! Debug statement to actually see whether
  ! QUIET is loaded into memory

  ! Precompute zodi lookups
  if (cpar%enable_tod_analysis .and. cpar%include_tod_zodi) then
     do i = 1, numband
          if (data(i)%tod_type == 'none') cycle
          if (.not. data(i)%tod%subtract_zodi) cycle        
          call data(i)%tod%precompute_zodi_lookups(cpar)
          call read_tod_zodi_params(cpar, zodi_model, data(i)%tod)
     end do
  end if

  ! initialize zodi samp mod
  if (cpar%sample_zodi .and. cpar%include_tod_zodi) call initialize_zodi_samp_mod(cpar)
!   call zodi_model_to_ascii(cpar, zodi_model, "/mn/stornext/u3/metins/dirbe/chains/chains_testing/init_zodi.dat")
!   call mpi_barrier(MPI_COMM_WORLD, ierr)
!   stop
! Example use case of zodi_to_ascii for HKE
!   call zodi_model_to_ascii(cpar, zodi_model, "/mn/stornext/u3/metins/dirbe/chains/chains_testing/init_zodi.dat")
<<<<<<< HEAD
!   call ascii_to_zodi_model(cpar, zodi_model, "/mn/stornext/d5/data/duncanwa/DIRBE/ascii_powell.dat")
!   call ascii_to_zodi_model(cpar, zodi_model, "/mn/stornext/u3/metins/dirbe/chains/chains_testing/ascii_K98.dat")
=======
!   call ascii_to_zodi_model(cpar, zodi_model, "/mn/stornext/u3/metins/dirbe/data/ascii_K98.dat")
>>>>>>> 2dca5de640afb7c5ae0a05579f75b9ba9a646116

  call initialize_signal_mod(cpar);         call update_status(status, "init_signal")
  call initialize_from_chain(cpar, handle, first_call=.true.); call update_status(status, "init_from_chain")


!write(*,*) 'Setting gain to 1'
!data(6)%gain = 1.d0

  ! Make sure TOD and BP modules agree on initial bandpass parameters
  ok = trim(cpar%cs_init_inst_hdf) /= 'none'
  if (ok) ok = trim(cpar%init_chain_prefix) /= 'none'
  if (cpar%enable_tod_analysis) call synchronize_bp_delta(ok)
  call update_mixing_matrices(update_F_int=.true.)       

  if (cpar%output_input_model) then
     if (cpar%myid == 0) write(*,*) 'Outputting input model to sample number 999999'
     call output_FITS_sample(cpar, 999999, .false.)
     call mpi_finalize(ierr)
     stop
  end if

  ! Output SEDs for each component
  if (cpar%output_debug_seds) then
     if (cpar%myid == cpar%root) call dump_components('sed.dat')
     call mpi_finalize(ierr)
     stop
  end if
  
  if (cpar%myid == cpar%root) call wall_time(t2)
  
  ! **************************************************************
  ! *                   Carry out computations                   *
  ! **************************************************************

  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) then 
     !write(*,*) '|'
     write(*,fmt='(a)') ' ---------------------------------------------------------------------'
     write(*,fmt='(a,f12.3,a)') ' |  Time to read data = ', t2-t1, ' sec'
     write(*,*) '|  Starting Gibbs sampling'
  end if


  ! Prepare chains 
  call init_chain_file(cpar, first_sample)
  !first_sample = 1
  if (first_sample == -1) then
     call output_FITS_sample(cpar, 0, .true.)  ! Output initial point to sample 0
     first_sample = 1
  else
     ! Re-initialise seeds and reinitialize
     call initialize_mpi_struct(cpar, handle, handle_noise, reinit_rng=first_sample)
     !first_sample = 10
     first_sample=first_sample-1 ! Reject last sample, which may be corrupt
     call initialize_from_chain(cpar, handle, init_samp=first_sample, init_from_output=.true., first_call=.true.)
     first_sample = first_sample+1
  end if
  call timer%stop(TOT_INIT)
  !data(1)%bp(0)%p%delta(1) = data(1)%bp(0)%p%delta(1) + 0.2
  !data(2)%bp(0)%p%delta(1) = data(1)%bp(0)%p%delta(1) + 0.2


  ! Run Gibbs loop
  iter  = first_sample
  first = .true.
  !----------------------------------------------------------------------------------
  ! Part of Simulation routine
  !----------------------------------------------------------------------------------
  ! Will make only one full gibbs loop to produce simulations
  !if (cpar%enable_tod_simulations) cpar%num_gibbs_iter = 2
  !----------------------------------------------------------------------------------
  do while (iter <= cpar%num_gibbs_iter)
     ok = .true.

     call timer%start(TOT_GIBBSSAMP)
     if (cpar%myid_chain == 0) then
        call wall_time(t1)
        write(*,fmt='(a)') ' ---------------------------------------------------------------------'
        write(*,fmt='(a,i4,a,i8)') ' |  Chain = ', cpar%mychain, ' -- Iteration = ', iter
     end if
     ! Initialize on existing sample if RESAMP_CMB = .true.
     if (cpar%resamp_CMB) then
        if (mod(iter-1,cpar%numsamp_per_resamp) == 0 .or. iter == first_sample) then
           curr_samp = mod((iter-1)/cpar%numsamp_per_resamp,cpar%last_samp_resamp-cpar%first_samp_resamp+1) + cpar%first_samp_resamp
           if (cpar%myid_chain == 0) write(*,*) '|  Re-initializing on sample ', curr_samp
           call initialize_from_chain(cpar, handle, init_samp=curr_samp)
           call update_mixing_matrices(update_F_int=.true.)       
        end if
     end if
     !----------------------------------------------------------------------------------
     ! Part of Simulation routine
     !----------------------------------------------------------------------------------
     ! If we are on 1st iteration and simulation was enabled,
     ! we copy real LFI data into specified folder.
     if ((iter == 1) .and. cpar%enable_tod_simulations) then
       call copy_tod(cpar, ierr)
       call write_filelists_to_disk(cpar, ierr)
     end if
     !----------------------------------------------------------------------------------
     ! Process TOD structures
     if (iter > 1 .and. cpar%enable_TOD_analysis .and. (iter <= 2 .or. mod(iter,cpar%tod_freq) == 0)) then

      ! Create zodi atlas
      ! if (.false.) then
      !    allocate(param_vec(base_zodi_model%N_PARAMETERS))
      !    do i = 1, base_zodi_model%N_PARAMETERS
      !       if (cpar%myid == cpar%root) print *, "Creating zodi atlas for parameter: ", i
      !       do j = 1, 3
      !          base_zodi_model = sampled_zodi_model
      !          base_zodi_model%param_i = i
      !          base_zodi_model%up_down_j = j
               
      !          param_vec = base_zodi_model%model_to_param_vec()
      !          if (j == 1) then
      !             param_vec(i) = param_vec(i) - 0.1*param_vec(i)
      !          else if (j == 3) then
      !             param_vec(i) = param_vec(i) + 0.1*param_vec(i)
      !          end if
      !          call base_zodi_model%param_vec_to_model(param_vec)
      !          call process_TOD(cpar, cpar%mychain, iter, handle)
      !       end do
      !    end do
      ! end if 
      
   !   if (iter == 1) then ! For faster component separation since we dont sample the cios

        ! First iteration should just be component separation, in case sky model
        ! is off
        call timer%start(TOT_TODPROC)
        call process_TOD(cpar, cpar%mychain, iter, handle)
        call timer%stop(TOT_TODPROC)
     end if


     if (cpar%enable_tod_simulations) then
        ! Skip other steps if TOD simulations
        exit
     end if

   if (iter > 1 .and. cpar%enable_TOD_analysis .and. cpar%sample_zodi) then
      call timer%start(TOT_ZODI_SAMP)
      if (iter == 2) call downsamp_invariant_structs(cpar) !downsample and cache tod and pointing
      call project_and_downsamp_sky(cpar) ! project skymodel down to downsampled pointing
      call sample_linear_zodi(cpar, handle, iter, zodi_model, verbose=.true.)
      call minimize_zodi_with_powell(cpar)
      ! call sample_zodi_group(cpar, handle, iter, zodi_model, verbose=.true.)
      call timer%stop(TOT_ZODI_SAMP)
   end if

     ! Sample non-linear parameters
     if (iter > 1 .and. cpar%sample_specind) then
        call timer%start(TOT_SPECIND)
        call sample_nonlin_params(cpar, iter, handle, handle_noise)
        call timer%stop(TOT_SPECIND)
     end if
     !if (mod(iter,cpar%thinning) == 0) call output_FITS_sample(cpar, 100+iter, .true.)

     ! Sample linear parameters with CG search; loop over CG sample groups
     !call output_FITS_sample(cpar, 1000+iter, .true.)
     if (cpar%sample_signal_amplitudes) then
        call timer%start(TOT_AMPSAMP)
        do samp_group = 1, cpar%cg_num_user_samp_groups
           if (cpar%myid_chain == 0) then
              write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
                   & samp_group, ' of ', cpar%cg_num_user_samp_groups
           end if
           call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

           if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

        end do
        call timer%stop(TOT_AMPSAMP)

        ! Perform joint alm-Cl Metropolis move
        call timer%start(TOT_CLS)
        do i = 1, 3
           if (cpar%resamp_CMB .and. cpar%sample_powspec) call sample_joint_alm_Cl(handle)
        end do
        call timer%stop(TOT_CLS)
     end if

     ! Sample power spectra
     call timer%start(TOT_CLS)
     if (cpar%sample_powspec) call sample_powspec(handle, ok)
     call timer%stop(TOT_CLS)

     ! Output sample to disk
     call timer%start(TOT_OUTPUT)
     if (mod(iter,cpar%thinning) == 0) call output_FITS_sample(cpar, iter, .true.)
     call timer%stop(TOT_OUTPUT)

     ! Sample partial-sky templates
     !call sample_partialsky_tempamps(cpar, handle)

     !call output_FITS_sample(cpar, 1000, .true.)

     ! Output zodi ipd and tod parameters to chain
     if (cpar%include_tod_zodi .and. cpar%enable_TOD_analysis) then
         call zodi_model%model_to_chain(cpar, iter)
         do i = 1, numband
            if (data(i)%tod_type == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle
            call output_tod_params_to_hd5(cpar, zodi_model, data(i)%tod, iter)
         end do
     end if

     call wall_time(t2)
     if (ok) then
        if (cpar%myid_chain == 0) then
           write(*,fmt='(a,i4,a,f12.3,a)') ' |  Chain = ', cpar%mychain, ' -- wall time = ', t2-t1, ' sec'
        end if
        iter = iter+1
     else
        if (cpar%myid_chain == 0) then
           write(*,fmt='(a,i4,a,f12.3,a)') ' |  Chain = ', cpar%mychain, ' -- wall time = ', t2-t1, ' sec'
           write(*,*) '|'
           write(*,*) '|  SAMPLE REJECTED'
        end if        
     end if
     
     first = .false.

     call timer%stop(TOT_GIBBSSAMP)
     call timer%incr_numsamp(0)
     !write(*,*) timer%numsamp
     call timer%dumpASCII(cpar%ds_label, trim(cpar%outdir)//"/comm_timing.txt")
  end do

  
  ! **************************************************************
  ! *                   Exit cleanly                             *
  ! **************************************************************

  ! Wait for everybody to exit
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Clean up
  if (cpar%myid == cpar%root .and. cpar%verbosity > 1) write(*,*) '|  Cleaning up and finalizing'

  ! And exit
  call free_status(status)
  call mpi_finalize(ierr)


contains

  subroutine process_TOD(cpar, chain, iter, handle)
    !
    ! Routine for TOD processing
    !
    implicit none
    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: chain, iter
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: i, j, k, l, ndet, ndelta, npar, ierr
    character(len=4)  :: ctext
    character(len=6)  :: samptext
    character(len=512)  :: prefix, postfix, filename
    real(dp)     :: t1, t2, dnu_prop, rms_EE2_prior
    real(dp),      allocatable, dimension(:)     :: eta
    real(dp),      allocatable, dimension(:,:,:) :: delta
    real(dp),      allocatable, dimension(:,:)   :: regnoise
    type(map_ptr), allocatable, dimension(:,:)   :: s_sky, s_gain
    class(comm_map),  pointer :: rms => null()
    class(comm_map),  pointer :: gainmap => null()
    class(comm_comp), pointer :: c => null()
    class(comm_N),    pointer :: N

    ndelta      = cpar%num_bp_prop + 1



    do i = 1, numband  
       if (trim(data(i)%tod_type) == 'none') cycle
       if (iter .ne. 2 .and. mod(iter-1, data(i)%tod_freq) .ne. 0) then
           if (cpar%myid == 0) then
             write(*,fmt='(a,i1,a)') '|  Only processing '//trim(data(i)%label)//' every ',& 
               & data(i)%tod_freq, ' Gibbs samples'
           end if
           cycle
       end if

       if (cpar%myid == 0) then
          write(*,*) '|  ++++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) '|  Processing TOD channel ', trim(data(i)%label)
       end if

       ! Set up EE l=2-subtracted CMB map for absolute gain calibration
       c => compList
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             if (trim(c%label) == 'cmb' .and. c%nmaps > 1) then
                rms_EE2_prior = sqrt(0.308827d-01 * 2*pi/(2.*3.)) / c%cg_scale(2) / c%RJ2unit(2) ! LCDM, Planck 2018 best-fit, uK_cmb^2
                gainmap        => comm_map(c%x)
             end if
          end select
          c => c%next()
       end do


       ! Compute current sky signal for default bandpass and MH proposal
       npar = data(i)%bp(1)%p%npar
       ndet = data(i)%tod%ndet
       allocate(s_sky(ndet,ndelta))
       allocate(s_gain(ndet,1))
       allocate(delta(0:ndet,npar,ndelta))
       allocate(eta(ndet))
       do k = 1, ndelta
          ! Propose new bandpass shifts, and compute mixing matrices
          if (k > 1) then
             if (data(i)%info%myid == 0) then
                do l = 1, npar
                   if (.not. data(i)%tod%sample_abs_bp .or. mod(iter,2) == 0) then
                   !if (.true. .or. mod(iter,2) == 0) then
                      !write(*,*) 'relative',  iter
                      ! Propose only relative changes between detectors, keeping the mean constant
                      delta(0,l,k) = data(i)%bp(0)%p%delta(l)
                      do j = 1, ndet
                         eta(j) = rand_gauss(handle)
                      end do
                      eta = matmul(data(i)%tod%prop_bp(:,:,l), eta)
                     !  write(*,*) "prop_bp: ", data(i)%tod%prop_bp(:,:,l)
                      do j = 1, ndet
                         delta(j,l,k) = data(i)%bp(j)%p%delta(l) + eta(j)
                      end do
                      delta(1:ndet,l,k) = delta(1:ndet,l,k) - mean(delta(1:ndet,l,k)) + &
                           & data(i)%bp(0)%p%delta(l)
                   else
                      !write(*,*) 'absolute',  iter
                      ! Propose only an overall shift in the total bandpass, keeping relative differences constant
                      dnu_prop = data(i)%tod%prop_bp_mean(l) * rand_gauss(handle)
                      do j = 0, ndet
                         delta(j,l,k) = delta(j,l,1) + dnu_prop
                      end do
                   end if
                end do
             end if
             call mpi_bcast(delta(:,:,k), (data(i)%tod%ndet+1)*npar, MPI_DOUBLE_PRECISION, 0, cpar%comm_chain, ierr)
          else
             do j = 0, ndet
                delta(j,:,k) = data(i)%bp(j)%p%delta
             end do
             do l = 1, npar
                delta(1:ndet,l,k) = delta(1:ndet,l,k) - mean(delta(1:ndet,l,k)) + data(i)%bp(0)%p%delta(l)
             end do
          end if

             do j = 0, ndet
                data(i)%bp(j)%p%delta = delta(j,:,k)

                !write(*,*) "delta, j, k: ", delta(j,:,k), j, k
                call data(i)%bp(j)%p%update_tau(data(i)%bp(j)%p%delta)
                if (j > 0 .and. cpar%enable_TOD_analysis .and. data(i)%tod%subtract_zodi) then
                   !write(*,*) 'alloc', i, j, allocated(data(i)%bp(j)%p%nu)
                   call update_zodi_splines(data(i)%tod, data(i)%bp(j), j, zodi_model)
                end if
             end do
             call update_mixing_matrices(i, update_F_int=.true.)

          ! Evaluate sky for each detector given current bandpass
          do j = 1, data(i)%tod%ndet
             !s_sky(j,k)%p => comm_map(data(i)%info)
             call get_sky_signal(i, j, s_sky(j,k)%p, mono=.false.)
             !s_sky(j,k)%p%map = s_sky(j,k)%p%map + 5.d0
             !0call s_sky(j,k)%p%smooth(0.d0, 180.d0)
          end do

          ! Evaluate sky for each detector for absolute gain calibration
          if (k == 1) then
             do j = 1, data(i)%tod%ndet
                if (associated(gainmap)) then
                   call get_sky_signal(i, j, s_gain(j,1)%p, mono=.false., &
                     & abscal_comps=data(i)%tod%abscal_comps, gainmap=gainmap) 
                else
                   call get_sky_signal(i, j, s_gain(j,1)%p, mono=.false.) 
                end if
             end do
          end if

       end do

       !       call s_sky(1,1)%p%writeFITS('sky.fits')

       rms => comm_map(data(i)%rmsinfo)
       call data(i)%tod%process_tod(cpar%outdir, chain, iter, handle, s_sky, delta, data(i)%map, rms, s_gain)

       ! Clear zodi cache for next band
       if (data(i)%subtract_zodi) call data(i)%tod%clear_zodi_cache()
       
       call timer%incr_numsamp(data(i)%id_abs)
       
       if (cpar%myid_chain == 0) then
         write(*,*) '|'
         write(*,*) '|  Finished processing ', trim(data(i)%label)
         write(*,fmt='(a)') '---------------------------------------------------------------------'
         !write(*,*) ''
       end if


       N => data(i)%N
       select type (N)
       class is (comm_N_lcut)
          ! Remove filtered modes
          call int2string(chain, ctext)
          call int2string(iter, samptext)
          prefix = trim(cpar%outdir) // '/tod_' // trim(data(i)%tod%freq) // '_'
          postfix = '_c' // ctext // '_k' // samptext // '.fits'
          call N%P(data(i)%map)
          call data(i)%map%writeFITS(trim(prefix)//'lcut'//trim(postfix))
       end select

       ! Update rms and data maps
       allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
       if (associated(data(i)%procmask)) then
          call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
       else
          call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, map=rms)
       end if
       if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
       !copy data map without regnoise, to write to chain file
       data(i)%map0%map = data(i)%map%map
       data(i)%map%map = data(i)%map%map + regnoise         ! Add regularization noise
       data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask
       deallocate(regnoise)
       call rms%dealloc

       ! Update mixing matrices based on new bandpasses
       do j = 0, data(i)%tod%ndet
          data(i)%bp(j)%p%delta = delta(j,:,1)
          call data(i)%bp(j)%p%update_tau(data(i)%bp(j)%p%delta)
          if (j > 0 .and. cpar%enable_TOD_analysis .and. data(i)%tod%subtract_zodi) call update_zodi_splines(data(i)%tod, data(i)%bp(j), j, zodi_model)
       end do
       call update_mixing_matrices(i, update_F_int=.true.)

       ! Clean up temporary data structures
       do j = 1, data(i)%tod%ndet
          do k = 1, ndelta
             call s_sky(j,k)%p%dealloc
          end do
          call s_gain(j,1)%p%dealloc
       end do
       deallocate(s_sky, s_gain, delta, eta)

       ! Set monopole component to zero, if active. Now part of n_corr
       call nullify_monopole_amp(data(i)%label)
       
    end do
    if (associated(gainmap)) call gainmap%dealloc()

  end subroutine process_TOD

  subroutine print_help()
    print '(a, /)', 'command-line options:'
    print '(a)',    '  -v, --version     print version information and exit'
    print '(a)',    '  -h, --help        print usage information and exit'
  end subroutine print_help

end program commander
