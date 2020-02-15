program commander
  use comm_param_mod
  use comm_data_mod
  use comm_signal_mod
  use comm_cr_mod
  use comm_chisq_mod
  use comm_output_mod
  use comm_comp_mod
  use comm_nonlin_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2015, all rights reserved                *
  ! *                                                                   *
  ! *                                                                   *
  ! *   NB! The code is provided as is, and *no* guarantees are given   *
  ! *       as far as either accuracy or correctness goes. Even though  *
  ! *       it is fairly well tested, there may be (and likely are)     *
  ! *       bugs in this code.                                          *
  ! *                                                                   *
  ! *  If used for published results, please cite these papers:         *
  ! *                                                                   *
  ! *      - Jewell et al. 2004, ApJ, 609, 1                            *
  ! *      - Wandelt et al. 2004, Phys. Rev. D, 70, 083511              *
  ! *      - Eriksen et al. 2004, ApJS, 155, 227 (Commander)            *
  ! *      - Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB)        *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b)        :: i, iargc, ierr, iter, stat, first_sample, samp_group, curr_samp, tod_freq
  real(dp)            :: t0, t1, t2, t3, dbp
  logical(lgt)        :: ok
  type(comm_params)   :: cpar
  type(planck_rng)    :: handle, handle_noise

  type(comm_mapinfo), pointer :: info => null()
  type(comm_map),     pointer :: m    => null()
  class(comm_comp),   pointer :: c1   => null()

  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************
  call wall_time(t0)
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, cpar%myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, cpar%numprocs, ierr)
  cpar%root = 0
    
  
  if (cpar%myid == cpar%root) call wall_time(t1)
  call read_comm_params(cpar)
  if (cpar%myid == cpar%root) call wall_time(t3)
  
  call initialize_mpi_struct(cpar, handle, handle_noise)
  call validate_params(cpar)  
  call init_status(status, trim(cpar%outdir)//'/comm_status.txt')
  status%active = cpar%myid == 0 !.false.
  
  if (iargc() == 0) then
     if (cpar%myid == cpar%root) write(*,*) 'Usage: commander [parfile] {sample restart}'
     call mpi_finalize(ierr)
     stop
  end if
  if (cpar%myid == cpar%root) call wall_time(t2)

  ! Output a little information to notify the user that something is happening
  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) then
     write(*,*) ''
     write(*,*) '       **********   Commander   *************'
     write(*,*) ''
     write(*,*) '   Number of chains                       = ', cpar%numchain
     write(*,*) '   Number of processors in first chain    = ', cpar%numprocs_chain
     write(*,*) ''
     write(*,fmt='(a,f12.3,a)') '   Time to initialize run = ', t2-t0, ' sec'
     write(*,fmt='(a,f12.3,a)') '   Time to read in parameters = ', t3-t1, ' sec'
     write(*,*) ''

  end if

  ! ************************************************
  ! *               Initialize modules             *
  ! ************************************************

  if (cpar%myid == cpar%root) call wall_time(t1)

  call update_status(status, "init")
  if (cpar%enable_tod_analysis) call initialize_tod_mod(cpar)
  call initialize_bp_mod(cpar);             call update_status(status, "init_bp")
  call initialize_data_mod(cpar, handle);   call update_status(status, "init_data")
  call initialize_signal_mod(cpar);         call update_status(status, "init_signal")
  call initialize_from_chain(cpar, handle); call update_status(status, "init_from_chain")

!write(*,*) 'Setting gain to 1'
!data(6)%gain = 1.d0

  ! Make sure TOD and BP modules agree on initial bandpass parameters
  if (cpar%enable_tod_analysis) call synchronize_bp_delta(trim(cpar%cs_init_inst_hdf) /= 'none')
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
     write(*,*) ''
     write(*,fmt='(a,f12.3,a)') '   Time to read data = ', t2-t1, ' sec'
     write(*,*) '   Starting Gibbs sampling'
  end if


  ! Prepare chains 
  call init_chain_file(cpar, first_sample)
!write(*,*) 'first', first_sample
  !first_sample = 1
  if (first_sample == -1) then
     call output_FITS_sample(cpar, 0, .true.)  ! Output initial point to sample 0
     first_sample = 1
  else
     ! Re-initialise seeds and reinitialize
     call initialize_mpi_struct(cpar, handle, handle_noise, reinit_rng=first_sample)
     !first_sample = 10
     call initialize_from_chain(cpar, handle, init_samp=first_sample, init_from_output=.true.)
     first_sample = first_sample+1
  end if

  !data(1)%bp(0)%p%delta(1) = data(1)%bp(0)%p%delta(1) + 0.2
  !data(2)%bp(0)%p%delta(1) = data(1)%bp(0)%p%delta(1) + 0.2

  ! Run Gibbs loop
  iter = first_sample
  do while (iter <= cpar%num_gibbs_iter)

     if (cpar%myid == 0) then
        call wall_time(t1)
        write(*,fmt='(a)') '---------------------------------------------------------------------'
        write(*,fmt='(a,i4,a,i8)') 'Chain = ', cpar%mychain, ' -- Iteration = ', iter
     end if

     ! Initialize on existing sample if RESAMP_CMB = .true.
     if (cpar%resamp_CMB) then
        if (mod(iter-1,cpar%numsamp_per_resamp) == 0 .or. iter == first_sample) then
           curr_samp = mod((iter-1)/cpar%numsamp_per_resamp,cpar%last_samp_resamp-cpar%first_samp_resamp+1) + cpar%first_samp_resamp
           if (cpar%myid == 0) write(*,*) 'Re-initializing on sample ', curr_samp
           call initialize_from_chain(cpar, handle, init_samp=curr_samp)
           call update_mixing_matrices(update_F_int=.true.)       
        end if
     end if

     ! Process TOD structures
     if (cpar%enable_TOD_analysis .and. (iter <= 2 .or. mod(iter,cpar%tod_freq) == 0)) then
        call process_TOD(cpar, cpar%mychain, iter, handle)
     end if

     ! Sample linear parameters with CG search; loop over CG sample groups
     if (cpar%sample_signal_amplitudes) then
        do samp_group = 1, cpar%cg_num_samp_groups
           if (cpar%myid == 0) then
              write(*,fmt='(a,i4,a,i4,a,i4)') '  Chain = ', cpar%mychain, ' -- CG sample group = ', &
                   & samp_group, ' of ', cpar%cg_num_samp_groups
           end if
           call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

           if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

        end do

        ! Perform joint alm-Cl Metropolis move
        do i = 1, 3
           if (cpar%resamp_CMB) call sample_joint_alm_Cl(handle)
        end do

     end if

     ! Sample power spectra
     call sample_powspec(handle, ok)

     ! Output sample to disk
     call output_FITS_sample(cpar, iter, .true.)

     ! Sample partial-sky templates
     !call sample_partialsky_tempamps(cpar, handle)

     !call output_FITS_sample(cpar, 1000, .true.)

     ! Sample non-linear parameters
     if (cpar%sample_specind) then
        do i = 1, cpar%num_ind_cycle
           call sample_nonlin_params(cpar, iter, handle)
        end do
     end if

     ! Sample instrumental parameters

     

     call wall_time(t2)
     if (ok) then
        if (cpar%myid == 0) then
           write(*,fmt='(a,i4,a,f12.3,a)') 'Chain = ', cpar%mychain, ' -- wall time = ', t2-t1, ' sec'
        end if
        iter = iter+1
     else
        if (cpar%myid == 0) then
           write(*,fmt='(a,i4,a,f12.3,a)') 'Chain = ', cpar%mychain, ' -- wall time = ', t2-t1, ' sec'
           write(*,*) 'SAMPLE REJECTED'
        end if        
     end if
     
  end do

  
  ! **************************************************************
  ! *                   Exit cleanly                             *
  ! **************************************************************

  ! Wait for everybody to exit
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Clean up
  if (cpar%myid == cpar%root .and. cpar%verbosity > 1) write(*,*) '     Cleaning up and finalizing'

  ! And exit
  call free_status(status)
  call mpi_finalize(ierr)


contains

  subroutine process_TOD(cpar, chain, iter, handle)
    implicit none
    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: chain, iter
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: i, j, k, l, ndet, ndelta, npar, ierr
    real(dp)     :: t1, t2, dnu_prop
    real(dp),      allocatable, dimension(:)     :: eta
    real(dp),      allocatable, dimension(:,:,:) :: delta
    real(dp),      allocatable, dimension(:,:)   :: regnoise
    type(map_ptr), allocatable, dimension(:,:)   :: s_sky
    class(comm_map), pointer :: rms => null()

    ndelta      = cpar%num_bp_prop + 1

    do i = 1, numband  
       if (trim(data(i)%tod_type) == 'none') cycle

       ! Compute current sky signal for default bandpass and MH proposal
       npar = data(i)%bp(1)%p%npar
       ndet = data(i)%tod%ndet
       allocate(s_sky(ndet,ndelta))
       allocate(delta(0:ndet,npar,ndelta))
       allocate(eta(ndet))
       do k = 1, ndelta
          ! Propose new bandpass shifts, and compute mixing matrices
          if (k > 1) then
             if (cpar%myid == 0) then
                do l = 1, npar
                   if (mod(iter,2) == 0) then
                      write(*,*) 'relative',  iter
                      ! Propose only relative changes between detectors, keeping the mean constant
                      delta(0,l,k) = data(i)%bp(0)%p%delta(l)
                      do j = 1, ndet
                         eta(j) = rand_gauss(handle)
                      end do
                      eta = matmul(data(i)%tod%prop_bp(:,:,l), eta)
                      do j = 1, ndet
                         delta(j,l,k) = data(i)%bp(j)%p%delta(l) + eta(j)
                      end do
                      delta(1:ndet,l,k) = delta(1:ndet,l,k) - mean(delta(1:ndet,l,k)) + &
                           & data(i)%bp(0)%p%delta(l)
                   else
                      write(*,*) 'absolute',  iter
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

          ! Update mixing matrices
          if (k > 1 .or. iter == 1) then
             do j = 0, ndet
                data(i)%bp(j)%p%delta = delta(j,:,k)
                call data(i)%bp(j)%p%update_tau(delta(j,:,k))
             end do
             call update_mixing_matrices(i, update_F_int=.true.)       
          end if

          ! Evaluate sky for each detector given current bandpass
          do j = 1, data(i)%tod%ndet
             !s_sky(j,k)%p => comm_map(data(i)%info)
             call get_sky_signal(i, j, s_sky(j,k)%p, mono=.false.) 
             !0call s_sky(j,k)%p%smooth(0.d0, 180.d0)
          end do

       end do

!       call s_sky(1,1)%p%writeFITS('sky.fits')

       ! Process TOD, get new map. TODO: update RMS of smoothed maps as well. 
       ! Needs in-code computation of smoothed RMS maps, so long-term..
       rms => comm_map(data(i)%info)
       call data(i)%tod%process_tod(cpar%outdir, chain, iter, handle, s_sky, delta, data(i)%map, rms)

       ! Update rms and data maps
       allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
       if (associated(data(i)%procmask)) then
          call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
       else
          call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, map=rms)
       end if
       if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
       data(i)%map%map = data(i)%map%map + regnoise         ! Add regularization noise
       data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask
       deallocate(regnoise)
       call rms%dealloc

       ! Update mixing matrices based on new bandpasses
       do j = 0, data(i)%tod%ndet
          data(i)%bp(j)%p%delta = delta(j,:,1)
          call data(i)%bp(j)%p%update_tau(delta(j,:,1))
       end do
       call update_mixing_matrices(i, update_F_int=.true.)       

       ! Clean up temporary data structures
       do k = 1, ndelta
          do j = 1, data(i)%tod%ndet
             call s_sky(j,k)%p%dealloc
          end do
       end do
       deallocate(s_sky, delta, eta)
    end do

  end subroutine process_TOD

end program commander
