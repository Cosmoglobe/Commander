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

  integer(i4b)        :: iargc, ierr, iter, stat, first_sample
  type(comm_params)   :: cpar
  type(planck_rng)    :: handle

  type(comm_mapinfo), pointer :: info
  type(comm_map), pointer :: m
  class(comm_comp), pointer :: c1

  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************
  call read_comm_params(cpar)
  call initialize_mpi_struct(cpar, handle)
  call init_status(status, 'comm_status.txt')
  
  if (iargc() == 0) then
     if (cpar%myid == cpar%root) write(*,*) 'Usage: commander [parfile] {sample restart}'
     call mpi_finalize(ierr)
     stop
  end if

  ! Output a little information to notify the user that something is happening
  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) then
     write(*,*) ''
     write(*,*) '       **********   Commander   *************'
     write(*,*) ''
     write(*,*) '   Number of chains                       = ', cpar%numchain
     write(*,*) '   Number of processors in first chain    = ', cpar%numprocs_chain
     write(*,*) ''
  end if

  ! ************************************************
  ! *               Initialize modules             *
  ! ************************************************

  call update_status(status, "init")
  call initialize_bp_mod(cpar);            call update_status(status, "init_bp")
  call initialize_data_mod(cpar, handle);  call update_status(status, "init_data")
  call initialize_signal_mod(cpar);        call update_status(status, "init_signal")
  call initialize_from_chain(cpar)

  if (cpar%output_input_model) then
     if (cpar%myid == 0) write(*,*) 'Outputting input model to sample number 999999'
     call output_FITS_sample(cpar, 999999, .false.)
     call mpi_finalize(ierr)
     stop
  end if

  ! Output SEDs for each component
  if (.false.) then
     if (cpar%myid == cpar%root) call dump_components('sed.dat')
     call mpi_finalize(ierr)
     stop
  end if
  
  ! **************************************************************
  ! *                   Carry out computations                   *
  ! **************************************************************

  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) write(*,*) '     Starting Gibbs sampling'

  ! Initialize output structures

  ! Run Gibbs loop
  first_sample = 1
  if (trim(cpar%init_chain_prefix) == trim(cpar%chain_prefix)) first_sample = cpar%init_samp+1
  do iter = first_sample, cpar%num_gibbs_iter

     ! Sample linear parameters with CG search
     call sample_amps_by_CG(cpar, handle)

     ! Sample non-linear parameters
     call sample_nonlin_params(cpar, handle)

     ! Sample instrumental parameters

     ! Sample power spectra

     ! Compute goodness-of-fit statistics
     
     ! Output sample to disk
     call output_FITS_sample(cpar, iter, .true.)
     
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

end program commander
