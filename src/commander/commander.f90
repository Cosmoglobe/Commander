program commander
  use comm_param_mod
  use comm_data_mod
  use comm_signal_mod
  use comm_cr_mod
  use comm_chisq_mod
  use comm_output_mod
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

  integer(i4b)        :: iargc, ierr, iter, stat
  type(comm_params)   :: cpar
  type(planck_rng)    :: handle

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
  call initialize_bp_mod(cpar);          call update_status(status, "init_bp")
  call initialize_data_mod(cpar);        call update_status(status, "init_data")
  call initialize_signal_mod(cpar);      call update_status(status, "init_signal")
  
  ! **************************************************************
  ! *                   Carry out computations                   *
  ! **************************************************************

  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) write(*,*) '     Starting Gibbs sampling'

  ! Initialize output structures

  ! Run Gibbs loop
  do iter = 1, cpar%num_gibbs_iter

     ! Sample linear parameters with CG search
     call sample_amps_by_CG(cpar, handle)

     ! Sample amplitude parameters with positivity prior

     ! Sample spectral indices

     ! Sample instrumental parameters

     ! Sample power spectra

     ! Compute goodness-of-fit statistics
     
     ! Output sample to disk
     call output_FITS_sample(cpar, iter)
     
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
