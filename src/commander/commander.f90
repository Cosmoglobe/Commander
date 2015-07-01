program commander
  use comm_param_mod
  use comm_data_mod
  use comm_signal_mod
  use comm_cr_mod
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
  type(comm_map), pointer      :: map

  real(dp), allocatable, dimension(:) :: q


  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************
  call read_comm_params(cpar)
  call initialize_mpi_struct(cpar)
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
  !call dump_components('test.dat')
  call dumpCompMaps('test', 'chains')

!!$  map => comm_map(data(1)%info)
!!$  call data(1)%map%writeFITS('in.fits')
!!$  call data(1)%map%YtW
!!$  call data(1)%B%conv(.false., data(1)%map, map)
!!$  call map%Y
!!$  call map%writeFITS('out.fits')

  ! **************************************************************
  ! *                   Carry out computations                   *
  ! **************************************************************

  if (cpar%myid == cpar%root .and. cpar%verbosity > 0) write(*,*) '     Starting Gibbs sampling'

  ! Initialize output structures

  ! Run Gibbs loop
  do iter = 1, cpar%num_gibbs_iter

     ! Sample linear parameters with CG search

     ! Sample amplitude parameters with positivity prior

     ! Sample spectral indices

     ! Sample instrumental parameters

     ! Sample power spectra

     ! Compute goodness-of-fit statistics
     
     ! Output sample to disk
     
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

  recursive function A(x, Nscale)
    use healpix_types
    implicit none
    real(dp), dimension(:),       intent(in)           :: x
    real(dp), dimension(size(x))                       :: A
    real(dp),                     intent(in), optional :: Nscale

    A(1) =  0.5*x(1) - 0.2*x(2)
    A(2) = -0.2*x(1) + 0.8*x(2)
  end function A
  
  recursive function invM(x, Nscale)
    use healpix_types
    implicit none
    real(dp), dimension(:),      intent(in)           :: x
    real(dp), dimension(size(x))                      :: invM
    real(dp),                    intent(in), optional :: Nscale

    invM = x
  end function invM
  
end program commander
