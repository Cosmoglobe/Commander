module comm_tod_simulations_mod
  !
  ! This module contains the collection of subroutines
  ! needed to simulate (e.g. LFI) tods.
  !
  use comm_hdf_mod
  use comm_fft_mod
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_param_mod
  use comm_utils
  !use comm_tod_LFI_mod
  implicit none

contains

   ! ************************************************
   !
   ! TODO: Put this into tod_noise_mod
   !> @brief This routine simulates correlated noise 
   !! component via FFTW. The formula to be used is:
   !!  \f[
   !!     P_n(f) = \sigma_0^2 \left[1 + \left(\frac{f_k}{f}\right)^\alpha\right]
   !!  \f]
   !! The algorithm is the following:
   !! 1. Do an FFT on sigma_0
   !! 2. Multiply this by the formula above
   !! 3. Do an inverse FFT
   !
   !> @author Maksym Brilenkov
   !
   !> @param[in]
   !> @param[out]
   !
   ! ************************************************
!   subroutine simulate_n_corr(self)
!     implicit none
!     class(comm_LFI_tod), intent(inout) :: self !< class instantiation variable
!     type(planck_rng),    intent(inout) :: handle
!     integer(i4b),        intent(in)    :: scan_id !< current PID
!     integer(i4b)                       :: ntod !< total amount of ODs
!     integer(i4b)                       :: ndet !< total amount of detectors

     ! Getting total number of tods and detectors
!     ntod = self%scans(scan_id)%ntod
!     ndet = self%ndet
!     samprate = self%samprate
!     alpha    = self%scans(scan_id)%d(i)%alpha
!     nu_knee  = self%scans(scan_id)%d(i)%fknee
!     N_wn = sigma_0 ** 2  ! white noise power spectrum

     ! Making FFTW plans with given number of threads,
     ! one plan is for forward FFTW and another is for
     ! backward FFTW.
!     nomp = omp_get_max_threads()
!     nfft = 2 * ntod
!     n = nfft / 2 + 1
!     call sfftw_init_threads(err)
!     call sfftw_plan_with_nthreads(nomp)

!     allocate(dt(nfft), dv(0:n-1))
!     call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
!     call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

!     do i = 1, ndet

!     end do

!   end subroutine simulate_n_corr


   subroutine copy_tod(cpar, ierr)
     !
     ! Routine which copies the original hdf5 TODs for
     ! subsequent processing.
     ! 
     ! It first reads-in values stored inside filelist*.txt
     ! to determine the total amount of ODs (i.e. files) to
     ! copy and then uses MPI to invoke multiple system calls
     ! to copy the files into predifined location.
     !
     ! Arguments:
     ! ----------
     !
     ! cpar:     derived type
     !           Object containing parameters from the parameterfile.
     !
     implicit none
     ! Parameter file variables
     type(comm_params), intent(in) :: cpar
     integer(i4b)                  :: id_abs   !< absolute ID of the channel which includes inactive bands
     character(len=512)            :: filelist !< file, which contains correspondance between PIDs and ODs
     character(len=512)            :: simsdir  !< directory where to output simulations 

     !class(comm_LFI_tod), intent(inout) :: self
     ! Simulation routine variables
     integer(i4b) :: unit    !< the current file list value
     integer(i4b) :: n_lines !< total number of raws in the, e.g. filelist_v15.txt file
     integer(i4b) :: n_elem  !< number of unique elements
     integer(i4b) :: val     !< dummy value
     integer(i4b) :: iostatus !< to indicate error status when opening a file
     integer(i4b) :: i, band     !< loop variables
     ! MPI variables
     integer(i4b), intent(in) :: ierr        !< MPI error status
     integer(i4b) :: nprocs !< number of cores
     integer(i4b) :: start_chunk !< Starting iteration value for processor of rank n
     integer(i4b) :: end_chunk   !< End iteration value for processor of rank n
     character(len=256), allocatable, dimension(:) :: input_array  !< array of input h5 file names
     character(len=256), allocatable, dimension(:) :: dummy_array
     character(len=256), allocatable, dimension(:) :: output_array !< array of output h5 file names

     nprocs = cpar%numprocs
     id_abs = cpar%numband
   end subroutine copy_tod


  subroutine write_filelists_to_disk(cpar, ierr)
    !
    ! Routine which copies and overwrites the original
    ! filelist.txt to point to the new (simulation) dir
    !
    ! Arguments:
    ! ----------
    !
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    !
    implicit none
    ! Parameter file variables
    type(comm_params), intent(in) :: cpar
    integer(i4b)                  :: id_abs   !< absolute ID of the channel which includes inactive bands
    character(len=512)            :: filelist !< file, which contains correspondance between PIDs and ODs
    character(len=512)            :: simsdir  !< directory where to output simulations 

    integer(i4b) :: unit    !< the current file list value
    integer(i4b) :: iostatus !< to indicate error status when opening a file
    integer(i4b) :: n_lines !< total number of raws in the, e.g. filelist_v15.txt file
    integer(i4b) :: n_elem  !< number of unique elements
    integer(i4b) :: val     !< dummy value
    integer(i4b) :: i, band     !< loop variables
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: sims_filelist
    character(len=512) :: freq !< central frequency label
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    character(len=256), allocatable, dimension(:) :: input_array  !< array of input h5 file names
    integer(i4b), allocatable, dimension(:) :: pid_array  !< array of input pid labels
    character(len=256), allocatable, dimension(:) :: output_array !< array of output h5 file names
    real(sp), allocatable, dimension(:) :: column3, column4, column5
    ! MPI variables
    integer(i4b), intent(in) :: ierr        !< MPI error status
    integer(i4b) :: nprocs !< number of cores


    nprocs = cpar%numprocs
    id_abs = cpar%numband
    ! looping through all the bands

  end subroutine write_filelists_to_disk

  

  ! ************************************************
  !
  !> @brief Modified version of routine taken from here:
  !!  https://rc.fas.harvard.edu/wp-content/uploads/2013/03/MPI_Plamen_Krastev.pdf
  !!  to calculate the loop range for each processor
  !
  !> @param[in]
  !> @param[out]
  !> @param[inout]
  !
  ! ************************************************
  subroutine split_workload(n1, n2, nprocs, rank, start_chunk, end_chunk)
    implicit none
!    class(comm_LFI_tod), intent(inout) :: self
    integer, intent(in)    :: n1          !< Lowest val of iteration variable
    integer, intent(in)    :: n2          !< Highest values of iteration variable
    integer, intent(inout) :: nprocs      !< # cores
    integer, intent(in)    :: rank        !< processor ID
    integer, intent(out)   :: start_chunk !< Starting iteration value for processor of rank n
    integer, intent(out)   :: end_chunk   !< End iteration value for processor of rank n

    integer :: quotient
    integer :: reminder

    quotient = (n2 - n1 + 1) / nprocs
    reminder = mod((n2 - n1 + 1), nprocs)
    start_chunk = rank * quotient + 1 + min(rank, reminder)
    end_chunk   = start_chunk + quotient - 1
    if(reminder > rank) end_chunk = end_chunk + 1

  end subroutine split_workload 

end module comm_tod_simulations_mod
