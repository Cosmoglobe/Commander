! ************************************************
!
!> @brief This module contains a collection of
!! subroutines to simulate (e.g. LFI) tods.
!
!> @author Maksym Brilenkov
!
!> @param[in]
!> @param[out]
!
! ************************************************
module comm_tod_simulations_mod
  use comm_hdf_mod
  use comm_fft_mod
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_param_mod
  use comm_utils
  use comm_tod_LFI_mod
  !use comm_tod_mod
  !use comm_map_mod
  !use comm_conviqt_mod
  !use pix_tools
  !use healpix_types
  !use comm_huffman_mod
  !use comm_4D_map_mod
  !use comm_zodi_mod
  !use comm_tod_mapmaking_mod
  !use comm_tod_pointing_mod
  !use comm_tod_gain_mod
  !use comm_tod_bandpass_mod
  !use comm_tod_orbdipole_mod
  implicit none

  !private
  !public comm_LFI_tod

  !integer(i4b), parameter :: N_test      = 19
  !integer(i4b), parameter :: samp_N      = 1
  !integer(i4b), parameter :: prep_G      = 15
  !integer(i4b), parameter :: samp_G      = 2
  !integer(i4b), parameter :: prep_acal   = 3
  !integer(i4b), parameter :: samp_acal   = 4
  !integer(i4b), parameter :: prep_rcal   = 18
  !integer(i4b), parameter :: samp_rcal   = 19
  !integer(i4b), parameter :: prep_relbp  = 5
  !integer(i4b), parameter :: prep_absbp  = 16
  !integer(i4b), parameter :: samp_bp     = 11
  !integer(i4b), parameter :: samp_sl     = 6
  !integer(i4b), parameter :: samp_N_par  = 7
  !integer(i4b), parameter :: sel_data    = 8
  !integer(i4b), parameter :: bin_map     = 9
  !integer(i4b), parameter :: calc_chisq  = 10
  !integer(i4b), parameter :: output_slist = 12
  !integer(i4b), parameter :: samp_mono   = 13
  !integer(i4b), parameter :: sub_sl      = 14
  !integer(i4b), parameter :: sub_zodi    = 17
  !logical(lgt), dimension(N_test) :: do_oper


  !type, extends(comm_tod) :: comm_LFI_tod
  !  class(orbdipole_pointer), allocatable :: orb_dp !orbital dipole calculator
  ! contains
  !   !procedure     :: process_tod        => process_LFI_tod
  !   !----------------------------------------------------------------------------------
  !   ! Simulation Routine
  !   procedure     :: simulate_LFI_tod
  !   procedure     :: copy_LFI_tod
  !   procedure     :: split_workload 
  !   !----------------------------------------------------------------------------------
  !end type comm_LFI_tod

  !interface comm_LFI_tod
  !   procedure constructor
  !end interface comm_LFI_tod

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

   ! ************************************************
   !
   !> @brief Subroutine to retrieve OD-PID dependence.
   !! Reason:
   !! Each OD has a set of PIDs associated with it, 
   !! each of which contains a set of tods for a given 
   !! detector. Each tod is uniquely identified via its
   !! location on the sky (i.e. via angles \theta, \phi
   !! & \psi); thus, we need to know exactly what tod 
   !! located inside what PID (or OD) to correctly
   !! simulate data. Unfortunately, Commander3 doesn't
   !! know/care about OD-PID dependence, so we need a
   !! separate routine for thet.
   !
   !> @author Maksym Brilenkov
   !
   !> @param[in]
   !> @param[out]
   !
   ! ************************************************
   subroutine get_od_pid_dependence(cpar, current_band, nprocs, ierr)
     implicit none
     ! Parameter file variables
     type(comm_params), intent(in) :: cpar
     integer(i4b),      intent(in) :: current_band !< current channel to work on (e.g. 30GHz) 
     character(len=512)            :: filelist !< file, which contains correspondance between PIDs and ODs
     character(len=512)            :: datadir  !< data directory, which contains all h5 files 
     character(len=512)            :: simsdir  !< directory where to output simulations 
     ! Simulation routine variables
     integer(i4b) :: unit    !< the current file list value
     integer(i4b) :: n_lines !< total number of raws in the, e.g. filelist_v15.txt file
     integer(i4b) :: n_elem  !< number of unique elements (i.e. total number of ODs)
     integer(i4b) :: val     !< dummy value
     integer(i4b) :: iostatus !< to indicate error status when opening a file
     integer(i4b) :: i       !< loop variables
     ! MPI variables
     integer(i4b) :: nprocs !< number of cores
     integer(i4b), intent(in) :: ierr   !< MPI error status
     integer(i4b) :: start_chunk !< Starting iteration value for processor of rank n
     integer(i4b) :: end_chunk   !< End iteration value for processor of rank n
     character(len=256), allocatable, dimension(:) :: input_array  !< array of input h5 file names
     character(len=256), allocatable, dimension(:) :: dummy_array
     character(len=256), allocatable, dimension(:) :: output_array !< array of output h5 file names

     simsdir = trim(cpar%sims_output_dir)//'/'
     datadir = trim(cpar%datadir)//'/'
     filelist = trim(datadir)//trim(cpar%ds_tod_filelist(current_band))

     n_lines = 0
     n_elem  = 0
     val     = 0
     ! processing files only with Master process
     if (cpar%myid == 0) then
       write(*,*) "   Starting copying files..."
       unit = getlun()
       ! open corresponding filelist, e.g. filelist_30_v15.txt
       open(unit, file=trim(filelist), action="read")
       ! we loop through the file until it reaches its end
       ! (iostatus will give positive number) to get the
       ! total number of lines in the file
       iostatus = 0
       do while(iostatus == 0)
         read(unit,*, iostat=iostatus) val
         n_lines = n_lines + 1
       end do
       close(unit)
       ! an input array of strings,
       ! which will store filenames
       allocate(input_array(1:n_lines))
       ! array which will store pid values
       !allocate(pid_array(1:n_lines))
       write(*,*) "--------------------------------------------------------------"
       ! once again open the same file to start reading
       ! it from the top to bottom
       open(unit, file=trim(filelist), action="read")
       ! we need to ignore the first line, otherwise it will appear inside an input array
       do i = 0, n_lines-2
         if (i == 0) then
           read(unit,*) val
         else
           read(unit,*) val, input_array(i)
           !pid_array(i) = val
           !write(*,*) "pid_array(i) is ", pid_array(i)
         end if
       end do
       close(unit)
       allocate(dummy_array(size(input_array)))
       write(*,*) "--------------------------------------------------------------"
       do i = 2, size(input_array)
         ! if the number already exists in result check next
         if (any(dummy_array == input_array(i))) cycle
         ! No match was found, so add it to the output
         n_elem = n_elem + 1
         dummy_array(n_elem) = input_array(i)
       end do
       deallocate(input_array)
       write(*,*) "--------------------------------------------------------------"
       ! reducing the size of output array of strings from 45000 to 1490
       allocate(output_array(1:n_elem))
       do i = 1, size(output_array)
         output_array(i) = dummy_array(i)
       end do
       deallocate(dummy_array)
     end if
     ! passing in the array length to all processors
     call MPI_BCAST(n_elem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     ! allocating an array which contains a list of OD names
     if (cpar%myid /= 0) allocate(output_array(n_elem))
     ! mpi passes not a string but each character value,
     ! which means we need to multiply the length of each
     ! path to a file on the value of string length
     call MPI_BCAST(output_array, n_elem * 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     !write(*,*) "n_elem", n_elem
     !write(*,*) "output_array", output_array(1490)
     ! dividing the task to equal (more or less) chunks to loop on
     call split_workload(1, size(output_array), nprocs, cpar%myid, start_chunk, end_chunk)
     ! synchronising processors
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     ! copying all the files with multiprocessing support
     ! each processor has its own chunk of data to work on
     do i = start_chunk, end_chunk
       call system("cp "//trim(output_array(i))//" "//trim(simsdir))
     end do


   end subroutine get_od_pid_dependence

   ! ************************************************
   !
   !> @brief Subroutine to copy original hdf5 files.
   !! It first reads in values stored inside 
   !! filelist*.txt to determine the total amount of
   !! ODs (i.e. files) to copy. And then uses MPI to 
   !! invoke multiple system calls to actually copy
   !! the files to predifined location.
   !
   !> @author Maksym Brilenkov
   !
   !> @param[in]
   !> @param[out]
   !
   ! ************************************************
   subroutine copy_LFI_tod(cpar, ierr)
     implicit none
     ! Parameter file variables
     type(comm_params), intent(in) :: cpar
     integer(i4b)                  :: id_abs   !< absolute ID of the channel which includes inactive bands
     character(len=512)            :: filelist !< file, which contains correspondance between PIDs and ODs
     character(len=512)            :: datadir  !< data directory, which contains all h5 files 
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
     ! looping through all the bands
     do band = 1, id_abs
       ! if the band is not included then skip it
       if (.not. cpar%ds_active(band)) cycle
       simsdir = trim(cpar%sims_output_dir)//'/'
       datadir = trim(cpar%datadir)//'/'
       filelist = trim(datadir)//trim(cpar%ds_tod_filelist(band))
       ! dump one simulation to disc and that is it
       n_lines = 0
       n_elem  = 0
       val     = 0
       ! processing files only with Master process
       if (cpar%myid == 0) then
         write(*,*) "   Starting copying files..."
         unit = getlun()
         ! open corresponding filelist, e.g. filelist_30_v15.txt
         open(unit, file=trim(filelist), action="read")
         ! we loop through the file until it reaches its end
         ! (iostatus will give positive number) to get the
         ! total number of lines in the file
         iostatus = 0
         do while(iostatus == 0)
           read(unit,*, iostat=iostatus) val
           n_lines = n_lines + 1
         end do
         close(unit)
         ! an input array of strings,
         ! which will store filenames
         allocate(input_array(1:n_lines))
         ! array which will store pid values
         !allocate(pid_array(1:n_lines))
         write(*,*) "--------------------------------------------------------------"
         ! once again open the same file to start reading
         ! it from the top to bottom
         open(unit, file=trim(filelist), action="read")
         ! we need to ignore the first line, otherwise it will appear inside an input array
         do i = 0, n_lines-2
           if (i == 0) then
             read(unit,*) val
           else
             read(unit,*) val, input_array(i)
           end if
         end do
         close(unit)
         allocate(dummy_array(size(input_array)))
         write(*,*) "--------------------------------------------------------------"
         do i = 2, size(input_array)
           ! if the number already exists in result check next
           if (any(dummy_array == input_array(i))) cycle
           ! No match was found, so add it to the output
           n_elem = n_elem + 1
           dummy_array(n_elem) = input_array(i)
         end do
         deallocate(input_array)
         write(*,*) "--------------------------------------------------------------"
         ! reducing the size of output array of strings from 45000 to 1490
         allocate(output_array(1:n_elem))
         do i = 1, size(output_array)
           output_array(i) = dummy_array(i)
         end do
         deallocate(dummy_array)
       end if
       ! passing in the array length to all processors
       call MPI_BCAST(n_elem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       ! allocating an array which contains a list of OD names
       if (cpar%myid /= 0) allocate(output_array(n_elem))
       ! mpi passes not a string but each character value,
       ! which means we need to multiply the legth of each
       ! path to a file on the value of string length
       call MPI_BCAST(output_array, n_elem * 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
       !write(*,*) "n_elem", n_elem
       !write(*,*) "output_array", output_array(1490)
       ! dividing the task to equal (more or less) chunks to loop on
       call split_workload(1, size(output_array), nprocs, cpar%myid, start_chunk, end_chunk)
       ! synchronising processors
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       ! copying all the files with multiprocessing support
       ! each processor has its own chunk of data to work on
       do i = start_chunk, end_chunk
         call system("cp "//trim(output_array(i))//" "//trim(simsdir))
       end do
       deallocate(output_array)
       if (cpar%myid == 0) write(*,*) "Finished copying files!"
       if (cpar%myid == 0) write(*,*) "--------------------------------------------------------------"
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     end do
     !call MPI_Finalize(ierr)
     !stop
   end subroutine copy_LFI_tod


   ! ************************************************
   !
   !> @brief Commander3 native simulation module. It
   !! simulates correlated noise and then rewrites
   !! the original timestreams inside the files.
   !
   !> @author Maksym Brilenkov
   !
   !> @param[in]
   !> @param[out]
   !
   ! ************************************************
   subroutine simulate_LFI_tod(cpar, scan_id, s_tot, handle)
     implicit none
     ! Parameter file variables
     type(comm_params),                     intent(in)    :: cpar
     ! Other input/output variables
     real(sp), allocatable, dimension(:,:), intent(in)    :: s_tot   !< total sky signal
     integer(i4b),                          intent(in)    :: scan_id !< current PID
     class(comm_LFI_tod), pointer                   :: ctod
     type(planck_rng),                      intent(inout) :: handle
     ! Simulation variables
     real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
     real(sp)                              :: gain   !< detector's gain value
     real(sp)                              :: sigma0
     integer(i4b)                          :: ntod !< total amount of ODs
     integer(i4b)                          :: ndet !< total amount of detectors
     ! Other vaiables
     integer(i4b)                          :: i, j !< loop variables
     integer(i4b)       :: mpi_err !< MPI error status

     ! shortcuts
     ntod = ctod%scans(scan_id)%ntod
     ndet = ctod%ndet

     ! Allocating main simulations' array
     allocate(tod_per_detector(ntod, ndet))       ! Simulated tod

     ! Main simulation loop
     do i = 1, ntod
       do j = 1, ndet
       ! skipping iteration if scan was not accepted
       if (.not. ctod%scans(scan_id)%d(j)%accept) cycle
       ! getting gain for each detector (units, V / K)
       ! (gain is assumed to be CONSTANT for EACH SCAN)
       gain   = ctod%scans(scan_id)%d(j)%gain
       write(*,*) "gain ", gain
       sigma0 = ctod%scans(scan_id)%d(j)%sigma0
       write(*,*) "sigma0 ", sigma0
       ! Simulating tods
       tod_per_detector(i,j) = gain * s_tot(i,j) + sigma0 * rand_gauss(handle)
       end do
     end do

     !----------------------------------------------------------------------------------
     ! Saving stuff to hdf file
     ! Getting the full path and name of the current hdf file to overwrite

     ! freeing memory up
     deallocate(tod_per_detector)
     call MPI_Finalize(mpi_err)
     stop

!     character(len=*),                      intent(in)    :: sims_output_dir !< output dir for simulated tods
!
!     ! FFTW variables
!     integer*8    :: plan_forward  !< FFT plan for forward transformation
!     integer*8    :: plan_backward !< FFT plan for backward transformation
!     integer(i4b) :: nfft 
!     real(sp),     allocatable, dimension(:) :: dt
!     complex(spc), allocatable, dimension(:) :: dv
!
!     real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
!     real(sp)           :: gain     !< detector's gain value
!     real(sp)           :: sigma0
!     real(sp)           :: progress !< percentage counter
!     character(len=50)  :: message  !< message to pass to progress bar
!     integer(i4b)       :: ntod !< total amount of ODs
!     integer(i4b)       :: ndet !< total amount of detectors
!     integer(i4b)       :: j, k, myindex !< loop variables
!     integer(i4b)       :: start_chunk !< Starting iteration value for processor of rank n
!     integer(i4b)       :: end_chunk   !< End iteration value for processor of rank n
!     integer(i4b)       :: mpi_err !< MPI error status
!     integer(i4b)       :: omp_err !< OMP error status
!
!     ! Doing small string manipulation - retreaving the name of the current file
!     character(len=512) :: mystring, mysubstring, cwd, currentHDFFile
!     character(len=6)   :: pidLabel
!     character(len=3)   :: detectorLabel
!     type(hdf_file)     :: file !< hdf5 file to work with
!     integer(i4b)       :: hdf5_error
!     integer(HID_T)     :: file_id       ! File identifier
!     integer(HID_T)     :: dset_id       ! Dataset identifier
!     integer(HSIZE_T), dimension(1) :: dims
!
!     ntod = self%scans(scan_id)%ntod
!     ndet = self%ndet
     ! Planning FFTW
!     nomp = omp_get_max_threads()
!     nfft = ntod !2 * ntod
!     n = nfft / 2 + 1
!     call sfftw_init_threads(omp_err)
!     call sfftw_plan_with_nthreads(nomp)
!     allocate(dt(nfft), dv(0:n-1))
!     call sfftw_plan_dft_r2c_1d(plan_forward,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
!     call sfftw_plan_dft_c2r_1d(plan_backward, nfft, dv, dt, fftw_estimate + fftw_unaligned)
     ! Need to run detenctor by detector
     ! Populating data arrays
!     do k = 1, ndet !<= here should be ntod or something like that
       ! skipping iteration if scan was not accepted
!       if (.not. self%scans(scan_id)%d(k)%accept) cycle
       ! getting gain for each detector (units, V / K)
       ! (gain is assumed to be CONSTANT for EACH SCAN)
!       gain   = self%scans(scan_id)%d(k)%gain
!       sigma0 = self%scans(scan_id)%d(k)%sigma0
       ! Starting to simulate 1/f for each detector 
!       dt(k) = sigma_0 * rand_gauss(handle)  
!     end do 
     ! Executing Forward FFT
!     call sfftw_execute_dft_r2c(plan_forward, dt, dv)
     ! Multiplying the result by [1 + (f_knee/f)^alpha]
!     do k = 1, ndet !<= here should be ntod or something like that
       ! skipping iteration if scan was not accepted
!       if (.not. self%scans(scan_id)%d(k)%accept) cycle
       ! Getting \alpha, f_knee
!       alpha = self%scans(scan_id)%d(k)%alpha
!       fknee = self%scans(scan_id)%d(k)%fknee

!     end do

!     dv(0) = sigma_0 * rand_gauss(handle)
     ! for all others
!     do k = 1, ndet
!       dv(k) = sigma_0 * cmplx(rand_gauss(handle), rand_gauss(handle)) / sqrt(2) * (1 + (f/fknee)**alpha)
!     end do
!     call sfftw_execute_dft_c2r(plan_forward, dv, dt)

     

!     allocate(tod_per_detector(ntod, ndet))       ! Simulated tod

     ! TODO:
     ! - Make this work with MPI => figure out how commander handles 2d arrays
     ! - Add correlated noise component with MPI support as well
     ! - Change (if needed) the output routine
     !call self%split_workload(1, size(output_array), self%numprocs, self%myid, start_chunk, end_chunk)
     !do i = start_chunk, end_chunk
     !  if (self%myid == 0) then
     !    message  = "Copying files: "
     !    progress = i * 100.0 / (end_chunk - start_chunk + 1)
     !    call pbar%run_progressbar(progress, message, 2)
     ! Dividing workload among processors and simulating tods
!     call self%split_workload(1, ntod, self%numprocs, self%myid, start_chunk, end_chunk)
     !do j = 1, ntod
!     do j = start_chunk, end_chunk
       !if (self%myid == 0) then
       !  message = " Simulating TODs per detector: "
       !  progress = j * 100.0 / (end_chunk - start_chunk + 1)
       !end if
!       do k = 1, ndet
          ! skipping iteration if scan was not accepted
!          if (.not. self%scans(scan_id)%d(k)%accept) cycle
          ! getting gain for each detector (units, V / K)
          ! (gain is assumed to be CONSTANT for EACH SCAN)
!          gain     = self%scans(scan_id)%d(k)%gain
!          sigma0   = self%scans(scan_id)%d(k)%sigma0
          ! Getting sample rate, \nu frequency, and alpha
          ! for each detector for each scan
!          samprate = self%samprate
!          alpha    = self%scans(scan_id)%d(k)%alpha
!          nu_knee  = self%scans(scan_id)%d(k)%fknee
!          N_wn     = sigma_0 ** 2  ! white noise power spectrum
          
!          tod_per_detector(j,k) = gain * s_tot(j,k) + sigma0 * rand_gauss(handle)
          !tod_per_detector(j,k) = gain * s_tot(k,j) + n_corr(k,j) + sigma0 * rand_gauss(handle)
!       end do
!     end do
     ! waiting for all processors to finish
     !call MPI_BARRIER(self%comm, ierr)
     !call MPI_BCAST(tod_per_detector, size(tod_per_detector), MPI_REAL, 0, self%comm, ierr)
     !call MPI_REDUCE(2darray, 2darray(1), 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

     !f_fill = constructor%nobs/(12.*constructor%nside**2)
     !call mpi_reduce(f_fill, f_fill_lim(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, constructor%info%comm, ierr)
     !call mpi_reduce(f_fill, f_fill_lim(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, constructor%info%comm, ierr)
     !call mpi_reduce(f_fill, f_fill_lim(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, constructor%info%comm, ierr)

     !----------------------------------------------------------------------------------
     ! Saving stuff to hdf file
     ! Getting the full path and name of the current hdf file to overwrite
!     mystring = trim(self%hdfname(scan_id))
!     mysubstring = 'LFI_0'
!     myindex = index(trim(mystring), trim(mysubstring))
!     call getcwd(cwd)
!     currentHDFFile = trim(cwd)//'/'//trim(sims_output_dir)//'/'//trim(mystring(myindex:))
!     write(*,*) "Write PID into "//trim(currentHDFFile)
!     ! Converting PID number into string value
!     call int2string(self%scanid(scan_id), pidLabel)
!
!     dims(1) = ntod
!     ! Initialize FORTRAN interface.
!     call h5open_f(hdf5_error)
!     ! Open an existing file - returns file_id
!     call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, file_id, hdf5_error)
!     do j = 1, ndet
!         detectorLabel = self%label(j)
!         ! Create new dataset inside "PID/Detector" Group
!
!         ! Delete group if it already exists
!         !call h5eset_auto_f(0, hdferr)
!         !call h5oget_info_by_name_f(file%filehandle, trim(adjustl(itext)), object_info, hdferr)
!         !if (hdferr == 0) call h5gunlink_f(file%filehandle, trim(adjustl(itext)), hdferr)
!         !write(*,*) 'group ', trim(adjustl(itext))
!         !call create_hdf_group(file, trim(adjustl(itext)))
!         !if (.not. cpar%resamp_CMB) call create_hdf_group(file, trim(adjustl(itext))//'/md')
!
!         ! Open an existing dataset.
!         call h5dopen_f(file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
!         ! Write tod data to a dataset
!         call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
!         ! Close the dataset.
!         call h5dclose_f(dset_id, hdf5_error)
!    end do
!    ! Close the file.
!    call h5fclose_f(file_id, hdf5_error)
!    ! Close FORTRAN interface.
!    call h5close_f(hdf5_error)
!
  end subroutine simulate_LFI_tod


  ! ************************************************
  !
  !> @brief Modified version of procedure, originally
  !!        written by Plamen Krastev 
  !!        (plamenkrastev@fas.harvard.edu), to calculate 
  !!        the loop range for each processor
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
