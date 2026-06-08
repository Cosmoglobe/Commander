module comm_tod_simulations_mod
  !
  ! This module contains the collection of subroutines
  ! needed to simulate (e.g. LFI) tods.
  !
  use comm_tod_mod
  use comm_hdf_mod
  use comm_fft_mod
!  use comm_shared_arr_mod
!  use spline_1D_mod
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
     integer(i4b), intent(inout) :: ierr        !< MPI error status
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
       filelist = trim(cpar%ds_tod_filelist(band))
       ! central frequency (label)
       !freq = cpar%ds_label(band)
       !write(*,*) "freq is "//trim(freq)
       n_lines = 0
       n_elem  = 0
       val     = 0
       ! processing files only with Master process
       if (cpar%myid == 0) then
         write(*,*) "|   Starting copying files..."
         ! copying an existing filelist into simulation folder
         !call system("cp "//trim(filelist)//" "//trim(simsdir))
         !mystring = filelist
         !mysubstring = ".txt"
         !myindex = index(trim(mystring), trim(mysubstring))
         unit = getlun()
         ! open corresponding filelist, e.g. filelist_30_v15.txt
         open(unit, file=trim(filelist), action="read")
         ! we loop through the file until it reaches its end
         ! (iostatus will give positive number) to get the
         ! total number of lines in the file
         iostatus = 0
         do while(iostatus == 0)
           read(unit,*, iostat=iostatus) val
           if (iostatus == 0) n_lines = n_lines + 1
         end do
         close(unit)
         ! an input array of strings,
         ! which will store filenames
         allocate(input_array(1:n_lines-1))
         ! array which will store pid values
         !allocate(pid_array(1:n_lines))
         write(*,*) "| --------------------------------------------------------------"
         ! once again open the same file to start reading
         ! it from the top to bottom
         open(unit, file=trim(filelist), action="read")
         ! we need to ignore the first line, otherwise it will appear inside an input array
         do i = 0, n_lines-1
           if (i == 0) then
             read(unit,*) val
           else
             read(unit,*) val, input_array(i)
           end if
         end do
         close(unit)
         allocate(dummy_array(size(input_array)))
         write(*,*) "| --------------------------------------------------------------"
         do i = 1, size(input_array)
           ! if the number already exists in result check next
           if (any(dummy_array == input_array(i))) cycle
           ! No match was found, so add it to the output
           n_elem = n_elem + 1
           dummy_array(n_elem) = input_array(i)
         end do
         deallocate(input_array)
         write(*,*) "| --------------------------------------------------------------"
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
       ! Last thing we need to do is to copy and modify filelist*.txt
       ! which points to the new (simulations) directory
       !if (cpar%myid == 0) then
       !  ! if the band is not included then skip it
       !  if (.not. cpar%ds_active(band)) cycle
       !  simsdir = trim(cpar%sims_output_dir)//'/'
       !  filelist = trim(cpar%ds_tod_filelist(band))
       !  ! central frequency (label)
       !  freq = cpar%ds_nu_c(band)!id_abs)
       !  constructor%freq          = cpar%ds_label(id_abs)
       !  constructor%central_freq  = cpar%ds_nu_c(id_abs)
       !  write(*,*) trim(cpar%ds_tod_filelist(band)), band
       !  !new_filelist = trim(simsdir)//trim()
       !  !call system("cp "//trim(filelist)//" "//trim(simsdir))
       !  !call 
       !end if 
       if (cpar%myid == 0) write(*,*) "| Finished copying files!"
       if (cpar%myid == 0) write(*,*) "| --------------------------------------------------------------"
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     end do
     !call MPI_Finalize(ierr)
     !stop
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
    do band = 1, id_abs
      ! if the band is not included then skip it
      if (.not. cpar%ds_active(band)) cycle
      simsdir = trim(cpar%sims_output_dir)//'/'
      filelist = trim(cpar%ds_tod_filelist(band))
      ! processing files only with Master process
      if (cpar%myid == 0) then
        ! GENERALIZE THIS TO NON-LFI DATA
        mysubstring = 'LFI_0'
        mysubstring = 'wmap_'
        n_lines = 0
        n_elem  = 0
        val     = 0
        ! copying an existing filelist and renaming it
        freq = cpar%ds_label(band)
        sims_filelist = trim(simsdir)//"filelist_"//trim(freq)//"_simulations.txt"
        write(*,*) "| filelist is "//trim(filelist)
        write(*,*) "| sims_filelist is "//trim(sims_filelist)

        call system("cp "//trim(filelist)//" "//trim(sims_filelist))
        ! Now, changing pointings inside the file
        unit = getlun()
        ! open corresponding filelist, e.g. filelist_30_v15.txt
        open(unit, file=trim(sims_filelist), action="read")
        ! we loop through the file until it reaches its end
        ! (iostatus will give positive number) to get the
        ! total number of lines in the file
        iostatus = 0
        do while(iostatus == 0)
          read(unit,*, iostat=iostatus) val
          if (iostatus == 0) n_lines = n_lines + 1
        end do
        close(unit)
        allocate(input_array(n_lines), pid_array(n_lines), output_array(n_lines))
        allocate(column3(n_lines), column4(n_lines), column5(n_lines))
        ! array which will store pid values
        !allocate(pid_array(1:n_lines))
        write(*,*) "| --------------------------------------------------------------"
        ! once again open the same file to start reading
        ! it from the top to bottom
        open(unit, file=trim(sims_filelist), action="read")
        ! we need to ignore the first line, otherwise it will appear inside an input array
        ! also, we will get 1 additional line (i.e. the number n_lines will be  1 value larger
        ! than it should be, because of the way we count lines); thus, we need to loop from
        ! 1 to (n_lines - 1), or from 0 to (n_lines - 2)
        do i = 0, n_lines-1
          if (i == 0) then
            read(unit,*) val
          else
            read(unit,*) pid_array(i), input_array(i), column3(i), column4(i), column5(i)
            ! making the new pointing
            mystring = trim(input_array(i))
            myindex = index(trim(mystring), trim(mysubstring))
            output_array(i) = trim(simsdir)//trim(mystring(myindex:))
          end if
        end do
        close(unit)
        ! opening this file again, now to (over)write new values 
        ! (recl=1024 is used so the last column won't be put into
        ! next row)
        open(unit, file=trim(sims_filelist), action="write", recl=1024)
        do i = 0, n_lines-1
          if (i == 0) then
            write(unit,*) val
          else
            write(unit,*) pid_array(i), '"'//trim(output_array(i))//'"', column3(i), column4(i), column5(i)
          end if
        end do
        close(unit)
        write(*,*) "| --------------------------------------------------------------"
        deallocate(input_array, pid_array, output_array)
        deallocate(column3, column4, column5)
      end if
    end do

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



  subroutine simulate_tod(self, scan_id, s_tot, n_corr, handle)
    !
    ! Commander3 native simulation routine. It simulates  correlated
    ! noise component, adds it to the commander-sampled total sky 
    ! signal (multiplied by gain factor for a given frequency) and 
    ! overwrites the original timestreams inside copied files.
    !
    !  Arguments:
    !  ----------
    !  s_tot:    real(sp), array(:,:)
    !            Total sky signal 
    !  scan_id:  integer(i4b)
    !            Local scan ID for the current core
    !  handle:   planck_rng derived type
    !            Healpix definition for random number generation
    !
    !  Returns:
    !  --------
    !
    implicit none
    class(comm_tod), intent(inout) :: self
    ! Parameter file variables
    !type(comm_params),                     intent(in)    :: cpar
    ! Other input/output variables
    real(sp),              dimension(:,:), intent(in)    :: s_tot   !< total sky signal
    real(sp),              dimension(:,:), intent(out)   :: n_corr  !< Correlated noise (output)
    integer(i4b),                          intent(in)    :: scan_id !< current PID
    type(planck_rng),                      intent(inout) :: handle
    ! Simulation variables
    real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
    real(sp)                              :: gain   !< detector's gain value
    real(sp)                              :: sigma0
    real(sp) :: N_c
    real(sp) :: samprate
    real(sp) :: fft_norm
    real(dp) :: chisq
    integer(i4b)                          :: ntod !< total amount of ODs
    integer(i4b)                          :: ndet !< total amount of detectors
    byte,    allocatable, dimension(:)    :: ztod 

    ! HDF5 variables
    character(len=6)   :: samptext, scantext
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=512) :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    type(hdf_file)     :: tod_file
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(hid_t)     :: dtype  ! hdf5 datatype
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err, errorcode !< MPI error status
    integer(i4b)       :: nomp !< Number of threads available
    integer(i4b)       :: omp_err !< OpenMP error status
    integer(i4b) :: omp_get_max_threads
    integer(i4b) :: n, nfft
    integer*8    :: plan_back
    real(sp) :: nu
    !real(sp), allocatable, dimension(:,:) :: n_corr
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    character(len=10) :: processor_label   !< to have a nice output to screen
    integer(i4b) :: ntoks
    character(len=512), dimension(100) :: toks

    !write(*,*) 'sim', self%scanid(scan_id), self%scans(scan_id)%d%accept

    ! shortcuts
    ntod = self%scans(scan_id)%ntod
    ndet = self%ndet

    ! Simulating 1/f noise
    !write(*,*) "Simulating correlated noise"
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i, j, k, dt, dv, sigma0, nu)
    allocate(dt(nfft), dv(0:n-1)) !, n_corr(ntod, ndet))
    !$OMP DO SCHEDULE(guided)
    do j = 1, ndet
      ! skipping iteration if scan was not accepted
      if (.not. self%scans(scan_id)%d(j)%accept) cycle
      ! getting gain for each detector (units, V / K)
      ! (gain is assumed to be CONSTANT for EACH SCAN)
      gain   = self%scans(scan_id)%d(j)%gain
      sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
      samprate = self%samprate
      ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
      fft_norm = sqrt(1.d0 * nfft)
      !
      !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      dv(0) = 0. ! fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) ! HKE: This expression is not correct for the monopole
      do k = 1, (n - 1)
        nu    = k * (samprate / 2) / (n - 1)
        N_c   = self%scans(scan_id)%d(j)%N_psd%eval_corr(nu)
        dv(k) = cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(N_c) /sqrt(2.0)
      end do
      ! Executing Backward FFT
      call timer%start(TOT_FFT)
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      call timer%stop(TOT_FFT)
      dt = dt / sqrt(1.d0*nfft)
      n_corr(:,j) = dt(1:ntod)
      !write(*,*) "n_corr ", n_corr(:, j)
    end do
    !$OMP END DO
    deallocate(dt, dv)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_back)

    ! Allocating main simulations' array
    allocate(tod_per_detector(ntod, ndet))       ! Simulated tod
    tod_per_detector = 1d30

    ! Main simulation loop
    do i = 1, ntod
      do j = 1, ndet
        ! skipping iteration if scan was not accepted
        if (.not. self%scans(scan_id)%d(j)%accept) cycle
        gain   = self%scans(scan_id)%d(j)%gain
        sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
      end do
    end do
    ! Digitizes the data to the nearest integer; probably should mimic the
    ! actual ADC conversion process
    if (self%compressed_tod) tod_per_detector = real(nint(tod_per_detector), kind=sp)

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = '/'

    myindex = index(trim(mystring), trim(mysubstring), back=.true.) + 1


    call get_tokens(trim(mystring), "/", toks=toks, num=ntoks)
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(toks(ntoks))
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    write(*,*) "!  Process:", self%myid, "started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) "!  "//trim(toks(ntoks))

    dims(1) = ntod
    call h5open_f(hdf5_error)
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    if (hdf5_error /= 0) call h5eprint_f(hdf5_error)

    ! Remake huffman, symbols for tod_per_detector
    ! decompress the zipped tods to remake the tod
    !do j = 1, 4
    !   if (.not. self%scans(scan_id)%d(j)%accept) cycle
    !   call self%decompress_tod(scan_id, j, tod_per_detector(:,j))
    !end do
    ! call hufmak(tod_per_detector, self%scans(scan_id)%todkey)
    ! Need to overwrite the keys in the simulated data

    if (self%compressed_tod) then
      call open_hdf_file(trim(self%sims_output_dir)//'/tod_'//pidLabel//'.h5', tod_file, 'w')
      do k = 1, self%ndet
        detectorLabel = self%label(k)

        call write_hdf(tod_file, '/'//trim(detectorLabel), tod_per_detector(:,k))
        call write_hdf(tod_file, '/xi_n_'//trim(detectorLabel), self%scans(scan_id)%d(k)%N_psd%xi_n)
        call write_hdf(tod_file, '/gain_'//trim(detectorLabel), self%scans(scan_id)%d(k)%gain)

      end do
      call write_hdf(tod_file, '/x_im', self%x_im)
      call close_hdf_file(tod_file)
    else
      do j = 1, ndet
        detectorLabel = self%label(j)
        call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      end do
    end if
    call h5dclose_f(dset_id, hdf5_error)
    call h5fclose_f(hdf5_file_id, hdf5_error)
    call h5close_f(hdf5_error)


    deallocate(tod_per_detector)
    write(*,*) "!  Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

  end subroutine simulate_tod


  subroutine simulate_tod_on_the_fly(self, sd, handle)
    ! on_the_fly simulation routine that generates simulated tod,
    ! consisting of sum of correlated and white noise generate here,
    ! and sky signal taken from sd%s_tot multiplied by gain,
    ! and puts it into self%tod (where it overwrites existing tod)
    implicit none
    class(comm_tod),          intent(inout) :: self   ! full data object with compressed tods
    class(comm_scandata),     intent(in)    :: sd     ! decompressed self%scans(scan)%d(det)%tod
    type(planck_rng),         intent(inout) :: handle ! for random number generator

    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp)                                :: gain, sigma0, N_c, samprate, fft_norm, chisq, nu, todsim
    integer(i4b)                            :: ntod, ndet, scan, i, j, k, n, nfft
    integer(i4b)                            :: nomp !< Number of threads available
    integer(i4b)                            :: omp_err !< OpenMP error status
    integer(i4b)                            :: omp_get_max_threads
    integer*8                               :: plan_back

    ! shortcuts
    scan   = sd%scan              ! scan number
    ndet   = sd%ndet              ! number of detectors                = self%ndet
    ntod   = sd%ntod              ! number of tod samples for my scan  = self%scans(scan)%ntod

    !write(*,*) 'sim', scan, self%scans(scan)%d%accept
    
    ! Simulating 1/f noise
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)
    allocate(dt(nfft), dv(0:n-1)) !, n_corr(ntod, ndet))
    do j = 1, ndet
       ! skipping iteration if scan was not accepted
       if (.not. self%scans(scan)%d(j)%accept) cycle
       ! getting gain for each detector (units, V / K)
       ! (gain is assumed to be CONSTANT for EACH SCAN)
       gain   = self%scans(scan)%d(j)%gain
       sigma0 = self%scans(scan)%d(j)%N_psd%sigma0
       samprate = self%samprate
       ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
       fft_norm = sqrt(1.d0 * nfft)
       !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
       ! HKE: expression above is not correct for the monopole
       dv(0) = 0. ! OBS spceial case for the monopole, better to not alter it
       do k = 1, (n - 1)
          nu    = k * (samprate / 2) / (n - 1)
          N_c   = self%scans(scan)%d(j)%N_psd%eval_corr(nu)
          dv(k) = cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(N_c) /sqrt(2.0)
       end do
       ! Executing Backward FFT
       call timer%start(TOT_FFT)
       call sfftw_execute_dft_c2r(plan_back, dv, dt)
       call timer%stop(TOT_FFT)
       dt = dt / sqrt(1.d0*nfft)
            
       ! getting parameters for tod simulations
       gain   = self%scans(scan)%d(j)%gain
       sigma0 = self%scans(scan)%d(j)%N_psd%sigma0
       !write(*,*) 'det,gain,sigma0', j, gain, sigma0
       !write(*,*) "tod before", j, self%scans(scan)%d(j)%tod(1), self%scans(scan)%d(j)%tod(ntod)
       ! make full tod simulations
       do i = 1, ntod
          ! simulated tod = gain * sky_signal + correlated noise + white noise
          ! sky signal taken from sd%s_tot with dimension (ntod, ndet, 0:hmax, nbp))
          todsim = gain * sd%s_tot(i,j,0,1) + dt(i) + sigma0 * rand_gauss(handle)
          ! Digitizes the data to the nearest integer; probably should mimic the
          ! actual ADC conversion process
          if (self%compressed_tod) todsim = real(nint(todsim), kind=sp)
          ! overwrite tod in self by simulation
          self%scans(scan)%d(j)%tod(i) = todsim
       end do
       !write(*,*) "tod after ", j, self%scans(scan)%d(j)%tod(1), self%scans(scan)%d(j)%tod(ntod)

!       if (self%myid==0 .and. j==1) then
!!          filename='dt_'//trim(int2string(j))//'.dat'
!          open(18, file='dt.dat')
!          do i = 1, ntod
!             write(18,*) i, gain *sd%s_tot(i,j,0,1), dt(i), self%scans(scan)%d(j)%tod(i) - gain *sd%s_tot(i,j,0,1) - dt(i),  self%scans(scan)%d(j)%tod(i)
!          end do
!          close(18)
!          write(*,*) 'written to file dt.dat'
!       end if
       
    end do
    deallocate(dv, dt)
    call sfftw_destroy_plan(plan_back)
    
  end subroutine simulate_tod_on_the_fly

end module comm_tod_simulations_mod
