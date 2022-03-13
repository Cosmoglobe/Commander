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
       !  datadir = trim(cpar%datadir)//'/'
       !  filelist = trim(datadir)//trim(cpar%ds_tod_filelist(band))
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
    character(len=512)            :: datadir  !< data directory, which contains all h5 files 
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
      datadir = trim(cpar%datadir)//'/'
      filelist = trim(datadir)//trim(cpar%ds_tod_filelist(band))
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

end module comm_tod_simulations_mod
