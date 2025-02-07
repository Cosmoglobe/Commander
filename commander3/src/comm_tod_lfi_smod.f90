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
submodule (comm_tod_lfi_mod) comm_tod_lfi_smod
contains

  !**************************************************
  !             Constructor
  !**************************************************
  
  
  module subroutine load_instrument_lfi(self, instfile, band)
    !
    ! Reads the LFI specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_LFI_tod
    !    the LFI tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_lfi_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band

    integer(i4b) :: i, j
    real(dp) :: weight
    integer(i4b), dimension(1) :: n_gmf
    character(len=1) :: id
    character(len=512) :: path
    real(dp), dimension(:), pointer :: gmf_s

    ! Read in mainbeam_eff
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'mbeam_eff', self%mb_eff(band))


    if(trim(self%level) == 'L1') then
      if(index(self%label(band), 'M') /= 0) then
       self%diode_names(band,:) = ['ref00','sky00','ref01','sky01']
       id = '0'
      else
       self%diode_names(band,:) = ['ref10','sky10','ref11','sky11']
       id = '1'
      end if

      ! read in the r checkpoints
      call get_size_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'gmfSplits', n_gmf)

      allocate(gmf_s(n_gmf(1)))

      call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'gmfSplits', gmf_s) 
      self%gmf_splits(band)%p => gmf_s

      !if(self%myid == 0) write(*,*) trim(self%label(band)), gmf_s
      ! read in the diode weights
      call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'diodeWeight', weight)
      self%diode_weights(band,1) = weight
      self%diode_weights(band,2) = 1.d0 - weight

      if (.not. self%L2_exist) then
        ! Read ADC corrections
        path = trim(adjustl(self%label(band)))//'/'//'adc91-'//id//'0'
      end if
   end if

  end subroutine load_instrument_lfi
 
  module subroutine initHDF_lfi(self, chainfile, path)
    ! 
    ! Initializes instrument-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
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
    class(comm_lfi_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
    
    real(dp), allocatable, dimension(:,:,:,:)           :: ref_filter
    integer(i4b)  :: i, j
    real(dp), dimension(2)                              :: boundary

    boundary = (0.d0, 0.d0) 


  end subroutine initHDF_lfi

 
  module subroutine diode2tod_lfi(self, scan, map_sky, procmask, tod)
    ! 
    ! Generates detector-coadded TOD from low-level diode data
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! scan:     int
    !           Scan ID number
    ! procmask: array of sp
    !           processing mask that cuts out the galaxy
    !
    ! Returns
    ! ----------
    ! tod:      ntod x ndet sp array
    !           Output detector TOD generated from raw diode data
    !
    implicit none
    class(comm_lfi_tod),                       intent(inout) :: self
    integer(i4b),                              intent(in)    :: scan
    real(sp),          dimension(0:,1:,1:,1:), intent(in)    :: map_sky
    real(sp),          dimension(0:),          intent(in)    :: procmask
    real(sp),          dimension(:,:),         intent(out)   :: tod

    integer(i4b) :: i,j,k,half,horn,n_mask, n_unmask, err, nsmooth
    real(sp), allocatable, dimension(:,:) :: diode_data, corrected_data, s_sky, mask, differenced_data
    real(dp), allocatable, dimension(:) :: nu_out
    real(dp), allocatable, dimension(:,:) ::  binned_corr
    integer(i4b), allocatable, dimension(:,:,:) :: pix, psi
    integer(i4b), allocatable, dimension(:,:)   :: flag
    real(dp) :: r1, r2, sum1, sum2, A(3,3,2), b(3,2), x(3,2), t1
    logical(lgt) :: gmf_split
    character(len=1024) :: filename

    nsmooth = self%get_nsmooth()

    allocate(diode_data(self%scans(scan)%ntod, self%ndiode))
    allocate(corrected_data(self%scans(scan)%ntod, self%ndiode))
    allocate(differenced_data(self%scans(scan)%ntod, self%ndiode/2))
    allocate(nu_out(nsmooth), binned_corr(1, nsmooth))
    allocate(pix(self%scans(scan)%ntod,self%ndet,self%nhorn), psi(self%scans(scan)%ntod,self%ndet,self%nhorn), flag(self%scans(scan)%ntod,self%ndet))
    allocate(s_sky(self%scans(scan)%ntod, self%ndet), mask(self%scans(scan)%ntod, self%ndet))

    call self%get_freq_bins(nu_out)

    diode_data = 0.0

    do i=1, self%ndet
      call self%decompress_pointing_and_flags(scan, i, pix(:,i,:), psi(:,i,:), flag(:,i))
    end do
    call project_sky(self, map_sky(:,:,:,1), pix(:,:,1), psi(:,:,1), flag, &
         & procmask, scan, s_sky, mask)

    do i=1, self%ndet

       if (.not. self%scans(scan)%d(i)%accept) cycle

       ! check if this is one of the weird chunks with 2 gain modulation factors
       gmf_split = .false.

       t1 = self%scans(scan)%t0(2) + 2**16 * self%scans(scan)%ntod / self%samprate
       do k = 1, size(self%gmf_splits(i)%p)
        !write(*,*) "Time", self%scans(scan)%t0(2), t1, self%gmf_splits(i)%p(k), self%samprate, self%scans(scan)%ntod, size(self%gmf_splits(i)%p)
        if (self%gmf_splits(i)%p(k) > self%scans(scan)%t0(2) .and. self%gmf_splits(i)%p(k) < t1) then
          gmf_split = .true.
          exit
        end if
       end do

       if(gmf_split) then 
        !self%scans(scan)%d(i)%accept = .false.
        write(*,*) trim(self%label(i)), "| Not cutting scan", self%scans(scan)%chunk_num, "because of gmf split", self%scans(scan)%t0(2), t1, self%gmf_splits(i)%p(k), k, size(self%gmf_splits(i)%p) 
        !cycle
       end if

        ! Decompress diode TOD for current scan
        call self%decompress_diodes(scan, i, diode_data)

        ! Apply ADC corrections
        ! Wiener-filter load data (do not apply load smoothing in BP10)
        !call self%filter_reference_load(i, corrected_data)

        ! Compute the gain modulation factors

        if(self%use_dpc_gain_modulation) then

          r1 = 0.d0
          r2 = 0.d0
          sum1 = 0.d0
          sum2 = 0.d0
          n_mask = 0
          n_unmask = 0
          
          do k = 1, size(corrected_data(:,1))
            if (mask(k,i) == 0.) cycle

            sum1 = sum1 + corrected_data(k,1)
            sum2 = sum2 + corrected_data(k,3)
            n_unmask = n_unmask + 1

            r1 = r1 + corrected_data(k,2)
            r2 = r2 + corrected_data(k,4)
            n_mask = n_mask + 1

          end do
!
          if (r1 == 0.d0 .or. r2 == 0.d0 .or. sum1 == 0.d0 .or. sum2 == 0.d0) then
             self%scans(scan)%d(i)%accept = .false.
             cycle
          end if
          self%R(scan,i,1) = (r1/n_mask)/(sum1/n_unmask)
          self%R(scan,i,2) = (r2/n_mask)/(sum2/n_unmask)
    
       else ! use fancy new gain modulation factor computation

        A = 0.d0
        b = 0.d0
        do k = 1, size(corrected_data(:,1))
          if (mask(k,i) == 0.) cycle

          A(1,1,1) = A(1,1,1) + 1.d0
          A(1,2,1) = A(1,2,1) + s_sky(k,i)
          A(1,3,1) = A(1,3,1) + corrected_data(k,1)
          A(2,2,1) = A(2,2,1) + s_sky(k,i)          * s_sky(k,i)
          A(2,3,1) = A(2,3,1) + s_sky(k,i)          * corrected_data(k,1)
          A(3,3,1) = A(3,3,1) + corrected_data(k,1) * corrected_data(k,1)
          b(1,1)   = b(1,1)   +                       corrected_data(k,2)
          b(2,1)   = b(2,1)   + s_sky(k,i)          * corrected_data(k,2)
          b(3,1)   = b(3,1)   + corrected_data(k,1) * corrected_data(k,2)

          A(1,1,2) = A(1,1,2) + 1.d0
          A(1,2,2) = A(1,2,2) + s_sky(k,i)
          A(1,3,2) = A(1,3,2) + corrected_data(k,3)
          A(2,2,2) = A(2,2,2) + s_sky(k,i)          * s_sky(k,i)
          A(2,3,2) = A(2,3,2) + s_sky(k,i)          * corrected_data(k,3)
          A(3,3,2) = A(3,3,2) + corrected_data(k,3) * corrected_data(k,3)
          b(1,2)   = b(1,2)   +                       corrected_data(k,4)
          b(2,2)   = b(2,2)   + s_sky(k,i)          * corrected_data(k,4)
          b(3,2)   = b(3,2)   + corrected_data(k,3) * corrected_data(k,4)
       end do
       !if (A(1,1,1) == 0.d0) then
       !    self%scans(scan)%d(i)%accept = .false.
       !    cycle
       ! end if        
        do j = 1, 3
           do k = j+1, 3
              A(k,j,:) = A(j,k,:)
           end do
        end do
        call solve_system_real(A(:,:,1), x(:,1), b(:,1))
        call solve_system_real(A(:,:,2), x(:,2), b(:,2))
        
        ! average sky value/average load value
        self%R(scan,i,1) = x(3,1)
        self%R(scan,i,2) = x(3,2)
        !if( self%scanid(scan) == 27676) then
        !  write(*,*) "new, old:", x(3,1), x(3,2), (r1/n_mask)/(sum1/n_unmask), (r2/n_mask)/(sum2/n_unmask) 
        !end if
      end if

        ! Compute output differenced TOD

        !w1(sky00 - R*ref00) + w2(sky01 - R*ref01)
        !if(self%myid == 0) write(*,*) r1, r2, n_mask, size(diode_data(:,1))
        !tod(:,i) = self%diode_weights(i,1) * (diode_data(:,2) - self%R(scan,i,1) * diode_data(:,1)) + self%diode_weights(i,2)*(diode_data(:,4) - self%R(scan,i,2) * diode_data(:,3))
!        do k = 0, 100
!           self%R(scan,i,2) = 0.70 + 0.003*k
           
            !determine cross corrlation between the two diodes
            !binned_corr = 0.d0
            !differenced_data(:,1) = corrected_data(:,2) - self%R(scan,i,1) * corrected_data(:,1)
            !differenced_data(:,2) = corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3)
            !call self%compute_ref_load_filter(differenced_data, binned_corr, nu_out, err)

            !write(filename, '(A12,I6.6,A3,A4)') trim('diode_xcorr_'), self%scans(scan)%chunk_num, trim(self%label(i)),  trim('.dat')
            !open(58 + self%myid, file=filename)
            !do j = 1, size(nu_out) 
            !  write(58 + self%myid, *) nu_out(j), binned_corr(1,j)
            !end do
            !close(58 + self%myid)

            tod(:,i) = self%diode_weights(i,1) * (corrected_data(:,2) - self%R(scan,i,1) * corrected_data(:,1)) + self%diode_weights(i,2)*( corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3))

            !tod(:,i) = self%diode_weights(i,2)*(corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3))

!           tod(:,i) = corrected_data(:,4) - self%R(scan,i,2) * corrected_data(:,3)
!           write(*,*) self%R(scan,i,2), sum(tod(:,i)**2), variance(1.d0*tod(:,i))
!        end do
        !stop

!!$    if (self%scanid(scan) == 3 .and. i == 1) then
!!$        open(58,file='res_L2fromL1_030_pid3.dat', recl=1024)
!!$        do j = 1, size(tod,1)
!!$           write(58,*) j, corrected_data(j,2) - x(1,1) - s_sky(j,1)*x(2,1) - corrected_data(j,1)*x(3,1)!, corrected_data(j,:)
!!$        end do
!!$        close(58)
!!$     end if


        !tod(:,i) = self%diode_weights(i,1) * (corrected_data(:,2) - filtered_data(:,1)) + self%diode_weights(i,2)*( corrected_data(:,4) - filtered_data(:,3))
        !tod(:,i) = (corrected_data(:,1) - corrected_data(:,3)) + (corrected_data(:,2) - corrected_data(:,4))
        

        !stop
        tod(:,i) = tod(:,i) - (sum(tod(:,i))/size(tod(:,i)))
 
    end do
!    stop

    !if (self%scanid(scan) == 3) then
    !    open(58,file='comm3_L2fromL1_030_pid3.dat', recl=1024)
    !    do j = 1, size(tod,1)
    !       write(58,*) j, tod(j,:), diode_data(j,:)!, corrected_data(j,:)
    !    end do
    !    close(58)
    ! end if

    deallocate(diode_data, corrected_data, pix, psi, flag, s_sky, mask)
    deallocate(differenced_data, nu_out, binned_corr)

!stop

!call mpi_finalize(i)
!stop

  end subroutine diode2tod_lfi

  module function get_nsmooth(self)
    implicit none
    class(comm_lfi_tod),  intent(in)   :: self
    integer(i4b)                       :: get_nsmooth  
    integer(i4b) :: j
    real(sp)     :: fbin, nu

    fbin         = 1.2 ! multiplicative bin scaling factor
    get_nsmooth  = 1
    nu           = 0.01
    do while (nu < self%samprate/2)
       get_nsmooth = get_nsmooth + 1
       nu          = nu * fbin
    end do
  end function get_nsmooth

  module subroutine get_freq_bins(self, freqs)
    implicit none
    class(comm_lfi_tod),    intent(in)    :: self
    real(dp), dimension(:), intent(inout) :: freqs  
 
    integer(i4b) :: j
    real(sp)     :: fbin, nu

    fbin         = 1.2 ! multiplicative bin scaling factor
    nu           = 0.01
    j            = 1
    do while (nu < self%samprate/2)
       freqs(j)    = nu
       j           = j + 1
       nu          = nu * fbin
    end do
    freqs(j) = self%samprate/2

  end subroutine get_freq_bins

  module subroutine compute_ref_load_filter(self, data_in, binned_out, nu_out, err)
    ! 
    ! Computes the binned weiner filter for the reference load
    !
    ! Arguments:
    ! ----------
    ! 
    ! self:     comm_tod_LFI object
    !           TOD processing class
    ! data_in:  float array (ntod, ndiode)
    !           input diode timestreams
    !
    ! Returns:
    ! --------
    !
    ! binned_out : float array
    !              array of filter transfer function for ref load
    ! nu_out     : float_array
    !              frequencies that index binned_out
    ! err        : error flag; 0 if OK, 1 if no data
    implicit none
    class(comm_lfi_tod),          intent(in)    :: self
    real(sp),     dimension(:,:), intent(in)    :: data_in
    real(dp),     dimension(:,:), intent(inout) :: binned_out
    real(dp),     dimension(:),   intent(in)    :: nu_out
    integer(i4b),                 intent(out)   :: err

    integer(i4b) :: i, j, k, nfft, n, n_bin
    real(dp)     :: num, denom, fsamp, fbin, nu, upper, subsum, nu_low, delta_nu, sum_ref, sum_sky
    integer*8    :: plan_fwd

    real(sp),     allocatable, dimension(:) :: dt_sky, dt_ref
    real(dp),     allocatable, dimension(:) :: filter
    complex(spc), allocatable, dimension(:) :: dv_sky, dv_ref

    ! This test should be replaced with something more fine-tuned
    if (all(data_in == 0.)) then
       err = 1
       return
    else
       err = 0
    end if

    n       = size(data_in(:,1))
    nfft    = n/2+1
    fsamp   = self%samprate

    allocate(dt_sky(n), dt_ref(n), dv_sky(0:nfft-1), dv_ref(0:nfft-1), filter(nfft-1))
   


    do i = 1, size(data_in(1,:))/2

       dt_ref = data_in(:, 2*i-1)
       dt_sky = data_in(:, 2*i)

       sum_ref = sum(dt_ref)
       sum_sky = sum(dt_sky)

      ! FFT of ref signal
      call sfftw_execute_dft_r2c(plan_fwd, dt_ref, dv_ref)

      ! FFT of sky signal
      call sfftw_execute_dft_r2c(plan_fwd, dt_sky, dv_sky)     

      ! Compute cross correlation
      do j = 1, nfft-1
         num =  real(dv_sky(j)*conjg(dv_ref(j)), dp)
         denom = real(abs(dv_sky(j)) * abs(dv_ref(j)), dp)
         if (denom < 1d-100) then
            filter(j) = 0.
         else 
            filter(j) = num/denom
         end if

      end do

      ! bin into bins defined by nu_out
      j = 1
      delta_nu     = ind2freq(2, fsamp, nfft)
      k            = nint(0.01d0/delta_nu)    ! First frequency to consider
      nu           = ind2freq(k, fsamp, nfft) ! Current frequency

      ! loop through all bins
      do j = 1, size(nu_out) - 1
        n_bin = 0
        subsum = 0.d0
        ! loop through all the frequencies in that bin
        do while (nu < nu_out(j+1) .and. k < size(filter))
          subsum = subsum + filter(k)
          k = k+1
          n_bin = n_bin + 1
          nu = nu + delta_nu
        end do
        if(n_bin > 0) then
          ! average the binned values
          binned_out(i, j) = binned_out(i, j) + subsum/n_bin
        end if
      end do 

    end do

    call sfftw_destroy_plan(plan_fwd)

    deallocate(dt_sky, dt_ref, dv_sky, dv_ref, filter)

  end subroutine compute_ref_load_filter


  module subroutine filter_reference_load(self, det, data)
    class(comm_lfi_tod),               intent(in)      :: self
    integer(i4b),                      intent(in)      :: det
    real(sp), dimension(:,:),          intent(inout)   :: data

    real(dp)     :: filt
    integer(i4b) :: i, j, nfft, n
    integer*8    :: plan_fwd, plan_back

    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv


    n       = size(data(:,1))
    nfft    = n/2+1

    allocate(dt(n), dv(0:nfft-1))

    call sfftw_plan_dft_r2c_1d(plan_fwd,  n, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, n, dv, dt, fftw_estimate + fftw_unaligned)

!!$    open(58,file='raw.dat')
!!$    do i = 1, n
!!$       write(58,*) data(i,:)
!!$    end do
!!$    close(58)

    do i = 1, self%ndiode/2

      ! Check if data is all zeros
      dt = data(:, 2*i -1)
      if(all(dt == 0)) cycle

      ! FFT of ref signalA

      ! Filter ref with cross correlation transfer function
!      open(58,file='filter.dat')
      do j=1, size(dv) -1
        if(ind2freq(j, self%samprate, nfft) < 7.d0) filt = 1.d0 ! removes regions where we don't want to filter because it actually adds noise somehow
        dv(j) = dv(j) * filt
      end do
!     close(58)

      ! IFFT ref signal
      
      ! Normalize
      data(:, 2*i-1) = dt/n

    end do

    call sfftw_destroy_plan(plan_fwd)
    call sfftw_destroy_plan(plan_back)

    deallocate(dt, dv)

!!$    open(58,file='filtered.dat')
!!$    do i = 1, n
!!$       write(58,*) data(i,:)
!!$    end do
!!$    close(58)
!!$    stop

  end subroutine filter_reference_load

  module subroutine dumpToHDF_lfi(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
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
    class(comm_lfi_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    character(len=10) :: diode_name
    integer(i4b) :: ierr, i, j
    real(dp), allocatable, dimension(:,:)   :: amp, amp_tot
    real(dp), allocatable, dimension(:,:,:) :: R, R_tot
    real(dp), allocatable, dimension(:,:,:,:) :: ref_filter, adc_corr

    allocate(amp(self%nscan_tot,self%ndet), amp_tot(self%nscan_tot,self%ndet))
    amp = 0.d0
    amp(self%scanid,:) = self%spike_amplitude
    call mpi_reduce(amp, amp_tot, size(amp), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

    if (trim(self%level) == 'L1') then
       allocate(R(self%nscan_tot,self%ndet,size(self%R,3)),R_tot(self%nscan_tot,self%ndet,size(self%R,3)))
       R = 0.d0
       R(self%scanid,:,:) = self%R
       call mpi_reduce(R, R_tot, size(R), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)
    end if

    if (self%myid == 0 .and. trim(self%level) == 'L1') then
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_temp', self%spike_templates)
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_ampl', amp_tot)
       call write_hdf(chainfile, trim(adjustl(path))//'R_factor', R_tot)
       call write_hdf(chainfile, trim(adjustl(path))//'w_diode', self%diode_weights)


    end if

    deallocate(amp, amp_tot)
    if (trim(self%level) == 'L1') deallocate(R, R_tot)

  end subroutine dumpToHDF_lfi

  module subroutine sample_1Hz_spikes(tod, handle, map_sky, m_gain, procmask, procmask2)
    !   Sample LFI specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    implicit none
    class(comm_lfi_tod),                          intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: m_gain
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, k, bin, ierr, nbin
    real(dp)     :: dt, t_tot, t, A, b, mval, eta
    real(dp)     :: t1, t2
    character(len=6) :: scantext
    real(dp), allocatable, dimension(:)     :: nval
    real(sp), allocatable, dimension(:)     :: res
    real(dp), allocatable, dimension(:,:)   :: s_sum
    real(dp), allocatable, dimension(:,:,:) :: s_bin
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '|    --> Sampling 1Hz spikes'

    dt    = 1.d0/tod%samprate   ! Sample time
    t_tot = 1.d0                ! Time range in sec
    nbin  = tod%nbin_spike      ! Number of bins 

    allocate(s_bin(0:nbin-1,tod%ndet,tod%nscan), s_sum(0:nbin-1,tod%ndet), nval(0:nbin-1))

    ! Compute template per scan
    s_bin = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       tod%apply_inst_corr = .false. ! Disable 1Hz correction for just this call
       call sd%init_singlehorn(tod, i, map_sky, m_gain, procmask, procmask2)
       tod%apply_inst_corr = .true.  ! Enable 1Hz correction again

       allocate(res(tod%scans(i)%ntod))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle

          !write(*,*) tod%scanid(i), j, maxval(abs(sd%tod(:,j))), tod%scans(i)%d(j)%gain, maxval(abs(sd%s_sky(:,j))), maxval(abs(sd%s_sl(:,j))), maxval(abs(sd%s_orb(:,j)))
          !res = sd%tod(:,j)/tod%scans(i)%d(j)%gain - (sd%s_sky(:,j) + &
          !     & sd%s_sl(:,j) + sd%s_orb(:,j))
          do k = 1, tod%scans(i)%ntod
             if (sd%tod(k,j) /= sd%tod(k,j)) then
                write(*,*) tod%scanid(i), j, k, sd%tod(k,j), tod%scans(i)%d(j)%gain, sd%s_sky(k,j), sd%s_sl(k,j), sd%s_orb(k,j)
             end if
             !res(k) = 1/tod%scans(i)%d(j)%gain - (sd%s_sky(k,j) + &
             !     & sd%s_sl(k,j) + sd%s_orb(k,j))
             res(k) = sd%tod(k,j)/tod%scans(i)%d(j)%gain - (sd%s_sky(k,j) + &
                  & sd%s_sl(k,j) + sd%s_orb(k,j))
          end do

          nval = 0.d0
          do k = 1, tod%scans(i)%ntod
             if (sd%mask(k,j) == 0.) cycle
             t = modulo(tod%scans(i)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
             bin = min(int(t*nbin),nbin-1)
             s_bin(bin,j,i) = s_bin(bin,j,i)  + res(k)
             nval(bin)      = nval(bin)       + 1.d0
          end do
          if (all(nval > 0)) then
             s_bin(:,j,i) = s_bin(:,j,i) / nval
             s_bin(:,j,i) = s_bin(:,j,i) - mean(s_bin(1:nbin/3,j,i))
          else
             s_bin(:,j,i) = 0.d0
             tod%scans(i)%d(j)%accept = .false.
          end if
       end do

!!$       if (trim(tod%freq) == '070') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_22S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,10,i)
!!$          end do
!!$          close(58)
!!$       end if
!!$
!!$       if (trim(tod%freq) == '044') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_26S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,6,i)
!!$          end do
!!$          close(58)
!!$       end if

       ! Clean up
        call sd%dealloc
        deallocate(res)
    end do

    ! Compute smoothed templates
    s_sum = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       do j = 1, tod%ndet
          s_sum(:,j) = s_sum(:,j) + s_bin(:,j,i)
       end do
    end do
    call mpi_allreduce(mpi_in_place, s_sum,  size(s_sum),  &
         & MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)

    ! Normalize to maximum of unity, and subtract mean
    do j = 1, tod%ndet
       !s_sum(:,j) = s_sum(:,j) - median(s_sum(:,j)) 
       mval = maxval(abs(s_sum(:,j))) 
       do k = 0, nbin-1
          s_sum(k,j) = s_sum(k,j) / mval
       end do
       tod%spike_templates(:,j) = s_sum(:,j) 

       tod%spike_templates(:,j) = tod%spike_templates(:,j) - &
            & sum(tod%spike_templates(:,j))/nbin
    end do

    ! Compute amplitudes per scan and detector
    tod%spike_amplitude = 0.
    do j = 1, tod%ndet
       A = 0.d0; b = 0.d0
       do i = 1, tod%nscan
          if (.not. tod%scans(i)%d(j)%accept) cycle
          b = b + sum(s_sum(:,j)*s_bin(:,j,i)) / tod%scans(i)%d(j)%N_psd%sigma0**2
          A = A + sum(s_sum(:,j)**2)           / tod%scans(i)%d(j)%N_psd%sigma0**2
       end do

       if (tod%info%myid == 0) eta = rand_gauss(handle)
       call mpi_bcast(eta, 1, MPI_DOUBLE_PRECISION, 0, tod%info%comm, ierr)
       call mpi_allreduce(mpi_in_place, A, 1,  &
            & MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)
       call mpi_allreduce(mpi_in_place, b, 1,  &
            & MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)

       if (A == 0.d0) then
          tod%spike_amplitude(:,j) = 0.d0
       else
          tod%spike_amplitude(:,j) = b/A
          if (trim(tod%operation) == 'sample') then
             tod%spike_amplitude(:,j) = tod%spike_amplitude(:,j) + eta / sqrt(A)
          end if
       end if
       !if (tod%info%myid == 0) write(*,*) 'Spike amplitude =', j, tod%spike_amplitude(1,j)
    end do

    ! Clean up
    deallocate(s_bin, s_sum, nval)

  end subroutine sample_1Hz_spikes

  module subroutine construct_corrtemp_lfi(self, scan, pix, psi, s)
    !  Construct an LFI instrument-specific correction template; for now contains 1Hz template only
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_lfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    dt    = 1.d0/self%samprate   ! Sample time
    t_tot = 1.d0                ! Time range in sec
    nbin  = self%nbin_spike      ! Number of bins 

    do j = 1, self%ndet
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t = modulo(self%scans(scan)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
          b = min(int(t*nbin),nbin-1)
          s(k,j) = self%spike_amplitude(scan,j) * self%spike_templates(b,j)
       end do
    end do

  end subroutine construct_corrtemp_lfi


  module subroutine preprocess_L1_to_L2(self, map_sky, procmask)
    implicit none
    class(comm_lfi_tod),                          intent(inout) :: self
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask

    integer(i4b) :: i, j, k, m, n, npix, unit, barrier, mpistat(MPI_STATUS_SIZE), ierr
    character(len=4)   :: id
    character(len=6)   :: scantext
    character(len=512) :: filename
    type(hdf_file) :: h5_file
    real(dp), allocatable, dimension(:,:)     :: m_buf
 !   real(sp), allocatable, dimension(:)       :: procmask
    real(sp), allocatable, dimension(:,:)     :: tod

!    npix = 12*self%nside**2

    ! Distribute processing masks
!!$    allocate(m_buf(0:npix-1,self%nmaps), procmask(0:npix-1))
!!$    call self%procmask%bcast_fullsky_map(m_buf); procmask  = m_buf(:,1)
!!$    deallocate(m_buf)
    
    !unit = getlun()
    !open(unit, file=trim(self%L2file), form='unformatted')


  end subroutine preprocess_L1_to_L2

  module subroutine remove_fixed_scans_lfi(self)
    ! 
    ! Sets accept = .false. for known bad scans
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_lfi_tod),                  intent(inout)  :: self

    integer(i4b) :: i, j, k

    do j = 1, self%ndet
       do i = 1, self%nscan
          k             = self%scanid(i)

          ! Chisquare excess in 70 GHz; unknown origim
          if ((k > 24900 .and. k <= 25300) .and. (trim(self%label(j)) == '18M' .or. trim(self%label(j)) == '18S')) self%scans(i)%d(j)%accept = .false.

          ! 44 GHz triple dot, with weaker effects in the other two channels
          if (k == 6144 .or. k == 6126) self%scans(i)%d(j)%accept = .false.

          ! The Day Planck Stood Still; 14389 has bad chisq
          if (k == 14389 .or. k == 14390) self%scans(i)%d(j)%accept = .false.
  
       end do
    end do


  end subroutine remove_fixed_scans_lfi

end submodule comm_tod_lfi_smod

