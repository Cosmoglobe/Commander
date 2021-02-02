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
module comm_tod_gain_mod
  use comm_tod_mod
  use comm_tod_noise_mod
  use comm_utils
  implicit none



contains


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Haavard: Get rid of explicit n_corr, and replace 1/sigma**2 with proper invN multiplication
!  subroutine sample_gain_per_scan(self, handle, det, scan_id, n_corr, mask, s_ref)
!   subroutine calculate_gain_mean_std_per_scan(self, det, scan_id, s_tot, invn, mask)
   subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invN, mask, s_ref, s_tot, handle)
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: s_invN, mask, s_ref, s_tot
    integer(i4b),                         intent(in)    :: scan_id
    type(planck_rng),                     intent(inout)  :: handle


    real(sp), allocatable, dimension(:,:) :: residual
    real(sp), allocatable, dimension(:)   :: r_fill
    integer(i4b) :: ext(2), j, i, ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_tot,1)

!    allocate(residual(size(stot_invN)))
      !g_old = tod%scans(scan_id)%d(det)%gain 

!   residual = (tod%scans(scan_id)%d(det)%tod - (tod%gain0(0) + &
!         & tod%gain0(det)) * s_tot) * mask

    allocate(r_fill(size(s_tot,1)))
    call tod%downsample_tod(s_tot(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet))
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          residual(:,j) = 0.d0
          cycle
       end if
       r_fill = tod%scans(scan_id)%d(j)%tod - tod%scans(scan_id)%d(j)%baseline & 
         & - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', real(tod%scans(scan_id)%d(j)%sigma0, sp), handle, tod%scans(scan_id)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do
    call multiply_inv_N(tod, scan_id, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    ! Get a proper invn multiplication from Haavard's routine here
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          tod%scans(scan_id)%d(j)%gain  = 0.d0
          tod%scans(scan_id)%d(j)%dgain = 0.d0
       else
          tod%scans(scan_id)%d(j)%dgain      = sum(s_invN(:,j) * residual(:,j))
          tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invN(:,j) * s_ref(:,j))
          if (tod%scans(scan_id)%d(j)%gain_invsigma < 0.d0) then
             write(*,*) 'Warning: Not positive definite invN = ', tod%scanid(scan_id), j, tod%scans(scan_id)%d(j)%gain_invsigma
          end if
       end if
    end do
!    tod%scans(scan_id)%d(det)%dgain      = sum(stot_invN * residual)
!    tod%scans(scan_id)%d(det)%gain_invsigma = sum(stot_invN * s_tot)
!    if (tod%scans(scan_id)%d(det)%gain_invsigma < 0) then
!       write(*,*) 's', sum(mask * stot_invN * s_tot), sum(stot_invN * s_tot)
!    end if

!!$    if
!(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_invsigma > 5.) then
!!$       open(58, file='tod.dat')
!!$       do i = p, q
!!$          if (mask(i) > 0.5) write(58,*) i, residual(i), s_tot(i)*0.001 
!!$       end do
!!$       close(58)
!!$       stop
!!$    end if

    !write(*,*) det, scan_id, real(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_invsigma,sp), real(tod%gain0(0),sp), real(tod%gain0(det),sp), real(g_old,sp)

   ! write(*,*) tod%scanid(scan_id), real(tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), real(tod%gain0(0) + tod%gain0(3) + tod%scans(scan_id)%d(3)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), '# deltagain'

    if (.false. .and. trim(tod%freq) == '030' .and. mod(tod%scanid(scan_id),100) == 0) then
       call int2string(tod%scanid(scan_id), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(1)%gain_invsigma
       open(58,file='gain_delta_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i, residual(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_ref(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_invN(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline &
          & - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline
       end do
       close(58)
    end if


    deallocate(residual, r_fill)

  end subroutine calculate_gain_mean_std_per_scan


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Eirik: Update this routine to sample time-dependent gains properly; results should be stored in self%scans(i)%d(j)%gain, with gain0(0) and gain0(i) included
  subroutine sample_smooth_gain(tod, handle, dipole_mods)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),   dimension(:, :),       intent(in)     :: dipole_mods


!    real(sp), dimension(:, :) :: inv_gain_covar ! To be replaced by proper matrix, and to be input as an argument
    integer(i4b) :: i, j, k, ndet, nscan_tot, ierr, ind(1)
    integer(i4b) :: currstart, currend, window, i1, i2, pid_id, range_end
    real(dp)     :: mu, denom, sum_inv_sigma_squared, sum_weighted_gain, g_tot, g_curr, sigma_curr
    real(dp), allocatable, dimension(:)     :: lhs, rhs, g_smooth
    real(dp), allocatable, dimension(:)     :: temp_gain, temp_invsigsquared
    real(dp), allocatable, dimension(:)     :: summed_invsigsquared, smoothed_gain
    real(dp), allocatable, dimension(:,:,:) :: g
    real(dp), allocatable, dimension(:,:)   :: pidrange_gainarr
    integer(i4b),   allocatable, dimension(:, :) :: pid_ranges
    integer(i4b),   allocatable, dimension(:, :) :: window_sizes
    integer(i4b), save :: count = 0
    character(len=128)  :: kernel_type

    ndet       = tod%ndet
    nscan_tot  = tod%nscan_tot

    ! Collect all gain estimates on the root processor
    allocate(g(nscan_tot,ndet,2))
    allocate(lhs(nscan_tot), rhs(nscan_tot))
    g = 0.d0
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          g(k,j,1) = tod%scans(i)%d(j)%dgain
          g(k,j,2) = tod%scans(i)%d(j)%gain_invsigma
       end do
    end do
    if (tod%myid == 0) then
       call mpi_reduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, tod%comm, ierr)
    else
       call mpi_reduce(g,            g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, tod%comm, ierr)
    end if

    if (tod%myid == 0) then
!!$       open(58,file='tmp.unf', form='unformatted')
!!$       read(58) g
!!$       close(58)

        count = count+1
        !write(*, *) "FREQ IS ", trim(tod%freq), count
       !nbin = nscan_tot / binsize + 1

        !open(58,file='gain_' // trim(tod%freq) // '.dat', recl=1024)
       do j = 1, ndet
          do k = 1, nscan_tot
             !if (g(k,j,2) /= 0) then
             !if (g(k,j,2) /= g(k,j,2)) write(*,*) j,k, real(g(k,j,1),sp), real(g(k,j,2),sp), real(g(k,j,1)/g(k,j,2),sp)
             if (g(k,j,2) > 0) then
                if (abs(g(k, j, 1)) > 1e10) then
                   write(*, *) 'G1'
                   write(*, *) g(k, j, 1)
                end if
                if (abs(g(k, j, 2)) > 1e10) then
                   write(*, *) 'G2'
                   write(*, *) g(k, j, 2)
                end if
                !if (abs(dipole_mods(k, j) > 1e10)) then
                !   write(*, *) 'DIPOLE_MODS'
                !   write(*, *) dipole_mods(k, j)
                !else
                  ! write(58,*) j, k, real(g(k,j,1)/g(k,j,2),sp), real(g(k,j,1),sp), real(g(k,j,2),sp), real(dipole_mods(k, j), sp)
                !end if
             else
                !write(58,*) j, k, 0., 0.0, 0., 0.
             end if
          end do
          !write(58,*)
       end do
       !close(58)

       allocate(window_sizes(tod%ndet, tod%nscan_tot))
       call get_smoothing_windows(tod, window_sizes, dipole_mods)
!!$       open(58, file='smoothing_windows' // trim(tod%freq) // '.dat', recl=1024)
!!$       do j = 1, ndet
!!$         do k = 1, size(window_sizes(j, :))
!!$            if (window_sizes(j, k) == 0) cycle
!!$            write(58, *) j, k, window_sizes(j, k)
!!$         end do
!!$       end do
!!$       close(58)

       allocate(pidrange_gainarr(nscan_tot, ndet))
       pidrange_gainarr = 0.d0
       do j = 1, ndet
         do k = 1, nscan_tot
            if (g(k, j, 2) > 0.d0) then
               pidrange_gainarr(k, j) = g(k, j, 1) / g(k, j, 2)
            end if
         end do
       end do
!       call get_pid_ranges(pid_ranges, tod, pidrange_gainarr, dipole_mods, window_sizes)
       call get_pid_ranges_tabulated(pid_ranges, tod)
       deallocate(pidrange_gainarr)
!!$       open(58, file='pid_ranges_' // trim(tod%freq) // '.dat', recl=1024)
!!$       do j = 1, ndet
!!$         do k = 1, size(pid_ranges(j, :))
!!$            if (pid_ranges(j, k) == 0) exit
!!$            write(58, *) j, k, pid_ranges(j, k)
!!$         end do
!!$       end do
!!$       close(58)

       do j = 1, ndet
         lhs = 0.d0
         rhs = 0.d0
         pid_id = 1
         k = 0
         !write(*,*) "PIDRANGE: ", pid_ranges(j, :)
         do while (pid_id < size(pid_ranges(j, :)))
            if (pid_ranges(j, pid_id) == 0) exit
            currstart = pid_ranges(j, pid_id)
            if (pid_ranges(j, pid_id+1) == 0) then
               currend = nscan_tot
            else
               !currend = pid_ranges(j, pid_id +1)
               currend = pid_ranges(j, pid_id +1) - 1
            end if
            sum_weighted_gain = 0.d0
            sum_inv_sigma_squared = 0.d0
            allocate(temp_gain(currend - currstart + 1))
            allocate(temp_invsigsquared(currend - currstart + 1))
            allocate(summed_invsigsquared(currend-currstart + 1))
            allocate(smoothed_gain(currend-currstart + 1))
            do k = currstart, currend
               !if (count == 12) write(*,*) j, k, g(k, j, 1) , g(k, j, 2)
               if (g(k,j,2) /= g(k,j,2)) then
                  write(*,*) 'GAIN IS NAN', k, j
                  temp_gain(k-currstart + 1) = 0.d0
               else if (g(k,j,2) > 0.d0) then
                  temp_gain(k-currstart + 1) = g(k, j, 1) / g(k, j, 2)
                  if (trim(tod%operation) == 'sample') then
                     temp_gain(k-currstart+1) = temp_gain(k-currstart+1) + rand_gauss(handle) / g(k, j, 2)
                  end if
               else
                  temp_gain(k-currstart + 1) = 0.d0
               end if
               temp_invsigsquared(k - currstart + 1) = max(g(k, j, 2),0.d0)
            end do
!            write(*, *) 'FREQ:', trim(tod%freq)
!            write(*, *) 'TEMP_INVSIGSQUARED:', temp_invsigsquared
!!$            write(*,*) 't', j, currstart, currend
!!$            if (j == 4 .and. currstart == 8684) then
!!$               write(*,*) temp_gain
!!$               write(*,*) temp_invsigsquared
!!$            end if
!!$            call moving_average_padded_variable_window(temp_gain, smoothed_gain, &
            kernel_type = 'boxcar'
            call moving_average_variable_window(temp_gain, smoothed_gain, &
               & window_sizes(j, currstart:currend), temp_invsigsquared, summed_invsigsquared, kernel_type)
            if (any(summed_invsigsquared < 0)) then
               write(*, *) 'WHOOOOPS'
               write(*, *) 'currstart', currstart
               write(*, *) 'currend', currend
               write(*, *) 'temp_invsigsquared', temp_invsigsquared
               stop
            end if
            !write(*, *) 'SMOOTHED_GAIN:', smoothed_gain
            !write(*, *) 'SUMMED_INVSIGSQUARED:', summed_invsigsquared
            do k = currstart, currend
               g(k, j, 1) = smoothed_gain(k - currstart + 1)
!               if (trim(tod%operation) == 'sample' .and. summed_invsigsquared(k-currstart+1) > 0.d0) then
!                  g(k, j, 1) = smoothed_gain(k - currstart + 1) + &
!                     & rand_gauss(handle) / &
!                     & sqrt(summed_invsigsquared(k - currstart + 1))
!               else
!                  g(k, j, 1) = smoothed_gain(k - currstart + 1)
!               end if
               if (summed_invsigsquared(k - currstart + 1) > 0) then
                  g(k, j, 2) = 1.d0 / sqrt(summed_invsigsquared(k - currstart + 1))
               else
                  g(k, j, 2) = 0.d0
               end if
            end do
            pid_id = pid_id + 1

            deallocate(temp_gain)
            deallocate(temp_invsigsquared)
            deallocate(summed_invsigsquared)
            deallocate(smoothed_gain)
         end do

         ! Doing binning for now, but should use proper covariance sampling
!         do i = 1, nbin
!            b1 = (i-1) * binsize + 1
!            b2 = min(i * binsize, nscan_tot)
!            if (count(g(b1:b2, j, 2) > 0) <= 1) cycle
!            n = 0
!            do k = b1, b2
!               if (g(k,j,2) <= 0.d0) then
!                  g(k,j,1) = 0.d0
!                  cycle
!               end if
!               if (g(k,j,1) == 0) cycle
!               n = n + 1
!               ! To be replaced by proper matrix multiplications
!               lhs(k) = g(k, j, 2)
!   !            rhs(k) = rhs(k) + g(k, j, 1) + sqrt(inv_gain_covar(j, k)) * rand_gauss(handle) + sqrt(g(k, j, 2)) * rand_gauss(handle) this is proper when we have the covariance matrix
!               ! HKE: Disabling fluctuations for now 
!               rhs(k) = g(k, j, 1) !+ sqrt(g(k, j, 2)) * rand_gauss(handle)
!!               g(k, j, 1) = rhs(k) / lhs(k)
!               vals(n) = rhs(k) / lhs(k)
!               if (abs(vals(n)) > 0.05d0) vals(n) = 0.d0
!            end do
!            where (g(b1:b2, j, 2) > 0)
!               g(b1:b2, j, 1) = median(vals(1:n))
!            end where
!         end do

         ! Apply running average smoothing
!         allocate(g_smooth(nscan_tot))
!         window = 500
!         do k = 1, nscan_tot
!            i1 = max(k-window,1)
!            i2 = min(k+window,nscan_tot)
!            g_smooth(k) = sum(g(i1:i2,j,1)) / (i2-i1)
!            !g_smooth(k) = median(g(i1:i2,j,1)) 
!         end do
!         g(:,j,1) = g_smooth - mean(g_smooth)
!         deallocate(g_smooth)
 
!!$         ! Apply running average smoothing
!!$         allocate(g_smooth(nscan_tot))
!!$         window = 10000
!!$         do k = 1, nscan_tot
!!$            i1 = max(k-window,1)
!!$            i2 = min(k+window,nscan_tot)
!!$            ind = minloc(g(i1:i2,j,2))
!!$            g_smooth(k) = g(ind(1),j,1)
!!$            !g_smooth(k) = median(g(i1:i2,j,1)) 
!!$         end do
!!$         g(:,j,1) = g_smooth - mean(g_smooth)
!!$         deallocate(g_smooth)

         mu  = 0.d0
         denom = 0.d0
         do k = 1, nscan_tot
            if (g(k, j, 2) <= 0.d0) cycle
            mu         = mu + g(k, j, 1)! * g(k,j,2)
            denom      = denom + 1.d0! * g(k,j,2)
!            write(*,*) j, k, g(k,j,1), rhs(k)/lhs(k)
         end do
         mu = mu / denom

         ! Make sure fluctuations sum up to zero
         if (tod%verbosity > 1) then
           write(*,*) 'mu = ', mu
         end if
         g(:,j,1) = g(:,j,1) - mu
       end do
!       open(58,file='gain_postsmooth' // trim(tod%freq) // '.dat', recl=1024)
       do j = 1, ndet
          do k = 1, nscan_tot
             !if (g(k,j,2) /= 0) then
             if (g(k,j,2) > 0) then
                if (abs(g(k, j, 1)) > 1e10) then
                   write(*, *) 'G1_postsmooth'
                   write(*, *) g(k, j, 1)
                end if
                if (abs(g(k, j, 2)) > 1e10) then
                   write(*, *) 'G2_postsmooth'
                   write(*, *) g(k, j, 2)
                end if
                if (abs(dipole_mods(k, j) > 1e10)) then
                   write(*, *) 'DIPOLE_MODS'
                   write(*, *) dipole_mods(k, j)
                else
!                   write(58,*) j, k, real(g(k,j,1)/g(k,j,2),sp), real(g(k,j,1),sp), real(g(k,j,2),sp), real(dipole_mods(k, j), sp)
                end if
             else
!                write(58,*) j, k, 0., 0.0, 0., 0.
             end if
          end do
!          write(58,*)
       end do
       close(58)

       deallocate(window_sizes, pid_ranges)
    end if

    ! Distribute and update results
    call mpi_bcast(g, size(g),  MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)    
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          !if (g(k, j, 2) <= 0.d0) cycle
          tod%scans(i)%d(j)%dgain = g(k,j,1)
          tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + g(k,j,1)
          !write(*,*) j, k,  tod%scans(i)%d(j)%gain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(g, lhs, rhs)

  end subroutine sample_smooth_gain


  ! Haavard: Remove current monopole fit, and replace inv_sigmasq with a proper invN(alpha,fknee) multiplication
!  subroutine accumulate_abscal(self, scan, det, mask, s_sub, &
!       & s_orb, A_abs, b_abs)
   ! This is implementing equation 16, adding up all the terms over all the sums
   ! the sum i is over the detector.
   subroutine accumulate_abscal(tod, scan, mask, s_sub, s_ref, s_invN, A_abs, b_abs, handle, out, s_highres)
    implicit none
    class(comm_tod),                   intent(in)     :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub, s_ref
    real(sp),          dimension(:,:), intent(in)     :: s_invN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    logical(lgt), intent(in) :: out
    real(sp),          dimension(:,:), intent(in), optional :: s_highres
 
    real(sp), allocatable, dimension(:,:)     :: residual
    real(sp), allocatable, dimension(:)       :: r_fill
    real(dp)     :: A, b, scale
    integer(i4b) :: i, j, ext(2), ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_sub,1)
!    p = 5000 ! Buffer width
!    q = size(mask)-p
    call tod%downsample_tod(s_sub(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet), r_fill(size(s_sub,1)))
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) then
          residual(:,j) = 0.
          cycle
       end if
       r_fill = tod%scans(scan)%d(j)%tod-s_sub(:,j) - tod%scans(scan)%d(j)%baseline
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', abs(real(tod%scans(scan)%d(j)%sigma0, sp)), handle, tod%scans(scan)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j))
       b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j))
!       if (out) write(*,*) tod%scanid(scan), real(sum(s_invN(:,j) * residual(:,j))/sum(s_invN(:,j) * s_ref(:,j)),sp), real(1/sqrt(sum(s_invN(:,j) * s_ref(:,j))),sp), '  # absK', j
    end do

!    if (trim(tod%freq) == '070') then
!       write(*,*) tod%scanid(scan), real(b/A,sp), real(1/sqrt(A),sp), '  # absK', det
 !   end if


    if (.false. .and. mod(tod%scanid(scan),1000) == 0 .and. out) then
       call int2string(tod%scanid(scan), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%gain0(0), tod%gain0(1)
       open(58,file='gainfit3_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i-ext(1)+1, residual(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_ref(i,1)
       end do
       write(58,*)
       scale = sum(s_invN(:,1) * residual(:,1)) / sum(s_invN(:,1) * s_ref(:,1))
       residual(:,1) = residual(:,1) - scale * s_ref(:,1)
       do i = 1, size(s_ref,1)
          write(58,*) i-ext(1)+1, residual(i,1) 
       end do
       close(58)

       open(58,file='gainfit4_'//itext//'.dat')       
       do i = 1, size(s_sub,1)
          write(58,*) i, tod%scans(scan)%d(4)%tod(i) - tod%scans(scan)%d(4)%baseline
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, r_fill(i)
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, s_highres(i,4)
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, s_sub(i,4)
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, mask(i,4)*10
       end do
       close(58)
    end if



    !if (det == 1) write(*,*) tod%scanid(scan), b/A, '  # absgain', real(A,sp), real(b,sp)

    deallocate(residual, r_fill)

  end subroutine accumulate_abscal

  ! Sample absolute gain from orbital dipole alone 
  subroutine sample_abscal_from_orbital(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%gain0(0) = sum(b)/sum(A)
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%gain0(0) = tod%gain0(0) + 1.d0/sqrt(sum(A)) * rand_gauss(handle)
       end if
       if (tod%verbosity > 1) then
         write(*,*) 'abscal = ', tod%gain0(0), sum(b), sum(A)
       end if
    end if
    call mpi_bcast(tod%gain0(0), 1,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(A, b)

  end subroutine sample_abscal_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  ! Eirik: Change this routine to sample the constant offset per radiometer; put the result in tod%gain0(i)
  subroutine sample_relcal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b, rhs, x
    real(dp), allocatable, dimension(:, :) :: coeff_matrix

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet), rhs(tod%ndet+1), x(tod%ndet+1))
    allocate(coeff_matrix(tod%ndet+1, tod%ndet+1))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    coeff_matrix = 0.d0
    rhs = 0.d0
    if (tod%myid == 0) then
       do j = 1, tod%ndet
         coeff_matrix(j, j) = A(j)
         rhs(j) = b(j) 
         if (trim(tod%operation) == 'sample') rhs(j) = rhs(j) + sqrt(A(j)) * rand_gauss(handle)
         coeff_matrix(j, tod%ndet+1) = 0.5d0
         coeff_matrix(tod%ndet+1, j) = 1
       end do
       coeff_matrix(tod%ndet+1, tod%ndet+1) = 0.d0
       rhs(tod%ndet+1) = 0.d0
       call solve_system_real(coeff_matrix, x, rhs)
       if (tod%verbosity > 1) then
         write(*,*) 'relcal = ', real(x,sp)
       end if
    end if
    call mpi_bcast(x, tod%ndet+1, MPI_DOUBLE_PRECISION, 0, &
       & tod%info%comm, ierr)

    tod%gain0(1:tod%ndet) = x(1:tod%ndet)
    deallocate(coeff_matrix, rhs, A, b, x)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

  end subroutine sample_relcal


  subroutine get_pid_ranges(pid_ranges, tod, dgains, dipole_mods, window_sizes)
     implicit none

     integer(i4b), allocatable, dimension(:, :), intent(out)     :: pid_ranges
     class(comm_tod),           intent(in)          :: tod
     real(dp), dimension(:, :), intent(in)          :: dipole_mods
     real(dp), dimension(:, :), intent(in)          :: dgains
     integer(i4b), dimension(:, :), intent(in)      :: window_sizes

     integer(i4b), allocatable, dimension(:)        :: jump_indices
     integer(i4b), allocatable, dimension(:)        :: sorted_indices
     real(dp),  allocatable,    dimension(:)        :: smoothed_data, smoothed_vars
     real(dp)                                       :: target_percentile
     real(dp)                                       :: quantile_quantum
     real(dp)                                       :: prev_jump_var
     real(sp)                                       :: jump_percentile
     integer(i4b)                                   :: i, j, k, percentile_index
     integer(i4b)                                   :: slow_smooth_window_size
     integer(i4b)                                   :: max_n_jumps
     integer(i4b)                                   :: nscan, range_idx
     integer(i4b)                                   :: start_idx, end_idx, pos
     logical(lgt)                                   :: in_high_var_region


     select case (trim(tod%freq))
         case ('030')
            slow_smooth_window_size = 300
            jump_percentile = 0.99
         case ('044')
            slow_smooth_window_size = 200
            jump_percentile = 0.999
         case ('070')
            slow_smooth_window_size = 150
            jump_percentile = 0.995
         ! Currently, the 70-GHz ds
         case default
            slow_smooth_window_size = 150
            jump_percentile = 0.995
      end select

      nscan = tod%nscan_tot

      quantile_quantum = 1.d0/nscan
      max_n_jumps = ceiling((1.d0 - jump_percentile) / quantile_quantum)
      percentile_index = floor(nscan * jump_percentile)
      allocate(jump_indices(nscan))
      allocate(sorted_indices(nscan))
      allocate(pid_ranges(tod%ndet, max_n_jumps+2))
      allocate(smoothed_data(nscan), smoothed_vars(nscan))

      do i = 1, nscan
         sorted_indices(i) = i
      end do

      do i = 1, tod%ndet
         smoothed_data = 0.d0
         smoothed_vars = 0.d0
!!$         open(58, file='gain_notmovaverage_' // trim(tod%label(i)) // '.dat', recl=1024)
!!$         do j = 1, size(dgains(:, i))
!!$            write(58, *) dgains(j, i)
!!$         end do
!!$         close(58)
         call moving_average(dgains(:, i), smoothed_data, slow_smooth_window_size, &
            weights=dipole_mods(:, i))
!!$         open(58, file='gain_movaverage_' // trim(tod%label(i)) // '.dat', recl=1024)
!!$         do j = 1, size(smoothed_data)
!!$            write(58, *) smoothed_data(j)
!!$         end do
!!$         close(58)
         call moving_variance(smoothed_data, smoothed_vars, slow_smooth_window_size)
!!$         open(58, file='gain_variance_' // trim(tod%label(i)) // '.dat', recl=1024)
!!$         do j = 1, size(smoothed_vars)
!!$            write(58, *) smoothed_vars(j)
!!$         end do
!!$         close(58)
         smoothed_vars = smoothed_vars * dipole_mods(:, i)
         ! Just reusing array as buffer instead of allocating a whole new one
         smoothed_data = smoothed_vars
         call Quicksort(sorted_indices, smoothed_data)
         target_percentile = smoothed_data(percentile_index)
         j = 1
         range_idx = 1
         pid_ranges(i, :) = 0
         pid_ranges(i, 1) = 1
         in_high_var_region = .false.
         prev_jump_var = 1d30
         do while (j <= nscan)
            if (smoothed_vars(j) > target_percentile .and. .not. in_high_var_region) then
               !write(*, *) 'In high var'
               in_high_var_region = .true.
               start_idx = j
               !write(*, *) 'start_idx: ', start_idx
            else if (in_high_var_region .and. & 
               & smoothed_vars(j) <= target_percentile .and. & 
               & smoothed_vars(j) /= 0) then
               !write(*, *) 'End high var'
               in_high_var_region = .false.
               end_idx = j
               !write(*, *) 'end_idx: ', end_idx
               pos = maxloc(smoothed_vars(start_idx:end_idx-1), dim=1) + start_idx
               !write(*, *) 'pos: ', pos
               !write(*, *) 'window_size: ', window_sizes(i, pos)
               !write(*, *) 'prev_pid_range: ', pid_ranges(i, range_idx)
               if ((nscan - pos) < window_sizes(i, pos)) then
                  j = j + 1
                  !write(*, *) 'Cycle 1'
                  cycle
               else if (pos - pid_ranges(i, range_idx) < window_sizes(i, pos)) then
                  ! If this proposed jump has a greater variance, choose it
                  ! instead of the previous one
                  if (prev_jump_var < smoothed_vars(pos)) then
                     pid_ranges(i, range_idx) = pos
                     prev_jump_var = smoothed_vars(pos)
                  end if
                  j = j + 1
                  !write(*, *) 'Cycle 2'
                  cycle
               end if
               !write(*, *) 'Not cycling'
               prev_jump_var = smoothed_vars(pos)
               range_idx = range_idx + 1
               pid_ranges(i, range_idx) = pos
            end if
            j = j + 1
         end do
         range_idx = range_idx + 1
      end do
      deallocate(jump_indices, sorted_indices, smoothed_data, smoothed_vars)

  end subroutine get_pid_ranges


  subroutine get_pid_ranges_tabulated(pid_ranges, tod)
     ! From NPIPE's gain jumps
     implicit none
     integer(i4b), allocatable, dimension(:, :), intent(out)     :: pid_ranges
     class(comm_tod),           intent(in)          :: tod

     integer(i4b)           :: n_jumps

!     n_jumps = 17 ! Npipe has 15 events + the beginning and end
     n_jumps = 26 ! Npipe has 15 events + the beginning and end (but two of them are too bunched up)

     allocate(pid_ranges(tod%ndet, n_jumps))
     pid_ranges(:, :) = 0

!!$     pid_ranges(:, 1) = 1
!!$     pid_ranges(:, 2) = 3352
!!$     pid_ranges(:, 3) = 5030
!!$     pid_ranges(:, 4) = 5484
!!$     pid_ranges(:, 5) = 10911
!!$     pid_ranges(:, 6) = 15957
!!$     pid_ranges(:, 7) = 16455
!!$     pid_ranges(:, 8) = 21484
!!$     pid_ranges(:, 9) = 25654
!!$     pid_ranges(:, 10) = 27110
!!$     pid_ranges(:, 11) = 27343
!!$     pid_ranges(:, 12) = 30387
!!$     pid_ranges(:, 13) = 32763
!!$     pid_ranges(:, 14) = 38591
!!$     pid_ranges(:, 15) = 43929
!!$!     pid_ranges(:, 16) = 44063
!!$     ! This last event is too close to the previous one
!!$     pid_ranges(:, 16) = 0

     pid_ranges(:, 1) = 1
     pid_ranges(:, 2) = 3352
     pid_ranges(:, 3) = 5030
     pid_ranges(:, 4) = 5484
     pid_ranges(:, 5) = 8309 ! 20K
     pid_ranges(:, 6) = 8503 ! 20K
     pid_ranges(:, 7) = 8606 ! 20K
     pid_ranges(:, 8) = 9613 ! 20K
     pid_ranges(:, 9) = 10117 ! 20K
     pid_ranges(:, 10) = 10512 ! 20K
     ! There's one more 20K at 10897 but that's very close to this one
     pid_ranges(:, 11) = 10911
     pid_ranges(:, 12) = 14061 ! 20K
     pid_ranges(:, 13) = 15957
     pid_ranges(:, 14) = 16204 ! 20K
     pid_ranges(:, 15) = 16455
     pid_ranges(:, 16) = 17024 ! 20K
     pid_ranges(:, 17) = 18338 ! 20K
     pid_ranges(:, 18) = 21484
     pid_ranges(:, 19) = 25654
     pid_ranges(:, 20) = 27110
     pid_ranges(:, 21) = 27343
     pid_ranges(:, 22) = 30387
     pid_ranges(:, 23) = 32763
     pid_ranges(:, 24) = 38591
     pid_ranges(:, 25) = 43929
     ! pid_ranges(:, 16) = 44063
     ! This last event is too close to the previous one
     pid_ranges(:, 26) = 0
     
  end subroutine get_pid_ranges_tabulated

  subroutine get_smoothing_windows(tod, windows, dipole_mods)
     implicit none

     integer(i4b), dimension(:, :), intent(out)         :: windows
     class(comm_tod),               intent(in)          :: tod
     real(dp), dimension(:, :),     intent(in)          :: dipole_mods

     real(sp), allocatable, dimension(:)                :: nonzero_mask
     integer(i4b)  :: window_size_dipole_minimum, window_size_dipole_maximum
     real(dp)   :: low_dipole_thresh, curr_dipole
     real(dp)   :: high_dipole_thresh

     real(dp)   :: reduce_fac 

     integer(i4b)   :: i, j, window_size


     high_dipole_thresh = 5d-6
     low_dipole_thresh = 1d-6
     reduce_fac = 1.d0

     select case (trim(tod%freq))
         case ('030')
            window_size_dipole_minimum = 1200
            window_size_dipole_maximum = 400
         case ('044')
            window_size_dipole_minimum = 1500
            window_size_dipole_maximum = 300
         case ('070')
            window_size_dipole_minimum = 1880
            window_size_dipole_maximum = 400
         case ('023-WMAP_K')
             ! I think this is the width of the number of scans to smooth the
             ! gain estimate over. Could be worth experimenting on this one.
             window_size_dipole_minimum = 120
             window_size_dipole_maximum = 40
         ! Currently, the 70-GHz ds
         case default
            window_size_dipole_minimum = 1800
            window_size_dipole_maximum = 400
      end select

      do i = 1, tod%nscan_tot
         do j = 1, tod%ndet
            curr_dipole = dipole_mods(i, j)
            if (curr_dipole == 0) then
               windows(j, i) = 0
            else if (curr_dipole < low_dipole_thresh) then
               windows(j, i) = window_size_dipole_minimum
            else if (curr_dipole > high_dipole_thresh) then
               windows(j, i) = window_size_dipole_maximum
            else
               ! Interpolate
               windows(j, i) = window_size_dipole_minimum + &
                  & int(real(window_size_dipole_maximum - & 
                  & window_size_dipole_minimum, dp) / (high_dipole_thresh - &
                  & low_dipole_thresh) * (curr_dipole - low_dipole_thresh))
            end if
            windows(j, i) = windows(j, i) * reduce_fac
         end do
      end do

  end subroutine get_smoothing_windows




end module comm_tod_gain_mod
