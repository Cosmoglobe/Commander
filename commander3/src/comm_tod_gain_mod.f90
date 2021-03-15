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
  !  subroutine sample_gain_per_scan(self, handle, det, scan_id, n_corr, mask, s_ref)
  !  subroutine calculate_gain_mean_std_per_scan(self, det, scan_id, s_tot, invn, mask)
   subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invN, mask, s_ref, s_tot, handle, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: s_invN, mask, s_ref, s_tot
    integer(i4b),                         intent(in)    :: scan_id
    type(planck_rng),                     intent(inout)  :: handle
    real(sp),             dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),             dimension(:,:), intent(in), optional :: tod_arr


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
       if (present(tod_arr)) then
         r_fill = tod_arr(:, j) - tod%scans(scan_id)%d(j)%baseline & 
           & - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       else
         r_fill = tod%scans(scan_id)%d(j)%tod - tod%scans(scan_id)%d(j)%baseline & 
           & - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       end if
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
          if (present(mask_lowres)) then
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invN(:,j) * residual(:,j) * mask_lowres(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invN(:,j) * s_ref(:,j)    * mask_lowres(:,j))
          else
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invN(:,j) * residual(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invN(:,j) * s_ref(:,j))
          end if
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

    !if (.false. .and. trim(tod%freq) == '030' .and. mod(tod%scanid(scan_id),100) == 0) then
    if (.false.) then
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
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 1) - tod%scans(scan_id)%d(1)%baseline &
            & - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
          else
            write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline &
            & - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
          end if
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 1) - tod%scans(scan_id)%d(1)%baseline
          else
            write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline
          end if
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

        open(58,file='gain_' // trim(tod%freq) // '.dat', recl=1024)
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
                   write(58,*) j, k, real(g(k,j,1)/g(k,j,2),sp), real(g(k,j,1),sp), real(g(k,j,2),sp), real(dipole_mods(k, j), sp)
                !end if
             else
                write(58,*) j, k, 0., 0.0, 0., 0.
             end if
          end do
          write(58,*)
       end do
       close(58)

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


       do j = 1, ndet
         lhs = 0.d0
         rhs = 0.d0
         pid_id = 1
         k = 0
         !write(*,*) "PIDRANGE: ", tod%jumplist(j, :)
         do while (pid_id < size(tod%jumplist(j, :)))
            if (tod%jumplist(j, pid_id) == 0) exit
            currstart = tod%jumplist(j, pid_id)
            if (tod%jumplist(j, pid_id+1) == 0) then
               currend = nscan_tot
            else
               !currend = tod%jumplist(j, pid_id +1)
               currend = tod%jumplist(j, pid_id +1) - 1
            end if
            sum_weighted_gain = 0.d0
            sum_inv_sigma_squared = 0.d0
            allocate(temp_gain(currend - currstart + 1))
            allocate(temp_invsigsquared(currend - currstart + 1))
            allocate(summed_invsigsquared(currend-currstart + 1))
            allocate(smoothed_gain(currend-currstart + 1))
            do k = currstart, currend
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
!         if (tod%verbosity > 1) then
!           write(*,*) 'mu = ', mu
!         end if
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

       deallocate(window_sizes)
    end if

    ! Distribute and update results
    call mpi_bcast(g, size(g),  MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)    
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          !if (g(k, j, 2) <= 0.d0) cycle
          tod%scans(i)%d(j)%dgain = g(k,j,1)
          tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + g(k,j,1)
          !write(*,*) j, k,  tod%gain0(0), tod%gain0(j), g(k,j,1), tod%scans(i)%d(j)%gain 
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
   subroutine accumulate_abscal(tod, scan, mask, s_sub, s_ref, s_invN, A_abs, b_abs, handle, out, s_highres, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                   intent(in)     :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub, s_ref
    real(sp),          dimension(:,:), intent(in)     :: s_invN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    logical(lgt), intent(in) :: out
    real(sp),          dimension(:,:), intent(in), optional :: s_highres
    real(sp),          dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),          dimension(:,:), intent(in), optional :: tod_arr
 
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
       if (present(tod_arr)) then
         r_fill = tod_arr(:,j)-s_sub(:,j)! - tod%scans(scan)%d(j)%baseline
         if (tod%scanid(scan) == 30 .and. out) write(*,*) tod%scanid(scan), sum(abs(tod_arr(:,j))), sum(abs(s_sub(:,j))), tod%scans(scan)%d(j)%baseline
       else
         r_fill = tod%scans(scan)%d(j)%tod - s_sub(:,j)! - tod%scans(scan)%d(j)%baseline
         if (tod%scanid(scan) == 30 .and. out) write(*,*) tod%scanid(scan), sum(abs(tod%scans(scan)%d(j)%tod)), sum(abs(s_sub(:,j))), tod%scans(scan)%d(j)%baseline
       end if
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', abs(real(tod%scans(scan)%d(j)%sigma0, sp)), handle, tod%scans(scan)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (present(mask_lowres)) then
          A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j)    * mask_lowres(:,j))
          b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j) * mask_lowres(:,j))
       else
          A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j))
          b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j))
       end if
       if (tod%scanid(scan) == 30 .and. out) then
         write(*,*) 'scan, s N^-1 r/s N^-1 s, sigma, absK, det, sum(abs(s_invN)), sum(abs(s_ref)), sum(abs(mask)), sum(abs(res), sum(s N^-1 ref*mask), sum(s N^-1 res*mask)'
         write(*,*) tod%scanid(scan), real(sum(s_invN(:,j) * residual(:,j))/sum(s_invN(:,j) * s_ref(:,j)),sp), real(1/sqrt(sum(s_invN(:,j) * s_ref(:,j))),sp), '  # absK', j, sum(abs(s_invN(:,j))), sum(abs(s_ref(:,j))), sum(abs( mask_lowres(:,j))), sum(abs(residual(:,j))), sum(s_invN(:,j) * s_ref(:,j)    * mask_lowres(:,j)), sum(s_invN(:,j) * residual(:,j) * mask_lowres(:,j))
       end if
    end do

    if (.false. .and. out) then
    !if (sum(s_invN(:,4) * residual(:,4) * mask_lowres(:,4))/sum(s_invN(:,4) * s_ref(:,4) * mask_lowres(:,4)) < -0.5d0) then
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
          write(58,*) i+ext(1)-1, residual(i+ext(1)-1,1) 
       end do
       if (present(mask_lowres)) then
          write(58,*)
          do i = 1, size(s_ref,1)
             write(58,*) i+ext(1)-1, mask_lowres(i,1) 
          end do
       end if
       close(58)

       open(58,file='gainfit4_'//itext//'.dat')       
       do i = 1, size(s_sub,1)
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 4) - tod%scans(scan)%d(4)%baseline
          else
            write(58,*) i, tod%scans(scan)%d(4)%tod(i) - tod%scans(scan)%d(4)%baseline
          end if
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, r_fill(i)
       end do
       write(58,*)
       if (present(s_highres)) then
          do i = 1, size(s_sub,1)
             write(58,*) i, s_highres(i,4)
          end do
          write(58,*)
       end if
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
         write(*,*) 'abscal = ', tod%gain0(0)
         write(*,*) 'sum(b), sum(A) = ', sum(b), sum(A)
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

  subroutine sample_imbal_cal(tod, handle, A_abs, b_abs)
    !  Subroutine to sample the transmission imbalance parameters, defined in
    !  the WMAP data model as the terms x_im; given the definition
    !  d_{A/B} = T_{A/B} \pm Q_{A/B} cos(2 gamma_{A/B}) \pm U_{A/B} sin(2 gamma_{A/B})
    !  we have
    !  d = g[(1+x_im)*d_A - (1-x_im)*d_B]
    !  Returns x_{im,1} for detectors 13/14, and x_{im,2} for detectors 23/24.
    !
    !
    !  Arguments (fixed):
    !  ------------------
    !  A_abs: real(dp)
    !     Accumulated A_abs = s_ref^T N^-1 s_ref for all scans
    !  b_abs: real(dp)
    !     Accumulated b_abs = s_ref^T N^-1 s_sub for all scans
    !
    !  
    !  Arguments (modified):
    !  ---------------------
    !  tod: comm_WMAP_tod
    !     The entire tod object. tod%x_im estimated and optionally sampled
    !  handle: planck_rng derived type 
    !     Healpix definition for random number generation
    !     so that the same sequence can be resumed later on from that same point
    !
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
       tod%x_im(1) = sum(b(1:2))/sum(A(1:2))
       tod%x_im(3) = sum(b(3:4))/sum(A(3:4))
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%x_im(1) = tod%x_im(1) + 1.d0/sqrt(sum(A(1:2))) * rand_gauss(handle)
          tod%x_im(3) = tod%x_im(3) + 1.d0/sqrt(sum(A(3:4))) * rand_gauss(handle)
       end if
       tod%x_im(2) = tod%x_im(1)
       tod%x_im(4) = tod%x_im(3)
       if (tod%verbosity > 1) then
         write(*,*) 'imbal =', tod%x_im(1), tod%x_im(3)
       end if
    end if
    call mpi_bcast(tod%x_im, 4,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)


    deallocate(A, b)

  end subroutine sample_imbal_cal

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
             window_size_dipole_minimum = 1 !120
             window_size_dipole_maximum = 1 ! 40
         case default
            window_size_dipole_minimum = 1 !1800
            window_size_dipole_maximum = 1 !400
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
