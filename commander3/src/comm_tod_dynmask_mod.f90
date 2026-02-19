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

module comm_tod_dynmask_mod
  use comm_status_mod
  use comm_tod_mod
  use comm_param_mod
  implicit none
  
  private
  public comm_dynmask

  integer(i4b), parameter :: N_rms    = 3
  integer(i4b), parameter :: flag_dyn = 2**30

  ! This module implements TOD-level dynamic masking, allowing the user to remove individual samples.
  ! The mask itself is stored in tod%scans%d%mask_dyn in the form a (2,ncut) array, where ncut is the
  ! number of flagged ranges, and mask_dyn(1,i) and mask_dyn(2,i) define the first and last excluded
  ! samples. The marked samples are flagged in comm_tod_mod->decompress_flags during decompression.
  !
  ! To use this feature:
  ! 1) Import comm_tod_dynmask_mod in comm_tod_{inst}_mod
  ! 2) Create the dynmask object in comm_tod_{inst}_smod->constructor_{inst}; choose which cuts to apply
  ! 3) Create the mask with dynmask->create in comm_tod_{inst}_smod->process_{inst}_tod
  ! 4) Report statistics to screen with dynmask->report at the end of comm_tod_{inst}_smod->process_{inst}_tod
  !
  ! Note that the cut definition is designed to be easily extended according to the needs of a given experiment
  ! by adding more cut criteria and processor routines; the behaviour of this module will evolve quickly in
  ! time.
  !
  ! Current list of cut procedures:
  ! 1) Pixel histogram  -- remove samples with extreme pixhist values; see comm_tod_pixhist_mod
  ! 2) Extreme outliers -- remove any samples with a residual more than N sigma(white noise)
  ! 3) Single samples   -- remove any samples with a residual more than N sigma(residual)
  ! 4) Excess RMS       -- remove any block of samples with a excess rms over a user-specified window length
  ! 5) Isolated samples -- remove samples that are sandwiched between two already flagged samples
  ! 6) Long chunks      -- remove any block of samples with more than x% flagged samples over a user-specified window
  ! 7) Object masks     -- remove any samples close to the Sun, Moon or Earth
  ! 8) Cosmic rays      -- remove any samples with a steep rise followed by slow decent; typical for CR impacts
  !
  type :: comm_dynmask
     class(comm_tod), pointer :: tod
     character(len=512) :: outdir
     integer(i4b)       :: output_scan, output_det    ! Write residuals after each cut to ascii for given detector and scan
     logical(lgt)       :: apply_pixhist              ! Pixel histogram outliers in units of sigma; tod%pixhist must be allocated
     real(sp)           :: threshold_extreme          ! Extreme outliers in white noise sigma; inside flagging mask; typically set to a high value
     real(sp)           :: threshold_singlesamp       ! Remove individual sample outliers above a given RMS threshold; only applied outside flagging mask
     real(sp)           :: threshold_excessRMS(N_rms) ! Excess variance in windows of N samples;
     integer(i4b)       :: window_excessRMS(N_rms)    ! Window length for excess variance
     logical(lgt)       :: remove_isolated_samples    ! Remove samples that have both masked neighbors
     real(sp)           :: threshold_longchunks       ! Maximum fraction of masked samples within long chunk
     integer(i4b)       :: window_longchunks          ! Length of long chunk to look for masked samples
     logical(lgt)       :: apply_solar_mask           ! Apply precomputed solar mask
     logical(lgt)       :: apply_moon_mask            ! Apply precomputed moon mask
     logical(lgt)       :: apply_earth_mask           ! Apply precomputed Earth mask
     real(sp)           :: threshold_cr               ! Cosmic ray threshold in sigma; defined by fast rise
     integer(i4b)       :: width_cr_mask              ! Number of samples to remove after CR hit
     integer(i8b), dimension(-1:9)  :: stats    ! Statistics for dynamic mask (ntod_tot, ncut_base, ncut_1, ncut_2,...)
   contains
     procedure :: create => create_dynamic_mask
     procedure :: remove_pixhist_outliers
     procedure :: remove_extreme_outliers
     procedure :: remove_single_outliers
     procedure :: remove_excessRMS
     procedure :: remove_isolated_samp
     procedure :: remove_longchunks
     procedure :: exclude_solar_mask
     procedure :: exclude_moon_mask
     procedure :: exclude_earth_mask
     procedure :: remove_cr_hits
     procedure :: report => report_dynamic_mask_stats
     procedure :: dump_residual
     procedure :: compress_mask
  end type comm_dynmask
    
  interface comm_dynmask
    procedure dynmask_constructor
  end interface comm_dynmask

contains
 
  !initializes a comm_dynmask class
  function dynmask_constructor(tod, cpar) result(c)
    implicit none
    class(comm_tod),     target,  intent(in) :: tod
    type(comm_params),            intent(in) :: cpar
    class(comm_dynmask), pointer            :: c
    
    integer(i4b) :: ierr
    
    ! Initialize internal parameters; disable all cuts by default
    allocate(c)
    c%outdir                  = cpar%outdir
    c%tod                     => tod
    
    c%output_scan             = -1
    c%apply_pixhist           = .false.
    c%threshold_extreme       = -1.
    c%threshold_singlesamp    = -1.
    c%threshold_excessRMS     = -1.
    c%window_excessRMS        = -1
    c%remove_isolated_samples = .false.
    c%threshold_longchunks    = -1.
    c%window_longchunks       = -1
    c%apply_solar_mask        = .false.
    c%apply_moon_mask         = .false.
    c%apply_earth_mask        = .false.
    c%threshold_cr            = -1.
    c%width_cr_mask           = -1
    c%stats                   =  0
  end function dynmask_constructor
  
  subroutine create_dynamic_mask(self, sd, det)
    implicit none
    class(comm_dynmask),  intent(inout) :: self
    class(comm_scandata), intent(inout) :: sd
    integer(i4b),         intent(in)    :: det
    
    integer(i4b) :: i, j, k, n, pix_nest, ntod, window, ntot, scan
    real(dp) :: rms0
    real(sp) :: var0, gain
    integer(i4b), allocatable, dimension(:,:) :: bad, buffer
    real(sp),     allocatable, dimension(:)   :: res, mask_dyn, var_window
    
    if (sum(sd%mask(:,det)) == 0) return 
    call timer%start(TOD_DYNMASK, self%tod%band)
    
    scan        = sd%scan
    ntod        = sd%ntod
    ntot        = count(iand(sd%flag(:,det),self%tod%flag0) .eq. 0)
    gain        = self%tod%scans(scan)%d(det)%gain
    self%stats(-1) = self%stats(-1) + ntod        ! Total number of samples
    self%stats( 0) = self%stats( 0) + ntod - ntot ! Number of samples removed by base flagging
    
    ! Compute residual
    allocate(res(ntod))
    res = (sd%tod(:,det)-real(self%tod%scans(scan)%d(det)%gain,sp)*sd%s_tot(:,det,0,1))/&
         & self%tod%scans(scan)%d(det)%N_psd%sigma0
    
    ! Generate dynamic mask
    allocate(mask_dyn(ntod))
    mask_dyn = 1.0
    
    ! Apply all cuts
    if (self%apply_pixhist)             call self%remove_pixhist_outliers(sd, det, res, mask_dyn)
    if (self%threshold_extreme > 0.)    call self%remove_extreme_outliers(sd, det, res, mask_dyn)
    if (self%threshold_cr > 0.)         call self%remove_cr_hits         (sd, det, res, mask_dyn)
    if (self%threshold_singlesamp > 0.) call self%remove_single_outliers (sd, det, res, mask_dyn)
    do i = 1, N_rms
       if (self%threshold_excessRMS(i) > 0.) then
          call self%remove_excessRMS(sd, det, self%threshold_excessRMS(i), self%window_excessRMS(i), res, mask_dyn)
       end if
    end do
    if (self%remove_isolated_samples)   call self%remove_isolated_samp   (sd, det, res, mask_dyn)
    if (self%threshold_longchunks > 0.) call self%remove_longchunks      (sd, det, res, mask_dyn)
    if (self%apply_solar_mask)          call self%exclude_solar_mask     (sd, det, res, mask_dyn)
    if (self%apply_moon_mask)           call self%exclude_moon_mask      (sd, det, res, mask_dyn)
    if (self%apply_earth_mask)          call self%exclude_earth_mask     (sd, det, res, mask_dyn)
    
    ! Compress and store final mask
    call self%compress_mask(scan, det, mask_dyn)
    
    ! Check if all samples have been excluded
    if (count(iand(sd%flag(:,det),self%tod%flag0) .eq. 0) == 0) then
       write(*,fmt='(a,a,i6,i4)') ' Dynamic mask, scan rejected due to no unflagged sample = ', trim(self%tod%freq), self%tod%scanid(scan), det
       self%tod%scans(scan)%d(det)%accept = .false.
    end if
    
    deallocate(mask_dyn, res)
    call timer%stop(TOD_DYNMASK, self%tod%band)
    
  end subroutine create_dynamic_mask
  
  
  subroutine remove_pixhist_outliers(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i, pix_nest
    
    ! Pixel histogram outliers
    q    = (self%tod%nside / self%tod%nside_pixhist)**2
    ncut = 0
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .eq. 0) then
          call ring2nest(self%tod%nside, sd%pix(i,det,1), pix_nest)
          pix_nest = pix_nest / q
          if (sd%tod(i,det) < self%tod%pixhist(4,pix_nest,det) .or. &
               & sd%tod(i,det) > self%tod%pixhist(5,pix_nest,det)) then
             mask_dyn(i)    = 0.
             sd%mask(i,det) = 0.
             sd%flag(i,det) = sd%flag(i,det) + flag_dyn
             ncut           = ncut + 1
          end if
       end if
    end do
    self%stats(1) = self%stats(1) + ncut
    
    if (self%output_scan == self%tod%scanid(sd%scan)) call self%dump_residual(sd, res, det, "pixhist")
    
  end subroutine remove_pixhist_outliers
  
  subroutine remove_extreme_outliers(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i, pix_nest
    
    ncut = 0
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .eq. 0 .and. abs(res(i)) > self%threshold_extreme) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut + 1
       end if
    end do
    self%stats(2) = self%stats(2) + ncut
    
    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "extreme")
    
  end subroutine remove_extreme_outliers
  
  subroutine remove_single_outliers(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i
    logical(lgt), allocatable, dimension(:) :: cut
    
    allocate(cut(sd%ntod))
    ncut = 0
    
    !if (self%output_scan == tod%scanid(scan)) open(58, file='flag_stage2.dat')
    
    ! Identify all samples above a given threshold and outside the mask
    do i = 1, sd%ntod
       cut(i) = (sd%mask(i,det) == 1. .and. abs(res(i)) > self%threshold_singlesamp)
       !if (self%output_scan == tod%scanid(scan) .and. sd%mask(i,det) == 1.) write(58,*) i, res(i), count(cut(i:i) == 1.)
    end do
    
    ! Check first sample manuall
    if (cut(1) .and. (.not. cut(2) .or. sd%mask(2,det) == 0.)) then
       mask_dyn(1)    = 0.
       sd%mask(1,det) = 0.
       sd%flag(1,det) = sd%flag(1,det) + flag_dyn
       ncut           = ncut + 1
    end if
    
    ! Check all intermediate samples
    do i = 2, sd%ntod-1
       if (cut(i) .and. (.not. cut(i-1) .or. sd%mask(i-1,det) == 0.) .and. (.not. cut(i+1) .or. sd%mask(i+1,det) == 0.)) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut + 1
       end if
    end do
    
    ! Check last sample manually
    if (cut(sd%ntod) .and. (.not. cut(sd%ntod-1) .or. sd%mask(sd%ntod-1,det) == 0.)) then
       mask_dyn(sd%ntod)    = 0.
       sd%mask(sd%ntod,det) = 0.
       sd%flag(sd%ntod,det) = sd%flag(sd%ntod,det) + flag_dyn
       ncut                 = ncut + 1
    end if
    !if (self%output_scan == tod%scanid(scan)) close(58)
    deallocate(cut)
    
    self%stats(3) = self%stats(3) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "singla")
    
  end subroutine remove_single_outliers

  subroutine remove_cr_hits(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    
  end subroutine remove_cr_hits

  subroutine remove_excessRMS(self, sd, det, threshold,  window, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),                           intent(in)    :: threshold
    integer(i4b),                       intent(in)    :: window
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: ncut, i, k
    real(sp)     :: var0
    real(sp), allocatable, dimension(:) :: var_window
    
    ! Look for excess variance excess in small windows; typically cosmic rays and other short glitches
    allocate(var_window(sd%ntod))
    call compute_running_variance(res, sd%mask(:,det), window, var_window, var_mean=var0, mean_full=.true.)
    var_window = sqrt(var_window)
    var0       = sqrt(var0)     
    ncut       = 0
    !if (self%output_scan == tod%scanid(scan)) open(58, file='flag_stage3.dat')
    do i = 1, sd%ntod
       !if (self%output_scan == tod%scanid(scan) .and. sd%mask(i,det) == 1.) write(58,*) i, res(i), var_window(i), var_window(i)/(threshold(4)*var0), threshold(4)*var0
       if (sd%mask(i,det) == 1. .and. var_window(i) > threshold*var0) then
          do k = max(i-window,1), min(i+window,sd%ntod)
             if (iand(sd%flag(k,det),self%tod%flag0) .eq. 0) then
                mask_dyn(k)    = 0.
                sd%mask(k,det) = 0.
                sd%flag(k,det) = sd%flag(k,det) + flag_dyn
                ncut           = ncut + 1
             end if
          end do
       end if
    end do
    !if (self%output_scan == tod%scanid(scan)) close(58)
    deallocate(var_window)
    
    self%stats(4) = self%stats(4) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "excessRMS")
    
  end subroutine remove_excessRMS
  
  subroutine remove_isolated_samp(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout)  :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i, pix_nest
    
    ncut       = 0
    
    ! Check first sample manually
    if (iand(sd%flag(1,det),self%tod%flag0) .eq. 0 .and. iand(sd%flag(2,det),self%tod%flag0) .ne. 0) then
       mask_dyn(1)    = 0.
       sd%mask(1,det) = 0.
       sd%flag(1,det) = sd%flag(1,det) + flag_dyn
       ncut           = ncut + 1
    end if
    
    ! Check intermediate samples
    do i = 2, sd%ntod-1
       if (iand(sd%flag(i-1,det),self%tod%flag0) .ne. 0 .and. iand(sd%flag(i,det),self%tod%flag0) .eq. 0 .and. iand(sd%flag(i+1,det),self%tod%flag0) .ne. 0) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut + 1
       end if
    end do
    
    ! Check last sample
    if (iand(sd%flag(sd%ntod,det),self%tod%flag0) .eq. 0 .and. iand(sd%flag(sd%ntod-1,det),self%tod%flag0) .ne. 0) then
       mask_dyn(sd%ntod)    = 0.
       sd%mask(sd%ntod,det) = 0.
       sd%flag(sd%ntod,det) = sd%flag(sd%ntod,det) + flag_dyn
       ncut              = ncut + 1
    end if
    
    self%stats(7) = self%stats(7) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "isolated")
    
  end subroutine remove_isolated_samp
  
  subroutine remove_longchunks(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: ncut, i, j, k, window
    
    ! Remove consecutive chunks with many flagged samples
    window = self%window_longchunks
    ncut       = 0
    !if (self%output_scan == self%tod%scanid(scan)) open(58, file='flag_stage7.dat')
    do i = 1, sd%ntod
       !if (self%output_scan == self%tod%scanid(scan)) write(58,*) i, res(i), iand(sd%flag(k,det),self%tod%flag0) .eq. 0
       j = max(i-window,1)
       k = min(i+window,sd%ntod)
       if (count(sd%flag(j:k,det) == flag_dyn)/real(k-j+1,sp) > self%threshold_longchunks) then
          do k = max(i-window,1), min(i+window,sd%ntod)
             if (iand(sd%flag(k,det),self%tod%flag0) .eq. 0) then
                mask_dyn(k)    = 0.
                sd%mask(k,det) = 0.
                sd%flag(k,det) = sd%flag(k,det) + flag_dyn
                ncut           = ncut + 1
             end if
          end do
       end if
    end do
    !if (self%output_scan == self%tod%scanid(scan)) close(58)
    !write(*,fmt='(a,a,i6,i4,a,f8.5,i8,i8)') ' Dynamic mask, consecutive  -- ', trim(self%tod%freq), self%tod%scanid(scan), det, ' = ', real(ncut,sp) / ntod, ncut
    
    self%stats(8) = self%stats(8) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "longchunks")
    
  end subroutine remove_longchunks
  
  subroutine exclude_solar_mask(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i, pix_nest
    
    ncut = 0
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .ne. 0) cycle
       if (self%tod%mask_solar(self%tod%scans(sd%scan)%d(det)%pix_sol(i,1),1) < 0.5) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut+1
       end if
    end do
    self%stats(9) = self%stats(9) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "solar")
    
  end subroutine exclude_solar_mask
  
  subroutine exclude_moon_mask(self, sd, det, res,  mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
    
    integer(i4b) :: q, ncut, i, pix_nest

    ncut = 0
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .ne. 0) cycle
       if (self%tod%mask_moon(self%tod%scans(sd%scan)%d(det)%pix_moon(i,1),1) < 0.5) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut+1
       end if
    end do
    
    self%stats(9) = self%stats(9) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "moon")
    
  end subroutine exclude_moon_mask
  
  subroutine exclude_earth_mask(self, sd, det, res, mask_dyn)
    implicit none
    class(comm_dynmask),                intent(inout) :: self
    class(comm_scandata),               intent(inout) :: sd
    integer(i4b),                       intent(in)    :: det
    real(sp),             dimension(:), intent(in)    :: res
    real(sp),             dimension(:), intent(inout) :: mask_dyn
     
    integer(i4b) :: q, ncut, i, pix_nest
    real(sp)     :: b_elon
    
    ncut = 0
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .ne. 0) cycle
       b_elon = max(min(int(self%tod%scans(sd%scan)%d(det)%earth_elon(i,1)/(pi/NBIN_EARTH_ELON)),NBIN_EARTH_ELON),1)
       if (self%tod%mask_earth(b_elon) < 0.5) then
          mask_dyn(i)    = 0.
          sd%mask(i,det) = 0.
          sd%flag(i,det) = sd%flag(i,det) + flag_dyn
          ncut           = ncut+1
       end if
    end do
    self%stats(9) = self%stats(9) + ncut

    if (self%output_scan == self%tod%scanid(sd%scan) .and. self%output_det == det) call self%dump_residual(sd, res, det, "earth")
    
  end subroutine exclude_earth_mask
  
  subroutine report_dynamic_mask_stats(self)
    implicit none
    class(comm_dynmask), intent(in) :: self
    
    integer(i4b) :: ierr
    integer(i8b) :: ntod
    
    ! Synchronize stats across cores
    call mpi_allreduce(MPI_IN_PLACE, self%stats, size(self%stats), MPI_INTEGER8, MPI_SUM, self%tod%comm, ierr)
    
    if (self%tod%myid == 0) then
       ntod = self%stats(-1)
       write(*,fmt='(a,a,a)')      'TOD flagging stats for ', trim(self%tod%freq), ' (      frac,          ntot     )'
       if (ntod == 0) then
          write(*,fmt='(a,f8.5,i16)') '  No non-flagged samples!'
       else
          write(*,fmt='(a,f8.5,i16)') '  Total number of samples     = ', real(self%stats(-1),dp)/ntod, ntod
          write(*,fmt='(a,f8.5,i16)') '  Base flagging               = ', real(self%stats( 0),dp)/ntod, self%stats( 0)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, pixhist       = ', real(self%stats( 1),dp)/ntod, self%stats( 1)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, extreme       = ', real(self%stats( 2),dp)/ntod, self%stats( 2)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, single spikes = ', real(self%stats( 3),dp)/ntod, self%stats( 3)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask,   5 window    = ', real(self%stats( 4),dp)/ntod, self%stats( 4)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask,  50 window    = ', real(self%stats( 5),dp)/ntod, self%stats( 5)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, 500 window    = ', real(self%stats( 6),dp)/ntod, self%stats( 6)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, single samp   = ', real(self%stats( 7),dp)/ntod, self%stats( 7)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, consecutive   = ', real(self%stats( 8),dp)/ntod, self%stats( 8)
          write(*,fmt='(a,f8.5,i16)') '  Dynamic mask, solar mask    = ', real(self%stats( 9),dp)/ntod, self%stats( 9)
          write(*,fmt='(a,f8.5,i16)') '  Final accept ratio          = ', real(ntod-sum(self%stats(0:9)),dp)/ntod, ntod-sum(self%stats(0:9))
       end if
    end if
    
  end subroutine report_dynamic_mask_stats
  
  subroutine dump_residual(self, sd, res, det, tag)
    implicit none
    class(comm_dynmask),               intent(in) :: self
    class(comm_scandata),              intent(in) :: sd
    real(sp),            dimension(:), intent(in) :: res
    integer(i4b),                      intent(in) :: det
    character(len=*),                  intent(in) :: tag

    integer(i4b) :: i
    character(len=6) :: scan_text
    character(len=4) :: det_text
    
    call int2string(sd%scan, scan_text)
    call int2string(det,  det_text)
    
    open(58, file=trim(self%outdir)//'/dynmask_'//trim(adjustl(tag))//'_'//scan_text//'_'//det_text//'.dat')
    do i = 1, sd%ntod
       if (iand(sd%flag(i,det),self%tod%flag0) .eq. 0) write(58,*) i, res(i)
    end do
    close(58)
    
  end subroutine dump_residual
  
  subroutine compress_mask(self, scan, det, mask_dyn)
    implicit none
    class(comm_dynmask),               intent(in) :: self
    integer(i4b),                      intent(in) :: scan, det
    real(sp),            dimension(:), intent(in) :: mask_dyn
    
    integer(i4b) :: i, n, nmax, ntod
    integer(i4b), allocatable, dimension(:,:) :: bad, buffer
    
    nmax = 1000 ! Maximum number of exclude ranges; extended when needed
    ntod = size(mask_dyn)
    
    ! Compress and store dynamic mask
    allocate(bad(2,nmax))
    bad = -1
    n   = 0
    do i = 1, ntod
       if (mask_dyn(i) == 0.) then
          ! Start new range if not already active
          if (bad(1,n+1) == -1) bad(1,n+1) = i
       else
          ! Close active range
          if (bad(1,n+1) /= -1 .and. bad(2,n+1) == -1) then
             bad(2,n+1) = i-1
             n          = n+1
          end if
       end if
       
       ! Increase array size if needed
       if (n == nmax) then
          nmax = 2*nmax
          allocate(buffer(2,nmax))
          buffer = -1
          buffer(:,1:nmax/2) = bad
          deallocate(bad)
          allocate(bad(2,nmax))
          bad = buffer
          deallocate(buffer)
       end if
    end do
    
    ! Close open range if needed at the end
    if (bad(1,n+1) /= -1 .and. bad(2,n+1) == -1) then
       bad(2,n+1) = ntod
       n          = n+1
    end if
    
    ! Store final array
    if (n > 0) then
       allocate(self%tod%scans(scan)%d(det)%mask_dyn(2,n))
       self%tod%scans(scan)%d(det)%mask_dyn = bad(:,1:n)
    end if
    
  end subroutine compress_mask
  
end module comm_tod_dynmask_mod
