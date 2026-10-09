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
module comm_tod_linemap_mod
   use comm_tod_driver_mod
   use comm_shared_arr_mod
   use comm_map_mod
   implicit none

   type comm_linemap
      integer(i4b)       :: nline, n_A, ndet, nobs, npix, numprocs_shared, chunk_size
      type(shared_3d_dp) :: sA_map
      type(shared_3d_dp) :: sb_map
      !real(dp), allocatable, dimension(:,:)   :: line_ratio
      real(dp), allocatable, dimension(:,:,:) :: A_map
      real(dp), allocatable, dimension(:,:,:) :: b_map
      real(dp), allocatable, dimension(:,:,:) :: A_det
      real(dp), allocatable, dimension(:,:)   :: b_det
    contains
      procedure :: synchronize => synchronize_linemap
   end type comm_linemap

   interface comm_linemap
      procedure constructor_linemap
   end interface comm_linemap
   
contains

  function constructor_linemap(tod) result(c)
    implicit none
    class(comm_tod),        intent(in)    :: tod
    class(comm_linemap),    pointer       :: c

    integer(i4b) :: i, ierr
    class(comm_mapinfo), pointer:: info => null()

    allocate(c)
    
    call timer%start(TOD_ALLOC, tod%band)
    c%nobs            = tod%pixcache%nobs
    c%npix            = tod%info%npix
    c%numprocs_shared = tod%numprocs_shared
    c%chunk_size      = c%npix/c%numprocs_shared

    c%nline = size(tod%line_ratio,1)
    c%n_A   = c%nline*(c%nline+1)/2
    c%ndet  = size(tod%line_ratio,2)
    info    => comm_mapinfo(tod%info%comm, tod%info%nside, 0, c%nline, .false.)

    allocate(c%A_map(c%ndet,c%n_A,c%nobs), c%b_map(c%ndet,c%nline,c%nobs))
    c%A_map = 0.d0; c%b_map = 0.d0
    call init_shared_3d_dp(tod%myid_shared, tod%comm_shared, &
         & tod%myid_inter, tod%comm_inter, [c%ndet,c%n_A,c%npix], c%sA_map)
    call mpi_win_fence(0, c%sA_map%win, ierr)
    if (c%sA_map%myid_shared == 0) c%sA_map%a = 0.d0
    call mpi_win_fence(0, c%sA_map%win, ierr)
    call init_shared_3d_dp(tod%myid_shared, tod%comm_shared, &
         & tod%myid_inter, tod%comm_inter, [c%ndet,c%nline,c%npix], c%sb_map)
    call mpi_win_fence(0, c%sb_map%win, ierr)
    if (c%sb_map%myid_shared == 0) c%sb_map%a = 0.d0
    call mpi_win_fence(0, c%sb_map%win, ierr)
    call timer%stop(TOD_ALLOC, tod%band)

    allocate(c%A_det(c%nline,c%nline,c%ndet))
    allocate(c%b_det(c%nline,c%ndet))
    
  end function constructor_linemap

  subroutine deallocate_linemap(self)
    implicit none
    class(comm_linemap), pointer, intent(inout) :: self

    integer(i4b) ::  i

    if (allocated(self%A_map))       deallocate(self%A_map, self%b_map)
    if (self%sA_map%init)            call dealloc_shared_3d_dp(self%sA_map)
    if (self%sb_map%init)            call dealloc_shared_3d_dp(self%sb_map)
    if (allocated(self%A_det))       deallocate(self%A_det, self%b_det)
    deallocate(self)
    nullify(self)

  end subroutine deallocate_linemap

  subroutine synchronize_linemap(self, tod)
    implicit none
    class(comm_linemap),  intent(inout) :: self
    class(comm_tod),      intent(in)    :: tod

    integer(i4b) :: i, j, start_chunk, end_chunk, ind1, ind2, ierr

    call timer%start(TOD_MAPSYN, tod%band)
    
    do i = 0, self%numprocs_shared-1
       start_chunk = mod(self%sA_map%myid_shared+i,self%numprocs_shared)*self%chunk_size
       end_chunk   = min(start_chunk+self%chunk_size-1,self%npix-1)
       call tod%pixcache%get_ind_range(start_chunk, end_chunk, ind1, ind2)
       !if (self%sA_map%myid_shared == 0) write(*,*) tod%pixcache%ind2pix(1), tod%pixcache%ind2pix(tod%pixcache%nobs), start_chunk, ind1, end_chunk, ind2

       call mpi_win_fence(0, self%sA_map%win, ierr)
       call mpi_win_fence(0, self%sb_map%win, ierr)
       if (ind1 > 0 .and. ind2 > 0) then
          do j = ind1, ind2
             self%sA_map%a(:,:,tod%pixcache%ind2pix(j)+1) = self%sA_map%a(:,:,tod%pixcache%ind2pix(j)+1) + &
                  & self%A_map(:,:,j)
             self%sb_map%a(:,:,tod%pixcache%ind2pix(j)+1) = self%sb_map%a(:,:,tod%pixcache%ind2pix(j)+1) + &
                  & self%b_map(:,:,j)
          end do
       end if
    end do
    call mpi_win_fence(0, self%sA_map%win, ierr)
    call mpi_win_fence(0, self%sb_map%win, ierr)

    ! Collect contributions from all nodes
    call mpi_win_fence(0, self%sA_map%win, ierr)
    if (self%sA_map%myid_shared == 0) then
       do i = 1, size(self%sA_map%a, 1)
          call mpi_allreduce(MPI_IN_PLACE, self%sA_map%a(i,:,:), size(self%sA_map%a(1,:,:)), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, self%sA_map%comm_inter, ierr)
       end do
    end if
    call mpi_win_fence(0, self%sA_map%win, ierr)
    call mpi_win_fence(0, self%sb_map%win, ierr)
    if (self%sb_map%myid_shared == 0) then
       do i = 1, size(self%sb_map%a, 1)
          call mpi_allreduce(mpi_in_place, self%sb_map%a(i, :, :), size(self%sb_map%a(1, :, :)), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, self%sb_map%comm_inter, ierr)
       end do
    end if
    call mpi_win_fence(0, self%sb_map%win, ierr)

    call timer%stop(TOD_MAPSYN, tod%band)
    
  end subroutine synchronize_linemap


    ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine bin_linemap(tod, scan, pix, flag, res, linemap)
    !        call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, linemap)
    ! Routine to bin time ordered data
    ! Assumes white noise after correctiom from correlated noise and calibrated data
    ! 
    ! Arguments:
    ! ----------
    ! tod:    
    !         
    ! scan:   integer
    !         scan number
    ! pix:    2-dimentional array
    !         Number of pixels from scandata
    ! psi:    2-dimentional array
    !         Pointing angle pr pixel
    ! flag:   2-dimentional array
    !         Flagged data to be excluded from the mapmaking
    ! data:   2-dim array
    !         Array of calibrated data
    !
    ! Returns:
    ! ----------
    ! binmap: pointer
    !         Pointer to array of binned map?
    ! 

    implicit none
    class(comm_tod),                             intent(in)    :: tod
    integer(i4b),                                intent(in)    :: scan
    integer(i4b),        dimension(1:,1:),       intent(in)    :: pix, flag
    real(sp),            dimension(1:,1:),       intent(in)    :: res
    type(comm_linemap),                           intent(inout) :: linemap

    integer(i4b) :: det, i, j, k, t, pix_, off, nline, psi_
    real(dp)     :: inv_sigmasq, eff

    nline = linemap%nline
 
    call timer%start(TOD_MAPBIN, tod%band)
    do det = 1, size(pix,2) ! loop over all the detectors
       if (.not. tod%scans(scan)%d(det)%accept) cycle
       inv_sigmasq = (tod%scans(scan)%d(det)%gain/tod%scans(scan)%d(det)%N_psd%sigma0)**2

       ! polarization efficiency
       do t = 1, size(pix,1)
          
          if (iand(flag(t,det),tod%flag0) .ne. 0) cycle ! leave out all flagged data
          
          pix_    = tod%pixcache%pix2ind(pix(t,det))  ! pixel index for pix t and detector det
          do i = 1, nline
             linemap%b_map(det,i,pix_) = linemap%b_map(det,i,pix_) + res(t,det) * inv_sigmasq
             do j = 1, i
                k = i*(i-1)/2 + j ! Position in A matrix
                linemap%A_map(det,k,pix_) = linemap%A_map(det,k,pix_) + inv_sigmasq
             end do
          end do
       end do
    end do
    call timer%stop(TOD_MAPBIN, tod%band)
    
  end subroutine bin_linemap

   subroutine finalize_linemap(tod, linemap)
    !
    ! Routine to finalize the binned maps
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! linemap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(in)    :: tod
    type(comm_linemap),                    intent(inout) :: linemap

    integer(i4b) :: i, j, k, l, nmaps, ierr, ndet, nline, n_A
    integer(i4b) :: det, np0, comm, myid, nprocs
    real(dp), allocatable, dimension(:,:)   :: A_inv
    real(dp), allocatable, dimension(:,:,:) :: b_tot
    real(dp), allocatable, dimension(:)     :: b_sum
    real(dp), allocatable, dimension(:,:,:) :: A_tot
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: smap 

    call timer%start(TOD_MAPSOLVE, tod%band)
    
    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    ndet  = linemap%ndet
    n_A   = size(linemap%sA_map%a,dim=2)
    nline = linemap%nline
    nmaps = nline
    
    allocate (A_tot(ndet, n_A, 0:np0-1), b_tot(ndet, nmaps, 0:np0-1))
    A_tot = linemap%sA_map%a(:,:, tod%info%pix + 1)
    b_tot = linemap%sb_map%a(:,:, tod%info%pix + 1)

    ! Solve for local map and rms
    allocate (A_inv(nline, nline), b_sum(nline))
    do i = 0, np0 - 1
       if (all(b_tot(1, :, i) == 0.d0)) then
          tod%rms_line%map(i,:) = 0.d0
          tod%map_line%map(i,:) = 0.d0
          cycle
       end if

       ! Build full system by summing over detectors, multiplying with line ratios
       A_inv = 0.d0
       b_sum = 0.d0
       do det = 1, ndet
          do j = 1, nline
             b_sum(j) = b_sum(j) + b_tot(det,j,i) * tod%line_ratio(j,det)
             do k = 1, j
                l = j*(j-1)/2 + k ! Position in A matrix
                A_inv(j,k) = A_inv(j,k) + A_tot(det,l,i) * &
                     & tod%line_ratio(j,det) * tod%line_ratio(k,det)
             end do
          end do
       end do

       ! Symmetrize A
       do j = 1, nline
          do k = 1, j
             A_inv(k,j) = A_inv(j,k)
          end do
       end do
       
       ! Solve system
       call invert_singular_matrix(A_inv, 1d-12)
       tod%map_line%map(i,:) = matmul(A_inv, b_sum)

       ! Diagonal matrix; store RMS
       do j = 1, nline
          tod%rms_line%map(i,j) = sqrt(A_inv(j, j))
       end do
    end do

    ! Subtract monopole and dipole
    do j = 1, nline
       call tod%map_line%subtract_mono_dipole(mask=tod%mask_line, col=j)
    end do

    ! Suppress noise
    !call tod%map_line%wiener_filter(tod%rms_line, spin0=.true.)
    
    deallocate (A_inv, A_tot, b_tot, b_sum)
    call timer%stop(TOD_MAPSOLVE, tod%band)
      
  end subroutine finalize_linemap

  subroutine sample_line_ratios(tod, linemap, update_line_ratio)
    !
    ! Routine to samole line ratios given existing template
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! linemap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(inout)    :: tod
    type(comm_linemap),                   intent(inout) :: linemap
    logical(lgt),                         intent(in)    :: update_line_ratio 

    integer(i4b) :: i, j, k, l, nmaps, ierr, ndet, nline, n_A
    integer(i4b) :: det, np0, comm, myid, nprocs
    logical(lgt) :: skip
    real(dp), allocatable, dimension(:,:)   :: A_inv, A, A_p
    real(dp), allocatable, dimension(:,:,:) :: b_tot
    real(dp), allocatable, dimension(:)     :: b_sum, b_p
    real(dp), allocatable, dimension(:,:,:) :: A_tot
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: smap 

    call timer%start(TOD_MAPSOLVE, tod%band)
    
    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    ndet  = linemap%ndet
    n_A   = size(linemap%sA_map%a,dim=2)
    nline = linemap%nline
    nmaps = nline
    
    allocate (A_tot(ndet, n_A, 0:np0-1), b_tot(ndet, nmaps, 0:np0-1))
    A_tot = linemap%sA_map%a(:,:, tod%info%pix + 1)
    b_tot = linemap%sb_map%a(:,:, tod%info%pix + 1)

    ! Build linear system for line ratios per detector
    allocate (A_inv(nline, nline), b_sum(nline), A_p(nline,nline), b_p(nline), A(nline,nline))
    do det = 1, ndet
       A_inv = 0.d0
       b_sum = 0.d0       
       loop_pix: do i = 0, np0 - 1

          if (all(b_tot(det,:,i) == 0.d0)) cycle
!!$          do j = 1, nline
!!$             if (tod%line_ref_map(j)%p%map(i,1) == 0.d0) cycle loop_pix
!!$          end do
          
          ! Build full system by summing over detectors, multiplying with line ratios
          do j = 1, nline
             if (tod%line_ref_map(j)%p%map(i,1) == 0.d0) cycle
             b_sum(j) = b_sum(j) + b_tot(det,j,i) * tod%line_ref_map(j)%p%map(i,1)
             do k = 1, j
                if (tod%line_ref_map(k)%p%map(i,1) == 0.d0) cycle
                l = j*(j-1)/2 + k ! Position in A matrix
                A_inv(j,k) = A_inv(j,k) + A_tot(det,l,i) * &
                     & tod%line_ref_map(j)%p%map(i,1) * tod%line_ref_map(k)%p%map(i,1)
             end do
          end do
       end do loop_pix

       ! Symmetrize A
       do j = 1, nline
          do k = 1, j
             A_inv(k,j) = A_inv(j,k)
          end do
       end do

       call mpi_allreduce(MPI_IN_PLACE, A_inv, size(A_inv), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
       call mpi_allreduce(MPI_IN_PLACE, b_sum, size(b_sum), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)

       ! Store system for later use
       linemap%A_det(:,:,det) = A_inv
       linemap%b_det(:,det)   = b_sum
       if (.not. update_line_ratio) cycle
       
       ! Solve data system 
       A = A_inv
       call invert_singular_matrix(A, 1d-12)
       b_sum = matmul(A, b_sum)
      
       ! Weighted sum with prior
!!$       A_p = 0.d0
!!$       do j = 1, nline
!!$          A_p(j,j) = 1.d0/tod%line_ratio_prior(2,j,det)**2
!!$          b_p(j)   =      tod%line_ratio_prior(1,j,det)
!!$          if (tod%myid == 0) write(*,*) 'd = ', A_inv(:,j), b_sum(j)
!!$          if (tod%myid == 0) write(*,*) 'p = ', A_p(:,j),   b_p(j)
!!$       end do
!!$       b_sum = matmul(A_inv,b_sum) + matmul(A_p,b_p)
!!$       A_inv = A_inv + A_p
!!$       call invert_singular_matrix(A_inv, 1d-12)
!!$       b_sum = matmul(A_inv, b_sum)
       
       tod%line_ratio(:,det) = b_sum
       if (tod%myid == 0) write(*,*) ' Line ratios = ', tod%line_ratio(:,det), det
       
       ! Add fluctuation term
       
    end do

    deallocate (A_inv, A_tot, b_tot, b_sum, A_p, b_p)
    call timer%stop(TOD_MAPSOLVE, tod%band)
      
  end subroutine sample_line_ratios

  subroutine sample_line_scaling(tod, linemap)
    !
    ! Routine to samole line ratios given existing template
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! linemap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(inout)    :: tod
    type(comm_linemap),                    intent(inout) :: linemap

    integer(i4b) :: i, j, k, l, nmaps, ierr, ndet, nline, n_A
    integer(i4b) :: det, np0, comm, myid, nprocs
    real(dp)     :: A, b, s, alpha
    real(dp), allocatable, dimension(:,:,:) :: b_tot
    real(dp), allocatable, dimension(:,:,:) :: A_tot

    call timer%start(TOD_MAPSOLVE, tod%band)
    
    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    ndet  = linemap%ndet
    n_A   = size(linemap%sA_map%a,dim=2)
    nline = linemap%nline
    nmaps = nline
    
    allocate (A_tot(ndet, n_A, 0:np0-1), b_tot(ndet, nmaps, 0:np0-1))
    A_tot = linemap%sA_map%a(:,:, tod%info%pix + 1)
    b_tot = linemap%sb_map%a(:,:, tod%info%pix + 1)

    ! Build linear system for line ratios per detector
    A = 0.d0
    b = 0.d0
    do det = 1, ndet
       do i = 0, np0 - 1
          if (all(b_tot(det,:,i) == 0.d0)) cycle
          
          ! Build full signal
          s = 0.d0
          do j = 1, nline
             s = s + tod%line_ratio(j,det) * tod%map_line%map(i,j)
          end do

          ! Build linear system
          b = b + b_tot(det,1,i) * s
          A = A + A_tot(det,1,i) * s**2
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE, A, 1, &
         & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
    call mpi_allreduce(MPI_IN_PLACE, b, 1, &
         & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)

    ! Solve data system
    alpha = b/A

    ! Add fluctuation term


    ! Rescale line ratios
    tod%line_ratio = alpha * tod%line_ratio
    
    if (tod%myid == 0) write(*,*) ' Line scaling = ', alpha
    do det = 1, ndet
       if (tod%myid == 0) write(*,fmt='(a,2f16.4,i5)') ' U [uK/(Kkms)] = ', tod%line_ratio(:,det), det
    end do
    
    deallocate (A_tot, b_tot)
    call timer%stop(TOD_MAPSOLVE, tod%band)
      
  end subroutine sample_line_scaling

  
  subroutine sample_line_ratios_amp(tod, linemap)
    !
    ! Routine to samole line ratios given existing template
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! linemap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(inout)    :: tod
    type(comm_linemap),                    intent(inout) :: linemap

    integer(i4b) :: i, j, k, l, ierr, ndet, nline,  det, np0
    real(dp)     :: chisq
    real(dp), allocatable, dimension(:)     :: p

    call timer%start(TOD_MAPSOLVE, tod%band)
    
    np0   = tod%info%np
    ndet  = linemap%ndet
    nline = linemap%nline

    allocate(p(nline*ndet-1))
    p = 1.d0
    call powell(p, powell_chisq_lineamp, ierr, tolerance=1d-4)

    !Make sure final solution is stored
    chisq = powell_chisq_lineamp(p)

    do det = 1, ndet
       if (tod%myid == 0) write(*,*) ' Line ratios = ', tod%line_ratio(:,det), det
    end do

    deallocate(p)

    ! Fit overall amplitudes
    
  contains

    function powell_chisq_lineamp(x)
      implicit none
      real(dp), dimension(:), intent(in),  optional :: x
      real(dp)                                      :: powell_chisq_lineamp

      integer(i4b) :: i, j, l, ierr
      real(dp)     :: alpha, A, b, sigma, chisq 

      ! Initialize line ratios, first detector has amplitude equal to 1
      if (tod%myid == 0) then
         tod%line_ratio(:,1) = 1.d0
         j = 1
         do i = 1, ndet
            do l = 1, nline
               if (i == 1 .and. l == 1) then
                  tod%line_ratio(l,i) = 1.d0
               else
                  tod%line_ratio(l,i) = x(j)
                  j                   = j+1
               end if
            end do
         end do
      end if
      call mpi_bcast(tod%line_ratio, nline*ndet, MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)

      if (any(tod%line_ratio < 0.d0)) then
         powell_chisq_lineamp = 1.d30
         return
      end if
      
      ! Solve for (unnormalized) amplitude maps
      call finalize_linemap(tod, linemap)

      ! Fit for scaling
      do l = 1, nline
         A = 0.d0; b = 0.d0
         do i = 0, np0 - 1
            if (tod%line_ref_map(l)%p%map(i,1) == 0.d0) cycle
            sigma = tod%rms_line%map(i,l)
            b = b + tod%map_line%map(i,l) * tod%line_ref_map(l)%p%map(i,1)    / sigma**2
            A = A +                         tod%line_ref_map(l)%p%map(i,1)**2 / sigma**2
          end do
          call mpi_allreduce(MPI_IN_PLACE, A, 1, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
          call mpi_allreduce(MPI_IN_PLACE, b, 1, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
          alpha = b/A
          tod%map_line%map(:,l) = tod%map_line%map(:,l)/alpha
          if (tod%myid == 0) write(*,*) ' l, alpha = ', l, alpha
      end do

      ! Compute chisq
      chisq = 0.d0
      do l = 1, nline
         do i = 0, np0 - 1
            if (tod%line_ref_map(l)%p%map(i,1) == 0.d0) cycle
            sigma = 1.d0 !tod%rms_line%map(i,l)
            chisq = chisq + (tod%map_line%map(i,l)-tod%line_ref_map(l)%p%map(i,1))**2 / &
                 & sigma**2
         end do
      end do
      call mpi_allreduce(MPI_IN_PLACE, chisq, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)

      powell_chisq_lineamp = chisq

      if (tod%myid == 0) write(*,*) ' chisq, p = ', chisq, real(x,sp)

      if (chisq < 100) then
         call tod%rms_line%writeFITS("rms.fits")
         call tod%map_line%writeFITS("map.fits")

         ! Compute chisq
         chisq = 0.d0
         do l = 1, nline
            do i = 0, np0 - 1
               if (tod%line_ref_map(l)%p%map(i,1) == 0.d0) cycle
               sigma = tod%rms_line%map(i,l)
               chisq = chisq + (tod%map_line%map(i,l)-tod%line_ref_map(l)%p%map(i,1))**2 / &
                    & sigma**2
               write(*,*) l, i, chisq, tod%map_line%map(i,l), tod%line_ref_map(l)%p%map(i,1), sigma
            end do
         end do

         call mpi_finalize(ierr)
         stop
      end if
      
    end function powell_chisq_lineamp
    
  end subroutine sample_line_ratios_amp

  subroutine sample_line_ratios_amp2(tod, linemap)
    !
    ! Routine to samole line ratios given existing template
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! linemap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(inout)    :: tod
    type(comm_linemap),                    intent(inout) :: linemap

    integer(i4b) :: i, j, k, l, ierr, ndet, nline,  det, np0
    real(dp)     :: chisq
    real(dp), allocatable, dimension(:)     :: p

    call timer%start(TOD_MAPSOLVE, tod%band)
    
    np0   = tod%info%np
    ndet  = linemap%ndet
    nline = linemap%nline

    allocate(p(nline*ndet))
    j = 1
    do i = 1, ndet
       do l = 1, nline
          p(j) = tod%line_ratio(l,i) 
          j    = j+1
       end do
    end do

    !p = sum(tod%line_ratio)/size(tod%line_ratio)
    if (tod%myid == 0) write(*,*) ' p = ', p
    call powell(p, powell_chisq_lineamp, ierr, tolerance=1d-4)

    !Make sure final solution is stored
    chisq = powell_chisq_lineamp(p)

    do det = 1, ndet
       if (tod%myid == 0) write(*,*) ' Line ratios = ', tod%line_ratio(:,det), det
    end do

    deallocate(p)

    ! Fit overall amplitudes
    
  contains

    function powell_chisq_lineamp(x)
      implicit none
      real(dp), dimension(:), intent(in),  optional :: x
      real(dp)                                      :: powell_chisq_lineamp

      integer(i4b) :: i, j, l, ierr
      real(dp)     :: alpha, A, b, sigma, chisq0, chisq1, lambda0, lambda1, res(nline)

      ! Initialize line ratios, first detector has amplitude equal to 1
      if (tod%myid == 0) then
         j = 1
         do i = 1, ndet
            do l = 1, nline
               tod%line_ratio(l,i) = x(j)
               j                   = j+1
            end do
         end do
      end if
      call mpi_bcast(tod%line_ratio, nline*ndet, MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)

      if (any(tod%line_ratio < 0.d0)) then
         powell_chisq_lineamp = 1.d30
         return
      end if
      
      ! Solve for (unnormalized) amplitude maps
      !write(*,*) ' line = ', tod%line_ratio
      call finalize_linemap(tod, linemap)
      !write(*,*) ' NaN = ', count(tod%map_line%map /= tod%map_line%map)

      ! Fit for scaling
      do l = 1, nline
         A = 0.d0; b = 0.d0
         do i = 0, np0 - 1
            sigma = 1.d0 !tod%rms_line%map(i,l)
            if (tod%line_ref_map(l)%p%map(i,1) == 0.d0 .or. sigma == 0.d0) cycle
            !if (sigma == 0.d0) write(*,*) 'sigma', l, i, tod%rms_line%map(i,l)
            b = b + tod%map_line%map(i,l) * tod%line_ref_map(l)%p%map(i,1)    / sigma**2
            A = A +                         tod%line_ref_map(l)%p%map(i,1)**2 / sigma**2
          end do
          call mpi_allreduce(MPI_IN_PLACE, A, 1, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
          call mpi_allreduce(MPI_IN_PLACE, b, 1, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
          !write(*,*) ' l,A,b = ', l, A, b
          alpha = b/A
          tod%map_line%map(:,l) = tod%map_line%map(:,l)/alpha
          if (tod%myid == 0) write(*,*) ' l, alpha = ', l, alpha
      end do

      ! Compute chisq
      chisq0 = 0.d0
      do l = 1, nline
         do i = 0, np0 - 1
            if (tod%line_ref_map(l)%p%map(i,1) == 0.d0) cycle
            sigma = 1.d0 !tod%rms_line%map(i,l)
            chisq0 = chisq0 + (tod%map_line%map(i,l)-tod%line_ref_map(l)%p%map(i,1))**2 / &
                 & sigma**2
         end do
      end do
      call mpi_allreduce(MPI_IN_PLACE, chisq0, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, linemap%sA_map%comm_shared, ierr)
      
      ! Add residual-based contribution
      chisq1 = 0.d0 
      do i = 1, ndet
         res = matmul(linemap%A_det(:,:,i),tod%line_ratio(:,i)) - linemap%b_det(:,i)
         chisq1 = chisq1 + sum(res**2)
      end do

      lambda0 = 0.d-3  ! Reference map difference
      lambda1 = 1.d0   ! TOD residual
      chisq = lambda0*chisq0 + lambda1*chisq1
      powell_chisq_lineamp = chisq

      if (tod%myid == 0) write(*,*) ' chisq, p = ', chisq0,chisq1, real(x,sp)

      if (.false. .and. chisq < 100) then
         call tod%rms_line%writeFITS("rms.fits")
         call tod%map_line%writeFITS("map.fits")

         ! Compute chisq
         chisq = 0.d0
         do l = 1, nline
            do i = 0, np0 - 1
               if (tod%line_ref_map(l)%p%map(i,1) == 0.d0) cycle
               sigma = tod%rms_line%map(i,l)
               chisq = chisq + (tod%map_line%map(i,l)-tod%line_ref_map(l)%p%map(i,1))**2 / &
                    & sigma**2
               write(*,*) l, i, chisq, tod%map_line%map(i,l), tod%line_ref_map(l)%p%map(i,1), sigma
            end do
         end do

         call mpi_finalize(ierr)
         stop
      end if
      
    end function powell_chisq_lineamp
    
  end subroutine sample_line_ratios_amp2

  
end module comm_tod_linemap_mod
