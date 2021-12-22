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
module comm_tod_mapmaking_mod
   use comm_tod_mod
   use comm_utils
   use comm_shared_arr_mod
   use comm_map_mod
   use comm_param_mod
   implicit none

   type comm_binmap
      integer(i4b)       :: ncol, n_A, nout, nobs, npix, numprocs_shared, chunk_size
      logical(lgt)       :: shared, solve_S
      type(shared_2d_dp) :: sA_map
      type(shared_3d_dp) :: sb_map
      class(map_ptr), allocatable, dimension(:)     :: outmaps
      real(dp),       allocatable, dimension(:,:)   :: A_map
      real(dp),       allocatable, dimension(:,:,:) :: b_map
    contains
      procedure :: init    => init_binmap
      procedure :: dealloc => dealloc_binmap
      procedure :: synchronize => synchronize_binmap
   end type comm_binmap



contains





  subroutine init_binmap(self, tod, shared, solve_S)
    implicit none
    class(comm_binmap),  intent(inout) :: self
    class(comm_tod),     intent(in)    :: tod
    logical(lgt),        intent(in)    :: shared, solve_S

    integer(i4b) :: i, ierr

    self%nobs            = tod%nobs
    self%shared          = shared
    self%solve_S         = solve_S
    self%npix            = tod%info%npix
    self%numprocs_shared = tod%numprocs_shared
    self%chunk_size      = self%npix/self%numprocs_shared
    if (self%chunk_size*self%numprocs_shared < self%npix) then
       self%chunk_size   = self%chunk_size + 1
    end if
    if (solve_S) then
       self%ncol = tod%nmaps + tod%ndet - 1
       self%n_A  = tod%nmaps*(tod%nmaps+1)/2 + 4*(tod%ndet-1)
       self%nout = tod%output_n_maps + tod%n_bp_prop
       !write(*,*) 'hei!', size(tod%bp_delta,2)
    else
       self%ncol = tod%nmaps
       self%n_A  = tod%nmaps*(tod%nmaps+1)/2
       self%nout = tod%output_n_maps
    end if
!!$    write(*,*) 'nout = ', tod%output_n_maps, self%nout
!!$    call mpi_finalize(ierr)
!!$    stop
    allocate(self%outmaps(self%nout))
    do i = 1, self%nout
       self%outmaps(i)%p => comm_map(tod%info)
    end do
    allocate(self%A_map(self%n_A,self%nobs), self%b_map(self%nout,self%ncol,self%nobs))
    self%A_map = 0.d0; self%b_map = 0.d0
    if (self%shared) then
       call init_shared_2d_dp(tod%myid_shared, tod%comm_shared, &
            & tod%myid_inter, tod%comm_inter, [self%n_A,self%npix], self%sA_map)
       call mpi_win_fence(0, self%sA_map%win, ierr)
       if (self%sA_map%myid_shared == 0) self%sA_map%a = 0.d0
       call mpi_win_fence(0, self%sA_map%win, ierr)
       call init_shared_3d_dp(tod%myid_shared, tod%comm_shared, &
               & tod%myid_inter, tod%comm_inter, [self%nout,self%ncol,self%npix], self%sb_map)
       call mpi_win_fence(0, self%sb_map%win, ierr)
       if (self%sb_map%myid_shared == 0) self%sb_map%a = 0.d0
       call mpi_win_fence(0, self%sb_map%win, ierr)
    else

    end if

  end subroutine init_binmap


  subroutine dealloc_binmap(self)
    implicit none
    class(comm_binmap), intent(inout) :: self

    integer(i4b) ::  i

    if (allocated(self%A_map)) deallocate(self%A_map, self%b_map)
    if (self%sA_map%init)  call dealloc_shared_2d_dp(self%sA_map)
    if (self%sb_map%init)  call dealloc_shared_3d_dp(self%sb_map)
    if (allocated(self%outmaps)) then
       do i = 1, self%nout
          call self%outmaps(i)%p%dealloc
       end do
       deallocate(self%outmaps)
    end if

  end subroutine dealloc_binmap

  subroutine synchronize_binmap(self, tod)
    implicit none
    class(comm_binmap),  intent(inout) :: self
    class(comm_tod),     intent(in)    :: tod

    integer(i4b) :: i, j, start_chunk, end_chunk, ierr

    if (.not. self%shared) return

    do i = 0, self%numprocs_shared-1
       start_chunk = mod(self%sA_map%myid_shared+i,self%numprocs_shared)*self%chunk_size
       end_chunk   = min(start_chunk+self%chunk_size-1,self%npix-1)
       !if (i == self%numprocs_shared-1) end_chunk = self%npix-1
       do while (start_chunk < self%npix)
          if (tod%pix2ind(start_chunk) /= -1) exit
          start_chunk = start_chunk+1
       end do
       do while (end_chunk >= start_chunk)
          if (tod%pix2ind(end_chunk) /= -1) exit
          end_chunk = end_chunk-1
       end do
       if (start_chunk < self%npix)  start_chunk = tod%pix2ind(start_chunk)
       if (end_chunk >= start_chunk) end_chunk   = tod%pix2ind(end_chunk)

       call mpi_win_fence(0, self%sA_map%win, ierr)
       call mpi_win_fence(0, self%sb_map%win, ierr)
       do j = start_chunk, end_chunk
          self%sA_map%a(:,tod%ind2pix(j)+1) = self%sA_map%a(:,tod%ind2pix(j)+1) + &
                  & self%A_map(:,j)
          self%sb_map%a(:,:,tod%ind2pix(j)+1) = self%sb_map%a(:,:,tod%ind2pix(j)+1) + &
                  & self%b_map(:,:,j)
       end do
      !  stop
    end do
    
    call mpi_win_fence(0, self%sA_map%win, ierr)
    call mpi_win_fence(0, self%sb_map%win, ierr)

  end subroutine synchronize_binmap

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine bin_TOD(tod, scan, pix, psi, flag, data, binmap)
    !        call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)
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
    integer(i4b),        dimension(1:,1:),       intent(in)    :: pix, psi, flag
    real(sp),            dimension(1:,1:,1:),    intent(in)    :: data
    type(comm_binmap),                           intent(inout) :: binmap

    integer(i4b) :: det, i, t, pix_, off, nout, psi_
    real(dp)     :: inv_sigmasq

    nout = size(data,1) 
    do det = 1, size(pix,2) ! loop over all the detectors
       if (.not. tod%scans(scan)%d(det)%accept) cycle
       off         = 6 + 4*(det-1)
       inv_sigmasq = (tod%scans(scan)%d(det)%gain/tod%scans(scan)%d(det)%N_psd%sigma0)**2

       
       do t = 1, size(pix,1)
          
          if (iand(flag(t,det),tod%flag0) .ne. 0) cycle ! leave out all flagged data
          
          pix_    = tod%pix2ind(pix(t,det))  ! pixel index for pix t and detector det
          psi_    = psi(t,det)

         
          binmap%A_map(1,pix_) = binmap%A_map(1,pix_) + 1.d0                                 * inv_sigmasq
          binmap%A_map(2,pix_) = binmap%A_map(2,pix_) + tod%cos2psi(psi_)                    * inv_sigmasq
          binmap%A_map(3,pix_) = binmap%A_map(3,pix_) + tod%cos2psi(psi_)**2                 * inv_sigmasq
          binmap%A_map(4,pix_) = binmap%A_map(4,pix_) + tod%sin2psi(psi_)                    * inv_sigmasq
          binmap%A_map(5,pix_) = binmap%A_map(5,pix_) + tod%cos2psi(psi_)*tod%sin2psi(psi_) * inv_sigmasq
          binmap%A_map(6,pix_) = binmap%A_map(6,pix_) + tod%sin2psi(psi_)**2                 * inv_sigmasq
          !binmap%A_map(1,pix_) = binmap%A_map(8,pix_) + 1.d0

          do i = 1, nout
             binmap%b_map(i,1,pix_) = binmap%b_map(i,1,pix_) + data(i,t,det)                      * inv_sigmasq
             binmap%b_map(i,2,pix_) = binmap%b_map(i,2,pix_) + data(i,t,det) * tod%cos2psi(psi_) * inv_sigmasq
             binmap%b_map(i,3,pix_) = binmap%b_map(i,3,pix_) + data(i,t,det) * tod%sin2psi(psi_) * inv_sigmasq
          end do
          
          if (binmap%solve_S .and. det < tod%ndet) then
             binmap%A_map(off+1,pix_) = binmap%A_map(off+1,pix_) + 1.d0               * inv_sigmasq 
             binmap%A_map(off+2,pix_) = binmap%A_map(off+2,pix_) + tod%cos2psi(psi_) * inv_sigmasq
             binmap%A_map(off+3,pix_) = binmap%A_map(off+3,pix_) + tod%sin2psi(psi_) * inv_sigmasq
             binmap%A_map(off+4,pix_) = binmap%A_map(off+4,pix_) + 1.d0               * inv_sigmasq
             do i = 1, nout
                binmap%b_map(i,det+3,pix_) = binmap%b_map(i,det+3,pix_) + data(i,t,det) * inv_sigmasq 
             end do
          end if

          
          
         end do
         
         
         
      end do
      



      ! stop

  end subroutine bin_TOD


   ! differential TOD computation, written with WMAP in mind.
   subroutine bin_differential_TOD(tod, data, pix, psi, flag, x_imarr, pmask, b, M_diag, scan, comp_S, b_mono)
    ! Routine to bin differential time ordered data
    ! Assumes white noise after correctiom from correlated noise and calibrated data
    ! 
    ! Arguments:
    ! ----------
    ! tod:    
    !         
    ! scan:   integer
    !         scan number
    ! pix:    2-dimensional array
    !         Number of pixels from scandata
    ! psi:    2-dimensional array
    !         Pointing angle pr pixel
    ! flag:   2-dimensional array
    !         Flagged data to be excluded from the mapmaking
    ! data:   2-dim array
    !         Array of calibrated data
    ! comp_S: logical
    !         Whether or not to explicitly solve for spurious component
    ! b_mono: 2-dim array
    !         not implemented
    ! pmask:  1-dim array
    !         Healpix map of good and bad values
    ! x_imarr: 1-dim array
    !          imbalance parameters, duplicated so it's (x_1, x_1, x_2, x_2)
    ! Returns:
    ! ----------
    ! binmap: pointer
    !         Pointer to array of binned map?
    ! 
      implicit none
      class(comm_tod), intent(in)                               :: tod
      integer(i4b), intent(in)                                  :: scan
      real(sp), dimension(1:, 1:, 1:), intent(in)               :: data
      integer(i4b), dimension(1:), intent(in)                   :: flag
      real(sp), dimension(0:), intent(in)                       :: pmask
      integer(i4b), dimension(1:, 1:), intent(in)               :: pix, psi
      real(dp), dimension(1:), intent(in)                       :: x_imarr
      real(dp), dimension(0:, 1:, 1:), intent(inout)            :: b
      real(dp), dimension(0:, 1:), intent(inout)                :: M_diag
      real(dp), dimension(0:, 1:, 1:), intent(inout), optional  :: b_mono
      logical(lgt), intent(in)                                  :: comp_S

      integer(i4b) :: det, i, t, nout
      real(dp)     :: inv_sigmasq, d, p, var, dx_im, x_im

      integer(i4b) :: lpix, rpix, lpsi, rpsi

      integer(i4b) :: f_A, f_B

      nout = size(b, dim=3)
      ! Note that x_imarr is duplicated
      dx_im = 0.5*(x_imarr(1) - x_imarr(3))
      x_im = 0.5*(x_imarr(1) + x_imarr(3))

      if (tod%scans(scan)%d(1)%accept) then
         !inv_sigmasq = 0.d0 
         var = 0
         do det = 1, 4
           var = var + (tod%scans(scan)%d(det)%N_psd%sigma0/tod%scans(scan)%d(det)%gain)**2/4
         end do
         inv_sigmasq = 1/var
         do t = 1, tod%scans(scan)%ntod
            if (iand(flag(t),tod%flag0) .ne. 0) cycle

            lpix = pix(t, 1)
            rpix = pix(t, 2)
            lpsi = psi(t, 1)
            rpsi = psi(t, 2)
            f_A = pmask(rpix)
            f_B = pmask(lpix)

            do i = 1, nout
               d = 0.d0
               p = 0.d0
               ! d = (d13 + d14 + d23 + d24)/4
               ! p = (d13 + d14 - d23 - d24)/4
               do det = 1, 4
                 d = d + data(i, t, det)/4
                 if (det < 3) then
                    p = p + data(i, t, det)/4
                 else
                    p = p - data(i, t, det)/4
                 end if
               end do
               ! T
               b(lpix, 1, i) = b(lpix, 1, i) + f_A*((1.d0+x_im)*d + dx_im*p)*inv_sigmasq
               b(rpix, 1, i) = b(rpix, 1, i) - f_B*((1.d0-x_im)*d - dx_im*p)*inv_sigmasq
               ! Q
               b(lpix, 2, i) = b(lpix, 2, i) + f_A*((1.d0+x_im)*p + dx_im*d)*tod%cos2psi(lpsi)*inv_sigmasq
               b(rpix, 2, i) = b(rpix, 2, i) - f_B*((1.d0-x_im)*p - dx_im*d)*tod%cos2psi(rpsi)*inv_sigmasq
               ! U
               b(lpix, 3, i) = b(lpix, 3, i) + f_A*((1.d0+x_im)*p + dx_im*d)*tod%sin2psi(lpsi)*inv_sigmasq
               b(rpix, 3, i) = b(rpix, 3, i) - f_B*((1.d0-x_im)*p - dx_im*d)*tod%sin2psi(rpsi)*inv_sigmasq
               ! S
               if (comp_S) then
                 b(lpix, 4, i) = b(lpix, 4, i) + f_A*((1.d0+x_im)*p + dx_im*d)*inv_sigmasq
                 b(rpix, 4, i) = b(rpix, 4, i) - f_B*((1.d0-x_im)*p - dx_im*d)*inv_sigmasq
               end if
            end do

            M_diag(lpix, 1) = M_diag(lpix, 1) + f_A*inv_sigmasq
            M_diag(rpix, 1) = M_diag(rpix, 1) + f_B*inv_sigmasq
            M_diag(lpix, 2) = M_diag(lpix, 2) + f_A*inv_sigmasq*tod%cos2psi(lpsi)**2
            M_diag(rpix, 2) = M_diag(rpix, 2) + f_B*inv_sigmasq*tod%cos2psi(rpsi)**2
            M_diag(lpix, 3) = M_diag(lpix, 3) + f_A*inv_sigmasq*tod%sin2psi(lpsi)**2
            M_diag(rpix, 3) = M_diag(rpix, 3) + f_B*inv_sigmasq*tod%sin2psi(rpsi)**2
            if (comp_S) then
              M_diag(lpix, 4) = M_diag(lpix, 4) + f_A*inv_sigmasq
              M_diag(rpix, 4) = M_diag(rpix, 4) + f_B*inv_sigmasq
            else
              ! Not a true diagonal term, just the off-diagonal estimate of the
              ! covariance for each pixel.
              M_diag(lpix, 4) = M_diag(lpix, 4)+f_A*inv_sigmasq*tod%sin2psi(lpsi)*tod%cos2psi(lpsi)
              M_diag(rpix, 4) = M_diag(rpix, 4)+f_B*inv_sigmasq*tod%sin2psi(rpsi)*tod%cos2psi(rpsi)
            end if

         end do
       end if

end subroutine bin_differential_TOD

   subroutine compute_Ax(tod, x_imarr, pmask, comp_S, M_diag, x_in, y_out)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Code to compute matrix product P^T N^-1 P m
      ! y = Ax
      ! Explicitly removes the monopole in temperature, since this mode is
      ! formally solvable for due to transmission imbalance, but takes the
      ! majority of the BICG-Stab iterations to solve for.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      class(comm_tod),                 intent(in)              :: tod
      real(dp),     dimension(1:),     intent(in)              :: x_imarr
      real(sp),     dimension(0:),     intent(in)              :: pmask
      logical(lgt), intent(in)                                 :: comp_S
      real(dp),                dimension(:,:), intent(in) :: M_diag
      real(dp),     dimension(0:, 1:), intent(in),    optional :: x_in
      real(dp),     dimension(0:, 1:), intent(inout), optional :: y_out

      integer(i4b), allocatable, dimension(:)         :: flag
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi

      logical(lgt) :: finished
      integer(i4b) :: j, k, ntod, ndet, lpix, rpix, lpsi, rpsi, ierr
      integer(i4b) :: nhorn, t, f_A, f_B, nside, npix, nmaps
      real(dp)     :: inv_sigmasq, var, iA, iB, sA, sB, d, p, x_im, dx_im, monopole
      real(dp), allocatable, dimension(:,:) :: x, y
      nhorn = tod%nhorn
      ndet  = tod%ndet
      nside = tod%nside
      nmaps = tod%nmaps
      npix  = 12*nside**2

      if (comp_S) then
         allocate(x(0:npix-1,nmaps+1), y(0:npix-1,nmaps+1))
      else
         allocate(x(0:npix-1,nmaps),   y(0:npix-1,nmaps))
      end if
      if (tod%myid == 0) then
         finished = .false.
         call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
         x = x_in
      end if
      call mpi_bcast(x, size(x),  MPI_DOUBLE_PRECISION, 0, tod%info%comm, ierr)

      x_im   = 0.5*(x_imarr(1) + x_imarr(3))
      dx_im  = 0.5*(x_imarr(1) - x_imarr(3))
      y      = 0.d0
      do j = 1, tod%nscan
         ntod = tod%scans(j)%ntod
         allocate (pix(ntod, nhorn))             ! Decompressed pointing
         allocate (psi(ntod, nhorn))             ! Decompressed pol angle
         allocate (flag(ntod))                   ! Decompressed flags
         !do k = 1, tod%ndet
         if (tod%scans(j)%d(1)%accept) then
            call tod%decompress_pointing_and_flags(j, 1, pix, &
                & psi, flag)

            !inv_sigmasq = 0.d0 
            var = 0
            do k = 1, 4
               var = var + (tod%scans(j)%d(k)%N_psd%sigma0/tod%scans(j)%d(k)%gain)**2/4
               !inv_sigmasq = inv_sigmasq  + (tod%scans(j)%d(k)%gain/tod%scans(j)%d(k)%sigma0)**2
            end do
            inv_sigmasq = 1.d0/var

            do t = 1, ntod

               if (iand(flag(t),tod%flag0) .ne. 0) cycle
               lpix = pix(t, 1)
               rpix = pix(t, 2)
               lpsi = psi(t, 1)
               rpsi = psi(t, 2)

               f_A = pmask(rpix)
               f_B = pmask(lpix)
               ! This is the model for each timestream
               iA = x(lpix, 1)
               iB = x(rpix, 1)
               if (comp_S) then
                  sA = x(lpix, 2)*tod%cos2psi(lpsi) + x(lpix, 3)*tod%sin2psi(lpsi) + x(lpix, 4)
                  sB = x(rpix, 2)*tod%cos2psi(rpsi) + x(rpix, 3)*tod%sin2psi(rpsi) + x(rpix, 4)
               else
                  sA = x(lpix, 2)*tod%cos2psi(lpsi) + x(lpix, 3)*tod%sin2psi(lpsi)
                  sB = x(rpix, 2)*tod%cos2psi(rpsi) + x(rpix, 3)*tod%sin2psi(rpsi)
               end if
               d  = (1.d0+x_im)*iA - (1.d0-x_im)*iB + dx_im*(sA + sB)
               p  = (1.d0+x_im)*sA - (1.d0-x_im)*sB + dx_im*(iA + iB)
               ! Temperature
               y(lpix, 1) = y(lpix, 1) + f_A*((1.d0 + x_im)*d + dx_im*p) * inv_sigmasq
               y(rpix, 1) = y(rpix, 1) - f_B*((1.d0 - x_im)*d - dx_im*p) * inv_sigmasq
               ! Q
               y(lpix, 2) = y(lpix, 2) + f_A*((1.d0 + x_im)*p + dx_im*d) * tod%cos2psi(lpsi)*inv_sigmasq
               y(rpix, 2) = y(rpix, 2) - f_B*((1.d0 - x_im)*p - dx_im*d) * tod%cos2psi(rpsi)*inv_sigmasq
               ! U
               y(lpix, 3) = y(lpix, 3) + f_A*((1.d0 + x_im)*p + dx_im*d) * tod%sin2psi(lpsi)*inv_sigmasq
               y(rpix, 3) = y(rpix, 3) - f_B*((1.d0 - x_im)*p - dx_im*d) * tod%sin2psi(rpsi)*inv_sigmasq
               ! S
               if (comp_S) then
                 y(lpix, 4) = y(lpix, 4) + f_A*((1.d0 + x_im)*p + dx_im*d) * inv_sigmasq
                 y(rpix, 4) = y(rpix, 4) - f_B*((1.d0 - x_im)*p - dx_im*d) * inv_sigmasq
               end if
            end do
         end if
         deallocate (pix, psi, flag)
      end do

      if (tod%myid == 0) then
         call mpi_reduce(y, y_out, size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
         monopole = sum(y_out(:,1)*M_diag(:,1)*pmask) &
                & / sum(M_diag(:,1)*pmask)
         y_out(:,1) = y_out(:,1) - monopole
      else
         call mpi_reduce(y, y,     size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
      end if

      deallocate(x, y)

   end subroutine compute_Ax

  subroutine finalize_binned_map(tod, binmap, handle, rms, scale, chisq_S, Sfilename, mask)
    !
    ! Routine to finalize the binned maps
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! binmap:
    ! handle:  planck_rng derived type
    !          Healpix definition for random number generation
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(in)    :: tod
    type(comm_binmap),                    intent(inout) :: binmap
    type(planck_rng),                     intent(inout) :: handle
    class(comm_map),                      intent(inout) :: rms
    real(dp),                             intent(in)    :: scale
    real(dp),        dimension(1:,1:),    intent(out),   optional :: chisq_S
    character(len=*),                     intent(in),    optional :: Sfilename
    real(sp),        dimension(0:),       intent(in),    optional :: mask

    integer(i4b) :: i, j, k, nmaps, ierr, ndet, ncol, n_A, off, ndelta
    integer(i4b) :: det, nout, np0, comm, myid, nprocs
    real(dp), allocatable, dimension(:,:)   :: A_inv, As_inv
    real(dp), allocatable, dimension(:,:,:) :: b_tot, bs_tot
    real(dp), allocatable, dimension(:)     :: W, eta
    real(dp), allocatable, dimension(:,:)   :: A_tot
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: smap 

    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    nout  = size(binmap%sb_map%a,dim=1)
!!$    write(*,*) 'nout = ', nout
!!$    call mpi_finalize(ierr)
!!$    stop
    nmaps = tod%info%nmaps
    ndet  = tod%ndet
    n_A   = size(binmap%sA_map%a,dim=1)
    ncol  = size(binmap%sb_map%a,dim=2)
    ndelta = 0; if (present(chisq_S)) ndelta = size(chisq_S,dim=2)

    if (present(Sfilename)) then
       info => comm_mapinfo(tod%comm, tod%info%nside, 0, ndet-1, .false.)
       smap => comm_map(info)
       smap%map = 0.d0
    end if

    ! Collect contributions from all nodes
    call mpi_win_fence(0, binmap%sA_map%win, ierr)
    if (binmap%sA_map%myid_shared == 0) then
       do i = 1, size(binmap%sA_map%a, 1)
          call mpi_allreduce(MPI_IN_PLACE, binmap%sA_map%a(i, :), size(binmap%sA_map%a, 2), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, binmap%sA_map%comm_inter, ierr)
       end do
    end if
      call mpi_win_fence(0, binmap%sA_map%win, ierr)
      call mpi_win_fence(0, binmap%sb_map%win, ierr)
      if (binmap%sb_map%myid_shared == 0) then
         do i = 1, size(binmap%sb_map%a, 1)
            call mpi_allreduce(mpi_in_place, binmap%sb_map%a(i, :, :), size(binmap%sb_map%a(1, :, :)), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, binmap%sb_map%comm_inter, ierr)
         end do
      end if
      call mpi_win_fence(0, binmap%sb_map%win, ierr)

      allocate (A_tot(n_A, 0:np0 - 1), b_tot(nout, nmaps, 0:np0 - 1), bs_tot(nout, ncol, 0:np0 - 1), W(nmaps), eta(nmaps))

      A_tot = binmap%sA_map%a(:, tod%info%pix + 1)
      b_tot = binmap%sb_map%a(:, 1:nmaps, tod%info%pix + 1)
      bs_tot = binmap%sb_map%a(:, :, tod%info%pix + 1)


      ! Solve for local map and rms
      allocate (A_inv(nmaps, nmaps), As_inv(ncol, ncol))
      if (present(chisq_S)) chisq_S = 0.d0

      do i = 0, np0 - 1
         if (all(b_tot(1, :, i) == 0.d0)) then
            if (.not. present(chisq_S)) then
               rms%map(i, :) = 0.d0
               do k = 1, nout
                  binmap%outmaps(k)%p%map(i, :) = 0.d0
               end do
            end if
            cycle
         end if

         A_inv = 0.d0
         A_inv(1, 1) = A_tot(1, i)
         A_inv(2, 1) = A_tot(2, i)
         A_inv(1, 2) = A_inv(2, 1)
         A_inv(2, 2) = A_tot(3, i)
         A_inv(3, 1) = A_tot(4, i)
         A_inv(1, 3) = A_inv(3, 1)
         A_inv(3, 2) = A_tot(5, i)
         A_inv(2, 3) = A_inv(3, 2)
         A_inv(3, 3) = A_tot(6, i)
         if (present(chisq_S)) then
            As_inv = 0.d0
            As_inv(1:nmaps, 1:nmaps) = A_inv
            do det = 1, ndet - 1
               off = 6 + 4*(det - 1)
               As_inv(1, 3 + det) = A_tot(off + 1, i)
               As_inv(3 + det, 1) = As_inv(1, 3 + det)
               As_inv(2, 3 + det) = A_tot(off + 2, i)
               As_inv(3 + det, 2) = As_inv(2, 3 + det)
               As_inv(3, 3 + det) = A_tot(off + 3, i)
               As_inv(3 + det, 3) = As_inv(3, 3 + det)
               As_inv(3 + det, 3 + det) = A_tot(off + 4, i)
            end do
         end if

         call invert_singular_matrix(A_inv, 1d-12)
         do k = 1, tod%output_n_maps
            b_tot(k, 1:nmaps, i) = matmul(A_inv, b_tot(k, 1:nmaps, i))
         end do
         if (present(chisq_S)) then
            call invert_singular_matrix(As_inv, 1d-12)
            bs_tot(1, 1:ncol, i) = matmul(As_inv, bs_tot(1, 1:ncol, i))
            do k = tod%output_n_maps + 1, nout
               bs_tot(k, 1:ncol, i) = matmul(As_inv, bs_tot(k, 1:ncol, i))
            end do
         end if
    

         if (present(chisq_S)) then
            do j = 1, ndet - 1
               !write(*,*) mask(tod%info%pix(i + 1)), As_inv(nmaps + j, nmaps +j)
               if (mask(tod%info%pix(i + 1)) == 0.) cycle
               if (As_inv(nmaps + j, nmaps + j) <= 0.d0) cycle
               chisq_S(j, 1) = chisq_S(j, 1) + bs_tot(1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               do k = 2, ndelta
                  chisq_S(j, k) = chisq_S(j, k) + bs_tot(tod%output_n_maps + k - 1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               end do
               if (present(Sfilename) .and. As_inv(nmaps+j,nmaps+j) > 0.d0) then
                  !smap%map(i,j) = bs_tot(tod%output_n_maps+1,nmaps+j,i) / sqrt(As_inv(nmaps+j,nmaps+j))
                  smap%map(i,j) = bs_tot(1,nmaps+j,i) / sqrt(As_inv(nmaps+j,nmaps+j))
               end if
            end do
         end if
         do j = 1, nmaps
            rms%map(i, j) = sqrt(A_inv(j, j))*scale
            do k = 1, tod%output_n_maps
               binmap%outmaps(k)%p%map(i, j) = b_tot(k, j, i)*scale
            end do
         end do
      end do

      if (present(chisq_S)) then
         if (myid == 0) then
            call mpi_reduce(mpi_in_place, chisq_S, size(chisq_S), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
         else
            call mpi_reduce(chisq_S, chisq_S, size(chisq_S), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
         end if
         if (present(Sfilename)) then
            call smap%writeFITS(Sfilename)
            call smap%dealloc
         end if
      end if

      


      deallocate (A_inv, As_inv, A_tot, b_tot, bs_tot, W, eta)

   end subroutine finalize_binned_map

   subroutine run_bicgstab(tod, handle, bicg_sol, npix, nmaps, num_cg_iters, epsil, procmask, map_full, M_diag, b_map, l, prefix, postfix, comp_S)
     !
     !
     !  Subroutine that runs the biconjugate gradient-stabilized mapmaking
     !  routine on differential data, solving the P_m^T N^-1 P x = P_m^T N^-1 d
     !  mapmaking equation, where P_m takes into account the asymmetric masking
     !  
     !  Explicitly removes the monople from P_m^T N^-1 d, so that it is not
     !  for in the routine, which wastes many CG iterations.
     !
     !  Arguments (fixed):
     !  ------------------
     !  tod: comm_tod
     !
     !  npix: int
     !
     !  nmaps: int
     !
     !  epsil: real (dp)
     !
     !  procmask: real(sp)
     !
     !  M_diag: real (dp)
     !
     !  b_map: real (dp)
     !
     !  l: int
     ! 
     ! comp_S: logical
     !         Whether or not to explicitly solve for spurious component
     !
     !  Arguments (modified):
     !  ---------------------
     !  handle: planck_rng
     !
     !  bicg_sol: real (dp)
     !
     !  num_cg_iters: int
     !
     !  map_full: real (dp)
     !
     implicit none
     class(comm_tod),                         intent(in) :: tod
     type(planck_rng),                     intent(inout) :: handle
     real(dp),         dimension(:, :, :), intent(inout) :: bicg_sol
     integer(i4b),                            intent(in) :: npix, nmaps
     integer(i4b),                         intent(inout) :: num_cg_iters
     real(dp),                                intent(in) :: epsil
     real(sp),                  dimension(:), intent(in) :: procmask
     real(dp),               dimension(:), intent(inout) :: map_full
     real(dp),                dimension(:,:), intent(in) :: M_diag
     real(dp),              dimension(:,:,:), intent(in) :: b_map
     integer(i4b),                            intent(in) :: l
     character(len=512),                      intent(in) :: prefix
     character(len=512),                      intent(in) :: postfix
     logical(lgt), intent(in)                            :: comp_S




     integer(i4b)                               :: recomp_freq = 10
     real(dp),     allocatable, dimension(:, :) :: m_buf
     integer(i4b)                               :: i_max, i_min, ierr, i
     real(dp)                                   :: delta_0
     real(dp)                                   :: alpha, beta, sigma_mono
     real(dp),     allocatable, dimension(:, :) :: r, s, q
     real(dp)                                   :: monopole
     logical(lgt)                               :: finished, write_cg
     real(dp)                                   :: rho_old, rho_new
     real(dp)                                   :: omega, delta_r, delta_s
     real(dp),     allocatable, dimension(:, :) :: rhat, r0, shat, p, phat, v
     real(dp),        allocatable, dimension(:) :: determ
     character(len=512)                         :: i_str, l_str

     write_cg = .false.
     !write_cg = tod%first_call

     if (tod%myid==0) then
        if (comp_S) then
           allocate (r     (0:npix-1, nmaps+1))
           allocate (rhat  (0:npix-1, nmaps+1))
           allocate (r0    (0:npix-1, nmaps+1))
           allocate (q     (0:npix-1, nmaps+1))
           allocate (p     (0:npix-1, nmaps+1))
           allocate (s     (0:npix-1, nmaps+1))
           allocate (shat  (0:npix-1, nmaps+1))
           allocate (m_buf (0:npix-1, nmaps+1))
           allocate (phat  (0:npix-1, nmaps+1))
           allocate (v     (0:npix-1, nmaps+1))
        else
           allocate (r     (0:npix-1, nmaps))
           allocate (rhat  (0:npix-1, nmaps))
           allocate (r0    (0:npix-1, nmaps))
           allocate (q     (0:npix-1, nmaps))
           allocate (p     (0:npix-1, nmaps))
           allocate (s     (0:npix-1, nmaps))
           allocate (shat  (0:npix-1, nmaps))
           allocate (m_buf (0:npix-1, nmaps))
           allocate (phat  (0:npix-1, nmaps))
           allocate (v     (0:npix-1, nmaps))
           allocate (determ(0:npix-1))
        end if

        i_max = 500
        i_min = 0

        if (.false. .and. l == 1) then
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, bicg_sol(:,:,1), v)
           r = b_map(:, :, l) - v 
        else
           r = b_map(:, :, l)
        end if
        monopole = sum(b_map(:,1,l)*M_diag(:,1)*procmask) &
               & / sum(M_diag(:,1)*procmask)
        r0 = b_map(:, :, l)
        r0(:, 1) = b_map(:, 1, l) - monopole
        ! This is if we are not solving for S
        ! In this case, M_diag(:,4) is the QU covariance term
        if (comp_S) then
          rhat =  r/M_diag
        else
          determ = M_diag(:,2)*M_diag(:,3) - M_diag(:,4)**2
          rhat(:,1) =  r(:,1)/M_diag(:,1)
          rhat(:,2) = (r(:,2)*M_diag(:,3)- r(:,2)*M_diag(:,4))/determ
          rhat(:,3) = (r(:,3)*M_diag(:,2)- r(:,3)*M_diag(:,4))/determ
        end if
        
        delta_r = sum(r*rhat)
        delta_0 = delta_r
        delta_s = delta_s

        omega = 1d0
        alpha = 1d0

        rho_new = sum(r0*r)
        i = 0
        bicg: do
           i = i + 1
           rho_old = rho_new
           call update_status(status, 'dot product')
           rho_new = sum(r0*r)
           call update_status(status, 'done dot product')
           if (rho_new == 0d0) then
             if (tod%verbosity > 1) write(*,*) '| Residual norm is zero'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if

           if (i==1) then
              p = r
           else
              beta = (rho_new/rho_old) * (alpha/omega)
              p = r + beta*(p - omega*v)
           end if

           if (comp_S) then
             phat =  p/M_diag
           else
             phat(:,1) =  p(:,1)/M_diag(:,1)
             phat(:,2) = (p(:,2)*M_diag(:,3)- p(:,2)*M_diag(:,4))/determ
             phat(:,3) = (p(:,3)*M_diag(:,2)- p(:,3)*M_diag(:,4))/determ
           end if
           
           call update_status(status, 'v=A phat')
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, phat, v)
           call update_status(status, 'done')
           num_cg_iters = num_cg_iters + 1

           alpha         = rho_new/sum(r0*v)
           s             = r - alpha*v
           if (comp_S) then
             shat =  s/M_diag
           else
             shat(:,1) =  s(:,1)/M_diag(:,1)
             shat(:,2) = (s(:,2)*M_diag(:,3)- s(:,2)*M_diag(:,4))/determ
             shat(:,3) = (s(:,3)*M_diag(:,2)- s(:,3)*M_diag(:,4))/determ
           end if
           delta_s       = sum(s*shat)

           if (tod%verbosity > 1) then 
              write(*,101) 2*i-1, delta_s/delta_0
101           format (' |', 6X, I4, ':   delta_s/delta_0:',  2X, ES9.2)
           end if

           bicg_sol(:,:,l) = bicg_sol(:,:,l) + alpha*phat

           if (write_cg) then
             write(i_str, '(I0.3)') 2*i-1
             write(l_str, '(I1)') l
             call write_map(trim(prefix)//'cgest_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & bicg_sol(:,:,l))
             call write_map(trim(prefix)//'cgres_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & r)
           end if

           if (delta_s .le. (delta_0*epsil) .and. 2*i-1 .ge. i_min) then
              if (tod%verbosity > 1) write(*,*) '| Reached bicg-stab tolerance'
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if

           call update_status(status, 'q=A shat')
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, shat, q)
           call update_status(status, 'done')

           omega         = sum(q*s)/sum(q*q)
           bicg_sol(:,:,l) = bicg_sol(:,:,l) + omega*shat


           if (omega == 0d0) then
             if (tod%verbosity > 1) write(*,*) '| omega is zero'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if

           if (mod(i, recomp_freq) == 1 .or. beta > 1.d8) then
              call update_status(status, 'A xhat')
              call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, bicg_sol(:,:,l), r)
              call update_status(status, 'done')
              r(:,1) = b_map(:,1,l) - monopole - r(:,1)
              r(:,2) = b_map(:,2,l)  - r(:,2)
              r(:,3) = b_map(:,3,l)  - r(:,3)
              if (comp_S)  r(:,4) = b_map(:,4,l)  - r(:,4)
           else
              r = s - omega*q
           end if

           if (write_cg) then
             write(i_str, '(I0.3)') 2*i
             write(l_str, '(I1)') l
             call write_map(trim(prefix)//'cgest_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & bicg_sol(:,:,l))
             call write_map(trim(prefix)//'cgres_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & r)
           end if

           ! If you are not solving for S
           ! If you are solving for S
           if (comp_S) then
             rhat = r/M_diag
           else
             rhat(:,1) =  r(:,1)/M_diag(:,1)
             rhat(:,2) = (r(:,2)*M_diag(:,3)- r(:,2)*M_diag(:,4))/determ
             rhat(:,3) = (r(:,3)*M_diag(:,2)- r(:,3)*M_diag(:,4))/determ
           end if
           delta_r      = sum(r*rhat)
           num_cg_iters = num_cg_iters + 1

           if (tod%verbosity > 1) then 
              write(*,102) 2*i, delta_r/delta_0
102           format (' |', 6X, I4, ':   delta_r/delta_0:',  2X, ES9.2)
           end if
           if (delta_r .le. delta_0*epsil .and. 2*i .ge. i_min) then
              if (tod%verbosity > 1) write(*,*) '| Reached bicg-stab tolerance'
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if
           if (i==i_max) then
             if (tod%verbosity > 1) write(*,*) '| Reached maximum number of iterations'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if
        end do bicg

        if (l == 1) then
           ! Maximum likelihood monopole
           monopole = sum((bicg_sol(:,1,1)-map_full)*M_diag(:,1)*procmask) &
                  & / sum(M_diag(:,1)*procmask)
           if (trim(tod%operation) == 'sample') then
              ! Add fluctuation term if requested
              sigma_mono = sum(M_diag(:,1) * procmask)
              if (sigma_mono > 0.d0) sigma_mono = 1.d0 / sqrt(sigma_mono)
              if (tod%verbosity > 1) then
                write(*,*) '| monopole, fluctuation sigma'
                write(*,*) '| ', monopole, sigma_mono
              end if
              monopole = monopole + sigma_mono * rand_gauss(handle)
           end if
           bicg_sol(:,1,1) = bicg_sol(:,1,1) - monopole
        end if
     else
        loop: do while (.true.) 
           call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
           if (finished) exit loop
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag)
        end do loop
     end if
     if (tod%myid == 0) deallocate (r, rhat, s, r0, q, shat, p, phat, v, m_buf)
     if (tod%myid == 0 .and. .not. comp_S) deallocate (determ)

   end subroutine run_bicgstab


end module comm_tod_mapmaking_mod
