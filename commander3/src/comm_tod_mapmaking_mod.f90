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

    call timer%start(TOD_ALLOC, tod%band)
    self%nobs            = tod%nobs
    self%shared          = shared
    self%solve_S         = solve_S
    self%npix            = tod%info%npix
    self%numprocs_shared = tod%numprocs_shared
    self%chunk_size      = self%npix/self%numprocs_shared
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
    !write(*,*) 'nout = ', tod%output_n_maps, self%nout
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
    call timer%stop(TOD_ALLOC, tod%band)

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
          if(tod%nmaps > 1) then
            binmap%A_map(2,pix_) = binmap%A_map(2,pix_) + tod%cos2psi(psi_)                    * inv_sigmasq
            binmap%A_map(3,pix_) = binmap%A_map(3,pix_) + tod%cos2psi(psi_)**2                 * inv_sigmasq
            binmap%A_map(4,pix_) = binmap%A_map(4,pix_) + tod%sin2psi(psi_)                    * inv_sigmasq
            binmap%A_map(5,pix_) = binmap%A_map(5,pix_) + tod%cos2psi(psi_)*tod%sin2psi(psi_) * inv_sigmasq
            binmap%A_map(6,pix_) = binmap%A_map(6,pix_) + tod%sin2psi(psi_)**2                 * inv_sigmasq
            !binmap%A_map(1,pix_) = binmap%A_map(8,pix_) + 1.d0
          end if 

          do i = 1, nout
             binmap%b_map(i,1,pix_) = binmap%b_map(i,1,pix_) + data(i,t,det)                      * inv_sigmasq
             if(tod%nmaps > 1) then
               binmap%b_map(i,2,pix_) = binmap%b_map(i,2,pix_) + data(i,t,det) * tod%cos2psi(psi_) * inv_sigmasq
               binmap%b_map(i,3,pix_) = binmap%b_map(i,3,pix_) + data(i,t,det) * tod%sin2psi(psi_) * inv_sigmasq
             end if
          end do
          
          if (binmap%solve_S .and. det < tod%ndet) then
             binmap%A_map(off+1,pix_) = binmap%A_map(off+1,pix_) + 1.d0               * inv_sigmasq 
             if(tod%nmaps > 1) then
               binmap%A_map(off+2,pix_) = binmap%A_map(off+2,pix_) + tod%cos2psi(psi_) * inv_sigmasq
               binmap%A_map(off+3,pix_) = binmap%A_map(off+3,pix_) + tod%sin2psi(psi_) * inv_sigmasq
               binmap%A_map(off+4,pix_) = binmap%A_map(off+4,pix_) + 1.d0               * inv_sigmasq
             end if 
            do i = 1, nout
                binmap%b_map(i,det+3,pix_) = binmap%b_map(i,det+3,pix_) + data(i,t,det) * inv_sigmasq 
             end do
          end if
          
       end do
    end do

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


      call timer%start(TOD_MAPBIN, tod%band)

      nout = size(b, dim=3)
      ! Note that x_imarr is duplicated
      dx_im = 0.5*(x_imarr(1) - x_imarr(3))
      x_im = 0.5*(x_imarr(1) + x_imarr(3))

      if (tod%scans(scan)%d(1)%accept) then
         !inv_sigmasq = 0.d0 
         var = 0
         ! 16 because each variable is divided by 4, variance goes as Var(aX) = a^2 Var(X)
         do det = 1, 4
           var = var + (tod%scans(scan)%d(det)%N_psd%sigma0/tod%scans(scan)%d(det)%gain)**2/16
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
      call timer%stop(TOD_MAPBIN, tod%band)

end subroutine bin_differential_TOD

   subroutine compute_Ax(tod, x_imarr, pmask, comp_S, M_diag, x, y, x_in, y_out)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Code to compute matrix product P^T N^-1 P m
      ! y = Ax
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      class(comm_tod),                 intent(in)              :: tod
      real(dp),     dimension(1:),     intent(in)              :: x_imarr
      real(sp),     dimension(0:),     intent(in)              :: pmask
      logical(lgt), intent(in)                                 :: comp_S
      real(dp),                dimension(:,:), intent(in) :: M_diag
      real(dp),                dimension(1:,0:), intent(inout), optional :: x
      real(dp),                dimension(1:,0:), intent(inout), optional :: y
      real(dp),     dimension(0:, 1:), intent(in),    optional :: x_in
      real(dp),     dimension(0:, 1:), intent(inout), optional :: y_out

      integer(i4b), allocatable, dimension(:)         :: flag
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi

      logical(lgt) :: finished
      integer(i4b) :: j, k, ntod, ndet, lpix, rpix, lpsi, rpsi, ierr
      integer(i4b) :: nhorn, t, f_A, f_B, nside, npix, nmaps
      real(dp)     :: inv_sigmasq, var, iA, iB, sA, sB, d, p, x_im, dx_im, x_im_pos, x_im_neg, sigT, sigP, lcos2psi, lsin2psi, rcos2psi, rsin2psi, monopole
      nhorn = tod%nhorn
      ndet  = tod%ndet
      nside = tod%nside
      nmaps = tod%nmaps
      npix  = 12*nside**2

      x      = 0.d0
      if (tod%myid == 0) then
         finished = .false.
         call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
         x = transpose(x_in)
      end if
      call mpi_bcast(x, size(x),  MPI_DOUBLE_PRECISION, 0, tod%info%comm, ierr)

      x_im     = 0.5*(x_imarr(1) + x_imarr(3))
      dx_im    = 0.5*(x_imarr(1) - x_imarr(3))
      x_im_pos = 1.d0 + x_im
      x_im_neg = 1.d0 - x_im

      y      = 0.d0

      do j = 1, tod%nscan
         ntod = tod%scans(j)%ntod
         allocate (pix(ntod, nhorn))             ! Decompressed pointing
         allocate (psi(ntod, nhorn))             ! Decompressed pol angle
         allocate (flag(ntod))                   ! Decompressed flags
         !do k = 1, tod%ndet
         if (tod%scans(j)%d(1)%accept) then
            !call update_status(status, 'decomp')
            call tod%decompress_pointing_and_flags(j, 1, pix, &
                & psi, flag)
            !call update_status(status, 'done')

            var = 0.d0
            do k = 1, 4
               var = var + (tod%scans(j)%d(k)%N_psd%sigma0/tod%scans(j)%d(k)%gain)**2/16
            end do
            inv_sigmasq = 1.d0/var

            !call update_status(status, 'loop')
            do t = 1, ntod

               if (iand(flag(t),tod%flag0) .ne. 0) cycle
               lpix = pix(t, 1)
               rpix = pix(t, 2)
               lcos2psi = tod%cos2psi(psi(t,1))
               lsin2psi = tod%sin2psi(psi(t,1))
               rcos2psi = tod%cos2psi(psi(t,2))
               rsin2psi = tod%sin2psi(psi(t,2))

               ! This is the model for each timestream
               iA = x(1,lpix)
               sA = x(2,lpix)*lcos2psi + x(3,lpix)*lsin2psi
               if (comp_S) sA = sA + x(4,lpix)
               iB = x(1,rpix)
               sB = x(2,rpix)*rcos2psi + x(3,rpix)*rsin2psi
               if (comp_S) sB = sB + x(4,rpix)

               d  = (x_im_pos*iA - x_im_neg*iB + dx_im*(sA + sB)) * inv_sigmasq
               p  = (x_im_pos*sA - x_im_neg*sB + dx_im*(iA + iB)) * inv_sigmasq

               if (pmask(rpix) > 0.5d0) then
                  sigT      = x_im_pos*d + dx_im*p 
                  sigP      = x_im_pos*p + dx_im*d
                  y(1,lpix) = y(1,lpix) + sigT 
                  y(2,lpix) = y(2,lpix) + sigP * lcos2psi
                  y(3,lpix) = y(3,lpix) + sigP * lsin2psi
                  if (comp_S) y(4,lpix) = y(4,lpix) + sigP
               end if

               if (pmask(lpix) > 0.5d0) then
                  sigT       = -(x_im_neg*d - dx_im*p)
                  sigP       = -(x_im_neg*p - dx_im*d)
                  y(1,rpix) = y(1,rpix) + sigT
                  y(2,rpix) = y(2,rpix) + sigP * rcos2psi
                  y(3,rpix) = y(3,rpix) + sigP * rsin2psi
                  if (comp_S) y(4,rpix) = y(4,rpix) + sigP
               end if

            end do
            !call update_status(status, 'done')
         end if
         deallocate (pix, psi, flag)
      end do

      if (tod%myid == 0) then
         call mpi_reduce(y, x, size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
         y_out    = transpose(x)
         monopole = sum(x(1,:)*M_diag(:,1)*pmask) &
                & / sum(M_diag(:,1)*pmask)
         y_out(:,1) = y_out(:,1) - monopole
      else
         call mpi_reduce(y, x,     size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
      end if

   end subroutine compute_Ax

  subroutine finalize_binned_map_unpol(tod, binmap, rms, scale, chisq_S, mask)
    !
    ! Routine to finalize temperature-only binned maps
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! binmap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(in)    :: tod
    type(comm_binmap),                    intent(inout) :: binmap
    class(comm_map),                      intent(inout) :: rms
    real(dp),                             intent(in)    :: scale
    real(dp),        dimension(1:,1:),    intent(out),   optional :: chisq_S
    real(sp),        dimension(0:),       intent(in),    optional :: mask


    integer(i4b) :: i, j, k, nmaps, ierr, ndet, ncol, n_A, off, ndelta
    integer(i4b) :: det, nout, np0, comm, myid, nprocs
    real(dp)     :: A_inv, As_inv
    real(dp), allocatable, dimension(:,:,:) :: b_tot, bs_tot
    real(dp), allocatable, dimension(:)     :: W, eta
    real(dp), allocatable, dimension(:,:)   :: A_tot

    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    nout  = size(binmap%sb_map%a,dim=1)
    ndet  = tod%ndet
    n_A   = size(binmap%sA_map%a,dim=1)
    ncol  = size(binmap%sb_map%a,dim=2)
    ndelta = 0; if (present(chisq_S)) ndelta = size(chisq_S,dim=2)

    ! Collect contributions from all nodes
    !TODO: figure out why this causes a crash
!    call mpi_win_fence(0, binmap%sA_map%win, ierr)
!    if (binmap%sA_map%myid_shared == 0) then
!       do i = 1, size(binmap%sA_map%a, 1)
!          write(*,*) "at point A, i=", i, binmap%sA_map%comm_inter
!          call mpi_allreduce(MPI_IN_PLACE, binmap%sA_map%a(i, :), size(binmap%sA_map%a, 2), size(binmap%sA_map%a, 2), &
!               & MPI_DOUBLE_PRECISION, MPI_SUM, binmap%sA_map%comm_inter, ierr)
!       end do
!    end if
!      call mpi_win_fence(0, binmap%sA_map%win, ierr)
!      call mpi_win_fence(0, binmap%sb_map%win, ierr)
!      if (binmap%sb_map%myid_shared == 0) then
!         do i = 1, size(binmap%sb_map%a, 1)
!            write(*,*) "at point B, i=", i, binmap%sb_map%comm_inter
!            call mpi_allreduce(mpi_in_place, binmap%sb_map%a(i, :, :), size(binmap%sb_map%a(1, :, :)), &
!                 & MPI_DOUBLE_PRECISION, MPI_SUM, binmap%sb_map%comm_inter,ierr)
!         end do
!      end if
!      call mpi_win_fence(0, binmap%sb_map%win, ierr)



      allocate (A_tot(n_A, 0:np0 - 1), b_tot(nout, 1, 0:np0 - 1), bs_tot(nout, ncol, 0:np0 - 1), W(1), eta(1))
      A_tot = binmap%sA_map%a(:, tod%info%pix + 1)
      b_tot = binmap%sb_map%a(:, 1:1, tod%info%pix + 1)
      bs_tot = binmap%sb_map%a(:, :, tod%info%pix + 1)

      ! Solve for local map and rms
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
         
         ! compute average
         A_inv = 1.d0/A_tot(1, i) 


         if (present(chisq_S)) then
            As_inv = A_inv
            write(*,*) "chisq_S not supported. TODO: whatever this is supposed to be" 
         end if


         if (present(chisq_S)) then
            ! TODO: compute inverse of chisq_S?
            do j = 1, ndet - 1
               if (mask(tod%info%pix(i + 1)) == 0.) cycle
               if (As_inv <= 0.d0) cycle
               chisq_S(j, 1) = chisq_S(j, 1) + bs_tot(1, 1 + j, i)**2/As_inv
               do k = 2, ndelta
                  chisq_S(j, k) = chisq_S(j, k) + bs_tot(tod%output_n_maps + k -1, 1 + j, i)**2/As_inv
               end do
            end do
         end if
         rms%map(i, 1) = sqrt(A_inv)*scale
         do k = 1, tod%output_n_maps
            binmap%outmaps(k)%p%map(i, 1) = b_tot(k, 1, i)/A_inv*scale
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
      end if

      deallocate (A_tot, b_tot, bs_tot, W, eta)


  end subroutine finalize_binned_map_unpol


  subroutine finalize_binned_map(tod, binmap, rms, scale, chisq_S, mask)
    !
    ! Routine to finalize the binned maps
    ! 
    ! Arguments:
    ! ----------
    ! tod:
    ! binmap:
    ! rms:
    ! scale
    ! chisq_S
    ! mask
    !
    implicit none
    class(comm_tod),                      intent(in)    :: tod
    type(comm_binmap),                    intent(inout) :: binmap
    class(comm_map),                      intent(inout) :: rms
    real(dp),                             intent(in)    :: scale
    real(dp),        dimension(1:,1:),    intent(out),   optional :: chisq_S
    real(sp),        dimension(0:),       intent(in),    optional :: mask

    integer(i4b) :: i, j, k, nmaps, ierr, ndet, ncol, n_A, off, ndelta
    integer(i4b) :: det, nout, np0, comm, myid, nprocs
    real(dp), allocatable, dimension(:,:)   :: A_inv, As_inv
    real(dp), allocatable, dimension(:,:,:) :: b_tot, bs_tot
    real(dp), allocatable, dimension(:)     :: W, eta
    real(dp), allocatable, dimension(:,:)   :: A_tot

    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    nout  = size(binmap%sb_map%a,dim=1)
    nmaps = tod%info%nmaps
    ndet  = tod%ndet
    n_A   = size(binmap%sA_map%a,dim=1)
    ncol  = size(binmap%sb_map%a,dim=2)
    ndelta = 0; if (present(chisq_S)) ndelta = size(chisq_S,dim=2)

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

         ! Can I return the condition number?
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
               if (mask(tod%info%pix(i + 1)) == 0.) cycle
               if (As_inv(nmaps + j, nmaps + j) <= 0.d0) cycle
               chisq_S(j, 1) = chisq_S(j, 1) + bs_tot(1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               do k = 2, ndelta
                  chisq_S(j, k) = chisq_S(j, k) + bs_tot(tod%output_n_maps + k - 1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               end do
            end do
         end if

         ! Store map in correct units
         do j = 1, nmaps
            do k = 1, tod%output_n_maps
               binmap%outmaps(k)%p%map(i, j) = b_tot(k, j, i)*scale
            end do
         end do

         ! Store N in correct units
         if (rms%info%nmaps == 3) then
            ! Diagonal matrix; store RMS
            do j = 1, nmaps
               rms%map(i,j) = sqrt(A_inv(j, j))*scale
            end do
         else if (rms%info%nmaps == 4) then           
            ! Block-diagonal T + QU; store N
            rms%map(i,1) = A_inv(1,1)*scale**2
            rms%map(i,2) = A_inv(2,2)*scale**2
            rms%map(i,3) = A_inv(3,3)*scale**2
            rms%map(i,4) = A_inv(2,3)*scale**2
         end if
      end do

      if (present(chisq_S)) then
         if (myid == 0) then
            call mpi_reduce(mpi_in_place, chisq_S, size(chisq_S), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
         else
            call mpi_reduce(chisq_S, chisq_S, size(chisq_S), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
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
     real(dp),         dimension(:, :),    intent(inout) :: bicg_sol
     integer(i4b),                            intent(in) :: npix, nmaps
     integer(i4b),                         intent(inout) :: num_cg_iters
     real(dp),                             intent(inout) :: epsil
     real(sp),                  dimension(:), intent(in) :: procmask
     real(dp),                dimension(:,:), intent(in) :: map_full
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
     real(dp)                                   :: alpha, beta
     real(dp),     allocatable, dimension(:, :) :: r, s, q
     logical(lgt)                               :: finished, write_cg
     real(dp)                                   :: rho_old, rho_new, monopole
     real(dp)                                   :: omega, delta_r, delta_s
     real(dp),     allocatable, dimension(:, :) :: rhat, r0, shat, p, phat, v, x_temp, y_temp
     real(dp),        allocatable, dimension(:) :: determ
     character(len=512)                         :: i_str, l_str

     call timer%start(TOD_MAPSOLVE, tod%band)
     write_cg = .false.
     !write_cg = .true.
     !write_cg = tod%first_call

     if (comp_S) then
        allocate (x_temp(nmaps+1,0:npix-1))
        allocate (y_temp(nmaps+1,0:npix-1))
     else
        allocate (x_temp(nmaps,0:npix-1))
        allocate (y_temp(nmaps,0:npix-1))
     end if
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


        call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, x_temp, y_temp, bicg_sol, r)
        r = b_map(:, :, l) - r
        monopole = sum(b_map(:,1,l)*M_diag(:,1)*procmask) &
               & / sum(M_diag(:,1)*procmask)
        !if (l == 1) then
        !  !bicg_sol = transpose(map_full)
        !else
        !  bicg_sol = 0d0
        !end if
        r0 = b_map(:, :, l) - monopole
        call tod%apply_map_precond(r0, rhat)
        
        delta_r = sum(r0*rhat)
        delta_0 = delta_r
        delta_s = delta_s

        omega = 1d0
        alpha = 1d0

        rho_new = sum(r0*r)
        i = 0

        if (write_cg) then
          write(i_str, '(I0.3)') 0
          write(l_str, '(I1)') l
          call write_map(trim(prefix)//'cgest_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                       & bicg_sol(:,1:3))
          call write_map(trim(prefix)//'cgres_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                       & r(:, 1:3))
        end if
        bicg: do
           i = i + 1
           rho_old = rho_new
           call update_status(status, 'dot product')
           rho_new = sum(r0*r)
           call update_status(status, 'done dot product')
           if (rho_new == 0d0) then
             if (tod%verbosity > 1) write(*,*) '|      Residual norm is zero'
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

           call tod%apply_map_precond(p, phat)
           
           call update_status(status, 'v=A phat')
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, x_temp, y_temp, phat, v)
           call update_status(status, 'done')
           num_cg_iters = num_cg_iters + 1

           alpha         = rho_new/sum(r0*v)
           s             = r - alpha*v
           call tod%apply_map_precond(s, shat)
           delta_s       = sum(s*shat)

           if (tod%verbosity > 1) then 
              write(*,101) 2*i-1, delta_s/delta_0
101           format (' |', 6X, I4, ':   delta_s/delta_0:',  2X, ES11.4)
           end if

           bicg_sol = bicg_sol + alpha*phat

           if (write_cg) then
             write(i_str, '(I0.3)') 2*i-1
             write(l_str, '(I1)') l
             call write_map(trim(prefix)//'cgest_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & bicg_sol(:,1:3))
             call write_map(trim(prefix)//'cgres_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & s(:, 1:3))
           end if

           if (abs(delta_s) .le. (delta_0*epsil) .and. 2*i-1 .ge. i_min) then
              if (tod%verbosity > 1) write(*,*) '|      Reached bicg-stab tolerance'
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if

           call update_status(status, 'q=A shat')
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, x_temp, y_temp, shat, q)
           call update_status(status, 'done')

           omega         = sum(q*s)/sum(q*q)
           bicg_sol = bicg_sol + omega*shat


           if (omega == 0d0) then
             if (tod%verbosity > 1) write(*,*) '|      omega is zero'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if

           if (mod(i, recomp_freq) == 1 .or. beta > 1.d8) then
              call update_status(status, 'A xhat')
              call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, x_temp, y_temp, bicg_sol, r)
              call update_status(status, 'done')
              r(:,1) = b_map(:,1,l)  - r(:,1) - monopole
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
                          & bicg_sol(:,1:3))
             call write_map(trim(prefix)//'cgres_'//trim(i_str)//'_'//trim(l_str)//trim(postfix), &
                          & r(:, 1:3))
           end if

           call tod%apply_map_precond(r, rhat)
           delta_r      = sum(r*rhat)
           num_cg_iters = num_cg_iters + 1

           if (tod%verbosity > 1) then 
              write(*,102) 2*i, delta_r/delta_0
102           format (' |', 6X, I4, ':   delta_r/delta_0:',  2X, ES11.4)
           end if
           if (abs(delta_r) .le. delta_0*epsil .and. 2*i .ge. i_min) then
              if (tod%verbosity > 1) write(*,*) '|      Reached bicg-stab tolerance'
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           else if (delta_r > delta_0*1000) then
              write(*,*) '|      Solution is diverging, killing search'
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if
           if (i==i_max) then
             if (tod%verbosity > 1) write(*,*) '|      Reached maximum number of iterations'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if
        end do bicg

     else
        loop: do while (.true.) 
           call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
           if (finished) exit loop
           call compute_Ax(tod, tod%x_im, procmask, comp_S, M_diag, x_temp, y_temp)
        end do loop
     end if
     if (tod%myid == 0) deallocate (r, rhat, s, r0, q, shat, p, phat, v, m_buf)
     deallocate (x_temp, y_temp)
     if (tod%myid == 0 .and. .not. comp_S) deallocate (determ)

     call timer%stop(TOD_MAPSOLVE, tod%band)

   end subroutine run_bicgstab


end module comm_tod_mapmaking_mod
