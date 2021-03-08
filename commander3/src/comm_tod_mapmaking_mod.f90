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
   implicit none

contains


  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine bin_TOD(tod, data, pix, psi, flag, A, b, scan, comp_S, b_mono)
    implicit none
    class(comm_tod),                             intent(in)    :: tod
    integer(i4b),                                intent(in)    :: scan
    real(sp),            dimension(1:,1:,1:),    intent(in)    :: data
    integer(i4b),        dimension(1:,1:),       intent(in)    :: pix, psi, flag
    real(dp),            dimension(1:,1:),       intent(inout) :: A
    real(dp),            dimension(1:,1:,1:),    intent(inout) :: b
    real(dp),            dimension(1:,1:,1:),    intent(inout), optional :: b_mono
    logical(lgt),                                intent(in)    :: comp_S

    integer(i4b) :: det, i, t, pix_, off, nout, psi_
    real(dp)     :: inv_sigmasq

    nout        = size(b,dim=1)

    do det = 1, tod%ndet
       if (.not. tod%scans(scan)%d(det)%accept) cycle
       !write(*,*) tod%scanid(scan), det, tod%scans(scan)%d(det)%sigma0
       off         = 6 + 4*(det-1)
       inv_sigmasq = (tod%scans(scan)%d(det)%gain/tod%scans(scan)%d(det)%sigma0)**2
       do t = 1, tod%scans(scan)%ntod
          
          if (iand(flag(t,det),tod%flag0) .ne. 0) cycle
         !if (flag(t,det)==1) cycle
          
          pix_    = tod%pix2ind(pix(t,det))
          psi_    = psi(t,det)
          
          A(1,pix_) = A(1,pix_) + 1.d0                                  * inv_sigmasq
          A(2,pix_) = A(2,pix_) + tod%cos2psi(psi_)                    * inv_sigmasq
          A(3,pix_) = A(3,pix_) + tod%cos2psi(psi_)**2                 * inv_sigmasq
          A(4,pix_) = A(4,pix_) + tod%sin2psi(psi_)                    * inv_sigmasq
          A(5,pix_) = A(5,pix_) + tod%cos2psi(psi_)*tod%sin2psi(psi_) * inv_sigmasq
          A(6,pix_) = A(6,pix_) + tod%sin2psi(psi_)**2                 * inv_sigmasq
          
          do i = 1, nout
             b(i,1,pix_) = b(i,1,pix_) + data(i,t,det)                      * inv_sigmasq
             b(i,2,pix_) = b(i,2,pix_) + data(i,t,det) * tod%cos2psi(psi_) * inv_sigmasq
             b(i,3,pix_) = b(i,3,pix_) + data(i,t,det) * tod%sin2psi(psi_) * inv_sigmasq
          end do
          
          if (present(b_mono)) then
             b_mono(1,pix_,det) = b_mono(1,pix_,det) +                      inv_sigmasq
             b_mono(2,pix_,det) = b_mono(2,pix_,det) + tod%cos2psi(psi_) * inv_sigmasq
             b_mono(3,pix_,det) = b_mono(3,pix_,det) + tod%sin2psi(psi_) * inv_sigmasq
          end if
          
          if (comp_S .and. det < tod%ndet) then
             A(off+1,pix_) = A(off+1,pix_) + 1.d0               * inv_sigmasq 
             A(off+2,pix_) = A(off+2,pix_) + tod%cos2psi(psi_) * inv_sigmasq
             A(off+3,pix_) = A(off+3,pix_) + tod%sin2psi(psi_) * inv_sigmasq
             A(off+4,pix_) = A(off+4,pix_) + 1.d0               * inv_sigmasq
             do i = 1, nout
                b(i,det+3,pix_) = b(i,det+3,pix_) + data(i,t,det) * inv_sigmasq 
             end do
          end if
          
       end do
    end do

  end subroutine bin_TOD


   ! differential TOD computation, written with WMAP in mind.
   subroutine bin_differential_TOD(tod, data, pix, psi, flag, x_imarr, pmask, b, M_diag, scan, comp_S, b_mono)
      implicit none
      class(comm_tod), intent(in)                               :: tod
      integer(i4b), intent(in)                                  :: scan
      real(sp), dimension(1:, 1:, 1:), intent(in)               :: data
      integer(i4b), dimension(1:), intent(in)                   :: flag
      integer(i4b), dimension(0:), intent(in)                   :: pmask
      integer(i4b), dimension(1:, 1:), intent(in)               :: pix, psi
      real(dp), dimension(1:), intent(in)                       :: x_imarr
      real(dp), dimension(0:, 1:, 1:), intent(inout)            :: b
      real(dp), dimension(0:, 1:), intent(inout)                :: M_diag
      real(dp), dimension(0:, 1:, 1:), intent(inout), optional  :: b_mono
      logical(lgt), intent(in)                                  :: comp_S

      integer(i4b) :: det, i, t, nout
      real(dp)     :: inv_sigmasq, d, p, var, dx_im, x_im

      integer(i4b) :: lpix, rpix, lpsi, rpsi, sgn

      integer(i4b) :: pA, pB, f_A, f_B

      nout = size(b, dim=3)
      dx_im = 0.5*(x_imarr(1) - x_imarr(2))
      x_im = 0.5*(x_imarr(1) + x_imarr(2))

      if (tod%scans(scan)%d(1)%accept) then
         !inv_sigmasq = 0.d0 
         var = 0
         do det = 1, 4
           var = var  + (tod%scans(scan)%d(det)%sigma0/tod%scans(scan)%d(det)%gain)**2/4
           !inv_sigmasq = inv_sigmasq  + (tod%scans(scan)%d(det)%gain/tod%scans(scan)%d(det)%sigma0)**2
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
               do det = 1, 4
                 !sgn = (-1)**((det + 1)/2 + 1) ! 1 for 13, 14, -1 for 23, 24
                 d = d + data(i, t, det)/4
                 p = p + data(i, t, det)/4*(-1)**((det + 1)/2 + 1)
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
            end do

            M_diag(lpix, 1) = M_diag(lpix, 1) + f_A*inv_sigmasq
            M_diag(rpix, 1) = M_diag(rpix, 1) + f_B*inv_sigmasq
            M_diag(lpix, 2) = M_diag(lpix, 2) + f_A*inv_sigmasq*tod%cos2psi(lpsi)**2
            M_diag(rpix, 2) = M_diag(rpix, 2) + f_B*inv_sigmasq*tod%cos2psi(rpsi)**2
            M_diag(lpix, 3) = M_diag(lpix, 3) + f_A*inv_sigmasq*tod%sin2psi(lpsi)**2
            M_diag(rpix, 3) = M_diag(rpix, 3) + f_B*inv_sigmasq*tod%sin2psi(rpsi)**2

            ! Not a true diagonal term, just the off-diagonal estimate of the
            ! covariance for each pixel.
            M_diag(lpix, 4) = M_diag(lpix, 4)+f_A*inv_sigmasq*tod%sin2psi(lpsi)*tod%cos2psi(lpsi)
            M_diag(rpix, 4) = M_diag(rpix, 4)+f_B*inv_sigmasq*tod%sin2psi(rpsi)*tod%cos2psi(rpsi)

         end do
       end if

end subroutine bin_differential_TOD

   subroutine compute_Ax(tod, x_imarr, pmask, x_in, y_out)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Code to compute matrix product P^T N^-1 P m
      ! y = Ax
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      class(comm_tod),                 intent(in)              :: tod
      real(dp),     dimension(1:),     intent(in)              :: x_imarr
      integer(i4b), dimension(0:),     intent(in)              :: pmask
      real(dp),     dimension(0:, 1:), intent(in),    optional :: x_in
      real(dp),     dimension(0:, 1:), intent(inout), optional :: y_out

      integer(i4b), allocatable, dimension(:)         :: flag
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi

      logical(lgt) :: finished
      integer(i4b) :: i, j, k, ntod, ndet, lpix, rpix, lpsi, rpsi, ierr
      integer(i4b) :: nhorn, t, sgn, pA, pB, f_A, f_B, nside, npix, nmaps
      real(dp)     :: inv_sigmasq, var, dA, dB, iA, iB, sA, sB, d, p, x_im, dx_im
      real(dp), allocatable, dimension(:,:) :: x, y
      nhorn = tod%nhorn
      ndet  = tod%ndet
      nside = tod%nside
      nmaps = tod%nmaps
      npix  = 12*nside**2

      allocate(x(0:npix-1,nmaps), y(0:npix-1,nmaps))
      if (tod%myid == 0) then
         finished = .false.
         call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
         x = x_in
      end if
      call mpi_bcast(x, size(x),  MPI_DOUBLE_PRECISION, 0, tod%info%comm, ierr)

      x_im   = 0.5*(x_imarr(1) + x_imarr(2))
      dx_im  = 0.5*(x_imarr(1) - x_imarr(2))
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
               var = var + (tod%scans(j)%d(k)%sigma0/tod%scans(j)%d(k)%gain)**2/4
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
               ! The sgn parameter is +1 for timestreams 13 and 14, -1
               ! for timestreams 23 and 24, and also is used to switch
               ! the sign of the polarization sensitive parts of the
               ! model
               iA = x(lpix, 1)
               iB = x(rpix, 1)
               sA = x(lpix, 2)*tod%cos2psi(lpsi) + x(lpix, 3)*tod%sin2psi(lpsi)
               sB = x(rpix, 2)*tod%cos2psi(rpsi) + x(rpix, 3)*tod%sin2psi(rpsi)
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
            end do
         end if
         deallocate (pix, psi, flag)
      end do

      if (tod%myid == 0) then
         call mpi_reduce(y, y_out, size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
      else
         call mpi_reduce(y, y,     size(y), MPI_DOUBLE_PRECISION,MPI_SUM,&
              & 0, tod%info%comm, ierr)
      end if

      deallocate(x, y)

   end subroutine compute_Ax

  subroutine finalize_binned_map(tod, handle, sA_map, sb_map, rms, outmaps, chisq_S, Sfile, mask, sb_mono, sys_mono, condmap)
    implicit none
    class(comm_tod),                      intent(in)    :: tod
    type(planck_rng),                     intent(inout) :: handle
    type(shared_2d_dp), intent(inout) :: sA_map
    type(shared_3d_dp), intent(inout) :: sb_map
    class(comm_map),                      intent(inout) :: rms
    class(map_ptr),  dimension(1:),       intent(inout), optional :: outmaps
    real(dp),        dimension(1:,1:),    intent(out),   optional :: chisq_S
    character(len=*),                     intent(in),    optional :: Sfile
    integer(i4b),    dimension(0:),       intent(in),    optional :: mask
    type(shared_3d_dp), intent(inout), optional :: sb_mono
    real(dp),        dimension(1:,1:,0:), intent(out),   optional :: sys_mono
    class(comm_map),                      intent(inout), optional :: condmap

    integer(i4b) :: i, j, k, nmaps, ierr, ndet, ncol, n_A, off, ndelta
    integer(i4b) :: det, nout, np0, comm, myid, nprocs
    real(dp), allocatable, dimension(:,:)   :: A_inv, As_inv, buff_2d
    real(dp), allocatable, dimension(:,:,:) :: b_tot, bs_tot, buff_3d
    real(dp), allocatable, dimension(:)     :: W, eta
    real(dp), allocatable, dimension(:,:)   :: A_tot
    class(comm_mapinfo), pointer :: info 
    class(comm_map), pointer :: smap 

    myid  = tod%myid
    nprocs= tod%numprocs
    comm  = tod%comm
    np0   = tod%info%np
    nout  = size(sb_map%a,dim=1)
    nmaps = tod%info%nmaps
    ndet  = tod%ndet
    n_A   = size(sA_map%a,dim=1)
    ncol  = size(sb_map%a,dim=2)
    ndelta = 0; if (present(chisq_S)) ndelta = size(chisq_S,dim=2)

    ! Collect contributions from all nodes
    !call update_status(status, "tod_final1")
    call mpi_win_fence(0, sA_map%win, ierr)
    if (sA_map%myid_shared == 0) then
       !allocate(buff_2d(size(sA_map%a,1),size(sA_map%a,2)))
       

!       call mpi_allreduce(sA_map%a, buff_2d, size(sA_map%a), &
!            & MPI_DOUBLE_PRECISION, MPI_SUM, sA_map%comm_inter, ierr)
         do i = 1, size(sA_map%a, 1)
            call mpi_allreduce(MPI_IN_PLACE, sA_map%a(i, :), size(sA_map%a, 2), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, sA_map%comm_inter, ierr)
         end do
         !sA_map%a = buff_2d
         !deallocate(buff_2d)
      end if
      call mpi_win_fence(0, sA_map%win, ierr)
      call mpi_win_fence(0, sb_map%win, ierr)
      if (sb_map%myid_shared == 0) then
         !allocate(buff_3d(size(sb_map%a,1),size(sb_map%a,2),size(sb_map%a,3)))
         !call mpi_allreduce(sb_map%a, buff_3d, size(sb_map%a), &
         !     & MPI_DOUBLE_PRECISION, MPI_SUM, sb_map%comm_inter, ierr)
         do i = 1, size(sb_map%a, 1)
            call mpi_allreduce(mpi_in_place, sb_map%a(i, :, :), size(sb_map%a(1, :, :)), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, sb_map%comm_inter, ierr)
         end do
         !sb_map%a = buff_3d
         !deallocate(buff_3d)
      end if
      call mpi_win_fence(0, sb_map%win, ierr)
      if (present(sb_mono)) then
         call mpi_win_fence(0, sb_mono%win, ierr)
         if (sb_mono%myid_shared == 0) then
            call mpi_allreduce(MPI_IN_PLACE, sb_mono%a, size(sb_mono%a), &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, sb_mono%comm_inter, ierr)
         end if
         call mpi_win_fence(0, sb_mono%win, ierr)
      end if

      allocate (A_tot(n_A, 0:np0 - 1), b_tot(nout, nmaps, 0:np0 - 1), bs_tot(nout, ncol, 0:np0 - 1), W(nmaps), eta(nmaps))
      A_tot = sA_map%a(:, tod%info%pix + 1)
      b_tot = sb_map%a(:, 1:nmaps, tod%info%pix + 1)
      bs_tot = sb_map%a(:, :, tod%info%pix + 1)
      if (present(sb_mono)) then
         do j = 1, ndet
            sys_mono(:, nmaps + j, :) = sb_mono%a(:, tod%info%pix + 1, j)
         end do
      end if
      !call update_status(status, "tod_final2")

      ! Solve for local map and rms
      if (present(Sfile)) then
         info => comm_mapinfo(tod%comm, tod%info%nside, 0, ndet - 1, .false.)
         smap => comm_map(info)
         smap%map = 0.d0
      end if
      allocate (A_inv(nmaps, nmaps), As_inv(ncol, ncol))
      !if (present(chisq_S)) smap%map = 0.d0
      if (present(chisq_S)) chisq_S = 0.d0
      do i = 0, np0 - 1
         if (all(b_tot(1, :, i) == 0.d0)) then
            !write(*,*) 'missing pixel',i, tod%myid
            if (present(sb_mono)) sys_mono(1:nmaps, 1:nmaps, i) = 0.d0
            if (.not. present(chisq_S)) then
               rms%map(i, :) = 0.d0
               if (present(outmaps)) then
                  do k = 1, nout
                     outmaps(k)%p%map(i, :) = 0.d0
                  end do
               end if
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

         if (present(condmap)) then
            call get_eigenvalues(A_inv(1:nmaps, 1:nmaps), W)
            if (minval(W) < 0.d0) then
!             if (maxval(W) == 0.d0 .or. minval(W) == 0.d0) then
               write (*, *) 'A_inv: ', A_inv
               write (*, *) 'W', W
!             end if
            end if
            if (minval(W) > 0) then
               condmap%map(i, 1) = log10(max(abs(maxval(W)/minval(W)), 1.d0))
            else
               condmap%map(i, 1) = 30.d0
            end if
         end if

         call invert_singular_matrix(A_inv, 1d-12)
         do k = 1, tod%output_n_maps
            b_tot(k, 1:nmaps, i) = matmul(A_inv, b_tot(k, 1:nmaps, i))
         end do
!!$       if (tod%info%pix(i) == 100) then
!!$          write(*,*) 'b', b_tot(:,:,i)
!!$       end if
         if (present(chisq_S)) then
            call invert_singular_matrix(As_inv, 1d-12)
            bs_tot(1, 1:ncol, i) = matmul(As_inv, bs_tot(1, 1:ncol, i))
            do k = tod%output_n_maps + 1, nout
               bs_tot(k, 1:ncol, i) = matmul(As_inv, bs_tot(k, 1:ncol, i))
            end do
         end if

         if (present(sb_mono)) sys_mono(1:nmaps, 1:nmaps, i) = A_inv(1:nmaps, 1:nmaps)
         if (present(chisq_S)) then
            do j = 1, ndet - 1
               if (mask(tod%info%pix(i + 1)) == 0) cycle
               if (As_inv(nmaps + j, nmaps + j) <= 0.d0) cycle
               chisq_S(j, 1) = chisq_S(j, 1) + bs_tot(1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               do k = 2, ndelta
                  chisq_S(j, k) = chisq_S(j, k) + bs_tot(tod%output_n_maps + k - 1, nmaps + j, i)**2/As_inv(nmaps + j, nmaps + j)
               end do
               if (present(Sfile) .and. As_inv(nmaps + j, nmaps + j) > 0.d0) then
                  smap%map(i, j) = bs_tot(1, nmaps + j, i)/sqrt(As_inv(nmaps + j, nmaps + j))
               end if
            end do
         end if
         do j = 1, nmaps
            rms%map(i, j) = sqrt(A_inv(j, j))*1.d6 ! uK
            if (present(outmaps)) then
               do k = 1, tod%output_n_maps
                  outmaps(k)%p%map(i, j) = b_tot(k, j, i)*1.d6 ! uK
               end do
            end if
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
         if (present(Sfile)) call smap%writeFITS(trim(Sfile))
      end if

      if (present(Sfile)) call smap%dealloc

      deallocate (A_inv, As_inv, A_tot, b_tot, bs_tot, W, eta)

   end subroutine finalize_binned_map


   subroutine run_bicgstab(tod, handle, cg_sol, npix, nmaps, num_cg_iters, epsil, procmask, map_full, M_diag, b_map, l)
     implicit none
     class(comm_tod), intent(in) :: tod
     type(planck_rng),                     intent(inout) :: handle
     real(dp), dimension(:, :, :), intent(inout) :: cg_sol
     integer(i4b), intent(in) :: npix, nmaps
     integer(i4b), intent(inout) :: num_cg_iters
     real(dp),  intent(in) :: epsil
     integer(i4b), dimension(:), intent(in)    :: procmask
     real(dp), dimension(:), intent(inout)     :: map_full
     real(dp), dimension(:,:), intent(in)     :: M_diag
     real(dp), dimension(:,:,:), intent(in)     :: b_map
     integer(i4b), intent(in)    :: l

     real(dp), allocatable, dimension(:, :)          :: m_buf
     integer(i4b) :: i_max, i_min, ierr, status, i
     real(dp) :: delta_0, delta_old, delta_new
     real(dp) :: alpha, beta, g, f_quad, sigma_mono
     real(dp), allocatable, dimension(:, :)    :: r, s, d, q
     real(dp) :: monopole
     logical(lgt) :: write_cg_iter=.false., finished
     real(dp) :: rho_old, rho_new
     real(dp) :: omega, delta_r, delta_s
     real(dp), allocatable, dimension(:, :) :: rhat, r0, shat, p, phat, v
     real(dp), allocatable, dimension(:)    :: determ

     if (tod%myid==0) then
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
         determ = M_diag(:,2)*M_diag(:,3) - M_diag(:,4)**2

         i_max = 100
         i_min = 0

        if (.false. .and. l == 1) then
           call compute_Ax(tod, tod%x_im, procmask, cg_sol(:,:,1), v)
           r = b_map(:, :, l) - v 
        else
           r  = b_map(:, :, l)
        end if
        r0 = b_map(:, :, l)
        rhat(:,1) =  r(:,1)/M_diag(:,1)
        rhat(:,2) = (r(:,2)*M_diag(:,3)- r(:,2)*M_diag(:,4))/determ
        rhat(:,3) = (r(:,3)*M_diag(:,2)- r(:,3)*M_diag(:,4))/determ
        delta_r = sum(r*rhat)
        delta_0 = delta_r
        delta_s = delta_s

        omega = 1d0
        alpha = 1d0

        rho_new = sum(r0*r)
        bicg: do i = 1, i_max
           rho_old = rho_new
           rho_new = sum(r0*r)
           if (rho_new == 0d0) then
             write(*,*) 'rho_i is zero'
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
           phat(:,1) =  p(:,1)/M_diag(:,1)
           phat(:,2) = (p(:,2)*M_diag(:,3)- p(:,2)*M_diag(:,4))/determ
           phat(:,3) = (p(:,3)*M_diag(:,2)- p(:,3)*M_diag(:,4))/determ
           
           call compute_Ax(tod, tod%x_im, procmask, phat, v)
           num_cg_iters = num_cg_iters + 1

           alpha         = rho_new/sum(r0*v)
           s             = r - alpha*v
           shat(:,1) =  s(:,1)/M_diag(:,1)
           shat(:,2) = (s(:,2)*M_diag(:,3)- s(:,2)*M_diag(:,4))/determ
           shat(:,3) = (s(:,3)*M_diag(:,2)- s(:,3)*M_diag(:,4))/determ
           delta_s       = sum(s*shat)

           if (tod%verbosity > 1) then 
              write(*,101) 2*i-1, delta_s/delta_0
101           format (6X, I4, ':   delta_s/delta_0:',  2X, ES9.2)
           end if

           cg_sol(:,:,l) = cg_sol(:,:,l) + alpha*phat

           if (delta_s .le. (delta_0*epsil) .and. 2*i-1 .ge. i_min) then
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if

           call compute_Ax(tod, tod%x_im, procmask, shat, q)

           omega         = sum(q*s)/sum(q*q)
           cg_sol(:,:,l) = cg_sol(:,:,l) + omega*shat
           if (omega == 0d0) then
             write(*,*) 'omega is zero'
             finished = .true.
             call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
             exit bicg
           end if

           if (mod(i, 10) == 1 .or. beta > 1.d8) then
              call compute_Ax(tod, tod%x_im, procmask, cg_sol(:,:,l), r)
              r = b_map(:, :, l) - r
           else
              r = s - omega*q
           end if

           rhat(:,1) =  r(:,1)/M_diag(:,1)
           rhat(:,2) = (r(:,2)*M_diag(:,3)- r(:,2)*M_diag(:,4))/determ
           rhat(:,3) = (r(:,3)*M_diag(:,2)- r(:,3)*M_diag(:,4))/determ
           delta_r      = sum(r*rhat)
           num_cg_iters = num_cg_iters + 1

           if (tod%verbosity > 1) then 
              write(*,102) 2*i, delta_r/delta_0
102           format (6X, I4, ':   delta_r/delta_0:',  2X, ES9.2)
           end if
           if (delta_r .le. delta_0*epsil .and. 2*i .ge. i_min) then
              finished = .true.
              call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
              exit bicg
           end if
        end do bicg

        if (l == 1) then
           ! Maximum likelihood monopole
           monopole = sum((cg_sol(:,1,1)-map_full)*M_diag(:,1)*procmask) &
                  & / sum(M_diag(:,1)*procmask)
           write(*,*) monopole
           if (trim(tod%operation) == 'sample') then
              ! Add fluctuation term if requested
              sigma_mono = sum(M_diag(:,1) * procmask)
              if (sigma_mono > 0.d0) sigma_mono = 1.d0 / sqrt(sigma_mono)
              if (tod%verbosity > 1) then
                write(*,*) 'monopole, fluctuation sigma'
                write(*,*) monopole, sigma_mono
              end if
              monopole = monopole + sigma_mono * rand_gauss(handle)
           end if
           write(*,*) 'sampled final monopole = ', monopole
           cg_sol(:,1,1) = cg_sol(:,1,1) - monopole
           write(*,*) 'cg_sol monopole', sum(cg_sol(:,1,1)*procmask)/sum(procmask)
        end if
     else
        loop: do while (.true.) 
           call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, tod%info%comm, ierr)
           if (finished) exit loop
           call compute_Ax(tod, tod%x_im, procmask)
        end do loop
     end if
     if (tod%myid == 0) deallocate (r, rhat, s, r0, q, shat, p, phat, v, m_buf, determ)

   end subroutine run_bicgstab

   subroutine sample_mono(tod, handle, sys_mono, res, rms, mask)
      implicit none
      class(comm_tod), intent(inout) :: tod
      type(planck_rng), intent(inout) :: handle
      real(dp), dimension(1:, 1:, 0:), intent(in)    :: sys_mono
      class(comm_map), intent(in)    :: res, rms, mask

      integer(i4b) :: i, j, k, nstep, nmaps, ndet, ierr
      real(dp)     :: t1, t2, s(3), b(3), alpha, sigma, chisq, chisq0, chisq_prop, naccept
      logical(lgt) :: accept, first_call
      real(dp), allocatable, dimension(:)   :: mono0, mono_prop, mu, eta
      real(dp), allocatable, dimension(:, :) :: samples

      call wall_time(t1)

      first_call = .not. allocated(tod%L_prop_mono)
      sigma = 0.03d-6 ! RMS proposal in K
      alpha = 1.2d0   ! Covariance matrix scaling factor
      nmaps = size(sys_mono, dim=1)
      ndet = size(sys_mono, dim=2) - nmaps
      allocate (mono0(ndet), mono_prop(ndet), eta(ndet), mu(ndet))

      ! Initialize monopoles on existing values
      if (tod%myid == 0) then
         do j = 1, ndet
            mono0(j) = tod%mono(j)
         end do

         if (first_call) then
            allocate (tod%L_prop_mono(ndet, ndet))
            tod%L_prop_mono = 0.d0
            do i = 1, ndet
               tod%L_prop_mono(i, i) = sigma
            end do
            nstep = 100000
            allocate (samples(ndet, nstep))
         else
            nstep = 1000
         end if
      end if
      call mpi_bcast(mono0, ndet, MPI_DOUBLE_PRECISION, 0, res%info%comm, ierr)
      call mpi_bcast(nstep, 1, MPI_INTEGER, 0, res%info%comm, ierr)

      ! Compute chisquare for old values; only include Q and U in evaluation;
      ! add old correction to output map
      chisq = 0.d0
      do k = 0, res%info%np - 1
         b = 0.d0
         do j = 1, ndet
            b = b + mono0(j)*sys_mono(:, nmaps + j, k)
         end do
         s = matmul(sys_mono(1:3, 1:3, k), b)*1.d6 ! uK
         if (mask%map(k, 1) == 1.d0) then
            chisq = chisq + sum((res%map(k, 2:3) - s(2:3))**2/rms%map(k, 2:3)**2)
         end if
      end do
      call mpi_reduce(chisq, chisq0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, res%info%comm, ierr)

      naccept = 0
      if (res%info%myid == 0) open (78, file='mono.dat', recl=1024)
      do i = 1, nstep

         ! Propose new monopoles; force zero average
         if (res%info%myid == 0) then
            do j = 1, ndet
               eta(j) = rand_gauss(handle)
            end do
            mono_prop = mono0 + matmul(tod%L_prop_mono, eta)
            mono_prop = mono_prop - mean(mono_prop)
         end if
         call mpi_bcast(mono_prop, ndet, MPI_DOUBLE_PRECISION, 0, res%info%comm, ierr)

         ! Compute chisquare for new values; only include Q and U in evaluation
         chisq = 0.d0
         do k = 0, res%info%np - 1
            b = 0.d0
            do j = 1, ndet
               b = b + mono_prop(j)*sys_mono(:, nmaps + j, k)
            end do
            s = matmul(sys_mono(1:3, 1:3, k), b)*1.d6 ! uK
            if (mask%map(k, 1) == 1) then
               chisq = chisq + sum((res%map(k, 2:3) - s(2:3))**2/rms%map(k, 2:3)**2)
            end if
         end do
         call mpi_reduce(chisq, chisq_prop, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, res%info%comm, ierr)

         ! Apply Metropolis rule or maximize
         if (res%info%myid == 0) then
            if (trim(tod%operation) == 'optimize') then
               accept = chisq_prop <= chisq0
            else
               accept = (rand_uni(handle) < exp(-0.5d0*(chisq_prop - chisq0)))
            end if
!!$          write(*,*) 'Mono chisq0 = ', chisq0, ', delta = ', chisq_prop-chisq0
!!$          write(*,*) 'm0 = ', real(mono0,sp)
!!$          write(*,*) 'm1 = ', real(mono_prop,sp)
         end if
         call mpi_bcast(accept, 1, MPI_LOGICAL, 0, res%info%comm, ierr)

         if (accept) then
            mono0 = mono_prop
            chisq0 = chisq_prop
            naccept = naccept + 1
         end if
         if (res%info%myid == 0 .and. mod(i, 100) == 0) write (78, *) i, chisq0, mono0
         if (first_call) samples(:, i) = mono0

      end do
      if (res%info%myid == 0) close (78)

      if (first_call .and. res%info%myid == 0) then
         do i = 1, ndet
            mu(i) = mean(samples(i, :))
            do j = 1, i
               tod%L_prop_mono(i, j) = mean((samples(i, :) - mu(i))*(samples(j, :) - mu(j)))
            end do
!          write(*,*) real(tod%L_prop_mono(i,:),sp)
         end do
         call compute_hermitian_root(tod%L_prop_mono, 0.5d0)
         !call cholesky_decompose_single(tod%L_prop_mono)
         tod%L_prop_mono = alpha*tod%L_prop_mono
         deallocate (samples)
      end if

      ! Update parameters
      tod%mono = mono0

      call wall_time(t2)
      if (res%info%myid == 0) then
         write (*, fmt='(a,f8.3,a,f8.3)') 'Time for monopole sampling = ', t2 - t1, ', accept = ', naccept/nstep
      end if

      deallocate (mono0, mono_prop, eta, mu)

   end subroutine sample_mono

end module comm_tod_mapmaking_mod
