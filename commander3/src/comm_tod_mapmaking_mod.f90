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
   use comm_shared_arr_mod
   use comm_map_mod
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


  end subroutine init_binmap


  subroutine dealloc_binmap(self)
    implicit none
    class(comm_binmap), intent(inout) :: self

    integer(i4b) ::  i


  end subroutine dealloc_binmap

  subroutine synchronize_binmap(self, tod)
    implicit none
    class(comm_binmap),  intent(inout) :: self
    class(comm_tod),     intent(in)    :: tod

    integer(i4b) :: i, j, start_chunk, end_chunk, ierr


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


end subroutine bin_differential_TOD

   subroutine compute_Ax(tod, x_imarr, pmask, comp_S, M_diag, split, x, y, x_in, y_out)
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
      integer(i4b), intent(in)         :: split
      real(dp),                dimension(1:,0:), intent(inout), optional :: x
      real(dp),                dimension(1:,0:), intent(inout), optional :: y
      real(dp),     dimension(0:, 1:), intent(in),    optional :: x_in
      real(dp),     dimension(0:, 1:), intent(inout), optional :: y_out

      integer(i4b), allocatable, dimension(:)         :: flag
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi
      !integer(i4b), dimension(10) :: t_arr=(/52131,52496,52861,53227,53592,53957,54322, &
      !                                    &  54688,55053,55418/)


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


   end subroutine finalize_binned_map

   subroutine run_bicgstab(tod, handle, bicg_sol, npix, nmaps, num_cg_iters, epsil, procmask, map_full, M_diag, b_map, l, prefix, postfix, comp_S, split)
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
     integer(i4b), intent(in)                            :: split





   end subroutine run_bicgstab


end module comm_tod_mapmaking_mod
