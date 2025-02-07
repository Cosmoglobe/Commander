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
module comm_N_rms_QUcov_mod
  use comm_N_mod
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This module handles the case of a white noise matrix where there is
  ! correlation between the Q and U Stokes parameters, but each pixel is
  ! independent. Explicitly, each matrix is represented mathematically as
  ! N^{-1}_{pp} = [\sum 1/\sigma^2, 0,                             0                               ]
  !               [0,              \sum \cos^2 2\psi/\sigma^2,     \sum \cos2\psi\sin2\psi/\sigma^2]
  !               [0,              \sum\cos2\psi\sin2\psi/\sigma^2 \sum \sin^2 2\psi/\sigma^2      ]
  ! If there were no off-diagonal term, the diagonal terms would be 1/\sigma_{T/Q/U}^2.
  !
  ! The N_map object contains four maps, namely the TT, QQ, UU, and QU elements of the per-pixel covariance
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private
  public comm_N_rms_QUcov, comm_N_rms_QUcov_ptr
  
  type, extends (comm_N) :: comm_N_rms_QUcov
     real(dp),        allocatable, dimension(:,:) :: siN, iN, iN_low
     class(comm_map), pointer                     :: N_map => null()
     class(comm_map), pointer                     :: N_low => null()
   contains
     ! Data procedures
     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map_lowres
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS_rms_QUcov
     procedure :: rms_pix     => returnRMS_rms_QUcov_pix
     procedure :: update_N    => update_N_rms_QUcov
  end type comm_N_rms_QUcov

  interface comm_N_rms_QUcov
     procedure constructor
  end interface comm_N_rms_QUcov

  type comm_N_rms_QUcov_ptr
     type(comm_N_rms_QUcov), pointer :: p => null()
  end type comm_N_rms_QUcov_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask, map)
    implicit none
    class(comm_N_rms_QUcov),                  pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask
    class(comm_map),                    pointer, intent(in), optional :: map

    integer(i4b)       :: i, ierr, nside_smooth
    real(dp)           :: sum_noise, npix
    type(comm_mapinfo), pointer :: info_smooth => null()
    type(comm_mapinfo), pointer :: info_cg => null()
    type(comm_mapinfo), pointer :: info_lowres => null()

    
    ! General parameters
    allocate(constructor)


  end function constructor

  subroutine update_N_rms_QUcov(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_rms_QUcov),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, j, ierr
    real(dp)     :: sum_tau, sum_tau2, sum_noise, npix
    class(comm_map),     pointer :: invW_tau => null(), N_tmp => null()
    class(comm_mapinfo), pointer :: info_lowres => null()
    class(comm_mapinfo), pointer :: info_pre => null()



  end subroutine update_N_rms_QUcov

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i

    if (self%pol) then
       do i = 0, self%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%iN(1,i) * map%map(i,1)
          map%map(i,2) = self%iN(2,i) * buff_Q + self%iN(4,i) * buff_U
          map%map(i,3) = self%iN(4,i) * buff_Q + self%iN(3,i) * buff_U
       end do
    else
       map%map(:,1) = map%map(:,1) * self%iN(1,:)
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_1map_lowres(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i

    if (self%pol) then
       do i = 0, self%N_low%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%iN_low(1,i) * map%map(i,1)
          map%map(i,2) = self%iN_low(2,i) * buff_Q + self%iN_low(4,i) * buff_U
          map%map(i,3) = self%iN_low(4,i) * buff_Q + self%iN_low(3,i) * buff_U
       end do
   else
       map%map(:,1) = self%iN_low(1,:)
   end if
!!$    if (present(samp_group)) then
!!$       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
!!$    end if
  end subroutine matmulInvN_1map_lowres

  ! Return map_out = N * map
  subroutine matmulN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i
    if (self%pol) then
       do i = 0, self%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%N_map%map(i,1) * map%map(i,1)
          map%map(i,2) = self%N_map%map(i,2) * buff_Q + self%N_map%map(i,4) * buff_U
          map%map(i,3) = self%N_map%map(i,4) * buff_Q + self%N_map%map(i,3) * buff_U
       end do
    else
       do i = 0, self%info%np-1
          map%map(i,1) = self%N_map%map(i,1) * map%map(i,1)
       end do
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b)  :: nmaps_band, nmaps_inp, nmaps
    nmaps_band = size(self%siN, dim=2)
    nmaps_inp  = size(map%map, dim=2)
    nmaps      = min(nmaps_inp,nmaps_band)
    map%map(:,1) = self%siN(1,:) * map%map(:,1)
    if (nmaps_inp > nmaps_band) then
      map%map(:,nmaps_band+1:nmaps_inp) = 0
    end if

    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) then
          map%map(:,1) = map%map(:,1) * self%samp_group_mask(samp_group)%p%map(:,1)
          map%map(:,nmaps+1:nmaps_inp) = 0.d0
       end if
    end if
  end subroutine matmulSqrtInvN_1map

!!$  ! Return map_out = invN * map
!!$  subroutine matmulInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_rms_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$    res%map = (self%siN%map)**2 * map%map
!!$  end subroutine matmulInvN_2map
!!$  
!!$  ! Return map_out = sqrtInvN * map
!!$  subroutine matmulSqrtInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_rms_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$    res%map = self%siN%map * map%map
!!$  end subroutine matmulSqrtInvN_2map

  ! Return RMS map
  subroutine returnRMS_rms_qucov(self, res, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: res
    integer(i4b),      intent(in),   optional  :: samp_group
    where (self%N_map%map(:,1:3) > 0.d0)
       res%map = sqrt(self%N_map%map(:,1:3))
    elsewhere
       res%map = 0
    end where
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) then
          where (self%samp_group_mask(samp_group)%p%map == 0.d0)
             res%map = 0
          end where
       end if
    end if
  end subroutine returnRMS_rms_qucov
  
  ! Return rms for single pixel
  function returnRMS_rms_qucov_pix(self, pix, pol, samp_group, ret_invN)
    implicit none
    class(comm_N_rms_QUcov),   intent(in)              :: self
    integer(i4b),        intent(in)              :: pix, pol
    real(dp)                                     :: returnRMS_rms_qucov_pix
    integer(i4b),        intent(in),   optional  :: samp_group
    logical(lgt),        intent(in),   optional  :: ret_invN
    if (self%N_map%map(pix,pol) > 0.d0) then
       returnRMS_rms_qucov_pix = sqrt(self%N_map%map(pix,pol))
    else
       returnRMS_rms_qucov_pix = 0
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) then
          if (self%samp_group_mask(samp_group)%p%map(pix,pol) == 0.d0) then
             returnRMS_rms_qucov_pix = 0
          end if
       end if
    end if
    if (present(ret_invN)) then
       if (ret_invN) returnRMS_rms_qucov_pix = 1d0/returnRMS_rms_qucov_pix**2
    end if
  end function returnRMS_rms_qucov_pix

  subroutine initialize_iN_siN(N, N_low, iN, siN, iN_low)
    implicit none
    class(comm_map),                   intent(in)    :: N
    class(comm_map),                   intent(in)    :: N_low
    real(dp),        dimension(1:,0:), intent(inout) :: iN
    real(dp),        dimension(1:,0:), intent(inout) :: siN
    real(dp),        dimension(1:,0:), intent(inout) :: iN_low

    integer(i4b) :: i
    real(dp) :: A(2,2)
    logical(lgt)        :: pol

    if (size(N%map, dim=2) .eq. 3) then
      pol = .true.
    else
      pol = .false.
    end if

    do i = 0, N%info%np-1
       ! T component
       if (N%map(i,1) > 0.) then
          iN(1,i)  = 1.d0 / N%map(i,1)
          siN(1,i) = sqrt(iN(1,i))
       else
          iN(1,i)  = 0.d0
          siN(1,i) = 0.d0
       end if

       if (pol) then
          ! QU block; check for positive definite matrix
          if (N%map(i,2)*N%map(i,3)-N%map(i,4)**2 > 0.) then
             A(1,1)  = N%map(i,2) ! QQ
             A(1,2)  = N%map(i,4) ! QU
             A(2,1)  = N%map(i,4) ! UQ
             A(2,2)  = N%map(i,3) ! UU

             siN(2,i) = A(1,1)    ! QQ
             siN(3,i) = A(2,2)    ! UU
             siN(4,i) = A(1,2)    ! QU = UQ
          else
             iN(2:4,i)  = 0.d0
             siN(2:4,i) = 0.d0
          end if
       end if
    end do

    do i = 0, N_low%info%np-1
       ! T component
       if (N_low%map(i,1) > 0.) then
          iN_low(1,i)  = 1.d0 / N_low%map(i,1)
       else
          iN_low(1,i)  = 0.d0
       end if

       if (pol) then
          ! QU block; check for positive definite matrix
          if (N_low%map(i,2)*N_low%map(i,3)-N_low%map(i,4)**2 > 0.) then
             A(1,1)  = N_low%map(i,2) ! QQ
             A(1,2)  = N_low%map(i,4) ! QU
             A(2,1)  = N_low%map(i,4) ! UQ
             A(2,2)  = N_low%map(i,3) ! UU

             iN_low(2,i)  = A(1,1)    ! QQ
             iN_low(3,i)  = A(2,2)    ! UU
             iN_low(4,i)  = A(1,2)    ! QU = UQ
          else
             iN_low(2:4,i)  = 0.d0
          end if
       end if
    end do

  end subroutine initialize_iN_siN

end module comm_N_rms_QUcov_mod
