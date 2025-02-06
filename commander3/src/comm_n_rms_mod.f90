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
module comm_N_rms_mod
  use comm_N_mod
  implicit none

  private
  public comm_N_rms, comm_N_rms_ptr
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: siN        => null()
     class(comm_map), pointer :: siN_lowres => null()
     class(comm_map), pointer :: rms0       => null()
   contains
     ! Data procedures

     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map_lowres
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS_rms
     procedure :: rms_pix     => returnRMS_rms_pix
     procedure :: update_N    => update_N_rms
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

  type comm_N_rms_ptr
     type(comm_N_rms), pointer :: p => null()
  end type comm_N_rms_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask, map)
    implicit none
    class(comm_N_rms),                  pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask
    class(comm_map),                    pointer, intent(in), optional :: map

    integer(i4b)       :: i, ierr, tmp, nside_smooth
    real(dp)           :: sum_noise, npix
    type(comm_mapinfo), pointer :: info_smooth => null()

    ! General parameters
    allocate(constructor)

    call constructor%update_N(info, handle, mask, regnoise, &
         & noisefile=trim(cpar%ds_noisefile(id_abs)))


  end function constructor

  subroutine update_N_rms(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_rms),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, j, ierr
    real(dp)     :: sum_tau, sum_tau2, sum_noise, npix
    class(comm_map),     pointer :: invW_tau => null(), iN => null()
    class(comm_mapinfo), pointer :: info_lowres => null()

    self%rms0     => comm_map(info, noisefile)
    self%siN     => comm_map(self%rms0)

    ! Set up diagonal covariance matrix
    self%invN_diag => comm_map(info)
    self%invN_diag%map = self%siN%map**2
    call compute_invN_lm(self%invN_diag)


  end subroutine update_N_rms

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    integer(i4b)  :: nmaps_band, nmaps_inp
  end subroutine matmulInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_1map_lowres(self, map, samp_group)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    integer(i4b)  :: nmaps_band, nmaps_inp
  end subroutine matmulInvN_1map_lowres

  ! Return map_out = N * map
  subroutine matmulN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    integer(i4b)  :: nmaps_band, nmaps_inp
  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    integer(i4b)  :: nmaps_band, nmaps_inp, nmaps
  end subroutine matmulSqrtInvN_1map


  ! Return RMS map
  subroutine returnRMS_rms(self, res, samp_group)
    implicit none
    class(comm_N_rms), intent(in)              :: self
    class(comm_map),   intent(inout)           :: res
    integer(i4b),      intent(in),   optional  :: samp_group
  end subroutine returnRMS_rms
  
  ! Return rms for single pixel
  function returnRMS_rms_pix(self, pix, pol, samp_group, ret_invN)
    implicit none
    class(comm_N_rms),   intent(in)              :: self
    integer(i4b),        intent(in)              :: pix, pol
    real(dp)                                     :: returnRMS_rms_pix
    integer(i4b),        intent(in),   optional  :: samp_group
    logical(lgt),        intent(in),   optional  :: ret_invN

    returnRMS_rms_pix = 0d0
  end function returnRMS_rms_pix

end module comm_N_rms_mod
