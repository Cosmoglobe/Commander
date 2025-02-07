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
module comm_N_QUcov_mod
  use comm_N_mod
  implicit none

  private
  public comm_N_QUcov, comm_N_QUcov_ptr
  
  type, extends (comm_N) :: comm_N_QUcov
     real(dp),        allocatable, dimension(:,:) :: Ncov
     real(dp),        allocatable, dimension(:,:) :: siN
     real(dp),        allocatable, dimension(:,:) :: iN
     class(comm_map), pointer                     :: siN_diag => null()
   contains
     ! Data procedures
     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS_QUcov
     procedure :: rms_pix     => returnRMS_QUcov_pix
     procedure :: update_N    => update_N_QUcov
  end type comm_N_QUcov

  interface comm_N_QUcov
     procedure constructor
  end interface comm_N_QUcov

  type comm_N_QUcov_ptr
     type(comm_N_QUcov), pointer :: p => null()
  end type comm_N_QUcov_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask)
    implicit none
    class(comm_N_QUcov),                pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask

    
    ! General parameters
    allocate(constructor)


  end function constructor

  subroutine update_N_QUcov(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_QUcov),                 intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map


  end subroutine update_N_QUcov



  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: map
    integer(i4b),        intent(in),   optional  :: samp_group
    
    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m


  end subroutine matmulInvN_1map


  ! Return map_out = N * map
  subroutine matmulN_1map(self, map, samp_group)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: map
    integer(i4b),        intent(in),   optional  :: samp_group

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m


  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: map
    integer(i4b),        intent(in),   optional  :: samp_group

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:) :: m, invN_m


  end subroutine matmulSqrtInvN_1map


  ! Return RMS map
  subroutine returnRMS_QUcov(self, res, samp_group)
    implicit none
    class(comm_N_QUcov), intent(in)              :: self
    class(comm_map),     intent(inout)           :: res
    integer(i4b),        intent(in),   optional  :: samp_group
    where (self%siN_diag%map > 0.d0)
       res%map = 1.d0/self%siN_diag%map
    elsewhere
       res%map = 0
    end where
  end subroutine returnRMS_QUcov
  
  ! Return rms for single pixel
  function returnRMS_QUcov_pix(self, pix, pol, samp_group, ret_invN)
    implicit none
    class(comm_N_QUcov),   intent(in)            :: self
    integer(i4b),          intent(in)            :: pix, pol
    real(dp)                                     :: returnRMS_QUcov_pix
    integer(i4b),        intent(in),   optional  :: samp_group
    logical(lgt),        intent(in),   optional  :: ret_invN

    if (self%siN_diag%map(pix,pol) > 0.d0) then
       returnRMS_QUcov_pix = 1.d0/self%siN_diag%map(pix,pol)
    else
       returnRMS_QUcov_pix = 0
    end if
    if (present(ret_invN)) then
       if (ret_invN) returnRMS_QUcov_pix = self%siN_diag%map(pix,pol)**2
    end if
  end function returnRMS_QUcov_pix

end module comm_N_QUcov_mod
