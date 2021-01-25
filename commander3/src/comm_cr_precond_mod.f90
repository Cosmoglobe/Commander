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
module comm_cr_precond_mod
  use comm_utils
  implicit none

  type invM
     integer(i4b)                              :: n
     integer(i4b), allocatable, dimension(:)   :: ind, comp2ind
     real(dp),     allocatable, dimension(:,:) :: M0, M
  end type invM

  type precond
     type(invM), allocatable, dimension(:,:)      :: invM_diff ! (0:nalm-1,nmaps)
     type(invM), allocatable, dimension(:,:)      :: invM_src  ! (1,nmaps)
     type(invM), allocatable, dimension(:,:)      :: invM_temp ! (1,1)
  end type precond

  type(precond) :: P_cr
  
end module comm_cr_precond_mod
