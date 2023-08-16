!===============================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
!===============================================================================

! This module handles cosmic rays correction to time ordered data

module comm_tod_cr_mod
  use comm_tod_mod
  use comm_utils
  implicit none

  private
  public comm_cr

  type :: comm_cr
    real(dp) :: amplitude, tau ! parameters to fit a template

  contains 
    procedure :: find_cosmic_rays
    procedure :: fit_template
    procedure :: cr_subtraction

  end type comm_cr

interface

  module subroutine find_cosmic_rays(self)
    implicit none
    class(comm_cr), intent(inout) :: self
  end subroutine find_cosmic_rays

  module subroutine fit_template(self)
    implicit none
    class(comm_cr), intent(inout) :: self
  end subroutine fit_template

  module subroutine cr_subtraction(self)
    implicit none
    class(comm_cr), intent(inout) :: self
  end subroutine cr_subtraction


end interface

end module comm_tod_cr_mod



