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
module comm_cr_utils
  use comm_utils
  implicit none

  integer(i4b)  :: ncr
  integer(i4b), allocatable, dimension(:,:) :: ind_comp    ! (start, length, nmaps)

  interface cr_insert_comp
     module procedure cr_insert_comp_1d, cr_insert_comp_2d
  end interface cr_insert_comp

  interface cr_extract_comp
     module procedure cr_extract_comp_1d, cr_extract_comp_2d
  end interface cr_extract_comp
  
contains

  subroutine cr_insert_comp_1d(id, add, x_in, x_out)
    implicit none

    integer(i4b),                  intent(in)    :: id
    logical(lgt),                  intent(in)    :: add
    real(dp),     dimension(1:),   intent(in)    :: x_in
    real(dp),     dimension(1:),   intent(inout) :: x_out

    integer(i4b) :: i, n, pos

    n   = size(x_in,1)
    pos = ind_comp(id,1)
    if (add) then
       x_out(pos:pos+n-1) = x_out(pos:pos+n-1) + x_in
    else
       x_out(pos:pos+n-1) = x_in
    end if
    
  end subroutine cr_insert_comp_1d

  subroutine cr_insert_comp_2d(id, add, x_in, x_out)
    implicit none

    integer(i4b),                   intent(in)    :: id
    logical(lgt),                   intent(in)    :: add
    real(dp),     dimension(0:,1:), intent(in)    :: x_in
    real(dp),     dimension(1:),    intent(inout) :: x_out

    integer(i4b) :: i, n, pos

    n   = size(x_in,1)
    pos = ind_comp(id,1)
    do i = 1, size(x_in,2)
       if (add) then
          x_out(pos:pos+n-1) = x_out(pos:pos+n-1) + x_in(:,i)
       else
          x_out(pos:pos+n-1) = x_in(:,i)
       end if
       pos                = pos+n
    end do
    
  end subroutine cr_insert_comp_2d


  subroutine cr_extract_comp_1d(id, x_in, x_out)
    implicit none

    integer(i4b),                            intent(in)  :: id
    real(dp),     dimension(:),              intent(in)  :: x_in
    real(dp),     dimension(:), allocatable, intent(out) :: x_out

    integer(i4b) :: i, n, pos, nmaps

    pos   = ind_comp(id,1)
    n     = ind_comp(id,2)
    allocate(x_out(n))
    x_out = x_in(pos:pos+n-1)
    
  end subroutine cr_extract_comp_1d

  subroutine cr_extract_comp_2d(id, x_in, x_out)
    implicit none

    integer(i4b),                              intent(in)  :: id
    real(dp),     dimension(:),                intent(in)  :: x_in
    real(dp),     dimension(:,:), allocatable, intent(out) :: x_out

    integer(i4b) :: i, n, pos, nmaps

    pos   = ind_comp(id,1)
    nmaps = ind_comp(id,3)
    n     = ind_comp(id,2)/nmaps
    allocate(x_out(0:n-1,nmaps))
    do i = 1, nmaps
       x_out(:,i) = x_in(pos:pos+n-1)
       pos        = pos + n
    end do
    
  end subroutine cr_extract_comp_2d

  subroutine cr_mask(id, mask_flag, mask)
    implicit none

    integer(i4b),                  intent(in)    :: id
    logical(lgt),                  intent(in)    :: mask_flag
    real(dp),     dimension(1:),   intent(inout) :: mask

    integer(i4b) :: i, n, pos, nmaps

    pos   = ind_comp(id,1)
    nmaps = ind_comp(id,3)
    n     = ind_comp(id,2)/nmaps
    do i = 1, nmaps
       if (mask_flag) then
          mask(pos:pos+n-1) = 1.d0
       else
          mask(pos:pos+n-1) = 0.d0
       end if
       pos        = pos + n
    end do
    
  end subroutine cr_mask

  
end module comm_cr_utils
