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

module comm_tod_crosstalk_mod
  use comm_tod_mod
  use comm_status_mod
  use comm_tod_driver_mod
  implicit none
  private

  public comm_crosstalk

  type :: comm_crosstalk
    real(dp), dimension (:,:), allocatable      :: crosstalk_matrix
  contains
    procedure :: estimate_crosstalk_matrix
    procedure :: remove_crosstalk_signal
  end type comm_crosstalk   

  interface comm_crosstalk
    procedure xtalk_constructor
  end interface comm_crosstalk

  type crosstalk_pointer
    class(comm_crosstalk), pointer :: p => null()
  end type crosstalk_pointer

contains
 
  !initializes a comm_crosstalk class
  function xtalk_constructor(correlations) result(c)
    implicit none
    logical(lgt), dimension(:,:),    intent(in)     :: correlations
    class(comm_crosstalk), pointer                  :: c

    allocate(c)
    allocate(c%crosstalk_matrix(size(correlations(1,:)), size(correlations(:,1))))
    c%crosstalk_matrix = 1.d0

  end function xtalk_constructor
 
  ! estimates the crosstalk coeficients between each detector
  subroutine estimate_crosstalk_matrix(self,sd)
    ! The (ndet,ndet) crosstalk matrix is computed for each scan by a least square estimation
    ! y = bx + n
    !
    ! y := d_i
    ! x := [d_1, ..., d_{i-1}, d_{i+1}, ..., d_{ndet}]
    ! b = (x.T * x)^-1 * x.T * y
    !

    implicit none
    class(comm_crosstalk),                            intent(inout)      :: self
    class(comm_scandata),                             intent(in)         :: sd
    integer(i4b) :: i, j, k, l, n
    real(dp), dimension(:),   allocatable   :: xTy
    real(dp), dimension(:,:), allocatable   :: sub_data, inv_xTx

    allocate(sub_data(sd%ntod,sd%ndet-1))
    allocate(inv_xTx(sd%ndet-1,sd%ndet-1))
    allocate(xTy(sd%ndet-1))
    inv_xTx = 0.0
    xTy = 0.0

    do i=1, sd%ndet
       
       ! collect x
       k=0
       do j=1, sd%ndet
          if (i/=j) then
             k = k+1
             sub_data(:,k) = sd%tod(:,j)
          end if
       end do

       ! compute (x.T * x)^-1
       do j=1, sd%ndet-1
          do k=1, sd%ndet-1
             do n=1, sd%ntod
                inv_xTx(j,k) = inv_xTx(k,j) + sub_data(n,j) * sub_data(n,k)
             end do
          end do
       end do
       call invert_matrix(inv_xTx)

       ! compute (x.T * y)
       do j=1, sd%ndet-1
          do n=1, sd%ntod
             xTy(j) = xTy(j) + sub_data(n,j) * sd%tod(n,i)
          end do
       end do

       ! solve for b
       k=0
       do j=1, sd%ndet
          if (i/=j) then
             k = k+1
             do l=1, sd%ndet-1
                self%crosstalk_matrix(j,i) = self%crosstalk_matrix(j,i) + inv_xTx(l,k) * xTy(l)
             end do
          end if
       end do
    end do

    deallocate(xTy, sub_data, inv_xTx)

  end subroutine estimate_crosstalk_matrix

  subroutine remove_crosstalk_signal(self, sd)
    implicit none
    class(comm_crosstalk),                      intent(inout)   :: self
    class(comm_scandata),                       intent(inout)   :: sd

    integer(i4b) :: i, j
    real(dp), dimension(:,:), allocatable :: corr_tod

    allocate(corr_tod(sd%ntod,sd%ndet))
    corr_tod = 0.d0
    call invert_matrix(self%crosstalk_matrix)

    do i=1, sd%ndet
       do j=1, sd%ndet
          corr_tod(:,i) = corr_tod(:,i) + self%crosstalk_matrix(j,i) * sd%tod(:,j)
       end do
    end do
    sd%tod = corr_tod

    deallocate(corr_tod)

  end subroutine remove_crosstalk_signal

end module comm_tod_crosstalk_mod
