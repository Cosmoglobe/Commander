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
module comm_F_mod
  use comm_comp_mod
  use comm_param_mod
  use comm_map_mod
  use comm_bp_mod
  use comm_F_int_mod
  use comm_F_int_0D_mod
  use comm_F_int_1D_mod
  use comm_F_int_2D_mod
  use comm_F_line_mod
  implicit none

!  private
!  public comm_F

  type :: comm_F
     ! Linked list variables
     class(comm_F), pointer :: nextLink => null()
     class(comm_F), pointer :: prevLink => null()

     ! Data variables
     class(comm_bp),    pointer :: bp
     class(comm_comp),  pointer :: c
     class(comm_map),   pointer :: F_diag
     class(comm_F_int), pointer :: F_int

     real(dp)                   :: checksum
   contains
     ! Linked list procedures
     procedure :: nextF    ! get the link after this link
     procedure :: prevF    ! get the link before this link
     procedure :: setNextF ! set the link after this link
     procedure :: addF     ! add new link at the end

     ! Data procedures
     procedure :: F          => matmulF
     procedure :: update     => updateTheta
     procedure :: write_F_to_FITS 
  end type comm_F

  interface comm_F
     procedure constructor_F
  end interface comm_F

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_F(info, comp, bp) result(c)
    implicit none
    class(comm_mapinfo), intent(in), target :: info
    class(comm_comp),    intent(in), target :: comp
    class(comm_bp),      intent(in), target :: bp
    class(comm_F),       pointer            :: c

    ! General parameters
    allocate(c)
    c%c        => comp
    c%bp       => bp
    c%F_diag   => comm_map(info)
    c%checksum =  -1.d30
    
  end function constructor_F

  ! Return map_out = F * map
  subroutine matmulF(self, map, res)
    implicit none
    class(comm_F),   intent(in)     :: self
    class(comm_map), intent(in)     :: map
    class(comm_map), intent(inout)  :: res
    res%map = self%F_diag%map * map%map
  end subroutine matmulF

  ! id = 1 => common {T,Q,U}
  ! id = 2 => T + common {Q,U}
  ! id = 3 => T + Q + U
  subroutine updateTheta(self, theta, id)
    implicit none
    class(comm_F),                 intent(inout) :: self
    class(comm_map), dimension(:), intent(in)    :: theta  ! Parameter maps; both alm and map
    integer(i4b),                  intent(in)    :: id

    integer(i4b) :: i, j, p
    real(dp)     :: checksum, t(self%c%npar)

    checksum = 0.d0
    do i = 1, self%c%npar
       checksum = checksum + sum(abs(theta(i)%alm))
    end do
    if (checksum == self%checksum) return

    ! Update with new parameter set
    self%checksum = checksum

    ! Temperature
    do p = 0, self%F_diag%info%np-1
       do j = 1, self%c%npar
          t(j) = theta(j)%map(p,1)
       end do
       self%F_diag%map(p,1) = self%F_int%eval(t)
    end do

    if (self%F_diag%info%nmaps == 3) then
       ! Stokes Q
       if (id > 1) then
          do p = 0, self%F_diag%info%np-1
             do j = 1, self%c%npar
                t(j) = theta(j)%map(p,2)
             end do
             self%F_diag%map(p,2) = self%F_int%eval(t)
          end do
       else
          self%F_diag%map(:,2) = self%F_diag%map(:,1)
       end if
       
       ! Stokes U
       if (id > 2) then
          do p = 0, self%F_diag%info%np-1
             do j = 1, self%c%npar
                t(j) = theta(j)%map(p,2)
             end do
             self%F_diag%map(p,3) = self%F_int%eval(t)
          end do
       else
          self%F_diag%map(:,3) = self%F_diag%map(:,2)
       end if
    end if

  end subroutine updateTheta
  
  subroutine write_F_to_FITS(self, filename)
    implicit none
    class(comm_F),    intent(in) :: self
    character(len=*), intent(in) :: filename
    call self%F_diag%writeFITS(filename)
  end subroutine write_F_to_FITS
  
  function nextF(self)
    class(comm_F) :: self
    class(comm_F), pointer :: nextF
    nextF => self%nextLink
  end function nextF

  function prevF(self)
    class(comm_F) :: self
    class(comm_F), pointer :: prevF
    prevF => self%prevLink
  end function prevF
  
  subroutine setNextF(self,next)
    class(comm_F) :: self
    class(comm_F), pointer :: next
    self%nextLink => next
  end subroutine setNextF

  subroutine addF(self,link)
    class(comm_F)         :: self
    class(comm_F), target :: link

    class(comm_F), pointer :: c
    
    c => self%nextLink
    do while (associated(c%nextLink))
       c => c%nextLink
    end do
    link%prevLink => c
    c%nextLink    => link
  end subroutine addF
  
end module comm_F_mod
