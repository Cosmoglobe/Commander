!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------

module utilities

  ! WARNING :most of the routines initially in the module utilities
  ! have been moved to module obsolete in version 1.2

  public :: die_alloc


contains

  !=======================================================================
  subroutine die_alloc(routine_name,array_name)
    !=======================================================================
    !     write a message and die
    !     to be used if status /= 0 after an allocate
    !=======================================================================
    implicit none
    character(len=*), intent(in) :: routine_name
    character(len=*), intent(in) :: array_name
    !=======================================================================

    print*,routine_name// &
         &     '> can not allocate memory for array : '//array_name
    stop 'program aborted'

    return
  end subroutine die_alloc


end module utilities
