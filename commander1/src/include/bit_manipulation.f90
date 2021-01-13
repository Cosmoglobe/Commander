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
!! This module provides functions that permute bits of the binary
!! representation of a number.
!!
!! @author Benjamin D. Wandelt October 1997
!! edited by E. Hivon, October 2001 to be 'F' compatible
module bit_manipulation

use healpix_types
implicit none
private

integer(i4b), parameter :: oddbits=89478485,evenbits=178956970

public :: swapLSBMSB, invswapLSBMSB, invLSB, invMSB

contains

!! Returns i with even and odd bit positions interchanged.
function swapLSBMSB(i)
  integer(i4b) :: swapLSBMSB
  integer(i4b), intent(in) :: i

  swapLSBMSB = IAND(i,evenbits)/2 + IAND(i,oddbits)*2
end function swapLSBMSB

!! Returns NOT(i) with even and odd bit positions interchanged.
function invswapLSBMSB(i)
  integer(i4b) :: invswapLSBMSB
  integer(i4b), intent(in) :: i

  invswapLSBMSB = NOT(swapLSBMSB(i))
end function invswapLSBMSB

!! Returns i with odd (1,3,5,...) bits inverted.
function invLSB(i)
  integer(i4b) :: invLSB
  integer(i4b), intent(in) :: i

  invLSB = IEOR(i,oddbits)
end function invLSB

!! Returns i with even (0,2,4,...) bits inverted.
function invMSB(i)
  integer(i4b) :: invMSB
  integer(i4b), intent(in) :: i

  invMSB = IEOR(i,evenbits)
end function invMSB

end module bit_manipulation
