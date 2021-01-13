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
module udgrade_nr

  ! --- in udgrade_template.f90 ---
  ! subroutine sub_udgrade_nest
  ! subroutine udgrade_ring
  ! subroutine udgrade_nest

  use healpix_types
  use misc_utils

  IMPLICIT none

  private

  public :: sub_udgrade_nest, udgrade_ring, udgrade_nest

  ! SP and DP interface
  interface sub_udgrade_nest
     module procedure sub_udgrade_nest_s, sub_udgrade_nest_d
  end interface

  ! SP and DP, 1-dim and N-dim interface
  interface udgrade_ring
     module procedure udgrade_ring_1d_s, udgrade_ring_1d_d, udgrade_ring_nd_s, udgrade_ring_nd_d
  end interface

  interface udgrade_nest
     module procedure udgrade_nest_1d_s, udgrade_nest_1d_d, udgrade_nest_nd_s, udgrade_nest_nd_d
  end interface
  !---------------------

contains

  ! include SP implementation of routines
  include 'udgrade_s_inc.f90'
  
  ! include DP implementation of routines
  include 'udgrade_d_inc.f90'


end module udgrade_nr
