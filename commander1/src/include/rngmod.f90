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
!
!  This file contains the random number generator
!  used by the Planck LevelS package.
!  The generator is an F90 port of an xorshift generator described
!  in Marsaglia, Journal of Statistical Software 2003, vol 8.
!  It has a period of 2^128 - 1.
!
!  Copyright (C) 2003, 2004 Max-Planck-Society
!  Author: Martin Reinecke & Eric Hivon
!
module rngmod
  use healpix_types
  implicit none
  private

  public planck_rng, rand_init, rand_uni, rand_gauss

  type planck_rng
     integer(i4b) :: x,y,z,w
     real(dp)     :: gset
     logical(lgt) :: empty
  end type planck_rng

  real(dp) :: small = 1.0_dp/4294967296._dp

contains

  subroutine twiddle (v)
    ! users generally choose small numbers as seed
    ! this function will spread the seed (!) on a wider range
    integer(i4b), intent(inout) :: v
    integer i

    do i=1,9
      v=ieor(v,ishft(v,13))
      v=ieor(v,ishft(v,-17))
      v=ieor(v,ishft(v,5))
    end do
  end subroutine

  subroutine rand_init(handle, seed1, seed2, seed3, seed4)
    type(planck_rng), intent(out) :: handle
    integer(i4b), intent(in), optional :: seed1, seed2, seed3, seed4
    integer(i4b) :: i
    real(dp) :: junk

    handle%empty = .true.

    ! default values
    handle%x = 123456789
    handle%y = 362436069
    handle%z = 521288629
    handle%w = 88675123

    if (present(seed1)) then
      if (seed1/=0) handle%x = seed1
      !!!seed1 = abs(seed1)
    endif
    call twiddle(handle%x)
    if (present(seed2)) then
      if (seed2/=0) handle%y = seed2
    endif
    call twiddle(handle%y)
    if (present(seed3)) then
      if (seed3/=0) handle%z = seed3
    endif
    call twiddle(handle%z)
    if (present(seed4)) then
      if (seed4/=0) handle%w = seed4
    endif
    call twiddle(handle%w)

    ! do some burn in so that user provided seed(s)
    ! interact with default ones
    do i=1, 16
       junk = rand_uni(handle)
    enddo

  end subroutine

  function rand_uni(handle)
    type(planck_rng), intent(inout) :: handle
    real(dp) :: rand_uni

    integer(i4b) tmp
    tmp = ieor(handle%x,ishft(handle%x,11))
    handle%x = handle%y
    handle%y = handle%z
    handle%z = handle%w
    handle%w = ieor(ieor(handle%w,ishft(handle%w,-19)),ieor(tmp,ishft(tmp,-8)))
    rand_uni = small*handle%w
    if (rand_uni<0.0_dp) rand_uni=rand_uni+1.0_dp
  end function rand_uni

  function rand_gauss(handle)
    type(planck_rng), intent(inout) :: handle
    real(dp) :: rand_gauss
    real(dp) :: fac,rsq,v1,v2

    if (handle%empty) then
1      v1 = 2.0_dp * rand_uni(handle) - 1.0_dp
       v2 = 2.0_dp * rand_uni(handle) - 1.0_dp
       rsq = v1**2 + v2**2
       if(rsq>=1.0_dp .or. rsq==0.0_dp) goto 1
       fac = sqrt(-2.0_dp * log(rsq)/rsq)
       handle%gset = v1*fac
       rand_gauss  = v2*fac
       handle%empty= .false.
    else
       rand_gauss   = handle%gset
       handle%empty = .true.
    endif
  end function rand_gauss

end module rngmod
