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
!! Generator for uniformly distributed and Gaussian random numbers
module ran_tools

use healpix_types
implicit none
private

public randgauss_boxmuller, ran_mwc

contains

!! Box-Muller method for converting uniform into Gaussian deviates.
function randgauss_boxmuller(iseed)
  integer(i4b), intent(inout) :: iseed  !! random number state
  real(sp) :: randgauss_boxmuller       !! result
  logical, save :: empty=.true.
  real(sp) :: fac,rsq,v1,v2
  real(sp),save :: gset !! test

  if (empty .or. iseed < 0) then ! bug correction, EH, March 13, 2003
    do
      v1 = 2._sp*ran_mwc(iseed) - 1._sp
      v2 = 2._sp*ran_mwc(iseed) - 1._sp
      rsq = v1**2+v2**2
      if((rsq<1._sp) .and. (rsq>0._sp)) exit
    end do

    fac = sqrt(-2._sp*log(rsq)/rsq)
    gset = v1*fac
    randgauss_boxmuller = v2*fac
    empty = .false.
  else
    randgauss_boxmuller = gset
    empty = .true.
  endif
end function randgauss_boxmuller


!! This routine implements the Marsaglia "multiply with carry method"
!! and adds a Marsaglia shift sequence to it.
!! (cf. references at http://www.cs.hku.hk/)
!!
!! You are welcome to replace this with your preferred generator
!! of uniform variates in the interval ]0,1[ (i.e. excluding 0 and 1).
!!
!! (Re-)initialise with setting iseed to a negative number.
!! Note that iseed=0 gives the same sequence as iseed=-1.
!! After initialisation iseed becomes abs(iseed) (or 1 if it was 0).

!! @author B. D. Wandelt, May 1999

function ran_mwc(iseed)
  integer(i4b), intent(inout):: iseed !! random number state
  real(sp) :: ran_mwc                 !! result

  integer(i4b) :: i,iseedl,iseedu,mwc,combined
  integer(i4b),save :: upper,lower,shifter
  integer(i4b),parameter :: mask16=65535,mask30=2147483647
  real(sp),save :: small
  logical, save :: first=.true.

  if (first .or. (iseed<=0)) then
    if (iseed==0) iseed=-1
    iseed = abs(iseed)
    small = nearest (1._sp,-1._sp)/mask30

    ! Users often enter small seeds - I spread them out using the
    ! Marsaglia shifter a few times.
    shifter=iseed
    do i=1,9
      shifter=ieor(shifter,ishft(shifter,13))
      shifter=ieor(shifter,ishft(shifter,-17))
      shifter=ieor(shifter,ishft(shifter,5))
    enddo

    iseedu=ishft(shifter,-16)
    upper=ishft(iseedu+8765,16)+iseedu !This avoids the fixed points.
    iseedl=iand(shifter,mask16)
    lower=ishft(iseedl+4321,16)+iseedl !This avoids the fixed points.

    first=.false.
  endif

  do
    shifter=ieor(shifter,ishft(shifter,13))
    shifter=ieor(shifter,ishft(shifter,-17))
    shifter=ieor(shifter,ishft(shifter,5))

    upper=36969*iand(upper,mask16)+ishft(upper,-16)
    lower=18000*iand(lower,mask16)+ishft(lower,-16)

    mwc=ishft(upper,16)+iand(lower,mask16)

    combined=iand(mwc,mask30)+iand(shifter,mask30)

    ran_mwc=small*iand(combined,mask30)
    if(ran_mwc/=0._sp) exit
  end do
end function ran_mwc

end module ran_tools
