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
module statistics
  !---------------------------------
  ! subroutine compute_statistics
  ! subroutine print_statistics
  ! function   median
  ! type       tstats
  !
  ! EH, IPAC, 2005-04-25
  !---------------------------------
  use healpix_types
  use misc_utils, only: assert
  use m_indmed, only: indmed
  implicit none

  private
  public :: median
  public :: tstats
  public :: compute_statistics, print_statistics

  interface median
     module procedure median_s, median_d
  end interface

  interface compute_statistics
     module procedure comp_stats_s, comp_stats_d
  end interface

  type tstats
     integer(i4b) :: ntot, nvalid
     real(dp)     :: mind, maxd
     real(dp)     :: average, absdev
     real(dp)     :: rms, var
     real(dp)     :: skew, kurt
  end type tstats


contains
  !
  ! compute_statistics(data, stats [, badval])
  !
  !===================================================================
  subroutine comp_stats_s(data, stats, badval)
  !===================================================================
    integer(i4b),     parameter :: KMAP = SP
    real(KMAP),       parameter :: ONE  = 1.0_KMAP
    character(len=*), parameter :: code = 'compute_statistics'
    !
    real(KMAP), dimension(0:), intent(in)  :: data
    type(tstats),              intent(out) :: stats
    real(KMAP),     optional,  intent(in)  :: badval
    !
    real(KMAP)   :: sentinel, precis
    integer(i4b) :: nvalid, i, n
    real(dp)     :: mind, maxd
    real(dp)     :: average, absdev
    real(dp)     :: rms, var
    real(dp)     :: skew, kurt
    real(dp)     :: x, dx, eps, p
    !-----------------------------------------------------------------
    precis = 10*epsilon(ONE)
    sentinel = -huge(ONE)
    if (present(badval)) then
       sentinel = badval
       call assert(sentinel /= 0, code//': BadValue should not be set to 0.0')
    endif

    n       = size(data)
    nvalid  = 0
    mind    =   huge(ONE)
    maxd    = - huge(ONE)
    average = 0.0_dp
    var     = 0.0_dp
    rms     = 0.0_dp
    skew    = 0.0_dp
    kurt    = 0.0_dp

    do i=0, n-1
       x = data(i)
       if (abs(x/sentinel-ONE) > precis) then
          mind = min(mind, x)
          maxd = max(maxd, x)
          average = average + x
          nvalid = nvalid + 1
       endif
    enddo

    if (nvalid == 0) goto 100
    average = average / nvalid

    do i=0, n-1
       x = data(i)
       if (abs(x/sentinel-ONE) > precis) then
          dx = x - average
          eps = eps + dx
          absdev = absdev + abs(dx)
          p = dx * dx
          var = var + p
          p = p * dx
          skew = skew + p
          p = p * dx
          kurt = kurt + p
       endif
    enddo

100 continue
    if (nvalid > 0) then
       absdev = absdev / nvalid
    else
       print*,'=================================='
       print*,'No valid data point for statistics'
       print*,'=================================='
    endif

    if (nvalid > 1) then
       var    = (var - eps * eps / nvalid) / (nvalid - 1)
       rms = sqrt(var)
    else
       print*,'============================================'
       print*,'Needs at least 2 valid points for statistics'
       print*,'============================================'
    endif

    if (var /= 0.0_DP) then
       skew = skew / (nvalid * rms**3)
       kurt = kurt / (nvalid * var**2) - 3.0_dp
    else
       print*,'=========================================='
       print*,'No skewness or kurtosis when zero variance'
       print*,'=========================================='
    endif

    stats%ntot    = n 
    stats%nvalid  = nvalid  
    stats%mind    = mind    
    stats%maxd    = maxd    
    stats%average = average 
    stats%absdev  = absdev
    stats%var     = var     
    stats%rms     = rms     
    stats%skew    = skew    
    stats%kurt    = kurt    

    return
  end subroutine comp_stats_s
  !===================================================================
  subroutine comp_stats_d(data, stats, badval)
  !===================================================================
    integer(i4b),     parameter :: KMAP = DP
    real(KMAP),       parameter :: ONE  = 1.0_KMAP
    character(len=*), parameter :: code = 'compute_statistics'
    !
    real(KMAP), dimension(0:), intent(in)  :: data
    type(tstats),              intent(out) :: stats
    real(KMAP),     optional,  intent(in)  :: badval
    !
    real(KMAP)   :: sentinel, precis
    integer(i4b) :: nvalid, i, n
    real(dp)     :: mind, maxd
    real(dp)     :: average, absdev
    real(dp)     :: rms, var
    real(dp)     :: skew, kurt
    real(dp)     :: x, dx, eps, p
    !-----------------------------------------------------------------
    precis = 10*epsilon(ONE)
    sentinel = -huge(ONE)
    if (present(badval)) then
       sentinel = badval
       call assert(sentinel /= 0, code//': BadValue should not be set to 0.0')
    endif

    n       = size(data)
    nvalid  = 0
    mind    =   huge(ONE)
    maxd    = - huge(ONE)
    average = 0.0_dp
    var     = 0.0_dp
    rms     = 0.0_dp
    skew    = 0.0_dp
    kurt    = 0.0_dp

    do i=0, n-1
       x = data(i)
       if (abs(x/sentinel-ONE) > precis) then
          mind = min(mind, x)
          maxd = max(maxd, x)
          average = average + x
          nvalid = nvalid + 1
       endif
    enddo

    if (nvalid == 0) goto 100
    average = average / nvalid

    do i=0, n-1
       x = data(i)
       if (abs(x/sentinel-ONE) > precis) then
          dx = x - average
          eps = eps + dx
          absdev = absdev + abs(dx)
          p = dx * dx
          var = var + p
          p = p * dx
          skew = skew + p
          p = p * dx
          kurt = kurt + p
       endif
    enddo

100 continue
    if (nvalid > 0) then
       absdev = absdev / nvalid
    else
       print*,'=================================='
       print*,'No valid data point for statistics'
       print*,'=================================='
    endif

    if (nvalid > 1) then
       var    = (var - eps * eps / nvalid) / (nvalid - 1)
       rms = sqrt(var)
    else
       print*,'============================================'
       print*,'Needs at least 2 valid points for statistics'
       print*,'============================================'
    endif

    if (var /= 0.0_DP) then
       skew = skew / (nvalid * rms**3)
       kurt = kurt / (nvalid * var**2) - 3.0_dp
    else
       print*,'=========================================='
       print*,'No skewness or kurtosis when zero variance'
       print*,'=========================================='
    endif

    stats%ntot    = n 
    stats%nvalid  = nvalid  
    stats%mind    = mind    
    stats%maxd    = maxd    
    stats%average = average 
    stats%absdev  = absdev
    stats%var     = var     
    stats%rms     = rms     
    stats%skew    = skew    
    stats%kurt    = kurt    

    return
  end subroutine comp_stats_d


  !--------------------------------------------------

  subroutine print_statistics(stats)

    type(tstats),        intent(in) :: stats
    integer(i4b) :: nmiss

    nmiss = stats%ntot-stats%nvalid
    print*,'Pixels  = ', stats%nvalid,' / ',stats%ntot
    write(*,9000) 'Missing = ',nmiss,(nmiss*100.)/stats%ntot
    write(*,9010) 'Average = ',stats%average
    write(*,9010) 'Abs dev = ',stats%absdev
    write(*,9010) 'Rms     = ',stats%rms
    write(*,9010) 'Min     = ',stats%mind
    write(*,9010) 'Max     = ',stats%maxd
    write(*,9010) 'Variance= ',stats%var
    write(*,9010) 'Skewness= ',stats%skew
    write(*,9010) 'Kurtosis= ',stats%kurt
    print*,'   '

9000 format(1x,a,i12,', (',f8.4,'  %)')
9010 format(1x,a,1pe14.6)
  end subroutine print_statistics

  !======================================================
  ! MEDIAN
  ! med = median(data, [badval, even]) 
  !======================================================

  function median_s(data, badval, even) result (med)
    integer(I4B), parameter :: KMAP = SP
    real(KMAP),   parameter :: ONE  = 1.0_KMAP
    character(len=*), parameter :: code = 'median'
    ! dummy variables
    real(KMAP)                                      :: med
    real(KMAP), dimension(:),  intent(in), target   :: data
    real(KMAP),                intent(in), optional :: badval
    logical(LGT),              intent(in), optional :: even
    ! local variables
    real(KMAP), dimension(:), pointer :: gdata
    real(KMAP)                        :: precis
    logical(LGT)                      :: do_even, do_bad
    integer(I4B)                      :: ndata, ngood, j, k
    !---------------------------------------------------------------------------
    precis = 10*epsilon(ONE)

    do_bad = present(badval)
    if (do_bad) call assert(badval /= 0, code//': BadValue should not be set to 0.0')
    do_even = .false.
    if (present(even)) do_even = even
    
    ! select valid data
    ndata = size(data)
    if (do_bad) then
       ngood = count(abs(data/badval-ONE) > precis)
       allocate(gdata(1:ngood))
       gdata = pack(data, mask= (abs(data/badval-ONE) > precis))
    else
       ngood = ndata
       gdata => data
    endif
    
    ! find median value
    if (do_even .and. mod(ngood,2) == 0) then
       call indmed( gdata, j)
       call indmed(-gdata, k)
       med = 0.5_KMAP * (gdata(j) + gdata(k))
    else
       call indmed(gdata, j)
       med = gdata(j)
    endif

    return
  end function median_s
  !=================================================================
  function median_d(data, badval, even) result (med)
    integer(I4B), parameter :: KMAP = DP
    real(KMAP),   parameter :: ONE  = 1.0_KMAP
    character(len=*), parameter :: code = 'median'
    ! dummy variables
    real(KMAP)                                      :: med
    real(KMAP), dimension(:),  intent(in), target   :: data
    real(KMAP),                intent(in), optional :: badval
    logical(LGT),              intent(in), optional :: even
    ! local variables
    real(KMAP), dimension(:), pointer :: gdata
    real(KMAP)                        :: precis
    logical(LGT)                      :: do_even, do_bad
    integer(I4B)                      :: ndata, ngood, j, k
    !---------------------------------------------------------------------------
    precis = 10*epsilon(ONE)

    do_bad = present(badval)
    if (do_bad) call assert(badval /= 0, code//': BadValue should not be set to 0.0')
    do_even = .false.
    if (present(even)) do_even = even
    
    ! select valid data
    ndata = size(data)
    if (do_bad) then
       ngood = count(abs(data/badval-ONE) > precis)
       allocate(gdata(1:ngood))
       gdata = pack(data, mask= (abs(data/badval-ONE) > precis))
    else
       ngood = ndata
       gdata => data
    endif
    
    ! find median value
    if (do_even .and. mod(ngood,2) == 0) then
       call indmed( gdata, j)
       call indmed(-gdata, k)
       med = 0.5_KMAP * (gdata(j) + gdata(k))
    else
       call indmed(gdata, j)
       med = gdata(j)
    endif

    return
  end function median_d
  !=================================================================

end module statistics

