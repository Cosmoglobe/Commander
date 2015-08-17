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
module coord_v_convert

! module provided by C.Burigana & D.Maino (U. Milano)
! edited by E.Hivon

  use healpix_types
  use misc_utils, only : strupcase

  Implicit NONE

  private
  real(DP), parameter, private :: DTOR = PI/180.0_dp

  public :: coordsys2euler_zyz ! produce Euler angles (ZYZ) for coordinate conversion
  public :: xcc_v_convert ! front end routine
  public :: xcc_DP_E_TO_Q ! Ecl -> Eq
  public :: xcc_DP_E_TO_G ! Ecl -> Gal
  public :: xcc_DP_G_TO_E ! Gal -> Ecl
  public :: xcc_DP_Q_TO_E ! Eq  -> Ecl
  public :: XCC_DP_PRECESS ! precess ecliptic coord

contains
  subroutine coordsys2euler_zyz(iepoch, oepoch, isys, osys, psi, theta, phi)
    !
    ! routine producing Euler angles psi, theta, phi 
    !  (for rotation around the unrotated Z, Y, Z axes respectively)
    ! corresponding to changes of coordinates from ISYS to OSYS
    ! for respective epochs IEPOCH and OEPOCH
    ! the angles definition match the one required for rotate_alm.
    ! EH, March 2005
    !     April 06 2005, corrected psi calculation
    !
    use pix_tools, only: angdist
    real(dp),         intent(in)  :: iepoch, oepoch
    character(len=*), intent(in)  :: isys, osys
    real(dp),         intent(out) :: psi, theta, phi

    real(dp), dimension(1:3) :: v1, v2, v3, v1p, v2p, v3p
    !-----------------------
 
    ! generate basis vector
    v1 = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    v2 = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    v3 = (/ 0.0_dp, 0.0_dp, 1.0_dp /)

    ! rotate them
    call xcc_V_CONVERT(v1,iepoch,oepoch,isys,osys,v1p)
    call xcc_V_CONVERT(v2,iepoch,oepoch,isys,osys,v2p)
    call xcc_V_CONVERT(v3,iepoch,oepoch,isys,osys,v3p)

    ! normalize vectors
    v1p = v1p / sqrt(dot_product(v1p,v1p))
    v2p = v2p / sqrt(dot_product(v2p,v2p))
    v3p = v3p / sqrt(dot_product(v3p,v3p))

    ! figure out Euler angles
    call angdist(v3,v3p, theta)
    psi = atan2(v2p(3),-v1p(3))
    phi = atan2(v3p(2), v3p(1))

    return
  end subroutine coordsys2euler_zyz

!
  Subroutine xcc_V_CONVERT(ivector,iepoch,oepoch,isys,osys,ovector)

!
!    Function to convert between standard coordinate systems, including
!  precession.
!  
!
!      Q,C:  Equatorial coordinates
!      E  :  Ecliptic coordinates
!      G  :  Galactic coordinates
!      SG :  Super Galactic coordinates  (TB implemented)
!      U  :  Universal Rest-frame coordinates (TB implemented)
!
!  03/14/88, a.c.raugh, STX.
!  11/15/01, C.Burigana & D.Maino, OAT
!
!
!  Arguments:
!
    real(dp) ::      IVECTOR(1:3)   ! Input coordinate vector
    real(dp) ::      IEPOCH         ! Input epoch
    real(dp) ::      OEPOCH         ! Output epoch        
    real(dp) ::      OVECTOR(1:3)   ! Output coordinate vector
    Character(len=*) ::    ISYS           ! Input system
    Character(len=*) ::    OSYS           ! Output system
!
!  Program variables:
!                 
    real(dp) ::      XVECTOR(1:3)     ! Temporary coordinate holder
    real(dp) ::      YVECTOR(1:3)     ! Temporary coordinate holder
    Character(len=20) ::    isysu     ! Input system, uppercase
    Character(len=20) ::    osysu     ! Output system, uppercase
!

    isysu = trim(strupcase(isys))
    osysu = trim(strupcase(osys))

!  If the input coordinate system was not ecliptic, convert to ecliptic
!  as intermediate form:
!
    if (trim(isysu) /= 'E') then
       if (trim(isysu) == 'Q' .or. trim(isysu) == 'C') then
          call xcc_dp_q_to_e(ivector,iepoch,xvector)
       elseif (trim(isysu) == 'G') then
          call xcc_dp_g_to_e(ivector,iepoch,xvector)
!            elseif (isys == 'SG') then
!              to be implemented
!            elseif (isys == 'U') then
!               to be implemented
       else
          print*,' unknown input coordinate system: '//trim(isysu)
       endif
    else
       xvector(1:3)=ivector(1:3)
    endif
!
!  If the begin and end epoch are different, precess:
!
    if (iepoch /= oepoch) then
       call xcc_dp_precess(xvector,iepoch,oepoch,yvector)
    else
       yvector(1:3)=xvector(1:3)
    endif
!
!  If output system is not ecliptic, convert from ecliptic to desired
!  system:
!
    if (trim(osysu) /= 'E') then 
       if (trim(osysu) == 'Q' .or. trim(osysu) == 'C') then
          call xcc_dp_e_to_q(yvector,oepoch,ovector)
       elseif (trim(osysu) == 'G') then
          call xcc_dp_e_to_g(yvector,oepoch,ovector)
       else
          print*,' unknown output coordinate system: '//trim(osysu)          
       endif
    else
       ovector(1:3)=yvector(1:3)
    endif
    return
  end subroutine xcc_v_convert
!===============================================================

  Subroutine xcc_DP_E_TO_Q(ivector,epoch,ovector)
!
!    Routine to convert from ecliptic coordinates to equatorial (celestial)
!  coordinates at the given epoch.  Adapted from the COBLIB routine by Dr.
!  E. Wright.
!
!  02/01/88, a.c.raugh, STX.
!  11/16/01, C.Burigana,D.Maino, OAT
!
!
!  Arguments:
!
    real(dp) ::     IVECTOR(3)     ! Input coordinate vector
    real(dp) ::     EPOCH          ! Equatorial coordinate epoch
    real(dp) ::     OVECTOR(3)     ! Output coordinate vector
!
!  Program variables:
!
    real(dp) ::     T              ! Julian centuries since 1900.0
    real(dp) ::     EPSILON        ! Obliquity of the ecliptic
    real(dp) ::     HVECTOR(3)     ! Temporary holding place for output vector
    real(dp) :: dc, ds
!
!  Set-up:
!
    T = (epoch - 1900.d0) / 100.d0 
    epsilon = 23.452294d0 - 0.0130125d0*T - 1.63889d-6*T**2 &
         &  + 5.02778d-7*T**3
!
!  Conversion:
!
    dc = cos(DTOR * epsilon)
    ds = sin(DTOR * epsilon)
    hvector(1) = ivector(1)
    hvector(2) = dc*ivector(2) - ds*ivector(3)
    hvector(3) = dc*ivector(3) + ds*ivector(2)
!
!  Move to output vector:
!
    ovector(1:3) = hvector(1:3)
!
    return
  end Subroutine xcc_DP_E_TO_Q
!===============================================================

  Subroutine xcc_DP_E_TO_G(ivector,epoch,ovector)
!
!    Routine to convert ecliptic (celestial) coordinates at the given
!  epoch to galactic coordinates.  The ecliptic coordinates are first
!  precessed to 2000.0, then converted. 
!
!                          --------------
!
!  Galactic Coordinate System Definition (from Zombeck, "Handbook of
!  Space Astronomy and Astrophysics", page 71):
!
!         North Pole:  12:49       hours right ascension (1950.0)
!                     +27.4        degrees declination   (1950.0)
!
!    Reference Point:  17:42.6     hours right ascension (1950.0)
!                     -28 55'      degrees declination   (1950.0)
!
!                          --------------
!
!  03/08/88, a.c.raugh, STX.
!  11/16/01, C.Burigana, D.Maino, OAT
!
!
!  Arguments:
!
    real(dp) ::    IVECTOR(3)      ! Input coordinate vector
    real(dp) ::    EPOCH           ! Input coordinate epoch
    real(dp) ::    OVECTOR(3)      ! Output coordinate vector
!
!  Program variables:
!
    real(dp) ::    HVECTOR(3)      ! Temporary holding place for coordinates
    real(dp) ::    T(3,3)          ! Galactic coordinate transformation matrix
    integer(i4b) ::  I             ! Loop variables
! 
    Data T / -0.054882486d0, -0.993821033d0, -0.096476249d0, & ! 1st column
         &    0.494116468d0, -0.110993846d0,  0.862281440d0, & ! 2nd column
         &   -0.867661702d0, -0.000346354d0,  0.497154957d0/ ! 3rd column
!
!     Precess coordinates to 2000.0 if necesssary:
!
    if (epoch /= 2000.d0) then
       call xcc_dp_precess(ivector,epoch,2000.d0,hvector)
    else
       hvector(1:3)=ivector(1:3)
    endif
!
!     Multiply vector by transpose of transformation matrix:
!
    DO I = 1,3
       ovector(i) = sum(hvector(1:3)*T(1:3,i))
    ENDDO
!
    RETURN
  END Subroutine xcc_DP_E_TO_G
!===============================================================

  Subroutine xcc_DP_G_TO_E(ivector,epoch,ovector)
!
!    Routine to convert galactic coordinates to ecliptic (celestial) 
!  coordinates at the given epoch.  First the conversion to ecliptic 
!  2000 is done, then if necessary the results are precessed.
!
!  03/08/88, a.c.raugh, STX.
!  11/16/01, C.Burigana, D.Maino, OAT
!
!
!  Arguments:
!
    real(dp) ::    IVECTOR(3)      ! Input coordinate vector
    real(dp) ::    EPOCH           ! Input coordinate epoch
    real(dp) ::    OVECTOR(3)      ! Output coordinate vector
!
!  Program variables:
!
    real(dp) ::    HVECTOR(3)      ! Temporary holding place for coordinates
    real(dp) ::    T(3,3)          ! Galactic coordinate transformation matrix
    integer(i4b) ::  I,J             ! Loop variables
!
    Data T / -0.054882486d0, -0.993821033d0, -0.096476249d0, &! 1st column
         &    0.494116468d0, -0.110993846d0,  0.862281440d0, &! 2nd column
         &   -0.867661702d0, -0.000346354d0,  0.497154957d0/ ! 3rd column
!
!     Multiply by transformation matrix:
!
    DO I = 1,3
       hvector(i) = sum(ivector(1:3)*T(i,1:3))
    ENDDO
!
!     Precess coordinates to epoch if necesssary:
!
    if (epoch /= 2000.d0) then
       call xcc_dp_precess(hvector,2000.d0,epoch,ovector)
    else
       ovector(1:3)=hvector(1:3)
    endif
!
    return
  end Subroutine xcc_DP_G_TO_E
!===============================================================

  Subroutine xcc_DP_Q_TO_E(ivector,epoch,ovector)
!
!    Routine to convert equatorial (celestial) coordinates at the given
!  epoch to ecliptic coordinates.  Adapted from the COBLIB routine by Dr.
!  E. Wright.
!
!  02/01/88, a.c.raugh, STX.
!  11/16/01, C.Burigana, D.Maino, OAT
!
!
!  Arguments:
!
    real(dp) ::    IVECTOR(3)     ! Input coordinate vector
    real(dp) ::    EPOCH          ! Epoch of equatorial coordinates
    real(dp) ::    OVECTOR(3)     ! Output coordinate vector
!
!  Program variables:
!
    real(dp) ::    T              ! Centuries since 1900.0
    real(dp) ::    EPSILON        ! Obliquity of the ecliptic, in degrees
    real(dp) ::    HVECTOR(3)     ! Temporary holding place for output vector
    real(dp) ::  dc, ds
!
!  Set-up:
!
    T = (epoch - 1900.d0) / 100.d0
    epsilon = 23.452294d0 - 0.0130125d0*T - 1.63889d-6*T**2 &
         &  + 5.02778d-7*T**3
!
!  Conversion:
!
    dc = cos(DTOR * epsilon)
    ds = sin(DTOR * epsilon)
    hvector(1) = ivector(1)
    hvector(2) = dc*ivector(2) + ds*ivector(3)
    hvector(3) = dc*ivector(3) - ds*ivector(2)
!
!  Move to output vector:
!
    ovector(1:3) = hvector(1:3)
!
    return
  end Subroutine xcc_DP_Q_TO_E
!===============================================================

  Subroutine XCC_DP_PRECESS(ivector,iepoch,oepoch,ovector)
!
!    Subroutine to precess the input ecliptic rectangluar coordinates from
!  the initial epoch to the final epoch.  First the precession around the
!  Z-axis is done by a simple rotation.  A second rotation about the 
!  X-axis is then performed to account for the changing obliquity of the
!  ecliptic.
!
!  02/25/88, a.c.raugh, STX.
!  03/22/88, acr: modified to include changing epsilon
!  11/16/01, Carlo Burigana & Davide Maino, OAT
!
!
!  Parameters:
!
    real(dp) ::   IVECTOR(3)  ! Input coordinates
    real(dp) ::   IEPOCH      ! Initial epoch, years (input)
    real(dp) ::   OEPOCH      ! Final epoch, years (input)
    real(dp) ::   OVECTOR(3)  ! Output coordinates
!
!  Routine variables:
!
    real(dp) ::   TEMP        ! Holding place for rotation results
    real(dp) ::   GP_LONG     ! General precession in longitude
    real(dp) ::   dE          ! Change in the obliquity of the ecliptic
    real(dp) ::   OBL_LONG    ! -Longitude about which the obliquity is changing
    real(dp) ::   dL          ! Longitude difference
    real(dp) ::   TM          ! Mid-epoch, in tropical centuries since 1900
    real(dp) ::   TVECTOR(3)  ! Temporary holding vector
    real(dp) :: dco, dso, dce, dse, dcl, dsl
!
    Tm = ((oepoch+iepoch)/2.d0 - 1900.d0) / 100.d0
    gp_long  = (oepoch-iepoch) * (50.2564d0+0.0222d0*Tm) / 3600.d0
    dE       = (oepoch-iepoch) * (0.4711d0-0.0007d0*Tm) / 3600.d0
    obl_long = 180.d0 - (173.d0 + (57.06d0+54.77d0*Tm) / 60.d0) &
         &   + gp_long / 2.d0
    dL       = gp_long - obl_long
!
!     Z-axis rotation by OBL_LONG:
!
    dco = cos(DTOR * obl_long)
    dso = sin(DTOR * obl_long)
    tvector(1) = ivector(1)*dco - ivector(2)*dso
    tvector(2) = ivector(1)*dso + ivector(2)*dco
    tvector(3) = ivector(3)
!
!     X-axis rotation by dE:
!
    dce = cos(DTOR * dE)
    dse = sin(DTOR * dE)
    temp       = tvector(2)*dce - tvector(3)*dse
    tvector(3) = tvector(2)*dse + tvector(3)*dce
    tvector(2) = temp
!
!     Z-axis rotation by GP_LONG - OBL_LONG:
!
    dcl = cos(DTOR * dL)
    dsl = sin(DTOR * dL)
    temp       = tvector(1)*dcl - tvector(2)*dsl
    tvector(2) = tvector(1)*dsl + tvector(2)*dcl
    tvector(1) = temp
!
!     Transfer results to output vector:
!
    ovector(1:3) = tvector(1:3)

    return
  end Subroutine XCC_DP_PRECESS
!

end module coord_v_convert


