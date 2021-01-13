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
module head_fits

  !-------------------------------------------------------------------------
  ! this module (head_fits) introduced in Healpix 1.2
  ! merges the module wrap_fits present in version 1.1
  ! and the modules strip_fits that appeared in some intermediate versions
  !
  ! It is designed to read and write header information in FITS files
  !
  ! EH, 2002-08-09
  !-------------------------------------------------------------------------

  ! subroutine add_card   [interface]
  ! subroutine write_hl
  ! subroutine merge_headers
  !
  ! subroutine get_card   [interface]
  !
  ! subroutine del_card

  USE healpix_types
  USE misc_utils
  interface add_card
     ! add_card generic routine
     ! ------------------------
     ! adds a card to a fits header
     ! syntax = call add_card(header, keyword, value [, comment])
     !   header  = character string array
     !   keyword = character string (NON BLANK)
     !   value   = can be either logical, integer, real, double, character
     !   comment = (optional) character string scalar
     !
     module procedure d_add_card, f_add_card, a_add_card, i_add_card, l_add_card, v_add_card
  end interface

  interface get_card
     ! get_card generic routine
     ! ------------------------
     ! read a keyword (case UNsensitive) from a FITS header
     ! syntax = call get_card(header, keyword, value [, comment, count])
     !   header  = (input) character string array
     !   keyword = (input) character string scalar
     !   value   = (output) can be either logical, integer, real, double, character
     !              contains on output value of keyword
     !   comment = (output, optional) character string scalar
     !              contains on output value of comment field
     !   count   = (output, optional) integer scalar,
     !       is set to 1 if the keyword is found and 0 otherwise
     !
     ! except if value is declared as a character string,
     ! its type should match that in the fits file
     !
     module procedure d_get_card, f_get_card, a_get_card, i_get_card, l_get_card
  end interface

  interface del_card
     ! del_card generic routine
     ! ------------------------
     ! deletes a card from a fits header
     ! syntax = call del_card(header, keywords) or
     !          call del_card(header, keyword)
     !   header   = character string array
     !   keywords = character string vector (NON BLANK)
     !   keywords = character string  (NON BLANK)
     !
     module procedure del_card1, del_cardn
  end interface

  private

  character (LEN=1 ), private :: dtype
  character (LEN=20), private :: stval
  character (LEN=80), private :: card, stcom
  integer(kind=I4B) , private :: nlh, statusfits, count_in = 0
  logical(kind=LGT) , private :: match, exact, casesen = .false., verbose = .true.

  public :: add_card, merge_headers

  public :: get_card

  public :: del_card

contains

  !=====================================================================
  subroutine d_add_card(header, kwd, value, comment) ! double precision
    real(kind=DP), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment

    write(stval,'(1pe20.12)') value
    call write_hl(header, kwd, stval, comment)
    RETURN
  END subroutine d_add_card
  !=====================================================================
  subroutine f_add_card(header, kwd, value, comment) ! single precision
    real(kind=SP), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment

    write(stval,'(1pe20.8)') value
    call write_hl(header, kwd, stval, comment)
    RETURN
  END subroutine f_add_card
  !=====================================================================
  subroutine a_add_card(header, kwd, value, comment) ! character
    character (LEN=*), intent(IN), OPTIONAL :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    character (LEN=240) :: st_value, st_comment

    st_value   = ''
    st_comment = ''
    if (present(value)  ) then
       write(st_value,'(a)')   value
       st_value = adjustl(st_value)
       if (st_value(1:1) /= "'" .and. &
            & trim(kwd) /= "COMMENT" .and. trim(kwd) /= "HISTORY") st_value = "'"//trim(st_value)//"'"
    endif
    if (present(comment)) write(st_comment,'(a)') comment
    call write_hl(header, kwd, st_value, st_comment)
    RETURN
  END subroutine a_add_card
  !=====================================================================
  subroutine i_add_card(header, kwd, value, comment) ! integer (i*4)
    integer(kind=I4B), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment

    write(stval,'(i20)') value
    call write_hl(header, kwd, stval, comment)
    RETURN
  END subroutine i_add_card
  !=====================================================================
  subroutine l_add_card(header, kwd, value, comment) ! logical
    logical(kind=LGT), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment

    write(stval,'(l7)') value
    call write_hl(header, kwd, stval, comment)
    RETURN
  END subroutine l_add_card
  !=====================================================================
  subroutine v_add_card(header) ! blank line
    character (LEN=80), dimension(:), intent(INOUT) :: header

    call write_hl(header, 'COMMENT', ' ', ' ')
!    call write_hl(header, ' ', ' ', ' ')
  END subroutine v_add_card
  !=====================================================================
  subroutine write_hl(header, kwd, st_value, comment)
    IMPLICIT none
    character (LEN=80), dimension(1:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL  :: st_value
    character (LEN=*), intent(IN), OPTIONAL  :: comment
    integer(kind=I4B) :: hdtype
    integer(kind=I4B) :: iw, j, j_low, j_hi

    character (LEN=240) :: headerline
    character (LEN=80)  :: buffheader
    character (LEN=10)  :: pad10
    !=====================================================================

    nlh = size(header)
    ! iw = first blank line
    iw = nlh + 1
    do
       if (iw == 1) exit
       if (trim(header(iw-1)) /= '') exit
       iw = iw - 1
    enddo

    if (iw > nlh) goto 10

    pad10='   '
    buffheader =''
    headerline = adjustl(kwd)
    if (present(st_value)) headerline = trim(headerline)//' '//adjustl(st_value)
    if (present(comment))  headerline = trim(headerline)//' '//TRIM(adjustl(comment))

    if (trim(headerline) == 'COMMENT') then ! COMMENT alone
       header(iw) = 'COMMENT'
       iw = iw + 1
       goto 20
    endif

    hdtype = 0
    statusfits = 0
    CALL ftgthd(headerline(1:79), buffheader, hdtype, statusfits)
    if (hdtype == -1) then ! delete keyword
       header(iw) = headerline(1:79)
    else ! append
       header(iw) = buffheader
    endif
    iw = iw + 1
    ! deal with long cards
    do j=0,2
       j_low = 80 + j*70
       j_hi  = min(j_low+69, len(headerline))
       if (len_trim(headerline) >= j_low) then
          if (iw > nlh) goto 10
          hdtype = 0
          statusfits = 0
          CALL ftgthd(pad10//headerline(j_low:j_hi), buffheader, hdtype, statusfits)
          header(iw) = buffheader
          iw = iw + 1
       else
          exit
       endif
    enddo
    goto 20

10  continue
    print*,'WARNING: Header is too short, card not written'

20  continue

    RETURN
  END subroutine write_hl
  !=====================================================================
  subroutine merge_headers( header1, header2)
    IMPLICIT none
    character (LEN=80), dimension(1:), intent(IN)    :: header1
    character (LEN=80), dimension(1:), intent(INOUT) :: header2
    integer(kind=I4B) :: iw1, iw2, s1, s2, ss

    s2 = size(header2)
    iw2 = s2
    do while(header2(iw2) == '' .and. iw2 > 1)
       iw2 = iw2 - 1
    enddo
    iw2 = iw2 + 1

    s1 = size(header1)
    iw1 = s1
    do while(header1(iw1) == '' .and. iw1 > 1)
       iw1 = iw1 - 1
    enddo
    iw1 = iw1 + 1

    ss = MIN(iw1, s2-iw2+1)

    if (ss < iw1) then
       print*,' Second header in merge_headers is not long enough'
       print*,' Should be ',iw1+iw2-2,' instead of ',s2
       print*,' Output will be truncated'
    endif

    header2(iw2:iw2+ss-1) = header1(1:ss)

    RETURN
  END subroutine merge_headers



  !===================================================================
  subroutine d_get_card(header, kwd, value, comment, count) ! double precision
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    real(kind=DP),                    intent(OUT) :: value
    character (LEN=*),                intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),                intent(OUT), OPTIONAL :: count
    integer :: i

    count_in = 0
    value = 0.0d0
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'F' .or. dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             count_in = 1
             if (present(comment)) comment = stcom
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected DOUBLE (F), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//kwd//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine d_get_card

  !===================================================================
  subroutine f_get_card(header, kwd, value, comment, count) ! single precision
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    REAL(kind=SP),                    intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer :: i

    count_in = 0
    value = 0.
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'F' .or. dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected REAL (F), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//kwd//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine f_get_card

  !===================================================================
  subroutine a_get_card(header, kwd, value, comment, count) ! ascii string
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN) :: kwd
    character (LEN=*),                intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count

    integer :: ifor, ibac, i

    count_in = 0
    value = ' '
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          stval = adjustl(stval)
          ! remove first and last quote
          ifor = index(stval,"'")
          ibac = index(stval,"'",back=.true.)
          if (ifor >= 1         ) stval(ifor:ifor) = " "
          if (ibac <= len(stval) .and. ibac > ifor) &
               &                  stval(ibac:ibac) = " "
          value = trim(adjustl(stval))
          count_in = 1
          if (present(comment)) comment = stcom
          if (present(count))   count = count_in
          return
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//kwd//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine a_get_card

  !===================================================================
  subroutine i_get_card(header, kwd, value, comment, count) ! integer
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN) :: kwd
    integer(kind=I4B),                          intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer(kind=I4B) :: i

    count_in = 0
    value = 0
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected integer (I), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//kwd//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine i_get_card

  !===================================================================
  subroutine l_get_card(header, kwd, value, comment, count) ! logical
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    logical(kind=LGT),                intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer(kind=I4B) :: i

    count_in = 0
    value = .false.
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'L') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected logical (L), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//kwd//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine l_get_card

  !===================================================================
  subroutine del_cardn(header, kwds) ! remove cards
    !===================================================================
    ! remove for header the card corresponding to keyword kwds
    ! kwds can be a vector
    !===================================================================
    character(len=80), dimension(1:), intent(INOUT) :: header
    character(len=*),  dimension(1:), intent(IN)    :: kwds
    integer(kind=I4B) :: i
    character(len=20) :: kwd
    !===================================================================
    nlh = size(kwds)
    do i = 1, nlh
       kwd = adjustl(kwds(i))
       if (trim(kwd) /= "") then
          kwd = "- "//kwd
          call write_hl(header, kwd)
       endif
    enddo
  end subroutine del_cardn
  !===================================================================
  subroutine del_card1(header, kwds) ! remove cards
    !===================================================================
    ! remove for header the card corresponding to keyword kwds
    ! kwds can be a vector
    !===================================================================
    character(len=80), dimension(1:), intent(INOUT) :: header
    character(len=*),                 intent(IN)    :: kwds
    character(len=20) :: kwd
    !===================================================================
       kwd = adjustl(kwds)
       if (trim(kwd) /= "") then
          kwd = "- "//kwd
          call write_hl(header, kwd)
       endif
  end subroutine del_card1

end module head_fits
