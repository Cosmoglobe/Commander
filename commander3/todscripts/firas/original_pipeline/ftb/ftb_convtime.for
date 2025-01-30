ch


	PROGRAM FTB_CONVTIME
!
!	Utility to convert GMT to Binary
!	and vice versa
!
!  	Include files used:
!	    $SSDEF
!	    CT$LIBRARY:CTUSER.INC
!
!-------------------------------------------------------------------------
!
! Changes:
!
!	R. Kummerer, December 3, 1986. Change input format from hex to
!	    decimal.
!
!  	Shirley M. Read, STX, August 26, 1988
!	    Reason: SPRs 1566 and 1763 requested that all FIRAS facilities 
!           should tell of successful completion and signal errors via
!	    $status for batch processing. Also added the interface to
!	    Fut_Error, the condition handler.
!
!-------------------------------------------------------------------------

	IMPLICIT	NONE

!	Include Files

	INCLUDE         '($SSDEF)'
	INCLUDE		'CT$LIBRARY:CTUSER.INC'

!	FTB messages

	EXTERNAL   	FTB_NORMAL
	EXTERNAL	FTB_ABERR

!	Condition Handler

	EXTERNAL        FUT_ERROR

!	Local Variables

	LOGICAL*1	CHOICE
	CHARACTER	GMT*14
	CHARACTER	ASC_TIME*23
	INTEGER*4	BINTIM(2)
	INTEGER*4	I
	INTEGER*2	TL
	INTEGER*4	CVTFLG
	INTEGER*4	STATUS

!	Functions

	INTEGER*4	STR$UPCASE
	INTEGER*4	SYS$BINTIM
	INTEGER*4	SYS$ASCTIM
!
!     Code begins here
!
!     Establish condition handler.            

	CALL LIB$ESTABLISH ( FUT_ERROR )

	STATUS = SS$_NORMAL

!	Process Loop

5	CONTINUE

	CALL LIB$ERASE_PAGE(1, 1)

	WRITE (6, 6001)
6001	FORMAT (//' Enter "B" to go from Binary to GMT (default)'/
	1	  '       "G" to go from GMT to Binary' /
	1	  '       "I" to go from Binary to ASCII' /
	1	  '       "T" to go from ASCII to GMT' //
	1	  ' 		===>> ', $)

	read (5, 5001) choice
5001	format (a)

	IF (CHOICE .EQ. 'G' .or. choice .eq. 'g') then
	  write (6, 6002)
6002	  format (//' Enter GMT as YYDDDHHMMSSTTT' /
	1	    '          --> ', $)
	  READ (5, 5001) GMT
	  CALL STR$TRANSLATE (GMT, GMT, '0', ' ')

	  CALL CT_GMT_TO_BINARY (GMT, BINTIM)
	  write (6, 6003) bintim
6003	  format (//' Binary time (1): ', I /
	1	    '             (2): ', I )
!6003	  format (//' Binary time (1): ', Z8.8 /
!	1	    '             (2): ', Z8.8 )
	else if (choice .eq. 'I' .or. choice .eq. 'i') then
	  write (6, 6004)
	  read (5, 5002) bintim(1)
	  write (6, 6005)
	  read (5, 5002) bintim(2)

	  status = sys$asctim(tl,asc_time,bintim,cvtflg)
	  if (status .ne. ss$_normal) then
	    call lib$signal(%val(status))
	  else 
	    write (6, 8001) asc_time
8001	    format (//' ASCII: ', a/)
	  endif
	else if (choice .eq. 'T' .or. choice .eq. 't') then
	  write (6, 7001)
7001      format (//' Enter ASCII time as DD-MON-YYYY HH:MM:SS.CC' /
	1	    '                 --> ', $)
	  read (5, 5001) asc_time
	  status = str$upcase(asc_time,asc_time)
	  if (status .ne. ss$_normal ) then
	    call lib$signal(%val(status))
	  else
	    status = sys$bintim(asc_time,bintim)
	    if (status .eq. ss$_normal) then
	      call ct_binary_to_gmt(bintim,gmt)
	      write (6, 6006) (gmt(I:I), I=1,14)
	    else
	      call lib$signal(%val(status))
	    endif
	  endif
	else
	  write (6, 6004)
6004	  format (//' Enter first(lower order) binary in decimal: ',$)
	  read (5, 5002) bintim(1)
5002	  format (I)
!5002	  format (Z8.8)
	  write (6, 6005)
6005	  format (' Enter second(higher order) binary in decimal: ',$)
	  read (5, 5002) bintim(2)

	  call ct_binary_to_gmt (bintim, gmt)
	  write (6, 6006) (gmt(I:I), I=1,14)
6006	  format (//' GMT: ', 2A, '-', 3A, '-',
	1	3(2A, '-'), 3A1)

	ENDIF

	WRITE (6, 6999)
6999	FORMAT (///' More times to convert (Y/[N])? ', $)

	read (5, 5001) choice
	if (choice .eq. 'Y' .or. choice .eq. 'y') goto 5

!       Exit the program.

	if ( status .eq. ss$_normal ) then
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
	endif
	end
	  
	end
