	integer*4 function fep_decode_dbname (name, inum, jchan)

C------------------------------------------------------------------------
C    PURPOSE: Decodes the database name fetched from FIRASCALRES.INP.
C	      Returns the decoded channel and calres group number.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: R. Isaacman
C            ARC
C
C    INVOCATION: STATUS = FEP_DECODE_DBNAME ( DBNAME, NUM, CHAN )
C
C    INPUT PARAMETERS:
C	DBNAME		CH*8		STOL database word.
C
C    OUTPUT PARAMETERS: 
C	NUM		I*4		Calres group number.
C	CHAN		I*4		Channel ID.
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: None
C
C----------------------------------------------------------------------

	implicit	none

	character	*8	name
	integer		*4	inum
	integer		*4	jchan
	integer		*2	istat	

	external	fep_normal
	external	fep_decoderead

	read (name(5:5),10,iostat=istat) inum		! inum=1,2,3,4
10	format (i1)
	if (istat .ne. 0) then
	   fep_decode_dbname = %loc(fep_decoderead)
	   call lib$signal(fep_decoderead,%val(1),%val(istat))
	   return
	end if

	jchan = 1					! A side, LO curr
	if (name(2:2) .eq. 'B') jchan = 3		! B side
	if (name(7:7) .eq. 'H') jchan = jchan + 1	! HI current

	fep_decode_dbname = %loc(fep_normal)

	return
	end
