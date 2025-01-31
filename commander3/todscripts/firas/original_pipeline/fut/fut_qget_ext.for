	INTEGER*4 FUNCTION FUT_QGET_EXT (EXT_LUN, EXT_REC)

C---------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FUT_QGET_EXT
C
C	PROGRAM DESCRIPTION:
C	  Program to read FIRAS Engineering Extrema archive
C
C	AUTHOR:
C	  R.Kummerer	(Adapted from work done by E.Fung)
C	  STX
C	  March 31, 1989
C
C	Input parameters:
C	  EXT_LUN	I*2	Logical unit for reading Engineering Extrema
C				Archive
C
C	Output parameters:
C	  EXT_REC	RECORD	Buffer containing engineering extrema
C
C	INPUT FILES:
C	  FIRAS Engineering Extrema
C
C	OUTPUT FILES:
C	  NONE
C
C	SUBROUTINES CALLED:
C	  CT_QUERY_GET
C
C	INCLUDE FILES USED:
C	  CT$LIBRARY:CTUSER.INC
C
C---------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'CT$LIBRARY:CTUSER.INC'

	DICTIONARY 'FXT_ENG_XTRM'
	RECORD /FXT_ENG_XTRM/ EXT_REC

	INTEGER		*2	EXT_LUN
	INTEGER		*2	ISTAT
	INTEGER		*2	CT_STAT(20)

	EXTERNAL		FUT_NORMAL
	EXTERNAL		FUT_EOF
	EXTERNAL		FUT_CTQGET_ERR


	CALL CT_QUERY_GET (, EXT_LUN, EXT_REC, CT_STAT)
	IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	  IF (CT_STAT(1) .EQ. CTP_ENDOFFILE) THEN
	    FUT_QGET_EXT = %LOC(FUT_EOF)
	  ELSE
	    CALL LIB$SIGNAL (FUT_CTQGET_ERR, %VAL(1), %VAL(CT_STAT(1)))
	    FUT_QGET_EXT = %LOC(FUT_CTQGET_ERR)
	  ENDIF
	  RETURN
	ENDIF

	FUT_QGET_EXT = %LOC(FUT_NORMAL)

	RETURN
	END
