	INTEGER*4 FUNCTION FUT_QGET_ANC (ANC_LUN, ANC_REC)

C---------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FUT_QGET_ANC
C
C	PROGRAM DESCRIPTION:
C	  Program to read FIRAS Ancillary Housekeeping Archive
C
C	AUTHOR:
C	  R.Kummerer	(Adapted from work done by E.Fung)
C	  STX
C	  FEBRUARY 13, 1989
C
C	Input parameters:
C	  ANC_LUN	I*2	Logical unit for reading Ancillary Housekeeping
C				Archive
C
C	Output parameters:
C	  ANC_REC	RECORD	Buffer containing Ancillary Housekeeping record
C
C	INPUT FILES:
C	  FIRAS ANCILLARY HOUSEKEEPING ARCHIVE
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

	DICTIONARY 'NFS_ANC'
	RECORD /NFS_ANC/ ANC_REC

	INTEGER		*2	ANC_LUN
	INTEGER		*2	ISTAT
	INTEGER		*2	CT_STAT(20)

	EXTERNAL		FUT_NORMAL
	EXTERNAL		FUT_EOF
	EXTERNAL		FUT_CTQGET_ERR


	CALL CT_QUERY_GET (, ANC_LUN, ANC_REC, CT_STAT)
	IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	  IF (CT_STAT(1) .EQ. CTP_ENDOFFILE) THEN
	    FUT_QGET_ANC = %LOC(FUT_EOF)
	  ELSE
	    CALL LIB$SIGNAL (FUT_CTQGET_ERR, %VAL(1), %VAL(CT_STAT(1)))
	    FUT_QGET_ANC = %LOC(FUT_CTQGET_ERR)
	  ENDIF
	  RETURN
	ENDIF

	FUT_QGET_ANC = %LOC(FUT_NORMAL)

	RETURN
	END
