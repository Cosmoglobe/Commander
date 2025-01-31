	INTEGER*4 FUNCTION FUT_QGET_EMF (EMF_LUN, EMF_REC)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FUT_QGET_EMF
C
C	PROGRAM DESCRIPTION:
C	  Program to read FIRAS Enginerring Mode Archive
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  November 14, 1987   (adapted from work done by E. Fung)
C
C	Input parameters:
C	  EMF_LUN	I*2	Logical unit for reading engineering mode archive
C
C	Output parameters:
C	  EMF_REC	RECORD	Buffer containing engineering mode record
C
C	INPUT FILES:
C	  FIRAS ENG MODE ARCHIVE
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
C------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'CT$LIBRARY:CTUSER.INC'

	DICTIONARY	'NFS_EMF'
	RECORD /NFS_EMF/ EMF_REC

	INTEGER		*2	EMF_LUN
	INTEGER		*2	ISTAT
	INTEGER		*2	CT_STAT(20)

	EXTERNAL		FUT_NORMAL
	EXTERNAL		FUT_EOF
	EXTERNAL		FUT_CTQGET_ERR


	CALL CT_QUERY_GET (, EMF_LUN, EMF_REC, CT_STAT)
	IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	  IF (CT_STAT(1) .EQ. CTP_ENDOFFILE) THEN
	    FUT_QGET_EMF = %LOC(FUT_EOF)
	  ELSE
	    CALL LIB$SIGNAL (FUT_CTQGET_ERR, %VAL(1), %VAL(CT_STAT(1)))
	    FUT_QGET_EMF = %LOC(FUT_CTQGET_ERR)
	  ENDIF
	  RETURN
	ENDIF

	FUT_QGET_EMF = %LOC(FUT_NORMAL)

	RETURN
	END
