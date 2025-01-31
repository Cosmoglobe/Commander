	SUBROUTINE FEP_READ_FIRAS_HKP (HKP_LUN, HKP_REC, CATALOG_ENTRY, 
	1	ISTAT)
C/
C/	PROGRAM NAME:
C/	  READ_FIRAS_HSK
C/
C/	PROGRAM DESCRIPTION:
C/	  Program to read FIRAS Housekeeping Archive
C/
C/	AUTHOR:
C/	  E.Fung
C/	  GSFC
C/	  OCTOBER 31, 1985
C/
C/	MODIFIED BY:
C/	  E.FUNG
C/	  GSFC
C/	  APRIL 21, 1986
C/	  REASON:	To report to calling pgm of catalog entry number.
C/
C/	Input parameters:
C/	  HKP_LUN	I*2	Logical unit for reading Housekeeping Archive
C/
C/	Output parameters:
C/	  HKP_REC		RECORD	Buffer containing Housekeeping record
C/	  CATALOG_ENTRY		I*4	CT Catalog Entry #
C/	  ISTAT			I*2	1 = Normal status
C/					-1 = No more HSK records to read
C/					Anything else: Bad status
C/
C/	INPUT FILES:
C/	  FIRAS HOUSEKEEPING ARCHIVE
C/
C/	OUTPUT FILES:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  CT_READ_ARCV
C/
C/	INCLUDE FILES USED:
C/	  HSKP.INC
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	ERROR HANDLING:
C/	  Through return status ISTAT.
C/

	IMPLICIT	NONE

	DICTIONARY 'NFS_HKP'
	RECORD /NFS_HKP/ HKP_REC

	INTEGER*2	HKP_LUN, ISTAT

	integer*4	catalog_entry

	INCLUDE		'CT$LIBRARY:CTUSER.INC'

	INTEGER*2	ct_stat(20)

	INTEGER*4	CAT_ENTRY_EQV

	EQUIVALENCE	(CAT_ENTRY_EQV,	CT_STAT(11))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ISTAT=1
	CALL CT_READ_ARCV (, HKP_LUN, HKP_REC, CT_STAT)
	IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	  IF (CT_STAT(1) .EQ. CTP_ENDOFFILE) THEN
	    ISTAT = -1
	    RETURN
	  ELSE
	    ISTAT = CT_STAT(1)
	    RETURN
	  ENDIF
	ENDIF

	CATALOG_ENTRY = CAT_ENTRY_EQV

	RETURN
	END
