	INTEGER*4 FUNCTION FUT_CHECKSUM ( SCI_REC, FLG_CKSM_ERR_ST)
C/
C/ PROGRAM NAME: 
C/	FUT_CHECKSUM
C/
C/ PROGRAM DESCRIPTION:
C/	This routine will perform a checksum of the 512 points of the Science
C/	interferogram data and compare it to the data checksum in the channel
C/	microprocessor header. 
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  December 11, 1987.
C/
C/ MODIFIED BY:
C/	None
C/
C/ CALLING SEQUENCE:
C/	Status = FUT_Checksum (Sci_Rec, Flg_Cksm_Err_St)
C/
C/ INPUT PARAMETERS:
C/	Sci_Rec -- Raw Science Record of dictionary structure Nfs_Sdf. 
C/	Flg_Cksm_Err_St  -- Data checksum quality flag 
C/
C/ OUTPUT PARAMETERS:
C/
C/ ENTRY/EXIT:
C/	Normal function entry and return.
C/
C/ ERROR HANDLING
C/	
C/ INPUT/OUTPUT DISK FILES: 
C/	None
C/
C/ PROCS/MACROS/INCLUDE FILES:
C/
C/ SUBROUTINES CALLED:
C/	Checksum from utility library.
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED: 
C/      Compute the data checksum on the 512 interferogram points of the
C/	Science Record. Discard the two higher order bytes. Compare this
C/      computed checksum to the Science Record Microprocessor_Header 
C/	Data_Checksum. If they do not match, set the input Data_Quality 
C/	flag to 1.
C/	Set function value to success for return.
C/
C/ PDL for FUT_CHECKSUM:
C/ 
C/	BEGIN
C/
C/	CALL CHECKSUM to PERFORM the following:
C/
C/	  SET the four byte integer checksum accumulator to zero
C/
C/	  DO from 1 to 512 interferogram points
C/
C/	    ADD the interferogram point value to the checksum accumulator
C/
C/	  ENDDO
C/
C/	  Discard the two higher order bytes for the computed checksum
C/
C/	Check the return status and if Ok complete the next steps
C/
C/	IF ( the computed checksum is not equal to the Science_Record
C/	     Microprocessor_Header Data_Checksum ) THEN
C/
C/	  SET the Data_Quality_Flag to 1
C/
C/	ENDIF
C/
C/	SET the return status for the function to success or failure.
C/
C/	RETURN
C/	END
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE
	
	DICTIONARY 'NFS_SDF'
	RECORD /NFS_SDF/ SCI_REC  ! Input Raw Science record for a channel

	EXTERNAL   	FUT_NORMAL
	EXTERNAL	FUT_ABERR
	EXTERNAL	FUT_LIBCHECKSM

	INCLUDE 	'($SSDEF)'

	INTEGER*4	CHECKSUM		! Utility Library function
 
	BYTE FLG_CKSM_ERR_ST    ! Input quality flag for checksum

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable

	INTEGER*2 	NUM_POINTS
	PARAMETER       ( NUM_POINTS = 512 )	! Number of points in each IFG

	INTEGER*2 	CHECKSUM_VALUE  	! Value returned from Checksum

C**	Set return status to success.

	RETSTAT = SUCCESS

C**	Call CHECKSUM to do the following:	
C**	Compute the checksum of all the 512 inteferogram points.
C**	Discard the two higher order bytes and  return Checksum_Value.

	STATUS = CHECKSUM ( NUM_POINTS, SCI_REC.IFG_DATA.IFG,
	1	 CHECKSUM_VALUE )
	IF ( STATUS .NE. SS$_NORMAL) THEN
	  RETSTAT = ERROR
	  CALL LIB$SIGNAL(FUT_LIBCHECKSM,%VAL(1),%VAL(STATUS))
	ENDIF

C**	Compare Checksum_Value to the Data Checksum in the Science Record 
C**	Micro Header Word 6. If not equal, then set the quality flag to 1.
 
	IF ( CHECKSUM_VALUE .NE. SCI_REC.SCI_HEAD.SC_HEAD6 ) THEN

	    FLG_CKSM_ERR_ST = 1

	ENDIF

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FUT_CHECKSUM = %loc(FUT_NORMAL)
	ELSE
	  FUT_CHECKSUM = %loc(FUT_ABERR)
	ENDIF

	RETURN
	END

