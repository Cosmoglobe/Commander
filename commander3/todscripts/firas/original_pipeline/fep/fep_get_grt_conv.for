	Integer*4 Function FEP_Get_GRT_Conv (  GRT_LUN,dbword,
     .						ohms, ctemps, npts )


C------------------------------------------------------------------------
C    PURPOSE: Retrieves GRT reference data  of converting GRT counts to
C	      engineering units from the archive
C
C    Arthour: Harte Wang
C              STX
C              Dec. 30, 1989        
C
C
C
C    INVOCATION: STATUS = FEP_GET_GRT_CONV ( GRT_LUN, DBWORD,
C					     OHMS, CTEMPS, NPTS )
C
C    INPUT PARAMETERS:
C       GRT_lun                 I*4             Fortran unit number for GRT REF
C                                               File
C	DBWORD			CH*8		STOL database word.
C
C    OUTPUT PARAMETERS: 
C	OHMS(*)			R*4		Counts to ohms convertion table.
C	CTEMPS(*)		R*4		Ohms to temps convertion table.
C	NPTS			I*4		Number of points in convertion
C						tables.
C
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	FUT_Error
C	$SSDef
C  
C  
C----------------------------------------------------------------------
c PDL: 
c  Set return status to success
c  Read GRT table reference file by data_base key word
C
C   Read reference data by GRT field database name from FDB_GRTRICH to buffer
C   If error
C   Then
C     Write error message to screen
C     Set status to error
C   Else 
C    								  
C     Fill in: the number of points in convertion table,
C              Counts to ohms convertion table
C              ohms to temps convertion table
C              with info read from FDB_GRTRICH file     
C								  
C   Endif   
C   IF  status equal error
C   Then
C    Set return status for this function to bad 
C    Write error message tothe screen
C   Else
C    Set return status for this function to success
C   Endif   
C
C   Return
C   End
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_grtpars)'
	Include		'(FUT_Error)'
	Include		'($SSDef)'

	Character	*8	dbword
	Integer		*4	npts
	Integer		*4	Grt_lun
	Real		*4	ohms(400)
	Real		*4	ctemps(400)

	Integer		*4	status


	External		FEP_Normal
	External		FEP_GRTREFREAD
	External		FUT_Normal
c
c Get the attributes of the field being converted.
c

c
c Open the STOL database word / GRT table reference file.

	INCLUDE		'SYS$LIBRARY:FORIOSDEF/NOLIST'
	BYTE		GRT_BUFF(GRT_RECSZ)
	INTEGER*4       NPTS_EQV,IST
	EQUIVALENCE	(NPTS_EQV,	GRT_BUFF(GRT_DB_KNCT))	! # points in lookup table


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							   !
!     The next will load the GRT lookup tables     !
!							   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   	  FEP_Get_GRT_Conv = %Loc(FEP_Normal)

	  READ (GRT_LUN, KEY=Dbword, IOSTAT=IST) GRT_BUFF
	  IF (IST .NE. 0) THEN
	        call lib$signal(FEP_GRTREFREAD,%val(2),%val(IST),
	1	'FDB_GRTRICH')
                FEP_GET_GRT_CONV=%Loc(FEP_GRTREFREAD)
                Return
          ELSE						! GOOD READ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								  !
!     Fill in  GRT_TBL with info read from DB file     !
!								  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    NPTS = NPTS_EQV		      
	    CALL LIB$MOVC3 (4*NPTS_EQV, 	
	1		GRT_BUFF(GRT_DB_KNTS),
	1		OHMS(1))

	    CALL LIB$MOVC3 (4*NPTS_EQV, 	! Move all temperatures values
	1		GRT_BUFF(GRT_DB_KNTS + 4 * NPTS_EQV),
	1		ctemps(1))

          Endif
	FEP_Get_GRT_Conv = %Loc(FEP_Normal)
	Return
	End
