	INTEGER*4 FUNCTION FDQ_MAKE_IDX (NEW_SEGMENTS, ENG_REC, SCI_REC,
	1	HSK_ANALOG, IPDU_STAT, stmon_sync, IDX_LUN, NCHAN,
	1	VALID_CHAN, IDXTOLS, IDXFLAGS, CUR_IDX_REC,TELM_FLAG,IDX_WRIT)
C/
C/	PROGRAM NAME:
C/	  FDQ_MAKE_IDX
C/
C/	PROGRAM DESCRIPTION:
C/	  This is the routine to construct the index record and determine
C/	  if a new index record ought to be written.
C/
C/	AUTHOR:
C/	  E. FUNG,
C/	  GSFC
C/	  March 30, 1986
C/
C/	MODIFIED BY:
C/	  E. FUNG,
C/	  GSFC
C/	  June 7, 1987
C/	  REASON:	To add new passed parameter NEW_SEGMENTS.  This 
C/			parameter is to be used to tell FDQ_NEW_IDX
C/			(a function called by FDQ_MAKE_IDX) to start with
C/			a brand new index record (disregard the old one,
C/			which would be the index record of the last
C/			processed segment).
C/			Another change done at this time is to change the
C/			record LAST_IDX_REC from a local record to a common
C/			record (in named COMMON-block LAST_IDX).  This is
C/			done so that LAST_IDX_REC is accessible to other
C/			routine, mainly FDQ, which has to write out the
C/			last IDX record when a segment is done.  This is
C/			incorrectly handled in the past in that the last index
C/			record will NOT be written out if there are no changes
C/			in any of the index variables (that exceed the tolerances).
C/
C/	D.WARD
C/	GSFC/STX
C/	OCT. 13,1987
C/	REASON:	     GET_IPDU_FIELDS was deleted.
C/		     Call GET_IPDU_FIELDS to get the IPDU individual statuses;
C/
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 14, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  June, 1988
C/	  REASON: 	Added input records IDXTOLS and IDXFLAGS which
C/			are now obtained from CCT_Get_Config.
C/
C/	MODIFIED BY:
C/	  H. Wang
C/	  STX
C/	  Jan. 20, 1991
C/	  REASON:	New requirements for FDQ
C/              The following new quantities will be tested for the generation
C/              of idex records
C/              * Bolometer CMD bias
C/              * IPDU relays
C/              * LVDT Status
C/              * MTM CAL Motor
C/              
C/
C/	  Larry P. Rosen, STX
C/	  30 August 1990
C/			Pass Report_File name and Report flag from call to
C/			FDQ_NEW_IDX for constructing its standardized file name.
C/
C/        SPR 7855, Missing gain data in fdq_idx
C/                  H. Wang, Dec. 12, 1990
C/
C/	CALLING SEQUENCE:
C/	  Status = FDQ_MAKE_IDX (ENG_REC, SCI_REC,HSK_ANALOG, IPDU_STAT,
C/		stmon_sync, IDX_LUN, NCHAN, VALID_CHAN, IDXTOLS, IDXFLAGS,
C/		CUR_IDX_REC,TELM_FLAG,IDX_WRIT)
C/
C/	INPUT PARAMETERS:
C/	  NEW_SEGMENTS		L*1	Flag to indicate whether we have
C/					just started a new set of segments
C/					or not.  This flag is not being consulted
C/					in FDQ_MAKE_IDX, but it should be passed
C/					down to the next level -- function
C/					FDQ_NEW_IDX for deciding to start
C/					with a brand new index record.
C/	  ENG_REC(1024) 	BYTE	Engineering record constructed under the current pass
C/
C/	  SCI_REC(1536, 4)	BYTE	The 4 Science file buffers
C/
C/	  HSK_ANALOG(52)	BYTE	HSK analog counts
C/
C/	  IPDU_STAT(8)		BYTE	IPDU statuses
C/
c/	  stmon_sync(2)		byte	Status monitor sync:
c/					stmon_sync(1) = Side A
c/					stmon_sync(2) = Side B
c/					if in sync stmon_sync = 1
c/					if not in sync, stmon_sync = 0
c/					(This is determined in routine FDQ_GET_BRK_HSK)
C/	  IDX_LUN		I*2	CT logical unit for index file
C/
C/	  NCHAN			I*2	# valid channel in this pass
C/
C/	  VALID_CHAN(4)		I*2	Array of valid channel ID's:
C/					1=RH, 2=RL, 3=LH, 4=LL
C/	  IDXTOLS               Rec     Index tolerances.
C/	  IDXFLAGS              Rec     Index check enable flags.
C/        TELM_FLAG             I*2     Telemetry quality of eng. record
C/                                      0 = good
C/                                      1 = bad   
C/	OUTPUT PARAMETER:
C/	  CUR_IDX_REC(512)		BYTE	Index record constructed in this pass
C/        IDX_WRIT              I*4     The number of new idx records
C/	INPUT/OUTPUT FILES:
C/	  FIRAS INDEX ARCHIVE (OUTPUT)
C/
C/	INCLUDE FILES USED:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  FDQ_GET_ENG_FLDS
c/	  FDQ_GET_HSK_FLDS
C/	  FDQ_GET_SCI_FLDS
C/	  FDQ_NEW_IDX (I*4 FUNCTION)
C/	  CT_WRITE_ARCV
C/
C/	ERROR HANDLING:
C/	  TBD
C/
C/	METHOD USED:
C/	  The following is the PDL:
C/
C/	  Call FDQ_GET_SCI_FLDS to get the fields of IDX record derived
C/		from the Science records;
C/	  Call FDQ_GET_ENG_FLDS to get the fields of IDX record derived
C/		from the Engineering record;
C/
C/	  If (FDQ_NEW_IDX (current IDX record, old IDX record, cur ENG time) then
C/	    Write out current IDX record;
C/	    Move current IDX record to "old" IDX record;
C/	  Endif;
C/
C/	  Return;
C/	  End.
C/

	IMPLICIT	NONE

	INCLUDE		'CT$LIBRARY:CTUSER.INC'

        dictionary 'fdq_idx'
        record /fdq_idx/    cur_idx_rec
	common /last_idx/   last_idx_rec
        record /fdq_idx	/    last_idx_rec

        dictionary 'fdq_eng'
        record /fdq_eng/    eng_rec

        dictionary 'nfs_sdf'
        record /nfs_sdf/    sci_rec(4)

        dictionary 'fex_idx_tols'
        record /fex_idx_tols/    idxtols
        dictionary 'fex_idx_flag'
        record /fex_idx_flag/    idxflags

	EXTERNAL	FDQ_NEWIDXER
	EXTERNAL	FDQ_CTWRITERR
	EXTERNAL	FDQ_NORMAL
	EXTERNAL	FDQ_WRITEIDX
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_GETENGERR
	EXTERNAL	FDQ_GETHSKERR
	EXTERNAL	FDQ_GETSCIERR
	
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable

	INTEGER*4	FDQ_NEW_IDX
	INTEGER*4	FDQ_GET_ENG_FLDS
	INTEGER*4	FDQ_GET_HSK_FLDS
	INTEGER*4	FDQ_GET_SCI_FLDS

        byte		hsk_analog(52), 
	1		IPDU_STAT(8), new_segments,
	1		STMON_SYNC(2)

	INTEGER*2	IDX_LUN, NCHAN, VALID_CHAN(4), TELM_FLAG, CT_STAT(20)

        integer*4       i              ! a counter
        INTEGER*4       IDX_WRIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C	Set return status to success.

	RETSTAT = SUCCESS
        IF (new_segments) IDX_WRIT = 0

	CUR_IDX_REC.idx_sync.side(1) = STMON_SYNC(1)
	CUR_IDX_REC.idx_sync.side(2) = STMON_SYNC(2)

	STATUS = FDQ_GET_SCI_FLDS (SCI_REC,NCHAN,
	1	VALID_CHAN, CUR_IDX_REC)
	IF ( STATUS .NE. %loc(FDQ_NORMAL) ) THEN
	  call lib$signal(FDQ_GETSCIERR,%val(1),%val(status))
	  RETSTAT = ERROR
	ENDIF

	IF ( RETSTAT .EQ. SUCCESS) THEN
	  STATUS = FDQ_GET_ENG_FLDS (ENG_REC, CUR_IDX_REC)
	  IF ( STATUS .NE. %loc(FDQ_NORMAL) ) THEN
	    call lib$signal(FDQ_GETENGERR,%val(1),%val(status))
	    RETSTAT = ERROR
	  ENDIF
	ENDIF

	IF ( RETSTAT .EQ. SUCCESS ) THEN 
	  STATUS = FDQ_GET_HSK_FLDS (HSK_ANALOG, IPDU_STAT, CUR_IDX_REC)       
	  IF ( STATUS .NE. %loc(FDQ_NORMAL) ) THEN
	    call lib$signal(FDQ_GETHSKERR,%val(1),%val(status))
	    RETSTAT = ERROR
	  ENDIF
	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									  !
!     At this point the index record CUR_IDX_REC should be all filled     !
!     in.  Compare this with the last IDX record written		  !
!									  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF ( RETSTAT .EQ. SUCCESS ) THEN

	  STATUS = FDQ_NEW_IDX( new_segments, CUR_IDX_REC, LAST_IDX_REC,
	1	eng_rec, idxtols, idxflags, telm_flag)

	  If ( STATUS .EQ. %loc(FDQ_WRITEIDX) ) THEN

!	    Set return status to success.
	    RETSTAT = SUCCESS
!
!     Fill in end time of LAST_IDX_REC
!
           last_idx_rec.header.gmt_end_time_ascii = eng_rec.ct_head.gmt 

           do i = 1,2
             last_idx_rec.header.binary_end_time(i) = eng_rec.ct_head.time(i)
           end do

!     Write to Index Archive

	   LAST_IDX_REC.HEADER.DATASET_ID=CTU_$FIR_IDX

	   Call CT_WRITE_ARCV (, IDX_LUN, LAST_IDX_REC, CT_STAT)
	   IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	      RETSTAT = ERROR
	      STATUS = CT_STAT(1)
	      CALL LIB$SIGNAL(FDQ_CTWRITERR,%VAL(2),%VAL(STATUS),
	1	'FDQ_IDX')
           ELSE
             IDX_WRIT = IDX_WRIT + 1
	   ENDIF

!     Fill in start time of CUR_IDX_REC

           cur_idx_rec.header.gmt_start_time_ascii = eng_rec.ct_head.gmt 

           do i = 1,2
             cur_idx_rec.header.binary_start_time(i) = eng_rec.ct_head.time(i)
           end do

!     Move buffer (that is just written) to LAST_IDX_REC

	   Call LIB$MOVC3 (512, CUR_IDX_REC, LAST_IDX_REC)

	  ELSEIF ( STATUS .EQ. %loc(FDQ_NORMAL)) THEN

	    RETSTAT = SUCCESS
	  ELSE
	    RETSTAT = ERROR
	  ENDIF
	ENDIF		! Retstat is succcess

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_MAKE_IDX = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_MAKE_IDX = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
