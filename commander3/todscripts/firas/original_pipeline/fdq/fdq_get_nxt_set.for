	INTEGER*4 FUNCTION FDQ_GET_NXT_SET (READ_FLAG, DONE, CT_LUN, 
	1	SCI_REC, NACTIVE, ACTIVE, NPROC, PROC_CHAN, 
	2	AVG_TIME, EARL,NSCI_READ, SET_SEC, IFG_COUNT)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_NXT_SET
C/
C/	PROGRAM DESCRIPTION:
C/	  The things this subroutine does were formerly in-line code in
C/	  FDQ.  It does the following:
C/	    (1) Get the next 4 SCI records to be processed, ready the
C/	        archives if necessary;
C/	    (2) Determine how many channels are "active" (this means that
C/	        there is at least one record to be processed; an inactive
C/	        channel means there are no more records in the current 
C/	        segment); this means variable NACTIVE and array ACTIVE;
C/	    (3) Determine which of the active channels are to be processed
C/	        in this pass (iteration of the do-loop in FDQ).  Channels
C/	        are considered to be in the same "group" if there times are
C/	        within 16 seconds of each other; this is variable NPROC and
C/	        array PROC_CHAN;
C/	    (4) Determine the average time of the Science channels that
C/	        are being processed.  This average time will be used for
C/	        the Engineering record.
C/
C/	AUTHOR:
C/	  Edwin Fung
C/	  GSFC
C/	  April, 1987
C/
C/	MODIFIED BY:
C/
C/         J. Durachta
C/             ARC
C/          10/22/87
C/           Reason:    Musnt tamper with the science record time (as 
C/                      FDQ_WITHIN has been doing). Local variable earl_time
C/                      added to accomplish this. FDQ_WITHIN function now
C/                      passes dummy varibles earl_time and bntm.
C/                       
C/	  Shirley M. Read
C/	      STX
C/	  January 12, 1988
C/	    Reason: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal. Added parameter, Set_Sec, to
C/	                calling sequence for setting the number of seconds for
C/	                the science record group.
C/
C/	  R. Kummerer
C/	      STX
C/	  October 28, 1988
C/	    Reason: 	SPR 2702, Convert AVG_TIME to CT "precision".
C/
CH	CHANGE LOG:
CH
CH	Version 4.2.1 02/09/89, SPR 2700, Shirley M. Read, STX
CH		FDQ falls over for out of order and bad time ranges.The FIRAS
CH	 	Stripper now writes the science records using the telemetry
CH		minor frame time at transmit for the primary key in order to
CH		ensure a valid time tag. A new FIRAS preprocessor reads the raw
CH		science records, computes the midpoint of collect time and flags
CH		any records with bad midpoint times. FDQ must now access the
CH		computed collect time from the new field in the science record
CH	        in order to group the four science channels with the bracketing
CH		housekeeping frames and set the data quality summary flag when
Ch		the badtime flag is set. The bad record is written back to the
CH		archive immediately so that the next record may be considered
CH	        for the group. The average midpoint of collect times is still
CH		used as the time tag of the FDQ engineering record and is thus
CH		stored in the science record as the cross reference. However,
CH		the transmit times of the science records are stored in the 
CH		FDQ engineering record. FDQ must now enter the science gain into
CH		each science record from the bracketing housekeeping frames. The
CH		two LS bits must now be checked in the certification mask when
CH		running FDQ. The bit must be set for FPP but not for FDQ.
CH
CH	Version 4.4.1 06/15/89, SPR 3961, Shirley M. Read, STX
CH              The attitude in the raw science record is left zero on IFGs
CH		failed for bad telemetry by FDQ. This causes problems in 
CH		quicklook programs which need the pixel number but do not
CH		check the telemetry quality or time tags.
CH
CH	Version 4.4.1 06/16/89, SPR 3964, Shirley M. Read, STX
CH              FDQ fails to make an entry in the IFG matrix for IFGs failed
CH		by FPP. All four FPR_SCCTL records are expected by the 
CH		FDQ_Write_Matrix function in the calling sequence. The
CH		FDQ_Get_Nxt_Set only passed one channel record at a time.
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after
CH		an FPP or FDQ abort. The FPR tracking system compounds the
CH		problems. Files with non-matching version numbers seem often
CH		to result from improper clean-up. Bad record times cause
CH		SEGCTL to abort and mess up the tracking system. It was
CH		decided to change the modify of the science records in FPP
CH		and FDQ to a simple COBETRIEVE read of the existing records
CH		from a dataset and write a modifed dataset with the same
CH		information which was entered on the modify. Two new science
CH		data sets will be required: a science dataset of raw science
CH		data plus FPP input and a science dataset with FPP science
CH		data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH		xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where
CH		xx is the channel id. The new datasets must be opened and
CH		processed in FPP and FDQ. 
CH	Version 4.4.1 08/23/89, SER 4210, R. Kummerer, STX
CH		Prevent raw science segment overlaps related changes.
CH
CH      Modified by:
CH                  H. Wang, STX, 1/29/91
CH                  New requirements for FDQ
CH                  FDQ_WRITE_MATRIX no longer be needed 
CH                  Because the tracking file no longer be required
CH                  
C/	CALLING SEQUENCE:
C/	  Return_Status = FDQ_GET_NXT_SET (READ_FLAG, DONE, CT_LUN, SCI_REC,
C/	    NACTIVE, ACTIVE, NPROC, PROC_CHAN, AVG_TIME, nsci_read, 
C/	    SET_SEC, ifg_count)
C/
C/	Input/Output parameters:
C/
C/	  READ_FLAG(4)  I,O	Flag(one for each channel) to indicate whether 
C/				one should call COBETRIEVE to read another 
C/				record. If true go ahead and read.
C/	  DONE(4)	I,O	Flag(one for each channel) to indicate whether
C/				the channel has more records to process.  If
C/				true means there are no more.
C/	  CT_LUN(22)	I	Logical unit numbers for COBETRIEVE access;
C/				CT_LUN(1) through CT_LUN(4) are for the Science
C/				archives:  RH, RL, LH and LL respectively.
C/	  SCI_REC(4)	I,O	Science record buffers(one for each channel)
C/				containing the last record read from archive.
C/	  NACTIVE	O	Number of active channels (active channels are
C/				those that are not "done") yet
C/	  ACTIVE(4)	O	Array of "active" channel ID's
C/	  NPROC		O	Number of channels "being processed".  This
C/				is the no. of channels within 16 seconds of
C/				the "earliest" channel including the earliest
C/				channel.
C/	  PROC_CHAN(4)	O	Array of "being processed" channel ID's
C/	  AVG_TIME	O	Average start time of all "being processed"
C/				channels in VAX binary format.
C/	  NSCI_READ(4)	I,O	These are counters for how many SCI records are
C/				read (for each channel) for a summary report at
C/				the end of DQ processing.
C/        EARL          O       The early science record
C/	  SET_SEC	I       Number of seconds for science set grouping.
C/	  IFG_COUNT(4)  I,O     Number of IFG records modified
C/
C/	INPUT/OUTPUT FILES:
C/	  COBETRIEVE FIRAS Science archives
C/
C/	SUBROUTINES CALLED:
C/	  CT_READ_ARCV (COBETRIEVE routine)
C/	  CT_WRITE_ARCV (COBETRIEVE routine)
C/	  FDQ_WITHIN
C/	  TIME_LT (COBETRIEVE function)
C/	  FUT_ATTITUDE
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	ERROR HANDLING:
C/	  Message to SYS$OUTPUT
C/
C/	METHOD USED:
C/
C/	  PDL for FDQ_GET_NXT_SET (April 1987, E.F.)
C/
C/	  Set NACTIVE to 0;
C/	  Set EARL to -1;
C/	  Do for channel = 1 to 4;
C/	    If (READ_FLAG(channel)) and (channel is not "done") then
C?	      DO WHILE Not EOF and Science record has a bad midpoint time
C/		  Read next SCI record for channel;
C/	        If EOF then
C/	          Mark channel as "done";
C/	        Else
C/		  If ( Midpoint_Time is good ) then
C/	            Note binary midpoint time;
C/	            Increment NACTIVE;
C/	            Set ACTIVE(NACTIVE) to channel;
C/		    Note the good value for the badtime flag;
C/		  Else
C/		    Increment the Ifg_Count for the channel;
C/		    Insert the Ifg_Count in the science record;
C/		    Set bit 5 of data quality summary flag in science record;
C?		    Call FUT_Attitude with attitude = none.
C/		    Write the modified science record back in the archive;
C/		  Endif
C/	        Endif;
C/	      ENDDO
C/	    Else (if channel is not "done") then
C/	       Increment NACTIVE;
C/	       Set ACTIVE(NACTIVE) to channel;
C/	    Endif;
C/	    If (Channel not "done") then
C/	      If (EARLISET = -1) or (channel time is earlier than "EARL" time) then
C/	        Set EARL to channel;
C/	      Endif;
C/	    Endif;
C/	  Enddo;		(Channel = 1 to 4)       
C/
C/
C/	  Set NBEING_PROC to 0;
C/	  Set DIFF_SUM to 0;
C/
C/	  If (more than one active channel) then
C/	    Do for all "active" channel
C/	      If (start time of channel is within 16 seconds of "EARL" channel) then
C/	        Mark channel as "being processed";
C/	        Save the time difference (between channel time and "EARL");
C/	      Endif;
C/	    Enddo;
C/	    Get "average channel time" from the time differences;
C/	  Elseif (only one active channel) then
C/	    Set "average channel time" to the one active channel time;
C/	  Endif;
C/
C/	  Return;
C/	  End.
C/
C/

	IMPLICIT	NONE

	INCLUDE 'CT$LIBRARY:CTUSER.INC'
	INCLUDE '($SSDEF)'
	INCLUDE '(FUT_PARAMS)'

	DICTIONARY	'NFS_SDF'
	LOGICAL		TIME_LT		! COBETRIEVE library function
	LOGICAL		FDQ_WITHIN	! DQ function
	INTEGER*4       FUT_ATTITUDE      ! FUT library attitude function

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ERADDX
	EXTERNAL	FDQ_WITHINER
	EXTERNAL	FDQ_WITHINTRUE
	EXTERNAL	FDQ_WITHINFALS
	EXTERNAL	FDQ_CTREADERR
	EXTERNAL	FDQ_CTWRITERR

	INTEGER*4       LIB$ADDX        ! System function - add quadword
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4       ZERO / 0 /
	integer*4	STATUS		! Dummy status variable
	Byte            BAD_TIME_VAL / 32 /

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			    !
!     Passed parameters     !
!			    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LOGICAL*1	DONE(4), READ_FLAG(4)
	INTEGER*2	ACTIVE(4), CT_LUN(22), NACTIVE, NPROC,
	1		nsci_read(4), PROC_CHAN(4)
	INTEGER*4	AVG_TIME(2)
	CHARACTER*14	GMT
	INTEGER*4       SET_SEC
	INTEGER*4       IFG_COUNT(4)
	RECORD		/NFS_SDF/ SCI_REC(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*2	BIGGER, CHAN, CT_STAT(20), EARL, I
	INTEGER*4	AVG_DIFF(2), BNTM(2), DIFF(2), 
	1		DIFF_SUM(2), GRP_THRES(2), earl_time(2)
	INTEGER*4       ATT_TYPE 		! Type of attitude = none

!	DATA		GRP_THRES(1) / 160000000 /	! Threshold is 16 seconds
							! for determining if 2 channels
							! are in the "same group"
	LOGICAL*1       FIRST / .TRUE. /, EOF, BADTIME

	Character*12	Datasets(22)
	Data Datasets / 'FPP_SDF_RH',
	1		'FPP_SDF_RL',
	1		'FPP_SDF_LH',
	1		'FPP_SDF_LL',
	1		'NFS_HKP',
	1		'FDQ_IDX',
	1		'FDQ_ENG',
	1		'FDQ_ETR',
	1		'FEX_AV_CALRS',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		'FDQ_SDF_RH',
	1		'FDQ_SDF_RL',
	1		'FDQ_SDF_LH',
	1		'FDQ_SDF_LL' /


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

C	First time set the number of second for group and attitude type.

	If ( First ) then

	  GRP_THRES(1) = SET_SEC * 10000000

	  ATT_TYPE = fac_none
	  First = .False.
	Endif

	NACTIVE = 0
	EARL = -1

	DO I = 1, 4
	 IF ( RETSTAT .EQ. SUCCESS ) THEN
	  IF (READ_FLAG(I)) THEN	! All DONE channel should have READ_FLAG set to false
	    EOF = .FALSE.
	    BADTIME = .TRUE.

	    DO WHILE ((.not. EOF) .and. (BADTIME) .and. (RETSTAT .EQ. SUCCESS))
	      CALL CT_READ_ARCV (, CT_LUN(I), SCI_REC(I), CT_STAT)
	
	      READ_FLAG(I) = .FALSE.

	      IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
		EOF = .TRUE.
	        IF (CT_STAT(1) .EQ. CTP_ENDOFFILE) THEN
	          DONE(I) = .TRUE.
	          EOF = .TRUE.
	        ELSE
	          RETSTAT = ERROR
	          STATUS = CT_STAT(1)
	          CALL LIB$SIGNAL(FDQ_CTREADERR,%VAL(2),%VAL(STATUS),
	1	    DATASETS(I))

	        ENDIF		! IF (CT_STAT .EQ. CTP_ENDOFFILE)
	      ELSE		! Good read from CT_READ_ARCV
	        nsci_read(i) = nsci_read(i) + 1
	        IF ( SCI_REC(I).COLLECT_TIME.BADTIME_FLAG .eq. Zero) Then
	          CALL LIB$MOVC3 (8, SCI_REC(I).COLLECT_TIME.MIDPOINT_TIME,
	1	      BNTM)
	          NACTIVE = NACTIVE + 1
	          ACTIVE(NACTIVE) = I
	 	  BADTIME = .FALSE.	
	        ELSE		! Bad midpoint of collect time
	          IFG_COUNT(I) = IFG_COUNT(I) + 1
	          SCI_REC(I).DQ_DATA.IFG_NO = IFG_COUNT(I)
	          SCI_REC(I).DQ_DATA.DATA_QUALITY(110) = BAD_TIME_VAL
		  STATUS = FUT_ATTITUDE( SCI_REC(I), ATT_TYPE,
	1	    SCI_REC(I).CT_HEAD.GMT, SCI_REC(I).CT_HEAD.GMT )
	          CALL CT_WRITE_ARCV(,CT_LUN(I+18),SCI_REC(I),CT_STAT)
	          IF (CT_STAT(1). ne. CTP_NORMAL) THEN
	            RETSTAT = ERROR
	            STATUS = CT_STAT(1)
	            CALL LIB$SIGNAL(FDQ_CTWRITERR,%VAL(2),%VAL(STATUS),
	1				DATASETS(I+18))
		  ENDIF		! CT_STAT(1)
	        ENDIF           ! Good or bad midpoint of collect time
	      ENDIF		! IF (CT_STAT .NE. CTP_NORMAL)
	    ENDDO
	  ELSE			! READ_FLAG is false
	    if (.not. DONE(i)) then
	      NACTIVE = NACTIVE + 1
	      ACTIVE(NACTIVE) = I
	    endif
	  ENDIF			! IF (READ_FLAG(I))

	  if (.not. DONE(i) .and. RETSTAT .EQ. SUCCESS) then
	    if ( (EARL .EQ. -1) .OR. 
	1	TIME_LT ( SCI_REC(I).COLLECT_TIME.MIDPOINT_TIME, 
	1		SCI_REC(EARL).COLLECT_TIME.MIDPOINT_TIME) ) THEN
	      EARL = I
	    endif
	  endif
	 ENDIF			! Return status is 'Success'.
	ENDDO			! I = 1, 4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								       !
!     Now that we have found out what is "active", let us find	       !
!     out what is "being processed", to be put in the PROC_CHAN array     !
!								       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF ( RETSTAT .EQ. SUCCESS ) THEN
	  NPROC = 0
	  DIFF_SUM(1) = 0
	  DIFF_SUM(2) = 0
          call lib$movc3(8,sci_rec(earl).collect_time.midpoint_time,earl_time)
	ENDIF

	IF ( RETSTAT .EQ. SUCCESS ) THEN
	  IF (NACTIVE .GT. 1) THEN
	    DO I = 1, NACTIVE
	     IF ( RETSTAT .EQ. SUCCESS ) THEN
	      CHAN = ACTIVE(I)
              call lib$movc3(8,sci_rec(chan).collect_time.midpoint_time,bntm)
	      STATUS = FDQ_WITHIN (bntm, earl_time, GRP_THRES, DIFF, 
	1	 BIGGER ) 
	      IF ( STATUS .eq. %LOC(FDQ_WITHINTRUE) ) then
	        NPROC = NPROC + 1
	        PROC_CHAN(NPROC) = CHAN
	        STATUS = LIB$ADDX (DIFF, DIFF_SUM, DIFF_SUM)
	        IF ( STATUS .NE. SS$_Normal ) THEN
	          CALL LIB$SIGNAL(FDQ_ERADDX,%VAL(1),%VAL(STATUS))
		  RETSTAT = ERROR
		ENDIF
	      ELSEIF ( STATUS .NE. %LOC(FDQ_WITHINFALS) ) THEN
 	        RETSTAT = ERROR
	        If ( status .eq. %LOC(FDQ_ABERR)) status = zero
	        CALL LIB$SIGNAL(FDQ_WITHINER,%VAL(1),%VAL(STATUS))
	      ENDIF		! IF (FDQ_WITHIN)
	     ENDIF		! Return status is 'Success'.
	    ENDDO

!     Get "average" time difference

            IF ( RETSTAT .EQ. SUCCESS ) THEN    
	      AVG_DIFF(1) = NINT ( FLOAT(DIFF_SUM(1)) / FLOAT(NPROC))  
	      STATUS = LIB$ADDX (SCI_REC(EARL).COLLECT_TIME.MIDPOINT_TIME, 
	1	AVG_DIFF, AVG_TIME)
	      IF ( STATUS .NE. SS$_Normal ) THEN
	        CALL LIB$SIGNAL(FDQ_ERADDX,%VAL(1),%VAL(STATUS))
		RETSTAT = ERROR
              ENDIF
	      CALL CT_BINARY_TO_GMT(AVG_TIME,GMT)
	      CALL CT_GMT_TO_BINARY(GMT,AVG_TIME)
	    ENDIF  	! Return status is 'Success'
	  ELSE
	    IF ( RETSTAT .EQ. SUCCESS ) THEN
	      NPROC  = 1
	      PROC_CHAN(1) = ACTIVE(1)
	      CALL LIB$MOVC3 (8, SCI_REC(EARL).COLLECT_TIME.MIDPOINT_TIME, 
	1	  AVG_TIME)
	    ENDIF
	  ENDIF		! If nactive gt 1
	ENDIF		! If return status is 'Success'

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_NXT_SET = %LOC(FDQ_NORMAL)
	ELSE
	  FDQ_GET_NXT_SET = %LOC(FDQ_ABERR)
	ENDIF

	RETURN
	END
