	INTEGER*4 FUNCTION FDQ_GET_BRK_HSK (HSK_LUN, SCI_BNTM, HSK_RECS,
	1	FIRST_MJ, stmon_sync, IPDU_RELAY, ISTAT, CT_STAT ,first_time)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_BRK_HSK
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine takes the Science start time as argument
C/	  and will successively read Housekeeping records until
C/	  one or two Housekeeping records is (or are) found
C/	  which contains 2 major frames that will "bracket" the
C/	  Science start time.
C/
C/	  This program assume that the next Housekeeping record
C/	  ready to be read when control comes to this program is
C/	  always BEFORE the Science time passed as input argument.
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  October 22, 1985
C/
C/	MODIFIED BY:
C/	  E. FUNG
C/	  GSFC
C/	  October 30, 1985
C/	  REASON:	Add WRITE statements for information and debug messages.
C/
C/	  E. FUNG
C/	  GSFC
C/	  NOVEMBER 5, 1985
C/	  REASON:	Add variable NREAD to count how many CT reads are
C/			performed for each entry into this subroutine.
C/
C/	  E.FUNG
C/	  November 6, 1985
C/	  REASON:	To add logic to check for no HSK records before a given
C/			Science record and flag it.
C/
C/	  E. FUNG
C/	  November 18,1985
C/	  REASON:	To change the status returns of EOF_HSK from 3 to 6
C/			and BRK_MJ_MISS from 4 to 7 to conform with calling
C/			program (DATA_QUALIFY).
C/
C/	  E. FUNG
C/	  December 7, 1985
C/	  REASON:	To change the writing of the log file (logical unit 16)
C/			to D-lines.
C/			Also, a new input parameter LOGGING, to see if any logging
C/			at all is desired (to logical unit 6).
C/
C/	  E.FUNG
C/	  GSFC
C/	  April 4, 1986
C/	  REASON:	To add the check for the toggling of the major frame
C/			ID in the status monitor status and pass the info to
C/			the calling program (which is main program DAQUALIFY);
C/			also pick up the IPDU power relay statuses and passed
C/			it back as an output parameter.
C/
C/	  E.FUNG
C/	  GSFC
C/	  May 7, 1986
C/	  REASON:	To take away the LOGGING parameter (no longer there
C/			in DQ 7.1) and put in d-lines for writing out the
C/			SCI and HSK times
C/
C/	  E.FUNG
C/	  GSFC
C/	  JUNE 11, 1986
C/	  REASON:	To add a write statement whenever a new HSK segment
C/			is being used.
C/
C/	  E.FUNG
C/	  GSFC
C/	  MAY 19, 1987
C/	  REASON:	Added variable NREAD_GOOD to count the no. of successful
C/			CT reads.
C/			Also correct a potential problem with HSK later than
C/			the first SCI time (quite common!):  In the check
C/			of DIFF look at DIFF(2) as well!
C/
C/        J. W. Durachta
C/        August 28, 1987
C/            ARC
C/        Reason:       First_time flags added to correct for problem when first
C/                      science record of segment fails to have hosekeeping data
c/                      associated with it. Upon such a condition, the major
C/                      frame counter is reset to "one", and processing cont-
C/                      inues as though that first science record had not been
C/                      present. Use of the CDD and record structures have also
C/                      been added. Housekeeping record is new format (ie. 576
C/                      bytes).
C/
C/	D.WARD
C/	GSFC/STX
C/	OCT. 13, 1987
C/	REASON:		Name was changed to FDQ_GET_BRK_HSK from 
C/			GET_BRACKETING_HSK.
C/
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 14, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	CALLING SEQUENCE:
C/	  Return Status =FDQ_GET_BRK_HSK (HSK_LUN, SCI_BNTM, HSK_RECS, FIRST_MJ,
C/		STMON_SYNC, IPDU_RELAY, ISTAT, CT_STAT, first_time )
C/
C/	INPUT PARAMETERS:
C/	  HSK_LUN	I*2	Logical unit number for reading the Housekeeping file
C/	  SCI_BNTM(2)	I*4	Binary time for Science time tag; when running
C/				with real instrument this should be the time
C/				calculated from the microprocessor headers.
C/	  FIRST_MJ	I*2	If FIRST_MJ = 0 it means FDQ_GET_BRK_HSK
C/				is being called the first time; otherwise
C/				it should point to one of the 4 major frames
C/				in HSK_RECS which is the first major frame
C/				for the LAST Science record process.  Note
C/				that this is also an output parameter.
C/        FIRST_TIME    L*1     Logical flag to indicate that this is the first
C/                              record in the current segment.
C/
C/	OUTPUT PARAMETERS:
C/	  HSK_RECS(1152) BYTE	This is an array that can hold 2 Housekeeping
C/				records.  Depending on SCI_BNTM and FIRST_MJ,
C/				either 1 or 2 Housekeeping record(s) will be
C/				read into it.
C/	  FIRST_MJ	I*2	On exit of this subroutine, this should hold
C/				the pointer to one of the 4 major frames in
C/				HSK_RECS which contains the first major of
C/				housekeeping data for the Science record that
C/				is currently being processed.
C/	  STMON_SYNC(2)	BYTE	Status monitor sync:
C/				0 = out of sync (FIRAS major frame ID not toggling)
C/				1 = in sync (FIRAS major frame ID toggling)
C/				STMON_SYNC(1) is for Side A
C/				STMON_SYNC(2) is for Side B
C/	  IPDU_RELAY(8)	BYTE	The IPDU power relay statuses
C/	  ISTAT		I*2	Return status from FDQ_GET_BRK_HSK:
C/				1 = GOOD
C/				2 = BAD
C/				6 = EOF on Housekeeping
C/				7 = Either 1 or 2 "bracketing" major frames missing
C/	  CT_STAT(16)	I*2	The status information returned from the
C/				last COBETRIEVE read.
C/
C/	INPUT FILES:
C/	  FIRAS Housekeeping Archive
C/
C/	OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	SUBROUTINES CALLED:
C/	  CT_READ_ARCV (COBETRIEVE subroutine)
C/	  GMT_TO_BINARY (COBETRIEVE subroutine)
C/	  TIME_LT (COBETRIEVE function)
C/	  LIB$SUBX
C/
C/	ERROR HANDLING:
C/	  Error info passed in ISTAT and CT_STAT.
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/
C/	  If (FIRST_MJ=0) then
C/	    Read HSK into first record of HSK_RECS;
C/	    Mark first record of HSK_RECS as "new";
C/	  Endif
C/	  Set BRKT_HSK_NOT_FOUND to true;
C/
C/	  Do while (BRKT_HSK_NOT_FOUND)
C/	    Increment FIRST_MJ by 1;
C/	    If FIRST_MJ=3 then
C/	      Mark first record of HSK_RECS as old;
C/	      If second record of HSK_RECS is old then
C/	        Read HSK into second record of HSK_RECS;
C/	      Endif
C/	    Elseif FIRST_MJ=5 then
C/	      Mark second record of HSK_RECS as old;
C/	      If first record of HSK_RECS is old then
C/	        Read HSK into first record of HSK_RECS;
C/	        Mark first record of HSK_RECS as "new";
C/	      Endif
C/	      Set FIRST_MJ to 1
C/	    Endif
C/	    Compare time tag of HSKP record pointed to by FIRST_MJ 
C/			with Science time;
C/	    If time is after Science time then
C/	      Set BRKT_HSK_NOT_FOUND to false (to exit main do-while loop);
C/	    Endif
C/	  Enddo
C/
C/	  Decrement FIRST_MJ by 1 (since it is pointing to the record
C/		that just comes after the Science time);
C/	  If FIRST_MJ=0 then
C/	    Set FIRST_MJ to 4;
C/	  Endif;
C/	  If bracketing major frames are across record boundaries (would be
C/		so if FIRST_MJ=2 or FIRST_MJ=4 then
C/	    Subtract the times between the bracketing major frames;
C/	    If difference is greater than 64 second then
C/	      Set status to "either 1 or 2 brking major frames missing";
C/	      Return;
C/	    Endif;
C/	  Endif;
C/
C/	  Determine status monitor sync by check status monitor major frame ID;
C/	  Move IPDU power relay from 1st major frame to output parameter IPDU_RELAY;
C/
C/	  Set status to good;
C/	  Return;
C/	  End.
C/

	IMPLICIT	NONE

        dictionary 'nfs_hkp'
        record/nfs_hkp/ hsk_recs(2)

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_CTREADERR
	EXTERNAL 	FDQ_ERSUBX
	EXTERNAL	FDQ_NHKPREAD
	EXTERNAL	FDQ_NEWHKPSEG
	INTEGER*4	Lib$Subx

	LOGICAL*1	IPDU_RELAY(8), STMON_SYNC(2),first_time

	INTEGER*2	CT_STAT(20), FIRST_MJ,
	1		HSK_LUN, ISTAT

	INTEGER*4	SCI_BNTM(2)

	PARAMETER	BAD = 2

	parameter	brk_mj_miss = 7

	parameter	eof_hsk = 6

	PARAMETER	GOOD = 1

	INCLUDE		'CT$LIBRARY:CTUSER.INC'
	INCLUDE		'($SSDEF)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LOGICAL*1	BRKT_HSK_NOT_FOUND,
	1		NEW(2), SCI_GMT(14),
	1		TIME_LT

	INTEGER*2	cur_cat_entry /0/, NEXT_MJ, nread_good /0/
	integer*2	nread			! This is a count of how many CT read
						! for each entry into this subroutine

	INTEGER*2	STMON_A(2), STMON_B(2)	! Temperory words to hold status monitor statueses

	INTEGER*4	DIFF(2), tot_reads,
	1		HSK_BNTM(2), HSK_BNTM1(2),
	1		SEC64 / 640000000 /, mj(2), bn(2)

	INTEGER*2	I, J, K

        character*14    hsk_gmt
        character*14    gmt1
        character*14    gmt2

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

        if(first_time)then
          tot_reads = 0   !Keep track of total number of reads per segment.
          first_time = .false.
        endif

	nread = 0			! Reset NREAD to 0 for each entry
	IF (FIRST_MJ .EQ. 0) THEN	! First time through for this subroutine

!
!     First time always read HSKP into first record in HSK_RECS
!

	  CALL CT_READ_ARCV (, HSK_LUN, HSK_RECS(1), CT_STAT)
          tot_reads = tot_reads + 1

	  IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	    call lib$signal(FDQ_NHKPREAD,%val(1),%val(nread_good))
	    if (ct_stat(1) .eq. ctp_endoffile) then
	      istat = eof_hsk

!	      Set function to return status

	      IF (RETSTAT.EQ.SUCCESS) THEN
	  	FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
	      ELSE
	  	FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
	      ENDIF

	      return
	    else
	      ISTAT = BAD
	      RETSTAT = ERROR
	      STATUS = CT_STAT(1)
	      CALL LIB$SIGNAL(FDQ_CTREADERR,%VAL(2),%VAL(STATUS),
	1	'NFS_HKP')

!	      Set function to return status

	      IF (RETSTAT.EQ.SUCCESS) THEN
	  	FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
	      ELSE
	  	FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
	      ENDIF

	      RETURN
	    endif
	  else
	    nread_good = nread_good + 1
	  ENDIF
	  IF (CT_STAT(11) .NE. CUR_CAT_ENTRY) THEN
	    call lib$signal(FDQ_NEWHKPSEG,%val(1),%val(ct_stat(11)))       
	    CUR_CAT_ENTRY = CT_STAT(11)
	  ENDIF
	  nread = nread + 1
	  NEW(1) = .TRUE.
	  NEW(2) = .FALSE.
	ENDIF

	BRKT_HSK_NOT_FOUND = .TRUE.

        if(tot_reads.eq.1 .and. istat.eq.brk_mj_miss .and.
     &              first_mj.eq.4)first_mj = 0       !No hsk for first science
                                                     !record. First hsk is NOT
                                                     !old. Try it for current
                                                     !science record.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			     !
!     Main DO WHILE loop     !
!			     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO WHILE (BRKT_HSK_NOT_FOUND)

	  FIRST_MJ = FIRST_MJ + 1

	  IF (FIRST_MJ .EQ. 3) THEN		! Previous good record is 2

	    NEW(1) = .FALSE.			! Mark first record as "old"
	    IF (.NOT. NEW(2)) THEN		! Second record is old,
						! so read HSK file into 2nd record
	      CALL CT_READ_ARCV (, HSK_LUN, HSK_RECS(2),CT_STAT)
              tot_reads = tot_reads + 1
	      IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	        call lib$signal(FDQ_NHKPREAD,%val(1),%val(nread_good))
!	        write (6, 6101) nread_good
	        if (ct_stat(1) .eq. ctp_endoffile) then
	          istat = eof_hsk

!		  Set function to return status

		  IF (RETSTAT.EQ.SUCCESS) THEN
	  		FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
		  ELSE
	  		FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
		  ENDIF

	          return
	        else
	          ISTAT = BAD
	          RETSTAT = ERROR
	          STATUS = CT_STAT(1)
	          CALL LIB$SIGNAL(FDQ_CTREADERR,%VAL(2),%VAL(STATUS),
	1       	'NFS_HKP')

!		  Set function to return status

		  IF (RETSTAT.EQ.SUCCESS) THEN
	  		FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
		  ELSE
	  		FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
		  ENDIF
	          
		  RETURN
	        endif
	      else
	        nread_good = nread_good + 1
	      ENDIF		! CT_STAT(1) .NE. CTP_NORMAL
	      IF (CT_STAT(11) .NE. CUR_CAT_ENTRY) THEN
		call lib$signal(FDQ_NEWHKPSEG,%val(1),%val(ct_stat(11)))  
	        CUR_CAT_ENTRY = CT_STAT(11)
	      ENDIF
	      nread = nread + 1
	      NEW(2) = .TRUE.

	    ENDIF		! .NOT. NEW(2)

	  ELSEIF (FIRST_MJ .EQ. 5) THEN

	    FIRST_MJ = 1			! HSK_RECS is a double circular buffer
	    NEW(2) = .FALSE.			! Mark 2nd record as old
	    IF (.NOT. NEW(1)) THEN		! If 1st record is old
						! then read HSKP into it
	      CALL CT_READ_ARCV (, HSK_LUN, HSK_RECS(1),CT_STAT)
              tot_reads = tot_reads + 1
	      IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	        call lib$signal(FDQ_NHKPREAD,%val(1),%val(nread_good))
!	        write (6, 6101) nread_good
	        if (ct_stat(1) .eq. ctp_endoffile) then
	          istat = eof_hsk

!		  Set function to return status

		  IF (RETSTAT.EQ.SUCCESS) THEN
	  		FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
		  ELSE
	  		FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
		  ENDIF

	          return
	        else
	          ISTAT = BAD
	          RETSTAT = ERROR
	          STATUS = CT_STAT(1)
	          CALL LIB$SIGNAL(FDQ_CTREADERR,%VAL(2),%VAL(STATUS),
	1	      'NFS_HKP')

!		  Set function to return status

		  IF (RETSTAT.EQ.SUCCESS) THEN
	  		FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
		  ELSE
	  		FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
		  ENDIF

	          RETURN
	        endif
	      else
	        nread_good = nread_good + 1
	      ENDIF		! CT_STAT(1) .NE. CTP_NORMAL
	      IF (CT_STAT(11) .NE. CUR_CAT_ENTRY) THEN
		call lib$signal(FDQ_NEWHKPSEG,%val(1),%val(ct_stat(11)))
!	        WRITE (6, 6401) CT_STAT(11)
	        CUR_CAT_ENTRY = CT_STAT(11)
	      ENDIF
	      nread = nread + 1
	      NEW(1) = .TRUE.

	    ENDIF		! .NOT. NEW(1)

	  ENDIF			! FIRST_MJ .EQ. 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									   !
!     Now compare time of HSKP mj fr pointed to by FIRST_MJ with	   !
!     Science time and declare done if Science time precedes HSKP time     !
!									   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(first_mj .eq. 1)then
            hsk_gmt = hsk_recs(1).ct_head.gmt
          elseif(first_mj .eq. 2)then
            hsk_gmt = hsk_recs(1).hskp_tail.gmt_mjf2
          elseif(first_mj .eq. 3)then
            hsk_gmt = hsk_recs(2).ct_head.gmt
          elseif(first_mj .eq. 4)then
            hsk_gmt = hsk_recs(2).hskp_tail.gmt_mjf2
          endif

	  CALL GMT_TO_BINARY (%ref(HSK_GMT), HSK_BNTM) 
	  IF (TIME_LT(SCI_BNTM, HSK_BNTM)) THEN
	    BRKT_HSK_NOT_FOUND = .FALSE.
	  ENDIF

	ENDDO			! (BRKT_HSK_NOT_FOUND) or main DO-WHILE loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								      !
!     We have just found the major frame which is just after the      !
!     Science time.  Therefore we should backtrack the pointer by     !
!     one to point to the 1st major frame which "brackets" the	      !
!     Science record.  Of course, backtracking from 1 means the	      !
!     previous major frame is in 4, as HSK_RECS is a circular	      !
!     double buffer.						      !
!								      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	NEXT_MJ = FIRST_MJ			! FIRST_MJ at this point is really the next major frame; remember it
	FIRST_MJ = FIRST_MJ - 1
	IF (FIRST_MJ .EQ. 0) THEN
	  FIRST_MJ = 4
	ENDIF

        if(first_mj .eq. 1)then
          gmt1 = hsk_recs(1).ct_head.gmt
        elseif(first_mj .eq. 2)then
          gmt1 = hsk_recs(1).hskp_tail.gmt_mjf2
        elseif(first_mj .eq. 3)then
          gmt1 = hsk_recs(2).ct_head.gmt
        elseif(first_mj .eq. 4)then
          gmt1 = hsk_recs(2).hskp_tail.gmt_mjf2
        endif

!
!     Now check if the 2 bracketing major frames are across record 
!     boundaries (will be so if FIRST_MJ=2 or FIRST_MJ=4); if so
!     check if the 2 major frames are more than 64 seconds apart;
!     (actually the tolerance should be a little over 32 seconds)
!     if more than 64 seconds apart then there is a gap, and we
!     really do not have 2 "bracketing" major frames so flag it
!     in ISTAT
!
!     Also, the test will also catch the situation when there is
!     no HSK record ahead of a Science record in the beginning of
!     the interval
!
	IF (FIRST_MJ .EQ. 2 .OR. FIRST_MJ .EQ. 4) THEN
	  CALL GMT_TO_BINARY (%ref(gmt1), HSK_BNTM1)	
!
!     HSK_BNTM1 is the time of the first major frame in binary; the second
!     major frame's binary time should still be in HSK_BNTM
!

	  STATUS = LIB$SUBX ( HSK_BNTM, HSK_BNTM1, DIFF )
	  IF ( STATUS .NE. SS$_NORMAL ) THEN
	    call lib$signal(FDQ_ERSUBX,%val(1),%val(status))
	    RETSTAT = ERROR
	  ENDIF
	  If (retstat.eq.success) then
	    IF (DIFF(1) .GT. SEC64 .or. diff(1) .lt. 0
	1	.or. diff(2) .gt. 0) THEN	! The 2nd test is needed because DIFF is
				! really a 64-bit integer and that DIFF(1) .lt. 0
				! means a very big integer!
				! 3rd test (checking DIFF(2)) added 5/20/87
	      ISTAT = BRK_MJ_MISS	
!	      Set function to return status

	      IF (RETSTAT.EQ.SUCCESS) THEN
	  	FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
	      ELSE
	 	FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
	      ENDIF

	      RETURN		
	    ENDIF	! Diff(1) .GT 64, etc.
	  endif		! Retstat is success
	ENDIF

!      Check return status before proceeding.
       if ( retstat .eq. success ) then
!
!     See if the FIRAS major frame ID in the status monitor is toggling
!     or not and set the sync flags STMON_SYNC
!
        if(first_mj .eq. 1)then       ! First major frame in buffer 1 of HSK_RECS
          bn(1) = 1
          bn(2) = 1
          mj(1) = 1
          mj(2) = 2
        elseif(first_mj .eq. 2)then   ! First major frame in buffer 2 of HSK_RECS
          bn(1) = 1
          bn(2) = 2
          mj(1) = 2
          mj(2) = 1
        elseif(first_mj .eq. 3)then   ! First major frame in buffer 3 of HSK_RECS
          bn(1) = 2
          bn(2) = 2
          mj(1) = 1
          mj(2) = 2
        else                          ! First major frame in buffer 4 of HSK_RECS
          bn(1) = 2
          bn(2) = 1
          mj(1) = 2
          mj(2) = 1
        endif
        stmon_a(1) = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(1)
        stmon_a(2) = hsk_recs(bn(2)).frame(mj(2)).hskp_head.stat_monitor_cmd(1)
        stmon_b(1) = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(5)
        stmon_b(2) = hsk_recs(bn(2)).frame(mj(2)).hskp_head.stat_monitor_cmd(5)

	if(ibits(stmon_a(1),14,1) .ne. ibits(stmon_a(2),14,1)) then	! Side A
	  stmon_sync(1) = 1		! Side A status monitor in sync
	else
	  stmon_sync(1) = 0		! Side A status monitor out of sync
	endif

	if (ibits(stmon_b(1),14,1) .ne. ibits(stmon_b(2),14,1)) then	! Side B
	  stmon_sync(2) = 1		! Side B status monitor in sync
	else
	  stmon_sync(2) = 0		! Side B status monitor out of sync
	endif
!
!     Move IPDU power relay statuses from FIRST_MJ (arbitrarily; can use the 2nd 
!     bracketing major frame as well.  This is done to save processing time)
!
	CALL LIB$MOVC3 (8, hsk_recs(bn(1)).frame(mj(1)).hskp_head.ipdu_stat,
     &                                                           IPDU_RELAY)

	ISTAT = GOOD
       endif 		! Retstat is success

!	Set function to return status
    
	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_BRK_HSK = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GET_BRK_HSK = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
