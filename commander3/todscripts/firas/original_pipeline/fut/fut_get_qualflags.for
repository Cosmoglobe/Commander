	INTEGER*4 FUNCTION FUT_GET_QUALFLAGS ( first_chan,
	1	END_SEGMENT, FIRST_MJ, SCI_REC, ENG_REC, HSK_RECS, 
 	1	SCILIM_REC, ENGLIM_REC, LIMFLAGS, LIMIT_ON, SC_SIDE, 
	1       PLOT, PLOT_TABLE, PLOT_START, PLOT_STOP, PLOT_DEVICE,
	1	MIN_OFFTIME, MAX_OFFTIME, REPORT_LUN, ENGLIM, SAMPRATE)

C------------------------------------------------------------------------------
C/
C/	PROGRAM NAME:
C/	  FUT_GET_QUALFLAGS
C/
C/	PROGRAM DESCRIPTION:
C/	  This program will determine the data quality of the Science record
C/	  based on the information in the Science record itself (being available
C/	  through named Record SCI_REC) and the housekeeping/engineering 
C/	  information (being available through named records ENG_REC and 
C/	  HSK_RECS) which brackets the Science record.
C/
C/	AUTHOR:
C/	  Shirley M. Read
C/	  STX
C/	  January 20, 1988
C/	  The entire list of the data quality flags was revised for Build 3.1.
C/	  Both the data items to be checked and the structure containing the 
C/	  limits for the data were modified to such an extent that the old
C/	  subroutine FUT_QUALFGS was unusable. The new FUT_QUALFLAGS was written
C/	  for Build 3.1 and will replace the old routine.
C/	
C/	MODIFIED BY:
C/	  Shirley M. Read
C/	  STX
C/	  March, 1988
C/	  After an analysis of the operation of the FIRAS instrument and an
C/	  examination of the circuit diagrams, a number of quality flags had
C/	  to be redefined and the limit enable/disable file will have to be 
C/	  redesigned. The redesign will be done for Build 4.0. As many of the 
C/	  data quality flags as possible were revised for Build 3.2. The 
C/	  requirements for many quality flags still need to be defined either 
C/	  as SPRs or enhancements for Build 4.0. 
C/
C/	  Shirley M.Read
C/	  STX
C/	  May, 1988
C/	  Modified the checking of the limits in accord with new citeria
C/	  developed as a result of analysis of the FIRAS instrument circuit
C/        diagram. Also included new limit checks based on decisions of the
C/	  FIRAS Team. 
C/	  Added the Limit_0n array, spacecraft side and report unit to the 
C/	  calling sequence.
C/	  
C/	  Shirley M. Read
C/	  STX
C/	  July, 1988
C/	  In discussions with the FIRAS instrument team about the commanding
C/	  of the instrument, we discovered that the calibrator resistors had
C/	  no limits defined in contrast to other analogs. This routine was
C/	  modified to mask out the limit checks on the calibrator resistors
C/	  for Build 4.1.
C/
C/	  H. Wang
C/	  STX
C/	  Jan., 1991
C/        Reason:     New requirements for FDQ
C/                    New requirements for reporting:
C/                     Banner announcement if red engineering
C/                     limits exceeded
C/
CH	CHANGE LOG:		New Format for Build 4.2 STX, 10/14/88
CH
CH      Version 4.1.1 10/15/88, SPR 2621, Shirley M. Read, STX
CH              Disable limit checking in FDQ for GRTs when in dwell mode.
CH      Version 4.1.1 10/15/88, SPR 2622, Shirley M. Read, STX
CH              Flag values for missing conversions need changes in FDQ.
CH	        The flag values for converted fields which are out of range
CH		need the same changes. GRTs in dwell mode need a flag value
CH              also. The SWG has selected a flag value of -9999 for all cases.
CH              The limit checking algorithms need to bypass setting the quality
CH              flags if the GRT or engineering analog has the flag value.
CH		Interpolation and averaging algorithms need to be modified to
CH              include the flag value. 
CH		Note: There are still some questions about setting the flags
CH		according to the FIRAS side which is powering the particular
CH	 	FIRAS instrument component. The post I and T analysis will
CH		determine the final algorithms. The analogs most likely to
CH		be affected are the IPDU temperatures, the drive box 
CH	 	temperatures and the MTM cal motor current. These analogs
CH		have been checked without regard to FIRAS A or B side since
CH		Build 3.2.
CH
CH	Version 4.2.1 02/09/89, SPR 3062, Shirley M. Read, STX
CH		Moon contaminated IFGs failed by data quality checking. FES and
CH		FCI will interpret the data as being unusable. It is needed for
CH		moon modeling. The moon flag bit will not be counted in the
CH		data quality summary flag.
CH	Version 4.2.1 02/09/89, SPR 2550, Shirley M. Read, STX
CH		Provide user option of setting the time interval for triggered
CH		plots.
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH	Version 4.2.1 03/05/89, SPR 3359, Shirley M. Read, STX
CH		FUT_Get_Qualflags has a data type error in the calling sequence
CH		for the logical unit number to write the report. The report
CH	        number is passed to the caller in an array of I*2 unit numbers.
Ch		It is declared as an I*4 in FUT_Get_Qualflags.
Ch
CH 	Version 4.5 09/07/89, SPR 4505, Shirley M. Read, STX
CH		FUT_Get_Qualflags has a tyo in the limit check of hot
CH		spot heaters. In a DO Loop the character 1 is used instead
CH		of ix.
CH
CH 	Version 5.4 1/25/90, SPR 5538, Harte Y. Wang, STX
CH		FUT_Get_Qualflags has to check bad telemetry from Housekeeping
CH	        Records and Science Records.	
CH
CH 	Version 5.8 3/7/90, SPR 6297, Harte Y. Wang, STX
CH		FUT_Get_Qualflags, If telemetry quality is bad, do not
CH              do the rest of limit checking.
CH
CH 	Version ?? 10/1/90, SPR 7495, Harte Y. Wang, STX
CH		FUT_Get_Qualflags, did not check the HKP telemetry quality 
CH              properly when the major fram pointer is 3 or 4.
CH
CH 	Version 7.1 10/31/90, SPR 7627, Harte Y. Wang, STX
CH              Checks on telemetry format.
CH
CH	Version 9.8 8/4/92, SER 9848, Steven Alexander, HSTX
CH		Add sampling rate as a passed parameter instead of a 
CH		hard-coded value.
C--------------------------------------------------------------------------
C/
C/	INPUT PARAMATERS:
C/
C/	  END_SEGMENT           L*1     Flag indicating segment ending
C/	  FIRST_MJ	        I*2     First major frame of HKP Rec set
C/	  SCI_REC(1536)  	BYTE	Buffer for current channel SCI record
C/	  ENG_REC(1024)         BYTE    Buffer for current Engineering record
C/	  HSK_RECS(576,2)       BYTE    Buffer for two bracketing HKP records
C/	  SCILIM_REC(512,2)     BYTE    Science/attitude limits: Red and Yellow
C/	  ENGLIM_REC(1024,4)    BYTE    Engineering/analog limits: Red Low,
C/					Yellow Low, Yellow High and Red High
C/	  LIMFLAGS(192)         BYTE    Limits On/Off flags
C/	  LIMIT_ON(110)         L*1     Limits enable/disable flags
C/	  SC_SIDE               I*2     Spacecraft side controlling FIRAS
C/	  PLOT		        I*2     Flag to generate plot. Values:
C/					0 = none, 1 = red limits only,
C/				        2 = yellow or red limits
C/	  PLOT_DEVICE           I*4     Plot device: line or laser printer
C/	  MIN_OFFTIME           I*4     Hours to offset plot start
C/	  MAX_OFFTIME           I*4     Hours to offset plot stop
C/	  REPORT_LUN            I*2     Unit number for report
C/	  ENGLIM		I*4	Per channel count of eng limit violations
C/	  SAMPRATE		R*4  	Sampling rate (I&T or Mission)
C/	
C/	OUTPUT PARAMETERS:
C/
C/	  SCI_REC(1534)         BYTE    Quality flags to be written to 
C/	    Record Field                Science record
C/	    DATA_QUALITY(110)		See INCLUDE file FUT_QUALFLAGS.TXT 
C/					for details
C/	  PLOT_TABLE(118)       L*1     Table of flags
C/	  PLOT_START(2)         I*4     Plot start time
C/	  PLOT_STOP(2)          I*4     Plot stop time
C/	
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/	
C/	INCLUDE FILES USED:
C/	  FUT_QUALFLAGS.TXT
C/	
C/	SUBROUTINES CALLED:
C/	  FUT_CHECKSUM
C/	  FUT_MJF_CHANGE
C/	  FUT_XTALK_ACCESS
C/	  FUT_XTALK_CHECK
C/	  FUT_SUM_FLG
C/	  FUT_QUAL_SUMMARY
C/        FUT_TRIGGER_PLOT
C/	
C/	ERROR HANDLING:
C/	  By calling Lib$Signal for error. Fut_Error has been established
C/	  to handle all errors and return control to the FDQ program.
C/	
C/	METHOD USED:
C/
!	FUT_GET_QUALFLAGS determines the quality of the data associated 
!	with 110 Data_Quality_Flags and sets these flags for each input
!	Science_Record. Whenever applicable, the red and yellow limits
!	are used for comparison. An overall summary flag is computed by
!	counting the red and yellow flags. The function is invoked for
!	each active science channel.
!
!	PDL for FUT_Get_Qualflags
!
!	BEGIN
!
!	C  For purposes of organizing the limit checking routines, the
!	C  Data_Quality_Flags have been grouped according to the type of
!	C  checking to be performed. Most of the flags are set to 0 if
!	C  the quality is good and 1 if the quality is bad.
!
!	C  Group 1 : FIRAS SCIENCE-TELEMETRY DATA QUALITY
!
!	IF ( LIMIT_ON for FLG_BADSCI ) THEN
!	  DO for Science_Record Data_Ready flags 1 to 8
!	    IF ( Science_Record Data_Ready is not equal 0 ) THEN
!	       SET FLG_BADSCI to 1
!              Banner announcement red limits exceeded 
! 	    ENDIF
!	  ENDDO
!	ENDIF
!	IF ( LIMIT_ON for FLG_BADHKP ) THEN
!	  DO for Housekeeping_Record Telemetry_Quality_Flags 1 to 2
!	    IF ( Housekeeping_Record Telemetry_Quality_Flags is not equal 0 )
!	      THEN
!	      SET FLG_BADHKP to 1
!              Banner announcement red limits exceeded 
!	    ENDIF
!	  ENDDO
!	ENDIF
!	IF ( LIMIT_ON for FLG_CAL ) THEN
!	  IF ( Science_Record External_Calibrator_Position is equal 0 or 3) THEN
!	    SET FLG_CAL to 1
!           Banner announcement red limits exceeded 
!	  ENDIF
!	ENDIF
!
!	C  GROUP 2 : FIRAS SCIENCE MICROPROCESSOR HEADER WORDS
!
!	  MASK out Bit 12 in the Science_Record Micro_Header Status_Word
!         AND the result with the two byte limits flags bit mask
!	  IF ( LIMIT_ON for FLG_MICRO ) THEN
!	    IF ( the masked Status_Word is not equal to 0 ) THEN
!	      SET FLG_MICRO to 1
!	      AND the result with the Limflags bit mask
!             Banner announcement red limits exceeded 
!	    ENDIF
!	  ENDIF 
!
!	IF ( LIMIT_ON for FLG_GLITCH_CT ) THEN
!	  IF ( Science_Record Micro_Header Deglitch_Overflow_Address is
!	    filled ) THEN
!	    SET FLG_GLITCH_CT to 1
!	  ELSEIF ( the Science_Record Micro_Header Glitch_Count > FEX_SCILIM 
!	    database limits Glitch_Count ) THEN
!	    SET FLG_GLITCH_CT to 1
!           Banner announcement red limits exceeded 
!	  ENDIF 
!	ENDIF
!	IF ( LIMIT_ON for FLG_SATURATE ) THEN
!	  IF ( Science_Record Micro_Header Saturated_Sample_Count > FEX_SCILIM
!	    database limits Saturated_Sample_Count ) THEN
!	    SET FLG_SATURATE to 1
!           Banner announcement red limits exceeded 
!	  ENDIF
!	ENDIF
!	IF ( LIMIT_ON for FLG_CKSM_ERR_ST ) THEN
!
!	  CALL FUT_CHECKSUM to compute the data checksum from the interferogram
!	     points and determine the value for FLG_CKSM_ERR_ST
!            Banner announcement red limits exceeded 
!	  
!	ENDIF
!
!	C  GROUP 3 : FIRAS ATTITUDE LIMITS
!
!	C  Set the bits for the red and yellow attitude limits flags   
!
!	IF ( LIMIT_ON for FLG ATT_ALT_RED ) THEN
!	  IF ( Science_Record Sun_Angle is outside FEX_SCILIM database red
!	    limits Sun_Angle ) THEN
!
!	    SET Bit 0 in FLG_ATT_ALT_RED
!
!	  ENDIF
!	  IF ( Science_Record Earth_Limb is outside FEX_SCILIM database red
!	    limits Earth_Limb ) THEN
!
!	    SET Bit 1 in FLG_ATT_ALT_RED
!
!	  ENDIF
!	  IF ( Science_Record Moon_Angle is outside FEX_SCILIM database red
!	    limits Moon_Angle ) THEN
!
!	    SET Bit 2 in FLG_ATT_ALT_RED
!
!	  ENDIF
!
!	IF ( LIMIT_ON for FLG_ATT_ALT_YEL ) THEN
!	  IF ( Science_Record Sun_Angle is outside FEX_SCILIM database yellow
!	    limits Sun_Angle ) THEN
!
!	    SET Bit 0 in FLG_ATT_ALT_YEL
!
!	  ENDIF
!	  IF ( Science_Record Earth_Limb is outside FEX_SCILIM database yellow
!	    limits Earth_Limb ) THEN
!
!	    SET Bit 1 in FLG_ATT_ALT_YEL
!
!	  ENDIF
!	  IF ( Science_Record Moon_Angle is outside FEX_SCILIM database yellow
!	    limits Moon_Angle ) THEN
!
!	    SET Bit 2 in FLG_ATT_ALT_YEL
!
!	  ENDIF
!
!
!
!	  ENDIF
!	  IF ( Science_Record No attitude found  
!	     ) THEN
!
!	    SET Flg_no_attitude 1
!
!	  ENDIF
!
!	C  GROUP 4 : GRT LIMITS
!
!	  DO for each GRT on A side
!	    IF ( LIMIT_ON for FLG_GRTRED_ST ) THEN
!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                red low limits for GRT OR > FEX_ENGLIM database
!		 red high limits for GRT ) THEN
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTRED_ST flags 
!               Banner announcement red limits exceeded 
!	      ENDIF 	
!	    ENDIF
!	    IF ( LIMIT_ON for FLG_GRTYEL_ST ) THEN
!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                yellow low limits for GRT OR > FEX_ENGLIM database
!		 yellow high limits for GRT ) THEN
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTYEL_ST flags
!	      ENDIF 	
!	    ENDIF
!	  ENDDO		    
!	  DO for each GRT on B side
!	    IF ( LIMIT_ON for FLG_GRTRED_ST ) THEN
!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                red low limits for GRT OR > FEX_ENGLIM database
!		 red high limits for GRT ) THEN
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTRED_ST flags
!                Banner announcement red limits exceeded 
!	      ENDIF 	
!	    ENDIF
!	    IF ( LIMIT_ON for FLG_GRTYEL_ST ) THEN
!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                yellow low limits for GRT OR > FEX_ENGLIM database
!		 yellow high limits for GRT ) THEN
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTYEL_ST flags
!	      ENDIF 	
!	    ENDIF
! 	  ENDDO		    
!	  AND each GRT flag with the bit mask in Limflags and mask out
!	      the calibrator resistors.
!
!	C  GROUP 5 : ENGINEERING LIMITS
!
!	C  Check all the channel specific words in the group and preamps which
!	C  depend on the spacecraft side controlling FIRAS. The flag setting
!       C  algoriths apply to the following set of data quality flags:
!	C  FLG_CHAN_TMP_ST, FLG_OPPA_TMP, FLG_CHPA_TMP and FLG_BOL_VOL_ST.
!
!	DO for all channel engineering flags 
!         IF ( LIMIT_ON for channel engineering flag ) THEN
!	    IF ( Engineering_Record Channel_Engineering_Word < FEX_ENGLIM
!	      database yellow low limits for Channel_Engineering_Word or
!	      > FEX_ENGLIM database yellow high limits for Channel_Engineering_
!             Word ) THEN
!	      SET channel engineering flag to 1
!	    ENDIF
!	    IF ( Engineering_Record Channel_Engineering_Word < FEX_ENGLIM
!	      database red low limits for Channel_Engineering_Word or
!	      > FEX_ENGLIM database red high limits for Channel_Engineering_
!             Word ) THEN
!	      SET channel engineering flag to 2
!             Banner announcement red limits exceeded 
!	    ENDIF
!	  ENDIF	   
!	ENDDO
!       DO for the Preamps
!         IF ( LIMIT_ON for Preamp on SC_Side ) THEN
!	    IF ( Engineering_Record Preamp_Engineering_Word < FEX_ENGLIM
!	      database yellow low limits for Preamp_Engineering_Word or
!	      > FEX_ENGLIM database yellow high limits for Preamp_Engineering_
!             Word ) THEN
!	      SET Preamp engineering flag to 1
!	    ENDIF
!	    IF ( Engineering_Record Preamp_Engineering_Word < FEX_ENGLIM
!	      database red low limits for Preamp_Engineering_Word or
!	      > FEX_ENGLIM database red high limits for Preamp_Engineering_
!             Word ) THEN
!	      SET Preamp engineering flag to 2
!             Banner announcement red limits exceeded 
!	    ENDIF
!	  ENDIF	   
!	ENDDO
!
!	C  The remainder of the engineering limits flags refer to engineering
!       C  words which have an A and B side. The side read may depend on the
!	C  spacecraft side controlling FIRAS. Some flags have two bytes. Byte 
!       C  one may refer to the active side and byte two for the other side,
!       C  or the bytes may depend on which side of the spacecraft is powering
!       C  FIRAS or the two bytes may refer to the FIRAS A and B sides. 
!	C  The next algorithm applies to all of the flags in this set.
!
!	  DO for all engineering words in the group
!	    IF ( LIMIT_ON for the engineering component ) THEN
!	      IF ( Engineering_Record Engineering_Word < FEX_ENGLIM
!	        database yellow low limits for Engineering_Word OR > 
!	        FEX_ENGLIM database yellow high limits for Engineering_Word )
!               THEN
!	        SET engineering field flag to 1
!	      ENDIF
!	      IF ( Engineering_Record Engineering_Word < FEX_ENGLIM
!	        database red low limits for Engineering_Word OR >
!	        FEX_ENGLIM database red high limits for Engineering_Word )
!               THEN
!	        SET engineering field flag to 2
!               Banner announcement red limits exceeded 
!	      ENDIF
!	    ENDIF	   
!	  ENDDO
!
!	C  GROUP 6 : MAJOR FRAME BOUNDARY CHANGES
!
!	C  Status changes across the Major Frame boundary, indicated by a 
!	C  change in state of a particular bit in the arrays of the status
!       C  monitor command in the two bracketing housekeeping records will
!	C  be used to set FLG_STCHG_MJ_ST. 
!
!	CALL FUT_MJF_CHANGE to extract the relevant bits from the status
!	     monitor command arrays in the two bracketing housekeeping
!	     records, check for a change in state and set the corresponding
!	     bit in the FLG_STCHG_MJ_ST to 1. AND with the bit masks in 
!	     Limflags.
!
!	C  GROUP 7 : STATUS FROM DIRBE, DMR or SPACECRAFT AFFECTING FIRAS
!
!	CALL FUT_Xtalk_access to access DIRBE,DMR,SPACECRAFT
!           states that will affect FIRAS 
!	CALL FUT_XTALK_CHECK to check DIRBE,DMR SPACECRAFT states that will affect FIRAS 
!
!	C  GROUP 8 : FIRAS DATA QUALITY SUMMARY FLAG
!
!	CALL FUT_SUM_FLG to count the red and yellow flags set for the
!	     entire Science_Record and determine a value for the overall
!	     summary quality flag, FLG_SUMMARY. The Moon_Angle flag will
!	     not be included in the count. 
!
!	C  A summary table will be kept for all individual flags set by
!	C  the routine for the entire segment. For each IFG the appropriate
!	C  data quality flag counters will be incremented is the flag indicates
!	C  bad quality. When the segment is completer, a summary will be sent
!	C  to the user via sys$out device in effect for the run. A message for
!	C  each quality flag will be included in FUT_MSG.MSG.
!
!	CALL FUT_QUAL_SUMMARY to increment the counters for the flags set for
!	     the channel Science Record and output the summary at the end of
!	     the segment.
!
!       IF the PLOT flag is set to red or red and yellow limits,
!         CALL FUT_TRIGGER_PLOT to set flags in a table to be used to write
!            a command file which will trigger Eng_Plots of all engineering
!            fields for which red or yellow limits were exceeded.
!	ENDIF
!
!	SET return status to success 
!	
!	RETURN
!       END
!
C-------------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'(FUT_QUALFLAGS)'

	EXTERNAL   	FUT_NORMAL
	EXTERNAL	FUT_ABERR
	EXTERNAL	FUT_CHECKSUMER
	EXTERNAL	FUT_MJFCHANGER
	EXTERNAL	FUT_SUMFLAGER
	EXTERNAL	FUT_QUALSUMER
	EXTERNAL	FUT_TRIGPLOTERR

!	Functions called

	INTEGER*4	FUT_Checksum
	INTEGER*4	FUT_Mjf_Change
	INTEGER*4	FUT_Sum_Flg
	INTEGER*4	FUT_Qual_Summary
	INTEGER*4       FUT_Trigger_Plot

!	The following functions will be added for a later build.	
	INTEGER*4	FUT_XTALK_ACCESS
	INTEGER*4	FUT_XTALK_CHECK
        BYTE            DUMMY1(10),DUMMY2(10),DUMMY3(10),DUMMY4(10)

!	Passed parameters
        Logical*1       First_chan      ! First channel flag 
  	logical*1	END_SEGMENT     ! End of segment flag
	integer*2       FIRST_MJ	! First major frame in HKP records
	logical*1       LIMIT_ON(110)   ! Limits enable/disable flags
	logical*1       LIMIT_ON_GRTRED(16)   ! GRT Limits enable/disable flags
	logical*1       LIMIT_ON_GRTYEL(16)   ! GRT Limits enable/disable flags
	integer*2	GRTFLAGS
	integer*2       SC_SIDE         ! Spacecraft side 1 or 2
	integer*2	PLOT		! Plot flag
        logical*1       Telm_flag       ! telemetry bad flag
	logical*1       PLOT_TABLE(118) ! Table of flags to trigger plots
	integer*4       PLOT_START(2)   ! Plot start time
	integer*4       PLOT_STOP(2)    ! Plot stop time
	integer*4       PLOT_DEVICE     ! Plot device: line or laser printer
	integer*4       MIN_OFFTIME     ! Hours to offset plot start
	integer*4       MAX_OFFTIME     ! Hours to offset plot stop
	integer*2       REPORT_LUN      ! Unit to write report
	integer*4	ENGLIM    	! engineering limits violation counts
	real*4		SAMPRATE	! sampling rate (I&T or Mission)
        dictionary 'nfs_sdf'
	record /nfs_sdf/  SCI_REC

        dictionary 'fdq_eng'
	record /fdq_eng/  ENG_REC

        dictionary 'nfs_hkp'
	record /nfs_hkp/  HSK_RECS(2)

        dictionary 'fex_scilim'
	record /fex_scilim/  SCILIM_REC

        dictionary 'fex_englim'
	record /fex_englim/  ENGLIM_REC

        dictionary 'fex_limflags'
	record /fex_limflags/  LIMFLAGS


!	Local variables

        Logical*1       First_redlim    
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS         ! Values for status
	parameter       ( SUCCESS = 1 )
        integer*4       ERROR
	parameter	( ERROR = 2 ) 
	integer*4	STATUS		! Dummy status variable
	integer*4	ZEROSTAT 
	parameter       ( ZEROSTAT = 0 )

	integer*2	CHAN, DIFF, Ix, Jx, side, curside, oppside
!	integer*2       all_set_mask
!	parameter       ( all_set_mask = 'FFFF'x )
!	integer*2       last_word_mask
!	parameter       ( last_word_mask = 3 )
	byte	        zero 
	parameter       ( zero =  0 )
	byte	        one 
	parameter       ( one =  1 )
	byte	        two 
	parameter       ( two =  2 )
	integer*2 	i2zero
	parameter	( i2zero = 0 )
	integer*2       i2one 
	parameter       ( i2one = 1 )
	integer*2       i2two
	parameter       ( i2two = 2 )
	integer*2       status_mask 
	parameter       ( status_mask = 'EFFF'x )
	integer*2       zeromask 
	parameter  	( zeromask = 0 )
	integer*2       masked_status_word
	byte		byte_status(2)
	equivalence     (masked_status_word, byte_status(1))
	integer*2       flag
	logical*1	proceed
	logical*1	skip
	byte	        Gmask(16) / '01'x, '02'x, '04'x, '08'x,
	1			    '10'x, '20'x, '40'x, '80'x,
	2	                    '01'x, '02'x, '04'x, '08'x,
	3			    '10'x, '20'x, '40'x, '80'x /

	integer*2       OEmask(10,2) 
	1		/ 1, 3, 5, 7, 9,  11, 13, 15, 17, 19,  
	2		  2, 4, 6, 8, 10, 12, 14, 16, 18, 20 / 
	
	byte		Calmask / 'C3'x /         ! Mask for cal resistors

	integer*2 	i2mask			  ! I*2 storage for masks
	byte		bytemask(2)               ! Byte storage for masks
	equivalence     (i2mask, bytemask(1))
	integer*2 	i2val			  ! I*2 storage for byte values
	byte		byteval(2)                ! Byte storage for byte values
	equivalence     (i2val, byteval(1))

	integer*2	Red / 1 /, Yellow / 2 /   ! Plot types
	logical*1	Statmon_Flag (2)	  ! Flag A and B sides for
						  ! exceeding stat. mon. limits
	logical*1	Volts_Flag (2,10)	  ! Flag IPDU component voltages
	logical*1	Cur_Flag   (2,6)          ! Flag IPDU component currents

	integer * 4	d_ready		!data ready dropout counter
        integer * 4	j,k             !counters
	integer * 4	start		!start counter
	integer * 4	stop		!stop counter
	real * 4        FV / -9999.0 /,gli_rate
	byte		All_Set / 'FF'x /  ! All bits set ( -1 )
        Logical*1       No_attitude
        Character*40    P_NAME
        CHARACTER*2     CHAN_NAME(4)/'RH','RL','LH','LL'/
        CHARACTER*8     SIDE1  
	character*15 GRT_Name(64) /  
	1 '  Xcal_Temp_Al ','  SkyHorn_Tp_Al','  RefHorn_Tp_Al',
	2 '  Ical_Temp_Al ','  Dihedrl_Tp_Al','  Bol_A_RH_T_Al',
	3 '  Bol_A_RL_T_Al','  Bol_A_LH_T_Al','  Bol_A_LL_T_Al',
	4 '  Mirror_Tmp_Al','  Calrs_RH_T_Al','  Calrs_RL_T_Al',
	5 '  Calrs_LH_T_Al','  Calrs_LL_T_Al','  Xcal_Tp_S5_Al',
	6 '  Colimatr_T_Al',  			
	1 '  Xcal_Temp_Ah ','  SkyHorn_Tp_Ah','  RefHorn_Tp_Ah', 
	2 '  Ical_Temp_Ah ','  Dihedrl_Tp_Ah','  Bol_A_RH_T_Ah',
	3 '  Bol_A_RL_T_Ah','  Bol_A_LH_T_Ah','  Bol_A_LL_T_Ah',
	4 '  Mirror_Tmp_Ah','  Calrs_RH_T_Ah','  Calrs_RL_T_Ah',
	5 '  Calrs_LH_T_Ah','  Calrs_LL_T_Ah','  Xcal_Tp_S5_Ah',
	6 '  Colimatr_T_Ah',   			
	1 '  Xcal_Temp_Bl ','  SkyHorn_Tp_Bl','  RefHorn_Tp_Bl',
	2 '  Ical_Temp_Bl ','  Dihedrl_Tp_Bl','  Bol_A_RH_T_Bl',
	3 '  Bol_A_RL_T_Bl','  Bol_A_LH_T_Bl','  Bol_A_LL_T_Bl',
	4 '  Mirror_Tmp_Bl','  Calrs_RH_T_Bl','  Calrs_RL_T_Bl',
	5 '  Calrs_LH_T_Bl','  Calrs_LL_T_Bl','  Xcal_Tp_S6_Bl',
	6 '  Colimatr_T_Bl',   			
	1 '  Xcal_Temp_Bh ','  SkyHorn_Tp_Bh','  RefHorn_Tp_Bh',
	2 '  Ical_Temp_Bh ','  Dihedrl_Tp_Bh','  Bol_A_RH_T_Bh',
	3 '  Bol_A_RL_T_Bh','  Bol_A_LH_T_Bh','  Bol_A_LL_T_Bh',
	4 '  Mirror_Tmp_Bh','  Calrs_RH_T_Bh','  Calrs_RL_T_Bh',
	5 '  Calrs_LH_T_Bh','  Calrs_LL_T_Bh','  Xcal_Tp_S5_Bh',
	6 '  Colimatr_T_Bh' /  			


	character*15 Ipdu_Volt_Name(10) /        ! IPDU Voltages
	1 '  Dig_Cnvt_N15v','  Dig_Cnvt_P15v','  Dig_Cnvt_P5v ',
	2 '  Anl_Cnvt_P15v','  Anl_Cnvt_N15v','  Bias_Pre_P25v',
	3 '  Int_Ps_P28v  ','  Int_Ps_P15v  ','  Int_Ps_N15v  ',
	4 '  Int_Ps_P5v   ' /

	character*15 Ipdu_Cur_Name(6) /		 ! IPDU Currents.
	1 '  Bias_Pre_Reg ','  Analog_Conv  ','  Digital_Conv ',
	2 '  Con_Current_H','  Con_Current_L','  Con_Int_Conv ' /

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !  
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.
        First_redlim = .false.
	retstat = success
        No_attitude = .false.
C	Initialize Science_Record Data_Quality flags to zero.

	DO Ix = 1, 110
	   Sci_Rec.DQ_Data.Data_Quality(Ix) = Zero
	ENDDO

	chan = Sci_Rec.Sci_Head.Chan_Id
					
!	C  For purposes of organizing the limit checking routines, the
!	C  Data_Quality_Flags have been grouped according to the type of
!	C  checking to be performed. Most of the flags are set to 0 if
!	C  the quality is good and 1 if the quality is bad.
!
!	C  Group 1 : FIRAS SCIENCE-TELEMETRY DATA QUALITY
!

	IF ( LIMIT_ON ( FLG_BADSCI) ) THEN
!	  If ( Science_Record Data Ready Count > FEX_SCILIM
!	    database limits Data_Ready Count) set the quality flag.
 	  d_ready = 0
	  do j = 1,6
            if(j.eq.1)then
              start = 6
              stop = 15
            else if(j.eq.6)then
              start = 0
              stop = 13
            else
              start = 0
              stop = 15
            end if
            do k = start,stop
              if(btest(sci_rec.sci_head.data_ready(j),k)
     .	                                  .eq. 0)d_ready = d_ready + 1
            end do
          end do
	  IF ( d_ready .GT. 
	1	Scilim_Rec.Lim(1).Sci_Head.Data_Ready(1)) THEN
	      Sci_Rec.DQ_Data.Data_Quality(FLG_BADSCI) = one
C              If (.not. first_redlim) then
C                Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                First_redlim=.true.
C              Endif
C              P_NAME='Missing data ready bits'
C              Write(report_lun,899)P_NAME
	  ENDIF 
	ENDIF

	
	IF ( LIMIT_ON (FLG_BADHKP) ) THEN
          Telm_flag = .false.
          Do ix=1, 60
           If(Sci_Rec.Sci_Head.Data_Qual(ix) .ne. 0) Telm_flag = .true.
          Enddo 
	  IF ( first_mj .eq. i2one ) then

	    DO Ix = 1, 2
	      IF (( Hsk_Recs(1).Mj_Frm(Ix).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_flag)   
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	    ENDDO
	    IF (( Hsk_Recs(1).ct_head.hskp1_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1   Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	    IF (( Hsk_Recs(1).ct_head.hskp2_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1   Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one

	  ELSEIF ( first_mj .eq. i2two ) THEN

	      IF ( (Hsk_Recs(1).Mj_Frm(2).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_flag)   
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	      IF ( (Hsk_Recs(2).Mj_Frm(1).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_flag)   
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
              IF (( Hsk_Recs(1).ct_head.hskp2_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1     Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	      IF (( Hsk_Recs(2).ct_head.hskp1_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1     Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one

	  ELSEIF ( first_mj .eq. 3 ) THEN
              	 

	    DO Ix = 1, 2
	      IF ( (Hsk_Recs(2).Mj_Frm(Ix).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_Flag)        
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	    ENDDO
	    IF (( Hsk_Recs(2).ct_head.hskp1_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1   Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	    IF (( Hsk_Recs(2).ct_head.hskp2_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1   Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one

	  ELSEIF ( first_mj .eq. 4 ) THEN
	 

	    
	      IF ( (Hsk_Recs(2).Mj_Frm(2).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_Flag)        
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	      IF ( (Hsk_Recs(1).Mj_Frm(1).Tlm_Qual_Maj_Frm .ne. zero )
	1	.or. Telm_Flag)        
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	      IF (( Hsk_Recs(2).ct_head.hskp2_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1   Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one
	      IF (( Hsk_Recs(1).ct_head.hskp1_tlm_fmt .eq. 0 )
	1	.or. Telm_flag)   
	1     Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP) = one

	  ENDIF
	ENDIF
C	IF ((Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP).eq. 1)) then
C          If (.not. first_redlim) then
C            Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C            First_redlim=.true.
C          Endif
C          P_NAME='Telemetry quality bad'
C          Write(report_lun,899)P_NAME
C        Endif
	IF (LIMIT_ON(FLG_CAL).and.(Sci_Rec.DQ_Data.Data_Quality
	1       (FLG_BADHKP).eq. 0)) then

	  IF (( Sci_Rec.DQ_Data.Xcal_Pos .ne. i2one ) .and.
	1     ( Sci_Rec.DQ_Data.Xcal_Pos .ne. i2two )) THEN

	      Sci_Rec.DQ_Data.Data_Quality(FLG_CAL) = one
C              If (.not. first_redlim) then
C                Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                First_redlim=.true.
C              Endif
C              P_NAME=' Sci or Cal. data Travel error '
C              Write(report_lun,899)P_NAME

	  ENDIF
	ENDIF
!
!	C  GROUP 2 : FIRAS SCIENCE MICROPROCESSOR HEADER WORDS
!
	IF ((LIMIT_ON(FLG_MICRO)).and.(Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP)
	1                     .eq. 0)) THEN
!	  MASK out Bit 12 in the Science_Record Micro_Header Status_Word.
!         AND with Limflags FLG_MICRO mask.

	  Masked_Status_Word = Sci_Rec.Sci_Head.SC_Head3 .AND.
	1		       Status_Mask
	  IF ( Masked_Status_Word .ne. zeromask ) THEN
	    if ( byte_status(1) .ne. zero ) then
	      byteval(1) = byte_status(1)
	      bytemask(1) = Limflags.Lim_Flags.FLG_MICRO(1)
	      i2val = i2val .AND. i2mask
	      Sci_Rec.DQ_Data.Data_Quality(FLG_MICRO) = byteval(1) 
	    endif
 	    if ( byte_status(2) .ne. zero ) then
	      byteval(1) = byte_status(2)
	      bytemask(1) = Limflags.Lim_Flags.FLG_MICRO(2)
	      i2val = i2val .AND. i2mask
	      Sci_Rec.DQ_Data.Data_Quality(FLG_MICRO + 1) = byteval(1) 
	    endif
C            If (byte_status(1) .ne. zero .or. byte_status(2) .ne. zero)
C	1     then
C                If (.not. first_redlim ) then
C                  Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                  First_redlim=.true.
C                Endif
C                P_NAME=' microprocessor status word '
C                Write(report_lun,899)P_NAME
C              endif
	  ENDIF 
	ENDIF

	IF(LIMIT_ON(FLG_GLITCH_CT).and.(Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP)
	1                     .eq. 0)) THEN

!	  If ( Science_Record Micro_Header Deglitch_Overflow_Address is
!	    filled ) set the quality flag.
	  IF ( Sci_Rec.Sci_Head.SC_Head22 .NE. i2zero ) THEN 
	      Sci_Rec.DQ_Data.Data_Quality(FLG_GLITCH_CT) = one

!	  Elseif ( the Science_Record Micro_Header Glitch_Count > FEX_SCILIM 
!	    database limits Glitch_Count ) set the quality flag.
	  ELSE
           If (sci_rec.sci_head.sc_head10 .ne. 0) then
            Gli_rate = sci_rec.sci_head.sc_head21/((sci_rec.sci_head.sc_head9*
	1	sci_rec.sci_head.sc_head11*sci_rec.sci_head.sc_head10) /
	2	samprate)
            If (sci_rec.sci_head.chan_id .eq. 1) then
              IF ( gli_rate .GT. 
	1	Scilim_Rec.Lim(1).Glitch_rate(1) ) THEN
	        Sci_Rec.DQ_Data.Data_Quality(FLG_GLITCH_CT) = one
              Endif
            ElseIf (sci_rec.sci_head.chan_id .eq. 2) then
              IF ( gli_rate .GT. 
	1	Scilim_Rec.Lim(1).Glitch_rate(2) ) THEN
	        Sci_Rec.DQ_Data.Data_Quality(FLG_GLITCH_CT) = one
              Endif
            ElseIf (sci_rec.sci_head.chan_id .eq. 3) then
              IF ( gli_rate .GT. 
	1	Scilim_Rec.Lim(1).Glitch_rate(3) ) THEN
	        Sci_Rec.DQ_Data.Data_Quality(FLG_GLITCH_CT) = one
              Endif
            ElseIf (sci_rec.sci_head.chan_id .eq. 4) then
              IF ( gli_rate .GT. 
	1	Scilim_Rec.Lim(1).Glitch_rate(4) ) THEN
	        Sci_Rec.DQ_Data.Data_Quality(FLG_GLITCH_CT) = one
              Endif
C              If (.not. first_redlim) then
C                Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                First_redlim=.true.
C              Endif
C              P_NAME=' Glitch count too high '
C              Write(report_lun,899)P_NAME
            Endif 
           Endif
	  ENDIF 
	ENDIF

	IF (LIMIT_ON(FLG_SATURATE).and.(Sci_Rec.DQ_Data.Data_Quality(FLG_BADHKP)
	1                     .eq. 0)) THEN 
!	  If ( Science_Record Micro_Header Saturated_Sample_Count > FEX_SCILIM
!	    database limits Saturated_Sample_Count ) set the quality flag.
	  IF ( Sci_Rec.Sci_Head.SC_Head20 .GT. 
	1	Scilim_Rec.Lim(1).Sci_Head.SC_Head20 ) THEN
	      Sci_Rec.DQ_Data.Data_Quality(FLG_SATURATE) = one
C              If (.not. first_redlim) then
C                Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                First_redlim=.true.
C              Endif
C              P_NAME='Saturated sample count high'
C              Write(report_lun,899)P_NAME
	  ENDIF 
	ENDIF

	IF(LIMIT_ON(FLG_CKSM_ERR_ST).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)) THEN
!
!	CALL FUT_CHECKSUM to compute the data checksum from the interferogram
!	     points and determine the value for FLG_CKSM_ERR_ST
!	  
	  STATUS = FUT_CHECKSUM ( SCI_REC,
	1	Sci_Rec.DQ_Data.Data_Quality(FLG_CKSM_ERR_ST))

	  if ( status .ne. %loc(FUT_NORMAL) ) then
	    if ( status .eq. %loc(FUT_ABERR)) status = zerostat
	    call lib$signal(FUT_CHECKSUMER,%val(1),%val(status))
	    retstat = error
C	  else
C           If (sci_rec.dq_data.data_quality(flg_cksm_err_st) .eq. 1) then
C            If (.not. first_redlim) then
C              Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C              First_redlim=.true.
C            Endif
C            P_NAME=' Checksum errors '
C            Write(report_lun,899)P_NAME
C           Endif
          endif  
	ENDIF
!
!	C  GROUP 3 : FIRAS ATTITUDE LIMITS
!
!	C  Set the bits for the red and yellow attitude limits flags   
!
        If (LIMIT_ON(FLG_NO_ATTITUDE).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)
	1	.AND.RETSTAT.EQ.SUCCESS ) then
          IF (sci_rec.attitude.solution .eq. 0) then
            sci_rec.dq_data.data_quality(FLG_NO_ATTITUDE) = 1
            no_attitude = .true.
          endif
        Endif   
	IF ( LIMIT_ON(FLG_ATT_ALT_RED).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0) .and. (.not. no_attitude)
	1	.AND.RETSTAT.EQ.SUCCESS ) THEN 

!	  If ( Science_Record Sun_Angle is outside FEX_SCILIM database red
!	    limits Sun_Angle ), SET Bit 0 in FLG_ATT_ALT_RED
!
	  IF ( Sci_Rec.Attitude.Sun_Angle .LT. 
	1      Scilim_Rec.Lim(1).Attitude.Sun_Angle ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED)
	    bytemask(1) = Gmask(1)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED) = byteval(1) 
	  ENDIF 

!	  If ( Science_Record Earth_Limb is outside FEX_SCILIM database red
!	    limits Earth_Limb ), SET Bit 1 in FLG_ATT_ALT_RED
	  IF ( Sci_Rec.Attitude.Earth_Limb .LT. 
	1      Scilim_Rec.Lim(1).Attitude.Earth_Limb ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED)
	    bytemask(1) = Gmask(2)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED) = byteval(1) 
	  ENDIF

!	  IF ( Science_Record Moon_Angle is outside FEX_SCILIM database red
!	    limits Moon_Angle ), SET Bit 2 in FLG_ATT_ALT_RED
	  IF ( Sci_Rec.Attitude.Moon_Angle .LT. 
	1      Scilim_Rec.Lim(1).Attitude.Moon_Angle ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED)
	    bytemask(1) = Gmask(3)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED) = byteval(1) 
	  ENDIF

!    The following routines will be added in a later build.
!	  CALL FUT_CHECK_SAA to determine whether Science_Record COBE latitude
!	  and longitude is within FEX_SCILIM database red limits COBE 
!	  latitude and longitude for the South Atlantic Anomaly and set 
!	  Bit 3 in FLG_ATT_ALT_RED
!
!	  CALL FUT_CHECK_VAB to determine whether Science_Record COBE 
!	  latitude and longitude is within FEX_SCILIM database red limits
!	  COBE latitude and longitude for the Van Allen Belts and set 
!	  Bit 4 in FLG_ATT_ALT_RED
!
!	AND with Limflags FLG_ATT_ALT_RED bit mask.

	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED)
	    bytemask(1) = Limflags.Lim_Flags.FLG_ATT_ALT_RED
	    i2val = i2val .AND. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_RED) = byteval(1) 

	ENDIF	 ! LIMIT_ON FLG_ATT_ALT_RED and Retstat is success

	IF ( LIMIT_ON ( FLG_ATT_ALT_YEL).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)  .and. (.not. no_attitude)
	1	.AND. 
	1	RETSTAT.EQ.SUCCESS ) THEN   
!	  IF ( Science_Record Sun_Angle is outside FEX_SCILIM database yellow
!	    limits Sun_Angle ), SET Bit 0 in FLG_ATT_ALT_YEL
	  IF ( Sci_Rec.Attitude.Sun_Angle .LT. 
	1      Scilim_Rec.Lim(2).Attitude.Sun_Angle ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL)
	    bytemask(1) = Gmask(1)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL) = byteval(1) 
	  ENDIF

!	  IF ( Science_Record Earth_Limb is outside FEX_SCILIM database yellow
!	    limits Earth_Limb ),SET Bit 1 in FLG_ATT_ALT_YEL
	  IF ( Sci_Rec.Attitude.Earth_Limb .LT. 
	1      Scilim_Rec.Lim(2).Attitude.Earth_Limb ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL)
	    bytemask(1) = Gmask(2)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL) = byteval(1) 
	  ENDIF

!	  IF ( Science_Record Moon_Angle is outside FEX_SCILIM database yellow
!	    limits Moon_Angle ), SET Bit 2 in FLG_ATT_ALT_YEL
	  IF ( Sci_Rec.Attitude.Moon_Angle .LT. 
	1      Scilim_Rec.Lim(2).Attitude.Moon_Angle ) THEN
	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL)
	    bytemask(1) = Gmask(3)
	    i2val = i2val .OR. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL) = byteval(1) 
	  ENDIF

!	  CALL FUT_CHECK_SAA to determine whether Science_Record COBE latitude
!	  and longitude is within FEX_SCILIM database yellow limits COBE 
!	  latitude and longitude for the South Atlantic Anomaly and set 
!	  Bit 3 in FLG_ATT_ALT_YEL
!
!	  CALL FUT_CHECK_VAB to determine whether Science_Record COBE 
!	  latitude and longitude is within FEX_SCILIM database yellow limits
!	  COBE latitude and longitude for the Van Allen Belts and set 
!	  Bit 4 in FLG_ATT_ALT_YEL
!
!	AND with Limflags FLG_ATT_ALT_YEL bit mask.

	    byteval(1) = Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL)
	    bytemask(1) = Limflags.Lim_Flags.FLG_ATT_ALT_YEL
	    i2val = i2val .AND. i2mask
	    Sci_Rec.DQ_Data.Data_Quality(FLG_ATT_ALT_YEL) = byteval(1) 

	ENDIF	 ! Lim_On FLG_ATT_ALT_YEL and Retstat is success

!       If the yellow or red moon angle limits are exceeded, set the Moon
!       flag to one or two, accordingly.

!	IF ( LIMIT_ON ( FLG_MOON) .and.(Sci_Rec.DQ_Data.Data_Quality(
!	1                    FLG_BADHKP) .eq. 0) .AND. 
!	1	RETSTAT.EQ.SUCCESS  .and. (.not. no_attitude)) THEN   
!	  IF ( Sci_Rec.Attitude.Moon_Angle .LT. 
!	1      Scilim_Rec.Lim(2).Attitude.Moon_Angle ) THEN
!	    Sci_Rec.DQ_Data.Data_Quality(FLG_MOON) = One
!	  ENDIF
!	  IF ( Sci_Rec.Attitude.Moon_Angle .LT. 
!	1      Scilim_Rec.Lim(1).Attitude.Moon_Angle ) THEN
!	    Sci_Rec.DQ_Data.Data_Quality(FLG_MOON) = Two
!	  ENDIF
!	ENDIF  ! FLG_MOON limit check is enabled
!
!	C  GROUP 4 : GRT Limits
!
!         Perform the following for Low Grts and then Hi Grts. A previous
!         routine has set all GRT red and yellow limits Limit_On flags to
!         either true or false depending on whether all red or yellow bits
!         are enabled or disabled. Thus only one Limit_On flag needs to be 
!         checked for red or yellow. The Limflags bit masks will remove the
!         disabled limit checks on the AND function. The calibrator resistors
!         are not to be checked for limits and will be disabled from checking.
!         If the dwell_stat(1) is equal to zero (1 = A side), do the A side.
!         If the dwell_stat(2) is equal to zero (1 = B side), do the B side.
!
	IF ( ( RETSTAT .EQ. SUCCESS ) .and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0) .and.
	1  ( Eng_Rec.En_Stat.Dwell_Stat(1) .EQ. Zero ) ) THEN

!	  DO for each GRT on A side
!	      If ( Engineering_Record GRT < FEX_ENGLIM database
!                red low limits for GRT OR > FEX_ENGLIM database
!		 red high limits for GRT ) 
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTRED_ST flags 

	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTRED_ST(Ix)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTRED((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTYEL_ST(Ix)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTYEL((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO

	  DO Ix = 1, 16
           If((Ix .lt. 11) .or. (Ix .gt. 14)) then  
	    IF ( LIMIT_ON_GRTRED (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.A_Lo_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(1).En_Analog.A_Lo_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.A_Lo_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(4).En_Analog.A_Lo_Grt(Ix)) .AND.
	1	  (Eng_Rec.En_Analog.A_Lo_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTRED_ST + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1)
                If(first_chan) then    
                  If (.not. first_redlim) then
                      Write(report_lun,900)eng_rec.ct_head.gmt
                      First_redlim=.true.
                  Endif
                  P_NAME=Grt_name(Ix)
                  Write(report_lun,901)P_NAME, Eng_rec.en_analog.a_lo_grt(ix)
                Endif
	      ENDIF 	
	    ENDIF

!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                yellow low limits for GRT OR > FEX_ENGLIM database
!		 yellow high limits for GRT ) 
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTYEL_ST flags
	
	    IF ( LIMIT_ON_GRTYEL (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.A_Lo_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(2).En_Analog.A_Lo_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.A_Lo_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(3).En_Analog.A_Lo_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.A_Lo_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTYEL_ST + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
	      ENDIF 	
             Endif
           Endif
	  ENDDO		    

	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTRED_ST(Ix+2)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTRED((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTYEL_ST(Ix+2)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTYEL((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 16
           If((Ix .lt. 11) .or. (Ix .gt. 14)) then  
	    IF ( LIMIT_ON_GRTRED (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.A_Hi_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(1).En_Analog.A_Hi_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.A_Hi_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(4).En_Analog.A_Hi_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.A_Hi_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTRED_ST + 2 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
                If(first_chan) then    
                  If (.not. first_redlim) then
                    Write(report_lun,900)eng_rec.ct_head.gmt
                    First_redlim=.true.
                  Endif
                  P_NAME=GRT_NAME(16+IX)
                  Write(report_lun,901)P_NAME, Eng_rec.en_analog.a_hi_grt(ix)
               Endif 
	      ENDIF 	
	    ENDIF
	    IF ( LIMIT_ON_GRTYEL (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.A_Hi_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(2).En_Analog.A_Hi_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.A_Hi_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(3).En_Analog.A_Hi_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.A_Hi_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTYEL_ST  + 2 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
	      ENDIF 	
	    ENDIF
           Endif
	  ENDDO
	ENDIF	! Retstat is success and not dwell mode -- A side GRTs

	IF ( ( RETSTAT .EQ. SUCCESS ) .AND.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0) .and.
	1  ( Eng_Rec.En_Stat.Dwell_Stat(2) .EQ. Zero ) ) THEN

!	  Perform the following for the Low Grts and then the Hi Grts.
!	  DO for each GRT on B side
!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                red low limits for GRT OR > FEX_ENGLIM database
!		 red high limits for GRT )
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTRED_ST flags

	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTRED_ST(Ix+4)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTRED((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTYEL_ST(Ix+4)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTYEL((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 16
           If((Ix .lt. 11) .or. (Ix .gt. 14)) then  
	    IF ( LIMIT_ON_GRTRED (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.B_Lo_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(1).En_Analog.B_Lo_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.B_Lo_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(4).En_Analog.B_Lo_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.B_Lo_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTRED_ST + 4 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
                If(first_chan) then    
                  If (.not. first_redlim) then
                    Write(report_lun,900)eng_rec.ct_head.gmt
                    First_redlim=.true. 
                  Endif
                  P_NAME=GRT_NAME(32+IX)
                  Write(report_lun,901)P_NAME, Eng_rec.en_analog.b_lo_grt(ix)
                endif
	      ENDIF 	
	    ENDIF

!	      IF ( Engineering_Record GRT < FEX_ENGLIM database
!                yellow low limits for GRT OR > FEX_ENGLIM database
!		 yellow high limits for GRT ) 
!	         SET the corresponding bit for the GRT in
!		 FLG_GRTYEL_ST flags.
	
	    IF ( LIMIT_ON_GRTYEL (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.B_Lo_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(2).En_Analog.B_Lo_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.B_Lo_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(3).En_Analog.B_Lo_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.B_Lo_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTYEL_ST + 4 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
	      ENDIF 	
	    ENDIF
           Endif
 	  ENDDO		    

	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTRED_ST(Ix+6)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTRED((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 2
	    grtflags = Limflags.Lim_Flags.FLG_GRTYEL_ST(Ix+6)
	    DO Jx = 1, 8
	      LIMIT_ON_GRTYEL((Ix-1)*8 + Jx) = BITest( grtflags, Jx-1)
	    ENDDO
	  ENDDO
	  DO Ix = 1, 16
           If((Ix .lt. 11) .or. (Ix .gt. 14)) then  
	    IF ( LIMIT_ON_GRTRED (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.B_Hi_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(1).En_Analog.B_Hi_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.B_Hi_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(4).En_Analog.B_Hi_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.B_Hi_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTRED_ST + 6 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
                If(first_chan) then    
                  If (.not. first_redlim) then
                    Write(report_lun,900)eng_rec.ct_head.gmt
                    First_redlim=.true.
                  Endif
                  P_NAME=GRT_NAME(48+IX)
                  Write(report_lun,901)P_NAME, Eng_rec.en_analog.b_hi_grt(ix)
                Endif
	      ENDIF 	
	    ENDIF
	    IF ( LIMIT_ON_GRTYEL (Ix )) THEN
	      IF ((Eng_Rec.En_Analog.B_Hi_Grt(Ix) .LT.
	1	  Englim_Rec.Lim(2).En_Analog.B_Hi_Grt(Ix) .OR.
	1	  Eng_Rec.En_Analog.A_Hi_Grt(Ix) .GT.
	1	  Englim_Rec.Lim(3).En_Analog.B_Hi_Grt(Ix) ) .AND.
	1	  (Eng_Rec.En_Analog.B_Hi_Grt(Ix) .ne. FV)) THEN
	          flag = FLG_GRTYEL_ST + 6 + ((Ix-1) / 8 )
	          byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	          bytemask(1) = Gmask(Ix)
	    	  i2val = i2val .OR. i2mask
		  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1) 
	      ENDIF 	
	    ENDIF
           Endif
	  ENDDO
	ENDIF	! Retstat is success and not dwell mode B side Grts

!       The 3.2 FDQ version performed the following algorithm:
!	AND the GRT data quality flags with the bit masks in Limflags.
!       The 4.1 version was changed to check the bits in the Limflags
!       first before doing the limit checking.

!	AND the even GRT data quality flags with the cal resistor mask.

	flag = FLG_GRTRED_ST - 1
	Do ix = 2, 8, 2
	  flag = flag + 2
	  byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	  bytemask(1) = calmask
	  i2val = i2val .AND. i2mask
	  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1)
	Enddo
	flag = FLG_GRTYEL_ST - 1
	Do ix = 2, 8, 2
	  flag = flag + 2
	  byteval(1) = Sci_Rec.DQ_Data.Data_Quality(flag)
	  bytemask(1) = calmask
	  i2val = i2val .AND. i2mask
	  Sci_Rec.DQ_Data.Data_Quality(flag) = byteval(1)
	Enddo

!	If both A and B sides are in dwell mode, set all of the GRT quality
!	flags because all GRTs will have the flag value and no temperatures
!       will be available from either side.

	IF (( RETSTAT .EQ. SUCCESS ).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0) .AND. 
	1	( Eng_Rec.En_Stat.Dwell_Stat(1) .NE. Zero ) .AND.
	1	( Eng_Rec.En_Stat.Dwell_Stat(2) .NE. Zero )) THEN

 	   Do ix = 11, 26
	      Sci_Rec.DQ_Data.Data_Quality(ix) = All_Set
	   Enddo
	ENDIF

!	Removed the checking of the 8 temperature controllers. The 12 bit words
!       in the Status-Monitor-Commands relating to commands for the reference 
!	source, reference horn, skyhorn and external calibrator for A and B 
!       sides cannot be converted to temperatures with a set of polynomial 
!	coefficients as other engineering database words. Furthermore, the
!       complete range of 4096 values ( 12 bits ) is a valid commanded word
!       and therefore no red and yellow limits can be defined. The temperature
!       controller (TEMP_CTRL) fields will not be used in limit checking. 
!
!	C  GROUP 5 : ENGINEERING LIMITS
!
!	C  Check channel related words and spacecraft side specific words in
!       C  the group. The flag setting algorithm applies to the following set:
!	C  FLG_CHAN_TMP_ST, FLG_CHPA_TMP, FLG_OPPA_TMP and FLG_BOL_VOL_ST.
!
!	DO for the four engineering flags, check yellow and then red limits.

	IF ( (RETSTAT .EQ. SUCCESS).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)) THEN 
          IF (( Limflags.Lim_Flags.FLG_CHAN_TMP_ST(chan)) .AND.
	1   ( Eng_Rec.En_Analog.Cna_Temp(chan) .ne. FV )) THEN
	    IF ( Eng_Rec.En_Analog.Cna_Temp(chan) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Cna_Temp(chan) .OR. 
	1	Eng_Rec.En_Analog.Cna_Temp(chan) .GT.
	1	Englim_Rec.Lim(3).En_Analog.Cna_Temp(chan) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_CHAN_TMP_ST) = one
	    ENDIF
	    IF ( Eng_Rec.En_Analog.Cna_Temp(chan) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Cna_Temp(chan) .OR. 
	1	Eng_Rec.En_Analog.Cna_Temp(chan) .GT.
	1	Englim_Rec.Lim(4).En_Analog.Cna_Temp(chan) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_CHAN_TMP_ST) = two
              If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='CHANNEL '//Chan_name(chan)// 'TEMP'
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.cna_temp(chan)
              Endif
	    ENDIF
	  ENDIF	   
          IF (( LIMIT_ON(FLG_CHPA_TMP) ) .AND.
	1   ( Eng_Rec.En_Analog.Pamp_Chan .ne. FV )) THEN
	    IF ( Eng_Rec.En_Analog.Pamp_Chan .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Pamp_Chan .OR. 
	1	Eng_Rec.En_Analog.Pamp_Chan .GT.
	1	Englim_Rec.Lim(3).En_Analog.Pamp_Chan ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_CHPA_TMP) = one
	    ENDIF
	    IF ( Eng_Rec.En_Analog.Pamp_Chan .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Pamp_Chan .OR. 
	1	Eng_Rec.En_Analog.Pamp_Chan .GT.
	1	Englim_Rec.Lim(4).En_Analog.Pamp_Chan ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_CHPA_TMP) = two
              If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='CHANNEL PREAMP '
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.PAMP_CHAN
              Endif
	    ENDIF
	  ENDIF	   

          IF (( LIMIT_ON(FLG_OPPA_TMP) ) .AND.
	1   ( Eng_Rec.En_Analog.Pamp_Op .ne. FV )) THEN
	    IF ( Eng_Rec.En_Analog.Pamp_Op .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Pamp_Op .OR. 
	1	Eng_Rec.En_Analog.Pamp_Op .GT.
	1	Englim_Rec.Lim(3).En_Analog.Pamp_Op ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_OPPA_TMP) = one
	    ENDIF
	    IF ( Eng_Rec.En_Analog.Pamp_Op .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Pamp_Op .OR. 
	1	Eng_Rec.En_Analog.Pamp_Op .GT.
	1	Englim_Rec.Lim(4).En_Analog.Pamp_Op ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_OPPA_TMP) = two
              If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='OPTICAL PREAMP'
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.PAMP_OP
              Endif
	    ENDIF
	  ENDIF	   

          IF (( Limflags.Lim_Flags.FLG_BOL_VOL_ST(chan) ) .AND.
	1   ( Eng_Rec.En_Analog.bol_volt(chan) .ne. FV )) THEN
	    IF ( Eng_Rec.En_Analog.bol_volt(chan) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.bol_volt(chan) .OR. 
	1	Eng_Rec.En_Analog.bol_volt(chan) .GT.
	1	Englim_Rec.Lim(3).En_Analog.bol_volt(chan) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_BOL_VOL_ST) = one
	    ENDIF
	    IF ( Eng_Rec.En_Analog.bol_volt(chan) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.bol_volt(chan) .OR. 
	1	Eng_Rec.En_Analog.bol_volt(chan) .GT.
	1	Englim_Rec.Lim(4).En_Analog.bol_volt(chan) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_BOL_VOL_ST) = two
              If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='BOLOMETER VOLTAGES '//Chan_name(chan)
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.bol_volt(chan)
              Endif
	    ENDIF
	  ENDIF	   
	ENDIF  		! Retstat is success
!
!	  The remainder of the engineering limits flags refer to engineering
!       words which have an A and B side. The side read in telemetry may
!	depend on the spacecraft side controlling FIRAS.
!	  Some flags have two bytes. The original requirements specified that
!       byte one referred to the side being checked and byte two to the other 
!	side. This requirement has been modified for Build 3.2 to a check for
!       the A side on byte 1 and B side on byte 2. 
!	  Both sides will be checked when the field values from both sides may
!       affect the IFGs. Many limit checks have been modified for Build 3.2 
!	and may be modified again for Build 4.0 after the requirements have
!       been clearly defined. Algorithms for Build 4.0 account for the fact 
!       that the spacecraft sides read only A or B FIRAS side components.
!       If the proper SC_Side is not on, the limit disable flas are set in a 
!       previous routine and the limit checks are therefore not done. Build
!       4.0 also includes the redefinition of the two IPDU temperature flags to
!       analog converter and digital converter and limit checking accordingly. 
!
	IF ( (RETSTAT .EQ. SUCCESS) .and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)) THEN
	  flag = FLG_IPDU_TMP - 1
	  DO Ix = 1, 2
	    flag = flag + 1 
	    IF (( LIMIT_ON(flag) ) .AND.
	1   ( Eng_Rec.En_Analog.IPDU_Temp(Ix) .NE. FV )) THEN
	      IF ( Eng_Rec.En_Analog.IPDU_Temp(Ix) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.IPDU_Temp(Ix) .OR. 
	1	Eng_Rec.En_Analog.IPDU_Temp(Ix).GT.
	1	Englim_Rec.Lim(3).En_Analog.IPDU_Temp(Ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = one
	      ENDIF
	      IF ( Eng_Rec.En_Analog.IPDU_Temp(Ix) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.IPDU_Temp(Ix) .OR. 
	1	Eng_Rec.En_Analog.IPDU_Temp(Ix) .GT.
	1	Englim_Rec.Lim(4).En_Analog.IPDU_Temp(Ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = two
              If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                If (IX .eq. 1 ) Side1 = '(SIDE A)'
                If (IX .eq. 2 ) Side1 = '(SIDE B)'
                P_NAME='IDPU TEMP '// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.IPDU_TEMP(IX)
              Endif
	      ENDIF
	    ENDIF  	! LIMIT_ON(flag)	   
	  ENDDO		! Ix = 1, 2

!	The drive box temperatures will not be checked for Build 3.2.
!	The thermistor being read depends on which side of the spacecraft
!	is controlling FIRAS. The DBX_TEMP(2) does not refer simply to FIRAS A 
!	and B, but refers to the fact that spacecraft side 1 reads A only and
!       spacecraft side 2 reads B only. The check is redone for Build 4.0.
!
	  flag = FLG_DBOX_TMP - 1
	  DO ix = 1, 2	
	    flag = flag + 1
	    IF (( LIMIT_ON(flag) ) .AND.
	1   ( Eng_Rec.En_Analog.Dbx_Temp(ix) .NE. FV )) THEN
	      IF ( Eng_Rec.En_Analog.Dbx_Temp(ix) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Dbx_Temp(ix) .OR. 
	1	Eng_Rec.En_Analog.Dbx_Temp(ix).GT.
	1	Englim_Rec.Lim(3).En_Analog.Dbx_Temp(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = one
	      ENDIF
	      IF ( Eng_Rec.En_Analog.Dbx_Temp(ix) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Dbx_Temp(ix) .OR. 
	1	Eng_Rec.En_Analog.Dbx_Temp(ix) .GT.
	1	Englim_Rec.Lim(4).En_Analog.Dbx_Temp(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                If (IX .eq. 1 ) Side1 = '(SIDE A)'
                If (IX .eq. 2 ) Side1 = '(SIDE B)'
                P_NAME='Dri. Box temp'// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.DBX_TEMP(IX)
               Endif
	      ENDIF
	    ENDIF
	  ENDDO	   
!
!	  For Build 3.2 both sides of the status monitor temperatures
!  	  will be checked. Control from the spacecraft side may have to be
!	  checked for future builds since there are thermistors for each of
!	  the spacecraft sides. For Build 4.0 both A and B FIRAS sides are
!         checked. There is only one flag here. The worst case is reflected.

	  DO Ix = 1, 2
	    Statmon_Flag(ix) = .false.
	    IF ((LIMIT_ON(FLG_STMON_TMP)) .AND. 
	1       ( Eng_Rec.En_Analog.Stat_Mon_Temp(ix) .NE. FV ) .AND.
	1	(Limflags.Lim_Flags.Flg_Stmon_Tmp(ix) .ne. zero) ) THEN
	      IF ( Eng_Rec.En_Analog.Stat_Mon_Temp(Ix) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Stat_Mon_Temp(Ix) .OR. 
	1	Eng_Rec.En_Analog.Stat_Mon_Temp(Ix) .GT.
	1	Englim_Rec.Lim(3).En_Analog.Stat_Mon_Temp(Ix) ) THEN 
	        IF (Sci_Rec.DQ_Data.Data_Quality(FLG_STMON_TMP) .ne. two)
	1	  Sci_Rec.DQ_Data.Data_Quality(FLG_STMON_TMP) = one
	        IF (Plot .eq. Yellow) Statmon_Flag(ix) = .true.
	      ENDIF
	      IF ( Eng_Rec.En_Analog.Stat_Mon_Temp(Ix) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Stat_Mon_Temp(Ix) .OR. 
	1	Eng_Rec.En_Analog.Stat_Mon_Temp(Ix) .GT.
	1	Englim_Rec.Lim(4).En_Analog.Stat_Mon_Temp(Ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_STMON_TMP) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                If (IX .eq. 1 ) Side1 = '(SIDE A)'
                If (IX .eq. 2 ) Side1 = '(SIDE B)'
                P_NAME='Stat. Mon. Temp'// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.STAT_MON_TEMP(IX)
               Endif
	       If ((Plot .eq. Red) .or. (Plot .eq. Yellow)) 
	1	   Statmon_Flag(ix) = .true. 
	      ENDIF
	    ENDIF	! Limit_On FLag is enabled and Limflags flag ne 0
	  ENDDO		! Ix = 1, 2
!
!	IPDU voltages and IPDU currents have odd and even array positions for
!	the A and B sides, respectively. Check both sides.

	  DO Ix = 1, 10
	    Jx = FLG_IPDU_VOL_ST + Ix - 1
	    DO side = 1, 2
	      Volts_Flag(side,ix) = .false.
	      IF (( LIMIT_ON(Jx) ) .AND.
	1       (Eng_Rec.En_Analog.IPDU_Volt(OEmask(Ix,side)) .NE. FV)) THEN
	        IF ( Eng_Rec.En_Analog.IPDU_Volt(OEmask(Ix,side)) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.IPDU_Volt(OEmask(Ix,side)) .OR. 
	1	Eng_Rec.En_Analog.IPDU_VOlt(OEmask(Ix,side)).GT.
	1	Englim_Rec.Lim(3).En_Analog.IPDU_Volt(OEmask(Ix,side))) THEN 
	    	  Sci_Rec.DQ_Data.Data_Quality(Jx) = one
	          If (Plot .eq. Yellow) Volts_Flag(side,ix) = .true.
	        ENDIF
	      ENDIF       ! LIMIT_ON(Jx)
	    ENDDO		! Side = 1, 2
	    DO side = 1, 2
	      IF (( LIMIT_ON(Jx) ) .AND.
	1       (Eng_Rec.En_Analog.IPDU_Volt(OEmask(Ix,side)) .NE. FV)) THEN
	        IF ( Eng_Rec.En_Analog.IPDU_Volt(OEmask(Ix,side)) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.IPDU_Volt(OEmask(Ix,side)) .OR. 
	1	Eng_Rec.En_Analog.IPDU_Volt(OEmask(Ix,side)) .GT.
	1	Englim_Rec.Lim(4).En_Analog.IPDU_Volt(OEmask(Ix,side))) THEN 
	    	  Sci_Rec.DQ_Data.Data_Quality(Jx) = two
                If(first_chan) then    
                  If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                    First_redlim=.true.
                  Endif
                  If (side .eq. 1 ) Side1 = '(SIDE A)'
                  If (Side .eq. 2 ) Side1 = '(SIDE B)'
                  P_NAME= Ipdu_Volt_Name(IX)// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.IPDU_VOLT
	1	     (OEmask(IX,side))
                 Endif
	          If ((Plot .eq. Red) .or. (Plot .eq. Yellow)) 
	1	     Volts_Flag(side,ix) = .true.
	        ENDIF
	      ENDIF	! LIMIT_ON(Jx)
	    ENDDO		! Side = 1, 2
 	  ENDDO		! Ix = 1, 10

	  DO Ix = 1, 6
	    Jx = FLG_IPDU_CUR_ST + Ix - 1
	    DO side = 1, 2
	      Cur_Flag(side,ix) = .false.
	      IF (( LIMIT_ON(Jx) ) .AND.
	1       (Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)) .NE. FV)) THEN
	        IF ( Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.IPDU_Amp(OEmask(Ix,side)) .OR. 
	1	Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)).GT.
	1	Englim_Rec.Lim(3).En_Analog.IPDU_Amp(OEmask(Ix,side)) ) THEN 
	    	  Sci_Rec.DQ_Data.Data_Quality(Jx) = one
	          If (Plot .eq. Yellow) Cur_Flag(side,ix) = .true.
	        ENDIF
	      ENDIF	! LIMIT_ON(Jx)
	    ENDDO		! Side = 1, 2
	    DO side = 1, 2
	      IF (( LIMIT_ON(Jx) ) .AND.
	1       (Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)) .NE. FV)) THEN
	        IF ( Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.IPDU_Amp(OEmask(Ix,side)) .OR. 
	1	Eng_Rec.En_Analog.IPDU_Amp(OEmask(Ix,side)) .GT.
	1	Englim_Rec.Lim(4).En_Analog.IPDU_Amp(OEmask(Ix,side)) ) THEN 
	    	  Sci_Rec.DQ_Data.Data_Quality(Jx) = two
                If(first_chan) then    
                  If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                    First_redlim=.true.
                  Endif
                  If (side .eq. 1 ) Side1 = '(SIDE A)'
                  If (Side .eq. 2 ) Side1 = '(SIDE B)'
                  P_NAME= Ipdu_Cur_Name(Ix)// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.IPDU_AMP
	1	     (OEmask(IX,side))
                 Endif
	          If ((Plot .eq. Red) .or. (Plot .eq. Yellow)) 
	1	     Cur_Flag(side,ix) = .true.
	        ENDIF
	      ENDIF	! LIMIT_ON(Jx)
	    ENDDO		! Side = 1, 2	   
 	  ENDDO		! Ix = 1, 6

!	Hot spot heaters are checked on both sides regardless of channel.

	  flag = FLG_HS_HEAT_A - 1
	  DO ix = 1, 2
	    flag = flag + 1
	    IF((LIMIT_ON(flag)) .AND. (RETSTAT .EQ. SUCCESS) .AND.
	1     ( Eng_Rec.En_Analog.Hot_Spot(Ix) .NE. FV )) THEN
	      IF ( Eng_Rec.En_Analog.Hot_Spot(ix) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Hot_Spot(ix) .OR. 
	1	Eng_Rec.En_Analog.Hot_Spot(ix).GT.
	1	Englim_Rec.Lim(3).En_Analog.Hot_Spot(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = one
	      ENDIF
	      IF ( Eng_Rec.En_Analog.Hot_Spot(ix).LT. 
	1	Englim_Rec.Lim(1).En_Analog.Hot_Spot(ix) .OR. 
	1	Eng_Rec.En_Analog.Hot_Spot(ix) .GT.
	1	Englim_Rec.Lim(4).En_Analog.Hot_Spot(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                If (IX .eq. 1 ) Side1 = '(SIDE A)'
                If (IX .eq. 2 ) Side1 = '(SIDE B)'
                P_NAME='Hot Spot'// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.hot_spot(IX)
               Endif
	      ENDIF
	    ENDIF	
	  ENDDO		! ix = 1, 2 for hot spot heaters   

!	The algorithm for checking the MTM cal motor is added for Build 4.0.
!       The currents for A and B sides will be checked if the Limit_On flag
!       is set.

	  flag = FLG_MTMCAL_MTR - 1
	  DO ix = 1, 2	
	    flag = flag + 1
	    IF (( LIMIT_ON(flag) ) .AND.
	1     ( Eng_Rec.En_Analog.Mtm_Cal_Mtr(Ix) .NE. FV )) THEN
	      IF ( Eng_Rec.En_Analog.Mtm_Cal_Mtr(ix) .LT. 
	1	Englim_Rec.Lim(2).En_Analog.Mtm_Cal_Mtr(ix) .OR. 
	1	Eng_Rec.En_Analog.Mtm_Cal_Mtr(ix).GT.
	1	Englim_Rec.Lim(3).En_Analog.Mtm_Cal_Mtr(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = one
	      ENDIF
	      IF ( Eng_Rec.En_Analog.Mtm_Cal_Mtr(ix) .LT. 
	1	Englim_Rec.Lim(1).En_Analog.Mtm_Cal_Mtr(ix) .OR. 
	1	Eng_Rec.En_Analog.Mtm_Cal_Mtr(ix) .GT.
	1	Englim_Rec.Lim(4).En_Analog.Mtm_Cal_Mtr(ix) ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(flag) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                If (IX .eq. 1 ) Side1 = '(SIDE A)'
                If (IX .eq. 2 ) Side1 = '(SIDE B)'
                P_NAME='MTM/CAL Motors '// Side1
                Write(report_lun,901)P_NAME, Eng_rec.en_analog.MTM_CAL_MTR(IX)
               Endif
	      ENDIF
	    ENDIF
	  ENDDO	   

!       New LMAC temperatures for analog and digital converters added for 
!        Build 4.0.	

          IF (( LIMIT_ON(FLG_ACLMAC_TMP) ) .AND.
	1     ( Eng_Rec.En_Tail.LMAC_Analog_Temp .NE. FV )) THEN
	      IF ( Eng_Rec.En_Tail.LMAC_Analog_Temp .LT. 
	1	Englim_Rec.Lim(2).En_Tail.LMAC_Analog_Temp .OR. 
	1	Eng_Rec.En_Tail.LMAC_Analog_Temp .GT.
	1	Englim_Rec.Lim(3).En_Tail.LMAC_Analog_TEmp ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_ACLMAC_TMP) = one
	      ENDIF
	      IF ( Eng_Rec.En_Tail.LMAC_Analog_Temp .LT. 
	1	Englim_Rec.Lim(1).En_Tail.LMAC_Analog_Temp .OR. 
	1	Eng_Rec.En_Tail.LMAC_Analog_Temp .GT.
	1	Englim_Rec.Lim(4).En_Tail.LMAC_Analog_Temp ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_ACLMAC_TMP) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='LMAC Analog Temp'
                Write(report_lun,901)P_NAME, Eng_rec.en_tail.LMAC_ANALOG_TEMP
               Endif
	      ENDIF
	  ENDIF	   
          IF (( LIMIT_ON(FLG_DCLMAC_TMP) ) .AND.
	1     ( Eng_Rec.En_Tail.LMAC_Digital_Temp .NE. FV )) THEN
	      IF ( Eng_Rec.En_Tail.LMAC_Digital_Temp .LT. 
	1	Englim_Rec.Lim(2).En_Tail.LMAC_Digital_Temp .OR. 
	1	Eng_Rec.En_Tail.LMAC_Digital_Temp .GT.
	1	Englim_Rec.Lim(3).En_Tail.LMAC_Digital_Temp ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_DCLMAC_TMP) = one
	      ENDIF
	      IF ( Eng_Rec.En_Tail.LMAC_Digital_Temp .LT. 
	1	Englim_Rec.Lim(1).En_Tail.LMAC_Digital_Temp .OR. 
	1	Eng_Rec.En_Tail.LMAC_Digital_Temp .GT.
	1	Englim_Rec.Lim(4).En_Tail.LMAC_Digital_Temp ) THEN 
	    	Sci_Rec.DQ_Data.Data_Quality(FLG_DCLMAC_TMP) = two
               If(first_chan) then    
                If (.not. first_redlim) then
                  Write(report_lun,900)eng_rec.ct_head.gmt
                  First_redlim=.true.
                Endif
                P_NAME='LMAC Digital Temp'
                Write(report_lun,901)P_NAME, Eng_rec.en_tail.LMAC_digital_TEMP
               endif 
	      ENDIF
	  ENDIF	   

	ENDIF		! Retstat is success
!
!	C  GROUP 6 : MAJOR FRAME BOUNDARY CHANGES
!
!	C  Status changes across the Major Frame boundary, indicated by a 
!	C  change in state of a particular bit in the arrays of the status
!       C  monitor command in the two bracketing housekeeping records will
!	C  be used to set FLG_STCHG_MJ_ST. 
!
!	CALL FUT_MJF_CHANGE to extract the relevant bits from the status
!	     monitor command arrays in the two bracketing housekeeping
!	     records, check for a change in state and set the corresponding
!	     bit in the FLG_STCHG_MJ_ST to 1.
!
	If ((Retstat .eq. Success).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)
	1                  .and. (Limit_On(FLG_STCHG_MJ_ST))) then
	  STATUS = FUT_MJF_CHANGE ( HSK_RECS, Limflags, First_MJ, Chan, 
	1	 Sci_Rec.DQ_Data.Data_Quality(FLG_STCHG_MJ_ST))
	  if ( status .ne. %loc(FUT_NORMAL) ) then
	    if ( status .eq. %loc(FUT_ABERR)) status = zerostat
	    call lib$signal(FUT_MJFCHANGER,%val(1),%val(status))
	    retstat = error
C          else
C           If (sci_rec.dq_data.data_quality(flg_stchg_mj_st) .gt. 0 .or. 
C	1	sci_rec.dq_data.data_quality(flg_stchg_mj_st+1) .gt. 0 .or.
C	1	sci_rec.dq_data.data_quality(flg_stchg_mj_st+2) .gt. 0 .or.
C	1	sci_rec.dq_data.data_quality(flg_stchg_mj_st+3) .gt. 0 )
C	1  then
C              If (.not. first_redlim) then
C                Write(report_lun,900)chan_name(chan),eng_rec.ct_head.gmt
C                First_redlim=.true.
C              Endif
C              P_NAME='status changes acrs. maj. frm. bdy.'
C              Write(report_lun,899)P_NAME
C           Endif
          
          endif  

	Endif 		! Retstat is success
!
!	C  GROUP 7 : STATUS FROM DIRBE, DMR or SPACECRAFT AFFECTING FIRAS
!
!	CALL FUT_Xtalk_access to access DIRBE,DMR,SPACECRAFT
!           states that will affect FIRAS 
!	CALL FUT_XTALK_CHECK to check DIRBE,DMR SPACECRAFT states that will affect FIRAS 
!
	If ((Retstat .eq. Success).and.(Sci_Rec.DQ_Data.Data_Quality(
	1                    FLG_BADHKP) .eq. 0)
	1                  .and. (Limit_On(FLG_DIRBE_ST))
	1                  .and. (Limit_On(FLG_DMR_ST))
	1                  .and. (Limit_On(FLG_SPACECR_ST))) then
           STATUS= FUT_XTALK_ACCESS(DUMMY1,DUMMY2)
           STATUS = FUT_XTALK_CHECK(DUMMY3,DUMMY3)
        Endif
!
!	C  GROUP 8 : FIRAS DATA QUALITY SUMMARY FLAG
!
!	CALL FUT_SUM_FLG to count the red and yellow flags set for the
!	     entire Science_Record and determine a value for the overall
!	     summary quality flag, FLG_SUMMARY. The Moon_Angle flag will
!	     not be included in the count. 

	If (( Retstat .eq. Success) .and. (Limit_On(FLG_SUMMARY)) ) then
	  STATUS = FUT_SUM_FLG (Chan, Sci_Rec.DQ_Data.Data_Quality)
	  if ( status .ne. %loc(FUT_NORMAL) ) then
	    if ( status .eq. %loc(FUT_ABERR)) status = zerostat
	    call lib$signal(FUT_SUMFLAGER,%val(1),%val(status))
	    retstat = error
	  endif
	Endif 		! Retstat is success
!
!	C  A summary table will be kept for all individual flags set by
!	C  the routine for the entire segment. For each IFG the appropriate
!	C  data quality flag counters will be incremented is the flag indicates
!	C  bad quality. When the segment is completer, a summary will be sent
!	C  to the user via sys$out device in effect for the run. A message for
!	C  each quality flag will be included in FUT_MSG.MSG.
!
!	CALL FUT_QUAL_SUMMARY to increment the counters for the flags set for
!	     the channel Science Record and output the summary at the end of
!	     the segment.
!
	If ( Retstat .eq. Success ) then
	  STATUS = FUT_QUAL_SUMMARY ( Chan, Sci_Rec.DQ_Data.Data_Quality,     
	1	 End_Segment, Report_Lun)
	  if ( status .ne. %loc(FUT_NORMAL) ) then
	    if ( status .eq. %loc(FUT_ABERR)) status = zerostat
	    call lib$signal(FUT_QUALSUMER,%val(1),%val(status))
	    retstat = error
	  endif
	Endif 		! Retstat is success
!
!	CALL FUT_TRIGGER_PLOT to keep track of all limits exceeded for which
!       Engplots can be generated and write a command text file to be invoked
!       at the completion of the run.
!
	If (( Retstat .eq. Success ) .AND. ( Plot .ne. Zero)) then
	  STATUS = FUT_TRIGGER_PLOT ( Sci_Rec, Chan, Plot, Statmon_Flag,     
	1   Volts_Flag, Cur_Flag, Min_Offtime, Max_Offtime, Plot_Device,
	2   Plot_Table, Plot_Start, Plot_Stop, englim )
	  if ( status .ne. %loc(FUT_NORMAL) ) then
	    if ( status .eq. %loc(FUT_ABERR)) status = zerostat
	    call lib$signal(FUT_TRIGPLOTERR,%val(1),%val(status))
	    retstat = error
	  endif
	Endif 		! Retstat is success
!	
!	Set function to return status

	if (retstat.eq.success) then
	  FUT_GET_QUALFLAGS = %loc(FUT_NORMAL)
	else
	  FUT_GET_QUALFLAGS = %loc(FUT_ABERR)
	endif
899     Format(5x,a40)
900     Format(//x,'***********************************************',
	1     /x,'*','       Red Engineering Limits Exceeded       *',
	1     /x,'*       Eng. Record GMT = ', A14,5x,' *',                   
	1     /x,'***********************************************')
901     Format(5x,a40,' = ',f12.5)
	RETURN
	END
