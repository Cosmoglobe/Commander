	INTEGER*4 FUNCTION FDQ_GET_OPTIONS ( OPTIONS, TIME_RANGE, 
	1	FILE_SEG, REF_CHAN, 
	2	SET_SEC, SEC_DEF, ATTITUDE_TYPE, LIMIT,
	3	PLOT_DEVICE, MIN_OFFTIME, MAX_OFFTIME, 
	4	REPORT_FILE, CURGMT, REPORT,Report_def)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_OPTIONS
C/
C/	PROGRAM DESCRIPTION:
C/
C/	FDQ_GET_OPTIONS This routine reads the CLD file for processing
C/	options for running DATA QUALIFY and checks consistency of choices.
C/
C/	AUTHOR:
C/	  Shirley M. Read
C/	  STX
C/	  December 1987
C/
C/	MODIFIED BY:
C/	  In December 1987 the original subroutine was replaced by the
C/	  current function. The new routine uses the Command Line 
C/	  Interpreter to obtain the user's options in accord with the
C/	  rest of the Firas Pipeline. User options of writing an Ascii
C/	  conversion coefficients file, a full listing or brief listing
C/	  of this file and number of seconds to be used as a criterion
C/	  for grouping science records into sets were also included in
C/	  this new version. These parameters will be passed as arguments
C/	  back to the appropriate FDQ functions.
C/
C/	  Shirley M. Read
C/	  STX
C/	  March 1988.
C/	    For Build 3.2 COBETRIEVE included modifications which made it 
C/	  necessary to change the interface of FDQ with COBETRIEVE. The
C/	  catalog number is no longer available. The segment for modification
C/	  must be specified by a filename. A time range selection is still
C/	  available. The CT Modify Open was also changed and new catalog 
C/	  functions were added to relace the old information routines. Thus
C/	  FDQ was redesigned to use the new COBETRIEVE features. The selected
C/	  file_segment and the channel referenced in the filename are added
C/	  to the passed parameters. The user option for Exclude segments is
C/	  removed.
CH
CH	CHANGE LOG:
CH	
CH		Version 4.2.1 01/19/89, SPR 3150, Shirley M. Read, STX
CH	 	        FDQ needs to access coarse attitude for the CSDR I&T
CH			Acceptance Test for Quick Look on January 19, 1989.
CH			The coarse and predicted options were not implemented
CH	  		at the time the real attitude access was implemented
CH			in FDQ. The keywords "COARSE" and "PREDICTED" must
CH			be picked up from the command line and the appropriate
CH			FAC Params set accordingly.
CH
CH	Version 4.2.1 02/09/89, SPR 2550, Shirley M. Read, STX
CH		Provide user option of setting the time interval for triggered
CH		plots.
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH	Version 4.2.1 02/09/89, SPR 2955, Shirley M. Read, STX
CH		FDQ should be given the correct COBE orbital period which is
CH 		used to form engineering trends. Make the orbital period a
CH		command line option.
CH	Version 4.2.1 03/08/89, SPR 3384, Shirley M. Read, STX
CH	        The start and stop strings[B need to be padded with standard
CH	 	CT type padding for subsequent CT calls. CT_GMT_to_Binary
CH		will not take the semicolon in the string.		
CH	Version 4.4.1 07/20/89, SPR 4085, Shirley M. Read, STX
CH	        FDQ should call the new FUT routine to get the orbital period.
CH		Get the orbital period from the command line if the qualifier
CH	        is present. Otherwise the FUT routine will be used when the
CH	        data stop time is determined.
CH	Version 4.4.1 08/20/89 SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments related changes.
CH	Version 5.2 12/29/89 SPR 5536, R. Kummerer, STX
CH		/NOTRACK erroneously selects IFG tracking.
CH
CH	Version 6.1 5/15/90 SPR 6726, H. Wang, STX
CH		FDQ must support talaris printer
CH
CH	Version     8/8/90 SER 4171, Larry P. Rosen, STX
CH		Standardize report file names.  
CH
CH      Modified by: H. WANG, STX, 1/29/91
CH         Reason: New requirements for FDQ
CH               The following command qualifiers no longer be required:
CH               SCSIDE
CH               XCALSET
CH               TRACK
CH               DELTIM
CH               BINTIM
CH               ASCFILE
CH               LSTFULL
CH               LSTSHORT
CH               ORBITAL_PD
CH                 
C/
C/	CALLING SEQUENCE:
C/	  STATUS = FDQ_GET_OPTIONS ( OPTIONS, TIME_RANGE, FILE_SEG, 
C/	           REF_CHAN,  SET_SEC,SEC_DEF, 
C/	           ATTITUDE_TYPE, LIMIT, PLOT_DEVICE, MIN_OFFTIME,
C/		   MAX_OFFTIME, REPORT_FILE, CURGMT, REPORT,report_def)
C/
C/	INPUT PARAMETERS:
C/	  NONE
C/
C/	OUTPUT PARAMETERS:
C/	  OPTIONS(10)	I*2	Array of various flags indicating processing
C/				options:
C/				OPTIONS(1)  --	Flag indicating whether a time
C/						range is used:
C/						0 = No time range used (default)
C/						1 = Binary start time used
C/						2 = Specific set of segments
C/						    to process
C/						3 = Specific set of science
C/						    segments to process; read
C/						    HKP by timerange
C/				OPTIONS(2)  --	Flag indicating whether all or
C/						one segment(s) are(is) to be
C/						processed:
C/						0 = Process only the next segment
C/						1 = Process all segments (default)
C/				OPTIONS(3)  --	Flag indicating how many segments
C/						are to be excluded from processing
C/						in this run.
C/						0 is the default (that mean no
C/						segment is excluded by default)
C/
C/				OPTIONS(4) - OPTIONS(10) -- Reserved for future use
C/
C/	  TIME_RANGE	   C*30	Time range for processing interval if OPTIONS(1)
C/				equals 1(BST).
C/	  FILE_SEG	   C*39 Name of raw science file representing a segment
C/			        to be modified.
C/	  REF_CHAN	   I*4  Channel referenced in file name.
C/	  SET_SEC          I*4  Number of seconds to be used as a criterion
C/	 			for grouping of science records into sets.
C/				Default is 16 seconds.
C/        SEC_DEF          L*1  FLAG indicating whether set_sec is default 
C/	  ATTITUDE_TYPE    I*4  Type of attitude requested. Values:
C/				0 = none, 1 = simulated, 2 = predicted,
C/			        3 = coarse, 4 = fine without Dirbe, 
C/				5 = fine with Dirbe, 6 = definitive
C/	  LIMIT            I*2  Flag for generating Engplots. Values:
C/			        0 = none, 1 = red only, 2 = red or yellow
C/	  PLOT_DEV         I*4	Plot device: line or laser printer
C/	  MIN_OFFTIME	   I*4	Offset in hours from plot start
C/	  MAX_OFFTIME	   I*4	Offset in hours from plot stop
C/	  REPORT_FILE	   C*33	Report file name
C/	  CURGMT	   C*14 Current run time in GMT
C/	  REPORT	   L*1	Flag whether to write report or not
C/        Report_def       L*1  Flag indicate default command qualifier
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C	  UPM_Stat_Msg -- UPM facility status and message file
C	  $SSdef     -- System status values
C/	  FDQ_OPTIONS.TXT
C/	  FUT_PARAMS.TXT
C/
C/	SUBROUTINES CALLED:
C/	  UPM_PRESENT
C/	  UPM_GET_VALUE
C/	  UPM_GET_WORD
C/	  UPM_GET_LONGWORD
C/	  UPM_GET_FLOAT
C/	  STR$TRANSLATE
C/	  SYS$GETTIM
C/
C/	ERROR HANDLING:
C/	  Call Lib$Signal with FDQ defined parameters. A condition handler
C/	was established in FDQ to signal the errors and return control to
C/	the FDQ Program.
C/
C/	METHOD USED:
C/
C/	PDL FOR	FDQ_GET_OPTIONS 
C/
C/	BEGIN
C/
C/	Set return status to success
C/
C/	IF (UPM_PRESENT('SEGMENT')) THEN
C/	  Set the flag in OPTIONS array for one specified segment
C/	  SET the flag in OPTIONS array for process one segment
C/	  CALL UPM_GET_VALUE to get the filename for the selected segment
C/	  IF ( Status is not normal ) THEN			
C/	    CALL LIB$SIGNAL ( Status )
C/          Set return status to failure
C/	  ELSE
C/	    EXTRACT the channel identifier in the filename. 
C/	    IF (present) THEN
C/	       Set the reference channel ID
C/	    ELSEIF (not present) THEN
C/	       CALL LIB$SIGNAL ( Badfile name )
C/	       Set return status to failure
C/	    ENDIF
C/	  ENDIF
C/	ELSE
C/	  Set the flag in OPTIONS array for select by binary data time
C/	ELSE
C/ 	  CALL LIB$SIGNAL ( No file nor time range selected )
C/	  Set the return status to failure
C/	ENDIF
C/
C/	  IF ( flag for select by binary start time or delivery time is set )
C/	    THEN
C/	    IF(UPM_PRESENT('JSTART')) THEN
C/	      CALL UPM_GET_VALUE to get the start time
C/	      IF ( Status is not normal ) THEN			
C/	        CALL LIB$SIGNAL ( Status )
C/              Set the return status to failure
C/	      ENDIF
C/	    ENDIF
C/	    IF(UPM_PRESENT('JSTOP')) THEN
C/	      CALL UPM_GET_VALUE to get the stop time
C/	      IF ( Status is not normal ) THEN			
C/	        CALL LIB$SIGNAL ( Status )
C/              Set the return status to failure
C/	      ENDIF
C/	    ENDIF
C/	  ENDIF
C/
C/	  IF ( Return status is success and flag for no time range 
C/	    selected is not set) THEN
C/	    Verify consistency of start and stop times If inconsistent,
C/          set the flag for no time range selected
C/	    IF (JSTOP lt JSTART) THEN
C/	      CALL LIB$SIGNAL ( stop less than start )
C/	      Set the flag in OPTIONS array for no time range selected
C/	    ELSE
C/	      CALL STR$TRANSLATE to convert the start and stop times to
C/                 a standard Ascii format time range to be used for calls
C/                 to Cobetrieve
C/	    ENDIF
C/	  ENDIF
C/
C/
C/	IF (UPM_PRESENT('SETSEC')) THEN
C/	  CALL UPM_GET_VALUE to get the number of seconds for the set group 
C/	  IF ( Status is not normal ) THEN			
C/	    CALL LIB$SIGNAL ( Status )
C/          Set return status to failure
C/	  ENDIF
C/	ELSE
C/	  Set the number of seconds for the set group to default of 16 seconds
C/	ENDIF
C/	IF (UPM_PRESENT('XCALSET')) THEN
C/	  Set the Xcal flag to True
C/	ELSE
C/	  Set the Xcal flag to False
C/	ENDIF
C/	IF (UPM_PRESENT('ATTITUDE')) THEN
C/	  CALL UPM_GET_VALUE to get the selected type
C/	  Set the Attitude type accordingly
C/	ELSE
C/	  Set the Attitude type to None
C/	ENDIF
C/	IF (UPM_PRESENT('SCSIDE')) THEN
C/	  CALL UPM_GET_WORD to get the spacecraft side
C/	ENDIF
C/	IF (UPM_PRESENT('LIMIT')) THEN
C/	  CALL UPM_GET_VALUE to get the selected type
C/	  Set the Plot parameter to red or yellow, accordingly
C/	ELSE
C/	  Set the Plot parameter to None
C/	ENDIF
C/      IF the Plot flag is greater than zero, get the Plot Device.
C/	IF the Plot flag is greater than zero,
C/	   Get the Minimum offset time.
C/	   Get the Maximum offset time.
C/	ENDIF
C/	Call SYS$Gettim to get current time.
C/	If status normal, convert binary time to GMT
C/	ELSE call lib$signal too signal error
C/	Set report default name
C/	STATUS=UPM_PRESENT('REPORT')
C/	IF(STATUS .EQ. UPM_PRES) THEN
C/	   REPORT = .TRUE.
C/	   RET_STATUS = UPM_GET_VALUE (to get the report file name)
C/	   IF (RET_STATUS .EQ. UPM_ABSENT) THEN use default name.
C/	ELSE IF(STATUS .EQ. UPM_DEFAULTED) THEN use default name.
C/	ELSE (negated)
C/	   REPORT = .FALSE. No report file to be opened
C/	ENDIF
C/	RETURN
C/	END

	IMPLICIT	NONE

C	Input/Output Passed Parameters

	INTEGER*2	OPTIONS(10)
	CHARACTER	TIME_RANGE*30
	CHARACTER*39    FILE_SEG
	INTEGER*4	REF_CHAN
	INTEGER*4	SET_SEC
	LOGICAL*1	SEC_DEF
	INTEGER*4  	ATTITUDE_TYPE
	INTEGER*2       LIMIT	        ! Plot type: none, red or yellow
	INTEGER*4       PLOT_DEVICE     ! Flag for plot device if an Engplots
                                        ! command file is requested
        INTEGER*4       MIN_OFFTIME     ! Hours to subtract from first limit
                                        ! exceeded time for plot start
        INTEGER*4       MAX_OFFTIME     ! Hours to add to first limit
                                        ! exceeded time for plot stop
	CHARACTER*33	REPORT_FILE	! Report file name to be opened
	CHARACTER*14	CURGMT		! Current run system time in GMT
	LOGICAL*1	REPORT		! Flag whether to write report or not
        Logical*1       Report_def
C	Include Files, Functions and External Parameters

	INCLUDE 	'(FDQ_OPTIONS)'
	INCLUDE	 	'(FUT_PARAMS)'
	INCLUDE		'($SSDEF)'
	INCLUDE	        '(UPM_STAT_MSG)'	! Command line parameters
						! & functions
	INTEGER*4	UPM_PRESENT, UPM_GET_VALUE
	INTEGER*4       UPM_GET_WORD, UPM_GET_LONGWORD, UPM_GET_FLOAT
	INTEGER*4	LIB$SKPC	! Skip characters
	INTEGER*4	LIB$LOCC	! Locate characters
	INTEGER*4	SYS$GETTIM
	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_SELFILE
	EXTERNAL	FDQ_REFCHAN
	EXTERNAL	FDQ_SELTIMINT
	EXTERNAL	FDQ_SELASCFUL
	EXTERNAL	FDQ_SELASCBRF
	EXTERNAL	FDQ_SELSEC
	EXTERNAL	FDQ_SELBST
	EXTERNAL	FDQ_SELDDT
	EXTERNAL        FDQ_DDTNOTIMP
	EXTERNAL        FDQ_SELNOTIM
	EXTERNAL	FDQ_ERSELFIL
	EXTERNAL	FDQ_FILEBLANK
	EXTERNAL	FDQ_BADFILE
	EXTERNAL	FDQ_BADREFCHAN
	EXTERNAL	FDQ_ERTIMINT
	EXTERNAL        FDQ_SELNTOBIG
	EXTERNAL	FDQ_NEXZERO
	EXTERNAL	FDQ_ERSELSEC
 	EXTERNAL	FDQ_BADTIMINV
	EXTERNAL	FDQ_SECTOOBIG
	EXTERNAL	FDQ_GETATTYPERR
	EXTERNAL	FDQ_INVALATTSEL
	EXTERNAL	FDQ_GETSCSIDERR
	EXTERNAL	FDQ_INVALSCSIDE
	EXTERNAL	FDQ_GETPLOTERR
	EXTERNAL	FDQ_INVALPLOT
	EXTERNAL	FDQ_TYPATTSEL
	EXTERNAL	FDQ_SCSIDESEL
	EXTERNAL	FDQ_PLOTSEL
	EXTERNAL	FDQ_SELORBPD
	EXTERNAL 	FUT_NORMAL

C	Local Declarations

	INTEGER*4	DEFAULT_SEC / 16 /  !Default of 16 
	INTEGER*4	I	! Index
	INTEGER*4	IPOS1, IPOS2	! Line positions
	INTEGER*4	START, STOP
	CHARACTER*2	CHAN_ID(4) / 'RH', 'RL', 'LH', 'LL' /
	CHARACTER*2	CHAN_ID_L(4) / 'rh', 'rl', 'lh', 'll' /
	CHARACTER	START_TIME*14, STOP_TIME*14
	INTEGER*2	SECOND_SEMICOLON
	INTEGER*4 	STATUS		! System return status
	INTEGER*4 	STATUS2		! System return status
	INTEGER*2       VALUE_SIZE	! Length of qualifier.
	INTEGER*4	RETSTAT		! Return status indicator for function
	INTEGER*4 	SUCCESS / 1 /
	INTEGER*4	ERROR / 2 /	
	CHARACTER*20    ATT_TYPE
	CHARACTER*20    PLOT_TYPE
	CHARACTER*32    TXTVAL		! Dummy storage varuable
	INTEGER*2       TXTLEN
	CHARACTER*33	REPORT_DEFAULT  ! Default report file name
	INTEGER*4	CURTIME(2)	! System binary time


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       -----------------------------------------------------------------------
C	Get values of qualifiers.
C		1.  FILE determines whether a filename will be used
C		    to get the Firas Cobetrieve files.
C		2.  BINTIM or DELTIM determines whether the query time will 
C		    be based on Binary Start and Stop Time of the data or 
C		    Data Delivery Time. If no segment and no times are 
C		    specified, an error is signalled.
C		3.  JSTART determines the start time for Cobetrieve.
C		4.  JSTOP determines the stop time for Cobetrieve.
C	        5.  ASCFILE determines if an Ascii conversions file is written.
C		    The default is no Ascii file written.
C	        6.  If the file is requested, LSTFULL or LSTSHORT determines
C	            if full or brief listing will be written. The default is 
C	            a brief listing.
C	        7.  SETSEC determines the number of seconds to be used as
C		    a criterion for grouping sets of science records. If not
C	            present, a value of 16 seconds is used.
C	        8.  If present, XCALSET sets the flag to set the Xcal bit 
C		    in the archive record. The default is False.
C	        9.  If present, ATTITUDE, get value of attitude type. 
C		    Defaultis simulated.
C	       10.  If present, SCSIDE, get value of spacecraft side.
C	       11.  If present, PLOT, get value for plot flag.
C	       12.  If present or default, REPORT, open report file.
C
C	If an error occurs in a return from a UPM call, the status is signaled,
C	and the program is exited with an error status.
C       -----------------------------------------------------------------------

C	Set return status to success

	RETSTAT = SUCCESS

	IF (UPM_PRESENT('FILE')) THEN	! A specific file is specified

C	  CALL UPM_GET_VALUE to get the filename for the selected segment

	  Status = UPM_GET_VALUE('FILE',File_Seg,Value_Size)	

	  IF (Status .EQ. SS$_Normal) THEN			

C	  CHECK validity of filename and extract the reference channel ID.

	    OPTIONS(OPT_TIME_RANGE) = OPTV_SPEC_SEG
	    OPTIONS(OPT_ALL_SEGMENTS) = OPTV_ONE_SEGMENT

	    ipos1 = Lib$Skpc( ' ', File_Seg )
	    if ( ipos1 .eq. 0 ) then
	      RETSTAT = ERROR
	      CALL LIB$SIGNAL(FDQ_FILEBLANK)
	    else
	      ipos2 = Lib$Locc( '.', File_Seg(ipos1:39))
	      if ( ipos2 .eq. 0 ) then
	        RETSTAT = ERROR
	        CALL LIB$SIGNAL(FDQ_BADFILE, %val(1), File_Seg )
	      else
		start = ipos1 + ipos2 - 3
		stop = start + 1
	        ref_chan = 0
	        do i = 1,4
		  if ( (File_Seg(start:stop) .eq. chan_id(i)) .or.
	1	       (File_Seg(start:stop) .eq. chan_id_l(i))) then
		    ref_chan = i
	          endif
		enddo
	        if ( ref_chan .eq. 0 ) then
		  RETSTAT = ERROR
		  CALL LIB$SIGNAL(FDQ_BADREFCHAN, %val(1), File_Seg )
	        else
	          CALL LIB$SIGNAL(FDQ_SELFILE, %VAL(1), File_Seg)
	          CALL LIB$SIGNAL(FDQ_REFCHAN, %val(1), chan_id(ref_chan))
	        endif
	      endif
	    endif
	   
	  ELSE IF (Status .Eq. UPM_Absent) Then
	    OPTIONS(OPT_TIME_RANGE) = OPTV_HKP_TIMERANGE
	    OPTIONS(OPT_ALL_SEGMENTS) = OPTV_ALL_SEGMENTS
	  ELSE
	    CALL LIB$SIGNAL(FDQ_ERSELFIL,%VAL(1),%VAL(Status))
	    RETSTAT = ERROR
	  ENDIF

	ELSE	      ! Binary start and stop requested	
C
C	    Set the flag in OPTIONS array for select by binary data time
C
	  OPTIONS(OPT_TIME_RANGE) = OPTV_BST_TIME_RANGE
	  OPTIONS(OPT_ALL_SEGMENTS) = OPTV_ALL_SEGMENTS
	  CALL LIB$SIGNAL(FDQ_SELBST)
	ENDIF

C	  IF ( flag for select by binary start time is set
C	       i.e. the flag for no time range selected is not set ) then

	IF (OPTIONS(OPT_TIME_RANGE).NE.OPTV_NO_TIME_RANGE. AND.
	1	RETSTAT .EQ. SUCCESS) THEN

	  IF(UPM_PRESENT('JSTART')) THEN
C	    CALL UPM_GET_VALUE to get the start time
	    Start_Time = FAC_JSTART_DEFAULT
	    Status = UPM_GET_VALUE('JSTART',Txtval,Txtlen)	
	    IF (Status .NE. SS$_Normal) THEN			
	      CALL LIB$SIGNAL(FDQ_ERTIMINT,%VAL(1),%VAL(status))
	      RETSTAT = ERROR
	    ELSE
	      Start_Time(1:Txtlen) = Txtval(1:Txtlen)
	    ENDIF
	  ENDIF
	  IF(UPM_PRESENT('JSTOP')) THEN
C	    CALL UPM_GET_VALUE to get the stop time
	    Stop_Time = FAC_JSTOP_DEFAULT
	    Status = UPM_GET_VALUE('JSTOP',Txtval,Txtlen)	
	    IF (Status .NE. SS$_Normal) THEN			
	      CALL LIB$SIGNAL(FDQ_ERTIMINT,%VAL(1),%VAL(status))
	      RETSTAT = ERROR
	    ELSE
	      Stop_Time(1:Txtlen) = Txtval(1:Txtlen)
	    ENDIF
	  ENDIF
	ENDIF

C	  IF ( Return status is success and flag for no time range 
C	    selected is not set) THEN
C
C 	    Verify consistency of start and stop times If inconsistent,
C           set the flag for no time range selected

	IF (OPTIONS(OPT_TIME_RANGE) .NE. OPTV_NO_TIME_RANGE .AND.
	1	OPTIONS(OPT_TIME_RANGE) .NE. OPTV_SPEC_SEG .AND.
	1	RETSTAT.EQ.SUCCESS ) THEN
	  IF (STOP_TIME .LT. START_TIME) THEN
	    CALL LIB$SIGNAL(FDQ_BADTIMINV)
	    OPTIONS(OPT_TIME_RANGE) = OPTV_NO_TIME_RANGE
	    CALL LIB$SIGNAL(FDQ_SELNOTIM)
	  ELSE
	    TIME_RANGE = START_TIME // ';' // STOP_TIME // ';'
	    CALL LIB$SIGNAL(FDQ_SELTIMINT,%VAL(2),Start_Time,Stop_Time)
	  ENDIF
	ENDIF

        
	IF (UPM_PRESENT('SETSEC').AND.RETSTAT.EQ.SUCCESS) THEN

C	  CALL UPM_GET_LONGWORD to get the number of seconds for the set group

	  Status = UPM_GET_LONGWORD('SETSEC',SET_SEC)	
          SEC_DEF = .false. 
	  IF (Status .NE. SS$_Normal) THEN			
	    CALL LIB$SIGNAL(FDQ_ERSELSEC,%VAL(1),%VAL(Status))
	    RETSTAT = ERROR

	  ELSE
	     If ( SET_SEC .LE. 16 .AND. SET_SEC .GE. 1 ) THEN

	       CALL LIB$SIGNAL(FDQ_SELSEC,%VAL(1),%VAL(SET_SEC))

	     ELSE

	       CALL LIB$SIGNAL(FDQ_SECTOOBIG,%VAL(1),%VAL(SET_SEC))
	     ENDIF
	  ENDIF
	ELSE
C/	  Set the number of seconds for the set group to default of 16 seconds
            SEC_DEF = .true.
	    SET_SEC = DEFAULT_SEC
	ENDIF

	IF (UPM_PRESENT('ATTITUDE').AND.RETSTAT.EQ.SUCCESS) THEN

C/	  Get the type of attitude.

	  Status = UPM_GET_VALUE('ATTITUDE',att_type,Value_Size)	

	  IF (Status .NE. SS$_Normal) THEN			
	    CALL LIB$SIGNAL(FDQ_GETATTYPERR,%VAL(1),%VAL(Status))
	    RETSTAT = ERROR

	  ELSE		! Get one of the 4 types for Build 4.0
	    IF ( att_type(1:3) .eq. 'NON' ) Then
	      ATTITUDE_TYPE = fac_none
	    ELSEIF ( att_type(1:3) .eq. 'SIM' ) Then
	      ATTITUDE_TYPE = fac_simulated
	    ELSEIF  (att_type(1:3) .eq. 'FIN') Then
	      ATTITUDE_TYPE = fac_fine_with_dirbe
	    ELSEIF ( att_type(1:3) .eq. 'DEF' ) Then
	      ATTITUDE_TYPE = fac_definitive
	    ELSEIF ( att_type(1:3) .eq. 'COA' ) Then
	      ATTITUDE_TYPE = fac_coarse
	    ELSE	 
	      CALL LIB$SIGNAL(FDQ_INVALATTSEL,%VAL(1), ATT_TYPE)
	      RETSTAT = ERROR
	    ENDIF
	  ENDIF
	  IF ( RETSTAT .eq. SUCCESS ) 
	1	CALL LIB$SIGNAL(FDQ_TYPATTSEL,%vaL(1),ATT_TYPE)	  

	ELSE
C/	  Set the attitude type  to default of none.

	  ATTITUDE_TYPE = fac_none
	ENDIF


	IF (UPM_PRESENT('LIMIT').AND.RETSTAT.EQ.SUCCESS) THEN

C/	  Get the plot flag.

	  Status = UPM_GET_VALUE('LIMIT',plot_type,Value_Size)	

	  IF (Status .NE. SS$_Normal) THEN			
	    CALL LIB$SIGNAL(FDQ_GETPLOTERR,%VAL(1),%VAL(Status))
	    RETSTAT = ERROR

	  ELSE		! Get one of the values for Build 4.0
	    IF ( plot_type(1:4) .eq. 'NONE' ) Then
	      LIMIT = fac_none
	    ELSEIF ( plot_type(1:3) .eq. 'RED' ) Then
	      LIMIT = 1
	    ELSEIF ( plot_type(1:3) .eq. 'YEL' ) Then
	      LIMIT = 2
	    ELSE	 
	      CALL LIB$SIGNAL(FDQ_INVALPLOT,%VAL(1), PLOT_TYPE)
	      RETSTAT = ERROR
	    ENDIF
	    IF ( RETSTAT .eq. SUCCESS ) 
	1	CALL LIB$SIGNAL(FDQ_PLOTSEL,%vaL(1),PLOT_TYPE)	  

	  ENDIF
	  
	ELSE
C/	  Set the plot type to default of none.

	  LIMIT  = fac_none
	ENDIF

!	If Engplots are requested, get the plot device qualifier.

	If ( LIMIT .Ne. fac_none ) Then
	  Status = UPM_GET_VALUE ('PLOT_DEVICE', Txtval, Txtlen)
	  If (Status .Ne. SS$_Normal) Then
	    RETSTAT = ERROR
	    CALL LIB$SIGNAL (%val(Status))
	  Else
	    If ( Txtval(1:2) .Eq. 'PR' ) Then
	      Plot_Device = fac_lineprinter
	    Elseif ( Txtval(1:3) .Eq. 'QMS' ) Then
	      Plot_Device = fac_laserqms
	    Elseif ( Txtval(1:2) .Eq. 'TF' ) Then
	      Plot_Device = fac_lasertf
	    Endif	 
	  Endif
	Endif

!	Get the minimum offset time for the plot start and maximum offset
!	time for the plot stop.

	If ( LIMIT .Ne. fac_none ) Then
	  If ((UPM_Present('MIN_OFFTIME') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('MIN_OFFTIME') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('MIN_OFFTIME') .Eq. SS$_Normal)) Then
	    Status = UPM_GET_LONGWORD ('MIN_OFFTIME', Min_Offtime)
	    If (Status .Ne. SS$_Normal) Then
	      RETSTAT = ERROR
	      CALL LIB$SIGNAL (%val(Status))
	    Endif
	  Endif
	  If ((UPM_Present('MAX_OFFTIME') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('MAX_OFFTIME') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('MAX_OFFTIME') .Eq. SS$_Normal)) Then

	    Status = UPM_GET_LONGWORD ('MAX_OFFTIME', Max_Offtime)
	    If (Status .Ne. SS$_Normal) Then
	      RETSTAT = ERROR
	      CALL LIB$SIGNAL (%val(Status))
	    Endif
	  Endif
	Endif

	

! Get system time, make report default name, get REPORT qualifier.

	Status = Sys$Gettim( Curtime )
	IF ( Status .Eq. SS$_Normal ) THEN
	    Call CT_Binary_To_Gmt( Curtime, Curgmt )
	ELSE
	    FDQ_GET_OPTIONS = %loc(FDQ_Aberr)
	    Call Lib$Signal(ERROR, %val(1), %val(Status))
	    RETSTAT = ERROR
	ENDIF
	Report_Default = 'FDQ_' // Start_Time(1:7) // '_' // 
	1     Stop_Time(1:7) // '.REP_' // CurGMT(1:9)
	Status = UPM_PRESENT('REPORT')
	IF (Status .EQ. upm_pres) THEN
            report_def =.false.
            Report = .True.
            Status2 = upm_get_value('REPORT',Report_File,Txtlen)
	    IF (Status2 .EQ. upm_absent) THEN
		Report_File = Report_Default
	    END IF
	ELSE IF (Status .EQ. upm_defaulted) THEN
            Report_def= .true.
	    Report = .True.
	    Report_File = Report_Default
	ELSE
            Report = .False.
	END IF

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_OPTIONS = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GET_OPTIONS = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
