	INTEGER*4 FUNCTION FDQ_PROC_CUR_SET (NEW_SEGMENTS, CT_LUN,
	1	HSK_RECS, FIRST_MJ, AVG_TIME, Eng_badtime,
	1	IPDU_RELAY, STMON_SYNC,
	1	NPROC, PROC_CHAN, READ_NXT, TOLS, SCI_REC, ENG_REC,
	1	XCAL_SET, Scilim_Rec, Englim_Rec, Limflags,
	1	Idxtols, Idxflags, Grtwt,grttrans,Attitude_Type, GMT_Start, 
	1	GMT_Stop,
	1	Sc_Side, LIMIT, Plot_Device, Min_Offtime, Max_Offtime,
	1	Plot_Table, Plot_Start, Plot_Stop, ifg_count, englim,IDX_WRIT,
	1	time_range, samprate)
C/-----------------------------------------------------------------------------
C/
C/	PROGRAM NAME:
C/	  FDQ_PROC_CUR_SET
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine does the actual processing of a set of SCI records
C/	  and bracketing housekeeping records to generate an ENG record 
C/	  and possibly an IDX (Index) record.
C/	  The original SCI records are modified to the effect that their
C/	  quality flags, attitude and some instrument related fields are 
C/	  filled in.
C/
C/	AUTHOR:
C/	  EDWIN FUNG
C/	  GSFC
C/	  APRIL 1987
C/
C/	MODIFIED BY:
C/	  D. T. WARD
C/	  GSFC/STX
C/	  October 2, 1987
C/	Reason:  The DATASET_ID was included for SCI_REC.
C/
C/	MODIFIED BY:
C/	  D. T. WARD
C/	  GSFC/STX
C/	  October 20, 1987
C/	Reason:  HEC and DEC were replaced with ETR.
C/
C/	MODIFIED BY:
C/	  Shirley M. Read
C/	  STX
C/	  January 4, 1988
C/	Reason:  The routine was redesigned in accord with a new set of
C/		 quality flags for the Science Record. The new quality 
C/	         flags are designed to relfect information pertaining to 
C/	         each channel. Thus each set of 110 data quality flags must 
C/	         be checked separately for each channel. Since some flags 
C/	         relate to the attitude, the FUT_Attitude must be called
C/	         before the FDQ_Get_Qualflags. Thus the design for 
C/	         FDQ_Proc_Cur_Set was changed to process all of this 
C/		 information in a channel loop.The quality flags are all set
C/		 at this point in the science record. Thus qualflags is
C/	         no longer needed in the calling sequence.
C/		 As routines are modified for FDQ, they are being converted
C/		 to functions. The status for the function will be checked
C/	 	 to determine whether processing should continue. All error
C/		 detection will now be done through the Fut_Error condition
C/	         handler via Lib$Signal.
C/	         Also added xcal_set, scilim_rec, englim_rec and limflags
C/		 to calling sequence. 
C/
C/	MODIFIED BY:
C/	  Shirley M. Read
C/	  STX
C/	  March, 1988
C/	Reason:  The FDQ_ATTITUDE function was moved to the FUT facility so
C/	         that the FIRAS simulator could use the same routine. Thus
C/		 FUT_ATTITUDE will now be invoked by FDQ.
C/
C/	MODIFIED BY:
C/	  Shirley M. Read
C/	  STX
C/	  May, 1988
C/	Reason:    The type of attitude solution, start GMT, stop GMT and 
C/		 spacecraft side were added to the calling sequence to be 
C/		 used in functions invoked by FDQ_Proc_Cur_Set. New functions
C/		 are added to get the spacecraft side and set the limits 
C/		 enable/disable flags according to the current Limflags 
C/		 Dataset. Added idxtols and idxflags as input, plot table, 
C/		 plot start and plot stop as output.
C/	           The logging option was never fully implemented. Other
C/		 reporting functions replaced logging. Passed the report
C/	         unit number through the lun array.
C/
C/	MODIFIED BY:
C/	  H. Wang
C/	  STX
C/	  Jan., 1991
C/	Reason:    New requirements for FDQ
C/              1) Gain and fakeit data will no longer be figured by FDQ
C/              2) The calibrator resistor counts used for the conversion
C/                 of temperature GRT's will be read from FEX_AV_CLARS
C/                 ref. file
C/              3) The instrument temperature changes between bracketing
C/                 major frames will be calculated and stored in FDQ_ENG
C/              4) New quantities will be tested for the generation of
C/                 Idex records
C/              5) The attitude info. - Terrestrial radiation region location
C/                 flag will be stored in each FDQ_SDF
C/ 
CH	Change Log: 	New Format for Build 4.2,  STX 10/15/88.
CH
CH	Version 4.1.1 10/15/88, SPR 2620, Shirley M. Read, STX
CH		Disable limit checking of attitude fields for /ATT=NONE.
CH	        The calling sequence to FDQ_Set_Limflags had to be changed.
CH      Version 4.1.1 10/15/88, SPR 2622, Shirley M. Read, STX
CH              Flag values for missing conversions need changes in FDQ.
CH	        The flag values for converted fields which are out of range
CH		need the same changes. GRTs in dwell mode need a flag value
CH              also. The SWG has selected a flag value of -9999 for all cases.
CH              The limit checking algorithms need to bypass setting the quality
CH              flags if the GRT or engineering analog has the flag value.
CH		Interpolation and averaging algorithms need to be modified to
CH              include the flag value. 
CH      Version 4.1.1 10/15/88, SPR 2625, Shirley M. Read, STX
CH		FDQ must call FUT_Attitude even with the /ATT=NONE qualifier.
CH	        The attitude routine will put a pixel number of -1 in the
CH		science record as a flag value. The attitude fields in the 
CH	        science record will be initialized also by the routine. 
CH	Version 4.1.1 10/20/88, SPR 2668, Shirley M. Read, STX
CH		The addition of the flag value of -9999.0 for engineering
CH		analogs has affected the calculation of the internal reference
CH		source temperature in FDQ. The FUT_Temperature_List has been
CH		updated to include this flag and return a better value for
CH		the averaged temperatures. This routine is common to other
CH		FIRAS pipeline facilities. Thus FDQ should use this routine.
CH	Version 4.2.1 11/30/88, SPR 2906, R. Kummerer, FDQ_GET_QUALFLAGS
Ch		becomes FUT_GET_QUALFLAGS.
CH	Version 4.2.1 1/12/89, SPR 2890, R. Kummerer, FUT_TEMPERATURE GET_CONFIG
CH		related changes.
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
CH	Version 4.2.1 02/09/89, SPR 2550, Shirley M. Read, STX
CH		Provide user option of setting the time interval for triggered
CH		plots.
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH	Version 4.4.1 06/15/89, SPR 3961, Shirley M. Read, STX
CH              the attitude in the raw science record is left zero on IFGs
CH		failed for bad telemetry by FDQ. This causes problems in 
CH		quicklook programs which need the pixel number but do not
CH		check the telemetry quality or time tags.
CH	Version 4.4.1 07/21/89, SPR 3963, Shirley M. Read, STX
CH	        Attitude=NONE filled when COARSE was selected.
CH	        This can happen either because, on the one hand the Att. or
CH	        Orbit Arcv Open fails, or on the other hand, because the
CH	        UAX/UOX calls fail from lack of data.  The two cases need to
CH	        be distinguished by the FDQ program, because in the latter
CH	        case there may be a temporary absence of data which will clear
CH	        up later in the timerange; but if the Arcvs can't be opened,
CH	        then there's not going to be any data, and we can abandon the
CH	        attempt.  FUT_ATTITUDE status flag FUT_OPENATTORB signals open
CH		error conditions.
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
CH	Version 5.8,  1990 Jan 31,  SPR 5278,  F. Shuman,  STX
CH		Separate open error conditions in FUT_ATTITUDE for Attitude
CH		Archive and Orbit Archive, so that when one of these archives
CH		fails to open, we can still recover quantities computed solely
CH		from the other.
CH
CH	Version ???, 30 August 1990, Larry P. Rosen, STX
CH		Pass Report_File name and Report logical flag from call to
CH		FDQ_MAKE_IDX for standard report file name.
CH
CH	Version 7.1, 25 Oct. 1990, H.Wang, STX
CH              SPR #7501, Out 0f order time-ordered-data(engineering)
CH
CH	Version 7.1, 25 Oct. 1990, H.Wang, STX
CH              SPR # 7524, Missing data in the fdq_eng
CH
CH	Version 7.3, 12 DEc. 1990, H.Wang, STX
CH              SPR # 7855, Missing Gain data in the fdq_idx
CH
CH	Version 9.8, 4 August 1992, S. Alexander, HSTX
CH		SER # 9848, Pass sampling rate as a parameter from FDQ to
CH		FUT_GET_QUALFLAGS
C------------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  STATUS = FDQ_PROC_CUR_SET (NEW_SEGMENTS, CT_LUN, HSK_RECS, FIRST_MJ,
C/		AVG_TIME,ENG_BADTIME, 
C/              IPDU_RELAY, STMON_SYNC, NPROC, PROC_CHAN, READ_NXT,
C/		TOLS, SCI_REC, ENG_REC, XCAL_SET, SCILIM_REC, ENGLIM_REC,
C/		LIMFLAGS, IDXTOLS, IDXFLAGS, ATTITUDE_TYPE, GMT_START,
C/		GMT_STOP, SC_SIDE, PLOT, PLOT_DEVICE, MIN_OFFTIME, MAX_OFFTIME, 
C/		PLOT_TABLE, PLOT_START, PLOT_STOP, IFG_COUNT, ENGLIM,IDX_WRIT,
C/	        time_range, samprate)
C/
C/	INPUT/OUTPUT PARAMETERS:
C/	  NEW_SEGMENTS	I	L*1	Flag to indicate whether this pass
C/					is the beginning of a new segment.  It
C/					is passed down to FDQ_MAKE_IDX (a subroutine
C/					on the next level down).
C/	  CT_LUN(22)	I	I*2	COBETRIEVE logical unit numbers, the 1st
C/					4 being SCI files (RH, RL, LH, LL respectively).
C/	  HSK_RECS(2)	I	REC	Array of 2 record buffers (4 major frames)
C/					filled by routine FDQ_GET_BRK_HSK
C/	  FIRST_MJ	I	I*2	Indicates which of the 4 major frames in
C/					HSK_RECS is the "first" major frame for the
C/					current set.
C/	  AVG_TIME(2)	I	I*4	Average start time (in binary format) for
C/					the channels being processed
C/        ENG_BADTIME   I       L*1     The flag to indicate bad enG. time
C/	  IPDU_RELAY(8)	I	BYTE	IPDU relay (to be passed down to FDQ_MAKE_IDX)
C/	  STMON_SYNC(2)	I	BYTE	Status monitor sync words (to be passed
C/					down to FDQ_MAKE_IDX)
C/	  NPROC		I	I*2	Number of channels being processed in this pass
C/	  PROC_CHAN(4)	I	I*2	List of channel ID's that are being processed
C/					in this pass
C/	  READ_NXT(4)	I,O	L*1	Array of flags for determining whether another
C/					SCI record should be read in from the archives
C/					This flag is set to "true" after CT_MODIFY_PUT
C/					has been called for a particular channel.
C/	  TOLS(10)	I	I*2	List of tolerances for setting quality
C/					flags.  
C/	  SCI_REC(4)	I,O	Rec	Record buffers of raw science files for
C/					each channel
C/	  ENG_REC	I,O	Rec	Record buffer for ENG archive.
C/	  XCAL_SET      I       L*1     User option to set xcal or not.
C/	  SCILIM_REC    I       Rec     Science Limits, Red and Yellow.
C/	  ENGLIM_REC    I       Rec     Engineering Limits, Red & Yellow, High
C/					& Low.
C/	  LIMFLAGS      I       Rec     Flags to switch limit checks ON/OFF.
C/	  IDXTOLS       I       Rec     Index tolerances.
C/	  IDXFLAGS      I       Rec     Index check enable flags.
C/	  GRTWT         I       Rec     GRT switches table for FUT_TEMPERATURE_LIST.
C/	  GRTTRANS      I       Rec     Recs for fex_grttrans
C/	  ATTITUDE_TYPE I       I*4     Type of attitude requested for run.
C/	  GMT_START     I       C*14    GMT start time.
C/	  GMT_STOP      I       C*14    GMT stop time.
C/	  SC_SIDE       I       I*2     Spacecraft side selected by user.
C/					This is temporary until SC Archive is
C/				        available to FDQ.
C/	  LIMIT		I       I*2     Flag indicating if Engplots to be made
C/					and if so, trigger on yellow/red or 
C/					only on red limits exceeded. 
C/	  PLOT_DEVICE   I       I*4     Plot device: line or laser printer
C/	  MIN_OFFTIME   I       I*4     Hours offset for plot start
C/	  MAX_OFFTIME   I       I*4     Hours offset for plot stop
C/	  PLOT_TABLE(118) O     L*1     Table of flags to trigger Engplots.
C/	  PLOT_START(2) O       I*4     Plot start time.
C/	  PLOT_STOP(2)  O       I*4     Plot stop time.
C/	  ENGLIM                I*4     Engineering limits exceeded flag
C/        IDX_WRIT      O       I*4     The number of new idx records  
C/	  SAMPRATE	I	R*4	Sampling rate (I&T or Mission)
C/
C/	INPUT/OUTPUT FILES:
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC
C/	  FUT_PARAMS.TXT
C/
C/	SUBROUTINES CALLED:
C/
C/	ERROR HANDLING:
C/	  Messages sent via Lib$Signal to Fut_Error condition handler.
C/
C/	METHOD USED:
C/
C/      PDL for FDQ_PROC_CUR_SET
C/
C/	BEGIN
C/
C/	SET return status to success. The return status will be set to failure
C/	    if any call results in an error status return.
C/	CALL FDQ_GET_SC_CONFIG to determine which side of the spacecraft
C/	    is controlling FIRAS.
C/	CALL FDQ_SET_LIMFLAGS to set the Limit_On flags corresponding to
C/	    the 110 data quality flags in the science record.
C/
C/	IF ( there is valid housekeeping data ) THEN
C/
C/	  CALL FDQ_CONVERT to convert the housekeping buffers to engineering
C/	       buffers.
C/	  CALL FDQ_INTERPOLATE to interpolate the engineering buffers produced
C/	       by the bracketing housekeeping records and produce the new
C/	       engineering record. It also calculates the internal reference
C/	       source temperature for the science record.
C/	  CALL FDQ_MAKE_INDEX to construct the index record and determine if
C/	       it should be written. If so, write the record in the archive.
C/	ENDIF
C/
C/	DO for each channel in the set
C/
C/	  Load the start and stop times from the channel science record into
C/	  the engineering record.
C/
C/	  IF ( there is valid housekeeping data ) THEN
C/
C/	    Load the external calibrator position from the index record
C/	    into the engineering record.
C/	    Load the external calibrator position from the index record
C/	    into the science record.
C/	    Load the internal reference source temperature into the science
C/	    record.
C/	  ENDIF
C/
C/	  IF ( attitude is requested and status from previous call OK ) THEN
C/	    CALL FUT_ATTITUDE to get the attitude information for the 
C/	       science record.
C/	  ENDIF
C/	  IF ( there is valid housekeeping data ) THEN
C/	    CALL FUT_GET_QUALFLAGS to check the quality of the science,       
C/	       attitude and engineering related data against the red and
C/	       yellow limits and set the 110 data quality flags in the
C/	       science record. The overall quality summary flag is computed
C/             for the interferogram and a count is kept of the number of
C/	       times each individual quality flag is bad for a given segment.
C/        ENDIF
C/
C/	  Load the average interferogram time into the engineering time
C/	  fields of the science record.
C/
C/	  Set the read next flag for this channel to True.
C/	ENDDO
C/	SET function to return status
C/	RETURN
C/	END  
C/-----------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE 'CT$LIBRARY:CTUSER.INC'
	include  '(fut_params)'

	dictionary	'fdq_eng'
	dictionary	'fdq_idx'
	dictionary	'nfs_hkp'
	dictionary	'nfs_sdf'
	dictionary      'fex_scilim'
	dictionary      'fex_englim'
	dictionary      'fex_limflags'
	dictionary      'fex_idx_tols'
	dictionary      'fex_idx_flag'
	dictionary      'fex_grtrawwt'
	dictionary      'fex_grttrans'
        dictionary      'FUT_enganlg'

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_CONVERTER
	EXTERNAL 	FDQ_INTERPOLER
	EXTERNAL	FDQ_MAKEIDXER
	EXTERNAL	FDQ_ATTITUDER
	EXTERNAL	FDQ_GETQUALFLG
	EXTERNAL	FDQ_CTWRITERR
	EXTERNAL   	FDQ_GETSCCONFERR
	EXTERNAL	FDQ_SETLIMFLGERR
	EXTERNAL        FDQ_WRITEMATRIX
	EXTERNAL        FUT_NORMAL
	EXTERNAL        FUT_ABERR
	EXTERNAL        FUT_ATTOPNERR
	EXTERNAL        FUT_ORBOPNERR
	EXTERNAL        FUT_ATTORBOPNERR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			    !
!     Passed parameters     !
!			    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LOGICAL*1	IPDU_RELAY(8), NEW_SEGMENTS,
	1		READ_NXT(4), STMON_SYNC(2)
        LOGICAL*1       ENG_BADTIME
	INTEGER*2	CT_LUN(22), FIRST_MJ, NPROC, PROC_CHAN(4),
	1		TOLS(10), TELM_FLAG
	INTEGER*4	AVG_TIME(2), IFG_COUNT(4), ENGLIM,IDX_WRIT
	LOGICAL*1       XCAL_SET        ! Flag to det xcal position.
	INTEGER*4       ATTITUDE_TYPE   ! Type of attitude requested.
	CHARACTER*14    GMT_START       ! GMT start time.
	CHARACTER*14    GMT_STOP        ! GMT stop time.
	INTEGER*2	SC_SIDE	  	! Spacecraft side controlling FIRAS.
	INTEGER*2	LIMIT		! Flag for Engplots, red and yellow.
	INTEGER*4       PLOT_DEVICE     ! Plot device: line or laser printer
	INTEGER*4	MIN_OFFTIME      ! Hours to offset for plot start
	INTEGER*4	MAX_OFFTIME      ! Hours to offset for plot stop
	LOGICAL*1       PLOT_TABLE(118) ! Table of flags to trigger Engplots.
	INTEGER*4       PLOT_START(2)   ! Plot start time.
	INTEGER*4       PLOT_STOP(2)    ! Plot stop time.
        Character*30    Time_range
	REAL*4		SAMPRATE	! Sampling rate (I&T or Mission)

	RECORD /FDQ_ENG/ ENG_REC
	RECORD /NFS_HKP/ HSK_RECS(2)
	RECORD /NFS_SDF/ SCI_REC(4)
	RECORD /FEX_SCILIM/ SCILIM_REC
	RECORD /FEX_ENGLIM/ ENGLIM_REC
	RECORD /FEX_LIMFLAGS/ LIMFLAGS
	RECORD /FEX_IDX_TOLS/ IDXTOLS
	RECORD /FEX_IDX_FLAG/ IDXFLAGS
	RECORD /FEX_GRTRAWWT/ GRTWT
	RECORD /FEX_GRTTRANS/ GRTTRANS
C	RECORD /FEX_av_clars/ CALRS
C	RECORD /FEX_saa_vab/ SAA
	RECORD /Fut_enganlg/ rawsigmas, correct_temp


!!!!!!!!!!!!!!!!!!!!!!
!		     !
!    "Constants"     !
!		     !
!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*2	BAD, ETR, ENG, GOOD, IDX,CAL,
	1		HSK, RS1, RS2, RS3, RS4, RPT,
	2		SM1, SM2, SM3, SM4,
	3		SC1, SC2, SC3, SC4,
	4		ORS1, ORS2, ORS3, ORS4

	PARAMETER	(GOOD = 1)
	PARAMETER	(BAD = 2)

	PARAMETER	(RS1 = 1)
	PARAMETER	(RS2 = 2)
	PARAMETER	(RS3 = 3)
	PARAMETER	(RS4 = 4)
	PARAMETER	(HSK = 5)
	PARAMETER	(IDX = 6)
	PARAMETER	(ENG = 7)
	PARAMETER	(ETR = 8)
	PARAMETER	(CAL = 9)
	PARAMETER       (RPT = 10)  ! Ct_Lun 10 is reserved for RMS report file
	PARAMETER	(SM1 = 11)
	PARAMETER	(SM2 = 12)
	PARAMETER	(SM3 = 13)
	PARAMETER	(SM4 = 14)
	PARAMETER	(SC1 = 15)
	PARAMETER	(SC2 = 16)
	PARAMETER	(SC3 = 17)
	PARAMETER	(SC4 = 18)
	PARAMETER	(ORS1 = 19)
	PARAMETER	(ORS2 = 20)
	PARAMETER	(ORS3 = 21)
	PARAMETER	(ORS4 = 22)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         !
!     Functions           !
!                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*4	FDQ_CONVERT,			! convert Hsk buffers
	1		FDQ_INTERPOLATE,		! interpolate Eng buffer
	1		FDQ_MAKE_IDX,			! make index record
	1		FUT_ATTITUDE,			! attitude info function
	1		FUT_GET_QUALFLAGS,		! get data quality flags
	1		FDQ_GET_SC_CONFIG,	        ! get SC configuration
	1		FDQ_SET_LIMFLAGS,		! set limit flags
	1		FUT_TEMP_CORRECT               ! correct GRT temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	BYTE		HSK_ANALOG(52)
	INTEGER*2	CHAN, CNT, CT_STAT(20), I, ISTAT,
	1		SCI_LUN(4)
	LOGICAL*1       LIMIT_ON(110)		    ! Limit check enable flags	
	INTEGER*4	R_STATUS,	            ! return status
	1	        RETSTAT
	INTEGER*4       STATUS				! I*4 dummy status var
	INTEGER*4	SUCCESS / 1 /, ERROR / 2 /	! value of return status
	INTEGER*4       ZERO / 0 /,combswitch/0/
	REAL*4		ical_avg, ical_hi_avg, ical_lo_avg,
	1		int_ref_src_temp
	LOGICAL*1	end_segment
        LOGICAL*4       singlifg/.true./
        Logical*1       First_chan
	INTEGER*4       att_type                    ! Type of attitude = none
	REAL*4		FV / -9999.0 /
	REAL*4		TEMP(10),outsigmas(10)		    ! Temperature list
        INteger*2       gain_val(0:7)/1,3,10,30,100,300,1000,3000/
	RECORD /FDQ_ENG/ ENG_BUFFS(2)
	RECORD /FDQ_IDX/ IDX_REC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to Success.

	RETSTAT = SUCCESS

	CALL LIB$MOVC3 (4*2, CT_LUN(ORS1), SCI_LUN(1))	! This is to achieve some kind
							! of independence of the indices
							! of SCI_LUN because they are
							! being used in a DO-loop

!	Increment the IFG count before any processing is done; counts
!	must keep in step with the IFG in the segment, as IFGs may fail
!	processing in this routine before the science record is modified
!	in the archive.

	do cnt = 1, nproc
	  chan = proc_chan(cnt)
	  ifg_count(chan) = ifg_count(chan) + 1
	end do

!	If new segment is starting, set the end_segment flag to print the
!	quality flags summary.

	If ( New_Segments ) then
	    end_segment = .true.
	Else
	    end_segment = .false.
	Endif

! 	Get the spacecraft configuration.

	R_STATUS = FDQ_GET_SC_CONFIG (SC_SIDE)
	IF ( R_STATUS .NE. %loc(FDQ_NORMAL)) THEN
	    if ( R_Status .eq. %loc(FDQ_ABERR)) R_status = zero
	    RETSTAT = ERROR
	    CALL LIB$SIGNAL( FDQ_GETSCCONFERR, %VAL(1), %VAL(R_STATUS))
	ENDIF

!	Call FDQ_SET_LIMFLAGS to enable/disable limit checking flags.

	IF ( RETSTAT .EQ. SUCCESS ) THEN
	  R_STATUS = FDQ_SET_LIMFLAGS (LIMFLAGS, SC_SIDE, 
	1	     ATTITUDE_TYPE, LIMIT_ON)
	  IF ( R_STATUS .NE. %loc(FDQ_NORMAL)) THEN
	    if ( R_Status .eq. %loc(FDQ_ABERR)) R_status = zero
	    RETSTAT = ERROR
	    CALL LIB$SIGNAL( FDQ_SETLIMFLGERR, %VAL(1), %VAL(R_STATUS))
	  ENDIF
	ENDIF

	IF ((eng_rec.en_tail.hskp_flag .EQ. 0) .and. 
	1  (retstat .eq. success)) THEN   ! This means we have valid HSK data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								    !
!     In the following 2 calls, FDQ_CONVERT will fill the 2 ENG	    !
!     buffers (ENG_BUFFS), and FDQ_INTERPOLATE will fill ENG_REC.   !
!								    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  R_STATUS = FDQ_CONVERT (HSK_RECS, FIRST_MJ, CT_LUN(cal),Time_range,
	1	ENG_BUFFS)
	  IF ( R_STATUS .NE. %loc(FDQ_NORMAL)) THEN
	    if ( R_Status .eq. %loc(FDQ_ABERR)) R_status = zero
	    RETSTAT = ERROR
	    CALL LIB$SIGNAL( FDQ_CONVERTER, %VAL(1), %VAL(R_STATUS))
	  ENDIF
 
	  IF ( RETSTAT.EQ.SUCCESS ) THEN
	    R_STATUS = FDQ_INTERPOLATE (HSK_RECS, FIRST_MJ, AVG_TIME, 
	1	ENG_BUFFS, ENG_REC, HSK_ANALOG, XCAL_SET,TELM_FLAG,
	1	grtwt,grttrans)
	    IF ( R_STATUS .NE. %loc(FDQ_NORMAL)) THEN
	      if ( R_Status .eq. %loc(FDQ_ABERR)) R_status = zero
	      RETSTAT = ERROR
	      CALL LIB$SIGNAL( FDQ_INTERPOLER, %VAL(1), %VAL(R_STATUS))
	    ENDIF
	  ENDIF	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									    !
!     Next calculate the Internal Reference Source Temperature from the     !
!     4 versions of the same quantity (A-high, A-low, B-high and B-low	    !
!									    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  IF ( RETSTAT.EQ.SUCCESS ) THEN 	! Invoke FUT routine 10/21/88
            R_STATUS  = FUT_TEMP_CORRECT(eng_rec.en_analog, correct_temp)
C   
C  SPR 8899, The IREF_TEMP value is no longer being calaulated by FDQ
C    
C	    R_STATUS = FUT_TEMP_LIST (eng_rec.en_analog, 
C	1	rawsigmas,grtswt,combswitch,singlifg, temp,outsigmas)
C
C	      int_ref_src_temp = temp(2)	! Ical temperature
	  ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						                       !
!     Invoke FDQ_MAKE_IDX to construct the index record and determine  !
!     if it should be written. If so, write the record in the archive. !
!						                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  IF ( RETSTAT.EQ.SUCCESS ) THEN
	    do	cnt = 1, nproc
	      chan = proc_chan(cnt)
	      Eng_rec.chan(chan).fakeit = sci_rec(chan).dq_data.fake
            enddo          
            If (.not. eng_badtime) then
	    R_STATUS = FDQ_MAKE_IDX (new_segments, eng_rec, sci_rec,
	1	 hsk_analog, ipdu_relay, stmon_sync, ct_lun(idx), NPROC,
	1	 PROC_CHAN, IDXTOLS, IDXFLAGS, idx_rec,telm_flag,IDX_WRIT)
	    IF ( R_STATUS .NE. %loc(FDQ_NORMAL)) THEN
	      if ( R_Status .eq. %loc(FDQ_ABERR)) R_status = zero
	      RETSTAT = ERROR
	      CALL LIB$SIGNAL( FDQ_MAKEIDXER, %VAL(1), %VAL(R_STATUS))
	    ENDIF
            Endif  
	  ENDIF
	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								     !
!     Now we have to do some processing on those Science records     !
!     that are marked "being processed".  The following DO loop	     !
!     handles this						     !
!								     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF ( RETSTAT.EQ.SUCCESS ) THEN
!
!     Science Mode Set-up (Ext. Cal. Position)
!
	   eng_rec.en_xcal.pos(1) = idx_rec.xcal.pos(1)
	   eng_rec.en_xcal.pos(2) = idx_rec.xcal.pos(2)
!     Hot spot command (on/off) store into Eng record
!
          eng_rec.en_stat.hot_spot_cmd(1) = idx_rec.ipdu_stat(1).
	1	                            hot_spot_heater
          eng_rec.en_stat.hot_spot_cmd(2) = idx_rec.ipdu_stat(2).
	1	                            hot_spot_heater
!
	  do	cnt = 1, nproc
	   IF ( RETSTAT.EQ.SUCCESS) THEN
	    chan = proc_chan(cnt)
            eng_rec.en_head.sci_time(chan).bin_time(1) = 
     &                                   sci_rec(chan).ct_head.time(1)
            eng_rec.en_head.sci_time(chan).bin_time(2) = 
     &                                   sci_rec(chan).ct_head.time(2)
	    eng_rec.en_head.ifg_no(chan) = ifg_count(chan)

	    IF (eng_rec.en_tail.hskp_flag .EQ. 0) THEN	! This means we have valid HSK data
!
!     Science Mode Set-up (Ext. Cal. Position)
!
!     Fill in Channel-Specific fields in ENG record
!     lib$movc3 used because 1st 6 fields are byte
!     while last is int*2
!
	      call lib$movc3(8,idx_rec.chan(chan),eng_rec.chan(chan))

              if(chan.eq.1 .or. chan.eq.2)then
                sci_rec(chan).dq_data.xcal_pos = idx_rec.xcal.pos(1)
              else
                sci_rec(chan).dq_data.xcal_pos = idx_rec.xcal.pos(2)
              endif

!     Fill in Internal Reference Temperature field
C   
C  SPR 8899, The IREF_TEMP value is no longer being calaulated by FDQ
C    
C
C              sci_rec(chan).dq_data.iref_temp = int_ref_src_temp
	    endif		! if (eng_rec.en_tail.hskp_flag .EQ. 0) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						       !
!     Invoke FUT_ATTITUDE to fill in attitude info     !
!						       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    R_STATUS = FUT_ATTITUDE ( SCI_REC(CHAN), ATTITUDE_TYPE,
	1                             GMT_START, GMT_STOP )

!  FUT_ATTITUDE now handles attit/orbit arcv open failures on its own,
!    computing whatever it is able.  There are a couple quantities it can
!    fill even if neither arcv opens -- those that depend only on positions
!    of Solar System objects (Moon phase, Sun-Moon distance).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						                        !
!     Invoke FUT_GET_QUALFLAGS to check the quality of the science,     !
!     attitude and engineering related data against the red and yellow  !
!     limits and set the 110 data quality flags in the science record.  !
!						                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    IF ( RETSTAT.EQ.SUCCESS ) THEN
	      IF (eng_rec.en_tail.hskp_flag .EQ. 0 .and. 
	1	  (.not. eng_badtime)) THEN  ! This means we 
          			             ! have valid HSK data
		If(Cnt .eq. 1) then
                  First_chan = .true.
                Else
                  First_chan = .false.
                Endif
	        R_STATUS = FUT_GET_QUALFLAGS (First_chan,
	1	END_SEGMENT, FIRST_MJ, SCI_REC(CHAN),    
	1	ENG_REC, HSK_RECS, SCILIM_REC, ENGLIM_REC, LIMFLAGS,
	1	LIMIT_ON, SC_SIDE, LIMIT, PLOT_TABLE, PLOT_START, PLOT_STOP, 
	1	PLOT_DEVICE, MIN_OFFTIME, MAX_OFFTIME, CT_LUN(10), ENGLIM,
	1	SAMPRATE)
	        IF ( R_STATUS .NE. %loc(FUT_NORMAL)) THEN
	          if ( R_Status .eq. %loc(FUT_ABERR)) R_status = zero
	          RETSTAT = ERROR
	          CALL LIB$SIGNAL(FDQ_GETQUALFLG,%VAL(1),%VAL(R_STATUS))   
	        ENDIF
	      ELSE	! Set the data quality summary flag to value indicating
			! there is no housekeeping data available.
	        sci_rec(chan).dq_data.data_quality(110) = fac_no_hskp
 	      ENDIF	! Housekeeping data is available.
	    ENDIF	! Retstat is success
!
!     Back to the associated science record
!
	    IF (RETSTAT.EQ.SUCCESS ) THEN
              sci_rec(chan).dq_data.eng_time(1) = Eng_rec.ct_head.time(1)
              sci_rec(chan).dq_data.eng_time(2) = Eng_rec.ct_head.time(2)
              sci_rec(chan).dq_data.ifg_no = ifg_count(chan)

!     Fill in the Data Quality Summary Flag in FDQ_ENG record.

	      eng_rec.en_head.dq_summary_flag(chan) = 
	1	sci_rec(chan).dq_data.data_quality(110)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							    !
!     Write out modified science record back to archive     !
!							    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      IF (CHAN .EQ. 1) 
	1	SCI_REC(CHAN).CT_HEAD.DATASET_ID=FAC_FDQ_RH
	      IF (CHAN .EQ. 2) 
	1	SCI_REC(CHAN).CT_HEAD.DATASET_ID=FAC_FDQ_RL
	      IF (CHAN .EQ. 3) 
	1	SCI_REC(CHAN).CT_HEAD.DATASET_ID=FAC_FDQ_LH
	      IF (CHAN .EQ. 4) 	
	1	SCI_REC(CHAN).CT_HEAD.DATASET_ID=FAC_FDQ_LL
              call lib$movc3(1,sci_rec(chan).sci_head.sc_head1a,
	1	  Eng_rec.chan(chan).up_sci_mode) 
              call lib$movc3(1,sci_rec(chan).sci_head.sc_head9,
	1	Eng_rec.chan(chan).up_adds_per_group) 
              call lib$movc3(1,sci_rec(chan).sci_head.sc_head11,
	1     Eng_rec.chan(chan).up_swps_per_ifg) 
              call lib$movc3(1,sci_rec(chan).sci_head.mtm_speed,
	1	Eng_rec.chan(chan).xmit_mtm_speed) 
              call lib$movc3(1,sci_rec(chan).sci_head.mtm_length,
	1	Eng_rec.chan(chan).xmit_mtm_len)
              If (sci_rec(chan).sci_head.gain .ge. 0) then 
              Eng_rec.chan(chan).sci_gain=gain_val(sci_rec(chan).sci_head.gain)
              else
              Eng_rec.chan(chan).sci_gain=-1
              endif
              Eng_rec.en_head.att_summary_flag(chan) = sci_rec(chan).
	1	dq_data.data_quality(109)
              Eng_rec.ct_head.orbit=sci_rec(chan).ct_head.orbit
	      CALL CT_WRITE_ARCV(,SCI_LUN(CHAN),SCI_REC(CHAN),CT_STAT)
	      IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	        RETSTAT = ERROR
	        STATUS = CT_STAT(1)
	        CALL LIB$SIGNAL( FDQ_CTWRITERR, %VAL(2), %VAL(STATUS),
	1	  DATASETS(CHAN+18))
	      ENDIF

	      read_nxt(chan) = .true.


	    ENDIF
	   ENDIF	! Retstat is success for previous loop

	  enddo		! do cnt = 1, nproc
	ENDIF

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_PROC_CUR_SET = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_PROC_CUR_SET = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
