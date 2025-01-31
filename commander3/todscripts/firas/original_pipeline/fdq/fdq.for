	PROGRAM FDQ
C/
C/	PROGRAM NAME:
C/	  FDQ	(FDQ_MAIN, DATA QUALIFY, or DQ)
C/
C/	PROGRAM DESCRIPTION:
C/	  This is the second of FIRAS' scientific "off-line"
C/	  processing program.  It will process each record in the 4 FIRAS Science
C/	  Files (All being FIRAS raw archive files) that falls within a
C/	  given time interval that is selected by the user.  The processing
C/	  includes the following:
C/		(1) Align the Science records with the corresponding
C/		    Housekeeping records (the other raw archive data set);
C/		(2) Convert the housekeeping data from counts to
C/		    engineering units (i.e., temperatures, voltages,
C/		    currents, etc.) and write out a corresponding
C/		    Engineering record, another FIRAS archive data set);
C/		(3) Check the Science record header for consistency of modes
C/		    and data quality and write back to the Science record
C/		    the quality flags thus calculated);
C/		(4) Check limits on housekeeping data and report any
C/		    out-of-limits condition to a report file as
C/		    well as setting appropriate quality flags mentioned
C/		    in item (3) above (this capability was not implemented
C/		    in either ITD1 or ITD2);
C/		(5) Compile engineering statistics on the hourly, daily
C/		    and mission statistics and write it to the Engineering
C/		    Statistics file, another FIRAS archive data set (this capability
C/		    was not available in the first ITD);
C/		(6) To write out the index record.  This function was performed
C/		    by the FIRAS stripper in ITD 1, but has been moved to this
C/		    piece of software for ITD 2 (and beyond) to alleviate the
C/		    real time requirements of the FIRAS stripper.
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  February 13, 1986
C/
C/	MODIFIED BY:
C/	  Edwin H. Fung
C/	  GSFC
C/	  April, 1987
C/	  REASON:	There is a new requirement for DATA QUALIFY for the
C/			3rd ITD:  That it should be relatively automated as
C/			compared to the previous versions.  Formerly the user
C/			is required to query the FIRAS archive catalog (using CATS)
C/			to get more or less "exact" start and stop times for
C/			the 4 Science files to use as input to DQ.  A human
C/			error in giving this start and stop times will cause
C/			difficulties in DQ processing.  This problem is especially
C/			acute because running DQ is essentially an irreversible
C/			process:  Science archives can only be modified once by
C/			DQ.  Aside from the ramifications of the possibility of
C/			human errors, it is also cumbersome for the users to
C/			manually get the start and stop times (via CATS) in
C/			order to run DQ.
C/			To address this requirement it is imperative that the
C/			former human operation of querying the catalog and
C/			picking out the start and stop times be done automatically
C/			by the software.  This is in effect one of the main
C/			reasons why we are having this new version of DQ.  In
C/			this new version it is no longer necessary for the
C/			human operator to type in the start and end times; DQ
C/			will read the catalog itself for the necessary information
C/			to start processing Science archive segments.
C/			Another feature changed in this release is that DQ will
C/			now write out one segment of ENG (engineering) and RSI
C/			(index) segment per set of SCI segments, as opposed to
C/			one segments for ALL SCI segments processed.  In other
C/			words, "time interval" is no longer the unit of archive
C/			for processing, but instead a "segment" is.   This is
C/			more in line with the COBETRIEVE philosophy of only
C/			modifying entire segments at a time (no partial segment
C/			modification) and is also inherently more elegant.
C/
C/      Modified by:
C/          J. W. Durachta
C/          August 28,1987
C/              ARC
C/            Reason:  First_time flag added to to correct for problem in get_
C/                     brk_hsk when first science record fails to have house-
C/                     keeping data.
C/
C/	Modified by:
C/		D.T.Ward
C/		October 2, 1987
C/		GSFC/STX
C/		Reason: DATASET_ID was modified for the eng and index records
C/			to CTU_$FIR_EDF & _IDX
C/
C/	Modified by:
C/		D.T.Ward
C/		October 20, 1987
C/		GSFC/STX
C/		Reason:  HES and DES were replace by ETR (Engineering TRends)
C/			 since only orbital statistics are now required.
C/
C/	Modified by:
C/		Shirley M. Read
C/		December 1987
C/		STX
C/		Reason: Redesign of error signaling process by establishing
C/			Fut_Error as the condition handler. A message file
C/			for FDQ will include all informational, warning and
C/			error messages to be output via Lib$Signal and the
C/			Fut_Error function. Error checking will be added for
C/			all functions and subroutine calls. All of the FDQ
C/			subroutines are being changed to functions. Success
C/			on completion of FDQ will also be signalled.
C/		        FDQ and all invoked FDQ functions will be redesigned
C/			to continue processing only if the status is 'success'.
C/			Since the Fut_Error function will also signal untrapped
C/			Fortran/System errors and return control to FDQ, the
C/		        existing FDQ_Exit will no longer be needed.
C/			Some calling sequences have been modified or corrected.
C/	                Write statement have been added and modified in
C/			prepraration for interface to FUT message utility in a
C/			later build.
C/
C/	Modified by:
C/		Shirley M. Read
C/		May 1988
C/		STX
C/		Reason: Modification of FDQ design to archive conversion and
C/			limits datasets and load them into buffers for FDQ
C/		        via a COBETRIEVE 'Get Configuration' routine to be
C/			called before processing every set of science records.
C/			The function would get all datasets needed by FDQ with
C/		        a time tag which covered the science set. For Build 4.0
C/		        only the FDQ_Load_Lim will be removed and the datasets
C/			will be read via the new routine, CCT_Get_Config_TOD.
C/			Modified calling sequences for some functions to pass
C/			parameters needed for FDQ enhancements.
C/			Added the capability to obtain real attitude in
C/			addition to the current simulated attitude. The user
C/			may select simulated, none or a type of real attitude.
C/			Added to capability to generate command line text
C/			files which will trigger Engplots of engineering
C/			fields which failed the limits check. Users may
C/			select switches to turn off Engplots, trigger on
C/			only exceeding red limits or trigger on exceeding
C/			either red or yellow limits. Added call to write
C/		        plot files, FDQ_Write_Engplot.
C/			Removed the logging function, which was never fully
C/			implemented, and replaced it with a report file.
C/
C/	Modified by:
C/		R. Kummerer
C/		June 27, 1988
C/		STX
C/		Reason: Interface change to FDQ_CLOSE_ARCV for Science
C/			Catalog / Matrix CT_LUN argument.
C/
C/	Modified by:
C/		Shirley M. Read
C/		July 1988
C/		STX
C/		Reason: The CCT_GET_CONFIG call was
C/		        exceeding the complete data time range for obtaining
C/			the complete set of configuration files instead of a
C/			time within the time tag of the dataset. To be safe,
C/			FDQ will determine the minimum and maximum data times
C/			from the catalog records of all science and housekeeping
C/			files to be read.
C/
C/	Modified by:
C/		Shirley M. Read
C/		August 1988
C/		STX
C/		Reason: The indexed configuration, time tagged datasets for
C/		        engineering conversions were entered in COBETRIEVE.
C/			FDQ has been modified to access these datasets using
C/			a new CCT_Get_Config routine for indexed datasets.
C/			This routine returns the unit number for the user
C/			routine to read the current time tagged segment with
C/			keyed access.
C/
CH	Version 4.1.1 10/15/88, SPR 2620, Shirley M. Read, STX
CH		Disable limit checking of attitude fields for /ATT=NONE.
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
CH      Version 4.1.1 10/15/88, SPR 2625, Shirley M. Read, STX
CH		FDQ must call FUT_Attitude even with the /ATT=NONE qualifier.
CH	        The attitude routine will put a pixel number of -1 in the
CH		science record as a flag value. The attitude fields in the
CH	        science record will be initialized also by the routine.
CH      Version 4.1.1 10/19/88, SPR 2657, Shirley M. Read, STX
CH		FDQ Convert needs to convert the new LMACs. The conversion
CH	        coefficients have now been defined and the include text file
CH	        has been expanded to include the LMACs.
CH      Version 4.1.1 10/20/88, SPR 2665, Shirley M. Read, STX
CH		FDQ needs the L3 FDQ_Save_Cat_Info to get common file versions.
CH		Processing of FDQ requires that the housekeeping file and
CH		four science files have the same file extension and version
CH	        numbers. This constraint occured when the new COBETRIEVE
CH	        allowed multiple version numbers for the same data.
CH      Version 4.1.1 10/20/88, SPR 2667, Shirley M. Read, STX
CH		FDQ has an error in counting the number of quality flags set
CH		for exceeding limits or encounters other bad data. The bit
CH		flags can have a negative value due to the msb being set.
CH	Version 4.1.1 10/20/88, SPR 2668, Shirley M. Read, STX
CH		The addition of the flag value of -9999.0 for engineering
CH		analogs has affected the calculation of the internal reference
CH		source temperature in FDQ. The FUT_Temperature_List has been
CH		updated to include this flag and return a better value for
CH		the averaged temperatures. This routine is common to other
CH		FIRAS pipeline facilities. Thus FDQ should use this routine.
CH	Version 4.2.1 11/30/88, SPR 2906, R. Kummerer, FDQ_QUAL_SUMMARY
CH		becomes FUT_QUAL_SUMMARY, FDQ_WRITE_ENGPLOT becomes
CH		FUT_WRITE_ENGPLOT.
CH
ch	version 4.2.1 12/01/88, ser 2379, J.T.Bonnell, STX@GSFC
ch		This program was modified to refer to the
ch		new firas archive logical names in response
ch		to SER 2379 (csdr$firas_in, _out, _raw,
ch		_ref, _uref, and _cal).
CH	Version 4.2.1 1/12/89, SPR 2890, R. Kummerer, FUT_TEMPERATURE_LIST
CH		GET_CONFIG related changes.
CH	Version 4.2.1 01/19/89, SPR 3150, Shirley M. Read, STX
CH		FDQ needs to access coarse attitude for the CSDR I&T
CH		Acceptance Test for Quick Look on January 19, 1989.
CH		The coarse and predicted options were not implemented
CH		at the time the real attitude access was implemented
CH		in FDQ. The keywords "COARSE" and "PREDICTED" must
CH		be picked up from the command line and the appropriate
CH		FAC Params set accordingly.
CH	Version 4.2.1 02/09/89, SPR 2700, Shirley M. Read, STX
CH		FDQ falls over for out of order and bad time ranges.The FIRAS
CH		Stripper now writes the science records using the telemetry
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
CH	Version 4.2.1 02/09/89, SPR 3062, Shirley M. Read, STX
CH		Moon contaminated IFGs failed by data quality checking. FES and
CH		FCI will interpret the data as being unusable. It is needed for
CH		moon modeling. The moon flag bit will not be counted in the
CH		data quality summary flag.
CH	Version 4.2.1 02/09/89, SPR 3139, Shirley M. Read, STX
CH		FDQ does not weed out the bad data flag value of -9999.0 when
CH		it forms the engineering trends records. Modules affected may be
CH		FDQ_Init_Stats, FDQ_Gen_Stats and FDQ_Form_Stats. In computing
CH		the statistics, the flag value should not be included in the
CH		number, sum or sum of squares.
CH	Version 4.2.1 02/09/89, SPR 2550, Shirley M. Read, STX
CH		Provide user option of setting the time interval for triggered
CH		plots.
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH	Version 4.2.1 02/09/89, SPR 2955, Shirley M. Read, STX
CH		FDQ should be given the correct COBE orbital period which is
CH		used to form engineering trends. Make the orbital period a
CH		command line option.
CH	Version 4.2.1 02/08/89, SPR 3297, Shirley M. Read, STX
CH		An attitude data quality summary flag should maintained apart
CH		from the existing over-all data quality summary flag. It will
CH		be assigned, just as the existing summary flag, a quality
CH		level dependent on the number of yellow and red limit
CH		violations of the attitude quantities checked. FDQ will fill
CH		the attitude summary flag.
CH	Version 4.2.1 03/08/89, SPR 3384, Shirley M. Read, STX
CH	        The start and stop strings need to be padded with standard
CH		CT type padding for subsequent CT calls. CT_GMT_to_Binary
CH		will not take the semicolon in the string.
CH	Version 4.4.1 07/21/89, SPR 4085, Shirley M. Read, STX
CH              FDQ should call the new FUT routine to get the orbital period.
CH		After the command line is checked for the presence of the
CH		orbital period qualifier and the stop time for the data is
CH		determined, if the period is still zero, call the FUT routine
CH		to read the NOR_ORB_HRD file from the CSDR$Space_Archive.
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
CH	Version 4.4.1 08/02/89, SER 3306, Q. C. Chung, STX
CH              Provide version number to track software update.
CH	Version 4.4.1 08/03/89, SPR 4251, Q. C. Chung, STX
CH              Exit with status in the log.
CH	Version 4.4.1 08/21/89, SER 4210, R. Kummerer STX
CH		Prevent overlaps in raw science segments. The FIRAS
CH		Project Pipeline has been designed with the assumption
CH		that the raw science data segments would not overlap
CH		in time. The EDIT data will not, but the PLAYBACK data
CH		will routinely contain about an hour overlap. The
CH		playback data may have to be used until 9 months into
CH		the mission. Thus in order to prevent problems and
CH		avoid processing IFGs in the overlap twice through FES
CH		and FCI, the segment overlaps must be removed. The FPP
CH		FIRAS Preprocessor is the most logical place to remove
CH		the overlaps before writing its intermediate science
CH		datasets for the four channels. A previous CCR was
CH		passed to change the modify feature in the FIRAS
CH		pipeline to a COBETRIEVE write. This required FPP and
CH		FDQ to write new science datasets. Since the removal
CH		of the overlap would be done before the segment catalog
Ch		is written (with segment start and stop times and IFG
CH		counts), the functions of SEGCTL would have to be done
CH		also in FPP. There is also an SER for a new facility to
CH		check for data gaps and reset the gap flag in the
CH		tracker records. This could be done also in FPP.
CH		It includes the following:
CH
CH			1. FPP must remove the segment overlaps.
CH
CH			2. FPP must now perform the functions
CH			   originally performed by SEGCTL. This
CH			   includes the writing and maintenance
CH			   of the segment catalog.
CH
CH			3. FDQ must query housekeeping.
CH
CH	Version 4.4.2 10/10/89, SPR 4705, R. Kummerer STX
CH		Exclude bad engineering data from the trends analysis.
CH		Bad engineering includes (1) engineering records produced
CH		when there were no housekeeping records at the science
CH		time and (2) engineering records produced when there were
CH		telemetry errors.  These conditions can be determined by
CH		checking the science data quality.  The engineering record
CH		is excluded from trending if any one of the science channels
CH		encounter either of the above conditions.
CH
CH	Version 4.4.3 11/14/89, SPR 5032, R. Kummerer STX
CH		Skip processing raw science segments from missing channels.
CH	Version 4.4.04 11/29/89, SPR 5207, R. Kummerer STX
CH		Skip closes on missing channel segments.
CH	Version 5.0 12/11/89, SPR 5313, R. Kummerer STX
CH		Perform appropriate raw science and IFG tracker archive closes.
CH	Version 5.2 12/29/89 SPR 5536, R. Kummerer, STX
CH		/NOTRACK erroneously selects IFG tracking.
CH	Version 5.5 1/25/90 SPR 5538, H. Wang, STX
CH		FDQ Does appear to be recognizing the bad telemetry.
CH
CH	Version 5.8 3/7/90 SPR 6297, H. Wang, STX
CH		FUT_get_qualflags routine do not need to do limit checking
CH              if the telemetry bad flag has been set
CH
CH	Version 6.1 5/15/90 SPR 6726, H. Wang, STX
CH		FDQ must support talaris 1590t printer
CH
CH	Version 6.7 8/7/90 SER 4171, Larry P. Rosen, STX
CH		Standardize report file name:
CH			Report option added to command line.
CH			Report file name and run time gotten in FDQ_GET_OPTIONS.
CH			FDQ_GET_OPTIONS is called before report file opened.
CH
CH	Version 6.7 8/20/90 SPR 7148,7129, H. Wang, STX
CH		FDQ plot command file need unique name.
CH              FDQ inapproprite error message signaled.
CH
CH	Version ??? 8/30/90 SER 4171, Larry P. Rosen, STX
CH		Standardize report file name:
CH			Report_File name gets passed to FDQ_NEW_IDX for
CH			constructing its report name.
CH
CH	Version 6.9 10/1/90 SPR 7495, H. Wang, STX
CH		FUT_GET_QUALFLAGS, did not check HKP telemtry quality 
CH              porperlt when major fram pointer are 3 or 4.
CH
CH	Version 7.1 10/25/90 SPR 7501, H. Wang, STX
CH		Out of order time-ordered-data (engineering)
CH
CH	Version 7.1 10/25/90 SPR 7524, H. Wang, STX
CH		Missing data in the FDQ_ENG record
CH
CH	Version 7.1 10/31/90 SPR 7627, H. Wang, STX
CH	 	Checks on telemetry format
CH
CH	Version 7.1 11/12/90,  L. Rosen, STX
CH		Succesful Completion if no science files found.  Not abort.
CH 
CH	Version 7.3 12/5/90 SPR 7829, H. Wang, STX
CH	 	Correct the POWER STATUS data in FDQ_IDX record
CH
CH	Version 7.3 12/11/90 SPR 7852, H. Wang, STX
CH	 	Glitch channel preamp
CH
CH	Version 7.3 12/12/90 SPR 7855, H. Wang, STX
CH	 	Missing gain data in fdq_idx
CH
CH     Modified by: Harte Wang, STX, 1/29/91
CH        Reason: New requirements for FDQ.
CH                1)FDQ new requirments for reporting:
CH                  * Account invoking facility
CH                  * Banner announcement if red engineering limits exceeded
CH                  * Full command line with defaults
CH                  * Logical name translation
CH                  * Error message
CH                2)FDQ requirements for general processing
CH                  * The Calibration resistor counts used for the conversion
CH                    of temperature GRTs will be read from FEX_AV_CALRS
CH                    ref. file
CH                  * The instrument temperature changes between bracketing
CH                    major frames will be calaulated and stored in engin.
CH                    files(FDQ_ENG)
CH                3) New requirements for idexing
CH                   * The DQ_IDX_REPORT is no longer required
CH                   * New quantities will be tested for the generation of
CH                     index records:
CH                      Bolometer CMD bias
CH                      IPDU relays
CH                      LVDT status
CH                      MTM CAL MOTOR
CH                 4) New Requirements for attitude
CH                     New attitude qualifiers
CH                      Terrestrial radiation region location flag
CH                 5) New requirements for Gain and fakeit
CH                      FDQ will no longer figure out the gains and fakeit
CH                      (FPP do it)
CH                 6) Tracking file no longer required
CH     
CH	Version 9.8, 4 August 1992, Steven Alexander, HSTX
CH		SER 9848, Retrieve sampling rate from reference dataset
CH		FEX_SAMPRATE and pass it to FUT_GET_QUALFLAGS via 
CH		FDQ_PROC_CUR_SET
C-------------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  This is a main program.
C/
C/	INPUT/OUTPUT PARAMETERS:
C/	  NONE
C/
C/	INPUT FILES:
C/	  FIRAS HOUSEKEEPING (ARCHIVE)
C/	  FIRAS SCIENCE FILE - RH (ARCHIVE)
C/	  FIRAS SCIENCE FILE - RL (ARCHIVE)
C/	  FIRAS SCIENCE FILE - LH (ARCHIVE)
C/	  FIRAS SCIENCE FILE - LL (ARCHIVE)
C/	  ATTITUDE FILE (ARCHIVE)
C/	  Curves and GRT files created by the IGSE database definition software
C/	  File for index tolerances
C/	  File for turning off certain quantities for indicing purposes
C/	  SCIENCE LIMITS FILE (ARCHIVE)
C/	  ENGINEERING LIMITS FILE (ARCHIVE)
C?	  LIMITS ENABLE/DISABLE FILE (ARCHIVE)
C/
C/	OUTPUT FILES:
C/	  FIRAS ENGINEERING DATA FILE (ARCHIVE)
C/	  FIRAS RAW INDEX FILE (ARCHIVE)
C/	  FIRAS ENGINEERING STATISTICS FILE (ARCHIVE)
C/	  FIRAS SCIENCE CATALOG FILE (ARCHIVE) (One for each channel)
C/	  FIRAS SCIENCE MATRIX FILE (ARCHIVE) (One for each channel)
C/	  ASCII LIMITS ALARM FILE
C/	  ASCII CONVERSIONS FILE
C/	  INDEX REPORT FILE (DQ_IDX_yydddhh_yydddhh.REP_yydddhhmm)
C/	  FDQ_YYDDDHH_YYDDDHH.REP_YYDDDHHMM	(report file)
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC (COBETRIEVE return status definitions)
C/	  CCT_Get_Config.Txt
C/	  FUT_Params.Txt
C/	  FDQ_Options.Txt
C/
C/	TEXT FILES USED:
C/	  FDQ_OPTIONS.TXT
C/
C/	SUBROUTINES CALLED:
C/	  CT_INIT (COBETRIEVE routine)
C/	  CT_WRITE_ARCV (COBETRIEVE routine)
C/	  CT_BINARY_TO_GMT
C/	  CT_GMT_TO_BINARY
C/	  FDQ_GET_OPTIONS
C/	  FDQ_LOAD_CONV
C/	  FDQ_GET_NXT_SEGMENT
C/	  FDQ_GET_NXT_SET
C/	  FDQ_GET_BRK_HSK
C/	  FDQ_PROC_CUR_SET
C/	  FDQ_OPEN_ARCV
C/	  FDQ_CLOSE_ARCV
C/	  FDQ_GEN_STAT
C/	  FDQ_FORM_STAT
C/	  FUT_QUAL_SUMMARY
C/	  LIB$GET_LUN
C/	  CCT_OPEN_CONFIG
C/	  CCT_GET_CONFIG_TOD
C/	  CCT_GET_CONFIG_IDX_TOD
C/	  CCT_CLOSE_CONFIG
C/
C/	ERROR HANDLING:
C/	  Error messages to SYS$OUTPUT (terminal if run interactively, log
C/	  file if run batch)
C/
C/	METHOD USED:
C/
C/	New PDL for FDQ (April 1987, E.F.)
C/
C/	  Call FDQ_LOAD_CONV to load buffer containing the
C/		conversions coefficients and the GRT lookup-tables;
C/	  Call FDQ_GET_OPTIONS to get process options, file or time range,
C/		and seconds for IFG sets.
C/	  Open the FDQ Report File.
C/        ** Write account invoking facility to the report file
C/        ** Write full command line with defaults to the report file
C/        ** Write logical name translation to the report file 
C/           (New requirements)-- H. Wang, STX, 12/17/90
C/	  Call CCT_OPEN_CONFIG to obtain access to all limits, flags and
C/	        tolerance datasets with time tags covering the FDQ run.
C/	  Do for all channels
C/	    Set MORE_SEGMENTS to true;
C/	  Enddo;
C/
C/	  Do while (MORE_SEGMENTS)
C/	    Call FDQ_GET_NXT_SEGMENT (OPTIONS, TIME_RANGE, FILE_SEG,
C/	         REF_CHAN, MORE_SEGMENTS, FILENAMES )
C/	    On the first pass if the orbital period was not entered on the
C/		 command line, call the FUT_ORBITAL_PERIOD to get the orbital
C/		 period in minutes.
C/	    If (MORE_SEGMENTS true for at least one channel) then
C/	      Call OPEN_ARCHIVES (FILENAMES, CT_LUN, TIME_RANGE, ATTITUDE_TYPE);
C/	      Do for CHAN = 1 to 4;
C/	        If (MORE_SEGMENTS for this channel) then
C/	          Set READ_FLAG for this channel to true;
C/	        Else
C/	          Set DONE_FLAG for this channel to true;
C/	        Endif;
C/	      Enddo;
C/	      Set NACTIVE to # channels with MORE_SEGMENTS true;
C/	      Do while (NACTIVE > 0)
C/
C/	        Zero out the ENG record buffer;
C/	        Call FDQ_GET_NXT_SET to do the following:
C/		* Read SCI archives if necessary to get next record (READ_FLAG
C/			will be cleared if SCI archive is read)
C/		* Determine new value of NACTIVE
C/		* Determine which channels should be marked "being processed"
C/		* Calculate AVG_TIME from the "being processed" channels' times
C/
C/	        If (NACTIVE > 0) then
C/	          Set ENG time (BNTM, GMT & PB5) to "average channel time" (AVG_TIME);
C/	          Call FDQ_GET_BRK_HSK to fill HSK buffers;
C/	          If (No HSK found) then
C/
C/	            Set "bad" quality flag on ENG record;
C/	            Increment NOHSK (for statistical purposes);
C/
C/	          Else
C/
C/	            Call CCT_Get_Config_TOD to get all of the archived
C/			configuration files with time tags covering the time
C/		        of the current science set. The files obtained will be
C/			science limits,engineering limits, limit check enable
C/			flags, index tolerances and index check enable flags.
C/	            Call CCT_Get_Config_Idx_TOD to get all of the archived
C/			configuration indexed files with time tags covering the
C/		        time of the current science set. The files obtained
C/		        will be conversion coefficients, GRT Tables and
C/			calibrator resistor information.
C/	            Call FDQ_PROC_CUR_SET to do the following:
C/	              * Convert counts in HSK buffers to eng. units in ENG buffers
C/	              * Interpolate the 2 converted ENG buffers
C/	              * Calculate internal reference source temperature
C/	              * Write out new index record if necessary
C/	              * Set quality flags on ENG record buffer
C/	              * Modify necessary fields in SCI record buffers and rewrite
C/			to archive for all "being processed" channels
C/		      * Write command files to trigger Engplots
C/	              * Set READ_FLAG for "being processed" channels
C/
C/	            Call FDQ_GEN_STAT to generate engineering statistics;
C/
C/	          Endif;		(No HSK found)
C/	          Write ENG record to archive;
C/	        Endif;		(NACTIVE > 0)
C/
C/	      Enddo;		(Do while (NACTIVE > 0))
C/
C/	    Endif;		(If (MORE_SEGMENTS true for at least 1 channel))
C/	    If (options = "one segment only") then
C/	      Set MORE_SEGMENTS of each channel to false;
C/	    Endif;
C/
C/	  Enddo;		(Do while (MORE_SEGMENTS))
C/
C/	  Call FUT_QUAL_SUMMARY to write the quality flag summary for the last
C/	       segment;
C/	  Call FDQ_FORM_STAT to form up engineering statistics;
C/	  call FUT_WRITE_ENGPLOT to write the command line text file to
C/	       trigger Engplots;
C/	  Close all archives;
C/	  Call CCT_CLOSE_CONFIG to close the configuration archives.
C/	  End.
C/
	IMPLICIT	NONE

	EXTERNAL	FUT_ERROR	! Condition handler

	EXTERNAL	FDQ_NORMAL
	EXTERNAL	FUT_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FUT_ABERR
	EXTERNAL	FDQ_CTINITERR
	EXTERNAL	FDQ_CTWRITERR
	EXTERNAL	FDQ_GETNXSETER
	EXTERNAL	FDQ_GETBRKHKER
	EXTERNAL	FDQ_OPENARCVER
	EXTERNAL	FDQ_PROCURSETR
	EXTERNAL	FDQ_CLOSARCVER
	EXTERNAL	FDQ_WRTPLOTERR
	EXTERNAL	FDQ_GETNXSEGER
	EXTERNAL	FDQ_GETOPTER
	EXTERNAL	FDQ_LOADCONVER
	EXTERNAL	FDQ_LOADLIMER
	EXTERNAL	FDQ_GENSTATER
	EXTERNAL	FDQ_FORMSTATER
	EXTERNAL	FDQ_MISSHKP
	EXTERNAL	FDQ_EOFHKP
	EXTERNAL	FDQ_HKPRFAIL
	EXTERNAL        FDQ_QUALSUMER
	EXTERNAL	FDQ_ERGETTIM
	EXTERNAL	FDQ_LUNGETERR
	EXTERNAL	FDQ_OPENERR
	EXTERNAL        FDQ_GETCONFIGERR
	EXTERNAL        FDQ_OPNCONFIGERR
	EXTERNAL        FDQ_CLSCONFIGERR
	EXTERNAL        FDQ_ERSUBX
	EXTERNAL	FDQ_GETORBPD
	EXTERNAL	FDQ_BADENG

	include		'ct$library:ctuser.inc'

	include		'($ssdef)'

	include		'(fdq_options)'

	include		'(fut_error)'

	include		'(fut_params)'

	include	        '(cct_get_config)'

	include	        '($jpidef)'

!!!!!!!!!!!!!!!!!!!!!!!
!		      !
!     "Constants"     !
!		      !
!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*2	BRK_MJ_MISS, ENG_CONV_INDX,
	1		EOF_HSK, GAIN_TOL, GLITCH_TOL,
	1		TLM_TOL

	PARAMETER	(BRK_MJ_MISS = 7)
	PARAMETER	(ENG_CONV_INDX = 88)
	PARAMETER	(EOF_HSK = 6)
	PARAMETER	(GAIN_TOL = 2)
	PARAMETER	(GLITCH_TOL = 1)
	PARAMETER	(TLM_TOL = 3)

	INTEGER*2	BAD, ETR, ENG, GOOD, IDX,CAL,
	1		HSK, RS1, RS2, RS3, RS4, RPT,
	2		SM1, SM2, SM3, SM4,
	3		SC1, SC2, SC3, SC4,
	4		ORS1, ORS2, ORS3, ORS4

	PARAMETER	(BAD = 2)
	PARAMETER	(ETR = 8)
	PARAMETER	(ENG = 7)
	PARAMETER	(GOOD = 1)
	PARAMETER	(CAL = 9)	
	PARAMETER	(HSK = 5)
	PARAMETER	(IDX = 6)
	PARAMETER	(RS1 = 1)
	PARAMETER	(RS2 = 2)
	PARAMETER	(RS3 = 3)
	PARAMETER	(RS4 = 4)
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

!!!!!!!!!!!!!!!!!!!!!
!		    !
!     Variables     !
!		    !
!!!!!!!!!!!!!!!!!!!!!

	dictionary	'fdq_idx'
	COMMON /LAST_IDX/	LAST_IDX_REC	! Common block for the last
	record /fdq_idx/	last_idx_rec	! IDX record; added 6-8-87

	dictionary 'fdq_eng'
	record /fdq_eng/ eng_rec

	dictionary 'nfs_sdf'
	record /nfs_sdf/ sci_rec(4)


	dictionary 'nfs_hkp'
	record /nfs_hkp/ hsk_recs(2)

	dictionary	'fex_scilim'

	dictionary	'fex_englim'

	dictionary	'fex_limflags'

	dictionary	'fex_idx_tols'

	dictionary	'fex_idx_flag'

	dictionary	'fex_grtrawwt'

	dictionary	'fex_grttrans'

	dictionary	'fex_samprate'

!	New structure for CCT_Get_Config_TOD

	structure / fdq_config /

	  record /fex_scilim/      scilim_rec

	  record /fex_englim/      englim_rec

	  record /fex_limflags/    limflags

	  record /fex_idx_tols/    idxtols

	  record /fex_idx_flag/    idxflags

	  record /fex_grtrawwt/      grtrawwt_rec

	  record /fex_grttrans/      grttrans_rec

	  record /fex_samprate/	     samprate

	endstructure	

	record / fdq_config / config

	record /config_status/ stat(8)
	record /config_status/ stat2(3)

!	Structured datasets

	integer*4	ndset /8/		! Number of datasets
	character*32    dset(8)			! Names of datasets
	data dset(1) / 'CSDR$FIRAS_REF:fex_scilim' /
	data dset(2) / 'CSDR$FIRAS_REF:fex_englim' /
	data dset(3) / 'CSDR$FIRAS_REF:fex_limflags' /
	data dset(4) / 'CSDR$FIRAS_REF:fex_idx_tols' /
	data dset(5) / 'CSDR$FIRAS_REF:fex_idx_flag' /
	data dset(6) / 'CSDR$FIRAS_REF:fex_grtrawwt' /
	data dset(7) / 'CSDR$FIRAS_REF:fex_grttrans' /
	data dset(8) / 'CSDR$FIRAS_REF:fex_samprate' /
	Character*32   Calrs_set/ 'CSDR$FIRAS_REF:fex_av_calrs' /
	integer*4	size(8) / 1024, 4096, 512, 1024, 512, 256, 384, 64/  
!				 Dataset sizes
	integer*4       con_lun(8)		! Logical unit numbers
	integer*4       start_time(2)		! Data start time for FDQ run
	integer*4       stop_time(2)		! Data stop time for FDQ run
	integer*4       ncache / 1 /		! Number of caches- not used
	integer*4       index(8)		! Initial cache pointers
	integer*4	ref_count		! Reference counter for cache
	logical*1       new_segment(8)		! Flag for new segments accessed
	character*1     access_mode / ' ' /     ! Access mode sequential-default
	character*19    dummy_space		! Dummy space pad
	character*14    con_time		! Data time for report
        character*14    conerr_time             ! Error time for CCT-Get_Config
	integer*4       sec_64(2)               ! Offset from start telemetry
	data            sec_64(1) / 640000000 / ! time for collect time
	integer*4       len / 2 /		! Lib$Subx length
	integer*4       tstart_time(2)          ! Dummy start time for data
					        ! start for CCT_Open_Config
	real*4		samprate		! Sampling rate of FEX_SAMPRATE

C	Indexed datasets.

	integer*4       size2(3) / 128, 7012, 20/
	character*5     access_mode2 / 'KEYED' /     ! Access mode sequential-default
	character*32    diset(3)			! Names of datasets
	data diset(1) / 'CSDR$FIRAS_REF:fdb_limcuvky' /
	data diset(2) / 'CSDR$FIRAS_REF:fdb_grtrich' /
	data diset(3) / 'CSDR$FIRAS_REF:fdb_grtcal' /
	integer*4	nidset / 3 /		! Number of datasets
	integer*4       config_lun(3)		! Logical unit numbers
	integer*4       index2(3)		! Initial cache pointers
	logical*1       new_segment2(3)		! Flag for new segments accessed

C	Functions

	INTEGER*4	FDQ_GET_OPTIONS
	INTEGER*4	FDQ_LOAD_CONV
	INTEGER*4	FDQ_GET_NXT_SEGMENT
	INTEGER*4	FDQ_OPEN_ARCV
	INTEGER*4	FDQ_GET_NXT_SET
	INTEGER*4	FDQ_GET_BRK_HSK
	INTEGER*4	FDQ_PROC_CUR_SET
	INTEGER*4	FDQ_GEN_STAT
	INTEGER*4	FDQ_FORM_STAT
	INTEGER*4	FDQ_CLOSE_ARCV
	INTEGER*4       FUT_QUAL_SUMMARY
	INTEGER*4       FUT_WRITE_ENGPLOT
	INTEGER*4	FUT_ORBITAL_PERIOD
	INTEGER*4	CCT_OPEN_CONFIG
	INTEGER*4       CCT_GET_CONFIG_TOD
	INTEGER*4       CCT_GET_CONFIG_IDX_TOD
	INTEGER*4	CCT_CLOSE_CONFIG
	INTEGER*4       LIB$GET_LUN     ! Get logical unit number
	INTEGER*4       LIB$SUBX        ! Extended precision subtract
        Integer*4       LIB$GETJPI
	integer*4	STATUS		! Status of FDQ processing
	integer*4	RETSTAT		! Return status from FDQ function
	integer*4	IOSTATUS        ! return status from I/O
	integer*4	SUCCESS / 1 /, ERROR_ST / 2 /  ! Values for status

	CHARACTER	TIME_RANGE*30
        Character*14    CONGMT_START  
	CHARACTER*14    GMT_START       ! Start GMT of data
	CHARACTER*14    GMT_STOP        ! Stop GMT of data
	CHARACTER*39    FILE_SEG	! User specified file-segment
	CHARACTER*39    FILENAMES(5)    ! Complete filenames for CT open
	CHARACTER*39    OUTNAMES(7)    ! 
	CHARACTER*39    BLANK           ! Blank characters
	CHARACTER*1     BLANKS(39) / 39 * ' ' /
	EQUIVALENCE     ( BLANK, BLANKS(1) )
	LOGICAL*1	SCI_OPEN(4)/.FALSE.,.FALSE.,.FALSE.,.FALSE./

	LOGICAL*1	DONE(4), IPDU_RELAY(8),
	1		MORE_SEGMENTS(4), NEW_SEGMENTS, first_time,
	1		READ_NXT(4), STMON_SYNC(2), XCAL_SET/.false./,SEC_DEF,
	1		WRITE_ASC_CONV, LIST_FULL, prt,report_def

	LOGICAL*1	NO_RECS / .FALSE. /	! Flag True if no science found.
	LOGICAL*1	FIRST / .TRUE. /
	LOGICAL*1	FIRST_ENG
	LOGICAL*1       TIME_LE,time_lt
	LOGICAL*1	ENG_badtime/.false./
        INTEGER*4       ENG_PREVTIME(2)
	INTEGER*2	ACTIVE(4), CT_LUN(22), CT_STAT(20), FIRST_MJ,
	1		I, ISTAT, NACTIVE, neng_writ, ix, chan,k,
	1		NNOHSK, NPROC, nsci_read(4), OPTIONS(10),
	1		PROC_CHAN(4), tols(10), SC_SIDE/0/, LIMIT
        Character*40    option_name(3)/'* Process all segments (default) ',
	1	                      '* No segment is excluded by default',
	1	                      '* Process daily segment '/
	INTEGER*4	LUN
	INTEGER*4       PLOT_START(2), PLOT_STOP(2)  ! Start and stop plot

	INTEGER*4	AVG_TIME(2),PreAVG_TIME(2),DELIVERY_TIME(2),ifg_count(4),
	1		ref_chan, set_sec, englim,
	1		ZEROES /0/

	INTEGER*4       ATTITUDE_TYPE	! Type of attitude to be used.
	LOGICAL*1       PLOT_TABLE(118) ! Table of flags to trigger Engplots

	INTEGER*4	PLOT_DEVICE	! Plot device: line or laser printer
	INTEGER*4	MIN_OFFTIME, MAX_OFFTIME  ! Offset in hours for plots
	REAL*4		ORBITAL_PD	! Time for one orbit in minutes
	LOGICAL*1	FORM_STAT	! Form trends from valid eng records;
					! skip those of FAC_NO_HSKP or
					! FAC_TLM_ERROR science data quality
	CHARACTER*33	REPORT_FILE	! Report filename
	CHARACTER*14	CURGMT		! Current system run time in GMT
	LOGICAL*1	REPORT/.FALSE./	! Flag to write report or not
	Character*24	Lim_Type(2)
	Data Lim_Type(1) / 'Red Limits only ' /
	Data Lim_Type(2) / 'Red and Yellow Limits' /
	Character*20    Dev
	Character*16    Att


        Integer*4       Cut_Register_version
        Integer*4       Cut_Display_Banner
        Integer*4       num_vol/80/
        Integer*4       rstatus
        Integer*4       lun_out/6/
        Character*6     version
        Character*8     owner
        Parameter       (version='9.8')
        Logical*1       first_ref/.true./  
        Logical*1       distim/.false./  
        integer*4       ETR_WRIT,IDX_WRIT
        Integer*4       CUT_TRANSLATE_ARCHIVE_ID
        Integer*4       Flen,Tlen
        Character*72    Log,Tlog,Flog
        Integer*4       UPM_GET_VALUE
        CHARACTER*256   COMMAND_LINE
        Integer*2       Clen
        Integer*2       EARL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set status for FDQ processing to success.

	STATUS = SUCCESS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			                      !
!     Declare exit handler                    !
!     Replace for Build 3.1 version of FDQ    !
!     with establish condition handler.       !
!			                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	Rstatus = Cut_Register_Version(version)
	Rstatus = CUT_Display_Banner(lun_out,num_vol,
	1		'FIRAS Facility FDQ_Data_Qualify')
	Write (lun_out,61)
61	Format (//)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						    !
!     Define tolerances for analog index values     !
!     to signify "change of state"		    !
!						    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	tols(glitch_tol) = 250
	tols(gain_tol) = 20
	tols(tlm_tol) = 10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				 !
!     Get processing options     !
!				 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	RETSTAT = FDQ_GET_OPTIONS ( OPTIONS, TIME_RANGE, FILE_SEG,
	1  REF_CHAN,  SET_SEC, SEC_DEF,
	2  ATTITUDE_TYPE, LIMIT, PLOT_DEVICE, MIN_OFFTIME,
	3  MAX_OFFTIME,  REPORT_FILE, CURGMT, REPORT,report_def)
	if ( retstat .ne. %loc(FDQ_NORMAL) ) then
	  if ( retstat .eq. %loc(FDQ_ABERR) ) retstat = zeroes
	  call lib$signal(FDQ_GETOPTER,%val(1),%val(retstat))
	  STATUS = ERROR_ST
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                !
!     Open report file.          !
!                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Get a unit number and open the report file.

	If ( (status .eq. success) .AND. REPORT) Then
		retstat = Lib$Get_Lun (lun)
		ct_lun(rpt) = lun
		If ( retstat .ne. SS$_Normal ) Then
		      status  = error_st
		      Call Lib$Signal(FDQ_LUNGETERR, %val(1), %val(retstat))
		Endif
	ENDIF

!	Check status before processing.

	If ( (status .eq. success) .AND. REPORT) Then

	      Open( Unit=Ct_Lun(rpt), File=Report_File,
	1	    Status='NEW', Form='FORMATTED', Access='SEQUENTIAL',
	2	    Organization='SEQUENTIAL', Iostat=iostatus )

	      If ( iostatus .ne. zeroes ) Then
		    status = error_st
	            Call Lib$Signal(FDQ_OPENERR,%val(1),%val(iostatus))
	      Endif
c
c ** call routines to put program version number on the report file
c
              rstatus = CUT_Register_Version(version)
              rstatus = CUT_Display_Banner(Ct_Lun(rpt),num_vol,
	1		'FIRAS Facility FDQ_Data_Qualify')
              Rstatus = Lib$getJPI(JPI$_UserName,,,,owner,)
              WRITE(CT_lun(rpt),62)owner, curgmt
 62     Format(//' Run By :    ',a,'                     TIME:    ',a)
             
	Endif	! status is success

	IF ( (STATUS .EQ. SUCCESS) .AND. REPORT ) THEN
          Write(ct_lun(rpt),64)
          log(1:14)='CSDR$FIRAS_RAW'
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),65) tlog(1:tlen)
          log(1:14)='CSDR$FIRAS_IN '
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),81) tlog(1:tlen)
          log(1:14)='CSDR$FIRAS_OUT'
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),82) tlog(1:tlen)
          log(1:14)='CSDR$FIRAS_REF'
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),66)tlog(1:tlen)
          log(1:18)='CSDR$ANCIL_ARCHIVE'
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),67)TLOG(1:tlen)
          log(1:18)='CSDR$SPACE_ARCHIVE'
          rstatus=cut_translate_archive_id(log,flog,flen,tlog,tlen)
          Write(ct_lun(rpt),68)Tlog(1:tlen)
 63     Format(/1x,'Options:')
 64     Format(//1x,'Logical names used:')
 65     Format(5x,'CSDR$FIRAS_RAW = ', a)
 81     Format(5x,'CSDR$FIRAS_IN = ', a)
 82     Format(5x,'CSDR$FIRAS_OUT = ', a)
 66     Format(5x,'CSDR$FIRAS_REF = ', a)
 67     Format(5x,'CSDR$ANCIL_ARCHIVE = ', a)
 68     Format(5x,'CSDR$SPACE_ARCHIVE = ', a)
 69     Format(//1x,'Command Line with defaults:'/5x,a/5x,a)
 70     format(//1x,'Report File Used: ', a)
          write(ct_lun(rpt),70)report_file
          Rstatus = UPM_GET_VALUE('$line',command_line,clen)
          If (sec_def) command_line(clen+1:clen+11)='/SETSEC=16'
          If (report_def) then
            clen = clen +12
            command_line(clen:clen+7)='/REPORT'
          Endif 
          Write( Ct_lun(rpt),69)command_line(1:50),command_line(51:clen+7)
          Write(CT_LUN(RPT),63)
          do i =1, 3
            If (options(i) .eq. 1) then
              WRITE (ct_lun(rpt), 6001) option_name(i)
            Endif 
          enddo
	    IF (options(opt_time_range) .ne. optv_no_time_range .and.
	1	options(opt_time_range) .ne. optv_spec_seg ) THEN
		WRITE (ct_lun(rpt), 6003) TIME_RANGE
	    ELSEIF (options(opt_time_range) .eq. optv_spec_seg) THEN
		write (ct_lun(rpt), 6002) file_seg
	    ENDIF
	    IF (attitude_type .eq. fac_none) THEN
		write (ct_lun(rpt), 7100)
	    ELSE
		IF (attitude_type .eq. fac_simulated) THEN
		    att = 'Simulated   '
		ELSEIF (attitude_type .eq. fac_fine_with_Dirbe ) THEN
		    att = 'Fine with Dirbe'
		ELSEIF (attitude_type .eq. fac_definitive ) THEN
		    att = 'Definitive '
		ELSEIF (attitude_type .eq. fac_coarse ) THEN
		    att = 'Coarse  '
		ENDIF
		write (ct_lun(rpt), 7200) att
	    ENDIF

	    IF (plot_device .eq. fac_none) THEN
		write (ct_lun(rpt), 7300)
	    ELSE
		IF (plot_device .eq. fac_lineprinter) THEN
		    dev = 'lineprinter '
		ELSE
                  If (Plot_device .eq. Fac_lasertf) then
		    dev = 'laserprinter("TF")'
                  Else
		    dev = 'laserprinter("QMS")'
                  Endif
		ENDIF
		write(ct_lun(rpt),7400)dev,lim_type(limit),min_offtime,max_offtime
	    ENDIF
          Fut_report_lun = ct_lun(rpt)
	  Call Lib$Establish ( FUT_ERROR )
	ENDIF		! Status is success

6001	FORMAT (5x,a)
6002    FORMAT (/1x,'One Segment Specified: ',a/)
6003	FORMAT (/1x,'Time Range: ' , A )
7000    format (/1x,'Seconds for science grouping =', I6,
	1	' Orbital Period =', F10.4/)
7100    format (/1x,'No attitude solution will be used for the run.')
7200    format (/1x,'Attitude: ',A)
7300    format (/1x,'No Engplots will be generated for the run.')
7400    format (/1x,'Engplots: ', A, ' for ', A/
	1	11x,'Offset Hours for Plot Start =',I6,
	2	'  for Plot Stop =', I6 )
!
!	Initialize variables
!
	IF ( STATUS .EQ. SUCCESS ) THEN
	  DO I = 1, 4
	    MORE_SEGMENTS(I) = .TRUE.
	  ENDDO

	  CALL CT_INIT (CT_STAT)
	  IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	      STATUS = CT_STAT(1)
	      CALL LIB$SIGNAL(FDQ_CTINITERR,%VAL(1),%VAL(STATUS))
	      STATUS = ERROR_ST
	      ISTAT = BAD				! BAD = 2
!	      STOP 'DATA QUALIFY aborts'
	   ENDIF
	ENDIF
        etr_writ=0

	DO WHILE (MORE_SEGMENTS(1) .OR. MORE_SEGMENTS(2) .OR.
	1	MORE_SEGMENTS(3) .OR. MORE_SEGMENTS(4))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				   !
!     Call FDQ_GET_NXT_SEGMENT     !
!				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  if ( status .eq. success ) then

	    RETSTAT = FDQ_GET_NXT_SEGMENT (OPTIONS, TIME_RANGE,
	1	FILE_SEG, REF_CHAN, MORE_SEGMENTS, FILENAMES)
	    if ( retstat .ne. %loc(FDQ_NORMAL) ) then
	      if ( retstat .eq. %loc(FDQ_ABERR) ) retstat = zeroes
	      if ( .NOT.MORE_SEGMENTS(1) .AND. .NOT.MORE_SEGMENTS(2) .AND.
	1	.NOT.MORE_SEGMENTS(3) .AND. .NOT.MORE_SEGMENTS(4) ) then
		NO_RECS = .true.
		if (REPORT) Write (ct_lun(rpt),500)
	      else
		call lib$signal(FDQ_GETNXSEGER,%val(1),%val(retstat))
	      endif
	      STATUS = ERROR_ST
	    endif
	
	    if ( status .eq. success ) then

	      if ( first ) then	 ! Get start and stop for attitude & config

	        gmt_start(1:14) =  time_range(1:14)
	        gmt_stop(1:14) = time_range(16:29)
	        first = .false.
		retstat = FUT_Orbital_Period ( gmt_start, orbital_pd )
		if (retstat .ne. %loc(FUT_Normal)) then
		  status = error_st
		  call Lib$Signal(FDQ_GETORBPD, %val(1),%val(retstat))
	        endif
		IF (report) write (ct_lun(rpt), 7000) set_sec, orbital_pd
	        if ((options(opt_time_range).eq.optv_spec_seg).AND.report)then
	          WRITE (ct_lun(rpt), 6003) TIME_RANGE, TIME_RANGE(1:14),
	1		TIME_RANGE(16:29)
		endif

!       Open configuration files.

!	Compute binary start and stop times for config file ranges.
                CONGMT_START = time_range(1:14)
	        call ct_gmt_to_binary(congmt_start,start_time)
	        call ct_gmt_to_binary(gmt_stop, stop_time)

	        if ( status .eq. success ) then
	          RETSTAT = CCT_OPEN_CONFIG ( start_time, stop_time,
	1	    ndset, dset, size, access_mode, ncache, con_lun,
	2	    index, stat, ref_count )
	          if ( .not. retstat ) then
	            STATUS = ERROR_ST
	            call lib$signal(FDQ_OPNCONFIGERR,%val(1),%val(retstat))
	          endif
	        endif

	        if ( status .eq. success ) then
	          RETSTAT = CCT_OPEN_CONFIG ( start_time, stop_time,
	1	    nidset, diset, size2, access_mode2, ncache, config_lun,
	2	    index2, stat2, ref_count )
	          if ( .not. retstat ) then
	            STATUS = ERROR_ST
	            call lib$signal(FDQ_OPNCONFIGERR,%val(1),%val(retstat))
	          endif
	        endif

	      endif	! First time in FDQ run

	      IF (MORE_SEGMENTS(1) .OR. MORE_SEGMENTS(2) .OR.
	1	MORE_SEGMENTS(3) .OR. MORE_SEGMENTS(4)) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					 !
!     Open CT archives for this pass     !
!					 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                FIRST_TIME = .TRUE.     ! Flag passed to FDQ_GET_BRK_HSK
                                        ! indicating that this is the first
                                        ! pass on the current segment
                FIRST_ENG = .TRUE.
	        NEW_SEGMENTS = .TRUE.	! This flag is to be passed
					! to FDQ_PROC_CUR_SET, then to
					! FDQ_MAKE_IDX and finally to
					! FDQ_NEW_IDX to signal it
					! to start over with a brand new
					! index record
					! (FDQ_NEW_IDX will set it
					! to false immediately)
		RETSTAT = FDQ_OPEN_ARCV ( FILENAMES, SCI_OPEN, TIME_RANGE,
	1				  CT_LUN,  ENGLIM,OUTNAMES )
		IF ( RETSTAT .NE. %LOC(FDQ_NORMAL) ) THEN
	          IF ( RETSTAT .EQ. %LOC(FDQ_ABERR) ) RETSTAT = ZEROES
		  CALL LIB$SIGNAL(FDQ_OPENARCVER,%VAL(1),%VAL(RETSTAT))
		  STATUS = ERROR_ST
	        ENDIF

	        IF ( STATUS .EQ. SUCCESS ) THEN
	          NACTIVE = 0
	          NENG_WRIT = 0
	          NNOHSK = 0
	          DO I = 1, 4
	            NSCI_READ(I) = 0
		    IFG_COUNT(I) = 0
		    ENGLIM = FAC_NOT_PRESENT
	            IF (FILENAMES(I) .NE. BLANK) THEN
	              DONE(I) = .FALSE.
	              READ_NXT(I) = .TRUE.
	              NACTIVE = NACTIVE + 1
	            ELSE
	              DONE(I) = .TRUE.
	              READ_NXT(I) = .FALSE.
	            ENDIF
	          ENDDO
	        ENDIF		! Status is success

	        DO WHILE ((NACTIVE .GT. 0) .AND. (STATUS .EQ. SUCCESS))
                    distim = .false.
	            RETSTAT = FDQ_GET_NXT_SET (READ_NXT, DONE,
	1		CT_LUN, SCI_REC, NACTIVE, ACTIVE, NPROC, PROC_CHAN,
	2		AVG_TIME, EARL,NSCI_READ, SET_SEC, IFG_COUNT)
		    if ( retstat .ne. %loc(FDQ_NORMAL) ) then
	              if ( retstat .eq. %loc(FDQ_ABERR) ) retstat = zeroes
		      call lib$signal(FDQ_GETNXSETER,%val(1),%val(retstat))
	              STATUS = ERROR_ST
		    endif
                    If(time_lt(avg_time, preavg_time)) distim=.true.   	
	            IF ((NACTIVE .GT. 0) .AND. (status.eq.success)) THEN

	              CALL LIB$MOVC5 (4, ZEROES, 0, 1024, ENG_REC)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									       !
!     Fill in GMT and binary times on the ENG record header using AVG_TIME     !
!									       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	              CALL LIB$MOVC3 (8, sci_rec(earl).collect_time.midpoint_
	1	            time, ENG_REC.CT_HEAD.TIME)
	              call CT_Binary_To_Gmt (sci_rec(earl).collect_time.midpoint
	1	                             _time,
	1			eng_rec.ct_head.gmt )
                  
	              If  (.not. First_eng) Then
                        If(Time_LE(Eng_rec.ct_head.time, eng_prevtime)) then
                          Eng_badtime = .true.
                        else
                          eng_badtime = .false.
                        Endif
                      Endif
                      Eng_prevtime(1) = Eng_rec.ct_head.time(1)
                      Eng_prevtime(2) = Eng_rec.ct_head.time(2)
                      If (.not. distim) then
	              RETSTAT = FDQ_GET_BRK_HSK (CT_LUN(HSK),
	1	        AVG_TIME,
	1		HSK_RECS, FIRST_MJ, STMON_SYNC, IPDU_RELAY,
	1		ISTAT, CT_STAT,first_time)
		      if (retstat .eq. %loc(FDQ_NORMAL)) THEN
	                IF (ISTAT .NE. GOOD) THEN
	                  IF (ISTAT .EQ. BRK_MJ_MISS) THEN	! Missing 1 or 2 MJ FR
			    call lib$signal(FDQ_MISSHKP)
	                  elseif (istat .eq. eof_hsk) then
			    call lib$signal(FDQ_EOFHKP)
	                  else
			    call lib$signal(FDQ_HKPRFAIL)
	                  endif
	                 eng_rec.en_tail.hskp_flag = 1  ! Signify no HSK in ENG buffer
	                ENDIF		! IF (ISTAT .NE. GOOD)
		      else
	                if ( retstat .eq. %loc(FDQ_ABERR) ) retstat = zeroes
		        call lib$signal(FDQ_GETBRKHKER,%val(1),%val(retstat))
	                STATUS = ERROR_ST
		      endif

                     
	              if (( eng_rec.en_tail.hskp_flag .eq. 1 )
	1		.and. ( status .eq. success )) then
	                nnohsk = nnohsk + 1
	              endif		! if (eng_rec.en_tail.hskp_flag .eq. 1)
                    endif ! if distim true
!	Call COBETRIEVE function to get all limits and enable flags datasets
!	with time tag ranges covering the Science set time.

		      if ( status .eq. success ) then
	                RETSTAT = CCT_GET_CONFIG_TOD( AVG_TIME, NDSET,
	1		  SIZE, CON_LUN, INDEX, CONFIG, NEW_SEGMENT, STAT)
	                if ( .not. retstat ) then
	                  status = error_st
	                  call lib$signal( FDQ_GETCONFIGERR,
	1		    %val(1), %val(retstat))
			  Call CT_Binary_To_Gmt(avg_time,conerr_time)
			  IF (REPORT) write(ct_lun(rpt),6009) conerr_time,
	1			gmt_start, gmt_stop
 6009	   format(1x/1x,'Error CCT_Get_Config_TOD. Input Time= ',a/
	1	     1x,'Time for CCT_Open_Config= ',a,' to ',a)
			else
                         IF ( first_ref) then
			  Call CT_Binary_To_Gmt(avg_time,con_time)
                          Write(CT_LUN(RPT),6008)con_time
                          first_ref=.false.
                         Endif  
 			  do i = 1, ndset
			    if (new_segment(i).and.report) then
			      write(ct_lun(rpt),6010) dset(i)
			    endif
			  enddo
	                endif
		      endif

 6008	FORMAT(1x,'Reference datasets Accessed at time: ',a )
 6010	FORMAT(5x,a )
!
!	Call COBETRIEVE function to get all conversions datasets.

		      if ( status .eq. success ) then
	                RETSTAT = CCT_GET_CONFIG_IDX_TOD( AVG_TIME, NIDSET,
	1		  CONFIG_LUN, INDEX2, NEW_SEGMENT2, STAT2)
	                if ( .not. retstat ) then
	                  status = error_st
	                  call lib$signal( FDQ_GETCONFIGERR,
	1		    %val(1), %val(retstat))
			  Call CT_Binary_To_Gmt(avg_time,conerr_time)
			  IF(REPORT) write(ct_lun(rpt),6011) conerr_time,
	1			gmt_start, gmt_stop
 6011	   format(1x/1x,'Error CCT_Get_Config_IDX_TOD. Input Time= ',a/
	1	     1x,'Time for CCT_Open_Config= ',a,' to ',a)
			else
			  Call CT_Binary_To_Gmt(avg_time,con_time)
			  do i = 1, nidset
			    if (new_segment2(i).and.report) then
			      write(ct_lun(rpt),6010)diset(i)
			    endif
			  enddo
			  if (new_segment2(i-1).and.report) then
			      write(ct_lun(rpt),6010) calrs_set
                          endif
	                endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				      !
!     Load conversions and limits     !
!				      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	                if ( status .eq. success) then
	                  RETSTAT = FDQ_LOAD_CONV(config_lun, new_segment2 )
	                  if ( retstat .ne. %loc(FDQ_NORMAL) ) then
	                    if ( retstat .eq. %loc(FDQ_ABERR) ) retstat = zeroes
	                      call lib$signal(FDQ_LOADCONVER,%val(1),%val(retstat))
	                      STATUS = ERROR_ST
	                    endif
	                endif
		      endif

!	Call the FDQ routine to process the current set of science records.

	              if ( status .eq. success ) then
	                RETSTAT = FDQ_PROC_CUR_SET (new_segments,
	1		  ct_lun, hsk_recs, first_mj, avg_time,Eng_badtime,
	1		  ipdu_relay, stmon_sync, nproc, proc_chan,
	1		  read_nxt, tols, sci_rec, eng_rec, xcal_set,
	1		  config.scilim_rec, config.englim_rec,
	1		  config.limflags, config.idxtols,
	1		  config.idxflags, config.grtrawwt_rec,
	1	          Config.grttrans_rec, 
	1		  attitude_type, gmt_start, gmt_stop,
	1		  sc_side, limit, plot_device, min_offtime, max_offtime,
	1		  plot_table, plot_start, plot_stop,
	1		  ifg_count, englim, idx_writ,time_range,
	1		  config.samprate.sampling_rate)
			  if ( retstat .ne. %loc(FDQ_NORMAL) ) then
			    if(retstat.eq.%loc(FDQ_ABERR)) retstat = zeroes
		            call lib$signal(FDQ_PROCURSETR,
	1			%val(1),%val(retstat))
	                    STATUS = ERROR_ST
		          endif
	              endif
		      form_stat = .TRUE.
		      do i=1,4
			 if (eng_rec.en_head.dq_summary_flag(i) .gt.
	1			fac_many_red) then
			    form_stat = .FALSE.
			 endif
		      enddo
	              if ( status .eq. success .and. form_stat) then
	                retstat = FDQ_GEN_STAT ( ct_lun(ETR), eng_rec,
	1		  orbital_pd,ETR_WRIT )
		        IF (retstat .NE. %loc(FDQ_NORMAL)) then
			  if ( retstat .eq. %loc(FDQ_ABERR) )
	1			retstat = zeroes
		          call lib$signal(FDQ_GENSTATER,
	1		    %val(1),%val(retstat))
		          STATUS = ERROR_ST
	        endif
                      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							    !
!	Write out Engineering record through COBETRIEVE     !
!							    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	        if ( status .eq. success ) then
                First_eng = .false.
                If ( .not.Eng_badtime) then
		ENG_REC.CT_HEAD.DATASET_ID=CTU_$FIR_EDF
	        CALL CT_WRITE_ARCV (, CT_LUN(ENG), ENG_REC, CT_STAT)
	        IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
		  status = ct_stat(1)
		  call lib$signal(FDQ_CTWRITERR,%val(2),%val(status),
	1		%val(ENG))
		  STATUS = ERROR_ST
!	          STOP 'DATA QUALIFY aborts'
	        ENDIF
                else
		  call lib$signal(FDQ_BadENG, %val(1),ENG_REC.ct_head.GMT)
                endif                
	        endif
	        if ( status .eq. success .and. (.not. eng_badtime)) then
	        NENG_WRIT = NENG_WRIT + 1
	        endif
	        ENDIF		! IF (NACTIVE .GT. 0 and Status is success)
              preavg_time(1) = avg_time(1)
              preavg_time(2) = avg_time(2)
	    ENDDO	! DO WHILE (NACTIVE .GT. 0 .AND. Status is success)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								       !
!     Write out last index record (if none has been written in the     !
!     last call to FDQ_MAKE_IDX)				       !
!								       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	   if ( status .eq. success .and. (.not. Eng_badtime)) then
	    if (last_idx_rec.header.gmt_end_time_ascii .ne.
	1		eng_rec.ct_head.gmt) then
	      last_idx_rec.header.gmt_end_time_ascii =
	1		eng_rec.ct_head.gmt
	      call lib$movc3 (8, eng_rec.ct_head.time,
	1		last_idx_rec.header.binary_end_time)
	      LAST_IDX_REC.HEADER.DATASET_ID=CTU_$FIR_IDX
	      call CT_WRITE_ARCV (, ct_lun(idx), LAST_IDX_REC, CT_STAT)
	      if (ct_stat(1) .ne. ctp_normal) then
	        status = ct_stat(1)
		call lib$signal( FDQ_CTWRITERR, %val(2),
	1		%val(status),%val(IDX))
	        STATUS = ERROR_ST
!	        stop 'DATA QUALIFY aborts'
              Else
                Idx_writ = Idx_writ + 1
	      endif
	    endif

C Eventually put summary into report routine.

	        If(report) then
                  write(ct_lun(rpt),6200)
                  Do k=1,4
                   write(ct_lun(rpt),6201)filenames(k),outnames(k)(16:39),
	1	                          nsci_read(k)
                  enddo
		  WRITE (ct_lun(rpt), 6202)FILENAMES(5),outnames(5)(16:39),
	1	                           NENG_WRIT
                  WRITE(ct_lun(rpt),6204)outnames(7)(16:39), IDX_WRIT    
                endif
	    endif		! Status is success
	   ENDIF		! IF (MORE_SEGMENTS(1) .OR. etc., etc.,..
	  endif			! Status is success after get_nxt_segment

	  endif		! Status is success at beginning of while more seg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								 !
!     If only processing one set of segments set switches so     !
!     that it will not go on after the first pass		 !
!								 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  if ((options(opt_time_range) .eq. optv_spec_seg) .or.
	1   (status .ne. success)) then
	    do i = 1, 4
	      more_segments(i) = .false.
	    enddo
	  endif

6200	  FORMAT (///1x,'Input Dataset', '             Output Dataset',
	1	 '         Records Written' )
C
6201      FORMAT(1x, a23,3x,a,2x,I5)
6202      FORMAT(1x, a23,3x,a,2x,I5)
6203      FORMAT(27x,a,2x,I5)
6204      FORMAT(27x,a,2x,I5)
6205      FORMAT(//1x,'# Dummy Eng. records (no HSK): ', I5//)

	ENDDO		! DO WHILE (MORE_SEGMENTS(1) .OR. etc., etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							  !
!     Call FDQ_FORM_STAT to form engin. statistics        !
!							  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! This forms and writes stats to the CT archive for the
	! last data set even though a full orbital period has not passed.
	if ( status .eq. success ) then
	    retstat = FDQ_FORM_STAT (ct_lun(ETR),ETR_WRIT)
	    IF (retstat .ne. %loc(FDQ_NORMAL)) then
	      if (retstat.eq.%loc(FDQ_ABERR)) retstat = zeroes
	      call lib$signal(FDQ_FORMSTATER,%val(1),%val(retstat))
	      STATUS = ERROR_ST
	    ELSE
              IF (REPORT) THEN
                WRITE(ct_lun(rpt),6203)outnames(6)(16:39), ETR_WRIT    
              ENDIF
            Endif
	endif		! Status is success
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!     Call FUT_QUAL_SUMMARY to print data quality flags   !
!	   for last segment.                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if ( (status .eq. success) .AND. report) then
          WRITE(ct_lun(rpt),6205) nnohsk    
	  prt = .true.
	  chan = 4
	  retstat = FUT_QUAL_SUMMARY ( chan,
	1   sci_rec(chan).dq_data.data_quality, prt, ct_lun(rpt))
	    IF (retstat .ne. %loc(FUT_NORMAL)) then
	      if (retstat.eq.%loc(FUT_ABERR)) retstat = zeroes
	      call lib$signal(FDQ_QUALSUMER,%val(1),%val(retstat))
	      STATUS = ERROR_ST
	    ENDIF
	endif		! Status is success



	RETSTAT = FUT_WRITE_ENGPLOT (Plot_Table, Plot_Start, Plot_Stop,
	1	  englim, Plot_Device)
	IF (retstat .ne. %loc(FUT_NORMAL)) THEN
	  if (retstat.eq.%loc(FUT_ABERR)) retstat = zeroes
	  call lib$signal(FDQ_WRTPLOTERR,%val(1),%val(retstat))
	  STATUS = ERROR_ST
	ENDIF

	RETSTAT = FDQ_CLOSE_ARCV (SCI_OPEN, CT_LUN,   ENGLIM)
	IF (retstat .ne. %loc(FDQ_NORMAL)) THEN
	  if (retstat.eq.%loc(FDQ_ABERR)) retstat = zeroes
	  call lib$signal(FDQ_CLOSARCVER,%val(1),%val(retstat))
	  STATUS = ERROR_ST
	ENDIF

	RETSTAT = CCT_CLOSE_CONFIG ( NDSET, CON_LUN, INDEX )
	if ( .not. retstat ) then
	    STATUS = ERROR_ST
	    call lib$signal(FDQ_CLSCONFIGERR, %val(1), %val(retstat))
	endif
	RETSTAT = CCT_CLOSE_CONFIG ( NIDSET, CONFIG_LUN, INDEX2 )
	if ( .not. retstat ) then
	    STATUS = ERROR_ST
	    call lib$signal(FDQ_CLSCONFIGERR, %val(1), %val(retstat))
	endif
	IF ( status .eq. success .OR. NO_RECS) THEN
	  call lib$signal(FDQ_NORMAL)
	  IF (REPORT) THEN
C	    write(unit=ct_lun(rpt),FMT=600)
	    Close(ct_lun(rpt))
	  ENDIF
	  call exit(ss$_normal)
	ELSE
	  call lib$signal(FDQ_ABERR)
          IF (REPORT) THEN
C	    write(unit=ct_lun(rpt),FMT=700)
	    Close(ct_lun(rpt))
	  ENDIF
	  call exit(ss$_abort)
	ENDIF
  500	format(/5x,' No catalog records found for any channel.')
C 600	format(/5x,' The FIRAS Data_Qualify Program Run Completed Successfully.')
C 700	format(/5x,' FIRAS_DATA_QUALIFY PROGRAM TERMINATED WITH ERROR. ')
	END
