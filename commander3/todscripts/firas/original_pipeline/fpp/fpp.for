C-----------------------------------------------------------------------------
	Program FPP
C-----------------------------------------------------------------------------
C	Program Name:
C		FPP_Pre_Processor
C
C	Program Description:
C	    The FIRAS Preprocessor, FPP, is the first facility in the pipeline
C	    processing of FIRAS data.  Its purpose is to compute the midpoint of
C	    collect time of the interferograms, to find the gain, fake-it
C	    status, and scan mode during the collection of the interferograms,
C	    and to write this information in th science record in the FIRAS
C	    archive.  It will also flag any invalid computed midpoint times,
C	    unknown gains, unknown fake-it statuses, data gaps, and bad
C	    telemetry quality.
C
C	Designer:
C	    Shirley M. Read, STX, January 1989.
C
C	Programmers:
C	    Shirley M. Read and Quoc C. Chung, STX, January 1989.
C
C	Major ReDesign and Programming:
C	    Larry P. Rosen, STX, January-April 1991.
C
C	Change Log since redesign: (changes prior to redesign are not relevant)
CH
CH	17 October 1991, Larry P. Rosen, Hughes STX, SER 9167.
CH	Add check for number of sweeps, and qualifier.  1.  IFG's are always
CH	flagged bad if # sweeps is < 0 or >16.  2. If qualifer true (default)
CH	then IFG's flagged bad if # sweeps not = 0, 1, 4, or 16.
CH
C	Calling Sequence:
C	    This is a main program.
C
C	Input/Output Parameters:
C	    None
C
C	Input Files:
C	    FIRAS Raw Science Data Archive Files : NFS_SDF_RH, NFS_SDF_RL,
C						   NFS_SDF_LH, NFS_SDF_LL
C	    FIRAS WORD 31 Data Archive File : NFS_ANC
C
C	    FIRAS Reference Data Set : FEX_GAIN.DAT, FEX_FAKEIT.DAT
C
C	Output Files:
C	    FIRAS Raw Science Data Archive Files : FPP_SDF_RH, FPP_SDF_RL,
C						   FPP_SDF_LH, FPP_SDF_LL
C	    Report file for run : FPP_Report.Lis_yydddhhmm
C
C	Include Files Used:
C	    CT$Library:CTUser.Inc (COBETRIEVE return status definitions)
C
C	Subroutines and Functions Called:
C	    FPP_Parse_Command
C	    FPP_Init_Report
C	    FPP_Cat_Info
C	    FPP_open_fex
C	    FPP_Open_Archive
C	    FPP_Collect_Time
C	    FPP_Gain
C	    FPP_FakeIt
C	    FPP_Close_Archive
C	    FPP_Check_Segment_Overlap
C	    FPP_Close_Report
C	    CT_Read_Arcv
C	    CT_Write_Arcv
C	    Lib$Establish
C	    Lib$Signal
C	    FUT_GMT_CONV
C
C	Error Handling:
C	    Establish the condition handler. Signal all errors via Lib$Signal.
C	    Error messages to SYS$OUTPUT (terminal if run interactively, log
C	    file if run batch) and report file if requested.
C
C	Method Used:  PDL for FPP_Pre_Processor
C
C (* indicates new feature)
C
C	Establish the condition handler.
C
C	Call FPP_Parse_Command to parse user options from the command line.
C
C	IF (report requested) THEN
C	    CALL FPP_Init_Report to open the report file and write initial
C	    information for the run:
C		Standard banner,
C *		Process invoking facility,
C *		Expansion of logical names,
C *		Invoking command line with defaults,
C		Current GMT time.
C	ENDIF
C
C	CALL FPP_Cat_Info to get the catalog entries for all NFS_SDF and NFS_ANC
C	segments for the run.  Datasets may be selected by time range or by a
C	channel filename.  A check is made for a valid NFS_SDF catalog file if a
C	single channel is selected by filename.  The complete data time range
C	for the run is determined from the catalog information.  Construct and
C	store the input and output data filenames of the channels to be
C	processed and the corresponding NFS_ANC files.
C
C *	Determine number of segments to write (one per day).
C
C *	DO for each segment (day):
C
C	    DO for each channel to be processed:
C
C *		OPEN gain and fake_it reference data files.
C
C		CALL FPP_Open_Archive to open the input data file NFS_SDF, the
C *		ancillary data file NFS_ANC, and the output FPP_SDF file.
C
C		Set the first record flag to true and EOF flag to false.
C
C		Do While ( not EOF ) (for each record in the input files)
C
C		    CALL CT_Read_Arcv to read the NFS_SDF science record.
C
C		    IF ( CT status is not end of file ) THEN
C
C			IF (report requested) THEN report any science data gap.
C
C			CALL FPP_Collect_Time to compute the midpoint of collect
C			time by determining and using the MTM scan speed and
C			length from the NFS_ANC file.  The midpoint time, the
C			badtime_flag, the MTM scan speed and the MTM scan length
C *			will be stored in the science record.  Start and end of
C			collect times will be returned.
C
C *			Call FPP_Gain to determine the gain for the record
C *			   from the info in the reference data file.
C *			Call FPP_FakeIt to determine the fake_it status
C *			    for the record from the info in the reference file.
C			CALL CT_Write_Arcv to write the FPP_SDF record.
C		    ELSE
C			Set EOF to True
C		    ENDIF
C *		    Increment BadTime flag summary.
C		ENDDO
C
C		CALL FPP_Close_Archive to close the NFS_SDF, NFS_ANC, FPP_SDF
C  *		archived files.  Close gain and fake_it reference data files.
C
C	    ENDDO (for each channel to be processed)
C
C	ENDDO (for each segment)
C
C	CALL FPP_Check_Segment_Overlap to report FPP overlaps in archive.
C
C	IF (Report) Then
C	    Write and close report file.
C	ENDIF
C
C	Signal the processing status and Exit.
C-----------------------------------------------------------------------------
	Implicit None
	
	External	Fut_Error	! condition handler
	
C  Include Files

	Include		'CT$Library:CTUser.Inc'
	Include		'($SSDef)'
	Include		'(FUT_Params)'
	Include 	'(FPP_Msg)'

C  Record Structures

	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ Sci_Rec
	
	Dictionary 'NFS_ANC'
	Record /NFS_ANC/ Anc_Rec

C  Functions

	Integer*4 FPP_Parse_Command
	Integer*4 FPP_Init_Report
	Integer*4 FPP_Cat_Info
	Integer*4 FPP_Open_Fex
	Integer*4 FPP_Open_Archive
	Integer*4 FPP_Collect_Time
	Integer*4 FPP_Gain
	Integer*4 FPP_Fakeit
	Integer*4 FPP_Close_Archive
	Integer*4 FPP_Check_Segment_Overlap
	Integer*4 FPP_Close_Report
	Integer*4 Cut_Register_Version
	Integer*4 Cut_Display_Banner
	Integer*4 FUT_GMT_Conv

C  Local Variables

	Character*6	Version
	Parameter	(Version='8.8')
	Character*15	Arch_In  /'CSDR$FIRAS_RAW:'/	! input archive
	Character*15	Arch_Out /'CSDR$FIRAS_OUT:'/	! output archive
	Character*15	Arch_Ref /'CSDR$FIRAS_REF:'/	! reference archive
	Integer*2	Maxn		! max number of output files
	Parameter	(Maxn=50)
	Integer*4	Status		! status of processing
	Integer*4	Rstatus		! return status from function call
	Integer*4	Success / 1 /, Error / 2 /  ! values for status
	Integer*4	Time_Present	! presence of /jstart, /jstop
	Character*14	Gmt_Start	! requested start gmt for data
	Character*14	Gmt_Stop	! requested stop gmt for data
	Integer*4	File_Present	! presence of /filename
	Character*39	File_Name	! user specified file name
	Logical*1	Report		! flag to enable writing a report file
	Integer*4	Max_Sec		! max # secs betw. mjfs for gap declared
	Integer*4	Anc_Gap		! max # secs betw. mjfs to report lg gap
	Integer*4	Sci_Gap		! max # secs betw. nfs_SDF recs for gap
	Logical*1	Sweep_Check	! Check the number of sweeps.
	Integer*4	Lun_Anc		! unit numbers for nfs_anc ct archive
	Integer*4	Lun_In_Sci	! unit numbers for nfs_SDF ct archive
	Integer*4	Lun_Out_Sci	! unit numbers for FPP_SDF ct archive
	Integer*4	Lun_Rpt		! report unit number
	Integer*4	Lun_Fake	! unit number for reference fakeit file
	Integer*4	Lun_Gain	! unit number for reference gain file
	Integer*2	Num_Chan_Tbl(Maxn)	! # chan segmnts to be processed
	Integer*2	Proc_Chan_Tbl(Maxn,4)	! table of channel segments
	Integer*2	Num_Chan	! number of channels to process
	Integer*2	Proc_Chan(4)	! channels to be processed
	Integer*2	Day /0/		! current daily output file number
	Integer*4	Num_Days	! number of daily output files
	Character*64	In_Files(Maxn,5)	! names of files to be processed
	Character*64	Out_Files(Maxn,4)	! names of files to be created
	Integer*4	Num_Rec(Maxn,4)	! number of records in each file
	Integer*4	Num_Good(Maxn,4)	! # of good records in each file
	Integer*4	Num_Bad(Maxn,4)	! number of bad records in each file
	Integer*4	Badflg(18,Maxn,4)	! counter for badtime flags
	Integer*2	Ct_Stat(20)	! cobetrieve status
	Integer*2	Chan_Num	! channel number
	Integer*2	Last_Chan / 0 / ! previous channel processed
	Integer*2	Chan		! channel index
	Logical*1	Eof		! flag for eof on nfs_SDF file
	Logical*1	Anc_Eof		! flag for eof on nfs_anc file
	Logical*1	First_Rec	! first reocrd in science file
	Integer*2	Bf		! current rec badtime flag (0-18)
	Integer*4	Check_Jstart(2)	! overlap check base start time
	Integer*4	Check_Jstop(2)	! overlap check base stop time
	Integer*4	Num_Vol/80/
        Integer*4	Lun_Out/6/
	Character*33	Report_File	! file name for report
	Character*79	Command_Line(3)	! command line with defaults
        Real*8		Tdiff
        Real*8		Prev_Time
        Real*8		Time_Sec
        Character*14	Prev_Gmt
        Character*2	Chan_Name(4)/'RH','RL','LH','LL'/
        Integer*2	Ihr, Imin		! gmt date/time parts
	Integer*4	Msec			! gmt date/time part
	Logical*1	Sync			! mtm in sync with micro
	Integer*4	Start_Collect(2)	! start time of science collect
	Integer*4	End_Collect(2)		! end time of science collect

C  Set status for FPP processing to Success.

	Status = Success

C  Establish condition handler.

	Call Lib$Establish ( Fut_Error )

	Rstatus = Cut_Register_Version (Version)
	Rstatus = Cut_Display_Banner (Lun_Out, Num_Vol,
	1   'Firas Facility FPP_Pre_Processor')
	Write ( Lun_Out, 10 )
  10	Format (//)

C  Get processing options from command line.

	Rstatus = FPP_Parse_Command ( Time_Present, Gmt_Start, Gmt_Stop,
	1   File_Present, File_Name, Proc_Chan, Num_Chan, Max_Sec, Report,
	2   Report_File, Command_Line, Sci_Gap, Anc_Gap, Sweep_Check )
	If ( Rstatus .Ne. %Loc(FPP_Normal) ) Status = Error

C  Open a report file for the run and write initial information.

	If ( (Status .Eq. Success) .And. (Report) ) Then
	  Rstatus = FPP_Init_Report ( Lun_Rpt, Report_File, Command_Line,
	1    Version, Arch_In, Arch_Out, Arch_Ref )
	  If ( Rstatus .Ne. %Loc(FPP_Normal) ) Status = Error
	Endif

C  Get the catalog information for NFS_SDF and NFS_ANC records to be processed.

	If ( Status .Eq. Success) Then
	  Rstatus = FPP_Cat_Info ( Time_Present, Gmt_Start, Gmt_Stop,
	1   File_Present, File_Name, Proc_Chan, Num_Chan, Report, Lun_Rpt,
	2   Maxn, Num_Chan_Tbl, Proc_Chan_Tbl, Num_Days, In_Files, Out_Files,
	3   Num_Rec, Check_Jstart, Check_Jstop, Arch_In, Arch_Out )
	  If ( Rstatus .Ne. %Loc(FPP_Normal) ) Status = Error
	Endif

C  DO for each output file (one per day per channel)

	Do While ((Day .Lt. Num_Days) .And. (Status .Eq. Success))
	  Day = Day + 1

C  Process channels, one by one.

	  Do Chan = 1, Num_Chan_Tbl(Day)

	    Chan_Num = Proc_Chan_Tbl(Day,Chan)

C  Open the NFS_ANC, NFS_SDF, and FPP_SDF files.

	    If ( Status .Eq. Success ) Then
	      Rstatus = FPP_Open_Archive ( Chan_Num, In_Files, Out_Files,
	1	Report, Lun_Rpt, Lun_Anc, Lun_In_Sci, Lun_Out_Sci, Day, Maxn )
	      If ( Rstatus .Ne. %Loc(FPP_Normal) ) Then
		Status = Error
		Eof = .True.
	      Else
		Eof = .False.
		Anc_Eof = .False.
		First_Rec = .True.
              Endif
	    Endif

C  Open the gain and fake-it reference data files.

	    If ( Status .Eq. Success) Then
	      Rstatus = FPP_Open_Fex ( Arch_Ref, Gmt_Start, Gmt_Stop, Lun_Rpt,
	1         Report, Lun_Fake, Lun_Gain, Chan_Num )
	      If ( Rstatus .Ne. %Loc(FPP_Normal) ) Status = Error
	    Endif

C  Process the science records one by one until all are done (Do for each rec).

	    Do While ((.Not. Eof) .And. (Status .Eq. Success))

	      Call Ct_Read_Arcv ( , Lun_In_Sci, Sci_Rec, Ct_Stat )

	      If (Ct_Stat(1) .Ne. Ctp_Normal) Then
		If ( Ct_Stat(1) .Eq. Ctp_Endoffile ) Then
		  Eof = .True.
		Else		! read error
		  Status = Ct_Stat(1)
		  Call Lib$Signal (FPP_Ctreaderr,%Val(2),%Val(Status),
	1	    In_Files(Day,Chan_Num))
		  Status = Error
		Endif		! read error
	      Else		! normal return from ct_read_arcv

	        If (Time_Present .EQ. Fac_Present) Num_Rec (Day, Chan_Num) =
	1          Num_Rec (Day, Chan_Num) + 1

		If (Report) Then

C  Check for large science gap and report it.

		  RStatus = FUT_Gmt_Conv (Time_Sec, Sci_Rec.Ct_Head.Gmt, 2)
		  If (.Not. First_Rec) Then
		    Tdiff = Time_Sec - Prev_Time
		    If (Tdiff .Gt. Sci_Gap) Then
	              Ihr = Int (Tdiff / Fac_Hour)
	              Tdiff = Tdiff - Dble(Ihr*Fac_Hour)
	              Imin = Int (Tdiff / Fac_Minute)
	              Tdiff = Tdiff - Dble(Imin*Fac_Minute)
	              Msec = Int (Tdiff)
		      Write (Lun_Rpt,30) Chan_Name(Chan_Num), Prev_Gmt,
	1		Sci_Rec.Ct_Head.Gmt, Ihr, Imin, Msec
  30		      Format (5X,'Gap: NFS_SDF_',A2,5X,A14,' To ',A14,' = '
	1                 I2.2,' Hr ',I2.2,' Min ', I2.2,' Sec')
		    Endif
		  Endif
		  Prev_Gmt = Sci_Rec.Ct_Head.Gmt
		  Prev_Time = Time_Sec
		Endif				! if report

C  Compute the midpoint of collect time for the science record.

		If (Status .Ne. Error) Then

	   	  Rstatus = FPP_Collect_Time ( Chan_Num, Max_Sec, Sci_Rec,
	1	    Lun_Anc, First_Rec, Report, Lun_Rpt, Anc_Gap, Anc_Eof,
	2	    Sync, Start_Collect, End_Collect, Chan, In_Files(Day,5),
	3	    Sweep_Check )
		  If ( Rstatus .Ne. %Loc(FPP_Normal)) Then
		    If (Rstatus .Eq. %Loc(FPP_Ctanceof)) Then
		      Eof = .True.
		    Else
		      Status = Error
		    Endif
	          Else

C  Get the gain for this channel and record.

		    Rstatus = FPP_Gain(Chan_Num, Sci_Rec, First_Rec, Lun_Gain,
	2	      Report, Lun_Rpt, Sync, Start_Collect, End_Collect )
	
C  Get the fake-it status for this channel and record.

		    Rstatus = FPP_Fakeit ( Chan_Num, Sci_Rec, First_Rec,
	2	      Lun_Fake, Report, Lun_Rpt, Start_Collect, End_Collect )
	
C  Write the FPP_SDF record.

		    Call Ct_Write_Arcv ( , Lun_Out_Sci, Sci_Rec, Ct_Stat )
		    If (Ct_Stat(1) .Ne. Ctp_Normal) Then
		      Status = Ct_Stat(1)
		      Call Lib$Signal (FPP_Ctwriterr, %Val(2), %Val(Status),
	1		'FPP_SDF_'//Fac_Channel_Ids(Chan_Num))
		      Status = Error
		    Endif		! write error
		  Endif			! good return from FPP_collect_time
		Endif			! status ok
	      Endif			! normal return from ct_read_arcv

C  Increment the the good or bad record counter and BadTime Flags counter.

	      If (Status .Ne. Error .And. .Not.Eof) Then
	         Bf = Sci_Rec.Collect_Time.Badtime_Flag
	         If ( Bf .Eq. 0 ) Then
		   Num_Good (Day, Chan_Num) = Num_Good (Day,Chan_Num) + 1
	         Else
		   Num_Bad (Day, Chan_Num) = Num_Bad (Day, Chan_Num) + 1
		   Badflg (Bf,Day,Chan_Num) = Badflg(Bf,Day,Chan_Num) + 1
	         Endif
	         First_Rec = .False.
	      Endif
	    Enddo		! for each record

C  Close the NFS_SDF and NFS_ANC Archive Files.

	    If ( Status .Eq. Success ) Then
	      Rstatus = FPP_Close_Archive ( Chan_Num, Report, Lun_Rpt, Lun_Anc,
	1	Lun_In_Sci, Lun_Out_Sci, Lun_Fake, Lun_Gain )
	      If ( Rstatus .Ne. %Loc(FPP_Normal)) Status = Error
	    Endif

	  Enddo					! for each channel
	Enddo 					! for each daily output file
	
C  Perform segment checks: FPP_SDF_xx overlaps.

	If ( Status .Eq. Success ) Then
	  Rstatus= FPP_Check_Segment_Overlap ( Proc_Chan, Num_Chan,
	1   Check_Jstart, Check_Jstop, Report, Lun_Rpt )
	  If ( Rstatus .Ne. %Loc(FPP_Normal)) Status = Error
	Endif

C  Close the FPP_Report File.

	If ( Report ) Then
	  Rstatus = FPP_Close_Report ( Maxn, Num_Chan_Tbl, Proc_Chan_Tbl,
	1   Num_Days, Out_Files, Num_Rec, Num_Good, Num_Bad, Badflg, Lun_Rpt,
	2   Status )
	  If ( Rstatus .Ne. %Loc(FPP_Normal)) Status = Error
	Else

C  Signal the processing status.

	   If (Status .Eq. Success) Then
	      Call Lib$Signal (FPP_Normal)
	   Else
	      Call Lib$Signal (FPP_Aberr)
	   Endif
	Endif

C  Exit with status.

	If (Status .Eq. Success) Then
	  Call Exit (Ss$_Normal)
	Else
	  Call Exit (Ss$_Abort)
	Endif

	End
