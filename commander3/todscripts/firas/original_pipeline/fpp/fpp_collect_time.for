C-------------------------------------------------------------------------------
	Integer*4 Function FPP_Collect_Time ( Chan_Num, Max_Sec, Sci_Rec,
	1   Lun_Anc, First_Rec, Report, Lun_Rpt, Anc_Gap, Anc_Eof, Sync,
	2   Start_collect, End_collect, Gap_Check, Ancfilename, Sweep_Check)
C-------------------------------------------------------------------------------
C	Purpose: To obtain the MTM scan speed and length from the NFS_ANC 
C		 record and call the routine to compute the midpoint of 
C	         collect time. If the midpoint time cannot be computed or
C		 the computed midpoint time is unvalid, the badtime flag
C		 is set with a code value indicating the reason for failure.
C
C	Author: Shirley M. Read
C		STX, January, 1989
C
C	Invocation: Status = FPP_Collect_Time ( Chan_Num, Max_Sec, Sci_Rec,
C	1   Lun_Anc, First_Rec, Report, Lun_Rpt, Anc_Gap, Anc_Eof, Sync,
C	2   Start_collect, End_collect, Gap_Check, Ancfilename, Sweep_Check)
C
CH  Change Log:
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH	There have been problems during test operations for FIRAS processing due
CH	to the required clean-up of the archives after an FPP or FDQ abort.  The
CH	FPR tracking system compounds the problems Files with non-matching
CH	version numbers seem often to result from improper clean-up. Bad record
CH	times cause SEGCTL to abort and mess up the tracking system. It was
CH	decided to change the modify of the science records in FPP and FDQ to a
CH	simple COBETRIEVE read of the existing records from a dataset and write
CH	a modifed dataset with the same information which was entered on the
CH	modify. Two new science data sets will be required: a science dataset of
CH	raw science data plus FPP input and a science dataset with FPP science
CH	data plus FDQ input. These datasets will be FPP_SDF_xx, where xx is the
CH	channel id (RH,RL, LH or LL) and FDQ_SDF_xx, where xx is the channel id.
CH	The new datasets must be opened and processed in FPP and FDQ.
CH 
CH	SPR 6738, H. Wang, STX, 6/29/90
CH         FPP collect time filled with Maj Frm Time when badtime flag set.
CH
CH	New Requirements, Larry Rosen, STX, 14 March 1991
CH	   Major changes as per new design.  See new PDL below.
CH
CH	17 October 1991, Larry P. Rosen, Hughes STX, SER 9167.
CH	Add check for number of sweeps, and qualifier.  1.  IFG's are always
CH	flagged bad if # sweeps is < 0 or >16.  2. If qualifer true (default)
CH	then IFG's flagged bad if # sweeps not = 0, 1, 4, or 16.
CH	Also, if the # sweeps = 0, then calculate the midpoint time as if it
CH	there were 1 sweep.
C ------------------------------------------------------------------------------
C	Input Parameters:
C	  Name		  Type		  Description
C 	  ----------------------------------------------------------------------
C	  Chan_Num        I*2             Channel number
C	  Max_Sec         I*4             Maximum number of seconds between 
C				          major frames before a gap is declared
C	  Sci_Rec    NFS_SDF Rec          Current FIRAS Raw Science Record
C	  Lun_Anc         I*4             Logical unit for NFS_ANC records
C	  First_Rec       L*1             Flag for first record in a file
C	  Report          L*1             Flag to write report
C	  Lun_Rpt         I*4             Logical unit for report
C	  Anc_gap	  I*4             Max number seconds for large gap reprt
C	  Gap_Check	  I*2		  Flag to check for sci and anc gaps.
C	  Ancfilename	  C*64		  Anc filename gives times for ref data.
C	  Sweep_Check     L*1		  Flag to check # sweeps
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Anc_Eof         L*1             End of file on NFS_ANC
C	  Sync            L*1             MTM sync flag
C	  Start_Collect   I*4(2)          Start time of science collection
C	  End_Collect     I*4(2)          End time of science collection
C	
C	Subroutines Called:
C	  COBETRIEVE: CT_GMT_to_Binary, Time_LT, Time_LE, Time_GE, Time_GT
C         Extended Precision Arithmetic: Lib$Emul, Lib$Addx, Lib$Subx, Lib$Ediv
C	  Lib$Signal
C
C	Include Files:
C	  FPP_Msg.Txt
C         $SSDef
C	  CT$Library:CTUser.Inc
C	  Fut_Params
C-------------------------------------------------------------------------------
C	New PDL for  FPP_Collect_Time               Larry P. Rosen, Feb,Apr 1991
C
C	(* indicates new feature)
C	(Note: Two buffers are used for ancillary records.  Each ancillary
C	record is read into buffer two, where the comparison to the science
C	record's reference time is made.  If the science reference time is
C	beyond the ancillary record in buffer two, that anc. record is copied
C	to buffer one, and the next anc. record is read into buffer two.)
C
C	IF (science record is "First Record" in the file) THEN
C
C  *	   CALL FPP_Read2 to read the first NFS_ANC record, where at least one
C	      frame has science telemetry format, into the current buffer (2).
C	   IF (error) Signal an error and set an error status for the return.
C	ENDIF
C
C  *	Initialize midpoint, start, and end of collect time as transmit time.
C
C  *	IF (# Sweeps in IFG < 0 or > 16) Then
C	   Set BadTime Flag = 18
C	   Save Transmit time for next call.
C	   Return
C  *	IF (Sweep_Check flag is True) Then
C  *	   IF (# Sweeps NE 0, 1, 4, or 16) Then
C	      Set BadTime Flag = 18
C	      Save Transmit time for next call.
C	      Return
C	   EndIf
C	EndIf
C
C  *	Determine MTM sync (in or out) from science record (status word 3,
C  *	   bit 10, vax bit 5); 0 = in sync.
C
C	Transmit time from science record = sci_rec.ct_head.time.
C
C  *	IF (in sync) THEN
C  *	   Get and check minor frame counters for collect and transmit.
C  *	   Difference of counters times sec_per_minor_frame gives Delta Time.
C  *	   Subtract Delta Time from Transmit time to get Start Collect.
C  *	   Add 1/2 second to Start Collect to get Ref_Time.
C  *	ELSE
C  *	   Subtract 1 sec_per_minor_frame from transmit time to get both
C  *	   End_Collect and Ref_Time.
C  *	ENDIF
C
C  *	Extract science record's adds per group from science header (sc_head9).
C  *	IF (adds_per_group .EQ. 0) THEN
C  *	   Set adds_per_group = 1
C  *	ELSE IF (adds_per_group .NE. 1 OR 2 OR 3 OR 8 OR 12) THEN
C  *	   Set BadTime flag in science record to seventeen.
C  *	ENDIF
C
C	IF (status ok) THEN
C	   Set the Compute flag to TRUE.
C	ELSE
C	   Set Compute flag to FALSE and Computation error to TRUE.
C	END IF
C
C  *	Determine ancillary major frame, minor frame, and extract scan mode:
C  *	See walkthrough document for figures of various cases and flow diagram.
C	Do (Until correlated: scan mode of ref_time found or badtime flag set)
C	   If (Telemetry Format Record 2, Frame 1 is good) Then
C	      If (Ref_Time LT Record 2, Frame 1 time) Then
C	         If (Telemetry Format Record 1, Frame 2 is good) Then
C	            If (Gap [Rec1,Frm2]->[Rec2,Frm1])
C	               Case 1. Ref_Time is in gap.
C	               Set Flag = 9
C	            Else
C	               Case 3. Ref_Time is between the last minor frame of the
C	               first ANC record and the start time of the second ANC
C	               record.  Index = 128.  Extract MTM speed (Vax bit 4,
C	               status word 31 (Index)).  Extract MTM length (Vax bit 5).
C	            EndIf
C	         ElseIf (Anc_Counter .LE. 1) Then
C	            Ref_Time is in gap before first frame.
C	            Set Flag = 1
C	         Else  (TF rec 1 Frm 2 is bad)
C	            Case 4. Ref_Time points to frame with non-science telemetry
C	            format.  Set flag = 16
C	         EndIf
C	      Else		! Ref time > MF1 time
C	         If (Ref_Time LT Last_MNF_Time of frame 1) Then
C	            Cases 5 and 6. Ref_Time points in frame 1. Compute the index
C	            to the NFS_ANC Word 31 minor frame status bit array.  Set
C	            MNF index to (Ref_Time - MF_1_Start_Time) * 128 frames/32.
C	            Since MNF is 0 to 127 and record array goes 1 to 128, add 1.
C	            Extract speed and length.
C	         Else
C	            If (Telemetry Format Record 2 Frame 2 is bad) Then
C	               Case 8. Ref_Time is after frame 1 and frame 2 is bad. Get
C	               next 2 frames where at least one frame is good Telemetry
C	               format.  Call function FPP_Read2.
C	            Else
C	               If (Ref_Time  LT  Anc_Rec(2).GMT_MJF2) ) Then
C	                  If (Gap) Then
C	                     Case 9. Ref_Time points to a gap between major
C	                     frames 1 and 2.  Set flag = 9.
C	                  Else
C	                     Case 12. Ref_Time is between the last minor frame
C	                     of major frame 1 and the the start time of major
C	                     frame 2 in the second ANC record.
C	                  EndIf
C	               Else
C	                  If (Ref_Time  LT  Last_MNF_Time2) Then
C	                     Case 10. Ref_Time points in frame 2. Compute the
C	                     index to the NFS_ANC Word 31 minor frame status
C	                     bit array.  Extract speed and length.
C	                  Else				! Ref > End of MJF2
C	                     Case 11.  Ref_Time points after end of Major frame
C	                     2.  Get next record that contains at least one
C	                     with science format.  Call FPP_Read2.
C	                  EndIf
C	               EndIf
C	            EndIf
C	         EndIf
C	      EndIf
C	   Else			! Telemetry format Frame 1 Rec 2 is bad
C	      If (Ref_Time  LT  Anc_Rec(2).GMT_MJF2) Then
C	         Cases 2 and 7.  Ref_Time points prior to Frame 2 and frame 1
C	         has non-science telemetry format. Check if this if 1st Anc rec.
C	         If (Anc_Counter = 1) Then
C	            Set flag = 1
C	         Else
C	            Set flag = 16
C	         EndIf
C	      Else
C	         If (Ref_Time  LT  Last_MNF_Time2) Then
C	            Case 10. Ref_Time points in frame 2. Compute the index to
C	            the NFS_ANC Word 31 minor frame status bit array.  Extract
C	            speed and length.
C	         Else				! Ref > End of MJF2
C	            Case 11.  Ref_Time points after end of Major frame 2.  Get
C	            next record that contains at least one with science format.
C	            Call FPP_Read2.
C	         EndIf
C	      EndIf
C	   EndIf
C	EndDo	! Do until correlated : when scan mode of ref_time found or flag
C
C	IF ("Compute Flag" is set AND status is good) THEN
C	   Call FPP_Midpoint_Time to compute the midpoint of collect time using
C	   the MTM_Scan_Speed and MTM_Scan_Length.  If the #sweeps = 0, use 1.
C  *	    IF ("in sync") THEN
C  *	       Also return End of Collect Time.
C  *	    EL("out of sync")
C  *	       Also return Start of Collect Time.
C  *	    ENDIF
C	    IF (return status = 1) THEN
C		Set the "Midpoint Time" in the science record to the returned
C		midpoint of collect time.
C	        Call FPP_Check_Time for the consistency check of the MTM scan
C	        speed and length for all minor frames in the collection.
C	    ENDIF
C	ENDIF
C
C	Save the transmit time from the current science record for comparison
C	with the next science record; sci_rec.ct_head.time
C
C	RETURN with normal or error status.
C-------------------------------------------------------------------------------
C			FPP BADTIME FLAG VALUES
C
C	(* indicates new flags)
C
C	    Zero  -- No error, normal.
C  a	    One   -- Telemetry gap before IFG record
C  a	    Two   -- IFG  transmitted before collected
C  a	    Three -- IFG transmit frame count is less than zero
C  a	    Four  -- IFG collect frame count is less than zero
C  a	    Five  -- IFG transmit frame equals the collect frame
C  c	    Six   -- Midpoint of collect > transmit time
C  c	    Seven -- Midpoint of collect < previous transmit time
C  c	    Eight -- Consistency check for scan mode failed
C  c	    Nine  -- Telemetry gap for consistency check
C  a	    Ten   -- Computational error during calculation of collect time
C  d*	    Eleven - Gain change and out of sync -> gain unknown
C  e*	    Twelve - Fake-it change -> unknown fake-it status
C  a*	    Thirteen - Science time-out (gap)
C  d&e*	    Fourteen - Housekeeping data missing (gap)
C  d&e*	    Fifteen  - Telemetry Quality bad
C  a*	    Sixteen  - Ancillary gap due to non-science telemetry format
C  a*	    Seventeen - Bad number Adds/Group in micro processor header
C  a*	    Eighteen - Bad number of Sweeps
C
C	a = set in FPP_Collect_Time
C	c = set in FPP_Check_Time
C	d = set in FPP_Gain
C	e = set in FPP_FakeIt
C*******************************************************************************
	Implicit None

C  Passed Parameters.

	Integer*2	Chan_Num        ! Channel number
	Integer*4	Max_Sec         ! Maximum number of seconds between 
				        ! major frames before a gap is declared
	Dictionary	'NFS_SDF'
	Record / NFS_SDF / Sci_Rec      ! Current FIRAS Raw Science Record

	Integer*4	Lun_Anc         ! logical unit for NFS_ANC records
	Logical*1	First_Rec       ! Flag for first record in a file
	Logical*1	Report		! Flag to write report
	Integer*4	Lun_Rpt		! Unit number for report
	Logical*1	Anc_Eof		! End of file on current NFS_ANC segment
	Integer*4	Anc_Gap		! max seconds for report of gaps in anc
	Logical*1	Sync		! MTM sync flag
	Integer*4	Start_Collect(2)  ! Start time of science collection
	Integer*4	End_Collect(2)	! End time of science collection
	Integer*2	Gap_Check	! Flag to look for gaps in sci and anc
	Character*64	AncFilename	! Gives times for reading ref data.
	Logical*1	Sweep_Check	! Flag to check # Sweeps.

C  Include files.

	Include		'(FPP_Msg)'
	Include		'(FUT_Params)'
	Include		'($SSDef)'
	Include		'CT$Library:CTUser.Inc'

C  Functions: Extended Precision Arithmetic, COBETRIEVE Time Comparison, FPP

	Integer*4	Lib$Emul
	Integer*4	Lib$Subx
	Integer*4	Lib$Addx
	Integer*4	Lib$Ediv
	Logical*1	Time_LT
	Logical*1	Time_LE
	Logical*1	Time_GE
	Logical*1	Time_GT
	Integer*4	FPP_Read2
	Integer*4	FPP_Check_Time
	Integer*4	FPP_Midpoint_Time

C  Local Parameters.

	Integer*4	Sec_per_MNF     ! Seconds per minor frame (0.25)
	Parameter	( Sec_per_MNF = 2500000 )
	Integer*4	Quart_Sec(2)
	Data		Quart_Sec(1) / 2500000 /
	Integer*4	One_Sec         ! One second
	Parameter	( One_Sec = 10000000 )
	Integer*4	Half_Sec(2)	! (0.5 sec)
	Data		Half_Sec(1) / 5000000 /

C  Local Declarations.

	Logical*1	Comperr			! Arithmetic computation error 
	Logical*1	First_Time/.True./	! First time
	Integer*4	Flag_Time(2)	! Quadword for earliest CT date known
	Integer*4	Addend / 0 /	! Used in Lib$Emul
	Integer*4	Status		! Return status for functions
	Integer*4	QMax_Sec(2)     ! Vax Quadword for Max_Sec input parm
	Logical*1	Compute		! Flag to proceed with computation of 
					! midpoint of collect
	Integer*4	Ntime/0/
	Integer*4	Prev_Time(2)	! Previous IFG quadword time
	Integer*2	Sweeps
	Logical*1	First_Anc	! Flag first read of Anc for the day
	Logical*1	Bad_Ancil		! telemetry formats of Ancil rec
	Dictionary	'NFS_ANC'
	Record / NFS_ANC / Anc_Rec(2)   ! Current FIRAS Word 31 Records
	Integer*2	CT_Stat(20)     ! COBETRIEVE return status
	Integer*4	Anc_Counter / 0 / ! NFS_ANC record counter
	Logical*1	TF_R1(2)	! Flag Telemetry format each frame, Rec1
	Logical*1	TF_R2(2)	! Flag Telemetry format each frame, Rec2
	Character*14	Last_gmt	! Last frame time with science tlm fmt
	Integer*4	Last_MNF_Time1(2)  ! Last minor frame time, Major Frame1
	Integer*4	Last_MNF_Time2(2)  ! Last minor frame time, Major Frame2
	Integer*4	MJF_Diff(2)     ! Difference between the 2 MJ F times.
	Integer*4	IFG_Time(2)	! Computed IFG quadword time
	Integer*4	Len / 2 /       ! Length for Lib$Subx & Lib$Addx array 
	Integer*4	I4_Bits		! Longword for storing status bits
	Logical*1	Correlated	! Flag indicating reference time is 
					! within an NFS_ANC frame time
	Integer*2	Adds			! Adds per group for IFG
	Integer*4	Transmit_Frame
	Integer*2	Transmit_Word(2)
	Equivalence	(Transmit_Frame, Transmit_Word(1))
	Integer*4	Collect_Frame
	Integer*2	Collect_Word(2)
	Equivalence	(Collect_Frame, Collect_Word(1))
	Integer*4	Delta_Time(2)	! Binary frame time difference
	Integer*4	Ref_Time(2)	! Time for obtaining MTM scan mode
	Integer*4	Del_Ref(2)      ! Difference between minor frame times
	Integer*4	Remainder	! Remainder in division function
	Integer*4	Index           ! Minor frame index to status bits array
	Integer*4	Zero / 0 /      ! Zero value
	Integer*4	One / 1 /       ! One value
	Integer*4	Two / 2 /       ! Two value
	Integer*4	Three / 3 /     ! Three value
	Integer*4	Four / 4 /      ! Four value
	Integer*4	Five / 5 /      ! Five value
	Integer*4	Nine / 9 /      ! Nine value
	Integer*4	Ten / 10 /      ! Ten value
	Integer*4	Sixteen / 16 /	! Sixteen value
	Integer*4	Seventeen / 17 /	! Seventeen value
	Integer*4	Eighteen / 18 / ! Eighteen

C  Set the function status to Normal.

	FPP_Collect_Time = %loc(FPP_Normal)
	Status = -1

C  Compute the earliest CT time for flag and quadword maximum secs for ANC gap.

	If ( First_Time ) Then
	   Call CT_GMT_to_Binary (fac_jstart_default, Flag_Time)
	   Status = Lib$Emul (Max_Sec, One_Sec, Addend, QMax_Sec)
	   If (.Not. Status) Then
	      Call Lib$Signal(FPP_EmulErr, %val(1), %val(Status))
	      FPP_Collect_Time = %loc(FPP_Aberr)
	   Endif
	   First_Time = .False.
	Endif	

C  Allow the program to perform the first record block even if there is a bad
C  processing status to catch additional errors for post mortem analysis.

	If (.Not. Status) Then
	   Comperr = .True.
	   Compute = .False.
	Else
	   Comperr = .False.
	   Compute = .True.
	Endif

C  If the science record is the "First Record" in the file, then call
C  CT_Read_Arcv to read the first NFS_ANC record into current buffer(2).
C  Also, initialize the previous record time to the CT flag time.

	If ( First_Rec ) Then
	   Prev_Time(1) = Flag_Time(1)
	   Prev_Time(2) = Flag_Time(2)
	   First_Anc = .True.
	   status = FPP_Read2 ( Lun_Anc, Anc_Rec, TF_R2, TF_R1, Gap_Check,
	1     Anc_Gap, Report, Lun_Rpt, Last_MNF_Time1, Last_MNF_Time2,
	2     First_Anc )
	   If (status .NE. %loc(FPP_Normal)) Then
	      FPP_Collect_Time = status
	      Correlated = .True.
	      Compute = .False.
	   EndIf
	   Anc_Counter = One
	   First_Anc = .False.
	EndIf 				! First_Rec

C  Initialize Midpoint, Start, and End of Collect time as the Transmit time.

	Sci_Rec.Collect_Time.Midpoint_Time(1) = Sci_rec.ct_head.time(1)
	Sci_Rec.Collect_Time.Midpoint_Time(2) = Sci_rec.ct_head.time(2)
	Start_Collect(1) = Sci_rec.ct_head.time(1)
	Start_Collect(2) = Sci_rec.ct_head.time(2)
	End_Collect(1) = Start_Collect(1)
	End_Collect(2) = Start_Collect(2)
	Sync = .True.

C  Check Number of sweeps.

	Sweeps = Sci_Rec.Sci_Head.SC_Head11
	IF (Sweeps .LT. 0 .OR. Sweeps .GT. 16) THEN
	   Sci_Rec.Collect_Time.Badtime_Flag = eighteen
	   Prev_Time(1) = Sci_Rec.CT_Head.Time(1)
	   Prev_Time(2) = Sci_Rec.CT_Head.Time(2)
	   Return
	EndIf
	IF (Sweep_Check .EQ. .True.) Then
	   IF (Sweeps .NE. 0 .AND. Sweeps .NE. 1 .AND. Sweeps .NE. 4 .AND.
	1     Sweeps .NE. 16) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = eighteen
	      Prev_Time(1) = Sci_Rec.CT_Head.Time(1)
	      Prev_Time(2) = Sci_Rec.CT_Head.Time(2)
	      Return
	   ENDIF
	ENDIF

C  Determine MTM Sync from science record (statwd3,bit10 (vaxbit5)) 0=in sync.
C  True = in sync

	Sync = .Not. Btest (sci_rec.sci_head.sc_head3, 5)

	If (Sync) Then				! In Sync

C  Get minor frame counters from collect to transit.

	   Transmit_Word(1) = Sci_Rec.Sci_Head.SC_Head4
	   Transmit_Word(2) = Sci_Rec.Sci_Head.SC_Head5
	   Collect_Word(1) = Sci_Rec.Sci_Head.SC_Head12
	   Collect_Word(2) = Sci_Rec.Sci_Head.SC_Head13
	   If (Transmit_Frame .Lt. Collect_Frame) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = two
	   ElseIf (Transmit_Frame .Lt. 0) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = three
	   ElseIf (Collect_Frame .Lt. 0) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = four
	   ElseIf (Transmit_Frame .Eq. Collect_Frame) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = five
	   EndIf

	   If (Sci_rec.Collect_Time.Badtime_Flag .Eq. zero) Then

C  The goal here is to compute the delta time between Collect_frame and
C  Transmit_frame in VAX quadword format.  Since the transmit frame always
C  follows the collect frame, the difference will always be non-negative.
C  Subtracting the delta time from the transmit time gives the reference time.
C  Transmit time of science record is sci_rec.ct_head.time

	      Status = Lib$Emul (Transmit_Frame - Collect_Frame, Sec_Per_MNF,
	1        Addend, Delta_Time)
	      If (.Not. Status) Then
	         Call Lib$Signal (FPP_EmulErr, %val(1), %val(Status))
	         Sci_Rec.Collect_Time.Badtime_Flag = ten
	         Comperr = .True.
	         Compute = .False.
	      Else
	         Status = Lib$Subx (sci_rec.ct_head.time, Delta_Time,
	1           Start_Collect, Len)
	         If (.Not. Status) Then
	            Call Lib$Signal(FPP_SubxErr, %val(1), %val(Status))
	            Sci_Rec.Collect_Time.Badtime_Flag = ten
	            Comperr = .True.
	            Compute = .False.
	         Else
	            Status = Lib$Addx (Start_Collect, Half_Sec, Ref_Time, Len)
	            If (.Not. Status) Then
	               Call Lib$Signal(FPP_AddxErr, %val(1), %val(Status))
	               Sci_Rec.Collect_Time.Badtime_Flag = ten
	               Comperr = .True.
	               Compute = .False.
	            EndIf
	         EndIf
	      EndIf
	   Else
	      Comperr = .True.
	      Compute = .False.
	   EndIf
	Else				! Out of sync
	   Status = Lib$Subx (sci_rec.ct_head.time, Quart_sec ,End_Collect,Len)
	   If (.Not. Status) Then
	      Call Lib$Signal(FPP_SubxErr, %val(1), %val(Status))
	      Sci_Rec.Collect_Time.Badtime_Flag = ten
	      Comperr = .True.
	      Compute = .False.
	   Else
	      Ref_Time(1) = End_Collect(1)
	      Ref_Time(2) = End_Collect(2)
	   EndIf
	EndIf

C  Extract science record's adds per group from science header.

	Adds = sci_rec.sci_head.sc_head9
	If (Adds .Eq. zero) Then
	   sci_rec.sci_head.sc_head9 = 1
	ElseIf (Adds.Ne.1 .And. Adds.Ne.2 .And. Adds.Ne.3 .And. Adds.Ne.8 .And.
	1  Adds.Ne.12 .And. Sci_Rec.Collect_Time.Badtime_Flag.EQ.0) Then
	   Sci_Rec.Collect_Time.Badtime_Flag = seventeen
	EndIf

C  Set the "Correlated Flag" to True or False.

	If ((FPP_Collect_Time .EQ. %loc(FPP_Normal)) .And. (.Not. Comperr)) Then
	  Correlated = .False.
	Else
	  Correlated = .True.
	Endif

C  Determine if the reference frame time is before the first MF of the NFS_ANC,
C  within the first MF, within the second MF or after the NFS_ANC record.
C  Process accordingly.  If telemetry format is bad for frame, its time is
C  not correct :  Care must be taken if format is bad (non-science).  SEE
C  FPP Design Walkthrough document for diagrams of the following cases and
C  flow chart.  Remember TF = true means Good (science) telemetry format.

	Do While ( .Not. Correlated )

	   If (TF_R2(1)) Then
	      If (Time_LT (Ref_Time, Anc_Rec(2).CT_Head.Time)) Then
	         If (TF_R1(2)) Then
	            Status = Lib$Subx (Anc_Rec(2).CT_Head.Time,
	1              Anc_Rec(1).GMT_MJF2, Delta_Time, Len)
	            If (.Not. Status) Then
	               Call Lib$Signal(FPP_SubxErr, %val(1), %val(Status))
	               Sci_Rec.Collect_Time.Badtime_Flag = ten
	               Comperr = .True.
	               Compute = .False.
	            Else
	               If (Time_GT (Delta_Time, QMax_Sec)) Then

C  Case 1. Ref_Time is in gap before first frame.

	                  Sci_Rec.Collect_Time.Badtime_Flag = nine
	                  Compute = .False.
	               Else

C  Case 3. Ref_Time is between the last minor frame of the first ANC record and
C  the start time of the second ANC record. Thus the .25 secs per minor frame is
C  not valid. The last index to the status bits array, 128, is used to get the
C  MTM mode.

	                  Index = 128
	                  I4_Bits =
	1                    Anc_Rec(1).Frame(2).Minor_Frame_Status_Bits(Index)
	                  Sci_Rec.Sci_Head.MTM_Speed = Jibits(I4_Bits,4,1)
	                  Sci_Rec.Sci_Head.MTM_Length = Jibits(I4_Bits,5,1)
	               EndIf				! If Gap
	            EndIf				! If status of SubX
	         ElseIf (Anc_Counter .LE. 1) Then

C Same as Case 1. Ref_Time is in gap before first frame.

	            Sci_Rec.Collect_Time.Badtime_Flag = one
	            Compute = .False.
	            If (report) Write (Lun_rpt,30)
  30	            Format (5X,'Science collection reference time occurred ',
	1             'before 1st ancillary record.')
	         Else
	            
C  Case 4.  Ref_Time points to frame with non-science telemetry format.

	            Sci_Rec.Collect_Time.Badtime_Flag = sixteen
	            Compute = .False.
	         EndIf					! Previous record TQ2
	         Correlated = .True.

	      Else					! Ref time > MF1 time

	         If (Time_LT (Ref_Time, Last_MNF_Time1)) Then

C  Cases 5 and 6. Ref_Time points in frame 1. Compute the index to the NFS_ANC
C  Word 31 minor frame status bit array.  Set MNF index to
C  (Ref_Time - MF_1_Start_Time) * 128 frames / 32.0. Since the MNF goes 0 to
C  127 and the record array goes 1 to 128, add 1.

	            Status =
	1              Lib$Subx (Ref_Time, Anc_Rec(2).CT_Head.Time, Del_Ref,Len)
	            If (.Not. Status) Then
	               Call Lib$Signal (FPP_SubxErr, %val(1), %val(Status))
	               Comperr = .True.
	               Compute = .False.
	            Else
	               Status = Lib$Ediv (Sec_per_MNF, Del_Ref, Index,Remainder)
	               If (.Not. Status) Then
	                  Call Lib$Signal(FPP_EdivErr, %val(1), %val(Status))
	                  Comperr = .True.
	                  Compute = .False.
	               Else
	                  If (Index .Lt. 128) Index = Index + 1
	                  I4_Bits =
	1                    Anc_Rec(2).Frame(1).Minor_Frame_Status_Bits(Index)
	                  Sci_Rec.Sci_Head.MTM_Speed = Jibits(I4_Bits,4,1)
	                  Sci_Rec.Sci_Head.MTM_Length = Jibits(I4_Bits,5,1)
	               EndIf
	            EndIf
	            Correlated = .True.
	         Else
	            If (.Not.TF_R2(2)) Then

C  Case 8.  Ref_Time is after frame 1 and frame 2 is bad. Get next 2 frames
C  where at least one frame is good Telemetry format.

	               status = FPP_Read2 ( Lun_Anc, Anc_Rec, TF_R2, TF_R1,
	1                 Gap_Check, Anc_Gap, Report, Lun_Rpt, Last_MNF_Time1,
	2                 Last_MNF_Time2, First_Anc )
	               If (status .NE. %loc(FPP_Normal)) Then
	                  FPP_Collect_Time = status
	                  Correlated = .True.
	                  Compute = .False.
	               EndIf
	            Else
	               If (Time_LT (Ref_Time, Anc_Rec(2).GMT_MJF2)) Then
	                  Status = Lib$Subx ( Anc_Rec(2).GMT_MJF2,
	1                    Anc_Rec(2).CT_Head.Time, Delta_Time, Len)
	                  If (.Not. Status) Then
	                     Call Lib$Signal(FPP_SubxErr, %val(1), %val(Status))
	                     Sci_Rec.Collect_Time.Badtime_Flag = ten
	                     Comperr = .True.
	                     Compute = .False.
	                  Else
	                     If (Time_GT (Delta_Time, QMax_Sec)) Then

C  Case 9.  Ref_Time points to a gap between major frames 1 and 2.

	                        Sci_Rec.Collect_Time.Badtime_Flag = nine
	                        Compute = .False.
	                        Correlated = .True.
	                     Else

C  Case 12.  Ref_Time is between the last minor frame of major frame 1 and the
C  the start time of major frame 2 in the second ANC record. Thus the .25 secs
C  per minor frame is not valid. The last index to the status bits array, 128,
C  is used to get the MTM mode.

	                        Index = 128
	                        I4_Bits =
	1                     Anc_Rec(2).Frame(1).Minor_Frame_Status_Bits(Index)
	                        Sci_Rec.Sci_Head.MTM_Speed = Jibits(I4_Bits,4,1)
	                        Sci_Rec.Sci_Head.MTM_Length =Jibits(I4_Bits,5,1)
	                        Correlated = .True.
	                     EndIf			! If Gap
	                  EndIf				! If status of SubX
	               Else				! Ref > MF2
	                  If (Time_LT (Ref_Time, Last_MNF_Time2)) Then

C  Case 10.  Ref_Time points in frame 2. Compute the index to the NFS_ANC
C  Word 31 minor frame status bit array.  Set MNF index to
C  (Ref_Time - MF_2_Start_Time) * 128 frames / 32.0. Since the MNF goes 0 to
C  127 and the record array goes 1 to 128, add 1.

	                     Status = Lib$Subx (Ref_Time,
	1                       Anc_Rec(2).GMT_MJF2, Del_Ref,Len)
	                     If (.Not. Status) Then
	                       Call Lib$Signal(FPP_SubxErr,%val(1),%val(Status))
	                       Comperr = .True.
	                       Compute = .False.
	                     Else
	                        Status = Lib$Ediv (Sec_per_MNF, Del_Ref, Index,
	1                          Remainder)
	                        If (.Not. Status) Then
	                           Call Lib$Signal (FPP_EdivErr, %val(1),
	1                             %val(Status))
	                           Comperr = .True.
	                           Compute = .False.
	                        Else
	                           If (Index .Lt. 128) Index = Index + 1
	                           I4_Bits =
	1                     Anc_Rec(2).Frame(2).Minor_Frame_Status_Bits(Index)
	                          Sci_Rec.Sci_Head.MTM_Speed=Jibits(I4_Bits,4,1)
	                         Sci_Rec.Sci_Head.MTM_Length=Jibits(I4_Bits,5,1)
	                        EndIf
	                     EndIf			! SubX status
	                     Correlated = .True.
	                  Else				! Ref > End of MJF2

C  Case 11.  Ref_Time points after end of Major frame 2.  Get next record that
C  contains at least one with science format.

	                     status = FPP_Read2 ( Lun_Anc, Anc_Rec, TF_R2,
	1                       TF_R1, Gap_Check, Anc_Gap, Report, Lun_Rpt,
	2                       Last_MNF_Time1, Last_MNF_Time2, First_Anc )
	                     If (status .NE. %loc(FPP_Normal)) Then
	                        FPP_Collect_Time = status
	                        Correlated = .True.
	                        Compute = .False.
	                     EndIf
	                  EndIf		! Ref < End of MJF2
	               EndIf		! Ref < MJF2
	            EndIf		! TF of Rec 2 Frm 2
	         EndIf			! Ref < End of MJF1
	      EndIf			! Ref < MJF1
	   Else				! Telemetry format Rec 2, Frm 1
	      If (Time_LT (Ref_Time, Anc_Rec(2).GMT_MJF2)) Then

C  Cases 2 and 7.  Ref_Time points prior to Frame 2 and frame 1 has non-science
C  telemetry format.

	         Compute = .False.
	         Correlated = .True.
	         If (Anc_Counter .LE. 1) Then

C  Ref_Time is in gap before first frame.

	            Sci_Rec.Collect_Time.Badtime_Flag = one
	            If (report) Write (Lun_rpt,30)
	         Else

	            Sci_Rec.Collect_Time.Badtime_Flag = sixteen

	         EndIf
	      Else			! Ref > MJF2
	         If (Time_LT (Ref_Time, Last_MNF_Time2)) Then

C  Case 10 (again).  Ref_Time points in frame 2. Compute the index to the
C  NFS_ANC Word 31 minor frame status bit array.  Set MNF index to
C  (Ref_Time - MF_2_Start_Time) * 128 frames / 32.0. Since the MNF goes 0 to
C  127 and the record array goes 1 to 128, add 1.

	            Status = Lib$Subx (Ref_Time,Anc_Rec(2).GMT_MJF2,Del_Ref,Len)
	            If (.Not. Status) Then
	               Call Lib$Signal(FPP_SubxErr,%val(1),%val(Status))
	               Comperr = .True.
	               Compute = .False.
	            Else
	               Status = Lib$Ediv (Sec_per_MNF, Del_Ref, Index,Remainder)
	               If (.Not. Status) Then
	                  Call Lib$Signal (FPP_EdivErr, %val(1), %val(Status))
	                  Comperr = .True.
	                  Compute = .False.
	               Else
	                  If (Index .Lt. 128) Index = Index + 1
	                  I4_Bits =
	1                    Anc_Rec(2).Frame(2).Minor_Frame_Status_Bits(Index)
	                  Sci_Rec.Sci_Head.MTM_Speed = Jibits(I4_Bits,4,1)
	                  Sci_Rec.Sci_Head.MTM_Length = Jibits(I4_Bits,5,1)
	               EndIf
	            EndIf			! SubX status
	            Correlated = .True.
	         Else				! Ref > End of MJF2

C  Case 11 (again).  Ref_Time points after end of Major frame 2.  Get next
C  record that contains at least one with science format.

	            status = FPP_Read2 ( Lun_Anc, Anc_Rec, TF_R2, TF_R1,
	1              Gap_Check, Anc_Gap, Report, Lun_Rpt, Last_MNF_Time1,
	2              Last_MNF_Time2, First_Anc )
	            If (status .NE. %loc(FPP_Normal)) Then
	               FPP_Collect_Time = status
	               Correlated = .True.
	               Compute = .False.
	            EndIf
	         EndIf			! Ref < End of MJF2
	      EndIf			! Ref < MJF2
	   EndIf			! Telemetry format Rec 2, Frm 1
	EndDo				! While Not Correlated

C  If the "Compute Flag" is set and status is good, then Call FPP_Midpoint_Time
C  to compute midpoint of collect time using MTM_Scan_Speed and MTM_Scan_Length.

	If (( Compute ) .And. (FPP_Collect_Time .Eq. %loc(FPP_Normal))) Then
	   If (Sweeps .EQ. 0) Sweeps = 1
	   Status = FPP_Midpoint_Time ( Sci_Rec.Sci_Head.MTM_Length,
	1     Sci_Rec.Sci_Head.MTM_Speed, Sweeps, Sci_Rec.CT_Head.Time, Sync,
	2     Start_Collect, End_Collect, IFG_Time, AncFilename )

	   If ( Status .Eq. One ) Then
	      Sci_Rec.Collect_Time.Midpoint_Time(1) = IFG_Time(1)
	      Sci_Rec.Collect_Time.Midpoint_Time(2) = IFG_Time(2)

C  Call FPP_Check_Time for the consistency check of the MTM scan speed
C  and length on the selected number of previous minor frames.

	      Status = FPP_Check_Time ( Prev_Time, Sci_Rec, Anc_Rec, TF_R1,
	1        TF_R2, Start_Collect, End_Collect, QMax_Sec, Lun_Anc,
	2        Last_MNF_Time1, Last_MNF_Time2, Gap_Check, Anc_Gap, Report,
	3        Lun_Rpt, First_Anc )

	   Else
	      Comperr = .True.
	   EndIf	! Status return of FPP_Midpoint_Time
	Endif 	 	! The "Compute Flag" is set and the status is good

C  If there was an arithmetic error, flag the record.
	
	If (Comperr .And. Sci_Rec.Collect_Time.Badtime_Flag.EQ.0) Then
	   Sci_Rec.Collect_Time.Badtime_Flag = Ten
	Endif

C  Save the minor frame transmit time from the current science record for
C  comparison with the next science record.

	Prev_Time(1) = Sci_Rec.CT_Head.Time(1)
	Prev_Time(2) = Sci_Rec.CT_Head.Time(2)

	Return
	End
