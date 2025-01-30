C*******************************************************************************
	Integer*4  Function  FPP_FAKEIT (Chan_Num, Sci_Rec, First_Rec, Lun_Fake,
	1		Report, Lun_Rpt, Start_Collect, End_Collect)
C*******************************************************************************
C	Author: Larry P. Rosen, STX, March 1991

C	Purpose: Extract fakeit status from Fex_Fakeit reference data file or
C	   flag housekeeping gap, telemetry quality, or unknown status.

C	Input Parameters:
C	  Name           Type	          Description
C 	  ----------------------------------------------------------------------
C	  Chan_Num        I*2		Channel number
C	  Sci_Rec   NFS_SDF Rec		Current FIRAS Raw Science Record.
C	  First_Rec       L*1		Flag 1st science record
C	  Lun_Fake        I*4		Logical unit for fakeit reference file.
C	  Report          L*1		Flag to write report.
C	  Lun_Rpt         I*4		Logical unit for report.
C	  Start_Collect   I*4(2)	Start time of science collection.
C	  End_Collect     I*4(2)	End time of science collection.
C
C	Output Parameters:
C	  Name		Type		  Description
C 	  ----------------------------------------------------------------------
C	  Sci_Rec   NFS_SDF Rec		With fakeit status or flags set.
C-------------------------------------------------------------------------------
C	(NEW)  PDL  for  FPP_FAKEIT	Larry P. Rosen, STX, 5 February 1991
C
C  NOTE:  Reference file was opened (JSTART - offset) to JSTOP.  Two buffers
C	will be used to hold new and old reference data records.
C
C	IF (FirstTime) Then
C	    READ FEX_FAKEIT.DAT record into "new buffer" (FEX_New).
C	    If (Collect time LT FEX_FAKEIT start time) signal error
C	ENDIF
C	IF (BadTime flag NE 0) Then
C	    Set fakeit value = -1 in science record
C	    Return
C	ENDIF
C	IF ( (Start_Collect_Time .GE. FEX_New.Start_Time)  AND
C	     (End_Collect_Time .LE. FEX_New.Stop_time) )
C	THEN (collection taken with this fakeit status):
C
C	    Copy fakeit status from FEX_New record to science record.
C
C	ELSEIF (End_Collect_Time .LE. FEX_New.Stop_time)
C	THEN (Collection started before current reference record):
C
C	    Check flags in the old buffer reference record to determine
C	    reason why record ended.
C
C	     IF (Reason .EQ. Telemetry Quality Bad) THEN
C		    Set housekeeping telemetry quality flag in science record.
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C	    ELSE IF (Reason .EQ. Housekeeping Data Gap) THEN
C		    Set housekeeping data gap flag in science record.
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C	    ELSE IF (Fakeit status changed for THIS channel) THEN
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C	    ELSE (Fakeit changed for a different channel)
C
C		    Copy fakeit status from FEX_Old  record to science_record.
C
C	    ENDIF  (reason)
C
C	ELSE  (find reference records that bracket collection):
C
C	    Bracket = FALSE.
C	    DO WHILE (.NOT. Bracket)
C		Move FEX_New record into "old" buffer, FEX_Old.
C		Read a new FEX record into "new" buffer.
C		IF (End_Collect_Time .LE. FEX_New.Stop_time) THEN Bracket= TRUE.
C	    END DO
C
C	    IF (Start_Collect_Time .GE. FEX_New.Start_Time)
C	    THEN (collection taken with this new record's fakeit status):
C
C		Copy fakeit status from FEX_New record to science record.
C
C	    ELSE  (collection spans break in reference records):
C
C		Check flags in the old buffer reference record to determine
C		reason why record ended.
C
C	 	IF (Reason .EQ. Telemetry Quality Bad) THEN
C		    Set housekeeping telemetry quality flag in science record.
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C		ELSE IF (Reason .EQ. Housekeeping Data Gap) THEN
C		    Set housekeeping data gap flag in science record.
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C		ELSE IF (Fakeit status changed for THIS channel) THEN
C		    Set BadTime flag in science record.
C		    Assign fakeit status "Unknown" (-1).
C
C		ELSE (Fakeit changed for a different channel)
C
C		    Copy fakeit status from FEX_Old  record to science_record.
C
C		ENDIF  (reason)
C	    ENDIF  (collection in new reference record)
C	ENDIF  (collection in current reference record)
C
C	RETURN with normal or error status.
C*******************************************************************************
	Implicit None

	Integer*2	Chan_Num	! Channel number
	Logical*1	First_Rec	! First Science record call to routine
	Integer*4	Lun_Fake	! Logical unit for fakeit refrence file.
	Logical*1	Report		! Flag to write report.
	Integer*4	Lun_Rpt		! Logical unit for report.
	Integer*4	Start_Collect(2)  ! Start time of science collection.
	Integer*4	End_Collect(2)    ! End time of science collection.
	Dictionary	'NFS_SDF'
	Record / NFS_SDF / Sci_Rec      ! Current FIRAS Raw Science Record

	Include		'(FPP_Msg)'
	Include		'CT$Library:CTUser.Inc'

C  Local variables

	Integer*2	Chan
	Dictionary 'fex_fakeit'
	Record / fex_fakeit / Fake_Rec(2)  ! Two records from Fakeit refrnc file
	Integer*2	CT_Stat(20)     ! COBETRIEVE return status
	Integer*4	Status		! Return from a function, 0 = good
	Logical*1	Bracket		! Flag collection bracketed by 2 Fex rec
	Logical*1	Time_LT
	Logical*1	Time_LE
	Logical*1	Time_GE
	Record / fex_fakeit / Blank_Rec  ! Empty record for initializing.

	FPP_Fakeit = %loc(FPP_Normal)
	Chan = 2 - mod ( chan_num, 2 )
	If (First_Rec) Then
	   Fake_Rec(1) = Blank_Rec
	   Call CT_Read_Arcv (, Lun_Fake, Fake_Rec(2), CT_Stat )
	   If (CT_Stat(1) .NE. CTP_Normal) Then
	      If (CT_Stat(1) .EQ. CTP_Endoffile) Then
	         Call Lib$Signal (FPP_CTRefEOF)
	      Else
	         Call Lib$Signal (FPP_CTReadErr, %val(1), %val(CT_Stat(1)),
	1           'Fex_FakeIt.DAt')
	         FPP_FakeIt = %loc(FPP_Aberr)
	      EndIf
	   Else					! CT_Stat(1) = CTP_Normal
	      If (Time_LT (Start_Collect, Fake_Rec(2).Start_Time) ) Then
	         Call Lib$Signal (FPP_LateRef)
	         FPP_FakeIt = %loc(FPP_Aberr)
	      EndIf
	   EndIf
	EndIf					! If first record
	If (Sci_Rec.Collect_Time.BadTime_Flag .NE. 0) Then
	   Sci_Rec.DQ_Data.Fake = -1
	   Return
	EndIf
	If (FPP_FakeIt .EQ. %loc(FPP_Normal)) Then
	   If (Time_GE (Start_Collect, Fake_Rec(2).Start_Time) .AND.
	2      Time_LE (End_Collect, Fake_Rec(2).Stop_Time) )  Then

C Collection taken with this fakeit status. Copy from FEX_New record.

	      Sci_Rec.DQ_Data.Fake = Fake_Rec(2).Fakeit(Chan)

	      ElseIf (Time_LE (End_Collect, Fake_Rec(2).Stop_Time) )  Then

C  Collection spans break in reference records.  Check flags in the old buffer
C  reference record to determine reason why record ended.

	         If (Fake_Rec(1).Tlm_Bad_Quality) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 15
	            Sci_Rec.DQ_Data.Fake = -1

	         ElseIf (Fake_Rec(1).FakeIt_Data_Gap) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 14
	            Sci_Rec.DQ_Data.Fake = -1

	         ElseIf (Fake_Rec(1).Fakeit_Change(Chan) ) Then

C Fakeit status changed for THIS channel.

	            Sci_Rec.Collect_Time.BadTime_Flag = 12
	            Sci_Rec.DQ_Data.Fake = -1

	         Else		! Fakeit status changed for a different channel

	            Sci_Rec.DQ_Data.Fake = Fake_Rec(1).Fakeit(Chan)

	         EndIf				! reason for reference gap

	   Else		! Find reference records that bracket collection

	      Bracket = .False.
	      Do While (.NOT. Bracket)
	         Fake_Rec(1) = Fake_Rec(2)
	         Call CT_Read_Arcv (, Lun_Fake, Fake_Rec(2), CT_Stat )
	         If (CT_Stat(1) .NE. CTP_Normal) Then
	            If (CT_Stat(1) .EQ. CTP_Endoffile) Then
	               Call Lib$Signal (FPP_CTRefEOF)
	               Bracket = .True.
	            Else
	               Call Lib$Signal(FPP_CTReadErr, %val(1), %val(CT_Stat(1)),
	1                 'Fex_FakeIt.DAt')
	               FPP_FakeIt = %loc(FPP_Aberr)
	            EndIf
	         Else				! CT_Stat(1) = CTP_Normal

C  IF (End_Collect_Time .LE. FEX_New.Stop_time) THEN Bracket = TRUE.

	            If (Time_LE (End_Collect, Fake_Rec(2).Stop_Time) ) Then
	               Bracket = .True.
	            EndIf
	         EndIf
	      EndDo			! Until Bracket = true

	      If (Time_GE (Start_Collect, Fake_Rec(2).Start_Time) .AND.
	2         Time_LE (End_Collect, Fake_Rec(2).Stop_Time) ) Then

C Collection taken with this fakeit status. Copy from FEX_New record.

	         Sci_Rec.DQ_Data.Fake = Fake_Rec(2).Fakeit(Chan)

	      Else

C  Collection spans break in reference records.  Check flags in the old buffer
C  reference record to determine reason why record ended.

	         If (Fake_Rec(1).Tlm_Bad_Quality) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 15
	            Sci_Rec.DQ_Data.Fake = -1

	         ElseIf (Fake_Rec(1).FakeIt_Data_Gap) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 14
	            Sci_Rec.DQ_Data.Fake = -1

	         ElseIf (Fake_Rec(1).Fakeit_Change(Chan) ) Then

C Fakeit status changed for THIS channel.

	            Sci_Rec.Collect_Time.BadTime_Flag = 12
	            Sci_Rec.DQ_Data.Fake = -1

	         Else		! Fakeit status changed for a different channel

	            Sci_Rec.DQ_Data.Fake = Fake_Rec(1).Fakeit(Chan)

	         EndIf				! reason for reference gap
	      EndIf			! collection in new reference record
	   EndIf			! collection in current reference record
	EndIf				! Status check

	Return
	End
