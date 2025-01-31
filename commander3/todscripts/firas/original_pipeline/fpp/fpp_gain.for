C*******************************************************************************
	Integer*4  Function  FPP_GAIN (Chan_Num, Sci_Rec, First_Rec, Lun_Gain,
	1		Report, Lun_Rpt, Sync, Start_Collect, End_Collect)
C*******************************************************************************
C	Author: Larry P. Rosen, STX, March 1991

C	Purpose: Extract gain from Fex_Gain reference data file or flag
C	   housekeeping gap, telemetry quality, or unknown gain.

C	Input Parameters:
C	  Name           Type	          Description
C 	  ----------------------------------------------------------------------
C	  Chan_Num        I*2		Channel number
C	  Sci_Rec   NFS_SDF Rec		Current FIRAS Raw Science Record.
C	  First_Rec       L*1		Flag 1st science record
C	  Lun_Gain        I*4		Logical unit for gain reference file.
C	  Report          L*1		Flag to write report.
C	  Lun_Rpt         I*4		Logical unit for report.
C	  Sync	          L*1		Flag MTM in or out of sync.
C	  Start_Collect   I*4(2)	Start time of science collection.
C	  End_Collect     I*4(2)	End time of science collection.
C
C	Output Parameters:
C	  Name		Type		  Description
C 	  ----------------------------------------------------------------------
C	  Sci_Rec   NFS_SDF Rec		With gains or flags set.
C-------------------------------------------------------------------------------
C	(NEW)  PDL  for  FPP_GAIN	Larry P. Rosen, STX, 5 February 1991
C
C  NOTE:  Reference file was opened (JSTART - offset) to JSTOP.  Two buffers
C	will be used to hold new and old reference data records.
C
C	IF (FirstTime) Then
C	    READ FEX_GAIN.DAT record into "new buffer" (FEX_New).
C	    If (Collect time LT FEX_Gain start time) signal error
C	ENDIF
C	IF (BadTime flag NE 0) Then
C	    Set gain value = -1
C	    Return
C	ENDIF
C	IF ( (Start_Collect_Time .GE. FEX_New.Start_Time)  AND
C	     (End_Collect_Time .LE. FEX_New.Stop_Time) )
C	THEN (collection taken at this gain):
C
C	    Copy gain from FEX_New record to science record.
C
C	ELSEIF ( End_Collect_Time .LE. FEX_New.Stop_Time)
C	THEN (collection spans break in reference records):
C
C	    Check flags in the old buffer reference record to determine
C	    reason why record ended.
C
C	    IF (Reason .EQ. Telemetry Quality Bad) THEN
C		    Set housekeeping telemetry quality flag in science record.
C		    Set BadTime flag in science record.
C		    Assign gain "Unknown" (-1).
C
C	    ELSE IF (Reason .EQ. Housekeeping Data Gap) THEN
C		    Set housekeeping data gap flag in science record.
C		    Set BadTime flag in science record.
C		    Assign gain "Unknown" (-1).
C
C	    ELSE IF (Gain didn't change for THIS channel) THEN
C		    (Gain changed for a different channel)
C		    Copy gain from FEX_Old record to science record.
C
C	    ELSE (Gain changed for this channel)
C
C		IF (MTM "Out of Sync") THEN
C			Set BadTime flag in science record.
C			Assign gain "Unknown" (-1).
C	    	ELSE
C		    IF ( (FEX_New.Start_Time .GE. Start_Collect_Time)
C			AND  (FEX_New.Start_Time .LE. End_Collect_Time) )
C		    THEN
C			(New gain sample time during collection)
C			Copy new gain from FEX_New record to science_record.
C
C		    ELSE IF ( (FEX_Old.Stop_Time .GE. Start_Collect_Time)
C			AND  (FEX_Old.Stop_Time .LE. End_Collect_Time) )
C		    THEN
C			(Old gain sample time during collection)
C			Copy old gain from FEX_Old record to science_record.
C		    ELSE
C			Assign gain "Unknown" (-1).
C			Set BadTime flag in science record.
C		    ENDIF  (sample during collect)
C		ENDIF  (sync)
C	    ENDIF  (reason)
C
C	ELSE  (find reference records that bracket collection):
C
C	    Bracket = FALSE.
C	    DO WHILE (.NOT. Bracket)
C		Move FEX_New record into "old" buffer, FEX_Old.
C		Read a new FEX record into "new" buffer.
C		IF (End_Collect_Time .LE. FEX_New.Stop_Time) THEN Bracket= TRUE.
C	    END DO
C
C	    IF (Start_Collect_Time .GE. FEX_New.Start_Time)
C	    THEN (collection taken at this new record's gain):
C
C		Copy gain from FEX_New record to science record.
C
C	    ELSE  (collection spans break in reference records):
C
C		Check flags in the old buffer reference record to determine
C		reason why record ended.
C
C	 	IF (Reason .EQ. Telemetry Quality Bad) THEN
C		    Set housekeeping telemetry quality flag in science record.
C		    Set BadTime flag in science record.
C		    Assign gain "Unknown" (-1).
C
C		ELSE IF (Reason .EQ. Housekeeping Data Gap) THEN
C		    Set housekeeping data gap flag in science record.
C		    Set BadTime flag in science record.
C		    Assign gain "Unknown" (-1).
C
C		ELSE IF (Gain didn't change for THIS channel) THEN
C		    (Gain changed for a different channel)
C		    Copy gain from FEX_Old record to science record.
C
C		ELSE (Gain changed for this channel)
C
C		    IF (MTM "Out of Sync") THEN
C			Set BadTime flag in science record.
C			Assign gain "Unknown" (-1).
C		    ELSE
C			IF ( (FEX_New.Start_Time .GE. Start_Collect_Time)
C			    AND  (FEX_New.Start_Time .LE. End_Collect_Time) )
C			THEN
C				(New gain sample time during collection)
C			    Copy new gain from FEX_New record to science_record.
C
C			ELSE IF ( (FEX_Old.Stop_Time .GE. Start_Collect_Time)
C			    AND  (FEX_Old.Stop_Time .LE. End_Collect_Time) )
C			THEN
C				(Old gain sample time during collection)
C			    Copy old gain from FEX_Old record to science_record.
C
C			ELSE
C			    Assign gain "Unknown" (-1).
C			ENDIF  (sample during collect)
C		    ENDIF  (sync)
C		ENDIF  (reason)
C	    ENDIF  (collection in new reference record)
C	ENDIF  (collection in current reference record)
C
C	RETURN with normal or error status.
C*******************************************************************************
	Implicit None

	Integer*2	Chan_Num	! Channel number
	Logical*1	First_Rec	! First Science record call to routine
	Integer*4	Lun_Gain	! Logical unit for gain reference file.
	Logical*1	Report		! Flag to write report.
	Integer*4	Lun_Rpt		! Logical unit for report.
	Logical*1	Sync		! Flag MTM in or out of sync.
	Integer*4	Start_Collect(2)  ! Start time of science collection.
	Integer*4	End_Collect(2)    ! End time of science collection.
	Dictionary	'NFS_SDF'
	Record / NFS_SDF / Sci_Rec      ! Current FIRAS Raw Science Record

	Include		'(FPP_Msg)'
	Include		'CT$Library:CTUser.Inc'

C  Local variables

	Integer*2	Chan
	Dictionary 'fex_gain'
	Record / fex_gain / Gain_Rec(2)	! Two records from Gain reference file
	Integer*2	CT_Stat(20)     ! COBETRIEVE return status
	Integer*4	Status		! Return from a function, 0 = good
	Logical*1	Bracket		! Flag collection bracketed by 2 Fex rec
	Logical*1	TIME_LT
	Logical*1	TIME_LE
	Logical*1	TIME_GE
	Record / fex_gain / Blank_Rec	! Empty record for initializing

	FPP_Gain = %loc(FPP_Normal)
	Chan = 2 - mod (chan_num , 2)
	If (First_Rec) Then
	   Gain_Rec(1) = Blank_Rec
	   Call CT_Read_Arcv (, Lun_Gain, Gain_Rec(2), CT_Stat )
	   If (CT_Stat(1) .NE. CTP_Normal) Then
	      If (CT_Stat(1) .EQ. CTP_Endoffile) Then
	         Call Lib$Signal (FPP_CTRefEOF)
	      Else
	         Call Lib$Signal (FPP_CTReadErr, %val(1), %val(CT_Stat(1)),
	1            'Fex_Gain.Dat')
	         FPP_Gain = %loc(FPP_Aberr)
	      EndIf
	   Else					! CT_Stat(1) = CTP_Normal
	      If (Time_LT (Start_Collect, Gain_Rec(2).Start_Time) )  Then
	         Call Lib$Signal (FPP_LateRef)
	         FPP_Gain = %loc(FPP_Aberr)
	      EndIf
	   EndIf
	EndIf					! If first record
	If (Sci_Rec.Collect_Time.BadTime_Flag .NE. 0) Then
	   Sci_Rec.Sci_Head.Gain = -1
	   Return
	EndIf
	If (FPP_Gain .EQ. %loc(FPP_Normal)) Then
	   If (Time_GE (Start_Collect, Gain_Rec(2).Start_Time) .AND.
	1      Time_LE (End_Collect, Gain_Rec(2).Stop_Time) ) Then

C Collection taken at this gain. Copy from FEX_New record to science record.

	      Sci_Rec.Sci_Head.Gain = Gain_Rec(2).Gain(Chan)

	   ElseIf (Time_LE (End_Collect, Gain_Rec(2).Stop_Time) ) Then

C  Collection spans break in reference records.  Check flags in the old buffer
C  reference record to determine reason why record ended.

	      If (Gain_Rec(1).Tlm_Bad_Quality) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 15
	            Sci_Rec.Sci_Head.Gain = -1

	      ElseIf (Gain_Rec(1).Gain_Data_Gap) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 14
	            Sci_Rec.Sci_Head.Gain = -1

	      ElseIf (.NOT. Gain_Rec(1).Gain_Change(Chan) ) Then

C Gain didn't change for THIS channel. Copy gain from FEX_Old record.

	            Sci_Rec.Sci_Head.Gain = Gain_Rec(1).Gain(Chan)

	      Else		! Gain changed for this channel

	         If (.NOT. Sync) Then		! MTM is Out of Sync
	               Sci_Rec.Collect_Time.BadTime_Flag = 11
	               Sci_Rec.Sci_Head.Gain = -1

	         Else

	            If (Time_GE (Gain_Rec(2).Start_Time, Start_Collect).AND.
	2               Time_LE (Gain_Rec(2).Start_Time, End_Collect) ) Then

C New gain was sampled during collection.  Copy gain from FEX_New record.

	               Sci_Rec.Sci_Head.Gain = Gain_Rec(2).Gain(Chan)

	            ElseIf (Time_GE (Gain_Rec(1).Stop_Time, Start_Collect).AND.
	2                   Time_LE (Gain_Rec(1).Stop_Time, End_Collect) ) Then

C Old gain was sampled during collection.  Copy gain from FEX_Old record.

	               Sci_Rec.Sci_Head.Gain = Gain_Rec(1).Gain(Chan)

	            Else
	               Sci_Rec.Collect_Time.BadTime_Flag = 11
	               Sci_Rec.Sci_Head.Gain = -1
	            EndIf			! sample during collect
	         EndIf				! sync
	      EndIf				! reason for reference gap

	   Else		! Find reference records that bracket collection

	      Bracket = .False.
	      Do While (.NOT. Bracket)
	         Gain_Rec(1) = Gain_Rec(2)
	         Call CT_Read_Arcv (, Lun_Gain, Gain_Rec(2), CT_Stat )
	         If (CT_Stat(1) .NE. CTP_Normal) Then
	            If (CT_Stat(1) .EQ. CTP_Endoffile) Then
	               Call Lib$Signal (FPP_CTRefEOF)
	               Bracket = .True.
	            Else
	               Call Lib$Signal(FPP_CTReadErr, %val(1), %val(CT_Stat(1)),
	1            'Fex_Gain.Dat')
	               FPP_Gain = %loc(FPP_Aberr)
	            EndIf
	         Else				! CT_Stat(1) = CTP_Normal

C  IF (End_Collect_Time .LE. FEX_New.Stop_Time) THEN Bracket = TRUE.

	            If (Time_LE (End_Collect, Gain_Rec(2).Stop_Time) ) Then
	               Bracket = .True.
	            EndIf
	         EndIf
	      EndDo			! Until Bracket = true

	      If (Time_GE (Start_Collect, Gain_Rec(2).Start_Time) .AND.
	2         Time_LE (End_Collect, Gain_Rec(2).Stop_Time) ) Then

C Collection taken at this gain. Copy from FEX_New record to science record.

	         Sci_Rec.Sci_Head.Gain = Gain_Rec(2).Gain(Chan)

	      Else

C  Collection spans break in reference records.  Check flags in the old buffer
C  reference record to determine reason why record ended.

	         If (Gain_Rec(1).Tlm_Bad_Quality) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 15
	            Sci_Rec.Sci_Head.Gain = -1

	         ElseIf (Gain_Rec(1).Gain_Data_Gap) Then
	            Sci_Rec.DQ_Data.Data_Quality(2) = 1
	            Sci_Rec.Collect_Time.BadTime_Flag = 14
	            Sci_Rec.Sci_Head.Gain = -1

	         ElseIf (.NOT. Gain_Rec(1).Gain_Change(Chan) ) Then

C Gain didn't change for THIS channel. Copy gain from FEX_Old record.

	            Sci_Rec.Sci_Head.Gain = Gain_Rec(1).Gain(Chan)

	         Else		! Gain changed for this channel

	            If (.NOT. Sync) Then		! MTM is Out of Sync
	               Sci_Rec.Collect_Time.BadTime_Flag = 11
	               Sci_Rec.Sci_Head.Gain = -1

	            Else

	               If (Time_GE (Gain_Rec(2).Start_Time, Start_Collect) .AND.
	2                  Time_LE (Gain_Rec(2).Start_Time, End_Collect) ) Then

C New gain was sampled during collection.  Copy gain from FEX_New record.

	                  Sci_Rec.Sci_Head.Gain = Gain_Rec(2).Gain(Chan)

	               ElseIf (Time_GE (Gain_Rec(1).Stop_Time, Start_Collect)
	2               .AND. Time_LE(Gain_Rec(1).Stop_Time, End_Collect) ) Then

C Old gain was sampled during collection.  Copy gain from FEX_Old record.

	                  Sci_Rec.Sci_Head.Gain = Gain_Rec(1).Gain(Chan)

	               Else
	                  Sci_Rec.Collect_Time.BadTime_Flag = 11
	                  Sci_Rec.Sci_Head.Gain = -1
	               EndIf			! sample during collect
	            EndIf			! sync
	         EndIf				! reason for reference gap
	      EndIf			! collection in new reference record
	   EndIf			! collection in current reference record
	EndIf				! Status check

	Return
	End
