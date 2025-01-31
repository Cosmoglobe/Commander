	Integer*4 Function FPP_Read2 (Lun_Anc, Anc_Rec, TF_R2, TF_R1, Gap_Check,
	1   Anc_Gap, Report, Lun_Rpt, Last_MNF_Time1, Last_MNF_Time2, First_Anc)
C*******************************************************************************
C  Larry P. Rosen, STX, 3 April 1991
C  Called by FPP_COLLECT_TIME.
C  Reads the next two ancillary frames where at least one frame has science
C  format.  If Gap_Check = 1 and telemetry format is ok then check for anc gaps.
C*******************************************************************************
	Implicit None

	Include		'(FPP_Msg)'
	Include		'(FUT_Params)'
	Include		'($SSDef)'
	Include		'CT$Library:CTUser.Inc'

C Passed.
	Integer*4	Lun_Anc			! Ancillary logical unit
	Dictionary	'NFS_ANC'
	Record / NFS_ANC / Anc_Rec(2)		! Current FIRAS Word 31 Records
	Logical*1	TF_R1(2), TF_R2(2)	! Telem. Format each record(frm)
	Integer*2	Gap_Check		! Flag to check for anc gaps.
	Integer*4	Anc_Gap			! Gap for Anc frames/ reporting
	Logical*1	Report			! Flag to write report
	Integer*4	Lun_Rpt			! Report logical unit
	Integer*4	Last_MNF_Time1(2)	! Quadword time end of MJF 1
	Integer*4	Last_MNF_Time2(2)	! Quadword time end of MJF 2
	Logical*1	First_Anc		! Flag 1st Anc record of day

C Functions.
	Integer*4	Lib$Emul
	Integer*4	Lib$Subx
	Integer*4	Lib$Addx
	Logical*1	Time_GT
	Integer*4	FUT_GMT_Conv

C Local.
	Logical*1	First_Time /.True./
	Integer*4	Status
	Integer*4	One /1/
	Integer*4	Last_MNF_Sec    ! Time increment for last minor frame
	Parameter	( Last_MNF_Sec = 317500000 )
	Integer*4	One_Sec
	Parameter	( One_Sec = 10000000 )
	Integer*4	QAnc_Gap(2)	! Quadword anc gap size.
	Integer*4	Addend /0/	! Additive for Emul
	Integer*4	Del_Last(2)	! Time to end of frame in quadword
	Record / NFS_ANC / Blank_rec	! Empty ANC record for clearing buffer
	Logical*1	Bad_Ancil	! Flag to keep reading
	Integer*2	CT_Stat(20)
        Character*14	GMT2
	Integer*4	Delta_Time(2)	! Difference of times between frames
	Integer*4	Len /2/		! Dimension for Lib$SubX
        Real*8		time_sec1, time_sec2, tjump
        Integer*2	ihr, imin
	Integer*4	msec
	Integer*4	Count /0/	! Count number of records at start

	FPP_Read2 = %loc(FPP_Normal)
	If (First_Time) Then
	   Status = Lib$Emul (One, Last_MNF_Sec, Addend, Del_Last)
	   If (.Not. Status) Then
	      Call Lib$Signal(FPP_EmulErr, %val(1), %val(Status))
	      FPP_Read2 = %loc(FPP_Aberr)
	   EndIf
	   Status = Lib$Emul (One_Sec, Anc_Gap, Addend, QAnc_Gap)
	   If (.Not. Status) Then
	      Call Lib$Signal(FPP_EmulErr, %val(1), %val(Status))
	      FPP_Read2 = %loc(FPP_Aberr)
	   EndIf
	   TF_R2(1) = .False.
	   TF_R2(2) = .False.
	   First_Time = .False.
	EndIf

C  Begin reading until at least 1 frame has science format.

	If (First_Anc) Then
	   Count = 0
	   Anc_Rec(2) = Blank_rec
	EndIf
	Bad_Ancil = .True.
	Do While (Bad_Ancil)
	   Anc_Rec(1) = Anc_Rec(2)
	   TF_R1(1) = TF_R2(1)
	   TF_R1(2) = TF_R2(2)
	   Call CT_Read_Arcv (, Lun_Anc, Anc_Rec(2), CT_Stat )
	   If (CT_Stat(1) .Ne. CTP_Normal) Then
	      If (CT_Stat(1) .Eq. CTP_Endoffile) Then
	         Call Lib$Signal(FPP_CTAncEof)
	         FPP_Read2 = %loc(FPP_CTAncEof)
	      Else
	         Call Lib$Signal(FPP_CTReadErr, %val(2), %val(CT_Stat(1)),
	1           'NFS_ANC')
	         FPP_Read2 = %loc(FPP_Aberr)
	      EndIf
	      Return
	   Else				! CT_Stat(1) = CTP_Normal

	      Count = Count + 1

C  Set Telemetry format flags.  1 means science.
C  In frame by frame check, only need 1 good frame initially.

	      If (Anc_Rec(2).CT_Head.Hskp1_Tlm_Fmt .EQ. 1) Then
	         TF_R2(1) = .True.
	      Else
	         TF_R2(1) = .False.
	      EndIf
	      If (Anc_Rec(2).CT_Head.Hskp2_Tlm_Fmt .EQ. 1) Then
	         TF_R2(2) = .True.
	      Else
	         TF_R2(2) = .False.
	      EndIf
	      If (TF_R2(1) .OR. TF_R2(2)) Bad_Ancil = .False.
	   EndIf				! read status
	EndDo					! Read loop

C  Check for gap in ancillary data for reporting.  If first anc record then
C  count = 1 and so don't check before frame 1.

	If (Gap_Check .EQ. 1 .AND. Report) Then
	   If (Count .GT. 1) Then
	      Call CT_binary_to_GMT (Anc_Rec(1).Gmt_Mjf2, GMT2)
	      If (TF_R2(1).AND.TF_R1(2)) Then
	         Status = Lib$SubX (Anc_Rec(2).CT_Head.Time,Anc_Rec(1).GMT_Mjf2,
	1           Delta_Time, Len)
	         If (.Not. Status) Then
	            Call Lib$Signal (FPP_SubxErr, %val(1), %val(Status))
	            FPP_Read2 = %loc(FPP_AbErr)
	            Return
	         EndIf
	         If (Time_GT (Delta_Time, QAnc_Gap)) Then
	            status = FUT_GMT_Conv (time_sec1,Anc_Rec(2).CT_Head.GMT,2)
	            status = FUT_GMT_Conv (time_sec2,GMT2,2)
	            tjump = time_sec1 - time_sec2
	            ihr = Int (tjump / fac_hour)
	            tjump = tjump - ihr * fac_hour
	            imin = Int (tjump / fac_minute)
	            tjump = tjump - imin * fac_minute
	            msec = Nint (tjump)
	            Write (lun_rpt,20) gmt2, anc_rec(2).ct_head.gmt, ihr, imin,
	1              msec
  20	            Format (5X,'Gap: NFS_ANC', 8X, A14,' To ', A14,' = ',
	1              I2.2, ' Hr ', I2.2, ' Min ', I2.2,' Sec')
	         EndIf				! gap
	      EndIf				! TF of Rec 2 Frm 1
	   EndIf				! Count > 1
	   Call CT_binary_to_GMT (Anc_Rec(2).Gmt_Mjf2, GMT2)
	   If (TF_R2(2).AND.TF_R2(1)) Then
	      Status = Lib$SubX (Anc_Rec(2).GMT_Mjf2, Anc_Rec(2).CT_Head.Time,
	1        Delta_Time, Len)
	      If (.Not. Status) Then
	         Call Lib$Signal (FPP_SubxErr, %val(1), %val(Status))
	         FPP_Read2 = %loc(FPP_AbErr)
	         Return
	      EndIf
	      If (Time_GT (Delta_Time, QAnc_Gap)) Then
	         status = FUT_GMT_Conv (time_sec1,Anc_Rec(2).CT_Head.GMT,2)
	         status = FUT_GMT_Conv (time_sec2,GMT2,2)
	         tjump = time_sec2 - time_sec1
	         ihr = Int (tjump / fac_hour)
	         tjump = tjump - ihr * fac_hour
	         imin = Int (tjump / fac_minute)
	         tjump = tjump - imin * fac_minute
	         msec = Nint (tjump)
	         Write (lun_rpt,20) anc_rec(2).ct_head.gmt, gmt2, ihr, imin,
	1           msec
	      EndIf				! gap
	   EndIf				! TF of Rec 2 Frm 2
	EndIf					! Gap_Check = 1 and Report

C  Compute the time of the last minor frame of each Major frame in the current
C  NFS_ANC record.  Note that if telemetry format is non-science that the time
C  calculated is incorrect since the stripper put in the last science format
C  times.

	Status = Lib$Addx (Anc_Rec(2).CT_Head.Time, Del_Last,Last_MNF_Time1,Len)
	If (.Not. Status) Then
	   Call Lib$Signal(FPP_AddxErr, %val(1), %val(Status))
	   FPP_Read2 = %loc(FPP_AbErr)
	Else
	   Status = Lib$Addx ( Anc_Rec(2).GMT_MJF2, Del_Last,Last_MNF_Time2,Len)
	   If (.Not. Status) Then
	      Call Lib$Signal(FPP_AddxErr, %val(1), %val(Status))
	      FPP_Read2 = %loc(FPP_AbErr)
	   EndIf
	EndIf
	Return
	End
