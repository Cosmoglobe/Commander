C*******************************************************************************
	Integer*4 Function FRD_GFG_Read2 ( Lun_hkp, Hkp_Rec, More_Data,
	1     Time_sec_1, Time_sec_2, Max_Sec, Frm1_Rec, Frm, Frm_Q,
	2     Lun_Rpt, Report, Act_MJ_A, Act_MJ_B )
C*******************************************************************************
C  Larry Rosen, STX, March 18 1991
C  If Frm = 2, then the examined frame was #2.  If Frm_Q = False (bad)
C  then frame 2 had bad quality and so we should begin reading housekeeping
C  records until 2 good consecutive frames are found.  They may be both in one
C  record (Frm1_Rec=1), or they may be in two different records (Frm1_Rec=2).
C  If Frm = 1, then the examined frame was #1.  If Frm_Q = False then
C  we should start with checking frame 2 and then read records.
C  If Frm_Q = true (good) then read next frame.  If it is good and no gap
C  then return that good record.  Else if it is good but gap, then get next
C  frame and check.  Else if bad, Then read records until 2 good frames are
C  found.
C  Determine actual major frame numbers from status monitor.  This information
C  is needed to extract gains, fakeit statuses and getting times.
C  The first time into this routine, Frm=2 and Frm_Q = False.  Check to see
C  that actual major frame A(1) .NE. A(2) and B(1) .NE. B(2) since if they are
C  then there is only data for one side of instrument.
C     PDL really wouldn't help anyone figure this code out (even me).  It is
C  much more useful to work through the flow chart (in desigh walkthrough),
C  while keeping in mind the above comments about frame quality Frm_Q, frame
C  number Frm, and telemetry quality TQ1 or TQ2.
C*******************************************************************************
	Implicit None

C Include.
      	Include 'CT$Library:CTUser.Inc'
C Functions.
	Integer*4       Lib$Signal  
	Integer*4	FUT_GMT_CONV
C External.
	External        FRD_Normal
	External        FRD_CTHkpEof
	External        FRD_Aberr
C Passed.
	Integer*4	LUN_HKP		! Unit number for NFS_HKP CT archive
        Dictionary 'NFS_HKP'
        Record /NFS_HKP/    HKP_REC(2)
	Logical*1       More_data	! Flag for more data
	Real*8		time_sec_2, time_sec_1
	Integer*4	Max_sec		! Max number secs before gap declared
	Integer*2	Frm1_Rec	! record number for frame 1
	Integer*2	Frm		! Current HKP frame examined (1 or 2)
	Logical*1	Frm_Q		! True = good, false = bad.
	Integer*4	LUN_Rpt		! Unit number for report file.
	Logical*1	Report		! Flag tells whether to write report.
        Integer*4       Act_Mj_A(2)     ! Actual a side major frame
        Integer*4       Act_mj_B(2)     ! Actual b side major frame
C Local.
 	Integer*2       TQ1, TQ2	! Telemetry quality frame 1 , 2
	Integer*2	CT_Stat(20)	! Cobetrieve status
	Real*8		Time_diff
	Real*4		Hour /3600./	! Large gap check
	Integer*4	status
	Integer*2	Frm2_Rec
	Parameter	(Frm2_Rec=1)	! record number for frame 2
	Integer*2       STAT_WD1A, STAT_WD1B, STAT_WD2A, STAT_WD2B
	Integer*4	ios		! i/o status
	Logical*1	First_Time /.True./	! First time needs A&B alternate

	FRD_GFG_Read2 = %loc(FRD_Normal)

	TQ1 = 1
	TQ2 = 1
	If (Frm .Eq. 1) Then
	   TQ2 = Hkp_Rec(2).Mj_Frm(2).Tlm_Qual_Maj_Frm
	   If (.Not. TQ2) Then				! 2nd frame was good
	      status = FUT_GMT_CONV(time_sec_2,Hkp_rec(2).Hskp_tail.gmt_mjf2,2)
	      time_diff = time_sec_2 - time_sec_1
	      If (time_diff .GT. hour .AND. Report) Then
	         Write (Lun_rpt,20,iostat=ios) time_diff,Hkp_rec(2).CT_Head.Gmt,
	1           Hkp_rec(2).hskp_tail.Gmt_mjf2
	         If (ios .NE. 0) Then
	            Write (6,25) ios
  25	            Format (1x, ' Error trying to write gap, status = ',I)
	            Write (6,*) ' Gap = ',time_diff, ' secs'
	         EndIf
	      EndIf
	      If (Frm_Q .And. time_diff .GT. 0 .AND. time_diff .LE. Max_sec)Then
	         Frm1_Rec = 1
	         Hkp_Rec(1) = Hkp_Rec(2)
	      Else
	         Do While (TQ1 .And. More_data)
	            Hkp_Rec(1) = Hkp_Rec(2)
	            Call CT_Read_Arcv (, Lun_Hkp, Hkp_Rec(2), CT_Stat )
	            If (CT_Stat(1) .Ne. CTP_Normal) Then
	               If (CT_Stat(1) .Eq. CTP_Endoffile) Then
	                  Call Lib$Signal(FRD_CTHkpEof)
	                  More_Data = .False.
	                  FRD_GFG_Read2 = %loc(FRD_CTHkpEof)
	                  Return
	               Else
	                  FRD_GFG_Read2 = %loc(FRD_Aberr)
	                  Return
	               Endif
	            Else
	               TQ1 = Hkp_Rec(2).Mj_Frm(1).Tlm_Qual_Maj_Frm
	            EndIf
	         EndDo
	         status = FUT_GMT_CONV(time_sec_1,Hkp_rec(2).Ct_Head.gmt,2)
	         TQ2 = Hkp_Rec(1).Mj_Frm(2).Tlm_Qual_Maj_Frm
	         status = FUT_GMT_CONV(time_sec_2,
	1                              Hkp_rec(1).Hskp_tail.gmt_mjf2,2)
	         time_diff = time_sec_1 - time_sec_2
	         If (time_diff .GT. hour .AND. Report) Then
	            Write (Lun_rpt,20,iostat=ios) time_diff,
	1              Hkp_rec(1).hskp_tail.Gmt_mjf2, Hkp_rec(2).CT_Head.Gmt
  20	            Format (1x,'Large housekeeping gap ',F10.0,' sec between ',
	1              A14, ' and ',A14)
	            If (ios .NE. 0) Then
	               Write (6,25) ios
	               Write (6,*) ' Gap = ',time_diff, ' secs'
	            EndIf
	         EndIf
	         If (time_diff .GT. 0 .AND. time_diff .LE. Max_sec .AND.
	1           .Not.TQ2) Then
	            Frm1_Rec = 2
	         Else
	            Hkp_Rec(1) = Hkp_Rec(2)
	            TQ2 = Hkp_Rec(1).Mj_Frm(2).Tlm_Qual_Maj_Frm
	            status =FUT_GMT_CONV(time_sec_2,
	1                              Hkp_rec(1).Hskp_tail.gmt_mjf2,2)
	            time_diff = time_sec_2 - time_sec_1
	            If (time_diff .GT. hour .AND. Report) Then
	               Write (Lun_Rpt,20,iostat=ios) 
	1              time_diff, Hkp_rec(2).CT_Head.Gmt,
	2              Hkp_rec(2).hskp_tail.Gmt_mjf2
	               If (ios .NE. 0) Then
	                  Write (6,25) ios
	                  Write (6,*) ' Gap = ',time_diff, ' secs'
	               EndIf
	            EndIf
	            If (time_diff .GT. 0 .AND. time_diff .LE. Max_sec .AND.
	1              .Not.TQ2) Then
	               Frm1_Rec = 1
	            Else
	               Frm = 2
	               Frm_Q = .False.
	               If (.Not. TQ2) Frm_Q = .True.
	            EndIf
	         EndIf
	      EndIf				! Frm_Q
	   Else
	      Frm = 2
	      Frm_Q = .False.
	   EndIf
	EndIf
	If (Frm .Eq. 2) Then
	   If (Frm_Q) GoTo 100
  200	   Do While (TQ2 .And. More_Data)
	      Call CT_Read_Arcv (, Lun_Hkp, Hkp_Rec(2), CT_Stat )
	      If (CT_Stat(1) .Ne. CTP_Normal) Then
	         If (CT_Stat(1) .Eq. CTP_Endoffile) Then
	            Call Lib$Signal(FRD_CTHkpEof)
	            More_Data = .False.
	            FRD_GFG_Read2 = %loc(FRD_CTHkpEof)
	            Return
	         Else
	            FRD_GFG_Read2 = %loc(FRD_Aberr)
	            Return
	         Endif
	      Else
	         TQ2 = Hkp_Rec(2).Mj_Frm(2).Tlm_Qual_Maj_Frm
	      EndIf
	      If (First_Time) Then

C Make sure both left and right sides are represented in the two frames.

	         STAT_WD1A=HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(1)
	         STAT_WD1B=HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(5)
	         STAT_WD2A=HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(1)
	         STAT_WD2B=HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(5)
	         ACT_MJ_A(1) = ABS(BTEST(STAT_WD1A,14)) + 1
	         ACT_MJ_B(1) = ABS(BTEST(STAT_WD1B,14)) + 1
	         ACT_MJ_A(2) = ABS(BTEST(STAT_WD2A,14)) + 1
	         ACT_MJ_B(2) = ABS(BTEST(STAT_WD2B,14)) + 1
	         If ( Act_MJ_A(1) .EQ. Act_MJ_A(2) .OR. Act_MJ_B(1) .EQ.
	1           Act_MJ_B(2) ) Then
	            TQ2 = 1				! read again
	         EndIf
	      EndIf
	   EndDo
	   TQ1 = Hkp_Rec(2).Mj_Frm(1).Tlm_Qual_Maj_Frm
	   status = FUT_GMT_CONV(time_sec_2,Hkp_rec(2).Hskp_tail.gmt_mjf2,2)
	   status = FUT_GMT_CONV(time_sec_1,Hkp_rec(2).Ct_Head.gmt,2)
	   time_diff = time_sec_2 - time_sec_1
	   If (time_diff .GT. hour .AND. Report) Then
	      Write (Lun_Rpt,20) time_diff,
	1        Hkp_rec(2).CT_Head.Gmt, Hkp_rec(2).hskp_tail.Gmt_mjf2
	      If (ios .NE. 0) Then
	         Write (6,25) ios
	         Write (6,*) ' Gap = ',time_diff, ' secs'
	      EndIf
	   EndIf
	   Hkp_Rec(1) = Hkp_Rec(2)
	   If (time_diff .GT. 0 .AND. time_diff .LE. Max_sec .AND. .Not.TQ1)Then
	      Frm1_Rec = 1
	   Else
  100	      Call CT_Read_Arcv (, Lun_Hkp, Hkp_Rec(2), CT_Stat )
	      If (CT_Stat(1) .Ne. CTP_Normal) Then
	         If (CT_Stat(1) .Eq. CTP_Endoffile) Then
	            Call Lib$Signal(FRD_CTHkpEof)
	            More_Data = .False.
	            FRD_GFG_Read2 = %loc(FRD_CTHkpEof)
	         Else
	            FRD_GFG_Read2 = %loc(FRD_Aberr)
	         Endif
	      Else
	         TQ1 = Hkp_Rec(2).Mj_Frm(1).Tlm_Qual_Maj_Frm
	         If (First_Time) Then
	            STAT_WD1A=HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(1)
	            STAT_WD1B=HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(5)
	            ACT_MJ_A(1) = ABS(BTEST(STAT_WD1A,14)) + 1
	            ACT_MJ_B(1) = ABS(BTEST(STAT_WD1B,14)) + 1
	            If ( Act_MJ_A(1) .EQ. Act_MJ_A(2) .OR. Act_MJ_B(1) .EQ.
	1              Act_MJ_B(2) ) Then
	               TQ1 = 1			! new frame is bad
	            EndIf
	         EndIf
	         status = FUT_GMT_CONV(time_sec_1,Hkp_rec(2).Ct_Head.gmt,2)
	         time_diff = time_sec_1 - time_sec_2
	         If (time_diff .GT. hour .AND. Report) Then
	            Write (Lun_Rpt,20,iostat=ios) 
	1              time_diff, Hkp_rec(1).hskp_tail.Gmt_mjf2,
	2              Hkp_rec(2).CT_Head.Gmt
	            If (ios .NE. 0) Then
	               Write (6,25) ios
	               Write (6,*) ' Gap = ',time_diff, ' secs'
	            EndIf
	         EndIf
	         If (time_diff .GT.0 .And. time_diff .LE. Max_sec .AND.
	1           .Not.TQ1) Then
	            Frm1_Rec = 2
	         Else
	            Hkp_Rec(1) = Hkp_Rec(2)
	            TQ2 = Hkp_Rec(1).Mj_Frm(2).Tlm_Qual_Maj_Frm
	            If (First_Time) Then
	             STAT_WD2A=HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(1)
	             STAT_WD2B=HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(5)
	               ACT_MJ_A(2) = ABS(BTEST(STAT_WD2A,14)) + 1
	               ACT_MJ_B(2) = ABS(BTEST(STAT_WD2B,14)) + 1
	               If ( Act_MJ_A(1) .EQ. Act_MJ_A(2) .OR. Act_MJ_B(1) .EQ.
	1                 Act_MJ_B(2) ) Then
	                  TQ2 = 1			! new frame is bad
	               EndIf
	            EndIf
	            status =FUT_GMT_CONV(time_sec_2,
	1                                Hkp_rec(1).Hskp_tail.gmt_mjf2,2)
	            time_diff = time_sec_2 - time_sec_1
	            If (time_diff .GT. hour .AND. Report) Then
	               Write (Lun_Rpt,20) 
	1                 time_diff, Hkp_rec(2).CT_Head.Gmt,
	2                 Hkp_rec(2).hskp_tail.Gmt_mjf2
	               If (ios .NE. 0) Then
	                  Write (6,25) ios
	                  Write (6,*) ' Gap = ',time_diff, ' secs'
	               EndIf
	            EndIf
	            If (time_diff .GT. 0 .And. time_diff .LE. Max_sec .AND.
	1              .Not.TQ2 .AND. .Not.TQ1) Then
	               Frm1_Rec = 1
	            ElseIf (.Not.TQ2) Then
	               GoTo 100			! get next record's 1st frame
	            Else
	               GoTo 200			! get new record with 2nd frame
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	EndIf
	If (Frm1_Rec .Eq. 1) Then
	   Frm = 1
	Else
	   Frm = 2
	EndIf

C  Determine actual major frame numbers from status monitor.  Word 1 = side
C  A, Word 5 = side B.  Bit test on bit 14  gives 0 for frm 1, 1 for frm 2.
C  These are needed to get gain and fakeit values and times.

	STAT_WD1A = HKP_REC(Frm1_Rec).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(1)
	STAT_WD1B = HKP_REC(Frm1_Rec).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(5)
	STAT_WD2A = HKP_REC(Frm2_Rec).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(1)
	STAT_WD2B = HKP_REC(Frm2_Rec).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(5)
	ACT_MJ_A(1) = ABS(BTEST(STAT_WD1A,14)) + 1
	ACT_MJ_B(1) = ABS(BTEST(STAT_WD1B,14)) + 1
	ACT_MJ_A(2) = ABS(BTEST(STAT_WD2A,14)) + 1
	ACT_MJ_B(2) = ABS(BTEST(STAT_WD2B,14)) + 1
	Frm_Q = .True.
	First_Time = .False.
	Return
	End
