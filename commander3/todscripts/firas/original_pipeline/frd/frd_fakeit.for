C------------------------------------------------------------------------------
	Integer*4 Function FRD_FakeIt (Hkp, Frm, Rec1, Rec2, Act_Mj_A, Act_Mj_B,
	1                              FEX_FakeIt_R_Rec, FEX_FakeIt_L_Rec, Both)
C------------------------------------------------------------------------------
C  Purpose:  Extract Fakeit status from housekeeping record and put into
C            reference record.
C
C  Larry Rosen & Nilo Gonzales, STX, 13-18 March 1991
C  Based on FRD_GAIN and FIRAS USER'S MANUAL (pages 17-21 & 17-22)
C-------------------------------------------------------------------------------
	Implicit None

C  Passed Parameters

	Dictionary 'NFS_HKP'
	Record /NFS_HKP/ HKP(2)

	Dictionary 'FEX_FAKEIT'
	Record /FEX_FAKEIT/ FEX_FAKEIT_R_REC, FEX_FAKEIT_L_REC

	Integer*4  Act_Mj_A(2)      ! Actual A side firas major frames
	Integer*4  Act_Mj_B(2)      ! Actual B side firas major frames

	External   FRD_Normal
	Integer*2  Rec1, Rec2       ! Record number for frame 1 , 2
	Integer*2  Frm              ! Current frame.
	Logical*1  Both             ! True means get both frames, else just frm

C  Local variables

	Integer*2  fake_word_a
	Integer*2  fake_word_b
	Integer*2  fakeit(4)        ! fakeit status for each channel
	Integer*2  i                ! channel number
C
	FRD_FAKEIT = %loc(FRD_Normal)

	If (Frm .Eq. 1 .Or. Both) Then
	   If (act_mj_a(1) .Eq. 1) Then
	      fake_word_a = HKP(rec1).frame(1).hskp_head.stat_monitor_cmd(1)
	      fakeit(1) = IBITS (fake_word_a, 9, 1)        ! RH
	      fakeit(2) = IBITS (fake_word_a, 8, 1)        ! RL
	   EndIf
	   If (act_mj_b(1) .Eq. 1) Then
	      fake_word_b = HKP(rec1).frame(1).hskp_head.stat_monitor_cmd(5)
	      fakeit(3) = IBITS (fake_word_b, 9, 1)        ! LH
	      fakeit(4) = IBITS (fake_word_b, 8, 1)        ! LL
	   EndIf
	EndIf
	If (Frm .Eq. 2 .Or. Both) Then
	   If (act_mj_a(2) .Eq. 1) Then
	      fake_word_a = HKP(rec2).frame(2).hskp_head.stat_monitor_cmd(1)
	      fakeit(1) = IBITS (fake_word_a, 9, 1)        ! RH
	      fakeit(2) = IBITS (fake_word_a, 8, 1)        ! RL
	   EndIf
	   If (act_mj_b(2) .Eq. 1) Then
	      fake_word_b = HKP(rec2).frame(2).hskp_head.stat_monitor_cmd(5)
	      fakeit(3) = IBITS (fake_word_b, 9, 1)        ! LH
	      fakeit(4) = IBITS (fake_word_b, 8, 1)        ! LL
	   EndIf
	EndIf

	Do i=1,2
	   FEX_FAKEIT_R_REC.fakeit(i) = fakeit(i)
	   FEX_FAKEIT_L_REC.fakeit(i) = fakeit(i+2)
	EndDo

	Return
	End
