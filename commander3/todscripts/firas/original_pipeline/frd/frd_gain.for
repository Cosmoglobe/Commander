C------------------------------------------------------------------------------
	Integer*4 Function FRD_Gain (Hkp, Frm, Rec1, Rec2, Act_Mj_A, Act_Mj_B,
	1                            FEX_Gain_R_Rec, FEX_Gain_L_Rec, Both)
C------------------------------------------------------------------------------
C  Purpose:  Extract Gains from housekeeping record and put into reference
C            record.
C
C  Larry Rosen & Nilo Gonzales, STX, 13-18 March 1991
C  Based on old FDQ_INTERPOLATE and FIRAS USER'S MANUAL (pages 17-21 & 17-22)
C-------------------------------------------------------------------------------
	Implicit None

C  Passed Parameters

	Dictionary 'NFS_HKP'
	Record /NFS_HKP/ HKP(2)

	Dictionary 'FEX_GAIN'
	Record /FEX_GAIN/ FEX_GAIN_R_REC, FEX_GAIN_L_REC

	External   FRD_Normal
	Integer*4  Act_Mj_A(2)      ! Actual A side firas major frames
	Integer*4  Act_Mj_B(2)      ! Actual B side firas major frames

	Integer*2  Rec1, Rec2       ! Record number for frame 1 , 2
	Integer*2  Frm              ! Current frame.
	Logical*1  Both             ! True means get both frames, else just frm

C  Local variables

	Integer*2  gain_word_a
	Integer*2  gain_word_b
	Integer*2  gain(4)          ! gain for each channel
	Integer*2  i                ! channel number
C
	FRD_GAIN = %loc(FRD_Normal)

	If (Frm .Eq. 1 .Or. Both) Then
	   If (act_mj_a(1) .Eq. 2) Then
	      gain_word_a = HKP(rec1).frame(1).hskp_head.stat_monitor_cmd(1)
	      call mvbits(gain_word_a,11,1,gain(1),2)		! RH channel
	      call mvbits(gain_word_a,9,1,gain(1),1)		! RH channel
	      call mvbits(gain_word_a,7,1,gain(1),0)		! RH channel
	      call mvbits(gain_word_a,10,1,gain(2),2)		! RL channel
	      call mvbits(gain_word_a,8,1,gain(2),1)		! RL channel
	      call mvbits(gain_word_a,6,1,gain(2),0)		! RL channel
	   EndIf
	   If (act_mj_b(1) .Eq. 2) Then
	      gain_word_b = HKP(rec1).frame(1).hskp_head.stat_monitor_cmd(5)
	      call mvbits(gain_word_b,11,1,gain(3),2)		! LH channel
	      call mvbits(gain_word_b,9,1,gain(3),1)		! LH channel
	      call mvbits(gain_word_b,7,1,gain(3),0)		! LH channel
	      call mvbits(gain_word_b,10,1,gain(4),2)		! LL channel
	      call mvbits(gain_word_b,8,1,gain(4),1)		! LL channel
	      call mvbits(gain_word_b,6,1,gain(4),0)		! LL channel
	   EndIf
	EndIf
	If (Frm .Eq. 2 .Or. Both) Then
	   If (act_mj_a(2) .Eq. 2) Then
	      gain_word_a = HKP(rec2).frame(2).hskp_head.stat_monitor_cmd(1)
	      call mvbits(gain_word_a,11,1,gain(1),2)		! RH channel
	      call mvbits(gain_word_a,9,1,gain(1),1)		! RH channel
	      call mvbits(gain_word_a,7,1,gain(1),0)		! RH channel
	      call mvbits(gain_word_a,10,1,gain(2),2)		! RL channel
	      call mvbits(gain_word_a,8,1,gain(2),1)		! RL channel
	      call mvbits(gain_word_a,6,1,gain(2),0)		! RL channel
	   EndIf
	   If (act_mj_b(2) .Eq. 2) Then
	      gain_word_b = HKP(rec2).frame(2).hskp_head.stat_monitor_cmd(5)
	      call mvbits(gain_word_b,11,1,gain(3),2)		! LH channel
	      call mvbits(gain_word_b,9,1,gain(3),1)		! LH channel
	      call mvbits(gain_word_b,7,1,gain(3),0)		! LH channel
	      call mvbits(gain_word_b,10,1,gain(4),2)		! LL channel
	      call mvbits(gain_word_b,8,1,gain(4),1)		! LL channel
	      call mvbits(gain_word_b,6,1,gain(4),0)		! LL channel
	   EndIf
	EndIf

	Do i=1,2
	   FEX_GAIN_R_REC.gain(i) = gain(i)
	   FEX_GAIN_L_REC.gain(i) = gain(i+2)
	EndDo

	Return
	End
