C--------------------------------------------------------------------------
	Integer*4 Function FRD_Gfg_Check ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
	1              LUN_R_Fakeit, LUN_L_Fakeit, LUN_Rpt, Report, Max_Sec )
C--------------------------------------------------------------------------
C       Purpose: This subroutine checks for telemetry quality, gain change,
C                fakeit change, time gaps. It also calls two subroutines
C                to extract fakeit statuses as well as gain values.
C       
C       Programmer:  Larry Rosen & Nilo Gonzales
C
C       Invocation: Status = FRD_Gfg_Check ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
C                      LUN_R_Fakeit, LUN_L_Fakeit, LUN_Rpt, Report, Max_Sec )
CH      Change Log:
CH		20 Sept. 1991, Larry Rosen, STX.  FGFG bombed when
CH		a 1991 HKP file was put in the read archive.  This caused
CH		a 121 DAY gap, which couldn't be written in the simple
CH		format.  To prevent things like this from bothering FGFG,
CH		we decided to add Jstart Jstop to the command line, and
CH		capture the write status.
C--------------------------------------------------------------------------
C       Input Parameters:
C         Name             Type        Description
C         -----------------------------------------------------------------
C 	  LUN_HKP  	    I*4         ! Unit numbers for NFS_HKP CT archive
C  	  LUN_R_Gain	    I*4         ! Unit numbers for FEX_GAIN_R CT archive
C  	  LUN_L_Gain	    I*4         ! Unit numbers for FEX_GAIN_L CT archive
C	  LUN_R_Fakeit      I*4         ! Unit numbers for FEX_FAKEIT_R archive 
C	  LUN_L_Fakeit      I*4         ! Unit numbers for FEX_FAKEIT_L archive 
C         LUN_Rpt           I*4         ! Report unit number
C         Report            L*1         ! Flag to enable writing a report file
C     	  Max_Sec           I*4         ! Maximum number of seconds allowed
C                                       ! between MJF before a gap is declared.
C       Include Files:
C       
C       CT$Library:CTUser.Inc
C       $SSdef              -- System status values
C
C      Externals:
C
C       FRD_Normal
C       FRD_CTHkpEof
C       FRD_CTReadErr
C
C      Processing Method: PDL for FRD_Gfg_Check
C   
C Begin PDL
C       First time get NFS_HKP record with at least frame 2 = good
C       Set Frame = 2
C       Set Frame Quality = False
C       Call FRD_GFG_Read2   ! Reads housekeeping record in the archive
C       ! and determine Actual major frame numbers from status monitor.
C       ! This information is needed to extract gains, fakeit statuses
C       ! and getting times.
C       Call FILL_Time_Buff  ! Compute offset times and put into time buffer
C       Call FRD_START_FAKE  ! Start a new Fakeit record with times from buffer.
C       Call FRD_START_GAIN  ! Start a new Gain record with times from buffer.
C       Call FRD_Fakeit      ! Extract fakeit status value from the HKP record
C       Call FRD_Gain        ! Extract Gain value from the HKP record
C 	Set First flag = False
C
C       Do While ( More_Data .And. Return status is Normal)
C          Set Next_Frm = False
C	   If (Frame .Eq. 1) Then
C	       Call CT_Read_Arcv         ! Read a new HKP record
C     	       Set Frame_1_Record = 2
C	       If (Read status is not Normal) Then
C    	          If (Read status is End of File) Then
C	             Signal end of file.
C  	             Set More_Data = False
C	             Set Gain Record End_of_data Flag
C	             Set Fakeit_Record End_of_dat Flag
C	             Call FRD_END_FAKE      ! End fakeit record
C	             If (FRD_END_FAKE Return status is Normal) Then
C                       Signal Return status to Normal
C	             Else
C                       Add 1 to Reference Fakeit record counter
C	       	        Call FRD_END_GAIN   ! End gain record
C	                If (FRD_END_GAIN Return status is Normal) Then
C	                   Signal Return status to Normal
C                       Else
C	                   Add 1 to Reference Gain record counter 
C               	   Signal End of file of the housekeeping record
C	                Endif
C	          Else
C	              Signal Cobetrieve read error
C                     Set More_data = False
C      	          Endif
C              Else
C      	           Actual major frame 1 from status monitor
C	           Call FUT_GMT_CONV     ! Convert GMT to time in seconds
C	       Endif
C	   Else                          ! Frame = 2
C	       Actual major frames 2 from status monitor
C     	       Call FUT_GMT_CONV         ! Convert GMT to time in seconds
C	   Endif                         ! Frame 1 or 2
C      	   If (FRD_GFG_Check status is Normal) Then
C	      Call Fill_Time_Buff        ! Compute offset times and put into
C                                        ! time buffer    
C     Check Telemetry Quality.
C             Telemetry Quality of the new major frame
C   	      If (Telemetry Quality of the current frame is bad) Then
C	         Set Next_Frame flag = True
C	         Set Fakeit Telemetry Quality Flag = 1 
C	         Set Gain   Telemetry Quality Flag = 1 
C	         Set Frame Quality = False
C                Call FRD_END_FAKE           ! End fakeit record
C                If (FRD_END_FAKE return status is not Normal) Then
C	             Signal FRD_END_FAKE return status to not Normal
C                Else
C	            Add 1 to Fakeit record counter 
C	            Store Fakeit GMT time into Temporary holding area
C	            Call FRD_END_GAIN        ! End gain record
C	            If (Return status is not Normal) Then
C	               Signal FRD_END_GAIN return status to not Normal
C	            Else
C	               Add 1 to Gain record counter
C	               Store Gain GMT time into Temporary holding area
C	               Add 1 to Telemetry counter 
C	            Endif
C    	          Endif
C     Begin reading HKP records until two good major frames are found.
C      	          If (FRD_GFG_Check status is Normal) Then
C	             Call FRD_GFG_Read2   ! Read new HKP frames and determine
C		     ! Actual major frame numbers from status monitor.  This
C       	     ! This information is needed to extract gains, fakeit
C		     ! statuses and getting times.
C	             If (Read status is not Normal) Then
C	                 Signal read status to not Normal
C	             Else
C      	                 Call FILL_Time_Buff  ! Compute offset times and
C	                                      ! put into time buffer.
C     Start the reference file records and reset.
C	               	 Fakeit record 1 = Fakeit record 2
C	               	 Gain record 1   = Gain record 2
C                        Call FRD_START_FAKE  ! Start a new Fakeit record with
C	                                      ! times from buffer.
C                        Call FRD_START_GAIN  ! Start a new Gain record with 
C	                                      ! times from buffer.
C                        Call FRD_Fakeit      ! Extract fakeit value from HKP
C                        Call FRD_Gain        ! Extract Gain value from HKP
C	              Endif                   ! return status from READ2
C              Endif                          ! Telemetry quality check
C     Check for Gap between Major Frames.
C     	       If (Next_Frame flag is false and FRD_GFG_Check status is Normal)
C	       Then
C                 If (Frame = 1) Then
C	             Calculate Time Difference between Frames
C	             If (Time difference is greater than 1 hour and Report is
C	                 selected) 
C      	                 Write 'Large Housekeeping Gap into the report file'
C            	  Else 
C	             Calculate Time Difference between Frames
C      	                 Write 'Large Housekeeping Gap into the report file'
C	          Endif
C	          If (Calculated Time Difference is equal Gap) Then 
C	              Set Next Frame = True
C	              Set Fakeit data gap flag = 1
C	              Set Gain data gap flag = 1
C	              Add 1 to Gap Fakeit Counter
C	              Add 1 to Gap Gain Counter     	              
C	              Call FRD_END_FAKE      ! End fakeit record
C                     If (FRD_END_FAKE return status is not Normal) Then
C	                  Signal FRD_END_FAKE return status to not Normal
C	              Else
C	                  Add 1 to Fakeit record counter
C	                  Store Fakeit GMT time into Temporary holding area
C	                  Call FRD_END_GAIN   ! End fakeit record
C                         If (FRD_END_GAIN return status is not Normal) Then
C	                     Signal FRD_END_GAIN return status to not Normal
C	                  Else
C	                      Add 1 to Gain record counter
C	                      Store Gain GMT time into Temporary holding area
C	                  Endif
C	              Endif
C
C     Current frame has good Telemetry quality, but was preceded by a gap.
C     Check next frame and if it is not good, then read until 2 good frames are
C     found. This is an option of Function FRD_GFG_Read2.
C 
C	               Call FRD_GFG_Read2   ! Read new HKP frames and determine
C		     ! Actual major frame numbers from status monitor.  This
C       	     ! This information is needed to extract gains, fakeit
C		     ! statuses and getting times.
C	               If (Read status is not Normal) Then
C	                   Signal read status to not Normal
C	               Else
C	                   Call Fill_Time_Buff  ! Compute offset times 
C                                               ! and put into time buffer
C     Start the reference file records and reset.
C	               	   Fakeit record 1 = Fakeit record 2
C	               	   Gain record 1   = Gain record 2
C                          Call FRD_START_FAKE  ! Start a new Fakeit
C	                                        ! record with times from buffer.
C                          Call FRD_START_GAIN  ! Start a new Gain record with 
C	                                        ! times from buffer.
C                          Call FRD_Fakeit      ! Extract fakeit value from HKP
C                          Call FRD_Gain        ! Extract Gain value from HKP
C	               Endif                    ! return status from READ2
C                 Endif                         ! If there is a gap.
C              Endif                            ! Not Next frame and status
C     If no gap and no bad telemetry quality for this frame. Both flag is set
C      to extract Gain and Fakeit frames. 
C	       If (Next Frame flag is False and FRD_GFG_Check is Normal) Then
C                  Call FRD_Fakeit            ! Extract fakeit value from HKP
C                  Call FRD_Gain              ! Extract Gain value from HKP
C     Check if Fakeit data changes. Set flags and increase counter for each 
C      Channel.
C              Do Channel 1 through 4
C                 If (Fakeit buffer2(Channel) .NE. Fakeit buffer1(Channel) Then 
C             	     Set Fakeit Change flag = 1
C             	     Add 1 to Channel fakeit Counter
C	          Endif
C	       EndDo
C              If (Any Fakeit Change flags are set) Then
C	       Call FRD_END_FAKE          ! End fakeit record
C	       If (FRD_END_FAKE Return status is not Normal) Then
C                  Signal Return status to not Normal
C	       Else
C	           Add 1 to Fakeit record counter
C	           Store Fakeit GMT time into Temporary holding area
C	           Fakeit record 1 = Fakeit record 2
C                  Call FRD_START_FAKE     ! Start a new Fakeit record with
C                                          ! times from buffer.
C	       Endif
C     Check if Gain data changes. Set flags and increase counter for each 
C      Channel.
C              Do Channel 1 through 4
C                 If (Gain buffer2(Channel) .NE. Gain  buffer1(Channel) Then 
C             	     Set Gain  Change flag = 1
C             	     Add 1 to Channel Gain  Counter
C	          Endif
C	       EndDo
C              If (Any Fakeit Change flags are set) Then
C	       Call FRD_END_GAIN            ! End gain record
C	       If (FRD_END_GAIN Return status is not Normal) Then
C                  Signal Return status to not Normal
C	       Else
C	           Add 1 to Gain record counter
C	           Store Gain GMT time into Temporary holding area
C	           Gain record 1 = Gain record 2
C                  Call FRD_START_GAIN  ! Start a new Gain record with
C                                       ! times from buffer.
C	       Endif                  
C              Time buffer1 = Time buffer2 
C	       Housekeeping rec 1 = Housekeeping rec 2
C              If (Frame is equal 2) Then
C     	          Set Frame = 1
C	       Else
C	          Set Frame = 2
C	       Endif
C          Endif                         ! If Next frame is false and status
C       Enddo                            ! More data
C	If (report) Then                 ! If Report qualifier is selected
C           Write ' Information to report file.'
C       Endif
C       Return
C       End
C ENDPDL
C--------------------------------------------------------------------------
	Implicit None

C	Passed Parameters:

	Integer*4	LUN_HKP  	! Unit numbers for NFS_HKP CT archive
	Integer*4	LUN_R_Gain	! Unit numbers for FEX_GAIN_R CT archive
	Integer*4	LUN_L_Gain	! Unit numbers for FEX_GAIN_L CT archive
	Integer*4       LUN_R_Fakeit    ! Unit numbers for FEX_FAKEIT_R archive
	Integer*4       LUN_L_Fakeit    ! Unit numbers for FEX_FAKEIT_L archive
	Character*21    Report_File     ! Report file name 
	Integer*4	LUN_Rpt         ! Report unit number
	Logical*1       Report          ! Flag to enable writing a report file
	Integer*4	Max_sec		! Max number secs before gap declared

C 	Include files:

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'
	Include '(FUT_Params)'
	Include '(FUT_Error)'

C      Externals:

	External        FRD_Normal
	External        FRD_CTHkpEof
	External        FRD_CTReadErr

C  	Functions:

	Integer*4       Lib$Signal  
	Integer*4       FRD_Fakeit
	Integer*4       FRD_Gain
	Integer*4	FUT_GMT_CONV
	Integer*4	FRD_START_FAKE          ! Start a new fakeit record
	Integer*4	FRD_START_GAIN          ! Start a new gain record
	Integer*4	FRD_END_FAKE		! end and write fakeit record
	Integer*4	FRD_END_GAIN		! end and write gain record
	Integer*4	Fill_Time_Buff		! fill time buffers
	Integer*4	FRD_GFG_Read2		! Reads HKP records until 2 good
						! frames are found (maybe not in
						! the same record).

C     	Local Variables:

        Dictionary 'NFS_HKP'
        Record /NFS_HKP/    HKP_REC(2)

	Dictionary 'FEX_FAKEIT'
	Record /FEX_FAKEIT/ Fex_Fakeit_R_Rec(2)
	Record /FEX_FAKEIT/ Fex_Fakeit_L_Rec(2)

	Dictionary 'FEX_GAIN'
	Record /FEX_GAIN/ Fex_Gain_R_Rec(2)
	Record /FEX_GAIN/ Fex_Gain_L_Rec(2)

	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ Time_Buff_R(2)  ! 1 = previous HKP, 2 = current, R
	Record /Time_Buffer/ Time_Buff_L(2)  ! 1 = previous HKP, 2 = current, L

	Real*4		Hour /3600./	! Secs in 1 hour for report large gaps.
	Integer*4 	RETSTAT		! Return status from function
	Integer*4	STATUS		! temporary status holder
	Integer*4	ios		! i/o status
        Integer*4       Fakeit_Rec_R_Ctr /0/  ! Fakeit record counter
        Integer*4       Fakeit_Rec_L_Ctr /0/  ! Fakeit record counter
        Integer*4       Gain_Rec_R_Ctr /0/    ! Gain record counter
        Integer*4       Gain_Rec_L_Ctr /0/    ! Gain record counter
        Integer*4       Fakeit_Chan_Ctr(4)  ! Number of Fakeit changes per chan
        Integer*4       Gain_Chan_Ctr(4)    ! Number of Gain changes per chan
        Integer*4       Gap_Fakeit_Ctr /0/  ! Gap Fakeit record counter
        Integer*4       Gap_Gain_Ctr /0/    ! Gap Gain record counter
	Integer*4       TQ_Ctr      /0/     ! Bad Telemetry counter
	Integer*2	Chan		    ! index, channel number
        Integer*4       Act_Mj_A(2)     ! Actual a side major frame
        Integer*4       Act_mj_B(2)     ! Actual b side major frame
	Logical*1       More_data /.True./      ! Flag for more data
	Integer*2       iyr,iday,ihr,imin
	Integer*4	msec
	Integer*2       STAT_WD1A, STAT_WD1B, STAT_WD2A, STAT_WD2B
	Logical*1	Next_Frm /.True./   ! Do main loop again
	Logical*1	Both                ! Get both frames gain and fakeit.
	Equivalence (Next_Frm, Both)        ! Always want both when flags set.
	Logical*1       First /.True./	    ! First time function called
	Integer*2	TQ		! Telemetry Quality current frame 0=good
	Real*8		time_sec_2, time_sec_1, time_diff
	Real*8          time_sec	 ! Temporary storage of time in secs
	Character*14    Prev_Gmt_Hold_F_R  ! Temp store previous FEX GMT Fakeit
	Character*14    Prev_Gmt_Hold_F_L  ! Temp store previous FEX GMT Fakeit
	Character*14    Prev_Gmt_Hold_G_R  ! Temp store previous FEX GMT Gain
	Character*14    Prev_Gmt_Hold_G_L  ! Temp store previous FEX GMT Gain
	Integer*2	CT_Stat(20)	! Cobetrieve status
	Integer*2	Frm		! Current HKP frame examined (1 or 2)
	Logical*1	Frm_Q		! Quality of current frame. True = good
	Integer*2	Frm1_Rec	! record number for frame 1
	Integer*2	Frm2_Rec
	Parameter	(Frm2_Rec=1)	! record number for frame 2
C
C  case 1:  Frm2_Rec = 1    [ 1 | 2 ]     case 2:  Frm2_Rec = 1     [   | 2 ]
C           Frm1_Rec = 1    	                   Frm1_Rec = 2     [ 1 |   ]

C	Set the function status to normal.

	FRD_GFG_Check = %loc(FRD_Normal)
	Retstat = %loc(FRD_Normal)

C  First time get first housekeeping record with at least frame 2 = good.

	Frm = 2
	Frm_Q = .False.
	Retstat = FRD_GFG_Read2 ( Lun_hkp, Hkp_Rec, More_Data, Time_sec_1,
	1  Time_sec_2, Max_Sec, Frm1_Rec, Frm, Frm_Q, Lun_Rpt, Report,
	2  ACT_MJ_A, ACT_MJ_B )

	If (Retstat .Ne. %loc(FRD_Normal)) Then
	   FRD_GFG_Check = Retstat
	Else

C  Write first frame time found to report file.

	   If (Report) Then
	      If (Frm1_Rec .eq. 1) Then
	         Write (Lun_Rpt,10) Hkp_Rec(1).ct_head.gmt
	      Else
	         Write (Lun_Rpt,10) Hkp_Rec(1).hskp_tail.gmt_mjf2
	      End If
  10	      Format(1x,/,1x,'Time of First Major Frame Found: ',A14)
	   End If	      

C  Based on actual frame number, put offset times into time buffer.

	   status = Fill_Time_Buff (act_mj_a, act_mj_b, time_buff_R(1),
	1     time_buff_L(1), time_sec_1, time_sec_2)

C  We now have all the times for gain and fakeit for this HKP record.
C  Start the reference file records with them.

	   status = FRD_START_FAKE ( Fex_Fakeit_R_rec(1),
	1     Time_buff_R(1), prev_gmt_hold_F_R, First )
	   status = FRD_START_FAKE ( Fex_Fakeit_L_rec(1),
	1     Time_buff_L(1), prev_gmt_hold_F_L, First )
	   status = FRD_START_GAIN ( Fex_Gain_R_rec(1),
	1     Time_buff_R(1), prev_gmt_hold_G_R, First )
	   status = FRD_START_GAIN ( Fex_Gain_L_rec(1),
	1     Time_buff_L(1), prev_gmt_hold_G_L, First )

C  Now get gain and fake_it data

	   status = FRD_Fakeit (HKP_Rec, Frm, Frm1_Rec, Frm2_Rec, Act_Mj_A,
	1     Act_Mj_B, Fex_Fakeit_R_Rec(1), Fex_Fakeit_L_Rec(1), Both)
	   status =  FRD_Gain (HKP_Rec, Frm, Frm1_Rec, Frm2_Rec, Act_Mj_A,
	1     Act_Mj_B, Fex_Gain_R_Rec(1), Fex_Gain_L_Rec(1), Both )

	   First = .False.

	EndIf			 ! Retstat from FRD_GFG_Read2

C  While there is more data: read the HKP records and get gain/fakeits as above.
C  MAIN LOOP

	Do While (More_Data .AND. (Retstat .EQ. %loc(FRD_Normal)))

	   Next_Frm = .False.				! also sets Both = false
	   If (Frm .Eq. 1) Then				! get new record
	      Call CT_Read_Arcv (, Lun_Hkp, Hkp_Rec(2), CT_Stat )
	      Frm1_Rec = 2
	      If (CT_Stat(1) .Ne. CTP_Normal) Then
	         If (CT_Stat(1) .Eq. CTP_Endoffile) Then
	            Call Lib$Signal(FRD_CTHkpEof)
	            More_data = .False.
	            Fex_Gain_R_Rec(1).End_of_data = .True.
	            Fex_Gain_L_Rec(1).End_of_data = .True.
	            Fex_Fakeit_R_Rec(1).End_of_data = .True.
	            Fex_Fakeit_L_Rec(1).End_of_data = .True.
	            Retstat = FRD_END_FAKE ( Fex_Fakeit_R_Rec(1),
	1             Time_buff_R(1), Lun_R_Fakeit )
	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Retstat = FRD_END_FAKE ( Fex_Fakeit_L_Rec(1),
	1                 Time_buff_L(1), Lun_L_Fakeit )
	               If (Retstat .NE. %loc(FRD_Normal)) Then
	                  FRD_GFG_Check = Retstat
	               Else
	                  Fakeit_Rec_R_Ctr = Fakeit_Rec_R_Ctr + 1
	                  Fakeit_Rec_L_Ctr = Fakeit_Rec_L_Ctr + 1
	                  Retstat = FRD_END_GAIN ( Fex_Gain_R_Rec(1),
	1                    Time_buff_R(1), Lun_R_Gain )
	                  If (Retstat .NE. %loc(FRD_Normal)) Then
	                     FRD_GFG_Check = Retstat
	                  Else
	                     Retstat = FRD_END_GAIN ( Fex_Gain_L_Rec(1),
	1                       Time_buff_L(1), Lun_L_Gain )
	                     If (Retstat .NE. %loc(FRD_Normal)) Then
	                        FRD_GFG_Check = Retstat
	                     Else
	                        Gain_Rec_R_Ctr = Gain_Rec_R_Ctr + 1
	                        Gain_Rec_L_Ctr = Gain_Rec_L_Ctr + 1
	                        FRD_GFG_Check = %loc(FRD_CTHkpEof)
	                     EndIf
	                  EndIf
	               EndIf
	            EndIf			! If end of file
	         Else
	            Call Lib$Signal(FRD_ctreaderr)
	            FRD_GFG_Check = %loc(FRD_CTREADERR)
	            More_data = .False.
	         Endif
	      Else				!  good read of HKP record

C  Determine actual major frame numbers from status monitor.
C  Get time in seconds, add offsets for gain and fakeit times, store in timebuff

	         STAT_WD1A = HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(1)
	         STAT_WD1B = HKP_REC(2).FRAME(1).HSKP_HEAD.STAT_MONITOR_CMD(5)
	         ACT_MJ_A(1) = ABS(BTEST(STAT_WD1A,14)) + 1
	         ACT_MJ_B(1) = ABS(BTEST(STAT_WD1B,14)) + 1
	         status = FUT_GMT_CONV (time_sec_1, Hkp_rec(2).Ct_Head.gmt, 2)
	      EndIf				! Read of HKP record
	   Else					! Frm = 2
	      STAT_WD2A = HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(1)
	      STAT_WD2B = HKP_REC(2).FRAME(2).HSKP_HEAD.STAT_MONITOR_CMD(5)
	      ACT_MJ_A(2) = ABS(BTEST(STAT_WD2A,14)) + 1
	      ACT_MJ_B(2) = ABS(BTEST(STAT_WD2B,14)) + 1
	      status = FUT_GMT_CONV (time_sec_2, 
	1                             Hkp_rec(2).hskp_tail.Gmt_mjf2,2)
	   EndIf				! Frm = 1 or 2

	   If (FRD_GFG_Check .Eq. %loc(FRD_Normal)) Then

C  Based on actual frame number, put offset times into time buffer.

	      status = Fill_Time_Buff (act_mj_a, act_mj_b, time_buff_R(2),
	1        time_buff_L(2), time_sec_1, time_sec_2)
   
C  CHECK TELEMETRY QUALITY of new frame.

	      TQ = Hkp_Rec(2).Mj_Frm(Frm).Tlm_Qual_Maj_Frm
	      If (TQ) Then
	         Next_Frm = .True.			! Both = true
	         Fex_Fakeit_R_Rec(1).Tlm_Bad_Quality = 1
	         Fex_Fakeit_L_Rec(1).Tlm_Bad_Quality = 1
	         Fex_Gain_R_Rec(1).Tlm_Bad_Quality = 1
	         Fex_Gain_L_Rec(1).Tlm_Bad_Quality = 1
	         Frm_Q = .False.
	         Retstat = FRD_END_FAKE ( Fex_Fakeit_R_Rec(1),
	1           Time_buff_R(1), Lun_R_Fakeit )
 	         If (Retstat .NE. %loc(FRD_Normal)) Then
	            FRD_GFG_Check = Retstat
	         Else
	            Retstat = FRD_END_FAKE ( Fex_Fakeit_L_Rec(1),
	1              Time_buff_L(1), Lun_L_Fakeit )
 	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Fakeit_Rec_R_Ctr = Fakeit_Rec_R_Ctr + 1
	               Fakeit_Rec_L_Ctr = Fakeit_Rec_L_Ctr + 1
	               Prev_Gmt_Hold_F_R = Fex_Fakeit_R_Rec(1).CT_head.Gmt
	               Prev_Gmt_Hold_F_L = Fex_Fakeit_L_Rec(1).CT_head.Gmt
	               Retstat = FRD_END_GAIN ( Fex_Gain_R_Rec(1),
	1                 Time_buff_R(1), Lun_R_Gain )
	               If (Retstat .NE. %loc(FRD_Normal)) Then
	                  FRD_GFG_Check = Retstat
	               Else
	                  Retstat = FRD_END_GAIN ( Fex_Gain_L_Rec(1),
	1                    Time_buff_L(1), Lun_L_Gain )
	                  If (Retstat .NE. %loc(FRD_Normal)) Then
	                     FRD_GFG_Check = Retstat
	                  Else
	                     Gain_Rec_R_Ctr = Gain_Rec_R_Ctr + 1
	                     Gain_Rec_L_Ctr = Gain_Rec_L_Ctr + 1
	                     Prev_Gmt_Hold_G_R = Fex_Gain_R_Rec(1).CT_head.Gmt
	                     Prev_Gmt_Hold_G_L = Fex_Gain_L_Rec(1).CT_head.Gmt
	                     TQ_Ctr = TQ_Ctr + 1
	                  EndIf
	               EndIf
	            EndIf
	         EndIf

C  Begin reading HKP records until we get two good consecutive frames.

	         If (FRD_GFG_Check .Eq. %loc(FRD_Normal))
	1           Retstat = FRD_GFG_Read2 ( Lun_hkp, Hkp_Rec, More_Data,
	2              Time_sec_1, Time_sec_2, Max_Sec, Frm1_Rec, Frm, Frm_Q,
	3              Lun_Rpt, Report, ACT_MJ_A, ACT_MJ_B )

	         If (Retstat .Ne. %loc(FRD_Normal)) Then
	            FRD_GFG_Check = Retstat
	         Else
	            status = Fill_Time_Buff (act_mj_a, act_mj_b, time_buff_R(1),
	1              time_buff_L(1), time_sec_1, time_sec_2)

C  Start the reference file records.  Clear first by setting Fex_rec1 = Fex_rec2

	            Fex_Fakeit_R_rec(1) = Fex_Fakeit_R_rec(2)
	            Fex_Fakeit_L_rec(1) = Fex_Fakeit_L_rec(2)
	            Fex_Gain_R_rec(1) = Fex_Gain_R_rec(2)
	            Fex_Gain_L_rec(1) = Fex_Gain_L_rec(2)
	            status = FRD_START_FAKE ( Fex_Fakeit_R_rec(1),
	1              Time_buff_R(1), prev_gmt_hold_F_R, First )
	            status = FRD_START_FAKE ( Fex_Fakeit_L_rec(1),
	1              Time_buff_L(1), prev_gmt_hold_F_L, First )
	            status = FRD_START_GAIN ( Fex_Gain_R_rec(1),
 	1              Time_buff_R(1), prev_gmt_hold_G_R, First )
	            status = FRD_START_GAIN ( Fex_Gain_L_rec(1),
 	1              Time_buff_L(1), prev_gmt_hold_G_L, First )
	            status = FRD_Fakeit (HKP_Rec, Frm, Frm1_Rec, Frm2_Rec,
	1              Act_Mj_A, Act_Mj_B, Fex_Fakeit_R_Rec(1),
	2              Fex_Fakeit_L_Rec(1), Both)
	            status =  FRD_Gain ( HKP_Rec, Frm, Frm1_Rec, Frm2_Rec,
	1              Act_Mj_A, Act_Mj_B, Fex_Gain_R_Rec(1),
	2              Fex_Gain_L_Rec(1), Both)

	         EndIf				! return status from read2
	      Endif				! TQ check

C  Check for Gap between Major Frames.

	      If (.Not.Next_Frm .And. FRD_GFG_Check .Eq. %loc (FRD_Normal))
	1     Then
	         If (Frm .Eq. 1) Then
	            time_diff = time_sec_1 - time_sec_2
	            If (time_diff .GT. hour .AND. Report) Then
	               Write (Lun_rpt,20,Iostat=ios) time_diff,
	2                    Hkp_rec(1).hskp_tail.Gmt_mjf2,
	3                    Hkp_rec(2).CT_Head.Gmt
  20	               Format (1x,'Large housekeeping gap ',F10.0,
	1                    ' sec between ', A14,' and ',A14)
	               If (ios .NE. 0) Then
	                  Write (6,25) ios
  25	                  Format (1x, ' Error trying to write gap, status = ',I)
	                  Write (6,*) ' Gap = ',time_diff, ' secs'
	               EndIf
	            EndIf
	         Else				! Frm = 2
	            time_diff = time_sec_2 - time_sec_1
	            If (time_diff .GT. hour .AND. Report) Then
	               Write (Lun_rpt,20,Iostat=ios) time_diff,
	2                 Hkp_rec(2).CT_Head.Gmt, Hkp_rec(2).hskp_tail.Gmt_mjf2
	               If (ios .NE. 0) Then
	                  Write (6,25) ios
	                  Write (6,*) ' Gap = ',time_diff, ' secs'
	               EndIf
                    EndIf
	         EndIf
	         If (FRD_GFG_Check .EQ. %loc(FRD_Normal)) Then
 	            If ( time_diff .LE. 0 .OR. time_diff .GT. Max_sec ) Then
	               Next_Frm = .True.
	               Fex_Fakeit_R_Rec(1).fakeit_data_gap = 1
	               Fex_Fakeit_L_Rec(1).fakeit_data_gap = 1
	               Fex_Gain_R_Rec(1).gain_data_gap = 1
	               Fex_Gain_L_Rec(1).gain_data_gap = 1
	               Gap_Fakeit_Ctr = Gap_Fakeit_Ctr + 1
	               Gap_Gain_Ctr = Gap_Gain_Ctr + 1
	               Retstat = FRD_END_FAKE ( Fex_Fakeit_R_Rec(1),
	1                 Time_buff_R(1), Lun_R_Fakeit )
	               If (Retstat .NE. %loc(FRD_Normal)) Then
	                  FRD_GFG_Check = Retstat
	               Else
	                  Retstat = FRD_END_FAKE ( Fex_Fakeit_L_Rec(1),
	1                    Time_buff_L(1), Lun_L_Fakeit )
	                  If (Retstat .NE. %loc(FRD_Normal)) Then
	                     FRD_GFG_Check = Retstat
	                  Else
	                     Fakeit_Rec_R_Ctr = Fakeit_Rec_R_Ctr + 1
	                     Fakeit_Rec_L_Ctr = Fakeit_Rec_L_Ctr + 1
	                     Prev_Gmt_Hold_F_R = Fex_Fakeit_R_Rec(1).CT_head.Gmt
	                     Prev_Gmt_Hold_F_L = Fex_Fakeit_L_Rec(1).CT_head.Gmt
	                     Retstat = FRD_END_GAIN ( Fex_Gain_R_Rec(1),
	1                       Time_buff_R(1), Lun_R_Gain )
	                     If (Retstat .NE. %loc(FRD_Normal)) Then
	                        FRD_GFG_Check = Retstat
	                     Else
	                        Retstat = FRD_END_GAIN ( Fex_Gain_L_Rec(1),
	1                          Time_buff_L(1), Lun_L_Gain )
	                        If (Retstat .NE. %loc(FRD_Normal)) Then
	                           FRD_GFG_Check = Retstat
	                        Else
	                           Gain_Rec_R_Ctr = Gain_Rec_R_Ctr + 1
	                           Gain_Rec_L_Ctr = Gain_Rec_L_Ctr + 1
	                           Prev_Gmt_Hold_G_R =
	1                             Fex_Gain_R_Rec(1).CT_head.Gmt
	                           Prev_Gmt_Hold_G_L =
	1                             Fex_Gain_L_Rec(1).CT_head.Gmt
	                        EndIf
	                     EndIf
	                  EndIf
	               EndIf

C  Current frame has good TQ, but was preceded by a gap.  So, check next frame
C  and if not ok then read records until 2 good frames. This is an option of
C  function FRD_GFG_Read2.

	               Retstat = FRD_GFG_Read2 ( Lun_hkp, Hkp_Rec, More_Data,
	1                 Time_sec_1, Time_sec_2, Max_Sec, Frm1_Rec, Frm, Frm_Q,
	2                 Lun_Rpt, Report, ACT_MJ_A, ACT_MJ_B )

	               If (Retstat .Ne. %loc(FRD_Normal)) Then
	                  FRD_GFG_Check = Retstat
	               Else
	                  status = Fill_Time_Buff (act_mj_a,act_mj_b,
	1                 time_buff_R(1), time_buff_L(1), time_sec_1,time_sec_2)

C  Start the reference file records.  Clear first by setting Fex_rec1 = Fex_rec2

	                  Fex_Fakeit_R_rec(1) = Fex_Fakeit_R_rec(2)
	                  Fex_Fakeit_L_rec(1) = Fex_Fakeit_L_rec(2)
	                  Fex_Gain_R_rec(1) = Fex_Gain_R_rec(2)
	                  Fex_Gain_L_rec(1) = Fex_Gain_L_rec(2)
	                  status = FRD_START_FAKE ( Fex_Fakeit_R_rec(1),
	1                    Time_buff_R(1), prev_gmt_hold_F_R, First )
	                  status = FRD_START_FAKE ( Fex_Fakeit_L_rec(1),
	1                    Time_buff_L(1), prev_gmt_hold_F_L, First )
	                  status = FRD_START_GAIN ( Fex_Gain_R_rec(1),
	1                    Time_buff_R(1), prev_gmt_hold_G_R, First )
	                  status = FRD_START_GAIN ( Fex_Gain_L_rec(1),
	1                    Time_buff_L(1), prev_gmt_hold_G_L, First )
	                  status = FRD_Fakeit (HKP_Rec, Frm, Frm1_Rec, Frm2_Rec,
	1                    Act_Mj_A, Act_Mj_B, Fex_Fakeit_R_Rec(1),
	2                    Fex_Fakeit_L_Rec(1), Both)
	                  status =  FRD_Gain ( HKP_Rec, Frm, Frm1_Rec, Frm2_Rec,
	1                    Act_Mj_A, Act_Mj_B, Fex_Gain_R_Rec(1),
	2                    Fex_Gain_L_Rec(1), Both)
	               EndIf			! return status from read2
	            EndIf			! If there is a gap
	         EndIf				! If FRD_GFG_Check .Eq. Normal
	      EndIf				! If .not.Next_Frm and status

C  If no gap and no bad TQ for this frame. Get gain and fakeit data. Both = F

	      If (.Not. Next_Frm .And. FRD_GFG_Check.Eq.%loc(FRD_Normal)) Then
	         status = FRD_Fakeit ( HKP_Rec, Frm, Frm1_Rec, Frm2_Rec,
	1           Act_Mj_A, Act_Mj_B, Fex_Fakeit_R_Rec(2),
	2           Fex_Fakeit_L_Rec(2), Both)
	         status = FRD_Gain ( HKP_Rec, Frm, Frm1_Rec, Frm2_Rec, Act_Mj_A,
	1           Act_Mj_B, Fex_Gain_R_Rec(2), Fex_Gain_L_Rec(2), Both)

C  Check if fakeit data changes. Set flags and increase counter for each chan.

	         Do Chan = 1,2
	            If ( Fex_Fakeit_R_Rec(2).Fakeit(chan)    .Ne.
	1                Fex_Fakeit_R_Rec(1).Fakeit(chan) )  Then
	               Fex_Fakeit_R_Rec(1).Fakeit_Change(chan) = 1
	               Fakeit_chan_ctr(chan) = Fakeit_chan_ctr(chan) + 1
	            EndIf
	            If ( Fex_Fakeit_L_Rec(2).Fakeit(chan)    .Ne.
	1                Fex_Fakeit_L_Rec(1).Fakeit(chan) )  Then
	               Fex_Fakeit_L_Rec(1).Fakeit_Change(chan) = 1
	               Fakeit_chan_ctr(chan+2) = Fakeit_chan_ctr(chan+2) + 1
	            EndIf
	         EndDo

C  If one or more fakeit changed.

	         If (Fex_Fakeit_R_Rec(1).Fakeit_Change(1) +
	1           Fex_Fakeit_R_Rec(1).Fakeit_Change(2) .GT. 0) Then
	            Retstat = FRD_END_FAKE ( Fex_Fakeit_R_Rec(1),
	1              Time_buff_R(1), Lun_R_Fakeit )
	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Fakeit_Rec_R_Ctr = Fakeit_Rec_R_Ctr + 1
	               Prev_Gmt_Hold_F_R = Fex_Fakeit_R_Rec(1).CT_head.Gmt
	               Fex_Fakeit_R_rec(1) = Fex_Fakeit_R_rec(2)
	               status = FRD_START_FAKE ( Fex_Fakeit_R_rec(1),
	1                  Time_buff_R(2), prev_gmt_hold_F_R, First )
	            EndIf
	         EndIf					! If Fakeit changes

	         If (Fex_Fakeit_L_Rec(1).Fakeit_Change(1) +
	1           Fex_Fakeit_L_Rec(1).Fakeit_Change(2) .GT. 0) Then
	            Retstat = FRD_END_FAKE ( Fex_Fakeit_L_Rec(1),
	1              Time_buff_L(1), Lun_L_Fakeit )
	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Fakeit_Rec_L_Ctr = Fakeit_Rec_L_Ctr + 1
	               Prev_Gmt_Hold_F_L = Fex_Fakeit_L_Rec(1).CT_head.Gmt
	               Fex_Fakeit_L_rec(1) = Fex_Fakeit_L_rec(2)
	               status = FRD_START_FAKE ( Fex_Fakeit_L_rec(1),
	1                  Time_buff_L(2), prev_gmt_hold_F_L, First )
	            EndIf
	         EndIf					! If Fakeit changes

C  Check if gain data changes. Set flags and increase counter for each channel.

	         Do Chan = 1,2
	            If ( Fex_Gain_R_Rec(2).Gain(chan)    .Ne.
	1              Fex_Gain_R_Rec(1).Gain(chan) )  Then
	               Fex_Gain_R_Rec(1).Gain_Change(chan) = 1
	               Gain_chan_ctr(chan) = Gain_chan_ctr(chan) + 1
	            EndIf
	            If ( Fex_Gain_L_Rec(2).Gain(chan)    .Ne.
	1              Fex_Gain_L_Rec(1).Gain(chan) )  Then
	               Fex_Gain_L_Rec(1).Gain_Change(chan) = 1
	               Gain_chan_ctr(chan+2) = Gain_chan_ctr(chan+2) + 1
	            EndIf
	         EndDo

C  If one or more Gain changed.

	         If (Fex_Gain_R_Rec(1).Gain_Change(1) +
	1           Fex_Gain_R_Rec(1).Gain_Change(2) .GT. 0) Then
	            Retstat = FRD_END_GAIN ( Fex_Gain_R_Rec(1), Time_buff_R(1),
	1              Lun_R_Gain )
	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Gain_Rec_R_Ctr = Gain_Rec_R_Ctr + 1
	               Prev_Gmt_Hold_G_R =Fex_Gain_R_Rec(1).CT_head.Gmt
	               Fex_Gain_R_rec(1) = Fex_Gain_R_rec(2)
	               status = FRD_START_GAIN ( Fex_Gain_R_rec(1),
	1                 Time_buff_R(2), prev_gmt_hold_G_R, First )
	            EndIf
	         EndIf				! If gain changed

	         If (Fex_Gain_L_Rec(1).Gain_Change(1) +
	1           Fex_Gain_L_Rec(1).Gain_Change(2) .GT. 0) Then
	            Retstat = FRD_END_GAIN ( Fex_Gain_L_Rec(1), Time_buff_L(1),
	1              Lun_L_Gain )
	            If (Retstat .NE. %loc(FRD_Normal)) Then
	               FRD_GFG_Check = Retstat
	            Else
	               Gain_Rec_L_Ctr = Gain_Rec_L_Ctr + 1
	               Prev_Gmt_Hold_G_L = Fex_Gain_L_Rec(1).CT_head.Gmt
	               Fex_Gain_L_rec(1) = Fex_Gain_L_rec(2)
	               status = FRD_START_GAIN ( Fex_Gain_L_rec(1),
	1                 Time_buff_L(2), prev_gmt_hold_G_L, First )
	            EndIf
	         EndIf				! If gain changed

	         Time_buff_R(1) = Time_buff_R(2)
	         Time_buff_L(1) = Time_buff_L(2)
	         Hkp_rec(1) = Hkp_rec(2)
	         If (Frm .Eq. 2) Then
	            Frm = 1
	         Else				! Frm = 1
	            Frm = 2
	         EndIf
	      EndIf				! If .Not.Next_Frm and status
	   EndIf				! read ok
	EndDo					! more_data

C  Write time of last frame found to report file

	If (Report) Then
	   If (Frm1_Rec .eq. 1) Then
	      Write (Lun_Rpt,30) Hkp_Rec(1).ct_head.gmt
	   Else
	      Write (Lun_Rpt,30) Hkp_Rec(1).hskp_tail.gmt_mjf2
	   End If
  30	   Format(1x,/,1x,'Time of Last Major Frame Found: ',A14)
	End If	      

C  Write record counters to report file

	If (Report) Write (lun_rpt,50) Fakeit_Rec_R_Ctr, Fakeit_Rec_L_Ctr,
	1   Gain_Rec_R_Ctr, Gain_Rec_L_Ctr, Gap_Fakeit_Ctr, Gap_Gain_Ctr,
	2   TQ_Ctr, Fakeit_Chan_Ctr, Gain_Chan_Ctr
  50	Format(1x,/,6x,'Number of FakeIt Records Written - Right Side:',I6,
	1      5x,'Left Side:',I6,/,
	1      6x,'Number of Gain Records Written - Right Side  :',I6,5x,
	1      'Left Side:',I6,/,
	1      6x,'Number of Housekeeping Gaps Causing FakeIt Records :',I6,/,
	1      6x,'Number of Housekeeping Gaps Causing Gain Records   :',I6,/,
	1      6x,'Number of Bad Telemetry Quality Reference Records  :',I6,/,
	1      6x,'Number of times Fakeit status changed :',/,25X,'RH:',I6,
	1      '   RL:',I6,'   LH:',I6,'   LL:',I6,/,
	1      6x,'Number of times Gain settings changed :',/,25X,'RH:',I6,
	1      '   RL:',I6,'   LH:',I6,'   LL:',I6)

	Return
	End
C******************************************************************************
	Integer*4 Function FRD_START_FAKE ( Fake_Rec, Time_buff, gmt, first )

C  Larry Rosen, STX, March 8 1991
C  Start a new write fake record

	Implicit None

	Character*14	gmt	! previous reference records gmt time.
	Logical*1	first	! true = first reference record, use self gmt.

	Dictionary 'FEX_FAKEIT'
	Record /FEX_FAKEIT/ Fake_Rec
	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ Time_Buff

	Integer*2	i
	Logical*1	TIME_LE

C  Find earliest time for CT_HEAD time of FEX reference records.

	Fake_Rec.CT_Head.GMT = Time_buff.FG
	Do i=1,2
	   Fake_Rec.CT_Head.Time(i) = Time_buff.FA(i)
	EndDo
	If (first) Then
	   Fake_Rec.Prev_GMT = Time_buff.FG
	Else
	   Fake_Rec.Prev_GMT = gmt
	EndIf
	Fake_Rec.Start_Gmt = Time_Buff.FG
	Do i=1,2
	    Fake_Rec.Start_Time(i)=Time_Buff.FA(i)
	EndDo

	Return
	End
C******************************************************************************
	Integer*4 Function FRD_START_GAIN ( Gain_Rec, Time_buff, gmt, first )

C  Larry Rosen, STX, March 8 1991
C  Start a new write gain record

	Implicit None

	Character*14	gmt	! previous reference records gmt time.
	Logical*1	first	! true = first reference record, use self gmt.

	Dictionary 'FEX_GAIN'
	Record /FEX_GAIN/ Gain_Rec
	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ Time_Buff

	Integer*2	i
	Logical*1	TIME_LE

C  Find earliest time for CT_HEAD time of FEX reference records.

	Gain_Rec.CT_Head.GMT = Time_buff.GG
	Do i=1,2
	   Gain_Rec.CT_Head.Time(i) = Time_buff.GA(i)
	EndDo
	If (first) Then
	   Gain_Rec.Prev_GMT = Time_buff.GG
	Else
	   Gain_Rec.Prev_GMT = gmt
	EndIf
	Gain_Rec.Start_Gmt = Time_Buff.GG
	Do i=1,2
	    Gain_Rec.Start_Time(i) = Time_Buff.GA(i)
	EndDo

	Return
	End
C******************************************************************************
	Integer*4 Function FRD_END_FAKE ( Fakeit_Rec, Time_buff, Lun )

C  Larry Rosen, STX, March 8 1991
C  Left and right independent, 18 June 1991
C  End and write fakeit record

	Implicit None

	Integer*4	LUN	! Unit number for FEX archive

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'

	External        FRD_Normal
	External        FRD_CTWritErr

	Dictionary 'FEX_FAKEIT'
	Record /FEX_FAKEIT/ Fakeit_Rec
	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ Time_Buff
	Integer*2	CT_Stat(20)	! Cobetrieve status
	Integer*4	status
	Integer*2	i

	FRD_END_FAKE = %loc(FRD_Normal)
	Fakeit_Rec.Stop_Gmt = Time_Buff.FG
	Do i=1,2
	    Fakeit_Rec.Stop_Time(i) = Time_Buff.FA(i)
	    Fakeit_Rec.Stop_Time(i) = Time_Buff.FA(i)
	EndDo
	Call CT_Write_arcv (,LUN, Fakeit_Rec, ct_stat)
	If ( CT_Stat(1) .NE. CTP_Normal ) Then
	    status = CT_Stat(1)
	    Call Lib$Signal (FRD_CTWRITERR, %val(2), %val(status),
	1       'FEX_Fakeit.Dat')
	    FRD_END_FAKE = %loc(FRD_CTWRITERR)
	EndIf
	Return
	End
C******************************************************************************
	Integer*4 Function FRD_END_GAIN ( Gain_Rec, Time_buff, Lun )

C  Larry Rosen, STX, March 8 1991
C  Left and right independent, 18 June 1991
C  End and write gain record

	Implicit None

	Integer*4	LUN	! Unit number for FEX archive

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'

	External        FRD_Normal
	External        FRD_CTWritErr

	Dictionary 'FEX_GAIN'
	Record /FEX_GAIN/ Gain_Rec
	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ Time_Buff
	Integer*2	CT_Stat(20)	! Cobetrieve status
	Integer*4	status
	Integer*2	i

	FRD_END_GAIN = %loc(FRD_Normal)
	Gain_Rec.Stop_Gmt = Time_Buff.GG
	Do i=1,2
	    Gain_Rec.Stop_Time(i) = Time_Buff.GA(i)
	    Gain_Rec.Stop_Time(i) = Time_Buff.GA(i)
	EndDo
	Call CT_Write_arcv (,LUN, Gain_Rec, ct_stat)
	If ( CT_Stat(1) .NE. CTP_Normal ) Then
	    status = CT_Stat(1)
	    Call Lib$Signal (FRD_CTWRITERR, %val(2), %val(status),
	1       'FEX_Gain.Dat')
	    FRD_END_GAIN = %loc(FRD_CTWRITERR)
	EndIf
	Return
	End
C******************************************************************************
	Integer*4 Function Fill_Time_Buff ( act_mj_a, act_mj_b, tbuff_R,
	1                                   tbuff_L, time_sec_1, time_sec_2)

C  Larry Rosen, STX, March 8, 1991
C  Based on actual frame number, put offset times into time buffer.
C  Get times of gain and fakeit data by converting gmt times to seconds,
C  adding offsets,  converting back to GMT and ADT times.
C  Gain fakeit times are frame gmt times plus offsets.
C   Take the major frame time, add 0.25 * (minor frame number - 4) seconds.
C  Fakeit statuses are in major frame 1, minor frame 12 for RH and RL, and
C  minor frame 14 for LH and LL.  This gives 2 times to start and 2 times
C  to end each reference record.
C   Gains are in major frame 2, minor frames 8, 12, and 16 for RH and RL,
C  and minor frames 10, 14, and 18 for LH and LL.  For gains, use the
C  first minor frame numbers (8 and 10) in this calculation.  This gives two
C  times to start and two times to end each reference record.
C   The time for a fakeit change is given by major frames 1 or 2 times (transmit
C  time in Ct_head) offset by the time into the frame that the fakeit status
C  change in the telemetry stream.

	Implicit None

        Integer*4       Act_Mj_A(2)     ! Actual a side major frame
        Integer*4       Act_mj_B(2)     ! Actual b side major frame
	Structure       /Time_Buffer/    ! Buffer to hold gain and fakeit times.
	  Integer*4     GA(2), FA(2)     ! Gain/Fakeit  ADT times.
	  Character*14  GG, FG           ! Gain/Fakeit  GMT times.
	End Structure
	Record /Time_Buffer/ TBuff_R, TBuff_L
	Real*8	time_sec_2, time_sec_1, time_sec
	Integer*4	status
      	Include 'CT$Library:CTUser.Inc'
	Integer*4	CT_GMT_to_Binary
	Integer*4	FUT_GMT_CONV

	If (act_mj_a(1) .eq. 1) Then
	    time_sec = time_sec_1 + 2
	    status = FUT_GMT_CONV(time_sec,TBUFF_R.FG,1)
	    Call CT_GMT_To_Binary (TBUFF_R.FG,TBUFF_R.FA)
	Else
	    time_sec = time_sec_1 + 1
	    status = FUT_GMT_CONV(time_sec,TBUFF_R.GG,1)
	    Call CT_GMT_To_Binary (TBUFF_R.GG,TBUFF_R.GA)
	Endif
	If (act_mj_a(2) .eq. 1) Then
	    time_sec = time_sec_2 + 2
	    status = FUT_GMT_CONV(time_sec,TBUFF_R.FG,1)
	    Call CT_GMT_To_Binary (TBUFF_R.FG,TBUFF_R.FA)
	Else
	    time_sec = time_sec_2 + 1
	    status = FUT_GMT_CONV(time_sec,TBUFF_R.GG,1)
	    Call CT_GMT_To_Binary (TBUFF_R.GG,TBUFF_R.GA)
	Endif
	If (act_mj_b(1) .eq. 1) Then
	    time_sec = time_sec_1 + 2.5
	    status = FUT_GMT_CONV(time_sec,TBUFF_L.FG,1)
	    Call CT_GMT_To_Binary (TBUFF_L.FG,TBUFF_L.FA)
	Else
	    time_sec = time_sec_1 + 1.5
	    status = FUT_GMT_CONV(time_sec,TBUFF_L.GG,1)
	    Call CT_GMT_To_Binary (TBUFF_L.GG,TBUFF_L.GA)
	Endif
	If (act_mj_b(2) .eq. 1) Then
	    time_sec = time_sec_2 + 2.5
	    status = FUT_GMT_CONV(time_sec,TBUFF_L.FG,1)
	    Call CT_GMT_To_Binary (TBUFF_L.FG,TBUFF_L.FA)
	Else
	    time_sec = time_sec_2 + 1.5
	    status = FUT_GMT_CONV(time_sec,TBUFF_L.GG,1)
	    Call CT_GMT_To_Binary (TBUFF_L.GG,TBUFF_L.GA)
	Endif
	Return
	End
