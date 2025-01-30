	INTEGER*4 FUNCTION FUT_QUAL_SUMMARY ( CHANNEL, QUALFLAGS, 
	1	END_SEGMENT, REPORT_LUN )
C/
C/ PROGRAM NAME: 
C/	FUT_QUAL_SUMMARY
C/
C/ PROGRAM DESCRIPTION:
C/	This routine accumulates counters for the data quality flags in the 
C/	channel science records. At the end of a segment the summary is printed.
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  January 1988.
C/
C/ MODIFIED BY:
C/	Shirley M. Read ( STX )  May 1988. 
C/	  Added report unit number to calling sequence.
C/	  Since FDQ_Qualflag_Names.Txt was moved to the FUT Facility, the
C/	  include file is now FUT_Qualflag_Names.
C/
C/ MODIFIED BY:
C/          H. Wang, STX, 1/29/91
C/          New requirements for FDQ
C/          Summary reports
C/
CH CHANGE LOG:			New Format for Build 4.2  STX, 10/15/88
CH
CH      Version 4.1.1 10/20/88, SPR 2667, Shirley M. Read, STX
CH		FDQ has an error in counting the number of quality flags set
CH		for exceeding limits or encounters other bad data. The bit
CH		flags can have a negative value due to the msb being set.
CH	Version 4.2.1 02/09/89, SPR 3062, Shirley M. Read, STX
CH		Moon contaminated IFGs failed by data quality checking. FES and
CH		FCI will interpret the data as being unusable. It is needed for
CH		moon modeling. The moon flag bit will not be counted in the
CH		data quality summary flag.
CH	Version 4.2.1 02/08/89, SPR 3297, Shirley M. Read, STX
CH		An attitude data quality summary flag should maintained apart
CH		from the existing over-all data quality summary flag. It will
CH		be assigned, just as the existing summary flag, a quality
CH		level dependent on the number of yellow and red limit 
CH		violations of the attitude quantities checked. FDQ will fill
CH		the attitude summary flag.
CH	Version 4.5 08/30/89, SER 4087, Shirley M. Read, STX
CH              FDQ needs additional limit check failure information for the
CH		Quicklook Analysis. The bit quality flags need to be 
CH		decommutated for verification of which engineering analog or
CH		which attitude limit or which microprocessor function caused
CH		the failure. This information needs to be available in the 
CH		report shortly after realtime passes as a check on the health
CH	        and safety of the instrument.
CH
C-------------------------------------------------------------------------------
C/
C/ CALLING SEQUENCE:
C/	Status = FUT_Qual_Summary( Channel, Qualflags, End_Segment, Report_Lun )
C/
C/ INPUT PARAMETERS:
C/	Channel     -- Current channel being processed
C/	Qualflags   -- 110 quality flags of current science record 
C/	End_Segment -- Flag indicating a segment has ended
C/	Report_Lun  -- Unit number for report
C/
C/ OUTPUT PARAMETERS:
C/
C/ ENTRY/EXIT:
C/	Normal function entry and return.
C/
C/ ERROR HANDLING
C/	Calls to Lib$Signal and interface with Fut_Error condition handler.
C/	
C/ INPUT/OUTPUT DISK FILES: 
C/	None
C/
C/ PROCS/MACROS/INCLUDE FILES:
C/	FUT_QUALFLAGS
C/	FUT_QUALFLAG_NAMES
C/
C/ SUBROUTINES CALLED:
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED: 
C/	This routine will clear output arrays for the report on the first
C/	call and start the accumulation of quality flag counts for the segment.
C/	On subsequent calls the flag counters will be accumulated. Only if
C/	it is a new segment will the report be written and the arrays reset.
C/
C/ PDL for FUT_QUAL_SUMMARY:
C/ 
C/	BEGIN
C/
C/	Set the status to success.
C/	If not the first time and end_segment flag is set, print the summary.
C/	Do for each of the non-summary quality flags in the channel array
C/	  If the quality flag is GT zero, increment the quality flag counter.
C/	Enddo
C/	Do for each of the bit quality flags in the channel array
C/        If the bit is set, increment the bit flag counter.
C/      Enddo
C/	If the quality summary flag is GT one, increment the summary quality
C/	   counter.
C/	Set function value for return.
C/
C/	RETURN
C/	END
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE
	
	EXTERNAL   	FUT_NORMAL
	EXTERNAL	FUT_ABERR
	Include 	'(FUT_PARAMS)'
	Include         '(FUT_QUALFLAGS)'
	include	        '(FUT_QUALFLAG_NAMES)'

!	Input/Output Passed Parameters

	integer*2	CHANNEL		! Micro-processor channel
	byte		QUALFLAGS(110)	! Science record data quality flags
	logical*1       END_SEGMENT     ! Segment ended
	integer*2	REPORT_LUN      ! Unit number for report

!	Local Variables
	
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable

	Integer*4	Qualflag_Counter(110,4)  ! Counters for 110 qual flags
						 ! for 4 channels
	logical*1       Write_Rpt /.false./, Decom /.false./ 
	logical*1 	First /.true./		 ! First time
	logical*1       P_Ipv /.false./, P_Ipc /.false./, 
	1	P_Micro /.false./, P_Grt_Red /.false./, 
	2	P_Grt_Yel /.false./, P_Mj_St /.false./,yel_tal
	INTEGER*4       CHANI4, FLAG
	integer*2 	RED, YELLOW     ! Counters for red and yellow flags
	byte            ZERO	        ! Zero value
	parameter       ( ZERO = 0 )
	byte		ONE
	parameter       ( ONE = 1 )	! Good quality flag
	byte	        two / 2 /       ! Some yellow, no red
	byte		three / 3 /     ! Many yellow, no red
	byte	        four / 4 /      ! No yellow, some red
	byte	        five / 5 /      ! Some yellow, some red
	byte	        six / 6 /	! Many yellow, some red
	byte	        seven / 7 /	! Many red
	byte		TLM_qual_failure  / 128 / ! TLM quality failure
	integer*2	Ix, Jx, Kx, Beg, Fin
	character*15 Micro_Name(0:15) / 	! Microprocessor status bits
	1 '  Grp Div Ovflo','  Grp Coad Ovfl','  Samp Div Ovfl',
	2 '  Samp Add Ovfl','  ADC Buf Ovflo','  Sync Error   ',
	3 '  Ill Instr Fet','  Ill Contex Sw','  ADC Pul Tmout',	
	4 '  Syn Pul Tmout','  Data Buf Clr ','  Cmd Overrun  ',
	5 '  Buffering    ','  Degl Math Ovf','               ',
	6 '               ' /
	integer*4 Micro_Flag(0:15,4) 
        Integer*4 Ipdu_TMP_A_RED_flag(4) 
        Integer*4 Ipdu_TMP_A_yel_flag(4)
        Integer*4 Ipdu_TMP_B_RED_flag(4) 
        Integer*4 Ipdu_TMP_B_Yel_flag(4) 
        Integer*4 CHAN_TMP_RED_flag(4) 
        Integer*4 CHAN_TMP_Yel_flag(4) 
        Integer*4 DBOX_TMP_A_RED_flag(4) 
        Integer*4 DBOX_TMP_A_yel_flag(4)
        Integer*4 DBOX_TMP_B_RED_flag(4) 
        Integer*4 DBOX_TMP_B_Yel_flag(4) 
        Integer*4 STMON_TMP_RED_flag(4) 
        Integer*4 STMON_TMP_Yel_flag(4) 
        Integer*4 CHPA_TMP_RED_flag(4) 
        Integer*4 CHPA_TMP_Yel_flag(4) 
        Integer*4 OPPA_TMP_RED_flag(4) 
        Integer*4 OPPA_TMP_Yel_flag(4) 
        Integer*4 HS_A_TMP_RED_flag(4) 
        Integer*4 HS_A_TMP_Yel_flag(4) 
        Integer*4 HS_B_TMP_RED_flag(4) 
        Integer*4 HS_B_TMP_Yel_flag(4) 
        Integer*4 MTM_A_RED_flag(4) 
        Integer*4 MTM_A_Yel_flag(4) 
        Integer*4 MTM_B_RED_flag(4) 
        Integer*4 MTM_B_Yel_flag(4) 
        Integer*4 BOL_RED_flag(4) 
        Integer*4 BOL_Yel_flag(4) 
        Integer*4 IPDU_vol_RED_flag(10,4) 
        Integer*4 IPDU_VOL_Yel_flag(10,4) 
        Integer*4 IPDU_cur_RED_flag(6,4) 
        Integer*4 IPDU_cur_Yel_flag(6,4) 
        
	character*15 Att_Name(0:7) / 	! Attitude limits
	1 '  Sun Angle    ','  Earth Limb   ','  Moon Angle   ',
	2 '  SAA          ','  VAB          ','               ',
	3 '               ','               ' /
	integer*4 Red_Att_Flag(0:7,4), Yel_Att_Flag(0:7,4)

C 	GRTs Temperatures: Side A, Low Current; Side A, High Current;
C			   Side B, Low Current; Side B, High Current.

	character*15 GRT_Name(0:63) /  
	1 '  Xcal_Temp_Al ','  SkyHorn_Tp_Al','  RefHorn_Tp_Al',
	2 '  Ical_Temp_Al ','  Dihedrl_Tp_Al','  Bol_A_RH_T_Al',
	3 '  Bol_A_RL_T_Al','  Bol_A_LH_T_Al','  Bol_A_LL_T_Al',
	4 '  Mirror_Tmp_Al','  Calrs_RH_T_Al','  Calrs_RL_T_Al',
	5 '  Calrs_LH_T_Al','  Calrs_LL_T_Al','  Xcal_Tp_S5_Al',
	6 '  Colimatr_T_Al',  			
	1 '  Xcal_Temp_Ah ','  SkyHorn_Tp_Ah','  RefHorn_Tp_Ah', 
	2 '  Ical_Temp_Ah ','  Dihedrl_Tp_Ah','  Bol_A_RH_T_Ah',
	3 '  Bol_A_RL_T_Ah','  Bol_A_LH_T_Ah','  Bol_A_LL_T_Ah',
	4 '  Mirror_Tmp_Ah','  Calrs_RH_T_Ah','  Calrs_RL_T_Ah',
	5 '  Calrs_LH_T_Ah','  Calrs_LL_T_Ah','  Xcal_Tp_S5_Ah',
	6 '  Colimatr_T_Ah',   			
	1 '  Xcal_Temp_Bl ','  SkyHorn_Tp_Bl','  RefHorn_Tp_Bl',
	2 '  Ical_Temp_Bl ','  Dihedrl_Tp_Bl','  Bol_A_RH_T_Bl',
	3 '  Bol_A_RL_T_Bl','  Bol_A_LH_T_Bl','  Bol_A_LL_T_Bl',
	4 '  Mirror_Tmp_Bl','  Calrs_RH_T_Bl','  Calrs_RL_T_Bl',
	5 '  Calrs_LH_T_Bl','  Calrs_LL_T_Bl','  Xcal_Tp_S6_Bl',
	6 '  Colimatr_T_Bl',   			
	1 '  Xcal_Temp_Bh ','  SkyHorn_Tp_Bh','  RefHorn_Tp_Bh',
	2 '  Ical_Temp_Bh ','  Dihedrl_Tp_Bh','  Bol_A_RH_T_Bh',
	3 '  Bol_A_RL_T_Bh','  Bol_A_LH_T_Bh','  Bol_A_LL_T_Bh',
	4 '  Mirror_Tmp_Bh','  Calrs_RH_T_Bh','  Calrs_RL_T_Bh',
	5 '  Calrs_LH_T_Bh','  Calrs_LL_T_Bh','  Xcal_Tp_S5_Bh',
	6 '  Colimatr_T_Bh' /  			

	integer*4 Red_GRT_Flag(0:63,4)
	integer*4 Yel_GRT_Flag(0:63,4)

!       Temperature Controllers are not converted as of this Build 4.5.
!       Thus no limits are checked.

	character*15 Ipdu_Volt_Name(10) /        ! IPDU Voltages
	1 '  Dig_Cnvt_N15v','  Dig_Cnvt_P15v','  Dig_Cnvt_P5v ',
	2 '  Anl_Cnvt_P15v','  Anl_Cnvt_N15v','  Bias_Pre_P25v',
	3 '  Int_Ps_P28v  ','  Int_Ps_P15v  ','  Int_Ps_N15v  ',
	4 '  Int_Ps_P5v   ' /

	character*15 Ipdu_Cur_Name(6) /		 ! IPDU Currents.
	1 '  Bias_Pre_Reg ','  Analog_Conv  ','  Digital_Conv ',
	2 '  Con_Current_H','  Con_Current_L','  Con_Int_Conv ' /

!	Limits are not checked in Build 4.5 for crosstalk with the spacecraft
!	or other instruments.
	character*15 DIRBE_Name(8) /		 ! DIRBE limits to be checked
	1 '  Anneal Htr1 Y','  Anneal Htr2 Y','  Anneal Htr1 R',
	2 '  Anneal Htr2 R','  IRS Sources Y','  IRS Sources R',
	3 '  Chopper      ','  Shutter      ' /
	character*15 DMR_Name(6) /               ! DMR limits to be checked
	1 '  Local Oscil Y','  Local Oscil R','  Noise Srce1 Y',
	2 '  Noise Srce2 Y','  Noise Srce1 R','  Noise Srce2 R' /
	character*15 SC_Name(8) /
	1 '  Torquer Bar Y','  Torquer Bar R','  Momentum Wh Y',
	2 '  Momentum Wh R','               ','               ',
	3 '               ','               ' /

	character*15 MJ_ST_Name(0:31) /    ! Status changes across the major frame
	1 '  GRT Dwell St ','  Micro Bus St ','  Bolo Bias St ',
	2 '  Temp Cont St1','  Temp Cont St2','  Temp Cont St3',
	3 '  Temp Cont St4','  Xcal Position','  Fakeit Bit St',
	4 '  Micro Sci Mod','  Micro Dat Typ','  Xcal Latch St',
	5 '  MTM Latch St ','  LVDT Status  ','  MTM Speed St ',
	6 '  MTM Length St','               ','               ',
	7 '               ','               ','               ',
	8 '               ','               ','               ',
	9 '               ','               ','               ',
	9 '               ','               ','               ',
	9 '               ','               ' /
	integer*4  Mj_St_Flag(0:31,4)
  
!	integer*2       Attitude_Flag   ! Temporary attitude flag value.
!	integer*2       Moon_Mask / '00FB'x /  ! Mask out moon bit in attitude
	byte 		Qual(2)
        integer*2       I2_Qual
        Integer*4       LX          
	equivalence	(Qual(1), I2_Qual)

!	Set return status to success.

	RETSTAT = SUCCESS

!       First Time: Initialize
 
	If ( First ) then
	  DO Ix = 1, 4
           Ipdu_TMP_A_RED_flag(Ix) = 0
           Ipdu_TMP_A_yel_flag(Ix)= 0
           Ipdu_TMP_B_RED_flag(IX) = 0
           Ipdu_TMP_B_Yel_flag(IX) = 0
           CHAN_TMP_RED_flag(IX) = 0
           CHAN_TMP_Yel_flag(IX) = 0
           DBOX_TMP_A_RED_flag(IX) = 0
           DBOX_TMP_A_yel_flag(IX)= 0
           DBOX_TMP_B_RED_flag(IX) = 0
           DBOX_TMP_B_Yel_flag(IX)= 0 
           STMON_TMP_RED_flag(IX) = 0
           STMON_TMP_Yel_flag(IX) = 0
           CHPA_TMP_RED_flag(IX) = 0
           CHPA_TMP_Yel_flag(IX) = 0
           OPPA_TMP_RED_flag(IX) = 0
           OPPA_TMP_Yel_flag(IX) = 0
           HS_A_TMP_RED_flag(IX) = 0
           HS_A_TMP_Yel_flag(IX) = 0
           HS_B_TMP_RED_flag(IX) = 0
           HS_B_TMP_Yel_flag(IX) = 0
           MTM_A_RED_flag(IX) = 0
           MTM_A_Yel_flag(IX) = 0
           MTM_B_RED_flag(IX) = 0
           MTM_B_Yel_flag(IX) = 0
           BOL_RED_flag(IX) = 0
           BOL_Yel_flag(IX) = 0
           DO KX = 1, 10
             IPDU_vol_RED_flag(KX,IX) = 0
             IPDU_VOL_Yel_flag(KX,IX) = 0
           ENDDO
           Do kx=1,6 
             IPDU_cur_RED_flag(kx,IX) =0
             IPDU_cur_Yel_flag(KX,IX) =0
           ENDDO    
	    DO Jx = 1, 110
	      Qualflag_Counter(Jx,Ix) = 0
	    ENDDO
	  ENDDO
	Endif
!
!	If a new segment is starting, print the summary and re-initialize.
!	Do not print on the first time.

	If ( End_Segment ) then
	   Write_Rpt = .true.
	Else
	   Write_Rpt = .false.
	Endif
	If ( First ) then
	   First = .false.
	   Write_Rpt = .false.
	   End_Segment = .false.
	Endif

	If ( Write_Rpt ) then
	  Write (Report_Lun, 100) 
 100	  Format (//1x,'Data quality summary:'//
	1   2x,'Quality Flag                   Number of Times Flag Set for '/,
	2   17x,'              RH          RL          LH         LL '/)
          Write (report_lun,149) ' Red Instrument '
	  DO Ix = 1, 85
	    If (( Qualflag_Counter(Ix,1) .GT. Zero) .OR.
	1       ( Qualflag_Counter(Ix,2) .GT. Zero) .OR.
	2       ( Qualflag_Counter(Ix,3) .GT. Zero) .OR.
	3       ( Qualflag_Counter(Ix,4) .GT. Zero)) Then
              If ((ix .ge. 1 .and. ix .le. 3) .or. (ix .ge. 8 .and.
	1	ix .le. 10)) then
                write(report_lun,200) qualflag_names(IX), 
	1	(Qualflag_Counter(Ix,Jx),Jx=1,4)
              endif
            If (ix .ge. 4 .and. ix .le. 5) then
	      If ( .not. P_Micro) Then
		Write (Report_Lun,150) Qualflag_Names(4)
	        Do Kx = 0, 15
	          If (( Micro_Flag(Kx,1) .gt. Zero ) .or.
	1	       ( Micro_Flag(Kx,2) .gt. Zero ) .or.
	2	       ( Micro_Flag(Kx,3) .gt. Zero ) .or.
	3	       ( Micro_Flag(Kx,4) .gt. Zero )) Then
	              Write (Report_Lun,200) Micro_Name(Kx), 
	1	      (Micro_Flag(Kx,Jx),Jx=1,4)
	            Endif
	        Enddo
	        P_Micro = .true.
	     Endif
	    ENDIF
            If (ix .ge. 11 .and. ix .le. 18) then
	      If ( .not. P_Grt_Red) Then
		Write (Report_Lun,150) Qualflag_Names(11)
	        Do Kx = 0, 63
	          If (( Red_GRT_Flag(Kx,1) .gt. Zero ) .or.
	1	      ( Red_GRT_Flag(Kx,2) .gt. Zero ) .or.
	2	      ( Red_GRT_Flag(Kx,3) .gt. Zero ) .or.
	3	      ( Red_GRT_Flag(Kx,4) .gt. Zero )) Then
	             Write (Report_Lun,200) GRT_Name(Kx), 
	1	      (Red_GRT_Flag(Kx,Jx),Jx=1,4)
		  Endif
	        Enddo
		P_Grt_Red = .true.
              Endif
	    ENDIF
             If ((ix .eq. 29) .and. 
	1	((IPDU_TMP_A_RED_FLAG(1) .ne. 0) .or.
	1	(IPDU_TMP_A_RED_FLAG(2) .ne. 0) .or.
	1	(IPDU_TMP_A_RED_FLAG(3) .ne. 0) .or.
	1	(IPDU_TMP_A_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(IPDU_TMP_A_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 30) .and. 
	1	((IPDU_TMP_B_RED_FLAG(1) .ne. 0) .or.
	1	(IPDU_TMP_B_RED_FLAG(2) .ne. 0) .or.
	1	(IPDU_TMP_B_RED_FLAG(3) .ne. 0) .or.
	1	(IPDU_TMP_B_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(IPDU_TMP_B_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 31) .and. 
	1	((CHAN_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(CHAN_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(CHAN_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(CHAN_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(CHAN_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 32) .and. 
	1	((DBOX_TMP_A_RED_FLAG(1) .ne. 0) .or.
	1	(DBOX_TMP_A_RED_FLAG(2) .ne. 0) .or.
	1	(DBOX_TMP_A_RED_FLAG(3) .ne. 0) .or.
	1	(DBOX_TMP_A_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(DBOX_TMP_A_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 33) .and. 
	1	((DBOX_TMP_B_RED_FLAG(1) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(2) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(3) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX-1), 
	1	(DBOX_TMP_B_RED_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 34) .and. 
	1	((STMON_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(STMON_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(STMON_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(STMON_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(STMON_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 35) .and. 
	1	((CHPA_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(CHPA_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(CHPA_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(CHPA_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(CHPA_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 36) .and. 
	1	((OPPA_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(OPPA_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(OPPA_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(OPPA_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(OPPA_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 37) .and. 
	1	((HS_A_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(HS_A_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(HS_A_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(HS_A_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(HS_A_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 38) .and. 
	1	((HS_B_TMP_RED_FLAG(1) .ne. 0) .or.
	1	(HS_B_TMP_RED_FLAG(2) .ne. 0) .or.
	1	(HS_B_TMP_RED_FLAG(3) .ne. 0) .or.
	1	(HS_B_TMP_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX-1), 
	1	(HS_B_TMP_RED_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 39) .and. 
	1	((MTM_A_RED_FLAG(1) .ne. 0) .or.
	1	(MTM_A_RED_FLAG(2) .ne. 0) .or.
	1	(MTM_A_RED_FLAG(3) .ne. 0) .or.
	1	(MTM_A_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(MTM_A_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 40) .and. 
	1	((MTM_B_RED_FLAG(1) .ne. 0) .or.
	1	(MTM_B_RED_FLAG(2) .ne. 0) .or.
	1	(MTM_B_RED_FLAG(3) .ne. 0) .or.
	1	(MTM_B_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX-1), 
	1	(MTM_B_RED_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 41) .and. 
	1	((BOL_RED_FLAG(1) .ne. 0) .or.
	1	(bol_RED_FLAG(2) .ne. 0) .or.
	1	(bol_RED_FLAG(3) .ne. 0) .or.
	1	(bol_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX), 
	1	(bol_RED_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 33) .and. 
	1	((DBOX_TMP_B_RED_FLAG(1) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(2) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(3) .ne. 0) .or.
	1	(DBOX_TMP_B_RED_FLAG(4) .ne. 0))) then
                write(report_lun,200) qualflag_names(IX-1), 
	1	(DBOX_TMP_B_RED_FLAG(Jx),Jx=1,4)
              endif
	      If ((Ix .ge. 42) .and. (Ix .le. 51)) Then
	        If ( .not. P_Ipv) Then
		   P_Ipv = .true.
	  	   Write (Report_Lun,150) Qualflag_Names(42)
 149		   Format(1x,a)
	        Endif
	        Kx = Ix - 41
                IF ( (IPDU_VOL_RED_FLAG(KX,1) .ne. 0) .or.
	1	 (IPDU_VOL_RED_FLAG(KX,2) .ne. 0) .or.
	1	 (IPDU_VOL_RED_FLAG(KX,3) .ne. 0) .or.
	1	 (IPDU_VOL_RED_FLAG(KX,4) .ne. 0)) then
	          Write (Report_Lun,200) IPDU_Volt_Name(Kx), 
	1	  (IPDU_VOL_RED_FLAG(KX,Jx),Jx=1,4)
                ENDIF
              ENDIF
	      if ((Ix .ge. 52) .and. (Ix .le. 57)) Then
	        If (.not. P_Ipc) Then
		   P_Ipc = .true.
	  	   Write (Report_Lun,150) Qualflag_Names(52)
	        Endif
	        Kx = Ix - 51
                IF ( (IPDU_CUR_RED_FLAG(KX,1) .ne. 0) .or.
	1	 (IPDU_CUR_RED_FLAG(KX,2) .ne. 0) .or.
	1	 (IPDU_CUR_RED_FLAG(KX,3) .ne. 0) .or.
	1	 (IPDU_CUR_RED_FLAG(KX,4) .ne. 0)) then
         	  Write (Report_Lun,200) IPDU_Cur_Name(Kx), 
	1	   (IPDU_CUR_RED_FLAG(kx,Jx),Jx=1,4)
                endif
	      ENDIF
	    IF(IX .ge. 58 .and. IX .le. 61) then
		If ( .not. P_Mj_St) Then
		  Write (Report_Lun,150) Qualflag_Names(58)
		  Do Kx = 0, 31
		    If ((Mj_St_Flag(Kx,1) .Gt. Zero) .or.
	1	        (Mj_St_Flag(Kx,2) .Gt. Zero) .or.
	2               (Mj_St_Flag(Kx,3) .Gt. Zero) .or.
	3	        (Mj_St_Flag(Kx,4) .Gt. Zero)) Then
	              Write (Report_Lun,200) Mj_St_Name(Kx), 
	1	      (Mj_St_Flag(Kx,Jx),Jx=1,4)
	            Endif
	          Enddo
		  P_Mj_St = .true.
	        Endif	
	    ENDIF
           ENDIF
	  ENDDO	
	  P_Micro = .false.
	  P_Grt_Red = .false.
	  P_Grt_Yel = .false.
	  P_Mj_St = .false.
	  P_Ipv = .false.
	  P_Ipc = .false.
          yel_tal=.false.
	  DO Ix = 1, 85
	    If (( Qualflag_Counter(Ix,1) .GT. Zero) .OR.
	1       ( Qualflag_Counter(Ix,2) .GT. Zero) .OR.
	2       ( Qualflag_Counter(Ix,3) .GT. Zero) .OR.
	3       ( Qualflag_Counter(Ix,4) .GT. Zero)) Then
              If (ix .ge. 19 .and. ix .le. 26) then
                If(.not. yel_tal) then
                  Write (report_lun,149) ' Yellow Instrument '
                  yel_tal = .true.
                Endif
	        If ( .not. P_Grt_Yel) Then
		  Write (Report_Lun,150) Qualflag_Names(19)
	          Do Kx = 0, 63
	            If (( Yel_GRT_Flag(Kx,1) .gt. Zero ) .or.
	1	      ( Yel_GRT_Flag(Kx,2) .gt. Zero ) .or.
	2	      ( Yel_GRT_Flag(Kx,3) .gt. Zero ) .or.
	3	      ( Yel_GRT_Flag(Kx,4) .gt. Zero )) Then
	              Write (Report_Lun,200) GRT_Name(Kx), 
	1	      (Yel_GRT_Flag(Kx,Jx),Jx=1,4)
		    Endif
	          Enddo
		  P_Grt_Yel = .true.
	        Endif
	    ENDIF
          If ((ix .eq. 29) .and. 
	1	((IPDU_TMP_A_YEL_FLAG(1) .ne. 0) .or.
	1	(IPDU_TMP_A_YEL_FLAG(2) .ne. 0) .or.
	1	(IPDU_TMP_A_YEL_FLAG(3) .ne. 0) .or.
	1	(IPDU_TMP_A_YEL_FLAG(4) .ne. 0))) then
                  If(.not. yel_tal) then
                   Write (report_lun,149)  ' Yellow Instrument '
                   yel_tal = .true.
                 Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(IPDU_TMP_A_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 30) .and. 
	1	((IPDU_TMP_B_YEL_FLAG(1) .ne. 0) .or.
	1	(IPDU_TMP_B_YEL_FLAG(2) .ne. 0) .or.
	1	(IPDU_TMP_B_YEL_FLAG(3) .ne. 0) .or.
	1	(IPDU_TMP_B_YEL_FLAG(4) .ne. 0))) then
                 If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                 Endif
                 write(report_lun,200) qualflag_names(IX), 
	1	(IPDU_TMP_B_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 31) .and. 
	1	((CHAN_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(CHAN_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(CHAN_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(CHAN_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(CHAN_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 32) .and. 
	1	((DBOX_TMP_A_YEL_FLAG(1) .ne. 0) .or.
	1	(DBOX_TMP_A_YEL_FLAG(2) .ne. 0) .or.
	1	(DBOX_TMP_A_YEL_FLAG(3) .ne. 0) .or.
	1	(DBOX_TMP_A_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(DBOX_TMP_A_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 33) .and. 
	1	((DBOX_TMP_B_YEL_FLAG(1) .ne. 0) .or.
	1	(DBOX_TMP_B_YEL_FLAG(2) .ne. 0) .or.
	1	(DBOX_TMP_B_YEL_FLAG(3) .ne. 0) .or.
	1	(DBOX_TMP_B_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX-1), 
	1	(DBOX_TMP_B_YEL_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 34) .and. 
	1	((STMON_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(STMON_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(STMON_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(STMON_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(STMON_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 35) .and. 
	1	((CHPA_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(CHPA_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(CHPA_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(CHPA_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(CHPA_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 36) .and. 
	1	((OPPA_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(OPPA_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(OPPA_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(OPPA_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(OPPA_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 37) .and. 
	1	((HS_A_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(HS_A_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(HS_A_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(HS_A_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(HS_A_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 38) .and. 
	1	((HS_B_TMP_YEL_FLAG(1) .ne. 0) .or.
	1	(HS_B_TMP_YEL_FLAG(2) .ne. 0) .or.
	1	(HS_B_TMP_YEL_FLAG(3) .ne. 0) .or.
	1	(HS_B_TMP_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149)  ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX-1), 
	1	(HS_B_TMP_YEL_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 39) .and. 
	1	((MTM_A_YEL_FLAG(1) .ne. 0) .or.
	1	(MTM_A_YEL_FLAG(2) .ne. 0) .or.
	1	(MTM_A_YEL_FLAG(3) .ne. 0) .or.
	1	(MTM_A_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149) ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(MTM_A_YEL_FLAG(Jx),Jx=1,4)
              endif
              If ((ix .eq. 40) .and. 
	1	((MTM_B_YEL_FLAG(1) .ne. 0) .or.
	1	(MTM_B_YEL_FLAG(2) .ne. 0) .or.
	1	(MTM_B_YEL_FLAG(3) .ne. 0) .or.
	1	(MTM_B_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149) ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX-1), 
	1	(MTM_B_YEL_FLAG(Jx),Jx=1,4)
              endif
          If ((ix .eq. 41) .and. 
	1	((BOL_YEL_FLAG(1) .ne. 0) .or.
	1	(bol_YEL_FLAG(2) .ne. 0) .or.
	1	(bol_YEL_FLAG(3) .ne. 0) .or.
	1	(bol_YEL_FLAG(4) .ne. 0))) then
                If(.not. yel_tal) then
                  Write (report_lun,149) ' Yellow Instrument '
                  yel_tal = .true.
                Endif
                write(report_lun,200) qualflag_names(IX), 
	1	(bol_YEL_FLAG(Jx),Jx=1,4)
              endif
	      If ((Ix .ge. 42) .and. (Ix .le. 51)) Then
	        Kx = Ix - 41
                IF ( (IPDU_VOL_YEL_FLAG(KX,1) .ne. 0) .or.
	1	 (IPDU_VOL_YEL_FLAG(KX,2) .ne. 0) .or.
	1	 (IPDU_VOL_YEL_FLAG(KX,3) .ne. 0) .or.
	1	 (IPDU_VOL_YEL_FLAG(KX,4) .ne. 0)) then
                 If(.not. yel_tal) then
                  Write (report_lun,149) ' Yellow Instrument '
                  yel_tal = .true.
                 Endif
	         If ( .not. P_Ipv) Then
		   P_Ipv = .true.
	  	   Write (Report_Lun,150) Qualflag_Names(42)
	         Endif
	         Write (Report_Lun,200) IPDU_Volt_Name(Kx), 
	1	 (IPDU_VOL_YEL_FLAG(KX,Jx),Jx=1,4)
                ENDIF
              ENDIF
	      if ((Ix .ge. 52) .and. (Ix .le. 57)) Then
	        Kx = Ix - 51
                IF ( (IPDU_CUR_YEL_FLAG(KX,1) .ne. 0) .or.
	1	 (IPDU_CUR_YEL_FLAG(KX,2) .ne. 0) .or.
	1	 (IPDU_CUR_YEL_FLAG(KX,3) .ne. 0) .or.
	1	 (IPDU_CUR_YEL_FLAG(KX,4) .ne. 0)) then
                  If(.not. yel_tal) then
                    Write (report_lun,149) ' Yellow Instrument '
                    yel_tal = .true.
                  Endif
	          If (.not. P_Ipc) Then
		   P_Ipc = .true.
	  	   Write (Report_Lun,150) Qualflag_Names(52)
	          Endif
         	  Write (Report_Lun,200) IPDU_Cur_Name(Kx), 
	1	   (IPDU_CUR_YEL_FLAG(kx,Jx),Jx=1,4)
                endif
	      ENDIF
	    ENDIF
	  ENDDO	
150	  Format(5x,a)
200       Format (5x,a,2x,i10,2x,i10,2x,i10,2x,i10)
205       Format (/2x,a,5x,i10,2x,i10,2x,i10,2x,i10)
	    IF (( Qualflag_Counter(6,1) .GT. Zero) .OR.
	1       ( Qualflag_Counter(6,2) .GT. Zero) .OR.
	2       ( Qualflag_Counter(6,3) .GT. Zero) .OR.
	3       ( Qualflag_Counter(6,4) .GT. Zero)) THEN
                  Write (report_lun,149) ' Red Attitude '
		  Write (Report_Lun,150) Qualflag_Names(6)
	          Do Kx = 0, 7
	            If (( Red_Att_Flag(Kx,1) .gt. Zero ) .or.
	1	      ( Red_Att_Flag(Kx,2) .gt. Zero ) .or.
	2	      ( Red_Att_Flag(Kx,3) .gt. Zero ) .or.
	3	      ( Red_Att_Flag(Kx,4) .gt. Zero )) Then
	              Write (Report_Lun,200) Att_Name(Kx), 
	1	      (Red_Att_Flag(Kx,Jx),Jx=1,4)
		    Endif
	          Enddo
            ENDIF
	    IF (( Qualflag_Counter(7,1) .GT. Zero) .OR.
	1       ( Qualflag_Counter(7,2) .GT. Zero) .OR.
	2       ( Qualflag_Counter(7,3) .GT. Zero) .OR.
	3       ( Qualflag_Counter(7,4) .GT. Zero)) THEN
                  Write (report_lun,149) ' Yellow Attitude '
		  Write (Report_Lun,150) Qualflag_Names(7)
	          Do Kx = 0, 7
	            If (( Yel_Att_Flag(Kx,1) .gt. Zero ) .or.
	1	      ( Yel_Att_Flag(Kx,2) .gt. Zero ) .or.
	2	      ( Yel_Att_Flag(Kx,3) .gt. Zero ) .or.
	3	      ( Yel_Att_Flag(Kx,4) .gt. Zero )) Then
	              Write (Report_Lun,200) Att_Name(Kx), 
	1	      (Yel_Att_Flag(Kx,Jx),Jx=1,4)
		    Endif
	          Enddo
            ENDIF
	    IF (( Qualflag_Counter(108,1) .GT. Zero) .OR.
	1       ( Qualflag_Counter(108,2) .GT. Zero) .OR.
	2       ( Qualflag_Counter(108,3) .GT. Zero) .OR.
	3       ( Qualflag_Counter(108,4) .GT. Zero)) THEN
   	        Write (Report_Lun,205) Qualflag_Names(108), 
	1	(Qualflag_Counter(108,Jx),Jx=1,4)
            Endif
	    DO Ix = 109, 110
	      Write (Report_Lun,205) Qualflag_Names(Ix), 
	1	(Qualflag_Counter(Ix,Jx),Jx=1,4)
	    ENDDO	
	  DO Ix = 1, 4
           Ipdu_TMP_A_RED_flag(Ix) = 0
           Ipdu_TMP_A_yel_flag(Ix)= 0
           Ipdu_TMP_B_RED_flag(IX) = 0
           Ipdu_TMP_B_Yel_flag(IX) = 0
           CHAN_TMP_RED_flag(IX) = 0
           CHAN_TMP_Yel_flag(IX) = 0
           DBOX_TMP_A_RED_flag(IX) = 0
           DBOX_TMP_A_yel_flag(IX)= 0
           DBOX_TMP_B_RED_flag(IX) = 0
           DBOX_TMP_B_Yel_flag(IX)= 0 
           STMON_TMP_RED_flag(IX) = 0
           STMON_TMP_Yel_flag(IX) = 0
           CHPA_TMP_RED_flag(IX) = 0
           CHPA_TMP_Yel_flag(IX) = 0
           OPPA_TMP_RED_flag(IX) = 0
           OPPA_TMP_Yel_flag(IX) = 0
           HS_A_TMP_RED_flag(IX) = 0
           HS_A_TMP_Yel_flag(IX) = 0
           HS_B_TMP_RED_flag(IX) = 0
           HS_B_TMP_Yel_flag(IX) = 0
           MTM_A_RED_flag(IX) = 0
           MTM_A_Yel_flag(IX) = 0
           MTM_B_RED_flag(IX) = 0
           MTM_B_Yel_flag(IX) = 0
           BOL_RED_flag(IX) = 0
           BOL_Yel_flag(IX) = 0
           DO KX = 1, 10
             IPDU_vol_RED_flag(KX,IX) = 0
             IPDU_VOL_Yel_flag(KX,IX) = 0
           ENDDO
           Do kx=1,6 
             IPDU_cur_RED_flag(kx,IX) = 0
             IPDU_cur_Yel_flag(KX,IX) = 0
           ENDDO    
             
	    DO Jx = 1, 110
	      Qualflag_Counter(Jx,Ix) = 0
	    ENDDO
	    DO Jx = 0, 63
	      Red_Grt_Flag(Jx,Ix) = 0
	      Yel_Grt_Flag(Jx,Ix) = 0
	    ENDDO
	    DO Jx = 0, 31
	      Mj_St_Flag(Jx,Ix) = 0
	    ENDDO
	    DO Jx = 0, 15
	      Micro_Flag(Jx,Ix) = 0
	    ENDDO
	    DO Jx = 0, 7
	      Red_Att_Flag(Jx,Ix) = 0
	      Yel_Att_Flag(Jx,Ix) = 0
	    ENDDO
	  ENDDO

	  Write_Rpt = .false.
	  End_Segment = .false.
	  Decom = .false.
	  P_Micro = .false.
	  P_Grt_Red = .false.
	  P_Grt_Yel = .false.
	  P_Mj_St = .false.
	  P_Ipv = .false.
	  P_Ipc = .false.

	Endif                               ! Write report

!	Accumulate counters. 

	Do Ix = 1, 85
	    if ( qualflags(Ix) .ne. zero ) then
	      Qualflag_Counter(Ix,Channel) = 
	1	Qualflag_Counter(Ix,Channel) + 1
              if (ix .eq. 29) then
	        If (qualflags(ix) .eq. 2) then
		 IPDU_TMP_A_RED_FLAG(CHANNEL) = IPDU_TMP_A_RED_FLAG(CHANNEL) + 1
                Else
                 If (qualflags(ix) .eq. 1) 
	1	 IPDU_TMP_A_YEL_FLAG(CHANNEL) = IPDU_TMP_A_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 30) then
                if (qualflags(ix) .eq. 2) then 
		 IPDU_TMP_B_RED_FLAG(CHANNEL) = IPDU_TMP_B_RED_FLAG(CHANNEL) + 1
                Else
                 if (qualflags(ix) .eq. 1) 
	1	IPDU_TMP_B_YEL_FLAG(CHANNEL) = IPDU_TMP_B_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 31) then
		if (qualflags(ix) .eq. 2) then 
		 CHAN_TMP_RED_FLAG(CHANNEL) = CHAN_TMP_RED_FLAG(CHANNEL) + 1
                else
		 if (qualflags(ix) .eq. 1) 
	1	CHAN_TMP_YEL_FLAG(CHANNEL) = CHAN_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 32) then
		If (qualflags(ix) .eq. 2) then 
	  	 DBOX_TMP_A_RED_FLAG(CHANNEL) = DBOX_TMP_A_RED_FLAG(CHANNEL) + 1
                else 
		  if(qualflags(ix) .eq. 1) 
	1	 DBOX_TMP_A_YEL_FLAG(CHANNEL) = DBOX_TMP_A_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 33) then
		if (qualflags(ix) .eq. 2) then 
		 DBOX_TMP_B_RED_FLAG(CHANNEL) = DBOX_TMP_B_RED_FLAG(CHANNEL) + 1
                else 
		  If(qualflags(ix) .eq. 1) 
	1	 DBOX_TMP_B_YEL_FLAG(CHANNEL) = DBOX_TMP_B_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 34) then
		If (qualflags(ix) .eq. 2) then 
		 STMON_TMP_RED_FLAG(CHANNEL) = STMON_TMP_RED_FLAG(CHANNEL) + 1
                else 
                 If(qualflags(ix) .eq. 1) 
	1	  STMON_TMP_YEL_FLAG(CHANNEL) = STMON_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 35) then
		If (qualflags(ix) .eq. 2) then 
		 CHPA_TMP_RED_FLAG(CHANNEL) = CHPA_TMP_RED_FLAG(CHANNEL) + 1
                else 
		 if(qualflags(ix) .eq. 1)  
	1	  CHPA_TMP_YEL_FLAG(CHANNEL) = CHPA_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif 
              if (ix .eq. 36) then
		if (qualflags(ix) .eq. 2) then 
		 OPPA_TMP_RED_FLAG(CHANNEL) = OPPA_TMP_RED_FLAG(CHANNEL) + 1
                else
		  if (qualflags(ix) .eq. 1) 
	1	   OPPA_TMP_YEL_FLAG(CHANNEL) = OPPA_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 37) then
		if (qualflags(ix) .eq. 2) then 
		 HS_A_TMP_RED_FLAG(CHANNEL) = HS_A_TMP_RED_FLAG(CHANNEL) + 1
                Else
		 if (qualflags(ix) .eq. 1) 
	1	 HS_A_TMP_YEL_FLAG(CHANNEL) = HS_A_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 38) then
		if (qualflags(ix) .eq. 2) then
	 	  HS_B_TMP_RED_FLAG(CHANNEL) = HS_B_TMP_RED_FLAG(CHANNEL) + 1
                Else
		  if (qualflags(ix) .eq. 1) 
	1	  HS_B_TMP_YEL_FLAG(CHANNEL) = HS_B_TMP_yel_FLAG(CHANNEL) + 1
                Endif
              Endif 
              if (ix .eq. 39) then
		if(qualflags(ix) .eq. 2) then 
		  MTM_A_RED_FLAG(CHANNEL) = MTM_A_RED_FLAG(CHANNEL) + 1
                else
		  if(qualflags(ix) .eq. 1) 
	1	   MTM_A_YEL_FLAG(CHANNEL) = MTM_A_yel_FLAG(CHANNEL) + 1
                Endif
              Endif
              if (ix .eq. 40) then
		if (qualflags(ix) .eq. 2)then 
		  MTM_B_RED_FLAG(CHANNEL) = MTM_B_RED_FLAG(CHANNEL) + 1
                else
		  if (qualflags(ix) .eq. 1) 
	1          MTM_B_YEL_FLAG(CHANNEL) = MTM_B_YEL_FLAG(CHANNEL) + 1
                endif
              endif 
              if (ix .eq. 41) then
		if (qualflags(ix) .eq. 2) then 
		  BOL_RED_FLAG(CHANNEL) = BOL_RED_FLAG(CHANNEL) + 1
                else
                  if (qualflags(ix) .eq. 1) 
	1	  BOL_YEL_FLAG(CHANNEL) = BOL_YEL_FLAG(CHANNEL) + 1
                endif
              endif
              if ((ix .ge. 42) .and. (ix .le. 51)) then
                LX=IX - 41
                if (qualflags(ix) .eq. 2) then 
		Ipdu_vol_red_FLAG(lx,CHANNEL) = Ipdu_vol_red_FLAG(lx,CHANNEL) + 1
                else
                 if (qualflags(ix) .eq. 1) 
	1	Ipdu_vol_yel_FLAG(lx,CHANNEL) = Ipdu_vol_yel_FLAG(lx,CHANNEL) + 1
                Endif
              endif
              if ((ix .ge. 52) .and. (ix .le. 57)) then
                lx=ix - 51
                if (qualflags(ix) .eq. 2) then
		Ipdu_cur_red_FLAG(lx,CHANNEL) = Ipdu_cur_red_FLAG(lx,CHANNEL) + 1
                else
		 if (qualflags(ix) .eq. 1) 
	1	Ipdu_cur_YEL_FLAG(lx,CHANNEL) = Ipdu_cur_yel_FLAG(lx,CHANNEL) + 1
                endif
              endif
C              if ((ix .ge. 58) .and. (ix .le. 61)) then
C                if (qualflags(ix) .ne. 0) 
C	1	stchg_red_FLAG(CHANNEL) = stchg_red_FLAG(CHANNEL) + 1
C              endif
	    endif
	Enddo
	If ( qualflags(108) .ne. zero ) Qualflag_Counter(108,Channel) =
	1	Qualflag_Counter(108,Channel) + 1

	If ( qualflags(109) .gt. fac_good_dq) Qualflag_Counter(109,Channel) =
	1	Qualflag_Counter(109,Channel) + 1
	If ( qualflags(110) .gt. fac_good_dq) Qualflag_Counter(110,Channel) =
	1	Qualflag_Counter(110,Channel) + 1
	
	Flag = FLG_MICRO - 1
	Do Ix = 1, 2
	  Flag = Flag + 1
	  Beg = 8 * ( Ix - 1 )
	  Fin = Beg + 7
          Qual(1) = Qualflags(flag)
	  Kx = -1
	  Do Jx = Beg, Fin
	    Kx = Kx + 1
            If ( Btest(I2_Qual,Kx) )      
	1	Micro_Flag(Jx,Channel) = Micro_Flag(Jx,Channel) + 1
	  Enddo
	Enddo
        Qual(1) = Qualflags(FLG_ATT_ALT_RED)
	Do Ix = 0, 7     
	  If ( Btest(I2_Qual,Ix) )
	1	Red_Att_Flag(Ix,Channel) = Red_Att_Flag(Ix,Channel) + 1
	Enddo
        Qual(1) = Qualflags(FLG_ATT_ALT_YEL)
	Do Ix = 0, 7     
	  If ( Btest(I2_Qual,Ix) )
	1	Yel_Att_Flag(Ix,Channel) = Yel_Att_Flag(Ix,Channel) + 1
	Enddo
	Flag = FLG_GRTRED_ST - 1
	Do Ix = 1, 8
	  Flag = Flag + 1
	  Beg = 8 * ( Ix - 1 )
	  Fin = Beg + 7
          Qual(1) = Qualflags(Flag)
	  Kx = -1
	  Do Jx = Beg, Fin
	    Kx = Kx + 1
	    If ( Btest(I2_Qual,Kx) )
	1     	Red_Grt_Flag(Jx,Channel) = Red_Grt_Flag(Jx,Channel) + 1
	  Enddo
	Enddo
	Flag = FLG_GRTYEL_ST - 1
	Do Ix = 1, 8
	  Flag = Flag + 1
	  Beg = 8 * ( Ix - 1 )
	  Fin = Beg + 7
          Qual(1) = Qualflags(Flag)
	  Kx = -1
	  Do Jx = Beg, Fin
	    Kx = Kx + 1
	    If (Btest(I2_Qual,Kx))  
	1     	Yel_Grt_Flag(Jx,Channel) = Yel_Grt_Flag(Jx,Channel) + 1
	  Enddo
	Enddo
	Flag = FLG_STCHG_MJ_ST - 1
	Do Ix = 1, 4
	  Flag = Flag + 1
	  Beg = 8 * ( Ix -1 )	
	  Fin = Beg + 7
          Qual(1) = Qualflags(Flag)
	  Kx = -1
	  Do Jx = Beg, Fin
            Kx = Kx + 1
	    If ( Btest(I2_Qual,Kx) ) 
	1	Mj_St_Flag(Jx,Channel) = Mj_St_Flag(Jx,Channel) + 1
	  Enddo
	Enddo
        
!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FUT_QUAL_SUMMARY = %loc(FUT_NORMAL)
	ELSE
	  FUT_QUAL_SUMMARY = %loc(FUT_ABERR)
	ENDIF

	RETURN
	END
