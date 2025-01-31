C-------------------------------------------------------------------------------
	Integer*4 Function FPP_Check_Time ( Prev_Time, Sci_Rec, Anc_Rec, TF_R1,
	1   TF_R2, Start_Collect, End_Collect, QMax_Sec, Lun_Anc,
	2   Last_MNF_Time1, Last_MNF_Time2, Gap_Check, Anc_Gap, Report,
	3   Lun_Rpt, First_Anc )
C-------------------------------------------------------------------------------
C	Purpose: To check the validity of the computed midpoint of collect time,
C	         set the Badtime Flag in the science record if indicated.
C	         Perform consistency check: check that all minor frames in
C		 science collection have same scan mode.
C
C	Author: Shirley M. Read
C		STX, January, 1989
C
CH	Change Log:
CH
CH       SPR 3898  Modified : QCC/STX, 05 25, 1989
CH      Purpose :  to take care the time gap > two major frames while the
CH                 micro-processor is not function properly for any reason.
CH
CH	New Requirements/Design, Larry P. Rosen, STX, April 1991
CH	Many changes.  See PDL.
C------------------------------------------------------------------------------
C	Input Parameters:
C	  Name		  Type		  Description
C 	  ----------------------------------------------------------------------
C	  Prev_Time(2)	I*4		Time of previous science record
C	  Sci_Rec	NFS_SDF struc	Current FIRAS science record
C	  Anc_Rec(2)	NFS_ANC struc	Current set of FIRAS Word 31 records
C	  TF_R1(2)	L*1		Flag science telemetry format each frame
C	  TF_R2(2)	L*1		Flag science telemetry format each frame
C	  Start_Collect	I*4(2)		Start time of science collection
C	  End_Collect	I*4(2)		End time of science collection
C	  QMax_Sec	I*4(2)		Max seconds allowed else gap.
C	  Lun_Anc	I*4		Logical unit of ancillary file
C	  Last_MNF_Time1(2)	I*4	Quadword time end of MJF 1
C	  Last_MNF_Time2(2)	I*4	Quadword time end of MJF 2
C	  Report          L*1             Flag to write report
C	  Lun_Rpt         I*4             Logical unit for report
C	  Anc_gap	  I*4             Max number seconds for large gap reprt
C	  Gap_Check	  I*2		  Flag to check for sci and anc gaps.
C	  First_Anc       L*1             Flag for first anc record
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Sci_Rec	NFS_SDF struc	Sci. Rec. with badtime flag set or not
C
C	Subroutines Called:
C
C	  Lib$Signal
C
C	Include Files:
C	  FPP_Msg.Txt
C	  $SSDef
C
C	Processing Method: PDL for FPP_Check_Time
C
C	IF (midpoint collect time .GE. transmit time) THEN
C	    Set the "Badtime Flag" in the science record to Six.
C	ENDIF
C
C	IF (midpoint collect time .LT. transmit time of previous record)
C	THEN
C	    Set the "Badtime Flag" in the science record to Seven.
C	ENDIF
C
C  *	Determine number of minor frames in the consistency check window by
C	subtracting the start collect time from end collect time and dividing
C	by the number of seconds per frames.
C
C  *	If end of collection goes beyond end of Anc Record 2 Frame 2, then
C	read another Anc record and swap buffers.
C
C  *	Determine Record and Frame numbers for Start of collect.  Compare start
C	collect time with frame times + 32 sec.  Bad telemetry format will be
C	caught by ELSEIF structure, since non-science format gives earlier
C	(wrong) times.
C
C  *	Now have frame and record for start of collect.  For all minor frames in
C	collection, extract the scan speed and length and compare with reference
C	time determined.  If inconsistent set BadTime flag = 8.
C
C  *	Do for each minor frame.  Extract the scan speed and length and compare
C	with those determined by the reference time.  If inconsistent set
C	BadTime flag = 8.  If collection crosses frame boundary (mnf=129 or 257)
C	, update pointers to correct frame and record.  Also check if there is
C	a gap, and if there is then set Flag to Nine.
C
C	RETURN
C------------------------------------------------------------------------------
	Implicit None

C  Passed Parameters.
	Integer*4  	Prev_Time(2)    ! Time of previous science record
	Dictionary      'NFS_SDF'
	Record         /NFS_SDF/ Sci_Rec    ! Current FIRAS science record
	Dictionary      'NFS_ANC'       
	Record	       /NFS_ANC/ Anc_Rec(2) ! Currect set of FIRAS Word 31 recs.
	Logical*1	TF_R1(2)	! Flag science telemetry format record 1
	Logical*1	TF_R2(2)	! Flag science telemetry format record 2
	Integer*4	Start_Collect(2)   ! Start time of science collection
	Integer*4	End_Collect(2)     ! End time of science collection
	Integer*4	QMax_Sec(2)	! Max seconds allowed else gap.
	Integer*4	Lun_Anc		! Logical unit of ancillary file
	Integer*4	Last_MNF_Time1(2)	! Quadword time end of MJF 1
	Integer*4	Last_MNF_Time2(2)	! Quadword time end of MJF 2
	Logical*1	Report		! Flag to write report
	Integer*4	Lun_Rpt		! Unit number for report
	Integer*4	Anc_Gap		! max seconds for report of gaps in anc
	Integer*2	Gap_Check	! Flag to look for gaps in sci and anc
	Logical*1	First_Anc	! Flag first read of Anc for the day

C  Include files.
	Include      	'(FPP_Msg)'
	Include      	'($SSDef)'
	Include		'CT$Library:CTUser.Inc'

C  Functions
	Logical*1	Time_LT		! CT time less than function
	Logical*1	Time_GE		! CT time greater than or equal function
	Logical*1	Time_GT		! CT time greater than function
	Integer*4       Lib$Addx
	Integer*4       Lib$Subx
	Integer*4	Lib$EMul
	Integer*4	Lib$EDiv
	Integer*4	FPP_Read2

C  Local Declarations.

	Logical*1	First_Time /.True./
	Integer*4	Last_MNF_Sec    ! Time increment for last minor frame
	Parameter	( Last_MNF_Sec = 317500000 )
	Integer*4	Addend / 0 /
	Integer*4	Del_Last(2)	! Time to end of frame in quadword
	Integer*4	One /1/, Six / 6 /, Seven / 7 /, Eight / 8 /, Nine / 9 /
	Integer*4	retstat		! function return statuses
	Integer*4	len /2/		! parameter for Lib$Subx
	Integer*4	Time_Diff(2)	! Difference of two times
	Integer*4	Sec_Per_Mnf /2500000/
	Integer*4	NMnFrames	! Number of minor frames in collection
	Integer*4	Sec_32(2)
	Data		Sec_32(1) /320000000/
	Integer*4	End_Time(2)	! End of frame 2 time (either record)
	Integer*2	CT_Stat(20)	! Cobetrieve read status
	Character*14	Stime		! A binary time as a gmt time
	Integer*2	Rec		! Record for performing check
	Integer*2	Frm		! Frame for performing check
	Integer*2	Index_Start	! Start index in Status Bits array for
					! consistency checking
	Integer*4	Remainder	! of Lib$Ediv.
	Integer*2	Index_Stop	! Stop index in Status Bits array for
					! consistency checking
	Integer*2	mnf		! Counter for minor frames in collection
	Integer*2	Index		! MNF index to Status Bits array
	Integer*4	Statbits	! I4 storage for status bits
	Integer*2	Length	 	! MTM scan length
	Integer*2	Speed		! MTM scan speed

C Set the function status to Normal.

	FPP_Check_Time = %loc(FPP_Normal)

	If (First_Time) Then
	   Retstat = Lib$Emul (One, Last_MNF_Sec, Addend, Del_Last)
	   If (.Not. Retstat) Then
	      Call Lib$Signal(FPP_EmulErr, %val(1), %val(Retstat))
	      FPP_Check_Time = %loc(FPP_Aberr)
	   EndIf
	   First_Time = .False.
	EndIf

C  Check if the midpoint of collect time for the science record is greater
C  or equal to the transmit frame time.

	If ( Time_GE ( Sci_Rec.Collect_Time.Midpoint_Time,
	1              Sci_Rec.CT_Head.Time )) Then
	   Sci_Rec.Collect_Time.Badtime_Flag = Six
	   Return
	EndIf

C  Check if the midpoint of collect time is less than the transmit frame 
C  time for the previous record. 

	If ( TIME_LT (Sci_Rec.Collect_Time.Midpoint_Time, Prev_Time)) Then
	   Sci_Rec.Collect_Time.Badtime_Flag = Seven
	   Return
	EndIf

C  Determine number of minor frames in the consistency check window (the entire
C  collection) by subtracting the start collect time from end collect time and
C  dividing by the number for seconds per frame.

	Retstat = Lib$Subx (End_Collect, Start_Collect, Time_Diff, len)
	If ( Retstat .Ne. ss$_normal) Then
	   Call Lib$Signal (fpp_subxerr,%val(1),%val(retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf
	Retstat = Lib$Ediv (Sec_Per_Mnf, Time_Diff, NMnFrames, Remainder)
	If ( Retstat .Ne. ss$_normal) Then
	   Call Lib$Signal (fpp_Ediverr,%val(1),%val(retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf

C  Need to check if collection crosses frame boundaries.  First check if end of
C  collect goes beyond end of Anc rec2, frame 2.  If yes, read one more record.

	Retstat = Lib$Addx (Anc_Rec(2).GMT_MJF2, Sec_32, End_Time, Len)
	If ( Retstat .Ne. ss$_normal) Then
	   Call Lib$Signal (fpp_Addxerr,%val(1),%val(retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf
	If (Time_GT (End_Collect, End_Time).AND.TF_R2(2)) Then
	   Retstat = FPP_Read2 ( Lun_Anc, Anc_Rec, TF_R2, TF_R1, Gap_Check,
	1     Anc_Gap, Report, Lun_Rpt, Last_MNF_Time1, Last_MNF_Time2,
	2     First_Anc )
	   If (Retstat .NE. %loc(FPP_Normal)) Then
	      FPP_Check_Time = %loc(FPP_Aberr)
	      Return
	   EndIf
	EndIf

C  Determine the Anc Record number and Frame number for the Start Collect time.

C    First check if start_collect before Anc record 1, frame 1.

	If (Time_LT (Start_Collect, Anc_Rec(1).CT_Head.Time)) Then
	   Call CT_Binary_To_GMT (Start_Collect, stime)
	   Call Lib$Signal (FPP_ColBefAnc, %val(1), stime)
	   Sci_Rec.Collect_Time.Badtime_Flag = eight
	   Return
	EndIf

C    Get end of frame 1 record 1, and see if start is less than it.

	Retstat = Lib$Addx (Anc_Rec(1).CT_Head.Time, Sec_32, End_Time, Len)
	If ( Retstat .Ne. ss$_normal) Then
	   Call Lib$Signal (fpp_Addxerr,%val(1),%val(retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf
	If (Time_LT (Start_Collect, End_Time)) Then
	   Rec = 1
	   Frm = 1
	Else

C    Get end of frame 2 record 1, and see if start is less than it.

	   Retstat = Lib$Addx (Anc_Rec(1).GMT_MJF2, Sec_32, End_Time, Len)
	   If ( Retstat .Ne. ss$_normal) Then
	      Call Lib$Signal (Fpp_Addxerr,%val(1),%val(retstat) )
	      Fpp_Check_Time = %loc(FPP_Aberr)
	      Return
	   EndIf
	   If (Time_LT (Start_Collect, End_Time)) Then
	      Rec = 1
	      If (Time_GE (Start_Collect, Anc_Rec(1).GMT_MJF2)) Then
	         Frm = 2
	      ElseIf (TF_R1(1)) Then	! end of frame 1 (check for gap later)
	         Frm = 1
	      Else			! Fails due to bad Telemetry format
	         Sci_Rec.Collect_Time.Badtime_Flag = eight
	         Return
	      EndIf
	   Else

C    Get end of frame 1 record 2, and see if start is less than it.

	      Retstat = Lib$Addx (Anc_Rec(2).CT_Head.Time, Sec_32, End_Time,Len)
	      If ( Retstat .NE. ss$_normal) Then
	         Call Lib$Signal (fpp_Addxerr,%val(1),%val(retstat) )
	         Fpp_Check_Time = %loc(FPP_Aberr)
	         Return
	      EndIf
	      If (Time_LT (Start_Collect, End_Time)) Then
	         If (Time_GE (Start_Collect,Anc_Rec(2).CT_Head.Time) ) Then
	            Rec = 2
	            Frm = 1
	         ElseIf (TF_R1(2)) Then		! end of frame 2 record 1
	            Rec = 1
	            Frm = 2
	         Else
	            Sci_Rec.Collect_Time.Badtime_Flag = eight
	            Return
	         EndIf
	      Else

C    Get end of frame 2 record 2, and see if start is less than it.

	         Retstat = Lib$Addx(Anc_Rec(2).GMT_MJF2, Sec_32, End_Time, Len)
	         If ( Retstat .Ne. ss$_normal) Then
	            Call Lib$Signal (fpp_Addxerr,%val(1),%val(retstat) )
	            Fpp_Check_Time = %loc(FPP_Aberr)
	            Return
	         EndIf
	         If (Time_LT (Start_Collect, End_Time)) Then
	            Rec = 2
	            If (Time_GE (Start_Collect, Anc_Rec(2).GMT_MJF2)) Then
	               Frm = 2
	            ElseIf (TF_R1(1)) Then	! end of frame 1
	               Frm = 1
	            Else		! Fails due to bad Telemetry format
	               Sci_Rec.Collect_Time.Badtime_Flag = eight
	               Return
	            EndIf
	         Else
	            Call CT_Binary_To_GMT (Start_Collect, stime)
	            Call Lib$Signal (FPP_ColAftAnc, %val(1), stime)
	            Sci_Rec.Collect_Time.Badtime_Flag = eight
	            Return
	         EndIf
	      EndIf
	   EndIf
	EndIf

C  Determine index of start_collect.

	If (Frm .EQ. 1) Then
	   Retstat = Lib$Subx (Start_Collect, Anc_Rec(Rec).CT_Head.Time,
	1     Time_Diff, Len)
	Else
	   Retstat =Lib$Subx (Start_Collect,Anc_Rec(Rec).GMT_MJF2,Time_Diff,Len)
	EndIf
	If ( Retstat .NE. SS$_Normal) Then
	   Call Lib$Signal (FPP_Subxerr, %val(1), %val(Retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf
	Retstat = Lib$Ediv (Sec_per_MNF, Time_Diff, Index_Start, Remainder)
	If ( Retstat .NE. SS$_Normal) Then
	   Call Lib$Signal (FPP_EDiverr, %val(1), %val(Retstat) )
	   Fpp_Check_Time = %loc(FPP_Aberr)
	   Return
	EndIf
	If (Index_Start .LT. 128) Index_Start = Index_Start + 1

C  Determine index of End_collect by adding number minor frames calcuted above

	Index_Stop = Index_Start + NMnFrames - 1

C  Do for each minor frame.  Extract the scan speed and length and compare with
C  reference time determined.  If inconsistent set BadTime flag = 8.
C  If collection crosses frame boundary (mnf=129 or 257), update pointers to
C  correct frame and record.  Also check if there is a gap, and if there is
C  then set Flag to Nine.

	Do mnf = Index_Start, Index_Stop
	   If (mnf .GT. 128) Then
	      If (mnf.EQ.129 .OR. mnf.EQ.257) Then
	         If (Frm .EQ. 1) Then
	            Retstat = Lib$Subx ( Anc_Rec(rec).GMT_MJF2,
	1              Anc_Rec(rec).CT_Head.Time, Time_Diff, Len)
	            If ( Retstat .NE. SS$_Normal) Then
	               Call Lib$Signal (FPP_Subxerr, %val(1), %val(Retstat) )
	               Fpp_Check_Time = %loc(FPP_Aberr)
	               Return
	            EndIf
	            Frm = 2
	         Else
	            Retstat = Lib$Subx ( Anc_Rec(2).CT_Head.Time,
	1              Anc_Rec(1).GMT_MJF2, Time_Diff, Len)
	            If ( Retstat .NE. SS$_Normal) Then
	               Call Lib$Signal (FPP_Subxerr, %val(1), %val(Retstat) )
	               Fpp_Check_Time = %loc(FPP_Aberr)
	               Return
	            EndIf
	            Rec = 2
	            Frm = 1
	         EndIf
	         If (Time_GT (Time_Diff, QMax_Sec)) Then	! Gap
	            Sci_Rec.Collect_Time.Badtime_Flag = Nine
	            Return
	         EndIf
	      EndIf
	      If (mnf.GT.256) Then
	         Index = mnf - 256
	      Else
	         Index = mnf - 128
	      EndIf
	   Else
	      Index = mnf
	   EndIf
	   StatBits = Anc_Rec(Rec).Frame(Frm).Minor_Frame_Status_Bits(Index)
	   Speed = Jibits(StatBits,4,1)
	   Length = Jibits(StatBits,5,1)
	   If (( Sci_Rec.Sci_Head.MTM_Speed .NE. Speed ) .OR.
	1     ( Sci_Rec.Sci_Head.MTM_Length .NE. Length )) Then
	      Sci_Rec.Collect_Time.Badtime_Flag = Eight
	      Return
	   EndIf
	EndDo

	Return
	End
