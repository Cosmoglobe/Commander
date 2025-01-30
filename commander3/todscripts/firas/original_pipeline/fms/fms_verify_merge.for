	Integer*4  Function  FMS_VERIFY_MERGE ( Combo, Order, Weight, Report,
	1                                       Lun_Rep, Rep_Default,
	2                                       Reportfile, Specif, Specif_Out,
	3                                       Tstart, Tstop, Current_Time,
	4                                       Combo_Type, Order_Num )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_VERIFY_MERGE.FOR
C
C  Purpose:  Check whether the combination is valid.  Determine output skymap
C            specification.  If report file name was defaulted, determine the
C            name, close the current report file, rename it, and open it
C            for appending.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FMS_Msg)'

C  Passed parameters.

	Character*9	Combo			! Combination type cmd line
	Character*6	Order			! Command combination order
	Character*8	Weight			! C_VECTOR or D_VECTOR
	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rep			! Report logical unit #
	Logical*1	Rep_Default		! Flag for default report name
	Character*45	Reportfile		! Filename for report
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*4	Specif_Out	! "chanscan" specification of output
	Character*14	Tstart			! Start time in GMT
	Character*14	Tstop			! Stop time in GMT
	Integer*4	Current_Time(2)		! Current system adt time
	Character*1	Combo_Type		! Combination type
			! Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other
	Integer*2	Order_Num		! Command order as an integer

C  Funciton

	Integer*4	Lib$Rename_File

C  Local

C  Lists of allowed skymap combinations.  If it isn't in a list, it won't be
C  allowed.  Each list element contains four strings.  They are the two
C  input file specifications, the one output file specification, and the
C  combination type: FREQUENCY, SCAN_MODE, SIDE, HYBRID.

	Integer*2	Num_In_List (4) / 4, 17, 14, 3 / ! # of combos in each
	Character*4	List0 (4,4) / 'RHSF','RHLF','RHFA','SCAN',
	1                             'RLSF','RLFL','RLFA','SCAN',
	2                             'LHSF','LHLF','LHFA','SCAN',
	3                             'LLSF','LLFL','LLFA','SCAN' /
	Character*4	List1 (4,17) / 'RHSS','LHSS','HISL','SIDE',
	1                              'RHFA','LHFA','HIFA','SIDE',
	2                              'RLSS','LLSS','LOSL','SIDE',
	3                              'RLFA','LLFA','LOFA','SIDE',
	4                              'RLLF','LLLF','HRES','SIDE',
	5                              'RHSS','RHFA','RTHI','SCAN',
	6                              'RLSS','RLFA','RTLO','SCAN',
	7                              'LHSS','LHFA','LTHI','SCAN',
	8                              'LLSS','LLFA','LTLO','SCAN',
	9                              'RHSS','RLSS','RTSL','FREQ',
	1                              'RHFA','RLFA','RTFA','FREQ',
	2                              'LHSS','LLSS','LTSL','FREQ',
	3                              'LHFA','LLFA','LTFA','FREQ',
	4                              'RHSS','LLSS','14SL','HYBR',
	5                              'RHFA','LLFA','14FA','HYBR',
	6                              'RLSS','LHSS','23SL','HYBR',
	7                              'RLFA','LHFA','23FA','HYBR' /
	Character*4	List2 (4,14) / 'HISL','HIFA','HIGH','SCAN',
	1                              'LOSL','LOFA','LOWF','SCAN',
	2                              'RTHI','LTHI','HIGH','SIDE',
	3                              'RTLO','LTLO','LOWF','SIDE',
	4                              'RTHI','RTLO','RIGT','FREQ',
	5                              'LTHI','LTLO','LEFT','FREQ',
	6                              'RTSL','RTFA','RIGT','SCAN',
	7                              'LTSL','LTFA','LEFT','SCAN',
	8                              'HISL','LOSL','SLOW','FREQ',
	9                              'HIFA','LOFA','FAST','FREQ',
	1                              'RTSL','LTSL','SLOW','SIDE',
	2                              'RTFA','LTFA','FAST','SIDE',
	3                              'RTHI','LTLO','RHLL','HYBR',
	4                              'RTLO','LTHI','RLLH','HYBR' /
	Character*4	List3 (4,3) / 'HIGH','LOWF','LRES','FREQ',
	1                             'RIGT','LEFT','LRES','SIDE',
	2                             'SLOW','FAST','LRES','SCAN' /

	Logical*1	Found			! Specif in list
	Integer*2	i			! Counters
	Integer*2	In_Order		! List number = input order
	Character*4	Combo_List		! Combination type from list
	Integer*4	Rstat			! Return status
	Character*45	Tempfile		! Temporary file name
	Character*14	Current_GMT		! Current time in GMT
	Integer*2	Slen			! Length of string

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Verify_Merge = %Loc (FMS_Normal)

C  Determine the Order of each input skymap by which list in section V of the
C  design walkthrough document the filename is in.

	Found = .FALSE.

C  Checking List0

	i = 1					! Combination number in list
	Do While (.NOT. Found .AND. i .LE. Num_In_List (1))
	   If ( ( Specif (1) .EQ. List0 (1,i) .AND.
	1         Specif (2) .EQ. List0 (2,i) ) .OR.
	1       ( Specif (1) .EQ. List0 (2,i) .AND.
	1         Specif (2) .EQ. List0 (1,i) ) ) Then
	      Found = .TRUE.
	      In_Order = 0					! Input order
	      Combo_List = List0 (4,i)
	      Specif_Out = List0 (3,i)
	   EndIf
	   i = i + 1
	EndDo

C  Checking, List1

	If (.NOT. Found) Then
	   i = 1				! Combination number in list
	   Do While (.NOT. Found .AND. i .LE. Num_In_List (2))
	      If ( ( Specif (1) .EQ. List1 (1,i) .AND.
	1            Specif (2) .EQ. List1 (2,i) ) .OR.
	1          ( Specif (1) .EQ. List1 (2,i) .AND.
	1            Specif (2) .EQ. List1 (1,i) ) ) Then
	         Found = .TRUE.
	         In_Order = 1					! Input order
	         Combo_List = List1 (4,i)
	         Specif_Out = List1 (3,i)
	      EndIf
	      i = i + 1
	   EndDo
	EndIf

C  Checking, List2

	If (.NOT. Found) Then
	   i = 1				! Combination number in list
	   Do While (.NOT. Found .AND. i .LE. Num_In_List (3))
	      If ( ( Specif (1) .EQ. List2 (1,i) .AND.
	1            Specif (2) .EQ. List2 (2,i) ) .OR.
	1          ( Specif (1) .EQ. List2 (2,i) .AND.
	1            Specif (2) .EQ. List2 (1,i) ) ) Then
	         Found = .TRUE.
	         In_Order = 2					! Input order
	         Combo_List = List2 (4,i)
	         Specif_Out = List2 (3,i)
	      EndIf
	      i = i + 1
	   EndDo
	EndIf

C  Checking, List3

	If (.NOT. Found) Then
	   i = 1				! Combination number in list
	   Do While (.NOT. Found .AND. i .LE. Num_In_List (4))
	      If ( ( Specif (1) .EQ. List3 (1,i) .AND.
	1            Specif (2) .EQ. List3 (2,i) ) .OR.
	1          ( Specif (1) .EQ. List3 (2,i) .AND.
	1            Specif (2) .EQ. List3 (1,i) ) ) Then
	         Found = .TRUE.
	         In_Order = 3					! Input order
	         Combo_List = List3 (4,i)
	         Specif_Out = List3 (3,i)
	      EndIf
	      i = i + 1
	   EndDo
	EndIf

C  Check if map was found.

	If (.NOT. Found) Then
	   FMS_VERIFY_MERGE = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Illegal_Merge, %Val(2), Specif(1), Specif(2))
	EndIf

C  Verify that the Order of the skymaps is equal to the command qualifier.

	If (Order .EQ. 'ZEROTH') Then
	   Order_Num = 0
	ElseIf (Order .EQ. 'FIRST') Then
	   Order_Num = 1
	ElseIf (Order .EQ. 'SECOND') Then
	   Order_Num = 2
	ElseIf (Order .EQ. 'THIRD') Then
	   Order_Num = 3
	EndIf
	If (In_Order .NE. Order_Num) Then
	   FMS_VERIFY_MERGE = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Wrong_Order)
	EndIf

C  Verify that the combination matches the command line COMBINATION type.

	If (Combo_List .NE. Combo (1:4)) Then
	   FMS_VERIFY_MERGE = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Wrong_Combo)
	EndIf

C  If Order is ZEROTH then D-Vector weighting should be used.  Give a
C  warning if not.

	If ((Order .EQ. 'ZEROTH') .AND. (Weight .NE. 'DVECTOR')) Then
	   Call Lib$Signal (FMS_Zeroth_Weight)
	EndIf

C  Determine Combo_Type: Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other

	If (Order .EQ. 'ZEROTH') Then
	   Combo_Type = 'Z'

	ElseIf (((Specif (1)(2:2) .EQ. 'L') .AND. (Specif (2)(2:2) .EQ. 'L'))
	1  .OR. ((Specif (1)(1:2) .EQ. 'LO') .AND. (Specif (2)(1:2) .EQ. 'LO'))
	2  .OR. ((Specif(1)(3:4) .EQ. 'LO') .AND. (Specif(2)(3:4) .EQ. 'LO')))
	3  Then

	   Combo_Type = 'L'

	ElseIf (((Specif(1)(2:2) .EQ. 'H') .AND. (Specif(2)(2:2) .EQ. 'H')).OR.
	1   ((Specif (1)(1:2) .EQ. 'HI') .AND. (Specif (2)(1:2) .EQ. 'HI')).OR.
	2   ((Specif(1)(3:4) .EQ. 'HI') .AND. (Specif(2)(3:4) .EQ. 'HI'))) Then

	   Combo_Type = 'H'

	ElseIf (((Specif(1)(2:2) .EQ. 'H') .AND. (Specif(2)(2:2) .EQ. 'L')).OR.
	1       ((Specif(1)(2:2) .EQ. 'L') .AND. (Specif(2)(2:2) .EQ. 'H')).OR.
	2   ((Specif (1)(1:2) .EQ. 'HI') .AND. (Specif (2)(1:2) .EQ. 'LO')).OR.
	3   ((Specif (1)(1:2) .EQ. 'LO') .AND. (Specif (2)(1:2) .EQ. 'HI')).OR.
	4   ((Specif (1)(3:4) .EQ. 'HI') .AND. (Specif (2)(3:4) .EQ. 'LO')).OR.
	5   ((Specif(1)(3:4) .EQ. 'LO') .AND. (Specif(2)(3:4) .EQ. 'HI'))) Then

	   Combo_Type = 'F'

	Else
	   Combo_Type = 'M'
	EndIf

C  Write combination to the report file.

	If (Report .AND. (FMS_Verify_Merge .EQ. %Loc (FMS_Normal))) Then
	   Write (Lun_Rep, 10) Order, Combo, Specif (1), Specif (2), Specif_Out
  10	   Format ( 2X, 'The ', A, ' order ', A, ' combination is:', /,
	1           5X, A, ' + ', A, ' ==> ', A, / )

C  If report name is defaulted, construct the new name, close the old file,
C  rename it, and open it for appending.

	   If (Rep_Default) Then
	      Close (Lun_Rep, Iostat=Rstat)
	      If (Rstat .NE. 0) Then
	         FMS_Verify_Merge = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Closerr, %Val(2), Reportfile,%Val(Rstat))
	      Else
	         Tempfile = Reportfile

C  Get current run time in GMT string form.

	         Call CT_Binary_to_GMT ( Current_Time, Current_GMT )

C  Get length of Order string excluding blank.

	         If (Order_Num .EQ. 1 .OR. Order_Num .EQ. 3) Then
	            Slen = 5
	         Else
	            Slen = 6
	         EndIf

	         Reportfile = 'FMS_' // Order (1:Slen) // '_' // Specif_Out //
	1                     '_' // Tstart (1:7) // '_' // Tstop (1:7) //
	2                     '.REP_' // Current_GMT (1:9)

	         Rstat = Lib$Rename_File (Tempfile, Reportfile)
	         If (Rstat .NE. SS$_Normal) Then
	            FMS_Verify_Merge = %Loc (FMS_Abort)
	            Call Lib$Signal ( FMS_Renamerr, %Val(3), Tempfile,
	1                             Reportfile, %Val (Rstat) )
	         Else
	            Open ( Unit=Lun_Rep, File=Reportfile, Status='Old',
	1                  Form='Formatted', Access='Append',
	2                  Organization='Sequential', Iostat=Rstat )

	            If (Rstat .NE. 0) Then
	               FMS_Verify_Merge = %Loc (FMS_Abort)
	               Call Lib$Signal ( FMS_Openerr, %Val(2), Reportfile,
	1                                %val (Rstat) )
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	EndIf
	Return
	End
