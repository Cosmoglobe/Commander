	Integer*4  Function  FMS_INIT_REPORT  ( Lun_Rep, Reportfile,
	1                                       Rep_Default, Rlen, Version,
	2                                       Cmdline, Script, Slen,
	3                                       Trans_Freq, In_Skymap,
	4                                       File_Ext, Arch_In, Arch_Out,
	5                                       Arch_Ref, Tstart, Tstop,
	6                                       Current_Time )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_INIT_REPORT.FOR
C
C  Purpose:  Open and initialize the report file.  If the report name was
C            defaulted, then the name must be constructed.  The name will
C            include the "chanscan" specification of the output file which
C            doesn't get determined until FMS_VERIFY_MERGE.  So, in the mean
C            time, open a temporary report file 'FMS_REPORT.TEMPORARY' to
C            catch any error messages.  Later close, rename, and reopen the
C            report file for appending.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FMS_Msg)'
	Include		'(FUT_Error)'
	Include		'($JPIDef)'
	Include		'(CCT_Query_Catalog_Record)'

C  Passed parameters.

	Integer*4	Lun_Rep			! Report logical unit #
	Character*45	Reportfile		! Filename for report
	Logical*1	Rep_Default		! Flag for default report name
	Integer*2	Rlen			! Length of report file name
	Character*6	Version
	Character*79	Cmdline (3)		! Command line with defaults
	Character*80	Script			! Script file name
	Integer*2	Slen			! Length of script name string
	Integer*2	Trans_Freq	! Transition Low-High Frequency in icm
	Character*33	In_Skymap (2)		! Input skymaps
	Character*20	File_Ext		! Output file extension
	Character*14	Arch_In			! Input data archive
	Character*15	Arch_Out		! Output data archive
	Character*15	Arch_Ref		! C or D Vector archive
	Character*14	Tstart			! Start time in GMT
	Character*14	Tstop			! Stop time in GMT
	Integer*4	Current_Time(2)		! Current system adt time

C  Functions

	Integer*4	Lib$Get_Lun
	Integer*4	CCT_Query_Catalog
	Integer*4	Time_LT
	Integer*4	CUT_Register_Version, CUT_Display_Banner
	Integer*4	Lib$GetJPI
	Integer*4	Sys$AscTim
	Integer*4 	CUT_Translate_Archive_ID

C  External

	External	Fut_Error
	External  	Cut_Translate_Archive_Id

C  Local

	Integer*4	Rstat				! Return status
	Character*20	TempRep /'FMS_REPORT.TEMPORARY'/
	Dictionary 'ccm_cme_catalog_entry'
	Record /query_catalog/		Query_Cat
	Record /ccm_cme_catalog_entry/	Cats (3)
	Character*8	Owner			! Invoking user name
	Integer*2	i			! Counter
	Character*72	Logn, Flogn, Tlogn, Blog ! Logical name and translated
	Integer*4	Flen, Tlen		! Length of translated logicals
	Integer*4	Larc, Larct		! Length of archive names
	Integer*2	Time_Len		! Length of time string
	Character*32	Time			! Current system time string

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Init_Report = %Loc (FMS_Normal)

C  Get the report logical unit number.  Establish the error handler.

	Rstat = Lib$Get_Lun (Lun_Rep)
	If (Rstat .NE. SS$_Normal) Then
	   FMS_Init_Report = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	Else
	   FUT_Report_Lun = Lun_Rep
	   Call Lib$Establish (FUT_Error)

C  If Rep_Default is false, the report name was entered.  Open the report.

	   If (.NOT. Rep_Default) Then

	      Open ( Unit=Lun_Rep, File=Reportfile (1:Rlen), Status='New',
	1            Form='Formatted', Access='Sequential',
	2            Organization='Sequential', Iostat=Rstat )

	      If (Rstat .NE. 0) Then
	         FMS_Init_Report = %Loc (FMS_Abort)
	         Call Lib$Signal ( FMS_Openerr, %Val(2), Reportfile (1:Rlen),
	1                          %Val(Rstat) )
	      EndIf
	   Else

C  Else Rep_Default is true.  Report name was defaulted and must be
C  constructed.  In this case the name will include the "chanscan"
C  specification of the output file which doesn't get determined until
C  FMS_VERIFY_MERGE.  So, in the mean time, open a temporary report file
C  'FMS_REPORT.TEMPORARY' to catch any error messages.  Later close, rename,
C  and reopen the report file for appending.

	      Reportfile = TempRep
	      Rlen = 20
	      Open ( Unit=Lun_Rep, File=Reportfile, Status='New',
	1            Form='Formatted', Access='Sequential',
	2            Organization='Sequential', Iostat=Rstat )

	      If (Rstat .NE. 0) Then
	         FMS_Init_Report = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Openerr, %Val(2), TempRep, %Val (Rstat))
	      EndIf
	   EndIf
	EndIf

C  Register the version number and display banner in report.

	If (FMS_Init_Report .EQ. %Loc (FMS_Normal)) Then
	   Rstat = CUT_Register_Version (Version)
	   Rstat  = CUT_Display_Banner ( Lun_Rep, 80,
	1                                'FIRAS  Facility  FMS_Merge_Skymap' )

C  Write user and current time to the report file.

	   Rstat = Lib$GetJPI (JPI$_Username,,,,Owner,)
	   Call Sys$GetTim ( Current_Time )
	   Rstat = Sys$AscTim ( Time_Len, Time, Current_Time, 0 )
	   Write (Lun_Rep, 10) Owner, Time (1:Time_Len)
  10	   Format (' Run by:   ', A, '   at  Time: ',A,/)

C  Write command line with defaults to the report file.

	   Do i = 1, 3
	      Write (Lun_Rep, 20) Cmdline (i)
	   EndDo
  20	   Format (1X, A)
	   Write (Lun_Rep,*)

C  Write translation of logical names.

	   Call Str$Upcase (Arch_In, Arch_In)
	   Larc = Len (Arch_In)
	   Logn (1:Larc) = Arch_In
	   Rstat = CUT_Translate_Archive_ID (Logn, Flogn, Flen, Tlogn, Tlen)
	   Write (Lun_Rep, 30) 'Input', Arch_In, Tlogn (1:Tlen)
  30	   Format (2X, 'Logical Translation for ', A, ' Archive:', /, 5X,
	1          A, ' = ', A)
	   Call Str$Upcase (Arch_Out, Arch_Out)
	   Larc = Len (Arch_Out)
	   Logn = Blog
	   Logn (1:Larc) = Arch_Out
	   Rstat = CUT_Translate_Archive_ID (Logn, Flogn, Flen, Tlogn, Tlen)
	   Write (Lun_Rep, 30) 'Output', Arch_Out, Tlogn (1:Tlen)
	   Call Str$Upcase (Arch_Ref, Arch_Ref)
	   Larc = Len (Arch_Ref)
	   Logn = Blog
	   Logn (1:Larc) = Arch_Ref
	   Rstat = CUT_Translate_Archive_ID (Logn, Flogn, Flen, Tlogn, Tlen)
	   Write (Lun_Rep, 30) 'Reference', Arch_Ref, Tlogn (1:Tlen)

C  Write contents of script file to report.

	   Write (Lun_Rep, 40) Script (1:Slen)
 40	   Format (/,2X,'Contents of Script file, ',A)
	   Do i=1,2
	      Write (Lun_Rep, 50) i, In_Skymap (i)
	   EndDo
 50	   Format (5X, 'Skymap(', I1, ') = ', A)
	   Write (Lun_Rep, 60) File_Ext
 60	   Format (5x, 'Filext = ', A, /)

C  Query the catalog to get the times of earliest and latest spectra.
C  This will be used in the report name if it is defaulted.

	   Query_Cat.Archive_Id = Arch_In
	   Query_Cat.Filename = In_Skymap (1)
	   Rstat = CCT_Query_Catalog (Query_Cat, Cats(1))
	   If (.NOT. Rstat) Then
	      FMS_INIT_REPORT = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Query_Cat, %Val(1), %Val(Rstat))
	   Else
	      Query_Cat.Filename = In_Skymap (2)
	      Rstat = Cct_Query_Catalog (Query_Cat, Cats(2))
	      If (.NOT. Rstat) Then
	         FMS_INIT_REPORT = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Query_Cat, %Val(1), %Val(Rstat))
	      Else
	         If (Time_LT (Cats(1).Initial_Time, Cats(2).Initial_Time)) Then
	            Call CT_Binary_To_GMT (Cats(1).Initial_Time, Tstart)
	         Else
	            Call CT_Binary_To_GMT (Cats(2).Initial_Time, Tstart)
	         EndIf
	         If (Time_LT (Cats(1).Final_Time, Cats(2).Final_Time)) Then
	            Call CT_Binary_To_GMT ( Cats(2).Final_Time, Tstop )
	         Else
	            Call CT_Binary_To_GMT ( Cats(1).Final_Time, Tstop )
	         EndIf
	         Write (Lun_Rep, 70) Tstart, Tstop
 70	         Format ( 2X, 'Earliest spectrum time: ', A14, /,
	1                 2X, 'Latest spectrum time: ', A14, / )
	      EndIf
	   EndIf
	EndIf

C  Write transition frequency:

	If (FMS_Init_Report .EQ. %Loc (FMS_Normal)) Then
	   Write (Lun_Rep, 80) Trans_Freq
 80	   Format (2X, 'High/Low Transition Frequency: ',I2,' icm',/)
	EndIf
	Return
	End
