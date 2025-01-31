	Integer*4  Function  FMS_PARSE  ( Script, Slen, Combo, Order,
	1                                 Model_Ext, Weight, Trans_Freq,
	2                                 Report, Reportfile, Rlen,
	3                                 Rep_Default, Specif, In_Skymap,
	4                                 File_Ext, Cmdline )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_PARSE.FOR
C
C  Purpose:  Parse the FMS command line, getting all the qualifier values and
C            constructing the full command line.  Open, read, and close the FMS
C            script file.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(UPM_Stat_Msg)'
	Include		'(FMS_Msg)'

C  Passed Parameters

	Character*80	Script			! Script file name
	Integer*2	Slen			! Length of script name string
	Character*9	Combo			! Combination type
	Character*6	Order			! Command combination order
	Character*12	Model_Ext		! Cal model soln extension
	Character*8	Weight			! C_VECTOR or D_VECTOR
	Integer*2	Trans_Freq	! Transition Low-High Frequency in icm
	Logical*1	Report			! Whether or not to report
	Character*45	Reportfile		! Filename for report
	Integer*2	Rlen			! Length of report file name
	Logical*1	Rep_Default		! Flag for default report name
	Character*4	Specif(2)		! "chanscan" of input skymaps
	Character*33	In_Skymap (2)		! Input skymaps
	Character*20	File_Ext		! Output file extension
	Character*79	Cmdline (3)		! Command line with defaults

C  Script file, Namelist definition.  Skymap is DIM(3) for error check.
C  Third must be blank.

	Namelist / Smaps / Skymap, Filext
	Character*33	Skymap (3)			! Input skymaps
	Character*20	Filext				! Output file extension

C  Functions

	Integer*4	UPM_Present
	Integer*4	UPM_Get_Value
	Integer*4	UPM_Get_Word
	Integer*4	Lib$Get_Lun

C  Local

	Character*79	Blankline
	Character*1	Blanks (79) / 79 * ' '/
	Equivalence	(Blankline, Blanks(1) )
	Integer*2	Cmdlines		! Number of lines in command
	Integer*4	Rstat			! Return status
	Integer*4	Txtlen			! Length of text string
	Integer*2	Cstart			! Next position in command line
        Integer*2	Cstop			! End of current Cmdline length
	Integer*2	Num_Len			! String length of number
	Character*10	Num_Str			! Number as a string
	Integer*2	Nstart		! Position of highest digit in string
	Integer*4	Lun_Script		! Logic unit # for script file
	Integer*2	i			! Counter
	Integer*2	Snum			! Number of skymaps input
	Integer*2	Period			! Location of . in filename

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Parse = %Loc (FMS_Normal)
	Cmdline(1) = Blankline
	Cmdline(2) = Blankline
	Cmdline(3) = Blankline
	Cmdline(1)(1:3) = 'FMS'
	Cmdlines = 1
	Cstop = 3
C
	Rstat = UPM_Present ('SCRIPT')
	If (Rstat .EQ. UPM_Pres) Then
	   Rstat = UPM_Get_Value ('SCRIPT', Script, Slen)
	   If (Rstat .EQ. SS$_Normal) Then
	      Call STR$Upcase (Script, Script)
	      Cstart = Cstop + 2
	      Cstop = Cstart + 7 + Slen
	      Cmdline(1)(Cstart:Cstop) = '/SCRIPT=' // Script (1:Slen)
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter script file name.' )
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	Else
	   Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                    'Must enter script file name.' )
	   FMS_Parse = %Loc (FMS_Abort)
	EndIf
C
	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('COMBINATION')
	   If (Rstat .EQ. UPM_Pres .OR. Rstat .EQ. UPM_Defaulted) Then
	      Rstat = UPM_Get_Value ('COMBINATION', Combo, Txtlen)
	      If (Rstat .EQ. SS$_Normal) Then
	         Call STR$Upcase (Combo, Combo)
	         Cstart = Cstop + 2
	         Cstop = Cstart + 12 + Txtlen
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 17 + Txtlen
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/COMBINATION=' //
	1           Combo (1:Txtlen)
	      Else
	         Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                          'Must enter combination type.')
	         FMS_Parse = %Loc (FMS_Abort)
	      EndIf
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter combination type.')
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	EndIf
C
	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('ORDER')
	   If (Rstat .EQ. UPM_Pres .OR. Rstat .EQ. UPM_Defaulted) Then
	      Rstat = UPM_Get_Value ('ORDER', Order, Txtlen)
	      If (Rstat .EQ. SS$_Normal) Then
	         Call STR$Upcase (Order, Order)
	         Cstart = Cstop + 2
	         Cstop = Cstart + 6 + Txtlen
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 11 + Txtlen
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/ORDER='// Order(1:Txtlen)
	      Else
	         Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                          'Must enter combination order.')
	         FMS_Parse = %Loc (FMS_Abort)
	      EndIf
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter combination order.')
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	EndIf
C
	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('MODEL_EXT')
	   If (Rstat .EQ. UPM_Pres .OR. Rstat .EQ. UPM_Defaulted) Then
	      Rstat = UPM_Get_Value ('MODEL_EXT', Model_Ext, Txtlen)
	      If (Rstat .EQ. SS$_Normal) Then
	         Call STR$Upcase (Model_Ext, Model_Ext)
	         Cstart = Cstop + 2
	         Cstop = Cstart + 12 + Txtlen
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 17 + Txtlen
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/MODEL_EXT=' //
	1           Model_Ext (1:Txtlen)
	      Else
	         Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                          'Must enter model extension.')
	         FMS_Parse = %Loc (FMS_Abort)
	      EndIf
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter model extension.')
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	EndIf
C
	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('WEIGHT')
	   If (Rstat .EQ. UPM_Pres .OR. Rstat .EQ. UPM_Defaulted) Then
	      Rstat = UPM_Get_Value ('WEIGHT', Weight, Txtlen)
	      If (Rstat .EQ. SS$_Normal) Then
	         Call STR$Upcase (Weight, Weight)
	         Cstart = Cstop + 2
	         Cstop = Cstart + 7 + Txtlen
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 12 + Txtlen
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop)='/WEIGHT='//Weight(1:Txtlen)
	      Else
	         Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                          'Must enter spectral weight.')
	         FMS_Parse = %Loc (FMS_Abort)
	      EndIf
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter spectral weight.')
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	EndIf
C
	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('TRANS_FREQ')
	   If (Rstat .EQ. UPM_Pres .OR. Rstat .EQ. UPM_Defaulted) Then
	      Rstat = UPM_Get_Word ('TRANS_FREQ', Trans_Freq)
	      If (Rstat .EQ. SS$_Normal) Then
	         Num_Len = Int (Log10 (Real (Trans_Freq))) + 1	! string length
	         Write (Num_Str,10) Trans_Freq		! convert # to string
  10	         Format (I10)
	         Nstart = 11 - Num_Len
	         Cstart = Cstop + 2
	         Cstop = Cstart + 11 + Num_Len
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 16 + Num_Len
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/TRANS_FREQ=' //
	1           Num_Str (Nstart:10)
	      Else
	         Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                          'Must enter transition frequency.')
	         FMS_Parse = %Loc (FMS_Abort)
	      EndIf
	   Else
	      Call Lib$Signal ( FMS_Parserr, %Val(1),
	1                       'Must enter transition frequency.')
	      FMS_Parse = %Loc (FMS_Abort)
	   EndIf
	EndIf

C  Get report file qualifier.  Note that if default report file name is
C  going to be used, that we do not yet know the file specification of
C  the output skymap and so cannot generate the complete name at this
C  time.  The trick is to use a dummy name and later rename it.

	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Rstat = UPM_Present ('REPORT')
	   If (Rstat .EQ. UPM_Pres) Then
	      Report = .TRUE.
	      Rstat = UPM_Get_Value ('REPORT', Reportfile, Rlen)
	      If (Rstat .EQ. UPM_Absent) Then
	         Rep_Default = .TRUE.
	         Cstart = Cstop + 2
	         Cstop = Cstart + 6
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 11
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/REPORT'
	      Else
	         Rep_Default = .FALSE.
	         Call STR$Upcase (Reportfile, Reportfile)
	         Cstart = Cstop + 2
	         Cstop = Cstart + 7 + Rlen
	         If (Cstop .GE. 80) Then
	            Cmdlines = Cmdlines + 1
	            Cstart = 5
	            Cstop = 12 + Txtlen
	         EndIf
	         Cmdline (Cmdlines)(Cstart:Cstop) = '/REPORT=' //
	1           Reportfile (1:Rlen)
	      EndIf
	   ElseIf (Rstat .EQ. UPM_Defaulted) Then
	      Report = .TRUE.
	      Rep_Default = .TRUE.
	      Cstart = Cstop + 2
	      Cstop = Cstart + 6
	      If (Cstop .GE. 80) Then
	         Cmdlines = Cmdlines + 1
	         Cstart = 5
	         Cstop = 11
	      EndIf
	      Cmdline (Cmdlines)(Cstart:Cstop) = '/REPORT'
	   Else
	      Report = .FALSE.
	      Rep_Default = .FALSE.
	   EndIf
	EndIf

C  Open, read, and close the FMS script file.

	If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	   Skymap (1) = ' '
	   Skymap (2) = ' '
	   Skymap (3) = ' '
	   Filext = ' '
	   Rstat = Lib$Get_Lun (Lun_Script)
	   If ( Rstat .NE. SS$_Normal ) Then
 	      FMS_Parse = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	   Else
	      Open ( Unit=Lun_Script, File=Script (1:Slen), Status='Old',
	1            Access='Sequential', Iostat=Rstat, READONLY )
	      If (Rstat .NE. 0) Then
	         FMS_Parse = %Loc (FMS_Abort)
	         Call Lib$Signal ( FMS_Openerr, %Val(2), Script (1:Slen),
	1                          %Val(Rstat) )
	      Else
	         Read (Lun_Script, Nml=Smaps, Iostat=Rstat)
	         If (Rstat .NE. 0) Then
	            FMS_Parse = %Loc (FMS_Abort)
	            Call Lib$Signal ( FMS_Readerr, %Val(2), Script (1:Slen),
	1                             %Val(Rstat) )
	         Else
	            Snum = 0
	            i = 1			! For each input skymap
	            Do While (i .LE. 3 .AND. FMS_Parse .EQ. %loc (FMS_Normal))
	               If (Skymap (i)(1:1) .NE. ' ') Then
	                  Snum = Snum + 1
	                  Call Str$Upcase ( In_Skymap (i), Skymap (i) )
	               EndIf
	               i = i + 1
	            EndDo

C  Make sure the number of skymaps is 2.

	            If (Snum .NE. 2) Then
	               Call Lib$Signal (FMS_Num_Maps)
	               FMS_Parse = %Loc (FMS_Abort)
	            EndIf

C  If the output file extension is given, use it, else use the first skymap's.

	            If (Filext(1:1) .NE. ' ') Then
	               File_Ext = Filext
	            Else
	               Period = Index ( In_Skymap (1), '.' )
	               File_Ext = In_Skymap (1) (Period+1 : Len (In_Skymap(1)))
	            EndIf
	            Close (Lun_Script, Iostat=Rstat)
	            If (Rstat .NE. 0) Then
	               FMS_Parse = %Loc (FMS_Abort)
	               Call Lib$Signal ( FMS_Closerr, %Val(2), Script (1:Slen),
	1                                %Val(Rstat) )
	            Else
	               Call Lib$Free_Lun (Lun_Script)
	            EndIf
	            If (FMS_Parse .EQ. %Loc (FMS_Normal)) Then
	               Specif (1) = In_Skymap (1) (9:12)
	               Specif (2) = In_Skymap (2) (9:12)
	            EndIf
	         Endif
	      EndIf
	   EndIf
	EndIf
        Return
        End
