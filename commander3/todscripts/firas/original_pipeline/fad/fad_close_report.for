	Integer*4  Function  FAD_Close_Report  ( lun_rpt, report_name )

c------------------------------------------------------------------------------
c   Purpose: Close the report file.
c
c   Input Parameters:
c      integer*4     lun_rpt       -- logical unit number for report file
c      character*45  report_name   -- name of the report to write
c
c   Output Parameters:
c      none
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c
c   Functions:
c      Lib$Signal
c
c   Author:  Larry Paris Rosen, Hughes STX, 20 April 1993
c------------------------------------------------------------------------------
	Implicit None

c Passed Parameters

	Integer*4	lun_rpt
	Character*45	report_name

c Include

	Include		'(FAD_msg)'

c Local

	Integer*4	rstat

c------------------------------------------------------------------------------

c Begin

	FAD_Close_Report = %Loc (FAD_Normal)

	Close (lun_rpt, Iostat=rstat)
	If (rstat .NE. 0) Then
	   FAD_Close_Report = %Loc (FAD_Abort)
	   Call Lib$Signal (fad_closerep, %Val(2), report_name, %Val(rstat))
	Endif
	Return
	End
