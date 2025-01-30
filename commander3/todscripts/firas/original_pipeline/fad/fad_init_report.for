	Integer*4  Function  FAD_Init_Report  (report_name, lun_rpt, cmdline,
	1                                      current_adt, version, arch_in,
	2                                      arch_cal, arch_out)

c------------------------------------------------------------------------------
c   Purpose: Open and initialize FAD report file.
c
c   Input Parameters:
c      character*45  report_name   -- name of the report to write
c      character*79  cmdline(2)    -- command line as a character string
c      integer*4     current_adt(2)  -- current system time in ADT format
c      character*6   version       -- version of FAD code.
c      character*14  arch_in       -- input archive name
c      character*15  arch_cal      -- reference file archive name
c      character*15  arch_out      -- output archive name
c
c   Output Parameters:
c      integer*4     lun_rpt  -- logical unit number for report file
c
c   Subroutines and Functions Called:
c
c      Lib$Signal
c      Lib$Establish
c      Cut_Register_Version
c      Cut_Display_Banner
c      Lib$Get_Lun
c      Time_LT
c      Lib$GetJpi
c      Sys$AscTim
c      Cut_Translate_Archive_Id
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      $SSdef        -- System status values
c      FUT_error     -- for error handling
c      $jpidef       -- for logical name translation
c
c   Author:  Larry Paris Rosen, Hughes STX, 19 April 1993
c------------------------------------------------------------------------------
	Implicit None

c Include files

	Include	 '($ssdef)'
	Include  '(fut_error)'
	Include  '($jpidef)'
	Include	 '(fad_msg)'

c Passed parameters

	Character*45	report_name
	Character*79	cmdline(2)			! Command line
	Integer*4	lun_rpt				! unit numbers
	Integer*4	current_adt(2)			! run time in adt
	Character*6	version
	Character*14	arch_in				! input data archive
	Character*15	arch_cal			! reference archive
	Character*15	arch_out			! output data archive

c Functions

	Integer*4  CUT_Register_Version, CUT_Display_Banner
	Integer*4  Lib$Get_Lun
	Integer*4  Time_LT
	Integer*4  Lib$GetJpi
	Integer*4  Sys$AscTim
	Integer*4  CUT_Translate_Archive_ID

c External

	External	fut_error
	External  	CUT_Translate_Archive_ID

c Local

	Integer*4	rstat			! return status
	Character*8	owner			! invoking user name
	Integer*2	time_len		! length of time string
	Character*32	current_time		! current system time string
	Character*72	logn, flogn, tlogn, blog ! logical name and translated
	Integer*4	flen, tlen		! length of translated logicals
	Integer*4	larc, larct		! length of archive names

c------------------------------------------------------------------------------

c Begin

	FAD_Init_Report = %Loc (FAD_Normal)

	rstat = Lib$Get_Lun (lun_rpt)
	If (rstat .NE. SS$_Normal) Then
	   FAD_Init_Report = %Loc (FAD_Abort)
	   CALL Lib$Signal (fad_lunerr, %Val(1), %Val(rstat))
	Else
	   FUT_Report_Lun = lun_rpt
	   Call Lib$Establish (fut_error)

	   Open (Unit=lun_rpt, File=report_name, Status='new',Form='formatted',
	1        Access='sequential', Organization='sequential', Iostat=rstat )

	   If (rstat .NE. 0) Then
	      FAD_Init_Report = %Loc (FAD_Abort)
	      Call Lib$Signal (FAD_Openrep, %Val(2), report_name, %Val(rstat))
	   Endif
	Endif
	If (FAD_Init_Report .EQ. %Loc (FAD_Normal)) Then
	   rstat = CUT_Register_Version (version)
	   rstat  = CUT_Display_Banner (lun_rpt, 80,
	1     'FIRAS  Facility  FAD_Apply_Destriper')

c Write user and current time to the report file.

	   rstat = Lib$GetJpi (jpi$_username,,,,owner,)
	   rstat = Sys$AscTim (time_len, current_time, current_adt, 0)
	   Write (lun_rpt, 10) owner, current_time (1:time_len)
  10	   Format (' Run by:   ', A, '   at  Time: ',A,/)

c Write command line with defaults to the report file.

	   Write (lun_rpt, 20) cmdline(1)
	   Write (lun_rpt, 20) cmdline(2)
  20	   Format (1X, A)
	   Write (lun_rpt, *)

c Write translation of logical names.

	   larc = Len (arch_in)
	   logn (1:larc) = arch_in
	   rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	   Write (lun_rpt, 30) 'Input', arch_in, tlogn(1:tlen)
  30	   Format (2X,'Logical Translation for ',A,' Archive:',/,5X,A,' = ',A)
	   larc = Len (arch_out)
	   logn = blog
	   logn (1:larc) = arch_cal
	   rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	   Write (lun_rpt, 30) 'Reference', arch_cal, tlogn(1:tlen)
	   larc = Len (arch_out)
	   logn = blog
	   logn (1:larc) = arch_out
	   rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	   Write (lun_rpt, 30) 'Output', arch_out, tlogn(1:tlen)
	Endif
	Return
	End
