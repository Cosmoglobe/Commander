C------------------------------------------------------------------------

	Integer*4 Function FTB_Overlap_Init ( jstart, jstop,
	1				      report, report_lun )

C------------------------------------------------------------------------
C    PURPOSE: 
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX Incorporated
C            September 5, 1989
C
C    INVOCATION: status = FTB_Overlap_Init ( jstart, jstop,
C					     report, report_lun )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	JSTART(2)	I*4		Binary GMT start time.
C	JSTOP(2)	I*4		Binary GMT stop time.
C	REPORT		L*1		Flag to generate a report.
C	REPORT_LUN	I*4		Report file LUN.
C
C    SUBROUTINES CALLED: 
C	FUT_GET_LUN
C	UPM_GET_VALUE
C	SYS$ASCTIM
C	SYS$GETTIM
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: FUT_Error.Txt,
C		    UPM_Stat_Msg,
C		    $SSDef,
C		    $JPIDef
C
C------------------------------------------------------------------------
C
C Changes:
C
C	SPR 4969.  R. Kummerer, November 9, 1989.  Subtract 24 hours from
C		the user supplied start time to ensure OVERLAP of seeing
C		all overlaps and properly selecting run times for FPP.
C
C------------------------------------------------------------------------

	Implicit	None

	Include         '($SSDef)'
	Include         '($JPIDef)'
	Include         '(FUT_Params)'
	Include         '(UPM_Stat_Msg)'

	Integer         *4	status
	Integer         *4	tstatus

	Integer		*4	jstart(2)
	Integer		*4	jstop(2)
	Logical		*1	report
	Integer		*4	report_lun
	Character	*64	report_file

	Character	*14	start
	Character	*14	stop
	Integer		*4	ilen

	Integer         *2	time_len
	Character       *32	time
	Integer		*4	current_time(2)

	Integer		*4	one_day(2)

	Integer		*4	i
	Character	*8	owner
	Character	*9	day
	Character	*256	command_line
	Integer		*2	clen
	Logical		*4	found

        Integer		*4      num_vol/80/
        Character	*5	version
        Parameter	(version='6.8')

	Integer		*4	SYS$ASCTim
	Integer         *4      SYS$GetTim
	Integer         *4      SYS$BinTim
	Integer		*4	LIB$Get_LUN
	Integer         *4      LIB$AddX
	Integer		*4	LIB$GetJPI
        Integer		*4      CUT_Register_Version
        Integer		*4      CUT_Display_Banner
	Integer		*4	UPM_Get_Value
	Integer		*4	UPM_Present

	External	FUT_Normal
	External	FTB_Normal
	External	FTB_RMSOpenRep

C Parse JSTART.

	status = UPM_Get_Value ( 'JSTART', start, ilen )

	If (status .Eq. SS$_Normal) Then
	   If (ilen .Lt. 14) start = start(1:ilen) //
	1				fac_jstart_default(ilen+1:)
	   Call CT_GMT_To_Binary ( start, jstart )
	End If

C Subtract 24 hours from the start time as a precaution for properly detecting
C overlapping segments.  For example if today's date where given as the start
C time, CCT_QUERY_CATALOG would only see today's segment and not be able to
C check the overlap with yesterday's segment.  Presumably, data older than 24
C hours has already been processed.

	status = SYS$BinTim ( '1 0:', one_day )
	status = LIB$AddX ( jstart, one_day, jstart )

C Parse JSTOP.

	status = UPM_Get_Value ( 'JSTOP', stop, ilen )

	If (status .Eq. SS$_Normal) Then
	   If (ilen .Lt. 14) stop = stop(1:ilen) //
	1				fac_jstop_default(ilen+1:)
	   Call CT_GMT_To_Binary ( stop, jstop )
	End If

C Parse REPORT.

	report = .True.
	status = UPM_Get_Value ( 'REPORT', report_file, )

C Open the report file.

	If (report) Then

	   status = LIB$Get_LUN ( report_lun )

	   If (status .Eq. SS$_Normal) Then

	      Open ( Unit=report_lun, File=report_file, Status='new',
	1               IOstat=tstatus )

	      If (tstatus .Eq. 0) Then

C
C Write the processing report.
C
	         tstatus = CUT_Register_Version(version)
                 tstatus = CUT_Display_Banner(report_lun,num_vol,
	1				'FIRAS Facility FTB_Overlap')
	         Write (report_lun,600)

	         tstatus = LIB$GetJPI (JPI$_UserName,,,,owner,)
	         Write (report_lun,200) owner

	         Call SYS$GetTim ( current_time )

	         tstatus = SYS$ASCTim ( time_len, time, current_time, 0 )
	         If (tstatus .Ne. SS$_Normal) Then
	            Call LIB$Signal(%Val(tstatus))
	         End If

	         Write (report_lun,300) time

	         tstatus = UPM_Get_Value ( '$LINE', command_line, clen )
	         If (tstatus .Ne. SS$_Normal) Then
	            Call LIB$Signal(%Val(tstatus))
	         End If

	         If (clen .Gt. 119) Then
		    i = 120
	 	    found = .False.
		    Do While (i .ge. 1. And. .Not. found)
		       i = i - 1
		       If (command_line(i:i) .Eq. '/') found = .True.
		    End Do
	  	    Write (report_lun,500) command_line(1:i-1), 
	1				      command_line(i:clen)
	         Else
		    Write (report_lun,400) command_line(1:clen)
	         End If

	      Else
                 status = %Loc(FTB_RMSOpenRep)
	         Call LIB$Signal(FTB_RMSOpenRep,%Val(2),%Val(tstatus),
	1	 		 report_file)
	      End If

	   Else
	      Call LIB$Signal(%Val(status))
	   End If

	End If

C The format statements in full glory.

200	Format (' Run by:     ',a)
300	Format (' Run Time:   ', a)
400	Format (' Invocation: ',a,/)
500	Format (' Invocation: ',a,'-',/,'             ',a,/)
600	Format(/)

	FTB_Overlap_Init = status

	Return
	End
