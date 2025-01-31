C-------------------------------------------------------------------------------
	INTEGER*4 FUNCTION fpp_init_report ( lun_rpt, report_file, command_line,
	1	version, arch_in, arch_out, arch_ref )
C-------------------------------------------------------------------------------
C	Purpose: To write initial information into the report file,
C		 primarily, the options selected from the command line.
C
C	Author: Shirley M. Read
C		STX, January, 1989
C
C	Invocation: Status = FPP_Init_Report ( Lun_rpt, Report_file,
C	1	Command_Line, Version, Arch_In, Arch_Out, Arch_Ref )
C
CH	Change Log:
CH
CH		Version 4.4.1 08/03/89 SER 3306, Q. C. Chung, STX
CH                     Provide version number to track software update.
CH
CH		Version 4.4.2 10/12/89, SPR 4711, R. Kummerer STX
CH			Correct record counts displayed in report file
CH			when processing by timerange.
CH
CH		Version 4.4.3 11/14/89, SPR 5034, R. Kummerer STX
CH			FPP must note when a channel is missing a raw science
CH			segment.
CH
CH		Version 4.4.4 11/29/89, SPR 5165, R. Kummerer STX
CH			Fails to skip missing segments.
CH
CH		Version 5.7 2/23/90, SPR 5876, H. Wang STX
CH			IFG TRACKING default not conform with rest of pipeline.
CH
CH		Version 7.4 7/24/90, SPR 4171, L. Rosen STX
CH			Standardize report filename.  The name is generated in
CH			FPP_PARSE_COMMAND.  The current GMT time is also
CH			gotten there, and passed here for use.
CH
CH		New Version  2/28/91, Larry Rosen, STX
CH			New requirements are to be met here.  Major rennovation.
CH			Command line has most of the required initial
CH			information in it.  Get user name. Get logicals.
C  ----------------------------------------------------------------------------
C	Input Parameters:
C	  Name			Type		Description
C	  ----------------------------------------------------------------------
C	  Report_File		C*33		Filename for report
C	  Command_Line(3)	C*79		Command line with defaults
C	  Version	 	C*6		Software version number from FPP
C	  Arch_In		C*(*)		Input data archive
C	  Arch_Out		C*(*)		Output data archive
C	  Arch_Ref		C*(*)		Reference data archive
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Lun_Rpt 	I*4		Logical unit for report file
C	
C	Subroutines Called:
C
C	  Lib$Get_Lun
C         Sys$Gettim
C	  CT_Binary_to_GMT
C	  Lib$Signal
C
C	Include Files:
C	  FPP_Msg.Txt
C         $SSDef
C	  FUT_Params.Txt
C
C	Processing Method: PDL for FPP_Init_Report
C
C	Get a logical unit number to open the report file.
C	Open the report file for write access.
C	Write banner, version number to report.
C	Establish error writer to report.
C	Get and write user name and system run time to report.
C	Write the command line with defaults to the report.
C	Write expansion of logical names.
C	Return with normal or error status.
C------------------------------------------------------------------------------	
	IMPLICIT NONE
	
C  Passed Parameters.

	INTEGER*4	lun_rpt		! logical unit for report file
	CHARACTER*33	report_file	! report filename
	CHARACTER*79	command_line(3)	! command line with defaults
	CHARACTER*6	version		! built version number from fpp
	CHARACTER*(*)	arch_in		! input data archive
	CHARACTER*(*)	arch_out	! output data archive
	CHARACTER*(*)	arch_ref	! reference data archive
	
C  Include files.

	INCLUDE		'(fut_params)'
	INCLUDE		'(fpp_msg)'
	INCLUDE		'($ssdef)'
	INCLUDE		'(fut_error)'
	INCLUDE		'($jpidef)'
	
C  Functions.

	INTEGER*4	lib$get_lun
	INTEGER*4	lib$establish
	INTEGER*4	lib$getjpi
	INTEGER*4	cut_register_version
	INTEGER*4	cut_display_banner
	INTEGER*4	sys$asctim
	INTEGER*4 	cut_translate_archive_id
	EXTERNAL  	cut_translate_archive_id
	
C  Local Declarations.

	EXTERNAL	fut_error
	INTEGER*4	status			! return status	
	INTEGER*4	rstatus			! return status
	INTEGER*4	zero / 0 /		! 0 value for fortran i/o status
	INTEGER*4	num_vol/80/
	CHARACTER*8	owner			! invoking user name
	INTEGER*4	current_time(2)		! current system adt time
	INTEGER*2	time_len		! length of time string
	CHARACTER*32	time			! current system time string
	CHARACTER*72	logn, flogn, tlogn	! logical name and translted
	INTEGER*4	flen, tlen		! length of translated logicals
	INTEGER*4	larc, larct		! length of archive names
	
C  Set the function status to Normal.

	fpp_init_report = %loc(fpp_normal)
	
C  Get a unit number and open the report file.

	status = lib$get_lun (lun_rpt)
	IF ( status .NE. ss$_normal ) THEN
	  fpp_init_report = %loc(fpp_aberr)
	  CALL lib$signal(fpp_getlunerr, %val(1), %val(status))
	ENDIF
	
C  Open the report file.

	IF ( fpp_init_report .EQ. %loc(fpp_normal) ) THEN
	  OPEN( unit=lun_rpt, file=report_file, status='new', form='formatted',
	1    access='sequential', organization='sequential', iostat=status )
	
	  IF ( status .NE. zero ) THEN
	    fpp_init_report = %loc(fpp_aberr)
	    CALL lib$signal(fpp_openerr,%val(1),%val(status))
	  ELSE
	    rstatus = cut_register_version(version)
            rstatus  = cut_display_banner(lun_rpt, num_vol,
	1     'Firas Facility FPP_Pre_Processor')
	    fut_report_lun = lun_rpt
	    status = lib$establish (fut_error)
	  ENDIF						! open ok
	ENDIF						! lun ok
	
C  Write user and current time to the report file.

	IF ( fpp_init_report .EQ. %loc(fpp_normal) ) THEN
	  rstatus = lib$getjpi (jpi$_username,,,,owner,)
	  CALL sys$gettim ( current_time )
	  rstatus = sys$asctim ( time_len, time, current_time, 0 )
	  WRITE (lun_rpt,20) owner, time
  20	  FORMAT (' Run by:   ', A, '   at  Time: ',A,/)
 
C  Write command line with defaults to the report file.

	  WRITE (lun_rpt,30) command_line(1), command_line(2), command_line(3)
  30	  FORMAT (3(1X,A,/))
 
	ENDIF		! function status is normal

C  Write translation of logical names.

	call str$upcase(arch_in,arch_in)
	larc = len (arch_in)
	logn(1:larc) = arch_in
	rstatus = cut_translate_archive_id(logn, flogn, flen, tlogn, tlen)
	WRITE (lun_rpt,40) arch_in, tlogn
  40	FORMAT(1X,'Logical Translation for Input Archive:',4X,A,' = ',/,25X,A)
	call str$upcase(arch_out,arch_out)
	larc = len (arch_out)
	logn(1:larc) = arch_out
	rstatus = cut_translate_archive_id(logn, flogn, flen, tlogn, tlen)
	WRITE (lun_rpt,50) arch_out, tlogn
  50	FORMAT(1X,'Logical Translation for Output Archive:',4X,A,' = ',/,25X,A)
	call str$upcase(arch_ref,arch_ref)
	larc = len (arch_ref)
	logn(1:larc) = arch_ref
	rstatus = cut_translate_archive_id(logn, flogn, flen, tlogn, tlen)
	WRITE (lun_rpt,60) arch_ref, tlogn
  60	FORMAT(1X,'Logical Translation for Reference Archive:',4X,A,' = ',/,25X,
	1    A,/)
	RETURN
	END
