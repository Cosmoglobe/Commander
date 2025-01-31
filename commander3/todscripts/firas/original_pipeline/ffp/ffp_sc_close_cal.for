	Integer*4  Function  FFP_SC_Close_Cal ( input, filexts, fnum, archin,
	1                                       inlun, outlun )
c------------------------------------------------------------------------------
c
c	Function FFP_SC_CLOSE_CAL
c
c	This function closes the FIRAS input fsl_cal data files and the ADB
c	output ffp data file.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		input	character*12			! Input filename base
c		filexts	Character*20 (fac_max_num)	! Input file extensions
c		fnum	integer*2		! Number of input files
c		archin  character*13		! input archive name
c		inlun	integer*4 (fac_max_num)	Input files luns
c		outlun	integer*4		Output file lun
c
c	Output:
c		none
c
c	Subroutines called:
c		lib$signal
c
c	Include files:
c		ffp_invoc_sky.txt
c		fut_params.txt
c		cct_query_catalog_record.txt
c
c------------------------------------------------------------------------------

	Implicit none

	Include '(fut_params)'
	Include '(ffp_invoc_sky)'
	Include 'CT$Library:CTUser.Inc'

c Passed parameters

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			!  Number of input maps
	Character*13	archin			!  Input archive name
	Integer*4	inlun (fac_max_num)	!  Input skymap luns
	Integer*4	outlun			!  Output skymap lun

c Local

	Integer*4	status			!  process status
	Integer*2	ct_stat(20)		!  return status
	Integer*4	fil			!  input file counter
	Character*47	infile
	Integer*4	io_status

C Externals

	External	ffp_normal
	External	ffp_closerr
C
C  Begin
C
	status = %loc(ffp_normal)
C
	Close ( outlun, iostat=io_status )
	If (io_status .NE. 0) Then
	   status = io_status
	   call Lib$Signal (ffp_closerr, %val(2), fcc_outfile(1:fcc_outlen),
	1              %loc(status))
	EndIf
	Do fil = 1, fnum
	   infile = archin // ':' // input // '.' // filexts (fil)
	   Call Ct_Close_Arcv ( , inlun (fil), ct_stat )
	   If (ct_stat(1) .NE. CTP_Normal) Then
	      status = ct_stat(1)
	      call Lib$Signal (ffp_closerr, %val(2), infile, %loc(status))
	   EndIf
	EndDo

	FFP_SC_Close_Cal = status

	Return
	End
