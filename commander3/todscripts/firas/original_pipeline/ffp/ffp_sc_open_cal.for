	Integer*4  Function  FFP_SC_Open_Cal ( input, filexts, fnum, archin,
	1                                      archout, inlun, outlun )
c------------------------------------------------------------------------------
c
c	Function FFP_SC_OPEN_CAL
c
c	This function opens the FIRAS input fsl_cal data files and the ADB
c	output ffp data file.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		input	character*12			! Input filename base
c		filexts	Character*20 (fac_max_num)	! Input file extensions
c		fnum	integer*2		! Number of input files
c		archin  character*13		! input archive name
c		archout character*14		! output archive name
c
c	Output:
c		inlun	integer*4 (fac_max_num)	Input files luns
c		outlun	integer*4		Output file lun
c
c	Subroutines called:
c		cct_query_catalog
c		ct_connect_read
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
	Include '(cct_query_catalog_record)'

c Passed parameters

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			!  Number of input maps
	Character*13	archin			!  Input archive name
	Character*14	archout			!  Output archive name
	Integer*4	inlun (fac_max_num)	!  Input skymap luns
	Integer*4	outlun			!  Output skymap lun

c Local

	Integer*4	status			!  return status
	Integer*4	rstatus			!  return status
	Integer*4	fil			!  input file counter
	Integer*4	cstatus			!  CT return status
	dictionary 'ccm_cme_catalog_entry'
	record /ccm_cme_catalog_entry/ cats(50)
	record /query_catalog/ query_cat
	Integer*4	io_stat			!  I/O return status
	Integer*2	len		!  Length of output record in longwords
	Character*47	infile
	Character*63	outfile

C Functions

	Integer*4	fut_get_lun
	Integer*4	cct_query_catalog
	Integer*4	ct_connect_read
	External	ct_connect_read

C Externals

	External	ffp_normal
	External	ffp_lunerr
	External	ffp_openerr
	External	ffp_ctnocatrecs
	External	ffp_ctquerycat
	External	fut_normal
	External	cct_q_no_cat_entry
	External	cct_normal
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin
C  Initialize the open.
C
	status = %loc(ffp_normal)
C
	If (fcc_coadd) then
	   len = fcc_ilen			! fil records
	Else
	   len = fcc_slen			! fsl records
	EndIf
c
c  Get the logical unit numbers for the files.
c
	rstatus = fut_get_lun (outlun)
	If (rstatus .NE. %loc(fut_normal)) Then
	   FFP_SC_Open_Cal = rstatus
	   Call lib$signal (ffp_lunerr, %val(1), %loc(rstatus))
	   Return
	EndIf
	Do fil=1,fnum
	   rstatus = fut_get_lun (inlun(fil))
	   If (rstatus .NE. %loc(fut_normal)) Then
	      FFP_SC_Open_Cal = rstatus
	      Call lib$signal (ffp_lunerr, %val(1), %loc(rstatus))
	      Return
	   EndIf
	EndDo
C
C  Open the output ADB file.
C
	outfile = archout // ':' // fcc_outfile (1:fcc_outlen)
	Open ( unit=outlun, file=outfile, iostat=io_stat, status='new',
	1      form='unformatted', recordtype='fixed', recl=len,
	2      access='sequential' )

	If (io_stat .NE. 0) Then
	   status = %loc(ffp_openerr)
	   Call lib$signal (ffp_openerr, %val(2), outfile, %val(io_stat))
	EndIf
	query_cat.archive_id = archin
	fil = 1
	Do While (fil .LE. fnum .AND. status .EQ. %loc(FFP_Normal))
C
C  Query the Cobetrieve catalog for input files.
C
	   query_cat.filename = input // '.' // filexts (fil)
	   cstatus = cct_query_catalog (query_cat, cats(1))
	   If (cstatus .EQ. %loc(cct_q_no_cat_entry)) Then
	      status = %loc(ffp_ctnocatrecs)
	      Call lib$signal (ffp_ctnocatrecs, %val(2), query_cat.archive_id,
	1                      query_cat.filename)
	   ElseIf (cstatus .NE. %loc(cct_normal)) Then
	      status = %loc(ffp_ctquerycat)
	      Call lib$signal (ffp_ctquerycat, %val(3), query_cat.archive_id,
	1                      query_cat.filename, %val(cstatus))
	   Else
C
C  Open the files.
C
	      infile = archin // ':' // input // '.' // filexts (fil)
	      Open (unit=inlun(fil), file=infile, iostat=io_stat, status='old',
	1           useropen=ct_connect_read)

	      If (io_stat .NE. 0) Then
	         status = %loc(ffp_openerr)
	         Call lib$signal (ffp_openerr, %val(2), infile, %val(io_stat))
	      EndIf
	   EndIf
	   fil = fil + 1
	EndDo

	FFP_SC_Open_Cal = status

	Return
	End
