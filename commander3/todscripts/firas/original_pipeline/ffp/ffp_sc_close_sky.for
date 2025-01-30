	Integer*4  Function  FFP_SC_Close_Sky ( inlun, outlun, fnum,
	1                                       input, filexts )
c------------------------------------------------------------------------------
c	Function FFP_SC_CLOSE_SKY
c
c	This function closes the FIRAS input spectrum skymaps and the ADB
c	output spectrum skymap.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		inlun	integer*4 (fac_max_num)	Input skymap lun's
c		outlun	integer*4		Output skymap lun
c		fnum	integer*2		Number of input skymaps
c		input	character*12			! Input filename base
c		filexts	Character*20 (fac_max_num)	! Input file extensions
c
c	Output:
c		none
c
c	Subroutines called:
c		csa_close_skymap
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		ffp_invoc_sky.txt
c		fut_params.txt
c
c------------------------------------------------------------------------------

	Implicit None

	Include '(fut_params)'
	Include '(ffp_invoc_sky)'

c Passed parameters

	Integer*4	inlun (fac_max_num)	!  Input skymap luns
	Integer*4	outlun			!  Output skymap lun
	Integer*2	fnum			!  Number of input maps
	Character*12	input			!  Input filename base
	Character*20	filexts (fac_max_num)	!  Input file extensions

c Local

	Integer*2	fil			!  File counter
	Integer*4	cstatus			!  CT or CSA return status
	Integer*4	io_stat			!  I/O return status
	Integer*4	rstatus			!  return status
	Integer*4	status			!  return status
	Character*33	infile			!  Input file name
	Character*1	blanks(33) /33*' '/
	Character*33	blankname		!  Blank input file name
	Equivalence	(blanks, blankname)

c Functions

	Integer*4	csa_close_skymap
	Integer*4	fut_free_lun

c External

	External	csa_normal
	External	ffp_csaclose
	External	ffp_normal
	External	fut_normal

	status = %loc(ffp_normal)
c
c  Close the input FIRAS skymap files and free the luns.
c
	fil = 1
	Do While (fil .LE. fnum .AND. status .EQ. %loc (ffp_normal))
	   infile = blankname
	   infile = input // '.' // filexts (fil)
	   cstatus = csa_close_skymap (inlun (fil), fac_skymap_no_levels)
	   If (cstatus .NE. %loc(csa_normal)) Then
	      status = %loc(ffp_csaclose)
	      Call lib$signal (ffp_csaclose, %val(2), infile, %val(cstatus))
	   Else
	      rstatus = fut_free_lun (inlun(fil))
	      If (rstatus .NE. %loc(fut_normal)) then
	         Call lib$signal ('Error freeing lun, status = ',%val(rstatus))
	      EndIf
	   EndIf
	   fil = fil + 1
	EndDo
c
c  Close the output ADB skymap file.
c
	cstatus = csa_close_skymap (outlun, fac_skymap_no_levels)
	If (cstatus .ne. %loc(csa_normal)) Then
	   status = %loc(ffp_csaclose)
	   Call lib$signal ( ffp_csaclose, %val(2), fcc_outfile(1:fcc_outlen),
	1                    %val(cstatus))
	Else
	   rstatus = fut_free_lun (outlun)
	   If (rstatus .NE. %loc(fut_normal)) then
	      Call lib$signal ('Error freeing outlun, status = ',%val(rstatus))
	   EndIf
	EndIf

	ffp_sc_close_sky = status

	Return
	End
