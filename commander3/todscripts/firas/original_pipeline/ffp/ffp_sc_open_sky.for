	Integer*4  Function  FFP_SC_Open_Sky ( input, filexts, fnum, archin,
	1                                      archout, inlun, outlun )
c------------------------------------------------------------------------------
c
c	Function FFP_SC_OPEN_SKY
c
c	This function opens the FIRAS input spectrum skymaps and the ADB output
c	spectrum skymap.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		input	character*12			! Input filename base
c		filexts	Character*20 (fac_max_num)	! Input file extensions
c		fnum	integer*2		! Number of input skymaps
c		archin  character*13		! Input archive name
c		archout character*14		! Output archive name
c
c	Output:
c		inlun	integer*4 (fac_max_num)	Input skymap lun's
c		outlun	integer*4		Output skymap lun
c
c	Subroutines called:
c		csa_field_offset_values
c		csa_open_skymap
c		fut_get_lun
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

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			!  Number of input maps
	Character*13	archin			!  Input archive name
	Character*14	archout			!  Output archive name
	Integer*4	inlun (fac_max_num)	!  Input skymap luns
	Integer*4	outlun			!  Output skymap lun

c Local

	Integer*4	rstatus			!  return status
	Integer*2	i			!  counter
	Integer*4	io_stat			!  I/O return status
	Integer*4	cstatus			!  CT or CSA return status
	Character*47	infile			!  Input file name
	Character*63	outfile
	Character*1	blanks(47) /47*' '/
	Character*47	blankname		!  Blank input file name
	Equivalence	(blanks, blankname)
	Integer*2	len		!  Length of output record in longwords

c Functions

	Integer*4	fut_get_lun
	Integer*4	csa_open_skymap
	Integer*4	csa_field_offset_values

c Externals

	External	csa_open_skymap
	External	csa_normal
	External	ffp_csafld
	External	ffp_csaopen
	External	ffp_normal
	External	ffp_lunerr
	External	fut_normal
C
C  Initialize the open.
C
	FFP_SC_Open_Sky = %loc(ffp_normal)
c
c  Get the logical unit numbers for the files.
c
	rstatus = fut_get_lun (outlun)
	If (rstatus .NE. %loc(fut_normal)) Then
	   FFP_SC_Open_Sky = rstatus
	   Call lib$signal (ffp_lunerr, %val(1), %loc(rstatus))
	   Return
	EndIf
	Do i=1,fnum
	   rstatus = fut_get_lun (inlun(i))
	   If (rstatus .NE. %loc(fut_normal)) Then
	      FFP_SC_Open_Sky = rstatus
	      Call lib$signal (ffp_lunerr, %val(1), %loc(rstatus))
	      Return
	   EndIf
	EndDo
C
	If (fcc_coadd) then
	   len = fcc_ilen			! fil records
	Else
	   len = fcc_slen			! fsl records
	EndIf
c
c  Open the output ADB skymap file.
c
	outfile = archout // ':' // fcc_outfile (1:fcc_outlen)
	Open ( unit=outlun, file=outfile, iostat=io_stat, status='new',
	1      form='unformatted', recordtype='fixed', recl=len,
	2      useropen=csa_open_skymap )

	If (io_stat .EQ. 0) Then
c
c  Set up the skymap offset values.
c
	   cstatus = csa_field_offset_values (fcc_pix_offset, -1, -1, outlun)
	   If (cstatus .NE. %loc(csa_normal)) Then
	      FFP_SC_Open_Sky = %loc(ffp_csafld)
	      Call lib$signal ( ffp_csafld, %val(2),
	1                       fcc_outfile(1:fcc_outlen), %val(cstatus) )
	   EndIf
	Else
 	   FFP_SC_Open_Sky = %loc(ffp_csaopen)
	   Call lib$signal ( ffp_csaopen, %val(2), fcc_outfile(1:fcc_outlen),
	1                    %val(io_stat) )
	EndIf
c
c  Open the input FIRAS skymap files.
c
	i = 1
	Do While (i .LE. fnum .AND. FFP_SC_Open_Sky .EQ. %loc(ffp_normal))
	   infile = blankname
	   infile = archin // ':' // input // '.' // filexts (i)
	   Open ( unit=inlun(i), file=infile, iostat=io_stat, status='old',
	1         form='unformatted', recordtype='fixed', readonly,
	2         useropen=csa_open_skymap )
	   If (io_stat .EQ. 0) Then
	      cstatus = csa_field_offset_values (fcc_pix_offset,-1,-1,inlun(i))
	      If (cstatus .NE. %loc(csa_normal)) Then
	         FFP_SC_Open_Sky = %loc(ffp_csafld)
	         Call lib$signal ( ffp_csafld, %val(2), infile, %val(cstatus) )
	      EndIf
	   Else
	      FFP_SC_Open_Sky = %loc(ffp_csaopen)
	      Call lib$signal ( ffp_csaopen, %val(2), infile, %val(io_stat) )
	   EndIf
	   i = i + 1
	EndDo

	Return
	End
