	Integer*4  Function  FIP_SC_Open_Sky ( input, filexts, fnum, archin,
	1                                      archout, inlun, outlun )

c------------------------------------------------------------------------------
c
c	Function FIP_SC_OPEN_SKY
c
c	This function opens the FIRAS input spectrum skymaps and the ADB output
c	spectrum skymap.
c
c	Author:  Larry P. Rosen, Hughes STX, 2 December 1994
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
c		fip_invoc_sky.txt
c		fut_params.txt
c
c------------------------------------------------------------------------------
	Implicit None

	Include '(fut_params)'
	Include '(fip_invoc_sky)'

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
	Integer*2	fic_alen		!  Length of output fic record
	parameter	(fic_alen=1088)		!  in longwords.
c					! Note that length of output fcf record
c					! is fcc_alen in fip_invoc_sky.txt
c Functions

	Integer*4	fut_get_lun
	Integer*4	csa_open_skymap
	Integer*4	csa_field_offset_values

c Externals

	External	csa_open_skymap
	External	csa_normal
	External	fip_csafld
	External	fip_csaopen
	External	fip_normal
	External	fip_lunerr
	External	fut_normal
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  Initialize the open.
C
	FIP_SC_Open_Sky = %loc(fip_normal)
c
c  Get the logical unit numbers for the files.
c
	rstatus = fut_get_lun (outlun)
	If (rstatus .NE. %loc(fut_normal)) Then
	   FIP_SC_Open_Sky = rstatus
	   Call lib$signal (fip_lunerr, %val(1), %loc(rstatus))
	   Return
	EndIf
	Do i=1,fnum
	   rstatus = fut_get_lun (inlun(i))
	   If (rstatus .NE. %loc(fut_normal)) Then
	      FIP_SC_Open_Sky = rstatus
	      Call lib$signal (fip_lunerr, %val(1), %loc(rstatus))
	      Return
	   EndIf
	EndDo
C
	If (fcc_coadd) then
	   len = fic_alen			! fic records
	Else
	   len = fcc_alen			! fcf records
	EndIf
C
C  Open the files.
C
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
	      FIP_SC_Open_Sky = %loc(fip_csafld)
	      Call lib$signal ( fip_csafld, %val(2),
	1                       fcc_outfile(1:fcc_outlen), %val(cstatus) )
	   EndIf
	Else
 	   FIP_SC_Open_Sky = %loc(fip_csaopen)
	   Call lib$signal ( fip_csaopen, %val(2), fcc_outfile(1:fcc_outlen),
	1                    %val(io_stat) )
	EndIf
c
c  Open the input FIRAS skymap files.
c
	i = 1
	Do While (i .LE. fnum .AND. FIP_SC_Open_Sky .EQ. %loc(fip_normal))
	   infile = blankname
	   infile = archin // ':' // input // '.' // filexts (i)
	   Open ( unit=inlun(i), file=infile, iostat=io_stat, status='old',
	1         form='unformatted', recordtype='fixed', readonly,
	2         useropen=csa_open_skymap )
	   If (io_stat .EQ. 0) Then
	      cstatus = csa_field_offset_values (fcc_pix_offset,-1,-1,inlun(i))
	      If (cstatus .NE. %loc(csa_normal)) Then
	         FIP_SC_Open_Sky = %loc(fip_csafld)
	         Call lib$signal ( fip_csafld, %val(2), infile, %val(cstatus) )
	      EndIf
	   Else
	      FIP_SC_Open_Sky = %loc(fip_csaopen)
	      Call lib$signal ( fip_csaopen, %val(2), infile, %val(io_stat) )
	   EndIf
	   i = i + 1
	EndDo
	Return
	End
