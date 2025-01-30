	Integer*4  Function  FFP_SC_Update_Report ( data_type, nrec, npix,
	1                                           input, filexts, fnum )
c------------------------------------------------------------------------------
c
c	Function FFP_SC_UPDATE_REPORT
c
c	This function writes the input and output file names, the data
c	selection cuts,  and the number of spectra and pixels processed to the
c	processing report.
c
c	Passed Parameters:
c	   Input:
c	      data_type		character*3	SKY or CAL
c	      nrec		integer*4	number of records processed
c	      npix		integer*4	number of pixels processed
c	      input		character*12	Input filename base
c	      filexts		character*20(fac_max_num) Input file extensions
c	      fnum		integer*2	Number of input files
c
c	   Output:
c	      none
c
c       Author:  S. Brodd, HSTX, 3/21/96
c
c	Subroutines called:
c		lib$signal
c		str$trim
c
c	Include files:
c		ffp_frequency.txt
c		ffp_invoc_sky.txt
c		fut_error.txt
c		fut_params.txt
c
c-----------------------------------------------------------------------------

	Implicit none

	Include '(fut_error)'
	Include '(fut_params)'
	Include '(ffp_invoc_sky)'
	Include '(ffp_frequency)'

c Passed Parameters

	Character*3	data_type	!  SKY or CAL
	Integer*4	npix		!  number of pixels containing spectra
	Integer*4	nrec		!  number of records processed
	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files

c Local

	Integer*4	io_stat		!  I/O return status
	Integer*2	fil		!  file counter
	Character*33	infile			!  Input file name
	Character*1	blanks(33) /33*' '/
	Character*33	blankname		!  Blank input file name
	Equivalence	(blanks, blankname)

c External

	External	ffp_normal
	External	ffp_repwrite

c Begin
c
c  Write the output file.
c
	Write (fut_report_lun,30,iostat=io_stat) fcc_outfile(1:fcc_outlen)
  30	Format (/,1x, 'Output File:  ', a,/)
	If (fcc_destriped .EQ. fac_present .AND. data_type .EQ. 'SKY') Then
	   Write (fut_report_lun,40,iostat=io_stat)
  40       Format (4x, 'Spectra were destriped.')
	EndIf
c
c  Write out the data selection cuts:
c
	If (fcc_galexc .eq. fac_present) Then
	   Write (fut_report_lun,50,iostat=io_stat) fcc_glat
  50	   Format (4x,'Galactic Latitude Cutoff:    ', 11x, F6.2,
	1          '     degrees')
	Else
	   Write (fut_report_lun,60,iostat=io_stat)
  60       Format (4x, 'All Galactic Latitudes Included.')
	EndIf
	If (fcc_freq .EQ. fac_present) Then
	   Write (fut_report_lun,80,iostat=io_stat) fcc_jlo, fcc_jhi
  80	   Format (4x, 'Frequency Indices:           ', 11x, I3, ' - ',I3)
	   Write (fut_report_lun,90,iostat=io_stat) fcc_nu0, fcc_dnu, fcc_nfreq
  90	   Format (4x, 'Initial Optical Frequency:   ', 11x, F7.3, '    GHz',/,
	1          4x, 'Optical Frequency Interval:  ', 11x, F7.3, '    GHz',/,
	2          4x, 'Number of Frequency Points:  ', 11x, I3)
	Else
	   Write (fut_report_lun,100,iostat=io_stat)
 100	   Format (4x, 'All Frequencies Included.')
	EndIf

	If (data_type .EQ. 'SKY') Then
c
c  Write the number of spectra processed:
c
	   Write (fut_report_lun,110,iostat=io_stat) nrec, npix
 110	   Format (4x, 'Number of spectra processed:      ', I5, /,
	1          4x, 'Number of pixels containing spectra:  ', I5, /)
	Else
	   Write (fut_report_lun,120,iostat=io_stat) nrec
 120	   Format (4x, 'Number of records processed:      ', I5)
	EndIf
c
c  Check the status of the writes:
c
	If (io_stat .NE. 0) Then
	   ffp_sc_update_report = %loc(ffp_repwrite)
	   Call lib$signal (ffp_repwrite, %val(2),
	1                   fcc_report_file(1:fcc_replen), %val(io_stat))
	Else
	   ffp_sc_update_report = %loc(ffp_normal)
	EndIf

	Return
	End
