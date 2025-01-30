	integer * 4 function  fcf_read_dvector ()

c-------------------------------------------------------------------------------
c
c	Function FCF_READ_DVECTOR
c
c	This function opens, reads, and closes the FIRAS D-Vector reference
c	dataset FEX_VAR_CCSS.VVV_XXXXXXXXXX from the directory CSDR$FIRAS_CAL 
c	(CC is channel, SS is scan mode, VVV is the version of the calibration
c	program that generated the solution, and XXXXXXXXXX is the model
c	solution label).  The file extension is specified by the command line
c	qualifier /MODEL_EXT.  The D-Vector is written into the include file
c	FCF_MODEL for use by FCF.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  25 October 1993
c		  SER 11397
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c		str$trim
c
c	Include files:
c		fcf_invoc.txt
c		fcf_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	var_lun		!  D-Vector file lun

	logical * 1	var_ver		!  D-Vector verification flag

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_var'
	record /fex_var/ fex_var

	external	fcf_varclose
	external	fcf_varinval
	external	fcf_varopen
	external	fcf_varread
	external	fcf_normal
	external	fut_normal

	status = %loc(fcf_normal)

C
C  Open, read, and close the D-Vector file.
C

c
c  Find the D-Vector filename.
c
	fcc_var_file = 'CSDR$FIRAS_CAL:FEX_VAR_' // fcc_scan_mode // '.' //
     .			fcc_model_ext
	call str$trim (fcc_var_file, fcc_var_file, fcc_varlen)

c
c  Open the D-Vector file.
c
	rstatus = fut_get_lun(var_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=var_lun, file=fcc_var_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_var_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the D-Vector.
c
	   read (var_lun, iostat=io_stat) fex_var
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_varread)
	      call lib$signal (fcf_varread, %val(2),
     .			       fcc_var_file(1:fcc_varlen), %val(io_stat))
	   endif

c
c  Close the D-Vector file.
c
	   close (unit=var_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_varclose)
	      call lib$signal (fcf_varclose, %val(2),
     .			       fcc_var_file(1:fcc_varlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(var_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fcf_varopen)
	   call lib$signal (fcf_varopen, %val(2),
     .			    fcc_var_file(1:fcc_varlen), %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fcf_normal)) then
C
C  Process the D-Vector.
C

c
c  Verify that the correct D-Vector has been read.
c
	   var_ver = .true.
	   if (fcc_chan .ne. fex_var.channel) var_ver = .false.
	   if (fcc_length .ne. fex_var.scan_length) var_ver = .false.
	   if (fcc_speed .ne. fex_var.scan_speed) var_ver = .false.
	   if (var_ver .eq. .false.) then
	      status = %loc(fcf_varinval)
	      call lib$signal (fcf_varinval, %val(1),
     .			       fcc_var_file(1:fcc_varlen))
	   endif

	endif

	if (status .eq. %loc(fcf_normal)) then
c
c  Insert the D-Vector into the common block.
c
	   call lib$movc5 (0,,0,2056,dvector)
	   do j = 1,257
	      dvector(j) = fex_var.variances(j)
	   enddo 
	endif


	fcf_read_dvector = status

	return
	end
