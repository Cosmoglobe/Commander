	integer * 4 function  fla_read_model ()

c-------------------------------------------------------------------------------
c
c	Function FLA_READ_MODEL
c
c	This function opens, reads, and closes the FIRAS calibration model
c	solution FEX_MOD_CCSS.VVV_XXXXXXXXXX from the directory CSDR$FIRAS_CAL 
c	(CC is channel, SS is scan mode, VVV is the version of the calibration
c	program that generated the solution, and XXXXXXXXXX is the model
c	solution label).  The file extension is specified by the command line
c	qualifier /MODEL_EXT.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  29 June 1993
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
c		fla_invoc_gain.txt
c		fla_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fla_invoc_gain)'
	include '(fla_model)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	mod_lun		!  model solution file lun
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	fla_normal
	external	fla_rmsclose
	external	fla_rmsopen
	external	fla_rmsread
	external	fut_normal

	status = %loc(fla_normal)

C
C  Open, read, and close the calibration model solution file.
C

c
c  Find the model solution file name.
c
	fcc_fex_file = 'CSDR$FIRAS_CAL:FEX_MOD_' // fcc_scan_mode // '.' //
     .			  fcc_model_ext
	call str$trim (fcc_fex_file, fcc_fex_file, fcc_fexlen)

c
c  Open the model solution file.
c
	rstatus = fut_get_lun(mod_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=mod_lun, file=fcc_fex_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_fex_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the model solution.
c
	   read (mod_lun, iostat=io_stat) fex_model
	   if (io_stat .ne. 0) then
	      status = %loc(fla_rmsread)
	      call lib$signal (fla_rmsread, %val(2),
     .			       fcc_fex_file(1:fcc_fexlen), %val(io_stat))
	   endif

c
c  Close the model solution file.
c
	   close (unit=mod_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fla_rmsclose)
	      call lib$signal (fla_rmsclose, %val(2),
     .			       fcc_fex_file(1:fcc_fexlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(mod_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fla_rmsopen)
	   call lib$signal (fla_rmsopen, %val(2),
     .			    fcc_fex_file(1:fcc_fexlen), %val(io_stat))
	endif	!  (open status


	fla_read_model = status

	return
	end
