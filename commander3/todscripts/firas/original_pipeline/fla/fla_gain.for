	program fla_gain

c-------------------------------------------------------------------------------
c
c	Program FLA_GAIN
c
c	This program computes the FIRAS gain function for a specific channel
c	and scan mode.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  30 June 1993
c                 SER 11413
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		fla_calc_responsivity
c		fla_init_gain_report
c		fla_parse_gain
c		fla_read_model
c		fla_read_reference
c		fla_update_gain_report
c		fla_write_gain
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		fla_invoc_gain.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_error)'
	include '(fut_params)'
	include '(fla_invoc_gain)'

	character * 14	current_gmt	!  GMT time of invocation

	integer * 4	io_stat		!  I/O return status
	integer * 4	lun_out /6/	!  terminal lun
	integer * 4	parse_status	!  return status
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	fla_calc_responsivity
	integer * 4	fla_init_gain_report
	integer * 4	fla_parse_gain
	integer * 4	fla_read_model
	integer * 4	fla_read_reference
	integer * 4	fla_update_gain_report
	integer * 4	fla_write_gain
	integer * 4	fut_free_lun

	external	fut_error

	external	fla_failure
	external	fla_normal
	external	fla_repclose
	external	fut_normal

c
c  Parse the command line.
c
	parse_status = fla_parse_gain (current_gmt)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .				     'FIRAS Facility FLA_GAIN')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fla_init_gain_report (current_gmt, parse_status)
	   if (status .eq. %loc(fla_normal)) call lib$establish (fut_error)
	else
	   status = %loc(fla_normal)
	endif


	if ((status .eq. %loc(fla_normal)) .and.
     .	    (parse_status .eq. %loc(fla_normal))) then
C
C  Read the model solution file.
C

c
c    The reference datsets.
c
	   status = fla_read_reference ()

c
c    The model solution file.
c
	   if (status .eq. %loc(fla_normal)) then
	      status = fla_read_model ()
	   endif


C
C  Compute the gain function.
C

c
c  Calculate the detector response.
c
	   if (status .eq. %loc(fla_normal)) then
	      status = fla_calc_responsivity ()
	   endif

c
c  Write the gain function.
c
	   if (status .eq. %loc(fla_normal)) then
	      status = fla_write_gain ()
	   endif

	endif	!  status from parse and report init


C
C  Close out the program.
C

c
c  Update the processing report
c
	if ((status .eq. %loc(fla_normal))  .and.
     .	    (fcc_report .eq. fac_present)) then
	   status = fla_update_gain_report ()
	endif

c
c  Signal the program completion status.
c
	if ((status .eq. %loc(fla_normal))  .and.
     .	    (parse_status .eq. %loc(fla_normal))) then
	   call lib$signal (fla_normal)
	else
	   call lib$signal (fla_failure)
	endif

c
c  Close the processing report file.
c
	if (fcc_report .eq. fac_present) then
	   close (unit=fut_report_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      call lib$signal (fla_repclose, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif
	   rstatus = fut_free_lun (fut_report_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	endif


	end
