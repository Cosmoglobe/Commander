	program fip_cvs

c-------------------------------------------------------------------------------
c
c	Program FIP_CVS
c
c	This program reformats the C-Vector files into the
c	Project Dataset Release format.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  17 November 1994
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		fip_init_cvs_report
c		fip_frequency_cut
c		fip_parse_cvs
c		fip_read_nyquist
c		fip_reformat_cvs
c		fip_update_cvs_report
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		fip_config_freq.txt
c		fip_invoc_cvs.txt
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
	include '(fip_invoc_cvs)'
	include '(fip_config_freq)'

	character * 14	current_gmt	!  GMT time of invocation

	integer * 4	io_stat		!  I/O return status
	integer * 4	lun_out /6/	!  terminal lun
	integer * 4	parse_status	!  return status
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	fip_frequency_cut
	integer * 4	fip_init_cvs_report
	integer * 4	fip_parse_cvs
	integer * 4	fip_read_nyquist
	integer * 4	fip_reformat_cvs
	integer * 4	fip_update_cvs_report
	integer * 4	fut_free_lun

	external	fut_error

	external	fip_failure
	external	fip_normal
	external	fip_repclose
	external	fut_normal

c
c  Parse the command line.
c
	parse_status = fip_parse_cvs (current_gmt)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .				     'FIRAS Facility FIP_CVS')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fip_init_cvs_report (current_gmt, parse_status)
	   if (status .eq. %loc(fip_normal)) call lib$establish (fut_error)
	else
	   status = %loc(fip_normal)
	endif


	if ((status .eq. %loc(fip_normal)) .and.
     .	    (parse_status .eq. %loc(fip_normal))) then
c
c  Get the Nyquist frequency.
c
	   status = fip_read_nyquist ()

c
c  Get the frequency cut.
c
	   if (status .eq. %loc(fip_normal)) then
	      status = fip_frequency_cut (fnyq_icm, fnyq_hz)

c
c  Process the C-Vector.
c
	      if (status .eq. %loc(fip_normal)) then
	         status = fip_reformat_cvs ()
	      endif

	   endif	!  status from read nyquist

	endif	!  status from parse and report init


c
c  Update the processing report
c
	if ((status .eq. %loc(fip_normal))  .and.
     .	    (fcc_report .eq. fac_present)) then
	   status = fip_update_cvs_report ()
	endif

c
c  Signal the program completion status.
c
	if ((status .eq. %loc(fip_normal))  .and.
     .	    (parse_status .eq. %loc(fip_normal))) then
	   call lib$signal (fip_normal)
	else
	   call lib$signal (fip_failure)
	endif

c
c  Close the processing report file.
c
	if (fcc_report .eq. fac_present) then
	   close (unit=fut_report_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      call lib$signal (fip_repclose, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif
	   rstatus = fut_free_lun (fut_report_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	endif


	end
