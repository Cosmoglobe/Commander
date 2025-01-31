	program fip_dust

c-------------------------------------------------------------------------------
c
c	Program FIP_DUST
c
c	This program reformats the FIRAS dust spectrum files into the
c	Project Dataset format.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  5 October 1994
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		fip_init_dust_report
c		fip_frequency_cut
c		fip_parse_dust
c		fip_read_nyquist
c		fip_reformat_dust
c		fip_update_dust_report
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c		str$trim
c
c	Include files:
c		fip_config_freq.txt
c		fip_invoc_dust.txt
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
	include '(fip_invoc_dust)'
	include '(fip_config_freq)'

	character * 14	current_gmt	!  GMT time of invocation

	integer * 4	fla_lun		!  input file lun
	integer * 4	fip_lun		!  output file lun
	integer * 4	io_stat		!  I/O return status
	integer * 4	lun_out /6/	!  terminal lun
	integer * 4	nspec		!  number of dust spectra
	integer * 4	parse_status	!  return status
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	fip_frequency_cut
	integer * 4	fip_init_dust_report
	integer * 4	fip_parse_dust
	integer * 4	fip_read_nyquist
	integer * 4	fip_reformat_dust
	integer * 4	fip_update_dust_report
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fla_dst'
	record /fla_dst/ in_rec
	dictionary 'fip_dst'
	record /fip_dst/ out_rec

	external	fut_error

	external	fip_failure
	external	fip_normal
	external	fip_repclose
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmsread
	external	fip_rmswrite
	external	fut_normal

c
c  Parse the command line.
c
	parse_status = fip_parse_dust (current_gmt)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .				     'FIRAS Facility FIP_DUST')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fip_init_dust_report (current_gmt, parse_status)
	   if (status .eq. %loc(fip_normal)) call lib$establish (fut_error)
	else
	   status = %loc(fip_normal)
	endif

	if ((status .eq. %loc(fip_normal))  .and.
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
	   endif
	endif


	if (status .eq. %loc(fip_normal)) then
C
C  Open the input and output spectrum files.
C

c
c  Find the file name.
c
	   fcc_fla_file = 'CSDR$FIRAS_IN:FLA_DST_' // fcc_scan_mode // '.' //
     .			   fcc_file_ext
	   fcc_fip_file = 'CSDR$FIRAS_OUT:FIP_DST_' // fcc_scan_mode // '.' //
     .			   fcc_file_ext
	   call str$trim (fcc_fla_file, fcc_fla_file, fcc_flalen)
	   call str$trim (fcc_fip_file, fcc_fip_file, fcc_fiplen)

c
c  Open the input file.
c
	   rstatus = fut_get_lun(fla_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=fla_lun, file=fcc_fla_file, status='old',
     .		 form='unformatted', access='sequential', readonly,
     .		 iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Open the output file.
c
	      rstatus = fut_get_lun(fip_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif

	      open (unit=fip_lun, file=fcc_fip_file, status='new',
     .		    form='unformatted', recordtype='fixed', recl=fcc_fip_size,
     .		    access='sequential', iostat=io_stat)

	      if (io_stat .eq. 0) then
C
C  Loop over the input records
C
	         nspec = 0
	         do while ((io_stat .eq. 0)  .and.
     .			   (status .eq. %loc(fip_normal)))
c
c  Read the input records.
c
	            read (fla_lun, iostat=io_stat) in_rec
	            if (io_stat .eq. 0) then
	               nspec = nspec + 1

c
c  Reformat the input record.
c
	               status = fip_reformat_dust (in_rec,out_rec)

	               if (status .eq. %loc(fip_normal)) then
c
c  Write the reformatted record to the output file.
c
	                  write (fip_lun, iostat=io_stat) out_rec
	                  if (io_stat .ne. 0) then
	                     status = %loc(fip_rmswrite)
	                     call lib$signal (fip_rmswrite, %val(2),
     .		   		     fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	                  endif
	               endif
	            elseif (io_stat .lt. 0) then
	               status = %loc(fip_normal)
	            else
	               status = %loc(fip_rmsread)
	               call lib$signal (fip_rmsread, %val(2),
     .		   	             fcc_fla_file(1:fcc_flalen), %val(io_stat))
	            endif  !  io_stat from read
	         enddo	!  while (io_stat


C
C  Close the spectrum files.
C

c
c  Close the output file.
c
	         close (unit=fip_lun, iostat=io_stat)
	         if (io_stat .ne. 0) then
	            status = %loc(fip_rmsclose)
	            call lib$signal (fip_rmsclose, %val(2),
     .	   			  fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	         endif

	         rstatus = fut_free_lun(fip_lun)
	         if (rstatus .ne. %loc(fut_normal)) then
	            call lib$signal (%val(rstatus))
	         endif

	      else
	         status = %loc(fip_rmsopen)
	         call lib$signal (fip_rmsopen, %val(2),
     .	   		       fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	      endif	!  (open status for output file

c
c  Close the input file.
c
	      close (unit=fla_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsclose)
	         call lib$signal (fip_rmsclose, %val(2),
     .	   		       fcc_fla_file(1:fcc_flalen), %val(io_stat))
	      endif

	      rstatus = fut_free_lun(fla_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif
	   else
	      status = %loc(fip_rmsopen)
	      call lib$signal (fip_rmsopen, %val(2),
     .	   		    fcc_fla_file(1:fcc_flalen), %val(io_stat))
	   endif	!  (open status for input file

	endif		!  (parse and report initialization


C
C  Close out the program.
C

c
c  Update the processing report with the number of spectra processed.
c
	if ((status .eq. %loc(fip_normal))  .and.
     .	    (fcc_report .eq. fac_present)) then
	   status = fip_update_dust_report (nspec)
	endif

c
c  Print out the number of spectra and pixels processed.
c
	type 20, nspec
 20	format (/, 4x, 'Number of dust spectra processed:  ', I5, /)

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
