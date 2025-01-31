	Program FFI

c-------------------------------------------------------------------------------
c
c	FFI_FISH_INPUT
c
c	This program drives the FISH spectrum generation routines.  The program
c	reads coadded calibration IFGs from the Cobetrieve archives and writes
c	FISH-format voltage spectra to an RMS file.  The program runs for a
c	single channel and scan mode.  The program can read coadds from two
c	input archives and write a hybrid output to the RMS file.
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  9 March 1993
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cct_close_config
c		ct_init
c		cut_display_banner
c		cut_register_version
c		ffi_get_reference
c		ffi_initialize_report
c		ffi_open_coadd
c		ffi_parse
c		ffi_produce_spectra
c		ffi_read_coadd
c		ffi_read_hybrid_coadd
c		ffi_update_report
c		ffi_write_spec
c		fut_free_lun
c		iwkin
c		lib$establish
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		ffi_config.txt
c		ffi_invoc.txt
c		ffi_spec.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_error)'
	include '(fut_params)'
	include '(ffi_config)'
	include '(ffi_invoc)'
	include '(ffi_spec)'

	character * 79	cmd_line(3)		!  command line invocation
	character * 14	current_gmt		!  GMT time of invocation

	integer *  2	ct_stat(20)		!  Cobetrieve return status

	integer *  4	cmdlen(3)		!  length of command lines in
						!    invocation
	integer *  4	cstatus			!  Close_Config return status
	integer *  4	ct_lun1			!  coadd file lun
	integer *  4	ct_lun2			!  coadd file lun
	integer *  4	io_stat			!  I/O return status
	integer *  4	j			!  a counter
	integer *  4	lun_out /6/		!  terminal lun
	integer *  4	ncmd			!  number of command lines in
						!    invocation
	integer *  4	num			!  number of current coadds read
	integer *  4	parse_status		!  command line parse return
						!    status
	integer *  4	read_status		!  coadd read return status
	integer *  4	rstatus			!  return status
	integer *  4	status			!  return status
	integer *  4	vs_lun			!  RMS file lun

	logical *  1	ref_open_dir		!  reference dataset open flag
	logical *  1	ref_open_seq		!  reference dataset open flag

	real	* 4	rwksp(2096)		!  IMSL work space
	common /worksp/ rwksp

	integer * 4	cct_close_config
	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	ffi_get_reference
	integer * 4	ffi_initialize_report
	integer * 4	ffi_open_coadd
	integer * 4	ffi_parse
	integer * 4	ffi_produce_spectra
	integer * 4	ffi_read_coadd
	integer * 4	ffi_read_hybrid_coadd
	integer * 4	ffi_update_report
	integer * 4	ffi_write_spec
	integer * 4	fut_free_lun

	dictionary 'fic_sky'
	record /fic_sky/ coadd_recs(fac_max_coad)
	record /fish_spec/ fish_recs(fac_max_coad)

	external	fut_error

	external	cct_normal
	external	ffi_chanerr
	external	ffi_cfgdirclose
	external	ffi_cfgseqclose
	external	ffi_ctinit
	external	ffi_eof
	external	ffi_failure
	external	ffi_normal
	external	ffi_repclose
	external	ffi_rmsclose
	external	fut_normal

C
C  Initialize the program.
C

c
c  Initialize the IMSL workspace.
c
	call iwkin(2096)

c
c  Parse the command line.
c
	parse_status = ffi_parse (current_gmt, ncmd, cmd_line, cmdlen)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .					   'FIRAS Facility FFI_Fish_Input')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = ffi_initialize_report (current_gmt, ncmd, cmd_line, cmdlen,
     .					   parse_status)
	   if (status .eq. %loc(ffi_normal)) call lib$establish (fut_error)
	else
	   status = %loc(ffi_normal)
	endif

	if ((parse_status .eq. %loc(ffi_normal))  .and.
     .	    (status .eq. %loc(ffi_normal))) then
c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)

	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Get the reference datasets.
c
	      status = ffi_get_reference (ref_open_dir, ref_open_seq)


	      if (status .eq. %loc(ffi_normal)) then
C
C  Process the spectra for the specified channel and scan mode.
C

c
c  Open the coadd file.
c
	         status = ffi_open_coadd (ct_lun1, ct_lun2, vs_lun)

	         if (status .eq. %loc(ffi_normal)) then
c
c  Initialize the coadd read.
c
	            fcc_nspec    = 0
	            fcc_nspec1   = 0
	            fcc_nspec2   = 0
	            num          = 0
	            read_status  = %loc(ffi_normal)

	            type 20, fac_channel_ids(fcc_chan),
     .			     fac_scan_mode_ids(fcc_smode)
 20	            format (x, 'Processing spectra for channel ', a,
     .			       ', scan mode ', a, /)

c
c  Read in a set of data for a particular scan mode.
c
	            do while (read_status .eq. %loc(ffi_normal))

	               if (fcc_hybrid .eq. fac_present) then
	                  read_status = ffi_read_hybrid_coadd (ct_lun1, ct_lun2,
     .							       num, coadd_recs)
	               else
	                  read_status = ffi_read_coadd (ct_lun1, num,
     .							coadd_recs)
	               endif

	               if ((read_status .eq. %loc(ffi_normal))  .or.
     .                     ((read_status .eq. %loc(ffi_eof)) .and.
     .			    (num .gt. 0))) then

	                  status = %loc(ffi_normal)
	                  j = 0
	                  fcc_nspec = fcc_nspec + num

	                  do while ((status .eq. %loc(ffi_normal))  .and.
     .                              (j .lt. num))
	                     j = j + 1
c
c  Produce the voltage spectra.
c
	                     status = ffi_produce_spectra (coadd_recs(j),
     .							   fish_recs(j))

	                  enddo		!  (while j .lt. num

c
c  Write the voltage spectra to the FISH-format RMS file.
c
	                  if (status .eq. %loc(ffi_normal)) then

	                     status = ffi_write_spec (vs_lun, num, fish_recs)

	                     if (read_status .eq. %loc(ffi_eof)) then
	                        close (unit=vs_lun, iostat=io_stat)
	                        if (io_stat .ne. 0) then
	                           status = %loc(ffi_rmsclose)
	                           call lib$signal (ffi_rmsclose, %val(2),
     .						    fcc_outfile(1:fcc_outlen),
     .						    %val(io_stat))
	                        endif
	                        rstatus = fut_free_lun (vs_lun)
	                        if (rstatus .ne. %loc(fut_normal)) then
	                           call lib$signal (%val(rstatus))
	                        endif
	                     endif

	                  endif		!  (status from processing spectra

c
c  Close the FISH-format files if no records were read the last time.
c
	               elseif ((read_status .eq. %loc(ffi_eof))  .and.
     .			       (num .eq. 0)) then
	                  close (unit=vs_lun, iostat=io_stat)
	                  if (io_stat .ne. 0) then
	                     status = %loc(ffi_rmsclose)
	                     call lib$signal (ffi_rmsclose, %val(2),
     .					      fcc_outfile(1:fcc_outlen),
     .					      %val(io_stat))
	                  endif
	                  rstatus = fut_free_lun (vs_lun)
	                  if (rstatus .ne. %loc(fut_normal)) then
	                     call lib$signal (%val(rstatus))
	                  endif
	               else
	                  status = %loc(ffi_chanerr)
	                  call lib$signal (ffi_chanerr)
	               endif	!  (read_status from ffi_read_coadd

	            enddo	!  (do while read_status


C
C  Close out the program.
C

c
c  Update the processing report with the number of spectra processed.
c
	            if ((status .eq. %loc(ffi_normal))  .and.
     .	                (fcc_report .eq. fac_present)) then
	               status = ffi_update_report ()
	            endif

c
c  Write out number of spectra processed.
c
	            if (status .eq. %loc(ffi_normal)) then
	               if (fcc_hybrid .eq. fac_present) then
	                  write (lun_out,30,iostat=io_stat) fcc_nspec,
     .							  fcc_nspec1, fcc_nspec2
	               else
	                  write (lun_out,40,iostat=io_stat) fcc_nspec
	               endif
	            endif
  30	format (/, 4x, 'Number of calibration spectra processed:  ', I5,
     .		/, 4x, '    Number from input 1:                  ', I5,
     .		/, 4x, '    Number from input 2:                  ', I5, //)
  40	format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)

	         endif	!	(status from ffi_open_coadd

	      endif	!  	(status from ffi_get_reference

c
c  Close the configuration files.
c
	      if (ref_open_seq) then
	         cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(ffi_cfgseqclose)
	            call lib$signal (ffi_cfgseqclose, %val(1), %val(cstatus))
	         endif
	      endif

	      if (ref_open_dir) then
	         cstatus = cct_close_config (ndset_dir, dir_lun, dir_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(ffi_cfgdirclose)
	            call lib$signal (ffi_cfgdirclose, %val(1), %val(cstatus))
	         endif
	      endif

	   else
	      status = %loc(ffi_ctinit)
	      call lib$signal (ffi_ctinit, %val(1), %val(ct_stat(1)))
	   endif		!  (ct_stat from ct_init

	endif			!  (status from ffi_parse and ffi_init_report

c
c  Signal the program completion status.
c
	if ((status .eq. %loc(ffi_normal))  .or.
     .	    (status .eq. %loc(ffi_chanerr))) then
	   call lib$signal (ffi_normal)
	else
	   call lib$signal (ffi_failure)
	endif

c
c  Close the processing report file.
c
	if (fcc_report .eq. fac_present) then
	   close (unit=fut_report_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      call lib$signal (ffi_repclose, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif
	   rstatus = fut_free_lun (fut_report_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	endif


	end
