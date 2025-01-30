	Program FCF

c-------------------------------------------------------------------------------
c
c	FCF_CALIBRATE_FIRAS
c
c	This program drives the Firas calibration routines.  The program reads
c	coadded IFGs from and writes calibrated spectra to the Cobetrieve
c	archives.  The program reads in the Firas calibration model solutions
c	from a Cobetrieve reference archive.
c
c	For a given channel and scan mode, FCF reads in coadded IFGs and
c	transforms them into voltage spectra.  It then applies the calibration
c	model to both the spectra and the spectrum variances.  FCF calibrates
c	both sky and calibration data.
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  23 September 1992
c                 SER 8292
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cct_close_config
c		cct_open_config
c		ct_init
c		cut_display_banner
c		cut_register_version
c		fcf_calibrate_spectra
c		fcf_display_spectra
c		fcf_initialize_report
c		fcf_open_coadd
c		fcf_parse
c		fcf_produce_spectra
c		fcf_read_coadd
c		fcf_update_report
c		fcf_write_cal_spec
c		fcf_write_sky_spec
c		fcf_write_temp_spec
c		fut_free_lun
c		iwkin
c		lib$establish
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fcf_config.txt
c		fcf_invoc.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Add DVECTOR command line switch and optional use of reference dataset
c	FEX_VAR in weighting the autophase corrector computation.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_error)'
	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_invoc)'

	character * 79	cmd_line(11)		!  command line invocation
	character * 14	current_gmt		!  GMT time of invocation

	complex	* 16	cspec(257)		!  calibrated spectrum
	complex	* 16	vspec(257)		!  voltage spectrum

	integer *  2	ct_stat(20)		!  Cobetrieve return status

	integer *  4	cmdlen(11)		!  length of command lines in
						!    invocation
	integer *  4	cs_lun			!  RMS file lun
	integer *  4	cstatus			!  Get_Config return status
	integer *  4	ct_lun			!  coadd file lun
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
	integer *  4	tstatus			!  return status
	integer *  4	vs_lun			!  RMS file lun

	logical *  1	ref_open_dir		!  reference dataset open flag
	logical *  1	ref_open_seq		!  reference dataset open flag

	real	*  8	cvar(3,257)		!  calibrated spectrum variances
	real	*  8	vvar(3,257)		!  voltage spectrum variances

	real	*  4	rwksp(2096)		!  IMSL work space
	common /worksp/ rwksp

	integer *  4	cct_close_config
	integer *  4	cct_open_config
	integer *  4	cut_display_banner
	integer *  4	cut_register_version
	integer *  4	fcf_calibrate_spectra
	integer *  4	fcf_display_spectra
	integer *  4	fcf_initialize_report
	integer *  4	fcf_open_coadd
	integer *  4	fcf_parse
	integer *  4	fcf_produce_spectra
	integer *  4	fcf_read_coadd
	integer *  4	fcf_update_report
	integer *  4	fcf_write_cal_spec
	integer *  4	fcf_write_sky_spec
	integer *  4	fcf_write_temp_spec
	integer *  4	fut_free_lun

	dictionary 'fic_sky'
	record /fic_sky/ coadd_recs(fcc_max_coadd)
	dictionary 'fcf_sky'
	record /fcf_sky/ cspec_recs(fcc_max_coadd)
	record /fcf_sky/ vspec_recs(fcc_max_coadd)

	external	fut_error

	external	cct_normal
	external	fcf_chanerr
	external	fcf_cfgdirclose
	external	fcf_cfgdiropen
	external	fcf_cfgseqclose
	external	fcf_cfgseqopen
	external	fcf_ctinit
	external	fcf_eof
	external	fcf_failure
	external	fcf_normal
	external	fcf_repclose
	external	fcf_tempclose
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
	parse_status = fcf_parse (current_gmt, ncmd, cmd_line, cmdlen)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .					   'FIRAS Facility FCF_Calibrate_Firas')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fcf_initialize_report (current_gmt, ncmd, cmd_line, cmdlen,
     .					   parse_status)
	   if (status .eq. %loc(fcf_normal)) call lib$establish (fut_error)
	else
	   status = %loc(fcf_normal)
	endif

	if ((parse_status .eq. %loc(fcf_normal))  .and.
     .	    (status .eq. %loc(fcf_normal))) then
c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)

	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Open the configuration files for sequential, direct, and model access.
c
	      ref_open_dir = .false.
	      ref_open_seq = .false.
	      call ct_gmt_to_binary(ref_gmt_start,ref_start)
	      call ct_gmt_to_binary(ref_gmt_stop,ref_stop)

	      cstatus = cct_open_config (ref_start, ref_stop, ndset_tod,
     .                                   dset_tod, size_tod, seq_access,
     .                                   ncache, tod_lun, tod_index,
     .                                   tod_stat, ref_count)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fcf_cfgseqopen)
	         call lib$signal (fcf_cfgseqopen, %val(1), %val(cstatus))
	      else
	         ref_open_seq = .true.

	         cstatus = cct_open_config (ref_start, ref_stop, ndset_dir,
     .                                      dset_dir, size_dir, dir_access, 
     .                                      ncache, dir_lun, dir_index,
     .                                      dir_stat, ref_count)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fcf_cfgdiropen)
	            call lib$signal (fcf_cfgdiropen, %val(1), %val(cstatus))
	         else
	            ref_open_dir = .true.
	         endif
	      endif


	      if (status .eq. %loc(fcf_normal)) then
C
C  Process the spectra for the specified channel and scan mode.
C

c
c  Open the coadd file.
c
	         status = fcf_open_coadd (ct_lun, vs_lun, cs_lun)

	         if (status .eq. %loc(fcf_normal)) then
c
c  Initialize the coadd read.
c
	            fcc_nspec      = 0
	            num            = 0
	            fcc_more_plots = .true.
	            read_status    = %loc(fcf_normal)

	            type 20, fac_channel_ids(fcc_chan),
     .			     fac_scan_mode_ids(fcc_smode)
 20	            format (x, 'Processing spectra for channel ', a,
     .			       ', scan mode ', a, /)

c
c  Read in a set of data for a particular scan mode.
c
	            do while (read_status .eq. %loc(fcf_normal))

	               read_status = fcf_read_coadd (ct_lun, num, coadd_recs)

	               if ((read_status .eq. %loc(fcf_normal))  .or.
     .                     ((read_status .eq. %loc(fcf_eof)) .and.
     .			    (num .gt. 0))) then

	                  status = %loc(fcf_normal)
	                  j = 0
	                  fcc_nspec = fcc_nspec + num

	                  do while ((status .eq. %loc(fcf_normal))  .and.
     .                              (j .lt. num))
	                     j = j + 1
c
c  Produce the voltage spectra.
c
	                     status = fcf_produce_spectra (coadd_recs(j), vspec,
     .							   vvar, vspec_recs(j),
     .							   cspec_recs(j))

c
c  Calibrate the spectra.
c
	                     if (fcc_calibrate .eq. fac_present) then
	                        status = fcf_calibrate_spectra (vspec, vvar,
     .								cspec, cvar,
     .								cspec_recs(j))
	                     endif

c
c  Plot the spectra.
c
	                     fcc_jump = fcc_jump - 1

	                     if ((status .eq. %loc(fcf_normal))  .and.
     .			         (fcc_plot .eq. fac_present)  .and.
     .			         (fcc_more_plots .eq. .true.)  .and.
     .			         (fcc_jump .le. 0)) then

	                        tstatus = fcf_display_spectra ()

	                        if (tstatus .eq. fac_not_present) then
	                           fcc_more_plots = .false.
	                        endif

	                     endif

	                  enddo		!  (while j .lt. num


c
c  Write the voltage spectra to a temporary RMS file.
c
	                  if ((status .eq. %loc(fcf_normal))  .and.
     .                        (fcc_write_vs .eq. fac_present)) then

	                     status = fcf_write_temp_spec (1, vs_lun, num,
     .							   vspec_recs)

	                     if (read_status .eq. %loc(fcf_eof)) then
	                        close (unit=vs_lun, iostat=io_stat)
	                        if (io_stat .ne. 0) then
	                           status = %loc(fcf_tempclose)
	                           call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_vs(1:fcc_toutlenvs),
     .					   %val(io_stat))
	                        endif
	                     endif

	                  endif		!  (status from processing spectra

c
c  Write the calibrated or differential spectra to a temporary RMS file.
c
	                  if (status .eq. %loc(fcf_normal)) then
	                     if (fcc_write_ds .eq. fac_present) then
	                        status = fcf_write_temp_spec (2, cs_lun, num,
     .							      cspec_recs)
	                     elseif (fcc_write_cs .eq. fac_present) then
	                        status = fcf_write_temp_spec (3, cs_lun, num,
     .							      cspec_recs)
	                     endif
	                     if (read_status .eq. %loc(fcf_eof)) then
	                        close (unit=cs_lun, iostat=io_stat)
	                        if (io_stat .ne. 0) then
	                           status = %loc(fcf_tempclose)
	                           if (fcc_write_ds .eq. fac_present) then
	                              call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_ds(1:fcc_toutlends),
     .					   %val(io_stat))
	                           elseif (fcc_write_cs .eq. fac_present) then
	                              call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_cs(1:fcc_toutlencs),
     .					   %val(io_stat))
	                           endif
	                        endif
	                     endif

	                  endif		!  (status from processing spectra

c
c  Close the temporary files if no records were read the last time.
c
	               elseif ((read_status .eq. %loc(fcf_eof))  .and.
     .			       (num .eq. 0)) then
	                  if (fcc_write_vs .eq. fac_present) then
	                     close (unit=vs_lun, iostat=io_stat)
	                     if (io_stat .ne. 0) then
	                        status = %loc(fcf_tempclose)
	                        call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_vs(1:fcc_toutlenvs),
     .					   %val(io_stat))
	                     endif
	                  endif
	                  if ((fcc_write_ds .eq. fac_present)  .or.
     .			      (fcc_write_cs .eq. fac_present)) then
	                     close (unit=cs_lun, iostat=io_stat)
	                     if (io_stat .ne. 0) then
	                        status = %loc(fcf_tempclose)
	                        if (fcc_write_ds .eq. fac_present) then
	                           call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_ds(1:fcc_toutlends),
     .					   %val(io_stat))
	                        elseif (fcc_write_cs .eq. fac_present) then
	                           call lib$signal (fcf_tempclose, %val(2),
     .					   fcc_temp_outfile_cs(1:fcc_toutlencs),
     .					   %val(io_stat))
	                        endif
	                     endif
	                  endif

	               else
	                  status = %loc(fcf_chanerr)
	                  call lib$signal (fcf_chanerr)
	               endif	!  (read_status from fcf_read_coadd

	            enddo	!  (do while read_status

c
c  Write the voltage spectra to the Cobetrieve archives.
c
	            if ((status .eq. %loc(fcf_normal))  .and.
     .	                (fcc_write_vs .eq. fac_present)) then
	               if (fcc_sky .eq. fac_present) then
	                  status = fcf_write_sky_spec (1, vs_lun)
	               else
	                  status = fcf_write_cal_spec (1, vs_lun)
	               endif
	            endif

c
c  Write the calibrated or differential spectra to the Cobetrieve archives.
c
	            if (status .eq. %loc(fcf_normal)) then
	               if (fcc_write_ds .eq. fac_present) then
	                  if (fcc_sky .eq. fac_present) then
	                     status = fcf_write_sky_spec (2, cs_lun)
	                  else
	                     status = fcf_write_cal_spec (2, cs_lun)
	                  endif
	               elseif (fcc_write_cs .eq. fac_present) then
	                  if (fcc_sky .eq. fac_present) then
	                     status = fcf_write_sky_spec (3, cs_lun)
	                  else
	                     status = fcf_write_cal_spec (3, cs_lun)
	                  endif
	               endif
	            endif

c
c  Update the processing report with the number of spectra processed.
c
	            if ((status .eq. %loc(fcf_normal))  .and.
     .	                (fcc_report .eq. fac_present)) then
	               status = fcf_update_report ()
	            endif

c
c  Write out number of spectra processed.
c
	            if (status .eq. %loc(fcf_normal)) then
	               if (fcc_sky .eq. fac_present) then
	                  write (lun_out,30,iostat=io_stat) fcc_nspec, fcc_npix
	               else
	                  write (lun_out,40,iostat=io_stat) fcc_nspec
	               endif
	            endif
  30	format (/, 4x, 'Number of sky spectra processed:      ', I5,
     .		/, 4x, 'Number of pixels containing spectra:  ', I5, //)
  40	format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)


C
C  Close all remaining open files.
C

c
c  Verify that the temporary RMS files are closed.
c
	            if (vs_lun .ne. 0) then
	               close (vs_lun)
	               rstatus = fut_free_lun(vs_lun)
	               if (rstatus .ne. %loc(fut_normal)) then
	                  call lib$signal (%val(rstatus))
	               endif
	            endif

	            if (cs_lun .ne. 0) then
	               close (cs_lun)
	               rstatus = fut_free_lun(cs_lun)
	               if (rstatus .ne. %loc(fut_normal)) then
	                  call lib$signal (%val(rstatus))
	               endif
	            endif

	         endif	!	(status from fcf_open_coadd

	      endif	!  	(cstatus from cct_open_config

c
c  Close the configuration files.
c
	      if (ref_open_seq) then
	         cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fcf_cfgseqclose)
	            call lib$signal (fcf_cfgseqclose, %val(1), %val(cstatus))
	         endif
	      endif

	      if (ref_open_dir) then
	         cstatus = cct_close_config (ndset_dir, dir_lun, dir_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fcf_cfgdirclose)
	            call lib$signal (fcf_cfgdirclose, %val(1), %val(cstatus))
	         endif
	      endif

	   else
	      status = %loc(fcf_ctinit)
	      call lib$signal (fcf_ctinit, %val(1), %val(ct_stat(1)))
	   endif		!  (ct_stat from ct_init

	endif			!  (status from fcf_parse and fcf_init_report

c
c  Signal the program completion status.
c
	if ((status .eq. %loc(fcf_normal))  .or.
     .	    (status .eq. %loc(fcf_chanerr))) then
	   call lib$signal (fcf_normal)
	else
	   call lib$signal (fcf_failure)
	endif

c
c  Close the processing report file.
c
	if (fcc_report .eq. fac_present) then
	   close (unit=fut_report_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      call lib$signal (fcf_repclose, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif
	   rstatus = fut_free_lun (fut_report_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	endif


	end
