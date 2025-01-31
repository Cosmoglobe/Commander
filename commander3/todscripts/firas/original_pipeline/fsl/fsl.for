	Program FSL

c-------------------------------------------------------------------------------
c
c	FSL_FIRAS_SPECTRA_LONG
c
c	This program drives the Firas calibration routines.  The program reads
c	coadded IFGs from and writes calibrated spectra to the Cobetrieve
c	archives.  The program reads in the Firas calibration model solutions
c	from a Cobetrieve reference archive.
c
c	For a given channel and scan mode, FSL reads in coadded IFGs and
c	transforms them into voltage spectra.  It then applies the calibration
c	model to both the spectra and the spectrum variances.  FSL calibrates
c	both sky and calibration data.
c
c	Author:   Shirley M. Read
c                 Hughes STX Corporation
c                 July 6, 1995
c                 SPR 12288
c
c                 FSL is a rewrite of the FIRAS Pipeline FCF Calibrate FIRAS
c                 to calibrate spectra of varying length. The FIRAS Pipeline
c                 was changed to process longer spectral arrays in order to
c                 improve the frequency resolution. The previous FCF was
c	          rewritten from the Pass 2 FFC for FIRAS Pass 3 processing by:
c                     Gene Eplee
c		      General Sciences Corporation
c		      23 September 1992
c                     SER 8292
c
c                 The original FFC was written for previous data processsing
c                 passes by:
c                     Robert Kummerer
c                     Hughes STX Corporation
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		ct_init
c		cut_display_banner
c		cut_register_version
c		fsl_calibrate_spectra
c		fsl_close_config
c               fsl_close_spectra
c		fsl_display_spectra
c		fsl_initialize_report
c		fsl_open_coadd
c		fsl_open_config
c               fsl_open_spectra
c		fsl_parse
c		fsl_produce_spectra
c		fsl_read_coadd
c		fsl_update_report
c		fsl_write_cal_spec
c		fsl_write_sky_spec
c		fut_free_lun
c		lib$establish
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c               fut_error
c               fut_params
c               fsl_config
c               fsl_invoc
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
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
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, July 12, 1995
c       Modified FCF to FSL for the new FIRAS pipeline which will process
c       long spectra to get improved frequency resolution.
c           1. Changed program, function, and status names for FSL.
c           2. Removed code for producing temporary files.
c           3. Redesigned the program to perform fucntions in a more logical
c              order. Restructured the if-else-endif blocking to decrease the
c              levels of nesting.
c           4. Added functions to open and close the spectra files. Transferred
c	       the open and close of configuration files to FSL functions.
c           5. Removed call to iwkin to allocate workspace and corresponding
c              workspace common block. By default, the current system
c              automatically allocates the correct amount of space.
c	    6. Removed obsolete External Symbol FSL_TempClose and declaration
c              of obsolete module FSL_Write_Temp_Spec.
c           7. Removed call to lib$establish for FUT_Error as the condition
c              handler from the "if block" after FSL_Initialize_Report and
c              moved the call before the FSL_Parse. It is now unconditional.
c
c	Fred Shuman, Hughes STX Corporation, 1995 Sep 27
c           8. Changed  fac_scan_mode_ids  to  fac_scan_mode_idsL.
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_error)'
	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_invoc)'

	character * 79	cmd_line(11)		!  command line invocation
	character * 14	current_gmt		!  GMT time of invocation

	complex	* 16	cspec(361)		!  calibrated spectrum
	complex	* 16	vspec(361)		!  voltage spectrum

	integer *  2	ct_stat(20)		!  Cobetrieve return status

	integer *  4	cmdlen(11)		!  length of command lines in
						!    invocation
	integer *  4	cs_lun			!  calibrated spectra file lun
	integer *  4	cstatus			!  Get_Config return status
	integer *  4	ct_lun			!  coadd file lun
	integer *  4	io_stat			!  I/O return status
	integer *  4	j			!  a counter
	integer *  4	lun_out /6/		!  terminal lun
	integer *  4	ncmd			!  number of command lines in
						!    invocation
	integer *  4    flag                    !  spectrum type flag
						!  1 = voltage
						!  2 = differential
						!  3 = calibrated
	integer *  4	num			!  number of current coadds read
	integer *  4	parse_status		!  command line parse return
						!    status
	integer *  4	read_status		!  coadd read return status
	integer *  4	rstatus			!  return status
	integer *  4	status			!  return status
	integer *  4	tstatus			!  return status
	integer *  4	vs_lun			!  voltage spectra file lun

	logical *  1	ref_open_dir /.false./	!  reference dataset open flag
	logical *  1	ref_open_seq /.false./	!  reference dataset open flag

	real	*  8	cvar(3,361)		!  calibrated spectrum variances
	real	*  8	vvar(3,361)		!  voltage spectrum variances

	character * 6   fsl_version
	parameter       (fsl_version = '13.9')	!  FSL version number for banner

	integer *  4	cut_display_banner
	integer *  4	cut_register_version
	integer *  4	fsl_calibrate_spectra
	integer *  4    fsl_close_config
	integer *  4	fsl_close_spectra
	integer *  4	fsl_display_spectra
	integer *  4	fsl_initialize_report
	integer *  4	fsl_open_coadd
	integer *  4    fsl_open_config
	integer *  4    fsl_open_spectra
	integer *  4	fsl_parse
	integer *  4	fsl_produce_spectra
	integer *  4	fsl_read_coadd
	integer *  4	fsl_update_report
	integer *  4	fsl_write_cal_spec
	integer *  4	fsl_write_sky_spec
	integer *  4	fut_free_lun

	dictionary 'fil_sky'
	record /fil_sky/ coadd_recs(fcc_max_coadd)
	dictionary 'fsl_sky'
	record /fsl_sky/ cspec_recs(fcc_max_coadd)
	record /fsl_sky/ vspec_recs(fcc_max_coadd)

	external	fut_error

	external	cct_normal
	external	fsl_chanerr
	external	fsl_ctinit
	external	fsl_eof
	external	fsl_failure
	external	fsl_normal
	external	fsl_repclose
	external	fut_normal

C
C  Initialize the program.
C
	call lib$establish (fut_error)
c
c  Parse the command line.
c
	parse_status = fsl_parse (current_gmt, ncmd, cmd_line, cmdlen)

c
c  Print the banner.
c
	rstatus = cut_register_version (fsl_version)
	rstatus = cut_display_banner (lun_out, 80,
     .		'FIRAS Facility FSL_Spectra_Long')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fsl_initialize_report (current_gmt, ncmd, cmd_line, cmdlen,
     .					   parse_status, fsl_version)
	else
	   status = %loc(fsl_normal)
	endif

	if ((parse_status .eq. %loc(fsl_normal))  .and.
     .	    (status .eq. %loc(fsl_normal))) then
c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)

	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Open the configuration files for sequential and direct access.
c
	      cstatus = fsl_open_config (ref_open_seq, ref_open_dir)
	      if (cstatus .ne. %loc(fsl_normal)) status = cstatus
	   else
	      status = %loc(fsl_ctinit)
	      call lib$signal (fsl_ctinit, %val(1), %val(ct_stat(1)))
	   endif		!  (ct_stat from ct_init

	endif			!  (status from fsl_parse and fsl_init_report)

C
C  Process the spectra for the specified channel and scan mode.
C
	if (status .eq. %loc(fsl_normal)) then	! (status from initialization)
c
c  Open the coadd file.
c
	     status = fsl_open_coadd (ct_lun)

	     if (status .eq. %loc(fsl_normal)) then
c
c  Initialize the coadd read.
c
	            fcc_nspec      = 0
	            num            = 0
	            fcc_more_plots = .true.
	            read_status    = %loc(fsl_normal)

	            write(lun_out,20) fac_channel_ids(fcc_chan),
     .			     fac_scan_mode_idsL(fcc_smode)
 20	            format (/, 4x, 'Processing spectra for channel ', a,
     .			           ', scan mode ', a, //)
c
c   Open the spectra file or files.
c
		    if ((fcc_write_vs .eq. fac_present) .or.
     .                  (fcc_write_cs .eq. fac_present) .or.
     .                  (fcc_write_ds .eq. fac_present)) then
	               status = fsl_open_spectra (vs_lun, cs_lun)
		    endif

	     endif     ! (status .eq. %loc(fsl_normal))  from fsl_open_coadd

	endif     ! (status .eq. %loc(fsl_normal))  from initialization

c
c   Read in a set of data for a particular scan mode. Calibration data are
c   read sequentially, one per call to Ct_Read_Arcv. Sky data are read by
c   pixel number through calls to CSA_Read_Pixels. Each call for a pixel
c   may return several records for the pixel. Thus an inner loop is required
c   to produce, calibrate and display spectra for each coadd.
c
	if (status .eq. %loc(fsl_normal)) then	  ! (status from open)

	   do while (read_status .eq. %loc(fsl_normal))

	       read_status = fsl_read_coadd (ct_lun, num, coadd_recs)

	       if (((read_status .eq. %loc(fsl_normal)) .or.
     .              (read_status .eq. %loc(fsl_eof))) .and.
     .		    (num .gt. 0)) then               ! data were read

	           status = %loc(fsl_normal)
	           j = 0
	           fcc_nspec = fcc_nspec + num

	           do while ((status .eq. %loc(fsl_normal)) .and. (j .lt. num))
	              j = j + 1
c
c  Produce the voltage spectra.
c
	              status = fsl_produce_spectra (coadd_recs(j), vspec, vvar,
     .						  vspec_recs(j), cspec_recs(j))
c
c  Calibrate the spectra.
c
	              if (fcc_calibrate .eq. fac_present) then
	                 status = fsl_calibrate_spectra (vspec, vvar, cspec,
     .                                                   cvar, cspec_recs(j))
	              endif
c
c  Plot the spectra.
c
	              fcc_jump = fcc_jump - 1

	              if ((status .eq. %loc(fsl_normal))  .and.
     .			  (fcc_plot .eq. fac_present)  .and.
     .			  (fcc_more_plots .eq. .true.)  .and.
     .			  (fcc_jump .le. 0)) then

	                 tstatus = fsl_display_spectra ()

	                 if (tstatus .eq. fac_not_present) then
	                     fcc_more_plots = .false.
	                 endif
	              endif
	           end do		!  (while j .lt. num
c
c  Write the voltage spectra to the Cobetrieve archives.
c
	           if ((status .eq. %loc(fsl_normal)) .and. (num .gt. 0)
     .	              .and. (fcc_write_vs .eq. fac_present)) then
	              flag = 1
	              if (fcc_sky .eq. fac_present) then

	                 status = fsl_write_sky_spec (flag, vs_lun, num,
     .                            vspec_recs)
	              else
	                 status = fsl_write_cal_spec (flag, vs_lun, num,
     .                            vspec_recs)
	              endif
	           endif		! write vspec_recs
c
c  Write the calibrated or differential spectra to the Cobetrieve archives.
c
	           if ((status .eq. %loc(fsl_normal)) .and. (num .gt. 0)) then
	              if (fcc_write_ds .eq. fac_present) then
	                 flag = 2
	                 if (fcc_sky .eq. fac_present) then

	                    status = fsl_write_sky_spec (flag, cs_lun, num,
     .                                                   cspec_recs)
	                 else
	                    status = fsl_write_cal_spec (flag, cs_lun, num,
     .                                                   cspec_recs)
	                 endif
	              elseif (fcc_write_cs .eq. fac_present) then
	                 flag = 3
	                 if (fcc_sky .eq. fac_present) then
	                    status = fsl_write_sky_spec (flag, cs_lun, num,
     .                                                   cspec_recs)
	                 else
	                    status = fsl_write_cal_spec (flag, cs_lun, num,
     .                                                   cspec_recs)
	                 endif
	              endif
	           endif		! write the spectra to archive

	       elseif ((read_status .ne. %loc(fsl_normal)) .and. ! error status
     .               (read_status .ne. %loc(fsl_eof))) then    ! fsl_read_coadd
	            status = %loc(fsl_chanerr)
	            call lib$signal (fsl_chanerr)
	       endif  !  read_status from fsl_read_coadd

C
C  If no more data, close all spectra archive files.
C
	       if (read_status .eq. %loc(fsl_eof)) then

	          status = fsl_close_spectra (vs_lun, cs_lun)

	       endif
	   enddo	!  (do while read_status .eq. %loc(fsl_normal))
	endif		!  (status .eq. %loc(fsl_normal) after open)

c
c  Update the processing report with the number of spectra processed.
c
	if ((status .eq. %loc(fsl_normal))  .and.
     .	    (fcc_report .eq. fac_present)) then
	       status = fsl_update_report ()
	endif
c
c  Write out number of spectra processed.
c
	 if (status .eq. %loc(fsl_normal)) then
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
C  Close the configuration files.
C
	if ((ref_open_seq) .or. (ref_open_dir)) then

	    cstatus = fsl_close_config (ref_open_seq, ref_open_dir)
	    if (cstatus .ne. %loc(fsl_normal)) status = cstatus
	endif
c
c  Signal the program completion status.
c
	if (status .eq. %loc(fsl_normal)) then
	   call lib$signal (fsl_normal)
	elseif (status .eq. %loc(fsl_chanerr)) then
	   call lib$signal (fsl_chanerr)
	else
	   call lib$signal (fsl_failure)
	endif
c
c  Close the processing report file.
c
	if (fcc_report .eq. fac_present) then
	   close (unit=fut_report_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      call lib$signal (fsl_repclose, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif
	   rstatus = fut_free_lun (fut_report_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	endif

	end
