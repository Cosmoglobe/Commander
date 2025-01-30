	integer * 4 function  fsl_open_cal_coadd (ct_lun)

c-------------------------------------------------------------------------------
c
c	Function FSL_OPEN_CAL_COADD
c
c	This function identifies and then opens the Cobetrieve calibration coadd
c	file(s) specified by the command line invocation.  The function uses
c	either the specified file name extension or the jstart and jstop times
c	to identify the files. It also builds the filenames for output spectra
c       files from the input coadd file names.
c
c	Author:
c                 FCF_Open_Cal_Coadd
c                 Gene Eplee
c		  General Sciences Corp.
c		  17 June 1992
c
c                FSL_Open_Cal_Coadd
c                Shirley M. Read
c                Hughes STX Corporation
c                July 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ct_lun		integer * 4		Cobetrieve logical unit
c							number
c
c	Output:
c		none
c
c	Subroutines called:
c		cct_query_tod_catalog
c		cct_query_catalog
c		ct_binary_to_gmt
c		ct_connect_read
c		lib$signal
c		str$trim
c
c	Include files:
c		cct_query_catalog_record.txt
c		cct_query_tod_catalog_record.txt
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
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
c       Shirley M. Read, Hughes STX Corporation, July 26, 1995
c       Modified FCF_Open_Cal_Coadd to FSL_Open_Cal_Coadd for the new FIRAS
c       pipeline which will process long spectra to get improved frequency
c       resolution.
c           1. Changed status, function and filenames.
c           2. Removed the build of temporary RMS file names.
c
c	Fred Shuman, Hughes STX Corporation, 1995 Sep 27
c           3. Changed  fac_scan_mode_ids  to  fac_scan_mode_idsL.
c-------------------------------------------------------------------------------

	implicit none

	include '(cct_query_catalog_record)'
	include '(cct_query_tod_catalog_record)'
	include '(fut_params)'
	include '(fsl_invoc)'

	character * 64	filebuff1	!  file name buffer
	character * 64	filebuff2	!  file name buffer

	integer * 2	fblen1		!  file name buffer length
	integer * 2	fblen2		!  file name buffer length
	integer * 2	num_cats	!  number of coadd file found

	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	smode		!  mtm scan mode
	integer * 4	status		!  return status

	integer * 4	cct_query_catalog
	integer * 4	cct_query_tod_catalog
	integer * 4	ct_connect_read
	external	ct_connect_read

	dictionary 'ccm_cme_catalog_entry'
	record /ccm_cme_catalog_entry/ cats(50)
	record /query_catalog/ query_cat
	record /query_tod_catalog/ query_tod

	external	cct_normal
	external	cct_q_no_cat_entry
	external	fsl_ctifgopen
	external	fsl_ctquerycat
	external	fsl_ctnocatrecs
	external	fsl_normal

	status = %loc(fsl_normal)

C
C  Set the scan mode from the FSL_Invoc include file.
C
	smode  = fcc_smode

	if (fcc_file .eq. fac_present) then
C
C  Define the input calibration file.
C
	   fcc_infile = 'CSDR$FIRAS_IN:FIL_CAL_' // fac_channel_ids(fcc_chan)
     .						 // fac_scan_mode_idsL(smode)
     .						 // '.' // fcc_file_ext
	   call str$trim (fcc_infile, fcc_infile, fcc_inlen)

	   query_cat.archive_id = 'CSDR$FIRAS_IN'
	   query_cat.filename =   'FIL_CAL_' // fac_channel_ids(fcc_chan)
     .					     // fac_scan_mode_idsL(smode)
     .					     // '.' // fcc_file_ext

c
c  Query the Cobetrieve catalog.
c
	   cstatus = cct_query_catalog (query_cat, cats(1))
	   call str$trim (filebuff1, query_cat.archive_id, fblen1)
	   call str$trim (filebuff2, query_cat.filename, fblen2)

	   if (cstatus .eq. %loc(cct_q_no_cat_entry)) then
	      fsl_open_cal_coadd = %loc(fsl_ctnocatrecs)
	      call lib$signal (fsl_ctnocatrecs, %val(2), filebuff1(1:fblen1),
     .							 filebuff2(1:fblen2))
	      return
	   elseif (cstatus .ne. %loc(cct_normal)) then
	      fsl_open_cal_coadd = %loc(fsl_ctquerycat)
	      call lib$signal (fsl_ctquerycat, %val(3), filebuff1(1:fblen1),
     .					 filebuff2(1:fblen2), %val(cstatus))
	      return
	   endif

c
c  Get the start and stop time for the coadd file.
c
	   fcc_jstart(1) = cats(1).initial_time(1)
	   fcc_jstart(2) = cats(1).initial_time(2)
	   fcc_jstop(1)  = cats(1).final_time(1)
	   fcc_jstop(2)  = cats(1).final_time(2)
	   call ct_binary_to_gmt (fcc_jstart,fcc_jstart_time)
	   call ct_binary_to_gmt (fcc_jstop,fcc_jstop_time)
	   fcc_time_range  = fcc_jstart_time // ';' //
     .			     fcc_jstop_time // ';'


	else
C
C  Identify the calibration coadd files from the timetags.
C
	   fcc_infile = 'CSDR$FIRAS_IN:FIL_CAL_' // fac_channel_ids(fcc_chan)
     .						 // fac_scan_mode_idsL(smode)
     .						 // '/' // fcc_time_range
	   call str$trim (fcc_infile, fcc_infile, fcc_inlen)

	   query_tod.archive_id    = 'CSDR$FIRAS_IN'
	   query_tod.dataset_name  = 'FIL_CAL_' // fac_channel_ids(fcc_chan)
     .						// fac_scan_mode_idsL(smode)
	   query_tod.start_time(1) = fcc_jstart(1)
	   query_tod.start_time(2) = fcc_jstart(2)
	   query_tod.stop_time(1)  = fcc_jstop(1)
	   query_tod.stop_time(2)  = fcc_jstop(2)

c
c  Query the Cobetrive catalog.
c
	   cstatus = cct_query_tod_catalog (query_tod, cats, num_cats)
	   call str$trim (filebuff1, query_tod.archive_id, fblen1)
	   call str$trim (filebuff2, query_tod.dataset_name, fblen2)

	   if (cstatus .eq. %loc(cct_q_no_cat_entry)) then
	      fsl_open_cal_coadd = %loc(fsl_ctnocatrecs)
	      call lib$signal (fsl_ctnocatrecs, %val(2), filebuff1(1:fblen1),
     .							 filebuff2(1:fblen2))
	      return
	   elseif (cstatus .ne. %loc(cct_normal)) then
	      fsl_open_cal_coadd = %loc(fsl_ctquerycat)
	      call lib$signal (fsl_ctquerycat, %val(3), filebuff1(1:fblen1),
     .					 filebuff2(1:fblen2), %val(cstatus))
	      return
	   endif

c
c  Get the file name extension.
c
	   fcc_file_ext = cats(1).filename_extension(1:3) //
     .			  fcc_jstart_time(1:7) // '_' // fcc_jstop_time(1:7)

	endif		!  (fcc_file

	if (status .eq. %loc(fsl_normal)) then
C
C  Determine the output spectra file names from input coadd file name. Only
C  command line selected spectra files will actually be needed for the FSL run.
C

	   fcc_outfile_vs = 'CSDR$FIRAS_OUT:FSL_VCL_' // fcc_scan_mode // '.' //
     .			     fcc_file_ext
	   fcc_outfile_ds = 'CSDR$FIRAS_OUT:FSL_DCL_' // fcc_scan_mode // '.' //
     .			     fcc_file_ext
	   fcc_outfile_cs = 'CSDR$FIRAS_OUT:FSL_CAL_' // fcc_scan_mode // '.' //
     .			     fcc_file_ext
	   call str$trim (fcc_outfile_vs, fcc_outfile_vs, fcc_outlenvs)
	   call str$trim (fcc_outfile_ds, fcc_outfile_ds, fcc_outlends)
	   call str$trim (fcc_outfile_cs, fcc_outfile_cs, fcc_outlencs)

C
C  Open the coadd files.
C
	   open (unit=ct_lun, file=fcc_infile, iostat=io_stat, status='old',
     .		 useropen=ct_connect_read)

	   if (io_stat .ne. 0) then
	      status = %loc(fsl_ctifgopen)
	      call lib$signal (fsl_ctifgopen, %val(2), fcc_infile(1:fcc_inlen),
     .					      %val(io_stat))
	   endif

	endif	!	(status from catalog query


	fsl_open_cal_coadd = status

	return
	end
