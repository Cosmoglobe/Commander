	integer * 4 function  fsl_open_sky_coadd (ct_lun)

c-------------------------------------------------------------------------------
c
c	Function FSL_OPEN_SKY_COADD
c
c	This function identifies and then opens the Cobetrieve coadd skymap
c	file specified by the command line invocation.  The function uses either
c	the specified file name extension or the jstart and jstop times to
c       identify the skymap. The function builds the filenames for output
c       spectra files from the input coadd file names.
c
c	Author:
c                FCF_Open_Sky_Coadd
c                Gene Eplee
c		 General Sciences Corp.
c		 17 June 1992
c
c                FSL_Open_Sky_Coadd
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
c		cct_query_catalog
c		cct_query_ttg_catalog
c		csa_open_skymap
c		ct_binary_to_gmt
c		lib$signal
c		str$trim
c		time_gt
c
c	Include files:
c		cct_query_catalog_record.txt
c		cct_query_ttg_catalog_record.txt
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
c       Modified FCF_Open_Sky_Coadd to FSL_Open_Sky_Coadd for the new FIRAS
c       pipeline which will process long spectra to get improved frequency
c       resolution.
c           1. Changed status, function and filenames.
c           2. Removed the build of temporary RMS file names.
c           3. Removed code related to pick up of scan mode for low frequency
c              long fast data. FIL already processes the 6 channel/scan modes.
c
c	Fred Shuman, Hughes STX Corporation, 1995 Sep 27
c           4. Changed  fac_scan_mode_ids  to  fac_scan_mode_idsL.
c-------------------------------------------------------------------------------

	implicit none

	include '(cct_query_catalog_record)'
	include '(cct_query_ttg_catalog_record)'
	include '(fut_params)'
	include '(fsl_invoc)'

	character * 64	filebuff1	!  file name buffer
	character * 64	filebuff2	!  file name buffer

	integer * 2	fblen1		!  file name buffer length
	integer * 2	fblen2		!  file name buffer length
	integer * 2	num_cats	!  number of CT skymap records

	integer * 4	activity(2)	!  binary skymap activity time
	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  Cobetrieve logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	cct_query_catalog
	integer * 4	cct_query_ttg_catalog
	integer * 4	csa_open_skymap
	external	csa_open_skymap
	logical * 1	time_gt

	dictionary 'ccm_cme_catalog_entry'
	record /ccm_cme_catalog_entry/ cats(50)
	record /query_catalog/ query_cat
	record /query_ttg_catalog/ query_ttg

	external	cct_normal
	external	cct_q_no_cat_entry
	external	fsl_csaifgopen
	external	fsl_ctquerycat
	external	fsl_ctnocatrecs
	external	fsl_normal

	status = %loc(fsl_normal)

	if (fcc_file .eq. fac_present) then
C
C  Define the skymap file names.
C

	   query_cat.archive_id = 'CSDR$FIRAS_IN'
	   query_cat.filename   = 'FIL_SKY_' // fac_channel_ids(fcc_chan)
     .					     // fac_scan_mode_idsL(fcc_smode)
     .					     // '.' // fcc_file_ext

c
c  Query the Cobetrieve catalog.
c
	   cstatus = cct_query_catalog (query_cat, cats(1))
	   call str$trim (filebuff1, query_cat.archive_id, fblen1)
	   call str$trim (filebuff2, query_cat.filename, fblen2)

	   if (cstatus .eq. %loc(cct_q_no_cat_entry)) then
	      fsl_open_sky_coadd = %loc(fsl_ctnocatrecs)
	      call lib$signal (fsl_ctnocatrecs, %val(2), filebuff1(1:fblen1),
     .							 filebuff2(1:fblen2))
	      return
	   elseif (cstatus .ne. %loc(cct_normal)) then
	      fsl_open_sky_coadd = %loc(fsl_ctquerycat)
	      call lib$signal (fsl_ctquerycat, %val(3), filebuff1(1:fblen1),
     .					 filebuff2(1:fblen2), %val(cstatus))
	      return
	   endif

c
c  Get the start and stop time for the skymap
c
	   fcc_jstart(1) = cats(1).initial_time(1)
	   fcc_jstart(2) = cats(1).initial_time(2)
	   fcc_jstop(1) = cats(1).final_time(1)
	   fcc_jstop(2) = cats(1).final_time(2)
	   call ct_binary_to_gmt (fcc_jstart,fcc_jstart_time)
	   call ct_binary_to_gmt (fcc_jstop,fcc_jstop_time)
	   fcc_time_range  = fcc_jstart_time // ';' //
     .			     fcc_jstop_time // ';'


	else
C
C  Identify the coadd skymap file from the timetags.
C
	   query_ttg.archive_id    = 'CSDR$FIRAS_IN'
	   query_ttg.dataset_name  = 'FIL_SKY_' // fac_channel_ids(fcc_chan)
     .						// fac_scan_mode_idsL(fcc_smode)
	   query_ttg.start_time(1) = fcc_jstart(1)
	   query_ttg.start_time(2) = fcc_jstart(2)
	   query_ttg.stop_time(1)  = fcc_jstop(1)
	   query_ttg.stop_time(2)  = fcc_jstop(2)

c
c  Query the Cobetrieve catalog.
c
	   cstatus = cct_query_ttg_catalog (query_ttg, cats, num_cats)
	   call str$trim (filebuff1, query_ttg.archive_id, fblen1)
	   call str$trim (filebuff2, query_ttg.dataset_name, fblen2)

	   if (cstatus .ne. %loc(cct_normal)) then
	      fsl_open_sky_coadd = %loc(fsl_ctquerycat)
	      call lib$signal (fsl_ctquerycat, %val(3), filebuff1(1:fblen1),
     .					 filebuff2(1:fblen2), %val(cstatus))
	      return
	   elseif (num_cats .eq. 0) then
	      fsl_open_sky_coadd = %loc(fsl_ctnocatrecs)
	      call lib$signal (fsl_ctnocatrecs, %val(2), filebuff1(1:fblen1),
     .							 filebuff2(1:fblen2))
	      return
	   endif

c
c  Determine which skymap has the last activity time.
c
	   activity(1) = 0
	   activity(2) = 0

	   do j = 1, num_cats
	      if (time_gt(cats(j).activity_time, activity)) then
	         activity(1) = cats(j).activity_time(1)
	         activity(2) = cats(j).activity_time(2)
	         fcc_file_ext = cats(j).filename_extension
	      endif
	   enddo

	endif	!	(fcc_file


	if (status .eq. %loc(fsl_normal)) then
C
C  Determine the file names.
C

c
c  Determine the input file name.
c
	   fcc_infile = 'CSDR$FIRAS_IN:FIL_SKY_'// fac_channel_ids(fcc_chan)
     .						// fac_scan_mode_idsL(fcc_smode)
     .						// '.' // fcc_file_ext
	   call str$trim (fcc_infile, fcc_infile, fcc_inlen)

c
c  Determine the output spectra file names from input coadd file name. Only
c  command line selected spectra files will actually be needed for the FSL run.
c

	   fcc_outfile_vs = 'CSDR$FIRAS_OUT:FSL_VSK_' // fcc_scan_mode // '.' //
     .			     fcc_file_ext
	   fcc_outfile_ds = 'CSDR$FIRAS_OUT:FSL_DSK_' // fcc_scan_mode // '.' //
     .			     fcc_file_ext
	   fcc_outfile_cs = 'CSDR$FIRAS_OUT:FSL_SKY_' // fcc_scan_mode // '.'//
     .			     fcc_file_ext
	   call str$trim (fcc_outfile_vs, fcc_outfile_vs, fcc_outlenvs)
	   call str$trim (fcc_outfile_ds, fcc_outfile_ds, fcc_outlends)
	   call str$trim (fcc_outfile_cs, fcc_outfile_cs, fcc_outlencs)

C
C  Open the coadd skymap.
C
	   open (unit=ct_lun, file=fcc_infile, iostat=io_stat, status='old',
     .	         form='unformatted', recordtype='fixed', readonly,
     .		 useropen=csa_open_skymap)

	   if (io_stat .ne. 0) then
	      status = %loc(fsl_csaifgopen)
	      call lib$signal (fsl_csaifgopen, %val(2), fcc_infile(1:fcc_inlen),
     .					       %val(io_stat))
	   endif


	endif		!  (status from catalog query


	fsl_open_sky_coadd = status

	return
	end
