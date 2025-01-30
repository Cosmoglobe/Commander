	integer*4  function  fip_open_cal (flun, alun)

c------------------------------------------------------------------------------
c
c	Function FIP_OPEN_CAL
c
c	This function opens the FIRAS input fcf_cal data skymap and the ADB
c	output fip data file.
c
c	Author:  Larry P. Rosen, Hughes STX, 5 October 1994
c
c------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		alun		integer * 4		ADB output file lun
c		flun		integer * 4		FIRAS input file lun
c
c	Subroutines called:
c		cct_query_catalog
c		ct_connect_read
c		lib$signal
c
c	Include files:
c		fip_invoc_sky.txt
c		fut_params.txt
c		cct_query_catalog_record.txt
c
c------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_sky)'
	include '(cct_query_catalog_record)'

	integer*4	alun			!  ADB skymap lun
	integer*4	flun			!  FIRAS skymap lun

	integer*4	status			!  return status
	integer*4	rstatus			!  return status
	integer*2	colon			!  position of colon
	integer*4	cstatus			!  CT or CSA return status
	dictionary 'ccm_cme_catalog_entry'
	record /ccm_cme_catalog_entry/ cats(50)
	record /query_catalog/ query_cat
	integer*4	io_stat			!  I/O return status

C Functions

	integer * 4	fut_get_lun
	integer * 4	cct_query_catalog
	integer * 4	ct_connect_read
	external	ct_connect_read

C Externals

	external	fip_normal
	external	fip_lunerr
	external	fip_openerr
	external	fip_ctnocatrecs
	external	fip_ctquerycat
	external	fut_normal
	external	cct_q_no_cat_entry
	external	cct_normal

C  Begin
C  Initialize the open.
C
	status = %loc(fip_normal)
C
C  Get the logical unit numbers for the files.
C
	rstatus = fut_get_lun (flun)
	rstatus = fut_get_lun (alun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal ( fip_lunerr, %val(1), %val (rstatus) )
	   status = %loc (fip_lunerr)
	else
	   query_cat.archive_id = 'CSDR$FIRAS_IN'
	   colon = index (fcc_infile, ':')
	   query_cat.filename = fcc_infile (colon+1:fcc_inlen)
C
C  Query the Cobetrieve catalog.
C
	   cstatus = cct_query_catalog (query_cat, cats(1))
	   if (cstatus .eq. %loc(cct_q_no_cat_entry)) then
	      status = %loc(fip_ctnocatrecs)
	      call lib$signal (fip_ctnocatrecs, %val(2), query_cat.archive_id,
     *	                       query_cat.filename)
	   elseif (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(fip_ctquerycat)
	      call lib$signal (fip_ctquerycat, %val(3), query_cat.archive_id,
     *	                       query_cat.filename, %val(cstatus))
	   else
C
C  Open the files.
C
	      open (unit=flun, file=fcc_infile, iostat=io_stat, status='old',
     *	            useropen=ct_connect_read)

	      if (io_stat .ne. 0) then
 	         status = %loc(fip_openerr)
	         call lib$signal (fip_openerr, %val(2),
     *	                          fcc_infile(1:fcc_inlen), %val(io_stat))
	      else
C
C  Open the output ADB skymap file.
C
	         open ( unit=alun, file=fcc_outfile(1:fcc_outlen),
     *	                iostat=io_stat, status='new', form='unformatted',
     *	                recordtype='fixed', recl=fcc_alen,
     *	                access='sequential')

	         if (io_stat .ne. 0) then
 	            status = %loc(fip_openerr)
	            call lib$signal (fip_openerr, %val(2),
     *	                             fcc_outfile(1:fcc_outlen), %val(io_stat))
	         endif
	      endif
	   endif
	endif
	fip_open_cal = status
	return
	end
