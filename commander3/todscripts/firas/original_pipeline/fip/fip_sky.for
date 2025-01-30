	program fip_sky

c-------------------------------------------------------------------------------
c
c	Program FIP_SKY
c
c	This program reformats a FIRAS Production Skymap into the Initial
c	Release ADB format.  FIP_SKY reformats FCF, FAD, or FCS skymaps into
c	FIP_SKY skymaps.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  25 May 1993
c		  SER 11058
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		csa_read_pixels
c		csa_write_pixels
c		ct_init
c		cut_display_banner
c		cut_register_version
c		fip_close_skymaps
c		fip_init_sky_report
c		fip_frequency_cut
c		fip_galactic_cut
c		fip_open_skymaps
c		fip_parse_sky
c		fip_reformat_spectrum
c		fip_update_sky_report
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		csa_pixel_input_rec.txt
c		csa_pixel_output_rec.txt
c		ct$library:ctuser.inc
c		fip_invoc_sky.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Change input record buffer size to work with FCF records.
c	Gene Eplee, GSC, 7 July 1994
c
c-------------------------------------------------------------------------------

	implicit none

	include '(csa_pixel_input_rec)'
	include '(csa_pixel_output_rec)'
	include 'ct$library:ctuser.inc'
	include '(fut_error)'
	include '(fut_params)'
	include '(fip_invoc_sky)'

	character * 14	current_gmt	!  GMT time of invocation

	integer * 2	ct_stat(20)	!  Cobetrieve return status

	integer * 4	alun		!  ADB skymap lun
	integer * 4	blocks		!  block count for CSA
	integer * 4	cstatus		!  CSA return status
	integer * 4	flun		!  FIRAS skymap lun
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	lun_out /6/	!  terminal lun
	integer * 4	max_in		!  maximun number of pixels to read
					!  in one call to CSA
	integer * 4	npix		!  number of pixels processed
	integer * 4	nspec		!  number of spectra reformatd
	integer * 4	num_in		!  number of pixels to be read
					!    in one call to CSA
	integer * 4	num_out		!  number of records to be written
					!    in one call to CSA
	integer * 4	num_read	!  number of pixels read by CSA
	integer * 4	pixel_no	!  pixel number to be read by CSA
	integer * 4	parse_status	!  return status
	integer * 4	pstatus		!  pixel status
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	logical * 1	first_time	!  first time flag for frequency cut

	real	* 4	fnyq_hz		!  Nyquist frequency in hz
	real	* 4	fnyq_icm	!  Nyquist frequency in icm

	integer * 4	csa_read_pixels
	integer * 4	csa_write_pixels
	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	fip_close_skymaps
	integer * 4	fip_frequency_cut
	integer * 4	fip_galactic_cut
	integer * 4	fip_init_sky_report
	integer * 4	fip_open_skymaps
	integer * 4	fip_parse_sky
	integer * 4	fip_reformat_spectrum
	integer * 4	fip_update_sky_report
	integer * 4	fut_free_lun

	dictionary 'fcf_sky'
	record /fcf_sky/ in_recs(fac_max_skymap_recs)
	dictionary 'fip_sky'
	record /fip_sky/ out_recs(fac_max_skymap_recs)

	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist

	external	fut_error

	external	cct_normal
	external	csa_normal
	external	fip_csaread
	external	fip_csawrite
	external	fip_ctinit
	external	fip_failure
	external	fip_maxrec
	external	fip_normal
	external	fip_repclose
	external	fut_normal

c
c  Parse the command line.
c
	parse_status = fip_parse_sky (current_gmt)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .				     'FIRAS Facility FIP_SKY')
	write(lun_out,10)
 10	format (/)

c
c  Initialize the processing report.
c
	if (fcc_report .eq. fac_present) then
	   status = fip_init_sky_report (current_gmt, parse_status)
	   if (status .eq. %loc(fip_normal)) call lib$establish (fut_error)
	else
	   status = %loc(fip_normal)
	endif


	if ((status .eq. %loc(fip_normal))  .and.
     .	    (parse_status .eq. %loc(fip_normal))) then
C
C  Open the skymap files.
C

c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)

	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Open the files
c
	      status = fip_open_skymaps (flun, alun)


	      if (status .eq. %loc(fip_normal)) then
C
C  Reformat the skymap.
C

c
c  Loop over the pixels.
c
	         cstatus = %loc(csa_normal)
	         pstatus = %loc(fip_normal)
	         blocks    = 0
	         first_time = .true.
	         npix  = 0
	         nspec = 0
	         inlist.level_no = fac_skymap_level
	         max_in  = fac_max_skymap_recs
	         num_in  = 1
	         num_out = 1
	         pixel_no = -1

	         do while ((cstatus .eq. %loc(csa_normal))  .and.
     .			   (status .eq. %loc(fip_normal))  .and.
     .			   (pstatus .eq. %loc(fip_normal))  .and.
     .			   (pixel_no .lt. 6143))

c
c  Read the pixels from the FIRAS skymap.
c
	            pixel_no = pixel_no + 1
	            if (fcc_galexc .eq. fac_present) then
	               pstatus = fip_galactic_cut (pixel_no)
	            endif
	            if (pstatus .eq. %loc(fip_normal)) then
	               inlist.pixel_no = pixel_no
	               cstatus = csa_read_pixels (flun, inlist, num_in, in_recs,
     .					      max_in, outlist, num_read, blocks)

	               if (outlist.no_records .gt. fac_max_skymap_recs) then
	                  status = %loc(fip_maxrec)
	                  call lib$signal (fip_maxrec, %val(2),
     .				       %val(outlist.no_records), %val(pixel_no))
	               elseif (outlist.no_records .ne. 0) then
	                  if (cstatus .eq. %loc(csa_normal)) then
	                     npix = npix + 1
	                     nspec = nspec + outlist.no_records
c
c  Set the frequency limits the first time through.
c
	                     if (first_time) then
	                        first_time = .false.
	                        fcc_destriped = in_recs(1).spec_data.destriped
	                        fnyq_icm = in_recs(1).coad_spec_data.nyquist_icm
	                        fnyq_hz =in_recs(1).coad_spec_data.nyquist_hertz
	                        status = fip_frequency_cut (fnyq_icm, fnyq_hz)
	                     endif

c
c  Reformat the spectrum records.
c
	                     do j = 1, outlist.no_records
	                        status = fip_reformat_spectrum (in_recs(j),
     .								out_recs(j))
	                     enddo

	                     if (status .eq. %loc(fip_normal)) then
c
c  Write out the reformated records to the ADB skymap.
c
	                        j = 0
	                        do while ((cstatus .eq. %loc(csa_normal))  .and.
     .					  (j .lt. outlist.no_records))
	                           j = j + 1
	                           cstatus = csa_write_pixels (alun,
     .							       out_recs(j),
     .							       num_out, blocks)
	                        enddo

	                        if (cstatus .ne. %loc(csa_normal)) then
	                           status =%loc(fip_csawrite)
	                           call lib$signal (fip_csawrite, %val(2),
     .						    fcc_outfile(1:fcc_outlen),
     .						    %val(cstatus))
	                        endif
	                     endif

	                  else
	                     status = %loc(fip_csaread)
	                     call lib$signal (fip_csaread, %val(2),
     .					      fcc_infile(1:fcc_inlen),
     .					      %val(cstatus))
	                  endif	!  return status from read pixels
	               endif	!  outlist ne 0
	            endif	!  pstatus ne eof
	         enddo		!  read pixels


C
C  Close the skymaps.
C
	         status = fip_close_skymaps (flun, alun)

	      endif		!  return status from open skymaps

	   else
	      status = %loc(fip_ctinit)
	      call lib$signal (fip_ctinit, %val(1), %val(ct_stat(1)))
	   endif		!  return status from CT_INIT

	endif			!  return status from parse and report init

c
c  Update the processing report with the number of spectra processed.
c
	if ((status .eq. %loc(fip_normal))  .and.
     .	    (fcc_report .eq. fac_present)) then
	   status = fip_update_sky_report (nspec, npix)
	endif

c
c  Print out the number of spectra and pixels processed.
c
	type 20, nspec, npix
 20	format (/, 4x, 'Number of sky spectra reformated:     ', I5,
     .		/, 4x, 'Number of pixels containing spectra:  ', I5, /)

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
