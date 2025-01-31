	program frd_extract_pixell
c-------------------------------------------------------------------------------
c
c	FRD_EXTRACT_PIXELL
c
c	Program to extract pixels of coadd records from multiple skymaps and
c	write them to a single output skymap.  The pixel numbers and skymaps
c	are specified in the file CSDR$FIRAS_REF:FRD_PIX_SCR_SS.TXT where
c       SS is the scan mode (SS, SF, or LF).
c
c	Author: S. Brodd, HSTX, 6/5/96, SER 12336
c-------------------------------------------------------------------------------

	implicit none

	include 'csdr$library:ctuser.inc'
	include '(csa_pixel_input_rec)'
	include '(csa_pixel_output_rec)'
	include '(fut_params)'
	include '(upm_stat_msg)'

	character * 33	file		!  RMS file name
	character * 18	in_ext(32)	!  skymap file ext
	character * 45	in_file		!  input skymap file name
	character * 44	out_file	!  output skymap file name

	integer * 2	ct_stat(20)	!  CT return status
	integer * 2	rec_size	!  coadd record size in
	parameter	(rec_size = fac_coad_spec_sizel / 4)  !  long words

	integer * 4	chan		!  channel number
	integer * 4	cstatus		!  CSA return status
	integer * 4	in_lun		!  input file logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	lun		!  RMS file lun
	integer * 4	npix		!  number of pixels to read
	integer * 4	num_in		!  number of pixels read by csa
	integer * 4	out_lun		!  output file logical unit number
	integer * 4	pix		!  pixel number
	integer * 4	pixel(32)	!  pixel numbers to be extracted
	integer * 4	pstatus		!  return status
	integer * 4	rstatus		!  return status
	integer * 4	smode		!  scan mode
	integer * 4	scr_smode	!  scan mode for script files
	integer * 4	status		!  return status

	integer * 4	csa_close_skymap
	integer * 4	csa_field_offset_values
	integer * 4	csa_open_skymap
	external	csa_open_skymap
	integer * 4	csa_read_pixels
	integer * 4	csa_write_pixels
	integer * 4	fut_get_lun
	integer * 4	fut_free_lun
	integer * 4	upm_present

	dictionary 'fil_sky'
	record /fil_sky/ coadd_recs(fac_max_coad)

	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist

	external	csa_normal
	external	frd_csaclose
	external	frd_csaopen
	external	frd_csaread
	external	frd_csawrite
	external	frd_ctinit
	external	frd_failure
	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
C
C  Parse the command line.
C
c
c  Set the channel specifier invocation flags.
c
	pstatus = upm_present ('CHANNEL')
	if (pstatus .eq. upm_pres) then
	   pstatus = upm_present ('CHANNEL.RH')
	   if (pstatus .eq. upm_pres) chan = 1
	   pstatus = upm_present ('CHANNEL.RL')
	   if (pstatus .eq. upm_pres) chan = 2
	   pstatus = upm_present ('CHANNEL.LH')
	   if (pstatus .eq. upm_pres) chan = 3
	   pstatus = upm_present ('CHANNEL.LL')
	   if (pstatus .eq. upm_pres) chan = 4
	else
	   chan = 1
	endif
c
c  Set the scan mode specifier invocation flags.
c
	pstatus = upm_present ('SCAN_MODE')
	if (pstatus .eq. upm_pres) then
	   pstatus = upm_present ('SCAN_MODE.SS')
	   if (pstatus .eq. upm_pres) then
	      smode = 1
              scr_smode = 1
	   endif
	   pstatus = upm_present ('SCAN_MODE.SF')
	   if (pstatus .eq. upm_pres) then
	      smode = 2
              scr_smode = 2
	   endif
	   pstatus = upm_present ('SCAN_MODE.LS')
	   if (pstatus .eq. upm_pres) then
	      smode = 3
              scr_smode = 3
	   endif
	   pstatus = upm_present ('SCAN_MODE.LF')
	   if (pstatus .eq. upm_pres) then
	      smode = 4
              scr_smode = 4
	   endif
	   pstatus = upm_present ('SCAN_MODE.FS')
	   if (pstatus .eq. upm_pres) then
	      smode = 5
              scr_smode = 2
	   endif
	   pstatus = upm_present ('SCAN_MODE.FL')
	   if (pstatus .eq. upm_pres) then
	      smode = 6
              scr_smode = 4
	   endif
	else
	      smode = 1
              scr_smode = 1
	endif
C
C  Read the pixel list file.
C
c
c  Open the file.
c
	status = %loc(frd_normal)
	file = 'csdr$firas_ref:frd_pix_scr_' // fac_scan_mode_idsl(scr_smode)
     .          // '.txt'
	rstatus = fut_get_lun(lun)
	open (unit=lun, name=file, status='old', form='formatted',
     .	      access='sequential', readonly, iostat=io_stat)
c
c  Read the list.
c
	if (io_stat .eq. 0) then
	   read (lun, *, iostat=io_stat) npix
	   do j = 1,npix
	      read (lun, 10, iostat=io_stat) pixel(j), in_ext(j)
	   enddo
   10	   format (x, i4, 2x, a18)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), file, %val(io_stat))
	   endif
c
c  Close the pixel list file.
c
	   close (lun, iostat=io_stat)
	   rstatus = fut_free_lun(lun)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal(frd_rmsclose, %val(2), file, %val(io_stat))
	   endif
	else
	   status = %loc(frd_rmsopen)
	   call lib$signal(frd_rmsopen, %val(2), file, %val(io_stat))
	endif

	if (status .eq. %loc(frd_normal)) then
C
C  Initialize Cobetrieve.
C
	   call ct_init(ct_stat)
	   if (ct_stat(1) .eq. ctp_normal) then
c
c Open the output skymap file.
c
	      rstatus = fut_get_lun(out_lun)
	      out_file = 'CSDR$FIRAS_OUT:FIL_SKY_' // fac_channel_ids(chan)
     .			  // fac_scan_mode_idsl(smode) // '.EP'
	      open (unit=out_lun, file=out_file, iostat=io_stat,
     .		    status='new', form='unformatted', recordtype='fixed',
     .		    recl=rec_size, useropen=csa_open_skymap)

	      if (io_stat .eq. 0) then
c
c  Set up the skymap offset values.
c
	         cstatus = csa_field_offset_values (fac_coad_spec_pix_offsetl,
     .						   fac_time_offset, -1, out_lun)
c
c  Loop over pixels and input skymaps.
c
	         pix = 0
	         do while ((status .eq. %loc(frd_normal))  .and.
     .			   (pix .lt. npix))
	            pix = pix + 1
	            inlist.level_no = fac_skymap_level
	            inlist.pixel_no = pixel(pix)
c
c  Open the input skymap file.
c
	            rstatus = fut_get_lun(in_lun)
	            in_file = 'CSDR$FIRAS_IN:FIL_SKY_' // fac_channel_ids(chan)
     .			       // fac_scan_mode_idsl(smode) // '.' //in_ext(pix)
	            open (unit=in_lun, file=in_file, iostat=io_stat,
     .			  status='old', form='unformatted',
     .			  recordtype='fixed', readonly,
     .			  useropen=csa_open_skymap)

	            if (io_stat .eq. 0) then
c
c  Read the coadd record for that pixel
c
	               cstatus = csa_read_pixels (in_lun, inlist, 1,
     .						     coadd_recs, fac_max_coad,
     .						     outlist, num_in, 0)

	               if (cstatus .eq. %loc(csa_normal)) then
	                  if (outlist.no_records .gt. 0) then
c
c  Write the coadd records to the output skymap file.
c
	                     cstatus = csa_write_pixels (out_lun, coadd_recs,
     .							 outlist.no_records, 0)
	                     if (cstatus .ne. %loc(csa_normal)) then
	                        status = %loc(frd_csawrite)
	                        call lib$signal (frd_csawrite, %val(2),
     .						 out_file, %val(cstatus))
	                     endif
	                  endif
	               else
	                  status = %loc(frd_csaread)
	                  call lib$signal (frd_csaread, %val(2), in_file,
     .							   %val(cstatus))
	               endif
c
c  Close the input skymap file.
c
	               cstatus = csa_close_skymap (in_lun, fac_skymap_no_levels)
	               if (cstatus .ne. %loc(csa_normal)) then
	                  status = %loc(frd_csaclose)
	                  call lib$signal (frd_csaclose, %val(2), in_file,
     .							    %val(cstatus))
	               endif
	            else
	               status = %loc(frd_csaopen)
	               call lib$signal (frd_csaopen, %val(2), in_file,
     .							%val(io_stat))
	            endif	!  (io_stat from input skymap open
	            rstatus = fut_free_lun(in_lun)
	         enddo		!  loop over pixels
c
c  Close the output skymap file.
c
	         cstatus = csa_close_skymap (out_lun, fac_skymap_no_levels)
	         if (cstatus .ne. %loc(csa_normal)) then
	            status = %loc(frd_csaclose)
	            call lib$signal (frd_csaclose, %val(2), out_file,
     .						       %val(cstatus))
	         endif
	      else
	         status = %loc(frd_csaopen)
	         call lib$signal (frd_csaopen, %val(2), out_file,
     .						   %val(io_stat))
	      endif	!  (io_stat from output skymap open
	      rstatus = fut_free_lun(out_lun)
	   else
	      status = %loc(frd_ctinit)
	      call lib$signal (frd_ctinit, %val(1), %val(ct_stat(1)))
	   endif	!  (ct_stat from ct_init
	endif		!  status from pixel list read
c
c  Signal the program completion status.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (frd_failure)
	endif

	end
