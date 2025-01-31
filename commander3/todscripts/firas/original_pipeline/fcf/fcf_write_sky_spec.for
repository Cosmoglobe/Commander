	integer * 4 function  fcf_write_sky_spec (flag, tmp_lun)

c-------------------------------------------------------------------------------
c
c	Function FCF_WRITE_SKY_SPEC
c
c	This function writes sky voltage spectrum records and calibrated
c	spectrum records out to the Cobetrieve skymap archives.  The spectra are
c	read in from the temporary RMS files written by FCF_Write_Temp_Spec.
c	The spectra are written out by the CSA routine CSA_Write_Pixels.
c
c	Author:	  Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  13 March 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		flag		integer * 4		spectrum type flag
c							1 = voltage
c							2 = differential
c							3 = calibrated
c		tmp_lun		integer * 4		RMS spectra file lun
c
c	Output:
c		none
c
c	Subroutines called:
c		csa_close_skymap
c		csa_field_offset_values
c		csa_open_skymap
c		csa_write_pixels
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_invoc)'

	character * 52	temp_file	!  RMS file name
	character * 48	out_file	!  Cobetrieve file name

	integer * 2	out_len		!  length of skymap file name
	integer * 2	temp_len	!  length of temporary file name

	integer * 4	blocks		!  block count for csa_write_pixels
	integer * 4	cstatus		!  CSA return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	flag		!  spectrum type flag
	integer * 4	io_stat		!  open return status
	integer * 4	num		!  number of spectra written out
	integer * 4	num_out		!  number of records to be written out
					!    in one call to csa
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	tmp_lun		!  temporary RMS spectra file lun

	integer * 4	csa_close_skymap
	integer * 4	csa_field_offset_values
	integer * 4	csa_open_skymap
	external	csa_open_skymap
	integer * 4	csa_write_pixels
	integer * 4	fut_get_lun
	integer * 4	fut_free_lun

	dictionary 'fcf_sky'
	record /fcf_sky/ spec_rec

	external	csa_normal
	external	fcf_csafldoffset
	external	fcf_csaspecclose
	external	fcf_csaspecnum
	external	fcf_csaspecopen
	external	fcf_csaspecwrite
	external	fcf_normal
	external	fcf_tempclose
	external	fcf_tempopen
	external	fcf_tempread
	external	fut_normal

C
C  Open the files.
C

	status = %loc(fcf_normal)

c
c  Determine the file names.
c
	if (flag .eq. 1) then
	   temp_file = fcc_temp_outfile_vs
	   temp_len  = fcc_toutlenvs
	   out_file  = fcc_outfile_vs
	   out_len   = fcc_outlenvs
	elseif (flag .eq. 2) then
	   temp_file = fcc_temp_outfile_ds
	   temp_len  = fcc_toutlends
	   out_file  = fcc_outfile_ds
	   out_len   = fcc_outlends
	elseif (flag .eq. 3) then
	   temp_file = fcc_temp_outfile_cs
	   temp_len  = fcc_toutlencs
	   out_file  = fcc_outfile_cs
	   out_len   = fcc_outlencs
	endif

c
c  Open the temporary RMS spectrum file.
c
	open (unit=tmp_lun, file=temp_file, access='sequential',
     .	     status='old', form='unformatted', dispose='delete', iostat=io_stat)

	if (io_stat .eq. 0) then
c
c Open the Cobetrieve skymap file.
c
	   rstatus = fut_get_lun(ct_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=ct_lun, file=out_file, iostat=io_stat, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_record_size,
     .		 useropen=csa_open_skymap)

	   if (io_stat .eq. 0) then
c
c  Set up the skymap offset values.
c
	      cstatus = csa_field_offset_values (fac_coad_spec_pix_offset,
     .						 fac_time_offset, -1, ct_lun)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fcf_csafldoffset)
	         call lib$signal (fcf_csafldoffset, %val(2),
     .				  out_file(1:out_len), %val(cstatus))
	      endif


C
C  Write the spectrum records to the Cobetrieve file.
C

c
c  Read spectrum records in from the RMS file and write them out to the
c	skymap file.
c
	      blocks = 0
	      num = 1
	      num_out=1

	      read (unit=tmp_lun, iostat=io_stat) spec_rec

	      do while ((cstatus .eq. %loc(csa_normal)) .and. (io_stat .eq. 0))

	         spec_rec.coad_spec_head.coadd_no = num
	         cstatus = csa_write_pixels (ct_lun, spec_rec, num_out, blocks)

	         num = num + 1
		 read (unit=tmp_lun, iostat=io_stat) spec_rec

	      enddo

c
c  Check the status of the reads and writes.
c
	      num = num - 1
	      if (fcc_nspec .ne. num) then
	         status = %loc(fcf_csaspecnum)
	         call lib$signal (fcf_csaspecnum, %val(3), %val(num),
     .				  %val(fcc_nspec), out_file(1:out_len))
	      endif

	      if ((io_stat .ne. 0)  .and.  (io_stat .ne. -1)) then
	         status = %loc(fcf_tempread)
	         call lib$signal (fcf_tempread, %val(2), temp_file(1:temp_len),
     .						%val(io_stat))
	      endif

	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fcf_csaspecwrite)
	         call lib$signal (fcf_csaspecwrite, %val(2),
     .				  out_file(1:out_len), %val(cstatus))
	      endif


C
C  Close the files.
C

c
c  Close the Cobetrieve spectrum file and the temporary RMS spectrum file.
c
	      cstatus = csa_close_skymap (ct_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fcf_csaspecclose)
	         call lib$signal (fcf_csaspecclose, %val(2),
     .				  out_file(1:out_len), %val(cstatus))
	      endif

	   else
	      status = %loc(fcf_csaspecopen)
	      call lib$signal (fcf_csaspecopen, %val(2), out_file(1:out_len),
     .						%val(io_stat))
	   endif	!  (io_stat for skymap file

	   close (unit=tmp_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_tempclose)
	      call lib$signal (fcf_tempclose, %val(2), temp_file(1:temp_len),
     .					      %val(io_stat))
	   endif

	else
	   status = %loc(fcf_tempopen)
	   call lib$signal (fcf_tempopen, %val(2), temp_file(1:temp_len),
     .					  %val(io_stat))
	endif		!  (io_stat for RMS file

c
c  Free the logical unit numbers.
c
	rstatus = fut_free_lun(ct_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	rstatus = fut_free_lun(tmp_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	tmp_lun = 0


	fcf_write_sky_spec = status

	return
	end
