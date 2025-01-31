	integer * 4 function  fcf_write_cal_spec (flag, tmp_lun)

c-------------------------------------------------------------------------------
c
c	Function FCF_WRITE_CAL_SPEC
c
c	This function writes calibration voltage spectrum records and calibrated
c	spectrum records out to the Cobetrieve time-ordered archives.  The
c	spectra are read in from the temporary RMS files written by
c	FCF_Write_Temp_Spec.  The spectra are written out by the Cobetrieve
c	routine CT_Connect_Write.
c
c	Author:	  Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  21 April 1992
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
c		ct_close_arcv
c		ct_connect_write
c		ct_write_arcv
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
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

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fcf_invoc)'

	character * 48	out_file	!  Cobetrieve file name
	character * 52	temp_file	!  RMS file name

	integer * 2	ct_stat(20)	!  CT return status
	integer * 2	out_len		!  length of output file name
	integer * 2	temp_len	!  length of temporary file name

	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	flag		!  spectrum type flag
	integer * 4	io_stat		!  open return status
	integer * 4	num		!  number of spectra written out
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	tmp_lun		!  temporary RMS spectra file lun

	integer * 4	ct_connect_write
	external	ct_connect_write
	integer * 4	fut_get_lun
	integer * 4	fut_free_lun

	dictionary 'fcf_sky'
	record /fcf_sky/ spec_rec

	external	fcf_ctspecclose
	external	fcf_ctspecnum
	external	fcf_ctspecopen
	external	fcf_ctspecwrite
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
c Open the Cobetrieve spectrum file.
c
	   rstatus = fut_get_lun(ct_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   end if

	   open (unit=ct_lun, file=out_file, iostat=io_stat,
     .		   status='new', useropen=ct_connect_write)


C
C  Write the spectrum records to the Cobetrieve file.
C

	   if (io_stat .eq. 0) then
c
c  Read spectrum records in from the RMS file and write them out to the
c	Cobetrieve file.
c
	      ct_stat(1) = ctp_normal
	      num = 1
	      read (unit=tmp_lun, iostat=io_stat) spec_rec

	      do while ((io_stat .eq. 0)  .and.  (ct_stat(1) .eq. ctp_normal))

	         spec_rec.coad_spec_head.coadd_no = num
	         call ct_write_arcv (, ct_lun, spec_rec, ct_stat)

	         num = num + 1
		 read (unit=tmp_lun, iostat=io_stat) spec_rec

	      end do

c
c  Check the status of the reads and writes.
c
	      num = num - 1
	      if (fcc_nspec .ne. num) then
	         status = %loc(fcf_ctspecnum)
	         call lib$signal (fcf_ctspecnum, %val(3), %val(num),
     .				  %val(fcc_nspec), out_file(1:out_len))
	      endif

	      if ((io_stat .ne. 0)  .and.  (io_stat .ne. -1)) then
	         status = %loc(fcf_tempread)
	         call lib$signal (fcf_tempread, %val(2), temp_file(1:temp_len),
     .						%val(io_stat))
	      end if

	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fcf_ctspecwrite)
	         call lib$signal (fcf_ctspecwrite, %val(2), out_file(1:out_len),
     .						   %val(ct_stat(1)))
	      end if


C
C  Close the files.
C

c
c  Close the Cobetrieve spectrum file and the temporary RMS spectrum file.
c
	      call ct_close_arcv (,ct_lun, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fcf_ctspecclose)
	         call lib$signal (fcf_ctspecclose, %val(2), out_file(1:out_len),
     .						   %val(ct_stat(1)))
	      end if

	   else
	      status = %loc(fcf_ctspecopen)
	      call lib$signal (fcf_ctspecopen, %val(2), out_file(1:out_len),
     .					       %val(io_stat))
	   end if	!  (io_stat for Cobetrieve file

	   close (unit=tmp_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_tempclose)
	      call lib$signal (fcf_tempclose, %val(2), temp_file(1:temp_len),
     .					      %val(io_stat))
	   end if

	else
	   status = %loc(fcf_tempopen)
	   call lib$signal (fcf_tempopen, %val(2), temp_file(1:temp_len),
     .					  %val(io_stat))
	end if		!  (io_stat for RMS file

c
c  Free the logical unit numbers.
c
	rstatus = fut_free_lun(ct_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	end if

	rstatus = fut_free_lun(tmp_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	end if

	tmp_lun = 0


	fcf_write_cal_spec = status

	return
	end
