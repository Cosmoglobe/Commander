	integer * 4 function  ffi_open_coadd (ct_lun1, ct_lun2, vs_lun)

c-------------------------------------------------------------------------------
c
c	Function FFI_OPEN_COADD
c
c	This function opens the coadd archives CSDR$FIRAS_IN1 and
C	CSDR$FIRAS_IN2 for cal data.  It then opens the FISH-format spectrum
c	RMS file.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 March 1993
c		 SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none	
c
c	Output:
c		ct_lun1		integer * 4		coadd CT file lun
c		ct_lun2		integer * 4		coadd CT file lun
c		vs_lun		integer * 4		voltage spectra RMS
c							   file lun
c
c	Subroutines called:
c		ct_connect_read
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		ffi_invoc.txt
c		ffi_spec.txt
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

	include '(fut_params)'
	include '(ffi_invoc)'
	include '(ffi_spec)'

	integer * 4	ct_lun1		!  CT logical unit number
	integer * 4	ct_lun2		!  CT logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	rstatus		!  return status
	integer * 4	smode		!  mtm scan mode
	integer * 4	status		!  return status
	integer * 4	vs_lun		!  voltage spectra RMS file lun

	integer * 4	ct_connect_read
	external	ct_connect_read
	integer * 4	fut_get_lun

	external	ffi_ctifgopen
	external	ffi_normal
	external	ffi_rmsopen
	external	fut_normal

	status = %loc(ffi_normal)


C
C  Pick up the correct scan mode for low frequency long fast data.
C
	if (fcc_xllf .eq. fac_present) then
	   smode  = 4
	else
	   smode  = fcc_smode
	endif


C
C  Open the primary coadd file.
C

c
c  Get the logical unit number for Cobetrieve.
c
	rstatus = fut_get_lun (ct_lun1)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

c
c  Detemine the coadd file name.
c
	fcc_infile1 = 'CSDR$FIRAS_IN1:FIC_CAL_' // fac_channel_ids(fcc_chan)
     .						// fac_scan_mode_ids(smode) 
     .						// '/' // fcc_time_range
	call str$trim (fcc_infile1, fcc_infile1, fcc_inlen1)

c
c  Open the files.
c
	open (unit=ct_lun1, file=fcc_infile1, iostat=io_stat, status='old',
     .	      useropen=ct_connect_read)
	if (io_stat .ne. 0) then
 	   status = %loc(ffi_ctifgopen)
	   call lib$signal (ffi_ctifgopen, %val(2), fcc_infile1(1:fcc_inlen1),
     .					   %val(io_stat))
	endif


	if ((fcc_hybrid .eq. fac_present)  .and.
     .      (status .eq. %loc(ffi_normal))) then
C
C  Open the secondary coadd file for hybrid input.
C

c
c  Get the logical unit number for Cobetrieve.
c
	   rstatus = fut_get_lun (ct_lun2)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

c
c  Detemine the coadd file name.
c
	   fcc_infile2 = 'CSDR$FIRAS_IN2:FIC_CAL_' // fac_channel_ids(fcc_chan)
     .						   // fac_scan_mode_ids(smode) 
     .						   // '/' // fcc_time_range
	   call str$trim (fcc_infile2, fcc_infile2, fcc_inlen2)

c
c  Open the files.
c
	   open (unit=ct_lun2, file=fcc_infile2, iostat=io_stat, status='old',
     .	         useropen=ct_connect_read)
	   if (io_stat .ne. 0) then
 	      status = %loc(ffi_ctifgopen)
	      call lib$signal (ffi_ctifgopen, %val(2),
     .			       fcc_infile2(1:fcc_inlen2), %val(io_stat))
	   endif

	endif


	if (status .eq. %loc(ffi_normal)) then
C
C  Open the FISH-format RMS spectrum file.
C

c
c  Get the logical unit number for the RMS file.
c
	   rstatus = fut_get_lun (vs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

c
c  Detemine the RMS file name.
c
	   fcc_outfile = 'CSDR$FIRAS_OUT:FFI_' // fcc_scan_mode // '.'
     .					       // fcc_file_ext
	   call str$trim (fcc_outfile, fcc_outfile, fcc_outlen)

c
c  Open the file.
c
	   open (unit=vs_lun, file=fcc_outfile, access='sequential',
     .		 status='new', form='unformatted', recordtype='fixed',
     .		 recl=fcc_vrec_size, iostat=io_stat)

	   if (io_stat .ne. 0) then
              status = %loc(ffi_rmsopen)
	      call lib$signal (ffi_rmsopen, %val(2), fcc_outfile(1:fcc_outlen),
     .					    %val(io_stat))
	   end if

	endif		!  (status from coadd open


	ffi_open_coadd = status

	return
	end
