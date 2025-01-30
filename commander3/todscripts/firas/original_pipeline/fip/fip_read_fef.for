	integer * 4 function   fip_read_fef (num, in_recs)

c-------------------------------------------------------------------------------
c
c	Function FIP_READ_FEF
c
c	This function reads the FIRAS Extra Factor records for FIP_FEF.
c
c	Author:	  Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  9 September 1994
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		num		integer * 4		number of records read
c		in_recs	FEF records		FEF records
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c		str$trim
c
c	Include files:
c		fip_invoc_fef.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_fef)'

	integer * 4	fef_lun		!  FEF logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	num		!  number of records read
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fef_spc'
	record /fef_spc/ in_rec
	record /fef_spc/ in_recs(fac_max_spec)

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmsread
	external	fut_normal

C
C  Read the records from the FEF file.
C

c
c  Find the file name.
c
	fcc_fef_file = 'CSDR$FIRAS_IN:FEF_SPC_' // fcc_scan_mode // '.' //
     .			fcc_file_ext
	call str$trim (fcc_fef_file, fcc_fef_file, fcc_feflen)

c
c  Open the file.
c
	rstatus = fut_get_lun(fef_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fef_lun, file=fcc_fef_file, status='old',
     .	      form='unformatted', access='sequential', readonly, iostat=io_stat)


	if (io_stat .eq. 0) then
C
C  Read the data.
C
	   num = 0
	   do while ((io_stat .eq. 0)  .and. (num .lt. fac_max_spec))

c
c  Read the records.
c
	      read (fef_lun, iostat=io_stat) in_rec

	      if (io_stat .eq. 0) then
	         num = num + 1
	         in_recs(num) = in_rec
	      endif

	   enddo	!  while (io_stat

c
c  Close the error file.
c
	   close (unit=fef_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fef_file(1:fcc_feflen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fef_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif


C
C  Check the status of the read.
C
	   if (io_stat .le. 0) then
	      status = %loc(fip_normal)
	   else
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2),
     .			       fcc_fef_file(1:fcc_feflen), %val(io_stat))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2),
     .			    fcc_fef_file(1:fcc_feflen), %val(io_stat))
	endif		!  status from open


	fip_read_fef = status

	return
	end
