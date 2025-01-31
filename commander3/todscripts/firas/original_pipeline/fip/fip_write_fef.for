	integer * 4 function   fip_write_fef (num, out_recs)

c-------------------------------------------------------------------------------
c
c	Function FIP_WRITE_FEF
c
c	This function write the FIRAS Extra Factor PDS records for FIP_FEF.
c
c	Author:	  Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  9 September 1994
c
c-------------------------------------------------------------------------------
c
c	Input:
c		num		integer * 4	number of records to write
c		out_recs	FIP records	FIP records
c
c	Output:
c		none
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

	integer * 4	fip_lun		!  FIP logical unit number
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	num		!  number of records to write
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fip_fef'
	record /fip_fef/ out_recs(fac_max_spec)

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmswrite
	external	fut_normal

C
C  Write the records to the FIP file.
C
	status = %loc(fip_normal)

c
c  Find the file name.
c
	fcc_fip_file = 'CSDR$FIRAS_OUT:FIP_FEF_' // fcc_scan_mode // '.' //
     .			fcc_file_ext
	call str$trim (fcc_fip_file, fcc_fip_file, fcc_fiplen)

c
c  Open the file.
c
	rstatus = fut_get_lun(fip_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fip_lun, file=fcc_fip_file, status='new',
     .	      form='unformatted', recordtype='fixed', recl=fcc_fip_size,
     .	      access='sequential', iostat=io_stat)


	if (io_stat .eq. 0) then
C
C  Write the data.
C
	   j = 0
	   do while ((j .lt. num)  .and.  (io_stat .eq. 0))
	      j = j+1
	      write (fip_lun, iostat=io_stat) out_recs(j)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmswrite)
	         call lib$signal (fip_rmswrite, %val(2),
     .				  fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	      endif
	   enddo



C
C  Close the error file.
C
	   close (unit=fip_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fip_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2),
     .			    fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	endif	!  (open status

	fip_write_fef = status

	return
	end
