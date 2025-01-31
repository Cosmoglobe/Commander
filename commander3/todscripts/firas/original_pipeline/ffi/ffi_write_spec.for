	integer * 4 function ffi_write_spec (lun, num, fish_recs)

c-------------------------------------------------------------------------------
c
c	Function FFI_WRITE_SPEC
c
c	This function writes spectrum records to the FISH-format RMS output
c	file.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  29 September 1992
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Inputs:
c		lun		integer * 4		logical unit number
c		num		integer * 4		number of spectra in
c							current ensemble
c		fish_recs(num)				spectrum records
c
c	Outputs:
c		none
c
c	Subroutines called:
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
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffi_invoc)'
	include '(ffi_spec)'

	integer * 4	io_stat		!  write return status
	integer * 4	j		!  a counter
	integer * 4	lun		!  RMS file logical unit number
	integer * 4	num		!  number of spectra to be written
	integer * 4	status		!  return status

	record /fish_spec/ fish_recs(num)

	external	ffi_normal
	external	ffi_rmswrite


c
c  Write the voltage spectra to the FISH output file.
c
	status = %loc(ffi_normal)
	j = 1

	do while ((status .eq. %loc(ffi_normal))  .and.  (j .le. num))

	   write (unit=lun, iostat=io_stat) fish_recs(j)

	   if (io_stat .ne. 0) then
	      status = %loc(ffi_rmswrite)
	      call lib$signal (ffi_rmswrite, %val(2), fcc_outfile(1:fcc_outlen),
     .					     %val(io_stat))
	   endif

	   j = j + 1

	end do


	ffi_write_spec = status

	return
	end
