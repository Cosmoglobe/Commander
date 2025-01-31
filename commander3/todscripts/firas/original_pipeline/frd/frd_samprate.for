	program frd_samprate

c-------------------------------------------------------------------------------
c
c	Program FRD_SAMPRATE
c
c	This program reads the Firas I&T or mission sampling rate from the ASCII
c       file FEX_SAMPRATE.TXT and writes it into the binary reference file 
c       FEX_SAMPRATE.DAT, using the RDL FEX_SAMPRATE.  This file can then be 
c       CCL'ed into the Firas reference archive.
c
c	Author:   Steven Alexander, HSTX, August 3, 1992, SER 9859.
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c
c	Include files:
c		$ssdef
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(upm_stat_msg)'

	character * 16	in_file  /'FEX_SAMPRATE.TXT'/	!  input file name
	character * 16	out_file /'FEX_SAMPRATE.DAT'/	!  output file name

	integer * 4	in_lun				!  input file lun
	integer * 4	io_stat				!  I/O return status
	integer * 4	out_lun				!  output file lun
	integer * 4	pstatus				!  parse return status
	integer * 4	status				!  return status
	integer * 4	rstatus				!  return status

	logical * 4	int /.false./			!  I&T data flag
	logical * 4	flight /.false./		!  flight data flag

	real	* 4	int_rate			!  I&T sampling rate
	real	* 4	flight_rate			!  flight sampling rate

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun
	integer * 4	upm_present

	dictionary 'fex_samprate'
	record /fex_samprate/ fex_samprate

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	frd_rmswrite

C
C  Parse the command line.
C

	status = %loc(frd_normal)
	pstatus = upm_present ('SAMPLE')
	if (pstatus .eq. upm_pres) then
           pstatus = upm_present ('SAMPLE.INT')
	   if (pstatus .eq. upm_pres) int = .true.
           pstatus = upm_present ('SAMPLE.MISSION')
	   if (pstatus .eq. upm_pres) flight = .true.
	else
	   flight = .true.
	endif	


C
C  Read the sampling rate file.
C

c
c  Open the sampling rate file.
c
	rstatus = fut_get_lun (in_lun)
	if (.not. rstatus) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, name=in_file, status='old', iostat=io_stat,
     .	      readonly, shared)

	if (io_stat .eq. 0) then
c
c  Read the sampling rates.
c
	   read (in_lun, *, iostat=io_stat) int_rate
	   read (in_lun, *, iostat=io_stat) flight_rate
	   if (io_stat .eq. 0) then
	      if (int .eq. .true.) fex_samprate.sampling_rate = int_rate
	      if (flight .eq. .true.) fex_samprate.sampling_rate = flight_rate
	   else
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), in_file, %val(io_stat))
	   endif

c
c  Close the sampling rate file.
c
	   close (in_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal(frd_rmsclose, %val(2), in_file, %val(io_stat))
	   endif

	   rstatus = fut_free_lun (in_lun)
	   if (.not. rstatus) then
	      call lib$signal(rstatus)
	   endif

	else
	   status = %loc(frd_rmsopen)
	   call lib$signal(frd_rmsopen, %val(2), in_file, %val(io_stat))
	endif

C
C  Write the sampling rate into the binary reference file.
C

c
c  Open the sampling rate reference file.
c
	rstatus = fut_get_lun (out_lun)
	if (.not. rstatus) then
	   call lib$signal(rstatus)
	endif

	open (unit=out_lun, file=out_file, status='new', recl=16,
     .	      recordtype='fixed', form='unformatted', access='sequential',
     .	      iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Write the sampling rate to the binary file.
c
	   write (out_lun, iostat=io_stat) fex_samprate
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmswrite)
	      call lib$signal(frd_rmswrite, %val(2), out_file, %val(io_stat))
	   endif

c
c  Close the binary file.
c
	   close (out_lun, iostat = io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal(frd_rmsclose, %val(2), out_file, %val(io_stat))
	   endif
 
           rstatus = fut_free_lun (out_lun)
           if (.not. rstatus) then
              call lib$signal(rstatus)
           endif

        else
           status = %loc(frd_rmsopen)
           call lib$signal(frd_rmsopen, %val(2), out_file, %val(io_stat))
        endif
c
c  Signal processing status of the program.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (%val(ss$_abort))
	endif

	end
