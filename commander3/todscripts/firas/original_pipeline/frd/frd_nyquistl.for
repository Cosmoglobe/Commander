	program frd_nyquistl

c-------------------------------------------------------------------------------
c
c	Program FRD_NYQUISTL
c
c	This program reads the Firas MTM sampling rates from the ASCII file
c       FEX_SAMPRATE.TXT and the optical Nyquist frequency correction from the 
C       ASCII file FEX_NYQUIST.TXT.  It then computes the Nyquist frequencies 
c       in icm and hz for all channels and scan modes for either I&T data or 
c       flight data.  These Nyquist frequencies are written into the binary 
c       reference file FEX_NYQUISTL.DAT, using the RDL FEX_NYQUISTL.  The file 
c	FEX_NYQUISTL.DAT can then be CCLed into the FIRAS reference archive.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  25 June 1992
c		  SER 9790
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
c	Changes:  Added scan modes 'FS' and 'FL'.  Changed determination of
c                 Nyquist frequency in Hz to use FFT length from FUT_PARAMS
c                 file.  THIS ROUTINE DOES NOT RETURN THE CORRECT NYQUIST
c                 FREQUENCY IN HZ FOR SCAN MODES "FS" AND "FL" FOR FAKEIT DATA.
c                 (The frequency in Hz returned should be divided by 4 in
c                 those special cases.)
c                 Alice Trenholme, GSC, 2/95, SER 12244.
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(upm_stat_msg)'
        include '(fut_params)'    !necessary to pull up spectrum length

	character * 16	in_file1  /'FEX_SAMPRATE.TXT'/	!  input file name
	character * 16	in_file2  /'FEX_NYQUIST.TXT'/	!  input file name
	character * 16	out_file  /'FEX_NYQUISTL.DAT'/	!  output file name

	integer * 4	in_lun				!  input file lun
	integer * 4	io_stat				!  I/O return status
	integer * 4	j				!  a counter
	integer * 4	mtm_speed(10)			!  commanded mtm scan
	data mtm_speed	/ 1, 2, 1, 2, 1, 2, 1, 2 , 2, 2/!    speeds + 1
	integer * 4	out_lun				!  output file lun
	integer * 4	pstatus				!  parse return status
	integer * 4	status				!  return status
	integer * 4	rstatus				!  return status

	logical * 4	int /.false./			!  I&T data flag
	logical * 4	flight /.false./		!  flight data flag

	real	* 4	fakeit(5)			!  fakeit adds per group
	data fakeit	/ 1.0, 2.0, 3.0, 8.0, 12.0 /
	real	* 4	fringes				!  mtm grating spacing
	parameter	(fringes = 20.00e-04)		!    in cm
	real	* 4	int_rate			!  I&T sampling rate
	real	* 4	flight_rate			!  flight sampling rate
        real    * 4     sampling_rate			!  sampling rate
	real	* 4	freq_shift			!  optical Nyquist
							!   frequency correction
	real	* 4	multiplier(2)			!  fringe multiplier for
	data multiplier	/ 6.0, 4.0 /			!    mtm scan speeds
	real	* 4	ngroup(10)			!  normal adds per group
	data ngroup	/ 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 12.0, 8.0, 8.0, 8.0 /
	real	* 4	norm				!  Nyquist frequency
							!    normalization
	real	* 4	optical_path			!  FIRAS optical path
	parameter	(optical_path = 3.4641016)	!    = 4*cos(pi/6)
	real	* 4	speed(2)			!  actual mtm scan
							!    speeds

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun
	integer * 4	upm_present

	dictionary 'fex_nyquistl'
	record /fex_nyquistl/ fex_nyquistl

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	frd_rmswrite
	external	fut_normal

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
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, name=in_file1, status='old', iostat=io_stat,
     .	      readonly, shared)

	if (io_stat .eq. 0) then
c
c  Read the sampling rates.
c
	   read (in_lun, *, iostat=io_stat) int_rate
	   read (in_lun, *, iostat=io_stat) flight_rate
	   if (io_stat .eq. 0) then
	      if (int .eq. .true.) sampling_rate = int_rate
	      if (flight .eq. .true.) sampling_rate = flight_rate
	   else
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), in_file1, %val(io_stat))
	   endif

c
c  Close the sampling rate file.
c
	   close (in_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal(frd_rmsclose, %val(2), in_file1, %val(io_stat))
	   endif

	   rstatus = fut_free_lun (in_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal(rstatus)
	   endif

	else
	   status = %loc(frd_rmsopen)
	   call lib$signal(frd_rmsopen, %val(2), in_file1, %val(io_stat))
	endif

C
C  Read the frequency correction file.
C

c
c  Open the frequency correction file.
c
	rstatus = fut_get_lun (in_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, name=in_file2, status='old', iostat=io_stat,
     .	      readonly, shared)

	if (io_stat .eq. 0) then
c
c  Read the frequency correction.
c
	   read (in_lun, *, iostat=io_stat) freq_shift

	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), in_file2, %val(io_stat))
	   endif

c
c  Close the frequency correction file.
c
	   close (in_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal(frd_rmsclose, %val(2), in_file2, %val(io_stat))
	   endif

	   rstatus = fut_free_lun (in_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal(rstatus)
	   endif

	else
	   status = %loc(frd_rmsopen)
	   call lib$signal(frd_rmsopen, %val(2), in_file2, %val(io_stat))
	endif

	if (status .eq. %loc(frd_normal)) then
C
C  Calculate the Nyquist frequencies.
C

c
c  Calculate the optical Nyquist frequencies.
c
	   norm = freq_shift / fringes / optical_path / 2.0
	   do j = 1,10
	      fex_nyquistl.icm(j) = norm * multiplier(mtm_speed(j)) / ngroup(j)
	   enddo

c
c  Calculate the scan speeds.
c
	   speed(1) = sampling_rate/fex_nyquistl.icm(1)/multiplier(1)
	   speed(2) = sampling_rate/fex_nyquistl.icm(1)/multiplier(2)

c
c  Calculate the electronic Nyquist frequencies.
c
	   do j = 1,10
	      fex_nyquistl.hz(j) = fex_nyquistl.icm(j) * speed(mtm_speed(j))
	   enddo

c
c  Calculate the fakeit Nyquist frequencies.
c
	   do j = 11, 15
	      fex_nyquistl.hz(j) = dble(fac_fft_length(1))/ fakeit(j-10)
	   enddo


C
C  Write the Nyquist frequencies into the binary reference file.
C

c
c  Open the Nyquist frequency reference file.
c
	   rstatus = fut_get_lun (out_lun)
	   if (.not. rstatus) then
	      call lib$signal(rstatus)
	   endif

	   open (unit=out_lun, file=out_file, status='new', recl=32,
     .		 recordtype='fixed', form='unformatted', access='sequential',
     .		 iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the Nyquist frequencies to the binary file.
c
	      write (out_lun, iostat=io_stat) fex_nyquistl
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

	endif	!  (io_stat from sampling rate file read


c
c  Signal processing status of the program.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (%val(ss$_abort))
	endif

	end
