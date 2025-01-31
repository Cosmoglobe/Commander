	program frd_vibcorrl

c-------------------------------------------------------------------------------
c
c	Program FRD_VIBCORRL
c
c	This program read the Firas vibration correction frequency offsets from
c	the ASCII file FEX_VIBCORRL.TXT and writes them into the binary reference
c	file FEX_VIBCORRL.DAT, using the RDL FEX_VIBCORRL.  The file 
c	FEX_VIBCORRL.DAT can then be CCL'ed into the Firas reference archive.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  16 January 1992
c                 SER 8292
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		ct_binary_to_gmt
c		fut_free_lun
c		fut_get_lun
c		sys$gettim
c
c	Include files:
c		$ssdef
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Added low frequency short fast offsets.
c	Gene Eplee, GSC, 17 August 1993
c
c	Changed for floating-point offsets for 341-point spectra.
c	Gene Eplee, GSC, 1 July 1994
c
c       Changed to read floating-point offsets for 361-point spectra and if
c       actual FFT length is less than 720, adjust the size of the offset.
c       Uses FFT length from FUT_PARAMS file.
c       Alice Trenholme, GSC, 3 February 1995, SPR 12287.
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
        include '(fut_params)'

	character * 14	gmt		!  current GMT
	character * 16	in_file  /'FEX_VIBCORRL.TXT'/	!  input file name
	character * 16	out_file /'FEX_VIBCORRL.DAT'/	!  output file name

	integer * 4	current_time(2)	!  current time
	integer * 4	in_lun		!  input file lun
	integer * 4	io_stat		!  I/O return status
	integer * 4	k		!  a counter
	integer * 4	out_lun		!  output file lun
	integer * 4	status		!  return status
	integer * 4	rstatus		!  return status

	real 	* 4	vibcorr(12,2)	!  vibration correction frequency offset
        real    * 4     corr_fact       !  adjustment for variable length FFT

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_vibcorrl'
	record /fex_vibcorrl/ fex_vib

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	frd_rmswrite

	status = %loc(frd_normal)

C
C  Read the vibration correction frequency offset file.
C

c
c  Open the vibration correction file.
c
	rstatus = fut_get_lun (in_lun)
	if (.not. rstatus) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, name=in_file, status='old', iostat=io_stat,
     .		readonly, shared)

	if (io_stat .eq. 0) then
c
c  Read the vibration correction frequency offsets.
c
	   read (in_lun, *, iostat=io_stat) (vibcorr(k,1),k=1,6)
	   read (in_lun, *, iostat=io_stat) (vibcorr(k,1),k=7,12)
	   read (in_lun, *, iostat=io_stat) (vibcorr(k,2),k=1,6)
	   read (in_lun, *, iostat=io_stat) (vibcorr(k,2),k=7,12)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), in_file, %val(io_stat))
	   endif

c
c  Close the offset file.
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


	if (status .eq. %loc(frd_normal)) then
C
C  Fill the vibration correction frequency offsets into the RDL.  
C
c
c  Since the vibration corrections in the file are based on an FFT length
c  of 720 compute a correction factor and adjust if the actual FFT length
c  is less.  This formula works for the special low FS and FL scan modes also
c  because the proportionality constant is the same regardless (i.e.,
c  fac_fft_length(5 or 6)/180.d0 = fac_fft_length(1)/720.d0.
c
           corr_fact = fac_fft_length(1)/720.d0
	   do k = 1,12
	      fex_vib.primary_offset(k)  = vibcorr(k,1)*corr_fact
	      fex_vib.secondary_offset(k) = vibcorr(k,2)*corr_fact
	   enddo

	   call sys$gettim(current_time)
	   call ct_binary_to_gmt(current_time,gmt)

	   fex_vib.ct_head.gmt = gmt
	   do k = 1,2
	      fex_vib.ct_head.time(k) = current_time(k)
	   enddo


C
C  Write the vibration correction RDL into the binary reference file.
C

c
c  Open the vibration correction reference file.
c
	   rstatus = fut_get_lun (out_lun)
	   if (.not. rstatus) then
	      call lib$signal(rstatus)
	   endif

	   open (unit=out_lun, file=out_file, status='new', recl=64,
     .		 recordtype='fixed', form='unformatted', access='sequential', 
     .           iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the RDL to the binary file.
c
	      write (out_lun, iostat=io_stat) fex_vib
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

	endif	!  (io_stat from vibration correction frequency file read


c
c  Signal processing status of the program.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (%val(ss$_abort))
	endif

	end
