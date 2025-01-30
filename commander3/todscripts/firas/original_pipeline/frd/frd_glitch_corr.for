	program frd_glitch_corr

c-------------------------------------------------------------------------------
c
c	Program FRD_GLITCH_CORR
c
c	This program read the Firas glitch rate correction slopes and
c	intercepts from the ASCII file FEX_GLTCHCOR.TXT and writes them into the
c	binary reference file FEX_GLTCHCOR.DAT, using the RDL FEX_GLTCHCOR.  The
c	file FEX_GLTCHCOR.DAT can then be CCL'ed into the Firas reference
c	archive.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  2 December 1993
c		  SER 11702
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
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'

	character * 14	gmt		!  current GMT
	character * 16	in_file  /'FEX_GLTCHCOR.TXT'/	!  input file name
	character * 16	out_file /'FEX_GLTCHCOR.DAT'/	!  output file name

	integer * 4	current_time(2)	!  current time
	integer * 4	in_lun		!  input file lun
	integer * 4	io_stat		!  I/O return status
	integer * 4	k		!  a counter
	integer * 4	out_lun		!  output file lun
	integer * 4	status		!  return status
	integer * 4	rstatus		!  return status

	real	* 8	glitchcor(12,2)	!  glitchrate correction slope and
					!    intercept

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_gltchcor'
	record /fex_gltchcor/ fex_glitch

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	frd_rmswrite

	status = %loc(frd_normal)

C
C  Read the glitchrate correction slope and intercept file.
C

c
c  Open the glitchrate correction file.
c
	rstatus = fut_get_lun (in_lun)
	if (.not. rstatus) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, name=in_file, status='old', iostat=io_stat,
     .		readonly, shared)

	if (io_stat .eq. 0) then
c
c  Read the glitchrate correction slopes and intercepts.
c
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,1),k=1,3)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,1),k=4,6)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,1),k=7,9)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,1),k=10,12)
	   read (in_lun, *, iostat=io_stat)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,2),k=1,3)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,2),k=4,6)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,2),k=7,9)
	   read (in_lun, *, iostat=io_stat) (glitchcor(k,2),k=10,12)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsread)
	      call lib$signal(frd_rmsread, %val(2), in_file, %val(io_stat))
	   endif

c
c  Close the correction file.
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
C  Fill the glitchrate correction slopes and intercepts into the RDL.
C
	   do k = 1,12
	      fex_glitch.slope(k)     = glitchcor(k,1)
	      fex_glitch.intercept(k) = glitchcor(k,2)
	   enddo

	   call sys$gettim(current_time)
	   call ct_binary_to_gmt(current_time,gmt)

	   do k = 1,2
	      fex_glitch.time(k) = current_time(k)
	   enddo
	   fex_glitch.gmt = gmt


C
C  Write the glitchrate correction RDL into the binary reference file.
C

c
c  Open the glitchrate correction reference file.
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
	      write (out_lun, iostat=io_stat) fex_glitch
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

	endif	!  (io_stat from glitchrate correction file read


c
c  Signal processing status of the program.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (%val(ss$_abort))
	endif

	end
