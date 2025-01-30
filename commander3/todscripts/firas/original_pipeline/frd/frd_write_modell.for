	integer * 4 function  frd_write_modell ()

c-------------------------------------------------------------------------------
c
c	Function FRD_WRITE_MODELL
c
c	This function writes the Firas calibration model solution contained in
c	the RDL FEX_MOD to the binary file FEX_MOD_CCSS.VVV_XXXXXXXXXX, where
c	CC is the channel, SS is the scan mode, VVV is the version of Dale
c	Fixsens program that generated the solution, and XXXXXXXXXX is the
c	label of the solution.
c
c	Author:
c		Gene Eplee
c		General Sciences Corp.
c		513-7768
c		2 October 1992
c		SER 8178
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		str$trim
c
c	Include files:
c		frd_model_invoc.txt
c		frd_model_solnl.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:  S. Brodd, HSTX, 12/12/95, SPR 12282.  Modifications for 
c                 improved FISH input.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(frd_model_invoc)'
	include '(frd_model_solnl)'

	character * 48	out_file		!  output file name

	integer *  2	out_len			!  output file name length

	integer *  4	io_stat			!  I/O status
	integer *  4	out_lun			!  output file lun
	integer *  4	rstatus			!  return status
	integer *  4	status			!  return status

	integer *  4	fut_free_lun
	integer *  4	fut_get_lun

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmswrite
	external	fut_normal

	status = %loc(frd_normal)

c
c  Determine the output model parameter file name.
c
	out_file = 'csdr$firas_out:fex_mod_' // fac_channel_ids(fcc_chan) //
     .				   fac_scan_mode_idsl(fcc_smode) // '.' //
     .				   fcc_infile_ext
	call str$trim (out_file, out_file, out_len)

c
c  Open the model parameter file.
c
	rstatus = fut_get_lun(out_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal(rstatus)
	endif

	open (unit=out_lun, file=out_file, status='new', form='unformatted',
     .        recl=7424, recordtype='fixed', access='sequential',
     .	      iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Write the calibration model to the model parameter file.
c
	   write (out_lun, iostat=io_stat) fex_model

	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmswrite)
	      call lib$signal (frd_rmswrite, %val(2), out_file(1:out_len),
     .					     %val(io_stat))
	   endif

c
c  Close the model parameter file.
c
	   close (out_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal (frd_rmsclose, %val(2), out_file(1:out_len),
     .					     %val(io_stat))
	   endif
	   rstatus = fut_free_lun(out_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal(rstatus)
	   endif

	else		!  (io_stat from opening of model
	   status = %loc(frd_rmsopen)
	   call lib$signal (frd_rmsopen, %val(2), out_file(1:out_len),
     .					 %val(io_stat))
	endif


	frd_write_modell = status

	return
	end
