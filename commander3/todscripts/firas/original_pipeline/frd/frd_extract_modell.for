	Program FRD_Extract_Modell

c-------------------------------------------------------------------------------
c
c	FRD_EXTRACT_MODELL
c
c	The program FRD_Extract_Modell reads Dale Fixsens Firas calibration
c	model solution from the RMS file FEX_MOD_CCSS.TXT, where CC is the
c	channel and SS is the scan mode.  The model solution is renormalized for
c	pipeline processing and inserted into the RDL FEX_MOD.  It is then
c	written to the binary file FEX_MOD_CCSS.VVV_XXXXXXXXXX, where VVV is the
c	version of Dale Fixsens program that generated the solution and
c	XXXXXXXXXX is the label of the solution.  The binary file can then be
C	CCLed into the Firas Calibration Reference Archive.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  30 July 1992
c		  SER 8178
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		cct_close_config
c		cct_get_config_tod
c		cct_open_config
c		ct_init
c		frd_parse_modell
c		frd_read_modell
c		frd_write_modell
c
c	Include files:
c		csdr$library:ctuser.inc
c		frd_model_configl.txt
c		frd_model_solnl.txt
c
c-------------------------------------------------------------------------------
c
c	Changes
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 22 October 1993
c	SER 11394
c
c       Modifications for 361 point spectra.
c       Steve Brodd, HSTX, 12 December 1995, SPR 12282.
c
c-------------------------------------------------------------------------------

	implicit none

	include 'csdr$library:ctuser.inc'
	include '(frd_model_configl)'
	include '(frd_model_solnl)'

	integer * 2	ct_stat(20)	!  Cobetrieve return status

	integer * 4	cstatus		!  Get_Config return status
	integer * 4	status		!  return status

	integer * 4	cct_close_config
	integer * 4	cct_get_config_tod
	integer * 4	cct_open_config
	integer * 4	frd_parse_modell
	integer * 4	frd_read_modell
	integer * 4	frd_write_modell

	external	cct_normal
	external	frd_closeconfig
	external	frd_ctinit
	external	frd_failure
	external	frd_getconfig
	external	frd_normal
	external	frd_openconfig

c
c  Parse the command line.
c
	status = frd_parse_modell ()

	if (status .eq. %loc(frd_normal)) then
c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)

	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Get the reference dataset.
c
	      cstatus = cct_open_config (ref_start, ref_stop, ndset_tod,
     .					 dset_tod, size_tod, seq_access, ncache,
     .					 tod_lun, tod_index, tod_stat,
     .					 ref_count)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(frd_openconfig)
	         call lib$signal (frd_openconfig, %val(1), %val(cstatus))
	      else
	         cstatus = cct_get_config_tod (ref_time, ndset_tod, size_tod,
     .					       tod_lun, tod_index, nyquistl,
     .				               new_tod_segment, tod_stat)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(frd_getconfig)
	            call lib$signal (frd_getconfig, %val(1), %val(cstatus))
	         endif
	         cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(frd_closeconfig)
	            call lib$signal (frd_closeconfig, %val(1), %val(cstatus))
	         endif
	      endif


	      if (status .eq. %loc(frd_normal)) then
c
c  Initialize the model solution RDL.
c
	         call lib$movc5(0,,0,29696,fex_model)

c
c  Read the calibration model solution from the RMS file.
c
	         status = frd_read_modell ()

	         if (status .eq. %loc(frd_normal)) then
c
c  Write the calibration model solution to the binary file.
c
	            status = frd_write_modell ()

	         endif

	      endif	!  status from get_config

	   else
	      status = %loc(frd_ctinit)
	      call lib$signal (frd_ctinit, %val(1), %val(ct_stat(1)))
	   endif	!  status from ct_init

	endif		!  status from command line parse


c
c  Terminate the program.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (frd_failure)
	endif

	end
