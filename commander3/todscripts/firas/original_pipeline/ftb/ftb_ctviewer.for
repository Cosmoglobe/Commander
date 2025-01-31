	Program FTB_CTVIEWER

c------------------------------------------------------------------------------
c    PURPOSE: Produce plots of interferograms and spectrum from any
c	      dataset in the archive.
c
c       REQUIREMENTS REFERENCE NUMBERS:
c
c       INPUT DATA:
c
c       OUTPUT DATA:
c
c    AUTHOR: Wes Young
c            SASC Technologies
c            August, 1986
c
c    INVOCATION:
c
c    INPUT PARAMETERS:
c
c    OUTPUT PARAMETERS:
c
c    SUBROUTINES CALLED:
c
c    COMMON VARIABLES USED:
c
c    INCLUDE FILES:
c
c    PROCESSING METHOD:
c
c----------------------------------------------------------------------
c
c     Change Log:
c
c	R. Kummerer, Nov 26, 1986. Completed coding so that CTVIEWER runs.
c
c	R. Kummerer, Dec 23, 1986. Fix CT open error complaint on exit.
c
c	R. Kummerer, Dec 24, 1986. Add feature to display coadd data.
c
c	R. Kummerer, Jan 8, 1987. Add display of spectra.
c
c	R. Kummerer, Jan 16, 1987. FUT_DISPLAY zero-point argument.
c
c	R. Kummerer, Jan 20, 1987. Zoom plotted wrong set of points.
c
c	R. Kummerer, Jan 22, 1987. CDD conversion.
c
c	R. Kummerer, Jan 29, 1987. Fix null channel plotting.
c
c	R. Kummerer, Feb 17, 1987. FUT_DISPLAY/LABEL device type argument.
c
c	R. Kummerer, Mar 13, 1987. Multi-device plot capability.
c
c	R. Kummerer, Mar 23, 1987. Plot label corrections.
c
c	R. Kummerer, March 31, 1987. Plot BB curve with cal spectra.
c
c	R. Kummerer, Apr 9, 1987. Update FUT_DISPLAY_MODEL argument list.
c
c	R. Kummerer, July 23, 1987. Display phase correctors.
c
c	R. Kummerer, August 6, 1987. Remove USTART,UEND, in FUT_DISPLAY now.
c
c	R. Kummerer, August 20, 1987. Display new model parameters...
c	                              model covariance and solution residuals.
c
c	F. Shuman, 1988 Feb  5.  Revise for coadd split.
c
c	R. Kummerer, July 21, 1988, SER 2064, Smooth capabilities.
c
c	F. Shuman, 1988 Aug  8.  SER 2064, Smooth capabilities.  Conversion
c	          to USEROPEN = CT_CONNECT_READ and CT_CONNECT_QUERY_NONRAW
c	          cause ctuser.inc is no longer maintained.
c
c	R. Kummerer, August 29, 1988, SER 2436, FUT_SETXAX error status.
c
c	Shirley M. Read, August 26, 1988, SPRs 1566 and 1763 requested that
c	          all FIRAS facilities should tell of successful completion
c		  and signal errors via $status for batch processing. Also
c		  added the interface with Fut_Error, the condition handler.
c
c	Shirley M. Read, September 15, 1988, Recoded labels for values which
c		  may be out of bounds, i.e., bad data values. Recoded the
c		  user entry for Boxcar function and window.
c
c	Version 4.1.1 10/25/88, SPR 2479, Shirley M. Read, STX
c		FTB_Ctviewer fails on unrecognized input. The smooth option
c		does not trap unrecognized input but bombs out of the program.
c		For example, the boxcar option, allows only the 'b', not
c		'boxcar' as a character string.
c	Version 4.1.1 10/25/88, SPR 2537, Shirley M. Read, STX
c		Teach FTB_Ctviewer to run in batch. All FIRAS facilities are
c		required to run in batch mode. The interactive/batch option,
c		a dataset qualifier and an Rse file qualifier must be added
c		to the command line. Code switches nust be implemented for
c		the two operational modes.
c	Version 4.1.1 10/25/88, SPR 2582, Shirley M. Read, STX
c		FTB_Ctviewer produces a stack dump when no FIRAS science data
c		records are found for a specified time range. This is a symptom
c		of the lack of proper error checking and condition handling in
c		the program.
c	Version 4.1.2 12/05/88, SPR 2846, Shirley M. Read, STX
c		FTB_Ctviewer does not fill the plot buffer when picking up
c		FNT_Noise spectra records. Zeros appear in the noise spectra.
c		The plot buffer is now updated to fill in the noise spectra.
c	Version 4.1.2 12/07/88, SPR 2963, R. Kummerer, STX
c		Correct smoothing of phase-corrected spectra: smooth real
c		and imaginary parts, instead of amplitude.
c
c	R. Kummerer, January 11, 1989, SPR 2890, FUT_TEMPERATURE GET_CONFIG
c		related changes.
c
c	R. Kummerer, March 24, 1989, SPR 3114, 3398.  Plot model variance,
c		covariance changes.
c
c	R. Kummerer, May 31, 1989, SPR 3931.  Fetch RSE from archive.
c
c	R. Kummerer, July 22, 1989, SER 4176.  Plot new datasets FPP_SDF_xx
c		and FDQ_SDF_xx.
c
c	R. Kummerer, August 18, 1989, SER 3493.  Use new PLT graphics
c		to display IFGs and spectra.
c
c	Q. Chung, STX, August 28, 1989, SER 3306.  Provide version number
c		to track software update.
c
c       D. Bouler , STX, Oct 24, 1989, SPR 4792.  upm_get_value for 'rse'
c               returns an error if only /rse is specified on the command
c               line. Remove code that signals error. Program will get
c               get configured rse if filename for rse is blank.
c
c	R. Kummerer, Dec 27, 1989, Version 5.2, SPR 5523.  Allow fetching
c		short and raw science from separate archives.
c
c	S. Alexander, March 6, 1990, SER's 5726, 5598.  Add capability to
c		pass in a PLT command file.
c
c	R. Kummerer, Mar 6, 1990, Version 5.8, SPR 5161.  Produce plots
c		with useful X-axis units.
c
c	R. Kummerer, Apr 10, 1990, Version 5.10, SPR 5616.  Allow RSE
c		data selection in interactive mode.
c
c       N. Gonzales, Aug 13, 1990, Version 6.7  , SPR 7246. If the
c               CT_KEYEDREAD_ARCV return status is not CTP_NORMAL,
c               signal a fatal error.
c
c	S. Alexander, November 27, 1990. SPR 7732, 7494.  Add OFFSET
c		argument to FUT_SETXAX call.
c
c       F. Shuman, 1991 May 15. SPR 8188. Replaced FUT_TEMPERATURE_LIST
c               with new FUT_TEMP_LIST.
c
c       N. Gonzales, Nov 14, 1991. SPR 9247. Updated the old FEX_GRTSwt file
c               into 3 reference files: FEX_GRTTrans, FEX_GRTRawWt,FEX_GRTCoaWt
c
c	F. Shuman, Hughes STX, 1991 Dec 11.  SPR 9335:  Replace obsolete logical
c	   name 'CSDR$FIRAS_URef' with 'CSDR$FIRAS_Ref'.
c
c       N. Gonzales, Hughes STX, 1992 Jan 23. SPR 9430: Rename FTB_GET_RSE to
c               FUT_GET_RSE.
c
c       N. Gonzales, Hughes STX, 1992 June 30, 1992. SPR 9798: Remove obsolete
c               datasets within the code from FCI, FES, FCS, FFC, FPR, FPS, FSF
c
c       N. Gonzales, Hughes STX, 1992 July 22, 1992. SPRs 9584, 9790: Modified
c               the code to read new FEX_Nyquist reference dataset and to access
c               new FUT_SETXAX subroutine.
c------------------------------------------------------------------------------

	implicit  none

c	Include Files

	include   'ct$library:ctuser.inc'
	include   '(cct_status_record)'
	include   '(cct_get_config)'
	include   '(fut_params)'
	include   '(upm_stat_msg)'
	include   '($ssdef)'

c       Local variables

	integer*4	pstatus		      ! Status of FTB processing
	integer*4	success/1/, error/2/  ! Local status values

	integer   *  4	npts		!number of points in ifg/spec = 512
	integer   *  4	nchan		!number of channels = 4
	integer   *  4	ios		!io status
	integer   *  4	i		! a counter
	integer   *  4	j
	integer   *  4	count
	integer   *  4	ct_nr
	integer   *  4	cv_nr
	integer   *  4	rd_nr
	integer   *  4	rs_nr
	integer   *  2	opens
	integer   *  2  database_id
	integer   *  2  dataset_id
	integer   *  2	rs_name		!dataset name for RSI
	integer   *  4	len
	integer   *  4	key(2)

	integer   *  4	ngroup
	integer   *  4	ngrp
	integer   *  4	sweeps
	integer   *  4	mtm_speed
	integer   *  4	mspd
	integer   *  4	mtm_length
	integer   *  4	chan
	integer   *  4	sci_mode
	integer   *  4	fakeit
	integer   *  4	gain
	integer   *  4  firstsamp /1/
	real      *  4  startx
	real      *  4  spacex

	integer	  *  4  interactive   !operational mode indicator
	character * 32  plot_device
	integer	  *  4  plt_com
	character * 64  plt_com_file
	character * 100 title(3)
	character * 100 ntitle(3)
	character * 100 rtitle(3)

	integer   *  4  pixel_no
	integer   *  4	bol_cmd
	real      *  4	bol_volt
	integer   *  4	resoln
	character * 14	gmt_time
	integer   *  4  time(2)
	character * 60	label
	character *  1	xcal_pos
	character *  3	scan_mode
	integer   *  4  ireftemp
	integer   *  4  pixno
	character * 64  plot_label

	real      *  4	fnyq
	real      *  4	espec(512)
	real      *  4	subdivid(0:1)	!Nyquist freq calc params
	real      *  4	fringes
	real      *  4	opt_path
	real      *  4	temps(10)	!avg instrument temps
	integer   *  4	combine /0/
	logical   *  4	singlifg /.True./
	real      *  4	dumsig(10)

	character * 12	dataset_name  !selected DAFS dataset name for batch
	character * 64  rse_file      !filename for rse file
	integer	  *  4  RSE_Present
	character * 64  filename      !filename string
	integer   *  2  namlen        !length in characters
	character *  3	datatyp

	logical   *  1  proceed       !flag to continue processing current
				      !time range and dataset
	logical   *  4  val_dataset   !dataset chosen in QryBldr is valid
	logical   *  4  eof           !end-of-file flag
	logical   *  4  more
	real      *  4  tmp(512)
	real      *  4  rtmp(512)
	real      *  4  itmp(512)
	real      *  4  tmpn(512)
	complex   *  8  ctmp(1024)
	complex   *  8  ctmp2(1024)
	complex   *  8  ctmpn(1024)

	character *  4	ans
	character *  4	nplt
	character *  4	rplt
	character *  8	inp
	character * 14	str
	character * 14  key_time

	integer   *  4  neg1 / 'FFFFFFFF'X /
	character * 14	starting_time
	character * 14	ending_time
	character * 30	time_range
	logical   *  1	prompt_times
	integer   *  4	lstart
	integer   *  4	lstop
	character *128	rse(16)
	logical   *1	oldrse/.false./

	integer   *4	sfcn
	integer   *4	window
	character *8	smth_fcn
	integer   *4	smth_wnd
	integer   *4	test_wnd
	character *8	smooth_types(1) / 'BOXCAR  ' /
	integer   *  4  ibeg, iend, pos1, pos2, numlen
	integer   *  4  mask / '00000010'x /
	character * 10  cnum
	character * 10  blank /'          '/
	character * 80  inline

c	GET_CONFIG variables.

	dictionary 'fex_nyquist'

	structure /CONFIG_DATA/
	   record /fex_nyquist/ fex_nyquist
	endstructure

	record /CONFIG_DATA/ CONFIG
	record /config_status/ stat(1)

	character *  1	access_mode/' '/	! data set access mode
	integer   *  4  number/1/		! number of data sets
	integer   *  4	size(1)/128/    	! size of data sets

	character * 32  name(1)			! names of data sets
	data name(1)/'csdr$firas_ref:fex_nyquist'/

	integer   *  4	lun(1)			! logical unit numbers
	integer   *  4	cindex(1)		! initial cache pointers
	logical   *  1	new_segment(1)		! flag for new segments
	integer   *  4	ncache/1/
	integer   *  4	ref_count

	integer   *  4  gmt_jstart(2)
	integer   *  4  gmt_jstop(2)

	integer   *  4  fut_query_rse
	integer   *  4  fut_get_rse
	integer   *  4	fut_error
	integer   *  4	fut_setxax
	integer   *  4  fut_smooth_data
	integer   *  4	fut_temp_list
	integer   *  4	fut_timerange
	integer   *  4	upm_get_value
	integer   *  4	upm_present
	integer   *  4  cut_register_version
	integer   *  4  cut_display_banner
	integer   *  4	cut_read_dafs_record_igse_key
	integer   *  4	cut_read_dafs_record_csdr_key
	integer   *  4  cct_open_config
	integer   *  4  cct_get_config_tod
	integer   *  4  cct_close_config
	integer   *  4	lib$get_lun
	integer   *  4	str$upcase
	integer   *  4	str$trim
	integer   *  4  lib$skpc
	integer   *  4  lib$locc
	integer   *  4  ots$cvt_ti_l

	integer   *  4	ct_connect_query_nonraw
	integer   *  4	ct_connect_read

	integer   *  4	status
	integer   *  4	io_stat
	integer   *  4	io_normal

	integer   *  4  num_vol/80/
	integer   *  4  lun_out/6/
	character *  6  version
	parameter       (version='9.8')

	data	subdivid /6., 4./	! Subdivide fringe grating
					! for slow, fast MTM speeds
	parameter	(io_normal = 0)
	parameter	(database_id = 1)
	parameter	(fringes = 20.00E-04)	! Grating spacing in cm
	parameter	(opt_path = 3.4641)	! = 4 cos(30 deg)  =  2 sqrt(3)

	external	ct_connect_query_nonraw
	external	ct_connect_read
	external	fut_error

	external	ftb_normal
	external        ftb_aberr
	external	ftb_ctinit
	external	ftb_ctqbld
	external	ftb_inv_choi1
	external	ftb_inv_choi2
	external	ftb_ctopen
	external	ftb_ctqarc
	external	ftb_ctqget
	external        ftb_ctqarceof
	external	ftb_ctclos
	external	ftb_dsnamlen
	external        ftb_missdsname
	external	ftb_opnconfigerr
	external	ftb_getconfigerr
	external	ftb_clsconfigerr
	external	ftb_ctrarc_err

	record /cct_status/ struc

	dictionary 'cdu_dafs'
	record /cdu_dafs/ dafs_rec

	dictionary 'nfs_sdf'
	record /nfs_sdf/ sci_rec
	dictionary 'fec_sscal'
	record /fec_sscal/ ss_rec
	dictionary 'fnt_noise'
	record /fnt_noise/ noise_rec


c
c Set status for FTB processing to success.
c
	pstatus = success
c
c Establish condition handler and display program banner.
c
	call lib$establish ( fut_error )

	ios = cut_register_version(version)
	ios = cut_display_banner(lun_out,num_vol,
	1		'FIRAS Facility FTB_CTViewer')
c
c Get the operational mode. the default is 'interactive'.
c
	interactive = fac_present
	if ((upm_present('INTERACTIVE') .eq. UPM_Pres) .or.
	1   (upm_present('INTERACTIVE') .eq. UPM_Defaulted)) then
	  interactive = fac_present
	elseif (upm_present('INTERACTIVE') .eq. UPM_Negated) then
	  interactive = fac_not_present
	endif
c
c If batch mode, get the DAFS dataset name.
c
	if (upm_present('DATASET')) then
	  ios = upm_get_value('DATASET', dataset_name, namlen)
	  if (namlen .gt. 12) then
	    pstatus = error
	    call lib$signal(ftb_dsnamlen,%val(1),dataset_name(1:12))
	  endif
	else if (upm_present('DATASET') .Ne. UPM_Absent) Then
	  pstatus = error
	  call lib$signal(ftb_missdsname)
	endif
c
c Get the plot type. The default for interactive mode is the terminal.
c The default for batch is the lineprinter.
c
	if (interactive .eq. fac_present) then
	   plot_device = 'PLT_DEVICE'
	else
	   plot_device = 'PLT_HARDCOPY'
	end if

	if (upm_present ( 'PLOTDEVICE' ) .eq. UPM_Pres) then
	   status = upm_get_value ( 'PLOTDEVICE', plot_device, namlen )
	end if
c
c Get the PLT command file if one is specified:
c
	plt_com = fac_not_present

	if (upm_present ( 'PLTFILE' ) .eq. UPM_Pres) then
	   plt_com = fac_present
	   status = UPM_Get_Value ( 'PLTFILE', plt_com_file, namlen )
	end if
c
c Get the JSTART and JSTOP times if running from an automation command file:
c
	starting_time = fac_jstart_default
	ending_time = fac_jstop_default
	prompt_times = .true.

	if (upm_present ( 'JSTART' ) .eq. UPM_Pres) then
	  ios = upm_get_value ( 'JSTART', starting_time, namlen )

	  if (ios .eq. ss$_normal) then
	    lstart = index(starting_time,' ') - 1
	    if (lstart .eq. -1) lstart=14
	    starting_time = starting_time(1:lstart) //
	2                   fac_jstart_default(lstart+1:)
	  end if

	  prompt_times = .false.
	end if

	if (upm_present ( 'JSTOP' ) .eq. UPM_Pres) then
	  ios = upm_get_value ( 'JSTOP', ending_time, namlen )

	  if (ios .eq. ss$_normal) then
	    lstop = index(ending_time,' ') - 1
	    if (lstop .eq. -1)lstop=14
	    ending_time = ending_time(1:lstop) //
	2                 fac_jstop_default(lstop+1:)
	  end if

	  prompt_times = .false.
	end if

	time_range = starting_time//';'//ending_time//';'

c
c Initialize:
c
	if (pstatus .eq. success) then
	 call ct_init(struc)
c
c If batch mode, get the name of the Rse file.
c
	 RSE_Present = fac_not_present

	 if (upm_present('RSE') .eq. UPM_Pres) then
	  RSE_Present = fac_present
	  ios = upm_get_value('RSE', rse_file, namlen)
	  if (rse_file .eq. ' ') then
	    call ct_gmt_to_binary(starting_time,gmt_jstart)
	    ios = fut_query_rse('FEX_FTBRSE', gmt_jstart, gmt_jstart, rse)
	    if (.not. ios) then
	      pstatus = error
	    endif
	  else
	    ios = fut_get_rse (rse_file, rse)
	    if (.not. ios) then
	      pstatus = error
	    endif
	  end if             ! (rse_file .eq. ' ')
	 elseif (upm_present('RSE') .eq. UPM_Negated) then
	  do i = 1,16
	    rse(i) = ' '
	  enddo
	  do i = 1,16
	    call lib$movc3(4,neg1,%ref(rse(i)(125:128)))
	  enddo
	 endif

	 sfcn = 1
	 window = 1
c
c ....If CobeTrieve initializes, proceed ....
c

	 if (struc.cterr .eq. ctp_normal) then
	  inp = 'Y'

	  do while (inp(1:1) .eq. 'Y')

	   proceed = .true.
	   val_dataset = .true.
c
c If Interactive mode, call the CT Query Builder, Else get the dataset_id.
c
	   if ( Interactive .eq. fac_present .And.
	1	RSE_Present .eq. fac_not_present ) then
c
c Call the CobeTrieve Query Builder for data type (retn'd in struc.dataset_id)
c    and time range
c
	     call ct_query_bld(ctu_$firas, time_range, rse, struc, oldrse)

	   else
	     ios = CUT_Read_Dafs_Record_CSDR_Key (dataset_name, dafs_rec)
	     if ( .not. ios ) then
	       pstatus = error
	       call lib$signal(%val(ios))
	     endif
	     struc.dataset_id = dafs_rec.igse_id.igse_dataset_id
	     inp(1:1) = 'N'
	   endif

c
c First parse TIME_RANGE to get the effective timerange for GET_CONFIG.
c
	   lstart = index(time_range(1:),';') - 1
	   lstop  = index(time_range(lstart+2:),';') - 1

	   starting_time = time_range(1:lstart)
	   ending_time = time_range(lstart+2:lstart+lstop+1)

	   starting_time = starting_time(1:lstart) //
	2                   fac_jstart_default(lstart+1:)
	   ending_time = ending_time(1:lstop) //
	2                 fac_jstop_default(lstop+1:)

	   call ct_gmt_to_binary(starting_time,gmt_jstart)
	   call ct_gmt_to_binary(ending_time,gmt_jstop)

c
c Get GRTs using configuration files (FEX_GRTTRANS, FEX_GRTRAWWT)
c
	   ios = cct_open_config ( gmt_jstart, gmt_jstop,
	1			        number, name, size, access_mode,
	2			        ncache, lun, cindex, stat, ref_count )
	   if (.not. ios) then
	      call lib$signal(ftb_opnconfigerr,%val(1),%val(ios))
	      pstatus = error
	   endif

	   if (((struc.cterr .eq. ctp_normal) .or.
	1       (interactive .eq. fac_not_present)) .and.
	2      (pstatus .eq.success) .and.
	3      (ios)) then

c
c	Encode the 'data type' ...
c
	    datatyp = '   '

	    if ((struc.dataset_id .ge. ctu_$fir_rs1 .and.
	2        struc.dataset_id .le. ctu_$fir_rs4)) then
	         datatyp = 'RSI'
	         plot_label = 'Raw Science Interferogram (NFS)'

	    else if ((struc.dataset_id .ge. fac_fpp_rh .and.
	2             struc.dataset_id .le. fac_fpp_ll)) then
	      datatyp = 'RSP'
	      plot_label = 'Raw Science Interferogram (FPP)'

	    else if ((struc.dataset_id .ge. fac_fdq_rh .and.
	2             struc.dataset_id .le. fac_fdq_ll)) then
	      datatyp = 'RS0'
	      plot_label = 'Raw Science Interferogram (FDQ)'

	    else if (struc.dataset_id .eq. ctu_$fir_ss1) then
	      datatyp = 'SS0'
	      rs_name = ctu_$fir_rs1
	      plot_label = 'Raw Science Interferogram (FEC_SSCAL)'

	    else if (struc.dataset_id .eq. ctu_$fir_ss2) then
	      datatyp = 'SS0'
	      rs_name = ctu_$fir_rs2
	      plot_label = 'Raw Science Interferogram (FEC_SSCAL)'

	    else if (struc.dataset_id .eq. ctu_$fir_ss3) then
	      datatyp = 'SS0'
	      rs_name = ctu_$fir_rs3
	      plot_label = 'Raw Science Interferogram (FEC_SSCAL)'

	    else if (struc.dataset_id .eq. ctu_$fir_ss4) then
	      datatyp = 'SS0'
	      rs_name = ctu_$fir_rs4
	      plot_label = 'Raw Science Interferogram (FEC_SSCAL)'

	    else if (struc.dataset_id .ge. ctu_$fir_nrh .and.
	2            struc.dataset_id .le. ctu_$fir_nll) then
	      datatyp = 'NRH'
	      plot_label = 'Noise Spectrum'

	    else
	      val_dataset = .false.
	    end if
c
c	Open the archive(s) (after converting the numeric dataset_id to the
c	12 (or less) char string name dafs_rec.dataset in interactive mode).
c
	    if ( interactive .eq. fac_present ) then

	      dataset_id = struc.dataset_id
	      status = cut_read_dafs_record_igse_key(database_id, dataset_id,
	2                                          dafs_rec)
	      if (.not. status) proceed = .false.
	    endif

	    status = STR$Trim(dafs_rec.dataset, dafs_rec.dataset, len)
	    if (.not. status) proceed = .false.
	    filename = 'CSDR$FIRAS_ARCHIVE:' // dafs_rec.dataset(1:len)
	2              // '/' // time_range

	    status = LIB$Get_LUN(ct_nr)
	    if (.not. status ) proceed = .false.
	    if ( proceed) then
	       open (unit = ct_nr,
	2            file = filename,
	3            status = 'old',
	4            iostat = io_stat,
	5            useropen = ct_connect_query_nonraw)

	       if (datatyp .eq. 'SS0' .and. io_stat .eq. io_normal) then
	          opens = 1
	          filename = 'CSDR$FIRAS_RAW:' // 'FDQ_SDF_' //
	2                dafs_rec.dataset(len-1:len) // '/' // time_range
	          status = LIB$Get_LUN(rs_nr)
	          open (unit = rs_nr,
	2               file = filename,
	3               status = 'old',
	4               iostat = io_stat,
	5               useropen = ct_connect_read)
	       end if
c
c	Failed to open archive.
c	If the io_stat is not normal, set the proceed flag to false. If the
c       mode is batch, set the processing status to error.
c
	       if (io_stat .ne. io_normal) then
	          call lib$signal(ftb_ctopen,%val(1),%val(io_stat))
	          proceed = .false.
	          if (interactive .eq. fac_not_present) then
	             pstatus = error
	          end if
	       end if

	    end if    !proceed
c
c	Form the collection:
c
	    if ((pstatus .eq. success) .and. ((struc.cterr .eq. ctp_normal)
	1    .or. (interactive .eq. fac_not_present)) .and. proceed) then

	     call ct_query_arcv(, ct_nr, rse, struc)

	     if (struc.cterr .eq. ctp_normal) then

	      more = .true.
c
c	.... Read a record from the proper archive:
c
	      if (datatyp .eq. 'RS0' .or.
	2        datatyp .eq. 'RSI' .or.
	3        datatyp .eq. 'RSP') then
	         call ct_query_get(,ct_nr,sci_rec,struc)
	      else if (datatyp .eq. 'SS0') then
	         call ct_query_get(,ct_nr,ss_rec,struc)
	         if (struc.cterr .eq. ctp_normal) then
	           key(1) = ss_rec.time(1)
	           key(2) = ss_rec.time(2)
	           call ct_keyedread_arcv(,rs_nr,sci_rec,key(1),struc)
	         endif
		 if (struc.cterr .ne. ctp_normal) then
		     call ct_binary_to_gmt(ss_rec.time, key_time)
		     call lib$signal (ftb_ctrarc_err, %val(1), key_time)
	         endif
	      else if (datatyp .eq. 'NRH') then
	         call ct_query_get(,ct_nr,noise_rec,struc)
	      end if

	      if (struc.cterr .ne. ctp_normal) then
	         if (struc.cterr .ne. ctp_endoffile) then
	            call lib$signal(ftb_ctqget,%val(1),%val(struc.cterr))
	            if (interactive .eq. fac_not_present) pstatus = error
		    proceed = .false.
	         end if
	         eof = .true.
	      else
	         eof = .false.
	      end if

	      if (interactive .eq. fac_not_present) ans(1:1) = 'P'

	      do while ( val_dataset .and. (.NOT. eof) )


	         if (datatyp .eq. 'RS0' .or.
	2            datatyp .eq. 'RSI' .or.
	3            datatyp .eq. 'RSP' .or.
	5            datatyp .eq. 'SS0') then

	            chan = sci_rec.sci_head.chan_id
	            ngroup = sci_rec.sci_head.sc_head9
	            sweeps = sci_rec.sci_head.sc_head11
	            mtm_speed = sci_rec.sci_head.mtm_speed
	            mtm_length = sci_rec.sci_head.mtm_length
	            gmt_time = sci_rec.ct_head.gmt
	            gain = fac_gains(sci_rec.sci_head.gain)
	            sci_mode = sci_rec.sci_head.sc_head1a
	            fakeit = sci_rec.dq_data.fake
		    pixel_no = sci_rec.attitude.pixel_no

	            if      (sci_rec.dq_data.xcal_pos .eq. fac_xcaltrans) then
	               xcal_pos = 'T'
	            else if (sci_rec.dq_data.xcal_pos .eq. fac_xcalout)   then
	               xcal_pos = 'O'
	            else if (sci_rec.dq_data.xcal_pos .eq. fac_xcalin)    then
	               xcal_pos = 'I'
	            else
	               xcal_pos = '?'
	            end if

	            if (mtm_speed .eq. 0) then
	               scan_mode(1:1) = 'S'
	            elseif (mtm_speed .eq. 1) then
	               scan_mode(1:1) = 'F'
	            else
		       scan_mode(1:1) = '?'
	            end if

	            if ((sci_rec.sci_head.sc_head9 .lt. 0) .or.
	1	        (sci_rec.sci_head.sc_head9 .gt. 99)) then

	              scan_mode(2:2) = '?'
		      scan_mode(3:3) = '?'

	            else

	              encode (2,60,scan_mode(2:3)) sci_rec.sci_head.sc_head9

	            endif

		    ireftemp = nint ( sci_rec.dq_data.iref_temp * 1000.0 )
		    if ((ireftemp .lt. 0) .or. (ireftemp .gt. 99999))
	1	       ireftemp = 99999

	            pixno = sci_rec.attitude.pixel_no
		    if ((pixno .lt. 0) .or. (pixno .gt. 9999))
	1	       pixno = 9999

	            encode (60,50,label)
	2              sci_rec.ct_head.gmt(1:11),
	3              scan_mode,
	4              ireftemp,
	5              xcal_pos,
	6              pixno

		    if (datatyp .eq. 'RS0') then
	               call ct_gmt_to_binary (gmt_time, time)
c
c Get data from FEX_Nyquist reference dataset using GET_CONFIG.
c
	               ios = cct_get_config_tod (time,number,size,lun,
	1                                     cindex, config, new_segment,stat)
	               if(.not. ios) then
	                 call lib$signal (ftb_opnconfigerr,%val(1),%val(ios))
	               end if

		       ios = fut_setxax(0,fakeit,mtm_speed,mtm_length,chan,
	1                             ngroup,sci_mode,firstsamp,config,
	2                             startx,spacex)
		    else
		       startx = -1.
		       spacex = 1.
		    end if

	         else if (datatyp .eq. 'NRH') then

	            chan = noise_rec.chan.chan_id
	            ngroup = noise_rec.chan.ngroup
	            sweeps = noise_rec.chan.sweeps
	            mtm_speed = noise_rec.chan.mtm_speed
	            mtm_length = noise_rec.chan.mtm_length
	            call ct_binary_to_gmt(
	2                    noise_rec.chan.av_collect_time,
	3                    gmt_time)
	            gain = noise_rec.chan.gain
	            sci_mode = noise_rec.chan.sci_mode
	            fakeit = noise_rec.chan.fakeit
		    pixel_no = noise_rec.attit.pixel_no
	            bol_cmd = noise_rec.chan.bol_cmd_bias

	            if (bol_cmd .lt. 0) then
	               bol_cmd = 256 + bol_cmd
	            end if

	            bol_volt = noise_rec.chan.bol_volt
	            label = noise_rec.coa_head.label

		    ios = fut_setxax(-1,fakeit,mtm_speed,mtm_length,chan,
	1                           ngroup,sci_mode,firstsamp,config,
	2                           startx,spacex)
	        end if

		if(interactive .eq. fac_present) then
c
c	.... Display the data "header":
c
	         type *
	         type *, 'Time tag is ',gmt_time
                 type *, 'Plot for channel ', fac_channel_ids(chan)

	         if (mtm_speed .eq. 0) then
	            type *, 'MTM speed is SLOW'
	         else if (mtm_speed .eq. 1) then
	            type *, 'MTM speed is FAST'
	         else
	            type *, 'MTM speed is UNKNOWN'
	         end if

                 if (mtm_length .eq. 0) then
                    type *, 'MTM length is SHORT'
                 else if (mtm_length .eq. 1) then
                    type *, 'MTM length is LONG'
                 else
                    type *, 'MTM length is UNKNOWN'
                 end if

                 if (gain .ne. -1) then
	            type *, 'Gain setting', gain
	         else
	             type *, 'Gain setting UNKNOWN'
	         end if

                 if (ngroup .ge. 0 .and. ngroup .le. 12) then
	            type *, 'Adds per group is ', ngroup
	         else
	            type *, 'Adds per group is UNKNOWN'
	         end if

                 if (sweeps .ge. 0 .and. ngroup .le. 99) then
	            type *, 'MTM sweeps is ', sweeps
	         else
                    type *, 'MTM sweeps is UNKNOWN.'
                 end if

	         if (sci_mode .ge. 0 .and. sci_mode .le. 4) then
	            type *, 'Science mode is', sci_mode
	         else
	            type *, 'Science mode is UNKNOWN'
	         end if

	         if (fakeit .eq. 1) then
	            type *, 'Fake-it mode is ON'
	         else
	            type *, 'Fake-it mode is OFF'
                 end if

                 type *, 'Pixel number is', pixel_no

                 if (datatyp .ne. 'RS0' .and.
	2               datatyp .eq. 'RSI' .or.
	3               datatyp .eq. 'RSP' .or.
	5               datatyp .ne. 'SS0') then
                     type *, 'Commanded bias is', bol_cmd
                     type *, 'Bolometer voltage is',bol_volt
                 end if
	        endif		!Interactive eq fac_present

c
c	Build the plot title from the information just extracted.
c
		 call fut_plot_title ( label, ngroup, sweeps, mtm_speed,
	2			       mtm_length, chan, sci_mode,
	3			       xcal_pos, fakeit, gain,
	4			       plot_label, title )

		 do i=1,3
		    ntitle(i) = title(i)
		 end do
		 ntitle(1) = ntitle(1)(1:index(ntitle(1),'\')-1) // ' Noise\'


c
c       Display the plot until further notice.
c
	         do while (more)

	          if (interactive .eq. fac_present) then
c
c	.... Prompt for an instruction:
c
	            type 40, smooth_types(sfcn), window
	            accept 12, ans
	            ios = str$upcase(ans,ans)
	           endif
c
c	.... Plot the data:
c
	            if (ans(1:1) .eq. 'P') then
	               more = .false.
	               count = 0
	               if (interactive .eq. fac_not_present) ans(1:1) = ' '
	               if (datatyp .eq. 'RS0' .or.
	2                  datatyp .eq. 'RSI' .or.
	3                  datatyp .eq. 'RSP' .or.
	5                  datatyp .eq. 'SS0') then
	                   do i=1,512
	                      tmp(i) = sci_rec.ifg_data.ifg(i)
	                   end do
		       else if (datatyp .eq. 'NRH') then
			    do i=1,512
			       tmp(i) = noise_rec.chan.spec(i)
			    end do
	               end if

c
c		Display the plot:
c
	                  ios = fut_smooth_data(tmp,512,tmp,sfcn,
	2                                               window)
	                  do i=1,512
	                     ctmp(i) = tmp(i)
	                  end do

	                  if      (datatyp .eq. 'RSI' .or.
	2                          datatyp .eq. 'RSP') then
	                     call fut_plot ( ctmp, 512, startx, spacex, title,
	2				     'Sample Number', 'Counts',
	3				     fac_not_present, interactive,
	4				     plot_device, plt_com,
	5				     plt_com_file )
	                  else if (datatyp .eq. 'RS0' .or.
	3                          datatyp .eq. 'SS0') then
	                     call fut_plot ( ctmp, 512, startx, spacex, title,
	2				     'Distance (cm)', 'Counts',
	3				     fac_not_present, interactive,
	4				     plot_device, plt_com,
	5				     plt_com_file )
	                  else if (datatyp .eq. 'NRH') then
	                     call fut_plot ( ctmp, 257, startx, spacex, title,
	2                                    'Frequency (Hz)', 'Volts/Hz**0.5',
	3				     fac_present, interactive,
	4				     plot_device, plt_com,
	5				     plt_com_file )
	                  end if

		       if(interactive .eq. fac_not_present) ans(1:1) =' '
c
c	    .... End  Display the plot
c
c	.... Jump plots:
c
		    else if (ans(1:1) .eq. 'J') then	            
	               type 100
	               accept *,count
	               more = .false.
	               eof = .false.
c
c	.... Smooth plots:
c
		    else if (ans(1:1) .eq. 'S') then
	               type 200
	               accept 12,inline
	               pos1 = lib$skpc ( ' ', inline )
	               ibeg = pos1
		       pos2 = lib$locc ( ' ', inline(pos1:80))
		       iend = ibeg + pos2 - 2
	               smth_fcn = inline(ibeg:iend)
		       ibeg = iend + 1
		       pos1 = lib$skpc ( ' ', inline(ibeg:80))
	               if ( pos1 .eq. 0 ) then
			 type *, ' No window specified.'
		       else	! Window specified
		         ibeg = ibeg + pos1 - 1
			 pos2 = lib$locc ( ' ', inline(ibeg:80))
		         if (pos2 .eq. 0) then
			   iend = 80
			   numlen = 81 - ibeg
		         else
		           iend = ibeg + pos2 - 2
			   numlen = iend - ibeg + 1
			 endif
		         cnum(1:10) = blank(1:10)
			 pos1 = 11 - numlen
			 cnum(pos1:10) = inline(ibeg:iend)
			 status = ots$cvt_TI_L ( cnum, smth_wnd,
	1			  %val(2), mask )
		         If ( status .ne. ss$_normal ) then
			   type *, ' Invalid window specified.'
		         else
	                   more = .false.
	                   eof = .false.
	                   ios = str$upcase(smth_fcn,smth_fcn)
	                   if (smth_fcn(1:1) .eq. 'B') then
	                     sfcn = fac_boxcar
	                     test_wnd = smth_wnd - 1
	                     if (mod(test_wnd,2) .eq. 0) then
	                       window = smth_wnd
	                     else
	                       type *, ' Illegal window size.'
	                     end if
	                   else
	                     type *, ' Unknown smooth function.'
	                   end if
			 endif  ! Invalid window specified
		       endif	! Window specified or not
c
c	.... Quit CTVIEWER:
c
	            else if (ans(1:1) .eq. 'Q') then	            
	               more = .false.
	               eof = .true.
c
c	.... Hit return:
c
		    else if (ans(1:1) .eq. ' ') then
	               more = .false.
	               count = 1
	               eof = .false.
	               if(interactive .eq. fac_not_present) ans(1:1) = 'P'
c
c	.... Avoid skipping to the next record
c
		    else
	               more = .false.
	               count = 0
	               eof = .false.
	            end if              !if (ans(1:1) .eq. 'P')
c
c	  .... :End Plot the data
c
	         end do                 !do while (more)

c	  .... :End Prompt for an instruction

c
c	Get next data record.
c
	         if (.not. eof) then
	            do while (count .ge. 1 .and. struc.cterr .eq. ctp_normal)
	               if (datatyp .eq. 'RS0' .or.
	2                  datatyp .eq. 'RSI' .or.
	3                  datatyp .eq. 'RSP') then
	                   call ct_query_get(,ct_nr,sci_rec,struc)
	               else if (datatyp .eq. 'SS0') then
	                  call ct_query_get(,ct_nr,ss_rec,struc)
	                  if (struc.cterr .eq. ctp_normal) then
	                    key(1) = ss_rec.time(1)
	                    key(2) = ss_rec.time(2)
	                    call ct_keyedread_arcv(,rs_nr,sci_rec,key(1),struc)
	                    if (struc.cterr .ne. ctp_normal) then
		                call ct_binary_to_gmt(ss_rec.time, key_time)
				call lib$signal (ftb_ctrarc_err, %val(1),
	1                                        key_time)
			    endif
			  end if
	               else if (datatyp .eq. 'NRH') then
	                  call ct_query_get(,ct_nr,noise_rec,struc)
	               end if
	               count = count - 1
	            end do
	            if (struc.cterr .eq. ctp_normal) then
	               more = .true.
	               eof = .false.
	            else
	               if (struc.cterr .ne. ctp_endoffile) then
	                  call lib$signal(ftb_ctqget,%val(1),%val(struc.cterr))
	               end if
	               more = .false.
	               eof = .true.
	            end if
	         end if                 !if (.not. eof)
c
c	  .... :End Get next data record
c
	      end do                    !do while (val_dataset .and. .not. eof)
c
c   .... Failed to form collection:
c
	     else
	      if (struc.cterr .eq. ctp_endoffile) then
	        call lib$signal(ftb_ctqarceof, %val(1), %val(struc.cterr))
	      else
	        call lib$signal(ftb_ctqarc, %val(1), %val(struc.cterr))
	        if (interactive .eq. fac_not_present) then
	          pstatus = error
	        endif
	      endif
	     end if
c
c
c   .... End of file; close the archive(s):
c
	     call ct_close_arcv(, ct_nr, struc)
c
c   .... (Note: This version stops trying after the 1st arcv close failure.)
c
	     if (struc.cterr .eq. ctp_normal) then
	      call lib$erase_page(1,1)
	     else
	      call lib$signal(ftb_ctclos, %val(1), %val(struc.cterr) )
	      pstatus = error
	     end if
c
c   .... Failed to open archive or failed to setup for open correctly.
c
	    else
	     if (interactive .eq. fac_not_present) pstatus = error
	    end if
c
c   .... Failure in CT Query Builder:
c
	   else
	    if (interactive .eq. fac_present) then
	      call lib$signal(ftb_ctqbld, %val(1), %val(struc.cterr) )
	    endif
	   end if

c
c Close configuration file.
c
	   ios = cct_close_config ( number, lun, cindex )
	   if (.not. ios) then
	      call lib$signal(ftb_clsconfigerr,%val(1),%val(ios))
	   endif

	   if ( .not. val_dataset ) then
	    call lib$signal( ftb_inv_choi1, %val(0), ftb_inv_choi2, %val(0) )
	   end if
	   if (interactive .eq. fac_present) then
	     type 5
	     accept 12,inp
	     ios = str$upcase(inp,inp)
	     type *
	   endif
	  end do                ! do while (inp(1:1) .eq. 'Y')

	 else
c
c .... if CobeTrieve won't initialize, then, Yer outa there!
c
	  call lib$signal ( ftb_ctinit, %val(1), %val(struc.cterr) )
	  pstatus = error
	 end if
	endif    !pstatus is success after input from command line

c
c	Format statements.
c
5	format (/,' Another query? (Y/[N]) ',/,'          => ',$)
12	format(a)
13	format(' Plot the noise (Y/[N])? ', $)
40	format(' P to see plot',/,
	2      ' S to smooth plots  Function = ',a,i,/,
	3      ' J to jump plots ',/,
	4      ' Q to quit displaying plots',/,
	5      ' Hit return for the next plot',/,
	6      10x,' => ',$)
50	format(a11,'_',a3,'_',i5,'_',a1,'_',i4)
60	format(i2)
100	format(' Enter number of plots to skip => ',$)
200	format(	' BOXCAR <window>',/,
	2       '   Enter Smooth Function => ',$)

c
c       Exit the program.
c
	if ( pstatus .eq. success ) THEN
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
	endif

	end
