	Integer * 4 Function FFL_Produce_Spectra (coadrec, fflc, ffli, ffls)

c-------------------------------------------------------------------------------
c
c	Function FFL_PRODUCE_SPECTRA
c
c	This function apodizes and rotates the coadded IFG, FFTs the IFG into
c	a spectrum, and converts the spectrum from counts to volts.  It then
c	apodizes the variances and produces a voltage sigma vector.
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  19 February 1993
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c	    coadrec record	coadd structure
c	    fflc record		configdata structure (defined in ffl_config.txt)
c	    ffli record		invoc structure (defined in ffl_invoc.txt)
c	    ffls record		fishspec structure (defined in ffl_spec.txt)
c
c	Output:
c	    fflc record		configdata structure (defined in ffl_config.txt)
c	    ffli record		invoc structure (defined in ffl_invoc.txt)
c	    ffls record		fishspec structure (defined in ffl_spec.txt)
c
c	Subroutines called:
c	    FUT_Get_Recnum
c	    FUT_Temp_List
c	    FUT_Apod_RecNumL
c	    FUT_Apod_RotL
c	    DFFTRF
c	    LIB$MovC5
c	    LIB$Signal
c
c	Include files:
c	    fut_params.txt
c	    ffl_config.txt
c	    ffl_invoc.txt
c	    ffl_spec.txt
c
c-------------------------------------------------------------------------------
c   Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c       Modifications to handle low frequency short fast data for
c       length 680 FFTs.      Alice Trenholme, GSC, 29 November 1994
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c	Record declaration removed from include file ffl_spec.txt and placed in
c	   here; use vrecsize from there in establishing memory space for fish
c	   record.  Fred Shuman, HSTX, 1995 June 30.
c
c	Routine FFL_Read_Reference deleted; only job remaining for it to do was
c	   to read the ETF.  Fred Shuman, HSTX, 1995 July 13.
c
c	Use vspecsiz to generate various array sizes and in LIB$MovC5 calls.
c	Fred Shuman, HSTX, 1995 July 21.
c-------------------------------------------------------------------------------

	Implicit None

	Include '(fut_params)'
	Include '(ffl_invoc)'		! defines record ffli (structure invoc)
c
c  Call arguments:
c
	Dictionary 'fil_sky'
	Record /fil_sky/ coadrec
	Include '(ffl_config)'		! defines record fflc (structure config)
		!  and certain other vbles; Includes CCT_Get_Config; Dictionary
		!  structures FEX_GRTCOAWT, FEX_GRTTRANS.
		!defines structure fishspec and parameters vrecsize and vspecsiz
	Record /fishspec/ ffls
c
c  All other variables and functions:
c
	Complex * 16	vspec(vspecsiz+1)!  double precision voltage spectrum

	Integer *  2	cmd_bias	!  integer valued commanded bias

	Integer *  4	j, k		!  counters
	Integer *  4	last_upmode /-1/!  microprocessor change flag
	Integer *  4	status		!  return status
	Integer *  4	speed		!  MTM speed
	Integer *  4	offgrts		!  return value from FUT_Temp_List
	Integer *  4	io_stat		!  I/O return status
	Integer *  4	upmode		!  microprocessor mode
	Integer *  4	fakeit		!  fakeit flag
	Integer *  4	erecno		!  etf record number
	Integer *  4	arecno		!  apodization function record number
	Integer *  4	lnrized /0/	!  'linearized' flag for call to
					!		 FUT_Apod_RecNumL
	Integer *  4	peakpos		!  peak position read from FIL record
	Integer *  4	fftlen		!  FFT length (from coadd record)
	Integer *  4	speclen		!  spectrum length computed from fftlen

	Real	*  4	temp(10)	!  instrument temperatures
	Real	*  4	tsigma(10)	!  temperature sigmas

	Real	*  4	aifg(512)	!  single precision ifg
	Real	*  8	daifg(512)	!  double precision ifg
	Real	*  8	apod(512)	!  apodization function
	Real	*  8	rifg(2*vspecsiz)!  apodized, rotated IFG
	Real	*  8	vsigma(vspecsiz+1)!  voltage sigmas
	Real	*  8	zifg(2*vspecsiz)!  double precision transformed ifg
	Complex	* 16	phase_shift	!  double prec. phase shift for HSF

	Save last_upmode

	Integer *  4	FFL_Read_Reference
	Integer *  4	FUT_Apod_RecNumL
	Integer *  4	FUT_Apod_RotL
	Integer *  4	FUT_Temp_List
	Real	*  8	UOE_ADT2T68

	External	ffl_etfread
	External	ffl_ftemplistgrt
	External	ffl_ftemplistcall
	External	ffl_fapodrecnuml
	External	ffl_apodread
	External	ffl_fapodrotl
	External	ffl_normal
	External	fut_normal

CC  Initialize the routine.
CC
	status = %Loc(ffl_normal)
c
c  Initialize the spectrum record.
c
	Call LIB$MovC5 (0, , 0, 4*vrecsize, ffls)
	fftlen = coadrec.coad_data.fft_length
	speclen = fftlen/2 + 1
	ffls.nifgs = coadrec.coad_spec_head.num_ifgs
	ffls.adj_nifgs = coadrec.coad_spec_head.adj_num_ifgs
	ffls.time=(uoe_adt2t68(coadrec.ct_head.time)-fflc.apco_eject_time) /
     &                 year_len
	ffls.volt = Dble(coadrec.coad_spec_data.bol_volt)
	cmd_bias  = coadrec.coad_spec_data.bol_cmd_bias		!  Two's
	if (cmd_bias .Lt. 0) cmd_bias = cmd_bias + 256		!  complement
	ffls.bias = Dble(cmd_bias)/25.5D0
	ffls.gain = coadrec.coad_spec_data.gain_sum /
     &			coadrec.coad_spec_head.num_ifgs
	upmode = coadrec.coad_spec_data.sci_mode
	fakeit = coadrec.coad_spec_data.fakeit
	speed = coadrec.coad_spec_data.mtm_speed
	peakpos = coadrec.coad_spec_data.peak_pos
c
c  Extract the IFG from the coadd record--for smode = [1..4] extract the full
c  512-point ifg; for smode = 5 (6) [xLSF (xLLF) data] the IFG has already been
c  decimated (truncated) by FIL, so extract the first 128 points of the coadd.
c
	Call LIB$MovC5 (0, , 0, 4096, daifg)
	Do j=1,ifglen(ffli.smode)
	   daifg(j) = Dble(coadrec.coad_data.ifg(j))
	EndDo
c
c  Since the variances are already apodized and have had their units fixed by
c  FIL (the new version of FIC), they may be converted to sigmas here without
c  fussing around.  Extract the variances from the coadd record; normalize by
c  Nyquist_icm*etendu.
c
	Call LIB$MovC5 (0, , 0, 8*(vspecsiz+1), vsigma)

	Do j = 1, speclen
	   vsigma(j) = DSqrt(Dble(coadrec.coad_data.real_var(j))) *
     &                  coadrec.coad_spec_data.Nyquist_icm * fac_etendu
	EndDo
c
c  Read the appropriate electronics transfer fn for a new microproc. mode --
c    note that for getting the ETF record number for scan modes FS & FL, adds/gp
c    of 2 (ffli.ngpetf), not 8 (ffli.ngroup) must be used; read the ETF.
c
	ffli.ngroup = coadrec.coad_spec_data.adds_per_group
	ffli.ngpetf = ffli.ngroup
	if (ffli.smode .Ge. 5) ffli.ngpetf = 2
	If (last_upmode .Ne. upmode) Then
	   last_upmode = upmode
	   Call FUT_Get_Recnum (fakeit, speed, ffli.chan, upmode,
     &                          ffli.ngpetf, erecno)
	   Read (fflc.dir_lun(1)''erecno, iostat=io_stat) fflc.ztrans

	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_etfread)
	      Call LIB$Signal(ffl_etfread,%Val(1),%Val(erecno),%Val(io_stat))
	   EndIf
	EndIf

	If (status .Ne. %Loc(ffl_normal)) Then
	   FFL_Produce_Spectra = status
	   Return
	EndIf
c
c  Combine the instrument temperatures.  If returned value of FUT_Temp_List is
c    .Eq. 0, it represents a normal return status.  If .Le. 10,  it represents
c    the lowest index of 'temp' for which all GRT switches were off.  It can be
c    .Gt. 10 only if FUT_Temp_List is called improperly--a matter for the
c    programmer, not the user, to fix.
c
	offgrts = FUT_Temp_List (coadrec.en_analog, coadrec.en_sigma,
     &				 fflc.config.grtcoawt, fflc.config.grttrans,
     &				 0, .False., temp, tsigma)
	If (offgrts .Gt. 0 .And. offgrts .Le. 10) Then
	   FFL_Produce_Spectra = %Loc(ffl_ftemplistgrt)
	   Call LIB$Signal (ffl_ftemplistgrt, %Val(2), %Val(offgrts),
     &			    coadrec.ct_head.gmt)
	   Return
	ElseIf (offgrts .Gt. 10) Then
	   FFL_Produce_Spectra = %Loc(ffl_ftemplistcall)
	   Call LIB$Signal (ffl_ftemplistcall, %Val(1), %Val(offgrts))
	   Return
	Else
	   Do j=1,6
	      ffls.temp(j) = temp(j)
	      ffls.tsigma(j) = tsigma(j)
	   EndDo
	   ffls.temp(7) = temp(ffli.chan+6)
	   ffls.tsigma(7) = tsigma(ffli.chan+6)
	EndIf
CC
CC  Generate the voltage spectrum.
CC
c  Apodize and rotate the IFG in preparation for the FFT.
c
	Do j=1,512
	   aifg(j) = Sngl(daifg(j))
	EndDo
	status = FUT_Apod_RecNumL (ffli.chan, ffli.smode, fakeit, upmode,
     &                             ffli.ngroup, lnrized, arecno)
	If (status .Ne. %Loc(fut_normal)) Then
	   FFL_Produce_Spectra = %Loc(ffl_fapodrecnuml)
	   Call LIB$Signal (ffl_fapodrecnuml, %Val(1), %Val(status))
	   Return
	Else
	   Read (fflc.dir_lun(2)''arecno, iostat=io_stat) apod
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_apodread)
	      Call LIB$Signal(ffl_apodread,%Val(2),%Val(arecno),%Val(io_stat))
	   EndIf

	   status = FUT_Apod_RotL (aifg, apod, fftlen, peakpos, rifg)
	   If (status .Ne. %Loc(fut_normal)) Then
	      FFL_Produce_Spectra = %Loc(ffl_fapodrotl)
	      Call LIB$Signal (ffl_fapodrotl, %Val(1), %Val(status))
	      Return
	   Else
c
c  Do the fft, repacking the complex spectrum from the real transformed array,
c    zifg, into the complex array, vspec.
c
	      status = %Loc(ffl_normal)
	      Call LIB$MovC5 (0, , 0, 16*(vspecsiz), zifg)
	      Call LIB$MovC5 (0, , 0, 16*(vspecsiz), vspec)

	      Call DFFTRF (fftlen, rifg, zifg)
	      vspec(1) = DCmplx(zifg(1), 0.0D0)
	      Do j=2,speclen-1
	         vspec(j) = DCmplx(zifg(2*j-2), zifg(2*j-1))
	      EndDo
	      vspec(speclen) = DCmplx(zifg(fftlen), 0.0D0)
c
c  Normalize the spectrum and convert it from counts to volts.
c    Fix the short fft for the low frequency fast spectra.
c    Shift the phase of the high frequency short fast spectra.
c
	      if ((ffli.chan .Eq. 1  .Or.  ffli.chan .Eq. 3)  .And.
     &                                                     ffli.smode .Eq. 2)
     &              phase_shift = DCmplx(0.0D0, fac_dpi/fftlen)

	      Do j = 2,speclen
	         vspec(j) = vspec(j) / fac_adc_scale / fflc.ztrans(j)
     &				    * CDExp(Dble(j-1)*phase_shift)
	      EndDo
CC
CC  Fill in the voltage spectrum record.
CC
	      Do j = 2,speclen
	         ffls.vspec(j-1) = vspec(j)
	         ffls.vsigma(j-1) = vsigma(j)
	      EndDo

	   EndIf		!	( status from FUT_Apod_RotL is normal
	EndIf			!	( status from FUT_Apod_RecNumL is normal

	FFL_Produce_Spectra = status
	Return
	End
