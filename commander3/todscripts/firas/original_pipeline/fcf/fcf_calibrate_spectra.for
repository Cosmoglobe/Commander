	Integer * 4 function  fcf_calibrate_spectra (vspec, vvar, cspec,
     .						     cvar, spec_rec)

c-------------------------------------------------------------------------------
c
c	Function FCF_CALIBRATE_SPECTRA
c
c	This function drives the FCF calibration routines by calling the
c	routines FCF_Temporal_Drift, FCF_Calc_Responsivity,
c	FCF_Pcalib_Variances, FCF_Apply_Model, and FCF_Doppler_Shift.
c	These routines produce calibrated spectra and calibrated spectrum
c	variances with the units of ergs/sec/cm**2/sr/icm.  The sky spectra are 
c	in the barycentric frame of reference.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  11 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		spec_rec			calibrated spectrum record
c		vspec		complex*16	voltage spectrum
c		vvar(3)		real   * 8	voltage variances
c
c	Output:
c		spec_rec			calibrated spectrum record
c		cspec		complex*16	calibrated spectrum
c		cvar(3)		real   * 8	calibrated variances
c
c	Subroutines called:
c		fcf_apply_model
c		fcf_calc_responsivity
c		fcf_doppler_shift
c		fcf_pcalib_variances
c		fcf_temporal_drift
c		fut_temp_list
c		lib$signal
c
c	Include files:
c		fcf_config.txt
c		fcf_display.txt
c		fcf_invoc.txt
c		fcf_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants:
c
c		VBLIM		40.0  km/sec		upper limit on projected
c							barycentric velocity
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Rearrange DISPLAY structure for differential and reference spectra.
c	Update all four bolometer temperatures with the operating point
c	temperature.
c	Gene Eplee, GSC, 11 July 1994
c	SER 11826
c
c	Write the updated bolometer temperatures from the internal temperature
c	array into the spectrum record.
c	Gene Eplee, GSC, 5 August 1994
c	SPR 11862
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_display)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	complex * 16	cspec(257)		!  calibrated spectrum
	complex * 16	vspec(257)		!  voltage spectrum

	integer * 2	bol_cmd_bias		!  commanded bolometer bias

	integer *  4	j			!  a counter
	integer *  4	k			!  a counter
	integer *  4	nifgs			!  number of ifgs in the
						!    spectrum
	integer *  4	rstatus			!  return status
	integer *  4	status			!  return status

	real	*  4	s_temp(10)		!  instrument temperatures
	real	*  4	s_tsig(10)		!  temperature sigmas

	real	*  8	bol_volt		!  detector voltage
	real	*  8	cmd_bias		!  commanded bolometer bias
	real	*  8	cvar(3,257)		!  calibrated variancess
	real	*  8	phase_corr		!  linear phase corrector
	real	*  8	power			!  IR power incident on detector
						!    from Planck functions
	real	*  8	primary_vib		!  primary (time-dependent)
						!    vibration correction
	real	*  8	Qrad			!  IR power incident on detector
						!    from bolometer model
	real	*  8	Tbol			!  actual detector temperature
	real	*  8	ticald			!  Ical temperature drift
						!    correction
	real	*  8	Tdet			!  measured detector temperature
	real	*  8	temp(10)		!  instrument temperatures
	real	*  8	tsig(10)		!  temperature sigmas
	real	*  8	vbary			!  proj barycentric velocity
	real	*  8	vblim			!  upper limit on vbary
	parameter	(vblim = 40.0D0)
	real	*  8	vvar(3,257)		!  voltage variances

	integer *  4	fcf_apply_model
	integer *  4	fcf_calc_responsivity
	integer *  4	fcf_doppler_shift
	integer *  4	fcf_pcalib_variances
	integer *  4	fcf_temporal_drift
	integer *  4	fut_temp_list

	dictionary 'fcf_sky'
	record /fcf_sky/	spec_rec

	external	fcf_ftemplist
	external	fcf_invalbaryvel
	external	fcf_normal
	external	fut_normal

C
C  Combine the instrument spectra.
C

	if (fcc_single .eq. fac_present) then
c
c  Combine the instrument temperatures for single ifg spectra, setting the
c	combined temperature sigmas to zero.
c
	   rstatus = fut_temp_list (spec_rec.en_analog, spec_rec.en_sigma,
     .				    config.grtcoawt, config.grttrans,
     .				    0, .true., s_temp, s_tsig)
	   if (rstatus .ne. 0) then
	      fcf_calibrate_spectra = %loc(fcf_ftemplist)
	      call lib$signal (fcf_ftemplist, %val(3), %val(status),
     .			 spec_rec.ct_head.gmt, %val(spec_rec.attitude.pixel_no))
	      return
	   else
	      status = %loc(fcf_normal)
	   endif

	   do j = 1,10
	      temp(j) = dble(s_temp(j))
	      tsig(j) = 0.0D0
	      spec_rec.coad_spec_data.temp(j) = s_temp(j)
	      spec_rec.coad_spec_data.temp_sigma(j) = s_tsig(j)
	   enddo

	else
c
c  Combine the instrument temperatures for ensemble ifg spectra.
c
	   rstatus = fut_temp_list (spec_rec.en_analog, spec_rec.en_sigma,
     .				    config.grtcoawt, config.grttrans,
     .				    0, .false., s_temp, s_tsig)
	   if (rstatus .ne. 0) then
	      fcf_calibrate_spectra = %loc(fcf_ftemplist)
	      call lib$signal (fcf_ftemplist, %val(3), %val(status),
     .			 spec_rec.ct_head.gmt, %val(spec_rec.attitude.pixel_no))
	      return
	   else
	      status = %loc(fcf_normal)
	   endif

	   do j = 1,10
	      temp(j) = dble(s_temp(j))
	      tsig(j) = dble(s_tsig(j))
	      spec_rec.coad_spec_data.temp(j) = s_temp(j)
	      spec_rec.coad_spec_data.temp_sigma(j) = s_tsig(j)
	   enddo

	endif	!	(fcc_single


C
C  Calibrate the variances and the spectrum.
C

c
c  Calculate the temporal drift corrections.
c
	status = fcf_temporal_drift (spec_rec.ct_head.time(2), ticald,
     .				     primary_vib)
	temp(2) = temp(2) - ticald

	if (status .eq. %loc(fcf_normal)) then
c
c  Calculate the actual detector temperature and the detector responsivity and
c	time constant.
c
	   bol_volt = dble(spec_rec.coad_spec_data.bol_volt)
	   bol_cmd_bias = spec_rec.coad_spec_data.bol_cmd_bias !Two's complement
	   if (bol_cmd_bias .lt. 0) bol_cmd_bias = bol_cmd_bias + 256
	   cmd_bias = dble(bol_cmd_bias)/25.5D0
	   Tdet     = temp(fcc_chan+6)
	   status = fcf_calc_responsivity (bol_volt, cmd_bias, Tdet, Tbol, Qrad)
	   do j = 1,4
	      temp(6+j) = Tbol
	   enddo

	   if (status .eq. %loc(fcf_normal)) then
c
c  Compute the bolometer delay.
c
	      call lib$movc5 (0,,0,4112,B)
	      do j = 2,257
	         if (cdabs(emiss(1,j)) .gt. trlim) then
	            B(j) = dcmplx(1.0D0, afreq(j)*tau)
	         endif
	      enddo

	      if (fcc_single .eq. fac_present) then
c
c  Set the calibrated variances to flag values for single ifg spectra.
c
	         do j = 1,3
	            do k = 1,257
	               cvar(j,k) = fcc_varflag
	            enddo
	         enddo
	      else
c
c  Apply the calibration model to the variances for coadded ifg spectra.
c
	         status = fcf_pcalib_variances (vvar, cvar)
	      endif

	      if (status .eq. %loc(fcf_normal)) then
c
c  Apply the calibration model to the spectra.
c
	         nifgs = spec_rec.coad_spec_head.num_ifgs
	         status = fcf_apply_model (primary_vib, nifgs, temp, tsig,
     .					  vspec, cvar, power, phase_corr, cspec)
	      endif


	      if ((status .eq. %loc(fcf_normal))  .and.
     .		  (fcc_sky .eq. fac_present)) then
c
c  Correct the sky spectra for the projected spacecraft barycentric velocity
c	along the Firas line of sight, checking for valid projected barycentric
c	velocities.
c
	         vbary =
     .	         dble(spec_rec.attitude.projected_barycentric_velocity)/100.0D0
	         if (dabs(vbary) .gt. vblim) then
	            status = %loc(fcf_invalbaryvel)
	            call lib$signal (fcf_invalbaryvel, %val(2),
     .				     spec_rec.ct_head.gmt,
     .				     %val(spec_rec.attitude.pixel_no))
	         elseif (vbary .ne. 0.0D0) then
	            status = fcf_doppler_shift (vbary, cspec)
	         endif
	      endif

	   endif	!  (status from fcf_calc_responsivity

	endif		!  (status from fcf_temporal_drift


C
C  Update the spectrum record.
C

	if (status .eq. %loc(fcf_normal)) then
c
c  Write the calibrated spectrum and variances into the spectrum record.
c
	   do j = 1,257
	      display.cspec(j) = cmplx(cspec(j))
	      spec_rec.spec_data.spec(j) = display.cspec(j)
	      display.cvar(1,j) = sngl(cvar(1,j))
	      display.cvar(2,j) = sngl(cvar(2,j))
	      display.cvar(3,j) = sngl(cvar(3,j))
	      spec_rec.spec_data.real_var(j)      = display.cvar(1,j)
	      spec_rec.spec_data.imag_var(j)      = display.cvar(2,j)
	      spec_rec.spec_data.real_imag_var(j) = display.cvar(3,j)
	   enddo

c
c  Fill in the model solution fields of the spectrum record.
c  Update the Ical and detector temperature fields in the spectrum record.
c
	   spec_rec.spec_data.responsivity          = sngl(S0)
	   spec_rec.spec_data.time_constant         = sngl(tau)
	   spec_rec.spec_data.phase_corr            = sngl(phase_corr)
	   spec_rec.spec_data.qrad                  = sngl(Qrad)
	   spec_rec.spec_data.ir_power              = sngl(power)
	   spec_rec.coad_spec_data.temp(2)          = sngl(temp(2))
	   do j = 1,4
	      spec_rec.coad_spec_data.temp(6+j)     = sngl(temp(6+j))
	   enddo
	   spec_rec.spec_data.model_ttag            = model_ttag
	   spec_rec.spec_data.model_label           = model_label

	endif	!  (status from fcf_apply_model


	fcf_calibrate_spectra = status

	return
	end
