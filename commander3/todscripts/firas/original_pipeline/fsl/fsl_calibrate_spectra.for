	Integer * 4 function  fsl_calibrate_spectra (vspec, vvar, cspec,
     .						     cvar, spec_rec)

c-------------------------------------------------------------------------------
c
c	Function FSL_CALIBRATE_SPECTRA
c
c	This function drives the FSL calibration routines by calling the
c	routines FSL_Temporal_Drift, FSL_Calc_Responsivity,
c	FSL_Pcalib_Variances, FSL_Apply_Model, and FSL_Doppler_Shift.
c	These routines produce calibrated spectra and calibrated spectrum
c	variances with the units of ergs/sec/cm**2/sr/icm.  The spectra and
c       variances are then converted to units of MJy/sr for FSL output records.
c       The sky spectra are in the barycentric frame of reference.
c
c	Author:   
c                FCF_Calibrate_Spectra
c                Gene Eplee
c	         General Sciences Corp.
c		 11 May 1993
c       
c                FSL_Calibrate_Spectra
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fsl_apply_model
c		fsl_calc_responsivity
c		fsl_doppler_shift
c		fsl_pcalib_variances
c		fsl_temporal_drift
c		fut_temp_list
c		lib$signal
c
c	Include files:
c		fsl_config.txt
c		fsl_display.txt
c		fsl_invoc.txt
c		fsl_model.txt
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
c	Changes for FCF:
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
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 9, 1995 
c       Modified FCF_Calibrate_Spectra to FSL_Calibrate_Spectra for the new 
c       FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.
c           1. Changed status, include file, record, and function names for FSL.
c           2. To accomodate longer spectral arrays, changed array sizes.
c           3. Converted units from ergs/sec/cm**2/sr/icm to MJy/sr.
c           4. Used start and stop frequency indices corresponding to model
c              solution to store calibrated spectrum and variances in output
c              spectrum record and display record.
c           5. Stored new fields in output spectrum record: fft_length,
c              lofreq_bin, and hifreq_bin. Stored number of points to plot
c              in display record.
c           6. Check for flag to set combined temperature sigmas from 
c              fut_temp_list to zero. If the flag is set or if the number of
c              IFGs in the coadd is less than or equal to 2, then set the
c              combined temperature sigmas to zero for internal processing and
c              in the ouput spectrum record.
c           7. If the flag to use FIL averaged variances for the autophase 
c              corrector computation is set, weight the FIL variance vectors
c              by the glitch rate adjusted number of IFGs for the coadd and
c              replace the individual voltage variances by the FIL average 
c              weighted variances. Then call FSL_Pcalib_Variances as before.
c           8. Set the flags in the output sprectrum record accoring to whether
c              the individual coadd variances, the FIL averaged voltage 
c              variances, or the Dvector computed by a previous FSL run was
c              used as a weight for the autopahse corrector computation.
c           9. Separate the storing of the spectrum and variances into the
c              display record and output spectum record. If the Dvector was
c              used for the autophase corrector computation, store the weighted
c              Dvector in the real-real variance vector and set the others to 
c              zero. 
c          10. If the individual coadd variances are used and there are less 
c              than 3 IFGs in the coadd, the coadded IFG cannot be calibrated.
c              Set the calibrated flag to zero in this case and to 1 for
c              more than 2 IFGs or use of the FIL or FSL variances.  
c          11. Since requirements are changed to weight the Dvector be the 
c              glitch rate adjusted number of IFGs, this number must be passed
c              to FSL_Apply _Model.
c
c       S. Brodd, HSTX, 2/11/97, SPR 12342.  Use FIL averaged variances for
c                                            autophase corrector, but put
c                                            individual calibrated variances
c                                            in output record.
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_display)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	complex * 16	cspec(361)		!  calibrated spectrum
	complex * 16	vspec(361)		!  voltage spectrum

	integer * 2	bol_cmd_bias		!  commanded bolometer bias

	integer *  4	j			!  a counter
	integer *  4	k			!  a counter
	integer *  4	nifgs			!  number of ifgs in the
						!    spectrum
	integer *  4	rstatus			!  return status
	integer *  4	status			!  return status

	real	*  4	s_temp(10)		!  instrument temperatures
	real	*  4	s_tsig(10)		!  temperature sigmas
	real    *  4    adjnifgs                !  glitch rate adjusted number
                                                !  of ifgs
	
	real	*  8	bol_volt		!  detector voltage
	real	*  8	cmd_bias		!  commanded bolometer bias
	real	*  8	cvar(3,361)		!  calibrated variances
	real	*  8	cvar_out(3,361)		!  calibrated variances
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
	real	*  8	vvar(3,361)		!  voltage variances
	real    *  8    erg_to_mjy_sq            !  square of erg to mjy factor

	integer *  4	fsl_apply_model
	integer *  4	fsl_calc_responsivity
	integer *  4	fsl_doppler_shift
	integer *  4	fsl_pcalib_variances
	integer *  4	fsl_temporal_drift
	integer *  4	fut_temp_list

	dictionary 'fsl_sky'
	record /fsl_sky/	spec_rec

	external	fsl_ftemplist
	external	fsl_invalbaryvel
	external	fsl_normal
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
	      fsl_calibrate_spectra = %loc(fsl_ftemplist)
	      call lib$signal (fsl_ftemplist, %val(3), %val(status),
     .			 spec_rec.ct_head.gmt, %val(spec_rec.attitude.pixel_no))
	      return
	   else
	      status = %loc(fsl_normal)
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
	      fsl_calibrate_spectra = %loc(fsl_ftemplist)
	      call lib$signal (fsl_ftemplist, %val(3), %val(status),
     .			 spec_rec.ct_head.gmt, %val(spec_rec.attitude.pixel_no))
	      return
	   else
	      status = %loc(fsl_normal)
	   endif
c
c  If the flag is on to set combined temperature sigmas to zero or if the number
c  of IFGs in the coadd is less than or equal to two, set the sigmas to zero.
c
	   nifgs = spec_rec.coad_spec_head.num_ifgs

	   if ((fcc_tsig0 .eq. fac_present) .or. (nifgs .le.2)) then
	      do j = 1,10
		 s_tsig(j) = 0.0
	      enddo
	   endif
c
c  Convert the temperature and sigmas to double precision and store the single
c  precision number on the output spectrum record.
c
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
	status = fsl_temporal_drift (spec_rec.ct_head.time(2), ticald,
     .				     primary_vib)
	temp(2) = temp(2) - ticald

	if (status .eq. %loc(fsl_normal)) then
c
c  Calculate the actual detector temperature and the detector responsivity and
c	time constant.
c
	   bol_volt = dble(spec_rec.coad_spec_data.bol_volt)
	   bol_cmd_bias = spec_rec.coad_spec_data.bol_cmd_bias !Two's complement
	   if (bol_cmd_bias .lt. 0) bol_cmd_bias = bol_cmd_bias + 256
	   cmd_bias = dble(bol_cmd_bias)/25.5D0
	   Tdet     = temp(fcc_chan+6)
	   status = fsl_calc_responsivity (bol_volt, cmd_bias, Tdet, Tbol, Qrad)
	   do j = 1,4
	      temp(6+j) = Tbol
	   enddo

	   if (status .eq. %loc(fsl_normal)) then
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
	            do k = 1,361
	               cvar(j,k) = fcc_varflag
	               cvar_out(j,k) = fcc_varflag
	            enddo
	         enddo
	      else
c
c  If FIL averaged variances are used, weight the FIL variance vectors by
c  the glitch rate adjusted number of IFGs for the coadd and replace the
c  individual coadd voltage variance vectors with the weighted FIL variances.
c
                 if ((nifgs .ge. 3) .or. ((nifgs .eq. 2) .and.
     .               (spec_rec.coad_spec_data.sec_template.subtracted .ne.
     .                fac_present))) then
   	            status = fsl_pcalib_variances (vvar, cvar_out)
                 else
	            do j = 1,3
	               do k = 1,361
	                  cvar_out(j,k) = fcc_varflag
	               enddo
	            enddo
                 end if
	         adjnifgs = spec_rec.coad_spec_head.adj_num_ifgs
                 if (fcc_flv .eq. fac_present) then
		    do j =1,3
	 	       do k = 1,fcc_spec_length
			  vvar(j,k) = fil_var(j,k) / adjnifgs
		       enddo
	            enddo
	         endif
c
c  Apply the calibration model to the variances for coadded ifg spectra.
c
	         status = fsl_pcalib_variances (vvar, cvar)
	      endif

	      if (status .eq. %loc(fsl_normal)) then
c
c  Apply the calibration model to the spectra. Add adjnifgs to call.
c
	         status = fsl_apply_model (primary_vib, nifgs, adjnifgs, temp, 
     .                        tsig, vspec, cvar, power, phase_corr, cspec)
	      endif


	      if ((status .eq. %loc(fsl_normal))  .and.
     .		  (fcc_sky .eq. fac_present)) then
c
c  Correct the sky spectra for the projected spacecraft barycentric velocity
c	along the Firas line of sight, checking for valid projected barycentric
c	velocities.
c
	         vbary =
     .	         dble(spec_rec.attitude.projected_barycentric_velocity)/100.0D0
	         if (dabs(vbary) .gt. vblim) then
	            status = %loc(fsl_invalbaryvel)
	            call lib$signal (fsl_invalbaryvel, %val(2),
     .				     spec_rec.ct_head.gmt,
     .				     %val(spec_rec.attitude.pixel_no))
	         elseif (vbary .ne. 0.0D0) then
	            status = fsl_doppler_shift (vbary, cspec)
	         endif
	      endif

	   endif	!  (status from fsl_calc_responsivity

	endif		!  (status from fsl_temporal_drift

C
C  Update the spectrum record.
C

	if (status .eq. %loc(fsl_normal)) then
c
c  Convert the spectrum and variances from ergs/sec/cm**2/sr/icm to MJy/sr.
c
	   do j = start_idx,stop_idx
	      cspec(j) = cspec(j) * fac_erg_to_mjy
	   enddo
	   erg_to_mjy_sq = fac_erg_to_mjy**2
	   do j = 1,3 
              do k = start_idx,stop_idx
	         cvar(j,k) = cvar(j,k) * erg_to_mjy_sq
	         cvar_out(j,k) = cvar_out(j,k) * erg_to_mjy_sq
	      enddo		  
	   enddo
c
c  Set the calibrated flag and the flag for one of the three types of variances
c  used for the autophase corrector computation.
c
           if ((fcc_flv .eq. fac_present) .or. (fcc_dvec .eq. fac_present) 
     .        .or. (nifgs .ge. 3)) then
	      spec_rec.spec_data.calibrated = 1
	   else
	      spec_rec.spec_data.calibrated = 0
           endif
	   if (fcc_flv .eq. fac_present) then
	      spec_rec.spec_data.fil_vars = 1
	   elseif (fcc_dvec .eq. fac_present) then
	      spec_rec.spec_data.fsl_vars = 1
	   else
	      spec_rec.spec_data.coadd_vars = 1
	   endif	   
c
c  Write the calibrated spectrum and variances into the spectrum record.
c  Use the start and stop frequency indices from the FSL_Model include file
c  to store only valid spectral values over the range where the model 
c  solution exists. 
c
	   do j = start_idx,stop_idx
	      display.cspec(j) = cmplx(cspec(j))
	      spec_rec.spec_data.spec(j) = display.cspec(j)
	   enddo

	   do j = start_idx,stop_idx
	      if (fcc_dvec .eq. fac_present) then
	         display.cvar(1,j) = sngl(dvector(j))
	         display.cvar(2,j) = 0.0
	         display.cvar(3,j) = 0.0
	      else
	         display.cvar(1,j) = sngl(cvar_out(1,j))
	         display.cvar(2,j) = sngl(cvar_out(2,j))
	         display.cvar(3,j) = sngl(cvar_out(3,j))
	      endif
              if ((nifgs .eq. 1) .or. ((nifgs .eq. 2) .and. 
     .            (spec_rec.coad_spec_data.sec_template.subtracted .eq. 
     .             fac_present))) then
	         spec_rec.spec_data.real_var(j)      = fcc_varflag
	         spec_rec.spec_data.imag_var(j)      = fcc_varflag
	         spec_rec.spec_data.real_imag_var(j) = fcc_varflag
              else
	         spec_rec.spec_data.real_var(j)      = display.cvar(1,j)
	         spec_rec.spec_data.imag_var(j)      = display.cvar(2,j)
	         spec_rec.spec_data.real_imag_var(j) = display.cvar(3,j)
              endif
	   enddo
c 
c  Store field values in display and output spectrum record.
c
	   display.numpts = stop_idx
	   spec_rec.spec_data.fft_length = fcc_fft_length
	   spec_rec.spec_data.lofreq_bin = start_idx
	   spec_rec.spec_data.hifreq_bin = stop_idx
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

	endif	!  (status from fsl_apply_model


	fsl_calibrate_spectra = status

	return
	end
