	integer * 4 function  fip_write_model ()

c-------------------------------------------------------------------------------
c
c	Function FIP_WRITE_MODEL
c
c	This function convrts the FIRAS calibration model solution and
c	reference datasets into the Initial Product format.  It then writes
c	the Initial Product model solution to the RMS binary file
c	CSDR$FIRAS_OUT:FIP_MOD_CCSS.VVV_XXXXXXXXXX, where CCSS refers to
c	channel/scan mode and VVV_XXXXXXXXXX is the calibration model solution
c	file name extension.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  9 June 1993
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
c		lib$signal
c		str$trim
c
c	Include files:
c		fip_config_model.txt
c		fip_invoc_model.txt
c		fip_frequency.txt
c		fip_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Added FEX_CVS reference dataset.  Gene Eplee, GSC, 29 September 1993.
c
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fip_invoc_model)'
	include '(fip_config_model)'
	include '(fip_frequency)'
	include '(fip_model)'
	include '(fut_params)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	mod_lun		!  model solution file lun
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmswrite
	external	fut_normal

C
C  Write the primary model solution fields to the Initial Product RDL.
C

c
c  The model solution identification fields.
c
	fip_model.model_label = fex_model.mod_head.label
	fip_model.chanscan    = fcc_scan_mode
	fip_model.nu_zero     = fcc_nu0
	fip_model.delta_nu    = fcc_dnu
	fip_model.omega_zero  = fcc_w0
	fip_model.delta_omega = fcc_dw
	fip_model.num_freq    = fcc_nfreq

c
c  The detector response function.
c
	fip_model.normalization = spec_norm
	fip_model.time_constant = sngl(tau)
	fip_model.dc_response   = sngl(S0)

c
c  The gain functions and the detector noise.
c	Take the square root of the C-Vector^2 and convert the result from
c	eplees to MJy/Sr.
c
	do j = fcc_jlo,fcc_jhi
	   fip_model.relex_gain(j-fcc_freq_offset) = real(cmplx(ztrans(j)))
	   fip_model.ielex_gain(j-fcc_freq_offset) = aimag(cmplx(ztrans(j)))
	   fip_model.rtransfer(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.transfer(j)))
	   fip_model.itransfer(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.transfer(j)))
	   fip_model.detector_noise(j-fcc_freq_offset) = 
     .          sngl(dsqrt(fex_cvs.cvector(j))) * fac_erg_to_mjy
	enddo

c
c  The apodization function.
c
	do j = 1,512
	   fip_model.apodization(j) = sngl(apod_fcn(j))
	enddo


C
C  Write the secondary model solution fields to the Initial Product RDL.
C

c
c  The bolometer parameters.
c
	fip_model.bolparm_R0   = sngl(fex_model.bolparm(1))
	fip_model.bolparm_T0   = sngl(fex_model.bolparm(2))
	fip_model.bolparm_G1   = sngl(fex_model.bolparm(3))
	fip_model.bolparm_beta = sngl(fex_model.bolparm(4))
	fip_model.bolparm_rho  = sngl(fex_model.bolparm(5))
	fip_model.bolparm_C1   = sngl(fex_model.bolparm(7))
	fip_model.bolparm_C3   = sngl(fex_model.bolparm(6))
	fip_model.bolparm_JO   = sngl(fex_model.bolparm(8))
	fip_model.bolparm_JG   = sngl(fex_model.bolparm(9))
	fip_model.bolparm_RL   = sngl(fac_load_resist)

c
c  The mission average bolometer state.
c
	fip_model.bolom_bias        = sngl(cmd_bias)
	fip_model.bolom_voltage     = sngl(bol_volt)
	fip_model.bolom_bathtemp    = sngl(Tdet)
	fip_model.bolom_temperature = sngl(Tbol)

c
c  The ical drift correction and xcal offset.
c
	fip_model.drift_amp       = sngl(fex_model.bolparm(12))
	fip_model.drift_tc        = 1.0/sngl(fex_model.bolparm(13))
	fip_model.drift_offset    = sngl(fex_model.bolparm(14))
	fip_model.xcal_correction = sngl(fex_model.bolparm(21))

c
c  The harmonic and vibration correction coefficients
c
	fip_model.optparm_2H  = sngl(fex_model.bolparm(11))
	fip_model.optparm_3H  = sngl(fex_model.bolparm(10))
	fip_model.optparm_VP0 = sngl(fex_model.bolparm(16))
	fip_model.optparm_VP1 = sngl(fex_model.bolparm(17))
	fip_model.optparm_VP2 = sngl(fex_model.bolparm(18))
	fip_model.optparm_VP3 = sngl(fex_model.bolparm(19))
	fip_model.optparm_VP4 = sngl(fex_model.bolparm(20))
	fip_model.optparm_VS  = sngl(fex_model.bolparm(15))

c
c  The optical model emissivities.
c
	do j = fcc_jlo,fcc_jhi
	   fip_model.rical(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.ical(j)))
	   fip_model.iical(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.ical(j)))
	   fip_model.rskyhorn(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.skyhorn(j)))
	   fip_model.iskyhorn(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.skyhorn(j)))
	   fip_model.rrefhorn(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.refhorn(j)))
	   fip_model.irefhorn(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.refhorn(j)))
	   fip_model.rdihedral(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.dihedral(j)))
	   fip_model.idihedral(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.dihedral(j)))
	   fip_model.rstructure(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.struct(j)))
	   fip_model.istructure(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.struct(j)))
	   fip_model.rbolometer(j-fcc_freq_offset) =
     .		real(cmplx(fex_model.bolometer(j)))
	   fip_model.ibolometer(j-fcc_freq_offset) =
     .		aimag(cmplx(fex_model.bolometer(j)))
	enddo


C
C  Write the Initial Product model solution to an RMS file.
C

	status = %loc(fip_normal)
c
c  Find the model solution file name.
c
	fcc_fip_file = 'CSDR$FIRAS_OUT:FIP_MOD_' // fcc_scan_mode // '.' //
     .			fcc_model_ext
	call str$trim (fcc_fip_file, fcc_fip_file, fcc_fiplen)

c
c  Open the model solution file.
c
	rstatus = fut_get_lun(mod_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=mod_lun, file=fcc_fip_file, status='new',
     .	      form='unformatted', recordtype='fixed', recl=fcc_fip_size,
     .	      access='sequential', iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Write the model solution.
c
	   write (mod_lun, iostat=io_stat) fip_model
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmswrite)
	      call lib$signal (fip_rmswrite, %val(2),
     .			       fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	   endif

c
c  Close the model solution file.
c
	   close (unit=mod_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(mod_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2),
     .			    fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	endif	!  (open status


	fip_write_model = status

	return
	end
