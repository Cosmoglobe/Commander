	integer * 4 function fcf_produce_spectra (coadd_rec, vspec, vvar,
     .						  vspec_rec, cspec_rec)

c---------------------------------------------------------------------------
c
c	Function FCF_PRODUCE_SPECTRA
c
c	This function apodizes and rotates the coadded IFG, FFTs the IFG into
c	a spectrum, and converts the spectrum from counts to
c	volts/cm**2/sr/icm.  It also apodizes the voltage variances.
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  21 October 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		coadd_rec			coadd record
c
c	Output:
c		vspec		complex * 16	voltage spectrum
c		vvar(3)		real * 8	apdoized voltage variances
c		vspec_rec			voltage spectrum record
c		cspec_rec			calibrated spectrum record
c
c	Subroutines called:
c		dfftrf
c		fcf_apodize_and_rotate
c		fcf_read_reference
c		lib$movc3
c		lib$movc5
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
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_display)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	complex	* 16	vspec(257)	!  voltage spectrum

	integer *  4	i		!  a counter
	integer *  4	j		!  a counter
	integer *  4	last_upmode /-1/!  microprocessor change flag
	integer *  4	status		!  return status
	integer *  4	t_stat		!  return status

	real	*  8	daifg(512)	!  ifg
	real	*  8	imag_var	!  imaginary-to-imaginary variance
	real	*  8	ivar(3,257)	!  unapodized voltage variances
	real	*  8	real_var	!  real-to-real variance
	real	*  8	reim_var	!  read-to-imaginary variance
	real	*  8	rifg(512)	!  apodized, rotated ifg
	real	*  8	sample		!  average ifg sample
	real	*  8	s2		!  sinc^2 function
	real	*  8	vvar(3,257)	!  apodized voltage variances
	real	*  8	zifg(512)	!  transformed ifg

	save last_upmode

	integer *  4	fcf_apodize_and_rotate
	integer *  4	fcf_read_reference

	dictionary 'fic_sky'
	record /fic_sky/ coadd_rec

	dictionary 'fcf_sky'
	record /fcf_sky/ cspec_rec
	record /fcf_sky/ vspec_rec

	external	fcf_normal

C
C  Initialize the routine.
C

c
c  Initialize the display block.
c
	call lib$movc5 (0,,0,21678,display)
	display.gmt         = coadd_rec.ct_head.gmt
	display.coadd_label = coadd_rec.coad_spec_head.label
	display.fakeit      = coadd_rec.coad_spec_data.fakeit
	display.gain        = coadd_rec.coad_spec_data.gain_sum /
     .			      coadd_rec.coad_spec_head.num_ifgs
	display.sweeps      = coadd_rec.coad_spec_data.sweeps
	display.upmode      = coadd_rec.coad_spec_data.sci_mode
	display.xcal_pos    = coadd_rec.coad_spec_data.xcal_pos
	display.pixel_no    = coadd_rec.attitude.pixel_no

c
c  Read the reference data for a new microprocessor mode.
c
	if (last_upmode .ne. display.upmode) then
	   last_upmode = display.upmode
	   status = fcf_read_reference (display.fakeit, display.upmode)
	   if (status .ne. %loc(fcf_normal)) then
	      fcf_produce_spectra = status
	      return
	   endif
	endif

c
c  Extract the IFG from the coadd record.
c
	call lib$movc5 (0,,0,4096,daifg)
	if (fcc_xlsf .eq. fac_present) then
c
c    Shift the low frequency short fast data left by three samples.
c    Average four adjacent samples to extract a 128-point ifg.
c    Point 128 = 0.0
c
	   do j=1,127
	      sample = 0.0D0
	      do i=0,3
	         sample = sample + dble(coadd_rec.coad_data.ifg(4*j+i))
	      enddo
	      daifg(j) = sample/4.0D0
	      display.aifg(j) = sngl(daifg(j))
	   enddo
	elseif (fcc_xllf .eq. fac_present) then
c
c    Extract the first 128 points of the ifg to match the xlsf data.
c    Point 128 = 0.0
c
	   do j=1,127
	      display.aifg(j) = coadd_rec.coad_data.ifg(j)
	      daifg(j) = dble(coadd_rec.coad_data.ifg(j))
	   enddo
	else
c
c  Extract the full 512-opoint ifg.
c
	   do j=1,512
	      display.aifg(j) = coadd_rec.coad_data.ifg(j)
	      daifg(j) = dble(coadd_rec.coad_data.ifg(j))
	   enddo
	endif

c
c  Extract the variances from the coadd record.
c
	call lib$movc5 (0,,0,6168,ivar)
	if (fcc_xlsf .eq. fac_present) then
c
c  Convolve the low frequency short fast variances with a sinc^2 function to
c    match the smoothed low frequency long fast data.
c
	   do j = 2,65
	      real_var = 0.0D0
	      imag_var = 0.0D0
	      reim_var = 0.0D0
	      do i = 2,257
	         if (i .eq. j) then
	            s2 = 1.0D0
	         else
	            s2 = (dsin(fac_dpi*dble(j-i)/4.0D0) /
     .			 (fac_dpi*dble(j-i)/4.0D0))**2
	         endif
	         real_var = real_var + s2*dble(coadd_rec.coad_data.real_var(i))
	         imag_var = imag_var + s2*dble(coadd_rec.coad_data.imag_var(i))
	         reim_var = reim_var +
     .			    s2*dble(coadd_rec.coad_data.real_imag_var(i))
	      enddo
	      ivar(1,j) = real_var/4.0D0
	      ivar(2,j) = imag_var/4.0D0
	      ivar(3,j) = reim_var/4.0D0
	      do i = 1,3
	         display.ivar(i,j) = sngl(ivar(i,j))
	      enddo
	   enddo
	elseif (fcc_xllf .eq. fac_present) then
c
c  Resample the low frequency long fast variances to smooth them to the low
c    frequency short fast spectral resolution.
c
	   do j = 2,65
	      display.ivar(1,j) = coadd_rec.coad_data.real_var(4*j-3)/4.0
	      display.ivar(2,j) = coadd_rec.coad_data.imag_var(4*j-3)/4.0
	      display.ivar(3,j) = coadd_rec.coad_data.real_imag_var(4*j-3)/4.0
	      do i = 1,3
	         ivar(i,j) = dble(display.ivar(i,j))
	      enddo
	   enddo
	else
c
c  Extract the full variances.
c
	   do j = 1,257
	      display.ivar(1,j) = coadd_rec.coad_data.real_var(j)
	      display.ivar(2,j) = coadd_rec.coad_data.imag_var(j)
	      display.ivar(3,j) = coadd_rec.coad_data.real_imag_var(j)
	      do i = 1,3
	         ivar(i,j) = dble(display.ivar(i,j))
	      enddo
	   enddo
	endif


C
C  Generate the voltage spectrum.
C

c
c  Apodize and rotate the IFG in preparation for the FFT.
c	Apodize the voltage variances.
c
	status = fcf_apodize_and_rotate (daifg, ivar, rifg, vvar)

	if (status .eq. %loc(fcf_normal)) then
c
c  Do the fft, unpacking the complex spectrum from the real transformed array.
c
	   call lib$movc5 (0,,0,4096,zifg)
	   call lib$movc5 (0,,0,4112,vspec)
	   if ((fcc_xlsf .eq. fac_present)  .or.
     .	       (fcc_xllf .eq. fac_present)) then
c
c  Do a 128-point fft for the low frequency fast ifgs.
c
	      call dfftrf (128,rifg,zifg)
	      vspec(1) = dcmplx(zifg(1),0.0D0)
	      display.spec(1) = cmplx(vspec(1))
	      do i=2,126,2
	         j=i/2+1
	         vspec(j) = dcmplx(zifg(i),zifg(i+1))
	         display.spec(j) = cmplx(vspec(j))
	      enddo
	      vspec(65) = dcmplx(zifg(128),0.0D0)
	      display.spec(65) = cmplx(vspec(65))
	   else
c
c  Do a 512-point fft for the rest of the data.
c
	      call dfftrf (512,rifg,zifg)
	      vspec(1) = dcmplx(zifg(1),0.0D0)
	      display.spec(1) = cmplx(vspec(1))
	      do i=2,510,2
	         j=i/2+1
	         vspec(j) = dcmplx(zifg(i),zifg(i+1))
	         display.spec(j) = cmplx(vspec(j))
	      enddo
	      vspec(257) = dcmplx(zifg(512),0.0D0)
	      display.spec(257) = cmplx(vspec(257))
	   endif

c
c  Normalize the spectrum and  convert it from counts to volts.
c    Fix the short fft for the low frequency fast spectra.
c    Shift the phase of the high frequency short fast spectra.
c    The units of the normalized spectrum are volts/cm**2/sr/icm.
c
	   vspec(1) = dcmplx(0.0D0,0.0D0)
	   if ((fcc_xlsf .eq. fac_present)  .or.
     .	       (fcc_xllf .eq. fac_present)) then
	      do j = 2,65
	         vspec(j) = 4.0D0 * vspec(j) / spec_norm / ztrans(j)
	      enddo
	   elseif (fcc_xhsf .eq. fac_present) then
	      do j = 2,257
	         vspec(j) = vspec(j) / spec_norm / ztrans(j)
     .						 * cdexp(dble(j-1)*phase_shift)
	      enddo
	   else
	      do j = 2,257
	         vspec(j) = vspec(j) / spec_norm / ztrans(j)
	      enddo
	   endif


C
C  Fill in the voltage spectrum record and the calibrated spectrum record.
C
	   call lib$movc3 (fac_coad_spec_size,coadd_rec,vspec_rec)
	   call lib$movc5 (0,,0,fac_spec_data_size,vspec_rec.spec_data)
	   if (fcc_xllf .eq. fac_present) then
	      vspec_rec.coad_spec_data.mtm_length     = fcc_length
	      vspec_rec.coad_spec_data.adds_per_group = fcc_ngroup
	      vspec_rec.coad_spec_data.nyquist_icm    = sngl(fnyq_icm)
	      vspec_rec.coad_spec_data.nyquist_hertz  = sngl(fnyq_hz)
	   endif
	   cspec_rec = vspec_rec

	   do j = 1,257
	      display.vspec(j)            = cmplx(vspec(j))
	      vspec_rec.spec_data.spec(j) = display.vspec(j)
	      display.vvar(1,j) = sngl(vvar(1,j))
	      display.vvar(2,j) = sngl(vvar(2,j))
	      display.vvar(3,j) = sngl(vvar(3,j))
	      vspec_rec.spec_data.real_var(j)      = display.vvar(1,j)
	      vspec_rec.spec_data.imag_var(j)      = display.vvar(2,j)
	      vspec_rec.spec_data.real_imag_var(j) = display.vvar(3,j)
	   enddo

	endif			!	(fcf_apodize_and_rotate


	fcf_produce_spectra = status

	return
	end
