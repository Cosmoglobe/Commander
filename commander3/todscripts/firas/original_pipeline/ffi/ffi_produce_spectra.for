	integer * 4 function ffi_produce_spectra (coadd_rec, fish_rec)

c---------------------------------------------------------------------------
c
c	Function FFI_PRODUCE_SPECTRA
c
c	This function apodizes and rotates the coadded IFG, FFT's the IFG into
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
c		coadd_rec			coadd record
c		fish_rec			spectrum record
c
c	Output:
c		none
c
c	Subroutines called:
c		dfftrf
c		ffi_apodize_and_rotate
c		ffi_read_reference
c		fut_temp_list
c		lib$movc5
c
c	Include files:
c		ffi_config.txt
c		ffi_invoc.txt
c		ffi_spec.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffi_config)'
	include '(ffi_invoc)'
	include '(ffi_spec)'

	complex * 16	vspec(257)	!  double precision voltage spectrum

	integer *  2	cmd_bias	!  integer valued commanded bias

	integer *  4	fakeit		!  fakeit flag
	integer *  4	i		!  a counter
	integer *  4	j		!  a counter
	integer *  4	last_upmode /-1/!  microprocessor change flag
	integer *  4	status		!  return status
	integer *  4	upmode 		!  microprocessor mode

	real	*  4	temp(10)	!  instrument temperatures
	real	*  4	tsigma(10)	!  temperature sigmas

	real	*  8	daifg(512)	!  double precision ifg
	real	*  8	imag_var	!  imaginary-to-imaginary variances
	real	*  8	real_var	!  real-to-real variances
	real	*  8	reim_var	!  real-to-imaginary variances
	real	*  8	rifg(512)	!  apodized, rotated IFG
	real	*  8	sample		!  average ifg sample
	real	*  8	s2		!  sinc^2 function
	real	*  8	vsigma(257)	!  voltage sigmas
	real	*  8	vvar(3,257)	!  apodized voltage variances
	real	*  8	zifg(512)	!  double precision transformed ifg

	save last_upmode

	integer *  4	ffi_apodize_and_rotate
	integer *  4	ffi_read_reference
	integer *  4	fut_temp_list

	dictionary 'fic_sky'
	record /fic_sky/ coadd_rec
	record /fish_spec/ fish_rec

	external	ffi_ftemplist
	external	ffi_normal

C
C  Initialize the routine.
C

c
c  Initialize the spectrum record.
c
	call lib$movc5 (0,,0,4256,fish_rec)
	fish_rec.nifgs = coadd_rec.coad_spec_head.num_ifgs
	fish_rec.time = sngl(dble(coadd_rec.ct_head.time(2) - fac_apco_date)
     .			/ fac_vax_year_len)
	fish_rec.volt = dble(coadd_rec.coad_spec_data.bol_volt)
	cmd_bias  = coadd_rec.coad_spec_data.bol_cmd_bias	!  Two's
	if (cmd_bias .lt. 0) cmd_bias = cmd_bias + 256		!  complement
	fish_rec.bias = dble(cmd_bias)/25.5D0
	fish_rec.gain = coadd_rec.coad_spec_data.gain_sum /
     .			coadd_rec.coad_spec_head.num_ifgs
	fakeit = coadd_rec.coad_spec_data.fakeit
	upmode = coadd_rec.coad_spec_data.sci_mode

c
c  Extract the IFG from the coadd record.
c
	call lib$movc5 (0,,0,4096,daifg)
	if (fcc_xlsf .eq. fac_present) then
c
c    Shift the low frequency short fast data left by three samples.
c    Average four adjacent samples to extract a 128-point ifg.
c    Point 128 = 0.0.
c
	   do j=1,127
	      sample = 0.0D0
	      do i=0,3
	         sample = sample + dble(coadd_rec.coad_data.ifg(4*j+i))
	      enddo
	      daifg(j) = sample/4.0D0
	   enddo
	elseif (fcc_xllf .eq. fac_present) then
c
c    Extract the first 128 points of the ifg to match the xlsf data.
c    Point 128 = 0.0
c
	   do j=1,127
	      daifg(j) = dble(coadd_rec.coad_data.ifg(j))
	   enddo
	else
c
c  Extract the full 512-point ifg.
c
	   do j=1,512
	      daifg(j) = dble(coadd_rec.coad_data.ifg(j))
	   enddo
	endif

c
c  Extract the variances from the coadd record.
c
	call lib$movc5 (0,,0,6168,vvar)
	if (fcc_xlsf .eq. fac_present) then
c
c  Convolve the low frequency short fast variances with a sinc**2 function to
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
	      vvar(1,j) = real_var/4.0D0
	      vvar(2,j) = imag_var/4.0D0
	      vvar(3,j) = reim_var/4.0D0
	   enddo
	elseif (fcc_xllf .eq. fac_present) then
c
c  Resample the low frequency long fast variances to smooth them to the low
c    frequency short fast spectral resolution.
c
	   do j = 2,65
	      vvar(1,j) = 4.0D0*dble(coadd_rec.coad_data.real_var(4*j-3))
	      vvar(2,j) = 4.0D0*dble(coadd_rec.coad_data.imag_var(4*j-3))
	      vvar(3,j) = 4.0D0*dble(coadd_rec.coad_data.real_imag_var(4*j-3))
	   enddo
	else
c
c  Extract the full variances.
c
	   do j = 1,257
	      vvar(1,j) = dble(coadd_rec.coad_data.real_var(j))
	      vvar(2,j) = dble(coadd_rec.coad_data.imag_var(j))
	      vvar(3,j) = dble(coadd_rec.coad_data.real_imag_var(j))
	   enddo
	endif

c
c  Read the reference data for a new microprocessor mode.
c
	if (last_upmode .ne. upmode) then
	   last_upmode = upmode
	   status = ffi_read_reference (fakeit, upmode)
	   if (status .ne. %loc(ffi_normal)) then
	      ffi_produce_spectra = status
	      return
	   endif
	endif

c
c  Combine the instrument temperatures.
c

	status = fut_temp_list (coadd_rec.en_analog, coadd_rec.en_sigma,
     .				config.grtcoawt, config.grttrans,
     .				0, .false., temp, tsigma)
	if (status .ne. 0) then
	   ffi_produce_spectra = %loc(ffi_ftemplist)
	   call lib$signal (ffi_ftemplist, %val(2), %val(status),
     .			    coadd_rec.ct_head.gmt)
	   return
	else
	   do j=1,6
	      fish_rec.temp(j) = temp(j)
	      fish_rec.tsigma(j) = tsigma(j)
	   enddo
	   fish_rec.temp(7) = temp(fcc_chan+6)
	   fish_rec.tsigma(7) = tsigma(fcc_chan+6)
	endif


C
C  Generate the voltage spectrum.
C

c
c  Apodize and rotate the IFG in preparation for the FFT.
c	Apodize the voltage variances and compute the voltage sigmas.
c
	status = ffi_apodize_and_rotate (daifg, vvar, rifg, vsigma)

	if (status .eq. %loc(ffi_normal)) then
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
	      do i=2,126,2
	         j=i/2+1
	         vspec(j) = dcmplx(zifg(i),zifg(i+1))
	      enddo
	      vspec(65) = dcmplx(zifg(128),0.0D0)
	   else
c
c  Do a 512-point fft for the rest of the data.
c
	      call dfftrf (512,rifg,zifg)
	      vspec(1) = dcmplx(zifg(1),0.0D0)
	      do i=2,510,2
	         j=i/2+1
	         vspec(j) = dcmplx(zifg(i),zifg(i+1))
	      enddo
	      vspec(257) = dcmplx(zifg(512),0.0D0)
	   endif

c
c  Normalize the spectrum and convert it from counts to volts.
c    Fix the short fft for the low frequency fast spectra.
c    Shift the phase of the high frequency short fast spectra.
c
	   vspec(1) = dcmplx(0.0D0,0.0D0)
	   if ((fcc_xlsf .eq. fac_present)  .or.
     .	       (fcc_xllf .eq. fac_present)) then
	      do j = 2,65
	         vspec(j) = 4.0D0 * vspec(j) / fac_adc_scale / ztrans(j)
	      enddo
	   elseif (fcc_xhsf .eq. fac_present) then
	      do j = 2,257
	         vspec(j) = vspec(j) / fac_adc_scale / ztrans(j)
     .						 * cdexp(dble(j-1)*phase_shift)
	      enddo
	   else
	      do j = 2,257
	         vspec(j) = vspec(j) / fac_adc_scale / ztrans(j)
	      enddo
	   endif


C
C  Fill in the voltage spectrum record.
C
	   do j = 2,176
	      fish_rec.vspec(j-1) = vspec(j)
	      fish_rec.vsigma(j-1) = vsigma(j)
	   enddo

	endif			!	(ffi_apodize_and_rotate


	ffi_produce_spectra = status

	return
	end
