	integer * 4 function fsl_produce_spectra (coadd_rec, vspec, vvar,
     .						  vspec_rec, cspec_rec)

c---------------------------------------------------------------------------
c
c	Function FSL_PRODUCE_SPECTRA
c
c	This function apodizes and rotates the coadded IFG, FFTs the IFG
c	into a spectrum, and converts the spectrum from counts to
c	volts/cm**2/sr/icm.
c
c	Author:
c                FCF_Produce_Spectra
c                Gene Eplee
c		 General Sciences Corporation
c		 21 October 1992
c
c                FSL_Produce_Spectra
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fut_apod_rotl
c		fsl_read_reference
c		lib$movc3
c		lib$movc5
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
c	Changes for FCF:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 8, 1995
c       Modified FCF_Produce_Spectra to FSL_Produce_Spectra for the new FIRAS
c       pipeline which will process long spectra to get improved frequency
c       resolution.
c           1. Changed status, include files, and function names for FSL.
c           2. Removed special processing for low frequency channel FS and FL
c              scan modes. FIL now performs this processing on the ifgs.
c           3. To accomodate the longer spectra, modified spectra array sizes
c              and ifg array sizes for use in FFT, apodization and rotation
c              functions. Used new FUT_Params long coadd/spectrum record size
c              and long spectral data size to fill arrays.
c           4. Replaced call to FCF_Apodize_and_Rotate with FUT_Apod_Rotl. The
c              new function requires input parameters of ifg, corresponding
c              apodization function, FFT length, and peak position from the
c              coadd record. It returns the apodized and rotated ifg.
c           5. Performed all FFT calls using the FFT length and spectra length
c              specific to each of the 6 scan modes. The values are stored in
c              FUT_Params.
c           6. Added a check on the fakeit status from the last record to
c              determine if new reference data should be read. The FCF routine
c              already contained a check for change of microprocessor mode for
c              this purpose.
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_display)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	complex	* 16	vspec(361)	!  voltage spectrum

	integer *  4	i		!  a counter
	integer *  4	j		!  a counter
	integer *  4	last_upmode /-1/!  microprocessor change flag
	integer *  4	last_fakeit /-1/!  fakeit change flag
	integer *  4	status		!  return status
	integer *  4	t_stat		!  return status
	integer *  4    peak_pos        !  ifg peak position from coadd record

	real	*  4	aifg(512)	!  ifg from coadd record
	real	*  8	imag_var	!  imaginary-to-imaginary variance
	real	*  8	real_var	!  real-to-real variance
	real	*  8	reim_var	!  read-to-imaginary variance
	real	*  8	rifg(720)	!  apodized, rotated ifg
	real	*  8	sample		!  average ifg sample
	real	*  8	vvar(3,361)	!  apodized voltage variances
	real	*  8	zifg(720)	!  transformed ifg

	save last_upmode
	save last_fakeit

	integer *  4	fut_apod_rotl
	integer *  4	fsl_read_reference

	dictionary 'fil_sky'
	record /fil_sky/ coadd_rec

	dictionary 'fsl_sky'
	record /fsl_sky/ cspec_rec
	record /fsl_sky/ vspec_rec

	external	fsl_normal
	external	fut_normal

C
C  Initialize the routine.
C

c
c  Initialize the display block.
c
	call lib$movc5 (0,,0,25254,display)
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
c  Read the reference data for a new microprocessor mode or fakeit status.
c
	if ((last_upmode .ne. display.upmode) .or.
     .      (last_fakeit .ne. display.fakeit)) then
	   last_upmode = display.upmode
	   last_fakeit = display.fakeit
	   status = fsl_read_reference (display.fakeit, display.upmode)
	   if (status .ne. %loc(fsl_normal)) then
	      fsl_produce_spectra = status
	      return
	   endif
	endif

c
c  Extract the full 512 point IFG from the coadd record.
c
	call lib$movc5 (0,,0,2048,aifg)

	do j= 1,512
	    display.aifg(j) = coadd_rec.coad_data.ifg(j)
	    aifg(j) = dble(coadd_rec.coad_data.ifg(j))
	enddo

c
c  Extract the full variances from the coadd record.
c
	call lib$movc5 (0,,0,8664,vvar)

	do j = 1,361
	     display.vvar(1,j) = coadd_rec.coad_data.real_var(j)
	     display.vvar(2,j) = coadd_rec.coad_data.imag_var(j)
	     display.vvar(3,j) = coadd_rec.coad_data.real_imag_var(j)
	     do i = 1,3
	        vvar(i,j) = dble(display.vvar(i,j))
	     enddo
	enddo


C
C  Generate the voltage spectrum.
C

c
c  Apodize and rotate the IFG in preparation for the FFT.
c
	peak_pos = coadd_rec.coad_spec_data.peak_pos

	status = fut_apod_rotl (aifg, apod_fcn, fcc_fft_length, peak_pos, rifg)

	if (status .eq. %loc(fut_normal)) then
c
c  Do the fft, unpacking the complex spectrum from the real transformed array.
c  The fft length and spectra length are specific to the 6 scan modes.
c
	   call lib$movc5 (0,,0,5760,zifg)
	   call lib$movc5 (0,,0,5776,vspec)

	   call dfftrf (fcc_fft_length, rifg, zifg)
	   vspec(1) = dcmplx(zifg(1),0.0D0)
	   display.spec(1) = cmplx(vspec(1))
	   do i=2,fcc_fft_length-2,2
	         j=i/2+1
	         vspec(j) = dcmplx(zifg(i),zifg(i+1))
	         display.spec(j) = cmplx(vspec(j))
	   enddo
	   vspec(fcc_spec_length) = dcmplx(zifg(fcc_fft_length),0.0D0)
	   display.spec(fcc_spec_length) = cmplx(vspec(fcc_spec_length))

c
c  Normalize the spectrum and  convert it from counts to volts.
c  Shift the phase of the high frequency short fast spectra.
c  The units of the normalized spectrum are volts/cm**2/sr/icm.
c
	   vspec(1) = dcmplx(0.0D0,0.0D0)
	   if (fcc_xhsf .eq. fac_present) then
	      do j = 2,fcc_spec_length
	         vspec(j) = vspec(j) / spec_norm / ztrans(j)
     .						 * cdexp(dble(j-1)*phase_shift)
	      enddo
	   else
	      do j = 2,fcc_spec_length
	         vspec(j) = vspec(j) / spec_norm / ztrans(j)
	      enddo
	   endif


C
C  Fill in the voltage spectrum record and the calibrated spectrum record.
C  The following lines of code copy the coadd data into the voltage spectrum
C  and calibrated spectrum records, initialize voltage and calibrated spectral
C  data group blocks to zero, and fill the voltage spectra record and display
C  include file with spectral data . The new long spectra sizes are used.
C
	   call lib$movc3 (fac_coad_spec_sizel,coadd_rec,vspec_rec)
	   call lib$movc5 (0,,0,fac_spec_data_sizel,vspec_rec.spec_data)
	   cspec_rec = vspec_rec

	   do j = 1,fcc_spec_length
	      display.vspec(j)            = cmplx(vspec(j))
	      vspec_rec.spec_data.spec(j) = display.vspec(j)
	      display.vvar(1,j) = sngl(vvar(1,j))
	      display.vvar(2,j) = sngl(vvar(2,j))
	      display.vvar(3,j) = sngl(vvar(3,j))
	      vspec_rec.spec_data.real_var(j)      = display.vvar(1,j)
	      vspec_rec.spec_data.imag_var(j)      = display.vvar(2,j)
	      vspec_rec.spec_data.real_imag_var(j) = display.vvar(3,j)
	   enddo

	endif			!	(fsl_apodize_and_rotate


	fsl_produce_spectra = status

	return
	end
