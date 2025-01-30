	integer * 4 function  fsl_display_spectra ()

c-------------------------------------------------------------------------------
c
c	Function FSL_DISPLAY_SPECTRA
c
c	This function plots the IFGs, spectra, and sigmass that are produced
c	by FSL.
c
c	Author:  
c                FCF_Display_Spectra
c                Gene Eplee
c		 General Sciences Corp.
c		 12 March 1993
c
c                FSL_Display_Spectra
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fut_plot
c		fut_plot_title
c		fut_setxaxl
c		lib$movc5
c		ots$cvt_ti_l
c		str$trim
c		str$upcase
c
c	Include files:
c		fsl_config.txt
c		fsl_display.txt
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Add reference spectrum plotting capability.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 10, 1995 
c       Modified FCF_Display_Spectra to FSL_Display_Spectra for the new FIRAS 
c       pipeline which will process long spectra to get improved frequency 
c       resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed call to new FUT_Setxaxl to handle longer spectra. 
c              Modified the calling sequence from old routine.
c           3. Used number of frequency points corresponding to model solution
c              to extract spectral data for plots and to 
c           4. Removed obsolete plots for unapodized variances and renumbered
c              plot options.
c           5. Changed units to MJy/sr for final calibrated spectra and
c              variances.

c------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_display)'
	include '(fsl_invoc)'

	character *  16  answer			!  query response
	character *  60  label			!  coadd label
	character *  32  plot_device		!  plot device
	character * 100  plot_label		!  plot label
	character *  64	 plt_com_file		!  plt command file name
	character * 100  title(3)		!  plot title

	complex *  8	dspec(512)		!  data array for plotting

	integer *  4	fakeit			!  fakeit flag
	integer *  4    first_sample / 1 /      !  first (start) sample of ifg
	integer *  4	gain			!  commanded gain
	integer *  4    npts                    !  number of points to display
	integer *  4	j			!  a counter
	integer *  4	nans			!  number of answer
	integer *  4	numlen			!  position of number in answer
	integer *  4	plt_com			!  plt command flag
	integer *  4    scale0 / 0 /            !  coadd type of X axis scale 
	integer *  4    scale1 / 1 /            !  spectrum type of X axis scale
	integer *  4	status			!  return status
	integer *  4	sweeps			!  number of mtm sweeps
	integer *  4	upmode			!  microprocessor mode
	integer *  4	xcal_pos		!  xcal position flag
	integer *  4	zp			!  zero point flag

	real	*  4	spacex			!  tickmark spacing for spectra
	real	*  4	spacexi			!  tickmark spacing for ifgs
	real	*  4	startx			!  starting tickmark for spectra
	real	*  4	startxi			!  starting tickmark for ifgs

	integer *  4	fut_setxaxl
	integer *  4	ots$cvt_ti_l

	external	fut_normal

C
C  Initialize the routine.
C

c
c  Initialize the plot.
c
	fakeit       = display.fakeit
	gain         = display.gain
	label        = display.coadd_label
	sweeps       = display.sweeps
	upmode       = display.upmode
	xcal_pos     = display.xcal_pos
	zp           = fac_present
	plot_device  = fcc_plot_device
	plt_com      = fcc_plt_com
	plt_com_file = fcc_plt_com_file
	npts = display.numpts

c
c  Call FUT_Setxaxl to get start point and spacing for plot. Scale is set to
c  zero for a coadded ifg and one for a spectrum. FSL processes only ifgs
c  with first (start) sample of 1.

	status = fut_setxaxl (fcc_chan, fcc_smode, fakeit, upmode, 
     .                        fcc_ngroup, config.nyquist, scale0, 
     .                        first_sample, startxi, spacexi)
	status = fut_setxaxl (fcc_chan, fcc_smode, fakeit, upmode, 
     .                        fcc_ngroup, config.nyquist, scale1, 
     .                        first_sample, startx, spacex)
	if (status .ne. %loc(fut_normal)) then
	   fsl_display_spectra = status
	   return
	endif

c
c  Get the user input.
c
	fsl_display_spectra = fac_present
	answer       = 'Y'
	do while (answer(1:1) .ne. 'Q')

	   type 10, fcc_scan_mode, display.gmt, display.pixel_no
  10	   format (/, 1x, 'Plots of ', a, ' spectrum ', a, ', pixel ', i4, /,
     .		    3x, 'Choose what to view (numbers or capital letters): ', /,
     .		    5x, '  (1) Interferogram', /,
     .		    5x, '  (2) Spectrum in counts', /,
     .		    5x, '  (3) Voltage spectrum', /,
     .		    5x, '  (4) Differential spectrum', /,
     .		    5x, '  (5) Reference spectrum', /,
     .		    5x, '  (6) Calibrated spectrum', /,
     .              5x, '  (7) Apodized voltage spectrum Sigmas', /,
     .		    5x, '  (8) Apodized voltage spectrum Cross Variances', /,
     .		    5x, '  (9) Calibrated spectrum Sigmas', /,
     .		    5x, ' (10) Calibrated spectrum Cross Variances', /,
     .		    5x, ' (11) Jump n spectra', /,
     .		    5x, ' (12) Quit displaying plots', /,
     .		    5x, 'Hit return to continue processing ==> ', $)

	   accept 20, answer
  20	   format(a)

	   call str$trim (answer, answer, numlen)
	   call str$upcase (answer, answer)
	   status = ots$cvt_ti_l (answer(1:numlen), nans, %val(4), 1)
	   if (.not. status) then
	      nans = 0
	   endif
	   call lib$movc5 (0,,0,4096,dspec)	!  Initialize the data array.


C
C  Do the plotting.
C

	   if (nans .eq. 1  .or.  answer(1:1) .eq. 'I') then
c
c  Plot the coadded IFG.
c
	      do j = 1, 512
	         dspec(j) = cmplx(display.aifg(j),0.0)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Coadded Calibration Interferogram'
	      else
	         plot_label = 'Coadded Sky Interferogram'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, 512, startxi, spacexi, title,
     .			     'Optical Path Difference (cm)', 'Counts',
     .			     fac_not_present, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 2  .or.  answer(1:1) .eq. 'S') then
c
c  Plot the spectrum in counts.
c
	      do j = 1, npts
	         dspec(j) = display.spec(j)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Calibration Spectrum in Counts'
	      else
	         plot_label = 'Sky Spectrum in Counts'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'Counts',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 3  .or.  answer(1:1) .eq. 'V') then
c
c  Plot the voltage spectrum.
c
	      do j = 1, npts
	         dspec(j) = display.vspec(j)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Voltage Calibration Spectrum'
	      else
	         plot_label = 'Voltage Sky Spectrum'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'Volts/Cm**2/Sr/Icm',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 4  .or.  answer(1:1) .eq. 'D') then
c
c  Plot the differential spectrum.
c
	      do j = 1, npts
	         dspec(j) = display.espec(j)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Differential Calibration Spectrum'
	      else
	         plot_label = 'Differential Sky Spectrum'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'Ergs/Sec/Cm**2/Sr/Icm',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 5  .or.  answer(1:1) .eq. 'R') then
c
c  Plot the reference spectrum.
c
	      do j = 1, npts
	         dspec(j) = -display.pspec(j)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Reference Calibration Spectrum'
	      else
	         plot_label = 'Reference Sky Spectrum'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'Ergs/Sec/Cm**2/Sr/Icm',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 6  .or.  answer(1:1) .eq. 'C') then
c
c  Plot the calibrated spectrum.
c
	      do j = 1, npts
	         dspec(j) = display.cspec(j)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Calibrated Calibration Spectrum'
	      else
	         plot_label = 'Calibrated Sky Spectrum'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'MJy/Sr',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 7  .or.  answer(1:2) .eq. 'AS') then
c
c  Plot the apodized voltage spectrum real-real and imaginary-imaginary
c	sigmas.
c
	      do j = 1, npts
	         dspec(j) =
     .		 cmplx(sqrt(display.vvar(1,j)),sqrt(display.vvar(2,j)))
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Apodized Voltage Calibration Spectrum Sigmas'
	      else
	         plot_label = 'Apodized Voltage Sky Spectrum Sigmas'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'Volts/Cm**2/Sr/Icm',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 8  .or.  answer(1:3) .eq. 'ACV') then
c
c  Plot the apodized voltage spectrum real-imaginary cross variances.
c
	      do j = 1, npts
	         dspec(j) = cmplx(display.vvar(3,j),0.0)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label =
     .		'Apodized Voltage Calibration Spectrum Cross Variances'
	      else
	         plot_label = 'Apodized Voltage Sky Spectrum Cross Variances'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', '(Volts/Cm**2/Sr/Icm)^2',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 9  .or.  answer(1:2) .eq. 'CS') then
c
c  Plot the calibrated spectrum real-real and imaginary-imaginary sigmas.
c
	      do j = 1, npts
	         dspec(j) =
     .		 cmplx(sqrt(display.cvar(1,j)),sqrt(display.cvar(2,j)))
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Calibrated Calibration Spectrum Sigmas'
	      else
	         plot_label = 'Calibrated Sky Spectrum Sigmas'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', 'MJy/Sr',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 10  .or.  answer(1:3) .eq. 'CCV') then
c
c  Plot the calibrated spectrum real-imaginary cross variances.
c
	      do j = 1, npts
	         dspec(j) = cmplx(display.cvar(3,j),0.0)
	      enddo

	      if (xcal_pos .eq. fac_xcalin) then
	         plot_label = 'Calibrated Calibration Spectrum Cross Variances'
	      else
	         plot_label = 'Calibrated Sky Spectrum Cross Variances'
	      endif

	      call fut_plot_title (label, fcc_ngroup, sweeps, fcc_speed,
     .				   fcc_length, fcc_chan, upmode, xcal_pos,
     .				   fakeit, gain, plot_label, title)
	      call fut_plot (dspec, npts, startx, spacex, title,
     .			     'Wave Number (icm)', '(MJy/Sr)^2',
     .			     zp, fac_present, plot_device,
     .			     plt_com, plt_com_file)


	   elseif (nans .eq. 11  .or.  answer(1:1) .eq. 'J') then
c
c  Jump over n spectra.
c
	      type 30
  30	      format (3x, 'Enter the number of spectra to jump ==> ', $)
	      accept *, fcc_jump
	      answer = 'Q'


	   elseif (nans .eq. 12  .or.  answer(1:1) .eq. 'Q') then
c
c  Disable plot flag.
c
	      fsl_display_spectra = fac_not_present
	      answer = 'Q'


	   elseif (answer(1:3) .eq. '   ') then
c
c  Go on to the next spectrum.
c
	      answer = 'Q'


	   else
c
c  Reply not recognized.
c
	      type 40
  40	      format (3x, 'Please try again.')

	   endif	!	(answer

	enddo		!	(while answer

	type 50
  50	format (/)


	return
	end
