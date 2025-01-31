	integer * 4 function  fsl_calc_responsivity (bol_volt, cmd_bias, Tdet,
     .						     Tbol, Qrad)

c-------------------------------------------------------------------------------
c
c	Function FSL_CALC_RESPONSIVITY
c
c	This function calculates the responsivity and time constant for the
c	Firas bolometers from the linear part of Steve Meyers detector model.
c	The input parameters are the detector voltage, commanded bias, and
c	measured temperature for a given voltage spectrum and the derived
c	bolometer parameters from the FIRAS calibration model solution.  The
c	output parameters are the actual bolometer temperature, the IR power
c	incident on the bolometer (in Watts), the detector responsivity
c	(in volts/erg), and the detector time constant for the spectrum.
c
c	Author:	 
c                FCF_Calc_Responsivity
c                Gene Eplee
c		 General Sciences Corp.
c		 10 March 1993
c
c                FSL_Calc_Responsivity
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		bol_volt	real * 8		detector voltage
c		cmd_bias	real * 8		commanded bolometer bias
c		Tdet		real * 8		measured detector
c							temperature
c
c	Output:
c		Tbol		real * 8		actual detector 
c							temperature
c		Qrad		real * 8		IR power incident on
c							bolometer
c
c	Subroutines called:
c		lib$signal
c
c	Include files:
c		fsl_display.txt
c		fsl_model.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants:
c
c		DELTA_T		1.0D-06  K		bolometer temperature
c							   convergence criterion
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 9, 1995 
c       Modified FCF_Calc_Responsivity to FSL_Calc_Responsivity for the new 
c       FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.Changed status, include file, and function names 
c       for FSL.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fsl_display)'
	include '(fsl_model)'

	integer * 4	j			!  a counter

	real	* 8	bol_volt		!  detector voltage
	real	* 8	BP			!  bias power on the bolometer
	real	* 8	C			!  total heat capacity
	real	* 8	cmd_bias		!  commanded bolometer bias
	real	* 8	delta_T			!  detector temperature
	parameter	(delta_t = 1.0D-06)	!    convergence tolerance
	real	* 8	DT			!  detector temperature
						!    derivative
	real	* 8	G			!  total thermal conductance
	real	* 8	P			!  power on the bolometer
	real	* 8	Qrad			!  IR power on the bolometer
	real	* 8	R			!  total detector resistance
	real	* 8	Tbol			!  actual detector temperature
	real	* 8	Tconv			!  detector temperature for
						!    convergence check
	real	* 8	Tdet			!  measured detector temperature
	real	* 8	V			!  total detector bias
	real	* 8	X			!  non-ideal electric field term
	real	* 8	Z			!  detector impedence

c
c  Bolometer parameters
c
	real	* 8	R0			!  detector resistance at
						!    infinite temperature
	real	* 8	T0			!  characteristic temperature
						!    for detector resistance
						!    function
	real	* 8	G1			!  coefficient of detector
						!    thermal conductance
	real	* 8	beta			!  index of temperature 
						!    dependence of detector
						!    thermal conductance
	real	* 8	rho			!  electric field dependence  of
						!    detector resistance
	real	* 8	C3			!  coefficient of cubic heat
						!    capacity term
	real	* 8	C1			!  coefficient of linear heat
						!    capacity term
	real	* 8	Jo			!  JFET offset
	real	* 8	Jg			!  JFET gain
	real	* 8	RL			!  detector load resistance

	external	fsl_dettempconv
	external	fsl_normal

c
c  Extract the bolometer parameters from the array bolparm.
c
	R0   =  bolparm(1)
	T0   =  bolparm(2)
	G1   =  bolparm(3)
	beta =  bolparm(4)
	rho  =  bolparm(5)
	C3   =  bolparm(6)
	C1   =  bolparm(7)
	Jo   =  bolparm(8)
	Jg   =  bolparm(9)
	RL   =  bolparm(10)

c
c  Find the bolometer state.
c
	V    =  (bol_volt - Jo) / Jg
	R    =  RL*V / (cmd_bias - V)
	X    =  V*rho
	Z    =  R/R0/X
	Tbol =  Tdet

c
c  Iterate to find the actual detector temperature.
c
	do j = 1, 7
	   Tbol = T0 / dlog(Z*Tbol*dsinh(X/Tbol))**2
	enddo
	Tconv = Tbol
	Tbol = T0 / dlog(Z*Tbol*dsinh(X/Tbol))**2
	if (dabs(Tbol-Tconv) .ge. delta_t) then
	   fsl_calc_responsivity = %loc(fsl_dettempconv)
	   call lib$signal (fsl_dettempconv, %val(2), display.gmt,
     .					     %val(display.pixel_no))
	   return
	endif

c
c  Find the IR power incident on the detector.
c
	BP = V * (cmd_bias-V) / RL
	P = (G1/(beta+1.0)) * (Tbol**(beta+1.0) - Tdet**(beta+1.0))
	Qrad = P - BP

c
c  Find the detector responsivity and time constant.
c
	X   =  Tbol/X * dtanh(X/Tbol)
	G   =  G1*Tbol**beta
	C   =  C3*Tbol**3 + C1*Tbol
	DT  =  1.0/X - 1.0 - 0.5*dsqrt(T0/Tbol)
	Z   =  (G*Tbol*R + DT*V**2) / (G*Tbol*R/X - DT*V**2)
	S0  =  rscale * R * (Z-X) / (V * (Z*R/RL + 1.0) * (X+1.0))
	tau =  C/G * (Z+1.0) * (R*X + RL) / ((Z*R + RL) * (X+1.0))


	fsl_calc_responsivity = %loc(fsl_normal)

	return
	end
