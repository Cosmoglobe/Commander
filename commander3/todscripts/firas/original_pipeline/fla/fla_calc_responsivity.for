	integer * 4 function  fla_calc_responsivity ()

c-------------------------------------------------------------------------------
c
c	Function FLA_CALC_RESPONSIVITY
c
c	This function calculates the responsivity and time constant for the
c	Firas bolometers from the linear part of Steve Meyer's detector model.
c	The input parameters are the detector voltage, commanded bias, and
c	measured temperature for a given voltage spectrum and the derived
c	bolometer parameters from the FIRAS calibration model solution.  The
c	output parameters are the actual bolometer temperature, the detector
c	responsivity (in volts/erg), and the detector time constant for the
c	spectrum.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 1 June 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Include files:
c		fla_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fla_model)'

	integer * 4	j			!  a counter

	real	* 8	C			!  total heat capacity
	real	* 8	DT			!  detector temperature
						!    derivative
	real	* 8	G			!  total thermal conductance
	real	* 8	R			!  total detector resistance
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

	external	fla_normal

c
c  Extract the bolometer parameters from the array bolparm.
c
	R0   =  fex_model.bolparm(1)
	T0   =  fex_model.bolparm(2)
	G1   =  fex_model.bolparm(3)
	beta =  fex_model.bolparm(4)
	rho  =  fex_model.bolparm(5)
	C3   =  fex_model.bolparm(6)
	C1   =  fex_model.bolparm(7)
	Jo   =  fex_model.bolparm(8)
	Jg   =  fex_model.bolparm(9)
	RL   =  fac_load_resist

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
	do j = 1, 8
	   Tbol = T0 / dlog(Z*Tbol*dsinh(X/Tbol))**2
	enddo

c
c  Find the detector responsivity and time constant.
c
	X   =  Tbol/X * dtanh(X/Tbol)
	G   =  G1*Tbol**beta
	C   =  C3*Tbol**3 + C1*Tbol
	DT  =  1.0/X - 1.0 - 0.5*dsqrt(T0/Tbol)
	Z   =  (G*Tbol*R + DT*V**2) / (G*Tbol*R/X - DT*V**2)
	S0  =  fac_erg_to_watt * R * (Z-X) / (V * (Z*R/RL + 1.0) * (X+1.0))
	tau =  C/G * (Z+1.0) * (R*X + RL) / ((Z*R + RL) * (X+1.0))


	fla_calc_responsivity = %loc(fla_normal)

	return
	end
