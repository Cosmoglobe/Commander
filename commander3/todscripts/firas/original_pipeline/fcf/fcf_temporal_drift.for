	integer * 4 function  fcf_temporal_drift (spec_time, ticald,
     .						  primary_vib)

c-------------------------------------------------------------------------------
c
c	Function FCF_TEMPORAL_DRIFT
c
c	This function computes the time in years sinc the aperture cover
c	ejection for a spectrum.  It then calculates the ical temporal drift
c	correction and the primary (time_dependent) vibration correction
c	for the spectrum.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  26 June 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		spec_time	integer * 4		VAX ADT time of spectrum
c
c	Output:
c		primary_vib	real * 8		primary vibration
c								correction
c		ticald		real * 8		Ical temporal drift
c								correction
c
c	Include files:
c		fcf_model.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants:
c
c		APCO_DATE	9626073 sec		time of aperture cover
c							   ejection in VAX ADT
c							   time
c
c		YEAR_LEN	73426.0D+00  sec	length of year VAX ADT
c							    time
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fcf_model)'

	integer	* 4	apco_date		!  aperature cover ejection
	parameter	(apco_date = 9626073)   !    date in VAX ADT time
	integer * 4	j			!  a counter
	integer * 4	spec_time		!  high-order word of VAX ADT
						!    time of spectrum

	real	* 8	primary_vib		!  primary (time-dependent)
						!    vibration correction
	real	* 8	ticald			!  Ical temperature drift
						!    correction
	real	* 8	time			!  time since aperture cover
						!    ejection in years
	real	* 8	year_len		!  length of a year in VAX ADT
	parameter	(year_len = 73426.0D0)	!    time

	external	fcf_normal


c
c  Calculate the time since aperature cover ejection for the spectrum.
c
	time = dble(spec_time - apco_date) / year_len

c
c  Calculate the Ical temporal drift correction.
c
	ticald = param(12) * dexp(-param(13)*time) - param(14)

c
c  Calculate the primary (time-dependent) vibration correction.
c
	primary_vib = 0.0D0
	do j = 4,0,-1
	   primary_vib = primary_vib * time + param(16+j)
	enddo


	fcf_temporal_drift = %loc(fcf_normal)

	return
	end
