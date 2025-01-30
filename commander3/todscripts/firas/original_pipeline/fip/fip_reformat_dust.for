	integer * 4 function  fip_reformat_dust (in_rec, out_rec)

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_DUST
c
c	This function reformats the FIRAS dust spectra into the Project Dataset
c	format.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  5 October 1994
c
c-------------------------------------------------------------------------------
c
c	Input:
c		in_rec		FLA_DST			input record
c
c	Output:
c		out_rec		FIP_DST			output record
c
c	Subroutines called:
c		lib$movc5
c
c	Include files:
c		fip_frequency.txt
c		fip_invoc_dust.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c		Imaginary part of input spectrum was being written into real
c		part of output spectrum.  This was corrected.  F. Shuman, HSTX,
c		1995 Jan 23.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_dust)'
	include '(fip_frequency)'

	integer * 4	is			!  scan mode pointer
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	start			!  lower profile index

	real	* 4	imag_spectrum(180)	!  imaginary part of spectrum
	real	* 4	real_spectrum(180)	!  real part of spectrum

	dictionary 'fla_dst'
	record /fla_dst/ in_rec
	dictionary 'fip_dst'
	record /fip_dst/ out_rec

	external	fip_normal

c
c  Initialize the routine.
c
	call lib$movc5 (0,,0,720,real_spectrum)
	call lib$movc5 (0,,0,720,imag_spectrum)
	call lib$movc5 (0,,0,1488,out_rec)

c
c  Fill in the non-spectral fields.
c
	out_rec.pixel    = in_rec.pixel
	out_rec.eclon    = in_rec.eclon
	out_rec.eclat    = in_rec.eclat
	out_rec.num_ifgs = in_rec.num_ifgs
	out_rec.chanscan = in_rec.chanscan
	out_rec.nu_zero  = fcc_nu0
	out_rec.delta_nu = fcc_dnu
	out_rec.num_freq = fcc_nfreq
	out_rec.galon    = in_rec.galon
	out_rec.galat    = in_rec.galat
	out_rec.ra       = in_rec.ra
	out_rec.dec      = in_rec.dec

c
c  Read the spectrum fields from the input record.
c
	is = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length
	start = lolim(is)-1
	do j=lolim(is),uplim(is)
	   real_spectrum(j) = in_rec.real_spectrum(j-start)
	   imag_spectrum(j) = in_rec.imag_spectrum(j-start)
	enddo

c
c  Write the spectrum fields to the output record.
c	Convert from ergs/sec/cm^2/sr/icm to MJy/sr.
c
	do j = fcc_jlo,fcc_jhi
	   out_rec.real_spectrum(j-fcc_freq_offset) =
     .			real_spectrum(j) * fac_erg_to_mjy
	   out_rec.imag_spectrum(j-fcc_freq_offset) =
     .			imag_spectrum(j) * fac_erg_to_mjy
	enddo

	fip_reformat_dust = %loc(fip_normal)

	return
	end
