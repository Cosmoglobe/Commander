	integer * 4 function  fip_reformat_fef (in_rec, out_rec)

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_FEF
c
c	This function reformats the FIRAS FEF records into the Project Dataset
c	format.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  15 September 1994
c
c-------------------------------------------------------------------------------
c
c	Input:
c		in_rec		FEF_SPC			input record
c
c	Output:
c		out_rec		FIP_FEF			output record
c
c	Subroutines called:
c		lib$movc5
c
c	Include files:
c		fip_frequency.txt
c		fip_invoc_fef.txt
c		fut_params.txt
c		uoe_constants.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_fef)'
	include '(fip_frequency)'
	include '(uoe_constants)'

	complex * 8	diff_spec(180)		!  difference spectrum

	integer * 4	atime(2)		!  ADT time
	integer * 4	btime(2)		!  ADT time
	integer * 4	ctime(2)		!  ADT time
	integer * 4	dtime(2)		!  ADT time
	integer * 4	is			!  scan mode pointer
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	start			!  lower profile index

	real	* 4	diff_spec_sig(180)	!  difference spectrum sigmas
	real	* 4	ratio_spec(180)		!  ratio spectrum
	real	* 4	ratio_spec_sig(180)	!  ratio spectrum sigmas

	real	* 8	uoe_adt2t68

	dictionary 'fef_spc'
	record /fef_spc/ in_rec
	dictionary 'fip_fef'
	record /fip_fef/ out_rec

	external	fip_normal

c
c  Initialize the routine.
c
	call lib$movc5 (0,,0,1440,diff_spec)
	call lib$movc5 (0,,0,720,diff_spec_sig)
	call lib$movc5 (0,,0,720,ratio_spec)
	call lib$movc5 (0,,0,720,ratio_spec_sig)
	call lib$movc5 (0,,0,3844,out_rec)

c
c  Fill in the identification fields.
c
	out_rec.descrip   = in_rec.descrip
	out_rec.chanscan1 = in_rec.chanscan1
	out_rec.chanscan2 = in_rec.chanscan2
	do k=1,4
	   do j=1,2
	      atime(j) = in_rec.start_times1(j,k)
	      btime(j) = in_rec.stop_times1(j,k)
	      ctime(j) = in_rec.start_times2(j,k)
	      dtime(j) = in_rec.stop_times2(j,k)
	   enddo
	   out_rec.start_times1(k) = uoe_adt2t68(atime) - T68IRASBase
	   out_rec.stop_times1(k)  = uoe_adt2t68(btime) - T68IRASBase
	   out_rec.start_times2(k) = uoe_adt2t68(ctime) - T68IRASBase
	   out_rec.stop_times2(k)  = uoe_adt2t68(dtime) - T68IRASBase
	enddo
	do j=1,2
	   out_rec.tau1(j) = in_rec.tau1(j)
	   out_rec.tau2(j) = in_rec.tau2(j)
	enddo
	out_rec.nu_zero   = fcc_nu0
	out_rec.delta_nu  = fcc_dnu
	out_rec.num_freq  = fcc_nfreq

c
c  Read the spectrum fields from the input record.
c
	is = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length
	start = lolim(is)-1
	do j=lolim(is),uplim(is)
	   diff_spec(j)      = in_rec.diff_spec(j-start)
	   diff_spec_sig(j)  = in_rec.diff_spec_sig(j-start)
	   ratio_spec(j)     = in_rec.ratio_spec(j-start)
	   ratio_spec_sig(j) = in_rec.ratio_spec_sig(j-start)
	enddo

c
c  Write the spectrum fields to the output record.
c	Convert from ergs/sec/cm^2/sr/icm to MJy/sr.
c
	do j = fcc_jlo,fcc_jhi
	   out_rec.rdiff_spec(j-fcc_freq_offset)     = 
     .			real(diff_spec(j)) * fac_erg_to_mjy
	   out_rec.idiff_spec(j-fcc_freq_offset)     = 
     .			aimag(diff_spec(j)) * fac_erg_to_mjy
	   out_rec.diff_spec_sig(j-fcc_freq_offset)  =
     .			diff_spec_sig(j) * fac_erg_to_mjy
	   out_rec.ratio_spec(j-fcc_freq_offset)     = ratio_spec(j)
	   out_rec.ratio_spec_sig(j-fcc_freq_offset) = ratio_spec_sig(j)
	enddo

	fip_reformat_fef = %loc(fip_normal)

	return
	end
