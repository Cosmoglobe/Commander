	Integer*4  Function  FCS_Spec  ( weight, fcs_rec, in_rec, sum,
	1                                rec1, first, dof )
c-----------------------------------------------------------------------------
c Function to store sum of weighted spectra, variance, bin vectors.
c Author: Larry P. Rosen, March 1992, Hughes STX
c Modified: Larry P. Rosen, February 1993, Hughes STX
c    Add high and low frequency cut-offs.  Determine index for low and high
c    frequencey cut-offs.  Set spectra and variances outside of the cut-off
c    indices to 0.
c Modified: Larry P. Rosen, October 1993.  Weight is now just number of ifgs.
c Modified: Larry P. Rosen, July 1994.  Weight is now glitch rate weighted.
c Modified: Larry P. Rosen, 11 August 1994. SER 11870.  FCS must copy
c    new attitude flag of record outside FDS galactic cut.  If first record
c    in pixel, copy flag.  If not, verify that the flag matches previous.
c-----------------------------------------------------------------------------
	Implicit None

	Include		'(fcs_msg)'
	Include		'(fut_fcs_include)'

c Passed Parameters

	Real*8		weight        ! glitchrate weighted nifgs
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec
	Dictionary	'FCS_SKY'
	Record / fcs_sky / fcs_rec
	Record / double_fcf / sum			! Record to store sums.
	Record / double_fcf / rec1	! Record to store first data record.
	Logical*1	first
	Real*8		dof	! degrees of freedom = weight - k

c local
	Integer*2	m	! counter
	Integer*2	k	! = 2 if 2nd template subtracted, 1 if NOT.
c-----------------------------------------------------------------------------
c Begin

	FCS_Spec = %loc (fcs_normal)

	k = in_rec.coad_spec_data.sec_template.subtracted + 1
        dof = weight - k

c Bin vector

	Do m = 1, 256
	   fcs_rec.coad_spec_data.bin_info(m) =
	1     fcs_rec.coad_spec_data.bin_info(m) +
	2     in_rec.coad_spec_data.bin_info(m)
	Enddo

c Sum weighted complex spectrum, sum variance.

	If (first) Then
	   first = .FALSE.
	   Do m = 1, 257
	      rec1.spec_data.spec(m) = in_rec.spec_data.spec(m)
	      sum.spec_data.real_var(m) = Dble (in_rec.spec_data.real_var(m))
	1        * Dble (weight * dof)
	      sum.spec_data.imag_var(m) = Dble (in_rec.spec_data.imag_var(m))
	1        * Dble (weight * dof)
	      sum.spec_data.real_imag_var(m) =
	1        Dble (in_rec.spec_data.real_imag_var(m))
	1        * Dble (weight * dof)
	   Enddo
	   fcs_rec.attitude.outside_galaxy_cut =
	1     in_rec.attitude.outside_galaxy_cut
	Else
	   Do m = 1, 257

c The following "If" is a fix to the vax stupid handling of floating point
c underflow.

	      If ( ((Abs (Real (in_rec.spec_data.spec(m))) .LE. 1.0e-30) .AND.
	1           (Abs (Dreal (rec1.spec_data.spec(m)))  .LE. 1.0e-30)) .OR.
	1          ((Abs (Aimag (in_rec.spec_data.spec(m))) .LE. 1.0e-30) .AND.
	1           (Abs (Dimag (rec1.spec_data.spec(m))) .LE. 1.0e-30)) ) Then
	         sum.spec_data.spec(m) =  (0.0,0.0)
	      Else
	         sum.spec_data.spec(m) = sum.spec_data.spec(m) + (Dcmplx
	1           (in_rec.spec_data.spec(m)) - rec1.spec_data.spec(m))
	1           * weight
              Endif
	      sum.spec_data.real_var(m) = sum.spec_data.real_var(m) +
	1        (Dble (in_rec.spec_data.real_var(m))
	1         * Dble (weight * dof))
	      sum.spec_data.imag_var(m) = sum.spec_data.imag_var(m) +
	1        (Dble (in_rec.spec_data.imag_var(m))
	1         * Dble (weight * dof))
	      sum.spec_data.real_imag_var(m) =
	1        sum.spec_data.real_imag_var(m) +
	1        (Dble (in_rec.spec_data.real_imag_var(m))
	1         * Dble (weight * dof))
	   Enddo
	   If ( fcs_rec.attitude.outside_galaxy_cut .NE.
	1       in_rec.attitude.outside_galaxy_cut ) Then
	      Call Lib$Signal ( fcs_galcut, %val(3),
	1                       in_rec.attitude.pixel_no,
	2                       fcs_rec.attitude.outside_galaxy_cut,
	3                       in_rec.attitude.outside_galaxy_cut )
	   Endif
	Endif
	Return
	End
