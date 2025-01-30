	Integer*4  Function  FCS_Non_Avg  ( fcs_rec, in_rec, nscan, scan,
	1                                   scan2, input )
c-----------------------------------------------------------------------------
c Function to store non-averaged data in output record.
c Author: Larry P. Rosen, Hughes STX, March 1992
c Note:  If two scan modes are being combined, then the speed/length fields
c are filled in according to the rule: if the fields differ, a -1 is put in,
c otherwise just the correct value is put into the fcs record.
c-----------------------------------------------------------------------------
	Implicit None

c Include

	Include		'(fcs_msg)'
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec
	Dictionary	'FCS_SKY'
	Record / fcs_sky / fcs_rec
	Integer*2	nscan				! Number scan modes
	Character*2	scan, scan2			! scan modes
	Character*3	input			! Input data type; FAD or FCF
	Integer*4	m

c Begin

	FCS_Non_Avg = %loc (fcs_normal)
	fcs_rec.spec_data.model_ttag = in_rec.spec_data.model_ttag
	fcs_rec.spec_data.model_label = in_rec.spec_data.model_label
	fcs_rec.coad_spec_data.chan_id = in_rec.coad_spec_data.chan_id
	If (nscan .EQ. 1) Then
	   fcs_rec.coad_spec_data.mtm_speed = in_rec.coad_spec_data.mtm_speed
	   fcs_rec.coad_spec_data.mtm_length = in_rec.coad_spec_data.mtm_length
	Else
	   If (scan(1:1) .EQ. scan2(1:1)) Then
	      fcs_rec.coad_spec_data.mtm_length=in_rec.coad_spec_data.mtm_length
	   Else
	      fcs_rec.coad_spec_data.mtm_length = -1
	   Endif
	   If (scan(2:2) .eq. scan2(2:2)) then
	      fcs_rec.coad_spec_data.mtm_speed = in_rec.coad_spec_data.mtm_speed
	   Else
	      fcs_rec.coad_spec_data.mtm_speed = -1
	   Endif
	Endif
	fcs_rec.coad_spec_data.fakeit = in_rec.coad_spec_data.fakeit
	fcs_rec.coad_spec_data.sci_mode = in_rec.coad_spec_data.sci_mode
	fcs_rec.coad_spec_data.xcal_pos = in_rec.coad_spec_data.xcal_pos
	fcs_rec.attitude.terr_rad_byte = -1
	If (input .EQ. 'FAD') Then
	   fcs_rec.spec_data.destriped = 1
	Endif
	Return
	End
