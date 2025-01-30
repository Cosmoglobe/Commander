	Integer*4  Function  FUT_FCS_Non_Avg  ( out_rec, in_rec, scan, scan2,
	1                                       chan, chan2, input )
c-----------------------------------------------------------------------------
c  Purpose:  Function to store non-averaged data of an FCS type record in the
c     output record.
c
c  Input:
c     in_rec - an FAD, FCF, FCS or FMS type record
c     scan    - Char*2 scan mode of skymap
c     scan2   - Char*2 scan mode of second skymap
c     chan    - Char*2 channel of skymap
c     chan2   - Char*2 channel of second skymap
c     input   - Char*3 input data type 'FAD', 'FCF', 'FCS', or 'FMS'
c  Output:
c     out_rec - an FCS or FMS type record
c
c  Author: Larry P. Rosen, Hughes STX, March 1992
c
c  Modified: LPR, October 1993.
c     Routine changed from FCS_Non_Avg to FUT_FCS_Non_Avg, an FUT routine to
c     be used by both FCS combining FCF or FAD records, and by FMS combining
c     FCS type records in various states of combination.
c  Modified: LPR, March 1994.
c     Pixel, Solution, Pixel_Definition, and Skymap_Index are now set here.
c
c  Notes:
c     When FCS is run, chan = chan2 and scan = scan2.
c     If the two channels are the same, then the output CHAN_ID is the same
c     as the input.  If the two scan modes are the same then the output
c     MTM_SPEED and MTM_LENGTH are the same as the input.
c     If the channels are different, then the output CHAN_ID will be -1.
c     If the speeds are different, then the output MTM_SPEED will be -1.
c     If the lengths are different, then the output MTM_LENGTH will be -1.
c-----------------------------------------------------------------------------
	Implicit None

c Include

	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec
	Dictionary	'FCS_SKY'
	Record / fcs_sky / out_rec
	Integer*2	nscan				! Number scan modes
	Character*2	scan, scan2			! scan modes
	Character*2	chan, chan2			! channels
	Character*3	input				! Input data type: FAD,
							!   FCF, FCS, or FMS
c Begin

	FUT_FCS_Non_Avg = 0
	out_rec.attitude.pixel_no = in_rec.attitude.pixel_no
	out_rec.attitude.solution = in_rec.attitude.solution
	out_rec.attitude.pixel_definition = in_rec.attitude.pixel_definition
	out_rec.attitude.skymap_index = in_rec.attitude.skymap_index
	out_rec.spec_data.model_ttag = in_rec.spec_data.model_ttag
	out_rec.spec_data.model_label = in_rec.spec_data.model_label

	If (chan .EQ. chan2) Then
	   out_rec.coad_spec_data.chan_id = in_rec.coad_spec_data.chan_id
	Else           ! Channels are different; the output CHAN_ID will be -1.
	   out_rec.coad_spec_data.chan_id = -1
	Endif
	If (scan(1:1) .EQ. scan2(1:1)) Then
	   out_rec.coad_spec_data.mtm_length = in_rec.coad_spec_data.mtm_length
	Else
	   out_rec.coad_spec_data.mtm_length = -1
	Endif
	If (scan(2:2) .EQ. scan2(2:2)) then
	   out_rec.coad_spec_data.mtm_speed = in_rec.coad_spec_data.mtm_speed
	Else
	   out_rec.coad_spec_data.mtm_speed = -1
	Endif
	out_rec.coad_spec_data.fakeit = in_rec.coad_spec_data.fakeit
	out_rec.coad_spec_data.sci_mode = in_rec.coad_spec_data.sci_mode
	out_rec.coad_spec_data.xcal_pos = in_rec.coad_spec_data.xcal_pos
	out_rec.attitude.terr_rad_byte = -1

c Set destripe flag to input.

	If (input .EQ. 'FAD') Then
	   out_rec.spec_data.destriped = 1
	Else
	   out_rec.spec_data.destriped = in_rec.spec_data.destriped
	Endif
	Return
	End
