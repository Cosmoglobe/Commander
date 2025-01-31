	Integer*4  Function  FCS_Minmax  (first, fcs_rec, in_rec)
c-----------------------------------------------------------------------------
c	Assign minimum and maximum values and time.
c	Larry P. Rosen, Hughes STX, March 1992
c-----------------------------------------------------------------------------
	Implicit None

	Include		'(fcs_msg)'

c Passed
	Logical*1	first
	Dictionary	'FCS_SKY'
	Record / fcs_sky / fcs_rec
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec

c Function
	Logical*1	Time_LT

c Local
	Integer*4	m			! counter
c-----------------------------------------------------------------------------
c Begin

	FCS_Minmax = %loc (fcs_normal)
	If (first) Then				! for min,max,1st,last values.
	   fcs_rec.coad_spec_head.first_gmt = in_rec.coad_spec_head.first_gmt
	   fcs_rec.coad_spec_head.first_time(1) =
	1     in_rec.coad_spec_head.first_time(1)
	   fcs_rec.coad_spec_head.first_time(2) =
	1     in_rec.coad_spec_head.first_time(2)
	   Do m = 1,6
	      fcs_rec.coad_spec_head.first_space_time(m) =
	1        in_rec.coad_spec_head.first_space_time(m)
	      fcs_rec.coad_spec_head.last_space_time(m) =
	1        in_rec.coad_spec_head.last_space_time(m)
	   Enddo
	   fcs_rec.coad_spec_head.last_gmt = in_rec.coad_spec_head.last_gmt
	   fcs_rec.coad_spec_head.last_time(1) =
	1     in_rec.coad_spec_head.last_time(1)
	   fcs_rec.coad_spec_head.last_time(2) =
	1     in_rec.coad_spec_head.last_time(2)
	   fcs_rec.coad_spec_head.first_mjr_frm_no =
	1     in_rec.coad_spec_head.first_mjr_frm_no
	   fcs_rec.coad_spec_head.last_mjr_frm_no =
	1     in_rec.coad_spec_head.last_mjr_frm_no
	   Do m = 1,10						! (2),(3)
	      fcs_rec.coad_spec_data.temp_min(m) =
	1        in_rec.coad_spec_data.temp_min(m)
	      fcs_rec.coad_spec_data.temp_max(m) =
	1        in_rec.coad_spec_data.temp_max(m)
	   Enddo
	   Do m = 1,4
	      fcs_rec.coad_spec_data.bol_volt_max(m) =
	1        in_rec.coad_spec_data.bol_volt_max(m)
	      fcs_rec.coad_spec_data.bol_volt_min(m) =
	1        in_rec.coad_spec_data.bol_volt_min(m)
	   Enddo
	   fcs_rec.attitude.terr_rad_byte = in_rec.attitude.terr_rad_byte
	   fcs_rec.coad_spec_data.dq_summary_flag =
	1     in_rec.coad_spec_data.dq_summary_flag
	   fcs_rec.coad_spec_data.att_summary_flag =
	1     in_rec.coad_spec_data.att_summary_flag
	Else
	   If ( Time_LT (in_rec.coad_spec_head.first_time,
	1                fcs_rec.coad_spec_head.first_time) ) Then
	      fcs_rec.coad_spec_head.first_gmt = in_rec.coad_spec_head.first_gmt
	      fcs_rec.coad_spec_head.first_time(1) =
	1        in_rec.coad_spec_head.first_time(1)
	      fcs_rec.coad_spec_head.first_time(2) =
	1        in_rec.coad_spec_head.first_time(2)
	      Do m = 1,6
	         fcs_rec.coad_spec_head.first_space_time(m) =
	1           in_rec.coad_spec_head.first_space_time(m)
	      Enddo
	   Endif
	   If ( Time_LT (fcs_rec.coad_spec_head.last_time,
	1                in_rec.coad_spec_head.last_time) ) Then
	      fcs_rec.coad_spec_head.last_gmt =
	1        in_rec.coad_spec_head.last_gmt
	      fcs_rec.coad_spec_head.last_time(1) =
	1        in_rec.coad_spec_head.last_time(1)
	      fcs_rec.coad_spec_head.last_time(2) =
	1        in_rec.coad_spec_head.last_time(2)
	      Do m = 1,6
	         fcs_rec.coad_spec_head.last_space_time(m) =
	1           in_rec.coad_spec_head.last_space_time(m)
	      Enddo
	   Endif
	   If (	in_rec.coad_spec_head.first_mjr_frm_no .LT.
	1       fcs_rec.coad_spec_head.first_mjr_frm_no) Then
	      fcs_rec.coad_spec_head.first_mjr_frm_no =
	1        in_rec.coad_spec_head.first_mjr_frm_no
	   Endif
	   If (	in_rec.coad_spec_head.last_mjr_frm_no .GT.
	1       fcs_rec.coad_spec_head.last_mjr_frm_no) Then
	      fcs_rec.coad_spec_head.last_mjr_frm_no =
	1        in_rec.coad_spec_head.last_mjr_frm_no
	   Endif
	   Do m = 1,10					! (2),(3)
	      fcs_rec.coad_spec_data.temp_min(m) =
	1        Min (fcs_rec.coad_spec_data.temp_min(m),
	2             in_rec.coad_spec_data.temp_min(m))
	      fcs_rec.coad_spec_data.temp_max(m) =
	1        Max (fcs_rec.coad_spec_data.temp_max(m),
	2             in_rec.coad_spec_data.temp_max(m))
	   Enddo
	   Do m = 1,4
	      fcs_rec.coad_spec_data.bol_volt_max(m) =
	1        Max (fcs_rec.coad_spec_data.bol_volt_max(m),
	2             in_rec.coad_spec_data.bol_volt_max(m))
	      fcs_rec.coad_spec_data.bol_volt_min(m) =
	1        Min (fcs_rec.coad_spec_data.bol_volt_min(m),
	2             in_rec.coad_spec_data.bol_volt_min(m))
	   Enddo
	   fcs_rec.attitude.terr_rad_byte = 
	1     Max (Zext (fcs_rec.attitude.terr_rad_byte),
	2          Zext (in_rec.attitude.terr_rad_byte ) )
	   fcs_rec.coad_spec_data.dq_summary_flag =
	1     Max (Zext (fcs_rec.coad_spec_data.dq_summary_flag),
	2          Zext (in_rec.coad_spec_data.dq_summary_flag ) )
	   fcs_rec.coad_spec_data.att_summary_flag =
	1     Max (Zext (fcs_rec.coad_spec_data.att_summary_flag),
	2          Zext (in_rec.coad_spec_data.att_summary_flag ) )

	Endif		! first; for min,max,1st,last values.
	Return
	End
