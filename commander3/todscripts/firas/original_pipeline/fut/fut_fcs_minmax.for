	Integer*4  Function  FUT_FCS_Minmax  ( first, out_rec, in_rec )
c-----------------------------------------------------------------------------
c	Assign minimum and maximum values and time.
c	Larry P. Rosen, Hughes STX, March 1992
c-----------------------------------------------------------------------------
	Implicit None

c Passed
	Logical*1	first
	Dictionary	'FCS_SKY'
	Record / fcs_sky / out_rec
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec

c Function
	Logical*1	Time_LT

c Local
	Integer*4	m			! counter
c-----------------------------------------------------------------------------
c Begin

	FUT_FCS_Minmax = 0
	If (first) Then				! for min,max,1st,last values.
	   out_rec.coad_spec_head.first_gmt = in_rec.coad_spec_head.first_gmt
	   out_rec.coad_spec_head.first_time(1) =
	1     in_rec.coad_spec_head.first_time(1)
	   out_rec.coad_spec_head.first_time(2) =
	1     in_rec.coad_spec_head.first_time(2)
	   Do m = 1,6
	      out_rec.coad_spec_head.first_space_time(m) =
	1        in_rec.coad_spec_head.first_space_time(m)
	      out_rec.coad_spec_head.last_space_time(m) =
	1        in_rec.coad_spec_head.last_space_time(m)
	   Enddo
	   out_rec.coad_spec_head.last_gmt = in_rec.coad_spec_head.last_gmt
	   out_rec.coad_spec_head.last_time(1) =
	1     in_rec.coad_spec_head.last_time(1)
	   out_rec.coad_spec_head.last_time(2) =
	1     in_rec.coad_spec_head.last_time(2)
	   out_rec.coad_spec_head.first_mjr_frm_no =
	1     in_rec.coad_spec_head.first_mjr_frm_no
	   out_rec.coad_spec_head.last_mjr_frm_no =
	1     in_rec.coad_spec_head.last_mjr_frm_no
	   Do m = 1,10						! (2),(3)
	      out_rec.coad_spec_data.temp_min(m) =
	1        in_rec.coad_spec_data.temp_min(m)
	      out_rec.coad_spec_data.temp_max(m) =
	1        in_rec.coad_spec_data.temp_max(m)
	   Enddo
	   Do m = 1,4
	      out_rec.coad_spec_data.bol_volt_max(m) =
	1        in_rec.coad_spec_data.bol_volt_max(m)
	      out_rec.coad_spec_data.bol_volt_min(m) =
	1        in_rec.coad_spec_data.bol_volt_min(m)
	   Enddo
	   out_rec.attitude.terr_rad_byte = in_rec.attitude.terr_rad_byte
	   out_rec.coad_spec_data.dq_summary_flag =
	1     in_rec.coad_spec_data.dq_summary_flag
	   out_rec.coad_spec_data.att_summary_flag =
	1     in_rec.coad_spec_data.att_summary_flag
	Else
	   If ( Time_LT (in_rec.coad_spec_head.first_time,
	1                out_rec.coad_spec_head.first_time) ) Then
	      out_rec.coad_spec_head.first_gmt = in_rec.coad_spec_head.first_gmt
	      out_rec.coad_spec_head.first_time(1) =
	1        in_rec.coad_spec_head.first_time(1)
	      out_rec.coad_spec_head.first_time(2) =
	1        in_rec.coad_spec_head.first_time(2)
	      Do m = 1,6
	         out_rec.coad_spec_head.first_space_time(m) =
	1           in_rec.coad_spec_head.first_space_time(m)
	      Enddo
	   Endif
	   If ( Time_LT (out_rec.coad_spec_head.last_time,
	1                in_rec.coad_spec_head.last_time) ) Then
	      out_rec.coad_spec_head.last_gmt =
	1        in_rec.coad_spec_head.last_gmt
	      out_rec.coad_spec_head.last_time(1) =
	1        in_rec.coad_spec_head.last_time(1)
	      out_rec.coad_spec_head.last_time(2) =
	1        in_rec.coad_spec_head.last_time(2)
	      Do m = 1,6
	         out_rec.coad_spec_head.last_space_time(m) =
	1           in_rec.coad_spec_head.last_space_time(m)
	      Enddo
	   Endif
	   If (	in_rec.coad_spec_head.first_mjr_frm_no .LT.
	1       out_rec.coad_spec_head.first_mjr_frm_no) Then
	      out_rec.coad_spec_head.first_mjr_frm_no =
	1        in_rec.coad_spec_head.first_mjr_frm_no
	   Endif
	   If (	in_rec.coad_spec_head.last_mjr_frm_no .GT.
	1       out_rec.coad_spec_head.last_mjr_frm_no) Then
	      out_rec.coad_spec_head.last_mjr_frm_no =
	1        in_rec.coad_spec_head.last_mjr_frm_no
	   Endif
	   Do m = 1,10					! (2),(3)
	      out_rec.coad_spec_data.temp_min(m) =
	1        Min (out_rec.coad_spec_data.temp_min(m),
	2             in_rec.coad_spec_data.temp_min(m))
	      out_rec.coad_spec_data.temp_max(m) =
	1        Max (out_rec.coad_spec_data.temp_max(m),
	2             in_rec.coad_spec_data.temp_max(m))
	   Enddo
	   Do m = 1,4
	      out_rec.coad_spec_data.bol_volt_max(m) =
	1        Max (out_rec.coad_spec_data.bol_volt_max(m),
	2             in_rec.coad_spec_data.bol_volt_max(m))
	      out_rec.coad_spec_data.bol_volt_min(m) =
	1        Min (out_rec.coad_spec_data.bol_volt_min(m),
	2             in_rec.coad_spec_data.bol_volt_min(m))
	   Enddo
	   out_rec.attitude.terr_rad_byte = 
	1     Max (Zext (out_rec.attitude.terr_rad_byte),
	2          Zext (in_rec.attitude.terr_rad_byte ) )
	   out_rec.coad_spec_data.dq_summary_flag =
	1     Max (Zext (out_rec.coad_spec_data.dq_summary_flag),
	2          Zext (in_rec.coad_spec_data.dq_summary_flag ) )
	   out_rec.coad_spec_data.att_summary_flag =
	1     Max (Zext (out_rec.coad_spec_data.att_summary_flag),
	2          Zext (in_rec.coad_spec_data.att_summary_flag ) )

	Endif		! first; for min,max,1st,last values.
	Return
	End
