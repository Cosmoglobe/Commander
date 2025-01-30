	Integer*4 Function FTB_Fill_TempC(data, num_recs, num_pts, time_flags,
	2                                 hkp, max, data_present, side)

C-----------------------------------------------------------------------------
C
C      Purpose:
C
C          To evaluate the temperature controller status bit fields and
C          place them in the R*4 array for PLT.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C---------------------------------------------------------------------------
C	Revisions:
C
C	   Version 4.4  1989 May 10, SPR 3045.  Fred Shuman, STX.
C			New facility BOZO to be brought into CSDR standards
C			and renamed FTB_TempC.
C
C	   Version 4.6  1989 Nov 21, SPR 5023.  Fred Shuman, STX.
C		Should not require the integrators to be on to plot.
C
C-----------------------------------------------------------------------------

	Implicit None

	Integer*4 num_recs	      ! Number of records possible in time range
	Real   *4 data(num_recs*2,1)  ! hkp bit data - output matrix
	Integer*4 num_pts	      ! Number of plot points in time range
	Integer*4 time_flags(1)       ! "Flags" indicating   \   0 = too short
C				      !  quality of the step  >  1 = OK
C				      !  from previous time  /   2 = too long
	  Dictionary 'nfs_hkp'	   
	  Record/nfs_hkp/ hkp(1)   
	Real   *4 max		   ! Max time - output parameter
	Logical*4 data_present
	Character*1 side

	Integer*4 i, fr, j, k
	Integer*4 LIB$SubX
	Integer*4 LIB$EDiv
	Integer*4 diff_time(2)
	Integer*4 diff_mj(2)
	Integer*4 diff_rec
	Integer*4 sec60/600000000/
	Integer*4 quo
	Integer*4 rem
	Integer*4 status
	Integer*4 on
	Integer*4 enable
	Integer*4 mjfm
	Integer*4 stat_mon
	Integer*4 start_bit
	Integer*4 mf_a1
	Integer*4 mf_a2
	Integer*4 mf_b1
	Integer*4 mf_b2
	Integer*4 bin_time(2)
	Integer*4 cli$get_value
	Integer*4 tc4a
	Integer*4 tc2d
	Integer*4 tcsl
	Integer*4 telem_i4_temp
	Integer*4 telem_i4_gain
	Integer*4 vax_i4_temp
	Integer*4 vax_i4_cur
	Integer*4 vax_i4_int
	Integer*4 vax_i4_pro
	Integer*4 LIB$MovC3
	Integer*4 CT_GMT_to_Binary

	Real   *4 no_data/-1.2E-34/

	Logical*4 goodrec		! Flag -- Are the 2 MFs of this record
C					!   different FIRAS MF numbers?

	Common /data_info/ diff_mj, diff_rec
	Common /tc_info/ on, enable, mjfm, stat_mon, start_bit

	FTB_Fill_TempC = 1
	data_present = .False.

	status = cli$get_value('side', side,)

	num_pts = 1
	Do i = 1,num_recs

C   If this time is flagged for following the previous time by too much, send
C     'no data' values to PLT to create a gap in the plot before this point.

	  If (time_flags(i) .Eq. 2) Then
	    data(num_pts, 1) = data(num_pts-1, 1)

	    Do j = 2,5
	      data(num_pts, j) = no_data
	    End Do

	    num_pts = num_pts + 1
	  End If

C   If data present, get it.

	  If (time_flags(i) .Eq. 1 .Or. time_flags(i) .Eq. 2) Then

	    If (side .Eq. 'A') Then

C   Temperature controllers, FIRAS side A --
C   ..Is the first S/C major frame FIRAS maj frm 2 or 1 ?
C
	      If (btest(hkp(i).frame(1).hskp_head.stat_monitor_cmd(1), 14)) Then
	        mf_a1 = 2
	      Else
	        mf_a1 = 1
	      End If
C
C   ..Is the 2nd S/C major frame FIRAS maj frm 2 or 1 ?
C
	      If (btest(hkp(i).frame(2).hskp_head.stat_monitor_cmd(1), 14)) Then
	        mf_a2 = 2
	      Else
	        mf_a2 = 1
	      End If
C
C.....This concludes this test of which major frame is what on side A.
C
	      If (mf_a1 .Eq. mf_a2) Then

	        goodrec = .False.

	      Else If (mf_a1 .Eq. mjfm) Then

	        goodrec = .True.
	        tc4a = hkp(i).frame(1).hskp_head.ipdu_stat(3)
	        tc2d = hkp(i).frame(1).hskp_head.ipdu_stat(5)
	        tcsl = hkp(i).frame(1).hskp_head.stat_monitor_cmd(1)

	        status = LIB$MovC3(8, hkp(i).ct_head.time, bin_time)
	        telem_i4_temp =
	2                   hkp(i).frame(1).hskp_head.stat_monitor_cmd(stat_mon)
	        telem_i4_gain =
	2                   hkp(i).frame(1).hskp_head.stat_monitor_cmd(4)

	      Else

	        goodrec = .True.
	        tc4a = hkp(i).frame(2).hskp_head.ipdu_stat(3)
	        tc2d = hkp(i).frame(2).hskp_head.ipdu_stat(5)
	        tcsl = hkp(i).frame(2).hskp_head.stat_monitor_cmd(1)

	        status = CT_GMT_to_Binary(hkp(i).hskp_tail.gmt_mjf2, bin_time)
	        telem_i4_temp =
	2                   hkp(i).frame(2).hskp_head.stat_monitor_cmd(stat_mon)
	        telem_i4_gain =
	2                   hkp(i).frame(2).hskp_head.stat_monitor_cmd(4)

	      End If      ! (mf_a1 .Eq. mf_a2)

	    Else

C   Temperature controllers, FIRAS side B --
C   ..Is the first S/C major frame FIRAS maj frm 2 or 1 ?
C
	      If (btest(hkp(i).frame(1).hskp_head.stat_monitor_cmd(5), 14)) Then
	        mf_b1 = 2
	      Else
	        mf_b1 = 1
	      End If
C
C   ..Is the 2nd S/C major frame FIRAS maj frm 2 or 1 ?
C
	      If (btest(hkp(i).frame(2).hskp_head.stat_monitor_cmd(5), 14)) Then
	        mf_b2 = 2
	      Else
	        mf_b2 = 1
	      End If
C
C.....This concludes this test of which major frame is what on side B.
C
	      If (mf_b1 .Eq. mf_b2) Then

	        goodrec = .False.

	      Else If (mf_b1 .Eq. mjfm) Then

	        goodrec = .True.
	        tc4a = hkp(i).frame(1).hskp_head.ipdu_stat(4)
	        tc2d = hkp(i).frame(1).hskp_head.ipdu_stat(6)
	        tcsl = hkp(i).frame(1).hskp_head.stat_monitor_cmd(5)

	        status = LIB$MovC3(8, hkp(i).ct_head.time, bin_time)
	        telem_i4_temp =
	2                 hkp(i).frame(1).hskp_head.stat_monitor_cmd(stat_mon+4)
	        telem_i4_gain =
	2                   hkp(i).frame(1).hskp_head.stat_monitor_cmd(8)

	      Else

	        goodrec = .True.
	        tc4a = hkp(i).frame(2).hskp_head.ipdu_stat(4)
	        tc2d = hkp(i).frame(2).hskp_head.ipdu_stat(6)
	        tcsl = hkp(i).frame(2).hskp_head.stat_monitor_cmd(5)

	        status = CT_GMT_to_Binary(hkp(i).hskp_tail.gmt_mjf2, bin_time)
	        telem_i4_temp =
	2                 hkp(i).frame(2).hskp_head.stat_monitor_cmd(stat_mon+4)
	        telem_i4_gain =
	2                   hkp(i).frame(2).hskp_head.stat_monitor_cmd(8)

	      End If      ! (mf_b1 .Eq. mf_b2)

	    End If        ! (side .Eq. 'A')

	    status = LIB$SubX(hkp(i).ct_head.time, hkp(1).ct_head.time,
	2                          diff_time)
	    status = LIB$EDiv(sec60, diff_time, quo, rem)
	    data(num_pts, 1) = Floatj(quo)/60. + Floatj(rem)/(sec60*60.)

C   Note that VAX bit position is reversed from telemetry position

	    If (bjtest(tc4a, 7-on) .And. bjtest(tc2d, 7-on) .And. goodrec) Then

	      data_present = .True.

	      vax_i4_temp = 0
	      Call mvbits(telem_i4_temp, 0, 12, vax_i4_temp, 0)
	      data(num_pts, 2) = Floatj(vax_i4_temp)

	      vax_i4_cur = 0
	      k = start_bit
	      Call mvbits(telem_i4_gain, 14-k, 2, vax_i4_cur, 0)

	      vax_i4_int = 0
	      vax_i4_pro = 0
	      Call mvbits(telem_i4_gain, 9-3*k, 3, vax_i4_int, 0)
	      Call mvbits(telem_i4_gain, 6-3*k, 3, vax_i4_pro, 0)

	      data(num_pts, 3) = Floatj(vax_i4_cur)
	      data(num_pts, 4) = Floatj(vax_i4_int)
	      data(num_pts, 5) = Floatj(vax_i4_pro)

C   Tempcontrollers are not on, so send 'no data' values to PLT.

	    Else

	      Do j = 2,5
	        data(num_pts, j) = no_data
	      End Do

	    End If    ! ( bjtest(tc4a, 7-on) .And. bjtest(tc2d, 7-on) )

	  Else If (time_flags(i) .Eq. 0) Then

C   Else, this time is flagged for following the previous time by too little,
C     so send 'no data' values to PLT to suppress plotting it and make a gap.

	    data(num_pts, 1) = data(num_pts-1, 1)
	    Do j = 2,5
	      data(num_pts, j) = no_data
	    End Do

	  End If      ! (time_flags(i) = 1 or 2)

	  num_pts = num_pts + 1

	End Do        ! i = 1,num_recs

	num_pts = num_pts - 1
	max = data(num_pts, 1)

	Return
	End
