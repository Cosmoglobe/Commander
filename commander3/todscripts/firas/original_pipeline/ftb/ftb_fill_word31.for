	Integer*4 Function FTB_Fill_Word31(data, step, num_recs, num_pts, bits,
	2                        num_fields, mult_fac, time_flags, word31, max)

C-----------------------------------------------------------------------------
C
C      Purpose:
C
C          To evaluate the bit fields and place them in the R*4 array for
C          PLT.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C-----------------------------------------------------------------------------
CH
CH	Change Log:
CH
CH	   Shirley M. Read	Corrected error in computation of record array 
CH	       STX              pointers which was causing an access of status
CH	   September, 1988      bits from Word31 records beyond the number of
CH				records which had been read. Also corrected an
CH			        error in filling the data array with no_data.
CH
CH	   Version 4.1.1 10/24/88, SPR 2678, Shirley M. Read, STX
CH		FTB_Word31 does not give the correct match between the plot
CH		labels and the actual use of the bit fields. The bits 
CH		corresponding to the labels (selected by named bits on the
CH		command line) are extracted by the program in VAX bit order
CH	        (lsb = bit 0, msb = bit 7) instead of the telemetry bit order
CH		(msb = bit 0, lsb = bit 7).
CH
CH	   Version 4.2.1  1989 May 10, SPRs 3045, 3799.  Fred Shuman, STX.
CH		3045:  Change name from FTB_Fill_Buffers to FTB_Fill_Word31 to
CH			distinguish from similar routine in the new facility
CH			FTB_TempC.
CH		3799:  Label the plot with actual start time of data plotted.
CH
CH	   Version 4.5.1  1989 Nov 06,  SPR 4784.  Fred Shuman, STX.
CH	      Word31 occasionally retraces a major frame's worth of
CH	      data.  This happens every time there is an overlap of an odd
CH	      number of MjF's between neighboring segments, a 50-50 occurrence
CH	      when a Playback Ingest is interrupted and re-started in a way
CH	      which writes overlapping segments to the archive.  CT_Read_Arcv
CH	      delivers the first NFS_ANC record (= 2 MjF's !) in the new
CH	      segment with a start-time later than that of the last record in
CH	      the old segment.  This segment has a 50% chance of being 1 MjF
CH	      later than the previous one, giving rise to a 1-MjF overlap in
CH	      Word31's plot data.  CT_Query_Get, which is used by FGA, for
CH	      example, apparently returns the 1st segment which does not
CH	      overlap the previous one.  The fix to this SPR is to mimic the
CH	      behavior of CT_Query_Get and skip any record that overlaps its
CH	      predecessor.  This can cause a 1-MjF gap in the plot.
CH	      FTB_Word31, FTB_Fetch_Data, FTB_Read_Word31, FTB_Graph_Word31,
CH	      and FTB_Fill_Word31.
C----------------------------------------------------------------------------

	Implicit None

	Dictionary 'nfs_anc'
	Record/nfs_anc/ word31(1)

	Integer*4 num_recs	! Number of records in time range
	Integer*4 num_pts	! Number of plot points in time range
	Integer*4 mult_fac	! Number of minor frames in a major frame / step
	Integer*4 num_fields	! Number of bit fields requested for plot
	Integer*4 step		! Number of word31 values to skip each iteration
	Integer*4 bits(num_fields)  ! Bits to extract from word31
	Integer*4 i4
	Integer*4 i, fr, j, k
	Integer*4 LIB$SubX
	Integer*4 LIB$EDiv
	Integer*4 diff_time(2)
	Integer*4 diff_mj(2)
	Integer*4 diff_rec
	Integer*4 sec60/600000000/
	Integer*4 quo
	Integer*4 rem
	Integer*4 cnt
	Integer*4 marker
	Integer*4 status
	Integer*4 bit_num

	Real*4 data(2*num_recs*(mult_fac+1),1)  ! Word31 bit data - output matrix
	Real*4 delta
	Real*4 no_data/-1.2E-34/
	Real*4 max			! Max time - output parameter

	Integer*4 time_flags(1)		! Flags indicating time gaps

	Common /data_info/ diff_mj, diff_rec

	FTB_Fill_Word31 = 1

C   Determine the time span between successive minor frames to be displayed.
C   Fill in the time as the abscissa, interpolating where necessary.

	delta = (floatj(diff_mj(1)) * floatj(step))/(128. * floatj(sec60))
	cnt = 0

	Do i = 1, num_recs

	  cnt = cnt + 1

C   If time difference from the previous record to this one was long, insert a
C   point with the previous record's time, whose plot values will later be
C   filled with 'no_data' to create a gap in the PLT plot.

	  If (time_flags(i) .Eq. 2) Then
	    data(cnt,1) = data(cnt-1,1)
	    cnt = cnt + 1
	  End If

C   Use major frame time as base.  If time difference was short, we will later
C   put 'no_data' into the plot values to create a gap in the PLT plot.  Convert
C   from seconds to minutes.

	  status = LIB$SubX(word31(i).ct_head.time, word31(1).ct_head.time,
	2                   diff_time)
	  status = LIB$EDiv(sec60, diff_time, quo, rem)
	  data(cnt,1) = floatj(quo) + floatj(rem)/sec60

C   Interpolate to minor frame time

	  If (time_flags(i) .Ne. 1) Then
	    marker = cnt
	    Do k = 1,2
	      Do j = 1+step, 128, step
	        cnt = cnt + 1
	        data(cnt,1) = data(cnt-1,1) + delta
	      End Do

	      If (k .Eq. 1) Then           ! First major frame
	        cnt = cnt + 1
	        data(cnt,1) = data(marker,1) + floatj(diff_mj(1))/sec60
	      End If
	    End Do   ! k = 1,2
	  End If

	End Do       ! i = 1, num_recs

	max = data(cnt,1)

C   Extract bits and fill data arrays.  Reset the word31 record counter for
C   each field as well as the data array pointers.

	Do k = 1, num_fields
	  cnt = 0
	  Do i = 1,num_recs

C   If time difference is wrong, use no-data flag.  Increment the data pointer.

	    If (time_flags(i) .Ne. 0) Then
	      cnt = cnt + 1
	      data(cnt, k+1) = no_data
	    End If

C   If data doesn't overlap previous record, get it.

	    If (time_flags(i) .Ne. 1) Then

	      Do fr = 1,2
	        Do j = 1,128,step
	          i4 = word31(i).frame(fr).minor_frame_status_bits(j)
	          cnt = cnt + 1
	          bit_num = 7 - bits(k)			! 10/24/88
	          data(cnt, k+1) = floatj(jibits(i4, bit_num, 1))
	        End Do
	      End Do

	    End If

	  End Do     ! i = 1,num_recs
	End Do       ! k = 1, num_fields

	num_pts = cnt

	Return
	End
