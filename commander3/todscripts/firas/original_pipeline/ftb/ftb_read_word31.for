	Integer*4 Function FTB_Read_Word31(word31, time_flags, num_recs,
	2                                  data_start, lun)

C---------------------------------------------------------------------------
C
C      Purpose:
C
C          To read the desired data and form the time flags array for the 
C          purpose of identifying time gaps.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C---------------------------------------------------------------------------
C      Modifications:
C
C	   Shirley M. Read
C	       STX
C	   August, 1988
C	   Reason:  Output messages via Lib$Signal.
C
C	   Version 4.5.1  1989 Nov 06,  SPR 4784.  Fred Shuman, STX.
C	      Word31 occasionally retraces a major frame's worth of
C	      data.  This happens every time there is an overlap of an odd
C	      number of MjF's between neighboring segments, a 50-50 occurrence
C	      when a Playback Ingest is interrupted and re-started in a way
C	      which writes overlapping segments to the archive.  CT_Read_Arcv
C	      delivers the first NFS_ANC record (= 2 MjF's !) in the new
C	      segment with a start-time later than that of the last record in
C	      the old segment.  This segment has a 50% chance of being 1 MjF
C	      later than the previous one, giving rise to a 1-MjF overlap in
C	      Word31's plot data.  CT_Query_Get, which is used by FGA, for
C	      example, apparently returns the 1st segment which does not
C	      overlap the previous one.  The fix to this SPR is to mimic the
C	      behavior of CT_Query_Get and skip any record that overlaps its
C	      predecessor.  This can cause a 1-MjF gap in the plot.
C	      FTB_Word31, FTB_Fetch_Data, FTB_Read_Word31, FTB_Graph_Word31,
C	      and FTB_Fill_Word31.
c
c           version 4.5.2 Spr 5162 Dec 1, 1989    D. Bouler STX
c             If telemetry is in dump or engineering mode, time for major
c             frame 2 is not filled in, but set to zero. When the two times
c             are subtracted and multiplied by two to get time between major
c             frames, integer overflow occurrs. Put in check while reading
c             ancillary data to discard data unless telemetry mode is science
c             (nfs_anc.ct_head.hskp(1 and 2)_tlm_fmt = 1) . Only this file
c             was modified. 
c
C---------------------------------------------------------------------------

	Implicit none

	Include 'ct$library:ctuser.inc'

	Dictionary 'nfs_anc'
	Record/nfs_anc/ word31(1)

	Integer*4	time_flags(1)
	Integer*4	num_recs
	Character*14	data_start
	Integer*4	lun

	Integer*4	i
	Integer*4	j
	Integer*4	status
	Integer*2	ct_stat(20)
	Integer*4	min_fram /2500000/
	Integer*4	sec80  /800000000/
	Integer*4	LIB$SubX
	Integer*4	LIB$EDiv
	Integer*4	LIB$EMul
	Integer*4	diff(2)
	Integer*4	diff_mj(2)
	Integer*4	diff_rec
	Integer*4	quo, rem
        integer*4       tformat_1   ! telemetry format for major frame 1
        integer*4       tformat_2   ! telemetry format for major frame 2
        integer*4       num_skip    ! number of records not in science format


	External	FTB_CTRead

	Common /data_info/ diff_mj, diff_rec

	FTB_Read_Word31 = 1
	ct_stat(1) = ctp_normal
	num_recs = 0
        num_skip = 0
        tformat_1 = 0
        tformat_2 = 0

C   Read the data; make sure telemetry format is science or don't count it

        do while (ct_stat(1) .eq. ctp_normal)  
          num_recs = num_recs + 1
          Call CT_Read_Arcv(, lun, word31(num_recs), ct_stat)
          tformat_1 = word31(num_recs).ct_head.hskp1_tlm_fmt
          tformat_2 = word31(num_recs).ct_head.hskp2_tlm_fmt
          if ((tformat_1 .ne. 1) .or. (tformat_2 .ne. 1)) then 
               num_recs = num_recs - 1
               num_skip = num_skip + 1
          endif
	End Do
        num_recs = num_recs - 1        

        if (ct_stat(1) .ne. ctp_endoffile) then
          call lib$signal(FTB_ctread, %val(1), %val(ct_stat(1)) )
          ftb_read_word31 = 0

        else if (num_recs .eq. 0) then

          ftb_read_word31 = 0

	Else if (num_recs .gt. 0) then

c   The start of data will be where telemetry format was 1 (science mode).

          data_start = word31(num_skip+1).ct_head.gmt

C   Find the time difference between rec(j) and rec(j-1).  If > 1.25 rec or
C    < 1 rec, set the time flags to indicate a gap or overlap, respectively.

	  time_flags(1) = 0
 	  status = LIB$SubX(word31(1).gmt_mjf2, word31(1).ct_head.time, diff_mj)
	  diff_rec = 2 * diff_mj(1)

	  Do j=2,num_recs

	    status = LIB$SubX(word31(j).ct_head.time,
	2                     word31(j-1).ct_head.time, diff)

	    If (diff(2) .Eq. 0  .And.
	2         diff(1) .Ge. 0  .And.  diff(1) .Lt. diff_rec-min_fram) Then
	      time_flags(j) = 1
	    Else If (diff(1).Gt.sec80 .Or. diff(1).Lt.0 .Or. diff(2).Gt.0) Then
	      time_flags(j) = 2
	    Else
	      time_flags(j) = 0
	    End If

	  End Do  !  (j=2,num_recs)

	End If    !  (ct_stat .eq. ctp_normal)

	Return
	End
