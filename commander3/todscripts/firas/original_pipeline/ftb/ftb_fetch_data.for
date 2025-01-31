	Integer*4 Function FTB_Fetch_Data(fetch_type, data_address, num_recs,
	2                      flag_address, start_time, end_time, data_start)

C---------------------------------------------------------------------------
C
C      Purpose:
C
C          To check the specified time span for errors and allocate
C          the necessary memory for the data
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
C	   Shirley M. Read
C              STX
C          August, 1988
C	   Reason: Output messages via Lib$Signal.
C
C	   Shirley M. Read
C              STX
C          September, 1988
C	   Reason: Corrected the number of bytes per NSF_ANC record from
C		   386 to 392 on the lib$get_vm and lib$free_vm.
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
C
C	   Version 4.6  1989 Nov 21,  SPR 5023.  Fred Shuman, STX.
C		Should not require the integrators to be on to plot.
C
C---------------------------------------------------------------------------

	Implicit None

	Include 'CT$Library:ctuser.inc'

C   FTB messages

	External FTB_Stop_Lt_Star
	External FTB_OpnFail_Anc
	External FTB_OpnFail_Hkp

C   External Function

	External  CT_Connect_Read

C   Local variables and functions

	Integer*4	fetch_type
	Integer*4	data_address
	Integer*4	num_recs
	Integer*4	flag_address
	Character*14	start_time
	Character*14	end_time
	Character*14	data_start

	Integer*2	ct_stat(20)
	Integer*4	lun
	Integer*4	ios
	Integer*4	st_time(2)
	Integer*4	ed_time(2)
	Integer*4	diff_time(2)

	Integer*4	FTB_Read_Word31
	Integer*4	FTB_Read_TempC
	Integer*4	CLI$Get_Value
	Integer*4	LIB$SubX
	Integer*4	LIB$Get_VM
	Integer*4	LIB$Free_VM
	Integer*4	LIB$Get_Lun
	Integer*4	LIB$Free_Lun
	Integer*4	num_allocate
	Integer*4	status
	Integer*4	len1, len2
	Integer*4	sec60/600000000/
	Integer*4	quo, rem
	Integer*4	LIB$EDiv
	Integer*4	CT_Connect_Read

	Character*14	zero/'00000000000000'/
	Character*30	time_span

C   FTB Parameters

	Integer*4 Word31
	Integer*4 TempC
	Parameter (Word31 = 1)
	Parameter (TempC = 2)

	FTB_Fetch_Data = 1
	time_span = ' '

C   Get the start and stop times and check that start < stop.

	status = CLI$Get_Value('jstart', start_time, len1)
	status = CLI$Get_Value('jstop', end_time, len2)
	start_time(len1+1:14) = zero(1:14-len1)
	end_time(len2+1:14) = zero(1:14-len2)

	time_span = start_time//';'//end_time//';'

	Call GMT_to_Binary(%ref(start_time), st_time)
	Call GMT_to_Binary(%ref(end_time), ed_time)

	status = LIB$SubX(ed_time, st_time, diff_time,)

	If (diff_time(2) .Lt. 0) Then

	  Call LIB$Signal(FTB_Stop_Lt_Star)
	  FTB_Fetch_Data = 0
	  Return

	Else

C   Calculate a maximal number of records based on the time span. Allocate the
C   memory space necessary. Read the data and close the archive.

	  status = LIB$EDiv(sec60, diff_time, quo, rem)
	  num_allocate = quo + 1

	  Call CT_Init(ct_stat)

	  status = LIB$Get_Lun(lun)

	  If (fetch_type .Eq. Word31) Then

	    Open(UNIT=lun,
	2        FILE='CSDR$FIRAS_Archive:nfs_anc/'//time_span,
	3        IOSTAT=ios, STATUS='old', USEROPEN=CT_Connect_Read)

	  Else If (fetch_type .Eq. TempC) Then

	    Open(UNIT=lun,
	2        FILE='CSDR$FIRAS_Archive:nfs_hkp/'//time_span,
	3        IOSTAT=ios, STATUS='old', USEROPEN=CT_Connect_Read)

	  End If

	  If (ios .Ne. 0) Then

	    If (fetch_type .Eq. Word31) Then
	      Call LIB$Signal(FTB_OpnFail_Anc,%val(1),%val(ios))
	    Else If (fetch_type .Eq. TempC) Then
	      Call LIB$Signal(FTB_OpnFail_Hkp,%val(1),%val(ios))
	    End If

	    FTB_Fetch_Data = 0
	    Return

	  Else

	    If (fetch_type .Eq. Word31) Then
	      status = LIB$Get_VM(392*num_allocate, data_address)
	      status = LIB$Get_VM(4*num_allocate, flag_address)
	      status = FTB_Read_Word31(%val(data_address), %val(flag_address),
	2                              num_recs, data_start, lun)
	    Else If (fetch_type .Eq. TempC) Then
	      status = LIB$Get_VM(576*num_allocate, data_address)
	      status = LIB$Get_VM(4*num_allocate, flag_address)
	      status = FTB_Read_TempC(%val(data_address), %val(flag_address),
	2                             num_recs, data_start, lun)
	    End If

	    If (status .Ne. 1) Then

	      If (fetch_type .Eq. Word31) Then
	        status = LIB$Free_VM(392*num_allocate, data_address)
	        status = LIB$Free_VM(4*num_allocate, flag_address)
	      Else If (fetch_type .Eq. TempC) Then
	        status = LIB$Free_VM(576*num_allocate, data_address)
	        status = LIB$Free_VM(4*num_allocate, flag_address)
	      End If
	      FTB_Fetch_Data = 0

	    Else
	      If (fetch_type .Eq. Word31) Then
	        status = LIB$Free_VM(392*(num_allocate-num_recs),
	2                            392*num_recs + data_address)
	        status = LIB$Free_VM(4*(num_allocate-num_recs),
	2                            4*num_recs + flag_address)
	      Else If (fetch_type .Eq. TempC) Then
	        status = LIB$Free_VM(576*(num_allocate-num_recs),
	2                            576*num_recs + data_address)
	        status = LIB$Free_VM(4*(num_allocate-num_recs),
	2                            4*num_recs + flag_address)
	      End If

	    End If  ! (status .Ne. 1)

	  End If    ! (ios .Ne. 0)

	End If      ! (diff_time(2) .Lt. 0)

	Call CT_Close_Arcv(, lun, ct_stat)
	status = LIB$Free_Lun(lun)

	If (ct_stat(1) .Ne. ctp_normal) Then
	  FTB_Fetch_Data = 0
	End If

	Return
	End
