	Integer*4 Function FTB_Read_TempC(hkp, time_flags, num_recs,
	2                                 data_start, lun)

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
C	Revisions:
C
C	   Version 4.4  1989 May 10, SPR 3045.  Fred Shuman, STX.
C			New facility BOZO to be brought into CSDR standards
C			and renamed FTB_TempC.
C
C	   Version 4.6  1989 Nov 21, SPR 5023.  Fred Shuman, STX.
C		Should not require the integrators to be on to plot.
C
C          Version 7.0  1990 Nov 01, SPR 7628. Nilo Gonzales, STX.
C               Add statement to check telemetry format.
C---------------------------------------------------------------------------

	Implicit None

	Include 'ct$library:ctuser.inc'

	Dictionary 'nfs_hkp'
	Record/nfs_hkp/ hkp(1)

	Integer*4	time_flags(1)  ! "Flags" indicating   \   0 = too short
C				       !  quality of the step  >  1 = OK
C				       !  from previous time  /   2 = too long
	Integer*4	num_recs
	Character*14	data_start
	Integer*4	lun

	Integer*4	j
	Integer*4	status
	Integer*2	ct_stat(20)
	Integer*4	min_fram /2500000/
	Integer*4	sec80  /800000000/
	Integer*4	LIB$SubX
	Integer*4	LIB$EDiv
	Integer*4	LIB$EMul
	Integer*4	LIB$Signal
	Integer*4	diff(2)
	Integer*4	diff_mj(2)
	Integer*4	diff_rec
	Integer*4	bin_mjf2(2)
        integer*4       tformat_1   ! telemetry format for major frame 1
        integer*4       tformat_2   ! telemetry format for major frame 2
        integer*4       num_skip    ! number of records not in science format
	Integer*4	CT_GMT_to_Binary

	External	FTB_CTRead

	Common /data_info/ diff_mj, diff_rec

	FTB_Read_TempC = 1
	ct_stat(1) = ctp_normal
	num_recs = 0
	num_skip = 0
	tformat_1 = 0
	tformat_2 = 0

C   Read the data; make sure that the telemetry format is science, or
C   don't count it   

	Do While (ct_stat(1) .Eq. ctp_normal)
	  num_recs = num_recs + 1
	  Call CT_Read_Arcv(, lun, hkp(num_recs), ct_stat)
	  tformat_1 = hkp(num_recs).ct_head.hskp1_tlm_fmt
	  tformat_2 = hkp(num_recs).ct_head.hskp2_tlm_fmt
	  if ((tformat_1 .ne. 1) .or. (tformat_2 .ne. 1)) then
	      num_recs = num_recs - 1
	      num_skip = num_skip + 1
	  endif
	End Do
	num_recs = num_recs - 1

	If (ct_stat(1) .Ne. ctp_endoffile) Then

	  Call LIB$Signal(FTB_CTRead,%val(1),%val(ct_stat(1)))
	  FTB_Read_TempC = 0

	Else

C   Find the time difference between rec(j) and rec(j-1).  If > 1.25 rec, or
C    < 1 rec, set the time flags to indicate a gap or overlap, respectively.

	  time_flags(1) = 1
	  status = CT_GMT_to_Binary(hkp(1).hskp_tail.gmt_mjf2, bin_mjf2)
	  status = LIB$SubX(bin_mjf2, hkp(1).ct_head.time, diff_mj)
	  diff_rec = 2 * diff_mj(1)

	  Do j=2,num_recs

	    status = LIB$SubX(hkp(j).ct_head.time,
	2                     hkp(j-1).ct_head.time, diff)

	    If (diff(2).Eq.0 .And. diff(1).Ge.0 .And.
	2                          diff(1).Lt.diff_rec-min_fram) Then
	      time_flags(j) = 0
	    Else If (diff(1).Gt.sec80 .Or. diff(1).Lt.0 .Or. diff(2).Gt.0) Then
	      time_flags(j) = 2
	    Else
	      time_flags(j) = 1
	    End If

	  End Do

	End If

	Return
	End

