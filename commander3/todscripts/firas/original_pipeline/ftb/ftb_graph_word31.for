	Integer*4 Function FTB_Graph_Word31(data_address, num_fields, bits,
     2                   flag_address, num_recs, step,
     3                   start_time, end_time, data_start, ylabels, yslab,
     4                   inter,plotfile,plt_com,plt_com_file)

C----------------------------------------------------------------------------
C
C      Purpose:
C
C          Fills the ftb_word31 arrays for display using ftb_fill_Word31 and
C          creates the command lines for PLT.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C----------------------------------------------------------------------------
CH
CH     Change Log:
CH
CH	   Shirley M. Read	Modified titles, ylabels and character size
CH	       STX		so that the titles and labels would fit on
CH	   September, 1988	the laser plots. If more than 3 bits are to
CH			        be plotted, the Y labels must be shortened.
CH			        A short Y label array is now passed to
CH				FTB_Graph_Data.
CH
CH	   Version 4.1.1 10/24/88, SPR 2678, Shirley M. Read, STX
CH		FTB_Word31 does not give the correct match between the plot
CH		labels and the actual use of the bit fields. The bits
CH		corresponding to the labels (selected by named bits on the
CH		command line) are extracted by the program in VAX bit order
CH		(lsb = bit 0, msb = bit 7) instead of the telemetry bit order
CH		(msb = bit 0, lsb = bit 7).
CH
CH	   Version 4.2.1  1989 Apr 26, SPR 3673.  Fred Shuman, STX.
CH		Y-axis numeric plot labels in FTB_Word31 overlap, rendering
CH		them unreadable.
CH
CH	   Version 4.2.1  1989 May 10, SPRs 3045, 3799.  Fred Shuman, STX.
CH		3045:  Change name from FTB_Graph_Data to FTB_Graph_Word31 to
CH			distinguish from similar routine in the new facility
CH			FTB_TempC.
CH		3799:  Label the plot with actual start time of data plotted.
CH
CH         Version 4.4.1    1989 August 3, SPR  3663,  Jon Smid, STX
CH              3663:  FTB  Word31 and Histo need to be run in operations
CH                     in noninteractive mode
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
CH
CH	SER 5726, Steven Alexander, STX, March 6, 1990.  Add capability
CH		to pass in a PLT command file.
CH---------------------------------------------------------------------------

	Implicit None

	Integer*4	data_address
	Integer*4	num_fields
	Integer*4	bits(num_fields)
	Integer*4	flag_address
	Integer*4	num_pts
	Integer*4	num_recs
	Integer*4	step
	Character*14	start_time
	Character*14	end_time
	Character*14	data_start
	Character*15	ylabels(num_fields)
	Character*6	yslab(num_fields)

	Integer*4	i
	Integer*4	status
	Integer*4	FUT_Erase
	Integer*4	LIB$Get_VM
	Integer*4	LIB$Free_VM
	Integer*4	bit_address
	Integer*4	mult_fac
	Integer*4	iery(9)
	Integer*4	ier, ncmd
	Integer*4	pos/45/
	Integer*4	LIB$Erase_Page
	Integer*4	FTB_Fill_Word31
	Integer*4	OTS$Cvt_L_TU
	Integer*4	num_fit / 3 /



	Real*4		yplotmin/-0.5/
	Real*4		yplotmax/1.5/
	Real*4		timeplotmin, timeplotmax
	Real*4		max

	Character*150	cmd(50)
	Character*70	title/' FIRAS WORD 31 TELEMETRY BITS (MSB -> LSB): '/
	Character*50	xlabel/'Time (minutes) beginning at '/
	Character*50    plotfile

	Logical*1	plt_com
	Character*64	plt_com_file

        Integer*4	xlen /28/
	Character*1	char
        logical*1       inter


C   Calculates the number of display points (minor frames) within a major
C   frame. Allocate the memory necessary for the status bit values. Fill
C   the input array to PLT.

	mult_fac = 128/step
	status = LIB$Get_VM(8*num_recs*(mult_fac+1)*(num_fields+1), bit_address)

	status = FTB_Fill_Word31(%val(bit_address), step, num_recs, num_pts,
     2                         bits, num_fields, mult_fac, %val(flag_address),
     3                         %val(data_address), max)

C   Create the PLT commands

	Do i=1,num_fields+1
	  iery(i) = 0.
	End Do

	status = OTS$Cvt_L_TU(bits(1), char,,)
	If (bits(1) .Eq. 0) Then
	  char = '0'
	End If
	title(pos:pos) = char

	xlabel(xlen+1:) = data_start

	Do i = 2,num_fields
	  pos = pos + 1
	  status = OTS$Cvt_L_TU(bits(i), char,,)
	  title(pos:pos+2) = ', '//char
	  pos = pos + 2
	End Do


        If (inter) Then
	   cmd(2) = 'D /VT'
        Else
           cmd(2) = 'D '//PLOTFILE
        Endif


        cmd(3) = 'LA OT '//title
	cmd(4) = 'LA T Timerange selected: '//start_time//' - '//end_time
	cmd(5) = 'LA F'
	cmd(6) = 'LA X '//xlabel
	cmd(7) = 'CS 1.0'
	cmd(8) = 'V .2 .1 .9 .8'

	ncmd = 8

	timeplotmin = 0.
	timeplotmax = max

C   Note that PLT has a bug such that one gets a floating point 0 divide
C   if one is plotting only 1 graph in VERTICAL mode and the values all
C   happen to be the same. Thus, if there is only one item (ie. the scan
C   direction), I have elected to use the OVERLAY mode. Only a certain
C   number of Y label fields with 15 characters will fit on the Y axis for
C   vertical plots. When this number is exceeded, the shortened labels must
C   be used.

	If (num_fields .Gt. 1) Then

	  cmd(1) = 'PLOT VERTICAL'

	  Do i = 1, num_fields
	    ncmd = ncmd + 1
	    Write (cmd(ncmd), '(a3, i1, 2(1x, f4.1))')
	2                         'R Y', i+1, yplotmin, yplotmax

	    ncmd = ncmd + 1
	    Write (cmd(ncmd), '(a3, i1, a4)')  'G Y', i+1, ' 2 1'

	    ncmd = ncmd + 1
	    If (num_fields .Le. num_fit) Then
	      Write (cmd(ncmd), '(a4, i1, 1x, a15)')  'LA Y', i+1, ylabels(i)
	    Else
	      Write (cmd(ncmd), '(a4, i1, 1x, a6)')  'LA Y', i+1, yslab(i)
	    End If
	  End Do

	Else

	  cmd(1) = 'PLOT OVERLAY'
	  Write (cmd(9), '(a3, 2(1x, f4.1))')  'R Y', yplotmin, yplotmax
	  cmd(10) = 'LA Y '//ylabels(1)
	  ncmd = 10

	End If

	If (plt_com) then
	   ncmd = ncmd + 1
	   cmd(ncmd) = '@' // plt_com_file
	end if

        If (.not.inter) then
               ncmd = ncmd + 2
               cmd(ncmd-1)  = 'P'
               cmd(ncmd)    = 'Q'
        Endif

	status = LIB$Erase_Page(1,1)

	Call PLT (%val(bit_address), iery, 2*num_recs*(mult_fac+1),
     2            num_pts, num_fields + 1, cmd, ncmd, ier)

	status = FUT_Erase()

C   Free the display data memory.

	status = LIB$Free_VM(8*num_recs*(mult_fac+1)*(num_fields+1), bit_address)

	Return
	End
