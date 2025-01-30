	Program FTB_Word31

C-----------------------------------------------------------------------------
C
C      Purpose:
C
C          To produce a graphic display of the status bits in FIRAS word 31
C          over time.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C
C       Include files used:
C	        $ssdef
C
C-----------------------------------------------------------------------------
CH
CH	Change Log:
CH
CH          Shirley M. Read
CH              STX
CH          August 26, 1988
CH
CH	Reason: SPRs 1566 and 1763 requested that all FIRAS facilities
CH		should tell of successful completion and signal errors
CH		via $status for batch processing. Also added interface
CH		to condition handler and new messages.
CH
CH          Shirley M. Read
CH              STX
CH          August 26, 1988
CH
CH	Reason: Modified the length of the Y labels for vertical plots.
CH		If there are too many fields for the vertical plot, the
CH		Y labels must be shortened. Added code for ALL bits.
CH
CH          Shirley M. Read
CH              STX
CH          September 20, 1988
CH
CH	Reason: Corrected the number of bytes per NFS_ANC record from
CH		386 to 392 on the lib$free_vm.
CH
CH	Version 4.1.1 10/24/88, SPR 2678, Shirley M. Read, STX
CH		FTB_Word31 does not give the correct match between the plot
CH		labels and the actual use of the bit fields. The bits
CH		corresponding to the labels (selected by named bits on the
CH		command line) are extracted by the program in VAX bit order
CH		(lsb = bit 0, msb = bit 7) instead of the telemetry bit order
CH		(msb = bit 0, lsb = bit 7).
CH
CH	Version 4.2.1  1989 Apr 26, SPR 3673.  Fred Shuman, STX.
CH		Y-axis numeric plot labels in FTB_Word31 overlap, rendering
CH		them unreadable.  Change made in module FTB_Graph_Word31.FOR.
CH
CH	Version 4.2.1  1989 May 10, SPR 3799.  Fred Shuman, STX.
CH		Label the plot with actual start time of data plotted.
CH
CH	Version 4.4.1 1989 August 4, SPR 3663. Jon Smid, STX
CH              Implementation of noninteractive mode
CH
CH	Version 4.5.1  1989 Nov 06,  SPR 4784.  Fred Shuman, STX.
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
CH	Version 5.1, 1989 Dec 22, SPR 5023 on FTB_TempControl. Fred Shuman, STX.
CH	      When TempControl was fixed (Qfix 769) for SPR 5023, num_pts was
CH	      removed from the calling list of FTB_Fetch_Data.  Since this
CH	      module is also called by FTB_Word31, modification was needed here
CH	      as well.
CH
CH	Version 5.2, 1990 Mar 6, SER 5726. Steven Alexander, STX.
CH	      Add capability to pass in a PLT command file.
C----------------------------------------------------------------------------

	Implicit None

c      FTB messages

	External	ftb_normal
	External	ftb_aberr
	External	ftb_bad_stepsize

c      Condition Handler

	External	fut_error

c      External Functions

	External	cli$_defaulted
	External	cli$_present

c      Include files

	Include		'($ssdef)'

c      Local variables and functions

	Integer*4	i, num_recs
	Integer*4	step
	Integer*4	status
	Integer*4	cli$present
	Integer*4	cli$get_value
	Integer*4	ftb_fetch_data
	Integer*4	bits(8)/8*0/
	Integer*4	FTB_Graph_Word31
	Integer*4	upm_get_longword
	Integer*4	data_address
	Integer*4	num_fields
	Integer*4	flag_address
	Integer*4	lib$free_vm
	Integer*4	fetch_type
	Integer*4	len
	Character*14	start_time
	Character*14	end_time
	Character*14	data_start

	Character*15	ylabels(8)
	Character*6	yslab(8)
	Character*50	plotfile
	Logical*1	plt_com /.false./
	Character*64	plt_com_file 
	Logical*1	proceed / .true. /
	Logical*1	inter/.true./

C   FTB Parameters

	Integer*4	Word31
	Integer*4	TempC
	Parameter	(Word31 = 1)
	Parameter	(TempC = 2)

!   Retrieve the data from COBETRIEVE.

       fetch_type = Word31
       status = FTB_Fetch_Data(fetch_type, data_address, num_recs,
     &                         flag_address, start_time, end_time, data_start)

       if(status .ne. 1) proceed = .false.

!   Determine the number of minor frames the user wishes to plot per time step.

       if ( proceed ) then

         status = upm_get_longword('step', step)

         if((step .le. 0) .or. (step .gt. 128)
     &                            .or. status .ne. ss$_normal) then
           proceed = .false.
           call lib$signal(ftb_bad_stepsize)
         endif
       endif

       num_fields = 0

!   Determine the fields for display

       if ( proceed ) then
	if (cli$present('ALL')) then
	  num_fields = 8
	  do i = 1,8
	    bits(i) = i - 1
	  enddo
	  yslab(1) = 'SCDIR'
	  yslab(2) = 'XCAL'
	  yslab(3) = 'LEN'
	  yslab(4) = 'SPEED'
	  yslab(5) = 'RH_DR'
	  yslab(6) = 'RL_DR'
	  yslab(7) = 'LH_DR'
	  yslab(8) = 'LL_DR'
	else
         if((cli$present('SCAN_DIRECTION') .eq. %loc(cli$_present))
     &   .or. (cli$present('SCAN_DIRECTION') .eq. %loc(cli$_defaulted)))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'SCAN_DIRECTION'
	   yslab(num_fields) = 'SCDIR'
           bits(num_fields) = 0
         endif

         if(cli$present('XCAL_STATUS') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'XCAL_STATUS'
	   yslab(num_fields) = 'XCAL'
           bits(num_fields) = 1
         endif

         if(cli$present('LENGTH') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'LENGTH'
	   yslab(num_fields) = 'LEN'
           bits(num_fields) = 2
         endif

         if(cli$present('SPEED') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'SPEED'
	   yslab(num_fields) = 'SPEED'
           bits(num_fields) = 3
         endif

         if(cli$present('RH_DATA_READY') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'RH_DATA_READY'
	   yslab(num_fields) = 'RH_DR'
           bits(num_fields) = 4
         endif

         if(cli$present('RL_DATA_READY') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'RL_DATA_READY'
	   yslab(num_fields) = 'RL_DR'
           bits(num_fields) = 5
         endif

         if(cli$present('LH_DATA_READY') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'LH_DATA_READY'
	   yslab(num_fields) = 'LH_DR'
           bits(num_fields) = 6
         endif

         if(cli$present('LL_DATA_READY') .eq. %loc(cli$_present))then
           num_fields = num_fields + 1
           ylabels(num_fields) = 'LL_DATA_READY'
	   yslab(num_fields) = 'LL_DR'
           bits(num_fields) = 7
         endif
	endif !All or some
       endif  !Proceed

! Find out if the process is interactive
         if ( cli$present('INTERACTIVE') .eq. %loc(cli$_present) .or.
     &        cli$present('INTERACTIVE') .eq. %loc(cli$_defaulted) ) then
                   INTER = .True.
         else
                   INTER = .False.
         endif
         status = cli$Get_Value('PLOTDEVICE',PLOTFILE,LEN)

! Get the PLT command file if specified

	if ( cli$present('PLTFILE') .eq. %loc(cli$_present) ) then
	   plt_com = .true.
	   status = cli$get_value('PLTFILE',plt_com_file,len)
	end if

!   Graph the status bits

       if ( proceed ) then
         status = FTB_Graph_Word31(data_address, num_fields, bits,
     &            flag_address, num_recs, step, start_time, end_time, data_start,
     &            ylabels, yslab, inter, plotfile, plt_com, plt_com_file )
       endif

!   Deallocate the WORD31 and time flag memory.

       status = lib$free_vm(392*num_recs, data_address)
       status = lib$free_vm(4*num_recs, flag_address)

!   Exit the program.

       if ( proceed ) then
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
       else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
       endif
       end
