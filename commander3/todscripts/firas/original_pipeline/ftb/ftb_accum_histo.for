	Integer*4 Function ftb_accum_histo(chan_id, num, yvals, inter,
     ^                                     starting_time, ending_time )

c--------------------------------------------------------------------------
c	Author: RAShafer
c		685
c		1988 Oct 5
c		Adapted from FTB_READ_HISTO by Bob Kummerer
c
c	Input:	Channel number
c		the accumulation histogram
c
c	Output:	number of records found
c		the accumulation histogram
c
c	Include files:
c		ct$library:ctuser.inc
c
c	Calling Sequence
c		status = ftb_accum_histo(channel number,
c					number of records found,
c					the time series,
c					fake-it flag)
c	Gist of routine
c		   QUERY for start and stop times
c		   PAD with zeros if necessary
c		CALL CT_OPEN to open science data for reading
c		IF return status is an error THEN
c		   SET return status to an error
c		ELSE
c		      DO UNTIL an error occurs or end of file
c			 CALL CT_READ_ARCV to read science records
c		         INCREMENT record counter
c			 Accumulate histogram
c		      END DO
c		      IF an error occured THEN
c		         SET return status to ERROR
c		      ELSE
c	                 SET return status to normal
c		      END IF
c                  END IF
c               END IF
c	        RETURN
c	        END
c___________________________________________________________________________
c Modification History
c
c D. Bouler, STX  april 24, 1989  SPS's 3664 and 3128
c    Divide ifg number by the number of sweeps so that bin number is
c    in the range 1 to 4096.
c
c J. Smid,  STX August 8, 1989  SPR 3883
c    The implementation of the noninteractive mode
c
c Version 4.5   1989 Sep 25,  SPR 4517,  Fred Shuman, STX.
c    Clear yvals() array each time a new timerange is requested.
c    FTB_HISTO and FTB_ACCUM_HISTO.
C
C A. Panitz, STX November 3, 1989 SPR 4900
C    Pass Start and End times back to FTB_HISTO so they can go on plot label
c
c D. Bouler STX Nov 8, 1989 spr 4797
c    ftb_accum_histo was initially being set to ftb_normal at the beginning
c    of the subroutine, but was returning an anonymous successful completion
c    message at the end. Ftb_histo, the main routine, did not signal
c    an error if it did not get ftb_normal return from ftb_accum_histo;
c    it just fell through and signalled ftb_aberr. Added ftb_status to
c    receive status return, then just before return set ftb_accum_histo
c    equal to ftb_status.
c
c
c--------------------------------------------------------------------------

	Implicit None

	Include '(FUT_Params)'
	Include 'CT$Library:CTUser.Inc'

C  Function Arguments
	Integer * 4	chan_id		!channel id
	Integer * 4	num		!number of records found
	real	* 4	yvals(*)	! The accumulation array
	logical	* 1	inter
	Character * 14	starting_time
	Character * 14	ending_time

C  Local Variables
	Integer * 4	start_chan	!first channel to be run
	Integer * 4	i		!a counter
	integer * 4	k		! another counter
	integer * 4	nbin		! The bin number of the current counts
	integer * 4     num_bad_records ! number of bad records found
	Integer	* 4     status
	Integer	* 4	start_time(2)
	Integer	* 4	end_time(2)
	Integer * 4	fake_it		!fake it flag
	Integer * 4	fak_tim(2)	!time to find fake it bit
	Integer	* 4	FUT_Timerange
	integer * 4     test_bin        !test number of ifg counts
	integer * 4     len1
	integer * 4     len2
	integer * 4     cli$get_value
	integer * 4     int
	integer * 4     ios
	integer * 4     ct_connect_read
	integer * 4     lib$get_lun

	integer * 2     sweeps          ! number of sweeps in data
	Integer * 2	ct_unit		!CT logical unit
	Integer * 2	ct_stat(20)	!CT return status
	Integer * 2	ct_type(4)	!CT data type
	Integer * 2	check_fake	!temporary variable

	Character * 14  zero         /'00000000000000'/
	Character * 30	time_span	!combined start-stop times
	Character *  3  input
	Character *  2  chan
	Character * 64  infile


	External	ftb_normal
	External	ftb_maxreccou
	External	ftb_ctinit
	External	ftb_ctopen
	External	ftb_ctread
	External	ftb_ctclos
	External        ct_connect_read

	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ sci_data	!science records

c   Data Testing Information

	logical*4 fut_goodsci
	logical*4 alltests/.true./
	logical*1 tests(2)
	logical*1 flags(2)
	logical*4 bad_bins           ! flag for bad bin numbers in a record

	integer         *4      str$upcase
	integer         *4      ftb_status

C
C Fill the science data buffer.
C
	ct_type(1) = ctu_$fir_rs1
	ct_type(2) = ctu_$fir_rs2
	ct_type(3) = ctu_$fir_rs3
	ct_type(4) = ctu_$fir_rs4

	ftb_status = %loc(ftb_normal)
	ct_stat(1) = ctp_normal
	num_bad_records = 0

C
C Fetch a timerange
C
	If ( inter ) then
	   status = FUT_Timerange ( starting_time, start_time,
     .                              ending_time, end_time )
	   Call CT_Binary_To_GMT(start_time,starting_time)
	   Call CT_Binary_To_GMT(end_time,ending_time)
	else
	   status = Cli$get_value ('JSTART', starting_time, len1 )
	   status = Cli$get_value ('JSTOP', ending_time, len2 )
	   starting_time(len1+1:14) =  zero(1:14-len1)
	   ending_time(len2+1:14) =  zero(1:14-len2)
	endif

	time_span = starting_time // ';' // ending_time // ';'

	If (inter) then
	   type 5, fac_channel_ids(chan_id)
5	   format (/,' Fetching data for channel ', a2)
	   chan = fac_channel_ids(chan_id)
	   type 55
55         format (' Enter data set type ( RAW, FPP, FDQ ):'$)
	   accept 56, input
56         format (a3)
	   status = str$upcase (input,input)

	else
	   status = Cli$get_value ('INPUT', input, len2 )
	   status = Cli$get_value ('CHANNEL', Chan, len2 )
	endif

	If ( input(1:1) .eq. 'R' ) then
	   infile = 'CSDR$FIRAS_ARCHIVE:NFS_SDF_' // Chan //'/'//Time_Span
	else if ( input(1:2) .eq. 'FP' ) then
	   infile = 'CSDR$FIRAS_ARCHIVE:FPP_SDF_' // Chan //'/'//Time_Span
	else
	   infile = 'CSDR$FIRAS_ARCHIVE:FDQ_SDF_' // Chan //'/'//Time_Span
	end if
	status = lib$get_lun(ct_unit)

	open  ( unit = ct_unit,
     .          file = infile,
     .          status = 'old',
     .          iostat = ios,
     .          useropen = ct_connect_read )

	If ( ios .ne. 0 ) then
	   call lib$signal(ftb_ctopen, %val(1), %val(ct_stat(1)))
	   ftb_status = %loc(ftb_ctopen)
	else
	   call ct_read_arcv(,ct_unit, sci_data, ct_stat)
	   num = 0
	   do while(ct_stat(1) .eq. ctp_normal )
	      if( fut_goodsci( sci_data, alltests, tests, flags) )then
	         bad_bins = .false.
	         num = num + 1
	         sweeps = sci_data.sci_head.sc_head11
c
c   Divide ifg counts by the number of sweeps.
c   Check record for bad bin numbers; if found, discard current record
c
	         do k=1, 512
	            test_bin = (sci_data.ifg_data.ifg(k) / sweeps) + 2049
	            if ((test_bin .lt. 0) .or. (test_bin .gt. 4096))
     *                      bad_bins = .true.
	         enddo
	         if ( bad_bins ) then
	            if ( inter ) then
	                write(6,*) 'Found bad bin numbers in record ',num
	                num_bad_records = num_bad_records + 1
	            endif
	         else
	            do k = 1, 512
	                nbin = (sci_data.ifg_data.ifg(k) / sweeps) + 2049
	                yvals(nbin) = yvals(nbin) + 1
	            end do
	         endif
	      else
	         if ( inter ) then
	            write(*,'(2a)')' **** Bad IFG at: ',sci_data.ct_head.gmt
	            if(.not.flags(1)) write(*,*)' Failed Data Ready Flag test'
	            if(.not.flags(2)) write(*,*)' Failed Data Checksum test'
	         endif
	      end if
	      call ct_read_arcv(,ct_unit, sci_data, ct_stat)
	   end do
	   if (ct_stat(1) .eq. ctp_endoffile ) then
	      call ct_close_arcv(,ct_unit,ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         ftb_status = %loc(ftb_ctclos)
	         call lib$signal(ftb_ctclos, %val(1), %val(ct_stat(1)))
	      end if

	   else
	      ftb_status = %loc(ftb_ctread)
	      call lib$signal(ftb_ctread, %val(1), %val(ct_stat(1)))
	   end if
	end if

	If ( inter ) then
	   write(6,*) 'Number of records with bad bin numbers: ',num_bad_records
	endif

	num = num - num_bad_records
	ftb_accum_histo = ftb_status

	return
	end
