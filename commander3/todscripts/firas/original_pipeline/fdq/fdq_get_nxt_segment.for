	INTEGER*4 FUNCTION FDQ_GET_NXT_SEGMENT ( OPTIONS, TIME_RANGE, 
	1	  FILE_SEG, REF_CHAN, MORE_SEGMENTS, FILENAMES )
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_NXT_SEGMENT
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine will query the COBETRIEVE catalog to find out what 
C/	  segment or segments are to be processed. it will set up the filenames
C/	  for the next set of segments to process, 4 raw science files and 1 
C/	  housekeeping file.  It will also check to see if the science data
C/	  had been processed by FPP (a must) or already by FDQ (a no-no).
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  April 16, 1987
C/
C/	MODIFIED BY:
C/	  Shirley M. Read
C/	  STX
C/	  January 25, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  March 1988
C/	  REASON:	COBETRIEVE has been modified. The catalog numbers are no
C/	                longer available. The Fortran open with filename and 
C/			user open connect key specifications now replace the 
C/			old CT_Open_Arcv and CT_Modify_Arcv. There are also
C/	                new CT Catalog query functions which will get all the
C/			information necessary to group science and housekeeping
C/			files in a segment and obtain all filenames for any
C/			valid specified time range. The CT_Cat_Info will be 
C/			replaced by the new query functions. The returned 
C/			catalog record contains information to be used by other
C/			FDQ functions and they will be modified accordingly.
C/
C/	  Shirley M. Read
C/	  STX
C/	  May 1988
C/	  REASON:	The start and stop times of the total data range 
C/			are needed for opening the attitude archive. If the
C/		        time range is input by the user, the start and stop
C/		        times can be extracted from it. When the user selects
C/			only an input file (segment), the start and stop must
C/			be obtained from the catalog. This routine was modified
C/			to load the times into the Time_Range  passed parameter
C/			after the catalog query.
C/
C/	  Shirley M. Read
C/	  STX
C/	  July 1988
C/	  REASON:	The start and stop times of the total data range 
C/			are needed for opening the configuration files, even
C/	                if the user specifies a data range covering most
C/			of the data or specifies a single file. The total
C/		        span will now cover all available channel files and
C/			the housekeeping files.
CH
CH	CHANGE LOG:
ch
ch	version 4.1.1 12/01/88, ser 2379, J.T.Bonnell, STX@GSFC
ch		This program was modified to refer to the
ch		new firas archive logical names in response
ch		to SER 2379 (csdr$firas_in, _out, _raw,
ch		_ref, _uref, and _cal).
ch
CH	Version 4.2.1 02/09/89, SPR 2700, Shirley M. Read, STX
CH		FDQ falls over for out of order and bad time ranges.The FIRAS
CH	 	Stripper now writes the science records using the telemetry
CH		minor frame time at transmit for the primary key in order to
CH		ensure a valid time tag. A new FIRAS preprocessor reads the raw
CH		science records, computes the midpoint of collect time and flags
CH		any records with bad midpoint times. FDQ must now access the
CH		computed collect time from the new field in the science record
CH	        in order to group the four science channels with the bracketing
CH		housekeeping frames and set the data quality summary flag when
Ch		the badtime flag is set. The bad record is written back to the
CH		archive immediately so that the next record may be considered
CH	        for the group. The average midpoint of collect times is still
CH		used as the time tag of the FDQ engineering record and is thus
CH		stored in the science record as the cross reference. However,
CH		the transmit times of the science records are stored in the 
CH		FDQ engineering record. FDQ must now enter the science gain into
CH		each science record from the bracketing housekeeping frames. The
CH		two LS bits must now be checked in the certification mask when
CH		running FDQ. The bit must be set for FPP but not for FDQ.
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after
CH		an FPP or FDQ abort. The FPR tracking system compounds the
CH		problems. Files with non-matching version numbers seem often
CH		to result from improper clean-up. Bad record times cause
CH		SEGCTL to abort and mess up the tracking system. It was
CH		decided to change the modify of the science records in FPP
CH		and FDQ to a simple COBETRIEVE read of the existing records
CH		from a dataset and write a modifed dataset with the same
CH		information which was entered on the modify. Two new science
CH		data sets will be required: a science dataset of raw science
CH		data plus FPP input and a science dataset with FPP science
CH		data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH		xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where
CH		xx is the channel id. The new datasets must be opened and
CH		processed in FPP and FDQ. 
CH
CH	Version 4.4.3 11/14/89, SPR 5032, R. Kummerer STX
CH		Skip processing raw science segments from missing channels.
CH
CH
CH	Version 6.7 8/20/90, SPR 7192, H. WANG STX
CH		Inappropriate error message signaled by FDQ.
CH
C/		
C/	CALLING SEQUENCE:
C/	  Return_Status = FDQ_GET_NXT_SEGMENT ( OPTIONS, TIME_RANGE, 
C/	                  FILE_SEG, REF_CHAN, MORE_SEGMENTS, FILENAMES )
C/
C/	INPUT PARAMETERS:
C/	  OPTIONS(10)	I*2	Array of various flags indicating processing
C/				options:
C/				OPTIONS(1)  --	Flag indicating whether a time
C/						range is used:
C/						0 = No time range used (default)
C/						1 = Binary start time used
C/						2 = Delivery time used
C/						3 = Specific segments are to
C/						    be processed
C/						4 = Specific segments are to
C/						    be processed; HKP read by
C/						    timerange
C/				OPTIONS(2)  --	Flag indicating whether all or
C/						one segment(s) are(is) to be
C/						processed:
C/						0 = Process only the next segment
C/						1 = Process all segments (default)
C/				OPTIONS(3)  --	Flag indicating how many segments
C/						are to be excluded from processing
C/						in this run.
C/						0 is the default (that mean no
C/						segment is excluded by default)
C/
C/				OPTIONS(4) - OPTIONS(10) -- Reserved for future use
C/
C/	  TIME_RANGE	   C*30	Time range for processing interval if OPTIONS(1)
C/				equals 1(BST) or 2(DDT)
C/	  FILE_SEG	   C*39 User specified file-segment name.
C/	  REF_CHAN	   I*4  Channel specified in filename.
C/
C/	OUTPUT PARAMETERS:
C/	  TIME_RANGE	   C*30	Time range for processing interval if OPTIONS(1)
C/				equals 3, specific segment selected and total
C/				data range.
C/	  MORE_SEGMENTS(4)  L*1	Flag(one for each channel) to indicate whether 
C/				there are more segments to process or not, set
C/				to "true" if there are more segments to process
C/	  FILENAMES(5)	   C*39 The filenames for the 4 raw science and 1 
C/				housekeeping files to process. Filenames(1 - 4)
C/			        are RH, RL, LH and LL, respectively. Filename(5)
C/				is the housekeeping.
C/				If channel I is "done", or if this channel
C/				has no segment to match a channel with an
C/				earlier segment, the filename will be set to 
C/				blank.
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE (yet)
C/
C/	INCLUDE FILE USED:
C/	  FDQ_OPTIONS.TXT
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	SUBROUTINES CALLED:
C/	  CCT_QUERY_CATALOG
C/	  CCT_QUERY_TOD_CATALOG
C/	  FDQ_GET_VALID_SEGMENTS
C/	  FDQ_SAVE_CAT_INFO
C/	  FDQ_GET_TIME_RANGE
C/	  CT_BINARY_TO_GMT
C/
C/	METHOD USED:
C/	  < PDL for this routine > :
C/
C/	If (first time around) then
C/	  Set FIRST to false;
C/	  If (specific segment requested) then
C/	    Call CCT_Query_Catalog_Record;
C/	    Check validity of input segment ID and if processed already then
C/	      Signal error and set status for return to error;
C/	    Else if OK then
C/	      Setup filenames for segment for CT user open;
C/	    Endif
C/	    Return;
C/	  Endif;
C/	  If (time range used (from OPTIONS)) then 
C/	    Do for each science channel
C/	      Call CCT_Query_TOD_Catalog with time range to get all raw science
C/		 catalog entries; 
C/	      Call FDQ_SAVE_CAT_INFO to pick up info from CCM_CME_CATALOG_ENTRY
C/	         and store into arrays NSEG (number of segments in the time 
C/		 range) and check for exceeding maximum allowed in the run;
C/	      If (NSEG(channel) > 0) then
C/	        Set MORE_SEGMENTS(channel) to true;
C/	      Else
C/	        Set MORE_SEGMENTS(channel) to false;
C/	      Endif;
C/	    Enddo;
C/	    Call CCT_Query_TOD_Catalog with time range to get the housekeeping
C/	         catalog entries;
C/	    If (number of entries = 0) then
C/	        Signal error and set return status to error;
C/	        Set MORE_SEGMENTS for each channel to false;
C/          Endif;
C/	  Endif;
C/	Endif first time;
C/      If (return status is still OK) then
C/	  Call FDQ_GET_VALID_SEGMENTS to pick out the next "valid" segments
C/	    i.e., the complete filenames for the next unprocessed segment and 
C/	    delivery time for later s/w build; Also obtain the
C/	    corresponding housekeeping filename.
C/	Endif;
C/	If (return status is still OK) then
C/	   Call FDQ_GET_Time_Range to get the complete data time range.
C/	Endif;
C/	Return;
C/	End.
C/
	IMPLICIT NONE

!	External Parameters

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ASCTIMER
	EXTERNAL	FDQ_SAVCATINF
	EXTERNAL	FDQ_GETVALSEG
	EXTERNAL	FDQ_CATNOENTRY
	EXTERNAL	FDQ_CLNCAT
	EXTERNAL        FDQ_GETCATENTF
	EXTERNAL	FDQ_GETCATINFO
	EXTERNAL	FDQ_CHECKENTRY
	EXTERNAL        FDQ_CTCATINFOM
	EXTERNAL	FDQ_CHKENTRYER
	EXTERNAL        FDQ_GETRANGE  
	EXTERNAL        FDQ_NOSCIENCE
	EXTERNAL	CCT_Q_NO_CAT_ENTRY
	EXTERNAL        FUT_NORMAL

!	Functions

	INTEGER*4	FDQ_Save_Cat_info
	INTEGER*4	FDQ_Get_Valid_Segments
	INTEGER*4       FDQ_Get_Time_Range
	INTEGER*4	CCT_Query_Catalog
	INTEGER*4	CCT_Query_TOD_Catalog
	INTEGER*4	FUT_Clean_Catalog
	INTEGER*4	SYS$Asctim
	integer*4	cct_q_no_cat_entry

!	Input/Output Parameters

	INTEGER*2	OPTIONS(10)
	CHARACTER	TIME_RANGE*30
	CHARACTER*39	FILE_SEG
	INTEGER*4	REF_CHAN
	LOGICAL*1	MORE_SEGMENTS(4)
	CHARACTER*39	FILENAMES(5)

!	Include files and dictionary references

	INCLUDE		'(FDQ_OPTIONS)'
	INCLUDE	        '($SSDef)'
	INCLUDE		'CT$LIBRARY:CTUSER.INC'
	INCLUDE		'(FUT_PARAMS)'

	integer*4	max_entries		! MAXimum # of segments
	parameter	(max_entries = 50)	! per channel; this affects
						! the dimension of records
	dictionary 'CCM_CME_Catalog_Entry'	!Info retrieved from CT
	record     /CCM_CME_Catalog_Entry/Catalog_Rec(max_entries)

	include '(CCT_Query_Catalog_Record)'
	record  /Query_Catalog/Query_Catalog_Rec

	include '(CCT_Query_TOD_Catalog_Record)'
	record  /Query_TOD_Catalog/Query_TOD_Catalog_Rec

!     Local storage     

	INTEGER*2	CATP_BAD_RSE
	PARAMETER	(CATP_BAD_RSE = 3)
	INTEGER*2	CATP_INSUFF_VM
	PARAMETER	(CATP_INSUFF_VM = 4)
	INTEGER*2	CATP_INVALID_DB
	PARAMETER	(CATP_INVALID_DB = 2)
	INTEGER*2	CATP_INVALID_ITEM
	PARAMETER	(CATP_INVALID_ITEM = 5)
	INTEGER*2	CATP_NO_ENTRY
	PARAMETER	(CATP_NO_ENTRY = 6)
	INTEGER*2	CATP_NORMAL
	PARAMETER	(CATP_NORMAL = 1)

	CHARACTER	CAT_INFO_MSG(2:6)*40 /
	1		'CAT_INFO status: Invalid Database ',
	1		'CAT_INFO status: Error processing RSE',
	1		'CAT_INFO status: Insufficient virt. mem.',
	1		'CAT_INFO status: Invalid Item Name ',
	1		'CAT_INFO status: No entries found' /,
	1		END_ASCTIM*23, ITEM*5, RSE*128,
	1		START_ASCTIM*23, TIME_OPT*3
	character	msg(0:3)*40 /
	1	        'Segment has not yet been processed by FPP',
	1	        'Segment is enabled for processing',
	1		'Segment has been processed by FDQ already',
	1		'Segment is not RS1.' /
	logical*1	first /.true./
	integer*4	check_entry	! catalog entry status from external fn.
	integer*2	first_colon, i, j, lvalid,
	1		nseg(5), second_colon
	integer*4	end_bntm(2),nentries,
	1		start_bntm(2), status
	integer*2       num_cat_rec     ! Number of catalog records from CT
	character*20    filex(max_entries,5)	! Extension names of files to 
					        ! be processed for each dataset
	logical*1	delflg(max_entries,5)   ! Deletion flag for each file
	integer*4	save_init(2,5)          ! Save array for initial times
	integer*4	save_fin(2,5)           ! Save array for final times
	integer*2       pnt		! Position of DOT in filename
	integer*2       id		! Channel 1 - 4 or HKP (5)
	character*20    file_extension  ! Copy of file extension
	integer*4 	retstat		! Return status
	integer*4 	success / 1 /, error / 2 /  ! Values for status
	integer*4       zero  / 0 /

	character*14     ARC_ID /'CSDR$FIRAS_IN '/
	character*14     ARC_ID_H /'CSDR$FIRAS_RAW'/
	character*10     ARCFILE(4)
	data ARCFILE(1) /'FPP_SDF_RH'/
	data ARCFILE(2) /'FPP_SDF_RL'/
	data ARCFILE(3) /'FPP_SDF_LH'/
	data ARCFILE(4) /'FPP_SDF_LL'/
	character*7 	 ARCHKP /'NFS_HKP'/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Begin Function       !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	retstat = success

	IF (FIRST) THEN
	  first = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							      !
!     A specific segment has been specified by the user.      !
!     Call the Query COBETRIEVE Catalog routine for the info. !
!     Check if selected segment is  there and if it is        !
!     processed already.				      !
!							      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  IF (OPTIONS(OPT_TIME_RANGE) .EQ. OPTV_SPEC_SEG) THEN
	
!	    First get catalog info for the reference channel.

	    query_catalog_rec.archive_id(1:14) = arc_id(1:14)
	    query_catalog_rec.filename(1:39) = file_seg(1:39)

	    status = CCT_query_catalog( query_catalog_rec,catalog_rec)

	    if (.not. status ) then		! Failed to get cat rec
	        retstat = error
	        call lib$signal(FDQ_GETCATENTF,%val(2),%val(status),file_seg)
	    else if (status .eq. %loc(cct_q_no_cat_entry)) then
	        retstat = error
	        call lib$signal( FDQ_CATNOENTRY, %val(1), file_seg)
	    else			! Call to CCT_query_catalog is OK

!	Call FDQ_SAVE_CAT_INFO to store the number of segments (files)    !
!	for the channel, the extension names of files, the deletion       !
!	flags for each file segment and the run flags for each segment.   !
!       After calling FDQ_SAVE_CAT_INFO, NSEG(I) will contain the number  !
!       of segments (in the channel) to be processed, up to a maximum of  !
!       MAX_ENTRIES entries.  FDQ_SAVE_CAT_INFO will also screen out all  !
!       the "excluded" nodes.						  !
									  !
	      nentries = 1
	      id = ref_chan
	      status = FDQ_Save_Cat_Info ( id, nentries, max_entries,
	1	 	   catalog_rec, nseg, filex, delflg )
	      if ( status .ne. %loc(FDQ_NORMAL) ) then
	         if (status .eq. %loc(FDQ_ABERR)) status = zero
	         call lib$signal(FDQ_SAVCATINF,%val(1),%val(status))
	         retstat = error
	      else
 	         more_segments(ref_chan) = .true.
		 file_extension = catalog_rec(1).filename_extension

!	Save the time range for the attitude and configuration files.

		 do i = 1, 2
		   save_init(i,ref_chan) = Catalog_Rec(1).Initial_Time(i)
	         enddo
	         call CT_GMT_To_Binary ( Catalog_Rec(1).Final_Time_Key,
	1	 		     save_fin(1,ref_chan))
	      endif

!	      If reference channel exists and is not processed already, 
!	      get catalog info for other science channels and housekeeping file.

	      If ( Retstat .eq.success ) then
		do  i = 1, 4
	          if ( i .ne. ref_chan ) then
	            query_catalog_rec.filename(1:39) =
	1	      			arcfile(i)//'.'//file_extension
	            status = CCT_query_catalog( query_catalog_rec,
	1		     catalog_rec)
	            if (.not. status ) then  ! Failed to get cat rec
	               retstat = error
	               call lib$signal(FDQ_GETCATENTF,%val(2),
	1		    %val(status),query_catalog_rec.filename)
	            else if (status .eq. %loc(cct_q_no_cat_entry)) then
	                Call Lib$Signal( FDQ_CATNOENTRY, %val(1), 
	1	          query_catalog_rec.filename)
			more_segments(i) = .false.
	            else		! Good status from CCT_query_catalog
		      nentries = 1
		      id = i
		      status = FDQ_Save_Cat_Info ( id, nentries, max_entries,
	1	 	   catalog_rec, nseg, filex, delflg )
	              if ( status .ne. %loc(FDQ_NORMAL) ) then
	                if (status .eq. %loc(FDQ_ABERR)) status = zero
	                call lib$signal(FDQ_SAVCATINF,%val(1),%val(status))
	                retstat = error
		      else
		        more_segments(i) = .true.

!	Save the time range for the attitude and configuration files.

		        do j = 1, 2
		          save_init(j,i) = Catalog_Rec(1).Initial_Time(j)
	                enddo
	                call CT_GMT_To_Binary ( 
	1		  Catalog_Rec(1).Final_Time_Key,save_fin(1,i))
		      endif		! Good status from sav cat info
	    	    endif		! Good status from CCT_query_catalog
		  endif			! Not reference channel
	        enddo			! Channel loop

!	Get the corresponding housekeeping files.

	        if ( retstat .eq. success ) then
	          query_catalog_rec.filename(1:39) =
	1	      			archkp//'.'//file_extension
	          status = CCT_query_catalog( query_catalog_rec,
	1		   catalog_rec)
	          if (.not. status ) then  ! Failed to get cat rec
	             retstat = error
	             call lib$signal(FDQ_GETCATENTF,%val(2),
	1	        %val(status),query_catalog_rec.filename)
	          else if (status .eq. %loc(cct_q_no_cat_entry)) then
	              Call Lib$Signal( FDQ_CATNOENTRY, %val(1), 
	1	          query_catalog_rec.filename)
		      more_segments(i) = .false.
	          else
		    nentries = 1
		    id = 5
		    status = FDQ_Save_Cat_Info ( id, nentries, max_entries,
	1	 	   catalog_rec, nseg, filex, delflg )
	            if ( status .ne. %loc(FDQ_NORMAL) ) then
	              if (status .eq. %loc(FDQ_ABERR)) status = zero
	              call lib$signal(FDQ_SAVCATINF,%val(1),%val(status))
	              retstat = error
	            else

!	Save the time range for the attitude and configuration files.

		      do i = 1, 2
		        save_init(i,5) = Catalog_Rec(1).Initial_Time(i)
	              enddo
	              call CT_GMT_To_Binary ( 
	1		Catalog_Rec(1).Final_Time_Key,save_fin(1,5))
	            endif			! Bad status from sav cat info
	    	  endif			! Not status from CCT_query_catalog
	        endif		        ! Retstat is success
	      Endif			! Retstat is success
	    endif			! Failed to get cat rec for ref chan

!	End of one segment selected.

	  ELSEIF (OPTIONS(OPT_TIME_RANGE) .NE. OPTV_SPEC_SEG) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							                !
!     A time range was selected by the user.                            !
!     Convert time range to start and stop VAX quadwords.               !
!     Query COBETRIEVE TOD catalog to get files in time range.          !
!     Call FDQ_SAVE_CAT_INFO to store the number of segments (files)    !
!     for the channel, the names of files, the number of records in     !
!     each file and the run flags for each segment.                     !
!							                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      first_colon = index(time_range, ';')
	      second_colon = first_colon + 
	1		index(time_range(first_colon+1:), ';')
	      start_asctim = time_range(1:first_colon-1)
	      end_asctim = time_range(first_colon+1:second_colon-1)
	      call gmt_to_binary (%ref(start_asctim), start_bntm)
	      call gmt_to_binary (%ref(end_asctim), end_bntm)

!	Query by data delivery time is not implemented for Build 3.2.
!	      if (options(opt_time_range) .eq. optv_ddt_time_range) then
!	            time_opt = 'ddt'
!	      else 			! this is bst
!	            time_opt = 'bst'
!	      endif

	      query_tod_catalog_rec.archive_id = arc_id 
	      do i = 1, 2
	        query_tod_catalog_rec.start_time(i)= start_bntm(i)
		query_tod_catalog_rec.stop_time(i) = end_bntm(i) 
	      enddo

	      do  i = 1, 4
	        query_tod_catalog_rec.dataset_name = 
	1		arcfile(i)
	        status = CCT_query_tod_catalog( query_tod_catalog_rec,
	1		      catalog_rec, num_cat_rec)
	        if (.not. status ) then  ! Failed to get cat rec
	            retstat = error
	            call lib$signal(FDQ_GETCATENTF,%val(2),%val(status),
	1				query_tod_catalog_rec.dataset_name)
	        else if (num_cat_rec .eq. 0) then
	            Call Lib$Signal( FDQ_CATNOENTRY, %val(1), 
	1	          		query_tod_catalog_rec.dataset_name)
	            more_segments(i) = .false.
	        else		! Good status from CCT_query_catalog
		  status = FUT_Clean_Catalog(catalog_rec,num_cat_rec,
	1					fac_present)
		  if (status .ne. %loc(FUT_NORMAL)) then
	            retstat = error
	            Call Lib$Signal( FDQ_CLNCAT, %val(1), status)
		  else
		    nentries = num_cat_rec
		    id = i
		    status = FDQ_Save_Cat_Info ( id, nentries, max_entries,
	1	 	   catalog_rec, nseg, filex, delflg )
	            if ( status .ne. %loc(FDQ_NORMAL) ) then
	              if (status .eq. %loc(FDQ_ABERR)) status = zero
	                call lib$signal(FDQ_SAVCATINF,%val(1),%val(status))
	                retstat = error
		      else
		        more_segments(i) = .true.

!	Save the time range for the attitude and configuration files.

		        do j = 1, 2
		          save_init(j,i) = Catalog_Rec(1).Initial_Time(j)
	                enddo
	                call CT_GMT_To_Binary ( 
	1		  Catalog_Rec(nentries).Final_Time_Key,save_fin(1,i))
		      endif		! Good status from sav cat info
	    	    endif		! Good status from CCT_query_catalog
		  endif			! Good status from FUT_clean_catalog
	        enddo			! Channel loop

! 	If no raw science channel files were found, signal error.

	        if ((.not. more_segments(1)) .and. (.not. more_segments(2))
	1	  .and. (.not. more_segments(3)) .and. 
	2	  (.not. more_segments(4))) then
		  retstat = error
	          call lib$signal(FDQ_NOSCIENCE)
	        endif

!	Get the corresponding housekeeping files.  We'll get the housekeeping
!	segments only if we are processing by specific segments.

	        If ( retstat .eq. success ) then
	          query_tod_catalog_rec.archive_id = arc_id_h
	          query_tod_catalog_rec.dataset_name = archkp
	          status = CCT_query_tod_catalog( query_tod_catalog_rec,
	1		      catalog_rec, num_cat_rec)
	          if (.not. status ) then  ! Failed to get cat rec
	             retstat = error
	             call lib$signal(FDQ_GETCATENTF,%val(2),%val(status),
	1				query_tod_catalog_rec.dataset_name)
	          else if (num_cat_rec .eq. 0) then
	             call Lib$Signal( FDQ_CATNOENTRY, %val(1), 
	1	          		query_tod_catalog_rec.dataset_name)
	             more_segments(i) = .false.
	          else                  ! Good status from CCT_query_tod_catalog
		    status = FUT_Clean_Catalog(catalog_rec,num_cat_rec,
	1					fac_present)
		    if (status .ne. %loc(FUT_NORMAL)) then
	              retstat = error
	              Call Lib$Signal( FDQ_CLNCAT, %val(1), status)
	  	    else
		      id = 5
	              nentries = num_cat_rec
		      status = FDQ_Save_Cat_Info ( id, nentries, max_entries,
	1	 	   catalog_rec, nseg, filex, delflg )
	              if ( status .ne. %loc(FDQ_NORMAL) ) then
	                if (status .eq. %loc(FDQ_ABERR)) status = zero
	                call lib$signal(FDQ_SAVCATINF,%val(1),%val(status))
	                retstat = error
		      else

!	Save the time range for the attitude and configuration files.

		        do i = 1, 2
		          save_init(i,5) = Catalog_Rec(1).Initial_Time(i)
	                enddo
	                call CT_GMT_To_Binary ( 
	1		  Catalog_Rec(nentries).Final_Time_Key,save_fin(1,5))
	              endif		! Bad status from FDQ_save_cat_info
		    endif		! Bad status from FUT_clean_catalog
	    	  endif			! Not status from CCT_query_catalog
	        Endif			! Retstat is success


	  ENDIF			! Time range specified.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									   !
!	Set the time range to cover the complete data range.               !
!									   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( Retstat .eq. success ) then
	    STATUS = FDQ_GET_TIME_RANGE ( SAVE_INIT, SAVE_FIN,
	1	MORE_SEGMENTS, TIME_RANGE )
	      if ( status .ne. %loc(FDQ_NORMAL) ) then
	        if (status .eq. %loc(FDQ_ABERR)) status = zero
	        call lib$signal(FDQ_GETRANGE,%val(1),%val(status))
	        retstat = error
	      endif
	  Endif	! Retstat is success

	ENDIF	! First

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
!    Call FDQ_GET_VALID_SEGMENTS to get the group of filenames for the     !
!    next valid segment. A segment with the run flag set will be           !
!    excluded. If there are no more segments for a dataset, the flag for   !
!    more_segments will be set to false and the filename to blanks. All    !
!    filenames in the group will have the same file extension name. If a   !
!    dataset does not have a file for the current qroup, the filename will !
!    be set to blank. The others in the group will be processed. If any    !
!    segment in the group has a delete flag set, a warning will be         !
!    signalled and the group will be bypassed.
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	If ( Retstat .eq. success ) then
	  STATUS = FDQ_GET_VALID_SEGMENTS (MAX_ENTRIES, FILEX, DELFLG, 
	1	NSEG, MORE_SEGMENTS, FILENAMES )
	    if ( status .ne. %loc(FDQ_NORMAL) ) then
	        if (status .eq. %loc(FDQ_ABERR)) status = zero
	        call lib$signal(FDQ_GETVALSEG,%val(1),%val(status))
	        retstat = error
	    endif
	Endif	! Retstat is success

!	Set function to return status

	if (retstat.eq.success) then
	  FDQ_GET_NXT_SEGMENT = %loc(FDQ_NORMAL)
	else
	  FDQ_GET_NXT_SEGMENT = %loc(FDQ_ABERR)
	endif

	return
	end
