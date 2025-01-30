C------------------------------------------------------------------------

	Integer*4 Function FTB_Overlap_Search( jstart, jstop,
	1				       report, report_lun )

C------------------------------------------------------------------------
C    PURPOSE: Search the NFS_SDF raw science archives for overlaps and
C	      produces a list of non-overlapping timeranges for FPP
C	      to run.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            October 13, 1989
C
C    INVOCATION: status = FTB_Overlap_Search ( jstart,
C					       jstop,
C					       report,
C					       report_lun )
C
C    INPUT PARAMETERS:
C
C	JSTART(2)	I*4	GMT binary start time of the period
C				to be checked for missing segment
C				catalog entries and data gaps.
C	JSTOP(2)	I*4	GMT binary stop time of the period
C				to be checked for missing segment
C				catalog entries and data gaps.
C	REPORT		L*1	Write a report?
C	REPORT_LUN	I*4	LUN to write report on.
C
C    OUTPUT PARAMETERS: 
C
C	STATUS		I*4	Success status.
C
C    SUBROUTINES CALLED: 
C	CCT_Query_Tod_Catalog
C	CCT_Get_FileSpec_Fields
C	FTB_Remove_Overlap
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	CTUser.Inc
C	CCT_Query_TOD_Catalog_Record
C	CCT_FileSpec_Fields_Record
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit None

	Include	'($SSDef)'
	Include 'CT$LIBRARY:CTUser.Inc'
        Include '(CCT_Query_Catalog_Record)'
	Include '(CCT_Query_TOD_Catalog_Record)'
	Include	'(FUT_Params)'
	Include	'(CUT_Params)'

	Dictionary 'CCM_CME_Catalog_Entry'
	Record /CCM_CME_CATALOG_ENTRY/ cats(CUT$_Max_Segments,4)

	Record /QUERY_TOD_CATALOG/ cat_query

	Integer		*4	jstart(2)
	Integer		*4	jstop(2)
	Logical		*1	report
	Integer		*4	report_lun

	Integer		*4	status
	Integer		*2	ct_status(20)
	Integer		*4	i
	Integer		*4	j

	Integer		*4	chan
	Integer		*2 	nrecs(4)
	Integer		*4	fpp_jstart(2,CUT$_Max_Segments)
	Integer		*4	fpp_jstop(2,CUT$_Max_Segments)
	Integer		*2 	ntimes
	Character	*19	gmt_jstart
	Character	*19	gmt_jstop
	Character	*32	file_ext

	Integer		*4 	CCT_Query_Catalog
	Integer		*4 	CCT_Query_TOD_Catalog
	Integer		*4	FTB_Remove_Overlap
	Integer		*4	FTB_Refine_Overlap

	External	FTB_Normal
	External	FTB_CTInit
	External	FTB_CTQueryCat
	External	FTB_NoCatRecs

C Initialize.

	Call CT_Init ( ct_status )

	If (ct_status(1) .Ne. CTP_Normal) Then
          FTB_Overlap_Search = %Loc(FTB_CTInit)
          Call LIB$Signal ( FTB_CTInit, %Val(1), %Val(ct_status(1)) )
	  If (report) Then
	     Write (report_lun, 80) ct_status(1)
	  End If
	  Return
	End If

C Process all channels.

	Do chan=1,4

C Check all raw science segments in each channel.  Query the archive catalog
C for raw science segments that are in the specified timerange.

	  cat_query.archive_id = 'CSDR$FIRAS_RAW'
	  cat_query.dataset_name = 'NFS_SDF_' // fac_channel_ids(chan)
	  cat_query.start_time(1) = jstart(1)
	  cat_query.start_time(2) = jstart(2)
	  cat_query.stop_time(1) = jstop(1)
	  cat_query.stop_time(2) = jstop(2)

	  status = CCT_Query_TOD_Catalog ( cat_query, cats(1,chan),
	1				   nrecs(chan) )

	  If (.Not. status) Then
	     nrecs(chan) = 0
             Call LIB$Signal ( FTB_CTQueryCat, %Val(1), %Val(status) )
	     If (report) Then
	        Write (report_lun, 90) status
	     End If
	  End If

	  If (nrecs(chan) .Gt. 0) Then

C Search the set of NFS_SDF raw science catalog records for overlaps.
C The returned list of catalog records contain non-overlapping timeranges.

	     status = FTB_Remove_Overlap ( cats(1,chan), nrecs(chan),
	1				   report, report_lun )

	  Else

	     Call LIB$Signal(FTB_NoCatRecs,%Val(1),%Val(chan))
	     If (report) Then
	        Write (report_lun, 5) fac_channel_ids(chan)
	     End If

	  End If

	End Do

C Form a list of timeranges on which FPP can be invoked.  The non-overlapping
C timeranges are currently for each channel.  However, FPP cannot be given a
C timerange to process each channel by, but a single timerange.  This routine
C merges times from each channel into a single list of channel independent
C times.

	status = FTB_Refine_Overlap ( cats, nrecs, fpp_jstart, fpp_jstop,
	1			      ntimes )

C List the non-overlapping timeranges on which FPP may be run.

	Write (6,200)
	If (report) Then
	   Write (report_lun, 200)
	End If

	Do i=1,ntimes

	   Call FUT_Format_Time ( fpp_jstart(1,i), gmt_jstart )
	   Call FUT_Format_Time ( fpp_jstop(1,i),  gmt_jstop  )

	   Write (6,300) i, gmt_jstart, gmt_jstop
	   If (report) Then
	      Write (report_lun, 300) i, gmt_jstart, gmt_jstop
	   End If

	End Do

	FTB_Overlap_Search = %Loc(FTB_Normal)

5	Format (//, 5x, '*** No NFS raw science catalog records found in channel ', a, '.')
80	Format (//, 5x, '*** CT_INIT error = ', i)
90	Format (//, 5x, '*** CCT_QUERY_TOD_CATALOG error = ', i)
200	Format (//, ' Candidate FPP run times:',/)
300	Format (    '     (', i3, ')  ', a, ' to ', a)

	Return
	End
