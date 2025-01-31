      Integer*4 Function FUT_Query_RSE(dataset,jstart,jstop,rse)

C------------------------------------------------------------------------
C    PURPOSE: Fetches a configured RSE from the specified RSE archive.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            May 17, 1989
C
C    INVOCATION: status = FUT_Query_RSE ( dataset, jstart, jstop, rse )
C
C    INPUT PARAMETERS:
C
C	DATASET		CH*12	Archive dataset name.
C	JSTART(2)	I*4	GMT binary start time.
C	JSTOP(2)	I*4	GMT binary stop time.
C
C    OUTPUT PARAMETERS: 
C	STATUS		I*4	Success status.
C	RSE(16)		CH*128	Record selection expression.
C
C    SUBROUTINES CALLED: 
C	CCT_Query_TTG_Catalog
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	CTUser.Inc
C	CCT_Query_TTG_Catalog_Record
C	$SSDef
C
C	SPR 7476 -- Wrong RSE retrieved from reference archive.
C		S. Alexander, Oct 17, 1990; STX.
C
C----------------------------------------------------------------------

	Implicit None

	Include	'($SSDef)'
	Include '(CCT_Query_TTG_Catalog_Record)'
	Include '(CUT_Params)'
	Include '(FUT_Params)'

	Dictionary 'CCM_CME_Catalog_Entry'
	Record /CCM_CME_CATALOG_ENTRY/ old_cats(cut$_max_segments),
     .			new_cats(cut$_max_segments),final_cats

	Record /QUERY_TTG_CATALOG/ ttg_query

	Character	*(*)	dataset
	Integer		*4	jstart(2)
	Integer		*4	jstop(2)
	Character	*64	rse_file
	Character	*128	rse(16)
	Integer		*2 	cats_num
	integer		*2	old_num, new_num
	Integer		*4	status
	Integer		*2	ct_status(20)
	Integer		*4	i,j

	Integer		*4 	CCT_Query_TTG_Catalog
	Integer		*4 	FUT_Read_RSE

	External	FUT_NoCatRecs
	External	FUT_CTQueryCat

C
C Query the archive catalog for record selection expressions that
C are in the specified timerange.
C
	ttg_query.archive_id = 'CSDR$FIRAS_REF'
	ttg_query.start_time(1) = jstart(1)
	ttg_query.start_time(2) = jstart(2)
	ttg_query.stop_time(1) = jstop(1)
	ttg_query.stop_time(2) = jstop(2)
	ttg_query.dataset_name = dataset

	status = CCT_Query_TTG_Catalog(ttg_query,old_cats,cats_num)

	If (status) Then

	   old_num = cats_num
	   j = 0

	   Do i = 1, old_num
	      If (.not. old_cats(i).deletion_flag) Then
		 j = j+1
		 new_cats(j) = old_cats(i)
	      End If
	   End Do

	   new_num = j

	   If (new_num .ne. 0) Then

	      final_cats = new_cats(1)

	      Do j = 2, new_num
	         If (new_cats(j).version_number .gt. final_cats.version_number)
     .		    Then
		    final_cats = new_cats(j)
		 End If
	      End Do

	      i = index(final_cats.filename_extension,' ') - 1
	      If (i .Eq. -1) Then
		 i = Len(final_cats.filename_extension)
	      End If

	      rse_file = 'CSDR$FIRAS_REF:' // dataset // '.' //
     .			 final_cats.filename_extension(1:i) // ';' //
     .			 final_cats.version_number

	      status = FUT_Read_RSE(rse_file, rse)

	   Else
	      status = %Loc(FUT_NoCatRecs)
	      Call LIB$Signal(FUT_NoCatRecs)
	   End If

	Else
	  FUT_Query_RSE = %Loc(FUT_CTQueryCat)
          Call LIB$Signal(FUT_CTQueryCat,%Val(1),%Val(status))
	End If


	FUT_Query_RSE = status

	Return
	End
