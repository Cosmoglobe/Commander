C-----------------------------------------------------------------------------
	Integer*4 Function FPP_Check_Segment_Overlap ( Proc_Chan, Num_Chan,
	1		Check_JStart, Check_JStop, Report, LUN_Rpt )
C-----------------------------------------------------------------------------
C	Purpose: To search for overlapping FPP_SDF datasets.
C
C       Author:  R. Kummerer / STX, August 25, 1989
C
C   	Invocation:  Status = FPP_Check_Segment_Overlap ( Proc_Num, Num_Chan,
C			Check_JStart, Check_JStop, Report, LUN_Rpt )
C
CH	Change Log:
CH
CH		Version 4.4.1, 09/04/89, SPR 4480, R. Kummerer, STX
CH			Limit overlap check timerange to four week window
CH			centered on segment being processed.
CH
CH		Version 4.4.04, 11/29/89, SPR 5206, R. Kummerer, STX
CH			Prevent call to FUT_CLEAN_CATALOG when NRECS is 0.
CH
CH		Version New, 4/9/91, New Requirements, Larry P. Rosen, STX
CH			Modify report format.
C-----------------------------------------------------------------------------
C  	Input Parameters:
C     
C	  Name           Type		Description
C         ---------------------------------------------------------------------
C         Proc_Chan(4)	  I*2		Processed channel number list
C	  Num_Chan	  I*2		Number of channels processed
C         Report          L*1		Flag for printing a report file or not
C         LUN_Rpt         I*4		Logical unit for report file
C	  Check_JStart(2) I*4		Overlap check start time
C	  Check_JStop(2)  I*4		Overlap check stop time
C-----------------------------------------------------------------------------
C	 Subroutines And Functions Used :
C	 ------------------------------
C        CCT_QUERY_TOD_CATALOG
C        FUT_CLEAN_CATALOG
C        LIB$SIGNAL
C-----------------------------------------------------------------------------
C	 Include Files :
C	 ---------------
C        FPP_MSG			External message params needed for FPP
C        CT$LIBRARY:CTUSER:INC	        COBETRIEVE return status definitions
C        $SSDEF				System status values
C	 FUT_PARAMS
C        CCT_QUERY_TOD_CATALOG_RECORD
C-----------------------------------------------------------------------------
C	 Processing Method : PDL for FPP_Check_Segment_Overlap
C
C			PDL for FPP_Check_Segment_Overlap
C			---------------------------------
C	 DO FOR each channel xx selected
C	   CALL CT_QUERY_TOD_CATALOG to fetch all FPP_SDF_xx catalog records
C		in the timerange 86001000000000 to 99365235959990
C	   DO FOR 2 to total number of FPP_SDF_xx catalog records
C	     IF (current FPP_SDF_xx catalog record overlaps
C		 last    FPP_SDF_xx catalog record) THEN
C		INCREMENT overlap counter
C		IF (writting a report) THEN
C		   DISPLAY the overlapping segments in the report
C		ENDIF
C	     ENDIF
C	   ENDDO
C	   IF (overlap counter <> 0) THEN
C	     CALL LIB$SIGNAL to issue an overlap warning
C	   ENDIF
C	 ENDDO
C-----------------------------------------------------------------------------
	Implicit None

C Include Files  

	Include	'($SSDef)'
	Include 'CT$Library:CTUser.Inc'
        Include '(CCT_Query_TOD_Catalog_Record)'
        Include '(FUT_Params)'
	Include	'(CUT_Params)'
        Include '(FPP_Msg)'

C Dictionaries
 
	Dictionary 'CCM_CME_Catalog_Entry'

C Records

	Record /QUERY_TOD_CATALOG/ Cat_Query
	Record /CCM_CME_CATALOG_ENTRY/ Catalog_Recs(CUT$_Max_Segments)

C Passed parameters.

	Integer*2    Proc_Chan(4) ! Channels processed list
	Integer*2    Num_Chan     ! Number of channels processed
	Logical*1    Report       ! Flag for printing output report
	Integer*4    Lun_Rpt      ! Logical unit for output report
	Integer*4    Check_JStart(2) ! Overlap check start time
	Integer*4    Check_JStop(2)  ! Overlap check start time

C Used as Function Values

	Integer*4    CCT_Query_Tod_Catalog
	Integer*4    CCT_Get_FileSpec_Fields
	Logical*1    Time_EQ
	Logical*1    Time_LT
	Logical*1    Time_GT
	Integer*4    FUT_Clean_Catalog

C Local Variables

	Integer*2    Chan_Num   ! Channel number
	Integer*2    Chan
	Integer*2    NRecs
	Integer*4    Status
	Integer*4    i
	Integer*4    j
	Integer*4    k
	Integer*4    l
	Character*32 Ext
	Character*64 First_File
	Character*64 Second_File
	Character*14 First_Initial
	Character*14 First_Final
	Character*14 Second_Initial
	Character*14 Second_Final
	Integer*4    Overlap_Count

	Character*1  Dashes(79)  / 79 * '-' /
	Character*79 Dash 
	Equivalence  ( Dashes(1), Dash )
	Integer*4    Iostatus	! Return IO status

C Check FPP_SDF_xx segments for overlap.

        FPP_Check_Segment_Overlap = %Loc(FPP_Normal)

	Do Chan = 1, Num_Chan

	   Chan_Num = Proc_Chan(Chan)

C Query the archive catalog for raw science segments that
C are in the standard CT data timerange.

	   Cat_Query.Archive_Id = 'CSDR$FIRAS_RAW'
	   Cat_Query.Dataset_Name = 'FPP_SDF_' // FAC_Channel_Ids(Chan_Num)
	   Cat_Query.Start_Time(1) = Check_JStart(1)
	   Cat_Query.Start_Time(2) = Check_JStart(2)
	   Cat_Query.Stop_Time(1) = Check_JStop(1)
  	   Cat_Query.Stop_Time(2) = Check_JStop(2)

	   Status = CCT_Query_TOD_Catalog ( Cat_Query, Catalog_Recs, NRecs )

	   If (Status .EQ. 0) Then
	      NRecs = 0
              FPP_Check_Segment_Overlap = %Loc(FPP_CTQueryCat)
              Call LIB$Signal ( FPP_CTQueryCat, %Val(2), %Val(Chan_Num),
	1			%Val(Status) )
	   Endif

C Remove duplicates and deleted segments from the list just queried.

	   If (NRecs .Gt. 0) Then

	      Status = FUT_Clean_Catalog ( Catalog_Recs, NRecs, FAC_Present )

	      If (Status .EQ. 0) Then
	         NRecs = 0
                 FPP_Check_Segment_Overlap = %Loc(FPP_ClnCat)
                 Call LIB$Signal ( FPP_ClnCat, %Val(1), %Val(Status) )
	      Endif

	   Endif

C Perform the overlap check using the catalog records just cleaned.

	   Overlap_Count = 0

	   Do i=1,NRecs-1

	      Do j=i+1,NRecs

	         If (Time_Eq(Catalog_Recs(i).Initial_Time,
	1			Catalog_Recs(j).Initial_Time)
	1		.Or.
	1	     Time_Eq(Catalog_Recs(i).Final_Time,
	1			Catalog_Recs(j).Final_Time)
	1		.Or.
	1	     Time_Eq(Catalog_Recs(i).Initial_Time,
	1			Catalog_Recs(j).Final_Time)
	1		.Or.
	1	     Time_Eq(Catalog_Recs(i).Final_Time,
	1			Catalog_Recs(j).Initial_Time)
	1		.Or.
	1	     (Time_Lt(Catalog_Recs(i).Initial_Time,
	1			Catalog_Recs(j).Initial_Time)
	1		.And.
	1	      Time_Gt(Catalog_Recs(i).Final_Time,
	1			Catalog_Recs(j).Initial_Time))
	1		.Or.
	1	     (Time_Gt(Catalog_Recs(i).Initial_Time,
	1			Catalog_Recs(j).Initial_Time)
	1		.And.
	1	      Time_Lt(Catalog_Recs(i).Initial_Time,
	1			Catalog_Recs(j).Final_Time)) )
	1           Then

C Form overlapping raw science segment catalog filenames.

	               k = Index(Catalog_Recs(i).Dataset_Name,' ')
	               l = Index(Catalog_Recs(i).Filename_Extension,' ')
	               Ext = Catalog_Recs(i).Filename_Extension(1:l-1) // ';' //
	1				Catalog_Recs(i).Version_Number
	               First_File = Catalog_Recs(i).Dataset_Name(1:k-1) //
	1				    '.' // Ext

		       Call CT_Binary_To_GMT(Catalog_Recs(i).Initial_Time,
	1					First_Initial)
		       Call CT_Binary_To_GMT(Catalog_Recs(i).Final_Time,
	1					First_Final)

	               k = Index(Catalog_Recs(j).Dataset_Name,' ')
	               l = Index(Catalog_Recs(j).Filename_Extension,' ')
	               Ext = Catalog_Recs(j).Filename_Extension(1:l-1) // ';' //
	1				Catalog_Recs(j).Version_Number
	               Second_File = Catalog_Recs(j).Dataset_Name(1:k-1) //
	1				    '.' // Ext

	               k = Index(Second_File,' ')
	               l = Index(First_File,' ')

		       Call CT_Binary_To_GMT(Catalog_Recs(j).Initial_Time,
	1					Second_Initial)
		       Call CT_Binary_To_GMT(Catalog_Recs(j).Final_Time,
	1					Second_Final)

	               If (Report) Then
	                  If (Overlap_Count .EQ. 0) Write (Unit=Lun_Rpt,
	1                    FMT=300, Iostat=Iostatus) FAC_Channel_Ids(Chan_Num)
	                  Write (LUN_Rpt, 600) First_File(1:l-1),
	1                    First_Initial, First_Final, Second_File(1:k-1),
	1                    Second_Initial, Second_Final
	                  EndIf
	              Overlap_Count = Overlap_Count + 1
	         Endif
	      Enddo
	   Enddo

	   If (Overlap_Count .Gt. 0) Call LIB$Signal
	1	(FPP_Overlap, %Val(1), FAC_Channel_Ids(Chan_Num))

	Enddo

300	Format (1x,/,20x, 'FPP segment overlaps encountered in channel ', A,':')
600	Format (1x,/, 5x, 'Segment ', A, ' (', A, '-', A,') overlaps', /,6x,
	1                  'segment ', A, ' (', A, '-', A, ')', /)

	Return
	End
