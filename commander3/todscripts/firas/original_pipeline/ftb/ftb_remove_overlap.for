C------------------------------------------------------------------------

      Integer*4 Function FTB_Remove_Overlap ( Cats, NRecs, Report, Report_LUN )

C------------------------------------------------------------------------
C    PURPOSE: Remove overlaps in the initial and final times found in
C	      catalog records input to FPP.  The initial and final times
C	      are adjusted to remove overlaps.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            July 31, 1989
C
C    CHANGE:
C           SPR 7020, FTB_OVERLAP CT read error during period:90175-90176
C           H.WANG, STX, 7/2/90
C
C           SPR 7346, FTB_OVERLAP CT read error during period:90230-90232
C             FTB_OVERLAP should handle the case:
C               |-----------|  < A segment delivered after B
C               |--------------| < B segment
C                 from A      B
C               |-----------||-| < create 2 new segments
C           H.WANG, STX, 8/22/90
C
C    INVOCATION: 
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: 
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES: 
C
C----------------------------------------------------------------------

	Implicit None

	Include	'(FUT_Params)'
	Include	'(CUT_Params)'

	Integer		*2	NRecs
	Logical		*1	Report
	Integer		*4	Report_LUN
	Integer		*4	status

	Dictionary 'CCM_CME_Catalog_Entry'
	Record /CCM_CME_Catalog_Entry/ Cats(CUT$_Max_Segments)

	Integer		*4	lst_rec		!Last record = potential over-
						!lapped record.
	Integer		*4	cur_rec		!Current record = record to be
						!checked for overlaps with last
						!record.
	Integer		*4	DeltaTime(2)	!Prevent start and stop times
						!of overlapping segments from
						!being exactly the same.
        Integer         *4      temp_s(2) 
        Integer         *4      temp_e(2) 
	Integer		*4	SYS$BinTim
	Integer		*4	LIB$Addx
	Integer		*4	LIB$Subx
	Integer		*4	FUT_Clean_Catalog
	Logical		*2	Time_GE		!CT time comparsion routines
	Logical		*2	Time_GT
	Logical		*2	Time_LE
	Logical		*2	Time_EQ
	Logical		*2	Time_LT
	Integer		*4	FTB_Sort_Catalog
	Integer		*4	FTB_Overlap_Set_Start
	Integer		*4	FTB_Overlap_Set_Stop

	External	FTB_Normal

C Initialize.

	status = SYS$BinTim ( '0 0:0:0.100', DeltaTime )

C Remove overlaps by comparing initial and final times of successive catalog
C records.  At this point the catalog records are sorted into ascending
C initial time and then by ascending activity time.

	status = FUT_Clean_Catalog ( Cats, NRecs, fac_present )

	status = FTB_Sort_Catalog ( Cats, NRecs )

	cur_rec = 2

	Do While (cur_rec .Le. NRecs)

	   lst_rec  = cur_rec - 1

C Adjust current record's initial time to exclude it's overlapping portion
C if the last record's final time is bracketted by the current record.

	   If ( Time_GE ( Cats(lst_rec).Initial_Time,
	1		  Cats(cur_rec).Initial_Time )  .And.
	2	Time_LE ( Cats(lst_rec).Final_Time,
	3		  Cats(cur_rec).Final_Time ) )  Then

C Completely overlapped segment.  This segment will be removed from the
C catalog list by the latter FUT_CLEAN_CATALOG call (FUT_CLEAN_CATALOG
C removes catalog records with the delete flag set).
C
C		|--------------------------------|  < A delivered after B
C			|---------------|           < B 
C	Time ...............................................>
C
C Segment A completely overlaps B.  Segment B entirely superceded by
C segment A; segment B, R.I.P..  Marking segment B as deleted here
C will cause the latter FUT_CLEAN_CATALOG call to remove it from the
C catalog list.

	      Cats(lst_rec).Deletion_Flag = .True.
	      Cats(lst_rec).Initial_Time(1) = -1
	      Cats(lst_rec).Initial_Time(2) = -1
	      Cats(lst_rec).Final_Time(1) = -1
	      Cats(lst_rec).Final_Time(2) = -1

	   Else If ( Time_GE ( Cats(lst_rec).Final_Time,
	1		       Cats(cur_rec).Initial_Time )  .And.
	2	     Time_LE ( Cats(lst_rec).Final_Time,
	3		       Cats(cur_rec).Final_Time )) then
C Bracketted final time.
C
C		        1	     2
C		|----------------||------|       < Create segments
C			|---------------|       < A delivered after B 
C		|----------------|		< B
C	Time ...............................>
C
C Segment A partially overlaps segment B.  One new segment is created.
C Segment 1 is the new segment, segment 2 is the original segment A.
              If(time_lt(cats(lst_rec).final_time,cats(cur_rec).final_time))
	1	then
               Cats(cur_rec).initial_time(1) = Cats(lst_rec).final_time(1) 
               Cats(cur_rec).initial_time(2) = Cats(lst_rec).final_time(2) 
	      status = LIB$subx ( Cats(cur_rec).Initial_Time,
	1			  DeltaTime,
	2			  Cats(cur_rec).Initial_Time )

	      status = FTB_Overlap_Set_Start ( Cats(cur_rec), Report,
	1				       Report_LUN )

C Delta times are negative.  Therefore, LIB$Subx performs an addition.
            else
               Cats(lst_rec).final_time(1) = Cats(cur_rec).initial_time(1) 
               Cats(lst_rec).final_time(2) = Cats(cur_rec).initial_time(2) 
	      status = LIB$addx ( Cats(lst_rec).final_Time,
	1			  DeltaTime,
	2			  Cats(lst_rec).final_Time )

	      status = FTB_Overlap_Set_Stop ( Cats(lst_rec), Report,
	1				       Report_LUN )

                
            endif
C Delta times are negative.  Therefore, LIB$Addx performs a subtraction.

 

	   Else If ( Time_LT ( Cats(lst_rec).Final_Time,
	1		       Cats(cur_rec).Initial_Time ) ) Then

C Segments do not overlap in any way.

	      Continue

	   Else If ( Time_GT ( Cats(lst_rec).Final_Time,
	1		       Cats(cur_rec).Final_Time ) ) Then

C More recent segment does not completely overlap latter segment final time. 
C Create a new segment.
C
C		    1		2	     3
C		|------||---------------||-------|  < Create segments
C			|---------------|           < A delivered after B
C		|--------------------------------|  < B
C	Time ...............................................>
C
C Segment A, delivered after segment B, overlaps a middle portion of B.
C The result are two "new segments", 1 and 3 in the diagram above.
C Segment 2, the original segment A, remains unchanged.

	     If ( Time_eq ( Cats(lst_rec).initial_Time,
	1		       Cats(cur_rec).initial_Time ) ) Then
C               |--------------------------------||---| create new segment
C		|--------------------------------|      < A delivered after B
C	        |-------------------------------------|           < B 
C	Time ...............................................>

              temp_s(1) = cats(lst_rec).initial_time(1)
              temp_s(2) = cats(lst_rec).initial_time(2)
              temp_e(1) = cats(lst_rec).final_time(1)
              temp_e(2) = cats(lst_rec).final_time(2)

	      Cats(lst_rec).Final_Time(1) =
	1			Cats(cur_rec).Final_Time(1)
	      Cats(lst_rec).final_Time(2) =
	1			Cats(cur_rec).final_Time(2)
	      status = LIB$addx ( Cats(cur_rec).final_Time,
	1			  DeltaTime,
	2			  Cats(cur_rec).Initial_Time )
	      Cats(cur_rec).Final_Time(1) = temp_e(1)
	      Cats(cur_rec).Final_Time(2) = temp_e(2)

	      status = FTB_Overlap_Set_Start ( Cats(cur_rec), Report,
	1				       Report_LUN )
	1			
            else
	      NRecs = NRecs + 1

C Create new segment (3).

	      Cats(NRecs) = Cats(lst_rec)

	      Cats(NRecs).Initial_Time(1) =
	1			Cats(cur_rec).Final_Time(1)
	      Cats(NRecs).Initial_Time(2) =
	1			Cats(cur_rec).Final_Time(2)
	      Cats(NRecs).Final_Time(1) =
	1			Cats(lst_rec).Final_Time(1)
	      Cats(NRecs).Final_Time(2) =
	1			Cats(lst_rec).Final_Time(2)

C Delta times are negative.  Therefore, LIB$Subx performs an addition.

	      status = LIB$Subx ( Cats(NRecs).Initial_Time,
	1			  DeltaTime,
	2			  Cats(NRecs).Initial_Time )

	      status = FTB_Overlap_Set_Start ( Cats(NRecs), Report,
	1				       Report_LUN )

C Create new segment (1).

	      Cats(lst_rec).Final_Time(1) =
	1			Cats(cur_rec).Initial_Time(1)
	      Cats(lst_rec).Final_Time(2) =
	1			Cats(cur_rec).Initial_Time(2)

C Delta times are negative.  Therefore, LIB$Addx performs a subtraction.

	      status = LIB$Addx ( Cats(lst_rec).Final_Time,
	1			  DeltaTime,
	2			  Cats(lst_rec).Final_Time )

	      status = FTB_Overlap_Set_Stop ( Cats(lst_rec), Report,
	1				      Report_LUN )

C Re-sort the catalog list so the new segment is properly inserted and
C start the search over again.

	      status = FTB_Sort_Catalog ( Cats, NRecs )

	      cur_rec = 1
            endif 
	   End If

	   cur_rec = cur_rec + 1

	End Do

C Remove completely overlapped segments (those just marked as deleted) from
C the catalog list.

	status = FUT_Clean_Catalog ( Cats, NRecs, fac_present )


	FTB_Remove_Overlap = %Loc(FTB_Normal)

	Return
	End
