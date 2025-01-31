C----------------------------------------------------------------------------

	Integer*4 Function FTB_Refine_Overlap(Cats, NRecs,
	1				      FPP_JStart, FPP_JStop,
	2				      NTimes)

C-----------------------------------------------------------------------------
C
C	Purpose: Refines the non-overlapping times from each raw science
C		 channel to a single list of channel independent times
C	         on which FPP will be run.
C
C       Author:  R. Kummerer / STX, October 16, 1989
C
C   	Invocation:  status = FTB_Refine_Overlap(Cats, NRecs,
C						 FPP_JStart, FPP_JStop,
C						 NTimes)
C
CH	Change Log:
CH
C-----------------------------------------------------------------------------------
C
C       Input Files :
C
C	Output Files :
C
C  	Input Parameters:
C     
C	  Name               Type	Description
C         ---------------------------------------------------------------------
C	  Cats(384,4)	     Record	List of channel dependent non-overlap
C					times.
C	  NRecs(4)	     I*2	Number of catalog records.
C-------------------------------------------------------------------------------------
C
C	Output Parameters:
C	  Name            Type		Description
C         ---------------------------------------------------------------------
C	  FPP_JStart(384,2)  I*4        Selected start time
C	  FPP_JStop(384,2)   I*4        Selected start time
C	  NTimes	     I*2	Number of FPP invocation times.
C------------------------------------------------------------------------------------
C
C	 Subroutines And Functions Used :
C	 ------------------------------
C
C	 Common Variables Used :
C
C	 Name		Type	 Use		Description
C------------------------------------------------------------------------------
C	 Include Files :
C	 ---------------
C---------------------------------------------------------------------------------
C
C	 Processing Method : PDL for FTB_Refine_Overlap
C
C			PDL for FTB_Refine_Overlap
C			--------------------------
C
C-----------------------------------------------------------------------------------

	Implicit None

C	Include Files  

	Include	'(CUT_Params)'

C       Dictionaries

	Dictionary 'CCM_CME_Catalog_Entry'
 
C       Records

C       Passed parameters.

	Record /CCM_CME_CATALOG_ENTRY/ Cats(CUT$_Max_Segments,4)
	Integer*2	NRecs(4)

C	Out Parameters.

	Integer*4	FPP_JStart(2,CUT$_Max_Segments)
	Integer*4	FPP_JStop(2,CUT$_Max_Segments)
	Integer*2	NTimes

C  	Used as Function Values

	Logical*2	Time_GT
	Logical*2	Time_LT
	Logical*2	Time_GE
	Logical*2	Time_LE

C       Local Variables

        Integer*4	Status  	! Function return status
        Integer*4	Chan    	! Index for channel number
        Integer*4	Chk_Chan    	! Index for checked channel
        Integer*4	I
        Integer*4	Cur_Seg
        Integer*4	Ref_Seg
	Logical*1	Matched
	Logical*1	Match(CUT$_Max_Segments,4)
	Integer*4	Seg_Start(2,4)
	Integer*4	Seg_Stop(2,4)
        Integer*4	Time_Count
        Integer*4	Chk_Time
        Integer*4	overlap_count

C       External used

        External	FTB_Normal
        External	FTB_FPPOvrLap
        
C Successively loop over each channel, know as the reference channel, and
C compare each segment extension within that channel with all in the remaining
C channels.  Attempt to find a segment that matches the reference segment in
C each of the remaining channels.  Mark those segments that have identical
C extensions as matching segments and save the time range of each matched
C segment.  Once a set of matched segments are found, their time ranges
C are searched for minimum start and maximum stop times.  Such a minimum start
C and maximum stop time pair is a time range on which FPP may be run.

C Initialize flags, counters and the match array.

	Do Chan=1,4
	   Do Cur_Seg=1,CUT$_Max_Segments
	      Match(Cur_Seg,Chan) = .False.
	   Enddo
	Enddo

	NTimes = 0

C Loop over the reference channels.

	Do Chan=1,4

C Attempt to match each segment in the the reference channel with a segment
C in any one of the remaining channels.

	   Do Ref_Seg=1,NRecs(Chan)

C Check if the reference segment has been already matched; if not, proceed.

	      If ( .Not. Match(Ref_Seg,Chan) ) Then

	         Match(Ref_Seg,Chan) = .True.

C An unmatched reference segment consists of one time range on which FPP
C could be run.  Initialize the time range buffer ("matched segments list").

	         Time_Count = 1
	         Seg_Start(1,Time_Count) = Cats(Ref_Seg,Chan).Initial_Time(1)
	         Seg_Start(2,Time_Count) = Cats(Ref_Seg,Chan).Initial_Time(2)
	         Seg_Stop(1,Time_Count)  = Cats(Ref_Seg,Chan).Final_Time(1)
	         Seg_Stop(2,Time_Count)  = Cats(Ref_Seg,Chan).Final_Time(2)

C Loop over the remaining channels.

	   	 Do Chk_Chan=Chan+1,4

		    Matched = .False.
	            Cur_Seg = 1

C Compare the reference segment with each segment in the current channel iff
C they had not already been previously matched.

	            Do While ( (.Not. Matched) .And.
	1		       (Cur_Seg .Le. NRecs(Chan)) )

		       If ( .Not. Match(Cur_Seg,Chk_Chan) ) Then

C If a match is found, add the time range of the just matched segment to
C the time range buffer and mark the segment as matched.

		          If ( Cats(Ref_Seg,Chan).Filename_Extension .Eq.
	1		       Cats(Cur_Seg,Chk_Chan).Filename_Extension ) Then
		             Matched = .True.
		             Match(Cur_Seg,Chk_Chan) = .True.
		             Time_Count = Time_Count + 1
		             Seg_Start(1,Time_Count) =
	1				Cats(Cur_Seg,Chk_Chan).Initial_Time(1)
		             Seg_Start(2,Time_Count) = 
	1				Cats(Cur_Seg,Chk_Chan).Initial_Time(2)
		             Seg_Stop(1,Time_Count)  =
	1				Cats(Cur_Seg,Chk_Chan).Final_Time(1)
		             Seg_Stop(2,Time_Count)  =
	1				Cats(Cur_Seg,Chk_Chan).Final_Time(2)
		          Endif		! If (Cats(Ref_Seg,Chan).File ...

		       Endif		! If ( .Not. Match ) ...

		       Cur_Seg = Cur_Seg + 1

                    Enddo		! Do While ...

	         Enddo		! Do Chk_Chan=Chan+1,4

C All remaining channel segments have been checked against the reference
C channel segment.  A time range buffer of matched segments has been built
C and can now be searched for minimum start and maximum stop times on which
C FPP can be run.  Add that FPP time range to the FPP time ranges run list.

		 NTimes = NTimes + 1

		 FPP_JStart(1,NTimes) = Seg_Start(1,1)
		 FPP_JStart(2,NTimes) = Seg_Start(2,1)
		 FPP_JStop(1,NTimes)  = Seg_Stop(1,1)
		 FPP_JStop(2,NTimes)  = Seg_Stop(2,1)

		 Do Chk_Time=2,Time_Count
		    If (Time_Lt(Seg_Start(1,Chk_Time),FPP_JStart(1,NTimes))) Then
		       FPP_JStart(1,NTimes) = Seg_Start(1,Chk_Time)
		       FPP_JStart(2,NTimes) = Seg_Start(2,Chk_Time)
		    Endif
		    If (Time_Gt(Seg_Stop(1,Chk_Time),FPP_JStop(1,NTimes))) Then
		       FPP_JStop(1,NTimes) = Seg_Stop(1,Chk_Time)
		       FPP_JStop(2,NTimes) = Seg_Stop(2,Chk_Time)
		    Endif
		 Enddo

	      Endif		! If ( .Not. Match(Ref_Seg,Chan) ) Then

	   Enddo		! Do Ref_Seg=1,NRecs(Chan)

	Enddo			! Do Chan=1,3

C As a precaution, check the candidate FPP run times for overlap.

	overlap_count = 0

	Do Chk_Time=1,NTimes-1

	   If (Time_GE(FPP_JStop(1,Chk_Time),FPP_JStart(1,Chk_Time+1)) .And.
	1      Time_LE(FPP_JStop(1,Chk_Time),FPP_JStop(1,Chk_Time+1))) Then

	      overlap_count = overlap_count + 1

	   Endif

	Enddo

	If (overlap_count .Ne. 0) Then
	   Call LIB$Signal(FTB_FPPOvrLap,%Val(1),%Val(overlap_count))
	Endif


        FTB_Refine_Overlap = %Loc(FTB_normal)

	Return
	End   
