C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Accum_Eng_Table ( Fieldval,
	1	  Dom_Mode, Nbins, Bins, Eng_Table )

C-------------------------------------------------------------------------------
C
C	Purpose: To accumulate the record counts for each scan mode in the
C	         appropriate engineering bin.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Accum_Eng_Table ( Fieldval, Dom_Mode, Nbins,
C			     Bins, Eng_Table )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		  Type		Description
C 	  ----------------------------------------------------------------------
C	  Fieldval        R*4           FIRAS engineering field value
C	  Dom_Mode        I*2           Predominant MTM scan mode on the record
C	  Nbins           I*2           Number of engineering bins
C	  Bins(5)         R*4           Engineering value bins
C	
C	Output Parameters:
C	  Name		  Type		Description
C 	  ----------------------------------------------------------------------
C	  Eng_Table(6,4)  I*4           Counts of number of records for each
C				        engineering value range by scan mode
C	
C	Subroutines Called:
C
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C	  Set function return.
C	  Increment the record count in the proper bin for engineering range
C	    and scan mode.
C	  If the value is above the top bin value, increment the bin above.
C	  If the count has not been incremented for any bin, signal an error
C	    and set the function return to error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

C	Passed Parameters.

	Real *4      Fieldval        ! Engineering field value
	Integer*2    Dom_Mode        ! Predominat scan mode
	Integer*2    Nbins           ! Number of bins
	Real*4       Bins(5)         ! Ranges for engineering values
	Integer*4    Eng_Table(6,4)  ! Table of record counts: array of 
				     ! Bin value versus scan mode
C	Local Declarations.

	Logical*1    Entered	     ! Count has been incremented for record
	Integer*2    Ix		     ! Index
        Integer*2    Ind             ! Index
C	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NoBinMatch

	FTB_Accum_Eng_Table = %loc(FTB_Normal)

C	Increment the record counter in the proper bin.
        Ind = Dom_Mode + 1
	Entered = .False.
	Ix = 1
	Do While ((.NOT. Entered) .AND. (Ix .LE. Nbins))
	   If ( Fieldval .LE. Bins(Ix) ) Then
	     Eng_Table(Ix,Ind) = Eng_Table(Ix,Ind) + 1
	     Entered = .True.
	   Endif
           Ix = Ix + 1
	Enddo 

	If ((.NOT. Entered) .AND. ( Fieldval .GT. Bins(Nbins))) Then
	   Entered = .True.
	   Ix = Nbins + 1
	   Eng_Table(Ix,Ind +1) = Eng_Table(Ix,Ind) + 1
	Endif
	If ( .NOT. Entered ) Then
	   Call Lib$Signal(FTB_NoBinMatch)
	   FTB_Accum_Eng_Table = %loc(FTB_Aberr)
        Endif

	Return
	End
