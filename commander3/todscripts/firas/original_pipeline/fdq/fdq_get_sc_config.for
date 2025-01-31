C-------------------------------------------------------------------------------

	Integer*4 Function FDQ_Get_SC_Config ( SC_Side )

C-------------------------------------------------------------------------------
C
C	Purpose: To determine which side of th spacecraft is controlling
C	         FIRAS and return the side to the caller.
C
C	Author: Shirley M. Read
C		STX, May 1988
C
C	Invocation:
C	        Status = FDQ_Get_Sc_Config ( SC_Side )
C
C	Modificaton History:
C
C	  Author	    Date	  Modification
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  SC_Side       I*2             Spacecraft side 1 or 2
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  SC_Side       I*2             Spacecraft side 1 or 2
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
C	     Until the spacecraft archive is available to the FIRAS pipeline, 
C	   the user may select the spacecraft side on the input command line.
C	   The default will be side 1 since this side is currently scheduled 
C	   to power FIRAS for the I & T. When the decommutated spacecraft 
C	   archive becomes available, the archive will be read using the 
C	   average time of the current science set as a key to obtain the 
C	   currect spacecraft side.
C	     The information needed to identify physical thermistors and
C	   other FIRAS components actually read in telemetry will be
C	   decommutated from Subcom 34. There will be two bytes worth of data
C	   per major frame. This information will determine the spacecraft
C	   telemetry unit that was selected. 
C	     Since the spacecraft archive is not yet designed, the real time
C	   database words are used in this documentation to aid in the
C	   future implementation of this FDQ function. The database words
C	   contain the relay status bits for the telemetry units on the 
C	   spacecraft. The words are the following:
C
C	     CTU1PWR -- Central Telemetry Unit 1 - power for spacecraft side 1
C	     CTU2PWR -- Central Telemetry Unit 2 - power for spacecraft side 2
C	     CITU1PWR -- Central Instrument Telemetry Unit 1 - power for 
C		         spacecraft side 1
C            CITU2PWR -- Central Instrument Telemetry Unit 2 - power for
C		         spacecraft side 2
C	     CSTU1PWR -- S/C Telemetry Unit 1 -- power for spacecraft side 1,
C			 needed for FIRAS LMACs only
C	     CSTU2PWR -- S/C Telemetry Unit 2 -- power for spacecraft side 2,
C			 needed for FIRAS LMACs only
C	     CICU1PWR -- Central Instrument Command Unit 1 - power for command
C	                 unit only for spacecraft side 1
C            CICU2PWR -- Central Instrument Command Unit 2 - power for command
C			 unit only for spacecraft side 2
C            
C	   Algorithm for determining the side:
C
C	    1. Check the status bits for the CTU1PWR and CTU2PWR first.
C	       If only one is set, then this side is powering FIRAS.
C	       If neither is set, then signal an error.
C	    2. If both are set, there may be a problem. In this case, the 
C	       status bit for CITU1PWR and CITU2PWR could decide if only
C	       one side is really on. 
C           3. As a last resort, the CSTU1PWR and CSTU2PWR relating to the 
C	       LMACs could be checked. The final decision in this case is TBD.
C
C	     The command units can not be used for determining the spacecraft
C	     side. They are listed here in order to describe the related 
C	     telemetry information in Subcom 34.
C
C------------------------------------------------------------------------------
C
	Implicit None

!       External Parameters.
	
	External   	FDQ_NORMAL
	External	FDQ_ABERR

!	Input/Output Passed Parameters

	Integer*2       SC_Side         ! Spacecraft side powering FIRAS

!	Local Varaibles.

	Integer*4	Status, Retstat		! Status variables	
	Integer*4       Success / 1 /, Error / 2 /
	Integer*2       Zero / 0 /, One / 1/ 
	Logical*1       First / .true. /        ! First time called.

!       Set return status to success.

	retstat = success

!	If the user has not set the spacecraft side, set it to the default.

	If ( First ) Then
	   First = .false.
	   If ( Sc_Side .eq. Zero ) Sc_Side = One
	Endif

	If (retstat.eq.success) Then
	  FDQ_Get_SC_Config = %loc(FDQ_NORMAL)
	Else
	  FDQ_Get_SC_Config = %loc(FDQ_ABERR)
	Endif

	Return
	End
	
