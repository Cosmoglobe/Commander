C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Dump_Hkp_Time ( Gmt, Numrec, 
	1	  Time_Flag, Rpt_Lun)

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS housekeeping records
C		 with invalid time tags. The records bracketing the bad time 
C		 periods are also included. 
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Dump_Hkp_Time 
C			     ( Gmt, Numrec, Time_Flag, Rpt_Lun )
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
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Gmt           C*14            Time tag of record
C	  Numrec         I*4            Number of records with repeated time tag
C	  Time_Flag      C*2            Flag indicating reason for dump
C	  Rpt_Lun        I*4            Report unit number
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
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
C	  Set function return status.
C	  Load the Gmt values into a string with the print format.
C	  If first time, print the table header.
C	  Dump the record information.
C	  If status for write is bad, reset the function return status and
C	    signal an error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

C	Passed Parameters.

	Character*14 Gmt             ! Time tag of record
	Integer*4 Numrec             ! Number of records with repeated time tag
	Character*2 Time_Flag        ! Flag indicating reason for dump
	Integer*4 Rpt_Lun            ! Report unit number

C	Local Declarations.

	Logical*1 First / .True. /
	Integer*4 Rstatus
	Character*19 Gmt_Char 
	Integer*4 Zero / 0 /
        Character*1 Dash(130)/ 130 * '_'/

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr

	FTB_Dump_Hkp_Time = %loc(FTB_Normal)

	If ( First ) Then
	  Write ( Unit=Rpt_lun, FMT=100, Iostat=Rstatus )
 100      Format (4x,'Transmit_Time',4x,'Num.Rec.',1x,'Flag')
	  Write ( Unit=Rpt_lun, FMT=150, Iostat=Rstatus ) Dash
 150      Format(1x, 130a)
          First = .False.
	  If (Rstatus .NE. Zero ) Then
	    FTB_Dump_Hkp_Time = %loc(FTB_Aberr)
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Gmt_Char= Gmt(1:2)//':'//Gmt(3:5)//':'//Gmt(6:7)//':'//
	1	  Gmt(8:9)//':'//Gmt(10:11)//'.'//Gmt(12:14)

	If ( FTB_Dump_Hkp_Time .EQ. %loc(FTB_Normal)) Then
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Gmt_Char, 
	1     Numrec, Time_Flag
 200      Format(1x,a,1x,i8,2x,a)
	  If ( Rstatus .NE. Zero ) Then
!	    Allow print even if overflow
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Return
	End
