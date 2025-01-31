C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Dump_Eng_Time ( Channel, Gmt, Status_Bits,
	1	  Transmit_Frame, Numrec, Time_Flag, Rpt_Lun)

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS engineering records with
C		 invalid time tags. The records bracketing the bad time periods
C		 are also included. 
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Dump_Eng_Time ( Channel, Gmt, Status_Bits,
C			     Transmit_Frame, Numrec, Time_Flag, Rpt_Lun )
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
C	  Channel       C*2             FIRAS channel
C	  Gmt           C*14            Time tag of record
C	  Status_Bits   I*2             Micro header status bit word
C	  Transmit_Frame I*4            Micro counter for transmit frame
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
C	  Determine which bits are set in the status word and load the 
C	    bit values into an array for printing.
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

	Character*2 Channel          ! FIRAS channel
	Character*14 Gmt             ! Time tag of record
	Integer*2 Status_Bits        ! Micro header status bit word
	Integer*4 Transmit_Frame     ! Micro counter for transmit frame
	Integer*4 Numrec             ! Number of records with repeated time tag
	Character*2 Time_Flag        ! Flag indicating reason for dump
	Integer*4 Rpt_Lun            ! Report unit number

C	Local Declarations.

	Logical*1 First / .True. /
	Integer*4 Rstatus
	Character*19 Gmt_Char 
	Character*1 Cbits(16)
	Integer*4 Zero / 0 /
	Integer*2 Ix
        Character*1 Dash(130)/ 130 * '_'/

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr

	FTB_Dump_Eng_Time = %loc(FTB_Normal)

	If ( First ) Then
	  Write ( Unit=Rpt_lun, FMT=100, Iostat=Rstatus )
 100      Format (4x,'Transmit_Time',4x,
	1	  'Status_Bits VAX Ms->Ls',2x,
	2  	  'Transmit_MNF',1x,'Num.Rec.',1x,'Flag')
          Write( Unit=Rpt_lun, Fmt=150, Iostat=Rstatus) Dash
 150      format(1x,130a)
          First = .False.
	  If (Rstatus .NE. Zero ) Then
	    FTB_Dump_Eng_Time = %loc(FTB_Aberr)
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Do Ix = 0, 15
	  Cbits(16 - Ix) = '0'
	  If ( Btest( Status_Bits, Ix )) Cbits(16 - Ix) = '1'
	Enddo
	Gmt_Char= Gmt(1:2)//':'//Gmt(3:5)//':'//Gmt(6:7)//':'//
	1	  Gmt(8:9)//':'//Gmt(10:11)//'.'//Gmt(12:14)

	If ( FTB_Dump_Eng_Time .EQ. %loc(FTB_Normal)) Then
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Gmt_Char, Cbits,
	1     Transmit_Frame, Numrec, Time_Flag
 200      Format(1x,a,2x,4(4a,1x),3x,i10,2x,i10,2x,a)
	  If ( Rstatus .NE. Zero ) Then
!           Allow print even if overflow
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Return
	End
