C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Get_MTM_Mode ( Eng_Rec, Scan_Mode, Dom_Mode )

C-------------------------------------------------------------------------------
C
C	Purpose: To determine the scan mode for the four channels and the
C	         predominating mode among the four.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Get_MTM_Mode ( Eng_Rec, Scan_Mode, Dom_Mode )
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
C	  Eng_Rec       FDQ_ENG         FIRAS engineering record
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Scan_Mode(4)  I*2             Four channel scan modes
C	  Dom_Mode      I*2             Predominant scan mode
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
C	  Set the status to success.
C         Get the MTM scan length and speed from the record for the four 
C	    channels.
C         Determine the scan mode (0 to 3) for each channel.
C	  Determine which scan mode predominates.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

C	Passed Parameters.

	Dictionary 'FDQ_ENG'
	Record /FDQ_ENG/ Eng_Rec
	Integer*2 Scan_Mode(4)       ! Four channel scan modes
	Integer*2 Dom_Mode           ! Predominant scan mode
	
C	Local Declarations.

	Integer*2 Length(4)	     ! MTM scan length for 4 channels
	Integer*2 Speed(4)	     ! MTM scan speed for 4 channels
	Integer*2 Mode_Cnt(4)        ! Count of channels in each mode for record
	Integer*2 Max_Mode           ! Maximum mode number on record
	Logical*1 Val_Chan(4)	     ! Channel science record associated to Eng
	Logical*1 Nosci              ! No channel science data associated
	Integer*2 Ix

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NoSci

	FTB_Get_MTM_Mode = %loc(FTB_Normal)

C	Get the MTM scan length and speed for 4 channels.

	Do Ix = 1, 4
	  Length(Ix) = Eng_Rec.Chan(Ix).Xmit_MTM_Len
	  Speed(Ix) = Eng_Rec.Chan(Ix).Xmit_MTM_Speed
	Enddo

C	Get the scan mode for each channel.

	Do Ix = 1, 4
	  If ((Length(Ix) .EQ. 0) .AND. (Speed(Ix) .EQ. 0)) Then
	    Scan_Mode(Ix) = 0
	  ElseIf ((Length(Ix) .EQ. 0) .AND. (Speed(Ix) .EQ. 1)) Then
	    Scan_Mode(Ix) = 1
	  ElseIf ((Length(Ix) .EQ. 1) .AND. (Speed(Ix) .EQ. 0)) Then
	    Scan_Mode(Ix) = 2
	  ElseIf ((Length(Ix) .EQ. 1) .AND. (Speed(Ix) .EQ. 1)) Then
	    Scan_Mode(Ix) = 3
	  Endif
	Enddo

C	Find the predominant mode among the 4 channels.

	Do Ix = 1, 4
	  Mode_Cnt(Ix) = 0
	  Val_Chan(Ix) = .False.
	Enddo

	Do Ix = 1, 4
	  If ((Eng_Rec.En_Head.Sci_Time(Ix).Bin_Time(1) .NE. 0) .OR.
	1     (Eng_Rec.En_Head.Sci_Time(Ix).Bin_Time(2) .NE. 0)) Then
	    Mode_Cnt(Scan_Mode(Ix)+1) = Mode_Cnt(Scan_Mode(Ix)+1) + 1
	    Val_Chan(Ix) = .True.
	  Endif
	Enddo	
	Max_Mode = Mode_Cnt(1)
	Dom_Mode = 0
	Do Ix = 2, 4
	  If (Mode_Cnt(Ix) .GT. Max_Mode) Then
	    Dom_Mode = Ix - 1
	    Max_Mode = Mode_Cnt(Ix)
	  Endif
	Enddo
	Nosci = .True.
	Do Ix = 1, 4
	  If (Val_Chan(Ix)) Nosci = .False.
	Enddo  	  
	If (Nosci) Then	   
	   Call Lib$Signal(FTB_NoSci)
	   FTB_Get_MTM_Mode = %loc(FTB_Aberr)
	Endif

	Return
	End
