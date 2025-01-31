	Integer*4 Function FEP_Get_GRT_Attr ( fr_lun, field, 
     .						length,
     .						offset, caloff, iside)

C------------------------------------------------------------------------
C    PURPOSE: Retrieves GRT Attributes: Cal resistor offset, GRT side 
C	      number.
C
C    INVOCATION: STATUS = FEP_Get_GRT_Attr ( FIELD, LENGTH,
C					     OFFSET, CALOFF, ISIDE)
C
C    INPUT PARAMETERS:
C	FIELD			CH*64		GRT field name.
C
C    OUTPUT PARAMETERS: 
C	LENGTH			I*2		Field length in bytes.
C	OFFSET			I*2		Field offset within record.
C	CALOFF			I*4		Cal resistor field offset.
C	ISIDE			I*4		AL,AH,BL,BH (1,2,3,4)
C
C  
C    AUTHOR: Harte Wang
C            STX
C            Oct. 27, 1987
C
C
C    SUBROUTINES CALLED: 
C	FUT_Field_Attributes
C	LIB$Get_LUN
C	LIB$Free_LUN
C	LIB$Signal
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	FUT_Error
C	$SSDef
C
C   PDL:
c      Call FUT_Field_Attributes to
c      Get the attributes of the field being converted.
c
C      IF Error then
C         return
C      Endif
C      IF Offset great or equal 36 and offset less or equal 66
C      Then
C        set calibration resistor field offset to 56
C        set side to  AL
C	Else If (offset .Ge. 68 .And. offset .Le. 98) Then
C        set calibration resistor field offset to 88
C        set side to  AH
C	Else If (offset .Ge. 100 .And. offset .Le. 130) Then
C        set calibration resistor field offset to 120
C        set side to  BL
C	Else If (offset .Ge. 132 .And. offset .Le. 162) Then
C        set calibration resistor field offset to 152
C        set side to  BH
C	End If
C	If status is error Then
C         Set the return code to error	
C	  set the message to the screen
C	Else
C         Set return code to good
C	Endif
C	Return
C	End
C  
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'($SSDef)'

	Character	*64	field
	Integer		*2	length
	Integer		*2	offset
	Integer		*4	fr_lun
	Integer		*4	caloff
	Integer		*4	iside
	Integer		*4	status

	Integer		*4	FUT_Field_Attributes

	External		FEP_Normal
	External		FUT_Normal

c
c Get the attributes of the field being converted.
c
	status = FUT_Field_Attributes ( fr_lun, field, length,
     .					offset )

	If (status .Ne. %Loc(FUT_Normal)) Then
	   FEP_Get_GRT_Attr = status
	   Return
	End If

	If (offset .Ge. 36 .And. offset .Le. 66) Then
	   caloff = 56
	   iside = 1
	Else If (offset .Ge. 68 .And. offset .Le. 98) Then
	   caloff = 88
	   iside = 2
	Else If (offset .Ge. 100 .And. offset .Le. 130) Then
	   caloff = 120
	   iside = 3
	Else If (offset .Ge. 132 .And. offset .Le. 162) Then
	   caloff = 152
	   iside = 4
	End If
        FEP_GET_GRT_ATTR=%loc(FUT_NORMAL)
	Return
	End
