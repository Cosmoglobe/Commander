	Integer*4 Function FUT_Field_Attributes ( lun, field_name,
     .						  length, offset )

C------------------------------------------------------------------------
C    PURPOSE: Retrieve a field's length and offset given the field name.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 2, 1987
C
C    INVOCATION: STATUS = FUT_Field_ATTRIBUTES ( FIELD_NAME, LENGTH,
C						 OFFSET )
C
C    INPUT PARAMETERS:
C	FIELD_NAME		CH*64		Fully qualified field name.
C
C    OUTPUT PARAMETERS: 
C	LENGTH			I*2		Field length in bytes (1 or 2).
C	OFFSET			I*2		Offset from beginning of record
C						(bytes; 0 based).
C
C    SUBROUTINES CALLED: 
C	FR_Init
C	FR_Open_Field
C	FR_Attribute
C	FR_Close
C	STR$UpCase
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	FUT_Error
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'($SSDef)'

	Integer		*4	i
	Integer		*4	j
	Logical		*4	status
	Integer		*4	lun
	Integer		*4	data_length
	Integer		*2	length
	Integer		*4	field_offset
	Integer		*2	offset
	Character	*64	field_name

	Logical		*4	FR_Open_Field
	Logical		*4	FR_Attribute
	Integer		*4	STR$UpCase

	External	FUT_Normal
	External	FUT_CFROpnFld
	External	FUT_CFRAttr

c
c Open the field.
c
	status = FR_Open_Field ( lun, field_name )

	If (.Not. status) Then
	   FUT_Field_Attributes = %Loc(FUT_CFROpnFld)
	   Call LIB$Signal(FUT_CFROpnFld,%Val(1),field_name)
	   Return
	End If

c
c Get the field length.
c
	status = FR_Attribute ( lun, field_name, 'DATA_LENGTH',
     .					%Descr(data_length) )
	length = data_length

	If (.Not. status) Then
	   FUT_Field_Attributes = %Loc(FUT_CFRAttr)
	   Call LIB$Signal(FUT_CFRAttr,%Val(1),field_name)
	   Return
	End If

c
c Get the field offset.
c
	status = FR_Attribute ( lun, field_name, 'FIELD_OFFSET',
     .					%Descr(field_offset) )
	offset = field_offset

	If (.Not. status) Then
	   FUT_Field_Attributes = %Loc(FUT_CFRAttr)
	   Call LIB$Signal(FUT_CFRAttr,%Val(1),field_name)
	   Return
	End If

	FUT_Field_Attributes = %Loc(FUT_Normal)

	Return
	End
