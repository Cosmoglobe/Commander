	Integer*4 Function FUT_Smooth_Data ( data, num, smooth_data,
     .					     sfcn, window )

C------------------------------------------------------------------------
C    PURPOSE: Performs simple smoothing functions.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            July 15, 1988
C
C    INVOCATION: status = FUT_SMOOTH_DATA ( data, num, smooth_data,
C					    sfcn, window )
C
C    INPUT PARAMETERS:
C	DATA(NUM)		R*4	Input data to be smoothed.
C	NUM			I*4	Number of data points.
C	SFCN			I*4	Smooth function.
C					   FAC_BOXCAR (=1) = boxcar
C	WINDOW			I*4	Boxcar window size.
C
C    OUTPUT PARAMETERS: 
C	STATUS			I*4	Success flag.
C	SMOOTH_DATA		R*4	Smoothed data.
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	FUT_PARAMS.TXT
C  
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Params)'

	Integer		*4	num
	Real		*4	data(num)
	Real		*4	smooth_data(num)
	Integer		*4	sfcn
	Integer		*4	window

	Integer		*4	np
	Integer		*4	ix
	Integer		*4	jx
	Integer		*4	kx
	Integer		*4	weight
	Real		*4	acc_data

	External	FUT_Normal
	External	FUT_UnkSmth

C
C Perform boxcar smoothing.
C
	weight = window

	np = weight - 1
	np = np / 2

	If (sfcn .Eq. fac_boxcar) Then
C
C Compute the averaged points of the smoothed IFGs using the 
C boxcar function and the user selected number of points.
C
	      Do jx = np + 1, num - np

	         acc_data = 0.0

	         Do kx = jx - np, jx + np
	            acc_data = acc_data + data(kx)
	         End Do

	         smooth_data(jx) = acc_data / weight

	      End Do	  
C
C Estimate the smoothed points at the end points of the IFGs.
C
	      Do jx = 1, np

	         weight = np + jx
	         acc_data = 0.0

	         Do kx = 1, jx + np
	            acc_data = acc_data + data(kx)
	         End Do

	         smooth_data(jx) = acc_data / weight

              End Do

	      Do jx = num - np + 1, num

	         weight = np + (num + 1) - jx
	         acc_data = 0.0

	         Do kx = jx - np, num
	            acc_data = acc_data + data(kx)
	         End Do

	         smooth_data(jx) = acc_data / weight

              End Do

	Else

	   Call LIB$Signal ( FUT_UnkSmth )
	   FUT_Smooth_Data = %Loc(FUT_UnkSmth)
	   Return

	End If

	FUT_Smooth_Data = %Loc(FUT_Normal)

	Return
	End
