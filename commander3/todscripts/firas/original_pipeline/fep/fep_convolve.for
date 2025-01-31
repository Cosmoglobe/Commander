	Subroutine FEP_Convolve (y, dwellbuff, npts, select, fieldnum,
	2                        crossgaps, nwindow)

C------------------------------------------------------------------------
C    PURPOSE: Smoothes engineering plots data.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: R. Isaacman
C            ARC
C
C    INVOCATION: CALL FEP_CONVOLVE ( Y, DWELLBUFF, NPTS, SELECT, FIELDNUM,
C	                             CROSSGAPS, NWINDOW )
C
C    INPUT PARAMETERS:
C	Y(*)		R*4	Data to be smoothed.
C	DWELLBUFF(*)	B*1	Dwell flag for each data point.
C	NPTS		I*4	Number of data points.
C	SELECT		I*4	Flag controlling smooth selection.
C	FIELDNUM	I*4	Sequential number of field being smoothed.
C	CROSSGAPS      Ch*1	Whether smoothing should cross data gaps.
C
C    OUTPUT PARAMETERS: 
C	Y(*)		R*4	Smoothed data.
C	NWINDOW		I*4	Width of smooth window.
C
C    SUBROUTINES CALLED:
C	STR$Trim
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES:
C	FEP_Invoc
C	FUT_Params
C  
C----------------------------------------------------------------------
C
C Changes:
C
C	/SMOOTH window select for batch. R. Kummerer, July 16, 1987.
C
C	Adapted for PLT plot package.  F. Shuman,  1988 Apr 28.
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FEP_Invoc)'
	Include		'(FUT_Params)'

	Real		*4	y(1)
	Byte			dwellbuff(1)
	Integer		*4	npts
	Integer		*4	select
	Integer		*4	fieldnum
	Logical		*1	crossgaps
	Integer		*4	nwindow

	Real		*4	buff(fac_max_hskp_recs)
	Real		*4	sum
	Integer		*4	index
	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	nsum
	Integer		*4	ihalf
	Integer		*4	irem
	Logical		*1	ans_ok
	Character	*1	ans
	Integer 	*4	status

	Integer		*4	STR$Trim


	If (select .Eq. 1) Then

	   ans_ok = .False.

	   Do While (.Not. ans_ok)

	      Type 20, fieldnum
20	      Format (1x, 'Enter window width for field', i3,
	2                 ' (odd integer only): ',$)
	      Accept *, nwindow

	      irem = Mod(nwindow,2)

	      If (irem .Eq. 1 .And. nwindow .Ge. 0) Then
	         ans_ok = .True.
	      Else
	         Type *, ' Please try again.'
	      End If

	   End Do

	Else

	   nwindow = fcc_window

	End If

	nwindow = IAbs (nwindow)
	If (Mod(nwindow,2) .Eq. 0) Then
	   nwindow = nwindow + 1
	End If
	ihalf = nwindow/2

	Do j=1,npts
	   If (y(j) .Eq. fac_no_data .Or. dwellbuff(j) .Eq. 1) Then
	      buff(j) = fac_no_data
	   Else
	      sum = y(j)
	      nsum = 1

	      If (crossgaps) Then
	         index = j - 1
	         Do k=1,ihalf
	            Do While (index .Ge. 1   .And.
	2             (y(index) .Eq. fac_no_data .Or. dwellbuff(index) .Eq. 1))
	               index = index - 1
	            End Do
	            If (index .Ge. 1) Then
	               sum = sum + y(index)
	               nsum = nsum + 1
	            End If
	            index = index - 1
	         End Do

	         index = j + 1
	         Do k=1,ihalf
	            Do While (index .Le. npts   .And.
	2             (y(index) .Eq. fac_no_data .Or. dwellbuff(index) .Eq. 1))
	               index = index + 1
	            End Do
	            If (index .Le. npts) Then
	               sum = sum + y(index)
	               nsum = nsum + 1
	            End If
	            index = index + 1
	         End Do
	      Else
	         index = j - 1
	         Do While (index .Ge. j-ihalf .And. index .Ge. 1 .And.
	2              y(index) .Ne. fac_no_data .And. dwellbuff(index) .Ne. 1)
	            sum = sum + y(index)
	            nsum = nsum + 1
	            index = index - 1
	         End Do

	         index = j + 1
	         Do While (index .Le. j+ihalf .And. index .Le. npts .And.
	2              y(index) .Ne. fac_no_data .And. dwellbuff(index) .Ne. 1)
	            sum = sum + y(index)
	            nsum = nsum + 1
	            index = index + 1
	         End Do
	      End If

	      buff(j) = sum/nsum
	   End If
	End Do

	Do j=1,npts
	   y(j) = buff(j)
	End Do

	Return
	End
