	Subroutine FSD_Astroplots_Label ( data_label, chan, att_soln,
	1				  plot_label, title )

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            April 26, 1990
C
C    INVOCATION: 
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: 
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES: 
C
C    PROCESSING METHOD:
C  
C----------------------------------------------------------------------
C Changes:
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C----------------------------------------------------------------------

	Implicit 	None

	Include 	'(FUT_Params)'

	Character 	*60	data_label	!label from data record
	Integer 	*4	chan      	!channel id
	Character 	*20	att_soln	!attitude solution
	Character 	*(*)	plot_label	!title of plot
	Character	*100	title(3)	!constructed plot title

	Integer		*4	len
	Integer		*4	status

	Integer		*4	STR$Trim

C
C Place the plot label in the first element of the plot title array.
C
	status = STR$Trim ( plot_label, plot_label, len )

	If (chan .Ge. 1 .And. chan .Le. 4) Then
	   title(1) = plot_label(1:len) // '  Ch ' // fac_channel_ids(chan) //
	1		'  Attitude=' // att_soln
	Else
	   title(1) = plot_label(1:len) // '  Ch ?' //
	1		'  Attitude=' // att_soln
	End If

C
C Place the data label in the second element of the plot title array.
C
	status = STR$Trim ( data_label, data_label, len )

	title(2) = data_label(1:len) // '\'

C
C Display the plot legend.
C
	title(3) = 'Gal Pln=__, Moon: 5d=O, 15d=C, 30d=(, VAB=__, SAA=__'

	Return
	End
