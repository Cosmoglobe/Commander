	Subroutine FUT_Plot_Title ( data_label, ngroup, nsweeps,
	1			    mtm_speed, mtm_length,
	2			    chan, micro, xcal_pos, fake_it, gain,
	3			    plot_label, title )

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
C            July 26, 1989
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
C
C Changes:
C
C	Version 4.4.1, R. Kummerer, STX, SPR 4570 Sept 15, 1989.
C		Conserve space available for the plot label by removing
C		text string 'LABEL:' from the plot label.
C
C----------------------------------------------------------------------

	Implicit 	None

	Include 	'($JPIDef)'
	Include 	'($SSDef)'
	Include 	'(FUT_Params)'

	Character 	*60	data_label	!label from data record
	Integer 	*4	ngroup		!number of adds per group
	Integer 	*4	nsweeps		!number of MTM sweeps
	Integer 	*4	mtm_speed      	!scan speed 0=slow 1=fast
	Integer 	*4	mtm_length	!scan length 0=short 1=long
	Integer 	*4	chan      	!channel id
	Integer 	*4	micro		!microprocessor mode
	Integer 	*4	xcal_pos	!XCal position
	Integer 	*4	fake_it		!fake-it mode
	Integer 	*4	gain		!gain factor
	Character 	*(*)	plot_label	!title of plot
	Character	*100	title(3)	!constructed plot title

	Character 	*9	day      	!date
	Character 	*8	whom      	!plot by
	Integer		*4	len
	Character	*2	adds
	Character	*2	sweeps
	Character	*1	umode
	Character	*4	gfact
	Integer		*4	status

	Integer		*4	LIB$GetJPI
	Integer		*4	STR$Trim

C
C Place the plot label in the first element of the plot title array.
C
	status = STR$Trim ( plot_label, plot_label, len )

	title(1) = plot_label(1:len) // '\'

C
C Place the instrument mode information in the second element of the
C plot title array.  First display the channel.
C
	If (chan .Ge. 1 .And. chan .Le. 4) Then
	   title(2) = 'Ch ' // fac_channel_ids(chan) // '\'
	Else
	   title(2) = 'Ch ?\'
	End If

C
C Display the adds per group and the number of MTM sweeps.
C
	If (ngroup .Ge. 0 .And. ngroup .Le. 12) Then
	   Write (adds, '(i2)') ngroup
	Else
	   adds = ' ?'
	End If

	If (nsweeps .Gt. 0 .And. nsweeps .Le. 99) Then
	   Write (sweeps, '(i2)') nsweeps
	Else
	   sweeps = '? '
	End If

	title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  Ad/Sw ' // adds // '/' // sweeps // '\'

C
C Display the MTM attributes.
C
	If (mtm_speed .Eq. 0) Then

	   If (mtm_length .Eq. 0) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Sh Sl\'
	   Else If (mtm_length .Eq. 1) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Lo Sl\'
	   Else 
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM ? Sl\'
	   End If

	Else If (mtm_speed .Eq. 1) Then

	   If (mtm_length .Eq. 0) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Sh Fa\'
	   Else If (mtm_length .Eq. 1) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Lo Fa\'
	   Else 
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM ? Fa\'
	   End If

	Else 

	   If (mtm_length .Eq. 0) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Sh ?\'
	   Else If (mtm_length .Eq. 1) Then
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM Lo ?\'
	   Else 
	      title(2) = title(2)(1:Index(title(2),'\')-1) //
	1			'  MTM ? ?\'
	   End If

	End If

C
C Display the microprocessor mode.
C
	If (micro .Ge. 0 .And. micro .Le. 4) Then
	   Write (umode, '(i1)') micro
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  uP ' // umode // '\'
	Else
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  uP ?\'
	End If

C
C Display the XCal position.
C
	If (xcal_pos .Eq. fac_xcalin) Then
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  XCal IN\'
	Else If (xcal_pos .Eq. fac_xcalout) Then
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  XCal OUT\'
	Else If (xcal_pos .Eq. fac_xcaltrans) Then
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  XCal TRN\'
	Else If (xcal_pos .Eq. fac_xcalout) Then
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  XCal ERR\'
	End If

C
C Display the gain factor.
C
	If (gain .Eq.    1 .Or.
	1   gain .Eq.    3 .Or.
	2   gain .Eq.   10 .Or.
	3   gain .Eq.   30 .Or.
	4   gain .Eq.  100 .Or.
	5   gain .Eq.  300 .Or.
	6   gain .Eq. 1000 .Or.
	7   gain .Eq. 3000) Then
	   Write (gfact, '(i4)') gain
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  Gn ' // gfact // '\'
	Else
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  Gn ?\'
	End If

C
C Display the fake-it mode.
C
	If (fake_it .Eq. 0) Then
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  Fk OFF\'
	Else
	   title(2) = title(2)(1:Index(title(2),'\')-1) //
	1		'  Fk ON\'
	End If

C
C Place the data label in the third element of the plot title array.
C
	status = STR$Trim ( data_label, data_label, len )

	title(3) = data_label(1:len) // '\'

	Return
	End
