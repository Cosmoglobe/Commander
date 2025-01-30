	Subroutine FUT_Scatter_Plot ( x, y, npoints,
	2			      title, xlabel, ylabel,
	3			      zp, interactive, device,
	4			      plt_com, plt_com_file )

C------------------------------------------------------------------------
C
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Incorporated
C            September 27, 1989
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
C    Changes:
C
C	S. Alexander, March 5, 1990.  Add capability to pass in a PLT
C	command file.
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'

	Integer		*4	npoints
	Real		*4	x(npoints)
	Real		*4	y(npoints)
	Character	*100	title(3)
	Character	*(*)	xlabel
	Character	*(*)	ylabel
	Integer		*4	zp
	Integer		*4	interactive
	Character	*32	device
	Integer		*4	plt_com
	Character	*64	plt_com_file

	Character	*100	title1
	Character	*100	title2
	Character	*100	title3

	Real		*4	p(2048, 2)
	Real		*4	sum
	Integer		*4	iery(10)
	Integer		*4	nveci / 2 /
	Character	*150	cmd(55)
	Integer		*4	ncmd
	Integer		*4	ier

	Integer		*4	j

	External	FUT_Normal

c
c Load plot data into P(*,*).
c
c	P(1,*) contains the X-axis values from array X
c	P(2,*) contains the Y-axis values from array Y
c
	sum = 0.

	Do j=1,npoints

	   p(j,1) = x(j)
	   p(j,2) = y(j)

	End Do

c
c Make the plot...
c
	Call LIB$Erase_Page ( 1, 1 )

	title1 = title(1)(1:Index(title(1),'\')-1)
	title2 = title(2)(1:Index(title(2),'\')-1)
	title3 = title(3)(1:Index(title(3),'\')-1)

	cmd(1) = 'D ' // device
	cmd(2) = 'LA OT ' // title1
	cmd(3) = 'LA T '  // title2
	cmd(4) = 'LA F '  // title3
	cmd(5) = 'LA X '  // xlabel
	cmd(6) = 'LA OY ' // ylabel
	cmd(7) = 'CS 1.35'
	cmd(8) = 'V .2 .1 .9 .8'
	cmd(9) = 'LINE OFF'
	cmd(10) = 'MARKER 1 ON 2'
	cmd(11) = 'COLOR 1 ON 2'
	ncmd = 11

        if (plt_com .eq. fac_present) then
            ncmd = ncmd + 1
            cmd(ncmd) = '@' // plt_com_file
        endif
c
c Tell PLT to plot if we are not in interactive mode.
c
	If (interactive .Eq. fac_not_present) Then
	   ncmd = ncmd + 2
	   cmd(ncmd-1) = 'P'
	   cmd(ncmd) = 'Q'
	End If

	Do j=1,nveci
	   iery(j) = 0.
	End Do

c
c Display the plot until the user finishes looking at it.
c
	Call PLT (p, iery, 2048, npoints, nveci, cmd, ncmd, ier)

	Call LIB$Erase_Page ( 1, 1 )

	Return
	End
