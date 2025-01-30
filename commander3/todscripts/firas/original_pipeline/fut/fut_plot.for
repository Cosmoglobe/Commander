	Subroutine FUT_Plot ( plot_buff, npoints,
	1		      startx, spacex,
	2		      title, xlabel, ylabel,
	3		      zp, interactive, device,
	4		      plt_com, plt_com_file )

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
C	R. Kummerer, December 3, 1986. SPACEX was not being calculated
C	since it always assumed that data being plotted was noise.
C
C	R. Kummerer, January 16, 1987. Allow plotting a zero-point line.
C
C	R. Kummerer, February 12, 1987. Use full screen for plot. Use
C	TEMPLATE routines to print plot text.
C
C	R. Kummerer, May 8, 1987. Set axis types on linear, semi-log and
C	log-log plots.  Fix semi-log and log-log plot problem.
C
C	S. Alexander, March 5, 1990. Add capability to pass in a PLT
C	command file.
C	
C	S. Alexander, August 5, 1992. FUT_SETXAX related changes; SPR 9514.
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'

	Integer		*4	npoints
	Complex		*8	plot_buff(npoints)
	Real		*4	startx
	Real		*4	spacex
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

	Real		*4	y(2048, 4)
	Real		*4	sum
	Integer		*4	iery(10)
	Integer		*4	nveci
	Character	*150	cmd(55)
	Integer		*4	ncmd
	Integer		*4	ier

	Integer		*4	j

	External	FUT_Normal

c
c Load plot data into Y(*,*).
c
c	Y(1,*) contains the X-axis values computed from STARTX and SPACEX
c	       set up by FUT_SETXAX.
c	Y(2,*) contains the data values stored in the REAL portion of the
c	       "data compressed" plot buffer, which is always a non-zero
c	       vector.
c	Y(3,*) contains the data values stored in the IMAGINARY portion of
c	       the "data compressed" plot buffer.  If it is a non-zero vector,
c	       there is data to be plotted.  This can also be used as the
c	       optional "zero-point" vector if it contains no actual data to
c	       be plotted.  The "zero-point" vector is used to draw the x-axis
c	       on the plot.
c	Y(4,*) contains the optional "zero-point" vector if it contains no
c	       actual data to be plotted.  The "zero-point" vector is used
c	       to draw the x-axis on the plot.
c
	sum = 0.

	Do j=1,npoints

	   y(j,1) = (j-1.)*spacex + startx
	   y(j,2) = real ( plot_buff(j) )
	   y(j,3) = aimag ( plot_buff(j) )

	   sum = sum + y(j,3)

	   If (zp .Eq. fac_present) Then
	      y(j,4) = 0.
	   End If

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
	cmd(9) = 'PLOT OVERLAY'
	ncmd = 9

c
c Determine the number of vectors to be plotted.  Set the line styles.
c
	If (sum .Eq. 0) Then
	   If (zp .Eq. fac_present) Then
	      nveci = 3
	      ncmd = ncmd + 2
	      cmd(ncmd-1) = 'LS 1 ON 2'
	      cmd(ncmd)   = 'LS 2 ON 3'
	   Else
	      nveci = 2
	      ncmd = ncmd + 1
	      cmd(ncmd) = 'LS 1 ON 2'
	   End If
	Else
	   If (zp .Eq. fac_present) Then
	      nveci = 4
	      ncmd = ncmd + 3
	      cmd(ncmd-2) = 'LS 1 ON 2'
	      cmd(ncmd-1) = 'LS 4 ON 3'
	      cmd(ncmd)   = 'LS 2 ON 4'
	   Else
	      nveci = 3
	      ncmd = ncmd + 2
	      cmd(ncmd-1) = 'LS 1 ON 2'
	      cmd(ncmd)   = 'LS 4 ON 3'
	   End If
	End If

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
	Call PLT (y, iery, 2048, npoints, nveci, cmd, ncmd, ier)

	Call LIB$Erase_Page ( 1, 1 )

	Return
	End
