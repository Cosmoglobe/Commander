	Subroutine FSD_Astroplots_Display ( title, xlabel,
	1				    plt_device,
	2				    plt_com, plt_com_file )

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
C Changes:
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FSD_Astroplots)'

	Character	*100	title(3)
	Character	*(*)	xlabel
	Integer		*4	zp
	Character	*32	plt_device
	Integer		*4	plt_com
	Character	*64	plt_com_file

	Character	*100	title1
	Character	*100	title2
	Character	*100	title3

	Real		*4	y(6001, 4)
	Real		*4	sum
	Integer		*4	iery(10)
	Integer		*4	nveci
	Character	*150	cmd(55)
	Integer		*4	ncmd
	Integer		*4	ier

	Real		*4	degl_noise_max
	Real		*4	glitch_rate_max
	Real		*4	pk_ht_max
	Real		*4	plot_max

	Integer		*4	i
	Integer		*4	j
	Integer		*4	inum
	Integer		*4	mc

	External	FUT_Normal

c
c Load plot data into Y(*,*).
c
c	Y(1,*) contains timetags of the data points.
c	Y(2,*) contains the IFG peak height and galactic plane and Moon
c	       crossings.
c	Y(3,*) contains the glitch rate and VAB passages.
c	Y(4,*) contains the deglitcher noise and SAA passages.
c
	sum = 0.
	degl_noise_max = -9999.
	glitch_rate_max = -9999.
	pk_ht_max = -9999.

	Do j=1,num

	   y(j,1) = timeaxis(j)

	   If (pk_ht(j) .Eq. -9999.) Then
	      y(j,2) = fac_no_data
	   Else
	      y(j,2) = pk_ht(j)
	      If (pk_ht(j) .Gt. pk_ht_max)
	1		pk_ht_max = pk_ht(j)
	   End If

	   If (glitch_rate(j) .Eq. -9999.) Then
	      y(j,3) = fac_no_data
	   Else
	      y(j,3) = glitch_rate(j)
	      If (glitch_rate(j) .Gt. glitch_rate_max)
	1		glitch_rate_max = glitch_rate(j)
	   End If

	   If (degl_noise(j) .Eq. -9999.) Then
	      y(j,4) = fac_no_data
	   Else
	      y(j,4) = degl_noise(j)
	      If (degl_noise(j) .Gt. degl_noise_max)
	1		degl_noise_max = degl_noise(j)
	   End If

	End Do

	plot_max = Max(degl_noise_max, glitch_rate_max, pk_ht_max) * 1.05

c
c Build PLT commands to make the plot.
c
	If (interactive .Eq. fac_present) Then
	   Call LIB$Erase_Page ( 1, 1 )
	End If

	title1 = title(1)(1:Index(title(1),'\')-1)
	title2 = title(2)(1:Index(title(2),'\')-1)
	title3 = title(3)(1:Index(title(3),'\')-1)

	cmd(1) = 'D ' // plt_device
	cmd(2) = 'LA OT ' // title1
	cmd(3) = 'LA T '  // title2
	cmd(4) = 'LA F '  // title3
	cmd(5) = 'LA X '  // xlabel
	cmd(6) = 'LA Y ""'
	cmd(7) = 'LA OY ""'
	cmd(8) = 'LA OY2 "Peak Height"'
	cmd(9) = 'LA OY3 "Glitch Rate"'
	cmd(10) = 'LA OY4 "Degl Noise"'
	cmd(11) = 'LA Y2 "Cnts (GP,Mn)"'
	cmd(12) = 'LA Y3 "Gl/Sec (VAB)"'
	cmd(13) = 'LA Y4 "Cnts (SAA)"'
	cmd(14) = 'CS 1.35'
	cmd(15) = 'V .2 .1 .9 .8'
	cmd(16) = 'PLOT VERTICAL'
	ncmd = 16

	nveci = 4

        If (plt_com .Eq. fac_present) Then
           ncmd = ncmd + 1
           cmd(ncmd) = '@' // plt_com_file
        End If

c
c Initialize plot traces to the no data value for SAA and VAB crossings.
c
	inum = num

	If (nsaa .Gt. 0 .Or. nvab .Gt. 0 .Or. ngalactic .Gt. 0) Then

	   inum = inum + 1

	   Do i=1,4
	      y(inum,i) = fac_no_data
	   End Do

	   Do j=1,num
	      Do i=1,4
	         y(inum+j,i) = fac_no_data
	      End Do
	   End Do

	End If

c
c Mark galactic plane crossings.
c
	If (ngalactic .Gt. 0) Then

	   Do i=1,num
	      Do j=1,ngalactic-1,2
		 If (timeaxis(i) .Ge. galactic_cross(j) .And.
	1	     timeaxis(i) .Le. galactic_cross(j+1)) Then
	            y(inum+i,1) = timeaxis(i)
	            y(inum+i,2) = pk_ht_max * 1.05
		 End If
	      End Do
	   End Do

	End If

c
c Mark VAB crossings.
c
	If (nvab .Gt. 0) Then

	   Do i=1,num
	      Do j=1,nvab-1,2
		 If (timeaxis(i) .Ge. vab(j) .And.
	1	     timeaxis(i) .Le. vab(j+1)) Then
	            y(inum+i,1) = timeaxis(i)
	            y(inum+i,3) = glitch_rate_max * 1.05
		 End If
	      End Do
	   End Do

	End If

c
c Mark SAA crossings.
c
	If (nsaa .Gt. 0) Then

	   Do i=1,num
	      Do j=1,nsaa-1,2
		 If (timeaxis(i) .Ge. saa(j) .And.
	1	     timeaxis(i) .Le. saa(j+1)) Then
	            y(inum+i,1) = timeaxis(i)
	            y(inum+i,4) = degl_noise_max * 1.05
		 End If
	      End Do
	   End Do

	End If

	If (nsaa .Gt. 0 .Or. nvab .Gt. 0 .Or. ngalactic .Gt. 0) Then
	   inum = inum + num
	End If

c
c Mark Moon crossings.
c
	If (nmoon .Gt. 0) Then

	   mc = 0

	   Do i=1,num-1
	      Do j=1,nmoon
c
c		 Moon angle 5d.
c
		 If (amoon(j) .Eq. 5) Then
		    If (timeaxis(i)   .Le. tmoon(j) .And.
	1	        timeaxis(i+1) .Ge. tmoon(j)) Then
		       ncmd = ncmd + 1
		       mc = mc + 1
		       write (cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), 3a1)')
	1		'LA ', mc, ' P ', tmoon(j), plot_max, '"', 'O', '"'
		    End If
		 End If
c
c		 Moon angle 15d.
c
		 If (amoon(j) .Eq. 15) Then
		    If (timeaxis(i)   .Le. tmoon(j) .And.
	1	        timeaxis(i+1) .Ge. tmoon(j)) Then
		       ncmd = ncmd + 1
		       mc = mc + 1
		       write (cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), 3a1)')
	1		'LA ', mc, ' P ', tmoon(j), plot_max, '"', 'C', '"'
		    End If
		 End If
c
c		 Moon angle 30d.
c
		 If (amoon(j) .Eq. 30) Then
		    If (timeaxis(i)   .Le. tmoon(j) .And.
	1	        timeaxis(i+1) .Ge. tmoon(j)) Then
		       ncmd = ncmd + 1
		       mc = mc + 1
		       write (cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), 3a1)')
	1		'LA ', mc, ' P ', tmoon(j), plot_max, '"', '(', '"'
		    End If
		 End If

	      End Do
	   End Do

	End If

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
	Call PLT (y, iery, 6001, inum, nveci, cmd, ncmd, ier)

	If (interactive .Eq. fac_present) Then
	   Call LIB$Erase_Page ( 1, 1 )
	End If

	Return
	End
