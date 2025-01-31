	Subroutine FSD_Astroplots_IntAng(time,g_lat,m_ang,t_lat,t_lon)

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
C            October 25, 1989
C
C    INVOCATION: Call FSD_Astroplots_IntAng(time,g_lat,m_ang,t_lat,t_lon)
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: None
C
C    PROCESSING METHOD:
C  
C----------------------------------------------------------------------

	Implicit None

	Include		'(FSD_Astroplots)'

	Real		*4	time
	Real		*4	g_lat
	Real		*4	m_ang
	Real		*4	t_lat
	Real		*4	t_lon

	Integer		*4	i
	Integer		*4	il
	Integer		*4	iu
	Integer		*4	lb
	Integer		*4	ub

	Logical		*1	found

C The time in question should be bracketted by times in the TIMEAXIS array.
C If not however, return appropriate angles for the situation anyway.

	If (time .Le. timeaxis(1)) Then

	   g_lat = gal_lat(1)
	   m_ang = moon_angle(1)
	   t_lat = terr_lat(1)
	   t_lon = terr_lon(1)

	Else If (time .Ge. timeaxis(num)) Then

	   g_lat = gal_lat(num)
	   m_ang = moon_angle(num)
	   t_lat = terr_lat(num)
	   t_lon = terr_lon(num)

	Else

C Perform a binary search through the TIMEAXIS time tags to find the two
C times that bracket the time of interest.

	   found = .False.

	   il = 1
	   iu = num
	   i = (iu-il)/2

 	   Do While (.Not. found)

	      If (time .Eq. timeaxis(i) .Or. iu-il .Eq. 1) Then
	         found = .True.
	         lb = il
	         ub = iu
	      Else If (time .Gt. timeaxis(i) .And. iu-il .Ne. 1) Then
		 il = i
	         i = (iu-il)/2 + il
	      Else If (time .Lt. timeaxis(i) .And. iu-il .Ne. 1) Then
		 iu = i
	         i = (iu-il)/2 + il
	      End If

	   End Do

	   g_lat = gal_lat(lb) + (gal_lat(ub)-gal_lat(lb)) *
	1		(time-timeaxis(lb)) / (timeaxis(ub)-timeaxis(lb))
	   m_ang = moon_angle(lb) + (moon_angle(ub)-moon_angle(lb)) *
	1		(time-timeaxis(lb)) / (timeaxis(ub)-timeaxis(lb))
	   t_lat = terr_lat(lb) + (terr_lat(ub)-terr_lat(lb)) *
	1		(time-timeaxis(lb)) / (timeaxis(ub)-timeaxis(lb))
	   t_lon = terr_lon(lb) + (terr_lon(ub)-terr_lon(lb)) *
	1		(time-timeaxis(lb)) / (timeaxis(ub)-timeaxis(lb))

	End If

	Return
	End
