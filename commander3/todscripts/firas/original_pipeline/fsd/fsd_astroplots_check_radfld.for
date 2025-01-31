	Logical*4 Function FSD_Astroplots_Check_Radfld(n_lat, s_lat,
	1					       n_lon, s_lon,
	2					       n_npts, s_npts,
	3					       terr_lat, terr_lon)

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
C    INVOCATION: flag = FSD_Astroplots_Check_Radfld(n_lat, s_lat,
C						    n_lon, s_lon,
C						    n_npts, s_npts,
C						    terr_lat, terr_lon)
C
C    INPUT PARAMETERS:
C
C	N_LAT(100)	R*4		"North" latitude field table
C	S_LAT(100)	R*4		"South" latitude field table
C	N_LON(100)	R*4		"North" longitude field table
C	S_LON(100)	R*4		"South" longitude field table
C	N_NPTS		I*4		Number of points in the "North" table
C	S_NPTS		I*4		Number of points in the "South" table
C	TERR_LAT	R*4		Terrestial latitude to be tested
C	TERR_LON	R*4		Terrestial longitude to be tested
C
C    OUTPUT PARAMETERS: 
C
C	FLAG		L*4		Flag stating if the terrestial position
C					is within the radiation field.
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

	Real		*4	n_lat(100)
	Real		*4	s_lat(100)
	Real		*4	n_lon(100)
	Real		*4	s_lon(100)
	Integer		*4	n_npts
	Integer		*4	s_npts
	Real		*4	terr_lat
	Real		*4	terr_lon

	Real		*4	in_lat
	Real		*4	is_lat

	Integer		*4	j
	Logical		*4	field_flag
	Logical		*4	flag


C Perform the check iff the terrestial longitude to be checked is bounded
C by the minimum and maximum longitudes in the radiation field table.
C Remember that the field longitudes are sorted in increasing order.

	field_flag = .False.

	If ((terr_lon .Ge. n_lon(1) .And. terr_lon .Le. n_lon(n_npts)) .And.
	1   (terr_lon .Ge. s_lon(1) .And. terr_lon .Le. s_lon(s_npts))) Then

C First find bracketting longitudes from the "North" table and interpolate
C the latitude at terrestial longitude point.

	   j = 1
	   flag = .False.

	   Do While (.Not. flag)

	      If (terr_lon .Eq. n_lon(j)) Then
	         in_lat = n_lat(j)
	         flag = .True.
	      Else If (terr_lon .Gt. n_lon(j)) Then
	         j = j + 1
	      Else If (terr_lon .Lt. n_lon(j)) Then
	         in_lat = n_lat(j-1) + (n_lat(j)-n_lat(j-1)) *
	1			(terr_lon-n_lon(j-1)) / (n_lon(j)-n_lon(j-1))
	         flag = .True.
	      End If
	   End Do

C Now find bracketting longitudes from the "South" table and interpolate
C the latitude at terrestial longitude point.

	   j = 1
	   flag = .False.

	   Do While (.Not. flag)

	      If (terr_lon .Eq. s_lon(j)) Then
	         is_lat = s_lat(j)
	         flag = .True.
	      Else If (terr_lon .Gt. s_lon(j)) Then
	         j = j + 1
	      Else If (terr_lon .Lt. s_lon(j)) Then
	         is_lat = s_lat(j-1) + (s_lat(j)-s_lat(j-1)) *
	1			(terr_lon-s_lon(j-1)) / (s_lon(j)-s_lon(j-1))
	         flag = .True.
	      End If
	   End Do

C Check if the terrestial latitude is bracketted by the "north" and "south"
C latitudes interpolated from the radiation field tables.

	   If (terr_lat .Ge. is_lat .And. terr_lat .Le. in_lat) Then
	      field_flag = .True.
	   End If

	End If

	FSD_Astroplots_Check_Radfld = field_flag

	End
