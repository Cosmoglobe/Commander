	Integer*4 Function FUT_LonLat(gmt, gmtstart, gmtstop,
	2                             terr_lon, terr_lat)

C-----------------------------------------------------------------------------
C   Purpose:  Callable function to call the orbit routine, UOX_Get_Orbit, and 
C      return Terrestrial position.  Call it once for each time (gmt) an orbit
C      record is wanted.  On the first call to this routine, give it start and
C      stop times (gmtstart, gmtstop) to open the orbit archive.
C
C
C   Argument list:
C	I/O	Name		Type		Description
C	---	---------	-----		---------------------------
C	 I	gmt		Ch*14		! Time for current call
C	 I	gmtstart	Ch*14   	! Start time for orbit data
C	 I	gmtstop		Ch*14   	! Stop time for orbit data
C	 O	terr_lon	 R* 4		! Terrestrial longitude
C	 O	terr_lat	 R* 4		! Terrestrial latitude
C
C
C   Function/Subroutine Calls:
C
C	FUT_Open_Orbit_Arcv
C	UOX_Get_Orbit
C	AUT_ADT2DFloat
C	GST			( filename UCB_GST.FOR )
C
C
C   Written by:
C
C            Fred Shuman,  ST Systems Corp.,  1990 Jan 5
C-----------------------------------------------------------------------------

	Implicit None

C    Passed arguments

	Character*14	gmt		! Time for current call
	Character*14    gmtstart	! Start time for orbit data
	Character*14    gmtstop		! Stop time for orbit data
	Real     * 4	terr_lon
	Real     * 4	terr_lat

C    Include files

	Include '(UOX_Time_Str)'
	Record /UOX_Time_Str/ Time_Rec
	Include '(UOX_Ephem_Str)'
	Record /Ephem/ Ephem_Rec

C    Functions

	Integer  * 4	FUT_Open_Orbit_Arcv
	Integer  * 4	UOX_Get_Orbit
	Real     * 8	AUT_ADT2DFloat
	Real     * 4	GST			! filename UCB_GST.FOR

C    Messages

	External	FUT_NORMAL
	External	FUT_ATTORBOPNERR
	External	UOX_SUCCESS

C    Local variables

	Real     * 4	Pi
	Parameter	( Pi = 3.1415926536 )
	Integer  * 4	i		! Loop counter
	Integer  * 4	status		! Dummy status variable
	Integer  * 4	adtstart(2)	! Binary start time
	Integer  * 4	adtstop(2)	! Binary stop time
	Integer  * 4	adt(2)		! Binary time
	Integer  * 4	orb_lun		! Logical unit for orbit archive
	Real     * 8	t		! Time (10^-7 sec) since reference of
					!    JD 2400000.5 = 1858 Nov 17 0h UT
	Real     * 4	sidtim		! Sidereal "time", in radians
C
C    Note: Vector u is represented in Cartesian coordinates.
C        Reference Frame origin, Frame orientation, Units are given.
C        "Terr", "Terrestrial" denote Earth-co-rotating frame.
C
	Real     * 4	COBE_Geo_Pos(3)	! COBE position--time t
	Real     * 4	u(3)		! COBE position unit-vector--
					!   Geocenter, (terr) equatl, Dim'less
	Real     * 4	temp		! Temporary storage variable
	Logical  * 4	first /.True./	! First time flag--open orb arcv 1st 
					!   time only

C   Open the orbit archive and call UOX_Get_Orbit to access the orbit records.

	Save first

	Call CT_GMT_to_Binary(gmtstart, adtstart)
	Call CT_GMT_to_Binary(gmtstop, adtstop)
	time_rec.number_of_points = 1
	If (first) Then
	   first = .False.
	   status = FUT_Open_Orbit_Arcv ( gmtstart, gmtstop, orb_lun )
	   If ( status .Ne. %Loc(FUT_NORMAL) ) Then
	      FUT_LonLat = %Loc(FUT_ATTORBOPNERR)
	      Return
	   End If
	End If

	Call CT_GMT_to_Binary(gmt, adt)
	Do i = 1, 2
	   time_rec.initial_time(i) = adt(i)
	End Do

C  Get the COBE position.  Because it is returned in m, we multiply by 1.E-3 to
C   convert to km.

	status = UOX_Get_Orbit ( time_rec, orb_lun, ephem_rec )
	If ( status .Eq. %Loc(UOX_Success) ) Then
	   COBE_Geo_Pos(1) = ephem_rec.x_pos * 1.E-3
	   COBE_Geo_Pos(2) = ephem_rec.y_pos * 1.E-3
	   COBE_Geo_Pos(3) = ephem_rec.z_pos * 1.E-3
	Else
	   Call LIB$Signal(status)
	   FUT_LonLat = status
	   Return
	End If

C  Convert the binary time to time since 5/24/68 in seconds.
C   E. L. "Ned" Wright routines in UCB use JD 2440000.5 (= 1968 May 24, 0h UT)
C   as reference, with sec as unit; AUT_ADT2DFloat uses the VAX ADT reference,
C   viz., JD 2400000.5 (=1858 Nov 17 00h UT), with 10^-7 sec as unit.

	t = 1.D-7 * AUT_ADT2DFloat( adt ) - 40000.D0 * 86400.D0

C  Convert COBE Celestial position to Terrestrial (frame co-rotates with Earth)
C   and normalize to form a unit vector.  Conversion calculation swiped from
C   E. L. Wright's routine, TERRESTIAL (sic), to be found in UCB.

	sidtim = GST(t)
	u(1) =  cos(sidtim)*COBE_Geo_Pos(1) + sin(sidtim)*COBE_Geo_Pos(2)
	u(2) = -sin(sidtim)*COBE_Geo_Pos(1) + cos(sidtim)*COBE_Geo_Pos(2)
	u(3) =              COBE_Geo_Pos(3)
	Call Norm ( u )

C  Compute TERR_LON and TERR_LAT from the Terrestrial position unit vector of
C   COBE by converting to spherical coordinates.  Use more robust atan2 for
C   computing latitude.

	terr_lon = atan2d(u(2), u(1))
	temp = sqrt(u(1)**2 + u(2)**2)
	terr_lat = atan2d(u(3), temp)

C  Done.  Exit.

	FUT_LonLat = %Loc(FUT_Normal)

	Return
	End
