	Integer*4 Function FUT_Attitude ( sci_rec, attitude_type,
     &                                    gmt_start, gmt_stop )
c-----------------------------------------------------------------------------
c
c      Subroutine to calculate the FIRAS attitude and attitude dependent
c      quantities.  These include the Sun-Moon and COBE-Moon distance, the
c      Earth limb, the Moon aspect and angle, the Sun angle, the altitude
c      and the barycentric velocity vector.
c
c      Started by:			Re-written and expanded by:
c
c            J. W. Durachta		Fred Shuman and Shirley Read
c        ST Systems Corporation		   ST Systems Corporation
c             April, 1987		     1988 Jan, May, Oct
c
c
c   PDL for FUT_Attitude
c
c	FUT_ATTITUDE calculates the FIRAS attitude and attitude-dependent
c   quantities.  These include the the FIRAS pixel number; FIRAS celestial
c   equatorial pointing put into Cartesian components, as well as right
c   acsension and declination; the terrestrial pixel number, and terrestrial
c   latitude and longitude, describing the geocentric COBE position; the Earth
c   limb and Earth limb azimuth angles, the Sun angle, the Moon and Moon
c   azimuth angles; the Moon phase; the Sun-Moon distance; the COBE-Moon
c   distance; the altitude of COBE; the component of barycentric velocity
c   along the COBE line-of-sight; the McIlwain "L" parameter, characterizing
c   the effective shielding of the geomagnetic field against incoming charged
c   particles at COBE's position; the galactic longitude and latitude; the
c   ecliptic longitude and latitude; the orbital phase (the angle COBE has
c   swept through in orbit since its last ascending node crossing); the
c   component of geocentric velocity along the COBE line-of-sight; the scan
c   angle (the 'azimuthal' component of FIRAS look direction about the
c   COBE-Sun line, measured from ascending passage through the ecliptic plane);
c   and the spacecraft rotation angle.
c	The program is organized into three sections.  In the first section, the
c   necessary COBE orbit and attitude information is gathered, the source of
c   which depends on whether simulated or real attitude was specified.  In the
c   second section, the required Solar system ephemerides are obtained from
c   XFM_Locate_Object, and all the quantities that will be needed in the final
c   section are prepared.  In the third section, the actual quantities to be
c   written in the Attitude record are computed and written.  You have probably
c   detected some ambiguity in whether a given intermediate quantity is
c   calculated in section 2 or 3.  This is largely accurate, but generally,
c   anything that is needed in computing more than one attitude quantity (AQ) is
c   found in section 2.  This is done to allow calculation of the individual AQs
c   to be re-ordered at will within section 3, without having to worry that some
c   value needed to find an AQ is not ready until a later AQ has been found.
c
c   BEGIN
c
c Section 1:
c   IF the attitude type is 'none' THEN
c      SET the Attitude Solution Field in the Science Record to 0.
c      RETURN.
c   ENDIF
c   SET ADT to Science_Record Cobetrieve_Header Binary Time
c   IF the attitude type is simulated THEN
c      CALL FUT_SIM_QUAT to determine the COBE celestial equatorial pointing
c      and UTime.
c         This function extracts the year month, day, hour, minute and
c         second and calls TTAG to convert the time to T, the number of
c         seconds since 05/24/68, 00:00 UT.  It then calls ASPECT with input
c         T to get the Attitude matrix A describing S/C orientation at time T.
c         Next it calls M2Q with input A to get Q, the associated quaternion.
c         The -A(I,1) is the FIRAS pointing in the CELESTIAL Equatorial frame.
c      CALL ORBIT to get COBE's geocentric position at time T.
c      CALL ORBIT twice more to get COBE's position at times T +/- dT.
c      COMPUTE COBE velocity, assuming circular orbit.
c   ELSEIF ( Attitude_Type is real attitude ) THEN
c      IF first time THEN
c         CALL FUT_OPEN_ATT_ARCV to open the attitude data set for the input
c              time range and pass back the unit number for reading data.
c         CALL FUT_OPEN_ORBIT_ARCV to open the orbit data set for the input time
c              range and pass back the unit number.
c      ENDIF
c      CALL UAX_GET_ATTITUDE to obtain the attitude quaternion for the current
c                            science record.
c      SET the Attitude Solution Field in the Science Record to the type of
c          attitude obtained.
c      CALL AUT_QTR_TO_ROT to convert the quaternion to the Attitude matrix.
c      CALL UOX_GET_ORBIT to get orbit nr, COBE position and velocity.
c      SET the orbit number in the science record to the orbit returned.
c      CONVERT position and velocity from m and m/sec into km and km/sec .
c   ENDIF
c
c Section 2:
c   CALL GST to get Greenwich Sidereal Time angle in radians.
c   COMPUTE COBE Terrestrial position from Celestial equatorial position by
c           rotation through GST angle.
c   CALL NORM to make the COBE Terrestrial position a unit vector.
c   COMPUTE the geocentric COBE distance and speed.
c   COMPUTE the COBE orbit-plane-normal unit vector.
c   COMPUTE the FIRAS beam unit vector in the celestial equatorial frame by
c           extracting the negative of column 1 of the Attitude matrix.
c   CALL XCC_Q_to_E to convert the FIRAS beam unit vector to the ecliptic frame
c           (retaining the equatorial version).
c   CALL XFM_Locate_Object to get Moon and Sun positions in the ecliptic frame
c           and convert to km for computations with COBE position.
c   CALL XFM_Locate_Object before and after current time and compute the
c           barycentric velocity of the Earth from these.
c   CONVERT Sun and Moon vectors to the ecliptic frame.
c   COMPUTE the COBE_Moon vector by (Moon_Position - COBE_Position).
c   COMPUTE the COBE_Sun vector by (Sun_Position - COBE_Position).
c   COMPUTE the COBE_Moon_Dist by finding the length of the COBE_Moon vector.
c
c Section 3:
c   CALL FIRAS_PIXNO to get the FIRAS Pixel Number from the ecliptic coords. of
c           the FIRAS beam vector.
c   WRITE the FIRAS beam vector.
c   CONVERT the FIRAS beam vector in the cel. eq. frame to spherical coordinates
c           to get RA and DEC.
c   CALL FIRAS_PIXNO to convert COBE Terrestrial position to TERR_PIXEL_NO.
c   COMPUTE TERR_LATITUDE and TERR_LONGITUDE from the COBE Terrestrial position
c           unit vector.
c   COMPUTE the EARTH_LIMB, EARTH_LIMB_AZIMUTH, SUN_ANGLE, MOON_ANGLE, and
c           MOON_AZIMUTH_ANGLE using the Attitude matrix, and computed distances
c           between bodies.
c   COMPUTE the MOON_PHASE by (moon_longitude - sun_longitude) Modulo 2*Pi.
c   COMPUTE SUN_MOON_DIST by finding length of vector (Sun - Moon).
c   WRITE COBE_MOON_DIST found in Section 2.
c   COMPUTE the COBE ALTITUDE, corrected to first order in Earth oblateness, by
c       (COBE_Geo_Dist -
c           Earth_radius*(1. - Earth_oblate*COBE_Geo_Pos(3)/COBE_Geo_Dist)) )).
c   ADD the vectors for barycentric velocity of Earth and geocentric velocity
c           of COBE, then take dot product with FIRAS LOS to obtain the
c           PROJECTED_BARYCENTRIC_VELOCITY.
c   COMPUTE the MCILWAIN_L_PARAM.
c   CALL XCC_E_to_G to convert FIRAS beam vector from eclip. to galactic frame,
c           then convert to spherical coordinates to obtain GALACTIC_LONGITUDE
c           and GALACTIC_LATITUDE.
c   CONVERT FIRAS beam vector (in ecliptic frame) to spherical coordinates to
c           obtain ECLIPTIC_LONGITUDE and ECLIPTIC_LATITUDE.
c   COMPUTE the ORBITAL_PHASE angle from COBE position and velocity.
c   CALL DOT on the FIRAS beam unit vector and COBE velocity vector to get the
c           PROJECTED_GEOCENTRIC_VELOCITY.
c   COMPUTE the SCAN_ANGLE, and SC_ROTATION_ANGLE from the vectors representing
c           Sun position; COBE position, velocity and orbit-plane-normal; and
c           FIRAS beam.
c   SET return status to success or failure.
c
c   RETURN
c   END
c
c-----------------------------------------------------------------------------
c      Modified by:
c
c	  Shirley M. Read
c	  STX
c	  January 18, 1988
c	  REASON: Added computation of all attitude parameters as described
c		  in the new PDL below.  Modified FUT_Sim_Quat to return T, the
c		  time since JD 2440000.5, and A, the rotation matrix which
c		  describes the S/C orientation.  These are now used for
c		  computations in FUT_Attitude.  Also, modified for interface
c		  with the FUT_Error condition handler.  Added error checking
c		  and calls to Lib$Signal.  Removed the FUT_Sim_Pixno since
c		  the celestial equatorial pointing can now be used to get the
c                 FIRAS pixel number.
c
c	  Shirley M. Read
c	  STX
c	  May, 1988
c	  REASON: Added the capability to obtain real attitude as well as
c		  simulated attitude by opening the attitude archive and
c		  calling UAX_Get_Attitude.  The attitude_type, start GMT
c		  and stop GMT are added to the calling sequence.  The type
c		  of attitude solution is entered in the science record.
c
c	  Shirley M. Read
c	  STX
c	  July, 1988
c	  REASON: Added the capability to obtain the orbit number by opening
c		  the orbit data archive and calling UOX_Get_Orbit.
c
c	  R. Kummerer
c	  STX
c	  Sept 19, 1988
c	  REASON: Moon phase calculation falling outside -PI to PI
c		  range.
c
c	  F. Shuman
c	  STX
c	  Oct  4, 1988
c	  REASON: Correct the method for keeping the moon phase in the range
c		  -PI to PI.
c
c	  F. Shuman
c	  STX
c	  Oct 11, 1988
c	  REASON: Write ten additional quantities per SPR 2025, SER 2374:
c	PROJECTED_BARYCENTRIC_VELOCITY  ECLIPTIC_LONGITUDE  PITCH_ANGLE
c	McILWAIN_L_PARAM                ECLIPTIC_LATITUDE   SCAN_ANGLE
c	GALACTIC_LONGITUDE              ORBITAL_PHASE       SC_ROTATION_ANGLE
c	GALACTIC_LATITUDE
c
c	  R. Kummerer
c	  STX
c	  Oct 15, 1988
c	  REASON: SPR 2624, zero the attitude section of raw science record
c	          prior to filling the attitude.
c
c	  1988 Oct 25   F. Shuman,  STX
c	  REASON: Replace PITCH_ANGLE with PROJECTED_GEOCENTRIC_VELOCITY.
c
c	  1989 Jul 21   F. Shuman,  STX
c	  REASON: SPR 3963, Attitude=NONE filled when COARSE was selected.
c	          This can happen either because, on the one hand the Att. or
c	          Orbit Arcv Open fails, or on the other hand, because the
c	          UAX/UOX calls fail from lack of data.  The two cases need to
c	          be distinguished by the calling program, because in the latter
c	          case there may be a temporary absence of data which will clear
c	          up later in the timerange; but if the Arcvs can't be opened,
c	          then there's not going to be any data, and we can abandon the
c	          attempt.  The status FUT_OPENATTORB was added to the MSG file
c	          in order to inform the calling program whether the Arcvs were
c	          successfully opened.
c
c	  1989 Oct 25   R. Kummerer,  STX
c	  REASON: SPR 4825, Improper signalling of error condition from
c		  UAX_GET_ATTITUDE.
c
c	  1989 Nov 29   R. Kummerer,  STX
c	  REASON: SPR 5177, Use collect time instead of transmit time as time
c		  from which attitude is calculated.
c
c	  1990 Jan 31   F. Shuman,  STX
c	  REASON: SPR 5278, Separate open error conditions for Attitude Archive
c		and Orbit Archive, so that when one of these archives fails to
c		open, we can still recover quantities computed solely from the
c		other.
c
c	  1990 Apr 4    R. Kummerer,  STX
c	  REASON: SPR 6629, Mark orbit and attitude archives as being
c		successfully opened so the simulated attitude gets calculated.
c
c	  1990 May 2    R. Kummerer,  STX
c	  REASON: SPR 6746, Remove obsolete and unused UAX message symbols.
c		Obsolete UAX message symbols were generating link errors
c		FIRAS facilities using FUT_ATTITUDE.
c
c	  1990 Sep 25    F. Shuman,  STX
c	  REASON: SPR 7449, Correct the "Celestial" unit vector, which is
c		supposed to be equatorial coordinates, but is in ecliptic
c		instead.  RA and dec fields are OK.
c
c         1990 Nov 29    N. Gonzales, STX
c         REASON: SPR 7811, Add status check for UAX_Success_Non_Definitive.
c
c	  1991 Jan 30    F. Shuman,  STX
c	  REASON: SPR 8037, Correct ~1/3 deg error in conversion to galactic
c		coordinates caused by epoch discrepancy (1950 was being used;
c		2000 now used).  Calls to UOE/UCB routines ETOC, CTOE, CTOG
c		replaced with corresponding XCC routines, XCC_E_to_Q, etc.
c
c	  1991 Apr 19    F. Shuman,  STX
c	    SPR 8013, Remove signalling of UAX msgs--they are useless.
c
c	  1991 Apr 19    F. Shuman,  STX
c	    SPR 8386, Call FUT_VABSAAFlag to fill new field TERR_RAD_BYTE.
c	    Add call to FUT_ATTOFFSET as stub for possible future FIRAS
c	    beam offset from spacecraft axis.
c
c	  1991 Jul 09    who? what? where?
c	    SPR 8764, Change of field name in attitude.rdl from CELESTIAL(3) to
c	    EQUATORIAL(3).
c
c	  1991 Aug 30    F. Shuman,  STX
c	    SPR 8971, 8976, Corrections and bringing up-to-date in response to
c	    code review by Fred Patt.
c
c	  1991 Oct 17    F. Shuman,  STX
c	    SPR 9168, 'First call' logic was faulty.
c
c	  1991 Dec 13    F. Shuman,  Hughes STX
c	    SPR 9345, Correct the Earth oblateness correction.
c
c	  1992 Jun 25    F. Shuman,  Hughes STX
c	    SER 9040, Correct the Earth Limb Angle and Azimuth for oblateness.
c-----------------------------------------------------------------------------

	Implicit None

c      Passed arguments and records.

	Dictionary 'nfs_sdf'
	Record/nfs_sdf/  sci_rec	! NFS science record

	Integer*4        attitude_type	! Type of attitude requested
	Character*14     gmt_start	! Start time for attitude
	Character*14     gmt_stop	! Stop time for attitude

c      Include files

	Include '(FUT_Params)'
	Include '(UAX_Inc)'
	Include '(UOX_Time_Str)'
	Record /UOX_Time_Str/ Time_Rec
	Include '(UOX_Ephem_Str)'
	Record /Ephem/ Ephem_Rec

c      Functions

	Integer*4	FUT_Sim_Quat
	Integer*4	FUT_Open_Att_Arcv
	Integer*4	UAX_Get_Attitude
	Integer*4	FUT_Open_Orbit_Arcv
	Integer*4	UOX_Get_Orbit
	Integer*4	FUT_AttOffset
	Integer*4	FUT_VABSAAFlag
	Real*8		AUT_ADT2T68
	Real*4		GST		! facility UOE; filename UCB_GST.FOR
	Integer*4	CT_GMT_To_Binary
	Integer*4	XFM_Locate_Object

	External	FUT_NORMAL
	External	FUT_ATTOPNERR
	External	FUT_ORBOPNERR
	External	UAX_SUCCESS
	External	UAX_SUCCESS_WITH_EXTRAPOLATION
	External	UAX_SUCCESS_NON_DIRBE
	External        UAX_SUCCESS_NON_DEFINITIVE
	External	UOX_SUCCESS

c	Constants and Parameters

	Real*4		Pi
	Parameter	  ( Pi = fac_pi )
	Real*4		Earth_radius		   !Equat. radius of Earth in km
	Parameter	  ( Earth_radius = 6378.164 )
	Real*4		Earth_oblate		   !Oblateness = 1 - (Polar/
	Parameter	  ( Earth_oblate = 0.0033529 )	!             Equat rad)
	Real*4		mjupiter		   !Jupiter mass in Solar masses
	Parameter	  ( mjupiter = 9.546E-4 )
	Real*4		sjup			   !Sun speed due to Jupiter
	Parameter	  ( sjup = 0.013 )	   ! (km/sec)
	Real*4		AU			   !Astronomical Unit in km
	Parameter	  ( AU = 1.49597870E+8 )
	Real*4		epoch			   !..for coordinate conversions
	Parameter	  ( epoch = 2000. )
	Integer*4	dadtsun			   !..for solar velocity calc.
	Parameter	  ( dadtsun = 5941 )	   ! =~ 1 synodic mo in ADT(hi)
	Real*4		dtsun			   !..for solar velocity calc.
	Parameter	  ( dtsun = dadtsun*4.294967296D+2 )	! multiplier =
c	2^32 * 10^-7 = value of a unit in the high-signif. half of an ADT
	Character*14	GMTcen			   !..for solar velocity calc.
	Parameter	  ( GMTcen = '90112000000000' )	!..Midpt time of FIRAS

c	Local variables

	Integer*4	status		!Dummy status variable
	Integer*4	zero / 0 /

	Integer*4	att_lun         !Logical unit for attitude archive
	Integer*4	orb_lun         !Logical unit for orbit archive
	Integer*4	type_code       !Type of attitude to be returned by UAX
	Real*4		attitude(4)     !Attitude quaternion returned from UAX
	Real*4		att_var(10)     !Attitude variance returned from UAX
	Real*4		error_statistics !Error statistics from UAX
	Logical*1	first / .True. /     !First call flag
	Logical*1	attopntried/.False./ !Have tried to open attitude
	Logical*1	attopen / .False. /  !Flag: Attitude arcv open success

	Integer*4	pix_no		!Pixel number
	Integer*4	adt(2)		!Binary time tag of IFG
	Integer*4	adtbak(2)	!Binary time - dadtsun
	Integer*4	adtfor(2)	!Binary time + dadtsun
	Character*14	adtcen		!Midpt time of FIRAS for solar vel calcn
	Integer*4	i, j		!Loop counters

	Real*8		t		!Time (10^-7 sec) since reference of
					! JD 2400000.5 = 1858 Nov 17 0h UT
	Real*8		dt		!Time difference (sec)
	Real*4		q(4)		!Attitude quaternion returned by Simquat
	Real*4		amat(3,3)	!SC orientation matrix.  If vector r(3)
		! is in the SC frame, then   r'(i) = (sum j):[Amat(j,i) r(j)]
		! is r in (celes) eqtl. frame.  Amat is orthogonal, so its
		! inverse = its transpose, and  r(i)=(sum j):[Amat(i,j) r'(j)].
		! -Amat(1,i) is the Firas pointing vector in the celestial
		! equatorial frame.
	Real*4		aoffset(3,3)	!SC orientation matrix rotated to adjust
					! its 'x' direction to the FIRAS beam.
	Real*4		sidtim		!Sidereal "time", in radians
c
c	Note: All vectors are represented in Cartesian coordinates.
c	  Reference Frame origin, Frame orientation, Units are given.
c	  "Terr", "Terrestrial" denote Earth-co-rotating frame.
c	  "Spacecraft frame" co-rotates with COBE.
c
	Real*4		beam_eq(3)	   !Firas beam unit-vector--
					   ! COBE, (cel) equatorial, Dimen'less
	Real*4		beam(3)		   !Firas beam unit-vector--
					   ! COBE, ecliptic, Dimensionless
	Real*4		COBE_Geo_Pos(3)	   !COBE position--time t
	Real*4		COBE_Pos_Plus(3)   !COBE position--time t + dt
	Real*4		COBE_Pos_Minus(3)  !COBE position--time t - dt
					   ! Geocenter, (cel) equatorial, km
	Real*4		u(3)		   !COBE position unit-vector--
					   ! Geocenter, (terr) equatl, Dim'less
	Real*4		lon		   !COBE terrestrial longitude (rad)
	Real*4		lat		   !COBE terrestrial latitude (rad)
	Real*4		londeg		   !COBE terrestrial longitude (deg)
	Real*4		latdeg		   !COBE terrestrial latitude (deg)
	Byte		flag		   !Terrestrial radiation flag--
					   !  0=Not checked, 1=Outside all belts
					   !  2=VAB(N), 4=VAB(S), 8=SAA
	Real*4		Celpos_Sun(3)      !Position of Sun--time t
					   ! Geocenter, (cel) equatorial, km
	Real*4		Celpos_Moon(3)     !Position of Moon--
					   ! Geocenter, (cel) equatorial, km
	Real*4		Eclpos_Sun(3)      !Position of Sun--
	Real*4		Eclpos_sun_plus(3) !                 time t+dt
	Real*4		Eclpos_sun_minus(3)!                 time t-dt
					   ! Geocenter, ecliptic, km
	Real*4		Eclpos_Moon(3)     !Position of Moon--
					   ! Geocenter, ecliptic, km
	Real*4		COBE_Sun(3)        !Vector from COBE to Sun--
					   ! COBE, (cel) equatorial, km
	Real*4		COBE_Moon(3)	   !Vector from COBE to Moon--
					   ! COBE, (cel) equatorial, km
	Real*4		Moon_SC(3)	   !Vector from COBE to Moon--
					   ! COBE, Spacecraft frame, km
	Real*4		Earth_SC(3)	   !Vector from COBE to Earth--
					   ! COBE, Spacecraft frame, km
	Real*4		COBE_Geo_Vel(3)	   !COBE velocity--
					   ! Geocenter, (cel) equatorial, km/sec
	Real*4		Orb_Norm(3)	   !COBE orbit normal unit vector--
					   ! Geocenter, (cel) equatl, Dimen'less
	Real*4		Earth_Helio_Vel(3) !Earth velocity--
					   ! Sun, (cel) equatorial, km/sec
	Real*4		eclnor(3)	   !Ecliptic North unit vector--
					   ! Sun, (cel) equatorial, Dimen'less
	Real*4		vjup(3)		   !Sun velocity due to Jupiter
					   ! Sun, ecliptic, km/sec
	Real*4		vect(3)		   !For temporary use
	Real*4		cvec(3)		   !Basis vectors of a Sun-
	Real*4		bvec(3)		   !      N eclip pole-
	Real*4		avec(3)		   !      FIRAS eclip plane ascent frame
					   ! COBE, (cel) equatorial, Dimen'less

	Real*4	        sumsq		!Sum of squares (temporary)
	Real*4		COBE_Geo_Dist	!COBE-Earth Center distance (km)
	Real*4		COBE_Geo_Speed	!COBE Geocentric speed (km/sec)
	Integer*4	N		!Planet number (in UOE/UCB lib)
	Real*4		COBE_Moon_Dist  !COBE-Moon distance (km)
	Real*4		COBE_Sun_Dist   !COBE-Sun distance (km)
	Real*4		Earth_Sun_Dist  !Earth-Sun distance (km)
	Real*4		distance	!Distance, general storage variable
	Real*4		temp, cros	!Various temporary variables
	Real*4		dot, bdott		!Various dot products
	Real*4		tau, stau, ctau		!Horizon tangent angle from COBE
						! to a sphere of rad=Earth eqtr,
						! & its sin & cos
	Real*4		Moon_lon	!Moon longitude in ecliptic frame
	Real*4		Sun_lon		!Sun longitude in ecliptic frame
	Real*4		Sun_Moon_londif !Difference of ecliptic longitudes,
					! Sun-Moon, brought into range -pi,+pi
	Real*4		Moon_phase	!Moon phase
	Integer*4	npix		!Terrestrial FIRAS pixel number

c1  Initialize, annihilating the attitude portion of the science at the start.
c1   If no attitude is wanted, then we're done here.

	Save first, attopntried, attopen, vjup, eclnor

	Call Lib$MovC5 ( 0, , 0, fac_attit_size, sci_rec.attitude )
	sci_rec.attitude.solution = fac_none
	sci_rec.attitude.pixel_no = -1
	sci_rec.attitude.terr_pixel_no = -1

	If ( attitude_type .Eq. fac_none ) Then
	  FUT_Attitude = %Loc(FUT_Normal)
	  Return
	End If

c1  Set the time tag to the science record time.

	If (sci_rec.collect_time.midpoint_time(1) .Ne. 0 .Or.
     &      sci_rec.collect_time.midpoint_time(2) .Ne. 0) Then

	  adt(1) = sci_rec.collect_time.midpoint_time(1)
	  adt(2) = sci_rec.collect_time.midpoint_time(2)

	Else

	  adt(1) = sci_rec.ct_head.time(1)
	  adt(2) = sci_rec.ct_head.time(2)

	End If

c1  Process according to the type of attitude selected.

	If ( attitude_type .Eq. fac_simulated ) Then

c1  Determine the COBE attitude quaternion. FUT_Sim_Quat uses the input adt to
c1   extract year, month, day, hour, minute and second.  It then calls Ttag to
c1   convert the time to T, the number of seconds since 05/24/68, 00:00 UT.  It
c1   next calls Aspect, which uses T to get the rotation matrix, Amat, giving
c1   the spacecraft orientation at time T.  Finally, it calls M2Q which uses
c1   Amat to get Q, the associated quaternion.  FUT_Sim_Quat returns Amat and
c1   T, which will be needed for later computations; note that the Amat it
c1   returns is in the vcel=Amat.vSC sense and so has to be transposed in order
c1   to bring it into agreement with the Amat returned by AUT_Qtr_to_Rot.
c1  CALL ORBIT to get the COBE position, COBE_Geo_Pos(3), in km, referred to
c1    Earth center in the Celestial Equatorial frame.

	  status = FUT_Sim_Quat( adt, Q, Amat, T )
	  If ( status .Ne. %Loc(FUT_NORMAL) ) Then
	    FUT_Attitude = status
	    Return
	  Else
	    temp = Amat(1,2)
	    Amat(2,1) = temp
	    Amat(1,2) = Amat(2,1)
	    temp = Amat(2,3)
	    Amat(3,2) = temp
	    Amat(2,3) = Amat(3,2)
	    temp = Amat(3,1)
	    Amat(1,3) = temp
	    Amat(3,1) = Amat(1,3)
	    sci_rec.attitude.solution = fac_simulated
	  End If

	  dt = 2.D2
	  Call Orbit ( t, COBE_Geo_Pos )
	  Call Orbit ( t + dt, COBE_Pos_plus )
	  Call Orbit ( t - dt, COBE_Pos_minus )
c1
c1  To compute the COBE velocity, assume circular orbit.  Use the orbital period
c1   taken by the ORBIT routine, 102.7 min = 6162 sec.
	  Do i=1,3
	    COBE_Geo_Vel(i) = ( COBE_Pos_plus(i) - COBE_Pos_minus(i)) * pi /
     &                         ( 6162. * Sin(2*pi*dt/6162.) )
	  End Do

	  attopen = .True.

	Else

c1  Real attitude is requested.  Open the attitude archive using the
c1   start and stop times.  Call UAX_Get_Attitude to get the attitude
c1   quaternion.  Only one set of vectors is needed for one science record.
c1   Open the orbit archive and call UOX_Get_Orbit to get the orbit number.

	  If ( attitude_type .Eq. fac_fine_with_Dirbe ) Then
	    type_code = UAX_C_Dirbe
	  Else If ( attitude_type .Eq. fac_definitive ) Then
	    type_code = UAX_C_Definitive
	  Else
	    type_code = UAX_C_Sunearth
	  End If

	  If ( .Not. attopntried ) Then
	    attopntried = .True.
	    status = FUT_Open_Att_Arcv (type_code, gmt_start, gmt_stop, att_lun)
	    If ( status .Eq. %Loc(FUT_NORMAL) ) Then
	      attopen = .True.
	    Else
	      attopen = .False.
	      FUT_Attitude = %Loc(FUT_ATTOPNERR)
	    End If

	    time_rec.number_of_points = 1
	    status = FUT_Open_Orbit_Arcv ( gmt_start, gmt_stop, orb_lun )
	    If ( status .Ne. %Loc(FUT_NORMAL) ) Then
	      FUT_Attitude = %Loc(FUT_ORBOPNERR)
	      Return
	    End If
	  End If

	  Do i = 1, 2
	    time_rec.initial_time(i) = adt(i)
	  End Do

	  If (attopen) Then
	    status = UAX_Get_Attitude ( time_rec, att_lun, type_code,
     &                                  attitude, att_var, error_statistics, )
	    If ( status .Eq. %Loc(UAX_Success) ) Then
	      sci_rec.attitude.solution = attitude_type
	    Else If(status .Eq. %Loc(UAX_Success_With_Extrapolation)) Then
	      sci_rec.attitude.solution = attitude_type
	    Else If(status .Eq. %Loc(UAX_Success_Non_Dirbe)) Then
	      sci_rec.attitude.solution = fac_fine_without_Dirbe
	    Else If(status .Eq. %Loc(UAX_Success_Non_Definitive)) Then
	      sci_rec.attitude.solution = fac_coarse
	    Else   ! Bad status return
	      sci_rec.attitude.solution = fac_none
	      FUT_Attitude = status
	      Return
	    End If

c1  Convert the quaternion to the attitude rotation matrix.  Offset the matrix
c1   to bring the spacecraft negative spin axis to the FIRAS beam direction.

	    Call AUT_Qtr_to_Rot ( attitude, amat )
	    status = FUT_AttOffset(amat, aoffset)
	  End If

c1  Convert the science record binary time to time since 5/24/68 in seconds;
c1   AUT_ADT2T68 uses JD 2440000.5 (= 1968 May 24, 0h UT) as reference, with
c1   sec as unit, while taking 'leap seconds' into account; this is consistent
c1   with the E. L. "Ned" Wright routines in UOE/UCB.

	  t = AUT_ADT2T68( adt )

c1  Get the orbit number and COBE position and velocity.  Because these are
c1   returned in m and m/sec, we multiply by 1.E-3 to convert to km and km/sec.

	  status = UOX_Get_Orbit ( time_rec, orb_lun, ephem_rec )
	  If ( status .Eq. %Loc(UOX_Success) ) Then
	    sci_rec.ct_head.orbit = ephem_rec.orbit_num
	    COBE_Geo_Pos(1) = ephem_rec.x_pos * 1.E-3
	    COBE_Geo_Pos(2) = ephem_rec.y_pos * 1.E-3
	    COBE_Geo_Pos(3) = ephem_rec.z_pos * 1.E-3
	    COBE_Geo_Vel(1) = ephem_rec.vx_vel * 1.E-3
	    COBE_Geo_Vel(2) = ephem_rec.vy_vel * 1.E-3
	    COBE_Geo_Vel(3) = ephem_rec.vz_vel * 1.E-3
	  Else
	    Call LIB$Signal(status)
	    FUT_Attitude = status
	    Return
	  End If
	End If

c2
c2 Preliminary calculations...
c2

c2  Convert COBE Celestial equatorial position to Terrestrial (frame co-rotates
c2   with Earth) by fetching the "Greenwich sidereal time," which GST returns
c2   as an angle in radians, and normalize to form a unit vector.
c2   Conversion calculation swiped from E. L. Wright's routine, TERRESTIAL
c2   (sic), to be found in UOE/UCB.  Routines GST and Norm are also in UOE/UCB.

	sidtim = GST(t)
	u(1) =  Cos(sidtim)*COBE_Geo_Pos(1) + Sin(sidtim)*COBE_Geo_Pos(2)
	u(2) = -Sin(sidtim)*COBE_Geo_Pos(1) + Cos(sidtim)*COBE_Geo_Pos(2)
	u(3) =              COBE_Geo_Pos(3)
	Call Norm ( u )

c2  Compute lon and lat, the terrestrial longitude and latitude, from the
c2   Terrestrial position unit vector of COBE by converting to spherical
c2   coordinates.  Rather than using asin(u(3)) for computing latitude, use
c2   more robust atan2 form.  When abs(u(3)) is close to 1, asin is unreliable.

	temp = SQRT(u(1)**2 + u(2)**2)
	lon = ATan2(u(2), u(1))
	lat = ATan2(u(3), temp)
	londeg = lon*180./Pi
	latdeg = lat*180./Pi

c2  Find the geocentric COBE distance and speed.

	sumsq = COBE_Geo_Pos(1)**2 + COBE_Geo_Pos(2)**2 + COBE_Geo_Pos(3)**2
	COBE_Geo_Dist = SQRT(sumsq)
	sumsq = COBE_Geo_Vel(1)**2 + COBE_Geo_Vel(2)**2 + COBE_Geo_Vel(3)**2
	COBE_Geo_Speed = SQRT(sumsq)

c2  Compute the COBE orbit plane normal, L = R x V, and normalize it.

	Call CROSS(COBE_Geo_Pos, COBE_Geo_Vel, Orb_Norm)
	Call Norm ( Orb_Norm )

c2  Get the FIRAS look direction in the (cel.) equatorial frame, given by the
c2   - Amat(1,i) row of the matrix, since this is just Amat applied to (-1,0,0)
c2   [in the v.Amat sense], the COBE-frame coordinates of FIRAS pointing.

	If (attopen) Then
	  Do i=1,3
	    beam_eq(i) = -aoffset(1,i)
	  End Do

c2  Convert the FIRAS beam vector from (celes.) equatorial frame ('Q'), epoch
c2   2000, to ecliptic ('E'), epoch 2000.

	  Call XCC_Q_to_E (beam_eq, epoch, beam)
	End If

c2  CALL XFM_Locate_Object to get moon and sun positions in the ecliptic frame;
c2   then convert to km for computations with COBE position.

	Call XFM_Locate_Object ( 'MOON', adt, eclpos_moon )

	adtbak(1) = adt(1)
	adtbak(2) = adt(2) - dadtsun
	adtfor(1) = adt(1)
	adtfor(2) = adt(2) + dadtsun

	Call XFM_Locate_Object ( 'SUN', adtbak, eclpos_sun_minus )
	Call XFM_Locate_Object ( 'SUN', adt,    eclpos_sun )
	Call XFM_Locate_Object ( 'SUN', adtfor, eclpos_sun_plus )

	Do i = 1, 3
	  eclpos_moon(i)       = AU * eclpos_moon(i)
	  eclpos_sun_minus(i)  = AU * eclpos_sun_minus(i)
	  eclpos_sun(i)        = AU * eclpos_sun(i)
	  eclpos_sun_plus(i)   = AU * eclpos_sun_plus(i)
	End Do

c2  Compute the heliocentric velocity of the Earth using Sun positions forward
c2   and backward in time and assuming circular orbit.
	Do i=1,3
	  Earth_Helio_Vel(i) = (-eclpos_sun_plus(i) + eclpos_sun_minus(i)) *
     &           pi / ( fac_year * Sin(2*pi*dtsun/fac_year) )
	End Do

c2  Compute the Earth-Sun distance.
	sumsq = 0.
	Do i = 1, 3
	  sumsq = sumsq + eclpos_Sun(i)**2
	End Do
	Earth_Sun_Dist = SQRT(sumsq)

c2  Convert Sun and Moon vectors from ecliptic ('E') to (cel.) equatorial ('Q')
c2    frame.
	Call XCC_E_to_Q (eclpos_Sun, epoch, celpos_Sun)
	Call XCC_E_to_Q (eclpos_Moon, epoch, celpos_Moon)

c2  Compute the COBE_Sun vector and COBE_Sun_Dist.
	sumsq = 0.
	Do i = 1, 3
	  COBE_Sun(i) = celpos_Sun(i) - COBE_Geo_Pos(i)
	  sumsq = sumsq + COBE_Sun(i)**2
	End Do
	COBE_Sun_Dist = SQRT(sumsq)

c2  Compute the COBE_Moon vector and COBE_Moon_Dist.
	sumsq = 0.
	Do i = 1, 3
	  COBE_Moon(i) = celpos_Moon(i) - COBE_Geo_Pos(i)
	  sumsq = sumsq + COBE_Moon(i)**2
	End Do
	COBE_Moon_Dist = SQRT(sumsq)

c2  Find the heliocentric position of Jupiter at 90112=center time of FIRAS
c2   mission, so we can estimate the velocity of the Sun due to Jupiter.
c2   Multiply the orbital speed of Jupiter by its mass ratio with the Sun to get
c2   the Jupiter-induced speed of the S.S. Barycenter.
	If (first) Then
	  first = .False.
	  status = CT_GMT_To_Binary(GMTcen, adtcen)
	  Call XFM_Locate_Object ( 'JUP', adtcen, avec )
	  Do i = 1, 3
	    avec(i) = AU * avec(i) - eclpos_sun(i)
	  End Do
	  Call Norm(avec)
	  vjup(1) =-sjup*avec(2)
	  vjup(2) = sjup*avec(1)
	  vjup(3) = 0.

c2  Convert the unit z-vector in the ecliptic system to equatorial frame.
	  vect(1) = 0.
	  vect(2) = 0.
	  vect(3) = 1.
	  Call XCC_E_to_Q (vect, epoch, eclnor)
	End If

c2
c2 ...Finished with preliminary calculations.
c2


c3  Determine FIRAS pixel number by calling FIRAS_Pixno.

	If (attopen) Then
	  Call Firas_pixno( beam, pix_no)
	  sci_rec.attitude.pixel_no = pix_no

c3  Store the FIRAS look direction (cel. equatorial frame)

	  Do i=1,3
	    sci_rec.attitude.equatorial(i) = beam_eq(i)
	  End Do

c3  Compute RA and DEC from the FIRAS beam unit vector by converting equatorial
c3   from Cartesian to spherical coordinates.  Rather than using asin(beam(3))
c3   for computing declination angle, use more robust atan2 form.  When
c3   Abs(beam(3)) is close to 1, asin is unreliable.

	  sci_rec.attitude.ra = Nint( 1.E4 * ATan2(beam_eq(2), beam_eq(1)) )
	  temp = SQRT(beam_eq(1)**2 + beam_eq(2)**2)
	  sci_rec.attitude.dec = Nint( 1.E4 * ATan2(beam_eq(3), temp) )
	End If

c3  CALL FIRAS_PIXNO to convert unit Terrestrial vector to TERR_PIXEL_NO

	Call FIRAS_Pixno ( u, npix )
	sci_rec.attitude.terr_pixel_no = npix

c3  Compute TERR_LONGITUDE and TERR_LATITUDE from the previously computed
c3   Terrestrial longitude and latitude of COBE.

	sci_rec.attitude.terr_longitude = Nint( 1.E4 * lon )
	sci_rec.attitude.terr_latitude = Nint( 1.E4 * lat )

c3  Earth limb computation.  Limb angle is minimum angle from FIRAS beam to
c3   horizon.  For the spherical approximation to Earth figure, this was
c3   180 deg - tau(=nadir-horizon angle) - (angle from FIRAS beam to zenith).
c3   Corrected to first order in Earth oblateness, it is as shown below,
c3   assuming that the beam vector is nearly aligned with COBE_Geo_Pos.

	If (attopen) Then
	  stau = Earth_radius / COBE_Geo_Dist
	  ctau = Sqrt(1 - stau**2)
	  tau = ASin ( stau )

	  dot = 0.
	  Do i = 1, 3
	    dot = dot + beam_eq(i)*COBE_Geo_Pos(i)
	  End Do
	  dot = dot/COBE_Geo_Dist
	  cros = SqRt(1 - dot**2)
	  bdott = stau*cros - ctau*dot - Earth_oblate*stau**4*u(3)**2/ctau

	  sci_rec.attitude.earth_limb = Nint(1.E4 * ACos(bdott))

c3  Computation of the EARTH_LIMB_AZIMUTH (azimuthal angle in SC y-z plane of
c3   the minimum FIRAS beam-to-Earth limb angle).
c3  Convert COBE position vector, COBE_Geo_Pos, from Celestial Equatorial frame
c3   to SC (COBE-oriented) frame by calling AUT_Vect_Rot, using the Amat matrix.
c3   (Minus this vector is taken, because COBE_Geo_Pos is from Earth to COBE.)
c3   The angle from z toward y in the y-z plane of the S/C frame is
c3   EARTH_LIMB_AZIMUTH.  Not corrected for Earth oblateness.

	  Call AUT_Vect_Rot(amat, COBE_Geo_Pos, Earth_SC)
	  Do i = 1, 3
	    Earth_SC(i) = -Earth_SC(i)
	  End Do
	  sci_rec.attitude.earth_limb_azimuth =
     &                    Nint(1.E4 * ATan2 ( Earth_SC(2), Earth_SC(3) ))
	End If

c3  Sun_Angle computation:
c3  Take dot product of beam unit vector with COBE-Sun vector and divide by
c3   COBE_Sun_Dist.  Ratio always far from +/- 1, so may as well use acos.
c3   If we can't get orbit data, we replace the COBE-Sun vector with the
c3   Earth-Sun vector, a change in direction of approx. 1/20000 radian.

	If (attopen) Then
	  dot = 0.
	  Do i = 1, 3
	    dot = dot + beam_eq(i) * COBE_Sun(i)
	  End Do
	  sci_rec.attitude.sun_angle = Nint( 1.E4 * ACos(dot/COBE_Sun_Dist) )
	End If

c3  Moon_Angle computation.
c3  Dot product of beam unit vector with COBE_Moon, divided by COBE_Moon_Dist.
c3   Ratio can be close to +/- 1, so use more robust atan2 form to get angle.

	If (attopen) Then
	  dot = 0.
	  Do i = 1, 3
	    dot = dot + beam_eq(i)*COBE_Moon(i)
	  End Do

	  temp = SQRT(Abs(COBE_Moon_Dist**2 - dot**2))
	  sci_rec.attitude.moon_angle = Nint( 1.E4 * ATan2(temp, dot) )

c3  Computation of the Moon azimuth angle.
c3  Transform Moon position vector from Earth to COBE-oriented frame.  Find the
c3   angle from z toward y in the Y-Z plane; this is the MOON_AZ_ANGLE.

	  Call AUT_Vect_Rot(amat, COBE_Moon, Moon_SC)
	  sci_rec.attitude.moon_az_angle =
     &                       Nint( 1.E4 * ATan2(Moon_SC(2), Moon_SC(3)) )
	End If

c3  Computation of the Moon Phase.
c3  Find the difference in Sun and Moon longitudes (which will range from
c3   -2*pi to +2*pi), and use AMOD to bring it into the -pi to +pi range.

	Sun_lon = ATan2 ( Eclpos_Sun(2), Eclpos_Sun(1) )
	Moon_lon = ATan2 ( Eclpos_Moon(2), Eclpos_Moon(1) )

	sci_rec.attitude.moon_phase =
     &             Nint( 1.E4*(amod( Moon_lon-Sun_lon+3.*Pi, 2.*Pi ) - Pi) )

c3  Compute SUN_MOON_DIST by finding length of (Sun minus Moon) vector.

	sumsq = 0.
	Do i = 1, 3
	  sumsq = sumsq + ( Celpos_Sun(i) - Celpos_Moon(i) )**2
	End Do

	sci_rec.attitude.sun_moon_dist = SQRT(sumsq)

c3  Write the COBE_MOON_DIST found previously.

	sci_rec.attitude.cobe_moon_dist = COBE_Moon_Dist

c3  Compute the altitude of COBE, corrected to first order in Earth oblateness.

	sci_rec.attitude.altitude = Nint( 10.*(COBE_Geo_Dist - Earth_radius*
     &                (1. - Earth_oblate*(COBE_Geo_Pos(3)/COBE_Geo_Dist)**2)) )

c3  Compute the solar system barycentric velocity of COBE.
c3  Add the heliocentric velocity of the Earth to the geocentric velocity of
c3   COBE to find the heliocentric velocity of COBE.  Subtract from this the
c3   heliocentric velocity of the Solar System Barycenter, approximated well
c3   enough for our purposes by the effect of Jupiter.  Perform a dot product
c3   with FIRAS beam unit vector to get "PROJECTED_BARYCENTRIC_VELOCITY".

	If (attopen) Then
	  dot = 0.
	  Do i=1,3
	    dot = dot +beam_eq(i)*(Earth_Helio_Vel(i) +COBE_Geo_Vel(i) -vjup(i))
	  End Do
	  sci_rec.attitude.projected_barycentric_velocity = Nint(100.*dot)
	End If

c3  CALL < McIlwain > to get the MCILWAIN_L_PARAM

c3  CALL XCC_E_to_G to convert FIRAS beam vector [ecliptic --> galactic frame].
c3   Convert to spherical coordinates to get galactic long. and lat.  Use the
c3   temporary array VECT, and, for the latitude, again the more robust atan2.

	If (attopen) Then
	  Call XCC_E_to_G (beam, epoch, vect)
	  sci_rec.attitude.galactic_longitude =
     &                          Nint( 1.E4 * ATan2(vect(2),vect(1)) )
	  temp = SQRT(vect(1)**2 + vect(2)**2)
	  sci_rec.attitude.galactic_latitude = Nint( 1.E4*ATan2(vect(3),temp) )

c3  Convert FIRAS beam vector in the ecliptic frame to spherical coordinates to
c3   get ecliptic long. and lat.  Atan2 again instead of asin for the latitude.

	  sci_rec.attitude.ecliptic_longitude =
     &                          Nint( 1.E4 * ATan2(beam(2),beam(1)) )
	  temp = SQRT(beam(1)**2 + beam(2)**2)
	  sci_rec.attitude.ecliptic_latitude = Nint( 1.E4*ATan2(beam(3),temp) )
	End If

c3  Find the geocentric angle from ascending node passage to present position.
c3   This is based on the assumption of a circular orbit.  The position, r, and
c3   velocity, v, vectors are perpendicular and rotate with the COBE orbital
c3   period, T.  They are expressed in the equatorial frame, so the "z" axis
c3   lies along the North pole.  Consider the frame in which r and v are fixed.
c3   In this frame, "z" rotates in a plane parallel to the r-v plane, and when
c3   projected into that plane, the angle it makes from v toward r gives the
c3   orbital phase.  The components of this projection are just dot-products of
c3   unit "z" with unit "r" and "v" vectors, and these are just z components of
c3   r,v unit vectors.  Note that because we are taking only "z" components, we
c3   don't care that u and COBE_Geo_Vel are given in frames rotating relative to
c3   each other about z.

	sci_rec.attitude.orbital_phase =
     &               Nint( 1.E4 * ATan2(u(3), COBE_Geo_Vel(3)/COBE_Geo_Speed) )

c3  Find the PROJECTED_GEOCENTRIC_VELOCITY of COBE by performing a dot product
c3   of COBE's geocentric velocity with FIRAS beam unit vector.

	If (attopen) Then
	  dot = 0.
	  Do i=1,3
	    dot = dot + beam_eq(i)*COBE_Geo_Vel(i)
	  End Do

	  sci_rec.attitude.projected_geocentric_velocity = Nint( 1.E3 * dot )
	End If

c3  Find the scan angle (azimuthal angle, about the Sun direction, of the
c3   skyhorn, with the ascending crossing of the ecliptic plane as reference).
c3   First construct an orthonormal frame A,B,C with C toward the Sun, B in the
c3   direction of the ecliptic z-axis (North eclip. pole), and A = B x C, this
c3   being on the same side of the BC plane as the ascending ecliptic crossing
c3   of the FIRAS beam (if COBE had been launched into a morning northbound,
c3   evening southbound orbit, we would take A = -B x C).  Thus, the scan will
c3   turn about the C axis, the scan angle increasing from 0 in the AC plane
c3   toward the BC plane; the scan angle being given by the A, B components of
c3   FIRAS beam.  Use the orbit plane normal vector, Orb_Norm, to find the
c3   direction of A.

	If (attopen) Then
	  Do i=1,3
	    bvec(i) = eclnor(i)
	    cvec(i) = COBE_Sun(i)
	  End Do
	  Call Norm(cvec)
	  Call CROSS(bvec, cvec, avec)

c3  Use the temporary array VECT to find the A and B components of beam.

	  Do i=1,3
	    vect(i) = 0.
	  End Do
	  Do i=1,3
	    vect(1) = vect(1) + beam_eq(i)*avec(i)
	    vect(2) = vect(2) + beam_eq(i)*bvec(i)
	  End Do

	  sci_rec.attitude.scan_angle = Nint( 1.E4 * ATan2(vect(2),vect(1)) )
	End If

c3  To find the COBE rotation angle, project the Sun vector onto the y-z plane
c3   of the spacecraft and take the angle formed with the z-axis, increasing
c3   toward the y-axis.  Do this by converting the Sun vector in equat. frame
c3   into the SC frame, using the temporary array VECT for the result.

	If (attopen) Then
	  Call AUT_Vect_Rot(amat, Celpos_Sun, vect)
	  sci_rec.attitude.sc_rotation_angle=Nint( 1.E4*ATan2(vect(2),vect(3)) )
	End If

c3  Call FUT_VABSAAFlag to find which radiation belt COBE is in, if any.

	status = FUT_VABSAAFlag(londeg, latdeg, flag)
	sci_rec.attitude.terr_rad_byte = flag

c3  FUT_Attitude always gives the following fields the values shown.

	sci_rec.attitude.pixel_definition = 'q'
	sci_rec.attitude.skymap_index =  5
	sci_rec.attitude.exc_galactic_lat = 0


c3  Done.  Exit.

	FUT_Attitude = %Loc(FUT_Normal)

	Return
	End
