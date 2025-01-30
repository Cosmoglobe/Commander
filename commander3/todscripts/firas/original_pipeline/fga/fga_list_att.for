	SUBROUTINE FGA_LIST_ATT (WRITE_LUN, ATT_REC)

C-----------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_ATT
C
C	PROGRAM DESCRIPTION:
C	  The program will get a formatted listing of all the attitude
C	  info in a given FIRAS IFG or spectrum record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  April 8, 1988
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ATT (WRITE_LUN, ATT_REC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file; if
C					not given assume 6.  Assume the calling
C					routine has opened this list file.
C	  ATT_REC		RECORD	FIRAS Science record
C
C	OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  NONE
C
C	OUTPUT FILES:
C	  Listing file or FOR006
C
C	INCLUDE FILES USED:
C	  FUT_PARAMS.TXT
C
C	SUBROUTINES CALLED:
C	  NONE
C
C-----------------------------------------------------------------------
C   Modifications:
C
C	Update for modified Attitude.rdl record; SER 2374.
C	  R. Kummerer, 1988 Aug 29.
C	  F. Shuman, 1988 Oct 25.
C
C	Alternate pixelization related changes: ALT_PIXEL_NO no longer
C	supported; SER 2370.
C	  R. Kummerer, 1989 Feb 14.
C
C       Alternate pixel information field is enlarged to prevent conversion
C       error during dumps for FCI_SKYCO records. SPR 7538
C         N. Gonzales, 1990 Oct. 27.
C
C       SPR 8372, Updated the code to handle new and revised rdl's.
C         N. Gonzales/STX, September, 12, 1991.              
C-----------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'(FUT_PARAMS)'

	INTEGER		*4	WRITE_LUN
	INTEGER		*4	STATUS

	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L
	REAL		*4	RA
	REAL		*4	DEC
	REAL		*4	TLAT
	REAL		*4	TLONG
	REAL		*4	ELIMB
	REAL		*4	ELIMBAZ
	REAL		*4	SUNANG
	REAL		*4	MOONANG
	REAL		*4	MOONAZ
	REAL		*4	MOONPHZ
	REAL		*4	ALT
	REAL		*4	BARYVEL
	REAL		*4	LPARAM
	REAL		*4	GLAT
	REAL		*4	GLONG
	REAL		*4	ELAT
	REAL		*4	ELONG
	REAL		*4	ORBPHZ
	REAL		*4	GEOVEL
	REAL		*4	SCANANG
	REAL		*4	SCROTANG
	REAL            *4      EXCGALAT
	CHARACTER	*29	PDEF
	CHARACTER	*25	INTRO
	INTEGER		*4	INDEX
	REAL		*4	ORIGIN
	

	DICTIONARY	'FUT_ATTIT'
	RECORD /FUT_ATTIT/ ATT_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

C
C Convert the attitude to recognizable units and display.
C
	RA = ATT_REC.RA * FAC_ATT_CONV
	DEC = ATT_REC.DEC * FAC_ATT_CONV
	TLAT = ATT_REC.TERR_LATITUDE * FAC_ATT_CONV
	TLONG = ATT_REC.TERR_LONGITUDE * FAC_ATT_CONV
	ELIMB = ATT_REC.EARTH_LIMB * FAC_ATT_CONV
	ELIMBAZ = ATT_REC.EARTH_LIMB_AZIMUTH * FAC_ATT_CONV
	SUNANG = ATT_REC.SUN_ANGLE * FAC_ATT_CONV
	MOONANG = ATT_REC.MOON_ANGLE * FAC_ATT_CONV
	MOONAZ = ATT_REC.MOON_AZ_ANGLE * FAC_ATT_CONV
	MOONPHZ = ATT_REC.MOON_PHASE * FAC_ATT_CONV
	ALT = ATT_REC.ALTITUDE * 0.1
	BARYVEL = ATT_REC.PROJECTED_BARYCENTRIC_VELOCITY * 0.01
	LPARAM = ATT_REC.MCILWAIN_L_PARAM
	GLAT = ATT_REC.GALACTIC_LATITUDE * FAC_ATT_CONV
	GLONG = ATT_REC.GALACTIC_LONGITUDE * FAC_ATT_CONV
	ELAT = ATT_REC.ECLIPTIC_LATITUDE * FAC_ATT_CONV
	ELONG = ATT_REC.ECLIPTIC_LONGITUDE * FAC_ATT_CONV
	ORBPHZ = ATT_REC.ORBITAL_PHASE * FAC_ATT_CONV
	GEOVEL = ATT_REC.PROJECTED_GEOCENTRIC_VELOCITY
	SCANANG = ATT_REC.SCAN_ANGLE * FAC_ATT_CONV
	SCROTANG = ATT_REC.SC_ROTATION_ANGLE * FAC_ATT_CONV
	EXCGALAT = ATT_REC.EXC_GALACTIC_LAT * FAC_ATT_CONV

        IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 6002, IOSTAT=STATUS)
	1		FAC_ATT_SOLN(ATT_REC.SOLUTION),
	1		ATT_REC.PIXEL_NO,
	1               ATT_REC.EQUATORIAL,
	1		RA,
	1		DEC,
	1		ATT_REC.TERR_PIXEL_NO,
	1		TLAT,
	1		TLONG,
	1		ELIMB,
	1		ELIMBAZ,
	1		SUNANG,
	1		MOONANG,
	1		MOONAZ,
	1		MOONPHZ,
	1		ATT_REC.SUN_MOON_DIST,
	1		ATT_REC.COBE_MOON_DIST,
	1		ALT,
	1		BARYVEL,
	1		LPARAM,
	1		GLAT,
	1		GLONG,
	1		ELAT,
	1		ELONG,
	1		ORBPHZ,
	1		GEOVEL,
	1		SCANANG,
	1		SCROTANG,
	1               EXCGALAT

6002	FORMAT (//'***** Attitude Information *****'/
	1	' Attitude Solution: ', T50, A /
	1	' Celestial Pixel Number: ', T50, I4 /
	1    
	1	' FIRAS Pointing (Celestial Equatorial Coordinates):   X=',
	1	     G13.6, ', Y=', G13.6, ', Z=', G13.6 /
	1	' Right Ascension (deg.), Declination (deg.): ', T50, G13.6,
	1	     ', ', G13.6 /
	1	' Terrestrial Pixel Number: ', T50, I4 /
	1	' Terrestrial Latitude, Longitude (degrees): ', T50, G13.6,
	1	     ', ', G13.6 /
	1	' Earth Limb Angle (degrees): ', T50, G13.6 /
	1	' Earth Limb Azimuth (degrees): ', T50, G13.6 /
	1	' Sun Angle (degrees): ', T50, G13.6 /
	1	' Moon Angle (degrees): ', T50, G13.6 /
	1	' Moon Azimuth (degrees): ', T50, G13.6 /
	1	' Moon Phase (degrees): ', T50, G13.6 /
	1	' Sun-Moon Distance (km): ', T50, G13.6 /
	1	' COBE-Moon Distance (km): ', T50, G13.6 /
	1	' COBE Altitude (km): ', T50, G13.6 /
	1	' Projected Barycentric Velocity (km/sec): ', T50, G13.6 /
	1	' McIlwain L Parameter (Earth Radii): ', T50, G13.6 /
	1	' Galactic Latitude, Longitude (degrees): ', T50, G13.6,
	1	     ', ', G13.6 /
	1	' Ecliptic Latitude, Longitude (degrees): ', T50, G13.6,
	1	     ', ', G13.6 /
	1	' Orbital Phase (degrees): ', T50, G13.6 /
	1	' Projected Geocentric Velocity (m/sec): ', T50, G13.6 /
	1	' Scan Angle (degrees): ', T50, G13.6 /
	1	' Spacecraft Rotation Angle (degrees): ', T50, G13.6 /
	1       ' Exclude data within Galactic Latitude: ', T50, G13.6)

        IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	IF (ATT_REC.PIXEL_DEFINITION .Eq. 'Q') Then
	  PDEF = 'Standard Quad-Cube Pixel:'
	  INTRO = 'Skymap Index:'
	  INDEX = ATT_REC.SKYMAP_INDEX
	ELSE IF (ATT_REC.PIXEL_DEFINITION .Eq. 'O') Then
	  PDEF = 'Orbit Average Pixel:'
	  INTRO = 'Orbit Origin (degrees):'
	  ORIGIN = 0.1 * ATT_REC.SKYMAP_INDEX
	ELSE IF (ATT_REC.PIXEL_DEFINITION .Eq. 'S') Then
	  PDEF = 'Scan Angle Pixel:'
	  INTRO = 'Number of Scan Pixels:'
	  INDEX = ATT_REC.SKYMAP_INDEX
	ELSE IF (ATT_REC.PIXEL_DEFINITION .Eq. 'E') Then
	  PDEF = 'Geocentric Quad-Cube Pixel:'
	  INTRO = 'Skymap Index:'
	  INDEX = ATT_REC.SKYMAP_INDEX
	ELSE
	  PDEF = ATT_REC.PIXEL_DEFINITION
	  INTRO = 'Skymap Index:'
	  INDEX = ATT_REC.SKYMAP_INDEX
	End IF

	IF (ATT_REC.PIXEL_DEFINITION .Eq. 'O') Then
	  WRITE (WRITE_LUN,6010,IOSTAT=STATUS)
	2          PDEF,INTRO,ORIGIN
	ELSE
	  WRITE (WRITE_LUN,6020,IOSTAT=STATUS)
	2          PDEF,INTRO,INDEX
	End IF

6010	FORMAT (/' Alternate Pixel Information: '/
	1	' ', T10, A30, T50, A26, F6.1)

        IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

6020	FORMAT (/' Alternate Pixel Information: '/
	1	' ', T10, A30, T50, A26, I4)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C    Terrestrial radiation region location flag

	WRITE (WRITE_LUN,6040,IOSTAT=STATUS) ATT_REC.TERR_RAD_BYTE
6040    FORMAT (10x,/' Terrestrial Radiation Bit Flag: ',I4)

        IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
