	Integer*4 Function FUT_VABSAAFlag(lon, lat, flag)

C------------------------------------------------------------------------------
C    PURPOSE:  Given terrestrial longitude and latitude in degrees, determine
C	whether it is located in the SAA, one of the VABelts, or none of these.
C
C    AUTHOR: Fred Shuman,  ST Systems Corporation,  1991 April 17
C
C    INVOCATION: status = FUT_VABSAAFlag(lon, lat)
C
C    INPUT ARGUMENTS:
C	LON	R*4	Terrestrial longitude in degrees
C	LAT	R*4	Terrestrial latitude in degrees
C
C    OUTPUT ARGUMENT:
C	FLAG		B	Flag stating whether the given terrestrial
C				position is within a radiation belt:
C					1 = Outside the radiation belts
C		(These values are	2 = Inside the North Van Allen Belt
C		 as specified in	4 = Inside the South Van Allen Belt
C		 ATTITUDE.RDL)		8 = Inside the South Atlantic Anomaly
C
C    SUBROUTINES CALLED:
C	CCT_Open_Config
C	CCT_Get_Config_TOD
C	CCT_Close_Config
C
C    INCLUDE FILES:
C	CCT_Get_Config
C	$SSDef
C
C    PROCESSING METHOD:
C	IF this is the first call to this routine THEN
C	  OPEN the reference dataset CSDR$FIRAS_Ref:FEX_VABSAA
C	  READ its contents into structure RAD (Dictionary element FEX_VABSAA)
C	ENDIF
C	SET FLAG to indicate that the given terr position is outside the
C	      radiation regions
C	DO for both Van Allen Belts
C	  IF the given latitude is between the upper and lower limits given in
C	        the reference file THEN (go ahead and check the detailed VAB...)
C	    COMPUTE the consecutive longitudes in the reference file that
C	          bracket the given longitude.
C	    INTERPOLATE the latitude at the upper and lower boundaries of VAB
C	    IF the given latitude is between the upper and lower VAB bounds THEN
C	      SET FLAG to indicate that the given terr pos is inside this VAB
C	    ENDIF
C	  ENDIF
C	ENDDO
C	IF the given long and lat are both within upper and lower limits given
C	      for the SAA in the reference file THEN (check the detailed SAA...)
C	  COMPUTE the consecutive longitudes in the reference file that
C	        bracket the given longitude.
C	  INTERPOLATE the latitude at the upper and lower boundaries of the SAA
C	  IF the given latitude is between the upper and lower SAA bounds THEN
C	    SET FLAG to indicate that the given terr pos is inside the SAA
C	  ENDIF
C	ENDIF
C------------------------------------------------------------------------------

	Implicit None

C  Include files

	Include  '(CCT_Get_Config)'
	Include  '($SSDef)'

	Record / Config_status / stat

	Dictionary 'FEX_VABSAA'
	Record / FEX_VABSAA / Rad

C  Calling arguments

	Real		*4	lon, lat
	Byte			flag

C  Other variables

	Logical		*4	first / .True. /
	Integer		*4	i, j, lun, status
	Real		*4	lon1, lon2, lats1, lats2, latn1, latn2
	Real		*4	step, prod, lonmin / -180. /

C  Excess baggage for CCT_*_Config*

	Character	*14	gmtstart / '86001000000000' /
	Character	*14	gmtstop / '99365000000000' /
	Character	*14	gmt / '89365235959999' /
	Integer		*4	ndset / 1 /, adtstart(2), adtstop(2), adt(2)
	Character	*25	dataset / 'CSDR$FIRAS_Ref:FEX_VABSAA' /
	Integer		*4	size / 1024 /
	Character	*1	accmode / ' ' /
	Integer		*4	ncache / 1 /, index, refct
	Logical		*1	newseg

C  Functions and Externals

	Integer		*4	CCT_Open_Config
	Integer		*4	CCT_Get_Config_TOD
	Integer		*4	CCT_Close_Config
	Integer		*4	Lib$Get_Lun, Lib$Free_Lun

	External		FUT_Normal

	Save first
C Set return status.

	FUT_VABSAAFlag= %loc(FUT_Normal)

C Get a unit number and open the reference arcv only on the first pass.

	If (first) Then

	   first = .False.
	   status = Lib$Get_Lun(lun)
	   If ( status .Ne. SS$_Normal ) Then
	      FUT_VABSAAFlag = status
	      Return
	   End If

C Open the radiation field reference file for read.

	   Call CT_GMT_to_Binary(gmtstart,adtstart)
	   Call CT_GMT_to_Binary(gmtstop, adtstop)
	   Call CT_GMT_to_Binary(gmt, adt)
	   status = CCT_Open_Config (adtstart, adtstop, ndset, dataset, size,
	2                            accmode, ncache, lun, index, stat, refct)
	   If (.Not. status) Then
	      FUT_VABSAAFlag = status
	   Else
	      status = CCT_Get_Config_TOD(adt, ndset, size, lun, index,
	2                                 rad, newseg, stat)
	      status = CCT_Close_Config ( ndset, lun, index )
	   End If
	End If
C
C Unless we find otherwise, the given position is outside the radiation belts.
C
	flag = 1
C
C Is the terrestrial position within a VAB?  Check for inclusion in each VAB iff
C   the latitude is between the minimum and maximum latitudes of that VAB.
C
	Do i=1,2
	   If (lat .Ge. Rad.VAB(i).latmin .And. lat .Le. Rad.VAB(i).latmax) Then
	      step = Rad.VAB(i).lonstep
	      j = (lon - lonmin)/step
	      lon1 = lonmin + j*step
	      lon2 = lon1 + step
	      latn1 = Rad.VAB(i).latn(j+1)
	      latn2 = Rad.VAB(i).latn(j+2)
	      lats1 = Rad.VAB(i).lats(j+1)
	      lats2 = Rad.VAB(i).lats(j+2)
	      prod = step*lat
	      If (prod .Ge. (lon2-lon)*lats1 + (lon-lon1)*lats2  .And.
	2         prod .Le. (lon2-lon)*latn1 + (lon-lon1)*latn2)    Then
	         flag = 2*i
	      End If
	   End If
	End Do
C
C Check for inclusion in the SAA iff the terrestrial position falls within the
C   longitude-latitude "box" that circumscribes the SAA region.
C   If it does, check whether the terrestrial latitude is between the "north"
C   and "south" latitudes interpolated from the SAA boundary tables.
C
	If (lon .Ge. Rad.SAA.lonmin .And. lon .Le. Rad.SAA.lonmax .And.
	2   lat .Ge. Rad.SAA.latmin .And. lat .Le. Rad.SAA.latmax )      Then
	   step = Rad.VAB(i).lonstep
	   j = (lon - Rad.SAA.lonmin)/step
	   lon1 = Rad.SAA.lonmin + j*step
	   lon2 = lon1 + step
	   latn1 = Rad.SAA.latn(j+1)
	   latn2 = Rad.SAA.latn(j+2)
	   lats1 = Rad.SAA.lats(j+1)
	   lats2 = Rad.SAA.lats(j+2)
	   prod = step*lat
	   If (prod .Ge. (lon2-lon)*lats1 + (lon-lon1)*lats2  .And.
	2      prod .Le. (lon2-lon)*latn1 + (lon-lon1)*latn2)    Then
	      flag = 8
	   End If
	End If

	Return
	End
