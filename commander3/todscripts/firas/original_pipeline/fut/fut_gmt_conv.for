	Integer*4 Function FUT_GMT_CONV (time_sec,gmt,flag)

C  Larry Rosen, STX, March 7 1991.
C  If flag = 1, convert time in seconds to GMT.
C  If flag = 2, convert gmt to time in seconds.

	Real*8		time_sec
	Character*14	gmt
	integer*2	flag
	Real*8		thold, tyr, tday, thr, tmin, tsec

	Integer*2	iyr, iday, ihr, imin
	Integer*4	msec
	Include		'(FUT_Params)'

	IF (flag .eq. 1) Then
	    iyr = 89 + int (time_sec/(365.*fac_day))
	    thold = time_sec - (iyr-89) * (365.*fac_day)
	    iday = int (thold / fac_day) + 1
   	    thold = thold - dble((iday-1)*fac_day)
	    ihr = int (thold / fac_hour)
	    thold = thold - dble(ihr*fac_hour)
	    imin = int (thold / fac_minute)
	    thold = thold - dble(imin*fac_minute)
	    msec = int (thold*1000)
	    write (gmt,10) iyr,iday,ihr,imin,msec
  10	    FORMAT (I2,I3.3,I2.2,I2.2,I5.5)
	Else
	    READ (Gmt,10) iyr,iday,ihr,imin,msec
	    tyr = dble((iyr-89)*365.)*fac_day
	    tday = dble(iday-1) * fac_day
	    thr = dble(ihr) * fac_hour
	    tmin = dble(imin) * fac_minute
	    tsec = dble(msec*0.001)
	    time_sec = tsec + tmin + thr + tday + tyr
	EndIf
	Return
	End
