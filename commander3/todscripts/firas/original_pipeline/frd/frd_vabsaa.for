	Program FRD_VABSAA
C------------------------------------------------------------------------------
C  PURPOSE:
C    Read FREF files FEX_VABN, FEX_VABS, and FEX_SAA.TXT, and rearrange them to
C    form and write FEX_VABSAA.DAT.  The 3 *.TXT files were made by rearranging
C    the 5 files CSDR$MPD:VABxn.DAT (x=S,N; n=1,2) and CSDR$MPD:SAA.DAT into
C    (longitude, south lat, north lat)-triplets describing the south and north
C    boundaries of one of the radiation belts.  VABN defines the north Van Allen
C    Belt; VABS, the south Van Allen Belt; SAA, the South Atlantic Anomaly.
C
C  CALLING SEQUENCE:
C    This is a main program.
C
C  AUTHOR:  Fred Shuman,  STX,  1991 Apr 23
C
C  INPUT FILES:
C    FEX_VABN, FEX_VABS, and FEX_SAA.TXT
C
C  OUTPUT FILE:
C    FEX_VABSAA.DAT
C
C  INCLUDE FILES USED: None
C
C  SUBROUTINES AND FUNCTIONS CALLED: None
C
C------------------------------------------------------------------------------

	Implicit	None

	Dictionary	'FEX_VABSAA'
	Record / FEX_VABSAA / Rad

	Logical		* 4	next
	Character	*14	gmt
	Character	*32	text, infile, outfile / 'FEX_VABSAA.DAT' /
	Character	*32	filename(2) / 'FEX_VABN.TXT', 'FEX_VABS.TXT' /
	Integer		* 4	j, k
	Real		* 4	lon, lats, latn
	Real		* 4	lonprev, latsprev, latnprev
	Real		* 4	step
	Integer		* 4	ios, status, lun, curradt(2)

	External	FRD_RadBadLon, FRD_RMSOpen, FRD_RMSWrite

	Call Lib$Get_Lun(lun)
C
C  Read the 3 'Earth radiation field' files.  First, the Van Allen Belt files:
C
	Call SYS$GetTim(curradt)
	Call CT_Binary_to_GMT(curradt,gmt)

	Rad.ct_head.gmt = gmt

	Rad.ct_head.time(1) = curradt(1)
	Rad.ct_head.time(2) = curradt(2)

	Do j=1,2
	   infile = filename(j)

	   Open (UNIT=lun,    FILE=infile,     FORM='FORMATTED',
	2        READONLY,    STATUS='OLD',    RECORDTYPE='VARIABLE',
	3        SHARED,      IOSTAT=ios)

	   If (ios .Ne. 0) Then
	      Call Lib$Signal(FRD_RMSOpen, %Val(2), infile, %Val(ios))
	   End If

C   First, read past the comments in the text file...

	   ios = 1
	   Do While (ios .Ne. 0)
	      Read (lun,*,IOSTAT=ios) lonprev, latsprev, latnprev
	   EndDo

C   Now read the real data...

	   Read (lun,*,IOSTAT=ios) lon, lats, latn
	   k = 1
	   Rad.vab(j).latmin  = Min(latsprev, lats)
	   Rad.vab(j).latmax  = Max(latnprev, latn)
	   Rad.vab(j).lonstep = lon - lonprev
	   Rad.vab(j).lats(1) = latsprev
	   Rad.vab(j).latn(1) = latnprev

	   Do While (ios .Eq. 0)
	      If (lon .Eq. lonprev+Rad.vab(j).lonstep) Then
	         k = k + 1
	         Rad.vab(j).lats(k) = lats
	         Rad.vab(j).latn(k) = latn
	      Else
	         Call LIB$Signal( FRD_RadBadLon, %Val(4), infile,
	2                  %Val(Rad.vab(j).lonstep), %Val(lonprev), %Val(lon) )
	      End If

	      Rad.vab(j).latmin = Min(Rad.vab(j).latmin, lats)
	      Rad.vab(j).latmax = Max(Rad.vab(j).latmax, latn)
	      lonprev = lon
	      latsprev = lats
	      latnprev = latn
	      Read (lun,*,IOSTAT=ios) lon, lats, latn
	   End Do

	   Close(UNIT=lun)
	End Do
C
C  Now read the S. Atlantic Anomaly:
C
	infile = 'FEX_SAA.TXT'

	Open (UNIT=lun,    FILE=infile,     FORM='FORMATTED',
	2     READONLY,    STATUS='OLD',    RECORDTYPE='VARIABLE',
	3     SHARED,      IOSTAT=ios)

	If (ios .Ne. 0) Then
	   Call Lib$Signal(FRD_RMSOpen, %Val(2), infile, %Val(ios))
	End If
C
C   The SAA file is a sequence of (long., S.lat., N.lat.) records defining the
C   South and North boundaries of the SAA loop.  It starts at the min long. of
C   the loop and goes forward in equal longitude steps along both boundaries,
C   ending at the max long.  We store this into 2 arrays, Rad.saa.latn=North &
C   lats=South boundary.  First, read past the comments in the text file...

	ios = 1
	Do While (ios .Ne. 0)
	   Read (lun,*,IOSTAT=ios) lonprev, latsprev, latnprev
	EndDo

C  Now read the real data...

	Read (lun,*,IOSTAT=ios) lon, lats, latn
	k = 1
	Rad.saa.lonmin  = lonprev
	Rad.saa.latmin  = Min(latsprev, lats)
	Rad.saa.latmax  = Max(latnprev, latn)
	Rad.saa.lonstep = lon - lonprev
	Rad.saa.lats(1) = latsprev
	Rad.saa.latn(1) = latnprev

	Do While (ios .Eq. 0)
	   If (lon .Eq. lonprev+Rad.saa.lonstep) Then
	      k = k + 1
	      Rad.saa.lats(k) = lats
	      Rad.saa.latn(k) = latn
	   Else
	      Call LIB$Signal( FRD_RadBadLon, %Val(4), infile,
	2                  %Val(Rad.vab(j).lonstep), %Val(lonprev), %Val(lon) )
	   End If

	   Rad.saa.latmin = Min(Rad.saa.latmin, lats)
	   Rad.saa.latmax = Max(Rad.saa.latmax, latn)
	   lonprev = lon
	   latsprev = lats
	   latnprev = latn
	   Read (lun,*,IOSTAT=ios) lon, lats, latn
	End Do
	Rad.saa.lonmax = lon

	Close(UNIT=lun)
C
C Finished reading the 3 'Earth rad field' files.  Write the output record.
C
	Open (UNIT=lun,           NAME=outfile, STATUS='new',
	2     FORM='unformatted', ACCESS='sequential',
	3     RECORDSIZE=256,     ORGANIZATION='sequential',
	4     RECORDTYPE='fixed', IOSTAT=ios,   SHARED)

	If (ios .Eq. 0) Then

	   Write (lun, IOSTAT=ios) Rad

	   If (ios .Ne. 0) Then
	      status = %Loc(FRD_RMSWrite)
	      Call LIB$Signal(FRD_RMSWrite, %Val(2), outfile, %Val(ios))
	   End If

	Else
	   status = %Loc(FRD_RMSOpen)
	   Call LIB$Signal(FRD_RMSOpen, %Val(2), outfile, %Val(ios))
	End If

	Close (lun, IOSTAT=ios)

	End
