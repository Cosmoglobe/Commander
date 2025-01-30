	subroutine fep_dplotter (nrec, temps, dwell_loc, timetags, jside,
	2                        ngrt, grtnums, icurrent)

C-----------------------------------------------------------------------------
C
C Changes:
C
C	SPR 3132, Dwellplot fails to plot XCAL S5 and S6. R. Kummerer, STX,
C		January 12, 1989.
C
C	SPR 2922, Spurious and misplaced values.  Fred Shuman, STX,
C		1989 Feb 10.
C
C	SPR 2910, Show all dwell data in the timerange.  Fred Shuman, STX,
C		1989 Feb 17.
C
C	SPR 2738, Add plotting of cal resistors (in counts).  Fred Shuman, STX,
C		1989 Feb 20.
C
C-----------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Params)'

	Character *32 fcc_plot_device
	Character *14 timetags(1)
	Character *80 title
	Character *40 xlabel, ylabel(9)
	Character *15 names(16), buffname
	Character *3  side(2)
	Character *4  curr(2)
	Real	  *4  temps(1), timebuff(1000)
	Real	  *8  hr
	Integer	  *2  irec, nrec, dwell_loc(1), oldloc
	Integer	  *4  i, j, k, m, n
	Integer	  *4  iyr, iday, ihr, imin, msec
	Integer	  *4  jside, ngrt, grtnums(1), igrt
	Integer	  *4  icurrent, dlen, len, ioffset, minor
	Logical	  *1  same_units

	Integer   *4	nplotpts, max_pts
	Parameter (max_pts = 64000)
	Real	  *4	maxtime
	Real	  *4	timeplotmin
	Real	  *4	timeplotmax
	Real	  *4	ymin
	Real	  *4	ymax
	Real	  *4	yplotmin
	Real	  *4	yplotmax
	Real	  *4	ymn(5)
	Real	  *4	ymx(5)
	Real	  *4	yplotmn(5)
	Real	  *4	yplotmx(5)
	Real	  *4	y(max_pts, 10)
	Integer   *4	iery(6)
	Integer   *4	nveci
	Character *150	cmd(30)
	Integer   *4	ncmd
	Integer   *4	nfields
	Integer   *4	ier

	data fcc_plot_device /'/VT'/
	data names /'XCAL#',      'Sky Horn#',  'Ref Horn#',  'ICAL#',
	2           'Dihedral#',  'Detec RH#',  'Detec RL#',  'Detec LH#',
	3           'Detec LL#',  'Mirror M1#', 'Cal Res RH#', 'Cal Res RL#',
	4           'Cal Res LH#', 'Cal Res LL#', 'XCAL Sn#',   'Collim C1#'/
	data side /' A#', ' B#'/
	data curr/' Lo#', ' Hi#'/
c
c  First create plot label by looking at dwell address and side flag (A or B).
c  Use these to stitch together the appropriate character string, then append
c  the time tag.
c
	xlabel = 'Time (Minutes after #'
	len = index (xlabel,'#') - 1
	xlabel(len+1:len+13) = timetags(1)(1:11) // ')'

	title = 'Dwell:'
	len = 6

	same_units = .True.
	Do i=1,ngrt
	   igrt = grtnums(i) - 16*jside
	   dlen = index (names(igrt),'#') - 1
	   title(len+1 : len+dlen+7) = ', ' // names(igrt)(1:dlen)
	2             // side(jside+1)(1:2) // curr(icurrent+1)(1:3)
	   len = len + dlen + 7
	
	   If (igrt .Lt. 11  .Or.  igrt .Gt. 14) Then
	      ylabel(i) = 'Degrees K'
	      If (igrt .Eq. 15) Then
	         If(jside .Eq. 0) Then
	            title(len-5:len-5) = '5'
	         Else
	            title(len-5:len-5) = '6'
	         End If
	      End If
	   Else
	      ylabel(i) = 'Counts'
	   End If

	   If (ylabel(i) .Ne. ylabel(1)) Then
	      same_units = .False.
	   End If
	End Do

	title(7:7) = ' '
c
c Convert the ASCII time tags into seconds measured from the start of 1986.
c  Shortly we will put that value into every 64th record of the time axis buffer
c     (since there are 64 GRT readings in every housekeeping record).
c  Note that the last 5 digits can be taken together as milliseconds.
c  Also, any shift of the "86" must be by multiples of 4.
c
	do j=1,nrec
	   Read (timetags(j),20) iyr, iday, ihr, imin, msec
20	   Format (i2,i3,i2,i2,i5)
	   iyr = iyr - 86
	   timebuff(j) = (INT(iyr*365.25d0 +.4d0) + iday)*fac_day
	2           + ihr*fac_hour + imin*fac_minute + msec*1.d-3
	enddo
c
c  Refer all times to that of the first record.
c
	do n=nrec,1,-1
	   timebuff(n) = timebuff(n) - timebuff(1)
	enddo
c
c  Fill in all the intermediate times between housekeeping record boundaries.
c  HSKP records are 64 seconds apart, so GRT readings are at 1 sec intervals.
c  Load times into y(*,1).
c
	do j=1,nrec
	   ioffset = 64 * (j-1) + 1
	   y(ioffset,1) = timebuff(j)
	   do minor=1,63
	      m = ioffset + minor
	      y(m,1) = y(m-1,1) + 1.
	   enddo
	enddo
c
c  Convert times y(*,1) from sec to min; load temps into y(*,2:nfields+1)
c
	nplotpts = 64*nrec
	igrt = 1
	oldloc = dwell_loc(1)
	Do irec=1,nrec
	   ioffset = 64 * (irec-1) + 1
	   Do minor=0,63
	      m = ioffset + minor
	      y(m,1) = y(m,1)/60.
	      Do k=1,ngrt
	         y(m,k+1) = fac_no_data
	      End Do
	      If (dwell_loc(irec) .Ne. oldloc) Then
	         igrt = igrt + 1
	         oldloc = dwell_loc(irec)
	      End If
	      If (temps(m) .Ge. 1) Then
	         y(m,igrt+1) = temps(m)
	      End If
	   End Do
	End Do
c
c Make the plot...
c
c     Find the left and right plot boundaries.
c
	maxtime = y(nplotpts,1)
	timeplotmin = - (fac_margin * maxtime)
	timeplotmax = (1. + fac_margin)* maxtime
c
c Set up the PLT calling arguments, iery(), nveci, cmd(), and ncmd and find the
c     top and bottom plot boundaries.  Skip over "no_data" to find the first
c     good value.
c
	j = 1
	Do While (y(j,2) .Eq. fac_no_data)
	   j = j + 1
	End Do
	ymin = y(j,2)
	ymax = y(j,2)

	nfields = ngrt
	Do k=1,nfields
	   ymn(k) = y(1,k+1)
	   ymx(k) = y(1,k+1)
	   Do j=1,nplotpts
	      If (y(j,k+1) .Ne. fac_no_data) Then
	         ymn(k) = amin1(ymn(k), y(j,k+1))
	         ymx(k) = amax1(ymx(k), y(j,k+1))
	      End If 
	   End Do
	   ymin = amin1(ymin, ymn(k))
	   ymax = amax1(ymax, ymx(k))

	   If (ymn(k) .Eq. ymx(k)) Then
	      If (ymn(k) .Gt. 0.) Then
	         yplotmn(k) = (1. - fac_margin)*ymn(k)
	         yplotmx(k) = (1. + fac_margin)*ymx(k)
	      Else if (ymn(k) .Eq. 0.) Then
	         yplotmn(k) = -1.
	         yplotmx(k) =  1.
	      Else
	         yplotmn(k) = (1. + fac_margin)*ymn(k)
	         yplotmx(k) = (1. - fac_margin)*ymx(k)
	      End If 
	   Else
	      yplotmn(k) = ymn(k) - fac_margin*(ymx(k) - ymn(k))
	      yplotmx(k) = ymx(k) + fac_margin*(ymx(k) - ymn(k))
	   End If 
	End Do

	If (ymin .Eq. ymax) Then
	   If (ymin .Gt. 0.) Then
	      yplotmin = (1. - fac_margin)*ymin
	      yplotmax = (1. + fac_margin)*ymax
	   Else if (ymin .Eq. 0.) Then
	      yplotmin = -1.
	      yplotmax =  1.
	   Else
	      yplotmin = (1. + fac_margin)*ymin
	      yplotmax = (1. - fac_margin)*ymax
	   End If 
	Else
	   yplotmin = ymin - fac_margin*(ymax - ymin)
	   yplotmax = ymax + fac_margin*(ymax - ymin)
	End If 

	If (nplotpts .Eq. 1) Then
	   nplotpts = 2
	   y(2,1) =  maxtime
	   Do k=1,nfields
	      y(2,k+1) = -1.
	   End Do
	End If

	nveci = nfields + 1
	Do k=1,nveci
	   iery(k) = 0.
	End Do
	cmd(1) = 'D ' // fcc_plot_device
	cmd(2) = 'V .2 .1 .9 .8'
	cmd(3) = 'CS 1.5'
	cmd(4) = 'LA OT '
	cmd(5) = 'LA T '  // title
	cmd(6) = 'LA X '  // xlabel
	cmd(7) = 'LA Y '  // ylabel(1)
	ncmd = 7

	Do k=1,ngrt
	   Write ( cmd(ncmd+k), '(a4, i1, a41)' )
	2         'LA Y', k+1, ' ' // ylabel(k)
	End Do

	ncmd = ncmd + ngrt + 2
	cmd(ncmd-1) = 'LA F'
	cmd(ncmd) = 'CO 1 ON 2'

	ncmd = ncmd + 1
	Write (cmd(ncmd), '(a, 4(1x, e14.4))')
	2      'R', timeplotmin, timeplotmax, yplotmin, yplotmax

	ncmd = ncmd + 1
	If (same_units) Then
	   cmd(ncmd) = 'PLOT OVERLAY'
	Else
	   cmd(ncmd) = 'PLOT VERTICAL'
	End If

	Call PLT (y, iery, max_pts, nplotpts, nveci, cmd, ncmd, ier)

	return
	end
