	integer*4 function fut_default_peak(mtmspeed, scanlength,
     *	  chan, ngroup, sci_mode, linearized, loc_peak)
c
c  Given the MTM speed and data compression, this function returns the
c  index number of the position of the zero path difference in the
c  interferogram.  
c
c  Input parameters:	mtmspeed	integer*4	0=slow,1=fast
c			scanlength	integer*4	0=short,1=long
c			chan		integer*4	1, 2, 3, or 4
c			sci_mode	integer*4	0, 1, 2, 3, or 4
c			ngroup		integer*4	adds per group
c
c...The following input parameter has been activated to work correctly for the
c.....most common scan modes
c			linearized	integer*4	0=nonlinearized ifg,
c							 1=linearized ifg
c				
c  Output parameter:	loc_peak	integer*4	peak index number
c
c  Return status:	 fut_normal = normal
c			 fut_invdefpeak = peak is outside range (1,512)
c
c  Written by:		Rich Isaacman
c			Applied Research Corp.
c			29 Dec 1988
c			286-5758
c
c  Revision:		Alice Trenholme
c			General Sciences Corp.
c			9 February 1990
c			286-5569
c
c	SER 6794. Increased variables used in finding peak by including channel, 
c	microprocessor mode, scan length, and linearization factor in input 
c	variable list.  Also employs a "fudge factor" based on empirical data to
c	correctly give the peak position in all commonly used configurations of
c	scan modes and adds per group.  Complete information follows:
c
c
c                    INTERFEROGRAM PEAK POSITION
c
c	Calculation of the peak position (in counts) of the interferogram has 
c  been enhanced by several factors.  These include the following:
c
c	1.  The raw peak location is initially calculated to be 2 counts plus
c  the distance to the peak, 1.2 cm, divided by adds per group times the
c  distance covered for each sample, which is .001152 cm or .001728 cm
c  respectively for slow or fast scan. 
c
c	2.  Information is delayed by a buffer associated with the onboard
c  deglitcher; the number of counts "behind" that the information is depends on
c  scan length, mtm speed, and adds per group; in the common configurations it
c  is 3  for the short slow scan mode, 2 for the short fast scan mode, and 8 for
c  long  scan modes.  This buffer factor is added to the peak position. 
c  
c  	3.  If the digital filters are on, the information is delayed by a 
c  factor which has been empirically determined by examination of on-orbit data.
c  This factor is divided by adds per group times 2/3 for slow speed, 1 for
c  fast.  The integer quotient of 4 divided by adds per group is added to this
c  factor, and the whole is then added to the peak location. 
c  
c  	4.  If the interferogram has been linearized, the peak location is
c  shifted to the left by a factor which again has been empirically determined
c  by examination of on-orbit and I_&T data.  The factor for the common
c  instrument states is 2 for long scans and for short slow scan in the high
c  channel.  It is 3 for short fast scans in the high channel, 5 for short slow
c  scans in the low channels, and 6 for short fast scans in the low channels. 
c  
c
c            EMPIRICALLY DETERMINED FACTORS USED IN FUT_DEFAULT_PEAK
c                  Common scan modes and adds per group
c
c    =======================================================================
c    |       Common Configurations       | Microprocessor |  Linearization |
c    | Scan Mode  | Channel | Adds/Group |   Shift        |    Shift       |
c    |============|=========|============|================|================|
c    | Short Slow |   Low   |     3      |   + 10         |      - 5       |
c    |------------|---------|------------|----------------|----------------|
c    | Short Slow |   High  |     3      |    + 7         |      - 2       |
c    |------------|---------|------------|----------------|----------------|
c    | Long Fast  |   Low   |     8      |    + 1         |      - 2       |
c    |------------|---------|------------|----------------|----------------|
c    | Long Fast  |   High  |     2      |    + 2         |      - 2       |
c    |------------|---------|------------|----------------|----------------|
c    | Short Fast |   Low   |     2      |   + 12         |      - 6       |
c    |------------|---------|------------|----------------|----------------|
c    | Short Fast |   High  |     2      |    + 9         |      - 3       |
c    |------------|---------|------------|----------------|----------------|
c    | Long Slow  |   Low   |    12      |    + 2         |      - 2       |
c    |------------|---------|------------|----------------|----------------|
c    | Long Slow  |   High  |     3      |    + 4         |      - 2       |
c    =======================================================================


	implicit none
	real*8		dx(2)			!point spacing in cm (slow/fast)
	real*8		offset			!ZPD offset from scan start (cm)
	integer*4	mtmspeed		! mtm speed
	integer*4	ngroup			! adds per group
	integer*4	loc_peak		! peak location
	integer*4	scanlength		! scan length
	integer*4	chan			! channel number 
	integer*4	ch			! 0=low, 1=high channels 
	integer*4	sci_mode		! science mode
	integer*4	sc_m_or_l	! 0=sci mode 1,3 and not linearized;
c					       1=sci mode 2,4 or linearized
	integer*4	buffactor		! buffer factor
	integer*4	data_indx		! index for fudge factor
	integer*4	linearized		! 0=non, 1=linearized
	integer*4	takeoff		       !counts to subtract if linearized

	real*8		shift			! shift due to microprocessors 
	real*8		addon(8)		!fudge factor based on empirical
						! data.  
	external	fut_normal
	external	fut_invdefpeak

	data 		dx /1.152e-03, 1.728e-03/

	data		offset/1.2000/
	data		addon/18.,12.,20.,14.,16.,6.,8.,0./
     *		

c...Set channel indicator to 1 for high channels and 0 for low

	if ((chan .eq. 1) .or. (chan .eq. 3)) then
		ch = 1
	else
		ch = 0
	endif


c...Set science mode/linearization indicator to 1 if digital filters are on OR 
c.....if linearization has occurred, 0 if not.  Since the peak position after
c.....linearization is not affected by whether the digital filters are on or 
c.....off, the shift factor is added in regardless when the interferogram is
c.....already linearized. 
 
	if (((sci_mode .eq. 1) .or. (sci_mode .eq. 3)) .and. 
     .            (.not.(linearized))) then 
		sc_m_or_l = 0
	else
		sc_m_or_l = 1
	endif


c...Index of fudge factor is determined by binary arithmetic

	data_indx = 1 + ch + 2*mtmspeed + 4*scanlength 


c...Initial peak location is 2 + distance to peak (1.2 cm.) divided by adds per 
c.....group times distance per sample

	loc_peak = offset/(ngroup*dx(mtmspeed+1)) + 2.


c...If digital filters are on or if thew interferogram has been linearized,
c.....shift is the "fudge factor" divided by adds per group times 2/3 for slow
c.....speed, 1 for fast speed, plus the integer c.....quotient of 4 divided by
c.....adds per group.  Shift is 0 if digital filters are off and the
c.....interferogram is not linearized. 

	shift =  sc_m_or_l*(3*addon(data_indx)/(ngroup*(2.+mtmspeed))) + 
     *		sc_m_or_l*(4/ngroup)

c...Buffer factor is 3 for short slow scan mode, 2 for short fast scan mode, 
c.....and 8 for long scan modes

	buffactor = (scanlength*(5 + mtmspeed)-mtmspeed+3)/ngroup

c...Amount to take off if linearized is 2 for long scan length or short slow, 
c.....high channel.  It is 3 for short fast, high channel, 5 for short slow, low
c.....channel, and 6 for short fast, low channel.

	takeoff = 
     *   linearized*(2*scanlength + (1-scanlength)*(mtmspeed + 5 - 3*ch))

c...peak location is the sum of the original plus the shift and buffer factor,
c.....less the amount subtracted after linearization

	loc_peak =  loc_peak + shift + buffactor - takeoff

c...peak is good if between 1 and 512; otherwise not good

	if (loc_peak.ge.1 .and. loc_peak.le.512) then
	   fut_default_peak = %loc(fut_normal)
	else
	   call lib$signal(fut_invdefpeak)
	   fut_default_peak = %loc(fut_invdefpeak)
	endif

	return
	end
