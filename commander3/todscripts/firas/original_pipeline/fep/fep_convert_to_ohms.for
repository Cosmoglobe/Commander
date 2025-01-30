	subroutine fep_convert_to_ohms (oldbuff, nrec, ohms,Config_lun, 
	1	                        con_lun,dwell_loc,
	2                               timetags, jside, rcode)
C-----------------------------------------------------------------------------
C    Converts a single 2-MjF (64 pts) record of counts to ohms for dwellplot.
C
C    Input Arguments:
C	oldbuff(*)	 B	One entire housekeeping record
C	nrec		 I* 2	HKP rec # within the current timerange
C	dwell_loc	 I* 2	Dwell address, telling which GRT is dwelled
C       Config_lun       I*4    Unit number of reference file(Fdb_Grtcal) 
C       Con_Lun          I*4    Unit number of ref file (Fex_Calres)  
C	jside		 I* 4	Spacecraft side:  0=A, 1=B
C
C    Output Arguments:
C	ohms(*)		 R* 4	Great big array (64000) of converted values
C				in ohms; we fill 64 of these each time this
C				routine is called.
C	timetags(*)	Ch*14	Time tags
C       Rcode           I*4     Error return code
C                               0 = no error
C                               1 = error
C-----------------------------------------------------------------------------
C Changes:
C	Last two values appearing in each Major Frame should be placed at the
C	beginning.  (SPR 2922)  Fred Shuman, STX, 1989 Feb 16.
C
C	Add the ability to plot cal resistors (in counts).
C	(SPR 2738)  Fred Shuman, STX,  1989 Feb 20.
C
C       Spr 3254, Dwellplot has to fetch the reference data from the
C                  Reference data archive.
C                  H. Wang, STX, 1990 Mar. 2. 
C-----------------------------------------------------------------------------
	Implicit None

	byte oldbuff(1), null(2), buff(450)
	integer*2 nrec, atempcnts_mf1(32), btempcnts_mf1(32),
     .			atempcnts_mf2(32), btempcnts_mf2(32),
     .			dwell_loc
	Integer*4 i, n, jpoly, Config_lun, Con_Lun, Btime(2), rcode
	integer*4 kounts_cal(4), kounts_grt(64), jside, jcur, index_offset
	real*4 ohms(1), resistance(64), polycoefs(3,4)
	character*14 time, timetags(1)

	equivalence (buff(1),  time),
     .		    (buff(101), atempcnts_mf1(1)),
     .		    (buff(165), btempcnts_mf1(1)),
     .		    (buff(317), atempcnts_mf2(1)),
     .		    (buff(381), btempcnts_mf2(1))
c
c  First calculate polynomials for counts-to-ohms conversion for each side
c  and current range. Do this on first pass only.
c
	do n=1,450
	   buff(n) = oldbuff(n)
	enddo
C
           call CT_GMT_TO_BINARY(time, Btime)
           rcode = 0  
	   call fep_grtcoeffs (polycoefs, Config_lun, Con_Lun,Btime,rcode)
           If (rcode .ne. 0) return
              
c
c  Put timetags into plotting buffer, one per record
c
	timetags(nrec) = time
c
c  Fill the 64-element array of counts with the appropriate part of the
c  housekeeping file (containing the GRT counts). Take either A or B side
c  depending on flag, "jside".  The last two values in each MjF are the 
c  first two values generated.  (SPR 2922)
c
	If (jside .eq. 0) Then
	   kounts_grt(1)  = atempcnts_mf1(31)
	   kounts_grt(2)  = atempcnts_mf1(32)
	   kounts_grt(33) = atempcnts_mf2(31)
	   kounts_grt(34) = atempcnts_mf2(32)
	   Do i=3,32
	      kounts_grt(i)    = atempcnts_mf1(i-2)
	      kounts_grt(i+32) = atempcnts_mf2(i-2)
	   End Do
	Else
	   kounts_grt(1)  = btempcnts_mf1(31)
	   kounts_grt(2)  = btempcnts_mf1(32)
	   kounts_grt(33) = btempcnts_mf2(31)
	   kounts_grt(34) = btempcnts_mf2(32)
	   Do i=3,32
	      kounts_grt(i)    = btempcnts_mf1(i-2)
	      kounts_grt(i+32) = btempcnts_mf2(i-2)
	   End Do
	End If
c
c  First determine which set of polynomial coefficients to use; jpoly=1 for
c  A/low current, 2 for A/high current, 3 for B/low current, etc.
c
	jcur = dwell_loc/16 + 1
	jpoly = jcur + 2*jside
c
c  Now convert the 64-point array of counts into ohms based on the polynomial
c  coefficients calculated the first time round.  (If this GRT is a cal 
c  resistor, just return it, converted from I*4 to R*4.  SPR 2738.)
c
	index_offset = (nrec - 1) * 64
	If ( dwell_loc .Lt. 10  .Or.
	2   (dwell_loc .Gt. 13  .And.  dwell_loc .Lt. 26)  .Or.
	3    dwell_loc .Gt. 29 ) Then
	   Do n=1,64
	      ohms(index_offset+n) = polycoefs(2,jpoly) + 
	2                            polycoefs(1,jpoly) * kounts_grt(n) + 
	3                            polycoefs(3,jpoly) * kounts_grt(n)**2 
	   End Do
	Else
	   Do n=1,64
	      ohms(index_offset+n) = kounts_grt(n)
	   End Do
	End If

	Return

	End
