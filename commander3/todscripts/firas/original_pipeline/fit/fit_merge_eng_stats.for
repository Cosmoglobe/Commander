	Integer*4 Function   FIT_Merge_Eng_Stats (stats_buff, stats_gmts,
	2                                         bin_size,
	3                                         starting_time, ending_time,
	4                                  nrecs, bin_buff, bin_gmts, nbins)

C------------------------------------------------------------------------
C    PURPOSE: Merge statistics on a given engineering field into bins, 
C		given a time-bin size
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman, STX
C            1987 Nov 6
C
C    INVOCATION:       STATUS = FIT_MERGE_ENG_STATS (STATS_BUFF, STATS_GMTS,
C					BIN_SIZE, STARTING_TIME, ENDING_TIME,
C					NRECS, BIN_BUFF, BIN_GMTS, NBINS)
C
C    INPUT PARAMETERS:
C	STATS_BUFF(max,
C	     num_stats_types)	 R*4		ENG STATS data to be merged.
C	STATS_GMTS(max,3)	Ch*14		GMTs associated with each orbit:
C						    1=1st, 2=Ave, 3=Last IFG.
C	BIN_SIZE		 R*4		User-chosen time bin (sec).
C	STARTING_TIME		Ch*14		GMT at start of plot range.
C	ENDING_TIME		Ch*14		GMT at end of plot range.
C	NRECS			 I*4		Number of orbit records.
C
C    OUTPUT PARAMETERS:
C	BIN_BUFF(max,
C	     num_stats_types)	 R*4		ENG STATS merged into bins.
C	BIN_GMTS(max, 3)	Ch*14		GMTs for each bin:
C	NBINS			 I*4		Number of bins formed.
C						    1=1st, 2=Mean, 3=Last.
C		Note: "max" is fac_max_stats_recs, declared in FUT_PARAMS.TXT
C				= 1000 as of 1988 Jun 6
C    SUBROUTINES CALLED:
C	CT_GMT_To_Binary
C	CT_Binary_To_GMT
C	AUT_Dfloat2ADT
C	AUT_ADT2Dfloat
C
C    COMMON VARIABLES USED:  None
C
C    INCLUDE FILES: 
C	FIT_Invoc
C	FIT_Menu
C	FUT_Error
C	FUT_Params
C	$SSDEF
C
C-----------------------------------------------------------------------
C	Revised:
C	                                        F. Shuman,  1987 Nov 19
C	   Replace Vecplt with PLT for plots.   F. Shuman,  1988 Mar 31
C	   SPR 3139.  Reject records with FDQ's "bad record" flag.
C	                                        F. Shuman,  1989 Feb 23
C
C       version 4.5.1 QFix 530, SPRs 4669, 4670.  F. Shuman,  STX,  1989 Oct 2.
C           Successive cropping of endtime; shuffling of data.
C           Changes to FIT.for, FIT_Extract_Eng_Stats.for,
C           FIT_Merge_Eng_Stats.for, and FIT_Data.txt.
C-----------------------------------------------------------------------

	Implicit	None

	Include		'(FIT_Invoc)'
	Include		'(FIT_Menu)'
	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FIT_Data)'
	Include		'($SSDEF)'

	Integer		*4	status

	Real		*4	stats_buff(fac_max_stats_recs, num_stats_types)
	Character	*14	stats_gmts(fac_max_stats_recs, 3)
	Real		*4	bin_size
	Character	*14	starting_time
	Character	*14	ending_time
	Integer		*4	nrecs
	Real		*4	bin_buff(fac_max_stats_recs, num_stats_types)
	Character	*14	bin_gmts(fac_max_stats_recs, 3)
	Integer		*4	nbins

	Real		*8	stats_time
	Real		*8	bin_time
	Real		*8	sum_time

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	ADTime(2)
	Real		*4	minval
	Real		*4	maxval
	Integer		*4	nsum_wts
	Real		*4	sum_vals
	Real		*4	sum_sqrs
	Real		*8	ave_time
	Integer		*4	nx0
	Real		*4	x1
	Real		*4	x2
	Real		*4	mean
	Real		*4	s
	Real		*4	xoffset

	Integer		*4	CT_GMT_To_Binary
	Integer		*4	CT_Binary_To_GMT
	Integer		*4	AUT_Dfloat2ADT
	Real		*8	AUT_ADT2Dfloat

	External	FIT_Normal
	External	FIT_NoGooDat

c i = bin counter,  j = orbit counter (of the incoming data)
c
	i = 1
	j = 1
c
c Skip to start of plot timerange, then to 1st orbit with any data
c
	Do While ((stats_gmts(j,1) .Lt. starting_time
	2          .Or. JNINT(stats_buff(j,5)) .Eq. 0)
	3         .And. j .Lt. nrecs)
	   j = j + 1
	End Do
c
c While within plot timerange
c
	Do While (stats_gmts(j,1) .Le. ending_time .And. j .Le. nrecs)
c
c If the value was written by FDQ as a "sentinel" (= -9999.), flagging a
c  bad orbit record, then reject it.  Stats_buff(j,3) is the orbit-record
c  minimum, so data written by FDQ before its corresponding fix (SPR 3139)
c  will be rejected if any of the values that went into it were flagged.
c  Also check (j,1), the mean value, and (j,4), the maximum, because trends
c  data has been found with mean =~ -3400, but min =~ 2.
c
	   If (stats_buff(j,1) .Gt. -100.  .And.
	2      stats_buff(j,3) .Gt. -100.  .And.
	3      stats_buff(j,4) .Gt. -100.) Then
c
c Initialize the statistics
c
	      nsum_wts = 0
	      sum_vals = 0.
	      sum_sqrs = 0.
	      minval   = stats_buff(j,3)
	      maxval   = stats_buff(j,4)
	      sum_time = 0.

	      bin_gmts(i,1) = stats_gmts(j,1)
	      bin_gmts(i,3) = stats_gmts(j,3)
c
c Process one bin of data
c
	      Call CT_GMT_To_Binary ( stats_gmts(j,1), ADTime )
	      stats_time = 1.D-7 * AUT_ADT2DFloat ( ADTime )
	      Call CT_GMT_To_Binary ( bin_gmts(i,1), ADTime )
	      bin_time = 1.D-7 * AUT_ADT2DFloat ( ADTime )
c
c While within bin interval and timerange ...
c    (first set up an offset in x to guard against roundoff err in 
c     later calculation of sample std. devn.)
c
	      xoffset = stats_buff(j,1)

	      Do While (stats_time .Le. bin_time + bin_size
	2               .And. stats_gmts(j,1) .Le. ending_time
	3               .And. j .Le. nrecs)
c
c !# of records that went into orbit j; JNINT=nearest int
c
	         nx0 = JNINT ( stats_buff(j,5) )
c
c Skip any orbit with no data
c
	         If (nx0 .Gt. 0) Then

	            Call CT_GMT_To_Binary ( stats_gmts(j,2), ADTime )
	            stats_time = 1.D-7 * AUT_ADT2DFloat ( ADTime )
	            sum_time = sum_time + nx0*stats_time
c
c Ave value during orbit, less the offset
c
	            x1 = stats_buff(j,1) - xoffset
c
c Sample std devn during orbit
c
	            If (nx0 .Gt. 1) Then
	               x2 = stats_buff(j,2)
	            Else
	               x2 = 0.
	            End If

	            nsum_wts = nsum_wts + nx0
	            sum_vals = sum_vals + nx0*x1
	            sum_sqrs = sum_sqrs + (nx0 - 1)*(x2**2) + nx0*(x1**2)
	            minval   = AMIN1 ( minval, stats_buff(j,3) )
	            maxval   = AMAX1 ( maxval, stats_buff(j,4) )

	         End If

	         j = j + 1

	         Call CT_GMT_To_Binary ( stats_gmts(j,1), ADTime )
	         stats_time = 1.D-7 * AUT_ADT2DFloat ( ADTime )

	      End Do

	      nx0 = nsum_wts
	      x1  = sum_vals
	      x2  = sum_sqrs

	      If (nx0 .Gt. 0) Then
	         mean = x1 / nx0  +  xoffset
	         ave_time = sum_time/nx0
	      Else
	         mean = 0.
	         ave_time = 0.
	      End If
c
c Check for variance being 'computationally' < 0 (roundoff err can cause this)
c
	      If (nx0 .Gt. 1 .And. nx0*x2 - x1*x1 .Ge. 0.) Then
	         s = SQRT ( (nx0*x2 - x1*x1)/nx0/(nx0 - 1) )
	      Else
	         s = 0.
	      End If

	      bin_buff(i,1) = minval
	      bin_buff(i,2) = mean
	      bin_buff(i,3) = s
	      bin_buff(i,4) = maxval
	      bin_buff(i,5) = nsum_wts

	      Call AUT_DFloat2ADT ( 1.D+7 * ave_time, ADTime )
	      Call CT_Binary_To_GMT ( ADTime, bin_gmts(i,2) )
c
c Increment the bin #
c
	      i = i + 1

	   Else

	      j = j + 1

	   End If         ! stats_buff(j,3) .Gt. -9990.
c
c Skip to next orbit with any data
c
	   Do While (JNINT(stats_buff(j,5)) .Eq. 0
	2            .And. j .Lt. fac_max_stats_recs)
	      j = j + 1
	   End Do


	End Do             ! While within plot timerange

	nbins = i - 1

	If (nbins .Gt. 0) Then
	   FIT_Merge_Eng_Stats = %Loc(FIT_Normal)
	Else
	   FIT_Merge_Eng_Stats = %Loc(FIT_NoGooDat)
	End If

	Return
	End
