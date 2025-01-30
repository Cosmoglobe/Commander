        subroutine fsd_astroplots_timeaxis (nbreak,xlabel,tstart,tend,
     .           deltat,t_unit,c_flag)

C------------------------------------------------------------------------------
C
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C------------------------------------------------------------------------------
C Changes:
C
C	SPR 5697, Remove dependency on MUT_YD_TO_YMD; reorganize time
C	    convertion.  R. Kummerer / STX, February 1, 1990.
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C------------------------------------------------------------------------------

        implicit none
c
	include '(fsd_astroplots)'

	Integer         *2      nbreak(50)
	Character       *32     xlabel
        Real            *8      tstart
	Real            *8      tend
	Real            *8      deltat
	Real            *8      t_unit
	Logical         *4      c_flag

	Integer         *4      i
	Integer         *4      ii
	real            *8      t
	Real		*4	gap_thresh
	real            *8      tjump

	Real            *8      aut_adt2t68

	if (.not.c_flag)then

c	GET start and end time

	   tstart = aut_adt2t68(gmt(1,1))
	   tend   = aut_adt2t68(gmt(1,num))
c
	    if ((tend-tstart) .ge. 1.08e+04)then
	       t_unit=3600.0
	       deltat=60.0
	       xlabel = 'Time (Hour)'
               gap_thresh=5./60.
	    else
	       t_unit=60.
	       deltat=1.0
	       xlabel = 'Time (Minute)'
	       gap_thresh=5.
	    endif
	   timeaxis(1)=0.0
	endif

	Do i=2,num
	   t = aut_adt2t68(gmt(1,i))
	   timeaxis(i)=(t-tstart)/t_unit
	End Do
c	
c	Check for breaks in timeaxis
c
	ii=0
	Do i=2,num
	   tjump = timeaxis(i) - timeaxis(i-1)
	   If (tjump .Gt. gap_thresh) Then
	      ii=ii+1
	      nbreak(ii)=i-1
	   End If
	End Do
	     NBREAK(II+1)=NUM
	Return
	End
