
	Integer*4 Function FRD_L_Store_GRT( GRTParm, Number, RLimits, GRTs)

!	Program Name : FRD_L_Store_GRT
!
!	Programmer: Shirley M. Read, STX, January 1988
!
!       Modified: Shirley M. Read, STX, May 1988
!                 Changed the storing of GRT limits into GRT array to
!	          skip the calibrator resistors which have no linits
!		  defined.
!
!	Program Description:
!	    This subroutine stores a set of engineering GRT limits 
!	  ( Red Low, Yellow Low, Yellow High and Red High ) into the GRTs 
!	  array if the number is within the bounds of the GRTs array. 
!
!----------------------------------------------------------------------------

	Implicit None

!	Include files.

	Include '(FUT_Convbuff)'

!	Passed Parameters.

	Character*8	GRTParm       ! Input GRT field parameter
	Integer*2       Number        ! GRTs Array position to be stored
	Real*4    	RLimits(4)    ! Set of Red and yellow limits
	Real*4		GRTs(ngrt,4)  ! GRT limits for the GRT
				      ! engineering database limits.
				      ! Order Red Low, Yellow Low,
		       		      ! Yellow High and Red High
!	Local variables

	Character*8	grt_field(4) / 'A_LO_GRT', 'A_HI_GRT',
	1			       'B_LO_GRT', 'B_HI_GRT' /
	Integer*2       grt_pos(4) / 0, 12, 24, 36 /
	Integer*2       position      ! Output GRT array position
	Integer*4 	retstat	      ! Program Processing status
	Integer*4 	success / 1 /, error / 2 /
	Integer*4 	ix            ! Index
	Logical*1 	proceed       ! Proceed flag
	Logical*1       found	      ! Flag for match found

!	Set status to success.

	retstat = success

!       If the input number parameter <= 12, set the proceed flag to true .

	If ( number .le. 12 ) then
	  proceed = .true.
	else
	  proceed = .false.
	  write(6,100) Number
	endif

!	Find the position of the GRT in the output array.

	If ( proceed ) then
	  found = .false.
	  Do ix = 1, 4
	    If ( grtparm .eq. grt_field(ix) ) then
	      found = .true.
	      position = grt_pos(ix) + number
	    Endif
	  Enddo
	  If ( .not. found ) then
	    proceed = .false.
	    retstat = error
	    Write(6,200) grtparm
	  Endif
	Endif	! Proceed
	
!	If proceed, store the limits in the GRTs array.

	If ( Proceed ) then
	    Do ix = 1, 4
	      GRTs(position,ix) = Rlimits(ix)
	    Enddo
	    Write (6,300) GRTParm, Number, Rlimits
	Endif

!	Formats
 100    format(1x,'GRT input array position is out of bounds. ',
	1	'Number= ',i10)	
 200    format(1x,'GRT field name ',a,' is invalid.')
 300    format(1x,'GRT field ',a,' has limits stored in position ',
	1	i3,' of GRTs array.'
	1	/1x,'Limits = ',4(f16.5,1x))

	FRD_L_Store_Grt = retstat

	Return
	End

