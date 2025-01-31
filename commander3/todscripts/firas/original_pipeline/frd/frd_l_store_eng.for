
	Integer*4 Function FRD_L_Store_Eng ( Parm, Lenparm, RLimits,
	1	  Analogs )

!	Program Name : FRD_L_Store_Eng
!
!	Programmer: Shirley M. Read, STX, January 1988

!	Program Description:
!	    This subroutine stores a set of engineering analog limits 
!	  ( Red Low, Yellow Low, Yellow High and Red High ) into the analogs 
!	  array if the parameter name matches one of the STOL engineering 
!	  database names in FDQ_Firnames.Txt.
!
!----------------------------------------------------------------------------

	Implicit None

!	Include files

	Include '(FUT_Convbuff)'
	Include '(FUT_Firnames)'
	Include '(FUT_Grtpars)'
	Include '(FUT_Lckeypars)'

!	Passed Parameters.

	Character*20	Parm	      ! First parameter
	Integer*4	Lenparm       ! Length of first character parm
	Real*4    	RLimits(4)    ! Set of Red and yellow limits
	Real*4		Analogs(npoly,4) ! Analog limits for the STOL
				      ! engineering database words.
				      ! Order Red Low, Yellow Low,
		       		      ! Yellow High and Red High
!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 success / 1 /, error / 2 /
	Real*4    templims(4)   ! Temporary storage for limits
	Integer*4 ix, jx        ! Indices
	Integer*4 match         ! Var to store index of match


!	Set status to success.

	retstat = success

!       Check the input parameter and set the match index if found.

	Match = 0
	Do ix = 1, npoly
	  If ( Parm(1:lenparm) .eq. Poly_Names(ix)(1:lenparm) ) then 
	    Match = ix
	  Endif
	Enddo

!	If match found, store the limits in the Analogs array.

	If ( Match .ne. 0 ) then

!	Check for inverse curves.

	  If ( ( Rlimits(1) .gt. Rlimits(2) ) .OR.
	1      ( Rlimits(2) .gt. Rlimits(3) ) .OR.
	2      ( Rlimits(3) .gt. Rlimits(4) ) ) Then

!	Check for error in order, i.e., limits must be monotonically 
!	increasing or monotonically decreasing.

	    If ( ( Rlimits(1) .lt. Rlimits(2) ) .OR.
	1        ( Rlimits(2) .lt. Rlimits(3) ) .OR.
	2        ( Rlimits(3) .lt. Rlimits(4) ) ) Then
	      Retstat = error
	      Write (6,200) Parm(1:lenparm), Rlimits
	    Else

!	Reverse the limits.

 	      Do ix = 1, 4
 		Templims(ix) = Rlimits(ix)
 	      Enddo
	      Do ix = 1, 4
 		jx = 5 - ix
		Rlimits(ix) = Templims(jx)
	      Enddo
	      Write (6,300) Parm(1:lenparm)
	    Endif
	  Endif		! Inverse curve detected.

!	Store the limits.

	  If ( Retstat .eq. success ) Then
	      Do ix = 1, 4
	        Analogs(match,ix) = Rlimits(ix)
	      Enddo
	      Write (6,100) Parm(1:lenparm), Rlimits
	  Endif
	Else
	  Write(6,400) Parm(1:lenparm)
	Endif

	FRD_L_Store_Eng = retstat

	return

!	Formats
 100    format(1x,'Parameter ',a,' stored in Analog array.'
	1	/1x,'Limits = ',4(f14.5,1x)) 
 200    format(1x,'Error: Parameter ',a,' does not have either ',
	1	'an increasing or a decreasing set of limits.'/
	2	1x,'Limits (in order of appearence on line ) =',
	3	4(1x,f14.5))
 300    format(1x,'Parameter ',a,' has an inverse curve. '/
	1      1x,'Program will reverse limits for limits database.')
 400	format(1x,'Parameter ',a,' not found in Match Table.') 
	
	end

