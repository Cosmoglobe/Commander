 
	Integer*4 Function FRD_L_EXTRACT_GPARMS ( Buffer,
	1	  Parm, Number, RLimits )
!
!	Program Name : FRD_L_EXTRACT_GPARMS.For
!
!	Programmer  : Shirley M. Read, STX, January 1988.
!
!	Description :   FRD_L_EXTRACT_GPARMS takes an input string of 80
!		      characters and extracts several parameters, delimited
!	              by a blank ( or blanks ) from the string. The maximum
!		      length of each parameter is 16 characters. 
!		        The first parameter is the GRT group field name for the 
!		      engineering limits record. It must have one of the four
!		      values: A_LO_GRT, A_HI_GRT, B_LO_GRT, B_HI_GRT. The
!		      second parameter is the GRT number within the group, i.e.,
!		      a number from 1 to 12 ( calibrator resistors are not
!		      included ) indicating the relative order of the GRT.
!		      Following are four real numbers ( decimal or 
!		      exponential format ) word limits: Red Low, Yellow Low, 
!		      Yellow High and Red High. 
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 	   ! Input character string.
	Character*8     Parm       ! GRT group field name.
	Integer*2       Number     ! Number of the GRT.
	Real*4          RLimits(4) ! Red and yellow, high and low limits.
	
!	Include files

	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Skpc
	Integer*4 Lib$Locc
	Integer*4 Ots$Cvt_TI_L  ! Convert numeric text to integer
	Integer*4 Ots$Cvt_T_F	! Convert numeric text to floating
	
!	Local variables

	Integer*4 lenparm       ! Length of number
	Integer*4 fstatus	! Function Processing status
	Integer*4 status 	! Returned system status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 ipos1, ipos2  ! Relative positions on string
	Integer*4 pstart, pstop  ! Positions on string
	Character*8 grt_field(4) / 'A_LO_GRT', 'A_HI_GRT',
	1			   'B_LO_GRT', 'B_HI_GRT' /
	Logical*1 found		 ! GRT name found in match table
	Character*17 char_string
	Character*20 blank /'                    '/
	Integer*4 numlen        ! Length of string containing a number
	Integer*4 ix            ! Index
	Integer*4 start         ! Start position on string arg for OTS call
	Integer*4 Mask  / '00000017'x /
	Integer*4 MaskI / '00000010'x /

!	Set status to success.

	fstatus = success

!	Get the first parm, the GRT group field name.
	
	ipos1 = lib$skpc ( ' ', buffer )
	if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,100) 
	endif
	if ( fstatus .eq. success ) then
	    pstart = ipos1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,110) buffer 
 	    endif
	endif
 
	if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 8) then
	      write(6,120) lenparm, buffer(pstart:pstop)
	      fstatus = error
	    else
	      parm = buffer(pstart:pstop)
D	      write(6,910) parm, lenparm, pstart, pstop
	      found = .false.
	      do ix = 1,4
		if ( parm .eq. grt_field(ix) ) found = .true.
	      enddo
	      if ( .not. found ) then
		fstatus = error
	        write(6, 130) parm
	      endif
	    endif
	endif     ! Fstatus is success

!	Get the second parm, the GRT number.
	
	if ( fstatus .eq. success ) then
	  char_string(1:17) = blank(1:17)
	  pstart = pstop + 1
	  ipos1 = lib$skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 16) then
	      write(6,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 18 - lenparm
	      char_string(start:17) = buffer(pstart:pstop)
D	      write(6,900) char_string, lenparm, start, pstart, pstop
	      status = ots$cvt_ti_l(char_string, number, %val(2), maski)   
	      if ( status .eq. SS$_normal) then
D	        write(6,960) number
	      else
	        write(6,425) status
	        fstatus = error
	      endif
	    endif
	  endif     ! Fstatus is success from parse of number on line

!	Check to see if GRT position number is valid, 1 to 12.

	  if ( fstatus .eq. success ) then
	    if( (number .lt. 1) .or. (number .gt. 12)) then
	       fstatus = error
	       write(6,430) number
	    endif
	  endif	    ! Fstatus is success after parse of number on line
	endif	    ! Fstatus is success on first parameter

!	Get the limits.

	  do ix = 1, 4
	    if ( fstatus .eq. success ) then
		pstart = pstop + 1
		ipos1 = lib$skpc ( ' ', buffer(pstart:80))
		if ( ipos1 .eq. 0 ) then
		    fstatus = error
      		    write(6,300)
		else
		    pstart = pstart + ipos1 - 1
		    ipos2 = lib$locc ( ' ', buffer(pstart:80))
		    if ( ipos2 .eq. 0 .and. ix .ne. 4 ) then
			fstatus = error
	                write(6,180)
		    elseif ( ipos2 .eq. 0 ) then
			pstop = 80
			numlen = 81 - pstart
		    else
			pstop = pstart + ipos2 - 2
			numlen = pstop - pstart + 1
		    endif		      
		endif
		if ( numlen .gt. 16 ) then
		   fstatus = error
	           write(6, 350) numlen
	        endif		
		if ( fstatus .eq. success ) then
	            char_string(1:17) = blank(1:17)
		    start = 18 - numlen 
D		    write(6,925) ipos1,ipos2,numlen,pstart,pstop,start
		    char_string(start:17) =  buffer(pstart:pstop)
D		    write(6,950) char_string
	 	    status = ots$cvt_t_f( char_string, rlimits(ix),
	1		     ,, mask )
	            if (status .eq. SS$_Normal ) then
			write(6,400) parm, number, ix, rlimits(ix)
		    else
			fstatus = error
	                write(6,450) status
		    endif
		endif
	    endif	! Fstatus is success
	  enddo		! do ix = 1, 2

!	Processing Formats:
 100	format (1x,'Error: Input string is blank')
 110	format (1x,'Error: Input first parameter runs to end of line'
	1	/ 1x, a )
 120    format (1x,'Error: GRT name parameter too long. Length= ',i4
	1	/1x, 'Parm= ', a)
 130    format (1x,'Error: Input GRT field name ',a,' is invalid')
 140    format (1x,'Error: No GRT position number. Rest of line blank')
 150    format (1x,'Error: Input second parameter runs to end of line')
 160	format (1x,'Error: Parameter length too long. Length= ',i4) 
 180    format (1x,'Error: Missing parameters on line.')
 300	format (1x,'Error: End of line on attempt to decode number')
 350    format (1x,'Error: Number of characters in number= ',i4,
	1	' is too long. Limit is 16.')
 400    format (1x,'GRT Field= ',a,' Number= ',i6,' Index= ',i2,
	1	' Limit Value = ',f16.6)
 425    format (1x,'Error: Bad return from OTS$CVT_TI_L, Status= ',z8.8)
 430    format (1x,'Error: Invalid GRT array position. Number= ',i10)
 450    format (1x,'Error: Bad return from OTS$CVT_T_F, Status= ',z8.8)

!	Debug Formats:
 900    format(1x,'Parm= ',a, ' Parm length= ',i4,' Outstring start=',i4,       
	1	' pstart= ',i4,' pstop= ',i4)
 910    format(1x,'Parm= ',a, ' Parm length= ',i4,
	1	' pstart= ',i4,' pstop= ',i4)
 925    format(1x,'ipos1=',i4,' ipos2=',i4,' numlen=',i4,
	1	   ' pstart=',i4,' pstop=',i4,' start=',i4)
 950    format(1x,'Char_String=',a)
 960    format(1x,'Number of Grt=',i8)

	FRD_L_Extract_GParms = fstatus

	Return
	End

