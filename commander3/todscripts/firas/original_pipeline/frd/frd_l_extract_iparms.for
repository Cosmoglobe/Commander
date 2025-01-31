
	Integer*4 Function FRD_L_EXTRACT_IPARMS ( Buffer,
	1	  Parm, Lenparm, Limits )
!
!	Program Name : FRD_L_EXTRACT_IPARMS.For
!
!	Programmer  : Shirley M. Read, STX, January 1988.
!
!	Description : FRD_L_EXTRACT_IPARMS takes an input string of 80
!		      characters and extracts parameters, delimited
!	              by a blank ( or blanks ). The first parameter
!	              is assumed to be the name of the RDL field name 
!		      of one of the science record's science data related 
!		      information words or attitude limits. Two integer 
!		      word limits must follow : Red and Yellow. The 
!		      maximum length allowed for the first parm is 20 
!		      characters. The maximum allowed characters per number
!		      is 14.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 	! Input character string
	Character*20	Parm	! First parameter
	Integer*4	Lenparm   ! Length of first character parm
	Integer*2       Limits(2) ! Red and yellow limits
	
!	Include files

	Include '($SSdef)'

!	Functions

	integer*4 lib$skpc
	integer*4 lib$locc
	integer*4 ots$cvt_ti_l	! Convert signed integer text to integer
	
!	Local variables

	Integer*4 fstatus	! Function Processing status
	Integer*4 status 	! Returned system status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 ipos1, ipos2  ! Relative positions on string
	Integer*4 pstart, pstop  ! Positions on string
	Character*15 char_string
	Character*20 blank /'                    '/
	Integer*4 numlen        ! Length of string containing a number
	Integer*4 ix            ! Index
	Integer*4 start         ! Start position
	Integer*4 mask / '00000010'x /

!	Set status to success.

	fstatus = success

!	Get the first parm.

	parm(1:20) = blank(1:20)
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
		write (6,150) 
 	    endif
	endif
	if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1 
	    parm = buffer(pstart:pstop)
D	    write(6,900) parm, lenparm, pstart, pstop
	endif

!	Get the limits.

	do ix = 1, 2
	    if ( fstatus .eq. success ) then
		pstart = pstop + 1
		ipos1 = lib$skpc ( ' ', buffer(pstart:80))
		if ( ipos1 .eq. 0 ) then
		    fstatus = error
      		    write(6,200)
		else
		    pstart = pstart + ipos1 - 1
		    ipos2 = lib$locc ( ' ', buffer(pstart:80))
		    if ( ipos2 .eq. 0 .and. ix .ne. 2 ) then
			fstatus = error
		        write(6,230) 
 		    elseif ( ipos2 .eq. 0 ) then
			pstop = 80
			numlen = 81 - pstart
		    else
			pstop = pstart + ipos2 - 2
			numlen = pstop - pstart + 1
		    endif		      
		endif
		if ( numlen .gt. 14 ) then
		   fstatus = error
	           write(6, 300) numlen
	        endif		
		if ( fstatus .eq. success ) then
	            char_string(1:15) = blank(1:15)
		    start = 16 - numlen 
D		    write(6,925) ipos1,ipos2,numlen,pstart,pstop,start
		    char_string(start:15) =  buffer(pstart:pstop)
D		    write(6,950) char_string
	 	    status = ots$cvt_TI_L( char_string, limits(ix),
	1		     %val(2), mask )
	            if (status .eq. SS$_Normal ) then
			write(6,350) parm, ix, limits(ix)
		    else
			fstatus = error
	                write(6,400) status
		    endif
		endif
	    endif	! fstatus is success
	enddo		! do ix = 1, 2

!	Production Formats:
 100	format (1x,'Error: Input string is blank')
 150	format (1x,'Error: Input first parameter runs to end of line')
 200	format (1x,'Error: End of line on attempt to decode number')
 230    format (1x,'Error: Missing parameters on line.')
 300    format (1x,'Error: Number of characters in number= ',i4,
	1	' is too long. Limit is 14.')
 350    format (1x,'Parm= ',a,' Index= ',i2,' Limit Value = ',i8)
 400    format (1x,'Error: Bad return from OTS$CVT_TI_L, Status= ',z8.8)

!	Debug Formats:
 900    format(1x,'Parm= ',a, ' Parm length= ',i4,
	1	' pstart= ',i4,' pstop= ',i4)
 925        format(1x,'ipos1=',i4,' ipos2=',i4,' numlen=',i4,
	1	   ' pstart=',i4,' pstop=',i4,' start=',i4)
 950        format(1x,'Char_String=',a)

	FRD_L_Extract_IParms = fstatus
	return
	end

