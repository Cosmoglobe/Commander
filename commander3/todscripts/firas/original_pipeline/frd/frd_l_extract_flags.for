
	Integer*4 Function FRD_L_Extract_Flags ( Buffer, Flag,
	1	  Flaglen, Number, Parm )
!
!	Program Name : FRD_L_Extract_Flags.For
!
!	Programmer  : Shirley M. Read, STX, January 1988.
!
!	Modified:     Shirley M. Read, STX, May 1988.
!		      Modified to parse new input Fex_Limflags file which
!		      was redesigned in accord with the redesign of the flags
!		      to enable/disable individual FIRAS instrument components.
!
!	Description :   FRD_L_Extract_Flags takes an input string of 80
!		      characters and extracts parameters, delimited
!	              by a blank ( or blanks ). The first parameter
!		      is assumed to be the name of the flag. It must match
!		      one of the names in the included FDQ_Qualflag_Names.Txt.
!		      These flags correspond to one of the 110 quality flags
!		      in the science record's data_quality flags.
!		        The second parameter is assumed to be the number of
!	 	      input quality flags for the flags array field. The
!		      remainder of parameters are the values of the flags in
!		      the flag field array. The values must be in order. 
!	                The maximum length allowed for each flag name is 15.
!	              The maximum length allowed for each number field is 4.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 		! Input character string
	Character*15    Flag    	! Flag name
	Integer*2	Flaglen 	! Length of flag name
	Integer*2       Number 		! Number of values expected
	Byte		Parm(20)	! Values of flags	
	
!	Include files

	Include '(FUT_Qualflag_Names)'
	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Skpc
	Integer*4 Lib$Locc
	Integer*4 Ots$Cvt_TI_L	! Convert signed integer text to integer
	Integer*4 Ots$Cvt_TZ_L	! Convert hexadecimal text to integer or byte
	
!	Local variables

	Integer*4 fstatus	 ! Function Processing status
	Integer*4 status 	 ! Returned system status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 ipos1, ipos2   ! Relative positions on string
	Integer*4 pstart, pstop  ! Positions on buffer string
	Integer*4 lenparm	 ! Length of parameter string
	Character*5 char_string
	Character*20 blank /'                    '/
	Integer*4 numlen         ! Length of string containing a number
	Integer*4 ix             ! Index
	Integer*4 start          ! Start position on char_string
	Integer*4 mask / '00000010'x /
	Integer*2 numbers(2)     ! Storage for numbers in buffer
	Logical*1 hex		 ! Flag for hex or decimal characters
	Logical*1 found          ! Flag for finding valid qualflag name

!	Set status to success.

	fstatus = success

!	Initialize.

	start = 0
	pstart = 0
	pstop = 0

!	Get the first parm, the Flag field name.
	
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
	    if (lenparm .gt. 15) then
	      write(6,120) lenparm, buffer(pstart:pstop)
	      fstatus = error
	    else
	      flag = buffer(pstart:pstop)
D	      write(6,910) flag, lenparm, pstart, pstop
	      flaglen = lenparm
	      found = .false.
	      do ix = 1,110
		if ( flag .eq. qualflag_names(ix) ) found = .true.
	      enddo
	      if ( .not. found ) then
		fstatus = error
	        write(6, 130) flag
	      endif
	    endif
	endif     ! Fstatus is success

!	Get the second parm, the number of values expected on the line.
	
	if ( fstatus .eq. success ) then
	  char_string(1:5) = blank(1:5)
	  pstart = pstop + 1
	  ipos1 = Lib$Skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = Lib$Locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 4) then
	      write(6,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 6 - lenparm
	      char_string(start:5) = buffer(pstart:pstop)
D	      write(6,900) char_string, lenparm, start, pstart, pstop
	      status = Ots$Cvt_TI_L(char_string, number, %val(2), mask)   
	      if ( status .eq. SS$_normal) then
D	        write(6,960) number
	      else
	        write(6,400) status
	        fstatus = error
	      endif
	    endif
	  endif     ! Fstatus is success from parse of number on line
	endif	    ! Fstatus is success on first parm

!	Get the flag value parameters from the buffer.

	do ix = 1, number
	    if ( fstatus .eq. success ) then
		pstart = pstop + 1
		ipos1 = Lib$Skpc ( ' ', buffer(pstart:80))
		if ( ipos1 .eq. 0 ) then
		    fstatus = error
      		    write(6,200)
		else
		    pstart = pstart + ipos1 - 1
		    ipos2 = Lib$Locc ( ' ', buffer(pstart:80))
		    if ( ipos2 .eq. 0 .and. ix .ne. number ) then
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
!	Check to see if parm was input in hex or decimal characters.

	        if ( fstatus .eq. success ) then
		  if (buffer(pstart:pstart) .eq. 'X' ) then
		    hex = .true.
	            pstart = pstart + 1
	            numlen = numlen - 1
	          else
	            hex = .false.
	          endif
	        endif

!	Check length of actual character number.

		if ( numlen .gt. 4 ) then
		   fstatus = error
	           write(6, 300) numlen
	        endif		
!	Convert the characters to a number.

		if ( fstatus .eq. success ) then
	          char_string(1:5) = blank(1:5)
		  start = 6 - numlen 
D		  write(6,925) ipos1,ipos2,numlen,pstart,pstop,start
		  char_string(start:5) =  buffer(pstart:pstop)
D		  write(6,950) char_string
		  if (hex) then
		    status = Ots$Cvt_TZ_L(char_string, parm(ix), %val(1),)
	            if (status .ne. SS$_Normal ) then
			fstatus = error
	                write(6,425) status
		    endif
		  else
	 	    status = Ots$Cvt_TI_L( char_string, parm(ix),
	1		     %val(1), mask )
	            if (status .ne. SS$_Normal ) then
		        fstatus = error
	                write(6,400) status
		    endif
		  endif
	          if (status .eq. SS$_Normal ) then
			write(6,350) flag, ix, parm(ix), parm(ix)
		  endif
	       endif    ! fstatus is success
	    endif	! fstatus is success
	enddo		! do ix = 1, number


!	Production Formats:
 100	format (1x,'Error: Input string is blank')
 110	format (1x,'Error: Input first parameter runs to end of line'
	1	/ 1x, a )
 120    format (1x,'Error: Flag name parameter too long. Length= ',i4
	1	/1x, 'Parm= ', a)
 130    format (1x,'Error: Input flag field name ',a,' is invalid')
 140    format (1x,'Error: No numbers on line. Rest of line blank')
 150    format (1x,'Error: Input second parameter runs to end of line')
 160	format (1x,'Error: Second parameter length too long. Length= ',i4) 

 200	format (1x,'Error: End of line on attempt to decode number')
 230    format (1x,'Error: Missing parameters on line.')
 300    format (1x,'Error: Too many characters in flag value: ',i4,
	1	' Limit is 4.')
 350    format (1x,'Flag= ',a,' Index= ',i4,' Value(dec)= ',i6,
	1	' Value(hex)= ',z2.2)
 400    format (1x,'Error: Bad return from OTS$CVT_TI_L, Status= ',z8.8)
 425    format (1x,'Error: Bad return from OTS$CVT_TZ_L, Status= ',z8.8)
 500    format (1x,'Error: Input array position is out of bounds. ',
	1 	'Array position= ',i10)
 550    format (1x,'Error: Invalid value for limit flag.'/
	1	1x,'Value must be 0 or -1 (false or true). ',
	2	'This record has value=',i8)

!	Debug Formats:
 900    format(1x,'Parm= ',a, ' Parm length= ',i4,' Outstring start=',i4,       
	1	' pstart= ',i4,' pstop= ',i4)

 910    format(1x,'Parm= ',a, ' Parm length= ',i4,
	1	' pstart= ',i4,' pstop= ',i4)
 
 925    format(1x,'Ipos1=',i4,' ipos2=',i4,' numlen=',i4,
	1	   ' pstart=',i4,' pstop=',i4,' start=',i4)
 950    format(1x,'Char_String=',a)
 960    format(1x,'Number of flags to be stored= ',i4)

	FRD_L_Extract_Flags = fstatus

	Return
	End
