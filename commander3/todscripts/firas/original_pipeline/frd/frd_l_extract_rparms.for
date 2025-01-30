
	Integer*4 Function FRD_L_EXTRACT_RPARMS ( Buffer,
	1	  Parm, Lenparm, RLimits, Name, Engr )
!
!	Program Name : FRD_L_EXTRACT_RPARMS.For
!
!	Programmer  : Shirley M. Read, STX, January 1988.
!
!	Description :   FRD_L_EXTRACT_RPARMS takes an input string of 80
!		      characters and extracts several parameters, delimited
!	              by a blank ( or blanks ). The parameters appear in
!		      sets of two records in the IGSE Firlims.Inp file. The
!		      first parameter on the first record of the set is 
!		      assumed to be the IGSE STOL database name of an 
!		      engineering word to be checked. These names are
!		      contained in the FDQ_Firnames.Txt file which provides
!		      a cross reference to the FDQ_Eng archive record fields
!		      starting with the En_Analog.Temp_Ctrl field. The database
!		      name will be followed by the word 'RANGE' and then by
!		      either 'ENGR' for engineering units ( real numbers) or
!		      'COUNTS' for unconverted count units ( integer numbers).
!		      The limits of database words in counts are not to be 
!		      checked under the current design. The flag 'Engr' will
!		      be set to false and the limits on the next record will
!		      be bypassed. If the units are engineering, the next 
!		      record must contain four real number ( decimal or 
!		      exponential format ) word limits: Red Low, Yellow Low, 
!		      Yellow High and Red High. The flag 'Name' will indicate
!		      which type of record is to be decoded. 
!	                On the first record, the maximum length allowed for
!		      each parameter is 20 characters. The maximum allowed 
!	              characters per number is 15. The database names
!		      are probably only 8 characters maximum.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 	   ! Input character string.
	Character*20	Parm	   ! First parameter.
	Integer*4	Lenparm    ! Length of first character parm.
	Real*4          RLimits(4) ! Red and yellow, high and low limits.
	Logical*1       Name       ! Flag for record with DB name if true,
			  	   ! four numbers if false.
	Logical*1       Engr       ! Flag set to indicate if Eng units or not.
	
!	Include files

	Include '($SSdef)'

!	Functions

	integer*4 lib$skpc
	integer*4 lib$locc
	integer*4 ots$cvt_t_f	! Convert numeric text to floating
	
!	Local variables

	Integer*4 fstatus	! Function Processing status
	Integer*4 status 	! Returned system status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 ipos1, ipos2  ! Relative positions on string
	Integer*4 pstart, pstop  ! Positions on string
	Character*16 char_string
	Character*20 arg        ! Buffer to contain char args 
	Character*20 blank /'                    '/
	Integer*4 numlen        ! Length of string containing a number
	Integer*4 ix            ! Index
	Integer*4 Start         ! Start position
	Integer*4 Mask / '00000017'x /

!	Set status to success.

	fstatus = success

!	If Name, get the first parm.

	IF ( Name) THEN

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
		write(6,50) pstart, buffer(1:40), buffer(41:80)
		write (6,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 20) then
	      write(6,160) lenparm
	      fstatus = error
	    else
	      parm = buffer(pstart:pstop)
D	      write(6,900) parm, lenparm, pstart, pstop
	    endif
	  endif     ! Fstatus is success

!	Get the type of units.

	  do ix = 1, 2
	     if ( fstatus .eq. success ) then
	        pstart = pstop + 1
	        arg(1:20) = blank(1:20)
	        ipos1 = lib$skpc ( ' ', buffer(pstart:80))
		if ( ipos1 .eq. 0 ) then
		    fstatus = error
      		    write(6,170)
		else
		    pstart = pstart + ipos1 - 1
		    ipos2 = lib$locc ( ' ', buffer(pstart:80))
		    if ( ipos2 .eq. 0 .and. ix .ne. 2 ) then
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
		if ( numlen .gt. 20 ) then
		   fstatus = error
	           write(6, 160) numlen
	        endif		
	     Endif  ! Fstatus is success
	  enddo
	  if ( fstatus .eq. success ) then
            arg = buffer(pstart:pstop)
	    if ( arg(1:4) .eq. 'ENGR' .or. arg(1:4) .eq. 'engr') then 
	      engr = .true.
	    else
	      engr = .false.
	    endif
          endif

	ELSE		! Limits record ( 4 numbers )

!	Get the limits.

	  pstop = 0
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
		if ( numlen .gt. 15 ) then
		   fstatus = error
	           write(6, 350) numlen
	        endif		
		if ( fstatus .eq. success ) then
	            char_string(1:16) = blank(1:16)
		    start = 17 - numlen 
D		    write(6,925) ipos1,ipos2,numlen,pstart,pstop,start
		    char_string(start:16) =  buffer(pstart:pstop)
D		    write(6,950) char_string
	 	    status = ots$cvt_t_f( char_string, rlimits(ix),
	1		     ,, mask )
	            if (status .eq. SS$_Normal ) then
			write(6,400) parm, ix, rlimits(ix)
		    else
			fstatus = error
	                write(6,450) status
		    endif
		endif
	    endif	! Fstatus is success
	  enddo		! do ix = 1, 2
	ENDIF		! Name record

!	Processing Formats:
  50    format(1x,'Parm start position= ',i8,' for Buffer:'/1x,a/1x,a)
 100	format (1x,'Error: Input string is blank')
 150	format (1x,'Error: Input first parameter runs to end of line')
 160	format(1x,'Error: Parameter length too long. Length= ',i4) 
 170    format(1x,'Error: Only one parameter on record.')
 180    format(1x,'Error: Missing parameters on line.')
 
 300	format (1x,'Error: End of line on attempt to decode number')
 350    format (1x,'Error: Number of characters in number= ',i4,
	1	' is too long. Limit is 16.')
 400    format (1x,'Parm= ',a,' Index= ',i2,' Limit Value = ',f16.6)
 450    format (1x,'Error: Bad return from OTS$CVT_T_F, Status= ',z8.8)

!	Debug Formats:
 900    format(1x,'Parm= ',a, ' Parm length= ',i4,
	1	' pstart= ',i4,' pstop= ',i4)
 925    format(1x,'ipos1=',i4,' ipos2=',i4,' numlen=',i4,
	1	   ' pstart=',i4,' pstop=',i4,' start=',i4)
 950    format(1x,'Char_String=',a)
 
!	Set function return status.

	FRD_L_Extract_RParms = fstatus

	return
	end

