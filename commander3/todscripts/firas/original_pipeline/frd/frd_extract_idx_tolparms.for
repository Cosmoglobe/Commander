 
	Integer*4 Function FRD_Extract_IDX_Tolparms ( Buffer , perc_number ,
	1                       rpt_lun,report,grp_names, number, found)
!
!	Program Name : FRD_Extract_Idx_Tolparms.For
!
!	Programmer  : Quoc C  chung, STX, April 7 1989.
!
!	Description : FRD_extract_idx_tolparms takes an input string of 80
!		      characters and extracts Three parameters, delimited
!	              by a blank ( or blanks ) from the string. The maximum
!		      length of each parameter is 32 characters. 
!		      The first parameter is the group field name for the 
!		      engineering record. It must have one of the twenty two
!		      Names.The second parameter is the class number within
!                     the group, i.e., a number from 1 to 22 
!		      Following is Integer number from 0 - 100 (percentage).
!
!        version 4.4 : Quoc C Chung, STX, May 10 1989.
!                      SPR 3773 FRD_DEFINE_IDX_PARAMS WRITES REPORT TO
!                      FOR000.DAT. 
!---------------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 	   ! Input character string.
        Integer*2       perc_number! converted percentage value
	Logical*1       found	   ! GRT name found in match table
        Logical*1       report
        Integer*4       rpt_lun
	
!	Include files

	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Skpc
	Integer*4 Lib$Locc
	Integer*4 Ots$Cvt_TI_L  ! Convert numeric text to integer
	Integer*4 Ots$Cvt_T_F	! Convert numeric text to floating
	
!	Local variables

	Character*32    Tols_name  ! Tols group field name.
        Character*32    grp_names(22)
	Integer*2       Number     ! Class Number .
	Integer*2       namelen    ! length of tols group field name
	Integer*4       lenparm    ! Length of number
	Integer*4       fstatus	   ! Function Processing status
	Integer*4       status 	   ! Returned system status
	Integer*4       success / 1 /, error / 2 /
	Integer*4       ipos1, ipos2 ! Relative positions on string
	Integer*4       pstart, pstop  ! Positions on string
!
	Character*2 char_string2
	Character*3 char_string3
	Character*20 blank20 /'                    '/
	Integer*4 numlen        ! Length of string containing a number
	Integer*4 ix            ! Index
	Integer*4 start         ! Start position on string arg for OTS call
	Integer*4 Mask  / '00000017'x /
	Integer*4 MaskI / '00000010'x /

!	Set status to success.

	fstatus = success

!	Get the first parm, the group field name.
	
	ipos1 = lib$skpc ( ' ', buffer )
	if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,100) 
            if (report) write(rpt_lun,100)
	endif
	if ( fstatus .eq. success ) then
	    pstart = ipos1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,110) buffer 
                if (report) write(rpt_lun,110) buffer
 	    endif
	endif
 
	if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 32) then
	      write(6,120) lenparm, buffer(pstart:pstop)
	      if (report) write(rpt_lun,120) lenparm, buffer(pstart:pstop)
	      fstatus = error
	    else
	      tols_name = buffer(pstart:pstop)
D	      write(6,910) tols_name, lenparm, pstart, pstop
	      found = .false.
	      do ix = 1,22
		if ( tols_name .eq. grp_names(ix) ) then
                found = .true.
                goto 13
                endif
	      enddo
            
 13    continue
	      if ( .not. found ) then
		fstatus = error
	        write(6, 130) tols_name
	        if (report) write(rpt_lun, 130) tols_name
	      endif ! .not. found
            endif   ! length
	endif     ! Fstatus is success

!	Get the second parm, the class number.
	
	if ( fstatus .eq. success ) then
	  char_string2(1:2) = blank20(1:2)
	  pstart = pstop + 1
	  ipos1 = lib$skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,140) 
	    if (report) write(rpt_lun,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,150) 
		if (report) write (rpt_lun,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 2) then
D	      write(6,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 3 - lenparm
	      char_string2(start:2) = buffer(pstart:pstop)
D	      write(6,900) char_string, lenparm, start, pstart, pstop
	      status = ots$cvt_ti_l(char_string2, number, %val(2), maski)   
	      if ( status .eq. SS$_normal) then
D	        write(6,960) number
	      else
	        write(6,425) status
	        if (report) write(rpt_lun,425) status
	        fstatus = error
	      endif  ! status
	    endif    ! lenparm
	  endif     ! status is Normal from parse of number on line

!	Check to see if Class number is valid, 1 to 22.

	  if ( fstatus .eq. success ) then
	    if( (number .lt. 1) .or. (number .gt. 22)) then
	       fstatus = error
	       write(6,430) number
	       if (report) write(rpt_lun,430) number
	    endif
          endif	    ! Fstatus is success on second parameter
        endif	    ! Fstatus is success after parse  number on line
!	Get the percentage

	if ( fstatus .eq. success ) then
	  char_string3(1:3) = blank20(1:3)
	  pstart = pstop + 1
	  ipos1 = lib$skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    write(6,140) 
	    if (report) write(rpt_lun,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		write (6,150) 
		if (report) write (rpt_lun,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 3) then
	      write(6,160) lenparm
	      if (report) write(rpt_lun,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 4 - lenparm
	      char_string3(start:3) = buffer(pstart:pstop)
D	      write(6,900) char_string3, lenparm, start, pstart, pstop
	      status = ots$cvt_ti_l(char_string3, Perc_number, %val(2), maski)   
	      if ( status .eq. SS$_normal) then
D	        write(6,960) perc_number
	      else
	        write(6,425) status
	        if (report) write(rpt_lun,425) status
	        fstatus = error
	      endif
	    endif
	  endif     ! Fstatus is success 

!	Check to see if percentage is valid, 0 to 100.

	  if ( fstatus .eq. success ) then
	    if( (perc_number .lt. 0) .or. (perc_number .gt. 100)) then
	       fstatus = error
	       write(6,430) perc_number
	       if (report) write(rpt_lun,430) perc_number
	    endif
	  endif	    ! Fstatus is success 

	  endif	    ! Fstatus is success after parse of percentage on line


!	Processing Formats:
 100	format (1x,'Error: Input string is blank')
 110	format (1x,'Error: Input first parameter runs to end of line'
	1	/ 1x, a )
 120    format (1x,'Error: Group name parameter too long. Length= ',i4
	1	/1x, 'Parm= ', a)
 130    format (1x,'Error: Input group field name ',a,' is invalid')
 140    format (1x,'Error: No position number. Rest of line blank')
 150    format (1x,'Error: Input second parameter runs to end of line')
 160	format (1x,'Error: Parameter length too long. Length= ',i4) 
 425    format (1x,'Error: Bad return from OTS$CVT_TI_L, Status= ',z8.8)
 430    format (1x,'Error: Invalid percentage Number= ',i10)

!	Debug Formats:
 900    format(1x,'Parm= ',a, ' Parm length= ',i4,' Outstring start=',i4,       
	1	' pstart= ',i4,' pstop= ',i4)
 910    format(1x,'Parm= ',a, ' Parm length= ',i4,
	1	' pstart= ',i4,' pstop= ',i4)
 960    format(1x,'percentage number=',i8)

	FRD_Extract_idx_tolParms = fstatus

	Return
	End

