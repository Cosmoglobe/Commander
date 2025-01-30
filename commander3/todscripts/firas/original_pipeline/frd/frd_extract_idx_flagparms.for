 
	Integer*4 Function FRD_Extract_IDX_flagparms ( Buffer , flag, report,
	1                                             rpt_lun, flag_text,index)
!
!	Program Name : FRD_Extract_Idx_flagsparms.For
!
!	Programmer  : Quoc C  chung, STX, April 14 1989.
!
!	Description : FRD_extract_idx_flagparms takes an input string of 80
!		      characters and extracts Three parameters, delimited
!	              by a blank ( or blanks ) from the string. The maximum
!		      length of first parameter is 30 characters. 
!		      The first parameter is the group field name for the 
!		      engineering record. It must have one of the 284
!		      Names.The second parameter is the index number for
!                     the group, i.e., a number from 1 to 284.  The last
!		      paramter is a FLag ( True / False).
!
!----------------------------------------------------------------------
!	Version 4.4   Quoc C Chung, STX, May 10 1989
!                     SPR 3773 Frd_Define_Idx_Params writes report to
!                     FOR000.dat.
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters

	Character*80	Buffer 	   ! Input character string.
        Integer*2       Index      ! converted index value
	Logical*1       Flag	   ! True/False value
        logical*1       report
        Character*30    Flag_text  ! individual item name
	
!	Include files

	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Skpc
	Integer*4 Lib$Locc
	Integer*4 Ots$Cvt_TI_L  ! Convert numeric text to integer
	Integer*4 Ots$Cvt_T_F	! Convert numeric text to floating
	
!	Local variables

	Integer*2       namelen    ! length of tols group field name
	Integer*4       lenparm    ! Length of number
	Integer*4       fstatus	   ! Function Processing status
	Integer*4       status 	   ! Returned system status
	Integer*4       success / 1 /, error / 2 /
	Integer*4       ipos1, ipos2 ! Relative positions on string
	Integer*4       pstart, pstop  ! Positions on string
        Integer*4       rpt_lun
!
	Character*1 Flag_string
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
	  if (report) write(rpt_lun,100) 
	endif
	if ( fstatus .eq. success ) then
	    pstart = ipos1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
	       if(report) write (rpt_lun,110) buffer 
 	    endif
	endif
 
	if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 30) then
	      if (report) write(rpt_lun,120) lenparm, buffer(pstart:pstop)
	      fstatus = error
	    else
	      Flag_text = buffer(pstart:pstop)
            endif   ! length
	endif     ! Fstatus is success

!	Get the second parm, the Index number.
	
	if ( fstatus .eq. success ) then
	  char_string3(1:3) = blank20(1:3)
	  pstart = pstop + 1
	  ipos1 = lib$skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    if (report) write(rpt_lun,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		if (report) write (rpt_lun,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 3) then
	      if (report) write(rpt_lun,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 4 - lenparm
	      char_string3(start:3) = buffer(pstart:pstop)
	      status = ots$cvt_ti_l(char_string3, Index,%val(2), maski)   
	      if ( status .eq. SS$_normal) then

	      else
	        if (report) write(rpt_lun,425) status
	        fstatus = error
	      endif  ! status
	    endif    ! lenparm
	  endif     ! status is Normal from parse of number on line

!	Check to see if Index number is valid, 1 to 284.

	  if ( fstatus .eq. success ) then
	    if( (index .lt. 1) .or. (index .gt. 284)) then
	       fstatus = error
	       if (report) write(rpt_lun,430) index
	    endif
          endif	    ! Fstatus is success on second parameter
        endif	    ! Fstatus is success after parse index on line

!	Get the Flag

	if ( fstatus .eq. success ) then
	  Flag_string(1:1) = blank20(1:1)
	  pstart = pstop + 1
	  ipos1 = lib$skpc ( ' ', buffer(pstart:80) )
	  if ( ipos1 .eq. 0 ) then
	    fstatus = error
	    if (report) write(rpt_lun,140) 
	  endif
	  if ( fstatus .eq. success ) then
	    pstart = pstart + ipos1 - 1
	    ipos2 = lib$locc ( ' ', buffer(pstart:80) )
	    if (ipos2 .eq. 0 ) then
		fstatus = error
		if (report) write (rpt_lun,150) 
 	    endif
	  endif
 
	  if (fstatus .eq. success ) then
	    pstop = pstart + ipos2 - 2
	    lenparm = pstop - pstart + 1
	    if (lenparm .gt. 1) then
	      if (report) write(rpt_lun,160) lenparm
	      fstatus = error
	    else	! Parm length is valid for a number
	      start = 2 - lenparm
	      Flag_string(start:1) = buffer(pstart:pstop)
	    endif   ! lenparm
	  endif     ! Fstatus is success 
        Endif       ! Fstatus is success after parse flag value
!	Check to see if Flag is valid True / False

	  if ( fstatus .eq. success ) then
            If ( Flag_string .eq. 'F' .or. Flag_string .eq. 'f' ) then
              Flag = .False.
            elseif (Flag_string .eq. 'T' .or. flag_string .eq. 't') then
              Flag = .True.
            Else
              if (report) Write(rpt_lun,600) Flag_string,index
              Fstatus = error
            endif   ! Flag_string
	  endif	    ! Fstatus is success
         


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
 600    Format (1x,' Either True or False Flag!!! something wrong ',l1,
	1       ' for this index = ',i4)

!	Debug Formats:

	FRD_Extract_idx_flagParms = fstatus

	Return
	End

