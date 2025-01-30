
	Integer*4 Function Frd_parse_idx_flags(flag_rec,Lun,report,rpt_lun,
	1                                      newlun)
C
C       This routine reads the Fex_idx_flag.txt (RMS file), and
C       get the idx name and values of each record, use the values  
C       to update the Frd_idx_flags.dat .
! 
!	PROGRAM NAME : FRD_Parse_idx_flags
! 
!           
! 	AUTHOR:
! 	  Quoc C Chung
! 	  STX
! 	  April 14 1989
! 
! 	MODIFIED BY:
! 	  N/A
! 
! 	INPUT PARAMETERS:
! 	  NONE
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffer Flags
! 
!       INPUT FILES:
! 	  FEX_Idx_flag.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Idx_flag.Dat
! 
! 	INCLUDE FILES USED:
!        None
! 
! 	SUBROUTINES USED:
!	  Lib$Skpc
!	  FRD_L functions:
!		FRD_L_Remtab
! 		FRD_Extract_Idx_flagsparm
! 
! 	ERROR HANDLING:
! 	  Passed back in output parameter Status
! 
! 	METHOD USED:
! 
!	Description : FRD_Parse_idx_flags reads FEX_idx_flag.Txt, containing 
!		      the Field name ,index number and the T/F Value, into 
!		      an 80 character buffer. It calls FRD_Extract_idx 
!		      _flagsparm to extract the field name , index number ,
!                     and the T/F value and store value in a record array.
!                     when the reading of the text file is completed. The
!                     Values are moved into the Fex_idx_flag record, finally,
!                     write the record to a binary file (RMS file),
!                     Fex_idx_flag.dat.
!
!----------------------------------------------------------------------
!       VERSION 4.4  QUOC C CHUNG, STX, MAY 09, 1989
!       SPR 3788  FRD_DEFINE_IDX_PARAMS CAN NOT BE INTERACTIVE.
!----------------------------------------------------------------------       
!       VERSION 4.4  QUOC C CHUNG, STX, MAY 09, 1989
!       SPR 3773  FRD_DEFINE_IDX_PARAMS WRITES REPORT TO FOR000.DAT.
!----------------------------------------------------------------------
	Implicit None
!
!	Include Files

	Include '($SSdef)'
        Include '(Fut_text_item)'
!	
!	Pass parameter

	dictionary 'Fex_idx_flag'

        Record    /Fex_Idx_Flag/Flag_rec

        Integer*4 lun,rpt_lun,newlun
        
!	Functions

	Integer*4 FRD_L_Remtab
	Integer*4 FRD_Extract_idx_flagparms
	Integer*4 Lib$Skpc	     ! Library function to skip characters
	
!	Local variables

	Character*80	buffer 	      ! Input character string
        Character*30    item_text(284)! Group field names
        Character*30    Flag_text     ! individual field name from buffer
        Character       Ans           ! user response
!
	Integer*4 lpos               ! Position on line
	Integer*4 retstat 	     ! Program Processing status
	Integer*4 status 	     ! Returned system status
	Integer*4 systatus           ! I/O status 
        Integer*4 Ix
	Integer*4 success / 1 /, error / 2 /
!
        Integer*2 index              ! converted index number
        Integer*2 index_no(284)      ! Group index number
!
	Logical*1 eof                ! End of file
        Logical*1 idx_flags(284)/284*.true./
        Logical*1 Flags/.False./,Flag
        Logical*1 report

!	Set status to success.

	retstat = success



!	Read the Fex_idx_flag.Txt file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  do while (( .not. eof ) .and. (retstat .eq. success ))
            Read (unit=lun,iostat=systatus,end=500,err=4000,fmt=300)
	1	  buffer
	    lpos = lib$skpc(' ',buffer)
	    if ( lpos .ne. 0 ) then		! Line is not blank
	      if ( buffer(1:1) .ne. '*'.and. buffer(1:1) .ne. '!' ) then 
	        status = FRD_l_Remtab ( buffer )
	        if ( status .ne. success ) then
	          retstat = error
	        endif
	        if ( retstat .eq. success ) then
	          status = FRD_Extract_idx_flagparms ( buffer, Flag,report,
	1                                      rpt_lun,flag_text, index)
	          if ( status .eq. success ) then
                     If ( index .gt. 0 .or. index .le. 284) then
                       item_text(index) = Flag_text
                       index_no(index) = index                  
                       idx_flags(index) = Flag
                     else
                      retstat = error
                      if (report) write (rpt_lun,10) index 
                     endif  ! index
		  else
	            retstat = error
	          endif  ! Status 
	        endif    ! Retstat is success
	      endif	! Buffer is not a comment
	    endif	! Buffer is not blank
	  enddo
       endif    ! Retstat is success

 4000     If ( systatus .ne. 0 ) then
	    retstat = error
	    if (report)write(rpt_lun,200) status
	  Endif
 500      If ( retstat .eq. success ) then
	    if ( report) write(rpt_lun,250) 
	  Endif
!
! Display the Fex_Idx_flag.txt for verification
!
          write(6,61)
          if (report) write(rpt_lun,61)
          Do ix = 1 , 284
            write(6,6100) item_text(ix),index_no(ix),idx_flags(ix)
            if (report)
	1   write(rpt_lun,6100) item_text(ix),index_no(ix),idx_flags(ix)
          Enddo
 61     Format(10x,'Item Name ',18x,'Index#',4x,'Flag'/80('-'))
 6100	Format(4x,a,2x,i4,10x,l1)

!       Open the output file, Fex_Idx_flag.Dat 

	  if ( retstat .eq. success ) then

	    Open (unit=newlun, file='Fex_Idx_flag.Dat',
	1 	  status='new', form='unformatted', 
	2	  organization='sequential', iostat=status,
	3	  recordsize=128, recordtype='fixed')

	    if (status .ne. 0) then
	      Write(Rpt_lun,100) status
	      retstat = error
	    endif


! 	Move the idx flags from the array to the idx flags Record.

	If ( retstat .eq. success ) then
	  do ix = 1 , 284
	     flag_rec.idx_flags.flags(ix) = idx_flags(ix)
          enddo
	Endif

!	Close Fex_idx_flag.Txt. If error, keep going.

	  close(unit=lun,iostat=systatus)
          if (systatus .ne. 0 ) then
            write (6,275) status
            if (report) write (rpt_lun,275) status
	  endif

	Endif  	! Retstat is success	  

!	If success,  write Fex_idx_flag.Dat.

	If ( retstat .eq. success ) then
             write(unit= newlun,iostat=systatus) flag_rec
	    if ( systatus .ne. 0 ) then
	      if (report) write(rpt_lun,350) 
	      retstat = error
	    else
              if (report) Write(rpt_lun,400)
	    endif
	Endif


!	Set function return status.

	FRD_Parse_idx_flags = retstat

	Return

!	Formats
!
  10	format(1x,' Invalid class number: ',i4)
 100	format(1x,'Error: Opening Fex_Idx_flag.Dat. Status= ',z8.8)
 200    format(1x,'Error : Reading Fex_idx_flag.Txt. Status= ',z8.8)
 250    format(1x,'Success : Reading of Fex_idx_flag.Txt completed.')
 275    format(1x,'Error : Closing Fex_idx_flag.Txt. Status= ',z8.8)
 300    format(a)
 350    format(1x,'fail: Fex_idx_flag.Dat file written.')
 400    format(1x,'success: Fex_idx_flag.Dat file written.')
	end
