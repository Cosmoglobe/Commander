
	Integer*4 Function FRD_L_Process_flags ( Limflags, Lun, Lun2 )
! 
!	PROGRAM NAME : FRD_L_Process_flags
! 
! 	PROGRAM DESCRIPTION:
! 	    This subroutine reads the Fex_Limflags text file (RMS disk file) 
! 	  and loads On/Off flags into a record structured buffer. Each record
!	  should contain the lim_flags array position and the value of 0 or -1
!	  ( false or true ). The routine will make an RMS output file 
!	  Fex_Limflags.Dat for use by FDQ.
!           
! 	AUTHOR:
! 	  Shirley M. Read
! 	  STX
! 	  January 1988
! 
! 	MODIFIED BY:
! 	  N/A
! 
! 	INPUT PARAMETERS:
! 	  NONE
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffer Limflags 
! 
!       INPUT FILES:
! 	  FEX_Limflags.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Limflags.Dat
! 
! 	INCLUDE FILES USED:
!        None
! 
! 	SUBROUTINES USED:
!	  Lib$Skpc
!	  FRD_L functions:
!		FRD_L_Remtab
! 		FRD_L_Extract_Flags
!		FRD_L_Store_Flags
!	        FRD_L_Write_Flags
! 
! 	ERROR HANDLING:
! 	  Passed back in output parameter Status
! 
! 	METHOD USED:
! 
!	Description :   FRD_L_Process_Flags reads FEX_Limflags.Txt, containing 
!		      the Lim_flags array position and the flag value, into 
!		      an 80 character buffer. It calls FRD_L_Extract_Flags to 
!		      extract the array position and flag value and store
!		      the flag into a logical array, Flags. When the reading
!		      of the file is complete, the flags are moved into the
!		      Fex_Limflags record.
!	              Finally it calls FRD_L_Write_Flags to create and write 
!		      a binary file, FEX_Limflags.Dat.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed Parameters 

	dictionary	'Fex_Limflags'

	record		/Fex_Limflags/limflags

	Integer*4 lun, lun2, frlun     ! Logical unit numbers

!	Functions

	Integer*4 FRD_L_Remtab
	Integer*4 FRD_L_Extract_Flags
	Integer*4 FRD_L_Store_Flags
	Integer*4 FRD_L_Write_Flags
	Integer*4 Lib$Skpc	     ! Library function to skip characters
	
!	Local variables

	Byte	        flagbuf(320) ! Array to match record limflags
	Character*80	buffer 	     ! Input character string
	Integer*4       ipos1        ! Position on line
	Character*15	flag         ! Limits flag name
	Integer*2       flaglen      ! Length of flag name
	Integer*2	number	     ! Number of values for flag array
	Byte		parm(20)     ! Flag values
	Integer*4 retstat 	     ! Program Processing status
	Integer*4 status 	     ! Returned system status
	Integer*4 systatus           ! I/O status 
	Logical*4 frstatus           ! CFR status
	Logical*4 FR_Init	     ! Field retriever utilities
	Logical*4 FR_Close_RDF
	Integer*4 success / 1 /, error / 2 /
	Logical*1 eof                ! End of file

!	Set status to success.

	retstat = success

! 	Open and read the Fex_Limflags.Txt file. Extract parameters from each
!	line. When EOF, close the file

	Open (unit=lun, file='Fex_Limflags.Txt',
	1 	status='old', form='formatted', access='sequential',
	2	organization='sequential', iostat=systatus)
	If ( systatus .ne. 0 ) then
	    retstat = error
            write(6,150) systatus
	Endif

!	Initialize the field retriever.

	frstatus = FR_Init ( frlun )
	If ( .Not. frstatus ) then
	    retstat = error
            write(6,180)
	Endif

!	Read the Fex_Limflags.Txt file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  do while (( .not. eof ) .and. (retstat .eq. success ))
            Read (unit=lun,iostat=systatus,end=500,err=400,fmt=300)
	1	  buffer
	    ipos1 = lib$skpc(' ',buffer)
	    if ( ipos1 .ne. 0 ) then		! Line is not blank
	      if ( buffer(1:1) .ne. '*'.and. buffer(1:1) .ne. '!' ) then 
	        status = FRD_L_Remtab ( buffer )
	        if ( status .ne. success ) then
	          retstat = error
	        endif
	        if ( retstat .eq. success ) then
	          status = FRD_L_Extract_Flags ( buffer, flag,
	1		   flaglen, number, parm ) 
	          if ( status .eq. success ) then
		    status = FRD_L_Store_Flags( frlun, flag, flaglen,
	1		     number, parm, flagbuf )
	            if ( status .ne. success ) then
	              retstat = error
                    endif
		  else
	            retstat = error
	          endif
	        endif     ! Retstat is success
	      endif	! Buffer is not a comment
	    endif	! Buffer is not blank
	  enddo
 400      If ( systatus .ne. 0 ) then
	    retstat = error
	    write(6,200) status
	  Endif
 500      If ( retstat .eq. success ) then
	    write(6,250) 
	  Endif
! 	Move the limit flags from the array to the Limflags Record.

	If ( retstat .eq. success ) then
	  Call Lib$Movc3( 256, flagbuf(65), 
	1	limflags.lim_flags.flg_badsci )
	Endif

!	Close Fex_Limflags.Txt. If error, keep going.

	  close(unit=lun,iostat=systatus)
          if (systatus .ne. 0 ) then
            write (6,275) status
	  endif

	Endif  	! Retstat is success	  

!	If success, call routine to write Fex_Limflags.Dat.

	If ( retstat .eq. success ) then
            status = FRD_L_Write_Flags ( Limflags, lun2 )
	    if ( status .eq. success ) then
	      write(6,350) 
	    else
	      retstat = error
	    endif
	Endif

!	Close up the field retriever.

	frstatus = FR_Close_RDF (0,0)

	If (.Not. frstatus) Then
            write (6,375)
	End If

!	Set function return status.

	FRD_L_Process_Flags = retstat

	Return

!	Formats
 150    format(1x,'Error: Failed to open Fex_Limflags.Txt. Status= ',z8.8)
 175	format(1x,'Error: Failed to get lun for field retriever. Status=',z8.8)
 180	format(1x,'Error: Failed to initialize field retriever.')
 200    format(1x,'Error : Reading Fex_Limflags.Txt. Status= ',z8.8)
 250    format(1x,'Success : Reading of Fex_Limflags.Txt completed.')
 275    format(1x,'Error : Closing Fex_Limflags.Txt. Status= ',z8.8)
 300    format(a)
 350    format(1x,'Success: Fex_Limflags.Dat file written.')
 375    format(1x,'Error : Failed closing field retriever.')

	end
