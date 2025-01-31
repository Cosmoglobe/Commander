      
	Integer*4 Function FRD_L_Process_Sci ( Scilim_Rec, Lun, Lun2 )
! 
!	PROGRAM NAME : FRD_L_Process_Sci
! 
! 	PROGRAM DESCRIPTION:
! 	    This subroutine reads the Science Limits text file (RMS disk file) 
! 	  and loads the red and yellow limits into record structured buffers.
! 	  These limits will be used for checking both science and attitude 
! 	  data and setting related quality flags in the science records. 
! 	  The red and yellow (when applicable) limiting values for each 
! 	  entity corresponding to one of 110 science record quality flags will
!         be compared to the values stored in the Raw Science Record or the 
! 	  Housekeeping Record. The quality flags will be set in the Science 
!	  record according to the results of the comparisons. 
! 	  The FEX_SCILIM is organized like the NFS_SDF, but omits the 
! 	  structures with fields that are not involved in the quality flag
!         checks. It contains all limits needed for comparison with data in  
! 	  the Science Record. There will be one record in the file with a 
!         structure dimensioned by two, one containing red limits and one 
!	  containing yellow limits. This file will be maintained offline by 
!	  the Firas Subsystem. 
!           FRD_L_Process_Sci will read an edit file, FEX_SciLim.Txt and 
!	  produce an RMS binary data file, FEX_SciLim.Dat, to be used by FDQ. 
!	  Eventually the file will be archived under COBETRIEVE. 
!           
! 	AUTHOR:
! 	  Shirley M. Read
! 	  STX
! 	  January 1988
! 
! 	MODIFIED BY:
!	  Shirley M. Read
!	  STX
!	  May 1988
!	    The science limits file was redesigned to facilitate the access 
!	    of the datasets via a new COBETRIEVE routine which updates all
!  	    buffers with the datasets whose time tags cover the data time.
!	    The file now contains only one record with the data structure
!	    dimensioned by two for red and yellow limits.
! 
! 	INPUT PARAMETERS:
! 	  NONE
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffers SCILIM_REC 
! 
!       INPUT FILES:
! 	  FEX_Scilim.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Scilim.Dat
! 
! 	INCLUDE FILES USED:
!        None
! 
! 	SUBROUTINES USED:
!	  Str$Upcase
!	  Lib$Skpc
!	  FRD_L functions:
!		FRD_L_Remtab
! 		FRD_L_Extract_IParms
! 	        FRD_L_Store_Sci
!	        FRD_L_Write_Sci
! 
! 	ERROR HANDLING:
! 	  Passed back in output parameter STATUS
! 
! 	METHOD USED:
! 
!	Description :   FRD_L_Define_Limits reads Fex_Scilim.Txt, containing 
!		      parameters and their Red and Yellow Limits, into an 80
!		      character buffer. It calls FRD_L_Extract_IParms to extract
!	              the parameter and limits from the buffer and convert the
!		      limits to integer values. 
!		      It then calls FRD_L_Store_Sci to store the limits in the
!		      corresponding field in the Scilim_Rec buffer. Finally, 
!	              it calls FRD_L_Write_Sci to create and Write the binary
!	              file, FEX_SciLim.Dat.
!		        For the science text files, the first parameter is
!		      the RDL field name of one of the science or attitude
!		      fields. After the science or attitude field name two 
!		      integer limits, separated by blanks, must follow on 
!		      the same line. The order of the limits is Red followed
!		      by Yellow. In some cases the red and yellow will be
!		      the same since the data quality flag has only a "red" 
!		      value. 
!
!----------------------------------------------------------------------

	Implicit None

!	Passed Parameters 

	Dictionary	'Fex_Scilim'

	Record		/Fex_Scilim/scilim_rec

	Integer*4 lun, lun2     ! Logical unit numbers

!	Include files

	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Skpc	   ! Library function to skip character
	Integer*4 FRD_L_Remtab
	Integer*4 Str$Upcase       ! Translate lower case letters to upper case
	Integer*4 FRD_L_Extract_IParms
	Integer*4 FRD_L_Store_Sci
	Integer*4 FRD_L_Write_Sci
	
!	Local variables

	Character*80	Buffer 	  ! Input character string
	Character*20	Parm	  ! First parameter
	Integer*4	Lenparm   ! Length of first character parm
	Integer*2	Limits(2) ! Red and yellow limits
	Integer*4       Ipos1     ! Position on line
	Integer*4 retstat	  ! Program Processing status
	Integer*4 status 	  ! Returned system status
	Integer*4 systatus        ! I/O status 
	Integer*4 success / 1 /, error / 2 /
	Logical*1 eof             ! End of file

!	Set status to success.

	retstat = success

! 	Open and read the Fex_Scilim.Txt file. Extract parameters from each
!	line. When EOF, close the file

	Open (unit=lun, file='Fex_Scilim.Txt',
	1 	status='old', form='formatted', access='sequential',
	2	organization='sequential', iostat=systatus)
	If ( systatus .ne. 0 ) then
	    retstat = error
            Write(6,150) systatus
	Endif

!	Read the Fex_Scilim.Txt file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  Do While ((.not.eof).and.(retstat .eq. success))
            Read (unit=lun,iostat=systatus,end=500,err=400,fmt=300)
	1	  buffer
	    ipos1 = Lib$Skpc(' ',buffer)
	    If ( ipos1 .ne. 0 ) then		! Not a blank line
	      if ( buffer(1:1) .ne. '*'.and. buffer(1:1) .ne. '!' ) then  
	        status = FRD_L_Remtab ( buffer )
	        if ( status .ne. success ) then
	          retstat = error
	        endif
	        if ( retstat .eq. success ) then
	          status = Str$Upcase (buffer, buffer)
	          if ( status .ne. success ) then
	            Write(6,160) status
	            retstat = error
	          endif
	        endif  ! Retstat is success
	        if ( retstat .eq. success ) then
	          status = FRD_L_Extract_IParms ( buffer, 
	1		parm, lenparm, limits ) 
	          if ( status .eq. success ) then
	            status = FRD_L_Store_Sci ( parm, lenparm, limits, 
	1		 scilim_rec )
	            if ( status .ne. success ) retstat = error
	          else
	            retstat = error
	          endif
	        endif   ! Retstat is success
	      endif  ! If buffer is not a comment
	    Endif    ! If buffer is not a blank line
	  Enddo      ! Read loop
 400      If ( systatus .ne. 0 ) then
	    retstat = error
	    Write(6,200) status
	  Endif
 500      If ( retstat .eq. success ) then
	    Write(6,250) 
	  Endif

!	Close Fex_Scilim.Txt. If error, keep going.

	  Close(unit=lun,iostat=systatus)
          If (systatus .ne. 0 ) then
            Write (6,275) status
	  Endif

	Endif  	! Retstat is success	  

!	If success, call routine to Write Fex_Scilim.Dat.

	If ( retstat .eq. success ) then
            status = FRD_L_Write_Sci ( Scilim_Rec, lun2 )
	    if ( status .eq. success ) then
	      Write(6,350) 
	    else
	      retstat = error
	    endif
	Endif

	FRD_L_Process_Sci = retstat

	Return

!	Formats
 150    format(1x,'Error: Failed to open Fex_Scilim.Txt. Status= ',z8.8)
 160    format(1x,'Error: Bad return from Str$Upcase. Status= ',z8.8)
 200    format(1x,'Error : Reading Fex_Scilim.Txt. Status= ',z8.8)
 250    format(1x,'Success : Reading of Fex_Scilim.Txt completed.')
 275    format(1x,'Error : Closing Fex_Scilim.Txt. Status= ',z8.8)
 300    format(a)
 350    format(1x,'Success: Fex_Scilim.Dat file written.')

	End
