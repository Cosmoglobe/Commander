
	Integer*4 Function FRD_L_Process_Eng ( Englim_Rec, Lun, Lun2 )
! 
!	PROGRAM NAME : FRD_L_Process_Eng
! 
! 	PROGRAM DESCRIPTION:
! 	    This subroutine reads the engineering limits text files 
!	  (RMS disk files) and loads the red and yellow limits into record
!	  structured buffers. These limits will be used for checking 
!	  converted housekeeping data and setting related quality flags in 
!	  the science records. The red and yellow limiting values for each 
! 	  entity corresponding to one of 110 science record quality flags will
!         be compared to the values stored in the Engineering Record. The
!	  quality flags will be set in the Science record according to the
! 	  results of the comparisons. 
!           The FEX_ENGLIM contains all RDL substructures in Fex_ENG. The 
!	  limits for the engineering analogs will be obtained from the IGSE 
!	  database file, DB$:Firlims.DB, Which currently contains the red and 
!	  yellow limits for these analog database points. The red and yellow
!	  limits for the GRT low and high currents will be read from an edit
!	  file, FEX_Grtlim.Txt, by this program.
! 	    The binary engineering limits file produced by the program is 
!	  FEX_Englim.Dat. There will be one record with a structure dimensioned
!	  by four categories of limits: red low, yellow low, yellow high and
!         red high. 
!           
! 	AUTHOR:
! 	  Shirley M. Read
! 	  STX
! 	  January 1988
! 
! 	MODIFIED BY:
! 	  Shirley M. Read
! 	  STX
! 	  January 1988
!           The structure of the output engineering file has been modified
!           to facilitate access by the new COBETRIEVE routines in FDQ. The
!           file now has one record with a structure array for the four
!           categories of limits. New database words have been added.
! 
! 	INPUT PARAMETERS:
! 	  Englim_Rec -- engineering limits records
!	  Lun and Lun2 --- Unit numbers for input and output
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffers in ENGLIM_REC are filled.
! 
!       INPUT FILES:
! 	  DB$:Firlims.Inp
!         Grtlims.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Englim.Dat
! 
! 	INCLUDE FILES USED:
! 
! 	SUBROUTINES USED:
!	  Str$Upcase
! 	  FRD_L Functions:
! 		FRD_L_Extract_RParms
! 	        FRD_L_Store_Eng
!	        FRD_L_Extract_GParms
! 		FRD_L_Store_Grt
!	        FRD_L_Write_Eng
!	        FRD_L_Remtab
! 
! 	ERROR HANDLING:
! 	  Passed back in output parameter status from functions. 
! 
! 	METHOD USED:
! 
!	Description :   FRD_L_Process_Eng reads the STOL database text file,
!	              DB$Firlims.Inp containing engineering analog database 
!		      parameters and their Red and Yellow Limits into an 80
!	              character buffer. It calls FRD_L_Extract_RParms to
!		      extract the parameter and limits and convert the limits
!	              to real*4 values. It then calls FRD_L_Store_Eng to store
!		      the limits into an analog buffer. 
!		        FRD_L_Process_Eng then reads the Grtlims.Txt file
!		      to obtain GRT red and yellow limits. It calls 
!		      DQ_Extract_GParms to extract the GRT field name, number
!		      of the GRT within the filed array and the four limits and 
!		      convert the limits to real*4 values. The next call to 
!		      FRD_L_Store_GRT stores these limits into the corresponding
!		      GRT array position for output.
!		        Finally FRD_L_Process_Eng calls FRD_L_Write_Eng to 
!		      create and Write the binary FEX_Englim.Dat limits file
!		      for use in FDQ. 
!			For both the engineering analog and GRT limits, the
!	              text files are read into a string of 80 characters. 
!		      For the engineering text files, the first parameter
!	              is assumed to be the name of an engineering database 
!		      word. After the database word are the key words 'RANGE'
!		      and either 'ENGR' or 'COUNTS', meaning engineering 
!		      units or unconverted counts, respectively. Four real 
!		      limits, separated by blanks, must follow on the next 
!		      line. The order is Red Low, Yellow Low, Yellow High 
!		      and Red High. There are 64 analog limits to be checked.
!	              The GRT text file contains records consist of a GRT
!	              number followed by the red and yellow limits. There
!		      are 64 GRT limits to be checked.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters 

	Dictionary	'Fex_Englim'

	Record		/Fex_Englim/englim_rec

	Integer*4 lun, lun2     ! Logical unit numbers

!	Include files

	Include '(FUT_Convbuff)'
	Include '($SSdef)'

!	Functions

	Integer*4 FRD_L_Remtab
	Integer*4 FRD_L_Extract_RParms
	Integer*4 FRD_L_Extract_GParms
	Integer*4 FRD_L_Store_Eng
	Integer*4 FRD_L_Store_GRT
	Integer*4 FRD_L_Write_Eng
	Integer*4 Lib$Skpc	   ! Library function to skip characters
	Integer*4 Str$Upcase       ! Translate lower case letters to upper case

!	Local variables

	Character*80	Buffer 	   ! Input character string
	Character*20	Parm	   ! First parameter
	Integer*4	Lenparm    ! Length of first character parm
	Character*8 	GRTParm    ! GRT array field parameter
	Real*4	        RLimits(4) ! Real number limits for eng or grt file
				   ! Red Low, Yellow Low, Yellow High, Red High
	Logical*1       Name       ! Flag for eng and GRT name record or limits
	Logical*1       Engr       ! If true, engineering units : process
	Integer*2       Number     ! Number of the GRT
	Integer*4 	retstat	   ! Program Processing status
	Integer*4 	status 	   ! Returned system status
	Integer*4 	systatus   ! I/O status 
	Integer*4 	success / 1 /, error / 2 /
	Integer*4 	ipos1      ! Position on line
	Integer*4 	ix,jx      ! Index
	Logical*1 	eof        ! End of file
	Integer*4 	rec	   ! Index for records
	Logical*1       Store      ! Flag to store eng limits 
	Real*4    	Analogs(npoly,4)  ! Analog array: eng analogs,
					  ! 4 limits each
	Real*4	        GRTs(ngrt,4)      ! GRT array: GRTs, 4 limits each

!	Set status to success.

	retstat = success

! 	Open and read the DB$:Firlims.Inp file. Extract parameters from each
!	line. When EOF, close the file

	Open (unit=lun, file='DB$:Firlims.Inp',
	1 	status='old', form='formatted', access='sequential',
	2	organization='sequential', iostat=systatus,
	3	shared, readonly)
	If ( systatus .ne. 0 ) then
	    retstat = error
            Write(6,150) systatus
	Endif

!	Read the DB$:Firlims.Inp file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  name = .true.
	  engr = .false.
	  store = .false.
	  Do While ((retstat .eq. success) .and. (.not.eof))
            Read (unit=lun,iostat=systatus,end=500,err=400,fmt=300)
	1	  buffer
	    ipos1 = lib$skpc(' ',buffer)
	    If ( ipos1 .ne. 0 ) then	! Line is not blank
	      if ( buffer(1:1) .ne. '*'.and. buffer(1:1) .ne. '!' 
	1	 .and. buffer(1:6) .ne. 'LIMITS' 
	2	 .and. buffer(1:5) .ne. 'ENDLM' ) then  
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
	          status = FRD_L_Extract_RParms ( buffer, 
	1		parm, lenparm, rlimits, name, engr )
	          if ( status .eq. success ) then	      	       
	            if ( (engr) .and. (.not. name) ) then
		      store = .true.
	            elseif ( engr .and. name ) then
	              name = .false.
	            endif
	          else
	            retstat = error
	          endif
	        endif   ! Retstat is success
	        if (( retstat .eq. success) .and. (store)) then
	          status = FRD_L_Store_Eng ( parm, lenparm, rlimits, 
	1		 analogs )
	          if ( status .ne. success ) then
		    retstat = error
	          else
		    name = .true.
	            engr = .false.
	            store = .false.
	          endif
	        endif    ! Retstat is success and store
	      endif      ! Buffer line not a comment
	    Endif	 ! Buffer line is not blank
	  Enddo
400       If ( systatus .ne. 0 ) then
	    retstat = error
	    Write(6,200) status
	  Endif

!	Move the engineering analog limits into the Englim_Rec.
!	The first set of 62 analog limits are stored sequentially
!	in the Englim_Rec. There are no limits defined for the temperature
!       controller command bit readbacks; however the record definition 
!       still retains the fields. These fields will not have values.
!	The LMAC limits were added later to the record definition and must
!	be stored in the LMAC positions.

500       If ( retstat .eq. success ) then
	    Do ix = 1, 4
	      Call Lib$Movc3( 248, analogs(1,ix),
	1	Englim_Rec.lim(ix).en_analog.temp_ctrl )
	      Call Lib$Movc3( 8, analogs(63,ix),
	1	Englim_Rec.lim(ix).en_tail.lmac_analog_temp )
	    Enddo
	    Write(6,250) 
	  Endif
	Endif	! Retstat is success

!	Close DB$:Firlims.Inp file. If error, keep going.

	  Close(unit=lun,iostat=systatus)
          If (systatus .ne. 0 ) then
            Write (6,275) status
	  Endif

!       Engineering Analog formats
 150    format(1x,'Error: Failed to open DB$:Firlims.Inp. Status= ',z8.8)
 160    format(1x,'Error: Bad return from Str$Upcase. Status= ',z8.8)
 200    format(1x,'Error : Reading DB$:Firlims.Inp. Status= ',z8.8)
 250    format(1x,'Success : Reading of DB$:Firlims.Inp completed.')
 275    format(1x,'Error : Closing DB$:Firlims.Inp. Status= ',z8.8)
 300    format(a)

! 	Open and read the GRTLIM.TXT file. Extract parameters from each
!	line. When EOF, close the file

	Open (unit=lun, file='Fex_Grtlim.Txt',
	1 	status='old', form='formatted', access='sequential',
	2	organization='sequential', iostat=systatus)
	If ( systatus .ne. 0 ) then
	    retstat = error
            Write(6,1150) systatus
	Endif

!	Read the Fex_Grtlim.Txt file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  Do While ((retstat .eq. success) .and. (.not.eof))
            Read (unit=lun,iostat=systatus,end=1500,err=1400,fmt=1300)
	1	  buffer
	    ipos1 = lib$skpc(' ',buffer)
	    If ( ipos1 .ne. 0 ) then	! Line is not blank
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
	          status = FRD_L_Extract_GParms ( buffer, grtparm,
	1		   number, rlimits )
	          if ( status .eq. success ) then	      	
	            status = FRD_L_Store_Grt ( grtparm, number,
	1		     rlimits, grts )
	            if ( status .ne. success ) then
		      retstat = error
	            endif
	          endif
	        endif    ! Retstat is success
	      endif      ! Buffer line not comment
	    Endif	       ! Buffer line is not blank
	  Enddo
1400      If ( systatus .ne. 0 ) then
	    retstat = error
	    Write(6,1200) status
	  Endif

!	Move the GRT limits into the Englim_Rec. Skip the calibrator
!	resistor positions. Conversions and limits are not defined 
!       for the calibrator resistors.

1500      If ( retstat .eq. success ) then
	    Do ix = 1, 4
	      Call Lib$Movc3( 40, grts(1,ix),
	1	Englim_Rec.lim(ix).en_analog.A_lo_grt(1) )
	      Call Lib$Movc3( 8, grts(11,ix),
	1	Englim_Rec.lim(ix).en_analog.A_lo_grt(15) )
	      Call Lib$Movc3( 40, grts(13,ix),
	1	Englim_Rec.lim(ix).en_analog.A_hi_grt(1) )
	      Call Lib$Movc3( 8, grts(23,ix),
	1	Englim_Rec.lim(ix).en_analog.A_hi_grt(15) )
	      Call Lib$Movc3( 40, grts(25,ix),
	1	Englim_Rec.lim(ix).en_analog.B_lo_grt(1) )
	      Call Lib$Movc3( 8, grts(35,ix),
	1	Englim_Rec.lim(ix).en_analog.B_lo_grt(15) )
	      Call Lib$Movc3( 40, grts(37,ix),
	1	Englim_Rec.lim(ix).en_analog.B_hi_grt(1) )
	      Call Lib$Movc3( 8, grts(47,ix),
	1	Englim_Rec.lim(ix).en_analog.B_hi_grt(15) )
	    Enddo
	    Write(6,1250) 
	  Endif
	Endif	! Retstat is success
 
!	Close Fex_Grtlim.Txt file. If error, keep going.

	  Close(unit=lun,iostat=systatus)
          If (systatus .ne. 0 ) then
            Write (6,1275) status
	  Endif

!	Grt Formats
1150    format(1x,'Error: failed to open Fex_Grtlim.Txt. Status= ',z8.8)
1200    format(1x,'Error : Reading Fex_Grtlim.Txt. Status= ',z8.8)
1250    format(1x,'Success : Reading of Fex_Grtlim.Txt completed.')
1275    format(1x,'Error : Closing Fex_Grtlim.Txt. Status= ',z8.8)
1300    format(a)

!	If success, call routine to Write Fex_Englim.Dat.

	If ( retstat .eq. success ) then
            status = FRD_L_Write_Eng ( Englim_Rec, lun2 )
	    if ( status .eq. success ) then
	      Write(6,1675) 
	    else
	      retstat = error
	    endif
	Endif

!	Write Format 
1675    format(1x,'Success: Fex_Englim.Dat file written.')

	FRD_L_Process_Eng = Retstat

	Return
	End
