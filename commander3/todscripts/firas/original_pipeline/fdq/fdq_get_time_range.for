	INTEGER*4 FUNCTION FDQ_GET_TIME_RANGE ( save_init, save_fin,
	1		more_segments, time_range )
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_TIME_RANGE
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine will compare the initial and final times for all
C/	  science channels which have catalog entries and the corresponding
C/	  housekeeping file and set the time range to the entire data range.
C/
C/	AUTHOR:
C/	  Shirley M. Read
C/	  STX
C/	  July 20, 1988
C/
C/	MODIFICATIONS:
CH
CH	  Version 4.4.1 08/20/89 SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments related changes.
CH
CH	  Version 4.5 NOV/27/89 SPR 5165 , H. WANG, STX
CH		FAIL TO PROCESS WHENEVER MISSING CHANNEL SEGMENTS
CH              ENCOUNTER
C/
C/	CALLING SEQUENCE:
C/	  Return_Status = FDQ_GET_TIME_RANGE ( SAVE_INIT, SAVE_FIN,
C/		MORE_SEGMENTS, TIME_RANGE )
C/
C/	INPUT PARAMETERS:
C/
C/	  SAVE_INIT(2,5) I*4    Intial data times
C/	  SAVE_FIN(2,5)  I*4    Final data times
C/	  MORE_SEGMENTS  L*1    More segments available for science channels
C/
C/	OUTPUT PARAMETERS:
C/
C/	  TIME_RANGE       C*30 Time range for data
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE 
C/
C/	INCLUDE FILE USED:
C/	  $SSDEF
C/
C/	SUBROUTINES CALLED:
C/
C/	  TIME_LT
C/	  TIME_GT
C/	  CT_BINARY_TO_GMT
C/
C/	METHOD USED:
C/	  < PDL for this routine > :
C/	  SET the return status to a good value ( 1 ).
C/	  SET the min and max times to the first array position with times.
C/	  DO for each other initial time and final time
C/	    If the current initial time is less than the reference min,
C/	      RESET the min
C/	    If the current final time is greater than the reference max,
C/	      RESET the max
C/	  ENDDO
C/	  CALL CT_BINARY_TO_GMT to convert the data min time to ACSII
C/	       and load into the Time_Range
C/	  CALL CT_BINARY_TO_GMT to convert the data max time to ACSII
C/	       and load into the Time_Range
C/	  SET the function value to success for return.

	IMPLICIT	NONE

!	External Parameters and Include Files

	INCLUDE		'($SSDEF)'

	EXTERNAL	FDQ_NORMAL
	EXTERNAL        FDQ_ERSUBX
	EXTERNAL        FDQ_ERADDX
	EXTERNAL        FDQ_ABERR

!	Passed Parameters

	integer*4	save_init(2,5)          ! Initial data times
	integer*4	save_fin(2,5)           ! Final data times
	logical*1       more_segments(4)        ! More channel data available
	character*30    time_range              ! Data time range 

!	Functions

	logical*1	Time_Lt		        ! Ct time less than 
	logical*1	Time_Gt		        ! Ct time greater than 
	integer*4       LIB$Subx        	! Extended precision subtract
	integer*4       LIB$Addx        	! Extended precision add

!	Local Declarations

	integer*4 	status, success /1/, error /2/ ! Status of processing
	integer*4 	retstat
	integer*2	min, max		! Min and max index pointers
	integer*2       ix, jx			! Indices
	logical*1       found                   ! Found initial min and max
	character*14    chartime		! Character time string
	integer*4       sec_64(2)               ! Offset from start telemetry
	data            sec_64(1) / 640000000 / ! time for collect time
	integer*4       len / 2 / 		! Lib$Subx length
	integer*4	ub

!	Set return status

	status = success

!	Check min and max times. This subroutine is not called unless one
!	channel has data.

	found = .false.
	do ix = 1, 4
	  if ((more_segments(ix)) .and. (.not. found)) then
	    min = ix
	    max = ix
	    found = .true.
	  endif
	enddo
	jx = min
	do ix = jx, 5
         If(more_segments(ix)) then  
	  if (Time_LT(save_init(1,ix),save_init(1,min))) min = ix
	  if (Time_GT(save_fin(1,ix),save_fin(1,max))) max = ix
         endif
	enddo

!	Set the output time range from min and max.	

!	Compute adjusted binary start and stop times for housekeeping and
!	config file ranges.  Subtract 64 seconds from the earliest start time
!	and add 64 seconds to the latest stop time.
        Chartime = time_range(1:14)
	call CT_gmt_to_binary ( chartime, save_init(1,min) )
	retstat = Lib$Subx (save_init(1,min),sec_64,save_init(1,min),len)
	retstat = Lib$Subx (save_init(1,min),sec_64,save_init(1,min),len)
	if (retstat .ne. ss$_normal) then
	  status = error
	  call Lib$Signal(FDQ_ERSUBX, %val(1),%val(retstat))
	endif

	retstat = Lib$Addx (save_fin(1,max),sec_64,save_fin(1,max),len)
	if (retstat .ne. ss$_normal) then
	  status = error
	  call Lib$Signal(FDQ_ERADDX, %val(1),%val(retstat))
	endif

	call CT_Binary_to_GMT ( save_init(1,min), chartime )
	time_range(1:14) = chartime(1:14)
	time_range(15:15) = ';'
	call CT_Binary_to_GMT ( save_fin(1,max), chartime )
	time_range(16:29) = chartime(1:14)
	time_range(30:30) = ';'

!	Set function to return status

	if (status .eq. success) then
	  FDQ_GET_TIME_RANGE = %loc(FDQ_NORMAL)
	else
	  FDQ_GET_TIME_RANGE = %loc(FDQ_ABERR)
	endif
	return
	end
