	Integer*4 Function FUT_Timerange ( start_time, bin_start_time,
     .					   stop_time, bin_stop_time )

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            May 7, 1987
C
C    INVOCATION:  STATUS = FUT_TIMERANGE ( START_TIME, BIN_START_TIME,
C					   STOP_TIME,  BIN_STOP_TIME )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	START_TIME		C*14	ASCII starttime in YYDDDHHMMSSCCC.
C	BIN_START_TIME(2)	I*4	Binary form of START_TIME.
C	STOP_TIME		C*14	ASCII stoptime in YYDDDHHMMSSCCC.
C	BIN_STOP_TIME(2)	I*4	Binary form of STOP_TIME.
C
C    SUBROUTINES CALLED: 
C	LIB$Get_Input
C	SYS$BinTim
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	CTUser.Inc
C	$SSDef
C
C----------------------------------------------------------------------
C
C Changes:
C
C	R. Kummerer, December 1, 1987. Use LIB$GET_INPUT to accept
C		input.  SPR 1750
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'CT$Library:CTUser.Inc'
	Include		'($SSDef)'

	Integer 	*2	ljstart
	Integer 	*2	ljstop
	Character 	*14  	jstart_default
	Character 	*14  	jstop_default

	Integer		*4	bin_start_time(2)
	Integer		*4	bin_stop_time(2)
	Integer		*4	current_time(2)
	Character	*14	start_time
	Character	*14	stop_time
	Integer		*4	status
	Logical		*1	ans_ok

	Integer		*4	LIB$Get_Input
	Integer		*4	SYS$BinTim
	Logical		*2	Time_LT
	Logical		*2	Time_GT

	External	FUT_Normal
	External 	FUT_InvTim

	Data jstart_default/'85001000000000'/
	Data jstop_default/'99365235959999'/

	ans_ok = .False.

	Do While (.Not. ans_ok)

C
C Accept a starttime.
C
	   Type 5
5	   Format ('                     YYDDDHHMMSSCCC')

	   status = LIB$Get_Input ( start_time, 'Enter a start time: ',
     .					ljstart )

	   If (status .Ne. SS$_Normal) Then
	     Call LIB$Signal ( %Val(status) )
	     FUT_Timerange = status
	     Return
	   End If

C
C Accept a stoptime.
C
	   Type 5

	   status = LIB$Get_Input ( stop_time, 'Enter a stop time:  ',
     .					ljstop )

	   If (status .Ne. SS$_Normal) Then
	     Call LIB$Signal ( %Val(status) )
	     FUT_Timerange = status
	     Return
	   End If

C
C Form the binary times.
C
	   start_time = start_time(1:ljstart) //
     .			     		jstart_default(ljstart+1:)
	   stop_time = stop_time(1:ljstop) //
     .			     		jstop_default(ljstop+1:)

	   Call CT_GMT_To_Binary(start_time,bin_start_time)
	   Call CT_GMT_To_Binary(stop_time,bin_stop_time)

	   If (Time_LT(bin_stop_time,bin_start_time)) Then
	      Call LIB$Signal(FUT_InvTim)
	   Else
	      ans_ok = .True.
	   End If

	End Do


	FUT_Timerange = %Loc(FUT_Normal)

	Return
	End
