	Integer*4 Function FTB_IF_Midpoint_Time(Collect_frame,
	1					Transmit_frame,
	2					Scan_Length,
	3					Scan_Speed,
	4					Sweeps_per_IFG,
	5					Transmit_Time,
	6					IFG_time)
C
C******************************************************************************
C
C	 IF_Midpoint_Time
C
C
C FILE:  IF_Midpoint_Time.for
C
C
C FUNCTION: This routine computes the midpoint time of an interferogram.
C	It returns an odd number if successful and an even number if the
C	computation could not be made.
C	
C
C CREATED BY: S Hilinski			DATE: 2-4-85
C
C
C CALLING SEQUENCE:
C argument		type		input/output		description
C Collect_frame		I*4		input		frame at which the
C							interferogram collection
C							was begun.
C Transmit_frame	I*4		input		frame at which the
C							interferogram transmission
C							was begun.
C Scan_Length		I*2		input		1=long, 0=short.
C Scan_Speed		I*2		input		1=fast, 0=slow.
C Sweeps_per_IFG	I*2		input		Mirror sweeps in collection
C							cycle.
C Transmit_Time(2)	I*4		input		Spacecraft time associated
C							with Transmit_frame.
C IFG_Time(2)		I*4		output		Midpoint time.
C
C ERROR HANDLING:
C none
C
C ASSUMPTIONS:
C none
C
C******************************************************************************
C CHANGE LOG:
C	Version 4.2 11/15/88, SPR 2617, Shirley M. Read, STX
C		Multiple records with same time tags appeared in the FIRAS
C		science archive during the Observatory I and T. A copy of the 
C	        FIRAS Stripper's IF_Midpoint_Time subroutine has been modified
C	        to analyze the problem of invalid time tags. Th checks on the
C	        two frames have been commented out so that the routine can
C	        be run with some arguments in reverse to analyze the problem.
C
C******************************************************************************
C
C
	Implicit None
	Integer*4 Collect_frame		!Frame in which the IFG collection was
					!begun.
	Integer*4 Transmit_frame	!Frame in which the IFG was transmitted.
	Integer*2 Scan_Length
	Integer*2 Scan_Speed
	Integer*2 Sweeps_per_IFG
	Integer*4 Transmit_Time(2)	!Transmission time in VAX quadword format.
	Integer*4 IFG_time(2)		!Interferogram time in VAX quadword format.
c
	Integer*4 Delta_Time(2)		!frame difference X 1/4 seconds.
	Integer*4 Collect_Time(2)	!IFG collection start time.
	Integer*4 Half_Duration(2)	!1/2 IFG duration time.
	Integer*4 I
c
	Integer*4 Quarter_Second
	Parameter (Quarter_Second = 2500000)
c
	Integer*4 half_sweep_flyback(0:1,0:1) ! (sweep_time+flyback_time)/2
					  ! i=scan_length, j=scan_speed.
	data half_sweep_flyback/15500000, ! short-slow 3.1 sec.
	1			53500000, ! long-slow 10.7 sec.
	2			11750000, ! short-fast 2.35 sec.
	3			38500000/ ! long-fast 7.7 sec.
c
c
c
c Here we make sure that the transmit frame counter and collect frame
c counter do not contain screwy numbers.  (We don't want underflows or
c overflows.)
	IFG_Time(1)=Transmit_Time(1)
	IFG_Time(2)=Transmit_Time(2)

c Comment out the checks on the two frames so that the routine can be run 
c with some arguments in reverse for the analysis of the invalid times in 
c the science records.

!	if(Transmit_frame.lt.Collect_frame) then
!		IF_Midpoint_Time=0
!		go to 10000
!	else if(Transmit_Frame.lt.0) then
!		IF_Midpoint_Time=2
!		go to 10000
!	else if(Collect_Frame.lt.0) then
!		IF_Midpoint_Time=4
!		go to 10000
!	endif
c
c
c
c The goal here is to compute the delta time between Collect_frame and
c Transmit_frame in VAX quadword delta format.  Since the transmit frame always
c follows the collect frame, the difference will always be non-negative.
	Call lib$emul(Transmit_frame-Collect_frame,
	1 Quarter_second,0,Delta_Time)
c
c
c
c We must subtract the delta time from the transmission time to obtain the
c start of collection time.  This is all quadword arithmetic.
	Call lib$subx(Transmit_Time,Delta_Time,Collect_Time)
c
c
c
c We compute 1/2 the IFG duration.
	I=Sweeps_per_IFG
	Call lib$emul(half_sweep_flyback(Scan_Length,Scan_Speed),I,0,
	1 Half_Duration)
c
c
c
c The IFG midpoint is the collection time plus the half duration time.
	Call lib$addx(Collect_Time,Half_Duration,IFG_Time)
	FTB_IF_Midpoint_Time=1
c
10000	continue
	return
	end
