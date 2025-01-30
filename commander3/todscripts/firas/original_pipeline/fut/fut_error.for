      Integer*4 Function FUT_Error ( sigargs, mechargs )

C------------------------------------------------------------------------
C    PURPOSE: Flag and log errors.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STI Incorporated
C            September 18, 1986
C
C	     Adapted from worked performed by E.Cheng.
C
C    INVOCATION: Call LIB$Establish ( FUT_ERROR )
C
C    INPUT PARAMETERS:
C      SIGARGS(*)		I*4	Signal array
C      MECHARGS(*)		I*4	Mechanism array
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED:  LIB$SIGNAL
C
C    COMMON VARIABLES USED:
C      ERROR(100)              I*4      List of error messages.
C      ERROR_COUNT             I*2      Number of errors recorded.
C
C    INCLUDE FILES:  FUT_ERROR.TXT
C
C----------------------------------------------------------------------
C
C Changes:
C
C	SPR 1855, Trap system errors such as divide by zero. R. Kummerer,
C		January 12, 1987.
C
C	SPR 2529, Convert severities worse than WARNING to WARNING so that
C		tracebacks are available, especially on integer overflows
C		and the like.  R. Kummerer, April 7, 1989.
C
C       SPR 2805, Add statement to spool errors to FIRAS report files
C               via SYS$PUTMSG system command in the FUT_ERROR. This
C               subroutine calls FUT_WRITE_ERROR to write error messages
C               into a facility report file.   
C                   Nilo G. Gonzales/STX, Dec. 20, 1990.
C----------------------------------------------------------------------

	Implicit	None

	Include 	'($SSDEF)'
	Include 	'($MTHDEF)'
	Include 	'($STSDEF)'
	Include		'(FUT_Error)'			!Error tracking

	Integer		*4	sigargs(*)		!Signal array
	Integer		*4	mechargs(*)		!Mechanism array
	Integer         *4      newsigargs(10)
	Integer         *4      element
	Integer         *2      I
	External		fut_write_error
	Integer         *4      fut_write_error
	Integer		*4	imatch
	Integer		*4	severity
	Real		*8 	log_flag
	Real		*8 	square_root_flag
	Integer		*4 	integer_flag
        Real		*8  	invalid_arg_flag

	Parameter (log_flag = -9999.0)
	Parameter (square_root_flag = -9999.0)
	Parameter (integer_flag = -9999)
        Parameter (invalid_arg_flag = -9999.0)

	Integer		*4	LIB$Match_Cond
	Integer		*4	LIB$Decode_Fault
	Integer		*4	LIB$Fixup_Flt

	External		FUT_Fix_FltDiv_F
	External		FUT_Fix_FltOvf_F
C
C Determine if this is an error of interest
C
	IMATCH = LIB$MATCH_COND(SIGARGS(2)
     1		,ss$_unwind
     2		,ss$_intdiv
     3		,mth$_logzerneg
     4		,ss$_fltdiv_f
     5		,ss$_roprand
     6		,mth$_squrooneg
     7		,ss$_intovf
     8		,mth$_invargmat
     9		,ss$_fltovf_f)

C
C This error is none of our business
C
C Get rid of last two elements in SIGARGS (the PC and PSL)
C then pad NEWSIGARGS with zeros.

	ELEMENT = 1
	NEWSIGARGS (ELEMENT) = 0

	DO I = 1, SIGARGS(1) - 2
	   ELEMENT = ELEMENT + 1
	   NEWSIGARGS (ELEMENT) = SIGARGS (ELEMENT)
	END DO
	DO I = ELEMENT + 1, 10
	   ELEMENT = ELEMENT + 1
	   NEWSIGARGS (ELEMENT) = 0
	END DO
	if (IMATCH .eq. 0) then
		call MVBITS (SIGARGS (2), 0, 3, SEVERITY, 0)
		if (SEVERITY .EQ. STS$K_ERROR .OR.
     1		    SEVERITY .EQ. STS$K_SEVERE) then
		  call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
		end if
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal
C
C Ignore unwinding calls
C
	else if (IMATCH .eq. 1) then
		FUT_Error = ss$_continue
C
C Convert integer divide by zeros into -9999 (or -255 for bytes)
C
	else if (IMATCH .eq. 2) then
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal

C
C Convert LOG(0) and LOG(-n) to -9999.0
C
	else if (IMATCH .eq. 3) then
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal

		call LIB$MOVC3 (4, log_flag, mechargs (4))
C
C We have encountered a floating divide by zero fault.  In order to
C proceed properly, we must convert the fault to a trap.  Note that
C 'LIB$SIM_TRAP' does not return if the signaled error is indeed a 
C fault and it can do something about it (i.e. complete the instruction
C and then simulate a trap).
C
C *** This method of handling the fault leaves a reserved operand in  ***
C *** as the result.  This may lead to reserved operand faults in the ***
C *** future if the caller does not know how to handle them.          ***
C
C	else if (IMATCH .eq. 4) then
C		FUT_Error = LIB$SIM_TRAP(SIGARGS,MECHARGS)
C
C Here is a more involved solution which writes a floating zero into the 
C output operand.  We call 'LIB$DECODE_FAULT' to actually muck with the
C divide instruction which generated the fault and write a zero into the 
C output operand.  We should never return from this call, so the action
C routine must signal a floating divide by zero trap so we output a message.
C
	else if (IMATCH .eq. 4) then
		FUT_Error = LIB$DECODE_FAULT(SIGARGS,MECHARGS
     1					,%descr(FUT_FIX_FLTDIV_F))
C
C Convert reserved operands into zeroes
C
	else if (IMATCH .eq. 5) then
		FUT_Error = LIB$FIXUP_FLT(SIGARGS,MECHARGS)
C
C Convert negative square roots to -9999.0
C
	else if (IMATCH .eq. 6) then
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal

		call LIB$MOVC3 (4, square_root_flag, mechargs (4))
C
C Convert integer overflows into the maximum integer value
C
	else if (IMATCH .eq. 7) then
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal

C
C Convert invalid arguments to  -9999.0
C
	else if (IMATCH .eq. 8) then
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal

		call LIB$MOVC3 (4, invalid_arg_flag, mechargs (4))
C
C Convert floating overflows to -9999.0
C
	else if (IMATCH .eq. 9) then
		FUT_Error = LIB$DECODE_FAULT(SIGARGS,MECHARGS
     1					,%descr(FUT_FIX_FLTOVF_F))
C
C Here we convert the error into a gentle informational message
C
	else
		call MVBITS (STS$K_WARNING, 0, 3, SIGARGS (2), 0)
	CALL SYS$PUTMSG (NEWSIGARGS, FUT_WRITE_ERROR,)
		FUT_Error = ss$_resignal
	endif
C
C Record the error.
C
	error_count = error_count + 1

	If (error_count .Gt. 100) Then
	  error_count = 100
	End If

	error(error_count) = sigargs(2)
C
C Done
C
	return
	end

	integer*4 function FUT_FIX_FLTDIV_F(OPCODE,INSTR_PC,PSL,REGISTERS
	1			,OP_COUNT,OP_TYPES,READ_OPS,WRITE_OPS
	2			,SIGARGS,SIGNAL_ROUT,CONTEXT
	3			,USER_ARG,ORIG_REGS)
C
	integer*4 OPCODE,INSTR_PC,PSL
	integer*4 REGISTERS(0:15)
	integer*4 OP_COUNT
	integer*4 OP_TYPES(*)
	integer*4 READ_OPS(*)
	integer*4 WRITE_OPS(*)
	integer*4 SIGARGS(*)
	integer*4 SIGNAL_ROUT
	integer*4 CONTEXT
	integer*4 USER_ARG
	integer*4 ORIG_REGS(0:15)
C
	external SIGNAL_ROUT
C
C Define codes
C
	include '($SSDEF)'
	include '($PSLDEF)'
	include '($LIBDCFDEF)'
C
C Local variables
C
	integer*4 OP_DTYPE
	integer*4 OP_SIZE
	integer*4 SIGNAL_ARRAY(2)
	real*8 FLTDIV_FLAG
	parameter (FLTDIV_FLAG=-9999.0)
C
C Check for FPD and generate a reserved operand signal if true ('SIGNAL_ROUT' 
C never returns)
C
	if (btest(PSL,psl$v_fpd)) then
		SIGNAL_ARRAY(1) = 1
		SIGNAL_ARRAY(2) = ss$_roprand
		call SIGNAL_ROUT(1,CONTEXT,SIGNAL_ARRAY)
	endif
C
C Get data type of the instruction
C
	OP_DTYPE = ibits(OP_TYPES(2),lib$v_dcftyp,lib$s_dcftyp)
C
C Move the appropriate 'zero' into the result operand (always the last
C operand)
C
C F_floating
C
	if (OP_DTYPE .eq. LIB$K_DCFTYP_F) then
		call LIB$MOVC3(4,FLTDIV_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C D_floating
C
	elseif (OP_DTYPE .eq. LIB$K_DCFTYP_D) then
		call LIB$MOVC3(8,FLTDIV_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C G_floating
C
	elseif (OP_DTYPE .eq. LIB$K_DCFTYP_G) then
		call LIB$MOVC3(8,FLTDIV_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C This should never happen
C
	else
		stop '*** FUT_FIX_FLTDIV_F Internal Error - undefined OP_DTYPE ***'
	endif
C
C Fixup PSL bits to indicate a zero result
C
	PSL = ibclr(PSL,psl$v_n)
	PSL = ibset(PSL,psl$v_z)
	PSL = ibclr(PSL,psl$v_v)
	PSL = ibclr(PSL,psl$v_c)
C
C Done with fixup, resignal <ss$_fltdiv>
C
	SIGNAL_ARRAY(1) = 1
	SIGNAL_ARRAY(2) = ss$_fltdiv
	call SIGNAL_ROUT(0,CONTEXT,SIGNAL_ARRAY)
C
C The previous call should never return
C
	stop '*** FUT_FIX_FLTDIV_F Internal Error - SIGNAL_ROUT returned ***'
	end

	integer*4 function FUT_FIX_FLTOVF_F(OPCODE,INSTR_PC,PSL,REGISTERS
	1			,OP_COUNT,OP_TYPES,READ_OPS,WRITE_OPS
	2			,SIGARGS,SIGNAL_ROUT,CONTEXT
	3			,USER_ARG,ORIG_REGS)
C
	integer*4 OPCODE,INSTR_PC,PSL
	integer*4 REGISTERS(0:15)
	integer*4 OP_COUNT
	integer*4 OP_TYPES(*)
	integer*4 READ_OPS(*)
	integer*4 WRITE_OPS(*)
	integer*4 SIGARGS(*)
	integer*4 SIGNAL_ROUT
	integer*4 CONTEXT
	integer*4 USER_ARG
	integer*4 ORIG_REGS(0:15)
C
	external SIGNAL_ROUT
C
C Define codes
C
	include '($SSDEF)'
	include '($PSLDEF)'
	include '($LIBDCFDEF)'
C
C Local variables
C
	integer*4 OP_DTYPE
	integer*4 OP_SIZE
	integer*4 SIGNAL_ARRAY(2)
	real*8 FLTOVF_FLAG
	parameter (FLTOVF_FLAG=-9999.0)
C
C Check for FPD and generate a reserved operand signal if true ('SIGNAL_ROUT' 
C never returns)
C
	if (btest(PSL,psl$v_fpd)) then
		SIGNAL_ARRAY(1) = 1
		SIGNAL_ARRAY(2) = ss$_roprand
		call SIGNAL_ROUT(1,CONTEXT,SIGNAL_ARRAY)
	endif
C
C Get data type of the instruction
C
	OP_DTYPE = ibits(OP_TYPES(2),lib$v_dcftyp,lib$s_dcftyp)
C
C Move the appropriate 'zero' into the result operand (always the last
C operand)

C
C F_floating
C
	if (OP_DTYPE .eq. LIB$K_DCFTYP_F) then
		call LIB$MOVC3(4,FLTOVF_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C D_floating
C
	elseif (OP_DTYPE .eq. LIB$K_DCFTYP_D) then
		call LIB$MOVC3(8,FLTOVF_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C G_floating
C
	elseif (OP_DTYPE .eq. LIB$K_DCFTYP_G) then
		call LIB$MOVC3(8,FLTOVF_FLAG,%val(WRITE_OPS(OP_COUNT)))
C
C This should never happen
C
	else
		stop '*** FUT_FIX_FLTOVF_F Internal Error - undefined OP_DTYPE ***'
	endif
C
C Fixup PSL bits to indicate a negative result
C
	PSL = ibclr(PSL,psl$v_n)
	PSL = ibset(PSL,psl$v_z)
	PSL = ibclr(PSL,psl$v_v)
	PSL = ibclr(PSL,psl$v_c)
C
C Done with fixup, resignal <ss$_fltovf>
C
	SIGNAL_ARRAY(1) = 1
	SIGNAL_ARRAY(2) = ss$_fltovf
	call SIGNAL_ROUT(0,CONTEXT,SIGNAL_ARRAY)
C
C The previous call should never return
C
	stop '*** FUT_FIX_FLTOVF_F Internal Error - SIGNAL_ROUT returned ***'
	end
