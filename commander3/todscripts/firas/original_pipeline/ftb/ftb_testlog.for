
	Program FTB_TestLog

C------------------------------------------------------------------------
C    PURPOSE:	Select information from the test conductor's log based
C		on test type, test number, start and stop times, and a
C		keyword from the test title.  The user is then given
C		the option of writing the selected start and stop times
C		to an RMS file for input into the automated pipeline.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 1, 1987
C
C    INVOCATION:  TESTLOG
C
C    INPUT PARAMETERS:  Query input select.
C
C    OUTPUT PARAMETERS:  User selected RMS file containing JSTART, JSTOP
C			 times.
C
C    SUBROUTINES CALLED: 
C
C		Time_Lt
C		STR$Compare
C		STR$UpCase
C		STR$Trim
C		DTR$Init
C		DTR$Command
C		DTR$Continue
C		DTR$DTR
C		DTR$Finish
C
C    COMMON VARIABLES USED:  None
C
C    INCLUDE FILES: 
C
C		$SSDEF
C		DTR$Library:DAB
C
C--------------------------------------------------------------------------
C
C Changes:
C
C	Use FUT_TIMERANGE to enter a timerange. R. Kummerer, July 24, 1987.
C
C	SPR 1542, Select multiple items. R. Kummerer, January 4, 1988.
C
C       Shirley M. Read, STX, August 26, 1988
C	   Reason: SPRs 1566 and 1763 requested that all FIRAS facilities 
C          should tell of successful completion and signal errors via
C	   $status for batch processing. Also added interface with Fut_Error,
C	   the condition handler.
C
C--------------------------------------------------------------------------

	Implicit None

	Include		'($SSDEF)'
	Include		'DTR$Library:DAB'

C	FTB messages

	External   	ftb_normal
	External	ftb_aberr

C	Condition Handler

	External        Fut_Error

C	Local Variables

	Integer		*4	status
	Character	*255	cmd
	Character	*64	test_type
	Character	*64	test_num
	Character	*64     test_title
	Integer		*4	title_len
	Character	*8	ttitle_len
	Character	*20	ans
	Logical		*1	ans_ok
	Character	*64	qfile
	Logical		*1	fetch
	Integer		*2	lcmd
	Integer		*2	k
	Integer		*2	l
	Character	*1	BL/' '/
	Character	*1	TAB/'	'/

	Integer		*4	start_time(2)
	Integer		*4	end_time(2)
	Character	*14	starting_time
	Character	*14	ending_time

	Integer		*4	FUT_Timerange
	Logical		*2	Time_Lt
	Integer		*4	STR$Compare
	Integer		*4	STR$UpCase
	Integer		*4	STR$Trim
	Integer		*4	STR$Find_First_Not_In_Set
	Integer		*4	STR$Find_First_In_Set
	Integer		*4	STR$Translate
	Integer		*4	DTR$Init
	Integer		*4	DTR$Command
	Integer		*4	DTR$Continue
	Integer		*4	DTR$DTR
	Integer		*4	DTR$Finish

c     Establish condition handler.            

	call lib$establish ( fut_error )
c
c Inititalize.
c
	qfile = '*.''a filename'''

	status = DTR$Init ( DAB, 100, msg_buff, aux_buff,
     .				DTR$K_Semi_Colon_Opt )

	If (status .Eq. SS$_Normal) Then

c
c		Ready the test log domain.
c
	   cmd = 'READY TLDFI SHARED READ;'
	   status = DTR$Command ( DAB, cmd )
	   status = DTR$Continue ( DAB )

	   If (status .Eq. SS$_Normal) Then

c
c		Fetch records from the test log until further notice.
c
	      fetch = .True.

	      Do While (fetch)

	         Call LIB$Erase_Page ( 1, 1 )

c
c			Select by timerange.
c
		 Type 50
50		 Format ( ' Selecting by timerange? (Y/[N]): ', $ )
		 Accept 20, ans

		 status = STR$UpCase ( ans, ans )

		 If (ans(1:1) .Eq. 'Y') Then

	            status = FUT_Timerange ( starting_time,
     .			       start_time, ending_time, end_time )

		 Else

		    starting_time = '86001000000000'
		    ending_time = '99365235959990'

		 End If

	         cmd = 'FIND B IN TLDFI WITH (STIME GE "' //
     .			 starting_time // '" AND ETIME LE "' //
     .			 ending_time // '") SORTED BY STIME;'

		 status = DTR$Command ( DAB, cmd )
	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Select by test type.
c
		 Type 10
10		 Format ( ' Enter the test type(s): ', $ )
		 Accept 20, test_type
20		 Format (a)

		 status = STR$UpCase ( test_type, test_type )
		 status = STR$Compare ( test_type, '        ' )

		 status = STR$Translate ( test_type, test_type,
     .						' ', ',' )

		 cmd = 'FIND C IN B'

		 If (test_type .Ne. '        ') Then

		    l = STR$Find_First_Not_In_Set ( test_type, BL//TAB )
		    test_type = test_type(l:)

		    k = STR$Find_First_In_Set ( test_type, BL )

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ' WITH (TTYPE EQ "' //
     .				test_type(1:k-1) // '"'
		    test_type = test_type(k:)

		    l = STR$Find_First_Not_In_Set ( test_type, BL//TAB )

		    Do While (l .Ne. 0)
	 	      test_type = test_type(l:)
		      k = STR$Find_First_In_Set ( test_type, BL )
		      status = STR$Trim ( cmd, cmd, lcmd )
	              cmd = cmd(1:lcmd) // ' OR TTYPE EQ "' //
     .				test_type(1:k-1) // '"'
		      test_type = test_type(k:)
		      l = STR$Find_First_Not_In_Set ( test_type,
     .						      BL//TAB )
		    End Do

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ')'

		 End If

		 status = STR$Trim ( cmd, cmd, lcmd )
		 cmd = cmd(1:lcmd) // ';'

		 status = DTR$Command ( DAB, cmd )
	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Select by test number.
c
		 Type 100
100		 Format ( ' Enter the test number(s): ', $ )
		 Accept 20, test_num

		 status = STR$Compare ( test_num, '        ' )

		 status = STR$Translate ( test_num, test_num, ' ', ',' )

		 cmd = 'FIND D IN C'

		 If (test_num .Ne. '        ') Then

		    l = STR$Find_First_Not_In_Set ( test_num, BL//TAB )
		    test_num = test_num(l:)

		    k = STR$Find_First_In_Set ( test_num, BL )

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ' WITH (TNUMB EQ ' //
     .				test_num(1:k-1)
		    test_num = test_num(k:)

		    l = STR$Find_First_Not_In_Set ( test_num, BL//TAB )

		    Do While (l .Ne. 0)
	 	      test_num = test_num(l:)
		      k = STR$Find_First_In_Set ( test_num, BL )
		      status = STR$Trim ( cmd, cmd, lcmd )
	              cmd = cmd(1:lcmd) // ' OR TNUMB EQ ' //
     .				test_num(1:k-1)
		      test_num = test_num(k:)
		      l = STR$Find_First_Not_In_Set ( test_num,
     .						      BL//TAB )
		    End Do

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ')'

		 End If

		 status = STR$Trim ( cmd, cmd, lcmd )
		 cmd = cmd(1:lcmd) // ';'

		 status = DTR$Command ( DAB, cmd )
	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Select by title substring.
c
		 Type 110
110		 Format ( ' Enter the test title substring(s): ', $ )
		 Accept 20, test_title

		 status = STR$Compare ( test_title, '        ' )

		 cmd = 'FIND E IN D'

		 If (test_title .Ne. '        ') Then

		    l = STR$Find_First_Not_In_Set ( test_title,
     .							BL//TAB )
		    test_title = test_title(l:)

		    k = STR$Find_First_In_Set ( test_title, ',' ) - 1
		    If (k .Eq. -1) Then
		      status = STR$Trim ( test_title, test_title, k )
		    End If

		    Write (ttitle_len, 120) k
120		    Format (i8)

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ' WITH (TITLE CONT ' //
     .			    'FN$STR_EXTRACT("' // test_title(1:k) //
     .			    '",1,' // ttitle_len // ')'
		    test_title = test_title(k+2:)

		    l = STR$Find_First_Not_In_Set ( test_title,
     .							BL//TAB )

		    Do While (l .Ne. 0)
	 	      test_title = test_title(l:)
		      k = STR$Find_First_In_Set ( test_title, ',' ) - 1
		      If (k .Eq. -1) Then
		        status = STR$Trim ( test_title, test_title, k )
		      End If
		      Write (ttitle_len, 120) k
		      status = STR$Trim ( cmd, cmd, lcmd )
	              cmd = cmd(1:lcmd) // ' OR TITLE CONT ' //
     .			      'FN$STR_EXTRACT("' // test_title(1:k) //
     .			      '",1,' // ttitle_len // ')'
		      test_title = test_title(k+2:)
		      l = STR$Find_First_Not_In_Set ( test_title,
     .						      BL//TAB )
		    End Do

		    status = STR$Trim ( cmd, cmd, lcmd )
	            cmd = cmd(1:lcmd) // ')'

		 End If

		 status = STR$Trim ( cmd, cmd, lcmd )
		 cmd = cmd(1:lcmd) // ';'

		 status = DTR$Command ( DAB, cmd )
	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Display what has been selected.
c
	         cmd = 'PRINT TNUMB ("Number"), TTYPE ("Type"), ' //
     .			  'STIME ("Start"), ETIME ("Stop"), ' //
     .			  'TITLE ("Log Entry") USING X(25) OF E;'
	         status = DTR$Command ( DAB, cmd )
	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Record the start and stop times in an RMS file.
c
		 Type 30
30		 Format ( ' Record times in a file? (Y/[N]): ', $ )
		 Accept 20, ans

		 status = STR$UpCase ( ans, ans )

		 If (ans(1:1) .Eq. 'Y') Then

     		    cmd = 'PRINT STIME (-), ETIME (-), TNUMB (-), ' //
     .				'TTYPE (-), TITLE (-) USING X(25) ' //
     .				'OF E ON !cmd;'
	            status = DTR$Command ( DAB, cmd, qfile )
	            status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

		 End If

c
c			Release the current collection so another
c			can be formed.
c
!	         cmd = 'RELEASE B AND C AND D AND E;'
!	         status = DTR$Command ( DAB, cmd )
!	         status = DTR$DTR ( DAB, DTR$M_Opt_Cmd )

c
c			Select from the test log again.
c
		 Type 40
40		 Format ( ' Select again from the test log? ([Y]/N): ',
     .				$ )
		 Accept 20, ans

		 status = STR$UpCase ( ans, ans )

		 If (ans(1:1) .Eq. 'N') Then
		    fetch = .False.
		 End If

	      End Do

	   Else

	      Call LIB$Signal ( %Val(status) )

	   End If

	Else

	   Call LIB$Signal ( %Val(status) )

	End If

	status = DTR$Finish ( DAB )

C       Exit the program.

	if ( status .eq. ss$_normal ) then
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
	endif
	end

	End
