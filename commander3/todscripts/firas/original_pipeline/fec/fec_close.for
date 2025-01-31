	integer*4 function FEC_CLOSE ( Report_flag,
     .                                 CT_Unit,
     .                                 LU_Report,
     .                                 LU_Short_Sci,
     .                                 Start_Chan,
     .                                 Stop_Chan,
     .                                 Skip_Chan,
     .                                 NoData_Flag,
     .                                 Write_Flag)
c
c      By K. Jensen, STX, 25-FEB-1991
c
c      Revision of FEC_CLOSE_FILES_A.FOR
c      Written by Connie Lau, SASC Technologies Inc., May 1986
c      Commented by Reid Wilson, SASC Technologies Inc.,  22-AUG-1986 
c
c      Revision made to distinguish final requirements FEC from older
c      versions of FEC.
c
c      Modification History:
c
c      Jan 23, 1989 D. Bouler, STX
c          Added report_flag to parameter list in order to signal
c          presence/absence of REPORT qualifier. If no REPORT 
c          qualifier (report_flag = 0) do not close lu_report.
c
c      June 6, 1989 R. Kummerer, STX
c          SPR 3972, Replace LIB$GET_LUN / LIB$FREE_LUN with analogous FUT
c	   routines.
c
c      01/17/92  K. Jensen, STX
c          SPR 9374, Add a flag to enable a graceful abort if there is
c                    no data.
c
c      Description:
c            This function closes the files opened in FEC_OPEN.
c            The input file closed here should be the archive file.
c
c      Format:
c            Tstatus =  FEC_Close ( Report_flag,
c                                   CT_Unit,
c                                   LU_Report,
c                                   LU_Short_Sci,
c                                   Start_Chan,
c                                   Stop_Chan,
c                                   Skip_Chan,
c                                   NoData_Flag,
c                                   Write_Flag)
c
c      Input Parameters:
c
c            Report_flag         = signals presence/absence of REPORT
c                                  integer*4
c            CT_Unit             = logical unit connected to archive file
c                                  integer*2
c            LU_Report           = logical unit connected to summary file
c                                  integer*2
c            LU_Short_Sci(4)     = logical unit connected to output records file
c                                  integer*2
c            LU_Index            = logical unit connected to index file
c                                  integer*2
c	     Start_Chan		 = start channel; integer*4
c	     Stop_Chan		 = stop channel; integer*4
c	     Skip_Chan		 = increment channel; integer*4
c	     NoData_Flag	 = If no data found, bypass close_archive
c	     Write_Flag		 = Switch controlling writes to short science
c				   archive;  integer*4
c      
c      Output Parameters:
c            Tstatus            = status of operation
c                                  integer*4
c
c      Modules Called:
c            ct_close_arcv               = close archive file
c            errsns                      = translate system error to system 
c                                          message code  
c            FUT_FREE_LUN                = FIRAS system routine 
c
c      Include Files:
c            CT$:LIBRARY:CTUSER.INC      = COBETRIEVE routines
c            FEC_INC                     = parameters used by FEC
c            FEC_MSG                     = message definition for FEC
c
c
c      *** BEGIN ***
c
        implicit none
c
	include  'CT$LIBRARY:ctuser.inc'      ! Access COBETRIEVE routines
        include  '(FEC_INC)'                  ! parameters used by FEC
        include  '(FEC_MSG)'                  ! message definition
        include  '(FUT_PARAMS)'               ! Global FIRAS parameters
        include  '(FUT_ERROR)'
        include  '($SSDEF)'
c
        integer*4       report_flag           ! presence/absence of REPORT
	integer*4	ct_unit               ! logical unit for archive file
	integer*2	ct_status(20)         ! COBETRIEVE return status holder
	integer*4	LU_Report             ! logical unit for summary file
	integer*4	LU_Short_Sci(4)       ! logical unit for records file
c
	integer*4	Start_Chan
	integer*4	Stop_Chan
	integer*4	Skip_Chan
        integer*4       NoData_Flag
        integer*4       Write_Flag 	      ! Were the short science archives
					      !   open?
        integer*4       ios                   ! return status of operations
        integer*4       counter               ! loop control variable

	integer*4	status
        integer*4       tstatus
	integer*4	FUT_Free_LUN

        tstatus = %loc(FEC_NORMAL)

 	Call CT_Close_Arcv (, ct_unit, ct_status)
	If (ct_status(1) .Ne. ctp_normal) Then
	  Call LIB$Signal(FEC_CTCLOSEENG)
          tstatus=%loc(FEC_CTCLOSEENG)
        Endif
	status = FUT_Free_LUN(ct_unit)
        If (.Not. status) Then
          Call LIB$Signal(%Val(status))
          Tstatus = status
        End If

	If ( Write_flag .Eq. 1 ) Then

           Do Counter = Start_Chan, Stop_Chan, Skip_Chan
              Call CT_Close_Arcv (, LU_Short_Sci(counter), ct_status)
              If ( NoData_Flag .eq. 0 ) Then
	        If (ct_status(1) .Ne. ctp_normal) Then
                   Call LIB$Signal(FEC_CTCLOSESS,%Val(1),%Val(Counter))
                   tstatus=%loc(FEC_CTCLOSESS)
                Endif
 	        status = FUT_Free_LUN(LU_Short_Sci(Counter))
                If (.Not. status) Then
                   Call LIB$Signal(%Val(status))
                   Tstatus = status
	        End If
              Endif
           End Do

	End If

c
c       Check for presence of REPORT qualifier; if absent, skip close
c       of lu_report. D. Bouler Jan 20,1989
c

        if (report_flag .eq. 1) then
           If ( Tstatus ) Call LIB$Signal(FEC_Normal)
           If ( .not. Tstatus ) Call LIB$Signal(%Val(SS$_Abort))
           Close (unit=LU_Report, iostat=ios)
           If (ios .Ne. 0) then
	      Call LIB$Signal(FEC_CLOSEREP)
              Tstatus = %loc(FEC_CLOSEREP)
           Endif
           FUT_Report_Lun = 0

           status = FUT_Free_LUN (LU_Report)
           If (.Not. status) Then
              Call LIB$Signal(%Val(status))
              Tstatus = status
	   End If
        end if

        If ( .not. Tstatus ) then
           FEC_Close = Tstatus
        Else
           FEC_Close = %loc(FEC_NORMAL)
        Endif

	RETURN
	END
