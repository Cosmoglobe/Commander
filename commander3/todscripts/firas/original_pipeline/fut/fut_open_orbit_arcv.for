	Integer*4 Function FUT_Open_Orbit_Arcv ( Gmt_Start, Gmt_Stop, Orb_Lun )

C-------------------------------------------------------------------------------
C
C	Purpose: To open the orbit archive for read access using
C	         a time range. The unit number is returned to the caller.
C
C	Author: Shirley M. Read
C		STX, July 1988
C
C	Invocation:
C
C	Modificaton History:
C
C	  Author		   Date		  Modification
C	------------------------------------------------------------------------
C	Fred Shuman, STX	1990 Feb 18	SPR 5278--FUT_Attitude must be
C						able to separate bad attitude
C						& bad orbit open conditions.  If
C						orbit open is bad, close it.
C
C	Input Parameters:
C	  Name		Type		Description
C	  ----------------------------------------------------------------------
C	  Gmt_Start     C*14            GMT start time for open.
C	  Gmt_Stop      C*14            GMT stop time for open.
C
C	Output Parameters:
C	  Name		Type		Description
C	  ----------------------------------------------------------------------
C	  Orb_Lun       I*4             Logical unit number.
C
C	Subroutines Called:
C	  Lib$Get_Lun
C	  Lib$Free_Lun
C         Lib%Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  $SSDef
C
C	Processing Method:
C
C	  Get the logical unit number with Lib$Get_Lun.
C	  Open the orbit data archive for read access using the dataset name
C	  and the time range.
C	  Return the unit number to the caller.
C
C------------------------------------------------------------------------------

	Implicit	None

C	Passed Parameters.

	character*14    gmt_start       ! Start time for orbit data
	character*14    gmt_stop        ! Stop time for orbit data
	integer*4       orb_lun         ! Logical unit number

C	Include files and external parameters.

	include	 '($SSDef)'

	integer*4 CT_Connect_Read
	external  CT_Connect_Read

	external  FUT_Normal

C	Functions

	integer*4  Lib$Get_Lun		! Get logical unit number
	integer*4  Lib$Free_Lun		! Free logical unit number

C	Local Declarations

	integer*4  status		! Function return status
	integer*2  ctstat(20)		! Cobetrieve return status block
	integer*4  zero / 0 /
	character*64 filename           ! Filename for open
	character*19 archive / 'CSDR$Space_Archive:' /
	character*1 slash / '/' /
	character*30 time/'00000000000000;00000000000000;'/
	character*12 dataset / 'NOR_ORB_DATA' /

C	Set return status.

	FUT_Open_Orbit_Arcv = %loc(FUT_Normal)

C	Get the unit number.

	status = lib$get_lun(orb_lun)
	if ( status .ne. SS$_Normal ) then
	  call Lib$Signal(status)
	  FUT_Open_Orbit_Arcv = status
	  return
	endif

C	Open the orbit archive for read.

	time(1:14) = gmt_start
	time(16:29) = gmt_stop
	filename = archive // dataset // slash // time

	open (  unit=orb_lun, access='sequential', status='old',
	1	file=filename, organization='sequential',
	1	iostat=status, shared, useropen=CT_Connect_Read)

	if (status .ne. zero) then
	  Call CT_Close_Arcv (, orb_lun, ctstat)
	  status = lib$free_lun(orb_lun)
	  FUT_Open_Orbit_Arcv = status
	endif

	return
	end
