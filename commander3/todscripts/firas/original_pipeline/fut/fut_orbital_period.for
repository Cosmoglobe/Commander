C-------------------------------------------------------------------------------

	Integer*4 Function FUT_Orbital_Period ( Gmt_Stop, Orbital_Period ) 

C-------------------------------------------------------------------------------
C
C	Purpose: To open the orbit header archive for read access using 
C	         a time range, read the appropriate header record and then
C		 determine the COBE orbital period.
C
C	Author: R. Kummerer
C		STX, July, 1989
C
C	Invocation:
C
C	Modificaton History:
C
C	  Author	    Date	  Modification
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Gmt_Stop      C*14            GMT stop time for open.
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Period        R*4             COBE orbital period in minutes.
C	
C	Subroutines Called:
C	  Lib$Get_Lun
C         Lib$Signal
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
C	  Open the orbit header archive for read access using the dataset name
C	  	and the time range.
C	  Read the latest orbit header record with CT_READ_ARCV_REV.
C	  Calculate the COBE orbital period from the Mean Motion.
C	  Close the orbital header archive.
C
C------------------------------------------------------------------------------
C	

	implicit	none

C	Passed Parameters.

	character*14    gmt_stop        ! Stop time for orbit data
	integer*4       orb_lun         ! Logical unit number

C	Include files and external parameters.

	include	 '($SSDef)'
	include  'CT$Library:CTUser.Inc'
	include  '(FUT_Params)'

	integer*4 CT_Connect_Read
	external  CT_Connect_Read

	external  FUT_Normal
	external  FUT_OpenOrbHdr
	external  FUT_ReadOrbHdr
	external  FUT_CloseOrbHdr
	
C 	Functions

	integer*4  Lib$Get_Lun  	! Get logical unit number
	integer*4  Lib$Free_Lun  	! Free logical unit number

C	Local Declarations
	
	dictionary 'NOR_ORB_HDR'
	record /NOR_ORB_HDR/ ORB_HDR

	integer*4  status    		! Function return status
	integer*2  ct_stat(20)
	integer*4  zero / 0 /
	character*64 filename           ! Filename for open
	character*19 archive / 'CSDR$Space_Archive:' /
	character*1 slash / '/' /
        character*30 time/'86001000000000;00000000000000;'/
	character*11 dataset / 'NOR_ORB_HDR' /
	real*4  Orbital_Period

C     	Set return status.

	FUT_Orbital_Period = %loc(FUT_Normal)
	Orbital_Period = 0.
						 
C 	Get the unit number.

	status = lib$get_lun(orb_lun)
	if ( status .ne. SS$_Normal ) then
	  call Lib$Signal(status)
          FUT_Orbital_Period = status
	  return
	endif

C     	Open the attitude archive for read.      

	time(16:29) = gmt_stop
	filename = archive // dataset // slash // time

 	open (  unit=orb_lun, access='sequential', status='old',
	1	file=filename, organization='sequential',
	1	iostat=status, shared, useropen=CT_Connect_Read)

	if (status .ne. zero) then
	  call Lib$Signal ( FUT_OpenOrbHdr, %val(1), %val(status) )
	  FUT_Orbital_Period = %Loc(FUT_OpenOrbHdr)
	  return
	endif

C	Read the latest orbit header record.

	call ct_read_arcv_rev ( , orb_lun, ORB_HDR, ct_stat )

	if (ct_stat(1) .ne. ctp_normal .and.
     .	    ct_stat(1) .ne. ctp_beginoffile) then
	  call Lib$Signal ( FUT_ReadOrbHdr, %val(1), %val(ct_stat(1)) )
	  FUT_Orbital_Period = %Loc(FUT_ReadOrbHdr)
	endif

C	Close the orbit header archive.

	call ct_close_arcv ( , orb_lun, ct_stat )

	if (ct_stat(1) .ne. ctp_normal) then
	  call Lib$Signal ( FUT_CloseOrbHdr, %val(1), %val(ct_stat(1)) )
	endif

C	Calculate the orbital period in minutes:
C
C		2PI   864 sec    1 min
C	    P = --- * ------- * ------   where W is COBE mean angular velocity
C		 W     1 DUT    60 sec

	if (FUT_Orbital_Period) then
	  Orbital_Period = ( 2. * fac_pi / ORB_HDR.Mean_Motion ) * 864. / 60.
	else
	  Orbital_Period = 0.
	endif

	return
	end
