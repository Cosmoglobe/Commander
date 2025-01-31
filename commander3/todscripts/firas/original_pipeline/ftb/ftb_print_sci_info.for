C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Print_Sci_Info ( Channel, Sci_Rec, Rpt_Lun,
	1					jstart, jstop )

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS science records with
C		 valid time tags. 
C
C	Author: Shirley M. Read
C		STX, November, 1989
C
C	Invocation: Status = FTB_Print_Sci_Info ( Channel, Sci_Rec, Rpt_Lun,
C			     Jstart, Jstop)
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Channel        C*2            FIRAS channel
C	  Sci_Rec        NFS_SDF        RDL
C	  Rpt_Lun        I*4            Report unit number
C	  Jstart         C*14           GMT start
C	  Jstop          C*14           GMT stop
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C	  Set function return status.
C	  Dump the record information: transmit time, collect time, badtime
C	     flag, telemetry quality flag and terrestrial latitude and
C	     longitude.
C	  If status for write is bad, reset the function return status and
C	    signal an error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None
	dictionary 'NFS_SDF'
        record /NFS_SDF/ sci_rec

	include 	'ct$library:ctuser.inc'
	include		'($ssdef)'
	include		'(fut_params)'
	external 	ct_connect_read	

C	Passed Parameters.

	Character*2 Channel          ! FIRAS channel
	Integer*4 Rpt_Lun            ! Report unit number
	Character*14 Jstart	     ! GMT start
	Character*14 Jstop	     ! GMT stop

C       Watch MTM declarations

	integer   * 2		tlat		!  RDL latitude
	integer   * 2		tlon		!  RDL longitude

	real      * 4		coord_conv /57.3E-04/  !  Coordinate conversion
	real      * 4		glat		!  geocentric latitude
	real      * 4           glon            !  geocentric longitude
	character * 79		message		!  Watch_Alert message string

C  SCI Time declarations

	Character*1  Pound1(79) / 79 * '#' /
	Character*79 Pound
	Equivalence  (Pound,Pound1(1))
	Character*14 Gmt             ! Time tag of record
	Character*14 Col_Time        ! Collect time
	Integer*4 Bin_Col(2)         ! Binary collect time

C	Local Declarations.

	Integer*4 Rstatus
	Integer*4 Zero / 0 /
	Integer*2 Ix
        Integer*2   Iflag	     ! Badtime flag
        Integer*2   Tflag	     ! Telemetry quality flag
        Character*1 Dash(79)/ 79 * '_'/
        Character*79 Dashes                                     
	Equivalence  (Dashes, Dash(1))
              

C	External Parameters.

	External FUT_Normal
	External FTB_Normal
	External FTB_Aberr

	FTB_Print_Sci_Info = %loc(FTB_Normal)
        GMT = Sci_Rec.Ct_Head.Gmt
	Iflag = Sci_Rec.Collect_Time.Badtime_Flag
	Tflag = Sci_Rec.DQ_Data.Data_Quality(2)
	Do Ix = 1,2
	   Bin_Col(ix) = Sci_Rec.Collect_Time.Midpoint_Time(ix)
	Enddo 
	Call Ct_Binary_to_Gmt (Bin_Col,Col_Time)

	Write (Unit=Rpt_Lun, FMT=160, Iostat=Rstatus) Dashes
 160    Format(1x,a)	
	Write (Unit=Rpt_Lun, FMT=175, Iostat=Rstatus ) Gmt, 
	1	Col_Time, Tflag, Iflag
 175    Format(1x,'Transmit GMT: ',a,' Collect_GMT: ',a,
	1	' Tlmq= ', I5, ' Timflg= ',I3) 

********************************************************************************
*       The output format of this science record dump is the same as the
*  Geocentric latitude and longitude for the WATCHMTM Program.
*	The subroutine WATCH_GEOCENTRIC takes a time for an MTM event, as 
*  determined by WATCHMTM, and finds the geocentric coordinates of the COBE
*  spacecraft at the time of the event.  The routine works by assuming a CT
*  start and stop time of the event time +/- 40 seconds.  It then opens the
*  FDQ_SDF_RH record within that time and reads the geocentric coordinates of 
*  the spacecraft.  The coordinates are converted to degrees and written to the
*  WATCHMTM output file by the subroutine WATCH_ALERT.
*
********************************************************************************
c
c  Decode the science records and extract the geocentric coordinates.
c
	  Tlat = Sci_Rec.Attitude.Terr_Latitude
	  Tlon = Sci_Rec.Attitude.Terr_Longitude
	  Glat = Float(Tlat)*Coord_Conv
	  Glon = Float(Tlon)*Coord_Conv
	  Write (Message,20) Glat, Glon
   20	  Format (' **** GEOCENTRIC COORDINATES ****  latitude = ',
     &			f6.1,' longitude = ',f6.1)
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Message
 200      Format(1x,a)

	Return
	End
