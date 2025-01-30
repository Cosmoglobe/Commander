C-------------------------------------------------------------------------------ftb_print

	Integer*4 Function FTB_Print_BadSci_Info ( Channel, Sci_Rec, Rpt_Lun,
	1					Jstart, Jstop )

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS science records with
C		 badtime tags. 
C
C	Author: Shirley M. Read
C		STX, November, 1989
C
C	Invocation: Status = FTB_Print_BadSci_Info ( Channel, Sci_Rec, Rpt_Lun,
C				Jstart, Jstop )
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
C	  Jstart         C*14           GMT data start
C	  Jstop          C*14           GMT data stop
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
C	  Dump the record information.
C	  If status for write is bad, reset the function return status and
C	    signal an error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None
	Dictionary 'NFS_SDF'
        Record /NFS_SDF/ Sci_Rec

	Include 	'Ct$Library:Ctuser.Inc'
	Include		'($Ssdef)'
	Include		'(Fut_Params)'
	External 	Ct_Connect_Read	

C	Passed Parameters.

	Character*2 Channel          ! FIRAS channel
	Integer*4 Rpt_Lun            ! Report unit number
	Character*14 Jstart          ! Data start for attitude
	Character*14 Jstop	     ! Data stop for attitude

C       Watch MTM declarations

	Integer   * 2		Tlat		!  RDL latitude
	Integer   * 2		Tlon		!  RDL longitude
	Real      * 4		Coord_conv /57.3E-04/  !  Coordinate conversion
	Real      * 4		Glat		!  geocentric latitude
	Real      * 4		Glon		!  geocentric longitude
	Character * 79		Message		!  Watch_Alert message string

C  SCI Time declarations

	Character*1  Dash1(79) / 79 * '_' /
	Character*79 Dash
	Equivalence  (Dash,Dash1(1))
	Character*1  Dash2(130) / 130 * '_' /
	Character*130 Dashlong
	Equivalence  (Dashlong,Dash2(1))
	Character*14 Gmt             ! Time tag of record
	Character*14 Col_Time        ! Collect time
	Integer*4 Bin_Col(2)         ! Binary collect time
	Integer*2 Status_Bits        ! Micro header status bit word
	Integer*2 Speed              ! MTM scan speed
	Integer*2 Length             ! MTM scan length
	Integer*2 Sweeps_per_Ifg     ! MTM Sweeps per Ifg
	Integer*2 Points             ! Points process for collect cycle
        Integer*2   Iflag     	! Badtime flag
        Integer*2   TfLag	! Telemetry quality flag
	Character*14 Transmit_Gmt    ! Computed transmit time
	Integer*4  Transmit_Frame       ! Microprocessor header transmit frame
	Integer*2  Transmit_Word(2)     ! counter
	Equivalence ( Transmit_Frame, Transmit_Word(1) )
	Integer*4  Collect_Frame        ! Microprocessor header collect frame
	Integer*2  Collect_Word(2)      ! counter
	Equivalence ( Collect_Frame, Collect_Word(1) )
	Integer*4 Rstatus
	Character*19 Gmt_Char 
	Character*19 Col_Char 
	Character*1 Cbits(16)
	Integer*4 Zero / 0 /
	Integer*2 Ix

c	FUT_ATTITUDE params

	Integer*4 FUT_Attitude
	Integer*4 Att
	integer*4 Status

C	External Parameters.

	External FUT_Normal
	External FTB_Normal
	External FTB_Aberr

	FTB_Print_BadSci_Info = %loc(FTB_Normal)

        GMT = Sci_Rec.Ct_Head.Gmt
	Iflag = Sci_Rec.Collect_Time.Badtime_Flag
	Tflag = Sci_Rec.Dq_Data.Data_Quality(2)

	Do Ix = 1,2
	   Bin_Col(ix) = Sci_Rec.Collect_Time.Midpoint_Time(ix)
	Enddo 
	Call Ct_Binary_to_Gmt (Bin_Col,Col_Time)
	
	Status_Bits = Sci_Rec.Sci_Head.SC_Head3
	Speed = Sci_Rec.Sci_Head.MTM_Speed
	Length = Sci_Rec.Sci_Head.MTM_Length
	Sweeps_per_IFG = Sci_Rec.Sci_Head.SC_Head11
	Points = Sci_Rec.Sci_Head.SC_Head8

	Transmit_Word(1) = Sci_Rec.Sci_Head.SC_Head4  !LSW
	Transmit_Word(2) = Sci_Rec.Sci_Head.SC_Head5  !MSW
	Collect_Word(1) = Sci_Rec.Sci_Head.SC_Head12  !LSW
	Collect_Word(2) = Sci_Rec.Sci_Head.SC_Head13  !MSW
c
c Get the attitude.
c
	Att= Fac_Coarse
	Status = FUT_Attitude ( Sci_Rec, Att,
     .				Jstart, Jstop )

	If (Status .Ne. %Loc(FUT_Normal)) Then
	   Call LIB$Signal(%val(status))
	   Return
	End If
      
	Do Ix = 0, 15
	  Cbits(16 - Ix) = '0'
	  If ( Btest( Status_Bits, Ix )) Cbits(16 - Ix) = '1'
	Enddo
	Gmt(1:14) = Sci_Rec.Ct_Head.Gmt(1:14)
	Gmt_Char= Gmt(1:2)//':'//Gmt(3:5)//':'//Gmt(6:7)//':'//
	1	  Gmt(8:9)//':'//Gmt(10:11)//'.'//Gmt(12:14)
	Col_Char= Col_Time(1:2)//':'//Col_Time(3:5)//':'//Col_Time(6:7)//':'//
	1	  Col_Time(8:9)//':'//Col_Time(10:11)//'.'//Col_Time(12:14)

C Write science information in report file.

        Write( Unit=Rpt_lun, Fmt=150, Iostat=Rstatus) Dashlong
 150    Format(1x,130a)
	Write (Unit=Rpt_Lun, FMT=175, Iostat=Rstatus ) Gmt, Col_Time
  175   Format(1x,'Transmit_Gmt: ',a,' Collect_Gmt: ',a )

********************************************************************************
*
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
	   Tlon = Sci_rec.Attitude.Terr_Longitude 
	   Glat = Float(Tlat)*Coord_Conv
	   Glon = Float(Tlon)*Coord_Conv
	   Write (Message,20) Glat, Glon
   20	   Format (' **** GEOCENTRIC COORDINATES ****  latitude = ',
     &			f6.1,' longitude = ',f6.1)
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Message
 200      Format(1x,a)

	If ( FTB_Print_BadSci_Info .EQ. %loc(FTB_Normal)) Then
	  Write ( Unit=Rpt_lun, FMT=100, Iostat=Rstatus )
 100      Format (1x,'Mid_Point_of_Collect',4x,
	1	'Status_Bits VAX Ms->Ls',1x,'Speed',1x,
	2	'Len',1x,'Transmit_MNF',1x,
	3	'Start_Col_MNF',1x,'Sw',2x,'#Pts.',1x,
	4	'Frame Transmit_Gmt',1x,'Tlmq', ' TimFlg')
          If (Rstatus .NE. Zero ) Then
	    FTB_Print_BadSci_Info = %loc(FTB_Aberr)
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	If ( FTB_Print_BadSci_Info .EQ. %loc(FTB_Normal)) Then
	  Write (Unit=Rpt_Lun, FMT=300, Iostat=Rstatus ) Col_Char, Cbits,
	1     Speed, Length, Transmit_Frame, Collect_Frame,
	2     Sweeps_per_Ifg, Points, Gmt_Char, TFlag, Iflag
 300      Format(1x,a,5x,4(4a,1x),2x,i4,1x,i3,1x,i12,1x,i12,
	1	 3x,i3,1x,i6,1x,a,1x,i3,1x,i4)
	    If ( Rstatus .NE. Zero ) Then
!	    Allow print even if overflow.
	      Call Lib$Signal(%val(Rstatus))
            Endif
	Endif

	Return
	End
