C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Print_Eng_Info ( Field, Eng_Rec, Fieldval,
	1	           Scan_Mode, Dom_Mode, ERpt_Lun, Start, Stop )

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS engineering records for
C		 the selected engineering field.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Print_Eng_Info ( Field, Eng_Rec, Fieldval,
C		             Scan_Mode, Dom_Mode, ERpt_Lun, Start, Stop )
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
C	  Field          C*30           FIRAS engineering field
C	  Eng_Rec        FDQ_ENG        RDL
C	  Fieldval       R*4            Engineering value
C	  Scan_Mode(4)   I*2            Four channel scan modes
C	  Dom_Mode       I*2            Predominant scan mode
C	  ERpt_Lun       I*4            Report unit number
C	  Start          C*14           Start of data for latitude/longitude
C	  Stop           C*14           Stop of data for latitude/longitude
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
C	  Set function return             
C	  Call FUT routine to get the geocentric latitude and longitude.
C	  Dump the record information: Lat/Long, Fieldname, Engineering Value, 
C	     MTM scan mode.
C	  If status for write is bad, reset the function return status and
C	    signal an error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

	Include 	'Ct$Library:Ctuser.Inc'
	Include		'($Ssdef)'
	Include		'(Fut_Params)'
	External 	Ct_Connect_Read	

C	Passed Parameters.

	Character*30 Field            ! Engineering field name
	Dictionary   'FDQ_ENG'
	Record       /FDQ_ENG/ Eng_Rec
	Real*4       Fieldval         ! Engineering field value
	Integer*2    Scan_Mode(4)     ! Channel scan modes
	Integer*2    Dom_Mode         ! Predominanrt scan mode
	Integer*4    ERpt_Lun         ! Report unit number
	Character*14 Start            ! Data start for attitude
	Character*14 Stop	      ! Data stop for attitude



C  Local declarations

	Character*1  Dash1(79) / 79 * '_' /
	Character*79 Dash
        Character*26 Dash26 /'--------------------------'/
        Character*9  dash9/'---------'/
	Equivalence  (Dash,Dash1(1))
	Character*1  Dash2(130) / 130 * '_' /
	Character*130 Dashlong
	Equivalence  (Dashlong,Dash2(1))
	Integer*4    Rstatus          		! Function return status
	Integer*4 Zero / 0 /
	Integer*2 Ix
        Integer*2 Count                         ! Counter
        Logical*1 First_time /.True./

c	FUT_LonLat params

	Integer*4 FUT_LonLat
	Integer*4 Status
	Integer   * 2		Tlat		!  RDL latitude
	Integer   * 2		Tlon		!  RDL longitude
	Real      * 4		Coord_conv /57.3E-04/  !  Coordinate conversion
	Real      * 4		Glat		!  geocentric latitude
	Real      * 4		Glon		!  geocentric longitude
	Character * 100		Message		!  Watch_Alert message string

C	External Parameters.

	External FUT_Normal
	External FTB_Normal
	External FTB_Aberr

	FTB_Print_Eng_Info = %loc(FTB_Normal)

C 	Print header page title and subtitle.

        If (First_time) then
	Write (Unit=ERpt_Lun, FMT=10, Iostat=Rstatus ) start, stop
 10         Format (1x//5x,'Record Dump from FIRAS Engineering Records for ',
	1	'Geocentric Latitude and Longitude'//5x,
	1	'Fieldname, Engineering values, MTM Scan mode, and ',
	1       'Predominant Scan Mode.'//,5x
   	1       'COBETRIEVE Time Selected is ',a,' to ',a//)

C	

	Write (Unit=ERpt_Lun, FMT=150, Iostat=Rstatus ) 
 150        Format (1x,'Geocentric Lat & Longitude      Fieldname',
	1              '          Eng_Value  Scan Mode   Pred_mode')

	Write (Unit=ERpt_Lun, FMT=175, Iostat=Rstatus )
	1     Dash26, Dash9, Dash9, Dash9, Dash9
 175	    Format (1x,a,6x,a,10x,a,2x,a,3x,a/)
        First_time = .False.
	Count = 0
        endif

c Get the latitude and longitude.
c
	Status = FUT_LonLat (Eng_Rec.Ct_Head.Gmt, Start, Stop, Glon, GLat)

	If (Status .Ne. %Loc(FUT_Normal)) Then
	   Call LIB$Signal(%val(status))
	   Return
	End If

	If (.Not.First_time) Then
           If ( Count .LE. 46) Then   
	      Write (Message,20) Glat,Glon,Field,Fieldval,Scan_Mode,Dom_Mode
  20	      Format ('* GEO LAT ',F6.1,' LON ',F6.1,1x,A21,2x,F7.3,
	1     3x,' S',4I2,2x,' P',I2)
	      Write (Unit=ERpt_Lun, FMT=200, Iostat=Rstatus ) Message
 200          Format(1x,a)
              Count = Count + 1
	   Else
               Count = 0
	       Write (Unit=ERpt_Lun, FMT=240, Iostat=Rstatus ) start, stop
 240           Format (1H1//15x,'COBETRIEVE Time Selected is ',a,' to ',a//)
	       Write (Unit=ERpt_Lun, FMT=250, Iostat=Rstatus ) 
 250           Format (1x,'Geocentric Lat & Longitude      Fieldname',
	1              '          Eng_Value  Scan Mode   Pred_mode')

	      Write (Unit=ERpt_Lun, FMT=275, Iostat=Rstatus )
	1     Dash26, Dash9, Dash9, Dash9, Dash9
 275	      Format (1x,a,6x,a,10x,a,2x,a,3x,a/)
	      Write (Message,300) Glat,Glon,Field,Fieldval,Scan_Mode,Dom_Mode
 300	      Format ('* GEO LAT ',F6.1,' LON ',F6.1,1x,A21,2x,F7.3,
	1     3x,' S',4I2,2x,' P',I2)
	      Write (Unit=ERpt_Lun, FMT=400, Iostat=Rstatus ) Message
 400          Format(1x,a)
              Count = Count + 1
	   Endif
  	Endif
	Return
	End
