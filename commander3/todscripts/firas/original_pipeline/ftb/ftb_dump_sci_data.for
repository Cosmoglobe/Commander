C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Dump_Sci_Data ( Channel, Gmt, Status_Bits,
	1	  Absum, Gain, Speed, Length, Transmit_Frame, Collect_Frame,
	2	  Sweeps_per_Ifg, Sampcol, Points, Sampergr, Collect_Gmt, 
	3	  Numrec, Tlm_Flag, Badtime, Rpt_Lun)

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS science records with
C		 invalid time tags. The records bracketing the bad time periods
C		 are also included. 
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Dump_Sci_Data ( Channel, Gmt, Status_Bits,
C			     Absum, Speed, Length, Transmit_Frame, 
C			     Collect_Frame, Sweeps_per_Ifg, Sampcol, Points,
C	                     Sampergr, Collect_Gmt, Numrec, Tlm_Flag, 
C			     Badtime, Rpt_Lun )
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
C	  Channel       C*2             FIRAS channel
C	  Gmt           C*14            Time tag of record
C	  Status_Bits   I*2             Micro header status bit word
C	  Absum         I*4             Sum of absolute values of Ifg points
C	  Gain          I*2             Science gain
C	  Speed         I*2             MTM scan speed
C	  Length        I*2             MTM scan length
C	  Transmit_Frame I*4            Micro counter for transmit frame
C	  Collect_Frame  I*4            Micro counter for start collect frame
C	  Sweeps_per_Ifg I*2            MTM Sweeps per Ifg
C	  Sampcol        I*2            ADC samples collected for cycle
C	  Points         I*2            ADC samples processed for collect cycle
C	  Sampergr       I*2            Samples averaged per group
C	  Collect_Gmt    C*14           Computed collect time
C	  Numrec         I*4            Number of records with repeated time tag
C	  Tlm_Flag       C*2            Flag indicating telemetry quality
C	  Badtime        C*1            Bad collect time flag
C	  Rpt_Lun        I*4            Report unit number
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
C	  Load the Gmt values into a string with the print format.
C	  Determine which bits are set in the status word and load the 
C	    bit values into an array for printing.
C	  If first time, print the table header.
C	  Dump the record information.
C	  If status for write is bad, reset the function return status and
C	    signal an error.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

C	Passed Parameters.

	Character*2 Channel          ! FIRAS channel
	Character*14 Gmt             ! Time tag of record
	Integer*2 Status_Bits        ! Micro header status bit word
	Integer*4 Absum              ! Sum of absolute values of Ifg points
	Integer*2 Gain		     ! Science gain
	Integer*2 Speed              ! MTM scan speed
	Integer*2 Length             ! MTM scan length
	Integer*4 Transmit_Frame     ! Micro counter for transmit frame
	Integer*4 Collect_Frame      ! Micro counter for start collect frame
	Integer*2 Sweeps_per_Ifg     ! MTM Sweeps per Ifg
	Integer*2 Sampcol            ! ADC samples collected for collect cycle
	Integer*2 Points             ! ADC samples processed for collect cycle
	Integer*2 Sampergr           ! Samples averaged for group
	Character*14 Collect_Gmt     ! Computed collect time
	Integer*4 Numrec             ! Number of records with repeated time tag
	Character*2 Tlm_Flag         ! Flag indicating telemetry quality
	Character*1 Badtime          ! Bad collect time flag
	Integer*4 Rpt_Lun            ! Report unit number

C	Local Declarations.

	Logical*1 First / .True. /
	Integer*4 Rstatus
	Character*19 Gmt_Char 
	Character*14 Coll
	Character*19 Coll_Char 
	Character*1 Cbits(16)
	Integer*4 Zero / 0 /
	Integer*2 Ix
        Character*1 Dash(131)/ 131 * '_'/


C	External Parameters.

	External FTB_Normal
	External FTB_Aberr

	FTB_Dump_Sci_Data = %loc(FTB_Normal)

	If ( First ) Then
	  Write ( Unit=Rpt_lun, FMT=100, Iostat=Rstatus )
 100      Format (1x,'Frame_Transmit_Time',1x,
	1	  'Stat_Bits VAX Ms->Ls',2x,'IFG_Sum',1x,'Gain',
	2  	  1x,'S',1x,'L',1x,'TransmitMF',1x,
	3	  'StartCol_MF',1x,'Sw',1x,'#SpCol',1x,'#SpPro',1x,'AD',
	4	  1x,'Midpoint_of_Collect',1x,'#Rec',1x,'B',1x,'TQ')
          Write( Unit=Rpt_lun, Fmt=150, Iostat=Rstatus) Dash
 150      Format(1x,131a)
          First = .False.
           If (Rstatus .NE. Zero ) Then
	    FTB_Dump_Sci_Data = %loc(FTB_Aberr)
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Do Ix = 0, 15
	  Cbits(16 - Ix) = '0'
	  If ( Btest( Status_Bits, Ix )) Cbits(16 - Ix) = '1'
	Enddo
	Coll(1:14) = Collect_Gmt(1:14)
	Gmt_Char= Gmt(1:2)//':'//Gmt(3:5)//':'//Gmt(6:7)//':'//
	1	  Gmt(8:9)//':'//Gmt(10:11)//'.'//Gmt(12:14)
	Coll_Char= Coll(1:2)//':'//Coll(3:5)//':'//Coll(6:7)//':'//
	1	  Coll(8:9)//':'//Coll(10:11)//'.'//Coll(12:14)

	If ( FTB_Dump_Sci_Data .EQ. %loc(FTB_Normal)) Then
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Gmt_Char, Cbits,
	1     Absum, Gain, Speed, Length, Transmit_Frame, Collect_Frame,
	2     Sweeps_per_Ifg, Sampcol, Points, Sampergr, Coll_Char, 
	3     Numrec, Badtime, Tlm_FLag
 200      Format(1x,a,1x,4(4a,1x),i9,I5,1x,i1,1x,i1,i11,i11,
	1	 1x,i3,1x,i6,1x,i6,i3,1x,a,i5,1x,a,1x,a)
	  If ( Rstatus .NE. Zero ) Then
!	    Allow print even if overflow.
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Return
	End
