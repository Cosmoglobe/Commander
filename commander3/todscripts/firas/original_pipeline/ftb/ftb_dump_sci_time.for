C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Dump_Sci_Time ( Channel, Gmt, Status_Bits,
	1	  Absum, Speed, Length, Transmit_Frame, Collect_Frame,
	2	  Sweeps_per_Ifg, Points, Transmit_Gmt, Numrec, 
	3	  Time_Flag, Rpt_Lun)

C-------------------------------------------------------------------------------
C
C	Purpose: To dump the information from the FIRAS science records with
C		 invalid time tags. The records bracketing the bad time periods
C		 are also included. 
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Dump_Sci_Time ( Channel, Gmt, Status_Bits,
C			     Absum, Speed, Length, Transmit_Frame, 
C			     Collect_Frame, Sweeps_per_Ifg, Points,
C	                     Transmit_Gmt, Numrec, Time_Flag, Rpt_Lun )
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
C	  Speed         I*2             MTM scan speed
C	  Length        I*2             MTM scan length
C	  Transmit_Frame I*4            Micro counter for transmit frame
C	  Collect_Frame  I*4            Micro counter for start collect frame
C	  Sweeps_per_Ifg I*2            MTM Sweeps per Ifg
C	  Points         I*2            Points process for collect cycle
C	  Transmit_Gmt   C*14           Computed transmit time
C	  Numrec         I*4            Number of records with repeated time tag
C	  Time_Flag      C*2            Flag indicating reason for dump
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
	dictionary 'NFS_SDF'
        record /NFS_SDF/ sci_rec

C	Passed Parameters.

	Character*2 Channel          ! FIRAS channel
	Character*14 Gmt             ! Time tag of record
	Integer*2 Status_Bits        ! Micro header status bit word
	Integer*4 Absum              ! Sum of absolute values of Ifg points
	Integer*2 Speed              ! MTM scan speed
	Integer*2 Length             ! MTM scan length
	Integer*4 Transmit_Frame     ! Micro counter for transmit frame
	Integer*4 Collect_Frame      ! Micro counter for start collect frame
	Integer*2 Sweeps_per_Ifg     ! MTM Sweeps per Ifg
	Integer*2 Points             ! Points process for collect cycle
	Character*14 Transmit_Gmt    ! Computed transmit time
	Integer*4 Numrec             ! Number of records with repeated time tag
	Character*2 Time_Flag        ! Flag indicating reason for dump
	Integer*4 Rpt_Lun            ! Report unit number

C	Local Declarations.

	Logical*1 First / .True. /
	Integer*4 Rstatus
	Character*19 Gmt_Char 
	Character*14 Tran
	Character*19 Tran_Char 
	Character*1 Cbits(16)
	Integer*4 Zero / 0 /
	Integer*2 Ix
        Character*1 Dash(130)/ 130 * '_'/
        Integer*4   iflag


C	External Parameters.

	External FTB_Normal
	External FTB_Aberr

	FTB_Dump_Sci_Time = %loc(FTB_Normal)
        iflag = sci_rec.collect_time.badtime_flag
       
	If ( First ) Then
	  Write ( Unit=Rpt_lun, FMT=100, Iostat=Rstatus )
 100      Format (1x,'Mid_Point_of_Collect',4x,
	1	  'Status_Bits VAX Ms->Ls',2x,'IFG_Sum',2x,
	2  	  'S',1x,'L',1x,'Transmit_MNF',1x,
	3	  'Start_Col_MNF',1x,'Sw',2x,'#Pts.',1x,
	4	'Comp. Transmit_Gmt',1x,'#Rec.',1x,'TF', ' Flag')
          Write( Unit=Rpt_lun, Fmt=150, Iostat=Rstatus) Dash
 150      Format(1x,130a)
          First = .False.
           If (Rstatus .NE. Zero ) Then
	    FTB_Dump_Sci_Time = %loc(FTB_Aberr)
	    Call Lib$Signal(%val(Rstatus))
	  Endif
	Endif
	Do Ix = 0, 15
	  Cbits(16 - Ix) = '0'
	  If ( Btest( Status_Bits, Ix )) Cbits(16 - Ix) = '1'
	Enddo
	Tran(1:14) = Transmit_Gmt(1:14)
	Gmt_Char= Gmt(1:2)//':'//Gmt(3:5)//':'//Gmt(6:7)//':'//
	1	  Gmt(8:9)//':'//Gmt(10:11)//'.'//Gmt(12:14)
	Tran_Char= Tran(1:2)//':'//Tran(3:5)//':'//Tran(6:7)//':'//
	1	  Tran(8:9)//':'//Tran(10:11)//'.'//Tran(12:14)

	If ( FTB_Dump_Sci_Time .EQ. %loc(FTB_Normal)) Then
          IF (TIME_FLAG .NE. 'LT' .OR. TIME_FLAG .NE. 'EQ') THEN
	  Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus ) Gmt_Char, Cbits,
	1     Absum, Speed, Length, Transmit_Frame, Collect_Frame,
	2     Sweeps_per_Ifg, Points, Tran_Char, Numrec, Time_Flag,iflag
 200      Format(1x,a,5x,4(4a,1x),1x,i10,2x,i1,1x,i1,1x,i12,1x,i12,
	1	 1x,i3,1x,i6,1x,a,1x,i4,1x,a,i3)
	    If ( Rstatus .NE. Zero ) Then
!	    Allow print even if overflow.
	      Call Lib$Signal(%val(Rstatus))
            Endif
	  ENDIF     ! TIME_FLAG
	Endif
	Return
	End
