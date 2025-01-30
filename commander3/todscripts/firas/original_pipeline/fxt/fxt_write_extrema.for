C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Write_Extrema ( Xtrm_Rec, Lun_Xtrm,
	1	  Report, Lun_Rpt, Xtrm_Table ) 

C-------------------------------------------------------------------------------
C
C	Purpose: To write the FXT_Eng_Xtrm record to the archive.
C
C	Author: Shirley M. Read
C		STX, December, 1988
C
C	Invocation: Status = FXT_Write_Extrema ( Xtrm_Rec, Lun_Xtrm,
C		    Report, Lun_Rpt, Xtrm_Table )
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
C	  Name		Type		   Description
C 	  ----------------------------------------------------------------------
C	  Xtrm_Rec      Structured Record  FXT_Eng_Xtrm Ct record
C         Lun_Xtrm      I*4                Logical unit to write FXT_Eng_Xtrm
C         Report        L*1                Flag to write report
C 	  Lun_Rpt       I*4                Logical unit to write report
C         Xtrm_Table(118,2) L*1            Table of flags indicating extrema
C				           are exceeded for the current run
C	Output Parameters:
C	  Name		Type		   Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	  CT_Write_Arcv
C	  Lib$Signal 
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	  FXT_Msg.Txt
C	  CT$Library:CTUser.Inc
C	  FUT_Plot_Names.Txt
C
C	Processing Method:
C	  
C	  Set the function status to FXT_Normal.
C	  Write the FXT_Eng_Xtrm record to the archive.
C	  If the return status is not success,
C	    Signal an error.
C     	    Set the function status to FXT_Aberr.
C	  Endif
C	  If the function status is normal and a report is requested,
C	    Write the field names of the engineering analogs for which
C           the minima or the maxima were exceeded in the report file.
C           If the write status is not success,
C	      Log the error and set the function status to FXT_Aberr.
C           Endif
C         Endif
C	  Return
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Dictionary 'FXT_Eng_Xtrm'
	Record /FXT_Eng_Xtrm/ Xtrm_Rec 	!FXT_Eng_Xtrm Ct record

	Integer*4   Lun_Xtrm            ! Logical unit to write FXT_Eng_Xtrm
	Logical*1   Report              ! Flag to write report
	Integer*4   Lun_Rpt             ! Logical unit to write report
	Logical*1   Xtrm_Table(118,2)   ! Table of flags indicating extrema

!	Include Files

	Include '(FXT_Msg)'
	Include 'CT$library:CTUser.Inc'
	Include '(FUT_Plot_Names)'

!	Local Declarations.

	Integer*4   Status		! Return status
	Integer*4   Zero / 0 /          ! Zero value
	Integer*2   CT_Status(20)       ! Cobetrieve return status
	Character*30 Buf1(118)          ! Write buffer
        Character*7  Buf2(118)          ! Write buffer
	Character*78 Dash
	Character*1  Dashes(78) / 78 * '-' /
	Equivalence  ( Dash, Dashes(1) )
	Integer*2   Ix, Jx, Kx          ! Indices
	Logical*1   Exceed(118)	        ! Flags for exceeding extrema
	Character*3 Flag(118,2)         ! Exceeded field table	
	Character*1 Blanks(30) / 30 * ' ' /
	Character*30 Blank		! Blank fill characters
	Equivalence ( Blanks(1), Blank )
	Integer*2   Ib, Ie              ! Begin and end pointers

!	Set the function status to normal.

	FXT_Write_Extrema = %loc(FXT_Normal)

!       Write the FXT_Eng_Xtrm record in the archive.

	Call CT_Write_Arcv (, Lun_Xtrm, Xtrm_Rec, CT_Status )
	
	If ( CT_Status(1) .Ne. CTP_Normal ) Then
	  Status = CT_Status(1)
	  Call Lib$Signal (FXT_CTWritErr, %val(1), %val(Status))
	  FXT_Write_Extrema = %loc(FXT_Aberr)
	Endif

!	If a report is requested, log the minima and maxima exceeded.

	If((FXT_Write_Extrema .Eq. %loc(FXT_Normal)) .And. (Report)) Then
	  Write(Unit=Lun_Rpt, FMT=100, Iostat=Status) Dash
	  If (Status .Ne. Zero) Then
	    FXT_Write_Extrema =%loc(FXT_Aberr)
	    Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	  Endif
	  If (FXT_Write_Extrema .Eq. %loc(FXT_Normal)) Then
	    Do Ix = 1, 118
	      Exceed(Ix) = .False.
	      Do Jx = 1, 2
	        Flag(Ix,Jx) = '   '
		If (Xtrm_Table(Ix,Jx)) Then
	          Exceed(Ix) = .True.
	          Flag(Ix,Jx) = 'New'
	        Endif
	      Enddo
	    Enddo	
	    Jx = 0
	    Do Ix = 1, 118
	      Buf1(Ix) = Blank(1:30)
	      Buf2(Ix) = Blank(1:7)
	      If ( Exceed(Ix) ) Then  
	        Jx = Jx + 1
	        Buf1(Jx) = Plot_Name(Ix)
	        Buf2(Jx) = Flag(Ix,1) // ' ' // Flag(Ix,2)
	      Endif
	    Enddo
	    If ( Mod(Jx,2) .Eq. Zero ) Then
              Kx = Jx / 2
	    Else
	      Kx = (Jx / 2) + 1
	    Endif
	    Ib = 1
	    Ie = 2 
	    Do Ix = 1, Kx
	      Write(Unit=Lun_Rpt, FMT=200, Iostat=Status) 
	1	Buf1(Ib), Buf2(Ib), Buf1(Ie), Buf2(Ie)
	      Ib = Ib + 2
	      Ie = Ie + 2
	    Enddo
	    If (Status .Ne. Zero) Then
	      FXT_Write_Extrema = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	    Endif
	  Endif
	Endif   	! Report

 100    Format(1x/10x,'Table of Engineering Analogs with New ',
	1   'Extrema for Current Run'// 6x, 
	2   ' Engineering Analog ',6x,'Min Max',8x,
	3   ' Engineering Analog ',6x,'Min Max' / 1x, A)
 200    Format(1x,A,1x,A,3x,A,1x,A)

	Return
	End

