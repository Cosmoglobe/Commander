	PROGRAM  FRD_MTMSWEEP
C*******************************************************************************
C	Purpose:  To read the text file, FEX_MTMSWEEP.TXT, containing the
C	timing of the MTM sweep and flyback and to write data as a record
C	structure to a reference file, FEX_MTMSWEEP.DAT.
C
C	Author:  Larry P. Rosen, STX, 1 May 1991
C	Calling Sequence:  This is a main program.
C	Input File:	FEX_MTMSWEEP.TXT
C	Output File:	FEX_MTMSWEEP.DAT
C*******************************************************************************
	Implicit None

	Dictionary	'FEX_MTMSWEEP'
	Record /FEX_MTMSWEEP/	FEX_Rec
	Integer*4	Lun
	Character*17	In_File  / 'FEX_MTMSWEEP.TXT;' /
	Integer*4	Istat
	External	FRD_RMSOpen
	Integer*2	I
	Character*17	Out_File / 'FEX_MTMSWEEP.DAT;' /
	Logical*1	Next /.True./
	Character*20	Text

	Call Lib$Get_Lun (Lun)
	Open (UNIT=Lun, FILE=In_File, STATUS='OLD', FORM='FORMATTED',
	1   IOSTAT=Istat)
	If (Istat) Then
	   Call Lib$Signal (FRD_RMSOpen, %Val(2), In_File, %Val(Istat))
	Else
	   Read (Lun,10)
  10	   Format (1X,17(/))
	   Do While (Next)
	      Read (Lun,15) Text
	      If (Text(1:10).Eq.'Start Here') Next = .False.
	   EndDo
  15	   Format(A20)
	   Read (Lun,15) Text
	   Read (Lun,20) FEX_Rec.CT_Head.GMT
  20	   Format (A14,/)
	   Call CT_GMT_to_Binary (FEX_Rec.CT_Head.GMT, FEX_Rec.CT_Head.Time)
	   Read (Lun,30) (FEX_Rec.Total_Sweep_Flyback(I),I=1,4)
  30	   Format (4(1X,I9,/))
	   Read (Lun,30) (FEX_Rec.Flyback(I),I=1,4)
	   Close (Lun)
	   Open (UNIT=Lun, FILE=Out_File, STATUS='NEW', FORM='UNFORMATTED',
	1     ACCESS='SEQUENTIAL', RECORDSIZE=32, ORGANIZATION='SEQUENTIAL',
	2     RECORDTYPE='FIXED', IOSTAT=Istat, SHARED)
	   If (Istat) Then
	      Call Lib$Signal (FRD_RMSOpen, %Val(2), Out_File, %Val(Istat))
	   Else
	      Write (Lun) FEX_Rec
	      Close (Lun)
	   EndIf
	EndIf
	End
