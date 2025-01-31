	PROGRAM  FRD_CMDGAIN
C*******************************************************************************
C	Purpose:  To read the text file, FEX_CMDGAIN.TXT, containing the
C       commanded gain values for the four channels and send them
C	to a reference file, FEX_CMDGAIN.DAT.
C
C	Author:  John Sims, STX, 11 June 1991, SER 7985
C	Calling Sequence:  This is a main program.
C	Input File:	FEX_CMDGAIN.TXT
C	Output File:	FEX_CMDGAIN.DAT
C*******************************************************************************
	Implicit None

	Dictionary	'FEX_CMDGAIN'
	Record /FEX_CMDGAIN/	FEX_Rec
	Integer*4	Lun
	Character*17	In_File  / 'FEX_CMDGAIN.TXT;' /
	Integer*4	Istat
	External	FRD_RMSOpen
	Integer*2	I
	Character*17	Out_File / 'FEX_CMDGAIN.DAT;' /
	Logical*1	Next /.True./
	Character*20	Text
        Real*4          in_val1,in_val2,in_val3,in_val4

	Call Lib$Get_Lun (Lun)
	Open (UNIT=Lun, FILE=In_File, STATUS='OLD', FORM='FORMATTED',
	1   IOSTAT=Istat)
	If (Istat) Then
	   Call Lib$Signal (FRD_RMSOpen, %Val(2), In_File, %Val(Istat))
	Else
	   Do While (Next)
	      Read (Lun,15) Text
	      If (Text(1:1).NE.'C') Next = .False.
	   EndDo
  15	   Format(A20)
      	   Call SYS$GetTim (FEX_Rec.CT_Head.Time)
 	   Call CT_Binary_to_GMT (FEX_Rec.CT_Head.Time, FEX_Rec.CT_Head.GMT)
           do I = 1,8
              Read(Lun,*) in_val1,in_val2,in_val3,in_val4
              FEX_rec.chan(1).cmdgain(I) = in_val1  
              FEX_rec.chan(2).cmdgain(I) = in_val2 
              FEX_rec.chan(3).cmdgain(I) = in_val3 
              FEX_rec.chan(4).cmdgain(I) = in_val4                           
	   enddo
	   Close (Lun)
	   Open (UNIT=Lun, FILE=Out_File, STATUS='NEW', FORM='UNFORMATTED',
	1     ACCESS='SEQUENTIAL', RECORDSIZE=64, ORGANIZATION='SEQUENTIAL',
	2     RECORDTYPE='FIXED', IOSTAT=Istat, SHARED)
	   If (Istat) Then
	      Call Lib$Signal (FRD_RMSOpen, %Val(2), Out_File, %Val(Istat))
	   Else
	      Write (Lun) FEX_Rec
	      Close (Lun)
	   EndIf
	EndIf
	End
