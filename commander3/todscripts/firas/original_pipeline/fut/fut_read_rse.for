      Integer*4 Function FUT_Read_RSE(rse_file,rse)

C------------------------------------------------------------------------
C    PURPOSE: Read the specified archived RSE.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            May 16, 1989
C
C    INVOCATION: status = FUT_Read_RSE ( rse_file, rse )
C
C    INPUT PARAMETERS:
C
C	RSE_FILE	CH*64		Name of the RSE file.
C
C    OUTPUT PARAMETERS: 
C
C	STATUS		I*4		Success status.
C	RSE(16)		CH*128		Record selection expression.
C
C    SUBROUTINES CALLED: 
C	CT_Connect_Read
C	LIB$Get_LUN
C	LIB$Free_LUN
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	CTUser.Inc
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit None

	Include	'($SSDef)'
	Include 'CT$LIBRARY:CTUser.Inc'
	Include '(FUT_Params)'

	Integer		*4	ct_lun
	Character	*64	rse_file
	Character	*128	rse(16)
	Character	*128	crse
	Byte			irse(128)
	Integer		*4	status
	Integer		*4	tstatus
	Integer		*2	ct_status(20)
	Integer		*4	i

	Integer		*4	CT_Connect_Read
	External		CT_Connect_Read
	Integer		*4	LIB$Get_LUN
	Integer		*4	LIB$Free_LUN

	External	FUT_Normal
	External	FUT_CTOpenRSE
	External	FUT_CTReadRSE
	External	FUT_CTCloseRSE

	Equivalence (irse, crse)

C
C Fetch the LUN and open the RSE archive file.
C
	status = %Loc(FUT_Normal)

	tstatus = LIB$Get_Lun ( ct_lun )

	If (tstatus .Eq. SS$_Normal) Then

	   Open ( Unit=ct_lun, File=rse_file, Status='Old', 
     .		  IOStat=tstatus, UserOpen=CT_Connect_Read )

	   If (tstatus .Eq. 0) Then
C
C Read the RSE.
C
	      ct_status(1) = CTP_Normal

	      Do i=1,16
		 If (ct_status(1) .Eq. CTP_Normal) Then
	            Call CT_Read_Arcv( , ct_lun, irse, ct_status)
		    rse(i) = crse
		 End If
	      End Do

	      If (ct_status(1) .Ne. CTP_Normal) Then
		 status = %Loc(FUT_CTReadRSE)
	         Call LIB$Signal(FUT_CTReadRSE,%Val(1),%Val(ct_status(1)))
	      End If

C
C Close the archive file.
C
	      Call CT_Close_Arcv ( , ct_lun, ct_status )

	      If (ct_status(1) .Ne. CTP_Normal) Then
	         Call LIB$Signal(FUT_CTCloseRSE,%Val(1),%Val(ct_status(1)))
	      End If

	      tstatus = LIB$Free_LUN ( ct_lun )

	      If (tstatus .Ne. SS$_Normal) Then
		 Call LIB$Signal ( %Val(tstatus) )
	      End If

	   Else

	      tstatus = LIB$Free_LUN ( ct_lun )

	      If (tstatus .Ne. SS$_Normal) Then
		 Call LIB$Signal ( %Val(tstatus) )
	      End If

	      status = %Loc(FUT_CTOpenRSE)
	      Call LIB$Signal ( FUT_CTOpenRSE, %Val(1), %Val(tstatus) )
	      Return

	   End If

	Else

	   status = tstatus
	   Call LIB$Signal ( %Val(tstatus) )

	End If

	FUT_Read_RSE = status

	Return
	End
