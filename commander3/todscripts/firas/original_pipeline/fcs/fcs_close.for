	Integer*4  Function  FCS_Close  ( lun_in, snum, report, lun_rpt,
	1                                 lun_out, skymap, fileout, numpix,
	2                                 status, ok )
c-----------------------------------------------------------------------------
c Function to close input, output, and report files.
c Author: Larry P. Rosen, Hughes STX, March 1992.
c-----------------------------------------------------------------------------
	Implicit None

c Include
	Include		'($ssdef)'
	Include		'(fcs_msg)'
	Include		'(fut_error)'
	Include		'(fut_params)'

c Functions
	Integer*4	CSA_Write_Pixels, CSA_Close_Skymap

c Passed Parameters
	Integer*4	lun_out, lun_in (fac_max_num), lun_rpt
	Integer*2	numpix, snum
	Logical*1	report
	Character*56	skymap (fac_max_num)		! Input skymaps
	Character*60	fileout
	Integer*4	status, ok			! Processing status

c External
	external	csa_normal

c Local
	Integer*4	rstat, i
c-----------------------------------------------------------------------------
c Begin

	FCS_Close = %loc (fcs_normal)
	If (numpix .EQ. 0) Then
	   Call Lib$Signal (fcs_zerospec)
	Else
	   If (report) Then
	      Write (lun_rpt,10) numpix
   10	      Format (1X,'The number of pixels with data = ',I4)
	   Else
	      Write (6,10) numpix
	   Endif
	Endif
	If (status .EQ. ok) Then
	   Do i = 1,snum
	      rstat = CSA_Close_Skymap (lun_in(i), fac_skymap_no_levels)
	      If (rstat .NE. %loc (csa_normal)) Then
	         FCS_Close = %loc (fcs_abort)
	         Call Lib$Signal (fcs_csaclose, %val(2), skymap(i), %val(rstat))
	      Endif
	   Enddo
	   rstat = CSA_Close_Skymap (lun_out, fac_skymap_no_levels)
	   If (rstat .NE. %loc (csa_normal)) Then
	      FCS_Close = %loc (fcs_abort)
	      Call Lib$Signal (fcs_csaclose, %val(2), fileout, %val(rstat))
	   Endif
	Endif

c Signal the processing status.

	If (status .EQ. ok .AND. fcs_close .EQ. %loc (fcs_normal) ) Then
	   Call Lib$Signal (fcs_normal)
	Else
	   Call Lib$Signal (fcs_abort)
	Endif

c Close report

	If (report) Then
	   Close (Unit=lun_rpt, Iostat=rstat)
	   fut_report_lun = 0
	   If (rstat .NE. 0) Then
	      FCS_Close = %loc (fcs_abort)
	      Call Lib$Signal (fcs_closerep, %val(1), %val(rstat))
	   Endif
	Endif
	Return
	End
