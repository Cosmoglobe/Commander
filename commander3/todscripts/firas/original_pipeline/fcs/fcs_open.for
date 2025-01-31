	Integer*4  Function  FCS_Open  ( in_skymap, arch_in, fileout, arch_out,
	1                                chan, scan, fext, snum, tstart, tstop,
	2                                report, lun_rpt, goodmap, lun_in,
	3                                lun_out )

c-----------------------------------------------------------------------------
c Function to open input and output data sets.
c Additionally, it updates "snum" to be just the number of input skymap files
c that are "good" (right scan & chan), and updates the in_skymap, tstart, and
c tstop arrays to only contain the information for "good" skymaps.  Thus these
c arrays are effectively reduced to the length "snum".
c Author: Larry P. Rosen, Hughes STX, March 1992
c-----------------------------------------------------------------------------
	Implicit None

c Include

	Include		'($ssdef)'
	Include		'(fut_error)'
	Include		'(fut_params)'
	Include		'(fcs_msg)'

c Parameters

	Character*56	in_skymap (fac_max_num)		! Input skymaps
	Character*(*)	arch_in				! input data archive
	Character*60	fileout
	Character*(*)	arch_out			! output data archive
	Character*2	chan, scan			! channel, scan mode
	Character*20	fext				! Output file extension
	Integer*2	snum				! # of input skymaps
	Character*14	tstart (fac_max_num)		! start times
	Character*14	tstop (fac_max_num)		! stop times
	Logical*1	report				! Flag report file
	Integer*4	lun_rpt				! Report logical unit #
	Logical*1	goodmap (fac_max_num)	! map has right chan & scan
	Integer*4	lun_out, lun_in (fac_max_num)	! Logical units

c Function

	Integer*4	Lib$Get_Lun
	Integer*4	CSA_Field_Offset_Values

c External

	External	csa_normal
	External	csa_open_skymap
	External	csa_field_offset_values

c Local

	Integer*4	rstat
	Character*60	filein
	Integer*2	skymaplen		!  skymap length in longwords
	Integer*4	f_pix_offset	! I*4 of fac_coad_spec_pix_offset
	Integer*4	f_time_offset	! I*4 of fac_time_offset
	Integer*2	i, j /0/
	Integer*4	lun

c Begin

	FCS_Open = %loc (fcs_normal)
	rstat = Lib$Get_Lun (lun_out)
	If (rstat .NE. ss$_normal) Then
	   FCS_Open = %loc (fcs_abort)
	   Call Lib$Signal (fcs_lunerr, %val(1), %val(rstat))
	Else
	   fileout = arch_out // 'FCS_SKY_' // chan // scan // '.' // fext
	   skymaplen = fac_coad_spec_size / 4
	   f_pix_offset = fac_coad_spec_pix_offset
	   f_time_offset = fac_time_offset
	   Open ( Unit=lun_out, File=fileout, Status='NEW', form='UNFORMATTED',
	1         Recordtype='FIXED', Recl=skymaplen, Useropen=CSA_Open_Skymap,
	2         Iostat=rstat )
	   If (rstat .NE. 0) Then
	      FCS_Open = %loc (fcs_abort)
	      Call Lib$Signal (fcs_csaopen, %val(2), fileout, %val(rstat))
	   Else
	      rstat = CSA_Field_Offset_Values (f_pix_offset, f_time_offset, -1,
	1                                      lun_out)
	      If (rstat .NE. %loc (csa_normal)) Then
	         FCS_Open = %loc (fcs_abort)
	         Call Lib$Signal (fcs_csaoffset, %val(1), %val(rstat))
	      Endif
	   Endif
	Endif
	If (FCS_Open .EQ. %loc (fcs_normal)) Then
	   Do i = 1, snum
	      If (goodmap(i)) Then
	         rstat = Lib$Get_Lun (lun)
	         If (rstat .NE. ss$_normal) Then
	            FCS_Open = %loc (fcs_abort)
	            Call Lib$Signal (fcs_lunerr, %val(1), %val(rstat))
	         Else
	            filein = arch_in // in_skymap(i)
	            Open ( Unit=lun, File=filein, Status='OLD',
	1                  Form='UNFORMATTED', Recordtype='FIXED', Readonly,
	2                  Recl=skymaplen, Useropen=csa_open_skymap,
	3                  Iostat=rstat )
	            If (rstat .NE. 0) Then
	               FCS_Open = %loc (fcs_abort)
	               Call Lib$Signal (fcs_csaopen, %val(2),filein,%val(rstat))
	            Else
	               j = j + 1
	               lun_in(j) = lun
	               rstat = CSA_Field_Offset_Values (f_pix_offset,
	1                                               f_time_offset, -1,
	2                                               lun_in(j))
	               If (rstat .NE. %loc (csa_normal)) Then
	                  FCS_Open = %loc (fcs_abort)
	                  Call Lib$Signal (fcs_csaoffset, %val(1), %val(rstat))
	               Else

c Update in_skymap, tstart & tstop arrays to only be for "good" skymaps.

	                  in_skymap(j) = in_skymap(i)
	                  tstart(j) = tstart(i)
	                  tstop(j) = tstop(i)
	               Endif
	            Endif
	         Endif
	      Endif
	   Enddo
	Endif
	If (FCS_Open .EQ. %loc (fcs_normal)) Then
	   snum = j			! snum is now the number of good maps
	   If (report) Then
	      Write (lun_rpt,*) ' Number of good input skymaps = ', snum
	      Write (lun_rpt,*) ' Output file name = ', fileout
	   Endif
	Endif
	Return
	End
