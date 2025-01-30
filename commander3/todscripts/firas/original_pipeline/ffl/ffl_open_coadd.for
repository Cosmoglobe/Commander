	Integer * 4 Function  FFL_Open_Coadd (ct_lun1, ct_lun2, vs_lun, ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_OPEN_COADD
c
c	This function opens the coadd archives CSDR$FIRAS_IN1 and
c	CSDR$FIRAS_IN2 for cal data.  It then opens the FISH-format spectrum
c	RMS file.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 March 1993
c		 SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c	    ffli record		invoc Structure (defined in ffl_invoc.txt)
c
c	Output:
c	    ct_lun1		Integer * 4	coadd CT file lun
c	    ct_lun2		Integer * 4	coadd CT file lun
c	    vs_lun		Integer * 4	voltage spectra RMS file lun
c	    ffli record		invoc Structure (defined in ffl_invoc.txt)
c	(NOTE: ffl_spec.txt is included solely to define parameter vrecsize,
c	    which is used in opening the FFL outfile. -- FGDS 1995 Jun 30)
c
c	Subroutines called:
c		fut_get_lun
c		ct_connect_read
c		LIB$Signal
c
c	Include files:
c		fut_params.txt
c		ffl_invoc.txt
c		ffl_spec.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c	Name of facility changed from FFI to FFL.
c	Fred Shuman, HSTX, 1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c	Subsumed the FFL_SPEC.TXT include file into FFL_CONFIG.TXT .
c	Fred Shuman, HSTX, 1995 July 21.
c-------------------------------------------------------------------------------

	Implicit none

	Include '(fut_params)'
	Include '(ffl_invoc)'		! defines record ffli (structure invoc)
	Include '(ffl_config)'		! defines record structure fflc and
		!  certain other vbles; Includes CCT_Get_Config; Dictionary
		!  structures FEX_GRTCOAWT, FEX_GRTTRANS.
		!defines structure fishspec and parameters vrecsize and vspecsiz
c
c  Call arguments:
c
	Integer * 4	ct_lun1		!  CT logical unit number
	Integer * 4	ct_lun2		!  CT logical unit number
	Integer * 4	vs_lun		!  voltage spectra RMS file lun
c
c  All other variables and functions:
c
	Integer * 4	io_stat		!  I/O return status
	Integer * 4	rstatus		!  return status
	Integer * 4	status		!  return status

	Integer * 4	fut_get_lun
	Integer * 4	ct_connect_read
	External	ct_connect_read

	External	ffl_ctifgopen
	External	ffl_normal
	External	ffl_rmsopen
	External	fut_normal

	status = %Loc(ffl_normal)
CC
CC  Open the primary coadd file.
c  Get the logical unit number for Cobetrieve.
c
	rstatus = fut_get_lun (ct_lun1)
	If (rstatus .Ne. %Loc(fut_normal)) Then
	   Call LIB$Signal (%Val(rstatus))
	EndIf
c
c  Determine the coadd file name.
c
	ffli.infile1 = 'CSDR$FIRAS_IN1:FIL_CAL_' // fac_channel_ids(ffli.chan)
     &               // fac_scan_mode_idsl(ffli.smode) // '/' // ffli.time_range
	Call STR$Trim (ffli.infile1, ffli.infile1, ffli.inlen1)
c
c  Open the files.
c
	Open (unit=ct_lun1, file=ffli.infile1, iostat=io_stat, status='old',
     &	      useropen=ct_connect_read)
	If (io_stat .Ne. 0) Then
	   status = %Loc(ffl_ctifgopen)
	   Call LIB$Signal (ffl_ctifgopen, %Val(2), ffli.infile1(1:ffli.inlen1),
     &					   %Val(io_stat))
	EndIf


	If ((ffli.hybrid .Eq. fac_present)  .And.
     &      (status .Eq. %Loc(ffl_normal))) Then
CC
CC  Open the secondary coadd file for hybrid input.
c  Get the logical unit number for Cobetrieve.
c
	   rstatus = fut_get_lun (ct_lun2)
	   If (rstatus .Ne. %Loc(fut_normal)) Then
	      Call LIB$Signal (%Val(rstatus))
	   EndIf
c
c  Determine the 2nd coadd file name.
c
	   ffli.infile2 = 'CSDR$FIRAS_IN2:FIL_CAL_'// fac_channel_ids(ffli.chan)
     &               // fac_scan_mode_idsl(ffli.smode) // '/' // ffli.time_range
	   Call STR$Trim (ffli.infile2, ffli.infile2, ffli.inlen2)
c
c  Open the files.
c
	   Open (unit=ct_lun2, file=ffli.infile2, iostat=io_stat, status='old',
     &	         useropen=ct_connect_read)
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_ctifgopen)
	      Call LIB$Signal (ffl_ctifgopen, %Val(2),
     &			       ffli.infile2(1:ffli.inlen2), %Val(io_stat))
	   EndIf

	EndIf


	If (status .Eq. %Loc(ffl_normal)) Then
CC
CC  Open the FISH-format RMS spectrum file.
c  Get the logical unit number for the RMS file.
c
	   rstatus = fut_get_lun (vs_lun)
	   If (rstatus .Ne. %Loc(fut_normal)) Then
	      Call LIB$Signal (%Val(rstatus))
	   EndIf
c
c  Determine the RMS file name.
c
	   ffli.outfile = 'CSDR$FIRAS_OUT:FFL_' // ffli.scan_mode // '.'
     &					        // ffli.file_ext
	   Call STR$Trim (ffli.outfile, ffli.outfile, ffli.outlen)
c
c  Open the file.
c
	   Open (unit=vs_lun, file=ffli.outfile, access='sequential',
     &		 status='new', form='unformatted', recordtype='fixed',
     &		 recl=vrecsize, iostat=io_stat)

	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_rmsopen)
	      Call LIB$Signal (ffl_rmsopen, %Val(2),
     &			       ffli.outfile(1:ffli.outlen), %Val(io_stat))
	   EndIf

	EndIf		!  (status from coadd open


	FFL_Open_Coadd = status

	Return
	End
