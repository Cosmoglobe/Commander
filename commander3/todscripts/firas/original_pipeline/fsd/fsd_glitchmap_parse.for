	Subroutine FSD_Glitchmap_Parse (scitype, lun1)

C------------------------------------------------------------------------------
C
C    Routine to parse information from the command line.
C
C    Written by Reid Wilson, STX Inc.,  22-NOV-1987
C
C      This routine has been added to the program to bring it up to CSDR
C      standards and satisfy the interface as defined in the FSD_GLITCHMAP
C      Walkthrough given by Jessy Dave' of ARC.
C------------------------------------------------------------------------------
C
C       version 4.2.1 12/02/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C       version 4.2.1 12/21/88, SPR 3010, QUOC CHUNG, STX
C       FIXED ERROR NOT CONVERTING THE CHARACTER VALUE OF CHAN_ID TO
C       INTEGER VALUE PROPERLY.
C
C	SPR 3945, 4178, Allow glitchmap to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Glitchmap, _Parse,
C	    _Init, _QryRead, and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C	version 4.4.1 SER 4567 R. Kummerer STX 09/15/89
C		Use PLT instead of TEMPLATE for graphics.
C
C	version 4.4.1 SPR 4137 R. Kummerer STX 10/06/89
C		Remove /OPTION as is does not and never will do
C		anything.
C
C       SPR 4133, Correct the Glitchmap options for /CHANNEL to be: ALL
C           RIGHT, LEFT, HIGH, LOW, RH, RL, LH, LL. Harte Wang STX 10/17/89   
C
C	SER 5728, Add capability to pass in a PLT command file.
C		Steven Alexander STX 03/07/90
C
C------------------------------------------------------------------------------

	Implicit None

	Include '(FSD_Glitchmap)' !contains declarations needed
	Include '($SSDef)'

	Integer   *  4   CLI$Present
	Integer   *  4   CLI$Get_Value
	Integer   *  4   ichar        ! convert character to ASCII code
	Integer   *  4   char         ! convert ASCII code to character

	Character *  3   scitype
	Character *255   value
	Integer   *  4   ret_code
	Integer   *  4   length
	Character *255   generic_string
	Integer   *  4   rstatus
	Integer   *  4   status
	Integer   *  4   iostat
	Integer   *  4   lun1
	Integer   *  4   ios

	External         FSD_FOpen
	External         CLI$_Present
	External         CLI$_Negated

	Rstatus = 0

	batch = .Not. (CLI$Present('INTERACTIVE'))

	ifg_flag = CLI$Present('IFG')

	ret_code = CLI$Get_Value ('SCIENCE', scitype, length)

	plot = fac_present
	If (Batch) Then
	   plot_device = 'PLT_HARDCOPY'
	Else
	   plot_device = 'PLT_DEVICE'
	End If
	If (CLI$Present('PLOTDEVICE') .Eq. %Loc(CLI$_Present)) Then
	   status = CLI$Get_Value ('PLOTDEVICE', plot_device)
	Else If (CLI$Present('PLOTDEVICE') .Eq. %Loc(CLI$_Negated)) Then
	   plot = fac_not_present
	End If

	plt_com = fac_not_present
	if (cli$present('PLTFILE') .eq. %loc(cli$_present)) then
	   plt_com = fac_present
	   status = cli$get_value('PLTFILE',plt_com_file)
	end if

	ret_code = CLI$Get_Value ('WRITE', out_file, length)

	generic_string = 'ARCHIVE'
	If (out_file(1:length) .Eq. generic_string(1:length)) Then
	   out_file = 'FSD_ARCHIVE.dat'
	End If

	Call LIB$Get_Lun(lun1)
	Open (UNIT=lun1, NAME=out_file, FORM='formatted', STATUS='new',
	2     IOSTAT=ios, RECORDTYPE='variable')
	If (ios .Ne. 0) Then
	   Call Lib$Signal(Fsd_fopen, %Val(2), %Val(lun1), %Val(ios))
	   status = SS$_Abort
	   Call Exit(status)
	End If

	n_chan = 0
	ret_code = CLI$Get_Value ('CHANNEL', value)

	If (ret_code) Then
	   If (value .Eq. 'ALL') Then
	      Do i=1, 4
	         nchan(i) = i
	      End Do
	      n_chan = 4
	   Else If (value .Eq. 'RH') Then
	      chan_id = 1
	      n_chan = 1
	   Else If (value .Eq. 'RL') Then
	      chan_id = 2
	      n_chan = 1
	   Else If (value .Eq. 'LH') Then
	      chan_id = 3
	      n_chan = 1
	   Else If (value .Eq. 'LL') Then
	      chan_id = 4
	      n_chan = 1
	   Else If (value .Eq. 'RIGHT') Then
              nchan(1) = 1
              nchan(2) = 2	     
	      n_chan = 2
	   Else If (value .Eq. 'LEFT') Then
	      nchan(1) = 3
              nchan(2) = 4
	      n_chan = 2
	   Else If (value .Eq. 'HIGH') Then
	      nchan(1) = 1
              nchan(2) = 3
	      n_chan = 2
	   Else If (value .Eq. 'LOW') Then
	      nchan(1)= 2
              nchan(2) = 4
	      n_chan = 2
	   Else
	      Print *, '****> Wrong chan_id: ', value, ' -- Please try again <****'
	      status = SS$_Abort
	      Call Exit(status)
	   End If
	Else
	   status = SS$_Abort
	   Call Exit(status)
	End If

	n_scan = 0
	ret_code = CLI$Get_Value ('SCAN_MODE', value)
	If (ret_code) Then
	   If (value .Eq. 'ALL') Then
	      nscan(1) = 'SF'
	      nscan(2) = 'SS'
	      nscan(3) = 'LF'
	      nscan(4) = 'LS'
	      n_scan = 4
	   Else
	      Do While (ret_code .And. (n_scan .Lt. 4))
	         n_scan = n_scan + 1
	         nscan(n_scan) = value
	         ret_code = CLI$Get_Value ('SCAN_MODE', value)
	      End Do
	   End If
	Else
	   status = SS$_Abort
	   Call Exit(status)
	End If

	Return
	End
