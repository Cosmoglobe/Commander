	Subroutine FSD_Poserr_Select (fgs, ctr, ctr_int)

C-----------------------------------------------------------------------------
C
C ***	READS INTEFEROGRAMS FROM RMS FILE(MINIPIPE OUTPUT SELECT.*) OR
C	      FROM ARCHIVE.
C
C-----------------------------------------------------------------------------
C ***   MODIFIED:
C	    Shirley M. Read STX, July 1988
C	    In the IFG display loop the user was asked to select an IFG
C	    from the array to plot. However only the last IFG read from the
C	    file or archive was passed to the FUT_RDisplay instead of the
C	    requested IFG. This error was corrected for Build 4.1. Also,
C	    the science header word for number of mirror sweeps was not
C	    checked to see if the IFG had already been divided by the number
C	    of sweeps. This check is needed for data collected before 1987.
C
C	    R. Kummerer, August 8, 1988, SPR 2061. Scale input ifgs by gain.
C
C       version 4.2.1 12/01/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C       VERSION 4.2.1 01/10/89, SPR 3109, R. KUMMERER, STX
C	INTERFACE CHANGE FOR FUT_FIND_IFG_CTR TO PASS NGROUP AND
C	MTM_SPEED.
C
C       VERSION 4.2.1 01/12/89, SPR 3113, R. KUMMERER, STX
C	PREVENT DIVIDE BY ZERO IN CASE OF 0 CURVATURE IN PARABOLIC
C	INTERPOLATION WHEN FINDING THE IFG PEAK.
C
C       VERSION 4.2.1 02/08/89, SPR 3135, R. KUMMERER, STX
C	REMOVE DITHER FROM IFGS PRIOR TO FUT_FIND_IFG_CTR.
C
C	SPR 4178, Allow PosErr to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_PosErr, _Select, _Spectrum,
C	    _ZPDDisp, _PosDisp, and FSD.CLD.  Fred Shuman, STX / 1989 Sep 8.
C
C       Version 4.4.1 9/25/89 SER 4568, R. Kummerer, STX
C		Use PLT instead of TEMPLATE for graphics.
C
C       Version 4.4.1 9/27/89 SER 4636, R. Kummerer, STX
C		Keep POSERR from complaining about EOF.
C
C	SER 5728, Steven Alexander, STX, March 7, 1990.  Add capability
C		to pass in a PLT command file.
C
C       Version 5.8 03/08/90 SPR 6407, R. Kummerer, STX
C		Eliminate FUT_FIND_IFG_CTR call; use IFG peak table instead.
C-----------------------------------------------------------------------------

	Implicit None

	Include   '($SSDef)'
	Include   '(FUT_Params)'
	Include   '(FSD_Poserr)'
	Include   'CT$Library:CTUser.inc'

	Real      *  4      fgs(512, 100)
	Real      *  4      ctr(100)
	Real      *  4      ctr_int(100)

	Character * 40      in_file
	Character *  1      reply
	Character *  1      answer
	Character *  2      channame
	Character *  1      arch
	Character *255      filespec 
	Real      *  4      fy(3)
	Real      *  4      y(513)
	Real      *  4      yin(513)
	Integer   *  2      fake
	Integer   *  2      isweep
	Integer   *  2      igroup
	Character * 14      gmt_time
	Character *  2      scan_mode
	Integer   *  4      status
	Integer   *  4      adds

	Logical   *  4      selected
	Logical   *  4      i_plot
	Logical   *  4      bad
	Logical   *  4      eof
	Logical   *  4      more
	Logical   *  1      fdq_eng(1024)

	Integer   *  2      ios         !io status
	Integer   *  2      ctstat(20)  !CT return status
	Integer   *  2      ct_nr       !CT unit number

	Integer   *  4      ns
	Integer   *  4      numpts
	Integer   *  4      ndis
	Integer   *  4      ictr, jctr  !Peak center
	Integer   *  2      i, j, k, ik
	Integer   *  2      deltax
	Integer   *  2      sngr, snsw
	Real      *  4      dither      !IFG dither

	Integer   *  4      rstatus
	Integer   *  4      iunit
	Character * 14      start_time
	Character * 14      stop_time
	Character * 30      time_range
	Integer   *  2      txtlen                          ! Length of text string
	Character * 14      txtval                          ! String to hold text
	Character *  3      scitype                         ! NFS, FPP, or FDQ
	Integer   *  2      nbad(100)
	Character * 14      no_onboard_div /'87001000000000'/ !Time past which
	                                                    !no onboard division
	                                                    !of the IFG is done
	Integer   *  4      UPM_Get_Value
	Integer   *  4      LIB$Get_Lun
	Integer   *  4      LIB$Free_Lun

	External            FSD_FOpen
	External            FSD_FRErr
	External            FSD_FREnd
	External            FSD_MaxRec
	External            FSD_CTClos
	External            FSD_CTInit
	External            FSD_CTOpen
	External            FSD_CTRead
	External            FSD_CTREOF
	External            FUT_Normal
	External            CT_Connect_Read

	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ sci_rec

	num = 0
	Write (6,100)
100	Format (/,' Enter file type: RMS or [ARCHIVE] > ', $)
	Accept 1,arch
1	Format(a)
	Call STR$UpCase(arch, arch)

	If (arch .Eq. 'R') Then
	   Write (6,200)
200	   Format (' Enter input file name >', $)
	   Accept 1,in_file
C
C ***   OPEN data from RMS file
C
	   ios = 0
	   iunit = 3

	   Open(UNIT=iunit, FILE=in_file, FORM='unformatted', READONLY,
	2       STATUS='old', RECORDTYPE='variable', IOSTAT=ios)

	   If (ios .Ne. 0) Then
	      Call LIB$Signal(FSD_FOpen, %Val(2), %Val(iunit), %Val(ios))
	      status = ss$_abort
	      eof = .True.
	      Call Exit(status)
	   Else
	      eof = .False.
	   End If

	Else
C
C ***  Read data from the archive
C
	   If (init) Then
	      Call CT_Init(ctstat)
	      init = .False.
	      If (ctstat(1) .Ne. CTP_Normal) Then
	         rstatus = ctstat(1)
	         Call LIB$Signal( FSD_CTInit, %Val(1), %Val(rstatus) )
	         status = SS$_Abort
	         Call Exit(status)
	      End If
	   End If

	   Type 10
10	   Format(' Enter start time > ', $)
	   Accept 1,start_time

	   Type 20
20	   Format(' Enter stop time > ', $)
	   Accept 1,stop_time

	   time_range = start_time // ';' // stop_time // ';'
	   Call STR$Translate(time_range, time_range, '0', ' ')
C
C   Open the Archive file.
C
	   selected = .False.
	   Do While (.Not. selected)
	      Write (6,300)
300	      Format (' Enter channel ID:  RH, RL, LH, LL > ', $)
	      Accept 1,channame
	      Call STR$UpCase(channame, channame)
	      chan_id = 0
	      Do While (.Not. selected)
	         chan_id = chan_id + 1
	         If ( channame .Eq. fac_channel_ids(chan_id) ) Then
	            selected = .True.
	         End If
	      End Do
	   End Do
C
C  Find out which type of science data was requested.
C
	   status = UPM_Get_Value('SCIENCE', txtval, txtlen)
	   scitype = txtval

	   filespec = 'CSDR$FIRAS_ARCHIVE:' // scitype // '_SDF_' // channame //
	2             '/' // time_range

	   status = LIB$Get_Lun(ct_nr)

	   Open(UNIT=ct_nr, FILE=filespec, STATUS='old',
	2       IOSTAT=ios, USEROPEN=CT_Connect_Read)

	   If (ctstat(1) .Ne. CTP_Normal) Then
	      rstatus = ctstat(1)
	      Call LIB$Signal(FSD_CTOpen, %Val(2), %Val(chan_id), %Val(rstatus))
	      status = SS$_Abort
	      eof = .True.
	      Call Exit(status)
	   Else
	      more = .True.
	      eof = .False.
	   End If
	End If   ! (arch .Eq. 'R')

	Do While (.Not. eof)
	   If (arch .Eq. 'R') Then
	      Read (iunit, IOSTAT=ios) sci_rec
	      If (ios .Eq. 0) Then
	         more = .True.
	         Read (iunit, IOSTAT=ios) fdq_eng
	      End If
	      If (ios .Lt. 0) Then
	         eof = .True.
	         more = .False.
	         Call LIB$Signal(FSD_FREnd, %Val(2), %Val(IUNIT), %Val(IOS))
	      Else If (ios .Eq. 0) Then
	         more = .True.
	      Else
	         Call LIB$Signal(FSD_FRErr, %Val(2), %Val(IUNIT), %Val(IOS))
	         status = SS$_Abort
	         more = .False.
	         eof = .True.
	         Call Exit(status)
	      End If

	   Else
C
C   Read a record from the proper archive.
C
	      Call CT_Read_Arcv(, ct_nr, sci_rec, ctstat)

	      If (ctstat(1) .Eq. CTP_Normal)  Then
	         more = .True.
	      Else If (ctstat(1) .Eq. CTP_EndofFile)  Then
	         eof = .True.
	      Else
	         rstatus = ctstat(1)
	         Call LIB$Signal(FSD_CTREOF, %Val(1), %Val(rstatus))
	         eof = .True.
	      End If
	   End If    ! (arch .Eq. 'R')

	   If (more) Then

	      chan_id = sci_rec.sci_head.chan_id
	      ngroup = sci_rec.sci_head.sc_head9
	      speed = sci_rec.sci_head.mtm_speed
	      length = sci_rec.sci_head.mtm_length
	      upmode = sci_rec.sci_head.sc_head1a
	      xcal_pos = sci_rec.dq_data.xcal_pos
	      gain = fac_gains(sci_rec.sci_head.gain)
C
C   Check to see if IFG has been divided by number of mirror sweeps.
C
	      If ((sci_rec.ct_head.gmt .Gt. no_onboard_div) .And.
	2         (.Not. sci_rec.sci_head.sc_head23 )) Then
	         nsweep = sci_rec.sci_head.sc_head11
	      Else
	         nsweep = 1
	      End If

	      gmt_time = sci_rec.ct_head.gmt
	      fake = sci_rec.dq_data.fake

	      If (speed .Eq. 0) Then
	         scan_mode(2:2) = 'S'
	      Else
	         scan_mode(2:2) = 'F'
	      End If

	      If (length .Eq. 0) Then
	         scan_mode(1:1) = 'S'
	      Else
	         scan_mode(1:1) = 'L'
	      End If

	      If (num .Eq. 0) Then
	         scan = scan_mode
	         fakeit = fake
	      End If

	      If (scan .Ne. scan_mode) Then
	         Print *, 'Change of scan_mode ', gmt_time, scan_mode, '.'
	         Print *, 'Data discarded.'
	         more = .False.
	      End If

	      If (fake .Ne. fakeit) Then
	         Print *, 'Change of fake-it mode', gmt_time, fake, '.'
	         Print *, 'Data discarded.'
	         more = .False.
	      End If

	      If (more) Then
	         num = num + 1
	         ngr(num) = ngroup
	         nsw(num) = nsweep
	         gmt(num) = gmt_time
	         dither = 0.

	         Do i=1,512
	            y(i) = floati(sci_rec.ifg_data.ifg(i))
	            dither = dither + y(i)
	         End Do

	         dither = dither / 512.

	         Do i=1,512
	            y(i) = y(i) - dither
	         End Do

	         adds = ngroup
		 ictr = fac_ifg_peak(length,speed,chan_id)

	         Do i=1,512
	            y(i) = y(i) / fac_gains(sci_rec.sci_head.gain)
	         End Do

	         Call LIB$MovC3(512*4, y(1), fgs(1,num))

	         Call LIB$MovC3(3*4, y(ictr-1), fy(1))
	         If ( (fy(3)-fy(2)) * (fy(2)-fy(1)) .Lt. 0) Then
	            ctr(num) = ictr - .5*(fy(3)-fy(1)) /
	2                                 (fy(3) - 2.*fy(2) + fy(1))
	         Else
	            ctr(num) = ictr
	         End If
	         jctr = ctr(num)
	         Call LIB$MovC3(12, y(jctr), fy(1))
	         deltax = ctr(num) - jctr
	         ctr_int(num) = fy(1) + deltax*( (fy(2)-fy(1)) +
	2                            .5*(deltax-1)*(fy(3) - 2.*fy(2) + fy(1)) )

	         If (num .Eq. 100) Then
	            more = .False.
	            eof = .True.
	            Call LIB$Signal(FSd_MaxRec, %Val(1), %Val(num))
	         End If   ! (num .Eq. 100)

	      End If   ! (more)
	   End If      ! (more)
C
C   End of file
C
	   If (eof) Then
	      If (arch .Eq. 'R') Then
	         close(unit=IUNIT)
	      Else
	         Call CT_Close_Arcv(, ct_nr, ctstat)
	         If ( ctstat(1) .Ne. CTP_Normal) Then
	            rstatus = ctstat(1)
	            Call LIB$Signal(FSD_CTClos, %Val(1), %Val(rstatus))
	         End If
	      End If  ! (arch .Eq. 'R')
	   End If    ! (eof)

	End Do     ! eof

	If (num .Gt. 0) Then
	   Write (6,700) num
700	   Format (/, ' Number of IFGs: ', i3, '.')
	   i_plot = .True.
	   Do While (i_plot)
	      Write (6,400)
400	      Format (' Do you want to see an IFG? (Y/[N]) > ', $)
	      Accept 1, ans
	      Call STR$UpCase (ans, ans)
	      If (ans(1:1) .Eq. 'Y') Then
	         Write (6,500)
500		 Format (' Enter serial number of IFG to plot > ', $)
	         Accept *,ns
	         If ((ns .Le. num) .And. (ns .Ge. 1)) Then

	            Write (label,50) gmt(ns), ctr(ns), ctr_int(ns), ns
50	            Format(1x, a14, ' ', 'IFG Ctr=', f7.3, ' ', 'Peak Ht=',
	2		   f8.2, ' ', 'IFG Num=', i3)

	            Call LIB$MovC3(512*4, fgs(1,ns), yin(1))

	            igroup = ngr(ns)
	            isweep = nsw(ns)

		    startx = -1
		    spacex = 1
		    plot_label = 'Raw Science Interferogram'

		    Do i=1,512
		      difg(i) = yin(i)
		    Enddo

		    Call FUT_Plot_Title(label, ngroup, nsweep, speed,
	2				length, chan_id, upmode, xcal_pos,
	3				fakeit, gain, plot_label, ptitle)
		    Call FUT_Plot(difg, 512, startx, spacex, ptitle,
	2			  'Sample Number', 'Counts', zp,
	3			  interactive, plot_device, plt_com,
	4			  plt_com_file )

	         Else
	            Print *, 'No data for this case.'
	         End If     ! ((ns .Le. num) .And. (ns .Ge. 1))

	      Else          ! (answer(1:1) .Eq. 'Y')  [see IFG]
C
C ****  Discard unwanted IFGs
C
	         Type 40, num
40		 Format(1x, 'Do you want to discard any of these ', i3,
	2		' IFGs? (Y/[N]) > ', $)

	         Accept 1,ans
	         Call STR$UpCase(ans, ans)
	         If (ans .Eq. 'Y') Then
	            Call LIB$MovC5(0,, 0, 100, nbad)
	            i = 0
	            ndis = 1
	            Do While (ndis .Ne. 0)
	               i = i + 1
	               Write (6,600)
600		       Format (' Enter serial # of IFG to discard (0 for none) > ', $)
	               Accept *,ndis
	               If (ndis .Ne. 0) Then
	                  nbad(i) = ndis
	               End If
	            End Do
	            ndis = i - 1

	            j = 0
	            k = 0
	            Do While (k .Lt. num)
	               k = k + 1
	               bad = .True.
	               If (ndis .Ne. 0) Then
	                  Do ik=1,ndis
	                     If (k .Eq. nbad(ik)) Then
	                        bad = .False.
	                     End If
	                  End Do
	               End If

	               If (bad) Then
	                  j = j + 1
	                  Call LIB$MovC3 (512*4, fgs(1,k), fgs(1,j))
	                  ngr(j) = ngr(k)
	                  nsw(j) = nsw(k)
	                  ctr(j) = ctr(k)
	                  ctr_int(j) = ctr_int(k)
	               End If  !(Bad)

	            End Do ! (k < num)

	            num = j
	            Print *, 'Number of IFGs: ' , num, '.'
	         Else
	            I_plot = .False.
	         End If   ! (ans .Eq. 'Y')  [discard any ifgs]
	      End If      ! (ans .Eq. 'Y')  [see IFG]
	   End Do         ! (i_plot)

	   sngr = ngr(1)
	   snsw = nsw(1)
	   If (num .Ge. 2) Then
	      Do i=2,num
	         sngr = sngr + ngr(i)
	         snsw = snsw + nsw(i)
	      End Do
	      ngroup = sngr/num
	      nsweep = snsw/num
	   End If   ! (num .Ge. 2)
	End If      ! (num .Ge. 0)

	Return
	End
