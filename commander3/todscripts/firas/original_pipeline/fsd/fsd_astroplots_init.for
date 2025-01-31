	Integer*4 Function FSD_Astroplots_Init(plt_device,
	2		plt_com,plt_com_file,file_seg,
	3               scitype,GMT_start,GMT_stop,numchans,sel_chan,outtype,
	4		instr_qual,attit_qual)

C----------------------------------------------------------------------------
C	
C	This subroutine reads the 5 'Earth radiation field' files and parses
C	   the command line.
C
C----------------------------------------------------------------------------
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C       version 4.2.1 11/24/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C       version 4.2.1 12/07/88, SPR 2898, QUOC CHUNG, STX
C       FIXED PROBLEM THAT DOES NOT ALLOW USER TO SELECT INDIVIDUAL
C       CHANNEL.
C
C	version 4.4  Enable plotting vs sky.  (SER 2373)
C	Fred Shuman, STX, 1989 Jan 26.
C
C	version 4.4  Add a /PLOTS qualifier.  (SPR 3967)
C	Fred Shuman, STX, 1989 Jun 15.
C
C	SPR 3945, 4178, Allow Astroplots to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Astroplots, _Init, _Catinfo,
C	    and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C       SPR 4133  Correct  options for /channel to be: ALL, RIGHT, LEFT,
C           HIGH, LOW, RH, RL, LH, LL now. Harte Wang, STX /1989 Oct. 17       
C
C       SPR 4142  ASTROPLOTS uses a simulated attitude if no attitude is
C	    contained in the raw science.  R. Kummerer, STX / 1989 Oct 24
C
C       SPR 4823  ASTROPLOTS signals EOF detection on all radiation field
C	    files.  R. Kummerer, STX / 1989 Oct 25
C
C	SPR 5752, Use character JSTART and JSTOP in call to FUT_ATTITUDE.
C	    R. Kummerer, STX / 1990 Jan 30
C
C	SPR 5777, Select instrument and attitude data quality.
C	    R. Kummerer, STX / 1990 Feb 2
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C       SPR 7171, Increase the value of the local variable TXTVAL to avoid
C           truncation of the qualifier "/PLOTDEVICE=FSD_AST_FDQ.QMS/QMS".
C           N. Gonzales, STX / 1990 Aug 01 
C
C       SPR 9847, Update FSD to use new FEX_Samprate reference file.
C           N. Gonzales, Hughes STX / 1992 July 30.
C----------------------------------------------------------------------------

	Implicit None

	Include '(FSD_Astroplots)'       ! includes record structure RadFld,
	                                 ! '(FUT_Params)', and '(UPM_Stat_Msg)'
	Include '(CCT_Status_Record)'
	Include 'CSDR$Library:CTUser.inc'
	Include '($SSDEF)'

	Record / CCT_Status / struc

C   Call list

	Character	*32	plt_device
	Integer		*4	plt_com
	Character	*64	plt_com_file
	Character	*39	file_seg	! File segment name
	Character	*14	GMT_start
	Character	*14	GMT_stop
	Integer		*4	JStart(2)
	Integer		*4	JStop(2)
	Integer		* 4	outtype		! Selected output type
	Integer		* 4	numchans	! Number of channels selected
	Integer		* 4	sel_chan(4)	! Array of selected channels
	Integer		* 4	instr_qual	! Instrument data quality
	Integer		* 4	attit_qual	! Attitude data quality

C   Local Declarations.

	Integer		* 4	status		! Return status
	Integer		* 4	plstatus	! Another return status
	Logical		* 1	plotpres	! Is /PLOTS qualfr. present?
	Integer		* 2	txtlen		! Length of text string
	Character	*64	txtval		! String to hold text
	Character	*32	infile 
	Character	* 3	scitype		! NFS, FPP, or FDQ
	Character	*64	plot_file
	Integer		* 4	cpy_lun
	Integer		* 4	i
	Integer		* 4	i1
	Integer		* 4	i2
	Character	* 3	str
	Real		* 4	lonprev
	Real		* 4	lon
	Real		* 4	lat
	Real		* 4	lontemp
	Real		* 4	lattemp
	Integer		* 4	reversals
	Integer		* 4	ios
	Integer		* 4	iunit
	Character	*39	Blank		! Blanks
	Character	* 1	Blanks(39) / 39 * ' ' /
	Equivalence  ( Blank, Blanks(1) )

C   Functions.

	Integer		* 4	FUT_SetDev
	Integer		* 4	OTS$Cvt_TU_L
	Integer		* 4	UPM_Present
	Integer		* 4	UPM_Get_Longword
	Integer		* 4	UPM_Get_Value
	Logical		* 1	Time_LT

C   External symbols

	External	FSD_Normal
	External	FSD_Fopen
	External	FSD_IllOrdFil
	External	FSD_FRerr
	External	FSD_Aberr1
	External	FSD_InvTim

	FSD_Astroplots_Init = %loc(FSD_Normal)
C
C  Read the 5 'Earth radiation field' files.  First, the S. Atlantic Anomaly:
C
	Call Lib$Get_Lun(iunit)
	infile = 'CSDR$MPD:MPD_SAA.DAT'

	Open (UNIT=iunit,    File=infile,     Form='FORMATTED',
	2     Readonly,      Status='OLD',    Recordtype='VARIABLE',
	3     Shared,        Iostat=ios)

	If (ios .Ne. 0) Then
	   Call Lib$signal(FSD_Fopen,%val(3),infile,%val(iunit),%val(ios))
	End If
c
c   The SAA file is a sequence of (long.,lat.) forming a closed loop.  It
c   starts at the min long. of the loop and goes forward along the North
c   boundary, then turns and follows the South boundary, ending at start pt.
c   We will split this into 2 arrays, SAA1=North & SAA2=South boundary.
c
	i1 = 0
	i2 = 0
	reversals = 0
	Read (IUNIT,*,iostat=ios) lon,lat
	lonprev = lon - 1.

	Do While (reversals .Lt. 2  .And.  ios .Eq. 0)
	   If (reversals .Eq. 0  .And.  lon .Ge. lonprev) Then
	      i1 = i1 + 1
	      RadFld.saa1_lon(i1) = lon
	      RadFld.saa1_lat(i1) = lat
	   Else If (reversals .Eq. 0  .And.  lon .Lt. lonprev) Then
	      reversals = 1
	      i2 = i2 + 1
	      RadFld.saa2_lon(i2) = RadFld.saa1_lon(i1)
	      RadFld.saa2_lat(i2) = RadFld.saa1_lat(i1)
	      i2 = i2 + 1
	      RadFld.saa2_lon(i2) = lon
	      RadFld.saa2_lat(i2) = lat
	   Else If (reversals .Eq. 1  .And.  lon .Lt. lonprev) Then
	      i2 = i2 + 1
	      RadFld.saa2_lon(i2) = lon
	      RadFld.saa2_lat(i2) = lat
	   Else
	      reversals = 2
	   End If

	   lonprev = lon
	   Read (IUNIT,*,iostat=ios) lon,lat
	End Do

	RadFld.saa1_num = i1
	RadFld.saa2_num = i2
c
c   Shouldn't see a second reversal in the longitude direction of travel.
c
	If (reversals .Eq. 2) Then
	   Call Lib$signal(FSD_IllOrdFil,%val(1),infile)
	End If

	If (ios .Lt. 0) Then
c
c   Shouldn't get an End-Of-File before closing the loop.
c
	   If (RadFld.saa1_lon(1) .Ne. RadFld.saa2_lon(i2) .Or.
	2      RadFld.saa1_lat(1) .Ne. RadFld.saa2_lat(i2)) Then
	      Call Lib$Signal(FSD_IllOrdFil,%val(1),infile)
	   End If
	Else If (ios .Gt. 0) Then
	   Call Lib$signal(FSD_FRerr,%val(2),%val(iunit),%val(ios))
	End If
c
c   Reverse the SAA South boundary to put it into increasing longitude order
c
	Do i=1,i2/2
	   lontemp = RadFld.saa2_lon(i)
	   lattemp = RadFld.saa2_lat(i)
	   RadFld.saa2_lon(i) = RadFld.saa2_lon(i2+1-i)
	   RadFld.saa2_lat(i) = RadFld.saa2_lat(i2+1-i)
	   RadFld.saa2_lon(i2+1-i) = lontemp
	   RadFld.saa2_lat(i2+1-i) = lattemp
	End Do

	Close(UNIT=iunit)
C
C  Now read the North Van Allen Belt (North boundary) file:
C
	infile = 'CSDR$MPD:MPD_VAB_NORTH1.DAT'

	Open (UNIT=iunit,    File=infile,     Form='FORMATTED',
	2     Readonly,      Status='OLD',    Recordtype='VARIABLE',
	3     Shared,        Iostat=ios)

	If (ios .Ne. 0) Then
	   Call Lib$signal(FSD_Fopen,%val(3),infile,%val(iunit),%val(ios))
	End If

	i=0
	Do While (ios .Eq. 0)
	   Read (IUNIT,*,iostat=ios) lon,lat
	   If (ios .Eq. 0) Then
	      i = i + 1
	      RadFld.vabn1_lon(i) = lon
	      RadFld.vabn1_lat(i) = lat
	   Else If (ios .Lt. 0) Then
	      Continue
	   Else
	      Call Lib$signal(FSD_FRerr,%val(2),%val(iunit),%val(ios))
	   End If
	End Do
	RadFld.vabn1_num = i
	Close(UNIT=iunit)
C
C  Now read the North Van Allen Belt (South boundary) file:
C
	infile = 'CSDR$MPD:MPD_VAB_NORTH2.DAT'

	Open (UNIT=iunit,    File=infile,     Form='FORMATTED',
	2     Readonly,      Status='OLD',    Recordtype='VARIABLE',
	3     Shared,        Iostat=ios)

	If (ios .Ne. 0) Then
	   Call Lib$signal(FSD_Fopen,%val(3),infile,%val(iunit),%val(ios))
	End If

	i=0
	Do While (ios .Eq. 0)
	   Read (IUNIT,*,iostat=ios) lon,lat
	   If (ios .Eq. 0) Then
	      i = i + 1
	      RadFld.vabn2_lon(i) = lon
	      RadFld.vabn2_lat(i) = lat
	   Else If (ios .Lt. 0) Then
	      Continue
	   Else
	      Call Lib$signal(FSD_FRerr,%val(2),%val(iunit),%val(ios))
	   End If
	End Do
	RadFld.vabn2_num = i
	Close(UNIT=iunit)
C
C  Now read the South Van Allen Belt (North boundary) file:
C
	infile = 'CSDR$MPD:MPD_VAB_SOUTH1.DAT'

	Open (UNIT=iunit,    File=infile,     Form='FORMATTED',
	2     Readonly,      Status='OLD',    Recordtype='VARIABLE',
	3     Shared,        Iostat=ios)

	If (ios .Ne. 0) Then
	   Call Lib$signal(FSD_Fopen,%val(3),infile,%val(iunit),%val(ios))
	End If

	i=0
	Do While (ios .Eq. 0)
	   Read (IUNIT,*,iostat=ios) lon,lat
	   If (ios .Eq. 0) Then
	      i = i + 1
	      RadFld.vabs1_lon(i) = lon
	      RadFld.vabs1_lat(i) = lat
	   Else If (ios .Lt. 0) Then
	      Continue
	   Else
	      Call Lib$signal(FSD_FRerr,%val(2),%val(iunit),%val(ios))
	   End If
	End Do
	RadFld.vabs1_num = i
	Close(UNIT=iunit)
C
C  Finally, read the South Van Allen Belt (South boundary) file:
C
	infile = 'CSDR$MPD:MPD_VAB_SOUTH2.DAT'

	Open (UNIT=iunit,    File=infile,     Form='FORMATTED',
	2     Readonly,      Status='OLD',    Recordtype='VARIABLE',
	3     Shared,        Iostat=ios)

	If (ios .Ne. 0) Then
	   Call Lib$signal(FSD_Fopen,%val(3),infile,%val(Iunit),%val(ios))
	End If

	i=0
	Do While (ios .Eq. 0)
	   Read (IUNIT,*,iostat=ios) lon,lat
	   If (ios .Eq. 0) Then
	      i = i + 1
	      RadFld.vabs2_lon(i) = lon
	      RadFld.vabs2_lat(i) = lat
	   Else If (ios .Lt. 0) Then
	      Continue
	   Else
	      Call Lib$signal(FSD_FRerr,%val(2),%val(Iunit),%val(ios))
	   End If
	End Do
	RadFld.vabs2_num = i
	Close(UNIT=iunit)
C
C Finished reading the 5 'Earth radiation field' files.
C
C  Now parse the command line.  First, are we Interactive or "Batch"?
C
	Status = UPM_Present('INTERACTIVE')
	If (Status .Eq. UPM_Negated) Then
	   interactive = fac_not_present
	Else
	   interactive = fac_present
	End If
C
C Get the plot device.
C
	If (interactive .Eq. fac_not_present) Then
	  plt_device = 'PLT_HARDCOPY'
	Else
	  plt_device = 'PLT_DEVICE'
	End If

	plstatus = UPM_Present('PLOTDEVICE')
	plotpres = plstatus .Eq. UPM_Pres
	If (plotpres) Then
	   status = UPM_Get_Value('PLOTDEVICE', txtval, txtlen)
	   plt_device = txtval
	End If
C
C  Get the PLT command file.
C
	plt_com = fac_not_present
	plstatus = UPM_Present('PLTFILE')
	plotpres = plstatus .Eq. UPM_Pres
	If (plotpres) Then
	   status = UPM_Get_Value('PLTFILE', txtval, txtlen)
	   plt_com_file = txtval
	   plt_com = fac_present
	End If
C
C  Get the file segment or the time range.
C
	File_Seg(1:39) = Blank(1:39)

	If (UPM_Present('INPUT')) Then
	   Status = UPM_Get_Value('INPUT', File_Seg, Txtlen)
	   If (Status .Ne. SS$_Normal) Then
	      FSD_Astroplots_Init = %loc(FSD_Aberr1)
	      Call Lib$Signal (%val(Status))
	   Endif
	Elseif (UPM_Present('JSTART')) Then
	   GMT_start = fac_jstart_default
	   Status = UPM_Get_Value('JSTART', Txtval, Txtlen)
	   If (Status .Eq. SS$_Normal) Then
	      GMT_start(1:Txtlen) = Txtval(1:Txtlen)
	   Else
	      FSD_Astroplots_Init = %loc(FSD_Aberr1)
	      Call Lib$Signal (%val(Status))
	   Endif
	   GMT_stop = fac_jstop_default
	   Status = UPM_Get_Value('JSTOP', Txtval, Txtlen)
	   If (Status .Eq. SS$_Normal) Then
	      GMT_stop(1:Txtlen) = Txtval(1:Txtlen)
	   Else
	      FSD_Astroplots_Init = %loc(FSD_Aberr1)
	      Call Lib$Signal (%val(Status))
	   Endif
	   Call CT_Gmt_to_Binary( GMT_start, JStart)
	   Call CT_Gmt_to_Binary( GMT_stop, JStop)
	   If (Time_LT(JStop, JStart)) Then
	      FSD_Astroplots_Init = %loc(FSD_InvTim)
	      Call Lib$Signal (FSD_InvTim)
	   Endif
	Endif
C
C  Find out which type of science data was requested.
C
	status = UPM_Get_Value('SCIENCE', txtval, txtlen)
	scitype = txtval
C
C  Get the channel selection.
C
	status = UPM_Get_Value('CHANNEL', Txtval, Txtlen)
	If (Status .Ne. SS$_Normal) Then
	   FSD_Astroplots_Init = %loc(FSD_Aberr1)
	   Call Lib$Signal (%val(Status))
	Endif
	numchans = 0
        if (txtval(1:2) .eq. 'RH' .or. txtval(1:2) .eq. 'RL' .OR.
     *    txtval(1:2) .eq. 'LH' .or. txtval(1:2) .eq. 'LL' .or.
     *    txtval(1:3) .eq. 'ALL') then 
     	Do i=1,4
	   If (txtval(1:2) .Eq. fac_channel_ids(i)) Then
	      numchans = 1
	      sel_chan(1) = i
	   Else If (txtval(1:1) .Eq. 'A') Then
	      numchans = 4
	      sel_chan(i) = i
	   End If	   
	End Do
        ELSE
          if (txtval(1:5) .eq. 'RIGHT') then
             sel_CHAN(1) = 1
             sel_chan(2) = 2
             numchans = 2
          endif
          if (txtval(1:4) .eq. 'LEFT') then
             sel_CHAN(1) = 3
             sel_chan(2) = 4
             numchans = 2
          endif
          if (txtval(1:4) .eq. 'HIGH') then
             sel_CHAN(1) = 1
             sel_chan(2) = 3
             numchans = 2
          endif
          if (txtval(1:3) .eq. 'LOW') then
             sel_CHAN(1) = 2
             sel_chan(2) = 4
             numchans = 2
          endif
  
        Endif  
C
C Get the attitude solution.
C
        user_att_soln = fac_simulated

        If (UPM_Present ( 'ATTITUDE' ) .Eq. UPM_Pres) Then

           If (UPM_Present('ATTITUDE.SIMULATED') .Eq. UPM_Pres) Then

	      user_att_soln = fac_simulated

           Else If (UPM_Present('ATTITUDE.PREDICTED') .Eq. UPM_Pres) Then

	      user_att_soln = fac_predicted

           Else If (UPM_Present('ATTITUDE.COARSE') .Eq. UPM_Pres) Then

	      user_att_soln = fac_coarse

	   Else If (UPM_Present('ATTITUDE.FINE') .Eq. UPM_Pres .Or.
	2           UPM_Present('ATTITUDE.FINE') .Eq. UPM_Defaulted) Then

	      user_att_soln = fac_fine_with_dirbe

	      If (UPM_Present('ATTITUDE.FINE.DIRBE') .Eq. UPM_Pres) Then
	         user_att_soln = fac_fine_with_dirbe
	      Else If (UPM_Present('ATTITUDE.FINE.WITHOUT_DIRBE') .Eq.
	2              UPM_Pres) Then
	         user_att_soln = fac_fine_without_dirbe
	      End If

           Else If (UPM_Present('ATTITUDE.DEFINITIVE') .Eq. UPM_Pres) Then

	      user_att_soln = fac_definitive

	   End If

        End If
C
C Retrieve /QUALITY, determine acceptable data quality
C
	Instr_Qual = fac_many_yellow

	If (UPM_Present('QUALITY.INSTRUMENT')) Then

	  status = UPM_Get_LongWord ( 'QUALITY.INSTRUMENT', Instr_Qual )
          If (status .Ne. SS$_Normal) Then
	    FSD_Astroplots_Init = %Loc(FSD_Aberr1)
	    Call Lib$Signal (%Val(Status))
          Endif

	Endif

	Attit_Qual = fac_many_yellow

	If (UPM_Present('QUALITY.ATTITUDE')) Then

	  status = UPM_Get_LongWord ( 'QUALITY.ATTITUDE', Attit_Qual )
          If (status .Ne. SS$_Normal) Then
	    FSD_Astroplots_Init = %Loc(FSD_Aberr1)
	    Call Lib$Signal (%Val(Status))
          Endif

	Endif
C
C  And finally, what type of output:  Time plot (outtype=1) or Skymap (=2)?
C
	status = UPM_Get_Value('OUTPUT', Txtval, Txtlen)
	If (txtval(1:1) .Eq. 'T') Then
	   outtype = 1
	Else If (txtval(1:1) .Eq. 'S') Then
	   outtype = 2
	End If

	Return
	End
