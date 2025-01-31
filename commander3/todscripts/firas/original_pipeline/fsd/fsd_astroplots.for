	Program FSD_ASTROPLOTS

C-----------------------------------------------------------------------------
C
C   PURPOSE:  To study peak height, glitch rate, & deglitcher noise
C	of FIRAS science data from the archive and
C	SAA, VAB, and galactic plane crossing, and moon approach angle
C
C   AUTHOR:
C	Jessy H. Dave'
C	ARC, September 1987
C
C   SUBROUTINES CALLED:
C	FSD_Astroplots_Init
C	FSD_Astroplots_Catinfo
C	FSD_Astroplots_Read
C	FSD_Astroplots_Data
C	FSD_Astroplots_Timeaxis
C	FSD_Astroplots_Param- - - - FSD_Astroplots_Sort
C	                          / FSD_Astroplots_Display
C	FSD_Astroplots_Plot - - -<- FSD_Astroplots_Label
C	                          \ FSD_Astroplots_Zoom
C	FSD_Astroplots_Skymap
C	Lib$Establish
C	CT_Binary_To_GMT
C	STR$Translate
C	LIB$MovC5
C	FUT_DeaDev
C	Uerase
C	LIB$Signal
C	Exit
C
C   COMMON VARIABLES USED:
C	/astro/ time_range,gmt,timeaxis,
C	        pk_ht,glitch_rate,degl_noise,saa,vab,galactic_cross,amoon,tmoon
C	/astro_data/ RadFld
C	/astro_num/ num,nsaa,nvab,ngalactic,nmoon
C	/flags/ begin,more
C	/qualifier/ interactive,in_file
C
C   INCLUDE FILES:
C	FSD_Astroplots.txt  . . . . Contains Common blocks,
C	                    Structure /Rad_Field_Data/, (Record RadFld),
C	                    Include FUT_Params,
C	                    Include UPM_Stat_Msg,
C	                            + a few general-purpose variables
C	FUT_Error
C	$SSDEF
C
C-----------------------------------------------------------------------------
C
C   Changes:
C
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
C       v4.2.2 Enable Astroplots to make skymap.  SPR 2373.  Fred Shuman, STX,
C       1989 Jan 31.
C
C	v4.4 Astroplots uses only simulated attitude & orbit.  If real ones
C	exist, use them instead.  SPR 3774.  Fred Shuman, STX.  1989 May 22.
C
C	v4.4 Suppress zooming (and prompting for it) in batch.
C	    SPR 3966.  Fred Shuman, STX.  1989 Jun 15.
C
C	v4.4 Add a /PLOTS qualifier; SPR 3967.  Fred Shuman, STX.  1989 Jun 15.
C
C	v4.4 Access violation following call to CCT_Query_TOD_Catalog.
C	    SPR 4007.  Fred Shuman, STX.  1989 Jun 15.
C
C	v4.4.1 SER 3306 Q. CHUNG STX 08/28/89
C	    Provide version number to track software updates.
C
C	SPR 3945, 4178, Allow Astroplots to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Astroplots, _Init, _Catinfo,
C	    and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C	SPR 4142, 4143, Correct determination of galactic plane crossings by
C	    providing an attitude solution and fetching this attitude
C	    information in FSD_ASTROPLOTS_DATA and pass it to
C	    FSD_ASTROPLOTS_PARAM.  R. Kummerer, STX / 1989 Oct 24.
C
C	SPR 4820, Correct checking of radiation fields.  R. Kummerer, STX /
C	    1989 Oct 25.
C
C	SPR 5085, Comment-out the call to CUT_Display_Banner in order to avoid
C	    conflict with FUT_SetDev, so that Astroplots can write a lineprinter
C	    plot file once again.  These lines should be reinstated when
C	    Astroplots is converted to PLT.  Fred Shuman, STX / 1989 Nov 17.
C
C	SPR 5197, Archive open failure due to using old interface, CT_OPEN_ARCV.
C	    CT_OPEN_ARCV does not know about archives FPP_SDF and FDQ_SDF.
C	    R. Kummerer, STX / 1990 January 30.
C
C	SPR 5777, Select instrument and attitude data quality.
C	    R. Kummerer, STX / 1990 Feb 2
C
C	SPR 6431, Increase the maximum data input limit and display the time
C	    range input.  R. Kummerer, STX / 1990 Mar 8
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C       SPR 9847, Update FSD to use new FEX_SAMPRATE reference file using
C           CCT_Get_Config.  N. Gonzales, Hughes STX / 1992 August 4.
C-----------------------------------------------------------------------------

	Implicit None

	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ sci_rec

	Include '(FSD_Astroplots)'
	Include '(FUT_Error)'
	Include '($SSDEF)'
	Include '(CCT_Get_Config)'
	Include 'CT$Library:CTuser.inc'

	Integer   *  4      CUT_Register_Version
	Integer   *  4      CUT_Display_Banner
	Integer   *  4      CCT_Open_Config
	Integer	  *  4	    CT_Init
	Integer   *  4      ct_stat(20)     ! cobetrieve status return array
	Integer   *  4      num_vol/80/
	Integer   *  4      lun_out/6/
	Integer   *  4      success /0/, err /1/
	Character *  6      version
	Parameter          (version='9.8')

	Integer   *  4      i
	Integer   *  4      j
	Character *  3      scitype          ! NFS, FPP, or FDQ
	Character * 32      plt_device
	Integer   *  4      plt_com
	Character * 64      plt_com_file
	Integer   *  4      cpy_lun
	Logical   *  4      c_flag
	Logical   *  4      s_flag
	Logical   *  4      eof
	Character * 39      file_seg
	Character *  4      str
	Character * 32      xlabel
	Character * 40      title
	Character * 80      label
	Integer   *  2      nbreak(50)
	Real      *  8      tstart
	Real      *  8      tend
	Real      *  8      deltat
	Real      *  8      t_unit
	Character * 14      GMT_start
	Character * 14      GMT_stop
	Character * 14      data_start
	Character * 14      data_stop
	Character * 20      extension     ! Skymap filename extension
	Integer   *  4      ich
	Integer   *  4      numchans
	Integer   *  4      sel_chan(4)
	Integer   *  4      chan
	Integer   *  4      r_status
	Integer   *  4      status        ! General Program Status
	Integer   *  4      outtype
	Character * 64      skyname
	Real      *  4      skydata(6144,3)
	Integer   *  4      numdata(6144)
	Integer   *  4      instr_qual
	Integer   *  4      attit_qual
c
c GET_CONFIG variables.
c
	Dictionary 'fex_samprate'

	Structure /config_data/
	   Record /fex_samprate/ fex_samprate
	Endstructure

	Record /config_data/ config
	Record /config_status/ stat(1)

	Character *  1	access_mode/' '/    ! data set access mode
	Integer   *  4  number/1/	    ! number of data sets
	Integer   *  4	size(1)/64/         ! size of data sets in bytes
	                                    !  (record length).
	Character * 32  name(1)		    ! names of data sets
	Character * 14  ref_jstart /'86001000000000'/
	Character * 14  ref_jstop /'99365235958990'/
	Integer   *  4  jstart(2)
	Integer   *  4  jstop(2)

	Data name(1)/'csdr$firas_ref:fex_samprate'/

	Integer   *  4	lun(1)		    ! logical unit numbers
	Integer   *  4	cindex(1)           ! initial cache pointers
	Logical   *  1	new_segment(1)	    ! flag for new segments
	Logical   *  1  first_time /.true./
	Integer   *  4	ncache/1/
	Integer   *  4	ref_count

	Integer   *  4  UPM_Present
	Integer   *  4  FSD_Astroplots_Init
	Integer   *  4  FSD_Astroplots_Catinfo
	Integer   *  4  FSD_Astroplots_Read
	Integer   *  4  FSD_Astroplots_Data
	Integer   *  4  FSD_Astroplots_Param
	Integer   *  4  FSD_Astroplots_Plot
	Integer   *  4  FSD_Astroplots_Skymap
	Integer   *  4  FUT_DeaDev
	Integer   *  4  OTS$Cvt_TU_L

	External        FUT_Error
	External        FSD_opnconfigerr
	External        FSD_Normal
	External        FSD_CTInit
	External        FSD_Aberr1

	Call Lib$Establish ( FUT_Error )
	
	r_status = success
	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out,num_vol,
	1                            'FIRAS Facility FSD_Astroplots' )
C
C Initialize Cobetrieve.
C
	If (r_status .eq. success) Then
	    Call CT_Init(ct_stat)
	    If (ct_stat(1) .Ne. ctp_normal) Then
	       Call Lib$signal(FSD_CTInit,%val(1), %val(ct_stat(1)))
	       r_status = err
	    End If
	End If

	call ct_gmt_to_binary (ref_jstart,jstart)
	call ct_gmt_to_binary (ref_jstop,jstop)
c
c Access Fex_Samprate reference dataset using cct_get_config.
c
	If (First_time) Then
	    status = cct_open_config (jstart, jstop, number, name,
	1   	                      size, access_mode, ncache, lun,
	2			      cindex, stat, ref_count)
	    first_time = .false.
	    If (.not. status) Then
	       call lib$signal(fsd_opnconfigerr,%val(1),%val(status))
	    End if
	End If

	status = FSD_Astroplots_Init(plt_device,plt_com,plt_com_file,
	1			  file_seg,scitype,
	2                         GMT_start,GMT_stop,numchans,sel_chan,outtype,
	3			  instr_qual,attit_qual)

	If (status .Eq. %Loc(FSD_Normal)) Then

	   status = FSD_Astroplots_Catinfo(file_seg,scitype,GMT_start,GMT_stop,
	2                                  numchans,sel_chan,extension)

	   time_range = GMT_start//';'//GMT_stop//';'
	   Call Str$Translate(time_range,time_range,'0',' ')
C
C   Clear the arrays
C
	   Call Lib$movC5(0,,0,3000*4,timeaxis)
	   Call Lib$movC5(0,,0,200*4,saa)
	   Call Lib$movC5(0,,0,200*4,vab)
	   Call Lib$movC5(0,,0,200*4,galactic_cross)
	   Call Lib$movC5(0,,0,200*4,tmoon)
	   Call Lib$movC5(0,,0,200*2,amoon)
	   nsaa=0
	   nvab=0
	   ngalactic=0
	   nmoon=0
c
c
c	Loop over current channel until further notice
c
	   c_flag = .False.
	   s_flag = .False.
	   begin = .True.
	   Do ich=1,numchans
	      Call Lib$movC5(0,,0,6144*3*4,skydata)
	      Call Lib$movC5(0,,0,6144*1*4,numdata)
	      chan = sel_chan(ich)
	      If (interactive .Eq. fac_present) Then
	         print*,'Processing data for channel ',fac_channel_ids(chan)
	      End If
	      eof = .False.
	      Call Lib$movC5(0,,0,3000*4,pk_ht)
	      Call Lib$movC5(0,,0,3000*4,glitch_rate)
	      Call Lib$movC5(0,,0,3000*4,degl_noise)
	      Call Lib$movC5(0,,0,3000*4,gal_lat)
	      Call Lib$movC5(0,,0,3000*4,moon_angle)
	      Call Lib$movC5(0,,0,3000*4,terr_lon)
	      Call Lib$movC5(0,,0,3000*4,terr_lat)

	      num = 0
	      more = .True.
	      Do While (more)
	         If (.Not. s_flag) Then
	            status = FSD_Astroplots_Read(scitype,chan,file_seg,
	2					 eof,sci_rec)
	         End If
	         If (.Not. eof) Then
	            status = FSD_Astroplots_Data(GMT_start,GMT_stop,
	2					 chan,sci_rec,tend,outtype,
	3					 instr_qual,attit_qual,
	4                                        c_flag,s_flag,skydata,numdata,
	5                                        number, size, lun, cindex,
	6                                        config, new_segment,stat)
	         End If
	      End Do  !(more)
	      If (s_flag) Then
	         If (numchans .Gt. 1) Then
	            status = FSD_Astroplots_Read(scitype,chan,file_seg,
	2					 eof,sci_rec)
	            begin=.True.
	            s_flag=.False.
	         Else
	            begin=.False.
	         End If
	      Else
	         begin=.True.
	      End If

	      Call CT_Binary_To_GMT ( gmt(1,1), data_start )
	      Call CT_Binary_To_GMT ( gmt(1,num), data_stop )

	      Type 100, fac_channel_ids(chan), data_start, data_stop
100	      Format (' Processing channel ', a, '  ', a, ' to ', a)

	      If (outtype .Eq. 1) Then
C
C   Output selection is a plot vs time.
C
	         If (num .Gt. 0) Then
	            Call FSD_Astroplots_Timeaxis(nbreak,xlabel,tstart,tend,
	2                                        deltat,t_unit,c_flag)
	            If (.Not. c_flag) Then
	               status = FSD_Astroplots_Param(tstart,tend,deltat,t_unit)
	               c_flag=.True.
	            End If

	            status = FSD_Astroplots_Plot(chan,plt_device,
	1					 plt_com,plt_com_file,
	2					 xlabel,
	3					 fac_att_soln(user_att_soln))
	         End If !num > 0

	      Else If (outtype .Eq. 2) Then
C
C   Output selection is a skymap.  Convert quantities to averages; make skymap.
C
	         Do i=1,6144
	            If (numdata(i) .Gt. 1) Then
	               Do j=1,3
	                  skydata(i,j) = skydata(i,j) / numdata(i)
	               End Do
	            End If
	         End Do
	         skyname = 'FSD_SKY_' // fac_channel_ids(chan) // '.' //
	2		   extension
	         status = FSD_Astroplots_Skymap(skydata,numdata,skyname)
	      End If

	   End Do      ! ich    (i.e., chan)

	   Call LIB$Erase_Page(1,1)

	End If
c
c *** exit program with message.
c
	If (status .Eq. %loc(FSD_Normal)) Then
	   Call LIB$Signal(FSD_Normal)
	   Call Exit(SS$_Normal)
	Else
	   Call LIB$Signal(FSD_Aberr1)
	   Call Exit(SS$_Abort)
	End If

	End
