        Integer*4 Function FSD_Astroplots_Data(GMT_start,GMT_stop,chan,sci_rec,
	2				  tend,outtype,instr_qual,attit_qual,
	3                                 c_flag,s_flag,skydata,numdata,number,
	4                                 size,lun,cindex,config,new_segment,
	5                                 stat)

C---------------------------------------------------------------------
C
C	This subroutine extracts data for astroplots from the already
C	filled science record.
C
C	INVOCATION:
C	    status = FSD_Astroplots_Data(chan,sci_rec,tend,outtype,c_flag,
C                                        s_flag,skydata,numdata)
C	INPUT:
C			   Ch *14	GMT_start
C			   Ch *14	GMT_stop
C	                    I * 4	chan
C	       Dictionary 'NFS_SDF'	sci_rec
C	                    R * 8	tend
C	                    I * 4	outtype, 1=plot vs Time, 2=Skymap
C	                    L * 4	c_flag
C
C	OUTPUT:
C	                    L * 4	more, Common/flags/
C	                    I * 4	num, Common/astro_num/
C	                   Ch *14	gmt(3000), Common/astro/
C	                    R * 4	pk_ht(3000), Common/astro/
C	                    R * 4	glitch_rate(3000), Common/astro/
C	                    R * 4	degl_noise(3000), Common/astro/
C			    R * 4       gal_lat(3000), Common/astro/
C			    R * 4	moon_angle(3000), Common/astro/
C			    R * 4	terr_lat(3000), Common/astro/
C			    R * 4	terr_lon(3000), Common/astro/
C	                    L * 4	s_flag
C	                    R * 4	skydata(6144,3)
C	                    I * 4	numdata(6144)
C
C	INCLUDE FILES:
C	    FSD_Astroplots.txt  . . . . Contains Common blocks,
C	                                Include FUT_Params,
C	                                Include UPM_Stat_Msg,
C	                                + a few general-purpose variables
C
C---------------------------------------------------------------------   
C        
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C       version 4.2.1 11/24/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C       
C       version 4.2.1 01/10/89, SPR 3109, R. KUMMERER, STX
C	INTERFACE CHANGE FOR FUT_FIND_IFG_CTR TO PASS NGROUP AND
C	MTM_SPEED.
C
C	SER 2373, F. Shuman, STX
C	Enable production of skymap.  1989 Jan 30.
C
C	SPR 3135, R. Kummerer, STX.  1989 Feb 8.
C	Remove dither from raw IFG prior to calling FUT_FIND_IFG_CTR.
C
C	SPR 4143, R. Kummerer, STX.  1989 Oct 24.
C	Save galactic latitude from science record for determining
C	galactic crossings.
C
C	SPR 4184, Use terrestrial pixel instead od celestial pixel number
C	    a primary key for pixelization. R. Kummerer, STX / 1990 Jan 30
C
C	SPR 5752, Use character JSTART and JSTOP in call to FUT_ATTITUDE.
C	    R. Kummerer, STX / 1990 Jan 30
C
C	SPR 5759, Initialize SAMP_RATE array so glitch rate is calculated
C	    correctly.  R. Kummerer, STX / 1990 Jan 30
C
C	SPR 5760, Prevent filling 1 based pixel number into ASTROPLOTS
C	    skymap.  Add one to terrestial pixel number. R. Kummerer, STX /
C	    1990 Jan 31.
C
C	SPR 5697, Remove dependency on MUT_YD_TO_YMD; reorganize time
C	    convertion.  R. Kummerer / STX, February 1, 1990.
C
C	SPR 5777, Select instrument and attitude data quality.
C	    R. Kummerer, STX / 1990 Feb 2
C
C	SPR 6406, Replace call to FUT_FIND_IFG_CTR with table of IFG
C	    peak positions.  R. Kummerer, STX / 1990 Mar 8
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C	SPR 5866, Change IFG sampling rate from 684.0 to 681.43.
C           N. Gonzales, STX / 1990 Dec. 10
C
C	SPR 9847, Update FSD to use new FEX_Samprate reference file using
C           cct_get_config.  N. Gonzales, Hughes STX / 1992 August 4.
C---------------------------------------------------------------------   

	Implicit None

	Include '($SSDEF)'
	Include '(FSD_Astroplots)'
	Include '(CCT_Get_Config)'
c
c Functions 
c
 	Integer * 4 CCT_Get_Config_Tod
c
c Records
c
	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ sci_rec

	Dictionary 'fex_samprate'
	Structure /config_data/
	   Record /fex_samprate/fex_samprate
	Endstructure

	Record /config_data/ config
	Record /config_status/ stat(1)

	Real      *  8      tend
	Integer   *  4      outtype
	Logical   *  4      c_flag
	Logical   *  4      s_flag
	Real      *  4      skydata(6144,3)
	Integer   *  4      numdata(6144)
	Character * 14      GMT_start
	Character * 14      GMT_stop
	Character * 14      start_time
	Integer   *  4      instr_qual
	Integer   *  4      attit_qual

	Integer   *  4      i
	Integer   *  4      gain

	Integer   *  4      fake
	Integer   *  4      gmt_time(2)
	Integer   *  4      jtime(2)

	Integer   *  4      ngroup
	Integer   *  4      mtm_speed
	Integer   *  4      mtm_length
	Integer   *  4      nsweep
	Integer   *  4      nglitch
	Integer   *  4      ndegl_noise
	Integer   *  4      status
	Integer   *  4      ictr
	Integer   *  4      pix
	Integer   *  4      chan
	Real	  *  4	    dither
	Real      *  8      t
	Real      *  4      rifg(512)
	Real      *  4      samp_rate(2) 

	Real      *  8      aut_adt2t68
	Integer   *  4      fut_attitude

	Integer * 4     number   	! number of data sets
	Integer * 4     size            ! size of data sets in bytes
	Integer * 4     lun(1)		! logical unit numbers
	Integer * 4     cindex(1)       ! initial cache pointers
	Logical * 1     new_segment(1)	! flag for new segments

	External FSD_Normal
	External FUT_Normal
	External FSD_FAERR
	External FSD_MAXREC
	External FSD_getconfigerr

C  Extract needed information from each raw science record for processing
C  and plotting.  Determine an attitude for each raw science time (if
C  necessary) and find the peak of each IFG also.

	gmt_time(1) = sci_rec.ct_head.time(1)
	gmt_time(2) = sci_rec.ct_head.time(2)
	if (c_flag)then
	  t = aut_adt2t68(gmt_time)
	  if (t .gt. tend)then
	    more=.false.
	    return
	  endif
	endif

	ngroup = sci_rec.sci_head.sc_head9
	mtm_speed = sci_rec.sci_head.mtm_speed
	mtm_length = sci_rec.sci_head.mtm_length
	nsweep = sci_rec.sci_head.sc_head11
	gain = fac_gains(sci_rec.sci_head.gain)
	nglitch = sci_rec.sci_head.sc_head21
	ndegl_noise=sci_rec.sci_head.sc_head15
	fake = sci_rec.dq_data.fake
	start_time = sci_rec.ct_head.gmt

C  Obtain individual IFGs.

	call lib$movc5(0,,0,2048,rifg)
	dither = 0.
	do i=1,512
	  rifg(i)=floati(sci_rec.ifg_data.ifg(i))
	  dither = dither + rifg(i)
	enddo
	dither = dither / 512.
	do i=1,512
	  rifg(i) = rifg(i) - dither  
	enddo

C  Obtain an attitude if one does not exist in the raw science record.
C  If FDQ_SDF raw science data are input, skip those IFGs flagged bad.

	status = %loc(FSD_Normal)

	If (sci_rec.dq_data.data_quality(110) .le. instr_qual .and.
	2   sci_rec.dq_data.data_quality(109) .le. attit_qual) Then

	   If (sci_rec.attitude.solution .Eq. fac_none) Then

	      status = fut_attitude(sci_rec,user_att_soln,GMT_start,GMT_stop)

	   End If

	   If (status) Then

	      Call ct_gmt_to_binary (start_time, jtime)
c
c Access FEX_Samprate reference dataset using CCT_Get_Config. 
c
	      status = cct_get_config_tod ( jtime, number, size,
	1                         lun, cindex, config, new_segment, stat )

	      if (.not. status) then
	         call lib$signal(fsd_getconfigerr,%val(1),%val(status))
	      endif

	      samp_rate(1) = config.fex_samprate.sampling_rate
	      ictr = fac_ifg_peak(mtm_length,mtm_speed,chan)

	      If (outtype .Eq. 1) Then

C  Will want to plot quantities vs time.

	          num=num+1
	          gmt(1,num)=gmt_time(1)
	          gmt(2,num)=gmt_time(2)
	          pk_ht(num)=-9999.
	          If (rifg(ictr) .Ne. 0.) pk_ht(num)=rifg(ictr)/gain
		  glitch_rate(num)=-9999.
		  If (ngroup .Gt. 0 .And. nsweep .Gt. 0) glitch_rate(num)=
	1		       nglitch*samp_rate(fake+1)/ngroup/nsweep/512.
		  degl_noise(num)=-9999.
		  If (ndegl_noise .Ne. 0) degl_noise(num)=
	1                      Floatj(ndegl_noise)/Floatj(gain)

	          gal_lat(num)=sci_rec.attitude.galactic_latitude * fac_att_conv
	          moon_angle(num)=sci_rec.attitude.moon_angle * fac_att_conv
	          terr_lat(num)=sci_rec.attitude.terr_latitude * fac_att_conv
	          terr_lon(num)=sci_rec.attitude.terr_longitude * fac_att_conv

	      Else If (outtype .Eq. 2) Then

C  Will want to make a skymap containing the quantities.

	           samp_rate(2) = 512.
	           num=num+1
	           gmt(1,num)=gmt_time(1)
	           gmt(2,num)=gmt_time(2)
	           pix = sci_rec.attitude.terr_pixel_no + 1
	           skydata(pix,1) = skydata(pix,1) + rifg(ictr)/gain
	           skydata(pix,2) = skydata(pix,2) +
	2                  nglitch*samp_rate(fake+1)/ngroup/nsweep/512.
	           skydata(pix,3) = skydata(pix,3) + 
	1                            Floatj(ndegl_noise)/Floatj(gain)
	           numdata(pix) = numdata(pix) + 1
!	gr = nglitch*samp_rate(fake+1)/ngroup/nsweep/512.
!	write (99,100) sci_rec.ct_head.gmt,' ',nglitch,' ',ngroup,' ',nsweep,
!	2		' ',fake,' ',sci_rec.dq_data.data_quality(110),
!	3		' ',gr
!100	format (a,5(a,i8),a,f15.9,/)
	      End If

	   Else
	      Call Lib$signal(fsd_FAERR,%val(1),%val(status))
	   End If

	   If (outtype .Eq. 1 .And. num .ge. 3000) Then

C  If we're plotting quantities vs time, limit the plot groups to 3000 points.

	     more=.false.
	     s_flag=.true.
	     If (interactive .eq. fac_present) Then
	        Call Lib$signal(fsd_maxrec,%val(1),%val(num))
	     End If
	   End If

	End If


	FSD_Astroplots_Data = %loc(FSD_Normal)

	Return
	End
