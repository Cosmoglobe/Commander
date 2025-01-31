	Program FSD_Glitchmap
C      -----------------------
C
C------------------------------------------------------------------------------
C    Routine to study glitch rate of FIRAS science data from the archive.
C
C	Written by
C	Jessy H. Dave'
C	ARC, June 1987
C
C	Subroutines called:
C	   FSD_Glitchmap_Parse
C	   FSD_Glitchmap_Init
C	   FSD_Glitchmap_Qryread
C	   FSD_Glitchmap_Data
C	   FSD_Glitchmap_Lin_Disp - - - - FSD_Glitchmap_Zoom
C	   FSD_Glitchmap_Terr
C	   FSD_Glitchmap_Terr_Data
C	   CUT_Register_Version
C	   CUT_Display_Banner
C	   Lib$MovC5
C	   STR$UpCase
C	   LIB$Erase_Page
C	   Lib$Signal
C
C------------------------------------------------------------------------------
C       version 4.2.1 12/02/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C       version 4.2.1 12/21/88, SER 2575, QUOC CHUNG, STX
C       BRING FSD_GLITCHMAP TO USE QUERY BUILDER TO ACCESS RAW SCIENCE
C       BY USING RSE.
C
C	version 4.4.1 SER 3306 Q. CHUNG STX 08/28/89
C                     PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C
C	SPR 3945, 4178, Allow glitchmap to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  This module, FSD_Glitchmap_Parse,
C	    _Init, _QryRead, and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C	version 4.4.1 SER 4567 R. Kummerer STX 09/15/89
C		      Use PLT instead of TEMPLATE for graphics.
C
C       SPR 4133, Correct the Glitchmap options for /CHANNEL to be: ALL, RIGHT, 
C           LEFT,HIGH, LOW, RH, RL, LH, LL. Harte Wang, STX 10/17/89 
C
C       SPR 9847, Update FSD to use new FEX_SAMPRATE reference file using
C           CCT_Get_Config.  N. Gonzales, Hughes STX / 1992 August 4.
C------------------------------------------------------------------------------

	Implicit None

	Include '(FSD_Glitchmap)'
	Include '($SSDEF)'
	Include '(CCT_Get_Config)'
	Include 'CT$Library:CTuser.inc'
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

	Integer   *  4  CUT_Register_Version
	Integer   *  4  CUT_Display_Banner
	Integer   *  4  CCT_Open_Config
	Integer   *  4  CT_Init
	Integer   *  4  ct_stat(20)
	Integer   *  4  success /0/, err /1/
	Integer   *  4  num_vol/80/
	Integer   *  4  rstatus
	Integer   *  4  lun_out/6/
	Character *  5  version
	Parameter       (version='9.8')

	Character *  3  scitype
	Logical   *  4  i_stop
	Logical   *  4  i_chan
	Logical   *  4  s_flag

	Integer   *  4  status, r_status
	Integer   *  4  speed
	Integer   *  4  length
	Integer   *  4  ich
	Integer   *  4  lun1
	Integer   *  4  ct_nr
	Integer   *  2  grp_num

	Real      *  4  gf_ifg(512)
	Real      *  4  gf_hist(512)
	Real      *  4  grp_thre
	Real      *  4  g1, g2

	Character * 14  grp_start
	Character * 14  grp_stop
	Character * 14  terr_start
	Character * 14  terr_stop
	Character *128  rse(16)
	Character *  1  response, ask

	Integer   *  4  STR$UpCase

	External    FSD_Normal
	External    FSD_AbErr3
	External    FSD_CTInit
	External    FSD_opnconfigerr


	Rstatus = CUT_Register_Version(version)
	Rstatus = CUT_Display_Banner(lun_out, num_vol,
	2                            ' FIRAS Facility FSD_Glitchmap')

	init = .True.
	status = SS$_Normal
	r_status = success
C
C Parse routine added  22-NOV-1987 by REW to get command line info
C
	Call FSD_Glitchmap_Parse(scitype, lun1)
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

C
C **** Call the Query Builder to select a timerange and RSEs
C
 50	Continue
	If (.Not. batch) Then
	   Call FSD_Glitchmap_Init(ct_nr, rse)
	End If
C
C *** Loop through selected channel
C
 100	Continue
	Do ich=1,n_chan
	   If (n_chan .EQ. 4) Then
	      chan_id = ich
	   End If
           IF (N_chan .EQ. 2 ) Then
             Chan_id = nchan(ich)
           ENDIF
	   If (.Not. batch) Then
	      Print *, 'Processing channel ', fac_channel_ids(chan_id)
	   End If

	   i_chan = .True.
	   initial = .True.
	   i_stop = .False.
	   Do While (.Not. i_stop)
	      begin = .True.
	      eof = .False.
	      s_flag = .False.
	      grp_max = .False.
	      Do While (.Not. eof)
	         If (.Not. grp_max) Then
	            Call LIB$MovC5(0,, 0, 2048, gf_hist)
	            Call LIB$MovC5(0,, 0, 2048, gf_ifg)
	            grp_num = 0
	            grp_thre = 0.0
	         End If
	         Call LIB$MovC5(0,, 0, 2048, g_hist)
	         Call LIB$MovC5(0,, 0, 2048, g_ifg)

	         more = .True.
	         num = 0
	         Do While (more)
	            If (.Not. s_flag) Then
	                  Call FSD_Glitchmap_Qryread(ct_nr,rse,sci_rec,scitype,
	2                                         i_stop,s_flag)
	            End If
	            If (s_flag .And. i_stop .And. (.Not. more)) Then
	               Goto 150
	            End If
	            If (.Not. i_stop) Then
	               Call FSD_Glitchmap_Data(sci_rec,s_flag,number,
	1                                      size,lun,cindex,config,
	2                                      new_segment,stat)
	            End If

	         End Do        !(more)

	         If (num .Gt. 0) Then
	            If (.Not. batch) Then
	               Call FSD_Glitchmap_Lin_Disp(ans)
	            End If

	            If (i_write) Then
	               If (grp_num .Eq. 0) Then
	                  grp_start = gmt(1)
	               End If
	               grp_num = grp_num + num
	               g1 = (grp_num - num)/grp_num
	               g2 = num/grp_num
	               Do i=1,512
	                  gf_hist(i) = gf_hist(i) + g_hist(i)
	                  If (ifg_flag) Then
	                     gf_ifg(i) = gf_ifg(i)*g1 + g_ifg(i)*g2
	                  End If
	               End Do
	               Do i=1,num
	                  grp_thre = grp_thre + degl_thre(i)
	               End Do

	               If (.Not. grp_max) Then
	                  grp_stop = gmt(num)
	                  grp_thre = grp_thre/grp_num

	                  If (scan(1:1) .Eq. 'S') Then
	                     length = 0
	                  Else
	                     length = 1
	                  End If

	                  If (scan(2:2) .Eq. 'S') Then
	                     speed = 0
	                  Else
	                     speed = 1
	                  End If
	                  Write(lun1,*) chan_id,ngroup,nsweep,speed,
	2                                  length,fakeit
	                  Write(lun1,*) grp_start,grp_stop,grp_num,
	2                                  grp_thre
	                  Write(lun1,110) (gf_hist(i), i=1,512)
	                  If (ifg_flag) Then
	                     Write(lun1,110) (gf_ifg(i), i=1,512)
 110	                     Format(1x, 8f8.2)
	                  End If
	               End If  !i_write
	            End If  !linear

	         End If !num > 0

	      End Do  !(.Not. eof)

	      If (ans(1:1) .Eq. 'Q') Then
	         i_stop = .True.
	      End If

	   End Do !i_stop

	End Do !ich
C
C *** Does the user want to continue?
C
	Write(6,10)
 10	Format(' Would you like to select another channel? Y/[N]: ', $)
	Read(5,11) response
 11	Format(a)
	status = STR$UpCase(response, response)
	If (response .Eq. 'Y') Then
	   Write(6,20)
 20	   Format(' Please enter desired channel : ', $)
	   Read(5,12) chan_id
 12	   Format(i4)
	   Goto 100
C
C **** Does the user want another time range?
C
	Else
 150	   Continue
	   Write(6,13)
 13	   Format(' Would you like to select another time? Y/[N]: ', $)
	   Read(5,14) ask
	   status = STR$UpCase(ask, ask)
 14	   Format(a)
	   If (ask .Eq. 'Y') Then
	      If (n_chan .Gt. 1) Then
	         Goto 50
	      Else
	         Write(6,20)
	         Read(5,12) chan_id
	         Goto 50
	      End If
	   End If
	End If
	Call LIB$Erase_Page(1,1)
	If (i_write) Then
	   Close(UNIT=lun1)
	End If
C
C *** EXIT WITH STATUS MESSAGE
C
	If (status .Eq. SS$_Normal) Then
	   Call LIB$Signal(FSD_Normal, %Val(1), %Val(status))
	Else
	   Call LIB$Signal(FSD_ABErr3, %Val(1), %Val(status))
	End If

	Stop
	End
